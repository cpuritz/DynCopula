import torch
import numpy as np
import math
from functools import lru_cache
from typing import Mapping, Union

###############################################################################

def fit_gaussian(
	par0: np.ndarray,
	x: np.ndarray,
	NX: np.ndarray,
	h: float,
	degree: int,
	control: Mapping[str, Union[float, int]],
	x0: float
) -> float:
	"""
	Fit a Gaussian copula model using local likelihood.

	Parameters
	----------
	par0 : np.ndarray
		Initial values for the calibration coefficients.
	x : np.ndarray
		Covariate values. Shape `(n,)`.
	NX : np.ndarray
		Normal-transformed pseudo-observations. Shape `(n, d)`. Rows
		correspond to `x`.
	h : float
		Kernel bandwidth.
	degree : int
        Degree of local polynomial approximation.
	control : Mapping[str, float | int]
		Optimization control parameters.
	x0 : float
		Covariate value to estimate copula parameters at.

	Returns
	-------
	np.ndarray
		A 1D array containing the estimates of the calibration coefficients.
	"""

	max_outer = int(control["max_outer"])
	max_iter = int(control["max_itr"])
	history_size = int(control["history_size"])
	tolerance_grad = float(control["tolerance_grad"])
	tolerance_change = float(control["tolerance_change"])
	degree = int(degree)
	par0 = np.atleast_1d(par0)
	
	# Set initial derivative values to zero
	npar = len(par0)
	par0 = np.concatenate([par0, np.zeros(degree * npar)])
	
	# Set up tensors
	eta = torch.tensor(
		par0.tolist(),
		dtype = torch.float64,
		requires_grad = True
	)
	x = torch.tensor(x, dtype = torch.float64)
	x0 = torch.tensor(x0, dtype = torch.float64)
	NX = torch.tensor(NX, dtype = torch.float64)
	
	optimizer = torch.optim.LBFGS(
		[eta],
		line_search_fn = "strong_wolfe",
		max_iter = max_iter,
		history_size = history_size,
		tolerance_grad = tolerance_grad,
		tolerance_change = tolerance_change
	)
	
	def closure():
		optimizer.zero_grad()
		loss = -1.0 * _local_loglik(eta, x, x0, NX, h, degree)
		loss.backward()
		return loss

	for _ in range(max_outer):
		loss = optimizer.step(closure)
	par_est = eta.detach().numpy()
	
	# Only interested in function estimates, not derivative estimates
	par_est = par_est[0:npar]
	
	return par_est

###############################################################################
	
def _local_loglik(
	eta: torch.Tensor,
	x: torch.Tensor,
	x0: torch.Tensor,
	NX: torch.Tensor,
	h: float,
	degree: int
) -> torch.Tensor:
	"""
	Compute the local log-likelihood for a Gaussian copula.

	Parameters
	----------
	eta : torch.Tensor
		A vector of unconstrained parameters of length `choose(d, 2)`
		parametrizing the correlation matrix.
	x : torch.Tensor
		A 1D tensor of covariate values.
	x0 : torch.Tensor
		A scalar specifying the covariate value to estimate copula parameters
		at.
	NX : torch.Tensor
		Normal-transformed pseudo-observations. Must be tensor of size
		`(n x d)`, where `n` is the length of `x` (rows correspond to values in
		`x`).
	h : float
		Kernel bandwidth.
	degree : int
	    Degree of local polynomial approximation.

	Returns
	-------
	torch.Tensor
		A scalar tensor representing the local log-likelihood of `eta`.
	"""

	# Epanechnikov kernel weights
	dx = x - x0
	wgt = 3 / (4 * h) * torch.clamp(1 - (dx / h)**2, min = 0)
	# Record where weights are nonzero to avoid unnecessary likelihood
	# evaluations
	mask = wgt > 0

	# Compute model log-likelihoods. Note that these are not the exact
	# log-likelihoods for a Gaussian copula, but the other terms don't depend
	# on the correlation matrix and thus can safely be ignored.
	powers = torch.arange(0, degree + 1).unsqueeze(0)
	dxp = dx.unsqueeze(1) ** powers
	eta_local = torch.matmul(dxp, eta.view(degree + 1, -1))
	if degree == 0:
	    # Coefficients are identical, no need to compute more than one
	    # correlation matrix
	    eta_local = eta_local[0, :]
	else:
	    eta_local = eta_local[mask, :].T
	P = _log_mvn_density(NX[mask, :], eta_local)

	# Local log-likelihood
	return torch.dot(wgt[mask], P)

###############################################################################

def _log_mvn_density(
    x: torch.Tensor,
    V: torch.Tensor
) -> torch.Tensor:
	"""
	Compute the log-density of a multivariate Gaussian distribution.

 	Parameters
	----------
	x : torch.Tensor
        Input samples of shape `(N, d)`, where `d` is the dimensionality.
    V : torch.Tensor
        Either of shape `(npar,)` or `(npar, N)`. In the former case, one
        covariance matrix is constructed for all samples. In the latter case,
        one covariance matrix is constructed for each sample. `npar` must equal
        `choose(d, 2)`.

	Returns
	-------
	torch.Tensor
		The log-density of each sample under the parameterized multivariate
		Gaussian.
	"""

	d = x.shape[-1]
	
	# Convert the unconstrained parameter vector to a Cholesky factor
	L = _vec2chol(V)
	
	# Compute m = L^(-1) x
	m = torch.linalg.solve_triangular(L, x.unsqueeze(-1), upper = False)
	
	# Mahalanobis distance between x and the Gaussian copula specified by L
	M = (m * m).sum(dim = -2).squeeze(-1)

	# Compute 0.5*log|R| for R=LL^T
	diag = L.diagonal(dim1 = -2, dim2 = -1)
	half_log_det = diag.clamp_min(torch.finfo(torch.float64).eps).log().sum(-1)

	log2pi = x.new_tensor(2.0 * math.pi).log()
	return -0.5 * (d * log2pi + M) - half_log_det

###############################################################################

def _vec2chol(
    V: torch.Tensor,
    scale: float = 0.5
) -> torch.Tensor:
    """
    Map a vector of unconstrained values to a valid Cholesky factor.

    Parameters
    ----------
    V : torch.Tensor
        Either of shape `(npar,)` or `(npar, N)`.
    scale : float, optional
		Scaling factor to control steepness of `tanh` transformation.

    Returns
    -------
    torch.Tensor
        If V was 1D, then a 2D tensor of shape `(d, d)` representing a Cholesky
        factor. Otherwise, a 3D tensor of shape `(N, d, d)`, which each batch
        representing a separate Cholesky factor.
    """
    
    # Add batch dimension
    if V.ndim == 1:
        V = V.unsqueeze(-1)
    npar, nbatch = V.shape
    
    d = (1 + math.isqrt(1 + 8 * npar)) // 2
    r, c, mask = _tril_col_major(d)

    # Base identity stacked for batch
    H = torch.eye(d, dtype = torch.float64).expand(nbatch, d, d).clone()

    # Fill strictly lower triangular entries
    H[:, r, c] = torch.tanh(scale * V.T)

    # Compute cumulative product term
    X = H[:, :, :-1].pow(2).clamp_max(1 - torch.finfo(torch.float64).eps)
    logS = torch.log1p(-X) * mask[:, :-1]
    sqrtcprod = torch.exp(0.5 * torch.cumsum(logS, dim = 2))

    # Build Cholesky factor
    L = torch.zeros((nbatch, d, d), dtype = torch.float64)
    L[:, :, 0] = H[:, :, 0]
    L[:, :, 1:] = H[:, :, 1:] * sqrtcprod
    
    # If initially unbatched, remove batch dimension
    L = L.squeeze(0)
    
    return L

###############################################################################

@lru_cache(maxsize = 1)
def _tril_col_major(d: int):
	"""
    Compute strictly lower-triangular matrix indices of a `d x d` matrix in
	column-major order.

    Parameters
    ----------
    d : int
        Dimension of the square matrix.

    Returns
    -------
    rows, cols : torch.Tensor
        1D tensors containing the lower-triangular row/column index pairs.
    mask : torch.Tensor
        Boolean mask where `mask[i, j]` is True if `j < i`.
    """
	rows, cols = torch.tril_indices(d, d, offset = -1)
	# Convert from row-major order to column-major order
	order = torch.argsort(cols * d + rows)
	rows = rows[order]
	cols = cols[order]
	
	mask = torch.zeros((d, d), dtype = torch.bool)
	mask[rows, cols] = True
	return rows, cols, mask

###############################################################################
