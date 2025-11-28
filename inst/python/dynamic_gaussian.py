import torch
import numpy as np
import math
from functools import lru_cache

###############################################################################

def fit_gaussian(
	par0: np.ndarray,
	x: np.ndarray,
	NX: np.ndarray,
	h: float,
	control: Mapping[str, float | int],
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
	control : Mapping[str, float | int]
		Optimization control parameters.
	x0 : float
		Covariate value to estimate copula parameters at.

	Returns
	-------
	np.ndarray
		A 1D array containing the estimates of the calibration coefficients.
	"""

	max_epoch = int(control["max_epoch"])
	max_iter = int(control["max_itr"])
	history_size = int(control["history_size"])
	tolerance_grad = float(control["tolerance_grad"])
	tolerance_change = float(control["tolerance_change"])
	
	eta = torch.tensor(
		np.atleast_1d(par0).tolist(),
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
		loss = -1.0 * _local_loglik(eta, x, x0, NX, h)
		loss.backward()
		return loss

	for _ in range(max_epoch):
		loss = optimizer.step(closure)
		
	return eta.detach().numpy()

###############################################################################
	
def _local_loglik(
	eta: torch.Tensor,
	x: torch.Tensor,
	x0: torch.Tensor,
	NX: torch.Tensor,
	h: float
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
		A scalar specifying the covariate value to estimate copula parameters at.
	NX : torch.Tensor
		Normal-transformed pseudo-observations. Must be tensor of size `(n x d)`,
		where `n` is the length of `x` (rows correspond to values in `x`).
	h : float
		Kernel bandwidth.

	Returns
	-------
	torch.Tensor
		A scalar tensor representing the local log-likelihood of `eta`.
	"""

	# Epanechnikov kernel weights
	u = (x - x0) / h
	wgt = 3 / (4 * h) * torch.clamp(1 - u * u, min = 0)
	# Record where weights are nonzero to avoid unnecessary likelihood evaluations
	mask = wgt > 0

	# Model log-likelihoods. Not the exact log-likelihoods for a Gaussian copula,
	# but the other terms don't depend on the correlation matrix and thus are
	# irrelevant during estimation.
	P = _log_mvn_density(NX[mask, :], eta)

	# Local log-likelihood
	return torch.dot(wgt[mask], P)

###############################################################################

def _log_mvn_density(x: torch.Tensor, v: torch.Tensor) -> torch.Tensor:
	"""
	Compute the log-density of a multivariate Gaussian distribution.

 	Parameters
	----------
	x : torch.Tensor
		Input samples of shape `(..., d)`, where `d` is the dimensionality.
	v : torch.Tensor
		A vector of unconstrained parameters of length `choose(d, 2)`
		parametrizing the covariance matrix.

	Returns
	-------
	torch.Tensor
		The log-density of each sample under the parameterized multivariate
		Gaussian. Output shape matches the batch dimensions of `x`.
	"""

	d = x.shape[-1]

	# Convert the unconstrained parameter vector to a Cholesky factor
	L = _vec2chol(v)
	# Compute m = L^(-1) x
	m = torch.linalg.solve_triangular(L, x.unsqueeze(-1), upper = False)
	# Mahalanobis distance between x and the Gaussian copula specified by L
	M = (m * m).sum(dim = -2)[..., 0]

	# Compute 0.5 * log|R| for R=LL^T
	diag = L.diagonal(dim1 = -2, dim2 = -1)
	half_log_det = diag.clamp_min(torch.finfo(torch.float64).eps).log().sum(-1)

	log2pi = x.new_tensor(2.0 * math.pi).log()
	return -0.5 * (d * log2pi + M) - half_log_det

###############################################################################

def _vec2chol(
	v: torch.Tensor,
	rho_max: float = 0.9999,
	scale: float = 0.5
) -> torch.Tensor:
	"""
    Map a vector of unconstrained values to a valid Cholesky factor.

    Parameters
    ----------
    v : torch.Tensor
        A 1D tensor.
    rho_max : float, optional
        Maximum allowed magnitude for off-diagonal correlations.
    scale : float, optional
		Scaling factor to control steepness of `tanh` transformation.

    Returns
    -------
    torch.Tensor
        A lower-triangular `(d x d)` matrix representing a Cholesky factor.
    """

	# Length of v determines the dimension of the matrix. The length of v
	# should equal choose(d, 2) for a positive integer d.
	d = (1 + math.isqrt(1 + 8 * v.numel())) // 2

	# Fill strictly lower-triangular entries of H in column-major order
	r, c, mask = _tril_col_major(d)
	H = torch.eye(d, dtype = torch.float64)
	H[r, c] = rho_max * torch.tanh(scale * v)

	eps = 1e-12
	X = H[:, :-1].pow(2).clamp_max(1 - eps)
	logS = torch.log1p(-X) * mask.to(torch.float64)
	sqrtcprod = torch.exp(0.5 * torch.cumsum(logS, dim = 1))

	L = torch.zeros((d, d), dtype = torch.float64)
	L[:, 0] = H[:, 0]
	L[:, 1:] = H[:, 1:] * sqrtcprod
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

	# Lower triangular indices in column major order
	rows, cols = torch.tril_indices(d, d, offset = -1)
	order = torch.argsort(cols * d + rows)

	r = torch.arange(d).unsqueeze(1)
	c = torch.arange(d - 1).unsqueeze(0)
	mask = (c < r)
	
	return rows[order], cols[order], mask

###############################################################################
