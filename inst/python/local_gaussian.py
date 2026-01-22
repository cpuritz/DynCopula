import torch
import numpy as np
from typing import Mapping, Union
from corr_utils import _log_mvn_density

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
