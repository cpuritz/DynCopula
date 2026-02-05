import torch
import numpy as np
from typing import Mapping, Union
from corr_utils import _log_mvn_density

###############################################################################

def fit_gaussian_spline(
	par0: np.ndarray,
	x: np.ndarray,
	NX: np.ndarray,
	B: np.ndarray,
    lam: float,
    control: Mapping[str, Union[float, int]]
) -> float:
	"""
	Fit a dynamic Gaussian copula model using smooth splines.
	
	Parameters
	----------
	par0 : np.ndarray
		Initial values for the spline coefficients. Shape `(K, p)`.
	x : np.ndarray
		Covariate values. Shape `(n,)`.
	NX : np.ndarray
		Normal-transformed pseudo-observations. Shape `(n, d)`. Rows
		correspond to `x`.
	B: np.ndarray
        Basis matrix. Shape `(n, K)`.
	lam: float
	    Smoothing parameter vector.
	control : Mapping[str, float | int]
		Optimization control parameters.

	Returns
	-------
	np.ndarray
		A 1D array containing the estimates of the spline coefficients.
	"""

	max_outer = int(control["max_outer"])
	max_iter = int(control["max_itr"])
	history_size = int(control["history_size"])
	tolerance_grad = float(control["tolerance_grad"])
	tolerance_change = float(control["tolerance_change"])
	lam = float(lam)

	# Set up tensors
	par0 = np.atleast_1d(par0)
	par0 = np.ascontiguousarray(par0)
	beta = torch.tensor(par0, dtype = torch.float64, requires_grad = True)
	
	x = torch.tensor(x, dtype = torch.float64)
	NX = torch.tensor(NX, dtype = torch.float64)
	B = torch.tensor(B, dtype = torch.float64)

	# Precompute fixed penalty matrix
	S = pen_mat(B.shape[1])
    
	optimizer = torch.optim.LBFGS(
		[beta],
		line_search_fn = "strong_wolfe",
		max_iter = max_iter,
		history_size = history_size,
		tolerance_grad = tolerance_grad,
		tolerance_change = tolerance_change
	)
    
	def closure():
		optimizer.zero_grad()
		loss = _spline_loss(beta, NX, B, S, lam)
		loss.backward()
		return loss
    
	for _ in range(max_outer):
		loss = optimizer.step(closure)

	return beta.detach().numpy()
    
###############################################################################
    
def _spline_loss(
	beta: torch.Tensor,
	NX: torch.tensor,
	B: torch.Tensor,
	S: torch.Tensor,
	lam: float
) -> torch.Tensor:
	"""
	Compute the spline loss for a Gaussian copula.
	
	Parameters
	----------
	beta : torch.Tensor
        Spline coefficients. Shape `(K, p)`.
	NX : torch.tensor
	    Normal-transformed pseudo-observations. Shape `(n, d)`.
	B : torch.Tensor
	    Basis matrix. Shape `(n, K)`.
	S : torch.Tensor
	    Penalty matrix. Shape `(K, K)`.
	lam : float
	    Smoothing parameter.
	
	Returns
	-------
	torch.Tensor
		A scalar tensor representing the spline loss.
	"""

	# Parameter estimates
	eta = B @ beta
	
	# Negative log likelihood
	nll = -torch.sum(_log_mvn_density(NX, eta.T.contiguous()))
	
	# Penalty term = lambda/2 * trace(beta.T * S * beta)
	pen = 0.5 * lam * (beta * (S @ beta)).sum()
	
	# Spline loss
	loss = nll + pen
	
	return loss

###############################################################################

def pen_mat(K: int) -> torch.Tensor:
	"""
	Compute penalty matrix.

	Parameters
	----------
	K : int
		Degrees of freedom.
		
	Returns
	-------
	torch.Tensor
		A tensor representing the penalty matrix.
	"""
	K = int(K)
	I = torch.eye(K, dtype = torch.float64)
	D2 = I[:-2] - 2 * I[1:-1] + I[2:]
	S = D2.T @ D2
	return S

###############################################################################
