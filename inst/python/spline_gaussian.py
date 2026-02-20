import torch
import math
import numpy as np
from typing import Mapping, Union
from corr_utils import _log_mvn_density

###############################################################################

def fit_gaussian_spline(
	NX: np.ndarray,
	B: np.ndarray,
    lam: float,
    control: Mapping[str, Union[float, int]]
) -> np.ndarray:
	"""
	Fit a dynamic Gaussian copula model using smooth splines.
	
	Parameters
	----------
	NX : np.ndarray
		Normal-transformed pseudo-observations. Shape `(n, d)`. Rows
		correspond to `x`.
	B : np.ndarray
        Basis matrix. Shape `(n, K)`.
	lam : float
	    Smoothing parameter.
	control :  Mapping[str, Union[float, int]]
		Optimization control parameters.

	Returns
	-------
	beta : np.ndarray
		Estimated spline coefficients. Shape `(K, p)`.
	"""
	
	max_iter = int(control["max_itr"])
	history_size = int(control["history_size"])
	tolerance_grad = float(control["tolerance_grad"])
	tolerance_change = float(control["tolerance_change"])
	dtype = control["precision"]

	if dtype == "float32":
	    dtype = torch.float32
	else:
	    dtype = torch.float64
	    
	K = B.shape[1]
	d = NX.shape[1]
	npar = d * (d - 1) // 2

	# Set up tensors
	NX = torch.tensor(NX, dtype = dtype)
	B = torch.tensor(B, dtype = dtype)
	
	# Set initial coefficients all to zero
	beta = torch.zeros(K, npar, dtype = dtype, requires_grad = True)

	# Precompute fixed penalty matrix
	S = pen_mat(K = B.shape[1], dtype = dtype)
    
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
		loss = _spline_loss(
		    beta = beta,
		    NX = NX,
		    B = B,
		    S = S,
		    lam = lam,
		    dtype = dtype
		)
		loss.backward()
		return loss
    
	loss = optimizer.step(closure)

	return beta.detach().numpy()
    
###############################################################################
    
def gaussian_spline_cv(
	NX: np.ndarray,
	B: np.ndarray,
    lam: float,
    control: Mapping[str, Union[float, int]],
    min_test_ix: int,
    max_test_ix: int
) -> float:
	"""
	Fit a dynamic Gaussian copula model using smooth splines.
	
	Parameters
	----------
	NX : np.ndarray
		Normal-transformed pseudo-observations. Shape `(N, d)`. Rows
		correspond to `x`.
	B : np.ndarray
        Basis matrix. Shape `(N, K)`.
	lam : float
	    Smoothing parameter.
	control :  Mapping[str, Union[float, int]]
		Optimization control parameters.
	min_test_ix : int
	    Minimum index for testing data.
	max_test_ix : int
	    Maximum index for testing data
	
	Returns
	-------
	ll : float
	    Cross-validated log-likelihood.
	"""
	
	max_iter = int(control["max_itr"])
	history_size = int(control["history_size"])
	tolerance_grad = float(control["tolerance_grad"])
	tolerance_change = float(control["tolerance_change"])
	dtype = control["precision"]

	if dtype == "float32":
	    dtype = torch.float32
	else:
	    dtype = torch.float64
	    
	K = B.shape[1]
	N, d = NX.shape
	npar = d * (d - 1) // 2
	
	# Indices for training and testing data
	test_ix = np.arange(min_test_ix, max_test_ix + 1).astype(int)
	train_mask = np.ones(N, dtype = bool)
	train_mask[test_ix] = False
	train_ix = np.arange(N)[train_mask]

	# Set up tensors
	NX_train = torch.tensor(NX[train_ix, :], dtype = dtype)
	NX_test = torch.tensor(NX[test_ix, :], dtype = dtype)
	B_train = torch.tensor(B[train_ix, :], dtype = dtype)
	B_test = torch.tensor(B[test_ix, :], dtype = dtype)
	
	# Set initial coefficients all to zero
	beta = torch.zeros(K, npar, dtype = dtype, requires_grad = True)

	# Precompute fixed penalty matrix
	S = pen_mat(K = B.shape[1], dtype = dtype)
    
    # Fit using training data
	optimizer = torch.optim.LBFGS(
		[beta],
		line_search_fn = "strong_wolfe",
		max_iter = max_iter,
		history_size = history_size,
		tolerance_grad = tolerance_grad,
		tolerance_change = tolerance_change
	)
	
	last_good = {"beta": beta.detach().clone(), "loss": None}
    
	def closure():
		optimizer.zero_grad()
		loss = _spline_loss(
		    beta = beta,
		    NX = NX_train,
		    B = B_train,
		    S = S,
		    lam = lam,
		    dtype = dtype
		)
		print("loss:", loss.item())
		
		if torch.isfinite(loss):
			last_good["beta"] = beta.detach().clone()
			last_good["loss"] = loss.detach()
		else:
			raise RuntimeError("Non-finite loss")
        
		loss.backward()
		return loss
    
	try:
		loss = optimizer.step(closure)
		beta_hat = beta.detach()
	except RuntimeError:
		print("Something bad happened, using last good value")
		loss = last_good["loss"]
		beta_hat = last_good["beta"]

	# Predicted coefficients for test data
	H_test = B_test @ beta_hat
	
	# Marginal log-likelihood
	margin_ll = (-0.5 * (math.log(2 * math.pi) + NX_test * NX_test)).sum(dim = 1)
	# Copula log-likelihood
	copula_ll = _log_mvn_density(NX_test, H_test, dtype)
	# Average model log-likelihood
	ll = torch.sum(copula_ll - margin_ll) / margin_ll.shape[0]
	
	return ll.numpy()
    
###############################################################################
    
def _spline_loss(
	beta: torch.Tensor,
	NX: torch.tensor,
	B: torch.Tensor,
	S: torch.Tensor,
	lam: float,
	dtype: torch.dtype
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
	dtype : torch.dtype
	    Floating-point precision.
	
	Returns
	-------
	torch.Tensor
		A scalar tensor representing the spline loss.
	"""
	
	# Parameter estimates
	eta = B @ beta
	
	# Negative log likelihood
	nll = -torch.sum(_log_mvn_density(NX, eta, dtype))

	# Penalty term = lambda/2 * trace(beta.T * S * beta)
	pen = 0.5 * lam * (beta * (S @ beta)).sum()

	# Spline loss
	loss = nll + pen
	
	return loss

###############################################################################

def pen_mat(K: int, dtype: torch.dtype) -> torch.Tensor:
	"""
	Compute penalty matrix.

	Parameters
	----------
	K : int
		Degrees of freedom.
	dtype : torch.dtype
	    Floating point precision.
		
	Returns
	-------
	torch.Tensor
		A tensor representing the penalty matrix.
	"""
	
	K = int(K)
	I = torch.eye(K, dtype = dtype)
	D2 = I[:-2] - 2 * I[1:-1] + I[2:]
	S = D2.T @ D2
	return S

###############################################################################
