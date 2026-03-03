import torch
import math
import numpy as np
from typing import Mapping, Union
from loss import _linear_loss, _spline_loss, _log_mvn_density, _pen_mat

###############################################################################

def fit_gaussian_spline(
	NX: np.ndarray,
	B: np.ndarray,
	Z: np.ndarray,
    lam: float,
    control: Mapping[str, Union[float, int]]
) -> np.ndarray:
	"""
	Fit a dynamic Gaussian copula model using smooth splines.
	
	Parameters
	----------
	NX : np.ndarray
		Normal-transformed pseudo-observations. Shape `(N, d)`. Rows
		correspond to `x`.
	B : np.ndarray
        Basis matrix. Shape `(N, K)`.
    Z : np.ndarray
        Design matrix of categorical covariates. Shape `(N, L)` or `None`.
	lam : float
	    Smoothing parameter.
	control :  Mapping[str, Union[float, int]]
		Optimization control parameters.

	Returns
	-------
	beta : np.ndarray
		Estimated coefficient matrix. Shape `(K + L, p)`.
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

	# Set up tensors
	NX_t = torch.tensor(NX, dtype = dtype)
	B_t = torch.tensor(B, dtype = dtype)
	
	if Z is None:
		L = 0
		Z_t = None
	else:
		L = Z.shape[1]
		Z_t = torch.tensor(Z, dtype = dtype)
    
	# Set initial coefficients all to zero
	beta = torch.zeros(
        size = (K + L, d * (d - 1) // 2),
        dtype = dtype,
        requires_grad = True
    )

	# Precompute fixed penalty matrix
	S = _pen_mat(K = K, dtype = dtype)
    
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
		    NX = NX_t,
		    B = B_t,
		    Z = Z_t,
		    S = S,
		    lam = lam,
		    dtype = dtype
		)
		loss.backward()
		return loss
    
	loss = optimizer.step(closure)
	return beta.detach().numpy()

###############################################################################

def fit_gaussian_linear(
	NX: np.ndarray,
	Z: np.ndarray,
    control: Mapping[str, Union[float, int]]
) -> np.ndarray:
	"""
	Fit a Gaussian GLM copula model.
	
	Parameters
	----------
	NX : np.ndarray
		Normal-transformed pseudo-observations. Shape `(N, d)`. Rows
		correspond to `x`.
    Z : np.ndarray
        Design matrix of categorical covariates. Shape `(N, L)` or `None`.
	control :  Mapping[str, Union[float, int]]
		Optimization control parameters.

	Returns
	-------
	beta : np.ndarray
		Estimated coefficient matrix. Shape `(L, p)`.
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
	    
	d = NX.shape[1]
	L = Z.shape[1]

	# Set up tensors
	NX_t = torch.tensor(NX, dtype = dtype)
	Z_t = torch.tensor(Z, dtype = dtype)
    
	# Set initial coefficients all to zero
	beta = torch.zeros(
        size = (L, d * (d - 1) // 2),
        dtype = dtype,
        requires_grad = True
    )
    
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
		loss = _linear_loss(
		    beta = beta,
		    NX = NX_t,
		    Z = Z_t,
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
	Z: np.ndarray,
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
    Z : np.ndarray
        Design matrix of categorical covariates. Shape `(N, L)` or `None`.
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
	
	# Indices for training and testing data
	test_ix = np.arange(min_test_ix, max_test_ix + 1).astype(int)
	train_mask = np.ones(N, dtype = bool)
	train_mask[test_ix] = False
	train_ix = np.arange(N)[train_mask]

	# Set up tensors
	NX_t = torch.tensor(NX, dtype = dtype)
	NX_train = NX_t[train_ix, :]
	NX_test = NX_t[test_ix, :]
	
	B_t = torch.tensor(B, dtype = dtype)
	B_train = B_t[train_ix, :]
	B_test = B_t[test_ix, :]
	
	if Z is None:
		L = 0
		Z_t = None
		Z_train = None
		Z_test = None
	else:
		L = Z.shape[1]
		Z_t = torch.tensor(Z, dtype = dtype)
		Z_train = Z_t[train_ix, :]
		Z_test = Z_t[test_ix, :]

	# Set initial coefficients all to zero
	beta = torch.zeros(
        size = (K + L, d * (d - 1) // 2),
        dtype = dtype,
        requires_grad = True
    )

	# Precompute fixed penalty matrix
	S = _pen_mat(K = K, dtype = dtype)
    
    # Fit using training data
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
		    NX = NX_train,
		    B = B_train,
		    Z = Z_train,
		    S = S,
		    lam = lam,
		    dtype = dtype
		)
		loss.backward()
		return loss
    
	try:
		loss = optimizer.step(closure)
		beta = beta.detach()
	except RuntimeError:
	    # Something went wrong, fall back to large negative log-likelihood
		return -1e10

	# Predicted coefficients for test data
	H_test = B_test @ beta[:K, :]
	if Z is not None:
		H_test = H_test + Z_test @ beta[K:, :]
	
	# Marginal log-likelihood
	log2pi = math.log(2 * math.pi)
	margin_ll = (-0.5 * (log2pi + NX_test * NX_test)).sum(dim = 1)
	# Copula log-likelihood
	copula_ll = _log_mvn_density(NX_test, H_test, dtype)
	# Average model log-likelihood
	ll = torch.sum(copula_ll - margin_ll) / margin_ll.shape[0]
	
	return ll.numpy()

###############################################################################
