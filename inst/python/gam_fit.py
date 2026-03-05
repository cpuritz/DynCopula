import torch
import math
import numpy as np
from typing import Mapping, Union
from loss import _gam_loss, _log_mvn_density, _pen_mat

###############################################################################

def fit_gaussian_gam(
	NX: np.ndarray,
	B: np.ndarray,
	Z: np.ndarray,
	M: np.ndarray,
    lam: np.ndarray,
    control: Mapping[str, Union[float, int]]
) -> np.ndarray:
	"""
	Fit a GAM Gaussian copula model.
	
	Parameters
	----------
	NX : np.ndarray
		Normal-transformed pseudo-observations. Shape `(N, d)`. Rows
		correspond to `x`.
	B : np.ndarray
        Basis matrix. Shape `(N, K)`.
    Z : np.ndarray
        Design matrix of categorical covariates. Shape `(N, L2)` or `None`.
    M : np.ndarray
        Design matrix for interaction between the smooth covariate and discrete
        covariates. Shape `(N, L1)`.
	lam : np.ndarray
	    Penalty parameters. Shape `(L1 + 1,)`.
	control :  Mapping[str, Union[float, int]]
		Optimization control parameters.

	Returns
	-------
	beta : np.ndarray
		Estimated coefficient matrix. Shape `(K * L1 + L2, p)`.
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
	    
	N, d = NX.shape
	K = B.shape[1]
	L1 = M.shape[1]
	L2 = Z.shape[1]
	p = d * (d - 1) // 2

	# Set up tensors
	NX_t = torch.tensor(NX, dtype = dtype)
	B_t = torch.tensor(B, dtype = dtype)
	M_t = torch.tensor(M, dtype = dtype)
	Z_t = torch.tensor(Z, dtype = dtype)
	lam = torch.tensor(lam, dtype = dtype)
    
	# Set initial coefficients all to zero
	beta = torch.zeros(
        size = (K * L1 + L2, p),
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
		loss = _gam_loss(
		    beta = beta,
		    NX = NX_t,
		    B = B_t,
		    Z = Z_t,
		    M = M_t,
		    S = S,
		    lam = lam,
		    dtype = dtype
		)
		loss.backward()
		return loss
    
	loss = optimizer.step(closure)
	return beta.detach().numpy()

###############################################################################
    
def gaussian_gam_cv(
	NX: np.ndarray,
	B: np.ndarray,
	Z: np.ndarray,
	M: np.ndarray,
    lam: np.ndarray,
    control: Mapping[str, Union[float, int]],
    min_test_ix: int,
    max_test_ix: int
) -> float:
	"""
	Compute cross-validated log-likelihood for a GLM Gaussian copula model.
	
	Parameters
	----------
	NX : np.ndarray
		Normal-transformed pseudo-observations. Shape `(N, d)`. Rows
		correspond to `x`.
	B : np.ndarray
        Basis matrix. Shape `(N, K)`.
    Z : np.ndarray
        Design matrix of categorical covariates. Shape `(N, L2)` or `None`.
    M : np.ndarray
        Design matrix for interaction between the smooth covariate and discrete
        covariates. Shape `(N, L1)`.
	lam : np.ndarray
	    Penalty parameters. Shape `(L1 + 1,)`.
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
	    
	N, d = NX.shape
	K = B.shape[1]
	L1 = M.shape[1]
	L2 = Z.shape[1]
	p = d * (d - 1) // 2
	
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
	
	M_t = torch.tensor(M, dtype = dtype)
	M_train = M_t[train_ix, :]
	M_test = M_t[test_ix, :]
	
	Z_t = torch.tensor(Z, dtype = dtype)
	Z_train = Z_t[train_ix, :]
	Z_test = Z_t[test_ix, :]
	
	lam = torch.tensor(lam, dtype = dtype)

	# Set initial coefficients all to zero
	beta = torch.zeros(
        size = (K * L1 + L2, p),
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
		loss = _gam_loss(
		    beta = beta,
		    NX = NX_train,
		    B = B_train,
		    Z = Z_train,
		    M = M_train,
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
	betas = beta[L2:, :].view(L1, K, p)              # (L1, K, p)
	s_per_j = B_test.unsqueeze(0) @ betas            # (L1, N, p)
	H_test = torch.einsum('ij,jip->ip', M_test, s_per_j)  # (N, p)
	if Z is not None:
		H_test = H_test + Z_test @ beta[:L2, :]
	
	# Marginal log-likelihood
	log2pi = math.log(2 * math.pi)
	margin_ll = (-0.5 * (log2pi + NX_test * NX_test)).sum(dim = 1)
	# Copula log-likelihood
	copula_ll = _log_mvn_density(NX_test, H_test, dtype)
	# Average model log-likelihood
	ll = torch.sum(copula_ll - margin_ll) / margin_ll.shape[0]
	
	return ll.numpy()

###############################################################################
