import torch
import torch.nn.functional as F
import math
import numpy as np
from typing import Mapping, Union
from spline_loss import _log_mvn_density, _spline_loss, _pen_mat

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
		Normal-transformed pseudo-observations. Shape `(n, d)`. Rows
		correspond to `x`.
	B : np.ndarray
        Basis matrix. Shape `(n, K)`.
    Z : np.ndarray
        Design matrix of categorical covariates. Shape `(n, L)` or `None`.
        Entries must be integers with a minimum value of `0` in each column.
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

	# Set up tensors
	NX = torch.tensor(NX, dtype = dtype)
	B = torch.tensor(B, dtype = dtype)
	
	# Categorical covariates
	if Z is not None:
	    Z_mat = one_hot_encode(Z.astype(int), dtype = dtype)
	    npar_cat = Z_mat.shape[1]
	else:
	    Z_mat = None
	    npar_cat = 0
    
	# Set initial coefficients all to zero
	beta = torch.zeros(
        size = (K + npar_cat, d * (d - 1) // 2),
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
		    NX = NX,
		    B = B,
		    Z = Z_mat,
		    S = S,
		    lam = lam,
		    dtype = dtype
		)
		loss.backward()
		return loss
    
	loss = optimizer.step(closure)
	beta = beta.detach()
	
	# Estimated model coefficients
	eta_hat = B @ beta[:K, :]
	if Z_mat is not None:
		eta_hat = eta_hat + Z_mat @ beta[K:, :]

	return eta_hat.numpy()
    
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
        Design matrix of categorical covariates. Shape `(n, L)` or `None`.
        Entries must be integers with a minimum value of `0` in each column.
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
	NX_train = torch.tensor(NX[train_ix, :], dtype = dtype)
	NX_test = torch.tensor(NX[test_ix, :], dtype = dtype)
	B_train = torch.tensor(B[train_ix, :], dtype = dtype)
	B_test = torch.tensor(B[test_ix, :], dtype = dtype)
	
	# Categorical covariates
	if Z is not None:
        Z_mat = one_hot_encode(Z.astype(int), dtype = dtype)
        npar_cat = Z_mat.shape[1]
	else:
    	Z_mat = None
        npar_cat = 0
	
	# Set initial coefficients all to zero
	beta = torch.zeros(
        size = (K + npar_cat, d * (d - 1) // 2),
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
		    Z = Z_mat,
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
	H_test = B_test @ beta
	
	# Estimated model coefficients
	H_test = B_test @ beta[:K, :]
	if Z_mat is not None:
		H_test = H_test + Z_mat @ beta[K:, :]
	
	# Marginal log-likelihood
	log2pi = math.log(2 * math.pi)
	margin_ll = (-0.5 * (log2pi + NX_test * NX_test)).sum(dim = 1)
	# Copula log-likelihood
	copula_ll = _log_mvn_density(NX_test, H_test, dtype)
	# Average model log-likelihood
	ll = torch.sum(copula_ll - margin_ll) / margin_ll.shape[0]
	
	return ll.numpy()

###############################################################################

def one_hot_encode(
    Z: np.ndarray,
    dtype: torch.dtype
) -> torch.Tensor:
    """
    Construct sparse one-hot encodings for multiple categorical covariates.

    Parameters
    ----------
    Z : np.ndarray
        Integer-coded categorical labels. Minimum must be `0`.
	dtype : torch.dtype
	    Floating-point precision.

    Returns
    -------
    Z_mat : torch.Tensor (sparse COO)
        Binary tensor indicating covariate values.
    """
    
    # Number of covariates
    N, ncat = Z.shape
    # Number of unique values per covariate
    Ls = np.array([(Z[:, k].max().item() + 1) for k in range(ncat)])
    
    labels = torch.as_tensor(Z, dtype = torch.int64)
    Z_list = []
    for k in range(ncat):
        lab = labels[:, k]
        # Use 0 as reference
        mask = (lab != 0)
        rows = torch.arange(N, dtype = torch.int64)[mask]
        cols = (lab[mask] - 1).to(torch.int64)
        
        Zk = torch.sparse_coo_tensor(
            indices = torch.stack([rows, cols]),
            values = torch.ones(rows.numel(), dtype = dtype),
            size = (N, Ls[k] - 1),
            dtype = dtype
        ).coalesce()
        Z_list.append(Zk)
        
    return torch.cat(Z_list, dim = 1).coalesce()

###############################################################################
