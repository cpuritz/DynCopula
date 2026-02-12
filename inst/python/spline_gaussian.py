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
    control: Mapping[str, Union[float, int]],
    compute_edf: bool
) -> Mapping[str, Union[np.ndarray, float]]:
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
	B : np.ndarray
        Basis matrix. Shape `(n, K)`.
	lam : float
	    Smoothing parameter.
	control :  Mapping[str, Union[float, int]]
		Optimization control parameters.
	compute_edf : bool
	    Whether to compute degrees of freedom.

	Returns
	-------
	A dictionary with the following elements:
	    
	beta : np.ndarray
		Estimated spline coefficients. Shape `(K, p)`.
	edf : float
	    Effective degrees of freedom. If `compute_edf` is `False`, then the
	    value returned is `np.nan`.
	"""

	max_outer = int(control["max_outer"])
	max_iter = int(control["max_itr"])
	history_size = int(control["history_size"])
	tolerance_grad = float(control["tolerance_grad"])
	tolerance_change = float(control["tolerance_change"])
	dtype = control["precision"]
	lam = float(lam)
	compute_edf = bool(compute_edf)
	
	if dtype == "float32":
	    dtype = torch.float32
	else:
	    dtype = torch.float64

	# Set up tensors
	par0 = np.atleast_1d(par0)
	par0 = np.ascontiguousarray(par0)
	beta = torch.tensor(par0, dtype = dtype, requires_grad = True)
	
	x = torch.tensor(x, dtype = dtype)
	NX = torch.tensor(NX, dtype = dtype)
	B = torch.tensor(B, dtype = dtype)

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
    
	for _ in range(max_outer):
		loss = optimizer.step(closure)
	beta_hat = beta.detach()
	
	if compute_edf:
	    edf = edf_pen(beta_hat, NX, B, lam, dtype)
	else:
	    edf = np.nan
	
	return {
	    "beta": beta_hat.numpy(),
	    "edf": edf
	}
    
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
	nll = -torch.sum(_log_mvn_density(
	    x = NX,
	    V = eta.T.contiguous(),
	    dtype = dtype
	))
	
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

def edf_pen(
    beta_hat: torch.tensor,
    NX: torch.tensor,
    B: torch.tensor,
    lam: float,
    dtype: torch.dtype
) -> float:
    """
	Compute effective degrees of freedom (EDF).

	Parameters
	----------
	beta_hat : torch.tensor
		Estimated model coefficients.
	NX : torch.tensor
	    Normal-transformed pseudo-observations.
	B : torch.tensor
	    Spline basis matrix.
	lam : float
	    Smoothing parameter.
	dtype : torch.dtype
	    Floating point precision.
		
	Returns
	-------
	float
		A float representing the EDF.
	"""

    K, p = beta_hat.shape

    # Penalty matrix S (K x K) and P = I_p ⊗ S  (Kp x Kp)
    S = pen_mat(K = K, dtype = dtype)
    P = torch.kron(torch.eye(p, dtype = dtype), S) # (Kp, Kp)

    # Flatten beta so autograd can take Hessians in R^{Kp}
    beta0 = beta_hat.reshape(-1).clone().detach().requires_grad_(True)

    def ll_from_vec(beta_vec: torch.Tensor) -> torch.Tensor:
        beta = beta_vec.view(K, p)
        eta = B @ beta
        ll = torch.sum(_log_mvn_density(
            x = NX,
            V = eta.T.contiguous(),
            dtype = dtype
        ))
        return ll

    # Observed information matrix
    Iobs = -torch.autograd.functional.hessian(ll_from_vec, beta0)  # (Kp, Kp)

    # EDF = tr((Iobs + lam * P)^{-1} Iobs)
    edf = torch.trace(torch.linalg.solve(Iobs + lam * P, Iobs))
    
    return edf.item()
