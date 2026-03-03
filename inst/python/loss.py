import torch
import math
from corr_utils import _vec2chol

###############################################################################

def _linear_loss(
	beta: torch.Tensor,
	NX: torch.tensor,
	Z: torch.Tensor,
	dtype: torch.dtype
) -> torch.Tensor:
	"""
	Compute the loss for a Gaussian GLM copula.

	Parameters
	----------
	beta : torch.Tensor
        Coefficients. Shape `(L, p)`.
	NX : torch.tensor
	    Normal-transformed pseudo-observations. Shape `(N, d)`.
    Z : torch.Tensor
        Design matrix. Shape `(N, L)`.
	dtype : torch.dtype
	    Floating-point precision.

	Returns
	-------
	torch.Tensor
		A scalar tensor representing the loss.
	"""
	
	eta = Z @ beta
	nll = -torch.sum(_log_mvn_density(X = NX, V = eta, dtype = dtype))
	return nll

###############################################################################

def _spline_loss(
	beta: torch.Tensor,
	NX: torch.tensor,
	B: torch.Tensor,
	Z: torch.Tensor,
	S: torch.Tensor,
	lam: float,
	dtype: torch.dtype
) -> torch.Tensor:
	"""
	Compute the spline loss for a Gaussian copula.

	Parameters
	----------
	beta : torch.Tensor
        Coefficients. Shape `(K + L, p)`.
	NX : torch.tensor
	    Normal-transformed pseudo-observations. Shape `(N, d)`.
	B : torch.Tensor
	    Basis matrix. Shape `(N, K)`.
    Z : torch.Tensor
        Design matrix of shape `(N, L)` or `None`.
	S : torch.Tensor
	    Penalty matrix. Shape `(K, K)`.
	lam : float
	    Smoothing parameter.
	dtype : torch.dtype
	    Floating-point precision.

	Returns
	-------
	torch.Tensor
		A scalar tensor representing the loss.
	"""
	
	K = B.shape[1]
	# First K rows are for continuous covariate, remaining rows are for
	# categorical covariate
	beta_cts = beta[:K, :]
	eta = B @ beta_cts
	if Z is not None:
		eta = eta + Z @ beta[K:, :]

	# Penalty term = lambda/2 * trace(beta.T * S * beta)
	pen = 0.5 * lam * (beta_cts * (S @ beta_cts)).sum()

	# Negative log likelihood
	nll = -torch.sum(_log_mvn_density(X = NX, V = eta, dtype = dtype))

	return nll + pen

###############################################################################

def _log_mvn_density(
    X: torch.Tensor,
    V: torch.Tensor,
    dtype: torch.dtype
) -> torch.Tensor:
	"""
	Compute the log-density of a multivariate Gaussian distribution.

 	Parameters
	----------
	X : torch.Tensor
        Input samples of shape `(N, d)`, where `d` is the dimensionality.
    V : torch.Tensor
        Either of shape `(npar,)` or `(npar, N)`. In the former case, one
        covariance matrix is constructed for all samples. In the latter case,
        one covariance matrix is constructed for each sample. `npar` must equal
        `choose(d, 2)`.
    dtype : torch.dtype
        Floating-point precision.

	Returns
	-------
	torch.Tensor
		The log-density of each sample under the parameterized multivariate
		Gaussian.
	"""
	
	d = X.shape[-1]
	
	# Convert the unconstrained parameter vector to a Cholesky factor
	L = _vec2chol(V, d, dtype)
	
	# Compute m = L^(-1) X
	m = torch.linalg.solve_triangular(L, X.unsqueeze(-1), upper = False)

	# Mahalanobis distance between X and the Gaussian copula specified by L
	M = (m * m).sum(dim = -2).squeeze(-1)
	
	# Compute 0.5 * log(det(LL^T))
	diag = L.diagonal(dim1 = -2, dim2 = -1)
	half_log_det = diag.clamp_min(torch.finfo(dtype).eps).log().sum(-1)
	
	return -0.5 * (d * math.log(2 * math.pi) + M) - half_log_det

###############################################################################

def _pen_mat(K: int, dtype: torch.dtype) -> torch.Tensor:
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
    S = torch.zeros((K, K), dtype = dtype)
    ix = torch.arange(K)
    
    # Main diagonal: [1, 5, 6, ..., 6, 5, 1]
    S[ix, ix] = 6
    S[0, 0] = 1
    S[1, 1] = 5
    S[-2, -2] = 5
    S[-1, -1] = 1
    
    # 1-off diagonals: [-2, -4, ..., -4, -2]
    S[ix[:-1], ix[1:]] = -4
    S[ix[1:], ix[:-1]] = -4
    S[0, 1] = -2
    S[1, 0] = -2
    S[-2, -1] = -2
    S[-1, -2] = -2
    
    # 2-off diagonals: [1, 1, ..., 1]
    S[ix[:-2], ix[2:]] = 1
    S[ix[2:], ix[:-2]] = 1
    
    return S

###############################################################################
