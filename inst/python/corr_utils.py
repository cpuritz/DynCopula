import torch
import math
from functools import lru_cache

###############################################################################

def _log_mvn_density(
    x: torch.Tensor,
    V: torch.Tensor,
    dtype: torch.dtype
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
    dtype : torch.dtype
        Floating-point precision.

	Returns
	-------
	torch.Tensor
		The log-density of each sample under the parameterized multivariate
		Gaussian.
	"""

	d = x.shape[-1]
	
	# Convert the unconstrained parameter vector to a Cholesky factor
	L = _vec2chol(V, dtype)
	
	# Compute m = L^(-1) x
	m = torch.linalg.solve_triangular(L, x.unsqueeze(-1), upper = False)
	
	# Mahalanobis distance between x and the Gaussian copula specified by L
	M = (m * m).sum(dim = -2).squeeze(-1)

	# Compute 0.5*log|R| for R=LL^T
	diag = L.diagonal(dim1 = -2, dim2 = -1)
	half_log_det = diag.clamp_min(torch.finfo(dtype).eps).log().sum(-1)

	log2pi = x.new_tensor(2.0 * math.pi).log()
	return -0.5 * (d * log2pi + M) - half_log_det

###############################################################################

def _gaussian_cop_loglik(
    x: torch.Tensor,
    V: torch.Tensor,
    dtype: torch.dtype
) -> torch.Tensor:
	"""
	Compute the log-density of a multivariate Gaussian copula.

 	Parameters
	----------
	x : torch.Tensor
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
		The total log-density over all samples under the parameterized
		multivariate Gaussian copula.
	"""
	# Copula log-likelihood
	copula_ll = _log_mvn_density(
	    x = x,
	    V = V,
	    dtype = dtype
	)
	
	# Marginal log-likelihood
	log2pi = x.new_tensor(2.0 * math.pi).log()
	margin_ll = (-0.5 * (x * x + log2pi)).sum(dim = 1)

	# Sum over rows
	return (copula_ll - margin_ll).sum()

###############################################################################

def _vec2chol(
    V: torch.Tensor,
    dtype: torch.dtype
) -> torch.Tensor:
    """
    Map a vector of unconstrained values to a valid Cholesky factor.

    Parameters
    ----------
    V : torch.Tensor
        Either of shape `(npar,)` or `(npar, N)`.
    dtype : torch.dtype
        Floating-point precision.

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
    H = torch.eye(d, dtype = dtype).expand(nbatch, d, d).clone()

    # Fill strictly lower triangular entries
    H[:, r, c] = torch.tanh(0.5 * V.T)

    # Compute cumulative product term
    X = H[:, :, :-1].pow(2).clamp_max(1 - torch.finfo(dtype).eps)
    logS = torch.log1p(-X) * mask[:, :-1]
    sqrtcprod = torch.exp(0.5 * torch.cumsum(logS, dim = 2))

    # Build Cholesky factor
    L = torch.zeros((nbatch, d, d), dtype = dtype)
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
