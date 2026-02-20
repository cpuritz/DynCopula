import torch
import math
from functools import lru_cache

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
	L = _vec2chol(V, dtype)
	
	# Compute m = L^(-1) X
	m = torch.linalg.solve_triangular(L, X.unsqueeze(-1), upper = False)
	print("min m:", m.min().item())
	print("max m:", m.max().item())
	
	# Mahalanobis distance between X and the Gaussian copula specified by L
	M = (m * m).sum(dim = -2).squeeze(-1)
	print("min M:", M.min().item())
	print("max M:", M.max().item())

	# Compute 0.5 * log(det(LL^T))
	diag = L.diagonal(dim1 = -2, dim2 = -1)
	print("min diag:", diag.min().item())
	print("max diag:", diag.max().item())
	half_log_det = diag.clamp_min(torch.finfo(dtype).eps).log().sum(-1)
	print("min hld:", half_log_det.min().item())
	print("max hld:", half_log_det.max().item())
	
	return -0.5 * (d * math.log(2 * math.pi) + M) - half_log_det

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
        Shape `(N, d(d-1)/2)`.
    dtype : torch.dtype
        Floating-point precision.

    Returns
    -------
    torch.Tensor
        A tensor of shape `(N, d, d)`, which each batch representing a separate
        Cholesky factor.
    """

    N, npar = V.shape

    d = (1 + math.isqrt(1 + 8 * npar)) // 2
    rows, cols, mask = _tril_col_major(d)

    # Base identity stacked for batch
    H = torch.eye(d, dtype = dtype).expand(N, d, d).clone()

    # Fill strictly lower triangular entries
    H[:, rows, cols] = torch.tanh(0.5 * V)

    # Compute cumulative product term
    X = H[:, :, :-1].pow(2).clamp_max(1 - torch.finfo(dtype).eps)
    csum = torch.cumsum(torch.log1p(-X) * mask[:, :-1], dim = 2)
    # Clamp min sqrtcprod at 1e-10
    sqrtcprod = torch.exp(0.5 * csum.clamp_min(math.log(1e-20)))
    
    # Build Cholesky factor
    L = torch.zeros((N, d, d), dtype = dtype)
    L[:, :, 0] = H[:, :, 0]
    L[:, :, 1:] = H[:, :, 1:] * sqrtcprod

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
