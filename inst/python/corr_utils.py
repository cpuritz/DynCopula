import torch
import math
from functools import lru_cache

###############################################################################

def _vec2chol(
    V: torch.Tensor,
    d: int,
    dtype: torch.dtype
) -> torch.Tensor:
    """
    Map a vector of unconstrained values to a valid Cholesky factor.

    Parameters
    ----------
    V : torch.Tensor
        Shape `(N, d(d-1)/2)`.
    d : int
        Dimension.
    dtype : torch.dtype
        Floating-point precision.

    Returns
    -------
    torch.Tensor
        A tensor of shape `(N, d, d)`, which each batch representing a separate
        Cholesky factor.
    """

    N, npar = V.shape
    rows, cols, mask = _tril_col_major(d)

    # Base identity stacked for batch
    H = torch.eye(d, dtype = dtype).expand(N, d, d).clone()

    # Fill strictly lower triangular entries
    H[:, rows, cols] = torch.tanh(0.5 * V)

    # Compute cumulative product term
    X = H[:, :, :-1].pow(2).clamp_max(1 - torch.finfo(dtype).eps)
    csum = torch.cumsum(torch.log1p(-X) * mask[:, :-1], dim = 2)
    # Clamp minimum at 1e-10
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
