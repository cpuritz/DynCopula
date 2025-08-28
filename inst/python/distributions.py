import torch
import math
from functools import lru_cache

###############################################################################
	
def _local_loglik(
    eta: torch.Tensor,
    x: torch.Tensor,
    x0: torch.Tensor,
    NX: torch.Tensor,
    h: float
) -> torch.Tensor:
	# Epanechnikov kernel weights
	u = (x0 - x) / h
	wgt = 3 / (4 * h) * torch.clamp(1 - u * u, min = 0)
	mask = wgt > 0

	# Log likelihoods
	P = _log_mvn_density(NX[mask, :], eta)
		
	# Local log likelihood
	return torch.dot(wgt[mask], P)

###############################################################################

def _log_mvn_density(x: torch.Tensor, v: torch.Tensor) -> torch.Tensor:
    d = x.shape[-1]
    L = _vec2chol(v)
    M = _mahalanobis(x, L)
    diag = L.diagonal(dim1 = -2, dim2 = -1)
    half_log_det = diag.clamp_min(torch.finfo(torch.float64).eps).log().sum(-1)
    log2pi = x.new_tensor(2.0 * math.pi).log()
    return -0.5 * (d * log2pi + M) - half_log_det

###############################################################################

def _mahalanobis(x: torch.Tensor, L: torch.Tensor) -> torch.Tensor:
    # x: (..., d), L: (..., d, d)
    m = torch.linalg.solve_triangular(L, x.unsqueeze(-1), upper = False)
    return (m * m).sum(dim = -2)[..., 0]

###############################################################################

def _vec2chol(
    v: torch.Tensor,
    rho_max: float = 0.99,
    scale: float = 0.5
) -> torch.Tensor:
    d = (1 + math.isqrt(1 + 8 * v.numel())) // 2

    # Fill strictly lower-triangular entries of H in column-major order
    r, c, mask = _tril_col_major(d)
    H = torch.eye(d, dtype = torch.float64)
    H[r, c] = rho_max * torch.tanh(scale * v)
    
    eps = 1e-12
    X = H[:, :-1].pow(2).clamp_max(1 - eps)
    logS = torch.log1p(-X) * mask.to(torch.float64)
    logcprod = torch.cumsum(logS, dim = 1)
    sqrtcprod = torch.exp(0.5 * logcprod)
    
    L = torch.zeros((d, d), dtype = torch.float64)
    L[:, 0] = H[:, 0]
    L[:, 1:] = H[:, 1:] * sqrtcprod
    return L

###############################################################################

@lru_cache(maxsize = 1)
def _tril_col_major(d: int):
    # Lower triangular indices in column major order
    rows, cols = torch.tril_indices(d, d, offset = -1)
    order = torch.argsort(cols * d + rows)
    
    r = torch.arange(d).unsqueeze(1)
    c = torch.arange(d - 1).unsqueeze(0)
    mask = (c < r)
    
    return rows[order], cols[order], mask

###############################################################################
