import torch
import numpy as np
import math

######################################################################

def _log_mvn_density(x, L):
    d = x.size(0)
    M = _mahalanobis(L, x)
    half_log_det = L.diagonal(dim1 = -2, dim2 = -1).log().sum(-1)
    return -0.5 * (d * math.log(2 * math.pi) + M) - half_log_det

######################################################################

def _log_mvt_density(x, L, nu):
    d = x.size(0)
    M = _mahalanobis(L, x)
    half_log_det = L.diagonal(dim1 = -2, dim2 = -1).log().sum(-1)
    df = torch.tensor(nu, dtype = torch.float64)
    return (
        torch.lgamma((df + d) / 2)
        - torch.lgamma(df / 2)
        - 0.5 * d * math.log(df * math.pi)
        - half_log_det
        - 0.5 * (df + d) * torch.log1p(M / df)
    )

######################################################################

def _mahalanobis(L, x):
    d = x.size(0)
    flat_L = L.reshape(-1, d, d)
    flat_x = x.view(1, d, 1)
    M = torch.linalg.solve_triangular(flat_L, flat_x, upper = False)
    return M.pow(2).sum(1).squeeze()

######################################################################

def _vec2chol(v):
    # Lower triangular indices in column major order
    d = int(1 + np.sqrt(1 + 8 * v.numel()) / 2)
    indices = np.stack(np.tril_indices(d, k = -1), axis = 1)
    sorted_indices = indices[np.lexsort((indices[:, 0], indices[:, 1]))]
    l_tri = (
        torch.tensor(sorted_indices[:, 0], dtype = torch.int32),
        torch.tensor(sorted_indices[:, 1], dtype = torch.int32)
    )

    # Transform back to constrained space
    H = torch.eye(d, dtype = v.dtype)
    H = H.index_put(l_tri, torch.tanh(v))

    # Reconstruct Cholesky factor
    L = torch.zeros(d, d, dtype = v.dtype)
    L[0, :] = H[0, :]
    L[:, 0] = H[:, 0]
    for j in range(1, d):
        cprod = torch.cumprod(1 - H[j, 0:j]**2, dim = 0)
        L[j, 1:(j + 1)] = H[j, 1:(j + 1)] * torch.sqrt(cprod)
    return L

######################################################################
