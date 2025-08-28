import torch
from distributions import (
    _log_mvn_density,
    _log_mvn_mass,
    _local_loglik_cts,
    _local_loglik_count
)

###############################################################################

def aic_cts(eta_i, x, NX, i, h):
    eta_i = eta_i.detach().clone().requires_grad_(True)

    # Local likelihood function
    def L_local(par):
    	return _local_loglik_cts(par, x, x[i], NX, h)

    # Model likelihood function
    def L_i(par):
        return _log_mvn_density(NX[i, :], par)

    # Degrees of freedom
    W0 = 3 / (4 * h)
    J = torch.autograd.functional.hessian(L_local, eta_i)
    H_i = torch.autograd.functional.hessian(L_i, eta_i)
    nu = W0 * torch.trace(torch.linalg.solve(-J, -H_i))

    dev = -2 * L_i(eta_i)
    df = 2 * nu
    
    return dev.detach().numpy(), df.detach().numpy()

###############################################################################

def aic_count(eta_i, x, NXm, NX, i, h):
    eta_i = eta_i.detach().clone().requires_grad_(True)
    
    # Local likelihood function
    def L_local(par):
    	return _local_loglik_count(par, x, x[i], NXm, NX, h)
    
    # Model likelihood function
    def L_i(par):
        return _log_mvn_mass(NXm[i, ].unsqueeze(0), NX[i, ].unsqueeze(0), par)
    
    u = (x - x[i]) / h
    w = (1.0 - u * u).clamp_min(0.0) * (3.0 / (4.0 * h))
    pos_ix = torch.nonzero(w > 0, as_tuple = False).flatten()
    
    p = eta_i.numel()
    K = torch.zeros((p, p), dtype = torch.float64)
    for j in pos_ix.tolist():
        def L_j(par):
            return _log_mvn_mass(NXm[j, ].unsqueeze(0), NX[j, ].unsqueeze(0), par)
        sj = torch.autograd.grad(L_j(eta_i), eta_i, retain_graph = True)[0].reshape(p, 1)
        K = K + (w[j] * w[j]) * (sj @ sj.T)
    J = -torch.autograd.functional.hessian(L_local, eta_i)
    W0 = 3 / (4 * h)
    c_i = W0 / ((w * w).sum() + 1e-12)
    nu1 = c_i * torch.trace(torch.linalg.solve(J, K))

    # Degrees of freedom
    # W0 = 3 / (4 * h)
    # J = torch.autograd.functional.hessian(L_local, eta_i)
    # H_i = torch.autograd.functional.hessian(L_i, eta_i)
    # nu1 = W0 * torch.trace(torch.linalg.solve(-J, -H_i))
    
    dev = (-2 * L_i(eta_i)).detach().numpy()
    df = (2 * nu1).detach().numpy()
    return dev, df

###############################################################################
