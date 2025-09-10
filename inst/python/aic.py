import torch
from distributions import _log_mvn_density, _local_loglik

###############################################################################

def compute_aic(eta_i, x, NX, i, h):
    eta_i = eta_i.detach().clone().requires_grad_(True)

    # Local likelihood function
    def L_local(par):
    	return _local_loglik(par, x, x[i], NX, h)

    # Model likelihood function
    def L_i(par):
        return _log_mvn_density(NX[i, :], par)
    
    # Deviance
    dev = L_i(eta_i)

    # Degrees of freedom
    W0 = 3 / (4 * h)
    J = torch.autograd.functional.hessian(L_local, eta_i)
    H_i = torch.autograd.functional.hessian(L_i, eta_i)
    nu = W0 * torch.trace(torch.linalg.solve(-J, -H_i))
    
    aic = -2 * dev + 2 * nu
    return aic.detach().numpy()

###############################################################################
