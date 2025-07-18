import torch
import botorch
import numpy as np
from distributions import _vec2chol, _log_mvt_density

######################################################################

def fit_continuous_t(par0, nu, dx, TX, band, control):
	max_it = int(control["max_it"])
	reltol = control["reltol"]
	patience = int(control["patience"])

	eta = torch.tensor(
		par0,
		dtype = torch.float64,
		requires_grad = True
	)
	dx = np.array(dx)
	TX = torch.tensor(TX, dtype = eta.dtype)

	optimizer = torch.optim.SGD(
	    [eta],
	    lr = control["lr"],
	    momentum = control["momentum"],
	    weight_decay = control["weight_decay"],
	    nesterov = True
	)

	# Record loss history
	hist = []
	# Track last good value in case of error
	eta_good = eta.detach().clone()
	# eta history
	eta_hist = []

	'''
	Exit codes:
	  0  = converged
	  1 = reached max iterations
	  2 = error occurred
	'''
	exit_code = 1

	for i in range(max_it):
		optimizer.zero_grad()
		loss = _loglik_cts_t(eta, dx, TX, band, nu)
		hist.append(loss.item())
		eta_hist.append(eta.detach().clone().unsqueeze(1))

		# Check for convergence
		if i >= patience:
			hp = hist[i - patience]
			lhist = hist[(i - patience + 1):(i + 1)]
			if all(abs(hp - h) <= reltol * abs(hp) for h in lhist):
				exit_code = 0
				break
		loss.backward()
		optimizer.step()
		if torch.isnan(eta).any():
			exit_code = 2
			break
		eta_good = eta.detach().clone()
		
    opt = eta_good.detach().numpy()
    eta_hist = torch.cat(eta_hist, dim = 1).numpy()
    return {
        "par": opt,
        "loss_hist": hist,
        "convergence": int(exit_code),
        "eta_hist": eta_hist
    }

######################################################################

def _loglik_cts_t(eta, dx, TX, band, nu):
	d = int(1 + np.sqrt(1 + 4 * eta.numel()) / 2)

	# Kernel weights
	wgt = 3 / (4 * band) * np.maximum(1 - (dx / band) ** 2, 0)
	pos_ix  = np.where(wgt > 0)[0]
	wgt = torch.tensor(wgt[pos_ix], dtype = eta.dtype)

	P = torch.zeros(len(pos_ix), dtype = eta.dtype)
	npar = int(eta.numel() / 2)
	loc = torch.zeros(d, dtype = eta.dtype)

	for i in range(len(pos_ix)):
		v = eta[0:npar] + eta[npar:(2 * npar)] * dx[pos_ix[i]]
		# Reconstruct Cholesky factor
		L = _vec2chol(v)

		# Log probability
		u = TX[pos_ix[i], :]
		P[i] = _log_mvt_density(u, L, nu)

	# Weighted negative log likelihood
	return -torch.sum(wgt * P)

######################################################################
