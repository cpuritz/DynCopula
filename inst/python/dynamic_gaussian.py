import torch
import botorch
import numpy as np
from distributions import _vec2chol, _log_mvn_density

######################################################################

def fit_continuous_gaussian(par0, dx, NX, band, control):
	max_it = int(control["max_it"])
	reltol = control["reltol"]
	patience = int(control["patience"])

	eta = torch.tensor(
		par0,
		dtype = torch.float64,
		requires_grad = True
	)
	dx = np.array(dx)
	NX = torch.tensor(NX, dtype = eta.dtype)

	optimizer = torch.optim.RMSprop(
		[eta],
		lr = control["lr"],
		weight_decay = control["weight_decay"]
	)
	
	# Record loss history
	hist = []
	# Track last good value in case of error
	eta_good = eta.detach().clone()
	
	# Exit codes
	#  0  = converged
	#  1 = reached max iterations
	#  2 = error occurred
	exit_code = 1
	
	for i in range(max_it):
		optimizer.zero_grad()
		loss = _loglik_cts_gaussian(eta, dx, NX, band)
		hist.append(loss.item())

		# Check for convergence
		if i >= patience:
			hp = hist[i - patience]
			lhist = hist[(i - patience + 1):(i + 1)]
			if all((hp - h) / abs(hp) < reltol for h in lhist):
				exit_code = 0
				break
		loss.backward()
		optimizer.step()
		if torch.isnan(eta).any():
			exit_code = 2
			break
		eta_good = eta.detach().clone()

	opt = eta_good.detach().numpy()
	return {"par": opt, "hist": hist, "convergence": int(exit_code)}
	
######################################################################
	
def _loglik_cts_gaussian(eta, dx, NX, band):
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
		
		# Reconstruct Cholesky factor of correlation matrix
		L = _vec2chol(v)
		
		# Log probability
		# Note: this is not the Gaussian copula density, but rather the only
		# part of it that depends on the correlation matrix. The additional
		# term is -ndist.log_prob(u).sum().
		ndist = torch.distributions.Normal(loc = 0, scale = 1)
		u = NX[pos_ix[i], :]
		P[i] = _log_mvn_density(u, L)

	# Weighted negative log likelihood
	return -torch.sum(wgt * P)

######################################################################

def fit_discrete_gaussian(par0, dx, NX, NXm, band, control):
	max_it = int(control["max_it"])
	reltol = control["reltol"]
	patience = int(control["patience"])

	eta = torch.tensor(
		par0,
		dtype = torch.float64,
		requires_grad = True
	)
	dx = np.array(dx)
	NX = torch.tensor(NX, dtype = eta.dtype)
	NXm = torch.tensor(NXm, dtype = eta.dtype)

	optimizer = torch.optim.RMSprop(
		[eta],
		lr = control["lr"],
		weight_decay = control["weight_decay"]
	)

	# Record loss history
	hist = []
	# Track last good value in case of error
	eta_good = eta.detach().clone()
	
	'''
	Exit codes:
	  0  = converged
	  1 = reached max iterations
	  2 = error occurred
	'''
	exit_code = 1
	
	for i in range(max_it):
		optimizer.zero_grad()
		loss = _loglik_discrete_gaussian(eta, dx, NXm, NX, band)
		hist.append(loss.item())

		# Check for convergence
		if i >= patience:
			hp = hist[i - patience]
			lhist = hist[(i - patience + 1):(i + 1)]
			if all(abs(h - hp) / hp < reltol for h in lhist):
				exit_code = 0
				break
		loss.backward()
		optimizer.step()
		if torch.isnan(eta).any():
			exit_code = 2
			break
		eta_good = eta.detach().clone()

	opt = eta_good.detach().numpy()
	return {"par": opt, "hist": hist, "convergence": int(exit_code)}
	
######################################################################

def _loglik_discrete_gaussian(eta, dx, NXm, NX, band):
	d = int(1 + np.sqrt(1 + 4 * eta.numel()) / 2)
	
	wgt = 3 / (4 * band) * np.maximum(1 - (dx / band) ** 2, 0)
	pos_ix  = np.where(wgt > 0)[0]

	P = torch.zeros(len(pos_ix), dtype = eta.dtype)
	npar = int(eta.numel() / 2)
	
	for i in range(len(pos_ix)):
		v = eta[0:npar] + eta[npar:(2 * npar)] * dx[pos_ix[i]]
		# Reconstruct Cholesky factor
		L = _vec2chol(v)
		# Reconstruct correlation matrix
		R = L.matmul(L.t())

		# Compute log probability
		bounds = np.column_stack((NXm[pos_ix[i], :], NX[pos_ix[i], :]))
		domain = torch.tensor(bounds, dtype = eta.dtype)
		P[i] = botorch.utils.probability.MVNXPB(R, domain).solve()

	# Weighted negative log likelihood
	wgt_tens = torch.tensor(wgt[pos_ix], dtype = eta.dtype)
	nll = -torch.sum(wgt_tens * P)
	return nll

######################################################################
