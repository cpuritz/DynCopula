import torch
import botorch
import numpy as np

######################################################################

def fit_continuous(par0, dx, NX, scale, band, control):
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
	
	d = int(1 + np.sqrt(1 + 4 * eta.numel()) / 2)
	indices = np.stack(np.tril_indices(d, k = -1), axis = 1)
	sorted_indices = indices[np.lexsort((indices[:, 0], indices[:, 1]))]
	l_tri = (torch.tensor(sorted_indices[:, 0], dtype = torch.int32),
		   torch.tensor(sorted_indices[:, 1], dtype = torch.int32))

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
		loss = _loglik_cts(eta, dx, NX, scale, band, l_tri)
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
	return {"opt": opt, "hist": hist, "convergence": int(exit_code)}
	
######################################################################
	
def _loglik_cts(eta, dx, NX, scale, band, l_tri):
	d = int(1 + np.sqrt(1 + 4 * eta.numel()) / 2)	

	wgt = 3 / (4 * band) * np.maximum(1 - (dx / band) ** 2, 0)
	pos_ix  = np.where(wgt > 0)[0]
	P = torch.zeros(len(pos_ix), dtype = eta.dtype)
	npar = int(eta.numel() / 2)

	for i in range(len(pos_ix)):
		v = eta[0:npar] + eta[npar:(2 * npar)] * dx[pos_ix[i]]

		# Transform back to constrained space
		H = torch.eye(d, dtype = eta.dtype)
		H = H.index_put(l_tri, torch.tanh(v / scale))
		
		# Reconstruct Cholesky factor
		L = torch.zeros(d, d, dtype = eta.dtype)
		L[0, :] = H[0, :]
		L[:, 0] = H[:, 0]
		for j in range(1, d):
			cprod = torch.cumprod(1 - H[j, 0:j]**2, dim = 0)
			L[j, 1:(j + 1)] = H[j, 1:(j + 1)] * torch.sqrt(cprod)
		
		# Reconstruct correlation matrix
		R = L.matmul(L.t())

		# Compute log probability
		mvn = torch.distributions.multivariate_normal.MultivariateNormal(
		    loc = torch.zeros(d, dtype = eta.dtype),
		    covariance_matrix = R,
		    validate_args = False
        )
		P[i] = mvn.log_prob(NX[pos_ix[i], :])

	# Weighted negative log likelihood
	wgt_tens = torch.tensor(wgt[pos_ix], dtype = eta.dtype)
	nll = -torch.sum(wgt_tens * P)
	return nll
	
######################################################################

def fit_discrete(par0, dx, NX, NXm, scale, band, control):
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
	
	d = int(1 + np.sqrt(1 + 4 * eta.numel()) / 2)
	indices = np.stack(np.tril_indices(d, k = -1), axis = 1)
	sorted_indices = indices[np.lexsort((indices[:, 0], indices[:, 1]))]
	l_tri = (torch.tensor(sorted_indices[:, 0], dtype = torch.int32),
		   torch.tensor(sorted_indices[:, 1], dtype = torch.int32))

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
		loss = _loglik_discrete(eta, dx, NXm, NX, scale, band, l_tri)
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
	return {"opt": opt, "hist": hist, "convergence": int(exit_code)}
	
######################################################################

def _loglik_discrete(eta, dx, NXm, NX, scale, band, l_tri):
	d = int(1 + np.sqrt(1 + 4 * eta.numel()) / 2)
	
	wgt = 3 / (4 * band) * np.maximum(1 - (dx / band) ** 2, 0)
	pos_ix  = np.where(wgt > 0)[0]

	P = torch.zeros(len(pos_ix), dtype = eta.dtype)
	npar = int(eta.numel() / 2)
	
	for i in range(len(pos_ix)):
		v = eta[0:npar] + eta[npar:(2 * npar)] * dx[pos_ix[i]]

		# Transform back to constrained space
		H = torch.eye(d, dtype = eta.dtype)
		H = H.index_put(l_tri, torch.tanh(v / scale))

		# Reconstruct Cholesky factor
		L = torch.zeros(d, d, dtype = eta.dtype)
		L[0, :] = H[0, :]
		L[:, 0] = H[:, 0]
		for j in range(1, d):
			cprod = torch.cumprod(1 - H[j, 0:j]**2, dim = 0)
			L[j, 1:(j + 1)] = H[j, 1:(j + 1)] * torch.sqrt(cprod)

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
