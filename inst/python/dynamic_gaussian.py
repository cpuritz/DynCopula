import torch
import numpy as np
from distributions import _local_loglik_cts, _local_loglik_count
from aic import aic_cts, aic_count

###############################################################################

def fit_gaussian_cts(par0, x, NX, h, control, x0 = None, i = None):
	max_itr = int(control["max_itr"])
	patience = int(control["patience"])
	reltol = float(control["reltol"])

	if i is not None:
	    i = int(i)
	    x0 = x[i]
	
	eta = torch.tensor(
		np.atleast_1d(par0).tolist(),
		dtype = torch.float64,
		requires_grad = True
	)
	x = torch.tensor(x, dtype = torch.float64)
	x0 = torch.tensor(x0, dtype = torch.float64)
	NX = torch.tensor(NX, dtype = torch.float64)

	nesterov = (control["momentum"] != 0)
	optimizer = torch.optim.SGD(
		[eta],
		lr = control["lr"],
		momentum = control["momentum"],
		nesterov = nesterov
	)
	
	# Record loss history
	hist = []
	# Track last good value in case of error
	eta_good = eta.detach().clone()
	# eta history
	eta_hist = []
	
	'''
	Exit codes:
	  0 = converged
	  1 = reached max iterations
	  2 = error occurred
	'''
	exit_code = 1
	
	for j in range(max_itr):
		optimizer.zero_grad()
		loss = -1.0 *  _local_loglik_cts(eta, x, x0, NX, h)
		with torch.no_grad():
			hist.append(loss.item())
			eta_hist.append(eta.detach().clone().unsqueeze(1))

		# Check for convergence
		if j >= patience:
			eps = 1e-12
			rel = [
				abs(hist[j - k] - hist[j - k - 1]) / (abs(hist[j - k - 1]) + eps)
				for k in range(patience)
			]
			if (all(r >= 0 and r <= reltol for r in rel)):
			    if len(hist) == patience + 1:
			        # No improvement has been made, increase learning rate
			        # before stopping
			        lr = optimizer.param_groups[0]["lr"]
			        optimizer.param_groups[0]["lr"] = 10 * lr
			    else:
				    exit_code = 0
				    break

		loss.backward()
		torch.nn.utils.clip_grad_norm_(eta, max_norm = control["max_grad"])
		optimizer.step()
		if torch.isnan(eta).any():
			exit_code = 2
			break
		eta_good = eta.clone()

	par_opt = eta_good.detach().numpy()
	eta_hist = torch.cat(eta_hist, dim = 1).numpy()
	
	# Estimate AIC
	if i is not None:
	    dev, df = aic_cts(eta_good, x, NX, i, h)
	else:
	    dev = None
	    df = None

	return {
        "par": par_opt,
        "loss_hist": hist,
        "convergence": int(exit_code),
        "eta_hist": eta_hist,
        "deviance": dev,
        "df": df
    }
    
###############################################################################

def fit_gaussian_count(par0, x, NXm, NX, h, control, x0 = None, i = None):
	max_epoch = int(control["max_epoch"])
	max_itr = int(control["max_itr"])
	history_size = int(control["history_size"])
	tolerance_grad = float(control["tolerance_grad"])
	tolerance_change = float(control["tolerance_change"])
	
	if i is not None:
	    i = int(i)
	    x0 = x[i]
	
	eta = torch.tensor(
		np.atleast_1d(par0).tolist(),
		dtype = torch.float64,
		requires_grad = True
	)
	x = torch.tensor(x, dtype = torch.float64)
	x0 = torch.tensor(x0, dtype = torch.float64)
	NX = torch.tensor(NX, dtype = torch.float64)
	NXm = torch.tensor(NXm, dtype = torch.float64)
	
	# Record loss history
	hist = []
	# Track last good value in case of error
	eta_good = eta.detach().clone()
	# eta history
	eta_hist = []
	
	'''
	Exit codes:
	  0 = converged
	  1 = reached max iterations
	  2 = error occurred
	'''
	exit_code = 1

	optimizer = torch.optim.LBFGS(
	    [eta],
	    line_search_fn = "strong_wolfe",
	    max_iter = max_itr,
	    history_size = history_size,
	    tolerance_grad = tolerance_grad,
	    tolerance_change = tolerance_change
	)
	
	def closure():
		optimizer.zero_grad()
		loss = -1.0 * _local_loglik_count(eta, x, x0, NXm, NX, h)
		loss.backward()
		return loss

	for _ in range(max_epoch):
		loss = optimizer.step(closure)
		with torch.no_grad():
			hist.append(loss.detach().item())
			eta_hist.append(eta.detach().clone().unsqueeze(1))
			
	par_opt = eta.detach().numpy()
	eta_hist = torch.cat(eta_hist, dim = 1).numpy()
	
	# Estimate AIC
	if i is not None:
	    dev, df = aic_count(eta_good, x, NXm, NX, i, h)
	else:
	    dev = None
	    df = None

	return {
        "par": par_opt,
        "loss_hist": hist,
        "convergence": int(exit_code),
        "eta_hist": eta_hist,
        "deviance": dev,
        "df": df
    }
	
###############################################################################
