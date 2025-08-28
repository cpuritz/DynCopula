import torch
import numpy as np
from distributions import _local_loglik
from aic import compute_aic

###############################################################################

def fit_gaussian(par0, x, NX, h, control, x0 = None, i = None):
	max_itr = int(control["max_itr"])
	patience = int(control["patience"])
	reltol = float(control["reltol"])
	max_grad = float(control["max_grad"])
	
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

    # MLE using SGD with momentum
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
		# Minimize negative log likelihood
		loss = -1.0 *  _local_loglik(eta, x, x0, NX, h)
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
			# Check for convergence
			if (all(r >= 0 and r <= reltol for r in rel)):
			    if len(hist) == patience + 1:
			        # No improvement has been made. Try increasing learning rate
			        # by a factor of 10 before stopping.
			        lr = optimizer.param_groups[0]["lr"]
			        optimizer.param_groups[0]["lr"] = 10 * lr
			    else:
				    exit_code = 0
				    break

		loss.backward()
		torch.nn.utils.clip_grad_norm_(eta, max_norm = max_grad)
		optimizer.step()
		
		# If an NaN's appeared, stop and return the previous value of eta
		if torch.isnan(eta).any():
			exit_code = 2
			break
		eta_good = eta.clone()

	par_opt = eta_good.detach().numpy()
	eta_hist = torch.cat(eta_hist, dim = 1).numpy()
	
	# Compute AIC if time points to perform inference at were also time points
	# the time series was sampled at
	if i is not None:
	    aic = compute_aic(eta_good, x, NX, i, h)
	else:
	    aic = None

	return {
        "par": par_opt,
        "loss_hist": hist,
        "convergence": int(exit_code),
        "eta_hist": eta_hist,
        "aic": aic
    }

###############################################################################
