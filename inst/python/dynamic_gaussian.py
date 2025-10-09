import torch
import numpy as np
from torch.func import hessian
from distributions import _local_loglik, _log_mvn_density

from concurrent.futures import ProcessPoolExecutor
import os
import multiprocessing as mp
import math

###############################################################################

def fit_gaussian(par0, x, NX, h, control, x0):
	max_epoch = int(control["max_epoch"])
	max_iter = int(control["max_itr"])
	history_size = int(control["history_size"])
	tolerance_grad = float(control["tolerance_grad"])
	tolerance_change = float(control["tolerance_change"])
	
	eta = torch.tensor(
		np.atleast_1d(par0).tolist(),
		dtype = torch.float64,
		requires_grad = True
	)
	x = torch.tensor(x, dtype = torch.float64)
	x0 = torch.tensor(x0, dtype = torch.float64)
	NX = torch.tensor(NX, dtype = torch.float64)
	
	optimizer = torch.optim.LBFGS(
		[eta],
		line_search_fn = "strong_wolfe",
		max_iter = max_iter,
		history_size = history_size,
		tolerance_grad = tolerance_grad,
		tolerance_change = tolerance_change
	)
	
	def closure():
		optimizer.zero_grad()
		loss = -1.0 * _local_loglik(eta, x, x0, NX, h)
		loss.backward()
		return loss

	for _ in range(max_epoch):
		loss = optimizer.step(closure)
		
	return eta.detach().numpy()

###############################################################################

def model_aic(eta_i, x, NX, h, i):
	eta_i = torch.tensor(
	    eta_i,
		dtype = torch.float64,
		requires_grad = True
	)
	x = torch.tensor(x, dtype = torch.float64)
	NX = torch.tensor(NX, dtype = torch.float64)
	i = int(i)

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
	J = torch.func.hessian(L_local)(eta_i)
	H_i = torch.func.hessian(L_i)(eta_i)
	nu = W0 * torch.trace(torch.linalg.solve(J, H_i))

	aic = -2 * dev + 2 * nu
	return aic.item()

###############################################################################
