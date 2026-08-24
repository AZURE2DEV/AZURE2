"""Sample an AZURE2 R-matrix model with emcee."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")

import emcee
import numpy as np
from scipy import stats
from multiprocess import Pool

from pyazr import azure2

NPROCS, NSTEPS = 8, 5000

# One engine per worker process.  Building it at module scope is what makes
# that happen either way a Pool starts: under "spawn" each worker re-imports
# this module and constructs its own, under "fork" each inherits a copy.  The
# engine is not thread-safe, so process-level parallelism is the only kind.
azr = azure2("13N.azr")
best = np.asarray(azr.params_rwa, float)
ndim = len(best)
nwalkers = 2 * ndim

# one prior per free parameter, chosen from its kind
sys_err = azr.datasets.sys_errors(vary_only=True)
priors, k = [], 0
for p in azr.parameters:
    if p.fixed:
        continue
    if p.kind == "energy":
        priors.append(stats.norm(p.value, 0.05))
    elif p.kind == "norm":
        se = sys_err[k]; k += 1
        priors.append(stats.norm(1, se) if se > 0 else stats.uniform(0, 10))
    elif p.kind == "shift":
        priors.append(stats.norm(0, 0.01))
    else:                                          # reduced width
        priors.append(stats.uniform(-1e6, 2e6))

offset = sum(np.sum(np.log(2 * np.pi * e ** 2)) for e in azr.cross_err)


def log_prob(theta):
    lp = np.sum([pr.logpdf(t) for pr, t in zip(priors, theta)])
    if not np.isfinite(lp):
        return -np.inf
    chi2 = np.sum(azr.calculate_chi2_rwa(theta))
    return -0.5 * (chi2 + offset) + lp


if __name__ == "__main__":
    p0 = best + 1e-4 * np.abs(best) * np.random.randn(nwalkers, ndim)
    with Pool(NPROCS) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)
        sampler.run_mcmc(p0, NSTEPS, progress=True)
    np.save("chain.npy", sampler.get_chain())
