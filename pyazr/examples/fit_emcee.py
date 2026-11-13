"""Sample the posterior of a project with emcee on the engine's objective, ``pyazr.Objective.log_prob``.

``log_prob(x) = -0.5 * objective(x)``: the data chi-squared and every prior and constraint of the project, the
same terms the program's fit minimises and the program's own sampler takes (docs/upgrade/objective.md section 6),
without its normalisation constant; ``-inf`` where the program's sampler puts zero probability (outside the bounds
of a parameter whose record has ``sampler_support``, a ``gamma_decay`` record's numbers; Park's J < 0).  Flat in x
elsewhere.  For walls at every bound, wrap it:
``lambda x: f.log_prob(x) if all(lo <= v <= hi for v, (lo, hi) in zip(x, f.bounds)) else -numpy.inf``.

A session cannot be pickled and is not shared between processes: with a pool, every worker opens its own session
(``_open`` as the pool's initializer) and the walkers' points travel as arrays.  One thread per worker
(``OMP_NUM_THREADS=1``) keeps processes x threads within the machine.

    python fit_emcee.py path/to/project.azr [walkers steps processes]
"""
import os
import sys

import numpy as np

import pyazr

_objective = None  # this process's Objective (a pool worker's, or the main process's)


def _open(path, threads):
    """A pool worker's initializer: its own session and Objective."""
    global _objective
    os.environ["OMP_NUM_THREADS"] = str(threads)
    _objective = pyazr.Session(path).fit_objective()


def log_prob(x):
    """ln P of this process's Objective (a module-level function: a pool sends its name, not a session)."""
    return _objective.log_prob(x)


def start(f, nwalkers, spread, rng):
    """Walkers around x0 (relative spread), each redrawn until log_prob is finite (inside the support)."""
    x0 = f.x0
    width = spread * np.where(x0 != 0.0, np.abs(x0), 1.0)
    p0 = np.empty((nwalkers, len(x0)))
    for k in range(nwalkers):
        for _ in range(100):
            p0[k] = x0 + width * rng.standard_normal(len(x0))
            if np.isfinite(f.log_prob(p0[k])):
                break
        else:
            raise RuntimeError(f"no start inside the support near x0 for walker {k}; make the spread smaller")
    return p0


def sample(path, nwalkers=None, nsteps=100, spread=1e-3, seed=1, processes=0, threads=1):
    """A chain of ``nsteps`` steps; ``processes`` > 0 samples in a pool of that many workers (one session each).
    Returns the emcee sampler (``get_chain()``, ``get_log_prob()``) and the parameters' names."""
    import emcee
    global _objective
    _objective = pyazr.Session(path).fit_objective()
    f = _objective
    ndim = len(f.x0)
    nwalkers = nwalkers or 2 * ndim + 2
    rng = np.random.default_rng(seed)
    p0 = start(f, nwalkers, spread, rng)
    state = np.random.RandomState(seed).get_state()  # emcee's own draws, fixed by the seed
    if processes:
        import multiprocessing
        context = multiprocessing.get_context("spawn")  # a fresh interpreter per worker, no OpenMP state forked
        with context.Pool(processes, initializer=_open, initargs=(os.path.abspath(path), threads)) as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)
            sampler.random_state = state
            sampler.run_mcmc(p0, nsteps)
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, f.log_prob)
        sampler.random_state = state
        sampler.run_mcmc(p0, nsteps)
    return sampler, f.names


def main(path, nwalkers=None, nsteps=200, processes=2):
    sampler, names = sample(path, nwalkers, nsteps, processes=processes)
    chain = sampler.get_chain(discard=nsteps // 2, flat=True)
    for name, column in zip(names, chain.T):
        print(f"  {name:<24} {np.mean(column):14.7g} +- {np.std(column):.3g}")
    print(f"acceptance {np.mean(sampler.acceptance_fraction):.2f}")
    np.save("chain.npy", sampler.get_chain())
    return sampler


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args[0] if args else "T2.azr", *(int(a) for a in args[1:4]))
