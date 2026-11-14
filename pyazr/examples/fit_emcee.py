"""Sample the posterior of a project with emcee on the engine's objective, ``pyazr.Objective.log_prob``.

``log_prob(x) = -0.5 * objective(x)``: the data chi-squared and every prior and constraint of the project, the
same terms the program's fit minimises and the program's own sampler takes (docs/upgrade/objective.md section 6),
without its normalisation constant; ``-inf`` where the program's sampler puts zero probability (outside
``f.bounds``, the fit's limits of a parameter of any kind; Park's J < 0).  Flat in x elsewhere.

A session cannot be pickled and is not shared between processes: with a pool, every worker opens its own session
(``_open`` as the pool's initializer) and the walkers' points travel as arrays.  One thread per worker
(``OMP_NUM_THREADS=1``) keeps processes x threads within the machine.

The script prints what the objective holds besides the data (its prior terms at the start and the parameters with
finite limits), samples, and prints the diagnostics (acceptance fraction; the integrated autocorrelation time tau
of each parameter beside the chain's length, which should exceed 50 tau) and the posterior mean and standard
deviation with units.  ``--rate ENTRANCE EXIT`` adds a quantity derived per sample: the reaction rate
N_A<sigma v> at ``--t9`` temperatures for up to ``--draws`` draws of the chain (``pyazr.rate.rate_band(...,
draws=)``), as the 16th, 50th and 84th percentiles.

Usage (from the repository's root):

    python3 pyazr/examples/fit_emcee.py tests/upgrade/T3/T3.azr --steps 200 [--walkers 20] [--processes 0] \\
        [--rate 2 1 --t9 0.5 1] [--out chain.npy]
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
    _objective = pyazr.Session(path).objective()


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
    _objective = pyazr.Session(path).objective()
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


def main():
    import argparse
    import emcee
    from pyazr import rate
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr)")
    ap.add_argument("--walkers", type=int, default=None, help="default: 2 x the free parameters + 2")
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--processes", type=int, default=0, help="worker processes (0: none)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--rate", type=int, nargs=2, metavar=("ENTRANCE", "EXIT"), default=None,
                    help="the rate of these pair keys per sample")
    ap.add_argument("--t9", type=float, nargs="+", default=[0.5, 1.0], help="temperatures for --rate (GK)")
    ap.add_argument("--draws", type=int, default=30, help="draws of the chain for --rate")
    ap.add_argument("--out", default=None, help="save the chain (steps x walkers x parameters) as .npy here")
    args = ap.parse_args()

    with pyazr.Session(args.project) as s:
        ev = s.evaluate()
        unit = {q.name: q.unit for q in s.parameters.free}
        f = s.objective()
        names, bounds = f.names, f.bounds
    priors = {k: v for k, v in ev.terms.items() if k != "data, normalizations and energy shifts" and v}
    print(f"objective at the start {ev.objective:.6g}: data {ev.chi2:.6g}, other terms {priors or 'none'}")
    limited = [(n, lo, hi) for n, (lo, hi) in zip(names, bounds) if np.isfinite(lo) or np.isfinite(hi)]
    print(f"{len(names)} free parameters, {len(limited)} with finite limits"
          + "".join(f"\n  {n:<24} [{lo:.4g}, {hi:.4g}] {unit[n]}" for n, lo, hi in limited))

    sampler, names = sample(args.project, args.walkers, args.steps, seed=args.seed, processes=args.processes)
    chain = sampler.get_chain()
    burn = args.steps // 2
    flat = sampler.get_chain(discard=burn, flat=True)
    tau = emcee.autocorr.integrated_time(chain[burn:], tol=0)
    print(f"\n{chain.shape[1]} walkers x {args.steps} steps; acceptance {np.mean(sampler.acceptance_fraction):.2f}; "
          f"after a burn-in of {burn} steps the chain is {args.steps - burn} steps long, "
          f"{(args.steps - burn) / tau.max():.1f} x the longest tau (50 x is the usual minimum)")
    print(f"{'parameter':<24} {'mean':>14} {'std':>10} {'unit':<8} {'tau (steps)':>11}")
    for name, column, t_ in zip(names, flat.T, tau):
        print(f"{name:<24} {np.mean(column):14.7g} {np.std(column):10.3g} {unit[name]:<8} {t_:11.1f}")

    if args.rate:
        with rate.RateSession(args.project, *args.rate) as r:
            band = rate.rate_band(r, draws=flat, names=names, temperatures=args.t9, max_draws=args.draws,
                                  seed=args.seed)
        print(f"\nN_A<sigma v> of pair {args.rate[0]} -> {args.rate[1]} ({band.units['rate']}), {band.how}:")
        for T, lo, med, hi in zip(band.T9, band.low, band.median, band.high):
            print(f"  T9 = {T:g} GK: {med:.5g} (16 %: {lo:.5g}, 84 %: {hi:.5g})")
    if args.out:
        np.save(args.out, chain)
        print(f"wrote {args.out}")
    return sampler


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    main()
