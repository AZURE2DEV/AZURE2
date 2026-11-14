#!/usr/bin/env python3
"""emcee's affine-invariant ensemble on the posterior of an AZURE2 fit (Objective.log_prob: -objective/2, -inf
outside the fit's limits and where Park's J < 0, as the program's own sampler).

The walkers start as draws from the Laplace approximation N(x_map, (J^T J)^-1) moved inside the limits.  The
ensemble's own autocorrelation time is reported, with its rule: a run shorter than 50 tau is not quoted.  In many
dimensions an ensemble sampler can report a small tau while it has not equilibrated; compare its widths with a
NUTS run (sample_nuts.py) before quoting it.

Reads:  the project; the start point from --x-json ("x") or --sav, else the project's start.
Writes: --out NPZ: x (walkers, steps after burn-in, parameters; the layout mcmc_stats.py reads as chains, draws),
        logp, names, acceptance, tau; a summary on standard output.

    python3 sample_emcee.py fit.azr --x-json fit.json --walkers 32 --steps 2000 --burn 500 --out emcee.npz

Exit status: 0 done and steps after burn-in >= 50 tau, 1 shorter, 2 bad arguments (emcee missing), 3 engine error.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mcmc_stats import tau   # noqa: E402


def laplace_starts(f, x, nwalkers, rng, spread=1.0):
    """Walkers drawn from N(x, spread^2 (J^T J)^-1), clipped inside the fit's limits."""
    _, J = f.residuals_and_jacobian(x)
    C = np.linalg.pinv(J.T @ J)
    w, V = np.linalg.eigh(0.5 * (C + C.T))
    Lc = V * np.sqrt(np.clip(w, 0, None))
    p = x + spread * rng.standard_normal((nwalkers, len(x))) @ Lc.T
    lo = np.array([b[0] for b in f.bounds])
    hi = np.array([b[1] for b in f.bounds])
    return np.clip(p, np.where(np.isfinite(lo), lo + 1e-12 * np.abs(lo) + 1e-300, lo),
                   np.where(np.isfinite(hi), hi - 1e-12 * np.abs(hi) - 1e-300, hi))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("project")
    ap.add_argument("--x-json")
    ap.add_argument("--sav")
    ap.add_argument("--walkers", type=int, default=0, help="default max(16, 2 n + 2)")
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--burn", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    try:
        import emcee
    except ImportError:
        print("sample_emcee: emcee is not installed", file=sys.stderr)
        return 2
    import pyazr
    try:
        s = pyazr.Session(a.project)
        if a.sav:
            s.load_sav(a.sav)
        f = s.objective()
        x = f.x0
        if a.x_json:
            with open(a.x_json) as fh:
                x = np.asarray(json.load(fh)["x"], float)
        n = len(x)
        nw = a.walkers or max(16, 2 * n + 2)
        rng = np.random.default_rng(a.seed)
        p0 = laplace_starts(f, x, nw, rng)
        sampler = emcee.EnsembleSampler(nw, n, f.log_prob)
        sampler.random_state = np.random.RandomState(a.seed).get_state()   # emcee's own draws, fixed by the seed
        sampler.run_mcmc(p0, a.steps, progress=False)
    except pyazr.Error as e:
        print(f"sample_emcee: {e}", file=sys.stderr)
        return 3
    chain = sampler.get_chain()[a.burn:]               # (steps, walkers, n)
    xs = np.transpose(chain, (1, 0, 2))                # (walkers, steps, n)
    taus = np.array([tau(xs[:, :, j]) for j in range(n)])
    acc = float(np.mean(sampler.acceptance_fraction))
    for j, name in enumerate(f.names):
        print(f"{name:28s} mean {xs[:, :, j].mean():12.6g}  sd {xs[:, :, j].std():10.4g}  tau {taus[j]:7.1f}")
    kept = xs.shape[1]
    print(f"{nw} walkers, {kept} steps kept, acceptance {acc:.3f}, largest tau {taus.max():.1f} "
          f"({kept / taus.max():.0f} tau kept; quote only above 50)")
    np.savez(a.out, x=xs, logp=sampler.get_log_prob()[a.burn:].T, names=np.array(f.names), acceptance=acc,
             tau=taus, x_start=x)
    s.close()
    return 0 if kept >= 50 * taus.max() else 1


if __name__ == "__main__":
    sys.exit(main())
