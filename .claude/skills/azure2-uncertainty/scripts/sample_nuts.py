#!/usr/bin/env python3
"""NUTS on the posterior of an AZURE2 fit, with the analytic gradient, in coordinates whitened at the MAP.

The target is the program's: log p(x) = Objective.log_prob(x) / T = -objective(x) / (2T), -inf outside the fit's
limits of every parameter kind and where Park's overlap J of a level is negative (the engine sampler's rule); its
gradient is -grad objective / (2T) (Objective.value_and_grad).  The priors are those of the project
(<parameterSettings>, nuisance priors, a prior covariance): a direction neither data nor priors close shows as a
non-positive eigenvalue of the Hessian and stops the run (add a prior; do not sample an open direction).

Whitening: x = x_map + L z with L = V diag(lambda^-1/2), H/2 = V diag(lambda) V^T the Hessian of -log p at x_map
(T = 1): --hessian exact (central differences of the analytic gradient, 2n gradients; the metric to sample with) or
gauss-newton (2 J^T J of the rows, one Jacobian; cheap, for a first look).  In z the posterior is close to N(0, I)
where the Laplace approximation holds, and NUTS (nuts.py) needs few gradients per draw.

Reads:  the project; the MAP (or best fit) from --x-json ("x") or --sav, else the project's start.
Writes: --out NPZ: x (chains, draws, parameters), z, logp, divergent, depth, ngrad, names, x_map, L, laplace_sigma;
        a summary on standard output: per parameter mean, sd, the Laplace sd, R-hat and ESS; divergences and
        gradients per draw.

    python3 sample_nuts.py fit.azr --x-json fit.json --chains 4 --warmup 300 --draws 500 --out nuts.npz

Run chains as separate one-thread processes for speed (--chains 1 --seed k each, then merge the files); one
process runs its chains one after the other.
Exit status: 0 done and every R-hat < 1.01, 1 done with R-hat >= 1.01 or divergences above 1 %, 2 bad arguments,
3 engine error or an open direction at x_map.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nuts import run_chain          # noqa: E402
from mcmc_stats import split_rhat, ess   # noqa: E402


def hessian(f, x, kind="exact", rel_step=1e-3):
    """Hessian of the objective at x: central differences of the analytic gradient, steps rel_step times the
    Gauss-Newton error of each parameter; or 2 J^T J."""
    _, J = f.residuals_and_jacobian(x)
    H_gn = 2 * J.T @ J
    if kind == "gauss-newton":
        return H_gn
    sig = 1 / np.sqrt(np.clip(np.diag(H_gn) / 2, 1e-300, None))
    H = np.zeros((len(x), len(x)))
    for j in range(len(x)):
        h = rel_step * sig[j]
        e = np.zeros(len(x))
        e[j] = h
        H[:, j] = (f.grad(x + e) - f.grad(x - e)) / (2 * h)
    return 0.5 * (H + H.T)


def whitening(H, temperature=1.0, floor=1e-8):
    """(L, eigenvalues of H/2) with x = x_map + L z; refuses a non-positive eigenvalue (ValueError)."""
    lam, V = np.linalg.eigh(H / 2)
    if lam.min() <= floor * lam.max():
        raise ValueError(f"the Hessian at x_map has eigenvalue {lam.min():.3g} (largest {lam.max():.3g}): a "
                         "direction the data and priors do not close, or x_map is not a minimum")
    return V / np.sqrt(lam) * np.sqrt(temperature), lam


def target(f, x_map, L, temperature=1.0):
    """logp_grad(z) of the whitened posterior; -inf where the engine's sampler excludes x or the engine fails."""
    import pyazr
    lo = np.array([b[0] for b in f.bounds])
    hi = np.array([b[1] for b in f.bounds])

    def logp_grad(z):
        x = x_map + L @ z
        if np.any(x < lo) or np.any(x > hi):
            return -np.inf, np.zeros_like(z)
        try:
            if not np.isfinite(f.log_prob(x)):
                return -np.inf, np.zeros_like(z)
            v, g = f.value_and_grad(x)
        except pyazr.Error:
            return -np.inf, np.zeros_like(z)
        return -v / (2 * temperature), L.T @ (-g / (2 * temperature))
    return logp_grad


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("project")
    ap.add_argument("--x-json")
    ap.add_argument("--sav")
    ap.add_argument("--hessian", choices=("exact", "gauss-newton"), default="exact")
    ap.add_argument("--temperature", type=float, default=1.0, help="T: 1, or chi2/nu to inflate as a Birge ratio")
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--warmup", type=int, default=300)
    ap.add_argument("--draws", type=int, default=500)
    ap.add_argument("--max-depth", type=int, default=8)
    ap.add_argument("--metric", choices=("auto", "dense", "diag", "identity"), default="auto")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True, help="NPZ of the draws")
    a = ap.parse_args(argv)
    import pyazr
    try:
        s = pyazr.Session(a.project, parameter_file=a.sav)
        f = s.objective()
        x_map = f.x0
        if a.x_json:
            with open(a.x_json) as fh:
                x_map = np.asarray(json.load(fh)["x"], float)
        H = hessian(f, x_map, a.hessian)
        L, lam = whitening(H, a.temperature)
    except (pyazr.Error, ValueError) as e:
        print(f"sample_nuts: {e}", file=sys.stderr)
        return 3
    logp_grad = target(f, x_map, L, a.temperature)
    chains = []
    for c in range(a.chains):
        print(f"chain {c}", flush=True)
        chains.append(run_chain(logp_grad, len(x_map), a.warmup, a.draws, seed=[a.seed, c], max_depth=a.max_depth,
                                metric=a.metric, log=print))
    z = np.array([c["z"] for c in chains])
    x = x_map + z @ L.T
    div = np.array([c["divergent"] for c in chains])
    ngrad = np.array([c["ngrad"] for c in chains])
    lap = np.sqrt(np.diag(L @ L.T))
    rh = [split_rhat(x[:, :, j]) for j in range(x.shape[2])]
    es = [ess(x[:, :, j]) for j in range(x.shape[2])]
    for j, n in enumerate(f.names):
        print(f"{n:28s} mean {x[:, :, j].mean():12.6g}  sd {x[:, :, j].std():10.4g}  Laplace sd {lap[j]:10.4g}  "
              f"R-hat {rh[j]:6.3f}  ESS {es[j]:7.0f}")
    print(f"divergences {div.sum()} of {div.size}; gradients per draw {ngrad.mean():.1f}")
    np.savez(a.out, x=x, z=z, logp=np.array([c["logp"] for c in chains]), divergent=div, ngrad=ngrad,
             depth=np.array([c["depth"] for c in chains]), names=np.array(f.names), x_map=x_map, L=L,
             laplace_sigma=lap)
    s.close()
    return 0 if max(rh) < 1.01 and div.mean() <= 0.01 else 1


if __name__ == "__main__":
    sys.exit(main())
