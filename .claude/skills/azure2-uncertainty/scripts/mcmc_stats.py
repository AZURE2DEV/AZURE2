#!/usr/bin/env python3
"""Convergence numbers of Markov chains: split R-hat, effective sample size, integrated autocorrelation time.

- ``split_rhat(x)``: x of shape (chains, draws); each chain is cut in two halves and the potential scale reduction
  of Gelman et al. (BDA3, eq. 11.4) computed over the 2m half chains.  Quote a posterior only below 1.01.
- ``ess(x)``: the effective sample size of x (chains, draws) by Geyer's initial monotone sequence of the
  autocorrelation averaged over the chains.
- ``tau(x)``: the integrated autocorrelation time draws / ESS per chain; emcee's rule: a chain shorter than about
  50 tau is not quoted.
Before R-hat on an R-matrix posterior, fold the exact sign symmetry of each level (all widths of one level change
sign together without changing any observable); a chain stuck in the mirror copy otherwise looks unconverged.

Reads:  --npz FILE with an array of draws (key --key, default "x"; shape (chains, draws, parameters)).
Writes: a table of R-hat, ESS and tau per parameter on standard output.

    python3 mcmc_stats.py --npz nuts.npz --key x

Exit status: 0 every R-hat < --rhat-max, 1 not, 2 bad arguments.
"""
import argparse
import sys

import numpy as np


def ess(x):
    x = np.atleast_2d(np.asarray(x, float))
    m, n = x.shape
    xc = x - x.mean()
    f = np.fft.rfft(xc, 2 * n, axis=1)
    ac = np.fft.irfft(f * np.conj(f), axis=1)[:, :n].mean(0) / n
    if ac[0] <= 0:
        return float(m * n)
    rho = ac / ac[0]
    s, prev = 0.0, np.inf
    for k in range(0, n - 1, 2):
        pair = rho[k] + rho[k + 1]
        if pair < 0:
            break
        pair = min(pair, prev)
        s += pair
        prev = pair
    t = max(2 * s - 1, 1e-3)
    return m * n / t


def tau(x):
    x = np.atleast_2d(np.asarray(x, float))
    return x.shape[1] * x.shape[0] / ess(x)


def split_rhat(x):
    x = np.atleast_2d(np.asarray(x, float))
    n = x.shape[1] // 2
    halves = np.concatenate([x[:, :n], x[:, n:2 * n]], axis=0)
    W = halves.var(axis=1, ddof=1).mean()
    B = n * halves.mean(axis=1).var(ddof=1)
    if W <= 0:
        return 1.0
    return float(np.sqrt(((n - 1) / n * W + B / n) / W))


def table(draws, names=None, out=print):
    """draws (chains, draws, parameters) -> list of (name, rhat, ess, tau)."""
    rows = []
    for j in range(draws.shape[2]):
        name = names[j] if names is not None else f"p{j}"
        rows.append((name, split_rhat(draws[:, :, j]), ess(draws[:, :, j]), tau(draws[:, :, j])))
        out(f"{name:28s} R-hat {rows[-1][1]:7.4f}  ESS {rows[-1][2]:9.1f}  tau {rows[-1][3]:7.2f}")
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--npz", required=True)
    ap.add_argument("--key", default="x")
    ap.add_argument("--rhat-max", type=float, default=1.01)
    a = ap.parse_args(argv)
    d = np.load(a.npz, allow_pickle=False)
    if a.key not in d.files or d[a.key].ndim != 3:
        print(f"mcmc_stats: {a.npz} has no 3-d array {a.key!r} (keys {d.files})", file=sys.stderr)
        return 2
    names = list(d["names"]) if "names" in d.files else None
    rows = table(d[a.key], names)
    return 0 if all(r[1] < a.rhat_max for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
