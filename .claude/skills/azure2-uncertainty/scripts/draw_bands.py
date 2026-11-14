#!/usr/bin/env python3
"""Bands on the model of AZURE2 data sets from parameter draws: posterior samples (sample_nuts.py,
sample_emcee.py) or draws from N(x, C) with C = (J^T J)^-1 at the best fit; the 16, 50 and 84 % quantiles per point
against the linear band sqrt(g^T C g) (g by central differences of the model).

A linear band and a band from draws agree where the model is close to linear in the parameters over their
uncertainty; where they differ (near a limit, an interference minimum, a poorly constrained width) the draws are
the band to quote.  Rates and extrapolations are evaluated the same way on a session of the extrapolations
(azure2-calculate owns those); each draw is one evaluation.

Reads:  the project; --draws NPZ (key x, any leading shape, last axis the free parameters) or --laplace N with the
        best fit from --x-json ("x") or --sav.
Writes: --out NPZ: per data set its energies, the quantiles, the linear band; a short table on standard output.

    python3 draw_bands.py fit.azr --x-json fit.json --laplace 400 --dataset sigma_pg --out bands.npz
    python3 draw_bands.py fit.azr --x-json fit.json --draws nuts.npz --thin 5 --out bands.npz

Exit status: 0 done, 2 bad arguments, 3 engine error.
"""
import argparse
import json
import sys

import numpy as np


def models(s, x, names):
    ev = s.evaluate(x)
    return [ev.dataset(n).model for n in names]


def linear_band(s, f, x, names, C, rel_step=1e-3):
    """sqrt(diag(G C G^T)) per data set, G by central differences in steps of rel_step x sqrt(C_jj)."""
    sig = np.sqrt(np.diag(C))
    cols = []
    for j in range(len(x)):
        e = np.zeros(len(x))
        e[j] = rel_step * sig[j]
        up, dn = models(s, x + e, names), models(s, x - e, names)
        cols.append([(u - d) / (2 * e[j]) for u, d in zip(up, dn)])
    out = []
    for k in range(len(names)):
        G = np.array([c[k] for c in cols]).T          # (points, parameters)
        out.append(np.sqrt(np.einsum("ij,jk,ik->i", G, C, G)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("project")
    ap.add_argument("--x-json")
    ap.add_argument("--sav")
    ap.add_argument("--draws", help="NPZ with draws under 'x'")
    ap.add_argument("--thin", type=int, default=1)
    ap.add_argument("--laplace", type=int, default=0, help="number of draws from N(x, (J^T J)^-1)")
    ap.add_argument("--dataset", action="append", default=[], help="data set name (repeatable; default all)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    if bool(a.draws) == bool(a.laplace):
        ap.error("give one of --draws and --laplace")
    import pyazr
    try:
        s = pyazr.Session(a.project, parameter_file=a.sav)
        f = s.objective()
        x = f.x0
        if a.x_json:
            with open(a.x_json) as fh:
                x = np.asarray(json.load(fh)["x"], float)
        ev = s.evaluate(x)
        names = a.dataset or [d.name for d in ev.data]
        _, J = f.residuals_and_jacobian(x)
        C = np.linalg.inv(J.T @ J)
        if a.draws:
            X = np.load(a.draws)["x"].reshape(-1, len(x))[::a.thin]
        else:
            X = np.random.default_rng(a.seed).multivariate_normal(x, C, size=a.laplace)
        M = [[] for _ in names]
        for xi in X:
            try:
                for k, m in enumerate(models(s, xi, names)):
                    M[k].append(m)
            except pyazr.Error:
                continue
        lin = linear_band(s, f, x, names, C)
    except (pyazr.Error, KeyError, np.linalg.LinAlgError) as e:
        print(f"draw_bands: {e}", file=sys.stderr)
        return 3
    out = {}
    for k, n in enumerate(names):
        q = np.percentile(np.array(M[k]), [16, 50, 84], axis=0)
        d = ev.dataset(n)
        half = 0.5 * (q[2] - q[0])
        ratio = np.median(half / np.where(lin[k] > 0, lin[k], np.nan))
        print(f"{n:24s} {d.n:5d} points  median (q84 - q16)/2 over the linear band: {ratio:.3f}")
        out.update({f"{n}__energies": d.energies, f"{n}__q16": q[0], f"{n}__q50": q[1], f"{n}__q84": q[2],
                    f"{n}__best": d.model, f"{n}__linear": lin[k]})
    if a.out:
        np.savez(a.out, draws_used=len(M[0]), **out)
    s.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
