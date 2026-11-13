#!/usr/bin/env python3
"""The five Legendre expansions of quantity legendre (docs/upgrade/data_kinds.md section 4.4) from the engine's own
dsigma/dOmega (pyazr, full precision; the projection is numpy's, no AZURE2 code).

With dsigma/dOmega_cm(theta) = sum_L A_L P_L(cos theta), A_L = (2L+1)/2 int_{-1}^{1} dsigma/dOmega P_L(x) dx, done by
24-point Gauss-Legendre quadrature on the extrapolation dsdo_nodes of models/extrap.azr (its angles are the nodes;
exact for a polynomial of degree <= 47 in cos theta, and the model's is of degree 2 l_max = 8).  The expansions:
a_over_a0 A_L/A_0, absolute A_L, rsl A_L/((2L+1) A_0), 2l2 2 A_L/(2L+1), 4pi 4 pi A_L.

  truth.py cli <output dir of a menu-3 run of models/extrap.azr>
      (a) below on the files the CLI writes (AZUREOut_aa=2_R=5.acoeff: the five expansions in record order; .extrap:
          dsigma/dOmega at the nodes; 11 significant digits): exit 1 beyond 1e-8.
  truth.py check <dir with m.azr = models/extrap.azr>  (pyazr, full precision)
      (a) every order of every expansion against the projection, relative to |truth| (relative to the expansion's
          order-0 value where the truth is below 1e-8 of it: the model's odd orders and L >= 4 vanish, and the
          quadrature leaves 1e-14 of A_0 there): exit 1 beyond 1e-8;
      (b) the five expansions of one run against each other by their defining relations (absolute = a_over_a0 x A_0,
          rsl = a_over_a0/(2L+1), 2l2 = 2 absolute/(2L+1), 4pi = 4 pi absolute): exit 1 beyond 1e-12.
  truth.py data <dir>
      writes models/data/<expansion>.dat: rows E_lab, L, value, error at the three energies, L = 0 ... 6, the truth
      times 1 + 0.02 sin(1.3 i) with 3 % errors (2e-3 of the order-0 value where the coefficient vanishes).
"""
import math
import os
import sys

import numpy as np

EXPANSIONS = ["a_over_a0", "absolute", "rsl", "2l2", "4pi"]
LAB = [5.9, 6.1, 6.3]


def engine(workdir):
    from pyazr import Session
    with Session(os.path.join(workdir, "m.azr"), cwd=workdir, data_mode=False) as m:
        grids = m.evaluate().data
    return ({e: np.array(grids[i].extra["legendre"]) for i, e in enumerate(EXPANSIONS)}, grids[5].model,
            grids[5].angles)


def projection(dsdo, angles, kmax):
    x, w = np.polynomial.legendre.leggauss(24)
    out = []
    for i in range(len(LAB)):
        f, th = dsdo[i::len(LAB)], angles[i::len(LAB)]  # the grid's points: angle outer, energy inner
        if np.max(np.abs(np.cos(np.radians(th)) - x)) > 1e-10:
            raise SystemExit("the dsdo angles are not the Gauss-Legendre nodes")
        A = [(2 * L + 1) / 2 * np.sum(w * f * np.polynomial.legendre.Legendre.basis(L)(x)) for L in range(kmax + 1)]
        out.append(np.array(A))
    return np.array(out)


def forms(A):
    L = np.arange(A.shape[1])
    a0 = A[:, :1]
    return {"a_over_a0": A / a0, "absolute": A, "rsl": A / ((2 * L + 1) * a0), "2l2": 2 * A / (2 * L + 1),
            "4pi": 4 * math.pi * A}


def cli(outdir):
    blocks, cur = [], []
    for line in open(os.path.join(outdir, "AZUREOut_aa=2_R=5.acoeff")):
        if line.strip():
            cur.append([float(v) for v in line.split()[1:]])
        elif cur:
            blocks.append(np.array(cur))
            cur = []
    if cur:
        blocks.append(np.array(cur))
    rows = [[float(v) for v in line.split()] for line in open(os.path.join(outdir, "AZUREOut_aa=2_R=5.extrap"))
            if line.strip()]
    rows.sort(key=lambda r: (-r[2], r[0]))  # angle outer (the nodes from 180 deg down), energy inner, as pyazr
    if len(blocks) != len(EXPANSIONS) or len(rows) != 24 * len(LAB):
        print("expected 5 coefficient blocks and 72 dsdo rows, got", len(blocks), len(rows))
        return 1
    coef = dict(zip(EXPANSIONS, blocks))
    worst = against_projection(coef, np.array([r[3] for r in rows]), np.array([r[2] for r in rows]))
    print("CLI output: every order of every expansion against the projection: %.1e" % worst)
    return 0 if worst <= 1e-8 else 1


def against_projection(coef, dsdo, angles):
    truth = forms(projection(dsdo, angles, coef["a_over_a0"].shape[1] - 1))
    worst = 0.0
    for e in EXPANSIONS:
        got, want = coef[e], truth[e]
        scale = np.where(np.abs(want) > 1e-8 * np.abs(want[:, :1]), np.abs(want), np.abs(want[:, :1]))
        worst = max(worst, float(np.max(np.abs(got - want) / scale)))
    return worst


def check(workdir):
    coef, dsdo, angles = engine(workdir)
    worst_t, worst_r = against_projection(coef, dsdo, angles), 0.0
    c = coef
    L = np.arange(c["a_over_a0"].shape[1])
    pairs = [(c["absolute"], c["a_over_a0"] * c["absolute"][:, :1]), (c["rsl"], c["a_over_a0"] / (2 * L + 1)),
             (c["2l2"], 2 * c["absolute"] / (2 * L + 1)), (c["4pi"], 4 * math.pi * c["absolute"])]
    for got, want in pairs:
        scale = np.maximum(np.abs(want), 1e-300)
        nz = np.abs(want) > 0
        worst_r = max(worst_r, float(np.max(np.abs(got - want)[nz] / scale[nz])) if nz.any() else 0.0)
        if np.any(got[~nz] != 0):
            worst_r = 1.0
    print("every order of every expansion against the projection: %.1e; the expansions against each other: %.1e"
          % (worst_t, worst_r))
    return 0 if worst_t <= 1e-8 and worst_r <= 1e-12 else 1


def data(workdir):
    coef, dsdo, angles = engine(workdir)
    truth = forms(projection(dsdo, angles, 6))
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models", "data")
    for e in EXPANSIONS:
        with open(os.path.join(here, e + ".dat"), "w") as f:
            f.write("# E_lab [MeV]  L  coefficient (%s)  error  (truth.py data)\n" % e)
            n = 0
            for i, E in enumerate(LAB):
                for L in range(7):
                    t = truth[e][i, L]
                    v = t * (1 + 0.02 * math.sin(1.3 * n))
                    err = abs(t) * 0.03 if abs(t) > 1e-12 * abs(truth[e][i, 0]) else 2e-3 * abs(truth[e][i, 0])
                    f.write("%r  %d  %.17e  %.17e\n" % (E, L, v, err))
                    n += 1
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("cli", "check", "data"):
        print(__doc__)
        sys.exit(2)
    sys.exit({"cli": cli, "check": check, "data": data}[sys.argv[1]](sys.argv[2]))
