#!/usr/bin/env python3
"""Independent reference for tests/q_attenuation (usage: verify.py WORK; exit status = number of failures).

(b) q.azr (40 coefficients of a 15-degree cone) against a brute-force detector average: the c.m. distribution W*(x*)
    of plain.azr at each energy is fitted with Legendre polynomials (degree 20), carried to the lab
    (W_lab = W* dOmega*/dOmega_lab, cos theta_lab = (x* + g)/sqrt(1 + g^2 + 2 g x*)), averaged over a uniform cone of
    half-angle 15 degrees around the lab direction of each c.m. angle by 2-d quadrature, and taken back with
    dOmega_lab/dOmega* -- no Legendre transform of the attenuation coefficients is used.
(c) q5.azr (Q_0..Q_4 of that cone, Q_K = 1 above) against W* + J sum_{K<5} (Q_K - 1) b_K P_K(x_lab), with the lab
    Legendre coefficients b_K projected in the LAB variable (the engine integrates in the c.m. one).
"""
import os
import sys

import numpy as np
from numpy.polynomial import legendre as L

WORK = sys.argv[1]
M1, M2, M3, M4 = 4.0026, 13.0033, 1.00866, 15.9949          # alpha + 13C -> n + 16O(6.13), pairs 2 and 5
QV = 6.35869 - 4.14313 - 6.12989
HALF = np.radians(15.0)
OUT = "output/AZUREOut_aa=2_R=5.extrap"
failures = 0


def load(name):
    a = np.loadtxt(os.path.join(WORK, name, OUT))
    return {round(e, 8): a[np.isclose(a[:, 0], e)][:, [2, 3]] for e in np.unique(a[:, 0])}


def P(k, x):
    return L.legval(x, np.eye(k + 1)[k])


def cone_q(n):
    u, w = L.leggauss(200)
    c = 0.5 * (1 - np.cos(HALF)) * u + 0.5 * (1 + np.cos(HALF))
    return np.array([np.sum(w * P(k, c)) / np.sum(w) for k in range(n)])


def report(ok, text):
    global failures
    print(("  ok    " if ok else "  FAIL  ") + text)
    failures += 0 if ok else 1


plain, q, q5 = load("plain"), load("q"), load("q5")
Q5 = cone_q(5)
u, wu = L.leggauss(160)
cp = 0.5 * (1 - np.cos(HALF)) * u + 0.5 * (1 + np.cos(HALF)); wcp = 0.5 * (1 - np.cos(HALF)) * wu
phi = (np.arange(320) + 0.5) * 2 * np.pi / 320
sp = np.sqrt(1 - cp**2)
gx, gw = L.leggauss(200)
errb, errc, err_naive = 0.0, 0.0, 0.0
for e, d in plain.items():
    ang, w = d[:, 0], d[:, 1]
    x = np.cos(np.radians(ang))
    c = L.legfit(x, w, 20)
    g = np.sqrt(M1 * M3 * e / (M2 * M4 * (e + QV)))
    lab = lambda xc: (xc + g) / np.sqrt(1 + g * g + 2 * g * xc)
    jac = lambda xc: (1 + g * xc) / (1 + g * g + 2 * g * xc) ** 1.5      # dOmega_lab / dOmega*
    def cm_of_lab(xl):
        th = np.arccos(np.clip(xl, -1, 1))
        return np.cos(th + np.arcsin(g * np.sin(th)))
    w_lab = lambda xl: L.legval(cm_of_lab(xl), c) / jac(cm_of_lab(xl))
    # (b) brute-force cone average around each lab direction
    for xc, val in zip(x, q[e][:, 1]):
        tl = np.arccos(np.clip(lab(xc), -1.0, 1.0))
        xl = np.cos(tl) * cp[:, None] - np.sin(tl) * sp[:, None] * np.cos(phi)[None, :]
        ref = np.sum(wcp[:, None] * w_lab(xl)) / (len(phi) * np.sum(wcp)) * jac(xc)
        errb = max(errb, abs(val / ref - 1)) if np.isfinite(ref) else np.inf
        naive = np.sum(c[:5] * cone_q(5) * np.array([P(k, xc) for k in range(5)])) + L.legval(xc, np.r_[0 * c[:5], c[5:]])
        err_naive = max(err_naive, abs(naive / ref - 1))
    # (c) lab-variable projection
    b = np.array([(2 * k + 1) / 2 * np.sum(gw * P(k, gx) * w_lab(gx)) for k in range(5)])
    for xc, val in zip(x, q5[e][:, 1]):
        ref = L.legval(xc, c) + jac(xc) * np.sum((Q5 - 1) * b * np.array([P(k, lab(xc)) for k in range(5)]))
        errc = max(errc, abs(val / ref - 1)) if np.isfinite(ref) else np.inf
report(errb < 2e-6, f"(b) 40 coefficients vs the brute-force cone average: max rel. deviation {errb:.1e} "
                    f"(Q applied in the c.m. instead: {err_naive:.1e})")
report(errc < 2e-6, f"(c) 5 coefficients vs the lab-variable projection: max rel. deviation {errc:.1e}")
sys.exit(failures)
