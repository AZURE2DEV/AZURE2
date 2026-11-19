#!/usr/bin/env python3
"""The nuclear phase shifts of the Gaussian 3He + 4He potential VD (Tursunov, Turakulov and Kadyrov, Phys. Rev. C 97
(2018) 035802, Tables I and II: V0 exp(-alpha r^2), alpha = 0.15747 fm^-2, depths per (l, J) in MeV s1/2 67.5,
p1/2 81.815179, p3/2 83.589554, d3/2 66.0, d5/2 69.0, f5/2 75.9, f7/2 84.8; Coulomb of a uniform sphere R_C = 3.095
fm), independent of the engine: scipy's solve_ivp (DOP853, rtol 1e-12) outward from r = 1e-4 fm with u = r^(l+1),
matched at 15 fm to mpmath's F_l, G_l (30 digits).  Masses: the atomic masses of the models (3.0160293 u, 4.0026032
u), e^2 = alpha hbar c.  Writes data/phases.csv: l, 2J, E_cm (MeV), delta (deg, reduced to (-90, 90]).

    python3 make_phases.py > data/phases.csv      (a few seconds; needs numpy, scipy, mpmath)
"""
import math

import mpmath as mp
from scipy.integrate import solve_ivp

mp.mp.dps = 30
HBARC, UCONV, ALPHA = 197.3269804, 931.49410242, 1.0 / 137.035999084
M1, M2 = 3.0160293, 4.0026032
MU = M1 * M2 / (M1 + M2) * UCONV
RC, AL = 3.095, 0.15747
DEPTH = {(0, 1): 67.5, (1, 1): 81.815179, (1, 3): 83.589554, (2, 3): 66.0, (2, 5): 69.0, (3, 5): 75.9, (3, 7): 84.8}


def V(r, l, v0):
    e2 = ALPHA * HBARC
    vc = 4 * e2 / (2.0 * RC) * (3.0 - (r / RC) ** 2) if r < RC else 4 * e2 / r
    return vc - v0 * math.exp(-AL * r * r) + HBARC ** 2 / (2.0 * MU) * l * (l + 1) / r ** 2


def phase(l, v0, E):
    k2 = 2.0 * MU / HBARC ** 2
    r0, rm = 1e-4, 15.0
    sol = solve_ivp(lambda r, y: [y[1], k2 * (V(r, l, v0) - E) * y[0]], (r0, rm), [r0 ** (l + 1), (l + 1) * r0 ** l],
                    method="DOP853", rtol=1e-12, atol=1e-14 * r0 ** (l + 1))
    u, du = sol.y[0, -1], sol.y[1, -1]
    k = math.sqrt(2.0 * MU * E) / HBARC
    eta = 4 * ALPHA * math.sqrt(MU / (2.0 * E))
    F, G = mp.coulombf(l, eta, k * rm), mp.coulombg(l, eta, k * rm)
    dF = mp.diff(lambda x: mp.coulombf(l, eta, x), k * rm)
    dG = mp.diff(lambda x: mp.coulombg(l, eta, x), k * rm)
    d = math.atan2(float(u * k * dF - du * F), float(du * G - u * k * dG))
    d = math.degrees(d)
    while d > 90.0:
        d -= 180.0
    while d <= -90.0:
        d += 180.0
    return d


print("l,two_j,E_cm,delta_deg")
for (l, tj), v0 in DEPTH.items():
    for E in (0.5, 1.0, 2.0, 4.0):
        print(f"{l},{tj},{E},{phase(l, v0, E):.12e}")
