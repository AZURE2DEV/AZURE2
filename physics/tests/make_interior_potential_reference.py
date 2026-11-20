#!/usr/bin/env python3
"""Reference values for physics/tests/interior_potential_test.cpp (docs/upgrade/interior_potential.md), independent of
the engine: scipy's DOP853 with mpmath's Coulomb functions.

  python3 physics/tests/make_interior_potential_reference.py > physics/tests/data/interior_potential.csv

alpha + 12C (masses 4.00260325413 and 12 u), the l = 1 Woods-Saxon of Koonin, Tombrello and Fox (Nucl. Phys. A 220
(1974) 221): V0 = 45.891 MeV, R = 3.55 fm, a = 0.60 fm, R_C = 3.55 fm, used here for l = 0, 1, 2.  Rows:
  f      a u'(a)/u(a) of the regular solution at the channel radius a_fm
  norm   int_0^a u^2 dr / u(a)^2 (fm)
  delta  the nuclear phase shift (rad, in (-pi/2, pi/2]) of the potential, the regular solution matched to the Coulomb
         functions at 15 fm: independent of any channel radius (a_fm is 0 in these rows)
  pole   the energy (MeV) of the pole of R_pot = 1/(f - b) with b = 0 at a_fm, l = 1, and in the next row its reduced
         width gamma^2 = -1/(df/dE) (MeV)
"""
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

HBARC = 197.3269804; ALPHA = 1 / 137.035999084; U = 931.49410242
MU = 4.00260325413 * 12.0 / 16.00260325413 * U
Z = 12; E2 = ALPHA * HBARC
V0, R, A, RC = 45.891, 3.55, 0.60, 3.55
K = 2 * MU / HBARC**2


def V(r, l):
    vc = Z * E2 / r if r >= RC else Z * E2 / (2 * RC) * (3 - (r / RC) ** 2)
    return vc - V0 / (1 + np.exp((r - R) / A)) + l * (l + 1) / (K * r * r)


def regular(E, a, l):
    r0 = 1e-3
    s = solve_ivp(lambda r, y: [y[1], K * (V(r, l) - E) * y[0]], (r0, a), [r0 ** (l + 1), (l + 1) * r0**l],
                  method="DOP853", rtol=1e-13, atol=1e-300, dense_output=True)
    return s


def f_reg(E, a, l):
    s = regular(E, a, l)
    return a * s.y[1, -1] / s.y[0, -1]


def norm(E, a, l):
    s = regular(E, a, l)
    val = quad(lambda r: s.sol(r)[0] ** 2, 1e-3, a, epsabs=0, epsrel=1e-11, limit=400)[0]
    return val / s.y[0, -1] ** 2


def delta(E, l, rm=15.0):
    f = f_reg(E, rm, l)
    k = np.sqrt(K * E); eta = Z * E2 * MU / (HBARC**2 * k); rho = k * rm
    F = mp.coulombf(l, eta, rho); G = mp.coulombg(l, eta, rho)
    dF = mp.diff(lambda x: mp.coulombf(l, eta, x), rho); dG = mp.diff(lambda x: mp.coulombg(l, eta, x), rho)
    d = float(mp.atan2(f * F - rho * dF, rho * dG - f * G))
    return (d + np.pi / 2) % np.pi - np.pi / 2


print("kind,l,a_fm,E_MeV,value")
for l in (0, 1, 2):
    for a in (4.0, 5.35, 6.5):
        for E in (-1.0, -0.2, 0.3, 1.0, 2.0, 3.5, 5.0):
            print(f"f,{l},{a},{E},{f_reg(E, a, l):.15e}")
            print(f"norm,{l},{a},{E},{norm(E, a, l):.15e}")
    for E in (0.3, 1.0, 2.0, 2.4, 2.6, 3.5, 5.0):
        print(f"delta,{l},0,{E},{delta(E, l):.15e}")
a = 5.35
E0 = brentq(lambda E: f_reg(E, a, 1), 1.0, 3.0, xtol=1e-12)
h = 1e-4
print(f"pole,1,{a},0,{E0:.12e}")
print(f"gamma2,1,{a},0,{-1 / ((f_reg(E0 + h, a, 1) - f_reg(E0 - h, a, 1)) / (2 * h)):.12e}")
