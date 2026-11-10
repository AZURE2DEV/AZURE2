#!/usr/bin/env python3
"""Reference values for physics_hybrid_potential (tests/hybrid_potential_test.cpp, parts 3-5): the channel functions of
the hybrid model below, at and above threshold, by an integration independent of the engine's.

System: p + 12C of tests/hybrid_potential (pair 1): mu = 1.00728 x 12 / 13.00728 u, Z1 Z2 = 6, Woods-Saxon V0 = 20 MeV,
R = 3.6 fm, a = 0.6 fm, Coulomb of a uniform sphere R_C = 3.5 fm with e^2 = alpha hbar c, matching radius 15 fm,
channel radius 3.6 fm.  Method: scipy's solve_ivp (DOP853, rtol 1e-12) inward from 15 fm with the outer solution and
its derivative from mpmath (25 digits): E > 0 the Coulomb H+ = G + iF; E < 0 the Whittaker function
W_{-eta_b, l+1/2}(2 kappa r) (unit ANC); E = 0 sqrt(r) K_{2l+1}(sqrt(8 c r)), c = Z1 Z2 alpha mu / hbar c.  Not the
engine's Numerov, grid, Whittaker or Bessel routines.

Rows (CSV): kind, l, E (MeV), value
  S       S(E) = a u'(a)/u(a) (E <= 0) or Re of rho H'/H (E > 0)
  logW    ln W~(a) for E < 0 (W~ = W beyond 15 fm): the ANC <-> gamma factor
  N       int_a^inf W~(r)^2 dr (fm), the bound-state integral of external capture (NF integral numerator), E < 0
  Sc      the undistorted S (V0 = 0) at the same E, for the V_N -> 0 check

Run: python3 make_hybrid_bound_reference.py > hybrid_bound.csv   (about a minute)
"""
import math

import mpmath as mp
from scipy.integrate import solve_ivp, quad

mp.mp.dps = 25
HBARC = 197.3269804
UCONV = 931.49410242
ALPHA = 1.0 / 137.035999084
MU = 1.00728 * 12.0 / 13.00728 * UCONV
Z1Z2 = 6
A_C, R_M, R_C = 3.6, 15.0, 3.5


def V(r, l, V0):
    e2 = ALPHA * HBARC
    vc = Z1Z2 * e2 / (2.0 * R_C) * (3.0 - (r / R_C) ** 2) if r < R_C else Z1Z2 * e2 / r
    vn = -V0 / (1.0 + math.exp((r - 3.6) / 0.6)) if V0 else 0.0
    return vc + vn + HBARC ** 2 / (2.0 * MU) * l * (l + 1) / r ** 2


def outer(l, E):
    """The outer solution as a function of r (mpmath), and its kind."""
    c = Z1Z2 * ALPHA * MU / HBARC
    if E > 0:
        k = math.sqrt(2.0 * MU * E) / HBARC
        eta = Z1Z2 * ALPHA * math.sqrt(MU / (2.0 * E))
        return lambda r: mp.coulombg(l, eta, k * r) + 1j * mp.coulombf(l, eta, k * r)
    if E < 0:
        kappa = math.sqrt(-2.0 * MU * E) / HBARC
        return lambda r: mp.whitw(-c / kappa, l + 0.5, 2.0 * kappa * r)
    return lambda r: mp.sqrt(r) * mp.besselk(2 * l + 1, mp.sqrt(8.0 * c * r))


def solve(l, E, V0, dense=False):
    w = outer(l, E)
    u0 = complex(w(R_M))
    du0 = complex(mp.diff(w, R_M))
    def f(r, y):
        g = 2.0 * MU / HBARC ** 2 * (V(r, l, V0) - E)
        return [y[2], y[3], g * y[0], g * y[1]]
    sol = solve_ivp(f, (R_M, A_C), [u0.real, u0.imag, du0.real, du0.imag], method="DOP853", rtol=1e-12,
                    atol=1e-14 * abs(u0), dense_output=dense)
    return sol, w


def S_of(l, E, V0):
    sol, _ = solve(l, E, V0)
    y = sol.y[:, -1]
    u, du = complex(y[0], y[1]), complex(y[2], y[3])
    return (A_C * du / u).real


def bound(l, E, V0):
    sol, w = solve(l, E, V0, dense=True)
    y = sol.y[:, -1]
    logW = math.log(abs(y[0]))
    inner = quad(lambda r: sol.sol(r)[0] ** 2, A_C, R_M, epsabs=0, epsrel=1e-12, limit=200)[0]
    tail = float(mp.quad(lambda r: w(r) ** 2, [R_M, 30, 60, mp.inf]))
    return logW, inner + tail


print("kind,l,E_MeV,value")
for l in (0, 1, 2):
    for E in (-2.0, -1.94351, -1.0, -0.5, -0.2, -0.05, -0.01, 0.0, 0.001, 0.01, 0.05, 0.1, 0.5):
        print(f"S,{l},{E},{S_of(l, E, 20.0):.15e}")
        print(f"Sc,{l},{E},{S_of(l, E, 0.0):.15e}")
    for E in (-1.94351, -0.5):
        logW, N = bound(l, E, 20.0)
        print(f"logW,{l},{E},{logW:.15e}")
        print(f"N,{l},{E},{N:.15e}")
