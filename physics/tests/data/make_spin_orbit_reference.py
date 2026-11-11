#!/usr/bin/env python3
"""Reference values for physics_spin_orbit (tests/spin_orbit_test.cpp, part 3) and tests/pyazr/spin_orbit_potential_test.py:
the channel functions of a hybrid pair with a spin-orbit term, by an integration independent of the engine's.

System: p + 12C of tests/hybrid_potential (pair 1): mu = 1.00728 x 12 / 13.00728 u, Z1 Z2 = 6, Woods-Saxon V0 = 20 MeV,
R = 3.6 fm, a = 0.6 fm, Coulomb of a uniform sphere R_C = 3.5 fm with e^2 = alpha hbar c, matching radius 15 fm,
channel radius 3.6 fm; spin-orbit term on the proton (i = 1/2), V_so = 6 MeV, R_so = 3.6 fm, a_so = 0.6 fm, in
FRESCO's convention (docs/upgrade/spin_orbit_potential.md Eq. (3)):

    v_so(r) <2 l.i>,  v_so = V_so (hbar/m_pi c)^2 (1/r) df/dr,  (hbar/m_pi c)^2 = 2.000 fm^2,
    <2 l.i> = j(j+1) - l(l+1) - 3/4  (l for j = l + 1/2, -(l + 1) for j = l - 1/2).

Method: scipy's solve_ivp (DOP853, rtol 1e-12) inward from 15 fm with the outer solution and its derivative from
mpmath (25 digits): E > 0 the Coulomb H+ = G + iF; E < 0 the Whittaker function W_{-eta_b, l+1/2}(2 kappa r) (unit
ANC).  <2 l.i> is computed here from j(j+1) - l(l+1) - i(i+1); df/dr from the closed form of the derivative of the
Fermi function.  Not the engine's Numerov, grid, Coulomb or Whittaker routines.

Rows (CSV): kind, l, 2j, E (MeV), value
  S      Re rho H'/H (E > 0) or a W~'/W~ (E < 0)
  P      Im rho H'/H (E > 0)
  phi    arg H(a) (rad; the engine's hard-sphere phase factor is exp(-i phi), pyazr's delta_hs = -phi)
  logW   ln |W~(a)| (E < 0): the ANC <-> gamma factor of a bound state

Run: python3 make_spin_orbit_reference.py > spin_orbit.csv   (about a minute)
"""
import math

import mpmath as mp
from scipy.integrate import solve_ivp

mp.mp.dps = 25
HBARC = 197.3269804
UCONV = 931.49410242
ALPHA = 1.0 / 137.035999084
MU = 1.00728 * 12.0 / 13.00728 * UCONV
Z1Z2 = 6
A_C, R_M, R_C = 3.6, 15.0, 3.5
V0, R, A = 20.0, 3.6, 0.6
VSO, RSO, ASO = 6.0, 3.6, 0.6
CONLS = 2.000


def two_l_dot_i(l, two_j, two_i=1):
    j, i = two_j / 2.0, two_i / 2.0
    return j * (j + 1) - l * (l + 1) - i * (i + 1)


def V(r, l, ldots):
    e2 = ALPHA * HBARC
    vc = Z1Z2 * e2 / (2.0 * R_C) * (3.0 - (r / R_C) ** 2) if r < R_C else Z1Z2 * e2 / r
    vn = -V0 / (1.0 + math.exp((r - R) / A))
    x = (r - RSO) / ASO
    dfdr = -math.exp(x) / (ASO * (1.0 + math.exp(x)) ** 2)
    vso = VSO * CONLS * dfdr / r * ldots
    return vc + vn + HBARC ** 2 / (2.0 * MU) * l * (l + 1) / r ** 2 + vso


def outer(l, E):
    c = Z1Z2 * ALPHA * MU / HBARC
    if E > 0:
        k = math.sqrt(2.0 * MU * E) / HBARC
        eta = Z1Z2 * ALPHA * math.sqrt(MU / (2.0 * E))
        return lambda r: mp.coulombg(l, eta, k * r) + 1j * mp.coulombf(l, eta, k * r)
    kappa = math.sqrt(-2.0 * MU * E) / HBARC
    return lambda r: mp.whitw(-c / kappa, l + 0.5, 2.0 * kappa * r)


def solve(l, two_j, E):
    w = outer(l, E)
    ldots = two_l_dot_i(l, two_j)
    u0 = complex(w(R_M))
    du0 = complex(mp.diff(w, R_M))
    def f(r, y):
        g = 2.0 * MU / HBARC ** 2 * (V(r, l, ldots) - E)
        return [y[2], y[3], g * y[0], g * y[1]]
    sol = solve_ivp(f, (R_M, A_C), [u0.real, u0.imag, du0.real, du0.imag], method="DOP853", rtol=1e-12,
                    atol=1e-14 * abs(u0))
    y = sol.y[:, -1]
    return complex(y[0], y[1]), complex(y[2], y[3])


print("kind,l,two_j,E_MeV,value")
for l in (0, 1, 2):
    for two_j in sorted({abs(2 * l - 1), 2 * l + 1}):
        for E in (0.5, 1.0, 2.0, 3.0):
            u, du = solve(l, two_j, E)
            L = A_C * du / u
            print(f"S,{l},{two_j},{E},{L.real:.15e}")
            print(f"P,{l},{two_j},{E},{L.imag:.15e}")
            print(f"phi,{l},{two_j},{E},{math.atan2(u.imag, u.real):.15e}")
        for E in (-1.94351, -0.5):
            u, du = solve(l, two_j, E)
            print(f"S,{l},{two_j},{E},{(A_C * du / u).real:.15e}")
            print(f"logW,{l},{two_j},{E},{math.log(abs(u.real)):.15e}")
