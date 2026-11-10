#!/usr/bin/env python3
"""Reference values for tests/pyazr/potential_record_test.py (pyazr_potential_record): the channel functions of a hybrid pair whose Coulomb term
is a uniformly charged sphere of radius R_C = 5.0 fm (the `coulomb_radius` field of a `potential` record), and of the
same pair with R_C = 3.5 fm (the value the old <potential> block implied), by an integration independent of the
engine's.

System: p + 12C of tests/hybrid_potential (pair 1): mu = 1.00728 x 12 / 13.00728 u, Z1 Z2 = 6, Woods-Saxon V0 = 20 MeV,
R = 3.6 fm, a = 0.6 fm, e^2 = alpha hbar c, matching radius 15 fm, channel radius 3.6 fm.  V_C(r) = Z1 Z2 e^2 /(2 R_C)
(3 - r^2/R_C^2) inside R_C, Z1 Z2 e^2 / r outside.  Method: scipy's solve_ivp (DOP853, rtol 1e-12) inward from 15 fm
with the outer solution and its derivative from mpmath (25 digits): E > 0 the Coulomb H+ = G + iF, E < 0 the
Whittaker function W_{-eta_b, l+1/2}(2 kappa r).  Not the engine's Numerov, grid or Coulomb routines.

Rows (CSV): R_C (fm), l, E (MeV), quantity, value
  P   the penetrability rho / |u(a)|^2, rho = k a, u the outgoing solution (E > 0)
  S   the shift function Re(a u'(a) / u(a)) (E > 0), a W~'(a) / W~(a) (E < 0)

Run: python3 make_coulomb_radius_reference.py > coulomb_radius.csv   (about half a minute)
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
A_C, R_M = 3.6, 15.0
V0, R_WS, A_WS = 20.0, 3.6, 0.6


def V(r, l, rc):
    e2 = ALPHA * HBARC
    vc = Z1Z2 * e2 / (2.0 * rc) * (3.0 - (r / rc) ** 2) if r < rc else Z1Z2 * e2 / r
    vn = -V0 / (1.0 + math.exp((r - R_WS) / A_WS))
    return vc + vn + HBARC ** 2 / (2.0 * MU) * l * (l + 1) / r ** 2


def outer(l, E):
    c = Z1Z2 * ALPHA * MU / HBARC
    if E > 0:
        k = math.sqrt(2.0 * MU * E) / HBARC
        eta = Z1Z2 * ALPHA * math.sqrt(MU / (2.0 * E))
        return lambda r: mp.coulombg(l, eta, k * r) + 1j * mp.coulombf(l, eta, k * r)
    kappa = math.sqrt(-2.0 * MU * E) / HBARC
    return lambda r: mp.whitw(-c / kappa, l + 0.5, 2.0 * kappa * r)


def at_channel_radius(l, E, rc):
    w = outer(l, E)
    u0 = complex(w(R_M))
    du0 = complex(mp.diff(w, R_M))

    def f(r, y):
        g = 2.0 * MU / HBARC ** 2 * (V(r, l, rc) - E)
        return [y[2], y[3], g * y[0], g * y[1]]

    sol = solve_ivp(f, (R_M, A_C), [u0.real, u0.imag, du0.real, du0.imag], method="DOP853", rtol=1e-12,
                    atol=1e-14 * abs(u0))
    y = sol.y[:, -1]
    return complex(y[0], y[1]), complex(y[2], y[3])


def main():
    print("R_C_fm,l,E_MeV,quantity,value")
    for rc in (5.0, 3.5):
        for l in (0, 1, 2):
            for E in (0.5, 1.0, 2.0, 3.0):
                u, du = at_channel_radius(l, E, rc)
                k = math.sqrt(2.0 * MU * E) / HBARC
                print(f"{rc},{l},{E},P,{k * A_C / abs(u) ** 2:.15e}")
                print(f"{rc},{l},{E},S,{(A_C * du / u).real:.15e}")
            E = -1.94351  # the 13N ground state of tests/hybrid_potential
            u, du = at_channel_radius(l, E, rc)
            print(f"{rc},{l},{E},S,{(A_C * du / u).real:.15e}")


if __name__ == "__main__":
    main()
