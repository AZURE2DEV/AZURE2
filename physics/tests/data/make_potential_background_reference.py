#!/usr/bin/env python3
"""Reference values for physics_potential_background (tests/potential_background_test.cpp), independent of the engine's
Numerov integration, grid and Coulomb routines (docs/upgrade/potential_background.md section 7.2).

System: 3He + 4He, mu = m(3He) m(4He) / (m(3He) + m(4He)) with the atomic masses 3.0160293 u and 4.0026032 u (those of
tests/scattering_length/he3he4.azr), Z1 Z2 = 4, e^2 = alpha hbar c.

Rows (CSV): kind, l, a (fm), E (MeV), value
  ellC    a F_l'(a)/F_l(a), F_l(eta, k r) the regular point-Coulomb function (mpmath coulombf, 30 digits): the
          interior log-derivative of the pure point-Coulomb problem (no nuclear term)
  delta   the nuclear phase shift delta_l (rad, reduced to (-pi/2, pi/2]) of the Gaussian potential
          V_N = -V0_l exp(-(r/r0)^2), r0 = 1/sqrt(0.15747) fm, with the Coulomb potential of a uniform sphere of radius
          3.095 fm (Dubovichenko, arXiv:1006.4944, Table 3.1, 4He+3He: S 67.5, P3/2 83.83, D5/2 69.0, F7/2 84.8 MeV,
          one depth per l here), by scipy's solve_ivp (DOP853, rtol 1e-12, atol 1e-14) outward from r = 1e-4 fm with the
          series start u = r^(l+1), matched at 15 fm to mpmath's F_l, G_l: u ~ cos(delta) F + sin(delta) G.

Run: python3 make_potential_background_reference.py > potential_background.csv   (about a minute)
"""
import math

import mpmath as mp
from scipy.integrate import solve_ivp

mp.mp.dps = 30
HBARC = 197.3269804
UCONV = 931.49410242
ALPHA = 1.0 / 137.035999084
M1, M2 = 3.0160293, 4.0026032
MU = M1 * M2 / (M1 + M2) * UCONV
Z1Z2 = 4
RC = 3.095
R0 = 1.0 / math.sqrt(0.15747)
DEPTH = {0: 67.5, 1: 83.83, 2: 69.0, 3: 84.8}
R_M = 15.0


def kin(E):
    k = math.sqrt(2.0 * MU * E) / HBARC
    eta = Z1Z2 * ALPHA * math.sqrt(MU / (2.0 * E))
    return k, eta


def coulomb_fg(l, eta, rho):
    F = mp.coulombf(l, eta, rho)
    G = mp.coulombg(l, eta, rho)
    dF = mp.diff(lambda x: mp.coulombf(l, eta, x), rho)
    dG = mp.diff(lambda x: mp.coulombg(l, eta, x), rho)
    return float(F), float(dF), float(G), float(dG)


def V(r, l):
    e2 = ALPHA * HBARC
    vc = Z1Z2 * e2 / (2.0 * RC) * (3.0 - (r / RC) ** 2) if r < RC else Z1Z2 * e2 / r
    return vc - DEPTH[l] * math.exp(-(r / R0) ** 2) + HBARC ** 2 / (2.0 * MU) * l * (l + 1) / r ** 2


def phase(l, E):
    k2 = 2.0 * MU / HBARC ** 2
    r0 = 1e-4
    sol = solve_ivp(lambda r, y: [y[1], k2 * (V(r, l) - E) * y[0]], (r0, R_M), [r0 ** (l + 1), (l + 1) * r0 ** l],
                    method="DOP853", rtol=1e-12, atol=1e-14 * r0 ** (l + 1))
    u, du = sol.y[0, -1], sol.y[1, -1]
    k, eta = kin(E)
    F, dF, G, dG = coulomb_fg(l, eta, k * R_M)
    # u = c (cos d F + sin d G), u' = c k (cos d F' + sin d G'): tan d = (u k F' - u' F) / (u' G - u k G')
    d = math.atan2(u * k * dF - du * F, du * G - u * k * dG)
    while d > math.pi / 2:
        d -= math.pi
    while d <= -math.pi / 2:
        d += math.pi
    return d


def main():
    print("kind,l,a,E,value")
    for l in range(4):
        for a in (3.0, 5.0):
            for E in (0.5, 2.0):
                k, eta = kin(E)
                rho = k * a
                F = mp.coulombf(l, eta, rho)
                dF = mp.diff(lambda x: mp.coulombf(l, eta, x), rho)
                print(f"ellC,{l},{a},{E},{float(rho * dF / F):.15e}")
    for l in range(4):
        for E in (0.5, 1.0, 2.0, 4.0):
            print(f"delta,{l},0,{E},{phase(l, E):.15e}")


if __name__ == "__main__":
    main()
