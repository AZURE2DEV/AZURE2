#!/usr/bin/env python3
"""Reference values for physics_potential_tuning (tests/potential_tuning_test.cpp), independent of the engine's Numerov
integration, Coulomb routines, digamma and fit (docs/upgrade/potential_background.md section 13.1).

Systems (the engine's model: the potential acts from 0 to r_m = 15 fm, point-Coulomb or free waves beyond):
  he3   3He + 4He, atomic masses 3.0160293 u and 4.0026032 u, Z1 Z2 = 4, e^2 = alpha hbar c, Coulomb potential of a
        uniform sphere of radius 3.095 fm, Gaussian V_N = -V0 exp(-(r/r0)^2), r0 = 1/sqrt(0.15747) fm: the potential VD
        of Tursunov et al. (s wave 67.5 MeV, p3/2 83.589554 MeV; potential_background.md section 10).
  nalpha  n + 4He, masses 1.00866491595 u and 4.0026032 u, no charge, Gaussian V0 = 40 MeV, r0 = 2.0 fm (made up).

Rows (CSV): kind, system, l, x, value
  delta   the nuclear phase shift (rad, in (-pi/2, pi/2]) at E = x MeV: scipy solve_ivp (DOP853, rtol 1e-12) outward from
          r = 1e-4 fm with u = r^(l+1), matched at 15 fm to mpmath's coulombf/coulombg (Riccati-Bessel for no charge)
  a, r    the scattering length (fm^(2l+1)) and effective range of the effective-range function K(k^2) (Coulomb-modified
          for he3 l = 0, C0^2 k cot d + 2 k eta h(eta), h = Re psi(i eta) - ln eta by mpmath's digamma; k^(2l+1) cot d
          for nalpha), fitted by numpy with a quartic in k^2 to 40 points on E = 0.05 ... 0.4 MeV (he3) or
          k = 0.02 ... 0.2 fm^-1 (nalpha): a window and degree other than the engine's
  h       h(eta) = Re psi(i eta) - ln eta by mpmath's digamma (30 digits) at eta = x
  bound   the number of bound states (E < -0.01 MeV): the sign changes of the regular solution at E = -0.01 MeV out to
          400 fm (Sturm), with the potential depth x (MeV)
  depth   the depth (MeV) of the s wave of he3 that gives a0 = x fm with VD's two bound s states: scipy's brentq on the
          a0 above
  dw_*    a depth and Gaussian width (dw_depth, dw_width) and their s-wave a and r (dw_a, dw_r) and bound states
          (dw_bound), computed here: the engine's depth-and-range solve from the start potential is to find them back

Run: python3 make_potential_tuning_reference.py > potential_tuning.csv   (a few minutes)
"""
import math

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

mp.mp.dps = 30
HBARC = 197.3269804
UCONV = 931.49410242
ALPHA = 1.0 / 137.035999084
R_M = 15.0
SYSTEMS = {
    "he3": dict(m1=3.0160293, m2=4.0026032, z=4, rc=3.095, r0=1.0 / math.sqrt(0.15747)),
    "nalpha": dict(m1=1.00866491595, m2=4.0026032, z=0, rc=0.0, r0=2.0),
}


def mu(s):
    return s["m1"] * s["m2"] / (s["m1"] + s["m2"]) * UCONV


def V(s, V0, l, r, r0=None):
    e2 = ALPHA * HBARC
    if s["z"] == 0:
        vc = 0.0
    elif r < s["rc"]:
        vc = s["z"] * e2 / (2.0 * s["rc"]) * (3.0 - (r / s["rc"]) ** 2)
    else:
        vc = s["z"] * e2 / r
    return vc - V0 * math.exp(-(r / (r0 or s["r0"])) ** 2) + HBARC ** 2 / (2.0 * mu(s)) * l * (l + 1) / r ** 2


def solve(s, V0, l, E, rmax, dense=False, width=None):
    k2 = 2.0 * mu(s) / HBARC ** 2
    r0 = 1e-4
    return solve_ivp(lambda r, y: [y[1], k2 * (V(s, V0, l, r, width) - E) * y[0]], (r0, rmax),
                     [r0 ** (l + 1), (l + 1) * r0 ** l], method="DOP853", rtol=1e-12, atol=1e-30,
                     dense_output=dense)


def kin(s, E):
    k = math.sqrt(2.0 * mu(s) * E) / HBARC
    eta = s["z"] * ALPHA * math.sqrt(mu(s) / (2.0 * E))
    return k, eta


def phase(s, V0, l, E, width=None):
    sol = solve(s, V0, l, E, R_M, width=width)
    u, du = sol.y[0, -1], sol.y[1, -1]
    k, eta = kin(s, E)
    rho = mp.mpf(k * R_M)
    F = mp.coulombf(l, eta, rho)
    G = mp.coulombg(l, eta, rho)
    dF = mp.diff(lambda x: mp.coulombf(l, eta, x), rho)
    dG = mp.diff(lambda x: mp.coulombg(l, eta, x), rho)
    d = math.atan2(float(du * F - k * u * dF), float(k * u * dG - du * G))
    while d > math.pi / 2:
        d -= math.pi
    while d <= -math.pi / 2:
        d += math.pi
    return d


def ere(s, V0, l, width=None):
    if s["z"]:
        energies = np.linspace(0.05, 0.4, 40)
    else:
        energies = (np.linspace(0.02, 0.2, 40) * HBARC) ** 2 / (2.0 * mu(s))
    xs, ys = [], []
    for E in energies:
        k, eta = kin(s, E)
        d = phase(s, V0, l, E, width)
        if s["z"]:
            c02 = 2.0 * math.pi * eta / math.expm1(2.0 * math.pi * eta)
            h = float(mp.re(mp.digamma(1j * eta)) - mp.log(eta))
            K = c02 * k / math.tan(d) + 2.0 * k * eta * h
        else:
            K = k ** (2 * l + 1) / math.tan(d)
        xs.append(k * k)
        ys.append(K)
    c = np.polynomial.polynomial.polyfit(np.array(xs), np.array(ys), 4)
    return -1.0 / c[0], 2.0 * c[1]


def bound(s, V0, l, width=None):
    E = -0.01
    sol = solve(s, V0, l, E, 400.0, dense=True, width=width)
    r = np.linspace(1e-3, 400.0, 400001)
    u = sol.sol(r)[0]
    return int(np.sum(np.sign(u[1:]) * np.sign(u[:-1]) < 0))


def main():
    he3, nal = SYSTEMS["he3"], SYSTEMS["nalpha"]
    print("kind,system,l,x,value")
    for eta in (0.05, 0.5, 2.0, 5.0, 11.9, 12.1, 30.0, 100.0):
        print(f"h,he3,0,{eta},{float(mp.re(mp.digamma(1j * eta)) - mp.log(eta)):.17e}")
    for l, V0 in ((0, 67.5), (1, 83.589554)):
        for E in (0.05, 0.5, 2.0):
            print(f"delta,he3,{l},{E},{phase(he3, V0, l, E):.15e}")
    for l in (0, 1):
        for E in (0.05, 1.0):
            print(f"delta,nalpha,{l},{E},{phase(nal, 40.0, l, E):.15e}")
    a, r = ere(he3, 67.5, 0)
    print(f"a,he3,0,67.5,{a:.12e}")
    print(f"r,he3,0,67.5,{r:.12e}")
    for l in (0, 1):
        a, r = ere(nal, 40.0, l)
        print(f"a,nalpha,{l},40.0,{a:.12e}")
        print(f"r,nalpha,{l},40.0,{r:.12e}")
    for l, V0 in ((0, 67.5), (0, 30.0), (0, 110.0), (1, 83.589554), (2, 69.0), (3, 84.8)):
        print(f"bound,he3,{l},{V0},{bound(he3, V0, l)}")
    for V0 in (40.0, 20.0, 5.0):
        print(f"bound,nalpha,0,{V0},{bound(nal, V0, 0)}")
    target = 40.0
    lo = 67.5  # a0 falls with the depth within a branch: scan towards the target, then bracket
    step = -0.5 if ere(he3, lo, 0)[0] < target else 0.5
    hi = lo + step
    while (ere(he3, hi, 0)[0] - target) * (ere(he3, lo, 0)[0] - target) > 0:
        lo, hi = hi, hi + step
    depth = brentq(lambda v: ere(he3, v, 0)[0] - target, min(lo, hi), max(lo, hi), xtol=1e-9)
    print(f"bound,he3,0,{depth},{bound(he3, depth, 0)}")
    # depth and width for (a, r): the a and r of a chosen (V0, r0) on the branch of the start potential, which the
    # engine is to find back from its start (VD's s wave; n + 4He's 40 MeV, 2 fm)
    for name, sysname, V0, r0 in (("he3", "he3", 70.0, 2.4), ("nalpha", "nalpha", 35.0, 2.2)):
        sy = SYSTEMS[sysname]
        a, r = ere(sy, V0, 0, r0)
        print(f"dw_depth,{name},0,0,{V0}")
        print(f"dw_width,{name},0,0,{r0}")
        print(f"dw_a,{name},0,0,{a:.12e}")
        print(f"dw_r,{name},0,0,{r:.12e}")
        print(f"dw_bound,{name},0,0,{bound(sy, V0, 0, r0)}")
    print(f"depth,he3,0,{target},{depth:.12e}")


if __name__ == "__main__":
    main()
