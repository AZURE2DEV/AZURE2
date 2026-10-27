#!/usr/bin/env python3
"""Capture analyzing power of models/twolevel.azr against a magnetic-substate calculation.

Truth.  The two 1+ levels share no particle channel, so (Brune parametrization) the capture
amplitudes are  T_s(E) = gamma_s gamma_gamma / D_lambda(E),
D = E_lambda - E - gamma_s^2 [S(E) - S(E_lambda)] - i gamma_s^2 P(E),  s = 0 from level A and s = 1
from level B, with the engine's own reduced-width amplitudes (parameters.out) and the closed-form
neutral l = 1 functions S = -1/(1 + rho^2), P = rho^3/(1 + rho^2), rho = k a.  Every other factor
(i^l, the hard-sphere phase, the photon penetrability, the overall normalisation) is common to
the two amplitudes and drops out of A_y.

A_y is then formed in magnetic substates: the channel-spin amplitude
  A(s nu, m_f, lambda; theta) = sum over pathways with that s of
     i sqrt(pi) sqrt(2l+1) <s nu l 0|J nu> <I_f m_f L mu|J nu> T sqrt((2L+1)/4pi) d^L_{mu lambda}(theta) c_p/sqrt2
(the convention under which the unpolarized Seyler-Weller coefficients reproduce AZURE2's own
Blatt-Biedenharn cross section, 02_angular_distributions.md of the 2026-10-25 audit), the projectile
spin taken out of the channel spin LIGHT PARTICLE FIRST, <j1 m1 j2 m2|s nu> (Lane and Thomas, as
PolarizationFunc), and  A_y = 2 Im sum (A_up A_down^*) / sum (|A_up|^2 + |A_down|^2)
summed over m_2, m_f and the photon helicity.  With the channel spin coupled target first the sign of
every s = 0 / s = 1 cross term flips, which here is A_y -> -A_y: that was the engine before 2026-10-26.

  python3 verify.py WORKDIR   (WORKDIR/output/AZUREOut_aa=1_R=2.extrap and WORKDIR/output/parameters.out)
"""
import os
import re
import sys
from math import cos, factorial, pi, sin, sqrt

import numpy as np
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_models as M  # noqa: E402

HBARC, UCONV = 197.32696310, 931.4940880   # include/Constants.h


def R(x):
    return Rational(int(round(2 * x)), 2)


def CG(j1, m1, j2, m2, J, Mm):
    try:
        return float(clebsch_gordan(R(j1), R(j2), R(J), R(m1), R(m2), R(Mm)))
    except ValueError:
        return 0.0


def dsmall(j, mp, m, beta):
    """Wigner d^j_{mp m}(beta)."""
    j2 = int(round(2 * j))
    jpm, jmm, jpmp, jmmp = (int(round(j + m)), int(round(j - m)), int(round(j + mp)), int(round(j - mp)))
    pref = sqrt(factorial(jpmp) * factorial(jmmp) * factorial(jpm) * factorial(jmm))
    mpm = int(round(mp - m))
    res = 0.0
    for k in range(0, j2 + 1):
        a, b, c, d = jpm - k, k, jmmp - k, k + mpm
        if min(a, b, c, d) < 0:
            continue
        res += ((-1) ** (k + mpm) * cos(beta / 2) ** (2 * j - 2 * k - mpm) * sin(beta / 2) ** (2 * k + mpm)
                / (factorial(a) * factorial(b) * factorial(c) * factorial(d)))
    return pref * res


def mrange(j):
    return [-j + i for i in range(int(round(2 * j)) + 1)]


def substate_ay(paths, i1, i2, If, theta, cM=-1.0):
    """A_y from channel-spin pathways {J, s, l, L, p, T}; projectile i1 first in <i1 m1 i2 m2|s nu>."""
    sins = [abs(i1 - i2) + k for k in range(int(round(i1 + i2 - abs(i1 - i2))) + 1)]

    def A(s, nu, mf, lam):
        out = 0j
        for p in paths:
            if p["s"] != s:
                continue
            mu = nu - mf
            J, l, L = p["J"], p["l"], p["L"]
            if abs(mu) > L or abs(nu) > J:
                continue
            cp = 1.0 if p["p"] == 1 else cM * lam
            out += (1j * sqrt(pi) * sqrt(2 * l + 1) * CG(s, nu, l, 0, J, nu) * CG(If, mf, L, mu, J, nu) * p["T"]
                    * sqrt((2 * L + 1) / (4 * pi)) * dsmall(L, mu, lam, theta) * cp / sqrt(2))
        return out

    inter, den = 0j, 0.0
    for mf in mrange(If):
        for lam in (1, -1):
            for m2 in mrange(i2):
                up = sum(CG(i1, 0.5, i2, m2, s, 0.5 + m2) * A(s, 0.5 + m2, mf, lam) for s in sins if abs(0.5 + m2) <= s)
                dn = sum(CG(i1, -0.5, i2, m2, s, -0.5 + m2) * A(s, -0.5 + m2, mf, lam) for s in sins if abs(-0.5 + m2) <= s)
                inter += up * np.conj(dn)
                den += abs(up) ** 2 + abs(dn) ** 2
    return 2 * inter.imag / den


def read_gammas(path):
    """{(E_level, R, s): g_int} of the 1+ levels from parameters.out."""
    out, level = {}, None
    for line in open(path):
        m = re.match(r"\s*J =\s*([0-9.]+)([+-])\s+E_level =\s*([0-9.]+)", line)
        if m:
            level = (float(m.group(1)), m.group(2), float(m.group(3)))
            continue
        m = re.match(r"\s*R =\s*(\d+)\s+[lL] =\s*(\S+)\s+s =\s*([0-9.]+).*?g_int =\s*(-?[0-9.]+)", line)
        if m and level and level[0] == 1.0 and level[1] == "+":
            out[(level[2], int(m.group(1)), float(m.group(3)))] = float(m.group(4))
    return out


def neutral_l1(e_cm):
    mu = M.MN * M.MT / (M.MN + M.MT)
    k = sqrt(2.0 * mu * UCONV * e_cm) / HBARC
    rho = k * M.RADIUS
    return -1.0 / (1.0 + rho * rho), rho ** 3 / (1.0 + rho * rho)


def amplitude(e_cm, e_level, g_n, g_gam):
    s_e, p_e = neutral_l1(e_cm)
    s_l, _ = neutral_l1(e_level)
    d = (e_level - e_cm) - g_n * g_n * (s_e - s_l) - 1j * g_n * g_n * p_e
    return g_n * g_gam / d


def main(work):
    failures = 0
    g = read_gammas(os.path.join(work, "output", "parameters.out"))
    ea, eb = round(M.EA, 4), round(M.EB, 4)
    try:
        gA0, gAg, gB1, gBg = g[(ea, 1, 0.0)], g[(ea, 2, 0.0)], g[(eb, 1, 1.0)], g[(eb, 2, 0.0)]
    except KeyError as e:
        print(f"  FAIL  parameters.out: reduced-width amplitude {e} not found; have {sorted(g)}")
        return 1
    print(f"  reduced-width amplitudes (MeV^1/2): A s=0 {gA0:.6f}, A gamma {gAg:.3e}; B s=1 {gB1:.6f}, B gamma {gBg:.3e}")
    rows = np.loadtxt(os.path.join(work, "output", "AZUREOut_aa=1_R=2.extrap"))
    # block 1 is the A_y segment; a blank line separates it from the cross-section segment (loadtxt skips blanks),
    # so take the first half
    n = rows.shape[0] // 2
    ay_rows = rows[:n]
    worst, biggest = 0.0, 0.0
    for e_cm, theta_deg, ay in ay_rows[:, [0, 2, 3]]:
        paths = [dict(J=1, s=0.0, l=1, L=1, p=0, T=amplitude(e_cm, M.EA - M.S_N, gA0, gAg)),
                 dict(J=1, s=1.0, l=1, L=1, p=0, T=amplitude(e_cm, M.EB - M.S_N, gB1, gBg))]
        truth = substate_ay(paths, 0.5, 0.5, 0.0, np.radians(theta_deg))
        worst = max(worst, abs(ay - truth))
        biggest = max(biggest, abs(truth))
    print(f"  {n} points, |A_y| up to {biggest:.4f}, worst |engine - truth| {worst:.2e}")
    if biggest < 0.05:
        print("  FAIL  the test is vacuous: |A_y| < 0.05 everywhere"); failures += 1
    if worst <= 1e-5:
        print(f"  ok    capture A_y = light-particle-first substate truth (worst {worst:.2e})")
    else:
        print(f"  FAIL  capture A_y vs substate truth: worst {worst:.2e} (> 1e-5)"); failures += 1
        # diagnose the coupling order
        flipped = max(abs(ay + substate_ay([dict(J=1, s=0.0, l=1, L=1, p=0, T=amplitude(e, M.EA - M.S_N, gA0, gAg)),
                                            dict(J=1, s=1.0, l=1, L=1, p=0, T=amplitude(e, M.EB - M.S_N, gB1, gBg))],
                                           0.5, 0.5, 0.0, np.radians(t))) for e, t, ay in ay_rows[:, [0, 2, 3]])
        if flipped <= 1e-5:
            print(f"        the engine's A_y is the NEGATIVE of the truth (target-first coupling, worst {flipped:.2e})")
    ends = ay_rows[(ay_rows[:, 2] == 0.0) | (ay_rows[:, 2] == 180.0), 3]
    if np.all(np.abs(ends) < 1e-12):
        print("  ok    A_y = 0 at 0 and 180 degrees")
    else:
        print(f"  FAIL  A_y at 0/180 degrees: {ends}"); failures += 1
    return failures


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
