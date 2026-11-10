#!/usr/bin/env python3
"""Reference Rutherford and Mott cross sections from mpmath (30 digits), for tests/coulomb_scattering_test.cpp.

Independent of physics/: eta and k from the definitions (non-relativistic k = sqrt(2 mu E)/hbar c,
eta = Z1 Z2 alpha sqrt(mu / 2E); relativistic: s = (m + T)^2, k = p/hbar c with p^2 = (s - m^2)(s - (m1 - m2)^2)/(4s),
eta = Z1 Z2 alpha / beta_rel, beta_rel = p sqrt(s) / (s/2 - (m1^2 + m2^2)/2)), the cross sections from the Coulomb
amplitude f_C(theta) = -eta / (2k sin^2(theta/2)) exp(-i eta ln sin^2(theta/2)):

  Rutherford  |f_C(theta)|^2
  Mott        sum_{S=0}^{2s} (2S + 1) |f_C(theta) + (-1)^S f_C(pi - theta)|^2 / (2s + 1)^2   (two identical spin-s particles)

and d/dE by mpmath's numerical derivative of the same expressions.  Constants: include/Constants.h (CODATA 2018).
Rows: z1 z2 m1(u) m2(u) two_s(-1 = Rutherford) relativistic E_cm(MeV) theta_cm(deg) value(b/sr) d_dE(b/sr/MeV).
Usage: make_coulomb_scattering.py OUTFILE
"""
import sys

import mpmath as mp

mp.mp.dps = 30
HBARC = mp.mpf("197.3269804")
UCONV = mp.mpf("931.49410242")
ALPHA = 1 / mp.mpf("137.035999084")


def eta_k(z1z2, m1u, m2u, rel, e):
    m1, m2 = m1u * UCONV, m2u * UCONV
    if not rel:
        mu = m1 * m2 / (m1 + m2)
        return z1z2 * ALPHA * mp.sqrt(mu / (2 * e)), mp.sqrt(2 * mu * e) / HBARC
    m = m1 + m2
    s = (m + e) ** 2
    p = mp.sqrt((s - m ** 2) * (s - (m1 - m2) ** 2) / (4 * s))
    beta = p * mp.sqrt(s) / ((s - m1 ** 2 - m2 ** 2) / 2)  # relative velocity of the pair (c = 1)
    return z1z2 * ALPHA / beta, p / HBARC


def fc(eta, k, th):
    s2 = mp.sin(th / 2) ** 2
    return -eta / (2 * k * s2) * mp.expj(-eta * mp.log(s2))


def cross_section(z1z2, m1u, m2u, two_s, rel, e, theta_deg):
    eta, k = eta_k(z1z2, m1u, m2u, rel, e)
    th = mp.radians(theta_deg)
    if two_s < 0:
        return abs(fc(eta, k, th)) ** 2 / 100
    direct, exchange = fc(eta, k, th), fc(eta, k, mp.pi - th)
    total = mp.mpf(0)
    for S in range(two_s + 1):
        total += (2 * S + 1) * abs(direct + (-1) ** S * exchange) ** 2
    return total / (two_s + 1) ** 2 / 100


CASES = [
    # z1 z2  m1          m2           two_s rel
    (1, 6, "1.00782503", "12.0", -1, 0),        # p + 12C, Rutherford
    (2, 6, "4.00260325", "12.0", -1, 0),        # alpha + 12C
    (2, 6, "4.00260325", "12.0", -1, 1),        # alpha + 12C, relativistic kinematics
    (2, 2, "4.00260325", "4.00260325", 0, 0),   # alpha + alpha, Mott (bosons, spin 0)
    (2, 2, "4.00260325", "4.00260325", 0, 1),
    (1, 1, "1.00782503", "1.00782503", 1, 0),   # p + p, Mott (fermions, spin 1/2)
    (1, 1, "2.01410178", "2.01410178", 2, 0),   # d + d, Mott (bosons, spin 1)
    (6, 6, "12.0", "12.0", 0, 0),               # 12C + 12C
]
ENERGIES = ["0.3", "2.5", "40"]
ANGLES = ["10", "37.5", "90", "143", "170"]

with open(sys.argv[1], "w") as out:
    out.write("# z1 z2 m1 m2 two_s(-1 Rutherford) relativistic E_cm theta_cm value(b/sr) d_dE(b/sr/MeV)  "
              "make_coulomb_scattering.py, mpmath %s\n" % mp.__version__)
    for z1, z2, m1, m2, two_s, rel in CASES:
        for es in ENERGIES:
            for ths in ANGLES:
                f = lambda e: cross_section(z1 * z2, mp.mpf(m1), mp.mpf(m2), two_s, rel, e, mp.mpf(ths))
                v = f(mp.mpf(es))
                d = mp.diff(f, mp.mpf(es))
                out.write(f"{z1} {z2} {m1} {m2} {two_s} {rel} {es} {ths} {mp.nstr(v, 20)} {mp.nstr(d, 20)}\n")
