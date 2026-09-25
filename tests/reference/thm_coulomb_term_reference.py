#!/usr/bin/env python3
"""Independent mpmath reference for ThmCoulombTerm (src/ThmFunc.cpp).

The THM entrance vertex carries, besides the surface part
(B - 1) j_l(rho) - rho j_l'(rho), the external Coulomb term

    C_l(E) = 2 eta k  Int_a^inf dr  O_l(k r) / O_l(k a) * j_l(p r)

(Tribble et al., Rep. Prog. Phys. 77 (2014) 106901, eq. 2.79;
Mukhamedzhanov, Shubhchintak, Bertulani, PRC 96 (2017) 024623, eq. 27;
Typel & Baur, Ann. Phys. 305 (2003) 228, eq. A.4), with 2 eta k =
2 Z1 Z2 alpha mu c^2 / (hbar c) and O_l = G_l + i F_l = H+_l.  For E < 0,
O_l -> W_{-eta_b, l+1/2}(2 kappa r); at E = 0 the common limit
sqrt(r) K_{2l+1}(sqrt(8 c r)), c = eta k.

The radial integral converges only conditionally for E > 0 (integrand
~ 1/r times oscillation).  Here it is evaluated WITHOUT any asymptotic
expansion: j_l = (h1_l + h2_l)/2 and the two pieces are integrated along the
rays r = a + i t (h1 part, decays like exp(-(p+k) t)) and r = a - i t (h2 part,
decays like exp(-(p-k) t)), which is exact by Cauchy/Jordan since the
integrand is analytic for Re r > 0.  O_l at complex r comes from
mpmath.whitw via DLMF 33.2.7, H+ ~ W_{-i eta, l+1/2}(-2 i rho) (the constant
prefactor cancels in the ratio).  As a consistency check a second evaluation
splits at a real R = a + 15 fm: real-axis quadrature a..R plus the two rays
from R; the two must agree.

The C++ implementation uses a different method (real-axis Gauss-Legendre
panels with O_l integrated inwards from its asymptotic expansion, rays only
beyond the asymptotic radius; WhitFunc for E < 0), so agreement is a genuine
cross-check.

Run:  python3 tests/reference/thm_coulomb_term_reference.py [--check] [pair ...]
It prints the table hard-coded in thm_coulomb_term_test.cpp (about half an
hour for all pairs, most of it 12C+12C; naming pairs restricts the run).
"""
import sys
import numpy
import mpmath as mp

mp.mp.dps = 25

# include/Constants.h
HBARC = mp.mpf('197.32696310')
UCONV = mp.mpf('931.4940880')
FSTRUC = 1 / mp.mpf('137.0359996790')

PAIRS = {
    # name: (Z1, Z2, m1 [u], m2 [u], a [fm], B [MeV])
    '7Li+p':   (3, 1, '7.016003', '1.007276', '4.0', '2.2246'),
    '6Li+d':   (3, 1, '6.015123', '2.014102', '4.5', '1.4735'),
    '12C+12C': (6, 6, '12.0', '12.0', '6.0', '10.27'),
    '17O+n':   (8, 0, '16.999131', '1.008665', '5.0', '4.143'),
}
ENERGIES = ['-1.0', '-0.1', '0.0', '0.05', '0.3', '1.0', '3.0', '6.0']
LVALUES = [0, 1, 2, 3]


def sph_h(l, x, sign):
    """Spherical Hankel h^(1) (sign=+1) or h^(2) (sign=-1), complex x."""
    s = mp.mpc(0)
    for m in range(l + 1):
        s += (sign * 1j) ** m * mp.factorial(l + m) / (mp.factorial(m) * mp.factorial(l - m) * (2 * x) ** m)
    return (-sign * 1j) ** (l + 1) * mp.exp(sign * 1j * x) / x * s


def make_O(l, E, mu, c):
    """Return O(r) up to an r-independent factor (analytic in Re r > 0)."""
    if E > 0:
        k = mp.sqrt(2 * mu * E) / HBARC
        eta = c / k
        return lambda r: mp.whitw(-1j * eta, l + mp.mpf(1) / 2, -2j * k * r)
    if E < 0:
        kap = mp.sqrt(2 * mu * (-E)) / HBARC
        etab = c / kap
        return lambda r: mp.whitw(-etab, l + mp.mpf(1) / 2, 2 * kap * r)
    return lambda r: mp.sqrt(r) * mp.besselk(2 * l + 1, mp.sqrt(8 * c * r))


_GL_X, _GL_W = numpy.polynomial.legendre.leggauss(24)


def gl_panels(f, pts):
    """Composite 24-point Gauss-Legendre over consecutive breakpoints."""
    total = mp.mpc(0)
    for lo, hi in zip(pts[:-1], pts[1:]):
        half, mid = (hi - lo) / 2, (hi + lo) / 2
        total += half * mp.fsum(w * f(mid + half * x) for x, w in zip(_GL_X, _GL_W))
    return total


def ray_breaks(rate, tend):
    """Breakpoints 0 < t1 < ... >= tend, panels growing geometrically."""
    pts, h = [mp.mpf(0)], min(mp.mpf('0.5'), mp.mpf('0.5') / rate)
    while pts[-1] < tend:
        pts.append(pts[-1] + h)
        h = min(h * mp.mpf('1.6'), 5 / rate)
    return pts


def coulomb_term(name, E, l, split=0):
    Z1, Z2, m1, m2, a, B = PAIRS[name]
    if Z1 * Z2 == 0:
        return mp.mpc(0)
    m1, m2, a, B, E = mp.mpf(m1), mp.mpf(m2), mp.mpf(a), mp.mpf(B), mp.mpf(E)
    mu = m1 * m2 / (m1 + m2) * UCONV
    c = Z1 * Z2 * FSTRUC * mu / HBARC            # eta k, fm^-1
    p = mp.sqrt(2 * mu * (E + B)) / HBARC
    O = make_O(l, E, mu, c)
    Oa = O(a)
    R = a + split
    total = mp.mpc(0)
    if split > 0:
        jl = lambda x: mp.besselj(l + mp.mpf(1) / 2, x) * mp.sqrt(mp.pi / (2 * x))
        total += gl_panels(lambda r: O(r) * jl(p * r), mp.linspace(a, R, 31))
    # j_l = (h1 + h2)/2; h1 part on r = R + i t, h2 part on r = R - i t.
    if E > 0:
        k = mp.sqrt(2 * mu * E) / HBARC
        eta = c / k
        rate_up, rate_dn = p + k, p - k
        tend_up = (60 + eta * mp.pi / 2) / rate_up
        tend_dn = 60 / rate_dn
    else:
        rate_up = rate_dn = p
        tend_up = tend_dn = 60 / p
    up = lambda t: O(R + 1j * t) * sph_h(l, p * (R + 1j * t), +1) * 1j
    dn = lambda t: O(R - 1j * t) * sph_h(l, p * (R - 1j * t), -1) * (-1j)
    total += (gl_panels(up, ray_breaks(rate_up, tend_up)) +
              gl_panels(dn, ray_breaks(rate_dn, tend_dn))) / 2
    return 2 * c * total / Oa


def surface(name, E, l):
    """(B-1) j_l(rho) - rho j_l'(rho) with B = 0, for the magnitude table."""
    Z1, Z2, m1, m2, a, B = PAIRS[name]
    m1, m2, a, B, E = mp.mpf(m1), mp.mpf(m2), mp.mpf(a), mp.mpf(B), mp.mpf(E)
    mu = m1 * m2 / (m1 + m2) * UCONV
    rho = mp.sqrt(2 * mu * (E + B)) / HBARC * a
    j = lambda x: mp.besselj(l + mp.mpf(1) / 2, x) * mp.sqrt(mp.pi / (2 * x))
    return -j(rho) - rho * mp.diff(j, rho)


if __name__ == '__main__':
    check = '--check' in sys.argv
    only = [a for a in sys.argv[1:] if a in PAIRS]  # optional: restrict to named pairs
    print('// name, E (MeV), l, Re C, Im C   (generated by thm_coulomb_term_reference.py)')
    worst = 0
    for name in (only or PAIRS):
        for E in ENERGIES:
            for l in LVALUES:
                v = coulomb_term(name, E, l)
                if check and PAIRS[name][0] * PAIRS[name][1] != 0:
                    w = coulomb_term(name, E, l, split=15)
                    d = abs(v - w) / max(abs(v), mp.mpf('1e-30'))
                    worst = max(worst, d)
                    print('//   split check rel diff %.2e' % float(d))
                print('    {"%s", %s, %d, %s, %s},' % (name, E, l,
                      mp.nstr(v.real, 17, min_fixed=-3, max_fixed=3),
                      mp.nstr(v.imag, 17, min_fixed=-3, max_fixed=3)))
                sys.stdout.flush()
    if check:
        print('// worst split-check relative difference: %.2e' % float(worst))
