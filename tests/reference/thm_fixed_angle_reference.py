#!/usr/bin/env python3
"""Independent reference for the fixed-angle THM observable (src/ThmAngular.cpp).

With the quantization axis along p_xA the THM entrance carries only m_l = 0.
The HOES amplitude for entrance channel spin projection nu and exit nu' is
the direct sum

    F_{nu nu'}(theta) = sum_{J l l'} sqrt(2l+1) <s nu l 0|J nu> <s' nu' l' m'|J nu>
                        x^J_{(s'l'),(sl)} Y_{l'}^{m'}(theta, 0),   m' = nu - nu',

and the observable is sum_{nu nu'} |F|^2 (per channel spins s, s'), whose
4 pi integral is sum (2J+1) |x|^2, the angle-integrated HOES cross section.
Here the sum is done literally: sympy's exact Clebsch-Gordan coefficients,
mpmath's spherical harmonics, and for the window averages mpmath.quad in
cos(theta).  The engine instead uses the Legendre (Blatt-Biedenharn Z)
form with GSL's 3j/6j symbols and Gauss-Legendre nodes; nothing is shared.

The amplitudes x come from a toy two-level R-matrix with the same assembly
as THMMatrixFunc::CalculateTHMCrossSection:

    x = sqrt(K 2 P_f) e^{i(omega_f - phi_f)} sum_{lam mu} gamma_{lam f} A_{lam mu} V_mu^{(s,l)},
    V_mu^{(s,l)} = gamma_{mu (s,l)} M_l,
    (A^-1)_{lam mu} = (E_lam - E) delta - sum_c gamma_{lam c} gamma_{mu c} (S_c - B_c + i P_c),

with made-up (but fixed) S, P, B, M_l and exit phases.  Two cases:

  "li7": 7Li + p -> a + a like: J = 2+ with two levels and entrance
         (s,l) = (1,1), (1,3), (2,1), (2,3), J = 0+ with one level (1,1);
         exit a + a, s' = 0, l' = J.  Identical exit: even L only.
  "asym": s = 1/2 entrance, s' = 1/2 exit, J = 1/2+ (two levels, l = 0,
         l' = 1), J = 1/2- (l = 1, l' = 0) and J = 3/2- (l = 1, l' = 2):
         odd L, forward-backward asymmetric.

Run:  python3 tests/reference/thm_fixed_angle_reference.py
It prints the table hard-coded in tests/reference/thm_fixed_angle_test.cpp.
"""
import cmath
import math

import mpmath as mp
from sympy import Rational, S, N
from sympy.physics.wigner import clebsch_gordan

mp.mp.dps = 30


def half(x):
    return Rational(int(round(2 * x)), 2)


def cg(j1, m1, j2, m2, j, m):
    return float(N(clebsch_gordan(half(j1), half(j2), half(j), half(m1), half(m2), half(m)), 25))


def spins(s):
    return [-s + k for k in range(int(round(2 * s)) + 1)]


def two_level(E, levels, channels, gam, S_, P_, B_):
    """A matrix of the levels; gam[lam][c]."""
    n = len(levels)
    M = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            v = (levels[a] - E) if a == b else 0
            for c in channels:
                v -= gam[a][c] * gam[b][c] * (S_[c] - B_[c] + 1j * P_[c])
            M[a, b] = v
    return M ** -1


def amplitudes(case):
    """[(J, s, l, sp, lp, x)] for the toy model at E = 1.1 MeV."""
    E = 1.1
    K = 0.37  # k_f/mu_f-like kinematic factor
    out = []
    if case == "li7":
        groups = [
            # J, levels, entrance (s,l) list, exit (sp, lp)
            (2.0, [0.8, 2.3], [(1, 1), (1, 3), (2, 1), (2, 3)], (0, 2),
             [{(1, 1): 0.61, (1, 3): -0.22, (2, 1): 0.35, (2, 3): 0.12, "f": 0.45},
              {(1, 1): -0.28, (1, 3): 0.41, (2, 1): 0.19, (2, 3): -0.33, "f": -0.52}]),
            (0.0, [1.6], [(1, 1)], (0, 0), [{(1, 1): 0.5, "f": 0.7}]),
        ]
        M = {1: 0.83, 3: -0.29}
        Sx = {1: -0.41, 3: -1.2}
        Px = {1: 0.05, 3: 0.002}
        exitP = {0: 0.9, 2: 0.55}
        exitPhase = {0: cmath.exp(0.31j), 2: cmath.exp(-0.47j)}
    else:
        groups = [
            (0.5, [0.7, 1.9], [(0.5, 0)], (0.5, 1),
             [{(0.5, 0): 0.52, "f": 0.3}, {(0.5, 0): -0.44, "f": 0.21}]),
            (0.5, [1.4], [(0.5, 1)], (0.5, 0), [{(0.5, 1): 0.47, "f": -0.36}]),
            (1.5, [1.25], [(0.5, 1)], (0.5, 2), [{(0.5, 1): 0.33, "f": 0.4}]),
        ]
        M = {0: 1.07, 1: 0.64}
        Sx = {0: -0.2, 1: -0.6}
        Px = {0: 0.3, 1: 0.08}
        exitP = {0: 1.1, 1: 0.8, 2: 0.35}
        exitPhase = {0: cmath.exp(0.2j), 1: cmath.exp(1.1j), 2: cmath.exp(-0.9j)}
    for J, lev, ent, (sp, lp), gam in groups:
        chans = list(ent) + ["f"]
        S_ = {c: (Sx[c[1]] if c != "f" else -0.3) for c in chans}
        P_ = {c: (Px[c[1]] if c != "f" else exitP[lp]) for c in chans}
        B_ = {c: (Sx[c[1]] - 0.05 if c != "f" else -0.3) for c in chans}
        A = two_level(E, lev, chans, gam, S_, P_, B_)
        n = len(lev)
        for (s, l) in ent:
            amp = 0
            for a in range(n):
                for b in range(n):
                    amp += gam[a]["f"] * A[a, b] * gam[b][(s, l)] * M[l]
            x = complex(mp.sqrt(K * 2 * exitP[lp])) * exitPhase[lp] * complex(amp)
            out.append((J, s, l, sp, lp, x))
    return out


def direct(amps, theta):
    """sum_{nu nu'} |F_{nu nu'}(theta)|^2 by the literal M-sum."""
    tot = mp.mpf(0)
    for s in sorted({a[1] for a in amps}):
        for sp in sorted({a[3] for a in amps}):
            for nu in spins(s):
                for nup in spins(sp):
                    F = mp.mpc(0)
                    for (J, s1, l, sp1, lp, x) in amps:
                        if s1 != s or sp1 != sp or abs(nu) > J:
                            continue
                        m = nu - nup
                        if abs(m) > lp:
                            continue
                        c = math.sqrt(2 * l + 1) * cg(s, nu, l, 0, J, nu) * cg(sp, nup, lp, m, J, nu)
                        if c == 0:
                            continue
                        F += c * x * mp.spherharm(lp, int(round(m)), theta, 0)
                    tot += abs(F) ** 2
    return tot


def axis(amps):
    """theta = 0: Y_l'^m'(0) = delta_m'0 sqrt((2l'+1)/4 pi), so nu' = nu --
    the m_l = 0 amplitude along the axis, summed over nu."""
    tot = mp.mpf(0)
    for s in sorted({a[1] for a in amps}):
        for sp in sorted({a[3] for a in amps}):
            for nu in spins(s):
                if abs(nu) > sp:
                    continue
                F = mp.mpc(0)
                for (J, s1, l, sp1, lp, x) in amps:
                    if s1 != s or sp1 != sp or abs(nu) > J:
                        continue
                    F += mp.sqrt((2 * l + 1) * (2 * lp + 1) / (4 * mp.pi)) * cg(s, nu, l, 0, J, nu) * \
                        cg(sp, nu, lp, 0, J, nu) * x
                tot += abs(F) ** 2
    return tot


def window(amps, t1, t2):
    x1, x2 = mp.cos(mp.radians(t2)), mp.cos(mp.radians(t1))
    f = lambda c: direct(amps, mp.acos(c))
    return mp.quad(f, [x1, (x1 + x2) / 2, x2]) / (x2 - x1)


def legendre(amps, L):
    """b_L with sum|F|^2 = (1/pi) sum_L b_L P_L: pi (2L+1)/2 Int sum|F|^2 P_L dx."""
    f = lambda c: direct(amps, mp.acos(c)) * mp.legendre(L, c)
    return mp.pi * (2 * L + 1) / 2 * mp.quad(f, [-1, 0, 1])


def main():
    for case in ("li7", "asym"):
        amps = amplitudes(case)
        print(f"// case {case}: waves {{J, s, l, s', l'}}, x")
        for (J, s, l, sp, lp, x) in amps:
            print(f"  {{{{{J}, {s}, {l}, {sp}, {lp}}}, {{{x.real:.17e}, {x.imag:.17e}}}}},")
        sigma = sum((2 * a[0] + 1) * abs(a[5]) ** 2 for a in amps)
        print(f"  sigma = sum (2J+1)|x|^2 = {sigma:.17e}")
        maxL = 2 * max(a[4] for a in amps)
        for L in range(maxL + 1):
            print(f"  b_{L} = {mp.nstr(legendre(amps, L), 17)}")
        for t in (0, 37, 90, 143, 180):
            print(f"  dsdO({t}) = {mp.nstr(direct(amps, mp.radians(t)), 17)}")
        print(f"  axis(0) = {mp.nstr(axis(amps), 17)}")
        for (t1, t2) in ((50, 70), (0, 180), (0, 0.5), (20, 160), (110, 130)):
            print(f"  window({t1}-{t2}) = {mp.nstr(window(amps, t1, t2), 17)}")


if __name__ == "__main__":
    main()
