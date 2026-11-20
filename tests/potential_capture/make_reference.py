#!/usr/bin/env python3
"""Direct potential-model capture 3He(alpha,gamma)7Be, computed without the engine (docs/upgrade/potential_background.md
section 5.6): the reference of tests/potential_capture.

The potential: V(r) = -V0 exp(-alpha r^2) + V_C(r), alpha = 0.15747 fm^-2 (Gaussian width r0 = 1/sqrt(alpha)), V_C of
a uniformly charged sphere of R_C = 3.095 fm, the depths VD of Tursunov, Turakulov and Kadyrov, Phys. Rev. C 97 (2018)
035802 (arXiv:1707.05948v3, Tables I and II): s1/2 67.5, p1/2 81.815179, p3/2 83.589554, d3/2 66.0, d5/2 69.0 MeV.
Masses 3.0160293 u (3He) and 4.0026032 u (4He), 931.49410242 MeV/u, hbar c = 197.3269804 MeV fm, alpha =
1/137.035999084 (the engine's constants and the projects' masses; not the paper's hbar^2/2u).

Bound states: p3/2 (7Be 3/2-) and p1/2 (1/2-), one node, by shooting (scipy DOP853 outward from the origin and inward
from 30 fm from mpmath's Whittaker function, matched at 3 fm, brentq on the energy); normalised on [0, inf) (mpmath
beyond 30 fm), so the ANC C is that of the normalised state, u_f = C W_{-eta_b, l+1/2}(2 kappa r) outside the
potential.

Scattering waves: DOP853 outward from the origin to 15 fm, matched there to cos(d) F_l + sin(d) G_l (mpmath coulombf,
coulombg), and normalised to that form (unit amplitude): u_E = cos d F + sin d G beyond 15 fm.

E1, long-wavelength operator e Z_eff r Y_1, Z_eff = Z1 m2/M - Z2 m1/M (m1 = 3He):

    sigma(E) = (16 pi / 9) k_g^3 / (hbar v) * 1/(2 s + 1) * (4 pi / k^2)
               * sum_{l_i J_i} e^2 Z_eff^2 |<(l_f s) J_f || Y_1 || (l_i s) J_i>|^2 |I_{l_i J_i}|^2,
    I = int_0^inf u_f(r) r u_E(r) dr,
    |<(l_f s)J_f||Y_L||(l_i s)J_i>|^2 = (2J_i+1)(2J_f+1) {l_f J_f s; J_i l_i L}^2 (2l_f+1)(2L+1)(2l_i+1)/(4 pi) (l_f L l_i; 0 0 0)^2,

(Fermi's golden rule, B(EL) rate 8 pi (L+1) k_g^(2L+1) / (hbar L ((2L+1)!!)^2), the incident wave normalised to unit
flux, averaged over the 2 s + 1 = 2 spin states of 3He + 4He), k_g = E_g / hbar c, E_g = E + B_f (no recoil, as the
engine), hbar v = (hbar c)^2 k / mu c^2.  S(E) = E sigma exp(2 pi eta).

Output (stdout, CSV): the bound states (B, ANC) as comment lines, then per energy sigma (b) and S (MeV b) for each final
state and their sum.  Energies are kept above 10 keV (mpmath's Coulomb functions slow down without bound as eta grows).

With --s-depth=V, the s1/2 depth V (MeV) in place of 67.5 (check.sh section 8: 100 MeV, a narrow s-wave resonance at
E_cm = 0.30 MeV; data/reference_narrow.csv).

With --paper, the paper's hbar^2/2u = 20.7343 MeV fm^2 and nuclear masses (4.001506179127 u, 3.0149322473 u, its
section III.A), to compare with its numbers (B = 1.5866 MeV, C(3/2-) = 4.34 fm^-1/2, Fig. 6).

Run: python3 -I make_reference.py > reference.csv   (a few minutes)
"""
import math
import sys
from fractions import Fraction

import mpmath as mp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

mp.mp.dps = 20
HBARC, UCONV, ALPHA = 197.3269804, 931.49410242, 1.0 / 137.035999084
M1, M2 = 3.0160293, 4.0026032
PAPER = "--paper" in sys.argv  # the paper's constants: hbar^2/2u = 20.7343 MeV fm^2, nuclear masses (its section III.A)
if PAPER:
    M1, M2 = 3.0149322473, 4.001506179127
MU = M1 * M2 / (M1 + M2) * UCONV
Z1, Z2 = 2, 2
E2 = ALPHA * HBARC
RC, GAUSS = 3.095, 0.15747
DEPTH = {(0, 1): 67.5, (1, 1): 81.815179, (1, 3): 83.589554, (2, 3): 66.0, (2, 5): 69.0}
# --s-depth=V: the s1/2 depth V (MeV) in place of 67.5 (the narrow-resonance case of check.sh section 8: V = 100 MeV
# puts an s-wave resonance at E_cm = 0.30 MeV, about 20 keV wide)
for _a in sys.argv[1:]:
    if _a.startswith("--s-depth="):
        DEPTH[(0, 1)] = float(_a.split("=", 1)[1])
R_MATCH, R_OUT = 15.0, 30.0
K2 = (M1 * M2 / (M1 + M2) / 20.7343) if PAPER else 2.0 * MU / HBARC ** 2  # 2 mu / hbar^2 (fm^-2 MeV^-1)


def V(r, l, two_j):
    vc = Z1 * Z2 * E2 / (2.0 * RC) * (3.0 - (r / RC) ** 2) if r < RC else Z1 * Z2 * E2 / r
    return -DEPTH[(l, two_j)] * math.exp(-GAUSS * r * r) + vc + l * (l + 1) / (K2 * r * r)


def rhs(l, two_j, E):
    return lambda r, y: [y[1], K2 * (V(r, l, two_j) - E) * y[0]]


def outward(l, two_j, E, r_end):
    r0 = 1e-4
    return solve_ivp(rhs(l, two_j, E), (r0, r_end), [r0 ** (l + 1), (l + 1) * r0 ** l], method="DOP853", rtol=1e-12,
                     atol=1e-30, dense_output=True)


def whittaker(l, B):
    kappa = math.sqrt(K2 * B)
    eta_b = Z1 * Z2 * E2 * K2 / (2.0 * kappa)
    return lambda r: mp.whitw(-eta_b, l + 0.5, 2.0 * kappa * r)


def bound_state(l, two_j, B_lo, B_hi):
    r_m = 3.0

    def mismatch(B):
        o = outward(l, two_j, -B, r_m)
        W = whittaker(l, B)
        i = solve_ivp(rhs(l, two_j, -B), (R_OUT, r_m), [float(W(R_OUT)), float(mp.diff(W, R_OUT))], method="DOP853",
                      rtol=1e-12, atol=1e-40)
        return o.y[1, -1] / o.y[0, -1] - i.y[1, -1] / i.y[0, -1]

    B = brentq(mismatch, B_lo, B_hi, xtol=1e-13)
    W = whittaker(l, B)
    o = outward(l, two_j, -B, r_m)
    i = solve_ivp(rhs(l, two_j, -B), (R_OUT, r_m), [float(W(R_OUT)), float(mp.diff(W, R_OUT))], method="DOP853",
                  rtol=1e-12, atol=1e-40, dense_output=True)
    scale = i.y[0, -1] / o.y[0, -1]  # the outward solution in units of the unit-ANC exterior

    def u(r):  # unit ANC
        if r <= r_m:
            return scale * o.sol(r)[0]
        if r <= R_OUT:
            return i.sol(r)[0]
        return float(W(r))

    norm = (quad(lambda r: u(r) ** 2, 0.0, r_m, epsrel=1e-12, limit=200)[0] +
            quad(lambda r: u(r) ** 2, r_m, R_OUT, epsrel=1e-12, limit=200)[0] +
            float(mp.quad(lambda r: W(r) ** 2, [R_OUT, 60, mp.inf])))
    C = 1.0 / math.sqrt(norm)
    return B, C, (lambda r: C * u(r)), (lambda r: C * W(r))


def scattering(l, two_j, E):
    k = math.sqrt(K2 * E)
    eta = Z1 * Z2 * E2 * K2 / (2.0 * k)
    o = outward(l, two_j, E, R_MATCH)
    rho = k * R_MATCH
    F, G = float(mp.coulombf(l, eta, rho)), float(mp.coulombg(l, eta, rho))
    dF = float(mp.diff(lambda x: mp.coulombf(l, eta, x), rho))
    dG = float(mp.diff(lambda x: mp.coulombg(l, eta, x), rho))
    u, du = o.y[0, -1], o.y[1, -1] / k  # du/drho
    # u = A (c F + s G), du = A (c F' + s G'), with F'G - FG' = 1: A c = du G - u G', A s = u F' - du F
    ac, as_ = du * G - u * dG, u * dF - du * F
    A = math.hypot(ac, as_)
    if ac < 0:
        A = -A
    c, s = ac / A, as_ / A
    interior = lambda r: o.sol(r)[0] / A
    exterior = lambda r: c * mp.coulombf(l, eta, k * r) + s * mp.coulombg(l, eta, k * r)
    return k, eta, interior, exterior


def three_j_000(a, b, c):
    J = a + b + c
    if J % 2 or c > a + b or c < abs(a - b):
        return 0.0
    g = J // 2
    f = math.factorial
    val = (-1) ** g * math.sqrt(f(J - 2 * a) * f(J - 2 * b) * f(J - 2 * c) / f(J + 1)) * f(g) / (f(g - a) * f(g - b) * f(g - c))
    return val


def six_j(j1, j2, j3, j4, j5, j6):
    """{j1 j2 j3; j4 j5 j6} by Racah's formula; arguments as Fractions (half-integers allowed)."""
    def tri(a, b, c):
        return (a + b - c >= 0 and a - b + c >= 0 and -a + b + c >= 0 and (a + b + c).denominator == 1)

    def delta(a, b, c):
        f = math.factorial
        return math.sqrt(f(int(a + b - c)) * f(int(a - b + c)) * f(int(-a + b + c)) / f(int(a + b + c + 1)))

    for t in ((j1, j2, j3), (j1, j5, j6), (j4, j2, j6), (j4, j5, j3)):
        if not tri(*t):
            return 0.0
    f = math.factorial
    pre = delta(j1, j2, j3) * delta(j1, j5, j6) * delta(j4, j2, j6) * delta(j4, j5, j3)
    a1, a2, a3, a4 = j1 + j2 + j3, j1 + j5 + j6, j4 + j2 + j6, j4 + j5 + j3
    b1, b2, b3 = j1 + j2 + j4 + j5, j2 + j3 + j5 + j6, j3 + j1 + j6 + j4
    total = 0.0
    for t in range(int(max(a1, a2, a3, a4)), int(min(b1, b2, b3)) + 1):
        total += (-1) ** t * f(t + 1) / (f(int(t - a1)) * f(int(t - a2)) * f(int(t - a3)) * f(int(t - a4)) *
                                         f(int(b1 - t)) * f(int(b2 - t)) * f(int(b3 - t)))
    return pre * total


def reduced2(lf, two_jf, li, two_ji, L=1):
    s = Fraction(1, 2)
    sj = six_j(Fraction(lf), Fraction(two_jf, 2), s, Fraction(two_ji, 2), Fraction(li), Fraction(L))
    return ((two_ji + 1) * (two_jf + 1) * sj ** 2 * (2 * lf + 1) * (2 * L + 1) * (2 * li + 1) / (4 * math.pi) *
            three_j_000(lf, L, li) ** 2)


def overlap(uf_in, uf_out, ui_in, ui_out):
    inner = quad(lambda r: uf_in(r) * r * ui_in(r), 0.0, R_MATCH, epsabs=0, epsrel=1e-11, limit=400)[0]
    outer = float(mp.quad(lambda r: uf_out(r) * r * ui_out(r), [R_MATCH, 30, 60, 120, mp.inf]))
    return inner + outer


def main():
    zeff = Z1 * M2 / (M1 + M2) - Z2 * M1 / (M1 + M2)
    finals = {}
    for two_jf in (3, 1):
        B, C, u_in, u_out = bound_state(1, two_jf, 0.3, 4.0)
        finals[two_jf] = (B, C, u_in, u_out)
        print(f"# bound p{two_jf}/2: B = {B:.10f} MeV, ANC = {C:.10f} fm^-1/2")
    energies = [float(x) for x in sys.argv[1:] if not x.startswith("--")] or [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0]
    print("E_cm_MeV,sigma_32_b,sigma_12_b,sigma_b,S_32_MeVb,S_12_MeVb,S_MeVb")
    for E in energies:
        waves = {(li, tji): scattering(li, tji, E) for (li, tji) in ((0, 1), (2, 3), (2, 5))}
        k, eta = waves[(0, 1)][0], waves[(0, 1)][1]
        hv = 2.0 * k / K2  # hbar v = hbar^2 k / mu (MeV fm)
        sig = {}
        for two_jf, (B, C, uf_in, uf_out) in finals.items():
            kg = (E + B) / HBARC
            total = 0.0
            for (li, tji), (_, _, ui_in, ui_out) in waves.items():
                r2 = reduced2(1, two_jf, li, tji)
                if r2 == 0.0:
                    continue
                total += r2 * overlap(uf_in, uf_out, ui_in, ui_out) ** 2
            sig[two_jf] = (16 * math.pi / 9) * kg ** 3 / hv / 2.0 * (4 * math.pi / k ** 2) * E2 * zeff ** 2 * total * 0.01
        g = math.exp(2 * math.pi * eta)
        s32, s12 = sig[3], sig[1]
        print(f"{E},{s32:.10e},{s12:.10e},{s32 + s12:.10e},{E * s32 * g:.10e},{E * s12 * g:.10e},{E * (s32 + s12) * g:.10e}",
              flush=True)


if __name__ == "__main__":
    main()
