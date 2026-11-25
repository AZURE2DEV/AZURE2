#!/usr/bin/env python3
"""E2 and M1 capture with interior capture, computed without the engine: the reference of check.sh section 11
(potential_background.md section 5.12).

The model (models/e2m1.azr.in): 3He + 4He with VD's s1/2 (67.5 MeV) and d3/2 (66.0 MeV) waves under
`background interior ... capture interior`, and the d5/2 wave (69.0 MeV) `only bound_state`; the final state is VD's
bound 0d5/2 state (B = 3.6804742644 MeV, scipy), at the threshold's E_x = 0, with its ANC C.  Its wave function: inside
a_c the potential's regular solution at E_f continued in value, C W(a) phi_f(r), outside the point-Coulomb Whittaker
function C W(r) (`only bound_state`, section 5.9).  The entrance waves:
  - s1/2 and d3/2: the potential's scattering waves from 0 to infinity (u = cos d F + sin d G beyond the potential);
  - d5/2 (the final state's own wave, point-Coulomb channel functions): the hard-sphere wave cos phi F - sin phi G
    beyond a_c (tan phi = F(a)/G(a)), zero inside (no background, no interior capture).
Pathways of different (l, J) do not interfere in the total cross section, so

    sigma_XL = 8 pi (L+1) / (L ((2L+1)!!)^2) k_g^(2L+1) / (hbar v) * 1/(2s+1) * 4 pi / k^2 * sum_{lJ} |<J_f||O_XL||J>|^2 |I|^2,

  E2: <(l_f s)J_f||e Z2 r^2 Y_2||(l s)J>, Z2 = Z1 (m2/M)^2 + Z2 (m1/M)^2 (long-wavelength), I = int u_f r^2 u_E;
  M1: orbital only (the pair's g factors are 0 in the model): O = sqrt(3/4pi) mu_N g_l l, g_l = mu m_p (Z1/m1^2 + Z2/m2^2),
      <(l s)J_f||l||(l s)J> = (-1)^(l+s+J+1) sqrt((2J+1)(2J_f+1)) {l J_f s; J l 1} sqrt(l(l+1)(2l+1)), I = int u_f u_E,
      mu_N^2 = alpha hbar c (hbar c / 2 m_p c^2)^2.
Output: E, sigma_E2 (b), sigma_M1 (b), and the pathways (E2 from s1/2, d3/2, d5/2; M1 from d3/2, d5/2) at the energies
given (MeV, c.m.; above 0.05 MeV).

Run: python3 -I make_reference_e2m1.py 4.0 0.5 1 2 3 4 5 > data/reference_e2m1.csv  (a_c first)
"""
import math
import sys
from fractions import Fraction

import mpmath as mp
from scipy.integrate import quad

import os

args = [a for a in sys.argv[1:]]
sys.argv = [sys.argv[0]]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_reference as m  # noqa: E402

m.DEPTH[(2, 5)] = 69.0
A_C = float(args[0])
ENERGIES = [float(x) for x in args[1:]]
MP_U = 938.27208816 / m.UCONV


def three_j_000(a, b, c):
    return m.three_j_000(a, b, c)


def main():
    B, C, _, _ = m.bound_state(2, 5, 1.0, 6.0)  # the normalised state; only C is used
    # the final state as the engine defines it with only bound_state: interior shape continued in value at a_c
    o = m.outward(2, 5, -B, A_C)
    W = m.whittaker(2, B)
    Wa = float(W(A_C))
    phi_a = o.y[0, -1]
    uf_in = lambda r: C * Wa * o.sol(r)[0] / phi_a
    uf_out = lambda r: C * W(r)
    zeff2 = m.Z1 * (m.M2 / (m.M1 + m.M2)) ** 2 + m.Z2 * (m.M1 / (m.M1 + m.M2)) ** 2
    gl = (m.M1 * m.M2 / (m.M1 + m.M2)) * MP_U * (m.Z1 / m.M1 ** 2 + m.Z2 / m.M2 ** 2)
    muN2 = m.E2 * (m.HBARC / (2.0 * 938.27208816)) ** 2
    s = Fraction(1, 2)
    Jf = Fraction(5, 2)
    print(f"# d5/2 bound state: B = {B!r} MeV, ANC = {C!r} fm^-1/2; a_c = {A_C} fm")
    print("E_cm_MeV,sigma_E2_b,sigma_M1_b,E2_s12,E2_d32,E2_d52,M1_d32,M1_d52")
    for E in ENERGIES:
        k = math.sqrt(m.K2 * E)
        eta = m.Z1 * m.Z2 * m.E2 * m.K2 / (2.0 * k)
        hv = 2.0 * k / m.K2
        kg = (E + B) / m.HBARC

        def full(l, two_j, power):
            _, _, ui, uo = m.scattering(l, two_j, E)
            inner = quad(lambda r: (uf_in(r) if r < A_C else float(uf_out(r))) * r ** power * ui(r), 0.0, m.R_MATCH,
                         epsabs=0, epsrel=1e-11, limit=400, points=[A_C])[0]
            outer = float(mp.quad(lambda r: uf_out(r) * r ** power * uo(r), [m.R_MATCH, 30, 60, 120, mp.inf]))
            return inner + outer

        # the d5/2 hard-sphere wave beyond a_c (point Coulomb)
        rho = k * A_C
        Fa, Ga = float(mp.coulombf(2, eta, rho)), float(mp.coulombg(2, eta, rho))
        phi = math.atan2(Fa, Ga)
        X = lambda r: math.cos(phi) * mp.coulombf(2, eta, k * r) - math.sin(phi) * mp.coulombg(2, eta, k * r)

        def hs(power):
            return float(mp.quad(lambda r: uf_out(r) * r ** power * X(r), [A_C, 15, 30, 60, 120, mp.inf]))

        # E2
        e2 = []
        for (l, two_j), integral in (((0, 1), full(0, 1, 2)), ((2, 3), full(2, 3, 2)), ((2, 5), hs(2))):
            J = Fraction(two_j, 2)
            sj = m.six_j(Fraction(2), Jf, s, J, Fraction(l), Fraction(2))
            red2 = (two_j + 1) * 6 * sj ** 2 * 5 * 5 * (2 * l + 1) / (4 * math.pi) * three_j_000(2, 2, l) ** 2
            e2.append(8 * math.pi * 3 / (2 * 15 ** 2) * kg ** 5 / hv / 2.0 * (4 * math.pi / k ** 2) * 0.01 *
                      m.E2 * zeff2 ** 2 * red2 * integral ** 2)
        # M1 (orbital)
        m1 = []
        for (l, two_j), integral in (((2, 3), full(2, 3, 0)), ((2, 5), hs(0))):
            J = Fraction(two_j, 2)
            sj = m.six_j(Fraction(l), Jf, s, J, Fraction(l), Fraction(1))
            red2 = (two_j + 1) * 6 * sj ** 2 * l * (l + 1) * (2 * l + 1) * 3 / (4 * math.pi)
            m1.append(16 * math.pi / 9 * kg ** 3 / hv / 2.0 * (4 * math.pi / k ** 2) * 0.01 * muN2 * gl ** 2 * red2 *
                      integral ** 2)
        print(f"{E},{sum(e2):.10e},{sum(m1):.10e}," + ",".join(f"{v:.10e}" for v in e2 + m1), flush=True)


if __name__ == "__main__":
    main()
