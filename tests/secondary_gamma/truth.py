"""Independent secondary-gamma angular distributions for A(a,b)B, B -> C + gamma (used by make_models.py).

Brune and deBoer, PRC 102, 024628 (2020).  The populations of the magnetic substates of B after the direction of b is
integrated over are built by brute force from the scattering amplitude, Eq. (17), with the exit channel spin coupled
light particle first, (I_b I_B m_b m_B | s' nu') -- AZURE2's convention (PolarizationFunc.cpp).  The photon angular
distribution of a pure multipole L from those populations then follows from |X_Lq|^2 (Jackson, Eq. 9.170 et seq.),
without Eq. (22) or the radiation parameters R_k; multipole mixing goes through Eqs. (4) and (7) of the paper.

T is a dict {(J, s_out, l_out, s_in, l_in): complex}; all T share J's energy dependence in the single-level models.
Requires numpy and sympy (make_models.py only; check.sh compares against the committed expected/ files).
"""
import itertools
import numpy as np
from sympy import S
from sympy.physics.wigner import clebsch_gordan as _cg, racah as _racah
from scipy.special import sph_harm_y


def cg(j1, j2, j3, m1, m2, m3):
    return float(_cg(S(j1), S(j2), S(j3), S(m1), S(m2), S(m3)))


def W(a, b, c, d, e, f):
    return float(_racah(S(a), S(b), S(c), S(d), S(e), S(f)))


def projections(j):
    j = S(j)
    return [j - i for i in range(int(2 * j) + 1)]


def populations(T, Ia, IA, Ib, IB):
    """P(m_B) proportional to sum_{m_b m_a m_A} int dOmega_b |M_{m_b m_B: m_a m_A}|^2 (orthonormal Y_{l'm'})."""
    P = {mB: 0.0 for mB in projections(IB)}
    for ma, mA, mb, mB in itertools.product(projections(Ia), projections(IA), projections(Ib), projections(IB)):
        coef = {}
        for (J, sp, lp, s, l), t in T.items():
            J, sp, s = S(J), S(sp), S(s)
            nu = ma + mA
            if abs(nu) > s or abs(nu) > J:
                continue
            c_in = cg(Ia, IA, s, ma, mA, nu) * np.sqrt(2 * l + 1) * cg(l, s, J, 0, nu, nu)
            nup = mb + mB
            mp = nu - nup
            if abs(nup) > sp or abs(mp) > lp:
                continue
            c_out = cg(Ib, IB, sp, mb, mB, nup) * cg(lp, sp, J, mp, nup, nu)
            coef[(lp, mp)] = coef.get((lp, mp), 0) + c_in * c_out * t
        P[mB] += sum(abs(v) ** 2 for v in coef.values())
    return P


def _x2(L, q, th):
    """|X_Lq(theta)|^2 (phi-independent), Jackson Eq. (9.170) and following."""
    def y2(m):
        return np.abs(sph_harm_y(L, m, th, 0.0)) ** 2 if abs(m) <= L else 0.0 * th
    return ((L - q) * (L + q + 1) * y2(q + 1) + (L + q) * (L - q + 1) * y2(q - 1) + 2 * q * q * y2(q)) / (2 * L * (L + 1))


def gamma_distribution_pure(P, IB, IC, L, th):
    """W(theta_gamma) for a pure multipole L, normalised so that its angle average is 1."""
    w = np.zeros_like(th)
    for mB, p in P.items():
        for mC in projections(IC):
            q = mB - mC
            if abs(q) > L:
                continue
            w += p * cg(IC, L, IB, mC, q, mB) ** 2 * _x2(L, int(q), th)
    # each |X_Lq|^2 integrates to 1 and sum_{m_C q} (I_C m_C L q | I_B m_B)^2 = 1
    return w * 4 * np.pi / sum(P.values())


def Rk(k, L, Lp, IB, IC):
    """Eq. (4)."""
    return (np.sqrt((2 * IB + 1) * (2 * L + 1) * (2 * Lp + 1)) * (-1) ** (IB - IC + L - Lp + k + 1)
            * cg(Lp, L, k, 1, -1, 0) * W(L, Lp, IB, IB, k, IC))


def gamma_distribution_mixed(P, IB, IC, L, delta, th):
    """Eq. (14) with the R_k of Eq. (7) (L and L+1 mixed with ratio delta), angle average 1."""
    tot = sum(P.values())
    w = np.zeros_like(th)
    for k in range(0, 2 * L + 3, 2):
        tk0 = np.sqrt(2 * k + 1) * sum(p * cg(IB, k, IB, mB, 0, mB) for mB, p in P.items()) / tot
        r = (Rk(k, L, L, IB, IC) + 2 * delta * Rk(k, L, L + 1, IB, IC)
             + delta ** 2 * Rk(k, L + 1, L + 1, IB, IC)) / (1 + delta ** 2)
        w += tk0 * r * np.polynomial.legendre.legval(np.cos(th), [0] * k + [1])
    return w
