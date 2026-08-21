"""Seyler & Weller, Phys. Rev. C 20 (1979) 453 -- Eqs. (20) and (21).

A python prototype of the capture-gamma Legendre coefficients in the channel-spin
representation, checked against the worked example the paper prints on p. 458.
Its job is to pin down the one convention the paper does not state explicitly:
the normalization of Fano's X coefficient.

    sigma_u(theta)  = N sum_k a_k P_k(cos theta)
    sigma (theta,phi) = N sum_k [ a_k P_k + b_k P_k^1 p_y + ... ]
so
    A_y(theta) = sum_k b_k P_k^1(cos theta) / sum_k a_k P_k(cos theta).

Channel label t = {p L b l s}: photon mode p (1 = electric, 0 = magnetic) and
multipolarity L, compound spin b, entrance orbital angular momentum l and
channel spin s.  x is the projectile spin, a the target spin, c the residual
(final) nuclear spin.
"""
import itertools
import numpy as np
from sympy import S, sqrt as ssqrt
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j


def hat(j):
    return np.sqrt(2.0 * float(j) + 1.0)


def cg(j1, j2, j3, m1, m2, m3):
    """Clebsch-Gordan <j1 m1 j2 m2 | j3 m3>, matching AngCoeff::ClebGord."""
    w = wigner_3j(S(j1), S(j2), S(j3), S(m1), S(m2), S(-m3))
    return float(w) * (-1.0) ** float(j1 - j2 + m3) * hat(j3)


def racahW(a, b, c, d, e, f):
    """Racah W(a b c d; e f), matching AngCoeff::Racah."""
    w = wigner_6j(S(a), S(b), S(e), S(d), S(c), S(f))
    return float(w) * (-1.0) ** float(a + b + c + d)


def fanoX(a, b, c, d, e, f, g, h, i, normalized=True):
    """Fano's X coefficient.

    `normalized=True` is the usual convention,
        X = [prod (2j+1)]^{1/2} { 9j },
    `normalized=False` the bare 9-j symbol.  Which one Seyler and Weller use is
    what the worked example decides.
    """
    nine = float(wigner_9j(S(a), S(b), S(c), S(d), S(e), S(f), S(g), S(h), S(i)))
    if not normalized:
        return nine
    pref = 1.0
    for j in (a, b, c, d, e, f, g, h, i):
        pref *= hat(j)
    return pref * nine


def parity_bracket(L, p, Lp, pp, k):
    """Eq. (15): [ ] = 1/2 [1 + (-1)^{L+p+L'+p'+k}]."""
    return 0.5 * (1.0 + (-1.0) ** (L + p + Lp + pp + k))


class Channel:
    """One capture pathway t = {p, L, b, l, s}."""

    def __init__(self, p, L, b, l, s, name=""):
        self.p, self.L, self.b, self.l, self.s = p, L, b, l, s
        self.name = name

    def __repr__(self):
        return self.name or f"(p{self.p} L{self.L} b{self.b} l{self.l} s{self.s})"


def a_coeff(k, chans, R, x, a, c):
    """Eq. (20).  a_0 must come out as sum_t (2b+1)|R_t|^2."""
    tot = 0.0
    for (i, t), (j, tp) in itertools.product(enumerate(chans), repeat=2):
        br = parity_bracket(t.L, t.p, tp.L, tp.p, k)
        if br == 0.0:
            continue
        term = ((-1.0) ** float(t.s - c + 1)
                * br
                * hat(t.l) * hat(tp.l) * hat(t.L) * hat(tp.L)
                * (2.0 * t.b + 1.0) * (2.0 * tp.b + 1.0)
                * cg(t.l, tp.l, k, 0, 0, 0)
                * racahW(t.l, t.b, tp.l, tp.b, t.s, k)
                * cg(t.L, tp.L, k, 1, -1, 0)
                * racahW(t.L, t.b, tp.L, tp.b, c, k))
        if term == 0.0:
            continue
        tot += term * (R[i] * np.conj(R[j])).real
    return tot


def b_coeff(k, chans, R, x, a, c, normalized_X=True):
    """Eq. (21)."""
    if k == 0:
        return 0.0
    pref = 3.0 * np.sqrt(float(x)) * hat(x) * hat(k) / np.sqrt(
        (float(x) + 1.0) * k * (k + 1.0))
    tot = 0.0
    for (i, t), (j, tp) in itertools.product(enumerate(chans), repeat=2):
        br = parity_bracket(t.L, t.p, tp.L, tp.p, k)
        if br == 0.0:
            continue
        term = (br
                * hat(t.s) * hat(tp.s) * hat(t.l) * hat(tp.l)
                * hat(t.L) * hat(tp.L)
                * (2.0 * t.b + 1.0) * (2.0 * tp.b + 1.0)
                * (-1.0) ** float(a - x + c - t.b - t.s + t.l)
                * cg(t.l, tp.l, k, 0, 0, 0)
                * racahW(x, t.s, x, tp.s, a, 1)
                * cg(t.L, tp.L, k, 1, -1, 0)
                * racahW(t.L, t.b, tp.L, tp.b, c, k)
                * fanoX(t.l, t.s, t.b, tp.l, tp.s, tp.b, k, 1, k,
                        normalized=normalized_X))
        if term == 0.0:
            continue
        tot += term * (1j * R[i] * np.conj(R[j])).real
    return pref * tot


# --------------------------------------------------------------- the test case
# Seyler & Weller p. 457: "1/2 + 1 -> (E1 or E2) + 1/2", ignoring the s = 3/2
# channel spin so that s = s' = c = 1/2.  Parity then forces l = L.
X_SPIN, A_SPIN, C_SPIN = 0.5, 1.0, 0.5
CHANS = [
    Channel(1, 1, 0.5, 1, 0.5, "R_1,1/2"),
    Channel(1, 1, 1.5, 1, 0.5, "R_1,3/2"),
    Channel(1, 2, 1.5, 2, 0.5, "R_2,3/2"),
    Channel(1, 2, 2.5, 2, 0.5, "R_2,5/2"),
]

# The published coefficients, as {k: {(i,j): value}} with (i,j) indexing CHANS.
# A term "A R_i R_j C" is a cos(phi_i - phi_j) term, "S" a sin one.
PAPER_A = {
    0: {(0, 0): 2.0, (1, 1): 4.0, (2, 2): 4.0, (3, 3): 6.0},
    1: {(2, 0): 6.92, (2, 1): 1.38, (3, 1): 12.48},
    2: {(1, 1): -2.0, (2, 2): 2.0, (3, 3): 3.42, (1, 0): -4.0, (3, 2): 1.72},
    3: {(3, 0): -6.92, (2, 1): -8.32, (3, 1): -5.54},
    4: {(3, 3): -3.42, (3, 2): -13.72},
}
PAPER_B = {
    1: {(2, 0): -1.16, (2, 1): -0.92, (3, 1): 2.08},
    2: {(1, 0): -0.66, (3, 2): 0.47},
    3: {(3, 0): -0.77, (2, 1): 0.92, (3, 1): -0.15},
    4: {(3, 2): -1.14},
}


def extract(k, kind, normalized_X=True):
    """Pull the per-pair coefficients out by finite differences on |R| and phase.

    A term is  A_ij R_i R_j cos(phi_i - phi_j)  for a_k and
               B_ij R_i R_j sin(phi_i - phi_j)  for b_k, so setting one pair
    to unit modulus with a relative phase of 0 or pi/2 isolates it.
    """
    out = {}
    n = len(CHANS)
    for i in range(n):
        for j in range(i, n):
            R = np.zeros(n, dtype=complex)
            if i == j:
                R[i] = 1.0
                v = (a_coeff(k, CHANS, R, X_SPIN, A_SPIN, C_SPIN) if kind == "a"
                     else b_coeff(k, CHANS, R, X_SPIN, A_SPIN, C_SPIN, normalized_X))
                if abs(v) > 5e-3:
                    out[(i, i)] = v
            else:
                # cos term: phases equal; sin term: phi_i - phi_j = pi/2
                R[:] = 0
                R[i], R[j] = 1.0, 1.0
                cval = (a_coeff(k, CHANS, R, X_SPIN, A_SPIN, C_SPIN) if kind == "a"
                        else b_coeff(k, CHANS, R, X_SPIN, A_SPIN, C_SPIN, normalized_X))
                Rd = np.zeros(n, dtype=complex)
                Rd[i], Rd[j] = 1.0, 1.0
                diag = 0.0
                for m in (i, j):
                    Rs = np.zeros(n, dtype=complex); Rs[m] = 1.0
                    diag += (a_coeff(k, CHANS, Rs, X_SPIN, A_SPIN, C_SPIN) if kind == "a"
                             else b_coeff(k, CHANS, Rs, X_SPIN, A_SPIN, C_SPIN, normalized_X))
                cross_cos = cval - diag
                R[:] = 0
                R[i], R[j] = 1.0j, 1.0          # phi_i - phi_j = pi/2
                sval = (a_coeff(k, CHANS, R, X_SPIN, A_SPIN, C_SPIN) if kind == "a"
                        else b_coeff(k, CHANS, R, X_SPIN, A_SPIN, C_SPIN, normalized_X))
                cross_sin = sval - diag
                v = cross_cos if kind == "a" else cross_sin
                if abs(v) > 5e-3:
                    # a_k terms carry cos(phi_i - phi_j), which is symmetric;
                    # b_k terms carry sin(phi_i - phi_j), which is not, so the
                    # transposed key gets the opposite sign.  The paper writes
                    # each term with the larger-L channel first.
                    out[(i, j)] = v
                    out[(j, i)] = v if kind == "a" else -v
    return out


def compare(kind, paper, normalized_X=True):
    print(f"\n--- {kind}_k   (X {'normalized' if normalized_X else 'bare 9-j'}) ---")
    allok = True
    for k in sorted(paper):
        got = extract(k, kind, normalized_X)
        want = paper[k]
        line = []
        ok = True
        for key in sorted(want):
            g = got.get(key, 0.0)
            w = want[key]
            if abs(g - w) > 0.02 + 0.01 * abs(w):
                ok = False
            line.append(f"{CHANS[key[0]].name}x{CHANS[key[1]].name}: "
                        f"{g:+7.3f} (paper {w:+6.2f})")
        allok &= ok
        print(f" k={k} {'OK ' if ok else 'BAD'} " + " | ".join(line))
    return allok


if __name__ == "__main__":
    ok_a = compare("a", PAPER_A)
    for norm in (True, False):
        ok_b = compare("b", PAPER_B, normalized_X=norm)
        print(f"   -> b_k with normalized_X={norm}: {'MATCH' if ok_b else 'no match'}")
