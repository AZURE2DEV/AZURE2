#!/usr/bin/env python3
"""Independent external-capture cross sections for the ec_m1_exact models (writes expected/sigma.dat, sigma_a12c.dat).

Hard-sphere external capture beyond the channel radius a, final state u_f(r) = C W_{-eta_f, l_f+1/2}(2 kappa r) with
ANC C.  Nothing here comes from AZURE2: the cross section is built from the golden rule with an energy-normalized
continuum,

  sigma = (pi/k^2) sum_J g_J |S_J|^2,   |S_J|^2 = 2 pi Gamma_J(E),
  Gamma_J = 8 pi (L+1) / (L [(2L+1)!!]^2) k_gamma^(2L+1) B_J,   B_J = |<f||O_L||i>|^2 / (2J+1),
  u_i(r) = sqrt(2 mu / (pi hbar^2 k)) [F_l(kr) cos d + G_l(kr) sin d],   d = -atan(F_l(ka)/G_l(ka)),

with one-body operators for two point clusters (cluster i at c_i r, c_1 = m2/M, c_2 = -m1/M, s_1 = 1, s_2 = (-1)^L):
  E_L, long wavelength   e e_L r^L Y_L,  e_L = Z1 (m2/M)^L + Z2 (-m1/M)^L
  E_L, finite k_gamma    e L (2L+1)!!/k^L sum_i Z_i s_i G_L(k c_i r) Y_L,   G_L(y) = int_0^y j_L(t)/t dt
                         (the charge term of the extended Siegert theorem, Friar and Fallieros, PRC 29, 1645 (1984))
  ... + convection       the transverse field's tangential part w_t(x) r grad Y_L against the convection current
                         e Z_i p_i/m_i of the relative motion; for u/r Y_l states it adds to the radial operator
                         L (2L+1)!!/k^L [-(lf(lf+1) - li(li+1))/sqrt(L(L+1))] sum_i Z_i s_i hbar/(2 m_i c) w_t(k c_i r)/r,
                         w_t(x) = [(L+1) j_{L-1}(x) - L j_{L+1}(x)]/((2L+1) sqrt(L(L+1))) - sqrt(L(L+1)) G_L(x)/x
                         (derived and checked against exact eigenstates of one Hamiltonian, where the current form of
                         the gradient part equals the charge form to 1e-9: R-matrix archive
                         12C+a_onefile/10-23-26_ec_tangential_current).  Spin-magnetization currents are not included.
  M1                     sqrt(3/(4 pi)) mu_N (g_l l + g_p s),  g_l = mu (Z1/m1^2 + Z2/m2^2) in proton masses
and e^2 = alpha hbar c, mu_N = e hbar c / (2 m_p c^2).  Reduced matrix elements in (l s)J coupling.

expected/sigma.dat (p + 16O -> 17F(0.495)): E_cm, sigma E1 LWA, E1 charge term, E1 charge + convection, M1 (barn).
expected/sigma_a12c.dat (alpha + 12C -> 16O g.s., E1 + E2): E_cm, sigma LWA, charge term, charge + convection (barn).
Run from this directory: python3 truth.py   (needs mpmath, numpy, scipy, sympy; about 20 minutes)
"""
import os
import mpmath as mp
import numpy as np
from scipy import integrate, special
from sympy.physics.wigner import wigner_3j, wigner_6j

import make_models as mm

HBARC, UCONV, ALPHA = 197.32696310, 931.4940880, 1 / 137.0359996790
MPROTON = 1.00727638 * UCONV
X, W = np.polynomial.legendre.leggauss(24)


def dfact(n):
    return float(np.prod(np.arange(n, 0, -2))) if n > 0 else 1.0


def G_L(L, y):
    """int_0^y j_L(t)/t dt (closed form for L = 1, checked against quadrature in __main__)."""
    if y <= 0:
        return 0.0
    if L == 1:
        return 0.5 * special.sici(y)[0] - np.sin(y) / (2 * y * y) + np.cos(y) / (2 * y)
    return integrate.quad(lambda t: special.spherical_jn(L, t) / t, 0, y, epsabs=0, epsrel=1e-13, limit=500)[0]


def w_t(L, x):
    ll = np.sqrt(L * (L + 1))
    return ((L + 1) * special.spherical_jn(L - 1, x) - L * special.spherical_jn(L + 1, x)) / ((2 * L + 1) * ll) \
        - ll * G_L(L, x) / x


class System:
    def __init__(self, z1, z2, m1, m2, a, B, lf, anc):
        self.z, self.m, self.a, self.B, self.lf, self.anc = (z1, z2), (m1, m2), a, B, lf, anc
        M = m1 + m2
        self.c = (m2 / M, m1 / M)
        self.mu = m1 * m2 / M * UCONV
        self.kappa = np.sqrt(2 * self.mu * B) / HBARC
        self.etaf = z1 * z2 * ALPHA * np.sqrt(self.mu / (2 * B))
        rmax = a + 40 / (2 * self.kappa) + 60
        edges = np.linspace(a, rmax, int((rmax - a) / 6) + 1)
        self.r = np.concatenate([0.5 * (hi - lo) * X + 0.5 * (hi + lo) for lo, hi in zip(edges[:-1], edges[1:])])
        self.w = np.concatenate([0.5 * (hi - lo) * W for lo, hi in zip(edges[:-1], edges[1:])])
        self.uf = anc * np.array([float(mp.whitw(-self.etaf, lf + 0.5, 2 * self.kappa * x)) for x in self.r])

    def ui(self, l, E):
        k = np.sqrt(2 * self.mu * E) / HBARC
        eta = self.z[0] * self.z[1] * ALPHA * np.sqrt(self.mu / (2 * E))
        F = lambda r: float(mp.coulombf(l, eta, k * r))
        G = lambda r: float(mp.coulombg(l, eta, k * r))
        d = -np.arctan(F(self.a) / G(self.a))
        norm = np.sqrt(2 * self.mu / (np.pi * HBARC ** 2 * k))
        return k, eta, norm * np.array([F(r) * np.cos(d) + G(r) * np.sin(d) for r in self.r])

    def radial_E(self, L, li, kg, u):
        """int u_f O_L(r) u_i dr for the three electric operators (LWA, charge term, charge + convection)."""
        s = (1, (-1) ** L)
        eL = sum(Z * si * c ** L for Z, si, c in zip(self.z, s, self.c))
        lwa = np.sum(self.w * self.uf * eL * self.r ** L * u)
        pref = L * dfact(2 * L + 1) / kg ** L
        chg = np.zeros_like(self.r)
        cnv = np.zeros_like(self.r)
        dl = self.lf * (self.lf + 1) - li * (li + 1)
        for Z, si, c, m in zip(self.z, s, self.c, self.m):
            x = kg * c * self.r
            chg += Z * si * np.array([G_L(L, xx) for xx in x])
            cnv += Z * si * HBARC / (2 * m * UCONV) * np.array([w_t(L, xx) for xx in x]) / self.r
        charge = np.sum(self.w * self.uf * pref * chg * u)
        conv = np.sum(self.w * self.uf * pref * (-dl / np.sqrt(L * (L + 1))) * cnv * u)
        return lwa, charge, charge + conv


def red_Y(lf, L, li):
    """<lf||Y_L||li>."""
    return float((-1) ** lf * np.sqrt((2 * lf + 1) * (2 * L + 1) * (2 * li + 1) / (4 * np.pi)) * wigner_3j(lf, L, li, 0, 0, 0))


def orbital(lf, li, s, Jf, J, L):
    """<(lf s)Jf||T_L(orbital)||(li s)J> / <lf||T_L||li>."""
    return float((-1) ** (lf + s + J + L) * np.sqrt((2 * J + 1) * (2 * Jf + 1)) * wigner_6j(lf, Jf, s, J, li, L))


def gamma_pref(L, kg):
    return 2 * np.pi * 8 * np.pi * (L + 1) / (L * dfact(2 * L + 1) ** 2) * kg ** (2 * L + 1)


def sigma_p16O(sysp, E):
    """p + 16O -> 17F(1/2+): E1 from l = 1 (J = 1/2-, 3/2-), M1 from l = 0 (J = 1/2+); channel spin 1/2."""
    kg = (E + sysp.B) / HBARC
    k, _, u1 = sysp.ui(1, E)
    _, _, u0 = sysp.ui(0, E)
    Jf, j1 = 0.5, 0.5
    out = []
    for rad in sysp.radial_E(1, 1, kg, u1):
        s = 0.0
        for J in (0.5, 1.5):
            red = orbital(0, 1, 0.5, Jf, J, 1) * red_Y(0, 1, 1) * rad
            s += (2 * J + 1) / (2 * j1 + 1) * gamma_pref(1, kg) * ALPHA * HBARC * red ** 2 / (2 * J + 1)
        out.append(np.pi / k ** 2 * s * 0.01)
    J, gp = 0.5, 5.5857
    s_red = float((-1) ** (0 + 0.5 + Jf + 1) * np.sqrt((2 * J + 1) * (2 * Jf + 1)) * wigner_6j(0.5, Jf, 0, J, 0.5, 1)) \
        * np.sqrt(j1 * (j1 + 1) * (2 * j1 + 1))
    muN2 = ALPHA * HBARC * (HBARC / (2 * MPROTON)) ** 2
    overlap = np.sum(sysp.w * sysp.uf * u0)
    Bm1 = 3 / (4 * np.pi) * muN2 * (gp * s_red * overlap) ** 2 / (2 * J + 1)
    out.append(np.pi / k ** 2 * (2 * J + 1) / (2 * j1 + 1) * gamma_pref(1, kg) * Bm1 * 0.01)
    return out


def sigma_a12c(sysa, E):
    """alpha + 12C -> 16O(0+): E1 from l = 1, E2 from l = 2 (spinless clusters, J = l)."""
    kg = (E + sysa.B) / HBARC
    tot = np.zeros(3)
    for L in (1, 2):
        k, _, u = sysa.ui(L, E)
        for n, rad in enumerate(sysa.radial_E(L, L, kg, u)):
            tot[n] += np.pi / k ** 2 * gamma_pref(L, kg) * ALPHA * HBARC * (red_Y(0, L, L) * rad) ** 2 * 0.01
    return list(tot)


if __name__ == "__main__":
    for L in (1, 2):
        for y in (0.01, 0.5, 3.0, 20.0):
            q = integrate.quad(lambda t: special.spherical_jn(L, t) / t, 0, y, epsabs=0, epsrel=1e-13, limit=500)[0]
            assert abs(G_L(L, y) - q) < 1e-12 * max(1.0, abs(q)), (L, y, G_L(L, y), q)
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(os.path.join(here, "expected"), exist_ok=True)
    sysp = System(1, 8, mm.MP, mm.M16, 5.0, 0.599607 - 0.4953, 0, 84.0)
    with open(os.path.join(here, "expected", "sigma.dat"), "w") as f:
        f.write("# E_cm(MeV)  sigma_E1_LWA  sigma_E1_charge  sigma_E1_charge+convection  sigma_M1   (barn)"
                "  -- truth.py, independent of AZURE2\n")
        for E in mm.ECM:
            row = sigma_p16O(sysp, E)
            line = f"{E:.4f}  " + "  ".join(f"{v:.10e}" for v in row)
            f.write(line + "\n"); print("p16O", line, flush=True)
    sysa = System(2, 6, mm.MA, mm.M12, 5.43, mm.SA, 0, 100.0)
    with open(os.path.join(here, "expected", "sigma_a12c.dat"), "w") as f:
        f.write("# E_cm(MeV)  sigma_LWA  sigma_charge  sigma_charge+convection   (barn, E1 + E2)"
                "  -- truth.py, independent of AZURE2\n")
        for E in mm.ECM_A:
            row = sigma_a12c(sysa, E)
            line = f"{E:.4f}  " + "  ".join(f"{v:.10e}" for v in row)
            f.write(line + "\n"); print("a12c", line, flush=True)
