#!/usr/bin/env python3
"""Independent external-capture cross sections for the ec_m1_exact models (writes expected/sigma.dat).

p + 16O -> 17F(0.495 MeV, 1/2+) + gamma, hard-sphere external capture beyond the channel radius a = 5 fm, final
state u_f(r) = C W_{-eta_f, 1/2}(2 kappa r) with ANC C = 84 fm^-1/2.  Nothing here comes from AZURE2: the cross
section is built from the golden rule with an energy-normalized continuum,

  sigma = (pi/k^2) sum_J g_J |S_J|^2,   |S_J|^2 = 2 pi Gamma_J(E),
  Gamma_J = 8 pi (L+1) / (L [(2L+1)!!]^2) k_gamma^(2L+1) B_J,   B_J = |<f||O_L||i>|^2 / (2J+1),
  u_i(r) = sqrt(2 mu / (pi hbar^2 k)) [F_l(kr) cos d + G_l(kr) sin d],   d = -atan(F_l(ka)/G_l(ka)),

with the textbook one-body operators:
  E1 (long wavelength)  e e_1 r Y_1,  e_1 = Z1 m2/M - Z2 m1/M
  E1 (finite k_gamma)   e 3/k_gamma [Z1 g(k_gamma r m2/M) - Z2 g(k_gamma r m1/M)] Y_1,  g(y) = int_0^y j_1(t)/t dt
                        (the charge term of the extended Siegert theorem, Friar and Fallieros, PRC 29, 1645 (1984))
  M1                    sqrt(3/(4 pi)) mu_N (g_l l + g_p s),  g_l = mu (Z1/m1^2 + Z2/m2^2) in proton masses
and e^2 = alpha hbar c, mu_N = e hbar c / (2 m_p c^2).  Reduced matrix elements in (l s)J coupling.

Columns of expected/sigma.dat: E_cm (MeV), sigma E1 LWA, sigma E1 exact, sigma M1 (barn).
Run from this directory: python3 truth.py   (needs mpmath, numpy, scipy, sympy; about a minute)
"""
import os
import mpmath as mp
import numpy as np
from scipy import integrate, special
from sympy.physics.wigner import wigner_3j, wigner_6j

import make_models as mm

HBARC, UCONV, ALPHA = 197.32696310, 931.4940880, 1 / 137.0359996790
MPROTON = 1.00727638 * UCONV
Z1, Z2, M1, M2 = 1, 8, mm.MP, mm.M16
J1, GP = 0.5, 5.5857
A, SEP, EF, ANC = 5.0, 0.599607, 0.4953, 84.0
M = M1 + M2
MU = M1 * M2 / M * UCONV                         # MeV
B = SEP - EF
KAPPA = np.sqrt(2 * MU * B) / HBARC
ETAF = Z1 * Z2 * ALPHA * np.sqrt(MU / (2 * B))
X, W = np.polynomial.legendre.leggauss(24)
RMAX = A + 40 / (2 * KAPPA) + 60
EDGES = np.linspace(A, RMAX, int((RMAX - A) / 6) + 1)
RGRID = np.concatenate([0.5 * (hi - lo) * X + 0.5 * (hi + lo) for lo, hi in zip(EDGES[:-1], EDGES[1:])])
WGRID = np.concatenate([0.5 * (hi - lo) * W for lo, hi in zip(EDGES[:-1], EDGES[1:])])
UF = ANC * np.array([float(mp.whitw(-ETAF, 0.5, 2 * KAPPA * r)) for r in RGRID])     # lf = 0


def g1(y):
    """int_0^y j_1(t)/t dt in closed form (checked against quadrature in __main__)."""
    return 0.5 * special.sici(y)[0] - np.sin(y) / (2 * y * y) + np.cos(y) / (2 * y)


def ui(l, E):
    k = np.sqrt(2 * MU * E) / HBARC
    eta = Z1 * Z2 * ALPHA * np.sqrt(MU / (2 * E))
    F = lambda r: float(mp.coulombf(l, eta, k * r))
    G = lambda r: float(mp.coulombg(l, eta, k * r))
    d = -np.arctan(F(A) / G(A))
    norm = np.sqrt(2 * MU / (np.pi * HBARC ** 2 * k))
    return k, norm * np.array([F(r) * np.cos(d) + G(r) * np.sin(d) for r in RGRID])


def red_Y(lf, L, li):
    """<lf||Y_L||li>."""
    return float((-1) ** lf * np.sqrt((2 * lf + 1) * (2 * L + 1) * (2 * li + 1) / (4 * np.pi)) * wigner_3j(lf, L, li, 0, 0, 0))


def orbital(lf, li, s, Jf, J, L):
    """<(lf s)Jf||T_L(orbital)||(li s)J> / <lf||T_L||li>."""
    return float((-1) ** (lf + s + J + L) * np.sqrt((2 * J + 1) * (2 * Jf + 1)) * wigner_6j(lf, Jf, s, J, li, L))


def sigma(E):
    kg = (E + B) / HBARC
    k, u1 = ui(1, E)
    _, u0 = ui(0, E)
    Jf = 0.5
    # E1: l = 1 -> lf = 0, s = 1/2, J = 1/2-, 3/2-
    e1 = Z1 * M2 / M - Z2 * M1 / M
    rad_lwa = np.sum(WGRID * UF * e1 * RGRID * u1)
    op = np.array([3 / kg * (Z1 * g1(kg * r * M2 / M) - Z2 * g1(kg * r * M1 / M)) for r in RGRID])
    rad_ex = np.sum(WGRID * UF * op * u1)
    pref = 2 * np.pi * (16 * np.pi / 9) * kg ** 3 * ALPHA * HBARC
    out = []
    for rad in (rad_lwa, rad_ex):
        s = 0.0
        for J in (0.5, 1.5):
            red = orbital(0, 1, 0.5, Jf, J, 1) * red_Y(0, 1, 1) * rad
            s += (2 * J + 1) / (2 * J1 + 1) * pref * red ** 2 / (2 * J + 1)
        out.append(np.pi / k ** 2 * s * 0.01)                       # fm^2 -> b
    # M1: l = 0 -> lf = 0, J = 1/2+, spin term only (l = 0)
    J = 0.5
    s_red = float((-1) ** (0 + 0.5 + Jf + 1) * np.sqrt((2 * J + 1) * (2 * Jf + 1)) * wigner_6j(0.5, Jf, 0, J, 0.5, 1)) \
        * np.sqrt(J1 * (J1 + 1) * (2 * J1 + 1))
    muN2 = ALPHA * HBARC * (HBARC / (2 * MPROTON)) ** 2
    overlap = np.sum(WGRID * UF * u0)
    Bm1 = 3 / (4 * np.pi) * muN2 * (GP * s_red * overlap) ** 2 / (2 * J + 1)
    out.append(np.pi / k ** 2 * (2 * J + 1) / (2 * J1 + 1) * 2 * np.pi * (16 * np.pi / 9) * kg ** 3 * Bm1 * 0.01)
    return out


if __name__ == "__main__":
    for y in (0.01, 0.5, 3.0, 20.0):
        q = integrate.quad(lambda t: special.spherical_jn(1, t) / t, 0, y, epsabs=0, epsrel=1e-13, limit=500)[0]
        assert abs(g1(y) - q) < 1e-12 * max(1.0, abs(q)), (y, g1(y), q)
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(os.path.join(here, "expected"), exist_ok=True)
    with open(os.path.join(here, "expected", "sigma.dat"), "w") as f:
        f.write("# E_cm(MeV)  sigma_E1_LWA(b)  sigma_E1_exact(b)  sigma_M1(b)   -- truth.py, independent of AZURE2\n")
        for E in mm.ECM:
            row = sigma(E)
            f.write(f"{E:.4f}  {row[0]:.10e}  {row[1]:.10e}  {row[2]:.10e}\n")
            print(f"{E:.4f}  {row[0]:.10e}  {row[1]:.10e}  {row[2]:.10e}", flush=True)
