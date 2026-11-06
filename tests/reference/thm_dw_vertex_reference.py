#!/usr/bin/env python3
"""Independent reference for the distorted-wave THM entrance vertex (src/ThmDwVertex.cpp).

The engine reduces the six-dimensional surface integral to reduced amplitudes
h^l_{Ls La} with the distorted waves in partial waves (Numerov + COUL, Lagrange
tables).  This script uses none of that.  The source

    S(r) = Int d^3u phi(u) chi(-)*_{k_sF}(alpha r + u) chi(+)_{k_aA}(r + beta u)

is integrated directly on a three-dimensional product grid for points r on
spheres around r = a, with the Coulomb wave in closed form,

    chi(+)_k(R) = exp(-pi eta/2) Gamma(1 + i eta) e^{i k z} 1F1(-i eta, 1, i k (R - z)),

1F1 from mpmath (tabulated in xi = R - z and splined), the s + F wave a plane
wave (neutral spectator), phi the Yukawa tail e^{-kappa u}/u (the Whittaker
tail with eta_b = 0), phi~(q) = 4 pi/(kappa^2 + q^2).  Axial symmetry (qf:
k_sF along k_aA = z) leaves only m = 0: S_l0(r) = 2 pi Int dcos Y_l0 S, and
S_l0'(a) by a five-point difference in r.  Then (s, d) = (S_l0(a), a
S_l0'(a))/(4 pi phi~) and G = (4 pi/(2l+1)) [[|s|^2, s* d], [d* s, |d|^2]].

Case: 19F(d,n) at 55 MeV (19F beam on the d target, d = p + n, spectator n),
point Coulomb in d + 19F, channel radius 5.136 fm; l = 0 at E = 0.3239 MeV and
l = 1 at E = 0.2127 MeV (the 324 and 213 keV resonances of 19F(p,a)).

    python3 thm_dw_vertex_reference.py          # ~15 min, two resolutions
"""
import sys
import numpy as np
import mpmath as mp
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline
from scipy.special import sph_harm_y

HBARC = 197.3269804
UCONV = 931.494
FSTRUC = 1.0 / 137.035999084
AMU = 931.49410242
kD, kP, kN, kF19 = 2.0135532134, 1.0072764675, 1.0086649159, 18.9934652

# Kinematics as the engine builds them (ThmDistortion::Setup, ThmDwVertex::Build).
ma, mA, ms, mx = kD, kF19, kN, kP
Ebeam = 55.0
bind = (mx + ms - ma) * AMU
eAA = Ebeam * kD / (kF19 + kD)            # horse (d) is the target
mu_aA = ma * mA / (ma + mA) * UCONV
mu_sF = ms * (mx + mA) / (ms + mx + mA) * UCONV
mu_sx = ms * mx / (ms + mx) * UCONV
ka = np.sqrt(2 * mu_aA * eAA) / HBARC
eta = 1 * 9 * FSTRUC * mu_aA / (HBARC * ka)
kappa = np.sqrt(2 * mu_sx * bind) / HBARC
alpha = mA / (mx + mA)
beta = ms / ma
A_RAD = 5.136


_TABLE = {}


def coulomb_table(xmax, step=0.005):
    """F(xi) = 1F1(-i eta, 1, i k xi) on [0, xmax], splined (real and imaginary parts)."""
    if xmax in _TABLE:
        return _TABLE[xmax]
    xs = np.arange(0.0, xmax + 10 * step, step)
    mp.mp.dps = 20
    vals = np.array([complex(mp.hyp1f1(-1j * eta, 1, 1j * ka * x)) for x in xs])
    _TABLE[xmax] = (CubicSpline(xs, vals.real), CubicSpline(xs, vals.imag))
    return _TABLE[xmax]


def source(rvals, thetas, ks, F, nu_panels, nth_u, nph_u, umax):
    """S(r, theta_r) for r in rvals, theta_r in thetas (r in the xz plane)."""
    norm = complex(mp.exp(-mp.pi * eta / 2) * mp.gamma(1 + 1j * eta))
    gx, gw = leggauss(16)
    edges = np.linspace(0.0, umax, nu_panels + 1)
    u = np.concatenate([0.5 * (edges[i + 1] - edges[i]) * gx + 0.5 * (edges[i + 1] + edges[i]) for i in range(nu_panels)])
    uw = np.concatenate([0.5 * (edges[i + 1] - edges[i]) * gw for i in range(nu_panels)])
    cx, cw = leggauss(nth_u)
    ph = 2 * np.pi * (np.arange(nph_u) + 0.5) / nph_u
    pw = np.full(nph_u, 2 * np.pi / nph_u)
    out = np.zeros((len(rvals), len(thetas)), complex)
    # One u panel (16 nodes) at a time keeps the arrays small.
    for p0 in range(0, len(u), 16):
        U, C, P = np.meshgrid(u[p0:p0 + 16], cx, ph, indexing="ij")
        W = (uw[p0:p0 + 16, None, None] * cw[None, :, None] * pw[None, None, :]) * U * np.exp(-kappa * U)
        S = np.sqrt(1 - C * C)
        ux, uy, uz = U * S * np.cos(P), U * S * np.sin(P), U * C
        for i, r in enumerate(rvals):
            for j, t in enumerate(thetas):
                rx, rz = r * np.sin(t), r * np.cos(t)
                # chi(-)*_{k_sF}(alpha r + u): plane wave e^{-i k_s z_s}, k_sF along z.
                zs = alpha * rz + uz
                # chi(+)_{k_aA}(r + beta u): Coulomb along z.
                Xa, Ya, Za = rx + beta * ux, beta * uy, rz + beta * uz
                xi = np.sqrt(Xa * Xa + Ya * Ya + Za * Za) - Za
                f = F[0](xi) + 1j * F[1](xi)
                out[i, j] += np.sum(W * np.exp(1j * (ka * Za - ks * zs)) * f) * norm
    return out


def gram(l, e, res):
    esf = eAA - bind - e
    ks = np.sqrt(2 * mu_sF * esf) / HBARC
    q = abs(ks - beta * ka)
    phit = 4 * np.pi / (kappa ** 2 + q ** 2)
    umax = 26.0 / kappa
    nu_panels, nth_u, nph_u, nth_r = res
    F = coulomb_table(2 * (A_RAD + 0.2 + beta * umax) + 1.0)
    h = 0.02
    rvals = A_RAD + h * np.arange(-2, 3)
    tx, tw = leggauss(nth_r)
    thetas = np.arccos(tx)
    S = source(rvals, thetas, ks, F, nu_panels, nth_u, nph_u, umax)
    Y = np.real(sph_harm_y(l, 0, thetas, 0.0))
    Sl = 2 * np.pi * (S * (tw * Y)[None, :]).sum(axis=1)
    s = Sl[2]
    d = A_RAD * (Sl[0] - 8 * Sl[1] + 8 * Sl[3] - Sl[4]) / (12 * h)
    s /= 4 * np.pi * phit
    d /= 4 * np.pi * phit
    f = 4 * np.pi / (2 * l + 1)
    return f * abs(s) ** 2, f * abs(d) ** 2, f * np.conj(s) * d


if __name__ == "__main__":
    print(f"# k_aA = {ka:.10f} fm^-1, eta_aA = {eta:.10f}, kappa = {kappa:.10f}, alpha = {alpha:.10f}, beta = {beta:.10f}")
    for l, e in ((0, 0.3239), (1, 0.2127)):
        rows = []
        for res in ((36, 96, 64, 16), (48, 128, 96, 20)):
            g11, g22, g12 = gram(l, e, res)
            rows.append((g11, g22, g12))
            print(f"# l={l} E={e} res={res}: G11={g11:.12e} G22={g22:.12e} G12={g12.real:.12e}{g12.imag:+.12e}i", flush=True)
        g11, g22, g12 = rows[-1]
        d = max(abs(rows[0][0] - g11) / g11, abs(rows[0][1] - g22) / g22, abs(rows[0][2] - g12) / abs(g12))
        print(f'    {{{l}, {e}, {g11:.12e}, {g22:.12e}, {g12.real:.12e}, {g12.imag:.12e}}},  // resolution change {d:.1e}')
