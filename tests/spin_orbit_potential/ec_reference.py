#!/usr/bin/env python3
"""The effect of the spin-orbit term on the external-capture radial integrals of models/capture.azr, by an integration
independent of the engine's (tests/pyazr/spin_orbit_potential_test.py compares it with pyazr.ec_integrals).

12C(p,g)13N into the ground state (l_f = 1, J_f = 1/2, bound by B = 1.94351 MeV) from the entrance s and d waves
(l_i = 0, J_i = 1/2 and l_i = 2, J_i = 3/2), E1, long-wavelength operator r^1; p + 12C (1.00728 u, 12 u), a_c = 3.6 fm,
Woods-Saxon 20 MeV / 3.6 fm / 0.6 fm, R_C = 3.5 fm, and either no spin-orbit term or V_so = 6 MeV / 3.6 fm / 0.6 fm on
the proton (FRESCO's convention, <2 l.i> = j(j+1) - l(l+1) - 3/4).  With H_i = G_i + i F_i the distorted outgoing
entrance wave (Coulomb H+ beyond 15 fm), phi_i = arg H_i(a) and W~_f the distorted bound state (unit ANC: the
Whittaker function beyond 15 fm), FW = int_a^inf F_i W~_f r dr and GW alike, the engine's integrals are a factor that
does not depend on the potential times
    direct:           cos(phi_i) FW - sin(phi_i) GW
    channel capture:  sqrt(P_i) exp(-i phi_i) (GW + i FW)
(src/ECIntegral.cpp, CalculateDirect).  Rows: the ratio (with term) / (without term) of each, at E_c = 0.5 and 1 MeV.
Method: scipy DOP853 (rtol 1e-12) inward from 15 fm with mpmath's Coulomb and Whittaker functions, quad on [a, 15],
mpmath beyond.  Not the engine's Numerov, quadrature or special functions.

Run: python3 ec_reference.py > expected/ec_ratio.csv   (about a minute)
"""
import math

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp

mp.mp.dps = 20
HBARC, UCONV, ALPHA = 197.3269804, 931.49410242, 1.0 / 137.035999084
MU = 1.00728 * 12.0 / 13.00728 * UCONV
Z1Z2, A_C, R_M, R_C, B = 6, 3.6, 15.0, 3.5, 1.94351


def V(r, l, ldots, vso):
    e2 = ALPHA * HBARC
    vc = Z1Z2 * e2 / (2.0 * R_C) * (3.0 - (r / R_C) ** 2) if r < R_C else Z1Z2 * e2 / r
    x = (r - 3.6) / 0.6
    v = vc - 20.0 / (1.0 + math.exp(x)) + HBARC ** 2 / (2.0 * MU) * l * (l + 1) / r ** 2
    return v + vso * 2.000 * (-math.exp(x) / (0.6 * (1.0 + math.exp(x)) ** 2)) / r * ldots


def ldots(l, two_j):
    j = two_j / 2.0
    return j * (j + 1) - l * (l + 1) - 0.75


def inward(l, two_j, E, w, vso):
    u0, du0 = complex(w(R_M)), complex(mp.diff(w, R_M))
    ld = ldots(l, two_j)
    def f(r, y):
        g = 2.0 * MU / HBARC ** 2 * (V(r, l, ld, vso) - E)
        return [y[2], y[3], g * y[0], g * y[1]]
    return solve_ivp(f, (R_M, A_C), [u0.real, u0.imag, du0.real, du0.imag], method="DOP853", rtol=1e-12,
                     atol=1e-14 * abs(u0), dense_output=True)


def amplitudes(li, two_ji, E, vso):
    k = math.sqrt(2.0 * MU * E) / HBARC
    eta = Z1Z2 * ALPHA * math.sqrt(MU / (2.0 * E))
    H = lambda r: mp.coulombg(li, eta, k * r) + 1j * mp.coulombf(li, eta, k * r)
    kappa = math.sqrt(2.0 * MU * B) / HBARC
    W = lambda r: mp.whitw(-Z1Z2 * ALPHA * MU / HBARC / kappa, 1.5, 2.0 * kappa * r)
    hi, wf = inward(li, two_ji, E, H, vso), inward(1, 1, -B, W, vso)
    def inner(part):
        return quad(lambda r: (wf.sol(r)[0] * r) * (hi.sol(r)[1] if part == "F" else hi.sol(r)[0]), A_C, R_M,
                    epsabs=0, epsrel=1e-12, limit=400)[0]
    def outer(part):
        g = (lambda r: mp.coulombf(li, eta, k * r)) if part == "F" else (lambda r: mp.coulombg(li, eta, k * r))
        return float(mp.quad(lambda r: g(r) * W(r) * r, [R_M, 25, 40, 70, 120]))
    FW, GW = inner("F") + outer("F"), inner("G") + outer("G")
    y = hi.sol(A_C)
    h, dh = complex(y[0], y[1]), complex(y[2], y[3])
    phi = math.atan2(h.imag, h.real)
    P = (k * A_C * dh / k / h).imag  # rho H'/H with H' = dH/drho
    direct = math.cos(phi) * FW - math.sin(phi) * GW
    channel = math.sqrt(P) * complex(math.cos(phi), -math.sin(phi)) * complex(GW, FW)
    return direct, channel


print("li,two_ji,E_MeV,direct_ratio,channel_ratio_re,channel_ratio_im")
for li, two_ji in ((0, 1), (2, 3)):
    for E in (0.5, 1.0):
        d1, c1 = amplitudes(li, two_ji, E, 6.0)
        d0, c0 = amplitudes(li, two_ji, E, 0.0)
        c = c1 / c0
        print(f"{li},{two_ji},{E},{d1 / d0:.12e},{c.real:.12e},{c.imag:.12e}")
