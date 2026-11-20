#!/usr/bin/env python3
"""Reference values for tests/upgrade/interior_capture_test.cpp (docs/upgrade/interior_potential.md section 6),
independent of the engine: scipy's DOP853 with mpmath's Coulomb and Whittaker functions.

  python3 tests/upgrade/interior_capture/make_reference.py > tests/upgrade/interior_capture/reference.csv

alpha + 12C (masses 4.0026 and 12 u, as the test's project) with the Woods-Saxon of Koonin, Tombrello and Fox in every l (V0 = 45.891 MeV, R = 3.55 fm, a = 0.60 fm,
R_C = 3.55 fm), E1 in the long-wavelength form (operator r), entrance l = 1, final state l = 0 with binding energy B:
the potential's own 0+ bound state (row "bound", found here) and B = 7.16 MeV, which is not a state of the potential.
With u the regular solution, W~ the decaying one (Whittaker beyond 60 fm), O the outgoing one (G + iF beyond 60 fm),
f = a u'/u, S = a W~'/W~, L = a O'/O at the channel radius a, b the boundary constant, and

  w_in(r) = W~_f(a) (S_f - b)/(f_f - b) u_f(r)/u_f(a)      the potential part of the final state inside a,
  I_in    = int_0^a [u_i(r)/u_i(a)] w_in(r) r dr,

each row is the ratio of a capture term with the potential's interior to the same term from outside a alone:

  direct   1 + u_i(a) I_in / int_a^inf u_i W~_f r dr                      (no level: the potential's own wave)
  channel  1 + O(a) (L - b)/(f - b) I_in / int_a^inf O W~_f r dr          (a level's channel term, E > 0; complex)
  closed   1 + (S_i - b)/(f_i - b) I_in / int_a^inf [W~_i/W~_i(a)] W~_f r dr   (the same below threshold, E < 0)

and row "total" is int_0^inf u_i w r dr / (the same at a = 5.35 fm), u_i normalized beyond the potential, with w = w_in
inside a and W~_f outside: 1 at every a for the potential's bound state (b = 0).  Before the rows the script asserts
the identity behind "channel": with one level, R = 1/(f - b) + g^2/(E_l - E) and psi = I - U O, the potential's part
of psi at a, R_pot (a psi' - b psi), minus its value without the level equals -(U - U_pot) O(a) (L - b)/(f - b).
"""
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

HBARC = 197.3269804; ALPHA = 1 / 137.035999084; U = 931.49410242
MU = 4.0026 * 12.0 / 16.0026 * U  # the masses of tests/upgrade/interior_bound/interior_bound.azr
Z = 12; E2 = ALPHA * HBARC
V0, R, A, RC = 45.891, 3.55, 0.60, 3.55
K = 2 * MU / HBARC**2
RMAX, RMIN = 60.0, 1e-3


def V(r, l):
    vc = Z * E2 / r if r >= RC else Z * E2 / (2 * RC) * (3 - (r / RC) ** 2)
    return vc - V0 / (1 + np.exp((r - R) / A)) + l * (l + 1) / (K * r * r)


def solve(E, l, r0, r1, y0):
    return solve_ivp(lambda r, y: [y[1], K * (V(r, l) - E) * y[0]], (r0, r1), y0, method="DOP853", rtol=1e-13,
                     atol=1e-300, dense_output=True)


def regular(E, l):
    return solve(E, l, RMIN, RMAX, [RMIN ** (l + 1), (l + 1) * RMIN**l]).sol


def decaying(B, l):
    kappa = np.sqrt(K * B); eta = Z * E2 * MU / (HBARC**2 * kappa)
    w = lambda r: mp.whitw(-eta, l + 0.5, 2 * kappa * r)
    return solve(-B, l, RMAX, 3.0, [float(w(RMAX)), float(mp.diff(w, RMAX))]).sol


def outgoing(E, l):
    k = np.sqrt(K * E); eta = Z * E2 * MU / (HBARC**2 * k)
    g = lambda r: mp.coulombg(l, eta, k * r)
    f = lambda r: mp.coulombf(l, eta, k * r)
    G = solve(E, l, RMAX, 3.0, [float(g(RMAX)), float(mp.diff(g, RMAX))]).sol
    F = solve(E, l, RMAX, 3.0, [float(f(RMAX)), float(mp.diff(f, RMAX))]).sol
    return lambda r: (G(r)[0] + 1j * F(r)[0], G(r)[1] + 1j * F(r)[1])


def integ(fn, lo, hi):
    return quad(fn, lo, hi, epsabs=0, epsrel=1e-11, limit=800)[0]


# the potential's 0+ bound state: the regular and the decaying solution match (at 6 fm)
def mismatch(B):
    u, w = regular(-B, 0)(6.0), decaying(B, 0)(6.0)
    return u[1] / u[0] - w[1] / w[0]


B0 = brentq(mismatch, 2.3, 2.6, xtol=1e-12)

# the identity behind "channel" (one level on top of R_pot), at a = 5.35 fm, b = 0.4
a, b, E, El, g2 = 5.35, 0.4, 2.0, 3.2, 0.15
u = regular(E, 1)(a); o, do = outgoing(E, 1)(a)
f = a * u[1] / u[0]; L = a * do / o; inc, dinc = np.conj(o), np.conj(do)
Rpot = 1 / (f - b)
def surface(Rf):  # U and the potential's part of psi at a for the R-function Rf
    Uc = (inc - Rf * (a * dinc - b * inc)) / (o - Rf * (a * do - b * o))
    return Uc, Rpot * (a * (dinc - Uc * do) - b * (inc - Uc * o))
U0, chi0 = surface(Rpot); U1, chi1 = surface(Rpot + g2 / (El - E))
assert abs((chi1 - chi0) / (-(U1 - U0) * o * (L - b) / (f - b)) - 1) < 1e-10
assert abs(chi0 / (inc - U0 * o) - 1) < 1e-10  # without levels the potential's part is the whole wave

print("kind,a_fm,b,B_MeV,E_MeV,re,im")
print(f"bound,0,0,{B0:.12f},0,0,0")
for B in (B0, 7.16):
    wf = decaying(B, 0); uf = regular(-B, 0)
    for b in (0.0, -1.5):
        totals = {}
        for a in (5.35, 4.0, 6.5):
            Sf = a * wf(a)[1] / wf(a)[0]; ff = a * uf(a)[1] / uf(a)[0]
            w_in = lambda r: wf(a)[0] * (Sf - b) / (ff - b) * uf(r)[0] / uf(a)[0]
            for E in (0.5, 1.0, 2.0, 2.4, 3.0, 5.0, -1.0, -3.0):
                if E + B <= 0:
                    continue  # no photon: the entrance channel lies below the final state
                ui = regular(E, 1)
                I_in = integ(lambda r: ui(r)[0] / ui(a)[0] * w_in(r) * r, RMIN, a)
                if E > 0:
                    out = integ(lambda r: ui(r)[0] * wf(r)[0] * r, a, RMAX)
                    print(f"direct,{a},{b},{B:.12f},{E},{1 + ui(a)[0] * I_in / out:.12e},0")
                    norm = abs(ui(40.0)[0] + 1j * ui(40.0)[1] / np.sqrt(K * E))  # ~ the amplitude beyond the potential
                    totals[E] = (ui(a)[0] * I_in + out) / norm
                    if a != 5.35:
                        print(f"total,{a},{b},{B:.12f},{E},{totals[E] / first[E]:.12e},0")
                    o = outgoing(E, 1); oa, doa = o(a); L = a * doa / oa; fi = a * ui(a)[1] / ui(a)[0]
                    re = integ(lambda r: o(r)[0].real * wf(r)[0] * r, a, RMAX)
                    im = integ(lambda r: o(r)[0].imag * wf(r)[0] * r, a, RMAX)
                    c = 1 + oa * (L - b) / (fi - b) * I_in / (re + 1j * im)
                    print(f"channel,{a},{b},{B:.12f},{E},{c.real:.12e},{c.imag:.12e}")
                else:
                    wi = decaying(-E, 1); Si = a * wi(a)[1] / wi(a)[0]; fi = a * ui(a)[1] / ui(a)[0]
                    out = integ(lambda r: wi(r)[0] / wi(a)[0] * wf(r)[0] * r, a, RMAX)
                    print(f"closed,{a},{b},{B:.12f},{E},{1 + (Si - b) / (fi - b) * I_in / out:.12e},0")
            if a == 5.35:
                first = dict(totals)
