#!/usr/bin/env python3
"""Independent reference for the THM distortion factor (src/ThmDistortion.cpp).

The zero-range prior-form DWBA amplitude of A + a(x+s) -> s + F(x+A)
(Mukhamedzhanov & Pang, PRC 99 (2019) 064618, eqs. 20-24; Mukhamedzhanov,
arXiv:2609.04498, eqs. 22-30), with r_sF = r_sx = r and r_aA = beta r,
beta = m_s/m_a:

    M = Int d^3r chi(-)*_{k_sF}(r) phi_sx(r) chi(+)_{k_aA}(beta r)
      = 4 pi/(k_sF beta k_aA) sum_l (2l+1) e^{i(sigma_l^aA + sigma_l^sF)} P_l(x)
        Int_rmin^inf dr phi(r) u_l^sF(k_sF r) u_l^aA(k_aA beta r),

phi = 2 kappa e^{-kappa r} U(1 + eta_b, 2, 2 kappa r) (the Whittaker tail
W_{-eta_b,1/2}(2 kappa r)/r) or e^{-kappa r}/r, and the plane-wave limit
M_PW = 4 pi Int r^2 j_0(q r) phi dr, q = |k_sF - beta k_aA|.

Point Coulomb: u_l = F_l(eta, rho) from mpmath.coulombf, the radial
integrals by mpmath.quad on the real axis (they converge absolutely), the
phases from mpmath.loggamma.  Nothing is shared with the C++ code, which
integrates the radial equation by Numerov's method, matches to COUL and
integrates by Simpson's rule.

Optical potential (case "optical"): u_l from scipy's DOP853 (rtol 1e-12)
started from the power series at r = 1e-3 fm, matched to mpmath's F_l, G_l
beyond the potential and the turning point, and integrated by scipy.quad on
the dense output.  Different integrator, different Coulomb functions,
different quadrature from the engine.

Run:  python3 tests/reference/thm_distortion_reference.py [case ...]
It prints the table hard-coded in thm_distortion_test.cpp (a few minutes).
"""
import sys
import mpmath as mp
import numpy as np

mp.mp.dps = 20

# include/Constants.h, and the amu of EData::BuildThmGroups (B from masses).
HBARC = 197.32696310
UCONV = 931.4940880
FSTRUC = 1 / 137.0359996790
AMU = 931.49410242
# ThmNuclide table (nuclear masses, u).
MASS = {"n": 1.0086649159, "p": 1.0072764675, "d": 2.0135532134, "3He": 3.0149322434,
        "12C": 11.9967096429, "14N": 13.9992355671, "18O": 17.9947732059}
CHARGE = {"n": 0, "p": 1, "d": 1, "3He": 2, "12C": 6, "14N": 7, "18O": 8}

# name: beam, target, spectator, x, Ebeam, horse ('beam'|'target'), angle
# ('qf' or c.m. degrees), bound ('whittaker'|'yukawa'), rmin, energies,
# optical (None or (aA 10 numbers, sF 10 numbers)).
WS_AA = (50.0, 5.5, 0.6, 20.0, 5.6, 0.6, 0.0, 0.0, 0.0, 5.5)
WS_SF = (90.0, 3.4, 0.75, 0.0, 0.0, 0.0, 10.0, 3.8, 0.65, 3.8)
CASES = {
    "c12_qf": ("14N", "12C", "d", "12C", 30.0, "beam", "qf", "whittaker", 0.0, (0.8, 1.6, 2.6), None),
    "c12_cm90": ("14N", "12C", "d", "12C", 30.0, "beam", 90.0, "yukawa", 0.0, (1.2,), None),
    "c12_rmin": ("14N", "12C", "d", "12C", 30.0, "beam", "qf", "whittaker", 3.0, (1.2,), None),
    "he3_qf": ("18O", "3He", "d", "p", 115.0, "target", "qf", "whittaker", 0.0, (0.6,), None),
    "he3_cm60": ("18O", "3He", "d", "p", 115.0, "target", 60.0, "whittaker", 0.0, (0.6,), None),
    "optical": ("14N", "12C", "d", "12C", 30.0, "beam", "qf", "whittaker", 0.0, (1.0, 2.5), (WS_AA, WS_SF)),
}


def setup(case):
    beam, target, spec, xn, ebeam, horse, angle, bound, rmin, energies, optical = CASES[case]
    a = beam if horse == "beam" else target
    A = target if horse == "beam" else beam
    ma, mA, ms, mx = MASS[a], MASS[A], MASS[spec], MASS[xn]
    Za, ZA, Zs, Zx = CHARGE[a], CHARGE[A], CHARGE[spec], CHARGE[xn]
    bind = (mx + ms - ma) * AMU
    eAA = ebeam * MASS[target] / (MASS[beam] + MASS[target])
    muAA = ma * mA / (ma + mA) * UCONV
    kAA = np.sqrt(2 * muAA * eAA) / HBARC
    etaAA = Za * ZA * FSTRUC * muAA / (HBARC * kAA)
    muSF = ms * (mx + mA) / (ms + mx + mA) * UCONV
    muSx = ms * mx / (ms + mx) * UCONV
    kappa = np.sqrt(2 * muSx * bind) / HBARC
    etaB = Zs * Zx * FSTRUC * muSx / (HBARC * kappa)
    beta = ms / ma
    return dict(eAA=eAA, bind=bind, muAA=muAA, kAA=kAA, etaAA=etaAA, muSF=muSF, ZsF=(Zs, Zx + ZA),
                ZaA=(Za, ZA), kappa=kappa, etaB=etaB, beta=beta, horse=horse, angle=angle,
                bound=bound, rmin=rmin, energies=energies, optical=optical)


def phi(r, s):
    r = mp.mpf(r)
    if s["bound"] == "yukawa" or s["etaB"] == 0:
        return mp.exp(-s["kappa"] * r) / r
    return 2 * s["kappa"] * mp.exp(-s["kappa"] * r) * mp.hyperu(1 + s["etaB"], 2, 2 * s["kappa"] * r)


def sigma(l, eta):
    return mp.im(mp.loggamma(l + 1 + 1j * eta)) if eta != 0 else mp.mpf(0)


def optical_wave(l, k, eta, mu, Z1Z2, p, rstore):
    """u_l normalized to F + T H+ on the points rstore (scipy DOP853)."""
    from scipy.integrate import solve_ivp
    fac = 2 * mu / HBARC**2
    zz = Z1Z2 * FSTRUC * HBARC
    V, R, a, W, RW, aW, WD, RD, aD, RC = p

    def ws(r, R, a):
        return 1 / (1 + np.exp((r - R) / a)) if a > 0 else 0.0

    def U(r):
        vc = (zz / r if r >= RC else zz * (3 - r * r / RC**2) / (2 * RC)) if RC > 0 else zz / r
        v = vc - V * ws(r, R, a) - 1j * W * ws(r, RW, aW)
        if WD:
            x = (r - RD) / aD
            v -= 1j * WD * 4 * np.exp(x) / (1 + np.exp(x))**2
        return fac * v

    def rhs(r, y):
        u = y[0] + 1j * y[1]
        du = y[2] + 1j * y[3]
        d2 = (l * (l + 1) / r**2 + U(r) - k * k) * u
        return [du.real, du.imag, d2.real, d2.imag]

    r0 = 1e-3
    # series u = r^{l+1}(1 + a1 r + a2 r^2), RC > 0: no 1/r term
    w0 = fac * (zz * 1.5 / RC - V * ws(0, R, a) - 1j * W * ws(0, RW, aW)
                - (1j * WD * 4 * np.exp(-RD / aD) / (1 + np.exp(-RD / aD))**2 if WD else 0))
    a2 = (w0 - k * k) / (4 * l + 6)
    u0 = r0**(l + 1) * (1 + a2 * r0**2)
    du0 = (l + 1) * r0**l * (1 + a2 * r0**2) + r0**(l + 1) * 2 * a2 * r0
    rtp = (eta + np.sqrt(eta**2 + l * (l + 1))) / k
    r1 = max(rstore[-1], R + 12 * a, RC, rtp) + max(np.pi / k, 1.0)
    r2 = r1 - 0.5 * np.pi / k
    sol = solve_ivp(rhs, (r0, r1), [u0.real, u0.imag, du0.real, du0.imag], method="DOP853",
                    rtol=1e-12, atol=1e-300, dense_output=True)
    y1, y2 = sol.sol(r1), sol.sol(r2)
    v1, v2 = y1[0] + 1j * y1[1], y2[0] + 1j * y2[1]
    F1, G1 = float(mp.coulombf(l, eta, k * r1)), float(mp.coulombg(l, eta, k * r1))
    F2, G2 = float(mp.coulombf(l, eta, k * r2)), float(mp.coulombg(l, eta, k * r2))
    det = F1 * G2 - F2 * G1
    cf = (v1 * G2 - v2 * G1) / det
    cg = (F1 * v2 - F2 * v1) / det
    norm = cf - 1j * cg

    def u(r):
        if r <= r0:
            return (r / r0)**(l + 1) * u0 / norm
        y = sol.sol(r)
        return (y[0] + 1j * y[1]) / norm
    return u


def amplitude(case, E):
    s = setup(case)
    esf = s["eAA"] - s["bind"] - E
    ksf = np.sqrt(2 * s["muSF"] * esf) / HBARC
    etasf = s["ZsF"][0] * s["ZsF"][1] * FSTRUC * s["muSF"] / (HBARC * ksf)
    kb = s["beta"] * s["kAA"]
    if s["angle"] == "qf":
        x = 1.0
    else:
        x = (1 if s["horse"] == "beam" else -1) * np.cos(np.radians(s["angle"]))
    q = np.sqrt(max(0.0, ksf**2 + kb**2 - 2 * ksf * kb * x))
    rmin, kap = s["rmin"], s["kappa"]
    pts = [rmin] + [rmin + t / kap for t in (1, 3, 6, 12, 24)] + [mp.inf]
    mpw = 4 * mp.pi * mp.quad(lambda r: r * r * (mp.sin(q * r) / (q * r) if q * r > 1e-12 else 1) * phi(r, s), pts)
    total = mp.mpc(0)
    small = 0
    l = 0
    while True:
        if s["optical"] is None:
            def f(r, l=l):
                return phi(r, s) * mp.coulombf(l, etasf, ksf * r) * mp.coulombf(l, s["etaAA"], kb * r)
            I = mp.quad(f, pts)
        else:
            from scipy.integrate import quad
            rend = rmin + 50 / kap
            usf = optical_wave(l, ksf, etasf, s["muSF"], s["ZsF"][0] * s["ZsF"][1], s["optical"][1], [rend])
            uaa = optical_wave(l, s["kAA"], s["etaAA"], s["muAA"], s["ZaA"][0] * s["ZaA"][1],
                               s["optical"][0], [s["beta"] * rend])
            ph = lambda r: float(phi(r, s))
            g = lambda r: ph(r) * usf(r) * uaa(s["beta"] * r)
            brk = [rmin + t / kap for t in (0.5, 1, 2, 4, 8, 16, 32)] + [rend]
            I = 0j
            lo = max(rmin, 1e-12)
            for hi in brk:
                re = quad(lambda r: g(r).real, lo, hi, epsabs=0, epsrel=1e-12, limit=400)[0]
                im = quad(lambda r: g(r).imag, lo, hi, epsabs=0, epsrel=1e-12, limit=400)[0]
                I += re + 1j * im
                lo = hi
            I = mp.mpc(I)
        term = (2 * l + 1) * mp.expj(sigma(l, s["etaAA"]) + sigma(l, etasf)) * mp.legendre(l, x) * I
        total += term
        small = small + 1 if abs(term) < 1e-14 * abs(total) else 0
        if l >= 5 and small >= 3:
            break
        l += 1
    M = 4 * mp.pi / (ksf * kb) * total
    return float(abs(M)**2), float(mpw), l


if __name__ == "__main__":
    names = sys.argv[1:] or list(CASES)
    print("// case, E (MeV), |M|^2, M_PW, highest l  -- thm_distortion_reference.py")
    for name in names:
        for E in CASES[name][9]:
            m2, mpw, l = amplitude(name, E)
            print(f'    {{"{name}", {E}, {m2:.12e}, {mpw:.12e}}},  // l <= {l}', flush=True)
