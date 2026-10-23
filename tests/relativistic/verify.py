"""Independent checks of AZURE2's relativistic kinematics (tests/relativistic/check.sh runs it).

usage: verify.py <work dir holding the runs> -> prints ok/FAIL lines, exit status = number of failures

  (a) na_hs[_rel]: n + 4He s-wave hard sphere, sigma = 4 pi sin^2(k a)/k^2 with k non-relativistic or from the
      invariant mass, at the c.m. energies AZURE2 reports, which must be the lab grid converted the same way.
  (b) pa_hs[_rel]: p + 4He s-wave hard sphere plus Coulomb at 90 and 150 deg, |f_C + f_N|^2 from mpmath's
      Coulomb functions at eta and rho of the same kinematics.
  (c) he5[_rel]: t(d,n)a and 4He(n,d)t at equal excitation energy obey detailed balance,
      g_dt k_dt^2 sigma(dt -> na) = g_na k_na^2 sigma(na -> dt), with k of the run's kinematics.
  (d) relativistic.out: E_s = S0 + T0 + T0^2/(2 m0) and the width factor 1 + T0/m0.
  (e) he5_ang_rel, pa_ang_rel: lab-frame differential data (value 1) -- the c.m. angle and the c.m. value
      dOmega_lab/dOmega_cm against an explicit four-vector boost of the detected particle.
"""
import os
import sys

import numpy as np

try:
    import mpmath
except ImportError:
    mpmath = None

U, HBARC, ALPHA = 931.4940880, 197.32696310, 1.0 / 137.0359996790
M = dict(n=1.00866491595, p=1.00782503207, d=2.01410177812, t=3.0160492779, a=4.00260325413)
work = sys.argv[1]
fails = 0


def check(cond, what):
    global fails
    print(("  ok    " if cond else "  FAIL  ") + what)
    fails += 0 if cond else 1


def k2(m1, m2, T, rel):
    a, b = m1 * U, m2 * U
    m = a + b
    if not rel:
        return 2.0 * a * b / m * T / HBARC ** 2
    return T * (T + 2 * m) * (T + 2 * a) * (T + 2 * b) / (4 * (m + T) ** 2) / HBARC ** 2


def eta(z, m1, m2, T, rel):
    a, b = m1 * U, m2 * U
    m = a + b
    if not rel:
        return z * ALPHA * np.sqrt(a * b / m / (2 * T))
    p = np.sqrt(k2(m1, m2, T, True)) * HBARC
    return z * ALPHA * (2 * a * b + 2 * m * T + T * T) / (2 * p * (m + T))


def cm_from_lab(m1, m2, e, rel):
    a, b = m1 * U, m2 * U
    m = a + b
    return e * m2 / (m1 + m2) if not rel else 2 * b * e / (np.sqrt(m * m + 2 * b * e) + m)


def load(run, name):
    return np.loadtxt(os.path.join(work, run, "output", name), ndmin=2)


# (a) neutral hard sphere
for rel in (False, True):
    run = "na_hs" + ("_rel" if rel else "")
    d = load(run, "AZUREOut_aa=1_R=1.extrap")
    elab = 0.5 + 0.5 * np.arange(d.shape[0])
    T = d[:, 0]
    dT = np.max(np.abs(T - cm_from_lab(M['n'], M['a'], elab, rel)) / T)
    check(dT < 1e-9, f"(a) {run}: c.m. energies = the lab grid converted ({'invariant mass' if rel else 'm2/m'}), rel {dT:.1e}")
    k = np.sqrt(k2(M['n'], M['a'], T, rel))
    sig = 4 * np.pi * np.sin(k * 3.0) ** 2 / k ** 2 / 100.0      # b
    r = np.max(np.abs(d[:, 3] - sig) / sig)
    check(r < 1e-8, f"(a) {run}: sigma = 4 pi sin^2(ka)/k^2 to {r:.1e} (E_cm {T[0]:.3f}-{T[-1]:.3f} MeV)")
na, nar = load("na_hs", "AZUREOut_aa=1_R=1.extrap"), load("na_hs_rel", "AZUREOut_aa=1_R=1.extrap")
kr = np.sqrt(k2(M['n'], M['a'], nar[-1, 0], True)) / np.sqrt(k2(M['n'], M['a'], na[-1, 0], False))
print(f"        at E_lab 30 MeV: E_cm {1e3 * (nar[-1, 0] - na[-1, 0]):+.2f} keV, k {kr - 1:+.3%} (relativistic - NR)")

# (b) charged hard sphere, Coulomb functions
if mpmath is None:
    print("  skip  (b) needs mpmath")
else:
    for rel in (False, True):
        run = "pa_hs" + ("_rel" if rel else "")
        d = load(run, "AZUREOut_aa=1_R=1.extrap")
        worst = 0.0
        for T, ang, s_az in d[:, [0, 2, 3]][::5]:
            k = np.sqrt(k2(M['p'], M['a'], T, rel))
            e = eta(2, M['p'], M['a'], T, rel)
            rho = k * 3.0
            F, G = float(mpmath.coulombf(0, e, rho)), float(mpmath.coulombg(0, e, rho))
            phi = np.arctan2(F, G)                                 # hard-sphere phase
            sigma0 = float(mpmath.arg(mpmath.gamma(1 + 1j * e)))
            th = np.radians(ang)
            s2 = np.sin(th / 2) ** 2
            fc = -e / (2 * k * s2) * np.exp(-1j * e * np.log(s2) + 2j * sigma0)
            fn = np.exp(2j * sigma0) * (np.exp(-2j * phi) - 1) / (2j * k)
            s_ref = abs(fc + fn) ** 2 / 100.0                      # b/sr
            worst = max(worst, abs(s_az - s_ref) / s_ref)
        check(worst < 2e-6, f"(b) {run}: dsigma/dOmega = |f_C + f_N|^2 (mpmath, eta and k {'invariant' if rel else 'NR'}) "
                            f"to {worst:.1e}")

# (c) detailed balance across partitions
for rel in (False, True):
    run = "he5" + ("_rel" if rel else "")
    f, b = load(run, "AZUREOut_aa=2_R=1.out"), load(run, "AZUREOut_aa=1_R=2.out")
    dex = np.max(np.abs(f[:, 1] - b[:, 1]))
    check(dex < 1e-9, f"(c) {run}: forward and reverse points at the same E_x (max diff {dex:.1e} MeV)")
    lhs = 6.0 * k2(M['d'], M['t'], f[:, 0], rel) * f[:, 3]
    rhs = 2.0 * k2(M['n'], M['a'], b[:, 0], rel) * b[:, 3]
    r = np.max(np.abs(lhs - rhs) / np.abs(rhs))
    check(r < 1e-8, f"(c) {run}: detailed balance g k^2 sigma(dt->na) = g k^2 sigma(na->dt) to {r:.1e}")

# (d) relativistic.out
lines = open(os.path.join(work, "he5_rel", "output", "relativistic.out")).read().splitlines()
lev = [l for l in lines if l.startswith("J =")][0].split()
ex, es = float(lev[lev.index("E_x") + 2]), float(lev[lev.index("E_s") + 2])
m0, s0 = (M['n'] + M['a']) * U, -0.735
t0 = ex - s0
want = s0 + t0 + t0 * t0 / (2 * m0)
check(abs(es - want) < 2e-6, f"(d) relativistic.out: E_s {es:.6f} = S0 + T0 + T0^2/(2 m0) = {want:.6f} MeV")
fac = float(lev[-1].rstrip(")"))
check(abs(fac - (1 + t0 / m0)) < 1e-8, f"(d) relativistic.out: dE_s/dE_x {fac:.8f} = 1 + T0/m0")

# (e) relativistic lab -> c.m. angles and dOmega_lab/dOmega_cm against an explicit four-vector boost
def boost(e1, m1, m2, m3, m4, q, cth):
    """lab cos(theta) of particle 3 emitted at c.m. cos(theta*) = cth (Q = q from the separation energies)"""
    a, b, c, d = m1 * U, m2 * U, m3 * U, m4 * U
    ti = 2 * b * e1 / (np.sqrt((a + b) ** 2 + 2 * b * e1) + a + b)
    tf = ti + q
    pf = np.sqrt(k2(m3, m4, tf, True)) * HBARC
    e3 = np.sqrt(c * c + pf * pf)
    etot = e1 + a + b
    sqs = a + b + ti
    beta, gam = np.sqrt(e1 * (e1 + 2 * a)) / etot, etot / sqs
    sth = np.sqrt(max(0.0, 1 - cth * cth))
    pz = gam * (pf * cth + beta * e3)
    return pz / np.hypot(pf * sth, pz)


def lab_to_cm(e1, m1, m2, m3, m4, q, lab):
    target = np.cos(np.radians(lab))
    lo, hi = -1.0, 1.0                      # cos(theta_lab) increases with cos(theta*) (g < 1 here)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if boost(e1, m1, m2, m3, m4, q, mid) < target: lo = mid
        else: hi = mid
    x = 0.5 * (lo + hi)
    h = 1e-6
    xa, xb = max(-1.0, x - h), min(1.0, x + h)
    jac = (boost(e1, m1, m2, m3, m4, q, xb) - boost(e1, m1, m2, m3, m4, q, xa)) / (xb - xa)
    return np.degrees(np.arccos(x)), jac       # dOmega_lab/dOmega_cm = d cos(lab)/d cos(cm)


for run, fname, datf, masses, q in (
        ("he5_ang_rel", "AZUREOut_aa=2_R=1.out", "dt_ang.dat", (M['d'], M['t'], M['n'], M['a']), 16.854 + 0.735),
        ("pa_ang_rel", "AZUREOut_aa=1_R=1.out", "pa_ang.dat", (M['p'], M['a'], M['p'], M['a']), 0.0)):
    d = load(run, fname)
    lab = np.loadtxt(os.path.join(work, run, "data", datf), ndmin=2)
    wa = wf = 0.0
    for (e1, al, _, _), row in zip(lab, d):
        ang, jac = lab_to_cm(e1, *masses, q, al)
        wa = max(wa, abs(row[2] - ang))
        wf = max(wf, abs(row[5] - jac) / jac)
    check(wa < 2e-5 and wf < 2e-6, f"(e) {run}: c.m. angles (max diff {wa:.1e} deg) and dOmega_lab/dOmega_cm "
                                   f"(rel {wf:.1e}) = an explicit four-vector boost")
    nr = load(run.replace("_rel", ""), fname)
    print(f"        {run}: relativistic - NR: c.m. angle up to {np.max(np.abs(d[:, 2] - nr[:, 2])):.4f} deg, "
          f"c.m. cross section up to {np.max(np.abs(d[:, 5] / nr[:, 5] - 1)):.3%}")

print(f"{fails} failure(s)")
sys.exit(min(fails, 100))
