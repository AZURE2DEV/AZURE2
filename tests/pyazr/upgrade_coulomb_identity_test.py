#!/usr/bin/env python3
"""Identity check (ARCHITECTURE.md section 10): the engine's channel functions against mpmath, and the Wronskian.

For a channel of orbital angular momentum l, Sommerfeld parameter eta = Z1 Z2 alpha sqrt(mu c^2 / (2 E_cm)) and
rho = k a (k = sqrt(2 mu c^2 E_cm) / hbar c, a the channel radius in fm), the regular and irregular Coulomb
functions F_l, G_l and their rho-derivatives must

  (a) equal mpmath.coulombf / coulombg (and their numerical derivatives at 30 digits) to 1e-8 relative;
  (b) satisfy the Wronskian F_l' G_l - F_l G_l' = 1 to 1e-10;
  (c) give the penetrability P = rho / (F^2 + G^2), the shift function S = rho (F F' + G G') / (F^2 + G^2) and the
      hard-sphere phase delta_hs = -atan(F / G) that the engine reports, to 1e-8 relative (P; S relative to
      max(|S|, 1e-3)) and 1e-9 rad.  The engine forms P, S and delta_hs from its own F, G, whose 5e-10 error
      propagates to 2.4e-9 in S at the worst point.

One exception is stated, not hidden: for a neutral channel (eta = 0) the engine's regular function F_l at small rho
and l >= 3 comes from the closed-form Riccati-Bessel expression, which loses digits to cancellation as
rho^-(2l+1): F_4(rho = 0.184) is 7.9e-6 relative from mpmath (n + 7Be at E_cm = 50 keV), F_3 4e-9.  P, S and
delta_hs are unaffected (F << G there).  Those points are checked separately at 1e-5, and the gap is a target for
the pycoulomb channel functions of Phase 1 (ROADMAP.md).

Channels: p + 12C of T1 (eta up to 2.6) and p + 7Li and n + 7Be of T3 (eta = 0: Riccati-Bessel functions),
l = 0..4, E_cm = 0.05, 0.3, 1, 3 MeV, at the channel radius.  The engine's constants are those of
include/Constants.h (hbar c = 197.3269804 MeV fm, 931.49410242 MeV/u, alpha = 1/137.035999084); masses are the
pair's atomic masses from the .azr.  An error in eta, in the mass unit or in the derivative convention (d/drho, not
d/dr) fails (a) at the 1e-3 level or worse.  Engine module via pyazr; skips without it.

    python3 tests/pyazr/upgrade_coulomb_identity_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
UPGRADE = os.path.join(ROOT, "tests", "upgrade")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    import mpmath as mp
    from pyazr import Session
    sys.path.insert(0, os.path.join(UPGRADE, "lib"))
    import upgrade_project as up
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)   # ctest SKIP_RETURN_CODE (tests/pyazr/CMakeLists.txt)

HBARC, UCONV, ALPHA = 197.3269804, 931.49410242, 1.0 / 137.035999084
mp.mp.dps = 30


def mp_functions(l, eta, rho):
    f = lambda r: mp.coulombf(l, eta, r)  # noqa: E731
    g = lambda r: mp.coulombg(l, eta, r)  # noqa: E731
    return (float(f(rho)), float(mp.diff(f, rho)), float(g(rho)), float(mp.diff(g, rho)))


def channel(project, pair_key, m1, m2, z1, z2, a):
    mu = m1 * m2 / (m1 + m2) * UCONV
    energies = np.array([0.05, 0.3, 1.0, 3.0])
    worst = {"F": 0.0, "Fsmall": 0.0, "W": 0.0, "PS": 0.0, "d": 0.0}
    work = tempfile.mkdtemp(prefix=f"upgrade_coulomb_{project}.")
    try:
        _channel(work, project, pair_key, m1, m2, z1, z2, a, mu, energies, worst)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    _report(project, pair_key, z1 * z2, worst)


def _channel(work, project, pair_key, m1, m2, z1, z2, a, mu, energies, worst):
    with Session(up.project_copy(project, work), data_mode=False) as m:
        for l in range(5):
            cf = m.coulomb_functions(pair_key, energies, L=l)
            for i, e in enumerate(energies):
                k = np.sqrt(2.0 * mu * e) / HBARC
                rho = k * a
                eta = z1 * z2 * ALPHA * np.sqrt(mu / (2.0 * e))
                F, dF, G, dG = mp_functions(l, eta, rho)
                got = (cf["F"][i], cf["dF"][i], cf["G"][i], cf["dG"][i])
                for n_, (x, y) in enumerate(zip(got, (F, dF, G, dG))):
                    key = "Fsmall" if (z1 * z2 == 0 and l >= 3 and n_ < 2) else "F"
                    worst[key] = max(worst[key], abs(x - y) / abs(y))
                worst["W"] = max(worst["W"], abs(got[1] * got[2] - got[0] * got[3] - 1.0))
                P = rho / (F * F + G * G)
                S = rho * (F * dF + G * dG) / (F * F + G * G)
                worst["PS"] = max(worst["PS"], abs(cf["P"][i] - P) / P, abs(cf["S"][i] - S) / max(abs(S), 1e-3))
                d = (cf["delta_hs"][i] + np.arctan2(F, G) + np.pi) % (2 * np.pi) - np.pi
                worst["d"] = max(worst["d"], abs(d))


def _report(project, pair_key, z1z2, worst):
    label = f"{project} pair {pair_key}"
    check(f"{label}: F, F', G, G' = mpmath (l = 0..4, 4 energies)", worst["F"] < 1e-8, f"largest {worst['F']:.1e}")
    check(f"{label}: Wronskian F'G - FG' = 1", worst["W"] < 1e-10, f"largest |W - 1| {worst['W']:.1e}")
    if z1z2 == 0:
        check(f"{label}: F, F' for l = 3, 4 (closed-form cancellation, see the docstring)", worst["Fsmall"] < 1e-5,
              f"largest {worst['Fsmall']:.1e}")
    check(f"{label}: P and S from F, G", worst["PS"] < 1e-8, f"largest {worst['PS']:.1e}")
    check(f"{label}: hard-sphere phase -atan(F/G)", worst["d"] < 1e-9, f"largest {worst['d']:.1e} rad")


channel("T1", 1, 1.00782503223, 12.0, 1, 6, 3.6)
channel("T3", 2, 1.00782503223, 7.0160034366, 1, 3, 4.0)
channel("T3", 3, 1.00866491595, 7.016928717, 0, 4, 4.0)

print("FAIL" if failures else "ok")
sys.exit(1 if failures else 0)
