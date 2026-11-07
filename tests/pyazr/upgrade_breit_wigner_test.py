#!/usr/bin/env python3
"""Identity check (ARCHITECTURE.md section 10): an isolated narrow level gives the Breit-Wigner area.

For one isolated level of spin J in the entrance channel p + 12C (spins j1 = 1/2, j2 = 0) with observed partial widths
Gamma_p and Gamma_gamma at the resonance energy E_r (c.m.), the capture cross section near E_r is

    sigma(E) = (pi / k^2) g Gamma_p Gamma_gamma / ((E - E_r)^2 + Gamma^2 / 4),   g = (2J + 1) / ((2 j1 + 1)(2 j2 + 1)),

k the c.m. wave number at E_r, and Gamma the total width in the level denominator.  AZURE2 keeps photon channels
out of the level matrix (docs/source/user_guide/segments.rst, "Neutron Total Cross Sections"): capture is computed
to first order in the photon amplitudes, so Gamma = Gamma_p, not Gamma_p + Gamma_gamma.  The approximation is valid
while Gamma_gamma << Gamma_p; here it changes the peak by (1.1)^2 and the area by 1.1, and this test pins it (an
engine that puts Gamma_gamma into the denominator fails by 10 %).  Over the window |E - E_r| <= W its integral is

    A_W = (2 pi^2 / k^2) g Gamma_gamma (2 / pi) atan(2 W / Gamma)        (Gamma = Gamma_p).

The Brune parameters of the .azr are these observed widths at E_r, so the R-matrix result must give A_W up to the
energy dependence of the penetrability and of the shift over the window (relative size W dP/dE / P ~ 1e-5 here).
A wrong factor of g, pi/k^2, or a width taken at the wrong energy fails at the 1e-2 level or worse.

Model (built here): p + 12C -> 13N, one 1/2+ level at E_x = 2.3650 MeV (E_r = E_x - S_p = 0.42150 MeV), s wave,
Gamma_p = 1.0 eV, E1 photon width Gamma_gamma = 0.1 eV to the 1/2- ground state, no external capture (mask 0),
a = 3.6 fm.  The capture cross section is computed (menu 3) on a lab grid of step 6e-8 MeV (Gamma/17 in the c.m.),
W = 100 eV, and integrated by the trapezoidal rule (discretization error (h/Gamma)^2/12 ~ 3e-4 of the peak region,
< 1e-5 of the area).  Tolerance 1e-5 relative on the area (observed 1.1e-7 with the c1cc92c7 build) and 2e-4 on
the largest sampled value against the peak (the grid misses the peak by up to h/2; observed 7.6e-6).

    python3 tests/pyazr/upgrade_breit_wigner_test.py     (binary: $AZURE2_BIN, $AZURE2_BINARY or the newest build)
"""
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests", "upgrade", "lib"))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


try:
    import numpy as np
    import upgrade_project as up
    binary = up.find_binary()
except (ImportError, SystemExit) as e:
    print(f"skip: {e}")
    sys.exit(77)   # ctest SKIP_RETURN_CODE (tests/pyazr/CMakeLists.txt)

HBARC, UCONV = 197.3269804, 931.49410242
M = up.MASS
S_P = up.separation(M["p"], M["12C"], M["13N"])
EX, GP, GG = 2.3650, 1.0, 0.1                 # MeV, eV, eV
ER = EX - S_P                                  # c.m., MeV
W = 100e-6                                     # MeV
mu = M["p"] * M["12C"] / (M["p"] + M["12C"]) * UCONV
lab = (M["p"] + M["12C"]) / M["12C"]           # E_lab / E_cm (non-relativistic)

p12c = up.Pair(1, 0.5, 1, 0.0, 1, M["p"], M["12C"], 1, 6, S_P, R=3.6)
g0 = up.Pair(2, 1.0, 1, 0.5, -1, 0.0, M["13N"], 0, 7, 0.0, ptype=10)
C = up.Channel
step = 6e-8
model = up.Project("bw", [p12c, g0], [up.Level(0.5, 1, EX, [C(p12c, 1, 0, GP), C(g0, 1, 2, GG)])], data=[],
                   tests=[up.TestSegment(1, 2, round((ER - 1.05 * W) * lab, 10), round((ER + 1.05 * W) * lab, 10), step)])

work = tempfile.mkdtemp(prefix="upgrade_bw.")
try:
    up.run(binary, work, model.write(work), up.EXTRAP)
    rows = up.read_blocks(os.path.join(work, "output", "AZUREOut_aa=1_R=2.extrap"))[0]
finally:
    shutil.rmtree(work, ignore_errors=True)

e, sigma = rows[:, 0], rows[:, 3]              # E_cm (MeV), sigma (b)
inside = np.abs(e - ER) <= W
e, sigma = e[inside], sigma[inside]
area = float(np.sum(0.5 * (sigma[1:] + sigma[:-1]) * np.diff(e)))          # b MeV
lo, hi = e[0] - ER, e[-1] - ER                                             # the window actually sampled
k2 = 2.0 * mu * ER / HBARC ** 2                                            # fm^-2
G = GP * 1e-6                                                              # photon width not in the denominator
expected = (math.pi / k2) / 100.0 * 1.0 * (GP * GG * 1e-12) * (2.0 / G) * (
    math.atan(2 * hi / G) - math.atan(2 * lo / G))                          # b MeV, g = 1
check(f"{inside.sum()} points within |E - E_r| <= {W * 1e6:.0f} eV", inside.sum() > 3000)
rel = area / expected - 1.0
check("integral of sigma(p,gamma) = Breit-Wigner area", abs(rel) <= 1e-5,
      f"{area:.6e} b MeV against {expected:.6e}, relative {rel:+.1e}")
peak = sigma.max()
want_peak = (math.pi / k2) / 100.0 * GP * GG / (GP / 2) ** 2
check("peak = (pi/k^2) g 4 Gamma_gamma / Gamma_p", abs(peak / want_peak - 1.0) <= 2e-4,
      f"{peak:.6e} b against {want_peak:.6e}")
print("FAIL" if failures else "ok")
sys.exit(1 if failures else 0)
