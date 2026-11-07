#!/usr/bin/env python3
"""Identity check (ARCHITECTURE.md section 10): with elastic scattering the only open channel, U is unitary.

For a single open channel each partial wave has U_J = exp(2 i delta_J) with real delta_J, so the total cross
section from the optical theorem, sigma_tot = (2 pi / k^2) sum_J g_J (1 - Re U_J), equals the angle-integrated
elastic cross section sigma_el = (pi / k^2) sum_J g_J |1 - U_J|^2 at every energy.  A collision matrix that loses or
gains flux (a wrong sign of the hard-sphere phase in one place, a level matrix that is not symmetric, a width
counted twice) breaks the equality at the size of the error.

Model: n + 16O (neutral, so the optical theorem applies to the elastic segment with "total 1"), built here with
tests/upgrade/lib: a 1/2+ s-wave level at E_x = 5.0 MeV (Gamma_n = 100 keV), a 3/2- p-wave level at 5.5 MeV
(50 keV), a 5/2+ d-wave level at 6.0 MeV (20 keV) and a 1/2+ background pole at 15 MeV (5 MeV); a = 4.0 fm, Brune
basis.  Test segments (menu 3) from E_lab = 0.05 to 3.0 MeV, step 0.01 MeV (296 energies), elastic and total.
Tolerance |sigma_tot / sigma_el - 1| <= 1e-9 (the outputs carry 11 significant digits).

    python3 tests/pyazr/upgrade_unitarity_test.py     (binary: $AZURE2_BIN, $AZURE2_BINARY or the newest build)
"""
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

M16O, M17O = 15.99491461957, 16.99913175650
n16o = up.Pair(1, 0.5, 1, 0.0, 1, up.MASS["n"], M16O, 0, 8, up.separation(up.MASS["n"], M16O, M17O), R=4.0)
C = up.Channel
model = up.Project(
    "unitarity", [n16o],
    [up.Level(0.5, 1, 5.0, [C(n16o, 1, 0, 1.0e5)]),
     up.Level(1.5, -1, 5.5, [C(n16o, 1, 2, 5.0e4)]),
     up.Level(2.5, 1, 6.0, [C(n16o, 1, 4, 2.0e4)]),
     up.Level(0.5, 1, 15.0, [C(n16o, 1, 0, 5.0e6)])],
    data=[],
    tests=[up.TestSegment(1, 1, 0.05, 3.0, 0.01), up.TestSegment(1, 1, 0.05, 3.0, 0.01, extra="total 1")])

work = tempfile.mkdtemp(prefix="upgrade_unitarity.")
try:
    up.run(binary, work, model.write(work), up.EXTRAP)
    el, tot = up.read_blocks(os.path.join(work, "output", "AZUREOut_aa=1_R=1.extrap"))
finally:
    shutil.rmtree(work, ignore_errors=True)

check("296 energies in each segment", len(el) == len(tot) == 296, f"{len(el)}, {len(tot)}")
check("same c.m. energies", np.array_equal(el[:, 0], tot[:, 0]))
ratio = tot[:, 3] / el[:, 3]
worst = np.max(np.abs(ratio - 1.0))
check("sigma_tot (optical theorem) = sigma_el, single open channel", worst <= 1e-9,
      f"largest |ratio - 1| {worst:.1e}; sigma_el {el[:, 3].min():.3g}-{el[:, 3].max():.3g} b")
print("FAIL" if failures else "ok")
sys.exit(1 if failures else 0)
