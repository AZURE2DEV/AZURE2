#!/usr/bin/env python3
"""Identity check (ARCHITECTURE.md section 10): the optical theorem with open reaction channels (T3, n + 7Be).

The total cross section of the neutral entrance pair n + 7Be from the forward elastic amplitude,
sigma_tot = (2 pi / k^2) g sum_J (2J + 1) sum_c (1 - Re U_cc), must equal the sum of the angle-integrated cross
sections into every particle pair the R-matrix couples it to: sigma(n,n) + sigma(n,p) + sigma(n,alpha alpha).  This
holds only if U is unitary over all channels and each exit cross section is normalized consistently, including the
identical-boson exit alpha + alpha (one event per decay, symmetrized amplitude).  A factor 2 in the identical-particle
normalization breaks it by sigma(n,alpha)/sigma_tot, here 7e-4 to 1.5e-2 of the total.

Uses the <segmentsTest> lines of tests/upgrade/T3/T3.azr (menu 3): the n + 7Be total ("total 1"), elastic, (n,p) and
(n,alpha) for E_lab = 0.01-0.39 MeV, step 0.01 MeV (the engine's
step loop stops below the upper bound 0.40).  The beta pair is not a channel of U and photon pairs are absent.
Tolerance 1e-9 relative (outputs carry 11 significant digits).

    python3 tests/pyazr/upgrade_optical_theorem_test.py     (binary: $AZURE2_BIN, $AZURE2_BINARY or the newest build)
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

work = tempfile.mkdtemp(prefix="upgrade_optical.")
try:
    up.run(binary, work, up.project_copy("T3", work), up.EXTRAP)
    out = os.path.join(work, "output")
    total, elastic = up.read_blocks(os.path.join(out, "AZUREOut_aa=3_R=3.extrap"))
    np_ = up.read_blocks(os.path.join(out, "AZUREOut_aa=3_R=2.extrap"))[0]
    na = up.read_blocks(os.path.join(out, "AZUREOut_aa=3_R=1.extrap"))[0]
finally:
    shutil.rmtree(work, ignore_errors=True)

check("one energy grid of 39 points in each of the four segments", all(len(b) == 39 for b in (total, elastic, np_, na)),
      ", ".join(str(len(b)) for b in (total, elastic, np_, na)))
check("same c.m. energies", all(np.array_equal(total[:, 0], b[:, 0]) for b in (elastic, np_, na)))
parts = elastic[:, 3] + np_[:, 3] + na[:, 3]
worst = np.max(np.abs(total[:, 3] / parts - 1.0))
check("sigma_tot = sigma(n,n) + sigma(n,p) + sigma(n,alpha)", worst <= 1e-9, f"largest relative difference {worst:.1e}")
share = na[:, 3] / total[:, 3]
check("the identical-boson exit carries a share above the tolerance", share.min() > 1e-6,
      f"sigma(n,alpha)/sigma_tot {share.min():.1e}-{share.max():.1e}")
print("FAIL" if failures else "ok")
sys.exit(1 if failures else 0)
