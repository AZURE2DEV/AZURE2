#!/usr/bin/env python3
"""The first evaluation of a session uses the shift it is given, not the .azr's.

A segment's energy shift starts at the value on its line, but its points start unshifted.  Up to 2026-10-29
the "last applied" shift started at 0, so a first request of exactly 0 (constant shift or sqrt(E) coefficient)
was taken as already applied and skipped; the sub-point grid re-placement (EData::ReplaceSubPointGrids) then
re-energised the segment's points with the shift still set from the file.  The first call returned the chi2 of
the FILE's shift (17O+a Kunz profile, 2026-10-29: 0 keV gave the +8 keV chi2 to every digit).  It needs a
nonzero shift on the line and a <targetInt> sub-point effect on the segment.

Model: tests/13N with only segment 3 (Meyer, 12C(p,p) at 84.3 deg lab) active, a Gaussian beam spread of 10 keV on
it (so the 0.457 MeV resonance lies in the sub-point windows of the points near it and the grid re-placement moves
their grids with the shift; a capture segment would not do, its grids are never re-placed), and a shift of 10 keV (case a) or a sqrt(E) coefficient of 0.01 (case b) on its line, varied.  In one session:
x = 0, 0, file value, 0.  Every call at 0 must agree to 1e-12 relative, and differ from the file value's.

Run from anywhere:  python3 tests/pyazr/shift_first_call_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import numpy as np                         # noqa: E402
    from pyazr import azure2                   # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(0)

SRC = os.path.join(ROOT, "tests", "13N", "13N.azr")


def make_model(work, sqrt):
    lines = open(SRC).read().split("\n")
    i0, i1 = lines.index("<segmentsData>"), lines.index("</segmentsData>")
    for i in range(i0 + 1, i1):
        f = lines[i].split()
        if not f:
            continue
        if i == i0 + 3:                        # segment 3: shift on, varied
            if sqrt:
                f += ["sqrtshift", "0.01", "0", "1"]
            else:
                f[11], f[13] = "0.01", "1"
        else:
            f[0] = "0"
        lines[i] = " ".join(f)
    j0, j1 = lines.index("<targetInt>"), lines.index("</targetInt>")
    lines[j0 + 1:j1] = ['1 "3" 50 1 0.01 0 0 "" 0 0 0']
    path = os.path.join(work, "first_call.azr")
    open(path, "w").write("\n".join(lines))
    return path


def run(label, sqrt):
    print(f"\n{label}")
    work = tempfile.mkdtemp(prefix="pyazr_first_call_")
    try:
        shutil.copytree(os.path.join(ROOT, "tests", "13N", "data"), os.path.join(work, "data"))
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        a = azure2(make_model(work, sqrt), cwd=work)
        x0 = np.array(a.params_rwa, float)
        names = [p.name for p in a.parameters if not p.fixed]
        kind = "shift_sqrt" if sqrt else "shift"
        free = [p for p in a.parameters if not p.fixed]
        k = [i for i, p in enumerate(free) if p.kind == kind][0]
        check(f"the file value {x0[k]} is in the vector", abs(x0[k] - 0.01) < 1e-12, (names[k], x0[k]))
        chi = []
        for v in (0.0, 0.0, x0[k], 0.0):
            x = x0.copy()
            x[k] = v
            chi.append(a.calculate_chi2_rwa(x)[0])
        print("    chi2 at 0, 0, file, 0: " + ", ".join(f"{c:.6f}" for c in chi))
        zero = chi[3]
        check("first call at 0 = a later call at 0", abs(chi[0] - zero) <= 1e-12 * abs(zero), (chi[0], zero))
        check("second call at 0 = a later call at 0", abs(chi[1] - zero) <= 1e-12 * abs(zero), (chi[1], zero))
        check("the file value gives another chi2", abs(chi[2] - zero) > 1e-6 * abs(zero), (chi[2], zero))
    finally:
        shutil.rmtree(work, ignore_errors=True)


run("a. constant shift 10 keV on the line, first request 0", sqrt=False)
run("b. sqrt(E) coefficient 0.01 on the line, first request 0", sqrt=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
