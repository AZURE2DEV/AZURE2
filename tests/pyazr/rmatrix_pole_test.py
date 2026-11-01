#!/usr/bin/env python3
"""The R-matrix route at a pole of the R-matrix: finite, and equal to the A-matrix route.

The R-matrix has a pole at the energy of every level with a width: the formal energy in the standard formalism, the
observed one in Brune's (there the level's row of the real level matrix vanishes).  The adaptive sub-point grid of a
target effect places a node exactly at a level's energy, so a thick-target capture segment evaluates the R-matrix
there.  Before 2026-11-01 the R-matrix route (<config> first line false) then gave NaN in the standard formalism
(R = gamma^2 / 0) and an ill-conditioned inverse with Brune: tests/13N's artemov capture (target integration on)
was NaN at every point above the 1/2+ resonance, and 1.2 % off at the first one with Brune.  The group is now solved
in the level form at a pole, gamma^T (M - gamma L gamma^T)^{-1} gamma, which is finite there.

Model: tests/13N with only the artemov capture segment active, target integration as in the file; A-matrix and
R-matrix route, Brune and standard, ignore_externals=False (the R-matrix route does not implement that option).
Every point must be finite and equal to the A-matrix route's (relative 1e-10).  A second pass with the target
integration switched off is the control (no sub-points: equal before the fix too).

Run from anywhere:  python3 tests/pyazr/rmatrix_pole_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
SRC = os.path.join(ROOT, "tests", "13N")


def write_model(path, amatrix, target_int):
    lines = open(os.path.join(SRC, "13N.azr")).read().split("\n")
    i0, i1 = lines.index("<segmentsData>"), lines.index("</segmentsData>")
    for k in range(i0 + 1, i1):
        if lines[k].split() and "artemov" not in lines[k]:
            lines[k] = "0" + lines[k].lstrip()[1:]
    if not target_int:
        t0, t1 = lines.index("<targetInt>"), lines.index("</targetInt>")
        for k in range(t0 + 1, t1):
            if lines[k].split():
                lines[k] = "0" + lines[k].lstrip()[1:]
    c = lines.index("<config>") + 1
    lines[c] = ("true" if amatrix else "false") + "  " + lines[c].split(None, 1)[1]
    open(path, "w").write("\n".join(lines))


def main():
    try:
        import numpy as np
        from pyazr import azure2
    except Exception as err:  # same skip as the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 0
    failures = []

    def check(name, ok, detail=""):
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
        if not ok:
            failures.append(name)

    work = tempfile.mkdtemp(prefix="pyazr_rmatrix_pole_")
    try:
        shutil.copytree(os.path.join(SRC, "data"), os.path.join(work, "data"))
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))

        def calc(amatrix, target_int, brune):
            path = os.path.join(work, f"m_{int(amatrix)}{int(target_int)}{int(brune)}.azr")
            write_model(path, amatrix, target_int)
            with azure2(path, cwd=work, use_brune=brune, ignore_externals=False) as a:
                return np.asarray(a.calculate_rwa(np.array(a.params_rwa, float))[0], float)

        for target_int in (True, False):
            for brune in (True, False):
                label = f"{'Brune' if brune else 'standard'}, target integration {'on' if target_int else 'off'}"
                va, vr = calc(True, target_int, brune), calc(False, target_int, brune)
                bad = int((~np.isfinite(vr)).sum())
                check(f"{label}: every R-matrix-route point finite", bad == 0, f"{bad} of {len(vr)} not finite")
                d = np.abs(vr - va) / np.abs(va)
                d[~np.isfinite(vr)] = np.inf
                k = int(np.argmax(d))
                check(f"{label}: R-matrix route = A-matrix route (worst {d.max():.1e})", d.max() <= 1e-10,
                      f"point {k}: A {va[k]:.6e}, R {vr[k]:.6e}")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
