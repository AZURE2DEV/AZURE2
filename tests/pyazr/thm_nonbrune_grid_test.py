#!/usr/bin/env python3
"""The sub-point grid of a folded THM point follows a narrow level moved by a
fit also when the formal parameters are fitted (Brune formalism off).

EPoint::RefreshSubPointGrid rebuilds a point's grid when the current fit
parameters anchor a narrow level differently from the ones it was built with
(tests/thm_narrow_fold, part B, checks it under Brune, from the CLI).  It did
so under the Brune formalism only: with formal parameters -- pyazr's
``use_brune=False``, the only way to turn Brune off -- a fit moving the 33 eV
5- level of 17O(n,a) left it between coarse sub-points of the grid built at the
input energy.  The anchors are now placed at the Thomas estimate of the
observed energy, E_lambda - sum gamma^2 (S - B)/(1 + sum gamma^2 dS/dE).

Model: tests/17O (n+17O THM, 21 keV lab Gaussian) on 38 points E_cm 20-131
keV, the narrow level (E_cm 74 keV) with a free energy, as in
tests/thm_narrow_fold.  With Brune off, a session built with the level 1.5 keV
higher gives formal parameters p_moved; a session built at the original energy
is evaluated at p_moved.  Its folded curve must match a fold, done here, of
its own unfolded curve on a dense grid (0.5 eV within 2 keV of the level, a
third session at p_moved) to 1e-3 pointwise -- the tolerance of
thm_narrow_fold.  (The two builds are not compared with each other: with
formal parameters the channel constants B_c = S_c(E_1) follow the input
energies, so the same p is a slightly different model in each.)  Before, the
grid stayed at the input energy and the fold was off by up to 7 %.

Run from anywhere:  python3 tests/pyazr/thm_nonbrune_grid_test.py
"""
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "17O")
LEVEL_EX = 1.891
MOVED_EX = LEVEL_EX + 0.0015
TARGETINT = '1 "1" 150 1 0.021 0 0 "" 0 0 0 0 "" 0'

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def project(dst, level_ex):
    """tests/17O with the evaluation points, the level at level_ex (free
    energy) and no <parameterSettings>, as thm_narrow_fold builds it."""
    for d in ("data", "output", "checks"):
        os.makedirs(os.path.join(dst, d))
    with open(os.path.join(dst, "data", "pts.dat"), "w") as f:
        e = 0.020
        while e < 0.1315:
            f.write(f"{e * 18 / 17:.9f} 0 1 0.1\n")
            e += 0.003
    out, skip = [], False
    for line in open(os.path.join(SOURCE, "17O.azr")).read().split("\n"):
        tag = line.strip()
        if tag == "<segmentsData>":
            out += [line, "1 1 2 -1 1 0 180 10 1.0 0 0 0 0 0 data/pts.dat 0 0"]
            skip = True
            continue
        if tag == "<targetInt>":
            out += [line, TARGETINT]
            skip = True
            continue
        if tag == "<parameterSettings>":
            out.append(line)
            skip = True
            continue
        if tag in ("</segmentsData>", "</targetInt>", "</parameterSettings>"):
            skip = False
        if skip:
            continue
        f = line.split()
        if len(f) > 30 and abs(float(f[2]) - LEVEL_EX) < 1e-9:
            f[2] = f"{level_ex:.6f}"
            f[3] = "0"
            line = " ".join(f)
        out.append(line)
    with open(os.path.join(dst, "run.azr"), "w") as f:
        f.write("\n".join(out))
    return os.path.join(dst, "run.azr")


SIGMA_CM = 0.021 * 17 / 18
TOL = 1e-3


def dense_points():
    """E_cm 0.5 eV within 2 keV of the level, 20 eV elsewhere, +-5.2 sigma
    beyond the evaluation range."""
    lo, hi, er = 0.020 - 5.2 * SIGMA_CM, 0.131 + 5.2 * SIGMA_CM, 0.074
    far = np.arange(lo, hi, 2.0e-5)
    far = far[(far < er - 0.002) | (far > er + 0.002)]
    return np.sort(np.concatenate([far, np.arange(er - 0.002, er + 0.002, 5.0e-7)]))


with tempfile.TemporaryDirectory() as tmp:
    moved = project(os.path.join(tmp, "moved"), MOVED_EX)
    orig = project(os.path.join(tmp, "orig"), LEVEL_EX)
    dense = project(os.path.join(tmp, "dense"), LEVEL_EX)
    # the dense copy: no folding, the dense grid as its data (lab energies)
    with open(dense) as f:
        text = f.read().replace(TARGETINT + "\n", "")
    with open(dense, "w") as f:
        f.write(text)
    with open(os.path.join(tmp, "dense", "data", "pts.dat"), "w") as f:
        for e in dense_points():
            f.write(f"{e * 18 / 17:.10f} 0 1 0.1\n")

    with azure2(moved, use_brune=False) as m:
        p_moved = np.asarray(m.params_rwa, float)
    with azure2(orig, use_brune=False) as o:
        p_orig = np.asarray(o.params_rwa, float)
        y = np.asarray(o.calculate_rwa(p_moved)[0], float)
        e_cm = np.asarray(o.energies[0], float)
    with azure2(dense, use_brune=False) as d:
        yd = np.asarray(d.calculate_rwa(p_moved)[0], float)
        ed = np.asarray(d.energies[0], float)

    moved_by = np.max(np.abs(p_moved - p_orig))
    print(f"  formal parameters of the two builds differ by up to {moved_by:.3g} MeV")
    check("the level moved in the formal parameters", moved_by > 1e-3, moved_by)
    ok = y.size == e_cm.size >= 30 and np.all(y > 0) and yd.size == ed.size > 20000
    check("every point evaluated", ok, (y.size, e_cm.size, yd.size))
    if ok:
        order = np.argsort(ed)
        ed, yd = ed[order], yd[order]
        worst, at = 0.0, 0.0
        for e0, y0 in zip(e_cm, y):
            k = np.exp(-0.5 * ((ed - e0) / SIGMA_CM) ** 2)
            ref = np.sum(0.5 * (k[1:] * yd[1:] + k[:-1] * yd[:-1]) * np.diff(ed)) / \
                np.sum(0.5 * (k[1:] + k[:-1]) * np.diff(ed))
            r = abs(y0 / ref - 1)
            if r > worst:
                worst, at = r, e0
        print(f"  fold at p_moved vs the dense reference: max |rel err| {worst:.2e} "
              f"(E_cm {1000 * at:.1f} keV)")
        check(f"the grid follows the level moved by the fit (rel {TOL:g})", worst < TOL, worst)

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: the THM sub-point grid follows a narrow level under formal parameters")
