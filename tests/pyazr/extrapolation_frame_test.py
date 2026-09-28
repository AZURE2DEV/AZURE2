#!/usr/bin/env python3
"""AzrModel.add_extrapolation says which frame its energy grid is in, and the
engine evaluates the grid asked for.

A <segmentsTest> line holds lab energies (the light particle of the entrance
pair on the heavy one at rest), as data files do; add_extrapolation's
docstring called its grid entrance-channel c.m. MeV, so a c.m. grid came out
stretched by (M1 + M2)/M2.  It now takes ``frame="lab"`` (the default, what
the engine reads) or ``frame="cm"`` (converted to lab with the pair's masses).

On tests/17O (entrance pair 1, n + 17O), in extrapolation mode:

  1. frame="cm", 0.1-0.5 MeV in 0.1 steps: calculate_energies reports those
     c.m. energies (1e-9);
  2. frame="lab" and the default, same numbers: c.m. = lab M2/(M1 + M2);
  3. an unknown frame raises.

Run from anywhere:  python3 tests/pyazr/extrapolation_frame_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "17O")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2, AzrModel
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)

grid = np.array([0.1, 0.2, 0.3, 0.4, 0.5])

with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "17O")
    shutil.copytree(os.path.join(SOURCE, "data"), os.path.join(work, "data"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(work, d))
    mdl = AzrModel.from_file(os.path.join(SOURCE, "17O.azr"))
    pair = mdl._pair_template(1)
    m1, m2 = pair.M1, pair.M2
    mdl.clear_extrapolations()
    mdl.add_extrapolation(1, 2, 0.1, 0.5, 0.1, frame="cm")
    mdl.add_extrapolation(1, 2, 0.1, 0.5, 0.1, frame="lab")
    mdl.add_extrapolation(1, 2, 0.1, 0.5, 0.1)
    path = os.path.join(work, "extrap.azr")
    mdl.write(path)
    try:
        mdl.add_extrapolation(1, 2, 0.1, 0.5, 0.1, frame="CM")
        check("3. an unknown frame raises", False, "accepted 'CM'")
    except ValueError:
        check("3. an unknown frame raises", True)

    with azure2(path, data_mode=False) as m:
        e = [np.asarray(v, float) for v in m.calculate_energies(np.asarray(m.params, float))]
    print(f"  pair 1: M1 = {m1}, M2 = {m2}; c.m. energies per segment:")
    for k, v in enumerate(e, 1):
        print(f"    segment {k}: {np.array2string(v, precision=6)}")
    ok = len(e) == 3 and all(v.size == grid.size for v in e)
    check("three extrapolation segments of five energies", ok, [v.size for v in e])
    if ok:
        check("1. frame='cm' is evaluated at the c.m. grid (1e-9)",
              np.max(np.abs(e[0] - grid)) < 1e-9, e[0] - grid)
        lab = grid * m2 / (m1 + m2)
        check("2. frame='lab' is converted like data (c.m. = lab M2/(M1+M2))",
              np.max(np.abs(e[1] - lab)) < 1e-9, e[1] - lab)
        check("2. the default is frame='lab'", np.array_equal(e[1], e[2]), e[2] - e[1])

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: add_extrapolation's frame= is what the engine evaluates")
