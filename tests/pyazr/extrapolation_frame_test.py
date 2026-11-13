#!/usr/bin/env python3
"""An extrapolation record's energy grid states its frame, and the engine evaluates the grid asked for.

The engine reads lab energies (the light particle of the entrance pair on the heavy one at rest), as data files
hold; a grid in entrance-channel c.m. MeV was once taken as lab and came out stretched by (M1 + M2)/M2.  A format-3
grid names its frame (`energies range lo hi MeV cm|lab step s`); `cm` is taken to lab with the pair's masses
(project_format3.md section 3, experiment/Grids.h).

On tests/17O (entrance pair 1, n + 17O; its extrapolation removed), in extrapolation mode:

  1. `cm`, 0.1-0.5 MeV in 0.1 steps: evaluate().data[i].energies reports those c.m. energies (1e-9);
  2. `lab`, same numbers: c.m. = lab M2/(M1 + M2);
  3. a grid without a frame, or with an unknown one (`CM`), is refused (FormatError, a ValueError).

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
    from pyazr import Session, Project
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)

grid = np.array([0.1, 0.2, 0.3, 0.4, 0.5])

with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "17O")
    shutil.copytree(os.path.join(SOURCE, "data"), os.path.join(work, "data"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(work, d))
    mdl = Project.open(os.path.join(SOURCE, "17O.azr"))
    pair = [p for p in mdl.pairs() if p.key == 1][0]
    m1, m2 = pair.light_mass, pair.heavy_mass
    for name, r in mdl.records().items():
        if r.record_kind == "extrapolation":
            mdl.remove(name)
    reaction = mdl.reaction(1, 2)

    def grid_record(name, frame):
        mdl.add("extrapolation", name, reaction=reaction, quantity="sigma", unit="b",
                energies=f"range 0.1 0.5 MeV {frame} step 0.1".replace("  ", " "))

    grid_record("in_cm", "cm")
    grid_record("in_lab", "lab")
    path = os.path.join(work, "extrap.azr")
    mdl.save(path)
    for frame in ("", "CM"):
        try:
            grid_record("refused", frame)
            check(f"3. a grid with frame {frame!r} is refused", False, "accepted")
        except ValueError:
            check(f"3. a grid with frame {frame!r} is refused", True)

    with Session(path, data_mode=False) as m:
        e = [np.asarray(v, float) for v in [d.energies for d in m.evaluate().data]]
    print(f"  pair 1: M1 = {m1}, M2 = {m2}; c.m. energies per segment:")
    for k, v in enumerate(e, 1):
        print(f"    segment {k}: {np.array2string(v, precision=6)}")
    ok = len(e) == 2 and all(v.size == grid.size for v in e)
    check("two extrapolations of five energies", ok, [v.size for v in e])
    if ok:
        check("1. a cm grid is evaluated at the c.m. energies (1e-9)",
              np.max(np.abs(e[0] - grid)) < 1e-9, e[0] - grid)
        lab = grid * m2 / (m1 + m2)
        check("2. a lab grid is converted like data (c.m. = lab M2/(M1+M2))",
              np.max(np.abs(e[1] - lab)) < 1e-9, e[1] - lab)

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: the frame of an extrapolation grid is what the engine evaluates")
