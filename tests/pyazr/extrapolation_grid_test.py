#!/usr/bin/env python3
"""pyazr.extrapolation_grid and pyazr.refine_grid: the engine's grids (docs/upgrade/project_format3.md section 3.8).

  1. extrapolation_grid gives the values written by hand of range step / points / points log / list / file, keV and
     MeV, and a c.m. grid with the entrance pair's masses (E_lab = E_cm (m1 + m2)/m2); a c.m. grid without masses, a
     grid that does not parse and an angle grid of kind adaptive raise ValueError with the engine's message;
  2. refine_grid on a Python Lorentzian (1 eV at 0.3001234 MeV, tolerance 1e-3) gives the grid of the C++ test
     (experiment_grids_test: fewer than 400 points, the peak and the integral to 1e-3), a vector-valued model the
     union grid, and max_points reached is reported;
  3. a pyazr session on a project with an adaptive grid (tests/extrapolation_grids/narrow.azr.in) in extrapolation mode
     is refused: the grid is made by the engine's extrapolation run at one parameter vector.
Needs the compiled engine; skips (exit 77) without it.
"""
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def near(a, b, rel=1e-14):
    return abs(a - b) <= rel * max(abs(a), abs(b))


try:
    sys.path.insert(0, ROOT)
    import pyazr
    from pyazr import _azure2  # noqa: F401
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

# ---- 1. static grids
g = pyazr.extrapolation_grid("range 0.1 0.4 MeV lab step 0.1")
check("range step: lo + k step", g == [0.1, 0.1 + 0.1, 0.1 + 2 * 0.1, 0.1 + 3 * 0.1], g)
g = pyazr.extrapolation_grid("range 100 300 keV lab points 3 log")
check("range points log in keV", len(g) == 3 and near(g[0], 0.1) and near(g[1], 0.17320508075688773) and near(g[2], 0.3), g)
g = pyazr.extrapolation_grid("list 1.2 0.6 MeV cm", masses=(1.00728, 12.0))
check("list in c.m. with masses: 1.300728, 0.650364 MeV lab", near(g[0], 1.300728) and near(g[1], 0.650364), g)
tmp = tempfile.mkdtemp()
try:
    with open(os.path.join(tmp, "grid.txt"), "w") as f:
        f.write("# keV\n100\n250  # comment\n")
    g = pyazr.extrapolation_grid('file "grid.txt" keV lab', project_dir=tmp)
    check("file grid relative to project_dir", len(g) == 2 and near(g[0], 0.1) and near(g[1], 0.25), g)
finally:
    shutil.rmtree(tmp)
g = pyazr.extrapolation_grid("list 0 90 180 deg cm", angles=True)
check("an angle grid", g == [0.0, 90.0, 180.0], g)
for text, kw, what in [("list 0.3 MeV cm", {}, "a c.m. grid without masses"),
                       ("range 0.1 MeV lab step 0.1", {}, "a grid that does not parse"),
                       ("adaptive 0 180 deg cm tolerance 1e-3 max_points 10 start_step 10", {"angles": True},
                        "an adaptive angle grid")]:
    try:
        pyazr.extrapolation_grid(text, **kw)
        check(f"{what}: ValueError", False, "no exception")
    except ValueError as e:
        check(f"{what}: ValueError ({str(e)[:60]}...)", len(str(e)) > 0)

# ---- 2. refinement
centre, width = 0.3001234, 1e-6


def lorentzian(x, c=centre, w=width):
    h = w / 2
    return h * h / ((x - c) ** 2 + h * h)


start = [0.2 + k * 0.01 for k in range(21)] + [centre + k * width for k in (0, -0.5, 0.5, -1, 1, -2, 2, -4, 4)]
r = pyazr.refine_grid(start, lambda xs: [lorentzian(x) for x in xs], 1e-3, 2000)
x, v = r["x"], [row[0] for row in r["values"]]
integral = sum(0.5 * (x[i] - x[i - 1]) * (v[i] + v[i - 1]) for i in range(1, len(x)))
h = width / 2
exact = h * (math.atan((0.4 - centre) / h) + math.atan((centre - 0.2) / h))
check(f"refine_grid: {len(x)} points < 400, peak and integral to 1e-3",
      len(x) < 400 and abs(max(v) - 1) < 1e-3 and abs(integral / exact - 1) < 1e-3 and not r["budget_reached"],
      (len(x), max(v), integral / exact - 1))
two = pyazr.refine_grid(start[:21], lambda xs: [(lorentzian(t, 0.25, 1e-4), lorentzian(t, 0.35, 1e-4)) for t in xs],
                        1e-3, 2000)
check("refine_grid, two components: refined at both lines",
      sum(abs(t - 0.25) < 1e-4 for t in two["x"]) > 10 and sum(abs(t - 0.35) < 1e-4 for t in two["x"]) > 10)
b = pyazr.refine_grid(start, lambda xs: [lorentzian(t) for t in xs], 1e-3, 100)
check("refine_grid, max_points 100: budget_reached, 100 points, the deviation left above the tolerance",
      b["budget_reached"] and len(b["x"]) == 100 and b["largest_deviation"] > 1e-3, b["largest_deviation"])

# ---- 3. a session refuses an adaptive grid
tmp = tempfile.mkdtemp()
try:
    with open(os.path.join(ROOT, "tests", "extrapolation_grids", "narrow.azr.in")) as f:
        text = f.read().replace("@GRID@", "adaptive 0.2 0.4 MeV cm tolerance 1e-3 max_points 2000 start_step 0.01")
    path = os.path.join(tmp, "narrow.azr")
    with open(path, "w") as f:
        f.write(text)
    os.makedirs(os.path.join(tmp, "output"))
    refused = False
    try:
        with pyazr.Session(path, cwd=tmp, extrapolations=True) as m:
            m.evaluate()
    except Exception as e:  # the session's set-up failure
        refused = True
        print(f"        ({type(e).__name__}: {str(e)[:100]})")
    check("a pyazr session in extrapolation mode refuses an adaptive grid", refused)
finally:
    shutil.rmtree(tmp)

print(f"{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
