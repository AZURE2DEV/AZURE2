#!/usr/bin/env python3
"""The model as a curve over a data set (``Session.model_curve``, ``pyazr.model_curve``).

1. At the data set's own points the curve is ``Session.evaluate``'s model: on T1 (elastic scattering at lab angles,
   analyzing powers at c.m. angles, a phase shift, capture integrated and at two lab angles) and on an alpha + 12C
   ratio to Rutherford (seven angular distributions) to rounding (1e-12); on T2 (a thin, a thick, an implanted and a
   gas target, a beam's spread with a detector's resolution, a sum, a ratio, a natural target, THM) an unfolded data
   set to rounding, a folded one without an energy shift to 1e-8 (the curves' folds on one mesh, as the data
   session's), a folded one with a shift and THM to 2e-4 (on grids of their own in the data session: the quadrature).
2. Between the points: T2's unfolded THM OES data set, T1's capture to the ground state (integrated and at a lab
   angle) and elastic scattering at a lab angle given data files at the curve's grid energies (a second project, the
   data path of the engine): its model equals the curve (1e-12).
3. An energy shift (a fixed one on an unfolded data set, a fitted one moved in x): the curve at the shifted lab
   energies equals the model of the shifted points.
4. The grid: ``points`` evenly spaced values per branch with the data set's own, sorted, within the first and last
   point; ``at_data`` marks the data set's points; branches per fixed angle or energy.
5. Refused, with the three parts: a data set whose columns give each point a window (T2's beam profile), an unknown
   data set (ParameterKeyError), a session of the extrapolations, points < 2; a record the engine refuses on a grid
   (an analyzing power at a lab angle) refuses its data set with the engine's words and the others are made.
6. ``Session.model_curve`` (one data set, its own session) equals the curve of ``model_curves`` at unfolded points,
   and ``pyazr.endf.Material.curves`` takes a curve as it takes a data set (a request at the curve's points).

Run from anywhere:  python3 tests/pyazr/model_curve_test.py
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
    import warnings

    import numpy as np                                  # noqa: E402
    import pyazr                                        # noqa: E402
    from pyazr.errors import ParameterError, ParameterKeyError  # noqa: E402
    from pyazr.model_curve import CurveSession, model_curves, plan_curves  # noqa: E402
    pyazr.kinds()
except Exception as err:
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)

warnings.simplefilter("ignore", pyazr.EngineWarning)
tmp = tempfile.mkdtemp(prefix="pyazr_model_curve_")


def copy(source, name):
    target = os.path.join(tmp, name)
    shutil.copytree(source, target, ignore=shutil.ignore_patterns("output", "checks", ".*.ui*"))
    os.makedirs(os.path.join(target, "output"), exist_ok=True)
    return target


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300))) if a.size else 0.0


def at_points(curve, data):
    """The curve's values at the data set's points (the grid entries marked at_data with the point's lab energy and,
    for an angle branch or a fixed angle, its c.m. angle)."""
    out = []
    for i in range(data.n):
        hit = curve.at_data & (curve.energy_lab == data.extra["energy_lab"][i])
        if curve.variable == "angle" or "angle_lab" in data.extra:
            hit &= np.isclose(curve.angles, data.angles[i], rtol=1e-9, atol=0)
        k = np.flatnonzero(hit)
        out.append(curve.model[k[0]] if k.size else np.nan)
    return np.asarray(out)


def folded(record):
    return any(k in record.fields for k in ("target", "beam", "detector"))


# -- 1. at the data points -------------------------------------------------------------------------------------------
print("at the data set's points")
t1 = copy(os.path.join(ROOT, "tests", "upgrade", "T1"), "T1")
s1 = pyazr.Session(os.path.join(t1, "T1.azr"))
ev1 = s1.evaluate()
c1 = model_curves(s1)
for d in ev1.data:
    c = c1[d.name]
    check(f"T1 {d.name}: a curve", not isinstance(c, Exception), str(c))
    if not isinstance(c, Exception):
        check(f"T1 {d.name}: at its points to rounding", rel(at_points(c, d), d.model) <= 1e-12,
              rel(at_points(c, d), d.model))

rr = copy(os.path.join(ROOT, "tests", "ratio_to_rutherford", "models"), "rutherford")
srr = pyazr.Session(os.path.join(rr, "aC_data.azr"))
d = srr.evaluate().data[0]
c = model_curves(srr)[d.name]
check("ratio to Rutherford: seven angular distributions", c.variable == "angle" and len(c.branches) == 7,
      f"{c.variable} {len(c.branches)}")
check("ratio to Rutherford: at its points to rounding", rel(at_points(c, d), d.model) <= 1e-12,
      rel(at_points(c, d), d.model))
check("ratio to Rutherford: an angle branch has one lab energy", all(
    np.all(c.energy_lab[b.start:b.stop] == b.value) and b.fixed == "energy_lab" for b in c.branches))

t2 = copy(os.path.join(ROOT, "tests", "upgrade", "T2"), "T2")
azr2 = os.path.join(t2, "T2.azr")
s2 = pyazr.Session(azr2)
ev2 = s2.evaluate()
records2 = pyazr.Project.open(azr2).records()
c2 = model_curves(s2)
kinds_seen = set()
for d in ev2.data:
    c, r = c2[d.name], records2[d.name]
    if d.name == "s07_profile":
        continue
    if isinstance(c, Exception):
        check(f"T2 {d.name}: a curve", False, str(c))
        continue
    if not folded(r):
        bound, kind = 1e-12, "unfolded"
    elif "shift" not in r.fields and "thm" not in r.fields:
        bound, kind = 1e-8, "folded"
    else:
        bound, kind = 2e-4, "folded with a shift or THM"
    kinds_seen.add(kind)
    got = rel(at_points(c, d), d.model)
    check(f"T2 {d.name} ({kind}): at its points to {bound:g}", got <= bound, got)
check("T2: unfolded, folded and shifted data sets among them", len(kinds_seen) == 3, kinds_seen)
check("T2: a convolved (s02: beam spread and detector), a target-integrated (s04 thick, s01 thin) data set made",
      all(not isinstance(c2[n], Exception) for n in ("s02_diff090", "s04_thick", "s01_int_thin")))

# -- 2. between the points: the same points read as data -------------------------------------------------------------
print("between the points")


def as_data(base, azr, name, curve, data, branch=0):
    """A copy of the project where the data set ``name`` reads a data file of the curve's grid energies (lab) of
    ``branch`` at the angle of the data set's first point (``data``: its EvaluatedData); its model through the data
    path."""
    target = copy(base, f"{name}_grid")
    b = curve.branches[branch]
    project = pyazr.Project.open(os.path.join(target, os.path.basename(azr)))
    record = project[name]
    columns = record.fields["columns"].split()
    lines = []
    for e in curve.energy_lab[b.start:b.stop]:
        cells = []
        for c in columns:
            role = c.split("[")[0]
            cells.append(repr(float(e)) if role == "E_lab" else
                         repr(float(data.extra["angle_lab"][0])) if role == "theta_lab" else
                         repr(float(data.angles[0])) if role == "theta_cm" else "1.0" if role == "value" else
                         "0.1" if role == "error" else "0")
        lines.append(" ".join(cells))
    with open(os.path.join(target, "grid.dat"), "w") as f:
        f.write("\n".join(lines) + "\n")
    fields = {"file": "grid.dat", "energy_window": None}
    project.set(name, **fields)
    project.save()
    s = pyazr.Session(project.path)
    return s.evaluate().dataset(name)


for base, azr, name, curves, ev in ((t2, azr2, "s13_thm_oes", c2, ev2),
                                    (t1, os.path.join(t1, "T1.azr"), "cap_g0", c1, ev1),
                                    (t1, os.path.join(t1, "T1.azr"), "el_lab090", c1, ev1),
                                    (t1, os.path.join(t1, "T1.azr"), "cap_g0_lab090", c1, ev1)):
    c = curves[name]
    d = as_data(base, azr, name, c, ev.dataset(name))
    b = c.branches[0]
    check(f"{name}: the curve's {b.stop - b.start} grid points read as data give the curve",
          rel(c.model[b.start:b.stop], d.model) <= 1e-12, rel(c.model[b.start:b.stop], d.model))

# -- 3. energy shifts ----------------------------------------------------------------------------------------------
print("energy shifts")
project = pyazr.Project.open(azr2)
project.set("s13_thm_oes", shift="offset 0.004 MeV fixed")
shifted = os.path.join(t2, "T2_shifted.azr")
project.save(shifted)
ss = pyazr.Session(shifted)
x = np.array(ss.parameters.x0, float)
moved = [p for p in ss.parameters.free if p.name in ("G.offset", "G.a0")]
check("T2: the shared energy scale G is free", len(moved) == 2, [p.name for p in moved])
for p in moved:
    x[p.free_index] += 0.002
evs = ss.evaluate(x)
cs = model_curves(ss, ["s13_thm_oes", "s02_diff090"], x=x)
d = evs.dataset("s13_thm_oes")
check("a fixed shift: the points' lab energies moved", abs(d.extra["energy_lab"][0] - ev2.dataset("s13_thm_oes")
                                                            .extra["energy_lab"][0] - 0.004) < 1e-9)
check("a fixed shift on an unfolded data set: at its points to rounding",
      rel(at_points(cs["s13_thm_oes"], d), d.model) <= 1e-12)
d = evs.dataset("s02_diff090")
check("a fitted shift moved in x on a folded data set: at its points to 1e-4",
      rel(at_points(cs["s02_diff090"], d), d.model) <= 1e-4, rel(at_points(cs["s02_diff090"], d), d.model))
check("the curve's x is the x given", np.array_equal(cs["s02_diff090"].x, x))

# -- 4. the grid -----------------------------------------------------------------------------------------------------
print("the grid")
c = model_curves(s2, ["s09_ratio"], points=50)["s09_ratio"]
d = ev2.dataset("s09_ratio")
own = np.unique(d.extra["energy_lab"])
check("points=50: 50 evenly spaced values and the data set's own", c.n == len(np.unique(np.concatenate(
    [own, np.linspace(own[0], own[-1], 50)]))), c.n)
check("the grid is sorted and spans the first to the last point", np.all(np.diff(c.energy_lab) > 0)
      and c.energy_lab[0] == own[0] and c.energy_lab[-1] == own[-1])
check("at_data marks the data set's points", np.array_equal(np.sort(c.energy_lab[c.at_data]), own))
check("one branch of an angle-integrated quantity", len(c.branches) == 1 and c.branches[0].fixed is None)
c = c1["cap_g0_lab090"]
check("capture at a lab angle (its c.m. angle varies with energy): fixed angle_lab 90",
      c.branches[0].fixed == "angle_lab" and c.branches[0].value == 90.0, c.branches)
c = c1["el_lab090"]
check("elastic scattering at a lab angle (one c.m. angle at every energy too): fixed angle_lab 90, the file's frame",
      c.branches[0].fixed == "angle_lab" and c.branches[0].value == 90.0, c.branches)
c = c1["ay_cm090"]
check("an analyzing power at a c.m. angle: fixed angle_cm 90", c.branches[0].fixed == "angle_cm"
      and np.allclose(c.angles, 90.0))

# -- 5. refusals -----------------------------------------------------------------------------------------------------
print("refusals")
e = c2["s07_profile"]
check("a beam-profile window per point: refused with three parts", isinstance(e, ParameterError)
      and "window_lo" in e.what and e.allowed and e.fix, repr(e))
for label, call, cls in (
        ("an unknown data set", lambda: s2.model_curve("s99"), ParameterKeyError),
        ("a profile data set by the one-name call", lambda: s2.model_curve("s07_profile"), ParameterError),
        ("points=1", lambda: s2.model_curve("s09_ratio", points=1), ParameterError)):
    try:
        call()
        check(f"{label}: refused", False, "no error")
    except cls as err:
        check(f"{label}: refused with three parts", bool(err.what and err.allowed and err.fix), repr(err))
with pyazr.Session(os.path.join(t1, "T1.azr"), extrapolations=True) as sx:
    try:
        model_curves(sx)
        check("a session of the extrapolations: refused", False)
    except ParameterError as err:
        check("a session of the extrapolations: refused with three parts", bool(err.what and err.allowed and err.fix))
plan = plan_curves(s1, ["el_lab090", "ay_cm090"])
plan.records["ay_cm090"][0]["angles"] = "list 90.0 deg lab"  # what the engine does not compute: A_y at a lab angle
with CurveSession(s1, plan) as cs1:
    out = cs1.evaluate()
check("an engine refusal of one record refuses its data set with the engine's words",
      isinstance(out["ay_cm090"], ParameterError) and "lab angle" in out["ay_cm090"].what, repr(out["ay_cm090"]))
check("... and the others are made", not isinstance(out["el_lab090"], Exception))

# -- 6. the one-name call; ENDF's curves take a curve ----------------------------------------------------------------
print("the one-name call")
c = s2.model_curve("s13_thm_oes")
check("Session.model_curve equals model_curves at an unfolded data set",
      np.array_equal(c.model, c2["s13_thm_oes"].model) and np.array_equal(c.energies, c2["s13_thm_oes"].energies))
c = s2.model_curve("s02_diff090")
d = ev2.dataset("s02_diff090")
check("Session.model_curve of a folded data set alone (its own mesh): at its points to 3e-3",
      rel(at_points(c, d), d.model) <= 3e-3, rel(at_points(c, d), d.model))
from pyazr.endf import CurveRequest, request_of  # noqa: E402
r = request_of(s2, c2["s13_thm_oes"])
check("pyazr.endf.request_of takes a curve: its points and lab energies", isinstance(r, CurveRequest)
      and len(r.energies) == c2["s13_thm_oes"].n and np.array_equal(r.energies_lab, c2["s13_thm_oes"].energy_lab))

for s in (s1, s2, srr, ss):
    s.close()
shutil.rmtree(tmp, ignore_errors=True)
print("ok" if not failures else f"FAILED: {len(failures)}: {failures}")
sys.exit(1 if failures else 0)
