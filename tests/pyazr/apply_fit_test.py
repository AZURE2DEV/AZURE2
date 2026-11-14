#!/usr/bin/env python3
"""save_fit must put every fitted value on the row it belongs to (pyazr.snapshot.place_level_values, the
engine's experiment::LevelValueEdits through Project).

Two identities have to line up for that, and both used to be got wrong:

  * **Which level.** A level at Ex = 0 comes back from the API with
    ``level_energy = None``, so a key rebuilt from (J, parity, energy) cannot
    tell it from anything else -- and a J-group holding a background pole then
    reloads scaled rather than misassigned. The engine's own ``(jgroup, level)``
    is unambiguous; Project.levels() gives it for each level of the file.

  * **Which pair.** A ``Parameter``'s ``pair`` is the engine's number, which
    counts pairs in the order ``<levels>`` first mentions them. That is *not*
    the pair key the file writes. On the 8Be model engine pair 1 is file key 2
    and file key 1 is engine pair 6, so matching one against the other wrote
    every width to the wrong channel: 146 of 202 parameters, silently.

The pair half is tested against synthetic Parameter/Pair objects rather than a
session, because every project in tests/ happens to have an identity mapping
and so cannot catch it. The end-to-end snapshot uses a real model.  Needs the
compiled module (Project reads the file through the engine).

Run from anywhere:  python3 tests/pyazr/apply_fit_test.py
"""
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    from pyazr import Project, FormatError
    from pyazr import parameters as _params
    from pyazr.snapshot import place_level_values
    import pyazr._azure2  # noqa: F401
except Exception as err:                                  # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def first_project():
    for path in sorted(glob.glob(os.path.join(ROOT, "tests", "*", "*.azr"))):
        try:
            return path, Project.open(path)
        except (FormatError, ValueError):
            continue
    raise SystemExit("no project of tests/ opens")


PATH, _ = first_project()

# --------------------------------------------------------------------------
print("1. Project.levels() numbers levels the way the engine does")
model = Project.open(PATH)
keys = {(lv.group, lv.level): lv for lv in model.levels() if lv.active}
check("every active level got a key", len(keys) == len([lv for lv in model.levels() if lv.active]),
      f"{len(keys)} keys")
check("both indices are 1-based",
      all(j >= 1 and l >= 1 for j, l in keys))
# a J-group's levels must be numbered consecutively from 1, in file order
groups = {}
for (j, l), lv in sorted(keys.items()):
    groups.setdefault(j, []).append((l, lv.jpi))
check("levels consecutive within each group",
      all([n for n, _ in v] == list(range(1, len(v) + 1)) for v in groups.values()))
check("one J^pi per group",
      all(len({jpi for _, jpi in v}) == 1 for v in groups.values()))

# --------------------------------------------------------------------------
print("\n2. a fitted width follows the engine's pair number to the file's key")


def parameter(**kw):
    """A Parameter with the fields save_fit reads, others defaulted."""
    base = dict(index=0, name="width_1_1", kind="width", fixed=False,
                value=0.0, free_index=0, jgroup=1, level=1)
    base.update(kw)
    return _params.Parameter(**base)


def pair(number, key):
    """A Pair carrying just the engine number and file key that matter here."""
    return _params.Pair(number=number, key=key, ptype=0, is_entrance=False,
                        J1=0.5, parity1=1, Z1=1, M1=1.0,
                        J2=0.5, parity2=1, Z2=1, M2=1.0,
                        sep_energy=0.0, excitation=0.0, channel_radius=4.0,
                        i1i2factor=1.0)


def widths(project):
    return {(c.pair, c.l, c.s): c.width for c in project.level(keys[(1, 1)].label).channels}


model = Project.open(PATH)
targets = [(c.pair, c.l, c.s) for c in model.level(keys[(1, 1)].label).channels]
if len({t[0] for t in targets}) < 2:
    check("the fixture has two channel pairs to tell apart", False, f"{targets}")
else:
    # Claim the engine numbered the pairs the other way round: engine pair 1 is
    # the file's second channel-pair, and vice versa.
    second = next(t for t in targets if t[0] != targets[0][0])
    a_pair, b_pair = targets[0][0], second[0]
    pairs = [pair(number=1, key=b_pair), pair(number=2, key=a_pair)]
    # the repr shows BOTH numbers, since engine messages and isospin.out name pairs by key (2026-10-27)
    check("Pair repr shows the engine number and the project's key",
          f"Pair(#1 key={b_pair}" in repr(pairs[0]), repr(pairs[0]))
    p = parameter(pair=1, L=second[1], S=second[2])
    place_level_values(model, [p], pairs, [4.2e3])
    written = widths(model)
    check("value landed on the pair the key names",
          written[second] == 4.2e3,
          f"got {written[second]} on {second}")
    check("and not on the pair the engine numbered",
          written[targets[0]] != 4.2e3)

# --------------------------------------------------------------------------
print("\n3. a parameter with no home is refused, not skipped")
model = Project.open(PATH)
before = model.to_bytes()
orphan = parameter(pair=99, L=7, S=0.5, jgroup=1, level=1)
try:
    place_level_values(model, [orphan], [], [1.0])
    check("refused", False, "it wrote a partial file instead")
except ValueError as err:
    check("refused", True)
    check("the message names the parameter", "width_1_1" in str(err))

# A channel the level does not list has input 0, so the engine fixes its width (ZeroInputWidth) and rebuilds it from
# that zero on re-read; its value (0, formal-basis rounding, or gamma^int with --no-ignore-externals) is not written
# (docs/upgrade/external_capture.md section 9).  A free parameter with no row is still refused (above).
print("  a fixed width on a channel the level does not list is skipped, whatever its value")
for v in (0.0, 5.0):
    model = Project.open(PATH)
    unlisted = parameter(pair=99, L=7, S=0.5, jgroup=1, level=1, fixed=True, free_index=None)
    try:
        place_level_values(model, [unlisted], [], [v])
        check(f"fixed {v:g}, no row: nothing to place", model.to_bytes() == before)
    except ValueError as err:
        check(f"fixed {v:g}, no row: nothing to place", False, str(err).splitlines()[0])

# --------------------------------------------------------------------------
print("\n4. a snapshot reads back as the fit it was made from")
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np  # noqa: E402
from pyazr import Session  # noqa: E402
project = os.path.join(ROOT, "tests", "hybrid_potential")
azr = os.path.join(project, "hybrid_potential.azr")
# Beside the original, not in a temp dir: a .azr resolves its data files
# and output directory relative to itself, so that is the only place a
# snapshot of it runs.
out = os.path.join(project, "_snapshot_test.azr")
sav = os.path.splitext(out)[0] + ".sav"
try:
    with Session(azr, cwd=project) as m:
        want = np.asarray(m.transform(m.parameters.x0)[m.parameters.free_mask], float)
        written, savfile = m.save_fit(out)
        check("snapshot written", os.path.exists(written))
        check("companion param.sav written", savfile and os.path.exists(savfile))
        # the first line is the `#parametrization` tag, not a parameter; as a
        # comment line it is invisible to numpy.loadtxt, so positional readers
        # still see exactly the parameters
        check("param.sav has every parameter",
              sum(1 for l in open(savfile) if l.split() and not l.split()[0].startswith("#")) == len(m.parameters))
        check("param.sav tagged with its parametrization",
              open(savfile).readline().split()[0] == "#parametrization")
        check("numpy.loadtxt sees only the parameters",
              np.loadtxt(savfile, usecols=(1,)).size == len(m.parameters))
    with Session(out, cwd=project) as back:
        got = np.asarray(back.transform(back.parameters.x0)[back.parameters.free_mask], float)
    check("every R-matrix value round-trips",
          got.shape == want.shape and np.allclose(got, want, rtol=1e-4),
          f"{int(np.sum(~np.isclose(got, want, rtol=1e-4)))} differ"
          if got.shape == want.shape else f"{got.shape} vs {want.shape}")
finally:
    for f in (out, sav):
        if os.path.exists(f):
            os.remove(f)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all save_fit placement checks passed")
