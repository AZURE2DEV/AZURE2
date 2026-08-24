#!/usr/bin/env python3
"""AzrModel.apply_fit must put every fitted value on the channel it belongs to.

Two identities have to line up for that, and both used to be got wrong:

  * **Which level.** A level at Ex = 0 comes back from the API with
    ``level_energy = None``, so a key rebuilt from (J, parity, energy) cannot
    tell it from anything else -- and a J-group holding a background pole then
    reloads scaled rather than misassigned. The engine's own ``(jgroup, level)``
    is unambiguous; ``engine_level_keys`` reproduces it from file order.

  * **Which pair.** A ``Parameter``'s ``pair`` is the engine's number, which
    counts pairs in the order ``<levels>`` first mentions them. That is *not*
    the pair key the file writes. On the 8Be model engine pair 1 is file key 2
    and file key 1 is engine pair 6, so matching one against the other wrote
    every width to the wrong channel: 146 of 202 parameters, silently.

The pair half is tested against synthetic Parameter/Pair objects rather than a
project file, because every .azr in tests/ happens to have an identity mapping
and so cannot catch it. The level half and the end-to-end snapshot use the real
models, and are skipped when the compiled engine is not built.

Run from anywhere:  python3 tests/pyazr/apply_fit_test.py
"""
import glob
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

spec = importlib.util.spec_from_file_location(
    "azrfile", os.path.join(ROOT, "pyazr", "azrfile.py"))
azrfile = importlib.util.module_from_spec(spec)
spec.loader.exec_module(azrfile)
AzrModel = azrfile.AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


PROJECTS = sorted(glob.glob(os.path.join(ROOT, "tests", "*", "*.azr")))

# --------------------------------------------------------------------------
print("1. engine_level_keys numbers levels the way the engine does")
model = AzrModel.from_file(PROJECTS[0])
keys = model.engine_level_keys()
check("every level got a key", len(keys) == len(model.levels),
      f"{len(keys)} keys for {len(model.levels)} levels")
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


_pspec = importlib.util.spec_from_file_location(
    "parameters", os.path.join(ROOT, "pyazr", "parameters.py"))
_params = importlib.util.module_from_spec(_pspec)
_pspec.loader.exec_module(_params)


def parameter(**kw):
    """A Parameter with the fields apply_fit reads, others defaulted."""
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


model = AzrModel.from_file(PROJECTS[0])
level = model.engine_level_keys()[(1, 1)]
targets = [(c.pair, c.L, c.S) for c in level.channels]
if len(targets) < 2:
    check("the fixture has two channels to tell apart", False,
          f"{len(targets)} channel(s)")
else:
    # Claim the engine numbered the pairs the other way round: engine pair 1 is
    # the file's second channel-pair, and vice versa.
    a_pair, b_pair = targets[0][0], targets[1][0]
    pairs = [pair(number=1, key=b_pair), pair(number=2, key=a_pair)]
    p = parameter(pair=1, L=targets[1][1], S=targets[1][2])
    model.apply_fit([p], [4.2e3], physical=True, pairs=pairs)
    written = {(c.pair, c.L, c.S): c.gamma for c in level.channels}
    check("value landed on the pair the key names",
          written[targets[1]] == 4.2e3,
          f"got {written[targets[1]]} on {targets[1]}")
    check("and not on the pair the engine numbered",
          written[targets[0]] != 4.2e3)

# --------------------------------------------------------------------------
print("\n3. a parameter with no home is refused, not skipped")
model = AzrModel.from_file(PROJECTS[0])
orphan = parameter(pair=99, L=7, S=0.5, jgroup=1, level=1)
try:
    model.apply_fit([orphan], [1.0], physical=True)
    check("strict raised", False, "it wrote a partial file instead")
except ValueError as err:
    check("strict raised", True)
    check("the message names the parameter", "width_1_1" in str(err))
try:
    model.apply_fit([orphan], [1.0], physical=True, strict=False)
    check("strict=False writes the rest", model.applied == 0
          and len(model.unplaced) == 1)
except ValueError:
    check("strict=False writes the rest", False, "it still raised")

# --------------------------------------------------------------------------
print("\n4. rwa values are refused unless the caller says they are physical")
model = AzrModel.from_file(PROJECTS[0])
try:
    model.apply_fit([parameter(pair=1, L=0, S=0.5)], [1.0])
    check("neither transform= nor physical= raises", False)
except ValueError:
    check("neither transform= nor physical= raises", True)

# --------------------------------------------------------------------------
print("\n5. a snapshot reads back as the fit it was made from")
try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2
except Exception as err:                                  # engine not built
    print(f"  skip  engine not available ({type(err).__name__})")
else:
    project = os.path.join(ROOT, "tests", "hybrid_potential")
    azr = os.path.join(project, "hybrid_potential.azr")
    # Beside the original, not in a temp dir: a .azr resolves its data files
    # and output directory relative to itself, so that is the only place a
    # snapshot of it runs.
    out = os.path.join(project, "_snapshot_test.azr")
    sav = os.path.splitext(out)[0] + ".sav"
    try:
        with azure2(azr, cwd=project) as m:
            want = np.asarray(m.transform_rwa(m.params_rwa), float)
            written, savfile = m.save_fit(out)
            check("snapshot written", os.path.exists(written))
            check("companion param.sav written", savfile and os.path.exists(savfile))
            check("param.sav has every parameter",
                  sum(1 for _ in open(savfile)) == len(m.parameters))
        with azure2(out, cwd=project) as back:
            got = np.asarray(back.transform_rwa(back.params_rwa), float)
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
print("all apply_fit checks passed")
