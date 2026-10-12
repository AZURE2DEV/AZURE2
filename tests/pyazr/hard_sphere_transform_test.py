#!/usr/bin/env python3
"""A project with hard-sphere J^pi groups (includeAllJPi=1) opened through pyazr.

pyazr's default is transform=True (AZURE2's --transform-params path), which runs
CNuc::TransformIn on every J-group.  A hard-sphere group has channels but no
levels, and the transform built a 0 x 0 eigenproblem for it, whose GSL workspace
is a null pointer: a segfault on opening the session.  The CLI checks never saw
it (tests/hard_sphere runs without the transform).

The check: tests/hard_sphere's no_dummies.azr (the option) and with_dummies.azr
(the same partial waves as zero-width dummy levels) must open and agree on
both routes, as they do bit for bit in the CLI.

Run from anywhere:  python3 tests/pyazr/hard_sphere_transform_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

try:
    import numpy as np
    from pyazr import azure2
    from pyazr.azrfile import AzrModel
    import pyazr._azure2  # noqa: F401
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(0)

MODELS = os.path.join(ROOT, "tests", "hard_sphere", "models")
OPTION = os.path.join(MODELS, "no_dummies.azr")
DUMMIES = os.path.join(MODELS, "with_dummies.azr")
assert "includeAllJPi=1" in open(OPTION).read()

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def evaluate(path, transform, where):
    os.makedirs(os.path.join(where, "output")); os.makedirs(os.path.join(where, "checks"))
    local = os.path.join(where, "m.azr")
    AzrModel.from_file(path).write(local)
    with azure2(local, data_mode=False, use_brune=False, ignore_externals=False, transform=transform) as m:
        x = np.array(m.params_rwa)
        return [np.array(seg) for seg in m.calculate_rwa(x)]


d = tempfile.mkdtemp(prefix="pyazr_hs_transform_")
try:
    for transform in (True, False):
        tag = "transform" if transform else "no transform"
        a = evaluate(OPTION, transform, os.path.join(d, f"opt_{int(transform)}"))
        b = evaluate(DUMMIES, transform, os.path.join(d, f"dum_{int(transform)}"))
        check(f"{tag}: the includeAllJPi model opens and evaluates",
              len(a) > 0 and all(np.all(np.isfinite(s)) and s.size > 0 for s in a), [s.size for s in a])
        same = len(a) == len(b) and all(s.shape == t.shape for s, t in zip(a, b))
        worst = max(np.max(np.abs(s / t - 1)) for s, t in zip(a, b)) if same else float("inf")
        check(f"{tag}: option == dummy levels (1e-10)", same and worst < 1e-10, worst)
finally:
    shutil.rmtree(d, ignore_errors=True)

print("\n" + ("FAILED: " + ", ".join(failures) if failures else "PASSED"))
sys.exit(1 if failures else 0)
