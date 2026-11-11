#!/usr/bin/env python3
"""pyazr.predictions(path): the expression trees of a project's data sets and extrapolations (docs/upgrade/data_model.md
section 3), the same lines as `AZURE2 --print-predictions` (AZURE2_BIN), on T2 (SUM, RATIO, sample, phase shift, ...),
natural_neon (sample SUMs) and the format-3 project ratio_to_rutherford/models/aC_story (a data set and an extrapolation
labelled by their record names, as every format-3 project's; an old project's by their numbers); a project that cannot
be set up raises ValueError.
Needs the compiled engine; skips (exit 77) without it.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    import pyazr
    from pyazr import _azure2  # noqa: F401
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

binary = os.environ.get("AZURE2_BIN")
NAMED = "tests/ratio_to_rutherford/models/aC_story.azr"
for rel in ("tests/upgrade/T2/T2.azr", "tests/natural_neon/natural_neon.azr", NAMED):
    d, name = os.path.split(os.path.join(ROOT, rel))
    cwd = os.getcwd()
    os.chdir(d)
    try:
        lines = pyazr.predictions(name)
    finally:
        os.chdir(cwd)
    check(f"{rel}: one line per data set and extrapolation", len(lines) > 0 and all(": " in l for l in lines), lines[:3])
    format3 = any(l.strip() == "format=3" for l in open(os.path.join(d, name)))
    if rel == NAMED:
        check(f"{rel}: the data set and the extrapolation by their record names",
              [l.split(": ")[0] for l in lines] == ["data aC_RTH", "extrapolation aC_RTH_curve"], lines)
    elif format3:
        check(f"{rel}: a format-3 project's data sets by their record names",
              all(not l.split(": ")[0].split(" ", 1)[1].isdigit() for l in lines), lines[:3])
    else:
        check(f"{rel}: an old project's data sets by their numbers",
              all(l.split(": ")[0].split(" ", 1)[1].isdigit() for l in lines), lines[:3])
    if binary:
        cli = subprocess.run([binary, "--print-predictions", name], cwd=d, stdout=subprocess.PIPE, check=True)
        check(f"{rel}: equal to AZURE2 --print-predictions", cli.stdout.decode().splitlines() == lines)
with tempfile.TemporaryDirectory() as t:
    bad = os.path.join(t, "missing.azr")
    try:
        pyazr.predictions(bad)
        check("a project that cannot be read raises ValueError", False)
    except ValueError:
        check("a project that cannot be read raises ValueError", True)
if failures:
    print(f"{len(failures)} FAILURE(S)")
    sys.exit(1)
print("all passed")
