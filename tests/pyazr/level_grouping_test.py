#!/usr/bin/env python3
"""AzrModel groups <levels> lines as the engine does (io::LevelGrouper; docs/upgrade/project.md P2).

On tests/13N, tests/17O and tests/upgrade/T2: with the lines of <levels> interleaved (every other line first, then
the rest, so that no two lines of a level stand next to each other), AzrModel finds the same levels as in the
original file -- the same (J, parity, energy, channel lines) -- and AZURE2 reads the interleaved file with the same
parameters.  (The former rule, consecutive lines with the same ID, doubles the level count of these files: 24, 8 and
10 levels instead of 12, 4 and 5.)
Needs the compiled engine; skips cleanly without it.  Run from anywhere:  python3 tests/pyazr/level_grouping_test.py
"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
CASES = [("13N", "13N.azr"), ("17O", "17O.azr"), ("upgrade/T2", "T2.azr")]
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    from pyazr import azure2, AzrModel
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(0)


def signature(model):
    return sorted((lv.J, lv.parity, lv.energy, tuple(sorted(c.to_line() for c in lv.channels)))
                  for lv in model.levels)


for project, azr in CASES:
    work = tempfile.mkdtemp(prefix="level_grouping.")
    try:
        d = os.path.join(work, project)
        shutil.copytree(os.path.join(ROOT, "tests", project), d)
        for sub in ("output", "checks"):
            os.makedirs(os.path.join(d, sub), exist_ok=True)
        orig = os.path.join(d, azr)
        text = open(orig).read().split("\n")
        a = next(i for i, l in enumerate(text) if l.strip() == "<levels>")
        b = next(i for i, l in enumerate(text) if l.strip() == "</levels>")
        body = [l for l in text[a + 1:b] if l.strip()]
        rev = os.path.join(d, "interleaved.azr")
        with open(rev, "w") as f:
            f.write("\n".join(text[:a + 1] + body[::2] + body[1::2] + text[b:]))
        m0, m1 = AzrModel.from_file(orig), AzrModel.from_file(rev)
        check(f"{project}: the same {len(m0.levels)} levels from the interleaved lines", signature(m0) == signature(m1),
              f"{len(m0.levels)} vs {len(m1.levels)} levels")
        s0, s1 = azure2(orig, cwd=d), azure2(rev, cwd=d)
        check(f"{project}: the engine reads the same parameters", sorted(np.asarray(s0.params_rwa).round(12)) ==
              sorted(np.asarray(s1.params_rwa).round(12)), "")
    finally:
        shutil.rmtree(work, ignore_errors=True)

print("PASS" if not failures else f"FAIL: {len(failures)}")
sys.exit(1 if failures else 0)
