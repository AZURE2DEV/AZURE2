#!/usr/bin/env python3
"""pyazr.Project groups <levels> rows as the engine does (io::LevelGrouper, through experiment::LevelEdits).

On tests/13N, tests/17O and tests/upgrade/T2: with the rows of <levels> interleaved (every other row first, then
the rest, so that no two rows of a level stand next to each other), Project.levels() finds the same levels as in the
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
    from pyazr import Session, Project
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)


def signature(project):
    return sorted((lv.J, lv.parity, lv.energy, tuple(sorted((c.pair, c.l, c.s, c.width_text, c.fixed) for c in lv.channels)))
                  for lv in project.levels())


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
        m0, m1 = Project.open(orig), Project.open(rev)
        check(f"{project}: the same {len(m0.levels())} levels from the interleaved lines", signature(m0) == signature(m1),
              f"{len(m0.levels())} vs {len(m1.levels())} levels")
        s0, s1 = Session(orig, cwd=d), Session(rev, cwd=d)
        check(f"{project}: the engine reads the same parameters", sorted(np.asarray(s0.parameters.x0).round(12)) ==
              sorted(np.asarray(s1.parameters.x0).round(12)), "")
    finally:
        shutil.rmtree(work, ignore_errors=True)

print("PASS" if not failures else f"FAIL: {len(failures)}")
sys.exit(1 if failures else 0)
