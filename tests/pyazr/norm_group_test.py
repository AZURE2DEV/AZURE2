#!/usr/bin/env python3
"""Normalization groups as pyazr sees them.

A <segmentsData> line may end with ``normgroup <label>`` (see include/SegLine.h,
tests/norm_groups): the segments carrying one label share the normalization of
the first of them (its value, vary flag and penalty).  A follower whose own line
varies its norm keeps that parameter, fixed, so the parameter layout does not
depend on the groups.  pyazr must read the block, write it, and -- with the
engine -- list the follower's norm as fixed, penalise the group once and
difference the shared norm exactly as AZURE2 does.

Run from anywhere:  python3 tests/pyazr/norm_group_test.py
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


try:                                           # pyazr itself needs numpy
    import numpy as np                         # noqa: E402
    from pyazr.datasets import SegmentSet      # noqa: E402
    from pyazr.azrfile import AzrModel         # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(0)

PROJECT = os.path.join(ROOT, "tests", "norm_groups")
AZR = os.path.join(PROJECT, "norm_groups.azr")

print("1. the normgroup block is read")
s = {x.key: x for x in SegmentSet.from_file(AZR)}
check("segments 1 and 2: group D", s[1].norm_group == "D" and s[2].norm_group == "D",
      (s[1].norm_group, s[2].norm_group))
check("no shift group", s[1].shift_group is None and s[2].shift_group is None)
check("data files still parsed", s[1].data_file == "data/xs_40_s08.dat" and s[2].data_file == "data/xs_90_s09.dat",
      (s[1].data_file, s[2].data_file))

print("\n2. the block round-trips through AzrModel")
m = AzrModel.from_file(AZR)
check("byte-identical re-write", m.to_text() == open(AZR).read())
m.add_data_segment("data/xs_90_s08.dat", 1, 1, observable="differential-cm", angle_min=90, angle_max=90,
                   norm_group="E", shift_group="H", natural=True)
lines = [l for l in m.to_text().splitlines() if l.strip()]
last = lines[[i for i, l in enumerate(lines) if l.strip() == "</segmentsData>"][0] - 1]
check("add_data_segment writes normgroup before the sample block",
      last.split()[-6:] == ["shiftgroup", "H", "normgroup", "E", "sample", "1"], last)
back = SegmentSet._parse(last.split(), 99)
check("... and reads back", back.norm_group == "E" and back.shift_group == "H" and back.natural,
      (back.norm_group, back.shift_group, back.natural))
try:
    SegmentSet._parse(last.split()[:-3], 99)     # "... normgroup" with no label
    check("a normgroup without a label is refused", False, "no error")
except ValueError as err:
    check("a normgroup without a label is refused", "normgroup without a label" in str(err), str(err))

print("\n3. with the engine: parameters, penalty, chi2 and Jacobian")
try:
    from pyazr import azure2
    have_engine = True
except Exception as e:                              # pragma: no cover
    print(f"  skip  engine not built ({e})")
    have_engine = False
if have_engine:
    work = tempfile.mkdtemp(prefix="pyazr_normgroup_")
    try:
        shutil.copytree(os.path.join(PROJECT, "data"), os.path.join(work, "data"))
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        # Leader: norm 1.2, free, 10 %.  Follower: its own line "1, free, 10 %" -- ignored.
        text = open(AZR).read() \
            .replace("4 1.2 0 0 0 0 0 data/xs_40_s08.dat", "4 1.2 1 10 0 0 0 data/xs_40_s08.dat") \
            .replace("4 1 0 0 0 0 0 data/xs_90_s09.dat", "4 1 1 10 0 0 0 data/xs_90_s09.dat")
        azr = os.path.join(work, "t.azr")
        open(azr, "w").write(text)
        expected = {}
        for l in open(os.path.join(PROJECT, "expected", "chiSquared.out")):
            t = l.strip().split(",")
            if len(t) == 5 and t[0].isdigit():
                expected[int(t[0])] = float(t[1])
        with azure2(azr, cwd=work) as mdl:
            norms = mdl.parameters.norms
            check("two norm parameters (segments 1, 2): the layout is the ungrouped one",
                  [p.segment_key for p in norms] == [1, 2], [(p.name, p.segment_key) for p in norms])
            free = [p for p in norms if not p.fixed]
            check("only the leader's is free", [p.segment_key for p in free] == [1], [(p.segment_key, p.fixed) for p in norms])
            idx = free[0].free_index
            x0 = np.asarray(mdl.params_rwa, float)
            c = np.asarray(mdl.segment_chi2(x0), float)
            check("segment chi2 = the CLI's at the shared 1.2 (expected/chiSquared.out)",
                  all(abs(c[k - 1] - expected[k]) < 1e-5 * max(1.0, expected[k]) for k in expected),
                  [(k, c[k - 1], expected[k]) for k in expected])
            x = x0.copy()
            x[idx] = 1.32                                   # 1 sigma above the leader's nominal 1.2
            pen = mdl.penalties(x)
            check("penalty ((n - 1.2)/0.12)^2 = 1, once, on the leader",
                  abs(float(pen["norm"][0]) - 1.0) < 1e-9 and float(pen["norm"][1]) == 0.0, pen["norm"])
            check("objective = chi2 + every penalty",
                  abs(mdl.objective(x) - (float(np.sum(mdl.calculate_chi2_rwa(x))) + sum(
                      float(np.sum(pen[k])) for k in ("norm", "shift", "shift_sqrt", "shift_fn")))) < 1e-6)
            c2 = np.asarray(mdl.segment_chi2(x), float)
            check("moving the leader's norm moves both segments",
                  abs(c2[0] - c[0]) > 1.0 and abs(c2[1] - c[1]) > 1.0, (c, c2))
            r, J = mdl.residual_jacobian(x0)
            h = 1e-6
            xp, xm = x0.copy(), x0.copy()
            xp[idx] += h
            xm[idx] -= h
            fd = (mdl.residual_jacobian(xp)[0] - mdl.residual_jacobian(xm)[0]) / (2 * h)
            col = J[:, idx]
            n1 = len(r) // 2
            scale = max(float(np.max(np.abs(fd))), 1e-12)
            check("Jacobian column of the shared norm matches finite differences",
                  float(np.max(np.abs(col - fd))) < 1e-4 * scale, (float(np.max(np.abs(col - fd))), scale))
            check("... and covers the follower's rows", float(np.max(np.abs(col[n1:]))) > 1e-3 * scale,
                  float(np.max(np.abs(col[n1:]))))
    finally:
        shutil.rmtree(work, ignore_errors=True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): " + ", ".join(failures))
    sys.exit(1)
print("all checks passed")
