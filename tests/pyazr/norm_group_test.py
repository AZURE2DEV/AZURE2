#!/usr/bin/env python3
"""Normalization groups as pyazr sees them.

A format-3 ``normalization`` record is shared by the data sets whose ``norm`` field names it (tests/norm_groups,
docs/upgrade/project_format3.md): the engine gives each data set its own norm parameter, the first data set's free
(the record's value, mode and prior), every later one fixed and tied to it, so the parameter layout does not depend
on the groups.  pyazr must read the records, refuse a data set naming no normalization, and -- with the engine --
list the follower's norm as fixed, penalise the group once and difference the shared norm exactly as AZURE2 does.

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
    from pyazr import Project, Nuisance, FormatError  # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(0)

PROJECT = os.path.join(ROOT, "tests", "norm_groups")
AZR = os.path.join(PROJECT, "norm_groups.azr")

print("1. the group as records")
proj = Project.open(AZR)
check("data sets xs_40_s08 and xs_90_s09 name normalization D",
      proj["xs_40_s08"]["norm"] == "D" and proj["xs_90_s09"]["norm"] == "D",
      (proj["xs_40_s08"]["norm"], proj["xs_90_s09"]["norm"]))
check("D is a normalization record, 1.2 fixed",
      proj["D"].record_kind == "normalization" and proj["D"]["value"] == "1.2 fixed", proj["D"].fields)
check("data files still read", proj["xs_40_s08"]["file"] == "data/xs_40_s08.dat"
      and proj["xs_90_s09"]["file"] == "data/xs_90_s09.dat")
try:
    proj.set("xs_90_s09", norm="E")              # names no normalization record
    check("a norm naming no normalization is refused", False, "no error")
except FormatError as err:
    check("a norm naming no normalization is refused", proj["xs_90_s09"]["norm"] == "D", str(err))

print("\n2. with the engine: parameters, penalty, chi2 and Jacobian")
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
        # The group's norm: 1.2, free, prior 10 %.
        proj.set("D", value=Nuisance(1.2, "fit", prior=(1.2, "10%")))
        azr = os.path.join(work, "t.azr")
        proj.save(azr)
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
