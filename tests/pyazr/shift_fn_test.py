#!/usr/bin/env python3
"""Energy-shift functions as pyazr sees them.

A <segmentsData> line may end with ``shiftfn "<f(x)>" N v0 s0 vary0 ...`` (see
include/SegLine.h, tests/shift_fn): E' = E + shift + f(x; a0, ...), x the lab
energy in MeV, declared on a shift group's first segment for the whole group.
pyazr must read the block, write it, and -- with the engine -- list the
parameters (kind "shift_fn"), penalise them and difference them exactly as
AZURE2 does.

Run from anywhere:  python3 tests/pyazr/shift_fn_test.py
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

PROJECT = os.path.join(ROOT, "tests", "shift_fn")
AZR = os.path.join(PROJECT, "shift_fn.azr")

print("1. the shiftfn block is read, with the group label")
s = {x.key: x for x in SegmentSet.from_file(AZR)}
check("segment 2: f = a0*sqrt(x), a0 = 0.003 fixed, no penalty",
      s[2].shift_fn == "a0*sqrt(x)" and s[2].shift_fn_params == ((0.003, 0.0, False),),
      (s[2].shift_fn, s[2].shift_fn_params))
check("segment 1: sqrtshift, no function", s[1].shift_fn is None and abs(s[1].energy_shift_sqrt - 0.003) < 1e-12)
check("segment 5: group G, f = a0*(x-1)", s[5].shift_group == "G" and s[5].shift_fn == "a0*(x-1)",
      (s[5].shift_group, s[5].shift_fn))
check("segment 6: group G, no function of its own", s[6].shift_group == "G" and s[6].shift_fn is None)
check("data files still parsed", s[2].data_file == "data/xs_40.dat" and s[6].data_file == "data/xs_90_lin_mislabelled.dat",
      (s[2].data_file, s[6].data_file))
check("no composite flagged", all(not x.composite for x in s.values()))

print("\n2. the block round-trips through AzrModel")
m = AzrModel.from_file(AZR)
check("byte-identical re-write", m.to_text() == open(AZR).read())
m.add_data_segment("data/xs_90.dat", 1, 1, observable="differential-cm", angle_min=90, angle_max=90,
                   shift_group="H", shift_fn="a0*x^1.5+a1", shift_fn_params=[(0.001, 0.002, True), (0, 0, False)])
lines = [l for l in m.to_text().splitlines() if l.strip()]
last = lines[[i for i, l in enumerate(lines) if l.strip() == "</segmentsData>"][0] - 1]
check("add_data_segment writes the keyword blocks",
      last.split()[-11:] == ["shiftgroup", "H", "shiftfn", '"a0*x^1.5+a1"', "2",
                             "0.001", "0.002", "1", "0", "0", "0"], last)
back = SegmentSet._parse(last.split(), 99)
check("... and reads back", back.shift_fn == "a0*x^1.5+a1" and back.shift_group == "H"
      and back.shift_fn_params == ((0.001, 0.002, True), (0.0, 0.0, False)), (back.shift_fn, back.shift_fn_params))

print("\n3. with the engine: parameters, penalty, chi2 and Jacobian")
try:
    from pyazr import azure2
    have_engine = True
except Exception as e:                              # pragma: no cover
    print(f"  skip  engine not built ({e})")
    have_engine = False
if have_engine:
    work = tempfile.mkdtemp(prefix="pyazr_shiftfn_")
    try:
        shutil.copytree(os.path.join(PROJECT, "data"), os.path.join(work, "data"))
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        # Free a0 of the group's function (penalty width 0.002 around 0.005).
        text = open(AZR).read().replace('shiftfn "a0*(x-1)" 1 0.005 0 0', 'shiftfn "a0*(x-1)" 1 0.005 0.002 1')
        azr = os.path.join(work, "t.azr")
        open(azr, "w").write(text)
        expected = {}
        for l in open(os.path.join(PROJECT, "expected", "chiSquared.out")):
            t = l.strip().split(",")
            if len(t) == 5 and t[0].isdigit():
                expected[int(t[0])] = float(t[1])
        with azure2(azr, cwd=work) as mdl:
            fn = mdl.parameters.fn_shifts
            check("three shift_fn parameters (segments 2, 3, 5)", [p.segment_key for p in fn] == [2, 3, 5],
                  [(p.name, p.segment_key) for p in fn])
            check("named segment_<key>_energy_shift_fn0, fn_term 0",
                  [p.name for p in fn] == ["segment_2_energy_shift_fn0", "segment_3_energy_shift_fn0",
                                           "segment_5_energy_shift_fn0"] and all(p.fn_term == 0 for p in fn),
                  [(p.name, p.fn_term) for p in fn])
            free = [p for p in fn if not p.fixed]
            check("only segment 5's is free, at 0.005", [p.segment_key for p in free] == [5]
                  and abs(free[0].value - 0.005) < 1e-12, [(p.segment_key, p.value) for p in free])
            check("listed after every sqrt(E) coefficient",
                  min(p.index for p in fn) > max(p.index for p in mdl.parameters.sqrt_shifts))
            idx = free[0].free_index
            x0 = np.asarray(mdl.params_rwa, float)
            c = np.asarray(mdl.segment_chi2(x0), float)
            check("segment chi2 = the CLI's (expected/chiSquared.out)",
                  all(abs(c[k - 1] - expected[k]) < 1e-5 * max(1.0, expected[k]) for k in expected),
                  [(k, c[k - 1], expected[k]) for k in expected])
            x = x0.copy()
            x[idx] = 0.009                        # 2 sigma off the nominal
            pen = mdl.penalties(x)
            check("penalty ((a0 - 0.005)/0.002)^2 = 4, on segment 5", abs(float(pen["shift_fn"][4]) - 4.0) < 1e-9
                  and abs(float(pen["shift_fn"].sum()) - 4.0) < 1e-9, pen["shift_fn"])
            check("objective = chi2 + every penalty",
                  abs(mdl.objective(x) - (float(np.sum(mdl.calculate_chi2_rwa(x))) + sum(
                      float(np.sum(pen[k])) for k in ("norm", "shift", "shift_sqrt", "shift_fn")))) < 1e-9)
            c2 = np.asarray(mdl.segment_chi2(x), float)
            check("moving a0 moves the group (segments 5 and 6) and nothing else",
                  abs(c2[4] - c[4]) > 1.0 and abs(c2[5] - c[5]) > 1.0
                  and all(abs(c2[k] - c[k]) < 1e-9 for k in (0, 1, 2, 3, 6)), (c, c2))
            r, J = mdl.residual_jacobian(x0)
            h = 1e-6
            xp, xm = x0.copy(), x0.copy()
            xp[idx] += h
            xm[idx] -= h
            fd = (mdl.residual_jacobian(xp)[0] - mdl.residual_jacobian(xm)[0]) / (2 * h)
            col = J[:, idx]
            scale = max(float(np.max(np.abs(fd))), 1e-12)
            check("Jacobian column of a0 matches finite differences",
                  float(np.max(np.abs(col - fd))) < 1e-4 * scale and scale > 1.0,
                  (float(np.max(np.abs(col - fd))), scale))
    finally:
        shutil.rmtree(work, ignore_errors=True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): " + ", ".join(failures))
    sys.exit(1)
print("all checks passed")
