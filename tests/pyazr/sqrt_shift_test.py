#!/usr/bin/env python3
"""The sqrt(E) energy-shift term as pyazr sees it.

A <segmentsData> line may end with ``sqrtshift b bError vary`` (see
include/SegLine.h): E' = E + shift + b*sqrt(E/MeV).  pyazr must read the
block, write it, and -- with the engine -- count its parameter and its penalty
exactly as AZURE2 does.

Needs numpy and the compiled engine (importing the pyazr package loads
both); exits 77 (skipped) without them.

Run from anywhere:  python3 tests/pyazr/sqrt_shift_test.py
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
    import numpy as np
    from pyazr.datasets import SegmentSet
    from pyazr.azrfile import AzrModel
except Exception as err:                       # numpy or the engine missing
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)

AZR = os.path.join(ROOT, "tests", "energy_shift_sqrt", "energy_shift_sqrt.azr")

print("1. the sqrtshift block is read, and only from the lines that carry it")
segs = SegmentSet.from_file(AZR)
s = {x.key: x for x in segs}
check("segment 1: a = 0.004", abs(s[1].energy_shift - 0.004) < 1e-12, s[1].energy_shift)
check("segment 1: b = 0.003, fixed", abs(s[1].energy_shift_sqrt - 0.003) < 1e-12
      and not s[1].vary_shift_sqrt, (s[1].energy_shift_sqrt, s[1].vary_shift_sqrt))
check("segment 3: b = -0.002", abs(s[3].energy_shift_sqrt + 0.002) < 1e-12, s[3].energy_shift_sqrt)
check("segments 2, 4, 5: no sqrt term",
      all(s[k].energy_shift_sqrt == 0.0 and not s[k].vary_shift_sqrt for k in (2, 4, 5)))
check("data file still parsed", s[1].data_file == "data/xs_40.dat", s[1].data_file)
check("no composite flagged", all(not x.composite for x in segs))

print("\n2. the block round-trips through AzrModel")
m = AzrModel.from_file(AZR)
check("byte-identical re-write", m.to_text() == open(AZR).read())
m.add_data_segment("data/xs_90.dat", 1, 1, observable="differential-cm",
                   angle_min=90, angle_max=90, energy_shift_sqrt=0.001,
                   energy_shift_sqrt_error=0.002, vary_shift_sqrt=True)
last = [l for l in m.to_text().splitlines() if l.strip()]
last = last[[i for i, l in enumerate(last) if l.strip() == "</segmentsData>"][0] - 1]
check("add_data_segment writes the keyword block", last.split()[-4:] == ["sqrtshift", "0.001", "0.002", "1"], last)
m2 = AzrModel.from_file(AZR)
m2.add_data_segment("data/xs_90.dat", 1, 1, observable="differential-cm")
last2 = [l for l in m2.to_text().splitlines() if l.strip()]
last2 = last2[[i for i, l in enumerate(last2) if l.strip() == "</segmentsData>"][0] - 1]
check("a plain segment line is unchanged", "sqrtshift" not in last2, last2)

print("\n3. with the engine: one parameter per segment, and the penalty")
try:
    from pyazr import azure2
    have_engine = True
except Exception as e:                              # pragma: no cover
    print(f"  skip  engine not built ({e})")
    have_engine = False
if have_engine:
    work = tempfile.mkdtemp(prefix="pyazr_sqrt_")
    try:
        shutil.copytree(os.path.join(ROOT, "tests", "energy_shift_sqrt", "data"), os.path.join(work, "data"))
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        # Free the sqrt coefficient of segment 1 with a 0.001 penalty, so the
        # objective differs from the data chi-squared by a known amount.
        text = open(AZR).read().replace("sqrtshift 0.003 0 0", "sqrtshift 0.003 0.001 1")
        azr = os.path.join(work, "t.azr")
        open(azr, "w").write(text)
        with azure2(azr, cwd=work) as mdl:
            sq = mdl.parameters.sqrt_shifts
            check("five sqrt_shift parameters (one per segment)", len(sq) == 5, len(sq))
            free = [p for p in sq if not p.fixed]
            check("exactly segment 1's is free", [p.segment_key for p in free] == [1], [p.segment_key for p in free])
            check("named segment_1_energy_shift_sqrt", free[0].name == "segment_1_energy_shift_sqrt", free[0].name)
            check("value read from the file", abs(free[0].value - 0.003) < 1e-12, free[0].value)
            idx = free[0].free_index
            x0 = np.asarray(mdl.params_rwa, float)
            # The data chi-squared must respond to b: at the file value the
            # shifted segment 1 equals the pre-shifted segment 2.
            c = mdl.segment_chi2(x0)
            check("shifted segment 1 == pre-shifted segment 2 (engine)",
                  abs(c[0] - c[1]) < 1e-6 * max(1.0, abs(c[1])), (c[0], c[1]))
            x = x0.copy()
            x[idx] = 0.005                        # 2 sigma off the nominal
            pen = mdl.penalties(x)
            check("penalty ((b - b0)/err)^2 = 4", abs(float(pen["shift_sqrt"].sum()) - 4.0) < 1e-9, pen["shift_sqrt"])
            check("objective = chi2 + penalties",
                  abs(mdl.objective(x) - (float(np.sum(mdl.calculate_chi2_rwa(x))) + float(pen["norm"].sum())
                                           + float(pen["shift"].sum()) + float(pen["shift_sqrt"].sum()))) < 1e-9)
            c2 = mdl.segment_chi2(x)
            check("moving b moves segment 1 only", abs(c2[0] - c[0]) > 1.0 and
                  all(abs(c2[k] - c[k]) < 1e-9 for k in range(1, 5)), (c[0], c2[0]))
            # The Jacobian column of b: the API finite-differences every free
            # shift term; before 2026-10 a sqrt(E) column did not exist at all.
            r, J = mdl.residual_jacobian(x0)
            h = 1e-6
            xp, xm = x0.copy(), x0.copy()
            xp[idx] += h
            xm[idx] -= h
            fd = (mdl.residual_jacobian(xp)[0] - mdl.residual_jacobian(xm)[0]) / (2 * h)
            col = J[:, idx]
            scale = max(float(np.max(np.abs(fd))), 1e-12)
            check("Jacobian column of b matches finite differences",
                  float(np.max(np.abs(col - fd))) < 1e-4 * scale and scale > 1.0,
                  (float(np.max(np.abs(col - fd))), scale))
    finally:
        shutil.rmtree(work, ignore_errors=True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): " + ", ".join(failures))
    sys.exit(1)
if not have_engine:
    print("sections 1-2 passed; section 3 skipped (no engine)")
    sys.exit(77)
print("all checks passed")
