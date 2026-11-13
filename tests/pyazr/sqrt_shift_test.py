#!/usr/bin/env python3
"""The sqrt(E) energy-shift term as pyazr sees it.

A format-3 data set may carry ``shift sqrt <b> MeV^1/2 <mode> [prior ...]`` (docs/upgrade/project_format3.md):
E' = E + offset + b*sqrt(E/MeV).  pyazr must read the field, write it on a new data set, and -- with the engine --
count its parameter and its penalty exactly as AZURE2 does.

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
    from pyazr import Project
except Exception as err:                       # numpy or the engine missing
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)

AZR = os.path.join(ROOT, "tests", "energy_shift_sqrt", "energy_shift_sqrt.azr")


def shifts(record):
    v = record.get("shift", [])
    return v if isinstance(v, list) else [v]


print("1. the sqrt term is read, and only from the data sets that carry it")
proj = Project.open(AZR)
names = [n for n, r in proj.records().items() if r.record_kind == "dataset"]
check("five data sets", names == ["xs_40", "xs_40_preshifted", "ay_60", "ay_60_preshifted", "xs_90"], names)
check("xs_40: offset 0.004 MeV, sqrt 0.003 MeV^1/2, fixed",
      shifts(proj["xs_40"]) == ["offset 0.004 MeV fixed", "sqrt 0.003 MeV^1/2 fixed"], shifts(proj["xs_40"]))
check("ay_60: sqrt -0.002 MeV^1/2", shifts(proj["ay_60"]) == ["sqrt -0.002 MeV^1/2 fixed"], shifts(proj["ay_60"]))
check("the other three: no sqrt term",
      all(not any(t.startswith("sqrt") for t in shifts(proj[n])) for n in ("xs_40_preshifted", "ay_60_preshifted", "xs_90")))
check("data file still read", proj["xs_40"]["file"] == "data/xs_40.dat", proj["xs_40"]["file"])

print("\n2. a new data set with the term, and one without")
spec = dict(reaction="1H(p,p)1H", quantity="dsigma_domega", unit="b/sr", file="data/xs_90.dat",
            columns="E_lab[MeV] theta_cm[deg] value error", energy_window="0 100 MeV", angle_window="90 90 deg",
            norm="1 fixed")
m = Project.open(AZR)
m.add("dataset", "xs_90_sqrt", shift="sqrt 0.001 MeV^1/2 fit prior 0 0.002", **spec)
check("the field is written on the new data set",
      "    shift          sqrt 0.001 MeV^1/2 fit prior 0 0.002" in m["xs_90_sqrt"].text, m["xs_90_sqrt"].text)
m.add("dataset", "xs_90_plain", **spec)
check("a plain data set has no shift", "shift" not in m["xs_90_plain"] and "shift" not in m["xs_90_plain"].text,
      m["xs_90_plain"].text)

print("\n3. with the engine: one parameter per data set, and the penalty")
try:
    from pyazr import Session
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
        proj.set("xs_40", shift=["offset 0.004 MeV fixed", "sqrt 0.003 MeV^1/2 fit prior 0 0.001"])
        azr = os.path.join(work, "t.azr")
        proj.save(azr)
        with Session(azr, cwd=work) as mdl:
            sq = mdl.parameters.sqrt_shifts
            check("five sqrt_shift parameters (one per data set)", len(sq) == 5, len(sq))
            free = [p for p in sq if not p.fixed]
            check("exactly segment 1's is free", [p.segment_key for p in free] == [1], [p.segment_key for p in free])
            check("named xs_40.shift.sqrt", free[0].name == "xs_40.shift.sqrt", free[0].name)
            check("value read from the file", abs(free[0].value - 0.003) < 1e-12, free[0].value)
            idx = free[0].free_index
            x0 = np.asarray(mdl.params_rwa, float)
            # The data chi-squared must respond to b: at the file value the
            # shifted segment 1 equals the pre-shifted segment 2.
            c = np.array([d.chi2 for d in mdl.evaluate(x0).data])
            check("shifted segment 1 == pre-shifted segment 2 (engine)",
                  abs(c[0] - c[1]) < 1e-6 * max(1.0, abs(c[1])), (c[0], c[1]))
            x = x0.copy()
            x[idx] = 0.005                        # 5 sigma from zero (2 from the file value)
            pen = mdl.evaluate(x)
            # measured from zero since 2026-10-25, the file value being only the start
            check("penalty (b/err)^2 = 25 (from zero)", abs(float(np.array([d.shift_sqrt_chi2 for d in pen.data]).sum()) - 25.0) < 1e-9, np.array([d.shift_sqrt_chi2 for d in pen.data]))
            check("objective = chi2 + penalties",
                  abs(mdl.objective()(x) - (mdl.evaluate(x).chi2 + float(np.array([d.norm_chi2 for d in pen.data]).sum())
                                           + float(np.array([d.shift_constant_chi2 for d in pen.data]).sum()) + float(np.array([d.shift_sqrt_chi2 for d in pen.data]).sum()))) < 1e-9)
            c2 = np.array([d.chi2 for d in mdl.evaluate(x).data])
            check("moving b moves segment 1 only", abs(c2[0] - c[0]) > 1.0 and
                  all(abs(c2[k] - c[k]) < 1e-9 for k in range(1, 5)), (c[0], c2[0]))
            # The Jacobian column of b: the API finite-differences every free
            # shift term; before 2026-10 a sqrt(E) column did not exist at all.
            r, J = mdl.objective().residuals_and_jacobian(x0)
            h = 1e-6
            xp, xm = x0.copy(), x0.copy()
            xp[idx] += h
            xm[idx] -= h
            fd = (mdl.objective().residuals_and_jacobian(xp)[0] - mdl.objective().residuals_and_jacobian(xm)[0]) / (2 * h)
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
