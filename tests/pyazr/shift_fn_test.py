#!/usr/bin/env python3
"""Energy-shift functions as pyazr sees them.

A format-3 ``energy_scale`` record may hold ``function "<f(x)>" MeV`` and one ``param`` per coefficient (tests/shift_fn,
docs/upgrade/project_format3.md): E' = E + offset + f(x; a0, ...), x the lab energy in MeV, for every data set whose
``shift`` names the record.  pyazr must read the records, write new ones, refuse a data set with two terms of one
kind, and -- with the engine -- list the parameters (kind "shift_fn"), penalise them and difference them exactly as
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
    from pyazr import Project, FormatError     # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(0)

PROJECT = os.path.join(ROOT, "tests", "shift_fn")
AZR = os.path.join(PROJECT, "shift_fn.azr")


def shifts(record):
    v = record.get("shift", [])
    return v if isinstance(v, list) else [v]


print("1. the functions are read, with the group's record")
proj = Project.open(AZR)
check("xs_40_2 names xs_40_2_scale: f = a0*sqrt(x), a0 = 0.003 fixed, no penalty",
      "xs_40_2_scale" in shifts(proj["xs_40_2"]) and proj["xs_40_2_scale"]["function"] == '"a0*sqrt(x)" MeV'
      and proj["xs_40_2_scale"]["param"] == "a0 0.003 fixed", proj["xs_40_2_scale"].fields)
check("xs_40: a sqrt term, no function", shifts(proj["xs_40"]) == ["offset 0.004 MeV fixed", "sqrt 0.003 MeV^1/2 fixed"],
      shifts(proj["xs_40"]))
check("xs_40_lin_mislabelled: group G, f = a0*(x-1)",
      shifts(proj["xs_40_lin_mislabelled"]) == ["G"] and proj["G"]["function"] == '"a0*(x-1)" MeV',
      (shifts(proj["xs_40_lin_mislabelled"]), proj["G"].fields))
check("xs_90_lin_mislabelled: group G, no function of its own", shifts(proj["xs_90_lin_mislabelled"]) == ["G"])
check("data files still read", proj["xs_40_2"]["file"] == "data/xs_40.dat"
      and proj["xs_90_lin_mislabelled"]["file"] == "data/xs_90_lin_mislabelled.dat")

print("\n2. new records: a scale with a function, a data set with its own offset and a shared function")
spec = dict(reaction="1H(p,p)1H", quantity="dsigma_domega", unit="b/sr", file="data/xs_90.dat",
            columns="E_lab[MeV] theta_cm[deg] value error", energy_window="0 100 MeV", angle_window="90 90 deg",
            norm="1 fixed")
m = Project.open(AZR)
m.add("energy_scale", "H", function='"a0*x^1.5+a1" MeV', param=["a0 0.001 fit prior 0 0.002", "a1 0 fixed"])
m.add("dataset", "xs_90_H", shift="H", **spec)
check("the function and its two coefficients are written",
      m["H"].text.split("\n")[1:] == ['    function       "a0*x^1.5+a1" MeV', "    param          a0 0.001 fit prior 0 0.002",
                                      "    param          a1 0 fixed"] and shifts(m["xs_90_H"]) == ["H"], m["H"].text)
m.add("energy_scale", "F", function='"a0*(sqrt(x)-sqrt(1.7))" MeV', param="a0 0 fit prior 0 0.02")
m.add("dataset", "xs_90_F", shift=["offset 0.001 MeV fit prior 0 0.01", "F"], **spec)
check("a data set with its own constant shift and the function of F",
      shifts(m["xs_90_F"]) == ["offset 0.001 MeV fit prior 0 0.01", "F"], shifts(m["xs_90_F"]))
try:
    m.add("dataset", "xs_90_GF", shift=["G", "F"], **spec)   # G holds an offset: two scales, two offsets
    check("two energy scales on one data set are refused", False, "no error")
except FormatError as err:
    check("two energy scales on one data set are refused", "E_FIELD_TWICE" in err.codes and "xs_90_GF" not in m,
          err.codes)

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
        proj.set("G", param="a0 0.005 fit prior 0 0.002")
        azr = os.path.join(work, "t.azr")
        proj.save(azr)
        expected = {}
        for l in open(os.path.join(PROJECT, "expected", "chiSquared.out")):
            t = l.strip().split(",")
            if len(t) == 5 and t[0].isdigit():
                expected[int(t[0])] = float(t[1])
        with azure2(azr, cwd=work) as mdl:
            fn = mdl.parameters.fn_shifts
            check("three shift_fn parameters (data sets 2, 3, 5)", [p.segment_key for p in fn] == [2, 3, 5],
                  [(p.name, p.segment_key) for p in fn])
            check("named <energy_scale>.a0, fn_term 0",
                  [p.name for p in fn] == ["xs_40_2_scale.a0", "xs_40_3_scale.a0", "G.a0"]
                  and all(p.fn_term == 0 for p in fn),
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
            x[idx] = 0.009                        # 4.5 sigma from zero (2 from the file value)
            pen = mdl.penalties(x)
            # measured from zero since 2026-10-25, the file value being only the start
            check("penalty (a0/0.002)^2 = 20.25 (from zero), on segment 5", abs(float(pen["shift_fn"][4]) - 20.25) < 1e-9
                  and abs(float(pen["shift_fn"].sum()) - 20.25) < 1e-9, pen["shift_fn"])
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
