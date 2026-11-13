#!/usr/bin/env python3
"""chi-squared per segment and per dataset, and AZURE2's own objective.

`evaluate().chi2` returns one number, and it is the *data* chi-squared only.
Two things follow that used to be left to the caller, with a recipe in the skill
for each:

  * A fit is judged per experiment, not per segment, and an experiment often
    owns several segments. `evaluate(x).data[i].chi2` does the split.

  * AZURECalc::operator() does not minimize the data chi-squared. It also adds
    ((norm - nominal) / (nominal/100 * norm_error))^2 per segment, and the
    same for a free energy shift. Minimize the bare chi-squared and the
    normalizations drift to absorb every discrepancy. `objective` is what a
    hand-rolled Minuit should be pointed at.

The penalty is checked against AZURE2's own `Total-Norm-Chi-Squared`, not
against a transcription of the formula: the reference value comes from the
engine, so a change to either has to move both.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/objective_test.py
"""
import os
import shutil
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
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import Session, Project, Nuisance
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)

with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(ROOT, "tests", "13N"), work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))

    # Free one normalization, so the penalty has something to act on: the
    # project as shipped fixes every norm at its nominal, where it is zero.
    model = Project.open(os.path.join(work, "13N.azr"))
    model.set("meyer_84", norm=Nuisance(1, "fit", prior=(1, "5%")))
    azr = os.path.join(work, "free_norm.azr")
    model.save(azr)

    print("1. the per-segment split adds up")
    with Session(azr, cwd=work) as m:
        x = np.asarray(m.params_rwa, float)
        seg = np.array([d.chi2 for d in m.evaluate(x).data])
        total = m.evaluate(x).chi2
        check("one entry per segment", len(seg) == m.nsegments,
              f"{len(seg)} for {m.nsegments}")
        check("sums to the data chi2", np.isclose(seg.sum(), total, rtol=1e-9),
              f"{seg.sum()} vs {total}")

        print("\n2. the data sets hold the points of the segments")
        data = m.evaluate(x).data
        check("no data set lost", np.isclose(sum(d.chi2 for d in data), total, rtol=1e-9))
        check("point counts add up",
              sum(d.n for d in data) == sum(len(m.energies[i]) for i in range(m.nsegments)))
        check("one data set per segment", len(data) == m.nsegments)

        print("\n3. at nominal normalizations the penalty is zero")
        check("no penalty", np.array([d.norm_chi2 for d in m.evaluate(x).data]).sum() == 0.0)
        check("objective equals chi-squared there",
              np.isclose(m.objective()(x), total, rtol=1e-12))

        print("\n4. a moved normalization is penalised as AZURE2 penalises it")
        norm = m.parameters.norms[0]
        moved = x.copy()
        moved[norm.free_index] = 1.10           # nominal 1.0, sys error 5%
        penalty = float(np.array([d.norm_chi2 for d in m.evaluate(moved).data]).sum())
        check("((1.10 - 1.0) / 0.05)^2 = 4", np.isclose(penalty, 4.0, rtol=1e-9),
              str(penalty))
        chi2_moved = m.evaluate(moved).chi2
        check("objective = chi-squared + penalty",
              np.isclose(m.objective()(moved), chi2_moved + penalty, rtol=1e-12))

        # Hand the same parameters to the engine, and let it report its own.
        allrwa = list(np.asarray(m._engine.params_all_rwa(), float))
        it = iter(range(len(moved)))
        sav = os.path.join(work, "moved.sav")
        with open(sav, "w") as fh:
            for i, p in enumerate(m.parameters):
                v = moved[next(it)] if not m.fixed_params[i] else allrwa[i]
                fh.write(f"{p.name:>28s} {float(v): .7e} {0.0: .7e}\n")

    import glob
    import subprocess
    binary = os.environ.get("AZURE2_BIN")  # ctest's: the binary this build made
    for cand in [] if binary else sorted(glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))):
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            binary = cand
            break
    if binary is None:
        print("  skip  no AZURE2 binary to cross-check against")
    else:
        subprocess.run([binary, "--no-gui", "--no-readline", "free_norm.azr"],
                       cwd=work, input="1\nmoved.sav\n\n7\n", text=True,
                       capture_output=True, timeout=900)
        line = ""
        out = os.path.join(work, "output", "chiSquared.out")
        if os.path.exists(out):
            for l in open(out):
                if l.startswith("Total-Chi-Squared:"):
                    line = l
        check("AZURE2 produced a total", bool(line), "no chiSquared.out line")
        if line:
            fields = line.split()
            engine_norm = float(fields[fields.index("Total-Norm-Chi-Squared:") + 1])
            print(f"        engine: {line.strip()}")
            print(f"        pyazr : penalty {penalty:.6f}")
            check("pyazr's penalty is the engine's Total-Norm-Chi-Squared",
                  np.isclose(engine_norm, penalty, rtol=1e-3, atol=1e-6),
                  f"{engine_norm} vs {penalty}")

print()
# ---- 5. a shiftgroup follower is not charged a second shift penalty ----------------------------------------------
print("\n5. shiftgroup followers are not penalised twice")
with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "sg")
    shutil.copytree(os.path.join(ROOT, "tests", "shift_groups"), work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    # the energy scale G of both data sets with an offset of 0.004 MeV, free, prior 0 +- 0.01 MeV:
    # the engine penalises the group's shift once, (0.004/0.01)^2 = 0.16, on the first member
    sg = Project.open(os.path.join(work, "shift_groups.azr"))
    sg.set("G", offset=Nuisance(0.004, "fit", unit="MeV", prior=(0, 0.01)))
    azr = os.path.join(work, "sg.azr")
    sg.save(azr)
    with Session(azr, cwd=work) as m:
        x = np.asarray(m.params_rwa, float)
        pen = np.array([d.shift_constant_chi2 for d in m.evaluate(x).data])
        chi2 = m.evaluate(x).chi2
        obj = float(m.objective()(x))
        check("the group's shift penalty is charged once, on the first member",
              np.allclose(pen, [0.16, 0.0], atol=1e-12), str(pen.tolist()))
        check("objective = chi-squared + 0.16", np.isclose(obj, chi2 + 0.16, rtol=1e-12),
              f"{obj} vs {chi2} + 0.16")
    binary = os.environ.get("AZURE2_BIN")  # ctest's: the binary this build made
    for cand in [] if binary else sorted(glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))):
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            binary = cand
            break
    if binary is None:
        print("  skip  no AZURE2 binary to cross-check against")
    else:
        r = subprocess.run([binary, "--no-gui", "--no-readline", "sg.azr"], cwd=work,
                           input="1\n\n\n7\n", text=True, capture_output=True, timeout=900)
        total = None
        for l in r.stdout.splitlines():
            if l.startswith("Total Chi-Squared:"):
                total = float(l.split()[2])
        check("AZURE2 printed its total (with the penalties)", total is not None)
        if total is not None:
            check("pyazr objective = AZURE2's total", np.isclose(obj, total, rtol=2e-5), f"{obj} vs {total}")

if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all objective checks passed")
