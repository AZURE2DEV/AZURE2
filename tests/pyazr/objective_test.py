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

This test reads the binding itself (`_azure2.Session`), so it uses the private handle `Session._engine`.
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
        x = np.asarray(m.parameters.x0, float)
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
        norm = m.parameters.of_kind("norm")[0]
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
                v = moved[next(it)] if not p.fixed else allrwa[i]
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
        subprocess.run([binary, "calculate", "--parameter-file=moved.sav", "--no-gui", "free_norm.azr"],
                       cwd=work, stdin=subprocess.DEVNULL, text=True,
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
        x = np.asarray(m.parameters.x0, float)
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
        r = subprocess.run([binary, "calculate", "--no-gui", "sg.azr"], cwd=work,
                           stdin=subprocess.DEVNULL, text=True, capture_output=True, timeout=900)
        total = None
        for l in r.stdout.splitlines():
            if l.startswith("Total Chi-Squared:"):
                total = float(l.split()[2])
        check("AZURE2 printed its total (with the penalties)", total is not None)
        if total is not None:
            check("pyazr objective = AZURE2's total", np.isclose(obj, total, rtol=2e-5), f"{obj} vs {total}")

# ---- 6. the engine's chi-squared call takes the data term from the objective --------------------------------------
# The binding's calculate_chi2_rwa(x) (AZUREAPI::CalculateChi2RWA) is the data chi-squared plus the terms of the model
# alone (Park's wall, the abundances, the model data, the gamma_decay priors).  Its data part must be the objective's
# own data term, not a second sum: with no segment prior active (every norm_chi2 and shift_chi2 is 0, asserted, so
# that the objective's first term holds the data alone), calculate_chi2_rwa(x) equals, bit for bit, the objective's
# terms "data, normalizations and energy shifts" + "data covariance", plus the model's terms added in the objective's
# order (each such term adds its value in one step, so the Python sum repeats the engine's).  Until 2026-11-16 the API
# summed the squares of the pulls, (d/e)^2, where the objective sums d^2/e^2: 11 of these 24 cases differed in the
# last bits (T1 Brune x0, 13N formal x0, ...).  T1-T3 (tests/upgrade) and 13N, three parametrizations, at the start
# and at a point with every free level energy and width moved by 1e-4 relative.
# Also: evaluate(x).objective = objective()(x) bit for bit.  Since 2026-11-17 the evaluation walks the terms once and
# takes its total from the breakdown's sum (each term from 0, then added in order), objective()(x) from the running
# sum; the two are one sequence of additions as long as every term after the first adds its value in one step.
print("\n6. calculate_chi2_rwa = the objective's data term (+ the model's terms), bit for bit")
MODEL_TERMS = ("Park J > 0", "abundances", "model data", "gamma decay priors")
with tempfile.TemporaryDirectory() as tmp:
    for folder, name in (("upgrade/T1", "T1.azr"), ("upgrade/T2", "T2.azr"), ("upgrade/T3", "T3.azr"),
                         ("13N", "13N.azr")):
        work = os.path.join(tmp, folder.replace("/", "_"))
        shutil.copytree(os.path.join(ROOT, "tests", folder), work)
        for junk in ("output", "checks"):
            shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
            os.makedirs(os.path.join(work, junk))
        for basis, opts in (("Brune", {}), ("formal", {"use_brune": False}), ("Park", {"use_park": True})):
            with Session(os.path.join(work, name), cwd=work, **opts) as m:
                x0 = np.asarray(m.parameters.x0, float)
                x1 = x0.copy()
                for p in m.parameters.free:
                    if p.kind in ("energy", "width"):
                        x1[p.free_index] = x0[p.free_index] * (1.0 + 1e-4)
                for where, x in (("x0", x0), ("moved", x1)):
                    ev = m.evaluate(x)
                    priors = sum(d.norm_chi2 + d.shift_chi2 for d in ev.data)
                    data = ev.terms["data, normalizations and energy shifts"] + ev.terms["data covariance"]
                    model = 0.0
                    for term, value in ev.terms.items():
                        if term in MODEL_TERMS:
                            model += value
                    chi2 = float(m._engine.calculate_chi2_rwa(x))
                    label = f"{folder} {basis} {where}"
                    check(f"{label}: no segment prior active", priors == 0.0, repr(priors))
                    check(f"{label}: calculate_chi2_rwa {chi2!r} = data {data!r} + model {model!r}",
                          chi2 == data + model, f"difference {chi2 - (data + model)!r}")
                    value = float(m.objective()(x))
                    check(f"{label}: evaluate(x).objective {ev.objective!r} = objective()(x) {value!r}",
                          ev.objective == value, f"difference {ev.objective - value!r}")

if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all objective checks passed")
