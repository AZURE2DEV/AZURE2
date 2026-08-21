#!/usr/bin/env python3
"""chi-squared per segment and per dataset, and AZURE2's own objective.

`calculate_chi2_rwa` returns one number, and it is the *data* chi-squared only.
Two things follow that used to be left to the caller, with a recipe in the skill
for each:

  * A fit is judged per experiment, not per segment, and an experiment often
    owns several segments. `segment_chi2` and `dataset_chi2` do the split.

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
    from pyazr import azure2, AzrModel
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)

with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(ROOT, "tests", "13N"), work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))

    # Free one normalization, so the penalty has something to act on: the
    # project as shipped fixes every norm at its nominal, where it is zero.
    model = AzrModel.from_file(os.path.join(work, "13N.azr"))
    model.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    azr = os.path.join(work, "free_norm.azr")
    model.write(azr)

    print("1. the per-segment split adds up")
    with azure2(azr, cwd=work) as m:
        x = np.asarray(m.params_rwa, float)
        seg = m.segment_chi2(x)
        total = float(np.sum(m.calculate_chi2_rwa(x)))
        check("one entry per segment", len(seg) == m.nsegments,
              f"{len(seg)} for {m.nsegments}")
        check("sums to calculate_chi2_rwa", np.isclose(seg.sum(), total, rtol=1e-9),
              f"{seg.sum()} vs {total}")

        print("\n2. datasets aggregate their segments")
        byfile = m.dataset_chi2(x)
        check("no dataset lost", np.isclose(sum(c for c, _, _ in byfile.values()),
                                            total, rtol=1e-9))
        check("point counts add up",
              sum(n for _, n, _ in byfile.values())
              == sum(len(m.energies[i]) for i in range(m.nsegments)))
        check("segment counts add up",
              sum(k for _, _, k in byfile.values()) == m.nsegments)

        print("\n3. at nominal normalizations the penalty is zero")
        check("no penalty", m.penalties(x)["norm"].sum() == 0.0)
        check("objective equals chi-squared there",
              np.isclose(m.objective(x), total, rtol=1e-12))

        print("\n4. a moved normalization is penalised as AZURE2 penalises it")
        norm = m.parameters.norms[0]
        moved = x.copy()
        moved[norm.free_index] = 1.10           # nominal 1.0, sys error 5%
        penalty = float(m.penalties(moved)["norm"].sum())
        check("((1.10 - 1.0) / 0.05)^2 = 4", np.isclose(penalty, 4.0, rtol=1e-9),
              str(penalty))
        chi2_moved = float(np.sum(m.calculate_chi2_rwa(moved)))
        check("objective = chi-squared + penalty",
              np.isclose(m.objective(moved), chi2_moved + penalty, rtol=1e-12))

        # Hand the same parameters to the engine, and let it report its own.
        allrwa = list(np.asarray(m.sess.params_all_rwa(), float))
        it = iter(range(len(moved)))
        sav = os.path.join(work, "moved.sav")
        with open(sav, "w") as fh:
            for i, p in enumerate(m.parameters):
                v = moved[next(it)] if not m.fixed_params[i] else allrwa[i]
                fh.write(f"{p.name:>28s} {float(v): .7e} {0.0: .7e}\n")

    import glob
    import subprocess
    binary = None
    for cand in sorted(glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))):
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
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all objective checks passed")
