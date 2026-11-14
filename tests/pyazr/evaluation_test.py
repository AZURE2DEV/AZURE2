#!/usr/bin/env python3
"""pyazr.Evaluation (Session.evaluate(x)) against the program's own output files, on T1, T2 and T3 (tests/upgrade).

- The program runs menu 1 (a calculation with the data) on a copy of the project; a session on another copy
  evaluates at its start x.
- chiSquared.out: one row per data set (chi2, N, norm, norm chi2) and its Total- line (Total-Chi-Squared is
  Evaluation.chi2; norm and shift totals the sums of data[i].norm_chi2, .shift_chi2; the other Total-* the terms of
  the engine's names) equal to the digits the program prints.
- AZUREOut_*: every data set's block (c.m. energy, excitation energy, c.m. angle, model, datum x norm, error x norm;
  the Legendre coefficients of an angular distribution) to the eleven digits printed.
- T2 at a displaced point (norms and shifts off their centres): the rows and totals of the chiSquared.out the
  session's writer (write_output_files) gives at that point.
- The objective is Session.objective(x) bit for bit and the sum of the terms; Evaluation.chi2 plus the model terms
  (Park, abundances, model data) is the binding's calculate_chi2_rwa(x); the model data are pyazr.model_data's (T3).
- The data sets' rows of the objective's rows (data[i].rows): consecutive from row 0, one per point, the squares
  adding up to chi2, and to each data set's chi2 where no covariance joins data sets (T1, T3); the three parts of
  the shift prior add up to shift_chi2.
- Extrapolations (T1 in extrapolation mode): the model, energies and S-factor conversion of every grid are the
  engine's forward pass bit for bit, and every chi-squared field is None.
- One evaluation per Evaluation: an evaluate(x) takes the time of one objective(x) (13N; the engine's
  Coulomb-function counters, cache_stats, stay at zero when the Coulomb table serves the functions).
- The self-description: pyazr.kinds() names Objective, Evaluation and EvaluatedData with their fields and the point
  outputs; help(pyazr) names objective and evaluate.

Run: python3 tests/pyazr/evaluation_test.py   (AZURE2_BIN: the program; ctest gives it)

This test reads the binding itself (`_azure2.Session`), so it uses the private handle `Session._engine`.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
TESTS = os.path.join(ROOT, "tests")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    import pyazr
    from pyazr import Session
    pyazr.kinds()
except Exception as err:  # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)
if not os.environ.get("AZURE2_BIN"):
    print("skip: AZURE2_BIN (the program) is not set")
    sys.exit(77)

REL = 1e-9  # AZUREOut_* columns: eleven significant digits
# chiSquared.out's Total- fields that are one term of the engine's objective
TOTAL_TERMS = {"Total-Park-Chi-Squared:": "Park J > 0", "Total-Abundance-Chi-Squared:": "abundances",
               "Total-Screening-Chi-Squared:": "screening priors", "Total-Density-Chi-Squared:": "density priors",
               "Total-Doppler-Chi-Squared:": "Doppler priors", "Total-Thickness-Chi-Squared:": "thickness priors",
               "Total-Nuisance-Chi-Squared:": "nuisance priors", "Total-Model-Data-Chi-Squared:": "model data",
               "Total-Prior-Covariance-Chi-Squared:": "prior covariance"}


def copy_project(name, dst):
    src = os.path.join(TESTS, "upgrade", name)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("output", "checks"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(dst, d), exist_ok=True)
    return os.path.join(dst, name + ".azr")


def run_program(azr):
    proc = subprocess.run([os.environ["AZURE2_BIN"], "--no-readline", os.path.basename(azr)], input="1\n\n\n\n7\n",
                          text=True, cwd=os.path.dirname(azr), capture_output=True, timeout=900)
    return os.path.join(os.path.dirname(azr), "output"), proc


def chi_squared_file(directory):
    """Rows (texts of chi2, N, norm, norm chi2) and the Total- line as {label: text}."""
    rows, totals = [], {}
    for line in open(os.path.join(directory, "chiSquared.out")):
        line = line.strip()
        if not line or line.startswith("Segment#"):
            continue
        if line.startswith("Total-"):
            tokens = line.split()
            totals = {tokens[k]: tokens[k + 1] for k in range(0, len(tokens) - 1, 2)}
            continue
        rows.append(line.split(",")[1:])
    return rows, totals


def printed(value, text, digits=6):
    """The program's text of a number (ostream precision `digits`, general format) against a value."""
    return f"{value:.{digits}g}" == text or float(f"{value:.{digits}g}") == float(text)


def blocks(directory):
    out = []
    for name in sorted(os.listdir(directory)):
        if not (name.startswith("AZUREOut_") and name.endswith(".out")) or "sample_seg" in name:
            continue
        rows = []
        for line in list(open(os.path.join(directory, name))) + [""]:
            cells = line.split()
            if not cells:
                if rows:
                    out.append((name, rows))
                rows = []
                continue
            rows.append([float(c) for c in cells])
    return out


def close(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return a.shape == b.shape and bool(np.all(np.abs(a - b) <= REL * np.abs(b) + 1e-300))


def matches(d, block):
    if len(block) != d.n:
        return False
    if "legendre" in d.extra:  # an angular distribution: E, then the coefficients of each point
        return all(close([d.energies[k]] + list(d.extra["legendre"][k]), block[k]) for k in range(d.n))
    b = np.array(block)
    if b.shape[1] < 7:
        return False
    return (close(d.energies, b[:, 0]) and close(d.excitation_energies, b[:, 1]) and close(d.angles, b[:, 2])
            and close(d.model, b[:, 3]) and close(d.values * d.norm, b[:, 5]) and close(d.errors * d.norm, b[:, 6]))


with tempfile.TemporaryDirectory() as tmp:
    for name in ("T1", "T2", "T3"):
        print(f"== {name}")
        program, proc = run_program(copy_project(name, os.path.join(tmp, name + "_cli")))
        rows, totals = chi_squared_file(program)
        s = Session(copy_project(name, os.path.join(tmp, name + "_api")))
        x = np.asarray(s.parameters.x0, float)
        ev = s.evaluate(x)
        check(f"{name}: one data set per row of chiSquared.out ({len(ev.data)})", len(rows) == len(ev.data),
              f"{len(rows)} rows")
        bad = [d.name for d, r in zip(ev.data, rows)
               if not (printed(d.chi2, r[0]) and int(r[1]) == d.n and printed(d.norm, r[2])
                       and printed(d.norm_chi2, r[3]))]
        check(f"{name}: every data set's chi2, N, norm, norm chi2 are the program's row", not bad, bad)
        check(f"{name}: Total-Chi-Squared {totals.get('Total-Chi-Squared:')} = Evaluation.chi2 ({ev.chi2:.10g})",
              printed(ev.chi2, totals["Total-Chi-Squared:"]))
        check(f"{name}: Total-Norm-Chi-Squared = the data sets' norm priors",
              printed(sum(d.norm_chi2 for d in ev.data), totals["Total-Norm-Chi-Squared:"]))
        check(f"{name}: Total-N", int(float(totals["Total-N:"])) == sum(d.n for d in ev.data))
        if "Total-Shift-Chi-Squared:" in totals:
            check(f"{name}: Total-Shift-Chi-Squared = the data sets' shift priors",
                  printed(sum(d.shift_chi2 for d in ev.data), totals["Total-Shift-Chi-Squared:"]))
        for label, term in TOTAL_TERMS.items():
            if label in totals:
                digits = 15 if "Prior-Covariance" in label else 6
                check(f"{name}: {label} {totals[label]} = terms[{term!r}]",
                      printed(ev.terms[term], totals[label], digits), ev.terms[term])
        found = blocks(program)
        missing = [d.name for d in ev.data if not any(matches(d, b) for _, b in found)]
        check(f"{name}: every data set's arrays are a block of the program's AZUREOut_* files", not missing, missing)
        check(f"{name}: objective = Session.objective()(x) bit for bit ({ev.objective:.10g})",
              ev.objective == s.objective()(x))
        total = sum(ev.terms.values())
        check(f"{name}: the terms add up to the objective", abs(total - ev.objective) <= 1e-12 * max(1, ev.objective))
        model_terms = ev.terms["Park J > 0"] + ev.terms["abundances"] + ev.terms["model data"]
        old = float(s._engine.calculate_chi2_rwa(x))  # the binding's data chi-squared with the model terms
        check(f"{name}: chi2 + the model terms = the binding's calculate_chi2_rwa(x) ({old:.10g})",
              abs(ev.chi2 + model_terms - old) <= 1e-12 * max(1, old))
        if name == "T2":
            # away from the start (norms and shifts off their centres): the writer's chiSquared.out at that x
            xd = np.array([v * (1.0 + 0.02 * (1 if k % 2 else -1)) for k, v in enumerate(x)])
            evd = s.evaluate(xd)
            out = s.write_output_files(xd)
            rows_d, totals_d = chi_squared_file(out)
            bad = [d.name for d, r in zip(evd.data, rows_d)
                   if not (printed(d.chi2, r[0]) and printed(d.norm, r[2]) and printed(d.norm_chi2, r[3]))]
            check("T2 displaced: every row is the writer's chiSquared.out at that x (norm priors non-zero: "
                  f"{sum(d.norm_chi2 for d in evd.data):.4g})", not bad and sum(d.norm_chi2 for d in evd.data) > 0, bad)
            check("T2 displaced: Total-Chi-Squared, Total-Norm-, Total-Shift-Chi-Squared",
                  printed(evd.chi2, totals_d["Total-Chi-Squared:"])
                  and printed(sum(d.norm_chi2 for d in evd.data), totals_d["Total-Norm-Chi-Squared:"])
                  and printed(sum(d.shift_chi2 for d in evd.data), totals_d["Total-Shift-Chi-Squared:"]))
            check("T2: data sets have the lab energy, conversion and no Legendre list",
                  all("energy_lab" in d.extra and "sfactor_conversion" in d.extra for d in ev.data))
        # the data sets' rows of the objective's rows, and the three parts of the shift prior
        rows = s.objective().residuals(x)
        spans = [(d.rows.start, d.rows.stop) for d in ev.data]
        check(f"{name}: the data sets' rows follow one another from row 0, one per point",
              spans[0][0] == 0 and all(a[1] == b[0] for a, b in zip(spans, spans[1:]))
              and all(d.rows.stop - d.rows.start == d.n for d in ev.data), spans[:3])
        whitened = [float(np.sum(rows[d.rows] ** 2)) for d in ev.data]
        check(f"{name}: the squares of the data rows add up to chi2 (1e-10)",
              abs(sum(whitened) - ev.chi2) <= 1e-10 * max(1.0, ev.chi2), (sum(whitened), ev.chi2))
        if name != "T2":  # T2 joins two data sets by a covariance: chi2 splits it as chiSquared.out, the rows otherwise
            check(f"{name}: per data set, the squares of its rows are its chi2 (1e-10)",
                  all(abs(w - d.chi2) <= 1e-10 * max(1.0, d.chi2) for w, d in zip(whitened, ev.data)))
        check(f"{name}: the three parts of the shift prior add up to shift_chi2",
              all(abs(d.shift_constant_chi2 + d.shift_sqrt_chi2 + d.shift_function_chi2 - d.shift_chi2)
                  <= 1e-12 * max(1.0, d.shift_chi2) for d in ev.data))
        if name == "T3":
            want = pyazr.model_data(s, x)
            check(f"T3: the model data are pyazr.model_data's ({len(want)})", want and ev.model_data == want)
        check(f"{name}: table() names every data set", all(str(d.name) in ev.table() for d in ev.data))
        s.close()

    print("== the extrapolations (T1 in extrapolation mode)")
    work = os.path.join(tmp, "T1x")
    shutil.copytree(os.path.join(TESTS, "upgrade", "T1"), work)
    os.makedirs(os.path.join(work, "output"), exist_ok=True)
    with Session(os.path.join(work, "T1.azr"), extrapolations=True) as s:
        x = np.asarray(s.parameters.x0, float)
        ev = s.evaluate(x)
        n = int(s._engine.update_segments_rwa(x))  # the engine's forward pass, through the binding
        same = n == len(ev.data) and all(
            np.array_equal(d.model, s._engine.calculated_segments(i))
            and np.array_equal(d.energies, s._engine.calculated_energies(i))
            and np.array_equal(d.extra["sfactor_conversion"], s._engine.calculated_conv(i)) for i, d in enumerate(ev.data))
        check(f"the {n} grids: model, energies, S-factor conversion are the forward pass's bit for bit", same)
        check("no chi-squared: chi2, objective, terms, data[i].chi2, .norm, .rows are None or empty",
              ev.chi2 is None and ev.objective is None and ev.terms == {} and not ev.has_data
              and all(d.chi2 is None and d.norm is None and d.rows is None and d.residuals is None for d in ev.data))
        check("the grids carry their names", all(d.name for d in ev.data), [d.name for d in ev.data])

    print("== one evaluation per Evaluation (13N)")
    work = os.path.join(tmp, "13N")
    os.makedirs(os.path.join(work, "output"))
    os.symlink(os.path.join(TESTS, "13N", "data"), os.path.join(work, "data"))
    shutil.copy(os.path.join(TESTS, "13N", "13N.azr"), work)
    with Session(os.path.join(work, "13N.azr")) as s:
        x = np.asarray(s.parameters.x0, float)
        s.objective()(x)
        s.evaluate(x)  # both paths warm

        t = time.perf_counter()
        for _ in range(3):
            s.objective()(x)
        t_obj = (time.perf_counter() - t) / 3
        t = time.perf_counter()
        for _ in range(3):
            s.evaluate(x)
        t_ev = (time.perf_counter() - t) / 3
        check(f"evaluate takes the time of one objective ({t_ev * 1e3:.1f} ms against {t_obj * 1e3:.1f} ms)",
              t_ev < 1.3 * t_obj)

    print("== the self-description")
    k = pyazr.kinds()
    check("kinds() names the objects and their fields",
          {"Objective", "Evaluation", "EvaluatedData"} <= set(k["objects"])
          and "log_prob(x)" in k["objects"]["Objective"]["fields"] and "terms" in k["objects"]["Evaluation"]["fields"])
    check("kinds()['data_outputs'] lists the point outputs with units",
          {"energy_lab", "angle_lab", "model_e1", "legendre"} <= set(k["data_outputs"])
          and k["data_outputs"]["energy_lab"]["unit"] == "MeV")
    check("help(pyazr) names objective and evaluate", "s.objective()" in pyazr.__doc__
          and "evaluate(x)" in pyazr.__doc__)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok: Evaluation")
