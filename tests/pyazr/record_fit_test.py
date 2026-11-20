#!/usr/bin/env python3
"""pyazr.record_fit.fit_record_settings: settings of a potential record fitted from outside the engine
(docs/upgrade/interior_potential.md section 10).

On tests/interior_capture/models/direct.azr (alpha + 12C, the potential model's E1 capture into the potential's own
0+ bound state) with `interior_charge_e1 0.02`: the model's own cross sections at 11 energies are written as data
(5 % errors), so the answer is known.
  (1) from interior_charge_e1 = 0.012 the fit returns 0.02 to 1e-6 with an objective below 1e-6;
  (2) the interval of delta = 1: at both ends an independent session gives the minimum + 1 to 1e-3, the ends lie
      on either side of 0.02, and the linear error from the Jacobian is half the interval to 5 % (the cross section is
      quadratic in the charge);
  (3) two settings, the interior charge and the depth of the potential (true values 0.02 and 45.891 MeV), from
      0.017 and 45.86.  The two are nearly degenerate in this one cross section, so the minimum is a long valley:
      the objective found is below 0.01, both values lie within 3 of their joint linear errors of the truth, and
      the charge's joint error is more than 1.5 times its error alone;
  (5) with a free parameter of the project's own (the data set's norm, prior 1 +- 5 %), so that every trial runs
      the engine's fit: from 0.016 the charge comes back to 0.02 within 1e-4 with the norm at 1 (1e-3), and the
      linear error is half the profile interval to 10 % (the fits leave scatter in the rows; the differences of
      the Jacobian must be far above it); each trial fits from the project's values and from the best vector so
      far and keeps the lower, and starts=("best",) alone gives the same charge here (one minimum);
  (6) the interval search keeps every profile value it has computed: with a profile that answers differently when
      asked twice at the same value (a second fit ending elsewhere), the crossing is still found (the root finder
      used to refuse the bracket), and a bound before the crossing gives None;
  (4) a start the engine refuses (a zero depth) raises with the engine's reason; the project file is unchanged
      and no trial file is left.
Needs the compiled engine and scipy; skips cleanly without them.

Run from anywhere:  python3 tests/pyazr/record_fit_test.py
"""
import hashlib
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "interior_capture", "models", "direct.azr")
REACTION = "12C(a,g)16O[4.734348545]"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import scipy  # noqa: F401
    from pyazr import Project, Session
    from pyazr.record_fit import Setting, fit_record_settings
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

work = tempfile.mkdtemp(prefix="record_fit_")
try:
    for d in ("output", "checks", "data"):
        os.makedirs(os.path.join(work, d))
    text = open(SRC).read().replace("step 0.1", "step 0.5")
    assert "    interior_boundary 0\n" in text and "<data>\n</data>" in text
    truth = os.path.join(work, "truth.azr")
    with open(truth, "w") as f:
        f.write(text.replace("    interior_boundary 0\n", "    interior_boundary 0\n    interior_charge_e1 0.02\n"))
    ev = Session(truth, cwd=work, extrapolations=True, ignore_externals=False).evaluate()
    curve = ev.data[0]
    with open(os.path.join(work, "data", "sigma.dat"), "w") as f:
        for e_cm, value in zip(curve.energies, curve.model):
            e_lab = e_cm * 16.0026 / 12.0  # the project's masses (4.0026 and 12 u), non-relativistic as the engine
            f.write(f"{e_lab:.10e} 0 {value:.10e} {0.05 * value:.10e}\n")
    dataset = ("<data>\ndataset synthetic\n    reaction       " + REACTION + "\n    quantity       sigma\n"
               "    unit           b\n    file           data/sigma.dat\n    columns        E_lab[MeV] _ value error\n"
               "    norm           1 fixed\n</data>")
    project = os.path.join(work, "fit.azr")
    with open(project, "w") as f:
        f.write(text.replace("    interior_boundary 0\n", "    interior_boundary 0\n    interior_charge_e1 0.012\n")
                .replace("<data>\n</data>", dataset))
    before = hashlib.sha256(open(project, "rb").read()).hexdigest()
    q = Setting("a12C", "interior_charge_e1", "{:.12g}", start=0.012, step=0.003)

    r = fit_record_settings(project, [q], intervals=True)
    check("(1) interior_charge_e1 from 0.012: 0.02 to 1e-6, objective below 1e-6",
          abs(r.values[0] - 0.02) < 1e-6 and r.objective < 1e-6 and r.converged,
          f"{r.values[0]!r} {r.objective!r} after {len(r.evaluations)} trials")
    lower, upper = r.intervals[0]
    at = []
    for end in (lower, upper):
        if end is None:
            at.append(float("nan"))
            continue
        ended = Project.open(project)
        ended.set("a12C", interior_charge_e1=f"{end:.12g}")
        f_end = os.path.join(work, "end.azr")
        ended.save(f_end)
        at.append(float(Session(f_end, cwd=work, ignore_externals=False).evaluate().objective))
    check("(2) the interval of delta = 1: the objective at both ends is the minimum + 1 (1e-3), 0.02 between them",
          lower is not None and upper is not None and lower < 0.02 < upper and all(abs(v - r.objective - 1.0) < 1e-3 for v in at),
          f"{lower!r} {upper!r} {at!r}")
    check("(2) the linear error of the setting is half the interval to 5 %",
          r.errors is not None and lower is not None and upper is not None and abs(r.errors[0] / (0.5 * (upper - lower)) - 1.0) < 0.05,
          f"{r.errors!r} {lower!r} {upper!r}")

    depth = Setting("a12C", "central", "woods_saxon depth {:.12g} MeV radius 3.55 fm diffuseness 0.6 fm",
                    start=45.86, step=0.02)
    q2 = Setting("a12C", "interior_charge_e1", "{:.12g}", start=0.017, step=0.002)
    r2 = fit_record_settings(project, [q2, depth], tolerance=1e-4)
    e2 = r2.errors
    check("(3) the interior charge and the depth together: objective below 0.01, both within 3 joint errors of the truth, "
          "the charge's joint error over 1.5 times its error alone",
          r2.objective < 1e-2 and e2 is not None and abs(r2.values[0] - 0.02) < 3 * e2[0] and abs(r2.values[1] - 45.891) < 3 * e2[1]
          and e2[0] > 1.5 * r.errors[0],
          f"{r2.values!r} {e2!r} {r2.objective!r} after {len(r2.evaluations)} trials")

    fitted = os.path.join(work, "fit_norm.azr")
    with open(fitted, "w") as f:
        f.write(open(project).read().replace("    norm           1 fixed\n", "    norm           1 fit prior 1 5%\n"))
    r5 = fit_record_settings(fitted, [Setting("a12C", "interior_charge_e1", "{:.12g}", start=0.016, step=0.002)], intervals=True)
    lo5, up5 = r5.intervals[0]
    warm_only = fit_record_settings(fitted, [Setting("a12C", "interior_charge_e1", "{:.12g}", start=0.016, step=0.002)],
                                    starts=("best",))
    check("(5) starts=('best',), half the fits: the same charge (1e-5)",
          abs(warm_only.values[0] - r5.values[0]) < 1e-5, f"{warm_only.values!r} {r5.values!r}")
    norm = float(r5.fit.x[0]) if r5.fit is not None and len(r5.fit.x) == 1 else float("nan")
    check("(5) with the data set's norm fitted in every trial: 0.02 (1e-4), norm 1 (1e-3), linear error = half the interval (10 %)",
          abs(r5.values[0] - 0.02) < 1e-4 and abs(norm - 1.0) < 1e-3 and r5.errors is not None and lo5 is not None and up5 is not None
          and abs(r5.errors[0] / (0.5 * (up5 - lo5)) - 1.0) < 0.10,
          f"{r5.values!r} norm {norm!r} errors {r5.errors!r} interval {lo5!r} {up5!r} after {len(r5.evaluations)} trials")

    from pyazr.record_fit import _crossing
    asked = {}

    def twice(v):  # the true rise the first time a value is asked, -1 (no rise) every time after
        asked[v] = asked.get(v, 0) + 1
        return (v - 1.0) ** 2 - 1.0 if asked[v] == 1 else -1.0

    up = _crossing(twice, 1.0, 0.3, -1.0, xtol=1e-9)
    down = _crossing(twice, 1.0, -0.3, -1.0, xtol=1e-9)
    check("(6) a profile that changes when asked twice: crossings at 2 and 0 (1e-6), no value asked twice, None at a bound",
          up is not None and down is not None and abs(up - 2.0) < 1e-6 and abs(down) < 1e-6 and max(asked.values()) == 1
          and _crossing(twice, 1.0, 0.3, -1.0, upper=1.5) is None, f"{up!r} {down!r} {max(asked.values())}")

    zero = Setting("a12C", "central", "woods_saxon depth {:.12g} MeV radius 3.55 fm diffuseness 0.6 fm", start=0.0, step=1.0)
    message = ""
    try:
        fit_record_settings(project, [zero])
    except RuntimeError as e:
        message = str(e)
    leftovers = [n for n in os.listdir(work) if ".recordfit." in n]
    check("(4) a zero depth raises with the engine's reason; the project is unchanged, no trial file is left",
          "refused" in message and hashlib.sha256(open(project, "rb").read()).hexdigest() == before and not leftovers,
          f"{message!r} {leftovers!r}")
finally:
    shutil.rmtree(work, ignore_errors=True)

print("FAILED: " + ", ".join(failures) if failures else "all checks passed")
sys.exit(1 if failures else 0)
