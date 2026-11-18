#!/usr/bin/env python3
"""Every script of pyazr/examples runs by the command its docstring states, and prints one pinned number or text.

One row per run: the script, its arguments (``{out}`` is a fresh scratch directory of the row), the text its standard
output must contain, whether the row is ``slow`` (over about 4 s at one thread: run by ``--slow``, the ctest
``pyazr_examples_slow`` with the label ``slow``) and, for a row that cannot run here, the reason it is skipped (a
package that is not installed, the network).  Each script runs as a new process from the repository's root, as its
docstring says, with ``PYTHONPATH`` set to the root.  After all rows, nothing may have been written into the root
or beside a project a row named (a new file there fails the row's script): the scripts write only where ``--out`` or
a temporary directory says.

A new example is a file in pyazr/examples, a row here and a row of the page docs/source/user_guide/pyazr_examples.rst
(consistency check 58 fails without either).

Needs the compiled module.  Run from anywhere:  python3 tests/pyazr/examples_test.py [--slow | --all] [NAME ...]
"""
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import namedtuple

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

Row = namedtuple("Row", "script args pin slow skip needs", defaults=(False, None, ()))

T1, T2, T3 = "tests/upgrade/T1/T1.azr", "tests/upgrade/T2/T2.azr", "tests/upgrade/T3/T3.azr"

ROWS = [
    # build and edit a project
    Row("print_scheme", [T1], "1/2+  E_x = 2.3650 MeV  (3 channels, 4 free)"),
    Row("project_from_python", ["{out}"], "chi2 = 16.6119 on 33 points"),
    Row("edit_model", [T1, "--out", "{out}"],
        "edited  : 10 data sets, 450 points, 11 free parameters, chi2 = 392735.9470", slow=True),
    Row("exfor_fetch", [T1, "--entry", "T0001003.1", "--local", "tests/pyazr/data/exfor_master", "--out", "{out}"],
        "exfor_t0001003_1: 2 points, chi2 = 799.1"),
    Row("exfor_fetch", [T1, "--entry", "O2599004", "--target", "C-13", "--out", "{out}"], "",
        skip="the IAEA EXFOR and ENSDF services: no network call from a test (the row above runs it offline)"),
    Row("refusals", [T3], "codes    ['E_REACTION_NOT_ALLOWED', 'E_UNIT_NOT_ALLOWED']"),
    Row("nuisances", [T1, "el_lab090", "el_lab125", "el_lab160", "--out", "{out}"],
        "N_shared.value           0.996788 +- 0.00171", slow=True),
    Row("natural_target", ["tests/natural_target/models/M.azr"],
        "largest |x3 - (0.6 x1 + 0.4 x2)| / x3 = 0.00e+00"),
    Row("gamma_decay", ["tests/gamma_decay/models/record.azr"],
        "level5.total_width                  0.0384894 +- 0.000177 eV"),
    Row("thm_dataset", [T2], "norm (profiled) = 0.00196953 +- 1.8e-05"),
    Row("covariance_builder", [T1, "el_lab090", "--out", "{out}/T1_cov.azr"],
        "chi2 of el_lab090: 72.8685301 without, 47.34670333 with the covariance", slow=True),
    Row("prior_builder", [T1, "--out", "{out}/T1_prior.azr"],
        "prior on 1 parameter(s), parametrization 1 (Brune), from 1 part(s)"),
    # calculate and compare
    Row("per_dataset_chi2", [T2], "T2.azr: chi2 = 10905.2 over 330 points"),
    Row("model_curve", [T2], "s01_int_thin                b     40    257            0.39751          0.467377"),
    Row("decompose_cross_section", [T1, "--level", "1/2+ 2.365"], "1/2+ 2.365       7.841e-01"),
    Row("angular_distribution", [T1, "-n", "3"], "largest |a_k>0| over the range: 0.1889", slow=True),
    Row("channel_radius_scan", [T3, "--pair", "2", "--radii", "3.5", "4.0", "4.5", "--out", "{out}"],
        "     4.00        509.7     1.032"),
    Row("nuclear_potential", ["tests/hybrid_potential/hybrid_potential.azr"],
        "this session's parameters.x0         : 2464106.791", slow=True),
    Row("tabulate_cross_sections", [T1, "--out", "{out}/sigma.dat"], "259 knots", slow=True),
    # fit
    Row("fit_engine", [T1, "--out", "{out}"],
        "objective of the saved project as it stands: 462.273511 (the fit: 462.273511)", slow=True),
    Row("fit_iminuit", [T3], "MIGRAD: objective 509.31", needs=("iminuit",)),
    Row("fit_scipy", [T3], "least_squares: objective 509.31248", needs=("scipy",)),
    # uncertainties
    Row("fit_emcee", [T3, "--steps", "60", "--rate", "2", "1", "--draws", "5", "--out", "{out}/chain.npy"],
        "objective at the start 518.713: data 509.684", needs=("emcee",)),
    Row("sensitivities", [T1, "--check"], "worst column: E(1/2+)  rel. error 3.34e-04"),
    Row("uncertainty_band", [T1], "    0.4059       150.19            1.461     0.97%", slow=True),
    # rates
    Row("reaction_rate", [T1, "--entrance", "1", "--exit", "2"], "    0.1   5.59562e-07    0.0847", slow=True),
    Row("sfactor_extrapolation", [T1, "--entrance", "1", "--exit", "2", "--out", "{out}"],
        "peak S = 263.5 keV b at 0.4190 MeV"),
    # evaluations
    Row("endf_evaluation", [], "format check by pyazr.endf: 0 errors, 0 warnings"),
    Row("endf_start_model", ["{out}"], "MF=6/MT=102 (z,g) capture: 1081 lines", slow=True),
    # physics functions
    Row("channel_functions", [T1], " 0      1.0284      0.18764       2.0968      0.17374"),
    Row("ec_integrals", [T1], "20 external-capture pathways from pair 1", slow=True),
    Row("dimensionless_widths", [T1], "theta^2 = 0.281  vs  theta^2_formal = 1.102"),
    Row("transform_widths", [T1], "2.39511          33000            129255"),
]


def listing(directory):
    """The names in a directory (not recursive), for 'nothing was written here'."""
    try:
        return set(os.listdir(directory))
    except OSError:
        return set()


def main(argv):
    mode = "slow" if "--slow" in argv else "all" if "--all" in argv else "quick"
    names = [a for a in argv if not a.startswith("--")]
    rows = [r for r in ROWS if (mode == "all" or r.slow == (mode == "slow")) and (not names or r.script in names)]
    watched = {ROOT} | {os.path.dirname(os.path.join(ROOT, a)) for r in rows for a in r.args
                        if a.endswith(".azr") and os.path.isfile(os.path.join(ROOT, a))}
    before = {d: listing(d) for d in watched}
    env = dict(os.environ, PYTHONPATH=ROOT + os.pathsep + os.environ.get("PYTHONPATH", ""))
    env.setdefault("OMP_NUM_THREADS", "1")
    failures, ran, skipped, t_all = [], 0, 0, time.time()
    for r in rows:
        missing = [m for m in r.needs if importlib.util.find_spec(m) is None]
        reason = r.skip or (f"needs {', '.join(missing)}, not installed" if missing else None)
        if reason:
            print(f"  skip  {r.script}: {reason}")
            skipped += 1
            continue
        out = tempfile.mkdtemp(prefix=f"examples_{r.script}.")
        args = [a.replace("{out}", out) for a in r.args]
        t0 = time.time()
        p = subprocess.run([sys.executable, os.path.join("pyazr", "examples", r.script + ".py")] + args, cwd=ROOT,
                           env=env, capture_output=True, text=True, timeout=900)
        dt = time.time() - t0
        ok = p.returncode == 0 and r.pin in p.stdout
        print(f"  {'ok  ' if ok else 'FAIL'}  {r.script} {' '.join(r.args)}  ({dt:.1f} s)")
        if not ok:
            failures.append(r.script)
            print(f"        exit {p.returncode}; expected in the output: {r.pin!r}")
            print("        " + "\n        ".join(p.stdout.splitlines()[:3] + ["..."] +
                                                (p.stdout[-1500:] + p.stderr[-1500:]).splitlines()))
        shutil.rmtree(out, ignore_errors=True)
        ran += 1
    for d, names_before in before.items():
        for new in sorted(listing(d) - names_before):
            failures.append(f"{os.path.relpath(os.path.join(d, new), ROOT)} was written")
            print(f"  FAIL  a script wrote {os.path.relpath(os.path.join(d, new), ROOT)}")
    print(f"\n{mode} rows: {ran} run, {skipped} skipped, {time.time() - t_all:.1f} s; "
          + (f"{len(failures)} failure(s): {failures}" if failures else "all ok"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
