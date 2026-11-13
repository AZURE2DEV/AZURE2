#!/usr/bin/env python3
"""Results read before they are written (docs/upgrade/data_model.md, step S3 and check 16).

An evaluation writes its values into a result buffer (include/Results.h); under AZURE_POISON_RESULTS=1 every slot
of that buffer is a quiet NaN on entry to each evaluation, so a value read before this evaluation writes it -- a
value kept from an earlier evaluation, or one no evaluation ever writes -- is NaN in the output.  This test runs
the same calculations with and without the variable and requires equal results:

  CLI    T1, T2, T3: menu 1 (calculate with data) and menu 3 (extrapolate); T1 also a short fit (menu 2,
         --max-calls), which goes through the pooled working copies and the analytic gradient.  The model pass
         (data_model.md S4: Results::model and the objective's values, poisoned with the leaf values): menu 1 on
         tests/natural_neon (sample sums; T2 above has a SUM and a RATIO), and tests/transmission_yield's conv (a transmission and a yield under a Gaussian spread: the
         transform at every sub-point) and natT (a sample transmission).  Every output file must be byte for byte
         the same.
  rate   tests/reaction_rate's p12C model, menu 5 (reaction rate): each rate point evaluates into a buffer of its
         own (Unit::UsePrivateRecord), poisoned as every other; reactionrates.out byte for byte.
  pyazr  T1, T2, T3: the values, E1 and E2 parts, angles, energies and angular-distribution coefficients of
         every segment (AZUREAPI::UpdateSegments), each in its own process (the variable is read once).  The
         session owns the buffer it evaluates into (AZUREAPI::results()).

Every owner of a buffer is covered: the CLI calculation (menus 1 and 3), the fit's pooled working copies and the
iteration snapshot (menu 2), the reaction rate (menu 5) and the Python session.

Known case, accepted and listed: K1 -- the API returns the E1 and E2 parts of every point
(AZUREAPI::UpdateSegments), but the evaluation writes them only for an angle-integrated cross section
(GenMatrixFunc::CalculateBareCrossSection; for a point with a sub-point fold, only on an E1/E2 data line,
Unit::IntegrateTargetEffectComponents).  Elsewhere they are the 0 of a new record; poisoned, NaN.  pyazr reads
them only for E1/E2 lines, which are written.

    poison_results_test.py [--binary AZURE2] [projects]      (AZURE2_BIN / AZURE2_BINARY, PYAZR_ENGINE_DIR)

This test reads the binding itself (`_azure2.Session`), so it uses the private handle `Session._engine`.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))  # pyazr from this checkout
import upgrade_project as up  # noqa: E402

PROJECTS = {"T1": "T1.azr", "T2": "T2.azr", "T3": "T3.azr"}
# The model-pass projects (menu 1 only): name -> (directory relative to tests/, .azr); "transmission_yield" is made by
# its make_models.py into the scratch directory.
MODEL_PROJECTS = {"natural_neon": ("natural_neon", "natural_neon.azr"), "transmission": ("transmission_yield", "conv.azr"),
                  "natural_transmission": ("transmission_yield", "natT.azr")}
TESTS = os.path.normpath(os.path.join(HERE, ".."))
MENUS = {"calc": "1\n\n\n7\n", "extrap": "3\n\n\n7\n"}
FIT = ("T1", "2\nn\n\n\n7\n", ["--max-calls=30"])
RATE_DIR = os.path.normpath(os.path.join(HERE, "..", "reaction_rate"))
RATE_MENU = "5\n\n1\n2\nyes\ntemperatures.dat\n"  # no parameter file, entrance 1, exit 2, temperatures from a file


def run_rate(binary, poison, tmp):
    """Menu 5 on tests/reaction_rate/models/p12C.azr: the output files, as run_cli."""
    d = tempfile.mkdtemp(dir=tmp)
    shutil.copy(os.path.join(RATE_DIR, "models", "p12C.azr"), os.path.join(d, "t.azr"))
    shutil.copy(os.path.join(RATE_DIR, "temperatures.dat"), d)
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(d, sub))
    env = dict(os.environ, AZURE_POISON_RESULTS="1" if poison else "0")
    r = subprocess.run([binary, "--no-gui", "--no-readline", "--no-coulomb-table", "t.azr"], input=RATE_MENU, cwd=d,
                       env=env, capture_output=True, text=True, timeout=1800)
    out = {}
    for f in sorted(os.listdir(os.path.join(d, "output"))):
        with open(os.path.join(d, "output", f), "rb") as fh:
            out[f] = fh.read()
    return r.returncode, out


def run_cli(binary, project, menu, extra, poison, tmp, src=None, azr=None):
    d = tempfile.mkdtemp(dir=tmp)
    src = src or os.path.join(HERE, project)
    for f in os.listdir(src):
        p = os.path.join(src, f)
        if f in ("output", "checks"):
            continue
        (shutil.copytree if os.path.isdir(p) else shutil.copy)(p, os.path.join(d, f))
    for sub in ("output", "checks"):
        shutil.rmtree(os.path.join(d, sub), ignore_errors=True)
        os.makedirs(os.path.join(d, sub))
    env = dict(os.environ, AZURE_POISON_RESULTS="1" if poison else "0")
    if any(a.startswith("--max-calls") for a in extra):
        # A fit bounded by --max-calls is reproducible at one thread only: Minuit2's call counter is not atomic, so
        # with several threads the same fit can run one more iteration (docs/upgrade/fit_threads.md, cause B).
        env["OMP_NUM_THREADS"] = "1"
    r = subprocess.run([binary, "--no-gui", "--no-readline", "--no-coulomb-table", *extra, azr or PROJECTS[project]],
                       input=menu, cwd=d, env=env, capture_output=True, text=True, timeout=1800)
    out = {}
    for f in sorted(os.listdir(os.path.join(d, "output"))):
        with open(os.path.join(d, "output", f), "rb") as fh:
            out[f] = fh.read()
    return r.returncode, out


def api_values(project):
    """The API's arrays of every segment, as JSON (NaN as the string "nan")."""
    from pyazr.session import Session
    m = Session(os.path.join(HERE, project, PROJECTS[project]))
    s = m._engine
    n = int(s.update_segments(s.params_values()))
    get = {"sigma": s.calculated_segments, "e1": s.calculated_segments_e1, "e2": s.calculated_segments_e2,
           "energy": s.calculated_energies, "angle": s.calculated_angles, "angular": s.calculated_angular_dists}
    res = []
    for i in range(n):
        res.append({k: ["nan" if math.isnan(v) else v for v in map(float, f(i))] for k, f in get.items()})
    return res


def compare_api(project, plain, poisoned, problems, known):
    if len(plain) != len(poisoned):
        problems.append(f"{project} API: {len(plain)} and {len(poisoned)} segments")
        return
    for i, (a, b) in enumerate(zip(plain, poisoned)):
        for k in a:
            if len(a[k]) != len(b[k]):
                problems.append(f"{project} API segment {i} {k}: {len(a[k])} and {len(b[k])} values")
                continue
            for j, (x, y) in enumerate(zip(a[k], b[k])):
                if x == y:
                    continue
                if k in ("e1", "e2") and x == 0.0 and y == "nan":
                    known.add(f"K1 {project} segment {i}")
                    continue
                problems.append(f"{project} API segment {i} {k}[{j}]: {x} without, {y} with poisoning")
                break


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--api":
        print(json.dumps(api_values(sys.argv[2])))
        return 0
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--binary")
    ap.add_argument("projects", nargs="*", default=list(PROJECTS))
    args = ap.parse_args()
    binary = up.find_binary(args.binary)
    problems, known = [], set()
    with tempfile.TemporaryDirectory() as tmp:
        runs = [(p, name, menu, []) for p in args.projects for name, menu in MENUS.items()
                if not (name == "extrap" and p == "T2")]  # T2 has no <segmentsTest>
        if FIT[0] in args.projects:
            runs.append((FIT[0], "fit", FIT[1], FIT[2]))
        for project, name, menu, extra in runs:
            rc0, plain = run_cli(binary, project, menu, extra, False, tmp)
            rc1, poisoned = run_cli(binary, project, menu, extra, True, tmp)
            if rc0 != rc1:
                problems.append(f"{project} {name}: exit {rc0} without, {rc1} with poisoning")
            if not plain:
                problems.append(f"{project} {name}: no output")
            for f in sorted(set(plain) | set(poisoned)):
                if plain.get(f) != poisoned.get(f):
                    problems.append(f"{project} {name}: output/{f} differs with poisoning")
            print(f"{project} {name}: {len(plain)} files compared")
        rc0, plain = run_rate(binary, False, tmp)
        rc1, poisoned = run_rate(binary, True, tmp)
        if rc0 != rc1 or rc0 != 0:
            problems.append(f"rate: exit {rc0} without, {rc1} with poisoning")
        if "reactionrates.out" not in plain:
            problems.append("rate: no reactionrates.out")
        for f in sorted(set(plain) | set(poisoned)):
            if plain.get(f) != poisoned.get(f):
                problems.append(f"rate: output/{f} differs with poisoning")
        print(f"rate p12C: {len(plain)} files compared")
        made = os.path.join(tmp, "transmission_yield")
        os.makedirs(made)
        mk = subprocess.run([sys.executable, "-I", os.path.join(TESTS, "transmission_yield", "make_models.py"),
                             os.path.join(TESTS, "reaction_rate", "models", "n12C.azr"),
                             os.path.join(TESTS, "natural_neon", "natural_neon.azr"), made],
                            capture_output=True, text=True)
        if mk.returncode != 0:
            problems.append(f"transmission_yield models: {mk.stderr.strip()[-300:]}")
        shutil.copytree(os.path.join(TESTS, "natural_neon", "data"), os.path.join(made, "data"), dirs_exist_ok=True)
        for name, (rel, azr) in MODEL_PROJECTS.items():
            src = made if rel == "transmission_yield" else os.path.join(TESTS, rel)
            rc0, plain = run_cli(binary, name, MENUS["calc"], [], False, tmp, src, azr)
            rc1, poisoned = run_cli(binary, name, MENUS["calc"], [], True, tmp, src, azr)
            if rc0 != rc1:
                problems.append(f"{name} calc: exit {rc0} without, {rc1} with poisoning")
            if "chiSquared.out" not in plain:
                problems.append(f"{name} calc: no chiSquared.out")
            for f in sorted(set(plain) | set(poisoned)):
                if plain.get(f) != poisoned.get(f):
                    problems.append(f"{name} calc: output/{f} differs with poisoning")
            print(f"{name} calc: {len(plain)} files compared")
    try:
        import pyazr.session  # noqa: F401
        have_api = True
    except Exception as e:  # the module is optional (pyazr tests skip without it)
        print(f"pyazr part skipped: {e}")
        have_api = False
    if have_api:
        for project in args.projects:
            vals = []
            for poison in ("0", "1"):
                env = dict(os.environ, AZURE_POISON_RESULTS=poison)
                r = subprocess.run([sys.executable, os.path.abspath(__file__), "--api", project], env=env,
                                   capture_output=True, text=True, timeout=1800)
                if r.returncode != 0:
                    problems.append(f"{project} API run failed: {r.stderr.strip()[-500:]}")
                    break
                vals.append(json.loads(r.stdout.strip().splitlines()[-1]))
            if len(vals) == 2:
                compare_api(project, vals[0], vals[1], problems, known)
                print(f"{project} API: {len(vals[0])} segments compared")
    if known:
        print("known: " + ", ".join(sorted(known)))
    for p in problems:
        print(f"FAIL: {p}")
    print(f"poisoned results: {'FAIL' if problems else 'ok'} ({len(problems)} problems, {len(known)} known)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
