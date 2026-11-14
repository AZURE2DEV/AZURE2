#!/usr/bin/env python3
"""Golden outputs of the test projects T1 and T3: write them, or compare a build against them.

    python3 tests/upgrade/golden/golden.py [--binary AZURE2] [--write [--engine-commit HASH]] [--rtol 1e-6] [T1 T3]

The binary is --binary, else $AZURE2_BIN / $AZURE2_BINARY, else the newest build*/src/AZURE2 of the checkout.
Each project is copied to a temporary directory and run three ways, with OMP_NUM_THREADS=2 and default flags
except --no-coulomb-table (exact Coulomb functions whether or not a table exists on the machine):

  calc    menu 1, calculate segments from data at the truth (T<n>.azr): chiSquared.out, AZUREOut_*.out,
          parameters.out, param.par, scatteringLengths.out
  extrap  menu 3, the extrapolations: AZUREOut_*.extrap
  rate    menu 5, the reaction rate of the pair in PROJECTS below over temperatures.dat: reactionrates.out

The projects are format 3; the golden files were written from their old originals (tests/convert/fixtures).  A
golden file is compared as the format-3 project writes it (golden_lines): param.par with the converter's parameter
names, and the classes tools/format3_gate.py lists for T3 (the scattering lengths as model data).

--write stores the files under golden/<project>/<mode>/ and the engine commit in golden/ENGINE.  Without it the
files are compared token by token: every non-numeric token equal, every number within

    |a - b| <= rtol * max(|a|, |b|)      (rtol = 1e-6 by default; both below 1e-200 count as equal)

which is about ten times the last printed digit of AZUREOut (11 significant digits) and above thread-count and
summation-order noise on one machine.  Across compilers (GCC, Clang, MinGW) differences near 1e-5 are expected
(root CLAUDE.md, "Determinism"); pass --rtol 2e-5 there.  Exit status 1 when anything differs.
"""
import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import upgrade_project as up  # noqa: E402
import format3_gate  # noqa: E402  (the accepted differences between an old file's outputs and its conversion's)

# Exact Coulomb functions: the universal table (~/.azure2) changes results at the 1e-6 level and may or may not
# exist on a machine, so the golden runs never use it.
FLAGS = ["--no-coulomb-table"]
PROJECTS = {"T1": (1, 2), "T3": (2, 1)}   # rate entrance, exit (the old originals' <lastRun>)
MODES = {
    "calc": (up.CALC, ["chiSquared.out", "AZUREOut_*.out", "parameters.out", "param.par",
                       "scatteringLengths.out"]),
    "extrap": (up.EXTRAP, ["AZUREOut_*.extrap"]),
    "rate": (None, ["reactionrates.out"]),
}


def run_project(binary, name, dest):
    src = os.path.join(HERE, "..", name)
    for mode, (command, patterns) in MODES.items():
        work = tempfile.mkdtemp(prefix=f"{name}_{mode}_")
        try:
            shutil.copy(os.path.join(src, f"{name}.azr"), work)
            shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
            shutil.copy(os.path.join(src, "temperatures.dat"), work)
            if command is None:
                command = up.rate_command(*PROJECTS[name])
            up.run(binary, work, f"{name}.azr", command, flags=FLAGS, threads=2)
            out = os.path.join(dest, name, mode)
            os.makedirs(out, exist_ok=True)
            for pat in patterns:
                files = glob.glob(os.path.join(work, "output", pat))
                if not files:
                    raise RuntimeError(f"{name} {mode}: no output/{pat}")
                for f in files:
                    shutil.copy(f, out)
        finally:
            shutil.rmtree(work, ignore_errors=True)


def parameter_names(binary, name):
    """{old parameter name: format-3 name}: the report of the converter run on the project's old original
    (tests/convert/fixtures/tests/upgrade/<name>/<name>.azr), in a copy beside the project's data files."""
    work = tempfile.mkdtemp(prefix=f"{name}_names_")
    try:
        src = os.path.join(HERE, "..", name)
        shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
        shutil.copy(os.path.join(ROOT, "tests", "convert", "fixtures", "tests", "upgrade", name, f"{name}.azr"),
                    os.path.join(work, "old.azr"))
        subprocess.run([binary, "--convert", "old.azr", "new.azr", "--to", "3", "--exact-grids", "--json-report",
                        "report.json"], cwd=work, capture_output=True, check=True,
                       env=dict(os.environ, OMP_NUM_THREADS="1", LC_ALL="C"))
        with open(os.path.join(work, "report.json")) as f:
            return json.load(f)["parameters"]
    finally:
        shutil.rmtree(work, ignore_errors=True)


def golden_lines(name, fname, path, names):
    """The golden file as the format-3 project writes it.  The golden files were written from the old project; the
    format-3 conversion (tests/upgrade/<name>/<name>.azr) writes the same numbers, with the parameter names of the
    converter's map in param.par and, where tools/format3_gate.py lists them for this project, the classes
    'scattering lengths as model data (label)' (chiSquared.out) and '(table)' (scatteringLengths.out: the old
    section's table of measured values is in modelData.out)."""
    with open(path) as f:
        lines = f.read().splitlines()
    project = f"upgrade/{name}/{name}.azr"
    if fname == "param.par":
        return [" ".join([names.get(t[0], t[0])] + t[1:]) if t else "" for t in (l.split() for l in lines)]
    if fname == "chiSquared.out" and format3_gate.listed("scattering lengths as model data (label)", fname, project):
        return [l.replace("Scattering-Length", "Model-Data") for l in lines]
    if fname == "scatteringLengths.out" and format3_gate.listed("scattering lengths as model data (table)", fname,
                                                                 project):
        cut = [i for i, l in enumerate(lines) if l.startswith("Constraints (<scatteringLengths>)")]
        return lines[:cut[0]] if cut else lines
    return lines


def tokens(lines):
    return [line.split() for line in lines]


def compare_file(a, golden, rtol):
    """Return (worst relative difference, list of problems) for an output file and the golden file's lines."""
    with open(a) as f:
        ta = tokens(f.read().splitlines())
    tb = tokens(golden)
    problems, worst = [], 0.0
    if len(ta) != len(tb):
        return float("inf"), [f"{len(ta)} lines, golden {len(tb)}"]
    for n, (la, lb) in enumerate(zip(ta, tb), 1):
        if len(la) != len(lb):
            problems.append(f"line {n}: {len(la)} fields, golden {len(lb)}")
            continue
        for x, y in zip(la, lb):
            try:
                fx, fy = float(x.rstrip(",")), float(y.rstrip(","))
            except ValueError:
                if x != y:
                    problems.append(f"line {n}: '{x}' vs golden '{y}'")
                continue
            scale = max(abs(fx), abs(fy))
            if scale < 1e-200:
                continue
            d = abs(fx - fy) / scale
            worst = max(worst, d)
            if d > rtol:
                problems.append(f"line {n}: {x} vs golden {y} (rel {d:.2e})")
    return worst, problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("projects", nargs="*", default=list(PROJECTS))
    ap.add_argument("--binary")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--engine-commit", help="commit the binary was built from (default: this checkout's HEAD)")
    ap.add_argument("--rtol", type=float, default=1e-6)
    args = ap.parse_args()
    binary = up.find_binary(args.binary)

    if args.write:
        commit = args.engine_commit or subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True,
                                                      text=True).stdout.strip()
        for name in args.projects:
            shutil.rmtree(os.path.join(HERE, name), ignore_errors=True)
            run_project(binary, name, HERE)
        with open(os.path.join(HERE, "ENGINE"), "w") as f:
            f.write(f"engine commit {commit}\nprojects {' '.join(args.projects)}\n"
                    f"run: OMP_NUM_THREADS=2, flags {' '.join(FLAGS)}, modes {' '.join(MODES)}\n")
        print(f"golden outputs written for {' '.join(args.projects)} (engine {commit})")
        return 0

    failures = 0
    tmp = tempfile.mkdtemp(prefix="golden_cmp_")
    try:
        for name in args.projects:
            names = parameter_names(binary, name)
            run_project(binary, name, tmp)
            for mode in MODES:
                gdir = os.path.join(HERE, name, mode)
                golden = sorted(os.listdir(gdir))
                got = sorted(os.listdir(os.path.join(tmp, name, mode)))
                if golden != got:
                    print(f"  FAIL  {name} {mode}: files {got}, golden {golden}")
                    failures += 1
                for fname in golden:
                    if fname not in got:
                        continue
                    worst, problems = compare_file(os.path.join(tmp, name, mode, fname),
                                                   golden_lines(name, fname, os.path.join(gdir, fname), names),
                                                   args.rtol)
                    if problems:
                        failures += 1
                        print(f"  FAIL  {name} {mode} {fname}: {len(problems)} differences, e.g. {problems[0]}")
                    else:
                        print(f"  ok    {name} {mode} {fname} (largest relative difference {worst:.1e})")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"golden comparison: {'FAIL' if failures else 'ok'} ({failures} problems, rtol {args.rtol:g})")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
