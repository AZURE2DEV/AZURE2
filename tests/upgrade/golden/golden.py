#!/usr/bin/env python3
"""Golden outputs of the test projects T1 and T3: write them, or compare a build against them.

    python3 tests/upgrade/golden/golden.py [--binary AZURE2] [--write [--engine-commit HASH]] [--rtol 1e-6] [T1 T3]

The binary is --binary, else $AZURE2_BIN / $AZURE2_BINARY, else the newest build*/src/AZURE2 of the checkout.
Each project is copied to a temporary directory and run three ways, with OMP_NUM_THREADS=2 and default flags
except --no-coulomb-table (exact Coulomb functions whether or not a table exists on the machine):

  calc    menu 1, calculate segments from data at the truth (T<n>.azr): chiSquared.out, AZUREOut_*.out,
          parameters.out, param.par, scatteringLengths.out
  extrap  menu 3, the <segmentsTest> lines: AZUREOut_*.extrap
  rate    menu 5, the reaction rate of the pair named in <lastRun> over temperatures.dat: reactionrates.out

--write stores the files under golden/<project>/<mode>/ and the engine commit in golden/ENGINE.  Without it the
files are compared token by token: every non-numeric token equal, every number within

    |a - b| <= rtol * max(|a|, |b|)      (rtol = 1e-6 by default; both below 1e-200 count as equal)

which is about ten times the last printed digit of AZUREOut (11 significant digits) and above thread-count and
summation-order noise on one machine.  Across compilers (GCC, Clang, MinGW) differences near 1e-5 are expected
(root CLAUDE.md, "Determinism"); pass --rtol 2e-5 there.  Exit status 1 when anything differs.
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import upgrade_project as up  # noqa: E402

# Exact Coulomb functions: the universal table (~/.azure2) changes results at the 1e-6 level and may or may not
# exist on a machine, so the golden runs never use it.
FLAGS = ["--no-coulomb-table"]
PROJECTS = {"T1": (1, 2), "T3": (2, 1)}   # rate entrance, exit (also in each project's <lastRun>)
MODES = {
    "calc": (up.CALC, ["chiSquared.out", "AZUREOut_*.out", "parameters.out", "param.par",
                       "scatteringLengths.out"]),
    "extrap": (up.EXTRAP, ["AZUREOut_*.extrap"]),
    "rate": (None, ["reactionrates.out"]),
}


def run_project(binary, name, dest):
    src = os.path.join(HERE, "..", name)
    for mode, (answers, patterns) in MODES.items():
        work = tempfile.mkdtemp(prefix=f"{name}_{mode}_")
        try:
            shutil.copy(os.path.join(src, f"{name}.azr"), work)
            shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
            shutil.copy(os.path.join(src, "temperatures.dat"), work)
            if answers is None:
                answers = up.rate_answers(*PROJECTS[name])
            up.run(binary, work, f"{name}.azr", answers, flags=FLAGS, threads=2)
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


def tokens(path):
    with open(path) as f:
        return [line.split() for line in f]


def compare_file(a, b, rtol):
    """Return (worst relative difference, list of problems) for two output files."""
    ta, tb = tokens(a), tokens(b)
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
                    worst, problems = compare_file(os.path.join(tmp, name, mode, fname), os.path.join(gdir, fname),
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
