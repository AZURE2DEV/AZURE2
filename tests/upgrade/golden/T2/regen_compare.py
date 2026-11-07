#!/usr/bin/env python3
"""Golden outputs of T2: regenerate them with an AZURE2 binary and compare with the committed set.

Three runs of tests/upgrade/T2, each from a clean copy of the project:
  calc_start  menu 1 on T2.azr       (the fit project at its start values: every observable and
                                       apparatus effect of T2 evaluated once, chi2 per segment)
  calc_truth  menu 1 on T2_true.azr  (the truth)
  fit         menu 2 on T2.azr with --use-lm (falls back to MIGRAD with the numerical gradient: the
                                       THM segments have no analytic Jacobian; see check_recovery.py)
All runs use --no-coulomb-table (exact Coulomb functions), so the result does not depend on whether a
Coulomb table exists on the machine, and OMP_NUM_THREADS from --threads (results must not depend on it).

Comparison, token by token (the files must have the same tokens; every number is compared):
  calc_*:  |a - b| <= 1e-5 max(|a|, |b|) + 1e-300   (the cross-compiler spread stated in CLAUDE.md)
  fit:     the same with 1e-3: a MIGRAD path ends where its stopping rule (EDM) says, and a change of
           rounding moves the end point by a fraction of the parameter error, not by 1e-5.
A difference outside these is a change of the engine's numbers, to be explained before re-baselining
(tests/CLAUDE.md, "Re-baselining").

    AZURE2_BINARY=/path/to/AZURE2 python3 tests/upgrade/golden/T2/regen_compare.py            # compare
    AZURE2_BINARY=/path/to/AZURE2 python3 tests/upgrade/golden/T2/regen_compare.py --regenerate --commit <hash>

--regenerate overwrites the golden files and ENGINE (the engine commit, given with --commit, and the
date).  Standard library only.  About 3 s for the calculations and 6 min for the fit at 8 threads.
"""
import argparse
import datetime
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
T2 = os.path.normpath(os.path.join(HERE, "..", "..", "T2"))
PROJECT = ["T2.azr", "T2_true.azr", "cov_diff.txt", "prior_cov.txt", "data"]
RUNS = {
    "calc_start": ("T2.azr", "1\n\n\n7\n", 1e-5, []),
    "calc_truth": ("T2_true.azr", "1\n\n\n7\n", 1e-5, []),
    "fit": ("T2.azr", "2\nn\n\n\n7\n", 1e-3, ["--use-lm"]),
}
KEEP = re.compile(r"^(chiSquared\.out|AZUREOut_.*\.out|thm_experiments\.out|parameters\.out|param\.sav|"
                  r"normalizations\.out|shifts\.out)$")
NUM = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


def tokens(path):
    return open(path).read().replace(",", " ").split()


def compare(a, b, rtol):
    ta, tb = tokens(a), tokens(b)
    if len(ta) != len(tb):
        return f"{len(ta)} tokens, golden {len(tb)}"
    worst, where = 0.0, None
    for i, (x, y) in enumerate(zip(ta, tb)):
        if NUM.match(x) and NUM.match(y):
            fx, fy = float(x), float(y)
            d = abs(fx - fy)
            scale = max(abs(fx), abs(fy))
            r = d / scale if scale > 0 else 0.0
            if d > rtol * scale + 1e-300 and r > worst:
                worst, where = r, (i, x, y)
        elif x != y:
            return f"token {i}: '{x}', golden '{y}'"
    if where:
        return f"token {where[0]}: {where[1]}, golden {where[2]} (rel {worst:.3g} > {rtol:g})"
    return None


def run(binary, threads, name, work):
    azr, menu, _, extra = RUNS[name]
    d = os.path.join(work, name)
    os.makedirs(d)
    for f in PROJECT:
        src = os.path.join(T2, f)
        (shutil.copytree if os.path.isdir(src) else shutil.copy)(src, os.path.join(d, f))
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(d, sub))
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), LC_ALL="C")
    r = subprocess.run([binary, "--no-gui", "--no-readline", "--no-coulomb-table", *extra, azr], input=menu, cwd=d,
                       env=env, capture_output=True, text=True, timeout=3600)
    if not os.path.isfile(os.path.join(d, "output", "chiSquared.out")):
        sys.exit(f"{name}: no chiSquared.out\n{(r.stdout + r.stderr)[-3000:]}")
    return os.path.join(d, "output")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--binary", default=os.environ.get("AZURE2_BINARY", ""))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--regenerate", action="store_true")
    ap.add_argument("--commit", default="", help="engine commit recorded in ENGINE (with --regenerate)")
    ap.add_argument("--only", choices=sorted(RUNS), action="append", help="run only these (repeatable)")
    args = ap.parse_args()
    if not args.binary or not os.access(args.binary, os.X_OK):
        sys.exit("set AZURE2_BINARY or pass --binary path/to/AZURE2")
    binary = os.path.abspath(args.binary)
    work = tempfile.mkdtemp(prefix="t2gold.")
    fail = 0
    try:
        for name in args.only or list(RUNS):
            out = run(binary, args.threads, name, work)
            produced = sorted(f for f in os.listdir(out) if KEEP.match(f))
            gdir = os.path.join(HERE, name)
            if args.regenerate:
                shutil.rmtree(gdir, ignore_errors=True)
                os.makedirs(gdir)
                for f in produced:
                    shutil.copy(os.path.join(out, f), gdir)
                print(f"{name}: wrote {len(produced)} files")
                continue
            golden = sorted(os.listdir(gdir)) if os.path.isdir(gdir) else []
            if produced != golden:
                print(f"FAIL {name}: files {produced}, golden {golden}")
                fail += 1
                continue
            for f in produced:
                msg = compare(os.path.join(out, f), os.path.join(gdir, f), RUNS[name][2])
                print(f"  {'FAIL' if msg else 'ok  '} {name}/{f}" + (f": {msg}" if msg else ""))
                fail += 1 if msg else 0
        if args.regenerate:
            with open(os.path.join(HERE, "ENGINE"), "w") as fh:
                fh.write(f"engine commit: {args.commit or 'unknown'}\n")
                fh.write(f"generated: {datetime.date.today().isoformat()} with --no-coulomb-table, "
                         f"OMP_NUM_THREADS={args.threads}\n")
                fh.write("runs: " + ", ".join(args.only or list(RUNS)) + " (regen_compare.py)\n")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if not args.regenerate:
        print("PASS" if fail == 0 else f"FAIL ({fail})")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
