#!/usr/bin/env python3
"""Compare the capture cross sections of engine runs with the reference (make_reference.py) or with each other.

    compare.py reference.csv RUN_DIR [RUN_DIR ...] --tol T            every run = the reference within T (relative)
    compare.py --spread RUN_DIR RUN_DIR [...] --tol T                  the runs agree with each other within T
    compare.py --spread RUN_DIR RUN_DIR [...] --at-least T             the runs differ somewhere by at least T

RUN_DIR holds output/AZUREOut_aa=1_R=2.extrap (3/2- ground state) and output/AZUREOut_aa=1_R=3.extrap (1/2-), whose
columns are E_cm (MeV), E_x, angle, sigma (b), S (MeV b).  Prints the largest relative difference and where; exits 1
when the condition fails or a file is missing.
"""
import csv
import os
import sys


def read_run(run):
    out = {}
    for state, R in (("32", 2), ("12", 3)):
        path = os.path.join(run, "output", f"AZUREOut_aa=1_R={R}.extrap")
        with open(path) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 5:
                    out[(state, round(float(parts[0]), 6))] = float(parts[3])
    for E in sorted({e for (_, e) in out}):
        if ("32", E) in out and ("12", E) in out:
            out[("sum", E)] = out[("32", E)] + out[("12", E)]
    return out


def read_reference(path):
    out = {}
    with open(path) as f:
        rows = [r for r in f if not r.startswith("#")]
    for r in csv.DictReader(rows):
        E = round(float(r["E_cm_MeV"]), 6)
        out[("32", E)] = float(r["sigma_32_b"])
        out[("12", E)] = float(r["sigma_12_b"])
        out[("sum", E)] = float(r["sigma_b"])
    return out


def worst(a, b):
    keys = sorted(set(a) & set(b))
    if not keys:
        raise SystemExit("no common points")
    w, where = -1.0, None
    for k in keys:
        d = abs(a[k] / b[k] - 1.0)
        if d > w:
            w, where = d, k
    return w, where, len(keys)


def main():
    args = sys.argv[1:]
    tol = at_least = None
    if "--tol" in args:
        i = args.index("--tol")
        tol = float(args[i + 1])
        del args[i:i + 2]
    if "--at-least" in args:
        i = args.index("--at-least")
        at_least = float(args[i + 1])
        del args[i:i + 2]
    try:
        if args[0] == "--spread":
            runs = [read_run(r) for r in args[1:]]
            w, where, n = -1.0, None, 0
            for r in runs[1:]:
                d, k, n = worst(r, runs[0])
                if d > w:
                    w, where = d, k
        else:
            ref = read_reference(args[0])
            w, where, n = -1.0, None, 0
            for r in args[1:]:
                d, k, n = worst(read_run(r), ref)
                if d > w:
                    w, where = d, k
    except (OSError, ValueError, KeyError, SystemExit) as e:
        print(f"error: {e}")
        return 1
    print(f"largest relative difference {w:.3e} ({where[0]}, E = {where[1]} MeV; {n} values per run)")
    if tol is not None and not w <= tol:
        return 1
    if at_least is not None and not w >= at_least:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
