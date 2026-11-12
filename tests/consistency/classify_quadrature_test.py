#!/usr/bin/env python3
"""tools/classify_quadrature.py: the per-file rule of docs/upgrade/energy_grid.md section 7.1 (largest relative distance
to the fine reference <= 1.2 x today's, and rms <= today's) and the fine reference taken per file from a list of tags.

Four runs laid out as tools/runall.sh writes them, in a temporary directory: ref (today's engine), fine8 and fine4
(refined runs; fine8 lacks one file), new (the engine under test).  Ten folded numbers per file, ref 1e-3 from the fine
value at one number and 4e-4 at the others (rms 4.9e-4):
  better.out  new half as far everywhere                                   -> quadrature
  spread.out  new as far as ref at the largest, 9e-4 at the others          -> failure (rms 9.1e-4), passes the maximum
  peak.out    new 1.3e-3 at the largest, 5e-5 at the others                 -> failure (maximum), passes the rms (4.1e-4)
  only4.out   as better.out, present in fine4 alone                         -> quadrature, against fine4
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "..", "tools", "classify_quadrature.py")

FINE = [1.0 + 0.1 * i for i in range(10)]


def values(first, others):
    """Ten numbers at relative distances `first` (number 0) and `others` above the fine values."""
    return [f * (1.0 + (first if i == 0 else others)) for i, f in enumerate(FINE)]


def write(cmp, tag, name, numbers):
    out = os.path.join(cmp, "brune", tag, ".res-t", "t.azr", "output")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, name), "w") as f:
        f.write("".join(f"{i} {v:.12e}\n" for i, v in enumerate(numbers)))


def main():
    ref = values(1e-3, 4e-4)
    new = {"better.out": values(5e-4, 2e-4), "spread.out": values(1e-3, 9e-4), "peak.out": values(1.3e-3, 5e-5),
           "only4.out": values(5e-4, 2e-4)}
    expected = {"better.out": ("quadrature", "fine8"), "spread.out": ("failure", "fine8"),
                "peak.out": ("failure", "fine8"), "only4.out": ("quadrature", "fine4")}
    failures = []
    with tempfile.TemporaryDirectory() as cmp:
        for name, numbers in new.items():
            write(cmp, "ref", name, ref)
            write(cmp, "new", name, numbers)
            write(cmp, "fine4", name, FINE)
            if name != "only4.out":
                write(cmp, "fine8", name, FINE)
        r = subprocess.run([sys.executable, TOOL, cmp, "ref", "fine8,fine4", "new", "brune"], capture_output=True,
                           text=True)
        print(r.stdout, end="")
        if r.returncode != 1:
            failures.append(f"exit {r.returncode} (expected 1: two files fail): {r.stderr}")
        lines = {os.path.basename(l.split("\t")[1]): l.split("\t") for l in r.stdout.splitlines() if not l.startswith("##")}
        for name, (cls, tag) in expected.items():
            row = lines.get(name)
            if row is None:
                failures.append(f"{name}: not listed")
                continue
            if row[8] != f"file: {cls}":
                failures.append(f"{name}: {row[8]}, expected file: {cls}")
            if row[-1] != tag:
                failures.append(f"{name}: compared with {row[-1]}, expected {tag}")
        if "## brune: 4 files, 0 identical, 2 quadrature, 2 failure" not in r.stdout:
            failures.append("summary: expected 4 files, 2 quadrature, 2 failure")
    for f in failures:
        print("FAIL:", f)
    print("classify_quadrature_test:", "FAILED" if failures else "passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
