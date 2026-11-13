#!/usr/bin/env python3
"""The identity of tests/shift_true_energy: every data set with a shift against its twin without one at the shifted
energies (README.md).

    python3 -I compare.py TOLERANCE MODEL OUTPUT_DIR [MODEL OUTPUT_DIR ...]

Each MODEL (models/identity.azr, or models/shifted.azr and models/twins.azr) gives the data sets' order and names, its
OUTPUT_DIR is its menu-1 output, whose file holds the data sets in the model's order, one line per point (c.m. energy,
..., the model value in column 4).  Prints, per case, the largest relative difference of the c.m. energies and of the
model values; exit 1 when a model value is off by more than TOLERANCE or a pair is incomplete.
"""
import glob
import re
import sys


def main():
    tol = float(sys.argv[1])
    names, block = [], {}
    for model, out in zip(sys.argv[2::2], sys.argv[3::2]):
        these = re.findall(r"^dataset (\S+)", open(model).read(), re.M)
        files = glob.glob(out + "/AZUREOut_aa=1_*.out")
        if len(files) != 1:
            print("%s: expected one output file, found %d" % (out, len(files)))
            return 1
        rows = [l.split() for l in open(files[0]) if len(l.split()) >= 5 and not l.startswith("#")]
        npts = len(rows) // len(these) if these else 0
        if npts * len(these) != len(rows) or npts == 0:
            print("%s: %d rows for %d data sets" % (out, len(rows), len(these)))
            return 1
        names += these
        block.update({n: rows[i * npts:(i + 1) * npts] for i, n in enumerate(these)})
    worst, bad = 0.0, 0
    for n in names:
        if n.endswith("_twin"):
            continue
        twin = block.get(n + "_twin")
        if twin is None:
            print("%s: no twin" % n)
            bad += 1
            continue
        de = max(abs(float(a[0]) - float(b[0])) / abs(float(b[0])) for a, b in zip(block[n], twin))
        dv = max(abs(float(a[3]) - float(b[3])) / max(abs(float(b[3])), 1e-300) for a, b in zip(block[n], twin))
        good = dv <= tol and de <= 1e-12  # NaN fails
        worst = max(worst, dv) if dv == dv else float("nan")
        flag = "ok  " if good else "FAIL"
        bad += not good
        print("%s %-40s energies %.1e  values %.2e  (%s)" % (flag, n, de, dv, " ".join("%.6e" % float(a[3]) for a in block[n])))
    print("worst relative difference %.3e over %d pairs, tolerance %.1e" % (worst, len(names) // 2, tol))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
