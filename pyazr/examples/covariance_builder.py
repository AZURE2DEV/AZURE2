#!/usr/bin/env python3
"""A data set's covariance file built from named parts (pyazr.covariance), checked as the engine checks it.

The parts here: the points' own errors (the error column), a 5 % scale uncertainty fully correlated over the data
set (on the model: a relvector term), and 3 % correlated over 0.2 MeV of lab energy.  The script prints what the
engine will make of it (the checks, the notes, the total error per point beside the file's), writes the file into
the project's directory with the parts on its #pyazr.covariance line, names it on the data set and saves the project
as a copy (the original file is not changed).  The chi-squared with and without the covariance is printed.

Usage (from the repository's root, on a copy of a project whose data files are beside it):

    python3 pyazr/examples/covariance_builder.py project.azr dataset_name [--out project_cov.azr]
"""

import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np  # noqa: E402

from pyazr import Project, Session, covariance  # noqa: E402

PARTS = [{"kind": "own"},
         {"kind": "scale", "percent": 5},
         {"kind": "length", "percent": 3, "length": 0.2}]


def chi2_of(path, name):
    with Session(path) as s:
        ev = s.evaluate(np.asarray(s.parameters.x0, float))
        return next(d.chi2 for d in ev.data if d.name == name)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr")
    ap.add_argument("dataset")
    ap.add_argument("--out", default=None, help="the project copy that names the file (default <name>_cov.azr)")
    args = ap.parse_args()

    p = Project.open(args.azr)
    c = covariance.build(p, [args.dataset], PARTS)
    print(f"{c.rows.n} rows ({c.fitted.size} fitted) in {c.rows.unit}; refusal: {c.refusal}")
    print("checks:", {k: v for k, v in c.checks.items()})
    for note in c.notes:
        print(" -", note)
    print(f"{'row':>4} {'E_lab':>8} {'value':>12} {'column':>12} {'sqrt(C_ii)':>12}")
    for i in c.fitted[:5]:
        print(f"{c.rows.record[i]:4d} {c.rows.energy[i]:8.4f} {c.rows.value[i]:12.5g} {c.rows.error[i]:12.5g} "
              f"{c.total_error[i]:12.5g}")
    name = f"cov_{args.dataset}.txt"
    c.write(os.path.join(p.directory, name))
    out = args.out or os.path.splitext(args.azr)[0] + "_cov.azr"
    p.set(args.dataset, covariance=name)
    p.save(out)
    back = covariance.read_parts(os.path.join(p.directory, name))
    print(f"written: {name} (parts read back: {len(back['parts'])}); project: {out}")
    print(f"chi2 of {args.dataset}: {chi2_of(args.azr, args.dataset):.10g} without, "
          f"{chi2_of(out, args.dataset):.10g} with the covariance")
    return 0


if __name__ == "__main__":
    sys.exit(main())
