#!/usr/bin/env python3
"""A natural target: several compound nuclei in one project, their cross sections summed with the abundances.

A project may hold the levels of several compound nuclei (here p + 12C -> 13N and p + 13C -> 14N).  A ``target``
record lists the isotopes of a sample with their abundances (``isotope 12C 0.6 fixed``); a data set or extrapolation
measured on it names the target, and its reaction may be a sum (``12C(p,p)12C + 13C(p,p)13C``): the engine weights
each nucleus's cross section with its abundance.  An abundance written ``fit`` is a parameter of the fit
(``<target>.isotope.<nucleus>``).

The script reads the target's isotopes from the project (``Project.records()``), evaluates the extrapolations in a
session made for them, and checks the sum against the single-nucleus curves (pure arithmetic of the engine's
numbers: the engine computes every curve).

Usage (from the repository's root):

    python3 pyazr/examples/natural_target.py tests/natural_target/models/M.azr
"""
import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np  # noqa: E402

from pyazr import Project, Session  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="a format-3 project with a target of several isotopes and extrapolations")
    args = ap.parse_args()

    p = Project.open(args.project)
    records = p.records()
    for name, r in records.items():
        if r.record_kind == "target" and r.get("isotope"):
            print(f"target {name}: isotopes {r['isotope']}")
    print(f"pairs: {[(q.key, q.label) for q in p.pairs()]}")

    with Session(args.project, extrapolations=True) as s:
        curves = {d.name: d for d in s.evaluate().data}
    for name, d in curves.items():
        r = records[name]
        on = f" on the target {r['target']}" if r.get("target") else ""
        unit = r["unit"]                    # the record's (evaluate(x).data[i].unit is None for an extrapolation)
        print(f"  {name}: {r['reaction']}{on}, {r['quantity']} in {unit}, {d.n} points, "
              f"E_cm {d.energies.min():.3f}-{d.energies.max():.3f} MeV, largest {d.model.max():.5g} {unit}")

    # the abundance-weighted sum, from the engine's single-nucleus curves (the grids are the same lab energies)
    if {"x1", "x2", "x3", "x4"} <= set(curves):
        x1, x2, x3, x4 = (curves[k].model for k in ("x1", "x2", "x3", "x4"))
        print(f"\nlargest |x3 - (0.6 x1 + 0.4 x2)| / x3 = {np.max(np.abs(x3 - (0.6 * x1 + 0.4 * x2)) / x3):.2e}")
        print(f"largest |x4 - 0.4 x2| / x4 = {np.max(np.abs(x4 - 0.4 * x2) / x4):.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
