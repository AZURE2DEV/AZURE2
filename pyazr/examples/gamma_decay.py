#!/usr/bin/env python3
"""A level's radiative decay written as a ``gamma_decay`` record (total radiative width, branchings, mixing ratios),
fitted to capture data to each final state.

Instead of one partial radiative width per final state in ``<levels>``, a ``gamma_decay`` record gives a level's
total radiative width Gamma_gamma (eV), its branching b_f to each final state (``g+15O[6.1763]`` is the photon pair to
the state at 6.1763 MeV) and one final state as the ``remainder`` (1 - the sum of the others); a mixing ratio
``mixing`` where two multipolarities feed one final state.  Each number written ``fit`` is a fit parameter of kind
``radiative`` whose owner fields name the record (``decay_record``, ``decay_part``, ``final_state``).

The script prints the record, its rows of the parameter table, fits with the engine and prints Gamma_gamma, the
branchings and the remainder with their errors, and the chi-squared of the data set of each final state.

Usage (from the repository's root):

    python3 pyazr/examples/gamma_decay.py tests/gamma_decay/models/record.azr [--record level5]
"""
import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np  # noqa: E402

from pyazr import Project, Session  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="a format-3 project with a gamma_decay record")
    ap.add_argument("--record", default=None, help="the gamma_decay record's name (default: the first)")
    args = ap.parse_args()

    decays = {n: r for n, r in Project.open(args.project).records().items() if r.record_kind == "gamma_decay"}
    name = args.record or next(iter(decays))
    print(f"gamma_decay {name}:")
    for key, value in decays[name].items():
        for v in (value if isinstance(value, list) else [value]):
            print(f"    {key:<12} {v}")

    with Session(args.project) as s:
        rows = [q for q in s.parameters.free if getattr(q, "decay_record", None) == name]
        start = s.evaluate()
        r = s.fit_with_engine()
        end = s.evaluate(r.x)
    print(f"\nfit: {r.status}; objective {start.objective:.4f} -> {r.objective:.4f}")
    branching = 0.0
    for q in rows:
        k = r.names.index(q.name)
        unit = "eV" if q.name.endswith("total_width") else "(fraction)"
        print(f"  {q.name:<32} {r.x[k]:12.6g} +- {r.errors[k]:.3g} {unit}")
        branching += r.x[k] if ".branching." in q.name else 0.0
    print(f"  {'remainder ' + decays[name]['remainder']:<32} {1.0 - branching:12.6g} (fraction, 1 - the branchings)")
    print("\nchi2 per data set at the fit:")
    for d in end.data:
        print(f"  {d.name:<10} {d.n:3d} points  chi2 = {d.chi2:8.3f}   ({d.quantity} in {d.unit}, "
              f"largest {np.max(d.model):.4g})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
