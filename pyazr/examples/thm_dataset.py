#!/usr/bin/env python3
"""THM data sets in a fit: the experiment that shares one profiled norm, its chi-squared beside the direct data's.

A Trojan-horse (THM) data set measures a quasi-free two-body cross section inside a three-body reaction
(A + a(x+s) -> s + b + B).  Its absolute scale is unknown, so its norm is ``profiled``: at every evaluation the engine
solves for the norm (and a background, where the experiment has one) that minimises the chi-squared, and the data sets
of one experiment (``<thm>`` ``experiment[<name>]``) share it.  ``Session.thm_experiments(x)`` reports each experiment
at ``x``: its data sets, points, chi2, the profiled norm (in the unit of the data per that of the model) with its
1-sigma error, the background coefficients.  The THM options and the data sets' kinematics are read with
``Project.thm()``; the user guide's "THM" pages describe them.

Usage (from the repository's root):

    python3 pyazr/examples/thm_dataset.py tests/upgrade/T2/T2.azr
"""
import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "2")

from pyazr import Project, Session  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="a format-3 project with THM data sets")
    args = ap.parse_args()

    p = Project.open(args.project)
    print(f"<thm> of {os.path.basename(args.project)}: {p.thm()}")
    records = p.records()
    with Session(args.project) as s:
        ev = s.evaluate()
        experiments = s.thm_experiments(ev.x)
    for name, e in experiments.items():
        members = [d for d in ev.data if d.segment_key in e["segments"]]
        print(f"\nexperiment {name}: data sets {[d.name for d in members]}, {e['points']} points, "
              f"background {e['background']}")
        print(f"  norm ({e['status']}) = {e['norm']:.6g} +- {e['sigma'][0]:.3g} (data unit per model unit), "
              f"chi2 = {e['chi2']:.3f}")
        for d in members:
            print(f"  {d.name:<16} {records[d.name]['quantity']:<12} {d.n:3d} points  chi2 = {d.chi2:9.3f}")
    direct = [d for d in ev.data if not any(d.segment_key in e["segments"] for e in experiments.values())]
    print(f"\nthe other data sets: {len(direct)}, chi2 = {sum(d.chi2 for d in direct):.3f}; "
          f"all data: chi2 = {ev.chi2:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
