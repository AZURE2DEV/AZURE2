#!/usr/bin/env python3
"""The hybrid Coulomb model, switched on for one particle pair at a time.

AZURE2 can replace the pure Coulomb radial wave functions of a channel with the
solution through a nuclear potential -- Woods-Saxon or Gaussian -- integrated
outward by Numerov and matched to the Coulomb functions outside.  A potential
belongs to a *pair*: it bends the wave functions of that channel and no other.
In a project it is a ``potential`` record of the ``<potentials>`` section, one
per pair; a change of it is an edit of the :class:`pyazr.Project` and a new
session, as a change of a channel radius is.

This example, on a project with a ``potential`` record for the pair (``tests/hybrid_potential`` is 13N cut down to
its three elastic p + 12C segments; a new session costs about two seconds at one thread):

  1. reads back what every pair resolves to,
  2. scans the Woods-Saxon depth on the pair and prints chi-squared, one edited copy of the project and one session
     per depth (the copies in ``--out``, default a new temporary directory, removed at the end; their sessions read
     the data files from the project's directory),
  3. shows the one thing that will bite -- the fit vector of a session without the potential does not describe the
     model with it.

Usage (from the repository's root):

    python3 pyazr/examples/nuclear_potential.py tests/hybrid_potential/hybrid_potential.azr [--pair 1]
"""
import argparse
import os
import shutil
import sys
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "2")   # before numpy

import numpy as np  # noqa: E402

from pyazr import Project, Session  # noqa: E402


def chi2(model, x=None):
    """Chi-squared at ``x``, by default the session's start ``parameters.x0``."""
    return model.evaluate(np.asarray(model.parameters.x0 if x is None else x, float)).chi2


def variant(project, record, work, name, depth):
    """A copy of the project with the pair's potential at ``depth`` MeV (no record at all for ``None``)."""
    p = Project.open(project)
    if depth is None:
        p.remove(record)
    else:
        p.set(record, central=f"woods_saxon depth {depth} MeV radius 3.6 fm diffuseness 0.6 fm")
    return p.save(os.path.join(work, name))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="a format-3 project with a potential record for the pair")
    ap.add_argument("--pair", type=int, default=1, help="the particle pair's key")
    ap.add_argument("--out", default=None, help="the directory of the edited copies (default: a temporary one)")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(args.project))
    project = Project.open(args.project)
    pair = next(q for q in project.pairs() if q.key == args.pair)
    record = next(name for name, r in project.records().items()       # a potential names its pair by key or label
                  if r.record_kind == "potential" and r["pair"] in (str(pair.key), pair.label))
    work = args.out or tempfile.mkdtemp(prefix="nuclear_potential.")
    os.makedirs(work, exist_ok=True)
    try:
        with Session(args.project) as m:
            print("what each pair resolves to as the file is written")
            print("-" * 62)
            for key, setting in m.nuclear_potentials().items():
                who = "default" if key == 0 else f"pair {key}"
                print(f"  {who:<10} {setting}")

        with Session(variant(args.project, record, work, "none.azr", None), cwd=here) as m:
            baseline = chi2(m)
            stale = np.asarray(m.parameters.x0, float)  # the fit vector of the model without a potential
        print(f"\n  chi-squared with no hybrid potential: {baseline:.3f}")

        print(f"\nWoods-Saxon depth scan on pair {args.pair} (R = 3.6 fm, a = 0.6 fm)")
        print("-" * 62)
        print(f"  {'V0 (MeV)':>10}  {'chi-squared':>16}  {'vs baseline':>12}")
        for V0 in (5.0, 10.0, 20.0, 30.0):
            with Session(variant(args.project, record, work, f"V{V0:g}.azr", V0), cwd=here) as m:
                c = chi2(m)
            print(f"  {V0:>10.1f}  {c:>16.3f}  {c / baseline:>11.4f}x")

        print("\nthe fit vector belongs to its session -- read it from the new one")
        print("-" * 62)
        with Session(variant(args.project, record, work, "V20.azr", 20.0), cwd=here) as m:
            with_stale = chi2(m, stale)
            with_fresh = chi2(m)
        print(f"  the vector of the session without it : {with_stale:.3f}   <- wrong")
        print(f"  this session's parameters.x0         : {with_fresh:.3f}   <- right")
        print("\n  The two differ because the project holds physical widths and a")
        print("  reduced-width amplitude means something different once the")
        print("  penetrabilities have moved.  Same caveat as a channel-radius change,")
        print("  and the same fix: take the new session's vector, then refit.")
    finally:
        if not args.out:
            shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
