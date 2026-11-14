#!/usr/bin/env python3
"""Scan a channel radius and watch what chi-squared does about it.

The channel radius is a boundary, not an observable, so the fit ought not to care much where it is put.  In practice
it does, because the external region is taken to be free of nuclear forces and it is not: a_c sets the hard-sphere
phase and the lower limit of the external-capture integral, so a scan of chi2(a_c) is the usual way of choosing it.

The radius lives in ``<levels>``, and a session reads its project once, so each radius is an edited copy of the
project (``Project.set_pair(pair, radius=)`` rewrites every row of the pair) and a new in-process session on it.
The copies are written to ``--out`` (default: a new temporary directory, printed); their sessions read the data files
from the original project's directory (``Session(..., cwd=)``).  The engine checks the signature of a saved
external-capture file and recomputes what belongs to another radius.  The widths are the project's physical widths
at every radius (no refit), so the scan shows how far the radius alone moves the model.

Usage (from the repository's root):

    python3 pyazr/examples/channel_radius_scan.py tests/upgrade/T3/T3.azr --pair 2 --radii 3.5 4.0 4.5 [--out DIR]
"""
import argparse
import os
import sys
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "2")

from pyazr import Project, Session  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr)")
    ap.add_argument("--pair", type=int, default=1, help="the particle pair's key")
    ap.add_argument("--radii", type=float, nargs="+", default=[3.5, 4.0, 4.5, 5.0], help="channel radii (fm)")
    ap.add_argument("--out", default=None, help="the directory of the edited copies (default: a new temporary one)")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(args.project))
    out_dir = args.out or tempfile.mkdtemp(prefix="channel_radius_scan.")
    os.makedirs(out_dir, exist_ok=True)
    print(f"the copies go to {out_dir}")

    radii = {q.key: q.radius for q in Project.open(args.project).pairs() if q.kind == 0}
    print(f"channel radii of {os.path.basename(args.project)} (fm): {radii}")
    print(f"\n{'a_c (fm)':>9} {'chi2':>12} {'chi2/N':>9}")

    rows = []
    for a in args.radii:
        p = Project.open(args.project)
        p.set_pair(args.pair, radius=a)
        path = p.save(os.path.join(out_dir, f"radius_{a:g}.azr"))
        with Session(path, cwd=here) as s:
            ev = s.evaluate()
            n = sum(d.n for d in ev.data)
        rows.append((a, ev.chi2, ev.chi2 / n))
        print(f"{a:9.2f} {ev.chi2:12.1f} {ev.chi2 / n:9.3f}")

    best = min(rows, key=lambda r: r[1])
    span = max(r[1] for r in rows) - best[1]
    print(f"\nminimum at a_c = {best[0]:.2f} fm, chi2 = {best[1]:.1f}; spread over the scan {span:.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
