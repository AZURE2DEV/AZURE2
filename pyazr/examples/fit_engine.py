#!/usr/bin/env python3
"""Fit a project with the engine's own fitter, MINOS intervals included, and save the fit as a project that reopens
as the fit.

``Session.fit_with_engine(method, minos=)`` runs the program's fit in-process: the same fitters as ``AZURE2 fit``
(MIGRAD by default, then HESSE; ``pyazr.kinds()["fitters"]`` lists them), ``<parameterSettings>`` applied, no file
written.  It returns the end point ``x``, the ``objective`` (the program's ``Total Chi-Squared``), ``status``,
``covariance``, ``errors`` and, with ``minos=True``, the asymmetric intervals ``minos`` with Minuit2's verdict per
side (``minos_valid``).  ``Session.save_fit(path, x=)`` writes the fit: the ``<levels>`` lines as physical values
(partial widths in eV, ANCs in fm^-1/2), the fitted nuisances into their records, and a ``param.sav`` beside it;
it reopens what it wrote and raises if that is not the fit.  ``--params`` saves a fit made elsewhere (a
``param.sav`` of ``AZURE2 fit``) instead of fitting.

The fit is written to ``--out`` (default: a new temporary directory, printed); ``save_fit`` writes the data files'
paths relative to the new place, so the saved project opens alone.

Usage (from the repository's root):

    python3 pyazr/examples/fit_engine.py tests/upgrade/T1/T1.azr [--method migrad] [--no-minos] [--out DIR]
    python3 pyazr/examples/fit_engine.py tests/upgrade/T1/T1.azr --params output/param.sav [--out DIR]
"""
import argparse
import os
import sys
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "2")

from pyazr import Session  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr)")
    ap.add_argument("--method", default="migrad", help="a fitter of pyazr.kinds()['fitters']")
    ap.add_argument("--no-minos", action="store_true", help="skip MINOS")
    ap.add_argument("--params", default=None, help="save this parameter file's fit instead of fitting")
    ap.add_argument("--out", default=None, help="the directory of the saved fit (default: a new temporary one)")
    args = ap.parse_args()
    out_dir = args.out or tempfile.mkdtemp(prefix="fit_engine.")
    os.makedirs(out_dir, exist_ok=True)

    with Session(args.project, parameter_file=args.params) as s:
        unit = {q.name: q.unit for q in s.parameters.free}
        start = s.evaluate()
        if args.params:
            x, total = s.parameters.x0, start.objective
            print(f"{args.params}: objective {total:.6f}")
        else:
            r = s.fit_with_engine(args.method, minos=not args.no_minos)
            x, total = r.x, r.objective
            print(f"{args.method}: {r.status} after {r.calls} calls; objective {start.objective:.6f} -> {total:.6f}")
            for note in r.notices:
                print("  note:", note)
            print(f"{'parameter':<24} {'value':>14} {'error':>10} {'MINOS lower':>12} {'upper':>10}  unit")
            for name, value, error in zip(r.names, r.x, r.errors):
                lo, hi = r.minos.get(name, (float("nan"), float("nan"))) if r.minos else (float("nan"),) * 2
                print(f"{name:<24} {value:14.7g} {error:10.3g} {lo:12.3g} {hi:10.3g}  {unit[name]}")
        path, sav = s.save_fit(os.path.join(out_dir, os.path.basename(args.project)), x=x)
    print(f"\nwrote {path}\nwrote {sav}")

    with Session(path) as t:          # the saved project alone (save_fit wrote its data paths relative to it)
        again = t.evaluate().objective
    print(f"objective of the saved project as it stands: {again:.6f} (the fit: {total:.6f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
