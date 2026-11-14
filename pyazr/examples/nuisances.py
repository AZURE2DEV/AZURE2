#!/usr/bin/env python3
"""Data sets that share one normalization with a prior and one energy shift with a prior, fitted by the engine.

A ``normalization`` record is one normalization shared by every data set whose ``norm`` names it; an
``energy_scale`` record is one shift of the lab energy shared the same way.  Each is written as a nuisance:
``value 1 fit prior 1 5%`` (fitted, with a Gaussian prior of 5 % about 1), ``offset 0 MeV fit prior 0 0.003``
(fitted, a prior of 3 keV about 0).  The priors' penalties are part of the objective (in its term "data,
normalizations and energy shifts", ``Evaluation.terms``; ``objective - chi2``) and the nuisances are rows of the
parameter table (kind ``norm``, ``shift``).  The same can be given on one data set alone:
``norm=Nuisance(1, "fit", prior=(1, "5%"))``.

The script adds the two records to a copy of the project, points the data sets at them, fits with the engine
(``fit_with_engine``) and prints the fitted nuisances with their errors and how far each moved from its prior's
centre, in units of the prior's width.  The copy goes to ``--out`` (default: a new temporary directory, printed);
its session reads the data files from the project's directory (``Session(..., cwd=)``).

Usage (from the repository's root):

    python3 pyazr/examples/nuisances.py tests/upgrade/T1/T1.azr el_lab090 el_lab125 el_lab160 [--out DIR]
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
    ap.add_argument("datasets", nargs="+", help="the data sets that share the normalization and the shift")
    ap.add_argument("--norm-prior", type=float, default=5.0, help="the normalization prior's width (%%)")
    ap.add_argument("--shift-prior", type=float, default=0.003, help="the shift prior's width (MeV, lab)")
    ap.add_argument("--out", default=None, help="the directory of the edited copy (default: a new temporary one)")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(args.project))
    out_dir = args.out or tempfile.mkdtemp(prefix="nuisances.")
    os.makedirs(out_dir, exist_ok=True)

    p = Project.open(args.project)
    p.add("normalization", "N_shared", value=f"1 fit prior 1 {args.norm_prior:g}%")
    p.add("energy_scale", "S_shared", offset=f"0 MeV fit prior 0 {args.shift_prior:g}")
    for name in args.datasets:
        p.set(name, norm="N_shared", shift="S_shared")
    path = p.save(os.path.join(out_dir, os.path.basename(args.project)))
    print(f"wrote {path}: {', '.join(args.datasets)} share N_shared and S_shared")

    with Session(path, cwd=here) as s:
        P = s.parameters
        rows = [q for q in P.free if q.kind in ("norm", "shift")]
        print(f"the nuisances among the free parameters: {[(q.name, q.kind, q.unit) for q in rows]}")
        start = s.evaluate()
        r = s.fit_with_engine()
        end = s.evaluate(r.x)
    print(f"fit: {r.status}; objective {start.objective:.4f} -> {r.objective:.4f}, data chi2 {end.chi2:.4f}, "
          f"the priors' part {end.objective - end.chi2:.4f}")
    centre_width = {"norm": (1.0, args.norm_prior / 100.0), "shift": (0.0, args.shift_prior)}
    for q in rows:
        k = r.names.index(q.name)
        c, w = centre_width[q.kind]
        print(f"  {q.name:<20} {r.x[k]:12.6g} +- {r.errors[k]:.3g} {q.unit:<4} ({(r.x[k] - c) / w:+.2f} prior widths)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
