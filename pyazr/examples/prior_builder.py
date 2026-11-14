#!/usr/bin/env python3
"""A prior on the fit parameters built from named parts (pyazr.priors.build), checked as the engine checks it.

The parts here: a standard deviation of 5 % on every free normalization (around its start value), 10 % on the
first free reduced-width amplitude, and a correlation of 0.5 between that amplitude and the first normalization.
With --fit DIR (an earlier fit's output directory: its param.sav and covariance.dat, written by a fit with
--covariance-band), its free level energies and widths are taken from that fit instead, with their correlations.
The script prints the record (names, centres, sigmas, the correlation matrix, the checks, the notes), writes the file
into the project's directory with the parts on its #pyazr.prior line, sets the project's priorCovariance key and saves
the project as a copy (the original file is not changed).  The objective at the start point is printed without and
with the prior.

Usage (from the repository's root, on a copy of a project whose data files are beside it):

    python3 pyazr/examples/prior_builder.py project.azr [--fit A/output] [--out project_prior.azr]
"""

import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np  # noqa: E402

from pyazr import Project, Session, priors  # noqa: E402


def parts_for(parameters, fit_dir, directory):
    """The parts: the free norms at 5 %, then an earlier fit or the first free width at 10 % with a correlation."""
    free = parameters.free
    norms = [p.name for p in free if p.kind == "norm"]
    parts = [{"kind": "width", "names": norms, "percent": 5}] if norms else []
    if fit_dir:
        rel = os.path.relpath(fit_dir, directory)
        parts.append({"kind": "fit", "sav": os.path.join(rel, "param.sav"),
                      "covariance": os.path.join(rel, "covariance.dat")})
        return parts
    widths = [p.name for p in free if p.kind == "width" and p.start != 0]
    if widths:
        parts.append({"kind": "width", "names": widths[:1], "percent": 10})
        if norms:
            parts.append({"kind": "correlation", "pairs": [[widths[0], norms[0], 0.5]]})
    return parts


def objective_at_start(path):
    with Session(path) as s:
        x = np.asarray(s.parameters.x0, float)
        ev = s.evaluate(x)
        return s.objective()(x), ev.terms.get("prior covariance", 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr")
    ap.add_argument("--fit", default=None, help="an earlier fit's output directory (param.sav, covariance.dat)")
    ap.add_argument("--out", default=None, help="the project copy that names the file (default <name>_prior.azr)")
    args = ap.parse_args()

    p = Project.open(args.azr)
    with Session(args.azr) as s:
        parts = parts_for(s.parameters, args.fit, p.directory)
        if not parts:
            print("the project has no free normalization and no free width: nothing to put a prior on")
            return 1
        prior = priors.build(s.parameters, parts, directory=p.directory)
    print(prior)
    if prior.refusal:
        return 1
    name = "prior_parameters.txt"
    prior.use(p, name)                                # the file, and priorCovariance=prior_parameters.txt
    out = args.out or os.path.splitext(args.azr)[0] + "_prior.azr"
    p.save(out)
    back = priors.read_parts(os.path.join(p.directory, name))
    print(f"\nwritten: {name} (parts read back: {len(back['parts'])}); project: {out}")
    before, _ = objective_at_start(args.azr)
    after, term = objective_at_start(out)
    print(f"objective at the start point: {before:.10g} without, {after:.10g} with the prior (its term {term:.6g})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
