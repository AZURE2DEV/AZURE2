#!/usr/bin/env python3
"""A prior on the fit parameters built from named parts (pyazr.priors.build), checked as the engine checks it.

The parts here: a standard deviation of 5 % on every free normalization (around its start value), 10 % on the
first free reduced-width amplitude, and a correlation of 0.5 between that amplitude and the first normalization.
With --fit DIR (an earlier fit's output directory: its param.sav and covariance.dat, written by a fit with
--covariance-band), its free level energies and widths are taken from that fit instead, with their correlations.
The script prints the record (names, centres, sigmas, the correlation matrix, the checks, the notes), writes the file
beside the project copy ``--out`` (default: in a new temporary directory, printed) with the parts on its #pyazr.prior
line, sets the copy's priorCovariance key and saves the copy (the original project is not changed; the copy's session
reads the data files from the original's directory).  The objective at the start point is printed without and with
the prior.  In a project of your own, ``prior.use(project, "prior_parameters.txt")`` writes the file into the
project's directory and sets the key in one call.

Usage (from the repository's root):

    python3 pyazr/examples/prior_builder.py tests/upgrade/T1/T1.azr [--fit A/output] [--out DIR/T1_prior.azr]
"""

import argparse
import os
import sys
import tempfile

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


def objective_at_start(path, cwd):
    with Session(path, cwd=cwd) as s:
        x = np.asarray(s.parameters.x0, float)
        ev = s.evaluate(x)
        return s.objective()(x), ev.terms.get("prior covariance", 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr")
    ap.add_argument("--fit", default=None, help="an earlier fit's output directory (param.sav, covariance.dat)")
    ap.add_argument("--out", default=None, help="the project copy that names the file (default: in a temporary "
                                                "directory)")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(args.azr))
    out = os.path.abspath(args.out or os.path.join(tempfile.mkdtemp(prefix="prior_builder."),
                                                   os.path.basename(args.azr)))
    os.makedirs(os.path.dirname(out), exist_ok=True)

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
    name = os.path.join(os.path.dirname(out), "prior_parameters.txt")
    prior.write(name)                                 # the file; then the key that names it, in the copy
    p.set_config(priorCovariance=name if os.path.dirname(name) != here else os.path.basename(name))
    p.save(out)
    back = priors.read_parts(name)
    print(f"\nwritten: {name} (parts read back: {len(back['parts'])}); project: {out}")
    before, _ = objective_at_start(args.azr, here)
    after, term = objective_at_start(out, here)
    print(f"objective at the start point: {before:.10g} without, {after:.10g} with the prior (its term {term:.6g})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
