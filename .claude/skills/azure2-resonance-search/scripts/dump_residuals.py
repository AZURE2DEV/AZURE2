#!/usr/bin/env python3
"""Point-by-point residuals of a fit, for the locators of the resonance search (rank_map.py, blur_map.py,
peakfind.py, blur_compare.py, blur_figure.py), from one evaluation of the project.

Reads:  the project (format 3); the fit from --sav (loaded by name) or --x-json ("x"), else the project's own values
        (a fitted project written by Session.save_fit carries them).
Writes: OUT (.npz), one entry per data point, in the engine's order: key (the data set's segment key), name, ecm
        (MeV, c.m., after the energy shift), ex (excitation energy, MeV), ang (c.m. angle, deg), dat and err (the
        datum and its error as the engine compares them), fit (n dat + res n err: the model in the same frame), res
        (the residual (f/J - n y)/(n dy) that enters chi-squared, not whitened; "model low" is res < 0).
        Prints the objective, the data chi-squared and sum(res^2); they agree unless a data covariance joins points
        (then sum(res^2) is the uncorrelated chi-squared, and the ranking ignores the correlations).

    python3 dump_residuals.py fit.azr resid.npz [--sav fit.sav]

Exit status: 0 written, 3 engine error.
"""
import argparse
import json
import sys

import numpy as np


def dump(s, x):
    """The columns (dict of arrays) of every data point at x."""
    ev = s.evaluate(x)
    cols = {k: [] for k in ("key", "name", "ecm", "ex", "ang", "dat", "err", "fit", "res")}
    for d in ev.data:
        n = d.n
        cols["key"] += [d.segment_key] * n
        cols["name"] += [d.name] * n
        cols["ecm"] += list(d.energies)
        cols["ex"] += list(d.excitation_energies)
        cols["ang"] += list(d.angles if len(d.angles) == n else np.full(n, np.nan))
        cols["dat"] += list(d.values)
        cols["err"] += list(d.errors)
        cols["fit"] += list(d.norm * (d.values + d.residuals * d.errors))
        cols["res"] += list(d.residuals)
    return ev, {k: np.asarray(v) for k, v in cols.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("project")
    ap.add_argument("out")
    ap.add_argument("--sav")
    ap.add_argument("--x-json")
    a = ap.parse_args(argv)
    import pyazr
    try:
        with pyazr.Session(a.project) as s:
            if a.sav:
                s.load_sav(a.sav)
            x = s.objective().x0
            if a.x_json:
                with open(a.x_json) as fh:
                    x = np.asarray(json.load(fh)["x"], float)
            ev, cols = dump(s, x)
    except pyazr.Error as e:
        print(f"dump_residuals: {e}", file=sys.stderr)
        return 3
    np.savez(a.out, **cols)
    print(f"objective {ev.objective:.6f}  data chi2 {ev.chi2:.6f}  rows {len(cols['res'])}  "
          f"sum res^2 {np.sum(cols['res'] ** 2):.6f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
