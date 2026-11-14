#!/usr/bin/env python3
"""Per-data-set chi-squared, pulls and their structure for an AZURE2 project at a parameter vector.

Reads a format-3 project (``.azr``) and, optionally, a parameter file (``param.sav``, matched by name); evaluates
it once with pyazr and prints, per data set: points N, chi2, chi2/N, the largest single-point chi2 and its share,
the fraction of pulls within 1 and 2 sigma, the median |pull| and the longest run of pulls of one sign. Then the
totals: the data chi2, the objective (data chi2 plus every prior), the free parameters and chi2/nu with
nu = N - (free parameters). With ``--chisquared-out`` it compares the data chi2 with the ``Total-Chi-Squared`` of a
CLI run's ``chiSquared.out``.

Pulls are ``Evaluation.data[i].residuals``: (model - datum*n)/(error*n), not whitened. A data set with a data
covariance has a chi2 that is not the sum of its squared pulls; the report marks it (``cov``).

Writes nothing unless ``--json PATH`` is given (the same numbers as JSON).

Exit status: 0 on success, 1 when ``--chisquared-out`` disagrees beyond ``--rtol``, 2 on a usage or set-up error.
"""
import argparse
import json
import math
import re
import sys

import numpy as np


def option_value(text):
    """'0'/'1'/'true'/'false' -> bool, a number -> float, else the text."""
    low = text.lower()
    if low in ("true", "yes", "on", "1"):
        return True
    if low in ("false", "no", "off", "0"):
        return False
    try:
        return float(text)
    except ValueError:
        return text


def open_session(azr, options):
    import pyazr
    opts = dict(o.split("=", 1) for o in options)
    return pyazr.Session(azr, **{k: option_value(v) for k, v in opts.items()})


def longest_run(signs):
    """Length of the longest run of equal nonzero signs."""
    best = run = 0
    prev = 0
    for s in signs:
        run = run + 1 if (s == prev and s != 0) else (1 if s != 0 else 0)
        prev = s
        best = max(best, run)
    return best


def report(session, x):
    ev = session.evaluate(x)
    rows = []
    for d in ev.data:
        r = np.asarray(d.residuals, float)
        chi2_pulls = float(np.sum(r ** 2))
        order = np.argsort(np.asarray(d.energies, float), kind="stable")
        pts = r ** 2
        rows.append(dict(
            name=d.name, quantity=d.quantity, n=int(d.n), chi2=float(d.chi2), norm=float(d.norm),
            max_point_chi2=float(pts.max()) if r.size else 0.0,
            max_point_energy_cm=float(np.asarray(d.energies)[int(np.argmax(pts))]) if r.size else None,
            within_1sigma=float(np.mean(np.abs(r) <= 1)) if r.size else None,
            within_2sigma=float(np.mean(np.abs(r) <= 2)) if r.size else None,
            median_abs_pull=float(np.median(np.abs(r))) if r.size else None,
            longest_same_sign_run=longest_run(np.sign(r[order])),
            covariance=not math.isclose(chi2_pulls, float(d.chi2), rel_tol=1e-6, abs_tol=1e-9),
        ))
    n_points = sum(r["n"] for r in rows)
    n_free = int(np.asarray(x).size)
    return dict(rows=rows, chi2=float(ev.chi2), objective=float(ev.objective), terms=dict(ev.terms),
                n_points=n_points, n_free=n_free,
                chi2_per_dof=float(ev.chi2) / (n_points - n_free) if n_points > n_free else None)


def total_from_chisquared_out(path):
    """The Total-Chi-Squared of a CLI run's chiSquared.out."""
    with open(path) as fh:
        for line in fh:
            m = re.search(r"Total-Chi-Squared:\s*(\S+)", line)
            if m:
                return float(m.group(1))
    raise ValueError(f"{path}: no Total-Chi-Squared line")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", help="the format-3 project")
    ap.add_argument("--params", help="parameter file (param.sav), applied by name; default: the project's values")
    ap.add_argument("--option", action="append", default=[], metavar="KEY=VALUE",
                    help="a session option (use_brune=0, use_long_wavelength=0, ...); the project's <config> keys "
                         "apply without it")
    ap.add_argument("--chisquared-out", help="a CLI run's chiSquared.out to compare the data chi2 with")
    ap.add_argument("--rtol", type=float, default=1e-6, help="relative tolerance of that comparison (default 1e-6)")
    ap.add_argument("--json", help="write the report as JSON to this path")
    a = ap.parse_args(argv)
    try:
        s = open_session(a.azr, a.option)
    except Exception as err:  # the engine's refusal names the problem
        print(f"error: {err}", file=sys.stderr)
        return 2
    with s:
        x = s.parameters.read(a.params) if a.params else s.parameters.x0
        rep = report(s, x)
    print(f"{'data set':24s} {'N':>5s} {'chi2':>11s} {'chi2/N':>8s} {'norm':>7s} {'max pt':>9s} {'share':>6s} "
          f"{'<1s':>5s} {'<2s':>5s} {'med|r|':>6s} {'run':>4s}")
    for r in rep["rows"]:
        share = r["max_point_chi2"] / r["chi2"] if r["chi2"] > 0 else 0.0
        print(f"{r['name'][:24]:24s} {r['n']:5d} {r['chi2']:11.4f} {r['chi2'] / max(r['n'], 1):8.3f} "
              f"{r['norm']:7.4f} {r['max_point_chi2']:9.3f} {share:6.2f} {r['within_1sigma']:5.2f} "
              f"{r['within_2sigma']:5.2f} {r['median_abs_pull']:6.2f} {r['longest_same_sign_run']:4d}"
              + ("  cov" if r["covariance"] else ""))
    print(f"\ndata chi2 {rep['chi2']:.6f}   objective {rep['objective']:.6f}   N {rep['n_points']}   "
          f"free {rep['n_free']}   chi2/nu " + ("-" if rep["chi2_per_dof"] is None else f"{rep['chi2_per_dof']:.4f}"))
    for k, v in rep["terms"].items():
        if v:
            print(f"  term {k}: {v:.6f}")
    status = 0
    if a.chisquared_out:
        ref = total_from_chisquared_out(a.chisquared_out)
        ok = math.isclose(rep["chi2"], ref, rel_tol=a.rtol, abs_tol=1e-9)
        print(f"chiSquared.out Total-Chi-Squared {ref:.6f}: {'agrees' if ok else 'DIFFERS'} (rtol {a.rtol:g})")
        status = 0 if ok else 1
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=1)
    return status


if __name__ == "__main__":
    sys.exit(main())
