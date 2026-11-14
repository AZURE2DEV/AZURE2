#!/usr/bin/env python3
"""Least-squares fit of an AZURE2 project on the objective's rows and their analytic Jacobian (SciPy).

The objective is the program's: the data chi-squared (whitened where a data covariance joins points) plus every
prior and constraint, as rows r with sum(r**2) == Objective(x) (`Session.objective().residuals_and_jacobian`).
SciPy's trust-region reflective method takes the fit's limits (`Objective.bounds`) as bounds.

Reads:  the project (format 3) and, with --sav, a parameter file whose values are loaded by name first.
Writes: --out JSON: names, start, end point, objective and data chi-squared at both, per data set chi-squared,
        status, one entry per start; with --save-fit PATH the fitted project (PATH and PATH's .sav, written and
        verified by Session.save_fit).
Staged fits: --free REGEX (repeatable) fits only the free parameters whose names match; the others are held at the
start (their columns are passed as `skip`, so no finite difference is spent on them).
Multiple starts: --starts N draws N-1 more starts, each fitted parameter moved by jitter * N(0, 1) times its
parabolic error at the project's start (sqrt of the diagonal of (J^T J)^-1, the scale the data give it), and keeps
the lowest objective.

    python3 lsq_fit.py model.azr --out fit.json --save-fit fitted.azr
    python3 lsq_fit.py model.azr --free '^width_' --free 'norm$' --starts 8 --jitter 3 --out stage1.json

Exit status: 0 the best start converged (SciPy status > 0), 1 no start converged, 2 bad arguments, 3 engine error.
"""
import argparse
import json
import re
import sys

import numpy as np


def load(project, sav=None, **options):
    """A session of `project` (parameter file `sav` loaded by name first) and its objective."""
    import pyazr
    s = pyazr.Session(project, **options, parameter_file=sav)
    return s, s.objective()


def select(names, patterns):
    """Indices of the free parameters to fit: all when no pattern, else those matching any regex."""
    if not patterns:
        return list(range(len(names)))
    rx = [re.compile(p) for p in patterns]
    return [i for i, n in enumerate(names) if any(r.search(n) for r in rx)]


def inside(x, lo, hi):
    """x moved strictly inside finite bounds (least_squares needs a feasible start)."""
    span = np.where(np.isfinite(hi - lo), hi - lo, 1.0)
    return np.clip(x, lo + 1e-9 * span, hi - 1e-9 * span)


def fit_subset(f, x_start, cols, max_nfev=None, tol=1e-10):
    """Minimise f over the columns `cols` of the free vector, the others held at x_start.  Returns (x, scipy
    result, the columns dropped because no data set and no prior sees them at the start).

    Two guards a hand-written fit needs: a column that is identically zero at the start (a width no active data set
    sees, an amplitude at exactly 0) is held, since least_squares with x_scale="jac" would stop at once on it; and
    an engine failure at a trial point (Coulomb overflow, a singular level matrix) returns large rows with a zero
    Jacobian, so the trust region shrinks instead of the fit ending."""
    import pyazr
    from scipy.optimize import least_squares
    x_full = np.array(x_start, float)
    cols = np.asarray(cols, int)
    _, J0 = f.residuals_and_jacobian(x_full)
    dead = [int(c) for c in cols if not np.any(J0[:, c])]
    cols = np.array([c for c in cols if c not in dead], int)
    n_rows = J0.shape[0]
    held = [i for i in range(len(x_full)) if i not in set(cols.tolist())]
    lo = np.array([b[0] for b in f.bounds])[cols]
    hi = np.array([b[1] for b in f.bounds])[cols]
    cache = {}

    def rows_and_jac(z):
        key = z.tobytes()
        if key not in cache:
            cache.clear()
            x = x_full.copy()
            x[cols] = z
            try:
                r, J = f.residuals_and_jacobian(x, skip=held)
                cache[key] = (r, J[:, cols])
            except pyazr.Error:
                cache[key] = (np.full(n_rows, 1e6 / np.sqrt(n_rows)), np.zeros((n_rows, len(cols))))
        return cache[key]

    z0 = inside(x_full[cols], lo, hi)
    res = least_squares(lambda z: rows_and_jac(z)[0], z0, jac=lambda z: rows_and_jac(z)[1], bounds=(lo, hi),
                        method="trf", x_scale="jac", xtol=tol, ftol=tol, gtol=tol, max_nfev=max_nfev)
    x = x_full.copy()
    x[cols] = res.x
    return x, res, dead


def per_dataset(s, f, x):
    """Per data set: the program's chi-squared split and the whitened rows' split (they differ only where a data
    covariance joins data sets; only their sum over the group has a meaning then)."""
    ev = s.evaluate(x)
    r = f.residuals(x)
    return ev, [{"name": d.name, "points": d.n, "chi2": d.chi2, "chi2_whitened_rows": float(np.sum(r[d.rows] ** 2)),
                 "norm": d.norm} for d in ev.data]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("project", help="the .azr (format 3)")
    ap.add_argument("--sav", help="parameter file loaded by name before the fit")
    ap.add_argument("--free", action="append", default=[], metavar="REGEX",
                    help="fit only the free parameters whose names match (repeatable); default: all")
    ap.add_argument("--starts", type=int, default=1, help="number of starts (the first is the project's)")
    ap.add_argument("--jitter", type=float, default=3.0,
                    help="spread of the extra starts in parabolic errors at the start")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-nfev", type=int, default=None, help="SciPy's evaluation limit per start")
    ap.add_argument("--out", help="JSON result file")
    ap.add_argument("--save-fit", metavar="PATH",
                    help="write the fitted project (and its .sav) here, in the project's directory (its data "
                         "paths stay relative)")
    a = ap.parse_args(argv)
    if a.starts < 1:
        ap.error("--starts must be >= 1")
    import pyazr
    try:
        s, f = load(a.project, a.sav)
    except (pyazr.Error, OSError) as e:
        print(f"lsq_fit: {e}", file=sys.stderr)
        return 3
    cols = select(f.names, a.free)
    if not cols:
        print(f"lsq_fit: no free parameter matches {a.free}; the free parameters: {f.names}", file=sys.stderr)
        return 2
    x0 = f.x0
    rng = np.random.default_rng(a.seed)
    scale = np.zeros(len(x0))
    if a.starts > 1:
        _, J = f.residuals_and_jacobian(x0)
        Jc = J[:, cols]
        scale[cols] = np.sqrt(np.clip(np.diag(np.linalg.pinv(Jc.T @ Jc)), 0, None))
    runs, best = [], None
    for k in range(a.starts):
        start = x0.copy()
        if k:
            start[cols] = x0[cols] + a.jitter * scale[cols] * rng.standard_normal(len(cols))
        try:
            x, res, dead = fit_subset(f, start, cols, a.max_nfev)
            value = f(x)
            run = {"start": k, "objective": value, "status": int(res.status), "message": res.message,
                   "nfev": int(res.nfev), "njev": int(res.njev or 0), "held_dead": [f.names[c] for c in dead]}
            if dead:
                print(f"start {k}: held (zero Jacobian column): {run['held_dead']}")
        except pyazr.Error as e:
            x, value, run = None, np.inf, {"start": k, "objective": None, "status": -9, "message": str(e)}
        runs.append(run)
        print(f"start {k}: objective {run['objective']}  status {run['status']}  {run['message']}", flush=True)
        if x is not None and (best is None or value < best[1]):
            best = (x, value, run)
    if best is None:
        print("lsq_fit: every start failed in the engine", file=sys.stderr)
        return 3
    x, value, run = best
    ev, sets = per_dataset(s, f, x)
    result = {"project": a.project, "names": f.names, "fitted": [f.names[i] for i in cols], "x0": x0.tolist(),
              "x": x.tolist(), "objective_start": f(x0), "objective": value, "chi2": ev.chi2, "terms": ev.terms,
              "datasets": sets, "best_start": run["start"], "starts": runs}
    print(f"best: objective {value:.10g} (data chi2 {ev.chi2:.10g}) from start {run['start']}")
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(result, fh, indent=1)
    if a.save_fit:
        s.save_fit(a.save_fit, x)
    s.close()
    return 0 if run["status"] > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
