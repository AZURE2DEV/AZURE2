#!/usr/bin/env python3
"""Fisher information of an AZURE2 fit: F = J^T J of the objective's rows (whitened data rows, then one row per
prior and constraint) at a point x, its eigen-spectrum, the sloppy directions, the information per data set, the
expected covariance F^-1 (against the engine's HESSE covariance on request), and what a better measurement of one
data set would add.

With r the rows (sum r^2 = chi^2 + priors) and J = dr/dx, the Gauss-Newton Hessian of the objective is 2 F, so
F^-1 is the covariance at Delta chi^2 = 1 (errordef 1) that Levenberg-Marquardt reports.  The data rows of data set
i are Evaluation.data[i].rows; F_i = J_i^T J_i sums to F with the prior rows' part F_prior.  Where a data
covariance joins data sets their whitened rows mix the sets: F_i is then the information of the rows, and only the
sum over the joined group has a meaning (the script says so).

Reads:  the project; the point from --x-json ("x", e.g. lsq_fit.py's output) or --sav, else the project's start.
Writes: a text report on standard output; --out NPZ with F, the per-set F_i, F_prior, the eigenvalues and vectors
        (of the scaled matrix), names; --json a summary.
        --compare-engine METHOD runs Session.fit_with_engine(METHOD) from x and compares its errors with sqrt(diag
        F^-1).  --improve SET=k shows the errors when the errors of data set SET are divided by k
        (F + (k^2 - 1) F_SET).

    python3 fisher.py model.azr --x-json fit.json --improve "sigma_pa=2" --compare-engine lm

Exit status: 0 done, 1 F is singular (an unconstrained direction; listed), 2 bad arguments, 3 engine error.
"""
import argparse
import json
import sys

import numpy as np


def fisher(s, f, x):
    """(F, {set name: F_i}, F_prior, rows info) at x."""
    r, J = f.residuals_and_jacobian(x)
    ev = s.evaluate(x)
    n_data = max((d.rows.stop for d in ev.data), default=0)
    per_set = {d.name: J[d.rows].T @ J[d.rows] for d in ev.data}
    F_prior = J[n_data:].T @ J[n_data:]
    joined = [d.name for d in ev.data if abs(np.sum(r[d.rows] ** 2) - d.chi2) > 1e-6 * max(1.0, d.chi2)]
    return J.T @ J, per_set, F_prior, {"n_data_rows": n_data, "n_rows": len(r), "joined": joined,
                                      "objective": float(np.sum(r ** 2))}


def scaled(F):
    """D^-1/2 F D^-1/2 with D = diag F: the spectrum without the parameters' units."""
    d = np.sqrt(np.clip(np.diag(F), 1e-300, None))
    return F / np.outer(d, d), d


def errors(F):
    """sqrt(diag F^-1), or None when F is singular."""
    try:
        C = np.linalg.inv(F)
    except np.linalg.LinAlgError:
        return None
    dg = np.diag(C)
    return np.sqrt(dg) if np.all(dg > 0) else None


def leave_out(F, F_i):
    """Errors without the information of one set (inf where F - F_i is singular in that parameter)."""
    e = errors(F - F_i)
    return e if e is not None else np.full(len(F), np.inf)


def report(names, F, per_set, F_prior, info, top=3, out=print):
    S, d = scaled(F)
    lam, V = np.linalg.eigh(S)
    out(f"{len(names)} free parameters, {info['n_data_rows']} data rows, {info['n_rows'] - info['n_data_rows']} "
        f"prior rows; objective {info['objective']:.8g}")
    out(f"scaled Fisher matrix: eigenvalues {lam.min():.3g} .. {lam.max():.3g}, condition {lam.max() / lam.min():.3g}")
    for k in list(range(min(top, len(lam)))):
        comp = np.argsort(-np.abs(V[:, k]))[:3]
        out(f"  sloppy {k}: eigenvalue {lam[k]:.3g}  " + ", ".join(f"{V[c, k]:+.2f} {names[c]}" for c in comp))
    e = errors(F)
    if e is None:
        out("F is singular: the sloppiest direction above is not constrained by data or priors")
        return lam, V, None
    out(f"{'parameter':28s} {'sigma':>11s} " + " ".join(f"{n[:10]:>10s}" for n in per_set) + f" {'priors':>10s}")
    for j, n in enumerate(names):
        share = [Fi[j, j] / F[j, j] for Fi in per_set.values()] + [F_prior[j, j] / F[j, j]]
        out(f"{n:28s} {e[j]:11.4g} " + " ".join(f"{v:10.3f}" for v in share))
    out("(columns: share of the parameter's information F_jj from each data set and from the priors)")
    if info["joined"]:
        out(f"data sets joined by a data covariance (their shares are those of the whitened rows): {info['joined']}")
    return lam, V, e


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("project")
    ap.add_argument("--x-json", help="JSON with the point under 'x' (free parameters' order)")
    ap.add_argument("--sav", help="parameter file loaded by name")
    ap.add_argument("--improve", action="append", default=[], metavar="SET=k",
                    help="errors with data set SET's errors divided by k (repeatable)")
    ap.add_argument("--compare-engine", metavar="METHOD", help="fit_with_engine(METHOD) from x; compare errors")
    ap.add_argument("--top", type=int, default=3, help="sloppy directions to print")
    ap.add_argument("--out", help="NPZ with the matrices")
    ap.add_argument("--json", help="JSON summary")
    a = ap.parse_args(argv)
    import pyazr
    try:
        s = pyazr.Session(a.project)
        if a.sav:
            s.load_sav(a.sav)
        f = s.objective()
        x = f.x0
        if a.x_json:
            with open(a.x_json) as fh:
                x = np.asarray(json.load(fh)["x"], float)
        if len(x) != len(f.names):
            ap.error(f"the point has {len(x)} values; the project has {len(f.names)} free parameters")
        F, per_set, F_prior, info = fisher(s, f, x)
    except (pyazr.Error, OSError) as e:
        print(f"fisher: {e}", file=sys.stderr)
        return 3
    lam, V, e = report(f.names, F, per_set, F_prior, info, a.top)
    summary = {"names": f.names, "x": x.tolist(), "sigma": None if e is None else e.tolist(),
               "eigenvalues_scaled": lam.tolist(), "joined": info["joined"],
               "leave_one_out_sigma": {k: leave_out(F, Fi).tolist() for k, Fi in per_set.items()}}
    for spec in a.improve:
        name, _, k = spec.rpartition("=")
        if name not in per_set:
            print(f"fisher: no data set {name!r}; the data sets: {list(per_set)}", file=sys.stderr)
            return 2
        e2 = errors(F + (float(k) ** 2 - 1) * per_set[name])
        summary.setdefault("improve", {})[spec] = None if e2 is None else e2.tolist()
        if e is not None and e2 is not None:
            print(f"errors of {name} / {k}: " + ", ".join(f"{n} {v:.3g}" for n, v in zip(f.names, e2 / e)))
    if a.compare_engine:
        r = s.fit_with_engine(a.compare_engine, x0=x)
        ee = None if r.errors is None else np.asarray(r.errors)
        summary["engine"] = {"method": a.compare_engine, "objective": r.objective, "status": r.status,
                             "errors": None if ee is None else ee.tolist()}
        if ee is not None and e is not None:
            print(f"engine {a.compare_engine} ({r.status}, objective {r.objective:.8g}) errors / Fisher errors: "
                  + ", ".join(f"{n} {v:.4f}" for n, v in zip(f.names, ee / e)))
    if a.out:
        np.savez(a.out, F=F, F_prior=F_prior, names=np.array(f.names), x=x, eigenvalues=lam, eigenvectors=V,
                 **{"F_set_" + k: v for k, v in per_set.items()})
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(summary, fh, indent=1)
    s.close()
    return 0 if e is not None else 1


if __name__ == "__main__":
    sys.exit(main())
