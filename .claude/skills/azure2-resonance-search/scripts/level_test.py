#!/usr/bin/env python3
"""Test changes of the level scheme of an AZURE2 fit, each refitted, with their significance: add a level (at one
energy or a scan of energies), remove one, switch one off, move one, flip the sign of a width.

Each variant is a `Project` edit of the fitted project (levels by label, "<J^pi> <E_x as written>"), saved under
--workdir and refitted from its own values (least squares on the objective's rows and Jacobian, all free
parameters) in a session whose relative paths are the base project's.  Structural edits renumber the level
parameters, so the base must carry its fit in its records: give a fitted project (Session.save_fit), or --x-json /
--sav and the script writes that fit into a copy first.

Significance against the base (N data points, k free parameters, objective = chi^2 + priors at the refit):
  Delta = objective_smaller_model - objective_larger_model, Delta k free parameters between them;
  AIC = objective + 2k, BIC = objective + k ln N (lower is better);
  F = (Delta / Delta k) / (objective_larger / (N - k_larger)), p from the F distribution (the nested-model test;
  valid when the smaller model is the larger with parameters at fixed values, e.g. widths at 0);
  --hold-out SET: every model is refitted with data set SET off, and the SET's chi-squared at that fit is the
  prediction score (cross-validation on held-out data; the set's norm stays at its start).

Variants (repeat --variant):
  add:<J^pi>:<E_x>:<pair>/<l>/<s>/<width>[,...]     widths in the row convention (eV for an open channel), free
  scan:<J^pi>:<E1>,<E2>,...:<pair>/<l>/<s>/<width>[,...]   one add per energy
  remove:<label>      off:<label>      move:<label>:<E_x>      sign:<label>:<pair>/<l>/<s>

Reads:  the base project (and its data, relative to its directory).
Writes: the variants under --workdir; --out JSON with one row per model; a table on standard output.

    python3 level_test.py fitted.azr --variant "scan:3-:13.2,13.27,13.4:1/2/1/1700,2/3/0/14000" --out scan.json
    python3 level_test.py fitted.azr --variant "remove:3- 13.2692" --hold-out sigma_pa --out remove.json

Exit status: 0 done, 2 bad arguments, 3 engine error in the base (a failing variant is reported, not fatal).
"""
import argparse
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import f as f_dist


def refit(f, max_nfev=200):
    """Least squares from the session's start over every free parameter; returns (x, objective)."""
    lo = np.array([b[0] for b in f.bounds])
    hi = np.array([b[1] for b in f.bounds])
    span = np.where(np.isfinite(hi - lo), hi - lo, 1.0)
    x0 = np.clip(f.x0, lo + 1e-9 * span, hi - 1e-9 * span)
    if len(x0) == 0:
        return x0, f(x0)
    r = least_squares(f.residuals, x0, jac=f.jacobian, bounds=(lo, hi), method="trf", x_scale="jac", xtol=1e-10,
                      ftol=1e-10, gtol=1e-10, max_nfev=max_nfev)
    return r.x, f(r.x)


def parse_channels(text):
    rows = []
    for item in text.split(","):
        pair, l, s, w = item.split("/")
        rows.append({"pair": int(pair), "l": int(l), "s": float(s), "width": float(w), "fixed": False})
    return rows


def jpi(text):
    """'3/2-' -> (1.5, -1)."""
    j, p = text[:-1], text[-1]
    num, _, den = j.partition("/")
    return float(num) / float(den or 1), (1 if p == "+" else -1)


def edits(spec):
    """[(name, function(project))] for one --variant."""
    kind, _, rest = spec.partition(":")
    if kind in ("add", "scan"):
        jp, energies, channels = rest.split(":")
        J, parity = jpi(jp)
        rows = parse_channels(channels)
        return [(f"add {jp} {float(e):g}", lambda p, e=float(e): p.add_level(J, parity, e, rows, fixed=False))
                for e in energies.split(",")]
    if kind == "remove":
        return [(f"remove {rest}", lambda p: p.remove_level(rest))]
    if kind == "off":
        return [(f"off {rest}", lambda p: p.set_level(rest, active=False))]
    if kind == "move":
        label, e = rest.rsplit(":", 1)
        return [(f"move {label} to {e}", lambda p: p.set_level(label, energy=float(e)))]
    if kind == "sign":
        label, ch = rest.rsplit(":", 1)
        pair, l, s = ch.split("/")

        def flip(p):
            w = p.level(label).channel(int(pair), int(l), float(s)).width
            p.set_channel(label, int(pair), int(l), float(s), width=-w)
        return [(f"sign {label} {ch}", flip)]
    raise ValueError(f"unknown variant kind {kind!r} in {spec!r}")


def score(path, cwd, hold_out=None):
    """Refit the project at `path`; with hold_out, refit without that data set and score it at that fit."""
    import pyazr
    from pyazr import Project
    out = {}
    with pyazr.Session(path, cwd=cwd) as s:
        f = s.objective()
        x, obj = refit(f)
        ev = s.evaluate(x)
        out.update(objective=obj, chi2=ev.chi2, k=len(x), N=sum(d.n for d in ev.data), x=x.tolist(), names=f.names)
    if hold_out:
        p = Project.open(path)
        p.set(hold_out, active="no")
        reduced = p.save(path[:-4] + ".holdout.azr")
        with pyazr.Session(reduced, cwd=cwd) as s:
            f = s.objective()
            x, _ = refit(f)
            fitted = dict(zip(f.names, x))
        with pyazr.Session(path, cwd=cwd) as s:
            f = s.objective()
            x_full = np.array([fitted.get(n, v) for n, v in zip(f.names, f.x0)])
            out["held_out_chi2"] = s.evaluate(x_full).dataset(hold_out).chi2
    return out


def compare(base, row):
    """Delta, AIC, BIC and the F test between the base and a variant (whichever has more parameters is larger)."""
    big, small = (row, base) if row["k"] > base["k"] else (base, row)
    dk = big["k"] - small["k"]
    delta = small["objective"] - big["objective"]
    out = {"delta_objective": base["objective"] - row["objective"], "delta_k": row["k"] - base["k"],
           "aic": row["objective"] + 2 * row["k"], "bic": row["objective"] + row["k"] * math.log(row["N"])}
    if dk > 0 and big["N"] > big["k"]:
        F = (delta / dk) / (big["objective"] / (big["N"] - big["k"]))
        out.update(F=F, p_value=float(f_dist.sf(F, dk, big["N"] - big["k"])))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("project", help="the fitted base project")
    ap.add_argument("--variant", action="append", default=[], required=True)
    ap.add_argument("--x-json", help="write this fit ('x') into a copy of the base first")
    ap.add_argument("--sav", help="load this parameter file and write it into a copy of the base first")
    ap.add_argument("--hold-out", metavar="SET", help="data set held out for cross-validation")
    ap.add_argument("--workdir", default="level_test", help="where the variant projects are written")
    ap.add_argument("--out", help="JSON result")
    a = ap.parse_args(argv)
    import pyazr
    from pyazr import Project
    os.makedirs(a.workdir, exist_ok=True)
    cwd = os.path.dirname(os.path.abspath(a.project))
    base_path = os.path.abspath(a.project)
    try:
        variants = [e for spec in a.variant for e in edits(spec)]
    except ValueError as e:
        print(f"level_test: {e}", file=sys.stderr)
        return 2
    try:
        if a.x_json or a.sav:
            with pyazr.Session(base_path, parameter_file=a.sav) as s:
                x = s.objective().x0
                if a.x_json:
                    with open(a.x_json) as fh:
                        x = np.asarray(json.load(fh)["x"], float)
                base_path = os.path.abspath(os.path.join(a.workdir, "base_fit.azr"))
                s.save_fit(base_path, x)
        base = score(base_path, cwd, a.hold_out)
    except pyazr.Error as e:
        print(f"level_test: the base project: {e}", file=sys.stderr)
        return 3
    rows = [dict(name="base", **base)]
    print(f"{'model':32s} {'objective':>12s} {'k':>4s} {'Delta':>10s} {'AIC':>12s} {'BIC':>12s} {'p(F)':>10s}"
          + (f" {'held-out chi2':>14s}" if a.hold_out else ""))
    print(f"{'base':32s} {base['objective']:12.6g} {base['k']:4d}" + " " * 59
          + (f" {base['held_out_chi2']:14.6g}" if a.hold_out else ""))
    for k, (name, edit) in enumerate(variants):
        p = Project.open(base_path)
        try:
            edit(p)
            path = p.save(os.path.abspath(os.path.join(a.workdir, f"variant_{k}.azr")))
            row = score(path, cwd, a.hold_out)
        except (pyazr.Error, KeyError) as e:
            print(f"{name:32s} failed: {e}")
            rows.append({"name": name, "error": str(e)})
            continue
        stats = compare(base, row)
        rows.append(dict(name=name, project=path, **row, **stats))
        print(f"{name:32s} {row['objective']:12.6g} {row['k']:4d} {stats['delta_objective']:10.4g} "
              f"{stats['aic']:12.6g} {stats['bic']:12.6g} {stats.get('p_value', float('nan')):10.3g}"
              + (f" {row['held_out_chi2']:14.6g}" if a.hold_out else ""))
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(rows, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
