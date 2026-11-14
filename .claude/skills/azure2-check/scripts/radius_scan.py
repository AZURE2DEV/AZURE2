#!/usr/bin/env python3
"""Channel-radius scan with a refit at every radius.

The channel radius a_c is a model choice: a reduced width means something else at another radius, so chi2 at a new
radius means nothing until the model is refitted there. For every radius this script writes a copy of the project
with ``Project.set_pair(pair, radius=a)`` (beside the original, so that relative data paths hold, with its own output
directory so that no external-capture cache of another radius is reused), opens a session on it, evaluates the
objective at the start, refits, and records the objective and the data chi2 at the minimum. The radii are walked
outward from the one closest to the project's own radius and each fit starts from its converged neighbour (by
parameter name), which keeps the ends of the scan from being biased by one distant start.

Fitters: ``scipy`` (default; ``scipy.optimize.least_squares`` on ``Objective.residuals`` and ``.jacobian``, trust
region, the fit's bounds) or ``engine-lm`` / ``engine-migrad`` (``Session.fit_with_engine``, the program's own fitters).

Writes: per radius ``<stem>_a<pair>_<radius>.azr`` and ``<stem>_a<pair>_<radius>.sav`` (the fitted vector, by name)
and the output directory ``output_a<pair>_<radius>/`` beside the project; ``--clean`` removes them. ``--table PATH``
writes the scan as text columns (a_c fm, objective at start, objective, data chi2, free parameters, fit status).

Exit status: 0 when every radius was fitted, 1 when a radius failed (the others are still reported), 2 on usage.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np


def option_value(text):
    low = text.lower()
    if low in ("true", "yes", "on", "1"):
        return True
    if low in ("false", "no", "off", "0"):
        return False
    try:
        return float(text)
    except ValueError:
        return text


def refit(session, x0, fitter, max_nfev):
    """(x, objective, status) of a fit of the session's objective from x0."""
    f = session.objective()
    if fitter.startswith("engine-"):
        r = session.fit_with_engine(fitter.split("-", 1)[1], x0=x0)
        return np.asarray(r.x, float), float(f(np.asarray(r.x, float))), str(r.status)
    from scipy.optimize import least_squares
    lo = np.array([b[0] for b in f.bounds], float)
    hi = np.array([b[1] for b in f.bounds], float)
    x0 = np.clip(np.asarray(x0, float), np.nextafter(lo, hi), np.nextafter(hi, lo))
    r = least_squares(f.residuals, x0, jac=f.jacobian, bounds=(lo, hi), x_scale="jac", method="trf",
                      xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=max_nfev)
    return r.x, float(f(r.x)), f"scipy status {r.status}"


def variant_paths(azr, pair, radius):
    tag = f"a{pair}_{radius:g}"
    return azr.with_name(f"{azr.stem}_{tag}.azr"), azr.with_name(f"{azr.stem}_{tag}.sav"), f"output_{tag}/"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--pair", type=int, required=True, help="pair key of <levels> whose radius is scanned")
    ap.add_argument("--radii", type=float, nargs="+", required=True, help="radii in fm")
    ap.add_argument("--params", dest="parameter_file", help="start parameter file (param.sav), by name")
    ap.add_argument("--option", action="append", default=[], metavar="KEY=VALUE", help="a session option")
    ap.add_argument("--fitter", default="scipy", choices=["scipy", "engine-lm", "engine-migrad"])
    ap.add_argument("--max-nfev", type=int, default=200, help="scipy: largest number of function evaluations")
    ap.add_argument("--no-refit", action="store_true", help="only the objective at the start (not a result)")
    ap.add_argument("--table", help="write the scan to this text file")
    ap.add_argument("--clean", action="store_true", help="remove the variant files and output directories")
    a = ap.parse_args(argv)
    import pyazr

    opts = {k: option_value(v) for k, v in (o.split("=", 1) for o in a.option)}
    azr = a.azr.resolve()
    base = pyazr.Project.open(azr)
    own = {q.key: q.radius for q in base.pairs()}
    if a.pair not in own:
        print(f"error: pair {a.pair} not in {azr.name} (pairs {sorted(own)})", file=sys.stderr)
        return 2
    centre = own[a.pair]
    order = sorted(a.radii, key=lambda r: abs(r - centre))
    start = {}
    with pyazr.Session(str(azr), **opts) as s0:
        x = s0.parameters.read(a.parameter_file) if a.parameter_file else s0.parameters.x0
        start = dict(zip(s0.parameters.free.names, x))
    fitted = {}  # radius -> {name: value}
    rows, failed = [], 0
    for r in order:
        path, sav, outdir = variant_paths(azr, a.pair, r)
        p = pyazr.Project.open(azr)
        p.set_pair(a.pair, radius=r)
        p.set_config(output_dir=outdir)
        p.save(path)
        (azr.parent / outdir).mkdir(exist_ok=True)
        done = [q for q in fitted if (q - centre) * (r - centre) >= 0] or list(fitted)
        seed = fitted[min(done, key=lambda q: abs(q - r))] if done else start
        try:
            with pyazr.Session(str(path), **opts) as s:
                names = s.parameters.free.names
                x0 = np.array([seed.get(n, v) for n, v in zip(names, s.parameters.x0)], float)
                f0 = float(s.objective()(x0))
                if a.no_refit:
                    xf, ff, status = x0, f0, "not refitted"
                else:
                    xf, ff, status = refit(s, x0, a.fitter, a.max_nfev)
                chi2 = float(s.evaluate(xf).chi2)
                s.parameters.write(sav, xf)
                fitted[r] = dict(zip(names, xf))
            rows.append((r, f0, ff, chi2, len(names), status))
        except Exception as err:
            failed += 1
            rows.append((r, float("nan"), float("nan"), float("nan"), 0, f"failed: {err}"))
        if a.clean:
            for q in (path, sav):
                if q.exists():
                    q.unlink()
            shutil.rmtree(azr.parent / outdir, ignore_errors=True)
    rows.sort()
    lines = [f"# pair {a.pair}, project radius {centre:g} fm, fitter {'none' if a.no_refit else a.fitter}",
             f"# {'a_c_fm':>8s} {'objective_start':>16s} {'objective':>14s} {'data_chi2':>14s} {'free':>5s}  status"]
    lines += [f"  {r:8.4f} {f0:16.6f} {ff:14.6f} {c:14.6f} {n:5d}  {st}" for r, f0, ff, c, n, st in rows]
    good = [row for row in rows if np.isfinite(row[2])]
    if good:
        best = min(good, key=lambda t: t[2])
        lines.append(f"# minimum objective {best[2]:.6f} at a_c = {best[0]:g} fm; spread over the scan "
                     f"{max(t[2] for t in good) - best[2]:.6f}")
    print("\n".join(lines))
    if a.table:
        Path(a.table).write_text("\n".join(lines) + "\n")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
