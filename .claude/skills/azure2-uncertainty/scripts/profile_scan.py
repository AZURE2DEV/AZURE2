#!/usr/bin/env python3
"""Profile likelihood of one free parameter of an AZURE2 fit: the parameter held at each value of a grid, every
other free parameter refitted (least squares on the objective's rows and analytic Jacobian, warm-started from the
neighbour), and the interval where Delta objective <= delta (1: 68 % for one parameter), against the parabolic
error and, on request, the engine's MINOS.

Reads:  the project; the best fit from --x-json ("x") or --sav, else the project's start (refitted first).
Writes: --out JSON: the grid, the profile, the interval, the parabolic error, MINOS's interval when asked;
        a table on standard output.

    python3 profile_scan.py fit.azr --x-json fit.json --parameter energy_1 --points 9 --span 3 --minos

The grid spans +-span parabolic errors (sqrt of the diagonal of (J^T J)^-1) around the best fit.  A profile that
is not parabolic, or an interval that touches a limit (a width at 0, theta^2 = 1), is the reason to quote the
profile instead of the parabolic error.
Exit status: 0 interval found on both sides, 1 open on a side (widen --span), 2 bad arguments, 3 engine error.
"""
import argparse
import json
import sys

import numpy as np
from scipy.optimize import least_squares


def refit(f, x, held, max_nfev=200):
    """Minimise the objective over the free columns other than `held` (held at x); returns the end point."""
    cols = [i for i in range(len(x)) if i not in held]
    lo = np.array([f.bounds[i][0] for i in cols])
    hi = np.array([f.bounds[i][1] for i in cols])
    base = np.array(x, float)

    def full(z):
        y = base.copy()
        y[cols] = z
        return y

    def fun(z):
        return f.residuals(full(z))

    def jac(z):
        return f.jacobian(full(z), skip=list(held))[:, cols]

    span = np.where(np.isfinite(hi - lo), hi - lo, 1.0)
    z0 = np.clip(base[cols], lo + 1e-9 * span, hi - 1e-9 * span)
    r = least_squares(fun, z0, jac=jac, bounds=(lo, hi), method="trf", x_scale="jac", xtol=1e-10, ftol=1e-10,
                      gtol=1e-10, max_nfev=max_nfev)
    return full(r.x)


def crossing(values, deltas, centre, level):
    """Where the profile crosses `level` on each side of `centre`, interpolated linearly in sqrt(Delta) (exact for a
    parabolic profile); None where it does not cross."""
    out = []
    for side in (-1, 1):
        pts = sorted([(v, np.sqrt(max(d, 0.0))) for v, d in zip(values, deltas) if side * (v - centre) >= 0],
                     key=lambda t: abs(t[0] - centre))
        if not pts or abs(pts[0][0] - centre) > 0:
            pts = [(centre, 0.0)] + pts
        hit, q = None, np.sqrt(level)
        for (v1, d1), (v2, d2) in zip(pts, pts[1:]):
            if d1 <= q < d2:
                hit = v1 + (q - d1) * (v2 - v1) / (d2 - d1)
                break
        out.append(hit)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("project")
    ap.add_argument("--parameter", required=True, help="a free parameter's name (Objective.names)")
    ap.add_argument("--x-json")
    ap.add_argument("--sav")
    ap.add_argument("--points", type=int, default=9)
    ap.add_argument("--span", type=float, default=3.0, help="half-width of the grid in parabolic errors")
    ap.add_argument("--delta", type=float, default=1.0, help="the Delta objective of the interval")
    ap.add_argument("--minos", action="store_true", help="also MINOS (fit_with_engine('migrad', minos=True))")
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    import pyazr
    try:
        s = pyazr.Session(a.project, parameter_file=a.sav)
        f = s.objective()
        if a.parameter not in f.names:
            print(f"profile_scan: {a.parameter!r} is not free; the free parameters: {f.names}", file=sys.stderr)
            return 2
        j = f.names.index(a.parameter)
        x = f.x0
        if a.x_json:
            with open(a.x_json) as fh:
                x = np.asarray(json.load(fh)["x"], float)
        x = refit(f, x, held=[])
        best = f(x)
        _, J = f.residuals_and_jacobian(x)
        sigma = float(np.sqrt(np.linalg.inv(J.T @ J)[j, j]))
        grid = x[j] + sigma * np.linspace(-a.span, a.span, a.points)
        prof = {}
        for side in (np.where(grid >= x[j])[0], np.where(grid < x[j])[0][::-1]):
            y = x.copy()
            for k in side:
                y[j] = grid[k]
                y = refit(f, y, held=[j])
                prof[k] = f(y) - best
    except (pyazr.Error, np.linalg.LinAlgError) as e:
        print(f"profile_scan: {e}", file=sys.stderr)
        return 3
    deltas = [prof[k] for k in range(len(grid))]
    for v, d in zip(grid, deltas):
        print(f"{a.parameter} {v:14.8g}  Delta objective {d:10.4f}")
    lo, hi = crossing(grid, deltas, x[j], a.delta)
    print(f"best {x[j]:.8g}; parabolic error {sigma:.4g}; profile interval "
          f"{'open' if lo is None else f'{lo - x[j]:+.4g}'} / {'open' if hi is None else f'{hi - x[j]:+.4g}'}")
    result = {"parameter": a.parameter, "best": float(x[j]), "objective": best, "sigma_parabolic": sigma,
              "grid": grid.tolist(), "delta": deltas, "interval": [lo, hi]}
    if a.minos:
        r = s.fit_with_engine("migrad", minos=True, x0=x)
        result["minos"] = r.minos.get(a.parameter)
        print(f"MINOS ({r.status}, objective {r.objective:.8g}): {result['minos']}")
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(result, fh, indent=1)
    s.close()
    return 0 if lo is not None and hi is not None else 1


if __name__ == "__main__":
    sys.exit(main())
