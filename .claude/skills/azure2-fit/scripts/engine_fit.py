#!/usr/bin/env python3
"""The program's own fitters on a project, through `Session.fit_with_engine` (MIGRAD, MIGRAD with the analytic
gradient, Levenberg-Marquardt, the GSL trust region, NLopt where built), with MINOS, and their comparison.

Reads:  the project (format 3); with --sav a parameter file loaded by name first; with --x-json the start vector
        ("x", in the free parameters' order) of a JSON file such as lsq_fit.py writes.
Writes: --out JSON: per method its status, convergence, objective, EDM, calls, errors, MINOS intervals, and the end
        point; with --save-fit PATH the end point of the first method (Session.save_fit: PATH and its .sav).
        Nothing else: fit_with_engine writes no file.
Each method starts from the same point.  After a MIGRAD row the engine runs HESSE; the covariance after LM and
GSL-LM is (J^T J + priors)^-1.  The objective is the number the program prints as Total Chi-Squared.

    python3 engine_fit.py model.azr --method lm --method migrad-gradient --out fits.json
    python3 engine_fit.py model.azr --method migrad --minos --max-calls 5000

Exit status: 0 every method converged, 1 one did not, 2 bad arguments (an unknown method), 3 engine error.
"""
import argparse
import json
import sys

import numpy as np


def run(s, method, x0=None, minos=False, max_calls=None, log=print):
    """One fit_with_engine call; returns a plain dict."""
    r = s.fit_with_engine(method, minos=minos, x0=x0, max_calls=max_calls)
    out = {"method": method, "status": r.status, "converged": bool(r.converged), "objective": float(r.objective),
           "edm": float(r.edm), "calls": int(r.calls), "names": list(r.names), "x": np.asarray(r.x).tolist(),
           "errors": None if r.errors is None else np.asarray(r.errors).tolist(),
           "covariance": None if r.covariance is None else np.asarray(r.covariance).tolist(),
           "minos": {k: list(v) for k, v in (r.minos or {}).items()}}
    log(f"{method:16s} {r.status:12s} objective {r.objective:.10g}  calls {r.calls}  edm {r.edm:.3g}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("project")
    ap.add_argument("--method", action="append", default=[],
                    help="a fitter of pyazr.kinds()['fitters'] (repeatable; default migrad)")
    ap.add_argument("--minos", action="store_true", help="MINOS after each MIGRAD row")
    ap.add_argument("--max-calls", type=int, default=None)
    ap.add_argument("--sav", help="parameter file loaded by name first")
    ap.add_argument("--x-json", help="JSON file with the start vector under 'x'")
    ap.add_argument("--out", help="JSON result file")
    ap.add_argument("--save-fit", metavar="PATH", help="write the first method's end point as a project")
    a = ap.parse_args(argv)
    import pyazr
    known = list(pyazr.kinds()["fitters"])
    methods = a.method or ["migrad"]
    for m in methods:
        if m not in known:
            print(f"engine_fit: unknown method {m!r}; the engine's fitters: {', '.join(known)}", file=sys.stderr)
            return 2
    try:
        s = pyazr.Session(a.project, parameter_file=a.sav)
        x0 = None
        if a.x_json:
            with open(a.x_json) as fh:
                x0 = np.asarray(json.load(fh)["x"], float)
        rows = [run(s, m, x0=x0, minos=a.minos and m.startswith("migrad"), max_calls=a.max_calls)
                for m in methods]
    except (pyazr.Error, OSError, ValueError) as e:
        print(f"engine_fit: {e}", file=sys.stderr)
        return 3
    if a.out:
        with open(a.out, "w") as fh:
            json.dump({"project": a.project, "fits": rows}, fh, indent=1)
    if a.save_fit:
        s.save_fit(a.save_fit, np.asarray(rows[0]["x"]))
    s.close()
    return 0 if all(r["converged"] for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
