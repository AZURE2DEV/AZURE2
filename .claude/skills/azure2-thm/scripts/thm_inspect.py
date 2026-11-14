#!/usr/bin/env python3
"""What the THM experiments of an AZURE2 project do at a parameter vector: profiled norm, background, vertex reach.

Reads a format-3 project with a ``<thm>`` block and, optionally, a parameter file (by name). For every THM experiment
(``experiment[<name>]`` of ``<thm>``) prints the data sets it groups, its points and chi2, the profiled norm n* with its
1 sigma (from the profile at fixed R-matrix parameters, not scaled by chi2/nu), the background coefficients b0, b1, b2
(model units, E the c.m. energy of the THM entrance pair in MeV) and the profile's status; then the entrance vertex
over the experiment's data energies (``Session.thm_vertex``): the window or model, and how many energies the window or
the DW vertex grid does not reach. With ``--distortion`` also the distortion factor R(E) at those energies when the
experiment has ``distortion=``. Equations: docs/source/theory/thm_implementation.rst.

Flags (lines ``FLAG ...``): a profile status other than "profiled", a vertex that does not reach every data energy, a
norm whose relative error exceeds ``--max-norm-error`` (default 0.5: the data do not fix the scale).
Exit status: 0 without a flag, 1 with one, 2 on a usage error or a project without THM experiments.
"""
import argparse
import sys
import warnings

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


def inspect(session, x, distortion=False):
    """{experiment: report} with the data energies, the data sets and the vertex's reach added."""
    ev = session.evaluate(x)
    by_segment = {d.segment: d for d in ev.data}
    out = {}
    for name, rep in session.thm_experiments(x).items():
        sets = [by_segment[k] for k in rep["segments"] if k in by_segment]
        energies = np.unique(np.concatenate([np.asarray(d.energies, float) for d in sets])) if sets else np.array([])
        rep = dict(rep, datasets=[d.name for d in sets], energies=energies,
                   dataset_norms=[float(d.norm) for d in sets], dataset_chi2=sum(float(d.chi2) for d in sets))
        if energies.size:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                v = session.thm_vertex(name, energies, x)
            rep["vertex_window"] = v.get("window")
            rep["vertex_model"] = v.get("model")
            rep["vertex_unreached"] = int(np.sum(~np.asarray(v.get("reached", np.ones(energies.size)), bool)))
            if distortion:
                try:
                    rep["distortion"] = np.asarray(session.thm_distortion(name, energies)["R"], float)
                except (KeyError, ValueError) as err:
                    rep["distortion"] = str(err)
        out[name] = rep
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr")
    ap.add_argument("--params", help="parameter file (param.sav), by name; default: the project's values")
    ap.add_argument("--option", action="append", default=[], metavar="KEY=VALUE", help="a session option")
    ap.add_argument("--distortion", action="store_true", help="also R(E) of experiments with distortion=")
    ap.add_argument("--max-norm-error", type=float, default=0.5)
    a = ap.parse_args(argv)
    import pyazr
    opts = {k: option_value(v) for k, v in (o.split("=", 1) for o in a.option)}
    with pyazr.Session(a.azr, **opts) as s:
        x = s.parameters.read(a.params) if a.params else s.parameters.x0
        reports = inspect(s, x, a.distortion)
    if not reports:
        print("error: the project has no THM experiment (<thm> experiment[...] lines)", file=sys.stderr)
        return 2
    flags = []
    for name, r in reports.items():
        print(f"experiment {name}: data sets {', '.join(r['datasets'])}; {r['points']} points, chi2 {r['chi2']:.4f}")
        print(f"  norm n* = {r['norm']:.6g} +- {r['sigma'][0]:.3g}; background {r['background']} "
              f"{' '.join(f'{b:.4g}' for b in r['b'])}; status {r['status']}")
        if r.get("cbkg"):
            print("  coherent background: " + ", ".join(f"{k} {v:.4g}" for k, v in r["cbkg"].items()))
        if "vertex_model" in r:
            print(f"  vertex {r['vertex_model']} {r['vertex_window']}: {r['vertex_unreached']} of {r['energies'].size}"
                  f" data energies not reached")
        if isinstance(r.get("distortion"), np.ndarray):
            print(f"  distortion R(E) from {r['distortion'].min():.4g} to {r['distortion'].max():.4g}")
        if r["status"] != "profiled":
            flags.append(f"{name}: profile status {r['status']}")
        if r.get("vertex_unreached"):
            flags.append(f"{name}: the vertex does not reach {r['vertex_unreached']} data energies")
        if r["norm"] and r["sigma"][0] / abs(r["norm"]) > a.max_norm_error:
            flags.append(f"{name}: n* = {r['norm']:.4g} +- {r['sigma'][0]:.3g}: the data do not fix the scale")
    for f in flags:
        print("FLAG", f)
    return 1 if flags else 0


if __name__ == "__main__":
    sys.exit(main())
