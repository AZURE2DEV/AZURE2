#!/usr/bin/env python3
"""Evaluate one AZURE2 project in the three level-parameter bases and compare what they compute.

A format-3 project holds observed level energies and partial widths (eV) or ANCs. The engine reads them into the basis
of the session: Brune (``use_brune=True``, the default), Park (``use_park=True``: Brune's observed amplitudes times
sqrt(J) of the level) or the standard Lane-Thomas basis (``use_brune=False``: formal energies and amplitudes with the
boundary condition at the lowest level of each J^pi). The three describe the same collision matrix, so every data set's
model must agree to rounding; a project with very broad background poles may have no standard equivalent (the session
then raises, naming the J^pi). This script opens a session per basis on the project's own values (a parameter file
carries its basis, so compare a fit after writing it into the project with ``Session.save_fit``) and prints the data
chi2 per basis and the largest relative difference of every data set's model from Brune's. It also prints each basis's
fit vector size and, for Park, the norms J of the levels.

Exit status: 0 when every basis agrees with Brune within ``--rtol`` (default 1e-6) or is refused by the engine as
having no equivalent (reported), 1 on a larger difference, 2 on a usage or set-up error.
"""
import argparse
import sys

import numpy as np

BASES = {"brune": dict(use_brune=True, use_park=False), "park": dict(use_brune=True, use_park=True),
         "standard": dict(use_brune=False, use_park=False)}


def evaluate(azr, options):
    """(data chi2, [model arrays], free-vector size, park norms or None) of the project's own values."""
    import pyazr
    with pyazr.Session(azr, **options) as s:
        x = s.parameters.x0
        ev = s.evaluate(x)
        norms = s.park_norms(x) if options.get("use_park") else None
        return float(ev.chi2), [np.asarray(d.model, float) for d in ev.data], [d.name for d in ev.data], x.size, norms


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--bases", nargs="+", default=list(BASES), choices=list(BASES))
    a = ap.parse_args(argv)
    import pyazr
    results = {}
    for name in a.bases:
        try:
            results[name] = evaluate(a.azr, BASES[name])
        except pyazr.Error as err:
            results[name] = err
    ref = results.get("brune")
    if ref is None or isinstance(ref, Exception):
        print(f"error: the Brune session failed: {ref}", file=sys.stderr)
        return 2
    status = 0
    for name, r in results.items():
        if isinstance(r, Exception):
            print(f"{name:9s} refused by the engine: {r}")
            continue
        chi2, models, names, nfree, norms = r
        worst, where = 0.0, ""
        for m0, m, label in zip(ref[1], models, names):
            scale = np.maximum(np.abs(m0), 1e-300)
            d = float(np.max(np.abs(m - m0) / scale)) if m0.size else 0.0
            if d > worst:
                worst, where = d, label
        print(f"{name:9s} data chi2 {chi2:.10g}  free {nfree}  largest rel. model difference from Brune {worst:.2e}"
              + (f" ({where})" if where else ""))
        if norms is not None:
            print("          Park norms J: " + " ".join(f"{j:.5f}" for j in norms))
        if worst > a.rtol:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
