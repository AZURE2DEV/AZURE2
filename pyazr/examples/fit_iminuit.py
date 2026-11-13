"""Fit a project with iminuit (MIGRAD, then HESSE) on the engine's objective, ``pyazr.Objective``.

The objective is what the program's own fit minimises (the data chi-squared and every prior and constraint), with
its analytic gradient; ``errordef`` is 1 (a chi-squared), which Minuit reads from the object.  On T2 MIGRAD reaches
the program's minimum (``tests/pyazr/fit_recipes_test.py``).

    python fit_iminuit.py path/to/project.azr
"""
import sys

import numpy as np

import pyazr


def fit(f, hesse=True):
    """MIGRAD (and HESSE) on the Objective ``f``; returns the ``iminuit.Minuit`` object."""
    import iminuit
    m = iminuit.Minuit(f, f.x0, grad=f.grad, name=f.names)
    m.limits = f.bounds
    m.migrad()
    if hesse:
        m.hesse()
    return m


def main(path):
    with pyazr.Session(path) as s:
        f = s.fit_objective()
        m = fit(f)
        x = np.array(m.values)
        print(f"MIGRAD: objective {m.fval:.10g}, valid {m.valid}, {m.nfcn} calls, {m.ngrad} gradients")
        for name, value, error in zip(f.names, x, m.errors):
            print(f"  {name:<24} {value:14.7g} +- {error:.3g}")
        ev = s.evaluate(x)
        print(ev.table())
        return m


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "T2.azr")
