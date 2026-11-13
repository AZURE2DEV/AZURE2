"""Fit a project with scipy on the engine's objective, ``pyazr.Objective``.

Two recipes: ``scipy.optimize.minimize`` with the analytic gradient (L-BFGS-B, the parameters' bounds), and
``scipy.optimize.least_squares`` on the objective's rows with their Jacobian (``sum(f.residuals(x)**2) == f(x)``:
the data rows, then one row per prior and constraint, so both recipes minimise the same objective as the program's
fit).  The parameters of a project differ in scale by many orders (a reduced width in MeV^1/2, an areal density in
atoms/cm^2): least_squares takes ``x_scale="jac"``, and minimize works on x divided by a scale per entry, the
inverse length of the Jacobian's column at the start (``scale_of``; |x0| where the engine has no Jacobian).  On T2
from its start: least_squares in 8 evaluations, L-BFGS-B in about 35 (tests/pyazr/fit_recipes_test.py).

    python fit_scipy.py path/to/project.azr
"""
import sys

import numpy as np

import pyazr


def scale_of(f):
    """A scale per entry of x, so that a unit step moves the objective alike in every direction: 1/|J_k| of the
    Jacobian's column k at x0 (Gauss-Newton), or |x0| (1 where x0 is 0) where the engine has no Jacobian."""
    x0 = f.x0
    try:
        _, J = f.residuals_and_jacobian(x0)
        norm = np.sqrt(np.sum(J ** 2, axis=0))
        if np.all(norm > 0):
            return 1.0 / norm
    except pyazr.EngineError:
        pass
    return np.where(x0 != 0.0, np.abs(x0), 1.0)


def minimize(f, **options):
    """scipy.optimize.minimize (L-BFGS-B) with the gradient, on x / scale; returns (result, x)."""
    from scipy import optimize
    c = scale_of(f)

    def value_and_grad(y):
        value, grad = f.value_and_grad(y * c)
        return value, grad * c

    bounds = [(lo / ci, hi / ci) for (lo, hi), ci in zip(f.bounds, c)]
    opts = {"maxiter": 1000, "ftol": 1e-13, "gtol": 1e-8}
    opts.update(options)
    r = optimize.minimize(value_and_grad, f.x0 / c, jac=True, bounds=bounds, method="L-BFGS-B", options=opts)
    return r, r.x * c


def least_squares(f, **options):
    """scipy.optimize.least_squares on the objective's rows with their Jacobian; returns the result (x in r.x)."""
    from scipy import optimize
    lower = np.array([lo for lo, _ in f.bounds])
    upper = np.array([hi for _, hi in f.bounds])
    opts = {"x_scale": "jac", "xtol": 1e-12, "ftol": 1e-12, "gtol": 1e-12}
    opts.update(options)
    return optimize.least_squares(f.residuals, f.x0, jac=f.jacobian, bounds=(lower, upper), **opts)


def main(path):
    with pyazr.Session(path) as s:
        f = s.fit_objective()
        r, x = minimize(f)
        print(f"minimize (L-BFGS-B): objective {r.fun:.10g}, {r.nfev} evaluations, {r.message}")
        q = least_squares(f)
        print(f"least_squares: objective {f(q.x):.10g}, {q.nfev} evaluations, {q.njev} Jacobians")
        print(s.evaluate(x).table())
        return r, q


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "T2.azr")
