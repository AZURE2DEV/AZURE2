"""The objective of a session as one callable for any minimiser or sampler (docs/upgrade/pyazr3.md section 7).

:class:`Objective` is what the engine's fit minimises: the data chi-squared and every prior and constraint
(docs/upgrade/objective.md), on the free vector x of the basis the session fits in.  Every number comes from the
engine; this module holds no formula.  Recipes (``pyazr/examples/fit_scipy.py``, ``fit_iminuit.py``,
``fit_emcee.py``)::

    f = s.fit_objective()
    r = scipy.optimize.minimize(f.value_and_grad, f.x0, jac=True, bounds=f.bounds, method="L-BFGS-B")
    r = scipy.optimize.least_squares(f.residuals, f.x0, jac=f.jacobian)
    m = iminuit.Minuit(f, f.x0, grad=f.grad, name=f.names); m.limits = f.bounds; m.migrad(); m.hesse()
    sampler = emcee.EnsembleSampler(nwalkers, len(f.x0), f.log_prob)
"""
import numpy as np

from .errors import EngineError
from .locking import session_lock


def require_data_mode(session, what):
    """Refuse a session in extrapolation mode: it has no data and no objective."""
    if session.mode != "data":
        raise EngineError(f"{what} needs the session's data (it is in extrapolation mode, which has no data and no "
                          "objective); call data_mode() first.")


class Objective:
    """The engine's fit objective as a callable: chi-squared of the data plus every prior and constraint.

    Built by ``Session.fit_objective()``.  Everything is on the free vector x (``x0``, ``names``, ``bounds``, in
    ``parameters.free`` order), in the basis the session fits in (Brune, Park or formal amplitudes; norms, shifts and
    the other nuisances as they are).

    - ``f(x)``: the objective, ``Session.objective(x)`` bit for bit; the same terms and order as the program's fit.
    - ``f.grad(x)``, ``f.value_and_grad(x)``: its gradient (analytic, finite differences for energy shifts).
    - ``f.residuals(x)``, ``f.jacobian(x)``, ``f.residuals_and_jacobian(x)``: the objective as least-squares rows,
      ``sum(r**2) == f(x)`` (to rounding): the data rows (whitened where a data covariance joins points), the model
      data, the prior covariance, Park's wall, the ``gamma_decay`` priors, the abundance sums, then one row
      ``(p - centre)/sigma`` per prior of a nuisance parameter.  The first rows are ``Session.residuals(x)`` and the
      rows of ``Session.residual_jacobian(x)``, bit for bit.
    - ``f.log_prob(x)``: ``-0.5 * f(x)`` (no constant), ``-inf`` where the engine's ensemble sampler puts zero
      probability: outside ``bounds``, the fit's limits of a free parameter of any kind, and where Park's overlap J
      of a level is negative.
    - ``x0`` (a copy of the start), ``names``, ``bounds`` (``(lower, upper)``, ``-inf``/``inf`` where the fit has
      none: the limits a parameter is created with, then those of a ``<parameterSettings>`` block, as the
      program's fit and sampler have them), ``errordef`` (1: a chi-squared; for iminuit).

    Arrays returned are copies the caller owns.  The object holds its session and these copies, no cache.  Calls from
    several threads on one session are safe: they take the session's lock and run one after the other.  A session
    cannot be pickled: a process pool opens one session per worker (``pyazr/examples/fit_emcee.py``).
    """

    errordef = 1.0
    BUILT_BY = "Session.fit_objective()"
    FIELDS = {
        "__call__(x)": "the objective at x: data chi-squared plus every prior and constraint",
        "grad(x)": "its gradient over the free vector",
        "value_and_grad(x)": "(value, gradient) from one engine call",
        "residuals(x)": "the objective as least-squares rows, sum(r**2) == f(x)",
        "jacobian(x)": "the rows' Jacobian, shape (rows, free parameters)",
        "residuals_and_jacobian(x)": "(rows, Jacobian) from one engine call",
        "log_prob(x)": "-0.5 f(x), -inf outside the engine's sampler's support (outside bounds, Park's J < 0)",
        "x0": "a copy of the start of the free vector",
        "names": "the free parameters' names",
        "bounds": "(lower, upper) per free parameter, the fit's limits, -inf / inf where none; walls of log_prob",
        "errordef": "1 (a chi-squared)",
    }

    def __init__(self, session):
        require_data_mode(session, "fit_objective()")
        free = session.parameters.free
        self._session = session
        self._lock = session_lock(session)
        self._x0 = np.array(session.params_rwa, float).ravel()
        self.names = list(free.names)
        self.bounds = [(float(lo), float(hi)) for lo, hi in free.bounds]

    @property
    def x0(self):
        """A copy of the start of the free vector (the session's ``params_rwa`` when the object was made)."""
        return self._x0.copy()

    def _engine(self):
        return self._session.sess  # SessionClosed on a closed session

    @staticmethod
    def _vector(x):
        return np.asarray(x, float).ravel()

    def __call__(self, x):
        """The objective at x (a float): chi-squared of the data plus every prior and constraint."""
        with self._lock:
            return float(self._engine().calculate_objective_rwa(self._vector(x)))

    def value_and_grad(self, x):
        """``(f(x), gradient)`` from one engine call (for ``scipy.optimize.minimize(..., jac=True)``)."""
        with self._lock:
            resp = np.asarray(self._engine().calculate_objective_grad_rwa(self._vector(x)), float)
        return float(resp[0]), resp[1:].copy()

    def grad(self, x):
        """The gradient of the objective at x, one entry per free parameter (per unit of the parameter)."""
        return self.value_and_grad(x)[1]

    def residuals(self, x):
        """The objective's least-squares rows at x from a forward pass (``sum(r**2) == f(x)``)."""
        with self._lock:
            out = self._engine().objective_rows_rwa(self._vector(x), False)
        return np.asarray(out["rows"], float)

    def residuals_and_jacobian(self, x):
        """``(rows, J)``: the rows of :meth:`residuals` and their Jacobian, shape (rows, free parameters)."""
        with self._lock:
            out = self._engine().objective_rows_rwa(self._vector(x), True)
        if not out["ok"]:
            raise EngineError("the engine has no analytic Jacobian for a data set of this project (residual_jacobian "
                              "returned none); use residuals() with a finite-difference Jacobian instead.")
        return np.asarray(out["rows"], float), np.asarray(out["jacobian"], float)

    def jacobian(self, x):
        """The Jacobian of :meth:`residuals` at x, shape (rows, free parameters) (``least_squares(jac=...)``)."""
        return self.residuals_and_jacobian(x)[1]

    def log_prob(self, x):
        """ln P = -0.5 * f(x), or ``-inf`` outside the engine's sampler's support (see the class docstring)."""
        with self._lock:
            return float(self._engine().log_prob_rwa(self._vector(x)))

    def describe(self):
        """Plain data: the free parameters (name, bounds, start), errordef, the terms' names."""
        rows = [{"name": n, "bounds": list(b), "start": float(v)}
                for n, b, v in zip(self.names, self.bounds, self._x0)]
        return {"parameters": rows, "errordef": self.errordef, "terms": list(self._engine().objective_term_names())}

    def __repr__(self):
        return f"Objective({len(self._x0)} free parameters of {self._session.file!r})"
