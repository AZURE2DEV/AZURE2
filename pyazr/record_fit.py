"""Fit settings of a project's records that are not fit parameters of the engine, by refitting at each trial value.

A potential record's depth, radius, diffuseness, ``interior_boundary`` and ``interior_charge_e1`` are inputs of the
engine: its channel functions and capture integrals are made from them at set-up (docs/upgrade/interior_potential.md
section 10).  :func:`fit_record_settings` minimises over such settings from outside: for each trial value it writes
the value into a copy of the project, opens a session, runs the engine's fit of the project's own free parameters
and takes the objective there.  The result is the profile of the objective over the settings, its minimum, and on
request the interval of each setting over which the profile rises by ``delta`` (1: one standard deviation of one
setting).  It holds no formula; every number is the engine's.

    from pyazr.record_fit import Setting, fit_record_settings
    r = fit_record_settings("c12ag.azr", [
            Setting("a12C_p", "interior_charge_e1", "{:.8g}", start=-0.03, step=0.005),
            Setting("a12C_s", "central", "woods_saxon depth {:.8g} MeV radius 3.55 fm diffuseness 0.6 fm",
                    start=54.1, step=0.2)],
        intervals=True)
    r.values, r.objective, r.intervals

Each trial costs one set-up and one fit, so this suits projects whose fit takes seconds.  Needs scipy.
"""
import math
import os
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

__all__ = ["Setting", "RecordFitResult", "fit_record_settings"]


@dataclass
class Setting:
    """One number in one field of one record: ``text.format(value)`` is written as the field ``field`` of the record
    ``record`` (the whole field, so a structured field names its other values in ``text``).  ``start`` is the first
    trial, ``step`` the size of the first move (default 5 % of the start, or 0.01 at 0), ``lower`` and ``upper``
    optional limits (a trial outside them is not evaluated)."""
    record: str
    field: str
    text: str
    start: float
    step: Optional[float] = None
    lower: Optional[float] = None
    upper: Optional[float] = None

    def first_step(self):
        if self.step:
            return abs(float(self.step))
        return 0.05 * abs(self.start) if self.start else 0.01


@dataclass
class RecordFitResult:
    """What :func:`fit_record_settings` returns: ``values`` of the settings at the lowest objective found (their
    order), ``objective`` there (the engine's, after its fit of the project's free parameters), ``fit`` that fit's
    :class:`pyazr.session_fit.EngineFitResult` (``None`` for a project without free parameters), ``intervals`` one
    ``(lower, upper)`` per setting where the profile over the other settings rises by ``delta`` (``None`` on a side
    without a crossing inside the limits or within 60 steps; empty unless asked for), ``errors`` and ``covariance``
    of the settings from the Jacobian of the profiled least-squares rows at the minimum (linear errors, by central
    differences of a tenth of each step; ``None`` when the rows are fewer than the settings), ``evaluations`` every trial
    as ``(values, objective, outcome)`` in order, ``converged`` whether the outer minimiser met its tolerance,
    ``project`` the path of the project written with the values (``write=``), else ``None``."""
    values: np.ndarray
    objective: float
    fit: object
    converged: bool
    errors: Optional[np.ndarray] = None
    covariance: Optional[np.ndarray] = None
    intervals: list = field(default_factory=list)
    evaluations: list = field(default_factory=list)
    project: Optional[str] = None


def fit_record_settings(project, settings, *, method="migrad", session=None, starts=("project", "best"),
                        tolerance=1e-3, max_evaluations=400,
                        intervals=False, delta=1.0, write=None, progress=None):
    """Minimise the fitted objective of the project ``project`` (a path) over ``settings`` (:class:`Setting`).

    Parameters
    ----------
    method : the engine fitter of each trial (:meth:`Session.fit_with_engine`).
    session : keyword arguments of :class:`pyazr.Session` for every trial (``use_brune``, ``ignore_externals``, ...);
        the default is the session's own defaults except ``ignore_externals=False``: a setting of a potential
        usually acts through the computed capture.
    starts : where each trial's fit starts: ``"project"`` (the project's own values) and/or ``"best"`` (the fitted
        vector of the lowest trial so far).  With both, the default, each trial runs both fits and keeps the lower:
        a fit with several minima otherwise follows the path of the trials (seen on 12C(a,g): 143 from the warm
        start alone against 114 from the project's values).  ``("best",)`` halves the cost where the fit has one
        minimum.
    tolerance : the outer minimiser (Gauss-Newton on the objective's least-squares rows after each trial's fit,
        ``scipy.optimize.least_squares``) stops when the settings, in units of their steps, move by less than this.
    max_evaluations : the limit of Gauss-Newton steps of the minimisation (each takes one trial and two more per
        setting for the Jacobian).
    intervals : also find, for each setting, where the profile over the others rises by ``delta``.
    write : a path; the project with the best values is written there (the project's own free parameters keep their
        written values: the fit's are in ``fit``).
    progress : called with ``(values, objective, outcome)`` after each trial.

    A trial whose value the engine refuses (a diagnostic of the record) or whose fit does not end at a minimum
    counts as +infinity, with the reason as its outcome.  The project file is not changed; trials are written beside
    it (``.<name>.recordfit.<pid>.azr``, removed at the end) so that its relative paths hold.
    """
    from scipy import optimize  # not a dependency of pyazr itself

    from . import Session
    from .errors import Error
    from .project import Project

    settings = list(settings)
    starts = tuple(starts)
    if not starts or any(k not in ("project", "best") for k in starts):
        raise ValueError('fit_record_settings: starts is "project" and/or "best".')
    if not settings:
        raise ValueError("fit_record_settings: no settings.")
    path = os.path.abspath(project)
    directory = os.path.dirname(path)
    trial_path = os.path.join(directory, f".{os.path.basename(path)}.recordfit.{os.getpid()}.azr")
    options = {"ignore_externals": False}
    options.update(session or {})
    steps = np.array([s.first_step() for s in settings])
    evaluations = []
    best = {"objective": math.inf, "values": None, "fit": None, "x": None}

    def apply(values, target):
        edited = Project.open(path)
        for s, v in zip(settings, values):
            edited.set(s.record, **{s.field: s.text.format(float(v))})
        edited.save(target)

    def trial(values):
        """The fitted objective at ``values`` and its least-squares rows there (``None`` for a failed trial)."""
        values = np.atleast_1d(np.asarray(values, float))
        outcome, objective, fit, rows = "ok", math.inf, None, None
        if any((s.lower is not None and v < s.lower) or (s.upper is not None and v > s.upper)
               for s, v in zip(settings, values)):
            outcome = "outside the limits"
        else:
            try:
                apply(values, trial_path)
                S = Session(trial_path, cwd=directory, **options)
                x = np.asarray(S.parameters.x0, float)
                if len(S.parameters.free.names):
                    warm = best["x"] if best["x"] is not None and len(best["x"]) == len(S.parameters.free.names) else None
                    vectors = [None] if "project" in starts or warm is None else []
                    if "best" in starts and warm is not None:
                        vectors.append(warm)
                    ended = []
                    for x0 in vectors:
                        candidate = S.fit_with_engine(method, x0=x0)
                        ended.append(candidate.outcome)
                        if candidate.outcome == "converged" and math.isfinite(candidate.objective) and candidate.objective < objective:
                            objective, fit, x = candidate.objective, candidate, candidate.x
                    if fit is None:
                        outcome = "fit " + ", ".join(ended)
                else:
                    objective = float(S.evaluate().objective)
                    if not math.isfinite(objective):
                        objective, outcome = math.inf, "objective not finite"
                if math.isfinite(objective):
                    rows = S.objective().residuals(x)
            except Error as e:
                outcome = f"refused: {str(e).splitlines()[0][:200]}"
        evaluations.append((values.copy(), objective, outcome))
        if objective < best["objective"]:
            best.update(objective=objective, values=values.copy(), fit=fit, x=None if fit is None else fit.x)
        if progress is not None:
            progress(values.copy(), objective, outcome)
        return objective, rows

    state = {"rows": None}

    def minimise(free, fixed_values, limit):
        """The minimum over the settings ``free`` (indices) with the others at ``fixed_values``: (values, objective,
        converged, the least_squares result or None).  The unknowns are (value - first value)/step."""
        fixed_values = np.array(fixed_values, float)

        def full(u):
            v = fixed_values.copy()
            v[free] = fixed_values[free] + np.atleast_1d(u) * steps[free]
            return v

        if len(free) == 0:
            objective, _ = trial(fixed_values)
            return fixed_values, objective, True, None

        def rows_of(u):
            objective, rows = trial(full(u))
            if rows is None:  # a failed trial: far worse than any other, the same length as the others
                n = len(state["rows"]) if state["rows"] is not None else max(len(free), 1)
                return np.full(n, 1e8)
            state["rows"] = rows
            return rows

        lower = np.array([-np.inf if settings[i].lower is None else (settings[i].lower - fixed_values[i]) / steps[i] for i in free])
        upper = np.array([np.inf if settings[i].upper is None else (settings[i].upper - fixed_values[i]) / steps[i] for i in free])
        first = rows_of(np.zeros(len(free)))
        if state["rows"] is None:
            return fixed_values, math.inf, False, None
        if len(first) < len(free):  # fewer rows than unknowns: no Gauss-Newton; the simplex on the objective
            r = optimize.minimize(lambda u: trial(full(u))[0], np.zeros(len(free)), method="Nelder-Mead",
                                  options={"xatol": tolerance, "fatol": tolerance, "maxfev": limit})
            return full(r.x), float(r.fun), bool(r.success), None
        def jacobian(u):
            # Central differences of a tenth of each setting's step (one-sided at a limit).  scipy's own differences
            # take a step of 1e-8 where the unknown is 0, which it is at the start: far below the scatter the
            # trials' fits leave in the rows.
            u = np.atleast_1d(np.asarray(u, float))
            columns = []
            for k in range(len(free)):
                up, down = u.copy(), u.copy()
                up[k] = min(u[k] + 0.1, upper[k])
                down[k] = max(u[k] - 0.1, lower[k])
                columns.append((rows_of(up) - rows_of(down)) / (up[k] - down[k]) if up[k] > down[k] else np.zeros(len(first)))
            return np.column_stack(columns)

        r = optimize.least_squares(rows_of, np.zeros(len(free)), jac=jacobian, bounds=(lower, upper), xtol=tolerance,
                                   ftol=1e-12, gtol=1e-12, max_nfev=limit)
        return full(r.x), float(2.0 * r.cost), bool(r.status > 0), r

    try:
        start = np.array([s.start for s in settings], float)
        _, _, converged, solution = minimise(list(range(len(settings))), start, max_evaluations)
        if best["values"] is None:
            reasons = sorted({e[2] for e in evaluations})
            raise RuntimeError("fit_record_settings: no trial gave a finite objective (" + "; ".join(reasons) + ").")
        errors = covariance = None
        if solution is not None and solution.jac.shape[0] >= len(settings):
            jtj = solution.jac.T @ solution.jac
            if np.linalg.matrix_rank(jtj) == len(settings):
                covariance = np.linalg.inv(jtj) * np.outer(steps, steps)
                errors = np.sqrt(np.diag(covariance))
        found = []
        if intervals:
            floor = best["objective"]
            centre = best["values"].copy()
            others = lambda i: [k for k in range(len(settings)) if k != i]

            def rise(i, v):
                fixed = centre.copy()
                fixed[i] = v
                return minimise(others(i), fixed, 200)[1] - floor - delta

            for i, s in enumerate(settings):
                sides = []
                for sign in (-1.0, 1.0):
                    inside, crossing = centre[i], None
                    for k in range(1, 61):
                        outside = centre[i] + sign * steps[i] * 0.5 * k
                        if (s.lower is not None and outside < s.lower) or (s.upper is not None and outside > s.upper):
                            break
                        if rise(i, outside) > 0.0:
                            crossing = optimize.brentq(lambda v: rise(i, v), inside, outside,
                                                       xtol=1e-3 * tolerance * steps[i], maxiter=60)
                            break
                        inside = outside
                    sides.append(None if crossing is None else float(crossing))
                found.append(tuple(sides))
        written = None
        if write is not None:
            written = os.path.abspath(write)
            apply(best["values"], written)
        return RecordFitResult(values=best["values"], objective=best["objective"], fit=best["fit"], converged=converged,
                               errors=errors, covariance=covariance, intervals=found, evaluations=evaluations, project=written)
    finally:
        if os.path.exists(trial_path):
            os.remove(trial_path)
