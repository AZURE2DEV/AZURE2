"""The adapter of the web UI to the pyazr session: the only file of ``pyazr/ui`` that names the class
:class:`pyazr.Session` and its methods (docs/upgrade/web_ui.md section 11).

It runs inside a worker process (:mod:`pyazr.ui.workers`), never in the server: a session build holds the GIL for its
whole length (12.8 s on 13N, web_ui.md section 2.1).  Every answer is plain JSON data (lists, dicts, floats); a
parameter vector ``x`` is the engine's free vector in its order and basis (the fit vector, pyazr3.md decision 10).

The chi-squared per data set is the program's: the rows of the ``chiSquared.out`` the session writes into the UI's
own output directory (``<config>``'s output directory of the session's file points there, see
:meth:`pyazr.ui.state.State.session_text`), which is ``Evaluation.data[i].chi2``; the sum of the squares of a data
set's rows of the objective (``Evaluation.data[i].rows``) splits data sets joined by a data covariance otherwise
(web_ui.md section 2.2, finding 5).  The file has the program's six significant digits.
"""
import contextlib
import dataclasses
import inspect
import os
import time

import numpy as np

from pyazr import errors, session_checks
from pyazr import parameters as pyazr_parameters
from pyazr.locking import session_lock

# The methods of pyazr.Session the UI does not call, by reason (the route-table consistency test,
# tests/pyazr/ui/routes_consistency_test.py, requires every public method to be called here or listed here).
SESSION_NOT_SERVED = {
    "the sensitivities of the model values: the band jobs use them through pyazr.bands.uncertainty_bands":
        "model_gradients",
    "dimensionless widths and level quantities, not a route of this step":
        "wigner_widths dimensionless_widths park_norms scattering_lengths",
    "session edits: edits of the Project and a new session (pyazr3.md decision 11)":
        "configure recalculate_external_capture "
        "refreeze_covariance nuclear_potential nuclear_potentials",
    "the server owns sessions: lifecycle": "close is_alive",
    "later steps of web_ui.md section 11 (the Nucleus and Physics tabs)":
        "level_scheme nuclei mass_number basis ec_integrals prior_covariance_residuals prior_covariance_names "
        "datasets cache_stats output_dir",
}



# The residual the payload carries for every point of a data set: the engine's standardized residual, the one whose
# square is its chi-squared term (the data rows of pyazr's ``Objective.residuals``).
RESIDUAL = {
    "definition": "(model - n*datum) / (n*error)",
    "doc": "the engine's standardized residual of each point: model minus the datum times the data set's "
           "normalization n, over the error times n, in the data set's own frame and unit (energy shifts applied). "
           "Its squares add up to the data set's chi-squared of chiSquared.out (to the file's six digits) unless a "
           "data covariance joins the data set with others: then the residuals are whitened rows w = W r "
           "(W^T W = M^-1), whose squares add up to the group's chi-squared but split it differently from the "
           "program's sum r_i (M^-1 r)_i per data set (web_ui.md finding 5); residual_chi2 is the sum of their "
           "squares (the data set's rows of pyazr's Objective.residuals, Evaluation.data[i].rows).",
}
# The dimension (pyazr.kinds()["units"]) of a data set's engine unit for which the S factor is offered, with a
# charged entrance pair (the S factor's own rule in the engine's quantity table); the unit of the converted values.
S_FACTOR_DIMENSION = "cross_section"
S_FACTOR_UNIT = "MeV*b"

def _plain(a):
    return np.asarray(a, float).tolist()


class Engine:
    """One :class:`pyazr.Session` and what the routes ask of it."""

    def __init__(self):
        self.s = None
        self._extrapolation_session = None
        self._curves = None  # the session of the model curves (pyazr.model_curve.CurveSession) of the open session
        self._curves_at = (None, None)  # (x as bytes, {name: ModelCurve or refusal}) of the last curves evaluated

    # -- the session ----------------------------------------------------------------------------------------------

    def open(self, file, cwd):
        """Build the session of ``file`` with relative paths resolved in ``cwd`` (the project's directory)."""
        import pyazr
        if self.s is not None:
            self.s.close()
            self.s = None
        if self._extrapolation_session is not None:
            self._extrapolation_session.close()
            self._extrapolation_session = None
        self._close_curves()
        s = pyazr.Session(file, cwd=cwd)
        os.makedirs(s.output_dir, exist_ok=True)  # the external-capture cache and chiSquared.out are written there
        self.s = s
        self.info = [_dataset_info(d) for d in s.active_datasets]  # in the engine's order (a copy per access)
        self.names = [d["name"] for d in self.info]
        kinds = pyazr.kinds()
        self.units, self.data_kinds = kinds["units"], kinds["data"]
        self.charged = {p.key: p.Z1 * p.Z2 > 0 for p in s.pairs}  # pair key -> both constituents charged
        return self._describe()

    def _describe(self):
        """The session: its options and where each came from, free-vector length, data sets (in the engine's order) and parameters."""
        s = self.s
        return {"options": dict(s.options), "option_sources": dict(s.option_sources), "n_free": len(s.parameters.x0),
                "datasets": self._datasets(), "pairs": s.describe()["pairs"],
                "parameters": self._parameters()}

    def _datasets(self):
        return [dict(d, points=len(self.s.energies[i])) for i, d in enumerate(self.info)]

    def _parameters(self):
        """One row per parameter of the engine (fixed and free), the session's own description (``describe()``): name,
        kind, unit, whether it is free and its index in x, ``x`` (fit basis), ``physical_value``, bounds, owner
        fields; ``x0`` is ``x`` of a free parameter."""
        common = ("name", "kind", "unit", "start", "value", "fixed", "free_index", "bounds")
        rows = []
        for p in self.s.describe()["parameters"]:  # the session's records
            row = {"name": p["name"], "kind": p["kind"], "unit": p["unit"], "x": p["start"],
                   "physical_value": p["value"], "free": not p["fixed"], "free_index": p["free_index"],
                   "bounds": p["bounds"]}
            row.update((k, v) for k, v in p.items() if k not in common and v is not None)
            row["x0"] = row["x"] if row["free"] else None
            rows.append(row)
        return rows

    def _free_names(self):
        """Names of the free parameters, in the order of x."""
        return [p.name for p in self.s.parameters if p.free_index is not None and not p.fixed]

    # -- evaluations ----------------------------------------------------------------------------------------------

    def evaluate(self, x):
        """The objective with its named terms, the data chi-squared (with the model's own terms: Park's wall, the
        abundance priors and sums, the measured model data; 0 in most projects), and the program's chi-squared per
        data set (``chiSquared.out``: chi2, N, norm, norm chi2 per row; its Total-
        line as ``totals``)."""
        return self._evaluate(np.asarray(x, float))[0]

    def _evaluate(self, x):
        """``(the answer of evaluate, the Evaluation at x)``: the engine's data chi-squared and the program's output
        files at x (through the binding, under the session's lock: the number the page has always shown, which the sum
        of ``Evaluation.chi2`` and the model's terms gives only to rounding where a data covariance joins data sets),
        then one evaluation."""
        s = self.s
        with session_lock(s):
            chi2 = float(s._engine.calculate_chi2_rwa(x))
            os.makedirs(s.output_dir, exist_ok=True)
            if not s._engine.write_output_files():
                raise RuntimeError("the engine wrote no output files (no data object)")
            session_checks.warn_engine_notices(s._engine, stacklevel=2)  # the engine's notices as warnings
        ev = s.evaluate(x)
        rows, totals = _read_chi_squared(os.path.join(s.output_dir, "chiSquared.out"))
        if len(rows) != len(self.names):
            raise RuntimeError(f"chiSquared.out has {len(rows)} rows, the session {len(self.names)} data sets")
        per_set = [dict(name=n, **r) for n, r in zip(self.names, rows)]
        return {"objective": ev.objective, "chi2": chi2, "terms": dict(ev.terms), "datasets": per_set,
                "totals": totals}, ev

    def data(self, name, x):
        """The arrays of the data set ``name`` at x (:meth:`evaluation`'s entry for it, without the chi-squared
        table)."""
        i = self._index(name)
        x = np.asarray(x, float)
        return self._dataset_arrays(x, self.s.evaluate(x))[i]

    def evaluation(self, x):
        """The Results payload at x: the program's chi-squared split (:meth:`evaluate`, which also writes the program's
        output files into the session's output directory), every data set's arrays (its chi-squared row merged in),
        every extrapolation's arrays and the model data.  Every number is the engine's or pyazr's (see
        :meth:`_dataset_arrays` for the two products the adapter forms)."""
        x = np.asarray(x, float)
        out, ev = self._evaluate(x)
        arrays = self._dataset_arrays(x, ev)
        out["datasets"] = [dict(a, **row) for a, row in zip(arrays, out["datasets"])]
        out["extrapolations"] = self._extrapolation_arrays(x)
        from pyazr import model_data
        out["model_data"] = [dataclasses.asdict(d) for d in model_data(self.s, x)]
        out["residual"] = RESIDUAL
        return out

    def _dataset_arrays(self, x, ev):
        """Every data set's arrays at x, in the engine's order, from the evaluation ``ev`` at x.

        The model at each point (the value its residual compares with the datum), the point's c.m. energy (MeV, after
        the data set's energy shift), excitation energy (MeV), c.m. angle (deg) and S-factor conversion
        (``Evaluation.data[i]``).  The datum and its error as the engine read them (c.m., engine unit; pyazr's
        ``cross`` and ``cross_err``), the normalization in use (``data[i].norm``: a THM norm at its profiled value),
        the standardized residuals (the data set's rows of ``Objective.residuals``: see :data:`RESIDUAL`) and the
        sum of their squares.  Two products are formed here, as the program forms them when it writes ``AZUREOut_*``:
        the datum times the norm (``value_norm``, ``error_norm``: what the model is compared with), and, where the S
        factor is offered (:meth:`_has_s_factor`), the model and the normalized datum times the conversion."""
        s = self.s
        residuals = s.objective().residuals(x)
        out = []
        for i, info in enumerate(self.info):
            d_ev = ev.data[i]
            model = d_ev.model
            value, error = np.asarray(s.cross[i], float), np.asarray(s.cross_err[i], float)
            norm = d_ev.norm
            rows = residuals[d_ev.rows]
            d = dict(info, e_cm=_plain(d_ev.energies), e_x=_plain(d_ev.excitation_energies),
                     angle_cm=_plain(d_ev.angles), value=_plain(value), error=_plain(error), norm_used=norm,
                     value_norm=_plain(value * norm), error_norm=_plain(error * norm), model=_plain(model),
                     residual=_plain(rows), residual_chi2=float(np.sum(rows ** 2)))
            if self._has_s_factor(info):
                conv = d_ev.extra["sfactor_conversion"]
                d["s_factor"] = {"model": _plain(model * conv), "value": _plain(value * norm * conv),
                                 "error": _plain(error * norm * conv), "unit": S_FACTOR_UNIT}
            out.append(d)
        return out

    def endf_curves(self, files, x):
        """ENDF-6 files' values at every data set's points at x (``Material.curves(data, session=s)`` of
        :mod:`pyazr.endf`), beside the data set's datum times its normalization and its model, and on the grid of the
        data set's model curve (:meth:`curve`) where it has one; ``files``: [(name, absolute path)], each checked
        inside the root by the server.  See :func:`pyazr.ui.evaluations.curves_answer`."""
        from .evaluations import curves_answer
        x = np.asarray(x, float)
        ev = self.s.evaluate(x)
        curves = self._model_curves(x)
        return curves_answer(files, self.s, ev.data, self._dataset_arrays(x, ev),
                             {n: (c if isinstance(c, Exception) else (c, self._curve_arrays(c)))
                              for n, c in curves.items()})

    def curve(self, name, x):
        """The model of the data set ``name`` on a fine grid at x (``Session.model_curve``; the arrays of
        :class:`pyazr.model_curve.ModelCurve`: c.m. and lab energies, excitation energies, c.m. angles, the model in
        the engine's unit, which grid values are the data set's own points, the branches; the S factor where offered)
        or its refusal (``pyazr.ParameterError``, three parts).  The curves of every active data set are made in one
        session of the extrapolations (``pyazr.model_curve.model_curves``: the folds of the curves are then planned
        on one mesh, as the data session's), kept while its records are the same (they change with an energy shift),
        and evaluated once per x."""
        self._index(name)
        c = self._model_curves(np.asarray(x, float))[name]
        if isinstance(c, Exception):
            raise c
        return self._curve_arrays(c)

    def endf_export(self, x):
        """What an ENDF-6 file of the session at x would lose, before it is written
        (:func:`pyazr.ui.evaluations.export_preview`: ``pyazr.endf.from_session``)."""
        from .evaluations import export_preview
        return export_preview(self.s, np.asarray(x, float))

    def _model_curves(self, x):
        """``{name: ModelCurve or refusal}`` of every active data set at x (see :meth:`curve`)."""
        from pyazr.model_curve import CurveSession, plan_curves
        key = x.tobytes()
        if self._curves_at[0] != key:
            plan = plan_curves(self.s, None, x)
            if self._curves is None or self._curves.plan.text != plan.text:
                self._close_curves()
                self._curves = CurveSession(self.s, plan)
            self._curves_at = (key, self._curves.evaluate(x))
        return self._curves_at[1]

    def _curve_arrays(self, c):
        """A :class:`pyazr.model_curve.ModelCurve` as plain data, with the S factor where it is offered."""
        info = self.info[self._index(c.name)]
        out = {"name": c.name, "quantity": c.quantity, "engine_unit": c.unit, "variable": c.variable,
               "e_cm": _plain(c.energies), "e_lab": _plain(c.energy_lab), "e_x": _plain(c.excitation_energies),
               "angle_cm": _plain(c.angles), "model": _plain(c.model), "at_data": [bool(v) for v in c.at_data],
               "branches": [dataclasses.asdict(b) for b in c.branches], "points": c.n}
        if self._has_s_factor(info):
            out["s_factor"] = {"model": _plain(c.model * c.sfactor_conversion), "unit": S_FACTOR_UNIT}
        return out

    def _close_curves(self):
        if self._curves is not None:
            self._curves.close()
        self._curves, self._curves_at = None, (None, None)

    def _extrapolations(self):
        """The session of the same project and options made for its extrapolations (``extrapolations=True``), made
        on first use and kept with the data session (its external-capture integrals are computed once)."""
        import pyazr
        if self._extrapolation_session is None:
            s = self.s
            self._extrapolation_session = pyazr.Session(s.file, cwd=s.cwd, **dict(s.options, extrapolations=True))
        return self._extrapolation_session

    def _extrapolation_vector(self, x):
        """The data session's free vector ``x`` as the extrapolation session's: its level parameters and abundances
        (the block both lead with), the extrapolations' nuisances at their start."""
        return pyazr_parameters.carry_model_parameters(x, self.s.parameters, self._extrapolations().parameters)

    def _extrapolation_arrays(self, x):
        """Every active extrapolation's arrays at x (the model on its grid: c.m. energy, excitation energy, c.m. angle,
        the value in the engine's unit, and the S factor where offered), from the extrapolation session."""
        infos = [self._extrapolation_info(e) for e in self.s.extrapolations if e.active]
        if not infos:
            return []
        grids = self._extrapolations().evaluate(self._extrapolation_vector(x)).data
        out = []
        for info, d_ev in zip(infos, grids):
            model = d_ev.model
            d = dict(info, e_cm=_plain(d_ev.energies), e_x=_plain(d_ev.excitation_energies),
                     angle_cm=_plain(d_ev.angles), model=_plain(model))
            if self._has_s_factor(info):
                d["s_factor"] = {"model": _plain(model * d_ev.extra["sfactor_conversion"]), "unit": S_FACTOR_UNIT}
            out.append(d)
        return out

    def _extrapolation_info(self, e):
        """What the engine's description says of an extrapolation; its values are in its quantity's engine unit."""
        engine_unit = self.data_kinds.get(e.quantity, {}).get("engine_unit", "")
        return {"name": e.name, "quantity": e.quantity, "unit": engine_unit, "engine_unit": engine_unit,
                "observable": e.observable, "entrance": e.entrance_key, "exit": e.exit_key, "adaptive": e.adaptive}

    def outputs(self, x):
        """Write the program's output files at x into the session's output directory (pyazr's
        ``write_outputs``: the chi-squared evaluation first, then the files, as the program writes them at the
        end of a run); the directory and the names of the files this call wrote."""
        start = time.time_ns()
        directory = self.s.write_outputs(np.asarray(x, float))
        names = sorted(n for n in os.listdir(directory) if os.path.isfile(os.path.join(directory, n))
                       and os.stat(os.path.join(directory, n)).st_mtime_ns >= start)
        return {"directory": directory, "files": names}

    def _has_s_factor(self, info):
        """Whether the S factor is offered for a data set or extrapolation: its engine unit is a cross section and
        its entrance pair is charged."""
        unit = self.units.get(info["engine_unit"]) or {}
        return unit.get("dimension") == S_FACTOR_DIMENSION and self.charged.get(info["entrance"], False)

    def _index(self, name):
        if name not in self.names:
            raise errors.ParameterKeyError(f"no active data set {name!r}{errors.close_names(name, self.names)} "
                                           f"(the data sets: {errors.listed(self.names)}); GET /api/session lists them")
        return self.names.index(name)

    def transform(self, x):
        """Physical values (MeV, eV, fm^-1/2; norms and shifts as they are) of every parameter at x, by name
        (``Session.transform``)."""
        s = self.s
        return {p.name: float(v) for p, v in zip(s.parameters, s.transform(np.asarray(x, float)))}

    def channel_functions(self, pair, e_cm, l=0, radius=0.0, J=None):
        """The session's channel functions of ``pair`` at c.m. energies ``e_cm`` (MeV; a number or a list): F, Fp, G,
        Gp, P, S, delta_hs and open as lists, at most :data:`pyazr.ui.physics.MAX_ELEMENTS` energies."""
        from .physics import check_size
        check_size(np.shape(e_cm), "channel_functions", "split the energies over several requests")
        out = self.s.channel_functions(pair, np.asarray(e_cm, float), l, radius, J)
        return {k: np.asarray(v).tolist() for k, v in out.items()}

    def session_function(self, fn, arguments, x=None):
        """A function of the session's table (``Session.describe()["functions"]``) by name, on ``arguments`` (checked
        by :func:`pyazr.ui.physics.session_request`), at the fit vector x where its row depends on it (``at_x``):
        the session's answer as plain data."""
        from pyazr.session_functions import SESSION_FUNCTIONS
        row = SESSION_FUNCTIONS[fn]
        kwargs = dict(arguments)
        if row["axis"] in kwargs:
            kwargs[row["axis"]] = np.asarray(kwargs[row["axis"]], float)
        if row["at_x"]:
            kwargs["params"] = np.asarray(x, float)
        return _tree(getattr(self.s, fn)(**kwargs))

    # -- the THM experiments (pyazr/session_thm.py) as raw routes ----------------------------------------------------

    def thm_experiments(self, x):
        """Every THM experiment as a chi-squared evaluation at x profiles it (``Session.thm_experiments``): by name,
        its data sets, background, points, chi2, status, norm, background coefficients b, their covariance and sigma,
        and the coherent background's parameters (cbkg) by name."""
        return _tree(self.s.thm_experiments(np.asarray(x, float)))

    def thm_background(self, experiment, x):
        """The shared norm and background of one THM experiment at x (``Session.thm_background``)."""
        return _tree(self.s.thm_background(experiment, np.asarray(x, float)))

    def thm_lineshape(self, experiment, energies, x):
        """The Coulomb line shape of a THM experiment at c.m. energies (MeV) of x + A, the level poles of x
        (``Session.thm_lineshape``); at most :data:`pyazr.ui.physics.MAX_ELEMENTS` energies."""
        return _tree(self.s.thm_lineshape(experiment, _energies(energies, "thm_lineshape"), np.asarray(x, float)))

    def thm_vertex(self, experiment, energies, x):
        """The THM entrance vertex of an experiment at c.m. energies (MeV) of x + A, the vertex boundaries of x
        (``Session.thm_vertex``; energies beyond the window's reach are marked in ``reached`` and warned of)."""
        return _tree(self.s.thm_vertex(experiment, _energies(energies, "thm_vertex"), np.asarray(x, float)))

    def thm_distortion(self, experiment, energies):
        """The distortion factor R(E) of a THM experiment at c.m. energies (MeV) of x + A
        (``Session.thm_distortion``; it does not depend on the R-matrix parameters)."""
        return _tree(self.s.thm_distortion(experiment, _energies(energies, "thm_distortion")))

    def snapshot(self, path, x):
        """The project's text with the values of x written into it (``save_fit``: level values, the free nuisances
        into their records; the result read back and checked), written to ``path``; returns the text."""
        written, _ = self.s.save_fit(path, np.asarray(x, float), param_sav=False, relative_to=self.s.cwd)
        with open(written, "rb") as f:
            return f.read().decode("utf-8", "surrogateescape")


def required_arguments(fn):
    """The arguments of the session's method ``fn`` that have no default (a row's ``default`` None says none)."""
    from pyazr import Session
    signature = inspect.signature(getattr(Session, fn))
    return {n for n, p in signature.parameters.items() if n != "self" and p.default is inspect.Parameter.empty}


def _energies(energies, what):
    from .physics import check_size
    check_size(np.shape(energies), what, "split the energies over several requests")
    return np.asarray(energies, float)


def _tree(value):
    """Plain data of a nested answer: dicts, lists, numbers and texts (numpy arrays and scalars as lists and
    numbers)."""
    if isinstance(value, dict):
        return {str(k): _tree(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_tree(v) for v in value]
    if isinstance(value, np.ndarray):
        return _tree(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    return value


def _dataset_info(d):
    """What the engine's description says of a data set (pyazr.described.DescribedDataset)."""
    return {"name": d.name, "quantity": d.quantity, "unit": d.unit, "engine_unit": d.engine_unit,
            "observable": d.observable, "entrance": d.entrance_key, "exit": d.exit_key}


def session_text(path):
    """The Python line that opens a session on the project file ``path`` (the Physics tab's Copy as Python)."""
    return f"session = pyazr.Session({path!r})  # the saved file; the page evaluates the working text"


# -- for a job's worker (pyazr/ui/jobs.py): its own session, the objective, the engine's fitters -----------------

def open_for_job(file, cwd):
    """For a job: ``(adapter, x0, names of x)`` of a new session of ``file`` (the adapter's ``s`` is the session)."""
    e = Engine()
    e.open(file, cwd)
    return e, e.s.parameters.x0, e._free_names()


@contextlib.contextmanager
def session_named(s, file_name):
    """A session of the text and options of ``s`` read from a copy named ``file_name`` in a directory of its own
    (relative paths still resolved in ``s.cwd``), closed and removed at the end: for a call that names the session's
    file in what it writes (``pyazr.endf.write``'s MF=1 text names the project, not the UI's working file)."""
    import shutil
    import tempfile
    import pyazr
    work = tempfile.mkdtemp(prefix="pyazr_ui_named_")
    path = os.path.join(work, file_name)
    shutil.copyfile(s.file, path)
    named = pyazr.Session(path, cwd=s.cwd, **s.options)
    try:
        yield named
    finally:
        named.close()
        shutil.rmtree(work, ignore_errors=True)


def objective_of(e):
    """The session's :class:`pyazr.Objective` (``value_and_grad``, ``log_prob``, ``bounds``, ``names``)."""
    return e.s.objective()


def fit_with_engine(e, method, minos, x0, max_calls, progress, stop):
    """The engine's fitter ``method`` (a row of ``pyazr.kinds()["fitters"]``) on the session from ``x0``
    (``Session.fit_with_engine``): ``progress`` gets each progress record (a dict), ``stop`` (a threading.Event) ends
    the fit at its best point.  Plain data: the end point by name, the objective there, the engine's status and
    words, the outcome in the words every fitter shares, its notices, EDM (None where the fitter has none), calls,
    errors, covariance and MINOS intervals by name where the fitter gives them (None otherwise) with each side's
    validity, the fitter's log."""
    r = e.s.fit_with_engine(method, minos=bool(minos), x0=np.asarray(x0, float), max_calls=max_calls,
                            progress=progress, stop=stop)
    names = list(r.names)
    covariance = None if r.covariance is None else np.asarray(r.covariance, float)
    return {"method": r.method, "x": dict(zip(names, np.asarray(r.x, float).tolist())),
            "objective": float(r.objective), "status": r.status, "converged": bool(r.converged),
            "edm": float(r.edm) if np.isfinite(r.edm) else None, "calls": int(r.calls), "names": names,
            "covariance": None if covariance is None else covariance.tolist(),
            "errors": None if r.errors is None else dict(zip(names, np.asarray(r.errors, float).tolist())),
            "minos": {n: [float(lo), float(hi)] for n, (lo, hi) in r.minos.items()},
            "covariance_from_hesse": bool(r.covariance_from_hesse),
            "fell_back_to_migrad": bool(r.fell_back_to_migrad), "outcome": r.outcome, "notices": list(r.notices),
            "minos_valid": {n: {side: dict(flags) for side, flags in sides.items()}
                            for n, sides in r.minos_valid.items()}, "log": r.log}


def band_points(e, x):
    """Where a band is drawn: every active data set, then every active extrapolation, each with its name, c.m.
    energy (MeV), excitation energy (MeV), c.m. angle (deg), the model at x (the engine's unit), and the S-factor
    conversion of each point where the S factor is offered (None otherwise), in the order of the Results payload."""
    x = np.asarray(x, float)
    ev = e.s.evaluate(x)
    out = {"datasets": [_band_entry(e, info, d) for info, d in zip(e.info, ev.data)], "extrapolations": []}
    infos = [e._extrapolation_info(x_) for x_ in e.s.extrapolations if x_.active]
    if infos:
        grids = e._extrapolations().evaluate(e._extrapolation_vector(x)).data
        out["extrapolations"] = [_band_entry(e, info, d) for info, d in zip(infos, grids)]
    return out


def _band_entry(e, info, d):
    conv = d.extra["sfactor_conversion"] if e._has_s_factor(info) else None
    return {"name": info["name"], "e_cm": _plain(d.energies), "e_x": _plain(d.excitation_energies),
            "angle_cm": _plain(d.angles), "model": _plain(d.model), "conv": None if conv is None else _plain(conv)}


def models_at(e, xs, progress):
    """The model of every active data set and extrapolation at each free vector of ``xs`` (rows): ``(data sets,
    extrapolations)``, each a list of arrays of shape (len(xs), points); ``progress(k, objective or None, x)`` after
    each evaluation (the objective in data mode; extrapolations have none)."""
    xs = np.asarray(xs, float)
    data = []
    for k, x in enumerate(xs):
        ev = e.s.evaluate(x)
        if not data:
            data = [np.empty((len(xs), len(d.model))) for d in ev.data]
        for i, d in enumerate(ev.data):
            data[i][k] = d.model
        progress(k + 1, ev.objective, x)
    grids = []
    if any(x_.active for x_ in e.s.extrapolations):
        ext = e._extrapolations()
        for k, x in enumerate(xs):
            ev = ext.evaluate(e._extrapolation_vector(x))
            if not grids:
                grids = [np.empty((len(xs), len(d.model))) for d in ev.data]
            for i, d in enumerate(ev.data):
                grids[i][k] = d.model
            progress(len(xs) + k + 1, None, x)
    return data, grids


def linear_sigma(e, x, covariance, names):
    """The 1-sigma of the model of every active data set and extrapolation at x from a fit's covariance over the free
    vector ``names`` (``pyazr.bands.uncertainty_bands``: sigma^2 = g^T C g with the engine's sensitivities g, over
    the columns a band spans, ``pyazr.bands.rmatrix_columns``: the nuisances do not enter).  ``(data sets,
    extrapolations, the names of the parameters the band varied)``: lists of arrays in the order of
    :func:`band_points`."""
    from pyazr import bands
    free = e._free_names()
    if list(names) != free:
        raise errors.ParameterError(f"the covariance is over {len(names)} parameters ({errors.listed(names)}), the "
                                    f"session's free vector over {len(free)} ({errors.listed(free)}): a band needs "
                                    "the fit of the same project text")
    x = np.asarray(x, float)
    cols = bands.rmatrix_columns(e.s)
    block = np.asarray(covariance, float)[np.ix_(cols, cols)]
    n_data = sum(1 for d in e.s.datasets if d.active)
    data = bands.uncertainty_bands(e.s, segments=list(range(n_data)), params=x, covariance=block)
    out_data = [data[key].sigma for key in sorted(data, key=lambda k: data[k].segment)]
    out_grids = []
    if any(x_.active for x_ in e.s.extrapolations):
        n_grids = sum(1 for x_ in e.s.extrapolations if x_.active)
        grids = bands.uncertainty_bands(e._extrapolations(), segments=list(range(n_grids)), params=x[cols],
                                        covariance=block)
        out_grids = [grids[key].sigma for key in sorted(grids, key=lambda k: grids[k].segment)]
    return out_data, out_grids, [free[c] for c in cols]


def _read_chi_squared(path):
    """The rows and the Total- line of the program's ``chiSquared.out``: ([{chi2, n, norm, norm_chi2}], {name: v})."""
    rows, totals = [], {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Segment#"):
                continue
            if line.startswith("Total-"):
                tokens = line.replace(":", " ").split()
                totals.update({tokens[k]: float(tokens[k + 1]) for k in range(0, len(tokens) - 1, 2)})
                continue
            cells = line.split(",")
            rows.append({"chi2": float(cells[1]), "n": int(cells[2]), "norm": float(cells[3]),
                         "norm_chi2": float(cells[4])})
    return rows, totals
