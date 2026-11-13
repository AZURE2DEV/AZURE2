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
import dataclasses
import os
import time

import numpy as np

from pyazr import errors, session_checks
from pyazr.locking import session_lock

# The methods of pyazr.Session the UI does not call, by reason (the route-table consistency test,
# tests/pyazr/ui/routes_consistency_test.py, requires every public method to be called here or listed here).
SESSION_NOT_SERVED = {
    "the sensitivities of the model values (uncertainty bands, not a route of this step)": "model_gradients",
    "dimensionless widths and level quantities, not a route of this step":
        "wigner_widths dimensionless_widths park_norms scattering_lengths",
    "session edits: edits of the Project and a new session (pyazr3.md decision 11)":
        "set_channel_radius configure recalculate_external_capture "
        "refreeze_covariance set_nuclear_potential clear_nuclear_potential nuclear_potential nuclear_potentials "
        "potential_records",
    "the server owns sessions: lifecycle": "close is_alive",
    "later steps of web_ui.md section 11 (Nucleus, Physics, THM tabs)":
        "level_scheme nuclei mass_number basis ec_integrals thm_experiments thm_background "
        "thm_lineshape thm_vertex thm_distortion prior_covariance_residuals prior_covariance_names "
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

    # -- the session ----------------------------------------------------------------------------------------------

    def open(self, file, cwd):
        """Build the session of ``file`` with relative paths resolved in ``cwd`` (the project's directory)."""
        import pyazr
        if self.s is not None:
            self.s.close()
            self.s = None
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
                "datasets": self._datasets(),
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

    def _extrapolation_arrays(self, x):
        """Every active extrapolation's arrays at x (the model on its grid: c.m. energy, excitation energy, c.m. angle,
        the value in the engine's unit, and the S factor where offered).  The session switches to its extrapolation
        grids and back (``extrap_mode``, ``data_mode``; the first switch computes the grids' external-capture
        integrals, later ones reuse them)."""
        s = self.s
        infos = [self._extrapolation_info(e) for e in s.extrapolations if e.active]
        if not infos:
            return []
        s.extrap_mode()
        try:
            grids = s.evaluate(x).data
            out = []
            for info, d_ev in zip(infos, grids):
                model = d_ev.model
                d = dict(info, e_cm=_plain(d_ev.energies), e_x=_plain(d_ev.excitation_energies),
                         angle_cm=_plain(d_ev.angles), model=_plain(model))
                if self._has_s_factor(info):
                    d["s_factor"] = {"model": _plain(model * d_ev.extra["sfactor_conversion"]), "unit": S_FACTOR_UNIT}
                out.append(d)
            return out
        finally:
            s.data_mode()

    def _extrapolation_info(self, e):
        """What the engine's description says of an extrapolation; its values are in its quantity's engine unit."""
        engine_unit = self.data_kinds.get(e.quantity, {}).get("engine_unit", "")
        return {"name": e.name, "quantity": e.quantity, "unit": engine_unit, "engine_unit": engine_unit,
                "observable": e.observable, "entrance": e.entrance_key, "exit": e.exit_key, "adaptive": e.adaptive}

    def outputs(self, x):
        """Write the program's output files at x into the session's output directory (pyazr's
        ``write_output_files``: the chi-squared evaluation first, then the files, as the program writes them at the
        end of a run); the directory and the names of the files this call wrote."""
        start = time.time_ns()
        directory = self.s.write_output_files(np.asarray(x, float))
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

    def snapshot(self, path, x):
        """The project's text with the values of x written into it (``save_fit``: level values, the free nuisances
        into their records; the result read back and checked), written to ``path``; returns the text."""
        written, _ = self.s.save_fit(path, np.asarray(x, float), param_sav=False)
        with open(written, "rb") as f:
            return f.read().decode("utf-8", "surrogateescape")


def _dataset_info(d):
    """What the engine's description says of a data set (pyazr.described.DescribedDataset)."""
    return {"name": d.name, "quantity": d.quantity, "unit": d.unit, "engine_unit": d.engine_unit,
            "observable": d.observable, "entrance": d.entrance_key, "exit": d.exit_key}


def session_text(path):
    """The Python line that opens a session on the project file ``path`` (the Physics tab's Copy as Python)."""
    return f"session = pyazr.Session({path!r})  # the saved file; the page evaluates the working text"


def open_objective(file, cwd):
    """For a job: ``(value_and_grad(x) -> (f, g), x0, names of x)`` of a new session of ``file`` (its
    :class:`pyazr.Objective`)."""
    e = Engine()
    e.open(file, cwd)
    return e.s.objective().value_and_grad, e.s.parameters.x0, e._free_names()


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
