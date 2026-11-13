"""The adapter of the web UI to today's pyazr session: the only file of ``pyazr/ui`` that names the class
:class:`pyazr.azure2` and its methods (docs/upgrade/web_ui.md section 11).  When pyazr 3.0's ``Session``,
``Parameters`` and ``Evaluation`` replace ``azure2`` (pyazr3.md step 5), this file is the one that changes.

It runs inside a worker process (:mod:`pyazr.ui.workers`), never in the server: a session build holds the GIL for its
whole length (12.8 s on 13N, web_ui.md section 2.1).  Every answer is plain JSON data (lists, dicts, floats); a
parameter vector ``x`` is the engine's free vector in its order and basis (the fit vector, pyazr3.md decision 10).

The chi-squared per data set is the program's: the rows of the ``chiSquared.out`` the session writes into the UI's
own output directory (``<config>``'s output directory of the session's file points there, see
:meth:`pyazr.ui.state.State.session_text`), not ``segment_chi2``'s split, which differs for data sets joined by a data
covariance (web_ui.md section 2.2, finding 5).  The file has the program's six significant digits.
"""
import dataclasses
import os
import time

import numpy as np

from pyazr import errors, session_checks

# Today's methods of pyazr.azure2 the UI does not call, by reason (the route-table consistency test,
# tests/pyazr/ui/routes_consistency_test.py, requires every public method to be called here or listed here).
AZURE2_NOT_SERVED = {
    "replaced by evaluate / evaluation (one call per evaluation, pyazr3.md 2.3)":
        "calculate calculate_all_rwa calculate_sfactor calculate_sfactor_rwa calculate_analyzing_power "
        "calculate_analyzing_power_rwa calculate_angular_dists calculate_angular_dists_rwa calculate_energies "
        "calculate_angles calculate_excitation_energy calculate_chi2 dataset_chi2 penalties "
        "scattering_length_chi2 chi2_and_grad residual_jacobian model_gradients",
    "physical values only through transform (pyazr3.md decision 10)":
        "transform_physical transform_rwa physical_levels wigner_widths dimensionless_widths park_norms "
        "scattering_lengths",
    "parameter-vector files: Parameters.read/write of pyazr 3.0, not a route of this step":
        "read_sav load_sav full_rwa_from_sav read_sav_basis update_rwa_params_from_sav update_sav_from_rwa_params",
    "per-kind index views: the parameter table carries kind and free index":
        "level_free_width_indices without_level only_level norm_indices shift_indices energy_indices",
    "session edits: edits of the Project and a new session (pyazr3.md decision 11)":
        "set_channel_radius refresh_parameters refresh_pairs configure recalculate_external_capture "
        "refreeze_covariance set_nuclear_potential clear_nuclear_potential nuclear_potential nuclear_potentials "
        "potential_records",
    "the objective with its terms in one call (objective_terms_rwa; the same number)": "objective",
    "the server owns sessions: lifecycle": "close is_alive",
    "later steps of web_ui.md section 11 (Nucleus, Physics, THM tabs)":
        "level_scheme nuclei mass_number basis ec_integrals thm_experiments thm_background "
        "thm_lineshape thm_vertex thm_distortion prior_covariance_residuals prior_covariance_names "
        "datasets n_rmatrix cache_stats output_dir",
    "the older name of channel_functions (served: POST /api/session/channel_functions)": "coulomb_functions",
}


# The residual the payload carries for every point of a data set: the engine's standardized residual, the one whose
# square is its chi-squared term (pyazr's ``residuals``, the program's chi-squared in Python).
RESIDUAL = {
    "definition": "(model - n*datum) / (n*error)",
    "doc": "the engine's standardized residual of each point: model minus the datum times the data set's "
           "normalization n, over the error times n, in the data set's own frame and unit (energy shifts applied). "
           "Its squares add up to the data set's chi-squared of chiSquared.out (to the file's six digits) unless a "
           "data covariance joins the data set with others: then the residuals are whitened rows w = W r "
           "(W^T W = M^-1), whose squares add up to the group's chi-squared but split it differently from the "
           "program's sum r_i (M^-1 r)_i per data set (web_ui.md finding 5); residual_chi2 is the sum of their "
           "squares (pyazr's segment_chi2).",
}
# The dimension (pyazr.kinds()["units"]) of a data set's engine unit for which the S factor is offered, with a
# charged entrance pair (the S factor's own rule in the engine's quantity table); the unit of the converted values.
S_FACTOR_DIMENSION = "cross_section"
S_FACTOR_UNIT = "MeV*b"

def _plain(a):
    return np.asarray(a, float).tolist()


class Engine:
    """One session of today's pyazr (``pyazr.azure2``) and what the routes ask of it."""

    def __init__(self):
        self.s = None

    # -- the session ----------------------------------------------------------------------------------------------

    def open(self, file, cwd):
        """Build the session of ``file`` with relative paths resolved in ``cwd`` (the project's directory)."""
        import pyazr
        if self.s is not None:
            self.s.close()
            self.s = None
        s = pyazr.azure2(file, cwd=cwd)
        os.makedirs(s.output_dir, exist_ok=True)  # the external-capture cache and chiSquared.out are written there
        self.s = s
        self.info = [_dataset_info(d) for d in s.active_datasets]  # in the engine's order (a copy per access)
        self.names = [d["name"] for d in self.info]
        kinds = pyazr.kinds()
        self.units, self.data_kinds = kinds["units"], kinds["data"]
        self.charged = {p.key: p.Z1 * p.Z2 > 0 for p in s.pairs}  # pair key -> both constituents charged
        return self._describe()

    def _describe(self):
        """The session: its options, free-vector length, data sets (in the engine's order) and parameters."""
        s = self.s
        return {"options": dict(s.options), "n_free": len(s.params_rwa), "datasets": self._datasets(),
                "parameters": self._parameters()}

    def _datasets(self):
        return [dict(d, points=len(self.s.energies[i])) for i, d in enumerate(self.info)]

    def _parameters(self):
        """One row per parameter of the engine (fixed and free), the session's own description (``describe()``): name,
        kind, unit, whether it is free and its index in x, ``x`` (fit basis), ``physical_value``, bounds, owner
        fields; ``x0`` is ``x`` of a free parameter."""
        rows = self.s.describe()["parameters"]
        for r in rows:
            r["x0"] = r["x"] if r["free"] else None
        return rows

    def _free_names(self):
        """Names of the free parameters, in the order of x."""
        return [p.name for p in self.s.parameters if p.free_index is not None and not p.fixed]

    # -- evaluations ----------------------------------------------------------------------------------------------

    def evaluate(self, x):
        """The objective with its named terms, the data chi-squared, and the program's chi-squared per data set
        (``chiSquared.out``: chi2, N, norm, norm chi2 per row; its Total- line as ``totals``)."""
        s = self.s
        x = np.asarray(x, float)
        terms = s.sess.objective_terms_rwa(x)
        # azure2.write_output_files(x) is these two calls; made here so that the data chi-squared it computes is kept
        chi2 = float(s.sess.calculate_chi2_rwa(x))
        os.makedirs(s.output_dir, exist_ok=True)
        if not s.sess.write_output_files():
            raise RuntimeError("the engine wrote no output files (no data object)")
        session_checks.warn_engine_notices(s.sess, stacklevel=2)  # the engine's notices of these calls as warnings
        rows, totals = _read_chi_squared(os.path.join(s.output_dir, "chiSquared.out"))
        if len(rows) != len(self.names):
            raise RuntimeError(f"chiSquared.out has {len(rows)} rows, the session {len(self.names)} data sets")
        per_set = [dict(name=n, **r) for n, r in zip(self.names, rows)]
        return {"objective": float(terms["total"]), "chi2": chi2,
                "terms": {str(k): float(v) for k, v in zip(terms["names"], terms["values"])},
                "datasets": per_set, "totals": totals}

    def data(self, name, x):
        """The arrays of the data set ``name`` at x (:meth:`evaluation`'s entry for it, without the chi-squared
        table)."""
        i = self._index(name)
        return self._dataset_arrays(np.asarray(x, float))[i]

    def evaluation(self, x):
        """The Results payload at x: the program's chi-squared split (:meth:`evaluate`, which also writes the program's
        output files into the session's output directory), every data set's arrays (its chi-squared row merged in),
        every extrapolation's arrays and the model data.  Every number is the engine's or pyazr's (see
        :meth:`_dataset_arrays` for the two products the adapter forms)."""
        x = np.asarray(x, float)
        out = self.evaluate(x)
        arrays = self._dataset_arrays(x)
        out["datasets"] = [dict(a, **row) for a, row in zip(arrays, out["datasets"])]
        out["extrapolations"] = self._extrapolation_arrays(x)
        from pyazr import model_data
        out["model_data"] = [dataclasses.asdict(d) for d in model_data(self.s, x)]
        out["residual"] = RESIDUAL
        return out

    def _dataset_arrays(self, x):
        """Every data set's arrays at x, in the engine's order.

        From one forward pass (``calculate_rwa``): the model at each point (the value its residual compares with the
        datum), the point's c.m. energy (MeV, after the data set's energy shift), excitation energy (MeV), c.m. angle
        (deg) and S-factor conversion.  The datum and its error as the engine read them (c.m., engine unit; pyazr's
        ``cross`` and ``cross_err``), the normalization in use (``segment_norms``: a THM norm at its profiled value),
        the standardized residuals (``residuals``: see :data:`RESIDUAL`) and the sum of their squares per data set
        (``segment_chi2``).  Two products are formed here, as the program forms them when it writes ``AZUREOut_*``:
        the datum times the norm (``value_norm``, ``error_norm``: what the model is compared with), and, where the S
        factor is offered (:meth:`_has_s_factor`), the model and the normalized datum times the conversion."""
        s, sess = self.s, self.s.sess
        norms = np.asarray(s.segment_norms(x), float)
        sums = np.asarray(s.segment_chi2(x), float)
        residuals = np.asarray(s.residuals(x), float)
        models = s.calculate_rwa(x)
        edges = np.cumsum([0] + [len(m) for m in models])
        out = []
        for i, info in enumerate(self.info):
            model = np.asarray(models[i], float)
            value, error = np.asarray(s.cross[i], float), np.asarray(s.cross_err[i], float)
            norm = float(norms[i])
            d = dict(info, e_cm=_plain(sess.calculated_energies(i)), e_x=_plain(sess.calculated_excitation_energies(i)),
                     angle_cm=_plain(sess.calculated_angles(i)), value=_plain(value), error=_plain(error), norm_used=norm,
                     value_norm=_plain(value * norm), error_norm=_plain(error * norm), model=_plain(model),
                     residual=_plain(residuals[edges[i]:edges[i + 1]]), residual_chi2=float(sums[i]))
            if self._has_s_factor(info):
                conv = np.asarray(sess.calculated_conv(i), float)
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
            models = s.calculate_rwa(x)
            sess = s.sess
            out = []
            for i, info in enumerate(infos[:len(models)]):
                model = np.asarray(models[i], float)
                d = dict(info, e_cm=_plain(sess.calculated_energies(i)),
                         e_x=_plain(sess.calculated_excitation_energies(i)), angle_cm=_plain(sess.calculated_angles(i)),
                         model=_plain(model))
                if self._has_s_factor(info):
                    d["s_factor"] = {"model": _plain(model * np.asarray(sess.calculated_conv(i), float)),
                                     "unit": S_FACTOR_UNIT}
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
        """Physical values (MeV, eV, fm^-1/2; norms and shifts as they are) of every parameter at x, by name."""
        s = self.s
        full = s._all_rwa(np.asarray(x, float))
        values = np.asarray(s.transform_all_rwa(full, include_fixed=True), float)
        return {p.name: float(values[i]) for i, p in enumerate(s.parameters) if i < values.size}

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
    return f"session = pyazr.azure2({path!r})  # the saved file; the page evaluates the working text"


def open_objective(file, cwd):
    """For a job: ``(objective_and_grad(x) -> (f, g), x0, names of x)`` of a new session of ``file``."""
    e = Engine()
    e.open(file, cwd)
    return e.s.objective_and_grad, np.asarray(e.s.params_rwa, float), e._free_names()


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
