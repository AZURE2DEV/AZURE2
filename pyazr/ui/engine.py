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
import os

import numpy as np

from pyazr import errors, session_checks

# Today's methods of pyazr.azure2 the UI does not call, by reason (the route-table consistency test,
# tests/pyazr/ui/routes_consistency_test.py, requires every public method to be called here or listed here).
AZURE2_NOT_SERVED = {
    "replaced by evaluate / evaluation (one call per evaluation, pyazr3.md 2.3)":
        "calculate calculate_all_rwa calculate_sfactor calculate_sfactor_rwa calculate_analyzing_power "
        "calculate_analyzing_power_rwa calculate_angular_dists calculate_angular_dists_rwa calculate_energies "
        "calculate_angles calculate_excitation_energy calculate_chi2 segment_chi2 dataset_chi2 penalties "
        "scattering_length_chi2 chi2_and_grad residual_jacobian model_gradients segment_norms",
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
    "the mode switch: extrapolations are a later step": "data_mode extrap_mode",
    "the server owns sessions: lifecycle": "close is_alive",
    "later steps of web_ui.md section 11 (Nucleus, Physics, THM tabs)":
        "level_scheme pairs nuclei mass_number basis ec_integrals thm_experiments thm_background "
        "thm_lineshape thm_vertex thm_distortion prior_covariance_residuals prior_covariance_names extrapolations "
        "datasets n_rmatrix cache_stats output_dir",
    "the program's output files: a later route (/api/outputs)": "write_output_files",
    "the older name of channel_functions (served: POST /api/session/channel_functions)": "coulomb_functions",
}


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

    def data(self, name, x, models=None):
        """The arrays of the data set ``name`` at x: c.m. energies (MeV) and angles (deg), the data values and errors
        in the engine's units, the model, and the engine's standardized residuals of its rows (whitened for a data set
        with a data covariance)."""
        s = self.s
        i = self._index(name)
        x = np.asarray(x, float)
        model = (models if models is not None else s.calculate_rwa(x))[i]
        edges = np.cumsum([0] + [len(e) for e in s.energies])
        residuals = np.asarray(s.residuals(x), float)[edges[i]:edges[i + 1]]
        return {**self.info[i],
                "e_cm": _plain(s.energies[i]), "e_x": _plain(s.excitation_energies[i]),
                "angle_cm": _plain(s.angles[i]), "value": _plain(s.cross[i]), "error": _plain(s.cross_err[i]),
                "model": _plain(model), "residual": residuals.tolist()}

    def evaluation(self, x):
        """Every data set's arrays at x (one model pass for all)."""
        models = self.s.calculate_rwa(np.asarray(x, float))
        return {"datasets": [self.data(name, x, models) for name in self.names]}

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
