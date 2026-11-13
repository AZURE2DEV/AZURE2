"""The numbers of a session at one parameter point (docs/upgrade/pyazr3.md section 7).

:func:`evaluate` (``Session.evaluate(x)``) runs one evaluation of the engine and returns an :class:`Evaluation`: the
objective, its terms by the engine's names, the data chi-squared, and every data set's arrays, read where the
program's output writer reads them (``chiSquared.out``, ``AZUREOut_*``).  Every number is the engine's; the objects
hold copies the caller owns, and no session.
"""
import numpy as np

from .model_data import ModelDatum
from .locking import session_lock


class EvaluatedData:
    """One data set at the evaluation's x (``Evaluation.data[i]``), in the engine's order.

    - ``name``, ``quantity``, ``unit`` (the engine's unit of ``model``, ``values``, ``errors``), ``frame``
      (``"c.m."``): from the session's description of the project; ``segment`` (the engine's segment number, from 1)
      and ``segment_key``.
    - ``chi2``: the data set's row of ``chiSquared.out`` (a data covariance by its per-set partition
      sum r_i (M^-1 r)_i); ``norm_chi2``, ``shift_chi2``: its normalization prior and its energy-shift priors (the
      three parts summed), as the objective adds them, and the three parts apart: ``shift_constant_chi2`` (the
      constant shift), ``shift_sqrt_chi2`` (the sqrt(E) coefficient), ``shift_function_chi2`` (the parameters of a
      shift function); ``norm``: the normalization in use (a profiled norm at n*).
    - ``rows``: the data set's rows of the objective's least-squares rows, a ``slice`` (``Objective.residuals(x)
      [d.rows]``; one row per point, whitened where a data covariance joins points).  With a covariance joining data
      sets only the sum over the group has a meaning: ``chi2`` splits it as ``chiSquared.out`` does,
      ``sum(r[d.rows]**2)`` by the whitened rows; both add up to the group's chi-squared.
    - In extrapolation mode (the session's ``<segmentsTest>`` grids) there are no data: every chi-squared field,
      ``norm``, ``values``, ``errors``, ``residuals`` and ``rows`` are ``None``; the arrays come from the engine's
      forward pass, and ``extra`` holds ``sfactor_conversion`` and, for an angular distribution, ``legendre``.
    - ``energies`` (MeV, c.m., after the energy shift), ``angles`` (deg, c.m.), ``excitation_energies`` (MeV),
      ``model`` (the model value of each point, c.m.), ``values``, ``errors`` (the datum and its error as the engine
      compares them, c.m.), ``residuals`` ((f/J - y n)/(dy n) in the datum's frame, not whitened; a THM experiment's
      rows of its shared profile).  The set-up arrays (energies, angles, values, errors) are copies here; they become
      read-only views in a later release.
    - ``extra``: ``{name: array}`` of the numbers the points carry beside the model value that apply to this data set
      (``pyazr.kinds()["data_outputs"]``: the lab energy and angle, the S-factor conversion, the E1 and E2 parts of a
      capture, a list of Legendre coefficients or of a sample's per-nucleus terms per point).
    """

    BUILT_BY = "Evaluation.data[i]"
    FIELDS = {
        "name": "the data set's record name", "quantity": "its quantity", "unit": "the engine unit of model, values",
        "frame": "c.m.", "segment": "the engine's segment number", "segment_key": "the segment's key",
        "chi2": "the data set's row of chiSquared.out", "norm_chi2": "its normalization prior",
        "shift_chi2": "its energy-shift priors", "shift_constant_chi2": "the constant shift's prior",
        "shift_sqrt_chi2": "the sqrt(E) coefficient's prior", "shift_function_chi2": "the shift function's priors",
        "rows": "a slice: its rows of Objective.residuals(x)", "norm": "the normalization in use",
        "energies": "MeV, c.m., after the shift", "angles": "deg, c.m.", "excitation_energies": "MeV",
        "model": "the model value of each point", "values": "the datum (c.m., engine unit)",
        "errors": "its error", "residuals": "(f/J - y n)/(dy n) per point", "n": "the number of points",
        "extra": "{name: array} of the point outputs that apply (kinds()['data_outputs'])",
        "table()": "the points as text",
    }

    CHI2_FIELDS = ("chi2", "norm_chi2", "shift_chi2", "shift_constant_chi2", "shift_sqrt_chi2", "shift_function_chi2",
                   "norm")

    def __init__(self, row, description, has_data=True):
        self.name = description.name if description is not None else None
        self.quantity = description.quantity if description is not None else None
        self.unit = getattr(description, "engine_unit", None) if description is not None else None
        self.frame = "c.m."
        self.segment = int(row["segment"])
        self.segment_key = int(row["segment_key"]) if has_data else getattr(description, "key", -1)
        for key in self.CHI2_FIELDS:
            setattr(self, key, float(row[key]) if has_data else None)
        self.rows = slice(int(row["first_row"]), int(row["first_row"]) + int(row["rows"])) if has_data else None
        for key in ("energies", "angles", "excitation_energies", "model", "values", "errors", "residuals"):
            setattr(self, key, np.asarray(row[key], float))
        if not has_data:
            self.values = self.errors = self.residuals = None
        self.extra = {k: (np.asarray(v, float) if not isinstance(v, list) else [np.asarray(a, float) for a in v])
                      for k, v in row["extra"].items()}

    @property
    def n(self):
        """The number of points."""
        return len(self.model)

    def table(self):
        """The points as text: c.m. energy, angle, model, datum, error, residual, one line each (energy, angle and
        model alone for an extrapolation)."""
        if self.chi2 is None:
            lines = [f"{self.name}  ({self.quantity}, {self.unit}, {self.frame}; {self.n} points)",
                     f"{'E_cm/MeV':>14} {'angle/deg':>10} {'model':>14}"]
            lines += [f"{self.energies[k]:14.7g} {self.angles[k]:10.5g} {self.model[k]:14.7g}" for k in range(self.n)]
            return "\n".join(lines)
        lines = [f"{self.name}  ({self.quantity}, {self.unit}, {self.frame}; chi2 {self.chi2:.6g}, "
                 f"norm {self.norm:.6g}, {self.n} points)",
                 f"{'E_cm/MeV':>14} {'angle/deg':>10} {'model':>14} {'datum':>14} {'error':>12} {'residual':>10}"]
        for k in range(self.n):
            lines.append(f"{self.energies[k]:14.7g} {self.angles[k]:10.5g} {self.model[k]:14.7g} "
                         f"{self.values[k]:14.7g} {self.errors[k]:12.5g} {self.residuals[k]:10.4g}")
        return "\n".join(lines)

    def __repr__(self):
        chi2 = "" if self.chi2 is None else f", chi2={self.chi2:.6g}"
        return f"EvaluatedData({self.name!r}, {self.quantity}, {self.n} points{chi2})"


class Evaluation:
    """The numbers of a session at one point x of the free vector, from one engine evaluation.

    - ``x``: a copy of the point.
    - ``objective``: what the fit minimises, ``Session.objective(x)`` and ``Objective(x)`` bit for bit.
    - ``chi2``: the data chi-squared, ``Total-Chi-Squared`` of the program's ``chiSquared.out`` (THM norms profiled,
      a data covariance included); ``objective - chi2`` are the priors and constraints.
    - ``terms``: ``{name: value}``, one entry per term of the engine's objective (its table, ``src/Objective.cpp``),
      in its order; they add up to ``objective`` (to rounding: each is summed on its own).  The first term holds the
      data, the normalization priors and the energy-shift priors; ``data[i].norm_chi2`` and ``.shift_chi2`` split the
      priors per data set, as ``chiSquared.out`` does.
    - ``data``: one :class:`EvaluatedData` per active data set, in the engine's order.
    - ``model_data``: the ``<modelData>`` rows (:class:`pyazr.model_data.ModelDatum`).
    - ``log_prob``: ``Objective.log_prob(x)``; ``excluded``: the point is outside the engine's sampler's support.
    - ``has_data``: False for an evaluation of the extrapolations (a session in extrapolation mode): ``data`` are the
      extrapolations, and ``objective``, ``chi2``, ``log_prob``, ``excluded`` are ``None``, ``terms`` is empty.
    """

    BUILT_BY = "Session.evaluate(x)"
    FIELDS = {
        "x": "a copy of the point", "objective": "what the fit minimises (= Objective(x))",
        "chi2": "the data chi-squared, Total-Chi-Squared of chiSquared.out",
        "terms": "{name: value} per term of the engine's objective", "data": "one EvaluatedData per data set",
        "model_data": "the <modelData> rows", "log_prob": "Objective.log_prob(x)",
        "excluded": "outside the engine's sampler's support", "has_data": "False for the extrapolations",
        "table(i)": "text", "dataset(name)": "by name",
    }

    def __init__(self, x, raw, descriptions, has_data=True):
        self.x = np.array(x, float)
        self.has_data = bool(has_data)
        self.objective = float(raw["objective"]) if has_data else None
        self.chi2 = float(raw["chi2"]) if has_data else None
        self.terms = {str(k): float(v) for k, v in zip(raw["term_names"], raw["term_values"])} if has_data else {}
        self.log_prob = float(raw["log_prob"]) if has_data else None
        self.excluded = bool(raw["excluded"]) if has_data else None
        rows = list(raw["data"])
        if descriptions is not None and len(descriptions) != len(rows):
            descriptions = None  # the description and the engine disagree: the arrays keep their segment numbers
        self.data = [EvaluatedData(r, descriptions[i] if descriptions is not None else None, has_data)
                     for i, r in enumerate(rows)]
        self.model_data = [ModelDatum(**{k: (float(v) if isinstance(v, float) else v) for k, v in r.items()})
                           for r in raw.get("model_data", [])]

    def dataset(self, name):
        """The :class:`EvaluatedData` of the data set ``name``."""
        for d in self.data:
            if d.name == name:
                return d
        raise KeyError(f"no data set {name!r} in this evaluation (the data sets: {[d.name for d in self.data]})")

    def table(self, i=None):
        """Text: without ``i`` the objective, the terms and one line per data set; with ``i`` (an index or a name)
        that data set's points."""
        if i is not None:
            return (self.dataset(i) if isinstance(i, str) else self.data[i]).table()
        if not self.has_data:
            return "\n".join([f"extrapolations ({len(self.data)})"] + [f"{str(d.name):<24} {str(d.quantity):<18} "
                                                                       f"{d.n:>6}" for d in self.data])
        lines = [f"objective {self.objective:.10g}   data chi2 {self.chi2:.10g}   ln P {self.log_prob:.10g}"]
        lines += [f"  {name:<40} {value:.10g}" for name, value in self.terms.items() if value != 0.0]
        lines.append(f"{'data set':<24} {'quantity':<18} {'points':>6} {'chi2':>12} {'norm':>10} {'norm chi2':>10}")
        for d in self.data:
            lines.append(f"{str(d.name):<24} {str(d.quantity):<18} {d.n:>6} {d.chi2:>12.6g} {d.norm:>10.6g} "
                         f"{d.norm_chi2:>10.4g}")
        return "\n".join(lines)

    def describe(self):
        """Plain data: objective, chi2, terms, and per data set its name, quantity, unit, points, chi2, extra names,
        its rows of the objective's rows (``[first, count]``)."""
        return {"objective": self.objective, "chi2": self.chi2, "terms": dict(self.terms), "log_prob": self.log_prob,
                "has_data": self.has_data,
                "data": [{"name": d.name, "quantity": d.quantity, "unit": d.unit, "frame": d.frame, "points": d.n,
                          "chi2": d.chi2, "norm": d.norm, "extra": sorted(d.extra),
                          "rows": None if d.rows is None else [d.rows.start, d.rows.stop - d.rows.start]}
                         for d in self.data]}

    def __repr__(self):
        if not self.has_data:
            return f"Evaluation({len(self.data)} extrapolations)"
        return f"Evaluation(objective={self.objective:.10g}, chi2={self.chi2:.10g}, {len(self.data)} data sets)"


def _legendre_rows(flat):
    """The engine's frame of Legendre coefficients of a segment (per point: a count, then that many coefficients)
    as one array per point."""
    rows, i = [], 0
    flat = np.asarray(flat, float)
    while i < flat.size:
        count = int(round(flat[i]))
        rows.append(flat[i + 1:i + 1 + count].copy())
        i += 1 + count
    return rows


def _forward_pass(engine, x):
    """The extrapolations at x from the engine's forward pass (what ``AZUREOut_*`` of extrapolation mode prints):
    per grid the model value, c.m. energy, angle and excitation energy of each point, its S-factor conversion and,
    for an angular distribution, its Legendre coefficients.  The engine's objective evaluation needs the data and
    is not run."""
    n = int(engine.update_segments_rwa(x))
    data = []
    for i in range(n):
        extra = {"sfactor_conversion": np.asarray(engine.calculated_conv(i), float)}
        legendre = _legendre_rows(engine.calculated_angular_dists(i))
        if any(r.size for r in legendre):
            extra["legendre"] = legendre
        data.append({"segment": i + 1, "segment_key": -1, "energies": engine.calculated_energies(i),
                     "angles": engine.calculated_angles(i),
                     "excitation_energies": engine.calculated_excitation_energies(i),
                     "model": engine.calculated_segments(i), "values": [], "errors": [], "residuals": [],
                     "extra": extra})
    return {"data": data}


def evaluate(session, x=None):
    """One engine evaluation of ``session`` (a :class:`pyazr.Session`) at the free vector ``x`` (default its start,
    ``params_rwa``): an :class:`Evaluation`.  In extrapolation mode its data sets are the extrapolations (a vector of
    the data's layout is taken for its free R-matrix block, as the session's other calls take it) and it has no
    chi-squared."""
    has_data = session.mode == "data"
    with session_lock(session):
        x = np.asarray(session.params_rwa if x is None else x, float).ravel()
        if has_data:
            raw = session.sess.evaluation_rwa(x)
        else:
            x = np.asarray(session._mode_free_vector(x), float).ravel()
            raw = _forward_pass(session.sess, x)
        try:
            if has_data:
                described = session._parsed("active")
            else:
                described = [d for d in session._parsed("extrapolations") if d.active]
        except RuntimeError:  # the description counts other data sets than the engine: arrays without names
            described = None
    return Evaluation(x, raw, described, has_data)
