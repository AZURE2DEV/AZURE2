"""The numbers of a session at one parameter point (docs/upgrade/pyazr3.md section 7).

:func:`evaluate` (``Session.evaluate(x)``) runs one evaluation of the engine and returns an :class:`Evaluation`: the
objective, its terms by the engine's names, the data chi-squared, and every data set's arrays, read where the
program's output writer reads them (``chiSquared.out``, ``AZUREOut_*``).  Every number is the engine's; the objects
hold copies the caller owns, and no session.
"""
import numpy as np

from .model_data import ModelDatum
from .locking import session_lock
from .objective import require_data_mode


class EvaluatedData:
    """One data set at the evaluation's x (``Evaluation.data[i]``), in the engine's order.

    - ``name``, ``quantity``, ``unit`` (the engine's unit of ``model``, ``values``, ``errors``), ``frame``
      (``"c.m."``): from the session's description of the project; ``segment`` (the engine's segment number, from 1)
      and ``segment_key``.
    - ``chi2``: the data set's row of ``chiSquared.out`` (a data covariance by its per-set partition
      sum r_i (M^-1 r)_i); ``norm_chi2``, ``shift_chi2``: its normalization prior and its energy-shift priors (the
      three parts summed), as the objective adds them; ``norm``: the normalization in use (a profiled norm at n*).
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
        "shift_chi2": "its energy-shift priors", "norm": "the normalization in use",
        "energies": "MeV, c.m., after the shift", "angles": "deg, c.m.", "excitation_energies": "MeV",
        "model": "the model value of each point", "values": "the datum (c.m., engine unit)",
        "errors": "its error", "residuals": "(f/J - y n)/(dy n) per point", "n": "the number of points",
        "extra": "{name: array} of the point outputs that apply (kinds()['data_outputs'])",
        "table()": "the points as text",
    }

    def __init__(self, row, description):
        self.name = description.name if description is not None else None
        self.quantity = description.quantity if description is not None else None
        self.unit = description.engine_unit if description is not None else None
        self.frame = "c.m."
        self.segment = int(row["segment"])
        self.segment_key = int(row["segment_key"])
        self.chi2 = float(row["chi2"])
        self.norm_chi2 = float(row["norm_chi2"])
        self.shift_chi2 = float(row["shift_chi2"])
        self.norm = float(row["norm"])
        for key in ("energies", "angles", "excitation_energies", "model", "values", "errors", "residuals"):
            setattr(self, key, np.asarray(row[key], float))
        self.extra = {k: (np.asarray(v, float) if not isinstance(v, list) else [np.asarray(a, float) for a in v])
                      for k, v in row["extra"].items()}

    @property
    def n(self):
        """The number of points."""
        return len(self.model)

    def table(self):
        """The points as text: c.m. energy, angle, model, datum, error, residual, one line each."""
        lines = [f"{self.name}  ({self.quantity}, {self.unit}, {self.frame}; chi2 {self.chi2:.6g}, "
                 f"norm {self.norm:.6g}, {self.n} points)",
                 f"{'E_cm/MeV':>14} {'angle/deg':>10} {'model':>14} {'datum':>14} {'error':>12} {'residual':>10}"]
        for k in range(self.n):
            lines.append(f"{self.energies[k]:14.7g} {self.angles[k]:10.5g} {self.model[k]:14.7g} "
                         f"{self.values[k]:14.7g} {self.errors[k]:12.5g} {self.residuals[k]:10.4g}")
        return "\n".join(lines)

    def __repr__(self):
        return f"EvaluatedData({self.name!r}, {self.quantity}, {self.n} points, chi2={self.chi2:.6g})"


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
    """

    BUILT_BY = "Session.evaluate(x)"
    FIELDS = {
        "x": "a copy of the point", "objective": "what the fit minimises (= Objective(x))",
        "chi2": "the data chi-squared, Total-Chi-Squared of chiSquared.out",
        "terms": "{name: value} per term of the engine's objective", "data": "one EvaluatedData per data set",
        "model_data": "the <modelData> rows", "log_prob": "Objective.log_prob(x)",
        "excluded": "outside the engine's sampler's support", "table(i)": "text", "dataset(name)": "by name",
    }

    def __init__(self, x, raw, descriptions):
        self.x = np.array(x, float)
        self.objective = float(raw["objective"])
        self.chi2 = float(raw["chi2"])
        self.terms = {str(k): float(v) for k, v in zip(raw["term_names"], raw["term_values"])}
        self.log_prob = float(raw["log_prob"])
        self.excluded = bool(raw["excluded"])
        rows = list(raw["data"])
        if descriptions is not None and len(descriptions) != len(rows):
            descriptions = None  # the description and the engine disagree: the arrays keep their segment numbers
        self.data = [EvaluatedData(r, descriptions[i] if descriptions is not None else None)
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
        lines = [f"objective {self.objective:.10g}   data chi2 {self.chi2:.10g}   ln P {self.log_prob:.10g}"]
        lines += [f"  {name:<40} {value:.10g}" for name, value in self.terms.items() if value != 0.0]
        lines.append(f"{'data set':<24} {'quantity':<18} {'points':>6} {'chi2':>12} {'norm':>10} {'norm chi2':>10}")
        for d in self.data:
            lines.append(f"{str(d.name):<24} {str(d.quantity):<18} {d.n:>6} {d.chi2:>12.6g} {d.norm:>10.6g} "
                         f"{d.norm_chi2:>10.4g}")
        return "\n".join(lines)

    def describe(self):
        """Plain data: objective, chi2, terms, and per data set its name, quantity, unit, points, chi2, extra names."""
        return {"objective": self.objective, "chi2": self.chi2, "terms": dict(self.terms), "log_prob": self.log_prob,
                "data": [{"name": d.name, "quantity": d.quantity, "unit": d.unit, "frame": d.frame, "points": d.n,
                          "chi2": d.chi2, "norm": d.norm, "extra": sorted(d.extra)} for d in self.data]}

    def __repr__(self):
        return f"Evaluation(objective={self.objective:.10g}, chi2={self.chi2:.10g}, {len(self.data)} data sets)"


def evaluate(session, x=None):
    """One engine evaluation of ``session`` (a :class:`pyazr.Session` in data mode) at the free vector ``x``
    (default its start, ``params_rwa``): an :class:`Evaluation`."""
    require_data_mode(session, "evaluate()")
    x = np.asarray(session.params_rwa if x is None else x, float).ravel()
    with session_lock(session):
        raw = session.sess.evaluation_rwa(x)
    try:
        described = session._parsed("active")
    except RuntimeError:  # the description counts other data sets than the engine: arrays without names
        described = None
    return Evaluation(x, raw, described)
