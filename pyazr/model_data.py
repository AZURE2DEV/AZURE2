"""The model data of a project: literature values of level properties against the model.

A ``<modelData>`` section holds format-3 ``datum`` records -- a resonance
strength, a resonance energy, a total or partial width, a branching ratio, an
ANC, a scattering length -- each compared with a function of the level
parameters (docs/upgrade/project_format3.md section 3.10).  The engine reads
and checks them (``AZURE2 --validate`` checks the records alone), evaluates
them from the observed energies and widths at every evaluation, and adds
((q - d)/e)^2 to its objective: the term "model data" of
``Session.evaluate(x).terms`` and ``Total-Model-Data-Chi-Squared`` of
``chiSquared.out``.

:func:`model_data` gives the model value of each datum at a parameter vector.
"""
from dataclasses import dataclass
from typing import List, Optional

import numpy as np


@dataclass(frozen=True)
class ModelDatum:
    """One ``datum`` of ``<modelData>`` at a parameter vector.

    ``model``, ``value`` and ``error`` are in ``unit``, the unit the datum is
    written in (``keV`` for ``value 2.63 keV``); ``residual`` is
    (model - value) / error and ``chi2`` the datum's share of the objective.
    A scattering length with an imaginary part also has ``model_imag``,
    ``imag``, ``imag_error`` and ``residual_imag`` (fm); otherwise they are
    ``None``.  Widths and energies are the observed ones at the observed
    level energy, as ``parameters.out`` reports them, in every basis.
    """
    name: str
    quantity: str
    model: float
    value: float
    error: float
    unit: str
    residual: float
    chi2: float
    source: str = ""
    model_imag: Optional[float] = None
    imag: Optional[float] = None
    imag_error: Optional[float] = None
    residual_imag: Optional[float] = None


def model_data(model, params=None) -> List[ModelDatum]:
    """The ``<modelData>`` data of ``model`` (a :class:`pyazr.Session`) at a
    free RWA vector (default ``model.params_rwa``), one :class:`ModelDatum`
    per datum in file order; empty without the section."""
    x = np.asarray(model.params_rwa if params is None else params, float)
    terms = model.sess.objective_terms_rwa(x)
    return [ModelDatum(**{k: (float(v) if isinstance(v, float) else v) for k, v in r.items()})
            for r in terms.get("model_data", [])]
