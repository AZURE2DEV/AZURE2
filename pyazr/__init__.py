"""pyazr: the Python interface of AZURE2, the R-matrix code for nuclear reactions.

The engine runs in this interpreter (the compiled module ``pyazr._azure2``); a project is a format-3 ``.azr`` file
(an old one: :func:`convert_azr`).  The objects of today's API:

- :class:`Session` -- a calculation session: one project set up by the engine; ``parameters`` (the fit
  vector's table, its start ``parameters.x0``), ``transform(x)`` (the physical values), data, ``objective()``,
  ``evaluate(x)``, output files, ``describe()``.
- :class:`Project` -- the project file as records: read, edit (``add_<kind>``, ``set``), ``validate()``,
  ``describe()``, ``save()``; ``Project.session()`` opens a :class:`Session` on it.
- :class:`Objective` (``Session.objective()``) -- the fit's objective as one callable for any minimiser or sampler
  (``f(x)``, ``value_and_grad``, ``residuals``, ``jacobian``, ``log_prob``, ``x0``, ``names``, ``bounds``); recipes for
  scipy, iminuit and emcee in ``pyazr/examples/fit_*.py``.
- :class:`Evaluation` (``Session.evaluate(x)``) -- every number at one x from one engine evaluation: ``objective``,
  ``chi2``, ``terms`` (by the engine's names), ``data[i]`` (model, residuals, chi2, energies, angles, ``extra``).
- :class:`Parameters` / :class:`Parameter` -- the parameter table: what each entry of the fit vector is (name, kind,
  unit, start, bounds, owner), from the engine's records (``Session.parameters``): ``of_kind``, ``of_level``,
  ``put``, ``take``, ``read``, ``write``.
- :mod:`pyazr.stopping`, :func:`whittaker_wave` -- physics functions; :mod:`pyazr.sources` -- the libraries
  (``nds``: EXFOR and ENSDF; ``masses``: AME2020; ``exfor_coulomb``, ``exfor_sfactor``: EXFOR rules).

Ask the package about itself (plain data, ``json.dumps`` takes it): :func:`describe` (version, build, Coulomb
table, data in use, threads), :func:`kinds` (record kinds, data kinds, parameter kinds, session options, bound
physics functions, error classes), :func:`schema` (the format-3 table), ``Session(...).describe()``.
Errors derive from :class:`Error` (:mod:`pyazr.errors`): :class:`FormatError`, :class:`ParameterError`,
:class:`EngineError` (also named ``AZURE2Error``), :class:`DomainError`, :class:`SessionClosed`, :class:`MissingFile`;
engine notices are :class:`EngineWarning`.  Each message says what was wrong, what is allowed, which call helps.

A session (run from the repository's root; ``tests/upgrade/T1`` is a test project)::

    import pyazr
    s = pyazr.Session("tests/upgrade/T1/T1.azr")    # read and set up the project
    x = s.parameters.x0                             # the fit vector at its start
    print(s.parameters.free.table())                # what each entry of x is
    f = s.objective()                               # the fit's objective: f(x), f.value_and_grad(x), ...
    P = s.parameters
    w = P.of_kind("width").free                     # the free reduced-width amplitudes
    x2 = P.put(x, w, P.take(x, w) * 1.1)            # a copy of x with them 10 % larger
    print(f(x), f(x2))                              # data chi-squared plus every prior, at both
    print(s.describe()["datasets"][0]["quantity"])  # what the session holds
    ev = s.evaluate(x)                              # every number at x: ev.objective, ev.chi2, ev.terms, ev.data
    print(ev.table())                               # the terms by name, one line per data set
    print(ev.data[0].model[:3])                     # the model values of the first data set
    s.close()
"""
import os as _os

# A build configured with -DPYAZR_OUTPUT_DIR=<dir> writes the engine module
# _azure2 to <dir>, not to this directory; PYAZR_ENGINE_DIR=<dir> makes the
# package find it there (ctest sets it for the tests in tests/pyazr).
if _os.environ.get("PYAZR_ENGINE_DIR"):
    __path__.append(_os.environ["PYAZR_ENGINE_DIR"])

from .errors import (Error, FormatError, ParameterError, ParameterKeyError, ParameterIndexError,
                     ParameterVectorError, EngineError, AZURE2Error, DomainError, DomainZeroDivision, SessionClosed,
                     MissingFile, EngineWarning)
from .session import Session
from .session_physics import ScatteringLength
from .objective import Objective
from .evaluate import Evaluation, EvaluatedData
from .parameters import Parameter, Parameters, NuclearPotential
from . import system
from .convert import convert_azr, convert_params
from . import sources
from . import stopping
from . import physics
from .physics import whittaker_wave
from .outputs import fitted_norms
from .sources.nds import (ExforData, ExforDataset, GroundState, Gamma,
                  Reference, search_exfor, fetch_exfor, fetch_levels,
                  fetch_gammas, fetch_ground_state, reference, resolve_doi,
                  aps_doi, sfactor_to_cross_section, cm_to_lab)
from .widths import ChannelWidth, WidthTable, weisskopf_width
from .angular import angular_distribution
from .tabulate import tabulate, TabulatedCrossSection
from .bands import (Band, load_covariance, rmatrix_columns, best_fit_params,
                    live_parameters, step_sizes, sensitivities, trimmed_model,
                    uncertainty_bands, extrapolation_bands)
from .modelavg import (Variant, Averaged, ModelAverage, model_average,
                       parameter_label, azr_parameter_slots,
                       write_averaged_azr)
from .priors import prior_from_fit, read_prior, write_prior, read_covariance_dat
from .schema import schema, validate, describe, predictions
from .described import kinds
from .grids import extrapolation_grid, refine_grid
from .model_data import ModelDatum, model_data
from .described import (DescribedDataset, DescribedDatasets, DescribedExtrapolation,
                        DescribedExtrapolations)
from .project import Project, Record, Nuisance, FormatError, records_from_csv
from .nuclear_part import ProjectLevel, ProjectChannel, ProjectPair

__all__ = ["Session", "convert_azr", "convert_params", "ScatteringLength", "whittaker_wave",
           "tabulate", "TabulatedCrossSection",
           "angular_distribution",
           "Parameter", "Parameters", "system",
           "NuclearPotential",
           "sources", "stopping", "physics",
           "fitted_norms",
           "ExforData", "ExforDataset", "GroundState", "Gamma",
           "Reference", "search_exfor", "fetch_exfor", "fetch_levels",
           "fetch_gammas", "fetch_ground_state", "reference", "resolve_doi",
           "aps_doi", "sfactor_to_cross_section", "cm_to_lab",
           "ChannelWidth", "WidthTable", "weisskopf_width",
           "Band", "load_covariance", "rmatrix_columns", "best_fit_params",
           "live_parameters", "step_sizes", "sensitivities", "trimmed_model",
           "uncertainty_bands", "extrapolation_bands",
           "Variant", "Averaged", "ModelAverage", "model_average",
           "parameter_label", "azr_parameter_slots", "write_averaged_azr",
           "prior_from_fit", "read_prior", "write_prior", "read_covariance_dat",
           "schema", "validate", "describe", "predictions", "ModelDatum", "model_data",
           "extrapolation_grid", "refine_grid",
           "Project", "Record", "Nuisance", "FormatError", "records_from_csv",
           "ProjectLevel", "ProjectChannel", "ProjectPair",
           "DescribedDataset", "DescribedDatasets", "DescribedExtrapolation", "DescribedExtrapolations",
           "kinds",
           "Error", "ParameterError", "ParameterKeyError", "ParameterIndexError", "ParameterVectorError",
           "EngineError", "AZURE2Error", "DomainError", "DomainZeroDivision", "SessionClosed", "MissingFile",
           "EngineWarning", "Objective", "Evaluation", "EvaluatedData"]
__version__ = '3.0.0.dev3'
