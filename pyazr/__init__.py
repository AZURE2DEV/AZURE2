import os as _os

# A build configured with -DPYAZR_OUTPUT_DIR=<dir> writes the engine module
# _azure2 to <dir>, not to this directory; PYAZR_ENGINE_DIR=<dir> makes the
# package find it there (ctest sets it for the tests in tests/pyazr).
if _os.environ.get("PYAZR_ENGINE_DIR"):
    __path__.append(_os.environ["PYAZR_ENGINE_DIR"])

from .azure2 import azure2, ScatteringLength
from .parameters import (Pair, PairSet, Parameter, ParameterSet, LevelKey,
                         NuclearPotential)
from .scheme import LevelScheme, SchemeLevel, SchemeChannel
from .convert import convert_azr, convert_params
from . import masses
from . import stopping
from .physics import whittaker_wave
from .outputs import fitted_norms
# nds also defines a Level, but transform.Level owns the name here; reach the
# other one as pyazr.nds.Level.
from .nds import (ExforData, ExforDataset, GroundState, Gamma,
                  Reference, search_exfor, fetch_exfor, fetch_levels,
                  fetch_gammas, fetch_ground_state, reference, resolve_doi,
                  aps_doi, sfactor_to_cross_section, cm_to_lab)
from .widths import (ChannelWidth, WidthTable, weisskopf_width,
                     teichmann_wigner, lane_thomas_wigner)
from .transform import (Channel, Level, TransformedChannel, TransformedLevel,
                        transform_out, partial_widths, LevelWidthTransformer,
                        levels_from_azr, levels_from_scheme,
                        penetrability, shift, shift_derivative, whittaker)
from .angular import angular_distribution
from .tabulate import tabulate, TabulatedCrossSection
from .nuclei import (Nucleus, MergeReport, merge_projects, nucleus_of,
                     nucleus_of_pair)
from .bands import (Band, load_covariance, rmatrix_columns, best_fit_params,
                    live_parameters, step_sizes, sensitivities, trimmed_model,
                    uncertainty_bands, extrapolation_bands)
from .modelavg import (Variant, Averaged, ModelAverage, model_average,
                       parameter_label, azr_parameter_slots,
                       write_averaged_azr)
from .priors import prior_from_fit, read_prior, write_prior, read_covariance_dat
from .schema import schema, validate, describe, predictions
from .grids import extrapolation_grid, refine_grid
from .model_data import ModelDatum, model_data
from .described import (DescribedDataset, DescribedDatasets, DescribedExtrapolation,
                        DescribedExtrapolations)
from .project import Project, Record, Nuisance, FormatError, records_from_csv
from .nuclear_part import ProjectLevel, ProjectChannel, ProjectPair

__all__ = ["azure2", "convert_azr", "convert_params", "ScatteringLength", "whittaker_wave",
           "tabulate", "TabulatedCrossSection",
           "angular_distribution",
           "Pair", "PairSet", "Parameter", "ParameterSet", "LevelKey",
           "NuclearPotential",
           "LevelScheme", "SchemeLevel", "SchemeChannel",
           "masses", "stopping",
           "fitted_norms",
           "ExforData", "ExforDataset", "GroundState", "Gamma",
           "Reference", "search_exfor", "fetch_exfor", "fetch_levels",
           "fetch_gammas", "fetch_ground_state", "reference", "resolve_doi",
           "aps_doi", "sfactor_to_cross_section", "cm_to_lab",
           "ChannelWidth", "WidthTable", "weisskopf_width", "teichmann_wigner", "lane_thomas_wigner",
           "Channel", "Level", "TransformedChannel", "TransformedLevel",
           "transform_out", "partial_widths", "LevelWidthTransformer",
           "levels_from_azr", "levels_from_scheme",
           "penetrability", "shift", "shift_derivative", "whittaker",
           "Nucleus", "MergeReport", "merge_projects", "nucleus_of",
           "nucleus_of_pair",
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
           "DescribedDataset", "DescribedDatasets", "DescribedExtrapolation", "DescribedExtrapolations"]
__version__ = '3.0.0.dev0'
