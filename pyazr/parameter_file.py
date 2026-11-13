"""A parameter file (``param.sav``, ``param.par``) applied to a session as the program applies its external parameter
file: ``Session(file, parameter_file=path)``.

The program's external parameter file (``AZUREParams::ReadUserParameters``) sets every parameter it names, the fixed
level energies and widths included.  A fit vector carries the free parameters only, so a session that read the file's
free values alone would compute another model wherever the file's fixed values differ from the project's (11B+alpha:
ten fixed widths up to 1.6 % apart, chi-squared 9925.03 against the program's 9868.04).  The session therefore takes
the whole file when it is made: every later call (``evaluate``, ``objective``, ``transform``, ``save_fit``) uses the
file's fixed values, and ``Session.parameters.x0`` is the file's free values.

Names are matched as the program matches them (``experiment::MatchParameterFileNames``).  The free / fixed split is the
project's (the program frees a parameter fixed at 0 that the file gives a value, and fixes a free one the file sets to
0; a value is applied either way, so a calculation agrees).  In an extrapolation run the file's entries of data-set
parameters (normalizations, energy shifts) are not applied: they belong to the data sets, not to the extrapolations.
"""

import json
import os
import warnings

import numpy as np

from .errors import ParameterError
from .parameters import read_parameter_file


def apply_parameter_file(engine, project_file, path, basis, extrapolations):
    """Put every value the parameter file ``path`` names into the engine (the binding's session object ``engine``,
    held under the session's lock by the caller).

    ``basis`` is the session's (0 standard, 1 Brune, 2 Park); a file whose ``#parametrization`` tag differs is refused
    (the program converts it; this does not).  ``extrapolations`` is True for an extrapolation run.  Raises
    :class:`pyazr.ParameterError` (a ``ValueError``) for a file of another basis or one that names no parameter of the
    run; warns for names the project lacks and for run parameters the file leaves out."""
    tag, values = read_parameter_file(path)
    if tag is not None and tag != basis:
        raise ParameterError(f"{path} holds basis {tag} amplitudes (0 standard, 1 Brune, 2 Park) but this session uses "
                             f"basis {basis}; open the session with the matching use_brune / use_park, or convert the "
                             f"file with the program.")
    from . import _azure2
    fixed = np.asarray(engine.params_fixed(), bool)
    full = np.asarray(engine.params_all_rwa(), float).copy()
    segment_field = {k["key"]: k["segment_field"] for k in _azure2.parameter_kinds()}
    records = engine.parameter_records()
    names = [engine.params_names(i) for i in range(len(records))]
    _match_names(project_file, path, list(values), names, fixed, records, segment_field)
    hit = skipped = 0
    for i, name in enumerate(names):
        if name not in values:
            continue
        if extrapolations and segment_field[records[i]["kind"]]:
            skipped += 1
            continue
        full[i] = values[name]
        hit += 1
    if skipped:
        warnings.warn(f"{path}: {skipped} entries of data-set parameters (normalizations and energy shifts) not "
                      f"applied in an extrapolation run; an extrapolation's grid is the one its record asks for.",
                      stacklevel=4)
    if hit == 0:
        raise ParameterError(f"{path}: no parameter name matches this project (a file of another project?); "
                             f"Session.parameters.names lists the names.")
    engine.set_params_all_rwa(full)


def _match_names(project_file, path, in_file, names, fixed, records, segment_field):
    """The file's names against the run's, by the engine's rule (``experiment::MatchParameterFileNames``): no name of
    the run raises ``ParameterError`` (E_PARAMETER_FILE_NO_MATCH); names the project lacks or run parameters the file
    leaves out warn (W_PARAMETER_FILE_NAMES).  Fixed data-set and extrapolation nuisances may be left out; names the
    project gives to another run's parameters (an extrapolation's in a data run) are neither applied nor unknown."""
    from . import _azure2
    optional = [names[i] for i, r in enumerate(records) if segment_field[r["kind"]] and fixed[i]]
    with open(project_file, "rb") as f:
        other = list(_azure2.parameter_names(f.read()))
    for text in _azure2.match_parameter_file_names(os.fspath(path), in_file, names, optional, other):
        d = json.loads(text)
        message = f"{d.get('code', '')}: {d.get('message', '')}"
        if d.get("severity") == "error":
            raise ParameterError(message)
        warnings.warn(message, stacklevel=5)
