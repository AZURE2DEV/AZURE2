"""The converter of old project files: ``AZURE2 --convert OLD NEW`` and ``--convert-params`` from Python.

pyazr reads project format 3 only (``pyazr.Project``, ``pyazr.Session``); a format-1 or format-2 file is refused with
the engine's message (:func:`pyazr.project.require_format3`), which names this conversion.  The converter is the
engine's (``convert/``, through ``_azure2.convert_format3_text`` and ``_azure2.convert_parameter_file``): this module
only hands the texts over and writes the results.
"""
import json
import os

from .described import _engine as _engine_module

__all__ = ["convert_azr", "convert_params"]


def convert_azr(old, new, overwrite=False, *, to=3, report=None, answers=None, exact_grids=False, strict=False,
                allow_questions=False):
    """Write the project ``old`` (format 1 or 2) in project format 3 to ``new``.

    The converter of ``AZURE2 --convert OLD NEW`` (convert/, docs/upgrade/project_format3.md section 8): every data
    set, extrapolation, target, beam and detector a named record; implicit meanings stated with ``# CONVERTED:``
    comments, unused values dropped, both listed in the report.  ``answers``: a dict ``{"<record>.<field>": value}``
    for the converter's questions; ``exact_grids``, ``strict``, ``allow_questions`` as the command-line options;
    ``report``: a path for the JSON report.  Side files (a prior-covariance copy ``<name>.f3<ext>``) are written beside
    ``new``.  Returns the report as a dict; raises ``ValueError`` (with the report's reason) when the project cannot be
    converted.  Data files are never written.  ``to`` must be 3 (the engine reads format 3 only).  A parameter file of
    the old project is renamed for the conversion by :func:`convert_params`.

    ``new`` (and the side files and the report) are not overwritten unless ``overwrite``.
    """
    if to != 3:
        raise ValueError(f"to={to!r}: the converter writes format 3 only")
    _azure2 = _engine_module()
    with open(old, "rb") as f:
        text = f.read()  # bytes, converted as bytes (the file's encoding is kept)
    if os.path.exists(new) and not overwrite:
        raise FileExistsError(f"{new} exists; pass overwrite=True to replace it")
    if os.path.abspath(old) == os.path.abspath(new):
        raise ValueError("the converter never overwrites its input")
    result = _azure2.convert_format3_text(text,
                                          os.path.dirname(os.path.abspath(old)),
                                          {str(k): str(v) for k, v in (answers or {}).items()},
                                          exact_grids=bool(exact_grids), strict=bool(strict),
                                          allow_questions=bool(allow_questions), input=str(old), output=str(new))
    rep = json.loads(result["report"])
    if report is not None:
        if os.path.exists(report) and not overwrite:
            raise FileExistsError(f"{report} exists; pass overwrite=True to replace it")
        with open(report, "w") as f:
            f.write(result["report"])
    if not result["ok"]:
        raise ValueError(f"{old} is not converted: {result['refusal']}")
    outdir = os.path.dirname(os.path.abspath(new))
    for name, data in result["side_files"]:
        path = os.path.join(outdir, name)
        if os.path.exists(path) and not overwrite:
            raise FileExistsError(f"{path} exists; pass overwrite=True to replace it")
        with open(path, "wb") as f:
            f.write(data)
    with open(new, "wb") as f:
        f.write(result["text"])
    return rep


def convert_params(old, new, project, overwrite=False):
    """Write the parameter file ``old`` (a ``param.sav``, ``param.par`` or ``param.fit`` written for the old project
    ``project``) to ``new`` with the names of the project's format-3 conversion: ``AZURE2 --convert-params OLD NEW
    --project PROJECT`` (convert::ConvertParameterFile, docs/upgrade/project_format3.md section 8.4).  The map is that
    of :func:`convert_azr` on ``project`` (``segment_3_norm`` -> ``s03.norm``); every other byte is kept.  Returns the
    list of ``(old name, new name)`` applied; raises ``ValueError`` when ``project`` is format 3 or is not converted.
    ``new`` is not overwritten unless ``overwrite``."""
    _azure2 = _engine_module()
    if os.path.exists(new) and not overwrite:
        raise FileExistsError(f"{new} exists; pass overwrite=True to replace it")
    with open(old, "rb") as f:
        text = f.read()
    with open(project, "rb") as f:
        project_text = f.read()
    refusal, renamed_text, renamed = _azure2.convert_parameter_file(text, project_text,
                                                                     os.path.dirname(os.path.abspath(project)))
    if refusal:
        raise ValueError(f"{project}: {refusal}")
    with open(new, "wb") as f:
        f.write(renamed_text)
    return [tuple(r) for r in renamed]
