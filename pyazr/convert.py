"""The converter of old project files: ``AZURE2 --convert OLD NEW [--to 3]`` from Python.

pyazr reads project format 3 only (``pyazr.Project``, ``pyazr.azure2``); a format-1 or format-2 file is refused with
the engine's message (:func:`pyazr.project.require_format3`), which names this conversion.  The converter is the
engine's (``convert/``, through ``_azure2.convert_format3_text`` and ``_azure2.convert_project_text``): this module
only hands the texts over and writes the results.
"""
import json
import os
import warnings

from .described import _engine as _engine_module

__all__ = ["convert_azr"]


def convert_azr(old, new, overwrite=False, *, to=3, report=None, answers=None, exact_grids=False, strict=False,
                allow_questions=False):
    """Write the project ``old`` in the updated file format to ``new``.

    ``to=3`` (the default): project format 3, the converter of ``AZURE2 --convert OLD NEW --to 3`` (convert/,
    docs/upgrade/project_format3.md section 8): every data set, extrapolation, target, beam and detector a named
    record; implicit meanings stated with ``# CONVERTED:`` comments, unused values dropped, both listed in the report.
    ``answers``: a dict ``{"<record>.<field>": value}`` for the converter's questions; ``exact_grids``,
    ``strict``, ``allow_questions`` as the command-line options; ``report``: a path for the JSON report.  Side files
    (a prior-covariance copy ``<name>.f3<ext>``) are written beside ``new``.  Returns the report as a dict; raises
    ``ValueError`` (with the report's reason) when the project cannot be converted.  Data files are never written.

    ``to=2``: the format-2 conversion (``io::ConvertToFormat2``, ``AZURE2 --convert`` without ``--to``):
    ``<lastRun>`` and ``<mcmc>``'s run settings are removed, ``useAdaptiveGrid`` moves from ``<potential>`` to
    ``<config>``, ``<potential>`` stays only with a hybrid potential, and ``format=2`` is added to ``<config>``; every
    other byte is kept.  Returns the list of changes (empty for a file already in format 2).  pyazr does not read a
    format-2 file: convert it to format 3.

    ``new`` (and the side files and the report) are not overwritten unless ``overwrite``.
    """
    _azure2 = _engine_module()
    with open(old, "rb") as f:
        text = f.read()  # bytes, converted as bytes (the file's encoding is kept)
    if os.path.exists(new) and not overwrite:
        raise FileExistsError(f"{new} exists; pass overwrite=True to replace it")
    if to == 3:
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
    if to != 2:
        raise ValueError(f"to={to!r}: the converter writes format 2 or 3")
    converted, changes, problems = _azure2.convert_project_text(text)
    for p in problems:
        warnings.warn(f"{old}: {p}")
    with open(new, "wb") as f:
        f.write(converted)
    return list(changes)
