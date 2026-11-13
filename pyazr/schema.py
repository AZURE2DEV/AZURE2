"""The table of project format 3 (docs/upgrade/project_format3.md sections 6.1, 7).

The engine holds the one table of record kinds, fields, quantities, units and diagnostic codes
(``experiment/src/Schema.cpp``, ``Quantities.cpp``, ``Units.cpp``); this module only hands it over.  pyazr, the web UI
and agents read field names, allowed values and units from it rather than spelling them.
"""

import warnings

__all__ = ["schema", "validate", "describe", "predictions"]


def schema():
    """The table of project format 3 as a dict: the JSON document that ``AZURE2 --schema`` prints.

    Top-level keys: ``schema_version`` ("3.0"), ``sections``, ``value_types``, ``records`` (each kind's fields with
    their type, grammar, conditions and allowed values), ``quantities`` (units with factors to the engine unit, column
    roles, reaction rules, required and allowed fields, norm rules, tree, EXFOR codes, refusals),
    ``general_norm_rules``, ``units``, ``dimensions``, ``column_roles``, ``lists`` and ``diagnostic_codes``.
    Needs the compiled engine module.
    """
    from .described import _engine as _engine_module
    return _engine_module().schema()


def validate(path):
    """Check the format-3 project ``path`` against the table: every problem, as ``AZURE2 --validate --json`` reports it.

    Returns a list of dicts, one per diagnostic (docs/upgrade/project_format3.md section 6.2): ``severity``
    ("error", "warning"), ``code``, ``file``, ``line``, ``column``, ``section``, ``record_kind``, ``record``,
    ``field``, ``found``, ``allowed``, ``message`` and a syntax-only ``fix``.  An empty list: the file is valid.
    The checks that need ``<levels>`` or the data files are made too (the engine's build of the project, model data
    included; data files are read relative to the working directory, as the engine opens them); when they cannot be
    made (an error of the validation, no readable ``<levels>``) their number is given as a warning.  A format-1 or
    format-2 file raises ``ValueError`` with the conversion command (section 8.4; ``AZURE2 --convert``).
    Needs the compiled engine module.
    """
    import json
    from .described import _engine as _engine_module
    with open(path, "rb") as f:
        text = f.read()  # bytes, as the engine reads a file
    refusal, diagnostics, skipped = _engine_module().validate_format3_text(text, str(path))
    if refusal:
        raise ValueError(refusal)
    if skipped:
        warnings.warn(f"{path}: {skipped} build check(s) were not run (they need a valid project and its <levels>)",
                      stacklevel=2)
    return json.loads(diagnostics)


def predictions(path):
    """What every data set and extrapolation of the project ``path`` predicts, as the engine's expression trees.

    Returns the lines ``AZURE2 --print-predictions path`` prints: ``"data <label>: <tree>"`` for each data set, then
    ``"extrapolation <label>: <tree>"``, the label the record's name in a format-3 project and the number (1, 2, ...)
    in an old one (docs/source/reference/command_line.rst gives the tree's text form).  The
    project is set up as a calculation would set it up, its data files read relative to the working directory.  A
    set-up that fails raises ``ValueError`` with the engine's messages.  Needs the compiled engine module.
    """
    from .described import _engine as _engine_module
    lines, messages = _engine_module().predictions(str(path))
    if any(line.endswith(": set-up failed") or line.startswith("project: ") for line in lines):
        raise ValueError(f"{path}: the project cannot be set up:\n{messages}")
    return list(lines)


def describe(path=None, as_json=False, record=None):
    """Without ``path``: the package and its engine as plain data (:func:`pyazr.described.describe_package`: version,
    engine build, Coulomb table, stopping and mass data, threads).  With ``path``: what the engine made of the format-3
    project ``path``, as ``AZURE2 --describe [--json]`` prints it.

    For every data set and extrapolation: the quantity, the observable computed, the reaction tree with the ``<levels>``
    pair keys each term resolved to, the windows, the data file as the engine opens it (the resolved path: relative to
    the working directory), the rows of the data file read and kept inside the windows with the lab-energy range of the
    kept rows (JSON: ``rows_read``, ``rows_kept``, ``kept_E_lab_MeV``, ``rows_error``), the target, beam and detector
    effects linked to it and its fit parameters with their names and kinds; then the model data (docs/upgrade/project_format3.md section 6.3).  Returns the text, or with
    ``as_json`` a dict (``datasets``, ``extrapolations``, ``model_data``).  A project with an error raises
    ``ValueError`` with its diagnostics; a format-1 or format-2 file raises ``ValueError`` with the conversion command.
    ``record``: only the data set, extrapolation, datum or potential of that name (``AZURE2 --describe --record``).
    Needs the compiled engine module.
    """
    import json
    from .described import _engine as _engine_module
    if path is None:
        from .described import describe_package
        return describe_package()
    with open(path, "rb") as f:
        text = f.read()
    refusal, diagnostics, description = _engine_module().describe_format3_text(text, str(path), bool(as_json), "",
                                                                                record or "")
    if refusal:
        raise ValueError(refusal)
    if not description:
        errors = [d["message"] for d in json.loads(diagnostics) if d["severity"] == "error"]
        raise ValueError(f"{path}: not described, the project has errors: " + "; ".join(errors))
    return json.loads(description) if as_json else description
