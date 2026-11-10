"""The table of project format 3 (docs/upgrade/project_format3.md sections 6.1, 7).

The engine holds the one table of record kinds, fields, quantities, units and diagnostic codes
(``experiment/src/Schema.cpp``, ``Quantities.cpp``, ``Units.cpp``); this module only hands it over.  pyazr, the web UI
and agents read field names, allowed values and units from it rather than spelling them.
"""

import warnings

__all__ = ["schema", "validate"]


def schema():
    """The table of project format 3 as a dict: the JSON document that ``AZURE2 --schema`` prints.

    Top-level keys: ``schema_version`` ("3.0"), ``sections``, ``value_types``, ``records`` (each kind's fields with
    their type, grammar, conditions and allowed values), ``quantities`` (units with factors to the engine unit, column
    roles, reaction rules, required and allowed fields, norm rules, tree, EXFOR codes, refusals),
    ``general_norm_rules``, ``units``, ``dimensions``, ``column_roles``, ``lists`` and ``diagnostic_codes``.
    Needs the compiled engine module.
    """
    from .azrfile import _engine_module
    return _engine_module().schema()


def validate(path):
    """Check the format-3 project ``path`` against the table: every problem, as ``AZURE2 --validate --json`` reports it.

    Returns a list of dicts, one per diagnostic (docs/upgrade/project_format3.md section 6.2): ``severity``
    ("error", "warning"), ``code``, ``file``, ``line``, ``column``, ``section``, ``record_kind``, ``record``,
    ``field``, ``found``, ``allowed``, ``message`` and a syntax-only ``fix``.  An empty list: the file is valid.
    The checks that need ``<levels>`` or the data files are left to the format-3 reader; their number is given as a
    warning.  A format-1 or format-2 file raises ``ValueError`` with the conversion command (section 8.4).
    Needs the compiled engine module.
    """
    import json
    from .azrfile import _engine_module
    with open(path, "rb") as f:
        text = f.read()  # bytes, as the engine reads a file
    refusal, diagnostics, skipped = _engine_module().validate_format3_text(text, str(path))
    if refusal:
        raise ValueError(refusal)
    if skipped:
        warnings.warn(f"{path}: {skipped} check(s) that need <levels> or the data files were not run; "
                      "the format-3 reader of step 3b runs them", stacklevel=2)
    return json.loads(diagnostics)
