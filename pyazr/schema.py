"""The table of project format 3 (docs/upgrade/project_format3.md sections 6.1, 7).

The engine holds the one table of record kinds, fields, quantities, units and diagnostic codes
(``experiment/src/Schema.cpp``, ``Quantities.cpp``, ``Units.cpp``); this module only hands it over.  pyazr, the web UI
and agents read field names, allowed values and units from it rather than spelling them.
"""

__all__ = ["schema"]


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
