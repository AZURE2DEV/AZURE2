# Writing a format-3 project from Python (`pyazr.Project`)

The recipe is `pyazr/examples/project_from_python.py` (run by `tests/pyazr/project_recipe_test.py`); the user page is
`docs/source/user_guide/project_from_python.rst`. Steps:

1. **Schema, not prose.** `pyazr.schema()["records"][kind]["fields"]` (grammar, `applies_when`, `units`, `doc`) and
   `pyazr.schema()["quantities"][q]` (units, columns, requires/allows, refusals). `help(pyazr.Project.add_dataset)`
   prints the fields. Never guess a field name or a unit token.
2. **Start.** `p = Project.new("<dir>/m.azr", nuclear_part="<format-3 file with the levels>")` (that file's
   `<config>`, `<levels>` and other record-free sections), or `Project.open(path)`. A format-1/2 file is refused with
   the conversion command: `pyazr.convert_azr(old, new, to=3)` first. Data files are relative to the project's dir.
3. **Add records**, a referenced record first (target, beam, detector, normalization, energy scale before the data
   set naming it): `p.add_target(...)`, `p.add_dataset(name, reaction=..., quantity=..., unit=..., file=...,
   columns=..., norm=..., ...)`. Values: the file's text, a number, a `(value, "unit")` tuple,
   `pyazr.Nuisance(value, "fit", unit=None, prior=(centre, "5%"))`, a list for a repeatable field.
4. **Read the refusal.** `pyazr.FormatError`: `.codes` (the engine's, e.g. `E_FIELD_REQUIRED`, `E_UNIT_NOT_ALLOWED`,
   `E_UNKNOWN_REFERENCE`, `E_DUPLICATE_NAME`), `.diagnostics` (record, field, found, allowed, message, fix). Fix the
   argument and call again; the project is unchanged after a refusal. `p.validate()` lists every problem.
5. **Describe before computing.** `p.describe(name)` (`--describe` of that record: observable, pair keys, units, rows
   read and kept, parameters). Compare it with what you meant.
6. **Calculate.** `with p.session() as m: ...` (the project is saved first; `m` is a `pyazr.azure2` session).

Edits: `p[name].set(field=value)` changes one line (`None` removes the field); `p.rename(old, new)` follows every
reference and parameter name (`<parameterSettings>`, the prior-covariance file written by `save`); `p.remove(name)` is
refused while referenced. Read fields with `p[name].field` or `p[name]["field"]`; the record's own kind is
`record_kind` (`kind` is a target's field). Fitted values back into the records: `m.save_fit(path, x,
norms="fitted")` writes norms, shift terms and isotope abundances into their nuisances.
