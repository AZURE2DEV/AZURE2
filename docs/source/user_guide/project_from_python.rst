Writing a project from Python
=============================

A format-3 project (:doc:`../reference/format3`) is a set of named records: data
sets, targets, beams, detectors, normalizations, energy scales, extrapolations,
model data and potentials.  :class:`pyazr.Project` writes and edits them from
Python, for a person or an agent, in five steps: read the schema, add the
records, read the diagnostics, look at the description, calculate.  Nothing in
pyazr spells a field, a quantity or a unit: the ``add_<kind>`` methods, their
signatures and their help are generated from the engine's table
(``pyazr.schema()``, the JSON of ``AZURE2 --schema``), and every record is
checked by the engine's validator as it is added.

The example builds the worked example of a thick-target yield curve,
:math:`^{15}\mathrm{N}(p,\alpha)^{12}\mathrm{C}` counted per µC of an H⁺ beam
(charge state 1) with an efficiency × branching of 0.0182, on an infinite
:math:`^{15}\mathrm{N}` target, on the levels of the test project T2.  It is
``pyazr/examples/project_from_python.py``; ``tests/pyazr/project_recipe_test.py``
runs it and checks that its χ² equals that of the committed project it
rebuilds (``tests/units_scales/models/yield_story.azr``).

.. literalinclude:: ../../../pyazr/examples/project_from_python.py
   :language: python
   :start-after: # [recipe]
   :end-before: # [/recipe]
   :dedent: 4

What each step gives:

1. **Schema.** ``pyazr.schema()["quantities"][q]`` holds a quantity's units,
   column roles, required and allowed fields and refusals;
   ``pyazr.schema()["records"][kind]["fields"]`` every field of a record kind
   with its grammar, the condition where it applies, its units and its help.
   ``help(Project.add_dataset)`` prints the same.
2. **A project.** ``Project.new(path, nuclear_part=file)`` takes another
   format-3 file's sections that hold no records (``<config>``, ``<levels>``,
   ``<thm>``, ``<parameterSettings>``, ...); ``Project.open(path)`` opens a
   whole project.  Data files are named relative to the project's directory.
3. **Records.** One call per record; the keyword arguments are the field keys.
   A value is the text the file holds, a number, a ``(value, "unit")`` tuple, a
   :class:`pyazr.Nuisance` (``Nuisance(1, "fit", prior=(1, "3%"))`` is
   ``1 fit prior 1 3%``) or a list for a field given more than once.  A record
   that names another (a data set's ``target``) is checked against the project
   as it stands: add the target first.
4. **Diagnostics.** A record or an edit the engine refuses raises
   :class:`pyazr.FormatError`: ``.codes`` are the engine's diagnostic codes,
   ``.diagnostics`` its dicts (``code``, ``record``, ``field``, ``found``,
   ``allowed``, ``message``, a syntax-only ``fix``), the message its text.  Every
   problem is reported at once and the project is left as it was.
   ``p.validate()`` returns every problem of the project as it stands
   (``AZURE2 --validate --json``).
5. **Describe, then calculate.** ``p.describe(name)`` is ``AZURE2 --describe``
   of that record: the observable, the reaction with the pair keys of
   ``<levels>``, the units, the data file and its rows, the effects and the fit
   parameters.  ``p.session(**options)`` saves the project and returns a
   :class:`pyazr.azure2` session on it.

Edits after the first save change only the lines of the fields they name:
``p[name].set(field=value)`` (``None`` removes a field), ``p.rename(old, new)``
(the record, every reference to it and the parameter names made from it, in
``<parameterSettings>`` and in the prior-covariance file), ``p.remove(name)``
(refused while another record names it, with the referrers).  A project opened
and saved without an edit is written back byte for byte.
