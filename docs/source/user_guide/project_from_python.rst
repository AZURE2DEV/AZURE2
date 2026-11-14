Writing a project from Python
=============================

A format-3 project (:doc:`../reference/format3`) is a set of named records: data
sets, targets, beams, detectors, normalizations, energy scales, extrapolations,
model data and potentials.  ``pyazr.Project`` writes and edits them from
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
   ``pyazr.Nuisance`` (``Nuisance(1, "fit", prior=(1, "3%"))`` is
   ``1 fit prior 1 3%``) or a list for a field given more than once.  A record
   that names another (a data set's ``target``) is checked against the project
   as it stands: add the target first.
4. **Diagnostics.** A record or an edit the engine refuses raises
   ``pyazr.FormatError``: ``.codes`` are the engine's diagnostic codes,
   ``.diagnostics`` its dicts (``code``, ``record``, ``field``, ``found``,
   ``allowed``, ``message``, a syntax-only ``fix``), the message its text.  Every
   problem is reported at once and the project is left as it was.
   ``p.validate()`` returns every problem of the project as it stands
   (``AZURE2 --validate --json``).
5. **Describe, then calculate.** ``p.describe(name)`` is ``AZURE2 --describe``
   of that record: the observable, the reaction with the pair keys of
   ``<levels>``, the units, the data file and its rows, the effects and the fit
   parameters.  ``p.session(**options)`` saves the project and returns a
   ``pyazr.Session`` session on it.

Edits after the first save change only the lines of the fields they name:
``p[name].set(field=value)`` (``None`` removes a field), ``p.rename(old, new)``
(the record, every reference to it and the parameter names made from it, in
``<parameterSettings>`` and in the prior-covariance file), ``p.remove(name)``
(refused while another record names it, with the referrers).  A project opened
and saved without an edit is written back byte for byte.
A whole-line comment between repeats of a field that ``set`` replaces stays
where it was.  ``rename`` takes the parameter names from the engine (the names
it gives the project before and after, compared whole).

**Data sets from a sheet.** ``pyazr.records_from_csv(path, kind="dataset",
project=p)`` reads a CSV file whose header row holds field keys (a ``name``
column among them; one row per data set: ``name, reaction, quantity, unit, file,
columns, norm, ...``; a value with a comma, as a reaction, in quotes).  Every row
is checked by the engine against the project as it would stand; all refusals
come at once in one ``pyazr.FormatError``, each diagnostic with its
``row``.  It returns ``(kind, name, fields)`` triples for ``p.add``.

**The session's data sets.** On a format-3 project ``m.datasets``,
``m.active_datasets`` and ``m.extrapolations`` are the engine's description of
the built project (``pyazr.DescribedDatasets``: record name, quantity,
reaction with pair keys, windows, data file, rows read and kept, parameters),
and ``m.evaluate().data[i].name`` are the data sets' record names.

The nuclear part
----------------

``<levels>``, ``<config>`` and ``<thm>`` are not records; ``Project`` reads and
edits them through the engine as well (the engine reads each line, makes each
new or edited line and refuses what a run would refuse; a row or line an edit
does not name keeps its bytes, and a number replaced in a ``<levels>`` row keeps
its column).  Every edit is validated on the project as it would stand, returns
the engine's new warnings and raises ``pyazr.FormatError`` on an error,
the project unchanged.  A level is named by its label, ``"<J^pi> <E_x as
written>"``, by J^pi and E_x within 1 keV, or by the engine's ``(group,
level)``; a channel by its pair key and (l, s).

``p.levels()`` and ``p.level(name)``
    The levels in the order of their first rows (``pyazr.ProjectLevel``:
    label, J^pi, E_x in MeV, fixed, active, the engine's group and level, the
    channels with their width columns)::

        lv = p.level(("3/2-", 3.5))          # the 3/2- level within 1 keV of 3.5 MeV
        print(lv.label, lv.channel(1, 1, 0.5).width)

``p.pairs()`` and ``p.reaction(entrance, exit)``
    The pairs of ``<levels>`` (``pyazr.ProjectPair``: particles, masses in
    u, separation energy in MeV, channel radius in fm) and the reaction label
    of two pair keys, as a data set names it::

        p.add_extrapolation("S_low", reaction=p.reaction(1, 2), quantity="sigma",
                            unit="b", energies="range 0.1 1 MeV cm step 0.01")

``p.set_level(level, energy=, fixed=, active=)``
    E_x (MeV), whether it is fixed in a fit, whether the engine reads the
    level, on every row of the level::

        p.set_level("3/2- 3.5032", energy=3.51, fixed=False)

``p.set_channel(level, pair, l, s, width=, fixed=, width_is_amplitude=)``
    The width column of one channel, in the row's convention (a partial width
    in eV for an open channel or a photon, an ANC in fm^-1/2 for a closed one,
    a reduced-width amplitude in MeV^1/2 when the row's amplitude flag is set);
    ``width_is_amplitude`` sets that flag (column 33, a particle channel only;
    the value is not converted, so give the width in the new convention in the
    same call)::

        p.set_channel("3/2- 3.51", 1, 1, 0.5, width=6.0e4)
        p.set_channel("3/2- 3.51", 1, 1, 0.5, width=0.4, width_is_amplitude=True)

``p.set_pair(pair, radius=, light_mass=, heavy_mass=, separation_energy=, light_excitation=, heavy_excitation=, binding_energy=)``
    A pair's columns on every row of the pair (refit after a new radius);
    ``binding_energy`` is the THM binding energy B of the Trojan horse whose
    transferred particle is the pair's light partner (MeV, column 32, a particle
    pair only, >= 0)::

        p.set_pair(1, radius=4.0)
        p.set_pair(1, binding_energy=2.224566)

``p.add_pair(kind=, light=, heavy=, separation_energy=, radius=, ...)`` and ``p.remove_pair(pair, remove_channels=False)``
    A pair by its partners (dicts ``J, parity, Z, mass`` in u, ``g``); a photon
    pair names the nucleus in its final level as ``heavy`` and its E_x as
    ``heavy_excitation``.  ``add_pair`` returns the key.  ``<levels>`` holds a
    pair in its rows only, so a pair no level names yet is kept by a row of its
    own (inactive, level ID 0), which ``levels()`` does not list and
    ``add_level`` removes once a level names the pair.  ``remove_pair`` refuses
    while a level has a channel of the pair, naming them, unless
    ``remove_channels=True`` (the channel rows go with it; a level is never
    left without a channel); a data set that names the pair is refused by the
    validator::

        k = p.add_pair(light=dict(J=0, parity=1, Z=2, mass=4.0026), heavy=dict(J=1.5, parity=-1, Z=4, mass=9.0122),
                       separation_energy=10.6, radius=5.0)
        p.remove_pair(k)

``p.allowed_channels(J, parity, l_max=2, multipolarity_max=2, multipolarities_max=2, pairs=None)``
    The channels a level of J^pi has with the pairs (dicts ``pair, l, s``), by
    the rule of the Qt GUI's Levels tab: a particle pair every channel spin and
    every l <= l_max of the right parity; a photon pair the lowest
    ``multipolarities_max`` multipolarities L <= ``multipolarity_max`` the
    triangle rule allows; a beta pair Gamow-Teller and Fermi by the parent's
    spin.

``p.add_level(J, parity, energy, channels, fixed=True)`` and ``p.remove_level(level)``
    A level of an existing J^pi gets the rows of its group (the channels not
    given with width 0, fixed); a new J^pi one row per channel given.  It takes
    the next free level ID and goes at the end of ``<levels>``.
    ``channels="allowed"`` gives the group's channels, or for a new J^pi those
    of ``allowed_channels`` with its default limits (width 0, fixed)::

        p.add_level(1.5, -1, 7.0, [dict(pair=1, l=1, s=0.5, width=100.0)])
        p.add_level(0.5, +1, 9.0, "allowed")
        p.remove_level("3/2- 7")

``p.config()`` and ``p.set_config(**entries)``
    The ``<config>`` entries by name as the engine reads them: ``amatrix``,
    ``output_dir``, ``checks_dir``, the eight checks (``check_data``, ...) and
    the keys (``relativistic``, ``priorCovariance``, ``maxOrbitalMomentum``, ...);
    ``None`` removes a key::

        p.set_config(output_dir="run2/", relativistic=True)

``p.thm()``, ``p.set_thm(**options)``, ``p.set_thm_experiment(name, **keys)``, ``p.remove_thm_experiment(name)``
    The ``<thm>`` options and experiments (their data sets name them in their
    ``thm`` field), the block checked as the engine reads it::

        p.set_thm(kinematics="triple")
        p.set_thm_experiment("LC10", Ebeam=55)

``pyazr.Session.save_fit()`` of a format-3 project writes the fitted level
values into the rows they come from with the same engine edits, the THM
coherent backgrounds into their experiments and the free nuisances into their
records.
