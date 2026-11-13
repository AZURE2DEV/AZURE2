The pyazr command line
======================

``python -m pyazr COMMAND ...`` runs what pyazr can say and do about a project from a shell: each command
is a call of one pyazr function.  ``--json`` gives one JSON document on standard output
(``{"command", "ok", "exit", "result"}``, and ``"error"`` with its parts ``what``, ``allowed``,
``fix`` when it fails); without it the answer is text.  Exit codes: 0 done, 1 the command ran and failed,
2 a wrong call or a refused input, 3 a file that is not a format-3 project (as the program's
:doc:`../reference/command_line`).  ``python -m pyazr commands --json`` gives this table as data.

This page is generated from the table of commands (``pyazr/cli.py``, ``reference_page()``);
``tests/consistency/repo_consistency_test.py`` checks that it is.

``commands``
------------

List the commands: name, arguments with their types, the pyazr function each calls, help.

.. code-block:: text

   python -m pyazr commands [--json]

Calls ``pyazr.cli.COMMANDS``.

- ``--json``: one JSON document on standard output.

``describe``
------------

What the engine makes of a project (as AZURE2 --describe), or the package and its engine.

.. code-block:: text

   python -m pyazr describe [FILE] [--record TEXT] [--json]

Calls ``pyazr.describe``.

- ``FILE``: the project; without it, the package and its engine (an existing format-3 project file (.azr)).
- ``--record``: only the data set, extrapolation, datum or potential of this name (a text).
- ``--json``: one JSON document on standard output.

``validate``
------------

Every problem of a project, as AZURE2 --validate (exit 2 when one is an error).

.. code-block:: text

   python -m pyazr validate FILE [--json]

Calls ``pyazr.validate``.

- ``FILE``: the project (an existing format-3 project file (.azr)).
- ``--json``: one JSON document on standard output.

``convert``
-----------

Write an old project (format 1 or 2) in format 3, as AZURE2 --convert; OLD is kept.

.. code-block:: text

   python -m pyazr convert OLD NEW [--overwrite] [--report FILE] [--answers FILE] [--exact-grids] [--strict] [--json]

Calls ``pyazr.convert_azr``.

- ``OLD``: the project in format 1 or 2 (an existing file).
- ``NEW``: the project in format 3 to write (a file to write, in an existing folder).
- ``--overwrite``: replace NEW, the report and the side files if they exist (an option without a value).
- ``--report``: write the converter's report (JSON) to this file (a file to write, in an existing folder).
- ``--answers``: answers to the converter's questions, {"<record>.<field>": value} (a JSON file holding one object).
- ``--exact-grids``: keep the old energy grids exactly (an option without a value).
- ``--strict``: refuse instead of noting an implicit meaning (an option without a value).
- ``--json``: one JSON document on standard output.

``schema``
----------

The table of project format 3 (as AZURE2 --schema), or one record kind or quantity of it.

.. code-block:: text

   python -m pyazr schema [KIND] [--json]

Calls ``pyazr.schema``.

- ``KIND``: a record kind or a quantity: only its entry (a text).
- ``--json``: one JSON document on standard output.

``kinds``
---------

The kinds of thing the engine and pyazr know (record, data, parameter kinds, options, physics functions, errors).

.. code-block:: text

   python -m pyazr kinds [TABLE] [--json]

Calls ``pyazr.kinds``.

- ``TABLE``: records, data, parameters, options, physics, units, dimensions or errors (a text).
- ``--json``: one JSON document on standard output.

``stopping``
------------

Stopping cross sections of a material for an ion (SRIM-2013 as ported from LunaSRIM; cite SRIM).

.. code-block:: text

   python -m pyazr stopping MATERIAL ION ENERGY [ENERGY ...] [--phase CHOICE] [--compound-correction X] [--json]

Calls ``pyazr.stopping.stopping``.

- ``MATERIAL``: the target's chemical formula (Ta2O5, C, CH4) (a text).
- ``ION``: the projectile (one of p, d, t, 3He, a).
- ``ENERGY``: lab kinetic energies of the projectile, MeV (a number).
- ``--phase``: solid or gas (default: SRIM's for the material) (one of solid, gas).
- ``--compound-correction``: SRIM's factor of Bragg's rule (a number; default 1.0).
- ``--json``: one JSON document on standard output.

``coulomb-table check``
-----------------------

Load the Coulomb table as a session does and give its path and content identity (exit 1: not readable).

.. code-block:: text

   python -m pyazr coulomb-table check [FILE] [--json]

Calls ``pyazr.described.coulomb_table``.

- ``FILE``: the table file (default: $AZURE2_COULOMB_TABLE, else ~/.azure2/coulomb_table.bin) (a text).
- ``--json``: one JSON document on standard output.

``coulomb-table build``
-----------------------

Build the Coulomb table (not in pyazr yet: the program builds it).

.. code-block:: text

   python -m pyazr coulomb-table build [FILE] [--json]

Not available yet: pyazr has no function that builds the Coulomb table (the engine's CoulombTable::BuildFile is not
bound).  Instead: the program builds it: AZURE2 --build-coulomb-table[=FILE] (about 5 min on 4 threads, 88 MB).

- ``FILE``: the table file (default: $AZURE2_COULOMB_TABLE, else ~/.azure2/coulomb_table.bin) (a text).
- ``--json``: one JSON document on standard output.

``coulomb-table fetch``
-----------------------

Download the Coulomb table (not possible yet: no published copy exists).

.. code-block:: text

   python -m pyazr coulomb-table fetch [FILE] [--json]

Not available yet: no published copy of the Coulomb table exists yet to download.  Instead: a table built here:
AZURE2 --build-coulomb-table[=FILE].

- ``FILE``: the table file (default: $AZURE2_COULOMB_TABLE, else ~/.azure2/coulomb_table.bin) (a text).
- ``--json``: one JSON document on standard output.

``ui``
------

Serve a project to the web UI on 127.0.0.1 until interrupted (--json: the URL as soon as it listens).

.. code-block:: text

   python -m pyazr ui FILE [--port N] [--threads N] [--token TEXT] [--json]

Calls ``pyazr.ui.server.main``.

- ``FILE``: the project to serve (an existing format-3 project file (.azr)).
- ``--port``: the TCP port (0: a free one above 20000) (a whole number; default 0).
- ``--threads``: the compute budget in threads (a whole number; default 1).
- ``--token``: a fixed access token (default: a random one) (a text).
- ``--json``: one JSON document on standard output.

``import exfor``
----------------

One EXFOR data set: reaction, author, year, columns with units, values.

.. code-block:: text

   python -m pyazr import exfor ID [--source CHOICE] [--plus N] [--json]

Calls ``pyazr.fetch_exfor``.

- ``ID``: the EXFOR data set (O2599004; a pointer as O2599004.1) (a text).
- ``--source``: the Web-API, the local EXFOR master files ($AZURE2_EXFOR_MASTER), or the first that answers (one of
  auto, web, local; default auto).
- ``--plus``: 0 original units, 1 computational units, 2 the universal grid (a whole number; default 0).
- ``--json``: one JSON document on standard output.

``import ensdf``
----------------

The ENSDF adopted levels of a nuclide (energies in keV, J^pi, half-lives).

.. code-block:: text

   python -m pyazr import ensdf NUCLIDE [--source CHOICE] [--json]

Calls ``pyazr.fetch_levels``.

- ``NUCLIDE``: the nuclide, as 14n or 20Ne (a text).
- ``--source``: IAEA LiveChart, the NNDC ENSDF API, or the first that answers (one of auto, livechart, nndc; default
  auto).
- ``--json``: one JSON document on standard output.

