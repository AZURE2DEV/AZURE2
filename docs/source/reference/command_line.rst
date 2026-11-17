Command-Line Usage
==================

AZURE2 runs a project from the command line, or from Python through ``pyazr``
(:doc:`../user_guide/pyazr`). This version has no graphical setup utility
(:doc:`../user_guide/graphical_interface`).

Command-Line Mode
-----------------

.. code-block:: bash

   AZURE2 <command> [options] input_file.azr

The command names the calculation; its options give every answer the run
needs, by name (see `Commands`_ below).  AZURE2 reads nothing from standard
input: a run without a command is refused with exit status 2 and the list of
the commands.  The paths a project names -- its output and checks
directories, data, reference, grid, covariance and table files -- and the
files the options name (a parameter file, an amplitude file, a temperature
file) are relative to the directory of ``input_file.azr``, wherever AZURE2 is
started; messages print them as the project writes them.  ``--validate`` and
``--describe`` read the data files by the same rule (relative to the project
file's directory), as ``pyazr.validate`` and ``pyazr.describe`` do.
``--no-gui`` (needed by builds with the Qt setup utility) is accepted and has
no effect; ``--no-readline`` (from the time of the text menu) is refused with
exit status 2: drop it.  Scripts written for the menu: `Scripts written for the
menu`_.

Commands
^^^^^^^^

An option that the command does not take, a missing required answer, a value
out of range, a file that cannot be read (relative to the directory of the
project file, as the project's paths), and an unknown command are refused
before anything is calculated, with exit status 2 and a message in three
parts: what is wrong, what the command takes, how to fix it.  A command
option given without a command is refused the same way.  ``AZURE2 --help``
lists the commands and their options first; ``AZURE2 --describe FILE --json``
gives them as data (the member ``commands``: per command its name and its
options with kind, range, default, whether it is required and, for one of two
ways, the other way in ``instead``).  Both are made from one table
(``src/RunRequest.cpp``); the test ``tests/cli_commands`` checks that the
help, the JSON and this page name the same commands and options.

.. list-table::
   :widths: 18 82
   :header-rows: 1

   * - Command
     - What it runs, and its options
   * - ``calculate``
     - The data sets at their data points from the parameters (no fit).
       ``--covariance-band``, ``--parameter-file=FILE``,
       ``--capture-amplitudes=FILE``.
   * - ``fit``
     - Fits the parameters to the data sets, then calculates them.
       ``--covariance-band``, ``--scale-covariance``,
       ``--parameter-file=FILE``, ``--capture-amplitudes=FILE``.
   * - ``extrapolate``
     - The extrapolations of the project (no data).  ``--covariance-band``,
       ``--parameter-file=FILE``, ``--capture-amplitudes=FILE``.
   * - ``minos``
     - A fit, then the MINOS interval of every free parameter.
       ``--delta-chi2=X`` (required), ``--covariance-band``,
       ``--scale-covariance``, ``--parameter-file=FILE``,
       ``--capture-amplitudes=FILE``.
   * - ``rate``
     - The reaction rate of one entrance and exit pair.
       ``--rate-entrance=N`` and ``--rate-exit=N`` (required), and one of two
       ways to give the temperatures (GK): ``--temperature-file=FILE``, or
       ``--temperature-min=T``, ``--temperature-max=T`` and
       ``--temperature-step=T``; ``--parameter-file=FILE``.
   * - ``mcmc``
     - The ensemble sampler (a build with ``USE_MCMC``).  ``--walkers=N``,
       ``--steps=N``, ``--sampler-threads=N`` and ``--basis=rwa|physical``
       (required); ``--parameter-spread=P``, ``--energy-spread=E``,
       ``--overwrite-samples``, ``--parameter-file=FILE``,
       ``--capture-amplitudes=FILE``.

The options of the commands:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Option
     - Answer
   * - ``--delta-chi2``
     - ``=X``: the rise of chi-squared that bounds each MINOS interval
       (Minuit's Up; 1 for one standard deviation of one parameter), X > 0.
   * - ``--covariance-band``
     - Writes the cross-section uncertainty band from the parameter
       covariance (this fit's, else the ``covariance.dat`` a fit saved).
   * - ``--scale-covariance``
     - With ``--covariance-band`` after a fit: inflates the covariance by the
       reduced chi-squared when it exceeds 1 (Birge ratio); refused without
       ``--covariance-band``.
   * - ``--parameter-file``
     - ``=FILE``: start from the parameters of FILE (a ``param.sav``),
       applied by name; without it, from the ``<levels>`` of the project.
       A file none of whose names is a parameter of the run is refused
       (``E_PARAMETER_FILE_NO_MATCH``, exit 1: nothing of it would be
       applied); a file with names the project has no parameter for, or
       without names of parameters the run has, is applied with
       ``W_PARAMETER_FILE_NAMES``, which lists both sides.  Names the project
       gives to the parameters of its other run (the data sets' in a run
       without data, the extrapolations' in a run with data) and fixed
       data-set or extrapolation parameters the file leaves out are not
       reported.  A file written for a project before its conversion to
       format 3 names data-set parameters by their old names:
       ``--convert-params`` renames it.  The log says ``Reading User
       Parameter File...`` when the file was applied.
   * - ``--capture-amplitudes``
     - ``=FILE``: read the external-capture amplitudes from FILE (an
       ``intEC.dat`` of an earlier run, checked against its signature file;
       FILE is read, not rewritten) instead of computing them into
       ``intEC.dat``.  Used only when the run has external capture; otherwise
       the run says so and ignores it.
   * - ``--rate-entrance``
     - ``=N``: the rate's entrance pair key, 1 or more.
   * - ``--rate-exit``
     - ``=N``: the rate's exit pair key (at least -1, not 0; -1 the total
       radiative capture), different from the entrance.
   * - ``--temperature-file``
     - ``=FILE``: the rate's temperatures, GK, one per line; required unless
       the three options below are given.
   * - ``--temperature-min``
     - ``=T``: the lowest temperature, GK (with ``--temperature-max`` and
       ``--temperature-step``, without ``--temperature-file``).
   * - ``--temperature-max``
     - ``=T``: the highest temperature, GK.
   * - ``--temperature-step``
     - ``=T``: the temperature step, GK.
   * - ``--walkers``
     - ``=N``: the sampler's walkers, 2 to 1000.
   * - ``--steps``
     - ``=N``: the sampler's steps, 100 to 10000000.
   * - ``--parameter-spread``
     - ``=P``: the walkers' initial spread, percent of each parameter, above 0
       and up to 50; default 5.
   * - ``--energy-spread``
     - ``=E``: the walkers' initial spread of a level energy, keV, above 0 and
       up to 1; default 1.
   * - ``--sampler-threads``
     - ``=N``: the sampler's threads, 1 or more (``--threads`` sets the
       OpenMP threads of the run).
   * - ``--basis``
     - ``=rwa`` samples reduced width amplitudes, ``=physical`` the project's
       physical parameters.
   * - ``--overwrite-samples``
     - Overwrites ``samples.mcmc``; without it the sampler continues the
       file.

Examples:

.. code-block:: bash

   AZURE2 calculate project.azr
   AZURE2 fit --parameter-file=fit.sav --covariance-band project.azr
   AZURE2 extrapolate --parameter-file=output/param.sav project.azr
   AZURE2 minos --delta-chi2=1 project.azr
   AZURE2 rate --rate-entrance=1 --rate-exit=2 --temperature-min=0.01 \
          --temperature-max=2 --temperature-step=0.01 project.azr
   AZURE2 mcmc --walkers=24 --steps=1000 --sampler-threads=4 --basis=rwa project.azr

The program's other options (the calculation's, the fitters', ``--threads``
and the switches) go with any command.

.. _options-as-config-keys:

Options of the calculation in the project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The options that change what a project computes -- the formalism, the
transformation of the input parameters, the Coulomb table, the external
capture -- are keys of the project's ``<config>`` section, so that the project
file says everything a calculation depends on.  Each key has the name of the
pyazr keyword and the value ``0`` or ``1`` (``true``/``false``):

=========================  ============================================================  =======
key (= pyazr keyword)      flags (on / off)                                              default
=========================  ============================================================  =======
``use_brune``              ``--use-brune`` / ``--no-brune``                              1
``use_park``               ``--use-park``                                                0
``transform``              ``--no-transform`` (off)                                      1
``ignore_externals``       ``--ignore-externals`` / ``--no-ignore-externals``            1
``use_long_wavelength``    ``--no-long-wavelength`` (off)                                1
``ec_convection_current``  ``--ec-convection-current`` / ``--no-ec-convection-current``  1
``use_rmc``                ``--use-rmc``                                                 0
``use_coulomb_table``      ``--coulomb-table`` / ``--no-coulomb-table``                  1
``use_ec_table``           ``--ec-table`` / ``--no-ec-table``                            1
``replace_grids``          ``--grid-replacement`` / ``--no-grid-replacement``            1
=========================  ============================================================  =======

A flag given on the command line (or a keyword given to a pyazr session) wins
over the key; the key wins over the default.  A project without the keys
computes what it computed before.  For example, a project in the standard
(Lane--Thomas) parametrization states it once::

   <config>
   true
   output/
   ...
   use_brune=0
   format=3
   </config>

and ``AZURE2 project.azr`` then runs it without ``--no-brune``;
``AZURE2 --use-brune project.azr`` runs the same project with Brune's level
matrix for this run only.  ``--use-park`` also gives ``use_brune`` on, and
``--no-brune`` also gives ``use_park`` off, unless the same
command line states that option itself.  ``AZURE2 --describe`` prints each
option's value and where it came from (``default``, ``project``, ``command
line``) for every project (``--json``: the member ``options``).  Refused, by
``--validate`` and by the run, each with its code: a value other than 0, 1,
true, false (``E_CONFIG_OPTION_VALUE``); the removed option ``use_gsl_coul``
with either value (``E_CONFIG_KEY_REMOVED``: delete the line; see
``--gsl-coul`` below); and a pair that cannot hold together,
whatever gave the two values -- a flag, a key or the default
(``E_CONFIG_OPTIONS_CONTRADICT``, naming both with their sources, what is
allowed and the fix): ``use_park=1`` with ``use_brune=0`` (so ``--use-park
--no-brune`` is refused in either order), ``use_rmc=1`` with ``use_brune=1``
(Brune is on by default: ``--use-rmc`` needs ``--no-brune``, a key
``use_rmc=1`` needs ``use_brune=0``), ``ec_convection_current=0`` with
``use_long_wavelength=1`` (``--no-ec-convection-current`` needs
``--no-long-wavelength``).  A pyazr session refuses the same keywords the
same way.  pyazr sets a key with
``Project.set_config(use_brune=False)``; the converter, ``Project.new`` and
``save_fit`` write no key they did not find.

Not keys (not a property of a project): the fit's and the job's settings
(``--use-gradient``, ``--use-lm``, ``--max-calls``, ...), the table file of
``--coulomb-table=FILE``, the command, threads, and the environment
switches below.

Scripts written for the menu
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before the commands, AZURE2 asked for the calculation and its answers in a
text menu, and scripts piped the answers to it
(``printf '2\nn\nparam.sav\n\n7\n' | AZURE2 project.azr``, often with the
option ``--no-readline``, which is refused now).
The menu is gone; such a script is refused with exit status 2 and the list of
the commands, and nothing is calculated.  Replace the pipe by the command:
``python3 tools/menu_to_command.py '<the answers>'`` prints it (it walks the
answers through the menu's questions as the menu read them; whether the menu
asked for the amplitude file depended on the project, and the tool says when
that matters).  The common answers:

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Menu answers
     - Command
   * - ``1\n\n\n7\n``
     - ``AZURE2 calculate project.azr``
   * - ``1\nfit.sav\n\n7\n``
     - ``AZURE2 calculate --parameter-file=fit.sav project.azr``
   * - ``2\nn\n\n\n7\n``
     - ``AZURE2 fit project.azr``
   * - ``2\ny\ny\nfit.sav\n\n7\n``
     - ``AZURE2 fit --covariance-band --scale-covariance
       --parameter-file=fit.sav project.azr``
   * - ``3\nn\nfit.sav\n\n7\n``
     - ``AZURE2 extrapolate --parameter-file=fit.sav project.azr``
   * - ``4\n1\nn\n\n\n7\n``
     - ``AZURE2 minos --delta-chi2=1 project.azr``
   * - ``5\n\n1\n-1\nyes\ntemps.dat\n``
     - ``AZURE2 rate --rate-entrance=1 --rate-exit=-1
       --temperature-file=temps.dat project.azr``
   * - ``5\n\n1\n2\nno\n0.01\n2\n0.01\n``
     - ``AZURE2 rate --rate-entrance=1 --rate-exit=2 --temperature-min=0.01
       --temperature-max=2 --temperature-step=0.01 project.azr``
   * - ``6\n\n\n24\n1000\n\n\n4\nyes\nyes\n``
     - ``AZURE2 mcmc --walkers=24 --steps=1000 --sampler-threads=4
       --basis=rwa --overwrite-samples project.azr``

The second blank line of the calculate, fit and extrapolate answers was the
external-capture amplitude file (blank: compute them), asked only of a
project with external capture; ``--capture-amplitudes=FILE`` gives a file.
A script that left out the band line of modes 2 and 3 gave the
parameter-file name to the y/n question: such a run used the ``<levels>`` of
the project (before November 2026 without an error).  The log of an older run
says which: ``Reading User Parameter File...`` when it read the file,
``Creating New param.par File...`` when it did not.

Available Options
-----------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Flag
     - Description
   * - ``--help``
     - List all available command-line options. This is the only flag that
       can be used without specifying an Input file.
   * - ``--no-gui``
     - Accepted and ignored: this build has no graphical setup utility. Builds
       with the Qt setup utility need it to run in text mode.
   * - ``--use-brune``
     - Use the Brune parameterization (equivalent to the GUI's
       "Use Brune formalism" option).  This is the default.
   * - ``--no-brune``
     - Use the standard Lane-Thomas parameterization (constant boundary
       conditions) instead of Brune's or Park's.  The input parameters are
       transformed to formal ones, with the boundary condition of each J^pi
       group set at its lowest level and each formal level oriented like its
       observed level (:doc:`../theory/conventions`); a project whose observed
       widths have no standard equivalent is refused with the J^pi groups
       concerned (see :doc:`troubleshooting`).  Builds before 2026-11-07 set
       the boundary condition at the first level listed; a note at start-up
       names each group where that differs.
   * - ``--use-park``
     - Use Park's level-dependent boundary conditions (Phys. Rev. C 104,
       064612) in place of Brune's alternative level matrix.  The fit
       parameters are then the observed reduced width amplitudes,
       :math:`\Gamma_c = 2P_c\gamma_c^2`.  Same cross sections as the
       default; see the note below.
   * - ``--ignore-externals``
     - A photon channel whose width is given as exactly zero (a 0 in
       ``<levels>``, or a channel the level does not list) has no external
       part: no external width and no channel capture through that level.
       This is the default.  The test is on the input, in the observed
       (Brune) form, so it is the same in every basis and does not change
       during a fit; builds before 2026-11-08 tested the current amplitude
       against 10⁻⁸ MeV^½, so the basis in use and a fit crossing 10⁻⁸
       changed it (``docs/upgrade/external_capture.md`` section 3).
   * - ``--no-ignore-externals``
     - Keep the external width of a level whose input width is zero (the
       GUI option unticked; pyazr ``ignore_externals=False``).  Needed for the
       R-matrix route, see below.  Before November 2026 the CLI had no way to
       turn the option off.
   * - ``--use-rmc``
     - Use the RMC capture formalism for neutron capture (equivalent to the
       GUI option).  Needs ``--no-brune`` (or ``use_brune=0`` in the
       project): with Brune on, the run is refused.
   * - ``--gsl-coul``
     - Removed (2026-11-16): the run is refused (exit 2) with the reason.
       It chose GSL's Coulomb functions, whose :math:`P` and :math:`S` are
       wrong by up to 5 % inside the outer turning point, below the Coulomb
       barrier (see :doc:`../developer/coulomb_table`).  AZURE2's own
       functions, which every run without the flag already used, are the only
       ones.  The ``<config>`` key ``use_gsl_coul`` is refused in the same
       words (``E_CONFIG_KEY_REMOVED``: delete the line), and so is the
       keyword of ``pyazr.Session`` (``TypeError``).
   * - ``--coulomb-table[=FILE]``
     - Take :math:`P`, :math:`S` and the hard-sphere phase of charged
       channels, open and closed, from the universal Coulomb table, with the
       exact routines outside it (the default).  ``FILE`` defaults to
       ``$AZURE2_COULOMB_TABLE``, else ``~/.azure2/coulomb_table.bin``; build
       it once with ``--build-coulomb-table`` (about 5 min on four
       threads, 88 MB), or with pyazr's ``pyazr coulomb-table build``, or
       download the published copy with ``pyazr coulomb-table fetch``
       (:doc:`../getting_started/installation`).  Without a table file the run says so once and uses
       the exact routines.  The natural-neon example runs about four times faster.
       Pairs with a nuclear potential and neutral pairs are not tabulated
       (see :doc:`../developer/coulomb_table`).
   * - ``--no-coulomb-table``
     - Exact Coulomb functions everywhere.
   * - ``--threads=N``
     - The OpenMP threads of the process (``N`` >= 1), over
       ``OMP_NUM_THREADS``; without it, ``OMP_NUM_THREADS`` or the OpenMP
       runtime's default decide.  The run prints the count in force and where
       it came from after its banner (``Threads: 2 (--threads, over
       OMP_NUM_THREADS=4)``), ``--describe`` in its last block.  Results depend
       on the count only as they depend on ``OMP_NUM_THREADS``.
       ``--build-coulomb-table`` and ``azure2-coulomb-table`` (whose
       ``--threads=N`` means the same) take it.  The sampler still asks for its
       own thread count.
   * - ``--build-coulomb-table[=FILE]``
     - Build the Coulomb table with the default specification on all OpenMP
       threads, write it to ``FILE`` (default as for ``--coulomb-table``),
       report progress in steps of 10 % and exit; no project file is needed.
       The same build as **Configure → Build Coulomb Table...** in the GUI.
   * - ``--convert``
     - ``--convert OLD NEW [--json-report R] [--answers A] [--exact-grids] [--strict]
       [--allow-questions]``: write the old project ``OLD`` (format 1 or 2) in project format 3
       (:doc:`format3`; the design is ``docs/upgrade/project_format3.md`` section 8) to ``NEW``, which
       must not exist, and exit; this AZURE2 reads format 3 only (:doc:`../user_guide/older_projects`).
       Every data set, extrapolation, target, beam, detector, shared normalization and energy scale
       becomes a named record.  Where the old file left a meaning implicit (a thick target's numbers
       per incident particle, a Gaussian convolution as a beam spread, an effect that also folded the
       extrapolation of the same index), the record states it with a ``# CONVERTED:`` comment and a
       note; values the old engine did not use (a group member's own norm or shift) are dropped with a
       note; where the old meaning cannot be written (two pairs with the same particles and states) the
       converter asks.  ``R``: the report as JSON (notes, questions, the parameter names old -> new
       such as ``segment_3_norm`` -> ``s03.norm``); ``A``: answers, a JSON object
       ``{"<record>.<field>": value}``; ``--exact-grids``: the extrapolation energies as the list the
       old engine accumulated instead of ``range lo hi step s`` with E_k = lo + k s; ``--strict``: ask
       instead of noting an implicit meaning; ``--allow-questions``: write ``?`` for an unanswered
       question (the file then does not validate).  A prior-covariance file is written beside ``NEW``
       as ``<name>.f3<ext>`` with renamed parameters; data files are never written; ``OLD`` is kept
       unchanged and nothing existing is overwritten.  The output is validated before it is written.
       ``--to 3`` is accepted and changes nothing.  ``pyazr.convert_azr(OLD, NEW)`` is the same
       converter.
   * - ``--convert-params``
     - ``--convert-params OLD.sav NEW.sav --project OLD.azr``: the parameter file ``OLD.sav``
       (``param.sav``, ``param.par`` or ``param.fit``) written for the old project ``OLD.azr``, with
       the names that ``OLD.azr``'s conversion gives its parameters (the report's map:
       ``segment_3_norm`` -> ``s03.norm``, ``density_1`` -> ``<target>.areal_density``), written to
       ``NEW.sav``, which must not exist; the renamed names are listed.  Level parameters keep their
       names; every other byte is kept.  ``--project`` is the old project (the converter keeps it):
       the map is that of its conversion, which a format-3 project alone does not determine; a
       format-3 project there is refused.  ``pyazr.convert_params(OLD, NEW, project)`` does the same.
   * - ``--validate``
     - ``--validate FILE [--json]``: read the format-3 project ``FILE``, check every record
       against the table (:doc:`format3`) and print every problem, one per line as
       ``file:line: record.field: severity code: message``, or with ``--json`` one JSON array of
       diagnostics (code, line, column, section, record, field, found, allowed, message, and a
       syntax-only fix).  The checks that need ``<levels>`` or the data files are made too (data
       files are read relative to the working directory, as the engine opens them); when the
       project has an error they are not, and their number is printed as a note.  Exit 0 without
       errors, 2 with errors,
       3 for a file that is not a format-3 project (with the conversion command), 1 when ``FILE``
       cannot be read.  ``pyazr.validate(path)`` returns the same diagnostics.
   * - ``--print-predictions``
     - ``--print-predictions FILE``: set the project ``FILE`` up as a calculation would (its
       data files are read relative to the project file's directory), print what every data set and
       then every extrapolation predicts, one line each as ``data <label>: <tree>`` or
       ``extrapolation <label>: <tree>``, and exit; the label is the data set's or extrapolation's
       name.  The tree is the expression the model pass
       evaluates at each point: leaves ``sigma(1->2)``, ``sigma(1->2; mult=E1)``,
       ``dsdo_cm(1->2)``, ``ay``, ``phase(1->1; J=..; l=..; s=..)``, ``ak``/``Ak``,
       ``angdist(..; k=K)``, ``sigma_total(1->1)`` (``thm_`` before a THM observable;
       ``angle=`` for a term at a fixed lab angle); nodes ``sum(a, b, ...)``,
       ``c*J0/Jk*a`` (a SUM term with its scaling c and the frame factor J_0/J_k into the
       frame of the first term), ``abundance(<sample>, <nucleus>)*a``, ``ratio(a, b)*J1/J0``
       (each term with its own lab-to-c.m. factor J_k), ``transmission(n=<n>, a)`` and
       ``yield(n=<n>, a, s_t)`` (n the nominal sample thickness, atoms/b).  Exit 0, or 1 with
       the set-up's messages on standard error when the project cannot be set up.
       ``pyazr.predictions(path)`` returns the same lines.
   * - ``--describe``
     - ``--describe FILE [--json] [--record NAME]``: print what the engine makes of the format-3 project
       ``FILE``: for every data set the quantity, the observable computed, the reaction tree
       with the ``<levels>`` pair key each term resolved to, the windows, the data file as the
       engine opens it (its resolved path), the rows it holds and the rows kept inside the
       windows with their lab-energy range (JSON ``rows_read``, ``rows_kept``,
       ``kept_E_lab_MeV``), the target, beam and detector effects linked to it
       and its fit parameters with names and kinds; then the extrapolations and the model data.
       The place to see how a file was understood.  A project with errors prints them (exit 2).
       With ``--record NAME`` only the data set, extrapolation, datum or potential ``NAME`` is
       described (exit 2 when the project has no such record).  A data set without an
       ``energy_window`` reads "window E_lab all energies (no energy_window)".
       ``pyazr.describe(path, as_json=False)`` returns the same; ``Project.describe(name)`` one record.
       After the records come the options of the calculation, each with its
       value and source, for every project.  The text form ends with the
       settings of this run and machine: the threads in force and every switch
       of kind *numbers* with its value and source, then the other switches
       that are set (``NAME=value (kind)``, see `Environment switches`_);
       ``--json`` has them as the members ``threads`` and ``switches``.
       ``pyazr.describe`` leaves this last block out.
   * - ``--schema``
     - Print the table of project format 3 as one JSON document and exit: the
       record kinds and their fields with value types and conditions, the
       quantities with their units, column roles, reaction rules and refusals,
       the units with their factors, and the diagnostic codes
       (:doc:`format3`).  No project file is needed.  The same document is
       ``pyazr.schema()``; ``tests/schema/schema_v3.json`` is its committed copy.
   * - ``--schema-out=FILE``
     - Write the text ``--schema`` prints (byte for byte) to ``FILE``, all or
       nothing: the program writes a temporary file ``FILE.tmp`` beside it and
       renames it over ``FILE`` only after the whole text was written and
       flushed without error.  On any error (an unwritable directory, a full
       disk, ``FILE`` a directory) ``FILE`` is untouched, the message names the
       cause on standard output and the exit status is 1; ``FILE`` empty is
       refused with exit 2.  Use it, not ``--schema > FILE``, to regenerate
       ``tests/schema/schema_v3.json``: a redirect empties the file before the
       program runs.  ``--schema`` itself is unchanged.
   * - ``--no-long-wavelength``
     - Electric external capture without the long-wavelength approximation:
       the radial operator :math:`\bar e_L r^L` becomes the charge term of the
       extended Siegert theorem (Friar and Fallieros, Phys. Rev. C 29, 1645
       (1984)),
       :math:`L(2L+1)!!\,k_\gamma^{-L}\,[Z_1 g_L(k_\gamma r\,m_2/M) +
       (-1)^L Z_2 g_L(k_\gamma r\,m_1/M)]`,
       :math:`g_L(y) = \int_0^y j_L(t)/t\,dt`, which reduces to
       :math:`\bar e_L r^L` term by term as :math:`k_\gamma \to 0`, plus the
       convection-current term of the relative motion (the tangential part of
       the transverse field; see :doc:`../theory/external_capture_operator`).
       Small for ordinary E1 and E2 (0.1--0.3 % for p + :sup:`16`\ O and
       :math:`\alpha` + :sup:`12`\ C E2), large where the long-wavelength E1
       charge nearly vanishes (:math:`\alpha` + :sup:`12`\ C E1: 3.5--4.7 times
       the amplitude, and +10 % from the convection current).  M1 terms are
       the same in both modes.  Internal (resonant) photon channels are not
       affected.  pyazr: ``use_long_wavelength=False``.  Checked against an
       independent calculation in ``tests/ec_m1_exact``.
   * - ``--no-ec-convection-current``
     - With ``--no-long-wavelength``: the charge term alone, without the
       convection-current term (the exact form before 2026-10-23).  pyazr:
       ``ec_convection_current=False``.
   * - ``--no-ec-table``
     - Compute every external-capture integral directly, at every data point
       and sub-point, and interpolate the amplitudes of shifted points
       linearly between the segment's own points (flat beyond its ends), as
       before the energy tables.  By default each capture term is tabulated in
       energy once, at nodes converged to 1e-10, and every point -- shifted or
       not, sub-point or rate energy -- is interpolated from the table
       (8-point Lagrange, 50 nodes per decade, error below 1e-9; see
       :doc:`../developer/ec_table`).  ``AZURE2_EC_TABLE_CHECK=1`` also computes
       every tabulated integral directly and reports the largest difference;
       ``AZURE2_EC_TABLE_NODES=N`` sets the node density.
   * - ``--no-grid-replacement``
     - Keep the adaptive sub-point grids of experimental effects
       (``<targetInt>``) where the input file's levels put them.  By default
       they are laid out again on every evaluation for the current levels
       and energy shifts, so a fit integrates each narrow resonance on a grid
       built for where the resonance is, and the result does not depend on
       which file the session started from.  Segments with external capture
       keep their grids either way.
   * - ``--no-transform``
     - Do not perform parameter transformations; input formal R-matrix
       parameters directly (equivalent to the GUI option).
   * - ``--ec-table``
     - The default, stated explicitly: undoes an earlier ``--no-ec-table``
       (options are applied in order).
   * - ``--ec-convection-current``
     - The default, stated explicitly: undoes an earlier
       ``--no-ec-convection-current``.
   * - ``--grid-replacement``
     - The default, stated explicitly: undoes an earlier
       ``--no-grid-replacement``.
   * - ``--use-gradient``
     - Fit with MIGRAD using the analytic gradient instead of numerical
       derivatives.
   * - ``--max-calls=N``
     - Stop a MIGRAD or NLopt fit after about N function calls (MIGRAD checks
       the limit between its steps, so it overshoots a little); the fit ends at
       the lowest objective it evaluated.  LM and GSL-LM ignore the limit and
       say so (a ``Note:`` line, ``call_limit N (--max-calls) ignored``).  Default: no
       limit, MIGRAD's EDM criterion ends the fit.  ``0`` also means no limit.
       How the fit ended is printed and written to ``fitStatus.out``.  Before
       2026-10-31 the limit was a fixed 50,000 and reaching it was silent.
       ``N`` is a whole number from 0 to 10\ :sup:`9`; anything else, or the
       flag without ``=N``, stops the run (exit 2).
   * - ``--refreeze-covariance[=N]``
     - For data covariances with ``relvector`` terms, whose model is frozen at
       the first evaluation of a run (the fit's starting point): rerun MIGRAD
       up to N times (default 5) with the model re-frozen at the fitted
       parameters, until the objective stops changing.  Without the option the
       model is still re-frozen once after the fit, so the outputs match a
       calculation from ``param.sav``, and a warning is printed when the fit
       and the re-frozen objective differ.  ``N`` is a whole number from 0 to
       1000; anything else stops the run (exit 2).
   * - ``--use-lm``
     - Fit with the Levenberg-Marquardt minimizer on the analytic Jacobian;
       falls back to MIGRAD.  Ignores ``--max-calls`` and says so.
   * - ``--use-gsl-lm``
     - Fit with the GSL trust-region least-squares minimizer with geodesic
       acceleration (analytic Jacobian); falls back to MIGRAD.  Ignores
       ``--max-calls`` and says so.
   * - ``--use-nlopt``
     - Fit with NLopt's Subplex (SBPLX) instead of Minuit2.  Only in builds
       configured with ``-DUSE_NLOPT=ON``; a build without NLopt refuses the
       switch (exit 2) and names the switches it has.  The fitters and their
       switches: :doc:`../user_guide/fitting`, "The Engine's Fitters"; the
       usage text (``--help``) lists the fitters of the build from that
       table.
   * - ``--covariance-band``
     - Write the cross-section uncertainty band from the parameter covariance.
       Without the flag, calculate and fit runs ask for it interactively; with
       it they do not, so scripted runs are deterministic.
   * - ``--scale-covariance``
     - With ``--covariance-band`` after a fit: inflate the covariance by the
       reduced :math:`\chi^2` when it exceeds 1 (Birge ratio).  Never deflates.
   * - ``--use-api``
     - Removed.  AZURE2 no longer serves a socket API and stops with an error;
       drive it from Python with ``from pyazr import Session``.

``AZURE2 --help`` prints this list as the binary was built
(``tests/cli_help/check.sh`` keeps the two, and this table, in step).

Options can also be read from a file, one per line, named by
``--options-file=FILE`` or the environment variable ``AZURE_OPTIONS_FILE``;
the flags of the file are applied after those on the command line (for an
option of the calculation the file's flag wins), while for a switch option
(`Environment switches`_) the command line wins over the file.  An unknown option (a mistyped ``--use-gradiant``) stops the run before
anything is read: ``ERROR: Unknown option --use-gradiant.``, the usage line, exit
status 2.  A flag that
takes a value (``--max-calls=N``, ``--refreeze-covariance=N``,
``--coulomb-table=FILE``, ``--build-coulomb-table=FILE``) given without its
value -- ``--max-calls`` alone, or ``--max-calls 50`` with a space -- or with a
value that does not parse stops the run before anything is read: ``ERROR:``
names the flag and the form it needs, the usage line follows, exit status 2
(``tests/cli_option_values/check.sh``).  Before November 2026 ``--max-calls 50``
was an unknown option and the fit ran without a limit, ``--max-calls=x`` printed
a warning and ``--refreeze-covariance=x`` was read as 0.

Multiple options can be combined:

.. code-block:: bash

   AZURE2 --no-gui --use-brune --ignore-externals input_file.azr

A-matrix or R-matrix route
--------------------------

The first line of the ``<config>`` block chooses how the cross sections are
computed: ``true`` the A-matrix (level matrix), ``false`` the R-matrix (channel
matrix :math:`[1-RL]^{-1}R`).  The two agree to about :math:`10^{-13}`, with
Brune's parametrization or without.  The R-matrix route implements neither
``--ignore-externals``, which is on by default, nor ``--use-park``, so with
either on an R-matrix request runs the A-matrix and says so::

   WARNING: the R-matrix route (<config> first line false) does not implement
   ignore-externals, which is on (...).  The A-matrix route will be used.

Add ``--no-ignore-externals`` (and leave out ``--use-park``) to run the
R-matrix route; the run log then reads "Performing R-Matrix Calculation...".
pyazr and the GUI's calculations follow the same rule (pyazr's
``ignore_externals`` is ``True`` by default).  Before November 2026 the CLI
forced the A-matrix whenever Brune or ignore-externals was on
("--ignore-externals and --use-brune options require A-Matrix") and so never
ran the R-matrix route, while pyazr and the GUI ran it with Park's
parametrization (reading the Park amplitudes as Brune's: ``tests/13N`` up to a
factor 9 off) and pyazr also with ``ignore_externals=True`` (keeping the terms
that option drops: up to 4.5 % off).

Exit status
-----------

AZURE2 exits with 0 when the run completed and 1 when it was aborted or
refused ("Calculation was aborted.", a project that AZURE2 refuses with
``ERROR: ...`` -- e.g. a malformed ``<thm>`` block or a ``<thm>`` weight that
does not fit the data --, a failed fill, initialization or parameter
transformation, a failed MCMC run), so that a script or a job chain can tell.
A flag that takes a value given without it or with a value that does not
parse, or a bad ``AZURE_REFINE_QUADRATURE``, exits with 2 before the project
is read (``ERROR: ...`` and the usage line).  A refused command line (an
unknown command, an option the command does not take, a missing or bad
answer, a file an option names that cannot be read) exits with 2 before
anything is calculated (`Commands`_).
A format-3 project refused when it is read, or in a run with data because a
data file cannot be read (``E_DATA_FILE``), prints its diagnostics and then,
as its last error line, ``ERROR: <file> (format 3) is refused; AZURE2
--validate <file> lists the same diagnostics.``  A file that is not a format-3
project (format 1 or 2, a format-3 header over an old section such as
``<segmentsData>``, an empty file, a file without ``<config>``) is refused
before anything of it is read, exit 3, with one message as the run's last
lines (:doc:`../user_guide/older_projects`)::

   ERROR: T2.azr is a format-1 project (<segmentsData> at line 53). This AZURE2 reads format 3 only.
          Convert it:  AZURE2 --convert T2.azr T2_f3.azr     (or pyazr.convert_azr("T2.azr", "T2_f3.azr"))
          The converter keeps T2.azr unchanged and reports every change it makes.

``--validate``, ``--describe`` and ``--print-predictions`` end with the same
message and exit 3.  A project file that cannot be
opened prints ``Could not open <file>.  Check that file exists.`` and exits
with -1 (255 on POSIX systems).  Before November 2026 it exited with 0 in every case on
``dev``, and with -1 (255 on POSIX systems) after a failed run on the ``thm``
branch.

Brune and Park parametrizations
-------------------------------

Both take the observed level energies as parameters and need no boundary
condition constants.  They differ only in how the reduced width amplitudes of a
level are normalized:

.. math::

   \gamma^{\rm Park}_{\lambda c} = \gamma^{\rm Brune}_{\lambda c}\,\sqrt{J_\lambda},
   \qquad
   J_\lambda = 1-\sum_c \left(\gamma^{\rm Park}_{\lambda c}\right)^2
   \left.\frac{dS_c}{dE}\right|_{E_\lambda}
   = \frac{1}{1+\sum_c \left(\gamma^{\rm Brune}_{\lambda c}\right)^2 (dS_c/dE)_{E_\lambda}} .

Park's level matrix is Brune's multiplied from both sides by
:math:`\mathrm{diag}(\sqrt{J_\lambda})`, which leaves the collision matrix
unchanged, so ``--use-park`` reproduces the default calculation
(``tests/park_formalism/check.sh``), THM (HOES) segments included: their
amplitude is bilinear in the reduced widths with level-diagonal factors only
(:doc:`../theory/thm_implementation`, "Brune and Park";
``tests/thm_park/check.sh``).  A channel entered as an amplitude
(``gammaIsRWA``, the 33rd field of a level line) holds Brune's amplitude in
either mode, so a file is the same model in both.  What changes is the
parameter space the minimizer works in:

* a Park amplitude is the observed width, :math:`\Gamma_c = 2P_c\gamma_c^2`
  (an ANC for a closed channel), with no factor depending on the other channels
  of the level.  Fixing, bounding or putting a prior on a width therefore acts
  on one parameter, and the parameters of a level are less correlated;
* every real Brune amplitude is a valid model, but Park amplitudes must satisfy
  :math:`J_\lambda > 0`; a level with :math:`J_\lambda \le 0` has no standard
  R-matrix counterpart (its widths exceed what the channel radii allow).  The
  fit objective therefore carries a penalty :math:`\sum_\lambda (J_\lambda/10^{-3})^2`
  over levels with :math:`J_\lambda < 0`, MCMC rejects such points (a chain
  started from such a point is reported, and refused if its whole starting
  ensemble is outside: :doc:`../user_guide/mcmc`), and a run that ends there
  prints a warning naming the level;
* analytic derivatives (``--use-lm``, the analytic-gradient minimizer, the
  covariance band, the gradient and Jacobian of pyazr's ``Session.objective()``) are
  available in both modes;
* the Wigner-limit bounds act on the Park amplitudes, i.e. on the observed
  :math:`\theta^2`.

Input and output files are the same in both modes: energies, partial widths and
ANCs in ``<levels>`` and ``parameters.out``.  ``param.par`` / ``param.sav`` /
``param.fit`` / ``param.mcmc`` hold the mode's own amplitudes; their first line, ``#parametrization`` (0 standard,
1 Brune, 2 Park), says which (a comment line, so ``numpy.loadtxt`` and other
positional readers are unaffected), and a file written in the other alternative mode
is converted on read (``gamma_Park = gamma_Brune sqrt(J)``), with a message.  A
file without the tag (written before it existed) is taken as Brune's.  Before
2026-11-04 ``param.fit`` (the snapshot a fit writes every 100 calls, which an
interrupted job restarts from) and ``param.mcmc`` were written without the tag,
so a ``--use-park`` fit resumed from its own ``param.fit`` had every amplitude
converted a second time; give such a file the first line
``#parametrization 2 0`` before reading it in Park mode.  The GUI
offers the mode as "Use Park parametrization" under Runtime Options.

Level-matrix solver and profiling
---------------------------------

At every energy each :math:`J^\pi` group is solved in *channel space* when it
has fewer coupled channels than levels: the bilinear forms of the level matrix
come from a :math:`c \times c` system for :math:`c` channels instead of the
:math:`n \times n` level matrix of :math:`n` levels (Woodbury identity, see
:doc:`../developer/level_matrix`).  The other groups, and the energies too close
to a pole of the formal R-matrix for the channel-space form to keep full
precision, use the level matrix as before.  Results agree with the level-matrix
solution to about :math:`10^{-12}` relative; the analytic gradient is computed
in the same way.  Two environment variables control and report this:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Variable
     - Effect
   * - ``AZURE_LEVEL_SPACE=1``
     - Solve every group with the level matrix (the route of earlier versions),
       for comparison.
   * - ``AZURE_PROFILE=1``
     - After every evaluation and gradient, write to standard error the time
       spent (summed over threads) in the matrix fill and solve, the T-matrix,
       the channel solves, the cross-section sums, the Coulomb functions and
       penetrabilities, and the stages of the gradient, with the number of
       channel-space and level-matrix solves and why a group fell back to the
       level matrix.  The set-up's stages (every point's Coulomb functions,
       penetrabilities and shifts, computed once) come first in a report of
       their own (``[profile] set-up``); a calculation (``calculate``) reports under
       ``calculation``, the evaluation after a fit under ``final evaluation``.
       No cost when unset.

``AZURE_CHANNEL_VERIFY=1`` (slow) compares every channel-space solve with the
level matrix and reports the largest difference seen so far;
``AZURE_CHANNEL_DEBUG=1`` reports the formal representation of each group
(and why a group has none).

Analytic-gradient diagnostics (``--use-gradient`` and the pyazr Jacobian):
``AZURE_GRAD_CHECK=1`` compares sampled analytic derivatives, and the residual
Jacobian, with central finite differences computed from the single-threaded
:math:`\chi^2` (``AZURE_GRAD_CHECK=exit`` stops after the check;
``AZURE_GRAD_CHECK_N``, ``AZURE_GRAD_CHECK_IDX``, ``AZURE_GRAD_CHECK_SEGS`` and
``AZURE_GRAD_CHECK_H`` choose the number of sampled parameters, which
parameters or segments, and the step); ``AZURE_GRAD_TIMING=1`` reports the time
per gradient call; ``AZURE_GRAD_DEBUG=1`` names the segment and the reason when
the adjoint gradient is unavailable and the fit falls back to finite
differences; ``AZURE_GRAD_ADJOINT=1`` forces the adjoint gradient instead of the
Jacobian path; ``AZURE_TE_WEIGHTS_CHECK=1`` and ``AZURE_TE_WEIGHTS_FD=1`` check
the target-effect combiner weights; ``AZURE_POOL_CHECK=1`` compares the pooled,
parallel :math:`\chi^2` with a fresh serial one; ``AZURE_SERIAL_CHI2=1``
evaluates every :math:`\chi^2` serially.

``AZURE_TRACE_FCN=1`` prints to standard error, for every :math:`\chi^2`
evaluation of a calculation or a fit, its number, its value and the parameter
vector, in hexadecimal floating point (``[fcn] N chi2 | p...``), so that two runs
(at different thread counts, say) can be compared value by value;
``AZURE_TRACE_FCN=fresh`` also evaluates each fit vector on a fresh copy of the
model and prints that value after ``|| fresh``, with the segments whose
:math:`\chi^2` differ from the pooled working copy's (a fresh copy per
evaluation: slow).

``AZURE_DUMP_NONFINITE=1`` (or ``=<threshold>``, default :math:`10^{10}`):
whenever a :math:`\chi^2` evaluation is non-finite or above the threshold,
append to ``nonfinite_dump.txt`` in the project file's directory (pyazr: the session's ``cwd``) the parameters that
moved most since the last gradient point, the segments with the largest or
non-finite contributions, and their offending points, including non-finite
target-effect sub-points (at most 20 dumps per run).  A non-finite
:math:`\chi^2`, or one above :math:`10^{100}`, is handed to Minuit2 as
:math:`10^{100}`.

``AZURE_POISON_RESULTS=1``: at the start of every evaluation (after the
parameters, norms and energy shifts are handed to the data) every value the
evaluation produces -- each point's and sub-point's calculated cross section,
its E1 and E2 parts, A_y, A_0, angular-distribution coefficients and sample
terms -- is set to a quiet NaN before it is calculated, so a value read before
the evaluation writes it, for instance one kept from an earlier evaluation,
appears as NaN in the output.  This holds for every buffer an evaluation writes:
the calculation's, each working copy's in a fit or an MCMC run, the Python
session's and each reaction-rate point's.  A test of the engine, not of a project: the
output files of a correct engine do not change.  In pyazr the E1 and E2 parts
(``calculated_segments_e1``, ``_e2``) of a line that is not an
angle-integrated cross section become NaN, since they are never calculated
there (0 without the variable).

``AZURE_REFINE_QUADRATURE=K``: a diagnostic of the quadrature gate
(``docs/upgrade/energy_grid.md`` section 7): every sub-point grid of an
experimental effect gets K times its points per resonance width and K times its
smooth intervals (the smooth step divided by K; the step of the geometric tails
per distance follows the points per width).  A folded value moves by the error
of the grid it was computed on; an unfolded one does not move.  K is a whole
number from 1 (the default, the variable unset) to 64; anything else, an empty
value included, stops the run with the usage line (exit 2), and in pyazr makes
opening a session raise ``ValueError``.

``AZURE_NODE_CHECK=1``: at the end of every hand-over of a parameter vector and
before every parallel point pass, every unit of every plan task (a point without
a fold, or a sub-point of a folded point) is compared bit for bit with the node
its slot names: the stage's node table (one entry per compound nucleus, route
and bits of :math:`E_x` of the units the set-up placed, with the channel rows
there) or the working copy's overflow table (the nodes the hand-over moved off
the stage table: re-placed sub-point grids, energy shifts).  The node's nucleus,
route, :math:`E_x` and channel rows must be the unit's.  The counts go to
standard error (``AZURE_NODE_CHECK: ...``); a unit off the tables, on a wrong
node or with other rows stops the run.  The output files do not change.

``AZURE_PLAN_CHECK=1``: the points of an evaluation are calculated by the
energy tasks of the data set's calculation plan, one task per compound energy
with every point and sub-point at it, then the folds of the points with
sub-points (one writer per value, checked when the plan is built).  With the
variable set, the parallel point pass runs the tasks and folds one after the
other and checks that each writes only the values of its own points and
sub-points and that no value is written by two of them (compared after each
one, or after each batch in a pass of more than 2,000); a breach stops the run
with a message naming the task and the slot.
Then every unit of the pass (a point, or a sub-point of a fold) is solved
alone -- level-matrix fill and factorisation, its T-matrix, its record -- and
compared bit for bit with the solve it shares with the other units at the same
compound energy (the same nucleus and the same bits of E_x; the bilinear forms
of every J-group's level matrix, the T-matrix of its entrance and exit pair, and
its record).  The counts go to standard error (``AZURE_PLAN_CHECK: ... units in
... energy tasks ...: 0 differ``); a unit that differs stops the run.
A test of the engine, slow, single-threaded; the output files do not change.

Environment switches
--------------------

Every environment variable AZURE2 reads, from one table in the code
(``support/src/Environment.cpp``; ``tests/consistency/repo_consistency_test.py``
keeps this list and that table the same).  The engine reads the environment
once, when the process starts (the command line) or when the first pyazr session
is made: a variable set later in the same process changes nothing.  The kind
says what setting it does: *numbers* (an output number can move -- other
physics, another input, or another route equal to rounding), *cost* (the same
bits in another time), *prints* (prints, dumps or checks; may exit), *test*
(injects a failure for a test).  The switches of kind *numbers* are copied into
the run's settings when the run starts.

Each switch of kind *numbers* is also an option of the program (the column
Option), the stated way to set it; the variable is its equivalent and keeps
working.  A flag option (``--level-space``) has the effect of the variable set to
``1``.  The value in force is the first of: the option on the command line, the
option in the options file (``--options-file=FILE`` or
``AZURE_OPTIONS_FILE``), the variable, the default.  ``--refine-quadrature=K``
and ``AZURE_REFINE_QUADRATURE`` refuse a value outside 1..64 alike (exit 2); an
option that takes a value refuses to be given without one.  ``AZURE2 --help``
lists the options.  ``AZURE2 --describe`` ends with a block ``settings of this
run and machine (value, source):``: the threads in force (``--threads=N``,
``OMP_NUM_THREADS`` or the OpenMP runtime's default) and every switch of kind
*numbers* with its value and source (``option``, ``options file``,
``environment``, ``default``); then the other switches that are set, under
``Environment switches set:``.  ``--describe --json`` has the members
``threads`` and ``switches``.

.. list-table::
   :header-rows: 1
   :widths: 25 25 10 40

   * - Variable
     - Option
     - Kind
     - What setting it does
   * - ``AZURE2_COULOMB_TABLE``
     - ``--coulomb-table=FILE``
     - numbers
     - the Coulomb table file when ``--coulomb-table`` names none (default ``~/.azure2/coulomb_table.bin``)
   * - ``AZURE2_SRIM_DATA_DIR``
     - ``--srim-data-dir=DIR``
     - numbers
     - the directory of SRIM-2013's ``SCOEF03.DAT`` and ``SNUC03.dat`` (default: the copy compiled in)
   * - ``AZURE_OPTIONS_FILE``
     - ``--options-file=FILE``
     - numbers
     - a file of command-line options, one per line, applied after those on the command line
   * - ``AZURE_REFINE_QUADRATURE``
     - ``--refine-quadrature=K``
     - numbers
     - ``K``: every sub-point grid gets K times its points (1..64; a bad value exits 2)
   * - ``AZURE2_EC_TABLE_NODES``
     - ``--ec-table-nodes=N``
     - numbers
     - nodes per decade of the external-capture tables, over ``ecTableNodesPerDecade``
   * - ``AZURE2_SWAP_COUPLING_ORDER``
     - ``--swap-coupling-order``
     - numbers
     - couples the two particle spins of a channel spin the other way round
       (a convention test, :doc:`../theory/conventions`)
   * - ``AZURE_BOUND_SHIFT_NUMERICAL``
     - ``--bound-shift-numerical``
     - numbers
     - the bound-state shift function from the numerical Whittaker derivative, not the U ratio
   * - ``AZURE_NO_NEUTRAL_CLOSED_FORM``
     - ``--no-neutral-closed-form``
     - numbers
     - the neutral bound-state shift without its closed form
   * - ``AZURE_SHIFT_DE_1EV``
     - ``--shift-derivative-1ev``
     - numbers
     - dS/dE by a 1 eV central difference of the exact S (the route before 2026-10-19)
   * - ``AZURE_LEVEL_SPACE``
     - ``--level-space``
     - numbers
     - every J group solved in level space; equal to rounding
   * - ``AZURE2_THM_FD``
     - ``--thm-by-differences``
     - numbers
     - the THM rows of the gradient and Jacobian by finite differences
   * - ``AZURE_GRAD_ADJOINT``
     - ``--gradient-adjoint``
     - numbers
     - the adjoint gradient instead of the Jacobian path; equal to rounding
   * - ``AZURE_GRAD_DENSE``
     - ``--gradient-dense``
     - numbers
     - the dense level-matrix adjoint instead of the rank-one terms; equal to rounding
   * - ``AZURE_TE_WEIGHTS_FD``
     - ``--fold-weights-by-differences``
     - numbers
     - the target-effect combiner weights by finite differences
   * - ``AZURE_SERIAL_CHI2``
     - none
     - cost
     - the points of an evaluation one by one (the same bits)
   * - ``AZURE_GRAD_CHECK``
     - none
     - prints
     - analytic derivatives against finite differences; ``exit`` stops after the check
   * - ``AZURE_GRAD_CHECK_H``
     - none
     - prints
     - the relative step of the gradient check (default 1e-6)
   * - ``AZURE_GRAD_CHECK_IDX``
     - none
     - prints
     - the parameters the gradient check differentiates (``i,j,a-b``)
   * - ``AZURE_GRAD_CHECK_N``
     - none
     - prints
     - how many parameters the gradient check samples (default 16)
   * - ``AZURE_GRAD_CHECK_ROWS``
     - none
     - prints
     - how many rows of the gradient check are printed (default 60)
   * - ``AZURE_GRAD_CHECK_SEGS``
     - none
     - prints
     - the segments whose norms and shifts the gradient check differentiates
   * - ``AZURE_SHIFT_GRAD_CHECK``
     - none
     - prints
     - the energy-shift Jacobian columns against full differences
   * - ``AZURE_TE_WEIGHTS_CHECK``
     - none
     - prints
     - ``N``: the combiner weights of N points against finite differences
   * - ``AZURE_PROFILE``
     - none
     - prints
     - per-stage times and counters on standard error
   * - ``AZURE_GRAD_TIMING``
     - none
     - prints
     - the time of every gradient call
   * - ``AZURE_GRID_TIMING``
     - none
     - prints
     - the time and the moves of every re-placement of the sub-point grids
   * - ``AZURE_PLAN_CHECK``
     - none
     - prints
     - checks the calculation plan and the node tables (stops on a difference)
   * - ``AZURE_NODE_CHECK``
     - none
     - prints
     - checks every unit's nodes and channel rows after each hand-over (stops on a difference)
   * - ``AZURE_POOL_CHECK``
     - none
     - prints
     - the pooled, parallel chi-squared against a fresh serial one (first five calls)
   * - ``AZURE_CHANNEL_VERIFY``
     - none
     - prints
     - every channel-space solve against the level-space solve (slow)
   * - ``AZURE2_EC_TABLE_CHECK``
     - none
     - prints
     - every external-capture integral also directly, with the largest difference
   * - ``AZURE_TRACE_FCN``
     - none
     - prints
     - every chi-squared evaluation of a fit; ``fresh`` also on a fresh copy
   * - ``AZURE_DUMP_NONFINITE``
     - none
     - prints
     - dumps an evaluation whose chi-squared is not finite or above the threshold
   * - ``AZURE_GRAD_DEBUG``
     - none
     - prints
     - why the gradient falls back to differences
   * - ``AZURE_CHANNEL_DEBUG``
     - none
     - prints
     - the formal representation of each J group
   * - ``AZR_DEBUG_GRID``
     - none
     - prints
     - the adaptive integration grids as they are built
   * - ``AZURE_POISON_RESULTS``
     - none
     - prints
     - reused result storage filled with NaN at the start of every evaluation
   * - ``AZURE_TEST_FAIL_E``
     - none
     - test
     - the calculation of a point at this lab energy (MeV) fails (``tests/pyazr/failed_point_test.py``)
   * - ``AZURE_TEST_FAIL_EXIT``
     - none
     - test
     - ... only for this exit pair key

Examples
--------

Run a calculation with Brune formalism:

.. code-block:: bash

   AZURE2 --no-gui --use-brune my_analysis.azr

Show all available options:

.. code-block:: bash

   AZURE2 --help
