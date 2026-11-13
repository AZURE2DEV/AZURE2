Extending AZURE2
================

A new feature is a row in the one table of its kind and a test, never a
special case in unrelated code. Each section below is the whole list of
places to touch for one kind of extension, with the test that fails when one
of them is forgotten. "Check N" is check N of
``tests/consistency/repo_consistency_test.py`` (it runs in under a second);
a name like ``experiment_schema`` is a CTest test. The recipes in the
repository's ``CLAUDE.md`` (and ``pyazr/CLAUDE.md`` for Python and the web UI)
carry the same lists with more detail; :doc:`contributing` says how a change
is checked and merged.

A data quantity (a kind of data set)
------------------------------------

1. One function in ``experiment/src/Quantities.cpp`` that fills its
   ``Quantity`` member by member (name, title, engine unit, units, energy
   roles, the reactions it allows, the fields it allows or requires, the tree,
   the refusals), and its entry in ``Quantities()``.
2. Its tree builder: the case in ``Schema::Build`` (``experiment/src/Build.cpp``,
   ``ReadObservable``) that turns the record into a ``DatasetSpec``
   (``experiment/include/azure2/experiment/SetupSpecs.h``). Until it has one,
   the build refuses the quantity with ``E_NOT_BUILT``.
3. If no existing tree fits: a row of ``TreeKinds()``
   (``experiment/src/Prediction.cpp``); a value the R-matrix does not compute
   is a leaf part whose value at a point is computed in ``PointInputs``
   (``src/ModelPass.cpp``) by a function of its own file
   (``src/CoulombReference.cpp`` is the pattern), and the part is added to
   ``IsReferenceLeaf`` (``experiment/include/azure2/experiment/Prediction.h``).
   A factor that depends on the point (an S factor) is a ``FactorRef`` under
   the tree kind ``Scaled``.
4. Tests: a record of it in the ``variants`` of
   ``experiment/tests/quantities_test.cpp`` (``experiment_quantities`` fails for
   a row without one), a case in ``experiment/tests/build_test.cpp``, a test
   project with an independent answer (``tests/ratio_to_rutherford/`` is the
   pattern), its lines in ``tests/upgrade/predictions.txt``
   (``tests/upgrade/predictions_test.py --write``, ``upgrade_predictions``).
5. ``AZURE2 --schema > tests/schema/schema_v3.json`` (``experiment_schema``,
   ``pyazr_schema``) and ``python3 tools/gen_schema_docs.py`` (check 18: the
   generated ``reference/format3.rst``); check 35 for a reference leaf.

An observable of a reaction at a mesh energy
--------------------------------------------

1. A kind in ``SiteObservable`` (``include/SiteObservable.h``) and its row in
   ``kSiteObservables`` (``src/SiteObservables.cpp``): the record's slots and
   the contraction with a reader's row (``Site::ValueFor``).
2. Its formula at the energy: one branch of
   ``GenMatrixFunc::CalculateBareCrossSection``; its adjoint in
   ``AMatrixFunc::PointTBar``.
3. What a reader needs per slot (angular row, Coulomb amplitude) in its row
   set (``include/RowSet.h``, built by ``WorkingCopy::BuildRowSets``).
4. Tests: the analytic gradient against finite differences
   (``AZURE_GRAD_CHECK``) and a test project of the quantity that uses it.

An experimental effect (an apparatus kernel)
--------------------------------------------

The engine is not edited: the fold (``src/PointFold.cpp``), the extent of the
mesh (``EffectSupport``, ``src/TargetSetup.cpp``) and the gradient's weights
(``FoldWeights``) read the table.

1. One file ``experiment/src/KernelGaussian.cpp``-like (``Kernel<Name>.cpp``):
   the fold as a template on the value type (``double`` gives the value,
   ``LinearForm`` the exact weights), its support, and a function returning
   its ``Kernel`` row; its parameters are ``KernelParameter`` rows (name,
   dimension, engine unit, ``fittable``).
2. Its line in ``Kernels()`` (``experiment/src/Kernels.cpp``) and the file in
   ``experiment/CMakeLists.txt``.
3. Tests: cases in ``experiment/tests/kernels_test.cpp`` against closed forms
   (``experiment_kernels``; the table's own checks run on every row), the
   derivative against central differences, a format-3 fixture in
   ``experiment/tests/data/``, the schema and its page as above, and a test
   project (``tests/rectangular_resolution/`` is the pattern).
4. A line in :doc:`data_model`, "Kernels".

A fittable parameter of a kernel needs no engine code: its
``KernelParameter`` says ``fittable = true`` and the project writes ``fit``
(``ApparatusRows()``, ``src/ApparatusParameters.cpp``).

A parameter kind
----------------

A nuisance that is not a number of an apparatus record (of a data set, of the
level scheme, of a THM experiment):

1. ``include/ParameterKind.h``: the value of ``ParamKind`` and its row of
   ``ParamKindTable()`` (key, unit, owner fields, prior and vary flags, doc);
   a ``ParameterRecord`` factory for new owner fields.
2. Create it with ``AZUREParams::Add`` in ``CNuc::FillMnParams`` or
   ``WorkingCopy::FillMnParams``; every later question asks the record, never
   the name (check 23).
3. Its place in ``BuildParamIndexMap`` (``src/AZUREGrad.cpp``); the fields
   that say what it belongs to are rows of ``ParameterFields()``
   (``api/src/pybind_records.cpp``), which ``pyazr.Parameter`` reads
   (``tests/pyazr/parameter_records_test.py`` fails for a kind no owner field
   describes).
4. Docs: the kinds table of ``reference/output_files.rst``; a project with
   the kind in ``tests/pyazr/parameter_kinds_test.py``.

A term of the objective
-----------------------

1. One ``ObjectiveTerm`` subclass in ``src/Objective.cpp``, appended to the
   term list in the ``Objective`` constructor (``Objective::Terms``); its
   gradient beside the others in ``src/AZUREGrad.cpp``, and its residual rows
   in ``AppendModelRows`` if the Levenberg-Marquardt path needs them. No term
   is added anywhere else (check 11).
2. Tests: the analytic gradient against finite differences
   (``AZURE_GRAD_CHECK``), the program and ``pyazr`` giving the same
   :math:`\chi^2` (``tests/pyazr/objective_test.py``).
3. Its total in ``Results::objective`` and from there in ``chiSquared.out``
   (``reference/output_files.rst``).

A record kind or a key of the project format
--------------------------------------------

A **key** (a field of a record kind):

1. One ``FieldDef`` row in its kind's function in ``experiment/src/Schema.cpp``
   (``DatasetKind()``, ``TargetKind()``, ...), made by the helper that says
   where it applies (``Req``, ``Opt``, ``RequiredWhen``, ``OptionalWhen``,
   ``ByQ``); a structured value is a ``FormSpec`` of ``StructuredValues()``.
   A new value type is a parser and a case in ``ParseValue``.
2. A fixture line in ``experiment/tests/data/`` that validates, and a case in
   ``experiment/tests/diagnostics_test.cpp`` for a refusal of its own; a
   ``ByQ`` field needs a value in the ``plausible`` map of
   ``experiment/tests/quantities_test.cpp``.
3. The schema (``experiment_schema``, ``pyazr_schema``) and its page
   (check 18); ``Schema::kVersion`` is raised if an existing project changes
   meaning, with the converter's support (``convert/``).

A **record kind**: a function ``<Name>Kind()`` in ``experiment/src/Schema.cpp``,
its entry in ``Schema::Kinds()`` and its section in ``Schema::Sections()``;
its build in ``Schema::Build``; a fixture, the schema and the page as for a
key. A **diagnostic code** is a constant in
``experiment/include/azure2/experiment/Diagnostics.h``, its row in
``DiagnosticCodes()`` and a case that triggers it alone in
``experiment/tests/diagnostics_test.cpp`` (check 24). A **unit** is a row of
``experiment/src/Units.cpp`` (``experiment_parsers``,
``experiment_quantities``).

``pyazr.Project`` edits every record through the binding of the engine's
grammar and needs no Python change; an unedited project is written back byte
for byte (``tests/pyazr/project_test.py``).

A session option (an option of the calculation)
-----------------------------------------------

A switch a project's result depends on is a ``<config>`` key, a program flag
and a ``pyazr`` keyword under one name.

1. One ``Option`` row of ``io::ConfigKeys()`` (``io/src/ConfigSection.cpp``);
   the program's parser and ``--help`` (``ApplyOptionFlag``,
   ``printOptionFlags``), ``--describe``, ``pyazr.kinds()["options"]`` and the
   web UI read it there (check 1; ``io_config_section``). A flag that implies
   another, or two that contradict, is a row of
   ``ConfigOptionImplications()`` or ``ConfigOptionContradictions()``.
2. Its place in ``Config``: a row of ``OptionStorages()``
   (``src/ProjectOptions.cpp``; ``upgrade_project_options`` fails for a row
   without one or a default unlike ``Config::Reset``).
3. Its keyword in ``pyazr.Session``'s signature and docstring
   (``pyazr_describe``).
4. Both values in ``tests/config_options/check.sh`` and
   ``tests/pyazr/config_options_test.py``; its line in
   ``reference/command_line.rst``.

An option that is not project state (a file of this machine, the menu's mode)
is a field of ``RuntimeOptions`` (``api/src/pybind_api.cpp``) with its row and
reason in ``NotProjectOptionRows()`` (``api/src/pybind_describe.cpp``).

A command-line flag that is not an option of the calculation:
``parseOptions`` and ``printHelp`` (``src/AZURE2.cpp``) and
``reference/command_line.rst`` (check 1).

A fitter
--------

A minimiser the program and a pyazr session can both run:

1. Its function, ``FitResult FitX(AZURECalc &, AZUREParams &, const
   FitControl &)``, in a file of its own (``src/FitX.cpp``, in
   ``AZURE2_CORE_SOURCES``): it starts from the values, limits and fixed flags
   of the ``AZUREParams``, leaves its end point there, writes its lines to
   ``control.log``, hands its steps to ``control.progress`` (an evaluation
   fitter puts a ``FitEvaluations`` in place while it runs, which also keeps
   the best point and puts the stop request in force; an iteration fitter
   sends a ``FitProgress`` per iteration and checks ``control.stop`` itself)
   and fills ``FitResult`` (status, calls, covariance where it has one; a
   stopped fit returns the best point seen with status ``stopped``).
2. Its row in ``Fitters()`` (``src/Fitters.cpp``): name, the program's switch,
   the ``Config::paramMask`` bit the switch sets, gradient / residual rows /
   bounds / covariance / MINOS, progress label, one line of doc. Its place in
   the table is its precedence in the program.
3. The switch in the program's parser (``src/AZURE2.cpp``) and its row in
   ``reference/command_line.rst``; its row in the table "The Engine's
   Fitters" of ``user_guide/fitting.rst`` (check 46 fails for a row missing
   from either, and for a minimizer bit of the parser that is no row).
4. Its case in ``tests/pyazr/fit_with_engine_test.py``: ``fit_with_engine``
   ends where the program ends, with the program's covariance.
   ``pyazr.kinds()["fitters"]`` and ``fit_with_engine(method=)`` read the
   table and need no edit.

A quantity of the progress record is a field of ``FitProgress``
(``include/FitControl.h``), set where the record is made (``FitEvaluations``
for evaluations, the fitter for its iterations) and copied into the
dictionary of ``FitWithEngineDict`` (``api/src/pybind_fit.cpp``); the
program's printer (``ProgramProgress::Print``) is unchanged unless it is to
print it.

An environment switch
---------------------

A diagnostic (an ``AZURE_*`` or ``AZURE2_*`` variable that prints, checks or
selects another route):

1. A value of ``azure2::support::Switch``
   (``support/include/azure2/support/Environment.h``) and its row, in the
   same order, in the table of ``support/src/Environment.cpp`` (name, type,
   value when unset, kind, what it does, its page). A switch that can move
   a number also names its option of the program in the row (``--name`` for
   a flag, ``--name=VALUE`` with its range for a value): the parser, ``--help``
   and ``--describe`` take the option from the row (check 43 fails for a
   row of that kind without one).
2. Its row in "Environment switches" of ``reference/command_line.rst``,
   with the option (checks 40 and 43 keep the page and the table the same).
3. Its reader: ``azure2::support::ProcessEnvironment()`` where it is used;
   a switch that can move a number is a field of
   ``Config::EnvironmentSettings``, filled in ``SettingsOf``
   (``src/Config.cpp``) and read from the run's ``Config``. ``getenv`` is
   called nowhere else (check 39).
4. A counter it prints is a member of ``RunDiagnostics``
   (``include/RunDiagnostics.h``), never a ``static`` (check 38).

A physics function of the binding
---------------------------------

1. The library computes it (``physics/``, ``experiment/``, ``stopping/``),
   with its unit test against an independent value; ``physics/`` includes no
   engine header and no GSL (check 8) and takes its constants from
   ``Constants.h`` (check 10).
2. One row in ``BuildTable()`` of ``api/src/pybind_physics.cpp``: name,
   summary, its group in the list, arguments with unit and domain, returns,
   and an adapter that only converts and calls the library. Each argument
   without a default carries an ``example``, and the argument a plot runs
   over a ``sweep``; ``tests/pyazr/physics_table_test.py`` fails without
   them, because the web UI's Physics tab draws every row from these values
   without typing. The row types (``ArgSpec``, ``OutSpec``, ``ArgKind``) are
   in ``api/src/pybind_rows.h``.
3. A ``row_test`` in ``tests/pyazr/physics_table_test.py`` against mpmath, a
   published table or the engine's own use (check 41 fails without one).
4. ``python3 tools/gen_physics_docs.py`` (the user-guide page;
   ``pyazr_physics_table`` fails until it is run). ``pyazr/physics.py`` builds
   the Python function from the row; nothing is written in Python.

A stopping function is one row of ``BuildStoppingTable()`` in
``api/src/pybind_stopping.cpp``, of the same row types, with its citation;
``tests/pyazr/ui/physics_tab_test.py`` plots every row of both tables and
runs the Python text the page offers for it. Neither kind of row needs
anything in ``pyazr/ui`` or in the page. A new kind of argument is a value
of ``ArgKind`` with its name in ``KindName`` (``api/src/pybind_rows.h``) and
one row of ``ARGUMENT_KINDS`` in ``pyazr/ui/static/physicsform.js``
(``tests/pyazr/ui/js/physics.test.mjs`` fails until it exists).

A command of ``python -m pyazr``
--------------------------------

1. One ``Command`` row of ``COMMANDS`` in ``pyazr/cli.py`` (name, ``Arg`` rows
   with a type of ``ARG_TYPES``, the ``pyazr`` function it calls, the help
   line) and its ``run`` function, which only calls that function.
2. Its case in ``JSON_CASES`` of ``tests/pyazr/cli_test.py`` (``pyazr_cli``).
3. The reference page ``user_guide/pyazr_command_line.rst``, regenerated with
   ``cli.write_reference_page()`` (check 42 fails for a row without a case, a
   case without a row, an unknown argument type or a stale page).

A run-time data file of pyazr
-----------------------------

1. The file under ``pyazr/`` (the wheel carries the package directory as it
   is), read by its module from next to itself (``__file__``), never through
   an environment variable.
2. Its row of ``DATA_FILES`` in ``pyazr/datafiles.py`` (path, what it is, who
   reads it); a file the wheel must not carry is a row of ``NOT_SHIPPED``
   instead, with the same pattern in ``wheel.exclude`` of ``pyproject.toml``.

Check 48 fails for a tracked non-Python file under ``pyazr/`` in neither list,
a row naming no file, or ``wheel.exclude`` differing from ``NOT_SHIPPED``;
``pyazr_wheel`` (slow) fails when the built wheel lacks a listed file or holds
an unlisted one. A new C++ source of the binding is one name in the
``pybind11_add_module(_azure2 ...)`` list of ``api/src/CMakeLists.txt``
(nothing for the wheel). A new optional dependency is one entry of
``[project.optional-dependencies]`` in ``pyproject.toml``.

A route and a tab of the web UI
-------------------------------

A **route** is one row of ``ROUTES`` in ``pyazr/ui/routes.py`` (a ``Route``:
method, path, the callable it serves and how a request maps onto its
arguments); there is no handler written per route.
``tests/pyazr/ui/routes_consistency_test.py`` fails for a row whose callable
does not exist, and for a public method of ``pyazr.Project`` or of the
session that is neither served nor listed in ``PROJECT_NOT_SERVED``
(``routes.py``) or ``SESSION_NOT_SERVED`` (``pyazr/ui/engine.py``). A **job
kind** (a recipe, a sampler) is a function ``(tracker, x0, settings, names)``
and one row of ``JOBS`` in ``pyazr/ui/jobs.py``
(``tests/pyazr/ui/jobs_test.py``). A fitter of the engine needs nothing in
the web UI: its job rows (a fit, and a MINOS job after it where the fitter's
row allows MINOS) are made from the engine's fitter table by
``engine_rows`` (``tests/pyazr/ui/routes_consistency_test.py``,
``engine_jobs_test.py``). A new type of a job's setting is a row of
``SETTING_TYPES`` (``jobs.py``) and of ``SETTING_WIDGETS``
(``pyazr/ui/static/tabs/fitjob.js``).

A **tab** is a module ``pyazr/ui/static/tabs/<name>.js`` that exports
``render(ctx, element, args)`` and one row of ``TABS`` in
``pyazr/ui/static/app.js``; add its view to ``VIEWS`` of
``tests/pyazr/ui/smoke_test.py``. No test fails when a tab module is left out
of ``TABS``; the smoke test (headless Chrome, ``ctest -L browser``) renders
the views it lists. No module of the page but ``layout.js`` names a record
kind, quantity or parameter kind (``tests/pyazr/ui/static_test.py``).

Other tables
------------

The repository's ``CLAUDE.md`` has the recipes for a model-only datum
(``ModelQuantities()``, check 27), a derived parametrisation of a level
(``gamma_decay``), a new optical potential term or global optical set
(``physics/src/transfer``), a material or stopping source (``stopping/``,
``TargetStopping``), a cache, and a regression check (``tests/<name>/check.sh``
with ``tests/lib/edit_guard.sh``).
