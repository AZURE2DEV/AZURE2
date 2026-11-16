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
2. Its formula at the energy, as a function of the T-matrix elements and
   plain numbers in ``physics/include/azure2/physics/Observables.h``, with
   its closed-form check in ``physics/tests/observables_test.cpp`` (no
   project file needed).
3. One branch of ``GenMatrixFunc::CalculateBareCrossSection`` that walks the
   pathways, calls the function and writes the record; its adjoint in
   ``AMatrixFunc::PointTBar`` (the derivative is still written there, not in
   ``physics/``).
4. What a reader needs per slot (angular row, Coulomb amplitude) in its row
   set (``include/RowSet.h``, built by ``WorkingCopy::BuildRowSets``).
5. Tests: the analytic gradient against finite differences
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

1. ``inference/include/azure2/inference/ParameterKind.h``: the value of ``ParamKind`` and its row of
   ``ParamKindTable()`` (key, unit, owner fields, prior and vary flags, doc);
   a ``ParameterRecord`` factory for new owner fields.
2. Create it with ``AZUREParams::Add`` in ``CompoundSystem::FillMnParams`` or
   ``WorkingCopy::FillMnParams``; every later question asks the record, never
   the name (check 23).
3. Its place in ``BuildParamIndexMap`` (``src/AZUREGrad.cpp``): one
   ``Append...`` function of ``ParamIndexMap``
   (``inference/include/azure2/inference/ParameterIndexMap.h``) with its
   lookup, called at the kind's place in the vector, and its case in
   ``inference/tests/parameter_index_map_test.cpp`` (a level's
   number takes its place from ``LevelParameterLayout``, above); the fields
   that say what it belongs to are rows of ``ParameterFields()``
   (``api/src/pybind_records.cpp``), which ``pyazr.Parameter`` reads
   (``tests/pyazr/parameter_records_test.py`` fails for a kind no owner field
   describes).
4. Docs: the kinds table of ``reference/output_files.rst``; a project with
   the kind in ``tests/pyazr/parameter_kinds_test.py``.

A level parameter
-----------------

A number per level (or per level and channel) that the fit varies beside the
energy and the amplitudes:

1. ``LevelParameterLayout.h`` / ``LevelParameterLayout.cpp`` of ``inference/``
   (``inference/include/azure2/inference/``, ``inference/src/``; the test
   ``inference/tests/level_parameter_layout_test.cpp``): a
   ``Slot::Kind`` and its slots in the constructor, which fixes where the
   number sits in the level block, and its name in ``Name``. Every walk of
   the level block reads the layout (``CompoundSystem::FillMnParams``, the fills from a
   fit vector and from an observed vector, ``GetTransformParams``,
   ``ConvertAmplitudeBasis``, ``SetFitLevelParameter``, the level part of
   ``BuildParamIndexMap``, the labels of ``src/ParameterLabel.cpp``, the
   ``gamma_decay`` slots, the model data's derivative columns), so none
   counts positions of its own.
2. The value at the vector: a member of ``LevelState``
   (``include/LevelState.h``) and its branch in the fills that switch on the
   slot's kind (``FillCompoundFromParams``, ``SetFitLevelParameter``), which
   write the ``LevelState`` they are given. What a fill computes from it
   beside the level (as a ``gamma_decay`` record's ``GammaDecayFill``) is a
   member of ``LevelState`` too. Its readers take the state as an argument
   (the compound first, then the state, then the inputs), as every reader of
   an evaluation's values does.
3. The input value: a member of ``LevelInputs`` (``include/LevelInputs.h``),
   made by ``CompoundSystem::Fill`` (in ``SetUpCompoundSystem``) from the ``<levels>`` rows, written by
   ``TransformIn`` and ``FillCompoundFromParamsPhysical`` (the observed
   vector), read by ``CompoundSystem::FillMnParams`` as the start value; its observed
   value in ``GetTransformParams``. What an evaluation reads of the input
   (the zero-width flags) is
   written into the ``LevelState`` where the input is read (``Fill``,
   ``TransformIn``, the fill from an observed vector), so that the
   evaluation reads the state alone. ``ALevel`` holds only what the
   ``<levels>`` rows say (``FileE``, ``FileGamma``, ``FileZeroWidth``).
4. A parameter kind (below) for its record, its pull-back in
   ``GradAccum`` (``src/AZUREGrad.cpp``) if the gradient is analytic, and its
   readers in the matrix functions.
5. Tests: ``tests/upgrade/level_parameter_layout_test.cpp`` (the layout
   against the names ``FillMnParams`` gives the block),
   ``tests/upgrade/level_inputs_test.cpp`` (the input values and what a fill
   leaves in the state) and the objective's gradient against finite
   differences.

An owner of the levels
----------------------

Something that evaluates the model on levels of its own (a second sampler, a
scan, a new front end) holds the compound as ``const CompoundSystem *`` (or ``const
CompoundSystem &``) of the run or session that set it up and outlives it, and its own
``LevelState`` (with a ``LevelInputs`` if it transforms from the observed
basis), copied from its caller's: ``LevelState state(callerState)``. It
passes the compound, the state and the inputs to the fills and the readers.
If evaluations of it run at the same time from ``const`` methods, the idle
copies are a ``StatePool`` (``include/StatePool.h``) of its caller, lent as
the ``Results`` are (``FitterStates`` of ``AZURECalc``,
``SamplerStates`` of ``AZURECalcMCMC``); so is anything else they must
keep between calls (the fit's last gradient is in ``FitterStates``): a
``mutable`` member is refused by the consistency test. No file of the engine
changes: the owner's own file, and its test.

A new front end that sets up a project itself (as the program's run, the
sampler and the session do) makes the compound and its data with the one
call, ``SetUpCompoundSystem`` (``include/CompoundSystemSetUp.h``), and owns
what it returns:

.. code-block:: cpp

   WorkingCopy data;
   CompoundSetUp setUp = SetUpCompoundSystem(config, [&](const CompoundSystem &structure,
                                                         const LevelInputs &inputs)
                                                         -> std::optional<CompoundSetUpRequests> {
     if (data.Fill(config, &structure, inputs) != 0) return std::nullopt;  // refused: it said why
     CompoundSetUpRequests requests;
     RequestDataEntrances(structure, config, requests);
     RequestCaptureAyTerms(structure, data, requests);
     return requests;
   });
   if (setUp.refusal != CompoundSetUp::Refusal::kNone) return;  // say why in the front end's words
   // setUp.compound (std::unique_ptr<const CompoundSystem>), setUp.state, setUp.inputs: the front end's

``tests/upgrade/SetUpForTests.h`` has this form for the tests. A front end of the program (one more command
beside the run and the sampler) calls ``SetUpProgramCompound`` (``include/ProgramRun.h``) instead: the
program's data set-up, its check files and its words on a refusal (``PrintSetUpRefusal``) come with it.

A set-up fact of the compound system
------------------------------------

Something the compound system learns once from the project and that no evaluation changes (a property of a
pair, a group or a level, a table derived from them) is a member of ``CompoundSystem``
(``include/CompoundSystem.h``) with a ``const`` accessor, written by its set-up only, which
``SetUpCompoundSystem`` (``src/CompoundSystemSetUp.cpp``) runs: ``src/CompoundSystemBuild.cpp`` (``Fill`` for
what the data's set-up may read, ``Initialize`` for what needs the data's requests or the transformed
inputs), or ``src/ExternalCaptureFinalStates.cpp`` for the external capture. A fact ``Initialize`` writes has
an accessor that refuses before it (``CompoundSystem::RefuseUnfinished``), since the data's set-up sees the
structure before the last part. A fact the last part needs from the data is a field of
``CompoundSetUpRequests``. Its readers take ``const CompoundSystem &``. If the check file of the compound
nucleus shows it, ``WriteCompoundNucleusCheck`` (``src/LevelReport.cpp``) prints it. A value that moves with
the parameter vector is not a set-up fact: it belongs in ``LevelState`` (above).

A term of the objective
-----------------------

1. One ``ObjectiveTerm`` subclass in ``src/Objective.cpp``, appended to the
   term list in the ``Objective`` constructor (``Objective::Terms``).
   ``ObjectiveTerm`` is ``azure2::inference::Term<ObjectiveContext>``
   (``inference/include/azure2/inference/TermList.h``): ``Add`` is the term,
   ``ModelOnly``, ``AddForSampler`` and ``ExcludesForSampler`` have defaults
   (not of the model alone; the sampler takes the whole term; no point is
   excluded). A Gaussian prior on one parameter is ``GaussianPriorChi2`` of
   ``GaussianPrior.h`` and a line of ``ForEachDiagonalPrior``
   (``include/Objective.h``), which gives it its gradient, Minuit's
   :math:`G_2` and its least-squares row. Its
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

A field a user sets rarely is wrapped in ``Rare(...)`` in the same row: the
schema writes ``"rare": true`` and the web page's first view folds it under
"More fields" while the record does not set it (``pyazr/ui/static/forms.js``,
``folded``; ``tests/pyazr/ui/js/firstview.test.mjs``). Its ``doc`` is the
meaning the first view writes under the field's title: a field needs a doc
that says what it is, not where it is described.

A **record kind**: a function ``<Name>Kind()`` in ``experiment/src/Schema.cpp``,
its entry in ``Schema::Kinds()`` and its section in ``Schema::Sections()``;
its build in ``Schema::Build``; a fixture, the schema and the page as for a
key. A **diagnostic code** is a constant in
``experiment/include/azure2/experiment/Diagnostics.h``, its row in
``DiagnosticCodes()`` and a case that triggers it alone in
``experiment/tests/diagnostics_test.cpp`` (check 24). A **unit** is a row of
``experiment/src/Units.cpp`` (``experiment_parsers``,
``experiment_quantities``).

A **section** (a block ``<name> ... </name>`` of the project file): a section
of records is the record kind's row of ``Schema::Sections()``; a section with a
grammar of its own (as ``<levels>``, ``<thm>``, ``<parameterSettings>``) is a
row of ``Schema::GrammarSections()`` (``experiment/src/Schema.cpp``: its name,
the code that reads it, the sections it may hold) and its reader, which opens it
with ``Config::SectionStream("<name>")``. The reader of the project
(``ReadFormat3Project``, ``experiment/src/ProjectReader.cpp``) refuses every
section in neither table, and a section inside one that does not list it
(``E_UNKNOWN_SECTION``): a block nothing reads would leave the run as if it were
absent. A ``<config>`` key that is not a row of ``io::ConfigKeys()``, a
``<config>`` line that is not ``key=value`` and a line of ``<mcmc>`` outside its
``<parameters>`` are refused the same way (``E_UNREAD_LINE``). The converter
refuses an old section it does not convert (``kConvertedSections`` in
``convert/src/Convert.cpp``); the sections of a newer AZURE2 that this build
does not compute (``<channelRadii>``, ``<isospin>``) are named with what they
state in ``Schema::NewerSections()``, whose words both refusals use. A section
that a newer program adds is therefore refused until its row and its reader
are written.

``pyazr.Project`` edits every record through the binding of the engine's
grammar and needs no Python change; an unedited project is written back byte
for byte (``tests/pyazr/project_test.py``).

A session option (an option of the calculation)
-----------------------------------------------

A switch a project's result depends on is a ``<config>`` key, a program flag
and a ``pyazr`` keyword under one name.

1. One ``Option`` row of ``io::ConfigKeys()`` (``io/src/ConfigSection.cpp``),
   wrapped in ``Told(...)`` with the words a form shows (a title, the old
   program's where it had the setting; the group; what the engine does
   without it; whether it is set rarely): ``pyazr.kinds()["config"]`` and the
   web page's run settings read them, so a new key needs nothing in the page
   (``tests/pyazr/ui/api_test.py``);
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

An option that is not project state (a file of this machine, the program's command)
is a field of ``RuntimeOptions`` (``api/src/pybind_api.cpp``) with its row and
reason in ``NotProjectOptionRows()`` (``api/src/pybind_describe.cpp``).

A command-line flag that is not an option of the calculation nor of a command:
``parseOptions`` and ``printHelp`` (``src/AZURE2.cpp``) and
``reference/command_line.rst`` (check 1).

A command of the program, or an option of a command
---------------------------------------------------

The program's commands (its run modes) and every answer a run needs are one
table, ``include/RunRequest.h`` and ``src/RunRequest.cpp``; the command line
(``ReadRunCommand``), ``--help`` (``RunCommandsHelp``) and ``--describe
--json`` (``RunCommandsJson``) read it, and ``PrepareRun``
(``src/RunPreparation.cpp``) turns the answers into the run's ``Config``.
The program reads nothing from standard input.

A new option of a command (a question a run needs answered):

1. Its member in ``RunRequest`` (with its unit in the comment).
2. Its row in ``MakeQuestions()`` (``src/RunRequest.cpp``), at the place
   ``--help`` lists it: option, kind (flag, word, number, whole number,
   file), the commands that take it, the condition (``AskedWhen``), required
   (``Required``, or ``RequiredOr`` for one of two ways) or its default, its
   range, help, and the member it fills.  The range and default are printed
   from the row, not written in the help text.  A rule between two answers
   is a row of ``RunRules()``.
3. Where the answer goes in ``Config``: ``ApplyAnswers``
   (``src/RunPreparation.cpp``).
4. Its row in the options table of ``reference/command_line.rst``
   (``tests/cli_help``, ``tests/cli_commands``) and, when it changes what a
   command writes, a case of ``tests/cli_commands/check.sh``.

A new command: a ``RunCommand`` value and its row of ``RunCommands()``
(name, help), its bits in ``SelectCalculation``
(``src/RunPreparation.cpp``), the ``Commands(...)`` of the rows it takes, what
``main`` (``src/AZURE2.cpp``) runs for it if it is not ``AZUREMain``, and its
row in the commands table of ``reference/command_line.rst``.

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
   and fills ``FitResult`` (its own status word, which ``FitOutcome`` maps
   to the shared words, ``converged``, calls, covariance where it has one,
   ``notices`` for what the caller should know; ``x`` and ``objective`` are
   the lowest objective it reached: a stopped fit, or one ended by its call
   limit, returns the best point seen).
2. Its row in ``Fitters()`` (``src/Fitters.cpp``): name, the program's switch,
   the ``Config::paramMask`` bit the switch sets, gradient / residual rows /
   bounds / covariance / MINOS, whether it honours a call limit
   (``callLimit``: ``RunFitter`` adds the notice for a fitter that does not),
   the row it falls back to without an analytic Jacobian (``fallback``), the
   progress label, one line of doc. Its place in the table is its precedence
   in the program. The program's parser and usage text take the switch and
   the doc line from the row (``ApplyFitterSwitch``, ``FitterSwitchesHelp``);
   a fitter behind a build option also gets a row in ``AbsentFitters()`` for
   the builds without it, so that its switch is refused with how to build it.
   ``fitStatus.out`` (one writer, ``WriteFitStatus``) needs no edit.
3. The switch's row in ``reference/command_line.rst`` and the fitter's row
   in the table "The Engine's Fitters" of ``user_guide/fitting.rst`` (check
   46 fails for a row missing from either, for a switch named in
   ``src/AZURE2.cpp`` outside the table, and for a minimizer bit of the
   parser that is no row).
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

A section (MF/MT) of ``pyazr.endf``
-----------------------------------

The format work of ENDF-6 is ``pyazr/endf/`` (user guide :doc:`../user_guide/endf`); every rule it codes cites
the printed page of ENDF-102 (2023) in a comment. A section that can be read and evaluated:

1. Its records: a ``read_*`` function and a record class in ``pyazr/endf/distributions.py`` (or a new file of the
   package for a new file MF), read with the ``Cursor`` of ``records.py`` and ending with ``_require_done``; a
   method of ``Material`` (``material.py``) that calls it.
2. Its words in ``Material.summary()``: a branch of ``_describe_section`` (``summary.py``); a new reaction or file
   name is a row of the tables of ``names.py``.
3. Its record counts in the own check: a branch of ``_check_records`` (``format_check.py``).
4. If it gives a quantity of a data set: a branch of ``curves`` (``curves.py``); physics (a kinematic or Coulomb
   formula) is a row of ``pyazr.physics``, never a formula in the package.
5. A hand-written file of the section (``tests/pyazr/endf_testfiles.py`` writes it) and checks against hand values
   from the manual's equations in ``tests/pyazr/endf_read_test.py``.

A section (MF/MT) of the ENDF writer
------------------------------------

``pyazr.endf.write`` (``pyazr/endf/writer.py``) writes what the engine computes; the module writes records only.
MF=6 (LAW=2 of the exits, LAW=5 of elastic scattering, the primary photons of MT=102) was added this way.

1. The reaction: a role of ``engine_values.Reaction`` made by ``_reactions`` (``writer.py``) from the session's
   pairs (``"elastic"``, ``"exit"``, ``"capture"``, ``"final_state"``).
2. The numbers: an engine quantity of an extrapolation in ``engine_values.Sweeps`` (``sigma``, ``legendre`` with
   the expansion that is the format's own normalisation, ``dsigma_domega`` at angles), on the one energy grid of
   every section (``Sweeps.energy_grid``: ``grid.refine`` from ``grid.start_grid``, the rule applied to the
   quantities the file tabulates; a second axis, as the cosines of LAW=5, refined by the same ``grid.refine``). If
   the engine has no quantity for it, a door of the binding (a new file of ``api/src/``, registered in
   ``BindPhysics``, or a row of the point outputs of ``api/src/pybind_evaluation.cpp`` when the numbers are a
   record the evaluation writes); never a physics formula in Python. LAW=2 needed no door: the ``legendre``
   quantity with the expansion ``rsl`` is its a_l. LAW=5 LTP=1 has one: the point output ``nuclear_amplitudes``
   (``api/src/pybind_elastic_amplitudes.cpp``) reads the elastic record of a data set's point and gives the
   manual's b_l and a_l (``Sweeps.nuclear_amplitudes``; the energies chosen by ``pyazr/endf/amplitudes.py``). A section whose reader does not interpolate linearly between energies (LAW=5
   LTP=12: sigma_NI of MF=3 times P_NI of MF=6, each interpolated) gets energies of its own, chosen by
   ``grid.thin_to`` with the reader's own reconstruction instead of the chord (``pyazr/endf/elastic.py``).
3. The records: a function of ``sections.py`` that lays the section out as its page of ENDF-102 says (cite it),
   with a ``RecordWriter`` (``records.py``), and the line of ``write`` that puts it into ``written[(MF, MT)]``.
   The directory of MF=1 and the sequence numbers follow from the sections.
4. What it cannot carry: a row of ``NOT_WRITTEN_REASONS`` for a section a reaction needs and the writer leaves out;
   a line of ``Losses`` (``model.py``) for what the project holds and the file does not.
5. Reading it back: the reader's ``read_*`` of the section and a branch of ``curves`` (the recipe above; MF=6 needed
   capture to one final state in ``curves``); checks in ``tests/pyazr/endf_write_test.py`` that the file's values
   equal the engine at the data's points and between the grid's, and the own format check clean.

A checking code of ``pyazr.endf``
---------------------------------

The codes are found, never installed (``pyazr/endf/checkers.py``):

1. A row of ``PROGRAMS`` (name as the program file is called, family, kind: check, rewrite or derive, what it does,
   what it does to the numbers of a file: the page says it before a run). The page needs no line: its ENDF tab lists
   the codes from ``GET /api/endf/checkers`` and runs them through the jobs ``endf_check`` and ``endf_standardize``.
2. Its deck: for a checking code the one of ``_run_check``, for a rewriting code the branch of ``standardize`` (a
   PREPRO code: ``_prepro_deck``); take the deck from the program's source, cite the lines.
3. Its report: a branch of ``parse_report`` that gives findings (severity, MAT, MF, MT, record, text) and the
   totals the report states; known messages with their meaning are rows of ``MEANINGS``.
4. The stand-in ``tests/pyazr/endf_standin/standin.py`` writes its report in the program's layout (lines copied from
   a real report, cited); cases in ``tests/pyazr/endf_check_test.py``, and a case with the real program in its last
   part (run when ``AZURE2_ENDF_UTILITIES`` / ``AZURE2_PREPRO`` are set).
5. The command's choices in ``pyazr/cli.py`` (``endf check --program``, ``endf standardize --program``).

An example script of pyazr
--------------------------

An example tells one task from start to end on a project of the repository, for an expert (or an agent working for
one) to copy:

1. The file ``pyazr/examples/<name>.py``: a docstring that says the task, the calls it shows and the command that runs
   it from the repository's root on a project of ``tests/``; the project as its first argument, required; files
   written only where ``--out`` (or a printed temporary directory) says; every printed number with its unit; pyazr's
   own calls for what pyazr computes.
2. Its row in ``ROWS`` of ``tests/pyazr/examples_test.py``: the arguments, one number or text of its output that is
   pinned, ``slow=True`` above about 4 s at one thread, ``needs=`` the optional packages, ``skip=`` the reason a run
   cannot be made here (the network).
3. Its row on the page of examples (:doc:`../user_guide/pyazr_examples`), under its task: the script, the task in one
   sentence, the command.

Check 58 fails for a script without both rows, for a row whose script is not there, and for a page of the
repository's skills (``.claude/skills``) that names an example that is not there: a script folded into another or
removed leaves the skills too.

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
``engine_jobs_test.py``). A job that starts from another job's result (a
band from a fit's covariance or a sampler's chain, a chain continued) names
that job's kind in ``follows`` and the key of its result in ``takes``; the
server hands it that job's record (``tracker.source``), and the page offers
it beside such a result (``tests/pyazr/ui/sampler_bands_test.py``). A new
type of a job's setting is a row of
``SETTING_TYPES`` (``job_settings.py``) and of ``SETTING_WIDGETS``
(``pyazr/ui/static/tabs/fitjob.js``).

A **tab** is a module ``pyazr/ui/static/tabs/<name>.js`` that exports
``render(ctx, element, args)`` (and ``shown(ctx)`` when it is about a part of
a project that may be absent, as the THM tab) and one row of ``TABS`` in
``pyazr/ui/static/app.js``; add its view to ``VIEWS`` of
``tests/pyazr/ui/smoke_test.py``. The page has two depths: the expert's shows ``TABS``
as they are; the first view shows the rows of ``FIRST_VIEW`` in
``pyazr/ui/static/layout.js`` (the old program's path and words), each a view
of one tab (``tab``, the route it opens, a ``part`` of the Nucleus tab or the
Data tab's ``sections``). A tab draws less in the first view where
``ctx.expert`` is false (the rest under a ``details.more``); a new tab is in
the expert's depth only until a row of ``FIRST_VIEW`` names it. A function of the session drawn by a tab is
a row of the session's functions (``Session.describe()["functions"]``,
``pyazr/session_functions.py``) and its entry in ``PAGE`` of
``pyazr/ui/physics.py`` or ``pyazr/ui/thm.py`` naming that tab in ``tab`` (a
row without an entry is drawn by the Physics tab); the server walks the
session's answer by the row's ``path`` and ``lines`` and the tab draws it with
``functionview.js``. A key of ``<thm>`` is one row of ``src/ThmKeys.cpp``,
served as ``pyazr.kinds()["thm"]``: nothing in the page. No test fails when a tab module is left out
of ``TABS``; the smoke test (headless Chrome, ``ctest -L browser``) renders
the views it lists. No module of the page but ``layout.js`` names a record
kind, quantity or parameter kind (``tests/pyazr/ui/static_test.py``).

A **column of a list** of the page comes from a table, never from the tab:
the lists of the Data tab (``pyazr/ui/static/lists.js``) show, per record
kind, the fields its records set with the schema's ``title`` (a new field is a
column with no line of the page; a field marked rare stays in the form), or,
for a kind the engine describes, the entries of
``pyazr.project.DESCRIPTION_FIELDS`` marked ``list`` (``kinds()["description"]``:
one row there for a new entry of ``describe``); the list of pairs shows
``ProjectPair.FIELDS`` marked ``list`` (``pyazr/nuclear_part.py``, a new
field of a pair is one row with its title and unit; the pair form takes
``PARTNER_FIELDS`` and ``ARGUMENTS`` the same way). Which of those columns the
first view keeps is one table, ``FIRST_VIEW_COLUMNS`` (``pyazr/ui/first_view.py``,
served as ``first_view`` of ``GET /api/describe``). A **control** of a form
is the row of its value type in ``WIDGETS`` (``pyazr/ui/static/forms.js``):
a ``path`` offers the project folder's files (``GET /api/files``) and is
previewed under the field (``tabs/records.js`` ``filePreview``, ``GET
/api/project/records/{name}/preview``), a ``reference`` offers the records and
"new …", a ``number_list`` a box per number and one that ``grows`` (the form
adds the next box as it is typed in). A **kind of job** the Fit tab's
first view offers, and whose button follows a finished job there, is a kind of
``FIRST_VIEW_KINDS`` (``pyazr/ui/jobs.py``); a tab lists the jobs of the rows
whose ``tab`` it is.
A **view of the first view** is a row of ``FIRST_VIEW`` (``layout.js``); a
row with ``when`` is shown under the condition of that name (``CONDITIONS``
in ``app.js``).

A part of a builder (a data covariance, a prior)
------------------------------------------------

The two builders make a file the engine reads from named parts; the page's forms come from their tables.

- **A data covariance** (``pyazr/covariance.py``): one row of ``PART_KINDS`` (its fields: name, type, default,
  unit, doc) and a branch of ``_terms`` that returns the part's terms (``diagonal``, ``vector``, ``relvector``,
  ``block``); a case in ``tests/pyazr/covariance_builder_test.py`` against a matrix worked by hand.
- **A prior on the parameters** (``pyazr/priors.py``): one row of ``PART_KINDS`` and one function of ``_BLOCKS``
  that returns the part's ``_Block`` (its parameter names, centres, covariance and notes); a case in
  ``tests/pyazr/prior_builder_test.py``.  A new type of field also needs its control in ``fieldInput`` of
  ``pyazr/ui/static/tabs/prior.js`` (the types: ``path:<suffixes>``, ``names``, ``number``, ``integer``, ``pairs``,
  ``matrix``).

Nothing else: the checks and the notes are made on the whole matrix, and the screens (``tabs/covariance.js``,
``tabs/prior.js``) read the tables through their routes (``pyazr/ui/utilities.py``, ``pyazr/ui/prior.py``).

Other tables
------------

The repository's ``CLAUDE.md`` has the recipes for a model-only datum
(``ModelQuantities()``, check 27), a derived parametrisation of a level
(``gamma_decay``), a new optical potential term or global optical set
(``physics/src/transfer``), a material or stopping source (``stopping/``,
``TargetStopping``), a cache, and a regression check (``tests/<name>/check.sh``
with ``tests/lib/edit_guard.sh``).
