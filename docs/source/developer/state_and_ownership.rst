State and Ownership
===================

Who owns each object of a run, which objects are immutable, and why a result
does not depend on the number of threads. The rule behind it: the only
process-wide objects are immutable tables built before the first
evaluation; everything an evaluation writes belongs to the caller and is
passed in. :doc:`classes` marks the immutable and the per-evaluation types in
its figures.

The objects of a run
--------------------

The project as read
   ``Config::ReadConfigFile`` reads the project once and keeps it as
   ``Config::project``, a ``std::shared_ptr<const Format3Project>``
   (``experiment/include/azure2/experiment/ProjectReader.h``): the validated
   records, the ``<levels>`` rows and the specs of ``BuiltProject``. Nothing
   writes it afterwards; a copy of the ``Config`` shares it. ``pyazr.Project``
   holds its own ``Document`` (``_azure2.Document``) and edits that one; a
   session reads the file again when it is built.

``Config`` -- the run
   One per program run or per ``pyazr`` session: the options, the paths, the
   environment's settings (``Config::environment``, read once from the
   process's environment through ``azure2::support::ProcessEnvironment``),
   and three owned places for what a run keeps besides its results:
   ``Config::scratch`` (``RunScratch``, ``include/RunScratch.h``: a slot per
   thread and tables the threads share, such as the external-capture term
   store), ``Config::diagnostics`` (``RunDiagnostics``,
   ``include/RunDiagnostics.h``: counters and the profile) and the once-only
   notes (``Config::FirstTime``). Two sessions in one process share none of
   them.

The set-up objects, immutable
   Built by the set-up, then only read, and shared by every copy of the
   working copy, every evaluation and every thread:

   - ``DataSet`` (``include/DataSet.h``): the segments, target effects,
     covariance groups, meshes, THM experiments and the parameter layout.
     The set-up writes it through ``WorkingCopy::Draft``;
     ``WorkingCopy::Freeze`` ends that, and every working copy then holds it
     as ``std::shared_ptr<const DataSet>``. A setting changed afterwards (a
     prior covariance given by the caller) makes a new ``DataSet``
     (``WorkingCopy::ReviseDataSet``). One part is written after the set-up
     on purpose and named in ``DataSet.h``: a covariance group's frozen model
     of its relative terms (``DataCovariance::Group``), frozen once under the
     group's lock at the first evaluation of a minimization and forgotten
     between two minimizations (``WorkingCopy::RefreezeCovariance``).
   - The tables of the stage: the ``NodeTable`` and the slots of every plan
     task (``NodeState``, ``include/NodeTable.h``), the shared angular rows
     (``RowSet``), the ``ChannelTable``; all held as
     ``std::shared_ptr<const ...>`` or ``SharedTable``.
   - ``CompoundSystem``, the set-up of the compound nuclei: the pairs, groups, channels
     and levels, the two ``SystemOfCompound`` of the structure, the
     ``PathwayTable``, the model data, the ``gamma_decay`` set and the
     parameter layout. One function makes it, ``SetUpCompoundSystem``
     (``include/CompoundSystemSetUp.h``): the structure from the ``<levels>``
     rows, then the caller's data set-up, which reads the structure through a
     ``const CompoundSystem &`` and hands back the entrance and capture-A_y
     requests (``CompoundSetUpRequests``), then the last part (B_c, the input
     transformation, the pathway table, the external-capture set-up). The
     set-up calls are private; the function returns the finished object as a
     ``std::unique_ptr<const CompoundSystem>`` with the owner's levels. Its one
     owner, the program's run (``ProgramRun::compound``), the sampler's run or
     the session (``AZUREAPI::compound_``), outlives every reader, which holds
     a ``const CompoundSystem *`` or ``const CompoundSystem &``. What the last
     part writes (B_c, the final-state energies, the pathway table, the
     external-capture blocks, the capture masks) refuses to be read before it
     (``std::logic_error``). It holds no value of an evaluation and cannot be
     copied.

The per-evaluation objects
   - ``WorkingCopy`` (``include/WorkingCopy.h``): what an evaluation changes
     -- the state of every segment and term (``SegmentState``: norms,
     shifts, frames), the present values of the apparatus parameters
     (``TargetEffectValues``), of the isotopic abundances
     (``WorkingCopy::AbundanceValues``) and of the THM coherent backgrounds
     (``WorkingCopy::ThmCoherentValues``, handed to the HOES model), the points (``Unit``) and slots (``Site``) with
     their energies after a shift, the plan's positions. ``WorkingCopy::Clone``
     is a plain copy that shares the ``DataSet``.
   - ``LevelState`` (``include/LevelState.h``): the levels at the present
     parameter vector (energies, amplitudes, shift functions, Park's
     :math:`J`, the external-capture parts, the ``gamma_decay`` fills), and
     ``LevelInputs`` (``include/LevelInputs.h``) where an owner transforms
     from the observed basis. Each owner holds its own and hands it to the
     readers (the compound first, then the state, then the inputs): the
     program's run (``ProgramRun::state``, ``inputs``, made by
     ``SetUpCompoundSystem``), the session (``AZUREAPI::state()``,
     ``inputs()``: a method of ``AZUREAPI`` that writes the session's levels, data, ``Results`` or its lazily
     read fit settings is not ``const``, so a ``const AZUREAPI &`` only reads), a fit's evaluations (copies of the
     caller's state in its
     ``FitterStates``, ``include/AZURECalc.h``), the sampler (its own
     ``WalkerLevels`` and a copy per walker and per snapshot), the reaction
     rate (``ReactionRate::State``, a copy per thread), a limit's conversion
     (``ParameterLimitsManager``), the sub-point grids
     (``WorkingCopy::GridLevels``). A pool of idle copies is a
     ``StatePool`` (``include/StatePool.h``) owned by the caller, as the
     ``Results`` are.
   - ``azure2::Results`` (``include/Results.h``): the evaluation's buffer --
     one record per point and slot at an offset fixed by the layout
     (``WorkingCopy::LayoutResults``), the model values, the contracted
     angular values, the objective's totals, and the scratch of the energy
     tasks (``Results::scratch``, one ``TaskScratch`` per thread: the
     level-matrix solvers and their buffers). It belongs to whoever starts
     the evaluation: the program's calculation, each pooled copy of a fit or
     of the sampler (``WorkingData``, a working copy with its buffer), the
     Python session. A point holds the offset of its record, never a pointer
     to a buffer; every function that reads or writes a result takes the
     buffer as an argument.

The rule and its check
----------------------

``tests/consistency/repo_consistency_test.py`` holds the rule:

- check 38: no ``thread_local``, no non-``const`` ``static``, no
  namespace-scope variable and no ``mutable`` member outside the allow-list
  ``tests/consistency/allowed_state.tsv`` (scanner
  ``tests/consistency/state_scan.py``). A row whose reason starts with
  ``table:`` is permanent: GSL's process-wide error handler, the Coulomb
  tables loaded from files, the cache of input files, the locks of
  process-wide output streams, the version generator of ``LevelState``. A row whose
  reason is ``pending`` names state still to be moved to its owner: the pools
  and the gradient memo of ``AZURECalc``, the sampler's lazily filled layout
  in ``AZURECalcMCMC``. A row is never added for new code.
- check 39: the process's environment is read in one place,
  ``support/src/Environment.cpp``; no other file calls ``getenv``.
- check 25: a point or a data set holds where its record lies, never a
  buffer.

New state that a calculation keeps is one of three owned places: a member of
``azure2::TaskScratch`` (read by energy tasks), of ``azure2::ThreadScratch``
(``Config::scratch``, read by set-up or hand-over code), or a table built once
in the set-up and held ``const``. A warning printed once is
``Config::FirstTime``. The recipe is "Add a cache" in the repository's
``CLAUDE.md``.

Process-wide tables
-------------------

Built once, reached through a function that returns a ``const`` reference,
never written: the constants (``physics/include/azure2/physics/Constants.h``),
the schema and its quantities, units and kernels (``Schema``,
``Quantities()``, ``Units.cpp``, ``Kernels()``), the parameter kinds
(``ParamKindTable()``), the options of the calculation (``io::ConfigKeys()``),
the objective's term list (``Objective::Engine()``), the environment switches
(``azure2::support::ProcessEnvironment``), the physics table of the binding
(``BuildTable()`` in ``api/src/pybind_physics.cpp``).

Threads
-------

The engine runs its loops with OpenMP. What runs in parallel and what each
thread may write:

- **The point pass** (``WorkingCopy::CalculatePointsInParallel``,
  ``src/EnergyPass.cpp``). The energy tasks run on a dynamic schedule; each
  writes the records of its own units and nothing else, and reads the
  set-up's tables and the evaluation's ``LevelState``. The folds run
  after every energy task has finished; each writes the records of its own
  points. A task that throws marks its points failed (their value is NaN);
  the failure is collected after the pass, not in the middle of it. Under
  ``AZURE_PLAN_CHECK=1`` the pass checks at run time that each step writes
  only its own slots and that the shared solve of a task gives every unit the
  bits its own solve would give (``CalculationPlan::RunChecked``,
  ``CheckEnergyTasks``).
- **The gradient** (``src/AZUREGrad.cpp``, ``src/GradEnergyTasks.cpp``): the
  forward pass as above, then one adjoint pass per energy task.
- **A fit's evaluations.** During a fit every evaluation of ``AZURECalc``
  takes its own ``LevelState`` and ``WorkingData`` from the pools
  (``AZURECalc::GetPooledState``, ``GetPooledData``) and returns them after
  use, so evaluations that run at the same time never share one. An
  evaluation already inside a parallel region (Minuit2's numerical
  gradient) runs the same point pass on its own thread alone
  (``CalculatePointsInParallel`` opens no second team): one solve per
  compound energy, written into the evaluation's own ``Results``. The
  sampler does the same with its own pools.
- **The set-up's loops** (energy-dependent values, external capture, node
  and channel tables) write one entry per index into storage sized before
  the loop.

Why a result does not depend on the thread count:

1. every value is written by exactly one task, from inputs that do not
   depend on which thread runs it or on what ran before (no memo keyed by
   the history of a run: the Coulomb functions are computed at every call,
   the scratch holds no value an evaluation reads back);
2. every sum that runs over tasks is formed in a fixed order: the
   objective after the pass, segment by segment; the gradient's sums in
   chunks fixed independently of the thread count, each in task order, and a
   slot's cotangent summed over its readers in plan order;
3. the set-up tables are immutable while a fit runs.

To compare two runs value by value, ``AZURE_TRACE_FCN`` prints each objective
value with its parameter vector in hexadecimal floating point;
``AZURE_TRACE_FCN=fresh`` also evaluates on a fresh copy of the master objects,
which must give the same bits as the pooled copy
(:doc:`../reference/command_line`, "Environment switches").

.. note:: **Planned.**
   The remaining pools and ``mutable`` members of ``AZURECalc``,
   ``AZURECalcMCMC`` go to their owners (the fitter's adapter, the run), as the
   pools of level states did (``StatePool``); the allow-list then keeps only
   its ``table:`` rows.
