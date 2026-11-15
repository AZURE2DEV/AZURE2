The Binding and pyazr
=====================

How Python reaches the engine. The user's side is in the user guide
(``pyazr``); this page is the structure behind it.

``AZUREAPI``
------------

``AZUREAPI`` (``api/include/AZUREAPI.h``, ``api/src/AZUREAPI.cpp``) is the
engine as a session: it owns a ``CompoundSystem`` and a ``WorkingCopy`` (made in
``AZUREAPI::ReplaceModel``), reads the run's ``Config``, and keeps the
session's ``Results`` buffer. ``AZUREAPI::Initialize`` makes the same set-up
calls as the program's ``AZUREMain`` (:doc:`architecture`, step 2), so a
session and a run of the program evaluate the same objects. Its entry points
take a parameter vector (in the basis of the fit, "RWA", or the physical
one) and return the objective and its terms
(``AZUREAPI::CalculateObjectiveRWA``), the gradient, the residuals and their
Jacobian, the model values and the observables of every point, the parameter
records (``AZUREAPI::GetParameterRecords``) and bounds. It evaluates with the
same functions as the program: the point pass
(``WorkingCopy::CalculatePointsInParallel``), the model pass and
``Objective::Engine()``. A session's model is the project's: another channel
radius or nuclear potential is an edit of the project and a new session.
``AZUREAPI::FitWithEngine`` (``api/src/AZUREAPIFit.cpp``) runs a row of the
fitter table on the session with the program's functions
(``Session.fit_with_engine``; the binding ``api/src/pybind_fit.cpp`` releases
the GIL and adds ``StopRequest`` and ``fitters()``).

The module ``_azure2``
----------------------

The pybind11 module (``api/src/CMakeLists.txt``, target ``_azure2``) links
``AZUREAPI`` and the converter and is built into the package directory
``pyazr/``. One file per subject:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - File
     - What it binds
   * - ``pybind_api.cpp``
     - ``_azure2.Session``: a ``Config`` and an ``AZUREAPI`` per session,
       built from the project file and the options of the call
       (``RuntimeOptions``, ``apply_options``); every evaluation entry point.
   * - ``pybind_project.cpp``
     - ``_azure2.Document``: the lossless document ``io::Document`` that
       ``pyazr.Project`` reads and edits, and the project's records.
   * - ``pybind_nuclear_part.cpp``
     - The ``<levels>`` rows as levels and pairs, the ``<config>`` entries
       and the ``<thm>`` section, through the engine's grammars.
   * - ``pybind_schema.cpp``
     - The table of format 3 as a dict (the same document as
       ``AZURE2 --schema``).
   * - ``pybind_convert.cpp``
     - Validation, description and conversion of a project, as
       ``AZURE2 --validate``, ``--describe``, ``--convert``.
   * - ``pybind_describe.cpp``
     - What ``pyazr.describe()`` and ``pyazr.kinds()`` read: the build
       options, the Coulomb table, the session options
       (``NotProjectOptionRows()``).
   * - ``pybind_records.cpp``
     - One record per fit parameter and per pair (``ParameterFields()``),
       from which ``pyazr.Parameter`` and the parameter table are built.
   * - ``pybind_evaluation.cpp``
     - One evaluation as arrays (``AZUREAPI::EvaluateObjectiveRWA``,
       ``api/src/AZUREAPIObjective.cpp``): the objective, its terms by name,
       and each data set's points with the outputs of the table
       ``PointOutputs()``. ``pyazr.Evaluation`` (``pyazr/evaluate.py``) reads
       it; ``pyazr.Objective`` (``pyazr/objective.py``) is the fit's function
       with its gradient, its least-squares rows
       (``AZUREAPI::ObjectiveRowsRWA``) and the sampler's log-probability
       (``AZUREAPI::LogProbabilityRWA``).
   * - ``pybind_physics.cpp``
     - The table of physics functions (``BuildTable()``), each a row with its
       arguments, units, domain and an adapter to ``physics/``.
   * - ``pybind_stopping.cpp``
     - The submodule ``_azure2.stopping`` (``pyazr.stopping``), with the
       table of stopping functions (``BuildStoppingTable()``), each row
       with SRIM's citation.
   * - ``pybind_rows.h``
     - The row types the two tables share (``ArgSpec``, ``OutSpec``,
       ``ArgKind``).
   * - ``pybind_potentials.cpp``, ``pybind_model_data.cpp``, ``pybind_thm.cpp``
     - The view of a session's hybrid potentials, the model data, the THM
       tables and grammars.
   * - ``pybind_notices.cpp``
     - The engine's standard output of a session as notices, which
       ``pyazr`` turns into warnings (``pyazr/session_checks.py``).

The binding holds no grammar and no formula of its own: Python hands over
names and numbers, and the engine's code reads, writes and computes.

pyazr
-----

``pyazr.Project`` (``pyazr/project.py``, ``pyazr/nuclear_part.py``) is the
project as a document: it opens a format-3 file, edits records, levels, pairs
and ``<config>`` keys through the engine's grammars, and saves it; an
unedited project is written back byte for byte. It does not evaluate.

``pyazr.Session`` (``pyazr/session.py``) is a session: it owns one
``_azure2.Session`` and offers the parameters (``Parameters``), the data
sets, the objective (``Session.objective()``: a ``pyazr.Objective``, ``pyazr/objective.py``) and the
numbers at one parameter point (``Session.evaluate(x)``: a ``pyazr.Evaluation``, ``pyazr/evaluate.py``: model
values, χ², the objective's terms, the data sets' rows of the objective, the extrapolations in extrapolation mode),
the output files, and the THM
views (``pyazr/session_thm.py``). Its methods that need the session's pairs
and options are mixins of one job each (``session_physics.py``,
``session_thm.py``, ``potentials.py``), listed one per line in the class
statement. Every public method that enters the engine holds the session's
lock (``pyazr/locking.py``: ``engine_call``), as ``Objective`` and
``Evaluation`` do, so two threads on one session run one after the other;
``tests/pyazr/session_lock_test.py`` fails for a member that reaches the
engine without it. Fits, bands and model averages are Python on top of it
(``pyazr/bands.py``, ``pyazr/modelavg.py``), as is the model of a data set on
a fine grid (``Session.model_curve``, ``pyazr/model_curve.py``: a
``CurveSession`` owns a session of the extrapolations made from the data
sets' records, trimmed by ``bands``; ``refine_plan`` refines their grid with
the ENDF writer's rule, ``pyazr.endf.grid.start_grid`` and ``refine``).

``python -m pyazr`` (``pyazr/cli.py``) is one table of commands
(``COMMANDS``), each calling one ``pyazr`` function.

The web UI (``pyazr/ui``) is a server (``Server``, ``server.py``) that owns
the open project and its edit history (``State``), the jobs (``JobTable``)
and a session worker (``SessionWorker``): the session lives in a separate
process, where ``Engine`` (``engine.py``) owns a ``pyazr.Session``. Every route
of its API is one row of ``ROUTES`` (``routes.py``); the functions of the
session are pyazr's rows (``Session.describe()["functions"]``) with the page's
own part of each (``PAGE`` in ``physics.py`` and ``thm.py``), and the keys of
``<thm>`` are ``pyazr.kinds()["thm"]``. The page
(``pyazr/ui/static``) holds no state of the project and computes no physical
number.
