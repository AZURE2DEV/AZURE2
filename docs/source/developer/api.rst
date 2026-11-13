The Binding and pyazr
=====================

How Python reaches the engine. The user's side is in the user guide
(``pyazr``); this page is the structure behind it.

``AZUREAPI``
------------

``AZUREAPI`` (``api/include/AZUREAPI.h``, ``api/src/AZUREAPI.cpp``) is the
engine as a session: it owns a ``CNuc`` and a ``WorkingCopy`` (made in
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
``Objective::Engine()``. A session can rebuild its model with another channel
radius (``AZUREAPI::SetRadius``) or another nuclear potential.

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
   * - ``pybind_physics.cpp``
     - The table of physics functions (``BuildTable()``), each a row with its
       arguments, units, domain and an adapter to ``physics/``.
   * - ``pybind_stopping.cpp``
     - The submodule ``_azure2.stopping`` (``pyazr.stopping``).
   * - ``pybind_potentials.cpp``, ``pybind_model_data.cpp``, ``pybind_thm.cpp``
     - The hybrid potentials of a session, the model data, the THM tables
       and grammars.
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

``pyazr.azure2`` (``pyazr/azure2.py``) is a session: it owns one
``_azure2.Session`` and offers the parameters (``ParameterSet``), the data
sets, the objective and its gradient, residuals, observables and
extrapolations, the output files, and the THM views. Fits, bands and model
averages are Python on top of it (``pyazr/bands.py``, ``pyazr/modelavg.py``).

``python -m pyazr`` (``pyazr/cli.py``) is one table of commands
(``COMMANDS``), each calling one ``pyazr`` function.

The web UI (``pyazr/ui``) is a server (``Server``, ``server.py``) that owns
the open project and its edit history (``State``), the jobs (``JobTable``)
and a session worker (``SessionWorker``): the session lives in a separate
process, where ``Engine`` (``engine.py``) owns a ``pyazr.azure2``. Every route
of its API is one row of ``ROUTES`` (``routes.py``). The page
(``pyazr/ui/static``) holds no state of the project and computes no physical
number.

.. note:: **Planned.**
   ``pyazr`` 3.0 replaces the class ``azure2`` by ``pyazr.Session``, with
   ``Objective`` and ``Evaluation`` objects in place of the ``calculate_*``
   methods and ``Parameters`` in place of ``ParameterSet``; ``engine.py`` is
   then the only file of the web UI that changes.
