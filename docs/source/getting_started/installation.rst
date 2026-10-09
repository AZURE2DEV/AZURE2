Installation
============

AZURE2 is built with CMake. The command-line engine, the Qt graphical setup
utility and the ``pyazr`` Python module come from one source tree; the
README of the repository has the platform-specific package lists and the
ready-built downloads.

Dependencies
------------

**Required:**

- a C++ compiler with **OpenMP** (GCC or Clang; on macOS Homebrew's
  ``libomp``)
- **CMake** 3.16 or later, and ``pkg-config`` (CMake finds GSL and Qwt
  through it)
- **GNU Scientific Library (GSL)** (https://www.gnu.org/software/gsl/)

Minuit2, the Coulomb-function library, the SRIM stopping-power utilities,
pugixml and the MCMC sampler are vendored under ``external/`` and
``numcmc/`` and built with AZURE2; nothing else has to be installed for
them.

**Optional:**

- **Qt5** (``Core``, ``Widgets``, ``Svg``, ``Script``) -- the graphical
  interface (``BUILD_GUI``)
- **Qwt** built for Qt5 -- the Plot tab (``USE_QWT``; the THM Workspace
  draws its plots without it)
- **Readline** -- line editing at the command-line prompts
  (``USE_READLINE``)
- **pybind11** -- the ``_azure2`` module that ``pyazr`` drives
  (``USE_API``); without it CMake warns and builds no module

Building from Source
--------------------

.. code-block:: bash

   git clone https://github.com/rdeboer1/AZURE2.git
   cd AZURE2
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build --parallel 2

The executable is ``build/src/AZURE2``; the ``_azure2`` module is placed in
``pyazr/``. They are separate targets: building ``AZURE2`` alone does not
rebuild ``_azure2``, and a stale module reproduces engine behaviour that the
current source no longer has. Build both (the default target does) after an
engine change. On a machine with little memory keep the parallel jobs low.

Platform scripts that detect Qwt and set the options live in ``scripts/``
(``build_linux.sh``, ``build_macos.sh``, ``build_windows.sh``,
``build_docker.sh``).

CMake Options
^^^^^^^^^^^^^

The following options can be passed to CMake as ``-D<OPTION>=ON|OFF``:

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Option
     - Default
     - Description
   * - ``BUILD_GUI``
     - ON
     - Build the graphical setup utility.
   * - ``USE_QWT``
     - ON
     - Build the Plot tab (needs Qwt).
   * - ``USE_API``
     - ON
     - Build the ``AZUREAPI`` library and the pybind11 ``_azure2`` module
       (``pyazr``; the THM Workspace's Diagnostics page needs it too).
   * - ``USE_MCMC``
     - ON
     - MCMC Bayesian sampling with the bundled ``numcmc``.
   * - ``USE_ERYA``
     - ON
     - SRIM stopping-power utilities.
   * - ``USE_READLINE``
     - ON
     - Readline for console input.
   * - ``USE_STAT``
     - ON
     - Use ``stat()`` for directory checks (off on Windows).
   * - ``USE_NLOPT``
     - OFF
     - NLopt as an alternative minimizer (expects an ``nlopt/`` tree).
   * - ``BUILD_LIBRARY``
     - OFF
     - Build AZURE2 as a library.
   * - ``BUILD_MACOS_BUNDLE``
     - OFF
     - Build a macOS ``.app`` bundle.
   * - ``CROSS_COMPILE_WINDOWS``
     - OFF
     - Cross-compile a Windows binary with MinGW.
   * - ``CODE_COVERAGE``, ``USE_GCOV``
     - OFF
     - Coverage builds.

If GSL or Qwt are installed in a non-standard location, point CMake at them
with ``-DCMAKE_PREFIX_PATH=/path/to/prefix``.

The Python package
------------------

``pyazr`` needs Python 3.8 or later with NumPy, SciPy and mpmath:

.. code-block:: bash

   pip install -e .          # from the repository root
   pip install -e ".[all]"   # plus what the examples use (matplotlib, emcee, zeus)

``import pyazr`` from the repository root loads the ``_azure2`` module that
CMake built into ``pyazr/``. See :doc:`../user_guide/pyazr`.

Tests
-----

``tests/run_tests.sh [path/to/AZURE2]`` runs every project in ``tests/``
against its recorded :math:`\chi^2` and every ``tests/*/check.sh``;
``ctest --test-dir build`` runs the unit, GUI and ``pyazr`` tests.
:doc:`../developer/contributing` describes the conventions.

Docker
------

``scripts/build_docker.sh`` builds an image from
``packaging/docker/Dockerfile.azure2``; ``examples/run_gui.sh`` starts the GUI
in it. For HPC systems the image can be converted to Singularity/Apptainer:

.. code-block:: bash

   sudo apptainer build AZURE2.sif docker-daemon://azure2:latest
