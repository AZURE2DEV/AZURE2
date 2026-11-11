Installation
============

AZURE2 is built with CMake. The command-line program and the ``pyazr``
Python module come from one source tree; the README of the repository has the
platform-specific package lists. This version has no graphical setup utility
(:doc:`../user_guide/graphical_interface`).

Dependencies
------------

**Required:**

- a C++ compiler with **OpenMP** (GCC or Clang; on macOS Homebrew's
  ``libomp``)
- **CMake** 3.16 or later, and ``pkg-config`` (CMake finds GSL through it)
- **GNU Scientific Library (GSL)** (https://www.gnu.org/software/gsl/)

Minuit2, the Coulomb-function library, the SRIM stopping-power utilities,
pugixml and the MCMC sampler are vendored under ``external/`` and
``numcmc/`` and built with AZURE2; nothing else has to be installed for
them.

**Optional:**

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

Platform scripts that set the options live in ``scripts/``
(``build_linux.sh``, ``build_windows.sh``, ``build_docker.sh``,
``build_crc.csh``).

5. Build the Coulomb-function table once (recommended):

   .. code-block:: bash

      ./src/AZURE2 --build-coulomb-table

   AZURE2 takes the Coulomb functions of charged-particle channels from this
   table when it finds it; the results are the same and a calculation is
   several times faster (about four times for the natural-neon example). The
   table is not distributed with AZURE2. The build takes about 5 minutes on
   four threads and writes 88 MB to ``~/.azure2/coulomb_table.bin`` (or to the
   file named by the environment variable ``AZURE2_COULOMB_TABLE``). Without
   it, AZURE2 says so once per run and computes the functions directly. See
   :doc:`../user_guide/configure_menu`.

CMake Options
^^^^^^^^^^^^^

The following options can be passed to CMake as ``-D<OPTION>=ON|OFF``:

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Option
     - Default
     - Description
   * - ``USE_API``
     - ON
     - Build the ``AZUREAPI`` library and the pybind11 ``_azure2`` module
       (``pyazr``).
   * - ``USE_MCMC``
     - ON
     - MCMC Bayesian sampling with the bundled ``numcmc``.
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
   * - ``CROSS_COMPILE_WINDOWS``
     - OFF
     - Cross-compile a Windows binary with MinGW.
   * - ``CODE_COVERAGE``, ``USE_GCOV``
     - OFF
     - Coverage builds.
   * - ``AZURE2_SANITIZE``
     - (empty)
     - ``thread``, ``address`` or ``address,undefined``: a sanitizer build of
       every target, for ``scripts/sanitize_upgrade.sh`` (developers only;
       ThreadSanitizer needs clang with LLVM's libomp and libarcher).

If GSL is installed in a non-standard location, point CMake at it
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
``ctest --test-dir build`` runs the unit, reference and ``pyazr`` tests.
:doc:`../developer/contributing` describes the conventions.

Docker
------

``scripts/build_docker.sh`` builds an image from
``packaging/docker/Dockerfile.azure2``; ``examples/run_mcmc.sh`` runs a
script in it. For HPC systems the image can be converted to Singularity/Apptainer:

.. code-block:: bash

   sudo apptainer build AZURE2.sif docker-daemon://azure2:latest
