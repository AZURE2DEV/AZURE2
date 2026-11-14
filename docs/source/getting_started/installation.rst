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
   * - ``BUILD_TESTING``
     - ON
     - Add the test directories under ``tests/`` (off in a wheel build).
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

``pyazr`` needs Python 3.10 or later and NumPy 1.24 or later; SciPy comes with the
extras ``web`` (the web page's scipy fit job) and ``examples``. One command
builds and installs it with its engine from a checkout:

.. code-block:: bash

   pip install .             # from the repository root: compiles the engine (a few minutes)
   pip install ".[all]"      # plus what the examples use (matplotlib, emcee, iminuit, scipy)
   pip install ".[web]"      # plus the web page's server (fastapi, uvicorn, scipy)

``pyproject.toml`` builds with scikit-build-core, which runs CMake with the
option set of a wheel (``BUILD_TESTING`` off) and
builds two targets. The compilers, CMake 3.16 or later, ``pkg-config`` and
GSL listed above must be installed; pip fetches
scikit-build-core and pybind11 itself. ``pip wheel . -w dist`` writes the
wheel instead of installing it. The wheel holds:

- the package ``pyazr`` with its engine module ``_azure2`` (the engine, Minuit2
  and the SRIM-2013 data compiled in);
- its data files: the AME2020 mass table and the web page's files (with the
  vendored plotly);
- two commands on the ``PATH``: ``pyazr`` (:doc:`../user_guide/pyazr_command_line`;
  the same as ``python -m pyazr``) and the program ``AZURE2`` (its commands,
  :doc:`../reference/command_line`);
- the licence files: AZURE2's (GPL-3.0), the notice of the SRIM-2013 data
  files, plotly's.

GSL and the OpenMP runtime stay shared libraries of the machine (a wheel built
for distribution has them copied in by ``auditwheel`` or ``delocate``). The
examples of ``pyazr/examples/`` are not installed: they read the repository's
test projects.

**SRIM.** The stopping powers of ``pyazr.stopping`` and of the engine's
target effects are SRIM-2013's (J. F. Ziegler, M. D. Ziegler and J. P.
Biersack, Nucl. Instrum. Methods Phys. Res. B 268 (2010) 1818), whose data
files are compiled into the module under their authors' terms, not AZURE2's
licence. A result that uses them cites SRIM; ``pyazr.describe()["stopping"]``
and ``pyazr stopping`` print the citation.

**The Coulomb table** is not in the wheel (88 MB). Get it once:

.. code-block:: bash

   pyazr coulomb-table where     # where the engine looks, and whether a file is there
   pyazr coulomb-table fetch     # download the published table, checked by its SHA-256
   pyazr coulomb-table build     # or build it here, offline (about 5 minutes on 4 threads)

Both write the file the engine reads by default (``$AZURE2_COULOMB_TABLE``,
else ``~/.azure2/coulomb_table.bin``), so every later session and run of the
program uses it with nothing set. ``fetch`` refuses and removes a file whose
SHA-256 is not the published one; until a table is published it says so and
names ``build``. Without a table the engine computes the Coulomb functions
with its exact routines (the same results, about four times slower).
In Python: ``pyazr.coulomb_table`` (``where``, ``fetch``, ``build``).

**Developers** keep the in-tree build of the previous section: ``cmake --build
build`` writes ``_azure2`` into ``pyazr/`` (or into ``PYAZR_OUTPUT_DIR``, found
through ``PYAZR_ENGINE_DIR``), and ``import pyazr`` from the repository's root
loads it; every test uses this module. ``pip install .`` refuses to start
while such a module lies in ``pyazr/`` (the wheel would carry it): remove it,
or let the in-tree build write elsewhere with ``-DPYAZR_OUTPUT_DIR``. See
:doc:`../user_guide/pyazr`.

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
