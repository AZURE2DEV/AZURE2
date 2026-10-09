Contributing
============

A change to AZURE2 comes with a test that fails without it, and the existing
suites must pass on all three CI platforms. This page describes the test
layout, the conventions the scripts share, and what CI requires.

Test suites
-----------

``tests/run_tests.sh [path/to/AZURE2]``
   The physics regression suite. Every directory ``tests/<name>/`` holding
   ``<name>.azr``, its ``data/`` and ``expected/chiSquared.out`` is run
   (menu 1, calculate with data, from a clean ``output/``) and its total
   :math:`\chi^2`, the :math:`\chi^2` of every segment and every segment's
   point count are compared with the expected file: the counts exactly, the
   :math:`\chi^2` to the relative tolerance ``TOL`` (default
   :math:`10^{-3}`, which absorbs the spread between compilers and libm; a
   project may set its own in ``tests/<name>/tolerance``, with the reason
   next to it). Without an argument the newest ``build*/src/AZURE2`` is
   used. Afterwards every executable ``tests/*/check.sh`` is run with the
   binary as its argument, its output kept in ``check.log``; a
   ``check.sh`` that is not executable is reported and not run. Cases are
   found automatically: a new project or check needs no change to the
   runner.
``ctest --test-dir build``
   ``tests/reference`` (C++ checks against closed forms or independent
   implementations, several with a Python reference script beside them;
   its ``README.md`` lists what each one checks), ``tests/gui`` (headless
   Qt tests, ``QT_QPA_PLATFORM=offscreen``, built with the GUI) and
   ``tests/pyazr`` (``pyazr_*``, the Python module against the engine and
   the CLI binary of the same build, passed as ``AZURE2_BIN``).

Shared shell helpers
--------------------

``tests/lib/guard.sh`` puts a time limit on every engine run, so that a hang
fails a check instead of stalling the suite. It is sourced by
``run_tests.sh`` and the ``check.sh`` scripts:

.. code-block:: bash

   . "$REPO/tests/lib/guard.sh"
   RUN="$(guard_command 600)"
   printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr

``guard_command SECONDS`` gives GNU ``timeout``, Homebrew's ``gtimeout`` or,
where neither exists, the shell watchdog ``run_guard``, so the limit holds on
every platform. ``run_tests.sh`` allows ``TEST_TIMEOUT`` seconds per project
(default 1800).

``tests/lib/check_common.sh`` holds the helpers the ``check.sh`` scripts
share. A script sets ``fail=0`` and, for ``ran`` and ``model``, ``$WORK``
(one directory per run) and ``$OUT`` (the output file of the segment):

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - helper
     - does
   * - ``ok MSG``
     - prints a passed line
   * - ``bad MSG``
     - prints a failed line and sets ``fail=1``
   * - ``same A B``
     - both files exist and are equal up to ``\r``
   * - ``ran NAME``
     - the run in ``$WORK/NAME`` ended with status 0 and wrote ``$OUT``; else a failure with the end of its log
   * - ``model NAME``
     - ``E model`` per point of ``$OUT`` (columns 1 and 4)
   * - ``counted N``
     - ``N`` is a positive number, i.e. a comparison compared something

A check must not pass for the wrong reason: test that a run succeeded
(``ran``) before comparing its output, and that a comparison found rows
(``counted``) before reporting it equal. Strip ``\r`` before counting fields
(Windows line ends), and use portable ``sed``/``awk`` (macOS has BSD tools).

Skipped tests: exit 77
----------------------

A ``pyazr`` test that cannot run its engine part -- numpy, scipy or mpmath
missing, the ``_azure2`` module not built (no pybind11 at configure time), no
``AZURE2_BIN`` -- runs its pure-Python checks and then exits with status 77.
``tests/pyazr/CMakeLists.txt`` sets ``SKIP_RETURN_CODE 77`` on every test
registered above that block, so ctest lists such a test as *Skipped*, not
*Passed*. The block must stay last in the file; a new test only has to
follow the convention. An exit of 0 there once hid that CI had never run a
single CLI-against-pyazr comparison.

What CI requires
----------------

``.github/workflows/build.yml`` builds on Linux, macOS and Windows with
pybind11, numpy, scipy and mpmath installed, so that ``_azure2`` is built
and every ``pyazr`` engine test runs. Each job runs a smoke test,
``tests/run_tests.sh`` and ``ctest``; every skipped test is listed in a
warning, and a skipped ``pyazr_*`` test fails the job on a platform that
must run them. ``.github/workflows/docs.yml`` builds this documentation with
``sphinx-build -W --keep-going``: a broken reference or a page missing from
a toctree fails it.

Pins and platforms
------------------

A recorded value (``expected/chiSquared.out``, a number in a ``check.sh``)
is taken from a run you trust and changed only with the reason in the commit.
A local build may use instructions the CI runners do not (FMA on an
x86-64-v3 build), so quantities sensitive to round-off can differ in
the last digits between a local run and CI; pin with a tolerance that covers
that, and compare results of several sessions in one process to
:math:`10^{-12}` rather than bit for bit (the Coulomb-function memo takes
energies within :math:`10^{-12}` MeV as equal). Run one engine process at a
time on a small machine: a realistic THM session takes 0.2-2 GB.
