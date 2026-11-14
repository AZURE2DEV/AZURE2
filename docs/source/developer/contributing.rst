Contributing
============

AZURE2's results are used in publications and nuclear-data evaluations, so the
standard for a change is that it can be *shown* to be correct, not only that
it compiles and runs. This page is for people. The same rules, written for AI
coding assistants, are in ``CLAUDE.md`` at the repository root, with shorter
files in ``tests/`` and ``pyazr/``.

Branches
--------

* ``main`` and ``qt5`` (the GitHub default) are maintained by people only.
* ``dev`` is the integration branch. The documentation site and the rolling
  pre-release downloads are built from it.
* Each feature or fix is developed on its own branch cut from ``dev``:
  ``feature/<topic>``, ``fix/<topic>``, ``docs/<topic>``, ``test/<topic>`` or
  ``chore/<topic>``. It is merged back into ``dev`` once the checks below pass.

Building and checking
---------------------

.. code-block:: bash

   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
         -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
   cmake --build build --parallel
   ./scripts/check.sh            # consistency, build, ctest, physics suite, docs
   ./scripts/check.sh --quick    # without the physics suite and docs

``scripts/check.sh`` stops at the first failing stage. ``JOBS=<n>`` sets both the
compile jobs of the build and the number of tests ``ctest`` runs at once
(default 4; most tests use two OpenMP threads each). CI runs the same tests
on Linux, macOS and Windows, and a change is not finished until all three
pass. The test suites are:

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Suite
     - What it checks
   * - ``tests/<name>/``
     - Physics regression: each project's total and per-segment
       :math:`\chi^2` and point counts against ``expected/chiSquared.out``
       (relative tolerance :math:`10^{-3}`), run by ``tests/run_tests.sh``. A
       ``check.sh`` beside it covers what one :math:`\chi^2` cannot.
   * - ``tests/reference/``
     - Quantities computed a second way, from closed forms or independent
       implementations. These are the tests that show a number is *right*,
       not only unchanged.
   * - ``tests/pyazr/``
     - The ``.azr`` round trip through pyazr, and pyazr against the
       engine. Tests labelled ``slow`` run with ``ctest -L slow``.
   * - ``tests/consistency/``
     - Command-line flags against ``--help`` and :doc:`../reference/command_line`,
       the segment keywords of the engine and pyazr, test registration,
       and ``check.sh`` permissions.

What a finished change contains
-------------------------------

1. **A test that fails without it.** A bug fix comes with the test that would
   have caught the bug. A new physics feature is checked against something
   independent of the code: a closed form, an independent calculation (for
   example ``mpmath`` or a magnetic-substate sum, as in
   ``tests/ec_m1_exact/truth.py``), a published value, or an exact identity.
2. **Every front end.** A user-visible feature reaches the engine, the
   command line and pyazr together (see the checklist below).
3. **Documentation.** The page of this site that describes the behaviour is
   updated in the same branch. Claude Code skills in ``.claude/skills`` that
   describe the old behaviour are corrected; skills teach use and link this
   site for behaviour, and ``tools/check_skills.py`` checks their rules
   (``docs/upgrade/RULES.md`` §3), among them that no skill names a pyazr
   name the migration table of :doc:`../user_guide/pyazr` lists as removed.
   A skill's scripts are run by its ``scripts/selftest.py`` on a project of
   ``tests/``; ctest runs each as ``skill_<name>``
   (``tests/skills/CMakeLists.txt``): a change of pyazr or the engine that
   breaks a skill fails there.
4. **A commit message that carries the evidence.** It says what was wrong, how
   it was tested, whether the test fails on the old build, and the numbers
   before and after where results moved.

Adding an ``.azr`` keyword or option
------------------------------------

The project is read once, through ``io/`` and the format-3 table of
``experiment/`` (:doc:`architecture`); a new key is a row of that table, not
parser code. :doc:`extending` lists the places for a key or a record kind, an
option of the calculation, a command-line flag, a parameter kind and an
objective term, each with the test that fails when one is forgotten.

Rules for test references
-------------------------

* An ``expected/`` file, a tolerance or a reference value is never changed to
  make a test pass. It is changed only for an intended and understood physics
  change, in a separate commit that gives the old and new values, the reason,
  and the independent evidence for the new value.
* Results must not depend on the number of threads, the platform beyond about
  :math:`10^{-5}`, or the order of levels and channels.
* If a change alters what a cached file would contain (external-capture
  integrals, energy tables, the Coulomb table), the cache's signature is
  changed so that old files are recomputed rather than reused.
* Input files written by earlier versions must keep loading and give the same
  results. New keywords are optional, and their defaults reproduce the old
  behaviour.

Style
-----

* C++: ``.clang-format`` (two-space indent, no automatic line wrapping).
  Comments explain *why*, citing the equation or paper. New features go in
  their own class and file. ``external/`` is upstream code and is not edited.
* Python: Python 3.10 or newer, with numpy, scipy and mpmath as the only core
  dependencies. Every public name in ``pyazr.__all__`` has a docstring giving
  units and frames. Removing or renaming a public name requires a deprecation
  period.
* Versions: pyazr's version is ``__version__`` in ``pyazr/__init__.py``
  (semantic versioning). The engine's version is set in ``CMakeLists.txt`` and
  bumped when a release is tagged.

Test infrastructure
-------------------

How the suites above are run, the conventions their scripts share, and what
CI requires.

Running the suites
^^^^^^^^^^^^^^^^^^

``tests/run_tests.sh [path/to/AZURE2]``
   The physics regression suite. Every directory ``tests/<name>/`` holding
   ``<name>.azr``, its ``data/`` and ``expected/chiSquared.out`` is run
   (``AZURE2 calculate``, from a clean ``output/``) and its total
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
   its ``README.md`` lists what each one checks) and
   ``tests/pyazr`` (``pyazr_*``, the Python module against the engine and
   the CLI binary of the same build, passed as ``AZURE2_BIN``).

Shared shell helpers
^^^^^^^^^^^^^^^^^^^^

``tests/lib/guard.sh`` puts a time limit on every engine run, so that a hang
fails a check instead of stalling the suite. It is sourced by
``run_tests.sh`` and the ``check.sh`` scripts:

.. code-block:: bash

   . "$REPO/tests/lib/guard.sh"
   RUN="$(guard_command 600)"
   $RUN "$AZURE2_BIN" calculate run.azr

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

A script that sets ``pipefail`` (most do) never pipes into a reader that
leaves before the end of its input: ``grep -q``, ``grep -m``, ``head -n``,
``sed`` with ``q``, ``awk`` with ``exit`` outside ``END``. The reader leaves
at its first match, the writer's next write gets ``SIGPIPE``, and the
pipeline's status becomes 141, so on a loaded machine a match reads as "not
found" (and, under ``!``, a miss as a pass) at random. Search a variable with
a here-string, ``grep -q PATTERN <<< "$text"``; capture a command's output
first, ``out="$(cmd)"``, or let the reader read to the end (``sed -n 1p``
for the first line, a flag instead of ``exit`` in ``awk``). ``head -c N`` is
a cap and allowed. Check 50 of ``tests/consistency/repo_consistency_test.py``
finds such a pipe.

A comparison's failure must reach the caller. An ``awk`` program that prints a
mismatch mark (``grid``, ``x``) prints the mark alone: ``END`` tests the flag
the body set before it prints a number (``END { if (e) exit; printf ... }``),
or the caller reads two lines and never sees the mark. A value handed to
``awk`` with ``-v`` and compared with an upper bound is first tested to be a
number (``v ~ /^ *[-+]?\.?[0-9]/ && v + 0 <= tol``): ``awk`` compares an
empty value or a mark as text, or after ``v + 0`` as 0, and it passes. Check
56 (``tests/consistency/comparison_scan.py``) finds both forms in
``tests/*/check.sh`` and ``tests/lib/``. A tolerance read from a file is
refused when it is not a number (``tests/run_tests.sh``).

Skipped tests: exit 77
^^^^^^^^^^^^^^^^^^^^^^

A ``pyazr`` test that cannot run its engine part -- numpy, scipy or mpmath
missing, the ``_azure2`` module not built (no pybind11 at configure time), no
``AZURE2_BIN`` -- runs its pure-Python checks and then exits with status 77.
``tests/pyazr/CMakeLists.txt`` sets ``SKIP_RETURN_CODE 77`` on every test
registered above that block, so ctest lists such a test as *Skipped*, not
*Passed*. The block must stay last in the file; a new test only has to
follow the convention. An exit of 0 there once hid that CI had never run a
single CLI-against-pyazr comparison.

What CI requires
^^^^^^^^^^^^^^^^

``.github/workflows/build.yml`` builds on Linux, macOS and Windows with
pybind11, numpy, scipy and mpmath installed, so that ``_azure2`` is built
and every ``pyazr`` engine test runs. Each job runs a smoke test,
``tests/run_tests.sh`` and ``ctest``; every skipped test is listed in a
warning, and a skipped ``pyazr_*`` test fails the job on a platform that
must run them. ``.github/workflows/docs.yml`` builds this documentation with
``sphinx-build -W --keep-going``: a broken reference or a page missing from
a toctree fails it.

Pins and platforms
^^^^^^^^^^^^^^^^^^

A recorded value (``expected/chiSquared.out``, a number in a ``check.sh``)
is taken from a run you trust and changed only with the reason in the commit.
A local build may use instructions the CI runners do not (FMA on an
x86-64-v3 build), so quantities sensitive to round-off can differ in
the last digits between a local run and CI; pin with a tolerance that covers
that. Several sessions in one process give the bits of a fresh process (there
is no Coulomb-function memo since ``docs/upgrade/coulomb.md``; builds before it
agreed only to :math:`10^{-12}`). Run one engine process at a
time on a small machine: a realistic THM session takes 0.2-2 GB.
