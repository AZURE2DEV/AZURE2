# AZURE2 — guidance for Claude Code

AZURE2 is a multi-level, multi-channel R-matrix code for nuclear reactions and
scattering. Physicists use its numbers in publications and evaluations, so a
wrong number that looks right is the worst possible bug. Optimise for
**correctness you can show**, then for keeping the GUI, CLI and Python API in
step, then for speed.

Human-readable versions of this file: `docs/source/developer/contributing.rst`
and `docs/source/developer/architecture.rst`. Directory-specific rules are in
`tests/CLAUDE.md`, `pyazr/CLAUDE.md` and `gui/CLAUDE.md`.

## Upgrade (branch `feature/upgrade`)

The engine computes one chain: level parameters → collision matrix → observables → apparatus fold → χ². The
upgrade rebuilds the code so that each stage of that chain is one module that depends only on the stages before it,
and so that every change can be checked against synthetic projects with known parameters and against the present
engine's results. Why, and in what order: `docs/upgrade/README.md`.

This branch carries the staged rebuild described in `docs/upgrade/` (start with `docs/upgrade/HANDOFF.md`). On it,
and on the phase branches cut from it (`feature/upgrade-<phase>-<topic>`), these rules apply on top of everything
below (full text and reasons in `docs/upgrade/RULES.md`):

- Code is shaped like the physics: `physics/`, `experiment/`, `inference/`, `io/`, plus `support/` for plumbing.
  Modules use only the modules below them; each is its own CMake target, so a wrong dependency does not compile.
- The engine computes three lines (`docs/upgrade/ARCHITECTURE.md` §1). New physics enters at one extension point,
  changes one symbol, and appears in one of the test projects T1–T3 (`tests/upgrade/`).
- Setup state is read-only once a fit stage starts; parallel code writes only its own slice or scratch; no caches
  beyond the per-stage tables; results independent of thread count.
- One implementation per formula and format, in C++ (Coulomb functions from pycoulomb). Library readers live in pyazr
  `sources/`; the core reads only the project. No network in the engine; imported numbers carry provenance.
- One internal angular-momentum and phase convention; conversions only in `level_parameters`.
- Every physics function is exposed in the binding with array inputs, so the web UI can plot it.
- Each phase passes the gate in `docs/upgrade/ROADMAP.md` and deletes the code it replaces. Update
  `docs/upgrade/HANDOFF.md` and the roadmap ticks at the end of every session.
- Documentation follows `docs/upgrade/STYLE.md`: a reaction-theory register for experimentalists; every symbol
  defined with its unit and frame; conventions and approximations stated with their validity; numbers, not adjectives.
- Golden outputs of the present engine live in `tests/upgrade/golden/` (generated files allowed there, like
  `tests/*/expected/`); regenerate them only with `tests/upgrade/golden/golden.py --write` and a commit that states why.
- Commit small and often (each completed step, at least every hour of work), with HANDOFF updated, so the work can be
  continued from another machine at any point.
- THM follows `docs/upgrade/THM.md`: reproduce published analyses first, HOES and OES data kept apart, resolution in
  c.m., no energy shifts to move THM data, no negative backgrounds, DWBA vertices validated against FRESCO.

## Branches (hard rule)

- **`main` and `qt5` (the GitHub default) are human-only.** Never commit to,
  merge into, rebase, or push either of them, and never open a PR that targets
  them unless the user says so in this conversation.
- **`dev` is the integration branch.** Development happens on one branch per
  feature or fix, cut from an up-to-date `origin/dev`, named
  `feature/<topic>`, `fix/<topic>`, `docs/<topic>`, `test/<topic>` or
  `chore/<topic>`. Merge it back into `dev` only once the checks under
  "Definition of done" pass. Use `git merge --no-ff` with a one-line summary,
  as in `Merge fix/park-dsde: <what it changes>`.
- Bring `dev` into a long-lived feature branch by merging, not rebasing. Other
  sessions build on pushed branches.
- One topic per branch. A fix discovered along the way gets its own branch, or
  at least its own commit, so it can be reverted or cherry-picked alone.
- Before starting, `git fetch` and check how far behind `origin/dev` you are;
  `dev` moves many times a day. Use a worktree (`git worktree add`) rather than
  switching branches under someone else's build.
- **CI is part of the check.** Linux passing locally does not mean macOS
  (Apple clang, libomp) and Windows (MinGW) pass. After a push, look at the
  run: `gh run list --branch <branch> --limit 3`, and
  `gh run view <id> --log-failed`. A red `dev` is the first thing to fix,
  before anything new is merged on top of it. The pre-release downloads users
  rely on are published only from a green Windows build of `dev`.

## Layout

| Path | What | Notes |
|---|---|---|
| `src/`, `include/` | C++ engine | `AZURE2.cpp` (main, CLI flags) → `AZUREMain` → `CNuc` (levels, channels) + `EData` (segments, data, experimental effects) → `AZURECalc` (χ²), `AZUREGrad` (gradients), `AZURECalcMCMC` |
| `api/` | `AZUREAPI` + pybind11 module `_azure2` | Builds `pyazr/_azure2*.so` in place |
| `pyazr/` | Python package | In-process engine plus a pure-Python `.azr` reader/writer (`azrfile.py`, `datasets.py`) |
| `gui/` | Qt5 + Qwt setup utility | Its own `.azr` reader/writer (`AZURESetup.cpp`, per-tab models) |
| `numcmc/` | Ensemble MCMC sampler | |
| `external/` | Vendored Minuit2, coul, erya/pugixml | **Do not edit**; it is upstream code |
| `tests/` | All tests | See `tests/CLAUDE.md` |
| `docs/source/` | Sphinx site, built with `-W` | The user-facing source of truth |
| `.claude/skills/` | How to *use* AZURE2 | See "Skills" below; `tools/check_skills.py` |

## Build and verify

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"   # without pybind11 the pyazr module is silently skipped
cmake --build build --parallel
./scripts/check.sh            # everything below, in order; --quick skips the physics suite
```

What `scripts/check.sh` runs, and what to run by hand when iterating:

```bash
ctest --test-dir build -LE slow --output-on-failure   # unit, reference, GUI, pyazr, consistency (minutes)
./tests/run_tests.sh build/src/AZURE2                 # physics regression: χ² of every tests/<name>/ (tens of minutes)
ctest --test-dir build -R <name> --output-on-failure  # a single test
ctest --test-dir build -L slow                        # long pyazr tests: run when you change pyazr fitting/bands
sphinx-build -W --keep-going -b html docs/source docs/_build/html   # docs, warnings are errors (as CI)
```

`tests/run_tests.sh` without an argument takes the *newest* `build*/src/AZURE2`.
Pass the binary explicitly if more than one build directory exists. The
Coulomb table (`--build-coulomb-table`, about 5 min, 88 MB in `~/.azure2/`)
makes runs about 4× faster and does not change results.

## Definition of done

A change is finished when all of these are true, and the commit message says
how each was checked:

1. **It builds** with the GUI and pyazr enabled, and adds no compiler warnings
   in the files it touches.
2. **The tests pass**: `ctest -LE slow` and `tests/run_tests.sh`. If a
   reference χ² moved, see the next section first.
3. **It is tested.** A bug fix comes with a test that fails on the old build
   and passes on the new one, and you say so. A new physics feature is checked
   against something independent of this code: a closed form, an
   `mpmath`/m-scheme/brute-force `truth.py`, a published number, or an
   identity such as a constant cross section surviving a fold. Agreeing with
   itself is not a check.
4. **Every front end agrees.** A user-visible feature reaches the engine, the
   GUI, pyazr and the docs in the same branch (see the checklist below). If one
   front end is deliberately left out, say so in the docs.
5. **Docs and skills are true.** The Sphinx page that describes the behaviour
   is updated in the same branch. A skill that states the old behaviour is
   corrected too (search `.claude/skills/` for the keyword).

## Rules that keep the results trustworthy

- **Never change an `expected/` file, a tolerance or a reference value to make
  a test pass.** Re-baseline only when the physics change is understood and
  intended. Do it in its own commit, which states the old and new values
  (`chi2 43,081.9 -> 43,065.4`), why they moved, and what independent check
  shows the new value is right. A per-test `tolerance` file needs a written
  reason inside it.
- **Do not weaken, skip or delete a failing test** to get green. Find out why
  it fails. If it cannot be fixed now, stop and tell the user.
- **Cached results carry a signature.** If a change alters what a cache
  (`intEC.dat`, `intEC.extrap`, EC tables, the Coulomb table, parameter files)
  would contain, bump its signature or version (e.g. `EData::ECSignature`) so
  old files are recomputed, not silently reused. Add a stale-file case to the
  test.
- **Old inputs keep working.** `.azr` files, `param.sav`/`param.par`/`param.fit`
  and data files from earlier versions must still load. A new keyword is
  optional and its default reproduces the old result. A format change gets a
  reader for the old form, a message, and a test with an old-format fixture.
- **No silent fallbacks.** If an input is invalid, or a requested method cannot
  be used, AZURE2 says so (and stops if the result would be wrong). It never
  continues with a guess. An unparseable line once spun forever at 100 % CPU,
  and a misread file once looped to `bad_alloc`.
- **Frames and units are explicit.** Input data are lab-frame, forward
  kinematics; outputs are centre-of-mass. Energies are in MeV, masses are atomic
  masses in u (not mass numbers). In `<levels>`, open channels carry Γ in
  eV, not a reduced-width amplitude, and closed channels carry an ANC in
  fm^-1/2. Name the
  frame and unit in every new variable, keyword, column header and docstring
  where it could be ambiguous. The working list of frames and units is in the
  `azure2-check` skill.
- **Determinism.** Results must not depend on thread count, on which file a
  session started from, or on the order of levels and channels. Tests run with
  `OMP_NUM_THREADS=2`. Differences across GCC, Clang and MinGW stay near 1e-5.

## A new `.azr` keyword or option touches all of these

The `.azr` file is parsed independently by the engine, the GUI and pyazr, so a
feature is not finished until each one reads it, writes it back unchanged and
reports it the same way:

1. Engine: `include/SegLine.h` / `NucLine.h` / `Config.cpp` (parse), the
   owning class (`ESegment`, `EData`, `CNuc`, ...), and start-up error
   messages for invalid combinations.
2. If it adds fit parameters or penalties, all four objective paths:
   `src/AZURECalc.cpp`, `src/AZUREGrad.cpp`, `src/AZUREMain.cpp` (output) and
   `api/src/AZUREAPI.cpp`. They must give the same χ², and the analytic
   gradient must agree with finite differences.
3. CLI: the flag in `parseOptions`, `printHelp`, and
   `docs/source/reference/command_line.rst`. The consistency test enforces all
   three.
4. GUI: read, write and edit it (`gui/include/SegmentKeywords.h` for segment
   tails), plus a round-trip test in `tests/gui/`.
5. pyazr: `datasets.py` / `azrfile.py` read and write, and `azure2.py` when it
   affects χ² or penalties, plus a test in `tests/pyazr/` registered in its
   `CMakeLists.txt`.
6. Docs: the user-guide page and `reference/data_formats.rst` or
   `output_files.rst`. Then the affected skills.
7. Tests: a `tests/<name>/` project (`expected/chiSquared.out` and/or
   `check.sh`) with a README that says what is checked and against what.

Reference implementations to copy: `git show --stat 062b1bfe` (segment keyword
`shiftfngroup`) and `3e13fb28` (`normgroup`).

## Code structure

- Put a new feature in its **own class and file**, as `ShftFunc`,
  `CoulombTable`, `DataCovariance` and `ScatteringLength` are. Do not grow the
  large files further: `EData.cpp` (4300 lines), `CNuc.cpp`, `AMatrixFunc.cpp`,
  `EPoint.cpp`, `AZURECalc.cpp`, `pyazr/azure2.py`, `pyazr/azrfile.py`.
  Moving code out of them is welcome, but do it in its own commit with no
  behaviour change and the full test suite run, never mixed into a feature.
- Use **one definition** of each formula and convention, and call it from both
  the engine and the gradient. When pyazr needs it, expose it through
  `_azure2` rather than re-implementing it in Python; where a Python copy is
  unavoidable, test it against the engine.
- Format C++ with `.clang-format` (2-space indent, `ColumnLimit: 0`, so you
  break lines yourself). Match the surrounding comment density. Comments say
  *why* and cite the equation or paper (`Edmonds 7.1.7`, `PRC 104, 064612`).
- Add no new dependencies (C++ or Python core) without asking. Optional Python
  extras go in `pyproject.toml` `[project.optional-dependencies]`.
- Diagnostics are environment variables named `AZURE_*` / `AZURE2_*` (for
  example `AZURE_PROFILE`, `AZURE_GRAD_CHECK`). Document new ones in
  `reference/command_line.rst`.

## Versions

- pyazr: `__version__` in `pyazr/__init__.py`, the only place it is stated.
  Bump the minor version for new public API, the patch version for fixes, and
  the major version for breaking changes, which also need a deprecation period.
  Keep `__all__` to deliberate public names.
- Engine: `AZURE2_VERSION_*` in `CMakeLists.txt` (2.1.0). It is independent of
  pyazr's version, and the docs take theirs from it.

## Commits

Follow the existing style. The subject line is `area: what changed`, or
`fix(area): ...`, `feat: ...`, `docs(skill): ...`, `test: ...`, `perf: ...`.
The body explains:

- what was wrong or missing, and the physics in one or two sentences;
- the test that shows it, including whether it fails on the old build;
- numbers before and after, where results moved;
- what was deliberately left out.

Commit generated outputs only under `tests/*/expected/`. Never commit `build*/`,
`output/`, `checks/`, `run.log`, `*.so`, or `.DS_Store`.

## Skills

The rules are `docs/upgrade/RULES.md` §3; `tools/check_skills.py` enforces items 4, 5, 6 and 8 and checks that
every docs link resolves (pre-commit and CI).

- **Three places, three jobs.** `.claude/skills/` teaches how to **use** AZURE2; the Sphinx docs (`docs/source/`)
  say how it **behaves**; this file says how to **change** it. When behaviour changes, the Sphinx page is updated
  first; a skill links to it and carries only workflow, recipes and pitfalls. A skill is never the only place a
  behaviour is documented.
- **Size.** A `SKILL.md` has at most 500 lines (hard limit 800). Detail goes in `references/*.md` of at most 300
  lines each.
- **No history in skills.** No commit hashes, job IDs, session narratives, quotes or "fixed on" notes. The only
  dated text is a version note for users of older builds, written as `builds before vX.Y` or
  `builds before YYYY-MM-DD`, with `expires: YYYY-MM-DD` at most 90 days ahead; the checker fails on an expired
  note, so it is removed or renewed.
- **One owner skill per task.** Other skills link to the owner by name instead of repeating it. Repository and
  user-level skills never share a name or a description.
- **Descriptions** are at most 600 characters, name their triggers, and say when to use a sibling skill instead.
- **No personal absolute paths** (`/groups/...`, `/Users/...`, home directories). Use environment variables
  (`$AZURE2_BINARY`, `$AZURE2_SRC`, `$RMFIT_HOME`, `$CAMPAIGN_ROOT`, rmfit's own
  `$AZURE2_ROOT`, `$FERDINAND_HOME`) and say what they hold.

The skill set: `azure2-setup`, `azure2-calculate`, `azure2-fit`, `azure2-check`, `azure2-uncertainty`,
`azure2-plot`, `azure2-targets`, `azure2-thm`, `azure2-resonance-search`, `azure2-exchange`, plus
`rmfit-campaign` (the table in `docs/upgrade/RULES.md` §3 says which task each owns). Run
`python3 tools/check_skills.py` before committing a skill change.

## When to stop and ask

- A reference value or test has to change, and you cannot point to an
  independent check that the new value is right.
- A change would alter results for existing users' projects by default.
- A change would alter a file format, a public pyazr name, or a CLI flag in a
  way that is not backwards compatible.
- A change would touch `external/`, CI release publishing, or `main`/`qt5`.
