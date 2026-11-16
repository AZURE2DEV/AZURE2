# AZURE2 — guidance for Claude Code

AZURE2 is a multi-level, multi-channel R-matrix code for nuclear reactions and
scattering. Physicists use its numbers in publications and evaluations, so a
wrong number that looks right is the worst possible bug. Optimise for
**correctness you can show**, then for keeping the CLI and the Python API in
step, then for speed.

Human-readable versions of this file: `docs/source/developer/contributing.rst`
and `docs/source/developer/architecture.rst`. Directory-specific rules are in
`tests/CLAUDE.md` and `pyazr/CLAUDE.md`.

## Upgrade (branch `feature/upgrade`)

The engine computes one chain: level parameters → collision matrix → observables → apparatus fold → χ². The
upgrade rebuilds the code so that each stage of that chain is one module that depends only on the stages before it,
and so that every change can be checked against synthetic projects with known parameters and against the present
engine's results. Why, and in what order: `docs/upgrade/README.md`.

This branch carries the staged rebuild described in `docs/upgrade/` (start with `docs/upgrade/FOR_THE_OWNER.md`:
its section 11 is the reading order for an agent, its section 10 what is not finished). On it,
and on the phase branches cut from it (`feature/upgrade-<phase>-<topic>`), these rules apply on top of everything
below (full text and reasons in `docs/upgrade/RULES.md`):

- Code is shaped like the physics: `physics/`, `io/`, `stopping/`, `experiment/`, `inference/`, `convert/`, then
  the engine (`src/`, `include/`) and the binding (`api/`), plus `support/` for plumbing. (`inference/` holds the
  parameter vector as plain data, the priors and the rule of the objective's term list; the terms themselves read
  the engine's objects and are `src/Objective.cpp`: `docs/upgrade/inference.md`.) Modules use only the modules below
  them; each is its own CMake target, so a wrong dependency does not compile.
- Design before code. Each piece of physics (a channel type, a basis or convention, an observable, a kernel, a
  prior) is a self-contained part with a stated interface, its inputs and outputs in physical units, and one test
  that would fail if it were wrong; parts fit together only through those interfaces (`ARCHITECTURE.md` §9), never by
  reaching into each other's state. A new convention or feature starts as a short design note in `docs/upgrade/`
  (what it is physically, where it plugs in, what it replaces, how it is tested), reviewed before the code. No
  special cases threaded through unrelated code, no flags that change meaning elsewhere, no copy of a formula.
- The engine computes three lines (`docs/upgrade/ARCHITECTURE.md` §1). New physics enters at one extension point,
  changes one symbol, and appears in one of the test projects T1–T3 (`tests/upgrade/`).
- Setup state is read-only once a fit stage starts; parallel code writes only its own slice or scratch; no caches
  beyond the per-stage tables; results independent of thread count.
- One implementation per formula and format, in C++ (Coulomb functions: `physics/coulomb`, taken from pycoulomb). Library readers live in pyazr
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
- **Exception, until `feature/upgrade` and `dev` are joined:** work on the
  upgrade is cut from `feature/upgrade` and merged back into it, never cut from
  or merged with `dev`. A plain merge of `dev` conflicts in 93 files and merges
  others silently wrong (`docs/upgrade/FOR_THE_OWNER.md` section 7); a feature
  of `dev` is re-written in its named places (`docs/upgrade/dev_features.md`).
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
| `src/`, `include/` | C++ engine | `AZURE2.cpp` (main, CLI flags) → `AZUREMain` → `CompoundSystem` (the compound nuclei as set up: pairs, J^π groups, channels, levels' set-up data; made with the data by the one set-up `SetUpCompoundSystem`, `include/CompoundSystemSetUp.h`, and shared constant, `include/CompoundSystem.h`; the values at a vector are each owner's `LevelState`) + `DataSet` (the data as set up, immutable: segments, measurements, target effects, covariance; `include/DataSet.h`) held by each `WorkingCopy` (what an evaluation changes: norms, shifts, frames, the data points, `Unit`, and the slots of their target-effect folds, each a `Site` (`include/Site.h`: energies, angle rows, channel rows, Coulomb amplitude, screening, stopping, EC amplitudes; the observables read a `Site`), shared by the points of a sharing group on the shared mesh (`include/SlotGroup.h`: a slot is written by its first reader, `WorkingCopy::WritesSlot`, read by all; never write a fold's intermediate values into a slot); `include/WorkingCopy.h`) → `AZURECalc` (χ²), `AZUREGrad` (gradients), `AZURECalcMCMC` |
| `api/` | `AZUREAPI` + pybind11 module `_azure2` | Builds `pyazr/_azure2*.so` in place |
| `io/` | The `.azr` document and its format conversion | Standalone (no engine headers); the engine and pyazr read the file through it |
| `stopping/` | Stopping and straggling: materials, SRIM-2013's electronic and nuclear stopping (ported from LunaSRIM; SRIM's data files in `stopping/data/srim/`, cite SRIM), energy loss, range, TRIM/Bohr straggling with Tschalär propagation, with uncertainty and provenance | Uses `physics/` only; the engine calls it through `TargetStopping` (`stopping evaluation`, `stopping table`), pyazr through `pyazr.stopping`; `docs/upgrade/stopping.md` |
| `experiment/` | Model expressions (what a data set measures), lab → c.m. frames, the table of project format 3 | Uses `physics/` and `io/`'s `Record.h` only (no engine headers); `experiment/README.md` |
| `inference/` | The fit's parameter vector as plain data (kinds and records, the level block's layout, the index map full ↔ packed, the posterior's support), one Gaussian prior, the prior covariance (grammar and arithmetic), the sums of a list of terms | Uses nothing (no engine headers, no Minuit2: check 8b, `inference_standalone`); the engine names its types through `include/InferenceTypes.h` and keeps the walks, the values of the vector (`AZUREParams`) and the terms (`src/Objective.cpp`); `docs/upgrade/inference.md` |
| `physics/transfer` (in `physics/`) | Optical potentials, bound states, distorted waves, DWBA transfer, the THM vertex in the plane-wave limit | A leaf: target `azure2_physics_transfer`, used by the engine's THM distortion (`src/ThmDistortion.cpp`) and the binding; checked against FRESCO (`docs/upgrade/transfer.md`) |
| `pyazr/` | Python package | In-process engine; reads project format 3 only: `Project` (`project.py`, `nuclear_part.py`) edits the engine's document with the engine's grammars; old files are refused with the conversion command (`convert.py`: the converter) |
| `numcmc/` | Ensemble MCMC sampler | |
| `external/` | Vendored Minuit2 (Coulomb functions: `physics/coulomb/`, from pycoulomb; stopping: `stopping/`) | **Do not edit**; it is upstream code |
| `tests/` | All tests | See `tests/CLAUDE.md` |
| `docs/source/` | Sphinx site, built with `-W` | The user-facing source of truth |
| `.claude/skills/` | How to *use* AZURE2 | See "Skills" below; `tools/check_skills.py` |

## Build and verify

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"   # without pybind11 the pyazr module is silently skipped
cmake --build build --parallel
./scripts/check.sh            # everything below, in order; --quick skips the physics suite
JOBS=2 ./scripts/check.sh --quick   # JOBS: compile jobs and parallel ctest tests (default 4; about 2 threads per test)
```

What `scripts/check.sh` runs, and what to run by hand when iterating:

```bash
ctest --test-dir build -LE slow --output-on-failure   # unit, reference, pyazr, consistency (minutes)
./tests/run_tests.sh build/src/AZURE2                 # physics regression: χ² of every tests/<name>/ (tens of minutes)
ctest --test-dir build -R <name> --output-on-failure  # a single test
ctest --test-dir build -L slow                        # long pyazr tests: run when you change pyazr fitting/bands
sphinx-build -W --keep-going -b html docs/source docs/_build/html   # docs, warnings are errors (as CI)
ctest --test-dir build -R docs_build                  # the same build as a test (skipped without Sphinx)
```

`tests/run_tests.sh` without an argument takes the *newest* `build*/src/AZURE2`.
Pass the binary explicitly if more than one build directory exists. The
Coulomb table (`--build-coulomb-table`, about 5 min, 88 MB in `~/.azure2/`)
makes runs about 4× faster and does not change results.

## Definition of done

A change is finished when all of these are true, and the commit message says
how each was checked:

1. **It builds** with pyazr enabled, and adds no compiler warnings in the files
   it touches. `feature/upgrade` has no Qt GUI (decision 3 in
   `docs/upgrade/HANDOFF.md`); the web UI of Phase 9 replaces it, checked against
   `docs/upgrade/qt_parity.md`.
2. **The tests pass**: `ctest -LE slow` and `tests/run_tests.sh`. If a
   reference χ² moved, see the next section first.
3. **It is tested.** A bug fix comes with a test that fails on the old build
   and passes on the new one, and you say so. A new physics feature is checked
   against something independent of this code: a closed form, an
   `mpmath`/m-scheme/brute-force `truth.py`, a published number, or an
   identity such as a constant cross section surviving a fold. Agreeing with
   itself is not a check.
4. **Every front end agrees.** A user-visible feature reaches the engine, the
   CLI, pyazr and the docs in the same branch (see the checklist below). If one
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
  would contain, bump its signature or version (e.g. `WorkingCopy::ECSignature`) so
  old files are recomputed, not silently reused. Add a stale-file case to the
  test.
- **Old inputs are converted, not read.** The engine reads project format 3
  only and refuses anything else with one message (`io::NotFormat3Message`,
  exit 3); an older `.azr` goes through the converter (`AZURE2 --convert`,
  `convert/`), a `param.sav`/`param.par` of an old project through
  `AZURE2 --convert-params`, and the converter keeps the old meaning or refuses
  with the old engine's reason. Data files of earlier versions must still load.
  A new field is optional and its default reproduces the old result; a change
  to format 3 that changes an existing project's meaning bumps
  `Schema::kVersion` and gets the converter's support and a test with an old
  fixture (`tests/convert/fixtures/`).
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
- **No hidden state.** Constant tables (masses, constants, the schema) are built
  once and reached through a function returning a `const` reference; nothing
  else is global. What a calculation writes belongs to the caller and is passed
  in (`Results`, per-stage tables, scratch); a point or a data set holds where
  its record lies, never a buffer [`repo_consistency_test.py` check 25]. Do not add a function-level
  `static` or `thread_local` that holds a computed value, or a cache keyed by
  the last call. State lives in one of three owned places: the energy task's
  scratch (`include/TaskScratch.h`, in the caller's `Results`), the run's
  per-thread scratch (`include/RunScratch.h`, `Config::scratch`), the run's
  once-only warnings (`Config::FirstTime`). What a fitter's or a sampler's
  `const` evaluations keep between calls is in the caller's `FitterStates` or
  `SamplerStates`. The allow-list holds permanent tables only
  [check 38, `tests/consistency/allowed_state.tsv`; recipe "Add a cache"].
- **Extension is a table row.** A new quantity, χ² term, kernel or stopping
  model is one row in the one table of its kind, with its test. No class that
  registers itself from a static object, no factory hierarchy.

## Recipes: where a new feature goes

The engine reads the `.azr` (format 3 only) once, through `io/` (`io::Document`, `Config::SectionStream`) and
`experiment/`'s table (`Config::project`: the records validated and built into specs), and pyazr reads it
through the same code in the binding; the χ² is assembled in one place, `src/Objective.cpp` (its terms; the rule of
their sums is `inference/`'s `TermList.h`). Each recipe below is
the whole list of places to touch; `tests/consistency/repo_consistency_test.py` fails when the rule in brackets is
bypassed. Run it first (under a second), then the build and the tests of "Definition of done".

**A data-set or extrapolation keyword**: a field of the record, recipe "A format-3 field" below (and, when an old
segment keyword means the same, its conversion in `convert/`, `convert_constructs`). No engine code reads an old
section (check 34).

**A `<config>` key or a new `.azr` section** (a section of records: the recipe "A format-3 field" and its kind in
`Schema.cpp`; the one place is the table).
1. Read it in its owner from `Config::SectionStream("<name>")` [check 12: the project file is opened only in
   `src/Config.cpp`]. A section's grammar belongs to one parser; pyazr calls it through `api/src/` (as
   `pybind_nuclear_part.cpp` does for `<config>` and `<thm>`), never a Python copy.
2. Its absence must reproduce the old result. If an old file states the same thing another way, the converter
   writes it (`convert/`, a case in `convert/tests/constructs_test.cpp`); `AZURE2 --convert` and
   `pyazr.convert_azr` then cover it.
3. pyazr: `pyazr.Project` edits it through a binding of the engine's grammar (`api/src/`; Python hands over
   names and numbers, the engine makes the line); an unedited project writes the file back byte for byte
   (`tests/pyazr/project_test.py`, `azr_layout_test.py`).
4. Docs: `reference/data_formats.rst` (records: `reference/format3.rst`, generated).

**A term of the objective** (a prior, a penalty, a constraint).
1. One `ObjectiveTerm` subclass in `src/Objective.cpp`, appended to the term list in the `Objective` constructor
   (`docs/upgrade/objective.md` §4; `ObjectiveTerm` is `azure2::inference::Term<ObjectiveContext>`, whose defaults
   say how the sampler takes a term); its gradient beside the others in `src/AZUREGrad.cpp`, and its residual rows in
   `AppendModelRows` if the Levenberg-Marquardt path needs them. A Gaussian prior on one parameter is not a new
   formula: `GaussianPriorChi2` / `GaussianPriorGradient` / `GaussianPriorCurvature`
   (`inference/include/azure2/inference/GaussianPrior.h`) and a line in `ForEachDiagonalPrior` (`include/Objective.h`).
   Arithmetic that needs no engine object goes in `inference/` with its test in `inference/tests/` [check 8b: only
   its own headers and the standard library]. Every path -- CLI fit, gradient, API, MCMC, writer -- asks the objective
   [check 11: no term is added outside `src/Objective.cpp`].
2. The analytic gradient agrees with finite differences (`AZURE_GRAD_CHECK`, a pyazr test); the CLI and pyazr give
   the same χ² (`tests/pyazr/objective_test.py`).
3. Its total goes into `Results::objective` (the evaluation's buffer) and from there to `chiSquared.out`; docs
   `output_files.rst`.

**A physics function or an observable.**
1. In `physics/` (no engine headers, no GSL: check 8), with its unit test in `physics/tests/` against an
   independent value (closed form, mpmath `truth.py`, a published number).
2. Exposed in the binding with array inputs (`api/src/pybind_physics.cpp`) so pyazr and the web UI call it.
3. Constants come from `physics/include/azure2/physics/Constants.h` [check 10: no copy with another value], masses
   from `pyazr/data/mass_1.mas20` (`pyazr.sources.masses`; the engine's THM nuclei are generated from it and checked by
   `pyazr_thm_nuclide_masses`).

**An observable of a reaction at a mesh energy** (what a site's record holds; `include/SiteObservable.h`).
1. A kind in `SiteObservable` and its row in `kSiteObservables` (`src/SiteObservables.cpp`: the record's slots and the
   contraction with a reader's row, `Site::ValueFor`); `ObservableOf` names the terms of that kind.
2. Its formula at E_n, one branch of `GenMatrixFunc::CalculateBareCrossSection` (the sum itself a function of
   `physics/` `Observables.h`); its adjoint a function of `ObservablesAdjoint.h` (T-matrix cotangents from the
   observable's; its test against finite differences of the forward function in
   `physics/tests/observables_adjoint_test.cpp`), called from the branch of `AMatrixFunc::PointTBar` that walks the
   same pathways (checked with `AZURE_GRAD_CHECK`).  A kind whose formula reads the level matrix's elements
   instead of a T-matrix (kHoes, THM's σ_HOES: `HoesCrossSection`) is a branch of `SolveAtOneEnergy` and of the
   gradient's task (`GradEnergyTasks.cpp`, `HoesTask`), asked by the kind.
3. Anything a reader needs per slot (angular row, Coulomb amplitude, A_y factors) in its row set (`include/RowSet.h`,
   built by `Unit::RowSetOver`; shared by the readers of one angle as the data give it, `WorkingCopy::BuildRowSets`:
   for a lab angle the c.m. angle and dΩ_lab/dΩ_c.m. of each slot, `experiment::LabAngleAt`, energy_grid.md §4.3),
   held per slot or once by the set's one rule (`RowSet::perSlot`, read through its accessors), checked by
   `WorkingCopy::CheckReaderRows` (`AZURE_NODE_CHECK`).  A kind with two values per slot (`values = 2`: the weight w and
   the element r, as A_y and its cross section) is folded as (K·(w r))/(K·w) with no further code (`Unit::FoldsRatio`,
   `RatioPartsAt`; the adjoint in `AZUREGrad.cpp`); only where the weight's cotangent goes is named
   (`Unit::RatioWeightRoute`).

**A new optical potential term** (`physics/transfer`, `docs/upgrade/transfer.md` §2).
1. A field in `OpticalPotential` (`physics/include/azure2/physics/transfer/OpticalPotential.h`, with its formula, unit
   and sign in the header's Eq. (2.1)), its `operator==`, its value in `TermsAt` (one member of `PotentialTerms` per
   FRESCO TYPE, so that the listing stays comparable with FRESCO's fort.34) and its line in `Describe`.
2. If it depends on l or j, its factor where `SolveDistortedWaves` and `SolveBoundState` build U_lj.
3. A FRESCO case in `physics/tests/data/make_fresco_benchmarks.py` that uses it (run the script, commit the data) and
   its check in `physics/tests/transfer_fresco_test.cpp`: the listing first, then S and the cross sections.

**Add a global optical set** (`opticalAA=<name>` / `opticalSF=<name>` of a THM experiment; `docs/upgrade/transfer.md`
§24.4). Never a second copy of the numbers in `src/`.
1. One row of `kSets` in `physics/src/transfer/GlobalOptical.cpp`: name, projectiles in words and as (Z, A),
   reference, target-mass and lab-energy range, and the function of (Z_p, A_p, Z_t, A_t, E_lab) that returns the
   `OpticalPotential` (R = r0 A_t^{1/3}; the helper `Central` for central sets). Spin-orbit only with its primary
   source or an independent code's table in front of you (name the file in the comment) and a test.
2. Its points in `physics/tests/transfer_global_optical_test.cpp` (at least two per row; the expected ten numbers from
   the paper's table or an independent code, `tests/reference/thm_optical_reference.py` for RIPL-3 sets: say which).
3. Its line in the table of `docs/source/theory/thm_implementation.rst` ("Global optical potentials"). The engine,
   the parser and pyazr's `thm_global_opticals()` read the table; nothing else changes.

**A model-only datum** (no beam energy: a resonance strength, a width, an ANC; `docs/upgrade/data_model.md` §4, M10).
The objective gets no new term: T9, `ModelDataTerm`, evaluates every row of `<modelData>`.
1. Its formula in `physics/` (`ResonanceQuantities.h` is the pattern) with a unit test against numbers by hand and its
   derivatives against central differences; the binding (`api/src/pybind_physics.cpp`).
2. Its row in the table (`experiment/src/Quantities.cpp`, `q.recordKinds = kDatum`; the recipe "A format-3 quantity")
   and in `ModelQuantities()` (`experiment/src/ModelDatum.cpp`), with what `BuildModelData` must find in `<levels>`
   (level, pairs, channel) and its refusals (a code each, recipe "A format-3 diagnostic code").
3. Its case in `LevelValue` (`src/ModelData.cpp`), from the table `CompoundSystem::ObservedOfFitValues` returns (never
   formal parameters); its gradient comes with it (central differences over the level's J^π group).
4. A function `test_<quantity>` in `tests/model_data/check.sh` against a value by hand [check 27: every datum
   quantity of the table has its `ModelQuantities()` row and a test named after it]; `reference/data_formats.rst`.

**A fitter** (a minimiser the program and a pyazr session both run; `docs/upgrade/fit.md`).
1. Its function `FitResult FitX(AZURECalc &, AZUREParams &, const FitControl &)` in `src/FitX.cpp`: start, limits and
   fixed flags from the `AZUREParams`, its end point left there; lines to `control.log`; steps to `control.progress`
   (an evaluation fitter keeps a `FitEvaluations` alive while it runs: count, best point, the stop request in force;
   an iteration fitter sends a `FitProgress` per iteration and checks `control.stop`); a stopped fit returns the best
   point seen with status `stopped`.
2. Its row in `Fitters()` (`src/Fitters.cpp`; the row's place is its precedence in the program); the switch in
   `src/AZURE2.cpp`, `reference/command_line.rst` and the table "The Engine's Fitters" of
   `docs/source/user_guide/fitting.rst` [check 46].
3. Its case in `tests/pyazr/fit_with_engine_test.py` (the session's fit ends where the program's ends).
   `fit_with_engine(method=)` and `pyazr.kinds()["fitters"]` read the table.

**A CLI flag.** `parseOptions`, `printHelp` and `docs/source/reference/command_line.rst` [check 1].

**An option of the calculation** (a switch a project's result depends on: a formalism, an approximation, a numerical
method; `docs/upgrade/project_format3.md` §3.12). It is a `<config>` key, a flag and a pyazr keyword under one name.
1. One `Option(name, flag, negation, default, effect, doc)` row of `io::ConfigKeys()` (`io/src/ConfigSection.cpp`);
   the program's parser and `--help` (`ApplyOptionFlag`, `printOptionFlags`), `--describe`, `set_config`,
   `pyazr.kinds()["options"]` and the web UI read it there [check 1; `io_config_section`]. A flag that sets another
   option too, or a pair that cannot hold together, is a row of `ConfigOptionImplications()` /
   `ConfigOptionContradictions()`.
2. Its place in `Config` (a bit of `paramMask` or a bool member): one row of `OptionStorages()` in
   `src/ProjectOptions.cpp` [`upgrade_project_options`: a row without one, or a default unlike `Config::Reset`, fails].
3. Its keyword in `pyazr.Session.__init__`'s signature and docstring (with the row's default) [`pyazr_describe`].
4. Its case, both values, in `tests/config_options/check.sh` (key = flag) and `tests/pyazr/config_options_test.py`
   (key = keyword); its line in the table of `docs/source/reference/command_line.rst`.
An option that is not project state (a file or directory of this machine, the run's command) is a field of
`RuntimeOptions` with its row and reason in `NotProjectOptionRows()` (`api/src/pybind_describe.cpp`).

**A new regression check** (`tests/<name>/check.sh`; `tests/CLAUDE.md` for the house pattern). A variant of a
project is made on its format-3 records, and every edit must be shown to bite: an edit that matches nothing leaves
the variant equal to the base, and the check then compares a run with itself and passes without testing anything.
1. Edit through the guard of `tests/lib/edit_guard.sh`: `edited SRC DEST awk '...'` (or `sed`, `perl`, a shell
   function) fails when DEST would equal SRC; `must_differ A B` for a variant made another way (a heredoc, pyazr).
   THM records through `tests/lib/thm_records.sh` (it fails on a data set that is not there).
2. Where the check's logic needs the variant's result to differ from the base, assert that it does (a number, an
   output file), not only that the run ended.
3. Inject the error the check protects against once (a one-line engine change in a scratch build), see it fail, and
   say so in the commit message.
4. `tools/edit_audit.sh <AZURE2> <dir> <name>` then `python3 tools/edit_audit_report.py <dir> --details` lists every
   edit that hit nothing and every pair of runs with different inputs and equal outputs
   (`docs/upgrade/format3_tree.md` section 7).

**A fit parameter of an apparatus effect** (a number of a target, beam or detector record that a fit varies;
`docs/upgrade/apparatus_parameters.md`). A table kernel's own parameter needs no engine code: it is kind 12,
`apparatus`, through the row `kernel` of `ApparatusRows()`.
1. In the kernel's row, the parameter's `KernelParameter` says `fittable = true` (its name, dimension, engine unit); the
   kernel reads its value from `KernelParameters::values[k]`, k its place in the row; the fold's grid is laid out at the
   parameter's layout value through `KernelParameters::layoutValues` (`EffectSupport`; nothing in the kernel).
2. The format-3 form is generated with the nuisance grammar `<value> <unit> fixed|fit [prior <c> <w>[%]] [grid_max
   <v>]` (`Schema.cpp`); `BuildEffects.cpp` makes an `ApparatusParameterSpec` when the value is `fit` or has a prior
   or grid_max: name `<record>.<field>.<parameter>`, shared by every data set naming the record.
3. Tests: the kernel's `derivative` against central differences (`kernels_test.cpp`); a fit of pseudo-data made at a
   known value, with a prior, shared by two data sets, and the refusals (`tests/rectangular_resolution/check.sh` (c)-(f)
   is the pattern); the schema (`AZURE2 --schema`, `tools/gen_schema_docs.py`).
4. Docs: the kernel's row in `energy_grid.md` section 3.3 and its user-guide paragraph.
A number that reaches the calculation some other way than a kernel's values (as the stopping scale re-sets the
sub-points' stopping) is a new *row* of `ApparatusRows()`: one function file with its declare / write / read / follow
functions (`src/StoppingScale.cpp`, `StoppingScaleRow`) and its line in the table; nothing else in the engine.
Derivative: central differences of the whole chi^2 (`AZURECalc.cpp`, as a density); an analytic one is not wired.

**A derived parametrisation of a level** (numbers a user writes and fits in place of some of a level's partial widths:
Γ_γ and branchings, `gamma_decay`, `docs/upgrade/gamma_parametrisation.md`; a spectroscopic factor or an ANC in place
of a particle width would be the next). The engine's physics is not edited: the numbers sit in the vector at the
positions of the amplitudes they replace, one number per amplitude, and a map turns them into amplitudes at every fill.
1. The record: a kind in `experiment/src/Schema.cpp` (its section in `Schema::Sections()`), its build against the
   `<levels>` tables in a file of its own (`experiment/src/GammaDecay.cpp`), its refusals as codes (recipe "A format-3
   diagnostic code"); `Schema::Build` calls the build for `--validate` and pyazr.
2. The map in one engine file (`src/GammaDecayMap.cpp`): reading the section (as `ModelData::Read`), the input values
   at set-up (so that TransformIn converts them as any width), `AddParameter` (called in
   `CompoundSystem::FillMnParams` for the slots it owns), `ToEngine` (the first line of
   `CompoundSystem::FillCompoundFromParams`) through the engine's own conversion (`physics/`
   `amplitudes_from_observed`), `PullBack` (called by `GradAccum::Scatter`: the chain rule onto the numbers and onto
   the level's other parameters), the priors as rows of an objective term, the report in `parameters.out`.
   A second map is a second member of `CompoundSystem` beside `gammaDecays_` and one more line in each of those four
   places.
3. A parameter kind (recipe "A new kind of fit parameter"; `IsLevelSlot` covers a number in a level's slot).
4. Tests: the same project written both ways gives identical calculations (`tests/gamma_decay/check.sh` is the
   pattern), fits in both parametrisations with equal minima and propagated errors (`tests/pyazr/gamma_decay_test.py`),
   `AZURE_GRAD_CHECK` on the new numbers and on the level's energy and amplitudes, each check shown to fail once.

**A new kind of fit parameter** (a nuisance that is not a number of an apparatus record: of a segment, of the level
scheme, of a THM experiment; `docs/upgrade/parameter_kinds.md`). An apparatus record's number is the recipe above,
not a kind. What a parameter is lives in its record, never in its name [check 23: no prefix or suffix test of a
parameter name].
1. `inference/include/azure2/inference/ParameterKind.h`: the enum value (the next code) and its `ParamKindTable()` row
   (key, unit, owner fields,
   `ownPrior` / `segmentField` / `varyFromProject`, doc); a `ParameterRecord` factory for new owner fields.
2. Create it with `AZUREParams::Add(name, value, step, ParameterRecord::...)` in `CompoundSystem::FillMnParams` or
   `WorkingCopy::FillMnParams`, the name built there by one function; every later question asks `params.Record(i)`.
3. Follow the layout in `BuildParamIndexMap` (`src/AZUREGrad.cpp`: one `Append...` of `ParamIndexMap`,
   `inference/`'s `ParameterIndexMap.h`, at the kind's place in the vector, with its lookup and its case in
   `inference/tests/parameter_index_map_test.cpp`; a number of a level: a slot of `LevelParameterLayout`, which every
   walk of the level block reads; Developer Guide, "A level parameter"); the fields that say what it belongs to are
   rows of the record table `ParameterFields()` (`api/src/pybind_records.cpp`: name, unit, the kinds it describes,
   owner flag, doc, reader), which pyazr's `Parameter` and `Parameters.table()` read with no Python edit
   (`tests/pyazr/parameter_records_test.py` fails for a kind no owner field describes).
4. Docs: the kinds table of `reference/output_files.rst` and `parameter_kinds.md` section 1 [check 23]; a project
   with the kind in `tests/pyazr/parameter_kinds_test.py` (its names rebuilt from the records).

**A format-3 field** (a key of a record kind: `docs/upgrade/project_format3.md`, `experiment/README.md` §2).
1. One `FieldDef` row in its kind's function in `experiment/src/Schema.cpp` (`DatasetKind()`, `TargetKind()`, ...),
   made by the helper that says where it applies and whether it is required: `Req`, `Opt`, `RequiredWhen`,
   `OptionalWhen`, `ByQ` (the quantity decides), wrapped in `Repeatable`/`Defaulted`/`RequiredIf` where needed, e.g.
   `RequiredWhen("charge_state", "Charge state", IntAtLeast(1), Condition::Dimension({"per_charge"}), "the beam's
   charge state q")`. A structured value is a `FormSpec` of `StructuredValues()` built from `Keyword(...)`/`Unkeyed(...)`
   alternatives and `Pos`/`Kw`/`OptionalKw` parts. If the field depends on the quantity, `ByQ(...)` there and an
   `Allows(...)`/`Requires(...)` in each quantity function of `experiment/src/Quantities.cpp` that uses it.
   No parser code: the value type's parser exists (`Parsers.h`); a new value type is a parser plus a case in `ParseValue`.
   A positional part cannot be optional: an optional part is a keyword part (`OptionalKw("unit", ...)`, written
   `reference <file> unit mb`).  A `ByQ` field needs a value in the `plausible` map of
   `experiment/tests/quantities_test.cpp` (the test tries it on every quantity that does not list it).
2. A fixture line in `experiment/tests/data/` that validates, and a case in `experiment/tests/diagnostics_test.cpp`
   if the field brings a refusal of its own.
3. `AZURE2 --schema > tests/schema/schema_v3.json`, bump `Schema::kVersion` if an existing project changes meaning
   [check 21: `experiment_schema`, `pyazr_schema`]; `python3 tools/gen_schema_docs.py` [check 18].

**A format-3 quantity.** One function in `experiment/src/Quantities.cpp` that fills its `Quantity` member by member,
and its entry in `Quantities()`; members it does not set keep the defaults of `Quantities.h`:
```cpp
Quantity SFactor() {
  Quantity q;
  q.name = "sfactor";
  q.title = "Astrophysical S factor";
  q.doc = "S(E) = E sigma(E) exp(2 pi eta) with the engine's eta and atomic masses (data_kinds.md section 4.1).";
  q.recordKinds = kDataKinds;
  q.engineUnit = "MeV*b";
  q.units = {AnyUnitOf("sfactor")};
  q.energyRoles = {"E_cm", "E_lab"};
  q.reaction = {R::kChargedEntrance, "charged entrance"};
  q.fields = {Allows("multipole", kGammaExit)};
  q.tree = "Scale_SF(E) o Sigma(r)";
  q.exfor = {",SIG,,SFC", "PAR,SIG,,SFC", "EP/PAR,SIG,,SFC"};
  q.refusals = {RefusedAtValidation("E_REACTION_NOT_ALLOWED", "neutral",
                                    "the S factor is defined for charged particles (Gamow factor)")};
  return q;
}
```
Then a record of it in the `variants` of `experiment/tests/quantities_test.cpp` (the test fails for a row without
one), the schema and docs as above; later its tree builder (step 3).

**A new kind of data set** (a quantity the engine can fit, step 3b: the engine runs format 3).
1. Its row in the table: the recipe "A format-3 quantity" [`experiment_quantities` fails without its variant; check 18
   when `format3.rst` is not regenerated; check 24 for a new diagnostic code].
2. Its tree builder: the case in `Schema::Build` (`experiment/src/Build.cpp`, `ReadObservable`, called by `BuildDataset`) that
   turns the record into a `DatasetSpec` (`SetupSpecs.h`: the observable, the reaction's pair keys, the windows,
   the nuisances with their `ParameterNames`); until it has one, `Build` refuses the quantity with `E_NOT_BUILT`, so
   a file never runs with a quantity the engine does not compute.  What the engine computes for the spec is the
   model pass (`include/ModelPass.h`; the expression tree, S10).  A spec member the old lines never set (as
   `ReactionSpec::coulombReference`) goes into the comparison of `convert/tests/spec_equivalence_test.cpp` too; where an
   old construct means the same thing (`absolute 1` is `ReactionSpec::legendre` Absolute), the old-line adaptor
   `convert/src/SpecsFromOldLines.cpp` sets it, or `convert_spec_equivalence` fails.  A convention of the data (a
   normalisation, an expansion) is one table of its own in `experiment/` (`LegendreExpansion.h`), read by Build and
   the engine alike.
   A refusal that depends on the target, beam or detector the data set names is made where the effects are built
   (`experiment/src/BuildEffects.cpp`, `CoulombRatioFolds` is the pattern).
3. If the spec needs a tree the existing kinds do not build: one row of `TreeKinds()` (`experiment/src/Prediction.cpp`,
   its `TreeKind` and builder; `TreeKindOf` decides) and, for a value the R-matrix does not compute, a `LeafPart`.
   Such a leaf's value at a point is computed in `PointInputs` (`src/ModelPass.cpp`) by a function of its own file
   (`src/CoulombReference.cpp`: the Rutherford / Mott leaf; `src/ReferenceTableLeaf.cpp`: a tabulated reference),
   folded on the point's own sub-points by `FoldReferenceAtPoint` (`CoulombReference.h`: the one quadrature and the
   edge blend, never a copy of them); add the part to `IsReferenceLeaf` (`Prediction.h`) [check 35], which `PointInputs`,
   `ComposeModelValue` (no point is calculated for it) and `ModelPointBackward` (it leaves the adjoint's leaves) all
   ask.  A table the record names is read once by Build and kept in the spec in engine units and frame
   (`ReferenceTable.h`: energies taken to the lab once); the engine's own build reads only a data file's first row
   (`CheckDataFile`), so a refusal that needs every row (a row outside a table) reads them for that quantity only.
   A factor on the value that depends on the point (E exp(2 pi eta + g E), a Legendre order's normalisation) is a
   `FactorRef` kind (`Prediction.h`; its value in `PointInputs`) under the tree kind `Scaled`, not a new tree.
   The coefficients of a Legendre extrapolation (`AngularDistribution`) are no tree value: the `.acoeff` writer
   (`src/DataOutputFiles.cpp`) and the API take them from `Segment::LegendreCoefficients`, as the sigma and S columns
   of `AZUREOut` take theirs from `SigmaAndS` there.  A data value that is the
   same in the lab and c.m. frames (a ratio) answers `Segment::IsDimensionlessRatio`, so it is not converted at load.
   Then `python3 tests/upgrade/predictions_test.py --binary <AZURE2> --write` (lines added for new projects only).
4. A test: a case in `experiment/tests/build_test.cpp` (the spec member by member, and its refusals), a case in
   `experiment/tests/prediction_test.cpp` for a new tree, and when an old
   construct converts to it, `convert_spec_equivalence` (the old line's spec equals the record's) and a run of the
   converted file (`tests/format3_engine/check.sh` is the pattern: outputs byte for byte the old file's).  A quantity
   no old construct has gets a test project of its own in format 3 with an independent answer
   (`tests/ratio_to_rutherford/`: `truth.py` with mpmath, `check.sh`, models under `models/` with their data in
   `models/data/`; io_convert leaves a format-3 `.azr` alone).
   `AZURE2 --describe` must show the observable and the pair keys it resolved to.

**A format-3 unit.** One `Unit` row in `experiment/src/Units.cpp`, `{token, dimension, factor to the dimension's
engine unit[, aliases]}`, e.g. `{"keV", "energy", 1e-3}` (a new dimension is also a row of `Dimensions()` there, its
base unit with factor 1, as `inverse_energy`, 1/MeV); `experiment_parsers` checks every row's dimension and factor and `experiment_quantities` tries the unit
on every quantity; then the schema and docs as above.  A factor is applied only by `Units.h` `Rescale` (through
`ToBase`, `ToUnit`, `Context::InBase`; the columns of a data file through `experiment/DataRows.h`, where the points
are made), never as a literal in `Build`.  A quantity that is a scale of an observable the engine computes (an S
factor, a yield per unit charge) is `ReactionSpec::scale`, one tree kind `Scaled`: a constant factor is set in
`experiment/src/BuildUnits.cpp` (units, the yield per charge), an energy-dependent one in `Build.cpp` (`SFactorOf`:
`sFactor` and the modified S factor's g); its per-point value is a `FactorRef` evaluated in `PointInputs`
(`src/ModelPass.cpp`; the conversion in `src/CoulombReference.cpp`), and `src/DataOutputFiles.cpp` (`SigmaAndS`)
writes the sigma and S columns of such a data set.

**A format-3 diagnostic code.** A constant in `experiment/include/azure2/experiment/Diagnostics.h` (`io/`'s reader:
`io/include/azure2/io/Record.h`), its row in `DiagnosticCodes()` and a case in `experiment/tests/diagnostics_test.cpp`
that triggers it alone [check 24].

**A new apparatus kernel** (a resolution or broadening of a folded point: `experiment/include/azure2/experiment/Kernels.h`,
`docs/upgrade/energy_grid.md` section 3.3). A kernel is plain functions of the point's nodes, the point and the effect's
parameters; the engine applies whatever the table holds.
1. One file `experiment/src/Kernel<Name>.cpp` (`KernelGaussian.cpp` is the pattern): the fold as a template on V
   (`double` gives the value, `LinearForm` the exact weights the gradient reads: the same arithmetic), its support (the
   c.m. interval the point's nodes must span; the set-up gives the point the shared mesh's nodes over it, whose two
   outer nodes lie up to one interval beyond it: the fold integrates over `nodes.Low()` .. `nodes.High()`, never over
   the end nodes' span) and `Kernel <Name>Kernel()` filling its
   row member by member. Its own parameters are `KernelParameter` rows (name, dimension, engine unit, fittable) whose
   values arrive in `KernelParameters::values`; with `tableForm = true` and the keyword `detector.energy_resolution
   <form>` the format-3 form is generated from those rows (`Schema.cpp`) and read by `BuildEffects.cpp`: no parser
   code. An analytic derivative in a parameter goes in `derivative`.
2. Its line in `Kernels()` (`experiment/src/Kernels.cpp`) and the file in `experiment/CMakeLists.txt`.
3. Tests: cases in `experiment/tests/kernels_test.cpp` against closed forms (the table's own checks run on every row),
   a derivative against central differences; a format-3 fixture that validates (`experiment/tests/data/`, listed in
   `schema_test.cpp`); `AZURE2 --schema > tests/schema/schema_v3.json`, `python3 tools/gen_schema_docs.py` [check 18,
   21]. A test project of its own (`tests/rectangular_resolution/` is the pattern: the folded run against an average of
   the unfolded one) also adds its lines to `tests/upgrade/predictions.txt` (`python3 tests/upgrade/predictions_test.py
   --binary <AZURE2> --write`, lines added for the new project only) [`upgrade_predictions`].
4. Docs: its row in `energy_grid.md` section 3.3, its form in `reference/data_formats.rst` (generated part) and a line in
   `developer/edata_and_target_effects.rst`.
The engine is not edited: the fold (`src/PointFold.cpp`), the grid's extent (`EffectSupport`, `src/TargetSetup.cpp`) and
the gradient's weights (`FoldWeights`) read the table. Limits (2026-11): a table kernel cannot be combined with a
target, a beam spread, a beam profile or a UDR on one data set (refused in `TargetEffect::SelectKernels`; kernels are
chained by the grid branch); a kernel that needs a new datum of the point or a non-scalar parameter adds it to
`KernelPoint`/`KernelParameters` and to `TargetEffect::KernelParametersFor`; a parameter with `fittable = true` is a fit parameter when the project says `fit` (recipe "A fit parameter of an
apparatus effect").

**A new material or a measured stopping table** (`stopping/`, `docs/upgrade/stopping.md`). A material is data, not
code: `MaterialFromFormula("Ta2O5")` or `MaterialFromFractions(...)`, then set its fields: `phase` (SRIM's default is
gas only for a single natural-gas element; a gaseous compound such as CH4 needs `Phase::kGas`), `compound_correction`
(SRIM's factor, 1 = Bragg), `mass_density_g_cm3`. Elements 1-92 and ions 1-92 work without any edit: SRIM's data files
cover them. A measured stopping, or a correction to SRIM's for one element or projectile, is a user table:
`Material::user_stopping` with `ReadStoppingTable` (a unit line `E[MeV] eps[eV/1e15atoms/cm2]`, `S[MeV*cm2/g]` or
`S[keV/um]`, a stated uncertainty and a provenance string) or a `formula`; it overrides the evaluation and
`StoppingProvenance` shows it. SRIM's files in `stopping/data/srim/` are never edited. After a change to the formulas,
regenerate `stopping/tests/data/srim_reference.csv` with `stopping/tools/make_srim_reference.py` (any Python with
`lunasrim`) only if LunaSRIM itself changed, and run `ctest -R stopping`. Every output that uses these stopping powers
cites SRIM (`kSrimCitation`). The library includes only `azure2/stopping/`, `azure2/physics/` and the standard library
[check 32].

**A new stopping source** (beside `formula`, `evaluation`, `table`; `docs/upgrade/stopping.md` §14). A target's stopping
is computed in one place, `TargetStopping::PerActiveNucleus` (`src/TargetStopping.cpp`), lab MeV cm² per active nucleus;
`TargetEffect::StoppingAt` multiplies the stopping scale once and every reader asks it.
1. The keyword: one `Keyword(...)` of `f.stopping` in `experiment/src/Schema.cpp` (its parts, `scale` among them) and
   its `StoppingSpec::Source` value (`experiment/.../SetupSpecs.h`), filled by `StoppingInto` and keyed for the shared
   scale by `ScaleOf` (`experiment/src/BuildEffects.cpp`).
2. Its case in `TargetStopping`: built in the constructor (a library object from the spec, or `Error()`), evaluated in
   `LibraryPerActiveNucleus`, refused in `Check` (codes `E_STOPPING_DATA`, `E_STOPPING_RANGE`), described in
   `Provenance` (a source that is a published evaluation cites it there).
3. A test with an answer from outside the engine (`tests/stopping_yield` is the pattern), the schema and docs as in
   "A format-3 field".

Reference commits: `git show --stat 062b1bfe` (segment keyword `shiftfngroup`), `3e13fb28` (`normgroup`).

**Add a cache** (a value recomputed often enough to matter, e.g. a table of penetrabilities per mesh energy).  First
measure that it buys something (user-space instructions, `perf stat -e instructions:u`, are load independent).
Keys known at set-up (mesh energies, channels, the l range) make it a table: build it once in the set-up, hold it
`const` in the problem (the stage's node table, the channel table) and pass it to the readers.  Keys that move with
the parameters make it scratch: a member of `azure2::TaskScratch` (`include/TaskScratch.h`) when energy tasks read
it, of `azure2::ThreadScratch` (`include/RunScratch.h`) when set-up or hand-over code without a task reads it.  Key
it by everything the value depends on (never an address alone), bound its size, and give it a test that a changed
input misses.  Never a `static`, `thread_local`, namespace-scope variable or `mutable` member: check 38 refuses them
outside `tests/consistency/allowed_state.tsv`, and a row is never added for new code.  A warning printed once is
`if (configure.FirstTime("key"))`.

**Add a diagnostic switch** (an environment variable that prints, checks or selects another route).  Four places:
(1) a value of `azure2::support::Switch` in `support/include/azure2/support/Environment.h` and its row in the table
of `support/src/Environment.cpp` (same order; name, type, value when unset, kind -- numbers, cost, prints, test --,
what it does, its page, and for a switch of kind *numbers* its option of the program, `--name` for a flag or
`--name=VALUE`, with the range of a whole number its reader checks; check 43 refuses a numbers row without an option
and an option on another kind); (2) its row in the list "Environment switches" of
`docs/source/reference/command_line.rst`, with the option in the column Option (checks 40 and 43 keep the two the
same); (3) its reader: a switch that prints or checks reads the record where it is used,
`azure2::support::ProcessEnvironment().Set(Switch::kX)` / `.On(...)` / `.Text(...)` / `.Number(...)` (a hot path
keeps the value in a `static const` local); a switch of kind *numbers* is a field of `Config::EnvironmentSettings`
filled in `SettingsOf` (`src/Config.cpp`) and read from the run's `Config` (deep code without a `Config`: from the
pair's `PairFunctionSettings`, set in `CompoundSystem::Fill`); (4) a counter it prints is a member of
`azure2::RunDiagnostics` (`include/RunDiagnostics.h`), reached through `configure.diagnostics`, never a `static`.  The option needs nothing
else: the program's parser, `--help` and `--describe` read it from the row, and the record holds the option's value
as it holds the variable's (`StartProcessEnvironment`, made by `main` first; the command line over the options file
over the variable).  Never `getenv` (check 39 refuses it outside `support/src/Environment.cpp`); never `setenv` to
change a switch inside a process (pass the choice as an argument or on a copy of the `Config`).  Its test sets the
variable before the process starts (a ctest `ENVIRONMENT`, `env=` of a subprocess) or passes the option: the record
is made once.  A pyazr option for a *numbers* switch would be one field of `RuntimeOptions` and one line of
`apply_options` (`api/src/pybind_api.cpp`; `pyazr3.md` section 4, item 9).

**Comparing the outputs of two engines** (a change that must not move any number; the full text is
`docs/upgrade/ARCHITECTURE.md` section 10). Run every project of a `tests/` tree with each engine, then compare:
```bash
export SRC_TESTS=$PWD/tests WORK=<scratch> JOBS=<runs at once, one thread each; default 1>
export AZURE2_COULOMB_TABLE=<the table file the reference used>  # required: without it every charged channel moves
tools/runall.sh <absolute path of the reference engine> tab-ref   # once per reference
tools/runall.sh <absolute path of the new engine> tab-new         # a fresh tag for every run (REPLACE_TAG=1 to reuse one)
tools/cls.sh tab-ref tab-new; grep '^## ' $WORK/cls-tab-new.tsv   # "## <basis>: N files, N identical" when nothing moved
```
`AZURE2 calculate`, three bases (Brune, `--no-brune`, `--use-park`); the outputs land in `$WORK/cmp/<basis>/<tag>/.res-<folder>/
<file>.azr/` (`output/`, `checks/`, `run.log`), the table in `$WORK/cls-<tag>.tsv`, one line per file that differs.
*Identical* (`tools/classify.py`): the same lines of tokens after commas become blanks and `=` a token of its own,
numbers compared by value (`1.0` equals `1`), progress lines of `run.log` dropped. The converted-tree gate of format 3
(`tools/format3_gate.sh`) runs the same `runall.sh`.

## Code structure

- Put a new feature in its **own class and file**, as `ShftFunc`,
  `CoulombTable`, `DataCovariance` and `ScatteringLength` are. Do not grow the
  large files further: every hand-written source file of `src/`, `include/`,
  `physics/`, `io/`, `experiment/`, `inference/`, `api/src/`, `pyazr/` and `tools/` above
  1,000 lines is listed with its line count in
  `tests/consistency/file_sizes.tsv` [check 26: a listed file may not grow, an
  unlisted one may not pass 1,000 lines; only the director adds an entry].
  After shortening a listed file, `python3
  tests/consistency/repo_consistency_test.py --tighten-sizes` lowers its record
  (it never raises one). Moving code out of the large files is welcome, but do
  it in its own commit with no behaviour change and the full test suite run,
  never mixed into a feature.
- Use **one definition** of each formula, format and convention, and call it from both
  the engine and the gradient. When pyazr needs it, expose it through
  `_azure2` rather than re-implementing it in Python (pyazr needs the compiled
  module to read a project); where a Python copy is unavoidable, test it
  against the engine.
- Format C++ with `.clang-format` (2-space indent, `ColumnLimit: 0`, so you
  break lines yourself). Match the surrounding comment density. Comments say
  *why* and cite the equation or paper (`Edmonds 7.1.7`, `PRC 104, 064612`).
- No path of one machine (a directory under /data0 or /home) in a tracked file: tools find each other beside
  themselves in `tools/`, scratch is an environment variable (`$WORK`), data files are relative to the project
  [check 30; `docs/upgrade/HANDOFF.md`, a log, is left out].
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

The rules are `docs/upgrade/RULES.md` §3; `tools/check_skills.py` enforces items 4, 5, 6, 8, 10 and 11 and checks that
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

The skill set: `azure2-evaluation` (the procedure of a whole evaluation, which names the owner of each step),
`azure2-setup`, `azure2-calculate`, `azure2-fit`, `azure2-check`, `azure2-uncertainty`, `azure2-plot`,
`azure2-targets`, `azure2-thm`, `azure2-resonance-search`, `azure2-exchange`, plus `rmfit-campaign` (the table in
`docs/upgrade/RULES.md` §3 says which task each owns). A recipe of more than five lines of Python is a script in the
skill's `scripts/`, run by its `scripts/selftest.py` on a project of `tests/`; ctest runs every self-test
(`tests/skills/CMakeLists.txt`, `skill_<name>`). Run `python3 tools/check_skills.py` before committing a skill
change, and the skill's self-test when a script or pyazr changed.

## When to stop and ask

- A reference value or test has to change, and you cannot point to an
  independent check that the new value is right.
- A change would alter results for existing users' projects by default.
- A change would alter a file format, a public pyazr name, or a CLI flag in a
  way that is not backwards compatible.
- A change would touch `external/`, CI release publishing, or `main`/`qt5`.
