# Evaluation without copies and without hidden state (Phase 6)

Status: survey and design (2026-11-12, branch `feature/upgrade-6-evaluation`, from the merged head 88d55f17).  The
state table (section 1) has a column "became" that is kept true as the work lands; section 7 records what is done.

The rule (`ARCHITECTURE.md` §9, `RULES.md` §2, the repository's `CLAUDE.md`; the PI, 2026-11-10): **the only
process-wide objects are immutable tables built before the first evaluation; every other state is owned by the caller
and passed down.**  Why it matters: the owner must be able to say, for any number the engine prints, which inputs it
depends on; an AI agent that drives the engine must be able to run two calculations in one process, or one calculation
at two thread counts, and get the same bytes without knowing in which order caches were filled.

How the inventory was made.  `tests/consistency/state_scan.py` (new; the scanner of check 38) reads every `.cpp` and
`.h` of `src/`, `include/`, `api/`, `physics/`, `experiment/`, `stopping/`, `convert/`, `io/` (not `external/`, not
the tests, not `physics/*/transfer` of another branch), blanks comments and literals, tracks namespace, class and
block scopes, and reports every `thread_local`, every non-`const` `static` in a block or a class, every non-`const`
variable at namespace scope (named or unnamed namespace, `static` or not; a pointer to const is not const) and every
`mutable` member.  Its list was checked against `nm -C` of the built objects (every data symbol in `.bss`/`.data` of
the engine's objects): the scanner misses nothing that `nm` shows except `const` objects with a dynamic initialiser
(tables, and the `static const bool on = getenv(...)` switches of section 2, which are listed there).  Members that
are caches but not `mutable` (written in non-`const` methods of objects an evaluation owns) were found by reading:
section 1.6.

Counts at the start (symbols; a `thread_local int a, b;` is two): **thread_local 26** (32 declarations in 13 files, as
`ROADMAP.md` counts them), **static 40**, **namespace 42**, **mutable 53**.  `NuclearPotentialManager` and
`NodeTableCache` named by the roadmap are already gone (85d34aa6: the hybrid potential is a value held by its pair;
energy_grid.md branch 8: `NodeTableCache.h` deleted).

## 1. The state table

Columns: **symbol** (file); **holds**; **written by / when**; **read by**; **thread safety today** (mutex, per thread,
single-threaded phase, atomic, none); **order** — can its content change an output number depending on the order of
calls or the thread count, and how that is known (*read*: by reading the code; *run*: by the A-B-A test of section
7.2 or another run named); **becomes** — deleted / table (immutable, built at set-up, held by the problem) / scratch
(owned by an energy task or the caller, passed by reference) / Results (a field of the caller's buffer) / notices (the
run's message collector, `RunNotices`) / diagnostics (a record read once, section 2) / allowed (stays, on the
allow-list with its reason).  "Done" in the last column means the change is in this branch.

### 1.1 `thread_local` (the roadmap's 32 lines in 13 files)

| Symbol (file) | Holds | Written / read | Safety | Order | Becomes |
|---|---|---|---|---|---|
| `memo` (`include/WhitFunc.h`, `ScaledMemo`) | scaled Whittaker W at (l, η_b, ρ), 2¹⁷ entries | EC integrals' final-state factor, set-up and shifted EC | per thread | no: exact key of a pure function (read) | scratch of the EC integrals (section 7) |
| `t_radialBlocks`, `t_radialSize` (`src/CoulFunc.cpp`) | Coulomb waves per (pair data, l, 2J, E) and radius, 2²⁰ entries | `CoulFunc::operator()` under `SetRadialScratch` (EC integrals, EC phases) | per thread | no: key holds every input (read) | scratch of the EC integrals |
| `g_fgMemo` (`src/ECIntegral.cpp`) | F, G, W at each radius of one FW/GW pair | `ECIntegral::Integrate`, cleared per call | per thread | no (read) | a member of the integral's parameters (scratch of one call) |
| `holder` (`src/ECIntegral.cpp`, `ec_workspace`) | a GSL qag workspace | every EC integral | per thread | no: no value survives a call (read) | scratch of the EC integrals |
| `cache` (`src/ECIntegral.cpp`, `getRateGrid`) | Gauss–Laguerre nodes and final W per final channel (rate mode) | EC integrals of a rate | per thread | no (read) | scratch of the EC integrals |
| `cache` ×3 (`src/ChannelSpace.cpp`: `Cached`, `CurrentLevelArrays`, `CurrentLevelDerivatives`) | per J group: the formal-basis solver, the level arrays E, γ, S, 1/(E_i−E_j), the level-matrix derivatives | the channel-space route of every level-matrix solve and gradient task | per thread | no: keyed by address and `ALevel::Stamp` (a process-wide counter), so a reused address misses (read; A-B-A run) | scratch of the energy task (section 7) |
| `lastPath`, `last`, `lastSet`, `lastGeneration` (`src/CoulombTable.cpp`, `Shared`) | the last table looked up by this thread | every `CoulFunc` construction | per thread | no (read) | the table is resolved once at set-up and held in `Config`; deleted |
| `amatrix`, `thm` (`src/EnergyPass.cpp`, `SolveAtOneEnergy`) | an `AMatrixFunc` and a `THMMatrixFunc` with their buffers, re-pointed per task (`Reset`) | the point pass's energy tasks | per thread | no: `Reset` re-points them (read) | scratch of the pass, owned by the caller's `Results` |
| `units` (`src/EnergyPass.cpp`) | the task's unit pointers | the point pass | per thread | no | scratch of the pass |
| `amf`, `thm` (`src/GradEnergyTasks.cpp`) | as above for the gradient's tasks | gradient pass | per thread | no | scratch of the pass |
| `reusable` (`src/AZUREGrad.cpp`) | a `THMMatrixFunc` for `ComputeTHMRows` | THM rows | per thread | no | scratch of the caller |
| `memo`, `filled`, `next` (`src/ThmFunc.cpp`, `ThmCoulombTerm`) | 64 last THM Coulomb terms, keyed by pair address and values | channel tables of THM points | per thread | no: key holds the values (read) | deleted (measured, section 3.3) |
| `memo`, `filled`, `next` (`src/ThmVertexBoundary.cpp`) | 32 last S_c(E_λ) of the vertex boundary | `ThmVertexBoundary` set-up | per thread | no (read; a session leak by address was fixed earlier) | deleted (measured) |
| `memo`, `next`, `key` (`src/ThmLineshape.cpp`, `ThmLevelWidth`) | 64 last level widths for the THM line shape | THM points | per thread | no (read) | deleted (measured) |
| `table` ×2 (`src/ThmAngular.cpp`) | the geometric Legendre coefficients of a list of waves | THM angular windows | per thread | no | a value built by the caller per point (`ThmAngularCouplingOf`) |
| `entry` ×2 (`src/AZURECalc.cpp`, Windows only) | the `NAME=value` string given to `_putenv` | `AZURE_GRAD_CHECK` exit handling | per thread | no | deleted (`_putenv` copies its argument) |

### 1.2 Process-wide objects at namespace scope

| Symbol (file) | Holds | Written / read | Safety | Order | Becomes |
|---|---|---|---|---|---|
| `g_config` (`src/AZURE2.cpp`, `api/src/pybind_api.cpp`, `src/coulomb_table_main.cpp`) | a pointer to the run's `Config` (the binding swaps it per call) | read by `CoulFunc`, `ShftFunc`, `Kinematics.h` (relativistic switch), `UnitExternalCapture` | none: two sessions on two threads would race | yes in principle: a session's call reads the other's config if two calls run at once (read); not seen in one thread (A-B-A run) | deleted (done, branch 2): a pair carries the run's settings of its functions (`PPair::FunctionSettings`, given to its constructor by `CNuc::Fill`); `Kinematics::Channel`, `CoulFunc`, `ShftFunc` and the frame rules read them from the pair; `Kinematics::Relativistic` deleted |
| `g_ecAmplitudeCache` (`src/ECAmplitudeCache.cpp`) | the EC amplitudes on the reaction-rate/extrapolation grid | the binding swaps it per session; `UnitExternalCapture` reads | none across threads | yes if two sessions run at once (read) | the run's (done, branch 2): `RunScratch::ecAmplitudes`, read through `RunECAmplitudes(configure)`; the binding's swap deleted |
| `g_shiftDerivMemo`, `g_shiftDerivMutex` (`src/CNuc.cpp`) | dS/dE and d²S/dE² at level energies (Park), keyed by every input | `CNuc::CalcShiftFunctions` | mutex | no: exact key (read) | the run's (done, branch 2): `RunScratch::shiftDerivatives` |
| `g_terms`, `g_termsMutex`, `g_nodes`, `g_lookups`, `g_direct`, `g_nodesReport`, `g_checkMutex`, `g_maxDeviation` (`src/ECTable.cpp`) | the EC-integral table: nodes per term, filled on first use (`shared_future`), counters | EC amplitudes on the table | mutex + futures | no: a node's value is a pure function of its key (read) | the term table is the run's (done, branch 2): `ECTermStore` in `RunScratch::ecTerms` (filled on first use behind its lock: the terms are not known before the first evaluation); counters → diagnostics (milestone 5) |
| `g_files`, `g_bytes`, `g_stats`, `g_mutex` (`src/InputFiles.cpp`) | the bytes of data files read with `cacheInputFiles` (pyazr), keyed by path and stamp | data readers | mutex | no: re-read when the stamp moves (read) | allowed (the director, 2026-11-12; branch 3): the entries are the bytes of a file keyed by absolute path and file stamp (size, modification time), immutable once read, behind a lock; not a hazard, because a file changed between two reads has a new stamp and is read again, and no value of a calculation is stored; `cache_stats()` counts the process |
| `g_sharedGeneration` (`src/CoulombTable.cpp`) | moves when a missing table file is forgotten | `Shared` | atomic | no | deleted with the per-thread lookup |
| `SharedTables()::t`, `SharedFailed()::f`, `SharedMutex()::m` (`src/CoulombTable.cpp`) | the Coulomb tables loaded from files | `Shared` | mutex | no: immutable after load (read) | allowed (an immutable table read from a file) |
| `CoulombTable::hits_`, `misses_` (`include/CoulombTable.h`), `report` (`Shared`) | lookup counters printed at exit | every table lookup | atomic | no (diagnostic) | diagnostics, allowed until section 2 is done |
| `g_ns`, `g_calls`, `g_count`, `g_depth` (`src/Profile.cpp`) | `AZURE_PROFILE` accounting | timers | atomic | no (diagnostic) | allowed (diagnostics; section 2) |
| `g_teCombinerNs` (`src/AZUREGrad.cpp`) | gradient timing | `AZURE_GRAD_TIMING` | atomic | no | diagnostics |
| `g_jacGradForCheck`, `g_deferCheckExit` (`src/AZURECalc.cpp`) | the analytic gradient kept for `AZURE_GRAD_CHECK`, a deferred exit | gradient check | none | no output number (diagnostic) | diagnostics |
| (done 2026-11-13) the 17 variables of the running MCMC in `src/AZURECalcMCMC.cpp` (`g_mcmc_calc`, `g_sample_file`, ... `g_using_rwa_parameters`, the two progress callbacks) | | | | | a run record of `RunMCMCSampling` and members of the sampler; numcmc's callbacks are `std::function` |
| `g_last_printed_step`, `g_start_time` (`src/AZURE2.cpp`) | CLI progress of the MCMC | CLI | none | no | the same run object |
| `kSection` (`src/ModelData.cpp`) | the string "modelData" through a mutable pointer | readers | – | no | `const char *const` (done) |

### 1.3 Function-level and class-level statics

| Symbol (file) | Holds | Becomes |
|---|---|---|
| `told` (`src/CNucExternalCapture.cpp`, `WarnOnce`), `warned` ×2 (`src/HandOver.cpp`), `warned` (`src/ApparatusParameters.cpp`), `ecNoted` (`src/TargetSetup.cpp`), `reported` (`src/ModelPass.cpp`, `src/PointPass.cpp`), `warned` (`src/AMatrixFunc.cpp`, `InvertMatrices`), `told` (`src/AZUREGrad.cpp`) | "print this warning once per process": a second run in the same process (pyazr, the web UI) stays silent | notices: `RunNotices` held by the run's `Config`, one key per message (done) |
| `nChecked`, `calls` (`src/AZUREGrad.cpp`), `nDump`, `evals`, `pc`, `jc`, `jChecked`, `gradCalls`, `gradChecked`, `g2Calls`, `rjCalls` (`src/AZURECalc.cpp`), `calls` (`src/TargetSetup.cpp`), `worst` (`src/AMatrixFunc.cpp`, `AZURE_CHANNEL_VERIFY`) | call counters and worst values of diagnostic prints (`AZURE_GRAD_CHECK`, `AZURE_GRAD_TIMING`, `AZURE_DUMP_NONFINITE`, …) | diagnostics record (section 2); allowed until then |
| `fileMutex`, `traceMutex`, `progressMutex` (`src/AZURECalc.cpp`) | serialise a dump file, a trace and the progress line on `std::cerr`/the run's stream | allowed: a lock of a process-wide stream holds no value |
| `m`, `n`, `h` (`include/GSLException.h`, `GslQuiet`) | GSL's process-wide error handler, saved and restored around a scope | allowed: GSL's handler is process-wide by GSL's design (external) |
| `counter` (`src/ALevel.cpp`) | the stamp generator of `ALevel` (keys of the channel-space caches) | deleted with those caches (milestone 1) or allowed (an atomic counter that names objects) |
| `dummy` (`src/PolarizationFunc.cpp`, `AmplitudeMatrix::At`) | a zero returned by reference for an index out of range, writable by the caller | a `const` zero (done) |
| `kinds`, `kindName`, `basisName`, `names`, `parts` (5 files) | arrays of string literals through mutable pointers | `const char *const` (done) |

### 1.4 `mutable` members

| Members (file) | Holds | Becomes |
|---|---|---|
| `cnuc_pool_`, `edata_pool_`, `pool_mutex_`, `pools_initialized_` (`include/AZURECalc.h`, `include/AZURECalcMCMC.h`) | `max(4, threads)` clones of `CNuc` and of the working copy with its `Results`, taken by each evaluation of a fit | deleted (milestones 3–4): an evaluation reads the problem and writes its `Results` and scratch |
| `output_compound_`, `output_data_`, `output_results_`, `output_mutex_` (`include/AZURECalc.h`) | a clone kept for the iteration snapshots of a fit | deleted with the pools (the snapshot is an evaluation into its own `Results`) |
| `fixedMask_`, `g2CacheP_`, `g2Cache_`, `lastGradP_`, `lastGrad_`, `lastGradMutex_` (`include/AZURECalc.h`) | the free mask; Minuit's G2 computed with the last gradient at its parameters; the last gradient for a dump | the fitter's state (owned by the minimiser adapter, non-`const` methods): milestone 4 |
| `all_physical_`, `physical_`, `rwa_`, `all_rwa_`, `all_indexes`, `fixed_`, `paramNames_`, `paramKinds_`, `priorMeans_`, `priorStds_`, `usePriors_`, `freeKinds_`, `priorCov_`, `useRWA_` (`include/AZURECalcMCMC.h`) | the sampler's parameter layout, filled lazily in `const` methods | built at construction (milestone 4) |
| `levelGammas_`, `levelEnergies_`, `levelDecoupled_`, `shiftFunctions_`, `cached_max_levels_`, `cached_max_channels_` (`include/AMatrixFunc.h`) | buffers of the level-matrix fill | done (branch 2): plain members, written by the non-`const` fill |
| `inverse_valid_`, `inverse_`, `work_` (`include/LevelMatrixSolver.h`) | the explicit inverse formed on first request from the LU factors | done (branch 2): `Inverse`, `Bilinear` non-`const` (and `AMatrixFunc::GetAMatrix`, whose `const_cast` went) |
| `ec_shifted_` (`include/Site.h`) | per site, the EC amplitude at the last shifted energy | a field of the hand-over's table (EC unit amplitudes per mesh energy, roadmap Phase 6 box 1) |
| `channelFunc_` (`include/ThmVertexBoundary.h`) | a `ChannelFunc` whose calls are non-`const` (its last-value cache) | done (branch 2): `OnShellAt` non-`const` |
| `warned` (`ThmDistortion.h`, `ThmDwVertex.h`, `ThmLineshape.h`, `ThmOptions.h`) | warn-once flags of option objects | notices |
| `cache_mutex_` (`include/ECAmplitudeCache.h`) | guards the session's amplitude cache | done (branch 2): a plain member, the reader non-`const` |
| `nuisanceParams_`, `nuisance_`, `nuisanceRead_`, `results_` (`api/include/AZUREAPI.h`) | the session's nuisance layout read lazily; the session's `Results` written by `const` evaluations | the session owns its `Results` and passes them (non-`const` evaluation methods): milestone 3 |

### 1.5 Constant tables (stay; the allow-list's permanent rows)

`const` tables reached through a function (`Quantities()`, `Kernels()`, `Units()`, `Schema::Sections()`,
`DiagnosticCodes()`, `ParamKindTable()`, `ApparatusRows()`, `Objective::Engine()`, the mass and element tables): not
reported by the scanner (they are `const`), nothing to do.  The Coulomb table loaded from a file, the SRIM data
(`stopping/`, read once into `const` tables), GSL's error handler: on the allow-list with a reason.

### 1.6 Caches that are members of objects an evaluation owns (not `mutable`, found by reading)

| Member | Holds | Order | Becomes |
|---|---|---|---|
| `CoulFunc::lLast_ … coulLast_`, `tLast*` | the last call's waves and P, S, φ | no: exact key (read) | part of the object, which is the caller's (an energy task's); stays as scratch |
| `ALevel::shift_cache_energy_` | S_c at the last level energy | no (read) | parameter part of the evaluation (`CNuc`), milestone 3 |
| `AMatrixFunc::chan_cap_cache_`, `lastInEnergy_` | channel-capture terms at the last energy | no | scratch of the task |
| `SegmentState::lastEnergyShift…`, `shiftFnValues` | the shift the points carry | yes in the past (fit_threads.md cause A), now hand-over state of a working copy | the evaluation's hand-over table (milestone 3) |
| `WorkingCopy` overflow nodes (`NodeOverflow`), shift tables | nodes and channel rows of moved points | no (`AZURE_NODE_CHECK`) | the evaluation's hand-over table (milestone 3) |
| `TargetSetup` `gridNuc_` | a `CNuc` clone for grid placement | no | read from the problem (milestone 3) |
| `ECAmplitudeCache::cache_` | EC amplitudes per grid energy | no | a table of the problem built at set-up (EC unit amplitudes per mesh energy) |

## 2. The environment switches

The survey as it was made at the start of Phase 6 (2026-11-12).  Since the fourth branch (section 8) the list is the
table of `support/src/Environment.cpp` (39 engine switches; the survey's 44 counted `HOME`, `USERPROFILE`, the two
variables of the converter's test programs and the stopping test's own reading), checked against the manual's list in
`docs/source/reference/command_line.rst` (check 40); the "Proposal" column below was the survey's, the decisions of
section 6 are the PI's and none is taken: every switch kept its name and its effect.  At the survey: 44 names were read
with `getenv` (plus `HOME`/`USERPROFILE` for the default Coulomb-table path), each an input that no file and no
`--describe` showed.  "Changes": **numbers** (an output number moves), **route** (another algorithm,
numbers equal to rounding or to a stated tolerance), **diag** (prints or checks only; may exit), **test** (injects a
failure for a test), **speed**.  "Used by": tests (count of files under `tests/`), the director's scripts, nobody.

| Name | Read in | Changes | Used by | Proposal |
|---|---|---|---|---|
| `AZURE2_COULOMB_TABLE` (+`HOME`, `USERPROFILE`) | `CoulombTable.cpp` | numbers (1e-7: table vs exact) | tests 7, director | `<config>` key / CLI `--coulomb-table PATH` shown by `--describe`; the variable stays as the default's source |
| `AZURE2_SRIM_DATA_DIR` | `stopping/src/Tables.cpp` | numbers (which SRIM files) | tests 3 | the same: a key with the variable as default |
| `AZURE_OPTIONS_FILE` | `AZURE2.cpp` | numbers (any option) | tests 1 | CLI `--options FILE` |
| `AZURE_REFINE_QUADRATURE` | `CommandLineValues.cpp` | numbers (fold quadrature ×N) | tests 6 | CLI option (exists as a value of the command line; the variable goes) |
| `AZURE2_EC_TABLE_NODES` | `ECTable.cpp` | numbers (EC table density) | nobody | goes (the `<config>` key `ecTableNodesPerDecade` exists) |
| `AZURE2_SWAP_COUPLING_ORDER` | `PolarizationFunc.cpp` | numbers (A_y coupling order) | nobody | goes (the convention record of `ARCHITECTURE.md` §5 owns this) |
| `AZURE_BOUND_SHIFT_NUMERICAL`, `AZURE_NO_NEUTRAL_CLOSED_FORM`, `AZURE_SHIFT_DE_1EV` | `ShftFunc.cpp`, `CoulFunc.cpp` | numbers (older routes of S, dS/dE) | nobody | go |
| `AZURE_LEVEL_SPACE` | `ChannelSpace.cpp` | route | tests 3 | diagnostics record |
| `AZURE2_THM_FD`, `AZURE_GRAD_ADJOINT`, `AZURE_GRAD_DENSE`, `AZURE_TE_WEIGHTS_FD` | `AZUREGrad.cpp`, `AZURECalc.cpp`, `GradEnergyTasks.cpp`, `AMatrixFunc.cpp` | route (gradient by differences or another adjoint; a fit's path can move) | tests 3 / 2 / 0 / 0 | diagnostics record (the first two used by tests); the last two go |
| `AZURE_SERIAL_CHI2` | `EnergyPass.cpp` | speed (serial pass; same bits) | tests 1 | diagnostics record |
| `AZURE_GRAD_CHECK`, `_H`, `_IDX`, `_N`, `_ROWS`, `_SEGS`, `AZURE_SHIFT_GRAD_CHECK`, `AZURE_TE_WEIGHTS_CHECK` | `AZURECalc.cpp`, `AZUREGrad.cpp` | diag (may exit) | tests 49 / 8 / 0 / 22 / 10 / 0 / 2 / 0 | diagnostics record (one `--check-gradient` option later) |
| `AZURE_PROFILE`, `AZURE_GRAD_TIMING`, `AZURE_GRID_TIMING` | `Profile.h`, `WorkingCopy.cpp`, `AZURECalc.cpp`, `AZUREGrad.cpp`, `TargetSetup.cpp` | diag (timing) | tests 5 / 0 / 1 | diagnostics record |
| `AZURE_PLAN_CHECK`, `AZURE_NODE_CHECK`, `AZURE_POOL_CHECK`, `AZURE_CHANNEL_VERIFY`, `AZURE2_EC_TABLE_CHECK` | `CalculationPlan.cpp`, `NodeTable.cpp`, `AZURECalc.cpp`, `AMatrixFunc.cpp`, `ECTable.cpp` | diag (checks) | tests 4 / 3 / 0 / 1 / 2 | diagnostics record (`AZURE_POOL_CHECK` goes with the pools) |
| `AZURE_TRACE_FCN`, `AZURE_DUMP_NONFINITE`, `AZURE_GRAD_DEBUG`, `AZURE_CHANNEL_DEBUG`, `AZR_DEBUG_GRID` | `AZURECalc.cpp`, `AZUREGrad.cpp`, `ChannelSpace.cpp`, `AdaptiveIntegrationGrid.cpp` | diag (prints) | tests 1 / 0 / 0 / 0 / 0 | the first two: diagnostics record; the last three go |
| `AZURE_POISON_RESULTS` | `Results.h` | diag (NaN-filled buffers) | tests 3 | diagnostics record |
| `AZURE_TEST_FAIL_E`, `AZURE_TEST_FAIL_EXIT` | `PointPass.cpp` | test (injected failure) | tests 3 / 2 | diagnostics record, marked test-only |
| `AZURE2_SPEC_INJECT`, `AZURE2_WRITE_CONVERT_FIXTURES` | `convert/tests/` | test | tests | stay (test programs, not the engine) |

The diagnostics record: one `struct Diagnostics` filled **once** at start from the environment (the CLI's `main`, the
binding's session constructor), held by `Config` and passed down; `--describe` prints the switches that are set;
the readers take it from their `Config`.  A new switch is one field and one line of the reader (section 6, decision
3).  Not moved in this branch: the switches read as `static const bool x = getenv(...)` are `const` after their first
read (no order effect), and the others are read on use; neither is state an evaluation writes.

## 3. What an evaluation costs today

Measured on the base engine (`work/director/bin/AZURE2-thm5`, the branch's start), with the Coulomb table, one run at a
time, `OMP_NUM_THREADS` = 1 or 3, in a fresh copy of the project; wall time and peak resident memory from
`/usr/bin/time` (10 ms resolution); the machine was shared (load average 8–10 on 24 cores during the runs, another
agent's jobs): medians of three where three were run.  Script `work/eval/measure_all.sh`, table
`work/eval/timing.tsv`.

Wall time in s (median of three at 1 thread; one run at 3 threads) and peak resident memory in MB:

| Project | calc 1 thread | calc 3 threads | MB at 1 / 3 | fit (`--use-gradient`) 1 thread | fit 3 threads | fit MB at 1 / 3 |
|---|---|---|---|---|---|---|
| T2 (`tests/upgrade/T2`) | 0.22 | 0.18 | 121 / 123 | 33.2 (454 calls, χ² 297.568) | 24.6 | 190 / 213 |
| 13N | 7.75 | 5.75 | 119 / 120 | 17.98 (70 calls) | 9.56 | 141 / 143 |
| natural_neon | 1.35 | 0.85 | 343 / 343 | 1109 (one run) | not run | 812 / – |
| 18O_p_a_thm | 0.07 | 0.07 | 96 / 96 | – | – | – |
| beam_profile_kernel | 16.3 (one run) | 14.5 | 131 / 131 | – | – | – |

A calculation's wall time below 0.3 s is mostly start-up and is below the resolution of a wall-clock comparison on a
loaded machine: for the speed gate the user-space instruction count (`perf stat -e instructions:u`, one run, load
independent to 0.1 %) is reported beside the wall time (section 7).  The natural_neon fit at 3 threads was not run
(18 minutes at 1 thread; the thread budget of this branch is 3).

### 3.1 Memory as a function of the thread count (the roadmap's gate)

Measured: a calculation's peak memory does not depend on the thread count (T2 121 → 123 MB, 13N 119 → 120 MB,
natural_neon 343 → 343 MB from 1 to 3 threads).  A fit holds more: T2 +69 MB, 13N +22 MB, natural_neon +469 MB over the
calculation at 1 thread, and T2 +23 MB more at 3 threads.  Read from the code (`AZURECalc::InitializePools`): a fit
makes `max(4, threads)` clones of `CNuc` and of the working copy with its `Results`, plus one pair for the iteration
snapshots (`output_compound_`, `output_data_`), so from 1 to 4 threads the copies are the same five and above 4
threads each thread adds one.  Inferred from the two measurements: one copy is about 14 MB for T2 and about 95 MB for
natural_neon, so a natural_neon fit at 12 threads would hold about 13 copies, 1.6 GB (not run).  The roadmap's gate
("memory one problem regardless of threads") is therefore not met by fits today; it is met by calculations.

### 3.2 Where a fit's time goes

`perf record` (200 Hz, DWARF call graphs) of the base engine's `--use-gradient` fits at 1 thread, inclusive shares
of the whole run (`work/eval/perffit-T2-0ILa`, `perffit-13N-0N07`):

| | T2 (33 s) | 13N (18 s) |
|---|---|---|
| making the pooled copies (`InitializePools`, `Clone`) | 0.02 % | 0.03 % |
| re-initialising a pooled copy per evaluation: parameters (`FillCompoundFromParams*`) | 0.02 % | < 0.01 % |
| … energy shifts and re-placed grids (`FillEnergyShiftsFromParams`) | 6.3 % | not seen |
| the mapping between copies | 0 (deleted 2026-11-11) | 0 |
| waiting on a mutex (`pthread_mutex_lock`, pool and output locks) | not seen (< 0.01 %) | not seen |
| iteration snapshots on the output copies (`WriteIterationOutput`) | 2.2 % | not seen |
| point passes (`CalculatePointsInParallel`, all callers) | 56 % | 7 % |
| gradient (`AZURECalc::Gradient`) | 30 % | 40 % |
| set-up (the rest of `AZUREMain`, external-capture integrals) | ≈ 10 % | ≈ 38 % |

`AZURE_GRAD_TIMING` on the same T2 fit: "parameters + energy shifts" 0.44–2.1 ms against 38 ms for the points of an
evaluation.  What could not be separated: the share of the energy-shift hand-over that a copy repeats only because
its history differs from another copy's (the work itself is needed once per evaluation in any design).  The pools
cost memory (section 3.1), not time: removing them buys memory and simplicity, and at most the 2 % of the snapshots.

### 3.3 What the per-thread memos buy

The per-thread memos were deleted all at once in a trial build (`work/eval/bin/AZURE2-v1`; outputs identical), and the
calculations measured; then the ones that matter came back as caller-owned scratch (this branch, `AZURE2-v2`).
User-space instructions of one calculation at 1 thread (wall time in parentheses, ms of CPU):

| Project | base | all memos deleted | memos as caller-owned scratch |
|---|---|---|---|
| T2 | 1.4157e9 (233) | 1.4363e9, +1.5 % (288) | 1.4281e9, +0.9 % (246) |
| 18O_p_a_thm | 6.147e7 (75) | 6.434e7, +4.7 % (86) | 6.254e7, +1.7 % (72) |
| 13N | 7.4466e10 (7752) | 1.6109e11, +116 % (18490) | 7.4510e10, +0.06 % (7730) |
| natural_neon | 9.824e9 (1374) | – | 9.855e9, +0.3 % (1362) |

The external-capture memos (entrance Coulomb waves at the radial nodes, final-state Whittaker values) carry 13N: they
stay, as the run's per-thread scratch.  Of the THM memos only the vertex boundary's S_c(E_min) shows in a profile
(1.6 % of the 18O run); it stays as scratch; the level-width memo, the THM Coulomb-term memo and the angular-coupling
table showed nothing and are deleted (the coupling is built by the caller per point).

## 4. The target

**The problem** (immutable, built once per fit stage, shared `const` by every evaluation and every thread): the
system (partitions, channels, Jπ blocks, the pathway table), the data set (`DataSet`: measurements, segments, target
effects, covariance groups, parameter layout), the calculation plan and its energy tasks, the mesh and the stage's node
table with every node's channel rows, the kernels' weights that do not depend on fitted parameters, the EC unit
amplitudes per mesh energy, the THM vertex tables, the Coulomb table, the stopping tables, the diagnostics record and
the thread budget.  What a fit stage changes between stages (re-planned grids) is a new problem, built by the stage
change, never written in place.

**An evaluation**: parameter vector θ (and ν) in, `Results` out.  It owns (i) the parameter-dependent part of the
model — the formal level parameters E_λ, γ_λc of every Jπ block after `level_parameters`, the boundary values and
S_c(E_λ) of Brune and Park, the external widths and the norms and shifts of the data sets; (ii) the hand-over table of
the energy-shifted data (the nodes and channel rows of moved points, `NodeOverflow` and the shift tables today); (iii)
its `Results` (unit records, model values, contracted values); (iv) per energy task, scratch: the level-matrix
buffers (`AMatrixFunc`, `LevelMatrixSolver`, `THMMatrixFunc`), the task's unit list, the channel-space arrays.  The
parallel loop is over the energy tasks (for a gradient, over the tasks of each pass), with a fixed-order reduction,
and nowhere else; the thread budget is one number passed in (`EvaluationContext::threads`), never read from OpenMP's
global state.

**`CNuc` and `WorkingCopy` in that picture.**  `CNuc` today mixes both: its pairs, J groups, channels, the pathway table,
the level *structure* and the convention are problem; the level energies and widths as fitted (`ALevel::fitE_`,
`fitGamma_`, the shift-function values, the external widths, `SqrtNF`) are evaluation.  The split is the remaining
step 7 of `cnuc_split.md` (the parameter vector and the set-up out of `CNuc`); until it is done, an evaluation needs
its own `CNuc` (the `cnuc_pool_` stays) — the first milestone below does not need it.  `WorkingCopy` is already split:
`DataSet` is problem (shared through `shared_ptr<const DataSet>`), `SegmentState` (norms, shifts, frames, shift
tables), the present target-effect values and the overflow nodes are evaluation, the `Unit`s are problem except
their slots' hand-over values.

**Minuit2's parallel gradient and the MCMC's walkers** each run several evaluations at once.  Each needs a
`Results`, its parameter-dependent model part (until `CNuc` is split: a `CNuc` copy) and its hand-over table — nothing
else, and nothing shared but the `const` problem.  The per-evaluation state is a value (`Evaluation`: the three parts
above) constructed from the problem in microseconds once the pools' deep copies of `Unit`s are gone; a pool of them is
then an optimisation of allocation, not of copying, and may stay as a plain `std::vector<Evaluation>` indexed by the
calling thread.

**pyazr and `api/`** hold today: the session's `Config` and the global pointer swapped per call (`g_config`), the
session's `ECAmplitudeCache` swapped into `g_ecAmplitudeCache`, the session's `Results` (`mutable results_`), its
`WorkingCopy`/`CNuc` and a lazily read nuisance layout.  In the target, a session holds a problem (`shared_ptr<const
Problem>`) and one `Evaluation`; the swap of globals goes (milestone 2).  The public Python functions return what
they return now.

### 4.1 What an evaluation of a fit writes today (the inventory for milestone 3, branch 3)

Read from the code of 55a7e2c9 (`WorkingCopy.h`, `HandOver.cpp`, `ShiftFrames.cpp`, `NodeTable.cpp`,
`UnitExternalCapture.cpp`, `AZURECalc.cpp`, `AZURECalcMCMC.cpp`, `AZUREAPI.h`).  An evaluation of a fit is
`FillCompoundFromParams` (the `CNuc`), `FillNormsFromParams` and `FillEnergyShiftsFromParams` (the hand-over into the
working copy), the point pass and the objective into a `Results`.  What each member of `WorkingCopy` is:

| Member of `WorkingCopy` | Problem or evaluation | Written by an evaluation when |
|---|---|---|
| `set_` (`DataSet`: segments, terms, target effects, covariance, parameter layout) | problem (shared `const` by the copies) | never |
| `plan_`, `capacity_`, `iterations_` (an atomic counter read through the master copy), `covFiles_`, `ec*_` (set-up of `intEC.dat`) | problem | never |
| `segmentStates_` (`SegmentState`: norm, shifts, last shifts, shift-function values, thickness; per term the frame rows at the present shift and the shift tables) | evaluation | every evaluation (norms); shifted segments (the rest) |
| `effectValues_` (`TargetEffectValues`: U_e, density, temperature, stopping scale, kernel values), `thmCoherentValues_` | evaluation | a fit of these parameters |
| `units_`, `slots_` (`Unit`, `Site`: the points and the sites of their folds) | **both**: the measurement's values (datum, error, original energy, angle) are problem; lab and c.m. energy and angle after the shift, Legendre rows, penetrabilities, shifts and phases, Coulomb amplitude, channel rows (`RecalcEDependentValues`), screening U_e (`SetScreening`, every evaluation with a screening effect), target thickness and stopping (density, stopping scale), the EC amplitude at the shifted energy (`ec_shifted_`), the slots of a re-placed fold are evaluation | shifted terms, fitted target-effect parameters, re-placed grids |
| `slotBlock_`, `slotWriter_`, `slotReaders_`, the result layout (`resultsSize_`, `modelSize_`, a unit's record offsets) | problem as laid out at set-up; evaluation when a re-placed grid lays the units out again (`PlaceSubPointGrids`, `LayoutResults`) | re-placed grids |
| `rowSets_`, `contractedSize_` | evaluation when slots move (rebuilt by the hand-over, `BuildRowSets`) | re-placed grids, shifts |
| `nodes_` (`NodeState`: the stage table, shared; the slots of the moved tasks and the overflow nodes) | the table is problem; slots and overflow are evaluation | shifts, re-placed grids |
| `meshAnchors_`, `meshSupports_`, `meshCentres_`, `gridNuc_` (a `CNuc` holding the observed level values the meshes are laid out for) | evaluation (functions of the level parameters and shifts) | a fit of a level that anchors a mesh, shifts |

**Finding.**  Section 4 above says "the `Unit`s are problem except their slots' hand-over values".  That is true of a
project without energy shifts, fitted target-effect parameters and re-placed grids; for the others the hand-over
rewrites most of a `Site` (its energies, angles, energy-dependent channel values, Coulomb amplitude, rows, stopping
and thickness) and may lay the slots out again.  Removing `WorkingCopy::Clone` therefore needs a split of the `Site`
layer: the measurement and the set-up's values stay in the problem's sites, and the values above become a table of the
evaluation, indexed by site, written by the hand-over and read by every reader of a site (the observables, the folds,
the gradient's cotangents, the node checks, the writers) — about 40 accessors of `Site.h` and their readers in
`ShiftFrames.cpp`, `TargetSetup.cpp`, `NodeTable.cpp`, `EnergyPass.cpp`, `GradEnergyTasks.cpp`, the observables and
the output.  It is the "hand-over table" of section 4 at the size of the sites, not of the overflow nodes alone, and two
of its files are in another branch's work tonight (`ShiftFrames.cpp`, `TargetSetup.cpp`).  It is not done in this
branch (see section 7.5); what this branch does is what does not need it:

1. The pools hold one evaluation per evaluation running at the same time, made when first needed, instead of
   `max(4, threads)` copies made at the first evaluation.  An evaluation is a pair (a `CNuc` copy, a working copy
   with its `Results`) taken by one caller at a time.  A fit with the analytic gradient (one evaluation at a time,
   its parallelism inside the point pass) then holds one, at any thread count; Minuit2's numerical gradient and the
   MCMC's walkers hold one per thread that evaluates.
2. (Planned, not done in branch 3: see section 7.5.)  `AZUREAPI` reads its nuisance layout once, at construction
   and at a rebuild; its evaluating methods write the session's `Results` as non-`const` methods.
3. Kept, with the reason: `WorkingCopy::Clone`, `WorkingData` and the snapshot copies (`output_*`) of a fit (a
   snapshot runs `TransformOut` on its `CNuc`, which a pooled one must not see), `AZURE_POOL_CHECK` and
   `AZURE_TRACE_FCN=fresh` (they compare a pooled copy with a fresh one, which still exist), `Site::ec_shifted_`
   (a value of the `Site` split).

## 5. Milestones

Each is a commit series that leaves the tree passing (all projects × 3 bases identical to the reference, unit tests,
`run_tests.sh`, the eight fits identical, the session test of section 7.2).

1. **No `thread_local`, no warn-once statics, no mutable pointers to literals** (this branch).  Gate: check 38 with
   the allow-list; the session test; speed within 2 % (section 3).  Deletes the memos that buy nothing, moves the others
   to caller-owned scratch.  Risk: a memo hid a recomputation whose rounding differs from the memo's value (a memo
   keyed by an address could return a value computed by another object); found by the full comparison — a moved
   number would show there.  About −250 / +200 lines.
2. **No process-wide pointers** (`g_config`, `g_ecAmplitudeCache`, the swap in the binding) and the shift-derivative
   memo of `CNuc`: the readers take their `Config` from their caller; the EC amplitude cache becomes a member of the
   session's problem.  Gate as 1, plus two sessions on two Python threads at once (the binding releases the GIL)
   giving the bytes of one.  Risk: a reader deep in `Kinematics.h` (the relativistic switch) must get `Config` through
   ~15 call chains; a forgotten one is a compile error, not a silent change.  About ±400 lines.
3. **An evaluation of a fit writes only its `Results`, its hand-over table and task scratch** (the brief's part C):
   `edata_pool_`, `WorkingData`, `WorkingCopy::Clone`, the output copies go; the working copy's evaluation state
   (`SegmentState`, target-effect values, overflow nodes) becomes an `Evaluation` value next to the `const` problem.
   Needs no decision of section 6.  Deletes `include/WorkingData.h`, `WorkingCopy::Clone`, `AZURE_POOL_CHECK`, the
   pool code of `AZURECalc.cpp` (~250 lines).  Risk: the hand-over (energy shifts, re-placed grids) writes the working
   copy today; a copy whose history differed gave other bits before (fit_threads.md, cause A); the eight fits
   identical file by file and `AZURE_TRACE_FCN=fresh` (every pooled value equal to a fresh one) find it.
4. **`CNuc` split, the parameter vector out** (`cnuc_split.md` step 7, another agent's area): `cnuc_pool_` goes; the
   sampler's and the fitter's state become non-`const` objects of the minimiser adapter; the MCMC globals become a run
   object.  Gate: memory of a fit at 1, 3 and 12 threads equal within the `Results` size per thread.
5. **The diagnostics record** (section 2): one `Diagnostics` struct read once; the counters of section 1.3 become
   its members (or go); `--describe` lists the set switches.  Can run in parallel with 2–4.  Done in the fourth
   branch (sections 7.6 and 8): one table and one reader in `support/`, `Config::environment`, `RunDiagnostics`.

## 6. Decisions for the director

1. **Environment switches that change numbers** (`AZURE2_COULOMB_TABLE`, `AZURE2_SRIM_DATA_DIR`,
   `AZURE_OPTIONS_FILE`, `AZURE_REFINE_QUADRATURE`): become command-line options and `<config>` keys shown by
   `--describe`, the variables kept as defaults of the first two (recommended), or stay variables.  Users see this.
2. **Switches that select an older route** (`AZURE2_SWAP_COUPLING_ORDER`, `AZURE_BOUND_SHIFT_NUMERICAL`,
   `AZURE_NO_NEUTRAL_CLOSED_FORM`, `AZURE_SHIFT_DE_1EV`, `AZURE2_EC_TABLE_NODES`, `AZURE_GRAD_DENSE`,
   `AZURE_TE_WEIGHTS_FD`, `AZR_DEBUG_GRID`, `AZURE_GRAD_DEBUG`, `AZURE_CHANNEL_DEBUG`): no test uses them;
   recommended: delete.
3. **One diagnostics record** read once at start and shown by `--describe` (recommended), versus keeping `getenv` at
   the point of use.  An agent then sees every switch that was on in a run's log.
4. **The meaning of the thread setting.**  Today `OMP_NUM_THREADS` sets both the energy-task loop and Minuit2's
   parallel gradient (nested: an evaluation inside the gradient runs its tasks serially), and the pools hold
   `max(4, threads)` copies.  Recommended: one `--threads N` (and `threads=` in pyazr) passed down as the budget, with
   `OMP_NUM_THREADS` its default; Minuit2's own loop stays OpenMP (Minuit2 is not patched) and gets the same number.
5. **Memory for speed.**  Memos that buy more than 2 % stay as caller-owned scratch (no slower engine); a memo whose
   removal costs less than 2 % on the gate projects is deleted.  If the director prefers fewer lines over the last
   percent, the remaining scratch (EC radial blocks, Whittaker values) can go at a measured cost (section 3.3).

## 7. Record

### 7.1 State counts

`state_scan.py` (as committed) on the start 88d55f17 and on this branch; the allow-list after:

| Kind | before | after | on the allow-list as `table:` | as `pending` |
|---|---|---|---|---|
| thread_local | 26 | 0 | 0 | 0 |
| static (block or class) | 40 | 26 | 10 | 16 |
| namespace scope | 44 | 42 | 0 | 42 |
| mutable member | 53 | 49 | 0 | 49 |
| total | 163 | 117 | 10 | 107 |

(The first count in the introduction, 42 namespace-scope variables, was made before the scanner caught `T a(0),
b(0)`; 44 is the committed scanner's.)  Pending rows by milestone: 1 (solver buffers written by `const` methods) 10,
2 (global pointers, process-wide memos) 13, 3 (pooled working copies, the API's lazy state) 11, 4 (`CNuc` pools, the
fitter's and the sampler's state, MCMC globals) 44, 5 (diagnostic counters) 29.

### 7.2 Two calculations in one process (`tests/pyazr/session_isolation_test.py`)

Run on the base engine (`work/b-thm5/pyazr-out`, log `work/eval/aba_base.log`): **no leak found**.  13N (A),
18O_p_a_thm (B), hybrid_potential (C) and 13N with its first pair's radius moved to 5.0 fm inside the session (R), in
the orders A B A, B A B, A C A, A R A, R A R, each at 1 and 3 threads: every evaluation of a project gives the bytes
of that project alone in a fresh process (cross sections of every point, χ² per segment, objective and analytic
gradient), and 1 thread gives the bytes of 3.  The address-keyed memos that leaked before (the THM S_c memo, the
external-capture amplitude cache) had been fixed by full keys and per-session objects; the state that is left is
safe in one thread of one session at a time, not by construction.

### 7.3 Done in this branch

**Part 2 (no hidden state), first milestone, partly.**  Every `thread_local` is gone; the warn-once statics print
once per run; the mutable pointers to literals are `const`; the rest is on the allow-list with its milestone.
- Energy-task scratch: `include/TaskScratch.h`, `src/TaskScratch.cpp` (`TaskScratch`: the point pass's and the
  gradient's `AMatrixFunc` and `THMMatrixFunc`, the task's unit list, `ChannelSpace::Memo`; `TaskScratchSet` in
  `Results::scratch`, one slot per OpenMP thread, empty on copy).  `ChannelSpace::Memo` (`include/ChannelSpace.h`)
  replaces the three per-thread maps; an `AMatrixFunc` uses its own unless given the task's
  (`ShareChannelMemo`).  `ThmPointModelGradient` takes a solver owned by its parallel region.
- The run's per-thread scratch: `include/RunScratch.h`, `src/RunScratch.cpp` (`RunScratch`, `ThreadScratch`: an
  `ECScratch` — GSL workspace, `CoulRadialMemo`, `WhitScaledMemo`, rate grids — and a `ThmShiftMemo`), held by
  `Config::scratch`.  `CoulFunc::SetRadialMemo`, `WhitFunc::SetScaledMemo` take the caller's tables; the FW/GW memo
  is a member of the integral's parameters.
- Deleted: the THM level-width memo, the THM Coulomb-term memo, the per-thread table of the THM angular coupling
  (`ThmAngularCouplingOf`, built by the caller), the per-thread Coulomb-table lookup and its generation counter.
- Notices: `Config::printedOnce`, `Config::FirstTime(key)`, `Config::FirstTimeIn(log, key)` (set-up records keep the
  log); 12 flags (CNuc external capture, hand-over ×2, apparatus, target set-up, model pass, point pass, singular
  level matrix, gradient debug, four THM option objects).
- Check 38 (`tests/consistency/repo_consistency_test.py`, `state_scan.py`, `allowed_state.tsv`); shown failing on a
  planted namespace variable, static, `thread_local` and `mutable` member (the check plants them itself) and on an
  injected process-wide EC amplitude memo in `src/UnitExternalCapture.cpp` (two findings).
- The session test (`tests/pyazr/session_isolation_test.py`, ctest `pyazr_session_isolation`): passes on the base
  and on this branch with the same digests; with the injected memo it fails (A R A and R A R, 1 and 3 threads).

**Speed after part 2** (one thread; the same build as the full checks).  User-space instructions, one run each
(load independent): T2 calculation +0.9 %, 18O_p_a_thm +1.7 %, 13N +0.06 %, natural_neon +0.3 %; T2
`--use-gradient` fit 3.25865e11 → 3.26007e11 (+0.04 %), 13N fit 1.874908e11 → 1.875266e11 (+0.02 %).  Wall time,
medians of three, base / branch, runs interleaved at load 12–16 on 24 cores (another group's builds and tests):
T2 0.47 / 0.60 s, 13N 10.11 / 8.87 s, natural_neon 2.00 / 2.20 s, 18O 0.16 / 0.14 s, T2 fit 44.0 / 43.2 s, 13N fit
21.0 / 21.0 s; the spread of the three runs of one engine (up to 50 % for T2's calculation) is larger than any
difference, so the wall times neither show nor exclude a 2 % change; the instruction counts are within it.
Peak memory: calculations within 0.1 MB; T2 fit 189.4 → 192.9 MB (+3.6 MB), 13N fit 140.6 → 140.6 MB.  Inferred
cause of the T2 increase: the solvers are now one set per `Results` (per pooled copy, five in a fit) instead of one
per thread; it goes with the pools (milestone 3).

Not done in this branch: the 107 pending rows (milestones 1–5); the `g_config` and `g_ecAmplitudeCache` pointers;
the diagnostics record; part 3 (section 5, milestone 3).

### 7.4 Done in the second branch (`feature/upgrade-6-evaluation-2`)

**Milestone 2, except the input-file bytes.**  `g_config` is deleted with its three definitions and 30 test copies: a
pair carries the run's settings of its functions (`PPair::FunctionSettings`: relativistic kinematics, the use and the
file of the Coulomb table), passed to its constructor by `CNuc::Fill`, and every former reader takes them from its pair
(`Kinematics::Channel`, `CoulFunc`, `ShftFunc`, the frame rule of `Segment.cpp`, `Site::QLabBoostAt`,
`Site::CalculateCrossSectionConversionFactor`, the extrapolation seeds, the channel table's signature, the
shift-derivative key).  The table is loaded on first use as before.  The binding's `ConfigScope` (the swap of
`g_config` and `g_ecAmplitudeCache` around every call) is deleted.  The run's scratch (`RunScratch`) holds, besides
its per-thread part, the external-capture term table (`ECTermStore`), Park's shift-derivative memo
(`ShiftDerivativeMemo`) and the `--no-ec-table` amplitude store (`ECAmplitudeCache`).  Within one run the Coulomb
table is fixed, so the keys that hold only whether a table is used are complete (said in `ECTable.h` and
`RunScratch.h`).  Not moved: `InputFiles` (the byte store is keyed by absolute path and file stamp, behind a lock;
moving it per session changes what `cache_stats()` counts with two sessions in a process).

**Milestone 1 completed.**  No `mutable` in `AMatrixFunc.h`, `LevelMatrixSolver.h`, `ThmVertexBoundary.h`; the methods
that write the buffers are non-`const`; the `const_cast` of `AMatrixFunc::GetAMatrix` is gone.

**Two sessions at the same time** (`tests/pyazr/concurrent_sessions_test.py`, ctest `pyazr_concurrent_sessions`,
label slow, about 11 min at 2 threads).  The evaluations release the GIL (`py::call_guard<py::gil_scoped_release>` or
an inner `gil_scoped_release` block in `api/src/pybind_api.cpp`); the `Session` constructor (the set-up) keeps it, so
one session's set-up runs beside the other's evaluation, and two evaluations run side by side.  Found while writing it:
the engine resolves a project's relative paths against the process's working directory, which pyazr sets per call
(`_in_dir`); two sessions with different directories on two threads race on it (the test uses absolute paths).

Verification of the branch head (build `work/b-eval2`, merged with `feature/upgrade-4-data`): whole tree x 3 bases
identical to `tab-gridm-1` (827 / 800 / 827, run.log included); the converted-tree gate 628 / 601 / 628; the eight fits
identical to `work/b-head3`; T1, T2, T3, 13N at 1 and 3 threads identical; `tests/run_tests.sh` 160 / 0.  User-space
instructions at 1 thread against `b-head3`: T2 1.42771e9 -> 1.42820e9 (+0.03 %), 13N 7.45089e10 -> 7.45087e10,
natural_neon 9.8539e9 -> 9.8642e9 (+0.10 %), 18O_p_a_thm 6.25213e7 -> 6.25218e7, T2 `--use-gradient` fit 3.25702e11
-> 3.25903e11 (+0.06 %), 13N fit 1.875305e11 -> 1.875930e11 (+0.03 %).  Allow-list: 107 pending rows -> 88 (milestone
1: 10 -> 0, milestone 2: 13 -> 4, the input-file bytes; 3: 11; 4: 44; 5: 29).

### 7.5 Done in the third branch (`feature/upgrade-6-evaluation-3`)

**Milestone 3, the part that does not need the `Site` split (section 4.1).**  The copies an evaluation of a fit or
of the MCMC runs on are made when an evaluation finds no idle one, not `max(4, threads)` (fit) or twice the hardware
threads (MCMC, 48 on the 24-core machine) at the first evaluation; `InitializePools` and `pools_initialized_` are
deleted.  A fit with the analytic gradient holds the master pair, one evaluation pair and (from the first snapshot)
the snapshot pair at any thread count.  Kept: `WorkingCopy::Clone`, `WorkingData`, the pools' stacks, the snapshot
copies, `AZURE_POOL_CHECK`, `Site::ec_shifted_` (section 4.1, item 3); `AZUREAPI`'s `mutable` members (not reached).

**Relative paths (done, branch 3, continued).**  A project's relative paths are resolved against one directory held
by the run, never against the process's working directory, and pyazr no longer changes the working directory.
- One rule, one function: `azure2::experiment::PathToOpen(directory, path)` (`experiment/.../ProjectPaths.h`; an
  absolute path, or no directory, as written; an empty path is the directory itself).  The engine's
  `Config::OpenPath(path)` applies it with `Config::baseDirectory`, which `Config::ReadConfigFile` sets to the
  project file's directory (made absolute once, at start, against the working directory: the file name is the
  caller's) unless the caller named one; pyazr names the session's `cwd` (`RuntimeOptions.base_directory`).
- Every opening site takes `OpenPath` and every message prints the path as the project writes it: data files
  (`Segment::Fill`), covariance files and the data files they list (`DataCovariance::Read`, the groups' keys), the
  user-defined resolution file (`TargetEffect`), a stopping table (`TargetStopping`, the directory passed in), the
  output and checks files (`DataOutputFiles`, `AZUREOutput`'s prefix, `DataCheckFiles`, `CNucReport`,
  `AZUREParams`, `AZURECalc`, `AZUREMain`, `ReactionRate`, `ThmReports`, `DataCovarianceGroups`), `intEC.dat` and
  its signature (`ExternalCaptureSetup`, `WorkingCopy`, `AZUREAPI`), the MCMC files (`AZURECalcMCMC`, the CLI's
  overwrite), `kinflag.dat`, `nonfinite_dump.txt`, and a file named at a menu prompt (parameter, temperature,
  integrals file: `AZURE2.cpp`).  The output-directory check refuses an empty name as before.
- The experiment library's `projectDir` is used for opening only (`BuildContext::ProjectPath`, `Grids.cpp`):
  messages, `GridValues::file`, `ExtrapolationSpec::gridFile`, `StoppingSpec::tableFile` keep the written name;
  `--describe` still prints the resolved data file, as before.  The engine passes its `baseDirectory` (so
  `--print-predictions` and `pyazr.predictions`, which set the project up through `ReadConfigFile` as a run does,
  read relative to the project file now).
- pyazr: `_in_dir` and the `_Engine` wrapper are deleted; the `Session` gets the directory; a relative
  `coulomb_table` is made absolute against it (it was the session's directory's before, through the `chdir`).
- Unchanged, deliberately: THM tables and the prior covariance were already relative to the project file
  (`ThmRelativePath`, `PriorCovariance::ResolvePath`); the project file and `--coulomb-table`/`AZURE_OPTIONS_FILE`
  given on the command line are the caller's, relative to where the command is given; `pyazr.validate` and
  `pyazr.describe`, with `--validate` and `--describe` (held to the same text by `tests/pyazr/validate_test.py`),
  keep reading data files relative to the process's working directory: changing them changes what a public pyazr
  function returns, a decision for the director.
- Tests: `tests/run_from_elsewhere` (the CLI started from the parent and from an unrelated directory writes the same
  output and checks files and reads an integrals file named at the prompt as the project's; fails on the base); `tests/pyazr/concurrent_sessions_test.py` pair R S (13N and 18O_p_a_thm
  with their relative paths, each session also writing its output files, two threads at once).

**Two defects, by a test first.**  (a) The signature of `intEC.dat` held neither whether the Coulomb table is used nor
which file: it now holds the table's content identity (below) or "exact Coulomb functions"
(`tests/ec_table_signature`).  Every signature of a run with a table changes, so an existing `intEC.dat` is computed
again once.  The table is named by what it holds, not by its file (`CoulombTable::ContentIdentity`: format,
interpolation order, l range, the two axes, eta*rho maximum, generator text and a 64-bit hash of every tabulated
value, 8 bytes a step; computed on the first request and fixed after it, behind `std::call_once`): a copy of the table
at another path, or the same table on another machine, reuses the file; a table with other values does not.  Cost: one
pass over the 88 MB table, 1.2e8 user-space instructions (T1 2.9467e10 -> 2.9589e10, +0.41 %; 13N 7.4509e10 ->
7.4631e10, +0.16 %), paid only by a run that sets up external-capture integrals (T2: 1.42829e9 -> 1.42839e9).  (b) `--validate` and `--describe` read `relativistic` already (the project reader's `RelativisticFlag`);
`tests/describe_relativistic` holds it (a THM data set under `relativistic=1` is refused by both).

**Input-file bytes**: the four rows are permanent (`table:`), with the reason of section 1.2.

**Diagnostics record (milestone 5)**: not started (done in the fourth branch, section 7.6).

### 7.6 Done in the fourth branch (`feature/upgrade-6-evaluation-4`)

**Milestone 5: the environment record and the run's diagnostics** (design: section 8).
- One table, one reader: `support/` (library `azure2_support`, uses nothing; the stopping library and the engine link
  it), 39 switches: 14 *numbers*, 1 *cost*, 22 *prints*, 2 *test*.  `ProcessEnvironment()` reads them once; check 39
  refuses `getenv` anywhere else in the engine, its libraries and the binding (test programs under `tests/` aside);
  check 40 keeps the table and the manual's list (`command_line.rst`, "Environment switches", all 39 documented) the
  same, kinds included.
- Through the run's `Config` (`Config::environment`, from the record by `Reset`): `AZURE2_COULOMB_TABLE` (the default
  table file; to the pairs through `PairFunctionSettings`), `AZURE2_SRIM_DATA_DIR` (given by `TargetStopping` to the new
  `Stopping`/`StoppingProvenance` overloads that take SRIM's data; the library's own default `SrimTables()` reads the
  record), `AZURE2_EC_TABLE_NODES`, `AZURE_BOUND_SHIFT_NUMERICAL`, `AZURE_NO_NEUTRAL_CLOSED_FORM`, `AZURE_SHIFT_DE_1EV`
  (pair settings), `AZURE_LEVEL_SPACE`, `AZURE2_THM_FD`, `AZURE_GRAD_ADJOINT`, `AZURE_GRAD_DENSE`,
  `AZURE_TE_WEIGHTS_FD`, `AZURE_SERIAL_CHI2`.  `AZURE_OPTIONS_FILE` and `AZURE_REFINE_QUADRATURE` are read from the
  record where the options are read (the CLI's `parseOptions`, the binding's session).  Through the record:
  `AZURE2_SWAP_COUPLING_ORDER` (the A_y amplitude matrix of a stored record has no run object), the 22 printing
  switches and the two test hooks.
- Findings (switches read at every call, or set inside a process): section 8.3.  The engine's own `setenv`
  (`AZURE_POOL_CHECK` -> `AZURE_SERIAL_CHI2`, `AZURE_GRAD_CHECK` -> `AZURE_GRAD_ADJOINT`) is gone: a `Config` copy
  (`AZURECalc::Chi2OnWith`) and a gradient route (`AZURECalc::GradientBy`).  `tests/pyazr/thm_adjoint_test.py` runs its
  difference half in a child process started with `AZURE2_THM_FD=1`.
- Counters: `RunDiagnostics` (`Config::diagnostics`, shared by copies of the `Config`): the profile account
  (`AZUREProfile::Account`; every `Timer`, `Count`, `ScopedReport` takes the run's, `Config::profile`, null without
  `AZURE_PROFILE`), the Coulomb table's lookups (`CoulombTableUse`, given to the lookups by the pairs' settings;
  printed by `RunDiagnostics::Finish` when the run ends -- the CLI's `Config` at the end of `main`, or `ExitRun` before
  an early `std::exit` -- where the process's exit printed it before), and the call counters and worst values of the
  gradient, grid, pool and channel checks.  The external-capture table's counters are members of the run's term store
  (`ECTermStore`).  The gradient kept for `AZURE_GRAD_CHECK` and its deferred exit are an argument (the route).
- Allow-list: 29 rows of milestone 5 removed (and the two globals of the gradient check); the scanner now also
  reports a `std::once_flag` member (kind `once`): `CoulombTable::ContentIdentity`'s is a permanent row.  Pending
  53 (milestone 3: 11, milestone 4: 42), permanent 15.
- `--describe` (text form) ends with the switches that are set ("Environment switches set:", one line
  `NAME=value (kind)` each); `pyazr.describe` returns the description without it.  `--validate`, `--describe`,
  `pyazr.validate`, `pyazr.describe` read data files relative to the project file's directory
  (`azure2::experiment::ProjectDirectoryOf`, also the run's rule in `Config::ReadConfigFile`).
- Speed (user-space instructions, one run at 1 thread, base `work/b-gecm` / this branch): T1 2.95891e10 / 2.95854e10,
  T2 1.43138e9 / 1.42572e9 (-0.4 %), 13N 7.46322e10 / 7.46164e10, natural_neon 9.8673e9 / 9.8505e9, T2
  `--use-gradient` fit 1.737368e11 / 1.736046e11, 13N fit 1.885643e11 / 1.885169e11: none slower.
- Not done (another agent's file this night, `api/src/AZUREAPI.cpp`): the unused `intEC_cache` file name and the
  lazily read nuisance layout (`mutable` members) of the session.

## 8. The environment record and the run's diagnostics (milestone 5): design

Status: design (2026-11-13, branch `feature/upgrade-6-evaluation-4`, from 7f85e0f7), written before the code; done as
designed, with the choices recorded in section 7.6.

### 8.1 One table, one reader

- **The table** (`support/include/azure2/support/Environment.h`, `support/src/Environment.cpp`, library
  `azure2_support`, which depends on nothing and which the stopping library and the engine link): one row per switch
  with its name, type (flag, text, number, integer), default, kind, one line of what it does and the page that
  documents it.  Kinds: **numbers** (an output number can move: other physics, another route equal to rounding, or
  another input), **cost** (the same bits, another time), **prints** (prints, dumps or checks; may exit), **test**
  (injects a failure for a test).  The table is a `constexpr` array, so it cannot change.
- **The reader**: `azure2::support::ReadEnvironment(lookup)` is the only call of `getenv` in the engine, the libraries
  and the binding; `ProcessEnvironment()` returns the record read once, on its first call, and never written after (a
  `const` object with a dynamic initialiser, like the constant tables: the environment is the process's by nature).
  The first call is the CLI's `main` before anything else, and in pyazr the first `Config` a session builds.  The
  record also keeps the home directory (`HOME`, else `USERPROFILE`) for the default Coulomb-table path.  The record
  keeps each value as text with whether it was set at all, and each reader applies the test it applied before (set at
  all; set, non-empty and not "0"; a number by `atof`/`atoi`), so every switch keeps its effect.
- **A consistency check** (number 39) refuses `getenv` in `src/`, `include/`, `api/`, `physics/`, `experiment/`,
  `stopping/`, `io/`, `convert/`, `support/` outside `support/src/Environment.cpp` (test programs under `*/tests/`
  are not the engine), and every row's name must appear in `docs/source/reference/command_line.rst` (number 40;
  the table there is checked against the code's table both ways).

### 8.2 Config or record, per switch

The switches that can move a number and that code holding a `Config` reads are copied into the run's `Config` when the
`Config` is made (`Config::Reset`, from the record), and the code reads the `Config`: making one of them a command-line
option or a pyazr option later is one assignment in `parseOptions` or in the binding's `RuntimeOptions`.  The others
(prints and checks, the test hooks, the cost switch with no run object at hand) read the record where they are used.

| Switch | Kind | Read through |
|---|---|---|
| `AZURE2_COULOMB_TABLE` | numbers | `Config::environment.coulombTableDefault`, to the pairs' settings |
| `AZURE2_SRIM_DATA_DIR` | numbers | `Config::environment.srimDataDirectory`, given by `TargetStopping` |
| `AZURE_OPTIONS_FILE`, `AZURE_REFINE_QUADRATURE` | numbers | the record, where the options are read |
| `AZURE2_EC_TABLE_NODES` | numbers | `Config::environment.ecTableNodes` |
| `AZURE2_SWAP_COUPLING_ORDER` | numbers | the record (`PolarizationFunc`: no run object there) |
| the three shift-function routes | numbers | `Config::environment`, to `PairFunctionSettings` |
| `AZURE_LEVEL_SPACE` and the four gradient routes | numbers | `Config::environment` |
| `AZURE_SERIAL_CHI2` | cost | `Config::environment.serialPointPass` (the pool check: a `Config` copy) |
| the 22 printing and checking switches, the two test hooks | prints, test | the record, where used |

The library's own default (`stopping::SrimTables()`, for callers without a run) reads the record.

### 8.3 Switches read more than once or set inside a process (findings)

Read at every call today: `AZURE_SERIAL_CHI2`, the gradient switches of `AZURECalc.cpp`/`AZUREGrad.cpp`,
`AZURE_GRID_TIMING`, `AZURE_PROFILE` in `WorkingCopy.cpp`, `AZURE2_COULOMB_TABLE` (each `CoulombTable::Shared("")`),
`AZURE_REFINE_QUADRATURE` (each pyazr session).  Set inside a process: (a) the engine itself, `AZURE_POOL_CHECK` sets
`AZURE_SERIAL_CHI2` around its reference and `AZURE_GRAD_CHECK` sets `AZURE_GRAD_ADJOINT` around the second gradient:
both become explicit (a `Config` copy for the first, the gradient's route as an argument for the second), and the
`setenv` helpers go; (b) `tests/pyazr/thm_adjoint_test.py` sets `AZURE2_THM_FD` between two calls of one session: the
test runs its difference half in a child process started with the variable; (c) `tests/upgrade/results_owner_test.cpp`
sets `AZURE_POISON_RESULTS` before its first `Config`, which stays valid.  No public pyazr function changes.

### 8.4 The counters: owned by the run

`Config::diagnostics` (`include/RunDiagnostics.h`: a `RunDiagnostics` held through a `shared_ptr`, made by the
`Config`'s constructor and shared by its copies, as `scratch` is) holds what the diagnostics count: the profile
accounting (`AZUREProfile::Account`: stage times, counters, the nesting depth of reports), the Coulomb table's lookup
counts (given to the pairs' functions through `PairFunctionSettings`, printed when the run's `RunDiagnostics` ends, as
the process's exit printed them before: the CLI's run ends with the process), and the call counters and worst values
of the gradient, grid and channel checks.  The external-capture table's counters go into the run's term store
(`ECTermStore`, already in `RunScratch`), beside the terms they count.  The analytic gradient kept for
`AZURE_GRAD_CHECK` and the deferred exit are arguments of the second gradient call, not state.  The locks of the dump
file, the trace and the progress line stay (locks of process-wide streams, permanent rows).

### 8.5 `--describe` and the manual

`--describe` ends with the switches that are set (name, value, kind), from the table and the record; `pyazr.describe`
returns the project's description as before (its text stays the CLI's up to that block).  The manual's list
(`command_line.rst`, a section "Environment switches") is checked against the table (check 40).

## 9. The switches that move numbers as options of the program, and `--threads`: design

Status: design (2026-11-13, branch `feature/upgrade-6-switches-options`), written before the code. The PI approved
(2026-11-13): the switches of kind numbers become command-line options shown by `--help` and `--describe`; one
`--threads=N`; nothing deleted, renamed or re-defaulted (the variables keep working with the same effect).

### 9.1 The options, one per row

The row of `SwitchTable()` gains the option's name (`SwitchRow::option`). A switch of type flag takes the bare option
(its value is the text "1", which every reader counts as set and on, as the variable set to 1); the others take
`--name=VALUE`. Names follow the existing flags (lower case, words joined by `-`, `=` before a value).

| Variable | Option | Project or run | Why |
|---|---|---|---|
| `AZURE2_COULOMB_TABLE` | `--coulomb-table=FILE` (exists) | run | a file of this machine (equal to exact to 1e-7) |
| `AZURE2_SRIM_DATA_DIR` | `--srim-data-dir=DIR` | run | a directory of this machine, not a choice of the project |
| `AZURE_OPTIONS_FILE` | `--options-file=FILE` | run | a file of more options |
| `AZURE_REFINE_QUADRATURE` | `--refine-quadrature=K` | run (doubt) | a convergence check of the fold quadrature |
| `AZURE2_EC_TABLE_NODES` | `--ec-table-nodes=N` | run (doubt) | a convergence check of the external-capture tables |
| `AZURE2_SWAP_COUPLING_ORDER` | `--swap-coupling-order` | run (doubt) | a test of convention invariance |
| `AZURE_BOUND_SHIFT_NUMERICAL` | `--bound-shift-numerical` | run | an older route of the bound-state S |
| `AZURE_NO_NEUTRAL_CLOSED_FORM` | `--no-neutral-closed-form` | run | an older route of the neutral bound S |
| `AZURE_SHIFT_DE_1EV` | `--shift-derivative-1ev` | run | the older dS/dE route |
| `AZURE_LEVEL_SPACE` | `--level-space` | run | a route equal to rounding |
| `AZURE2_THM_FD` | `--thm-by-differences` | run | a gradient route (moves a fit's path, not a calculation) |
| `AZURE_GRAD_ADJOINT` | `--gradient-adjoint` | run | a gradient route equal to rounding |
| `AZURE_GRAD_DENSE` | `--gradient-dense` | run | a gradient route equal to rounding |
| `AZURE_TE_WEIGHTS_FD` | `--fold-weights-by-differences` | run | a gradient route |

No switch becomes a `<config>` key: each is a file or directory of the machine, a check of a numerical tolerance, an
older or a diagnostic route. Doubts for the PI: the two convergence checks (a project that needs K = 4 to converge
could say so in its file) and the coupling-order test (`ARCHITECTURE.md` §5 gives conventions to a record). None of the
names collides with a flag of another meaning (`--ec-table` / `--no-ec-table` and `--use-gradient` stay what they are;
the parser matches the whole option, or the text before `=`).

`AZURE_SERIAL_CHI2` (kind cost) stays a variable only: it changes no bit (it is checked by `AZURE_POOL_CHECK` to give
the same bits), so a user or an agent reading `--describe` to know what can move a number does not need it; it is a
tool to localise a defect of the parallel pass, like the printing switches. That is right as long as the parallel pass
is bit-identical to the serial one, which the tree's 1-against-N-threads comparison checks.

### 9.2 One parser, one record

- **The record carries the options.** `ProcessEnvironment()` is already read first thing by the program's `main`
  (nothing reads it during static initialisation: every reader is a function or a function-level static, all reached
  after `main` starts). `main` now makes it with `StartProcessEnvironment(arguments)`: the switch options of the
  command line, then those of the options file, then the variables. Every reader keeps reading the record (or the
  run's `Config::environment` copied from it), so an option has exactly the effect of its variable, also where no run
  object is at hand (`PolarizationFunc`'s coupling order, the stopping library's default SRIM tables). The record stays
  immutable after it is made; a second `StartProcessEnvironment` with options after the record exists is a programming
  error (thrown, never silently ignored). A pyazr session makes the record from the environment alone, as today.
- **One parser, generated from the table** (`ParseSwitchOption` in `support/src/Environment.cpp`): an argument is a
  switch option when it is a row's option (a flag row: the whole argument; a value row: the text before `=`). The
  program's `parseOptions` skips them (it has no `else if` per switch); an unknown option stays an error.
- **Value checks in one place**: the row says the range of a whole number that today's reader checks (only
  `AZURE_REFINE_QUADRATURE`, 1..64); `SwitchWholeNumber(record, id)` checks the value in force, whatever its source,
  and a bad value exits 2 naming the option or the variable it came from. The rows whose reader does not check today
  (`AZURE2_EC_TABLE_NODES`: `atoi`, a bad value is 0 = the default) keep that reading: checking it would change what
  the switch does (not approved). A value row given without a value (`--srim-data-dir=`) is refused (exit 2).
- **Precedence** (stated and tested): the option on the command line, then the options file, then the environment
  variable, then the default. No switch has a `<config>` key (9.1). A flag row given twice: the last wins. Where this
  differs from the flags of section 3.12 (`project_format3.md`): the options file's lines are applied **after** the
  command line's flags, so for the options of the calculation a flag in the file wins over one on the command line
  (as before; not changed here); for the switches the command line wins, as the brief asks. Listed for the PI.
- `--coulomb-table=FILE` is one row, not two mechanisms: it is the option of `AZURE2_COULOMB_TABLE` (the record) and,
  as before, turns `use_coulomb_table` on and names the run's table.

### 9.3 `--describe`, `--help`, the manual

- `--help` lists the switch options from the table (option, variable, what it does), after the flags of the options.
- `--describe` (text) ends, after the project's description, its options block (both of which `pyazr.describe`
  returns as well) and the stopping of its targets, with a block headed `settings of this run and machine (value,
  source):`: `  threads = N (source)` and one line per switch of kind numbers,
  `  --refine-quadrature=K (AZURE_REFINE_QUADRATURE) = 4 (option)`, with its value in force (the row's "unset" text
  when not set) and its source (`option`, `options file`, `environment`, `default`); then the other switches that are
  set (`Environment switches set:`, kinds cost, prints, test). The tests that compare a committed description cut the
  text at that header: it depends on the machine. `--json`: the members `threads` and `switches`
  (`{"AZURE_REFINE_QUADRATURE": {"option": "--refine-quadrature=K", "value": "4", "source": "option"}, ...}`).
- The manual's list of switches gains the option column; check 40 (the list against the table) checks it; the flag
  count check counts the switch options; a row of kind numbers without an option, or an option of another kind, fails
  the consistency test (check 43).

### 9.4 `--threads=N`

`--threads=N` (N ≥ 1, the checked whole-number parse of `--max-calls`; a bad value exits 2) sets the OpenMP thread
count of the process (`omp_set_num_threads`) in `main` before any parallel region and before the project commands.
Without it: today's behaviour (`OMP_NUM_THREADS`, else the runtime's default). With both: the option wins and one line
says so. The value in force and its source (`--threads`, `OMP_NUM_THREADS`, the runtime's default) are printed by
`--describe` and by the run after the banner. `OMP_NUM_THREADS` is read through the record (as `HOME`), not by a new
`getenv`. `--build-coulomb-table` takes `omp_get_max_threads()` and so honours it. The Coulomb-table tool
(`azure2-coulomb-table`) takes `--threads=N` with the same check and the same meaning (its older `--threads N` form
stays). The sampler's menu question for its thread count (`mcmcParams.nthreads`, `src/AZURECalcMCMC.cpp`) belongs to
the sampler: it should take `--threads` as its default and stop asking (listed, not done here).

### 9.5 pyazr (open items for pyazr 3.0)

A session reads the record from the environment, as before. For a later step: (1) a session keyword per switch of
kind numbers read through `Config::environment` (one assignment after `Config::Reset`, from the same table); (2) the
two switches read from the process record where no run object exists (the coupling order, the stopping library's
default tables) stay process-wide: a session cannot change them, and should say so; (3) `pyazr.describe()` should
return the same switch block (value, source); (4) the threads of a session (`OMP_NUM_THREADS` today) as one keyword
with the same meaning as `--threads`.
