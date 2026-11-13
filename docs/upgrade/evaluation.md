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
| `g_files`, `g_bytes`, `g_stats`, `g_mutex` (`src/InputFiles.cpp`) | the bytes of data files read with `cacheInputFiles` (pyazr), keyed by path and stamp | data readers | mutex | no: re-read when the stamp moves (read) | not moved in branch 2: a session-owned store would make `cache_stats()` (documented "for the whole process") count one session; left for a decision (section 6) |
| `g_sharedGeneration` (`src/CoulombTable.cpp`) | moves when a missing table file is forgotten | `Shared` | atomic | no | deleted with the per-thread lookup |
| `SharedTables()::t`, `SharedFailed()::f`, `SharedMutex()::m` (`src/CoulombTable.cpp`) | the Coulomb tables loaded from files | `Shared` | mutex | no: immutable after load (read) | allowed (an immutable table read from a file) |
| `CoulombTable::hits_`, `misses_` (`include/CoulombTable.h`), `report` (`Shared`) | lookup counters printed at exit | every table lookup | atomic | no (diagnostic) | diagnostics, allowed until section 2 is done |
| `g_ns`, `g_calls`, `g_count`, `g_depth` (`src/Profile.cpp`) | `AZURE_PROFILE` accounting | timers | atomic | no (diagnostic) | allowed (diagnostics; section 2) |
| `g_teCombinerNs` (`src/AZUREGrad.cpp`) | gradient timing | `AZURE_GRAD_TIMING` | atomic | no | diagnostics |
| `g_jacGradForCheck`, `g_deferCheckExit` (`src/AZURECalc.cpp`) | the analytic gradient kept for `AZURE_GRAD_CHECK`, a deferred exit | gradient check | none | no output number (diagnostic) | diagnostics |
| `g_mcmc_calc`, `g_sample_file`, `g_sample_count`, `g_nwalkers`, `g_current_step`, `g_samples_written_to_file`, `g_stop_requested`, `g_step_offset`, `g_last_*` (5), `g_param_names`, `g_mcmc_for_callbacks`, `g_using_rwa_parameters` (`src/AZURECalcMCMC.cpp`) | the running MCMC: sample file, step, walkers, last values for progress callbacks | the ensemble sampler's callbacks (a C interface without a user pointer) | none | no output of a calculation; the chain's file depends on them (read) | a run object passed through the callbacks' user pointer (milestone 4, with numcmc's retirement) |
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
| `all_physical_`, `physical_`, `rwa_`, `all_rwa_`, `all_indexes`, `fixed_`, `paramNames_`, `paramKinds_`, `parametersInitialized_`, `priorMeans_`, `priorStds_`, `usePriors_`, `freeKinds_`, `priorCov_`, `useRWA_` (`include/AZURECalcMCMC.h`) | the sampler's parameter layout, filled lazily in `const` methods | built at construction (milestone 4) |
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

44 names are read with `getenv` (plus `HOME`/`USERPROFILE` for the default Coulomb-table path).  Each is an input that
no file and no `--describe` shows.  "Changes": **numbers** (an output number moves), **route** (another algorithm,
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
   its members (or go); `--describe` lists the set switches.  Can run in parallel with 2–4.

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

