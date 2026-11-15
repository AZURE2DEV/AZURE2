# Sanitizer runs of the upgrade test projects: setup and findings

This document records what the ThreadSanitizer (TSan) and AddressSanitizer (ASan) runs of the Phase 0 CI item check
(`ROADMAP.md`, Phase 0; `RULES.md` §2.4), what each run needs from the toolchain and why, and every report the runs
produced on `feature/upgrade` with its diagnosis and a proposed fix. The runs check the parallel and memory behaviour
of the present engine on the test projects T1–T3 and on the tests of point mapping, energy shifts and target
integration; a report in AZURE2's own code is a finding, fixed on its own `fix/` branch, never suppressed.

## 1. What runs

`scripts/sanitize_upgrade.sh thread|address [build-dir]` configures a build with `-DAZURE2_SANITIZE=<kind>`
(`CMakeLists.txt`: every target, the vendored Minuit2 and coul included, built with `-fsanitize=<kind>`, frame
pointers and line tables; no GUI), then runs at `OMP_NUM_THREADS=8` (`SAN_THREADS`):

| section | runs | what it exercises |
|---|---|---|
| `upgrade` | T1, T3: calculation (menu 1), extrapolation (3), reaction rate (5), fit from the start model with the analytic gradient (2, `--use-gradient`); T2: calculation and a fit bounded to `SAN_FIT_CALLS` = 30 function calls (`--use-lm`, Minuit2's numerical gradient, as `scripts/bench_upgrade.sh`) | all three lines of `ARCHITECTURE.md` §1; Minuit2's OpenMP gradient calling the objective from eight threads |
| `checks` | 26 `tests/<name>/check.sh`: `angle_mapping`, `component_mapping`, `component_subpoint_angle`, `shared_point_tables`, `shared_tables_race`, `shared_effect_pairs`, `upos_pair_index`; `shift_fn`, `shift_fn_groups`, `shift_groups`, `shift_penalty_zero`, `shift_sqrt_subpoints`, `shift_subpoint_angle`, `energy_shift_sqrt`; `target_auto_windows`, `target_density`, `target_effect_shared`, `target_thickness_exact`, `target_window_fallback`, `thick_target`, `convolution_unify`, `convolution_window`, `straggling_bohr`, `natural_target`, `natural_neon` | point mapping across segments and angles and the shared point tables; energy shifts (groups, functions, √E, sub-points, penalties, the parallel re-energising of `EData::FillEnergyShiftsFromParams`); target integration (thickness, density, windows, thick targets, straggling, convolution, natural targets) |
| `regress` | menu 1 of `beam_profile_kernel`, `beam_spread_ay`, `thick_target_resonance` | the target and profile folds that have no `check.sh` |
| `reference` | ctest `beam_profile_reference`, `udr_reference`, `adaptive_grid_reference`, `compound_stopping_reference` | the integrators on their own |
| `pyazr` | 12 tests in-process: `shift_fn`, `sqrt_shift`, `shift_first_call`, `jacobian_shift_columns`, `target_density`, `thickness`, `natural_target`, `doppler`, `data_covariance`, `composite_terms`, `parallel_points`, `shared_state` | the `_azure2` module: sessions, the parallel point pass of the API, shifts and targets through pyazr |

Each sanitizer report is written to `<report-dir>/<step>/<kind>.<pid>` (`log_path`), so a report from a child
process (a `check.sh` run, a python test) is counted too. The script exits 1 if any report exists, 2 if a step
failed without a report, 3 if the toolchain or build is unusable. Leak detection is off (`detect_leaks=0`; §4.5).

## 2. Toolchain, and why

**OpenMP under TSan.** TSan sees synchronization only through instrumented atomics and the pthread functions it
intercepts. GCC's libgomp implements barriers, the end of a work-sharing loop and reductions with raw futex system
calls, so a GCC TSan build reports every value written inside a parallel region and read after it. The runs use clang
with LLVM's libomp and its OMPT tool libarcher (`OMP_TOOL_LIBRARIES`), which declares to TSan the happens-before
order of barriers, tasks, critical sections and reductions. The script refuses a TSan run on libgomp, finds
`libarcher.so` next to the libomp the binary loads (or `SAN_ARCHER`), and checks with `ARCHER_OPTIONS=verbose=1`
that the tool attached before it runs anything.

**Third-party code is instrumented too.** A race is reported only between instrumented accesses, and synchronization
done in uninstrumented code is invisible; Minuit2 and coul are therefore built with the same flags (the option is set
before `add_subdirectory`). libomp itself is not instrumented; its own bookkeeping produces one report (§4.4).

**The python module.** `_azure2` is a shared object loaded by an uninstrumented python. Clang does not link the
sanitizer runtime into shared objects, so python runs with the runtime in `LD_PRELOAD`
(`libclang_rt.<kind>-x86_64.so` or `libclang_rt.<kind>.so`, found with `-print-file-name`; for GCC the module's own
`DT_NEEDED`). `PYTHONMALLOC=malloc` lets ASan see the allocations, and `ignore_noninstrumented_modules=1` keeps TSan
from reporting python's own threads. A CLI started by a python test gets `LD_PRELOAD` removed (the executable
carries its own runtime; two runtimes in one process abort). The tests import pyazr from the directory two levels
above their own file, so they run from a mirror of the checkout whose `pyazr/` holds the sanitizer module; the build
writes the module to `<build>/pyazr-module` (`PYAZR_OUTPUT_DIR`), never into the checkout.

**Suppressions.** `scripts/sanitizers/tsan.supp` (run time) and `scripts/sanitizers/tsan.ignorelist` (compile time,
clang) hold only reports traced to third-party code, each with its reason. A run-time `race:` entry matches any frame
of either stack; an entry naming a Minuit2 frame would therefore also hide every AZURE2 race in the objective that
Minuit2 calls. Third-party races whose stacks contain AZURE2 frames are for that reason removed at compile time, one
function at a time. ASan can suppress only reports raised inside interceptors; `asan.supp` is empty.

**CI.** `.github/workflows/upgrade-sanitizers.yml`: Ubuntu 24.04, clang 18 with `libclang-rt-18-dev` (the sanitizer
runtimes) and `libomp-18-dev` (libomp and libarcher), a TSan and an ASan job, nightly at 04:23 UTC and on manual
dispatch, and on a push that changes the workflow, the script or the suppressions. The TSan job sets
`vm.mmap_rnd_bits=28`: the 32-bit mmap randomization of recent runner kernels collides with TSan's fixed shadow
mapping. Reports and logs are kept as a job artifact for 14 days. A step failure without a sanitizer report (exit 2)
is a warning there, since `build.yml` gates test failures; reports and setup errors fail the job. GitHub runs a
schedule only from the workflow file on the default branch: until the file is there, the nightly does not fire.

**The runs recorded here** used the CRC cluster (HTCondor, 8 cores per job): clang 23.1.4, with compiler-rt (ASan,
TSan) and libomp/libarcher built from the LLVM 23.1.4 sources because the cluster's clang installation carries
neither; and GCC 15.2 (ASan only, libgomp). A clang installation without compiler-rt and libomp needs both built from
the same LLVM release and passed with `-resource-dir` and `-L`/`-rpath` (`CFLAGS`, `CXXFLAGS`, `LDFLAGS`) and
`-DOpenMP_omp_LIBRARY=` (`SAN_CMAKE_ARGS`).

## 3. Results on `feature/upgrade` (engine sources of 356c6770), 2026-11-06

| run | steps | steps failed | sanitizer reports |
|---|---|---|---|
| clang 23.1.4 ASan (24 min with the build) | 54 | 1 (`adaptive_grid_reference`, §4.6; no report) | none |
| GCC 15.2 ASan (43 min) | 54 | 1 (the same) | none |
| clang 23.1.4 TSan, no suppressions (50 min) | 54 | 13 | 10 distinct report sites in 34 report files: §4.1–§4.4 |
| clang 23.1.4 TSan, with `tsan.supp` and `tsan.ignorelist` (49 min) | 54 | 5: T3-fit, T2-fit30, `shared_tables_race`, pyazr `shift_fn` (reports), `adaptive_grid_reference` (§4.6) | 8 distinct report sites in 25 report files: §4.1 (48 reports) and §4.2 (80) only |
| clang 23.1.4 TSan, as above, on 5b021003 with the fixes of §4.1 (fecf4924) and §4.2 (c556f052), 2026-11-07 (54 min) | 54 | 2: `shared_point_tables` (peak RSS limit (d) under TSan, §4.2; no report), `adaptive_grid_reference` (§4.6) | none |
| clang 23.1.4 ASan + LeakSanitizer, `upgrade` section | 10 | 10 (exit 1 from LSan) | one direct leak per run, §4.5 |

ASan is silent on every step with both compilers: no heap or stack overflow, use after free, or double free in T1–T3,
the 26 checks, the regression projects, the reference tests or the pyazr tests at 8 threads.

## 4. Findings

### 4.1 Fit iteration counter and progress line written from every thread (AZURE2; real race; fixed)

TSan, in every fit that calls the objective from Minuit2's numerical gradient (T3-fit, in its Hesse step; T2-fit30;
the fits of `shared_tables_race`, `shift_fn`, `shift_fn_groups`, `energy_shift_sqrt`, `natural_target`):

```
Write of size 4 by thread T2:   EData::Iterate()            src/EData.cpp:2023   iterations_++
                                AZURECalc::operator()       src/AZURECalc.cpp:566
                                MnUserFcn::operator()       external/minuit2/src/MnUserFcn.cxx:41
                                Numerical2PGradientCalculator::operator() (.omp_outlined)  Numerical2PGradientCalculator.cxx:163
Previous write by thread T4:    the same stack
Read of size 4 by thread T1:    EData::Iterations() const   src/EData.cpp:970    (AZURECalc.cpp:565)
Location: heap, the EData of AZUREMain (include/AZUREMain.h:23)

Write of size 8 by thread T3:   std::ios_base::width        AZURECalc::operator()  src/AZURECalc.cpp:709
Previous write by thread T2:    the same                    Location: global std::cout
```

`AZURECalc::operator()` reads the shared fit counter (`EData::iterations_`, a plain `int`), increments it, and every
tenth value prints a progress line to `configure().outStream` (`std::cout`); every hundredth writes a snapshot
(`WriteIterationOutput`, which is mutex-protected). Minuit2 calls the objective from all OpenMP threads during the
numerical gradient (`Numerical2PGradientCalculator.cxx:116`), so the read-increment is not atomic and the stream
state (`width`) and buffer are written concurrently.

Effect: increments are lost and two threads can read the same counter value, so the iteration number shown, the
choice of which evaluations print a progress line and which write the snapshot are not deterministic; concurrent
writes to `std::cout` may interleave or garble the progress line. No χ², parameter or output file of the fit result
depends on the counter.

**Fixed** by fecf4924 (`fix/fit-iteration-counter`). `iterations_` is a `std::atomic<int>`; `EData::Iterate()`
returns the value it replaced (`fetch_add`, relaxed: the counter orders no other memory), and `AZURECalc::operator()`
takes its number from that one call. The progress line and the flush after it are written under a mutex local to
`operator()`; `WriteIterationOutput` keeps its own. The clang 23.1.4 TSan run of the branch merged with §4.2's fix
(54 steps, 8 threads; §3) reports nothing: the 48 reports at `EData.cpp:970`, `EData.cpp:2023` and
`AZURECalc.cpp:709` of the run before are gone, and T3-fit and T2-fit30 pass. A T1 fit from the start model with
Minuit2's numerical gradient (73 progress lines) prints the same lines, the same iteration numbers and the same χ²
as the build before at one thread; at eight threads the numbers printed are exactly 10, 20, …, 730, each once.

Proposed fix (`fix/fit-iteration-counter`): make `iterations_` a `std::atomic<int>`, let `EData::Iterate()` return
the value it replaced (`fetch_add(1, std::memory_order_relaxed)`), so that each evaluation gets a unique number from
one atomic operation, and replace the two calls in `AZURECalc::operator()` by that one; copy the value with `load()`
in `EData::Clone` (`EData.cpp:5527`). Print the progress line under a mutex (the one of `WriteIterationOutput` or a
second), or only from the thread that holds `omp_get_thread_num() == 0`. Test: TSan on T3-fit reports nothing at
`EData.cpp` or `AZURECalc.cpp:709`; the progress line count of a T2 fit of N calls is ⌊N/10⌋.

### 4.2 Copy-on-write ownership test of the point tables (AZURE2; not a race under the C++ memory model, invisible to TSan; fixed)

TSan, in every run that re-energises points in parallel (`shared_tables_race`, `shift_fn`, `shift_fn_groups`,
`shift_groups`, `pyazr shift_fn`; five reports per occurrence, one per member of `EDependentTables`):

```
Write of size 8 by thread T4:   std::vector<int>::assign    EPoint::EDependentTables::Clear()  include/EPoint.h:160-164
                                EPoint::ClearEDependentRows()           src/EPoint.cpp:1705   (SoleOwner(edep_) true)
                                EPoint::RecalcEDependentValues          src/EPoint.cpp:1726
                                ESegment::UpdatePointEnergyWithShift    src/ESegment.cpp:1136
                                EData::FillEnergyShiftsFromParams (.omp_outlined)  src/EData.cpp:5426
Previous read of size 8 by main thread:
                                std::vector<int>::size()    EPoint::ClearEDependentRows()  src/EPoint.cpp:1697   (SoleOwner(edep_) false)
                                ... the same call chain ...
Location: heap, an EDependentTables made by EPoint::MutableEDep (src/EPoint.cpp:1760) during EData::Initialize
```

Two points hold the same table (`edep_`, a `std::shared_ptr` shared after `EData::MapData` or a clone). The thread of
the first point finds the table shared, reads its row sizes to reserve a table of its own and drops its reference;
the thread of the second then finds itself the sole owner (`SoleOwner`, `src/EPoint.cpp:62`) and empties the table
in place. The order is guaranteed: the reference is dropped by libstdc++'s `_M_release`, an acquire–release
fetch-add on the use count (`bits/shared_ptr_base.h:165`), and `SoleOwner` reads the count with a relaxed load
(`use_count()`, `shared_ptr_base.h:216`) followed by `std::atomic_thread_fence(std::memory_order_acquire)`
(3fba1192). A relaxed load that reads the value written by a release operation, followed by an acquire fence,
synchronizes with it; the read of the first thread happens before the write of the second. TSan does not model
standalone fences (`compiler-rt/lib/tsan/rtl/tsan_interface_atomic.cpp`, `OpFence::Atomic`: "FIXME(dvyukov): not
implemented"), so it sees a relaxed load with no acquire and reports the pair. The libc++ `use_count()` (macOS) is
also a relaxed load and its release an acquire–release decrement; the argument holds there too.

No new reference is taken while the parallel loop runs: a point mapped onto a host copies the host's `edep_`
(`CopyEDependentValuesFrom`) in the host's own work item, after the host, on the same thread
(`ESegment::UpdatePointGroupEnergyWithShift`, `EData.cpp:5405-5426`), and `EData::MapData` runs after the loop.

The report is not suppressible (AZURE2 code), so the TSan job reports it until the ownership test is one TSan can
see. Proposed fix (`fix/edep-ownership-tsan`): replace the `use_count()` + fence test by an ownership count TSan sees,
with the same orderings: an `std::atomic<int> holders` in `EDependentTables` (and in the shared Legendre table),
incremented (relaxed) wherever `edep_` is copied, decremented (`memory_order_acq_rel`) wherever a point drops it,
and loaded with `memory_order_acquire` in `SoleOwner`; or hold the tables through a small handle class that does
this, so the copies and drops cannot be missed. Phase 6 removes the copy-on-write tables altogether
(`ROADMAP.md`); the fix is worth making before that only if the TSan job is to be green earlier. Test: TSan on
`shared_tables_race` and `shift_fn` reports nothing at `EPoint.h:160-164`; the race test still passes at 8 threads.
A check that the fix branch has not made a real race invisible: remove the acquire on the new counter and see TSan
report the pair again.

**Fixed** by c556f052 (`fix/edep-ownership-tsan`). `edep_` and `legendreP_` are held through `SharedTable<T>`
(`include/SharedTable.h`): a `std::shared_ptr` for the lifetime, and a count of holders in the table itself. A copy
of a handle increments the count (relaxed), a dropped handle decrements it with release order after its last access
to the table, and `SoleOwner()` loads it with acquire order. A load that reads 1 reads the last read-modify-write in
the count's modification order; each earlier release decrement heads a release sequence that continues through every
later read-modify-write on the count ([intro.races]), so the load synchronizes with the decrement of every former
holder, and the reads of those holders happen before the sole holder's write. These are release and acquire
operations on one atomic object, which TSan models. The weak references of `NodeTableCache` are
`SharedTable<T>::Weak`; `lock()` counts a new holder, as `weak_ptr::lock` did, and a published table is still never
written again.

Results, clang 23.1.4 TSan at 8 threads, on a temporary branch with this fix and §4.1's (`feature/upgrade` 5b021003
+ fecf4924 + c556f052): 54 steps, no report file (before: 25 files, 128 reports at 8 distinct AZURE2 sites, §3).
The control asked for above: with the acquire load replaced by a relaxed one, the `checks` section (25 steps) gives
5 report files with 75 reports, all at `EDependentTables::Clear` (`shift_fn`, `shift_fn_groups`, `shift_groups`),
the pair of this section; TSan therefore sees the ordering through the new count and only through it.

The change adds 8 bytes per table (the count) and one relaxed increment per copied handle. The peak RSS of the
4-thread fit of `shared_point_tables` (d), a sixteenth of `natural_neon`, in plain GCC builds: 306 and 311 MB on
`feature/upgrade` (no GUI), 324 and 326 MB on this branch, 329 and 332 MB on `fix/fit-iteration-counter`, which does
not touch the tables (both with the GUI); this branch adds no measurable memory. The check's limit, 3.0 times the
calculation's peak RSS, is near the measured ratio: the calculation took 103 MB in one run and 119 MB in the next on
the same binary, and the 332 MB run failed (d) once in plain builds. Under TSan the fit used 1133 MB before and 1162
and 1164 MB in two runs after (TSan keeps a synchronization record for each atomic object on which a release is done,
one more per table), above 3.0 × 383 MB = 1149 MB; that step fails under TSan without a report. In plain builds both
branches pass `tests/run_tests.sh` (131 of 131, `shared_point_tables` included).

### 4.3 Minuit2 call counter (third party; real race; compile-time ignorelist)

TSan, in every fit with the numerical gradient (38 reports in the first run):
`ROOT::Minuit2::MnUserFcn::operator()` increments `mutable int MnFcn::fNumCall` (`MnUserFcn.cxx:21`) from every
thread of `Numerical2PGradientCalculator::operator()` (`.omp_outlined`, `Numerical2PGradientCalculator.cxx:163`),
write against write, location the `MnUserFcn` on the main thread's stack (`AZUREMain.cpp:531`). Increments are lost:
the call count Minuit2 reports, and the `--max-calls` bound that uses it, are approximate under threads; no fitted
value depends on it. Minuit2 is vendored and not edited here (root `CLAUDE.md`). The entry in
`scripts/sanitizers/tsan.ignorelist` removes the instrumentation of that one function; the objective it calls stays
instrumented (§4.1 is reported through it). If upstream Minuit2 makes the counter atomic, the entry goes.

### 4.4 libomp suspend/resume bookkeeping (third party; run-time suppression)

TSan, in every run (24 reports in the first run, two per process): `pthread_mutex_init` in
`__kmp_suspend_initialize_thread` (a worker initializing its own suspend mutex) against `pthread_mutex_lock` in
`__kmp_resume_template` (the thread waking it), heap block from `___kmp_allocate`. Both stacks are inside libomp;
AZURE2 frames appear only as the code that opened the parallel region (`AZUREMain.cpp:175`,
`ReactionRate.cpp:631`, `EData.cpp:4648`). libomp orders the two through its own atomic flags, which TSan does not
see because libomp is not instrumented and libarcher annotates OpenMP constructs, not the runtime's internals.
Suppressed in `tsan.supp` (`race:__kmp_suspend_initialize_thread`, `race:__kmp_resume_template`); both names are
libomp-internal, so the entries cannot hide an AZURE2 race.

**On the CI (2026-11-15).** The runner's libomp (ubuntu-24.04, clang 18, `libomp.so.5`) has no symbols: its frames
print `<null>`, the two names match nothing, and every step of the thread job failed with this report (read from
the job's annotation: `pthread_mutex_lock` from libomp on the main thread against `pthread_mutex_init` from libomp
on a worker, the block allocated by libomp, AZURE2 only as the opener of the region in `ChannelRowFill.cpp`). A
third entry, `race:^pthread_mutex_init$`, names the interceptor of the initialization itself: AZURE2 never calls
it (its locks are `std::mutex`), so the entry cannot hide a race between AZURE2's accesses. It was not tried on
this machine (no clang here on that date); whether other reports stand behind it shows in the first CI run with it.

### 4.5 Memory not freed at exit (AZURE2; deliberate; not a leak that grows)

With `SAN_LEAKS=1`, LeakSanitizer reports on every CLI run exactly one direct leak, the `AZUREMain` that
`src/AZURE2.cpp:1740` allocates with `new` and never deletes; every other byte (up to 60 MB in 258,296 allocations
for the T2 fit) is reachable only from it (indirect leaks). The comment above that line says why: the teardown of the
engine objects took 2.3–3.7 s at the end of a `natural_neon` calculation, and no engine object holds a resource the
operating system does not reclaim. Not a finding; the runs keep `detect_leaks=0`. A memory growth during a fit is a
different question, measured by the peak RSS of `scripts/bench_upgrade.sh`.

### 4.6 `adaptive_grid_reference` fails without any sanitizer (pre-existing; not a sanitizer finding)

The case "eV level crossing the window edge" (Γ = 0.8 eV, 10 keV window, 50 points) gives a swing of the integration
error of 1.685 × 10⁻⁴ of the area as the level moves, against the tolerance 10⁻⁴ (`tests/reference/
adaptive_grid_reference_test.cpp`). The same number comes from a plain GCC 11.5 Release build of
`feature/upgrade`, from GCC 15.2 ASan and from clang 23.1.4 ASan. It is a
failure of the present engine's sub-point grid (`d216517d`, graded tails), not a memory or threading error. It needs
its own branch (`fix/adaptive-grid-window-edge`): either the grid near a window edge or the test's tolerance with a
written reason, decided by the author of `d216517d`'s test.

Resolved by decision 8 (a), commit 47c119e4 on `fix/adaptive-grid-criterion`: the swing was the grid's tail bias
leaving the window with the level (a one-way drift), so for the two crossing cases the test now bounds the variation of
the error beyond its net change, by √3/(2π ppw²) + 10⁻⁴ = 2.10 × 10⁻⁴ (the Euler–Maclaurin end term of the core
lattice at a fixed edge plus the in-window smoothness tolerance; derivation in the test). The grid is unchanged.

## 5. Running it

```bash
# clang with libomp + libarcher and the sanitizer runtimes installed (CI: clang-18, libclang-rt-18-dev, libomp-18-dev)
CC=clang CXX=clang++ scripts/sanitize_upgrade.sh thread            # build-san-thread/, about 50 min with the build on 8 cores
CC=clang CXX=clang++ scripts/sanitize_upgrade.sh address           # build-san-address/, about 25 min
SAN_SECTIONS="upgrade" SAN_NO_BUILD=1 scripts/sanitize_upgrade.sh thread   # T1-T3 only, existing build
SAN_REPORT_DIR=$PWD/san-reports scripts/sanitize_upgrade.sh address        # keep reports and logs
```

`scripts/sanitize_upgrade.sh --help` lists the variables. The build directory must not be shared between kinds.
