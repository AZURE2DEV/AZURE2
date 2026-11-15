# Fits that depend on the number of threads: causes, experiments, fix

A MIGRAD fit of `tests/upgrade/T2` with Minuit's numerical gradient, bounded to 200 function calls, ended at a total
χ² of 691.608 with 1 thread and 595.877 with 8, and the prior-covariance term of `chiSquared.out` differed in the
ninth digit between 1, 2 and 4 threads (0.306025731622686 / …729007276 / …733669139); at 8 threads the result was not
even repeatable (`HANDOFF.md`, 2026-11-09). Calculations without a fit (menu 1) were identical at every thread count.
This note records the two causes found, the experiment that shows each, the fix of the first, and what the second
needs. Thread counts are `OMP_NUM_THREADS`. "Verified" below means measured on this host with the builds named;
"inferred" means read from the code and not measured.

## 1. How the fit is run, and the diagnostic

The benchmark fit (`scripts/bench_upgrade.sh`, `work/fitrep.sh`) is menu 2 with `--use-lm --max-calls=200
--no-coulomb-table`. T2 has THM segments, which have no analytic Jacobian, so the Levenberg–Marquardt path declines
and MIGRAD runs with Minuit2's numerical gradient (`src/AZUREMain.cpp:463-475`). That gradient is an OpenMP loop over
the free parameters (`external/minuit2/src/Numerical2PGradientCalculator.cxx:116`), each thread calling the objective
`AZURECalc::operator()` at its own parameter vectors p. Each call takes a working copy of the model (`CNuc`, `EData`)
from a pool (`AZURECalc::GetPooledCNuc/EData`, `src/AZURECalc.cpp:2040-2095`), a stack of `max(4, cores)` clones made
once; which copy a call receives depends on the order in which threads return theirs.

`AZURE_TRACE_FCN=1` (added here, `src/AZURECalc.cpp`, after the objective is computed in `operator()`) prints every
evaluation as `[fcn] N chi2 | p...` with χ² and p in hexadecimal floating point. Two runs can then be compared value by
value: if every vector p evaluated in both runs has the same χ² bits, the objective is a function of p alone and any
divergence of the fits comes from the minimizer; if not, the objective carries state. `AZURE_TRACE_FCN=fresh` also
evaluates each fit vector on a fresh clone of the master objects and prints that value; a pooled copy must give the
same bits. Reproduction: the shortest case is T1 (`tests/upgrade/T1`, 13 free parameters, 3 s), 200 calls, 1 against
4 threads, for cause B; T2 is needed for cause A (T1 has no re-placed sub-point grids), where it appears within the
first 30 calls.

## 2. Cause A: re-placed sub-point grids kept a node set within 10⁻¹² of the present one

**The code.** At every evaluation `EData::ReplaceSubPointGrids` (`src/EData.cpp:1228`) places each adaptive
sub-point grid again from the current observed level parameters and the point's energy shift. A grid whose nodes all
lay within 10⁻¹² (relative) of the nodes it already held counted as unchanged, and the old nodes were kept
(`src/EData.cpp:1318-1321` before the fix). The grid a working copy held therefore depended on the parameter vectors
that copy had evaluated before: two copies at the same p could integrate on node sets that differ in their last bits,
and give χ² values that differ by about 10⁻¹⁵ relative. Minuit's parallel gradient distributes the vectors over the
copies by thread and by timing, so the result depended on the thread count and, through the timing, on the run.

**Experiments (verified, T2, builds of `fix/upgrade-fit-threads` with the trace).**

| run | result |
|---|---|
| 1 vs 2 threads, 200 calls, value at the same p | 22 of the 67 common vectors differ (first at evaluation 45 of the 1-thread fit), e.g. 10905.111934917488 vs …428 |
| 1 vs 4 threads, 200 calls | 34 of 59 common vectors differ (first at evaluation 25) |
| pooled vs fresh copy, 60 calls (`fresh`) | 60 of 79 evaluations differ at 1 thread, 6 of 79 at 4 threads; segments 1, 3, 4, 5, 8 (the target-integrated segments), first segment 4 (thick target, √E shift) |
| same, `--no-grid-replacement` | 0 of 79 differ at 4 threads |
| same, nodes compared exactly (the fix) | 0 of 79 differ at 4 threads |
| 8 threads, 3 runs, before the fix | prior-covariance term 0.454833808889256, …808793684, …808697643: three different results |
| 8 threads, 3 runs, after the fix | 0.45483380618971 three times, `param.sav` identical |

The same pooled copy evaluated twice at the same p gives the same bits; only a copy with another history differs.
Without grid re-placement the differences vanish, and with exact comparison they vanish too: the tolerance is the
state. This is the repeatable dependence on the thread count (symptom a) and, at 8 threads, where the order in which
copies return to the pool varies from run to run, a part of the non-repeatability (symptom b).

**Fix (commit on this branch).** A grid is kept only if its nodes equal the new nodes bit for bit
(`src/EData.cpp`, `ReplaceSubPointGrids`). The grid a copy holds is then a function of the current parameters alone,
as `AdaptiveIntegrationGrid::Quantize` (`src/AdaptiveIntegrationGrid.cpp:23-33`) already intended. Nodes that do not
move still keep their tables (`EPoint::AdoptEDependentValues` requires the exact energy, so nothing else changes).

**After the fix (verified, T2, 200 calls).** At 1, 2 and 4 threads the fits evaluate the same 216 vectors with the
same χ² bits, and `param.sav` and `chiSquared.out` are identical (prior-covariance term 0.306025733674053 at each).
The 1-thread result moves in its last digits, because the 1-thread history was itself one arbitrary choice of kept
grids: prior-covariance term 0.306025731622686 → 0.306025733674053, segment 2 χ² 60.0218 → 60.0219, total
691.608 unchanged.

## 3. Cause B: Minuit2's call counter loses counts under threads

**The code.** `MnUserFcn::operator()` increments `mutable int fNumCall` without synchronization
(`external/minuit2/src/MnUserFcn.cxx:21`, also `MnFcn.cxx:26`), from every thread of the numerical gradient. Concurrent
increments are lost, so Minuit counts fewer calls than were made. The count decides only when MIGRAD stops at
`--max-calls` (the iteration loop runs while the count is below the limit, `external/minuit2/src/VariableMetricBuilder.cxx:344`; the count also enters the call budgets of lines 105-181) and what
`fitStatus.out` reports as `calls`. This race is listed in `scripts/sanitizers/tsan.ignorelist`, where it was judged
not to touch a value the fit computes; it does change where a bounded fit ends.

**Experiments (verified).**

| run | evaluations made (trace) | calls Minuit reports | final total χ² |
|---|---|---|---|
| T1, 200 calls, 1 thread | 205 | 204 | 462.736 |
| T1, 200 calls, 4 threads | 234 | 220 | 462.477 |
| T2 after fix A, 1 thread | 216 | 215 | 691.608 |
| T2 after fix A, 8 threads, 3 runs | 256 each | 229, 230, 228 | 595.877 each |
| T2 after fix A, Minuit counter atomic (scratch build), 8 threads, 2 runs | 216 | 215 | 691.608 |
| T1, Minuit counter atomic, 4 threads | — | — | 462.736 |

In T1 every one of the 204 vectors the 1-thread fit evaluated is evaluated by the 4-thread fit with the same χ²
bits: the objective is not involved. At 4 threads Minuit lost 13 counts, was below 200 when the 1-thread fit stopped,
and ran one more iteration. In T2 at 8 threads it lost 26–28 counts and ran one more iteration (40 evaluations): that
iteration is the whole of the 691.608 / 595.877 difference. With an atomic increment (`#pragma omp atomic` before
both increments, a build in a scratch copy of the source, not committed) the 8-thread T2 fit and the 4-thread T1 fit
end exactly where the 1-thread fits end; every output agrees except `param.fit` (below). Whether the lost counts
cross the limit varies with timing, which is the rest of symptom (b): one of three 8-thread runs of the merged head
had ended at 691.608.

**Not fixed here.** The fix is two lines in `external/minuit2` (an atomic increment), and the root `CLAUDE.md` forbids
editing `external/` without the maintainer's decision. Until then: a fit that runs to convergence (no `--max-calls`)
is not affected except for the `calls` line of `fitStatus.out`; a fit bounded by `--max-calls` depends on the thread
count by up to one MIGRAD iteration. An alternative that leaves Minuit untouched would be AZURE2's own numerical
gradient (an `FCNGradientBase` with the same step rules), which is larger than this note's scope.

## 4. A third, cosmetic dependence: the `param.fit` snapshot

`operator()` writes the snapshot `param.fit` (and the intermediate outputs) at the evaluation numbered 100, 200, …
(`src/AZURECalc.cpp`, `thisIteration % kOutputInterval`). Evaluations are numbered by an atomic counter in the order
threads reach it, so inside a parallel gradient the vector that gets number 100 depends on the schedule. After both
fixes `param.fit` is the only output that still differs between thread counts (verified: T1 4 threads, T2 8 threads).
The final outputs overwrite the other intermediate files. A deterministic rule (the first evaluation numbered ≥ 100k
outside a parallel region) would change which vector `tests/park_formalism/check.sh` reads back with its logged χ²,
so it is left as an open item.

## 5. What was ruled out

- **Reductions whose order depends on the thread count.** The objective at a given p is bit-identical at 1, 2, 4 and
  8 threads once grids are exact (216 of 216 vectors, T2); calculations were already identical. Minuit's gradient loop
  writes only its own element per parameter (`Numerical2PGradientCalculator.cxx:126-163`).
- **The analytic gradient (`src/AZUREGrad.cpp`).** Not on the path: these fits use Minuit's numerical gradient.
  Fits with `--use-gradient` call the objective serially and parallelize inside it; not measured here.
- **Other pool state.** With exact grids, pooled and fresh copies agree at all 79 fit evaluations of a 60-call T2 fit
  at 4 threads (`AZURE_TRACE_FCN=fresh`). Other carried state (the energy-shift "unchanged" test, which compares
  exactly; `NodeTableCache`, keyed by exact energy and angle) did not show. Inferred, not measured: the shifted-point
  external-capture memo (`EPoint::GetECAmplitudeWithShift`, `src/EPoint.cpp:674`) records the final level's energy at
  a copy's first shifted evaluation, so a copy keeps the level energy of the first vector it saw; T1–T3 have no
  shifted capture segment, so no test here covers it.

## 6. Test and cost

`tests/fit_threads/check.sh` (run by `tests/run_tests.sh` and registered in ctest as `fit_threads`, about 1 min):
(1) a 30-call T2 fit at 2 threads with `AZURE_TRACE_FCN=fresh`, every pooled value equal to the fresh value;
(2) 30-call T2 fits at 1 and 4 threads, the same vectors with the same values (the number of evaluations is not
compared, because of cause B). With the old tolerance restored it fails both: 26 of 79 evaluations differ from a fresh
copy (first at evaluation 15), and 34 of 59 common vectors differ between 1 and 4 threads.

Cost (`scripts/bench_upgrade.sh`, this host, three interleaved runs each, `work/b-merged3` against this branch): fit200
at 1 thread 63.6 / 62.6 / 63.9 s before, 64.6 / 61.3 / 65.6 s after (means 63.4 and 63.8 s, +0.7 %, within the
spread); at 8 threads 34.2 / 34.5 / 34.8 s before, 29.1 / 33.5 / 34.9 s after; peak RSS 276–277 MB and 643–653 MB
both. The 29.1 s run is the cause-B case: Minuit lost too few counts to run the extra iteration, and that run ended at
691.608 instead of 595.877 (the other five 8-thread runs: 595.877). The calculation (menu 1) is unchanged:
1.0 s and 0.23 s.

Calculations are unchanged: all 179 projects under `tests/` in the three bases (menu 1, one thread, Coulomb table)
give identical outputs with `work/b-merged3` and with this branch (689 + 662 + 689 files, all identical). ctest
`-LE slow` 131/131 (`pyazr_park_dsde` needs `AZURE2_COULOMB_TABLE` set, as in the earlier verification runs); `tests/run_tests.sh` 134/134 at 8 threads (the 133 before and `fit_threads`).

Phase 6 (`ROADMAP.md`: evaluation without copies) removes the pools and with them the class of state that cause A
belonged to; the exact comparison is the containment until then.

## 7. The point pass inside a caller's parallel region (2026-11-15)

**What was seen.** A plain MIGRAD fit of `tests/13N` (`fit --max-calls=300`) took 13.6 s at 1 thread, 42.0 s at 2 and
25.7 s at 4 (34.0 / 26.8 / 25.7 s in three repeats); the Levenberg-Marquardt fit of T2 36.1 s at 1 thread, 84.2 s at 2
and 48.5 s at 4. More threads made these fits slower; the processor time at 2 threads was five times the time at 1
(13N: 74.3 s against 13.6 s).

**Cause (verified).** `WorkingCopy::CalculatePointsInParallel` (`src/EnergyPass.cpp`) began with
`if (omp_in_parallel()) return -1;`, and returned -1 always in a build without OpenMP. An evaluation called from
inside a parallel region (Minuit2's numerical gradient, the residual groups of the Levenberg-Marquardt fit) therefore
skipped the pass and calculated every point alone: one solve of the level matrix per point instead of one per compound
energy. With `OMP_NUM_THREADS=1` no region is ever parallel, so the fit at 1 thread kept the pass and was the fastest.
The sampler's likelihood never called the pass at all (`AZURECalcMCMC::LogLikelihoodOfFilled` evaluated point by
point), at any thread count.

**Repair.** The pass always runs. Inside a caller's region it runs on the calling thread alone (the three loops carry
`if (ownTeam)`, with `ownTeam = !omp_in_parallel()`); without OpenMP the same, serially. `AZURE_SERIAL_CHI2` still
turns it off. This is safe because every record the pass writes belongs to one evaluation: the `Results` (its scratch
slots are taken by `omp_get_thread_num()`, and an evaluation inside a region stays on one thread), the pooled
`WorkingCopy` and `LevelState`; the run's scratch (`RunScratch::ForThisThread`) is kept per thread under a mutex. The
sampler's likelihood now calls the pass as a fit's evaluation does.

**Measured** (this host, one run each unless said; elapsed seconds; `work/b-trial79` before, the branch's build after;
scripts and logs in `work/director/speed80/`):

| case | threads | before | after |
|---|---|---|---|
| 13N, plain MIGRAD, 300 calls | 1 | 13.6 | 13.9 |
| | 2 | 42.0 | 10.8 |
| | 4 | 25.7 (34.0, 26.8) | 8.2 (9.6, 9.3, 7.6) |
| T2, plain MIGRAD, 600 calls | 1 | 17.0 | 17.8 |
| | 2 | 16.2 | 10.3 |
| | 4 | 10.3 | 6.9 |
| T2, Levenberg-Marquardt | 1 | 36.1 | 36.4 |
| | 2 | 84.2 | 20.6 |
| | 4 | 48.5 | 16.2 |
| T2, sampler, 16 walkers, 100 steps, physical | 1 | 326.5 | 133.5 |
| | 4 | 166.0 | 48.4 |
| T2, sampler, the same, RWA | 4 | 195.4 | 73.4 |
| 13N, sampler, 16 walkers, 100 steps, physical | 4 | 187.8 | 30.6 |

The differences at 1 thread of the fits (+0.3 to +0.8 s) are single runs and not separated from the scatter. Peak
memory is unchanged (T2 sampler at 4 threads 296 MB before, 288 MB after).

**The numbers are the same.** Every output file of the thirteen pairs of runs above is identical before and after
(the chains `samples.mcmc` included, and the sampler at 1 against 4 threads), with one exception: `output/param.fit`
in some fits at more than one thread. That file is the snapshot of section 4, whose vector depends on the thread
schedule: it differs between two runs of one build as well (before: run 1 against run 3; after: run 1 against runs 2
and 3), so it is not a difference of the repair.

**The diagnostic `AZURE_NODE_CHECK=1`.** The check of the node table runs before each pass, so it now also runs on
the threads of a numerical gradient at once. Its lines were written piece by piece to standard error and ran into each
other (`tests/upgrade/node_check.sh`, the T2 fit at 2 threads: 154 whole check lines against 161 of each of the two
other kinds; no line reported a breach). The lines of one check are now collected and written at once
(`WorkingCopy::CheckNodeTable`); the test passes, with about twice as many checks as before in that case (83 before:
the gradient's evaluations had no pass and so no check).

**Not done.** A thread-sanitizer run (no clang on this host; the nightly job of the continuous integration covers T1
to T3 at 8 threads once the branch is there). No test fails when the pass is declined inside a region again (a
timing is not a test): the node check's line count in the T2 fit case would show it, and is not asserted. The fit
with the analytic gradient (`--use-gradient`) was already
outside a region and is unchanged (not timed again here).
