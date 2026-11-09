# AZURE2 upgrade: handoff

Update this file at the end of every work session: what changed, what is next, what is open. Newest first.

## 2026-11-09 (night) — fits that depend on the thread count: two causes (branch `fix/upgrade-fit-threads`)

Record: `fit_threads.md`. Diagnostic `AZURE_TRACE_FCN=1|fresh` (every χ² with its parameter vector in hex floats;
`fresh` also on a fresh model copy).
- **Cause A, fixed:** `EData::ReplaceSubPointGrids` kept a sub-point grid whose nodes lay within 10⁻¹² of the new
  ones, so a pooled copy's grid, and its χ² in the last bits, depended on the vectors it had evaluated before. Shown:
  pooled vs fresh copy differ at 60 of 79 T2 evaluations (none with `--no-grid-replacement` or exact comparison);
  1 vs 4 threads, 34 of 59 common vectors with different χ². Fix: exact node comparison. After it, T2 at 1, 2, 4
  threads evaluates the same 216 vectors with the same bits; 8 threads repeatable (3 runs). The 1-thread T2 fit moves
  in its last digits (prior-covariance term 0.306025731622686 → 0.306025733674053; total 691.608). Cost within noise
  (fit200 63.4 → 63.8 s at 1 thread). Test `tests/fit_threads/check.sh` (ctest `fit_threads`), fails with the
  tolerance restored.
- **Cause B, not fixed (needs `external/`):** Minuit2's call counter (`MnUserFcn.cxx:21`) is not atomic; at 4–8
  threads it loses 13–28 counts in 200, so a `--max-calls` fit runs one more MIGRAD iteration, timing-dependent. This
  is the whole 691.608 / 595.877 difference. An atomic increment (scratch build, not committed) makes T2 at 8 threads
  and T1 at 4 threads end where 1 thread ends. Decision for the maintainers: patch the vendored Minuit2 (two lines).
- **Open:** `param.fit` (snapshot at evaluation 100, 200, …) holds a schedule-dependent vector under the parallel
  gradient; the shifted-capture memo `EPoint::GetECAmplitudeWithShift` keeps the level energy of a copy's first
  shifted evaluation (inferred, not covered by T1–T3).

## Phase 3 merged; guardrails (written 2026-11-09, after the plan entry below)

**Merged into `feature/upgrade`:** `feature/upgrade-3-project`, **Phase 3** (design and record in `project.md`;
decisions P1–P5). 31 commits; pyazr −460 lines net (−980 / +517), `io/` +937 lines.
- `io/` (`azure2_io`, standalone like `physics/`): the lossless `io::Document`, the format conversion
  (`io::ConvertToFormat2`, format 2 in `<config>`), the engine's level grouping. Round trip byte for byte of every
  `.azr` under `tests/` and `examples/` in 945 variants (P1: the 48,922-file archive is on another machine).
- The engine reads every section from the document (`Config::SectionStream`); a format-1 file is converted in memory
  with a note; `AZURE2 --convert OLD NEW` and `pyazr.convert_azr` write it (P5: old files keep working).
- pyazr reads through the binding: document, segment / extrapolation / target-effect lines, level grouping, `<thm>`
  (`read_thm_block`; the 650-line Python port is deleted). The places where the two readers disagreed are decided
  as P2/P3 (`project.md` §5). pyazr needs the compiled module to read a project; 2.12.0.
- Step 4: `AzrModel` writes an unedited file back byte for byte in any layout or line ending; an edit changes only
  its lines (`azr_layout_test.py`).
- Step 5: one mass table. The engine's THM nuclei are generated from `pyazr/data/mass_1.mas20` (they were up to
  1 × 10⁻⁸ u off it; no output of the 179 projects moved) and checked by `pyazr_thm_nuclide_masses`; u, ħc, α, m_e
  have one definition (`physics/.../Constants.h`).
- Phase 2 leftover closed: the API's Park / abundance / scattering-length terms come from `Objective::ModelChi2`.

**Guardrails (`CLAUDE.md` "Recipes", `repo_consistency_test.py` 10–12):** a recipe per kind of addition (segment
keyword, `<config>` key or section, objective term, physics function, CLI flag), each with the check that fails when
it is bypassed: a constant copied with another value (10), a χ² term added outside `Objective.cpp` (11), the project
file opened outside `Config.cpp` (12). Each was shown to fail on an injected bypass. The fast pre-merge check is
`scripts/check.sh --quick` (consistency first, under a second).

**Verified** (branch head bba7f965 for the engine; the later commits are API, tests and docs): 179 × 3 — every output
identical except run.log, whose one changed line is the format-1 note; ctest 130/130 and the slow tests 10/10 on the
final head; `run_tests.sh` 133/133.

**Open, for the maintainers (they move results):** `EPoint.cpp`'s m_e is CODATA 2002 (0.510998903 MeV) in the
Fermi-function argument; `ECIntegral.cpp`'s proton mass 1.00727638 u (CODATA 2018: 1.007276466621). The first is a listed exception of check 10 (which covers u, ħc, α, m_e).

**Open, found by the benchmark (not caused by Phase 3): fits depend on the thread count.** Same binary, T2, MIGRAD
with 200 calls: the final χ² is 0.306025731622686 / …729007276 / …733669139 (prior-covariance term) at 1 / 2 / 4
threads, each reproducible; the writer's total ends at 691.608 (1 thread) and 595.877 (8 threads; one of three
runs of the merged head ended at 691.608, so at 8 threads not even repeatable), because the fit
is unconverged and the 9th-digit differences grow over 200 steps. Present on the Phase 1 build (b-pathways) and the
Phase 2 head too; the 2026-11-07 entry deferred it to Phase 2 ("fixed-order reductions in the numerical gradient"),
which did not address it. Calculations (menu 1) stay identical at 1 and 8 threads. Suspect: the pooled working
copies (`AZURECalc::GetPooledCNuc/EData`) carry state between evaluations, and Minuit2's parallel numerical gradient
hands each copy a thread-dependent subset of calls. Next step: evaluate one parameter vector on a fresh copy and on a
pool copy after a different history, and diff. Belongs with Phase 6 (no pools) unless a smaller cause is found first.

**Benchmark** (this host, three interleaved runs each, rows in `tests/upgrade/bench/T2.tsv`; before = 93f8cbc5, the
start of Phase 1): calc 0.93 → 0.90 s (1 thread), 0.42 → 0.24 s (8); fit200 63.2 → 56.2 s (1), 36.4 → 30.7 s (8);
fit peak RSS 241 → 276 MB (1), 488 → 643 MB (8; the interim shift tables). The pre-upgrade `dev` engine (aab1de0d)
cannot run T2 (χ² 1.2 × 10⁸: keywords it does not know), so it is not a baseline. The 93f8cbc5 build also ends the
8-thread fit at 595.877: the thread dependence above predates Phase 1.

**Document for the owner:** `docs/upgrade/FOR_THE_OWNER.md` (two pages).

Thread check on the merged head (calculations, menu 1): 12/12 identical at 1 and 8 threads.

**Next.** The fit thread dependence above (diagnose before new features); then the first feature through the
recipes, extrapolation segments on a user or adaptive grid.

## 2026-11-10 — plan to the PI's meeting with the owner (next week)

**PI (maintainer):** the goal is to hand the owner a code that is efficient, easy to maintain and easy to extend, while
the owner's AI keeps adding features. Agreed:
- finish Phase 3 (`project.md`; P1–P5 decided);
- the `CLAUDE.md` guardrails: recipes for adding a feature, consistency checks that fail when one is bypassed, a fast
  pre-merge check, a speed benchmark;
- a short document for the owner.

**PI, 2026-11-10: no merge with `origin/dev` during the reconstruction;** it is planned only when everything is
complete. The assessment below is kept for then.

**Requested feature (PI, 2026-11-10), after Phase 3, as the first feature through the new recipes:** extrapolation
segments on a non-uniform grid. Proposed design (format 2, a `<segmentsTest>` keyword):
- `energies <file>`, or `energies E1 E2 ...` inline, replaces min/max/step;
- `adaptive <tol>` bisects, from the coarse uniform grid, every interval where the curve departs from linear
  interpolation by more than `tol` (relative), up to a point budget, seeded with the level energies and widths;
- reuse `AdaptiveIntegrationGrid` where it fits.

**`origin/dev` since the upgrade branched (2026-11-08, 18 commits, measured 2026-11-10):**
- documentation (skills), 10 commits;
- portable as they are: target-effect auto windows (442610f3, 3a104e06), UDR time-of-flight kernels at low energy
  (a68474f0), test and CI fixes (5a08c048);
- to redo in the new structure:
  - the heavy-ion Coulomb threshold fallback (56512d6a): it patches `external/coul`, now `physics/coulomb`; to be
    checked against mpmath;
  - the coincident-level analytic limits (6d3be80a), now in `physics/level_parameters`;
  - pyazr NUTS sampling (bf09666a, 7de653af), which should use `objective_and_grad()`;
- obsolete: the Qt logo (b4be0500; decision 3) and the `libcoul` GSL link (63f0c093; library gone).

## 2026-11-09 (late) — Phase 2 (one objective) merged; band and bound channel functions merged

**Merged into `feature/upgrade`:**
- `shift-band`: ShftFunc interpolates only beyond η_b = 1000; one number moved, 1.2 × 10⁻¹¹ (`coulomb.md` §1.4).
- `bound-functions`: S(0) from its own Bessel K (2 × 10⁻¹⁵ against mpmath); GSL's central derivative ported bit for
  bit; outputs identical.
- `upgrade-2-objective`: **Phase 2**, design and record in `objective.md`. One ordered list of terms (`Objective`)
  for the FCN, the finite-difference and LM cost, the MCMC likelihood (data part) and the API.
  - One enumeration of the diagonal priors serves the gradients, G2, the least-squares rows and the MCMC; no
    parameter is matched by name any more.
  - pyazr computes no prior (`objective()`, `penalties()` and the new `objective_and_grad()` call the engine;
    pyazr 2.11.0).
  - The PI accepted **D1–D7** and decided that **the MCMC shares the fit's configuration**. As built:
    - D1: the MCMC profiles THM norms;
    - D2, D7: Park's wall and the abundance sum-to-one are least-squares rows;
    - D3: the prior covariance replaces its parameters' own priors;
    - D4: no norm prior for a profiled THM norm;
    - D5: a fixed parameter pays no prior;
    - D6: nuisance priors in the API, pyazr, the writer and the MCMC.
  - Every changed output and test is classified in `objective.md` §4. Of note, D1 corrects a diverging THM chain
    (level energy −3.2 × 10⁷ MeV before, 6.29 MeV after).

**Merged head 8f015075 verified as a whole:** 179 × 3 identical to the objective branch's head; ctest 132/132 with the
slow tests; `run_tests.sh` 133/133; 1 vs 8 threads identical.

**Mistake and fix this session:** a new pyazr test was written over the existing `tests/pyazr/objective_test.py`
(86ac2970); it is restored byte-identical (e8468cdc), and the new test is `objective_terms_test.py`.

**Next.** (1) Phase 2 leftovers (`objective.md` §4): the writer from the terms instead of stored totals; `inference/` as
its own library once Phase 4 makes the data values. (2) Phase 3 (`io/Project`: one `.azr` reader/writer). (3) Step 6b
with the observables. (4) Lecture 2.

## 2026-11-09 (night) — six branches merged; step 6a merged; band step under verification

**PI (2026-11-09): "go on with the steps and you can merge yourself, no need for my review."** Merged into
`feature/upgrade` (not pushed): `external-capture` (step 4), `zero-width-fixed`, `whittaker`, `brune-iteration`,
`whittaker-engine`, `coupling` (step 5), in that order. Conflicts: test lists in `tests/pyazr/CMakeLists.txt` (all
kept) and pyazr's version (2.10.0). **After this merge everyone rebuilds the Coulomb table (format 2).** The merged head
db972c84 was verified as a whole:
- 179 × 3 against the `whittaker-engine` run (format-2 table): the only changes are those the branches recorded
  (step 5's two rounding files, `screening_6li` formal, the `brune-iteration` formal run.log and 10⁻¹¹ changes);
- ctest 127/127 including the slow tests; `run_tests.sh` 133/133; 1 vs 8 threads identical.

**Step 6a** (`feature/upgrade-1-pathways`, merged 87721b0c; `cnuc_split.md` §10 and §6). Step 6 as planned would
rewrite about 730 accessor calls in 23 files, which is Phase 2's work, so it is split. 6a moves the pathway
enumeration (`SortPathways`) into `physics/Pathways` (`reaction_pathways` of the System and the external-capture
final states); the engine's trees are filled from it. Evidence:
- the trees match a verbatim copy of the old code on all 179 `.azr` in two bases (51524 nodes);
- 179 × 3 bit-identical; ctest 119/119, run_tests 133/133, threads identical.

`CNuc.cpp` is down to 2913 lines. 6b (the consumers of the trees) goes with Phase 2's observables.

**ShftFunc band** (`feature/upgrade-1-shift-band`, not merged; design in `coulomb.md` §1.4, "The interpolation band
near threshold"). Measured: `whittaker_scaled` is fine near threshold but scatters at 10⁻⁹ beyond η_b ≈ 10⁴. So the
η_b ≤ 1000 limit stays, and the GSL-only condition (η_b > 50 with z < 0.05) is dropped. B_min for p+¹²C goes from
3.3 × 10⁻⁴ to 8.2 × 10⁻⁷ MeV. Verification is running (179 × 3 with and without a table, ctest, threads, run_tests).

**Next.** (1) Band: classify, merge. (2) The rest of `CNuc`'s physics: `CalcShiftFunctions` (Park's numerical S′ and
S″, `g_shiftDerivMemo`) onto `channel_functions`. (3) Phase 2 (`inference/` objective) and 6b.

## 2026-11-09 (evening) — split `CNuc`, step 5: coupling coefficients of the observables

Branch `feature/upgrade-1-coupling` from 9be284f0 (step 4), not merged. Commits: design (`cnuc_split.md` §9,
da3703d2, alone), physics (96917f36), engine (563a7172), record (9143239d). `CalcAngularDists`,
`CalcCaptureAnalyzingPower` and `ThmAngular`'s Z̄ take their coefficients from `physics/ObservableCoupling` (no
GSL; new `wigner_9j`, `racah_w`); `AngCoeff::Wigner9j` deleted; `CNuc.cpp` 3135 → 3061 lines. `SortPathways`
moves to step 6 with the pathway trees.

**Evidence.** The physics tests check the 9j and, independently of the code under test, the particle and photon
coefficients against m-scheme angular distributions (237k checks, ≤ 1.4 × 10⁻¹⁴). The engine-side test checks
every coefficient against the old GSL expression (≤ 4 × 10⁻¹⁵ relative, same terms kept). 179 × 3: 2063/2069
identical; the six files left (`ec_e1m1_phase` s_1, `gamma_lab_frame` cm, in three bases) differ by at most
3.2 × 10⁻¹⁵ of the column maximum. The sign flips in `ec_e1m1_phase` are 10⁻¹⁹ residues of a cross section that
is exactly zero at 0° and 180°. ctest 113/113, run_tests 133/133, threads identical, timing unchanged. **Ready
for PI review/merge**; independent of the four branches waiting.

**Next.** (1) PI: merge order of the five branches (`zero-width-fixed`, `whittaker`, `brune-iteration`,
`whittaker-engine`, `coupling`). (2) Step 6: data accessors and the pathway trees (`SortPathways`) on the System.
(3) `ShftFunc` band; lecture 2.

## 2026-11-09 (day) — PI decisions; the Brune blocker fixed; step 3b (ii)

**PI decisions.** (1) `brune-iteration`: option (i), the iteration must converge where the root exists.
(2) afd8ee6a (`apply_fit_test.py` to §9's rule) confirmed. (3) The Whittaker engine branch waits for step 3b (ii)
— now done, see below.

| branch | state |
|---|---|
| `feature/upgrade-1-brune-iteration` | **verified, blocker fixed** (`level_parameters.md` §10, design a3bb203a + d6afeaeb + ff484603, code 5d8dc1e1, verified e4ca5d9d). E − e(E) is strictly increasing (Weyl, S′ ≥ 0), so each level has one Brune energy; the fixed-point map repels where e′ < −1 (−2.88 on the blocker: 2-cycle −0.973 / 4.462 MeV around the root 2.36820 MeV). A level that ends unconverged is now the root by bracket + bisection (status `Bracketed`). ctest 111/111, run_tests 133/133, threads identical; 179 × 3: only `screening_6li` formal changes (three background levels now transformed, equal to the Brune input digit for digit) |
| `feature/upgrade-1-whittaker-engine` | step 3b (ii) done (design cf5ddcce, code 82ca0985, verified b70ffdf9): final-state W memoized on its exact arguments, bit-identical (2040/2040), 99.8 % hits; T1 5.14 → 2.77 s, 13N 9.71 → 7.08 s (GSL build: 2.89 / 7.60 s). **Ready for PI review/merge**; after the merge everyone rebuilds the Coulomb table (format 2) |

**Next.** (1) PI: merge order (`zero-width-fixed`, `whittaker`, `brune-iteration`, then `whittaker-engine`).
(2) Narrowing `ShftFunc`'s band near threshold (own comparison). (3) Lecture 2.

## 2026-11-09 (night) — the three branches built and verified; the Whittaker engine part on a new branch; lecture 1

Work moved to a 24-core node (local builds, no Condor); at most 12 threads. Environment: conda env with Python 3.12,
numpy, scipy, mpmath, pybind11; build `cmake ... -DBUILD_GUI=OFF -DUSE_API=ON -Dpybind11_DIR=...` (34 s on 12 jobs).
The 179 × 3 comparison takes 51 s with a Coulomb table and ~8 min without; `run_tests.sh` ~20 min on 11 threads.

| branch | state |
|---|---|
| `feature/upgrade-1-external-capture` | `run_tests.sh` rerun: 133/133 (pending item closed); `ctest -LE slow` 111/111 |
| `feature/upgrade-1-zero-width-fixed` | **verified** (`external_capture.md` §9, "Verified"): 2040/2040 identical, ctest 112/112, run_tests 133/133, threads identical. Two test commits: check (c) of the new test called `apply_fit` in a way it refuses (fixed); `apply_fit_test.py` "fixed nonzero, no line: still refused" changed to §9's rule — **PI to review** (afd8ee6a) |
| `feature/upgrade-1-whittaker` | **verified** (`coulomb.md` §1.4, "Part 1 verified"): 2040/2040 identical, ctest 114/114 with the binding test |
| `feature/upgrade-1-brune-iteration` | verified except one test — **blocker for the PI** (`level_parameters.md` §9, "Blocker"): `pyazr_formal_capture_mask` (c) is NaN because `brune_from_formal` 2-cycles (E = −0.973 / 4.462 MeV) at ±10⁻⁶ MeV on `energy_3`; a18fe89e hid the same non-convergence by keeping formal values. Otherwise 179 × 3: Brune, Park identical, formal 657/662 (1.7 × 10⁻¹¹); run_tests 133/133 |
| `feature/upgrade-1-whittaker-engine` (new, from `-whittaker`) | step 3b (i) part 2: design addendum (a2593c96), engine switch (48fb454f), EC signature v6 with a stale-file case (5ed85322). **Changes outputs; PI to review before merging.** Numbers below |
| `docs/lectures-1` (new, from 7704e48f) | lecture 1 text written (`docs/lectures/01-external-region/lecture.md`), figures redone (no label boxes over curves, ranges fitted) |

**Whittaker engine part (`feature/upgrade-1-whittaker-engine`).** `WhitFunc`, `CoulombTable::ExactBound` and
`ThmBoundStateTail` call `whittaker` / `whittaker_scaled`; no engine file includes `gsl_sf_hyperg.h`. The Coulomb
table format is 2 (its bound sheet changed): a format-1 file is refused with a warning, so **everyone must rebuild
their table after the merge** (`AZURE2 --build-coulomb-table`, 112 s on 12 threads). Comparison without a table
against the `-whittaker` head: Brune 595/689, formal 559/662, Park 605/689 identical; every observable within
2.2 × 10⁻⁹ of its column's largest value; parameter files at most 9 × 10⁻⁸ relative (last printed digits), the
apparent large differences being values at rounding (10⁻¹⁸–10⁻³⁷ formal zero-input widths; 10⁻²² cross sections
at 0°/180° that vanish), one 0.99999999 eV width printed as 1000 meV, and −0.000000 signs; the 20 `intEC.dat`
caches change (hence v6). With tables built by each build: `c12c12` identical in the three bases. Note: the
format-1 table in `~/.azure2` on this node (26 Oct) moves `c12c12` by 1.9 × 10⁻⁵ against a format-1 table
built by today's reference build — the table can go stale within one format number. ctest 114/114 (with a format-2
table; `pyazr_park_dsde` needs a usable table, it fails when the only table is refused); threads identical;
run_tests 133/133. **Cost:** one-thread calc 13N 7.10 → 9.14 s, T1 2.62 → 4.73 s, all of it in
set-up (the external-capture table, 1176 integrals at 1320 nodes in T1, recomputed every run; evaluation unchanged,
~1 ms) — the 13× slower W per call. Step 3b (ii), one W table per integral, is the remedy; do it before merging.

**Next.** (1) PI: the brune-iteration blocker (safeguarded iteration, or move the test model); afd8ee6a; whether the
engine switch may merge with its cost or waits for 3b (ii). (2) Step 3b (ii): the per-integral W table in
`ECIntegral`'s data path, then the timing again. (3) Narrowing `ShftFunc`'s band near threshold (own comparison).
(4) Lecture 2 (channels, level and collision matrices, the three bases).

## 2026-11-08, end of session — lectures added to the plan; work stopped here

The PI asked for R-matrix lectures (Notion Markdown, figures, later slides): plan and lecture 1's figures in
`docs/lectures/` (`PLAN.md`). Lecture 1's text is not written. Everything else is as in the entry below: four
branches, the three on top of step 4 BUILD-PENDING (Condor schedd on crcfe02 down).

## 2026-11-08, later — three branches on top of `feature/upgrade-1-external-capture`; builds pending

The PI accepted step 4's fixes and its two proposed defaults (2026-11-08). Each is a branch cut from
`feature/upgrade-1-external-capture` a18fe89e, design note first, then code; they do not depend on each other.
Merge order: this branch, then the three in any order. **Not built yet**: the login node's Condor schedd has been
down since the node rebooted (~07:40); the code is syntax-checked (`g++ -fsyntax-only` with the build's flags, no
warnings in the touched lines) and each code commit says "BUILD AND TESTS PENDING".

- `feature/upgrade-1-brune-iteration` (worktree `azure2-up-ec`): `brune_from_formal` inside an evaluation converges
  at max(10⁻¹² MeV, 10⁻¹²|E|), takes its last iterate at the rounding floor (≤ 10⁻⁹ MeV, status `RoundingLimited`),
  and gives NaN instead of the formal values otherwise; setup failures stop the run (`level_parameters.md` §9).
  Tests: `level_parameters_unit_test` (new cases), `pyazr_formal_brune_iteration` (slow; χ² of
  `beam_profile_kernel` smooth in `width_9_2`, which jumped 40.44 → 42.93).
- `feature/upgrade-1-zero-width-fixed` (worktree `azure2-up-zero`): a width given as 0 is fixed in every basis
  (`ZeroInputWidth`, not the transformed value); `apply_fit` skips fixed unlisted channels; pyazr 2.9.1
  (`external_capture.md` §9). Test `pyazr_zero_width_fixed`. Expected: formal free set = Brune free set (134 → 67 on
  `beam_profile_kernel`), `--no-ignore-externals` input-0 photon widths fixed.
- `feature/upgrade-1-whittaker` (worktree `azure2-up-whit`), step 3b (i) part 1: `physics/coulomb` `whittaker`,
  `whittaker_scaled` without GSL (integral representation of U, double-exponential quadrature; `coulomb.md` §1.4).
  Against mpmath on 1162 points: W 2.1 × 10⁻¹², ρW′/W 3.3 × 10⁻¹³ (GSL as the engine calls it: 6.9 × 10⁻⁹, 10 points
  without a value); COULCC at imaginary argument dropped (0/0 at η_b = l+1, l+2, …). 4.8 µs per call against GSL's
  0.36 µs. `physics_whittaker` and `whittaker_gsl` were run standalone (pass; four mutations each fail it); the
  binding `_azure2.whittaker` and `pyazr.whittaker_wave` (pyazr 2.10.0) wait for the build.

**To do when Condor is back** (per branch): build (`upgrade-work/build.sh SRC BUILD OFF`), `ctest -LE slow`, the new
slow test, `run_tests.sh`, 1 vs 8 threads, the 179-project comparison against a18fe89e (`ec-work/runall.sh`,
`classify.py`), each change classified; then a follow-up commit with the numbers. Also still pending:
`run_tests.sh` on this branch's 709ed040 build (`ec-work/verify.sh`).

**Next (step 3b (i) part 2, design in `coulomb.md` §1.4):** `WhitFunc` over `whittaker`/`whittaker_scaled`,
`ShftFunc`'s bound branch on ρW′/W (whether the interpolated band near threshold narrows is a separate comparison),
`CoulombTable`, `ThmNumerics`, `ThmFunc`; measure `ECIntegral`'s data path (13× per W call) before and after, and
route `pyazr.whittaker` (transform.py) through the engine. Bound-state outputs move in their last digits; classify.

## 2026-11-08 (`feature/upgrade-1-external-capture`) — split `CNuc`, step 4: `physics/external_capture`; the mask on the input

Branch `feature/upgrade-1-external-capture` from `feature/upgrade` 5055bec0, not merged. Commits: the design
(`external_capture.md`, alone), the physics part, the engine through the adapter, the mask on the input, the
lowest-level final state with a free energy (step 2 open item 2); then the addendum (section 7, alone), one adjoint
test for every gradient path with decoupled masked levels on P, `tests/pyazr/formal_capture_mask_test.py`, these docs.

**State.** `physics/external_capture`: N_f, the multipole rule, γ^ext (Eq. 4, radial integral injected), the
projector Π = T⁻¹PT; the adapter `src/CNucExternalCapture.cpp` replaces `CalcExternalWidth` and six copies of
N_f^−½ (`CNuc.cpp` 3349 → 3132 lines). `--ignore-externals` masks a photon channel whose input width is exactly zero
(`ALevel::ZeroInputWidth`), the same in every basis and continuous in a fit; the formal basis carries the Brune mask
with Π. The adjoint does not differentiate Π or the step-2 final-state amplitudes: `FormalPartsOutsideAdjoint`
(`src/AZUREGrad.cpp`) sends the χ² gradient, the residual Jacobian and the model gradients to finite differences
when a slot has a projector or a masked level has a free particle amplitude, or a formal final-state block exists;
a masked level with no particle width is decoupled and uses P (exact). How the parts stay apart: the physics part
gets arrays and returns arrays; the adapter only reads `CNuc` and calls it; the gradient asks `CNuc` one question.

**Evidence** (`external_capture.md` §6). 179 `.azr` × 3 bases against 5055bec0: Brune 689/689, Park 689/689
identical; formal 650/662, every change classified (`upgrade/T1` capture and χ², `beam_profile_kernel` capture,
χ² and three input-0 E2 widths), each changed formal output now equal to the Brune one. The guard commit:
2040/2040 files byte-identical to 3d08b7f0. χ² gradient vs central differences, formal basis, 13N with a masked
coupled 1/2⁺ level: 1.7 × 10⁻³ → 1.4 × 10⁻⁷ (new test; fails on e205038e; other checks fail on 5055bec0). No
warnings; check.sh PASS; `ctest -LE slow` 111/111 (`pyazr_zero_width_level` passes again); `run_tests.sh` 133/133 on 3d08b7f0 (on the final build cut short by the 2026-11-08 node reboot, rerun pending; its
menu-1 outputs are byte-identical to 3d08b7f0's);
1 vs 8 threads byte-identical (T1–T3, 13N, three bases); T2 bench (Brune; mean of three, same host, 5055bec0 → 3d08b7f0; the
guard is a no-op in Brune) fit200 61.58 → 62.16 s (1 thread, +0.9 %), 26.61 → 26.65 s (8); calc 1.03 → 0.98 s,
0.21 → 0.20 s.

**Open (proposed default in brackets).**
- **Blocker for merging the formal basis as "continuous":** `brune_from_formal` at 10⁻¹² MeV is at the rounding
  limit (§8). On `beam_profile_kernel` a 10⁻⁶ change of one width makes the 20 MeV 1⁻ level miss the tolerance
  (last step 7.7 × 10⁻¹² MeV after 1001 iterations); the formal values are then kept silently, its row of T becomes
  a unit vector and χ² jumps 40.44 → 42.93. The merged step-2 final-state amplitudes use the same call. [A
  tolerance relative to the rounding of E, e.g. max(10⁻¹² MeV, 64 ε |E|), in both callers; a level that still does
  not converge is an error of the evaluation, not formal values. Needs the PI: it changes merged formal-basis
  numbers at ~10⁻¹¹ and the no-silent-fallback rule decides the second half.]
- Width parameters with value 0 are fixed at setup (`CNuc.cpp:2259`), tested on the transformed value, so in the
  formal basis rounding decides: `beam_profile_kernel` frees 35 parameters with |value| < 10⁻¹² on 5055bec0 and 48
  here (10 new on the lowest 1⁻, 2⁺, 3⁻ levels, from 8df95d22's reordered sums; 3 masked photon channels now exactly
  γ^int = 0 plus rounding). [Fix on the input, as the mask: a width given as 0 is fixed in every basis. Changes which
  parameters a formal-basis fit varies; PI decision.]
- Formal basis with a projector or a final-state block: finite-difference gradients. [Accept; the derivative of T
  is a Phase 2 item.]
- From step 2, unchanged: `TransformOut`'s 10⁻⁶ MeV report iteration; `ThmLevelWidth`'s own total width.

**Next.** PI decisions on the two items above; then merge; then step 3b (i) (Whittaker functions without GSL).

## 2026-11-08 (`feature/upgrade-1-cnuc-split-2`) — split `CNuc`, step 2: observed quantities; formal-basis external capture

Branch `feature/upgrade-1-cnuc-split-2` from `feature/upgrade` 9d50667f, not merged. Commits: the design
(`cnuc_split.md` §8, alone), `observed_levels` and its callers, the formal-basis external-capture fix with its test,
the restriction to final states above the lowest level, these docs.

**State.** `level_parameters::observed_levels` is the one definition of the observed quantities from total amplitudes
γ^int + γ^ext: W = s 2P|γ^tot|²/N and the `<levels>` value (Γ eV, C fm^−½, μ, Q, B, or the Brune amplitude of a
`gammaIsRWA` column) with its kind. `TransformOut` = `TransformToBrune` + `ObservedWidths`; `GetTransformParams`,
`PrintTransformParams`, `CheckRadiativeWidths` read the stored values (the last now with N and the external part);
`ConvertAmplitudeBasis` calls `park_from_brune`/`brune_from_park`. The five γ → Γ copies are deleted. External capture
takes the final state's Brune amplitude and energy in every basis (`ALevel::ECGamma`, `ECEnergy`; formal basis: a
read-only setup record of the blocks with a final state above their lowest level, `brune_from_formal` to 10⁻¹² MeV in
`FillCompoundFromParams`; γ^ext subtracted before any group is made formal); the intEC signature hashes the energy
used. How the parts stay apart: `level_parameters` gets arrays (E, γ, γ^ext, flags, channel functions) and returns
arrays; the engine only gathers and stores; external capture still supplies γ^ext from `CNuc` (step 4).

**Evidence** (`cnuc_split.md` §6, step 2). 179 `.azr` × 3 bases against 9d50667f: Brune 689/689, Park 689/689
identical; formal 649/662, every change classified (8 `parameters.out` μ/Q instead of G = 0, 3 `screening_6li` logs
lose a warning built on formal amplitudes, `beam_profile_kernel` `param.par` and `parameters.out` whose photon widths
to the 6.049 MeV state now reproduce the input, e.g. 0.000053 → 0.000330 meV for 0.33). No re-baseline. Formal ↔ Brune
capture on `beam_profile_kernel` (111 energies, `--no-ignore-externals`): 6.049 MeV state 0.808 → 1.5 × 10⁻¹¹; the other
four final states 0 before and after. New `tests/beam_profile_kernel/check.sh` (fails on 9d50667f, 0.448). No warnings;
`check.sh` PASS; `ctest -LE slow` 109/109; `run_tests.sh` 133/133; 1 vs 8 threads byte-identical (T1–T3, three bases);
T2 bench (mean of three, same host) fit200 133.1 → 131.8 s (1 thread), 55.8 → 56.2 s (8, +0.7 %); calc 2.45 → 2.22 s,
0.42 → 0.44 s.

**Open (proposed default in brackets).**
- `--ignore-externals` (the default) drops the external resonant capture of levels whose internal photon amplitude is
  < 10⁻⁸ in the basis in use, so formal and Brune still differ by default (9.7 × 10⁻² on the `beam_profile_kernel`
  data). [Step 4: define the mask on the Brune amplitudes in every basis.]
- A formal-basis fit that moves the energy of a lowest-level final state: its Brune amplitude then differs from the
  formal one at second order and is not recomputed. [Accept until step 4; or record every block with a final state
  whose energy is free.]
- Formal basis with a final state above the lowest level: the analytic gradient falls back to finite differences.
  [Accept; the derivative of `brune_from_formal` is a Phase 2 item.]
- The reports (`parameters.out`) of the formal basis keep `TransformOut`'s 10⁻⁶ MeV iteration; deeply bound channels
  amplify it (an ANC input 0 reported as 38.65 fm^−½ in `beam_profile_kernel`, as before). [Tighten to 10⁻¹⁰ MeV in a
  later step, with its own comparison.]
- `ThmLevelWidth` (THM line shape) has its own total-width formula. [Phase 5 with THM.]

**Next.** Merge; then step 4 (`physics/external_capture`) or step 3b (i).

## 2026-11-07 (director) — step 3b merged (4a8ecf02)

Per-stage channel table and the Coulomb memo's deletion are on `feature/upgrade`. Accepted as interim, with the
reason: 13N +1–9 % per pair (inside the base's spread; revisit in step 4, `physics/external_capture`); T2 fit peak
memory +30 % at 8 threads from the three shift-state tables per segment (gone with Phase 2's analytic shift
derivative); the χ² of a fit differs in the 8th digit between 1 and 8 threads, already so on the base (Phase 2:
fixed-order reductions in the numerical gradient). Next: step 2 (observed widths from the total amplitude), step 4
(external capture), then the Whittaker functions without GSL and `--coulomb-method`.

## 2026-11-07 (`feature/upgrade-1-channel-table`) — step 3b: the per-stage channel table; memo deleted

Branch `feature/upgrade-1-channel-table` from `feature/upgrade` 815d9b5f, not merged. Commits: the design
(`channel_table.md`, alone, with the `perf` measurement of where `-nomemo` spent its time), the table, the radial
scratch, the `-nomemo` commits re-cut (c3dc22e5, 37a3637b), the entrance phase through the scratch, keys unique per
point and the table released after setup, the first docs; then (director's answers) the free-shift tables, the
capture amplitude from the table, distinct recorded rows, these docs.

**State.** `ChannelTable` (src/): exact keys (pair α, l, E_c), rows P, S, e^{−iφ}, e^{i(σ_l − σ_0)} from one function
(`ChannelTable::Compute`), built by `EData` before the point loops (keys from `EPoint::ChannelKeys`), one thread per
row, read-only, released after setup; at setup the capture amplitude reads √P and the entrance phase at a_c from it.
Free-shift refills: each segment keeps the tables of its last three shift states (key: the shift values bitwise;
`BeginShiftRefill`/`CommitShiftTable` single-threaded around the parallel loop; replaced, never edited) — interim,
deleted when the objective has the analytic shift derivative (Phase 2). `CoulFunc::SetRadialScratch`: thread-private
exact-argument scratch for the capture integrals (accepted under RULES §2.3, `channel_table.md` §3).
`CoulFuncCache` deleted. Stay: the per-slot rows of `EDependentTables` (Phase 4), `g_shiftDerivMemo` (until
`level_parameters` computes S, S′, S″ at E_λ).

**Evidence** (`channel_table.md` §6). 179 `.azr` × 3 bases against the memo base: observables ≤ 7.0 × 10⁻¹¹ (one
exact zero's residue), `intEC.dat` ≤ 4.0 × 10⁻¹¹; against `-nomemo` 670/689 byte-identical, 19 `intEC.dat`
≤ 2.8 × 10⁻¹³ (one ulp of E_c). T2 fit200 (1 thread) byte-identical to `-nomemo`. No re-baseline. No warnings;
`channel_table_unit`; check.sh PASS; `ctest -LE slow` 109/109; `run_tests.sh` 132/132; THM
byte checks 5/5; 1 vs 8 threads byte-identical on T1, T2, T3, 13N (menu 1); TSan (upgrade section) 0 reports (first
build).

**Speed** (qa-xp-003, alternating, three repeats; base → branch): T2 calc 2.29 → 2.30 s (1 thread), 0.68 → 0.37 s
(8); fit200 140.1 → 121.7 s (1), 68.8 → 53.3 s (8); 13N 15.52 → 15.76 s; natural_neon 21.48 → 18.66 s; `tgt_fixed`
57.65 → 57.41 s. 13N is slower in each of five pairs (+1 to +9 %) but inside the base's spread. Peak RSS of a T2 fit
+13 % (1 thread), +30 % (8): the shift tables.

**Open (proposed default in brackets).**
- 13N +1–9 % per pair. [Accept; the rest goes with step 4 (`physics/external_capture`, unit amplitudes per grid
  energy).]
- Fit peak RSS +30 % at 8 threads (three shift tables per segment and working copy). [Accept until Phase 2 removes
  them; or keep two states, which loses the return to ν₀ after ν₀ ± h.]
- A fit at 1 and 8 threads differs in χ²'s 8th digit — on the memo base and `-nomemo` too. [Separate item for the
  Phase 2 objective (fixed-order reductions in the numerical gradient).]

**Next.** Merge, then step 3b (i), (iii), (iv) or step 2.

## 2026-11-07 (`feature/upgrade-1-coulomb`) — split `CNuc`, step 3: Coulomb functions from pycoulomb (memo kept)

Branch `feature/upgrade-1-coulomb` from `feature/upgrade` f549fa23, not merged. Commits: the design (`coulomb.md`,
alone), `physics/coulomb`, the GSL cross-check, the engine routed and `external/coul` deleted, `channel_functions`
in physics/, the engine calling it, these docs. The memo's deletion (c3dc22e5, with its docs) is on
`feature/upgrade-1-coulomb-nomemo`, to be re-cut on the per-stage table (director's decision; `cnuc_split.md` §3,
step 3b ii′).

**State.** `physics/coulomb`: `coulomb(l, η, ρ, method)` → F, F′, G, G′, σ_l, pure, NaN with a status on failure;
methods coulsm (default; AZURE2's former `external/coul`), coul, coulfg, coulcc, one adaptor file each, sources
byte for byte from pycoulomb azure2-prereqs 7adef98 and coulcc-cpp 37f0f9b (`physics/coulomb/README.md`).
`physics/ChannelFunctions`: P, S, e^{−iφ}, S′ (Richardson), σ_l of an open channel and η, ρ (non-relativistic and
relativistic) with the constants; the engine's `Kinematics`, `Constants.h` and `CoulFunc` call them, and every
Coulomb evaluation goes through `CoulombWavesAt` → `coulomb(…, coulsm)`, behind the unchanged `CoulFuncCache`.
Deleted: `external/coul` (8 files, 4053 lines). Not done (step 3b): Whittaker functions without GSL (`WhitFunc`,
`ShftFunc`, `ChannelFunc`'s bound branch stay), the per-stage table, then the memo's deletion; hybrid channels,
`CoulombTable`, `--gsl-coul`.

**Evidence** (`cnuc_split.md` §6, step 3). 179 `.azr` × 3 bases. Routing: identical except the hybrid-potential
projects (cwfcomp's 10⁻¹⁵ imaginary parts no longer enter the hybrid H⁺): observables ≤ 9.9 × 10⁻¹¹, one formal
width 5.0 × 10⁻⁸. `channel_functions`: bit-identical (2040 files). No re-baseline. physics_coulomb (mpmath, 1296
points): coulsm/coul ≤ 8.3 × 10⁻¹⁰ where the functions oscillate, coulfg 6.3 × 10⁻¹⁴, coulcc 1.5 × 10⁻¹²
everywhere; identity 8.3 × 10⁻¹⁰; σ_l 7.1 × 10⁻¹⁵; 8 threads bit-identical. GSL cross-check 4.6 × 10⁻¹⁰ (448
points). physics_channel_functions: neutral closed forms 6.0 × 10⁻¹⁰, P/S/φ against mpmath ≤ 7.5 × 10⁻⁹, S′ against
coulcc 1.8 × 10⁻⁷.
Checks (GCC 11.5, GUI off, Condor, build of 08516cfb): no compiler warnings; `tests/upgrade/check.sh` PASS (goldens
unchanged); `ctest -LE slow` 110/110 — `pyazr_thm_profiled` and `pyazr_thm_fixed_angle` included, and both also pass
on f549fa23 in 2 of 2 repeats on the same host, so their failure in the director's run is intermittent (the memo is
per thread and call-order dependent); `tests/run_tests.sh` 131/132, the failure the known `shared_point_tables` (d)
peak-RSS ratio (316 MB > 3 × 105 MB). 1 vs 8 threads: T2 and T3 outputs byte-identical; T1 and 13N identical except
`intEC.dat` (≤ 3.4 × 10⁻¹⁴ and ≤ 5.1 × 10⁻¹² of the line's largest value): the memo's per-thread caches, gone on
`-nomemo`, where all four were byte-identical.

**Speed** (T2.tsv rows; `coulomb.md` §7). host d32cepyc086, mean of two, before (f549fa23) → this branch: T2 calc 1.56 → 1.36 s (1 thread; one
f549fa23 run at 1.89 s), 0.44 → 0.42 s (8); T2 fit200 80.5 → 81.9 s (1 thread), 36.5 → 37.0 s (8); menu-1
calculation, one thread, `13N` 10.05 → 10.16 s, `natural_neon` 12.35 → 12.45 s; `convolution_unify` `tgt_fixed`
(2 threads) 29.9 → 30.6 s. The routing costs at most 2 % (fit and `tgt_fixed`), within the run-to-run spread of the
calculations; peak RSS within 2 %.

**Measured on `-nomemo`** (`coulomb.md` §7): outputs ≤ 7.0 × 10⁻¹¹ (one exact zero's residue), the two THM checks
pass (their failure is the memo's call-order dependence), but 3–17× slower (T2 calc 1.37 → 6.16 s, `natural_neon`
12.3 → 105.6 s) and `convolution_unify` `tgt_fixed` 84 → 1414 s, past its 600 s limit.

**Open (proposed default in brackets).**
- The THM identity checks and the 1-vs-8-thread `intEC.dat` differences are intermittent memo effects. [Accept until
  `-nomemo` is re-cut on the table; there both are deterministic by construction.]
- Whittaker functions without GSL. [COULCC at imaginary argument in `physics/coulomb`, against mpmath and GSL; bound
  outputs will move in their last digits and are classified like this step.]
- `--gsl-coul`, `CoulombTable`. [Add `--coulomb-method=coulsm|coul|coulfg|coulcc` and retire `--gsl-coul` with a
  message; keep `CoulombTable` until the per-stage table is measured against it.]
- coulfg's BSD-3 text (the LLNL COPYRIGHT file) is not in pycoulomb. [Jakub adds it in pycoulomb; copied here then.]

**Next.** Step 3b (ii), the per-stage table, then `-nomemo` re-cut on it; or step 2, which does not depend on it.
## 2026-11-07 (director, end of session) — where to continue

**Branches (local, bundled in `azure2-upgrade.bundle`; pycoulomb branches in `pycoulomb.bundle`).**
- `feature/upgrade`: Phases 0 and 1 so far — test projects T1–T3, goldens, sanitizers, `physics/` (System, JpiBlock,
  Convention, level_parameters), CNuc split steps 1 and 1b, Qt removed, decisions 1–12 recorded below.
- `feature/upgrade-1-coulomb` (not merged): pycoulomb's C++ in `physics/coulomb/` (coulsm default, coul, coulfg,
  coulcc), `channel_functions`, engine routed through them, `external/coul` deleted; output-identical except the
  hybrid-potential rounding (≤ 9.9 × 10⁻¹¹). Its docs/timing commit may be missing if the session ended first: check
  `git log`; if the head is 08516cfb, add the docs (coulomb.md describes the memo as still present) and time T2,
  13N, natural_neon, convolution_unify against f549fa23 before merging.
- `feature/upgrade-1-coulomb-nomemo` (hold): deletes `CoulFuncCache`; correct (≤ 7 × 10⁻¹¹, fixes the two pyazr THM
  byte checks) but 3–17× slower until the per-stage table exists. Re-cut on top of step 3b.
- pycoulomb (Jakub's repo): `azure2-prereqs` (libcoulomb, tests, LICENSING_NOTES.md) and `coulcc-cpp` (COULCC in
  C++, bit-identical to the Fortran). Jakub: add the LLNL BSD-3 COPYRIGHT for coulfg, choose the licence, push, tag.

**Next steps, in order.**
1. ~~Merge `feature/upgrade-1-coulomb`~~ Done (a9cd42c3): routing costs ≤ 2 % (within spread); outputs identical
   except hybrid rounding; `ctest -LE slow` 110/110. The two pyazr THM byte checks are intermittent (the per-thread
   memo's call order), not a merge defect; they become deterministic with step 3b.
2. Step 3b (`cnuc_split.md` §3, `coulomb.md` §6): the per-stage channel-function table — P, S, S′, φ, σ_l per
   (α, l, channel energy) over the evaluation energies (decision 11, exact energies), read-only during an evaluation,
   replaced when a free shift moves; replaces the `EDependentTables` Coulomb rows, `g_shiftDerivMemo` and the memo.
   Gate: outputs ≤ 10⁻⁹ or explained; T2 and `convolution_unify` faster than with the memo. Then re-cut and merge
   `-nomemo`.
3. Whittaker functions without GSL (COULCC at imaginary argument); `--coulomb-method`; retire `--gsl-coul`.
4. Remaining CNuc split steps (`cnuc_split.md` §3), then Phase 2.

**Known open items.** `shared_point_tables` (d) peak-RSS limit has no margin; `light_excitation` (a) external-capture
integrals depend on thread scheduling at 5 × 10⁻¹¹; formal basis with external capture disagrees with Brune (step-1b
entry, item 2); Sphinx never built on CRC (`sphinx-build -W docs/source` before pushing). Agent brief:
`~/AI-R-scratch/upgrade-work/AGENT_BRIEF.md` on CRC (rules are also in CLAUDE.md).

## 2026-11-07 (director) — state after the step-1b and Qt-removal merges: two pyazr THM checks fail

`feature/upgrade` 160439ba, fresh build, GCC 11.5, no Qt (Condor 16205, 16206): no warnings; `tests/upgrade/check.sh`
ok; `tests/run_tests.sh` 131/132, the failure the known `shared_point_tables` (d) peak-RSS ratio (305 MB > 3 × 90 MB;
the limit has no margin). `ctest -LE slow` 105/107: `pyazr_thm_profiled` and `pyazr_thm_fixed_angle` each fail one
check, "write_output_files: AZUREOut_aa=5_R=4.out identical to the CLI's" (χ² agrees to 10⁻⁹; the other output
files are identical). Both passed on each parent of the merges. Under diagnosis on `fix/thm-pyazr-cli-identity`;
the Coulomb-function step waits for it. Decisions 11 (exact energies) and 12 (data as model expressions) added.

## 2026-11-07 (`feature/upgrade-1-cnuc-split-1b`) — split `CNuc`, step 1b: decision 10 on, one solver, one channel energy

Branch `feature/upgrade-1-cnuc-split-1b` from `feature/upgrade` 4bf7aa06, not merged. Commits: the design
(`cnuc_split.md` §7, alone), decision 10 with the options deleted, the channel energy and adapter, the sign rule
corrected, a strongly mixed test, the docs, the benchmark rows, this entry.

**State.** The formal basis has two rules and no options: B_c = S_c at the lowest level of each J^π group, and each
formal level oriented like the Brune level it comes from (`level_parameters.md` §3, §4; `conventions.rst`).
`level_parameters` uses its own eigensolvers only; `BoundaryRule`, `FormalSignRule`, the solver options and
`GslGeneralizedEigen`/`GslSymmetricEigen` are deleted. The adapter forms E_c = (E_x − S) − e, as at a data point, and
returns S, P, S′ at every call (`ChannelEnergyOrder`, `ChannelQuantities` deleted). In the formal basis the engine
prints a note for each group whose first listed level is not its lowest.

**Evidence** (`cnuc_split.md` §6, step 1b). 179 `.azr` × 3 bases, 4bf7aa06 vs 7f4af9d7: Brune 670/689 files
identical, Park 685/689, formal 595/662. Class (a), formal basis: 21 `param.par` (sign, B_c in `natural_neon` and
`be7p_fit`), 25 `parameters.out`, observables ≤ 8.5 × 10⁻¹¹. Class (b), last digits from the E_c order: Brune
observables ≤ 9.8 × 10⁻¹¹ printed (3.6 × 10⁻¹¹ full precision), Park ≤ 6.7 × 10⁻¹¹; Brune amplitudes up to 6.3 × 10⁻⁹ in
the ⁷Li + p models' 0⁺ background level (negative Brune denominator). No class (c). Goldens and `expected/` unchanged
(T2 fit within its own run-to-run spread); no re-baseline. GUI off, Condor: no warnings; `check.sh` PASS;
`ctest -LE slow` 107/107; `run_tests.sh` 131/132 (`shared_point_tables` (d), the flaky peak-RSS ratio of the step-1
entry: the calculation's RSS varies 67–105 MB between runs; fit RSS 301–308 MB here vs 310–312 MB before; 3 repeats:
1 failure here, 0 before). T2 bench, same host, mean of two: fit200 78.54 vs 78.53 s
(1 thread), 34.30 vs 34.36 s (8); calc 1.16 vs 1.21 s, 0.39 vs 0.39 s.

**Found.** (1) Decision 10's "largest |γ| positive" changes observables: it reversed the 6.917 MeV 2⁺ external-capture
final state of `beam_profile_kernel` (largest formal amplitude a photon one) and the formal-basis capture to it changed
by up to 4.51 relative. A final state's sign is shared with the photon amplitudes that end on it and with the
external-capture amplitude, so it is not free. Built instead: the orientation rule, also solver-independent. (2) The
formal and Brune bases disagree in capture with external capture, before and after: 1.5 × 10⁻² to 9.1 relative on
`beam_profile_kernel`'s five final states (the cross section takes the final state's formal amplitude, Γ_γ its Brune
one; with the group order of `TransformIn`, `level_parameters.md` §6); reported formal-basis photon widths with external
capture do not reproduce the input there. (3) The engine's S (about 10⁻¹⁰) and S′ (about 10⁻⁸) are not smooth in their
last digits, so a one-ulp change of E_c moves Brune amplitudes by 10⁻¹¹, more where a level's Brune denominator is
near 0 or negative; the 10⁻¹² expectation for Brune/Park was not met. (4) The step-1 harness ran `X` and `X/models` in
one work tree concurrently; one run of `shift_subpoint_angle` was lost here (fixed in the step-1b harness).

**Open (proposed default in brackets).**
- Sign rule: orientation of the Brune level instead of decision 10's "largest |γ| positive". [Keep the orientation
  rule; PI to confirm. The alternative needs the sign of every final state carried to the photon channels of all
  levels that decay to it, across blocks.]
- Formal-basis external capture differs from the Brune basis (Found 2). [Decide in step 2 or 4: compute γ^ext and the
  capture normalization from the final state's Brune (observed) amplitude in every basis.]
- Formal input written by earlier builds (`param.sav`, `--no-transform`) for out-of-order groups now means another
  R-matrix; the engine notes it. [Accept; no converter. A project that needs explicit B_c values gets a keyword with
  `io/Project`.]
- Brune/Park last digits depend on E_c order and the memo at 10⁻¹¹ (10⁻⁹ in ill-conditioned levels). [Accept until
  `channel_functions` (step 3) gives smooth S, S′.]

**Next.** Step 2 (observed widths from total amplitudes), with the external-capture question above.
## 2026-11-07 (`feature/upgrade-remove-qt`) — the Qt GUI removed (decision 3); parity checklist

Branch `feature/upgrade-remove-qt` from `feature/upgrade` 4bf7aa06, not merged. Commits: the parity checklist
(alone, before the removal), the build removal, the GUI-only engine entry points, packaging/CI, the Sphinx docs, the
skills, CLAUDE.md files and README, this entry with ROADMAP.

**State.** `docs/upgrade/qt_parity.md` lists what the Qt setup utility did at 4bf7aa06 (201 items in 15 sections:
.azr blocks and columns, keyword whitelist, round-trip guarantees, each tab, the THM workspace, network fetches,
environment hooks; 32 lines, one per `tests/gui` file); the Phase 9 web UI is checked against it. Removed: `gui/`
(124 files) and `tests/gui/` (33), 47,651 lines; `BUILD_GUI`, `USE_QWT`, `BUILD_MACOS_BUNDLE`, every Qt/Qwt
find_package, `FindQwt.cmake`, `GUI_BUILD`/`start_gui`/`gui_stub.cpp`, the Qt DLL copies, Qt/Qwt in CI (apt, brew,
MSYS2, windeployqt, linuxdeploy's Qt plugin), both Dockerfiles, `build_linux.sh`, `build_windows.sh`,
`build_crc.csh`; `build_macos.sh` (it built the Qt .app/.dmg) and `examples/run_gui.sh` deleted. Engine entry points
only the GUI called, each with no caller left: `AZURECalcMCMC::RequestStop/ClearStop/SetGUIResultsCallback/
SetGUIMessageCallback` (and the pointers only they fed), `NuclearPotentialManager` getters (7 functions,
`currentType_`), `AdaptiveIntegrationGrid::GetConfig`. Docs: new `user_guide/graphical_interface.rst`; `file_menu`,
`plot`, `gui_structure` pages deleted; tab-framed user-guide pages keep their content (it documents the file) under a
note. Skills, CLAUDE.md (two .azr parsers, not three; no GUI step), tests/CLAUDE.md, pyazr/CLAUDE.md, README,
repo_consistency (no tests/gui, no SegmentKeywords.h) follow.

**Kept, and why.** `--no-gui` is accepted and ignored (scripts, `tests/upgrade` and the goldens' runners pass it).
`ThmDistortion::GlobalEnds` (unused now; the THM Diagnostics page drew it): a physics function the Phase 9 Physics/THM
tab needs through the binding. The Windows icon and the AppImage icon moved to `packaging/icons/`. The rolling
pre-release job (dev only) is unchanged but for the Windows BUILD-INFO text; macOS CI now uploads a tar.gz of the bare
executables and the job keeps the previous .dmg. `docs/upgrade` history (`cnuc_split.md`, `skills_migration.md`,
`THM.md`, `MERGE_NOTES.md`) and the "files written by the Qt GUI" notes in skills are left as written.

**Evidence.** Condor (cluster 16179), modules python + gsl only (no qt, no qwt), fresh configure of 5f09a5ce: the
`--trace-expand` log has 0 find_package/find_library calls for Qt or Qwt, CMakeCache 0 Qt/Qwt entries, `ldd AZURE2`
0 Qt libraries; the 7 CMake warnings are all in `external/` and `numcmc/`, which this branch does not touch; build: 0 compiler warnings.
`tests/upgrade/check.sh` PASS; `ctest -LE slow` 107/107 (check_skills failed in that run because the skills were not
yet edited; 3/3 of check_skills and the consistency tests pass at a90474d4); `tests/run_tests.sh` 131/132, the one
failure `light_excitation` (a), see Found; `tools/check_skills.py` 0 errors. Docs: every :doc:, toctree entry, :ref:
and image of the 37 pages resolves (script; Sphinx is not installed). The engine's sources changed only by deletions
and the `--no-gui` help text, so no output is expected to move.

**Found.** (1) `tests/light_excitation` (a) is flaky on `feature/upgrade` itself: `intEC.dat` of two runs differs at
5 × 10⁻¹¹ relative (tolerance 10⁻¹²). Same host, 4bf7aa06 build vs this branch: 8 threads 3/3 vs 3/3 failures,
2 threads 4/10 vs 4/10. The EC integrals depend on thread scheduling at that level (determinism rule). (2) pyazr's
`AzrModel` cannot create a particle pair (skill azr-levels §1); without the GUI a new compound nucleus is set up by
writing `<levels>` lines. (3) From the parity review: the Pairs table hid the light-E_x column (column 15) and Ctrl+D
help indices predate two tabs; recorded in qt_parity.md, not fixed (Qt stays on dev).

**Open (proposed default in brackets).**
- A self-contained macOS package of the CLI. [Leave the tar.gz until Phase 9's wheels; no .dmg.]
- If this merges into `dev` before the upgrade is released, the dev pre-release stops refreshing the macOS .dmg and
  the GUI downloads. [The upgrade does not merge into dev before Phase 9; the director confirms.]
- `light_excitation` flake. [A fix/ branch: fixed reduction order in the EC integrals or a tolerance with a written
  reason, after finding which loop is scheduled.]
- `ThmDistortion::GlobalEnds` in the binding. [Phase 9, with the Physics tab.]

**Next.** The director merges; Phase 9 checks the web UI against `qt_parity.md`.
## 2026-11-07 (director) — COULCC in C++, bit-identical to the Fortran (pycoulomb, local branch `coulcc-cpp`)

COULCC (Thompson and Barnett, Comput. Phys. Commun. 36 (1985) 363) translated to structured C++17 in Jakub's
pycoulomb clone, branch `coulcc-cpp` (from 9fee957; e378a33 … 37f0f9b; not pushed). One function per Fortran routine,
COMMON blocks as members of a per-call object, no mutable global state. Against gfortran (-O2/-O3, x86-64, GCC 11.5)
on 12 814 cases (real, complex ρ/η/λ, Bessel kinds, all MODE1 × KFN): 618 988 values, every IFAIL, RERR and
diagnostic counter bit-identical. Bit identity needed gfortran's complex semantics (mixed real/complex promotion,
Smith division) and REAL*16 as `__float128`; elsewhere a double-double, equal to rounding. The literal goto port it
replaces had 39 IFAIL mismatches, sign errors of J_ν, Y_ν at negative real ρ, and a data race.

**Algorithm limits found (same in the Fortran; to report to I. J. Thompson).** (1) 2λ non-integer, ρ < ½: G from
the reflection formula divides by sin χ, which cancels (−2.2 × 10⁻¹⁴ at η = 5, λ = 2.7, ρ = 0.1), G 1.3 % wrong with
IFAIL 0; the Wronskian cannot see it. (2) 2λ integer, ρ < ½: the logarithmic origin series loses accuracy to
8 × 10⁻¹¹, not reported in RERR. (3) Re λ_min ≤ −1 (outside the stated domain): Wronskian wrong by up to 54 % with
IFAIL 0. Real-argument accuracy against mpmath otherwise ≤ 7 × 10⁻¹² (CF2), Bessel ≤ 3.4 × 10⁻¹⁴.

**Next.** Merge `coulcc-cpp` with `azure2-prereqs` and add COULCC as a libcoulomb method (one adaptor; complex
arguments through `coulomb_complex()`); then it comes into `physics/coulomb/` with the rest. Not used as AZURE2's
default (`coulsm` stays). Pushing both branches is Jakub's.

## 2026-11-07 (`feature/upgrade-1-cnuc-split`) — split `CNuc`, step 1: the engine calls `level_parameters`

Branch `feature/upgrade-1-cnuc-split` from `feature/upgrade` 93f8cbc5, not merged. Commits: the plan
(`cnuc_split.md`, alone), the physics diagnostics and eigensolver option, the adapters, then one commit per moved
function (`CalcBoundaryConditions`, `TransformIn`, `TransformOut`'s formal iteration), the docs, this entry.

**State.** `cnuc_split.md` classifies all CNuc members (nucleus data, physics functions, parsing, reporting,
parameter vector), states the interface of each target part, and orders the moves (steps 1, 1b, 2–7). Step 1 is done:
B_c, observed → Brune/Park amplitudes, the external-capture subtraction, Brune → formal (GSL's eigensolver passed as
an option) and the formal → Brune iteration of `parameters.out` are `level_parameters`' code; the old bodies are
deleted (`CNuc.cpp` 3628 → 3311 lines). `CNuc` keeps the wording of its warnings, Park's J on the level and N_f^½.

**Evidence.** Engine before (93f8cbc5) vs after, all 179 `.azr` under `tests/` (T1–T3 included) × Brune, `--no-brune`,
`--use-park`: 1503 output files byte-identical; logs identical except `screening_6li` (`--no-brune`), whose warning now
gives the last energy step (118.45, 40.9187, 85.4476 MeV) instead of 0. With the own eigensolver the formal basis
would change 53 files (698 sign flips in `param.par`, observables ≤ 2.2e-11), so it stays on GSL until step 1b.
GUI off, Condor: no warnings; `check.sh` PASS; `ctest -LE slow` 106/107 (`adaptive_grid_reference`, decision 8);
`run_tests.sh` 132/132 (GUI on earlier: 131/132, `shared_point_tables` (d), see Found). T2 bench same host: fit200 79.32 vs 80.10 s (1 thread), 36.16 vs 36.28 s (8); calc 1.17
vs 1.19 s, 0.42 vs 0.41 s; no slowdown beyond one 10 ms tick.

**Found.** (1) The Coulomb memo matches energies within 1e-12 MeV, so the engine's numbers depend on the order of
its Coulomb-function calls at 1e-13; with the three orders of E_c (E − e − S, E − (S + e), E − S − e) this decided
bit-identity (`cnuc_split.md` §5). The adapter reproduces each caller's order and call sequence. (2) A group whose first
listed level is not an R-matrix level had B_c uninitialized; it now takes the first R-matrix level (none in tests/).
(3) `tests/shared_point_tables` (d) peak-RSS ratio sits at 2.9–3.0 against its limit 3.0 for both engines: it failed
once in `run_tests.sh` (GUI build) and once in 3 repeats of the new GUI-off engine (316 > 315 MB), the old engine
passed 3/3 at 305–311 MB; same RSS ranges, so a flaky margin, not this change.

**Open (proposed default in brackets).**
- The adapter options `ChannelEnergyOrder` and `ChannelQuantities` exist only to reproduce the memo history. [Keep
  until step 3; or drop them in step 1b, whose formal re-baseline could absorb a 1e-13 change of the other bases too —
  director's choice.]
- The FACTOR 3.0 of `shared_point_tables` (d) has no margin. [Measure the ratio on 10 runs and set the limit with a
  written reason, separate commit.]
- `CNuc::TransformOut`'s observed widths, `GetTransformParams`, `PrintTransformParams`, `CheckRadiativeWidths` still
  each convert γ → Γ. [Step 2: one observed-width function with γ^ext.]

**Next.** Step 1b (decision 10 on, own eigensolvers, GSL adapters deleted, formal goldens re-baselined), then step 2.
## 2026-11-07 (director) — pycoulomb prerequisites done (in Jakub's repo, local branch)

Decision 8 (a) merged (409fef62): `ctest -LE slow` 107/107 with GUI off. pycoulomb, local branch `azure2-prereqs` of
`~/AI-R-scratch/pycoulomb-src/pycoulomb` (from 9fee957; commits 7e3c834 … 7adef98; not pushed, not tagged):
build products removed; `DESIGN.md`; `libcoulomb/`, one CMake target `pycoulomb::coulomb` without Python, the
vendorable set; pure, thread-safe API `coulomb(l, η, ρ, Method)` → {F, F′, G, G′, σ_l, status}, `coulomb_range`,
`coulomb_phase`, `whittaker`, `whittaker_scaled`; methods `coul` (AZURE2 before 2025-09), `coulsm` (AZURE2 now,
default), `coulfg`, `gsl`. Parity with AZURE2 on 7150 points: 0 bit differences (`coulsm` vs `external/coul`, gsl vs
`CoulFunc`, Whittaker vs `WhitFunc`/`CoulombTable`). Threads: 0 bit mismatches at 8 threads; TSan clean.

**Accuracy against mpmath (40 digits; 1296 points, l 0–8, η 0–20, ρ 10⁻³–50).** `coul`/`coulsm`: ≤ 10⁻⁹ for η > 0
(their built-in precision 10⁻¹⁰). `coulsm` at η = 0: F_l from unstable upward recurrence is wrong for l ≥ 2 at
small ρ (G exact), but P, S, φ depend on F only through F/G, absolute error ≤ 1.4 × 10⁻¹⁶ on the grid; this is the
open neutral-F item of the T1 session. `coulfg`: JWKB inside the turning point, wrong where G > 300. `gsl`: errors to
2.5 × 10⁻² inside the turning point for η ≥ 2. All four agree to 8.3 × 10⁻¹⁰ in the oscillating region.

**For Jakub.** Add the LLNL BSD-3 COPYRIGHT text next to `coulfg2.c` (its header refers to a missing file); choose
the repository licence; push and tag. Whether to fix `coulsm` at η = 0 [keep: AZURE2 stays bit-identical; P, S, φ
unaffected]. **Next.** After the CNuc split: take `libcoulomb/`'s sources into `physics/coulomb/` (not `external/`:
 Jakub, 2026-11-07, the Coulomb functions are physics) and build
`channel_functions` on it (replacing `external/coul`, `CoulFunc`, `CoulFuncCache`, `WhitFunc`), bit-identical with
`coulsm`; COULCC added as a method when `coulcc-cpp` is done. The host must keep GSL's error handler off (`GslQuiet`).

## 2026-11-07 — `fix/adaptive-grid-criterion`: decision 8 (a) in `adaptive_grid_reference`

Branch `fix/adaptive-grid-criterion` from `feature/upgrade` e5dd8ae7, not merged. Commits: the test (47c119e4), this
entry and the notes in `ROADMAP.md` and `sanitizer_findings.md` §4.6. No engine change; the grid is thm's, unchanged.

**Criterion.** The two cases "eV/keV level crossing the window edge" check the total variation of the trapezoid
error over the scan minus |net change| (relative to the area πΓ/2) instead of its swing; the accuracy check (2 × 10⁻⁴)
and the five other cases are unchanged. A one-way drift, here the tail bias (+1.69 × 10⁻⁴) leaving the window with
the level, costs nothing; each back-and-forth excursion of size d costs 2d. Tolerance, both cases:
√3/(2π ppw²) + 10⁻⁴ = 2.10 × 10⁻⁴ at 50 points per width.
- √3/(2π ppw²) = 1.10 × 10⁻⁴: a lattice of pitch h = Γ/ppw cut by a fixed edge has the Euler–Maclaurin end term
  (h²/12) f′(edge); as the peak crosses, f′(edge) runs 0 → −max → +max → 0 (max |f′| = 3√3/(4Γ) at Γ/(2√3) from
  E_R), a variation of 4 max|f′| with no net change. Measured 1.102 × 10⁻⁴ (eV level, scan of 0.00125 Γ steps across
  the edge); in the keV case the cells from 20 Γ below E_R upward give 1.102 × 10⁻⁴. It changes the slope of the
  in-window area in E_R by at most 2/(3 ppw²) = 3 × 10⁻⁴ of itself, so it cannot reverse a gradient.
- 10⁻⁴: the smoothness tolerance of the in-window cases, here for the sum of all other reversals. keV case:
  7.3 × 10⁻⁵ from the cells below the 20-width core (graded tail meeting the fixed 5 keV coarse lattice), a ripple of
  6.8 × 10⁻⁶ peak to peak; the in-window keV case has the same ripple (swing 2.4 × 10⁻⁵).
A sampled scan shows at most the excess of the continuous curve, so the bound holds for any scan step. Values
(2000 / 200 000 steps): eV 1.96 × 10⁻⁵ / 1.14 × 10⁻⁴ (its 1 eV step, 1.25 Γ, samples the edge term coarsely; net
change −1.685 × 10⁻⁴), keV 1.79 × 10⁻⁴ / 1.83 × 10⁻⁴ (net 4 × 10⁻⁷).

**Detection.** The test built against the grid before the graded tails (dev 07581f03, through a scratch shim that
feeds it the resonance list) reproduces the 17O+α failure: in-window eV case max error 7.85 × 10⁻², swing
3.53 × 10⁻². The new crossing check fails on it: variation beyond the net change 2.87 × 10⁻¹ against 2.10 × 10⁻⁴
(net −7.49 × 10⁻², swing 7.87 × 10⁻²); 7 checks fail there, none on the current grid.

**Verified.** Condor 16150/16152, GCC 11.5, `-DBUILD_GUI=OFF`, pyazr on: no warnings in the
rebuild; `adaptive_grid_reference` all checks pass; `ctest -LE slow -j8` 107/107 pass; `tests/upgrade/check.sh` PASS.

## 2026-11-07 (`feature/upgrade-1-level-parameters`) — `Convention` record and `level_parameters` (Phase 1, code part)

Branch `feature/upgrade-1-level-parameters` from `feature/upgrade` e54e6fcc, not merged. Commits: the record
(fbda8bf5), 3j/6j and eigensolvers (becc8bfe), `level_parameters` with its physics tests (45443582), the `CNuc`
adapter and T1–T3 test (7bd0e1e9), the docs (89061935), this entry.

**State.** `physics/` has the `Convention` record (`Convention.h`; the System carries one, `basis()` reads it) and
`level_parameters` (`LevelParameters.h`, `level_parameters.md`): observed ↔ Brune/Park, B_c (first listed or lowest
level), external-capture subtraction, Brune ↔ formal (own Cholesky + Jacobi; the engine's fixed-point inverse),
Brune ↔ Park, and the convention matrix T on the γ columns (LS orders, jj by 6j, iˡ as (−1)^⌊l/2⌋, E/M signs), with
`recoupling_gaps` reporting missing channels. Channel functions are an input (`ChannelFunctions`); the adapter
`src/CNucLevelParameters.cpp` supplies them from `ChannelFunc`/`WhitFunc`. Nothing in the engine calls the new code.
Checked on Condor (GCC 11.5, GUI on, pyazr on): no warnings; `tests/upgrade/check.sh` PASS (goldens unchanged);
`ctest -LE slow -j8` 137/138, the one failure `adaptive_grid_reference` (decision 8). Against `CNuc` on T1–T3: Brune
and Park amplitudes, Park's J and B_c bit-identical; formal E ≤ 6.6e-16, formal γ ≤ 1.0e-14 (up to the level's sign,
which differs from GSL's on 3 of 19 levels); formal → Brune vs `TransformOut` ≤ 6.9e-14; six conventions give the
same formal parameters to 1.5e-15 (parameter level; the engine reads no convention). Physics tests: 10 512 recoupling
coefficients vs m-scheme to 8.9e-16; formal vs Brune level matrix identity 4.4e-15. Six mutations each fail at least
one test (`level_parameters.md` §7).

**Found.** T3's 2⁺ block lists only s = 2 of the two channel spins J = 2 allows for p + ⁷Li and n + ⁷Be at l = 1, so
T3 has no jj form for that block (reported, as intended). In the formal basis `CNuc::TransformIn` evaluates external
capture after earlier groups were made formal, so a final state in an earlier multi-level group would enter with
formal E, γ; on T1–T3 this changes nothing (T1's 14 non-zero external amplitudes go to the single-level ground-state
group). It does not explain the 54 % T1 formal-basis capture difference of `tests/upgrade/README.md`.

**Open (proposed default in brackets).**
- The two proposals (B_c from the lowest level; largest |γ| positive) are options, off. [PI decides; keep off until
  the engine switches to `level_parameters`, then turn both on together with a golden re-baseline of formal-basis
  outputs only.]
- iˡ is represented by the real signs (−1)^⌊l/2⌋ (both i^l and i^−l); no code with that convention was compared.
  [Compare with one external code's amplitudes (SAMMY or FRESCO output for one level) before a project uses it.]
- jj means j = l + j_x with x = light by default (`spin_order` picks the particle); the final coupling J = j + j_y is
  fixed. [Enough for nucleon + target; add a field only when a source needs (j_y, j).]
- External-capture order in the formal basis (above). [Decide whether to compute every γ^ext from Brune final states
  when the engine switches to `level_parameters`; not changed here.]
- E2 static moments scale the input as v × (100 (αħc)^½), the engine as (v × 100) × (αħc)^½: possibly 1 ulp apart.
  [None in T1–T3; accept.]
- No `.azr` keyword, binding or GUI for the record. [Add with `io/Project` and the binding of `level_parameters`
  once `channel_functions` exists.]

**Next.** Split `CNuc`: the engine's TransformIn/TransformOut/CalcBoundaryConditions call `level_parameters`
(bit-identical for Brune/Park shown here; formal to 1e-14, so formal-basis goldens move at that level).

## 2026-11-07 (director) — Phase 1 merges; ignored headers recovered

**Merged into `feature/upgrade`.** `feature/upgrade-1-conventions` (655c5634), `feature/upgrade-1-physics-system`
(f31d2027), `docs/upgrade-1-fix-contradictions` (7e3e5643); ARCHITECTURE §3 amended (B_c per Jπ block, 1af53f98).

**Defect found and fixed (9da89268).** The `.gitignore` pattern `azure2` also matched `physics/include/azure2/`, so
`System.h` and `JpiBlock.h` were never committed and a fresh checkout of f31d2027 did not build. The pattern is now
`/azure2`; the headers were restored unchanged from the session that wrote them. Rule for every branch: run
`git status --short --ignored` before hand-back.

**Verification.** `feature/upgrade` 7e3e5643, GCC 11.5 with GUI, Condor 16128: no compiler warnings;
`tests/upgrade/check.sh` ok; `ctest -LE slow` 133/134, the one failure `adaptive_grid_reference` (decision 8).
Sphinx not built (not installed on CRC): run `sphinx-build -W docs/source` before pushing; the cross-page `:eq:`
references from `level_matrix.rst` into `conventions.rst` are untested.

**Next (Phase 1, without pycoulomb).** `Convention` record and `level_parameters` (formal ↔ Brune/Park, B_c from
the lowest level, sign rule for formal levels), checked by the convention-invariance test on T1–T3; then the `CNuc`
split. `channel_functions` waits for the pycoulomb prerequisites (decision 6).

## 2026-11-07 (`feature/upgrade-1-physics-system`) — `physics/` target, `System` and `JpiBlock` (Phase 1, first step)

Branch `feature/upgrade-1-physics-system` from `feature/upgrade` f8b7bf22, not merged. Commits: the physics target,
the CNuc adapter with its test, the docs.

**State.** `physics/` is its own CMake target `azure2_physics` (headers `physics/include/azure2/physics/`, sources
`physics/src/`, tests `physics/tests/`); it links nothing and is added before the global GSL/OpenMP/Minuit2 include
directories. `System` (nuclei, partitions α, channels c = (α, l, s) with a_c, Jπ blocks, the level basis) and
`JpiBlock` (2J, π, nucleus, channel index list, B_c, E[λ], γ[λ n_c + k]) are immutable values with 0-based indices
and doubled spins; the choices and their reasons are in the header comments and `physics_system.md`. The adapter
`BuildSystemFromCNuc` (`include/CNucSystem.h`, `src/CNucSystem.cpp`, in AZURE2Core) builds a System from a CNuc
after `Initialize`; nothing in the engine calls it. Checked on Condor (GCC 11.5, GUI off, pyazr on):
`physics_system_from_cnuc` compares 2540 numbers of T1, T2, T3 in the Brune, Park and formal bases with the CNuc
accessors, all equal, and the structure with the `<levels>` lines (5 mutations of the adapter each fail it);
`physics_system_unit` and `physics_standalone` pass (an `#include "CNuc.h"` in physics/ fails the latter);
`repo_consistency` now also refuses non-standard includes in physics/ (a `<gsl/gsl_math.h>` fails it);
`tests/upgrade/check.sh` PASS (golden T1, T3 largest difference 0.0, T2 calc_start/calc_truth ok, recovery ok);
`ctest -LE slow -j8` 102/103, the one failure `adaptive_grid_reference`, as before.

**Open (proposed default in brackets).**
- B_c is stored per block because the engine sets it per J-group from the first level; ARCHITECTURE §3 lists it with
  the channel. [Keep it in the block; `level_parameters` may later want B_c = S_c(E_λ) per level, which is a block
  quantity too.]
- The blocks hold the parameters of the project's basis (Brune by default), not formal ones; `System::basis()` says
  which. [Keep until `level_parameters` exists; it then fills the blocks in the internal convention and the basis
  tag moves to the `Convention` record.]
- External-capture final states that are not R-matrix levels, hybrid (α, l, j) channels and the THM binding energy
  are not in the System. [Add with the items that need them: external capture with `channel_functions`, hybrid
  channels with pycoulomb, B_THM with the THM vertex.]
- The adapter reads `ALevel::GetE/GetGamma` (the loaded parameters), not the moving fit values. [Right for setup
  state; per-evaluation E, γ come from `level_parameters`.]

**Next.** `channel_functions(α, l, E)` on the System's partitions; then split `CNuc` (ROADMAP Phase 1).
## 2026-11-07 (`feature/upgrade-1-conventions`) — conventions page, documentation part (ROADMAP Phase 1)

Branch `feature/upgrade-1-conventions` from `feature/upgrade` f8b7bf22, merged 655c5634. Commits: 754ac7d7 (test
`tests/coulomb_nuclear_sign`), 8faa430d (`docs/source/theory/conventions.rst` + Theory toctree), and this
HANDOFF/ROADMAP commit. No C++ changed.

**State.** The page states the present engine's conventions, each with file:line: s = j1 + j2 light first (LT
III.2a), J = s + l in the spin-amplitude route, Condon–Shortley CG/Racah; no i^l on particle channels (only the
E external-capture amplitude carries −i i^(l+L−lf), M1 +i); φ = arctan(F/G), ω = σ_l − σ_0, U = ΩΩ[1 + 2iP^½γAγP^½],
T = e^{2iω} − U; B_c = S_c(E of the first listed level) in the formal basis; sign of γ = sign of the input Γ/ANC/
amplitude (Brune, Park; photon total width with EC; static-moment sign flip for half-integer J); photon widths
never in the level matrix except --use-rmc; Blatt–Biedenharn z = (−1)^(s'−s) Z Z / 4 (photon: (−1)^(1+s−jf));
elastic Coulomb term with C_α = +(4π)^−½ η csc² e^{−2iη ln sin} and Re[i C T*]; identical-particle factors; the
planned Convention record (marked planned). `tests/coulomb_nuclear_sign/check.sh`: p+12C hard sphere vs an
mpmath partial-wave sum, 30 points, A- and R-matrix, Brune and --no-brune, max rel diff 8.4e-11 (Condor 16122);
the interference reversed moves the truth by 7–74 %. Sphinx not built (not installed on CRC; labels and :doc:
targets checked by script).

**Existing docs that contradict the code (not changed here).**
- `docs/source/developer/level_matrix.rst`: "L_c = L_0(c) is the shift function" (code: L⁰ = S − B + iP); U
  written without Ω_cΩ_c'; T_cc = e^{2iω} − e^{iω}e^{iφ}(1+U) and T_cc' = −e^{iω_c}e^{iφ_c'}U (code: Ω_cΩ_c' with
  e^{−iφ} on both channels); "ω_c = σ_l" (code: σ_l − σ_0); [1−RL] = δ − iRP (code: 1 − RL⁰ in the standard
  formalism); interference "(i/√π) g f_C* T P_l + c.c." (code: Re[(i/√π)(π/k²) g C_α T* P_l], opposite sign and
  half the size if f_C means C_α).
- `docs/source/reference/data_formats.rst` (phase-shift data): "−90° < θ_lab < 90°" names the wrong quantity; the
  phase shift is in (−90°, 90°], and in [0°, 180°) for identical pairs (`GenMatrixFunc.cpp:555-564`).

**Open (proposed default in brackets).**
- Sign of internal relative to external capture: fixed by the code's formulas, checked by no independent
  calculation. [An m-scheme or potential-model truth for one resonance + EC, as `tests/ec_e1m1_phase`.]
- The orbital–spin order is pinned only relative to the Blatt–Biedenharn route (polarization Gate 1). [Accept;
  the invariance test of the record will cover it.]
- LT and DB were not re-read (no PDF reader on CRC); the page says which comparisons rest on the code's citations.
  [Re-read LT VII–VIII and DB §2–3 on a machine with a PDF reader and add equation numbers.]
- Formal basis: B_c = S_c(first listed level) makes the formal parameters depend on level order, and GSL's
  eigenvector sign makes a formal level's overall sign arbitrary. Observables are unaffected. [Fix both in
  `level_parameters`: B_c from the lowest level, sign so that the largest |γ| of a level is positive.]
- Fix the two contradicting pages. [Separate `docs/` commit, by the director's choice.]
  Done by 57bcd03d on `docs/upgrade-1-fix-contradictions` (not merged): all listed contradictions confirmed
  and fixed; also the Reich–Moore sign on the A⁻¹ diagonal (`+= iγ²`, code `−= iγ²`).

**Next.** Convention record and `level_parameters` (code part of the ROADMAP item).

## 2026-11-07 (`fix/edep-ownership-tsan`) — copy-on-write ownership TSan can see (sanitizer_findings.md §4.2)

Branch `fix/edep-ownership-tsan` from `feature/upgrade` 5b021003, not merged. Commits: c556f052 (the fix), and the
findings/HANDOFF commit.

**State.** `EPoint::edep_` and `legendreP_` are `SharedTable<T>` (`include/SharedTable.h`): a shared_ptr plus an
atomic holder count (relaxed increment, release decrement, acquire load in `SoleOwner`); the use_count() + fence test
is gone. Checked on Condor: GCC 11.5 build with GUI and pyazr; `tests/upgrade/check.sh` PASS at 8 threads (golden
outputs unchanged); `tests/run_tests.sh` 131/131; `ctest -j8` 135/140, the five failures inherited from origin/thm.
clang TSan on a temporary branch with this fix and `fix/fit-iteration-counter` (`feature/upgrade` + fecf4924 +
c556f052): 54 steps, no report file (before: 25 files, 128 reports, 8 AZURE2 sites). Control: the acquire load made
relaxed brings back 75 reports at `EDependentTables::Clear` (checks section).

**Open (proposed default in brackets).**
- Under TSan, `shared_point_tables` (d) now fails without a report: peak RSS 1162/1164 MB against the limit 1149 MB
  (1133 MB before), from TSan's per-atomic synchronization records; plain builds are unaffected and the ratio is
  near the limit in plain builds too (one failure in six runs, on the branch that does not touch the tables).
  [Run that check under TSan with `FACTOR=3.5` from `sanitize_upgrade.sh`, with the reason in the script; the
  plain-build limit stays 3.0. Not done here: it changes a test limit, which is the director's call.]
- `pyazr_upgrade_thread_independence` failed once in a full `ctest -j8` on this branch (T3 fit `param.fit`, 2.3e-5
  relative, 1 vs 4 threads; parameters equal); four isolated repeats on this branch, on `fix/fit-iteration-counter`
  and on `feature/upgrade` all pass. [Cause not found; not reproduced. Rerun the test
  a few times under load before the next merge.]

**Next.** Merge with `fix/fit-iteration-counter` (both add a HANDOFF entry here: keep both).
## 2026-11-07 (`fix/fit-iteration-counter`) — fit iteration counter atomic (sanitizer_findings.md §4.1)

Branch `fix/fit-iteration-counter` from `feature/upgrade` 5b021003, not merged. Commits: fecf4924 (the fix), and the
findings/HANDOFF commit.

**State.** `EData::iterations_` is `std::atomic<int>`; `EData::Iterate()` returns the number of this evaluation
(`fetch_add`), `AZURECalc::operator()` uses it, and the progress line is written under a mutex. Checked on Condor:
GCC 11.5 build with GUI and pyazr; `tests/upgrade/check.sh` PASS at 8 threads (golden outputs unchanged);
`tests/run_tests.sh` 131/131; `ctest -j8` 135/140, the five failures inherited from origin/thm (fourth session,
below). A T1 fit with Minuit2's numerical gradient prints the same 73 progress lines as `feature/upgrade` at one
thread, and the numbers 10, 20, …, 730 once each at eight. clang TSan on a temporary branch with this fix and
`fix/edep-ownership-tsan`: 54 steps, no report (before: 25 report files; 48 reports at this finding's three sites).

**Next.** Merge with `fix/edep-ownership-tsan` (both add a HANDOFF entry here: keep both).
## 2026-11-07 — `fix/adaptive-grid-edge-swing`: adaptive_grid_reference, analysis only (decision needed)

**What fails.** "eV level crossing the window edge" (Γ = 0.8 eV, 10 keV window, 50 points, level moved from
0.9 to 1.1 of the window in 2000 steps): swing of the error 1.685 × 10⁻⁴ of the area πΓ/2, tolerance 10⁻⁴. The
check is from dev (d216517); the grid is thm's (02ae97b), kept at the merge of dev 2d97b1f0, where MERGE_NOTES.md
already lists this failure as a decision.

**Cause.** thm's grid has geometric tails whose step is r times the distance from the level,
r = TailFactor() = 2.5/50 × √(20/2) = 0.158 for a 20-width core. A chord over a 1/x² wing over-counts by about
r²/2, so the area comes out high by r²Γ/(2π · core) ≈ 2 × 10⁻⁴; measured +1.69 × 10⁻⁴ inside the window, constant
to 2.2 × 10⁻⁷ as the level moves (the "eV level, 10 keV window" case). In the crossing scan the error starts at
that bias and ends at ~0 with the level outside, so the swing equals the net change: swing 1.685 × 10⁻⁴, net
change −1.685 × 10⁻⁴, total variation beyond |net| 2.0 × 10⁻⁵, largest jump between neighbouring positions
(1 eV) 9.4 × 10⁻⁶. For any grid whose bias b passes the case's accuracy check (|b| < 2 × 10⁻⁴), this swing is at
least |b|: here the swing check is a 10⁻⁴ accuracy bound, not the smoothness check it is described as. It is not
an oscillation (the 17O+a failure it was written for had 3.5 × 10⁻² with a 10 eV period). dev's grid passes
because its bias is 7.2 × 10⁻⁵.

**Change.** a test commit that only prints the net change and the variation beyond it for every case; no check,
tolerance or case changed. The test still fails.

**Options.** (a) Keep thm's grid; for the crossing cases check the variation beyond the net change (2.0 × 10⁻⁵
eV case, 1.8 × 10⁻⁴ keV case at present, so its tolerance needs its own argument) and leave the bias to the
accuracy check. (b) Lower tailRatio 2.5 → ~1.9 (b ∝ r², ~30 % more tail nodes): every fold moves by ~10⁻⁴, so the
regression pins and the upgrade goldens (1e-6) move; not done. (c) dev's grid: moves THM outputs. [Proposed
default: (a), decided by the PI, since it changes a reference test's criterion.]
## 2026-11-07 — `fix/c12c12-kernel-pin`: pair_identity and the THM examples' chi-squared

**Cause.** `pair_identity` pinned the engine total of `examples/c12c12_tumino2018` at 112.0376 (data 111.613 +
priors 0.4248), the value of the merge of dev 64ed575. The merge of dev 2d97b1f0 brought dev 7d24c88: Gaussian
resolution kernels to 5σ instead of 3σ. MERGE_NOTES.md measured the move (data 111.613 → 111.661, direct segments
5-8 only; convolutionRange = 3 restores thm's files byte for byte) but the GUI test was not built and the
documented totals were not updated.

**Why the new value is right.** A Gaussian cut at ±kσ and renormalised has variance
σ²[1 − 2kφ(k)/(2Φ(k) − 1)]: 0.9733 σ² at k = 3 (the fold used a resolution 1.3 % narrower than stated),
1 − 1.5 × 10⁻⁵ at k = 5.

**Change.** f0981b53: the pin is 112.0856 (tolerance 0.01 unchanged). bb1bd4f0: examples.rst,
thm_implementation.rst and the azure2-thm calibration table carry the measured totals (build-upgrade, 4 threads):
c12c12 112.086 (111.661; THM segments 1-4 61.02, unchanged), li7_tumino2006 448.620 (439.149), li6_pizzone2011
158.967 (157.041). The other five THM examples give their documented values.
## 2026-11-07 — `fix/mcmc-shift-prior-centre`: mcmc_priors and the MCMC tab's shift-prior centre

**Cause.** `mcmc_priors` (thm, 2026-10-09) expected the energy-shift prior of segment 1 at its start value
(0.004 ± 0.002 MeV) and the sqrt(E) coefficient's at 0.003 ± 0.001. dev 0dd28e5e (2026-10-25), merged later,
centres both on no shift: the fit's penalty (`ESegment::AddShiftPenalty`) and the sampler's automatic prior
(`AZURECalcMCMC::BuildAutoPriors`). The tab follows it; the test was out of date. Reading the code showed a second
difference: the engine centres the shift prior on `GetNominalEnergyShift()`, which a
`segment_N_energy_shift prior_centre c` row of `<parameterSettings>` moves; the tab used 0 regardless.

**Change.** dc4cbaa6: the test expects 0 ± 0.002 MeV and 0 ± 0.001. 2a137b4c: the tab takes
`GetNominalEnergyShift()` (as it takes `GetNominalNorm()`); test section 2b, with the row 0.004, expects 0.004 ±
0.002 MeV for the shift and 0 ± 0.001 for the sqrt(E) term. The old line gives 0 there, so 2b fails on it.

**Overlap.** None with the `EData`/`AZURECalc` counter or `EPoint` ownership work: only `gui/src/MCMCTab.cpp`.

**Open.** `ESegment.h` describes `GetNominalEnergyShiftSqrt()` as "the penalty's reference", but the penalty uses
0. [Correct the comment on the next branch that touches `ESegment`.]
## 2026-11-07 — `fix/gui-pair-columns`: thm_opt_in and light_excitation_gui after the merge of dev 3264896d

**Cause.** The merge of dev 3264896d into thm (made without Qt; MERGE_NOTES.md) added dev's light-particle
excitation e1 to `PairsData` and kept thm's THM columns: pairs model column 15 is e1, 16 the THM binding energy
B (`PairsModel::BindingEnergyColumn`, 15 before), 17 the compound nucleus; the GUI always writes B and the width
flag (columns 32, 33) and e1 as column 34 when it is nonzero; NucLine tells 1, 2 and 3 trailing fields apart. Two
GUI tests were never built against that layout. `thm_opt_in` wrote B into column 15, so it set e1 = 2.2246 MeV and
the project rightly opened with THM off. `light_excitation_gui` (from dev) expected 16 fields and 31/32 columns.
The code is right in both cases; the tests were out of date.

**Change.** `thm_opt_in` uses `BindingEnergyColumn`. `light_excitation_gui` checks SIZE 17, the three named
columns, 33 columns written for a pair without e1 and 34 with it, a 31-column dev file written with 33, and (new)
the 33/34-column output read back with e1 on pair 2 only. `reference/data_formats.rst` now states the
trailing-field rule (it described dev's 32nd column only).

**Verified.** The four `fix/` branches of this session merged on 5b021003 (scratch, detached), GUI on, GCC 11.5,
Condor 16112/16113: `ctest -LE slow` 130 of 131 pass, only `adaptive_grid_reference` fails; `tests/upgrade/check.sh`
PASS (no engine code changed). The failing tests of this branch failed on origin/thm c6ddb741 (plain GUI build).

**Open.** The GUI turns a 31-column classic file into 33 columns on save; the engine and pyazr read both. Whether
the THM columns move behind e1 when thm goes into dev is the decision MERGE_NOTES.md already lists. [Keep the
present layout until that merge.]

- With `build.sh` the pyazr module lands in `BUILD/pyazr-out/` (`PYAZR_OUTPUT_DIR`), but ctest imports
  `pyazr` from the source tree: 31 pyazr tests print `skip` and `pyazr_angular_coefficients` fails at import (it
  has no skip guard). The counts above are with the module copied into the scratch tree's `pyazr/`. [Have ctest
  put `PYAZR_OUTPUT_DIR` on `PYTHONPATH`, and give `angular_coefficients_test.py` the skip guard; own `fix/`
  branch.]

## 2026-11-07 (fourth session) — Phase 0 merged; failures inherited from origin/thm

**State.** `feature/upgrade` carries all Phase 0 branches: rules and skills (09379840), T1/T3 (a80352e6), T2
(356c6770), sanitizers (73290ee2). `tests/upgrade/check.sh` runs T1, T3 and the two T2 calculations against the
golden outputs (e4cadf51); Condor 16086 at 8 threads: PASS. Both race fixes are on `feature/upgrade` (ROADMAP item
ticked).

**Five ctest failures come from origin/thm, not from the upgrade.** A plain build of origin/thm c6ddb741 (GCC 11.5,
Qt GUI on) fails the same five tests as `thm-review` 1174fdb2:

| test | failing check |
|---|---|
| `thm_opt_in` | "binding energy: THM on" |
| `pair_identity` | 112.0856 against the expected 112.0376 |
| `mcmc_priors` | energy shift 0.004 ± 0.002 MeV and the sqrt(E) coefficient prior 0.003 ± 0.001 |
| `light_excitation_gui` | `PairsData::SIZE` 16 and the 31/32-column pair lines |
| `adaptive_grid_reference` | "eV level crossing the window edge": swing 1.685 × 10⁻⁴ of the area, tolerance 10⁻⁴ |

The log is in the scratch build directory of the analysis machine; the commands are `ctest -j8` on a GUI build.
None of the five touches the THM form-factor GUI added on `thm-review` (7c543690). Each needs a `fix/` branch.

**Next.** `fix/fit-iteration-counter` and `fix/edep-ownership-tsan` (TSan green); `fix/` branches for the five
inherited failures; then Phase 1.

**Later the same session.** Four of the five inherited failures fixed and merged (05e08e06, 1d35e1d4, 22214e14;
analysis of `adaptive_grid_reference` in c51f973e, decision 8). `fix/pyazr-tests-engine-dir`: a build with
`PYAZR_OUTPUT_DIR` outside `pyazr/` (as made by the scratch build helper and the sanitizer builds) skipped 31 pyazr
tests; ctest now sets `PYAZR_ENGINE_DIR`, which the package adds to its search path. `ctest -R pyazr -LE slow`:
74/74 passed, none skipped (Condor 16114). The four tests outside `tests/pyazr` that import pyazr
(`beta_delayed_alpha`, `consistency`, `electron_screening`, `hybrid_potential`) do not get the variable yet.

**TSan fixes merged.** `fix/fit-iteration-counter` (a2392e40) and `fix/edep-ownership-tsan` (f8b7bf22): clang TSan on
both, 54 steps at 8 threads, 0 report files (before: 25 files, 128 reports, 8 sites). Merged `feature/upgrade`
f8b7bf22, GCC 11.5 with GUI, Condor 16120: `tests/upgrade/check.sh` PASS; `ctest -LE slow` 130/131, the one failure
`adaptive_grid_reference` (decision 8).

**Open (proposed default in brackets).**
- `shared_point_tables` check (d), a 4-thread fit's peak memory below 3.0 × the calculation's, fails under TSan
  (1163 MB against 1149 MB) with no report, and failed once in six plain runs. [Factor 3.5 under
  `sanitize_upgrade.sh` only, with the reason in the script; PI to confirm, since it is a test limit.]
- `pyazr_upgrade_thread_independence` failed once: the T3 fit's χ² differed by 2.3 × 10⁻⁵ relative between 1 and 4
  threads, fitted parameters equal; twelve reruns on three branches passed. Cause not found. [Watch; investigate if
  it recurs.]

## 2026-11-06 (sanitizer session) — Phase 0: TSan and ASan at 8 threads, script, CI, findings

Branch `feature/upgrade-0-sanitizers` (from `feature/upgrade` 356c6770), local, not pushed. Commits: c692bb8b
(cherry-pick of 189f3023 from `ai-r/dev`: `PYAZR_OUTPUT_DIR`), 557bb0ba (`AZURE2_SANITIZE`), cfd281ab
(`scripts/sanitize_upgrade.sh`, `scripts/sanitizers/`), 16309cf8 (TSan ignorelist), 8572b21c (workflow), 635b0a1a
(`sanitizer_findings.md`, installation option row), and the roadmap/HANDOFF commit.

**State.**
- `-DAZURE2_SANITIZE=thread|address|address,undefined` instruments every target, Minuit2 and coul included.
- `scripts/sanitize_upgrade.sh thread|address`: 54 steps at 8 threads (T1/T3 calc, extrap, rate, gradient fit; T2
  calc and a 30-call fit; 26 mapping/shift/target `check.sh`; 3 regression projects; 4 integrator reference tests;
  12 pyazr tests in-process). Exit 1 on a report, 2 on a failure without one, 3 on setup.
- `.github/workflows/upgrade-sanitizers.yml`: clang 18 + libomp/libarcher on Ubuntu 24.04, nightly, dispatch, and on
  a push touching the workflow/script/suppressions.
- Runs on Condor (`sanitizer_findings.md` §3): ASan clean with clang 23.1.4 and GCC 15.2. TSan, after the two
  third-party entries (libomp suspend/resume, suppressed; Minuit2 call counter, ignorelisted): 25 report files,
  all at two AZURE2 sites.

**Findings for `fix/` branches (not fixed here; details, stacks and proposed fixes in `sanitizer_findings.md` §4).**
1. §4.1 `fix/fit-iteration-counter` — real race: `EData::iterations_` read and incremented, and the progress line
   written to `std::cout`, from every thread of Minuit2's numerical gradient (`AZURECalc.cpp:565-713`). No fit value
   depends on it. Fix: atomic counter with `fetch_add` returning the iteration, progress line under a mutex.
2. §4.2 `fix/edep-ownership-tsan` — the copy-on-write ownership test `SoleOwner` (relaxed `use_count()` + acquire
   fence, 3fba1192) is ordered under the C++ memory model; TSan does not model standalone fences and reports
   `EDependentTables::Clear` against `ClearEDependentRows` (`shared_tables_race`, `shift_fn*`, pyazr `shift_fn`).
   Fix: an ownership count TSan can see (atomic holders count, acquire load), or wait for Phase 6, which removes the
   tables. The TSan job stays red until one of the two lands.
3. Third party, recorded only: Minuit2's `fNumCall++` from the parallel gradient is a real race (call count and
   `--max-calls` approximate under threads).
4. Not a sanitizer finding: `adaptive_grid_reference` fails on every build of 356c6770, plain GCC 11.5 included
   ("eV level crossing the window edge", swing 1.685e-4 > 1e-4). Needs its own branch.

**Not done / not verified.** The workflow has never run on GitHub (no push from here): the apt layout (runtime names,
`libarcher.so` and `omp.h` locations) is handled by lookups, not tested. The cluster's clang module carries no
compiler-rt and no libomp; the runs used both built from the LLVM 23.1.4 sources in scratch (outside the repository).
GCC TSan was not run: with libgomp it cannot be correct (§2 of the findings).

**Open (proposed default in brackets).**
- GitHub fires the schedule only from the default branch (`qt5`, human-only). [Leave the file on `feature/upgrade`;
  rely on the push trigger and manual dispatch until the upgrade merges, or the PI adds the file to the default
  branch.]
- Should the TSan job be green before §4.2 is fixed? [No suppression of AZURE2 code: take `fix/edep-ownership-tsan`
  next, a small change in `EPoint`; it also removes the last fence-based argument from the copy-on-write code.]
- The ASan job warns, not fails, on `adaptive_grid_reference` (exit 2 is a warning in the workflow). [Keep: build.yml
  gates test failures.]

## 2026-11-07 (T2 session) — test project T2, its golden outputs and the benchmark

Branch `feature/upgrade-0-testprojects-t2` (from `feature/upgrade` ff75fc28), commits 3184b1ee … 8870504f, not
merged.

**State.**
- `tests/upgrade/T2/`: p + ¹⁵N → ¹⁶O through an apparatus, 13 segments, 330 points, 19 free parameters (6 level,
  13 nuisance): thin, thick, implanted (fitted density) and gas (Doppler, E-dependent resolution) targets with stopping
  and straggling; beam profile; Q attenuation; constant, √E and f(E) shifts with groups; a norm group; norm,
  shift, density, abundance and covariance priors; a data covariance; SUM and RATIO through α₁; a natural-target
  sample; THM HOES (two runs, one profiled scale, 20 keV c.m. fold) and OES points. `make_t2.py` holds the truth and
  makes the data with the engine (seed 20261106); README states coverage, truth and TODOs.
- Recovery (`check_recovery.py`, HTCondor 8 threads): noisy χ² 294.17/330, max |pull| 2.15; Asimov χ² 0.17,
  max |pull| 0.38. Golden outputs in `tests/upgrade/golden/T2` (1e-5 calculations, 1e-3 fit;
  `regen_compare.py`). Benchmark `scripts/bench_upgrade.sh`, first row in `tests/upgrade/bench/T2.tsv`:
  calc 2.67 s / 1.77 s, 200-call fit 170 s / 96.7 s, peak RSS 238 / 475 MB at 1 / 8 threads.

**Findings for later phases.**
- With THM segments the CLI forces the analytic-gradient MIGRAD; on T2 it spends > 3 min per gradient call and
  did not converge in 10 min. THM has no analytic residual Jacobian, so `--use-lm` falls back to numerical MIGRAD
  (1082 calls, 585 s). The T2 scripts use `--use-lm`. Phase 2/5: a THM Jacobian.
- Exact Coulomb functions at re-placed sub-points are about 2/3 of a fit evaluation; a Coulomb table makes the T2
  fit 8× faster (72 s) with the same parameters to 8 digits.
- Fit speed-up from 1 to 8 threads is 1.76 and peak memory doubles (the 200-call benchmark).
- In the Brune basis the fit fixes γ, but the `.azr` enters Γ at the level energy: moving a level's energy or a
  free width in the file changes the γ of its fixed channels. A synthetic project must write the fixed Γ that give
  the true γ at the start values (`make_t2.py` `match_fixed`); T1 and T3 need the same care.

**Not done.** S(E) fold of the OES points (no engine support; README TODO); a fitted Doppler temperature (¹⁶O has
no resonance narrow enough); a two-nucleus natural target (¹⁴N(p,α) closed); TSan/ASan runs; T2 is not wired into
CTest or `run_tests.sh` (its fit takes 10 min).

**Next.** Merge the branch. Open: T2 has its own scripts (`T2/check_recovery.py`, `golden/T2/regen_compare.py`),
written in parallel with T1/T3's `recover.py` and `golden/golden.py`; folding T2 into their `PROJECTS` lists is a
small follow-up (proposed default: keep the T2 scripts, and let `tests/upgrade/check.sh` call
`golden/T2/regen_compare.py --only calc_start --only calc_truth`, 3 s, but not the 10-min fit). Decide whether the
T2 recovery runs in nightly CI (10 min at 8 threads, or 72 s with a Coulomb table built in the job).

## 2026-11-07 (third session) — Phase 0: rules into CLAUDE.md; ten-skill set; skill checker

**Branch** `feature/upgrade-0-rules-skills` (from `feature/upgrade` ff75fc28), local, not pushed.

**Done.**
- `CLAUDE.md` "Skills" section replaced by `RULES.md` §3 (use / behaviour / change; sizes; no history except
  expiring version notes; one owner; descriptions ≤ 600; no personal paths; the checker). The "date an entry" rule
  is gone; the frames-and-units pointer now names `azure2-check`.
- `.claude/skills/` is the ten-skill set plus `rmfit-campaign` (2,280 → 358 lines + three references). The ten old
  skills are deleted. `docs/upgrade/skills_migration.md` maps every section of every old skill (by line range in the
  ff75fc28 files) to its new location or a reason for dropping it (498 rows; coverage checked by script).
- `tools/check_skills.py` (R4, R5, R6, R7, R8, docs/reference/script links), its unit test
  `tests/consistency/check_skills_test.py` (12 cases), ctest `check_skills_unit` and `check_skills`, a local
  `.pre-commit-config.yaml` hook and `.github/workflows/skills.yml` (runs on push/PR touching skills or docs). The
  checker exits 0 on the branch, also with `--strict`.
- Developer content from the skills moved to `tests/CLAUDE.md` (CI portability, regression vs reference, m-scheme
  traps) and `gui/CLAUDE.md` (headless dialog checks). `reference/command_line.rst` gained the gradient-diagnostic
  and `AZURE_DUMP_NONFINITE` variables.

**Not checked.** `sphinx-build -W` (Sphinx is not installed on the analysis machine; the two RST edits were read by
eye). The CI workflow has not run (nothing is pushed).

**Open.**
- Version notes in the skills expire on 2027-02-01; the checker then fails until each is removed or renewed.
- The user-level `~/.claude/skills/azure2` must narrow its description so it does not load next to the repo skills
  (`RULES.md` §3, outside the repository).
- `pyazr/azrfile.py` (a docstring near line 3255) still says "see the azure2-eval skill"; left alone because this
  branch changes no engine or pyazr code. The owner is now `azure2-calculate`.
- Corrections made while moving content (reviewer please confirm): masses in the old project builder were called
  nuclear (they are atomic); `nuclear_mass()` and `zero_all_widths()` cited by old skills do not exist in pyazr; the
  CLI band prompt is asked in modes 2, 3 and 4; the level-merge rule follows `levels_channels.rst`; in
  `azure2-resonance-search` the gain is "control minus add" (positive = improvement).
## 2026-11-06 (evening, EST) — Phase 0: T1, T3, golden outputs, identity checks (branch `feature/upgrade-0-testprojects-t1t3`)

**State.** Branch `feature/upgrade-0-testprojects-t1t3`, cut from `feature/upgrade` c1cc92c7 (feature/upgrade ff75fc28
merged in for STYLE.md). Ready for review and merge; nothing pushed.
- 63c08471 `tests/upgrade/T1` (p + ¹²C → ¹³N, 442 points, 13 free parameters) and `tests/upgrade/T3` (⁸Be, 494
  points, 8 free parameters): generators with the true parameters, synthetic data (seeds 13 and 8), `.azr`, start
  models, READMEs with truth, data tables and coverage; `tests/upgrade/lib/upgrade_project.py`; `tests/upgrade/recover.py`.
- 3b48a101 golden outputs of engine c1cc92c7 (`tests/upgrade/golden/`, calc/extrap/rate, `--no-coulomb-table`,
  2 threads) with `golden.py` (compare at rtol 1e-6, `--write` to regenerate); `tests/upgrade/check.sh` (golden +
  recovery, run by `tests/run_tests.sh`, ~50 s).
- fbad9b6a, bf5f1713 seven identity tests `tests/pyazr/upgrade_*_test.py` (ctest `pyazr_upgrade_*`; the thread test is
  `LABELS slow`); `tests/upgrade/README.md` lists every §10 identity with result, and the pending ones.

**Verified (feature/upgrade build c1cc92c7).** Recovery: T1 χ² 462.3 at the fit vs 482.3 at the truth (N = 442), largest
pull 3.64 (5/2⁺ energy), noise-free fit ≤ 0.007 σ; T3 χ² 498.7 vs 509.7 (N = 494), largest pull 1.82, noise-free
≤ 0.036 σ; eight extra noise seeds for T1: pull means within ±0.7, standard deviations 0.7–1.5. Golden: two runs agree
bit for bit; an edited golden value fails. Identities: Coulomb F, G vs mpmath 4.5e-10, Wronskian 4e-16; unitarity to
print precision; optical theorem 2.6e-11; Breit–Wigner area 1.1e-7; Park = Brune exactly, formal = Brune 2.2e-11
(see open 1); gradient vs Richardson differences 1.2e-8 / 2.7e-6 of max|∇χ²|; threads 1 vs 4 equal to 3.6e-11 (open 2).
A build of this branch (CLI + pyazr, no GUI) ran `ctest -R "pyazr_upgrade|consistency"`: 8/8 passed (none
skipped), and `tests/upgrade/check.sh` passed. The full `ctest -LE slow` and `tests/run_tests.sh` were not run: the
branch adds tests and changes no engine, pyazr or GUI code.

**Open (proposed default in brackets).**
1. T1 capture **with external capture** in the formal basis (`--no-brune`) differs from Brune/Park by up to 54 %
   (E_c.m. = 1.568 MeV); without external capture it agrees exactly. The EC amplitude depends on B_c and the transformation
   of observed photon widths to the formal basis does not carry it. [Treat Brune/Park as the reference; resolve in Phase 1
   with `level_parameters`; meanwhile a warning for `--no-brune` with an EC mask would protect users — a `fix/` branch
   on dev, after James's group agrees.]
2. Thread dependence at 3.6e-11 relative in T1 capture to the excited (3/2⁻) final state, in the calculation and after a
   fit, also between two 4-thread runs; fitted parameters identical. [Accept until Phase 6; its gate tightens
   `upgrade_thread_independence` to bit for bit. Worth a TSan run on T1 in the Phase 0 CI box.]
3. Neutral channel F_l at small ρ and l ≥ 3 is computed with cancellation (7.9e-6 relative at l = 4, ρ = 0.18); P, S,
   δ_hs unaffected. [Fixed by the pycoulomb functions of Phase 1; the test then tightens to 1e-8.]
4. Photon widths are not in the level denominator (Γ = Γ_particle); `upgrade_breit_wigner` pins it. [Keep; state it in
   `docs/source/theory/conventions.rst` in Phase 1.]
5. `tests/beta_delayed_alpha` labels its energy column "E_cm(12C+a)", but the engine reads a β-delayed data energy as the
   lab energy of the light particle from the decay at rest (E_c.m. = E·16/12 there). The data are consistent with the
   engine, only the label is wrong. [Fix the label on dev.]
6. Golden files are generated outputs outside `tests/*/expected/`, which the root CLAUDE.md "Commits" rule allows only
   there. [Name `tests/upgrade/golden/` as the second allowed place in the Upgrade section.]
7. T1's convention round trip is a TODO for Phase 1 (no `Convention` record yet).

**Next.** Merge this branch; add T2's golden files and recovery to `golden.py` / `recover.py` (`PROJECTS`, the project
list) when T2 lands; the Phase 0 CI box should run `tests/upgrade/check.sh` and the `pyazr_upgrade_*` tests.

## 2026-11-07 (second session) — branches imported; THM merged in; Phase 0 started

**How to continue on any machine.** Pull `feature/upgrade` and any `feature/upgrade-<phase>-<topic>` branch listed
below, read this file, then `ROADMAP.md`. Every phase branch is cut from `feature/upgrade` and merged back with
`--no-ff` when its gate passes. Nothing in the repository depends on the machine it was written on: scripts take the
engine binary and build directories from arguments or environment variables (`AZURE2_BINARY`, `PYTHONPATH` for the
pyazr module) and document their defaults. Commits are local until Jakub pushes them; an agent never pushes.

**State.**
- `feature/upgrade` now contains: the plan (f99661cf); `origin/thm` c6ddb741 merged in (d20718f8), so the rebuild
  starts from the THM line (which already carries `dev` up to 69826014; the last five `dev` commits, aab1de0d and
  its four parents, come through the plan commit's base); `fix/edep-cow-acquire` 3fba1192 merged (37978592); an
  introduction in `README.md` and the `CLAUDE.md` Upgrade section (c1cc92c7).
- `thm-review` (4882ad04): the uncompiled `vertexModel=ff` GUI code from the analysis machine is applied in the
  working tree and being built and tested before it is committed; then `origin/thm` is merged into it.
- Phase 0 work in progress on `feature/upgrade-0-testprojects-t1t3` (T1, T3, identity checks) and
  `feature/upgrade-0-testprojects-t2` (T2, `scripts/bench_upgrade.sh`). Next: `feature/upgrade-0-rules-skills`.

**Open.** The decisions listed in the entry below are unchanged. The race-fix status of `thm` (a6bc6b9) is resolved on
this branch by the merge; `origin/thm` itself still lacks it.

## 2026-11-07 — plan written; branch created; nothing implemented yet

**State.**
- Branch `feature/upgrade`, cut from `origin/dev` aab1de0d. Contains only this plan (`docs/upgrade/`) and the "Upgrade"
  section of `CLAUDE.md`. No code changed.
- The plan: `ARCHITECTURE.md` (target design: three lines, modules, one-way rule, stopping, Coulomb via pycoulomb,
  conventions, physics explorer, tests, front ends), `THM.md` (simple valid THM with DWBA vertices), `ROADMAP.md`
  (phases 0–9 with checklists and gates), `RULES.md` (rules and the ten-skill set), `figures/structure.html` (open in
  a browser; generated by `figures/make_figures.py`).
- The plan came from six read-only surveys of `dev` (data/point model, experimental effects, fit loop and threading,
  fix history since 2026-08-01, classes vs physics with pyazr, skills), done on 2026-11-06.

**Work that exists only on the analysis machine (almalinux), shipped as patches in the handoff package.**
- `thm-review` (base `origin/thm` 30d3bd9, local, never pushed): 073fe48 `backgroundPositive=on` (engine, pyazr, docs,
  tests), 08ccea3 its GUI box, 4882ad0 `vertexModel=ff` full form factor (engine, pyazr, docs, tests; validation:
  Green identity 2.6e-10, R values vs formula sheet ≤ 1e-6, `ffVolume=off` byte-identical to `perlevel`, light suite
  69/69, THM ctests 30/30). Plus uncommitted GUI code for `ff` in 6 files (`gui/include/ThmExperimentsPage.h`,
  `ThmSettings.h`, `gui/src/ThmExperimentsPage.cpp`, `ThmSettings.cpp`, `ThmWorkspace.cpp`,
  `tests/gui/thm_workspace_test.cpp`): written, **not compiled** (the machine had < 2.7 GB free). Build `AZURE2` and
  `thm_workspace_test` in the GUI build, run the test, then commit in the style of 08ccea3. Note: `origin/thm` has
  since merged dev (c6ddb741); rebase or merge accordingly.
- `fix/edep-cow-acquire` (base `origin/dev` e8727574): 3fba1192, acquire fence on the copy-on-write ownership test of
  the shared point tables (`SoleOwner()` in `src/EPoint.cpp`) and the "a published table is never written again"
  invariant in `include/NodeTableCache.h`. Checked: race test passes at 4 threads; ctest 57/59 and physics 90/96 with
  exactly the same failures on unpatched dev (memory limit of that machine, plus `component_mapping` (c),
  `total_optical` (a), `pyazr_failed_point`, `pyazr_composite_terms`, `shared_point_tables` (d) — worth checking on a
  bigger machine; they may be real dev regressions). Not checked: clang TSan, macOS.
- `origin/thm` lacks the use-after-free fix a6bc6b9 that is on `dev`; merge dev into thm.

**Not in any repository.** `~/THM/` on the analysis machine (THM notes, evaluations, Marco's review in
`~/THM/review/marco_2026-11/`: `MARCO_REVIEW.md`, `IMPLEMENTATION_PLAN.md`, the formula sheet, meeting notes; the
15N D3 refit in `~/THM/evaluations/15N_p_a0/review2026/`). Copy it separately if needed (private 24Mg data must stay
out of any repository).

**Decisions open (Jakub, then James deBoer's group).**
1. Agree the structure and the phase order; who owns which phase.
2. The Phase 0 freeze on perf work in the point model and on new Qt features.
3. Retiring Qt for the web UI after one release with both. **Decided 2026-11-07 (Jakub): no Qt GUI in the upgrade.**
   Builds and checks run with `-DBUILD_GUI=OFF`; `gui/` and `tests/gui/` are removed on their own branch; the web UI
   (Phase 9) is the only graphical front end. Qt stays on `dev`/`qt5` for released versions.
4. The quadrature tolerance rule for Phase 5.
5. Whether pyazr keeps a pure-Python `.azr` reader after Phase 3 (plan: no).
6. pycoulomb licence and the origin of `coulsm` (Jakub, owner). 2026-11-07: Jakub says to take pycoulomb from
   github.com/skowrons94/pycoulomb (his); cloned at 9fee957. Still open: its LICENSE and the origin and terms of
   `coulsm` (SAMMY) and `coulcc` (CPC program library), before either enters this repository. **Decided 2026-11-07
   (Jakub): pycoulomb's C++ sources may be taken directly into AZURE2; permission for them obtained (incl. from Carl).
   Their purpose here is to switch between Coulomb algorithms.** The Fortran `coulcc` is not taken (C++ only).
   Checked by diff: pycoulomb `coulsm/` is AZURE2's present `external/coul` (the SAMMY-side patch, `cwfcomp_accurate`,
   in AZURE2 since 17fd9ded, 2025-09-05) up to a namespace and comments; pycoulomb `coul/` is AZURE2's coul before
   that patch. The default algorithm is therefore `coulsm`, and it must reproduce today's numbers bit for bit.
7. Redistribution terms of any stopping table before shipping one.
8. (Added 2026-11-07.) The criterion of `adaptive_grid_reference` for a level crossing the window edge: thm's grid
   overcounts a level's area by 1.69 × 10⁻⁴ by design, so the swing equals that offset and exceeds the 10⁻⁴
   tolerance. Options (a) test the variation beyond the net change, (b) finer tail steps (moves every fold by about
   10⁻⁴, regenerates the goldens), (c) take dev's grid. Analysis in the `fix/adaptive-grid-edge-swing` entry.
   [Proposed: (a).] **Decided 2026-11-07 (Jakub): (a).** **Implemented by 47c119e4** (`fix/adaptive-grid-criterion`): variation
   beyond the net change, tolerance √3/(2π ppw²) + 10⁻⁴.
9. (Added 2026-11-07.) The nightly sanitizer schedule fires only from the default branch (`qt5`). [Keep the file on
   `feature/upgrade`; push trigger and manual dispatch until the PI decides.]
10. (Added 2026-11-07.) Formal basis: B_c from the lowest level instead of the first listed, and each formal
   level's largest |γ| positive. Both exist as `Convention` options, off (`level_parameters.md`). [Turn both on when
   the engine switches to `level_parameters`, with a formal-basis golden re-baseline.] **Decided 2026-11-07 (Jakub):
   as proposed.** Implemented 2026-11-07 (`feature/upgrade-1-cnuc-split-1b`) with B_c from the lowest level; the sign
   rule as built is the orientation of the Brune level, because "largest |γ| positive" changed external-capture
   observables (step-1b entry, Found 1). **Confirmed 2026-11-07 (Jakub): keep the Brune orientation.**
11. (Added and decided 2026-11-07, Jakub.) Exact energies by default: observables are calculated at the data
   energies, shifted energies and fold nodes, as AZURE2 does today; interpolation on a fixed mesh only as an opt-in,
   checked mode (`ARCHITECTURE.md` §1, §3, §4).
12. (Added and decided 2026-11-07, Jakub.) Data must be easy to adapt (sums, ratios, yields, E1/E2 parts, …): a
   measurement's prediction is a model expression of observable selectors and combining nodes (`ARCHITECTURE.md`
   §3 experiment/, §9). Designed in Phase 4 with `io/Project` reading today's segment types into it.

**Next actions.**
1. Push `feature/upgrade` (and, if wanted, `thm-review` and `fix/edep-cow-acquire`) from a machine with GitHub access.
2. Phase 0, first tasks: T1–T3 generators and projects, golden outputs, identity tests, benchmark script.
3. pycoulomb prerequisites (Phase 1, first box) can start in parallel in the pycoulomb repository.
4. Build and commit the `ff` GUI code on `thm-review`; the THM re-analyses (D1–D6 of the Marco review plan) continue
   on `thm-review` until Phase 5.
