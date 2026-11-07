# AZURE2 upgrade: handoff

Update this file at the end of every work session: what changed, what is next, what is open. Newest first.

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
3. Retiring Qt for the web UI after one release with both.
4. The quadrature tolerance rule for Phase 5.
5. Whether pyazr keeps a pure-Python `.azr` reader after Phase 3 (plan: no).
6. pycoulomb licence and the origin of `coulsm` (Jakub, owner).
7. Redistribution terms of any stopping table before shipping one.
8. (Added 2026-11-07.) The criterion of `adaptive_grid_reference` for a level crossing the window edge: thm's grid
   overcounts a level's area by 1.69 × 10⁻⁴ by design, so the swing equals that offset and exceeds the 10⁻⁴
   tolerance. Options (a) test the variation beyond the net change, (b) finer tail steps (moves every fold by about
   10⁻⁴, regenerates the goldens), (c) take dev's grid. Analysis in the `fix/adaptive-grid-edge-swing` entry.
   [Proposed: (a).]
9. (Added 2026-11-07.) The nightly sanitizer schedule fires only from the default branch (`qt5`). [Keep the file on
   `feature/upgrade`; push trigger and manual dispatch until the PI decides.]
10. (Added 2026-11-07.) Formal basis: B_c from the lowest level instead of the first listed, and each formal
   level's largest |γ| positive. Both exist as `Convention` options, off (`level_parameters.md`). [Turn both on when
   the engine switches to `level_parameters`, with a formal-basis golden re-baseline.]

**Next actions.**
1. Push `feature/upgrade` (and, if wanted, `thm-review` and `fix/edep-cow-acquire`) from a machine with GitHub access.
2. Phase 0, first tasks: T1–T3 generators and projects, golden outputs, identity tests, benchmark script.
3. pycoulomb prerequisites (Phase 1, first box) can start in parallel in the pycoulomb repository.
4. Build and commit the `ff` GUI code on `thm-review`; the THM re-analyses (D1–D6 of the Marco review plan) continue
   on `thm-review` until Phase 5.
