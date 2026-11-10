# AZURE2 upgrade: roadmap

Each phase is done on its own branch cut from `feature/upgrade` (named `feature/upgrade-<phase>-<topic>`), merged back
into `feature/upgrade` with `--no-ff` when its gate passes. `feature/upgrade` merges `origin/dev` regularly (merge, not
rebase). Nothing reaches `dev` until James deBoer's group has agreed the direction (`HANDOFF.md`, decisions).

**Gate for every phase.** CI green on Linux, macOS and Windows; every reference χ² and output unchanged — or, where a
phase deliberately changes a grid, different within a quadrature tolerance stated per test in the commit; the T2
benchmark (time, peak memory, 1 and 8 threads) not worse; TSan and ASan clean at 8 threads on T1–T3; the code the phase
replaces is deleted in the same phase. Update `HANDOFF.md` at the end of every work session.

Tick boxes as work lands; put the commit hash after each tick.

## Phase 0 — Stabilise, measure, test projects, rules, skills

- [ ] Freeze: no perf work on `EData`/`EPoint`/`ESegment`; no new Qt features (announce in `dev`'s CLAUDE.md).
      Drafted on branch `docs/upgrade-freeze` (a7345b80, cut from `origin/dev` aab1de0d); not merged: decision 2 in
      `HANDOFF.md` is open.
- [x] Test projects T1, T2, T3 (`ARCHITECTURE.md` §10) under `tests/upgrade/T1..T3/`: generator script with the true
      parameters, synthetic data with noise, `.azr`, README stating what each covers and the truth.
      T1 and T3 done: 63c08471 (branch `feature/upgrade-0-testprojects-t1t3`); recovery check `tests/upgrade/recover.py`.
      T2 done (branch `feature/upgrade-0-testprojects-t2`): 3184b1ee, 10ec535a, dc898924, 52340da3, 425fbc8c
      (generator, data, fixes), 7ece1045 (README); recovery check 345bbb69, 599e6bee
      (`tests/upgrade/T2/check_recovery.py`): max |pull| 2.15 (noisy), 0.38 (Asimov). Not seconds-scale: the T2
      fit takes 585 s at 8 threads (72 s with a Coulomb table); its calculation 1.8 s.
- [x] Golden outputs of today's engine for T1–T3 (`tests/upgrade/golden/`), with the engine commit recorded.
      T1 and T3: 3b48a101 (engine c1cc92c7, `golden/golden.py`, rtol 1e-6). T2: 8870504f (`golden/T2`, engine
      sources 37978592, `golden/T2/regen_compare.py`, rtol 1e-5 calculations, 1e-3 fit).
- [x] Identity checks (no data) as unit tests: list in `ARCHITECTURE.md` §10. fbad9b6a, bf5f1713: every identity
      the present engine can express (Coulomb/mpmath/Wronskian, unitarity, optical theorem, Breit–Wigner area,
      Brune/formal/Park, gradient, threads). Pending on Phases 1, 5, 7 (kernels, stopping, THM, conventions): listed in
      `tests/upgrade/README.md`, each to be added with its module.
- [x] Benchmark: `scripts/bench_upgrade.sh` → time and peak RSS for T2 at 1 and 8 threads; stored per commit.
      345bbb69, 5a9426fd (script), f51f8d18 (first entry, `tests/upgrade/bench/T2.tsv`).
- [x] CI: nightly clang TSan and ASan at 8 threads on T1–T3 and the mapping/shift/target tests.
      Branch `feature/upgrade-0-sanitizers`: 557bb0ba, 16309cf8 (`AZURE2_SANITIZE`, TSan ignorelist), cfd281ab
      (`scripts/sanitize_upgrade.sh`, suppressions), 8572b21c (`.github/workflows/upgrade-sanitizers.yml`), 635b0a1a
      (`sanitizer_findings.md`). ASan clean (clang, GCC); TSan reports two AZURE2 sites, §4.1 (iteration counter,
      real) and §4.2 (copy-on-write ownership test, TSan-blind fence), to be fixed on `fix/` branches; the TSan job is
      red until then. The workflow has not run on GitHub; the nightly fires only once the file is on the default
      branch.
- [x] Rules into `CLAUDE.md` (`RULES.md`), replacing the "date an entry (`dev 0dd28e5`)" skill rule.
- [x] Skills rewritten into the ten-skill set (`RULES.md` §3) and `tools/check_skills.py` in CI.
- [x] Merge the race fixes everywhere: `a6bc6b9e` (via origin/thm c6ddb741, merged d20718f8) and the acquire-fence
      fix 3fba1192 (`fix/edep-cow-acquire`, merged 37978592) are both on `feature/upgrade`. Pushing them to `thm` and
      `dev` is the PI's.
Gate: T1–T3 run in seconds and recover their truth with today's engine; golden files committed.
Gate status (2026-11-07): met except (a) the T2 fit, which takes minutes, not seconds (its calculation, 3 s, is in
`tests/upgrade/check.sh`; the fit is for nightly CI), (b) [fixed: §4.1 and §4.2 merged a2392e40, f8b7bf22; TSan 0
reports] and (c) [fixed: `adaptive_grid_reference`, inherited from origin/thm c6ddb741, takes decision 8 (a)
in 47c119e4, `fix/adaptive-grid-criterion`]; the other four inherited failures are fixed (merges 05e08e06, 1d35e1d4, 22214e14;
`ctest -LE slow` 130/131). The freeze awaits decision 2.

## Phase 1 — Physics core as arrays; pycoulomb; conventions

- [ ] pycoulomb (owner Jakub) prerequisites: LICENSE file (GPL-3 compatible; state origin and terms of the
      SAMMY-derived `coulsm`); CMake library target and install without Python; stateless, thread-safe
      `coulomb(l, eta, rho) -> F, G, Fp, Gp, sigma_l`; Whittaker functions (move AZURE2's `WhitFunc`); remove committed
      build products (`.so`, `libcoulfg.a`, `cmake_install.cmake`, `.DS_Store`); tests against mpmath; tag a release.
- [ ] `physics/` library target. `System`, `JpiBlock` (E[λ], γ[λ, c] arrays), partitions, channels.
      First step done on `feature/upgrade-1-physics-system` (merged): the `azure2_physics` target, `System` and
      `JpiBlock` as values (`physics_system.md`), the `CNuc` → System adapter `src/CNucSystem.cpp`, tests
      `physics_system_unit`, `physics_standalone`, `physics_system_from_cnuc`. The engine does not use them yet.
- [ ] `channel_functions(α, l, E)` from pycoulomb, tables per (α, l, E); replaces `external/coul`, `CoulFunc`,
      `CoulFuncCache`, `WhitFunc`. Cross-method identity test (coul, coulsm, coulfg, GSL).
      *First part done on `feature/upgrade-1-coulomb` (2026-11-07, merged; design `coulomb.md`):
      `physics/coulomb` (coulsm, coul, coulfg, coulcc; mpmath, Wronskian, identity, threads; GSL cross-check),
      `external/coul` deleted, open-channel `channel_functions` and the kinematics in physics/ called by the engine;
      bit-identical outputs except hybrid channels (≤ 9.9 × 10⁻¹¹, explained). `CoulFuncCache` kept: its deletion
      (branch `-nomemo`) costs 3–17× and waits for the per-stage table. Open (step 3b): Whittaker functions without
      GSL, the per-stage table, then the memo; hybrid channels, `CoulombTable`, `--gsl-coul`.*
      *Step 3b (ii) on `feature/upgrade-1-channel-table` (2026-11-07, merged; design `channel_table.md`): the
      per-stage channel table and the radial scratch, `CoulFuncCache` deleted; outputs as on `-nomemo`, THM checks
      and thread identity deterministic; free-shift refills from the segment's last three shift tables (interim until
      the analytic shift derivative); T2 fit200 −13 % / −23 % (1 / 8 threads), 13N +1.5 % (within spread),
      fit peak RSS +13 % / +30 % (`channel_table.md` §6).*
- [ ] Split `CNuc`: nucleus data → `System`; boundary conditions, shift functions, Brune/Park transforms, angular
      coefficients → functions of the System.
      *Step 1 done on `feature/upgrade-1-cnuc-split` (2026-11-07, merged): plan `cnuc_split.md`;
      `CalcBoundaryConditions`, `TransformIn` (all bases) and `TransformOut`'s formal iteration call `level_parameters`;
      1503 output files bit-identical in three bases.*
      *Step 1b done on `feature/upgrade-1-cnuc-split-1b` (2026-11-07, merged): B_c from the lowest level, formal
      levels oriented like their Brune levels (not "largest |γ| positive", which changed capture observables), own
      eigensolvers, one channel energy; every changed output of 179 projects × 3 bases classified, no re-baseline.*
      *Step 2 done on `feature/upgrade-1-cnuc-split-2` (2026-11-08, merged): `observed_levels` is the one
      definition of Γ, C, μ, Q, B from γ^int + γ^ext (five copies deleted); `ConvertAmplitudeBasis` through
      `park_from_brune`; external capture takes the final state's Brune amplitude and energy in every basis (formal ↔
      Brune capture to `beam_profile_kernel`'s 6.049 MeV state 0.81 → agree; `tests/beam_profile_kernel/check.sh`);
      Brune and Park outputs identical, formal changes classified.*
      *Step 4 done on `feature/upgrade-1-external-capture` (2026-11-08, merged): `physics/external_capture` (N_f,
      γ^ext, the multipole rule, the projector Π); `--ignore-externals` on the input, the same in every basis; one
      adjoint test for the formal-basis parts; Brune and Park outputs identical, formal changes classified. Open before
      merge: `brune_from_formal`'s 10⁻¹² MeV tolerance at the rounding limit (`external_capture.md` §8). Next: 3b (i),
      5–7.*
      *Lectures (PI, 2026-11-08): an R-matrix lecture series in Notion-ready Markdown, later slides, from the code's
      own definitions (`docs/lectures/PLAN.md`); lecture 1 (external region) has its figures, text to write.*
- [ ] Internal convention written in `docs/source/theory/conventions.rst`; `Convention` record; recoupling (6j) and
      phase conversions in `level_parameters`; convention-invariance test on T1–T3.
      *Documentation part done on `feature/upgrade-1-conventions` (2026-11-07, merged): the page states the
      present engine's conventions with file:line references, lists where existing docs contradict the code, and
      sketches the record (planned); `tests/coulomb_nuclear_sign` pins the Coulomb–nuclear interference sign. Open:
      the record, `level_parameters`, the invariance test, and an independent check of the internal–external
      capture sign.*
      *Code part done on `feature/upgrade-1-level-parameters` (2026-11-07, merged): the `Convention` record in the
      System; `level_parameters` (observed ↔ Brune/Park, Brune ↔ formal, B_c, external-capture subtraction, LS ↔ jj
      by 6j, order and iˡ signs, an orthogonal T with its Jacobian), own 3j/6j and eigensolvers
      (`level_parameters.md`). Against `CNuc` on T1–T3: Brune and Park amplitudes and B_c bit-identical, formal to
      1e-14; invariance in six conventions at the parameter level. Open: the engine does not call it (with the `CNuc`
      split), no `.azr` keyword states a convention, the binding, the internal–external capture sign.*
Gate: bit-identical outputs (except the deliberate convention test cases); Coulomb functions agree with mpmath.

## Phase 2 — One objective

- [ ] `inference/` target: `ParameterVector`, `Objective` with one list of terms (value, gradient, residuals,
      Jacobian); every penalty of `EData`, `CNuc`, `ESegment` and pyazr's `penalties()` becomes a term.
- [ ] CLI fit, gradient, API, output writer all call it; the MCMC path keeps working until numcmc is removed.
Gate: bit-identical χ², gradients = finite differences; pyazr no longer computes penalties in Python.
*Done on `feature/upgrade-2-objective` (2026-11-09, merged; `objective.md`): one ordered list of terms (`Objective`)
for the FCN, the LM/FD cost, the MCMC likelihood and the API; one enumeration of the diagonal priors for gradients,
G2, least-squares rows and the MCMC; pyazr computes no prior (2.11.0). χ² bit-identical where no decision applies;
PI decisions D1–D7 and "the MCMC shares the fit's configuration" applied, each change classified. Open: the writer
from the terms rather than stored totals; `inference/` as a separate library (after Phase 4).*

## Phase 3 — One project model; library readers

- [x] `io/Project`: one `.azr` reader/writer; engine and pyazr read through it; old files keep loading.
- [x] Round trip of the archive of 48,922 `.azr` files (the sweep used for the parser hardening). *P1: the archive is
      on another machine; done on every `.azr` of the repository in 945 variants (layouts, line endings, removals).*
- [ ] pyazr `sources/`: AME/NUBASE, ENSDF, EXFOR (X4 master files), IAEA measured stopping; provenance records;
      fixtures of a few kB per library; `pyazr libraries update` as the only networked command.
Gate: byte-identical round trips; one mass table in the whole code.
*Done on `feature/upgrade-3-project` (2026-11-09, merged; `project.md`): `io/` (document, format 2 and `--convert`,
level grouping); the engine and pyazr read every section through it (pyazr's `<thm>` port deleted); unedited files
written back byte for byte, edits touch only their lines; the engine's THM nuclei generated from the one AME2020
file; u, ħc, α, m_e defined once; `CLAUDE.md` recipes with consistency checks 10–12. 179 × 3 identical except the
format-1 note in run.log. Open: `sources/` (the item above), `<scatteringLengths>` tokens still checked by pyazr.*

## Phase 4 — Data are data

Design: `data_model.md` (survey, types, decisions M1–M13, steps S1–S16), `data_kinds.md` (what the engine accepts,
what EXFOR holds, the proposal for named quantities, targets, Rutherford ratios, resonance strengths),
`project_format3.md` (the format the engine reads), `fit_threads.md` (the fit's dependence on the thread count).
Started 2026-11-09; integration branch `feature/upgrade-4-data`.
Done by 2026-11-10 (merged, verified): S3 (a) `Results` behind the old accessors and `AZURE_POISON_RESULTS`; S6 frames;
S8 expression nodes; the first cause of thread-dependent fits (S2, exact grid comparison); format 3 step 2 (record
reader, the table, `--schema`); S3 (b) the buffer owned by the caller; parameter kinds as records
(`parameter_kinds.md`); format 3 step 3a (`convert/`: the converter old → format 3, `--validate`); S1 (history
test), the S2 remainder and S4 (one model pass and one residual, `src/ModelPass.cpp`, check 28); S15 (`<modelData>`:
ωγ, resonance energy, widths, branching ratio, ANC, scattering length; `ModelDataTerm`, check 27); S9 and format 3
step 3b-1 as far as the specs (the engine's set-up consumes specs; `Schema::Build` makes them from format-3
records; spec-equivalence test, 183 of 186 converted files equal, 3 recorded). In work: S10 and S11 (the model pass
evaluates the expression tree, the gradient calls its `backward`); S5 (the calculation plan) with the follow-ups of
the model pass; the engine runs a format-3 file, `--describe`, the converted-tree gate.

- [ ] State of an evaluation in one place (S1–S5): a `Results` buffer owned by the caller, with slots per calculation
      unit; one model pass and one residual function; the mapping as a plan; no writes from one object into another and
      nothing carried from one evaluation to the next. `AZURE_POISON_RESULTS=1` fails on any read before write.
- [ ] `experiment/` target (S6–S8), a standalone library that links only `physics/`: frames (lab ↔ c.m., one
      implementation), `Measurement` (immutable, shared between threads) and `MeasurementFrame` (what depends on fitted
      shifts), and the model-expression nodes `Leaf`, `Sum`, `Scale`, `Weighted`, `Ratio`, `Transform` with `value` and
      `backward`.
- [ ] Model expressions in the engine (S9–S13): one keyword → tree table; the tree pass replaces the 19 loops and the
      two compositions (forward in `ESegment`, adjoint in `AZUREGrad`); components as leaves; amplitude filters
      (multipole, Jπ, l, s, internal/external) on any cross section. Test: every test project's χ² unchanged; an
      E1-only and an E2-only fit of T1 against a closed form.
- [ ] Format 3 (maintainer's decision 2026-11-09: no backward compatibility in the engine, a converter instead): a data
      set states its quantity, unit, frame, reference (Rutherford or Mott), target by name and normalisation in words;
      `<targets>` with names replaces `<targetInt>` index lists; `<modelData>` holds model-level data (resonance
      strengths ωγ, scattering lengths, widths). Old files are read only by `AZURE2 --convert` / `pyazr.convert_azr`.
      One table of quantities feeds the file validation, Python, the GUI and the documentation; `AZURE2 --schema` and
      `--validate` give them to programs and AI agents.
- [ ] New data kinds, in this order (`data_kinds.md` §5.3): ratio to Rutherford and to Mott; S factor as an engine
      quantity (with it, the fix of pyazr's S-factor conversion, which uses mass numbers: σ +2.9 % for ³He+α at 93 keV,
      +6.7 % for p+¹²C at 100 keV); named targets and thick-target yields per incident particle or per charge;
      ωγ; Legendre coefficients with their normalisations. Extrapolation segments on a list, a file or an adaptive grid.
- [ ] Recipe "a new kind of data" in `CLAUDE.md` and its checks (S16).
Gate: outputs of the 179 projects × 3 bases bit-identical through the restructuring (S1–S13); for format 3, the outputs
of the converted projects identical to the originals' on the previous engine.

**`CNuc` split, to finish** (own track; Phase 1 left it at steps 1–5 and 6a): `src/CNuc.cpp` was 2,907 lines and is
1,954 after the moves of 2026-11-10 (printing, checks, transformation messages, the `<levels>` grammar to `io/`).
Remaining (`cnuc_split.md`, "Remaining steps"): 3b (iii) hybrid-potential channels, 3b (iv) `CoulombTable` /
`--gsl-coul`, 6b consumers of the pathway trees, 7 the parameter vector and the set-up out, `CNuc` deleted. The
file-size ratchet is consistency check 26 (`tests/consistency/file_sizes.tsv`): no listed file grows, no other file
passes 1,000 lines, and the record only goes down.

**`EData`, `EPoint`, `ESegment` split** (own track beside the `CNuc` split; plan `edata_split.md`): at 070756a1
`src/EData.cpp` has 5,571 lines, `EDataThm.cpp` 1,233, `EPoint.cpp` 3,006, `ESegment.cpp` 1,446. Each step is a
verbatim move, a split in place or a deletion, bit-identical on 179 × 3. Now, beside the running branches: the output
writer (490 lines), the check files (240) and the THM reports and tables (461) into files of their own, and the writer
in four functions (N1–N4). After the current wave: dead code (W1) and the data covariance (W2). After 3d, S7 and S12:
the load-time checks, the parameter layout, the hand-over and the shift code into their own files; S5, S7, S10–S12
delete mapping, frames and composition. The parts Phases 5 and 6 rewrite (target effects and grids, external capture,
folds, channel tables) wait for the director's decision D1 (`edata_split.md` §8).

## Phase 5 — Energy grid, kernels, stopping, THM observable

- [ ] `EnergyGrid` on compound-nucleus energies; one collision-matrix solve per grid energy for all reactions; grid
      fixed per fit stage, planned from data, windows and level widths.
- [ ] Observables O_n: integrated σ, Legendre coefficients, nuclear amplitudes with f_C(θ) in kernel rows, A_y as a
      separate numerator, phase shifts, β-delayed spectra, optical theorem, identical particles.
- [ ] Kernel builders: target, resolution, beam profile, UDR, Doppler, energy shifts (interpolation weights), SUM,
      RATIO, natural targets, S-factor diagonal; analytic ∂K/∂ν.
- [ ] Lab data with a fitted shift compared in the frame of the data (decision M3 (b) of `data_model.md`, taken by
      the PI 2026-11-10): the model goes to the lab, the measured values and errors stay fixed. Every fit it moves
      is listed with its size.
- [ ] Stopping module (`ARCHITECTURE.md` §7): `Material`, ZBL nuclear, Bethe with corrections, straggling with
      Tschalär propagation, user tables; then the fitted low-energy evaluation for p and α (IAEA data), gases; then
      compounds and heavier ions; ERYA and `SRIM2013.xml` removed; licence check of any shipped table.
- [ ] THM (`THM.md` §1–2): σ_HOES observable with plane-wave and full-form-factor vertices, THM kernel (c.m.
      resolution, QF window, OES S-fold), free or double-ratio normalization, legacy K options marked.
Gate: outputs within the stated quadrature tolerance, each difference explained; smooth χ² in shifts; gradient = FD.

## Phase 6 — Evaluation without copies

- [ ] Parallel over grid-energy blocks with per-task scratch; explicit thread budget; remove pools, `EData::Clone`,
      point mapping, copy-on-write tables, `NodeTableCache`, process-wide memos; EC unit amplitudes per grid energy;
      rates and extrapolation on the grid (+ analytic narrow-resonance term).
- [ ] No hidden state (`ARCHITECTURE.md` §9): the mutable function-level statics of `src/` (about 50), the
      `NuclearPotentialManager` singleton and the `thread_local` caches of `ECIntegral.cpp` become caller-owned tables
      or go; a consistency check refuses new ones.
Gate: thread-count independence by construction; memory one problem regardless of threads; T2 benchmark recorded.

## Phase 7 — DWBA for THM vertices (leaf module; may start in parallel with Phases 4–6)

- [ ] `physics/potentials` (Woods–Saxon, Coulomb, global optical sets) and `physics/transfer` (bound states, level
      functions, distorted waves via pycoulomb's Numerov, zero-range then finite-range post-form DWBA vertex).
- [ ] Validation: plane-wave limit to 1e-6; FRESCO benchmarks for d, ³He, ⁶Li, ¹⁴N Trojan horses to 1 %.
Gate: benchmarks committed with their FRESCO inputs and outputs.

## Phase 8 — pyazr 3.0

- [ ] Object binding (`ARCHITECTURE.md` §11): modules, numpy views, `physics` explorer functions, `sources`, `cli`.
- [ ] Wheels via scikit-build-core + cibuildwheel (Linux, macOS, Windows) with the core and the `azure2` CLI.
- [ ] `azure2()` session shim with deprecation warnings for one release; docs recipes for iminuit, scipy, emcee.
Gate: ported pyazr tests pass; wheels install and run T1–T3 on three platforms.

## Phase 9 — Web UI; retire Qt, ERYA, numcmc

- [ ] `pyazr ui`: FastAPI + uvicorn, static page, plotly.js vendored; `/api` one-to-one with pyazr; jobs with progress
      and stop; tabs Project, Nucleus, Data, Fit, Results, Physics, THM (opt-in).
- [ ] Physics tab plots every exposed physics function (`ARCHITECTURE.md` §8, `THM.md` §4).
- [ ] pytest on `/api` with T1–T3; browser smoke tests; every item of the Qt parity checklist
      (`docs/upgrade/qt_parity.md`) ticked or dropped with a written reason.
- [x] Qt parity checklist written from the Qt code before its removal (`docs/upgrade/qt_parity.md`, 201 items and 32
      test lines, from `4bf7aa06:gui/` and `tests/gui/`).
- [x] `gui/`, `tests/gui/`, `BUILD_GUI`/`USE_QWT`/`BUILD_MACOS_BUNDLE`, Qt/Qwt in CMake, CI, Docker and build scripts
      removed (decision 3, 2026-11-07; branch `feature/upgrade-remove-qt`); `--no-gui` kept, accepted and ignored.
- [ ] Remove `numcmc/`, `external/erya`, readline menus, the GSL fitter.
Gate: parity checklist signed off by Jakub; release notes.
