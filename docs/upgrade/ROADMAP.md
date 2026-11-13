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

Design: `data_model.md` (survey, types, decisions), `data_kinds.md` (what the engine accepts, what EXFOR holds),
`project_format3.md` (the format the engine reads), `data_structure.md` (the data classes replaced), `cnuc_split.md`
(the pathway table), `fit_threads.md` (the fit's dependence on the thread count). Started 2026-11-09; integration
branch `feature/upgrade-4-data`. The dated record is `HANDOFF.md`; the summary for the owner is `FOR_THE_OWNER.md`.
State on 2026-11-11:

- [x] State of an evaluation in one place: a `Results` buffer owned by the caller with slots per calculation unit,
      one model pass and one residual function (`src/ModelPass.cpp`), the mapping as a plan
      (`src/CalculationPlan.cpp`), nothing carried from one evaluation to the next (`upgrade_history`,
      `AZURE_POISON_RESULTS=1`).
- [x] `experiment/`, a library that links only `physics/`: frames (lab ↔ c.m., one implementation), `Measurement`
      (what was measured, immutable, held once and shared by every copy and thread), `MeasurementFrame` (what depends
      on fitted shifts), the model-expression nodes with `value` and `backward`.
- [x] One keyword → tree table; the tree pass replaces the composition loops, the gradient calls its `backward`.
      Open: components as leaves, amplitude filters (multipole, Jπ, l, s, internal/external) on any cross section,
      `predict`.
- [x] Format 3 (no backward compatibility in the engine, a converter instead): a data set is a named record that
      states its quantity, unit, frame, reference, target by name and normalisation; 9 record kinds, 21 quantities,
      one table that feeds the file validation, Python and the documentation. `AZURE2 --convert`, `--validate`,
      `--describe`, `--schema` (for programs and AI agents); `pyazr.Project` edits a project as records. The engine
      runs format 3; every converted test project gives the old outputs (the gate, `tools/format3_gate.sh`).
- [ ] The engine reads **only** format 3. Done: all of the test tree and of the examples but the THM projects are
      format 3 (176 of 183 project files, the originals kept as the converter's fixtures; `format3_tree.md`); `--validate` reports the
      data rows, covariance files and final states the run would refuse; pyazr's session reads a format-3 project
      from the engine (`pyazr_on_records.md`). In work: the rest of the tree. Then: old files refused with the
      converter's command, the old reading path and `pyazr/azrfile.py` (3,756 lines) deleted, a second round of
      validation (the level rows, kinematic limits, combinations of target effects).
- [x] New data kinds: ratio to Rutherford and to Mott; the S factor as a quantity; units with dimensions checked;
      named targets and thick-target yields per particle or per charge; model-level data (ωγ, resonance energy,
      widths, branching ratio, ANC, scattering length); Legendre coefficients with their expansions; extrapolations
      on a list, a file, a logarithmic range or an adaptive grid; the spin–orbit term of the hybrid potential.
      The yield of a layer of finite thickness with the user's stopping formula (Phase 5, branch 4,
      `tests/layer_yield`); with the evaluated stopping when the stopping module is wired (Phase 5, branch 5).
- [ ] Recipe "a new kind of data" in `CLAUDE.md` with its checks.
Gate: outputs of all test projects × 3 bases bit-identical through the restructuring; for format 3, the outputs of
the converted projects identical to the originals' but for the classified differences of the gate.

**The data classes replaced** (`data_structure.md`; the PI, 2026-11-10: replace, do not carve). `EData`, `EPoint`,
`ESegment` were 10,665 lines in three files. Milestone 1 (merged): the measurement exists once. Milestone 2 (merged
2026-11-11): a term's settings are one shared immutable record, every point and sub-point is a unit of one flat
vector, the plan names units by index, a copy of the working copy is a plain copy. Milestone 3 (merged
2026-11-11): an immutable `DataSet` held as `shared_ptr<const DataSet>` and a small `WorkingCopy` apart; the unit
is 480 bytes (was 608); the three class names and their files are deleted and consistency check 31 keeps them out.

**`CNuc`, to finish** (`cnuc_split.md`): `src/CNuc.cpp` was 2,907 lines and is 1,774; printing, checks, messages, the
level grammar, external-capture set-up and the bridge to `physics/` are files of their own; the pathway table is
immutable, built at set-up and read by everything outside the data classes; the capture analyzing-power terms are
built at set-up. Remaining: external-capture amplitudes at set-up (with the decay-position defect of the data
classes), the pathway trees deleted, the parameter vector and the set-up out, `CNuc` deleted. The file-size ratchet
(`tests/consistency/file_sizes.tsv`): no listed file grows, no other file passes 1,000 lines.

## Phase 5 — Energy grid, kernels, stopping, THM observable

- [ ] `EnergyGrid` on compound-nucleus energies; one collision-matrix solve per grid energy for all reactions; grid
      fixed per fit stage, planned from data, windows and level widths.  Every folded point on the shared mesh (done
      2026-11-12, `feature/upgrade-5-one-grid`: the uniform grid of a point's own and the key `useAdaptiveGrid` removed).
- [ ] Observables O_n: integrated σ, Legendre coefficients, nuclear amplitudes with f_C(θ) in kernel rows, A_y as a
      separate numerator, phase shifts, β-delayed spectra, optical theorem, identical particles.
- [ ] Kernel builders: target, resolution, beam profile, UDR, Doppler, energy shifts (interpolation weights), SUM,
      RATIO, natural targets, S-factor diagonal; analytic ∂K/∂ν.
- [ ] Lab data with a fitted shift compared in the frame of the data (decision M3 (b) of `data_model.md`, taken by
      the PI 2026-11-10): the model goes to the lab, the measured values and errors stay fixed. Every fit it moves
      is listed with its size.
- [ ] Stopping module (`ARCHITECTURE.md` §7): the leaf library `stopping/` with SRIM-2013's stopping powers
      ported from LunaSRIM (any ion, any element, solid or gas, compound correction; decision of 2026-11-11),
      TRIM and Bohr straggling with Tschalär propagation, user tables (done, `stopping.md`); then wired into the
      engine with the stopping-scale nuisance; ERYA and `SRIM2013.xml` removed (both done, branch `stopping-wire`,
      2026-11-11: `stopping evaluation`, `stopping table`, kind `stopping_scale`, `pyazr.stopping`); the owner decides
      whether SRIM's two data files are shipped.
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
      Status 2026-11-12 (`evaluation.md` §7): milestone 1 done (no `thread_local`, no warn-once static, no `mutable`
      solver buffer); milestone 2 done but the input-file bytes (a decision: `cache_stats()`); milestones 3 (the pooled
      working copies), 4 (`CNuc` pools, fitter and sampler state, MCMC globals) and 5 (diagnostics) pending: 88 rows.
      Status 2026-11-13: milestone 5 done (one table of 39 environment switches, one reader in `support/`, read once;
      the switches that move numbers in `Config::environment`; the counters in `RunDiagnostics`); milestone 3 partly
      (copies when needed); 53 pending rows (milestones 3 and 4), 15 permanent.
Gate: thread-count independence by construction; memory one problem regardless of threads; T2 benchmark recorded.

- [x] A level's radiative decay as Γ_γ, branchings and mixing ratios (`gamma_parametrisation.md`; record
  `gamma_decay`, parameter kind 13; branch `feature/upgrade-8-gamma-branchings`); with external capture on its
  channels (note section 2.7; branch `feature/upgrade-8-gamma-external-capture`). Open: the formal basis (refused), a
  Rose-Brink δ (not offered).

## Phase 7 — DWBA for THM vertices (leaf module; may start in parallel with Phases 4–6)

- [x] `physics/transfer`, first part (`transfer.md`; branch `feature/upgrade-7-transfer-1`): optical potentials with
      explicit numbers, bound states, distorted waves (spin 0, ½, 1), zero-range and finite-range post-form DWBA to a
      bound state, the plane-wave THM vertex.
- [x] `physics/transfer`, second part (`transfer.md` §14–§18; branch `feature/upgrade-7-transfer-2`): scattering
      states and FRESCO's bins of x + A, the transfer to them (finite and zero range), d²σ/dΩ dE and the narrow-bin
      limit, the interior part against the level function, the prior form and the remnants.
- [x] Global optical sets: one table of the library (`GlobalOptical.h`; branch `feature/upgrade-7-thm-wiring`).
- [x] The wiring into the THM observable: the distortion factor R(E) and the waves of the DW vertex call the library
      (`transfer.md` §24); the DW vertex's source projection S_lm(a) is a library function and its quadrature is the
      library's (`SourceProjection.h`, §25; branch `feature/upgrade-7-thm-dw-vertex`).
- [x] Validation, first part: plane-wave limit to 1e-6; FRESCO benchmarks for the d and ⁶Li Trojan horses to a bound
      state (¹⁵N(d,n)¹⁶O, ¹²C(⁶Li,d)¹⁶O: < 0.2 %), committed with inputs and outputs (`physics/tests/data/fresco`).
- [x] Validation, continuum: FRESCO bins for ¹⁵N(d,n) (non-resonant, narrow and wide resonance, finite range: ≤ 0.5 %),
      committed; post = prior with remnants (10⁻³). Not reached: zero range and ⁶Li bins, prior form and remnant
      against FRESCO (FRESCO's own numerics, `transfer.md` §15, §17).
- [x] Spin-orbit in the distorted waves of a transfer, the projectile's iT11 (`transfer.md` §19; branch
  `feature/upgrade-7-transfer-3`): ¹⁵N(d,n) with spin-orbit in both channels against FRESCO 2e-4 (iT11 1e-4).
- [x] Validation: ³He and ¹⁴N Trojan horses to bound states (¹⁸O(³He,d)¹⁹F, g.s. and 5/2⁺, FR and ZR, with and without
  spin-orbit: ≤ 8.9e-4; ¹²C(¹⁴N,d)²⁴Mg, 0⁺ and 2⁺, FR: ≤ 8.3e-4; `transfer.md` §20–§21).
- [ ] Validation: the ⁶Li bins (FRESCO 3.5 and Frescox 7.4 alike lose accuracy there, `transfer.md` §23); the ³He and
  ¹⁴N Trojan horses to the continuum.
Gate: benchmarks committed with their FRESCO inputs and outputs.

## Phase 8 — pyazr 3.0

Design: `pyazr3.md` (the measured state of pyazr, the objects of 3.0, the order of work, the PI's decisions).

- [ ] Object binding (`ARCHITECTURE.md` §11): modules, numpy views, `physics` explorer functions, `sources`, `cli`.
- [ ] Wheels via scikit-build-core + cibuildwheel (Linux, macOS, Windows) with the core and the `azure2` CLI.
- [ ] Tested recipes for iminuit, scipy, emcee (no `azure2()` shim: struck by the PI, `pyazr3.md` section 5, 9).
Gate: ported pyazr tests pass; wheels install and run T1–T3 on three platforms.

## Phase 9 — Web UI; retire Qt, ERYA, numcmc

Design: `web_ui.md` (process model, `/api`, state rule, tabs, order of work, decisions for the PI; a measured spike).

- [ ] `pyazr ui`: FastAPI + uvicorn, static page, plotly.js vendored; `/api` one-to-one with pyazr; jobs with progress
      and stop; tabs Project, Nucleus, Data, Fit, Results, Physics, THM (opt-in).
- [x] The server (`pyazr/ui`: route table, state and history, workers, jobs; 2026-11-13).
- [x] The page: header, tabs Project, Nucleus, Data with generated forms, the engine's descriptions and a preview plot;
      plotly.js 4.1.2 vendored; node tests and a headless-Chrome smoke test (2026-11-13). Fit, Results, Physics, THM
      tabs to come.
- [x] The options of the calculation are `<config>` keys (decision 7; `project_format3.md` §3.12): flag or keyword >
      key > default, sources in `--describe` and `azure2.describe()`, contradictions refused, one table; the server's
      provisional `SESSION_OPTIONS` gone (2026-11-13).
- [ ] Physics tab plots every exposed physics function (`ARCHITECTURE.md` §8, `THM.md` §4).
- [ ] pytest on `/api` with T1–T3; browser smoke tests; every item of the Qt parity checklist
      (`docs/upgrade/qt_parity.md`) ticked or dropped with a written reason.
- [x] Qt parity checklist written from the Qt code before its removal (`docs/upgrade/qt_parity.md`, 201 items and 32
      test lines, from `4bf7aa06:gui/` and `tests/gui/`).
- [x] `gui/`, `tests/gui/`, `BUILD_GUI`/`USE_QWT`/`BUILD_MACOS_BUNDLE`, Qt/Qwt in CMake, CI, Docker and build scripts
      removed (decision 3, 2026-11-07; branch `feature/upgrade-remove-qt`); `--no-gui` kept, accepted and ignored.
- [ ] Remove `numcmc/`, readline menus, the GSL fitter (`external/erya`: removed in Phase 5, branch `stopping-wire`).
Gate: parity checklist signed off by Jakub; release notes.
