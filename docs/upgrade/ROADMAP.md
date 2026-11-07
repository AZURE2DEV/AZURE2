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
- [ ] CI: nightly clang TSan and ASan at 8 threads on T1–T3 and the mapping/shift/target tests.
- [x] Rules into `CLAUDE.md` (`RULES.md`), replacing the "date an entry (`dev 0dd28e5`)" skill rule.
- [x] Skills rewritten into the ten-skill set (`RULES.md` §3) and `tools/check_skills.py` in CI.
- [ ] Merge the race fixes everywhere: `a6bc6b9` is on `dev` only (missing on `thm`); the acquire-fence fix
      (local branch `fix/edep-cow-acquire`, 3fba1192, analysis machine) still to push and merge.
Gate: T1–T3 run in seconds and recover their truth with today's engine; golden files committed.

## Phase 1 — Physics core as arrays; pycoulomb; conventions

- [ ] pycoulomb (owner Jakub) prerequisites: LICENSE file (GPL-3 compatible; state origin and terms of the
      SAMMY-derived `coulsm`); CMake library target and install without Python; stateless, thread-safe
      `coulomb(l, eta, rho) -> F, G, Fp, Gp, sigma_l`; Whittaker functions (move AZURE2's `WhitFunc`); remove committed
      build products (`.so`, `libcoulfg.a`, `cmake_install.cmake`, `.DS_Store`); tests against mpmath; tag a release.
- [ ] `physics/` library target. `System`, `JpiBlock` (E[λ], γ[λ, c] arrays), partitions, channels.
- [ ] `channel_functions(α, l, E)` from pycoulomb, tables per (α, l, E); replaces `external/coul`, `CoulFunc`,
      `CoulFuncCache`, `WhitFunc`. Cross-method identity test (coul, coulsm, coulfg, GSL).
- [ ] Split `CNuc`: nucleus data → `System`; boundary conditions, shift functions, Brune/Park transforms, angular
      coefficients → functions of the System.
- [ ] Internal convention written in `docs/source/theory/conventions.rst`; `Convention` record; recoupling (6j) and
      phase conversions in `level_parameters`; convention-invariance test on T1–T3.
Gate: bit-identical outputs (except the deliberate convention test cases); Coulomb functions agree with mpmath.

## Phase 2 — One objective

- [ ] `inference/` target: `ParameterVector`, `Objective` with one list of terms (value, gradient, residuals,
      Jacobian); every penalty of `EData`, `CNuc`, `ESegment` and pyazr's `penalties()` becomes a term.
- [ ] CLI fit, gradient, API, output writer all call it; the MCMC path keeps working until numcmc is removed.
Gate: bit-identical χ², gradients = finite differences; pyazr no longer computes penalties in Python.

## Phase 3 — One project model; library readers

- [ ] `io/Project`: one `.azr` reader/writer; engine, Qt GUI and pyazr read through it; old files keep loading.
- [ ] Round trip of the archive of 48,922 `.azr` files (the sweep used for the parser hardening).
- [ ] pyazr `sources/`: AME/NUBASE, ENSDF, EXFOR (X4 master files), IAEA measured stopping; provenance records;
      fixtures of a few kB per library; `pyazr libraries update` as the only networked command.
Gate: byte-identical round trips; one mass table in the whole code.

## Phase 4 — Data are data

- [ ] `experiment/` target: `Measurement` arrays converted to c.m. at load; results in per-evaluation arrays, not in
      points; hosts no longer write into peers.
Gate: bit-identical.

## Phase 5 — Energy grid, kernels, stopping, THM observable

- [ ] `EnergyGrid` on compound-nucleus energies; one collision-matrix solve per grid energy for all reactions; grid
      fixed per fit stage, planned from data, windows and level widths.
- [ ] Observables O_n: integrated σ, Legendre coefficients, nuclear amplitudes with f_C(θ) in kernel rows, A_y as a
      separate numerator, phase shifts, β-delayed spectra, optical theorem, identical particles.
- [ ] Kernel builders: target, resolution, beam profile, UDR, Doppler, energy shifts (interpolation weights), SUM,
      RATIO, natural targets, S-factor diagonal; analytic ∂K/∂ν.
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
- [ ] pytest on `/api` with T1–T3; browser smoke tests; Qt parity checklist from today's Qt features and
      the GUI review in `azure2-plot` and `gui/CLAUDE.md`.
- [ ] One release with both UIs; then remove `gui/`, `numcmc/`, `external/erya`, readline menus, the GSL fitter.
Gate: parity checklist signed off by Jakub; release notes.
