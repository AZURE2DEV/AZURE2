# AZURE2 upgrade: handoff

Update this file at the end of every work session: what changed, what is next, what is open. Newest first.

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

**Next actions.**
1. Push `feature/upgrade` (and, if wanted, `thm-review` and `fix/edep-cow-acquire`) from a machine with GitHub access.
2. Phase 0, first tasks: T1–T3 generators and projects, golden outputs, identity tests, benchmark script.
3. pycoulomb prerequisites (Phase 1, first box) can start in parallel in the pycoulomb repository.
4. Build and commit the `ff` GUI code on `thm-review`; the THM re-analyses (D1–D6 of the Marco review plan) continue
   on `thm-review` until Phase 5.
