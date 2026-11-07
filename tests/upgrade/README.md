# Upgrade tests: synthetic projects, golden outputs, identity checks

This directory holds the tests every phase of the upgrade (`docs/upgrade/ROADMAP.md`) must pass: synthetic R-matrix
projects generated from known parameters, which a fit must recover; the outputs of the present engine on them, which
later builds must reproduce until a phase changes a result on purpose; and the list of identity checks of
`docs/upgrade/ARCHITECTURE.md` §10, the tests that need no data because the quantity is known exactly.

| path | what |
|---|---|
| `T1/` | p + ¹²C → ¹³N: elastic, A_y, phase shifts, capture to two bound states with external capture, S factor, rate, extrapolation (`T1/README.md`) |
| `T2/` | p + ¹⁵N → ¹⁶O through an apparatus (branch `feature/upgrade-0-testprojects-t2`) |
| `T3/` | ⁸Be: α + α identical bosons, p + ⁷Li, n + ⁷Be, β-delayed α, optical theorem, scattering lengths, hybrid potential (`T3/README.md`) |
| `lib/upgrade_project.py` | the project model and `.azr` writer, the CLI runner, the synthetic-data generator |
| `recover.py` | the recovery check: fit the start model, compare with the truth |
| `golden/` | outputs of the present engine and `golden.py` to regenerate and compare (`golden/README.md`) |
| `check.sh` | golden comparison and recovery check; run by `tests/run_tests.sh` |

All scripts need only python3 with numpy (and mpmath for the Coulomb identity). They take the engine from an
argument, `$AZURE2_BIN` or `$AZURE2_BINARY`, else the newest `build*/src/AZURE2` of the checkout, and the pyazr module
from the checkout's `pyazr/`. They work on temporary copies and never write into the source tree.

```bash
tests/upgrade/check.sh build/src/AZURE2                       # golden + recovery, ~50 s on two threads
ctest --test-dir build -R pyazr_upgrade --output-on-failure   # the identity checks below (-L slow: threads)
```

## Identity checks (no data)

Each check compares the engine with a closed form, mpmath or an exact identity, in `tests/pyazr/upgrade_*_test.py`
(registered in `tests/pyazr/CMakeLists.txt`). Results with the `feature/upgrade` build c1cc92c7 (2026-11-06):

| identity (§10) | test | what is compared | result |
|---|---|---|---|
| Coulomb functions vs mpmath, Wronskian | `upgrade_coulomb_identity` | F, G, F′, G′ (d/dρ) at a_c for p + ¹²C, p + ⁷Li, n + ⁷Be, l = 0–4, E_c.m. = 0.05–3 MeV; P, S, δ_hs from them | 4.5 × 10⁻¹⁰; Wronskian 4 × 10⁻¹⁶; P, S 2.4 × 10⁻⁹. Neutral F_l at small ρ, l ≥ 3: 7.9 × 10⁻⁶ (closed-form cancellation; Phase 1 target) |
| U unitary with elastic only | `upgrade_unitarity` | n + ¹⁶O, one open channel: σ_tot (optical theorem) = σ_el, 296 energies | equal to 11 printed digits |
| optical theorem | `upgrade_optical_theorem` | T3: σ_tot(n + ⁷Be) = σ(n,n) + σ(n,p) + σ(n,αα) | 2.6 × 10⁻¹¹ |
| isolated narrow level → Breit–Wigner area | `upgrade_breit_wigner` | one 1 eV ¹²C(p,γ) level: ∫σ dE over ±100 eV and the peak | area 1.1 × 10⁻⁷, peak 7.6 × 10⁻⁶; photon width outside the level denominator (stated in the test) |
| Brune ↔ formal ↔ Park round trip | `upgrade_parametrization_roundtrip` | T1, T3: every data point and extrapolation in the three bases | Park 0, formal 2.2 × 10⁻¹¹ — **except** T1 capture with external capture in the formal basis: up to 54 % (open question, HANDOFF) |
| gradient = finite differences | `upgrade_gradient` | T1, T3 at the truth and the start: adjoint ∇χ² vs Richardson central differences | 1.2 × 10⁻⁸ (T1), 2.7 × 10⁻⁶ (T3) of max \|∇χ²\| |
| thread-count independence | `upgrade_thread_independence` (label `slow`) | T1, T3 menus 1, 3, 5 and a fit, 1 vs 4 threads | equal to 3.6 × 10⁻¹¹; one T1 capture file not bit-identical (Phase 6 target: bit for bit) |

### Pending: they need modules that do not exist yet

These are listed, not stubbed. Where a present-engine test covers part of the identity, it is named; it stays the
oracle until the new module replaces the code it tests.

| identity (§10) | needs | present-engine coverage |
|---|---|---|
| constant σ survives every kernel | Phase 5 kernel builders (target, resolution, profile, UDR, Doppler, shifts, SUM/RATIO) | `tests/reference/beam_profile_reference_test.cpp`, `udr_reference_test.cpp`, `tests/convolution_unify` |
| Gaussian ∘ Gaussian = Gaussian (σ² add) | Phase 5 kernels | `tests/convolution_window` (partial) |
| probability conservation of kernel weights | Phase 5 kernels | none as a stand-alone identity |
| convention invariance | Phase 1 `Convention` record and `level_parameters` conversions | none; T1 README lists the convention round trip as TODO |
| stopping: thin-layer limit | Phase 5 stopping module | `tests/reference/compound_stopping_reference_test.cpp` (SRIM tables, to be removed) |
| straggling: Bohr limit | Phase 5 stopping module | `tests/straggling_bohr` |
| THM Green identity | Phase 5 THM observable (`THM.md`) | the THM tests of `thm-review` (Green identity 2.6 × 10⁻¹⁰, HANDOFF) |
| DWBA vertex vs FRESCO | Phase 7 `physics/transfer` | `tests/reference/thm_dw_vertex_test.cpp` (zero-range against mpmath, not FRESCO) |
