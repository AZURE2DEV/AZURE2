# T1 — capture and scattering, p + ¹²C → ¹³N

T1 is a synthetic R-matrix evaluation of the ¹³N compound nucleus between E_x = 2 and 3.6 MeV: proton elastic
scattering, the vector analyzing power A_y, the s-wave phase shift and radiative capture to two bound states, all
generated from the parameters listed below. A fit with the present engine must return those parameters
(`../recover.py`), and the present engine's outputs on it are kept in `../golden/T1/`. Each later phase of the
upgrade must reproduce both. The project runs in about 3 s (calculation) and 3 s (fit) on two threads.

## The system

Two particle partitions enter: p + ¹²C (pair 1) and γ + ¹³N (pairs 2 and 3).

| pair | particles | J^π of the heavy particle | E_x of the heavy particle (MeV) | S (MeV) | a_c (fm) |
|---|---|---|---|---|---|
| 1 | p + ¹²C | 0⁺ | 0 | S_p = 1.94350 | 3.6 |
| 2 | γ + ¹³N | 1/2⁻ | 0 | — | — |
| 3 | γ + ¹³N | 3/2⁻ | 1.000 | — | — |

Atomic masses are AME2020 (m_p = 1.00782503223 u, m(¹²C) = 12 u, m(¹³N) = 13.00573861 u); S_p = m_p + m(¹²C) − m(¹³N)
with 931.49410242 MeV/u. The g factors are g_p = 5.5856946893 and g(¹²C) = 0. The 3/2⁻ state at E_x = 1.000 MeV is
**fictitious**: ¹³N has no bound excited state. It exists so that capture to an excited final state, with an
asymptotic normalization coefficient (ANC) of its own, is part of the test. Both photon pairs allow external capture
in E1, M1 and E2 (multipolarity mask 7 in column 31).

## The true parameters

Level parameters are in the Brune basis (Brune, Phys. Rev. C 66, 044611 (2002)): E_x is the observed excitation
energy of ¹³N (MeV), Γ the observed partial width in the c.m. frame (eV), C the ANC of a bound channel (fm^−1/2). The
boundary condition is the shift function at each level energy. Free parameters (fitted in the recovery check) are
marked *.

| J^π | E_x (MeV) | channel | value |
|---|---|---|---|
| 1/2⁻ | 0 | p + ¹²C, l = 1, s = 1/2 | C = 1.85 * |
| 3/2⁻ | 1.000 | p + ¹²C, l = 1, s = 1/2 | C = 1.20 * |
| 1/2⁺ | 2.3650 * | p, l = 0 | Γ_p = 33.0 keV * |
| | | γ → g.s., E1 | Γ_γ0 = 0.45 eV * |
| | | γ → 3/2⁻ (1.0), E1 | Γ_γ1 = 0.10 eV * |
| 3/2⁻ | 3.5020 * | p, l = 1 | Γ_p = 55.0 keV * |
| | | γ → g.s., M1 | Γ_γ0 = 0.50 eV * |
| | | γ → g.s., E2 | 0 (E2 enters through external capture only) |
| | | γ → 3/2⁻ (1.0), M1 | Γ_γ1 = 0.20 eV * |
| 5/2⁺ | 3.5470 * | p, l = 2 | Γ_p = 47.0 keV * |
| | | γ → 3/2⁻ (1.0), E1 | Γ_γ1 = 0.050 eV * |
| 1/2⁺ | 15.0 | p, l = 0 | Γ_p = 3.0 MeV (background pole) |
| 3/2⁻ | 15.0 | p, l = 1 | Γ_p = 2.0 MeV (background pole) |

Every nonzero photon width is free. A photon width in the `.azr` is an observed width, and its external part depends
on the ANC of the final state; a fixed nonzero photon width would therefore move with a fitted ANC, and the start
model would not be the truth with only the free parameters moved. The 13 free parameters are the 1/2⁻ and 3/2⁻
ANCs, the energies of the three resonances, their proton widths and the five nonzero photon widths.

`generate.py` holds these numbers (`model()`); `T1.azr` is written from them, and `T1_start.azr` is the same model
with every free energy raised by 4 keV and every free width and ANC multiplied by 1.15.

## The data

The data are the engine's own observables at the true parameters (`upgrade_project.truth_values`: each data file is
first written with value 1, so the c.m. datum the engine reports is the lab-to-c.m. factor of that point) plus
Gaussian noise of the stated standard deviation, from `numpy.random.default_rng(13)`. The uncertainty column is
that standard deviation. All energies are E_lab of the proton (MeV) on a ¹²C target at rest.

| file | observable (isDiff) | frame | E_lab grid (MeV) | angle | noise | points |
|---|---|---|---|---|---|---|
| `el_lab090/125/160` | dσ/dΩ elastic (1) | lab | 0.300–2.000, step 0.025 | θ_lab = 90°, 125°, 160° | 2 % | 3 × 69 |
| `ay_cm090/150` | A_y elastic (7) | c.m. | 1.55–1.90, step 0.01 | θ_c.m. = 90°, 150° | 0.02 absolute | 2 × 36 |
| `phase_s12` | δ(J = 1/2, l = 0) (2) | — | 0.30–2.00, step 0.05, \|δ\| < 75° | — | 0.5° absolute | 35 |
| `cap_g0` | σ(p,γ₀) (0) | — | 0.12–2.00, step 0.04 | — | 5 % | 48 |
| `cap_g1` | σ(p,γ₁) to the 3/2⁻ state (0) | — | 0.15–2.00, step 0.05 | — | 5 % | 38 |
| `cap_g0_lab000/090` | dσ/dΩ(p,γ₀) (1) | lab | 1.40–2.00, step 0.03 | θ_lab = 0°, 90° | 5 % | 2 × 21 |

Phase shifts follow the engine convention −90° < δ < 90°; points with \|δ\| ≥ 75° at the truth are dropped, because
noise on a point near ±90° folds across the branch cut. χ² of the true parameters on these data is 482.3 for
N = 442 points (expected 442 ± 30).

The `<segmentsTest>` lines extrapolate the capture cross sections and S factors to E_lab = 0.02–0.60 MeV, give the
s1/2 and p3/2 phase shifts from 0.30 to 2.00 MeV and the capture angular distribution at 1.5–1.9 MeV. The reaction
rate of ¹²C(p,γ₀) is computed at T₉ = 0.01–3 (`temperatures.dat`).

## What T1 covers, against `docs/upgrade/ARCHITECTURE.md` §10

| item | status |
|---|---|
| elastic dσ/dΩ at lab angles | data, three angles |
| (p,γ) to the ground and an excited state | data (excited state fictitious, see above) |
| external capture E1/E2/M1 | mask 7 on both photon pairs; the resonant E2 width is zero |
| A_y | data, two c.m. angles |
| phase shifts | data (s1/2); test segments (s1/2, p3/2) |
| S factor | column 5 of the `.extrap` and `.out` files; golden outputs |
| reaction rate | menu 5; golden output `rate/reactionrates.out` |
| Brune ↔ formal ↔ Park round trip | identity test `tests/pyazr/upgrade_parametrization_roundtrip_test.py`: Park equal; formal equal except capture with external capture (up to 54 %, open question in `docs/upgrade/HANDOFF.md`) |
| convention round trip | **TODO (Phase 1)**: the `Convention` record and the conversions in `level_parameters` do not exist yet |
| extrapolation | test segments to E_lab = 20 keV |

## Recovery (present engine, 2026-11-06)

`python3 ../recover.py T1` with the build of `feature/upgrade` c1cc92c7, two threads, MIGRAD on the analytic
gradient: valid minimum after 49 calls in 2.7 s; χ² = 462.3 at the fit against 482.3 at the truth. Pulls
(p_fit − p_true)/σ_p, with σ_p the MIGRAD uncertainty of the engine's reduced-width parameter: largest \|pull\| 3.64 for
the 5/2⁺ energy (E_x = 3.54664 ± 0.00010 MeV against 3.5470 MeV), all others ≤ 2.03. On noise-free data the same fit
returns every parameter within 0.007 σ_p. An ensemble of eight other noise seeds gave pull means within ±0.7 and
standard deviations 0.7–1.5, so σ_p has the right scale; the tolerance of the noisy check is \|pull\| ≤ 4, that of the
noise-free check 0.05.

## Running it

```bash
export AZURE2_BIN=/path/to/build/src/AZURE2           # or AZURE2_BINARY; default: newest build*/src/AZURE2
python3 tests/upgrade/T1/generate.py                  # rewrite T1.azr, T1_start.azr, data/ (needs numpy)
python3 tests/upgrade/recover.py T1                   # the recovery check
python3 tests/upgrade/golden/golden.py T1             # compare with the golden outputs
```

Regenerating the data with another engine changes the data files at the level that engine differs from the one
named above; the committed files are the reference and no test regenerates them.
