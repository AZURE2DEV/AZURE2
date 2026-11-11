# T3 — identical bosons, neutrons and special channels, ⁸Be

The project files are format 3 (`docs/upgrade/format3_tree.md`); the old lines quoted on this page are those of the
originals in `tests/convert/fixtures/tests/upgrade/T3/`.

T3 is a synthetic R-matrix evaluation of ⁸Be with four partitions: α + α (identical spin-0 bosons), p + ⁷Li,
n + ⁷Be and the β⁻ decay of ⁸Li. It is generated from the parameters below; a fit with the present engine must
return them (`../recover.py`), and the present engine's outputs on it are kept in `../golden/T3/`. The level scheme
borrows the ⁸Be states near the nucleon thresholds, but the widths and background poles are chosen for coverage, not
taken from an evaluation. The project runs in 0.1 s (calculation) and about 1 s (fit) on two threads.

## The system

| pair | particles | J^π (light, heavy) | S (MeV) | a_c (fm) | notes |
|---|---|---|---|---|---|
| 1 | α + α | 0⁺, 0⁺ | S_α = −0.09184 | 5.5 | identical bosons; hybrid Woods–Saxon potential |
| 2 | p + ⁷Li | 1/2⁺, 3/2⁻ | S_p = 17.25444 | 4.0 | |
| 3 | n + ⁷Be | 1/2⁺, 3/2⁻ | S_n = 18.89955 | 4.0 | neutral: total cross section, scattering lengths |
| 4 | e⁻ + ⁸Li (β⁻ decay) | 1/2⁺, 2⁺ | Q_β = 16.00445 | 5.5 | pair type 20; Gamow–Teller channel (l = 1) |

Atomic masses are AME2020 (`../lib/upgrade_project.py`, `MASS`), S = m₁ + m₂ − m(⁸Be) and Q_β = m(⁸Li) − m(⁸Be) with
931.49410242 MeV/u. For a β pair the engine reads Q_β in the S column and computes the endpoint energy
Q_β − E_x; the electron radius a = 5.5 fm enters only the integrated Fermi function.

The α + α channel functions are those of the hybrid model (`<potential>`, pair 1): Coulomb plus a Woods–Saxon well
V₀ = 10 MeV, R = 3.0 fm, a = 0.6 fm, integrated by Numerov to the channel radius. The other pairs use the plain Coulomb
(or, for n + ⁷Be, spherical Bessel) functions. Changing V₀ to 30 MeV moves χ² of the α + α data from 86.7 to 90.3, so
the potential is active.

## The true parameters

Brune basis: E_x of ⁸Be (MeV), observed partial widths Γ (eV, c.m.), β amplitudes B (dimensionless, the engine's
Gamow–Teller matrix element in its own normalization). Free parameters are marked *.

| J^π | E_x (MeV) | channel | value |
|---|---|---|---|
| 2⁺ | 3.03 * | α + α, l = 2 | Γ_α = 1.50 MeV * |
| | | β (⁸Li 2⁺ → 2⁺, GT) | B = 1.0 * |
| 4⁺ | 11.35 | α + α, l = 4 | Γ_α = 3.5 MeV |
| 1⁺ | 17.64 * | p + ⁷Li, l = 1, s = 1 | Γ_p = 11.0 keV * |
| | | p + ⁷Li, l = 1, s = 2 | Γ_p = 1.0 keV |
| 2⁻ | 19.00 * | p + ⁷Li, l = 0, s = 2 | Γ_p = 120 keV * |
| | | n + ⁷Be, l = 0, s = 2 | Γ_n = 60 keV * |
| 0⁺ | 25.0 | α + α, l = 0 | Γ_α = 10 MeV (background) |
| 2⁺ | 22.0 | α + α, l = 2; p, l = 1, s = 2; n, l = 1, s = 2; β | 4.0 MeV; 1.0 MeV; 1.0 MeV; B = −0.5 (background) |
| 1⁻ | 22.0 | p, l = 0, s = 1; n, l = 0, s = 1 | 2.0 MeV each (background) |

Unnatural-parity levels (1⁺, 2⁻, 1⁻) have no α + α channel: two identical spin-0 bosons couple only to J^π = even⁺.
The 2⁺ background pole carries the ⁷Li(p,α) and ⁷Be(n,α) strength. `T3_start.azr` moves every free energy by
+10 keV and every free width and amplitude by a factor 1.15.

## The data

Engine truth plus Gaussian noise from `numpy.random.default_rng(8)` (standard deviation as stated, also written as
the uncertainty column). Energies are E_lab of the projectile (α on α, p on ⁷Li, n on ⁷Be), except for the
β-delayed spectrum, whose energy column is the laboratory energy E_α of one α particle from the decay at rest
(E_c.m.(α + α) = 2 E_α); the engine reads it so (`docs/source/user_guide/particle_pairs.rst`, "Relativistic
kinematics").

| file | observable (isDiff) | frame | E grid (MeV) | angle | noise | points |
|---|---|---|---|---|---|---|
| `aa_cm050/090` | dσ/dΩ α + α (4) | c.m. | E_lab 1.0–8.0, step 0.1 | θ_c.m. = 50°, 90° | 3 % | 2 × 71 |
| `aa_phase_d` | δ(J = 2, l = 2) α + α (2) | — | E_lab 1.0–8.0, step 0.25, \|δ\| < 80° | — | 1° absolute | 22 |
| `beta_alpha` | β-delayed α spectrum (0) | — | E_α 0.40–6.00, step 0.08 | — | 3 % | 71 |
| `p7li_pp_lab150` | dσ/dΩ ⁷Li(p,p) (1) | lab | E_lab 0.30–2.40, step 0.025 | θ_lab = 150° | 3 % | 85 |
| `p7li_pn` | σ ⁷Li(p,n)⁷Be (0) | — | E_lab 1.90–2.40, step 0.01 | — | 5 % | 51 |
| `p7li_pa` | σ ⁷Li(p,α)α (0) | — | E_lab 0.30–2.40, step 0.05 | — | 5 % | 43 |
| `n7be_total` | σ_tot n + ⁷Be (0, `total 1`) | — | E_lab 0.01–0.40, step 0.01 | — | 3 % | 40 |
| `n7be_np` | σ ⁷Be(n,p)⁷Li (0) | — | E_lab 0.01–0.40, step 0.01 | — | 5 % | 40 |

`total 1` makes the n + ⁷Be elastic segment the total cross section by the optical theorem,
σ_tot = (2π/k²) Σ_J g_J Σ_c∈n+⁷Be (1 − Re U_cc), with g_J = (2J + 1)/((2j₁ + 1)(2j₂ + 1)).

The `<scatteringLengths>` block constrains the complex s-wave scattering lengths of n + ⁷Be in the channel spins
s = 1 and s = 2 (J^π = 1⁻, 2⁻). The truth is a₀(s = 1) = 3.137375 − 0.131465 i fm and a₀(s = 2) = 0.566232 − 1.935850 i fm
(imaginary because n + ⁷Be lies above the open p + ⁷Li and α + α channels); the constraints are these plus Gaussian
noise of 0.1 fm (seed 1008), with σ = 0.1 fm on both parts.

χ² of the true parameters: 509.7 for N = 494 data points, plus 9.0 from the four scattering-length terms.

The `<segmentsTest>` lines compute, for E_lab(n) = 0.01–0.39 MeV (step 0.01 MeV), the n + ⁷Be total cross section and its three
parts (elastic, (n,p), (n,α)) for the optical-theorem identity test, the α + α phase shifts for l = 0, 2, 4, an α + α
angular distribution at E_lab = 4 MeV (10°–170°, symmetric about 90°), and the β-delayed α spectrum for E_α = 0.1–7.5
MeV. The rate of ⁷Li(p,α)α is computed at T₉ = 0.01–3.

## What T3 covers, against `docs/upgrade/ARCHITECTURE.md` §10

| item | status |
|---|---|
| α + α, identical bosons | c.m. cross sections, d-wave phase shift, symmetric test angular distribution |
| p + ⁷Li | elastic at θ_lab = 150°, (p,n), (p,α) to identical bosons |
| n + ⁷Be | total cross section, (n,p); (n,α) in the test segments |
| β-delayed α | data and test segment |
| optical theorem | `total 1` data segment; identity test `tests/pyazr/upgrade_optical_theorem_test.py` |
| scattering length | complex constraints in the fit; `scatteringLengths.out` in the golden outputs |
| hybrid potential channel | α + α channel functions from the Woods–Saxon well |

Not covered: a neutron capture channel, and the β-delayed proton or neutron branches (no ⁸Li branch exists
energetically). These need nothing beyond what T1 and T3 already exercise.

## Recovery (present engine, 2026-11-06)

`python3 ../recover.py T3` with the build of `feature/upgrade` c1cc92c7, two threads: valid minimum after 129 calls
in 1.0 s; χ² = 498.7 at the fit against 509.7 at the truth. Largest \|pull\| 1.82 (Γ_p of the 1⁺ level, s = 1: the
reduced width 0.6674 ± 0.0203 MeV^½ against 0.6306); all others ≤ 1.53. On noise-free data (data values and the
scattering lengths at the truth) the fit returns every parameter within 0.036 σ_p.

## Running it

```bash
export AZURE2_BIN=/path/to/build/src/AZURE2           # or AZURE2_BINARY; default: newest build*/src/AZURE2
python3 tests/upgrade/T3/generate.py                  # rewrite T3.azr, T3_start.azr, data/, temperatures.dat
python3 tests/upgrade/recover.py T3
python3 tests/upgrade/golden/golden.py T3
```
