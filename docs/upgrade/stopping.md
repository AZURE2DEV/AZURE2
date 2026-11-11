# Stopping and straggling: AZURE2's own module (`stopping/`)

`docs/upgrade/ARCHITECTURE.md` §7 decides that AZURE2 computes stopping and straggling itself, with no SRIM, ERYA or
LunaSRIM, every function with a stated formula, reference and uncertainty. This note describes the first branch of
that module: a leaf library that nothing in the engine links yet. The wiring into the target integration, the fit
to the IAEA stopping data base and the removal of ERYA are the next branches (§13).

## 1. What and where

| Path | What |
|---|---|
| `stopping/include/azure2/stopping/Tables.h`, `src/Tables.cpp` | the two data tables, compiled in; parsers |
| `Material.h`, `src/Material.cpp` | `Projectile`, `Material` (formula or fractions, phase, density), user tables |
| `Electronic.h`, `src/Electronic.cpp` | electronic stopping: low-energy fits, Bethe with corrections, the join, scaling |
| `Nuclear.h`, `src/Nuclear.cpp` | ZBL nuclear stopping |
| `Stopping.h`, `src/Stopping.cpp` | ε(E) of a material: Bragg additivity, user override, uncertainty, provenance text |
| `EnergyLoss.h`, `src/EnergyLoss.cpp` | N(E_in, E_out), E_out(E_in, N), range, stopping-scale derivatives |
| `Straggling.h`, `src/Straggling.cpp` | Bohr, Chu (table), Yang's charge-state term, Tschalär propagation, beam spread |
| `stopping/data/elements.tsv` | Z = 1–92: standard atomic weight, mean excitation energy I, phase |
| `stopping/data/electronic_coefficients.tsv` | per element, p and He: fit coefficients, range, uncertainty, source |
| `stopping/tests/` | six unit tests, `data/make_reference.py` → `data/reference.csv` |
| `stopping/tools/azure2_stopping.cpp` | the command-line tool `azure2-stopping` |
| `stopping/tools/erya_compare.cpp`, `compare_erya.sh` | the comparison with ERYA (§9) |

The library uses `physics/` (constants, element symbols) and the standard library only (consistency check 32). All
functions are free functions on plain records; nothing holds state but the two constant tables, parsed once and
reached by `ElementTable()` and `CoefficientTable()`.

## 2. Units

Stated once, in `Material.h`: energies in MeV, kinetic, laboratory frame, of the projectile; stopping cross
sections ε in eV/(10¹⁵ atoms/cm²) **per atom of the material** (all atoms counted); areal densities in atoms/cm²
(all atoms); straggling variances in MeV². An energy loss is ΔE [MeV] = ε N × 10⁻²¹. A stopping per atom of one
element (the "active" atom of a compound target) is ε / f, with f = `AtomFraction(material, Z)`.

## 3. The element table

`elements.tsv`, one row per element Z = 1–92: the standard atomic weight A_r and its stated interval (CIAAW abridged
values, T. Prohaska et al., Pure Appl. Chem. 94 (2022) 573; Tc, Pm, Po, At, Rn, Fr, Ra, Ac have none and carry the
mass number of a long-lived isotope), the mean excitation energy I (the ICRU Report 37 (1984) elemental values that
ICRU Report 49 (1993) and NIST PSTAR/ASTAR use; gas values for H, He, N, O, F, Ne, Cl, Ar, Kr, Xe, Rn) and the phase
of the element at room temperature. Symbols are checked against `physics/`'s `ElementSymbol`.

## 4. Electronic stopping of protons and α

### 4.1 Low energies: the fit of ICRU 49

`electronic_coefficients.tsv`, one row per element and projectile, columns `Z projectile form A1..A5 fit_max_MeV_u
rel_unc source note`. Two forms (both give eV/(10¹⁵ atoms/cm²)):

- `AZ_p`, protons (H. H. Andersen and J. F. Ziegler, *Hydrogen Stopping Powers and Ranges in All Elements*,
  Pergamon 1977, as revised and adopted in ICRU Report 49): with T the proton energy per mass unit in keV/u,
  S_low = A2 T^0.45, S_high = (A3/T) ln(1 + A4/T + A5 T), S = S_low S_high / (S_low + S_high) for T ≥ 10 keV/u;
  below 10 keV/u, S(10) (T/10)^½ (proportional to the velocity). A1, the published coefficient of A1 T^½ below
  10 keV, is kept in the table as published and not used: the continuation from S(10) reproduces PSTAR at
  1–10 keV to 0.1 % (§9), A1 T^½ is discontinuous at 10 keV.
- `Z_He`, helium-4 (J. F. Ziegler, *Helium Stopping Powers and Ranges in All Elements*, Pergamon 1977, as adopted
  in ICRU Report 49): with T the He-4 energy in MeV, S_low = A1 (1000 T)^A2, S_high = (A3/T) ln(1 + A4/T + A5 T),
  S = S_low S_high / (S_low + S_high) for T ≥ 1 keV; velocity-proportional below.

The variable of each form is taken at the projectile's velocity: T/u = E/M (M the projectile's mass in u), so d and t
use the proton row and ³He the helium row at the same velocity.

### 4.2 High energies: Bethe with corrections

ε = K z₁² Z₂/β² · L, K = 4π r_e² m_e c² (5.099 × 10⁻⁴ eV/(10¹⁵ atoms/cm²) from `physics/` constants),
L = L₀ − C/Z₂ + z₁ L₁ + z₁² L₂:

- L₀ = ln(2 m_e c² β²γ² / I) − β² (W_max ≈ 2 m_e c² β²γ² for an ion; no density effect: below 10 MeV/u it is
  negligible).
- Shell correction C (H. Bichsel's notation; formula of W. H. Barkas and M. J. Berger, NASA SP-3013 (1964), as given
  in ICRU Report 37): C = (0.422377 η⁻² + 0.0304043 η⁻⁴ − 0.00038106 η⁻⁶) 10⁻⁶ I² + (3.858019 η⁻² − 0.1667989 η⁻⁴
  + 0.00157955 η⁻⁶) 10⁻⁹ I³, η = βγ, I in eV. Valid for η ≥ 0.13 (7.9 MeV/u); below, C is held at its η = 0.13
  value. This is the largest approximation of the Bethe branch: between 3 and 8 MeV/u in heavy elements it is what
  leaves ≤ 3 % against ICRU 49 (§9).
- Barkas term L₁ = F(b/x^½) / (Z₂^½ x^{3/2}), x = β²/(α² Z₂) (J. C. Ashley, R. H. Ritchie and W. Brandt, Phys. Rev.
  B 5 (1972) 2393); F is the 47-point table and b the per-element screening parameters of Geant4's
  `G4EmCorrections` (after ICRU 49), linear interpolation in F, F ∝ 1/W above the table. Geant4 multiplies the term
  by an empirical 1.29; this module does not.
- Bloch term L₂ = −y² Σₙ 1/(n(n² + y²)) = ψ(1) − Re ψ(1 + iy), y = z₁α/β (F. Bloch, Ann. Phys. 16 (1933) 285).

### 4.3 The join

The fit alone up to T₁ = `fit_max_MeV_u` of the row (1 MeV/u for every row), Bethe alone from T₂ = 3 MeV/u
(`kBetheFrom_MeV_u`), and between them linear interpolation in ln S – ln T from S_fit(T₁) to S_Bethe(T₂). The choice
was made by comparing joins at 0.5/1 MeV/u and 3/4/6 MeV/u with PSTAR and ASTAR for H, He, N, O, Ne, Al, Ti, Au:
T₁ = 1, T₂ = 3 gives the smallest largest deviation (≤ 3 % except the rows of §9.1). ICRU 49 itself joins its fits
to Bethe by a smooth interpolation over a similar interval.

### 4.4 Heavier ions

z₁ ≥ 3: the proton value at the same velocity times (q z₁)², q/z₁ = 1 − exp(−0.95 v_r), v_r = v/(v₀ z₁^{2/3})
(J. M. Pierce and M. Blann, Phys. Rev. 173 (1968) 390). This is a scaling, not an evaluation: its stated
uncertainty is 20 % for v_r ≥ 1 and 50 % below, where the equilibrium charge underestimates the effective charge
of stopping. The Brandt–Kitagawa effective charge (as in ZBL 1985) is the better choice for the next branch.

## 5. Nuclear stopping

ZBL universal (J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping and Range of Ions in Solids*, Pergamon 1985):
ε_r = 32.53 M₂ E / (Z₁Z₂(M₁+M₂)(Z₁^0.23+Z₂^0.23)) (E in keV), s_n = ln(1 + 1.1383 ε_r) / (2(ε_r + 0.01321 ε_r^0.21226
+ 0.19593 ε_r^½)) for ε_r ≤ 30, ln ε_r/(2ε_r) above, S_n = 8.462 Z₁Z₂M₁ s_n / ((M₁+M₂)(Z₁^0.23+Z₂^0.23)). M₂ is the
standard atomic weight. The test computes the nuclear stopping of the universal potential itself (classical
scattering integral, `make_reference.py`): the fit agrees to 0.4 % for ε_r = 0.003–0.3, 1.04 % at 3, and is 4.3 %
low at 30, where it switches to its high-energy form. Assigned uncertainty 10 % (`kNuclearRelUncertainty`): ICRU 49's
nuclear stopping differs from ZBL by up to that amount at keV energies (§9.1).

## 6. Compounds and user stopping

Bragg additivity: ε = Σᵢ fᵢ (S_e,i + S_n,i), fᵢ the atom fractions, times `CompoundCorrectionFactor`, which is 1 in
this branch. That function is the place of the core-and-bond correction (J. F. Ziegler and J. M. Manoyan, Nucl.
Instrum. Methods B 35 (1988) 215) or of a correction fitted to measured compound data; the next branch fills it.

The phase of the material is recorded; the table rows are for one phase per element (gas for H, N, O, ...), and the
provenance text says when they differ (no gas–solid correction in this branch).

A user table or formula overrides the evaluation for one projectile (matched by charge and by mass within 0.01 u):
`UserStopping` with `energy_MeV`/`stopping` (interpolated linearly in ln E – ln ε, refused outside its range) or a
`formula` (any callable ε(E [MeV]); the engine will pass its expression parser), a stated relative uncertainty and a
provenance string. `ReadStoppingTable` reads two columns with a unit line, `E[MeV]` or `E[keV]`, and
`eps[eV/1e15atoms/cm2]`, `S[MeV*cm2/g]` (converted with the material's mean atomic weight) or `S[keV/um]` (needs the
mass density). `StoppingProvenance(material, projectile)` prints where every number comes from, down to the table
line.

## 7. Uncertainty

Every evaluation returns its relative 1σ uncertainty (`StoppingValue::rel_uncertainty`), the prior width of the
stopping-scale nuisance in Phase 5. In this branch the values are **assigned**, not fitted: 4 % for a proton row and
5 % for a helium row below 1 MeV/u (the scatter of measured data about ICRU 49 is of this order; the stated
uncertainties of ICRU 49 could not be consulted for this branch), 2 % in the Bethe region, the larger of the two in
the join, 20 % / 50 % for effective-charge scaling, 10 % for the nuclear part. A compound adds its constituents'
electronic uncertainties linearly (fully correlated, since they come from one evaluation) and the nuclear part in
quadrature. The fit to the IAEA data base (§13) replaces these numbers by a covariance.

## 8. What a target needs

On any ε(E) (`StoppingFunction`) and on a material:

- N(E_in, E_out) = ∫_{E_out}^{E_in} dE/ε, computed in u = √E (dE = 2u du, which makes the integrand finite at E = 0
  for a velocity-proportional ε), adaptive Gauss–Kronrod 7–15 to a relative tolerance of 10⁻¹⁰;
- E_out(E_in, N), the root of N(E_in, E) = N by safeguarded Newton steps (dN/dE_out = −1/ε); refused if N ≥ R;
- R(E) = N(E, 0);
- with ε → s ε, at s = 1: ∂N/∂s = −N, ∂R/∂s = −R, ∂E_out/∂s = −N ε(E_out).

The tests check them against closed forms for ε = c E^p (p = 0.5, 0, 0.3, −0.8) to 10⁻⁹, the derivatives against
central differences, and E_out and N against each other on Ta₂O₅ with all the joins inside the interval. The
quadrature refuses an integrand it cannot converge on (a discontinuous user formula, for instance) instead of
returning an inaccurate number.

## 9. Comparisons (not gates)

### 9.1 With ICRU 49 (NIST PSTAR / ASTAR)

Downloaded for comparison only, never committed (§12). Largest |ours/ICRU 49 − 1| in each energy band, electronic
stopping (E in MeV, lab; "C (gr.)" is graphite, "C (am.)" amorphous carbon):

| | 1–10 keV | 10–100 keV | 0.1–1 MeV | 1–3 MeV | 3–10 MeV |
|---|---|---|---|---|---|
| p H | 0.1 % | 0.1 % | +2.1 % (0.45) | +1.1 % | 0.2 % |
| p He | 0.0 % | 0.0 % | −0.9 % | −1.2 % | −0.5 % |
| p C (am.) | **−17.1 %** (1 keV) | **−9.0 %** (10 keV) | +1.9 % | +2.2 % | +1.3 % |
| p C (gr.) | **−20.7 %** | **−14.3 %** | **−8.1 %** (0.1) | +0.7 % | −0.9 % |
| p N | 0.0 % | 0.0 % | −0.8 % | +1.3 % | +1.7 % |
| p O | 0.1 % | 0.0 % | −1.4 % | +1.3 % | +1.8 % |
| p Ne | 0.1 % | 0.1 % | +1.9 % | −1.4 % | +1.4 % |
| p Al | 0.1 % | 0.0 % | +2.3 % | +2.4 % | +1.7 % |
| p Ti | 0.0 % | 0.0 % | −0.7 % | +1.7 % | +2.3 % |
| p Au | 0.1 % | 0.1 % | +2.3 % | +1.8 % | +2.9 % (4.5) |
| α H | 0.0 % | 0.0 % | 0.0 % | −1.0 % | −2.2 % |
| α He | 0.0 % | 0.1 % | 0.0 % | −0.6 % | −2.4 % |
| α C (gr.) | 0.0 % | 0.1 % | 0.0 % | +1.1 % | −2.3 % |
| α C (am.) | +4.0 % | +7.5 % | **+8.9 %** | +5.3 % | +2.4 % |
| α N | 0.0 % | 0.0 % | 0.0 % | +1.6 % | −1.2 % |
| α O | 0.0 % | 0.0 % | 0.0 % | −2.4 % | −4.5 % (4.5) |
| α Ne | 0.0 % | 0.0 % | 0.0 % | +2.3 % | +1.9 % |
| α Al | 0.1 % | 0.0 % | 0.0 % | −0.9 % | −2.4 % |
| α Ti | 0.0 % | 0.0 % | 0.1 % | −0.6 % | −3.0 % |
| α Au | **−23.1 %** (1 keV) | **−15.8 %** | **−7.2 %** (0.1) | −2.8 % | −2.1 % |

Where the same coefficients are used the agreement is 0.1 %; the 1–3 % above 0.3 MeV/u is the join (§4.3). Three
rows do not reproduce ICRU 49 within the assigned uncertainty, and the reason is the coefficients, not the form:
- **p in carbon below 40 keV** (−17 % amorphous at 1 keV, −9 % at 10 keV; graphite −8 % up to 0.1 MeV): ICRU 49 has
  phase-specific carbon curves at low energy that are not in the transcribed row.
- **α in carbon**: the row is ICRU 49's graphite (0.0 % against ASTAR graphite); amorphous carbon differs by up to
  9 % (a phase choice, §6).
- **α in gold below 0.5 MeV** (−23 % at 1 keV, −7 % at 0.1 MeV): ASTAR's gold is neither the transcribed Ziegler
  1977 row nor the alternative row the source names as ICRU's (which is worse, −11 % to −14 %).
Total stopping (electronic + nuclear) differs further at 1–10 keV by the nuclear part: ZBL against ICRU 49's nuclear
stopping, up to −6 % of the total for p in H at 1 keV.

### 9.2 With ERYA

`stopping/tools/compare_erya.sh <build> <out.tsv>` (a build with `USE_ERYA`) evaluates ERYA's `SRIMUtilities` (its
SRIM-2013 proton coefficients, read only, through the code the GUI used) and this module's electronic stopping for p
and α on H, He, C, N, O, Ne, Al, Ti, Ta, Au, Ta₂O₅, TiN from 20 keV to 10 MeV. ERYA has no helium data: its α column
is 4 × ERYA(p) at the same velocity, the bare-charge scaling its equation would give. Relative difference
ours/ERYA − 1 of the electronic stopping:

| E (MeV) | 0.02 | 0.1 | 0.5 | 1 | 2 | 5 | 10 |
|---|---|---|---|---|---|---|---|
| p H | +0.632 | +0.140 | +0.024 | −0.027 | −0.045 | +0.008 | +0.086 |
| p He | +1.399 | +0.289 | −0.117 | −0.121 | −0.090 | −0.014 | +0.050 |
| p C | −0.101 | −0.006 | −0.012 | +0.012 | −0.018 | −0.028 | −0.028 |
| p N | +0.433 | +0.183 | +0.001 | −0.009 | −0.003 | +0.013 | +0.011 |
| p O | +0.432 | +0.242 | +0.009 | −0.009 | −0.014 | −0.003 | −0.005 |
| p Ne | +0.077 | −0.019 | −0.004 | −0.014 | −0.011 | +0.002 | −0.006 |
| p Al | +0.019 | −0.015 | +0.015 | +0.023 | +0.006 | +0.009 | +0.002 |
| p Ti | −0.023 | +0.014 | −0.032 | −0.013 | +0.001 | +0.024 | +0.017 |
| p Ta | −0.007 | −0.002 | +0.015 | −0.008 | +0.006 | +0.039 | +0.011 |
| p Au | −0.044 | −0.024 | −0.007 | −0.018 | −0.010 | +0.024 | +0.012 |
| p Ta₂O₅ | +0.208 | +0.116 | +0.012 | −0.009 | −0.002 | +0.025 | +0.006 |
| p TiN | +0.110 | +0.069 | −0.021 | −0.012 | −0.000 | +0.021 | +0.016 |
| α H | −0.605 | −0.520 | −0.290 | −0.059 | −0.011 | −0.058 | −0.060 |
| α He | +0.011 | −0.159 | −0.246 | −0.166 | −0.129 | −0.106 | −0.087 |
| α C | −0.523 | −0.414 | −0.267 | −0.132 | −0.006 | −0.010 | −0.042 |
| α N | −0.515 | −0.365 | −0.196 | −0.029 | +0.016 | −0.010 | −0.014 |
| α O | −0.314 | −0.284 | −0.174 | −0.031 | +0.012 | −0.052 | −0.035 |
| α Ne | −0.325 | −0.385 | −0.308 | −0.142 | −0.006 | −0.003 | −0.021 |
| α Al | −0.563 | −0.453 | −0.246 | −0.108 | −0.025 | −0.004 | −0.014 |
| α Ti | −0.538 | −0.477 | −0.280 | −0.128 | −0.031 | −0.002 | −0.005 |
| α Ta | −0.416 | −0.339 | −0.251 | −0.091 | −0.027 | −0.028 | −0.003 |
| α Au | −0.627 | −0.472 | −0.227 | −0.107 | −0.016 | +0.012 | −0.012 |
| α Ta₂O₅ | −0.369 | −0.312 | −0.214 | −0.062 | −0.010 | −0.038 | −0.015 |
| α TiN | −0.531 | −0.444 | −0.253 | −0.095 | −0.016 | −0.004 | −0.008 |

Why, where it is large:
- **Protons in solids** (C, Al, Ti, Ta, Au): within ±4 %, the size of two evaluations of the same data.
- **Protons in the gases H, He, N, O below 0.2 MeV**: ERYA is 30–60 % lower. This module reproduces PSTAR there to
  0.1 % (§9.1), so the difference is ERYA's: its single eight-parameter formula, evaluated at all energies (its
  comment: "use high-energy formula for all energies"), or SRIM's coefficients for these elements (not inspected:
  the file is only run, never read for reuse); which of the two, this comparison cannot tell. Compounds with O and N (Ta₂O₅, TiN) inherit +7 to +21 %.
- **α everywhere below 2 MeV**: ERYA has no helium stopping; 4 × ERYA(p) assumes a bare He²⁺ at all velocities,
  while helium's effective charge falls well below 2 at these velocities. 30–60 % at 20 keV, 10–30 % at 0.5–1 MeV,
  ≤ 5 % above 2 MeV. In the engine an α target with ERYA's equation was therefore wrong by these amounts unless the
  user scaled it.

## 10. Straggling

### 10.1 Bohr

dΩ_B²/dN = 4π z₁² e⁴ Z₂ (N. Bohr, Phil. Mag. 30 (1915) 581), e² = α ħc, MeV² cm²; a compound uses the atom-fraction
weighted Z₂ (Bragg additivity of the variance). 0.2606 keV² per 10¹⁸ atoms/cm² for z₁ = Z₂ = 1 (the rounded "0.26"
of the literature, checked to 0.3 %).

### 10.2 Chu's correction

Ω² = H(E/M₁, Z₂) Ω_B² (W. K. Chu, Phys. Rev. A 13 (1976) 2057; Hartree–Fock–Slater). H is data: a `ChuTable`
(E/M₁ in keV/u × Z₂, linear in E/M₁, 1 above the table, its first value below). **No table is shipped**: Chu's
values are published as a table and figures (reproduced in the IBA handbooks and distributed with SIMNRA as
`CHU_CORR.DAT`, whose terms were not checked), and Yang et al.'s fit of H could not be obtained as numbers from an open
source for this branch. The function, the model switch and the tests (with a test table) are in place.

### 10.3 The charge-state term of Yang, O'Connor and Wang

Ω²/Ω_B² = γ² H + (ΔΩ/Ω_B)² (Q. Yang, D. J. O'Connor, Z. Wang, Nucl. Instrum. Methods B 61 (1991) 149), with, for
hydrogen ions (γ = 1), (ΔΩ/Ω_B)² = B₁Γ/((E − B₂)² + Γ²), Γ = B₃(1 − exp(−B₄E)); for helium and heavier,
Z₁^{4/3}/Z₂^{1/3} · C₁Γ/((ε − C₂)² + Γ²), Γ = C₃(1 − exp(−C₄ε)), ε = E/(Z₁^{3/2}Z₂^{1/2}); E in MeV/u. γ² is 1 for
hydrogen, the ZBL (1985) helium effective charge for helium, the Pierce–Blann fraction squared above. **The
constants B₁–B₄, C₁–C₄ are not shipped**: the paper is not open access and no open source with the numbers was found
(the SIMNRA guide gives the form, not the constants). `YangCoefficients` takes them; the tests check the form with
given constants.

### 10.4 Propagation and beam spread

σ²(x) = S(E(x))² [σ₀²/S(E₀)² + ∫₀ˣ (dΩ²/dx′)/S(E(x′))² dx′] (B. Tschalär, Nucl. Instrum. Methods 61 (1968) 141 and
64 (1968) 237), computed as S(E_out)² ∫_{E_out}^{E_in} (dΩ²/dN)/S³ dE in √E. The beam's entrance spread σ₀ is carried
through the layer by S(E_out)/S(E_in) and added in quadrature to the straggling (`StragglingWidths`: straggling,
beam, total). Tested: constant S gives Ω² + σ₀² exactly, a thin layer gives Bohr, and ε = c E^p against the closed
form to 10⁻⁸.

## 11. Command line and Python

`azure2-stopping` (built in `<build>/stopping/`):

```
azure2-stopping stopping   --projectile a --material Ta2O5 --from 0.02 --to 10 --points 50
azure2-stopping range      --projectile p --material TiN --from 0.1 --to 5
azure2-stopping layer      --projectile p --material Ta2O5 --energy 1.0 --areal-density 1e18 --beam-sigma 0.0005
azure2-stopping provenance --projectile a --material Ta2O5
```

Python: not bound in this branch. Binding the functions as a `pyazr.stopping` submodule needs one new file
`api/src/pybind_stopping.cpp` (array inputs, numpy outputs, as `pybind_physics.cpp`), its line in
`api/src/CMakeLists.txt` with `azure2::stopping` linked, and two lines in `api/src/pybind_api.cpp` (the declaration
`void BindStopping(py::module_ &m);` and the call `BindStopping(m);`), a file other branches are editing; it is left
for the wiring branch.

## 12. Where the numbers come from, and their terms

- The fit coefficients are those of ICRU Report 49 (1993) for the Andersen–Ziegler (1977) and Ziegler (1977) forms,
  **transcribed from Geant4** (`G4BraggModel.cc`, `G4BraggIonModel.cc`; Geant4 Software License 1.0, which permits
  redistribution with its notice), and checked against PSTAR/ASTAR to 0.1 % (§9.1). Neither SRIM's nor ERYA's
  coefficient file was used. The Barkas F table and b parameters are from Geant4's `G4EmCorrections.cc` (same
  licence). Whether the Geant4 notice must accompany these numbers in AZURE2 is a decision for the maintainers.
- I values: ICRU Report 37 / 49 (checked against Geant4's `G4NistMaterialBuilder.cc`, which lists the same values).
- Atomic weights: CIAAW, 2021 abridged table.
- NIST PSTAR/ASTAR tables and the SIMNRA user's guide were downloaded for comparison and reading only, into the
  scratch directory; nothing of them is committed.

## 13. Next branches

1. The fit to the IAEA stopping-power data base: per element, gas and solid separately where data exist, the same
   functional form with coefficients and their covariance fitted to measured data (after checking the data base's
   redistribution terms), replacing the transcribed rows and the assigned uncertainties; carbon and gold first.
2. Compounds: the core-and-bond correction (`CompoundCorrectionFactor`) where measured compound data show one.
3. Straggling data: a Chu table and Yang's constants with their terms checked.
4. Heavy ions: the Brandt–Kitagawa effective charge.
5. Wiring into the engine (Phase 5): `TargetEffect`'s stopping expression, its range-integral energy loss and Bohr
   straggling replaced by this module (a material and projectile per target in format 3), the binding, and the
   removal of ERYA and `SRIM2013.xml`.
