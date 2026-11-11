# Stopping and straggling: AZURE2's module (`stopping/`), SRIM-2013 through LunaSRIM

`docs/upgrade/ARCHITECTURE.md` §7: AZURE2 computes stopping and straggling itself, with **SRIM-2013's stopping
powers**, ported to C++ from **LunaSRIM** (J. Skowronski, github.com/skowrons94/LunaSRIM, MIT), the pure-Python
re-implementation of SRIM-2013 that reproduces SRIM's tables. The PI decided this on 2026-11-11: "I would use the SRIM
stopping powers from my LunaSRIM repository. They are free to use as long as SRIM is cited." An AZURE2 fit and a LUNA
analysis made with LunaSRIM then use the same stopping power. The module is a leaf library that nothing in the engine
links yet; the wiring and the removal of ERYA are the next branch (§13).

**Cite SRIM** in every result that uses these stopping powers: J. F. Ziegler, M. D. Ziegler and J. P. Biersack,
*SRIM – The stopping and range of ions in matter (2010)*, Nucl. Instrum. Methods Phys. Res. B 268 (2010) 1818,
doi:10.1016/j.nimb.2010.02.091.

## 1. What and where

| Path | What |
|---|---|
| `stopping/include/azure2/stopping/Tables.h`, `src/Tables.cpp` | SRIM's data files: reader, compiled-in copy, directory override |
| `Material.h`, `src/Material.cpp` | `Projectile`, `Material` (formula or fractions, phase, compound correction, density), user tables |
| `Electronic.h`, `src/Electronic.cpp` | SRIM's electronic stopping, phase factor, fade-out (§4) |
| `Nuclear.h`, `src/Nuclear.cpp` | SRIM's nuclear stopping (§5) |
| `Stopping.h`, `src/Stopping.cpp` | ε(E) of a material: Bragg sum, user override, uncertainty, break points, provenance text, the citation |
| `EnergyLoss.h`, `src/EnergyLoss.cpp` | N(E_in, E_out), E_out(E_in, N), range, stopping-scale derivatives |
| `Straggling.h`, `src/Straggling.cpp` | TRIM's and Bohr's straggling, Tschalär propagation, beam spread |
| `stopping/data/srim/` | `SCOEF03.DAT`, `SNUC03.dat` (SRIM-2013, unchanged) and their `README.md` (source, SHA-256, terms, citation) |
| `stopping/tests/` | six unit tests (seven ctest entries), `data/srim_reference.csv` (LunaSRIM's values) |
| `stopping/tools/azure2_stopping.cpp` | the command-line tool `azure2-stopping` (§11) |
| `stopping/tools/make_srim_reference.py` | writes `srim_reference.csv` with LunaSRIM (any Python that has `lunasrim`) |

The library uses `physics/` (constants, element symbols) and the standard library only (consistency check 32). All
functions are free functions on plain records. The only state is the constant table of SRIM's data, parsed once and
reached by `SrimTables()`; a caller may also read the files from any directory into a `SrimData` it owns
(`ReadSrimData`) and pass it to the element-level functions.

## 2. Units

Stated once, in `Material.h`: energies in MeV, kinetic, laboratory frame, of the projectile; stopping cross
sections ε in eV/(10¹⁵ atoms/cm²) **per atom of the material** (all atoms counted); areal densities in atoms/cm²
(all atoms); straggling variances in MeV². An energy loss is ΔE [MeV] = ε N × 10⁻²¹. A stopping per atom of one
element (the "active" atom of a compound target) is ε / f, with f = `AtomFraction(material, Z)`.

Inside the SRIM formulas the variable is SRIM's: the energy per mass unit e = E [keV] / M₁ [u], in keV/u. The
projectile's mass enters only there and in the nuclear kinematics; it is **the caller's** (`Projectile::mass_u`;
`Proton()`, `Alpha()`, ... carry AME2020 atomic masses). LunaSRIM's default, when the caller gives none, is the mass
of the most abundant isotope (`SCOEF03.DAT` column 3); a caller that wants to reproduce such a LunaSRIM run passes
that mass.

## 3. SRIM's data files

`stopping/data/srim/SCOEF03.DAT` (one row per element Z = 1–92: masses, densities, the stopping fits of H, He and Li
ions, the heavy-ion scaling polynomials, the mean excitation energies of the solid and the gas, two shell-correction
polynomials, the high-energy effective-charge parameters) and `SNUC03.dat` (four coefficients of the nuclear stopping
fit for every pair Z₁ ≤ Z₂). Column layout: LunaSRIM's `data_files.py`; here `Tables.h` (`SrimElement` names the
column of every field) and `ParseSrimData`. The files are copied unchanged from LunaSRIM 1.0.0 (SHA-256 in their
README) and compiled into the library by CMake (raw string pieces of 30 000 characters, `src/TableText.cpp.in`), so a
built library needs no file at run time. `AZURE2_SRIM_DATA_DIR`, if set, names a directory to read them from instead;
without the files in the source tree, the build embeds nothing and the variable is required. The test
`stopping_srim_directory` runs the tables test with the variable set.

The element masses used are SRIM's isotope-averaged atomic weights (column 4): for the nuclear stopping's M₂ and for
`MeanAtomicWeight` (unit conversions of user tables). The first branch's `elements.tsv` (CIAAW weights, ICRU 37 I
values, phases) is gone: SRIM's file gives the weights and I values it uses, and the natural phase is SRIM's list of
gaseous elements (`IsNaturalGas`); nothing else of that table was used.

## 4. Electronic stopping

Ported from LunaSRIM's `stopping.py`; each C++ function names the LunaSRIM function it comes from. All in keV/u and
eV/(10¹⁵ atoms/cm²):

- **SRIM's fit** (`SrimFit`, from `_fit`), coefficients c₁…c₁₀ of the element and ion class:
  S_low = c₁ e^c₂ + c₃ e^c₄ + c₉ e^½, S_high = c₅/e^c₆ · (ln(c₇/e + c₈ e) + c₁₀/e^0.2), S = S_low S_high/(S_low + S_high).
- **Bethe** (`SrimBethe`, from `bethe`): β² = 1 − (1 + e/931494)⁻², L = ln(1022000 β²/(1 − β²)) − β² − D − C − ln I
  + Barkas + Bloch, S = Z₁² Z₂ L / (1961.16 β²); C and D the two shell polynomials in ln(e/1000) (D = 0 for a gas),
  I the solid's or the gas's (gas: one value for H ions, one for heavier), Barkas = Z₁ (e/1000) q / (e/1000 + q) with
  q = 1.5 e^−0.4 + 45000/(Z₂ e^1.6) (solid) or 1.35 e^−0.4 (gas), Bloch = −y (1.20206 − y (1.042 − 0.8549 y + 0.343 y²)),
  y = (Z₁/137.036)²/β². The numbers are SRIM's own.
- **H, He, Li** (`SrimLightIon`, from `_light_ion`): the fit between 2 keV/u and E_max (1000 keV/u for H, 2000 for He
  and Li), linear in e from the fit at E_max to Bethe at 1.2 E_max, Bethe above, S(2 keV/u) (e/2)^½ below 2 keV/u.
  d and t are H ions, ³He a He ion: the same stopping at the same e.
- **Z₁ ≥ 4** (`SrimElementStopping`, from `electronic_stopping_element`): helium's stopping at the same e (clamped to
  2…30 000 keV/u), times the ion's polynomial and the target's polynomial in ln e and Z₁²/4; velocity-proportional
  below 2 keV/u; above 30 MeV/u scaled by the ratio of SRIM's high-energy effective charge (`SrimHighEnergyCharge`,
  from `_hi_charge`) times helium's Bethe value at e and at 30 MeV/u.
- **Phase** (`SrimPhaseFactor`, `SrimFadeOut`, from `phase_factor`, `fade_correction`): a natural gas (H, He, N, O,
  F, Ne, Cl, Ar, Br, Kr, I, Xe, At, Rn) in a solid target gets 0.92, a natural solid in a gas target 1/0.92, noble gases
  and elements in their natural phase 1; the factor F (times the compound correction C, §6) is applied fully at low
  velocity and faded out: 1 + (F C − 1)/(1 + exp(1.48 (√(e/25) − c))), c = 7 for H ions, 4 for heavier.
- `ElectronicStopping` (from `electronic_stopping`) puts them together at E [MeV] for a projectile, an element, a
  phase and a compound correction. Refused: E ≤ 0, E above 10⁷ keV/u (SRIM returns 0 there), Z outside 1–92, C ≤ 0
  (LunaSRIM reads C = 0 as 1, the convention of SRIM's input file; here it must be said).

The branches meet at 2 keV/u, E_max, 1.2 E_max and 30 MeV/u: the stopping is continuous there (tested to 10⁻⁹ on both
sides), with a jump in its slope. LunaSRIM's `vbmath` (Visual Basic's x87 arithmetic, which LunaSRIM needs to
reproduce SRIM's output files digit for digit) is not ported: the C++ uses `std::pow`, `std::log`, `std::exp`, and
agrees with LunaSRIM to rounding (§9.1). Sums of polynomials are taken in LunaSRIM's order.

## 5. Nuclear stopping

SRIM-2013's pair-specific nuclear stopping (`NuclearStopping`, from `nuclear_stopping`): ZBL's reduced energy
ε = 32.53 M₂ E / (Z₁ Z₂ (M₁ + M₂)(Z₁^0.23 + Z₂^0.23)) (E in keV), s_n = ½ ln(1 + c₁ε)/(ε + c₂ ε^c₃ + c₄ ε^½) with the
four coefficients of the pair from `SNUC03.dat` for ε ≤ 30, ln ε/(2ε) above, and S_n = s_n 8.462 Z₁ Z₂ M₁ /
((M₁ + M₂)(Z₁^0.23 + Z₂^0.23)). The test checks the constants 32.53 and 8.462 against the definition
s_n = ε S_n/(π a² γ E) (to 0.2 %, their rounding). The fit and the high-energy form **do not meet at ε = 30**: S_n jumps
there, by +18 % for H on H, 5.7–9 % in size for the other pairs checked (p on O and Ta, α on C, ¹²C on C). In the
total stopping of p, α and ¹²C in H, C and Ta₂O₅ the jump is at most 1.2 % (α in H, at 20 keV), because S_n is small
there. `StoppingBreakpoints` gives the energy of each element's ε = 30 to the integrators (§8). ZBL's universal formula, which the first branch used, is
not kept (nothing uses it).

One deviation from LunaSRIM's arithmetic: `std::log1p(c₁ε)` where LunaSRIM writes `Log(1 + c₁ε)`. Both are the same
to 10⁻¹⁶/(c₁ε) relative (≤ 1.1 × 10⁻¹² in the reference rows, ε ≥ 10⁻⁴); with `log(1 + x)` the rounding noise at
ε < 10⁻¹² kept the range integral down to E = 0 from converging.

## 6. Compounds, phase and user stopping

**Bragg's rule as SRIM's tables apply it** (`Stopping`): S_e = Σᵢ fᵢ S_e,i · FadeOut(Fᵢ C), S_n = Σᵢ fᵢ S_n,i(Mᵢ),
fᵢ the atom fractions, Fᵢ the phase factor of element i in the material's phase, C the material's compound
correction. **C enters once.** LunaSRIM's `compound_stopping` applies C a second time below 10 MeV
(S_e × FadeOut(C) with c = 7), because it follows the list mode of SRIM's SR module, which does so; LunaSRIM's own
documentation (`web/static/physics.md`) says that the table mode of the SR module, TRIM and its energy-loss
calculation apply C once. This module follows the latter, so that its stopping, its energy loss and LunaSRIM's
`straggling.energy_loss` agree; the reference rows of compounds with C ≠ 1 are the once-applied sum. With C = 1
the two are identical.

**Phase.** A `Material` has `phase` (solid or gas; a liquid is a solid for SRIM) and `compound_correction`
(default 1 = Bragg). `MaterialFromFormula` and `MaterialFromFractions` set `DefaultPhase`: gas for a material of
one natural-gas element, solid for everything else, compounds included (CH₄ or CO₂ as a gas must be said). The phase
used is in the provenance text and in the header of every `azure2-stopping` output.

**Compound correction.** A number the user gives. SRIM lists corrections for its compounds in its compound
dictionary (`Compound.dat` of SRIM-2013, shown in SRIM's "Stopping in Compounds" window; LunaSRIM reads its densities
and stoichiometries but not the corrections); the look-up by name was not implemented.

**User stopping.** A user table or formula overrides the evaluation for one projectile (matched by charge and by
mass within 0.01 u): `UserStopping` with `energy_MeV`/`stopping` (interpolated linearly in ln E – ln ε, refused
outside its range) or a `formula`, a stated relative uncertainty and a provenance string. `ReadStoppingTable` reads
two columns with a unit line, `E[MeV]` or `E[keV]`, and `eps[eV/1e15atoms/cm2]`, `S[MeV*cm2/g]` (converted with the
material's mean atomic weight) or `S[keV/um]` (needs the mass density). This is also how a measured correction to
SRIM's stopping for one element or projectile enters.

## 7. Uncertainty

Every evaluation returns its relative 1σ uncertainty (`StoppingValue::rel_uncertainty`), the prior width of the
stopping-scale nuisance in Phase 5.

- **Electronic**: SRIM's published accuracy, from the table "Accuracy of SRIM Stopping Power Calculations" on the
  SRIM-2013 web page www.srim.org/SRIM/SRIM2011.htm (read 2026-11-11), column SRIM-2010 (the page calls it the current
  stopping calculation), the average deviation of all experimental data points from SRIM: **H ions 4.0 %** (≈ 8300
  points), **He ions 3.9 %** (≈ 6500), **Li ions 4.8 %** (≈ 1400), **Be–U ions 5.8 %** (≈ 9000); 74 %, 76 %, 72 %,
  58 % of the points lie within 5 % of SRIM. The page adds that without the points more than 25 % off, all numbers
  would be about 20 % lower. The average deviation is used as the 1σ (an interpretation: SRIM states no confidence
  level). `ElectronicRelUncertainty(Z₁)`.
- **Nuclear**: 10 %, **provisional**, not a quoted number (`kNuclearRelUncertainty`): SRIM's table is for the total
  stopping; no source for the nuclear part alone was found.
- A compound adds its constituents' electronic uncertainties linearly (fully correlated, one evaluation) and the
  nuclear part in quadrature. A user table carries its own.

## 8. What a target needs

On any ε(E) (`StoppingFunction`) and on a material:

- N(E_in, E_out) = ∫_{E_out}^{E_in} dE/ε, computed in u = √E (dE = 2u du, which makes the integrand finite at E = 0
  for a velocity-proportional ε), adaptive Gauss–Kronrod 7–15 to a relative tolerance of 10⁻¹⁰;
- E_out(E_in, N), the root of N(E_in, E) = N by safeguarded Newton steps (dN/dE_out = −1/ε); refused if N ≥ R;
- R(E) = N(E, 0);
- with ε → s ε, at s = 1: ∂N/∂s = −N, ∂R/∂s = −R, ∂E_out/∂s = −N ε(E_out).

**Break points.** Every function takes an optional list of energies where ε changes form; the integral is split
there. On a material the list is `StoppingBreakpoints`: 2 keV/u, the end of the light-ion fit and 1.2 times it
(helium's, inside the heavy-ion scaling, and 30 MeV/u for Z₁ ≥ 4), and each element's nuclear ε = 30. Measured: the
adaptive quadrature also converges across all of them without the break points, to the same value within
5 × 10⁻¹² (p, α, ¹²C in H, C, Ta₂O₅); the break points are kept because the nuclear stopping jumps (§5) and a jump is
where an adaptive rule may misjudge its error. Tested: the integral across each break equals the sum of its two
sides, with and without break points, to 10⁻⁹.

The tests check N, R and E_out against closed forms for ε = c E^p (p = 0.5, 0, 0.3, −0.8) to 10⁻⁹, the derivatives
against central differences, E_out and N against each other on Ta₂O₅, and the energy after a layer against
LunaSRIM (§9.1).

**Speed** (one thread, `-O2`, p in Ta₂O₅): one ε(E) 0.8–1.1 µs (first branch 3.1–3.5 µs); one E_out of a
10¹⁸ atoms/cm² layer at 0.3–1 MeV 0.6–0.9 ms (first branch 1.6–2.3 ms). E_out recomputes the range integral at each
Newton step; a cumulative table is the obvious speed-up when the engine needs it.

## 9. Comparisons

### 9.1 With LunaSRIM (the tests)

`stopping/tools/make_srim_reference.py` writes `stopping/tests/data/srim_reference.csv` (8889 rows) with LunaSRIM
1.0.0 itself: electronic stopping of p, d, ³He, α, ⁷Li, ¹²C, ¹⁴N, ¹⁶O in H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg,
Al, Si, Ca, Ti, Fe, Ni, Cu, Zr, Ag, Ta, W, Au, Pb (the gases in both phases), at 0.3 keV/u to 100 MeV/u covering every
branch, each boundary exactly and 10⁻³ on each side; nuclear stopping of the same pairs at reduced energies 10⁻⁴ to
10³ and 10⁻⁶ on each side of 30; TiN, Ta₂O₅, TiH₂, LiF, CaF₂, H₂O, CH₄ (gas), CO₂ (gas) with C = 1, 0.93, 1.08;
eight layers with `lunasrim.straggling.energy_loss`. The C++ tests read the file and need no Python. Largest
relative deviation C++ / LunaSRIM − 1 (tolerance 10⁻¹⁰ for the stopping powers):

| | p | d | ³He | α | ⁷Li | ¹²C | ¹⁴N | ¹⁶O |
|---|---|---|---|---|---|---|---|---|
| electronic, 25 elements (5921 values) | 4.3e-16 | 4.8e-16 | 0 | 0 | 4.5e-16 | 1.6e-11 | 3.3e-16 | 2.4e-15 |
| nuclear, 25 elements (2000 values) | 1.1e-12 | 1.1e-12 | 7.7e-13 | 7.7e-13 | 9.8e-13 | 8.5e-13 | 9.5e-13 | 1.0e-12 |

Compounds (960 rows): electronic 0 (bit-identical), nuclear ≤ 4.2 × 10⁻¹⁶. The nuclear 10⁻¹² is `log1p` (§5).

Layers (energy loss ΔE = E_in − E_out, and the straggling σ of both models), tolerance 10⁻⁵:

| layer | ΔE | σ TRIM | σ Bohr |
|---|---|---|---|
| p 300 keV, 10¹⁷ Ta | 2.2e-6 | 9.2e-7 | 1.3e-6 |
| p 1 MeV, 2 × 10¹⁸ Ta | 1.3e-6 | 1.0e-6 | 1.4e-6 |
| p 200 keV, 5 × 10¹⁷ TiN | 4.3e-6 | 8.4e-7 | 4.6e-7 |
| p 2 MeV, 5 × 10¹⁸ TiN | 1.2e-6 | 4.0e-7 | 2.2e-8 |
| α 1 MeV, 10¹⁸ C | 4.6e-6 | 2.9e-6 | 2.5e-6 |
| α 5 MeV, 10¹⁹ C | 1.8e-6 | 8.2e-7 | 1.2e-6 |
| ¹²C 12 MeV, 5 × 10¹⁷ C | 1.7e-6 | 4.8e-7 | 8.6e-7 |
| ¹²C 30 MeV, 3 × 10¹⁸ C | 2.3e-6 | 1.7e-6 | 2.1e-6 |

The layer differences are LunaSRIM's integrator, not the stopping: LunaSRIM tabulates S on 200 points per decade,
interpolates it log-log, takes S′ as a numerical difference and steps the energy and the variance with RK4 (2000
steps in the reference); this module integrates the exact ε to 10⁻¹⁰. Hence 10⁻⁵, not 10⁻¹⁰.

A coefficient column shifted by one in the reader (the proton fit read from column 12 instead of 11) fails 1001
electronic rows; a phase factor always 1, the compound correction applied twice, TRIM's mean Z rounded half up, or
no break points each fail their tests (shown once while writing them).

### 9.2 With ICRU 49 (NIST PSTAR / ASTAR), an independent check

NIST's tables (downloaded for comparison only, not committed) against SRIM-2013's electronic stopping, the gases as
gases. Largest SRIM/ICRU 49 − 1 in each energy band, with the energy (MeV) where it occurs when it matters; "C (am.)"
amorphous carbon, "C (gr.)" graphite (SRIM has one carbon). No test: the two are different evaluations, and these
numbers are their difference, not an error of the port.

| | 1–10 keV | 10–100 keV | 0.1–1 MeV | 1–3 MeV | 3–10 MeV |
|---|---|---|---|---|---|
| p H | −3.2 % (0.01) | −3.2 % (0.01) | −3.4 % (0.175) | +0.5 % | +0.4 % |
| p He | −84 % (0.002) | −27 % (0.01) | −2.2 % | −0.5 % | +0.2 % |
| p C (am.) | +3.7 % | +8.6 % (0.03) | +2.3 % | +2.4 % | +1.6 % |
| p C (gr.) | −5.1 % | −6.0 % (0.1) | −7.1 % (0.15) | +1.1 % | +0.7 % |
| p N | −5.1 % | −5.1 % | −2.4 % | −2.4 % | −1.0 % |
| p O | −14 % (0.01) | −14 % (0.01) | −1.2 % | −0.5 % | +0.1 % |
| p Ne | −33 % (0.0015) | −16 % (0.01) | +2.7 % | +0.5 % | +0.5 % |
| p Al | −4.2 % | −4.2 % | −1.8 % | +1.8 % | +1.1 % |
| p Ti | +2.0 % | −1.8 % | −2.5 % | −3.0 % | −1.5 % |
| p Au | +11.6 % (0.0015) | −2.5 % | +2.6 % | +1.1 % | +0.7 % |
| α H | +60 % (0.001) | +15 % (0.01) | −1.8 % | −2.1 % | +1.3 % |
| α He | +25 % (0.001) | +2.6 % | −4.9 % (0.175) | +1.1 % | +3.1 % |
| α C (am.) | −26 % (0.001) | −10 % (0.01) | −3.1 % | +4.1 % | +4.9 % |
| α C (gr.) | −27 % (0.001) | −14 % (0.01) | −9.9 % (0.1) | −4.2 % | +3.7 % |
| α N | +26 % (0.001) | +11 % (0.01) | −3.6 % | −1.4 % | +1.2 % |
| α O | −15 % (0.001) | −10 % (0.01) | −4.0 % | −1.9 % | +1.7 % |
| α Ne | −12 % (0.001) | −8.5 % (0.01) | −1.3 % | +2.8 % | +2.8 % |
| α Al | +43 % (0.001) | +8.8 % (0.01) | −3.4 % | +1.9 % | +2.3 % |
| α Ti | −5.9 % | −5.8 % | +3.2 % | +0.8 % | −2.8 % |
| α Au | +32 % (0.001) | +10 % (0.01) | −5.9 % | −3.0 % | +1.1 % |

Above about 0.1 MeV/u the two evaluations agree within a few per cent (graphite and α in Au up to 6–10 % at
0.1–0.2 MeV). At the lowest energies (below about 10 keV/u, the first one or two columns) they differ by tens of per
cent, most of all in the gases: these are the energies where both rest on few measurements; which one is closer to
the data there was not examined. The first branch's ICRU-49-form evaluation agreed with PSTAR/ASTAR to 0.1 % where it
used their coefficients (its table is in the history of this file).

## 10. Straggling

Two models (`StragglingModel`), as LunaSRIM's `straggling.py` has them; `kTrim` is the default.

- **TRIM** (`TrimSpecificVariance`, from `w_trim`): per unit path W = w²/12 with w = min(k E, w_B),
  w_B = 177 Z₁ (Z₂ N)^½ eV/Å^½ × ((1/β² − ½)/(1/β² − 1)) (Bohr with TRIM's relativistic factor; 177²/12 = 2610.75 is
  TRIM's rounding of 4πe⁴ = 2605.6 eV² Å², +0.2 %), k² = 3.2 × 12 N (Z₁ + Z₂)^{8/3}/(1000 M₁) with E in eV (the
  low-energy limit; it acts only below a few hundred eV per u), N in atoms/Å³, and Z₂ TRIM's mean target charge, an
  integer: Z = round(fᵢ Zᵢ + Z) over the constituents in order, rounding half to even (`TrimMeanTargetCharge`; TiN:
  11, then round(14.5) = 14). Per atom, dΩ²/dN = W/N (independent of N), in MeV² cm² (× 10⁻²⁸ from eV² Å²).
- **Bohr** (`BohrSpecificVariance`): dΩ²/dN = 4π z₁² e⁴ Z₂ (N. Bohr, Phil. Mag. 30 (1915) 581), e² = αħc, with the
  exact atom-fraction-weighted Z₂ (`MeanTargetCharge`). 0.2606 keV² per 10¹⁸ atoms/cm² for z₁ = Z₂ = 1.

The first branch's Chu and Yang–O'Connor–Wang forms, which had no data, are removed.

**Propagation** (B. Tschalär, Nucl. Instrum. Methods 61 (1968) 141 and 64 (1968) 237): the solution of
dΩ²/dx = W(E) − 2 S′(E) Ω², which LunaSRIM integrates step by step, is σ²(x) = S(E_out)² [σ₀²/S(E₀)² +
∫_{E_out}^{E₀} (dΩ²/dx)/S³ dE], computed in √E with the break points of §8. The mean energy uses the total stopping,
the spread the electronic straggling only (nuclear collisions give a low-energy tail). The beam's entrance spread σ₀ is
carried through the layer by S(E_out)/S(E₀) and added in quadrature (`StragglingWidths`: straggling, beam, total).
Tested: constant S gives Ω² + σ₀² exactly, a thin layer gives dΩ²/dN N, ε = c E^p against the closed form to 10⁻⁸,
and eight layers against LunaSRIM (§9.1).

## 11. Command line and Python

`azure2-stopping` (built in `<build>/stopping/`):

```
azure2-stopping stopping   --projectile a --material Ta2O5 --from 0.02 --to 10 --points 50
azure2-stopping range      --projectile p --material TiN --from 0.1 --to 5
azure2-stopping layer      --projectile p --material Ta2O5 --energy 1.0 --areal-density 1e18 --beam-sigma 0.0005
azure2-stopping layer      --ion-z 6 --ion-mass 12 --material C --energy 12 --areal-density 5e17 --straggling bohr
azure2-stopping provenance --projectile a --material CO2 --phase gas --compound-correction 0.93
```

`--phase solid|gas` overrides the default phase, `--compound-correction C` sets SRIM's factor, `--ion-z`/`--ion-mass`
give any ion. Every table starts with comment lines naming SRIM-2013 and LunaSRIM, the citation to give, and the
material, phase and compound correction.

Python: not bound in this branch. Binding the functions as a `pyazr.stopping` submodule needs one new file
`api/src/pybind_stopping.cpp` (array inputs, numpy outputs, as `pybind_physics.cpp`), its line in
`api/src/CMakeLists.txt` with `azure2::stopping` linked, and two lines in `api/src/pybind_api.cpp` (the declaration
`void BindStopping(py::module_ &m);` and the call `BindStopping(m);`); it is left for the wiring branch.

## 12. Provenance and licence

- **Provenance.** `StoppingProvenance(material, projectile)` and `azure2-stopping provenance` print: SRIM-2013 with
  its citation (`kSrimCitation`), LunaSRIM 1.0.0 as the implementation ported (`kSrimImplementation`), the data source
  (compiled in or the directory), the phase and compound correction used, per element the `SCOEF03.DAT` row, I, M,
  phase factor and the ε = 30 energy of its nuclear stopping, and the uncertainties with their sources.
- **The code** of `stopping/` is AZURE2's; it is a port of LunaSRIM's `stopping.py` and `straggling.py`, whose MIT
  licence permits it (copyright notice: LunaSRIM's LICENSE, J. Skowronski).
- **The data files** `SCOEF03.DAT` and `SNUC03.dat` are SRIM-2013's (J. F. Ziegler, J. P. Biersack, M. D. Ziegler),
  **not covered by AZURE2's licence**; they stay under the terms of their authors. They are in the repository on the
  PI's decision of 2026-11-11; whether the public repository ships them is the owner's decision
  (`stopping/data/srim/README.md`, which also records LunaSRIM's reading of SRIM's terms: free for personal,
  non-commercial use). Without them the library reads them from `AZURE2_SRIM_DATA_DIR`.
- The reference values in `stopping/tests/data/srim_reference.csv` are LunaSRIM's output for those data.
- NIST PSTAR/ASTAR tables were downloaded for the comparison of §9.2 only, into the scratch directory; nothing of them
  is committed.

## 13. Next

1. Wiring into the engine (Phase 5, `energy_grid.md` §6): `TargetEffect`'s stopping expression, its range-integral
   energy loss and Bohr straggling replaced by this module (a material, phase, compound correction and projectile per
   target in format 3), the stopping-scale nuisance with the prior of §7, the Python binding (§11), and the removal of
   ERYA and `SRIM2013.xml`.
2. Optional: the compound correction by name from SRIM's `Compound.dat` (it would be a third SRIM file to ship).
3. A source for the uncertainty of the nuclear stopping, to replace the provisional 10 %.
