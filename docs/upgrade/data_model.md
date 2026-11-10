# Phase 4: data are data — design

A fit compares, point by point, a measured value y_i ± δy_i at a beam energy E_lab,i (MeV) and angle θ_lab,i (deg)
with a model value m_i(θ, ν) computed from the level parameters θ and the apparatus nuisances ν (norms, energy shifts,
densities). Today the measured value, its kinematics, the channel functions at its energy, the model value, the
values of other points it computes and pointers to them all live in one object, `EPoint`, and a "data type" is a set
of flags threaded through `EPoint`, `ESegment`, `EData`, the matrix functions and the gradient. This note describes
what the code does now, with file and line, and designs the Phase 4 target of `ROADMAP.md`: a `Measurement` that holds
data only, per-evaluation result arrays, and a model expression that says what each data set measures. The gate is
bit-identical outputs on the 179 test projects × 3 bases and unchanged fits; each step below says why it keeps the
bits, or which it may move and how that change is classified.

Abbreviations: `EP` = `src/EPoint.cpp`, `EPh` = `include/EPoint.h`, `ES` = `src/ESegment.cpp`, `ESh` =
`include/ESegment.h`, `ED` = `src/EData.cpp`, `EDh` = `include/EData.h`, `AC` = `src/AZURECalc.cpp`, `AG` =
`src/AZUREGrad.cpp`, `GM` = `src/GenMatrixFunc.cpp`. Line numbers are those of `feature/upgrade` at eb79d95a.
"Verified" means read in the code; "inferred" means concluded from the code without a run.

## 1. What exists

### 1.1 What a point, a segment and the data set hold

`EPoint` (`EP` 3132 lines, `EPh` 820 lines; about 96 data members, `EPh` 655–817) mixes six kinds of state:

| kind | members (`EPh`) | written |
|---|---|---|
| input datum | `lab_energy_`, `original_energy_`, `lab_angle_`, `lab_crosssection_`, `lab_dcrosssection_`, `bin_low_lab_`, `bin_high_lab_`, `fileRecord_`, `coefficient_order_`, `coefficient_absolute_` | load (`EPoint(DataLine, ESegment*)`, `ES` 461–583) |
| segment settings copied into every point | `entrance_key_`, `exit_key_`, `segment_key_`, `is_differential_`, `is_phase_`, `is_thm_`, `is_ang_dist_`, `is_analyzing_power_`, `is_total_`, `j_value_`, `l_value_`, `phase_s_value_`, `isUPOS_`, `secondaryDecayL_`, `Ic_`, `delta_`, `max_ang_dist_order_`, `crossSectionComponent_`, seven `thm_*` pointers | load (`EP` 130–165, `ES` 1278–1287) |
| kinematics derived from the datum | `cm_energy_`, `excitation_energy_`, `cm_angle_`, `cm_crosssection_`, `cm_dcrosssection_`, `angleKinFactor_`, `crossSectionKinFactor_`, `bin_low_cm_`, `bin_high_cm_`, `cm_factor_`, `photo_qgamma_`, `photo_mcn_`, mutable `labAngleMax_`, `labAngleWhy_` | load, **and again on every change of an energy shift** (`ES` 1100–1275) |
| energy-dependent tables (depend on E_c.m., θ_c.m., radii only) | `edep_` (L⁰, √P, e^{iω}, e^{−iφ} per (J group, channel); THM vertex pieces), `legendreP_`, `thm_ps_`, `geofactor_` (π/k²), `sfactorconv_`, `screen_k_`, `coulombamplitude_`, `coulombdirect_`, `coulombexchange_`, `ec_amplitudes_`, `ec_energies_`, `eDependentEnergy_`, `eDependentValid_` | load; on a shift; on a re-placed sub-point grid; shared copy-on-write between clones and mapped points (`SharedTable`) |
| apparatus | `targetEffectNum_`, `stoppingPower_`, `targetThickness_` (per ν: fitted density), `targetLabEnergy_`, `sub_grid_spec_`, `udr_kernel_`, `udr_mass_`, `target_blend_weight_`, `target_fallback_only_`, `target_fallback_depth_`, `screen_ue_` (per ν), `screen_param_`, `bohr_z1_`, `bohr_z2_`, `integrationPoints_` (the sub-points, themselves `EPoint`s, re-placed per θ), `parentData_` | load; per ν in the hand-over (`ED` 5329–5340) |
| topology | `is_mapped_`, `energy_map_`, `local_mapped_points_` (raw pointers to other points), `own_e_dependent_` | `MapData` at load and after every shift change (`ED` 5506) |
| **results of an evaluation** | `fitcrosssection_`, `fitE1crosssection_`, `fitE2crosssection_`, `natural_terms_`, `analyzing_power_`, `ang_dist_a0_`, `angularDists_`, mutable `ec_shifted_` (memo) | every evaluation |
| scoped flags | `doppler_applied_`, `target_fallback_now_` | inside one call (RAII scopes, `EPh` 333–345, 608–620) |

`ESegment` (`ES` 1852, `ESh` 609 lines; about 100 members, `ESh` 496–606) holds the segment line's settings (keys,
ranges, the isDiff-derived flags `isdifferential_`, `iscmdifferential_`, `isphase_`, `isTHM_`, `isAnalyzingPower_`,
`isTotalCapture_`, `crossSectionComponent_`, `isCoefficients_`, the keyword values), its points by value (`points_`),
raw pointers to its component segments (`componentSegments_`), and per-evaluation state: `dataNorm_` (fit value, and
the profiled THM norm written during the χ², `ES` 336, 342), `energyShift_`, `energyShiftSqrt_`, `shiftFnValues_`,
`thickness_` (ν), `lastEnergyShift_*` (the shift the points carry), `shiftTables_` (the channel tables of the last
three shift states, `ESh` 494), `shiftRows_` (scratch), `segment_chi_squared_` (result).

`EData` (`ED` 5819, `EDh` 506 lines) owns the segments (`segments_`), the component segments
(`componentSegments_`, a `std::deque<ESegment>` of full copies of their base segment, points included, `ED`
3916–3977), the target effects (with the fitted densities, Doppler temperatures and screening potentials inside),
the parameter offsets, the penalty totals of the last stored evaluation (`parkPenalty_` … `priorCovChi2_`), and
three objects shared between all its clones: `covGroups_`, `nodeCache_`, `priorCov_` (`ED` 5618–5620). It also owns
`gridNuc_`, a private `CNuc` copy built lazily for the sub-point grids (`ED` 1279).

A component of a SUM, RATIO or total-capture segment is therefore a complete second segment holding a second copy of
every datum, with another exit pair (`ED` 3920–3940); its data values are never used.

### 1.2 Writes from one object into another during an evaluation

| # | who writes | what, into whom | when | under threads |
|---|---|---|---|---|
| W1 | a host point (`EPoint::Calculate`, `EP` 2148–2225) | the cross section of every point mapped onto it (`CalculateCrossSection(mapped)`; a new T-matrix when the exit differs, `EP` 2208–2225); for THM, `CalculateTHMCrossSection(mapped)` (`EP` 2177–2183) | each evaluation | safe by construction of the mapping: only unmapped points are work items (`ED` 5778–5785); a mapped point has one host and is never a host itself (`ED` 4180–4195, 4237–4246), so the hosts write disjoint sets |
| W2 | a host with a sub-point effect | sub-point j of each mapped point (`Calculate(cnuc, cfg, this, j)` → `parent->GetLocalMappedPoint(i)->GetSubPoint(j)`, `EP` 2213–2218), then the mapped point's integral (`EP` 2329–2341) | each evaluation | as W1 |
| W3 | `ESegment::CalculateTheoreticalCrossSection` (`ES` 1698–1822) and its 19 callers | the SUM, RATIO or abundance-weighted value **over the base point's own cross section** (`AC` 558–560 and the loops in `AC` 112, 264, 297, 634, 692, 1761; `api/src/AZUREAPI.cpp` 372, 482, 593, 1002; `src/AZURECalcMCMC.cpp` 102, 189, 900; `AG` 705, 927, 979, 1803; `ED` 4877); the per-nucleus terms into the base point (`SetNaturalTerms`, `ES` 1633, 1742, 1773) | each evaluation, serial loop after the parallel pass | safe (serial); but the base point's leaf value is destroyed, so a reader cannot tell a leaf value from a model value |
| W4 | `ESegment::TransformedValue` (`ES` 1605–1696) | the transformed values into the base point's sub-point slots, then the integral, then the old values back (`ES` 1681–1691) | each evaluation | safe only because nothing else reads the base point meanwhile |
| W5 | the shift step, `ESegment::UpdatePointEnergyWithShift` (`ES` 1100–1275) | a mapped point takes its host's tables (`CopyEDependentValuesFrom`, `ES` 1172–1175, 1262–1265; `EP` 1813–1831) | when a shift changes | safe by grouping: a host and its same-segment mapped points are one work item, host first (`ES` 1077–1098; `ED` 5466–5500); tables are copy-on-write (`include/SharedTable.h`) |
| W6 | `EData::FillEnergyShiftsFromParams` (`ED` 5329–5568) | the base segment's shift into its component segments (`ED` 5408–5418, 5428–5446) | when a shift changes | serial |
| W7 | `EData::ReplaceSubPointGrids` (`ED` 1228–1454) | a host lays out the sub-point grid of every point mapped onto it (`layOut(mp, …)`, `ED` 1394–1404); tables published to and adopted from `nodeCache_`, shared by every clone (`ED` 1362–1371) | each evaluation (`ED` 5495) | parallel by host; the cache is sharded with a mutex, and a published table is never written again (`include/NodeTableCache.h`) |
| W8 | `AZURECalc::GetPooledCNuc/EData` (`AC` 2012–2087) | a pooled working copy is taken from a LIFO stack of max(4, threads) clones, refilled from p, evaluated, and returned; it keeps every applied shift, mapping, grid and table | each evaluation of a fit | one evaluation per copy at a time; which copy, and therefore which history, depends on scheduling |
| W9 | `EData::Clone` (`ED` 5598–5724) | copies segments and points (tables by pointer), re-points component pointers and mapped-point lists into the copy, shares `covGroups_`, `nodeCache_`, `priorCov_`; throws if a copied host still points outside the copy (`ED` 5700–5716) | pool set-up, output snapshots, `AZURE_POOL_CHECK` | serial |
| W10 | the objective | the profiled THM norm into the segment (`ES` 336, 342; `src/EDataThm.cpp` 698), the segment χ² and the penalty totals into `EData` when `storeTotals` (`src/Objective.cpp` 79) | each stored evaluation | serial |
| W11 | the ν hand-over | fitted densities into each point's thickness, U_e into each point, Doppler T and sample thickness into effects and segments (`ED` 5329–5340) | each evaluation | serial, before the points |

Probe copies are made inside the evaluation as well: with edge blending or an automatic tolerance a target-effect
point copies itself, sub-point vector included, and calculates the copy (`EP` 2238–2245, 2285–2293, 2331–2340).

### 1.3 State carried from one evaluation to the next

The fit depends on the thread count (HANDOFF, top entry); calculations do not. A calculation evaluates `data()`
itself once; a fit evaluates pooled copies (W8), whose state is the result of their own history. Every item below is
state a copy carries between evaluations. Each is harmless only if what it holds is a function of the current (θ, ν)
alone.

| item | where | a function of (θ, ν) alone? |
|---|---|---|
| H1 shift "unchanged since the last call": points re-energised only if the shift differs from `lastEnergyShift_` | `ED` 5390–5397, 5427–5432 | yes, if re-energising is a pure function of the shift (H2–H4) |
| H2 the channel tables of the last three shift states, reused by exact key | `ES` 1022–1050, `ESh` 494 | yes (exact keys) |
| H3 **sub-point grids kept when every node of the new grid lies within 10⁻¹² (relative) of the old one** | `ED` 1316–1321 | **no**: the node energies of a copy are those of the last θ′ for which the grid moved by more than 10⁻¹², so two copies at the same θ can hold grids that differ in the last bits. Inferred, not run. Active by default (`replaceSubPointGrids = true`, `src/Config.cpp` 44) on every project with `useAdaptiveGrid=1` and an integrated target, T2 included |
| H4 a lab angle pushed beyond the kinematic maximum by a shift keeps "its last valid value", and the relativistic cross-section factor its last valid factor | `ES` 1143–1152, `EP` 1026–1031 | **no**: the last valid value is the copy's history. Only when the warning is printed |
| H5 the energy-dependent values are not recomputed while E_c.m. is unchanged | `EP` 1702 | yes |
| H6 the external-capture amplitude at the last shifted energy | `EP` 664–683 | yes (exact energy; the level energy is set at set-up, `src/CNucExternalCapture.cpp` 131–186) |
| H7 tables adopted from the shared `nodeCache_` or from the old sub-points at bit-identical (segment, E, θ) | `ED` 1339–1371 | yes, if `CalcEDependentValues` is pure |
| H8 `gridNuc_`, cloned on first use and refilled from the current θ | `ED` 1279–1289 | yes |
| H9 the relvector covariance model frozen at the first evaluation that reaches `EnsureFactored`, shared by all copies | `ED` 4802, `src/DataCovariance.cpp` 520–527; seeds prepared in `AC` 827 | yes in a CLI fit (the first evaluation is the start); in a pyazr session it is whatever was evaluated first |
| H10 component segment keys from two process-wide counters, both starting at −1000 | `ED` 3930, 3970 | not thread-related: a fixed-angle component and a plain one can share a key; keys continue across pyazr sessions; keys are part of the `nodeCache_` and `ECAmplitudeCache` keys |

Outside the data path, seen in passing and not examined: the thread-local A-matrix and THM functions reused across
points (`EP` 2171, 2197), the level-space caches keyed by `JGroup*` and level stamps (`src/ChannelSpace.cpp` 36–120),
the process-wide shift-derivative memo (`src/CNuc.cpp` 2593).

H3 is the one candidate in the data path for the fit's thread dependence that I can name; H4 can contribute only in
fits that print its warning. A test that decides it: evaluate θ_A on a fresh clone, and on a pooled clone after
evaluating θ_B ≠ θ_A, and compare the per-segment χ² bit for bit (`AZURE_POOL_CHECK`, `AC` 704–745, does half of
this already).

### 1.4 Where lab → c.m. happens

- Energies: `EPoint::ConvertLabEnergy` (non-relativistic E_c.m. = E_lab m₂/(m₁+m₂); relativistic through
  `Kinematics::CMFromLab`), `ConvertDecayEnergy` (β-delayed: the light particle's energy, c.m. of the exit pair),
  `ConvertExcitationEnergy` (`EP` 768–828). Called at load (`ES` 535–543) and on every shift change (`ES` 1123–1130).
- Angles: `ConvertLabAngle` (elastic and inelastic, relativistic branch `RelativisticLabToCM`), `ConvertLabAngleGammas`
  (photon aberration) (`EP` 830–1021). Lab-frame differential data only (isDiff 1); not for c.m. data (4, 7), phase
  shifts, coefficients.
- Cross sections and their errors: `ConvertCrossSection`, `ConvertCrossSectionGammas` multiply y and δy by the
  Jacobian dΩ_lab/dΩ_c.m. (`EP` 1023–1110); skipped for RATIO data (`ES` 560, 571, 1154, 1162). The same Jacobian is written a
  second time in `CalculateCrossSectionConversionFactor` (`EP` 1061–1089), used by RATIO and SUM frame factors.
- The conversion sequence itself is written twice: in `ESegment::Fill` (`ES` 535–575) and in
  `UpdatePointEnergyWithShift` (`ES` 1123–1170).
- **With a fitted energy shift the c.m. datum depends on ν.** The shift moves E_lab, hence θ_c.m. and the Jacobian, and
  `cm_crosssection_` and `cm_dcrosssection_` are rewritten (`ES` 1132–1165). The residual (f − n y_c.m.)/(n δy_c.m.)
  (`src/Objective.cpp` 58–59) therefore has a weight 1/δy_c.m.(ν) that moves with the shift. Sub-points move rigidly
  by the parent's c.m. shift (`ES` 1188–1230) and are never angle-converted (`ES` 1236–1247).
- The S-factor conversion (`sfactorconv_`, `EP` 1441) is used only by the writer (`ED` 2666–2722); data are always
  cross sections in b or b/sr.

### 1.5 The data kinds the engine accepts

isDiff decoding: data lines `ES` 25–150, test lines `ES` 155–250 (the two tables differ, see the table); refusals
`ED` 374–410 (data), 726–742 (test). Keywords: `include/SegLine.h` 105–390, `include/ExtrapLine.h` 72–130.
Evaluation of observables: `GM` (σ 128–258, total 103–119, dσ/dΩ 270–425, coefficients 427–436, Coulomb–nuclear
interference 516, phase shift 566, A_y 573–589, screening 25–36), `src/THMMatrixFunc.cpp` (THM). Composition: `ES`
1605–1822; its adjoint a second time in `AG` 697–1100.

| kind | written as | parsed | predicted where | note |
|---|---|---|---|---|
| σ(E), angle-integrated | data isDiff 0, test 0 | `ES` 25 | `GM` 128–258 | identical entrance pair ×2 (`GM` 230–257) |
| dσ/dΩ, lab | data 1, test 1 | `ES` 46 | `GM` 270–425, 516 | data converted to c.m. (§1.4) |
| dσ/dΩ, c.m. | data **4**, test **5** | `ES` 50, 174 | as above | no conversion |
| phase shift δ_{J l s} | 2 (+ J, l; `spin s`) | `SegLine.h` 105, `ES` 54 | `GM` 566 | checked against the level scheme (`ED` 2042–2125) |
| total capture | data **3**, test **4** (exit −1) | `ES` 63, `ED` 531–556 | base segment on the first capture pair + one component segment per other capture pair, SUM | |
| σ_E1, σ_E2 (multipolarity) | data 5, 6 | `ES` 66–71 | `GM` 154–196 (E1/E2 T-matrix sums, particle and external capture) | angle-integrated only; the only "multipolarity mask" today |
| A_y (vector analysing power), c.m. | 7 | `ES` 44 | `GM` 573–589 | the target average is σ-weighted |
| Legendre coefficient a_k = A_k/A₀ (or A_k with `absolute 1`) | data 9 (angle column = k) | `ES` 81, 500–512 | `GM` 427–436 | test segments use 3 instead |
| angular distribution (all A_k), test only | test **3** + max order | `ExtrapLine.h` 26 | `GM` 223–425 | |
| THM (half-off-shell) | 10 + any code above (19 refused) | `ES` 39–41 | `src/THMMatrixFunc.cpp` | profiled norm and group background in the χ² (`ES` 312–357, `src/EDataThm.cpp`) |
| polarisation × σ | 8, 18 | refused (`ED` 397) | — | |
| SUM of reactions | `1 0 n e x θ [s]` | `include/ComponentList.h` | `ES` 1747–1777 | terms × scaling × frame factor J_base/J_k (`ES` 1487–1493) |
| RATIO of two reactions | `1 1 n e x θ [s]` | as above | `ES` 1778–1819 | data dimensionless (no Jacobian); each term's own Jacobian (`ES` 1501–1522); a component angle fixes the denominator's θ (`ED` 3948–3966) |
| neutron total (optical theorem) | `total 1` | `SegLine.h` 328 | `GM` 103–119 | exit = entrance |
| natural / isotopic sample | `sample <name>` (`natural`) | `SegLine.h` 316 | `ES` 1730–1777 | every term × abundance a_{s,n}, a fit parameter (T7) |
| transmission T = e^{−nσ_tot} | `transmission n [vary σ]` | `SegLine.h` 334 | `ES` 1605–1696 | transform at each sub-point, then the kernel (`ED` 3819–3870) |
| self-shielded yield Y₀ = (s/s_t)(1 − e^{−n s_t}) | `yield n [vary σ]` | as above | as above | a SAMMY neutron-sample yield, with an added total term (`ED` 3856); not a thick-target charged-particle yield |
| unobserved primary, observed secondary γ | positional UPOS block (L, I_c, δ) | `SegLine.h` 198–209 | `GM` 270–425 | |
| β-delayed spectrum | entrance pair of type 20 | `<pairs>` | `GM` 140, 172 | decay energies (`EP` 801) |
| S factor | — | — | output only (`ED` 2666–2722) | EXFOR S-factor data are converted to σ in Python (`pyazr/nds.py`) |
| electron screening | `<targetInt>` screening block | `src/TargetEffect.cpp` 17 | `GM` 25–36 | a factor on σ at each sub-point |
| target integration, Gaussian, beam profile, UDR, Doppler, Q-coefficients | `<targetInt>` | `src/TargetEffect.cpp` | `EP` 2350–3040 | apparatus; folded inside the point |
| shifts, norms, groups, covariance | segment keywords | `SegLine.h` | ν and C, not the prediction | |
| model-only: scattering lengths, abundance and Park constraints | `<scatteringLengths>`, `<abundances>` | `src/ScatteringLength.cpp`, `src/CNuc.cpp` | `ModelOnly()` terms (`src/Objective.cpp` 109–169) | no beam energy |

### 1.6 The same logic written more than once

- The model loop "calculate the segment's value, store it in the point": 19 call sites (W3).
- The residual (f − n y)/(n δy): `src/Objective.cpp` 58, `AC` 116, 265, 298, 459, 476, 832, 1762, `ES` 326–347, the
  writer.
- Composition (SUM, RATIO, abundances, transforms): forward in `ES` 1605–1822, adjoint in `AG` 697–1100 and
  1108–1380.
- The lab → c.m. sequence (§1.4) twice; the particle Jacobian twice; the Rutherford cross section and the S-factor
  conversion again in Python (`pyazr/nds.py` 386 and its unit conversions).

## 2. The target

### 2.1 Three kinds of state, three owners

| state | type | owner, lifetime | written |
|---|---|---|---|
| what was measured, and what the model must predict | `Measurement` | built once at load, `std::shared_ptr<const Measurement>`, shared by every working copy (never cloned) | never after load |
| what depends on ν | `MeasurementFrame` (E_c.m., E_x, θ_c.m., Jacobian, y_c.m., δy_c.m. per point) and, in Phase 4 still, the calculation units (`EPoint`s) with their tables and sub-point grids | the working copy (Phase 4), the stage's tables (Phase 6) | in the single-threaded hand-over when ν or θ changes, never during the point pass |
| what an evaluation produces | `Results` (leaf values per calculation slot, model value per measured point, per-nucleus terms, profiled norms, χ² per measurement) | the caller of the evaluation; reused storage is poisoned on entry under `AZURE_POISON_RESULTS=1` | by the evaluation only, each slot by exactly one task |

The rule that follows: no object writes into another object; a calculation task writes only the result slots the plan
assigns it, and reads anything set up before.

### 2.2 Types

All in `experiment/` (namespace `azure2::experiment`), plain values, no engine headers:

```
struct Reaction      { int entranceKey; std::vector<int> exitKeys; };     // exitKeys: one, or every capture pair
enum class EnergyIn  { BeamLab, DecayLab, CM };                              // what column 1 is
enum class AngleIn   { Lab, CM, LegendreOrder, None };                       // what column 2 is

struct Measurement {
  int key;                                  // segment number in the file (1-based), the label of its outputs
  Reaction reaction;
  EnergyIn energyIn; AngleIn angleIn;
  bool convertValue;                        // y, dy carry dOmega_lab/dOmega_cm (lab differential, not a ratio)
  // input, as read (MeV, deg, b or b/sr or the quantity's unit), one entry per point inside the segment window
  std::vector<double> E, theta, y, dy, binLo, binHi;
  std::vector<int> record;                  // data-file row (covariance files index by it)
  // derived at load at the nominal shift (MeV, deg, 1)
  MeasurementFrame frame0;
  Prediction prediction;                    // the model expression, §2.3
  ApparatusRef apparatus;                   // <targetInt> line, shift and norm specification, covariance group
  Provenance provenance;                    // file, segment line, checksum
};

struct MeasurementFrame { std::vector<double> Ecm, Ex, thetaCm, jacobian, ycm, dycm; };
```

`MeasurementFrame frame(const Measurement&, const ShiftState&, const Masses&)` is the one function that applies a
shift and converts lab to c.m. (§1.4's two sequences become this one). `frame0` is its value at the nominal shift.
A measurement without an active shift uses `frame0` for the whole run; that is "converted to c.m. at load". A
measurement with a fitted shift gets a new frame in the hand-over whenever its shift changes (decision M3).

```
struct Results {
  std::vector<double> slot;                 // leaf values: one per (calculation unit, quantity) -- sigma, sigma_E1, sigma_E2, A_y, ...
  std::vector<double> model;                // m_i, one per measured point, all measurements concatenated
  std::vector<double> natural;              // per point and nucleus, for the decomposition output
  std::vector<double> norm, chi2;           // per measurement: the norm used (profiled for THM) and its chi^2
};
```

### 2.3 Model expressions

`Prediction` is a tree whose leaves are observable selectors and whose nodes combine their values point by point.

```
struct AmplitudeFilter { unsigned multipoles; /* E1, E2, M1, ... bits; 0 = all */
                         std::optional<double> J; int parity; int l; double s; enum { All, Internal, External } part; };
struct ObservableSelector {
  ObservableKind kind;   // Sigma, DSigmaCM, Ay, PhaseShift, LegendreA, AngularDistribution, TotalOptical,
                         // BetaDelayed, UposGamma, Thm, (new) Rutherford
  int entranceKey, exitKey;
  AmplitudeFilter filter;
  KindOptions options;   // J, l, s of a phase shift; Legendre order k and absolute; UPOS L, I_c, delta; THM options;
                         // a fixed angle (deg, c.m.) when the term does not use the point's
};
```

Nodes, each one class with `value` and `backward` (the cotangents of its children from its own), tested alone:

| node | value at point i | parameters |
|---|---|---|
| `Leaf` | the leaf's slot in `Results` for the calculation unit of point i | — |
| `Sum` | v₀ + v₁ + … in the order of the children (left to right, as `ES` 1749–1774) | — |
| `Scale` | (c · J_i) · v, with c a constant and J_i an optional per-point array (the frame factor J_base/J_k) | — |
| `Weighted` | a_{s,n} · v | abundance a_{s,n} (θ) |
| `Ratio` | v_num / v_den · (J_den / J_num), NaN when v_den = 0 (`ES` 1793–1816) | — |
| `Transform` | T = e^{−n s} or Y₀ = (s/s_t)(1 − e^{−n s_t}) applied at every sub-point of the base unit, then the base unit's kernel (`ES` 1660–1691) | sample thickness n (ν) |

The normalization stays outside the tree, on the data side of the residual (decision M5).

**Leaves and the apparatus in Phase 4.** A leaf's value is today's point value after its own fold: the calculation
unit (`EPoint` with its sub-points) computes the observable at the sub-points and folds them with its
`<targetInt>` line (`EP` 2297–2346). The fold therefore sits inside the leaf (decision M4). The one exception is
`Transform`, which must act between the sub-point values and the kernel: it reads the sub-point slots of its children
and calls the base unit's kernel on a local vector (`EPoint::IntegrateKernelsImpl` already takes a value reader,
`EP` 2793); no slot is overwritten (W4 goes). Phase 5 lifts every fold into an explicit `Fold(kernel)` node above a
bare leaf; the leaf selectors and the combining nodes stay as they are.

**What a filter means.** A filter selects amplitudes before they are squared: σ with "E1" is |Σ_{E1} T|² summed over
the incoherent labels, so a sum of filtered leaves is not the unfiltered σ where filtered and unfiltered amplitudes
interfere. For angle-integrated σ there is no interference between different (J, l, l′) or multipolarities after the
angular integration, and E1 + E2 + … equals σ; for dσ/dΩ it does not. Phase 4 keeps filters on angle-integrated leaves
(as `GM` 154–196 does now), decision M8.

### 2.4 How the present evaluator computes through them

An evaluation in Phase 4:

1. Hand-over (serial, as now): `FillCompoundFromParams`, norms, `FillEnergyShiftsFromParams` (new frames for the
   measurements whose shift changed, re-energised units, re-placed grids), and the **calculation plan**.
2. The plan: a list of tasks, each a host unit and the units mapped onto it (today's `MapData` and
   `MapComponentSegments` rules, moved unchanged), each unit owning a fixed range of `Results::slot`. The plan is
   built when the mapping is built today (load, and after a shift change), and checked when built: every slot has
   exactly one writer.
3. Point pass (parallel, `schedule(dynamic,1)`, largest first, as `ED` 5765–5819): each task runs today's
   `EPoint::Calculate`, `GenMatrixFunc`, `THMMatrixFunc` and the fold, writing into its slots.
4. Tree pass: for every measurement and point, m_i = prediction.value(i). Each m_i depends only on slots, so the
   order of points does not matter and the pass may run in parallel by measurement.
5. The objective reads `Results::model`, the frame's y_c.m. and δy_c.m., and the norms (the residual becomes one
   function used by every path, §1.6).

Nothing in `GenMatrixFunc`, `AMatrixFunc` or the folds changes except where a value is stored: `point->SetFitCrossSection(x)`
becomes a write to the unit's slot in the `Results` the task was given. The lazy path used inside an enclosing
parallel region (`pre = false`) runs the plan's tasks serially and then the tree pass; `ScopedPointMappingOff` is a
plan with no mapping.

### 2.5 Library placement

`experiment/` becomes a CMake target `azure2_experiment` that links only `azure2_physics` (masses, constants,
kinematics), like `physics/` and `io/`, from its first commit (decision M1). (Relaxed by decision F11 of `project_format3.md`: it may
also include `io/`'s plain record structs, `azure2/io/Record.h`, the header-only target `azure2_io_record`.) It holds `Measurement`, the frames, the
expression nodes, the keyword → tree table and `Results`. The engine (`src/`) holds the adapter: the calculation
units and the plan, i.e. today's `EPoint`/`EData` machinery, which a leaf reads through its slot. `io/` parses the
segment line into a plain `SegmentSpec` (today's `SegLine`, moved and made pure as `project.md` §2 intended); the
tree is built from it in `experiment/` (decision M13).

### 2.6 What Phase 4 does not do

- No energy grid, no kernels as sparse rows, no explicit `Fold` nodes, no stopping module (Phase 5).
- No removal of pools, `EData::Clone`, mapping, copy-on-write tables or `NodeTableCache` (Phase 6). Phase 4 makes the
  pooled copies hold only ν-state and makes that state a function of (θ, ν) (step S2), which is what Phase 6 needs to
  delete them.
- No new comparison frame for shifted lab data (decision M3).
- The THM experiment machinery (profiled norms, coherent backgrounds) keeps its code; only its stored results move to
  `Results`.

### 2.7 Why the bits stay

- Moving a value from a member to a slot stores the same double; reading it reads the same double.
- The plan performs the same calculations as the mapping (the same host computes the same mapped points with the same
  matrices); only who holds the pointers changes.
- Each node's arithmetic is written in today's order: `Sum` adds the base first and the terms in component order;
  `Scale` forms (c_k · J_base/J_k) · v_k and then multiplies by the abundance (`ES` 1762–1767); `Ratio` forms
  v_b / v_d · (J_d / J_b) (`ES` 1814). Reordering any of these can move the last bit, so each node has a test against
  the old expression on recorded values.
- The frame function is the old sequence with the duplicate removed; the relativistic branch keeps its fallback until
  S2 decides it.

Where a step can move bits it says so in §6 and names the classification.

## 3. Every present kind as a tree, and the text form

*The text form of this section is superseded by `project_format3.md` (2026-11-09); the trees stand.*

Notation: `σ(a→b)` is a leaf, `[E1]` a filter, `·` a `Scale`, `w(s)` a `Weighted` by the abundances of sample s, J_k the
term's frame factor.

| kind | tree |
|---|---|
| isDiff 0 | `σ(a→b)` |
| isDiff 1 | `dσ/dΩ_cm(a→b)` with y converted by the frame |
| isDiff 4 (data), 5 (test) | `dσ/dΩ_cm(a→b)` |
| isDiff 2 | `δ(a→a; J, l, s)` |
| isDiff 3 (data), 4 (test) | `Sum(σ(a→γ₁), σ(a→γ₂), …)` over the capture pairs of the entrance nucleus, in pair order |
| isDiff 5, 6 | `σ(a→b)[E1]`, `σ(a→b)[E2]` |
| isDiff 7 | `A_y(a→b)` |
| isDiff 9 | `a_k(a→b)` or, with `absolute 1`, `A_k(a→b)` |
| test isDiff 3 | `A_{0..K}(a→b)` (angular distribution) |
| isDiff 10 + c | `THM(code c, THM options)(a→b)` |
| SUM | `Sum(base, c₁J₁·t₁, c₂J₂·t₂, …)`, t_k the component leaf (with its own filter, a fixed angle if given) |
| RATIO | `Ratio(base, t₁)` with the two frame factors |
| `total 1` | `σ_tot(a→a)` |
| `sample s` | `Sum(w_s(n₀)·base, w_s(n₁)·c₁J₁·t₁, …)` |
| `transmission n` | `Transform_T(n; Sum(terms))` |
| `yield n` | `Transform_Y(n; Sum(reaction terms), σ_tot(a→a))` |
| UPOS | `γ_UPOS(a→b; L, I_c, δ)` |
| β-delayed entrance | the same leaves with `EnergyIn = DecayLab` |

**Text form (decision M9).** Today's keywords are enough to read every existing file, and they remain the canonical
form for these kinds. A text form is added at the end of Phase 4 so that a new combination needs no code:

```
predict "<expr>"
expr   := leaf | "sum(" expr {"," expr} ")" | "ratio(" expr "," expr ")" | number "*" expr
        | "abundance(" sample ")*" expr | "transmission(" n "," expr ")" | "yield(" n "," expr "," expr ")"
leaf   := kind "(" pair "->" pair {";" option} ")"
kind   := sigma | dsdo_cm | ay | phase | ak | Ak | sigma_total | upos | rutherford
option := mult=E1|E2|M1|... | J=1/2+ | l=0 | s=1/2 | part=internal|external | angle=<deg c.m.> | k=<order>
```

Examples: `predict "ratio(dsdo_cm(1->1), rutherford(1->1))"`, `predict "sum(sigma(1->2; mult=E1), 0.5*sigma(1->3))"`.
The `predict` keyword replaces isDiff's meaning for the line (the isDiff field still says the frame of the data
columns); a line with both `predict` and a composite block is refused. `io/` parses the string into a plain spec and
reports errors with line and column; `experiment/` builds the tree; pyazr reads it through the binding (no Python
parser); the writer writes keywords where the tree has a keyword form, `predict` otherwise. Folds stay out of the
grammar until Phase 5.

## 4. A new kind of datum after Phase 4

**Energy-dependent: ratio to Rutherford.** A new leaf `rutherford(a→a)` = dσ_C/dΩ_c.m., |f_C(θ)|² (Mott-symmetrised
for identical particles), from the Coulomb amplitude each unit already holds (`coulombamplitude_`, `EP` 1892–1940).
The datum is dimensionless and frame-independent for elastic scattering (both cross sections carry the same Jacobian),
so `convertValue = false`, as for RATIO. The tree is `Ratio(dsdo_cm(a→a), rutherford(a→a))`, written with `predict`
or a keyword `rutherford` on an isDiff 1 or 4 line. It joins the χ² as every measurement does (T1). The Python
conversion in `pyazr/nds.py` then writes the ratio as measured instead of multiplying it out.

**Folded: a thick-target yield tied to a named target.** The measured Y(E_beam) (reactions per incident particle) is
∫ σ(E) / ε(E) dE over the target, with ε the stopping cross section of a named `Material`. It needs the stopping module
and an explicit fold, so it is a Phase 5 `Fold(kind = yield, target = name)` above a `σ` leaf, with the target named
in the measurement's `ApparatusRef` (a `target <name>` keyword, parsed in `io/`, resolved against a `<targets>`
section). Phase 4 reserves the field and the node name; it does not implement it. The name `yield` is already taken
by the SAMMY self-shielded neutron yield (§1.5): the folded kind is `target_yield` (decision M12).

**Model-only: a resonance strength ωγ.** No beam energy and no point: ωγ = (2J+1)/((2j₁+1)(2j₂+1)) (1+δ₁₂) Γ_a Γ_b / Γ
(eV, c.m.), from the observed partial widths of one level (`physics/level_parameters`). It enters as a
`ModelDatum{selector, value, error, provenance}` read from one `<modelData>` section in its owner through
`Config::SectionStream`; one `ModelDataTerm` in `src/Objective.cpp`, `ModelOnly() = true`, sums ((q_k(θ) − d_k)/δ_k)²
over the rows; the measured scattering lengths become rows of the same term with their section unchanged (decision
M10). A new model quantity is then a physics function and a selector, not a new term.

## 5. Decisions

**Taken 2026-11-09 (director's review):** M1 (a), M2 (c), M3 (a) with (b) carried to Phase 5 as an open decision for
the maintainer, M4 (a), M5 (a), M6 (b), M7 (b), M8 (b), M13 (b). **Taken after `data_kinds.md` (same day):** M9 — neither (a) nor (b) as written: the user-facing form is the named
`quantity <name> unit <u> frame <lab|cm> …` keyword block of `data_kinds.md` §5.1, one closed table of quantities
mapped to trees in `SegmentKinds`; a free `predict "<expr>"` form waits until a project needs a combination the table
and the SUM/RATIO blocks cannot write (after Phase 5). M10 (b), section name `<modelData>`. M12 as recommended; the
folded kind is `quantity yield … target <name>`. **Pending:** M11 until the thread-dependence diagnosis reports
(`fit_threads.md`).

**From the maintainer, 2026-11-09 (overrides M12 and the "old projects load unchanged" clauses of this note and of
`data_kinds.md` §5.1):** backward compatibility of the `.azr` is not required; a converter from the old format is.
The engine therefore reads one grammar only (format 3, `project_format3.md`): named quantities, named targets,
`<modelData>`. The integer isDiff codes, the positional segment columns and the `<targetInt>` index lists are known
to the converter (`io/src/Convert.cpp`, `AZURE2 --convert`, `pyazr.convert_azr`) and to nothing else; an old file
given to the engine is refused with a message naming the converter. The gate changes accordingly: the outputs of the
converted test projects are bit-identical to those of the originals on the previous engine. The format must read
naturally in the GUI and in Python: the same field names in the file, in the pyazr call and in the form, all taken
from one table of quantities exported by the engine.

**M1. Where `experiment/` lives.** (a) A standalone CMake target on `physics/` from the first step, with an engine
adapter; (b) inside `src/` first, split later. Recommendation (a): the one-way rule is then enforced by the linker
from the start, as for `io/` and `physics/`, and the nodes and frames are testable without an engine.

**M2. Where results live.** (a) In the points, as now; (b) per-thread scratch; (c) a `Results` owned by the caller,
slots assigned by the plan. Recommendation (c): it is the only option in which a value from an earlier evaluation
cannot be read without being written first (the poisoned run tests this), and it is what Phase 6's per-task scratch
becomes.

**M3. Lab data with a fitted shift.** (a) Keep today's comparison in c.m., y_c.m.(ν) and δy_c.m.(ν) recomputed in the
hand-over into a per-ν `MeasurementFrame`; (b) compare in the frame of the data: m_lab = m_c.m./J(ν), weights 1/δy_lab
fixed. Recommendation: (a) in Phase 4 (bit-identical); (b) proposed for Phase 5, since a χ² weight that moves with a
nuisance is a model choice nobody made; it moves fits with shifted lab-differential data and is classified there.

**M4. Folds in Phase 4.** (a) Inside the leaf (today's per-point fold), `Transform` as the only node that reads
sub-point values; (b) explicit `Fold` nodes now. Recommendation (a): the folds are entangled with mapping, probes,
automatic tolerances and windows, which Phase 5 rebuilds as kernels; moving them twice is work Phase 5 would undo.

**M5. The normalization.** (a) On the data side, f − n y, outside the tree; (b) a `Scale(n)` node, f/n − y. Recommendation
(a): (b) changes the arithmetic of every residual.

**M6. Component segments.** (a) Keep them as full `ESegment` copies; (b) components become leaves of their base
measurement; their calculation units belong to its plan, and the duplicated data go. Recommendation (b), with the term
order, the component keys and the `intEC.dat` layout (`ED` 3256–3470, `ECSignature`) kept; the layout is checked byte
for byte.

**M7. The gradient of the composition.** (a) Leave `AG` 697–1380 as the adjoint of SUM/RATIO/abundances/transforms;
(b) the nodes' `backward` replaces it in Phase 4. Recommendation (b): otherwise the composition stays written twice and
a new node would need a second hand-written adjoint. If a reordering moves gradient bits, the change is classified by
the archived fit runs and gradient = finite differences.

**M8. Amplitude filters.** (a) E1/E2 only, as now; (b) the general filter (multipoles, Jπ, l, s, internal/external) on
angle-integrated leaves in Phase 4, on dσ/dΩ in Phase 5 with the Legendre coefficients; (c) everything now.
Recommendation (b): the ROADMAP's T1 E1-only/E2-only test needs it, and on σ the filtered pieces add up exactly, so it
is testable by a sum rule; on dσ/dΩ the interference terms need the observable rework of Phase 5.

**M9. A text form for trees.** (a) Keywords only; (b) `predict "<expr>"` as the last step of Phase 4, without folds;
(c) after Phase 5. Recommendation (b): every new combination is otherwise a new keyword, i.e. new parsing code that an
assistant writes in a place of its choosing; with (b) it is a line in the project file.

**M10. Model-only data.** (a) A section and a term per quantity (as `<scatteringLengths>`); (b) one `<modelData>`
section and one term with a selector per quantity, the old section read into it. Recommendation (b), with the term's
position in the objective kept where T9 is now, so the χ² sum keeps its order.

**M11. History in working copies.** (a) Leave it to Phase 6, which removes the copies; (b) make a copy's state a function
of (θ, ν) now: compare re-placed grids exactly (H3), take the load-time angle and factor where a shift leaves the
kinematic range (H4), one component-key counter per `EData` (H10). Recommendation (b), coordinated with the agent
diagnosing the thread dependence (if its test shows H3 is the cause, this step is its fix). Calculations do not move
(one evaluation, no history); fits may, and each moved fit is classified "thread/history dependence removed" with the
1-thread and 8-thread results shown equal.

**M12. Names and codes.** The data and test isDiff tables disagree (3, 4, 5); `yield` means a neutron-sample yield.
Recommendation: keep both tables (old files), map both in one place (`SegmentKinds`), and give new kinds names in the
text form, never new isDiff codes; name the folded kind `target_yield`.

**M13. Who builds the trees.** (a) `io/` builds experiment objects (`ARCHITECTURE.md` §3); (b) `io/` parses into plain
specs and `experiment/` builds. Recommendation (b): `io/` stays dependent on nothing (Phase 3, check 9), and the
keyword → tree table is in one place next to the nodes.

## 6. Steps

Verification for every step unless stated otherwise: 179 × 3 identical (`/data0/skowrons/AZURE2-upgrade/tools/runall.sh`,
`classify.py`), `ctest -LE slow`, `tests/run_tests.sh` 133/133, calculations identical at 1 and 8 threads, and the
archived fit runs of `objective.md` §4 identical. Steps S2–S5 come first because they remove cross-object writes and
history.

| step | what | files | risk; bits | test that fails if it is wrong |
|---|---|---|---|---|
| S1 | history test (**done** 2026-11-10: `tests/upgrade/history_test.cpp`, `upgrade_history`): θ_A on a fresh clone vs on a pooled clone after θ_B, per-segment χ² bit for bit, on T2, `natural_neon`, `shift_sqrt_subpoints`, one THM project (share with, or take from, the thread-defect agent) | `AC` 704–745, `tests/upgrade/` | none (a test) | it is the oracle for S2 |
| S2 | copy state a function of (θ, ν) (**done** 2026-11-10 but H3, which `fit_threads.md` closed: the load-time frame, one key counter per `EData`, the β-delayed refusals): exact grid comparison (`ED` 1316–1321), load-time angle and factor beyond the kinematic range (`ES` 1143–1152, `EP` 1026–1031), one component-key counter per `EData` (`ED` 3930, 3970) | `ED`, `ES`, `EP` | calculations unchanged; fits may move (M11) | S1 passes; fits at 1 and 8 threads equal on T2 |
| S3 | `Results` and slots: the result members of `EPoint` (§1.1) go; `GenMatrixFunc`, `AMatrixFunc`, `THMMatrixFunc`, the folds write slots. In two commits: (a) slots behind the old accessors, (b) the buffer passed explicitly. **Done** ((a) 2026-11-09, (b) 2026-11-10: `EPoint` holds an offset, `EData` the layout, every owner its buffer; `include/Results.h`) | `EP`, `EPh`, `GM`, `src/AMatrixFunc.cpp`, `src/THMMatrixFunc.cpp`, `AG`, `api/src/AZUREAPI.cpp`, writer | large, mechanical; none | `AZURE_POISON_RESULTS=1` on T1–T3 and the 179 calculations: any slot read before written gives NaN |
| S4 | one model pass (**done** 2026-11-10: `include/ModelPass.h`, `src/ModelPass.cpp`; check 28): `Results::model` apart from leaf slots; the 19 loops become one `EvaluateModel`; the residual one function; `SetNaturalTerms`, the profiled norm and the segment χ² into `Results`; `TransformedValue` folds a local vector (W3, W4, W10 go) | `ES`, `AC`, `AG`, `src/Objective.cpp`, `src/EDataThm.cpp`, `src/AZURECalcMCMC.cpp`, API, writer | medium; none if the residual keeps (f − y n)/(δy n) | 179 × 3; `natural_target`, transmission and yield tests; MCMC tests |
| S5 | the plan (**done** 2026-11-10: `include/CalculationPlan.h`, `src/CalculationPlan.cpp`; tests `upgrade_plan`, `upgrade_plan_runtime`, `AZURE_PLAN_CHECK=1`; check 29): mapping as tasks and slot ranges; `local_mapped_points_` and `Clone`'s pointer re-mapping (`ED` 5640–5716) go; single writer per slot asserted (W1, W2, W7 become plan entries) | `ED` 4099–4250, 5598–5819, `EP` 2148–2348 | medium; none | threads 1 vs 8; ThreadSanitizer on T2 (`scripts/sanitize_upgrade.sh`); the assertion on all 179 projects |
| S6 | `experiment/` target and frames: the lab → c.m. functions moved out of `EPoint` with unit tests; the two conversion sequences and the two Jacobians become one | `experiment/`, `EP` 748–1120, `ES` 535–575, 1123–1170, CMake, check 8-style include check | low; none | closed forms: elastic θ_c.m., photon aberration, Jacobian = numerical dΩ_lab/dΩ_c.m. to 10⁻⁹; relativistic against `Kinematics.h` |
| S7 | `Measurement` and `MeasurementFrame`: input arrays shared, not cloned; the frame recomputed in the hand-over when a shift moves; `EPoint` input fields read from them | `experiment/`, `ES`, `ED` | medium; none | 179 × 3; T2 fit peak RSS at 8 threads recorded (expected to drop) |
| S8 | expression nodes with `value` and `backward` | `experiment/` | low; none (not wired) | per node: value against the old expression on recorded values bit for bit, `backward` against finite differences |
| S9 | `SegLine`/`ExtrapLine` to `io/` as `SegmentSpec`; `SegmentKinds`: the one keyword → tree table; `ESegment`'s isDiff flags derived from the tree; `--print-predictions` | `io/`, `experiment/`, `include/SegLine.h`, `ES` 25–250, `ED` 342–609, 688–970, `api/` | medium; none (set-up) | the printed tree of every segment line of `tests/` and `examples/` equals a committed table; `_azure2.parse_segment_line` returns the same fields |
| S10 | the tree pass replaces `CalculateTheoreticalCrossSection` and `TransformedValue` | `ES` 1605–1822, `AC` | the arithmetic order; none if S8's tests hold | 179 × 3; archived fits |
| S11 | the gradient through `backward` (M7); `AG` 697–1380 composition code deleted | `AG` | gradient bits may move | gradient vectors on T1–T3 bit for bit against the build before, else classified; `AZURE_GRAD_CHECK`; archived fits |
| S12 | components as leaves (M6): `componentSegments_` and `CreateComponentSegment` go | `ED` 3472–3977, 4218–4250 | `intEC.dat` layout, keys | `intEC.dat` byte for byte on the capture projects; 179 × 3 |
| S13 | general amplitude filters on σ (M8) | `GM` 128–258, `experiment/` | new; old outputs unchanged | T1: Σ over a complete set of filters = σ to 10⁻¹²; E1-only and E2-only fits of T1 recover the generator's widths; isDiff 5/6 unchanged |
| S14 | `predict` text form and the `rutherford` leaf (M9) | `io/`, `experiment/`, `GM`, `pyazr/datasets.py` `_KEYWORDS`, docs | new | ratio to Rutherford = 1 with all nuclear widths 0; round trip of a `predict` line; pyazr reads it through the binding |
| S15 | `<modelData>` and `ModelDataTerm`, ωγ first (M10); scattering lengths moved into it | `src/Objective.cpp`, `AG`, `physics/`, `api/` | χ² order kept; none for old projects | ωγ of an isolated level against Γ_aΓ_b/Γ by hand; the scattering-length tests unchanged |
| S16 | recipe in `CLAUDE.md`, checks 13–16 (§7) | `CLAUDE.md`, `tests/consistency/repo_consistency_test.py` | none | each check fails on an injected bypass |

Parallel branches: S6 and S8 can start at once, apart from the rest; S1–S2 is one branch; S3 → S4 → S5 is a chain;
S7 needs S3 and S6; S9 needs S8; S10 needs S4 and S9; S11, S12, S13, S14 need S10 and are independent of each other;
S15 needs only Phase 2 and can run any time; S16 last, though check 13 can land with S9.

## 7. The recipe for a new kind of data, and its checks

For `CLAUDE.md` "Recipes":

**A kind of data** (a measured quantity at a beam energy, as a ratio to Rutherford or an S factor).
1. If it is a combination of existing leaves and nodes, write it with `predict` in the project; no code. A keyword for
   it is optional sugar: one row in `experiment/src/SegmentKinds.cpp` and in `_KEYWORDS` (check 6).
2. A new observable: its formula in `physics/` with a test against an independent value (recipe "A physics function
   or an observable"); one `ObservableKind` value, its text-form name and its row in `SegmentKinds.cpp`; the
   calculation in the one `switch (kind)` of the leaf evaluation (compiled with `-Werror=switch`, so a missing case
   does not build).
3. A new node: one class in `experiment/` with `value` and `backward`, tested alone (value against a closed form,
   `backward` against finite differences).
4. Its data join the χ² as every measurement does; nothing is added to `src/Objective.cpp`.
5. A test project using it; `reference/data_formats.rst`.

**A model-only datum** (no beam energy: ωγ, a width, an ANC). Its formula in `physics/`, a selector in
`experiment/ModelDatum`, a row type in the `<modelData>` reader; `ModelDataTerm` evaluates every selector, so the
objective gets no new term (check 11).

Checks (`tests/consistency/repo_consistency_test.py`, and one run-time check):
- 13. Kinds are decided in one place: comparisons of isDiff codes, `IsAdvanced`, `GetLegacyOperationType` and the
  prediction keywords appear only in `SegmentKinds.cpp` and `io/`. An `if (seg->IsMyKind())` in `GenMatrixFunc` or
  `ESegment` fails it; after S9 `ESegment` and `EPoint` have no kind flags left to test.
- 14. No result is stored in a set-up object: `SetFitCrossSection`, `SetNaturalTerms`, `SetSegmentChiSquared` do not
  exist; `Measurement` is used only through `const Measurement&` outside its builder.
- 15. Every `ObservableKind` value and every node class has a unit test named after it in `experiment/tests/` and a
  row in the kinds table of `reference/data_formats.rst`.
- 16. Run time: `AZURE_POISON_RESULTS=1` fills `Results` with NaN at the start of each evaluation; ctest runs T2
  with it. A value kept from an earlier evaluation, wherever it was hidden, makes the χ² NaN.

## 8. Not established

- Whether H3 (or H4) causes the fit's thread dependence: inferred from the code, not run. S1 decides it.
- Whether every node can reproduce today's arithmetic bit for bit: designed so (§2.7), not yet shown; S8's tests show it
  before anything is wired.
- The THM experiment groups (`src/EDataThm.cpp`) were read only where they write norms; their fit and profile
  arithmetic is assumed to be unaffected by moving the stored values.
- The consumers of point fields in pyazr (beyond the binding's arrays) and in `EData::WriteOutputFiles` were not listed
  one by one; S3 (a) keeps the accessors so they keep working while they are moved.
- The kinds surveyed from EXFOR (`data_kinds.md`, another agent) were not available; §4 covers the three cases asked
  for.
