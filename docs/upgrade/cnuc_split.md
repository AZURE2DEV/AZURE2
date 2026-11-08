# Splitting `CNuc`: plan and record

`CNuc` (`include/CNuc.h`, `src/CNuc.cpp`, 3628 lines at 93f8cbc5) holds the compound-nucleus model of a project: the
particle pairs (partitions α), the J^π groups with their channels c = (α, l, s) and levels λ, and the level parameters
E_λ (MeV, E_x of the compound nucleus) and γ_λc (MeV^½). It also computes from them: the boundary values B_c, the
conversions between observed and formal parameters, the shift functions S_c(E_λ), the reaction pathways and angular
coefficients, and the external-capture amplitudes; and it parses the `<levels>`, `<abundances>` and external-capture
input, fills the Minuit parameter vector, and writes `parameters.out`. ROADMAP Phase 1 ("Split CNuc") moves the
nucleus data into the `System` (`physics_system.md`) and the computations into functions of the System in `physics/`;
the parsing goes to `io/` (Phase 3), the parameter vector and penalties to `inference/` (Phase 2). This page classifies
every member function, gives its callers and what it reads and writes, fixes the order of the moves, and records what
has moved and the evidence that the engine's outputs did not change.

## 1. Rule for each move

A move replaces the body of a `CNuc` member by a call into `physics/` (through the adapters `src/CNucSystem.cpp`,
`src/CNucLevelParameters.cpp`) and deletes the old body in the same commit. It is accepted when the engine's output
files are bit-identical to those of the engine before the move: every project of `tests/` and the test projects T1–T3
run with the binary before and after, in the Brune basis (default), the formal basis (`--no-brune`) and Park's
(`--use-park`), and every file of `output/` compared byte for byte; plus `tests/upgrade/check.sh`, `ctest -LE slow` and
`tests/run_tests.sh`. A move that cannot meet this (a different eigensolver, a different order of operations) is
either made with the engine's arithmetic injected (section 4) or left for a deliberate re-baseline with its own
commit. Step 1 kept the `Convention` options of HANDOFF decision 10 off; step 1b turned the formal basis to its rules
and deliberately changed outputs, each change classified (sections 6, 7).

## 2. Classification

Classes: **(a)** nucleus data → `System`; **(b)** a physics function of the System → `physics/`; **(c)** input parsing
→ `io/`; **(d)** output and reporting; **(e)** other (parameter vector and penalties → `inference/`; lifetime and
plumbing). "Callers" counts call sites outside `CNuc.cpp` by file; "int." are callers inside `CNuc.cpp`. Line numbers
are those of 93f8cbc5.

### 2.1 (a) Nucleus data

| member | callers | reads / writes | target |
|---|---|---|---|
| `NumNuclei`, `GetNucleus` (2), `NucleusLabel`, `FindOrAddNucleus` (private) | EData, ESegment, ParameterLabel, AZUREAPI, CNucSystem; int. Fill, PrintNuc, PrintTransformParams, AddPair | `nuclei_` (Z, A); FindOrAdd appends | `System::nuclei()` |
| `NucleusOfPair` | none | `GetPair(n)->GetNucleus()` | delete (no caller) |
| `NumPairs`, `GetPair`, `IsPairKey`, `GetPairNumFromKey`, `GetPairFromKey`, `IsPair`, `AddPair` | ~250 sites: AMatrixFunc, GenMatrixFunc, RMatrixFunc, EData, EPoint, ESegment, AZUREGrad, AZUREAPI, GUI, … | `pairs_`; Add/IsPair only from Fill | `System::partitions()` with a key → α map (the file's pair key is an io/ label) |
| `NumJGroups`, `GetJGroup`, `IsJGroup`, `AddJGroup` | ~180 sites, same files | `jgroups_` (groups own channels, levels, pathways) | `System::blocks()`; levels as `JpiBlock` arrays |
| `GetMaxLValue`, `SetMaxLValue` | int. Fill, ParseExternalCapture | `maxLValue_` (never initialized before Fill) | `System` (largest l) |
| `CapturePairKeys`, `RatePairsProblem` | ReactionRate, AZUREMain, GUI MCMC thread | photon partitions of a nucleus | `System` queries (capture partitions of a nucleus) |
| `IsTransformedIn`, `TransformError` | EData (3), CNucSystem; int. Initialize | `transformedIn_`, `transformError_` | disappear: level_parameters returns values, not a state |
| `ThmCoherentValues`, `SetThmCoherentValues` | THMMatrixFunc, EDataThm | `thmCoherentValues_` (THM coherent backgrounds, fit parameters) | (e) `inference/` parameters handed to the THM observable |
| `Clone` | AZURECalc, AZURECalcMCMC, EData, ReactionRate, ParameterLimitsManager (13) | deep copy | (e) disappears: the System is immutable and shared; per-evaluation parameters are arrays |

### 2.2 (b) Physics functions of the System

| member | callers | reads | writes | target |
|---|---|---|---|---|
| `CalcBoundaryConditions` (1975–2003) | int. Initialize | E of a group's first level, pair S, e; `ChannelFunc::Shift` | `AChannel::SetBoundaryCondition` | `boundary_values` — **moved (step 1)** |
| `TransformIn` (1200–1610) | AZURECalcMCMC (3), ParameterLimitsManager, AZUREAPI (2), GUI MCMCTab; int. Initialize | input E, Γ/ANC/γ, flags; ChannelFunc, WhitFunc; B_c; `CalcExternalWidth` | `SetGamma`, `SetE` (formal), `SetParkNorm`, `SetSqrtNFFactor`; warnings | `amplitudes_from_observed`, `subtract_external_capture`, `formal_from_brune` — **moved (step 1)** except `SqrtNFFactor` (external capture, step 4) |
| `TransformOut` (2540–2820) | AZUREMain, AZURECalc, AZURECalcMCMC (2), AZURE2, EData, ReactionRate, ParameterLimitsManager, AZUREAPI (5), GUI (3) | fit E, γ; B_c; ChannelFunc, WhitFunc; `CalcExternalWidth(false)` | `SetTransformE/Gamma/Iterations`, `SetSqrtNFFactor`, `SetExternalGamma`, `SetBigGamma` | formal → Brune iteration: `brune_from_formal` — **moved (step 1)**; observed widths with external capture: step 2 |
| `ConvertAmplitudeBasis` (3422–3470) | AZUREParams::ReconcileBasis | Minuit values, `ShiftDerivative` | Minuit values and errors | `park_from_brune`, `brune_from_park` (step 2) |
| `CalcShiftFunctions` (3237–3410), `InvalidateShiftFunctions` | ~40 sites (AZURECalc, AZURECalcMCMC, AZUREAPI, GUI, …); int. TransformOut | fit E, γ; Coulomb functions; process-wide memo `g_shiftDerivMemo` | level S_c, S'_c, S''_c, Park J, `SqrtNFFactor` | `channel_functions(α, l, E)` (needs pycoulomb) + Park norm from `park_from_brune`; the memo becomes the per-stage table (RULES §2.3) (step 3) |
| `ParkNormPenalty`, `kParkNormScale` | AZURECalc (3), AZURECalcMCMC (2), AZUREAPI (3), AZUREGrad | level Park J | — | (e) a term of the `Objective` (Phase 2) |
| `CalcExternalWidth` (3513–3595) | CNucLevelParameters; int. TransformIn, TransformOut | final state E, γ, N_f, conversion factor; `ECIntegral` | — | `physics/external_capture` (step 4) |
| `SortPathways` (1616–1859) | int. Initialize | pairs, groups, channels, EC final states | Decay → KGroup → MGroup/ECMGroup trees on `PPair` | coupling coefficients per reaction in `observables` (step 5) |
| `CalcAngularDists` (2051–2199) | int. Initialize | pathway trees; `AngCoeff::ClebGord`, `Racah` | KLGroup, Interference (Z̄₁Z̄₂), diagonal pathways | `observables` coefficients from `AngularCoupling` (step 5) |
| `CalcCaptureAnalyzingPower` (2248–2379) | AMatrixFunc::PointTBar, GenMatrixFunc (lazy, mutex) | pathway trees; Racah, 9j | `Decay::AddCaptureAyTerm` | `observables` (step 5); the lazy build under a mutex becomes setup state |

### 2.3 (c) Input parsing → `io/Project` (Phase 3)

| member | callers | reads | writes |
|---|---|---|---|
| `Fill` (188–427) | AZUREMain, AZURE2, AZUREAPI (2), GUI (6), tests (6) | `<levels>` via `NucLine`, radius overrides, `includeAllJPi` | pairs, groups, channels, levels; calls the five below |
| `AddHardSphereGroups` (606–653) | int. Fill | particle pairs | level-less groups (a System rule: every J^π the pairs allow up to l_max) |
| `CheckIdenticalChannels` (669–695) | int. Fill | identical pairs with spin | refusal |
| `ParseExternalCapture` (728–858), static `ExitPairKeysInFile` (707–726) | int. Fill | EC final states, segment exit keys (re-reads the file) | `SetECParams`, N_f integrals, conversion factors; 784–793 and 811–855 unreachable |
| `ReadAbundances` (438–514, private) | int. Fill | `<abundances>` | `samples_` |
| `CheckMassSeparationConsistency` (881–943) | int. Fill | masses, separation energies | warnings (a check of the System at load) |

### 2.4 (d) Output and reporting

| member | callers | writes |
|---|---|---|
| `PrintTransformParams` (2948–3036), `PrintRelativisticParams` (3048–3110) | AZUREMain, AZURECalcMCMC, AZURECalc | `parameters.out`, `relativistic.out` |
| `GetTransformParams` (3114–3174) | AZURE2, AZURECalcMCMC, ParameterLimitsManager, AZUREAPI (5), GUI (3) | the observed parameters as a vector |
| `CheckRadiativeWidths` (2841–2942) | AZUREMain (2) | warning Γ_γ/Γ_p > 0.1 |
| `PrintNuc`, `PrintPathways`, `PrintBoundaryConditions`, `PrintAngularDists` | AZUREMain, AZUREAPI; int. Initialize | `checks/*.chk` |
| `PrintCompoundFromParams` | none | `std::cout`; delete |

The γ → Γ conversion (penetrability, Whittaker factor, static moments, (E_γ/ħc)^(2L+1)) is written five times: in
`TransformIn`, `TransformOut`, `CheckRadiativeWidths`, `GetTransformParams`, `PrintTransformParams`. After step 2 the
reports read one function's output.

### 2.5 (e) Parameter vector, penalties, lifetime

| member | callers | target |
|---|---|---|
| `FillMnParams` (2439–2506) | 30 sites (AZURECalc, AZURECalcMCMC, AZUREAPI, ParameterLimitsManager, GUI, tests) | `inference/ParameterVector` |
| `FillCompoundFromParams` (2512–2534), `FillCompoundFromParamsPhysical` (1136–1158) | ~40 sites | disappear: level_parameters takes θ |
| abundance members: `HasAbundanceBlock`, `NumSamples`, `GetSample`, `FindSample`, `DefaultSample`, `NumAbundanceParams`, `AbundanceOrdinal`, `AbundanceOf`, `AbundanceParamName`, `GetAbundance`, `AbundanceSum`, `SampleSumSigma`, `AbundancePenalty`, `AbundanceSumPenalty` | AZURECalc, AZUREGrad, AZURECalcMCMC, EData, ESegment, ParameterLabel, AZUREAPI, GUI | samples are `experiment/` (target composition); penalties are `Objective` terms |
| `ScatteringLengthConstraints`, `HasScatteringLengthConstraints` | AZURECalc (5), AZUREMain, AZURECalcMCMC, AZUREAPI | `Objective` term; constraints read by `io/` |
| `Initialize` (945–1021) | AZUREMain, AZURE2, AZUREAPI (2), GUI (4), tests (3) | the setup stage: System + level_parameters + coefficients; disappears with the class |

`DoubleFactorial` is declared in `CNuc.h` but defined in `src/DoubleFactorial.cpp` and not used by `CNuc`; its
declaration moves with its users (`ECIntegral`, `EffectiveCharge`).

### 2.6 GUI callers

The Qt GUI is removed from the upgrade (HANDOFF decision 3). Its calls into `CNuc` (`gui/src/AZURESetup.cpp`,
`AZUREMCMCThread.cpp`, `MCMCTab.cpp`, `ThmDiagnostics.cpp`; `tests/gui/`) go with it and constrain no move. No member is
used only by `gui/`, so none is classified "GUI, to be deleted"; `RatePairsProblem` and `Fill` keep their CLI and API
callers.

### 2.7 Interfaces of the target parts

Each part below is self-contained: it receives the System (or arrays from it) and plain arrays in physical units, and
returns arrays; it reads no other part's state. `CNuc` reaches them only through the two adapters, whose reach into
`CNuc` is fixed here.

| part (module) | inputs | outputs | callers |
|---|---|---|---|
| `System` (`physics/`) | partitions α (Z, mass in u, 2j, π, g, e₁, e₂ and S in MeV), channels (α, kind, l, 2s, a_c in fm), blocks (2J, π, nucleus, channel list) | immutable value | everything below |
| `boundary_values` | System, block b, E_λ (MeV, E_x, listed order), rule, channel functions | B_c (dimensionless) per block channel | `CNuc::CalcBoundaryConditions` (step 1) |
| `amplitudes_from_observed` | System, b, columns, Γ (eV) / C (fm^−½) / γ̃ (MeV^½) with flags, basis, channel functions | γ̃ or γ^P (MeV^½), Park J, diagnostics D/2, terms, notes | `CNuc::TransformIn` (step 1) |
| `subtract_external_capture` | System, b, total photon γ (MeV^½), γ^ext (complex, MeV^½) | internal γ, entries where \|Im γ^ext\| > γ | `CNuc::TransformIn` (step 1) |
| `formal_from_brune` | System, b, Ẽ_λ, γ̃, B_c, channel functions | formal E_λ (MeV), γ (MeV^½), each level oriented like its Brune level | `CNuc::TransformIn` (step 1) |
| `brune_from_formal` | System, b, formal E, γ, B_c, channel functions, tolerance (MeV) | Ẽ, γ̃, iterations, status, last step (MeV), last ΔB_c | `CNuc::TransformOut` (step 1) |
| observed widths (step 2) | System, b, E, internal γ, γ^ext, basis, channel functions | Γ (eV), C (fm^−½), μ, Q, B per (λ, c) | `parameters.out`, `GetTransformParams`, `CheckRadiativeWidths` |
| `park_from_brune`, `brune_from_park` | System, b, E, γ, channel functions | γ in the other basis, J | `ConvertAmplitudeBasis` (step 2), Park norm (step 3) |
| `channel_functions` (step 3, open channels; 3b all E_c and the table) | channel kinematics, l, a_c (fm), E_c (MeV, c.m.) | P, S, S′ (MeV^−1), e^{−iφ}, σ_l | `CoulFunc` (step 3); the level-parameter functions, collision matrix (3b) |
| `external_capture` (step 4) | System, final state (E, γ, N_f), photon channel, E grid | γ^ext (complex, MeV^½), N_f^½ | level parameters, collision matrix |
| coupling coefficients (step 5) | System, reaction (entrance α, exit α′), L_max | Z̄₁Z̄₂ and analyzing-power coefficients per pathway pair | observables |
| `io/Project` (Phase 3) | `.azr` text | System, project levels, measurements, priors | front ends |
| `inference/` terms (Phase 2) | parameter vector θ | χ² terms (abundance, Park wall, scattering lengths) | the Objective |

Channel functions are an input of every level-parameter function until `channel_functions` exists; the engine supplies
them through `ChannelFunctionsFromCNuc`.

**What the adapters read from `CNuc`.** `BuildSystemFromCNuc` reads the public accessors of pairs, J-groups, channels
and levels (spins, masses, energies, radii, flags; level E and γ only for `CNucSystemContent::Parameters`).
`ChannelFunctionsFromCNuc` reads a channel's pair and l and evaluates `ChannelFunc`/`WhitFunc`;
`ProjectLevelsFromCNuc` reads the `<levels>` values of a filled `CNuc`; `ExternalWidthsFromCNuc` calls
`CNuc::CalcExternalWidth`. Nothing else. Since step 1b the adapter has no options: one channel energy, and every call
returns S, S′ and P (sections 5, 7).

## 3. Order of the moves

1. **Level-parameter conversions in the engine** (this branch): `CalcBoundaryConditions`, `TransformIn` (all bases),
   the formal → Brune iteration of `TransformOut`. They call `level_parameters`, which was shown bit-identical for the
   Brune and Park amplitudes and B_c, so the engine's numbers can stay identical; the formal basis keeps the engine's
   GSL eigensolver (section 4).
1b. **Decision 10 on** (formal basis; done, section 6): B_c from the lowest level and one sign rule as the only rules
   of `level_parameters` (the `Convention` options deleted), its own eigensolvers (the GSL adapters deleted), one
   channel energy and one adapter call for S, S′, P. The sign rule is the orientation of the Brune level, not "largest
   |γ| positive", which changes external-capture observables (section 6).
2. **Observed widths and the basis switch.** One physics function for Γ_λc (and C_λc, μ, Q) from the total amplitude
   γ^int + γ^ext; `TransformOut`'s second loop, `GetTransformParams`, `PrintTransformParams`, `CheckRadiativeWidths`
   read it. `ConvertAmplitudeBasis` calls `park_from_brune`/`brune_from_park` on the parameter vector.
3. **Channel functions** (design `coulomb.md`; done, section 6): the Coulomb functions from `physics/coulomb`
   (pycoulomb's algorithms, `external/coul` deleted); P, S, φ, S′ and η, ρ of an open channel as
   `channel_functions` in physics/, which `CoulFunc` and `Kinematics` call. The process-wide memo (`CoulFuncCache`,
   10⁻¹² MeV tolerance) stays: deleting it was built and measured (3–17× slower, a test past its time limit) and moved
   to step 3b.
3b. **Channel functions, rest** (`coulomb.md` §1.4, §2.3, §6): (i) the Whittaker functions in `physics/coulomb`
   without GSL (proposed: COULCC at imaginary argument, against mpmath and GSL), then `ShftFunc`, `WhitFunc` and
   `ChannelFunc`'s bound branch replaced by `channel_functions` for all E_c, and the adapter
   `ChannelFunctionsFromCNuc` calling it; (ii) the per-stage table of P, S, S′, φ, σ_l per (α, l, a_c, E_n) in place
   of the Coulomb rows of `EPoint::EDependentTables` and of `g_shiftDerivMemo`, and a per-integral table of F, G at
   the quadrature radii for the external-capture integrals — this gives back the speed the memo bought (section 6,
   step 3); then (ii′) `CoulFuncCache` deleted: the branch `feature/upgrade-1-coulomb-nomemo` (commit c3dc22e5 and
   its docs, on top of step 3) is re-cut on the table, its output comparison repeated and its timings required near
   those of step 3; (iii) the hybrid-potential channels as (α, l, j) entries; (iv) `CoulombTable` (an opt-in interpolation
   of the same functions) and `--gsl-coul`: keep, remap to a `coulomb::Method`, or retire (user-visible; to decide).
   The non-smooth last digits of S and S′ (section 5) remain until (i)–(ii): `coulsm` is good to 10⁻¹⁰ only.
4. **External capture**: `CalcExternalWidth` and the bound-state normalization N_f^½ (`SqrtNFFactor`, now in five
   places) into `physics/external_capture`.
5. **Coupling coefficients**: `SortPathways`, `CalcAngularDists`, `CalcCaptureAnalyzingPower` into `observables`,
   computed once per reaction with `AngularCoupling` (3j/6j/9j against `AngCoeff`: differences of one ulp are possible,
   so this step is checked to 10⁻¹⁴ and re-baselined by its own commit if an output moves).
6. **Data accessors**: callers of `GetPair`/`GetJGroup` read the System (partition α, block b) through the adapter
   map; the pathway trees stop living on `PPair`.
7. Phase 2 and 3 take (e) and (c). `CNuc` is then empty and deleted.

## 4. The formal basis and the eigensolver

`formal_from_brune` and `brune_from_formal` solve N a = E M a and the symmetric problem of the iteration with
`SymmetricEigen.h`'s solvers. In step 1 they took the eigensolver as an option and the engine passed GSL's, so that the
formal basis stayed bit-identical; with the own solvers the formal energies agree with GSL's to 6.6 × 10⁻¹⁶ relative
and the amplitudes to 1.0 × 10⁻¹⁴, but the overall sign of a formal level could differ. Step 1b removed the option: the
sign of each formal level is fixed by a rule, so the solver's choice no longer reaches an output. Decision 10 proposed
"largest |γ| positive"; the comparison of step 1b showed that this rule changes external-capture observables in the
formal basis (section 6), and the rule taken is the orientation of the Brune level (`level_parameters.md` §4).

## 5. Channel energy and the Coulomb memo

Before step 1b the engine formed the channel energy in three orders: E − e − S (`TransformIn`), E − (S + e)
(`CalcBoundaryConditions`), E − S − e (`TransformOut`), with e = e₁ + e₂ the excitation energies of the pair and S its
separation energy. They agree whenever e = 0 and may differ by one ulp otherwise. Since step 1b there is one, E_c =
(E_x − S) − e, the order of the collision matrix at a data point (section 7). Two properties of the engine's channel
functions make the last digits of its numbers depend on such details. (1) The Coulomb memo (`CoulFuncCache`) returns a
stored value for any energy within 10⁻¹² MeV of a stored one, so the order of the Coulomb-function calls decides
which of two nearby energies sets the value. (2) S from `coulsm` is accurate to about 10⁻¹⁰ and S′, a Richardson
difference of S with step 10⁻³ E_c, to about 10⁻⁸; neither is smooth below that level, so moving E_c by one ulp
re-draws those digits. Step 1 kept every caller's order and call sequence and stayed bit-identical; step 1b gives up
both, and the outputs move at the level these two properties allow (section 6). Both go when `channel_functions`
replaces the adapter (step 3).

## 6. Record

### Step 1 (branch `feature/upgrade-1-cnuc-split`, 2026-11-07)

Moved, each in its own commit, old code deleted:

- `CalcBoundaryConditions` → `boundary_values` (rule `FirstListedLevel`, order E − (S + e)). A group whose first
  listed level is not an R-matrix level (an external-capture final state listed first) used to keep B_c uninitialized;
  it now takes the first R-matrix level, or 0. No project in `tests/` has one.
- `TransformIn` → `amplitudes_from_observed`, `subtract_external_capture`, `formal_from_brune` (GSL's solver), in all
  three bases. `CNuc` keeps the wording of its warnings (from the diagnostics the functions now return: D/2, each
  column's term, column notes, the entries whose |Im γ^ext| exceeds the total), Park's J on the level, and the
  bound-state normalization N_f^½ (step 4). Groups are still made formal one after the other (external capture in
  the engine's order, `level_parameters.md` §6). Two levels of a group at one energy (`--no-brune`) are now refused
  with a message instead of giving NaN.
- `TransformOut`, the formal → Brune iteration → `brune_from_formal` (GSL's solver, order E − S − e). The warning for a
  level that does not converge now reports the last energy step; the old code printed the difference after the update,
  always 0. The observed-width loop stays (step 2).

`CNuc.cpp` went from 3628 to 3311 lines.

**Evidence.** The engine before (93f8cbc5) and after (this branch), same compiler (GCC 11.5, Release), were run on
every `.azr` file under `tests/` (179 projects, including T1, T1_start, T2, T2_true, T3, T3_start), menu 1, one
thread, in each basis, and every file written was compared byte for byte (progress-bar lines of the log excluded):

| basis | runs | output files | identical | differing |
|---|---|---|---|---|
| Brune (default) | 179 | 510 (80 `parameters.out`) | 510 | 0 |
| formal (`--no-brune`) | 179 | 483 (75 `parameters.out`) | 483 | 0 |
| Park (`--use-park`) | 179 | 510 (80 `parameters.out`) | 510 | 0 |

The logs (`run.log`) are identical except the three of `tests/screening_6li` (`--no-brune`), whose non-converging
levels now report the steps 118.45, 40.9187 and 85.4476 MeV instead of 0. The same comparison with
`level_parameters`' own eigensolver in place of GSL's (formal basis only) changes 53 output files: 20 `param.par` (698
formal amplitudes change sign), 25 `parameters.out` (last printed digit of widths), 7 observable files (up to
2.2 × 10⁻¹¹ relative, T1) and `covariance.out` of T2; that switch therefore waits for decision 10.

GUI off (decision 3), GCC 11.5, Condor: no compiler warnings; `tests/upgrade/check.sh` passes (goldens unchanged);
`ctest -LE slow` 106/107, the one failure `adaptive_grid_reference` (known, decision 8); `tests/run_tests.sh` 132/132. T2 benchmark, both engines on one host (d32cepyc138), wall time after/before: calc 1.17/1.19 s (1 thread),
0.42/0.41 s (8 threads); fit200 79.32/80.10 s, 36.16/36.28 s; peak RSS within 1 % (`tests/upgrade/bench/T2.tsv`).

### Step 1b (branch `feature/upgrade-1-cnuc-split-1b`, 2026-11-07)

Done, each in its own commit after the design (section 7): B_c at the lowest level and the formal sign as rules of
`level_parameters`, the `Convention` options `BoundaryRule`, `FormalSignRule` and the eigensolver options with the
GSL adapters deleted; one channel energy (E_x − S) − e and one adapter call returning S, S′, P, the options
`ChannelEnergyOrder`, `ChannelQuantities` deleted; the sign rule changed from "largest |γ| positive" to the orientation
of the Brune level after the comparison below; a strongly mixed test block. In the formal basis the engine prints a
note for each J^π group whose first listed level is not its lowest (formal input from earlier builds does not apply).

**Evidence.** The engine before (4bf7aa06) and after (7f4af9d7; the later commits change tests and documents only),
GCC 11.5, Release, GUI off, were run on every `.azr` under `tests/` (179 projects), menu 1, one thread, in each basis,
and every file written was compared number by number (`classify.py`; each nested test folder in a private copy: the
step-1 harness let `X` and `X/models` share a work tree, which lost one run of `shift_subpoint_angle` here).
Printed numbers carry 11 significant digits in `AZUREOut`, 6 decimals in `parameters.out`, so a printed difference
of 10⁻¹¹–10⁻¹⁰ can come from a smaller change; the full-precision values of five projects were compared through pyazr.
Classes: (a) formal basis, B_c or sign; (b) last digits, from the order of E_c and the call sequence (section 5);
(c) anything else.

| basis | files | identical | changed, by class (largest relative change) |
|---|---|---|---|
| Brune | 689 | 670 | (b) 8 observable files, printed ≤ 9.8 × 10⁻¹¹, full precision ≤ 3.6 × 10⁻¹¹ (`channel_space`), 2.0 × 10⁻¹¹ (`7Li_p_ay`); 5 `parameters.out` (last printed digit of the 0⁺ background level of the ⁷Li + p models, full precision 6.3 × 10⁻⁹); 6 `run.log` (numbers in warnings, ≤ 1.7 × 10⁻¹¹) |
| Park | 689 | 685 | (b) 3 observable files, printed ≤ 6.7 × 10⁻¹¹ (`secondary_gamma` twolevel; full precision 3.9 × 10⁻¹¹ with the amplitudes bit-identical: the call sequence through the memo); 1 `run.log` (1.1 × 10⁻¹²) |
| formal | 662 | 595 | (a) 21 `param.par`: formal levels reversed (GSL's sign → orientation) with magnitudes unchanged to 1.3 × 10⁻¹⁵ of the level's largest amplitude (exact zeros become ≤ 2 × 10⁻¹⁷), new B_c in `natural_neon` (5 groups) and `be7p_fit` (1), an energy 0 → −8.8 × 10⁻³³ MeV in `screening_6li`; 25 `parameters.out`: reported signs, B_c (iterations, energies), the formal → Brune iteration's last digits (≤ 5.8 × 10⁻⁹; 4.5 × 10⁻⁷ for the non-converging levels of `screening_6li`), and in `beam_profile_kernel` one photon width with external capture, 0.000149 → 0.000053 meV, because `TransformOut` now gives its final state the orientation of its input; 5 observable files, ≤ 8.5 × 10⁻¹¹ (`natural_neon`, `T1`, `T2_true`); `covariance.out` of T2, 5.0 × 10⁻¹³; 10 `run.log` (the notes, 16 groups in 7 projects, and warning numbers). (b) 3 observable files of `channel_space`, ≤ 9.1 × 10⁻¹¹ (as in the Brune basis) |

No file is in class (c) with the rule as built. With decision 10's "largest |γ| positive" one was: on
`tests/beam_profile_kernel` the rule reversed the 6.917 MeV 2⁺ final state (its largest formal amplitude is the E2
photon one), and the formal-basis ¹²C(α,γ) capture to it, extrapolated over 111 energies, changed by up to 4.51
relative (E_c.m. = 4.50 MeV), the Brune basis unchanged; the menu-1 runs of `tests/` contain no such capture output, so
the classification above did not see it. With the orientation rule the five captures of that project are
identical to those of the engine before in both bases (5 × 111 printed values). That project also shows that the formal and Brune bases disagree in capture with external capture,
before and after this step: by 1.5 × 10⁻² to 9.1 relative over its five final states (`level_parameters.md` §6).

The Brune and Park outputs move by more than the 10⁻¹² expected. A variant with the old channel-energy order and the
new call sequence leaves the Brune and Park amplitudes of the five probed projects bit-identical and their observables
within 5 × 10⁻¹⁴, except the Park `secondary_gamma` case (3.9 × 10⁻¹¹); the order of E_c therefore causes the amplitude
changes, through the non-smooth last digits of S and S′ (section 5), and the 0⁺ level of the ⁷Li + p models, whose
Brune denominator is negative (the engine warns that its widths are out of reach), turns them into 6.3 × 10⁻⁹.

Goldens: the T1, T3 (`golden.py --write`) and T2 calculation goldens (`regen_compare.py --regenerate`), regenerated
with both engines, are byte-identical; the T2 fit differs by ≤ 3.0 × 10⁻⁷ (`param.sav`), less than the engine before
differs from its own committed fit (6.3 × 10⁻⁷, a run-to-run spread of the 8-thread MIGRAD end point). No golden and no
`tests/*/expected/` file is re-baselined. `transform_level_order` (b) states the new rule (commit of decision 10).

Checks (GCC 11.5, GUI off, Condor): no compiler warnings; `tests/upgrade/check.sh` PASS; `ctest -LE slow` 107/107;
`tests/run_tests.sh` 131/132, the one failure `shared_point_tables` (d), the peak-RSS ratio without margin
(HANDOFF, step 1): fit at 4 threads 316 MB > 3 × the calculation's; in three repeats of each engine the step-1b
engine failed once (308 MB > 3 × 67 MB) and the engine before passed 3/3, with the fit at 301–308 MB (step 1b) and
310–312 MB (before): the calculation's peak RSS varies from 67 to 105 MB between runs, not the fit's. T2 benchmark, both engines alternately on one host (d32cepyc147), mean of two:
calc 1.16 / 1.21 s (1 thread), 0.39 / 0.39 s (8); fit200 78.54 / 78.53 s, 34.30 / 34.36 s; peak RSS within 1 %
(`tests/upgrade/bench/T2.tsv`).

### Step 3 (branch `feature/upgrade-1-coulomb`, 2026-11-07)

Done, each in its own commit after the design (`coulomb.md`, alone): `physics/coulomb` (coulsm, coul, coulfg, coulcc
from pycoulomb azure2-prereqs 7adef98 and coulcc-cpp 37f0f9b, tests against mpmath and GSL); the engine's five
direct calls of `cwfcomp` replaced by `CoulombWavesAt` → `coulomb(l, η, ρ, coulsm)` and `external/coul` deleted
(4053 lines); `channel_functions` and the channel kinematics in physics/, called by `CoulFunc` and `Kinematics`.
Built and moved to `feature/upgrade-1-coulomb-nomemo` (step 3b): `CoulFuncCache` deleted (457 lines, two reference
tests, every `SetUseGlobalCache`, the `TransientScope`).

**Evidence.** The engine before (f549fa23) and after steps 3 (cd1dc9ad), 3-kinematics (the `channel_functions`
routing) and, on `-nomemo`, 3-memo, GCC 11.5, Release, GUI off, every `.azr` under `tests/` (179), menu 1, one thread, three bases,
every file compared number by number (`classify.py`, step-1b harness):

| comparison | Brune | formal | Park | changed |
|---|---|---|---|---|
| before → routing | 682/689 identical | 652/662 | 682/689 | hybrid-potential projects only (`T3`, `T3_start`, `hybrid_potential`): observables ≤ 9.9 × 10⁻¹¹; one formal-basis width 5000.000255 → 5000.000006 keV (input 5000 keV; the formal → Brune iteration); 1 log 1.9 × 10⁻¹² |
| routing → `channel_functions` | 689/689 | 662/662 | 689/689 | none |
| `channel_functions` → memo deleted (`-nomemo` only) | 486/689 | 458/662 | 486/689 | observables ≤ 7.0 × 10⁻¹¹ (`energy_shift_sqrt`, `relativistic` he5, T1), except an exact zero's rounding residue (`ec_e1m1_phase` s_1, E1 at 0° and 180°, 1.5 → 3.0 × 10⁻²⁵ against 1.2 × 10⁻⁹ at 60°); 20 `intEC.dat` ≤ 4.0 × 10⁻¹¹ of the line's largest value; χ² of `energy_shift_sqrt` and `shift_fn` in the 12th digit; every log without the "Cleared Coulomb function cache" line |

Why the routing moved the hybrid channels: `cwfcomp` returns complex F, G whose imaginary parts are not exactly 0
for real arguments (a probe: 7 of 125 points, up to 1.0 × 10⁻¹⁵ relative); the hybrid model built H⁺ = G + iF from
the complex values, `coulomb()` returns real parts. On `-nomemo` the memo's deletion moved the others: it returned the
functions of an energy within 10⁻¹² MeV of the one asked for. Every change is ≤ 10⁻¹⁰ relative on observables or the explained
residue; no golden and no `tests/*/expected/` file is re-baselined. `intEC.dat` changes at 4 × 10⁻¹¹, below its
use; its signature is not bumped (an old file stays valid to that level).

Checks and speed: HANDOFF, 2026-11-07, `feature/upgrade-1-coulomb`; speed in `coulomb.md` §7.

### Step 3b (ii), (ii′) (branch `feature/upgrade-1-channel-table`, 2026-11-07)

Done, design first (`channel_table.md`, alone): the per-stage channel table (P, S, e^{−iφ}, e^{i(σ_l − σ_0)} per
exact (α, l, E_c), built by `EData` before the point loops, read-only, released after setup), a thread-private
exact-argument scratch of the radial Coulomb waves of external capture, and `CoulFuncCache` deleted (the `-nomemo`
commits re-cut). Outputs against the memo base: as on `-nomemo` (≤ 7.0 × 10⁻¹¹ on observables, one exact zero's
residue, `intEC.dat` ≤ 4.0 × 10⁻¹¹), and byte-identical to `-nomemo` in all 689 Brune-basis files. The THM byte checks
pass 5/5; 1 vs 8 threads byte-identical (T1, T2, T3, 13N). Then, on the director's answers: free-shift refills read
the segment's tables of its last three shift states (interim until the analytic shift derivative), and the capture
amplitude reads √P and the entrance phase at a_c from the table. Speed against the memo base on one host: T2 fit200
−13 % (1 thread), −23 % (8), calc equal or faster, natural_neon −13 %, `tgt_fixed` equal, 13N +1.5 % (+1 to +9 % over
five pairs, inside the base's spread); peak RSS of a fit +13 % (1 thread), +30 % (8) (`channel_table.md` §6). Not done: (i) the Whittaker functions, (iii) hybrid (α, l, j) rows, (iv) `CoulombTable`/`--gsl-coul`; the
per-slot rows of `EDependentTables` and `g_shiftDerivMemo` stay (`channel_table.md` §2).

### Remaining steps

Steps 2, 3b and 4–7 of section 3. Next: the observed widths from total amplitudes (step 2), which removes the five copies of the
γ → Γ conversion and lets `ConvertAmplitudeBasis` call `park_from_brune`. Step 2 or 4 should also settle how the formal
basis evaluates external capture (the final state's formal or Brune amplitude, and the order of `TransformIn`), which
the comparison above shows differs from the Brune basis.

## 7. Step 1b: design

Step 1b makes the formal basis follow one rule and the adapter compute the channel functions in one way. It changes
no observable beyond rounding; it changes the formal parameters of a block (E_λ in MeV, γ_λc in MeV^½) by the choice of
B_c and by the overall sign of a level. What each part becomes:

**Boundary values** (`boundary_values`). B_c = S_c(Ẽ_low), Ẽ_low the lowest Brune energy of the block's R-matrix
levels (MeV, E_x). The rule is not an option: `BoundaryRule` and `Convention::boundary` are deleted. With it the formal
parameters of a block no longer depend on the order in which the file lists its levels: a permutation of the listing
permutes the formal levels and changes nothing else (to rounding). The lowest level's formal energy equals its Brune
energy, since Ẽ_λ is an eigenvalue of the formal level matrix with B_c = S_c(Ẽ_λ) (Brune 2002).

**Sign of a formal level** (`formal_from_brune`). *As designed:* the amplitude of largest |γ_λc| of each formal level
positive (decision 10). *As built, after the comparison (section 6):* each formal level keeps the orientation of the
Brune level it comes from, a_λ,k(λ) > 0 for the eigenvector of rank k(λ), and `brune_from_formal` orients the
eigenvectors of its iteration the same way. The design's premise, that observables are even in a level's sign, holds
within the R-matrix but not for an external-capture final state f: the photon amplitudes of the levels that decay to f
and the external-capture amplitude ∝ γ_f carry f's sign together, and "largest |γ| positive" reversed one such state
of `tests/beam_profile_kernel` and changed its capture by a factor of up to 4.5. Both rules are independent of the
eigensolver. `FormalSignRule` and `Convention::formal_sign` are deleted.

**Eigensolvers.** `formal_from_brune` and `brune_from_formal` use `SymmetricEigen.h` only. The solver options
(`GeneralizedEigenSolver`, `SymmetricEigenSolver`, `FormalToBruneOptions::solver`) and the adapter's
`GslGeneralizedEigen`, `GslSymmetricEigen` are deleted.

**Channel energy.** One formula, E_c = (E_x − S_α) − e (MeV, c.m.), with S_α the separation energy of partition α and
e = e₁ + e₂ the excitation energies of the pair: E_x − S_α is the energy above the partition's ground-state
threshold, and e is then the energy the pair holds in excitation. It is the order of the collision matrix at a data
point (`EPoint.cpp`, `inEnergy - GetSepE() - GetExE()`) and of `CNuc::CalcShiftFunctions` and `TransformOut`'s widths,
so a level energy and a data point at the same E_x give the same E_c to the bit. `level_parameters`' own E_c (which
branch, above or below threshold, a column takes) uses the same formula. `ChannelEnergyOrder` is deleted.

**Channel functions** (`ChannelFunctionsFromCNuc`). Every call returns S_c, S′_c and P_c of a particle channel,
evaluated on one `ChannelFunc` in the order S, P, S′ (P from `WhitFunc` below threshold). `ChannelQuantities` is deleted.
The extra evaluations (S where only P and S′ are read, at setup; S′ where only S is read, in the formal iteration)
cost nothing measurable in a Brune fit, which calls the adapter only at setup; the T2 benchmark checks it. Because the
Coulomb memo (`CoulFuncCache`, energies matched to 10⁻¹² MeV) now sees another sequence of calls, Brune and Park
outputs may move at the 10⁻¹³ level (section 5); that is accepted and reported.

**Compatibility.** Formal parameters given as input — `--no-transform` with `--no-brune`, or a parameter file
(`param.sav`, `param.par`) written by an earlier build in the formal basis — were defined with B_c at the first listed
level. For a block whose first listed level is not its lowest they now describe another R-matrix. At load, in the
formal basis, the engine prints one note per such block with both energies. No project in `tests/` reads formal
input for such a block except `transform_level_order`, which writes its own parameter file with the same build.

**Tests.** `physics_level_parameters_unit`: B_c at the lowest level; the formal parameters of the block listed in two
orders agree to 10⁻¹³; the lowest level's formal energy equals its Brune energy; the sign rule (as built: each formal
level points the way of its Brune level, also in a strongly mixed block where the eigensolver's sign does not); widths
too large for the radii are refused. `physics_level_parameters_from_cnuc`: `level_parameters` against
the engine without the level-sign allowance (same solver, same rule), so a sign or B_c that differs fails it.
`transform_level_order` (b): the formal parameters of `asc` and `desc` are the same levels. Engine before (4bf7aa06)
and after on all 179 `.azr` under `tests/` in three bases, every changed output file classified: (a) formal parameters
(B_c, sign), observables ≤ 10⁻¹⁰ relative; (b) last digits from the channel-energy order and the call sequence; (c)
anything else stops the step.

## 8. Step 2: design

Step 2 gives the observed quantities of a level one definition, computed from its total amplitudes, and settles which
amplitude of an external-capture final state the formal basis uses. It changes no observable in the Brune and Park
bases beyond rounding; in the formal basis it changes the external-capture observables of projects whose final state
shares its J^π group with other levels (below).

**The quantities.** For a level λ of block b with Brune (or Park) energy Ẽ_λ (MeV, E_x) and, in each channel c, the
total amplitude γ^tot_λc = γ^int_λc + γ^ext_λc (MeV^½; γ^ext complex, non-zero only for a photon channel with external
capture, from `CNuc::CalcExternalWidth` until step 4), with N_λ = 1 + Σ_{c particle} (γ^int_λc)² S′_c(Ẽ_λ) (Brune) or
N_λ = 1 (Park):

    W_λc = s_λc · 2 P_c(Ẽ_λ) |γ^tot_λc|² / N_λ,      s_λc = sign Re γ^tot_λc                                       (8.1)

with P_c the width factor of `ChannelValues` (`level_parameters.md` §1). What W is, and the value reported in the
`<levels>` form (the form `amplitudes_from_observed` reads back):

| channel | W | reported value |
|---|---|---|
| particle, E_c > 0 | Γ_λc (MeV) | Γ (eV) = 10⁶ W, signed |
| particle, E_c < 0 | C²_λc (fm⁻¹), P = P^b | C (fm^−½) = s |W|^½ |
| particle, E_c = 0 or P = 0 | 0 | 0 |
| particle with the amplitude flag (gammaIsRWA) | as above | γ̃ (MeV^½): γ^int (Brune), γ^int / J^½ (Park) |
| photon | Γ_γ (MeV), P = (E_γ/ħc)^{2L+1} (1 with `--use-rmc`) | Γ_γ (eV), signed |
| photon to the level's own state, M1 / E2 (static moment) | μ² (μ_N²) / (Q · 100 (αħc)^½)², sign reversed for half-integer J | μ (μ_N) / Q (b), signed |
| beta (Fermi, Gamow–Teller) | B = Re γ^tot (Brune), Re γ^tot / J^½ (Park) | B |

This is Eq. (2.1)–(2.2) of `level_parameters.md` read backwards, with the external part added inside the square.
Function: `observed_levels(System, b, Ẽ, γ^int, γ^ext, is_amplitude, basis, cf)` in `level_parameters` → per (λ, c)
the kind of the row, W, the reported value; per level N_λ. `observed_from_amplitudes` (the inverse used by the tests)
becomes the same function with γ^ext = 0. Validity as in §2: N_λ > 0 (Brune), J_λ > 0 (Park).

**The five copies and where each goes.**

| copy (`CNuc.cpp` at 9d50667f) | what it computes | after step 2 |
|---|---|---|
| `TransformIn` | Γ, C, μ, Q → γ̃ | already `amplitudes_from_observed` (step 1); unchanged |
| `TransformOut`, second loop (2423–2516) | W per channel (`SetBigGamma`), γ^ext (`SetExternalGamma`) | `CNuc::ObservedWidths` calls `observed_levels`; stores W and the reported value |
| `GetTransformParams` (2796–2881) | the reported values from W (√, signs, μ, Q, 10⁶) | reads the stored reported values |
| `PrintTransformParams` (2630–2711) | the same, printed | reads W and the stored values; prints as today (C unsigned) |
| `CheckRadiativeWidths` (2523–2624) | Γ_p and Γ_γ from the fit amplitudes (formal amplitudes read as Brune; no N in Γ_γ; no external part) | the Brune values of the start parameters (as `TransformOut`, without its messages), then `ObservedWidths`; the ratio from the stored W |

`ConvertAmplitudeBasis` (3100–3150) computes Park's J a sixth time; it calls `park_from_brune` / `brune_from_park` on
each block's (E, γ) of the parameter vector and scales the errors by the J^½ these return. `TransformOut` is split into
the Brune values (`brune_from_formal`, with the messages) and `ObservedWidths`. Levels outside the R-matrix (an
external-capture final state in a group of its own) have no block of the R-matrix System; the adapter gets a content
`CNucSystemContent::AllGroups` that lists every group and level, used only for the observed quantities.

**External capture in the formal basis.** The capture amplitude to a final state f is proportional to f's asymptotic
normalization, C_f ∝ γ_fc N_f^{−½} with N_f the bound-state normalization (`SqrtNFFactor`), because outside a_c the
bound state is C_f W_l(2κr). That relation holds for the amplitude of f's own pole: Brune's γ̃_fc (B_c = S_c(Ẽ_f)) is
the residue of the full R-matrix at Ẽ_f, so C_f follows from γ̃_f in every basis. A formal amplitude γ_fc (constant
B_c = S_c(Ẽ_low)) equals γ̃_fc only when f is the only level of its block (then the eigenproblem is 1 × 1); in a block of
several levels it is a mixture of the block's poles, and using it for C_f describes a different bound state. The
engine today takes γ_f as the cross section's final-state amplitude in each basis: the formal one in the formal basis
(`AMatrixFunc.cpp` 560, 1193; `RMatrixFunc.cpp` 392; N_f from the same amplitude), and in `TransformIn` γ^ext of a later
group sees the formal amplitude of a final state in an earlier one (`level_parameters.md` §6), while `TransformOut`
uses the Brune one. On `tests/beam_profile_kernel` the 0⁺ ground state and 6.049 MeV state share a group of four
levels, and the 2⁺ 6.917 MeV state one of several: the two bases disagree by 1.5 × 10⁻² to 9.1 in those captures.

Decision: **the final state enters external capture with its Brune (observed) amplitude in every basis**, which is what
the Brune and Park bases already do. In the formal basis this needs Brune's γ̃_f at the current formal parameters of
f's block, i.e. `brune_from_formal` of that block at each evaluation; its own amplitude γ̃_f and N_f computed from it
replace γ_f and its N_f in the three cross-section sites and in `CalcExternalWidth`. It is stored on the final state
(`ALevel::ECGamma`, set where the fit parameters are set, `FillCompoundFromParams`; equal to the fit amplitude in the
Brune and Park bases and for a final state alone in its block, so those stay bit-identical), from a read-only setup
record of the blocks that hold a final state (built in `Initialize`, shared by the clones). `TransformIn` computes γ^ext
of every group before any group is made formal. The formal parameters of a block (E, γ, B) are unchanged; only what the
capture term reads from the final state changes.

**Tests.** (1) `physics_level_parameters_unit`: `observed_levels` against hand values for each row of the table (a
level with a bound channel, a static M1 and E2 moment at half-integer J, a beta channel, an external part with
Re and Im), and observed → Brune → observed and observed → Park → observed with γ^ext = 0 (as before). A wrong N, a
dropped |·|² of the complex amplitude or a wrong sign rule fails it at order 1. (2) Brune vs formal: on
`tests/beam_profile_kernel` the five captures (111 energies each) and the `parameters.out` photon widths agree between
`--no-brune` and the default to ≤ 10⁻⁹ relative after the change (they differ by up to 9.1 before) — the independent
check, since the two bases are different parametrizations of one R-matrix; a ctest runs it (two engine runs, 41
energies). (3) Engine before (9d50667f) and after on all 179 `.azr` × three bases: Brune and Park ≤ 10⁻⁹ or explained;
formal: (a) external capture to a final state that shares its group, expected to move, each listed; (b) rounding;
(c) anything else stops the step. `run.log` of `CheckRadiativeWidths` changes where a flagged level's ratio now
includes N_λ and γ^ext (listed).

**Not in step 2.** N_f^½ itself and γ^ext stay engine code (step 4); the formal-basis evaluation pays one
`brune_from_formal` per final-state block per evaluation (measured on T1–T3, which have no such block, and on
`beam_profile_kernel`).
