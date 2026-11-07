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
commit. The `Convention` options of HANDOFF decision 10 (B_c from the lowest level, largest |γ| positive) stay off.

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

## 3. Order of the moves

1. **Level-parameter conversions in the engine** (this branch): `CalcBoundaryConditions`, `TransformIn` (all bases),
   the formal → Brune iteration of `TransformOut`. They call `level_parameters`, which was shown bit-identical for the
   Brune and Park amplitudes and B_c, so the engine's numbers can stay identical; the formal basis keeps the engine's
   GSL eigensolver (section 4).
2. **Observed widths and the basis switch.** One physics function for Γ_λc (and C_λc, μ, Q) from the total amplitude
   γ^int + γ^ext; `TransformOut`'s second loop, `GetTransformParams`, `PrintTransformParams`, `CheckRadiativeWidths`
   read it. `ConvertAmplitudeBasis` calls `park_from_brune`/`brune_from_park` on the parameter vector.
3. **Channel functions** (after the pycoulomb prerequisites, ROADMAP Phase 1): `CalcShiftFunctions` becomes tables
   of `channel_functions(α, l, E)` per fit stage; the process-wide memo goes; the adapters' `ChannelFunc` and `WhitFunc`
   calls are replaced. Here the three orders of the channel energy (section 5) are unified, with a re-baseline at the
   10⁻¹⁶ level if any output moves.
4. **External capture**: `CalcExternalWidth` and the bound-state normalization N_f^½ (`SqrtNFFactor`, now in five
   places) into `physics/external_capture`.
5. **Coupling coefficients**: `SortPathways`, `CalcAngularDists`, `CalcCaptureAnalyzingPower` into `observables`,
   computed once per reaction with `AngularCoupling` (3j/6j/9j against `AngCoeff`: differences of one ulp are possible,
   so this step is checked to 10⁻¹⁴ and re-baselined by its own commit if an output moves).
6. **Data accessors**: callers of `GetPair`/`GetJGroup` read the System (partition α, block b) through the adapter
   map; the pathway trees stop living on `PPair`.
7. Phase 2 and 3 take (e) and (c). `CNuc` is then empty and deleted.

## 4. The formal basis and the eigensolver

`formal_from_brune` and `brune_from_formal` solve N a = E M a and the symmetric problem of the iteration. With
`SymmetricEigen.h`'s solver the formal energies agree with GSL's to 6.6 × 10⁻¹⁶ relative and the amplitudes to
1.0 × 10⁻¹⁴, and the overall sign of a formal level differs from GSL's on 3 of 19 levels of T1–T3
(`level_parameters.md` §7). Both functions take the eigensolver as an option (`GeneralizedEigenSolver`,
`FormalToBruneOptions::solver`); the engine passes GSL's (`GslGeneralizedEigen`, `GslSymmetricEigen` in
`src/CNucLevelParameters.cpp`), so the formal basis runs through `level_parameters` with the engine's arithmetic. The
switch to the own solver goes together with decision 10 (the sign rule makes the solver's sign irrelevant) and a
re-baseline of formal-basis outputs; section 6 records what it would change.

## 5. Channel energy

The engine forms the channel energy in three orders: E − e − S (`TransformIn`), E − (S + e) (`CalcBoundaryConditions`),
E − S − e (`TransformOut`), with e = e₁ + e₂ the excitation energies of the pair and S its separation energy. They
agree whenever e = 0 and may differ by one ulp otherwise. The adapter `ChannelFunctionsFromCNuc` uses E − e − S. Section
6 states whether any output moved.

## 6. Record

(Filled in as the moves are made.)
