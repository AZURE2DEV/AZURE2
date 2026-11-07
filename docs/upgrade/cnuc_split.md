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
| `formal_from_brune` | System, b, Ẽ_λ, γ̃, B_c, channel functions, sign rule, eigensolver | formal E_λ (MeV), γ (MeV^½) | `CNuc::TransformIn` (step 1) |
| `brune_from_formal` | System, b, formal E, γ, B_c, channel functions, tolerance (MeV), eigensolver | Ẽ, γ̃, iterations, status, last step (MeV), last ΔB_c | `CNuc::TransformOut` (step 1) |
| observed widths (step 2) | System, b, E, internal γ, γ^ext, basis, channel functions | Γ (eV), C (fm^−½), μ, Q, B per (λ, c) | `parameters.out`, `GetTransformParams`, `CheckRadiativeWidths` |
| `park_from_brune`, `brune_from_park` | System, b, E, γ, channel functions | γ in the other basis, J | `ConvertAmplitudeBasis` (step 2), Park norm (step 3) |
| `channel_functions` (step 3) | partition α, l, E_c (MeV, c.m.) arrays | P, S, S′ (MeV^−1), φ, σ_l | the level-parameter functions, collision matrix |
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
`CNuc::CalcExternalWidth`. Nothing else. The adapter carries two options that exist only to reproduce the engine's
history (section 5): the order of E_c and the quantities to compute. They change no meaning elsewhere and go with the
adapter in step 3.

## 3. Order of the moves

1. **Level-parameter conversions in the engine** (this branch): `CalcBoundaryConditions`, `TransformIn` (all bases),
   the formal → Brune iteration of `TransformOut`. They call `level_parameters`, which was shown bit-identical for the
   Brune and Park amplitudes and B_c, so the engine's numbers can stay identical; the formal basis keeps the engine's
   GSL eigensolver (section 4).
1b. **Decision 10 on** (formal basis): `BoundaryRule::LowestLevel` and `FormalSignRule::LargestPositive` in the
   engine's calls, `level_parameters`' own eigensolvers, the GSL adapters deleted; formal-basis goldens re-baselined
   in their own commit (observables unchanged to ~10⁻¹¹, formal parameters change by B_c and sign).
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
1.0 × 10⁻¹⁴, and the overall sign of a formal level can differ from GSL's (`level_parameters.md` §7). Both functions
take the eigensolver as an option (`GeneralizedEigenSolver`, `FormalToBruneOptions::solver`); the engine passes GSL's
(`GslGeneralizedEigen`, `GslSymmetricEigen` in `src/CNucLevelParameters.cpp`), so the formal basis runs through
`level_parameters` with the engine's arithmetic and its outputs stay bit-identical. With the own solver they do not
(section 6). Decision 10 is taken: B_c from the lowest level and the largest |γ| of each formal level positive, both
on when the engine runs through `level_parameters`, with a re-baseline of the formal-basis goldens. That is the step
after this one (step 1b): the sign rule makes the solver's sign irrelevant, so the GSL option and its adapter go in
the same commit.

## 5. Channel energy and the Coulomb memo

The engine forms the channel energy in three orders: E − e − S (`TransformIn`), E − (S + e) (`CalcBoundaryConditions`),
E − S − e (`TransformOut`), with e = e₁ + e₂ the excitation energies of the pair and S its separation energy. They agree
whenever e = 0 and may differ by one ulp otherwise. The Coulomb memo (`CoulFuncCache`) returns a stored value for any
energy within 10⁻¹² MeV of a stored one, so which of two energies one ulp apart is evaluated first decides the value
both receive: the engine's numbers depend on the order of its Coulomb-function calls at the 10⁻¹³ level. The adapter
therefore takes the order of E_c and the quantities to compute from its caller (`ChannelEnergyOrder`,
`ChannelQuantities`) and evaluates each on its own `ChannelFunc` in the engine's call sequence. A first version that
computed S, dS/dE and P at every call, in one order, moved the Brune amplitudes of `tests/channel_space` by 10⁻¹³,
flipped the GSL sign of 6 of its formal levels, and moved two `secondary_gamma` outputs in the last printed digit. Both
the three orders and the memo's tolerance go when `channel_functions` replaces the adapter (step 3), with a re-baseline
at that level if an output moves.

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

### Remaining steps

Steps 2–7 of section 3. Next: the observed widths from total amplitudes (step 2), which removes the five copies of the
γ → Γ conversion and lets `ConvertAmplitudeBasis` call `park_from_brune`.

## 7. Step 1b: design

Step 1b makes the formal basis follow one rule and the adapter compute the channel functions in one way. It changes
no observable beyond rounding; it changes the formal parameters of a block (E_λ in MeV, γ_λc in MeV^½) by the choice of
B_c and by the overall sign of a level. What each part becomes:

**Boundary values** (`boundary_values`). B_c = S_c(Ẽ_low), Ẽ_low the lowest Brune energy of the block's R-matrix
levels (MeV, E_x). The rule is not an option: `BoundaryRule` and `Convention::boundary` are deleted. With it the formal
parameters of a block no longer depend on the order in which the file lists its levels: a permutation of the listing
permutes the formal levels and changes nothing else (to rounding). The lowest level's formal energy equals its Brune
energy, since Ẽ_λ is an eigenvalue of the formal level matrix with B_c = S_c(Ẽ_λ) (Brune 2002).

**Sign of a formal level** (`formal_from_brune`). The amplitude of largest |γ_λc| of each formal level is positive; of
two equal magnitudes the first channel of the block decides. The eigenproblem leaves the sign free and observables are
even in it, so the rule replaces the eigensolver's choice. `FormalSignRule` and `Convention::formal_sign` are deleted.
The Brune amplitudes that `brune_from_formal` returns carry the formal level's sign: each iteration's eigenvector is
dominated by the level itself, whose component `symmetric_eigen` makes positive.

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
orders agree to 10⁻¹³; the lowest level's formal energy equals its Brune energy; the largest |γ| of each formal level is
positive; widths too large for the radii are refused. `physics_level_parameters_from_cnuc`: `level_parameters` against
the engine without the level-sign allowance (same solver, same rule), so a sign or B_c that differs fails it.
`transform_level_order` (b): the formal parameters of `asc` and `desc` are the same levels. Engine before (4bf7aa06)
and after on all 179 `.azr` under `tests/` in three bases, every changed output file classified: (a) formal parameters
(B_c, sign), observables ≤ 10⁻¹⁰ relative; (b) last digits from the channel-energy order and the call sequence; (c)
anything else stops the step.
