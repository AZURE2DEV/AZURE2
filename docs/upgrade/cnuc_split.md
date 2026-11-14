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
plumbing). The table is the inventory of `src/CNuc.cpp` (2907 lines) and `include/CNuc.h` (306 lines) at a16f2a69,
the start of `feature/upgrade-cnuc-split-3`, with the members defined in `src/CNucExternalCapture.cpp` (step 4).
"Lines" is the definition's length (1 for a member defined in the header). "Callers" counts call sites outside
`src/CNuc*.cpp` (ext.) and inside (int.), by a name search over `src/`, `api/src/`, `include/`, `tests/` (script
`work/cnuc3/inventory.py`): names that other classes also use (`GetPair`, `GetNucleus`, `NumPairs`, `Fill`,
`Initialize`, `Clone`) are over-counted (≈). "Destination" is where the member goes and the step of section 3 that
moves it; "done" names the branch of this note's record (section 6) where the move is recorded.

### 2.1 Inventory at a16f2a69

| member | lines | class | callers ext. / int. (main files) | destination |
|---|---|---|---|---|
| `NumPairs`, `GetPair` | 3, 4 | (a) | ≈26 / 21, ≈310 / 110 (every observable, EData, EPoint, API, tests) | `System::partitions()` (step 6b) |
| `IsPairKey`, `GetPairNumFromKey`, `GetPairFromKey` | 9, 12, 4 | (a) | 26 / 0, 114 / 1, 7 / 3 | key → α map of the System (the key is an `io/` label) |
| `IsPair`, `AddPair` | 12, 4 | (a)/(c) | 0 / 3, 0 / 2 (Fill, ParseExternalCapture) | built from the records (step (c)) |
| `NumJGroups`, `GetJGroup`, `IsJGroup`, `AddJGroup` | 3, 4, 16, 3 | (a) | 63 / 28, ≈235 / 136, 1 / 6, 0 / 3 | `System::blocks()`, levels as `JpiBlock` arrays (step 6b) |
| `NumNuclei`, `GetNucleus` (2), `NucleusLabel`, `FindOrAddNucleus` | 3, 3+3, 5, 9 | (a) | 6 / 5, ≈51 / 24, 8 / 5, 0 / 1 | `System::nuclei()` |
| `NucleusOfPair` | 3 | (a) | 0 / 0 | delete (no caller) |
| `GetMaxLValue`, `SetMaxLValue` | 3, 3 | (a) | 0 / 1, 0 / 2 | `System` (largest l) |
| `CapturePairKeys`, `RatePairsProblem` | 9, 15 | (a) | 2 / 2 (ReactionRate), 1 (AZUREMain) | System queries |
| `IsTransformedIn`, `TransformError` | 1, 1 | (a) | 3 / 2 (EData), 0 / 1 | disappear: level_parameters returns values, not a state |
| `CalcBoundaryConditions` | 35 | (b) | 0 / 1 (Initialize) | glue of `boundary_values` (done, step 1); goes with the class |
| `TransformIn` | 206 | (b)+(d) | 6 / 1 (API, MCMC, ParameterLimitsManager) | glue of `amplitudes_from_observed`, `subtract_external_capture`, `formal_from_brune` (done, step 1); 95 lines are its warnings (d) |
| `TransformOut`, `TransformToBrune`, `ObservedWidths` | 4, 104, 42 | (b)+(d) | 15 / 1, 0 / 2, 0 / 2 | glue of `brune_from_formal`, `observed_levels` (done, steps 1, 2); 35 lines of `TransformToBrune` are its warnings (d) |
| `ConvertAmplitudeBasis` | 46 | (b)/(e) | 2 (AZUREParams) | glue of `park_from_brune`, `brune_from_park` on the parameter vector (done, step 2); goes with the parameter vector (step 7) |
| `CalcShiftFunctions`, `InvalidateShiftFunctions` | 167, 9 | (b) | 34 / 1, 3 / 1 | `channel_functions` + Park's J; the memo `g_shiftDerivMemo` → the per-stage table (step 3b (ii) remainder, `channel_table.md` §2) |
| `ParkNormPenalty`, `kParkNormScale` | 14 | (e) | 4 (Objective, MCMC) | a term of the `Objective` (it already asks it) |
| `SortPathways` | 85 | (b) | 0 / 1 | fills the trees from `reaction_pathways` (done, step 6a); disappears with the trees (step 6b) |
| `CalcAngularDists` | 144 | (b) | 0 / 1 | reads the trees, coefficients from `ObservableCoupling` (done, step 5); a consumer of step 6b |
| `CalcCaptureAnalyzingPower` | 78 | (b) | 2 (AMatrixFunc, GenMatrixFunc; lazy, mutex) | the same; its lazy build becomes setup state (step 6b) |
| `SetUpECFinalBlocks`, `ECFinalStatesFromFormal`, `SetUpFormalCaptureMasks`, `FormalCaptureMasksFromParams`, `FormalCaptureMasksOutsideAdjoint`, `MaskedLevelCoupled`, `FormalCaptureMask`, `HasFormalECFinalBlocks` | 52, 25, 37, 34, 12, 5, 5, 1 (`CNucExternalCapture.cpp`) | (b) | 0–1 each (AMatrixFunc, AZUREGrad) | glue of `external_capture` and `brune_from_formal` per evaluation (done, step 4); level_parameters of θ (step 7) |
| `Fill` | 225 | (c) | ≈18 (AZUREMain, AZURE2, API, tests) | builds itself from plain `<levels>` records read by `io/` (step (c) of this branch) |
| `ReadAbundances` | 77 | (c) | 0 / 1 | the `<abundances>` grammar into `io/`; the samples into `experiment/` (with the parameter vector, step 7) |
| static `ExitPairKeysInFile` | 27 | (c) | 0 / 1 | a second reader of segment lines: the segment reader's |
| `ParseExternalCapture` | 131 | (c)/(b) | 0 / 1 | final-state set-up (N_f integral, K_c); lines 784–797 and 815–859 unreachable: delete |
| `AddHardSphereGroups` | 48 | (c)/(a) | 0 / 1 | a System rule (every J^π the pairs allow up to l_max) |
| `CheckIdenticalChannels`, `CheckMassSeparationConsistency` | 27, 63 | (d) | 0 / 1 each (Fill) | checks file, const access (step (d)) |
| `Initialize` | 80 | (e) | ≈24 | the setup stage; 37 of its lines are two warnings (d); disappears with the class |
| `PrintNuc`, `PrintPathways`, `PrintBoundaryConditions`, `PrintAngularDists` | 86, 104, 37, 49 | (d) | 2 / 0, 0 / 1, 0 / 1, 0 / 1 | report file, const access (step (d)) |
| `PrintTransformParams`, `PrintRelativisticParams` | 81, 63 | (d) | 3 (AZUREMain, AZURECalc, MCMC), 0 / 1 | report file (step (d)) |
| `CheckRadiativeWidths` | 93 | (d) | 2 (AZUREMain) | 75 lines of message → checks file; the observed widths at a vector stay (`ObservedWidthsAt`) |
| `PrintCompoundFromParams` | 14 | (d) | 0 / 0 | delete (no caller) |
| `GetTransformParams` | 19 | (e) | 10 (API, AZURE2, MCMC, EData, ParameterLimitsManager, ReactionRate) | the observed parameter vector: with the parameter vector (step 7) |
| `FillMnParams`, `FillCompoundFromParams`, `FillCompoundFromParamsPhysical` | 71, 23, 20 | (e) | 43, 36 / 3, 8 | `inference/ParameterVector`; level_parameters takes θ (step 7; another agent changes the vector now) |
| abundance members (`HasAbundanceBlock`, `NumSamples`, `GetSample` (2), `FindSample`, `DefaultSample`, `NumAbundanceParams`, `AbundanceOrdinal`, `AbundanceOf`, `AbundanceParamName`, `GetAbundance`, `AbundanceSum`, `SampleSumSigma`, `AbundancePenalty`, `AbundanceSumPenalty`) | 1–13 each, 78 in all | (e) | 1–10 each (Objective, AZUREGrad, ESegment, API, ParameterLabel, MCMC) | samples → `experiment/` (target composition); penalties → `Objective` terms (step 7) |
| `ScatteringLengthConstraints`, `HasScatteringLengthConstraints` | 1, 1 | (e) | 9, 8 | the constraints read by `io/`, the term in `Objective` (step 7) |
| `ThmCoherentValues`, `SetThmCoherentValues` | 1, 1 | (e) | 2, 2 (THM) | parameters handed to the THM observable (step 7) |
| `Clone` | 4 | (e) | ≈30 (AZURECalc, MCMC, ReactionRate) | disappears: the System is immutable and shared (Phase 6) |
| the shift-derivative memo (`ShiftDerivKey`, `g_shiftDerivMemo`) | 35 | (b) | CalcShiftFunctions, InvalidateShiftFunctions | the per-stage table (step 3b (ii) remainder) |

`DoubleFactorial` is declared in `CNuc.h` but defined in `src/DoubleFactorial.cpp` and not used by `CNuc`; its
declaration moves with its users (`ECIntegral`, `EffectiveCharge`).

By class, of the 2907 lines (definitions only): (a) 144, (b) 955 (of which ≈ 130 are the warnings of `TransformIn` and
`TransformToBrune`), (c) 508, (d) 617, (e) 290 (of which 80 are `Initialize`); the other 393 are includes and the
comments and blank lines between definitions.

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
2. **Observed widths and the basis switch** (done, section 6). One physics function for Γ_λc (and C_λc, μ, Q) from the total amplitude
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
5. **Coupling coefficients** (done, section 6; design section 9; `SortPathways` moved to step 6 with the trees it
   builds): `SortPathways`, `CalcAngularDists`, `CalcCaptureAnalyzingPower` into `observables`,
   computed once per reaction with `AngularCoupling` (3j/6j/9j against `AngCoeff`: differences of one ulp are possible,
   so this step is checked to 10⁻¹⁴ and re-baselined by its own commit if an output moves).
6. **Data accessors** (section 10): callers of `GetPair`/`GetJGroup` read the System (partition α, block b) through the adapter
   map; the pathway trees stop living on `PPair`. In two parts: 6a the enumeration (`SortPathways`) as a physics
   function of the System, which the engine's trees are filled from; 6b the consumers of the trees, one by one
   (design: section 11).
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

### Step 2 (branch `feature/upgrade-1-cnuc-split-2`, 2026-11-08)

Done after the design (section 8, alone), each in its own commit: `observed_levels` in `level_parameters` and the
engine reading it (`TransformOut` split into `TransformToBrune` and `ObservedWidths`; `GetTransformParams`,
`PrintTransformParams`, `CheckRadiativeWidths` read the stored W, value and kind; `ConvertAmplitudeBasis` through
`park_from_brune`/`brune_from_park`); the formal-basis external capture (final state's Brune amplitude, N_f from it,
Brune energy; γ^ext subtracted before any group is made formal); then the final states restricted to those above the
lowest level of their block. Deleted: the γ → Γ code of `TransformOut`'s second loop, `GetTransformParams`,
`PrintTransformParams` and `CheckRadiativeWidths`, and the J of `ConvertAmplitudeBasis` (CNuc.cpp 3306 → 3349 lines
with the two new members; the formulas now live once, in `LevelParameters.cpp`).

As built, beyond the design: (1) the final state's **energy** in the capture integrals and the photon energy
(`EPoint`, formerly `GetE()`, the formal energy in the formal basis) is its Brune energy (`ALevel::ECEnergy`), and the
intEC signature hashes it (formal-basis files are recomputed once; Brune and Park signatures unchanged); without it the
6.049 MeV capture still differed by 0.48. (2) Only final states above the lowest level of their block are recomputed:
the lowest level's formal parameters are its Brune ones (B_c is taken there), and recomputing them sent `tests/13N`'s
ground state through the gradient fall-back. (3) The analytic gradient falls back to finite differences in the formal
basis when such a final state exists (`AZUREGrad.cpp`, `FormalPartsOutsideAdjoint`). (4) The default `--ignore-externals` still makes
the bases differ (below).

**Evidence.** Engine before (9d50667f) and after (09e7b3d0), GCC 11.5, Release, GUI off, all 179 `.azr` under `tests/`,
menu 1, one thread, three bases, every file compared number by number (`classify.py`):

| basis | files | identical | changed |
|---|---|---|---|
| Brune | 689 | 689 | none |
| Park | 689 | 689 | none |
| formal | 662 | 649 | (a) 8 `parameters.out` (`thick_target_resonance`, `shift_sqrt_subpoints/models` × 7): a photon channel to the level's own state printed as μ or Q = 0 instead of G = 0 (the static-moment test now uses the Brune energy; the formal one is not the state's); 3 `run.log` (`screening_6li`): the radiative-width warning for two non-converging background levels at 1100 and 694 MeV is gone (it had read formal amplitudes as Brune ones); `beam_profile_kernel` `param.par` (22 formal photon amplitudes, up to 0.53 relative: γ^int from Brune final states) and `parameters.out`, whose photon widths to the 6.049 MeV state now reproduce the input: 0.000053 → 0.000330 meV (input 0.33 meV), 77.52 → 52.70 meV (52.7), 189.5 → 244.0 meV (244), 0.0233 → 0.0270 meV (0.027), 1.737 → 1.900 meV (1.9), 59.96 → 49.00 meV (49). No observable file changed: no menu-1 run of `tests/` captures to a final state above the lowest of its group. |

Brune against formal on `beam_profile_kernel`, ¹²C(α,γ) at 90° to the five final states (exit keys 4–8), 111 lab
energies E_α = 1.0–6.5 MeV, `--no-ignore-externals`, largest relative difference over the energies:

| final state | 9d50667f | 09e7b3d0 |
|---|---|---|
| 0⁺ g.s. (R = 4), 7.117 1⁻ (5), 6.917 2⁺ (6), 6.130 3⁻ (7): lowest of their groups | 0 (11 printed digits) | 0 (unchanged in both bases) |
| 0⁺ 6.049 MeV (R = 8), second of four 0⁺ levels | 0.808 | 1.5 × 10⁻¹¹ (Brune basis unchanged) |

With the default `--ignore-externals` the bases still differ: the `beam_profile_kernel` data (ground-state capture,
menu 1) by 9.7 × 10⁻² before and after; with `--no-ignore-externals` by 2.4 × 10⁻¹¹. That option drops, at a data point,
the external resonant capture of every level whose internal photon amplitude is below 10⁻⁸ MeV^½ (`AMatrixFunc.cpp`,
channel capture); which levels those are depends on the basis, since a formal level mixes the Brune ones. Not changed
here (HANDOFF, open).

Checks (GCC 11.5, GUI off, Condor): no compiler warnings; `tests/upgrade/check.sh` PASS (goldens unchanged);
`ctest -LE slow` 109/109; `tests/run_tests.sh` 133/133 (with the new `beam_profile_kernel/check.sh`, which fails on
9d50667f: 0.448 for R = 8); 1 vs 8 threads byte-identical on T1, T2, T3 in each basis; `physics_level_parameters_unit`
section 7 (observed_levels against hand values, 1.1 × 10⁻¹⁶). T2 benchmark, both engines alternately on one host
(d12chas084), mean of three, before → after: calc 2.45 → 2.22 s (1 thread), 0.42 → 0.44 s (8); fit200 133.1 → 131.8 s
(1), 55.8 → 56.2 s (8); peak RSS equal (`tests/upgrade/bench/T2.tsv`).

### Step 5 (branch `feature/upgrade-1-coupling`, 2026-11-09)

Design section 9 (da3703d2); physics functions and tests (96917f36); the engine switch (563a7172).
`CalcAngularDists`, `CalcCaptureAnalyzingPower` and `ThmAngular`'s Z̄ call `ObservableCoupling` (3j, 6j and a
new 9j of `AngularCoupling`, no GSL); `AngCoeff::Wigner9j` deleted. `CNuc.cpp` 3135 → 3061 lines.

**Evidence.** (a) `physics_observable_coupling`: the 9j against its zero reduction (117649 checks), the six-3j sum
(74088, 9002 non-zero) and its orthogonality (7415), to 6.7 × 10⁻¹⁶. The particle coefficients match the m-scheme
angular distribution to 1.1 × 10⁻¹⁵ (6750 cases, = (1/π) Σ_k c_k P_k). The photon coefficients match to
1.4 × 10⁻¹⁴ (31008 cases, = 8 Σ_k c_k P_k, with the helicity amplitude L̂ (L μ j_f m_f | J ν) d^L_{μλ} times λ
for M). (b) `physics_observable_coupling_gsl`: every coefficient (3.7 × 10⁶ over seven kinds) matches the old GSL
expression to 4 × 10⁻¹⁵ of its kind's largest value, with no keep-test mismatch. (c) 179 `.azr` × 3 bases against
the step-4 engine (a18fe89e code; Coulomb table format 1): 2063/2069 files identical. The six changed files are
`ec_e1m1_phase` `s_1` R = 3 and `gamma_lab_frame` `cm` R = 2, in each basis. Every changed number is within
3.2 × 10⁻¹⁵ of its column's largest value: classified (b). `ec_e1m1_phase`'s sign flips are residues of
10⁻²⁵–10⁻¹⁹ in a cross section that vanishes at 0° and 180°. ctest 113/113; `run_tests.sh` 133/133 (no
expected value moved, so there was no re-baselining commit); 1 vs 8 threads identical (12/12). One-thread
calculation: 13N 7.24 → 7.26 s, 13N_capture_ay 1.94 → 2.06 s (mean of three).

### Step 6a (branch `feature/upgrade-1-pathways`, 2026-11-09)

Design section 10 (365c325f); physics (8f92d046); engine (416bd630). `SortPathways` builds the System, calls
`reaction_pathways` per entrance pair and fills the trees from its result; `CNuc.cpp` 3061 → 2913 lines.

**Evidence.** (a) `physics_pathways_unit`: 1011 checks on a synthetic System (61 internal, 30 external pathways; 47 and
0 with R-matrix capture). (b) `physics_pathways_from_cnuc`: the trees against a verbatim copy of the old
`SortPathways` on every `.azr` under `tests/` (179), Brune and formal, every particle and beta pair an entrance:
51524 nodes identical over 347 loads (11 skipped: models `Fill` refuses, by design, and formal transformations that
do not exist); nothing is written. (c) 179 × 3 against the merged head db972c84: 2040/2040 identical. ctest
119/119, `run_tests.sh` 133/133, 1 vs 8 threads identical.

### Steps (d), (c) (branch `feature/upgrade-cnuc-split-3`, 2026-11-10)

Done, each in its own commit, every one bit-identical (179 `.azr` × 3 bases against `tab-p3f`: 689 / 662 / 689 files
identical after each), with `scripts/check.sh --quick` (ctest −LE slow 142/142, then 144/144 with the two io tests):

- **(d) output and checks**: `src/CNucReport.cpp` (the four check files, `parameters.out`, `relativistic.out`) and
  `src/CNucChecks.cpp` (the identical-channel refusal and warning, the Legendre-order warning, the
  separation-energy/mass check, the radiative-width warning), reading the nucleus through const overloads of the
  accessors; `CheckRadiativeWidths` split into `ObservedWidthsAt` (stays) and `WarnLargeRadiativeWidths`;
  `PrintCompoundFromParams`, `NucleusOfPair` deleted.  The check files are written by no run of the 179: T1, T3, 13N
  and `beam_profile_kernel` were run with the four checks on (file and screen), base and branch, identical.
- **(d) transformation messages**: `src/CNucTransformMessages.cpp`, the warnings of `TransformIn` and
  `TransformToBrune` (same wording).  Three 13N variants trigger the three column notes, identical to the base;
  the external-capture, eigen-failure and not-converged messages are reached by no input.
- **(c) input**: `NucLine` is `io/`'s `LevelRow` (`LevelRows.h`) with `ReadLevelSection`.  `Fill` builds from the
  records.  ctest `io_level_rows` (10554 rows of 189 files against a verbatim copy of the old reader, every field);
  six malformed variants give the base build's messages byte for byte.  The `<abundances>` line is read by
  `convert/`'s `ReadAbundanceLine` (`azure2/convert/AbundanceLine.h`), the parser the converter to format 3 uses
  too; the `io::ReadAbundanceSection` this step first wrote was a second parser of the same line and was removed
  when the converter was merged (the old grammar goes with step 3d).
- **(c) dead code**: the two unreachable branches of `ParseExternalCapture` and the write-only max-l members.

`CNuc.cpp` 2907 → 1942 lines, `CNuc.h` 306 → 281.  Still in `CNuc.cpp` and its stream: `Fill` obtains the
`<levels>` stream from `Config` and hands it to `io/` unread; `ExitPairKeysInFile` (a second reader of the first three
tokens of segment lines) is not moved, see "Remaining steps".

### Remaining steps

- **Step 4** is done (section 6 of `external_capture.md`, merged before this branch): N_f^½ and γ^ext come from
  `physics/external_capture` through `src/CNucExternalCapture.cpp`, the final state's Brune amplitude and energy are
  set there, and the `--ignore-externals` mask is defined on the input (Brune/Park) width in every basis
  (`ALevel::ZeroInputWidth`). What is left in `CNuc.cpp` is one adapter call `SqrtNFFactor(...)` in each of
  `TransformIn`, `TransformToBrune`, `CalcShiftFunctions`, `FillCompoundFromParams`,
  `FillCompoundFromParamsPhysical`; the last two belong to the parameter vector (step 7).
- **`ExitPairKeysInFile`**: `ParseExternalCapture` needs the exit keys of the active segments before the data are
  read.  Either `Fill`'s callers pass them (EData knows the segments, but reads them after `Fill`), or the segment
  readers (`SegLine`, `ExtrapLine`) provide them; both change behaviour on a malformed segment line, so they wait for a
  decision.
- **Step 6b** (consumers of the pathway trees read the physics result): design in section 11 — the inventory of the
  readers (320 tree accesses outside `CNuc.cpp`, 86 in it), the immutable `PathwayTable` built by `CNuc::Initialize`
  from `reaction_pathways` and the coupling coefficients and shared by every working copy, and eight bit-identical
  steps ending with the deletion of the trees.  Done: steps 1, 2, 4, 5, 6 and 7 (§11.12); open: steps 3 and 8
  (they edit `EData*`/`EPoint*`).
- **Step 3b (iii)** (hybrid channels as (α, l, j)): today's hybrid potential is central and keyed by the pair
  (`NuclearPotentialManager::tagFor(pairKey)`), so the channel functions do not depend on j; (α, l, j) entries need a
  spin-orbit term first (a new physics feature: `<potential>` grammar, `NuclearPotentialManager`, `CoulFunc`, the
  table key in `ChannelTable`), with its design note.
  The design note is `spin_orbit_potential.md` (key (α, l, 2J) for a pair with one spinning particle; steps in its §5).
- **Step 3b (iv)**, **step 7** (class (e): `FillMnParams`, `FillCompoundFromParams*`, `GetTransformParams`,
  `ConvertAmplitudeBasis`, the abundance members, penalties, `Initialize`, `Clone`) and the deletion of `CNuc`: not
  started (another agent changes the parameter vector).

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

## 9. Step 5: design

Step 5 moves the coupling coefficients of the observables into physics/: the Legendre coefficients of the angular
distributions (`CNuc::CalcAngularDists`), the capture analyzing-power coefficients (`CNuc::CalcCaptureAnalyzingPower`)
and the THM angular coefficient (`ThmAngular.cpp`, the same Z̄). They are then computed without GSL, from the 3j and 6j
of `AngularCoupling` and a new 9j. `SortPathways` is not in this step; it moves with the pathway trees (step 6), because
its output is those trees on `PPair`. The pathway enumeration, the thresholds that decide which terms are kept, the
`KLGroup`/`Interference`/`CaptureAyTerm` objects and the lazy build of the capture table under its mutex stay as they are.

**The functions** (`physics/include/azure2/physics/ObservableCoupling.h`; all angular momenta doubled, as in physics/).
With x̂ = (2x + 1)^½ and the Racah W(abcd; ef) = (−1)^(a+b+c+d) {a b e; d c f}:

    Z̄(l₁ j₁ l₂ j₂; s L)    = l̂₁ l̂₂ ĵ₁ ĵ₂ (l₁ 0 l₂ 0 | L 0) W(l₁ j₁ l₂ j₂; s L)                                      (9.1)
    Z̄γ(L₁ j₁ L₂ j₂; j_f L) = L̂₁ L̂₂ ĵ₁ ĵ₂ (L₁ 1 L₂ −1 | L 0) W(L₁ j₁ L₂ j₂; j_f L)                                 (9.2)

and the interference coefficients the engine stores. A particle exit is (−1)^(s′−s) Z̄(l₁ j₁ l₂ j₂; s L) Z̄(l′₁ j₁ l′₂ j₂; s′ L) / 4.
A photon exit is (−1)^(1+s−j_f) Z̄(l₁ j₁ l₂ j₂; s L) Z̄γ(L₁ j₁ L₂ j₂; j_f L) / 4. The particle exit with an unobserved primary
(Brune and deBoer, PRC 102, 024628, Eq. 22, in AZURE2's light-first channel-spin order, as the comment at
`CNuc.cpp` 1930 derives) is a product of three Racah W. For capture A_y, a_k and b_k are Seyler and Weller,
PRC 20, 453, Eqs. (20) and (21), with the (−1)^(s−s′) of the light-first coupling (audit item 9). The pruning of
`CalcCaptureAnalyzingPower` (a factor below 10⁻¹² or a 9j below 10⁻¹⁴ gives zero) goes into the function, so it
returns the same zeros. `wigner_9j` is added to `AngularCoupling` as the sum over x of three 6j (Edmonds 6.4.3).

**Not in step 5.** The remaining `AngCoeff` callers evaluate at each point, not at setup, and stay on GSL until the
observables move (Phase 2): `PolarizationFunc` (Clebsch-Gordan, Y_lm), `ECIntegral` (Clebsch-Gordan, Racah of the
multipole operators), and `LegendreP1` in `AMatrixFunc`/`GenMatrixFunc`. `AngCoeff`'s 9j then has no caller and is
deleted.

**Numbers.** GSL's 3j and 6j and those of `AngularCoupling` differ in the last bits. For 2j ≤ 12, 69 % of the
non-zero 6j differ, by at most 7.3 × 10⁻¹⁵ relative; for 2j ≤ 20, 97 % of the 3j differ. `AngularCoupling` agrees
with exact values to 10⁻¹⁵, and GSL returns 6.9 × 10⁻¹⁸ for some 6j that are exactly zero. The step therefore
changes the coefficients by about 10⁻¹⁵ relative, and each observable by the same amount. A fit can amplify that to
its tolerance. It is checked, and if needed re-baselined, as in §3 item 5.

**Verification.**
(a) physics tests, on things that do not use the code under test. The 9j is checked against exact values, against
its reduction to a 6j when one argument is zero (Edmonds 6.4.14), against its definition as a sum over all
projections of six 3j for 2j ≤ 4, and against its orthogonality. (9.1) and (9.2) are checked against the m-scheme
sum of Clebsch-Gordan products that defines them, for l, L ≤ 3 and 2j ≤ 6.
(b) an engine-side ctest compares every new function with the engine's old expression (`AngCoeff`, GSL) for all
quantum numbers up to l, L ≤ 4, 2j, 2s ≤ 8 and k ≤ 6. The limit is 10⁻¹⁴ of the largest coefficient of its kind,
and the zeros must agree.
(c) the engine before (9be284f0) and after, on all 179 `.azr` in three bases, every changed file classified. An
observable within 10⁻¹² relative of its column's largest value is (b) rounding. A fitted value within the fit's
tolerance is (b). Anything else stops the step. Then `ctest`, `run_tests.sh` and 1 vs 8 threads. A test whose
expected value moves only by the amounts in (c) is re-baselined in its own commit, with the classification.

## 10. Step 6a: design

The pathway trees on `PPair` (`Decay` → `KGroup` → `MGroup`/`ECMGroup`) are read by every observable (`AMatrixFunc`,
`GenMatrixFunc`, `RMatrixFunc`, `EPoint`, `EData`, `AZUREAPI`, `CNucExternalCapture`: about 340 calls of
`GetDecay`, `GetKGroup`, `GetMGroup`, `GetECMGroup`, `GetKLGroup`). `GetPair`/`GetJGroup` are called about 730 times
in 23 files. Moving all of them at once is the observables move of Phase 2, so step 6 is split. Step 6a moves what
decides the trees, the enumeration of `SortPathways`, into physics/. The trees, their consumers and
`CalcAngularDists` (which reads the trees) stay. Step 6b then moves the consumers one at a time onto the physics
result and deletes the trees when the last one has moved.

**The function** (`physics/include/azure2/physics/Pathways.h`):

    reaction_pathways(System, entrance partition α, final states, options) → [ExitPathways]

The System is the R-matrix structure (`CNucSystemContent::Structure`: blocks are the J-groups in the R-matrix). A
final state is described by its 2J, parity, photon partition, multipole mask (E1, M1, E2) and particle channels
(partition, l, 2s). Final states can lie in J-groups outside the R-matrix, so they are not System blocks. The options
are the R-matrix-capture flag (`USE_RMC_FORMALISM`: no photon exit pathways). The result has one entry per exit
partition α′, in order of first appearance. Each entry holds the channel-spin groups (2s, 2s′ and, for the
unobserved-primary groups, 2s₂′), each with its internal pathways and its external pathways:
- internal pathway: block b, entrance channel k, exit channel k′, statistical spin factor (2J + 1)·I₁I₂ factor;
- external pathway: multipole, L, entrance l, compound 2J, final state, final channel and, for channel capture,
  the channel-capture group and pathway it continues and the compound block's photon channel.

Every list is in the order `SortPathways` appends: its loops (entrance pair, exit pair, s, s′, s₂′, group, channel,
channel′) run over System indices in the same order, since the System lists the partitions, R-matrix groups and their
channels in `CNuc` order (`physics_system.md`). The rules are those of `SortPathways`: the four exit kinds (beta
entrance with l ≤ 1, elastic, particle exit with s₂′, photon), the deduplication of channel-capture groups on the
pathway triple, the radiation type from the parities, the multipole mask and J triangle, and the coupling rule
`couples` of `ExternalCapture.h` (the same expression).

**The engine.** `SortPathways` builds the System, calls `reaction_pathways` for each entrance pair, and fills
`Decay`/`KGroup`/`MGroup`/`ECMGroup` from the result through the adapter map (block → J-group number, block channel
→ group channel number, final state → J-group and level number). The trees are therefore the same objects in the same
order with the same numbers, so every output is expected bit-identical.

**Verification.**
(a) A physics test on synthetic Systems checks properties the enumeration must have, without reusing its loops:
- every internal pathway joins two channels of one block with the stated channel spins;
- every pair of block channels with entrance partition α appears in exactly one group per exit kind;
- every external pathway satisfies the parity, mask, J and coupling rules;
- the count of each kind equals a direct count.
(b) An engine-side ctest keeps a verbatim copy of the old `SortPathways` (at db972c84). It builds both trees on
T1–T3 and on the projects with UPOS, external capture, beta-delayed and capture-A_y pathways
(`upos_pair_index`, `light_excitation`, `13N`, `13N_capture_ay`, `ec_e1m1_phase`, `beta_delayed_alpha`,
`capture_ay_coupling`). Every field and every position must be equal.
(c) 179 `.azr` × 3 bases against the merged head before 6a (db972c84): every file identical; anything else stops the
step. Then ctest, `run_tests.sh` and 1 vs 8 threads.

## 11. Step 6b: design

Step 6b moves every reader of the pathway trees on `PPair` (`Decay` → `KGroup` → `MGroup`/`ECMGroup`, and
`KLGroup` → `Interference`) onto one immutable table built from the result of `reaction_pathways` (§10) and the
coupling coefficients (§9), and deletes the trees, their classes and their adapter with the last reader. It changes
no number: every step keeps each output file of the 183 projects × 3 bases byte-identical to `tab-trial8` (704 /
677 / 704 files, Brune / formal / Park), and 1 thread against several. Line numbers below are those of the start
commit cb8f5912 (the integration branch `feature/upgrade-4-data`). `src/EData.cpp`, `src/EPoint.cpp` and
`src/ESegment.cpp` are being cut into smaller files by verbatim moves (in work, `edata_split.md`): consumers are named
by function, and the file given is the file at cb8f5912.

### 11.1 Inventory of the readers

Counted with `work/pathnote/inventory.py` (a name search attributing each occurrence to the enclosing function
definition). "Accesses" counts the container accessors of the trees: `GetDecay`, `NumDecays`, `IsDecay`,
`GetKGroup`, `NumKGroups`, `IsKGroup`, `GetMGroup`, `NumMGroups`, `IsMGroup`, `GetECMGroup`, `NumECMGroups`,
`GetKLGroup`, `NumKLGroups`, `IsKLGroup`, `GetInterference`, `NumInterferences`, `FlatInterferences`,
`InterferencePathways`, `NumInterferenceOrders`, `DiagonalPathways`, `TablesBuilt`, `GetCaptureAyTerm`,
`NumCaptureAyTerms`, `IsCaptureAyBuilt` and the writers (`Add*`, `Set*`, `BuildFlatInterferences`, `ClearDecays`).
Counting only the five `Get…Group`/`GetDecay` getters gives 211 in the readers and 47 in `CNuc.cpp` (258, the
"about 270" of "Remaining steps"). The node fields read (`GetJNum`, `GetChNum`, …) are listed by name; `GetS`, `GetL`,
`GetJ` and `GetRadType` are not counted, since `JGroup` and `AChannel` have getters of the same names. "Hot" is inside
the point loop of an evaluation (forward or gradient); "OMP" says whether it runs inside an OpenMP region.

| file | function | def. / accesses at | accesses | node fields read | hot | OMP |
|---|---|---|---|---|---|---|
| `src/GenMatrixFunc.cpp` | `GenMatrixFunc::CalculateBareCrossSection` | 38 / 49–536 | 96 | Decay pair, flat interferences, `TablesBuilt`; KGroup s, s′, s₂′, cross-spin; MGroup group, ch, ch′, stat. factor; ECMGroup l, L, J, E/M; KLGroup k, L; Interference m₁, m₂, type, Z̄₁Z̄₂, Z̄₁Z̄₂ UPOS | yes | yes (point pass, `EData.cpp` 4444) |
| `src/GenMatrixFunc.cpp` | `GenMatrixFunc::CalculateAmplitudeMatrix` | 801 / 814–828 | 10 | KGroup cross-spin; MGroup group, ch, ch′ | yes (particle A_y) | yes |
| `src/GenMatrixFunc.cpp` | `GenMatrixFunc::CalculateCaptureAnalyzingPower` | 845 / 857–877 | 6 | capture-A_y terms; calls `CNuc::CalcCaptureAnalyzingPower` (a write, §11.2) | yes (capture A_y) | yes |
| `src/GenMatrixFunc.cpp` | `GenMatrixFunc::BuildCoulombInterference` | 923 / 927–947 | 7 | diagonal pathways, `TablesBuilt`; MGroup group, ch, ch′, stat. factor | yes (elastic) | yes |
| `src/AMatrixFunc.cpp` | `AMatrixFunc::CalculateTMatrix` | 513 / 528–599 | 16 | MGroup group, ch, ch′; ECMGroup final group, level, channel, channel-capture flag and back-reference, photon channel | yes | yes |
| `src/AMatrixFunc.cpp` | `AMatrixFunc::PointTBar` | 672 / 703–1113 | 62 | as `CalculateBareCrossSection`, plus interference pathways and orders, capture-A_y terms (calls `CalcCaptureAnalyzingPower`) | yes (gradient) | yes (`AZUREGrad.cpp` 1280, 1387) |
| `src/AMatrixFunc.cpp` | `AMatrixFunc::TBarAdjoint` | 1129 / 1147–1315 | 18 | as `CalculateTMatrix` | yes (gradient) | yes |
| `src/AMatrixFunc.cpp` | `AMatrixFunc::CollectBilinearTerms` | 1852 / 1859–1875 | 11 | MGroup group, ch, ch′; ECMGroup count | yes (gradient) | yes |
| `src/AMatrixFunc.cpp` | `AMatrixFunc::TBarProducts` | 2690 / 2692–2695 | 4 | flat interferences, pathways, orders, `TablesBuilt` | yes (gradient) | yes |
| `include/AMatrixFunc.h` | `TBarProducts`, `tbar_decay_` | 350–355 | 1 (+ `Decay *` in two declarations) | — | — | — |
| `include/GenMatrixFunc.h` | `BuildCoulombInterference`, `rt_decay_` | 178, 188 | 0 (`Decay *` in two declarations) | — | — | — |
| `src/RMatrixFunc.cpp` | `RMatrixFunc::CalculateTMatrix` | 340 / 348–413 | 16 | as `AMatrixFunc::CalculateTMatrix` | yes (the R-matrix route: `<config>` first line `false`, `EPoint.cpp` 1924–1929) | yes |
| `src/EData.cpp` | `WalkECAmplitudes` (template, the intEC signature and count) | 1997 / 2011–2052 | 11 | KGroup s; ECMGroup l, L, J, E/M, final channel, channel capture and back-reference | no (set-up) | no |
| `src/EData.cpp` | `EData::CalculateECAmplitudes` (writes `intEC.dat`) | 2231 / 2394–2395 | 5 | group and pathway counts | no | its loop at 2361 is a parallel region; these lines are after it |
| `src/EData.cpp` | `EData::InitializeComponentSegments` (component intEC lines) | 2447 / 2680–2738 | 6 | group and pathway counts | no | as above (2642) |
| `src/EPoint.cpp` | `EPoint::GetECAmplitudeWithShift` | 594 / 611–615 | 3 | ECMGroup final group, level; then `ECPathwayAmplitude` | yes (with external capture, on an energy the point has not stored) | yes |
| `src/EPoint.cpp` | `ECPathwayAmplitude` (static) | 1679 / 1711–1741 | 3 | KGroup s; ECMGroup l, L, J, E/M, final channel, channel capture and back-reference | set-up, and as above | yes |
| `src/EPoint.cpp` | `EPoint::CalculateECAmplitudes` | 1754 / 1765–1768 | 6 | group and pathway counts | no (set-up) | yes (`EData.cpp` 2361, 2642) |
| `src/CNucExternalCapture.cpp` | `CNuc::SetUpFormalCaptureMasks` | 230 / 241–248 | 10 | ECMGroup channel capture, back-reference, photon channel; MGroup group | no (set-up) | no |
| `api/src/AZUREAPI.cpp` | `AZUREAPI::GetECIntegrals` | 1727 / 1751–1797 | 9 | KGroup s; ECMGroup l, L, J, E/M, final channel, channel capture and back-reference | no (API call) | no |
| `src/CNucReport.cpp` | `WritePathwaysCheck` | 126 / 144–212 | 14 | every field of KGroup, MGroup, ECMGroup | no (check file) | no |
| `src/CNucReport.cpp` | `WriteAngularDistributionsCheck` | 246 / 262–274 | 6 | KLGroup k, L; Interference m₁, m₂, type, Z̄₁Z̄₂ | no (check file) | no |
| `src/CNuc.cpp` | `CNuc::SortPathways` (fills the trees) | 890 / 928–969 | 18 (writers) | — | no | no |
| `src/CNuc.cpp` | `CNuc::CalcAngularDists` (KLGroups, flat and diagonal lists) | 1020 / 1024–1160 | 58 (reads and writes) | MGroup group, ch, ch′, stat. factor; ECMGroup l, L, J, E/M; KGroup s, s′, s₂′, cross-spin | no | no |
| `src/CNuc.cpp` | `CNuc::CalcCaptureAnalyzingPower` (capture-A_y terms, lazily) | 1190 / 1196–1263 | 10 (reads and writes) | MGroup group, ch, ch′; ECMGroup l, L, J, E/M; KGroup s | **yes** (first capture-A_y point) | **yes** |

Readers outside `CNuc.cpp`: 320 accesses (`GenMatrixFunc.cpp` 119, `AMatrixFunc.cpp` 111, `CNucReport.cpp` 20,
`EData.cpp` 22, `RMatrixFunc.cpp` 16, `EPoint.cpp` 12, `CNucExternalCapture.cpp` 10, `AZUREAPI.cpp` 9,
`AMatrixFunc.h` 1); in `CNuc.cpp` 86. Not readers: `src/AZUREGrad.cpp` (0 accesses; it passes `TMatrixCotangent`
buffers shaped by `PointTBar`), the THM matrix functions, `PolarizationFunc` (it takes group and channel numbers),
and the other `src/EData*.cpp` files at cb8f5912. The class files themselves (`Decay`, `KGroup`, `KLGroup`, `PPair`,
`MGroup`, `ECMGroup`: 109 occurrences) and `tests/upgrade/pathways_from_cnuc_test.cpp` (114) go with the trees.

**The numbering is a contract.** The pathway positions (k, m) — 1-based spin group within its decay, and pathway
within its group, separately for internal and external — key per-evaluation and per-point state outside the trees:
`GenMatrixFunc::tmatrix_[decay][k][m]` and `ec_tmatrix_[k][m]` (`AddTMatrixElement` asserts the order), the
cotangents `TMatrixCotangent::t`, `ec` (`AMatrixFunc.h` 20), the point's `ec_amplitudes_`, `ec_energies_`,
`ec_shifted_` (`EPoint::GetECAmplitude(k, m)`), the keys of `ECAmplitudeCache::AmplitudeKey` (k, m, entrance, exit,
segment), the line order of `intEC.dat` and its signature (`WalkECAmplitudes` hashes the pathways in tree order), and
the indices stored in the interference and capture-A_y terms. The table keeps these numbers exactly; renumbering
them is a separate change (§11.9).

### 11.2 Writes to the trees after set-up

The trees are written by `CNuc::SortPathways` and `CNuc::CalcAngularDists`, both called once from `CNuc::Initialize`
(`CNuc.cpp` 700, 709), and by the test `pathways_from_cnuc_test.cpp` (it clears and rebuilds them). One writer runs
after set-up, during the point pass:

- **`CNuc::CalcCaptureAnalyzingPower`** (`CNuc.cpp` 1190–1267), called from
  `GenMatrixFunc::CalculateCaptureAnalyzingPower` (`GenMatrixFunc.cpp` 864) and `AMatrixFunc::PointTBar`
  (`AMatrixFunc.cpp` 865) by the first capture-A_y point of a decay, on whichever OpenMP thread computes it. It
  appends to `Decay::captureAy_` and sets `captureAyBuilt_`, under a function-level `static std::mutex` (`CNuc.cpp`
  1199), a mutable static that ARCHITECTURE §9 forbids. The threads of the point pass and of the gradient
  (`AZUREGrad.cpp` 1282: one `AMatrixFunc` per thread on the same `CNuc`) share the trees, so the table is built once
  per `CNuc` object; a clone made before the build has none and builds its own. The terms are built for the Legendre
  order `point->GetMaxLOrder()` of the first caller. For a differential point that is `configure.maxLOrder`
  (`EPoint::Initialize`, `EPoint.cpp` 677), so all callers ask for the same order; it is not established here that no
  caller (a target-integration sub-point of an A_y point) has `GetMaxLOrder() = −1`, which would build an empty table
  and mark it built.

Mutable fields of the nodes and who writes them: `PPair::decays_` (`AddDecay`, `ClearDecays`); `Decay::kgroups_`,
`klgroups_` (`AddKGroup`, `AddKLGroup`), `flat_`, `flatPathways_`, `flatOrders_` (`BuildFlatInterferences`), `diag_`
(`SetDiagonalPathways`), `tables_built_` (`SetTablesBuilt`) — all `CalcAngularDists` or `SortPathways`;
`Decay::captureAy_`, `captureAyBuilt_` — `CalcCaptureAnalyzingPower`, after set-up; `KGroup::mgroups_`,
`ec_mgroups_`, `MGroup::statspinfactor_` — `SortPathways`; `KLGroup::interferences_`, `interferenceIndex_` —
`CalcAngularDists`; `ECMGroup::statspinfactor_` — no writer (stays 0.0; `ECMGroup::GetStatSpinFactor` has no caller);
`MGroup::tmatrix_`, `ECMGroup::tmatrix_` — never initialised, written or read.

### 11.3 What the nodes hold beyond the physics result

Every field is structural: none depends on the parameters θ, the nuisances ν or the energy, so the whole content is
set-up state and fits an immutable table. Some fields depend on the configuration (`maxLOrder`, R-matrix capture,
which pairs are entrances), which is fixed before `CNuc::Initialize`.

| node, field | in `reaction_pathways`' result | derived at set-up | depends on θ or E |
|---|---|---|---|
| `Decay::pair_` | `ExitPathways::partition` (+1: `pairOfPartition` is the identity) | — | no |
| `KGroup::s_`, `sp_`, `sp2_`, `hasSp2_` | `SpinGroup::two_s`, `two_sp`, `two_sp2`, `has_sp2` (halved; exact) | — | no |
| `MGroup::jnum_`, `ch_`, `chp_` | `InternalPathway::block`, `k_in`, `k_out` through `CNucSystemMap` (`groupOfBlock`; `channelOfBlockChannel[b][k] = k + 1` as built, `CNucSystem.cpp` 102–146) | — | no |
| `MGroup::statspinfactor_` | `stat_spin_factor`, 0.0 when `has_stat_spin_factor` is false (beta entrance) | — | no |
| `ECMGroup::radtype_`, `mult_`, `li_`, `ji_`, `chf_` | `X`, `L`, `l`, `two_J / 2.`, `final_channel + 1` | — | no |
| `ECMGroup::jGroupNum_`, `levelNum_` | `final_state` through the final-state list of `SortPathways` (J-group, level) | — | no |
| `ECMGroup::ischancap_`, `chdecay_`, `chkgroup_`, `chmgroup_` | `channel_capture`, `cc_exit + 1`, `cc_group + 1`, `cc_pathway + 1` (the tree stores 0 for direct capture; the result carries the triple for both) | — | no |
| `ECMGroup::internalChannel_` | `photon_channel` through the map, 0 when −1 | — | no |
| `KLGroup` (k, L), `Interference` (m₁, m₂, type, Z̄₁Z̄₂, Z̄₁Z̄₂ UPOS) | — | `CalcAngularDists` from the pathways' J, l, l′, s′, E/M and the exit partition's spins, by `particle_interference`, `photon_interference`, `upos_interference` (§9), kept when \|Z\| > 10⁻¹⁰, for L = 0 … `maxLOrder` | no |
| `Decay::flat_`, `flatPathways_`, `flatOrders_` | — | `BuildFlatInterferences`, from the KLGroups in their order | no |
| `Decay::diag_` | — | the internal pathways with entrance = exit channel (s, l, stat. factor) | no |
| `Decay::captureAy_` | — | `capture_ay_coefficients` (§9) over all pathway pairs of a photon decay, kept when \|a\| or \|b\| ≥ 10⁻¹⁰, for k = 0 … max L | no (built lazily, §11.2) |
| `Decay::tables_built_`, `captureAyBuilt_`; `KLGroup::interferenceIndex_` | — | flags and a set-up search index | — |
| `MGroup::tmatrix_`, `ECMGroup::tmatrix_`, `ECMGroup::statspinfactor_` | — | unused | — |

So the table is the physics result in the engine's numbering plus three derived lists per decay (interferences,
diagonal pathways, capture-A_y terms). Nothing of it belongs in the per-evaluation state (`Results`, the plan's
tasks) or in the per-stage energy tables.

### 11.4 The table

**Two layouts.**

(a) *The physics result itself.* Keep `std::vector<physics::ExitPathways>` per entrance (nested: exits → spin groups →
internal and external pathway vectors, System indices, 0-based), add the three derived lists as physics types, and
let the engine translate a block index to its J-group number through `CNucSystemMap` in each loop. One representation,
in physics/. The hot loops pay one map lookup per pathway (`groupOfBlock[block]`, a vector index) and keep three
levels of vectors.

(b) *A flat engine table of the same result.* Per entrance, one array per kind of row (exits, spin groups, internal
pathways, external pathways, interferences, interference pathways, diagonal pathways, capture-A_y terms), each in
the order of the trees, with offsets and counts linking them; rows carry the engine's numbers (J-group, channel,
level) beside the System indices. Declared in the engine, because it carries engine numbers.

**Recommendation: (b).** It removes the map lookups and the nested vectors from the hottest loops (§11.8), and its
row index is a global pathway index per entrance, which is what Phase 5 needs to store T-matrix elements in one flat
buffer instead of `tmatrix_[decay][k][m]` (§11.9). The cost is that the engine numbers are a second copy of the
System's indices; they go when the per-point tables and the level matrix are keyed by System block (`channel_table.md`
§2: Phase 4/5), and the table then moves to physics/ without them. What would change the recommendation: if the
per-point channel tables (`EDependentTables::Slot(jGroup, channel)`) and `AMatrixFunc`'s per-group matrices were keyed
by System block before 6b starts, the engine numbers would be unneeded and (a), read directly, would be the simpler
choice.

**Types** (`include/PathwayTable.h`, `src/PathwayTable.cpp`; plain structs, indices only, no pointers):

```cpp
namespace azure2 {
struct InternalRow {           // an MGroup
  int group, ch, chp;          // J-group number and its entrance and exit channel numbers (CNuc's, 1-based)
  double stat_spin_factor;     // (2J + 1) / ((2 I1' + 1)(2 I2' + 1)) of the exit pair; 0.0 for a beta entrance
  int block, k_in, k_out;      // physics::InternalPathway (System indices, 0-based)
};
struct ExternalRow {           // an ECMGroup
  char rad_type;               // 'E' or 'M'
  int mult, l;                 // multipolarity L; entrance orbital angular momentum l
  double J;                    // compound J (two_J / 2.)
  int final_group, final_level, final_channel;   // the final state's J-group, level and channel (1-based)
  bool channel_capture;
  int cc_exit, cc_group, cc_pathway;   // channel capture: the internal pathway it passes through, 1-based as the
                                        // check file prints them; 0 for direct capture (as the tree)
  int cc_row;                  // the same pathway as an index into EntrancePathways::internal; -1 for direct capture
  int photon_channel;          // the J-group channel of the block's photon channel (1-based), 0 if none
};
struct SpinGroupRow {          // a KGroup
  double s, sp, sp2;           // channel spins; sp2 = 0.0 unless has_sp2
  bool has_sp2;
  int first_internal, num_internal, first_external, num_external;
  bool CrossSpin() const { return has_sp2 && sp != sp2; }
};
struct ExitRow {               // a Decay
  int pair;                    // exit pair number (1-based)
  int first_group, num_groups;
  bool has_interferences;      // Decay::TablesBuilt: false for a beta entrance
  int first_interference, num_interferences, num_orders;
  int first_interference_pathway, num_interference_pathways;
  int first_diagonal, num_diagonal;
  int first_capture_ay, num_capture_ay;
};
struct EntrancePathways {
  std::vector<ExitRow> exits;
  std::vector<SpinGroupRow> groups;
  std::vector<InternalRow> internal;
  std::vector<ExternalRow> external;
  std::vector<FlatInterference> interferences;          // with a new field z_upos (Interference::GetZ1Z2_UPOS)
  std::vector<InterferencePathway> interference_pathways;
  std::vector<DiagonalPathway> diagonal;
  std::vector<CaptureAyTerm> capture_ay;
  int ExitOf(int exitPair) const;                       // 1-based exit position, 0 if the pair is no exit
};
class PathwayTable {
 public:
  const EntrancePathways &Entrance(int pair) const;     // pair 1-based; empty for a pair that is not an entrance
 private:
  std::vector<EntrancePathways> entrances_;             // one per pair, in CNuc order
};
}
```

`FlatInterference`, `InterferencePathway`, `DiagonalPathway` (now in `Decay.h`) and `CaptureAyTerm` keep their
fields and their 1-based (k, m) and move to `PathwayTable.h`; `FlatInterference` gains `z_upos`, so the UPOS and
angular-coefficient loops that walk `KLGroup`s today walk the flat list (it is built in the same order: KLGroups in
order of creation, interferences in order within each, and no (k, L) group is split, since `CalcAngularDists` loops
k → L → m₁ → m₂ and a (m₁, m₂, type) never repeats within one (k, L)).

**Who builds it, when.** One function, `BuildPathwayTable(const CNucSystem &, final states, options)`, called by
`CNuc::Initialize` where `SortPathways` and `CalcAngularDists` are called now, after the input transformation and
`SetUpECFinalBlocks`. It takes the System and its map (`BuildSystemFromCNuc(…, Structure)`), the final-state list
that `SortPathways` builds (`physics::FinalState` and its J-group and level), and options: R-matrix capture,
`maxLOrder`, the entrance pairs (`PPair::IsEntrance`, set while the data are read, before `CNuc::Initialize`) and the
(entrance, exit) pairs that need capture-A_y terms (§11.5). It calls `reaction_pathways` per entrance, appends the
rows in the order of the result, then the interferences (the loop of `CalcAngularDists`, reading J, l, l′, s′ and
E/M from the System: `two_J`, channel `l`, `two_s`, `kind`, the same integers and the same halved doubles that
`JGroup`/`AChannel` return), the diagonal pathways and the requested capture-A_y terms. Nothing writes it afterwards.

**Who owns it.** Not the System (it depends on the final states, the configuration and the entrances, none of which
are the R-matrix structure) and not the calculation plan (the plan is per data set and rebuilt when a shift changes;
the table is per run and is also read by the reaction rate, the extrapolations and the API). It is set-up state of the
run: `CNuc::Initialize` produces it as a `std::shared_ptr<const PathwayTable>`, and during 6b `CNuc` carries that
handle (`const PathwayTable &CNuc::Pathways() const`); `CNuc::Clone` copies the handle, so every working copy shares
one table instead of deep-copying the trees. In step 7, when `CNuc` is deleted, the handle moves to the object that
replaces it as the holder of the run's immutable set-up (the System and the tables built from it); the readers do not
change, only the object `Pathways()` is called on. `PPair` holds nothing.

**How a point and a thread reach it.** Through the objects they are given, never a global. Each matrix-function object
already receives its `CNuc *` (the point pass `Reset`s its per-thread `AMatrixFunc` with the point's `CNuc` at every
point, `EPoint.cpp` 1922–1926; the gradient builds one per thread, `AZUREGrad.cpp` 1282); a reader takes `const
EntrancePathways &P = compound()->Pathways().Entrance(aa)` once per point, before its loops. `EPoint` and `EData`
set-up functions already receive the `CNuc *` they read; new signatures take `const PathwayTable &` directly. Threads
only read.

**The hottest loop.** `AMatrixFunc::CalculateTMatrix` runs for every point of every evaluation (forward) and its
loop shape is shared by `RMatrixFunc::CalculateTMatrix`, `TBarAdjoint` and `CollectBilinearTerms`. Today
(`AMatrixFunc.cpp` 539–598, abridged):

```cpp
for (int ir = irStart; ir <= irEnd; ir++) {
  Decay *theDecay = compound()->GetPair(aa)->GetDecay(ir);
  for (int k = 1; k <= theDecay->NumKGroups(); k++) {
    for (int m = 1; m <= theDecay->GetKGroup(k)->NumMGroups(); m++) {
      MGroup *theMGroup = theDecay->GetKGroup(k)->GetMGroup(m);
      int jNum = theMGroup->GetJNum();
      int chNum = theMGroup->GetChNum();
      int chpNum = theMGroup->GetChpNum();
      JGroup *theJGroup = compound()->GetJGroup(jNum);          // the three pointers are not used below
      AChannel *entranceChannel = theJGroup->GetChannel(chNum);
      AChannel *exitChannel = theJGroup->GetChannel(chpNum);
      complex coulombPhaseEn = point->GetExpCoulombPhase(jNum, chNum);
      ...                                                         // phases, sqrt(P), GetUBilinear(jNum, chNum, chpNum)
      this->AddTMatrixElement(k, m, tmatrix);                     // or (k, m, ., ir) with R-matrix capture
    }
    for (int m = 1; m <= theDecay->GetKGroup(k)->NumECMGroups(); m++) {
      ECMGroup *theECMGroup = theDecay->GetKGroup(k)->GetECMGroup(m);
      ALevel *finalLevel = compound()->GetJGroup(theECMGroup->GetJGroupNum())->GetLevel(theECMGroup->GetLevelNum());
      ...
      if (theECMGroup->IsChannelCapture()) {
        MGroup *chanMGroup = compound()->GetPair(aa)->GetDecay(theECMGroup->GetChanCapDecay())
                                 ->GetKGroup(theECMGroup->GetChanCapKGroup())->GetMGroup(theECMGroup->GetChanCapMGroup());
        ...                                                       // chanMGroup->GetJNum(), GetChNum(), GetChpNum()
      }
      this->AddECTMatrixElement(k, m, tmatrix);
    }
  }
}
```

On the table, with every arithmetic line unchanged and in the same order:

```cpp
const EntrancePathways &P = compound()->Pathways().Entrance(aa);
for (int ir = irStart; ir <= irEnd; ir++) {
  const ExitRow &exit = P.exits[ir - 1];
  for (int k = 1; k <= exit.num_groups; k++) {
    const SpinGroupRow &g = P.groups[exit.first_group + k - 1];
    for (int m = 1; m <= g.num_internal; m++) {
      const InternalRow &r = P.internal[g.first_internal + m - 1];
      const int jNum = r.group, chNum = r.ch, chpNum = r.chp;
      complex coulombPhaseEn = point->GetExpCoulombPhase(jNum, chNum);
      ...                                                         // unchanged
      this->AddTMatrixElement(k, m, tmatrix);
    }
    for (int m = 1; m <= g.num_external; m++) {
      const ExternalRow &x = P.external[g.first_external + m - 1];
      ALevel *finalLevel = compound()->GetJGroup(x.final_group)->GetLevel(x.final_level);
      ...
      if (x.channel_capture) {
        const InternalRow &cc = P.internal[x.cc_row];
        ...                                                       // cc.group, cc.ch, cc.chp
      }
      this->AddECTMatrixElement(k, m, tmatrix);
    }
  }
}
```

**What stays identical, so that the bits stay.** The loops visit exits, spin groups and pathways in the order of the
trees (the rows are appended in that order, and the step-1 test below proves it), so every sum over pathways
(`CalculateBareCrossSection`'s temporary T-matrices grouped by (J, l, l′) in order of first appearance, the flat
interference sum into `rt_wl_` by Legendre order, the Coulomb-interference terms, the capture-A_y sums, the
cotangent sums of the gradient) adds the same terms in the same order. The operands are the same: integers, and
doubles that are exact halves (s, s′, J) or the same `stat_spin_factor` and coefficient doubles computed by the same
functions with the same arguments. The (k, m) passed to `AddTMatrixElement`, `GetECAmplitude`, the cotangent
buffers and the EC cache keys are the same numbers.

**Tempting, and not to be taken in 6b.**
- Skipping zero rows or reordering rows "for locality" (by block, by l): it reorders the sums.
- Merging the internal and external rows into one array, or renumbering m across groups: (k, m) is the contract of
  §11.1.
- Replacing the linear search of `IsTempTMatrix` by a precomputed (J, l, l′) group index per row: a sound later
  speed-up, but the group index must reproduce the order of first appearance; it is not part of 6b.
- Sharing the rows of a cross-spin UPOS group with the equal-spin groups whose pathways it repeats: each copy is
  its own (k, m) slot of the T-matrix storage.
- Fixing a defect of §11.9 inside a move: each fix is its own commit with its own classification.

**The gradient.** The gradient pass (`AZUREGrad.cpp`, `PointTBar`, `TBarAdjoint`, `CollectBilinearTerms`,
`TBarProducts`) needs from the table exactly what the forward pass reads, read-only: the shapes of the cotangent
buffers (groups × internal and external counts of one exit), the (J, l, l′) grouping, the flat interferences with
their pathway positions (s₁, s₂) and orders, the diagonal pathways, the capture-A_y terms, and for each external row
its final state and channel-capture row. The table has no derivative (nothing in it depends on θ), and the
cotangent buffers stay per evaluation and per thread as now. `TBarProducts(Decay *)` and the cached pointer
`tbar_decay_` (and `GenMatrixFunc::rt_decay_`) key a per-object scratch by the decay's address; they become a pointer
to the `ExitRow` (or the pair (entrance, exit)), with the same invalidation (`AddTMatrixElement` resets
`rt_decay_`).

**What the layout must not make harder: the coupled spin–orbit treatment** (`spin_orbit_potential.md` §2.3, option
(b), §5 step 6). There the channel functions of the channel spins of one (α, l, J) become a matrix block, and the
collision matrix of a pathway (c, c′) needs the block's matrices instead of products of per-channel phases and √P.
The table must therefore name a pathway by its block and block channels (it keeps `block`, `k_in`, `k_out` beside the
engine numbers), must hold no channel-function value or per-channel scalar (it holds none: they depend on the energy),
and must not assume that a channel spin identifies one channel of an (α, l, J) (the spin groups are a property of the
pathway pairs, not of the channels). Option (b) changes the arithmetic of the T-matrix in the matrix functions, not
the pathway list.

### 11.5 The other pathway shapes

- **Unobserved primary, observed secondary (UPOS).** A particle exit α′ ≠ α has spin groups (s, s′, s₂′); a group with
  s′ ≠ s₂′ holds the pathways of both exit spins. `SpinGroupRow` keeps `sp2`, `has_sp2` and `CrossSpin()`; the
  interference rows keep `z` (zero for a cross-spin group) and the new `z_upos`, and `sp_equal`. The readers that
  skip cross-spin groups (`CalculateAmplitudeMatrix`, the RMC branch and the angle-integrated sum of
  `CalculateBareCrossSection`, `PointTBar`'s guard) test `g.CrossSpin()` or `g.sp != g.sp2` exactly as they test
  the KGroup now.
- **Beta-delayed.** A beta entrance has groups (2l, 0) for l = 0, 1, internal pathways only, `stat_spin_factor` 0.0,
  and no interference lists (`CalcAngularDists` skips a pair of type 20). Its exits have `has_interferences = false`,
  which is `TablesBuilt() == false` today; the readers' fall-back branches for that case walk empty `KLGroup` lists
  and contribute nothing (§11.9).
- **Channel capture.** An external row continues an internal pathway of another exit of the same entrance; the tree
  stores the position triple (decay, k, m), the table also stores `cc_row`, the internal row itself, which is all
  that `CalculateTMatrix`, `TBarAdjoint`, `RMatrixFunc::CalculateTMatrix`, `ECPathwayAmplitude`, `WalkECAmplitudes`,
  `SetUpFormalCaptureMasks` and `GetECIntegrals` use it for (group, channel, exit channel, and the exit channel's l
  and s). The triple stays for the check file. `photon_channel` is the J-group channel number the tree calls
  `GetIntChannelNum`, 0 when the block has no photon channel to that final state.
- **Capture A_y.** Today the terms are built lazily under a mutex (§11.2). In the table they are built at set-up for
  each (entrance, exit) that a capture-A_y segment uses (an analyzing-power segment with a photon exit: the
  condition under which `GenMatrixFunc::CalculateCrossSection` calls `CalculateCaptureAnalyzingPower`,
  `GenMatrixFunc.cpp` 573–581), at `maxLOrder`, with the same early returns (photon exit only; projectile spin > 0).
  The set-up of the run knows these segments before `CNuc::Initialize` (the data are read first). A capture-A_y point
  whose (entrance, exit) was not requested must stop with a message, not compute A_y = 0 silently. Building the terms
  for every photon exit instead would remove the request list at the cost of 9j symbols for every capture project with
  a spinning projectile (13N, T1, …); that cost was not measured (§11.10).

### 11.6 Steps

Each step is one branch, mergeable alone, bit-identical: 183 projects × 3 bases against `tab-trial8` (704 / 677 /
704 identical), `ctest -LE slow`, `tests/run_tests.sh`, and calculations at 1 against several threads on T1, T2, T3,
13N and the projects named in the step, plus the measurement of §11.8. Counts are accessor occurrences as in §11.1,
before → after, in the files the step touches.

1. **The table beside the trees.** `include/PathwayTable.h`, `src/PathwayTable.cpp` (`BuildPathwayTable`), the handle
   on `CNuc` and its build in `CNuc::Initialize`; the trees are still built and read. No reader changes (accesses
   406 → 406). Test (new ctest, in `tests/upgrade/pathways_from_cnuc_test.cpp`, which already loads every `.azr` of
   `tests/` in the Brune and formal bases with every particle and beta pair an entrance): table against trees, every
   field and every position — exits, groups, internal and external rows (with `cc_row` resolving to the tree's
   back-reference), each interference (k, m₁, m₂, L, types, Z̄₁Z̄₂, Z̄₁Z̄₂ UPOS) in flat order, interference pathways,
   orders, diagonal pathways, and the capture-A_y terms of every photon decay (forcing the tree's lazy build with
   `CalcCaptureAnalyzingPower` at `maxLOrder`); doubles compared bit for bit. Plus a unit test of properties that do
   not reuse the builder: every row's (group, ch, chp) is a channel pair of one J-group with the group's spins; every
   interference names rows that exist; for a photon exit, a_k = 4 Z̄₁Z̄₂ for the pathway pairs of one spin group (the
   relation `CNuc.cpp` 1175–1176 and `GenMatrixFunc.cpp` 895–897 state; if it holds only summed over the pairs of one
   order, the test compares the sums). Outputs identical by construction (nothing reads the table). Memory: one table
   per run.
2. **Check files and set-up readers** (disjoint files, may run as three parallel branches after step 1):
   (2a) `WritePathwaysCheck`, `WriteAngularDistributionsCheck` (`CNucReport.cpp` 20 → 0); (2b)
   `CNuc::SetUpFormalCaptureMasks` (`CNucExternalCapture.cpp` 10 → 0); (2c) `AZUREAPI::GetECIntegrals`
   (`AZUREAPI.cpp` 9 → 0). Test: the four check files of T1, T3, 13N and `beam_profile_kernel` written with the checks
   on are byte-identical to the base's (as in the record of steps (d)); for 2c, `tests/ec_integrals_signature/check.sh`
   and `tests/ec_relativistic_signature/check.sh` (they reach `GetECIntegrals` through pyazr) and the arrays of
   `pyazr/examples/ec_integrals.py` on 13N compared number for number, base against step.
3. **External-capture amplitudes at set-up**: `WalkECAmplitudes`, `EData::CalculateECAmplitudes`,
   `EData::InitializeComponentSegments` (`EData.cpp` 22 → 0), `EPoint::CalculateECAmplitudes`,
   `EPoint::GetECAmplitudeWithShift`, `ECPathwayAmplitude` (`EPoint.cpp` 12 → 0). After the file split of
   `edata_split.md` is merged (these functions move to `EDataExternalCapture.cpp` and `EPointExternalCapture.cpp`, in
   work). The exit is looked up as the tree's numbering does today (`GetDecay(exit pair number)`, §11.9 item 1): the
   table's `ExitOf` gives the same position in every probed project, and the fix is a separate commit. Test: the 21
   `intEC.dat` files of the Brune basis identical, and the signature (`EData::ECSignature`, `EData.cpp` 2117, which
   hashes the pathways in tree order through `WalkECAmplitudes`) unchanged, so that a file written by the base is
   accepted; 13N and `beam_profile_kernel` byte-identical at 1 and 8 threads.
4. **The R-matrix route**: `RMatrixFunc::CalculateTMatrix` (`RMatrixFunc.cpp` 16 → 0). Independent of 3 and 5.
   The 183 × 3 comparison does not reach it: no `.azr` under `tests/` has `false` on the first line of `<config>`
   (counted: 0), and `Config::ResolveMatrixRoute` switches to the A-matrix whenever ignore-externals (the default) or
   Park's basis is on. The step therefore adds its own comparison, base against step, byte for byte: T1, T3, 13N and
   `beam_profile_kernel` copied with that line `false`, run with `--no-ignore-externals` in the Brune and formal
   bases.
5. **Forward hot path**: (5a) `AMatrixFunc::CalculateTMatrix` (`AMatrixFunc.cpp` 111 → 95); (5b)
   `GenMatrixFunc::CalculateBareCrossSection`, `CalculateAmplitudeMatrix`, `BuildCoulombInterference` and the
   `rt_decay_` key (`GenMatrixFunc.cpp` 119 → 6, `GenMatrixFunc.h` `Decay *` → `ExitRow *`). 5a and 5b touch disjoint
   files and may run in parallel. The R-matrix-capture branches (`--use-rmc`: `CalculateTMatrix`'s `isRMC`, the
   first branch of `CalculateBareCrossSection`) are not reached by the 183 × 3 comparison either (`--use-rmc` occurs
   only in `tests/thm_lineshape/check.sh`): the step compares T1 and 13N with `--use-rmc`, base against step. The
   fall-back branches for `TablesBuilt() == false` either go (they add nothing) in a separate commit with its own
   comparison, or are kept reading the empty list. Test: as for every step, and
   T1–T3, `upos_pair_index`, `light_excitation`, `beta_delayed_alpha`, `ec_e1m1_phase`, `identical_entrance_reaction`
   at 1 and 8 threads.
6. **Capture A_y at set-up**: the request list from the set-up of the run;
   `GenMatrixFunc::CalculateCaptureAnalyzingPower` and the capture-A_y block of `PointTBar` read the table;
   `CNuc::CalcCaptureAnalyzingPower` (78 lines and its 24-line comment), `captureAyMutex`, `Decay::captureAy_`,
   `captureAyBuilt_` and their accessors are deleted (`GenMatrixFunc.cpp` 6 → 0, `AMatrixFunc.cpp` 95 → 92). After 5b
   (same file). Test: `13N_capture_ay`, `capture_ay_coupling` byte-identical at 1 and 8 threads; a ctest that a
   capture-A_y point of an unrequested reaction stops with its message.
7. **Gradient**: `PointTBar` (rest), `TBarAdjoint`, `CollectBilinearTerms`, `TBarProducts`, `tbar_decay_`
   (`AMatrixFunc.cpp` 92 → 0, `AMatrixFunc.h` 1 → 0). After 5a and 6. Test: the gradient tests of ctest (analytic
   against finite differences), the T2 fit (200 calls) at 1 thread identical in every output file, at 8 threads as the
   base (the base's own 1-vs-8-thread fit differs in the 8th digit of χ², `channel_table.md` §6).
8. **Delete the trees.** `SortPathways` (85 lines) and `CalcAngularDists` (144) go (`CNuc::Initialize` calls
   `BuildPathwayTable` only); `PPair::decays_`, `AddDecay`, `ClearDecays`, `NumDecays`, `IsDecay` (2), `GetDecay` (2);
   the class files `Decay` (117 + 206 lines), `KGroup` (57 + 124), `MGroup` (35 + 52), `ECMGroup` (70 + 156),
   `KLGroup` (47 + 67), `Interference` (42 + 67) — 1040 lines — and `CaptureAyTerm.h` (43, its struct moved in step 1);
   `pathways_from_cnuc_test.cpp` (383 lines, ctest `physics_pathways_from_cnuc`: both its comparisons need the trees);
   the six sources from `src/CMakeLists.txt`. `CNuc.cpp` loses lines 886–1163 except `CalcBoundaryConditions`
   (`SortPathways` with its comment, 886–974, and `CalcAngularDists` with its comment, 1016–1163: 237 lines);
   `CalcCaptureAnalyzingPower` went in step 6. The property tests of step 1 and `physics_pathways_unit` stay. Test:
   the full set of the step rule.

Order: 1, then 2a, 2b, 2c, 3, 4, 5a, 5b in any order (3 after the file split), then 6 (after 5b), 7 (after 5a and 6),
8 last. Steps 2a–2c, 3, 4, 5a and 5b touch disjoint files.

### 11.7 Why the steps are bit-identical

Each step replaces where an integer or a double comes from (a node's field → a row's field holding the same value),
never an arithmetic expression or the order of a sum. The only places that compute anything are the builder's
coefficient loops, which are `CalcAngularDists`' and `CalcCaptureAnalyzingPower`'s loops moved with their arguments
taken from the System; step 1's test compares their doubles bit for bit with the trees' on every project before any
reader switches. A difference in a later step is therefore a reader that iterates or indexes differently, and the
step stops until it is found.

### 11.8 Cost

**Indirections per inner iteration.** Today an internal pathway of `CalculateTMatrix` costs, in the loop condition and
body, the calls `Decay::GetKGroup` (twice), `KGroup::NumMGroups`, `KGroup::GetMGroup`, `MGroup::GetJNum`, `GetChNum`,
`GetChpNum`, `CNuc::GetJGroup` and `JGroup::GetChannel` (twice): eleven out-of-line calls (the getters are defined in
`.cpp` files; the build is `-O3` without link-time optimisation, `CMakeCache.txt` of the base build), each one or two
dependent loads. On the table it is one indexed load of the row (the bases of `P.internal` and `g` are loop-invariant)
and three field reads, no call. A channel-capture row adds today `PPair::GetDecay`, `Decay::GetKGroup`,
`KGroup::GetMGroup` and three getters; on the table one indexed load. `CalculateBareCrossSection`'s angle-integrated
loop reaches `CNuc::GetJGroup` and `JGroup::GetChannel` for l, l′ and J of every pathway; those stay in 6b (they are
channel data, the other half of step 6), so its gain is smaller. Indexing the table cannot be slower than chasing the
trees unless the compiler fails to hoist the row bases; inlining changes when code moves between files (the file split
stopped once on the T2 calculation's inlining, in work), which is why each step is measured.

**Measured on the base build** (`work/b-trial9`, this host, one thread, `--no-coulomb-table`):
- T2 calculation (menu 1), `perf stat -e cycles:u,instructions:u`, five runs: 3.308–3.317 × 10⁹ cycles (spread
  0.3 %), 11 688 402 464–11 688 402 966 instructions (spread 4 × 10⁻⁸).
- `perf record` of the T2 calculation (1060 samples): the pathway readers are below 1 % (`CalculateTMatrix` 0.37 %,
  `CalculateBareCrossSection` 0.28 % self); 94 % of the run is the Coulomb functions of the set-up.
- `perf record` of a T2 fit of 30 calls (`--use-lm --max-calls=30`, 22 133 samples): `CalculateBareCrossSection`
  1.39 % self (2.1 % with callees), `CalculateTMatrix` 0.58 %, `JGroup::GetChannel` 0.37 %, `CNuc::GetJGroup` 0.24 %,
  `Decay::GetKGroup` 0.09 %, `KGroup::GetMGroup` 0.04 %, `PPair::GetDecay` 0.02 % (self). The readers and the getters
  they call are 3–4 % of a short fit (`GetChannel` and `GetJGroup` also serve the level matrix); a step that only
  changes how rows are reached moves a fraction of that, so a 2 % change of the fit's cycles means it added work.

**What each step measures** (base = the head before the step, same host, alternating base and step, `OMP_NUM_THREADS=1`):

    perf stat -r 5 -e cycles:u,instructions:u <AZURE2> --no-gui --no-readline --no-coulomb-table T2.azr < menu
    # T2 calculation: menu '1\n\n\n7\n' in a copy of tests/upgrade/T2 (as scripts/bench_upgrade.sh)
    perf stat -r 3 -e cycles:u,instructions:u <AZURE2> --no-gui --no-readline --no-coulomb-table \
        --use-lm --max-calls=200 T2.azr < menu
    # T2 fit of 200 calls: menu '2\nn\n\n\n7\n'

and, for the steps that touch them, 13N (external capture at set-up: step 3), `13N_capture_ay` (step 6) and the T2
fit with the analytic gradient path (step 7). **Stop rule:** the mean cycles of the step exceed the base's by more
than 2 % on the T2 calculation or on the 200-call fit: stop, profile (`perf record -g`), and report. Instructions are
reproducible to 10⁻⁷ and show any change in the work done; a rise in instructions with cycles inside 2 % is reported,
not a stop. Peak RSS (`/usr/bin/time -f %M`) of the fit at 1 and 8 threads is reported; it is expected to fall, since
the working copies share one table instead of copying the trees.

### 11.9 What 6b does not solve, and what was met while reading

**What remains.** After 6b `CNuc` still holds the pairs, J-groups, channels and levels (the `GetPair`/`GetJGroup`
half of step 6: about 730 calls, §10), the boundary values, the external-capture set-up, the parameter vector, the
abundances and penalties (step 7), and the handle to the table. The readers still call `CNuc::GetJGroup` and
`JGroup::GetChannel` for channel data (l, J, E/M of a pathway's channels) and `EPoint`'s per-(J-group, channel)
tables; the engine numbers in the rows exist only for these.

**What Phases 5 and 6 will want differently.**
- One collision-matrix solve per grid energy for all reactions (Phase 5): the T-matrix of a pathway becomes a slot of
  a flat per-energy buffer indexed by the row number of the table, not `tmatrix_[decay][k][m]`; the (k, m) contract
  of §11.1 (and the intEC line order, the cotangent buffers, the EC cache key) then changes in one commit, with the
  intEC signature bumped.
- Observables as functions of U with "coupling coefficients precomputed per reaction" (ARCHITECTURE §3): the flat
  interference and capture-A_y lists are those coefficients; they move to physics/ with System indices when the
  engine numbers go.
- Evaluation without copies (Phase 6): the table is already shared and immutable; the per-object scratches keyed by
  the last decay (`rt_decay_`, `tbar_decay_`) become per-task scratch.
- The total-capture and SUM components loop over several exits of one entrance (the R-matrix-capture branch):
  Phase 5's reaction-indexed observables would rather index exits by (entrance, exit) than by decay position.

**Defects and oddities (not fixed; file and line at cb8f5912).**
1. *An exit pair number used as a decay position.* `EData.cpp` 2011, 2394–2395, 2680–2681; `EPoint.cpp` 611,
   1765–1766; `AZUREAPI.cpp` 1751–1752 call `entrancePair->GetDecay(ir)` with `ir` the exit pair's number. The decays
   of an entrance are in pair order but skip a beta pair and any pair without a pathway from that entrance (another
   compound nucleus, an unconnected pair), so whenever such a pair is numbered below the photon pair these lines read
   another decay or index past the end. The other readers search the decay by `GetPairNum()`. Probed: 14 of the 21
   projects of `tests/` with external capture (one `.azr` of each folder) were run with the pathways check file on;
   in every one each decay's position equals its pair number, so `tests/` does not exercise the defect. A test project
   with a beta pair or a second nucleus listed before the photon pair is needed.
2. *The not-found check of the decay search is unreachable.* The search `while (ir < N) { ir++; if (match) break; }`
   (`GenMatrixFunc.cpp` 49, 814, 857; `AMatrixFunc.cpp` 533, 703, 1147, 1859) ends with `ir = N` when no decay
   matches, so the following `if (ir > NumDecays()) return false` (`GenMatrixFunc.cpp` 818, 861) never fires and the
   last decay is used; with no decay at all `GetDecay(0)` indexes −1. `ExitOf` returns 0 for "no
   exit", and the readers must then refuse instead of reading.
3. *Capture-A_y terms built in the point pass* under a function-level static mutex, for the first caller's Legendre
   order (§11.2); removed by step 6.
4. *`AZUREAPI::GetECIntegrals` uses the final state's formal energy* (`ecLevel->GetE()`, `AZUREAPI.cpp` 1786); the
   engine's integrals use the Brune energy `GetECEnergy()` since step 2 (`EPoint.cpp` 615, 1769). In the formal basis
   the API's integrals are at another energy than the engine's.
5. *Dead fields.* `MGroup::tmatrix_`, `ECMGroup::tmatrix_` (never initialised or read), `ECMGroup::statspinfactor_`
   (never set; its getter has no caller). The three unused pointers `theJGroup`, `entranceChannel`, `exitChannel` in
   `AMatrixFunc::CalculateTMatrix` (`AMatrixFunc.cpp` 549–551) cost three calls per pathway; they go with step 5a.
6. *Fall-back branches for `TablesBuilt() == false`* (`GenMatrixFunc.cpp` 287–316 and the Coulomb-interference
   branch 936–948): reached only by a decay without interference lists (beta entrance), where the `KLGroup` list they
   walk is empty.

### 11.10 What this design rests on

Counted (scripts and outputs in `work/pathnote/`): the accessor counts of §11.1 (`inventory.py`, `grep -o`); the
writers of §11.2 (`grep` of every `Add…`, `Set…`, `Clear…` of the tree classes over `src/`, `api/`, `include/`,
`tests/`); the line counts of §11.6 (`wc -l`); the decay numbering of 14 external-capture projects (`ecprobe.sh`,
base build, one thread). Measured: the `perf stat` and `perf record` figures of §11.8 (base build `work/b-trial9`, this
host, one thread). Inferred by reading, not run: that every row reproduces the tree (step 1's test establishes it);
that every capture-A_y caller asks for `maxLOrder` (§11.2); that `channelOfBlockChannel` is the identity
(`CNucSystem.cpp`); the indirection counts of §11.8; the cost of building capture-A_y terms for every photon exit
(not measured).

### 11.11 Decisions of the director (2026-11-10)

1. The flat engine table of §11.4 is the layout. Its handle (`shared_ptr<const PathwayTable>`) sits on `CNuc`
   while the readers move and goes to the object that replaces `CNuc` with step 7.
2. Capture A_y: the terms are built at set-up for the exits that the data sets with capture-A_y points name. The
   set-up collects them from the segments; a user writes nothing. A capture-A_y point whose exit was not collected
   stops with a message that names the data set (it cannot happen through the set-up; it guards the API). The static
   mutex goes.
3. Defect 1 (an exit pair number used as a decay position) is fixed before the external-capture readers move
   (step 3), in its own commit with a test project that has a pair without a pathway numbered below the photon
   pair; the outputs of every present project stay identical.
4. Defect 2 (the unreachable not-found check) is fixed with the reader it belongs to: no exit refuses with a
   message, it never reads the last decay.
5. The fall-back branches for `TablesBuilt() == false` are deleted in a commit of their own, after a run of every
   project with an abort placed in them shows that they are not reached with a non-empty list.
6. The step-1 comparison of table and trees goes into `pathways_from_cnuc_test.cpp`, which step 8 deletes.
7. Defect 4 (`AZUREAPI::GetECIntegrals` at the formal energy) is fixed with step 2c and gets a pyazr test in the
   formal basis.
8. The steps start after the cut of `EData`, `EPoint` and `ESegment` into smaller files is merged. The R-matrix
   route and `--use-rmc` are compared separately at the steps that touch them (§11.6), since no project of
   `tests/` selects them.

### 11.12 Status (2026-11-11, evening)

Done: steps 1, 2, 4, 5 and 7 (`feature/upgrade-4-pathway-table-1`), step 6 (`feature/upgrade-4-capture-ay`), steps 3
and 8 (`feature/upgrade-4-pathway-table-2`, below).
Capture A_y at set-up: `RequestCaptureAyTerms` (`src/CaptureAyRequests.cpp`) names, after the data are filled and
before `CNuc::Initialize`, every term (own and leaf) of an analyzing-power data set or extrapolation with a photon
exit; it runs in the four set-up paths (`AZUREMain`, `runMCMC`, `AZUREAPI`'s session set-up and `RebuildImpl`).
`BuildPathwayTable` builds the requested terms in the table's order (independent of the order of the requests) and
none for a beta entrance.  A capture-A_y point without its terms is refused with a message naming the data set.

Step 3 and step 8 (`feature/upgrade-4-pathway-table-2`, 2026-11-11), with decision 3 first:
- Defect 1 is fixed (its own commit): the external-capture code finds the exit to the photon pair by the pair's number
  (`ExternalCaptureExit`, `include/ExternalCaptureExit.h`, in `src/ExternalCaptureSetup.cpp`), and an entrance without
  a pathway to it is refused at the signature walk, naming the data segment and both pair keys.  Test project
  `tests/ec_exit_position` (13N with a 13C level whose n + 12C pair is numbered below or above the photon pair): the
  base wrote an empty `intEC.dat` and died with SIGSEGV on the first order; now both orders give the same outputs byte
  for byte; ctest `upgrade_external_capture_exit` checks the position and the refusal.
- The memo of the shifted amplitude (`Unit::GetECAmplitudeWithShift`) recorded the final level's energy at the first
  shifted evaluation of a copy.  That energy cannot depend on the copy's history: `ALevel::GetECEnergy` is written only
  by `CNuc::SetUpECFinalBlocks` in `CNuc::Initialize`.  The memo keeps (energy, amplitude) and reads the set-up energy;
  `history_test` has a project `ec` (external capture, a fitted shift, the final state's energy free) and a sequence
  B -> C against fresh C on every project (shown failing with the energy recorded per copy from the fit energy).
- Step 3: the walk of the intEC count and signature, the amplitude loops of the data segments and of the leaves (the
  latter moved verbatim from `WorkingCopy.cpp` to `WorkingCopy::LeafECAmplitudes`), `Unit::CalculateECAmplitudes`,
  `Unit::GetECAmplitudeWithShift` and `ECPathwayAmplitude` (`src/UnitExternalCapture.cpp`) read the table.  The
  signature is unchanged: files written by the base and by the table build are identical and each reads the other's.
- Step 8: one commit ran every project in the three bases with the trees not built and their accessors aborting
  (identical to `tab-trial18`); the next deleted `SortPathways`, `CalcAngularDists`, `PPair::decays_` and its
  accessors, the classes `Decay`, `KGroup`, `MGroup`, `ECMGroup`, `KLGroup`, `Interference` (12 files, 971 lines) and
  `tests/upgrade/pathways_from_cnuc_test.cpp` (531 lines).  `CNuc.cpp`: 1,774 -> 1,548 lines.  The pathways check
  file is written from the table (byte-identical).
- The double final state (two photon pairs with one final state): refused at set-up (`FinalStateOfAnotherPair`);
  `tests/ec_final_state_ambiguous` (d) and (e) are exactly the case met in step 6 (13N with a copy of its photon
  pair), so no test was added.

The pathways now exist once, in the table.  What remains of step 6 of §10 is the other half of `CNuc`'s readers:
`GetPair`/`GetJGroup` for channel data in the matrix functions and the per-(J-group, channel) tables of the points
(§11.9 "What remains"); step 7 (the parameter vector, the abundances and penalties, the set-up from the level rows,
the handle to the table) is open.

Met while doing step 6:
- The trees built capture-A_y terms for a beta entrance; both readers return before the terms for a non-particle
  entrance, so no output read them (every output identical without them).
- A photon pair names one final state (its J^pi group and the level within 1 keV of its E_x; a second level in the
  window is refused, `tests/ec_final_state_ambiguous`), so two final states cannot share one photon pair, and the
  loops of `GetECIntegrals` and `EPoint::CalculateECAmplitudes` over the EC levels of an exit's photon pair see one
  level.  The reverse loads: two photon pairs with the same final state.  `ParseExternalCapture` then sets the level's
  capture pair twice (`ALevel::SetECParams`; the last wins) and appends its NF integrals and conversion factors twice,
  so capture into the first pair runs without external capture and without a message (13N with a copy of its photon
  pair: 1.95e-9 against 4.32e-9 at the first point, fit column of AZUREOut; `GetECIntegrals` returns the second
  pair's 8 pathways, as the engine computes them).  The defect is in `ParseExternalCapture`, not in `GetECIntegrals`; not fixed.

## 12. The end of `CNuc`: design (2026-11-13)

This section lists again what `CNuc` still holds, measured on the code at 6e8d000b (the data head; the engine sources
are those of `work/b-trial50`), and designs the last steps so that agents can carry them out without a second design
round. It replaces "Remaining steps" (section 6), §11.9 "What remains" and the remaining items of §11.12; where those
and the code disagree, the code is taken (§12.1.5). No code was changed for this section.

### 12.1 Inventory, measured

Scripts and outputs in `work/cnucdesign/` (run from the worktree root):
`defs.py CNuc src/CNuc*.cpp src/*.cpp` (every `CNuc::` definition with its length: `cnuc_defs.tsv`);
`members.py . <own files> <names>` (calls `.M(`/`->M(` of each member per file outside the class's files:
`cnuc_callers.tsv`, `alevel_members.tsv`); `naming.py`, `group.py` (the files naming `CNuc`: `naming.tsv`,
`groups.tsv`); `clone_probe.cpp` (bytes and time of a copy, below); `grep` counts as quoted. Counts are name searches:
a member name another class also has (`GetNucleus`, `Fill`, `Initialize`, `FillMnParams`, `GetE`) is over-counted (≈).

**Size.** `src/CNuc.cpp` 1,527 lines (1,237 in definitions), `include/CNuc.h` 327; with `src/CNuc*.cpp` (seven files,
1,463 lines) and their headers, 3,647 lines. About 69 member functions (62 defined in `.cpp` files, the rest in the
header) and 17 data members.

#### 12.1.1 The members of `CNuc`

Classes: **(a)** immutable after set-up; **(a′)** set-up input that a re-set-up rewrites (a physical vector handed in
again: `FillCompoundFromParamsPhysical` + `TransformIn`, per evaluation in the sampler's physical mode); **(b)** the
fit's working state (changes with the parameter vector); **(c)** scratch of one evaluation or one report; **(d)** dead
or reachable only inside the class. "Ext." is calls outside `src/CNuc*.cpp` / files.

- `pairs_`, `NumPairs`, `GetPair` ×2, `GetPairNumFromKey`, `GetPairFromKey`, `IsPairKey` (3+4+3+12+4+9) — written when:
  `Fill`; `PPair::SetEntrance` by the data set-up (`WorkingCopy.cpp` ×4, `AZUREMain.cpp`, `AZUREAPI.cpp`) after `Fill`;
  ext. calls / files (main readers): 18/9, 190/40, 114/31, 8/6, 34/14 (WorkingCopy 76, matrix functions 108); class: (a)
  but for the entrance flag; goes to: `CompoundSystem`; the entrance set an input of the build (step 1).
- `jgroups_`, `NumJGroups`, `GetJGroup` ×2, `IsJGroup` (3+4+3+16) — written when: `Fill`, `AddHardSphereGroups`, B_c in
  `Initialize`; ext. calls / files (main readers): 60/23, 192/26 (RMatrixFunc 41, AMatrixFunc 32, GenMatrixFunc 18,
  GammaDecayMap 15), 1/1; class: (a), holding the (b) levels; goes to: `CompoundSystem` (groups, channels, level set-up
  data); the levels' (b) part to `LevelState`.
- `nuclei_`, `NumNuclei`, `GetNucleus` ×2, `NucleusLabel`, `FindOrAddNucleus` (3+3+3+5+9) — written when: `Fill`; ext.
  calls / files (main readers): 4/2, ≈42/12, 8/2; class: (a); goes to: `CompoundSystem`.
- `CapturePairKeys`, `RatePairsProblem` (9, 15) — ext. calls / files (main readers): 2/1, 1/1; class: (a) queries; goes
  to: `CompoundSystem`.
- `IsPair`, `AddPair`, `AddJGroup`, `AddHardSphereGroups`, `ParseExternalCapture`, `Fill` (12, 4, 3, 48, 71, 215) —
  written when: set-up; ext. calls / files (main readers): 0, 0, 0, 0, 0, ≈12/6; class: set-up; goes to:
  `BuildCompoundSystem` (step 5); `IsPair`/`Add*` private.
- `Initialize`, `CalcBoundaryConditions`, `BuildPathways`, `RequestCaptureAy`, `captureAyRequests_`, `pathways_`,
  `Pathways` (35, 35, 10, 6, —, —, 4) — written when: set-up; ext. calls / files (main readers): ≈15/8, 0, 0, 1/1, —, —,
  14/7; class: (a); goes to: `BuildCompoundSystem`; the table handle on `CompoundSystem` (§11.11 decision 1).
- `SetUpECFinalBlocks`, `ecFinalBlocks_`, `SetUpFormalCaptureMasks`, `captureMaskSetup_` (53, 34) — written when:
  `Initialize`; ext. calls / files (main readers): 0; class: (a), shared by copies; goes to: `CompoundSystem`
  (external-capture set-up).
- `ECFinalStatesFromFormal`, `FormalCaptureMasksFromParams`, `captureMasks_`, `FormalCaptureMask`,
  `FormalCaptureMasksOutsideAdjoint`, `MaskedLevelCoupled`, `HasFormalECFinalBlocks` (25, 34, —, 5, 12, 5, 1) — written
  when: per vector (`FillCompoundFromParams`); ext. calls / files (main readers): 0, 0, —, 1/1, 1/1, 0, 1/1; class: (b);
  goes to: `LevelState` (step 3).
- `TransformIn`, `transformedIn_`, `IsTransformedIn`, `transformError_`, `TransformError` (117, 1, 1) — written when:
  `Initialize`; again after a physical vector; ext. calls / files (main readers): 6/3 (MCMC 3, API 2, limits 1), 2/1, 0;
  class: (a′); the error (c); goes to: `LevelInputs` (a value; step 4).
- `CalcShiftFunctions`, `InvalidateShiftFunctions` (168, 10) — written when: per vector; ext. calls / files (main
  readers): 32/10, 1/1; class: (b); goes to: `LevelState` (step 3).
- `FillCompoundFromParams`, `SetFitLevelParameter` (24, 16) — written when: per vector; ext. calls / files (main
  readers): 35/8, 3/1 (ModelData); class: (b); goes to: `FillLevelState` (steps 3–4).
- `FillCompoundFromParamsPhysical` (22) — written when: re-set-up; ext. calls / files (main readers): 8/5; class: (a′);
  goes to: `LevelInputs` from an observed vector (step 4).
- `FillMnParams`, `ConvertAmplitudeBasis` (65, 46) — written when: creates θ; ext. calls / files (main readers): ≈36/8,
  1/1; class: parameter vector; goes to: `LevelParameterLayout` (step 4).
- `TransformOut`, `TransformToBrune`, `ObservedWidths`, `ObservedWidthsAt`, `ObservedOfFitValues`, `GetTransformParams`
  (4, 62, 42, 12, 21, 19) — written when: per report; `ObservedOfFitValues` inside an evaluation (model data); ext.
  calls / files (main readers): 13/8, —, —, 3/2, 1/1, 10/6; class: (c), writing the levels' report slots; goes to:
  `ObservedLevelTable`, a value returned (step 3).
- `ParkNormPenalty`, `kParkNormScale` (14) — written when: per vector; ext. calls / files (main readers): 3/2
  (Objective, MCMC); class: (b) reader; goes to: a term of `Objective` reading `LevelState` (step 2).
- `samples_`, `hasAbundanceBlock_`, `SetAbundances`, `ReadAbundances` and 14 abundance accessors (43, 3, ≈60) — written
  when: set-up; the values per vector; ext. calls / files (main readers): 1–8 each in 6 files (Objective, AZUREGrad,
  ParameterLabel, AZUREAPI, ModelPass, WorkingCopy); class: spec (a), values (b); goes to: spec in `DataSet`, values in
  the working copy, penalties as terms (step 2).
- `AbundancePenalty`, `AbundanceSumPenalty` (13, 9) — written when: per vector; ext. calls / files (main readers): 1/1,
  2/1; class: (b) reader; goes to: `Objective` terms (step 2).
- `modelData_`, `GetModelData` — written when: `Fill`; ext. calls / files (main readers): 10/5; class: (a), shared; goes
  to: `CompoundSystem`.
- `gammaDecays_`, `GetGammaDecays` — written when: `Fill`; ext. calls / files (main readers): 9/3; class: (a), shared;
  goes to: `CompoundSystem`.
- `gammaDecayFills_`, `GammaDecayFillOf`, `GammaDecayValues`, `HasGammaDecayValues`, `SetGammaDecayFill` — written when:
  per vector (`GammaDecaySet::ToEngine`); ext. calls / files (main readers): 10/1 (GammaDecayMap); class: (b); goes to:
  `LevelState`.
- `gammaDecayInputs_` — written when: `Fill`, `FillCompoundFromParamsPhysical`; ext. calls / files (main readers): 0;
  class: (a′); goes to: `LevelInputs`.
- `thmCoherentValues_`, `ThmCoherentValues`, `SetThmCoherentValues` — written when: per vector, by the data
  (`ThmNorm.cpp` ×2); ext. calls / files (main readers): 2/1, 2/1; class: (b), written by another object; goes to: the
  working copy, which owns the parameters (step 2).
- `Clone` (4) — ext. calls / files (main readers): 17/5 (AZURECalc 8, MCMC 4, ReactionRate 3); goes to: deleted (step
  4).
- `GammaDecayInputs()` (getter), `IsPair`, `AddPair`, `AddJGroup`, `TransformError`, `AbundanceParamName`,
  `BuildPathways`, `CalcBoundaryConditions` — ext. calls / files (main readers): 0 outside; class: (d) or internal only;
  goes to: private or deleted (step 1).

By class, of the 1,237 lines in definitions of `CNuc.cpp`: set-up and structure (a) ≈ 590 (`Fill` 215,
`ParseExternalCapture` 71, `AddHardSphereGroups` 48, accessors ≈ 110, `Initialize` and B_c 70, abundances' set-up 46,
pathways 30); (a′) `TransformIn` and the physical refill 139; (b) ≈ 270 (`CalcShiftFunctions` 168, fills 40, Park and
abundance penalties 36, abundance values ≈ 25); (c) reports 140; the parameter vector 111 (`FillMnParams`,
`ConvertAmplitudeBasis`). Another 168 lines of (b) are in `src/CNucExternalCapture.cpp` and 37 in `CNucObserved.cpp`.

#### 12.1.2 The classes that hang from it

**`ALevel`** (`include/ALevel.h` 217, `src/ALevel.cpp` 504; 544 bytes plus 19 vectors). Three parallel sets of values
and the set-up data, by writer (calls outside `ALevel.cpp`, `alevel_members.tsv`):

- `ids_`, `isinrmatrix_`, `energyfixed_`, `channelfixed_`, `gamma_is_rwa_`, `isECLevel_`, `ecPairNum_`, `ecMultMask_`,
  `nf_integrals_`, `ec_conv_factors_`, `ec_energy_` — written by: `Fill`, `ParseExternalCapture`, `SetUpECFinalBlocks`;
  read by (calls/files): `IsInRMatrix` 87/20, EC data 30/8; class: (a).
- `zero_input_`, `external_gammas_`, `level_e_`, `gammas_` (the input set) — written by: `Fill`, `TransformIn`,
  `FillCompoundFromParamsPhysical`, `GammaDecayMap` (one `SetGamma`); read by (calls/files): `GetE` ≈33/13, `GetGamma`
  12/6, `ZeroInputWidth` 10/4; class: (a′).
- `fitlevel_e_`, `fitgammas_`, `ec_gammas_`, `sqrt_nf_factor_`, `shifts_`, `shift_cache_energy_`, `park_norm_`,
  `shift_derivs_`, `shift_derivs2_`, `shift_d2_energy_`, `stamp_` (the fit set) — written by: `FillCompoundFromParams`,
  `CalcShiftFunctions`, `ECFinalStatesFromFormal`, `SetFitLevelParameter`, and three copies outside (`THMMatrixFunc.cpp`
  446–463 on a copy of the level, `TargetSetup.cpp` 199–200, `ReactionRate.cpp` 112–114); read by (calls/files):
  `GetFitGamma` 57/14, `GetFitE` 41/15, `GetShiftFunction` 13/5, `GetParkNorm` 13/8, `GetShiftDerivative` 10/5,
  `GetSqrtNFFactor` 5/4, `GetECGamma` 4/2, `Stamp` 3/1 (`ChannelSpace` keys its memo on the level's address and stamp);
  class: (b).
- `transform_e_`, `transform_gammas_`, `big_gammas_`, `observed_values_`, `observed_kinds_`, `transform_iter_` (the
  report set) — written by: `TransformToBrune`, `ObservedWidths`; read by (calls/files): `CNucReport` 12, `ModelData` 5,
  `AdaptiveExtrapolation` 4, `CNucChecks` 3, `GammaDecayMap` 1; class: (c).
- `NumGammas`, `SetShiftDerivatives` — read by (calls/files): 0; class: (d).

**`PPair`** (183 + 283 lines; 368 bytes): immutable after `Fill` except `entrance_` (`SetEntrance`, six writers in the
data set-up, read by `BuildPathways` and the external-capture set-up). `binding_energy_`'s setter `SetBindingEnergy` has
no caller; `ec_entrance_` is written (false) and never read: (d). It carries the run's `PairFunctionSettings` (two
strings and a pointer to the run's Coulomb-table counters) and the hybrid potential, which the System does not hold.
**`JGroup`** (81 + 136; 80 bytes): immutable after `Fill` (`SetInRMatrix`, `SetNucleusLabelled` in `Fill`).
**`AChannel`** (56 + 129; 40 bytes): immutable after `Initialize` (B_c written by `CalcBoundaryConditions`; the Wigner
limit in `Fill`). Read: `GetBoundaryCondition` 20 calls in 11 files; the channel kind is tested inline,
`GetRadType()` 62 times in 20 files and `PPair::GetPType()` 93 times in 30 files (155 tests in 37 files).

**The external-capture members.** Set-up: the final state of a photon pair (`ParseExternalCapture`: N_f integrals,
Whittaker conversion factors, `SetECParams`; `SetUpECFinalBlocks`: the Brune energy), immutable. Per vector:
`ec_gammas_` and `sqrt_nf_factor_` (`FillCompoundFromParams`; `CalcShiftFunctions` again in Park's basis with J;
`ECFinalStatesFromFormal` in the formal basis), and `captureMasks_`. Per report: `external_gammas_` written by
`ObservedWidths` (and at `TransformIn`). The adapter `SqrtNFFactor(...)` is called in five places
(`TransformIn`, `CalcShiftFunctions`, `FillCompoundFromParams`, `FillCompoundFromParamsPhysical`, `TransformToBrune`)
plus `CNucObserved.cpp` and `CNucExternalCapture.cpp`.

#### 12.1.3 The files that name `CNuc`

134 files of `src/`, `include/`, `api/` name the class at 6e8d000b (`grep -lw CNuc`; 133 when this work was
planned), and 75 files elsewhere (tests, documents, `physics/`/`io/`/`experiment/` comments). Grouped by the first of:
the class; printing; calls to the parameter vector or set-up; reads of level values; the pathway table; channel data;
nothing but the type (`group.py`):

- **the class (`CNuc*.h/.cpp`)** — files: 17.
- **type only: forward declarations, signatures that pass `CNuc *` on, comments (64 files, e.g. `Unit.h`, `Site.h`,
  `NodeTable.h`, `PointPass.cpp`, `EnergyPass.cpp`)** — files: 64; what they want: nothing of their own: a parameter
  type.
- **parameter vector and set-up, the owners (`AZURE2.cpp`, `AZUREMain.cpp`, `AZURECalc.cpp`, `AZURECalcMCMC.cpp`,
  `AZUREAPI.cpp`, `AZUREAPIModelData.cpp`, `pybind_api.cpp`, `AZUREGrad.cpp`, `Objective.cpp/.h`, `CovarianceBand.cpp`,
  `ParameterLimitsManager.cpp`, `ReactionRate.cpp`, `TargetSetup.cpp`, `AdaptiveExtrapolation.cpp`, `Predictions.cpp`,
  `ProjectSetup.cpp`, `ModelPass.cpp`, `WorkingCopy.cpp`, `GammaDecayMap.cpp`, `THMMatrixFunc.cpp`, `ThmNorm.cpp`)** —
  files: 22; what they want: `Fill`, `Initialize`, `Clone`, `FillMnParams`, `FillCompoundFromParams*`,
  `TransformIn/Out`, `CalcShiftFunctions`, abundances; and channel data.
- **channel data only (`ChannelTables.cpp`, `ChannelTable.cpp`, `NodeTable.cpp`, `Unit.cpp`, `UnitLegendre.cpp`,
  `Segment.cpp`, `ShiftFrames.cpp`, `MeshSlots.cpp`, `WorkingCopyUnits.cpp`, `PolarizationFunc.cpp`,
  `ScatteringLength.cpp`, `CoulombReference.cpp`, `GradEnergyTasks.cpp`, `AdaptiveIntegrationGrid.cpp`,
  `AZUREAPIChannelFunctions.cpp`, `pybind_potentials.cpp`)** — files: 16; what they want: pairs, groups, channels.
- **level values (`AMatrixFunc.cpp` 45 reads of the fit set, `RMatrixFunc.cpp` 24, `ModelData.cpp` 7,
  `ThmVertexBoundary.cpp` 3, `ThmLineshape.cpp` 2, `AZUREParams.cpp`)** — files: 6; what they want: the fit set of the
  levels, B_c, channel data.
- **the pathway table and channel data (`GenMatrixFunc.cpp`, `ExternalCaptureSetup.cpp`, `UnitExternalCapture.cpp`,
  `SiteObservables.cpp`, `CaptureAyRequests.cpp`)** — files: 5; what they want: `Pathways()`, pairs, groups.
- **printing (`DataCheckFiles.cpp`, `ParameterLabel.cpp/.h`, `ThmReports.cpp`)** — files: 4; what they want: channel
  data, labels.

The fit set of the levels is also read by files that do not name `CNuc`: `ChannelSpace.cpp` (keys and values),
`THMMatrixFunc.cpp` (19). In all, 907 calls of the structure accessors and 539 of `GetLevel`/`GetChannel`/`NumLevels`/
`NumChannels` outside the class.

#### 12.1.4 What `Clone` copies, and what a copy needs

`Clone` is the implicit copy: every pair (with its two strings), group, channel and level (19 vectors each), the
samples, the per-vector members; the set-up objects behind `shared_ptr` (pathway table, model data, `gamma_decay` set,
external-capture set-up) are shared. Measured with `clone_probe.cpp` (linked against `work/b-trial50`'s libraries, so
`Clone`, `FillCompoundFromParams` and `CalcShiftFunctions` are the build's; the probe loads a project as
`tests/upgrade/system_from_cnuc_test.cpp` does, Brune basis, external capture on, one thread; heap bytes counted from
the vectors' capacities; times are means over 2,000 copies, 100 for the shift functions, on this shared host):

| project | pairs/groups/levels/level-channels/θ | copy: bytes | (b) | (c) | time | refill | refill + S_c |
|---|---|---|---|---|---|---|---|
| T2 | 3 / 2 / 5 / 20 / 26 | 7.8 kB | 0.8 kB | 0.5 kB | 5.2 µs | 0.9 µs | 16.5 µs |
| 13N | 2 / 8 / 12 / 30 / 42 | 12.9 kB | 1.4 kB | 0.9 kB | 7.5 µs | 1.1 µs | 7.2 µs |
| `natural_neon` | 7 / 41 / 86 / 416 / 505 | 130.6 kB | 21.8 kB | 13.6 kB | 56 µs | 9.9 µs | 81 µs |

"Copy" is all the bytes `Clone` copies, "(b)" and "(c)" the fit set and the report set within them, "time" the time
of one `Clone`, "refill" one `FillCompoundFromParams`, "refill + S_c" a refill after `InvalidateShiftFunctions` with
every shift function recomputed (`CalcShiftFunctions` without second derivatives). `natural_neon` has the largest
`<levels>` of the tree (502 rows).

A copy needs only the (b) part, 10–17 % of the bytes; the rest is set-up the copy could share. A calculation makes
almost no copies: the program at one thread (`printf '1\n\n\n7\n' | AZURE2 --no-gui --no-readline`, under `gdb` with a
counting breakpoint) cloned 1 (T2), 0 (13N) and 1 (`natural_neon`) times, against 0.23 s, 8.0 s and 1.35 s of wall
time; `FillCompoundFromParams` ran 3 times, `CalcShiftFunctions` once. A fit keeps one copy per thread in a pool
(`cnuc_pool_`), and the sampler one per concurrent walker plus one per snapshot (every 10 steps). The end of `CNuc`
therefore gains no measurable time in a calculation or a fit; it gains a copy that cannot write the shared model
(state owned by the evaluation) and 85–90 % less memory per copy, which matters only for many sessions of a large
project (the web UI: 130 kB per session and per pooled copy for `natural_neon`, 22 kB after).

#### 12.1.5 Where the notes and the code disagree (the code is taken)

1. `ROADMAP.md` ("`CNuc`, to finish") gives `CNuc.cpp` 1,774 lines (1,527) and lists as remaining the external-capture
   amplitudes at set-up with the decay-position defect, and the deletion of the pathway trees: all three are done
   (§11.12; `ExternalCaptureExit`, `include/ExternalCaptureExit.h`; the trees' classes are gone).
2. Section 6, "Remaining steps": `ExitPairKeysInFile` no longer exists. `ParseExternalCapture` reads the exit keys from
   the format-3 specs (`ProjectExitPairKeys`, `src/ProjectSetup.cpp` 66–75), so a malformed segment line is the format-3
   reader's matter; nothing is left to decide. Step 6b is complete (steps 3 and 8 are done, §11.12).
3. §11.12, last paragraph: "`ParseExternalCapture` sets the level's capture pair twice ... not fixed" contradicts the
   bullet above it and the code: two photon pairs naming one final state are refused (`FinalStateOfAnotherPair`,
   `CNuc.cpp` 649–652, commit a629a391 "fix(cnuc): refuse two photon pairs that name the same final state";
   `tests/ec_final_state_ambiguous` (d), (e)).
4. `ROADMAP.md` (Phase 1, `Convention`) and `physics_system.md`: "the engine does not call it". The engine calls eleven
   functions of `level_parameters` (`amplitudes_from_observed`, `subtract_external_capture`, `formal_from_brune`,
   `brune_from_formal`, `park_from_brune`, `brune_from_park`, `boundary_values`, `observed_levels`,
   `internal_channels`, `external_reduced_width`, `capture_mask`, through `src/CNuc*.cpp` and `GammaDecayMap.cpp`).
   It does not call the convention part (`level_parameters()`, `convention_transform`, `project_channels`,
   `recoupling_coefficient`): no project states a convention, so there is nothing to convert. It rebuilds the
   `System` from `CNuc` at every use (twelve calls of `BuildSystemFromCNuc` counted here, nine in the engine at
   e9e37a18; one per `TransformIn`, `TransformToBrune`, `ObservedWidths`, `ConvertAmplitudeBasis`, set-up), since
   the structure it reads can still be written.
5. §2.1 lists `ScatteringLengthConstraints` and `HasScatteringLengthConstraints` on `CNuc`: gone. `CNuc.h` still calls
   `ReadAbundances` "TEMPORARY (step 3d)".

### 12.2 The end state

When `CNuc` is gone, the compound-nucleus side of the engine is six objects, one of which exists already:

- `azure2::physics::System` (exists) — one job: partitions, channels, J^π blocks of the R-matrix; owner: held by
  `CompoundSystem`, built once at set-up (two: R-matrix blocks, and all groups for the reports); immutable: yes; file:
  `physics/include/azure2/physics/System.h`.
- `CompoundSystem` — one job: the project's compound nuclei as set up: nuclei, pairs (`PPair`), J^π groups (`JGroup`)
  with channels (`AChannel`, B_c) and the levels' set-up data (`ALevel`: IDs, flags, external-capture final-state data),
  the pair-key map, the entrance pairs, the `System` with its index map, the external-capture set-up (formal-basis
  final-state blocks, capture-mask set-up), the model data, the `gamma_decay` set, the pathway table, the
  level-parameter layout; owner: `shared_ptr<const CompoundSystem>`, made by `BuildCompoundSystem`, held by the run or
  the session and shared by every evaluation, thread and walker; immutable: yes; file: `include/CompoundSystem.h`,
  `src/CompoundSystem.cpp` (queries), `src/CompoundSystemBuild.cpp` (from `Fill`, `AddHardSphereGroups`, `Initialize`),
  `src/ExternalCaptureFinalStates.cpp` (from `ParseExternalCapture`, `SetUpECFinalBlocks`, `SetUpFormalCaptureMasks`).
- `LevelInputs` — one job: the levels' starting values in the engine's basis (what `TransformIn` makes of `<levels>` or
  of an observed vector): E, γ, zero-input flags, the input external amplitudes and Park J, the `gamma_decay` numbers;
  owner: a value, held by the run or the session; the sampler's physical mode makes one per evaluation; immutable: yes
  (a value; a new physical vector makes a new one); file: `include/LevelInputs.h`, `src/LevelInputs.cpp`
  (`TransformIn`), `src/LevelInputMessages.cpp` (now `CNucTransformMessages.cpp`).
- `LevelState` — one job: the levels at one parameter vector and what follows from it: E_λ, γ_λc per group as flat
  arrays laid out like `JpiBlock` (λ n_c + k), the external-capture amplitude and N_f^½, S_c(E_λ) with the energy it is
  valid at, Park's J, dS/dE, d²S/dE², the formal capture-mask projectors, the `gamma_decay` fills, a version per level
  (the memo key `ChannelSpace` takes now from `ALevel::Stamp`); owner: the evaluation, like `Results` (one per pooled
  evaluation, thread or walker); immutable: no; file: `include/LevelState.h`, `src/LevelState.cpp`
  (`FillCompoundFromParams`, `CalcShiftFunctions`), `src/LevelStateExternalCapture.cpp` (`ECFinalStatesFromFormal`,
  `FormalCaptureMasksFromParams`).
- `ObservedLevelTable` — one job: the Brune (Park) and observed values of every level at a `LevelState`: E, γ, Γ, C, μ,
  Q, B, kinds, iterations, external amplitudes; the observed vector (`GetTransformParams`); owner: returned by
  `ObservedLevelsOf(system, inputs, state)`, held by its caller (reports, model data, the rate and sub-point grids, the
  sampler's output); immutable: a value; file: `include/ObservedLevelTable.h`, `src/ObservedLevelTable.cpp`
  (`TransformToBrune`, `ObservedWidths`).
- `LevelParameterLayout` — one job: where each level's energy and amplitudes sit in θ, with their records and kinds; the
  one walk of the level block; owner: member of `CompoundSystem`; immutable: yes; file:
  `include/LevelParameterLayout.h`, `src/LevelParameterVector.cpp` (`FillMnParams`, `ConvertAmplitudeBasis`, the level
  part of `BuildParamIndexMap`, `SetFitLevelParameter`).

The rest goes to objects that exist: the abundance samples (names, nominal values, priors, sum widths) to the
`DataSet` (they describe targets), their current values and the THM coherent values to the working copy (the data own
those parameters already: `WorkingCopy::FillThmCoherentFromParams`); the Park wall and the abundance priors and sums to
rows of the `Objective` (`src/Objective.cpp`), reading `LevelState` and the working copy; the reports
(`CNucReport.cpp`, `CNucChecks.cpp`) become `LevelReport.cpp` and `CompoundSystemChecks.cpp` reading
`CompoundSystem` and an `ObservedLevelTable`. An evaluation then reads `(const CompoundSystem &, const LevelState &)`
where it reads `CNuc *` today; set-up code reads `const CompoundSystem &`.

**Where each of today's members goes**: the "goes to" of each item of §12.1.1, and for `ALevel` (§12.1.2): (a) stays
on `ALevel` inside `CompoundSystem`; (a′) to `LevelInputs`; (b) to `LevelState`; (c) to `ObservedLevelTable`; (d)
deleted.

**The physics objects instead of second copies.**
- *Level-parameter transforms.* The engine already calls them (§12.1.5, item 4); what remains is glue that copies the
  fit set out of the levels into arrays for each call. `LevelState` holds the arrays in `JpiBlock`'s layout, so the
  calls take them directly. The convention part (`level_parameters()`) waits for a convention keyword in format 3
  (a new feature, not part of this work); after this work it has one caller to gain, `LevelInputs`' constructor.
- *`System`.* What stops the engine from holding one today: the structure it is built from can still be written
  after `Fill` (the entrance flags by the data set-up, B_c in `Initialize`), so each use rebuilds it. Step 1 freezes the
  structure and builds the `System` once; that closes the second copy in time. The engine's `PPair`/`JGroup`/
  `AChannel` remain as the set-up view the `System` is built from.
- *`JpiBlock`.* Its E and γ are values of a `System` (immutable); the per-vector values are `LevelState`'s arrays in
  the same layout, passed to the physics functions beside the `System`, as they already take them.

**Two designs for the structure; chosen: (1).**
1. *Keep `PPair`, `JGroup`, `AChannel`, `ALevel` (stripped to set-up data) as the engine's immutable view inside
   `CompoundSystem`, with the `System` built once from them.* The 907 + 539 reader calls change only from `CNuc *` to
   `const CompoundSystem &` and `const` pointers (mechanical; numbers read are the same objects). Cost: the structure
   is described twice (the view and the `System`), but the `System` is derived from the view once and both are
   immutable, so they cannot drift.
2. *The `System` is the only structure; the four classes are deleted*, the engine's extras (pair key, function
   settings, hybrid potential, Wigner limit, external-capture final-state data, level flags) become tables indexed by
   partition, block and level. Every reader changes from 1-based pointer accessors to 0-based indices: about 1,450
   calls in about 50 files, with an off-by-one risk at each, and the `System` gains fields `physics_system.md` §5
   keeps out on purpose (fit flags, the potential). Cost ≈ 3–4 agent-days and a conflict with every open branch that
   touches a matrix function.

(1) now; (2) when Phase 5 rewrites the readers for one collision-matrix solve per grid energy (§11.9), since those
readers change then anyway. Both give the state rule: an immutable shared structure and caller-owned state.

### 12.3 The steps

Each step is one branch an agent can do in a day or less. Steps 1–5 are bit-identical (every output file of the 183
projects × 3 bases against the reference tag; 1 against several threads); step 2b changes numbers on purpose. Each
bit-identical step ends with the pattern of §11.12 step 8: one commit in which what is moved away aborts if reached
(an accessor, a setter, a copy), the whole tree run identical, then the deletion in the next commit. Line counts are
estimates from §12.1.

**Step 1. The structure frozen after set-up; the `System` built once.** (≈ 300 lines changed, −80 net; 14 files:
`CNuc.h/.cpp`, `CNucSystem.*`, `CNucExternalCapture.cpp`, `CNucPathways.cpp`, `GammaDecayMap.cpp`, `WorkingCopy.cpp`,
`AZUREMain.cpp`, `AZUREAPI.cpp`, `PPair.*`, `ALevel.*`.)
- The entrance pairs become an input of the set-up, collected from the specs before `Initialize` (as the capture-A_y
  requests are, `RequestCaptureAyTerms`; `ProjectReactionPairs` and the rate's entrance), instead of six
  `SetEntrance` calls in the data set-up. First commit: both sets computed and compared, abort on a difference, the
  tree run; then the old calls go.
- `BuildSystemFromCNuc` runs once after `Fill` (the content `Structure` and `AllGroups`; neither holds B_c nor levels,
  `CNucSystem.h`), held as `shared_ptr<const CNucSystem>`; the twelve calls read it.
- Deleted: `PPair::SetBindingEnergy`, `ec_entrance_`, `ALevel::NumGammas`, `ALevel::SetShiftDerivatives`,
  `CNuc::GammaDecayInputs()`; `IsPair`, `AddPair`, `AddJGroup` private.
- The check: the non-`const` `GetPair`/`GetJGroup` abort after set-up in one commit (the tree run), then readers that
  only read take `const` pointers.
- Byte-identical because: the `System` is a function of `Fill`'s structure only (radii, hybrid potential, pairs,
  groups, channels, R-matrix flags), all fixed when `Fill` returns, so one build equals each of the twelve; the entrance
  set, if equal (the comparison commit shows it), gives the same pathway table and external-capture set-up.

*Step 1 as built* (branch `feature/upgrade-5-cnuc-frozen`). Measured first on the start commit (e9e37a18):

| what | design (6e8d000b) | code at the start | after |
|---|---|---|---|
| `SetEntrance` calls in the engine | 6 | 6 (and 2 in tests) | 0 |
| `BuildSystemFromCNuc` calls in the engine | 12 | 9 (and 5 in tests) | 2, in `CNuc::Fill` |
| dead members of §12.4 item 5 | 5 | all 5 present | deleted |
| `CNuc.cpp` / `CNuc.h` lines | 1,527 / 327 | 1,527 / 327 | 1,526 / 350 |

The six `SetEntrance` calls: `WorkingCopy.cpp` 4, `AZUREMain.cpp` 1, `AZUREAPI.cpp` 1; the nine builds: `CNuc.cpp` 5,
`CNucExternalCapture.cpp` 2, `CNucPathways.cpp` 1, `GammaDecayMap.cpp` 1.

- The entrances: `RequestDataEntrances` and `RequestExtrapolationEntrances` (`src/ProjectSetup.cpp`) request the
  entrance pair of every term (leading and components) of the active data sets or extrapolations with
  `CNuc::RequestEntrance`; a rate run requests its entrance pair. The comparison commit (37ac9885) kept the old flags
  and aborted on a difference at the pathway build and at the end of `WorkingCopy::Initialize`: all projects × 3 bases
  ran without an abort and identical. The program runs only menu 1 there: the extrapolation path was compared by
  the C++ tests that load projects without data (`upgrade_capture_ay_requests`, no abort); the rate path and the
  binding's sessions were not compared (the binding's tests of that run failed on an unrelated mismatch of the
  Python sources).
- The `System`: two are built once, at the end of `Fill` (`SystemStructure`, `SystemAllGroups`); a failure of the build
  is kept and rethrown by each reader, so a refused project stops where it did.
- The structure read-only: the design's check (the non-`const` `GetPair`/`GetJGroup` abort after set-up, the tree run)
  cannot be a run-time check: a reader on a non-`const` `CNuc *` calls the non-`const` overload, so it would abort in
  every reader. The compiler does it instead: the public accessors return `const` pointers only, the set-up writes
  through the private `Pair`/`Group`, and the level values that move after set-up (fit values, shift functions,
  report slots; the `gamma_decay` input widths at set-up) through `CNuc::MutableLevel`, whose callers outside the class
  are `TargetSetup.cpp`, `ReactionRate.cpp` (the sub-point-grid and rate copies), `GammaDecayMap.cpp` and two tests.
  The design had missed that the levels are reached through `GetJGroup` and written after set-up.
- Found on the way: `WorkingCopy::InitializeLeaves` flagged entrance pairs during `WorkingCopy::Initialize`, after
  the pathway table was built (the comparison showed it never added one); `tests/upgrade/channel_table_test.cpp`
  changed a channel radius after set-up (it fills a second compound with the radius override now).

**Step 2. Abundances, THM coherent values and the two penalties out of `CNuc`.** (≈ 350 lines moved; `CNuc.h/.cpp`,
`ProjectSetup.cpp`, `Objective.cpp/.h`, `AZUREGrad.cpp`, `ParameterLabel.cpp`, `ModelPass.cpp`, `WorkingCopy.*`,
`DataSet.*`, `ThmNorm.cpp`, `THMMatrixFunc.cpp`, `AZUREAPI.cpp`, `AZURECalcMCMC.cpp`.)
- The samples (`AbundanceSample`, `SampleAbundance`, `SetAbundances`) into the `DataSet`; the current abundances and
  the THM coherent values into the working copy; `WorkingCopy::FillMnParams` creates the abundance parameters first,
  at the position `CNuc::FillMnParams` gives them now (right after the level block), so θ's layout is unchanged.
- `ParkNormPenalty`, `AbundancePenalty`, `AbundanceSumPenalty` become the bodies of their `Objective` terms (the terms
  already exist and call these members).
- Byte-identical because the same numbers are summed in the same order by the same expressions; the parameter names
  and records are made by the same functions.
- Can run side by side with step 1 (different functions of `CNuc.cpp` and different lines of `AZUREAPI.cpp` and
  `WorkingCopy.cpp`; merge conflicts at most in includes).

*Step 2 as built* (branch `feature/upgrade-5-cnuc-passengers`, start adfb0118). Measured first:

- Calls outside `CNuc` that put the abundance block into a vector or read it back (the design did not count them):
  17 at the start (`GetTransformParams` / `FillCompoundFromParamsPhysical`: `AZUREAPI.cpp` 7, `AZURECalcMCMC.cpp` 4,
  `ParameterLimitsManager.cpp` 2, `SamplerSupport.cpp`, `AZURE2.cpp`, `TargetSetup.cpp`, `ReactionRate.cpp`); after:
  one place, `AbundanceSamples` (`AddParameters`, `AppendTo`, `ReadFrom`).
- THM coherent values: held twice at the start (`WorkingCopy`, `CNuc`; `THMMatrixFunc.cpp` read the `CNuc` copy at
  two places); after: once, in the working copy.
- `CNuc.cpp` / `CNuc.h` lines: 1,527 / 327 in the design, 1,526 / 350 at the start, 1,355 / 272 after.

- The design missed that the abundances travel through both vectors (the fit's and the observed one). They have one
  place now: `AbundanceSamples` (`include/AbundanceSamples.h`) holds the samples (names, the file's values, priors,
  sum widths) and writes the block into the fit's vector (`AddParameters`, called by `WorkingCopy::FillMnParams`
  first, at the position `CNuc::FillMnParams` gave it), into an observed vector (`AppendTo`: the block of the fit's
  vector it was made from, called by the seven places that hand out an observed vector) and reads it from either
  (`ReadFrom`, through `WorkingCopy::FillAbundancesFromVector`, called by `FillEnergyShiftsFromParams` and after a
  compound is filled from an observed vector).  The other five callers of `GetTransformParams` round-trip the level
  block only and read no abundance.
- The samples are the `DataSet`'s (`DataSet::Abundances`, set by `WorkingCopy::SetAbundanceSamples` at the start of
  the data set-up, `src/ProjectSetup.cpp`); the present values and the THM coherent values are the working copy's
  (`AbundanceValues`, `ThmCoherentValues`).  The HOES model gets the coherent values as an argument
  (`SolveAtOneEnergy`, `HoesCrossSection`, `ThmAdjointContext::coherentValues`).  `GradAccum`'s abundance array
  grows with the ordinals added, as its screening array does (the compound no longer knows the count).
- The comparison commit (2e17e299) kept the compound's copies and aborted wherever the two differed.  It found the
  places where only the compound had been filled from a vector and something read the value later: the session's
  `SetParamsAllRWA` and its three `Transform...` functions (read by the pair record's abundance), and the sampler's
  reference test (`tests/reference/sampler_objective_test.cpp`).  These read the block into the working copy too.
- The penalties are the bodies of their terms: `ParkTerm` (T6) computes the wall, `AbundanceTerm` (T7) the priors
  and sums; `kParkNormScale` is a constant of `Objective.h`.
- One visible change: the reader's messages of the `<abundances>` block (a row naming a nucleus no pair forms) print
  at the data set-up, after "Filling Data Structures...", instead of in `CNuc::Fill`; a rate run (no data) no longer
  reads the block.  No project of the tree prints one.

**Step 2b. Defect 1 fixed (changes numbers; its own step).** (≈ 30 lines in `TargetSetup.cpp` and a test project.)
The sub-point-grid copy takes the evaluation's N_f^½ and external-capture amplitudes with E and γ (§12.4, defect 1).
The acceptance rule for a step that changes a number: per output file, the largest distance to the refined reference
at most 1.2 times the old engine's and the rms at most 1.03 times, or at most 10⁻⁷; every file of the tree that moves
named with its size and cause. Side by side with steps 1 and 2 (one file of its own); before step 3, so that steps 3
and 4 need no exception to bit-identity.

**Step 3. `LevelState` and `ObservedLevelTable`: the fit and report sets out of `ALevel`.** (≈ 800 lines changed:
about 200 reader calls; `ALevel` −250 lines; new `LevelState` ≈ 350 and `ObservedLevelTable` ≈ 200 lines; ≈ 25 files:
`ALevel.*`, `CNuc.*`, `CNucExternalCapture.*`, `CNucObserved.cpp` (deleted), `CNucReport.cpp`, `CNucChecks.cpp`,
`CNucTransformMessages.cpp`, `AMatrixFunc.cpp`, `RMatrixFunc.cpp`, `THMMatrixFunc.cpp`, `ChannelSpace.cpp`,
`AZUREGrad.cpp`, `Objective.cpp`, `ModelData.cpp`, `GammaDecayMap.cpp`, `TargetSetup.cpp`, `ReactionRate.cpp`,
`ThmLineshape.cpp`, `ThmVertexBoundary.cpp`, `ThmReports.cpp`, `AdaptiveExtrapolation.cpp`,
`AdaptiveIntegrationGrid.cpp`, `AZUREAPI.cpp`, `AZURECalcMCMC.cpp`.) After step 1.
- `CNuc` holds one `LevelState` while the owners do not change (its copy still copies it); `FillCompoundFromParams`,
  `CalcShiftFunctions`, `ECFinalStatesFromFormal`, `FormalCaptureMasksFromParams`, `SetFitLevelParameter` write it;
  every reader of the fit set reads it (`state.E(j, λ)`, `state.Gamma(j, λ, c)`, ...). `ChannelSpace` keys its memo on
  (state, level, version) instead of (level address, stamp).
- `TransformOut` returns an `ObservedLevelTable`; `GetTransformParams`, the reports, `ModelData`, `ObservedWidthsAt`
  read it. `ObservedOfFitValues` no longer saves and restores N_f^½ (it computes a table and changes nothing).
  `THMMatrixFunc`'s numerical width derivative perturbs a copy of the level's values, not of the `ALevel`.
- Byte-identical because every value is computed by the same function from the same inputs in the same order; the
  memos (`ChannelSpace`, the levels' shift-function energies, `RunScratch::shiftDerivatives`) are keyed exactly, so a
  hit returns what the computation gives and a different key changes only whether it is recomputed (with step 2b
  merged, no copy reads a stale value).
- The check: the moved `ALevel` accessors abort (one commit, the tree run), then deleted; `ALevel` keeps (a) only.

**Step 4. The parameter vector; the copies replaced.** (≈ 600 lines changed; `CNuc.cpp` −400 into
`LevelParameterVector.cpp` and `LevelInputs.cpp`; files: `CNuc.*`, `AZUREGrad.cpp`, `AZUREParams.cpp`,
`ParameterLabel.cpp`, `AZURECalc.*`, `AZURECalcMCMC.*`, `AZUREAPI.*`, `AZUREAPIModelData.cpp`, `AZURE2.cpp`,
`AZUREMain.cpp`, `ReactionRate.cpp`, `TargetSetup.cpp`, `ParameterLimitsManager.cpp`, `CovarianceBand.cpp`,
`AdaptiveExtrapolation.cpp`, `ModelData.cpp`, `tests/consistency/allowed_state.tsv`.) After steps 2 and 3.
- `LevelParameterLayout` built at set-up; `FillMnParams`, `FillCompoundFromParams` (now `FillLevelState`),
  `FillCompoundFromParamsPhysical` (now a `LevelInputs` from an observed vector), `GetTransformParams`,
  `ConvertAmplitudeBasis`, `SetFitLevelParameter` and the level part of `BuildParamIndexMap` read it (one walk instead
  of eight).
- `TransformIn` returns a `LevelInputs`; the sampler's physical mode makes one per evaluation.
- The owners hold the `CNuc` as `shared_ptr<const CNuc>` plus their `LevelInputs`; `cnuc_pool_` (fitter and sampler)
  becomes a pool of `LevelState`; `Clone` is deleted; `ReactionRate`'s and `TargetSetup`'s copies become
  `ObservedLevelsOf` on the evaluation's state (the same values since step 2b). The two `cnuc_pool_` rows
  and `pool_mutex_` of `allowed_state.tsv` (milestone 4) go.
- Byte-identical because θ's order and records are unchanged (the layout lists the slots in the order of today's
  loops, `BuildParamIndexMap`'s comment says it mirrors them), and each evaluation reads the same values from a
  `LevelState` as from its copy; the check: `Clone` aborts in one commit (the tree, `ctest -LE slow`, the sampler's
  chains of the check scripts identical), then deleted.

**Step 5. The set-up out and the name gone.** (≈ 900 lines touched, mostly signatures in ≈ 110 files; `CNuc.cpp`
→ `CompoundSystemBuild.cpp` and `ExternalCaptureFinalStates.cpp`; the `CNuc*` files renamed: `CNucReport` →
`LevelReport`, `CNucChecks` → `CompoundSystemChecks`, `CNucSystem` → `SystemOfCompound`, `CNucLevelParameters` →
`EngineChannelFunctions`, `CNucPathways` into the build, `CNucTransformMessages` → `LevelInputMessages`.) Last.
- `Fill`, `AddHardSphereGroups`, `ParseExternalCapture`, `ReadAbundances`/`SetAbundances` (already moved in step 2),
  `Initialize` become `BuildCompoundSystem(const Config &, const SetupRequests &)` returning
  `shared_ptr<const CompoundSystem>`; `SetupRequests` carries what the data set-up collects first (entrance pairs,
  capture-A_y terms). Every `CNuc *` becomes `const CompoundSystem &` (set-up and readers) and, in evaluation
  signatures, `const LevelState &` beside it.
- A consistency check (the next free number after the highest in `repo_consistency_test.py`) refuses the names `CNuc`,
  `Clone()` of the compound and `cnuc_pool_` in `src/`, `include/`, `api/`, as check 31 does for the data classes; the
  `file_sizes.tsv` row of `CNuc.cpp` is removed by the director (or lowered with `--tighten-sizes` on the way).
- Byte-identical because it renames and moves code that steps 1–4 left without state: no function body changes but for
  the receiver of its accessors. The check is the build itself plus the tree.

**Side by side.** Steps 1, 2 and 2b together; then 3; then 4; then 5. Step 3 cannot run beside 1 (both edit every
external-capture function and `CNuc.cpp`'s fills), nor 4 beside 3 (the same owners' lines). Step 5 is the conflict
magnet (≈ 110 files): it should start when no other branch edits the engine's signatures (the sampler, the records,
the session options and the physics table, §12.6, merged first) and be merged within the day.

**Total.** ≈ 2,980 lines changed over six branches; `CNuc.cpp` (1,527) and `CNuc.h` (327) deleted, ≈ 1,700 lines of
new files (each under 1,000), net ≈ −500 lines; `ALevel` from 721 to ≈ 450 lines.

**What would change a number.** Nothing in steps 1–5 by design; step 2b by design. Named risks, each with the check
above: the entrance set from the specs differing from the data set-up's (step 1: then the pathway table differs; the
comparison commit finds it before anything moves); a reader that read a fit value after a report had overwritten a slot
it shares with the calculation (step 3: only N_f^½ is shared, by `ObservedOfFitValues` and `TransformToBrune` in the
formal basis, and both restore or recompute it; a difference would show in the formal-basis external-capture projects);
the order of the abundance parameters (step 2: the records' names and `param.par` would differ). The fix of defect 1
changes numbers and is its own step, 2b.

### 12.4 Defects in the way

1. **The sub-point-grid copy reads a stale N_f^½** (found while reading; not run). `WorkingCopy::ReplaceSubPointGrids`
   (`TargetSetup.cpp` 190–202) keeps one clone of the compound for the whole run (`gridNuc_`), copies only E and γ of
   the fit into it before each `TransformOut`, and the observed photon widths with external capture are formed with the
   final state's N_f^½ the clone had when it was made (`ExternalReducedWidth`, `initial = false`, reads
   `GetSqrtNFFactor()`; `TransformToBrune` refreshes it only in the formal and Park bases). In the Brune basis, with
   external capture and a fitted final-state amplitude, the grids are laid out from widths at the start values of N_f.
   `ReactionRate`'s copy is made fresh at each call, so it does not have the defect. Fixed before step 3, in its own
   step (2b) with the acceptance rule; a test project with a fitted final state and a target effect, and the run of
   the tree naming every file that moves (in a menu-1 calculation the copy is made at the vector it is used at, so
   only fits and the API's repeated evaluations can move).
2. **Decay position of the external-capture amplitudes**: fixed (§11.12, decision 3 of §11.11). Nothing to do.
3. **Two photon pairs with one final state** (`ParseExternalCapture`): refused since a629a391. Nothing to do; the
   contradicting paragraph of §11.12 is stale (§12.1.5).
4. **`ExitPairKeysInFile` and a malformed segment line**: the function is gone; the exit keys come from the format-3
   specs (§12.1.5). Nothing to do.
5. **Write-only and dead members** (§12.1.2): `PPair::ec_entrance_`, `PPair::SetBindingEnergy`, `ALevel::NumGammas`,
   `ALevel::SetShiftDerivatives`, `CNuc::GammaDecayInputs()`. Deleted in step 1.
6. **A report that writes the calculation's state.** `ObservedOfFitValues` saves and restores N_f^½ around
   `TransformToBrune`, and `TransformToBrune` in Park's basis calls `CalcShiftFunctions` (it fills the second
   derivatives into the evaluation's copy as a side effect). Not wrong today; removed by step 3.
7. **A positional boolean** `ExternalReducedWidth(..., bool initial, ...)` (`CNucExternalCapture.cpp` 142) chooses
   between the input set and the report set; step 3 replaces it by the two sources it reads from.
8. **The parameter layout written eight times** (`FillMnParams`, `FillCompoundFromParams`,
   `FillCompoundFromParamsPhysical`, `GetTransformParams`, `ConvertAmplitudeBasis`, `BuildParamIndexMap`,
   `ParameterLabel.cpp`'s two walks); kept in lock-step by comments. Step 4.

Radiative widths stay out of a level's total width (the PI's decision); nothing here touches it.

### 12.5 The PI's measure, before and after

- **(i) a new kind of level parameter (e.g. a reduced-width amplitude in another convention, or a new map like
  `gamma_decay`)** — places today: the record (`Schema.cpp`, its build in `experiment/`), the map file, a member of
  `CNuc` beside `gammaDecays_` and its lines in `Fill`, `FillMnParams`, `FillCompoundFromParams`,
  `FillCompoundFromParamsPhysical`, `GetTransformParams`; `ParameterKind.h`; `BuildParamIndexMap`; `AZUREGrad`'s
  pull-back; the objective's prior rows; `CNucReport.cpp`; the binding: ≈ 12 files; after step 5: the record and its
  build, the map file, one entry in the list of maps `CompoundSystem` holds (which `FillLevelState`, `LevelInputs` and
  the layout walk), `ParameterKind.h`, the pull-back, the report, the binding: ≈ 8 files. A convention keyword
  specifically: the record, the converter, `LevelInputs`' constructor calling `level_parameters()`.
- **(ii) a new channel type (the spin–orbit (α, l, j) channel of `spin_orbit_potential.md` option (b))** — places today:
  the `System` (`ChannelKind`), `AChannel`, `Fill`, the adapter `CNucSystem.cpp`, the channel table key (done for 2J),
  and 155 inline tests of the pair or radiation type in 37 files, of which the matrix functions hold 39; after step 5:
  the same, but for `Fill` → `CompoundSystemBuild.cpp`: the end of `CNuc` does not touch the inline tests. Removing them
  needs a channel-kind table (P, S, φ, the observed-width factor, predicates) read by the matrix functions: a Phase 5
  item with the readers' rewrite.
- **(iii) a new particle-pair property read from the project** — places today: as the light-particle excitation
  (`c49acd17`): the `<levels>` column (`io/LevelRows`), `LevelTables`/`LevelEdits` in `experiment/`, `PPair.h/.cpp`,
  `CNuc::Fill` (its consistency check), `CNucSystem.cpp`, `CNucReport.cpp`, the binding's pair record and pyazr
  (`nuclear_part.py`, `scheme.py`, `parameters.py`, `azure2.py`), its readers: ≈ 20 files; after step 5: the same less
  one (`CNuc.cpp` → the build file); the engine side stays: `PPair`, the build, `SystemOfCompound.cpp` if physics reads
  it, the report. No gain worth claiming.
- **an AI agent changing the level transforms safely** — places today: physics `LevelParameters`, `ExternalCapture`
  (1,353 lines), the adapters `CNucLevelParameters`, `CNucSystem` (381), the transform members of `CNuc.cpp` (532),
  `CNucExternalCapture` (463), `CNucObserved` (49), `CNucTransformMessages` (222), `ALevel` (721, to know which of its
  three sets a function reads), `GammaDecayMap` (930): ≈ 4,650 lines in 18 files, and the rule which object a function
  may write is nowhere written down; after step 5: physics (1,353, unchanged), `LevelInputs` with its messages (≈ 450),
  `LevelState` with its external capture (≈ 550), `ObservedLevelTable` (≈ 200), the channel-functions adapter (≈ 160),
  `GammaDecayMap` (930): ≈ 3,650 lines in 12 files, and the rule in the types: inputs and the structure are `const`,
  only `FillLevelState` writes a `LevelState`, a report returns a table.

Honestly: the end of `CNuc` improves (i) and the level transforms, and the state rule everywhere; it does not improve
(ii) or (iii) much, because their cost is in the inline channel-type tests and in the number of front ends, not in
`CNuc`.

### 12.6 Costs and risks

**Lines.** Step 1: ≈ 300 changed, −80; step 2: ≈ 350 moved, ±0; step 2b: ≈ 30 and a test; step 3: ≈ 800 changed,
−150 (`ALevel` −250, new files +550, readers ±0); step 4: ≈ 600 changed, −250 (eight walks become one); step 5:
≈ 900 touched, ±0 (renames).

**The file-size ratchet** (`tests/consistency/file_sizes.tsv`): `CNuc.cpp` (1,548 recorded, 1,527 now) only shrinks
and is deleted; `AMatrixFunc.cpp` (2,766), `AZUREAPI.cpp` (1,779), `AZURECalc.cpp` (2,018), `AZURECalcMCMC.cpp`
(1,433), `AZUREGrad.cpp` (1,273), `WorkingCopy.cpp` (1,372) may not grow: steps 2–5 replace calls one for one or
remove walks, so none grows if the agent does not add helpers there; a new file stays under 1,000 lines.

**Copy size and time.** A copy of today's `CNuc`: 7.8 / 12.9 / 131 kB and 5 / 7.5 / 56 µs (T2 / 13N / `natural_neon`);
after: a `LevelState` of 0.8 / 1.4 / 22 kB, made by a refill (0.9 / 1.1 / 9.9 µs) instead of a copy. **The sampler**
gains a walker's copy that holds only its state and cannot write the shared model, and loses the snapshot's clone
every 10 steps; no measurable time. **The web UI's sessions** may share one `CompoundSystem` when they open the same
project and hold a `LevelState` each; the gain is memory (≈ 110 kB per session for `natural_neon`) and the guarantee
that two sessions evaluating at once do not see each other (today they do not share a `CNuc`, so they do not either).

**Other open work on the same files** (read in their worktrees; not merged at 6e8d000b; what is assumed):
- *The sampler* (`feature/upgrade-6-sampler-objective`, HANDOFF entry on its branch): rewrites `AZURECalcMCMC.cpp`
  (the posterior from the objective's terms, the 17 globals into a run record) and `Objective.cpp/.h`,
  `ParameterKind.h`, `AZURE2.cpp`. Assumed merged before steps 2 and 4: step 2 then edits only `Objective.cpp`'s terms
  for the penalties (the sampler's own calls of `ParkNormPenalty`/`AbundanceSumPenalty` at lines 72–175 go with that
  branch), step 4 replaces its `cnuc_pool_`.
- *The parameter records* (`feature/upgrade-8-pyazr3-records`): deletes `AZUREAPI::GetParameterInfo`/`GetPairsInfo`
  and adds `pybind_records.cpp` (a table of the fields of a parameter record and a pair record, with "the reader of the
  engine object"). Assumed merged before step 2: the abundance and pair readers of that table then follow the moves
  (its pair rows read `PPair`, unchanged by these steps; its abundance rows read `CNuc::GetAbundance` today).
- *The session options* (`feature/upgrade-4-session-options`): touches `AZURE2.cpp` and `ProjectOptions` only; no
  conflict beyond step 5's signatures.
- *The physics table* (`feature/upgrade-8-pyazr3-physics`): its two engine calls (`AChannel::SetWignerLimit` through
  `physics/Widths.h`, `ChannelFunctionsFromCNuc`'s photon penetrability in `src/CNucLevelParameters.cpp`) stay where
  they are through steps 1–4; step 5 renames the second file. Assumed merged before step 5.
- *Phase 6, milestone 4* (`evaluation.md` §5): "`CNuc` split, the parameter vector out" is step 4 here; its gate (memory
  of a fit at 1, 3 and 12 threads equal within the `Results` size per thread) is step 4's speed and memory measurement.

### 12.7 What this design rests on, and decisions

**Taken from notes or inferred, not checked by a run:**
- that the six `SetEntrance` calls mark exactly the entrances of the active data sets, leaves, extrapolations and the
  rate (step 1's comparison commit decides it);
- that `boundary_values` and the other users of the `Structure` content read no B_c from the `System` (the adapter's
  comment says the content holds B_c = 0);
- defect 1 (read in `TargetSetup.cpp`, `CNucExternalCapture.cpp`, `CNuc.cpp`; no project was run to show it);
- that `ChannelSpace`'s memo returns only exact values (its key holds every input, by reading of `ChannelSpace.cpp`);
- the line estimates of §12.3 and §12.6 (from the counts of §12.1, not from a trial cut);
- the open branches' content (from their diffs and HANDOFF entries, not built).
**Measured**: the counts of §12.1 (scripts named there); the bytes and times of §12.1.4 (`clone_probe`, one thread);
the clone and refill counts of one calculation of T2, 13N and `natural_neon` (`gdb` breakpoints on `b-trial50`).

**Decisions for the director or the PI** (each with the choice taken here):
1. *Structure: keep the engine's view classes inside `CompoundSystem` (design 1) or replace them by the `System` now
   (design 2)?* Chosen: design 1; the `System` replaces them in Phase 5 with the readers' rewrite.
2. *Names.* `CompoundSystem`, `LevelInputs`, `LevelState`, `ObservedLevelTable`, `LevelParameterLayout`; `PPair`,
   `JGroup`, `AChannel`, `ALevel` keep their names (renaming them touches ≈ 1,450 lines of every open branch for no
   change of meaning). Alternative: rename them in step 5.
3. *Abundances: owned by the data (samples in the `DataSet`, values in the working copy) or by `CompoundSystem` and
   `LevelState`?* Chosen: the data, since a sample describes a target and a segment chooses it; θ's order is kept.
4. *A `LevelState` per evaluation or a pool?* Chosen: a pool of `LevelState` in step 4 (the shift-function memo
   keeps its hits, the order of evaluations as today); a fresh state per evaluation (`evaluation.md`'s `Evaluation`)
   costs one refill with all shift functions (16–81 µs here) and is Phase 6's choice.
5. *Defect 1*: fix it in its own step before step 3 with the acceptance rule (chosen), or keep it and name it (then
   steps 3 and 4 must copy it on purpose to stay bit-identical).
6. *Order against other work*: the sampler and the records branches merged before step 2, the physics table before
   step 5 (chosen), so that no step of this list rebases over them.
