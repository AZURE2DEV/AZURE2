# Splitting `EData`, `EPoint` and `ESegment`: plan

`EData` (`include/EData.h`, `src/EData.cpp`, `src/EDataThm.cpp`, `src/OldLineSetup.cpp`) is the data set of a
project: its segments, their component segments, the target effects, the parameter layout of the data-side
nuisances (norms, energy shifts, densities, Doppler temperatures, screening potentials, sample thicknesses), the THM
experiments and the data covariance; it also lays the sub-point grids out, runs the per-stage tables, hands a parameter
vector over to the points, and writes the output files. `ESegment` (`include/ESegment.h`, `src/ESegment.cpp`) is one
data or extrapolation set: its settings, its points by value, its normalization and energy-shift state, the channel
tables of its last shift states, and its composition (SUM, RATIO, total capture, sample, transform). `EPoint`
(`include/EPoint.h`, `src/EPoint.cpp`) is one measured point and one calculation unit: the datum, its c.m. frame, the
energy-dependent tables (L⁰, √P, Coulomb and hard-sphere phases, Legendre polynomials, Coulomb amplitude), the
external-capture amplitudes, its sub-points and the apparatus folds over them, the mapping onto other points and the
offsets of its result record. At 070756a1 the files have 5,571 (`EData.cpp`), 1,233 (`EDataThm.cpp`), 3,006
(`EPoint.cpp`), 1,446 (`ESegment.cpp`), 546, 905 and 602 lines (the three headers). The aim of the upgrade is files
with one job, small enough that one file is read to change one thing (no file above about 800 lines). This page lists
every member of the three classes with its responsibility, callers and destination, fixes which parts move now, which
wait for the step that deletes or rewrites them, and orders the moves so that each is small, bit-identical and
mergeable alone. It follows the form of `cnuc_split.md`; it does not design the calculation plan (S5) or
`Measurement` (S7), which are in `data_model.md`.

## 1. Rule for each move

A **move** cuts member definitions (or a file-local helper together with every user of it) from one `.cpp` file and
pastes them, unchanged, into another `.cpp` file of the same class; the old text is deleted in the same commit, the new
file is added to `AZURE2_CORE_SOURCES` (`src/CMakeLists.txt`, on a line of its own: the running branches edit line 10),
and nothing else changes: no signature, no evaluation order, no inline/out-of-line status, no `static` or
`thread_local` object split from its users. A **split in place** turns one function into several in the same file with
the statements in the same order. A **deletion** removes a member that has no caller or no reader. Each is accepted
when (a) the build has 0 warnings; (b) the 179 projects × 3 bases are identical to `tab-p3f`
(`tools/runall.sh`, `tools/cls.sh`: 689 / 662 / 689 files); (c) `ctest -LE slow` and `tests/run_tests.sh` pass;
(d) calculations at 1 and 8 threads are identical (the director's run); (e)
`python3 tests/consistency/repo_consistency_test.py` passes and `--tighten-sizes` lowers the entries of the files that
shrank; (f) the step's own check of §5 passes — every step names the outputs that (b) does not reach.

## 2. Inventory at 070756a1

### 2.1 How the numbers were obtained

All numbers of this page were counted at commit 070756a1 by scratch scripts, not committed
(`defs.py`, `inventory.py`, `spans.py`, `tables.py`, `members.py`):
top-level definitions are found by brace matching outside comments, strings and character literals; a header's members
by the same matching inside the class body. "Lines" of a function is its definition from signature to closing brace;
the **span** of a move (§2.3, §5) adds the comment and blank lines between it and the definition before it, which move
with it. "Callers int/data/ext" counts call sites `Name(` (and `Class::Name(`) in `src/`, `include/`, `api/src/`,
`api/include/` and `tests/` (`*.cpp`, `*.h`), outside the definition and its declaration: *int* in the class's own
files, *data* in the files of the other two classes, *ext* everywhere else, with the four files of most ext calls in
brackets. A name that another class also defines (`Fill`, `GetJ`, `Initialize`, …) is over-counted and marked ≈. Calls
written `Name<V>(` are not found: the three templated folds `IntegrateTargetEffectImpl`, `IntegrateKernelsImpl`,
`TargetFallbackImpl` show 0 and are called at `EPoint.cpp` 2456, 2470, 2561, 2606, 2700, 2843, 2848, 2854.
"State named" is the number of the class's own data members named in the body, and how many of them are written (an
assignment, `++`/`--` or a mutating container call is seen; reads and writes through accessors of other objects are
not counted).

Totals: `EData` has 139 member definitions (overloads counted apart; 105 out of line with 5,053 lines, 34 in the header
with 50) and 43 data members; `EPoint` 226 (143 out of line, 2,038 lines; 81 in the header, 110 lines; 2 declared and
never defined) and 93 data members; `ESegment` 187 (132 out of line, 906 lines; 55 in the header, 78 lines) and 98 data
members. None of the three declares a friend.

### 2.2 Responsibility classes

| code | responsibility |
|---|---|
| S | reading and set-up: building segments and points from specs, load-time refusals and checks |
| K | the point's kinematics and frames: lab → c.m. energy, angle, Jacobian, the kept load-time frame |
| T | energy-dependent tables of a point and their orchestration per stage: channel rows, Legendre polynomials, Coulomb amplitude, the per-stage channel table |
| X | external-capture amplitudes and the `intEC` file (signature, read/write, progress) |
| F | apparatus folds: target effects, sub-point grids (placement and re-placement), the kernels that fold sub-point values |
| P | the per-point physics pass: `EPoint::Calculate` (dispatch to the R/A/T-matrix functions in `GenMatrixFunc`, `AMatrixFunc`, `RMatrixFunc`, `THMMatrixFunc`), probes, the parallel pass |
| H | THM experiments: set-up of their physics objects, the profiled norm and background, reports and tables |
| B | β-delayed data: no member of its own; the decay-energy conversion (`EPoint::ConvertDecayEnergy`), one refusal (`BetaDelayedAllowed`) and `GetPType() == 20` branches in S, F and X code |
| O | outputs and printing: `AZUREOut_*`, `.band`, `chiSquared.out`, `normalizations.out`, `shifts.out`, check files |
| V | fit plumbing: parameter vector layout and names, nuisance blocks, norm/shift/thickness groups, priors and penalties, the hand-over of a parameter vector to segments and points |
| E | an energy shift applied to the points (re-energising, refill tables) |
| C | composition: component segments, SUM, RATIO, total capture, sample weights, transmission and yield |
| M | point mapping and cloning |
| R | result records: the layout of `Results` and a point's offsets and accessors into it |
| D | data covariance groups and seeds |
| I | containers, iteration and run flags |
| A | plain accessors of the datum and of the segment settings, including the kind flags |

### 2.3 Partition of the source files

Each `.cpp` file divides into these groups without remainder (spans of §2.1; the first group of a file carries its
`#include` block).

| file (lines) | group: lines |
|---|---|
| `src/EData.cpp` (5,571) | S refusals with the includes (1–230, 266–325): 290 · S `Fill`, `ApplyPriorCentres`, `MakePoints`: 352 · I containers, flags, target-effect list: 227 · F target windows, grids, `SetUpTargetEffects`, `ReplaceSubPointGrids`: 938 · S/T `Initialize`: 137 · T stage orchestration (`CalcLegendreP`, `BuildChannelTable`, `CalcEDependentValues`, `CalcCoulombAmplitude`): 101 · O writer (`WriteNaturalDecomposition`, `WriteBandLine`, `WriteOutputFiles`): 490 · O check files (`Print*`): 240 · X external capture: 501 · C components: 457 · M mapping: 267 · V parameters and groups: 660 · V/E hand-over: 407 · D covariance: 259 · R result layout: 52 · M/P `Clone`, `CalculatePointsInParallel`: 193 |
| `src/EDataThm.cpp` (1,233) | H set-up (helpers, `BuildThmGroups`, `FillThmCoherentFromParams`): 651 · H profile (`ThmGroupOf` … `ThmBackgroundAt`) and `ThmExperimentReports`: 121 · H/O `WriteThmExperiments` and the three tables: 461 |
| `src/EPoint.cpp` (3,006) | constructors with the includes: 201 · A/S accessors (207–602): 331 · T table getters: 86 · X external capture: 248 · `Initialize`: 23 · K kinematics: 209 · setters (934–1036): 108 · T Legendre and Q coefficients: 133 · T channel tables and Coulomb amplitude: 571 · P calculation (`ScopedPointMappingOff`, `MarkFailed`, `Calculate`, thin target, observable integration): 394 · M mapping: 55 · H/F THM setters, `SetTargetEffectNum`, `AddSubPoint`: 56 · F folds: 463 · R and sub-point accessors (2863–3006): 128 |
| `src/ESegment.cpp` (1,446) | S construction (`SetReaction`, two constructors): 182 · A settings accessors: 296 · H `ProfileNormChiSquared`: 33 · S/K `Fill` and conversion: 136 · V shift parameters and penalties: 233 · E shift application and refill tables: 270 · I setters and containers: 118 · C composition: 148 · H THM attachments: 30 |

`src/EDataIterator.cpp` (113 lines) and `include/EDataIterator.h` (38) are iteration only (I) and already apart;
`src/OldLineSetup.cpp` (53) is the old-format entry (`EData::Fill`/`MakePoints` overloads), deleted by step 3d.

### 2.4 Member functions

Rows of one function of 15 lines or more; smaller members of one class and destination are grouped (lines summed,
callers summed). Destinations name the step of §5 (N = now, W = after the current wave, S5/S7/S12/3d = after that step,
D1 = the moves that wait for the director's decision D1, Phase 5/6).

#### `EData`

| member(s) | lines | where | class | callers int/data/ext (main ext. files) | state named | destination |
|---|---|---|---|---|---|---|
| `NumScreeningParams`, `GetScreeningEffect`, `GetScreeningTargetEffect`, `GetScreeningTargetEffect`, `NumDensityParams`, `GetDensityEffect`, `GetDensityTargetEffect`, `GetDensityTargetEffect`, `NumDopplerParams`, `GetDopplerEffect`, `GetDopplerTargetEffect`, `GetDopplerTargetEffect`, `NumThicknessParams`, `GetThicknessSegment`, `SetPriorCovariance`, `GetPriorCovariance`, `GetNormParamOffset`, `GetEnergyShiftParamOffset`, `GetEnergyShiftSqrtParamOffset`, `GetEnergyShiftFnParamOffset`, `SetNormParamOffset`, `SetEnergyShiftParamOffset`, `SetEnergyShiftSqrtParamOffset`, `SetEnergyShiftFnParamOffset`, `SyncShiftGroupParams`, `NormOwner`, `SegmentParameterName`, `ThicknessParamOfSegment`, `DopplerParamOfEffect` | 85 in 29 | EData.cpp, EData.h | V | 81/0/81 (sum) | — | after S7: `src/EDataParameters.cpp` |
| `GetScreeningParamOffset`, `GetDensityParamOffset`, `DensityParamOfEffect`, `GetDopplerParamOffset`, `GetThicknessParamOffset` | 7 in 5 | EData.h | V | 0/0/0 (sum) | — | delete (no caller), W1 |
| `NormUsed`, `SizeResults`, `ForEachResultUnit`, `BeginEvaluation` | 35 in 4 | EData.cpp, EData.h | R | 10/0/4 (sum) | — | S5 (plan) |
| `SetCovarianceError`, `HasCovariance`, `NumCovarianceGroups`, `GetCovarianceGroup`, `HasRelativeCovariance`, `RefreezeCovariance`≈ | 12 in 6 | EData.cpp, EData.h | D | 2/3/19 (sum) | — | W2 `src/EDataCovariance.cpp` (after brief I) |
| `SetKinematicsError`, `EData`, `NumSegments`, `NumTargetEffects`, `IsFit`, `IsErrorAnalysis`, `IsSegmentKey`, `SetFit`, `SetErrorAnalysis`, `Iterate`, `AddSegment`, `AddTargetEffect`, `DeleteLastSegment`, `GetSegment`, `GetSegmentFromKey`, `GetTargetEffect`, `begin`, `end`, `GetSegments` | 85 in 19 | EData.cpp, EData.h | I | 430/35/516 (sum) | — | stays `src/EData.cpp` |
| `ResultsSize` | 1 in 1 | EData.h | R | 0/0/0 (sum) | — | delete (no caller), W1 |
| `NumThmCoherentParams`, `GetThmCoherentParam`, `NumThmGroups`, `GetThmGroup`, `ThmGroupOf`, `IsLastOfThmGroup`, `ThmGroupResiduals`, `ThmBackgroundAt` | 34 in 8 | EData.h, EDataThm.cpp | H | 9/0/18 (sum) | — | stays `src/EDataThm.cpp` (profile) |
| `GetThmCoherentParamOffset` | 1 in 1 | EData.h | H | 0/0/0 (sum) | — | delete (no caller), W1 |
| `Fill` | 150 | EData.cpp 353–502 | S | ≈3/2/22 (AZUREAPI.cpp 4, history_test.cpp 4, AZURE2.cpp 2, AZUREMain.cpp 2) | 2 named, 0 written | stays `src/EData.cpp`; Fill/MakePoints one function after 3d |
| `ApplyPriorCentres` | 28 | EData.cpp 512–539 | V | 1/0/0 | — | after S7: `src/EDataParameters.cpp` |
| `MakePoints` | 151 | EData.cpp 547–697 | S | 2/0/4 (AZUREAPI.cpp 2, AZURE2.cpp 1, AZUREMain.cpp 1) | — | stays `src/EData.cpp`; Fill/MakePoints one function after 3d |
| `Iterations`, `ResetIterations` | 6 in 2 | EData.cpp | I | 0/0/0 (sum) | — | delete (no caller), W1 |
| `ReplaceSubPointGrids` | 228 | EData.cpp 961–1188 | F | 1/0/0 | 4 named, 0 written | Phase 5 (kernels); `src/EDataTargetEffects.cpp` if D1 |
| `SetUpTargetEffects` | 489 | EData.cpp 1190–1678 | F | 2/0/0 | 4 named, 4 written: `screeningEffects_`, `screeningCopies_`, `dopplerEffects_`, `dopplerCopies_` | Phase 5 (kernels); `src/EDataTargetEffects.cpp` if D1 |
| `Initialize` | 131 | EData.cpp 1759–1889 | S/T | ≈0/2/32 (AZUREAPI.cpp 5, history_test.cpp 4, pybind_api.cpp 2, AZURE2.cpp 2) | 3 named, 3 written: `nucleusLabels_`, `channelTable_`, `numNuclei_` | stays `src/EData.cpp`; its phase-shift checks to `EDataChecks.cpp` after 3d |
| `PrintData` | 86 | EData.cpp 1903–1988 | O | 0/0/2 (AZUREAPI.cpp 1, AZUREMain.cpp 1) | — | N2 `src/EDataCheckFiles.cpp` |
| `CalcLegendreP`≈, `CalcCoulombAmplitude`≈ | 22 in 2 | EData.cpp | T | 10/13/0 (sum) | — | stays `src/EData.cpp` (stage tables; Phase 6) |
| `PrintLegendreP` | 37 | EData.cpp 2012–2048 | O | 1/0/0 | — | N2 `src/EDataCheckFiles.cpp` |
| `BuildChannelTable` | 20 | EData.cpp 2058–2077 | T | 2/0/0 | 1 named, 0 written | stays `src/EData.cpp` (stage tables; Phase 6) |
| `CalcEDependentValues` | 39 | EData.cpp 2079–2117 | T | ≈4/5/1 (ScatteringLength.cpp 1) | 1 named, 1 written: `channelTable_` | stays `src/EData.cpp` (stage tables; Phase 6) |
| `PrintEDependentValues` | 56 | EData.cpp 2123–2178 | O | 1/0/0 | — | N2 `src/EDataCheckFiles.cpp` |
| `PrintCoulombAmplitude` | 41 | EData.cpp 2198–2238 | O | 1/0/0 | — | N2 `src/EDataCheckFiles.cpp` |
| `WriteOutputFiles` | 435 | EData.cpp 2259–2693 | O | ≈0/0/7 (AZUREAPI.cpp 2, AZUREAPI.h 1, pybind_api.cpp 1, AZURECalc.cpp 1) | 3 named, 0 written | N1 `src/EDataOutputFiles.cpp` |
| `CountECAmplitudes`, `ECSignaturePath` | 7 in 2 | EData.cpp | X | 4/0/0 (sum) | — | Phase 6; `src/EDataExternalCapture.cpp` if D1 |
| `ECSignature` | 71 | EData.cpp 2862–2932 | X | 1/0/0 | 1 named, 0 written | Phase 6; `src/EDataExternalCapture.cpp` if D1 |
| `CalculateECAmplitudes` | 210 | EData.cpp 2985–3194 | X | ≈3/4/1 (AZUREAPI.cpp 1) | 5 named, 4 written: `ecUsePrevious_`, `ecOutputFile_`, `ecSignature_`, `ecReadPos_` | Phase 6; `src/EDataExternalCapture.cpp` if D1 |
| `InitializeComponentSegments` | 287 | EData.cpp 3201–3487 | C | 1/0/0 | 6 named, 1 written: `channelTable_` | S12 deletes or rewrites |
| `CheckComponentPairs` | 21 | EData.cpp 3502–3522 | C | 2/0/0 | 1 named, 1 written: `componentError_` | S12 deletes or rewrites |
| `FinishTransformedSegment` | 43 | EData.cpp 3529–3571 | C | 2/0/0 | — | S12 deletes or rewrites |
| `CheckTransformedSegments` | 45 | EData.cpp 3580–3624 | C | 2/0/0 | — | S12 deletes or rewrites |
| `CreateComponentSegment` | 21 | EData.cpp 3631–3651 | C | 7/0/0 | 2 named, 2 written: `componentSegments_`, `nextComponentKey_` | S12 deletes or rewrites |
| `MapData` | 106 | EData.cpp 3771–3876 | M/P | 4/0/0 | 1 named, 0 written | S5 (plan; brief I, in work) |
| `MapComponentSegments` | 29 | EData.cpp 3890–3918 | M/P | 1/0/0 | — | S5 (plan; brief I, in work) |
| `ResolveShiftGroups` | 37 | EData.cpp 3960–3996 | V | 1/0/0 | 1 named, 0 written | after S7: `src/EDataParameters.cpp` |
| `ResolveShiftFunctions` | 192 | EData.cpp 4031–4222 | V | 1/0/0 | 1 named, 0 written | after S7: `src/EDataParameters.cpp` |
| `ResolveNormGroups` | 35 | EData.cpp 4248–4282 | V | 1/0/0 | 1 named, 0 written | after S7: `src/EDataParameters.cpp` |
| `CovarianceFileFor` | 18 | EData.cpp 4290–4307 | D | 0/1/0 | 2 named, 2 written: `covarianceError_`, `covFiles_` | W2 `src/EDataCovariance.cpp` (after brief I) |
| `ResolveCovarianceGroups` | 105 | EData.cpp 4316–4420 | D | 1/0/0 | 2 named, 2 written: `covGroups_`, `covFiles_` | W2 `src/EDataCovariance.cpp` (after brief I) |
| `CovarianceChi2Correction` | 34 | EData.cpp 4422–4455 | D | 1/0/3 (AZUREAPI.cpp 2, Objective.cpp 1) | 1 named, 0 written | W2 `src/EDataCovariance.cpp` (after brief I) |
| `WhitenResidualRows` | 31 | EData.cpp 4457–4487 | D | 0/0/4 (AZUREAPI.cpp 2, AZURECalc.cpp 2) | 1 named, 0 written | W2 `src/EDataCovariance.cpp` (after brief I) |
| `WriteCovarianceSummary` | 26 | EData.cpp 4494–4519 | D | 1/0/0 | 1 named, 0 written | W2 `src/EDataCovariance.cpp` (after brief I) |
| `PrepareCovarianceSeeds` | 17 | EData.cpp 4531–4547 | D | 0/0/2 (AZUREAPI.cpp 1, AZURECalc.cpp 1) | 1 named, 0 written | W2 `src/EDataCovariance.cpp` (after brief I) |
| `SyncNormGroupParams` | 27 | EData.cpp 4556–4582 | V | 0/0/5 (AZUREMain.cpp 4, AZURECalc.cpp 1) | 3 named, 0 written | after S7: `src/EDataParameters.cpp` |
| `FillMnParams` | 121 | EData.cpp 4600–4720 | V | ≈0/0/44 (AZURECalc.cpp 10, AZURECalcMCMC.cpp 9, AZUREAPI.cpp 8, history_test.cpp 4) | 6 named, 5 written | after S7: `src/EDataParameters.cpp` |
| `ResolveThicknessGroups` | 46 | EData.cpp 4722–4767 | V | 1/0/0 | 2 named, 1 written: `thicknessSegments_` | after S7: `src/EDataParameters.cpp` |
| `FillThicknessFromParams`, `ThicknessPenalty`, `DopplerPenalty`, `DensityPenalty`, `FillScreeningFromParams`, `ScreeningPenalty`, `FillNormsFromParams`, `FillEnergyShiftsFromParams` | 63 in 8 | EData.cpp | V/E | 10/1/53 (sum) | — | after S7: `src/EDataHandOver.cpp` |
| `FillDopplerFromParams` | 23 | EData.cpp 4805–4827 | V/E | 1/0/0 | 3 named, 0 written | after S7: `src/EDataHandOver.cpp` |
| `SetUpDensityParams` | 55 | EData.cpp 4842–4896 | V | 1/0/0 | 3 named, 3 written: `densityParamOfEffect_`, `densityLeadCopies_`, `densityEffects_` | after S7: `src/EDataParameters.cpp` |
| `FillDensitiesFromParams` | 47 | EData.cpp 4898–4944 | V/E | 1/0/0 | 5 named, 0 written | after S7: `src/EDataHandOver.cpp` |
| `ApplyScreening` | 18 | EData.cpp 4953–4970 | V/E | 2/0/0 | 4 named, 0 written | after S7: `src/EDataHandOver.cpp` |
| `LayoutResults` | 24 | EData.cpp 5040–5063 | R | 3/0/0 | 2 named, 1 written: `resultsSize_` | S5 (plan) |
| `FillEnergyShiftsFromParamsImpl` | 236 | EData.cpp 5078–5313 | V/E | 1/0/0 | — | after S7: `src/EDataHandOver.cpp` |
| `Clone` | 126 | EData.cpp 5347–5472 | M/P | ≈1/0/28 (AZURECalc.cpp 10, AZURECalcMCMC.cpp 6, ReactionRate.cpp 3, park_chi2_fallback_test.cpp 2) | 29 named, 29 written | S5 (plan; brief I, in work) |
| `CalculatePointsInParallel` | 55 | EData.cpp 5517–5571 | M/P | 0/0/4 (AZURECalc.cpp 2, AZUREAPI.cpp 1, history_test.cpp 1) | — | S5 (plan; brief I, in work) |
| `FillThmCoherentFromParams` | 8 in 1 | EDataThm.cpp | H | 1/0/3 (sum) | — | stays `src/EDataThm.cpp` (set-up; after 3d) |
| `BuildThmGroups` | 201 | EDataThm.cpp 451–651 | H | 1/0/0 | 3 named, 3 written: `thmCoherentParamOffset_`, `thmCoherentParams_`, `thmGroups_` | stays `src/EDataThm.cpp` (set-up; after 3d) |
| `ProfileThmGroup` | 24 | EDataThm.cpp 692–715 | H | 1/0/3 (AZUREAPI.cpp 1, AZUREGrad.cpp 1, Objective.cpp 1) | 1 named, 0 written | stays `src/EDataThm.cpp` (profile) |
| `ThmExperimentReports` | 38 | EDataThm.cpp 735–772 | H | 1/0/1 (AZUREAPI.cpp 1) | 2 named, 0 written | stays `src/EDataThm.cpp` (profile) |
| `WriteThmExperiments` | 156 | EDataThm.cpp 774–929 | H/O | 1/0/0 | 1 named, 0 written | N3 `src/EDataThmReports.cpp` |
| `ThmLineshapeTable` | 70 | EDataThm.cpp 931–1000 | H/O | 0/0/1 (AZUREAPI.cpp 1) | 1 named, 0 written | N3 `src/EDataThmReports.cpp` |
| `ThmDistortionTable` | 54 | EDataThm.cpp 1002–1055 | H/O | 0/0/1 (AZUREAPI.cpp 1) | 1 named, 0 written | N3 `src/EDataThmReports.cpp` |
| `ThmVertexTable` | 177 | EDataThm.cpp 1057–1233 | H/O | 0/0/1 (AZUREAPI.cpp 1) | 1 named, 0 written | N3 `src/EDataThmReports.cpp` |
| `Fill`≈, `MakePoints` | 20 in 2 | OldLineSetup.cpp | S | 5/2/26 (sum) | — | 3d deletes (`src/OldLineSetup.cpp`, the old-format path) |

#### `EPoint`

| member(s) | lines | where | class | callers int/data/ext (main ext. files) | state named | destination |
|---|---|---|---|---|---|---|
| `GetThmWeight`, `GetThmLineshape`≈, `GetThmSpectatorWindow`, `GetThmDistortion`≈, `GetThmDwVertex`, `GetThmAngleWindow`, `GetThmCoherent`, `IsTHM`≈, `SetThmLineshape`≈, `SetThmSpectatorWindow`≈, `SetThmDistortion`≈, `SetThmDwVertex`≈, `SetThmCoherent`≈, `SetThmAngleWindow`≈ | 34 in 14 | EPoint.cpp, EPoint.h | H | 22/48/31 (sum) | — | stays `src/EPoint.cpp` (Phase 5 THM observable) |
| `IsAnalyzingPower`, `SetIsAnalyzingPower`, `IsTotal`, `SetIsTotal`, `IsCoefficient`, `GetCoefficientOrder`, `GetCoefficientAbsolute`, `SetCoefficientOrder`, `SetCoefficientAbsolute`, `GetPhaseS`, `IsDifferential`≈, `IsPhase`≈, `IsAngularDist`≈, `IsUPOS`≈, `GetSecondaryDecayL`≈, `GetIc`≈, `GetDelta`≈, `GetL`≈, `GetMaxAngDistOrder`≈, `GetJ`≈ | 40 in 20 | EPoint.cpp, EPoint.h | A | 41/130/356 (sum) | — | stays; kind flags derived from the tree (S9/S10, check 13) |
| `GetAnalyzingPower`, `SetAnalyzingPower`, `GetAngDistA0`, `SetAngDistA0`, `GetFitCrossSection`, `GetFitE1CrossSection`, `GetFitE2CrossSection`, `SetFitCrossSection`, `SetFitE1CrossSection`, `SetFitE2CrossSection`, `ResultSlot`, `ResultAngularCapacity`, `ResultAngularSlots`, `ResultNaturalSlots`, `SetModelRecord`, `GetModelValue`, `SetModelValue`, `GetNumModelAngularDists`, `GetModelAngularDist`, `CopyAngularDistsToModel`, `Rec`, `Rec`, `ModelRecord`, `GetNumAngularDists`, `GetAngularDist`, `SetAngularDists`, `GetNaturalTerms`, `SetNaturalTerms`, `SetResultRecord`, `UsePrivateRecord`, `UsePrivateRecord`, `NoResultRecord`, `NoModelRecord` | 116 in 33 | EPoint.cpp, EPoint.h, ModelPass.cpp | R | 115/25/148 (sum) | — | stays `src/EPoint.cpp` (S5 may renumber) |
| `SetUdrKernel`, `GetUdrKernel`, `GetUdrProjectileMass`, `SetSubGridSpec`, `GetSubGridSpec`, `SetTargetKinematics`, `GetCMFactor`, `IsSubPoint`, `SetScreening`, `GetScreeningParam`, `GetScreeningUe`, `ScreeningFactor`, `GetTargetBlendWeight`, `SetTargetBlendWeight`, `SetTargetFallbackOnly`, `GetTargetFallbackDepth`, `SetTargetFallbackDepth`, `GetTargetLabEnergy`, `SetTargetLabEnergy`, `SetPhotoKinematics`, `IsTargetEffect`≈, `NumSubPoints`≈, `GetTargetEffectNum`≈, `GetStoppingPower`, `GetTargetThickness`, `SetTargetEffectNum`≈, `AddSubPoint`, `SetParentData`, `SetStoppingPower`, `SetTargetThickness`, `GetParentData`, `GetSubPoint`, `GetSubPoints` | 78 in 33 | EPoint.cpp, EPoint.h | F | 133/149/85 (sum) | — | stays `src/EPoint.cpp` until Phase 5 (sub-points become kernels) |
| `GetScreeningK`, `NumThmPsNodes`, `GetThmPsWeight`, `ThmPsTableBytes`, `OwnsEDependentValues`, `SetOwnsEDependentValues`, `GetGeometricalFactor`, `GetSFactorConversion`, `GetSqrtPenetrability`, `GetThmFormFactor`, `GetThmFormFactor`, `GetLoElement`, `GetExpCoulombPhase`, `GetExpHardSpherePhase`, `GetCoulombAmplitude`, `GetCoulombAmplitude`, `SetGeometricalFactor`, `SetSFactorConversion`, `CompoundEnergy`, `MutableEDep`, `MutableLegendre`, `AdoptEDependentValues`, `CopyEDependentValuesFrom`, `AddLoElement`, `AddSqrtPenetrability`, `AddThmFormFactor`, `AddExpCoulombPhase`, `AddExpHardSpherePhase`, `SetCoulombAmplitude` | 141 in 29 | EPoint.cpp, EPoint.h | T | 38/24/108 (sum) | — | D1 (after S5): `src/EPointChannelTables.cpp`; per-stage tables in Phase 6 |
| `IsTargetFallbackOnly`, `ModelSlot`, `ModelAngularSlots`, `StoreSubpointOffsets`, `ApplySubpointShift`, `GetAngleKinFactor`, `SetExcitationEnergy`, `SetAngleKinFactor` | 12 in 8 | EPoint.cpp, EPoint.h | A | 3/0/0 (sum) | — | delete (no caller or no definition), W1 |
| `GetFileRecord`, `SetFileRecord`, `SetDataError`, `GetSegmentKey`≈, `GetEntranceKey`≈, `GetExitKey`≈, `GetLabAngle`, `GetCMAngle`, `GetLabEnergy`, `GetCMEnergy`, `GetExcitationEnergy`, `GetOriginalEnergy`, `GetLabCrossSection`, `GetCMCrossSection`, `GetLabCrossSectionError`, `GetCMCrossSectionError`, `SetLabEnergy`, `SetCMEnergy`, `SetLabAngle`, `SetCMAngle`, `SetExitKey`≈, `SetEntranceKey`≈ | 58 in 22 | EPoint.cpp, EPoint.h | A | 135/425/173 (sum) | — | stays `src/EPoint.cpp`; read from `Measurement`/frame after S7 |
| `GetLabAngleMax`, `GetLabAngleWhy`, `KeepLoadFrame`, `RestoreLoadFrame`, `HasBinWindow`, `SetBinWindow`, `GetBinLowCM`, `GetBinHighCM`, `GetCrossSectionKinFactor`, `ApplyAngleLimit`, `ConvertLabEnergy`, `ConvertExcitationEnergy`, `RelativisticLabToCM`, `ConvertLabAngle`, `ConvertLabAngleNR`, `ConvertLabAngleGammas`, `CalculateCrossSectionConversionFactor`, `ConvertCrossSectionGammas`, `CalculateCrossSectionGammaConversionFactor`, `SetCrossSectionKinFactor` | 101 in 20 | EPoint.cpp, EPoint.h | K | 16/25/1 (sum) | — | S7 (`MeasurementFrame`) |
| `EPoint` | 54 | EPoint.cpp 57–110 | S/A | 5/0/1 (ReactionRate.cpp 1) | 45 named, 45 written | stays `src/EPoint.cpp`; input fields read from `Measurement` after S7 |
| `EPoint` | 43 | EPoint.cpp 118–160 | S/A | 5/0/1 (ReactionRate.cpp 1) | 41 named, 41 written | stays `src/EPoint.cpp`; input fields read from `Measurement` after S7 |
| `EPoint` | 33 | EPoint.cpp 169–201 | S/A | 5/0/1 (ReactionRate.cpp 1) | 30 named, 30 written | stays `src/EPoint.cpp`; input fields read from `Measurement` after S7 |
| `IsMapped`, `NumLocalMappedPoints`, `GetMap`, `SetMap`, `ClearMapping`, `AddLocalMappedPoint`, `ClearLocalMappedPoints`, `GetLocalMappedPoint`, `GetMappedPoints` | 40 in 9 | EPoint.cpp | M | 23/45/19 (sum) | — | S5 deletes (plan; in work) |
| `GetMaxLOrder`, `GetLegendreP`, `AddLegendreP`, `ClearLegendrePolynomials` | 16 in 4 | EPoint.cpp | T | 1/4/16 (sum) | — | D1 (after S5): `src/EPointLegendre.cpp`; angle weights in Phase 5 |
| `GetECAmplitude`, `AddECAmplitude`, `ClearECAmplitudes` | 13 in 3 | EPoint.cpp | X | 4/14/3 (sum) | — | Phase 6; `src/EPointExternalCapture.cpp` if D1 |
| `GetECAmplitudeWithShift` | 49 | EPoint.cpp 632–680 | X | 0/0/3 (AMatrixFunc.cpp 2, RMatrixFunc.cpp 1) | 4 named, 1 written: `ec_shifted_` | Phase 6; `src/EPointExternalCapture.cpp` if D1 |
| `Initialize` | 17 | EPoint.cpp 703–719 | T | ≈0/2/32 (AZUREAPI.cpp 5, history_test.cpp 4, pybind_api.cpp 2, AZURE2.cpp 2) | — | stays `src/EPoint.cpp` |
| `ConvertDecayEnergy` | 5 in 1 | EPoint.cpp | K/B | 0/1/0 (sum) | — | S7 (`MeasurementFrame`, `EnergyIn::DecayLab`) |
| `ConvertLabAngle` | 20 | EPoint.cpp 813–832 | K | 1/2/0 | 2 named, 2 written: `angleKinFactor_`, `cm_angle_` | S7 (`MeasurementFrame`) |
| `ConvertCrossSection` | 20 | EPoint.cpp 867–886 | K | 0/1/0 | 3 named, 3 written: `crossSectionKinFactor_`, `cm_dcrosssection_`, `cm_crosssection_` | S7 (`MeasurementFrame`) |
| `QLabBoost` | 19 | EPoint.cpp 1086–1104 | T | 1/0/0 | — | D1 (after S5): `src/EPointLegendre.cpp`; angle weights in Phase 5 |
| `CalcLegendreP` | 46 | EPoint.cpp 1124–1169 | T | ≈2/8/0 | — | D1 (after S5): `src/EPointLegendre.cpp`; angle weights in Phase 5 |
| `ChannelKeys` | 21 | EPoint.cpp 1199–1219 | T | 1/1/0 | — | D1 (after S5): `src/EPointChannelTables.cpp`; per-stage tables in Phase 6 |
| `CalcEDependentValues` | 244 | EPoint.cpp 1221–1464 | T | ≈3/6/1 (ScatteringLength.cpp 1) | 8 named, 7 written | D1 (after S5): `src/EPointChannelTables.cpp`; per-stage tables in Phase 6 |
| `ClearEDependentRows` | 17 | EPoint.cpp 1469–1485 | T | 3/0/0 | 1 named, 1 written: `edep_` | D1 (after S5): `src/EPointChannelTables.cpp`; per-stage tables in Phase 6 |
| `RecalcEDependentValues` | 37 | EPoint.cpp 1491–1527 | T | 3/2/0 | 2 named, 0 written | D1 (after S5): `src/EPointChannelTables.cpp`; per-stage tables in Phase 6 |
| `ShareEDependentValues` | 16 | EPoint.cpp 1562–1577 | T | 0/2/0 | 8 named, 0 written | D1 (after S5): `src/EPointChannelTables.cpp`; per-stage tables in Phase 6 |
| `CalcCoulombAmplitude` | 45 | EPoint.cpp 1688–1732 | T | ≈3/10/0 | 3 named, 3 written: `spinDependentCoulomb_`, `coulombexchange_`, `coulombdirect_` | D1 (after S5): `src/EPointChannelTables.cpp`; per-stage tables in Phase 6 |
| `CalculateECAmplitudes` | 29 | EPoint.cpp 1825–1853 | X | ≈2/5/1 (AZUREAPI.cpp 1) | — | Phase 6; `src/EPointExternalCapture.cpp` if D1 |
| `AddECAmplitude` | 37 | EPoint.cpp 1859–1895 | X | 3/6/0 | 3 named, 2 written: `ec_amplitudes_`, `ec_energies_` | Phase 6; `src/EPointExternalCapture.cpp` if D1 |
| `MarkFailed` | 15 | EPoint.cpp 1924–1938 | P | 0/2/2 (ModelPass.cpp 2) | 1 named, 0 written | D1 (after S5): `src/EPointCalculate.cpp` |
| `Calculate` | 202 | EPoint.cpp 1940–2141 | P | 9/1/7 (AZUREAPI.cpp 3, ReactionRate.cpp 2, AZURECalcMCMC.cpp 1, ModelPass.cpp 1) | 1 named, 1 written: `integrationPoints_` | D1 (after S5): `src/EPointCalculate.cpp` |
| `CalculateThinTarget` | 60 | EPoint.cpp 2154–2213 | P | 1/0/0 | 2 named, 0 written | D1 (after S5): `src/EPointCalculate.cpp` |
| `IntegrateTargetEffectComponents` | 21 | EPoint.cpp 2223–2243 | P | 1/0/0 | 1 named, 0 written | D1 (after S5): `src/EPointCalculate.cpp` |
| `IntegrateTargetEffectForObservable` | 44 | EPoint.cpp 2261–2304 | P | 4/0/0 | — | D1 (after S5): `src/EPointCalculate.cpp` |
| `IntegrateTargetEffectImpl` | 19 | EPoint.cpp 2453–2471 | F | 0/0/0 | 1 named, 0 written | Phase 5 (kernels); `src/EPointFolds.cpp` if D1 |
| `DopplerWeights` | 71 | EPoint.cpp 2473–2543 | F | 2/0/2 (AZUREGrad.cpp 1, ModelPass.cpp 1) | 1 named, 0 written | Phase 5 (kernels); `src/EPointFolds.cpp` if D1 |
| `DopplerTemperatureDerivative` | 18 | EPoint.cpp 2545–2562 | F | 0/0/3 (AZUREGrad.cpp 3) | 1 named, 0 written | Phase 5 (kernels); `src/EPointFolds.cpp` if D1 |
| `TargetFallbackImpl` | 33 | EPoint.cpp 2564–2596 | F | 0/0/0 | 2 named, 0 written | Phase 5 (kernels); `src/EPointFolds.cpp` if D1 |
| `IntegrateKernelsImpl` | 244 | EPoint.cpp 2598–2841 | F | 0/0/0 | 13 named, 0 written | Phase 5 (kernels); `src/EPointFolds.cpp` if D1 |
| `IntegrateTargetEffect`, `IntegrateSubPointValues`, `TargetEffectWeights` | 20 in 3 | EPoint.cpp | F | 7/0/14 (sum) | — | Phase 5 (kernels); `src/EPointFolds.cpp` if D1 |

#### `ESegment`

| member(s) | lines | where | class | callers int/data/ext (main ext. files) | state named | destination |
|---|---|---|---|---|---|---|
| `ReplaceNormPrior`, `IsNormPriorReplaced`, `SetNominalNorm`, `IsProfiledNorm`, `IsVaryNorm`, `GetNorm`, `GetNominalNorm`, `GetNormError`, `SetNorm` | 21 in 9 | ESegment.cpp, ESegment.h | V | 1/38/47 (sum) | — | stays `src/ESegment.cpp`; norm groups end in S7 |
| `ReplaceShiftPrior`, `IsShiftPriorReplaced`, `ReplaceShiftSqrtPrior`, `IsShiftSqrtPriorReplaced`, `ReplaceShiftFnPrior`, `IsShiftFnPriorReplaced`, `SetNominalEnergyShift`, `GetEnergyShift`, `GetLastEnergyShift`, `GetNominalEnergyShift`, `GetEnergyShiftError`, `IsVaryEnergyShift`, `SetEnergyShift`, `SetLastEnergyShift`, `GetEnergyShiftSqrt`, `GetLastEnergyShiftSqrt`, `GetEnergyShiftSqrtError`, `GetShiftGroup`, `GetShiftLeaderPos`, `IsShiftGroupFollower`, `GetShiftFnGroup`, `SetShiftGroupFollower`, `GetNormGroup`, `GetNormLeaderPos`, `IsNormGroupFollower`, `SetNormGroupFollower`, `IsVaryEnergyShiftSqrt`, `TotalEnergyShift`, `ShiftFnTerm`, `OwnsShiftFn`, `GetDeclaredShiftFn`, `GetNominalShiftFnValues`, `GetShiftFnErrors`, `IsVaryShiftFn`, `HasSqrtShiftBlock`, `HasShiftFn`, `GetShiftFnExpression`, `NumShiftFnParams`, `GetShiftFnSourcePos`, `GetShiftFnParamStart`, `SetShiftFnParamStart`, `GetShiftFnValues`, `GetLastShiftFnValues`, `SetLastShiftFnValues`, `SetShiftFnValue`, `InstallShiftFn`, `CopyShiftFnFrom`, `AddShiftPenalty`, `ShiftPenaltyParts`, `SetEnergyShiftSqrt`, `SetLastEnergyShiftSqrt` | 144 in 51 | ESegment.cpp, ESegment.h | V | 13/159/52 (sum) | — | after S7: `src/ESegmentParameters.cpp` |
| `ReplaceThicknessPrior`, `IsThicknessPriorReplaced`, `GetSampleThickness`, `SetSampleThickness`, `GetNominalThickness`, `HasThicknessFit`, `IsThicknessVaried`, `GetThicknessSigma`, `GetThicknessGroup`, `GetThicknessLeaderPos`, `IsThicknessGroupFollower`, `SetThicknessGroupFollower`, `ThicknessPenalty`≈ | 20 in 13 | ESegment.h | V | 3/27/16 (sum) | — | stays; thickness groups end in S7 |
| `IsAnalyzingPower`, `GetPhaseS`, `IsTotal`, `IsCoefficients`, `IsCoefficientAbsolute`, `IsDifferential`≈, `IsCMDifferential`, `IsPhase`≈, `GetCrossSectionComponent`, `IsAngularDist`≈, `IsUPOS`≈, `GetSecondaryDecayL`≈, `GetIc`≈, `GetDelta`≈, `GetL`≈, `GetMaxAngDistOrder`≈, `GetJ`≈ | 41 in 17 | ESegment.cpp, ESegment.h | A | 21/144/355 (sum) | — | stays; derived from the tree (S9/S10, check 13) |
| `GetCovarianceFile`, `GetCovarianceGroup`, `SetCovarianceGroup`, `GetCovarianceSeed`, `SetCovarianceSeed` | 10 in 5 | ESegment.h | D | 0/7/6 (sum) | — | brief I moves the seeds to `Results` (in work) |
| `GetComponents`, `IsNatural`, `IsTransmission`, `IsYield`, `IsTransformed`, `ClearTransform`, `IsYieldTotalTerm`, `SetYieldTotalTerm`, `NeedsAssembly`, `IsTotalTerm`, `GetNaturalSample`, `SetNaturalSample`, `GetNaturalToken`, `SetNatural`, `IsTotalCapture`, `IsAdvanced`, `GetLegacyOperationType`, `ApplyTotalToPoints`, `SetIsTotalCapture`, `AddComponentSegment`, `SetOperationType`, `GetOperationType`, `HasComponents`, `GetComponentSegments`, `ClearComponents`, `SetComponentScaling`, `GetComponentScaling`, `TermFrameFactor`, `RatioTermFrameFactor`, `Transmission`, `YieldY0`, `HasSubPointIntegration` | 92 in 32 | ESegment.cpp, ESegment.h | C | 17/95/66 (sum) | — | S10/S12 (tree; components as leaves) |
| `SetThmWeight`, `GetThmWeight`, `GetThmLineshape`≈, `GetThmSpectatorWindow`, `GetThmDistortion`≈, `GetThmDwVertex`, `GetThmAngleWindow`, `GetThmCoherent`, `IsTHM`≈, `SetThmLineshape`≈, `SetThmDistortion`≈, `SetThmDwVertex`≈, `SetThmCoherent`≈, `SetThmAngleWindow`≈, `SetThmSpectatorWindow`≈ | 35 in 15 | ESegment.cpp, ESegment.h | H | 7/65/31 (sum) | — | stays `src/ESegment.cpp` |
| `SetReaction` | 40 | ESegment.cpp 31–70 | S | 2/0/0 | 23 named, 23 written | stays `src/ESegment.cpp` (G3 edits it now) |
| `ESegment` | 55 | ESegment.cpp 76–130 | S | 3/0/1 (old_line_specs_test.cpp 1) | 50 named, 50 written | stays `src/ESegment.cpp` (G3 edits it now) |
| `ESegment` | 47 | ESegment.cpp 136–182 | S | 3/0/1 (old_line_specs_test.cpp 1) | 43 named, 42 written | stays `src/ESegment.cpp` (G3 edits it now) |
| `IsInSegment`, `CMEnergyFactor`, `ConvertPointEnergyToCM` | 22 in 3 | ESegment.cpp | S/K | 4/4/0 (sum) | — | S7 (`Measurement`, `frame`) |
| `ProfileNormChiSquared` | 24 | ESegment.cpp 256–279 | H | 0/1/3 (AZUREAPI.cpp 2, Objective.cpp 1) | 1 named, 0 written | stays `src/ESegment.cpp` |
| `IsTargetEffect`≈, `NumPoints`, `GetEntranceKey`≈, `GetExitKey`≈, `GetTargetEffectNum`≈, `GetSegmentKey`, `GetMinEnergy`, `GetMaxEnergy`, `GetMinAngle`, `GetMaxAngle`, `GetDataFile`, `AddPoint`, `SetTargetEffectNum`≈, `SetSegmentKey`, `SetExitKey`≈, `SetEntranceKey`≈, `SetMinAngle`, `SetMaxAngle`, `GetPoint`, `GetPoints` | 70 in 20 | ESegment.cpp | I/A | 49/476/193 (sum) | — | stays `src/ESegment.cpp` |
| `Fill` | 93 | ESegment.cpp 389–481 | S/K | ≈0/6/22 (AZUREAPI.cpp 4, history_test.cpp 4, AZURE2.cpp 2, AZUREMain.cpp 2) | 1 named, 0 written | S7 (`Measurement`, `frame`) |
| `ConvertPointToCM` | 23 | ESegment.cpp 495–517 | S/K | 2/2/0 | — | S7 (`Measurement`, `frame`) |
| `GetNominalEnergyShiftSqrt`, `HasShiftGroupFollowers`, `SetHasShiftGroupFollowers`, `HasNormGroupFollowers`, `SetHasNormGroupFollowers`, `SetVaryNorm`, `AddComponent` | 14 in 7 | ESegment.cpp | V | 0/2/0 (sum) | — | delete (no reader), W1 |
| `UpdatePointEnergiesWithShift`, `ShiftState`, `BeginShiftRefill`, `ShiftRowSource`, `ShiftHost` | 34 in 5 | ESegment.cpp | E | 5/5/3 (sum) | — | after S7: `src/ESegmentShifts.cpp` |
| `CommitShiftTable` | 15 | ESegment.cpp 934–948 | E | 0/1/0 | 2 named, 2 written: `shiftTables_`, `shiftRows_` | after S7: `src/ESegmentShifts.cpp` |
| `UpdatePointGroupEnergyWithShift` | 15 | ESegment.cpp 983–997 | E | 0/2/1 (AZURECalc.cpp 1) | 1 named, 0 written | after S7: `src/ESegmentShifts.cpp` |
| `UpdatePointEnergyWithShift` | 145 | ESegment.cpp 1006–1150 | E | 3/0/0 | — | after S7: `src/ESegmentShifts.cpp` |
| `TransformTerms` | 24 | ESegment.cpp 1393–1416 | C | 0/0/2 (AZUREGrad.cpp 1, ModelPass.cpp 1) | 3 named, 0 written | S10/S12 (tree; components as leaves) |

### 2.5 Data members

| class | EData data members (header line) | count | destination |
|---|---|---|---|
| I containers and run state | `targetEffects_` 466, `segments_` 467, `kinematicsError_` 477, `iterations_` 478, `isFit_` 497, `isErrorAnalysis_` 498, `nucleusLabels_` 513, `numNuclei_` 514 | 8 | stay |
| C composition | `componentSegments_` 468, `nextComponentKey_` 472, `componentError_` 473 | 3 | S12 (components as leaves); `componentError_` is write-only: delete in W1 |
| D data covariance | `covFiles_` 474, `covGroups_` 475, `covarianceError_` 476 | 3 | stay (W2 moves the functions) |
| V parameter layout and nuisances | `normParamOffset_` 479, `energyShiftParamOffset_` 480, `energyShiftSqrtParamOffset_` 481, `energyShiftFnParamOffset_` 482, `screeningParamOffset_` 483, `screeningEffects_` 484, `screeningCopies_` 485, `densityParamOffset_` 486, `densityEffects_` 487, `densityLeadCopies_` 488, `densityParamOfEffect_` 489, `dopplerParamOffset_` 490, `dopplerEffects_` 491, `dopplerCopies_` 492, `thicknessParamOffset_` 493, `thicknessSegments_` 494, `priorCov_` 512 | 17 | S7 / Phase 2 `inference/` |
| X external-capture file state | `ecReadPos_` 499, `ecUsePrevious_` 502, `ecSignature_` 503, `ecOutputFile_` 504 | 4 | Phase 6 |
| H THM experiments | `thmGroups_` 505, `thmCoherentParams_` 506, `thmCoherentParamOffset_` 507 | 3 | stay |
| R result layout | `resultsSize_` 516, `modelSize_` 516 | 2 | S5 (plan) |
| F sub-point grids | `gridNuc_` 537, `nodeCache_` 539 | 2 | Phase 5/6 (`ARCHITECTURE.md` §12) |
| T stage channel table | `channelTable_` 542 | 1 | Phase 6 (per-stage tables) |

| class | EPoint data members (header line) | count | destination |
|---|---|---|---|
| A input datum | `fileRecord_` 754, `lab_angle_` 761, `original_energy_` 762, `lab_energy_` 764, `lab_crosssection_` 770, `lab_dcrosssection_` 771, `bin_low_lab_` 792, `bin_high_lab_` 793, `coefficient_order_` 813, `coefficient_absolute_` 814 | 10 | S7 (`Measurement`) |
| A segment settings copied into the point | `is_differential_` 731, `is_phase_` 732, `is_thm_` 733, `thm_weight_` 734, `thm_lineshape_` 735, `thm_window_` 736, `thm_distortion_` 737, `thm_dw_` 738, `thm_angle_` 739, `thm_coherent_` 740, `is_ang_dist_` 742, `is_analyzing_power_` 743, `entrance_key_` 749, `exit_key_` 750, `segment_key_` 751, `l_value_` 752, `max_ang_dist_order_` 759, `j_value_` 777, `phase_s_value_` 778, `isUPOS_` 798, `secondaryDecayL_` 799, `Ic_` 800, `delta_` 801, `is_total_` 857, `crossSectionComponent_` 876 | 25 | S9/S10 (tree), S7 |
| K kinematics derived from the datum | `cm_angle_` 760, `cm_energy_` 763, `excitation_energy_` 765, `cm_crosssection_` 768, `cm_dcrosssection_` 769, `angleKinFactor_` 782, `labAngleMax_` 783, `labAngleWhy_` 784, `crossSectionKinFactor_` 785, `bin_low_cm_` 794, `bin_high_cm_` 795, `loadFrame_` 788 | 12 | S7 (`MeasurementFrame`); `angleKinFactor_` is read only by a dead getter: W1 |
| T energy-dependent tables | `geofactor_` 772, `sfactorconv_` 773, `screen_k_` 774, `coulombamplitude_` 803, `spinDependentCoulomb_` 806, `coulombdirect_` 807, `coulombexchange_` 808, `legendreP_` 812, `edep_` 823, `thm_ps_` 850, `eDependentEnergy_` 854, `eDependentValid_` 855, `own_e_dependent_` 856 | 13 | Phase 5/6 (grid tables) |
| X external capture | `ec_amplitudes_` 858, `ec_energies_` 859, `ec_shifted_` 868 | 3 | Phase 6 |
| F apparatus | `is_sub_point_` 744, `cm_factor_` 745, `bohr_z1_` 746, `bohr_z2_` 747, `sub_grid_spec_` 748, `targetEffectNum_` 753, `target_blend_weight_` 755, `target_fallback_only_` 756, `target_fallback_depth_` 758, `udr_kernel_` 766, `udr_mass_` 767, `screen_ue_` 775, `screen_param_` 776, `stoppingPower_` 779, `targetThickness_` 780, `targetLabEnergy_` 781, `photo_qgamma_` 796, `photo_mcn_` 797, `integrationPoints_` 870, `parentData_` 871 | 20 | Phase 5 (kernels) |
| F scoped flags set during an evaluation | `doppler_applied_` 730, `target_fallback_now_` 757 | 2 | Phase 5 / 6 (no hidden state) |
| M mapping | `is_mapped_` 741, `energy_map_` 802, `local_mapped_points_` 869 | 3 | S5 (plan; in work) |
| R result record offsets | `slot_` 878, `num_ang_slots_` 879, `num_natural_slots_` 880, `num_model_ang_` 881, `model_slot_` 896 | 5 | stay |

| class | ESegment data members (header line) | count | destination |
|---|---|---|---|
| A segment line settings and kind flags | `isdifferential_` 501, `iscmdifferential_` 502, `isphase_` 503, `isTHM_` 504, `phaseS_` 512, `crossSectionComponent_` 514, `isTargetEffect_` 515, `isAngDist_` 517, `isAnalyzingPower_` 518, `isUPOS_` 519, `secondaryDecayL_` 520, `Ic_` 521, `delta_` 522, `entrancekey_` 524, `exitkey_` 525, `l_` 526, `targetEffectNum_` 527, `segmentKey_` 528, `maxAngDistOrder_` 529, `min_e_` 530, `max_e_` 531, `min_a_` 532, `max_a_` 533, `j_` 534, `isTotal_` 579, `isCoefficients_` 589, `coefficientsAbsolute_` 590, `datafile_` 592 | 28 | S9/S10 (tree), S7 (`Measurement`) |
| H THM attachments | `thmWeight_` 505, `thmLineshape_` 506, `thmWindow_` 507, `thmDistortion_` 508, `thmDwVertex_` 509, `thmAngle_` 510, `thmCoherent_` 511 | 7 | stay |
| V normalization | `varyNorm_` 516, `dataNorm_` 535, `dataNormNominal_` 536, `dataNormError_` 537, `normGroup_` 552, `normLeaderPos_` 553, `hasNormGroupFollowers_` 554, `normPriorReplaced_` 564 | 8 | S7; `hasNormGroupFollowers_` has no reader: W1 |
| V energy shift parameters | `energyShift_` 538, `lastEnergyShift_` 539, `energyShiftNominal_` 540, `energyShiftError_` 541, `varyEnergyShift_` 542, `energyShiftSqrt_` 543, `lastEnergyShiftSqrt_` 544, `energyShiftSqrtNominal_` 545, `energyShiftSqrtError_` 546, `varyEnergyShiftSqrt_` 547, `shiftGroup_` 548, `shiftLeaderPos_` 549, `hasShiftGroupFollowers_` 550, `shiftFnGroup_` 551, `hasSqrtShiftBlock_` 558, `shiftFnDeclared_` 560, `shiftFnNominal_` 561, `shiftFnError_` 562, `shiftFnVary_` 563, `shiftPriorReplaced_` 564, `shiftSqrtPriorReplaced_` 564, `shiftFnPriorReplaced_` 566, `hasShiftFn_` 568, `shiftFn_` 569, `shiftFnValues_` 570, `lastShiftFnValues_` 571, `shiftFnConfig_` 572, `shiftFnSourcePos_` 573, `shiftFnParamStart_` 574 | 29 | after S7 `ESegmentParameters.cpp`; `hasShiftGroupFollowers_`, `energyShiftSqrtNominal_` read only by dead getters: W1 |
| E channel tables of the last shift states | `shiftTables_` 497, `shiftRows_` 498 | 2 | after S7 `ESegmentShifts.cpp` |
| D data covariance | `covarianceFile_` 555, `covGroup_` 556, `covSeed_` 557 | 3 | `covSeed_` to `Results` (brief I, in work) |
| C composition, sample, transform | `isTotalCapture_` 523, `thicknessPriorReplaced_` 565, `isAdvanced_` 575, `operationType_` 576, `naturalToken_` 577, `naturalSample_` 578, `transmission_` 580, `yield_` 581, `thickness_` 582, `thicknessFit_` 583, `thicknessVary_` 584, `thicknessSigma_` 585, `thicknessGroup_` 586, `thicknessLeaderPos_` 587, `yieldTotalTerm_` 588, `components_` 591, `componentSegments_` 596, `segmentOperationType_` 597, `componentScaling_` 598, `componentCalculationMutex_` 599 | 20 | S10/S12; `componentCalculationMutex_` is never locked: W1 |
| I points | `points_` 593 | 1 | stay (S7: inputs shared) |

Members declared `mutable`: `EPoint::labAngleMax_`, `labAngleWhy_` (written by the `const` `ApplyAngleLimit`),
`ec_shifted_` (written by the `const` `GetECAmplitudeWithShift`), `ESegment::componentCalculationMutex_` (never locked). Members that point at another
object: `EPoint::parentData_` (the `EData` that holds the point; every fold reaches its target effect through it),
`local_mapped_points_` (S5 removes it), `ESegment::componentSegments_` (into `EData::componentSegments_`, a `std::deque`
so that the pointers stay valid), `shiftFnConfig_` (a `Config`).

### 2.6 Dead code and duplicates

Dead (no caller found, or no reader of what is written):
- `EData`: `GetScreeningParamOffset` (`EData.h` 96), `GetDensityParamOffset` (119), `DensityParamOfEffect` (121–123),
  `GetDopplerParamOffset` (143), `GetThicknessParamOffset` (160), `ResultsSize` (311), `GetThmCoherentParamOffset`
  (390), `Iterations` (`EData.cpp` 703–705), `ResetIterations` (1750–1752); `componentError_` (`EData.h` 473) is written
  at `EData.cpp` 3507 and 3517 and never read.
- `EPoint`: `StoreSubpointOffsets`, `ApplySubpointShift` (`EPoint.h` 716–717: declared, defined nowhere),
  `IsTargetFallbackOnly` (334), `ModelSlot` (551), `ModelAngularSlots` (553), `GetAngleKinFactor` (`EPoint.cpp`
  545–547; with it `angleKinFactor_`, written at 108, 158, 824, 829, 857 and by `SetAngleKinFactor` 1026–1028, read only
  by the getter), `SetExcitationEnergy` (987–989).
- `ESegment`: `GetNominalEnergyShiftSqrt` (723–725; with it `energyShiftSqrtNominal_`), `HasShiftGroupFollowers` (745)
  and `hasShiftGroupFollowers_`, `HasNormGroupFollowers` (763) and `hasNormGroupFollowers_` (their setters are called at
  `EData.cpp` 3974 and 4262), `SetVaryNorm` (1235–1237), `AddComponent` (1275–1278), `componentCalculationMutex_`
  (`ESegment.h` 599: made at `ESegment.cpp` 129 and 181, never locked).
- Unreachable: `EData::FillEnergyShiftsFromParamsImpl` lines 5261–5312 (52 lines), the branch for `data == nullptr`;
  all 24 calls of `FillEnergyShiftsFromParams` in `src/`, `api/src/` and `tests/` pass the data set itself.
- A loop that never runs: `EData.cpp` 3452–3458 (`mappedPoint < point->GetMappedPoints().begin()`); the S5 branch
  deletes it (in work).

Two implementations of one thing:
1. The sub-point grid of a point is laid out twice in `SetUpTargetEffects`: for segments (`EData.cpp` 1288–1464) and for
   component segments (1465–1646, 182 lines); the two differ in the variable names, two messages, and the Doppler and
   β-delayed refusals (1264–1278), which only the first has. A sub-point is built in two more places with the same
   four statements: `AddTargetFallbackNodes` (768–774) and the `layOut` lambda of `ReplaceSubPointGrids` (1083–1092).
2. `EData::Fill` (353–502) and `EData::MakePoints` (547–697) repeat the segment set-up (refusals, the total-capture
   pair check, the THM weight, the components); `diff` of the two bodies, indentation ignored, has 109 differing lines
   of 301 (both sides counted).
3. The external-capture amplitudes are computed, written and read back by one loop for the segments
   (`CalculateECAmplitudes`, 2985–3194) and by a second one for the component segments (in
   `InitializeComponentSegments`, 3345–3471); `WalkECAmplitudes` (2741–2817) enumerates the same layout a third
   time to count and hash it (by design: its comment says it mirrors the first).
4. The stage tables are run twice: `EData::CalcEDependentValues`, `CalcLegendreP`, `CalcCoulombAmplitude` (1994–2192)
   for segments, and the same three passes inside `InitializeComponentSegments` for components.
5. A base segment's shift is copied into its components by two blocks with different fields:
   `FillEnergyShiftsFromParamsImpl` 5156–5166 (SUM/RATIO: `SetEnergyShift`, `SetEnergyShiftSqrt`, `CopyShiftFnFrom`,
   no `SetLast*`) and 5169–5200 (total capture: also `SetLast*`); the unreachable branch holds two more copies
   (5277–5305).
6. The energy-shift prior is written twice in `ESegment`: `AddShiftPenalty` (851–860) and `ShiftPenaltyParts` (862–872);
   a third time as indices in `include/Objective.h` 110–125 (`ForEachShiftPenalty`).
7. The normalization prior width n₀σ/100 is formed in four places: `src/Objective.cpp` 62, `include/Objective.h` 191,
   `src/AZURECalc.cpp` 882, and the writer, `EData.cpp` 2461.
8. The bare probe of a target-effect point (copy, private record, effect 0, no mapped points, no sub-points, calculate)
   is written three times in `EPoint::Calculate`: 2039–2047, 2094–2102, 2129–2138.
9. `ESegment` holds its operation type twice: `operationType_` (`ESegment.h` 576, read by `GetLegacyOperationType`) and
   `segmentOperationType_` (597, read by `GetOperationType`).
10. `EPoint::IntegrateTargetEffectComponents` (2223–2243) writes the E1 and E2 values into the sub-points' σ slots,
    integrates and restores them — the pattern step S4 removed from the transforms — while `IntegrateSubPointValues`
    (2847–2849) integrates given values without writing.
11. The particle Jacobian with two relativistic fall-backs: `ConvertCrossSection` (867–886) keeps the load-time factor,
    `CalculateCrossSectionConversionFactor` (894–902) falls back to the non-relativistic one (known: HANDOFF, director
    entry of 2026-11-10, open findings).

## 3. Target files

All in `src/` (engine adapter code, `data_model.md` §2.5), named after the class whose members they define and their
job, as the `CNuc` split did (`CNucReport.cpp`, `CNucChecks.cpp`); free functions on plain records go to `experiment/`
when they need no engine header. Sizes are the spans counted in §2.3 at 070756a1, before the code later steps add, and
without the `#include` block of a new file (`src/CNucReport.cpp` has 33 lines before its first definition).

| file | job | taken from (lines) | size | includes | form | when |
|---|---|---|---|---|---|---|
| `src/EDataOutputFiles.cpp` | the run's output files | `EData.cpp` 231–265, 2239–2693 | 490 | `EData.h`, `AZUREOutput.h`, `CovarianceBand.h`, `CNuc.h`, `Config.h` | members (read `segments_`, `nucleusLabels_`, `thmGroups_`) | N1, N4 |
| `src/EDataCheckFiles.cpp` | the four data check files | `EData.cpp` 1898–1988, 2007–2048, 2118–2178, 2193–2238 | 240 | `EData.h`, `CNuc.h`, `Config.h` | members | N2 |
| `src/EDataThmReports.cpp` | `thm_experiments.out` and the THM tables of the API | `EDataThm.cpp` 773–1233 | 461 | `EData.h`, `ThmReports.h`, `ThmLineshape.h`, `ThmDistortion.h`, `ThmDwVertex.h`, `ThmVertexBoundary.h`, `CNuc.h` | members | N3 |
| `src/EDataThm.cpp` (stays) | THM set-up and profile | — | 772 | as now | members | after N3 |
| `src/EDataCovariance.cpp` | data covariance groups, χ² correction, whitening, summary | `EData.cpp` 4289–4547 | 259 | `EData.h`, `DataCovariance.h` | members (`covGroups_`, `covFiles_`) | W2 |
| `src/EDataChecks.cpp` | load-time refusals and checks that need the level scheme | `EData.cpp` 44–230, 266–325 less what 3b/3d puts in `Build`; `Initialize` 1764–1848 | ≤ 332 | `EData.h`, `CNuc.h`, `Config.h` | free functions on `const ESegment&`, `CNuc&` (they are today) | after 3d |
| `src/EDataParameters.cpp` | layout and names of the data-side parameters, groups, prior centres, nuisance blocks | `EData.cpp` 503–539, 3927–4288, 4548–4776, 4797–4803, 4835–4896 | ≤ 697 | `EData.h`, `ParameterKind.h`, Minuit2 | members | after S7 (`inference/` in Phase 2's remainder) |
| `src/EDataHandOver.cpp` | a parameter vector into segments, effects and points | `EData.cpp` 4777–4796, 4804–4834, 4897–4987, 4997–5024, 5077–5313 less 5261–5312 | ≤ 355 | `EData.h`, `CNuc.h` | members (they write the points) | after S7 |
| `src/EData.cpp` (stays) | the data set: construction from specs, containers, flags, the set-up sequence, the stage tables | — | 722 (+ what S5 leaves of `Clone` and the parallel pass, ≤ 193) | | members | end of Phase 4 |
| `src/EDataTargetEffects.cpp` | target effects, sub-point grids and their re-placement | `EData.cpp` 741–1678 less 1465–1646 | 756 | `EData.h`, `TargetEffect.h`, `AdaptiveIntegrationGrid.h`, `NodeTableCache.h` | members (`gridNuc_`, `nodeCache_`) | D1, after S12; Phase 5 replaces it |
| `src/EDataExternalCapture.cpp` | external-capture amplitudes, `intEC` signature and I/O | `EData.cpp` 2694–3194 | 501 | `EData.h`, `ECTable.h`, `ECAmplitudeCache.h` | members (`ec*_`) | D1; Phase 6 replaces it |
| `src/EPoint.cpp` (stays) | construction, accessors, setters, result records | — | 823 (S7 removes the input and frame accessors) | | members | end of Phase 4 |
| `src/EPointCalculate.cpp` | `Calculate`, thin target, observable integration, failure marking | `EPoint.cpp` 1911–2304 | 394 less S5's deletions | `EPoint.h`, `AMatrixFunc.h`, `RMatrixFunc.h`, `THMMatrixFunc.h`, `TargetEffect.h` | members | D1, after S5 |
| `src/EPointChannelTables.cpp` | channel rows, their sharing, Coulomb amplitude, table getters | `EPoint.cpp` 468–506, 556–602, 1170–1740 | 657 | `EPoint.h`, `ChannelTable.h`, `CNuc.h` | members (`edep_`, `MutableEDep`) | D1; Phase 6 replaces it |
| `src/EPointLegendre.cpp` | Legendre polynomials with Q coefficients | `EPoint.cpp` 1037–1169 | 133 | `EPoint.h`, `TargetEffect.h` | file-local helpers + members | D1 |
| `src/EPointExternalCapture.cpp` | the point's external-capture amplitudes | `EPoint.cpp` 603–680, 1741–1910 | 248 | `EPoint.h`, `ECIntegral.h`, `ECTable.h` | members | D1; Phase 6 replaces it |
| `src/EPointFolds.cpp` | the apparatus kernels on sub-point values | `EPoint.cpp` 2400–2862 | 463 | `EPoint.h`, `TELinearForm.h`, `TargetEffect.h` | members today; free functions of arrays as the first step of Phase 5 | D1; Phase 5 replaces it |
| `include/TELinearForm.h` | the linear form the folds are evaluated on | `EPoint.h` 66–105 | 40 | `<vector>` | class | with `EPointFolds.cpp` |
| `include/EPointTables.h` | `EDependentTables`, `EDependentShare`, `ThmPsTable` as namespace-level types, aliased in `EPoint` | `EPoint.h` 137–206, 828–849 | 92 | `SharedTable.h`, `Constants.h` | plain records | with `EPointChannelTables.cpp` |
| `src/ESegment.cpp` (stays) | construction from specs, settings, THM attachments, profiled norm | — | 641 | | members | end of Phase 4 |
| `src/ESegmentParameters.cpp` | norm-group and energy-shift parameters, shift function, shift priors | `ESegment.cpp` 632–670, 687–880 | 217 (233 less W1) | `ESegment.h`, `Equation.h` | members | after S7 |
| `src/ESegmentShifts.cpp` | re-energising the points after a shift, refill tables | `ESegment.cpp` 881–1150 | 270 | `ESegment.h`, `ChannelTable.h`, `CNuc.h` | members | after S7 |

What stays a member and why: everything that reads or writes private state that a planned step removes (the writer
reads `segments_`, `nucleusLabels_`, `thmGroups_`; the hand-over writes into points; the target-effect set-up owns
`gridNuc_` and `nodeCache_`). Making them free functions on `const EData&` needs `const` overloads of `GetSegments`,
`GetSegment`, `GetTargetEffect`, `NormUsed`, `ThmGroupOf`, `begin`/`end` (all non-`const` today) for code that S5, S7
and Phases 5–6 change anyway. Free functions on plain records, preferred where the inputs allow: the load-time refusals
(already free functions of a segment and the nucleus; with step 3d they become `Build` checks in `experiment/` or stay
in `EDataChecks.cpp`), the Q-coefficient helpers (`GaussLegendreRule`, `QTransformRule`, `LegendreSeries`, already
file-local), the folds (the first step of Phase 5: functions of sub-point energies, stopping powers and values that the
kernel builders reproduce), `ShiftedEnergy`, `ECSignatureHash`. No new mutable `static` or `thread_local`, no
registration from static objects: a moved function keeps the statics it has (listed in §6), none is added.

What is left of the three classes at the end of Phase 4 (after N, W, S5, S7, S10–S12, 3d and D1): `EData` is the
measurement set with its set-up sequence, stage tables, hand-over, parameters and outputs in eleven files, none above
772 lines but `EData.cpp` (722 plus what S5 leaves of the 193 lines of `Clone` and the parallel pass); `ESegment` is a segment's settings and shift state in three files; `EPoint` is the calculation unit in six
files, `EPoint.cpp` the largest at 823 lines before S7 takes the datum out. Phase 5 deletes the folds, the sub-point
grids and `EDataTargetEffects.cpp` (kernels on the energy grid); Phase 6 deletes the per-point tables, the
external-capture amplitudes per point, `Clone`, `nodeCache_` and `gridNuc_` (`ARCHITECTURE.md` §12).

## 4. What waits for what

### 4.1 The branches in work

Checked at 2026-11-10 10:42 with `git -C <wt> diff --stat 070756a1` (the three branches start from cdf58d41; the
difference to 070756a1 in `include/Constants.h`, `src/ECIntegral.cpp` and one line of `src/EPoint.cpp` is 070756a1's own
commits, not theirs). They are unfinished and unreviewed; what they change is "in work".

| branch | brief | edits in the data classes (hunks against cdf58d41) |
|---|---|---|
| `wt-expr` (H: S10, S11) | `ModelPass.cpp`, `Predictions.cpp` (new), `experiment/Expression*`, later `AZUREGrad.cpp` | `EData.h` 238 (`BuildPredictions`); `ESegment.h` 10, 269, 591 (`reaction_`, `prediction_` replace `components_`); `EData.cpp` 497, 692 (`Fill`, `MakePoints`); `ESegment.cpp` 64 |
| `wt-plan` (I: follow-ups, S5) | `CalculationPlan.{h,cpp}` (new), `AZUREAPI.cpp` | `EData.h` 3, 240–243, 326, 468; `EPoint.h` 318, 555, 640, 712, 715, 869; `ESegment.h` 212–220, 557 (covariance seeds out); `EData.cpp` 1043–1158 (`ReplaceSubPointGrids`), 3452 (dead loop), 3773–3914 (mapping), 4442 (seeds), 5132–5232 (hand-over), 5372–5569 (`Clone`, parallel pass); `EPoint.cpp` 326, 1444, 1926–2202 (`Calculate`), 2329–2343 (mapping), 2965–3006; `ESegment.cpp` 970–992 |
| `wt-fmt3` (G3: format 3 runs) | `Config.cpp`, `CNuc`, `experiment/Build*`, `convert/` | `EData.h` 83–153, 296–300; `ESegment.h` 392, 586 (`parameterNames_`); `EData.cpp` 4588–4759 (`SegmentParameterName`, `FillMnParams`); `ESegment.cpp` 76, 136 (constructors) |

All three also edit `src/CMakeLists.txt` (line 10, and 92–100) and `tests/consistency/repo_consistency_test.py`. No
branch edits `src/EDataThm.cpp`, `EData.cpp` 231–265 or 1898–2693, which steps N1–N3 take.

### 4.2 Parts that a planned step deletes or rewrites

Moving these before that step is waste and a source of merge conflicts.

| part (lines counted) | deleted or rewritten by |
|---|---|
| mapping, result layout, `Clone`, parallel pass (`EData.cpp` 267 + 52 + 193; `EPoint.cpp` 55; `ScopedPointMappingOff`; `ESegment::ShiftHost`, `UpdatePointGroupEnergyWithShift`) | S5 (brief I, in work) |
| covariance seeds (`ESegment.h` 212–220, `EData.cpp` 4422–4547) | brief I part 1 (in work) |
| the composition in `ESegment` (148) and `ModelPass.cpp` | S10, S11 (brief H, in work) |
| the input datum and its frame (`EPoint.cpp` kinematics 209, the constructors' input part, `ESegment::Fill` and conversion 136), norm and shift groups with their follower parameters (`EData.cpp` 3960–4288, 4556–4582), the hand-over of shifts (`FillEnergyShiftsFromParamsImpl`) | S7 (`Measurement`, `MeasurementFrame`) |
| component segments (`EData.cpp` 3195–3651, 457; the component grid block 1465–1646, 182; `ESegment::componentSegments_`) | S12 |
| the old-format entry (`src/OldLineSetup.cpp`), the load-time refusals' second home, `FillMnParams`' names | 3b/3d (G3 in work) |
| target effects and grids (938 less 182), folds (`EPoint.cpp` 463), sub-points, Q-coefficient angle weights | Phase 5 (energy grid, kernels) |
| per-point channel tables (571 + 86), external capture (`EData.cpp` 501, `EPoint.cpp` 248), `Clone`, pools, `nodeCache_`, `gridNuc_`, mutable statics | Phase 6 |

Stable now and edited by nobody: the output writer and the check files of `EData`, and the THM reports and tables of
`EDataThm.cpp`. Those are the first steps.

## 5. Ordered steps

Each step is one commit (N4: two) on its own branch, mergeable alone, with the evidence of §1. "Sizes after" are counted
at 070756a1; a step merged after other branches applies to their sizes the same differences.

### 5.1 Now, in parallel with H, I, G3

**N1 — the output writer to `src/EDataOutputFiles.cpp`.** Moves `WriteNaturalDecomposition` (`EData.cpp` 231–265),
`WriteBandLine` and `WriteOutputFiles` (2239–2693): 490 lines, no other file-local helper is used by them (checked).
`EData.cpp` 5,571 → 5,081. Check 11 (one objective) allows `ScreeningPenalty`, `DensityPenalty`, `DopplerPenalty` and
`ThicknessPenalty` only in `DEFINERS`, which lists `src/EData.cpp`: the step adds four `ALLOWED_PENALTY_USES` rows for
the new file with the reason (`chiSquared.out` prints the terms `Results::objective` does not hold; decision D2). The
179 × 3 comparison reads every `AZUREOut_*`, `chiSquared.out` and `normalizations.out` the menu-1 runs write; not
reached by it and compared byte for byte, base against branch: the `.band` files (`tests/thm_band`,
`tests/prior_covariance`, `tests/fixed_param_par` as their `check.sh` run them), the lab-frame output of `kinflag.dat = 1`
(`tests/gamma_lab_frame`, `tests/component_subpoint_angle`), and `shifts.out` of one project with a varied shift and a shift function
(`tests/shift_fn`, `tests/energy_shift_sqrt`).

**N2 — the check files to `src/EDataCheckFiles.cpp`.** Moves `PrintData` (1898–1988), `PrintLegendreP` (2007–2048),
`PrintEDependentValues` (2118–2178), `PrintCoulombAmplitude` (2193–2238): 240 lines. `EData.cpp` 5,081 → 4,841. No
project under `tests/` or `examples/` sets the four data checks (lines 6–9 of `<config>`) to `file` or `screen` (0 of
the 189 `.azr` files, counted), so the 179 runs do not reach this code: T1, T3, 13N, `natural_neon` and `thm_lineshape` are
run base against branch with the four checks `file`, and once `screen`, and the check files and the screen output
compared byte for byte. To add: a ctest that runs T1 with the four checks on and fails when a file is missing or its
line count changes (`tests/check_files_data/`), so these writers are reached by the suite from now on.

**N3 — the THM reports and tables to `src/EDataThmReports.cpp`.** Moves `WriteThmExperiments` (`EDataThm.cpp`
773–929), `ThmLineshapeTable` (930–1000), `ThmDistortionTable` (1001–1055), `ThmVertexTable` (1056–1233): 461 lines.
`ThmGroupEnergyRange` (`EDataThm.cpp` 52–61, file-local) is used by the set-up (280, 325) and by the writer (857, 880,
904): it becomes a private member of `EData` (one declaration in `EData.h`, no arithmetic). `ThmExperimentReports`
stays with the profile, which shares `GatherThmGroup` with it. `EDataThm.cpp` 1,233 → 772. Checks: 179 × 3
(`thm_experiments.out` of the THM projects); ctest `pyazr_thm_lineshape`, `pyazr_thm_distortion`,
`pyazr_thm_vertex_reach`, `pyazr_thm_dw_vertex`, `pyazr_thm_experiment` and `tests/thm_dw_vertex/check.sh`; the three tables of `tests/thm_lineshape`,
`tests/thm_distortion`, `tests/thm_dw_vertex` printed through pyazr with `repr` from both builds and compared.

**N4 — `WriteOutputFiles` in four functions** (split in place, in `EDataOutputFiles.cpp`): the per-segment files with
their band and sample files (today 2259–2491), `chiSquared.out` (the per-segment rows now written inside the segment
loop at 2453–2478, and the totals 2494–2568), `normalizations.out` (2571–2604), `shifts.out` (2607–2692); the THM
writer is called where it is now. The per-segment χ² rows move to a loop of their own over the same segments in the
same order, so the sums and the file are the same. Checks: those of N1. Size unchanged (about +20 lines of signatures).

### 5.2 After the current wave (H, I, G3) is merged

**W1 — dead code** (§2.6): the 24 dead member functions (two of them declared only), the write-only members `componentError_`, `angleKinFactor_`,
`hasShiftGroupFollowers_`, `hasNormGroupFollowers_`, `energyShiftSqrtNominal_`, `componentCalculationMutex_` with their
writes, the two never-defined declarations, and the unreachable 52 lines of `FillEnergyShiftsFromParamsImpl` (5261–5312;
brief I edits that function now). `EData.cpp` −70, `EPoint.cpp` −24 and the eight statements that write `angleKinFactor_`, `ESegment.cpp` −34, and their
declarations in the headers. Check: a removed member that still has a caller does not compile; 179 × 3.

**W2 — data covariance to `src/EDataCovariance.cpp`**: `EData.cpp` 4289–4547, 259 lines, after brief I has moved the
seeds into `Results` (it edits 4442 now). Checks: 179 × 3 (`covariance.out` of the covariance projects), ctest
`pyazr_objective`, `upgrade_history`, and the poison test, which reads the seeds.

### 5.3 After 3d

**F1 — `src/EDataChecks.cpp`**: the refusals that step 3b/3d does not turn into `Build` checks (today `EData.cpp`
44–230 and 266–325, 247 lines without the includes) and the phase-shift checks of `Initialize` (1764–1848, 85 lines);
`NaturalAllowed` also resolves the sample (item 11 of §7), which must stay in the set-up wherever the check goes, and
the positional `bool isTest` of seven helpers becomes a named argument.
**F2 — one segment builder**: `Fill` and `MakePoints` share their set-up (duplicate 2 of §2.6) once the old-format
overloads are gone. Checks: the refusal tests (`tests/*/models/bad_*.azr` through their `check.sh`), the five refusals
of `project_format3.md` §8.5, 179 × 3.

### 5.4 After S7

**M1 — `src/EDataParameters.cpp`** (layout, names, groups, prior centres, nuisance blocks; ≤ 697 lines), **M2 —
`src/EDataHandOver.cpp`** (the `Fill*FromParams` functions, penalties, the shift hand-over; ≤ 355), **M3 —
`src/ESegmentParameters.cpp`** (norm-group and shift parameters, shift function, shift priors; 217) and
`src/ESegmentShifts.cpp` (re-energising, refill tables; 270). S7 itself
deletes the frame code of `EPoint` (209) and `ESegment` (136). Checks: 179 × 3; the archived fits of `objective.md`
§4; `upgrade_history` (a pooled copy equals a fresh one after shifts move); `fit_threads`.

### 5.5 After S10–S12

The composition leaves the classes: S12 deletes `EData.cpp` 3195–3651 (457) and 1465–1646 (182), and the component
copies of duplicates 3–5 of §2.6; S10/S11 leave `ESegment`'s composition helpers (148) to the tree. No move of this plan
touches them before.

### 5.6 Decision D1: the parts Phases 5 and 6 rewrite

Under the rule of §4.2 these stay where they are until Phase 5 and 6, and the end of Phase 4 then leaves `EData.cpp` at
1,979 lines (722 + target effects 756 + external capture 501, plus what S5 leaves of `Clone` and the parallel pass) and
`EPoint.cpp` at 2,718 (3,006 less dead code 24, mapping 55 and kinematics 209), less what S7 takes out. If the director decides to cut them mechanically before those phases (recommended in
§8), the moves are, in this order and each alone: `EDataTargetEffects.cpp` (756, after S12), `EDataExternalCapture.cpp`
(501), `EPointCalculate.cpp` (394 less S5's deletions, after S5), `EPointFolds.cpp` with `include/TELinearForm.h` (463
and 40), `EPointChannelTables.cpp` with `include/EPointTables.h` (657 and 92), `EPointLegendre.cpp` (133),
`EPointExternalCapture.cpp` (248). Checks for each: 179 × 3, `ctest -LE slow` (target-effect, Doppler, UDR, beam-profile
and THM tests), 1 against 8 threads, and the T2 benchmark (`scripts/bench_upgrade.sh`) against the base: a fold or a
table accessor that stopped being inlined shows there first.

### 5.7 End state

| file | lines (counted at 070756a1) |
|---|---|
| `src/EData.cpp` | 722, plus what S5 leaves of `Clone` and the parallel pass (≤ 193) |
| `src/EDataOutputFiles.cpp`, `EDataCheckFiles.cpp` | 490, 240 |
| `src/EDataThm.cpp`, `EDataThmReports.cpp` | 772, 461 |
| `src/EDataCovariance.cpp`, `EDataChecks.cpp` | 259, ≤ 332 |
| `src/EDataParameters.cpp`, `EDataHandOver.cpp` | ≤ 697, ≤ 355 |
| `src/EDataTargetEffects.cpp`, `EDataExternalCapture.cpp` (D1) | 756, 501 |
| `src/EPoint.cpp` | 823, less the datum S7 takes out |
| `src/EPointCalculate.cpp`, `EPointFolds.cpp`, `EPointChannelTables.cpp`, `EPointLegendre.cpp`, `EPointExternalCapture.cpp` (D1) | 394, 463, 657, 133, 248 |
| `src/ESegment.cpp`, `ESegmentParameters.cpp`, `ESegmentShifts.cpp` | 641, 217, 270 |
| `include/EData.h`, `EPoint.h`, `ESegment.h` | 546, 905, 602 less W1's declarations; `EPoint.h` also less 40 (`TELinearForm.h`), 92 (`EPointTables.h`) and 24 (`ScopedPointMappingOff`, S5); not counted further |
| deleted by the steps (not moved) | `EData.cpp` 1,026 (components 457, component grids 182, mapping 267, result layout 52, dead 68); `EPoint.cpp` 288 (kinematics 209, mapping 55, dead 24); `ESegment.cpp` 318 (frame 136, composition 148, dead 34) |

The rows sum back to the files: `EData.cpp` 722 + 193 + 490 + 240 + 259 + 332 + 697 + 355 + 756 + 501 + 1,026 =
5,571; `EPoint.cpp` 823 + 394 + 463 + 657 + 133 + 248 + 288 = 3,006; `ESegment.cpp` 641 + 217 + 270 + 318 = 1,446;
`EDataThm.cpp` 772 + 461 = 1,233.

## 6. Risks

**Floating-point bits.** The build has no link-time optimization and no `-march` (`CMakeLists.txt`, `src/CMakeLists.txt`):
a function in another translation unit is not inlined into its callers now, and a verbatim move leaves every
expression as it is. What can move bits: (i) changing a function from out-of-line to inline or back — the table getters
`GetLoElement`, `GetSqrtPenetrability`, `GetExpCoulombPhase`, `GetExpHardSpherePhase` (3 lines each in `EPoint.cpp`,
called per channel from `AMatrixFunc.cpp`, `RMatrixFunc.cpp`, `GenMatrixFunc.cpp`, `THMMatrixFunc.cpp`) are loads and do
not change bits, but `GetCoulombAmplitude(s)` and `GetThmFormFactor` do complex arithmetic that an inlined copy may
contract into fused multiply–adds where the compiler contracts across statements (GCC's default
`-ffp-contract=fast` on a target with FMA, e.g. Linux aarch64; x86-64 without `-march` has no FMA, and Apple clang
contracts within a statement only) — so a move keeps the inline status, and a change of it is its own commit, checked
also on macOS CI; (ii) reordering a sum (N4 keeps the per-segment χ² order); (iii) splitting `Calculate` so that its
`static thread_local` matrix functions (`EPoint.cpp` 1964, 1990) become two objects (memory, not bits).

**Hot loops.** The point pass calls the table getters and the result accessors (`Rec`, `GetFitCrossSection`,
`SetFitCrossSection`, inline in `EPoint.h`) per channel and per sub-point. None of N1–N4, W1, W2 touches them. The D1
moves of the tables and folds are measured with the T2 benchmark.

**OpenMP.** Regions in the data files: `EData.cpp` 1137 (grid re-placement, `schedule(dynamic, 4)`, atomics and a mutex
captured by reference), 2000 and 2186 (Legendre, Coulomb amplitude), 2097 (`shared(localStop, configure)`, unnamed
`critical` at 2105), 3115 and 3396 (external capture, named `critical(ec_progress)` also at 2964, 3126, 3404, 3410),
3302 (components, unnamed `critical` at 3310), 4477 and 4541 (covariance; the latter `if (!omp_in_parallel())`), 5247
(re-energising), 5556 (`reduction(+ : failed)`); `EPoint.cpp` 1933 and `ModelPass.cpp` 219 share
`critical(azure_failed_point)`. Named and unnamed critical sections are global by name across translation units, so
moving a region keeps its lock. N1–N4 move no region; W2 moves 4477 and 4541; D1 moves 1137, 3115, 2964.

**State that must move with its users.** `tlPointMappingOff` (`EPoint.cpp` 1917, `thread_local` in an anonymous
namespace) and `ScopedPointMappingOff` (1920–1922) must stay in one file: split, each file would get its own flag and
mapping would silently stay on in the other. Function-level statics move with their function: `EData.cpp` 970
(`ecNoted`), 1169 (timing counter), 4816, 4910 (`warned`); `EPoint.cpp` 1931 (`reported`), 1964, 1990; and the
immutable `static const` test hooks at 1945–1946.

**Private access.** The classes declare no friends; moved functions stay members of the same class and keep their access.
A function turned into a free function needs public `const` accessors that do not exist today (§3).

**Consistency checks.** Check 11 lists `src/EData.cpp` and `src/ESegment.cpp` as files that may call the penalty
functions (N1, M2: allowance rows); check 25 names `include/EPoint.h` and `include/EData.h` (a header split must add its
new headers to the list, or a `Results` member could be declared there unchecked); check 26 records the sizes (lower
with `--tighten-sizes` after each step; a new file above 1,000 lines fails); check 28 allows composition only in
`ModelPass.cpp`.

**API and pyazr.** `api/src/AZUREAPI.cpp` calls 40 member functions of `EData` (94 call sites), 15 of `EPoint` (71) and 22
of `ESegment` (54) by unambiguous name, and `src/AZUREGrad.cpp` 18 (63), 24 (144) and 25 (83); pyazr reaches them only
through the binding (`api/src/pybind_api.cpp` calls none directly). A move inside a class changes none of these; S5 and
S7 change signatures and are where the API is adapted.

**Copies.** `EData::Clone` copies member by member: a step that adds a member (none of N1–N4) must add it to `Clone` and
to the `GridNucleus`-style "not copied" list where it is a cache.

**Merges.** All three running branches edit `src/CMakeLists.txt` line 10 and `repo_consistency_test.py`; the new files go
on a line of their own and the check-11 rows at the end of `ALLOWED_PENALTY_USES`.

## 7. Defects and oddities met while reading (not fixed)

1. `EData.h` 476: the comment of `covarianceError_` ("a component failed CheckComponentPairs") describes
   `componentError_` (473), which nothing reads.
2. `EData.cpp` 955–960: the documentation of `ReplaceSubPointGrids` says a grid within 10⁻¹² (relative) is kept; since
   the fix of cause A (`fit_threads.md` §2) only a bit-identical grid is kept (1049).
3. `ESegment.h` 171: `GetNominalEnergyShiftSqrt` is documented as "the penalty's reference"; the penalty is centred on 0
   (`ESegment.cpp` 855) and the function has no caller.
4. `EPoint.cpp` 2121–2122: edge blending of a target-effect point blends σ only; σ_E1 and σ_E2 keep the integrated values,
   while the automatic-tolerance branch sets all three to the bare values (2090–2091). An E1-only segment with a blended
   edge therefore compares an unblended value. Inferred, not run.
5. `EData.cpp` 1264–1278: the Doppler (particle entrance) and β-delayed refusals of a target effect are applied to base
   segments only; a component with another entrance pair is laid out (1465–1646) without them. Inferred, not run.
6. `EPoint.cpp` 824–830, 857–858, 876–885, 917–918: `angleKinFactor_` and `crossSectionKinFactor_` are assigned and then
   assigned again through their setters with the same value.
7. `EData.cpp` 3232: an invalid exit key of a component is a "Warning" and the run continues.
8. `EData.cpp` 2459–2464, 2509–2541: the writer recomputes six terms of the objective (norm prior, shift, screening,
   density, Doppler, thickness) instead of reading what the objective computed; `Results::objective` holds the other
   terms (Phase 2, open item "the writer from the terms").
9. `EPoint.cpp` 1964, 1990: `static thread_local` matrix functions constructed with the first call's `CNuc` and `Config`
   and re-pointed by `Reset` on every call; they persist across evaluations and sessions on a thread (Phase 6,
   no hidden state).
10. `ESegment.cpp` 916–917, `EData.cpp` 5254, 5277, 5289, 5309: comments that tell the history of a change ("after fixing
    EData::Clone", "CRITICAL FIX") rather than what the code does.
11. `EData.cpp` 218: `NaturalAllowed`, a refusal, also writes the resolved sample into the segment
    (`SetNaturalSample`); a caller that skips the check skips the resolution.
12. `EData.cpp` 93, 116, 131, 155, 166, 191, 200: positional `bool isTest` arguments (data or test segment), against the
    upgrade's rule of no positional booleans; `TotalAllowed`, `BetaDelayedAllowed`, `NaturalAllowed` take a non-`const`
    segment.

## 8. Decisions for the director

**Decided by the director, 2026-11-10.** D1 yes: the parts Phases 5 and 6 rewrite are cut into files after S5 and
S12, one mechanical commit per file with the benchmark. D2: the four allowance rows now; the six totals go into
`Results::objective` as a separate step. D3: flat `src/` with class-prefixed names. D4: members now; free functions
on `const` records when S7 provides the records. D5: the duplicates are fixed on their own branches; this track moves
code and deletes dead code. Spot-checked in review: the line ranges of N1, the loop that never runs (3452), the
merge of this note. N1–N4 are launched as brief P.

- **D1** — Cut the parts that Phases 5 and 6 rewrite (target effects 756, external capture 501 + 248, folds 463, channel
  tables 657, Legendre 133, `Calculate` 394) into files before those phases, in a wave in which no other branch edits
  them? Recommended: yes, after S5 and S12, because otherwise `EData.cpp` (1,979) and `EPoint.cpp` (2,718) stay far
  above the 800-line aim for the whole of Phase 5, and each phase then replaces one file instead of carving code out of a
  large one; the cost is a mechanical commit per file with the checks of §1.
- **D2** — N1's check-11 allowance (four rows) now, or first store the six recomputed totals in `Results::objective`
  (a behaviour-neutral change of the writer, Phase 2's open item)? Recommended: allowance now, totals as a separate
  step.
- **D3** — Placement: flat `src/` with class-prefixed names (chosen, as the `CNuc` split) or a `src/data/` directory?
- **D4** — The outputs as members (chosen) or as free functions over `const EData&` with new `const` accessors?
- **D5** — W1 deletes the duplicate state (`operationType_` / `segmentOperationType_`) and unifies the shift prior
  (`AddShiftPenalty` / `ShiftPenaltyParts`), or are those fixes of §2.6 left to their own `fix/` branches? Recommended:
  own branches; this track only moves and deletes dead code.

## 9. Not determined

- The sizes after S5, S7, S12 and 3d: those steps add code whose size is not known; §5.7 counts only what they remove.
- Which load-time refusals step 3b/3d turns into `Build` checks (G3 is deciding it now); `EDataChecks.cpp` is an upper
  bound.
- Whether `shifts.out` is written by any of the 179 menu-1 runs: not checked; N1 names two projects to compare by hand.
- Items 4 and 5 of §7 were read in the code, not run.
- No build, no run and no test of this note's steps was made: they are a plan.
