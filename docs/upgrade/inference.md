# The `inference/` library: the parameter vector, the priors, the rule of the term list

Status: built, 2026-11-16, on `feature/upgrade-2-inference`, cut from `feature/upgrade` at `f252ba30`. Section 7,
"Done", says what exists and what was run; sections 1 to 6 are the design as written before the code.
`ROADMAP.md` Phase 2 and `ARCHITECTURE.md` section 3 planned a library holding the fit's parameter vector θ and the
objective χ²(θ). This note says which part of the two can be a library today, which part cannot without a redesign,
and how each move is tested. Every output file stays identical: code is moved with its order of operations, not
rewritten.

## 1. What a library may be

`inference/` is the CMake target `azure2_inference`, namespace `azure2::inference`, public headers
`inference/include/azure2/inference/`. It includes its own headers and the C++ standard library, nothing else: no
header of `src/`, `include/` or `api/`, no Minuit2, no GSL. Three things enforce it, as for `physics/`: the target
links nothing; the ctest `inference_standalone` configures and builds the directory alone; a check of
`tests/consistency/repo_consistency_test.py` refuses any other include. The plan let it use `physics/` and
`experiment/`; nothing that moves needs them.

## 2. What the objective depends on today

`src/Objective.cpp` (495 lines) and `include/Objective.h` (264 lines) include 14 engine headers. Each of the
thirteen terms reads engine objects: T1-T4 walk the `Segment`s of a `WorkingCopy` and call the model pass
(`azure2::PointResidual`, `ProfileThmGroup`, `ProfileNormChiSquared`); T2 is `WorkingCopy::CovarianceChi2Correction`;
T6 walks `JGroup` and `ALevel` and reads `LevelState::ParkNorm`; T7 reads `AbundanceSamples`; T8 calls the working
copy's sums; T9 and T11 are `ModelData` and `GammaDecaySet`. The data χ² is the end of the model pass: its inputs
are the caller's `Results` buffer and the working copy, not arrays. To give the terms arrays, every evaluation
would copy the residuals, norms and prior tables out of those objects first; that is a second data path, which the
upgrade removed elsewhere, and a redesign of the evaluation, not a move.

So the cut is: **the library owns what is defined without an engine object; the engine keeps the walk.**

## 3. What moves

Each entry: what it is, where it is today, its header in `inference/include/azure2/inference/`, and how the engine
dependency is cut.

1. **The kinds of parameter, their table, the record of one parameter.** `include/ParameterKind.h` →
   `ParameterKind.h`. Nothing to cut: plain data already.
2. **The level block of the vector.** `include/LevelParameterLayout.h` and its source → `LevelParameterLayout.h`,
   `.cpp`. Nothing to cut: it is built from two arrays of sizes.
3. **The full vector's index map, packed ↔ full.** `ParamDesc` and `ParamIndexMap` of `include/AZUREGrad.h` and
   `src/AZUREGrad.cpp` → `ParameterIndexMap.h`, `.cpp`. The two `friend` builders go: the map is built by public
   appends, one per kind, and `Pack(fixed)`. `BuildParamIndexMap` (the walk of the compound and the data) stays in
   the engine and calls them.
4. **The support of the posterior.** `PosteriorSupport` (`Limit`, `Contains`) → `PosteriorSupport.h`. The two
   makers that read Minuit's parameters and `<parameterSettings>` become engine functions
   (`include/PosteriorSupport.h`).
5. **One Gaussian prior** ((x − x₀)/σ)², its gradient 2(x − x₀)/σ² and curvature 2/σ², written out at each use
   today → `GaussianPrior.h`: functions of three numbers. The engine's enumerations (`ForEachDiagonalPrior`) stay
   and call them.
6. **The prior covariance** (p − μ)ᵀC⁻¹(p − μ): the file's grammar, the marginalisation of fixed parameters,
   C = LLᵀ, the whitened residuals, the gradient, the rows, the logarithm of the normalisation.
   `include/PriorCovariance.h` and its source → `PriorCovariance.h`, `.cpp`. The grammar reads a `std::istream`;
   opening the file through the cache of input files, the path relative to the project and the messages to the
   run's stream stay in the engine (`src/PriorCovarianceFile.cpp`).
7. **The rule of the term list**: the sum in order, the sum with a breakdown, the sampler's sum and support, the
   terms of the model alone. `Objective` in `src/Objective.cpp` → `TermList.h` (`Term<Context>`,
   `TermList<Context>`), a template over the context a term reads; the engine instantiates it with
   `ObjectiveContext`.

`x` and `x₀` carry the unit of the parameter's kind (the column `unit` of `ParamKindTable()`): MeV for a level energy
and an energy shift, MeV^½ for a reduced-width amplitude, eV for U_e, atoms cm⁻² for an areal density, K for a
Doppler temperature. The vector's order is the one of the parameter file: the level block (J^π group after group,
level after level, the energy then one amplitude per channel), the abundances, the normalizations, the three kinds of
energy shift, screening, density, Doppler, thickness, the THM coherent background, the apparatus parameters.

## 4. What stays in the engine, and why

- **The terms themselves** (`src/Objective.cpp`), `ForEach…Prior`, `AppendModelRows`, `ReplaceDiagonalPriors`:
  section 2. A new term is still one class there (`CLAUDE.md`, recipe "A term of the objective").
- **The values, limits and fixed flags of the vector**: `AZUREParams` holds them as Minuit2's `MnUserParameters`,
  which MIGRAD, HESSE and MINOS take as they are. A plain `ParameterVector` class would be a second holder of the
  same numbers, kept equal by hand. The library therefore describes the vector (kinds, records, layout, index map,
  support) and does not hold it.
- **Physical ↔ fit transforms** (Γ in eV ↔ γ in MeV^½, Brune and Park ↔ formal): `physics/`'s `LevelParameters.h`
  with the engine's `LevelTransforms.cpp`; they need the channel functions.
- **`ParameterLimitsManager`** (reads `<parameterSettings>`, converts a width limit with the penetrability), the
  fitters (`src/Fit*.cpp`, on `AZURECalc`), the gradient (`src/AZUREGrad.cpp`), the sampler.

## 5. Slices, each a commit with every output file identical

1. **The parameter vector.** The target, entries 1-4 of section 3, the engine's names of the library's types
   (`include/InferenceTypes.h`, as `ChannelSpace.h` names `physics/`'s). Tests, none needing a project file:
   `inference_parameter_kinds` (one row per kind in the order of the codes; a row moved or a key repeated fails);
   `inference_level_parameter_layout` (three groups of unequal size: every slot's index, name and level number by
   hand; an amplitude placed before its energy fails); `inference_parameter_index_map` (a vector with one parameter
   of every kind and two fixed: each lookup, the offsets, packed ↔ full; a kind appended out of the file's order
   moves an index and fails); `inference_posterior_support` (inside, on and outside a limit, a vector too short, a
   NaN). `upgrade_level_parameter_layout` keeps checking the layout against the engine's projects.
2. **The priors and the list.** Entries 5-7. Tests: `inference_gaussian_prior` (value, gradient and curvature against
   central differences); `inference_prior_covariance` (a 2 × 2 and a 3 × 3 matrix by hand: χ² to 10⁻¹², the gradient
   against central differences, the rows' Jacobian = L⁻¹, a fixed parameter marginalised = the Gaussian of the
   remaining block, ln of the normalisation, each refusal of the grammar by its message); `inference_term_list`
   (terms that add known numbers: the order of the floating-point sum, breakdown = parts, the sampler's default
   takes the whole term, the model-only subset, a term that excludes). `sampler_objective`, `pyazr_objective` and
   `tests/prior_covariance` keep checking the engine's use.
3. **Documents**: `CLAUDE.md`, `FOR_THE_OWNER.md`, `ROADMAP.md`, `ARCHITECTURE.md`, the Developer Guide and its
   class tables, the consistency tables.

Gate of each code slice: build with 0 warnings; every project in three parametrizations against the reference
outputs (887, 851, 887 files identical); the eight fits identical to the base build; `ctest -LE slow`;
`tests/run_tests.sh`. Instructions at one thread of a T2 gradient fit within 1 % of the base build.

## 6. The output writer (audit of `WorkingCopy::WriteChiSquaredFile`, `src/DataOutputFiles.cpp`)

`chiSquared.out` takes its numbers from two places. (a) What the terms of the last evaluation stored in the
caller's `Results::objective`: each segment's χ² and the normalization used, Park's wall, the abundances, the
nuisance priors, the model data, the prior covariance. (b) Sums the writer forms again from the working copy with
the gathers and the arithmetic the terms use (since `objective_arrays.md` slices 1 and 2: `AppendNormPrior`,
`AppendShiftPriors`, `AppendTargetPriors` of `include/DiagonalPriors.h` and `PriorsChi2`): each segment's
normalization prior, the shift total and the target and apparatus totals, equal to the terms' while the working copy
is the one the last evaluation filled.  (Until slice 2 the normalization prior was written out a second time in the
writer, with its own rule.)
Found: T11 (the priors of the `gamma_decay` records) had no field, so for such a project the file's totals did not
add up to the minimised χ² (since 2026-11-16 the term stores its value in `Results::objective.gammaDecay` and the
writer prints `Total-Gamma-Decay-Chi-Squared` whenever the project has a record; the PI's word); a session's `write_outputs` refreshes the model-only terms (`ModelChi2`: T6, T7, T9, T11)
but neither the nuisance total nor the prior covariance's (that one is set by the session's last `CalculateChi2RWA`). A writer that asks the terms needs the compound, the level state and the nuisance function,
which a `WorkingCopy` method does not have, and a new field changes an output file: not done here.

## 7. Done (2026-11-16)

Slices 1 and 2 are built as designed, in two commits; slice 3 is the documents.

**The library.** `inference/`: seven headers (743 lines), four sources (421 lines), seven test programs and their
helper (553 lines). It includes `azure2/inference/` headers and the standard library, and links nothing. The
engine names its types through `include/InferenceTypes.h` (35 lines).

| header | holds | engine side that stays |
|---|---|---|
| `ParameterKind.h` (234) | `ParamKind`, `ParamKindTable()`, `ParameterRecord` | `AZUREParams` (the values) |
| `LevelParameterLayout.h` (68) | the level block | `CompoundSystem::ParameterLayout()` |
| `ParameterIndexMap.h` (153) | `ParamDesc`, `ParamIndexMap` | `BuildParamIndexMap` (`src/AZUREGrad.cpp`) |
| `PosteriorSupport.h` (35) | the limits, `Contains` | `PosteriorSupportOfFit`, `PosteriorSupportOfParameters` |
| `GaussianPrior.h` (34) | value, gradient, curvature | `ForEachDiagonalPrior` and its kin (`Objective.h`) |
| `PriorCovariance.h` (113) | grammar, set-up, arithmetic | `src/PriorCovarianceFile.cpp` (path, cache, messages) |
| `TermList.h` (106) | `Term<Context>`, `TermList<Context>` | the thirteen terms (`src/Objective.cpp`) |

In the engine: `src/Objective.cpp` 495 → 468 lines, `include/Objective.h` 264 → 252, `src/AZUREGrad.cpp` 1259 → 1108,
`include/AZUREGrad.h` 476 → 348; `include/ParameterKind.h`, `include/LevelParameterLayout.h`,
`src/LevelParameterLayout.cpp`, `include/PriorCovariance.h` and `src/PriorCovariance.cpp` are gone. The Gaussian
prior's three functions replace 30 written-out copies (the second commit's message says 25). Removed with no caller
left: `ParamIndexMap::EIndex`, `GammaIndex` and the four block offsets.

**Differences from the design.** `GaussianPriorChi2` forms the value as r·r, r = (x − x₀)/σ. Of the copies it
replaces, 11 were written `pow(r, 2.0)` and 7 `r * r`. GCC at `-O3` compiles both to one multiplication (read in the
assembly of a two-line file), so no number moved in the gate; a build without optimisation called the library's
`pow` at those 11 places and no longer does, which can move a last digit there.

**What was run.** `work/director/verify_branch.sh wt-infer infer tab-t74-1 trial92 6`, once per code commit.

| commit | warnings | files identical | differ, new | fits identical | quick tests | regression scripts |
|---|---|---|---|---|---|---|
| slice 1 | 0 | 887, 851, 887 | 0, 0 | 8 of 8 | 319 of 319 | 172 of 172 |
| slice 2 | 0 | 887, 851, 887 | 0, 0 | 8 of 8 | 322 of 322 | 172 of 172 |

The three numbers of files are the Brune, the formal (`--no-brune`) and the Park parametrization, against the
reference outputs `tab-t74-1`; the fits are compared with the base build `work/b-trial92`.

The quick set had 314 tests at the start; the library adds eight (seven programs and the stand-alone build).
User-space instructions at one thread, base build `work/b-trial92` against slice 2
(`work/director/infer_speed.sh`): a calculation of T2 1,242,597,477 → 1,242,582,263 (−0.001 %); a gradient fit of
T2 169,174,329,302 → 169,559,125,470 (+0.23 %); a gradient fit of 13N 185,403,197,700 → 185,402,404,063
(−0.0004 %); no output file differs. The 0.23 % was not traced to a function.

**Not run.** The slow tests (`ctest -L slow`), the sanitizers, the thread comparison, any other compiler or
platform (the CI was not run: nothing was pushed), a build without optimisation. No test was shown to fail by an
injected error.

**Not done, and why.**

- The terms, `ForEachDiagonalPrior`, `AppendModelRows` and `ReplaceDiagonalPriors` are engine code (section 2).
  Roadmap box "`inference/` target" stays open for them.
- No class `ParameterVector` (section 4): the values are `AZUREParams` on Minuit2's parameters.
- The fitters and the gradient are engine code.
- The writer is audited (section 6) and unchanged apart from calling `GaussianPriorChi2`. T11's field in
  `chiSquared.out` (`Total-Gamma-Decay-Chi-Squared`) was added afterwards, with the PI's word.
- `ObjectiveContext` still holds the two flags of a request (`thmOnly`, `segmentPriors`) beside what a term reads.

## 8. The diagonal priors on one list (2026-11-16, `objective_arrays.md` slices 1-3)

Slices 1, 2 and 3 of `objective_arrays.md` section 5 are built, one commit each.

- `inference/`'s `DiagonalPriors.h`: the list of diagonal priors (full index, the value the objective reads, centre,
  width, the source of the prior, a group: the segment or the apparatus row) and its arithmetic, `AddPriorsChi2`
  (into a running sum), `PriorsChi2` (a subtotal from 0, by source and group), `AddPriorsGradient` (on `p[index]`),
  `AddPriorsCurvature`; test `inference_diagonal_priors`.
- The engine's gathers, `src/DiagonalPriors.cpp` (`include/DiagonalPriors.h`): `AppendNormPrior`,
  `AppendShiftPriors`, `AppendTargetPriors`, `AppendAbundancePriors`, `DiagonalPriorsOf` (every prior with its
  index, in the order the API's rows expect) and, for T5, `NuisancePriorOf` / `NuisancePriorsOf`.  Each rule of
  eligibility is written there once; the terms (T3, T4, T5, T8), the gradients, G2, the least-squares rows, the
  sampler (`ForEachDiagonalPrior`, a loop over `DiagonalPriorsOf`) and the `chiSquared.out` writer read it.
- Gone: the five `WorkingCopy::*Penalty` sums, the per-object penalties of `TargetEffect`, `Segment` and
  `ApparatusParam`, `Segment::AddShiftPenalty` and `ShiftPenaltyParts`, the six `ForEach...Prior` helpers of
  `Objective.h`, the two written-out copies of the nuisance rule (`AZURECalcNuisance.cpp`, `AZUREAPIChi2.cpp`).
- Before slice 2, a scratch build compared the old enumeration's priors with the gather's on 343 runs (every project
  of `tests/` with data, and the runs of 26 regression scripts): equal in all, element by element (the commit message
  of slice 2 has the counts).

Not moved then: T7's value and gradient, T6, T1a and T1b (slices 4-6 of the note; section 9).

## 9. The abundances, Park's wall and the profiled norms on arrays (2026-11-17, `objective_arrays.md` slices 4-6)

Slices 4, 5 and 6 of `objective_arrays.md` section 5 are built, one commit each.  Each term keeps its class in
`src/Objective.cpp`; its walk of the engine's objects is one gather in the engine, its arithmetic one header of
`inference/`, called by every path.

| term | `inference/` | engine gather | callers |
|---|---|---|---|
| T7, abundances | `AbundancePriors.h`: `AbundanceBlock` (w, prior centre, width, active per ordinal; sum centre and width per sample), `AbundanceChi2` (with or without the priors: the value, the sampler's part), `AbundanceGradient`, `AbundanceSumRows` | `AbundanceBlockOf` (`src/DiagonalPriors.cpp`) | `AbundanceTerm`, `AddAbundancePenaltyGradient`, `AppendModelRows`, `AppendAbundancePriors` |
| T6, Park's wall | `ParkWall.h`: `ParkLevel` (J; the particle channels' γ, S', S''), `OnParkWall` (J < 0; false for NaN), `ParkWallChi2`, `ParkWallGradient` (the gradient's form), `ParkWallRow` (the row's form) | `ParkLevelsOf` (`src/Objective.cpp`) | `ParkTerm` (value, sampler's support), `AddParkPenaltyGradient`, `AppendModelRows` |
| T1a, T1b | `LinearProfile.h`, `src/LinearProfile.cpp`: `ProfiledNormChi2`, `ThmProfile`, `SolveThmProfile`, `ThmProfileDerivative` (moved unchanged from `src/ThmExperiment.cpp` and `include/ThmProfile.h`, which goes) | `Segment::ProfileNormChiSquared`, `GatherThmGroup` (`src/ThmNorm.cpp`) | as before |

- The gradient and the rows of Park's wall round differently (factor 2J/s² times dJ against J/s and dJ/s); both forms
  stay, side by side in `ParkWall.h`, each the expression of its path.  The rule "on the wall" is one function for
  the value, the gradient and the rows.
- `AbundanceSamples::Sum` went (its loop is `AbundanceSum`); `src/AZUREGrad.cpp` 1108 → 1083 lines.
- New tests: `inference_abundance_priors`, `inference_park_wall`, `inference_linear_profile` (the first unit test of
  `SolveThmProfile`); each shown to fail once on an injected arithmetic error.
- Before the walks were switched, a scratch build (not committed) computed the old walks verbatim beside the arrays
  at every call and compared bit for bit, over every project of `tests/` in the three bases, 24 regression scripts,
  25 ctests and ten fits (468 processes): T7 5288 evaluations (value, sampler's part, gradient, rows, prior list),
  T6 value 2767, gradient 62, rows 111 (233 levels on the wall), `ProfiledNormChi2` 8629, `SolveThmProfile` 4066,
  `ThmProfileDerivative` 89; no difference.

**Where the note was wrong in detail.** `ParkLevelsOf`'s flag says whether the channels are gathered (the value
reads J only), and they are gathered for the levels on the wall only (`OnParkWall`), so no walk reads γ, S' or S''
of a level off the wall, as before.  The note's `AbundanceSumRows` signature
(a vector of rows and one of sparse derivatives) became a list of `AbundanceSumRow` (sample, r, 1/s); the engine maps
the sample's ordinals to the vector.  `SolveThmProfile`'s names were kept (`ThmProfile` in a header called
`LinearProfile.h`), as the note says "move unchanged".

**What was run** is in the three commit messages and the hand-over entry; the full check of the three slices is
`work/director/verify_obj2b.out`.

**After (B4, 2026-11-17).** The gradients and rows of T6, T7, T9 and T11 and every total and field of
`chiSquared.out` are now members of the terms' classes (`AddGradient`, `AppendRows`, `WriteTotal`), walked by
`Objective::AddGradient`, `AppendModelRows`, `WriteTotals`; `AddParkPenaltyGradient`,
`AddAbundancePenaltyGradient`, `Objective::ModelDataGradient` and the free `AppendModelRows` went
(`objective_arrays.md` section 8).
