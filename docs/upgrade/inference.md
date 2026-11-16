# The `inference/` library: the parameter vector, the priors, the rule of the term list

Status: design, written 2026-11-16 on `feature/upgrade-2-inference`, cut from `feature/upgrade` at `f252ba30`.
Section 6, "Done", says what of it exists; until a line stands there, the rest is a plan. `ROADMAP.md` Phase 2 and
`ARCHITECTURE.md` section 3 planned a library holding the fit's parameter vector θ and the objective χ²(θ). This note
says which part of the two can be a library today, which part cannot without a redesign, and how each move is tested.
Every output file stays identical: code is moved with its order of operations, not rewritten.

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

`chiSquared.out` takes its numbers from three places. (a) What the terms of the last evaluation stored in the
caller's `Results::objective`: each segment's χ² and the normalization used, Park's wall, the abundances, the
nuisance priors, the model data, the prior covariance. (b) Sums the writer forms again from the working copy with
the functions the terms call (`Segment::AddShiftPenalty`, `WorkingCopy::ScreeningPenalty`, `DensityPenalty`,
`DopplerPenalty`, `ThicknessPenalty`, `ApparatusPenalty`): equal to the terms' while the working copy is the one
the last evaluation filled. (c) The normalization prior of each segment, written out a second time in the writer.
Found: T11 (the priors of the `gamma_decay` records) has no field, so for such a project the file's totals do not
add up to the minimised χ²; a session's `write_outputs` refreshes the model-only terms and the prior covariance but
not the nuisance total. A writer that asks the terms needs the compound, the level state and the nuisance function,
which a `WorkingCopy` method does not have, and a new field changes an output file: not done here. Slice 2 removes
(c): the writer and the term call the one `GaussianPrior`.

## 7. Done

Nothing yet.
