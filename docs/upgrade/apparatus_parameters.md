# A fitted parameter of an apparatus effect: the path, and the one kind that carries it

What a fit parameter of an apparatus effect is (a number of a target, beam or detector record that a fit may vary:
the stopping scale s of a material and projectile, the full width w of a rectangular energy resolution), every place
it passes from the record to the outputs, and the generic parameter kind that makes the next one cost no engine code.
Context: `energy_grid.md` §3.3 (kernels and their parameters), §6.2 (the stopping scale); `parameter_kinds.md` (the
kinds of the parameter vector).

## 1. The path of a fitted apparatus parameter before the generic kind (2026-11-11, head 50e9359a)

Four parameters of apparatus records were fitted, each a parameter kind of its own: the areal density N of a target
(kind 8), the free-gas Doppler temperature T (kind 9), the sample thickness n of a transmission or yield (kind 10), the
stopping scale s (kind 12, a day old).  Each passes the same stations; each station named every kind once more.

| # | Station | What it does | File (density, Doppler, thickness, stopping scale alike unless noted) |
|---|---|---|---|
| 1 | Record → spec | the format-3 field (`areal_density <v> <unit> fit prior ...`, `doppler free_gas <T> K fit ...`, `stopping evaluation scale <s> fit ...`) read into the effect's spec | `experiment/src/BuildEffects.cpp` (`TargetInto`, `StoppingInto`, `ScaleOf`, `SharedScales`), `SetupSpecs.h` |
| 2 | Spec → effect | the nominal value, vary flag, prior width, layout value (the largest value the sub-point grids span) kept on the effect | `include/TargetEffect.h`, `src/TargetEffect.cpp` (thickness: `include/Segment.h`) |
| 3 | Grouping | effects (or segments) that share one parameter: density groups, the stopping scale's (material, projectile) key; the list of groups | `src/HandOver.cpp`/`src/TargetSetup.cpp` (density), `src/DataSet.cpp` (Doppler, thickness), `src/StoppingScale.cpp`; the lists in `include/DataSet.h` |
| 4 | Working-copy value | the value an evaluation uses, per effect | `include/SegmentState.h` (`TargetEffectValues`), `src/WorkingCopy.cpp` (`AddTargetEffect`), `src/TargetSetup.cpp` (start values) |
| 5 | Vector slot | name, start value, step, limits, fixed flag, record (`ParameterRecord`), block offset | `src/DataParameters.cpp` (`FillMnParams`, `Layout`), `include/ParameterKind.h` (enum, row, record fields, factory) |
| 6 | Hand-over | vector → working copy, and what must follow (the stopping scale: every sub-point's stopping, every point's thickness) | `src/HandOver.cpp` (`Fill*FromParams`), `src/StoppingScale.cpp`, `include/WorkingCopy.h` |
| 7 | Prior | the penalty ((x − x0)/σ)², its term of the objective, the diagonal-prior list (least squares, covariance), the prior covariance's replacement of it | `src/Objective.cpp` (`TargetPriorTerm`, `ReplacedPriors`), `include/Objective.h` (`ForEach*Prior`), `include/WorkingCopy.h` |
| 8 | Derivative | analytic (Doppler: `FoldTemperatureDerivative`; thickness: the transform's backward), else central differences of the whole χ² | `include/AZUREGrad.h` (`ParamDesc`, `ParamIndexMap`), `src/AZUREGrad.cpp` (`BuildParamIndexMap`, `GradAccum`), `src/GradEnergyTasks.cpp` (cotangents `'d'`, `'t'`), `src/AZURECalc.cpp` (the kinds differenced: `FillShiftColumns`, `Gradient`), `api/src/AZUREAPI.cpp` (the same list with THM) |
| 9 | Outputs | `chiSquared.out` line (`Total-…-Chi-Squared`), labels of the covariance and correlation files, the MCMC table's name | `src/DataOutputFiles.cpp`, `src/ParameterLabel.cpp` (count and text), `src/AZURE2.cpp` (`PrettyParameterName`), `src/CovarianceBand.cpp` (`IdentityKey`) |
| 10 | API and pyazr | the 14-field listing (`GetParameterInfo`), the records (`parameter_records`), the kinds table | `api/src/AZUREAPI.cpp`, `api/src/pybind_api.cpp`, `pyazr/parameters.py` (type codes) |
| 11 | Documents and tests | the kinds table, the record grammar | `docs/source/reference/output_files.rst`, `docs/upgrade/parameter_kinds.md`, `tests/pyazr/parameter_kinds_test.py` |

The stopping scale's commit (fdd11148) edited 20 existing engine-side files for stations 2–10 (`ParameterKind.h`,
`DataSet.h`, `TargetEffect.h/.cpp`, `SegmentState.h`, `WorkingCopy.h/.cpp`, `HandOver.cpp`, `DataParameters.cpp`,
`Objective.h/.cpp`, `AZUREGrad.h/.cpp`, `AZURECalc.cpp`, `CovarianceBand.cpp`, `DataOutputFiles.cpp`,
`ParameterLabel.cpp`, `AZURE2.cpp`, `AZUREAPI.cpp`, `pybind_api.cpp`) plus `src/StoppingScale.cpp` new.  A kernel's
own parameter (`KernelParameter::fittable`) reached none of these stations: the rectangular width was a constant of
the effect.

## 2. The generic kind

Kind 12 is the **apparatus parameter** (`ParamKind::Apparatus`, API key `apparatus`): a number of an apparatus record
that a fit may vary.  It replaces the stopping-scale kind (one day old and unreleased: no project or parameter file
outside this tree names kind 12, so the code keeps its value and its place, last in the vector).

- **Declaration (the record says it is fitted).**  One grammar: the nuisance of format 3,
  `<value> [<unit>] fixed|fit [prior <centre> <width>[%]] [grid_max <value>]`, which the areal density, the Doppler
  temperature and the stopping scale already use.  A table kernel's parameter whose `KernelParameter` row says
  `fittable = true` is generated with this grammar by the schema (`Schema.cpp`), and `BuildEffects.cpp` turns the value
  into an `ApparatusParameterSpec` of the effect (`SetupSpecs.h`): row `kernel`, slot (the parameter's place in the
  row), sharing key (record, field, parameter), name `<record>.<field>.<parameter>`, unit, value, vary, prior,
  grid_max.  No line of either names a kernel.
- **Engine rows** (`include/ApparatusParameters.h`, table `ApparatusRows()` in `src/ApparatusParameters.cpp`): one row
  per *way a value reaches the calculation*, not per parameter.  Row `stopping_scale` (`src/StoppingScale.cpp`): the
  value is `TargetEffectValues::stoppingScale`, and the sub-points' stopping and the points' thickness follow it.  Row
  `kernel`: the value is slot k of `TargetEffectValues::kernelValues`, which the fold hands to the kernel as
  `KernelParameters::values`; every table kernel's parameter uses this row.  A row says how its effects declare their
  parameters, how a value is written into the working copy, and what must follow a change.
- **Groups.**  Declarations with one sharing key are one parameter (one detector record named by two data sets: one
  width; one material and projectile: one scale).  `WorkingCopy::SetUpApparatusParams` makes the groups in row order,
  then in effect order; `DataSet::ApparatusParam(k)` holds each group's declaration, its effects and the prior's
  replacement flag.
- **Engine side, one file** (`src/ApparatusParameters.cpp`): the groups, the vector slots (name, start, step, limits:
  [v0/10, the layout value]), the hand-over, the penalty and its `chiSquared.out` line (one per row with a fitted prior:
  `Total-Stopping-Scale-Chi-Squared`, `Total-Kernel-Parameter-Chi-Squared`), the labels, the API's listing, the
  index map.  Every `switch` over the kinds has one case `Apparatus` that asks the group.
- **Layout.**  The grid of a fitted kernel parameter is laid out at its layout value (grid_max, else
  max(v0 + 3σ, 1.1 v0) with a prior, else 1.5 v0, as the stopping scale's): `KernelParameters::layoutValues`, read by
  `EffectSupport` for the support only.  The fold reads the working copy's value.
- **Derivative.**  Central differences of the whole χ² (the route of the density and the stopping scale): nothing is
  needed from the row or the kernel.

## 3. What the density, the Doppler temperature and the thickness would need to move onto it

Not moved (they have analytic routes and parameter names that every project's `param.par` carries).
- **Areal density** (kind 8): a row whose write sets `TargetEffectValues::density` and whose follow-up is the
  points' thickness and fallback depth (`FillDensitiesFromParams`); its groups are the density groups of
  `TargetSetup.cpp` (a target named by several data sets), the same rule as the sharing key; its names
  (`density_<n>`, `<target>.areal_density`) would become declaration names.  Blocked by the vector order: the
  density's block comes before the Doppler, thickness and THM blocks, the apparatus block last, so every project with
  a fitted density would see its parameters reordered unless the generic block took the density's place.
- **Doppler temperature** (kind 9): a row writing `TargetEffectValues::dopplerTemperature`; the analytic derivative
  `FoldTemperatureDerivative` is the broadening kernel's own and would become the row's derivative through the
  gradient's cotangent list (`GradEnergyTasks.cpp`, `'d'`); same vector-order issue.
- **Sample thickness** (kind 10): a parameter of a segment's transform, not of an effect: a row keyed on the segment
  (owner field `segment`), with the transform's backward as its derivative (`'t'`).  It is the furthest from an
  apparatus record of the three.
