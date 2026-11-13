# Parameter kinds: what a fit parameter is, without reading its name

Format 3 (`project_format3.md`, decision F3, §7.4 check 23) lets the user name data sets, targets and
normalizations, so a parameter's name stops being a fixed pattern. The engine used to decide what a parameter is
from its name (`segment_…_norm`, `density_…`, `…_thickness`); with user names that breaks silently, and every new
nuisance had to be added to each prefix list. Now each parameter carries a **record** — its kind and what it belongs
to — written where it is created and held beside the Minuit parameters in `AZUREParams`
(`include/ParameterKind.h`). Names are labels for people and files. Check 23 of
`tests/consistency/repo_consistency_test.py` fails on a prefix or suffix test of a parameter name.

## 1. The kinds

One enumeration, `ParamKind`, and one table, `ParamKindTable()` (`include/ParameterKind.h`). The code is the `type`
of the API's parameter listing (`AZUREAPI::GetParameterInfo`) and never changes; the key is pyazr's
`Parameter.kind` (`pyazr.parameters.parameter_kinds()` returns the table through the binding). The vector holds the
kinds in the order of their codes, except that the abundances (6) follow the level block (0, 1).

| code | kind (key) | unit | owner fields of the record | created in | name (label) |
|---|---|---|---|---|---|
| 0 | `LevelEnergy` (`energy`) | MeV | jGroup, level, levelNumber | `CNuc::FillMnParams` | `energy_<levelNumber>` |
| 1 | `Gamma` (`width`) | MeV^1/2 (reduced-width amplitude, run's basis) | jGroup, level, channel, levelNumber | `CNuc::FillMnParams` | `width_<levelNumber>_<channel>` |
| 2 | `Norm` (`norm`) | 1 | segment, segmentKey | `EData::FillMnParams` | `EData::SegmentParameterName` |
| 3 | `EnergyShift` (`shift`) | MeV (lab) | segment, segmentKey | `EData::FillMnParams` | ″ |
| 4 | `EnergyShiftSqrt` (`shift_sqrt`) | MeV^1/2 | segment, segmentKey | `EData::FillMnParams` | ″ |
| 5 | `EnergyShiftFn` (`shift_fn`) | the shift function's | segment, segmentKey, fnTerm | `EData::FillMnParams` | ″ |
| 6 | `Abundance` (`abundance`) | 1 | sample, sampleName, nucleus, nucleusZ, nucleusA | `CNuc::FillMnParams` | `CNuc::AbundanceParamName` |
| 7 | `Screening` (`screening`) | eV | effect | `EData::FillMnParams` | `screening_<effect>` |
| 8 | `Density` (`density`) | atoms/cm² | effect | `EData::FillMnParams` | `density_<effect>` |
| 9 | `Doppler` (`doppler`) | K | effect | `EData::FillMnParams` | `doppler_<effect>` |
| 10 | `Thickness` (`thickness`) | atoms/barn | segment, segmentKey | `EData::FillMnParams` | `EData::SegmentParameterName` |
| 11 | `ThmCoherent` (`cbkg`) | HOES amplitude (c1 per MeV) | thmExperiment, jGroup, channel, exitChannel, thmTerm | owners `EData::BuildThmGroups` (`EDataThm.cpp`), added in `EData::FillMnParams` | `ThmCoherentParamStem` + `re0`… |
| 12 | `Apparatus` (`apparatus`) | the row's (`ParameterRecord::apparatusUnit`) | effect (the first target effect), `apparatus` (the row of `ApparatusRows()`), `owner` | grouped by sharing key in `WorkingCopy::SetUpApparatusParams` (`src/ApparatusParameters.cpp`), added last by `AddApparatusParams`; rows: the stopping scale (`src/StoppingScale.cpp`), a table kernel's own parameter (any `KernelParameter` with `fittable`) | the owner's values joined by `.`: `stopping_scale.<material>.<projectile>`, `<record>.<field>.<parameter>` (`apparatus_parameters.md`) |
| 13 | `RadiativeDecay` (`radiative`) | eV (`total_width`), 1 (`branching`, `mixing`) | jGroup, level, channel (the photon channel whose slot it holds), levelNumber, `decayRecord`, `decayPart`, `finalState` | `GammaDecaySet::AddParameter` (`src/GammaDecayMap.cpp`), called by `CNuc::FillMnParams` at the slots of the level's photon amplitudes; `CNuc::FillCompoundFromParams` maps them to amplitudes (`gamma_parametrisation.md`) | `<record>.total_width`, `<record>.branching.<final state>`, `<record>.mixing.<final state>` |

Every record also has `ordinal`, its place among the parameters of its kind (the `k` of `GetScreeningEffect(k)`,
`GetThicknessSegment(k)`, …), set by `AZUREParams::Add`. The table's three flags answer the questions the name
tests used to answer:

- `ownPrior` — the parameter has a prior of its own (segment error, `<abundances>`, a target-effect or transmission
  block); a `<parameterSettings>` nuisance flag adds nothing (kinds 2–10).
- `segmentField` — it belongs to a data segment; a run without data does not take it from a parameter file, and the
  GUI's `<parameterSettings>` numbering counts it also when fixed (2, 3, 4, 5, 10).
- `varyFromProject` — its project line's vary flag decides whether it varies; a parameter file's zero or nonzero
  value does not fix or free it (3, 4, 5, 7, 8, 9, 10).

## 2. The sites that recognised a kind by name (before this change), and what they decide now

| site (line numbers before) | what it decided | now |
|---|---|---|
| `src/AZUREParams.cpp` 9–15, 122–126 | a param-file value does not free / fix a shift, screening, density, Doppler or thickness | `varyFromProject` |
| `src/AZUREParams.cpp` 98 | a run without data skips the file's `segment_*` rows | `segmentField` of the matched parameter |
| `include/ParameterLabel.h` `IsNormName`, `IsEnergyShiftName`, `IsEnergyShiftFnName`, `EndsWith` | norm / shift by form | deleted |
| `api/src/AZUREAPI.cpp` `GetNormalizationIndices`, `GetEnergyShiftIndices` | the free norms, the free shift terms | record kind |
| `src/ParameterLimitsManager.cpp` 209–213 | `energy_N` / `width_N_C` → the GUI's `Level N …` setting name | `levelNumber`, `channel` |
| ″ 245–254 | which unmatched settings may bind by their stored index | one-token (parameter-name) entries never; `Level N` entries to that level's kind |
| ″ 262, 321 | the GUI numbering counts fixed `segment` parameters | `segmentField` |
| ″ 304–306 | no nuisance prior on a parameter with its own prior | `ownPrior` |
| ″ 399–461, 657–722 | a width's settings are converted eV → reduced amplitude | kind `Gamma` (`priorIsWidth_`) |
| ″ `ApplyParameterSetting`, `IsNuisanceParameter`, `GetParameterError`, `GetConverted…` | (no callers) | deleted |
| `src/AZURECalc.cpp` 1305 | the limits manager's numbering | `segmentField` |
| `src/AZURECalc.cpp` 1451 | Levenberg–Marquardt gradient scale of a density | kind `Density` |
| `src/AZUREMain.cpp` 784 | the band covariance's column names: free level parameters | kinds `LevelEnergy`, `Gamma` |
| `src/AZURECalcMCMC.cpp` 572; its `PARAM_*` enum | a prior covariance on level parameters refused for a physical chain; prior rules | kinds through `FullIndex`; `ParamKind`, `ownPrior` |
| `src/AZURE2.cpp` 1102–1228 `createPrettyParameterName` | the MCMC table's names; matching `<mcmc>` prior rows | `PrettyParameterName(record)`; `FindMcmcPriorParameter` (a row naming no parameter is an error) |
| `src/ParameterLimitsManager.cpp` nuisance conversion | converted every non-segment nuisance row through the level parameters (out of bounds for screening, cbkg, …) | widths only (kind `Gamma`) |
| `src/EData.cpp` 641–650 | `segment_N_norm` and `segment_N_energy_shift prior_centre` rows | whole name against `SegmentParameterName` of every line key: active applied, inactive warned, any other name an error |
| `src/ThmExperiment.cpp` `ApplyThmCoherentValues` | (no callers) | deleted |
| `pyazr/parameters.py` 179; `_KINDS` | an abundance's sample name; the kind keys | record `sample_name`; `parameter_kinds()` |
| `pyazr/azure2.py` 862, 1161–1180 | extrapolation skips segment rows of a `.sav`; `save_fit` of `cbkg_*` | `segment_field`; `thm_*` fields |
| `include/AZUREGrad.h` `ParamKind` | the gradient's kinds | moved to `ParameterKind.h` |

Where a **file** names a parameter (`param.sav` / `param.par`, a prior covariance, `<parameterSettings>`, the
`<mcmc>` prior rows) the name is looked up whole in the vector (`AZUREParams::IndexOfName`, `PriorCovariance::Load`)
or compared whole with the name its creation site gives it; the record then says what it is. One file grammar is
read without a parameter vector and stays as it is, listed in check 23: the example
`pyazr/examples/transform_widths.py`.  (pyazr's reading of the old `prior_centre` rows went with the old file class,
pyazr 3.0.)

Still a second walk of the layout (positions, not names; each must follow the creation order):
`BuildParamIndexMap` (`src/AZUREGrad.cpp`), `AZUREAPI::GetParameterInfo` (which now throws if its kind disagrees with
the record), `AZURELabel::Parameter` (`src/ParameterLabel.cpp`), `CNuc::GetTransformParams` /
`ConvertAmplitudeBasis`. Building them from the records is the next step (the inference/ `ParameterVector` of
`ARCHITECTURE.md` §3).

## 3. Recipe: a new kind of fit parameter

1. `include/ParameterKind.h`: the enum value (the next code) and its `ParamKindTable()` row — key, unit, owner
   fields, `ownPrior` / `segmentField` / `varyFromProject`, one line of doc; a `ParameterRecord` factory if its owner
   fields are new [check 23: one row per value, in code order].
2. Create it with `AZUREParams::Add(name, value, step, ParameterRecord::…)` in `CNuc::FillMnParams` or
   `EData::FillMnParams`, its name built there by one function. Nothing else may test the name.
3. Follow the layout in `BuildParamIndexMap` and in `AZUREAPI::GetParameterInfo` (and its fields in pyazr's
   `Parameter.from_record`); the owner fields in the binding's `parameter_records()`.
4. Docs: the table of `docs/source/reference/output_files.rst` and §1 here [check 23]; in
   `tests/pyazr/parameter_kinds_test.py` a project with the kind, whose names are rebuilt from the records.

## 4. The prior of a nuisance parameter: one rule for every kind

Format 3 writes a nuisance as `<value> [<unit>] fit|fixed [prior <centre> <width>[%]]` (`Values.h`, `Nuisance`). The
rule, the same for every kind that can carry a prior (norm, the three shift terms, abundance, `isotope_sum`,
screening, areal density, Doppler temperature, sample thickness, stopping scale, a kernel's parameter):

- the penalty is ((p − centre) / width)², with the **centre as written** and the width as written (a percentage is
  per cent of the centre); it is paid only while the parameter is free (`objective.md` D5);
- the **start** of the fit is the value written; the start does not move the prior;
- a prior without a centre (only an old line read by the old-line reader) is centred as before: on the start value,
  on 0 for the three shift terms, on 1 for the sum of a sample's abundances;
- the sub-point grids and the limits of a fitted density, temperature, stopping scale or kernel parameter are laid out
  from the **centre** where one is written and the start lies inside the range so laid out; otherwise from the start
  (below, "The layout of a fitted value").

Why it matters: a fit continued from its fitted values (the program's `param.sav` reloaded, or the project `pyazr`'s
`save_fit` writes: the fitted value with the original prior) must pay the same prior as the fit it continues.

**The layout of a fitted value.** Four kinds lay out a grid from their value: the areal density (the sub-point grids
span the layout density, the depth integral is cut at the present one), the Doppler temperature (the grids and the
integration windows span the layout temperature, the Gauss-Legendre pieces are u(b)/4 wide), the stopping scale (the
grids span the target at the layout scale; a thick target's Bohr straggling support takes the stopping at the base
scale) and a table kernel's parameter (its grid spans the layout value). Their limits are [f·b, H(b)], b the base:

| kind | lowest f·b | highest H(b) without `grid_max` (σ the prior's width) | no prior |
|---|---|---|---|
| areal density | b/100 | max(b + 3σ, 1.1 b) | 1.5 b |
| Doppler temperature | b/4 (and the pieces) | max(b + 3σ, 1.21 b) | 2 b |
| stopping scale, kernel parameter | b/10 | max(b + 3σ, 1.1 b) | 1.5 b |

`grid_max` (the Doppler fit form's Tmax) replaces H(b). Until this rule, b was the start value: a project written by
`save_fit` (the fitted value as the start, the prior unchanged) was computed on other grids and with other limits than
the fit that produced it. Measured on the base build (the fit, its `param.sav` read by pyazr, written back with
`save_fit`, both sessions evaluated at the same vector):

| kind (project) | objective, fit vs written (rel.) | per data set (rel.) | layout and limits, fit → written |
|---|---|---|---|
| areal density (`tests/target_density/fit_prior.azr`) | 2.1e-14 | meyer_84 1.2e-13 | 1.1e18 → 1.19e18; lower 1e16 → 1.08e16 |
| Doppler temperature (`doppler_fgm` `fit_prior`) | 5.5e-6 | t 6.0e-6, g 5.6e-6 | 363 → 379.6 K; lower 75 → 78.4 K; pieces for 300 → 313.7 K |
| stopping scale (`stopping_yield` `fit.azr`, densities with 10 % priors) | 3.4e-8 | scan1 3e-16, scan2 5.8e-7 | scale 1.12 → 1.126, densities 6.5e17 → 6.92e17, 1.3e18 → 1.23e18 |
| kernel parameter (`rectangular_resolution` `fit.azr`, `prior 14 3`) | 0 | 0, 0 | 23 → 29.0 keV; lower 1.4 → 2.0 keV |
| abundance (`natural_target` `M_data.azr`) | 0 | 0 | no grid |
| sample thickness (`transmission_yield` `fitP`) | 0 | 0 | no grid |

The Doppler difference shows in `chiSquared.out` (4.85481 → 4.85478) and is above `save_fit`'s own check (1e-6), which
then refuses its snapshot; the density's and the kernel's χ² reproduce, but their limits move, so a continued fit may
go where the fit could not.

The rule (`azure2/experiment/LayoutBase.h`, `ChooseLayoutBase`, the one function of the four kinds): the base is
the **centre** when the parameter is fitted, a centre is written, and the start lies in [f·c, H(c)] (`grid_max`
included as written); otherwise the **start**, and when a written centre was passed over (a start outside its range,
on purpose) one line of output says so. What it guarantees:

- a fit started at value = centre and the fit continued from its result share the base (the centre), so the same
  grids, pieces, windows and limits: the written project computes the fit's χ² and its second fit starts at the
  minimum, as long as the fitted value lies inside the limits, which it does (they are the fit's limits);
- value = centre, a prior without a written centre, no prior, or a fixed parameter: the base is the start itself (the
  same double), so every formula computes what it computed before; this is why the whole tree stays identical;
- the program's `param.sav` read back leaves the project's values (= centres) as the start: unchanged.

What it does not cover: a parameter **without a prior** has no centre to lay out from, so a fit continued from a
written project lays its grid out from the fitted value (1.5 × the value, not 1.5 × the start). Writing `grid_max`
(or a prior) into such a project fixes its layout; whether `save_fit` should write one is the PI's decision.

**A new nuisance kind whose value lays out a grid** needs: its `LayoutSpan` (a row of three numbers in
`LayoutBase.h`, or one of the three there), one `ChooseLayoutBase` call where the engine object reads its spec (the
centre from its `PriorCentre`), and the base, not the start, in every place that lays out from it (its grid, its
limits, its hand-over clamp). Before this rule the layout was written per kind (four formulas, the kernel's twice).

**One carrier.** `PriorCentre` (`SetupSpecs.h`: written, value, source; `Or(default)`) is the centre of every kind in
the specs, the one the norms and energy shifts had (`DatasetSpec::normCentre`, `shiftCentre`) made general. The build
fills it from the nuisance (`Build.cpp`, `BuildEffects.cpp`), the engine object of the kind keeps it beside its width
and uses it in its penalty, and the one enumeration of diagonal priors (`ForEachDiagonalPrior`, `Objective.h`)
hands it to the gradient, the least-squares rows and the MCMC priors. The API listing reports it
(`Parameter.prior_nominal`, now also for an abundance and an apparatus parameter).

**Writing back.** The program writes fitted values only to `param.sav` (no prior); reading it back keeps the
project's centres. pyazr's `save_fit` and `Project` write a fitted value into its record through the build's
`ParameterSource` (`experiment::ParameterValueEdit`), which changes the value token only, so the centre and width stay
as written. Every kind with a prior has a source now (the areal density, the Doppler temperature, the screening
potential and the sample thickness had none, so their fitted values were not written back); a source names the
engine's unit where it is not the base unit (`engineUnit`: eV for screening, atoms/b for a thickness).

**A new nuisance kind with a prior** needs: its width and `PriorCentre` in its spec, filled by `PriorCentreOf` where
the build reads the nuisance; the centre kept beside the width by its engine object and used in its penalty; one
line in `ForEachDiagonalPrior` (`Objective.h`, from which the gradient, the least-squares rows and the MCMC take it);
a `ParameterSource` for its write-back; its field in the API listing.

Before (measured by calculations, value ≠ centre against value = centre, at 441a0907): the norm and the constant shift
used the centre; the sqrt(E) shift, the shift function's parameters, the abundance, `isotope_sum`, screening, Doppler
temperature and sample thickness dropped it without a word (the penalty measured from the start, or from 0 / 1); the
areal density, the stopping scale and a kernel parameter refused a centre away from the value.
