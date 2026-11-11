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
