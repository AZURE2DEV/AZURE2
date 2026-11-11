# pyazr on records: what still uses `pyazr/azrfile.py`

`pyazr/azrfile.py` (`AzrModel`, 3756 lines) is a second grammar of the old project formats in Python. When the engine
reads format 3 only it is deleted; everything pyazr needs to know about a project then comes from the engine: the
format-3 records (`pyazr.Project`, the engine's `Document` and validator), the description of the built project
(`describe`, `AZURE2 --describe --json`), the table (`schema`) and the calculation session (`pyazr.azure2`).
This page lists every remaining use (2026-11-11, pyazr 2.31.0), what it needs and what replaces it. The deletion
follows the conversion of the test tree to format 3.

## Done

- `azure2.datasets`, `active_datasets`, `extrapolations`, `dataset_chi2` on a format-3 session: from the engine's
  description (`pyazr/described.py`, `DescribedDatasets` / `DescribedExtrapolations`: record names, quantity,
  reaction with pair keys, windows, data file, rows read and kept, parameters); `dataset_chi2` keyed by record name.
  Older formats keep `pyazr.datasets.SegmentSet` (it reads lines through the engine's segment-line tokenizer).
- `Project.rename`: the parameter names come from the engine (`_azure2.parameter_names`,
  `experiment::ParameterNamesOf`), matched whole; no pattern on names in Python.
- `pyazr/project.py`, `schema.py`, `grids.py`, `potentials.py`: the engine module without `azrfile._engine_module`.

## Remaining uses in `pyazr/`

| Where | What it uses `AzrModel` for | Replacement |
|---|---|---|
| `azure2.save_fit` | `apply_fit` (fitted level energies and widths into `<levels>`), `_save_cbkg` (THM coherent background), `set_abundance`, `set_segment_values` / `set_prior_centres` (old formats), `set_output_dir` of the verification copy | format 3: `Project` with an engine edit of `<levels>` rows (a binding that writes the fitted values of the level parameters, as `parameter_value_edits` does for nuisances); the THM background as a record field; output directory a `<config>` edit of `Project` |
| `azure2.scattering_length_chi2` | `scattering_lengths` (active `<scatteringLengths>` lines) | the engine's objective term (`penalties()["scattering_length"]` or a binding returning the term's rows) |
| `azure2._parsed("model")` | the two uses above | goes with them |
| `_HybridPotentials._write_potential_records` (`potentials.py`) | the potential records written into the model before `save_fit` writes it | `Project.set` on the `potential` records |
| `bands.trimmed_model`, the worker set-up of `bands.sensitivities` (`nprocs > 1`) | `keep_extrapolations`, `set_extrapolations`, `_tag_index`; a model edited for a worker | `Project.remove` / `add_extrapolation` on a copy |
| `angular.angular_distribution` | `clear_extrapolations`, `add_extrapolation`, `write` on a copy | `Project.new(nuclear_part=...)` + `add_extrapolation` |
| `tabulate` | `AzrModel.from_file`, extrapolation edits, `write` | as `angular` |
| `modelavg.write_averaged_azr`, `azr_parameter_slots` | level setters (`engine_level_keys`) | a `<levels>` edit binding (as for `save_fit`) |
| `transform.levels_from_azr` | the levels and channels of `<levels>` | the engine's level tables (`LevelTablesOfRows`; a binding) or `azure2.level_scheme` |
| `nuclei.merge_projects` | levels, pairs, `<config>`, `_tag_index` | records and the engine's level rows; merging `<levels>` needs a row writer binding |
| `nds` (docstring example) | `add_data_segment` | `Project.add_dataset` with the fetched keywords |
| `datasets._engine_module` | its own loader (works when loaded by path) | stays with `datasets.py` (old formats only) |
| `__init__` | exports `AzrModel`, `AzrLevel`, `AzrChannel`, `ScatteringLengthLine`, `PairMassCheck`, `convert_azr` | `convert_azr` moves to the converter's binding (`_azure2.convert_format3_text`); the rest are deleted with the file |

## Remaining uses in `tests/pyazr/` (58 files)

| What the tests do with it | Files (examples) | Replacement |
|---|---|---|
| open a project and edit data-segment lines (`add_data_segment`, `set_segment_norm`, `clear_data_segments`) | `angular_coefficients`, `data_covariance`, `norm_group`, `shift_fn`, `sqrt_shift`, `thickness`, `objective`, `physical_norm`, `natural_target`, `inactive_lines` | the converted tree's format-3 files and `Project.add_dataset` / `set` |
| edit extrapolations (`add_extrapolation`, `clear_extrapolations`) | `extrap_layout`, `extrapolation_frame`, `geometry_no_transform`, `shared_state`, `roundtrip` | `Project.add_extrapolation` / `remove` |
| edit levels (`add_level`, `set_energy`, `set_fixed`, `set_active`, `engine_level_keys`) | `ec_integrals_formal`, `formal_capture_mask`, `light_excitation`, `zero_width_level`, `model_average`, `azr_layout`, `inactive_lines`, `thm_profiled` | a `<levels>` row edit in `Project` (to be built: the rows are not records) |
| edit target effects (`add_target_effect`, `set_doppler`, `set_screening`, `set_sample_thickness`, `set_abundance`) | `doppler`, `screening`, `thickness`, `natural_target` | `Project.add_target` / `set` |
| THM options and experiments (`set_thm_*`, `clear_thm_*`, `_thm_read`) | the twelve `thm_*` tests | the `<thm>` section is read by the engine (`pybind_thm.cpp`); an edit binding of its options |
| `<config>` edits (`set_output_dir`, `set_channel_radius`, `relativistic`, `fill_masses`, `mass_check`, `prior_covariance`, `set_nuclear_potential`) | `output_dir`, `channel_order`, `radius_override`, `rebuild_resources`, `relativistic`, `masses`, `thm_nuclide_masses`, `prior_covariance`, `potential_record`, `save_fit_format3` | `Project` edits of `<config>` keys (not built), potential records |
| the old segment fields and priors (`segment_values`, `prior_centres`, `CLASSIC`) | `save_fit_norms`, `apply_fit` | gone with formats 1-2 |
| `convert_azr` | `convert_azr`, `potential_record`, `project`, `save_fit_format3`, `validate` | `_azure2.convert_format3_text` (stays) |
| byte-for-byte layout of the old writer | `roundtrip`, `azr_layout`, `sav_fixed_params`, `zero_width_fixed` | deleted with the old formats; `Project.open(f).save(g)` keeps the bytes of a format-3 file (`project_test`) |

## What the engine still has to offer before the deletion

1. An edit of `<levels>` rows from Python (fitted values, a level added or removed, active flags), the counterpart of
   `parameter_value_edits` for level parameters: `save_fit`, `modelavg`, the level tests.
2. `<config>` keys as edits of `Project` (output directory, channel radius, relativistic flag, prior covariance).
3. The scattering-length term's rows from the objective.
4. The `<thm>` options as an edit binding.
