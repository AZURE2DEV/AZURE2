# pyazr on records: done (pyazr 3.0.0.dev0)

pyazr reads project format 3 only. `pyazr/azrfile.py` (`AzrModel`, a second grammar of the old formats in Python),
`pyazr/datasets.py` (`SegmentSet`, `TestSegmentSet`), `tests/pyazr/old_tree.py` and the old model generators of
`tests/doppler_fgm` and `tests/transmission_yield` (kept under `tests/convert/fixtures/`) are deleted. Everything
pyazr knows about a project comes from the engine: the records (`pyazr.Project`, the engine's `Document` and
validator), the nuclear part (`<levels>`, `<config>`, `<thm>`: the engine's edits, section "Edits of the nuclear
part" below), the description of the built project (`describe`) and the calculation session (`pyazr.azure2`).
The engine still reads old files (its refusal switch is the next branch, `format3_tree.md` §5).

## What replaced what

| before | now |
|---|---|
| `AzrModel.from_file`, `.write` | `Project.open`, `.save` (the lossless document: an unedited project gives its bytes back) |
| level and channel setters, `add_level`, `remove_level`, `deactivate_level`, `set_channel_radius`, masses | `Project.levels()`, `level()`, `pairs()`, `set_level`, `set_channel`, `set_pair(radius=, light_mass=, heavy_mass=, separation_energy=, light_excitation=, heavy_excitation=)`, `add_level`, `remove_level` (engine: `experiment/LevelEdits.h`) |
| `AzrChannel.light_excitation` (the 32nd column) | `ProjectPair.light_excitation`, `set_pair(light_excitation=)` (the engine writes the column; a row without it gets it only for a value other than 0) |
| `mass_check`, `fill_masses`, `PairMassCheck` | `pyazr.masses.mass_check(project)`, `masses.fill_masses(project, pairs=, separation=)`, `masses.PairMassCheck` (through `Project.pairs()` and `set_pair`) |
| `engine_level_keys`, `apply_fit` | `ProjectLevel.group` / `.level` (the engine's numbering); `pyazr.snapshot.place_level_values` (save_fit: a level by (group, level), a channel by pair key, l, 2s) |
| segment lines, `add_data_segment`, `<targetInt>`, `<abundances>`, `<scatteringLengths>` | records: `dataset`, `extrapolation`, `target` (isotopes, density, screening, Doppler), `normalization`, `energy_scale`, `datum` (`Project.add/set/remove/rename`) |
| `<thm>` helpers of the old class | `Project.thm()`, `set_thm`, `set_thm_experiment(name, datasets=, ...)`, `remove_thm_experiment`; a data set's `thm` and `thm_weight` fields |
| `relativistic`, `prior_covariance`, `set_output_dir` | `Project.config()`, `set_config(relativistic=, priorCovariance=, output_dir=)` |
| `nuclei.merge_projects` on old files | `merge_projects(projects, out, *, isotopes=, target=, isotope_sum=, shared=)` on format 3: the `<levels>` rows merged by the engine (`experiment::MergeLevelRows`, binding `merge_level_rows`), records copied (a name taken by an earlier project renamed `<name>_<k>`), a `target` with `isotope` fields, shared data sets as one sum on that target; returns `MergeReport(pair_map, record_map, warnings)` |
| `convert_azr` in `azrfile.py` (default `to=2`) | `pyazr/convert.py`, default `to=3` (the engine's converter) |
| `datasets.fitted_norms` | `pyazr/outputs.py` |
| the old-format branch of `save_fit`, `scattering_length_chi2`, `trimmed_model`, `angular_distribution`, `tabulate`, `write_averaged_azr`, `levels_from_azr`, `datasets` / `extrapolations` | removed: each works on the format-3 project only |
| an old file handed to pyazr | one refusal, `pyazr.project.require_format3`: `FormatError` with the engine's message (`convert::NotFormat3Message`, binding `not_format3_message`), naming `AZURE2 --convert OLD NEW --to 3`; called by a session and `Project.open`, and through them by every function that opens a project |
| `load_sav`'s own name check | the engine's rule (`experiment::MatchParameterFileNames`, binding `match_parameter_file_names`): W_PARAMETER_FILE_NAMES warns, E_PARAMETER_FILE_NO_MATCH refuses |

## The tests of `tests/pyazr/` that used the old class

Moved: the test reads the tree's format-3 project with `Project`, makes the same variants by record and nuclear-part
edits and asserts the same committed numbers. Converter: it calls `pyazr.convert_azr` on an original kept in
`tests/convert/fixtures/`. Deleted: it tested the old reader or writer.

| test | state | notes |
|---|---|---|
| `roundtrip` | deleted | counterparts `project_test` (every format-3 file gives its bytes; CRLF), `project_nuclear_part_test` (an edit changes only its tokens), `azr_layout` |
| `azr_layout` | moved | every format-3 file of `tests/` byte for byte; CRLF, no-blank and interleaved variants; `set_level` changes only the level's rows; `remove_level` keeps the line ending |
| `level_grouping` | moved | `Project.levels()` on interleaved rows |
| `channel_order` | moved | the variants made by moving and deleting whole rows in the test (Project offers no reordering) |
| `inactive_lines` | moved | `add_level`, `set_level(active=False)`, rows moved in the test; `remove_level` of the inactive level |
| `masses` | moved | `masses.mass_check` / `fill_masses` on a Project |
| `light_excitation` | moved | `set_pair(light_excitation=)` on the nuclear part alone; the 35-token refusal is `io/tests/level_row_text_test.cpp`'s |
| `hard_sphere_transform`, `rebuild_resources`, `sav_fixed_params`, `zero_width_fixed` | moved | `sav_fixed_params` also checks the engine's name match; `zero_width_fixed` (c) is save_fit's placement |
| `apply_fit` | moved | `place_level_values` with synthetic Parameter/Pair (engine number vs key, orphan refused, fixed unlisted skipped) and the snapshot round trip |
| `model_average` | moved | `azr_parameter_slots`, `write_averaged_azr` on a Project |
| `natural_target` | moved | `merge_projects` on format 3, two targets, abundances, Jacobian, save_fit |
| `doppler`, `thickness`, `parameter_kinds` | moved | the tree's format-3 generators on the tree's n12C; the record grammar is their checks' |
| `norm_group`, `sqrt_shift`, `shift_fn`, `angular_coefficients`, `data_covariance` | moved | the reader parts deleted (counterpart: the record fields, `project_test`) |
| `prior_covariance`, `objective`, `screening`, `relativistic`, `radius_override`, `output_dir`, `scattering_length` | moved | |
| `target_density`, `composite_gradient` | moved | the director's decisions (below) |
| `save_fit_norms` | moved | the old half deleted (counterpart `save_fit_format3`) |
| `thm_experiment`, `thm_lineshape`, `thm_spectator_window`, `thm_distortion`, `thm_fixed_angle`, `thm_dw_vertex`, `thm_spectator_angles`, `thm_coherent_background`, `thm_model_average` (and `scripts/thm_model_average.py`), `thm_vertex_reach`, `tabulate_thm`, `extrapolation_frame`, `geometry_no_transform`, `thm_profiled`, `thm_vertex_memo_sessions`, `thm_nuclide_masses` | moved | `Project.thm()`, `set_thm_experiment(datasets=)`, `tests/pyazr/thm_records.py` |
| `convert_azr` | converter | to=2 and to=3 on the originals; the refusal of an old file by every entry point |
| `validate`, `potential_record`, `save_fit_format3`, `session_records`, `composite_terms` (2b) | converter | they read the originals directly; the parts that ran a session on an old file run it on the conversion |
| `azr_parse_cache` | moved | the tree's `channel_space` plus an E1 data set |

## Coverage lost

- `target_density`: two different targets tied to one density, and the group refusals of the old lines (format 3 has
  one target record per density).
- `composite_gradient`: the thick-target yield of a sum (E1 part + full cross section), which format 3 cannot state; a
  `yield` of a sum on an `infinite` target is refused by the validator as written today, so it is not tested.
- `scattering_length`: an inactive constraint line (format 3 has no inactive datum).
- `norm_group`, `prior_covariance`: the follower line's own ignored norm (not writable in format 3).
- `save_fit_norms`: the `<parameterSettings>` prior_centre rows (gone; a prior lives in its nuisance).
- The byte layout of old segment lines, `<thm>` aliases and the old writer's number format (the grammar is gone).
- THM: some refusals come when a session starts, not at the edit (`Project` accepts the text, the engine refuses
  it); an `extrapolation` without a frame is refused rather than defaulting to lab.
- `light_excitation`: setting the excitation back to 0 writes `0` in the column rather than dropping it.

## Edits of the nuclear part

The nuclear part of a format-3 project is what is not a record: `<levels>` (one row per channel of a level),
`<config>` (ten positional lines, then `key=value` lines) and `<thm>` (global `key=value` lines and
`experiment[<name>]` lines). The engine owns each grammar: it reads the section, says what each line is, makes the
edited line and the new line, and refuses what it would refuse on a run. Python holds no column position and no
grammar: it hands over names and numbers, applies the line replacements the engine returns to the lossless document
(`Document`), and the edit goes through `Project._apply` (the validator on the project as it would stand; a new
error refuses the edit, the project unchanged).

| what | engine (C++) | binding | `Project` |
|---|---|---|---|
| a `<levels>` row's text | `io/LevelRowText.h` | — | — |
| the levels | `experiment/LevelEdits.h` (`ReadLevelRowTable`) | `level_table(doc)` | `levels()`, `level()`, `pairs()` |
| level edits | `LevelEdits.h`: `SetLevel`, `SetChannel`, `SetPair` (radius, masses, separation energy, excitation energies), `AddLevel`, `RemoveLevel`, `LevelValueEdits`, `MergeLevelRows` | `set_level_edits`, `set_channel_edits`, `set_pair_edits`, `add_level_edits`, `remove_level_edits`, `level_value_edits`, `merge_level_rows` | `set_level`, `set_channel`, `set_pair`, `add_level`, `remove_level`; `_write_level_values` (save_fit); `nuclei.merge_projects` |
| `<config>` | `io/ConfigSection.h` | `config_entries`, `config_edits` | `config()`, `set_config(**entries)` |
| `<thm>` | `include/ThmBlockEdits.h` | `thm_entries`, `thm_option_edits`, `thm_experiment_edits`, `check_format3_thm` | `thm()`, `set_thm`, `set_thm_experiment`, `remove_thm_experiment` |

**Addressing.** A level is named by its label, `"<J^π> <E_x as written>"` (`"3/2- 3.5032"`; `#<k>` appended when
two levels share it), by `(J^π, E_x)` with E_x within 1 keV, by the `ProjectLevel` that `levels()` returned, or by
the engine's `(group, level)`. A channel by its pair key and (l, s); a pair by its key; a `<config>` entry by its
name; a THM option by its key; an experiment by its name. An edit returns the validator's new warnings; an error
raises `FormatError` with the engine's diagnostics.

**Not offered** (by design): reordering the rows of a level, adding a pair, removing a channel from a group,
renumbering level IDs, `<parameterSettings>` and `<mcmc>` edits (their rows name parameters; `Project.rename` maps
them).

**Bindings with no Python caller since pyazr 3** (kept; candidates for removal with the engine's old reading path):
`check_thm_block`, `read_thm_block`, `thm_option`, `thm_spin_text`, `thm_global_optical_for`, `read_thm_weight_table`
(`api/src/pybind_thm.cpp`), `parse_segment_line`, `parse_extrapolation_line`, `parse_target_effect_line`,
`LevelGrouper` (`api/src/pybind_api.cpp`).
