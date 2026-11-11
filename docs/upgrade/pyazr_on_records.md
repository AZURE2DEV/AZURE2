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

## Remaining uses in `pyazr/` (pyazr 2.32.0)

Every module below works on a format-3 project without `AzrModel` (the engine's edits of the nuclear part, the next
section, and `Project`'s record edits). On an older file each keeps calling `AzrModel` behind one test of the format
at its top (`project._is_format3`); that branch goes with the file.

| Where | Format 3 | Old formats (deleted with the file) |
|---|---|---|
| `azure2.save_fit` (`pyazr/snapshot.py`) | `write_format3`: level values by the engine's (group, level) and (pair key, l, 2s) (`level_value_edits`), THM coherent backgrounds by `set_thm_experiment(cbackground=)`, the potentials as records (`Project._replace_record_section`), the nuisances (`_write_parameter_values`); the verification copy's output directory by `set_config` | `write_old_format` |
| `azure2.scattering_length_chi2` | the model-data term's `scattering_length` rows (`model_data`) | `<scatteringLengths>` lines |
| `_HybridPotentials._write_potential_records` | the `<potentials>` body replaced by the session's records | refused (unchanged) |
| `bands.trimmed_model`, `bands._worker_setup` | extrapolation records removed and re-added in the order asked (`_trimmed_project`), `set_config(output_dir=)` | `_trimmed_old_format` |
| `angular.angular_distribution` | one `legendre` extrapolation per energy (`expansion a_over_a0`, `orders 0 n`), the reaction from `Project.reaction` | `clear_extrapolations` / `add_extrapolation` |
| `tabulate` | data sets, extrapolations, normalizations, energy scales and THM experiments removed (`_strip_format3`), one `sigma` data set per pair combination | `_strip_old_format` |
| `modelavg.write_averaged_azr`, `azr_parameter_slots` | the levels of `Project.levels()`, values by `_write_level_values` | `engine_level_keys` |
| `transform.levels_from_azr` | `Project.levels()` and `pairs()` (`_levels_from_project`; groups and levels as the engine numbers them) | the old level list |
| `nds` docstring | `Project.add_dataset` | — |
| `nuclei.merge_projects` | **not moved**: merging projects of several compound nuclei is the natural-target construct of group 6 (`natural_target` is its only test and still old); its format-3 form (records of several nuclei, `isotope` targets) waits for that group | `AzrModel` |
| `datasets._engine_module` | its own loader (works when loaded by path) | stays with `datasets.py` (old formats only) |
| `__init__` | exports `AzrModel`, `AzrLevel`, `AzrChannel`, `ScatteringLengthLine`, `PairMassCheck`, `convert_azr` | `convert_azr` moves to the converter's binding (`_azure2.convert_format3_text`); the rest are deleted with the file |
| `azure2._parsed("model")` | not used | `scattering_length_chi2` of an old file |

## Tests in `tests/pyazr/` (pyazr 2.32.0)

**Moved to the format-3 file** (they read the tree's project with `Project`, edit records and the nuclear part,
and run the same runs with the same expectations; none reads an original through `old_tree.py` any more):
`jacobian_skip`, `shift_first_call`, `jacobian_shift_columns`, `zero_width_level`, `formal_capture_mask`,
`physical_norm` (its scattering length a `datum` record: `scattering_length_chi2` from the model-data term),
`extrap_layout` (five `sigma` extrapolation records), `ec_integrals_formal` (`add_level`), `ec_relativistic`
(`set_config(relativistic=)`).  New: `project_nuclear_part` (every edit of the next section, the refusals, the bytes
kept, `save_fit` on 13N, `<thm>` on `examples/o18_lacognata2010`).

**Still on an original through `old_tree.py`** (format3_tree.md section 5), and why:

| test | what it tests | state |
|---|---|---|
| `roundtrip`, `sqrt_shift` (parts 1-2), `norm_group` (parts 1-2), `shift_fn` (reader part), `angular_coefficients` (reader part), `data_covariance` (reader part) | the old reader and writer of segment lines (`SegmentSet`, `add_data_segment`, byte-identical re-write) | deleted with the file; the format-3 counterpart is `project_test` (a file opened and saved keeps its bytes) and the record edits |
| `norm_group` (part 3), `sqrt_shift` (part 3), `shift_fn` (engine part), `prior_covariance`, `objective`, `screening`, `parameter_kinds` (screening part), `data_covariance` (engine part), `angular_coefficients` (engine part) | the engine on variants written as old segment / `<targetInt>` lines (norm groups, shift blocks, prior-centre rows, screening blocks) | to rewrite on records: the variants exist in format 3 (`normalization`, `energy_scale`, target `screening`; the checks of `tests/norm_groups`, `shift_fn`, `electron_screening` were rewritten so in groups 4a, 4b, 5), but the parameter layouts and names the tests assert are the old ones (a norm group lists one parameter in format 3, the old layout one per line) |
| `save_fit_norms` | `save_fit(norms=)` on old segment fields and `prior_centre` rows | the old-format half is deleted with the file; the format-3 half (`norms="fitted"` writes the nuisances with their priors) is `save_fit_format3` |
| `target_density`, `composite_gradient` | densityfit groups; a SUM of an E1 part and the full cross section | variants the converter refuses (`composite_gradient`'s SUM is a thick-target yield through the old `<targetInt>` line, which format 3 cannot state: "a thick target on isDiff 5 data"); a model decision is needed before they move |

**Using `AzrModel` on a format-3 file of the tree** (it still reads `<levels>` and `<config>` there): `channel_order`
(reorders the channel rows of a level: no `Project` call does that, by design), `inactive_lines`, `masses`
(`mass_check` / `fill_masses`: the table comparison would move to `Project.pairs` and `set_pair`), `hard_sphere_transform`,
`rebuild_resources`, `light_excitation` (the 32nd column read and written by the old class), `level_grouping`, `azr_layout`,
`sav_fixed_params`, `zero_width_fixed` (the old writer's bytes): to move or to delete with the file, one by one.

**Groups 6-8, still old in this tree** (left as they are): `natural_target`, `doppler`, `thickness` (group 6); of group 7
(converted 2026-11-11, `format3_tree.md` §5) the tests of AzrModel's THM API read the originals through `old_tree.py`
(`thm_experiment`, `thm_lineshape`, `thm_spectator_window`, `thm_distortion`, `thm_fixed_angle`, `thm_dw_vertex`,
`thm_spectator_angles`, `thm_coherent_background`, `thm_model_average`, `thm_vertex_reach`, `tabulate_thm`,
`extrapolation_frame`, `geometry_no_transform`, the 7Li_p_a parts of `roundtrip` and `save_fit_norms`), while
`thm_park`, `thm_adjoint`, `thm_nonbrune_grid` edit records (`tests/pyazr/thm_records.py`); `apply_fit`, `output_dir`, `radius_override`, `relativistic`, `scattering_length` (group 8 and the
rest).

## Edits of the nuclear part

The nuclear part of a format-3 project is what is not a record: `<levels>` (one row per channel of a level),
`<config>` (ten positional lines, then `key=value` lines), `<thm>` (global `key=value` lines and `experiment[<name>]`
lines) and, among the records, the `scattering_length` data of `<modelData>`. A person and an agent edit it through
`pyazr.Project`, as they edit records; the same calls are the operations of a GUI table (a row of the level table, a
cell of it, a key of a settings form).

**Who owns what.** The engine owns each grammar: it reads the section, says what each line is, makes the edited line
and the new line, and refuses what it would refuse on a run. Python holds no column position and no grammar: it hands
over names and numbers, applies the line replacements the engine returns to the lossless document (`Document`), and
the edit goes through `Project._apply` (the format-3 validator on the project as it would stand; a new error refuses
the edit, the project unchanged). A row, a `<config>` line or a `<thm>` line that the edit does not touch keeps its
bytes; a touched row keeps its layout (each new number right-aligned in the width of the token it replaces, as the
old `AzrModel` did, so a single edit does not reflow the block).

| what | engine (C++) | binding | `Project` |
|---|---|---|---|
| a `<levels>` row's text | `io/LevelRowText.h`: the names of the 34 columns, a row split into named tokens, a row with tokens replaced (beside `LevelRows.h`, the reader) | — | — |
| the levels | `experiment/LevelEdits.h`: the rows with their lines, grouped as `CNuc::Fill` groups them (`LevelTablesOfRows`, extended with the slot of each row), each level with its label; inactive rows grouped by consecutive ID | `level_table(doc)` | `levels()`, `pairs()` |
| level edits | `LevelEdits.h`: set a level's energy / fixed / active flags, a channel's width / fixed flag, a pair's radius, masses, separation energy (every row of the pair), add a level (the rows cloned from its J^π group, or from a row of each pair for a new group), remove a level, write fitted values by the engine's (group, level, channel) | `level_edit(doc, op, ...)`, `level_value_edits(doc, values)` | `set_level`, `set_channel`, `set_pair`, `add_level`, `remove_level`; `_write_level_values` (`save_fit`) |
| `<config>` | `io/ConfigSection.h`: the reader `Config::ReadConfigFile` uses (moved there from `src/Config.cpp`), each entry with its line, and the line of an edited entry; the table of keys | `config_entries(doc)`, `config_edit(doc, name, value)` | `config()`, `set_config(**entries)` |
| `<thm>` | `include/ThmBlockEdits.h` beside the reader (`ThmExperiment.h`): the line of a global option, an experiment's key set or removed in its line, an experiment removed; the block checked as the engine reads a format-3 one (members from the records, `Format3ThmLines`) | `thm_edit(doc, ...)`, `check_format3_thm(text)` | `thm()`, `set_thm(**options)`, `set_thm_experiment(name, **keys)`, `remove_thm_experiment(name)` |
| scattering lengths | the model-data term (`ModelData`): its rows with quantity `scattering_length` | `objective_terms_rwa()["model_data"]` (exists) | `azure2.scattering_length_chi2` sums them |

**Addressing.** A level is named by its label, `"<J^π> <E_x as written>"` (`"3/2- 3.5032"`; `#<k>` appended when two
levels share it), by `(J^π, E_x)` with E_x within 1 keV (the engine's rule for one level, `io::LevelGrouper`), or by
the `Level` that `levels()` returned. A channel by its pair key and (l, s); a pair by its key. `save_fit` names
levels by the engine's (group, level) indices of the fit parameters (`Parameter.jgroup`, `.level`) and channels by the
pair key, l and 2s, so a fitted value lands where the engine read it from (inactive rows take no index, as in
`CNuc::Fill`). A `<config>` entry by its name (`output_dir`, `checks_dir`, `amatrix`, the eight checks, and the keys as
the engine reads them: `relativistic`, `priorCovariance`, ...); a THM option by its key; an experiment by its name.

**What an edit returns.** The new warnings of the validator (a list of the engine's diagnostics, empty when there are
none); an error raises `FormatError` with the engine's diagnostics.  Values are numbers (written as the shortest
decimal that reads back as the same double) or, for `<config>` and `<thm>`, the words of the engine's grammar.

**Not offered** (no caller needs it): adding a pair (a new pair is a new row of a level, cloned from nothing: write it
by hand or with the converter), removing a channel from a group, renumbering level IDs, `<parameterSettings>` and
`<mcmc>` edits (their rows name parameters; `Project.rename` maps them).

## What the engine offers (built, pyazr 2.32.0)

1. Edits of `<levels>` rows (`experiment/LevelEdits.h`, `io/LevelRowText.h`): `Project.set_level`, `set_channel`,
   `set_pair`, `add_level`, `remove_level`, and `_write_level_values` for `save_fit` and `modelavg`.
2. `<config>` entries as edits (`io/ConfigSection.h`, the one reader, which `Config::ReadConfigFile` and
   `experiment::RelativisticFlag` now use): `Project.config`, `set_config`.
3. The scattering-length rows: the model-data term already returns them (`objective_terms_rwa()["model_data"]`).
4. `<thm>` options and experiments as edits (`include/ThmBlockEdits.h`, beside the reader; `Format3ThmLines`, which
   the engine's reading of a format-3 `<thm>` now uses): `Project.thm`, `set_thm`, `set_thm_experiment`,
   `remove_thm_experiment`.
5. The reaction label of two pair keys (`experiment::ReactionLabel`, which the converter now writes its labels with):
   `Project.reaction`.
