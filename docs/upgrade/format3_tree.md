# The test tree and the examples in format 3

Every project file of `tests/` and `examples/` is replaced, under its own name, by its format-3 conversion
(`project_format3.md` §8.5 item 6), and every `check.sh` that makes variants of a project edits format-3 records
instead of old keyword lines. This page is the plan and the state of that work: a successor continues from the
first group whose state is not "converted". The engine's old reading path stays intact until the whole tree is
converted, so the tree may be mixed at any commit.

## 1. Rules

- **Converting a file.** `tools/format3_tree_convert.sh <AZURE2> <path>...` (paths from the repository root) runs
  `AZURE2 --convert <path> <tmp> --to 3 --exact-grids --json-report …` at one thread in the file's folder (the
  converter reads the data files there), moves the original with `git mv` to `tests/convert/fixtures/<path>` and
  puts the conversion under the original name; side files `<stem>.f3<ext>` are added beside it. The reports go to
  `$REPORTS` (default `$WORK/format3_tree_reports`), not into the tree.
- **The originals** (`tests/convert/fixtures/<path from the repository root>`) are the converter's fixtures. They have
  no data files of their own: whatever reads one uses the data files of its home folder `<path>`'s directory.
  - `tools/convert_tree.sh` (and so the converted-tree gate `tools/format3_gate.sh`) copies each original back to its
    home in its copy of the tree before converting, and drops `tests/convert/` from the copy: the gate keeps
    converting the old files and comparing their outputs with the reference tag of the old files.
  - `convert_spec_equivalence` reads an original in place of the committed conversion of the same path;
    `convert_fixtures` reads its projects the same way (`Original` in `convert/tests/fixtures_test.cpp`).
  - `tools/runall.sh`, `tests/upgrade/plan_test.sh` and `tests/upgrade/predictions_test.py` leave `tests/convert/` out
    (no data files there). `tests/run_tests.sh` never looks below `tests/<name>/`.
  - pyazr tests that test the old file class (`pyazr/azrfile.py`, `pyazr/datasets.py`) on a tree project read its
    original through `tests/pyazr/old_tree.py` (`original(path)`, `copytree(src, dst)` with the originals laid over the
    copy). They are the old reader's tests; they go when that reader goes (§5).
- **A check that edited old lines** is rewritten to edit format-3 records, in this order of preference: (a) `sed` on
  one `key value` line of a record; (b) a committed variant file (format 3) when the edit is long; (c) the engine's
  tools (`AZURE2 --convert` of a variant that only the old grammar can state, `--describe`, `--validate`). It runs the
  same runs and expects the same numbers. An expected text that names old parameters (`segment_5_norm`) is updated in
  the same commit, with the class of `tools/format3_gate.py` (`CLASSES`) that explains the change.
- **Names.** The converter's names are kept (data sets after their data files, `tgt_targetInt<n>` for shared
  targets, `x<k>` for extrapolations); renaming is left for a later pass (it would move output and parameter names).
- **Per commit:** the touched directories' regression projects and `check.sh` (`tests/run_tests.sh` restricted to
  them), the ctests that read the converted files, `tests/upgrade/predictions.txt` regenerated where only data-set
  labels moved (`data 3:` → `data meyer_144:`; the tree of every line otherwise unchanged, checked with the diff
  filtered of labels), the consistency test, the build without warnings when code changed.

## 2. Measured before the work (all 185 convertible files converted at once, nothing else changed)

- Regression projects (`expected/chiSquared.out`): all pass on the converted files.
- `check.sh`: 51 of 121 pass unchanged, 70 fail (they edit old lines, read old sections, or copy a project without
  its side file).
- `ctest -LE slow`: 57 of 182 fail: 50 pyazr tests (the old file class on tree projects), `convert_fixtures`,
  `convert_spec_equivalence` (they read the tree), `park_chi2_fallback` (adds a `<targetInt>` line),
  `upgrade_plan_runtime` (edits old lines of `tests/transmission_yield`'s models), `fit_threads` (copies T2 without
  `prior_cov.f3.txt`), `upgrade_poison_results`, `upgrade_predictions` (data-set labels).
- Not converted by the converter (stay old, §4): `tests/angular_coefficients/models/bad_norm.azr`,
  `bad_test.azr`, `tests/electron_screening/bad.azr` (refusals of the old grammar), `tests/natural_target/models/
  M_thick.azr` (4 questions). Format 3 already: the 8 files of `ratio_to_rutherford`, `spin_orbit_potential`,
  `units_scales`.

## 3. Groups and state

A group is converted in one commit with every script and test that reads its files. Shared projects tie
directories together (`tests/13N/13N.azr` is read by 15 scripts and 9 pyazr tests, `tests/18O_p_a_thm` by 13 scripts
and 12 pyazr tests), so the groups follow the projects, not the directory names.

| group | projects converted | scripts and tests that read them | state |
|---|---|---|---|
| 1. self-contained | regression and model files whose checks edit nothing old: `13N_capture_ay`, `angle_mapping`, `capture_ay_coupling`, `convolution_unify` (not `bare`), `coulomb_nuclear_sign`, `data_file_comments`, `ec_e1m1_phase`, `ec_final_state_listed`, `ec_m1_exact`, `gamma_lab_frame`, `identical_dd_mott`, `identical_pp_mott`, `q_attenuation` (not `plain`), `ratio_frame` (not `ratio_data`), `reaction_rate` (not `n12C`, `p12C`), `secondary_gamma` (not `grad`, `mixing`), `shift_sqrt_subpoints`, `shift_subpoint_angle`, `sum_lab_frame` (`plain_data`, `sum_free`), `target_effect_ranges`, `transform_level_order`, the examples but `c12c12_tumino2018`, `natural_neon`, `o17_guardo2017_fit` (81 files) | their own `check.sh` (no edit of old lines; their `sed`/`awk` touch `<config>` lines or outputs only), `upos_pair_index`, `target_effect_shared`, `thm_beam_profile`, `thm_threshold_grid`, `thm_fixed_angle`, `failed_point`, `thm_park` (its example part) | converted |
| 2. 13N | `13N`, `mcmc_physical_norm/13N_fast`, `q_attenuation/models/plain`, `ec_threshold/models/ng`, `secondary_gamma/models/mixing`, `ec_m1_exact/models/m1`, `7Li_p_ay` | rewritten: `component_mapping` (new composites, see §4), `ec_final_state_ambiguous`, `ec_integrals_signature`, `ec_relativistic_signature`, `ec_shift`, `excitation_energy_column`, `kinematic_limits`, `mcmc_physical_norm` (and its `make_model.py`), `mcmc_snapshot_components`, `mixing_guard_lorder`, `normalizations_out` (part a), `progress_bars`, `straggling_bohr` (part d), `thick_target`, `total_capture_code` (§4); through the converter: `segline_formats`; unchanged and passing: `13N`, `level_ids`, `level_channel_order`, `nuisance_level_index`, `parameter_settings_nan`, `ec_m1_exact`, `ec_threshold`, `q_attenuation`, `secondary_gamma`; pyazr tests of the old reader on the original of 13N: `roundtrip`, `physical_norm`, `shift_first_call`, `target_density`, `jacobian_skip`, `jacobian_shift_columns`, `zero_width_level`, `formal_capture_mask`, `composite_gradient` | converted |
| 3. thick target | `thick_target_resonance` | rewritten (each model's records written by the check): `straggling_bohr`, `target_auto_windows` (but (f), §4), `target_thickness_exact`, `target_window_fallback`; unchanged: `shift_sqrt_subpoints` | converted |
| 4a. norm groups, prior centres | `norm_groups`, `15N_p_a` | rewritten: `norm_groups` (the follower's own flags cannot be written: the fixed-group case checks that no norm parameter exists), `normalizations_out` (part b), `prior_covariance` (data sets named `s1`, `s2`; parameters `s1.norm`, `s1.shift.offset`), `prior_centre` (centres in the `prior` of the field; the row refusals through the converter, §4); unchanged: `nuisance_prior`; pyazr tests on the originals: `norm_group`, `prior_covariance`, `save_fit_norms`, `extrap_layout`, `objective` (13N and shift_groups copies) | converted |
| 4b. shifts | `shift_fn`, `shift_groups`, `energy_shift_sqrt` (`7Li_p_a` moved to group 7: its readers are THM checks) | `shift_fn` (old `shiftfn` blocks, 7 variants and 6 refusals of the block grammar), `shift_fn_groups`, `shift_groups`, `shift_penalty_zero`, `shared_tables_race`, `energy_shift_sqrt`, `mcmc_classification` (part), pyazr `shift_fn`, `sqrt_shift`, `objective_terms` | to do: the checks build their variants from old segment lines with `shiftfn`/`shiftgroup`/`sqrtshift` blocks; each becomes `energy_scale` records (`function`, `param`) and `shift` fields; the refusals of the block grammar (short block, follower with its own function, sqrt and function on one line) are the converter's; the parameter names move (`segment_1_energy_shift_fn0` -> `<scale>.a0`) |
| 5. targets and apparatus | `electron_screening`, `screening_6li`, `target_density`, `convolution_window`, `convolution_unify/models/bare`, `data_covariance`, `reaction_rate/models/n12C`, `p12C` (read by `transmission_yield/make_models.py`), `secondary_gamma/models/grad`, `angular_coefficients`, `beam_profile_kernel`, `hard_sphere` | `electron_screening`, `screening_6li`, `target_density`, `data_covariance`, `doppler_fgm`, `transmission_yield`, `clone_composite`, `angular_coefficients`, `beam_profile_kernel`, `hard_sphere`, `cli_route_exit`, `extrap_sav_segments`; `park_chi2_fallback`, `upgrade_plan_runtime`; pyazr tests | to do |
| 6. natural targets | `natural_neon` (tests and example), `natural_target` | `natural_neon`, `segline_formats`, `shared_effect_pairs`, `shared_point_tables`, `total_optical`, `transmission_yield`, `natural_target`, `component_subpoint_angle`; pyazr tests | to do |
| 7. THM | `18O_p_a_thm`, `17O`, `7Li_p_a`, `6Li_d`, examples `c12c12_tumino2018`, `o17_guardo2017_fit`, `thm_coherent_background/toy` | the `thm_*` scripts, `18O_p_a_thm`, `level_at_threshold`, `excitation_energy_column`, `mcmc_classification`; pyazr THM tests | to do |
| 8. the rest | `upgrade/T1–T3`, `relativistic`, `scattering_length`, `channel_space`, `beta_delayed_alpha`, `identical_entrance_reaction`, `identical_pp_res`, `hybrid_potential`, `sum_lab_frame` (`sum_data`, `extrap`), `ratio_frame/models/ratio_data`, `ec_m1_exact/models/m1` | `upgrade`, `format3_engine`, `model_data`, `fit_threads`, `relativistic`, `scattering_length`, `phase_channel_spin`, `channel_space`, `beta_delayed_alpha`, `fixed_param_par`, `mass_separation_check`; `upgrade_poison_results`; pyazr tests | to do |
| 9. the refusal switch | — | §8.4 of `project_format3.md` | not started; only when nothing loads an old file |

## 3a. State at the last commit (verification)

92 of the 185 convertible files are converted (85 projects of `tests/`, 7 examples); 4 are not converted by the
converter and 8 are format 3 already (§2). The converted-tree gate on this tree (`REF=tab-trial8`, tag `gtree3a`)
gives the base's counts: 483 (456 formal) identical, 38 classified, 0 differ, 0 missing (the originals laid back by
`convert_tree.sh`). The tree's own run, `tab-tree3` (all projects × 3 bases), is the reference for the next agent: every
output of a converted project equals the gate's conversion of its original (`classify.py tab-gtree3a tab-tree3`: the 39
projects that differ are all unconverted ones, old against converted), and every output of an unconverted project
equals `tab-trial13` (the 41 projects that differ from `tab-trial13` are all converted ones, by the gate's classes).

## 4. Kept old on purpose, and checks that changed their model

| file or check | what | why |
|---|---|---|
| `total_capture_code` (a), (b) | an old line with exit key −1 and a data type that is not total capture, refused by the old engine | format 3 writes the two halves as one token (`[total]`), so the mismatch cannot be written. The converter converts such a line to `[total]` without a note (the old engine refused the file): a decision for the director (refuse in the converter, or accept). Part (d), the total-capture type with a real exit key, is checked through the converter. |
| `segline_formats` | the layouts of an old segment line (no shift, one shift, three fields; keyword blocks before `sample`) | the old segment-line reader is the converter's (`convert/include/azure2/convert/SegLine.h`): the check converts the original 13N and natural_neon with the line replaced and runs the conversion, the same expectations. |
| `component_mapping` | composites of the 13N model | two composites of the old model cannot be written in format 3 (a sum on a thick-target yield, `E_REACTION_NOT_ALLOWED`; an angle-integrated elastic term of a charged pair in a sum, no quantity). The new model: the composites in the check's header; the same assertions (outputs identical with and without the mapping, fewer level-matrix solutions, gradient = differences). |
| `progress_bars` | the SUM model | data set 6 (σ) instead of data set 1 (a thick-target yield, which takes one reaction in format 3). |
| `target_auto_windows` (f) | automatic windows and explicit ranges on one line (`"auto,0.48-0.50"`) | format 3's `fold_windows` is `auto[:tol[:margin]]` or a list of ranges, not both; the converter writes `auto` and drops the ranges with a note. The case is left out of the check: a decision for the director (a grammar for both, or none). |
| `target_window_fallback` (e) | the automatic per-point tolerance with no windows (old: an empty ranges token, then the tolerance) | format 3 states `auto_tolerance` only inside `fold_windows`, after a range or `auto`; the check writes one window over every point (`0-10 MeV auto_tolerance 0.001`), which gives the old outputs (compared). The converter drops such a tolerance **without a note** (`convert/src/BuildApparatus.cpp`, the windows branch requires ranges or auto): a defect for the director; `tests/reference/park_chi2_fallback_test.cpp` writes such a line (tolerance 1e-12). |
| `prior_centre` (unknown, typo, bad, nonpos) | `prior_centre` rows of `<parameterSettings>` that name no line, are misspelt or malformed, or put a norm's centre at 0 | the rows are the old grammar's; the converter refuses each with the old engine's reason (the check runs it). The format-3 counterpart of the last, `norm … fit prior 0 15%`, is **accepted** by `--validate` and the run (a prior of zero width, silently not applied: the chi2 is the classic one): a defect for the director (validation, out of this branch's files). |
| `norm_groups` (fixed leader) | a fixed group whose follower's own line says "free" | format 3 has one value per normalization record; the check now asserts that a fixed group has no norm parameter at all (the old file kept a fixed follower row at the leader's value). |
| `kinematic_limits` (c) | pair 1 with its particles swapped | the labels become `1H(12C,12C)1H`, `1H(12C,g)13N`: `p` is not accepted as a target nucleus in a label. |

Checks whose numbers were compared with the base's (`work/tree3/cmpbase.sh`: the base check on the base build and
the rewritten one on this branch, outputs diffed): identical for `thick_target`, `straggling_bohr`,
`segline_formats`, `normalizations_out`, `ec_relativistic_signature`, `ec_integrals_signature`, `ec_shift`,
`ec_final_state_ambiguous`, `kinematic_limits`, `excitation_energy_column`, `mixing_guard_lorder`,
`mcmc_physical_norm`, `mcmc_snapshot_components`, `target_thickness_exact`, `target_window_fallback`,
`normalizations_out`; identical but for the case left out for `target_auto_windows`, the names in the refusal
messages of `prior_covariance` (and two gradient agreements at the 1e-15 level), the fixed-leader line of
`norm_groups`, the converter's lines of `prior_centre`.

## 5. What still reads old files (blocks the refusal switch)

- The groups not yet converted (§3).
- `AZURE2 --convert` on purpose: `segline_formats`, `total_capture_code` (d), `ec_final_state_ambiguous` (e) convert an
  old model at run time; they keep working after the switch (the converter reads old files).
- pyazr tests of the old file class, on the originals (`tests/pyazr/old_tree.py`): those listed in §3. Each loads the
  original through `pyazr/azrfile.py` and runs the engine on it: they block the switch until pyazr's session reads
  format 3 and the old class is deleted (another agent's work).
- The model generators `tests/*/make_models.py` (29 of them; `make_model.py` of `mcmc_physical_norm`) write old files.
  Their committed outputs are converted; those run by a check (`transmission_yield`, `light_excitation`,
  `natural_target`, ...) are rewritten with their group. The others are developer tools that read a tree project and
  write old lines: to regenerate, point them at the original (`tests/convert/fixtures/...`, as
  `mcmc_physical_norm/make_model.py` now does) and convert the output; they are to be rewritten for format 3 or
  deleted (open).

## 6. User-guide pages that show old keyword lines (for a later documentation pass)

`docs/source/user_guide/segments.rst`, `experimental_effects.rst`, `chi_squared.rst`, `pyazr.rst`;
`docs/source/reference/data_formats.rst`, `command_line.rst`; `docs/source/examples/natural_neon.rst`;
`docs/source/theory/thm_implementation.rst`, `polarization_implementation.rst`; `docs/source/developer/
architecture.rst`, `coulomb_table.rst`, `edata_and_target_effects.rst`; `docs/source/getting_started/whats_new.rst`.
No test runs an `.azr` snippet of `docs/source/`.
