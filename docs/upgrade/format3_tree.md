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
| 1. self-contained | regression and model files whose checks edit nothing old: `13N_capture_ay`, `6Li_d`, `angle_mapping`, `capture_ay_coupling`, `convolution_unify` (not `bare`), `coulomb_nuclear_sign`, `data_file_comments`, `ec_e1m1_phase`, `ec_final_state_listed`, `ec_m1_exact`, `gamma_lab_frame`, `identical_dd_mott`, `identical_pp_mott`, `q_attenuation` (not `plain`), `ratio_frame` (not `ratio_data`), `reaction_rate` (not `n12C`, `p12C`), `secondary_gamma` (not `grad`, `mixing`), `shift_sqrt_subpoints`, `shift_subpoint_angle`, `sum_lab_frame` (`plain_data`, `sum_free`), `target_effect_ranges`, `transform_level_order`, the examples but `c12c12_tumino2018`, `natural_neon`, `o17_guardo2017_fit` (82 files) | their own `check.sh` (no edit of old lines; their `sed`/`awk` touch `<config>` lines or outputs only), `upos_pair_index`, `target_effect_shared`, `thm_beam_profile`, `thm_threshold_grid`, `thm_fixed_angle`, `failed_point`, `thm_park` (its example part) | converted |
| 2. 13N | `13N`, `mcmc_physical_norm/13N_fast`, `q_attenuation/models/plain`, `ec_threshold`, `secondary_gamma/models/mixing` | `13N`, `component_mapping`, `ec_final_state_ambiguous`, `ec_integrals_signature`, `ec_relativistic_signature`, `ec_shift`, `kinematic_limits`, `mixing_guard_lorder`, `normalizations_out`, `progress_bars`, `segline_formats`, `straggling_bohr`, `thick_target`, `total_capture_code`, `level_ids`, `level_channel_order`, `nuisance_level_index`, `parameter_settings_nan`, `mcmc_*`; pyazr tests on 13N | to do |
| 3. thick target | `thick_target_resonance` | `straggling_bohr`, `target_auto_windows`, `target_thickness_exact`, `target_window_fallback` | to do |
| 4. groups, shifts, priors | `norm_groups`, `shift_fn`, `shift_groups`, `energy_shift_sqrt`, `15N_p_a`, `7Li_p_a` | `norm_groups`, `normalizations_out`, `prior_covariance`, `prior_centre`, `nuisance_prior`, `shift_*`, `shared_tables_race`, `energy_shift_sqrt`; pyazr tests | to do |
| 5. targets and apparatus | `electron_screening`, `screening_6li`, `target_density`, `convolution_window`, `convolution_unify/models/bare`, `data_covariance`, `reaction_rate/models/n12C`, `p12C` (read by `transmission_yield/make_models.py`), `secondary_gamma/models/grad`, `angular_coefficients`, `beam_profile_kernel`, `hard_sphere` | `electron_screening`, `screening_6li`, `target_density`, `data_covariance`, `doppler_fgm`, `transmission_yield`, `clone_composite`, `angular_coefficients`, `beam_profile_kernel`, `hard_sphere`, `cli_route_exit`, `extrap_sav_segments`; `park_chi2_fallback`, `upgrade_plan_runtime`; pyazr tests | to do |
| 6. natural targets | `natural_neon` (tests and example), `natural_target` | `natural_neon`, `segline_formats`, `shared_effect_pairs`, `shared_point_tables`, `total_optical`, `transmission_yield`, `natural_target`, `component_subpoint_angle`; pyazr tests | to do |
| 7. THM | `18O_p_a_thm`, `17O`, `7Li_p_ay`, examples `c12c12_tumino2018`, `o17_guardo2017_fit`, `thm_coherent_background/toy` | the `thm_*` scripts, `18O_p_a_thm`, `level_at_threshold`, `excitation_energy_column`, `mcmc_classification`; pyazr THM tests | to do |
| 8. the rest | `upgrade/T1–T3`, `relativistic`, `scattering_length`, `channel_space`, `beta_delayed_alpha`, `identical_entrance_reaction`, `identical_pp_res`, `hybrid_potential`, `sum_lab_frame` (`sum_data`, `extrap`), `ratio_frame/models/ratio_data`, `ec_m1_exact/models/m1` | `upgrade`, `format3_engine`, `model_data`, `fit_threads`, `relativistic`, `scattering_length`, `phase_channel_spin`, `channel_space`, `beta_delayed_alpha`, `fixed_param_par`, `mass_separation_check`; `upgrade_poison_results`; pyazr tests | to do |
| 9. the refusal switch | — | §8.4 of `project_format3.md` | not started; only when nothing loads an old file |

## 4. Kept old on purpose

| file | why | where it is tested |
|---|---|---|
| (none yet) | | |

## 5. What still reads old files (blocks the refusal switch)

To be filled as the groups are converted.

## 6. User-guide pages that show old keyword lines (for a later documentation pass)

`docs/source/user_guide/segments.rst`, `experimental_effects.rst`, `chi_squared.rst`, `pyazr.rst`;
`docs/source/reference/data_formats.rst`, `command_line.rst`; `docs/source/examples/natural_neon.rst`;
`docs/source/theory/thm_implementation.rst`, `polarization_implementation.rst`; `docs/source/developer/
architecture.rst`, `coulomb_table.rst`, `edata_and_target_effects.rst`; `docs/source/getting_started/whats_new.rst`.
No test runs an `.azr` snippet of `docs/source/`.
