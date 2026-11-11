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
| 4b. shifts | `shift_fn`, `shift_groups`, `energy_shift_sqrt`, `pyazr/data/shift_fn_groups` (`7Li_p_a` moved to group 7: its readers are THM checks) | rewritten (each variant's `<data>` records written by the check, `energy_scale` records with `offset`, `function`, `param`; parameters `G.offset`, `G.offset[<data set>]`, `G.a0`, `<data set>.shift.offset`, `<data set>.shift.sqrt`): `shift_fn`, `shift_fn_groups`, `shift_groups`, `shift_penalty_zero`, `shared_tables_race`, `energy_shift_sqrt` (fit part), `mcmc_classification` (its energy-shift runs); the refusals of the old block grammar through the converter (§4); pyazr tests of the old reader on the originals: `shift_fn`, `sqrt_shift` | converted |
| 5. targets and apparatus | `electron_screening` (not `bad`), `screening_6li`, `target_density` (not `bad`), `convolution_window`, `convolution_unify/models/bare`, `data_covariance`, `secondary_gamma/models/grad`, `angular_coefficients` (not `bad_norm`, `bad_test`), `beam_profile_kernel`, `hard_sphere` | rewritten: `electron_screening` (fit models' `screening` field; part 5, the old block without its numbers, through the converter: `bad.azr` moved to the fixtures), `screening_6li`, `target_density` (`<target>.areal_density`; the group member's own density and `bad.azr`, densityfit without target integration, through the converter), `data_covariance` (refusals on records; `verify.py` reads the data-set records), `clone_composite` (its `make_models.py` writes format 3), `beam_profile_kernel` (its two extrapolations in `extrapolations.txt`, the records the converter writes for the old test lines), `hard_sphere`; unchanged and passing: `convolution_window`, `convolution_unify`, `secondary_gamma`, `angular_coefficients`, `cli_route_exit`, `extrap_sav_segments`; pyazr tests of the old reader on the originals: `data_covariance`, `angular_coefficients`, `screening`, `parameter_kinds` (screening_6li part). **Left old:** `reaction_rate/models/n12C.azr`, `p12C.azr`: the generators `doppler_fgm/make_models.py`, `transmission_yield/make_models.py` (and through them `pyazr` `doppler`, `thickness`, `parameter_kinds`, `upgrade_plan_runtime`, `upgrade_poison_results`) build their models from these files' old lines; they move with the rewrite of those two generators (group 6, with `natural_neon`, which `transmission_yield` also reads) | converted but n12C, p12C |
| 6. natural targets | `natural_neon` (tests and example; the example's `natural_neon_fit.sav` renamed to the format-3 parameter names, outputs equal), `natural_target` (not `M_thick`; `bad_key`, `bad_sample` to the fixtures: the converter refuses them), `reaction_rate/models/n12C.azr`, `p12C.azr` (14 files) | rewritten: `natural_neon` (the extrapolations' `beam`/`quadrature` removed instead of the `<targetInt>` lines switched off), `natural_target` (parameters `natural.isotope.12C` ...; the shared-norm case writes a `normalization` record; the prior on the `isotope` field; `bad_key`, `bad_sample` through the converter), `shared_effect_pairs` (beam records), `component_subpoint_angle` (extrapolation records), `doppler_fgm` and `transmission_yield` with their generators `make_models.py` (format-3 records; the old generators kept in `tests/convert/fixtures/tests/<dir>/make_models.py` for the pyazr tests of the old reader), `natural_target/make_models.py` deleted (its README says how the models were made); unchanged and passing: `segline_formats`, `shared_point_tables`, `reaction_rate`, `angle_mapping`, `upgrade_plan_runtime`, `upgrade_poison_results`; `total_optical` left old here, format 3 since the gaps branch (§4); pyazr tests of the old reader on the originals: `natural_target`, `doppler`, `thickness`, `parameter_kinds` (the old generators); `pyazr_predictions` labels a format-3 project's data sets by name | converted (`M_thick` and `total_optical` since the gaps branch, §4) |
| 7. THM | `18O_p_a_thm`, `17O`, `7Li_p_a`, `6Li_d`, examples `c12c12_tumino2018`, `o17_guardo2017_fit`, `thm_coherent_background/toy` | the `thm_*` scripts, `18O_p_a_thm`, `level_at_threshold`, `excitation_energy_column`, `mcmc_classification`; pyazr THM tests | to do. Measured (2026-11-11, converted and reverted): `17O` alone converts with 3 notes and its regression passes, but besides `pyazr_azr_layout`, `pyazr_convert_azr` (both edit or convert the old file: point them at the original) its readers `thm_narrow_fold`, `thm_energy_shift`, `thm_rwa_grid` (checks) and `pyazr_thm_nonbrune_grid`, `pyazr_extrapolation_frame` fail on the conversion and must be rewritten with it (`thm_park`, `pyazr_level_grouping` pass); predictions.txt moves its two lines to `data o17_guardo`, `extrapolation x1` |
| 8. the rest | `upgrade/T1–T3` (`T1`, `T1_start`, `T2`, `T2_true`, `T3`, `T3_start`), `relativistic` (14 models), `scattering_length` (5), `channel_space`, `beta_delayed_alpha/16N`, `identical_entrance_reaction` (`be8`, `c12c12`, `dd_pt`), `identical_pp_res`, `hybrid_potential`, `sum_lab_frame` (`sum_data`, `extrap`), `ratio_frame/models` (`ratio_data`, `fixed_angle`) (36 files; `ec_m1_exact/models/m1` was converted with group 2) | rewritten: `upgrade` (`golden/golden.py` compares the golden files of the old originals as the format-3 project writes them: param.par through the converter's name map, the T3 classes of the scattering lengths; `golden/T2/regen_compare.py` the class "follower listed"; `recover.py` sets the `<modelData>` scattering lengths), `format3_engine` (the old T2 from the fixtures; the committed T2 = its conversion), `model_data`, `fit_threads`, `scattering_length` (the `<modelData>` data, refusals on records), `phase_channel_spin` (`s` required; the old segment without "spin" through the converter), `channel_space`, `beta_delayed_alpha` (but (e), §4), `fixed_param_par`, `relativistic` (g), `cli_route_exit`, `extrap_sav_segments`, `level_at_threshold` (its c12c12 part), `hybrid_potential/check_per_pair.py`, `scripts/bench_upgrade.sh`, `tests/reference/park_chi2_fallback_test.cpp` (the beam spread and window as records, project text only); the generators of T1–T3 (`upgrade/lib/upgrade_project.py`, `T1/generate.py`, `T3/generate.py`, `T2/make_t2.py`, `T2/check_recovery.py`) write format 3 through the converter and give the committed files byte for byte; pyazr tests that run the projects edit records: `hybrid_derivative`, `potential_sessions`, `api_model_value`, `composite_terms` (but its part 2b), `sav_readers`, `output_penalties`, `predictions`, `upgrade_unitarity`, `upgrade_breit_wigner`; pyazr tests of the old reader or the converter on the originals: `scattering_length` (its block part), `azr_parse_cache`, `potential_record`, `save_fit_format3` (part 2), `validate`, `convert_azr`, `session_records` (its old T2), `composite_terms` (2b); unchanged and passing: `cli_eof`, `fit_call_limit`, `mcmc_park_start`, `park_formalism`, `identical_entrance_reaction`, `mass_separation_check`, `check_files_data`, `failed_point`, `sum_lab_frame`, `ratio_frame`, `spin_orbit_potential` | converted |
| 9. the refusal switch | — | §8.4 of `project_format3.md` | not started; only when nothing loads an old file |

## 3a. State at the last commit (verification)

176 of the 183 convertible files are converted (groups 1, 2, 3, 4a, 4b, 5, 6, 8). Left: group 7 (THM: `18O_p_a_thm`,
`17O`, `7Li_p_a`, `6Li_d`, `thm_coherent_background/toy`, examples `c12c12_tumino2018`, `o17_guardo2017_fit`). Not
converted on purpose: `angular_coefficients/models/bad_norm.azr`, `bad_test.azr` (old, in place; `natural_target/models/
M_thick.azr` is format 3 since the gaps branch, §4), and the old-grammar refusals `electron_screening/bad.azr`, `target_density/bad.azr`,
`natural_target/models/bad_key.azr`, `bad_sample.azr` (in the fixtures only; the converter refuses each, the last two
since group 6, so they no longer count as convertible: the spec equivalence gives 198 files, 183 converted, 183 equal,
0 recorded differences, 0 differ, 15 not converted). `total_optical` builds its models in format 3
(§4, since the gaps branch). 8 files were format 3 already. Groups 6 and 8 were converted on two branches and merged; the verification of the
merged tree is in the director's entry of `HANDOFF.md`.

Verification of group 8 (branch `feature/upgrade-4-tree-format3-4`): the converted-tree gate (`REF=tab-trial8`, tag
`gtree3d`) gives the base's counts, 483 (456 formal) identical, 38 classified, 0 differ, 0 missing (after
`tools/convert_tree.sh` drops the committed `prior_cov.f3.txt` before converting T2's original; the first run had T2
missing). The tree's own run, `tab-tree3d` (all projects x 3 bases), is the reference for the next agent: against
`tab-tree3b` every project not converted here is identical; the converted ones differ only by listed classes
(param.par names, the step and follower rows of `hybrid_potential` and T2, T2's follower listed in normalizations.out,
the scattering lengths as model data in chiSquared.out / scatteringLengths.out and the new modelData.out of
`be7p_fit`, `he3he4_fit`, T3, T3_start) and run.log. `ctest -LE slow` 183/183, `ctest -L slow` (OMP 2) 10/10; T1-T3,
13N, identical_pp_res and two relativistic models at 1 and 4 threads byte-identical; the T2 benchmark at 1 thread:
calculation chi2 10899.8, fit of 200 calls 691.608 (the base's), its output files identical to the base's fit of the
original but for param.par/param.sav names, the follower rows (`D.value[s03_diff150]`, `G.offset[s03_diff150]`) and
normalizations.out's follower row (classes "follower", "follower listed"). `tests/run_tests.sh` whole: 145 passed,
1 failed (`spin_orbit_potential`: its `unitarity.py` called the T1–T3 helper's old `write`; fixed in the last commit and
the directory rerun: passed).

Verification of group 6 (branch `feature/upgrade-4-tree-format3-3`): the converted-tree gate (`REF=tab-trial8`, tag `gtree3c`) gives the base's counts: 483 (456 formal) identical, 38
classified, 0 differ, 0 missing. The tree's own run, `tab-tree3c-1` (all projects × 3 bases), is the reference for
the next agent: against `tab-tree3b` every file is identical but natural_neon's `param.par` (parameter names, values
equal) and its `AZUREOut_sample_seg=*.out` headers (the sample's name for "1"), and two files fewer (the `run.log` of
`bad_key`, `bad_sample`, moved to the fixtures; `run.log` is not compared).

## 4. Kept old on purpose, and checks that changed their model

| file or check | what | why |
|---|---|---|
| `total_capture_code` (a), (b) | an old line with exit key −1 and a data type that is not total capture, refused by the old engine | format 3 writes the two halves as one token (`[total]`), so the mismatch cannot be written. Director's decision: the converter **refuses** such an active line with the old engine's reason and a non-zero exit (`convert/src/BuildData.cpp`; `convert_constructs`; the check gives both old lines to `AZURE2 --convert`). Part (d), the total-capture type with a real exit key, converts to `[total]`. |
| `segline_formats` | the layouts of an old segment line (no shift, one shift, three fields; keyword blocks before `sample`) | the old segment-line reader is the converter's (`convert/include/azure2/convert/SegLine.h`): the check converts the original 13N and natural_neon with the line replaced and runs the conversion, the same expectations. |
| `component_mapping` | composites of the 13N model | two composites of the old model cannot be written in format 3 (a sum on a thick-target yield, `E_REACTION_NOT_ALLOWED`; an angle-integrated elastic term of a charged pair in a sum, no quantity). The new model: the composites in the check's header; the same assertions (outputs identical with and without the mapping, fewer level-matrix solutions, gradient = differences). |
| `progress_bars` | the SUM model | data set 6 (σ) instead of data set 1 (a thick-target yield, which takes one reaction in format 3). |
| `target_auto_windows` (f) | automatic windows and explicit ranges on one line (`"auto,0.48-0.50"`) | a construct format 3 does not have (director's decision: no new grammar). `fold_windows` is `auto[:tol[:margin]]` or a list of ranges; the converter writes `auto`, drops the ranges and says so in a note of the conversion report. The case stays out of the check. |
| `target_window_fallback` (e) | the automatic per-point tolerance with no windows (old: an empty ranges token, then the tolerance) | format 3 states `auto_tolerance` only inside `fold_windows`; the check writes one window over every point (`0-10 MeV auto_tolerance 0.001`). The converter now writes the same for such an old line (`0-1e+06 MeV auto_tolerance <tol>`, with a note) instead of dropping the tolerance; its outputs equal the old line's on `identical_pp_res` with a convolution at tolerances 1e-12, 1e-3, 3e-2 (`convert_constructs` case). |
| `shift_fn` (6), `shift_fn_groups` (4) | refusals of the old block grammar: a follower of a shift group or of a shift-function group with a shiftfn of its own, a line in both kinds of group, a shift-function group whose first line has no shiftfn, a short block, sqrtshift with shiftfn | the check gives the old lines to `AZURE2 --convert`: the converter refuses each with the old engine's reason (the follower cases and the line in both groups were converted with the block dropped, or with a note, before this branch: a file the old engine refused ran), except sqrtshift with shiftfn, which converts to a data set with a sqrt term and a scale with a function and is refused by the run with the old message ("two descriptions of one energy calibration"). The refusals the format-3 records can state (function plus sqrt term, undeclared or unused parameter, a term with no effect) are run on records. |
| `target_density` (5), `electron_screening` (5) | old `<targetInt>` blocks the old engine refused (densityfit on a line without target integration; a screening block without its numbers) | kept as fixtures only (`tests/convert/fixtures/tests/<dir>/bad.azr`); the check converts them and expects the converter's refusal with the old reason. The converter used to convert the first one (a refused file ran with the density fixed): it now refuses a densityfit / densitygroup the old engine refused (`BuildApparatus.cpp`, `convert_constructs`). |
| `shift_groups` (1), `shift_fn` (3) | a follower whose own shift value is ignored | a data set naming an `energy_scale` has no value of its own: the check compares the record's two data sets with the same shift given directly. |
| `prior_centre` (unknown, typo, bad, nonpos) | `prior_centre` rows of `<parameterSettings>` that name no line, are misspelt or malformed, or put a norm's centre at 0 | the rows are the old grammar's; the converter refuses each with the old engine's reason (the check runs it). The format-3 counterpart of the last, `norm … fit prior 0 15%` (a centre at 0 or below), is refused by validation and so by the run with `E_PRIOR_CENTRE` (director's decision; the check's case "zero"). |
| `norm_groups` (fixed leader) | a fixed group whose follower's own line says "free" | format 3 has one value per normalization record; the check now asserts that a fixed group has no norm parameter at all (the old file kept a fixed follower row at the leader's value). |
| `beta_delayed_alpha` (e) | a Gaussian resolution (`<targetInt>` line with a convolution, 0.02 MeV) on a β-delayed spectrum; the old engine ran it | **format 3 since `feature/upgrade-4-format3-gaps`** (director's option 1): a detector's `energy_resolution` on a β-delayed spectrum (no beam) is the engine's Gaussian fold on the measured particle energy (`BuildEffects.cpp` `DetectorInto`; other data without a beam profile stay `E_NOT_BUILT`, the c.m. frame too); the converter writes the old line on β-delayed lines as that detector. The part runs on records; its outputs equal the converted old line's file for file but param.par (chi2 758.095, the old run's). A beam (`spread`) and a target with energy loss on β-delayed data stay refused (parts (b)–(d)). |
| `phase_channel_spin` (d) | a phase-shift segment without "spin" whose (J, l) has two channel spins (the old engine refused it) | format 3 requires `s` of a phase shift (`E_FIELD_REQUIRED`, checked); the old line through `AZURE2 --convert` is refused with the old engine's reason (checked on the original fixture). |
| `scattering_length` (7), `model_data` | the old `<scatteringLengths>` block's refusals and its own chi² line | the format-3 data (`datum … quantity scattering_length`) are refused on records (`E_NOT_A_CHANNEL`, `E_UNKNOWN_PAIR`, `E_FIELD_REQUIRED` for an error or an `imag_error`); the totals are `Total-Model-Data-Chi-Squared` and `modelData.out` (class "scattering lengths as model data"); `model_data` no longer appends an old block beside its `<modelData>` and checks the datum's own χ² in `modelData.out` instead of `Total-Scattering-Length-Chi-Squared` (rows printed to 5 digits: tolerance 1e-4). |
| `hybrid_potential/check_per_pair.py` (2) | a pair without a section of its own inherits the old block's default | format 3 has no default: one record per particle pair. The inheritance is the converter's and `convert_spec_equivalence` compares the potential of every particle pair; the script checks the record of pair 1, no default, and the other parts as before. |
| `pyazr/composite_terms` (2b) | the component layouts of an old segment line (two numbers per component of old GUI files, too few numbers) | old grammar only: the part edits the old original of `sum_data.azr` and runs it (blocks the refusal switch); parts 1 and 2 run on records (the ratio's denominator reaction, a scaled term `s*…`). |
| `cli_route_exit` (badscale) | a SUM component scaling that is not a number: refused by the set-up, exit 1 | in format 3 the scale of a term is refused when the file is read (`E_REACTION_SYNTAX`), where every refused project exits 255 (`AZURE2.cpp`, with the misleading line "Could not open … Check that file exists."; listed as a defect): the check asks for a non-zero exit. |
| `T2` truth runs of `make_t2.py`, `model_data` (fit project) | a THM experiment whose data sets are all inactive (the old engine ran it) | refused in format 3 ("'segments=' is not key=value": `ConfigThm.cpp` writes `segments=` with no key); the generator and the check leave the experiment line out of such runs. Listed as a defect. |
| `kinematic_limits` (c) | pair 1 with its particles swapped | the labels become `1H(12C,12C)1H`, `1H(12C,g)13N`: `p` is not accepted as a target nucleus in a label. |
| `total_optical` (all of it) | a natural sum of totals: "natural 1 total 1" on a line with components, its gradient check on that line, and "total 1" on a differential line (refused) | **format 3 since `feature/upgrade-4-format3-gaps`** (director's option 1): `sigma_total` takes a sum of elastic terms, one per nucleus, on a target with isotopes (`BuildEffects.cpp` `NaturalTotals` refuses a sum without isotopes, two terms on one pair, a factor on a term, a term whose exit pair is not its entrance); the converter converts the old line (`convert_constructs`). The check builds its models from `tests/natural_neon` (records); its outputs equal those of the old models on the base build file by file but for param.par names and the sample name in the `AZUREOut_sample_seg` headers; the analytic gradient of (b) agrees with the old run's to 5.6e-15 of the largest component (the finite differences and the Jacobian path to every digit; the conversion of the old model gives the same 5.5e-15, a property of the reading path). The refusals are on records (charged pair, inelastic label, an angle column, and the three set-up refusals). |
| `natural_target/models/M_thick.azr` | one `<targetInt>` line (a layer, nothing fitted) over x1, x2 (no sample) and x3, x4 (the natural sample): the converter's four questions "its target is shared with x1/x2, which has not the same sample" | **format 3 since `feature/upgrade-4-format3-gaps`**. What the old engine does with the line: it integrates every listed line over the one target and weights only the sample's lines by the abundances (`NaturalAllowed`, the sample token of the line). So the four questions have one answer: the sample's lines get a copy of the target with the isotopes, the others the target without (`BuildApparatus.cpp`, `BuildSamples`; a target with a fitted density, Doppler temperature or screening stays a question, since the copies would be two parameters). Converted in place, the original to the fixtures; outputs equal to the old file's but param.par names and the sample name in the `AZUREOut_sample_seg` headers (`convert_constructs` case; `predictions.txt` lines relabelled). |
| converter: a transmission or yield line under a Doppler-only `<targetInt>` line | the question "a sample thickness and a `<targetInt>` target: two areal densities" | since the gaps branch the converter writes one sample target with the `doppler` field and a note (`MoveDopplerToSamples`; outputs byte-identical to the old line's on the old `transmission_yield` grad_fit model but param.par); a fitted temperature shared by two lines stays a question (it would become two parameters). |
| `natural_target` (d) `bad_key`, `bad_sample` | a component naming a pair key `<levels>` lacks; a sample `<abundances>` lacks | the old engine refused both; the converter dropped the line (`bad_key`) or kept it with a note (`bad_sample`), so the conversions ran. The converter now refuses an active line of either kind with the old engine's reason (`convert_constructs`); the old files are fixtures only and the check converts them and expects the refusal. |
| `doppler_fgm` (g) | one fitted temperature (`doppler_1`) on a `<targetInt>` line shared by plain cross-section lines and transmission/yield lines | a cross section takes a target of kind none or layer, a transmission a target of kind sample (`E_TARGET_KIND`), so the plain and the transmission data sets have two targets and two temperatures (D0, D1; D2 as before): the gradient check covers three temperatures instead of two (16 -> 17 parameters; the gradient agreements equal the base's to the printed digits but the Jacobian path, 3.9e-14 -> 2.7e-13). |
| `doppler_fgm` (i) `vary 2`, `transmission_yield` (e) both keywords, transmission of a differential line, (l) `vary 2`, `thicknessgroup` without a label or without its own block | old-grammar forms | no format-3 form. The converter refuses each with the old engine's reason (`convert_constructs`; the converter used to convert `doppler ... 2 0` to a fixed temperature, `... 1 -5` to a fitted one, `... 1 0 200` to `grid_max 200`, both keywords to a transmission and a member without its own block to a plain line — files the old engine refused ran). The format-3 counterparts that can be written are refused by the run: T <= 0 and a sample thickness of 0 (new `E_OUT_OF_RANGE` checks in `experiment/src/BuildEffects.cpp`: format 3 ran both, the first without broadening, the second as a plain cross section), a prior without its width or with a negative one (`E_PRIOR_SYNTAX`), `grid_max` <= T (the run's message), a cross section on a sample (`E_TARGET_KIND`). |

Group 8: `beta_delayed_alpha` identical to the base check's output on the base build; `model_data` identical but for the
two lines of the old `<scatteringLengths>` total; `upgrade` identical but for the T3 noise-free recovery (d): the base's
`recover.py` wrote the edited old block with `</scatteringLengths>` on the line of its last constraint, which the old
engine then dropped (one of two constraints active); the rewrite keeps both (fit parameters move by < 0.01 sigma).

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
- pyazr tests of the old file class, on the originals (`tests/pyazr/old_tree.py`): those listed in §3 (groups 2, 4a,
  4b, 5: `roundtrip`, `physical_norm`, `shift_first_call`, `target_density`, `jacobian_skip`, `jacobian_shift_columns`,
  `zero_width_level`, `formal_capture_mask`, `composite_gradient`, `norm_group`, `prior_covariance`, `save_fit_norms`,
  `extrap_layout`, `objective`, `shift_fn`, `sqrt_shift`, `data_covariance`, `angular_coefficients`, `screening`,
  `parameter_kinds`; group 6: `natural_target`, `doppler`, `thickness`, `parameter_kinds` (the old generators)). `docs/upgrade/pyazr_on_records.md` lists their replacements as "gone with formats 1-2" (old
  segment fields, `SegmentSet`) or "to be built" (level-row edits); none was moved to the format-3 session for that
  reason. `pyazr_inactive_lines` (part 2), `pyazr_shared_state` and `pyazr_rmatrix_pole` edit format-3 records. Each loads the
  original through `pyazr/azrfile.py` and runs the engine on it: they block the switch until pyazr's session reads
  format 3 and the old class is deleted (another agent's work).
- Group 8 (this page, §3, §4): pyazr tests of the old reader or the converter on the originals (`scattering_length`
  block part, `azr_parse_cache`, `potential_record`, `save_fit_format3` part 2, `validate`, `convert_azr`,
  `session_records` old T2, `composite_terms` 2b); `beta_delayed_alpha` (e) converts its old line with `AZURE2 --convert` (on purpose, as `segline_formats`); the generator of
  `light_excitation` (run by its check; models written at run time in the old grammar).
- The model generators. Run by a check and writing format 3: `clone_composite/make_models.py`,
  `light_excitation/make_models.py` (its models as records). Writing format 3 through the converter (they compose the
  old grammar and call `AZURE2 --convert --to 3 --exact-grids`; they give the committed files byte for byte): the
  generators of T1–T3 (`tests/upgrade/lib/upgrade_project.py`, `T1/generate.py`, `T3/generate.py`, `T2/make_t2.py`).
  Rewritten with group 6 to write format 3: `doppler_fgm`, `transmission_yield`; `natural_target/make_models.py` is deleted.
  **Deleted** (branch `feature/upgrade-4-format3-gaps`; they wrote old files, their committed output, converted, is
  what the tests use): the generators of `angular_coefficients`, `capture_ay_coupling`, `convolution_unify`,
  `convolution_window`, `data_covariance`, `ec_e1m1_phase`, `ec_m1_exact`, `ec_threshold`, `electron_screening`,
  `gamma_lab_frame`, `hard_sphere`, `mcmc_physical_norm` (`make_model.py`), `q_attenuation`, `ratio_frame`,
  `relativistic`, `scattering_length`, `screening_6li`, `secondary_gamma`, `shift_sqrt_subpoints`,
  `shift_subpoint_angle`, `sum_lab_frame`, `target_density`, `transform_level_order` and `reaction_rate`
  (`make_truth.sh` reads the committed `*_truth.azr`); each README has a section "How the models were made" with the
  generator's description and the commit that holds it.  Left: `18O_p_a_thm/make_inputs.py` (group 7, another
  branch).
- The old generators `tests/convert/fixtures/tests/doppler_fgm/make_models.py`, `.../transmission_yield/make_models.py`
  (on the original n12C), for the pyazr tests `doppler`, `thickness`, `parameter_kinds` of the old file class.

## 6. User-guide pages that show old keyword lines (for a later documentation pass)

`docs/source/user_guide/segments.rst`, `experimental_effects.rst`, `chi_squared.rst`, `pyazr.rst`;
`docs/source/reference/data_formats.rst`, `command_line.rst`; `docs/source/examples/natural_neon.rst`;
`docs/source/theory/thm_implementation.rst`, `polarization_implementation.rst`; `docs/source/developer/
architecture.rst`, `coulomb_table.rst`, `edata_and_target_effects.rst`; `docs/source/getting_started/whats_new.rst`.
No test runs an `.azr` snippet of `docs/source/`.
