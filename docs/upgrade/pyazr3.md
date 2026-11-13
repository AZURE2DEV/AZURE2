# pyazr 3.0: where the binding stands, and the design of the release

This note lays the Python binding `pyazr` (3.0.0.dev0) against the target of `ARCHITECTURE.md` §8 and §11 and
`ROADMAP.md` Phase 8, from measurements on the tree at commit 441a0907 (`feature/upgrade-4-data`), and gives the
design of pyazr 3.0: its objects, one rule for arrays, one for errors and self-description, the fit, the physics
explorer, the command line, the packaging, the order of work and the decisions that are the PI's. The PI's measures:
pyazr must be intuitive in Python and in the GUI, usable by an AI agent without the manual, easy to maintain, and a new
kind of data or physics must be cheap to add; no backward compatibility is owed (a converter of old files only). The
web UI of Phase 9 will be a thin layer over pyazr objects, so what pyazr exposes is what the UI can show.

Every number of section 1 comes from a script of `work/pyazr3/` (outside the repository, listed in section 1.7 with
how to run it), run on this commit with the build `work/b-gecm` of the same commit, two threads.

## 1. Where pyazr stands (measured)

### 1.1 Public names

148 public names: the 99 of `pyazr.__all__` (97 functions and classes, two modules: `masses`, `stopping`) and 49 more
that submodules make public (their `__all__`, or names without a leading underscore) and the package does not
re-export.

| measure | `__all__` | submodules | all |
|---|---|---|---|
| names | 99 | 49 | 148 |
| reach the engine (`_azure2`), static call graph | 30 | 12 | 42 |
| compute in Python only (records, readers, formulas) | 67 | 31 | 98 |
| modules and constants | 2 | 6 | 8 |
| named by a test of `tests/pyazr/` (imported or `pyazr.<name>`) | 36 | 22 | 58 |
| entered by a test (function-level trace of 115 of the 120 tests), of the callables | 53 of 88 | 35 of 41 | 88 of 129 |
| named in `docs/source` | 49 | 17 | 66 |

The session class `azure2` has 85 public members; 18 are entered by no test: `calculate_analyzing_power`,
`calculate_analyzing_power_rwa`, `calculate_angular_dists`, `calculate_sfactor`, `dimensionless_widths`,
`energy_indices`, `level_free_width_indices`, `level_scheme`, `mass_number`, `nuclear_potentials`, `only_level`,
`without_level`, `recalculate_external_capture`, `refresh_pairs`, `refresh_parameters`, `transform_physical`,
`wigner_widths`, `SHIFT_PENALTY_FROM_ZERO`. Of `ParameterSet`'s 24 members (the list's own excluded) 13 are entered
by no test (`by_level`, `by_name`, `put`, `take`, `indices`, `table`, `to_records`, `values`, `fixed`, `shifts`,
`dopplers`, `thicknesses`, `wigner_limits`). Of `pyazr.bands` only `best_fit_params` and `live_parameters` are entered
(`sensitivities`, `uncertainty_bands`, `extrapolation_bands`, `trimmed_model` are not), and none of
`pyazr.transform`'s formulas is.

Every name, by module. Flags: E reaches the engine, P Python only, - a module or a constant; t named in a test; x
entered by a test (- for a constant or a record class without methods, where the trace cannot tell); d named in
`docs/source`. (`Project.add_<kind>` are generated closures: the trace does not attribute them; seven tests call them,
`project_test` among them.)

- `azure2`: azure2 Etxd, ScatteringLength P.x.
- `convert`: convert_azr Etxd, convert_params Etxd
- `physics`: whittaker_wave Etxd
- `tabulate`: tabulate Etx., TabulatedCrossSection P.x.
- `angular`: angular_distribution Etxd
- `parameters`: Pair Ptxd, PairSet P.x., Parameter Etx., ParameterSet E.x., LevelKey P.xd, NuclearPotential P.x.
- `scheme`: LevelScheme P..., SchemeLevel P..., SchemeChannel P...
- `masses`: masses -..d
- `stopping`: stopping -t.d
- `outputs`: fitted_norms P...
- `nds`: ExforData P.x., ExforDataset P.-., GroundState P.-., Gamma P.-d, Reference P.-., search_exfor P..d,
  fetch_exfor P.xd, fetch_levels P.xd, fetch_gammas P..., fetch_ground_state P..., reference P.xd, resolve_doi P..d,
  aps_doi P..., sfactor_to_cross_section Etxd, cm_to_lab P...
- `widths`: ChannelWidth P..., WidthTable P..., weisskopf_width P..., teichmann_wigner Ptx., lane_thomas_wigner Ptxd
- `transform`: Channel P..., Level P.x., TransformedChannel P..., TransformedLevel P..., transform_out P...,
  partial_widths P..., LevelWidthTransformer P..., levels_from_azr Etx., levels_from_scheme P..., penetrability P..d,
  shift P..d, shift_derivative P..d, whittaker P..d
- `nuclei`: Nucleus P..., MergeReport P.-., merge_projects Etxd, nucleus_of P..., nucleus_of_pair P...
- `bands`: Band P..., load_covariance P..d, rmatrix_columns E..., best_fit_params Etxd, live_parameters E.x.,
  step_sizes P..., sensitivities E..d, trimmed_model E..., uncertainty_bands E..d, extrapolation_bands E..d
- `modelavg`: Variant Ptxd, Averaged P.x., ModelAverage P.x., model_average Ptxd, parameter_label Ptxd,
  azr_parameter_slots Ptx., write_averaged_azr Etxd
- `priors`: prior_from_fit Ptxd, read_prior Ptxd, write_prior Ptxd, read_covariance_dat Ptxd
- `schema`: schema Etxd, validate Etxd, describe Etxd, predictions Etxd
- `model_data`: ModelDatum Pt-., model_data Etxd
- `grids`: extrapolation_grid Etxd, refine_grid Etxd
- `project`: Project Etxd, Record E.x., Nuisance Et-d, FormatError Ptxd, records_from_csv Etxd
- `nuclear_part`: ProjectLevel P.xd, ProjectChannel P.-., ProjectPair P.-d
- `described`: DescribedDataset P.x., DescribedDatasets P.xd, DescribedExtrapolation P.x., DescribedExtrapolations
  P.x.
- `bands` (not in `__all__`): active_indices P...
- `exfor_coulomb` (not in `__all__`): quantity_code P.x., is_elastic P.x., quantity P.xd, chosen_quantity P.x.,
  reaction_label P.x., record P.xd
- `exfor_sfactor` (not in `__all__`): reaction_label P.x., dataset P.xd
- `masses` (not in `__all__`): AME_VERSION -t-., MassEntry P.x., PairMasses P.-., lookup Ptx., parse_nuclide Ptx.,
  atomic_mass P..d, mass_excess P..., separation_energy Ptxd, pair_masses Ptxd, identify P.x., U_MEV -t-.,
  PairMassCheck P.x., mass_check Ptxd, fill_masses Ptxd
- `modelavg` (not in `__all__`): level_label Ptx., width_label Ptx., WEIGHT_METHODS -.-.
- `nds` (not in `__all__`): ExforColumn P.-., exfor_master_path P.xd, elastic_lab_to_cm P.x., Level P.x.
- `nuclear_part` (not in `__all__`): NuclearPart E.x.
- `nuclei` (not in `__all__`): SYMBOLS -t-.
- `parameters` (not in `__all__`): parameter_kinds Etxd
- `snapshot` (not in `__all__`): write_format3 E.x., place_level_values Ptxd, write_cbkg E.x.
- `stopping` (not in `__all__`): CITATION -t-d, IMPLEMENTATION -.-., Material Etxd, projectile Etx., energy_loss Etx.,
  energy_after Etxd, areal_density_between Etx., range_of Etxd, straggling Etxd, provenance Etxd
- `transform` (not in `__all__`): sommerfeld_eta P..., coulomb_rho P..., coulomb_waves P...

The full table, with the first line of each docstring, is `work/pyazr3/names.tsv` and `exercised.tsv`.

What `help()` gives an agent today: `pyazr.__doc__` is `None`; `help(pyazr.azure2)` shows the docstring of the mixin
`_HybridPotentials` ("Mixed into :class:`pyazr.azure2`: uses its ``sess``, ``pairs`` and ``configure``"), because
the session class has none of its own; the constructor's docstring describes `use_rmc` as "the
R-matrix-with-channels (RMC) formalism" where the program's `--use-rmc` is the Reich-Moore approximation; the free
vector `params_rwa`, named by the docstrings of `chi2_and_grad` and `save_fit`, is an instance attribute that no
`help()` lists; `ParameterSet.densities`, `dopplers`, `thicknesses` and `save_fit` still speak of `<targetInt>`,
`<segmentsData>` lines and `<parameterSettings>` `prior_centre` rows of the old formats.

### 1.2 The binding `_azure2` and `AZUREAPI`

155 bound names in eight files: `pybind_api.cpp` 89 (the `Session` class 71, `RuntimeOptions` 14 fields, four
functions), `pybind_project.cpp` 18, `pybind_nuclear_part.cpp` 17, `pybind_physics.cpp` 9, `pybind_stopping.cpp` 8,
`pybind_convert.cpp` 6, `pybind_schema.cpp` 5, `pybind_thm.cpp` 3; `pybind_model_data.cpp` and
`pybind_potentials.cpp` hold helpers.

- Used by pyazr: 135 (the 14 option fields through `setattr` in `azure2.__init__`). Not used by pyazr: 20. Of these,
  11 are used by tests only (`project_format`, `Session.norms`, `calculated_segments_e1/e2`, `statistical_factor`,
  `resonance_strength`, `branching_ratio`, `rutherford_cross_section`, `mott_cross_section`, `schema_json`,
  `thm_nuclides`) and 9 by nothing (`Session.configured_potential_pairs`, `update_norms`, `update_parameters`,
  `params_all`, `norms_errors`, `transform_all_rwa_full`, `modified_sfactor_conversion`, `Document.common_ending`,
  `thm_global_opticals`).
- Copies or views: every array is a copy. Of the `Session` methods 39 return a `py::array_t` built from a
  `std::vector` (`to_array`), 7 lists or tuples, 4 dicts, 15 scalars, 6 nothing. No method returns a view.
- State: the `Session` object holds the engine (`AZUREAPI`, its `Config` and run scratch) and nothing process-wide;
  the evaluations release the GIL, the constructor keeps it. `Document` holds a project's text. The free functions
  (physics, stopping, schema, convert, nuclear-part edits) hold no state.
- `AZUREAPI` (`api/include/AZUREAPI.h`): 72 public methods; 69 are called by the binding; 3 are not
  (`RadiusOverrides`, `configure`, `results`). Two of them hand over positional numeric layouts kept "in sync" by
  hand with pyazr: `GetParameterInfo` (16 doubles per parameter, `Parameter._NFIELDS`) and `GetPairsInfo` (21 per
  pair, `Pair._NFIELDS`); a third, `parameter_records()`, already gives each parameter as a record with its kind and
  owner fields.  (Since the second step of section 4 both layouts are deleted: `parameter_records()` and
  `pair_records()` are the only path, their fields one table of the binding, `api/src/pybind_records.cpp`.)

### 1.3 Physics computed in Python

17 formulas; 2 have a test against the engine; the penetrability is compared here (use case (c)).

| formula (pyazr) | replaced by (`pyazr.physics`, library) | largest rel. difference (`physics_table_test`) |
|---|---|---|
| `transform.sommerfeld_eta`, `coulomb_rho` | `sommerfeld_eta`, `wave_number` (`ChannelFunctions.h`) | 0; 2.2e-16 |
| `transform.coulomb_waves` (mpmath) | `coulomb_waves` (`physics/coulomb`) | 5.2e-10 for η > 0 (η = 0: see 1.6) |
| `transform.penetrability` | `penetrability` (`ChannelFunctions.h`) | 9.7e-10 (four pairs, l = 0–3) |
| `transform.shift` | `shift` (same) | 2.6e-8 (3e-9 absolute; Python right) |
| `transform.shift_derivative` (central difference) | `shift_derivative` (same) | 9.5e-6 (Python right, 1.6) |
| `transform.whittaker` (scipy, mpmath) | `whittaker_wave` (`Whittaker.h`) | 3.6e-10 |
| photon penetrability (E_γ/ħc)^(2L+1) (`transform`) | `photon_penetrability` (`Widths.h`) | 0 |
| `transform_out`, `partial_widths`, `LevelWidthTransformer` | the session's `transform_rwa` | not a row |
| `transform.ancs` | the session's ANCs (`parameters.out`) | not a row |
| `widths.weisskopf_width` | none: the engine has no Weisskopf estimate | not bound |
| `widths.teichmann_wigner`, `lane_thomas_wigner` | `wigner_limit` (`Widths.h`, `AChannel`) | 0 |
| `nds.cm_to_lab` | `cm_to_lab_energy` (`experiment/Frames.h`) | 2.2e-16 |
| `nds.elastic_lab_to_cm` (angle and Jacobian) | `elastic_lab_to_cm_angle` (`Frames.h`) | 2.5e-16; J 1.6e-14 |
| `bands.sensitivities`, `uncertainty_bands` | `CovarianceBand` (not bound) | not a row |
| the rate integral of `examples/reaction_rate.py` | `ReactionRate.cpp` (not bound) | not a row |

`azure2.wigner_widths` obtains the penetrability at the level energy by setting every reduced width to a small value
and reading the engine's widths back (a workaround for a function the binding did not have).
`nds.sfactor_to_cross_section` calls the engine's conversion (tested in `units_scales_test`). Not counted as
physics: `masses` (reads the AME2020 table, which the engine does not have), `modelavg` (statistics), `priors` (file
formats), `tabulate` (interpolation of engine values).

### 1.4 What the engine can do that pyazr cannot reach

**Command-line options** (35 in `AZURE2 --help`): 4 concern the program itself (`--help`, `--no-gui`, `--no-readline`,
`--use-api`); 17 are session options (`transform`, `use_long_wavelength`, `use_brune`, `use_park`,
`ignore_externals`, `use_rmc`, `use_gsl_coul`, `use_coulomb_table`, `coulomb_table`, `use_ec_table`,
`ec_convection_current`, `replace_grids`, with their negations); 6 are functions (`convert_azr`, `convert_params`,
`validate`, `describe`, `predictions`, `schema`); `--refreeze-covariance` is reached in part (one
`refreeze_covariance()`, the refits are the caller's); 7 are not reachable: `--build-coulomb-table`,
`--use-gradient`, `--max-calls`, `--use-lm`, `--use-gsl-lm`, `--covariance-band`, `--scale-covariance` (and
`--use-nlopt`, which the help does not list).

**Menu modes**: calculate with data and extrapolate are reachable (a session, `extrap_mode`); fit (the engine's
MIGRAD, LM, GSL-LM, NLopt), MINOS, the reaction rate and the MCMC are not.

**Record kinds, data kinds, parameter kinds**: all reachable. The 10 record kinds through `Project.add_<kind>`
(generated from the schema), the 23 quantities by any data set or datum a project states (the session evaluates
them; use case (b) added one), the 14 parameter kinds through `parameters` (`Parameter.kind`); 12 of the 14 have a
view of their own in `ParameterSet` (`radiative` and `apparatus` have none).

**Output files**: on 13N the program's menu 1 writes 9 files, pyazr's `write_output_files` 3 of them
(`AZUREOut_aa=1_R=1.out`, `R=2`, `chiSquared.out`, byte for byte the program's); on T2 it writes also
`normalizations.out`, `shifts.out`, `thm_experiments.out`, `covariance.out`. Never written by pyazr: `modelData.out`,
`param.par`, `parameters.out`, `scatteringLengths.out`, `relativistic.out`, `fitStatus.out`, `covariance.dat` and the
band, `reactionrates.out`, `cross_section_grid.dat`, the MCMC chain. Some of their numbers are reachable as values
(`model_data`, `scattering_lengths`, `physical_levels`, `transform_rwa`); `param.sav` is written by `save_fit` or by
`update_sav_from_rwa_params` (which needs an existing `param.sav`).

**Library functions** (the explorer of `ARCHITECTURE.md` §8): 218 functions are declared in the public headers of
`physics/`, `physics/coulomb`, `physics/transfer` and `stopping/`; the binding calls 29 of them (stopping 17 of 55,
`physics/` 10 of 78, Coulomb 1 of 7, transfer 1 of 78). Against the list of §8:

- Coulomb F, G, F′, G′: through a session only (`azure2.coulomb_functions`, at a pair's radius, on energies, not
  on ρ); the methods of pycoulomb side by side: no (only `use_gsl_coul` per session).
- P_l, S_l, φ_l: through a session (`coulomb_functions`: `P`, `S`, `delta_hs`); σ_l: no (the library has
  `coulomb_phase` in `Coulomb.h`, not bound; use case (c) used scipy).
- Whittaker functions: yes (`whittaker_wave`, bound; and `transform.whittaker` in Python).
- Wigner limits: per channel (`Parameter.wigner_limit`); single-particle widths: no.
- Level-matrix poles and eigenphases, R and U matrix elements: no. Phase shifts: only as a quantity of a project.
- THM vertex W_l(p), momentum distributions, distortion, line shape: per experiment through a session
  (`thm_vertex`, `thm_distortion`, `thm_lineshape`); the global optical sets: bound, no pyazr name.
- DWBA bound states and distorted waves: no.
- Stopping and straggling versus E and depth: yes (`pyazr.stopping`).
- Kernel rows K_i(E_n) and the grid E_n of a data set: no.
- Resonance quantities, Rutherford and Mott cross sections, S-factor conversions: bound, only the S-factor one has a
  pyazr name (`nds.sfactor_to_cross_section`).
- External-capture integrals: yes (`azure2.ec_integrals`).

### 1.5 Four use cases, as a user and as an agent

Each was written from the docstrings and `help()`, with every pyazr call counted and every wrong call tried once.

**(a) T2: load, χ², change a width, fit with scipy, write the parameters** (`uc_a_t2.py`, log `uc_a.log`). 23
calls. iminuit and emcee are not in the environment; scipy's L-BFGS-B on `objective_and_grad`. The session's χ²
equals the program's on the same project: `calculate_chi2_rwa` 10905.2313 against `Total-Chi-Squared: 10905.2`, the
objective 10910.1338 = data + `Total-Prior-Covariance-Chi-Squared` 4.9025, and `chiSquared.out` written by the
session is byte for byte the program's. Changing a width: in the vector (`ParameterSet([w]).put`) and in the file
(`Project.set_channel` then `session()`). The fit: 25 iterations, 31 evaluations, 23 s at two threads, objective
10910.1 → 689.26. Writing: `save_fit` refused the snapshot (below); with `verify=False` it wrote `fit.azr` and
`fit.sav`, and the reopened snapshot's objective is 689.2479 against the fit's 689.2579.
Looked up outside the docstrings: that the session's docstring is on `__init__`; that the start vector is the
attribute `params_rwa` (found by `dir()`); how to get a `ParameterSet` of one parameter (a slice is a plain `list`:
`AttributeError: 'list' object has no attribute 'put'`); which output file holds which term.
Messages that do not say what to do next: `save_fit`'s "the snapshot run on its own gives chi2 = 683.7333853, the
fit 683.7433904; it has been removed." (no cause, no `verify=False`); `objective` with a NaN: "Coulomb functions
failed at l = 0, eta = nan, rho = nan (physics/coulomb, coulsm)" (names no parameter); `parameters.by_name` of an
unknown name: `KeyError: 'nonexistent'`; `update_sav_from_rwa_params` without a `param.sav`: a bare
`FileNotFoundError`; `set_channel` on a channel the level lacks: right at its end, but 870 characters long because it
starts with the whole `ProjectLevel` repr. Good messages: a vector of the wrong length ("expected the free parameter
vector (19 values, params_rwa / params order), got 18"), an unknown level label (lists the levels), a record without a
required field (the code, the field, its grammar).
Late or silent: `calculate_chi2_rwa` returns a one-element list; it adds the abundance prior to the data χ² (0.010036
at the fitted point) although its docstring says "Total data chi-squared" and `segment_chi2`'s says their sum equals
it exactly: the sum is 683.7333544, `calculate_chi2_rwa` 683.7433904 (`uc_a_state.py`). This difference is what
`save_fit` compares, and the cause of its refusal (section 1.6).

**(b) A resonance strength added to a copy of 13N** (`uc_b_datum.py`, `uc_b.log`). 10 calls. `add_datum` with the
help of `help(Project.add_datum)` (generated from the engine's table: fields, units, which are required); the
engine's warning `W_BROAD_LEVEL` came at once. The model value 0.4799932236 eV is the program's (`modelData.out`),
and `chiSquared.out` and both `AZUREOut` files are byte for byte the program's. Looked up: that the datum's "level ID
of <levels> (column 9)" is `ProjectLevel.ids`; that the program needs the `checks/` directory the project names
(a session does not). `write_output_files` does not write `modelData.out` (nor `param.par`, `parameters.out`,
`scatteringLengths.out`). Every wrong call (J that does not match, unknown level, unit of the wrong dimension, no
unit, a misspelt field) was refused at once with the engine's code and a "did you mean".

**(c) A penetrability and a Coulomb phase for p + ¹²C of 13N** (`uc_c_channel.py`, `uc_c.log`, table and plot in
`uc_c/`). 5 calls. `coulomb_functions(1, E, L)` gives F, G, P, S and the hard-sphere phase on an energy grid;
pyazr's Python `penetrability` agrees to 3.8e-10 (l = 0), 2.7e-10 (1), 2.6e-11 (2) at 20 ms per point (scalar
only). The Coulomb phase σ_l is in no pyazr name: computed with `scipy.special.loggamma`. Looked up: the units of
`penetrability`'s `red_mass` (u, from the source). Silent: L = −1, a negative energy and the photon pair all return
zeros; pair 0 or 9 raise `IndexError: index 0 is out of bounds for axis 0 with size 0`. Printed to standard output
by every session whose output directory does not exist: "Could not write to EC Amplitude File." (no file named, no
consequence said; the run is right).

**(d) Two sessions in two directories on two threads** (`uc_d_threads.py`, `uc_d.log`). 11 calls per pair of
sessions. Two copies of 13N, the second with pair 1's radius at 4.2 fm, each thread building its session, evaluating
the objective and gradient and writing its files; the bytes of the objective and gradient and the three output files
equal those of each session alone, and the two differ from each other. Wall time 24.2 s on two threads against 23.1 s
one after the other: no gain, since the set-up (6 s per session) holds the GIL and dominates. The thread count of a
session can only be set through `OMP_NUM_THREADS` for the process.

### 1.6 Wrong numbers and defects found

pyazr against the program on the same project and the same parameters: **no wrong number**. On T2 at the fitted
point, the program on the snapshot and a fresh session on it write identical files (`uc_a_state.py`); the session's
`chiSquared.out` at x0 is the program's byte for byte (use case (a)); on 13N with a datum all common files are
identical (use case (b)).

Defects found (not fixed here):

1. **The engine drops the prior of a target `isotope` whose start value differs from its prior centre** (program and
   pyazr alike). T2 with `isotope 15N 0.0037 fit prior 0.0036 0.0002` (prior term 0.25 expected at the start) gives
   no `Total-Abundance-Chi-Squared` and the same `Total-Chi-Squared` as with `prior 0.0037`; with `isotope 15N
   0.00368 fit prior 0.0037 0.0002` (0.010 expected) none either; `Parameter.prior_nominal` is `None` for it. The
   term appears only when the parameter moves away from its start value (0.010036 at the fitted point). The prior's
   centre is the start value, not the field's `<centre>`; the norms of the same file keep their centre
   (`work/pyazr3/iso_check/`).
2. Through 1, `save_fit` on a fit with a free abundance refuses its own snapshot: the snapshot's start value is the
   fitted abundance, so its prior term vanishes; the message names a χ² difference, not the cause.
3. `calculate_chi2_rwa` includes the abundance prior; its docstring and `segment_chi2`'s say otherwise (1.5 (a)).
4. `write_output_files` writes only a part of the program's output files (3 of 9 on 13N; 1.4).
5. The documentation of `use_rmc` contradicts the program's (1.1).
6. Two tests of `tests/pyazr/` run the engine at four threads whatever the caller's budget:
   `parallel_points_test` (`OMP_NUM_THREADS=4` in its environment) and `upgrade_thread_independence_test`.

Found when the physics table was compared with mpmath (2026-11-13, `tests/pyazr/physics_table_test.py`), all in the
library, none changed: the engine's shift function S is good to about 3e-9 absolute (2.6e-8 relative at S = -0.06,
α + 12C, l = 1, 20 MeV, a = 5.5 fm), from the 5e-10 precision of coulsm's F and G; the engine's dS/dE (Richardson's
difference of S with h = 1e-3 E) carries that noise to 9.5e-6 relative (p + 12C, l = 0, 3.42 MeV, a = 5 fm), against
the 1e-8 that `ChannelFunctions.h` states; the Python formula (mpmath, a five-point difference with h = 1e-3 MeV) is
right there to 1e-11.  coulsm and coul return a wrong F for η = 0, l ≥ 2 at small ρ (factor 4 at l = 6, ρ = 0.2),
as `physics/coulomb/README.md` records; P, S and φ take F only through F/G and are not affected.

### 1.7 The scripts (in `work/pyazr3/`, outside the repository)

Environment: `env.sh` (conda environment, `PYAZR_ENGINE_DIR=work/b-gecm/pyazr-out`, the worktree on `PYTHONPATH`,
two threads). Each inventory regenerates with the tree as its argument.

- `inv_names.py <tree> > names.tsv`: public names, engine or Python (static call graph), tests naming them, docs.
- `trace/run.sh`, `trace/site/sitecustomize.py`: the tests of `tests/pyazr/` under a function-level tracer (115 of
  120; not run: `thickness`, `concurrent_sessions`, `doppler`, `session_isolation` for time,
  `upgrade_thread_independence` stopped at four threads); `inv_trace.py names.tsv trace/out > exercised.tsv`.
- `inv_bindings.py <tree> bindings.tsv api.tsv`: bound names, return kinds, GIL, users in pyazr; `AZUREAPI` methods.
- `inv_physics.py <tree> > physics.tsv`: functions of the library headers and the binding files that call them.
- Use cases: `uc_a_t2.py`, `uc_a_history.py`, `uc_a_cli_at_x.py`, `uc_a_state.py`, `uc_b_datum.py`,
  `uc_c_channel.py`, `uc_d_threads.py` (logs `*.log` beside them). Quick enough to become tests: (b), (c), (d)
  (6–25 s); (a) takes 40 s.

## 2. The design of pyazr 3.0

### 2.1 Principles

1. **One object per thing the core has**, named as the core names it (`ARCHITECTURE.md` §11): the project, the
   session (one built model), the parameters (the fit vector and what each entry is), the objective, an evaluation
   (the numbers at one parameter point), the physics functions. Everything else is a function or a plain record.
2. **The engine's tables are the API's tables.** Record kinds, data kinds, parameter kinds, options, physics
   functions, diagnostics: pyazr reads each from the engine at import (as `Project.add_<kind>` is generated today
   from the schema) and keeps no list of its own. A new row in the engine is a new name in Python with no Python edit.
3. **One parameter vector**: the fit vector `x` (the engine's free parameters, in its order and basis). Physical
   values are a transformation of `x` asked for explicitly, never a second family of methods.
4. **State owned by the caller**: a session is an object, two sessions share nothing, the only process-wide objects
   are immutable tables (`ARCHITECTURE.md` §9). The thread count is an argument of the session.
5. **No physics in Python** but library parsing and presentation; any Python formula left has a test against the
   engine.
6. **Self-describing**: every object answers `describe()` with plain, JSON-serialisable data; every error says what
   was wrong, what is allowed and which call shows or fixes it.

### 2.2 Modules and objects

- `pyazr.project`: `Project`, `Record`, `Nuisance`, `ProjectLevel`, `ProjectChannel`, `ProjectPair` (as today).
- `pyazr.session`: `Session` (replaces the class `azure2`, `azure2.py`, 2316 lines).
- `pyazr.parameters`: `Parameters` (the table of the fit vector), `Parameter` (replace `ParameterSet`, `LevelKey`,
  the `*_indices` methods).
- `pyazr.objective`: `Objective` (replaces `objective`, `objective_and_grad`, `chi2_and_grad`,
  `residual_jacobian`, `residuals`, `penalties`).
- `pyazr.evaluate`: `Evaluation` (model values, residuals and χ² terms at one `x`), `DataSet` (the set-up of one data
  set, read-only views) (replace the 14 `calculate_*` methods, `segment_chi2`, `dataset_chi2`, `datasets`).
- `pyazr.system`: `Pairs`, `Pair`, `Levels`, `Level` (the built system: masses, radii, channels, physical levels)
  (replace `Pair`, `PairSet`, `LevelScheme`, `Scheme*`, `physical_levels`, `nuclei`).
- `pyazr.physics`: free functions on arrays (section 2.7) (replace `transform`'s and `widths`' formulas).
- `pyazr.stopping`: as today.
- `pyazr.sources`: `exfor`, `ensdf`, `masses` (library readers; today `nds`, `exfor_coulomb`, `exfor_sfactor`,
  `masses`).
- `pyazr.bands`: `uncertainty_bands`, `sensitivities`, `extrapolation_bands`, `model_average` (today `bands`,
  `modelavg`, `tabulate`).
- `pyazr.plot`: plots of an `Evaluation` and of explorer arrays (optional, matplotlib).
- `pyazr.cli`: the command `pyazr` (section 2.8).
- top level: `open(path)` (a `Project`), `describe()`, `kinds()`, `schema()`, the exceptions.

**Project**: `open`, `new`, `save`, `records`, `add`/`add_<kind>`, `set`, `remove`, `rename`, `levels`, `set_level`,
`set_channel`, `set_pair`, `config`, `set_config`, `validate`, `describe`, `session(**options)`.
```python
p = pyazr.open("13N.azr")
p.set_channel("3/2- 3.5032", pair=1, l=2, s=1, width=60e3)          # eV
p.add_datum("wg", quantity="resonance_strength", level=5, J="3/2-", reaction="12C(p,g)13N",
            value=(0.25, "eV"), error=(0.02, "eV"))
p.save("13N_edit.azr"); s = p.session(threads=2)
```

**Session** (one built model): `parameters`, `objective()`, `evaluate(x)`, `data` (`DataSet`s), `extrapolations`,
`system`, `transform(x)` (physical values), `channel_functions(pair, e_cm, l)`, `level_matrix(e_cm, x)`,
`write_outputs(x, directory=None)` (every file the program writes for a calculation), `save_fit(path, x)`,
`describe()`, `close()`, a context manager.
```python
with pyazr.Session("T2.azr", threads=2) as s:
    x = s.parameters.x0                       # a copy of the start vector
    ev = s.evaluate(x)                        # an Evaluation
    print(ev.chi2, ev.terms)                  # the data chi^2 and every named term, as chiSquared.out
    s.save_fit("fit.azr", x)
```

**Parameters** (one row per parameter; the columns are the engine's parameter record): `x0`, `names`, `kinds`,
`units`, `free`, `bounds`, `of_kind(kind)`, `of_level(label)`, `by_name(name)`, `put(x, rows, values)`,
`take(x, rows)`, `read(path)`, `write(path, x)`, `table()`, `to_records()`. A selection is again a `Parameters`.
```python
P = s.parameters
w = P.of_level("1- 12.425800").of_kind("width")
x1 = P.put(P.x0, w[0], 2 * P.take(P.x0, w[0]))
print(P.table())                              # name, kind, unit, value, free, bounds, owner
```

**Objective** (a callable for any minimiser or sampler): `__call__(x)`, `grad(x)`, `value_and_grad(x)`,
`residuals(x)`, `jacobian(x)`, `log_prob(x)`, `x0`, `names`, `bounds`, `errordef` (= 1).
```python
f = s.objective()
r = scipy.optimize.minimize(f.value_and_grad, f.x0, jac=True, bounds=f.bounds, method="L-BFGS-B")
print(f(r.x), dict(zip(f.names, r.x)))
```

**Evaluation** (the numbers at one `x`; owns copies): `chi2` (data only, a float), `terms` (name → value, the
`Total-*` of `chiSquared.out`), `objective`, `data[i].model`, `.residuals`, `.chi2`, `.energies`, `.angles`,
`.extra` (the named columns of the data kind: A_y, Legendre coefficients, E1 and E2 parts), `table(i)`.
```python
ev = s.evaluate(x)
for d in ev.data:
    print(d.name, d.kind, d.chi2, d.model[:3])
assert ev.objective == f(x)
```

**DataSet** (set-up of one data set; read-only views): `name`, `kind`, `unit`, `frame`, `energies`, `angles`,
`values`, `errors`, `grid()` (the sub-point grid E_n), `kernel()` (the rows K_i(E_n)).

### 2.3 What stays, what is renamed, what goes

What a script written for 3.0.0.dev0 has to change (no deprecation shims are owed; the converter of old project files
stays):

- `pyazr.azure2(file, **opts)` → `pyazr.Session(file, threads=, **opts)` or `Project.session()`; option names kept.
- `s.params_rwa` (214 uses in `tests/`) → `s.parameters.x0`.
- `calculate_rwa`, `calculate`, `calculate_all_rwa`, `calculate_sfactor[_rwa]`, `calculate_analyzing_power[_rwa]`,
  `calculate_angular_dists[_rwa]`, `calculate_energies`, `calculate_angles`, `calculate_excitation_energy` →
  `s.evaluate(x).data[i]` (`model`, `extra`, `energies`, ...).
- `calculate_chi2_rwa` (a one-element list, data plus abundance prior), `calculate_chi2`, `segment_chi2`,
  `dataset_chi2`, `penalties`, `scattering_length_chi2` → `Evaluation.chi2`, `.terms`, `.data[i].chi2`.
- `objective`, `objective_and_grad`, `chi2_and_grad`, `residual_jacobian`, `residuals` → `Objective`.
- `transform_rwa`, `transform_all_rwa`, `transform_physical`, `physical_levels`, `wigner_widths`,
  `dimensionless_widths` → `s.transform(x)` and `s.system.levels`.
- `norm_indices`, `shift_indices`, `energy_indices`, `level_free_width_indices`, `only_level`, `without_level`, the
  12 per-kind views of `ParameterSet` → `Parameters.of_kind`, `of_level`, `put`.
- `read_sav`, `load_sav`, `update_sav_from_rwa_params`, `update_rwa_params_from_sav`, `full_rwa_from_sav`,
  `read_sav_basis` → `Parameters.read(path)` and `Parameters.write(path, x)`.
- `write_output_files` → `write_outputs`; `extrap_mode`, `data_mode` → `s.extrapolations` beside `s.data` (no mode
  switch of a live object).
- `set_channel_radius`, `set_nuclear_potential`, `clear_nuclear_potential` on a session → edits of the `Project`
  (radius, `potential` record) and a new session (decision 11).
- `Pair`, `PairSet`, `LevelScheme`, `SchemeLevel`, `SchemeChannel` → `pyazr.system`.
- `transform.*` and `widths.*` → `pyazr.physics` and `Session.transform` (the engine's numbers).
- `nds`, `exfor_coulomb`, `exfor_sfactor`, `masses` → `pyazr.sources.*`.
- `tabulate`, `angular_distribution` (both write temporary projects and open new sessions) → `Session.evaluate` on
  energies given as an argument.
- `modelavg`, `priors`, `bands` → `pyazr.bands` and `Parameters.read/write`.
- The 20 bindings pyazr does not use (section 1.2) and the old-format words of the docstrings go.

### 2.4 Arrays: one rule

**Read-only numpy views where the engine owns the memory and its lifetime is the session's; copies otherwise.** A
view has `flags.writeable == False`, and its `base` keeps the engine's buffer alive (a `py::capsule` holding a
`shared_ptr` to the set-up arrays), so it stays valid and unchanged after `close()` or a rebuild; a rebuild makes new
buffers and never writes into old ones. The results of an evaluation are written by the next evaluation into the same
engine buffers, so they are copied into the `Evaluation`: a view there would change under its reader without a sign.

- Views: the energies and angles of a data set (as the engine uses them: c.m., after unit conversion), the data
  values and errors (in the engine's units), the sub-point grid E_n and the kernel rows of a data set.
- Copies: model values, residuals, χ² per data set (in the `Evaluation`); the parameter vector `x0` and any `x`
  handed back; the gradient and the Jacobian; the covariance of a fit; level and pair tables (they depend on `x`).
- Built once per process from the engine's immutable tables: the schema, the kinds, the options, the physics
  functions' descriptions (Python dicts and lists).

Today every array is a copy (section 1.2).

### 2.5 Errors and self-description

**Errors**: one base class `pyazr.Error` and five subclasses that also derive from the builtin a user catches.

- `FormatError` (`ValueError`): the project text or a record edit is refused. The message: one line per diagnostic,
  `file:line: record.field: CODE: message`; `.diagnostics` holds them (as today).
- `ParameterError` (`ValueError`): a vector of the wrong length, a NaN, an unknown name or kind. The message: the
  length expected and given, the first wrong entry by its parameter name, the call that gives a valid vector.
- `EngineError` (`RuntimeError`): the engine refuses a set-up or an evaluation (today `AZURE2Error`). The message:
  the engine's, and the record, level or pair concerned.
- `DomainError` (`ValueError`): an argument of a physics function outside its domain: the function, the element's
  index and value, the domain.
- `SessionClosed` (`RuntimeError`): a call on a closed session; which call; that a new session is needed.

Every message says (1) what was wrong with the value given, (2) what is allowed, (3) which call shows or fixes it, and
names objects by their label, not by their repr. An engine notice (printed on standard output today, as "Could not
write to EC Amplitude File.") becomes `warnings.warn(..., pyazr.EngineWarning)` with its cause and consequence. A
test runs a table of wrong calls (those of section 1.5 among them) and asserts the class and the three parts.

**Self-description**: what an agent can ask, all plain data.

- `pyazr.describe()`: version, engine commit, build options (MCMC, NLopt, readline), the Coulomb table in use (path
  and content hash, or "exact"), the SRIM data in use, the mass table, the thread count.
- `pyazr.kinds()`: record kinds with fields and units; data kinds with units, frames and extra outputs; parameter
  kinds with units and owner fields; session options with type, default, CLI flag and doc; physics functions with
  arguments and units.
- `pyazr.schema()` (as `AZURE2 --schema`), `Project.describe()`, `Project.validate()` (as today).
- `Session.describe()`: pairs, levels, data sets (kind, points, units, frame), extrapolations, parameters (name,
  kind, unit, value, free, bounds), options in effect, output directory.
- `Parameters.table()`, `Evaluation.table(i)`: readable text.
- `help(pyazr)`: a module docstring of about thirty lines naming the five objects, with a ten-line session.

### 2.6 The fit

pyazr holds no fitting framework (`ARCHITECTURE.md` §11). It hands out an `Objective`, a callable with its gradient,
bounds and names:

```python
f = s.objective()
r = scipy.optimize.minimize(f.value_and_grad, f.x0, jac=True, bounds=f.bounds, method="L-BFGS-B")      # scipy
r = scipy.optimize.least_squares(f.residuals, f.x0, jac=f.jacobian)                                  # scipy, LSQ
m = iminuit.Minuit(f, f.x0, grad=f.grad, name=f.names); m.errordef = f.errordef; m.migrad(); m.hesse()  # iminuit
sampler = emcee.EnsembleSampler(nwalkers, len(f.x0), f.log_prob, pool=pool)    # one session per worker process
```

The three recipes are examples under `pyazr/examples/` and tests that skip without the package (iminuit and emcee are
not installed on this machine; the scipy recipe ran here, use case (a)). `bounds` come from the engine (the limits of
`ParameterLimitsManager`, the nuisances' `grid_max`); a parameter without one is `(-inf, inf)`. A `Session` is not
pickled: a pool's initializer opens one per process. The engine's own fitters (MIGRAD with or without the analytic
gradient, LM, GSL-LM, NLopt; MINOS; the covariance band) stay the reference behind the program; whether pyazr offers
them as one call is decision 8. The MCMC of `numcmc` stays in the program until Phase 9; pyazr's sampler is emcee.

### 2.7 The physics explorer

`pyazr.physics`: free functions on numpy arrays (broadcasting, units in the argument names), each calling one
function of `physics/` or `stopping/`, with no session.

- Kinematics: `sommerfeld_eta`, `wave_number`, `lab_to_cm_energy`, `cm_to_lab_energy`, `lab_to_cm_angle` (with the
  Jacobian). Library: `ChannelFunctions.h`, the `EPoint` conversions. Bound today: no.
- Coulomb: `coulomb_waves` (F, G, F′, G′ on ρ, every method of pycoulomb by name), `coulomb_phase` (σ_l),
  `whittaker_wave`. Library: `physics/coulomb` (`coulomb`, `coulomb_phase`, `Whittaker.h`). Today: Whittaker only.
- Channel functions: `penetrability`, `shift`, `hard_sphere_phase`, `shift_derivative`. Library:
  `ChannelFunctions.h`. Today: through a session only.
- Widths: `wigner_limit`, `single_particle_width`, `weisskopf_width`. Today: Python (`widths.py`).
- Resonance quantities, Coulomb scattering, S factor: `statistical_factor`, `resonance_strength`, `branching_ratio`,
  `rutherford_cross_section`, `mott_cross_section`, `sfactor_conversion`, `modified_sfactor_conversion`. Today:
  bound, without pyazr names.
- Transfer and THM: `bound_state`, `distorted_wave`, `thm_vertex_function` (W_l(p)), `momentum_distribution`,
  `global_optical`. Library: `physics/transfer`. Today: the global sets only.
- Stopping: `pyazr.stopping` as today.
- Per session, where the level scheme is needed: `Session.channel_functions(pair, e_cm, l)` (today
  `coulomb_functions`), `Session.level_matrix(e_cm, x)` (R and U matrix elements, eigenphases, poles),
  `DataSet.grid()`, `DataSet.kernel()`.

**The rule that makes a new function one place and a test.** The binding holds one table of physics functions in
`api/src/pybind_physics.cpp`: a row has the Python name, the argument names with their units and domains, the
docstring, and an adapter to the C++ function on arrays (one template that loops element-wise and returns a status per
element). `_azure2.physics_functions()` returns the rows' descriptions; `pyazr/physics.py` builds one Python function
per row at import, with its signature, docstring and the `DomainError` of a status, as `Project.add_<kind>` is built
from the schema today. `pyazr.kinds()["physics"]` and the documentation table come from the same rows. A consistency
check requires a test naming each row (against mpmath, a published table, or the engine's own use on a project).

### 2.8 The command line `pyazr`

Worth having in 3.0 (each a thin call of a pyazr function; `--json` for agents): `pyazr describe FILE`,
`pyazr validate FILE`, `pyazr convert OLD NEW`, `pyazr schema`, `pyazr kinds`, `pyazr import exfor ID [--to FILE]`,
`pyazr import ensdf NUCLIDE`, `pyazr stopping MATERIAL ION ENERGY...`, `pyazr coulomb-table build|fetch|check`. The
program `azure2` itself is in the wheel. Waiting for the UI (Phase 9): `pyazr ui`, `pyazr libraries update`.

Built (`pyazr/cli.py`, `python -m pyazr` until packaging gives an entry point): one table of commands, `COMMANDS`
(name, arguments with types, the pyazr function called, help line), from which the parser, every `--help`, `pyazr
commands --json` and the reference page (`docs/source/user_guide/pyazr_command_line.rst`) are made. `describe [FILE]`
(`pyazr.describe`; without FILE the package), `validate`, `convert`, `schema [KIND]`, `kinds [TABLE]`, `stopping`,
`coulomb-table check [FILE]`, `ui FILE`, `import exfor ID` (`fetch_exfor`, the local EXFOR master files with
`--source local`), `import ensdf NUCLIDE` (`fetch_levels`), `commands`. `--json` gives one document on standard output
(`command`, `ok`, `exit`, `result`, and `error` with `what`, `allowed`, `fix`); the engine's output goes to standard
error meanwhile. Exit codes those of the program's project commands: 0, 1 a failed run, 2 a wrong call or a refused
input, 3 an old project. `describe` and `validate` call the engine functions the program calls
(`experiment::ReadFormat3Text`, `DescribeJson`/`DescribeText`) and give the same JSON on T2; the program's
`--validate` and `--describe` also report the stopping of each target (`ReportTargetStoppings` in
`src/ProjectCommands.cpp`), which the binding does not. Not built, as rows that exit 2 with what to do instead:
`coulomb-table build` (the engine's `CoulombTable::BuildFile` is not bound; the program builds the table) and
`coulomb-table fetch` (no published copy). Not here: `pyazr libraries update`; `import exfor --to FILE` (pyazr has
no single call that writes any EXFOR entry as a data set: `to_dataset` takes ratios to Rutherford or Mott,
`to_sfactor_dataset` S factors, `to_azr` the old format). A `pyazr evaluate FILE` belongs here once `Evaluation`
exists (a call and a printer: the numbers at the project's parameters, for an agent without a session).

### 2.9 Packaging

A wheel carries the module `_azure2` (the core linked in), the program `azure2`, the mass table
(`pyazr/data/mass_1.mas20`, 464 kB) and the SRIM-2013 data (two files, 388 kB, compiled into the library today) under
decision 6; not the Coulomb table (85 MB, decision 7). GSL is linked statically (GPL, as AZURE2); the OpenMP runtime is
vendored by the platform's repair tool (auditwheel, delocate, delvewheel); readline is not needed once the menus are
gone (decision 3).

What scikit-build-core and cibuildwheel need from the CMake files as they are on this commit:

1. An `install(TARGETS _azure2 LIBRARY DESTINATION pyazr)`: there is none; the module is written into the source
   tree (`LIBRARY_OUTPUT_DIRECTORY ${PYAZR_OUTPUT_DIR}`, `api/src/CMakeLists.txt`) and shipped as setuptools package
   data (`pyproject.toml`).
2. `install(TARGETS AZURE2 ...)` exists (`src/CMakeLists.txt`, `${CMAKE_INSTALL_DIR}`); a wheel installs it into
   `${SKBUILD_SCRIPTS_DIR}`.
3. `USE_MCMC` defaults to ON and stops the configuration when `numcmc/` is missing, `USE_READLINE` to ON: a wheel
   build sets both OFF.
4. The test directories are added unconditionally (`add_subdirectory(tests/...)` in the top `CMakeLists.txt`):
   guard them with `BUILD_TESTING`.
5. `pyproject.toml` builds with setuptools and reads the version from `pyazr/__init__.py`: scikit-build-core's
   backend with its regex version provider on the same line.
6. cibuildwheel: GSL before the build (`gsl-devel` in manylinux_2_28, Homebrew on macOS). On Windows the code is
   built with MinGW (MSYS2, `.github/workflows/build.yml`), while cibuildwheel uses MSVC: either a MinGW wheel (a
   build step of its own, delvewheel for the runtime DLLs) or an MSVC port of the core; this is the largest unknown
   of the packaging. The wheel's test command runs T1–T3 on the installed wheel.

Not testable on this machine: the macOS and Windows wheels, and installing any wheel (no installation is allowed
here). They are tested in a workflow of their own on the three runners of `build.yml`, which builds and tests wheels
and publishes nothing; the release files are not touched.

## 3. Extensibility of the binding, today and after this design

Measured by searching `api/` and `pyazr/` for the names of existing kinds (`neutron_yield`,
`beta_delayed_spectrum`, `resonance_strength`; `gamma_decay`, `detector`, `energy_scale`; `radiative`, `cbkg`,
`thickness`, `doppler`, `apparatus`, `shift_fn`) and reading the places found.

**A new data kind (quantity).** Today: none, when its model value is one number per point. `Project.add_dataset`
takes the quantity from the engine's schema (`quantities`, the quantity table of `experiment/`, `project_format3.md`
§7), the session evaluates it, `DescribedDataset` reads the engine's description; `neutron_yield` and
`beta_delayed_spectrum` are named nowhere in `api/` or `pyazr/`. A kind with more than one number per point needs a
method per shape: `calculate_angular_dists[_rwa]`, `calculate_analyzing_power[_rwa]`, the binding's
`calculated_segments_e1/e2` and `calculated_angular_dists`. After: none; the quantity's row declares its extra
outputs and `Evaluation.data[i].extra` returns them by name.

**A new record kind.** Today: none for the record (`Project.add_<kind>` is generated from the schema at import,
`project.py`; `validate` and `describe` are the engine's); if the record carries fit parameters, the places of a
parameter kind. After: none.

**A new fit-parameter kind.** Today: the engine's row (`ParameterKind.h`; bound generically by `parameter_kinds()`
and `parameter_records()`) and five places of the binding: (1) `Parameter` (`pyazr/parameters.py`) has one attribute
per owner field (`decay_record`, `decay_part`, `final_state` were added for `radiative`); (2) `ParameterSet` has a
hand-written view per kind (12; none for `radiative`, `apparatus`); (3) `Parameter.__repr__` branches per kind; (4)
the numeric layout of `GetParameterInfo` (16 doubles, "keep in sync" with `Parameter._NFIELDS`) with type codes read
by number in `azure2.py` (`kind <= 1` or `kind == 6`); (5) per-kind branches in `snapshot.py`, `modelavg.py`,
`bands.py` for the kinds they treat. After: the engine's row only; `Parameter` is the engine's record,
`Parameters.of_kind(key)` serves every kind, the row's label gives the text; the numeric layouts of
`GetParameterInfo` and `GetPairsInfo` are deleted.

**A new physics function.** Today: five places and a test: `api/src/pybind_physics.cpp` (an `m.def` and an array
loop), `pyazr/physics.py` (a wrapper: broadcasting, status to exception, docstring), its `__all__`,
`pyazr/__init__.py` for a top-level name, `docs/source/user_guide/pyazr.rst`. Six of the nine functions bound so far
have no pyazr name. After: one row of the binding's function table and a test (section 2.7).

**A new session option.** Today: `RuntimeOptions` in `pybind_api.cpp` (the struct, `apply_options`,
`def_readwrite`), the constructor's signature, docstring and `self.options` in `azure2.py`, and the program's flag in
a separate list (`src/AZURE2.cpp`). After: one row of an options table read by the program's parser and by the
binding (`pyazr.kinds()["options"]`).  Built 2026-11-13 (`project_format3.md` §3.12): an option of the calculation
is one row of `io::ConfigKeys()` (key = keyword, flags, default, effect: the program's parser and help, `--describe`,
`kinds()` and the UI read it), one row of its place in `Config` (`src/ProjectOptions.cpp`), and its keyword in the
constructor's signature and docstring; `RuntimeOptions.options` takes it by name.

## 4. The order of work

Each step leaves `ctest`, `tests/run_tests.sh` and every project's outputs as they are: the binding changes, the
engine's physics does not. Sizes are estimates from the files of this commit.

1. **Errors and self-description.** `pyazr/errors.py` (the hierarchy; `AZURE2Error` registered as `EngineError`),
   the package docstring, `pyazr.describe()`, `pyazr.kinds()` (adding the options table and `physics_functions()`
   to the binding), `Session.describe()` on today's session. New `errors.py`, `describe.py`; `__init__.py`;
   `pybind_api.cpp` (+60). About 500 lines; tests `errors_test` (the wrong calls of 1.5), `describe_test`.
   **Done** (branch `feature/upgrade-8-pyazr3-errors`, additions only).  Built: `pyazr/errors.py` (`Error`;
   `FormatError` moved there from `project.py`, one class; `ParameterError` with `ParameterKeyError`,
   `ParameterIndexError`, `ParameterVectorError` for the places that raised `KeyError`, `IndexError` or
   `AZURE2Error` before; `AZURE2Error`, also named `EngineError`, the class `_azure2.AZURE2Error` is, keeping the
   name `AZURE2Error` that two tests compare; `DomainError`, `SessionClosed`, `MissingFile` (a
   `FileNotFoundError`), `DomainZeroDivision` (also a `ZeroDivisionError`, where a Python formula divided by zero
   before), `EngineWarning`); the binding raises the Python classes (`BindErrors`) and checks every
   free vector for a NaN or an infinity, naming the parameter; `pyazr/session_checks.py` (the argument checks and
   messages of the session's calls); the engine's notices pass through a line buffer per session
   (`api/src/pybind_notices.cpp`, the session's `Config::outStream`), whose table of four notices (EC amplitude
   file written or appended, normalization and shift files) become `EngineWarning`s with the file and the
   consequence; `pyazr.describe()` without a path, `pyazr.kinds()`, `azure2.describe()` (in `described.py`, which
   was extended rather than doubled by a `describe.py`: a module `pyazr.describe` would hide the function of that
   name); `api/src/pybind_describe.cpp` (`build_info`, `coulomb_table_info`, `session_options`) and
   `_azure2.physics_functions` (the names `BindPhysics` registers); the package docstring with a tested session;
   the stale docstrings of 1.1 and 1.6 corrected.  Tests `errors_test` (27 wrong calls), `describe_test`.
   Left: the options table is the binding's (the program's parser and `--help` in `src/AZURE2.cpp` are not a
   table; a test ties the three together); `kinds()["data"]` has no extra outputs (the engine's quantity table
   does not declare them); `azure2.describe()` gives no parameter bounds (the binding does not reach the engine's
   limits); notices printed outside `Config::outStream`, and the engine's `WARNING:` lines (which tests read on
   standard output), stay on standard output; the physics table is names and pybind signatures (step 4).
2. **Parameter records instead of numeric layouts.** `parameter_records()` and a `pair_records()` as dicts;
   `Parameter` and `Parameters` built on them; `of_kind`; delete `GetParameterInfo`, `GetPairsInfo`, the two
   `_NFIELDS` and the type codes. `AZUREAPI.h/.cpp` (−200), `pybind_api.cpp`, `parameters.py` (795 → about 400),
   the readers of the layout in `azure2.py`, `snapshot.py`, `bands.py`, `modelavg.py`. About 900 lines changed;
   `parameter_kinds_test` and the tests that use `parameters`.
   **Done** (branch `feature/upgrade-8-pyazr3-records`).  Built: the record table `api/src/pybind_records.cpp` (50
   parameter fields, 22 pair fields: name, unit, the kinds a field describes, whether the table shows it as the
   owner, doc, the reader of the engine object), bound as `Session.parameter_records()`, `Session.pair_records()`,
   `_azure2.parameter_record_fields()`, `pair_record_fields()` and in `pyazr.kinds()["parameter_fields"]`,
   `["pair_fields"]`; each kind's owner fields in `parameter_kinds()` come from it; the parameters' limits as created
   (`AZUREAPI::GetParameterBounds`).  `pyazr.Parameters` (`x0`, `names`, `kinds`, `units`, `free_mask`, `bounds`,
   `free`, `of_kind`, `of_level`, `by_name`, `put`, `take`, `read`, `write`, `table`, `to_records`; a selection is a
   `Parameters` with the same `x0`), `Parameter` and `Pair` the records themselves (a field of another kind reads
   `None`), `ParameterSet` = `Parameters` plus the views of earlier releases;
   `pyazr.snapshot.write_values(session, project, x)` (a `Project` edited in memory, web UI wish 8).
   Removed: `GetParameterInfo`, `GetPairsInfo`, `kParamInfoFields`,
   `kPairInfoFields` (301 lines of `AZUREAPI.h/.cpp`), `parameter_info`, `pairs_info`, the two `_NFIELDS`, the
   numeric type codes (in the records, `parameter_kinds()`, `kinds()`, `bands.py`, `azure2.py`), the fall-backs for
   engines with 19 or 20 pair fields, `Parameter.from_record`, `Pair.from_record`.  Differs from the plan: `free` is
   the selection of the free rows (as before; the per-row flag is `free_mask`); `put(x, values)` and `take(x)` of
   earlier releases are kept beside `put(x, rows, values)` and `take(x, rows)`; five names of earlier releases read
   their field (`prior_nominal` → `prior_centre`, `input_is_rwa` → `width_is_amplitude`, `density_max` and
   `doppler_max` → `layout_max`, `Parameter.sample` → `sample_name`); `bounds` are the limits the engine creates a
   parameter with (a `<parameterSettings>` block's limits act in the program's fit only); the prior of every kind
   with one (norms and the three shift terms too) is the declared prior, read where the objective reads it;
   `parameters.py` 804 → 620 lines (it keeps `Pair` and `NuclearPotential`).
3. **`Objective` and `Evaluation` beside today's methods** (added, nothing removed), and the three fit recipes as
   examples and tests. New `objective.py`, `evaluate.py` (about 600 lines), three examples. About 850 lines; tests:
   the object against `objective_and_grad`, `Evaluation` against `chiSquared.out` on T1–T3, the recipes (scipy runs;
   iminuit and emcee skip without the package). Needs 2.  Design and as built: section 7.
4. **The physics table** and `pyazr.physics` generated from it; bind the channel functions, the Coulomb waves and
   phase, kinematics, widths, the transfer functions; `Session.channel_functions` and `level_matrix`.
   `pybind_physics.cpp` (240 → about 700), `physics.py` (generated, about 150). About 900 lines; one test per row
   (about 20) and the consistency check. Needs 1 (for `DomainError`).
   **Done** (branch `feature/upgrade-8-pyazr3-physics`).  `api/src/pybind_physics.cpp` (about 760 lines): 24 rows
   (kinematics 6, Coulomb 3, channel functions 4, THM plane-wave vertex 1, widths 3, resonance quantities, Coulomb
   scattering and S factor 7;
   the seven bound one by one before moved into the table), one template with a status per element (0; argument k
   outside its domain; -2 the arguments together; -1 not converged); `_azure2.physics_functions()` returns the rows
   and `_azure2.physics_evaluate` evaluates one.  `pyazr/physics.py` (170 lines) builds a function per row;
   `errors="mask"` is the keyword for NaN and a status array.  `kinds()["physics"]` is the rows (with `axis`, the
   argument a plot puts on its horizontal axis, for the web UI), `kinds()["stopping"]` the stopping functions.
   `physics/Widths.h` adds `wigner_limit`, `single_particle_width`, `photon_penetrability` (the engine calls the
   first and the last).  `Session.channel_functions` in `pyazr/session_physics.py` (a mixin; the rename step moves
   it into `session.py`).  Differs from the plan: `level_matrix` not written (the engine hands out no R or U matrix
   without new code in `AZUREAPI`: the collision matrices live inside the evaluation), the transfer functions on
   potentials not bound (they take potentials and return functions of r, not one number per element), no
   Weisskopf width (no engine function: new physics, a decision); the Python formulas are untouched (the rename
   step deletes them).
5. **Rename and remove.** `Session` (`azure2.py` split into `session.py`, `data.py`, `system.py`), one parameter
   basis, the `_rwa`/physical pairs, `calculate_*`, the mode switch, the per-session edits; delete the formulas of
   `transform.py` and `widths.py`, the 20 unused bindings, the old-format words. Port the tests (105 files of
   `tests/pyazr/` call `azure2(`; 28 scripts elsewhere in `tests/` import pyazr), the 19 examples and the user guide
   (`pyazr.rst`, 1303 lines). The largest step: about 4000 lines of `pyazr/` and 2500 of `tests/` touched; done
   by one agent in one stride. Needs 2, 3, 4.
   **Begun** (branch `feature/upgrade-8-pyazr3-session`; the measured table of every name is its HANDOFF entry).
   Built: `pyazr/session.py` with the class `Session` (no `pyazr.azure2`, no alias; options keyword-only), the THM
   methods in the mixin `session_thm.py`, the session's lock on every public member that reaches the engine
   (`pyazr/locking.py`, `tests/pyazr/session_lock_test.py`); tests, examples, the web UI and the pages ported.
   Differs from the plan: no `threads=` (the engine has no per-session thread count); `azure2.py` is not split into
   `data.py` and `system.py` yet.  Not yet done: the evaluation family, one parameter vector and table, the Python
   physics, the libraries.  Measured on T2 for the ports to come: `segment_chi2` and `dataset_chi2` split a
   data-covariance group otherwise than `Evaluation.data[i].chi2` (equal sums); `chi2_and_grad` and
   `residual_jacobian`'s `skip` have no replacement; the engine's evaluation runs in extrapolation mode (models bit
   for bit those of `calculate_rwa`), so `evaluate` can serve the extrapolations without engine code.
6. **Read-only views of set-up arrays** (`DataSet.energies`, `values`, `errors`, `grid`, `kernel`): `shared_ptr`
   set-up buffers in the engine's data, capsule-based arrays in the binding. `AZUREAPI.cpp`, `pybind_api.cpp`, the
   data access of `experiment/`. About 400 lines; a view test (read-only, alive after `close()`, unchanged after a
   rebuild). Needs 5.
7. **The `pyazr` command line** (`pyazr/cli.py`, an entry point). About 350 lines; a test running each command once
   with `--json`. Needs 1.
   **Done** (branch `feature/upgrade-8-pyazr3-cli`; section 2.8 lists what is built).  `pyazr/cli.py` (about 720
   lines: the table, a call per command, the printer, the parser), `pyazr/__main__.py`; no entry point yet (step 8).
   Tests: `tests/pyazr/cli_test.py` (ctest `pyazr_cli`: one case per row with `--json`, equality with the program on
   T2, wrong calls, `--help`), consistency check 42 (a row without a case, a stale reference page).  Differs from the
   plan: `coulomb-table build` is not available (no bound build function: decision 7's "built on first use" is not
   implemented either); a missing file exits 2 where the program's `--validate FILE` exits 1 for an unreadable file.
8. **Packaging**: scikit-build-core `pyproject.toml`, the CMake install rules and options of 2.9, a cibuildwheel
   configuration and a workflow of its own building and testing wheels on Linux, macOS, Windows (publishing
   nothing). `pyproject.toml`, three `CMakeLists.txt` (about 60 lines), the workflow (about 120). Needs decisions 5,
   6, 7.

Side by side: 1, 2, 4, 7 and 8 touch different files and can run at once (1 and 2 both touch `pybind_api.cpp`, in
different functions: merge 1 first). 3 needs 2; 5 needs 2, 3, 4; 6 needs 5. The critical path is 2 → 3 → 5 → 6.
The engine's defect of section 1.6 (the isotope prior) is a branch of its own on the engine, before step 3's test
of `Evaluation` against `chiSquared.out`.

## 5. Decisions for the PI

Each with the options, the recommendation and what it costs to change later. Nothing a user can set is removed by
this note; the removals of section 2.3 wait for these answers.

**Taken by the PI on 2026-11-13** ("I accept your recommendations for pyazr", on the director's list): 1 (a), 2 (a),
3 (a) for the menus and (b) for `numcmc`, 4 all once bound, 5 (a), **6 (a)** (the SRIM data compiled in, with the
citation notice; the PI's word on the terms: free to use as long as SRIM is cited; the owner is to be told),
7 (c) with (b) as the offline path, 8 (b), 9 (a), 10 (a), 11 (a). Decision 6 differs from this note's
recommendation below; the others are the note's. The removal of the readline menus (3) is a change for the owner's
users and stands on the list of things he must see.

1. **Names of the top-level objects.** (a) `Session`, `Project`, `Parameters`, `Objective`, `Evaluation` (the
   binding's class is already `_azure2.Session`); (b) `Model` for the session; (c) keep `azure2`, a lower-case class
   named after the program, whose `help()` shows a mixin's docstring today. *Recommendation: (a).* Later: a rename
   touches every script (105 test files call `azure2(`); cheap before 3.0, dear after.
2. **Views against copies.** (a) read-only views of set-up arrays, copies of results (2.4); (b) copies everywhere
   (today); (c) views of results too. *Recommendation: (a).* Later: (a) and (b) have the same signatures and swap
   freely; (c) cannot be taken back without breaking scripts that keep a result.
3. **The readline menus and `numcmc`: with 3.0 or with Phase 9.** (a) 3.0: the program takes subcommands and flags
   (`azure2 calculate|fit|extrapolate|minos|rate|mcmc FILE`), the menus go; (b) both wait for Phase 9 (the roadmap's
   order). *Recommendation: (a) for the menus* (an agent cannot drive a menu safely: `src/AZURE2.cpp` documents the
   y/n trap of scripted answers; the wheel then needs no readline), *(b) for `numcmc`* (its replacement is the UI's
   sampler and the emcee recipe). Later: small either way.
4. **Which Python physics duplicates are deleted for bound functions.** The 17 formulas of 1.3, the penetrability
   workaround of `wigner_widths` and the rate integral of `examples/reaction_rate.py`. *Recommendation: all, in step
   5, once step 4 binds the engine's functions*; keep `masses`, `nds` parsing, `modelavg`, `priors` (library
   reading, statistics, file formats). Later: none; a name can come back as a bound function.
5. **Minimum Python and numpy.** Today `>=3.8`, `numpy>=1.20`; this machine has 3.12 and numpy 2.5. (a) Python 3.10,
   numpy 1.24 (3.8 and 3.9 are out of support); (b) keep 3.8. *Recommendation: (a)*, the wheel built against numpy 2
   (a pybind11 module so built runs with numpy 1.x), CI on the oldest and newest. Later: raising is cheap, lowering
   is not.
6. **The SRIM data's terms in the wheel.** The two SRIM-2013 files are not under AZURE2's licence
   (`stopping/data/srim/README.md`). (a) compiled in, as in the source build, with the notice; (b) a wheel without
   them, `stopping evaluation` refused with a message naming `AZURE2_SRIM_DATA_DIR`; (c) a separate data package with
   its own licence text, found by pyazr when installed. *Recommendation: (c) until the authors' permission for
   binary redistribution is in writing, then (a).* Later: none for the API (`pyazr.describe()` names the data used).
7. **The Coulomb table (85 MB).** (a) in the wheel (close to PyPI's file-size limit); (b) built on first use (about
   5 min on 4 threads); (c) fetched on request from a release asset, checked by its content hash (the table is named
   by content already); without it the exact routines. *Recommendation: (c), with (b) as the offline path.* Later:
   small.
8. **The engine's own fitters from Python.** (a) none: pyazr gives the objective, the program fits; (b) one call,
   `Session.fit_with_engine(method="migrad"|"lm"|"gsl-lm", minos=False)`, giving the program's numbers (MIGRAD's
   covariance, MINOS intervals). *Recommendation: (b)*: not a framework, and the GUI and agents get the reference
   fit and MINOS without a subprocess. Later: the same cost as now.
9. **The roadmap's "`azure2()` session shim with deprecation warnings for one release"** conflicts with "no backward
   compatibility". (a) strike it; (b) keep a shim for one release. *Recommendation: (a)*; section 2.3 is the
   migration list. Later: none.
10. **One parameter vector in the API.** (a) the fit vector only, physical values through `Session.transform(x)`;
    (b) keep both families (`*_rwa` and physical). *Recommendation: (a)*: half the session's methods go, and the
    basis (Brune, Park, formal) stays where the engine settles it. Later: dear (every script).
11. **The per-session edits** (`set_channel_radius`, `set_nuclear_potential`, the mode switch). (a) edits of the
    `Project` and a new session (the project is the state, as the web UI will have it); (b) keep them on the
    session. *Recommendation: (a)*; a radius scan costs one session per radius, as with the program. Later: small.

## 6. Not looked at

The Qt GUI and the web UI's design (Phase 9); the Sphinx build of `docs/source`; the four slowest tests of
`tests/pyazr/` under the tracer (`thickness`, `concurrent_sessions`, `doppler`, `session_isolation`) and
`upgrade_thread_independence` (stopped: four threads); the 28 scripts outside `tests/pyazr/` that import pyazr (only
counted); the examples (not run, `fit_emcee.py` and `fit_zeus.py` cannot run here); `nds.py` against the network;
Python 3.8–3.11 and numpy 1.x (not installed); any wheel build; macOS and Windows; the MCMC and MINOS paths.

## 7. The objective and the evaluation (step 3 of section 4)

### 7.1 As designed for step 3

Added beside today's methods; nothing is removed or renamed.  Two entry methods of today's session class:
`azure2.fit_objective()` (an `Objective`; the name `objective` is today's method `objective(x)`, so the rename step
gives the object its final name) and `azure2.evaluate(x=None)` (an `Evaluation`; the name is free and is the final
one).  Both refuse a session in extrapolation mode (there are no data and no objective there).  Everything is on the
free vector x of the basis the session fits in (decision 10).

**Which engine call gives each field.**  No formula in Python: a number the program writes comes from the code that
writes it, and the binding reads it where the writer reads it (the results buffer of the evaluation).

- `Objective.x0`, `names`, `bounds`: the parameter records (`start`, `name`, `bounds` of the free rows).
- `Objective.errordef`: 1 (the objective is a χ²; ln L = −½ χ²).
- `Objective(x)`: `AZUREAPI::CalculateObjectiveRWA` (`Objective::Engine().Chi2`).
- `grad(x)`, `value_and_grad(x)`: `CalculateObjectiveGradRWA` (one call gives both).
- `residuals(x)`, `jacobian(x)`: the rows of `CalculateResidualJacobianRWA` (data, whitened by a data covariance;
  model data; prior covariance), then the rows the program's least-squares drivers carry: `AppendModelRows` (Park's
  wall, the `gamma_decay` priors and wall, the abundance sums) and one row (p − c)/σ per diagonal prior
  (`ForEachDiagonalPrior`) and per nuisance prior, so that Σ r² is the objective.  New: `AZUREAPI::ObjectiveRowsRWA`.
- `log_prob(x)`: −½ × the objective, −∞ outside the sampler's support (below).  New: `Session.log_prob_rwa`.
- `Evaluation`, all of it: one model evaluation (new `AZUREAPI::EvaluateObjectiveRWA`: fill, the model of every
  point, the objective's sum and its breakdown by term from the same model values), then reads of the points and of
  the results buffer, as the output writer reads them.  New: `Session.evaluation_rwa`.
  - `objective`: the sum `Objective(x)` returns, bit for bit (the plain sum; the breakdown is summed apart).
  - `chi2`: Σ of the per-segment χ² the evaluation stores, in segment order: `Total-Chi-Squared` of `chiSquared.out`
    (a data covariance by its per-segment partition).
  - `terms`: name → value, one per `ObjectiveTerm` (`ObjectiveTerm::Name`).
  - `data[i].chi2`, `.norm_chi2`, `.shift_chi2`: the row of `chiSquared.out` (`Results::objective.segmentChi2`), the
    segment's normalization prior and its three shift-prior parts summed (`ObjectiveContext::SegmentPriors`).
  - `data[i].energies`, `.angles`, `.excitation_energies`: the points (`GetCMEnergy` after the shift, `GetCMAngle`,
    `GetExcitationEnergy`).
  - `data[i].model`: the model value of each point, engine unit, c.m. (`Unit::GetModelValue`, what `AZUREOut_*`
    prints).
  - `data[i].values`, `.errors`, `.norm`: the datum and its error as the engine compares them (c.m., engine unit,
    `Segment::CMValue`, `CMError`) and the norm in use (`WorkingCopy::NormUsed`: a profiled norm at n*).
  - `data[i].residuals`: (f/J − y n)/(dy n) of each point in the datum's frame (`azure2::PointResidual`), not
    whitened.
  - `data[i].extra`: the rows of the binding's table of point outputs that apply to the data set (below).  New: the
    table.
  - `model_data`: the `<modelData>` rows (`AZUREAPI::ModelDataRWA`, as `pyazr.model_data`).
  - `log_prob`, `excluded`: as `Objective.log_prob`; whether a term puts x outside the sampler's support
    (`Objective::ExcludesForSampler`).

The names and kinds of the data sets of `data[i]` (name, quantity, engine unit) come from the session's description
of the project (`active_datasets`), in the engine's order; the arrays from the engine.  The data sets are the
engine's reported segments (the first of each run of equal keys, as `calculate_rwa`); a data set's χ² is the sum of
the rows of its segments.

**The term names.**  `Evaluation.terms` is a dict whose keys are the names of the engine's table of terms
(`Objective::Terms()`, `ObjectiveTerm::Name`, bound once through `ObjectiveTermNames`), in its order: a new term of
`src/Objective.cpp` appears in `Evaluation.terms` with no edit of `pyazr/` (section 7.2 shows it).  The engine's
first term holds the data, the normalization priors and the energy-shift priors together; the per-data-set split of
the latter two is `data[i].norm_chi2` and `.shift_chi2`, as `chiSquared.out` splits them.

**The columns of a data set.**  The engine's quantity table (`experiment/Quantities.cpp`) does not see the points
(it is the `experiment/` layer, without engine headers), and which extra numbers a point carries depends on the data
set (a capture reaction has its E1 and E2 parts, a sample its per-nucleus terms, an angular distribution its Legendre
coefficients), not on the quantity alone.  So the declaration is a table of the binding, as the parameter records are
(`api/src/pybind_evaluation.cpp`): one row per point output, with its name, unit, doc, the condition on the data set
(a predicate of the segment) and the reader of the point.  `pyazr.kinds()["data_outputs"]` lists the rows;
`kinds()["data"][q]["extra_outputs"]` stays `None` (the quantity alone does not decide).  No engine change is needed
for it.  Rows of this step: the lab energy and the lab angle of a point (web UI finding 22), the S-factor conversion,
the E1 and E2 parts, the Legendre coefficients of an angular distribution, a sample's per-nucleus terms.  The σ and S
columns of `AZUREOut_*` (finding 23) are formed by the writer (`SigmaAndS` of `DataOutputFiles.cpp`); serving them
means moving that struct where the writer and the binding both call it (an engine change, outputs identical):
a later commit if time remains.

**`log_prob`, one rule with the engine's sampler** (`objective.md` section 6): ln P(x) = −½ × the objective, with no
constant (the sampler's `LogPriorNormalisation` is a constant and is left out), and −∞ where the engine's sampler
excludes the point: a term says so (`Objective::ExcludesForSampler`: Park's J < 0), or x lies outside the limits of a
parameter whose kind's row has `samplerSupport` (today the `gamma_decay` record's numbers).  The binding makes both
tests (`Session.log_prob_rwa`); the column reaches Python as the parameter record's field `sampler_support`, so a
user sees which bounds are walls.  The limits of the other kinds are not turned on (the PI's open decision).  A user
who wants walls at every bound writes one line:
`lp = lambda x: f.log_prob(x) if all(lo <= v <= hi for v, (lo, hi) in zip(x, f.bounds)) else -numpy.inf`.

**Arrays** (section 2.4): every array of an `Evaluation` and every array an `Objective` returns is a copy the caller
owns.  The set-up arrays of a data set (energies, angles, data values and errors) are copies too in this step, marked
so in their docstrings; they become read-only views in step 6.

**Threads and state.**  An `Objective` holds its session and copies of x0, names and bounds, nothing else (no cache:
the engine keeps its own state of the last evaluation, and that is the session's).  An `Evaluation` holds copies and
no session.  Calls from two threads on one session are safe: every engine call of `Objective` and `Evaluation`
takes the session's lock (a `threading.RLock` of the session object), so they run one after the other and each
returns what it would return alone (the engine evaluates on the session's one compound and data, so two calls at once
would race).  Today's methods do not take the lock; the rename step makes every engine call of `Session` take it.
Two sessions share nothing and run side by side.  A session is not pickled: a process pool opens one per worker.

**The defect 293.94 against 297.86** (fitted T2, `work/director/fits-t50-1/T2-migrad`): 297.860308551 is the
program's printed `Total Chi-Squared` at the end of the fit, the objective; `chiSquared.out` writes
`Total-Chi-Squared: 293.94`, the data χ², and beside it the priors: norm 1.20857, shift 1.15434, abundance 2.1e-7,
density 0.166164, prior covariance 1.39169; 293.94 + 1.20857 + 1.15434 + 0.166164 + 1.39169 = 297.861.
`calculate_chi2_rwa` is the data χ² plus the model terms (Park, abundances, model data), 293.94 here.  So no number
is wrong; `Evaluation.chi2` is the 293.94 (`Total-Chi-Squared`), `Evaluation.objective` the 297.86, and
`Evaluation.terms` holds the difference by name.

**Tests** (each shown to fail once with an injected error): `tests/pyazr/objective_object_test.py` (the object
against today's `objective_and_grad`, `residuals`, `residual_jacobian` bit for bit on T1, T2, T3, 13N_fast, the
nuisance-prior projects and the `gamma_decay` record with priors; Σ r² against the objective; the gradient against
a central difference; `log_prob` inside and at the excluded points, and against the engine's sampler; two threads),
`tests/pyazr/evaluation_test.py` (against the program's `chiSquared.out` and `AZUREOut_*` on T1, T2, T3; one
evaluation per `Evaluation`), `tests/pyazr/fit_recipes_test.py` (the three recipes; iminuit and emcee skip without
the package).

### 7.2 As built

As designed (7.1), with these differences and numbers.

- Files: `api/src/AZUREAPIObjective.cpp` (`EvaluateObjectiveRWA`, `ObjectiveRowsRWA`, `OutsideSamplerSupport`,
  `LogProbabilityRWA`), `api/src/pybind_evaluation.cpp` (`Session.evaluation_rwa`, `objective_rows_rwa`,
  `log_prob_rwa`, `objective_term_names`, `sampler_log_probability`; `_azure2.point_output_fields`), the record field
  `sampler_support` (`pybind_records.cpp`), `pyazr/objective.py`, `pyazr/evaluate.py`, two entry methods of
  `azure2` (`fit_objective`, `evaluate`), `kinds()["data_outputs"]` and `kinds()["objects"]`, the examples
  `fit_scipy.py`, `fit_iminuit.py` and `fit_emcee.py` (rewritten on `log_prob`).  pyazr 3.0.0.dev3.
- `Session.sampler_log_probability(x0, x)` binds the engine's sampler (`AZURECalcMCMC`, set up as the program sets
  it up at x0) in about 20 lines, so the test ties `log_prob` to it: equal to the sampler's ln P minus its constant
  on the `gamma_decay` record with priors and on T2, at the start and at a displaced point (1e-10).
- Without the session's lock, two threads calling one `Objective` crash the process (segmentation fault, shown
  with the lock removed); with it, each gets its values alone bit for bit.  Today's methods do not take it.
- One evaluation per `Evaluation`, measured on 13N: 20.4 ms against 20.0 ms for one `objective(x)` (a second
  evaluation injected: 42.0 ms).  The engine's Coulomb-function counters (`cache_stats`) stay at zero when the
  Coulomb table serves the functions, so they cannot count evaluations.
- The recipes on T2 from its start (one thread): iminuit's MIGRAD 297.86032061 in 113 calls and 19 gradients
  (HESSE brings the calls to 417), the program's MIGRAD 297.860308551; scipy's least_squares on the rows 297.859396915
  in 8 evaluations; L-BFGS-B 297.859396915 in 35 evaluations, on x divided by the inverse length of each Jacobian
  column (on x / |x0| it stood at 521 after 250 evaluations).  The program's MIGRAD stops 9.1e-4 (3.1e-6 relative)
  above the minimum both scipy recipes agree on to 2e-11; its gradient at the program's point is not small (91 in an
  energy, 65 in a width).  The test asks iminuit for the program's minimum to 1e-6 and scipy for that or lower.
- emcee: 8 walkers x 60 steps on the record (three free branchings, one under a prior centred at -0.01) stay in
  [0, 1]; the same seed gives the same chain, and a pool of two worker processes (spawned, one session each) the
  same chain as one process.

**Extensibility, before and after this step** (the places that must know):

- A new term of the objective.  Before: `src/Objective.cpp` (the class and its place in the list), the key of
  `azure2.penalties()` (`pyazr/azure2.py`, a hand-written map of names), the writer's `Total-` field
  (`src/DataOutputFiles.cpp`), its least-squares rows where the drivers need them (`AppendModelRows`), its gradient.
  After, for pyazr: `src/Objective.cpp` only; it appears in `Evaluation.terms` and in `Objective(x)`.  Its rows, if it
  is a sum of squares, are one call in `ObjectiveRowsRWA` (`api/src/AZUREAPIObjective.cpp`; the test of
  `sum(r**2) == f(x)` fails without them), its gradient in the engine's gradient paths as before.
- A new data kind with two extra output columns.  Before: the quantity row (`experiment/Quantities.cpp`), arrays in
  `AZUREAPI::ReportSegments`, two bindings (`calculated_*`), a method per shape in `azure2.py`.  After: the quantity
  row and two rows of `PointOutputs()` (`api/src/pybind_evaluation.cpp`); nothing in `pyazr/`.
- A new field of `Evaluation`.  Before: a new session method (`azure2.py`) and its binding.  After: one entry of the
  dict of `EvaluationDict` (`pybind_evaluation.cpp`), its attribute and its `FIELDS` line in `pyazr/evaluate.py`.

Not shown by a build: a term added to `src/Objective.cpp` reaching `Evaluation.terms` (the names come from
`ObjectiveTermNames`, the engine's list, read at each evaluation; the evaluation test reads every term by its
engine name).  Not done (left for later): the dense model line of a data set (finding 21), the σ and S columns
as the writer forms them (`SigmaAndS`), read-only views.
