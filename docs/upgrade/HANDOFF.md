# AZURE2 upgrade: handoff

Update this file at the end of every work session: what changed, what is next, what is open. Newest first.

## 2026-11-13 (agent, branch `feature/upgrade-8-pyazr3-objective`) — pyazr: `Objective`, `Evaluation`, three recipes

**Done** (beside today's methods; nothing removed; design and as built: `pyazr3.md` section 7).
- `azure2.fit_objective()` → `pyazr.Objective`: `f(x)` (= `objective(x)` bit for bit), `grad`, `value_and_grad`,
  `residuals`/`jacobian` (the objective as rows, `sum(r**2) == f(x)`: today's residual rows, then Park's wall, the
  `gamma_decay` rows, the abundance sums, the diagonal and nuisance priors; `AZUREAPI::ObjectiveRowsRWA`),
  `log_prob` (the engine sampler's rule: -f/2, -inf outside the limits of a kind with `samplerSupport` and where Park's
  J < 0; `AZUREAPI::LogProbabilityRWA`), `x0`, `names`, `bounds`, `sampler_support` (a new record field), `errordef`.
- `azure2.evaluate(x)` → `pyazr.Evaluation` from one evaluation (`AZUREAPI::EvaluateObjectiveRWA`,
  `api/src/pybind_evaluation.cpp`): `objective`, `chi2` (`Total-Chi-Squared`), `terms` (the engine's term names),
  `data[i]` (the row of `chiSquared.out`, norm and shift priors, energies, angles, model, datum, residuals, `extra`
  from the binding's table of point outputs, `kinds()["data_outputs"]`), `model_data`, `log_prob`, `table()`.
- Calls of both objects take a per-session lock (without it two threads crash the process: segmentation fault).
- Examples `fit_scipy.py`, `fit_iminuit.py`, `fit_emcee.py` (rewritten on `log_prob`, a spawned pool, one session per
  worker).  User guide section "The objective and the evaluation"; `mcmc.rst`; `pyazr/CLAUDE.md` recipe.
- 293.94 against 297.86 on the fitted T2: the data χ² (`Total-Chi-Squared`, `Evaluation.chi2`) against the
  objective (the program's printed `Total Chi-Squared`, `Evaluation.objective`); the difference is the priors.

**Verified** (build `work/b-pyobj`, 0 warnings): all projects x 3 bases against `tab-t50-1` 887 / 851 / 887
identical (tag `tab-pyobj-1`); `ctest -LE slow` 244 of 245 (the failure: `pyazr_ui_routes_consistency`, below); the
slow tests 12 of 12; the consistency test; the new tests `pyazr_objective_object`, `pyazr_evaluation`,
`pyazr_fit_recipes`, each shown to fail with an injected error (list in the report).  Not run: `tests/run_tests.sh`
and the converted-tree gate (no engine code outside the API library changed), Sphinx (not installed).

**Open.**
- `pyazr_ui_routes_consistency` fails: `fit_objective` is a new public name of the session and `pyazr/ui` (not this
  branch's) neither calls it nor lists it in `AZURE2_NOT_SERVED`; `web_ui.md` section 12 lists the calls to move to.
- The program's MIGRAD on T2 stops 9.1e-4 above the minimum scipy's least_squares and L-BFGS-B find (297.8594).
- Today's session methods are not safe from two threads on one session (they take no lock).
- Not done: the dense model line (finding 21), the σ and S columns as the writer forms them, read-only views.

## 2026-11-13 (agent, branch `feature/upgrade-docs-developer-guide`) — the Developer Guide (documents only)

**Done.** `docs/source/developer/` rewritten for the code at ad6235e6: `architecture.rst` (the libraries from the
CMake link lines, a run step by step for the program and for a pyazr session), `classes.rst` (the scheme of the
classes: five figures, libraries and four layers), `state_and_ownership.rst`, `extending.rst` (each kind of extension
with its places and the test that fails), `data_model.rst` (replaces `edata_and_target_effects.rst`, which stays as an
orphan page pointing to it), `level_matrix.rst` (history, timings and source line numbers removed; under 500 lines),
`physics_core.rst` and `api.rst` (were stubs), `index.rst` in reading order; one paragraph in `whats_new.rst`;
`contributing.rst`'s keyword section points to `extending.rst`; `theory/conventions.rst` names sections of
`level_matrix.rst` instead of its line numbers. The scheme is drawn from `docs/source/developer/figures/classes.tsv`
and `arrows.tsv` by `make_figures.py` (standard library, Graphviz `dot`); the SVGs are committed. Consistency check
44: every type of `classes.tsv` is declared in the file it names, every arrow joins two boxes of its figure.

**Verified how.** `python3 tests/consistency/repo_consistency_test.py` passes; check 44 shown to fail with a renamed
row, an arrow left without its box and a class renamed in `include/Segment.h` (restored). A name check (scratch
script) finds every file path and code name of the pages in the sources but the planned names and include-path
notation. Every figure rendered with headless Chrome at 800 px and looked at. The site was not built (no Sphinx on
this machine); the markup was checked by script (underlines, list tables, inline literals, `:doc:` targets, figure
paths) and by reading.

**Open.** `coulomb_table.rst` and `ec_table.rst` keep dates, status lines and paths of the fit archive. When `CNuc`
is split or pyazr 3.0 lands, the rows naming `CNuc`, `azure2`, `ParameterSet` in `classes.tsv`/`arrows.tsv` change
(then `make_figures.py`), and the "Planned" boxes of the pages go.

## 2026-11-13 (agent, branch `feature/upgrade-9-webui-physics`) — the web UI's Physics tab

**Done, first part** (the tab on the physics routes; design in `web_ui.md` section 6.6, as built in 11.5). The tab
lists the rows of `/api/physics` and the session's functions (`/api/session/functions`, new) grouped by the row's
`group`, with a search box; each row becomes a form of the existing generator (`physicsform.js`: the table
`ARGUMENT_KINDS` maps an argument kind to a value type of `forms.js`); the axis is sent as `{from, to, points,
spacing}` and made by the server with numpy; several values of an argument are a broadcast dimension each (one line
per combination); the plot has one trace per field and line with gaps (null) where the status is not 0; the status
line counts each status with its words; CSV of the arrays as returned; the answer's `python` (made by the server from
the request) shown, copied, saved. The binding's rows gained `group` and, per argument, `example` and `sweep`
(binding code only, `api/src/pybind_physics.cpp`; no engine number); the row check in `physics_table_test.py`.

**Done, second part** (stopping). The row types moved into `api/src/pybind_rows.h` (shared by the physics and
stopping tables; `ArgKind::text` added); `pybind_stopping.cpp` has a table of five rows (`stopping_power`, `range`,
`energy_after`, `areal_density_between`, `energy_straggling`) with SRIM's citation, `_azure2.stopping.functions()`
and `evaluate()`; `pyazr.stopping.table()` and `evaluate()`; `kinds()["stopping"]` is that table (it was a signature
and a doc per binding function); `GET`/`POST /api/stopping` served by the physics route's handler; the page lists the
rows under "stopping" and draws the citation under every stopping plot, in the CSV and the Python text.

**Verified how** (build `work/b-uiphys`, 0 warnings): T1, T2, T3, 13N x 3 bases identical to `tab-t50-1` (tag
`tab-uiphys-1`, 78 files per basis); the consistency test; `physics_table_test.py`; the UI tests as scripts
(`physics_tab_test` 4, `physics_routes_test` 6, `routes_consistency_test` 8, `static_test` 9, `node_test` 46 node
tests); the smoke test on T1 and T2 with screenshots (`work/uiphys/shots*`); each new check shown failing once with
an injected error (14 cases, `tests/pyazr/ui/README.md`). Second part: T1, T2, T3, 13N x 3 bases identical to
`tab-t50-1` (`tab-uiphys-2`); 4 more injected errors. At the end: `ctest -R pyazr_ -LE slow` 125 of 125 (one at a
time); ten further smoke runs 10 of 10 (`work/uiphys/smoke10/`). Not run: all projects x 3 bases (no C++ outside
`api/src/pybind_*.cpp` and the new `api/src/pybind_rows.h` changed), the slow tests, `tests/run_tests.sh`, sphinx.

**Open.** The level matrix, the transfer functions on potentials and the reaction rate (no engine call); the THM
tab; the session's functions described by the binding (finding 24); stopping of the project's own targets (26); the
stopping provenance per element not shown; `J` of `channel_functions` not offered.

## 2026-11-13 (agent, branch `feature/upgrade-4-cnuc-finish-design`) — the end of `CNuc`: design (documents only)

**Done.** `cnuc_split.md` §12: the inventory of `CNuc` measured on the code at 6e8d000b (every member with its writers,
readers and class; the fit, input and report sets of `ALevel`; `PPair`, `JGroup`, `AChannel`; the 134 files naming the
class, grouped by what they want), the end state (`CompoundSystem` immutable and shared, holding the `System` built
once; `LevelInputs` a value; `LevelState` owned by the evaluation; `ObservedLevelTable` returned by a function;
`LevelParameterLayout` the one walk of the level block), six steps (five bit-identical, one that fixes a defect and
changes numbers), the defects in the way, the PI's measure before and after, costs, other open work, decisions. One
pointer line in `ROADMAP.md`. No code, test or build file changed.

**Measured.** The counts by scripts in `work/cnucdesign/` (`defs.py`, `members.py`, `naming.py`, `group.py`); the copy:
`clone_probe.cpp` linked against `work/b-trial50`'s libraries, one thread: `Clone` copies 7.8 / 12.9 / 131 kB in
5 / 7.5 / 56 µs (T2 / 13N / `natural_neon`), of which the per-vector part is 0.8 / 1.4 / 22 kB; one calculation of each
under `gdb` cloned 1 / 0 / 1 times. Nothing else was run (no tree comparison: no code changed).

**Found.** The notes are behind the code in five places (§12.1.5: the decay-position defect, the trees and
`ExitPairKeysInFile` are done or gone, two photon pairs of one final state are refused, the engine does call eleven
`level_parameters` functions). A defect by reading, not run: the sub-point-grid copy of `TargetSetup.cpp` keeps the
final state's N_f^½ of the moment it was made (Brune basis, external capture, fitted final-state amplitude); proposed
as its own step before the state moves. Dead: `PPair::SetBindingEnergy`, `PPair::ec_entrance_`, `ALevel::NumGammas`,
`ALevel::SetShiftDerivatives`, `CNuc::GammaDecayInputs()`.

**Open.** The decisions of §12.7 (structure classes kept in design 1 or replaced by the `System` now; names;
abundances owned by the data; a pool of `LevelState`; the defect's step; the order against the sampler, records and
physics-table branches).

## 2026-11-13 (agent, branch `feature/upgrade-8-pyazr3-cli`) — the pyazr command line

**Done.** `python -m pyazr COMMAND ... [--json]` (`pyazr/cli.py`, `pyazr/__main__.py`; `pyazr3.md` 2.8 and step 7 of
section 4). One table, `COMMANDS` (name, arguments with types, the pyazr function called, help line), makes the parser,
every `--help`, `pyazr commands --json` and the reference page `docs/source/user_guide/pyazr_command_line.rst`.
Commands: describe [FILE], validate, convert, schema [KIND], kinds [TABLE], stopping, coulomb-table check, ui,
import exfor, import ensdf, commands; coulomb-table build and fetch are rows that exit 2 with what to do instead.
`--json`: one document on standard output (`command`, `ok`, `exit`, `result`, `error` with `what`, `allowed`, `fix`),
the engine's output sent to standard error meanwhile. Exit codes those of the program's project commands (0; 1 a
failed run; 2 a wrong call or a refused input; 3 an old project). Recipe "A command of the command line" in
`pyazr/CLAUDE.md`. No C++ changed.

**How.** `tests/pyazr/cli_test.py` (ctest `pyazr_cli`, about 15 s at one thread): one case per row with `--json`, the
JSON read; describe, validate and schema equal to `AZURE2 --describe/--validate --json` and `--schema` on T2 and on a
project with errors; ten wrong calls in both modes (exit 2 or 3, the three parts, nothing on standard output without
`--json`); `--help` of every command names every argument. Consistency check 42: a row without a case, a stale page.
Each shown failing with an injected error (a case removed, validate's status dropped, a print on standard output in
`--json` mode, a flag's help suppressed, the page changed). Extensibility drill: a command `predictions` was three
places (a row and its four-line call in `cli.py`, a case in the test, the regenerated page; a scratch commit, not kept).

**Open.** `coulomb-table build` needs the engine's `CoulombTable::BuildFile` bound (and decision 7's "built on first
use" is not implemented); `fetch` needs a published copy. `pyazr.describe(path, record=NAME)` with an unknown NAME
says "the project has errors:" with no error (the program says the record does not exist). The program's
`--validate`/`--describe` report the stopping of each target, the binding's functions do not. `python -m pyazr`
needs the compiled engine even for `--help` (`pyazr/__init__.py` imports it). `ui --json` catches the URL line
`pyazr.ui.server.main` prints, for want of a function that returns it. No `import exfor --to FILE` (no single pyazr
call writes any EXFOR entry as a data set). The program's `--validate FILE` exits 1 for an unreadable file; the
command line exits 2 for a missing one (a wrong call, as the brief asks).

## 2026-11-13 (agent, branch `feature/upgrade-4-session-options`) — the options of the calculation are `<config>` keys

**Measured** (the table of `project_format3.md` §3.12): eleven options of `RuntimeOptions` and the program's flags
change what a project computes (`use_brune`, `use_park`, `transform` — conventions; `ignore_externals`,
`use_long_wavelength`, `ec_convection_current`, `use_rmc` — physics; `use_gsl_coul`, `use_coulomb_table`,
`use_ec_table`, `replace_grids` — numerics); three are not project state (`data_mode`, the Coulomb table file,
`cwd`); the fit's, output, program and command flags and the environment switches are left out. The defaults of the
program and a session are the same for all eleven; two rules of combination differ: `--use-rmc` drops Brune and Park
with a warning where a session keeps them, and the program applies `--use-park`/`--no-brune` in the given order where
a session lets Park win. Not unified (decision for the PI; recommendation in §3.12).

**Done.** The options are bool keys of `<config>` under the pyazr keyword names; a flag or keyword given wins over the
key, the key over the default. One table: the option rows of `io::ConfigKeys()` (flags, default, effect) with the
implications and contradictions beside it; `src/ProjectOptions.cpp` keeps where `Config` stores each and the sources;
the program's parser and `--help` read the rows; the binding's option table is derived (plus the three non-project
options with their reasons). `--validate`, the run and `Project.set_config` refuse a bad value
(`E_CONFIG_OPTION_VALUE`) and a contradiction with a key (`E_CONFIG_OPTIONS_CONTRADICT`, both named with their
sources). `--describe` (text and JSON) and `azure2.describe()["option_sources"]` say each option's source. pyazr's
signature is unchanged: a decorator records the keywords a call named (`pyazr/session_options.py`). The web UI's
`SESSION_OPTIONS` is gone: the options shown are the project's with a `source` field; edited through the config route.

**Verified** (build `work/b-options`, 2 threads, 0 warnings). On the code commit: all projects × 3 bases against
`tab-t50-1` 887 / 851 / 887 identical (`tab-options-1`); the eight fits identical to `b-trial50`; the converted-tree
gate 679 / 643 / 679 identical; `ctest -LE slow` 236 of 237 (`pyazr_validate`: the describe texts, fixed in the next
commit); the slow tests 11 of 11; `tests/run_tests.sh` 165 passed, 2 failed (`ratio_to_rutherford`,
`units_scales`: their committed `--describe` text, fixed in the last commit by printing the options only when one is
not at its default). The later commits change only the description text: the failed tests, `pyazr_describe`,
`pyazr_config_options`, `pyazr_ui_api`, `format3_engine`, `describe_relativistic`, `cli_help` rerun and pass.
`tests/config_options/check.sh` (key = flag, 22 cases) and `pyazr_config_options` (key = keyword) pass; each new
check failed once with an injected error.

**Open.** The two combination rules above (PI). `save_fit` on T1 in the formal basis fails its own verification and
`use_rmc` on T1 fails a session's set-up but not the program's run (both as before; not of this branch).

## 2026-11-13 (agent, branch `feature/upgrade-6-sampler-objective`) — the sampler on the fit's objective, its state

**Measured first** (`docs/upgrade/objective.md` section 6, with the table term by term and line numbers): at the
branch start, every term of the fit was in the sampler's posterior up to the Gaussian priors' normalisation (the
difference sampler − (−½ χ²) constant at three points: T2 0.881047, 13N_fast 2.076794, `nuisance_prior` 3.686232,
0 for the record without priors, 13N, `energy_shift_sqrt`, the THM toy, `identical_pp_res` under Park), in both
bases, except (1) the `gamma_decay` record's priors and wall (difference 3.06 / 0.90 / 81.1 at the three points of
the record with priors; 20161.8 with Σb > 1) and (2) the limits of every parameter (never read by the sampler:
b = −0.1, Γ_γ = −0.01, an abundance of −0.001 all finite).

**Done.** The rule: ln P = −½ × the objective from the same table of terms + C (the Gaussian priors'
normalisation), −∞ outside the support. Each `ObjectiveTerm` says how the sampler takes it (`AddForSampler`, default
its whole value, so a new term reaches the sampler with no edit of `AZURECalcMCMC.cpp`; `ExcludesForSampler` for
Park); `Objective::SamplerChi2` replaces the hand-written sum (bit for bit the same where the posterior did not
change). The limits: −∞ outside the MINUIT limits of a kind whose row says `samplerSupport` (the record's numbers
only). The state: the 17 process-wide variables of `AZURECalcMCMC.cpp` are gone: the run's file, steps and best
walker are a record owned by `RunMCMCSampling`; the callbacks are `std::function` members, the stop request an
atomic member (`RequestStop`); `g_param_names` was a copy of `paramNames_`; the two write-only counters are dropped.
numcmc takes `std::function` callbacks (source compatible). The CLI's two progress globals in `AZURE2.cpp` became
captures of a lambda. 18 rows of `allowed_state.tsv` removed.

**Verified how** (build `work/b-sampler`, 0 warnings): the new ctest `sampler_objective` (identity at three points in
Brune/Park × both sampling bases, −∞ outside the limits, a chain pulled below 0 stays inside, the record's priors move
ln P by −½ × their term, T2's identity with its constant, two samplers in two threads = two runs alone), each check
shown to fail with one injected error (limits off; T11 not taken by the sampler; the run record made static); the MCMC
check scripts (`mcmc_classification`, `mcmc_park_start`, `mcmc_physical_norm`, `mcmc_snapshot_components`,
`nuisance_prior`, `prior_covariance`, `thm_park`, `clone_composite`, `cli_eof`) pass; chains of 8 walkers × 100 steps
before (`b-trial45`, and again `b-trial50` after the merge) and after: identical `samples.mcmc`, `walkers.mcmc`,
`param.mcmc` for the record without priors (Brune RWA and observed, Park RWA), 13N_fast (both bases),
`energy_shift_sqrt`, T2; changed for the record with priors (mean ln P over steps 50–99 −25.72 → −29.43, Γ_γ 0.038374 →
0.038389 eV): the record's priors, as intended. On the merged state (the data branch at 3298f372 merged in; build 0
warnings): all projects × 3 bases against `tab-t50-1` 887 / 851 / 887 identical; `ctest -LE slow` 236 of 236; the slow
tests 11 of 11; `tests/run_tests.sh` 166 passed, 0 failed; the converted-tree gate 679 / 643 / 679 identical; T1, T2,
T3, 13N × 3 bases at 1 against 2 threads identical; `tests/gamma_decay/check.sh` 63 ok; the MCMC check scripts and
`prior_centre` pass; the eight fits against `b-trial50`: 8 of 8 identical. The scratch measurement re-run on the merged
engine: sampler − (−½ χ²) equals the documented constant `LogPriorNormalisation()` at every point (0 for the record with
priors, 0.881047 for T2), −∞ outside the record's limits.

**Open, for the director.** The limits of the other kinds (abundances, thickness, densities, Doppler temperatures,
apparatus parameters, reduced widths, `<parameterSettings>` limits) are not applied by the sampler; doing it is one
column of `ParamKindTable()` but changes those posteriors (objective.md section 6). Outside a limit the fit reads a
clamped density for its prior, the sampler the raw value (0.496 in T2's −½ χ²). The `loglikelihood` column now holds
the record's priors (it is `logprob − logprior`; documented in the user guide).

## 2026-11-13 (agent, branch `feature/upgrade-8-pyazr3-physics`) — pyazr: the physics table and `pyazr.physics`

**Done** (no engine number changed; the Python formulas untouched).
- `api/src/pybind_physics.cpp`: one table of physics functions, 24 rows (name, arguments with unit, kind and domain,
  returns with units, doc, library function, plot axis, adapter); one template evaluates a row element-wise on
  arrays with a status per element (0; k > 0 argument k outside its domain; -2 the arguments together; -1 not
  converged); `_azure2.physics_functions()` (the rows, strict JSON) and `_azure2.physics_evaluate`.  The seven
  functions bound one by one before (`sfactor_conversion`, `modified_sfactor_conversion`, `statistical_factor`,
  `resonance_strength`, `branching_ratio`, `rutherford_cross_section`, `mott_cross_section`; `whittaker` is the
  row `whittaker_wave`) are rows now, with the THM plane-wave vertex `thm_vertex_function`; their callers
  (`nds.sfactor_to_cross_section`, four tests) follow.
- `physics/Widths.h` (`wigner_limit`, `single_particle_width`, `photon_penetrability`) with `physics_widths`; the
  engine's `AChannel::SetWignerLimit` and the photon branch of the level parameters call the first and the last.
- `pyazr/physics.py`: one function per row, built at import (signature, numpydoc docstring with units and domains,
  broadcasting, scalars stay scalars, `DomainError` naming argument, value, element and domain; `errors="mask"`
  for NaN and the status array; `status_text`, `functions`, `table()`).  `pyazr.kinds()["physics"]` = the rows,
  `kinds()["stopping"]` the stopping functions (they were mixed into "physics").
- `pyazr/session_physics.py`: `channel_functions(pair, e_cm, l, radius, J)` mixed into the session class (two
  edited lines of `azure2.py`); `coulomb_functions` stays as the old name.  The rename step should move the method
  into `session.py` (or keep the mixin as a base of `Session`) and drop `coulomb_functions`.
- `tools/gen_physics_docs.py` writes `docs/source/user_guide/physics_functions.rst` from the rows (the test checks
  it is current); consistency check 41: every row is named by a `row_test` (shown to fail without one).
- Test `pyazr_physics_table`: every row against mpmath, the engine's own use on 13N (P, S, phi bit for bit with the
  session's exact Coulomb routine, `Parameter.wigner_limit` bit for bit) and the Python formula it replaces.

**Measured first: the functions of `pyazr3.md` 2.7, before this branch.**

| function | library function | bound before | domain | outside it before |
|---|---|---|---|---|
| sommerfeld_eta, wave_number | `ChannelFunctions.h` `eta`, `rho` | no | E > 0 | inf at 0, NaN below |
| lab/cm energies | `experiment/Frames.h` `BeamToCM`, `BeamLabFromCM` | no | E >= 0 | a number of the formula |
| lab_to_cm angle, J | `Frames.h` `ElasticLabToCM`, `ReactionLabToCM` | no | forward branch | false + reason |
| coulomb_waves | `coulomb/Coulomb.h` `coulomb` | through a session, at the radius | l 0-10000, rho > 0 | status, NaN |
| coulomb_phase | `Coulomb.h` `coulomb_phase` | no | l >= 0 | NaN |
| whittaker_wave | `coulomb/Whittaker.h` `whittaker_scaled` | `_azure2.whittaker` | eta_b >= 0, rho > 0 | status, NaN |
| P, S, phi, dS/dE | `ChannelFunctions.h` `penetrability`, `shift`, ... | session only | E > 0 | session: zeros |
| wigner_limit | `AChannel::SetWignerLimit` (engine member) | per channel | mu, a > 0 | 0 |
| single_particle_width | none (Python: `wigner_widths` probe) | no | E > 0 | -- |
| weisskopf_width | none | no | -- | -- |
| resonance, Coulomb, S factor (7) | `ResonanceQuantities.h`, ... | one by one | E > 0 | inf, NaN |
| THM vertex W_l(p) | `transfer/ThmVertex.h` `PlaneWaveSurfaceVertex` | no | p, a > 0 | NaN |
| bound_state, distorted_wave, momentum distr. | `transfer/BoundState.h`, `DistortedWaves.h` | no | potentials | -- |
| global optical sets | `transfer/GlobalOptical.h` | `thm_global_opticals` (list) | -- | -- |

The 17 Python formulas and their rows: `pyazr3.md` section 1.3 (the largest difference per formula).

**Added after review.** The web UI's two physics routes served (`pyazr/ui/physics.py`, one server method, one adapter
method, two rows of `routes.py`): `POST /api/physics/{fn}` in the server's process (no session; the C++ loop releases
the GIL), `POST /api/session/channel_functions` in the session worker; both at most 20000 elements per request
(`physics.MAX_ELEMENTS`). `GET /api/physics` is `kinds()["physics"]`. A refusal of the server raised in a worker
(`Refused`) is now 422 (`errors.STATUS`). `coulomb_functions` is listed in `engine.AZURE2_NOT_SERVED` as the older name.
Test `pyazr_ui_physics_routes` (registered above the block that must stay last).
Verified on the merge with d1b193db (no C++ changed after 94eb4828): build 0 warnings; the consistency test; ctest -R
pyazr_ui -j1 8 of 8; ctest -LE slow 237 of 238, the failure pyazr_ui_smoke on T2 with the two known "Failed to fetch"
requests of the Data tab, passing on three reruns of it alone; T1, T2, T3, 13N x 3 bases identical to tab-t50-1; the
routes-consistency test shown to fail with a missing row naming an existing call.

**Verified** (after merging the data branch at 3298f372): build 0 warnings; all projects × 3 bases identical to
`tab-t50-1` (887 / 851 / 887, tag `tab-physics-2`); the eight fits identical to `trial50`; `ctest -LE slow` 236 of
237 (the failure: `pyazr_ui_routes_consistency`, below); slow 11 of 11; `tests/run_tests.sh` 166 passed, 0 failed.
- `pyazr_ui_routes_consistency` failed by its own design once the calls existed; it passes since the two routes are
  served (below).

**Open.**
- Not bound: `level_matrix` (the engine hands out no R or U matrix without new code in `AZUREAPI`), the transfer
  functions that take potentials and return functions of r (bound states, distorted waves, momentum distributions),
  the Weisskopf width (no engine function: new physics, a decision for the director), the reaction rate and the
  covariance band (session quantities).
- Found, not changed (library accuracy, `pyazr3.md` 1.6): S good to about 3e-9 absolute; the engine's dS/dE 9.5e-6
  relative at p + 12C, l = 0, 3.42 MeV, against the 1e-8 `ChannelFunctions.h` states.

## 2026-11-13 (agent, branch `feature/upgrade-8-pyazr3-records`) — pyazr: parameters and pairs as the engine's records

**Done** (no engine number changed; the binding and pyazr only, and one comment of `include/ParameterKind.h`).
- `api/src/pybind_records.cpp`: one table of the fields of a parameter record (50 rows: name, unit, the kinds a field
  describes, whether `Parameters.table()` shows it as the owner, doc, the reader of the engine object) and of a pair
  record (22 rows).  `Session.parameter_records()`, `Session.pair_records()`, `_azure2.parameter_record_fields()`,
  `pair_record_fields()`; `parameter_kinds()`'s `owner` from the table, its `code` gone; `pyazr.kinds()` gains
  `parameter_fields` and `pair_fields`.  `AZUREAPI::GetParameterBounds` (limits as created).
- Deleted: `AZUREAPI::GetParameterInfo`, `GetPairsInfo`, `kParamInfoFields`, `kPairInfoFields` (301 lines of
  `AZUREAPI.h/.cpp`), the bindings `parameter_info`, `pairs_info`, the two `_NFIELDS`, the numeric type codes, the
  fall-backs for engines with 19 or 20 pair fields, `Parameter.from_record`, `Pair.from_record`.
- `pyazr/parameters.py` (804 → 620 lines): `Parameter` and `Pair` are the records; `Parameters` (the table of the fit
  vector: `x0`, `names`, `kinds`, `units`, `free_mask`, `bounds`, `free`, `of_kind`, `of_level`, `by_name`, `put`,
  `take`, `read`, `write`, `table`, `to_records`); `ParameterSet` = `Parameters` and the old views.
  `pyazr.snapshot.write_values(session, project, x)`: a `Project` edited in memory (`snapshot.py`, web UI wish 8; a
  function, not a session method: a new public name of the session would need a line in `pyazr/ui/engine.py`).
- Tests: `tests/pyazr/parameter_records_test.py` (new); the tests of the removed layouts read the records
  (`parameter_kinds`, `extrap_layout`, `thm_profiled`, `light_excitation`).

**Measured before the change** (start commit da9a3b20).  Readers of the two layouts:

| reader | fields read |
|---|---|
| `pyazr/azure2.py:421` `_build_pairs` | all 21 (`Pair.from_record`; the 19- and 20-field fall-backs) |
| `pyazr/azure2.py:473` `_build_parameters` | all 16 (`Parameter.from_record`, `parameters.py:185`) |
| `pyazr/azure2.py:2253` `_mode_free_vector`, `:2342` `_sync_layout` | 0 type, 11 fixed |
| `pyazr/bands.py:168` `live_parameters` | 0, 1, 4, 6 |
| `tests/pyazr/extrap_layout_test.py:61`, `thm_profiled_test.py:123` | the row count |
| `tests/pyazr/parameter_kinds_test.py:116` | 0, 6, 13 against the records |
| `tests/pyazr/light_excitation_test.py:111` | `Pair.from_record` of 19, 20, 21 fields |

`snapshot.py`, `modelavg.py`, `scheme.py`, `described.py` read `Parameter` attributes, not the layout.  No reader
outside `pyazr/`, `api/`, `tests/` (the skill `azure2-fit` quoted an error message of the layout's reader).

| field of the 16 | meaning | in the records before | record field now |
|---|---|---|---|
| 0 type | kind code | `code`, `kind` | `kind` |
| 1 jgroup (abundance: Z) | J group | yes (`nucleus_Z`) | `jgroup`, `nucleus_Z` |
| 2 J, 3 parity | the group's J^pi | no | `J`, `parity` |
| 4 level (abundance: A) | level in its group | yes | `level`, `nucleus_A` |
| 5 level energy; prior centre | (screening, density, Doppler, thickness) | no | `level_energy`, `prior_centre` |
| 6 channel; sample; j; effect | (and transmission/yield) | yes but the last | `channel`, `sample_index`, `fn_term`, |
| | | | `effect`, `sample_kind` |
| 7 L; thickness leader | | no | `L`, `thickness_leader` |
| 8 S; layout maximum | | no | `S`, `layout_max` |
| 9 pair (cbkg: exit pair) | | no | `pair` |
| 10 radiation type | | no | `radiation_type` |
| 11 fixed, 12 value | | no | `fixed`, `value` (and `start`, `free_index`, `bounds`) |
| 13 segment key | | yes | `segment_key` |
| 14 Wigner limit; prior width | | no | `wigner_limit`, `prior_sigma` |
| 15 input is an amplitude | | no | `width_is_amplitude` |

The 21 pair fields (number, key, type, entrance, J, parity, Z, M of each particle, separation energy, excitation,
radius, spin factor, nucleus Z and A, abundance, light excitation, THM binding energy) were in no record; all are
fields of the pair record now, plus `pair_kind`.  Read by nobody outside `parameters.py` itself: `input_is_rwa`,
`sample_kind`, `thickness_leader`, the apparatus value of field 5, the pair's `i1i2factor` (repr only).

**Verified** on the merged state (build `work/b-records`, 0 warnings; 2 threads): all projects in three bases
against `tab-t50-1` (tag `tab-records-2`): 887 / 851 / 887 identical; the converted-tree gate (`grecords2`): 679 /
643 / 679 identical; the eight fits of `work/director/fits_vs_base2.sh` against `work/b-trial50`: 8 identical;
`ctest -LE slow`: 235 of 236, then `ui_routes_consistency` (a new public name of the session) and the slow
`formal_brune_iteration` (it finds a parameter by its repr) after two Python fixes, with the other tests concerned:
10 of 10; the slow tests 10 of 11 before that fix; `tests/run_tests.sh`: 166 passed, 0 failed;
`repo_consistency_test.py`: consistent.  Before the merge: the same tree and gate against `tab-t45-1` identical, the
eight fits against `b-trial46` identical.  Injected errors: a kind without an owner field (C++) and `start` read from
the physical values (C++) each made `parameter_records_test` fail.  Speed on 13N (before / after): session
construction 0.59 / 0.62 s (noise), the first `parameters` 1.5 / 0.6 ms, cached accesses 0.1 µs each, `pairs`
0.2 / 0.03 ms.

**Merged** `feature/upgrade-4-data` at 3298f372 (prior centres, diagnostics record, web UI): `PriorOf` reads the
centres where the objective reads them (`include/Objective.h`); the conflicts were the deleted `GetParameterInfo` /
`GetPairsInfo` and `Parameter.from_record` (whose new `prior_nominal` of an abundance and an apparatus number is the
record's `prior_centre`).  The web UI (`pyazr/ui/engine.py`) reads `describe()["parameters"]` and
`parameters` (`free_index`, `fixed`, `name`), both from the records now; its text of a fit (`save_fit` to a file,
read back) is what `pyazr.snapshot.write_values` gives in memory.

**Open.** (1) `bounds` are the limits as created; a `<parameterSettings>` block's limits act only in the
program's fit.  (3) `_sync_layout` and `_mode_free_vector` count energies and widths as the free R-matrix block, not
the `radiative` numbers in a level's slot (unchanged behaviour; a project with a gamma_decay record and
extrapolations may carry a wrong vector over a mode switch).  (4) Web UI wish 11: the pair record carries the THM
binding energy, the width record `width_is_amplitude`; `Project.set_pair` / `set_channel` still cannot edit them.

## 2026-11-13 (agent, branch `feature/upgrade-4-layout-centre`) — grids and limits laid out from the prior's centre

**Done.** Measured first (base build `work/b-prior`; fit, `param.sav` read by pyazr, `save_fit`, both sessions at the
same vector; `parameter_kinds.md` section 4 has the table): the Doppler temperature's χ² moved by 5.5e-6 (6.0e-6 and
5.6e-6 per data set, 4.85481 → 4.85478 in `chiSquared.out`, above `save_fit`'s own 1e-6, which refused its snapshot),
the stopping scale's by 5.8e-7 on one data set; the areal density's and the kernel width's χ² reproduced, but every
layout and lower limit moved with the start (density 1.1e18 → 1.19e18, 363 → 379.6 K, scale 1.12 → 1.126, width
23 → 29.0 keV); abundance and thickness (no grid) reproduced exactly. Rule (design committed first, section 4 "The
layout of a fitted value"): one function, `ChooseLayoutBase` (`experiment/.../LayoutBase.h`), with the three spans as
rows; the base of the layout, the lower limit and the Doppler pieces is the written centre when the parameter is fitted
and the start lies in the centre's range (`grid_max` as written), else the start with one `Note:` line. Used in
`TargetEffect` (density, Doppler, stopping scale, the kernel's parameters: the one place, the declaration reads it);
`ApparatusLayout` removed; the fields that were the base renamed `layoutBase`. After: the four written projects at the
fit's vector give the fit's objective and every data set's χ² exactly (difference 0), the same layouts, and
`save_fit(verify=True)` passes for the Doppler case.
The web page's "accept into project" (`pyazr/ui/server.py` `Server.apply` → `engine.py` `Engine.snapshot`) writes
through the same `save_fit(path, x, param_sav=False)` (verify on) that `continued.py` tests; its text merge
(`state.merge_snapshot`) is not covered here. The data branch (3298f372) merged in.

**Verified how** (build `work/b-layout`, 0 warnings), on the merged head: all projects x 3 bases against
`tab-t50-1` 887 / 851 / 887 identical (also against `tab-t45-1` before the merge); the eight fits of
`fits_vs_base2.sh` identical to `b-trial50`; the converted-tree gate 679 / 643 / 679 identical; `ctest -LE slow`
237/237; the slow tests 11/11; `tests/run_tests.sh` 166 passed, 0 failed; `tests/nuisance_prior_centre/check.sh`
with and without `PYAZR_ENGINE_DIR` (new check (h) = `continued.sh`). Injected errors (inclusive lower end; the
kernel's layout from the start) fail the new checks; the base build fails the Doppler χ², the note, `save_fit` and
the four layout checks. Not run: 1 against several threads, Sphinx.

**Open.** A parameter without a prior has no centre: a continued fit lays its grid out from the fitted value (1.5 ×);
whether `save_fit` should write `grid_max` there is the PI's call. `save_fit` writes a stopping scale shared by two
targets on the first only, and the written project is refused (`E_STOPPING_SCALE_CONFLICT`; pyazr write-back, not
fixed). The Doppler pieces from the start did not change χ² in the test case (piece counts equal), so only the window
and layout are covered there; the fit limits are not reported by the API, so the lower limits are checked by reading.

## 2026-11-13 (director, 16:49) — merged: the web UI's Results tab

**Merged** `feature/upgrade-9-webui-results` (the entry below; it carries the Fit tab). No C++ changed, so the
reference tag and the base build stay `tab-t50-1` and `work/b-trial50`.

**Checked by the director** on the branch's last commit with the verified build's engine, one thread: the
consistency test; the UI tests api, jobs, server, routes_consistency, static, node, results and the browser smoke
test, all ok; five further smoke runs, no failure.

**The smoke test's failed fetches** (the first open item of the entry of 16:16) are closed as a property of the
test's Chrome, not of the page: in each run's fresh profile Chrome loads its component extensions a few seconds
after start and cancels requests of the page that have not left the renderer (the server never saw them; 4 runs
in 10 on the page before either tab). The test's Chrome now starts without extensions, component updates and
background requests; no retry, no longer time budget. Not shown: which extension, and that a browser in daily use
never meets it.

**Open.**
- Nobody has used the page in a real browser by hand (the export buttons, the output-file buttons, the keyboard).
- The model line of a data set joins the engine's values at the data points; a dense line needs pyazr's
  `Session.evaluate` on a grid (finding 21 of `web_ui.md`); the lab energy of each point and the program's σ and S
  columns are findings 22 and 23.
- T3's extrapolation at small angles reaches 1e63 on a log axis (the engine's number, not looked into).

## 2026-11-13 (agent, branch `feature/upgrade-9-webui-results`) — the web UI's Results tab

**Done.** The director's corrections to the page (bounded file preview, wrapping tables, repeated groups as small
tables with number boxes as wide as their text, the full χ² in the status line, paths broken at `/`, the UI tests
registered above the last block of `tests/pyazr/CMakeLists.txt`). The Results tab: the payload
`GET /api/evaluation` (the program's χ² split, every data set's and extrapolation's arrays, residuals, S factor,
model data), `POST`/`GET /api/outputs`, `tabs/results.js`, plot rows in `plots.js` (residual axis, level markers,
log axes, CSV with 17 digits, SVG/PNG). `web_ui.md` 11.4 and findings 21–23; user guide; `pyazr/CLAUDE.md`.

**How.** `ctest -R pyazr_ui` equivalents run as scripts at one thread (api, jobs, server, routes, static, node,
results, smoke: all ok); the smoke test's "Failed to fetch" (1 in 10 runs on T2, 4 in 10 on the page's merge commit)
came from Chrome loading component extensions into each run's fresh profile, which cancels requests not yet sent;
the test's Chrome now runs without extensions, component updates and background requests: 0 failures in 20 runs
(`tests/pyazr/ui/README.md`);
`results_test.py` compares the payload with the program's `AZUREOut_*`, `.extrap` and `chiSquared.out` on T1 and
T2; each new check shown failing with an injected error (`tests/pyazr/ui/README.md`); screenshots in both themes
and sizes looked at (`work/uiresults/shots/`). No C++ changed.

**Open.** The model of a data set is drawn at its points only; lab energy, Legendre orders, the program's σ of an
S-factor data set, rate and bands wait for pyazr (findings 6, 14, 21–23). Nobody has used the page in a real browser
by hand.

## 2026-11-13 (director, 16:16) — merged: the web UI's Fit tab

**Merged** `feature/upgrade-9-webui-fit` (the entry below). No C++ changed (the diff against the data branch in
`src include api physics experiment stopping support` is empty), so the reference tag and the base build stay
`tab-t50-1` and `work/b-trial50`.

**Checked by the director** on the branch's last commit with the verified build's engine, one thread: the
consistency test; the UI tests api, jobs, server, routes_consistency, static, node and the browser smoke test, all
ok; the screenshots of the idle, running, finished, MINOS and accepted states looked at.

**Open.**
- The browser smoke test fails now and then on T2 with two requests of the Data tab that never reach the server
  ("Failed to fetch"): one run in two for the director on the Results branch, about one in eight for the Fit
  agent. It is older than both tabs. The Results agent is finding the cause; no retry may hide it.
- After Accept the page shows its note about the two χ² also when they differ in the last digit only (T1:
  462.27351124742853 against 462.2735112474378). The note should appear only above rounding (a relative difference
  of about 1e-12), and Accept should be disabled once the job's point is the project's. With the second part of the
  Fit tab.
- T2 from the page: iminuit's MIGRAD needs 1.7 s per call on one thread and HESSE turns its valid minimum invalid
  (objective 298.9167 after 416 calls). The engine's own MIGRAD is not offered yet (pyazr's `fit_with_engine` with
  progress and stop); compare the two on T2 when it is.
- A fitted areal density is laid out from its start value, so T2's χ² as written differs from the fit's by 2e-7
  relative: `feature/upgrade-4-layout-centre`.
- `GET /api/parameters`, `GET /api/x` and `PATCH /api/x` now answer while a job holds the thread budget (the header's
  χ² then waits for the job's end and says so).
- Nobody has used the page in a real browser by hand.

## 2026-11-13 (agent, branch `feature/upgrade-9-webui-fit`) — the web UI's Fit tab: parameters and recipe-fit jobs

**Done.** The Fit tab (`pyazr/ui/static/tabs/fit.js`, `tabs/fitjob.js`, `fitting.js`, `selftestfit.js`; a row of
`PLOTS`): the job form drawn from the server's job table (settings with types and defaults; the rows that wait for
pyazr greyed with their reason), a running job followed by a long poll (calls, objective, best, time, the objective
against call on a log scale, the log's last lines, Stop, Kill; found again after a reload), the result (the recipe's
words, χ² and N at the best point, Accept into project with Undo, the χ² of the project as written beside the fit's,
MINOS offered after an iminuit fit), the table of every parameter (kind filter from the engine's kinds, free only,
text; the working x edited in place; the job's best point, errors, MINOS intervals, start point), covariance and
correlation as tables and CSV. Server (`pyazr/ui` only): settings as typed rows (`Setting`, `SETTING_TYPES`),
iminuit's MINOS as the job kind `minos/iminuit`, covariance and correlation from HESSE, the job's trace and log tail,
`GET /api/jobs/{id}?wait=s`, a job's revision named `text_rev`, the accepted revision kept in the job's record, the
parameter table and x readable while a job holds the budget. Docs: `web_ui.md` 11.3 and findings 17–20, user guide,
`pyazr/CLAUDE.md`, `FOR_THE_OWNER.md`, `ROADMAP.md`, the tests' README. No C++ changed.

**How.** The seven UI test files at one thread on the base build's binding (`work/b-trial50`), all ok; each new check
shown failing with an injected error (`tests/pyazr/ui/README.md`); the smoke test (headless Chrome) on T1 and T2 with
the new walk through the tab (Run, progress, Stop, Accept, Undo: the header's χ² back exactly). A walk over HTTP on T2
against a server process at one thread (`work/uifit/walk.py`): 62 parameters, 19 free; χ² 10905.2313020616 (N 330) at
the start; an iminuit fit stopped after 7 calls (the long poll answered 0.88 s after the stop request; best objective
7087.924126); started again from its best point, it ended after 416 calls (779 s): objective 298.9166933, χ²
295.0675134652589 at the best point; accepted, the project as written gives χ² 295.06757397066707 (2 × 10⁻⁷ more; likely
because T2 fits a target density, the known open point; on T1, which does not, the two differ by 3 × 10⁻¹⁵); Undo gives
10905.2313020616 again. A job GET costs 2.5 ms (2.1 kB). Screenshots in `work/uifit/walk/`. A build of the branch with 0
warnings and `ctest -R pyazr_ui`: see the report.

**Open.** The engine's fitters, the sampler and bands wait for pyazr (`fit_with_engine`, `Objective`); a fixed
parameter is edited in its record or level, not from the table (finding 17); the physical value's unit is not sent
(finding 18). No event stream: a job is long-polled, the revision polled. Nobody has used the Fit tab in a real
browser by hand. A smoke-test run on T2 once saw two Data-tab requests fail in Chrome before reaching the server
(one run in about eight; not reproduced since).

## 2026-11-13 (director, 14:20) — merged: diagnostics record, prior centres, the web UI's server and page

**Merged** into the data branch after one verification of the four branches together (with pyazr's errors, merged
earlier): `feature/upgrade-6-evaluation-4` (the diagnostics record of an evaluation; the environment switches read
once per process by one reader, listed in one table and in the manual; `--describe` prints the switches that are
set), `feature/upgrade-4-prior-centre` (the centre of every nuisance prior honoured: areal density, Doppler
temperature, stopping scale, screening, thickness, abundance, apparatus; `save_fit` writes the fitted target values
back; a rebuild of the model drops the nuisance manager of the old model), `feature/upgrade-9-webui-server`
(`pyazr/ui`: server, state with history and revision, workers, jobs, the route table) and
`feature/upgrade-9-webui-page` (the page: header, tabs Project, Nucleus, Data, generated forms, preview plot;
plotly.js 4.1.2 vendored under MIT; no engine change).

**Verified** on the four engine branches together (`work/b-trial50`, logs `work/director/vt50/`): all projects in
three bases identical to the reference (887 / 851 / 887 files), 1 against 3 threads identical; the eight fits
identical to the base build; the switch spot check; `tests/gamma_decay/check.sh` at 3 threads;
`tests/nuisance_prior_centre/check.sh` without and with pyazr; `ctest -LE slow` 232 of 232; the slow tests 11 of 11;
`tests/run_tests.sh` 166 passed, 0 failed; the specification equivalence 181 equal, 0 differ; the converted-tree
gate 679 / 643 / 679 identical. The page (Python and static files only) on the merged tree: the consistency test
and its seven tests (api, jobs, server, routes, static, node, the headless-Chrome smoke test) pass.

**Reviewed by the director.** The prior-centre branch in full (no objection); the page's server edits: the file
route stays inside the root (resolved path, tested), the Content-Security-Policy allows inline styles and `data:` /
`blob:` images for plotly and nothing else beyond the server itself, no server text is inserted as HTML.

**Open.**
- A fit continued from its written-back values is computed on other grids than the fit (density, temperature,
  stopping scale and kernel parameter are laid out from the start value): branch `feature/upgrade-4-layout-centre`.
- The page's corrections (bounds of the file preview, cut text in the session-options table, the layout of repeated
  groups, the web UI's tests registered below the block that must stay last in `tests/pyazr/CMakeLists.txt`) go
  with the Results tab.
- The stopping scale's default width is relative to the start, not the centre.
- Nobody has used the page in a real browser by hand.

**Reference** for the next branches: outputs `tab-t50-1`, base build `work/b-trial50` (outside the repository).

## 2026-11-13 (agent, branch `feature/upgrade-9-webui-page`) — the web UI's page: Project, Nucleus, Data

**Done.** `pyazr/ui/static/`: the header (project, unsaved, revision, Undo, Redo, Save, the engine's χ², status
line), the tabs Project, Nucleus and Data, placeholders for Fit, Results, Physics, THM; forms generated from
`pyazr.schema()` (one widget per value type in `forms.js`, layout hints in `layout.js` only); the engine's description
beside each record; a preview plot of data and model with plotly.js 4.1.2 vendored (MIT, hash checked), a table
without it, CSV export. Server edits: static rows and handler, the token-to-cookie redirect, a CSP header,
`GET /api/files/{path}`. Twenty Qt items marked dropped in `qt_parity.md` (the PI, 2026-11-13). Docs: `web_ui.md`
11.2 and findings 13–16, user guide "The web interface", `pyazr/CLAUDE.md`, `FOR_THE_OWNER.md`, `ROADMAP.md`.

**How.** `ctest -R pyazr_ui` (7 tests, node on PATH: all pass); the smoke test on T1, T2, `tests/gamma_decay/
models/record.azr`, `examples/li7_tumino2006` (THM) and `tests/reaction_rate/models/moved_azr_truth.azr` (no data
sets: the page shows the engine's refusal), screenshots in both themes and sizes; a walk through the meeting's list
against a server (χ² after a width edit equal to a direct pyazr session's bit for bit; Undo restores it exactly).
No C++ changed; build 0 warnings.

**Open.** Nobody has used the page in a real browser by hand. Session options read-only until they are `<config>`
keys. No event stream (polling). Pairs add/remove, channel units, Legendre orders, a `<config>` schema wait for
pyazr (findings 10, 13–15). The static files are not yet package data of the wheel (step 6).

## 2026-11-13 (agent, branch `feature/upgrade-9-webui-server`) — the web UI server on pyazr's errors and description

**Done.** Merged `feature/upgrade-8-pyazr3-errors` into the server branch and adapted `pyazr/ui`: the status of an
exception is the first row of one ordered table naming its class or any base (`pyazr/ui/errors.py`; every class of
`pyazr.errors` maps to 404, 409 or 422 by it), and a worker sends the names of the class and its bases; the
engine's notices (`EngineWarning`) of a request in the session worker come back in the answer's `warnings`;
`/api/kinds` and `/api/physics` are served by `pyazr.kinds()`; `/api/describe` carries `pyazr.describe()` and the
session options (the engine's table with the value in use); the adapter's parameter rows are `azure2.describe()`'s;
a request whose `Origin` is not the server's own is refused, and a Host `[::1]` is parsed. `web_ui.md` 11.1 and 12
say what the errors branch delivered and what is still missing.

**How.** Build `work/b-uiserver` at two jobs (0 warnings); the four UI test files at one thread (api 18, jobs 5,
server 2, consistency 8: all pass); `errors_test.py` and `describe_test.py` of the merged branch (all ok); each new
rule shown failing with an injected error (`tests/pyazr/ui/README.md`); `ctest -R pyazr_ui`.

**Open.** The three parts of a message are one text, not attributes: `allowed` and `fix` of the JSON body stay empty
outside the engine's diagnostics. A job's notices go to its log, not to the client.

## 2026-11-13 (director, 13:50) — merged: pyazr's errors, engine notices and self-description (first step of pyazr 3.0)

**Merged** `feature/upgrade-8-pyazr3-errors` (entry below). Verified by the director on a rehearsal of this merge
(`work/b-trial46`): all projects × 3 bases identical to the reference (887 / 851 / 887), the converted-tree gate
identical (679 / 643 / 679), T1/T2/T3/13N at 1 against 3 threads identical, `ctest -LE slow` 226/226, slow 11/11,
`tests/run_tests.sh` 165 passed and 0 failed, spec equivalence without a difference.

**Behaviour a user of pyazr meets.** A parameter vector with NaN or infinity raises `ParameterVectorError` (it was
computed before); `coulomb_functions` outside its domain raises `DomainError`; `pyazr.__version__` is `3.0.0.dev1`.

**Next on pyazr.** Parameter records in place of the two numeric layouts (second step), the physics table (fourth),
the command line (seventh); the briefs are written. `pyazr/azure2.py` is one line below its size limit: the
rename step splits it.

## 2026-11-13 (agent, branch `feature/upgrade-8-pyazr3-errors`) — pyazr: errors, engine notices, self-description

**Done** (additions only: no public name renamed or removed, no engine number changed).
- `pyazr/errors.py`: `Error` and its classes, each also the builtin the call raised before (`FormatError` moved here
  from `project.py`; `ParameterError`, `ParameterKeyError`, `ParameterIndexError`, `ParameterVectorError`;
  `AZURE2Error` = `EngineError` = `_azure2.AZURE2Error`, one class that keeps the name `AZURE2Error`;
  `DomainError`, `DomainZeroDivision`, `SessionClosed`, `MissingFile`; `EngineWarning`).  The binding raises the
  Python classes and refuses a NaN or infinity in a parameter vector, naming the parameter
  (`api/src/pybind_api.cpp`, `BindErrors`, `expect_finite`).  `pyazr/session_checks.py`: the session's checks.
- Engine notices: the session's `Config::outStream` passes through `api/src/pybind_notices.cpp`; four notices (EC
  amplitude file written / appended, normalization and shift files) become `EngineWarning`s with file and consequence.
- `pyazr.describe()` (no path), `pyazr.kinds()`, `azure2.describe()` (`pyazr/described.py`);
  `api/src/pybind_describe.cpp`; `_azure2.physics_functions`; package docstring with a tested session; stale
  docstrings corrected; user guide sections "Asking the package about itself" and "Errors"; `pyazr/CLAUDE.md`.
- Tests `tests/pyazr/errors_test.py` (27 wrong calls), `tests/pyazr/describe_test.py`.

**Verified** (build `work/b-pyerr`, 2 threads).  Build: 0 warnings.  `repo_consistency_test.py`: consistent.
`JOBS=1 scripts/check.sh --quick`: 225 of 225 tests without the slow label passed (111 pyazr among them).  The slow
pyazr tests: 9 of 9 at two threads, then `parallel_points` and `upgrade_thread_independence` alone: passed.
`tests/run_tests.sh`: 165 passed, 0 failed.  T1, T2, T3 and 13N in the three bases (program, menu 1, one thread)
against `tab-gecm-1`: every file identical (tag `tab-pyerr-1`).  Injected errors (Python: the old KeyError, a
RuntimeError for a closed session, no L check, no missing-file check, a message without its call; C++: no NaN
check, an option without a row) each made `errors_test` or `describe_test` fail.  The whole tree is not compared.

**Open.** (1) A new session option still needs a row in four places of the binding and Python and the program's
parser (`src/AZURE2.cpp` has no option table). (2) `AZURE2 --describe` names parameter kinds with words of its own
("energy shift" for `shift`, `experiment/src/Describe.cpp`). (3) A session with `use_rmc=True` keeps `use_brune`
(the program drops it). (4) `kinds()["data"]` has no extra outputs; `azure2.describe()` no parameter bounds.

## 2026-11-13 (agent, branch `feature/upgrade-9-webui-server`) — the web UI's server and process model (no page)

**Done.** `pyazr/ui/` (new subpackage, not imported by `import pyazr`): `server.py` (the app factory, `Server`,
`main`; `python -m pyazr.ui project.azr`), `routes.py` (the route table, 70 rows: 53 served, 9 answering 501 with the
missing pyazr call, 8 with the later step), `state.py` (the text, history, revisions, working file and crash copy,
the root, the session's text), `workers.py` (session and job workers as processes, the budget, death with the
server), `jobs.py` (the job table: scipy and iminuit fits as rows, progress, stop, kill, checkpoints), `errors.py`,
`engine.py` (the adapter, the only file naming `azure2`). Tests `tests/pyazr/ui/` (29, registered in
`tests/pyazr/CMakeLists.txt`); `web_ui.md` section 11.1 (what was built, the measured numbers, eleven departures);
the user guide page `web_ui_api.rst`; `pyproject.toml` extra `web` and the package `pyazr.ui`. No existing pyazr
file, test or engine file changed; no number of the engine.

**How.** The four test files at one thread on the base build's binding (all pass; each rule's test shown failing with
an injected error, `tests/pyazr/ui/README.md`); evaluate equals a direct session bit for bit and the program's
`chiSquared.out` row for row on T1 and T2; a real server process over HTTP, killed during a fit: its workers end in
0.05 s; measured opening T2 and 13N, evaluation overhead, memory, stop times (`web_ui.md` 11.1).

**Open.** (1) The engine resolves `priorCovariance` against the directory of the file it reads, not the session's
base directory: the UI's session files must lie beside the project (departure 1). (2) The χ² per data set has the
six digits of `chiSquared.out` until pyazr gives the program's split (finding 5). (3) The event stream, the file
browser, a job queue and budgets of two or more threads (written, not run). (4) The decisions of `web_ui.md` section
13 stay the PI's; each choice made here is in one named place (section 11.1).

## 2026-11-13 (director, 11:40) — merged: the DW vertex on the transfer library; the `gamma_decay` record in the observed-parameter paths

Two branches, rehearsed together (`rehearsal/trial45`), read by the director (no objection) and verified once on the
rehearsal; the merged code is the rehearsal's (the difference is this file and one line of `web_ui.md`).

**Verified by the director** (build `work/b-trial45`, 0 warnings): all projects in the three bases against the
reference run `tab-gecm-1`: 887 / 851 / 887 files identical, none new; the eight fits identical to the base build
(`work/b-gecm`); `tests/gamma_decay/check.sh` three times at 3 threads (63 ok lines, 0 failures each);
`tests/pyazr/gamma_decay_test.py`; `ctest -LE slow` 224 of 224; the slow tests 11 of 11; `tests/run_tests.sh`
165 passed, 0 failed; spec equivalence 181 of 181 converted equal; the converted-tree gate 679 / 643 / 679 identical;
T1, T2, T3, 13N in the three bases at 1 against 3 threads identical. New reference run: `tab-t45-1`; base build
`work/b-trial45`.

**Accepted choices.** The THM check is 7.4 % slower (inside the 10 % the task allowed); the expansion phases live in
the library; `CNuc::GammaDecayInputs` is a member filled where the level parameters are filled.

**Open, from the two agents (not fixed here).** The engine's sampler has neither a record's priors and wall nor its
limits, in both parameter bases, and keeps its run in process-wide variables (a task is written: the sampler samples
the fit's objective and owns its state). `param.fit` is not reproducible run to run with the numerical gradient at
3 threads. A record and its twin differ by about 1e-10 in a fit. The THM tests print 13 digits.

**Next.** The pyazr errors and self-description branch (rehearsed on top, all projects identical, in verification);
the web UI's server (adapted to that branch, its 33 tests pass on the rehearsal); the prior-centre fix and the
diagnostics record (agents still at work); then pyazr's parameter records and the sampler.

## 2026-11-13 (agent, branch `feature/upgrade-8-gamma-physical-paths`) — the `gamma_decay` record in the observed-parameter paths; E_f in the pull-back

**Done.** Reproduced first (scratch `work/gammapp/`, build of the branch's start): every caller that hands
`GetTransformParams`' vector (the record's numbers at the record's slots) to `FillCompoundFromParamsPhysical`. Wrong:
the paths that set a vector in observed form, pyazr's `calculate_chi2` / `update_segments` and the observed-parameter
MCMC (`CalculateLogLikelihoodPhysical`, its snapshot, `ParkNonPositiveLevels`): `FillMnParams` took the file's
numbers, so the vector's were ignored (χ² 44.79 instead of 3889.3 for Γ_γ × 1.2; the MCMC's likelihood constant).
Read wrongly but in no output: the grid copies in `ReactionRate.cpp` and `TargetSetup.cpp` (0.870 eV of photon
width for Γ_γ = 0.0385 eV; the anchors and the tabulated rate use the particle width; the sum is used only by RMC,
refused for a record, and the rate's quadrature fallback). Harmless: `ConvertPhysicalLimitToReduced` (the converted
prior was the twin's bit for bit), and everything that returns the vector. Rule (note section 2.8, committed
first): the observed vector carries the record's numbers in both directions; `GammaDecaySet::WithPartialWidths` is
the one place that turns them into the levels' widths; `FillCompoundFromParamsPhysical` keeps them as the
compound's input numbers (`CNuc::GammaDecayInputs`), which `FillMnParams` reads. No call site changed. The pull-back
has the term in a free final-state energy E_f (central difference, step of E_λ), `InputIndices` names it; note 2.5's
step corrected to 10⁻⁴ MeV (gradient check at 10⁻⁴ and 10⁻⁶ MeV: both meet the claim), 2.7's sentence on E_f.

**Verified how** (build `work/b-gammapp`, 0 warnings, on c9ee6dd5): the whole tree × 3 bases against `tab-gecm-1`
(tag `tab-gpp-final1`): 887 / 851 / 887 identical, none added; `ctest -LE slow` 223/223; the slow tests 11/11;
`tests/run_tests.sh` 165 passed, 0 failed; the eight fits of `fits_vs_base2.sh` identical to `b-gecm`; T1-T3 and
13N calculations at 1 against 3 threads identical (12 of 12); the converted-tree gate 679 / 643 / 679 identical;
`tests/gamma_decay/check.sh` three times at 3 threads (63 checks, all passed each time, the new (i) and (j)
included); `tests/pyazr/gamma_decay_test.py` all passed (part 5 new). Five injected errors fail the new checks
(README). Not run: Sphinx, TSan/ASan.

**Open**: the MCMC (both bases) has neither the record's priors and wall (term T11 is not in its likelihood) nor its
limits b ∈ [0, 1], Γ_γ ≥ 0 (by reading; samples of the start build reached Γ_γ = −6.6×10¹³); `param.fit`'s
snapshot depends on thread scheduling in a fit with the numerical gradient (two runs at 3 threads differed); the
check (i) reads the grid's line of `AZR_DEBUG_GRID`.

## 2026-11-13 (director, 11:10) — the PI's decision on plots: plotly.js; the bundle fetched (outside the repository)

The PI ("use plotly", after asking whether something else would be better; the alternatives uPlot, Chart.js, ECharts
and D3 were laid out with what each would leave us to write) took decision 11 of `web_ui.md` section 13 as
recommended: the cartesian partial bundle of plotly.js (MIT). Fetched from the npm registry and checked against the
registry's sha512: `plotly.js-cartesian-dist-min` 4.1.2 (published 2026-10-05; 4.3.0 of 2026-11-10 was the newest and
is deliberately not taken), `plotly-cartesian.min.js` 1 507 990 bytes, sha256
`aaf12ff1d4671f7e93f6df4291651ab4ca6b76c0ea5f1576ec1030e3c2819a20`, with its `LICENSE`, under
`work/vendor/plotly/unpacked-4.1.2/package/` (outside the repository). The page step that draws the first plot
copies the file and the licence into `pyazr/ui/static/vendor/` and records version and hash beside them; every plot
goes through one module of the page, so that the library is named in one place. Decisions 5 (drop 20 parity items)
and 7 (session options as `<config>` keys) of that section are still open; the others are built on the note's
recommendations, each in one named place.

## 2026-11-13 (agent, branch `feature/upgrade-4-prior-centre`) — the centre of every nuisance prior is the number written

**Done.** One rule for every nuisance kind with a prior (`parameter_kinds.md` section 4): the penalty is
((p − centre)/width)² with the centre written, the start is the value written, a percentage width is per cent of the
centre; grids and limits stay laid out from the start. Before (calculations, value ≠ centre, base build): norm and
constant shift used the centre; sqrt shift, shift function, abundance, `isotope_sum`, screening, Doppler and thickness
dropped it; areal density, stopping scale and a kernel's parameter refused it. One carrier `PriorCentre`
(`SetupSpecs.h`, the norm/shift centre made general), `build::PriorCentreOf`; the engine objects keep the centre beside
the width; `ForEachDiagonalPrior` hands it to gradient, least-squares rows and MCMC; the API reports it
(`prior_nominal`, now also abundance and apparatus). Write-back: areal density, Doppler, screening and thickness had no
`ParameterSource` (save_fit kept their old values); now they have one, with `engineUnit`. API: a rebuild or mode
switch kept the nuisance priors' manager of the old model (freed objects, old radius): `AZUREAPI::ReplaceModel`.
Tests: `tests/nuisance_prior_centre/check.sh`, `save_fit_format3_test.py` (abundance and density moved),
`tests/pyazr/rebuild_nuisance_test.py`; all fail on the base build.

**Verified.** Build 0 warnings; all projects × 3 bases identical to `tab-gecm-1` (887/851/887, tag `tab-prior-1`);
converted-tree gate 679/643/679 identical; the eight fits identical to `b-gecm`; T1/T2/T3/13N 1 vs 3 threads identical;
`ctest -LE slow` 223/223, slow 9/9 and the two four-thread tests (capped at 3) 2/2; `tests/run_tests.sh` 165 passed,
1 failed (rectangular_resolution expected the lifted refusal; check updated, passes).

**Open.** Grids/limits of a fitted density, temperature, stopping scale or kernel parameter follow the start value: a
fit continued from a written project without `grid_max` lays them out anew (the Doppler pieces too). (Closed by the
branch `feature/upgrade-4-layout-centre`, entry above, for every parameter with a prior.)

## 2026-11-13 (director, 09:20) — merged: the design note of the web UI (documents only)

`docs/upgrade/web_ui.md` is a proposal with eleven decisions for the PI (its section 13), none taken; its twelve
findings on pyazr 3.0 (section 12) go into the briefs of the pyazr steps. Checked by the director: the spike's own
tests rerun at one thread (5 passed); by reading, that `Config::stopFlag` is set by nobody and makes
`AZURECalc::operator()` return 0 in a fit, that the snapshot interval is 100 calls, and that the difference of the
per-data-set chi-squared of T2's s02/s03 between `segment_chi2` and `chiSquared.out` is the convention stated in
`segment_chi2`'s docstring (whitened rows against r_i (M^-1 r)_i; the sums agree), so no wrong number. Not
repeated: the timings, the memory figures, the interrupted program fit, the parity classification of the 201 items.

## 2026-11-13 (agent, branch `feature/upgrade-9-webui-design`) — the web UI: design and a measured spike (no code)

**Done.** `docs/upgrade/web_ui.md` (new): the 201 Qt parity items by tab and state (115 served by today's pyazr, 17
designed in `pyazr3.md`, 20 missing from its design, 29 the page alone, 20 dropped with reasons); the process model
(server with the project text and its history, a session worker, one worker per job, a compute budget); the `/api`
as a route table over pyazr's objects with errors mapped to statuses; the state rule and the cost of each kind of
edit; seven tabs with wireframes; the page (vanilla JS, plotly.js vendored: described, not fetched); security
(loopback, a token); tests; extensibility; six steps of work; twelve findings for pyazr 3.0; eleven decisions for
the PI. One line in `ROADMAP.md` Phase 9. No code of the repository changed.

**How.** A FastAPI spike on today's pyazr in `work/webui/spike/` (outside the repository), run at one thread on
copies of T1, T2, 13N with the build `work/b-gecm`: over HTTP on 127.0.0.1 (server PIDs recorded and checked gone)
and with the test client (5 tests pass). Measured: open, builds cold and warm, the GIL during a build, evaluation
overhead, JSON sizes, a recipe fit's cooperative stop and kill, the program's MIGRAD under SIGINT, a server killed
during a fit, two clients, memory. T2's χ² through the API equals the program's (10905.2 + prior 4.9025 =
10910.1338). `python3 tests/consistency/repo_consistency_test.py` after the commit.

**Open.** (1) The findings of `web_ui.md` section 12 for the pyazr steps (session from the text, progress and stop of
engine fits, writing x into the text, pairs with channel generation, the reaction rate, ...). (2) The decisions of
its section 13, among them vendoring plotly.js (a download) and where session options live. (3) The per-data-set χ²
of `segment_chi2` differs from `chiSquared.out` for data sets joined by a covariance (T2 s02/s03; documented in its
docstring); the UI is to show the program's. (4) A job worker outlives a killed server today (the spike); the design
makes workers die with it.

## 2026-11-13 (director, 09:00) — the PI's decisions on pyazr 3.0; packages installed

The PI accepted the director's list of the eleven decisions of `pyazr3.md` section 5: the note's recommendations,
except 6 (SRIM data compiled into the wheel with the citation notice, on his word that the data are free to use
when SRIM is cited; the owner is to be told). Recorded at the head of that section; the roadmap's and the
architecture's shim sentences struck (decision 9). With the PI's yes, installed into the `azure2-up` environment by
conda (`--freeze-installed`, nothing present changed): fastapi 0.162.5, uvicorn 0.54.0, httpx 0.28.1, iminuit
2.33.0, emcee 3.1.6. Not asked for and not present: plotly.js (to be vendored for the web UI), scikit-build-core.

## 2026-11-13 (director, 08:20) — merged: the design note of pyazr 3.0 (documents only); a defect of prior centres

`docs/upgrade/pyazr3.md` is a proposal: its section 5 lists the decisions that are the PI's, none is taken. Its
inventory was measured at 441a0907 with the scripts of `work/pyazr3/` (outside the repository). The director read
the design sections and did not repeat the inventory's counts.

The defect of its section 1.6, item 1, repeated by the director (`work/director/iso_prior/`, T2, a calculation,
build `work/b-gecm`): `isotope 15N 0.0037 fit prior 0.0037 0.0002` and `... prior 0.0036 0.0002` give the same
total chi-squared 10910.1338428 (0.25 expected between them); start 0.0036 with centre 0.0036 or 0.0037 give
10917.025275 both. The cause, by reading: `AbundanceRowSpec` has no centre, `BuildEffects.cpp` copies the prior's
width only, and `CNuc` sets the nominal value to the start value. The same reading shows three kinds: norms and
energy shifts use the centre written; the stopping scale, the areal density and the apparatus parameters refuse a
centre away from the value ("this engine centres it on the value"); abundances, the Doppler temperature, the
screening potential and a data set's thickness or transmission drop the centre without a word. Consequence: a fit
started again from fitted values loses (or, for the refusing kinds, cannot state) its prior. Not fixed here; a
branch of its own.

## 2026-11-13 (agent, branch `feature/upgrade-8-pyazr3-design`) — pyazr 3.0: measured state and design (no code)

**Done.** `docs/upgrade/pyazr3.md` (new): the inventory of pyazr on 441a0907 measured by scripts (public names, the
binding, `AZUREAPI`, Python formulas, what the engine can do that pyazr cannot reach), four use cases run as a user
and an agent would, the design of 3.0 (objects, arrays, errors, self-description, fit, physics explorer, command
line, packaging), the extensibility of the binding today and after, eight steps of work with sizes and order, eleven
decisions for the PI. One line in `ROADMAP.md` Phase 8. No code changed.

**How.** Scripts and logs in `work/pyazr3/` (listed in the note's section 1.7), with the build `work/b-gecm` at two
threads: 115 of the 120 tests of `tests/pyazr/` under a function-level tracer; the use cases compared with the
program on the same project (identical files where both write them).
`python3 tests/consistency/repo_consistency_test.py` after the commit.

**Open.** (1) Engine defect, not fixed: a target `isotope` whose start value differs from its prior centre loses its
prior (no `Total-Abundance-Chi-Squared`, program and pyazr alike; note section 1.6); through it `save_fit` refuses its
own snapshot after a fit with a free abundance. (2) `calculate_chi2_rwa` includes the abundance prior although its
docstring says data only. (3) `parallel_points_test` and `upgrade_thread_independence_test` run the engine at four
threads whatever the caller's budget. (4) The PI's decisions of the note's section 5 before step 5 of its plan.

## 2026-11-13 (agent, branch `feature/upgrade-6-evaluation-4`) — the environment record and the run's diagnostics

**Done.** Milestone 5 of `evaluation.md` (design section 8, record section 7.6).
- One table of the 39 environment switches (`support/include/azure2/support/Environment.h`,
  `support/src/Environment.cpp`, library `azure2_support`): name, type, value when unset, kind (numbers 14, cost 1,
  prints 22, test 2), what it does, its page; one reader, `ProcessEnvironment()`, read once.  No other `getenv` in the
  engine, its libraries or the binding (check 39); the manual's list (`command_line.rst`, "Environment switches") is
  the table's (check 40).
- The switches that can move a number are copied into the run's `Config` (`Config::environment`, by `Reset`) and read
  there or from the pair's settings; SRIM's directory reaches the stopping library through `TargetStopping`.  The
  printing switches and the test hooks read the record.  The engine no longer sets variables inside the process
  (pool check: a `Config` copy; gradient check: `AZURECalc::GradientBy`).  `thm_adjoint_test.py` runs its
  `AZURE2_THM_FD` half in a child process.
- The counters are the run's: `include/RunDiagnostics.h` (`Config::diagnostics`): profile account (`Config::profile`),
  Coulomb table lookups (printed when the run ends; `ExitRun` before an early exit), the call counters of the checks;
  the EC table's counters in `ECTermStore`.  29 allow-list rows removed; the scanner reports `std::once_flag` members
  (kind `once`), the Coulomb table's content identity is a permanent row.  Pending 53, permanent 15.
- `--describe` ends with the switches that are set.  `--validate`, `--describe`, `pyazr.validate`, `pyazr.describe`
  read data files relative to the project file's directory (`experiment::ProjectDirectoryOf`); the validate test
  converts T2 beside its data.  `pyazr.rst` says what two sessions may do at once and what is still shared.  The
  recipe "Add a diagnostic switch" is in `CLAUDE.md`.
- Merged `feature/upgrade-4-data` at 441a0907 (no conflict).

**Verified** (build `work/b-eval4`, 0 warnings; base `work/b-gecm`, reference `tab-gecm-1`): all projects x 3 bases
identical (tag `tab-ev4-full2`: 887 / 851 / 887); the eight fits identical to the base build; T1, T2, T3, 13N at 3
threads identical to the 1-thread reference; the converted-tree gate (tag `gev4-1`); every switch set, base against
branch, the same outputs and the same log but for times and progress bars (`work/eval4/switch_results.tsv`);
user-space instructions not higher on T1, T2, 13N, natural_neon calculations and T2, 13N gradient fits; `ctest -LE
slow` 224 of 224; the slow tests 11 of 11 (`pyazr_concurrent_sessions` included; `pyazr_session_isolation` is in the
first set); `tests/run_tests.sh` 165 passed after two tests that compare the whole `--describe` text, or ran it on a
project away from its data, were fixed (`ratio_to_rutherford`, `units_scales`).

**Open.** `api/src/AZUREAPI.cpp` (another agent's file tonight): the unused `intEC_cache` file name and the lazily
read nuisance layout (`mutable`) stay.  The PI's decisions of `evaluation.md` section 6 (options for the *numbers*
switches, removal of the unused ones) are not taken.

## 2026-11-13 (director, 07:32) — merged: the `gamma_decay` record on a level whose transitions have external capture

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (nothing pushed):
`feature/upgrade-8-gamma-external-capture` (the agent's entry below says what was done, by file). No output file of
a project of the tree changes; reference tag `tab-gecm-1`, base build `work/b-gecm`.

**What it is.** A record's partial width on a channel with external capture is the observed width, the external
part included, as a width in `<levels>` is. At every fill the map evaluates γ^ext at the current parameters and
subtracts it as `TransformIn` does, so the record's widths stay what the record says while a fit moves a proton
width, the level energy or an ANC of the final state (`gamma_parametrisation.md` section 2.7). The formal basis
stays refused.

**Verified by the director** (rehearsal merge on the third Phase 6 branch, build `work/b-gecm`, reference
`tab-eval3m-1`).
Verified by the director on the rehearsal merge (this branch on the head after the third Phase 6 branch; build with 0
warnings): all projects in three bases against the tree of the head before (887 / 851 / 887 files) identical, none
new; the eight fits identical to the base build's; tests/gamma_decay/check.sh three times at 3 threads, all checks
passed each time; tests/pyazr/gamma_decay_test.py passes; calculations at 1 against 3 threads identical; ctest without
the slow label 223 of 223; the slow tests 11 of 11; tests/run_tests.sh 165 passed, 0 failed; spec equivalence 181 of
181 converted files equal; the converted-tree gate 679 / 643 / 679 identical. The director's own gradient check
against central differences on a 15O project with external capture on two transitions (not the agent's parameter
points; Brune and Park, the Park run reaching the final state's spin term): particle width 1.6e-8, final-state ANCs at
most 7.9e-9, the record's numbers at most 3.3e-11, the level energy 8.8e-9 at a step of 1e-8 (2.4e-4 at the default
step, the known noise of the mesh placement). Not run: TSan or ASan (not available), the Sphinx pages, MINOS and the
MCMC on such a record.

**The director's own gradient check** (`work/director/vgec/own/`), at a point not in the test (Γ_p 800 eV,
E_λ 7.5565 MeV, ground-state ANC 0.40, Γ_γ 0.02 eV, b(6.1763) 0.2, b(g.s.) 0.1), every free parameter, against
central differences. Brune: the proton width 1.6e-8, the three ANCs of the final states at most 7.0e-9, the
record's numbers at most 2.4e-11, the level energy 8.8e-9 (step 1e-8). Park, whose term of N_f in the final
state's amplitudes (J_f) the test's gradient points do not reach: the ANCs at most 7.9e-9, the proton width
7.5e-10, the record's numbers at most 3.3e-11, the level energy 3.7e-9.

**Reviewed** by reading: `src/GammaDecayMap.cpp` (the map, the floor, the pull-back's factors R and B against
Eq. 2.8), `ExternalReducedWidthAt` and `ExternalCaptureFinalState` (the old search for the final state and its
mask, unchanged in what it finds), the gradient's mask. Notes: (1) the pull-back has no term in the final state's
energy; the note's reason ("a bound state's energy is never free") holds only for a level below every particle
threshold of its nucleus (`CNuc::FillMnParams`), not for a final state above a second pair's threshold or an
energy freed through pyazr: next branch; (2) the set keeps a pointer to the run's `Config` (for the integrals'
scratch), which must outlive the compound and its clones: true for the program and for a pyazr session by
reading; (3) a record fit costs about 14 % more than its twin (agent's figure).

**Defect found, not fixed here** (agent's report, confirmed by the director by reading, not by a run):
`CNuc::GetTransformParams` returns the record's numbers at the record's slots, and `ReactionRate.cpp`,
`TargetSetup.cpp`, `ParameterLimitsManager.cpp`, the MCMC's observed-parameter path and the API hand that vector
to `FillCompoundFromParamsPhysical` as observed widths. At a narrow resonance a branching 0.5 is then read as a
width of 0.5 eV where the integration grid is placed. Present since the record's merge (04:03 today). Next
branch, with a reproduction first: `work/briefs/BX_gamma_record_physical_paths.md`.

**Not run by anyone**: Sphinx; MINOS and the MCMC on a record with external capture.

## 2026-11-13 (agent, branch `feature/upgrade-8-gamma-external-capture`) — the `gamma_decay` record on a channel with external capture

**Done.** The refusal of a record whose photon channel has external capture in the run is lifted for the Brune and
Park bases (the formal basis stays refused). Design: `gamma_parametrisation.md` section 2.7 (committed first). The
record's partial width of such a channel is the observed total, as in `<levels>`; at every fill the map forms the
total amplitude as before, evaluates γ^ext at the current parameters (the level's energy and particle amplitudes, the
final state's ANCs and N_f; Park: n0 = J_f) and subtracts it with `physics::subtract_external_capture`, as
`TransformIn` does. A width below 2P(Im γ^ext)²/n follows `TransformIn` (γ_int = −Re γ^ext, a warning once per run, a
note in `parameters.out`; no wall). Code: `ExternalCaptureFinalState` and `ExternalReducedWidthAt` (values in, exact
slopes out) in `src/CNucExternalCapture.cpp`, of which `ExternalReducedWidth` is now a call (bit-identical);
`src/GammaDecayMap.cpp` (per-slot final states, `Map`, the fill's `GammaDecayFill` in the compound, the pull-back's new
terms in the level's particle amplitudes and the final states' amplitudes, the factor dγ_int/dg, the energy's central
difference through the whole map and only when the energy is free, `InputIndices` for the gradient's mask in
`src/AZUREGrad.cpp`). The radial integrals are the direct ones with the run's per-thread scratch (as `TransformIn` and
`TransformOut`), not the energy tables.

**Verified how** (build `work/b-gammaec`, 0 warnings). `tests/gamma_decay/check.sh` (h): the ¹⁴N(p,γ)¹⁵O twin with
capture into the ground state and 6.1763 MeV (E1, g_ext −0.0049 and +0.114 MeV^½) and 6.7931 MeV (M1, g_ext 0):
every output file identical to the twin in both bases, also after moving Γ_p, E_λ and the ground state's ANC; the
gradient in all free parameters (three ANCs of final states) at three points, largest relative difference 3.6e-8
(level energy 2.4e-7 with step 1e-8); gradient MIGRAD and LM from displaced starts: χ² 237.111 as the twin, the LM
branchings equal the twin's Γ_f/Γ_γ to 1e-6; 1 = 3 threads; the floor case. Six injected errors each fail (README).
pyazr: `tests/pyazr/gamma_decay_test.py` part 4. Whole tree × 3 bases against `tab-gammam-1` (tag `tab-gec-1`):
887 / 851 / 887 identical, no project added; the eight fits identical to `b-gammam`; the converted-tree gate
(`ggec-1`): 679 / 643 / 679 identical; `ctest -LE slow` 221/221; the slow tests 11/11; `tests/run_tests.sh` 162
passed, 0 failed; check.sh three times
at 3 threads: all checks passed each time; calculations at 1 against 3 threads: 74 files identical (T1, T2, T3, 13N,
13N_capture_ay, thick_target_resonance and the external-capture twin and record). Cost (instructions:u, 1 thread):
calculation +0.25 %, gradient MIGRAD fit +14 % (89 calls against 85).

**Open.** Not run: Sphinx (no sphinx in the environment), MINOS and MCMC on a record with external capture. Found
outside the scope (not fixed): `CNuc::GetTransformParams` returns the record's numbers at its slots, and
`ReactionRate.cpp` (rate grid) and `TargetSetup.cpp` (grid copy) feed it to `FillCompoundFromParamsPhysical` as
observed widths, so a record level's photon "widths" there are Γ_γ, b_f, δ_f (read for grid placement; not checked
by a run); note section 2.5 says the energy step is 10⁻⁶ MeV, the code uses 10⁻⁴.

## 2026-11-13 (agent, branch `feature/upgrade-7-thm-dw-vertex`) — the DW vertex's source projection is the library's

**Done.** (1) Design in `transfer.md` §25 (what the projection returns, how `EntranceVertex`'s surface term is built
from it, what moves and what stays, the named settings, tests, extension points), committed before code.
(2) `physics/transfer/SourceProjection.h/.cpp`: `PrepareSourceProjection` (the plan: u and cos θ_u quadratures, φ at
the nodes, the triples (l; L_s, L_a) with their coefficients, the entrance waves), `SolveSourceExit`, the sum over
the points in blocks whose order the plan fixes (`SourceSumBlocks`, `SourceBlockSize`, `SumSourceBlock`,
`CombineSourceBlocks`; `ProjectSource` does it on one thread), `SourceMultipolesAt` (S_lm(a), S′_lm(a) at a
direction), `SourceMomentumAmplitude` (φ̃(q) on the same nodes). Pieces added to the library and used by it:
`SixPointLagrange` (value and derivative; `InterpolateOnGrid` now uses it), `SphericalHarmonicsWithDerivative`
(any sign of θ; the recurrence shared with `SphericalHarmonicsAtPhiZero`), `clebsch_gordan_integer`
(`AngularCoupling.h`, integer j without the j ≤ 40 table). (3) `src/ThmDwVertex.cpp` calls it: the engine keeps the
energy grid and its passes (OpenMP over the exit waves and over the blocks), the directions, the normalisation by
4πφ̃(q), the Gram matrix, its factorisation and interpolation; its quadrature, `Lagrange6`, `YTable`, `WaveTable` and
`ThmCG` are gone (src/ and include/: −355, +68 lines). `EntranceVertex` is unchanged (§25.2: it takes any wave; making
it call the projection would replace the exact plane wave of `PlaneWaveEntranceVertex` by a quadrature); a test builds
it from the projection instead. (4) The two zero-range implementations are tied by a direct test (3e-15).
(5) `CoulombSphere` is finite at r = 0 (one function `CoulombSphereTerm` in `HybridPotential.h` for the sphere,
`RadialPotential` unchanged bit for bit; `TabulateChannel`'s value at the origin from `TermsAt` at r = 0 instead of its
own formula beside a call at 1e-300). Not done: the optional spin-orbit terms of `kd03`.

**Verified** on the last code commit (build `work/b-dwvertex`, 0 warnings), against `tab-gecm-1` and `work/b-gecm`:
the whole tree in three bases 887 / 851 / 887 files identical (tag `tab-dwv1`); the converted-tree gate 679 / 643 /
679 identical; `ctest -LE slow` 224 of 224; the slow tests 11 of 11; `tests/run_tests.sh` 165 passed, 0 failed; the
eight fits identical to the base build file by file; T1, T2, 13N, 18O_p_a_thm at 3 threads identical to the tree run at
1; the checks thm_dw_vertex and thm_fit pass and their outputs differ from the base build's by at most 2.0e-10 relative
(the fit's EDM, a diagnostic of size 1e-8, by 8e-4); `thm_dw_vertex_test -v` 644 quantities, largest change 2.7e-11;
`thm_distortion_test`, `thm_optical_test` unchanged at their printed digits. Each new test failed once with an
injected error (named in the report). Speed: `transfer.md` §25.6.

**Open.** The +7 % instructions of the vertex check (per-point overhead of the library's pieces, largest for l = 0);
`kd03` spin-orbit; a spin-orbit distortion in the vertex (waves per j, spins in the triples, §25.5).

## 2026-11-13 (director, 07:10) — merged: Phase 6, third branch (a fit's copies made when needed; the Coulomb table named by content; relative paths are the project's)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (nothing pushed): `feature/upgrade-6-evaluation-3` (the
agent's two entries "Phase 6 branch 3" below say what was done, by file). New reference tag `tab-eval3m-1`, base
build `work/b-eval3m`.

**What it is.** (1) The fitter and the sampler no longer build a pool of working copies up front: a copy is made
when no idle one exists, so a fit holds one copy of the problem less. (2) The signature of `intEC.dat` names the
Coulomb table by its content (`CoulombTable::ContentIdentity`: format, order, l_max, the axes, the generator and a
64-bit hash of every sheet's values) or says "exact Coulomb functions"; an `intEC.dat` written before is recomputed
once. (3) A project's relative paths are opened against the project file's directory (the command-line program) or
the directory the pyazr session names, by one function (`azure2::experiment::PathToOpen`, `Config::OpenPath`); pyazr
no longer changes the process's working directory. The project file itself, `--coulomb-table` and
`AZURE_OPTIONS_FILE` on the command line stay relative to the working directory.

**What a user sees.** `AZURE2 /some/where/project.azr` started in another directory now reads and writes where the
project is (before: relative to where it was started, so it failed or wrote elsewhere); `nonfinite_dump.txt` is
written in the project's directory; the 23 `intEC.dat.sig` files per basis of the tree have a new last line.
`--validate` and `--describe` (and `pyazr.validate`, `pyazr.describe`) still resolve against the working directory:
decided by the director that they follow the project's directory too (next branch).

**Verified by the director** (rehearsal merge on the THM merge, build `work/b-eval3m`, reference `tab-gammam-1`).
Verified by the director on the rehearsal merge (this branch on the head after the THM merge; build with 0 warnings):
all projects in three bases against the tree of the head before (887 / 851 / 887 files): 864 / 828 / 864 identical;
the 23 that differ per basis are all `intEC.dat.sig` (the signature now names the Coulomb table by its content), no
other file differs and none is new; the eight fits identical to the base build's except the `intEC.dat.sig` of the two
13N fits (same reason, the hash line only); the checks ec_table_signature, run_from_elsewhere, layer_yield,
extrap_sav_segments and describe_relativistic pass at 3 threads; calculations at 1 against 3 threads identical; ctest
without the slow label 223 of 223; the slow tests 11 of 11; tests/run_tests.sh 165 passed, 0 failed; spec equivalence
181 of 181 converted files equal; the converted-tree gate 679 / 643 / 679 files with the same 23 signature files per
basis as the only differences. Not run: TSan or ASan (not available), the Sphinx pages, the agent's memory figures
(not measured again).

**The director's own check of the path rule.** Every project folder of the tree (73, Brune basis) run from an empty
directory elsewhere with the project's absolute path: all 832 output and checks files byte-identical to the run from
the project's directory, no file left in the foreign directory; 55 `run.log` differ only in the project's name as
typed and in progress lines.

**Reviewed** by reading: `ProjectPaths.h`, `Config`, `CoulombTable`, `ExternalCaptureSetup.cpp`, `pyazr/azure2.py`,
`AZURECalc` and `AZURECalcMCMC`, the allow-list, the API's diff, and every remaining raw file open of `src/`,
`api/src` and `experiment/src` (each goes through `OpenPath` / `PathToOpen` or was project-relative before). Notes:
(1) the four rows of `InputFiles.cpp` (the bytes of input files, keyed by absolute path, size and modification time)
were moved by the agent from pending to permanent "table:"; accepted by the director as a content-keyed cache, for
the owner to see (`cache_stats()` counts it per process); (2) `ContentIdentity` is computed lazily under
`std::call_once` behind a pointer, which the state scanner does not see: to be listed or made eager in the next
branch; (3) the memory figures in the agent's entry (a T2 gradient fit 192 to 175 MB at one thread, natural_neon
802 to 711 MB) are the agent's, measured before the last merges, not repeated by the director.

**Not done** (brief `work/briefs/BW_phase6_diagnostics_record.md`, next branch): the diagnostics record (0 of 44
switches; the 29 milestone-5 rows of the allow-list), the `mutable` members of `AZUREAPI`, a sentence of
`docs/source/user_guide/pyazr.rst` that still tells users to drive sessions from one thread, an unused `intEC_cache`
name in `api/src/AZUREAPI.cpp`. The allow-list now has 82 pending rows (milestone 3: 11, 4: 42, 5: 29) and 14
permanent ones.

## 2026-11-13 (director, 06:57) — merged: the THM distortions and the DW vertex's waves from the transfer library; global optical sets one table

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (nothing pushed): `feature/upgrade-7-thm-wiring` (the
agent's entry below says what was done, by file). No output file of a project of the tree changes (no project of the
tree uses `distortion=` or `vertexModel=dw`); the quantities of the THM unit tests move by at most 1e-8 relative.

**New reference tag of the whole tree: `tab-thmwm-1`** (the director's run of the merged tree, build
`work/b-thmwm`; the files are those of `tab-gammam-1`).

**What it is.** `src/ThmDistortion.cpp` keeps the THM kinematics, the window nodes, the reference energy and the
tables per stage, and has no solver: waves, matching, the Trojan-horse tail, the zero-range weight's partial-wave sum
and its plane-wave limit are `physics/transfer` (`DistortedWaves.h`: `SolveStoredWaves`, `MatchingFor`;
`ZeroRangeWeight.h`: `PrepareZeroRangeEntrance`, `ZeroRangeWeightSum`, `ZeroRangeMomentumAmplitude`; `Overlap`'s
`r_min_fm`). The seven global optical sets are one table (`GlobalOptical.h/.cpp`, a row and a test per set; recipe
in `CLAUDE.md`); `src/ThmOptical.cpp` and its header are gone. `ThmDwVertex` takes its waves from the library and
keeps its own 2D quadrature and Gram matrix.

**Verified by the director** (rehearsal `wt-trial42`, build `work/b-thmwm`, against `tab-gammam-1` and `b-gammam`).
Whole tree x 3 bases (tag tab-thmwm-1) against tab-gammam-1: 887 / 851 / 887 files identical, none new (no project of
the tree has a distortion). The eight fits identical to the base build, file by file. The checks thm_distortion,
thm_dw_vertex, thm_fit, thm_spectator_angles and thm_coulomb_consistency pass at 3 threads. ctest -LE slow 223/223;
the slow tests 11/11; tests/run_tests.sh 162 passed, 0 failed; spec equivalence 181 of 181 converted projects equal;
the converted-tree gate 679 / 643 / 679 identical; calculations at 1 against 3 threads identical. The director's own
comparison of the unit tests' quantities, old build against new: thm_distortion_test 82 quantities, 47 moved, largest
2.9e-9 relative; thm_optical_test 14 common quantities, largest 9.1e-9; thm_dw_vertex_test 644 quantities, 83 moved,
largest 6.7e-10; all below the 1e-7 of the accepted rule.

**The moved numbers, measured by the director** (`work/director/vthm/unit/`: the three unit tests run with `-v` in
`work/b-gammam` and `work/b-thmwm` at 1 thread, every printed quantity compared): `thm_distortion_test` 82
quantities, 47 moved, largest 2.9e-9; `thm_optical_test` 14 common quantities, all moved, largest 9.1e-9 (the other
230 of the old test are the parameter values, now 23 points of ten numbers in
`physics/tests/transfer_global_optical_test.cpp`, passing); `thm_dw_vertex_test` 644 quantities, 83 moved, largest
6.7e-10. All are below 1e-7, the limit of the accepted rule for a moved number (the PI's rule for the fold
quadrature, applied here by the director's decision; the owner has to see it). No golden value was re-baselined.

**Reviewed by reading**: `src/ThmDistortion.cpp`, the diff of `src/ThmDwVertex.cpp`, `ZeroRangeWeight.h/.cpp`, the
changes of `BoundState.cpp`, `Overlap.cpp`, `RadialNumerics.cpp`, `Transfer.cpp`. Notes: (1) the library now has two
implementations of the zero-range orbital integral for an l = 0 overlap (`ZeroRangeOrbitalIntegral`, with
Clebsch-Gordan coefficients and j ≤ 40, and `ZeroRangeWeightSum`, by the Legendre recurrence without that limit),
tied to each other only through the FRESCO case at 1 %: a direct equality test, or one calling the other, is owed;
(2) `CoulombSphere` (`physics/src/transfer/OpticalPotential.cpp`) returns NaN at r → 0 for a positive Coulomb
radius; the wiring steps around it, the library has to be fixed; (3) the |M|² of an optical case is converged to
about 5e-7 only, in the old code and in the new (radial step and matching radius each about 3e-6, of opposite
sign); (4) the global sets are central terms only.

**Not done: the DW vertex's quadrature in the library.** The engine's vertex is the surface bracket
(B − 1) S_lm(a) − a S′_lm(a) of the source's projection; the library's `EntranceVertex` returns that term together
with volume and interior terms that hold for a solution of the free equation and needs a level function, and no
public function returns S_lm(a). The roadmap's box stays open; the continuation is a public source projection in
`physics/transfer` that both `EntranceVertex` and the engine call (`work/briefs/BV_thm_dw_vertex_projection.md`).

**Not run by anyone**: a refined calculation of the old engine (the moved quantities are under 1e-7, where the rule
does not ask for one); a THM project with a distortion in the tree's comparison (there is none; the checks of
`tests/thm_distortion`, `tests/thm_dw_vertex` and `tests/thm_fit` cover it).

## 2026-11-13 (agent, branch `feature/upgrade-7-thm-wiring`) — the THM distortions call physics/transfer; global optical sets one table

**Done.** (1) Design in `transfer.md` §24 (what each THM quantity comes from, what the engine keeps, the DW vertex,
the named settings). (2) The seven global optical sets are one table of the library (`GlobalOptical.h/.cpp`), test
`physics_transfer_global_optical` (RIPL-3/FRONT21 values); `src/ThmOptical.{h,cpp}` deleted; recipe "Add a global
optical set" in `CLAUDE.md`. Tree identical after this part (860/833/860 against `tab-truem-1`). (3) The library gains
waves one l at a time (`TabulateChannel`, `SolvePartialWave`, matching at the last node or beyond the turning point),
exact plane waves, r_min on overlaps, the zero-range THM weight with its bound (`ZeroRangeWeight.h`), the a_3 term of
the Numerov start; test `physics_transfer_zero_range_weight` (closed forms, the two matchings, FRESCO's zero-range
¹²C(⁶Li,d) through the new sum to 6.3e-5). `ThmDistortion` keeps kinematics, windows, E_ref and the ln R grid and has
no Numerov, matching or partial-wave sum; `ThmDwVertex` takes its waves from the library. Unit-test quantities move
≤ 9e-9 (record and refinement table in `transfer.md` §24.5); no golden changed. Speed: the THM checks −80 % / −71 % /
−7 % instructions (distortion / fit / dw vertex). Merged `feature/upgrade-4-data` at e4fb58f1.

**Not done (stopped by the brief's rule).** The DW vertex itself: `EntranceVertex` is not the engine's surface-only
vertex for a distorted source, and the library has no function for S_lm(a) (§24.2). Spin-orbit terms of the global
sets (no primary source at hand; Frescox's `opticals/KDParamCode.f` has Koning–Delaroche's, not built here).

**Verified how** (last commit, merged onto e4fb58f1, build 0 warnings): whole tree against `tab-gammam-1` 887/851/887
identical; `ctest -LE slow` 223/223, slow 11/11; `tests/run_tests.sh` 162 passed; the eight fits identical to
`work/b-gammam`; converted-tree gate 679/643/679 identical; 1 against 3 threads identical (T1, T2, T3, 13N,
18O_p_a_thm and two distortion projects × 3 bases); the THM fits of `tests/thm_fit` same χ², parameters, calls.
**Open**: the source projection S_lm(a) in the library; matching of optical waves just beyond R + 12a costs ~3e-6
in |M|² (old and new alike, §24.5).

## 2026-11-13 (director, 04:03) — merged: a level's radiative decay as total width, branchings and mixing ratios (record `gamma_decay`)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (nothing pushed):
`feature/upgrade-8-gamma-branchings` (the agent's four entries below say what was done, by file). Asked for by the PI
on 2026-11-12 ("I would also add a possibility of fitting branchings directly"). No output number of an existing
project changes; three projects are added (`tests/gamma_decay/models`).

**New reference tag of the whole tree: `tab-gammam-1`** (the director's run of the merged tree, build
`work/b-gammam`); against `tab-truem-1` no file of an existing project differs.

**Verified by the director** (own build of the branch merged onto the head, 0 warnings).
Whole tree at three bases against `tab-truem-1`: 860 / 833 / 860 files identical, no file of an existing project
differs; the 27 / 18 / 27 added files are the three projects of `tests/gamma_decay/models`. The eight fits (T2 with
Levenberg-Marquardt, gradient, MIGRAD; 13N gradient and Levenberg-Marquardt; the layer and the two stopping fits)
identical file by file to the base build `work/b-truem`, `fitStatus.out` included. `tests/gamma_decay/check.sh` three
times at 3 threads: all checks passed each time. `ctest -LE slow` 221/221; the slow tests 11/11; `tests/run_tests.sh`
162 passed, 0 failed; the specification equivalence 181 equal; the converted-tree gate against `tab-truem-1` 656 / 629
/ 656 identical, none differs; calculations at 1 against 3 threads identical (12 of 12, and the model-data projects).

**Reviewed by reading**: the map `src/GammaDecayMap.cpp` (partial widths, photon amplitudes, the pull-back of the
gradient, the penalty rows), its call sites in the fill of the compound nucleus and in the gradient's scatter, the
call counter. Two findings went back to the agent and are fixed in the branch: (1) the gradient in a mixing ratio was
exactly zero at δ = 0 (the pull-back divided by the L + 1 width), so a Levenberg-Marquardt fit started there did not
move: now the derivative is taken in amplitude form, test (g) of `tests/gamma_decay/check.sh`; (2) the MIGRAD call
count in `fitStatus.out` differed by one between 1 and 3 threads, because Minuit2's counter is a plain integer
incremented from its own OpenMP gradient loop: Minuit2 is not patched (the PI's decision), AZURE2 counts the calls
of its objective itself, test (f).

**Decided by the director** (the note's five items, each as the agent recommended): the record's and keys' names
stay; the formal basis (`--no-brune`) stays refused with a message until someone needs it; external capture on a
record's channel stays refused in this branch and is lifted with the next work on the capture term; the remainder's
wall (width 1e-3 on a sum of branchings above 1) stays; Rose-Brink's sign convention for δ is not offered unless the
owner asks (the record is Krane-Steffen's, as ENSDF).

**Kept by the PI's decision (2026-11-12), not to be raised again**: the radiative widths stay out of a level's total
width; a `resonance_strength` datum uses the width including the photons, a difference of Γ_γ/ΣΓ_particle, which the
user guide states.

**Not run by anyone**: the Sphinx build of the new user-guide page (no sphinx in the environment); MINOS and the MCMC
with a branching at its boundary; `--max-calls` still acts on Minuit2's counter (approximate by a few calls under
threads). The volume and page of the Krane-Steffen-Wheeler reference in the note were not checked against the
paper.

## 2026-11-13 (agent, branch `feature/upgrade-8-gamma-branchings`) — the gradient in a mixing ratio at δ = 0

**Done.** `GammaDecaySet::PullBack` took ∂p/∂δ as p/(2Γ)·∂Γ/∂δ, which is 0/0 for the L+1 channel at δ = 0 and was
skipped: the gradient in δ was exactly 0 there (AZURE_GRAD_CHECK on an angular distribution: analytic 0, finite
difference −1.29e6), so an LM fit started at δ = 0 did not move. Now ∂p_L/∂δ = −p_L δ/(1+δ²) and ∂p_L+1/∂δ =
K/(1+δ²)^{3/2}, K = 2^½ p_L+1 from the map at δ = 1. Test `tests/gamma_decay/check.sh` (g) on
`models/sign_ad.azr` (the sign case's dσ/dΩ as data): the gradient at δ = 0 and 0.3 (2e-11, 4e-11) and MIGRAD and LM
fits from 0 ending at +0.5; on the old code the gradient check at 0 and the LM fit failed (MIGRAD with the gradient
still reached 0.5 through Minuit's own line search). Note §2.5/2.6 and the user guide say δ = 0 is regular and that
angle-integrated data alone leave the sign of δ open.

## 2026-11-13 (agent, branch `feature/upgrade-8-gamma-branchings`) — MIGRAD's call count is AZURE2's own

**Done.** `fitStatus.out` and the MIGRAD log lines report the calls counted by the objective
(`AZURECalc::CountCallsInto`,
a `std::atomic<long>` owned by the fit run in `AZUREMain`, reset before each MIGRAD), not Minuit2's `NFcn()`, whose
`MnFcn::fNumCall` is a plain int incremented from the threads of Minuit2's OpenMP numerical gradient (lost increments:
65 or 66 for the gamma_decay record fit at 3 threads, 62 or 63 for its partial-width twin, also with the base build).
`external/` is not touched; `--max-calls` still acts on Minuit2's counter (approximate by a few calls under threads).
`tests/gamma_decay/check.sh` (f): the record's fit at 1 and 3 threads, every output file identical.

**Verified how.** At 1 thread the count equals `NFcn()` (record 66, twin 63 = `b-head3`); ten 3-thread runs each: 66
and 63 every time; check (f) failed on its first run with `NFcn()` put back; the eight fits of `fits_vs_base2.sh`, the
tree x 3 bases against `tab-gridm-1` (827/800/827 identical), `ctest -LE slow` 221/221, consistency.

## 2026-11-12 (agent, branch `feature/upgrade-8-gamma-branchings`) — implementation: the `gamma_decay` record

**Done.** A level's radiative widths can be written as Γ_γ, branchings b_f (n − 1 of them and a named remainder) and
mixing ratios δ_f (Krane-Steffen) in a record `gamma_decay` of `<gammaDecays>`: `experiment/` (`GammaDecay.h/.cpp`,
the kind in `Schema.cpp`, five codes), the engine's map in `include/GammaDecayMap.h`, `src/GammaDecayMap.cpp` (set-up
widths, parameters of kind 13 `radiative` at the photon slots, `ToEngine` in `FillCompoundFromParams` through
`amplitudes_from_observed`, the gradient's pull-back in `GradAccum::Scatter`, objective term T11 for the priors and the
remainder's wall, the block of `parameters.out`), pyazr's `Parameter.decay_*`; docs (`format3.rst` generated,
`user_guide/gamma_decay.rst`, `output_files.rst`, `parameter_kinds.md`, `objective.md`, `FOR_THE_OWNER.md`,
`ROADMAP.md`, the recipe in `CLAUDE.md`); tests `tests/gamma_decay/check.sh`, `tests/pyazr/gamma_decay_test.py`, cases
in `experiment/tests/diagnostics_test.cpp`, `parameter_kinds_test.py`.

**Verified how.** Per commit: build without warnings, repo_consistency, the experiment and pyazr tests concerned, T1,
T2, T3, 13N, 13N_capture_ay, thick_target_resonance, model_data x 3 bases identical to `tab-gridm-1`. The new checks
pass and each was shown to fail with one injected error (listed in `tests/gamma_decay/README.md`). End of branch
(2026-11-13, b-gamma, 3 threads): build 0 warnings; all projects x 3 bases against `tab-gridm-1` (tag
`tab-gammab-final`): 827 / 800 / 827 of the reference identical, the other files are the three new projects
(only-new); `ctest -LE slow` 221/221; `tests/run_tests.sh` 161 passed, 0 failed; the converted-tree gate (`ggammab1`)
628 / 601 / 628 identical; the eight fits of `fits_vs_base2.sh` identical to `b-head3`; the record's fit at 1 and 3
threads: every output identical but `fitStatus.out`'s MIGRAD call count (66 against 65; with partial widths 66 and 66),
cause not found.

**Open.** The formal basis and external capture on a record channel are refused (decisions 2 and 3 of the note);
Sphinx was not built here (no sphinx in the environment).

## 2026-11-12 (agent, branch `feature/upgrade-8-gamma-branchings`) — design note: total radiative width, branchings, mixing ratios

**Done.** `docs/upgrade/gamma_parametrisation.md`: what exists (photon widths in `<levels>`, their amplitudes in the
three bases, external capture, priors, what the datum kinds use for the total width), the parametrisation Γ_γ, b_f,
δ_f and its map to the engine's amplitudes (`amplitudes_from_observed`, the engine's own conversion, at every fill),
the Jacobian, the constraint (n − 1 branchings and a named remainder, so the vector keeps one parameter per photon
channel), the record `gamma_decay` of a new section `<gammaDecays>`, and the decisions for the director.

**Verified how.** The sign of δ by the engine's angular distributions (base engine `b-head3`): seven cases (M1/E2 and
E1/M2, J = 3/2 and 5/2, J up and down, both parities), the engine against Krane-Steffen-Wheeler to 2e-7 in a2/a0 and
a4/a0; widths of equal sign are δ_KSW < 0 in all seven (note section 2.3; scripts in `work/gammab/sign/`).

**Open.** The implementation (Part B of the brief); item 4 (b) of the brief dropped by the PI (2026-11-12): the
radiative widths stay out of the total width.

## 2026-11-13 (director, 03:20) — merged: an energy shift is the true beam energy (stopping, thickness and window of a shifted point at E + s)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (nothing pushed):
`feature/upgrade-5-true-energy-stopping` (the agent's entry below says what was done, by file). Decided by the PI on
2026-11-12 ("I would fix it"). **Numbers change**: every folded data set that has both an energy shift and a target
integration with an energy-dependent stopping, by (d ln ε/d ln E) × shift/E on a yield on a resonance.

**New reference tag of the whole tree: `tab-truem-1`** (the director's run of the merged tree; `tab-gridm-1` is the
tree before). Against `tab-gridm-1` only the files of T2 and of `tests/shift_sqrt_subpoints` differ, and the files of
the new models are added.

**Verified by the director** (own build `work/b-truem` of the branch merged onto the head, 0 warnings).
Whole tree × 3 bases (tag `tab-truem-1`) against `tab-gridm-1`: 805 / 778 / 805 files identical; the 22 that differ
per basis are T2's (T2.azr 1.2e-6, T2_true 2.07e-3) and those of five `shift_sqrt_subpoints` models (0.96 %), and
33 files of the new models are added. T2.azr (zero shifts, fitted: the planned support is wider) against the eightfold
refined run of the old engine: largest distance 1.000 and rms 0.999 times the old engine's (the acceptance rule).
The eight fits against the build of the head before: T2's three routes change to χ² 297.860, the five others
(13N by two routes, layer yield, stopping yield by two) are identical file by file. `ctest -LE slow` 220 of 220, the
slow tests 11 of 11, `run_tests.sh` 161 of 161 (`shift_true_energy` among them), spec equivalence 181 equal, the
converted-tree gate 613 / 586 / 613 identical with 15 files per basis that differ, all of T2 and
`shift_sqrt_subpoints`; calculations at 1 against 3 threads identical (T1, T2, T3, 13N × 3 bases).

**An independent check by the director** (not one of the branch's tests): `tests/upgrade/T2/T2_true.azr` against a
twin written by hand, the shift lines of data sets 1 (thin TiN target, offset 2 keV) and 4 (thick target,
0.0015 √E) removed and their data energies moved to E + s. With the engine before the branch the shifted project and
the twin differ by 2.07e-3 (data set 1, at 0.189 MeV c.m.) and 1.45e-3 (data set 4, at 0.235 MeV); with the branch
they agree in every printed digit, and the twin's numbers are the same on both engines. The two differences are the
changes the agent reports for T2_true (2.08e-3 and 1.45e-3).

**Reviewed.** `ShiftedFold.h` (the rule, questions 1, 3, 5), `ShiftedFold.cpp` (`MoveTargetToShift`,
`SupportThicknessChange`, `PlannedShiftSupport`, `MovedSupport`), `TargetEffect::SupportThickness` and
`KernelThickness`, the call in `Segment::UpdatePointEnergyWithShift`, `PointFold.cpp`, the set-up's and the re-plan's
call of the planner, the test's README and `check.sh`. Notes: (1) the thickness is recomputed at the true energy with
the working density and stopping scale, and the lab energy it was computed at is stored, so a fitted density or scale
recomputes it at the same energy; (2) the window's growth is 0 bit for bit without a shift, for a constant stopping
and for a kernel without a target, which is why no other project moves; (3) the bound on the planned support's low
edge rests on three samples of the thickness change over the reach plus their spread: an argument, not a proof, for a
stopping that turns inside the reach; a support that falls outside its planned nodes is caught at run time
(`ShiftedFoldNodes`), not silently wrong.

**Decided by the director.** Left as they are and listed for the owner: the beam-profile kernel reads the nominal
energy for its reconstructed-energy window and its detailed-balance norm; a user-defined resolution kernel is built at
the nominal energy; the support of a Gaussian beam spread with an energy-dependent width is 5 σ of the nominal energy
(a tail effect of about 1e-7). None reads the stopping. Whether the first two should follow the rule depends on what
the shift means for those data (a beam calibration, or a calibration of the reconstructed energy): the owner's and
the PI's call.

**T2's fits moved**: χ² 297.568 → 297.860 by the three routes, the fitted shifts by less than 0.1 σ; the gradient
route needs 117 calls where it needed 454, MIGRAD without a gradient converges in 117 calls where it stopped at the
call limit of 200 before.

## 2026-11-13 (agent, branch `feature/upgrade-5-true-energy-stopping`) — an energy shift is the true beam energy

**Done.** The PI's rule of 2026-11-12: the fold of a point measured at E with the shift s is the fold of the same point
without a shift at E + s.  Stated once in `include/ShiftedFold.h` (head of the comment, questions 1, 3, 5) and in
`energy_grid.md` section 2.5.  Before, the target kernel of a shifted point read eps(E_n - Delta) (`TranslatedStopping`,
deleted) and the thickness and window of the nominal energy.  Now:
- `MoveTargetToShift` (ShiftedFold.cpp), called by the shift step (`Segment::UpdatePointEnergyWithShift`), sets the
  point's target lab energy and thickness at the true energy (`TargetEffect::KernelThickness`) and
  `Unit::SupportThicknessChange`, dT, the growth of the support's target thickness (`TargetEffect::SupportThickness`,
the
  one formula the set-up also lays the support out with; `SubGridSpec::supportThickness` holds the set-up's value);
- the moved support is [a + Delta - dT, b + Delta] (`MovedSupport`); the kernel reads each node's own stopping
  (`PointFold.cpp`);
- the planner's support of a shifted point is widened below by G, the largest dT over the reach plus the spread of
  three samples (`PlannedShiftSupport`, `SupportThicknessChange`; the bound is derived in the header).  dT = 0 bit for
  bit for a constant stopping and for kernels without a target.
- The gradient's shift column (a central difference of the fold on the fixed nodes) carries the kernel's E_0 dependence,
  since the move goes through `MoveTargetToShift`.

**Tests.** New `tests/shift_true_energy` (identity of 36 shifted data sets with twins at the shifted energies, one
mesh 0 / tol 1e-9, two meshes 1.04e-6 / tol 3e-6; gradient of a fitted offset on a 72 keV layer against a central
difference of chi^2, 3.9e-7 / tol 1e-5).  `tests/shift_sqrt_subpoints` check (g): the analytic plateau on narrow_sqrt
(2.06e-4 / tol 5e-4; refined x8 8.3e-5) with the new model `narrow_sqrt_nores`.  Injected errors shown failing (README
of shift_true_energy).  `tests/upgrade/predictions.txt` regenerated for the new projects (additions only).

**Numbers that moved** (all projects x 3 bases against `tab-gridm-1`, tag `truestop-full1`): only T2 (T2.azr 1.2e-6,
quadrature of the wider planned support; T2_true.azr 0.21 %, chi2 330.425 -> 330.835) and shift_sqrt_subpoints const,
narrow_const, narrow_sqrt, shiftfn, sqrt (-0.89 to -0.97 %); new files of the new models.  T2 fits: chi2 297.568 ->
297.860, shifts within 0.1 sigma, fewer evaluations.  Instruction counts within 0.2 % (T2_true -1.4 %).

**Re-baselined** (commit of its own): `tests/upgrade/golden/T2/calc_truth` (AZUREOut R=2, chiSquared) and
`tests/convert_params/expected/chi2_renamed` 293.675 -> 293.977, both T2 at non-zero shifts.

**Verified by the agent** (build `work/b-truestop`, 0 warnings): all projects x 3 bases against `tab-gridm-1` (tag
`truestop-full1`, 805 / 778 / 805 identical; the differing files are those of T2 and shift_sqrt_subpoints named above,
plus the new projects' files); converted-tree gate `gtruestop1` 613 / 586 / 613 identical, 15 per basis differ, all of
T2 and shift_sqrt_subpoints; `ctest -LE slow` 219/220 on the first commit (`upgrade_predictions`, regenerated in the
second commit, then passes); `tests/run_tests.sh` 159/161 before the re-baselining (convert_params, upgrade T1-T3's
calc_truth; both pass after it); 1 against 3 threads identical (T1, T2, T2_true, T3, 13N, narrow_sqrt, shiftfn, the
identity model; 3 bases); fits (T2 by lm, gradient, migrad 200; 13N, layer and stopping yields: identical);
`perf stat -e instructions:u` within 0.2 % except T2_true (-1.4 %).

**Open.** The beam profile reads the nominal energy for its window and detailed-balance norm; the Gaussian beam spread
with an energy-dependent sigma has its support clipped at 5 sigma of the nominal energy; a UDR kernel is built at the
nominal energy.  None depends on the stopping; each is a candidate for the same rule (the director's call).  A thin
target whose set-up support was clipped at the 1 keV floor keeps the floor clip (header, question 3).

## 2026-11-13 (director, 01:50) — merged: Phase 6, second branch (no process-wide pointers or memos; the solvers' buffers non-const)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (nothing pushed): `feature/upgrade-6-evaluation-2` (the
agent's entry below says what was done, by file). No output file of a project changed.

**What it is.** `g_config` is gone (its three definitions, 30 copies in test programs, `Kinematics::Relativistic`):
a pair is built with the run's settings of its functions (`PairFunctionSettings` in `include/PPair.h`: relativistic
kinematics, the Coulomb table and its file), a required constructor argument, and every former reader takes them
from its pair. The Python binding no longer swaps a global `Config` and a global amplitude cache per call; the
external-capture amplitudes, the external-capture term table and Park's shift-derivative memo belong to the run
(`include/RunScratch.h`). No `mutable` is left in `AMatrixFunc.h`, `LevelMatrixSolver.h`, `ThmVertexBoundary.h`.
`tests/pyazr/concurrent_sessions_test.py` (slow): pairs of sessions evaluating at the same moment on two Python
threads, and two Coulomb table files in one process, each give the bytes of the session alone in a fresh process;
on the engine before the branch 20 of 40 and 36 of 40 evaluations differed.

**Not done by the branch** (its brief asked for them): the pooled copies of a fit and of the MCMC, and the
diagnostics record. "No hidden state" is not reached and the roadmap's box is not ticked: 88 pending rows of
`tests/consistency/allowed_state.tsv` (milestone 2: 4, milestone 3: 11, milestone 4: 44, milestone 5: 29).

**Verified by the director** (own build `work/b-eval2m` of the branch head, reference `tab-gridm-1`, base build
`work/b-head3`).
Own build of the branch head (`work/b-eval2m`, 0 warnings). Whole tree × 3 bases (tag `tab-eval2m-1`) against
`tab-gridm-1`:
827 / 800 / 827 files identical. The eight fits (T2 by three routes, 13N by two, layer yield, stopping yield by two):
identical file by file to the build of the merged head before (`work/b-head3`). `ctest -LE slow` 220 of 220, the slow
tests 11 of 11 (`pyazr_concurrent_sessions` 515 s), `run_tests.sh` 160 of 160, spec equivalence 181 equal, the
converted-tree gate 628 / 601 / 628 identical, calculations at 1 against 3 threads identical (T1, T2, T3, 13N × 3
bases). Not run by the director: the agent's instruction counts (its table: within 0.1 % of the base).

**Reviewed.** `PairFunctionSettings` and its readers (`Kinematics::Channel`, `CoulFunc`, `ShftFunc`, the key of the
shift-derivative memo); `RunScratch`'s three shared members and their locks; the binding's diff (every `ConfigScope`
removed, `InitializeECAmplitudeCache` and `InvalidateShiftFunctions` given the session's `Config`). Notes: (1) the
settings are copied into the pairs at `CNuc::Fill`; nothing changes them afterwards (the three are set from the
command line, the `<config>` block or the session's options before the set-up; searched); (2) one expected line of
`tests/upgrade/predictions.txt` changed: `--print-predictions` ran before `main` set the global and so ignored
`relativistic=1`; it now refuses the 44.9 degree point of `li7pn_edge_rel` as the calculation does. Accepted: it is
the correction of a defect, and no output file moves.

**Decided by the director.** The bytes of input files stay in one store per process (`src/InputFiles.cpp`): entries
are keyed by absolute path and file stamp, immutable once read, behind a lock, and `cache_stats()` is documented per
process; its rows of the allow-list become permanent in the next branch.

**Found by the agent, in work in the next branch** (`feature/upgrade-6-evaluation-3`): the engine resolves a
project's relative paths against the process's working directory, which pyazr changes per call, so two sessions in
two directories on two threads race (the new test uses absolute paths); the signature of `intEC.dat` did not hold
the Coulomb table. A third finding of that branch, for a later one: the pooled copies cannot be deleted without
splitting the `Site` layer into problem and per-evaluation values (`evaluation.md` section 4.1 there).

**Decided by the PI since the entry of 22:55 below** (its items 4 and 5): the stopping of a shifted, folded point is
to be evaluated at the shifted energy (in work, `feature/upgrade-5-true-energy-stopping`); the radiative widths
stay out of the level's total width: the engine's resonance strength, omega Gamma_gamma Gamma_p over the sum of the
particle widths, is the convention of AZURE2 and is not to be changed or made optional. Also in work at the PI's
wish: the total radiative width, the gamma branchings and the mixing ratios of a level as fit parameters
(`feature/upgrade-8-gamma-branchings`).

## 2026-11-13 (agent, Phase 6 branch 3 continued) — the Coulomb table named by content; relative paths without the working directory; the branch's gates

**Done.**  (1) The `intEC.dat` signature names the Coulomb table by what it holds (`CoulombTable::ContentIdentity`:
header and a 64-bit hash of every value, computed once per loaded table on first request behind `std::call_once`),
not by its path and time: a copy elsewhere reuses the file, other values refuse it (`tests/ec_table_signature`).
Cost 1.2e8 instructions per run that sets up external capture (T1 +0.41 %, 13N +0.16 %), none otherwise.
(2) A project's relative paths are resolved against `Config::baseDirectory` (the project file's directory, or the
pyazr session's `cwd`) through one function (`azure2::experiment::PathToOpen`, `Config::OpenPath`) at every opening
site; messages print the project's own string; pyazr's `_in_dir` and `_Engine` are deleted (`evaluation.md` 7.5;
`tests/run_from_elsewhere`, `concurrent_sessions_test.py` pair R S).  (3) Merged `feature/upgrade-4-data` at
2ed0dba9 and e4fb58f1.

**Verified** on the merged head (build `work/b-eval3`, 0 warnings, 3 threads): whole tree x 3 bases against
`tab-gammam-1` (tag `tab-eval3-j`): 887 / 851 / 887 files, every number and every `run.log` equal, the 23
`intEC.dat.sig` per basis differ (the new signature); before the merge a second whole run with the table copied
elsewhere (`tab-eval3-g` against `tab-eval3-f`) gave the same bytes but for the `run.log` line naming the table.
Converted-tree gate `geval3b` (REF `tab-gammam-1`): 679 / 643 / 679, only the signature files differ.  `ctest -LE
slow` 221 / 221; `run_tests.sh` 164 / 165 at first (`layer_yield`: its refusal cases ran a project written outside
its data's directory and relied on the working directory; the fixture now copies the project beside the data:
165 / 165 by rerunning that test on this build and on `b-gammam`); `concurrent_sessions_test.py` PASS (pairs A B,
N M, A T, R S, 40 evaluations each; the R S pair with the base's pyazr and engine fails: a session cannot set up);
1 against 3 threads identical (T1, T2, T3, 13N x 3 bases); the eight fits identical to `b-gammam` file by file
but 13N's `intEC.dat.sig`.  Instructions against `b-gammam` (1 thread): T2 -0.001 %, T1 +0.41 %, 13N +0.16 %,
18O_p_a_thm +0.06 %, natural_neon +0.001 %, T2 gradient fit +0.36 %, 13N fit +0.02 %.  Memory against `b-eval2`
before the merges (`work/eval3/mem_t1.txt`, `mem_t3.txt`).

**Not done.**  The diagnostics record (`evaluation.md` section 2); `AZUREAPI`'s `mutable` members.
**Open for the director.**  `--validate`, `--describe`, `pyazr.validate`, `pyazr.describe` still read data files
relative to the process's working directory (held to one text by `tests/pyazr/validate_test.py`): moving them changes
what public pyazr functions return.

## 2026-11-13 (agent, Phase 6 branch 3, `feature/upgrade-6-evaluation-3`) — copies made when needed; the Coulomb table in the intEC signature

**Done.**  (1) `evaluation.md` 4.1: what an evaluation of a fit writes, member by member.  Finding: for shifted
terms, fitted target-effect parameters and re-placed grids the hand-over rewrites most of a `Site` (energies, angles,
channel values, Coulomb amplitude, rows, stopping, thickness, slot layout), so removing `WorkingCopy::Clone` needs a
split of the `Site` layer, not only of `SegmentState` and the overflow nodes.  (2) The pools of `AZURECalc` and
`AZURECalcMCMC` start empty and clone a copy only when no idle one is left (was `max(4, threads)` for a fit, twice
the hardware threads for the MCMC); `InitializePools`, `pools_initialized_` deleted.  (3) The `intEC.dat` signature
holds the Coulomb table's file (path, size, modification time) or "exact Coulomb functions"
(`tests/ec_table_signature`, fails on the base).  (4) `tests/describe_relativistic`: `--validate`/`--describe` read
`relativistic` (already true on the base).  (5) The input-file byte store's four allow-list rows are permanent.
Allow-list pending rows 88 -> 82 (milestone 2: 4 -> 0; 4: 44 -> 42; 3: 11; 5: 29).

**Verified** (build `work/b-eval3`, 0 warnings, one thread): T1, T2, T3, 13N, natural_neon, 18O_p_a_thm x 3 bases
against `tab-gridm-1`: 96 / 96 / 96 files, every number equal, only the 3 `intEC.dat.sig` per basis differ (13N, T1,
T1_start: the new signature); ctest `upgrade_history`, `pyazr_session_isolation`, `repo_consistency` pass.

**Not done** (the thread budget stayed at 1): the whole tree, the eight fits, `AZURE_TRACE_FCN=fresh`, the memory
table, instruction counts, `ctest -LE slow`, `run_tests.sh`, 1 against 3 threads, the converted-tree gate.  Not
started: relative paths without the working directory (plan in `evaluation.md` 7.5), the diagnostics record,
`AZUREAPI`'s `mutable` members.

## 2026-11-12 (director, 22:55) — merged: one grid (the key `useAdaptiveGrid` and the uniform grid of a point's own removed)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (nothing pushed): `feature/upgrade-5-one-grid` (the
agent's entry below says what was done, by file). The PI allowed the removal on 2026-11-12 on the condition that
narrow resonances stay covered.

**What it is.** Every folded point takes its nodes from the shared mesh of its nucleus; the second path (a uniform
grid of the point's own, moved rigidly with an energy shift) is deleted with its readers. `SetUpTargetEffects` lays
out the base term and the components with one function. The reaction rate has one path, the tabulated grid.
`<config> useAdaptiveGrid` is refused by the reader with `E_CONFIG_KEY_REMOVED`; the converter drops it from an old
file and says what changes when it was 0.

**What changes for a user.** Projects with the key at 1 (the default): no number. Projects with the key at 0: folded
values move to the converged fold. In the tree that is `examples/c12c12_tumino2018`, χ² 111.661 → 117.516.

**Verified by the director** (the branch merged onto the first Phase 6 branch in `rehearsal/trial38`, build
`work/b-gridm`; references `tab-eval-d1` and the agent's `tab-grid-1`).
Whole tree × 3 bases (tag `tab-gridm-1`, 827 / 800 / 827 files): against `tab-eval-d1` no number differs; the 17
differing files per basis are the 8 outputs of the three removed models (`shift_sqrt_subpoints/models/dense`,
`reaction_rate/models/n12C_pert`, `p12C_2g_pert`) and 9 `run.log` whose `E_DATA_FILE` messages name a line one lower
(the `useAdaptiveGrid` line is gone from the files); against the agent's `tab-grid-1` every file but one `run.log`
(the Coulomb-table line of `thm_coherent_background/toy`, from the Phase 6 branch). The eight fits (T2 by three
routes, 13N by two, layer yield, stopping yield by two): identical file by file to the build of the Phase 6 branch.
`ctest -LE slow` 217 of 217, the slow tests 10 of 10, `run_tests.sh` 160 of 160, spec equivalence 181 equal, the
converted-tree gate 628 / 601 / 628 identical, calculations at 1 against 3 threads identical (T1, T2, T3, 13N × 3
bases). Build: 0 warnings. By hand: `examples/c12c12_tumino2018` χ² on the engine before the branch with the key at 0:
111.660780 (117.328303 at 8 times the points, 117.515866 at 64 times), with the key at 1: 117.515816932; after the
branch 117.515816932 (117.515832232 refined eightfold), the bits of the key at 1. `narrow_sqrt`, first point: after the
branch 2.1802211257e-06 (refined 8 and 64 times 2.1798404627e-06, 2.1798336905e-06), the bits of the old key at 1;
the old uniform grid 2.06e-08 (5.09e-07 at 64 times the points). A file with the line: exit 1 with
`E_CONFIG_KEY_REMOVED` and the line number; `--convert` of the old `c12c12_tumino2018`: exit 0, the note printed, the
key absent. **Against the analytic narrow-resonance yield** (the PI's question; `narrow_sqrt` with a constant stopping
and no shift, the non-resonant part 1.42e-9 b taken from a run with the level's proton width at zero, refined 64
times): (λ_r²/2) ω Γ_γ (m+M)/(M ε N) with the arctan edge factor gives 2.160633e-06 b against the engine's
2.160405e-06 (1.1e-4) for Γ_p = 9.749 eV, and 2.161543e-06 against 2.161505e-06 (2e-5) for Γ_p = Γ_γ = 38.5 meV, where
the default mesh and the refined one differ by 1.9e-4. See items (4) and (5) below for what the comparison needed.

**Reviewed.** The diff of `TargetSetup.cpp` line by line against the two loops it replaces (same statements, the
component's pairs and factor as arguments; the error text of a non-positive thickness now names the component in both
cases); `ShiftFrames.cpp`, `Site.h`, `GradEnergyTasks.cpp`, `WorkingCopyUnits.cpp` (the readers of a row set without
stored values are gone with the own grid); `ReactionRate.cpp`; the reader's diagnostic and the converter's note. "Never
coarser than the uniform grid of the same number of points" holds per support: the mesh's smooth step is
window / ceil(window / uniform step) (`MeshPlanner.h`).

**Open, for the PI.** (1) The per-temperature quadrature of the rate is now reached only when fewer than three grid
points have σ > 0, and no test reaches it (the two test models that used it selected it with the key): keep it as an
untested fallback, give it a test, or replace it by a refusal. (2) A folded point whose pairs give no excitation
energy (a β-delayed entrance without an exit pair) is refused at set-up; no project of the tree has one. (3) A
format-3 file that still has the line is refused, not silently accepted: deliberate (format freedom), but a user of
the last weeks' files must delete one line. (4) Found with the narrow-resonance check, not caused by this branch (the
engine before it gives the same bits): with an energy shift on a folded data set the stopping is evaluated at the
nominal energies, not the shifted ones. On `narrow_sqrt` (a shift of 10 keV at 272 keV) the yield is 0.90 to 0.97 %
above the same true beam energy without a shift, equal to ε(E_r)/ε(E_r − shift) to five digits at the four points; in
general (dlnε/dlnE) × shift/E. A correction changes numbers of shifted, folded data sets with an energy-dependent
stopping: the PI decides. Not checked: whether `origin/dev` does the same. (5) The same check shows what the engine's
resonance strength is: ω Γ_γ Γ_p / (sum of the particle widths), because radiative widths do not enter the level
matrix (`docs/source/theory/conventions.rst`; the loop is the same in the original engine; a warning is printed above
Γ_γ = 10 % of the particle width). 0.39 % above the textbook ωγ for the test level, a factor 2 for Γ_p = Γ_γ.

## 2026-11-12 (agent, branch `feature/upgrade-5-one-grid`) — one grid: `useAdaptiveGrid` removed, every folded point on the mesh

**Done** (the PI's decision of the evening: the key may go, provided narrow resonances are covered by the mesh).
- Engine: the uniform grid of a point's own is gone with everything only it used: the own-grid branches of
  `SetUpTargetEffects` (its two near-copy loops, base segment and components, are one function `layOutTerm` now; a
  folded point whose pairs give no E_x frame is refused there, it threw in `ReplaceSubPointGrids` before),
  `AddGridSlots`, `Segment::MoveOwnGridSlots`, `Unit::CalcNodeRows`, the readers' branches for a row set with
  `valueOffset` −1 (`Site.h`, `WorkingCopyUnits.cpp`, `GradEnergyTasks.cpp`), the key's early return in
  `ReplaceSubPointGrids`, `Config::useAdaptiveGrid`; dead after it: `AdaptiveIntegrationGrid::GenerateGrid(CNuc*)`,
  `GenerateGridWithInfo`, `GetExpectedPointCount`, `SetConfig`, `GridPoint`.
- Rate: always the tabulated grid; the per-temperature quadrature (`RunPerTemperatureRateCalculation`) is only the
  fallback of a grid with fewer than three points of σ > 0. `tests/reaction_rate` (c) and (h) and their models went.
- Format: a file with the line is refused (`E_CONFIG_KEY_REMOVED`, "delete the line" and why); `io::ConfigKeys()`
  no longer lists it; the converter drops it and notes what changes when it was 0. The line removed from every project
  of the tree, the converter's expected outputs and the generators (the old originals in `tests/convert/fixtures` keep
it).
- Tests: `shift_sqrt_subpoints` without `dense`; (c) `sqrt` against itself refined eightfold (9.3e-5), (f)
  `narrow_sqrt` against itself refined eightfold (5.4e-4, measured 1.78e-4): fails at 0.98 with the level anchors
  taken out of the mesh plan.
- Docs: `ShiftedFold.h`, `energy_grid.md` (status), `project_format3.md` §3.2, `FOR_THE_OWNER.md`, `ROADMAP.md`,
  `qt_parity.md`, `project.md`, the `azure2-calculate` skill, `project_from_python.rst`.

**Verified** (build `work/b-grid` of the last code commit). Whole tree × 3 bases against `tab-thm5` (tag `tab-grid-1`):
every output file of a project with the key on identical; run.log of `electron_screening`, `target_density`,
`mcmc_physical_norm/13N_fast` differ only in the line numbers of their data-file errors (one line fewer); files 835 /
808 / 835 → 827 / 800 / 827 (the 6 files of `shift_sqrt_subpoints/models/dense` and the run.log of the two deleted
rate models). `ctest -LE slow` 215/215 (without `pyazr_park_gradient` and `gsl_quiet`, which set 4 threads). 1 against
3 threads, T1, T2, T3, 13N × 3 bases: outputs identical (run.log: progress bars, Coulomb lookup counts). Eight fits
identical to `b-thm5`. `tests/reaction_rate`: rates and `cross_section_grid.dat` of p12C, n12C, p12C_2g, moved_azr
byte-identical to `b-thm5`. `examples/c12c12_tumino2018` (the key at 0): with the key at 1 on `b-thm5` byte-identical
to this branch; against this branch refined eightfold the old uniform grid was off by up to 57 % in a folded value,
the mesh by 2.7e-5; χ² 111.661 → 117.516 (refined 117.516). Converted-tree gate with `tab-grid-1` as its reference:
628 / 601 / 628 files identical, none classified, differing or missing. `tests/run_tests.sh`: 158 of 159 passed; the one
failure, `cli_help` ("--answers is accepted but --help does not print it"), passed when run again alone, on this build
and on `b-thm5` (no CLI code is touched); `thm_dw_vertex` was not run (it sets 4 threads, above the budget of 3).

**Open**: the per-temperature quadrature of the rate has no test now (no model reaches its fallback). The refusal of a
folded point without an E_x frame (a beta-delayed spectrum without an exit pair) is reached by no project.

## 2026-11-12 (agent, Phase 6 branch 2) — no process-wide pointers or memos; the solvers' buffers non-const (`feature/upgrade-6-evaluation-2`)

**Done** (`evaluation.md` section 7.4).  `g_config` is deleted: each particle pair is built with the run's settings of
its functions (`PPair::FunctionSettings`: relativistic kinematics, use and file of the Coulomb table), and every
former reader takes them from its pair.  The binding's `ConfigScope` swap is deleted.  The EC term table, Park's
shift-derivative memo and the `--no-ec-table` amplitude store are members of the run's scratch (`RunScratch`).  No
`mutable` in `AMatrixFunc.h`, `LevelMatrixSolver.h`, `ThmVertexBoundary.h` (the writers are non-`const`; the
`const_cast` of `GetAMatrix` is gone).  New test `tests/pyazr/concurrent_sessions_test.py` (ctest
`pyazr_concurrent_sessions`, label slow, ~11 min): fails on the base (20/40 and 36/40 evaluations differ, one pair
crashes at set-up; the two-table sequences give other bytes), passes here.  Merged `feature/upgrade-4-data` (3c6d0e9b).
Allow-list 107 -> 88 pending rows.  One expected line changed: `tests/upgrade/predictions.txt` for li7pn_edge_rel
(`--print-predictions` read "non-relativistic" through the null g_config; now refuses the point as the calculation
does).
**Verified** (build `work/b-eval2` of the merged head): whole tree x 3 bases against `tab-gridm-1` (tag `tab-eval2-a`):
827 / 800 / 827 identical, run.log included; T1, T2, T3, 13N x 3 bases at 1 and 3 threads identical;
`ctest -LE slow` 219 of 220 before the predictions line was regenerated (the one failure was that line; it passes
after); `tests/run_tests.sh` 160 passed, 0 failed; the eight fits identical file by file to `work/b-head3`; the
converted-tree gate (`geval2-a`) 628 / 601 / 628 identical; user-space instructions against `b-head3` at 1 thread:
T2 +0.03 %, 13N -0.0002 %, natural_neon +0.10 %, 18O_p_a_thm +0.001 %, T2 fit +0.06 %, 13N fit +0.03 %.
**Not done**: the input-file bytes (a decision: `cache_stats()` counts the process); part 3 (the pooled working copies),
part 4 (diagnostics record).

## 2026-11-12 (director, 22:13) — merged: Phase 6, first branch (the engine's state surveyed, the per-thread caches and warn-once flags owned by the caller, a check against new state)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (nothing pushed): `feature/upgrade-6-evaluation` (the
agent's two entries below say what was done, by file). No output number of a project changed.

**What it is.** `docs/upgrade/evaluation.md`: every object that outlives an evaluation, the 44 environment switches,
what an evaluation costs, the target (an immutable problem, a parameter vector in, a `Results` out) and five
milestones. In code: no `thread_local` is left (26) and no warn-once static (12); what was in them is scratch of the
energy tasks (`include/TaskScratch.h`, in the caller's `Results`), scratch of the run (`include/RunScratch.h`,
`Config::scratch`) or the run's once-only notices (`Config::FirstTime`), or was deleted where it bought nothing.
Check 38 refuses a new `static`, `thread_local`, namespace-scope variable or `mutable` member outside
`tests/consistency/allowed_state.tsv`; `tests/pyazr/session_isolation_test.py` runs sessions of three projects (and
one project at two channel radii) one after another in one process, in both orders, at 1 and 3 threads.

**"No hidden state" is not reached and the roadmap's box is not ticked.** The allow-list holds 10 permanent rows
(tables, locks of output streams, GSL's error handler) and 107 pending ones. Of the pending, those that hold values
of a calculation process-wide: `g_config` and `g_ecAmplitudeCache` (swapped per call by the binding: two sessions
evaluating at the same time on two threads race; by reading, no test), the shift-derivative memo of `CNuc.cpp`, the
external-capture term table of `ECTable.cpp`, the input-file bytes. The keys of the last two memos hold
`useCoulombTable` but not the table's path (two sessions with two table files would share entries; by reading). The
rest: the fitter's and sampler's pools and `mutable` members, the solvers' `mutable` buffers, diagnostic counters,
numcmc's callback globals.

**Verified by the director** (build `work/b-eval` of the branch head, reference `tab-thm5`).
Whole tree × 3 bases (tag `tab-eval-d1`): every output file identical to `tab-thm5` (835 / 808 / 835). One `run.log`
per basis differs, `thm_coherent_background/toy`: the line "Coulomb table: … (l = 0..10, …)" is no longer printed; the
old shift memo of `ThmVertexBoundary` loaded the table only to form its key, and that project never looks the table
up. The eight fits (T2 by three routes, 13N by two, layer yield, stopping yield by two): identical file by file to the
base build. `ctest -LE slow` 217 of 217 (with `pyazr_session_isolation`), the slow tests 10 of 10, `run_tests.sh` 160
of 160, spec equivalence 184 equal, the converted-tree gate 633 / 606 / 633 identical, calculations at 1 against 3
threads identical (T1, T2, T3, 13N × 3 bases). Build: 0 warnings. After the merge with the third part of the transfer
library, `repo_consistency_test.py` (with check 38 over the new `physics/transfer` files) passes on the merged tree;
the merged tree itself was not rebuilt here (the two branches share no source file; it is built with the next branch).

**Reviewed.** `TaskScratch`, `RunScratch`, the `Config` and `Results` changes and every call site of the two
`ForThisThread`; the diffs of the passes, the gradient, the hand-over, `ECIntegral`, `CoulFunc`, `ThmVertexBoundary`;
the scanner and the allow-list. Notes: (1) `TaskScratchSet` is indexed by the OpenMP thread number, which is safe
because every concurrent evaluation owns its `Results` (true at the call sites read); (2) `RunScratch` is a map per
OS thread behind a mutex taken once per `ECIntegral` object; its contention above 3 threads is not measured; (3)
`Config` now carries shared mutable run state (`printedOnce`, `scratch`): an interim home until a run object exists
(milestone 2); (4) the THM shift memo's key lost the Coulomb table pointer, correct because the memo now belongs to
one run; (5) the scanner was tried on planted statics in lambdas, templates, nested blocks and class members: all
found (a static inside a lambda without a parameter list is reported under the lambda's name).

**Corrected by the director in the documents.** The owner note said that two sessions "never share a value"; they
share the process-wide memos listed above, and two sessions at the same time are not safe yet. Rewritten.

**Decided.** The continuation (`work/briefs/BO_phase6_no_pointers_no_pools.md`): the process-wide pointers and memos,
the solvers' `mutable` buffers, then an evaluation of a fit that writes only its `Results` (the pools of working
copies go), then a diagnostics record read once, with every switch keeping its name and effect.

**For the PI** (not decided by the director; several of the switches are in the user manual): whether the
number-changing environment switches (`AZURE2_COULOMB_TABLE`, `AZURE2_SRIM_DATA_DIR`, `AZURE_OPTIONS_FILE`,
`AZURE_REFINE_QUADRATURE`) become command-line options shown by `--describe`; whether ten switches no test uses are
deleted (`AZURE2_SWAP_COUPLING_ORDER`, `AZURE_BOUND_SHIFT_NUMERICAL`, `AZURE_NO_NEUTRAL_CLOSED_FORM`,
`AZURE_SHIFT_DE_1EV`, `AZURE2_EC_TABLE_NODES`, `AZURE_GRAD_DENSE`, `AZURE_TE_WEIGHTS_FD`, `AZR_DEBUG_GRID`,
`AZURE_GRAD_DEBUG`, `AZURE_CHANNEL_DEBUG`; eight are named in `docs/source`); whether one `--threads N` replaces
`OMP_NUM_THREADS` as the stated budget.

**Not run by anyone**: ThreadSanitizer or AddressSanitizer (the compiler here has none usable); a test that a
once-only warning is printed again by a second run in one process.

## 2026-11-12 (agent, Phase 6 branch 1, part 2) — no thread_local, caller-owned scratch, check 38 (`feature/upgrade-6-evaluation`)

**Done** (`evaluation.md` section 7.3 by file).  Every `thread_local` is gone: the energy tasks' solvers and the
channel-space values are `TaskScratch` in the caller's `Results`; the external-capture integrals' radial values,
Whittaker values, workspace and rate grids and the THM vertex boundary's S_c are the run's per-thread scratch
(`RunScratch`, `Config::scratch`); the THM level-width, Coulomb-term and angular-coupling memos and the per-thread
Coulomb-table lookup are deleted (measured: they bought nothing).  Warnings printed once per process print once per
run (`Config::FirstTime`).  Check 38 refuses new state outside `tests/consistency/allowed_state.tsv` (117 rows: 10
tables, 107 pending with their milestone).  `tests/pyazr/session_isolation_test.py` (A B A in one process, 1 and 3
threads) passes on the base and here with the same digests, fails with an injected process-wide memo.
**Verified** (build of 41896d9f, `work/eval/bin/AZURE2-b1`): every project × 3 bases against `tab-thm5`
(`tab-eval-b1`): 0 output files differ; one `run.log` differs (thm_coherent_background/toy: the line "Coulomb
table: … loaded" is no longer printed, because the deleted THM memos loaded the table only to put its address in
their keys; the project's numbers do not read the table).  `ctest -LE slow` 217/217.  1 against 3 threads (T1, T2,
T3, 13N × 3 bases) identical.  The eight fits identical to the base file by file.  Instructions within 2 % (largest
+1.7 % on 18O's calculation).  `tests/run_tests.sh` 160 passed, 0 failed.  The converted-tree gate (`geval-b1`): 633 /
606 / 633 files identical, 0 differ.
**Open.** The pending rows (milestones 1–5 of `evaluation.md`); part 3 of the brief (the pools) not started.

## 2026-11-12 (director, 21:47) — merged: physics/transfer, third part (spin-orbit in the distorted waves, the ³He and ¹⁴N Trojan horses, the ⁶Li bins against Frescox)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (nothing pushed): `feature/upgrade-7-transfer-3` (the
agent's four entries below say what was done, by file). A leaf library: nothing of the engine calls it, the diff
outside `physics/` and `docs/` is empty, no output number of a project can change.

**What it is.** The kernel integrals are indexed by pairs of partial waves (l, j); one function does the angular and
spin sums with or without spin-orbit, in finite and zero range (`TransferAmplitudes.cpp`, `transfer.md` §19), and
gives the cross section and the projectile's iT11. Against FRESCO 3.5: ¹⁵N(d,n) with spin-orbit 2e-4; ¹⁸O(³He,d)¹⁹F
(two states, finite and zero range, with and without spin-orbit) ≤ 8.9e-4; ¹²C(¹⁴N,d)²⁴Mg to two ¹²C + ¹²C states,
finite range, ≤ 8.3e-4 (asked: 1 %).

**Verified by the director** (fresh build `work/b-dwba3v` of the branch head).
Fresh build of the branch head (`work/b-dwba3v`): 0 warnings. `ctest -R physics`: 34 of 34.
`repo_consistency_test.py`: consistent. The three new tests run by hand print the accuracies the agent reported:
¹⁵N(d,n) with spin-orbit 2.0e-4 (iT11 1.0e-4), zero range 6e-6; ¹⁸O(³He,d)¹⁹F worst 8.9e-4 (5/2⁺ state, finite range),
FRESCO's own reference against its refined run up to 1.8e-3 for the ground state; ¹²C(¹⁴N,d)²⁴Mg 8.3e-4 and 3.3e-4
(FRESCO's own refinement 1.0e-3 and 3.1e-4); without spin-orbit the earlier cross sections are reproduced to 7e-16.
Four benchmark cases regenerated with the director's own runs of FRESCO 3.5 and Frescox 7.4
(`he3d_gs_finite_range_so`, `n14d_2p_finite_range`, `dn_finite_range_so`, `li6d_bin_wide_frx`): identical to the
committed files. Not run by the director: the whole tree, the fits, the gate (no file they read is touched); the
narrow ⁶Li bin with Frescox (the agent did not run it either).

**Reviewed.** `TransferAmplitudes.{h,cpp}` read: the index arithmetic of the projection sets, the vertex couplings
of (19.6), the normalisation of iT11 (√3/2 A_y for spin 1, A_y/√2 for spin ½) and the trace Tr(T S_y T†); no mutable
state. The FRESCO comparisons, with the reduction test (waves flagged with a vanishing spin-orbit term against waves
without, 7e-16), are the evidence for the algebra; it was not re-derived by the director.

**Limits to keep in mind.** FRESCO's own step refinement moves its ³He finite-range results by up to 1.8e-3 (0.2 %
was asked of it); the ¹⁴N cases use test potentials and unsymmetrised ¹²C + ¹²C (the exchange term is described in
`transfer.md` §21, not implemented); zero range is meaningless for ¹⁴N; two ANC checks have their own stated
tolerances (3e-5 for p + d, 1e-4 for ¹²C + ¹²C). The transferred cluster is one overlap component (no deuteron
d-state); spins above 1 are refused by the distorted waves.

**The ⁶Li bins (the PI's request, with Frescox 7.4, approved by the PI on 2026-11-12).** Frescox was installed under
`opt/frescox` (outside the repository; built with the system gfortran 11.5, as FRESCO 3.5 was, because the conda
environment has no Fortran compiler; accepted). It reproduces FRESCO 3.5 to 5e-4 on the wide ⁶Li bin and has the same
loss of accuracy at high L, as expected of the same lineage. At a cut of 12 fm our extrapolation in the step is within
0.004–0.4 % of Frescox's refined run from 0° to 60° and 1.2 % at 90°; neither code is converged to 0.2 % on a bin cut
where it is not small. So 1 % is not shown for the ⁶Li bins and the remaining 1 % at 90° is not attributed. The narrow
⁶Li bin was not run. The agent's proposal, a bin whose envelope has a node at the cut radius (no first-order step
error in either code), is the next benchmark object if this is pursued.

**Not done.** ³He and ¹⁴N to the continuum (item 4 of the brief).

**Open** (carried). The wiring into the THM observable (`transfer.md` §12, after the Phase 6 branches); a spin-dependent
THM weight needs the amplitude table public (`AmplitudesAt`); the engine's `ThmDistortion` is spinless.

## 2026-11-12 (agent, branch `feature/upgrade-7-transfer-3`) — physics/transfer, third part, item 3b: the ⁶Li bins against Frescox

**Done** (`transfer.md` §23): Frescox 7.4 cloned from github.com/LLNL/Frescox (eca88f3) to `opt/frescox-src`, built
with the system gfortran 11.5 (3 jobs; the conda environment has no Fortran compiler), binary `opt/frescox/bin/frescox`,
README beside it; five distribution test jobs equal its gfortran reference outputs. Generator: `--frescox` and
cases marked `code="frescox"`; data `li6d_bin_wide_frx` (not compared by a test).

**Found**: Frescox reproduces FRESCO 3.5 to 5e-4 on part 2's ⁶Li bin and shares its ERWIN accuracy loss (L = 26–35);
at RSP 20–30 fm it does not converge (results change by factors with rmatch), so it cannot settle the RSP question.
At RSP 12 fm neither code is converged to 0.2 % (Frescox 0.3–0.6 % in its step; ours first order in h because the
bin is cut where it is not small); ours extrapolated is 0.004–0.4 % from Frescox's refined run at 0–60° and 1.2 % at
90°. Asked 1 %: not shown. Not done: a narrow ⁶Li bin with Frescox.

**End of branch** (build `work/b-dwba3`, 0 warnings): `scripts/check.sh --quick` with JOBS=1: 219/219 (the
standalone physics build included); T1, T2, T3 and 13N in three bases (`tools/runall.sh`, tag `dwba3-sub1`, a tree of
those four folders): every output and check file identical to `tab-thm5`; `git diff --stat 0800944c -- src include
api pyazr experiment convert stopping external tools tests`: empty. Not run: the whole tree, `run_tests.sh`, the gate,
fits (nothing they test is touched); item 4 (³He and ¹⁴N to the continuum).

## 2026-11-12 (agent, branch `feature/upgrade-7-transfer-3`) — physics/transfer, third part, item 3: the ¹⁴N Trojan horse

**Done** (`transfer.md` §21, §22): ¹²C(¹⁴N,d)²⁴Mg at E_lab = 30 MeV to the 0⁺ g.s. (l = 0, six nodes) and the 2⁺
(l = 2, five nodes) ¹²C + ¹²C states, finite range, unsymmetrised; FRESCO data `n14d_*` (generator: `N14_*`, a
`cutr` setting; the bound-state listing parser takes FRESCO's `*******`); new test `physics_transfer_n14_fresco`.
Zero range shown meaningless by numbers (ZR/FR 1.2–211); what symmetrising the two ¹²C would add is in the note.

**Verified**: ours against FRESCO 8.3e-4 (g.s.), 3.3e-4 (2⁺); ours converged to 2.5e-4 (0.9 and 2.8 s at one thread).
FRESCO converges only with rnl ≥ 13 fm and its integration started 12 fm inside the Coulomb turning point (cutr);
with cuts ∝ L it is 1–3 % off, without them ERWIN loses 14 digits (table in §21). The ¹²C + ¹²C ANC: FRESCO's u/W
still rises at 12 fm; the test allows 1e-4 there. Injected error: the Jacobian (1 − αβ)² instead of ³ (both cases
fail by 14 %).

## 2026-11-12 (agent, branch `feature/upgrade-7-transfer-3`) — physics/transfer, third part, item 2: the ³He Trojan horse

**Done** (`transfer.md` §20, §22): ¹⁸O(³He,d)¹⁹F at E_lab = 25 MeV to the ½⁺ ground state and the 5/2⁺ state, finite
and zero range, with and without spin-orbit in both distortions: eight FRESCO cases `he3d_*` (generator: `HE3_*`;
FRESCO needs cutl here, ERWIN at L ≥ 35 otherwise spoils its transfer), new test `physics_transfer_he3_fresco`. The
generator reads FRESCO's bound-state listing by number pattern (columns can run together); part 1's regenerated
`dn_finite_range` and `li6d_finite_range` data are identical.

**Verified**: the new test passes: ours against FRESCO ≤ 8.9e-4 (FR), ≤ 9.2e-5 (ZR), iT11 ≤ 4.3e-4. FRESCO's own
refinement moves its FR results by up to 1.8e-3 (asked 0.2 %: not reached by FRESCO; ours converged to 2e-4). The
p + d ANC differs from FRESCO's by 2.2e-5 (unresolved; the test allows 3e-5 for it). Injected error: the sign of the
S_− term of S_y (iT11 of ³He fails in three cases).

## 2026-11-12 (agent, branch `feature/upgrade-7-transfer-3`) — physics/transfer, third part, item 1: spin-orbit in the distorted waves of a transfer

**Done** (leaf library, nothing of the engine touched; `transfer.md` §19, §22):
- The transfer amplitude with (l, j)-dependent distorted waves in both channels (projectile spin J_a, ejectile spin
  J_s; target spin a spectator): the kernel integrals (6.4) per pair of waves (`TransferKernel.cpp`), the angular and
  spin sums in one function for all cases (`TransferAmplitudes.{h,cpp}`, new): without spin-orbit the j sums are the
  identity (19.7) and part 1's numbers return (8e-15). Zero range goes through the same sums.
- The projectile's vector analysing power iT11 (Madison) from the same amplitudes: `FiniteRangeObservables`,
  `ZeroRangeObservables` (`Transfer.h`).
- FRESCO data `dn_finite_range_so`, `dn_zero_range_so` (generator: `pol` lines, `drop` of unread lines); the shared
  FRESCO checks moved into `transfer_fresco_common.h`; new test `physics_transfer_spin_orbit`.

**Verified** (build `work/b-dwba3`, 0 warnings): `ctest -R physics` 32/32; the printed lines of the seven earlier
transfer tests identical before/after; new test: reduction 6.7e-16, reflection 2e-16, FRESCO FR 2.0e-4, ZR 6.2e-6,
iT11 ≤ 1.0e-4. Injected errors, each failed: j sum without spin-orbit missing j = l + s (reduction), Y_{l,|n|} for
negative n (reflection and iT11), the kernel taking the first j wave of l_i (FRESCO FR), the sign of m_i in the
entrance coupling (FRESCO FR and ZR).

**Open**: items 2 (³He), 3 (¹⁴N), 3b (⁶Li bins against Frescox), 4.

## 2026-11-12 (director, 19:05) — merged: physics/transfer, second part (transfer to the x + A continuum, FRESCO bins, prior form and remnants)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (nothing pushed): `feature/upgrade-7-transfer-2` (the
agent's entry below says what was done, by file). A leaf library: nothing of the engine calls it yet, no output number
of a project can change and none did.

**What it is in four lines.** Scattering states and FRESCO's continuum bins of x + A (`Continuum.h`); the residual
overlap is one row of a table of kinds (bound state, bin, scattering state: `Overlap.h`) and the interaction a sum of
terms (post, prior, remnants: `Interaction.h`), so the kernel has one code path; d²σ/dΩ dE to the continuum as the
narrow-bin limit; the interior part of the continuum amplitude against the R-matrix level function (`transfer.md` §16).

**Verified by the director** (rehearsal of the merge on the merged head, build `work/b-trial37`).
Unit tests 216/216 (213 and the three new ones; the 31 of `physics/` include the standalone build of the library).
Every project on three bases identical to the merged engine (`tab-thm5`: 835, 808, 835 files). `git diff` against the
base outside `physics/` and `docs/`: the consistency test only. One FRESCO case (`dn_bin_narrow_fr`: reference,
refined and three variants) regenerated by the director with `opt/fresco/bin/fresco`: byte-identical to the committed
`reference.txt`. The numbers the new tests print equal the agent's report: ours against FRESCO 4.3e-3, 1.8e-3, 5.2e-3
at NK = 1600 (1.1e-3, 8.7e-4, 2.6e-3 at NK = 3200); narrow-bin limit 6.3e-4 with error ratio 4.00; post = prior
1.1e-3. Not re-run by the director: the injected-error runs, the ⁶Li and remnant FRESCO cases, `run_tests.sh`, the
gate and the fits (nothing they test is touched).

**Reviewed.** The headers and `Continuum.cpp`, `Interaction.cpp`, `TransferKernel.cpp` read: the normalisation
factors of (14.2) and the bin's norm (14.3) re-derived; the coordinates of the remnant terms (r_sA = r + ρ,
r_aA = r + βρ, r_sF = αr + ρ) checked against the kernel's (R, ρ); no mutable state (two `static const` tables of
kinds, one of powers of i). The one change outside `physics/` and `docs/`: `tests/consistency/repo_consistency_test.py`
lets a test of `physics/tests/` include a helper header of that folder by bare name.

**Corrected by the director in the documents.** "Within 0.3 % of FRESCO" (owner note, roadmap, `THM.md`) held only
with the finer k rule (NK = 3200); at the benchmark's NK = 1600 the worst is 0.5 %. Written as 0.5 % (asked: 1 %).

**What it says about the engine** (`transfer.md` §16; the agent's table, reproduced by a test, not re-derived by the
director): inside the channel radius the continuum amplitude is the level's vertex times the level's energy
dependence, exactly at the resonance energy and to 0.1–1.9 % one width away; the nuclear part outside the channel
radius is 9 % of the interior at a = 5 fm and 0.5 % at 7 fm for p + ¹⁵N, l = 2. So the surface/full form factor of
the THM fit is the right limit when the channel radius lies beyond the nuclear potential.

**Not reached, each with its cause named by the agent** (FRESCO's own numerics; data committed, not compared):
zero range to a bin (FRESCO's zero range disagrees with its own finite range by ×7); ⁶Li bins (FRESCO loses accuracy
at high L; at RSP = 12 fm 0.1–1.5 %); prior form and remnants against FRESCO (FRESCO does not converge in its
non-local range; ours: post = prior to 1.1 × 10⁻³). FRESCO's own step refinement is 0.3–0.5 % at two minima, above
the 0.2 % asked. Not done: the ³He and ¹⁴N Trojan horses.

**For the PI.** Two questions of the agent: whether bins cut at a fixed radius are the right benchmark object (the
far tail of a bin is physics: 88 against 58 mb/sr at 0° for a cut at 40 against 60 fm, in both codes), and whether
the ⁶Li bins deserve another code than FRESCO. Neither blocks the wiring.

**Open.** The wiring into the THM observable (`transfer.md` §12: the global optical sets move from `ThmOptical` then,
so that there is one copy of the numbers); ³He and ¹⁴N.

## 2026-11-12 (agent, branch `feature/upgrade-7-transfer-2`) — physics/transfer, second part: the x + A continuum

**Done** (leaf library, nothing of the engine touched; `docs/upgrade/transfer.md` §11, §12, §14–§18):
- Scattering states of x + A (`Continuum.h`: amplitude 1, δ(k − k'), δ(E − E') normalisations) and FRESCO's bins
  (ISC 2, 4, 12, 14; FRESCO's k rule or Simpson's; FRESCO's listing).
- Overlap kinds as one table (`Overlap.h`: bound state, bin, scattering state → one `RadialOverlap` for the kernel);
  the interaction as a sum of terms (`Interaction.h`: post, prior, core-core and optical remnants); the finite-range
  kernel in its own file (`TransferKernel.cpp`), complex G.
- d²σ/dΩ dE to the continuum and the narrow-bin limit; the interior part against the level function (§16 table); post
  = prior with remnants.
- FRESCO data: `dn_bin_{nonresonant,narrow,wide}_{fr,zr}`, `li6d_bin_wide_fr`, `dn_prior`, `dn_prior_remnant`,
  `dn_post_remnant` (generator updated; `rsp` and NK settings).

**Verified** (build `work/b-dwba2`): new tests `physics_transfer_continuum`, `_continuum_transfer`,
`_continuum_fresco` pass; `physics_transfer_fresco` gives part 1's numbers unchanged. Injected errors, each failed its
test and was removed: bin scale √(2.1/πN) (continuum: narrow-bin factor and norm), sign of −U_aA (continuum_transfer:
post ≠ prior by ×180), cos δ for the ISC = 4 weight (continuum_fresco: listing and cross sections).

**Open**: ⁶Li bins, zero range to bins, prior/remnant against FRESCO (FRESCO's numerics, §15, §17); ³He and ¹⁴N
Trojan horses; global optical sets; the wiring (§12).

## 2026-11-12 (agent, Phase 6 branch 1, part 1) — the state of the engine surveyed (`feature/upgrade-6-evaluation`)

**Done.** `docs/upgrade/evaluation.md`: the state table (every `thread_local`, mutable static, namespace-scope
variable, `mutable` member and member cache of `src/`, `include/`, `api/`, `physics/`, `experiment/`, `stopping/`,
`convert/`, `io/`, with what it becomes), the 44 environment switches with a proposal each, the cost of an evaluation
today (time and memory at 1 and 3 threads, the shares of a fit), the target, five milestones, decisions for the
director.  Measured: pool copying is 0.02 % of a fit's time, its cost is memory (a natural_neon fit holds 812 MB
against 343 MB for the calculation: about five copies of 95 MB, `max(4, threads)` plus the snapshot pair).  The
A-B-A test on the base engine found no leak.
**Open.** Parts 2 (no hidden state) and 3 of the brief follow on this branch.

## 2026-11-12 (director, 18:27) — merged: Trojan-horse data on the shared mesh, their gradient through the energy tasks, the shift column of a profiled THM data set

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (the branch held the data head, so the merge's tree is
the branch's plus the director's entries; nothing pushed): `feature/upgrade-5-thm-mesh` (four entries of the agent
below say what was done, by file).

**What it is, in four lines.**  σ_HOES is a kind of the observables' table (`SiteObservable::kHoes`: a scalar without
an angle), so a THM point is folded like any other datum: its supports are planned on the mesh of its compound
nucleus, a site is computed once per mesh energy and experiment settings and shared by every reader, and where direct
data read the same energy the level matrix is solved once for both.  The gradient of THM data goes through the energy
tasks (the model cotangent of a profiled scale by the envelope theorem, the HOES adjoint once per site).  No data kind
has an energy grid of its own any more, except a point under `useAdaptiveGrid` off.

**A defect of the merged engine, found by the director and fixed here.**  The analytic gradient's column of a fitted
energy shift on a THM data set with a profiled norm could hold only the prior's derivative (the local difference
compared the models at the start norm 1 with data that need n* of about 3e4).  A project with THM data takes the
analytic gradient in every MIGRAD fit ("THM segments present: enabling the analytic gradient"), so plain MIGRAD was
affected as `--use-gradient` was.  The director's runs on 18O(p,α) THM with a prior on the shift
(`work/director/thmshift/`, `work/director/thmshift5/`), merged engine against this branch:
- the gradient check with the level fixed: analytic −1973.6 against −1459.3 by differences; now −1241.977222 against
  −1241.977216;
- the fit with the level fixed and the shift free: both end at −6.698 keV and converge, in 31 calls before and 16 now;
- the fit with the level energy and the shift free: before, MIGRAD ended **invalid** (EDM 0.16 against a goal of
  2e-4, 47 calls; the warning is printed) with an error of the level energy of 1.8 keV; now it converges (EDM 1e-7,
  36 calls), error 0.73 keV; the end points agree (shift −5.650 keV, E 8.61657 MeV).
`tests/thm_fit` checks the column with and without a prior (run by the director: fails on the merged engine, passes
here).

**Verified by the director** (build `work/b-thm5` of the branch's last commit, reference `tab-trial36` = the merged
engine before it; the new reference of the merged engine is `tab-thm5`).
Every project on three bases against the merged engine (`tab-trial36`): only T2, T2_true, 18O_p_a_thm, 6Li_d, 7Li_p_a
and 17O differ. Against the director's refined references (quadrature ×8 and ×16, consistent): 0 failures by the
accepted per-file rule. Eight fits: the five without THM data identical, the three T2 fits at χ² 297.601 → 297.568
with the largest parameter move 0.006 σ and equal statuses. Unit tests 213/213, slow tests 10/10, converter
equivalence 184/184, the converted-tree gate without a difference, 1 against 3 threads identical. `run_tests.sh`: 159
of 160 on the agent's head; the one failure was the recorded χ² of `tests/convert_params` (T2 at the old fit's
parameters), which the agent had not re-run: merged engine 293.705, this branch 293.675, both refined eightfold
293.647; re-pinned by the director (293.703 → 293.675), the check passes. Gradient checks: T2 worst 2.6e-6, 18O
4.9e-7; `tests/thm_fit/check.sh` fails on the merged engine and passes here. Notes and runs:
`work/director/thm_review_notes.md`, `work/director/thm5/`, `work/director/thmshift5/`.

**Reviewed, with the director's own runs** (`work/director/thm_review_notes.md`, `work/director/thm5/`).
- *Numbers.*  Whole tree × 3 bases: only the six projects with a THM fold move (17O, 18O(p,α) THM, 6Li+d, 7Li(p,α),
  T2 and its truth file).  Against the merged engine refined 8-fold and 16-fold (the director's runs, tags `dir-fine8`,
  `dir-fine16`) every changed file is closer to the refined numbers than before, in the largest and in the rms
  distance (the accepted per-file rule, 0 failures): T2 1.2e-3 → 4.9e-4, T2's THM table 9.5e-4 → 1.5e-4, 18O 4.2e-4 →
  2.7e-4, 7Li(p,α) 7.1e-5 → 3.9e-5, 6Li+d 1.9e-5 → 1.4e-5, 17O 1.06e-4 → 1.04e-4.
- *T2's direct data move too*: the thin and thick target segments 1 and 4 share the mesh of 16O with the THM
  supports, which adds nodes to their folds (χ² 697.326 → 697.517, refined 697.74; 286.444 → 286.473, refined 286.59;
  THM segments 92.0754 → 92.148, refined 92.163; 121.489 → 121.586, refined 121.606; the profiled norm 1.96991e-3 →
  1.96953e-3).  The T2 expectation was re-pinned by the agent for these files; its `normalizations.out` now holds the
  follower's row (segment 3) that the comparison used to add.
- *Gradient* (`AZURE_GRAD_CHECK`, the director's runs): T2 2.6e-6 (19 components, the same as before), 18O with the
  level free 4.9e-7.  `AZURE_NODE_CHECK` on T2: 0 differences.  Gradient time at one thread, 18O: 291 → 23 ms, one χ²
  31 → 3.4 ms; T2 unchanged within 5 %.
- *One stored number changed*: `tests/thm_options`, `kinematics=lambda32`, 5604.98 → 5561.64 (the agent's refined
  runs: 5510 and 5507): 1/k_i at the entrance threshold inside a fold, the χ² rests on the quadrature there.
- No new static, thread-local or mutable state; +636 −320 lines in `src` and `include`.

**For the PI and the owner** (behaviour a user sees).  (1) Folded THM numbers move by the quadrature error they had
(table above); direct data that share a compound nucleus with THM data can move in the 4th digit of a folded value
too.  (2) A MIGRAD fit with a free energy shift on a THM data set with a profiled norm could use a wrong gradient
component until this merge (seen as an invalid minimum and a wrong error).  (3) THM fits are faster (18O: the gradient
13 times at one thread).

**Decision asked of the PI** (not taken here): the `<config>` key `useAdaptiveGrid`.  Off, a folded point keeps a
uniform grid of its own, the last path beside the mesh (`ShiftedFold.h` holds the agent's proposal); the same key
selects the reaction-rate grid (`ReactionRate.cpp`).  Either a mesh setting for a uniform lattice (the key stays), or
the key goes for target integration with a note in the converter.  Until decided the path stays.

**Open** (carried).  Not deleted by this branch, for one later small branch after a check on the tree: the leaves'
angle copying and `SeparateSlotsByAngle` (no reader of the angle was found), the central-difference fold weights and
the per-(point, slot) adjoint of `ComputeTHMRows` (it stays for rows: band, least squares, the API, and for
components, a coherent background, RMC).  Still by differences: THM norm columns, shift columns, the coherent
background.  "The shift on thin kernels" (the entry of 17:01).  Phase 7 part 2 (`feature/upgrade-7-transfer-2`, the
transfer to the continuum against FRESCO) is in work in `physics/` only.

## 2026-11-12 (agent, Phase 5 branch 10, after the merge of follow-ups 2) — the shift column of a profiled THM data set (`feature/upgrade-5-thm-mesh`)

**Merged** `feature/upgrade-4-data` (f952c444): HANDOFF kept both sides; THM supports pass through `QueueMeshPoint`
(label, kernel width for the node budget), THM shift columns take `InsideShiftStep`.
**Cause** (the director's finding): the local central difference of a shift column (`SegmentLocalChi2`) compared the
shifted models of a profiled THM data set with its data at the segment's start norm (`State().norm`, 1), not at its
profile n* (about 3e4 here), so the models were negligible against the data and the column was the prior's derivative
alone -- in the base engine and on this branch alike.
**Fix.**  `SegmentLocalChi2` moved to `ShiftedFold.cpp` (out of the size-ratcheted `AZURECalc.cpp`) with the norm as an
argument; `HeldNorms` gives each moved segment its own norm, or for a profiled one its profile n* at the present shift
(`Segment::ProfileNormChiSquared`, the objective's function), held fixed across the difference (envelope theorem).  An
experiment's segments stay outside the local path (whole chi^2 differences, as before).  `ThmModelCotangents` takes the
profiled scale of a data set without an experiment from `ProfileNormChiSquared` too (the second copy of S_md/S_mm is
gone).
**Checks.**  `AZURE_GRAD_CHECK`, folded 18O with the 8.6126 level free and a fitted shift: with a prior, shift column
-1494.8 against central differences 4269.1 before (rel 0.9998), 4269.13 = 4269.13 after (rel 1.3e-9); without a prior
-1.45 against 4388.2 before, rel 1.3e-11 after.  `tests/thm_fit` checks both (shown failing on the merge commit's
build).  The least-squares path: a project with THM data has no analytic Jacobian, so `--use-lm` falls back to MIGRAD
with whole-chi^2 differences (as the T2 golden notes say), which profile: its shift column was right before and after.
Fit of the director's case (level fixed at 8.6126, shift fit with prior 0.002 +- 0.003, 3 threads): `--use-gradient`,
`--use-lm`, plain MIGRAD all end at x = -6.6982 keV +- 0.49 keV, chi^2 1507.02 (to 1e-10 in x); the scan of chi^2 at
fixed shifts plus the prior (-8.0 ... -6.0 keV: 1522.43, 1518.09, 1515.80, 1515.59, 1517.44) has its parabola's minimum
at -6.70 keV.
**Ran**: build 0 warnings; consistency; quick set x 3 identical to milestone 1 (the merge's build and the fix's);
`tests/mesh_node_budget`, `thin_layer_limit`, `ay_particle_entrance`, `tests/upgrade/check.sh`, `tests/thm_fit` pass.

## 2026-11-12 (agent, Phase 5 branch 10, milestones 3 and 4) — the THM gradient through the energy tasks; the documents (`feature/upgrade-5-thm-mesh`)

**Gradient.**  The scalar pass (`EnergyTaskGradient`, CLI `AZURECalc::Gradient`) takes THM data when
`ThmThroughEnergyTasks` holds (no components, no coherent-background parameter, not RMC, no AZURE2_THM_FD): a THM
point's model cotangent is 2 r s f/e from its profile (`ThmModelCotangents`, `include/ThmCotangents.h`: at a profiled
norm and background d chi^2/d s = 0, the envelope theorem), its fold's weights take it to its sites as for any datum,
and `HoesTask` (`GradEnergyTasks.cpp`) sums the cotangents of all readers of a kHoes site and takes the HOES adjoint
once per site on one level-space fill per task.  T-matrix groups and fills skip kHoes units (`ReadsTMatrix`).
`ComputeTHMRows` stays for rows (band, LM, API) and for what the predicate refuses; its fold weights are still central
differences, one adjoint per (point, slot).  Still finite-difference: THM norm columns (a fitted or prior norm, through
the THM-only chi^2), shift columns, the coherent background (the old path), AZURE2_THM_FD=1.
**Checks** (`AZURE_GRAD_CHECK=exit`, `--use-gradient`, base `b-obs5` against this branch, worst relative disagreement):
T2 2.62e-6 = 2.62e-6 (19 components, 0 above 1e-2 in both); folded 18O variants with the 8.6026 level free (thm_fit's
edit, fold kept) plain 4.9e-7 =, experiment 4.8e-7 =, lineshape 2.9e-7 =, ps window 3.7e-6 =, DW vertex 3.6e-8 =.
Gradient time (`[grad-timing]`, 3 threads) 18O 88 -> 8.7 ms, T2 372 -> 316 ms.  1 against 3 threads: T2 and 18O-ps
gradients identical to the last digit.  A fitted THM shift (18O plain, offset fit prior 0.002 0.003): the shift column
disagrees with central differences in the base build and here alike (-1496 vs 4229; rel 0.9998) -- a defect of the
shift column of a profiled THM segment, not of this branch (listed in the report, not fixed).  New test:
`tests/thm_fit` (folded experiment, ps window: rel < 1e-5); shown failing with the profiled scale left out of the
cotangent (rel 0.9987, 0.989).
**Ran**: build 0 warnings; consistency; ctest subset 58/58; quick set x 3 identical to milestone 1; `tests/thm_fit`.
**Documents** (milestone 4): `energy_grid.md` section 8 (row 10 and its status paragraph), `THM.md` section 5,
`FOR_THE_OWNER.md` (a paragraph), the recipe "An observable of a reaction at a mesh energy" in `CLAUDE.md` (a kind that
reads the level matrix's elements).  Not deleted (see the report): the leaves' angle copying and
`SeparateSlotsByAngle` (they act on every scalar kind's leaves; no reader of the angle found, deletion not verified on
the tree), the uniform grid's own path (a user-visible key; proposal in `ShiftedFold.h`).

## 2026-11-12 (agent, Phase 5 branch 10, milestone 2) — sigma_HOES a kind of the table; one solve per node for THM and direct readers (`feature/upgrade-5-thm-mesh`)

**Done.**  `SiteObservable::kHoes` (`SiteObservable.h`, its row in `kSiteObservables`): the half-off-shell cross
section, a scalar without an angle (the record's sigma, read as it is; no row, no contraction).  `Site::AngleFree()` is
now "the kind's readers bring a row" (`ReadThroughRow`, the row's `contract`), so kValue and kHoes are both scalars.
The HOES formula is a function of a filled level matrix (`HoesCrossSection(const AMatrixFunc &, Site *, Results &)`,
`THMMatrixFunc::CalculateTHMCrossSection` calls it).  The plan has no THM route: an energy task is (nucleus, route of
the configuration, E_x bits), whatever its kinds; `SolveAtOneEnergy` fills with the first unit that reads a T-matrix
and gives each kHoes unit `HoesCrossSection` on that fill when no J-group took the channel route
(`AMatrixFunc::LevelSpaceThroughout`), else on a level-space fill of the task's THM units; the T-matrix groups skip
the THM units.  `CheckEnergyTasks` (`AZURE_PLAN_CHECK`) compares every unit's record from the task's shared solve with
its own solve, of its own kind, and prints the THM tasks and how many share one fill.  `NodeKey`'s THM flag stays: it
keeps the rows a THM site holds beside the channel rows (j_l, rho j_l', the Coulomb term, the window's table) on a node
of their own; the solve is the task's.  The sharing key compares a data set's THM weight table by value
(`SameThmWeight`; check 37 reads `SameX(a.f, b.f)` as a comparison); the experiment's attachments by identity (one
object per experiment).
**Numbers.**  Quick set x 3 bases against milestone 1 (`thm-m1a`): 360 of 360 files identical, run.log included.
**Solves.**  T2: 4,227 energy tasks -> 3,055; 1,172 of its 1,192 THM nodes are solved with the direct readers at
the node, every one by one fill (all J-groups of 16O take the level route at those nodes).  `18O_p_a_thm` (THM only):
656 tasks, unchanged (on milestone 1 its 14,688 units of the base became 656 slots of one block).
**Tests.**  `upgrade_plan_runtime` asserts T2's shared fills; injected once a THM task reading a level matrix it did not
fill: 18O "656 differ", abort.  `upgrade_node_check`: T2's s11/s12 weighted by two tables with equal rows keep 9
blocks, with other rows 10; with the tables compared by pointer (injected once) 10 and 10.  `upgrade_plan`'s route
check: a task is on the configuration's route.
**Ran**: build 0 warnings; consistency; ctest subset (58; `upgrade_plan` failed on the old route rule, fixed, then the
four plan/node/consistency tests again: pass); AZURE_PLAN_CHECK and AZURE_NODE_CHECK on T2 and 18O_p_a_thm: 0 differ.

## 2026-11-12 (agent, Phase 5 branch 10, milestone 1) — THM data on the shared mesh (`feature/upgrade-5-thm-mesh`)

**Design** (the six questions, answered in the headers that own them, commit 485badd7): `SiteObservable.h` (what a THM
site holds: sigma_HOES at E_n is a function of the node, the reaction and the THM experiment only -- nothing reads the
data point's angle or row, the spectator angles are the experiment's), `TermSettings.h` (the key), `MeshSlots.h`
(support, floor, nodes), `THMMatrixFunc.h` (the kernel is the table's; the level matrix shared per node),
`ShiftedFold.h`
(shifts; the uniform grid's proposal), `AZUREGrad.h` (the gradient).
**Done.**  A THM point is a folded point of the shared mesh of its compound nucleus (`OnSharedMesh` asks no THM
question): its support is the kernel's (Gaussian or beam profile, 5 sigma) clipped at the data kind's floor,
`Segment::IntegrationFloor` (E_x = 1 keV of the compound system for half-off-shell data, so E_cm of a node may be
negative), a hard edge of the planned support there; a THM shift is the shift kernel's case (`ShiftedFold.h`), not
floored (`Segment::ShiftIsFloored`, the one question the shift step asks).  Deleted: the own-grid re-placement of
`ReplaceSubPointGrids` (its `replace` path, the reshift list it returned to the hand-over, the `shifted` branch of its
site layout) and `AdaptiveIntegrationGrid::GridConfig::resonanceOffset` (their only user was THM); a folded point off
the
mesh with `useAdaptiveGrid` on is now a `logic_error`.  `Segment::MoveOwnGridSlots` and `Unit::CalcNodeRows` are left to
the uniform grid (`useAdaptiveGrid` off) alone; `CalcNodeRows` never served THM (a scalar kind has no row set).  Lines:
src +56 -121, include +21 -23, tests +5 -2.
**Moves** (`work/thm/quick.sh`, tags `thm-m1a` against `tab-obs5`, three bases): T1, T3, 13N, natural_neon,
`thm_coherent_background/toy` (unfolded THM) identical; the THM projects 17O, `18O_p_a_thm`, 6Li_d, 7Li_p_a, T2, T2_true
moved.  T2's direct segments move too (thin target 697.326 -> 697.517, thick 286.444 -> 286.473): the THM supports join
the mesh of 16O and the planner's nodes change where they overlap.  Classified (`tools/classify_quadrature.py`, base
`tab-obs5`, the base engine refined eightfold `thm-fine8`, made here): Brune 16 quadrature, no-brune 12, Park 15, 0
failure; every file closer to the refined run than before (largest relative distance new / before: T2 `AZUREOut` 4.9e-4
/ 1.2e-3, `thm_experiments.out` 1.5e-4 / 9.5e-4; 18O 2.7e-4 / 4.2e-4; 17O 1.04e-4 / 1.06e-4; 6Li_d 1.4e-5 / 1.9e-5;
7Li_p_a 3.9e-5 / 7.1e-5).  One pin regenerated: `tests/thm_options` kinematics=lambda32, chi^2 5604.98 -> 5561.64
(refined 8x 5510.14, 16x 5506.51; the THM curve's largest distance to the 8x run 0.158 -> 0.082): its 1/k_i reaches
1/sqrt(1e-6 MeV) at the threshold inside the folds.
**Ran**: build 0 warnings; consistency test; the quick set above; the check.sh of `18O_p_a_thm`, `thm_band`,
`thm_beam_profile`, `thm_coherent_background`, `thm_coulomb_consistency`, `thm_energy_shift`, `thm_experiment`,
`thm_lineshape`, `thm_narrow_fold`, `thm_options` (with the new pin), `thm_park`, `thm_rwa_grid`,
`thm_spectator_window`, `thm_threshold_grid` pass.
## 2026-11-12 (director, 17:01) — merged: the second batch of follow-ups (thin layers with a shift, unknown options, the gradient's gates, a node budget, the runner)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (the source of `rehearsal/trial36` = the observables
head with the batch; nothing pushed): `feature/upgrade-5-followups-2` (the agent's two entries below say what was
done, by file).

**What a user sees.**
- A layer thinner than 1 eV with an energy shift on its data set is calculated (the thin limit at the shifted
  energy); it was a failed point.
- A gradient or least-squares fit with a free energy shift on a thin layer no longer ends the process: the shift
  columns keep their difference step inside the planned nodes (`InsideShiftStep`, `src/ShiftedFold.cpp`).
- An unknown command-line option is an error (exit 2 with the usage line); `--use-nlopt` in a build without NLopt is
  one of them.  `--no-readline` exists in every build.
- An analysing power on a β-delayed or photon entrance is refused when the file is read (it was calculated as 0, and
  its gradient was wrong by 1e-3).
- A mesh that would need more than 1e6 nodes (`kMeshNodeBudget`, `include/MeshPlanner.h`) is refused with a message
  that names the data set, the point, the kernel's width, the planned width and what to do; before, such a project
  ran out of memory.
- `tools/runall.sh` refuses a tag that exists (it mixed the results of two runs silently) and ends with exit 1 when a
  project failed.

**Verified by the director** (build `work/b-trial36` of the rehearsal, reference `tab-obs5`).
- Every project on the three bases (835 / 808 / 835 output files, fresh tag `tab-trial36`): identical to the merged
  engine (`tab-obs5`) to the byte.  Eight fits against it: identical end points, χ² and times (T2 297.6010).
- Unit tests 213/213 (two new: `shift_column_inside`, `tools_runall`), slow tests 10/10, `run_tests.sh` 160/160 (two
  new: `ay_particle_entrance`, `mesh_node_budget`; the T2 golden comparison passes with the re-pin of the observables
  merge), converter equivalence 184/184, the converted-tree gate without a difference, 1 against 3 threads identical.
- The director's own cases (`work/director/batch2_checks.sh`): a 1e14 layer without a shift, with a fixed shift of 0
  and of 20 keV, and 1e16 with 20 keV all run; a free shift without a prior under `--use-lm` at 3e14, 1e15, 1e16 and
  under `--use-gradient` at 3e14, 1e15 ends with exit 0, where the same 3e14 fit on the build before the batch ends in
  `std::terminate` ("the upper edge of the moved support ... is outside the point's planned nodes"); on the test's
  flat χ² the gradient fit at 1e16 walks the shift away and is stopped by the node budget with the message (exit 1);
  an unknown option and `--use-nlopt` exit 2 with the usage line; the runner refuses an existing tag (exit 2, nothing
  touched).

**Reviewed; conflicts with the observables branch resolved by the director.**  `AMatrixFunc::PointTBar`: the batch's
assertions on the observables branch's row (the coefficient flag and the Legendre row are the reader's).
`GradTargetEffectAdjoint`: the batch's `logic_error` for a multipole part of a ratio, on the one `ratio` flag
(`Unit::FoldsRatio`).  Fixed beside it: the `AZURE_SHIFT_GRAD_CHECK` self-check divided by the nominal step, not by
the step it took.  Notes and runs: `work/director/batch2_review_notes.md`, `work/director/batch2_checks.out`.

**Decided.**  Not done in this batch and left for one later small branch, "the shift on thin kernels":
- the planned support of a shifted point is filled at the kernel's own step, so the cost grows as prior width over
  kernel width (the agent's two proposals are in its entry below: a floor on the reach with a step that grows, or a
  capped reach with re-plans);
- a thin-limit point (kernel below 1 eV) with a varied shift: χ²(ν) is a staircase, MIGRAD's error of the shift is
  too small there; the remedy is the unfolded path with the thin-limit factor;
- 6.2e-6 at one energy of `thin_layer_limit`'s target-averaged σ under a −20 keV shift, unexplained;
- the node-budget error inside a fit is thrown in Minuit2's parallel gradient and ends the process by
  `std::terminate` (the message is printed, exit 134); a clean end needs an abort flag through the minimizer (Minuit2
  is not patched).

**Open** (carried).  THM on the mesh (`feature/upgrade-5-thm-mesh`).  Small: `MeshSupport` carries a label string for
the budget's message (it could be made from the data set's index when the message is written); a window beyond the
kinematic maximum without a shift could be refused at set-up; the two near-copy loops of `SetUpTargetEffects`.

## 2026-11-12 (agent, follow-ups 2, the director's fixes) — shift columns, A_y entrance, node budget (`feature/upgrade-5-followups-2`)

**Done** (one commit each).
- *Shift columns inside the planned nodes.*  The local shift columns (gradient, residual Jacobian) moved the points by
  +-h without a re-plan; with |Delta - c| + h > R the refold threw and the process ended (3e14: R 0.55 eV < h 0.96 eV;
  1e15, 1e16: a shift within h of its window's edge; a 0.1 eV beam spread too).  `InsideShiftStep` (ShiftedFold.cpp)
  halves the step until both sides stay inside, else the column takes the path that plans again.  Test
  `shift_column_inside`.
- *A_y needs a particle entrance*: refused in `experiment/src/Build.cpp`; the PointTBar gate is an assertion.  Test
  `tests/ay_particle_entrance`.
- *Node budget*: PlanMesh refuses a plan beyond `kMeshNodeBudget` = 1e6 nodes (50 x natural_neon's 21,070) before
  making it, naming the data set, point, kernel width, planned width and step.  Set-up: ERROR, exit 1; inside a fit
  (Minuit's parallel gradient): std::terminate with the message, exit 134.  Test `tests/mesh_node_budget`.

**Open, design.**  The planned support of a shifted point is filled at the kernel's own step, so the cost grows as prior
width / kernel width (1e14 layer, kernel 6e-8 MeV: a 10 eV prior 1.3 s, 0.1 keV 10.6 s and 2.4 GB, 3 keV ~3e7 nodes,
refused now).
- A floor on the reach (a few eV) and a planned step that grows with the planned width (the kernel's number of points
  across it): bounded cost; moves the nodes of every shifted point whose planned width exceeds its kernel's, i.e. every
  varied shift with a prior or without (their outputs at the quadrature level).
- A reach capped at a multiple of the kernel width, with re-plans as the shift walks: nodes as without a shift; moves
the
  same projects' fit paths (more re-plans, chi^2 steps by the quadrature difference at each).
- Thin-limit points (kernel below 1 eV) with a varied shift: chi^2(nu) is a staircase (MIGRAD's error 10x too small at
  1e14); routing them through the unfolded path with the thin-limit factor is not done (cost in the report).
- At -20 keV with a 10 eV prior, thin_limit's target-averaged sigma at 0.6 MeV is 6.2e-6 off (0.7, 0.8 MeV and the yield
  3e-7, as expected): the size of linear interpolation over a keV interval, not the planned 2e-9 MeV step; not resolved.

## 2026-11-12 (agent, follow-ups 2) — thin shifted layers, unknown options, the gradient's gates, runall tags (`feature/upgrade-5-followups-2`)

**Done** (one commit each).
- *A layer below 1 eV on a data set with an energy shift.*  The point was NaN ("the lower edge of the moved support ...
  is outside the point's planned nodes"): its planned support is still below `MeshZeroWidth()` (fixed shift: no reach;
  varied without a prior: reach 3 support widths), the planner gives one node (its upper edge), and `ShiftedFoldNodes`
  compared the moved support's edges with it to 1e-9 MeV.  For one node it now checks that the moved support lies
  within `MeshZeroWidth()` below the node and the thin limit folds it (`src/ShiftedFold.cpp`; header paragraph on
  single nodes).  `tests/thin_layer_limit/check.sh` gained a shift scan at 1e14 and the 1e14–1e16 scan with +20 keV
  fixed and varied, against the bare data set at the shifted energies; fails on the base.
- *An unknown command-line option is an error* (exit 2, the usage line), also in `AZURE_OPTIONS_FILE`; `--no-readline`
  is accepted in every build.  `tests/cli_option_values`.
- *The gradient's gates.*  In `AMatrixFunc::PointTBar` the gates no input reaches are assertions (`std::logic_error`
  naming the broken invariant; the reason at each), one deleted; `GradTargetEffectAdjoint`'s E1/E2 gates likewise.
  `park_chi2_fallback` variant D: `--use-rmc` (the fallback against central differences).
- *`tools/runall.sh` / `run_one.sh`.*  A file's results are replaced as a whole; a failed move is an error (run_one
  exit 1, runall "FAILED", exit 1); an existing tag is refused unless `REPLACE_TAG=1`. 
`tests/consistency/runall_test.py`.

**Verified** (scratch `work/small2/`): builds 0 warnings; consistency test; `cli_option_values`, `cli_help`,
`cli_route_exit`, `thin_layer_limit`, `park_chi2_fallback`, `tools_runall`; T1–T3, 13N, 13N_capture_ay,
thick_target_resonance, layer_yield/models, shift_sqrt_subpoints/models, shift_subpoint_angle(+/models),
thin_layer_limit/models in three bases byte-identical to `tab-small1`; each new test shown failing on the base or with
an injection; `scripts/check.sh --quick` (JOBS=1) at the end.

**Not run**: the whole tree × 3, `run_tests.sh`, the converted-tree gate, the fits of the tree.

**Open** (details in the report).
- A varied shift on a thin layer: the analytic gradient's shift column (`--use-gradient`, `--use-lm`, `--use-gsl-lm`)
  ends the process (`std::logic_error` from `ShiftedFoldNodes` inside the parallel refold) when the reach is below the
  difference step (~1 eV): seen at 3e14, 1e15 and 1e16 atoms/cm² on the base too.  MIGRAD's numerical gradient works.
- A prior of keV on such a layer: the planned support keeps the layer grid's smooth step (millions of nodes, 28 GB).
- A varied shift without a prior on a single node: chi^2(nu) is a staircase of sigma' x 2R; MIGRAD's error of the
  shift comes out 10x too small at 1e14.
- An analyzing-power data set on a beta-delayed reaction is accepted and calculated as 0; its analytic gradient is
  wrong by ~1e-3 (PointTBar treats the point as the beta-delayed spectrum).
- Gates of PointTBar not settled (kept as refusals): beta-delayed with an E1/E2 part, the A_y sub-point gate
  (`ayCrossSection`), the capture-A_y term gates, a vanishing denominator, an empty amplitude matrix, an A_y exit that
  is neither particle nor photon.

## 2026-11-12 (director, 16:14) — merged: observables per mesh energy and reaction, the angle in the reader's row, the lab-frame fold per node (D5)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (the branch started from the follow-ups' head, so the
merge's tree is the branch's plus the director's entries; nothing pushed): `feature/upgrade-5-observables` (five
entries of the agent below, one per milestone, say what a site holds, what the key is and what went, by file).

**What it is, in four lines.**  A site of the mesh holds the observables of a reaction at one energy without an angle
(the integrated cross section and Legendre coefficients, the elastic amplitudes, the analysing power's coefficients:
one table of kinds, `include/SiteObservable.h`); a data point brings its angle in a row (`include/RowSet.h`), and
points of one block with the same angular identity share the row and the contracted value.  Two data sets of one
reaction on one mesh share their sites (the key has no segment, no angle, no Legendre order).  A point given at a lab
angle is folded in the lab: θ_c.m.(E_n; θ_lab) and dΩ_lab/dΩ_c.m.(E_n; θ_lab) at every node (decision D5).

**Verified by the director** (build `work/b-obs5` of the branch's last commit, references `tab-small1c` = the merged
engine before it and `tab-obs4d` = milestone 4; the new reference of the merged engine is `tab-obs5`).
- Every project on the three bases (835 / 808 / 835 output files, fresh tag `tab-obs5`): identical to milestone 4's run
  (`tab-obs4d`) to the byte, the run logs included; against the merged engine before the branch the known set only
  (T2's `covariance.out` 1.5e-5; six observable files, the largest the `ec_e1m1_phase` case below).
- Eight fits against the merged engine before the branch: T2 ends at χ² 297.6010 (297.5990 before D5) with
  `--use-lm`, `--use-gradient` and plain MIGRAD, every parameter within 0.004 of its error, in 64 / 44 / 26 s
  (88 / 65 / 38 s before); the two 13N fits at the same χ² and parameters to 8 digits; `layer_yield` and the two
  stopping fits identical.
- Unit tests 211/211, slow tests 10/10, converter equivalence 184/184, the converted-tree gate without a difference
  (633 / 606 / 633 files), 1 against 3 threads identical (12 runs and the Markdown projects).
- `run_tests.sh`: 157 of 158 at the branch's last commit; the one failure was `tests/upgrade`'s comparison of T2 with
  its golden files (1e-5), which D5 moves by 2.56e-5 and the branch had not re-pinned.  Re-pinned by the director
  (58c8e96c: segments 2 and 3, the two data sets at a lab angle; χ² 99.168 → 99.1682, 414.736 → 414.728 at the start
  values, 40.3817 → 40.3825 at the truth; nothing else moves); `tests/upgrade/check.sh` passes with it.
- Speed (one thread, side by side): `beam_profile_kernel` 5 min 27 s → 16 s and 913 → 134 MB, outputs identical.

**Reviewed, with the director's own runs** (`work/director/obs_review_notes.md`, `work/director/grad_obs5/`,
`work/director/speed_obs5/`, `work/director/labfold/`).
- *Milestones 1–3 (rounding only)*, whole tree × 3 bases at each: every number within the rounding class but one file,
  `ec_e1m1_phase/models/s_1` at 0° and 180°, where an interference of 1e-25 b/sr is the difference of terms of 1e-20
  (0.43 relative on a number 1e-5 of its neighbours; the same bytes at every later milestone).
- *Milestone 4 (D5), the one physics change.*  Only T2 and its truth file move: at most 2.6e-5 of a model value,
  1/1900 of the datum's error.  `tests/lab_frame_fold` (7Li(p,n) at 30° lab, a 300 keV window 0.4 MeV above
  threshold) run by the director with both builds: the new fold agrees with a brute-force integral of the unfolded
  engine to 2.2e-5 (the quadrature), the fold before D5 was 2.2e-2 off.  **The agent stopped on the brief's stop
  condition** (six of T2's 60 points moved by more than ten times the brief's estimate); accepted: the brief's estimate
  lacked the first-order term s²(ln σ)′(ln f)′, and with it the moves are explained (`energy_grid.md` §4.3, T2's
  table).
- *After milestone 4*: the ratio fold and its adjoint (A_y, a normalised coefficient) are one branch each, read from
  the table (`Unit::FoldsRatio`, `RatioPartsAt`, `RatioWeightRoute`); the director's question about a node beyond the
  kinematic maximum of a lab angle found an abort (an uncaught exception at set-up when the reach margin of a fitted
  shift held such a node), now a failed point only for a fold that reads the node (`tests/lab_frame_fold` case (d),
  run by the director: at offset 0 the folded value is the unshifted one bit for bit; the milestone-4 build fails the
  case).  Milestone 5 deleted the key's angle and Legendre order and the code that served sites with an angle.
- *Gradient* (`AZURE_GRAD_CHECK`, the director's runs on the last commit, largest relative difference to central
  differences): T2 2.6e-6, 13N 2.2e-7, the folded analysing power (`beam_spread_ay`, 202 parameters) 1.7e-7, capture
  A_y 1.2e-9, the lab-angle fold 3.6e-6, `layer_yield` 2.0e-9.  `AZURE_NODE_CHECK` on T2: 60 readers, 2 row sets, both
  at a lab angle, 0 differences.

**Corrected here.**  `FOR_THE_OWNER.md` said a lab angle beyond the kinematic maximum inside the window "is refused at
set-up": since the reach fix it is a failed point (NaN, χ² 1e100, a warning that names the cause and the remedy).

**For the PI and the owner** (behaviour a user sees).  (1) D5: a folded point given at a lab angle changes by the
window's width times the slopes of the cross section and of the angular factor; small for the tree's data (2.6e-5),
2 % for a wide window near a threshold.  (2) A window or a shift that reaches an energy where the lab angle does not
exist gives a failed point with the reason, not a refusal of the file.  (3) An analysing power given at a lab angle
is still refused by the format (its fold adjoint exists and is untested).  (4) Speed: see the verification paragraph
(`beam_profile_kernel`, the project with the most angles per energy).

**Open** (carried).  THM on the mesh (`feature/upgrade-5-thm-mesh`, started from this tree): the last data kind with
a grid of its own; its site must become angle-free first (the key has no angle).  The second batch of follow-ups
(`feature/upgrade-5-followups-2`): a thin layer with a shift, unknown options, the fallback's gates, a crash of
gradient fits with a free shift without prior on thin layers, a node budget.  Small: a window that reaches beyond the
kinematic maximum **without** a shift could be refused at set-up again (clearer than a failed point); a double-valued
lab angle (two c.m. solutions) takes one branch, as before; a ratio of two coefficient leaves needs a grammar field
(`experiment/` refuses it).

## 2026-11-12 (agent, Phase 5 branch 9, milestone 5) — the ratio fold from the table, a lab angle beyond the kinematic maximum, the deletions (`feature/upgrade-5-observables`)

**Decision on milestone 4** (the director): the D5 moves are accepted; the estimate of `energy_grid.md` 4.3 now has both
terms, s²(ln σ)′(ln f)′ + ½ s² f″/f, with T2's table.  `FOR_THE_OWNER.md` says it in plain words.
**Remark 3, done.**  A ratio observable is folded by one branch, (K·(w r))/(K·w), and differentiated by one branch:
`Unit::FoldsRatio` (a kind with `values = 2` in the table, or a normalised coefficient), `RatioPartsAt` (weight and
element of a slot), `RatioWeightRoute` (where the weight's cotangent goes).  Same products in the same order: no number
moves; the folded A_y variant's gradient lines are identical to milestone 4's.
**A lab angle beyond the kinematic maximum at a node** (the director's question).  With a fitted shift the block holds
the shift's reach margin; a margin node beyond the kinematic maximum made `Unit::RowSetOver` throw -- uncaught, the run
aborted (7Li(p,n) at 30° lab, offset prior 0.3 MeV).  Now such a slot holds NaN angle and J; a fold that reads one is a
failed point (NaN, the engine's warning: "the lab angle is beyond the kinematic maximum at an energy its fold reads (an
energy shift or the window reaches there); give the point at a c.m. angle, or narrow the window or the shift's
prior."), a fold that does not is unaffected.  `tests/lab_frame_fold` case (d) (offset 0: finite; -0.32 MeV: NaN with
the warning; the milestone-4 build aborts).
**Untested branch.**  The A_y branch of the fold adjoint at a lab angle (the weight's cotangent divided by J_n) is
unreached: format 3 refuses an analysing power at a lab angle (`E_NOT_BUILT`), so it is untested.
**Milestone 5, deleted** (nothing else served the sites with an angle): `SlotGroupKey::legendreOrder` and `angle` (a
differential site on the mesh has no angle, a scalar site has none, THM -- keyed by angle -- is off the mesh and never
keyed); `GenMatrixFunc::CalculateAmplitudeMatrix` with its `AZURE2_POL_DEBUG` hook, `CalculateCaptureAnalyzingPower`,
the
local-record path of a differential `kValue` site (now a logic_error), `Site::ContractRecord`, `Site::QLabBoost`.  Kept:
the points' own Legendre rows and Coulomb amplitude (unfolded points), the slots' entrance symmetry, the leaves' angles
and `SeparateSlotsByAngle` (THM leaves), `Polarization::AmplitudeMatrix` with an angle (the A_y adjoint).  Lines
(c8f53adf -> here): src +86 -232, include +28 -49, tests +36 -3, docs +40 -18; the deletion commit alone src +8 -152,
include +12 -44.  `GenMatrixFunc.cpp` 890 -> 756 lines, `UnitLegendre.cpp` 206 -> 202; the ten largest engine files are
unchanged but `AZUREGrad.cpp` 1253 -> 1234.
**Checks.**  Whole tree, Brune, against `tab-obs4d` (tag `obs4-tree-m5`): 835 of 835 files identical; quick set, the
lab-frame projects, `beam_spread_ay`, `lab_frame_fold` x 3 bases identical to `tab-obs4d`; ctest subset 50/50;
`AZURE_NODE_CHECK` 13N, natural_neon 0 breaches; `scripts/check.sh --quick` on the last commit (see the report).
**Speed.**  thick_target_resonance again, 5 alternated rounds (load 15-17): evaluation median 300.3 ms (`b-trial35`)
-> 302.0 ms (+0.6 %); cross sections 97.5 -> 93.9 ms, matrix fill 22.4 -> 23.4, T-matrix ~150 both: the +4 % of
milestone 4 was noise.
**Open for THM on the mesh** (branch 10).  THM sites are `kValue` keyed by their angle (they hold σ_HOES at the point's
angle), off the mesh with a grid of their own that moves with a shift (the last user of `Unit::CalcNodeRows`' own set
and of the leaves' angle copying).  On the mesh they would need an angle-free record (the HOES amplitudes per node) and
a reader row for the THM angular window, as A_y has; until then the key has no angle and must not be asked for THM.

## 2026-11-12 (agent, Phase 5 branch 9, milestone 4) — a point given at a lab angle is folded in the lab (D5); stopped before milestone 5 on a stop condition (`feature/upgrade-5-observables`)

**Merged first** `feature/upgrade-5-followups-1` (only this file conflicted: both sides kept, newest first by commit
time); quick set x 3 identical to `tab-obs3`, `beam_spread_ay` byte-identical to `tab-obs3`'s `target_effect_ranges`,
`thin_layer_limit` to `tab-small1c`.
**Review remarks.**  (1) Check 37 of the consistency test: every field of `TermSettings` is compared in
`SameObservable` or named in `kReaderOnlyTermSettings` beside it (shown failing with a field added and with `ic`
dropped).  (2) Adjoint solves per gradient (`AZURE_PROFILE=1`, first gradient of a `--use-gradient` fit, 1 thread,
`b-trial35` -> this branch): back-propagations T2 2,959 =, 13N 5,327 = (scalar route; 48,531 = on the Jacobian route),
hybrid_potential 2,472 =; `PointTBar` calls on the scalar route hybrid_potential 112,051 -> 4,884.  (3) not done (the
A_y and normalised-coefficient folds are still named in `Unit::FoldSubPoints` and `GradTargetEffectAdjoint`).
(4) Done: one rule for per-slot factors (`RowSet::perSlot`: angle, P_L, Y, P^1; f_C and J always per slot).
**D5** (`energy_grid.md` 4.3, status there).  A reader whose term's angle is given in the lab
(`experiment::AngleIsLab`) shares a row set with every reader of its lab angle; the set holds theta_cm(E_n; theta_lab)
and J(E_n; theta_lab) per slot (`experiment::LabAngleAt`: `ConvertRow` at `BeamLabFromCM` of the node -- a slot's
"lab energy" field holds its c.m. energy, so it is not read), the contracted value is divided by J_n, the fold is
multiplied by the point's own frame Jacobian (`Unit::LabFoldFactor`).  Reference leaves the same
(`FoldReferenceAtPoint`); fold adjoint and A_y adjoint follow.  `Unit::CalcNodeRows` (a set of one's own) is left to
points with a grid of their own whose slots move with a shift.  A node beyond the kinematic maximum is refused.
**Moved**: whole tree, Brune, against `tab-obs3` (`work/obs4/cls-tree-m4.tsv`): T2 only (7 files), up to 2.56e-5 of a
model value; quick set + lab-frame projects x 3 bases the same.  Row sets: T2 60 -> 2, hybrid_potential 2 = 2.
**Brute force** `tests/lab_frame_fold`: per node 2.2e-5 (quadrature), frozen angle 2.2e-2; milestone 3 fails it.
**Stop condition met** (not decided here): at 6 of T2's 60 points the move exceeds ten times 1/2 Delta^2 |(ln JP)''|
(at most 508x at a point where (ln f)'' nearly vanishes, 2.9e-8 against 5.8e-11); the formula lacks the term
sigma^2 (ln sigma)'(ln f)', which with the second-order term reproduces the 150 deg moves to 0.65-1.12.  Every move is
below 1/1900 of the datum's error.  Milestone 5 not started because of it.
**Checks**: build 0 warnings; consistency; ctest subset 50/50 (predictions.txt: the new project's 4 lines);
check.sh of angle_mapping, convolution_unify, data_covariance, ratio_to_rutherford, relativistic,
shift_subpoint_angle, target_density, lab_frame_fold; AZURE_NODE_CHECK 0 differ (T2, hybrid_potential, lab_fold);
AZURE_GRAD_CHECK as milestone 3 (T2 5.23e-6, hybrid_potential 7.25e-7), lab_fold with a free level 2.2e-6, 2.95e-2 with
the adjoint's J factor dropped (injected once); 1 vs 3 threads identical (T2, 13N, hybrid_potential, lab_fold).
**Measured** (1 thread, 3 alternated runs, medians, menu 1, `b-trial35` -> this branch; units / evaluation ms / wall s
/ peak MB; the machine was loaded, load average 10-13): T2 23,673 -> 18,313 / 50.9 -> 48.7 / 0.30 -> 0.26 / 136 -> 129;
13N 3,825 -> 2,972 / 25.8 -> 23.0 / 8.02 -> 7.99 / 121 -> 119; natural_neon 50,813 = / 418.7 -> 419.4 / 1.34 -> 1.39 /
340 -> 342; thick_target_resonance 28,924 = / 219.6 -> 228.5 (+4.1 %, alone; +5.0 % beside another run; the milestone-3
binary measured 372 against 247 in a third round: noise of this load) / 0.83 -> 0.84 / 271 -> 272; hybrid_potential
4,886 -> 2,473 / 26.7 -> 23.7 / 1.18 -> 1.16 / 106 -> 102; beam_profile_kernel 83,764 -> 3,196 / 467 -> 89 / 369 -> 16.8
/ 891 -> 131; lab_frame_fold 872 = / 1.1 -> 1.2 / 0.14 -> 0.15 / 95 =.
**Open.**  The director's decision on the estimate; milestone 5 (deletions) and remark 3; an A_y at a lab angle is
refused by format 3, so D5's A_y path is not reached by any project.

## 2026-11-12 (agent, Phase 5 branch 9, milestone 3) — shared angular rows, A_y without an angle, sites shared across data sets (`feature/upgrade-5-observables`)

**Row sets** (`include/RowSet.h`).  The readers of one block of angle-free slots that see them under one angle (same
kind, c.m. and lab angle bits, Q-coefficient effect, entrance pair) share one set of rows over the union of their
slots: P_L (one row, or one per slot where the Q transform depends on the energy), f_C per slot for elastic data,
Y_l'^mu and P_k^1 at the angle for A_y.  Built once by `WorkingCopy::BuildRowSets` (set-up; the hand-over when slots
are placed again); a point whose angle moves with a shift gets a set of its own (`Unit::CalcNodeRows`, not cached).
The values a reader folds are contracted once per evaluation into `Results::contracted` (the caller's buffer), in one
pass between the energy pass and the folds (`WorkingCopy::ContractRowSets`, one task per set; `Unit::Calculate` for a
point on its own), and read by every reader (`Site::SlotValue`, `SlotAnalyzingPower`).  The scalar gradient sums the
cotangents of the reads of one slot through one set before calling the point adjoint once (`MergeSharedReads`).
`Site::node_rows_` / `node_coulomb_` are gone.  `CheckReaderRows` checks every reader's set (range, rows, f_C, A_y
factors) bit for bit; injected once (f_C of every slot at the set's first energy): hybrid_potential "112030 differ".
**Records sized by the reaction**: the angular capacity is a function of the compound's largest Legendre order and
amplitude count (`RecordCapacityOf`, `WorkingCopy::Capacity`); the elastic record is packed after its counts; scale and
poison are entries of the table of kinds (no branch per kind left in `ScaleObservable` / `PoisonObservable`).
**A_y** (`SiteObservable::kAnalyzingPower`, layout `ay_record` in `SiteObservable.h`): the site holds the cross
section's record (elastic or Legendre, the formula of every other point) and A_y's coefficients -- for a particle exit
C(a, l') = sum i sqrt(2l+1) cg cg T identical (the AddPathway loop accumulated per exit l',
`AmplitudeMatrix::AddPathwayCoefficients`), for capture the Seyler-Weller d_k, n_k summed per order.  The contraction
assembles the amplitude matrix from C and the row's Y and calls the unchanged `AnalyzingPowerAy()`; the adjoint
builds the matrix at the reader's angle with the reader's f_C at the slot (`AnalyzingPowerBar`, `PathwayAdjoint` with
Y from a `HarmonicTable`).  One kind with two forms (particle, capture), chosen by the exit pair at the writer: the
term settings do not know the exit's type.
**The key** (`SlotGroupKey`): the term settings the energy pass reads (`SameObservable`, TermSettings.h), the mesh, the
target's line (`TargetEffect::GetLineNumber`: the per-segment copies of one target line share its parameters, hence the
slots' stopping), the stopping conversion, frame shift, screening; no segment, no nominal energy.  The hand-over
re-places a shared block once with the readers of every term.  Deleted: `SlotGroupKey::nominalEnergy`,
`Segment::ReangleSlots`, the lab-frame exception of `FoldOnlyShift` (no shift column recomputes node values now).
**What else is per (E_n, r) already**: phase shifts and total sigma are scalars of the reaction (kValue; the key's
settings hold J, l, s); beta-delayed spectra are angle-integrated (kValue); external-capture amplitudes are per site
energy and reaction.  THM stays kValue keyed by angle, off the mesh (branch 10).  kValue now means a scalar at E_n or
THM.
**Numbers.**  Quick set + 13 neighbours x 3 bases against `tab-trial35`: identical as printed except the known
`sum_lab_frame` chi^2 noise, `ec_e1m1_phase` (as milestone 2), run.log counts (`beam_profile_kernel` as before; 13N
external-capture integrals 28,384 -> 21,560: its two capture sets share sites), and `angle_mapping` A_y of 8e-8 at
low energy, 163 numbers by at most 2.2e-10 relative: A_y there is a 1e-7 difference of O(1) amplitude products, so
rounding of 1e-16 in the amplitudes moves it by up to ~1e-9 of itself; the sum over pathways now runs before the
multiplication by Y.  A folded A_y variant of 13N (target on two baumann sets; `work/obs/m3/ayfold`): every output
file byte-identical to `b-trial35`, 49 readers' rows checked.
**Sharing.**  hybrid_potential: meyer_84 and meyer_114 read one block of 2413 slots (before 2 blocks, 4826), numbers
identical; `tests/upgrade/node_check.sh` asserts it (fails with the copy number in the key: "2 blocks, 4826 slots").
**Measured** (1 thread, 3 alternated runs, medians, menu 1, base `b-trial35`; units / evaluation ms / peak MB):
hybrid_potential 4,886 -> 2,473 / 23.3 -> 22.6 / 106 -> 101; T2 23,673 -> 18,313 / 58.9 -> 55.0 / 136 -> 129; 13N
3,825 -> 2,972 / 22.2 -> 21.2 / 121 -> 119; natural_neon 50,813 = / 453 -> 444 / 340 -> 342; thick_target_resonance
28,924 = / 181.7 -> 180.1 / 271 -> 272; beam_profile_kernel (16 data sets of one reaction now share their sites)
83,764 -> 3,196 / 355 -> 82 / 891 -> 131, wall 331 -> 16.5 s; 13N_capture_ay 108 = / 0.9 -> 0.7 / 99 =.  1 against 3
threads: hybrid_potential, T2, 13N, beam_profile_kernel identical (every output file).  (Row sets alone, before the
key change: hybrid 21.9 -> 22.7 ms, 106 = MB.)
The project with the most angles per energy is 13N (38 A_y angles at 1.738 MeV, unfolded; script
`work/obs/m3/scripts/angles_per_energy.py`); unfolded points gain nothing from row sets.
**Gradient** (`AZURE_GRAD_CHECK`, same agreement as `b-trial35`): 13N_capture_ay fd 7.2e-10, folded-A_y 13N 6.7e-8,
hybrid_potential 7.25e-7 (with shared sites), T2 5.23e-6.
**Extensibility.**  (a) A new observable of a reaction (a tensor analysing power): a kind and its layout
(`SiteObservable.h`), its table row (slots, contraction, values, scale, poison) and `ObservableOf`
(`SiteObservables.cpp`), its writer branch (`GenMatrixFunc::CalculateBareCrossSection` + a coefficients function),
its adjoint branch (`AMatrixFunc::PointTBar`); new angular factors in `RowSet.h` / `Unit::RowSetOver`.  More than one
row because a ratio observable's fold (`Unit::FoldSubPoints`) and its fold adjoint (`AZUREGrad.cpp`) are written for
A_y by name, and a term flag needs a field in TermSettings and a quantity in `experiment/`.  (b) A new angular
dependence (a detector's solid angle): the row builder `Site::LegendreRowAt` (with its effect) and the row identity
(`RowIdentity`, WorkingCopyUnits.cpp).  (c) A ratio of two coefficient leaves: the engine reads coefficient leaves
already (`Site::ValueFor`), but `experiment/` refuses it: "quantity 'legendre' with Ratio(...): takes one reaction, not
a ratio" (E_REACTION_NOT_ALLOWED), and one order column serves both leaves, so B_2(R1)/B_0(R2) cannot be written --
a grammar change (a new field), outside this branch.
**Ran**: build 0 warnings; consistency test; `ctest -R "experiment|kernels|grad|upgrade_plan|upgrade_node|upgrade_mesh|
consistency|results_owner|history|ay|polariz|shift"` 48/48; `AZURE_NODE_CHECK` (node_check, 13N) and
`AZURE_PLAN_CHECK` (13N: 0 differ); `scripts/check.sh --quick` (JOBS=2) on the last code commit: 210/210.  Not run:
the whole tree, `run_tests.sh`, the gate, the eight fits.
**Open.**  Milestones 4 (D5) and 5 not started.  The A_y and coefficient-ratio folds are named per observable
(see (a)).

## 2026-11-12 (director, 13:35) — merged: six follow-ups (thin layers, the local shift column, value flags, a rename, two tests)

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (the branch started from the integration head, so the
merge's tree is the branch's plus this entry; nothing pushed): `feature/upgrade-5-followups-1` (the agent's entry
below says what each item does, by file).

**Verified by the director** (build `work/b-small`, reference `tab-trial35` = the shift kernel as merged; the new
reference of the merged engine is `tab-small1c`).
Every project on three bases under a fresh tag (`tab-small1c`; see the note on the runner below): 817 / 795 / 817 of
835 / 803 / 835 files identical to `tab-trial35`, and every other entry belongs to the renamed project (only in the
reference under its old name, only in the new run under its new one; its outputs equal the old folder's, `diff -rq`,
three bases) or to the new project `thin_layer_limit`.  Eight fits (T2 LM, gradient, MIGRAD; 13N gradient, LM;
`layer_yield` gradient; stopping-scale gradient, LM) end at the χ² of `b-trial35` to all printed digits.  Unit tests
211/211, slow tests 10/10, `run_tests.sh` 157/157, converter equivalence 184/184, the converted-tree gate 628 / 601 /
628 identical with none missing, 1 against 3 threads identical on T1, T2, T3, 13N in three bases.
**A defect of the runner met on the way**: `tools/run_one.sh` does not replace an existing output folder of a tag
(`mv` onto a non-empty directory fails and the error is discarded), so the first whole-tree run, under a tag that
had been used at 10:31 for seven folders, compared old files for those and reported them identical; the gate's five
"missing" files exposed it.  Everything above is from the fresh tag.  The fix is in the second batch of follow-ups.

**Reviewed, with the director's own runs** (`work/director/small1_review_notes.md`, `work/director/thin_shift/`).
- *A layer below 1 eV on a data set with an energy shift is still a failed point.*  `thin_limit.azr` at 1e14
  atoms/cm² with `shift offset 0.02 MeV fixed`, `offset 0.02 MeV fit` or `offset 0 MeV fit` on its data set: a warning
  ("shifted fold: the lower edge of the moved support … is outside the point's planned nodes [E, E]") and NaN in the
  output; with `offset 0 MeV fixed` it gives the unshifted bytes, and at 1e16 with a shift it works.  The engine before
  this branch refuses the file altogether (the false "beam does not leave the layer"), so nothing got worse, but the
  item is complete only for unshifted data.  In work: `feature/upgrade-5-followups-2`.
- *The old local shift column moved only the first term of a data set.*  In the Jacobian path (Levenberg–Marquardt)
  the shift column of a segment moved the points of term 0 only; `MoveSegmentsToShift` moves every term, as the
  hand-over does.  A summed data set (several reactions) with a shift of its own had a wrong shift column for its other
  terms before.  No project of the tree has that combination, so no number moved and the old error is not measured.
- *Two slips of the agent, both reported by itself and without consequence*: about six threads for three minutes
  (budget four); a flag test run against the base build nearly wrote `~/.azure2/coulomb_table.bin` (stopped by its
  60 s guard; the file's date is unchanged, 2026-10-26; the test now redirects `HOME`).

**For the PI and the owner** (behaviour a user sees): a value flag without its value (`--max-calls`, an empty
`--coulomb-table=`, a bad `AZURE_REFINE_QUADRATURE`) is an error with the usage line and exit code 2, where it was a
warning or was ignored; a yield on a layer thinner than 1 eV is σ(E)·N, where the file was refused; the project
`tests/target_effect_ranges` is `tests/beam_spread_ay`.

**Open** (carried).  The second batch of follow-ups (`feature/upgrade-5-followups-2`): the thin layer with a shift, an
unknown option as an error, the finite-difference fallback's gates one by one.  The observables branch
(`feature/upgrade-5-observables`): milestones 1–3 checked by the director on the whole tree (`tab-obs2`, `tab-obs3`),
milestone 4 (the lab-frame fold per node) and 5 (deletions) in work.

## 2026-11-12 (agent, follow-ups) — six small items after the shared mesh and the shift kernel (`feature/upgrade-5-followups-1`)

**Done** (one commit each, on 7ba01de7).
- *A target layer thinner than 1 eV.*  A support narrower than `MeshZeroWidth()` is one node, and every kernel returned
  0 for one node (a yield on such a layer was even refused at set-up as a layer "the beam does not leave").  The fold
  of one node is now the thin limit, in one place (`Compose` in `experiment/src/Kernels.cpp`): the node's value times
  each kernel's `Kernel::thinLimit` (1 for the averages; for the target kernel 1 / 1e-24 N / 1e-24 T/ε for a thin
  target / a layer yield / a thick target; none for the beam profile, whose support is at least 1 keV: one node would
  be NaN).  `src/TargetLayer.cpp` accepts one node.  Tests: `layer_kernel` case 6; `tests/thin_layer_limit/check.sh`
  (33 areal densities 1e13–1e17: Y/(1e-24 N) and the averaged σ follow σ(E0) within T|σ'|/σ, equal it to 1e-9 with one
  node, continuous across 1 eV).  New project in the tree: `tests/thin_layer_limit/models/thin_limit.azr`.
- *`pyazr_shift_fn` under `AZURE_NODE_CHECK=1`.*  A defect of the Jacobian's local shift column (`FillShiftColumns`,
  `src/AZURECalc.cpp`): it moved the points of its segments without looking their slots up again, so the point pass
  grouped a moved unit with the energy it had left.  Numbers were right only because the units sharing a node moved
  together.  `WorkingCopy::MoveSegmentsToShift` (`src/NodeTable.cpp`) does what the hand-over does.  Test
  `pyazr_shift_fn_node_check`.  Jacobians bit-identical to the base build.
- *`pyazr_shift_fold_ec`* compares the χ² of a shift moved from 0 to ±0.2, 0.05, 0.02 MeV with the same shift written
  in the file, fixed and free (1e-12).
- *Value-taking flags.*  `--max-calls` without `=N`, a bad `--max-calls=`, `--refreeze-covariance=`, an empty
  `--coulomb-table=` or `--build-coulomb-table=`, and an empty or bad `AZURE_REFINE_QUADRATURE` are errors with the
  usage line, exit 2 (`include/CommandLineValues.h`; pyazr raises ValueError for the variable instead of ending the
  process).  Consistency check 36.  Test `tests/cli_option_values/check.sh`; `tests/fit_call_limit` case 4 changed
  (it pinned the old warning).
- *Rename.*  `tests/target_effect_ranges` is `tests/beam_spread_ay` (whole-tree folder
  `.res-target_effect_ranges/target_effect_ranges.azr` → `.res-beam_spread_ay/beam_spread_ay.azr`).
- *`park_chi2_fallback`* has a third variant that takes the finite-difference fallback of the gradient: a folded
  capture analyzing power (13N_capture_ay with a beam spread), which the adjoint refuses.

**Verified** (scratch `work/small/`): every commit builds with 0 warnings; consistency test; the unit tests of each
item; T1, T2, T3, 13N and each item's projects in three bases byte-identical to `tab-trial35`; each new test shown
failing (base build or an injection).  `scripts/check.sh --quick` (JOBS=1) at the end: consistency, build, `ctest -LE
slow` 211/211 passed (2185 s).

**Not run**: the whole tree × 3, `run_tests.sh`, the converted-tree gate, the fits (the director's).

**Open.**  `tests/upgrade/predictions.txt` gained the new project's entry and the renamed path.  The other routes of
the gradient's finite-difference fallback are listed in the report (formal-basis external-capture blocks and
`--ignore-externals` projectors, `--use-rmc`, a few observable gates in `AMatrixFunc::PointTBar`); not tested one by
one.
An unknown option is still only a warning.

## 2026-11-12 (agent, Phase 5 branch 9, milestone 2) — elastic scattering: a site holds the elastic amplitudes, the reader applies f_C and its row; one table of observable kinds (`feature/upgrade-5-observables`)

**Merged first** `feature/upgrade-5-shift-kernel` (e97e06c4: the translated stopping's clamp, the external-capture
branch of `ReplaceSubPointGrids`); quick set against `tab-trial35` identical as printed but the `sum_lab_frame` noise.
**One table.**  `include/SiteObservable.h`: `enum class SiteObservable {kValue, kLegendre, kElastic}`, `ObservableOf`
(the one classification of a term, replacing the predicate `Site::HoldsLegendreCoefficients`) and `kSiteObservables`
(`src/SiteObservables.cpp`): per kind the slots of the record's angular part and the contraction a reader applies
(`Site::ValueFor` is one line through the table).  `Site::AngleFree()` (kind != kValue) is what the sharing key, the
slots' set-up, the leaves, the shift and the gradient ask.  The record functions moved there from `Unit.cpp`.
A new kind: its row, its branch in `GenMatrixFunc::CalculateBareCrossSection`, its adjoint in `PointTBar` (recipe in
`CLAUDE.md`, "An observable of a reaction at a mesh energy").
**What an elastic site holds** (kElastic: differential elastic data; elastic angular distributions stay kLegendre):
B_L^NN (b/sr, s = 4 for identical pairs, from `LegendreCoefficientsOf`), c_0 = geom/100, c_x = 2 c_0 (-1)^(2j)/(2j+1)
(identical pairs with spin), and D_l, X_l (complex) = (i geom f_it / (100 sqrt pi)) sum_s [1, (-1)^s] c_{s,l} from the
Coulomb-nuclear sums `rt_it_` (`GenMatrixFunc::ElasticCoefficientsOf`; layout `elastic_record` in the header; 108
slots).  dsigma/dOmega = c_0 (|f_d|^2 + |f_x|^2) + c_x Re(f_d f_x*) + sum B_L P_L + Re sum (f_d D_l + f_x X_l) P_l.
No angle.  **The reader** holds f_C: its own (`Site::CoulombAt(0)`), or per slot `Site::node_coulomb_` beside
`node_rows_` (`Unit::CalcNodeRows`, at the reader's angle and each slot's energy; recomputed where `node_rows_` is: set-
up, the hand-over's moved tasks only, an angle moved by a shift); `Site::CoulombAmplitudesAt` is the one formula (also
`CalcCoulombAmplitude`'s).  A slot without an angle computes no Coulomb amplitude and no Legendre rows.  The Rutherford/
Mott leaf takes the point's angle at each slot's energy (`CoulombReference.cpp`).  kValue differential sites (A_y,
THM) contract the same coefficients with their own row on a local record (`Site::ContractRecord`): one formula.
**Deleted** in `GenMatrixFunc.cpp` (-199/+61 lines with the new writer): the per-site differential sum — the
per-interference loop with the unreachable `angularCoeff` accumulation (118 lines), the per-site CT/RT/IT and its
uncached Coulomb branch, `rt_it_built_`.  `Unit.cpp` -56 (moved).  Folders: src +353/-317, include +129/-22 (new files
included).  **Still keyed by angle**: A_y (and THM, off the mesh).  Their users keep `SlotGroupKey::nominalEnergy`,
`Segment::ReangleSlots` and the solving shift column (`FoldOnlyShift` false): lab-frame folded A_y with a shift
(13N's baumann, electron_screening).  `legendreOrder` of the key is now always -1 (coefficient sites are kLegendre).
**A test that catches a wrong row**: `AZURE_NODE_CHECK` now checks every reader of angle-free slots: its row and
Coulomb amplitudes at each slot equal `LegendreRowAt`/`CoulombAmplitudesAt` at the slot's energy, bit for bit
(`WorkingCopy::CheckReaderRows`, a line per check; `tests/upgrade/node_check.sh` adds hybrid_potential, folded elastic).
Injected once (rows and f_C at the point's energy in `CalcNodeRows`): T2 fit "60 readers, 8132 slots, 8132 differ",
hybrid_potential "112032 differ", both abort; removed.
**Ran** (build 0 warnings; `b-obs`): consistency test; `ctest -R "experiment|kernels|grad|upgrade_plan|upgrade_node|
upgrade_mesh|consistency"` 35/35 (node check as above, AZURE_PLAN_CHECK in upgrade_plan_runtime); quick set + target_
density, angle_mapping, hybrid_potential, electron_screening, data_covariance, ec_e1m1_phase, beam_profile_kernel x 3
bases against `tab-trial35`: identical as printed except `sum_lab_frame` (chi^2 4e-17 noise) and `ec_e1m1_phase`
(the named cancellation: 12 numbers, 8 Park, -1.5e-25 -> -2.1e-25 b/sr at 0/180 deg beside 1.6e-9, as milestone 1);
`beam_profile_kernel` run.log: 1,089,452 -> 42,068 capture integrals (milestone 1's gain).  AZURE_GRAD_CHECK against
`b-trial35`, identical chi^2 and agreement: 13N elastic only (meyer 84/114/144, 14 free; fd 5.86e-7, Jacobian 8.9e-16),
hybrid_potential (folded elastic, 14 free; 7.25e-7, 5.2e-16), T2 (5.23e-6), 13N (6.16e-8), layer_yield/models fit
(3.98e-9), 13N_capture_ay (7.19e-10).  1 vs 2 threads: T2, 13N, hybrid_potential, ratio_to_rutherford identical.
**Measured** (1 thread, 3 alternated runs, medians, menu 1, base `b-trial35`; units / evaluation ms / wall s / peak MB):
T2 23,673 -> 21,085 / 46.4 -> 46.2 / 0.27 -> 0.27 / 136 -> 133; 13N 3,825 = / 19.4 -> 20.3 / 6.92 -> 6.89 / 121 = ;
natural_neon 50,813 = / 375.8 -> 368.2 / 1.21 -> 1.20 / 340 -> 342; thick_target_resonance 28,924 = / 151.7 -> 151.6 /
0.56 = / 271 -> 272; ratio_to_rutherford aC_data 56 = / 0.1 -> 0.2 / 0.06 = / 95 -> 93; hybrid_potential 4,886 = /
21.6 -> 26.0 (+20 %; first round 20.0 -> 26.2) / 1.01 = / 106 -> 112.  The hybrid_potential loss: its 76 readers
already shared one angle (4,826 slots, 112,032 reads), so a slot held one finished number; now each read contracts
the record (B_L, D_l with the reader's row and f_C) in the fold and again in the model pass, and the record (108
slots per site) and the per-slot f_C (3.6 MB) are new.  13N gains nothing: its 38-angle sets are A_y (milestone 3),
its elastic data are excitation functions at three angles, unfolded.
**Open.**  A memo of the contraction per (block, reader row) would remove the hybrid_potential loss where readers
share an angle.  Milestone 3 not started (see the report): capture A_y has the Legendre cut (Seyler-Weller a_k, b_k);
particle A_y (the amplitude matrix, `Polarization::AmplitudeMatrix`, with f_C) has no per-order form in the code.
`scripts/check.sh --quick` (JOBS=1) on 879f5f68: 210/210 tests pass.  Not run: the whole tree, `run_tests.sh`, the gate,
the eight fits, 1 vs 3 threads (the budget is 2).
**Next (milestone 3).**  (1) Key without the segment: `MeshGroupKey` holds `term.segment`; what the writer reads from
the settings must replace it (reaction pair keys, UPOS secondary L, I_c, delta, the E1/E2 component, THM/phase/total are
kValue).  (2) Capture A_y: a kind holding B_L and the Seyler-Weller numerator per order (`CaptureAyTerm` ak, bk), reader
P_L^1 rows.  (3) Particle A_y (`Polarization::AmplitudeMatrix`, Y_lmu at the point's angle, f_C) has no per-order
form: a P_L^1 expansion would be a second derivation (a stop condition of the brief); it stays kValue and keeps
`nominalEnergy`, `ReangleSlots` and the solving shift column for its lab-frame shifted segments.

## 2026-11-12 (director, 10:20) — merged: folded data with an energy shift on the fixed nodes of the shared mesh

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (one merge commit, tree of `rehearsal/trial35`; nothing
pushed): `feature/upgrade-5-shift-kernel` (the agent's entry below says what it does, by file).  Every folded point is
now on the shared mesh except THM and the uniform own grids; `NodeTableCache.h` and its locks are gone.

**Verified by the director** (build `work/b-trial35`, reference `tab-trial33` = the D3 deletions as merged).
All projects × 3 bases: 802 / 775 / 802 of 823 / 796 / 823 files identical; the others are T2, T2_true and the five
models of `shift_sqrt_subpoints`; classifier per basis 13 quadrature, 1 named (T2 `covariance.out`), 0 failure
(`work/director/clsq_trial35.txt`).  ctest 210/210 and slow 10/10; `run_tests.sh` 155 passed, 0 failed; spec
equivalence 184 of 184 converted equal; converted-tree gate 623 / 596 / 623 identical; 1 against 3 threads identical.
Eight fits at 1 thread: 13N (two routes), the layer yield and the stopping yield (two routes) byte-identical to the
engine before; T2 ends at χ² 297.598987 (LM), 297.598987 (gradient), 297.599950 (MIGRAD, 200 calls), byte-identical
before and after the review fixes, with no "unsupported imaginary number" line in the logs (14,711 / 9,046 / 4,790
before).  The review fixes alone move the five `shift_sqrt_subpoints` models by at most 4.5e-6 (one χ²) and 8e-9
(yields); everything else of the tree is byte-identical before and after them.

**Reviewed, with the director's own measurements** (`work/director/shift_review/`).
- *The fit of T2 moves by 0.13 of an error* (sqrt shift of the thick target +0.131 σ, offset of the thin target
  +0.11 σ; the brief's stop was 0.1).  At equal stopping the new fit agrees with the fit of the old engine refined
  eightfold within 0.010 σ, the old engine within 0.131 σ: the move is toward the converged quadrature.  Accepted.
- *T2's `covariance.out`* (the correlated χ² of the two lab-frame segments) fails the per-file rule: 2.63e-4 from the
  refined value against 2.28e-4 before.  Per point both segments are 3–4 times closer than before (largest 2.9e-4 /
  3.0e-4 → 7.5e-5 / 1.2e-4); the total is dominated by one point (E 1.1228 MeV, pull +15), and the old total was
  nearer only because other points' errors cancelled part of it.  A named row of the classifier (factor 1.2).
- *The T2 gradient fit took 456 calls where it took 118.*  Not a slower convergence: the gradient agrees with central
  differences at moved shifts and at the minimum, χ² is smooth along the three shift parameters, and with the start
  value of one normalisation changed (1.002, 0.997, 1.01, 0.99) both engines end either after 117–118 calls (EDM
  about 1e-4) or after 452–456 (EDM about 7e-8); at 0.99 the old engine takes 454 and this one 118.  The fit sits on
  the edge of Minuit's validity check after HESSE.  **For the PI and the owner**: the call count of this fit is not a
  measure of the engine.
- *A step at a re-plan*: when a shift leaves the reach its nodes were planned for, the mesh is planned again and χ²
  steps by the quadrature difference: 1.3e-2 in a χ² of 10,500 at +9.0e-3 MeV on T2's offset with prior 0.003 (three
  prior widths).  None inside the reach.  The old engine had steps of the same size inside the range a fit visits
  (1.7e-2 at 1.16e-3 MeV on the same parameter).
- *Two defects found in the review, fixed by the agent before the merge, each with a test that fails before the fix*
  (repeated by the director on the builds before and after): the stopping formula was called at a negative energy at
  the floor node of a thick target with a shift above 1 keV (thousands of "unsupported imaginary number" lines in a
  T2 fit; the argument is now held at the support's lower edge; T2's outputs do not move); a folded data set with
  external capture whose shift left the reach stopped the fit with an error (its nodes are now re-placed and take their
  capture amplitudes from the tables).  For the second the director compared, on 13N `skowronski_fels` with a beam
  spread, the χ² at a shift moved from 0 to +0.2, −0.2, +0.05, +0.02 MeV with the same shift written in the file
  (fixed and free): equal to all printed digits in the four cases (`work/director/shift_review/ec_far.out`).

**What it costs and gives** (the agent's measurements, T2, one thread): level-matrix solves of one gradient 158,929 →
48,526; peak memory 300 → 259 MB; one calculation 0.23 → 0.28 s (the reach adds nodes).  Engine lines: `src/` +408
−165, `include/` +193 −144 (the engine grew: the fold of a shifted point is a named file, `ShiftedFold.h/.cpp`, 312
lines, against the deleted cache of 111).

**Open** (in the order they will be met).  The lab-frame differential data with a shift still carry the point's angle
in the site (`SlotGroupKey::nominalEnergy`, `Segment::ReangleSlots`, a shift column that solves): branch 9 removes
them.  THM and uniform own grids keep the rigid move (`Segment::MoveOwnGridSlots`): branch 10.  The target kernel of a
shifted point reads the stopping at the node's energy minus the shift (the old rigid move's physics, kept so that
numbers stay): whether the stopping should be read at the true node energy is a physics question for the owner.  The
shift columns of the gradient are central differences of the fold on fixed nodes; an analytic ∂K/∂E_0 per kernel is
not written.  `pyazr_shift_fn` fails under `AZURE_NODE_CHECK=1` (`residual_jacobian` raises) on the engine before this
branch too; not investigated.  `pyazr_shift_fold_ec` checks that χ² is finite and returns; the comparison with the
exactly planned shift (above) should become one of its checks.  Carried: the silent zero of a layer below 1 eV, the
near-copy loops of `SetUpTargetEffects`, `park_chi2_fallback`, the name of `tests/target_effect_ranges`.

## 2026-11-12 (agent, Phase 5 branch 9, milestone 1) — non-elastic angular data: a site holds Legendre coefficients, the point's row applies the angle (`feature/upgrade-5-observables`)

**What a site holds now.**  A site of a differential term of a reaction whose exit pair is not its entrance pair,
or of any angular-distribution / Legendre-coefficient term (`Site::HoldsLegendreCoefficients`), holds in the angular
part of its record the coefficients B_0 .. B_{K-1} (b/sr) of dsigma/dOmega_c.m. = sum_L B_L P_L(cos theta) at its
energy (K = the exit's `num_orders`, at most `Config::maxLOrder` + 1 = 21 slots), and for an angular distribution
A_0 = sigma_int / 4 pi.  No angle.  The cut is where it was in today's sums: `GenMatrixFunc::BuildOrderProducts`
(the per-order products Z T1 T2* of the pathway table that the differential sum already summed once per T-matrix,
`rt_wl_`), scaled by geom I1I2 s / (100 pi) in `GenMatrixFunc::LegendreCoefficientsOf`; the UPOS weights
(Z_upos sqrt(2L+1)/4 R_L) are per order too.  No second formula: the external-capture interferences (ER, RE, EE)
are in the same per-order list.  s = `LegendreSymmetryFactor` (4 elastic identical, 2 another exit of an identical
entrance, else 1); a_k A_0 = B_k / s (`Site::CoefficientOfOrder`).
**The row.**  A point reads a site through `Site::ValueFor(reader, slot, results)`: sum_L B_L RowP(slot, L) for
differential data, B_k / s (absolute) or B_k / (s A_0) (normalised) for a coefficient datum of order k.  The row is the
point's own Legendre rows (with its Q coefficients), or, where the detector's lab transform of the Q coefficients
makes the row depend on the energy (particle exit, not UPOS / A_y / charged elastic: `Site::RowDependsOnEnergy`), the
point's rows at the energies of its slots (`Unit::CalcNodeRows`, `Site::node_rows_`, computed from
`Site::LegendreRowAt`, the former body of `CalcLegendreP`).  They are recomputed when the point's angle moves with a
shift (`ShiftFrames.cpp`) and for the points whose slots were placed again (`HandOver.cpp`,
`WorkingCopy::CalcNodeRows`).  Readers: the fold (`PointFold.cpp`), the coefficient and A_y folds (`PointPass.cpp`),
the transform terms (`ModelPass.cpp` `LeafValue`), the model record's a_k (`Unit::CopyAngularDistsToModel`), the
gradient (`UnitCotangent::rows`, `rowSlot`; `AMatrixFunc::PointTBar(..., rows, rowSlot)` takes the Legendre rows and
coefficient order from the reader).
**The key** (`MeshGroupKey`) has angle 0, Legendre order -1 and no nominal energy for these groups; slots get no
Legendre rows and no order (`BelongToRowOf`); `ReangleSlots` is not called for them and the shift column of a
lab-frame segment of them refolds on the fixed nodes (`FoldOnlyShift`); leaves of them are not separated by angle.
`CheckSlotGroups` knows the angle-free key.  Electron screening scales the coefficients (`ScaleLegendreCoefficients`).
A failed site has NaN coefficients (`MarkFailed`).
**Deleted / unreachable.**  The record write of normalised a_k (`SetAngularDists`) is gone; the angular-distribution
part of the old interference loop (`angularCoeff`) is unreachable and goes with milestone 5.  Model records size
from `ModelAngularCapacity` (unchanged), site records from `ResultAngularCapacity` (the coefficients).
**Largest difference.**  Quick set (T1, T2, T3, 13N, natural_neon, thick_target_resonance, layer_yield/models,
angular_coefficients, 13N_capture_ay, ratio_to_rutherford, sum_lab_frame; three bases) against `tab-trial34`: every
file identical as printed except `sum_lab_frame/models/{sum_data,sum_free}` chiSquared.out / run.log, chi^2 =
3.99e-17 (data equal to the model, rounding noise) moving by up to 5.0e-6 of itself (2e-22 absolute), the case the
lab-frame branch listed.  A first version that read the point's row at its own energy instead of the node's moved
T2 s02/s03 (Q attenuation on 15N(p,a)) by 5.1e-7 relative in 119 numbers: the per-node row above removes it.
**Ran.**  Build 0 warnings; consistency test; `ctest -R "experiment|kernels|grad|upgrade_plan|upgrade_node|
upgrade_mesh|consistency"` (all pass, upgrade_plan_runtime and upgrade_node_check run AZURE_PLAN_CHECK=1 and
AZURE_NODE_CHECK=1); `check.sh` of angular_coefficients, ratio_to_rutherford, sum_lab_frame pass; AZURE_GRAD_CHECK on
T1, T2, T3, 13N, layer_yield/models/fit, 13N_capture_ay: the same chi^2, fd and Jacobian agreement as the base build.
Injected once: every coefficient reader reading order 0 -> angular_coefficients fails 5 checks (data, shift, conv,
gradient); removed.
**Measured** (1 thread, 3 alternated runs, medians, menu 1, against `b-trial34`): T2 units of the energy pass
23,673 -> 21,085, evaluation 60.4 -> 57.5 ms, wall 0.36 -> 0.34 s, peak 136 -> 133 MB; 13N, natural_neon,
thick_target_resonance: units unchanged, wall and peak within 1 %.
**Not run.**  The whole tree, `run_tests.sh`, the gate, the eight fits, 1-vs-3-thread comparison (1 vs 2 threads,
the budget: T2, 13N, angular_coefficients menu 1 identical).
**Stopped here** (context of the session): milestones 2-5 not started.  Still keyed by angle: elastic differential
data (Coulomb amplitude and Legendre rows on the site), A_y data, phase shifts, total, THM.  Their users keep
`SlotGroupKey::nominalEnergy`, `Segment::ReangleSlots` and the recomputing shift column (`FoldOnlyShift` false):
lab-frame differential elastic segments with a shift.  Next: milestone 2 (elastic: B_L^NN with s = 4 for identical
pairs is already what `LegendreCoefficientsOf` gives for an elastic angular distribution; the Coulomb-nuclear terms
rt_it_ per (s, l) go in the record, f_C per reader and slot beside `node_rows_`).

## 2026-11-12 (agent, Phase 5 branch 8) — folded data with an energy shift on the fixed nodes of the shared mesh (`feature/upgrade-5-shift-kernel`)

**Done.**  A folded point of a data set with an energy shift folds on the shared mesh like every other folded point:
the nodes stay, the kernel moves (design: the comment of `include/ShiftedFold.h`, summary in `energy_grid.md` §2.5).
The planner gets the support moved by a centre and widened by a reach (3 prior widths of the shift parameters; 0 for a
fixed shift); a shift leaving the reach re-plans the nucleus in the hand-over (a function of the shift alone, a Note in
`run.log`).  The fold reads the sub-range of the point's planned slots enclosing the moved support, clips there, and
gives a target kernel ε(E_n − Δ) (the old rigid move's physics).  Lab-frame differential groups take the point's moved
c.m. angle (Legendre rows, Coulomb amplitude).  Shift columns of the gradient and residual Jacobian fold again on the
fixed nodes (no solve) where no site depends on the shift.  Deleted: `NodeTableCache.h`, `SharedTable::Weak`, the
shift conditions of `OnSharedMesh`, the rigid move of mesh slots (kept for THM and uniform grids:
`Segment::MoveOwnGridSlots`).  Fixed on the way: `Unit::CalcLegendreP` on a shifted point appended to its slots' rows
at every move.

**Verified** (scripts `work/shift/`): build 0 warnings; consistency test; ctest of experiment, kernels, grad, upgrade,
consistency, shift, shared_point, parallel_points, jacobian (39/39); `upgrade_node_check` (adapted: the shifted
calculation uses no overflow node); T1, T3, 13N, natural_neon, thick_target_resonance, layer_yield/models byte-identical
to `tab-trial33` in three bases; T2 and shift_sqrt_subpoints classified (`energy_grid.md` §2.5: all per-file quadrature
but T2 `covariance.out`); `check.sh` of convert_params, energy_shift_sqrt, shift_fn, shift_groups, shift_sqrt_subpoints,
shift_subpoint_angle ((b) back to 1e-6, measured 0), target_density, tests/upgrade (golden re-pinned) pass;
`AZURE_GRAD_CHECK=exit` on T2 and on shift_sqrt_subpoints sqrt and shiftfn with the shift varied: shift components
agree with the whole-χ² difference to 1e-12 or better.  Fits (T2, 200 calls, 1 thread, 3 alternated runs): the
fitted shifts move 0.131 and 0.114 of their errors (more than 0.1: a stop condition of the brief), to the minimum of
the engine before refined eightfold (sqrt shift: base 1.8638e-3, new 1.9533e-3, refined 1.9538e-3; offset 2.1977e-3,
2.2275e-3, 2.2300e-3).

**After the director's review (same day).**  Both stop conditions decided by the director: T2's `covariance.out` is a
named row of the classifier (factor 1.2; reason in `energy_grid.md` §7.1), and the 0.13 σ move of the fit is the
quadrature's improvement (the new fit agrees with the refined engine's within 0.010 σ).  The 118 against 204-456
calls of the T2 gradient fit are not a slower convergence: with one normalisation's start value changed (1.002, 0.997,
1.01, 0.99) both engines end either after 117-118 calls (EDM about 1e-4) or after 452-456 (EDM about 7e-8), and at
0.99 the base takes 454 and this engine 118; the fit sits on the edge of Minuit's check after HESSE.
Fixed: (A) the translated stopping of a thick target with a positive shift was evaluated below the floor (E_n − Δ < 0
at the outer node of s04_thick: 4790 "unsupported imaginary number" warnings in a T2 gradient fit, 2275 in the LM fit,
17 in one calculation of s04 at sqrt shift 0.002); its argument is held at the nominal support's lower edge, 0
warnings, every output file of both T2 fits byte-identical to before the fix (`pyazr_shift_fold_floor`, fails before
the fix).  (B) a folded segment with external capture whose shift left the reach threw mid-fit; its blocks are now
re-placed and take their amplitudes from the tables (`pyazr_shift_fold_ec`: 13N skowronski_fels with a beam spread and
a free shift without prior, +0.2 MeV = 3.4 reaches: χ² 4874462.5, base 4874562.0, base refined eightfold 4874427.9;
back at the file value to 1e-12).

**Open.**  T2 `covariance.out` and the 0.13 σ fit move are for the director (both explained, both stop conditions of
the brief).  Analytic ∂K/∂E_0 per kernel; the lab-frame differential column still solves (branch 9); THM keeps the
rigid move (branch 10).  Not run: the whole tree × 3 bases, `run_tests.sh`, the gate, the eight fits, other thread
counts.

## 2026-11-12 (director, 08:10) — merged: fold windows, edge blending, the automatic tolerance and the σ(Ē) fallback are deleted

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (one merge commit, tree of `rehearsal/trial33`; nothing
pushed): `feature/upgrade-5-d3-deletions` (the agent's entry below says what went, by file).  A folded data set now
names its effect and, optionally, the settings of the mesh; every point of it is folded with the whole kernel.

**Verified by the director** (build `work/b-trial33`, reference `tab-trial31` = the shared mesh as merged).
All projects in three bases: 820 of 823, 796 of 796, 820 of 823 files identical; the three that differ per basis are
`target_effect_ranges`' observable file, `chiSquared.out` and `run.log`.  `ctest -LE slow` 207/207 (1010 s), slow
10/10; `tests/run_tests.sh` 155 passed, 0 failed (2817 s; the two deleted checks gone); spec equivalence 215 files,
184 converted, 184 equal; the converted-tree gate against the rehearsal's own run 623/596/623 identical; 1 against 3
threads identical; the eight fits at 1 thread byte-identical to the shared mesh as merged
(`work/director/fits_t33_t31`).
Speed and memory are the agent's measurement (equal to a base built the same way), not repeated by the director.

**Reviewed.**  The engine's diff (`src/PointPass.cpp`, `src/AZUREGrad.cpp`, `src/TargetSetup.cpp`,
`experiment/src/KernelTarget.cpp`) is deletion only; no branch of the fold or of the adjoint was rewritten.

**The one moved project, accepted.**  `tests/target_effect_ranges` (a Gaussian beam spread folded into A_y data inside
a window 1.95–2.55 MeV, blend 0.12 MeV, automatic tolerance 0.002): χ² 3838.66 → 3760.95; outside the old window the
largest change of A_y is 0.046 (1.5 of the datum's σ), inside 0.004.  This is larger than the brief's stop condition
(the old tolerance), so the agent stopped and asked.  Repeated by the director: the engine **before** the deletion on
the old file with the window widened to 0–1000 MeV gives both output files of the new engine bit for bit
(`diff -rq`), so the new numbers are the old engine's full fold and the move is entirely the points the window left
unfolded.  The window was a restriction of the physics (a point outside it was compared without the beam spread), not
an approximation with a bound; decision D3 removes it.  **For the PI and the owner**: an old file with windows
converts with a note that names the line, and its numbers outside the windows change by the size of the fold there.

**Found beside it, not caused by it** (present in the engine before the mesh, `work/b-trial30`, with the same
numbers): a layer whose energy loss is below 1 eV (`MeshZeroWidth()`; `thick_target_resonance` with
`areal_density 3.2e+14`) is a single node, and every kernel returns 0 for fewer than two nodes — a silent zero where
the thin limit σ(E)·N is meant.  At 3.2e+15 (1 eV) and above the engine before the mesh and this one agree to 3e-5
or better.  To do: the single-node case returns the node's value times the kernel's norm, or the project is refused;
a test of the thin limit.

**Corrections.**  The entry below was headed 07:15; its commits are of 07:03 (header corrected).  `FOR_THE_OWNER.md`
said `dev`'s
automatic target-effect windows "port as they are": with D3 the two commits (442610f3, 3a104e06) are obsolete, and
the `azure2-eval` skill's notes on windowed refits will need a sentence at the merge with `dev`.

**Open** (carried, in the order they will be met).  The shift kernel (`feature/upgrade-5-shift-kernel`, started from
this tree): shifted folded segments on the fixed nodes, `shift_subpoint_angle` (b) back to 1e-6, `NodeTableCache.h`
and the rigid move deleted.  The two loops of `SetUpTargetEffects` are near copies (merge them when the placement of
slots is next touched).  `park_chi2_fallback` no longer exercises the finite-difference fallback of the gradient: a
construct the adjoint does not cover is needed for that test, or the fallback goes.  `tests/target_effect_ranges` has a
historical name (it is now "a beam spread folded into A_y").  The acceptance rule of changed folded numbers
(`energy_grid.md` §7.1)
was **confirmed by the PI on 2026-11-12**; the owner has still to see it (`FOR_THE_OWNER.md`).

## 2026-11-12 (director, 07:03) — merged: the engine on the shared mesh, the shared sites, the expectations the mesh moves

**Merged** into `feature/upgrade-4-data` and `feature/upgrade` (one merge commit, tree of `rehearsal/trial32`; nothing
pushed): `feature/upgrade-5-shared-mesh` (the engine places a folded point's nodes from one plan per compound nucleus;
the per-file rule of the classifier; the regenerated expectations) and `feature/upgrade-5-shared-sites` (the folded
points of a sharing group read one `Site` per mesh node).  **This is the first merge of the upgrade that changes
numbers on purpose.**

**Verified by the director** (build `work/b-trial31`; whole-tree tag `tab-trial31`, **the reference for the next
branches**):
- all projects x 3 bases against the mesh engine alone (`tab-mesh2`): every output identical (823 / 796 / 823 files;
  the one difference per basis is `13N`'s `run.log`: 28384 external-capture integrals instead of 367096) — the shared
  sites change no number;
- against the engine before the mesh (`tab-trial30`) with the refined references (`tab-fine8,tab-fine4`): Brune 577
  files, 542 identical, 33 quadrature, 2 named, 0 failures; no-brune 550 / 517 / 31 / 2 / 0; use-park 577 / 542 / 33
  / 2 / 0 (`work/director/clsq_trial31.txt`; per-project table now in `energy_grid.md` §7.1);
- `ctest -LE slow` 207 / 207, `ctest -L slow` 10 / 10; `tests/run_tests.sh` with the regenerated expectations:
  157 passed, 0 failed (with the old expectations: the six failures the expectations agent was given, no other); spec
equivalence
  184 of 184 converted equal; the converted-tree gate against the engine's own whole-tree run 623 / 596 / 623
  identical, 0 classified, 0 differ, 0 missing (**the gate's reference is now the tag of the engine under test**,
  `work/director/full_verify5.sh WT TAG RUNTAG [THR] [GREF]`: the converted tree must reproduce the originals with
  the same engine; the comparison with the engine before is the classifier's); 1 against 3 threads identical.

**Review of the shared sites** (director, the whole of `include/SlotGroup.h` and the diffs of `MeshSlots`,
`CalculationPlan`, `EnergyPass`, `PointPass`, `Unit`, `WorkingCopyUnits`, `AZUREGrad`, `GradEnergyTasks`): no defect
found.  The key of a sharing group holds everything a site's value depends on besides its energy; a site is written
by its first reader in data order and by nobody else; a failed site fails every reader; the gradient sums a shared
site's cotangents in the readers' plan order, so the sum does not depend on the thread count.

**Fits** (eight fits, 1 thread, `work/director/fits_t31_mesh2`, `fits_t31_t30`): byte-identical to the mesh engine
alone.  Against the engine before the mesh the fitted parameters move, in units of their fitted error: T2 0.01, 13N
0.000, `layer_yield` 0.03, `stopping_yield` 0.14 (gradient route; 0.06 by the least-squares route's error; areal
density 9.00435e17 -> 9.00498e17, 7e-5 relative).  Direction checked on that fit with the engine before the mesh
refined eightfold, started near the minimum: areal densities 6.00387e17 / 9.00435e17 (before), 6.00426e17 / 9.00498e17
(mesh), 6.00449e17 / 9.00533e17 (refined); stopping scale 1.0292742, 1.0292093, 1.0291722; chi2 72.04, 71.67, 71.41
— all three parameters move towards the refined fit.  **Side finding**: the same refined fit from the test's far
start ended invalid (chi2 154618 after 67 calls, "MIGRAD did NOT converge"); a far start's path depends on the
quadrature, the minimum does not beyond the numbers above.  Not a matter of the mesh; for the fit-reproducibility
statement to the owner.

**Speed and memory**, one calculation, 1 thread, best of three (`work/director/speed_trial31.out`):

| project | before the mesh | mesh alone | mesh with shared sites |
|---|---|---|---|
| `natural_neon` | 6.21 s, 1390 MB | 7.01 s, 1716 MB | 1.32 s, 340 MB |
| `thick_target_resonance` | 0.91 s, 388 MB | 1.30 s, 548 MB | 0.59 s, 272 MB |
| T2 | 0.26 s, 138 MB | 0.31 s, 150 MB | 0.26 s, 134 MB |
| 13N | 8.31 s, 225 MB | 8.32 s, 225 MB | 7.73 s, 122 MB |
| `stopping_yield` | 0.68 s, 333 MB | 0.66 s, 338 MB | 0.11 s, 105 MB |
| `layer_yield` | 0.30 s, 200 MB | 0.30 s, 205 MB | 0.10 s, 113 MB |
| `hybrid_potential` | 1.73 s, 338 MB | 1.73 s, 339 MB | 1.07 s, 106 MB |

Fits, mesh alone -> with shared sites: `stopping_yield` gradient 43 s -> 8 s, least squares 42 s -> 10 s; 13N
gradient 27 s -> 18 s; T2 gradient 36 s -> 34 s.  The mesh alone is slower than the engine before it (every point
of a data set held its own copy of the mesh's nodes in its range); the shared sites are what makes the mesh pay.

**Decisions of the director on what the expectations agent left.**
- `13N/expected` (814.349 pinned; engine before 814.351, now 814.346): left, within the check's tolerance; the pin
  already differed from the engine in other segments.
- `target_effect_ranges/expected`: regenerated by the branch that deletes the fold windows (it changes that project).
- `golden/T2/fit` (chi2 294.173, the fit of the engine before the mesh): left; `convert_params` needs its `param.sav`
  in the old names and now compares with `expected/chi2_renamed`.  To regenerate when that check has an input of its
  own.
- Historical records that quote old numbers (`tests/upgrade/bench/T2.tsv`, `format3_tree.md`, `MERGE_NOTES.md`) stay
  as records.

**Open.**
- **For the PI to confirm: the acceptance rule** of a changed folded number (`energy_grid.md` §7.1: per file, the
  largest distance to the refined reference at most 1.2 x the old engine's and the rms at most 1.03 x; converged
  files below 1e-7 pass; `hybrid_potential` named with 1.25).  It is the director's, not the PI's, until confirmed.
- `hybrid_potential` at its own settings is 5–6 % from its converged value in both engines: a finding about the
  project's settings for the owner.
- `shift_subpoint_angle` (b): tolerance 1e-5 (measured 3.5e-6) until the shift kernel puts a shifted segment on the
  mesh (branch 8), then 1e-6 again.
- A shared site's cotangent is one `PointTBar` call per reader (summing per site first changes the floating-point
  order: a rounding-class change, to do with the gradient's next branch); the serial `Unit::Calculate` path still
  solves all of a point's slots; `NodeTableCache` stays until shifted points are on fixed nodes.
- In work: the deletion of the four mechanisms of decision D3 (`feature/upgrade-5-d3-deletions`).

## 2026-11-12 (agent, Phase 5 branch 7 part 2, milestone 3) — fold windows, edge blending, the automatic tolerance and the σ(Ē) fallback deleted (`feature/upgrade-5-d3-deletions`)

**Done (decision D3).** Every point of a folded data set is folded with the whole kernel on the shared mesh; a user
sets for a folded data set the effect (`target`, `beam`, `detector`) and, optionally, the mesh's settings
(`quadrature`), nothing else.
- Engine: `src/TargetSetup.cpp` loses `AddTargetFallbackNodes`, `CheckTargetWindows`, `TargetWindowsFor` and the
  window and fallback branches of both loops of `SetUpTargetEffects` (968 -> 762 lines; the function 405 -> 365 lines,
  54 -> 46 `if`, 8 -> 4 `continue`).  `Unit` loses the blend weight, the fallback flags, depth and scope,
  `FoldThinTarget`, the probe (`DetachFromPlan`, the two-argument `Site::UsePrivateRecord`); `Unit::FoldSubPoints` is
  the one fold (it was `IntegrateTargetEffectForObservable`, behind a probe/tolerance front).  `TargetEffect` loses
  ranges, blend width, both tolerances, the automatic-window settings, `BlendWeight` (two), `GetRanges`,
  `GetTransitionWidth`, `GetAutoTolerance`, `IsAutoWindows`, `GetAutoWindowTolerance`, `GetAutoWindowMargin`,
  `HasWindows`, `MeanEnergyLoss`.  The hand-over and the stopping scale no longer track E_bar; the Coulomb reference no
  longer blends; `AZUREGrad.cpp` loses the two guards of `GradTargetEffectAdjoint` and the THM blend; the transform no
  longer refuses a blend (`WorkingCopy::CheckTransformedSegments`).
- Kernels: `KernelPoint` loses `blendWeight`, `fallbackOnly`, `fallbackNow`, `fallbackDepth`; `KernelTarget.cpp` loses
  `TargetFallback` and the blend.  A thin target inside one mesh interval is an ordinary folded point (the kernel's
  comment and `energy_grid.md` §2.2 say what it gives: the interpolated σ at the depth-averaged energy).
- Format 3: `fold_windows` (its form, the `Windows` value type, `FoldWindows`, `ParseWindows`, `SplitRange`,
  `E_WINDOW_SYNTAX`, the layer and Rutherford-ratio refusals of it) deleted.  A file that writes it is refused with the
  new `E_FIELD_DELETED` from a per-kind list `RecordKind::deleted` (exported in the schema JSON, printed in
  `format3.rst`).  The `quadrature` line never had a tolerance (`points`, `width_multiplier`, `per_width`): nothing
  of it is deleted; the `tolerance` of an extrapolation's adaptive grid is another thing and stays.
- Converter: an old `<targetInt>` line with ranges, `auto` windows or a tolerance converts with a note naming the
  line, what was dropped and what moves; a blend without windows (no effect in the old engine) is no note.
  `convert_spec_equivalence` requires one such note per windowed active line (shown failing without the note).
- Tests deleted: `tests/target_auto_windows` (automatic windows from the level scheme: placement, Note, threshold
  windows, refusals), `tests/target_window_fallback` (σ(Ē) outside the windows against a truth table, the
  |σ''|Δ²/24 error, blend in the kernel, refusals).  `tests/target_effect_ranges` is kept as a beam spread folded into
  A_y data (README rewritten), re-pinned in its own commit; `target_effect_shared` no longer strips the field;
  `layer_yield` expects `E_FIELD_DELETED`; `park_chi2_fallback` variant B is a convolved A_y on the adjoint.

**Moved numbers.** Only `tests/target_effect_ranges` (the one project of the tree, and of the old originals under
`tests/convert/fixtures/` and `examples/`, that set windows, a blend or a tolerance): χ² 3838.66 -> 3760.95; outside
its old window the largest A_y change is 0.046 (1.5 data σ), inside 0.004 (0.13 σ); 31 of 72 points keep their bits.
The base engine on the file without windows, and on the old file with the window widened to 0-1000 MeV, gives the
new outputs bit for bit.  This is larger than the old tolerance (5e-3 on χ², 0.002 per point) and than the data's σ:
the windows of a beam spread were a restriction of the physics, not an approximation.  Decision for the director.

**Verified how** (build `work/b-d3`, scripts `work/d3/`: `qset.sh`, `measure.sh`, `cmp_ter.py`).  Build 0 warnings;
consistency test; `ctest -R "experiment|convert|kernels|grad|upgrade_plan|upgrade_node|upgrade_mesh|consistency|schema"`
41/41 and `park_chi2|target_effect_line|parsers|schema` 5/5; `check.sh` of `target_effect_shared` and `layer_yield`
pass; T1, T2, T3, 13N, natural_neon, thick_target_resonance, layer_yield/models byte-identical to `tab-trial31`
(Brune, 109 files, run.log without timing lines), also under `AZURE_PLAN_CHECK=1 AZURE_NODE_CHECK=1` (outputs
identical, every check line "0 differ" / "0 off the tables" / "0 slots with a wrong writer").  The refusal
(`E_FIELD_DELETED`, one per field) and the conversion note shown on the old `target_effect_ranges` files; all 189 old
originals converted: one note, that file's.  Injected once each (then removed): the dataset kind's `deleted` list
emptied (`experiment_diagnostics` fails 3 checks) and the converter's note suppressed (`convert_constructs` 3,
`convert_spec_equivalence` 1).  `AZURE_GRAD_CHECK=exit` on T2 and `layer_yield/models/fit.azr` (1 thread): every
printed component equal to the base's (only the timing line differs).  Speed and memory, 1 thread, 7 alternated runs,
medians, base = commit ab54094a built the same way (`work/b-d3-base`): T2 0.25 / 0.25 s, 134 / 134 MB;
natural_neon 1.26 / 1.27 s (user 1.12 / 1.10 s), 340 / 340 MB; thick_target_resonance 0.57 / 0.57 s, 272 / 272 MB;
target_effect_ranges 0.59 -> 0.72 s, 104 -> 106 MB (the old file left 34 of 72 points unfolded; on the same windowless
file the two engines take 0.73 / 0.72 s, 109 / 109 MB).  Not run (the director's, on the last commit): the whole tree x
3 bases, `tests/run_tests.sh`, the gate, the fits, other thread counts.  `scripts/check.sh --quick` (JOBS=1,
OMP_NUM_THREADS=2) on commit e28c9812: 207 tests, 0 failed.

**Open.** The base and component loops of `SetUpTargetEffects` are now near copies (placement of slots: not in this
step's scope).  `park_chi2_fallback` no longer exercises the finite-difference fallback of the gradient (no p+p
construct the adjoint misses is left).  `.claude/skills/rmfit-campaign/references/policy-and-objective.md` describes
rmfit's own copy of the windows (another program; not edited).

## 2026-11-12 (agent, Phase 5 branch 7 part 2) — the stored expectations the shared mesh moves

**Done** (tests and docs only; no engine file changed).  For each of the six checks that `tests/run_tests.sh` failed on
024221f7, the new number was shown to be the mesh's and the closer one to today's engine refined eightfold
(`AZURE_REFINE_QUADRATURE=8`, `work/b-trial30`), then the expectation was regenerated by the project's own means:
- `thick_target_resonance/expected/chiSquared.out`: 30.9094 (15.3329 / 15.5765) -> 31.3524 (15.6762 / 15.6762); refined
  31.4508 (15.7233 / 15.7275).  Both segments read the same points and now share the mesh (step = the finer of the two
  multipliers), hence one value.  README re-pin paragraph with the yields.
- `convert_params`: (b) compares with `expected/chi2_renamed` (294.164, this engine at the old fit's parameters) instead
  of the old fit's 294.173 (`golden/T2/fit`, not regenerated: its `param.sav` in the old names is the check's input);
  refined 293.649; segments 5-8 (the folded ones) each closer to the refined value.
- `ec_relativistic_signature`: (b) accepts either refusal (count or signature): on the mesh the stale file holds 71240
  amplitudes where 71232 are needed, so the count refuses it first.  New (d): the stale file cut to 71232 with its own
  signature, refused by the signature alone (shown failing: with a matching signature it is read and chi2 is 697.463
  instead of 810.243).
- `shared_point_tables/expected/chi2_fit`: 109.250542525 -> 109.263831704.  At the old capped fit's `param.sav`:
  109.250994 (before), 109.252912 (mesh), 109.259623 (refined); the capped fit (5 calls) of the refined engine ends at
  109.219622, so the fitted value follows the path, not the rule.
- `shift_subpoint_angle` (b): tolerance 1e-6 -> 1e-5, measured 3.5e-6 (a0 13484.5150511 on the mesh, av 13484.5627097
  on its own grids, refined a0 13484.5302622).
- `upgrade/golden/T2`: `calc_start` and `calc_truth`, files `AZUREOut_aa=1_R=2.out` and `chiSquared.out` regenerated
  (`regen_compare.py --regenerate --only calc_start --only calc_truth`), ENGINE note appended; `normalizations.out` of
  `calc_start` kept (the regenerated one adds the format-3 follower row, no number moves).  Per-file rule passes
  (rms 5.49e-4 vs 5.57e-4, 5.17e-4 vs 5.27e-4); strict per-number rule does not on 2-4 numbers per file, e.g.
  `T2_true` segment 5 chi2 27.0364 vs refined 27.0713 (before 27.0864).

- Also re-pinned (test expectations whose pin was today's engine's bits, passing within 1e-3 but moved by the mesh):
  `natural_neon` 3342.51 -> 3342.7 (refined 3342.82), `hybrid_potential` 2.46614e6 -> 2.46411e6 (the named project of
  §7.1), `beam_profile_kernel` segments 7.0552 / 33.3841 -> 7.05514 / 33.3842; READMEs carry a re-pin line.

**Verified** (build `work/b-mesh3` of 024221f7, scratch `work/meshx`, 2 CPUs): `tests/run_tests.sh` 157 passed, 0
failed (before the three extra re-pins, which copy that run's own outputs); `repo_consistency_test.py`; `ctest -R
"upgrade|consistency"`; `pyazr_ec_relativistic`; the variants of convert_params, shared_point_tables, a0 and T2
(no Coulomb table) classified with `tools/classify_quadrature.py` against old / k = 8 runs: all files per-file
quadrature.  Not run: the whole tree in three bases (the director's `tab-mesh2`), `ctest -LE slow`.

**Open.**
- Branch 8 (shift kernel): restore 1e-6 in `shift_subpoint_angle` (b) once shifted segments are on the mesh.
- `golden/T2/fit` is from the engine before the mesh (its chi2 294.173 is no longer this engine's at its parameters);
  regenerating it rewrites `param.sav` in format-3 names, which `convert_params` cannot use as its old-name input.
- Not regenerated, within their tolerance: `13N/expected` (814.349 pinned; before 814.351, mesh 814.346; the pin already
  differs from today's engine in other segments) and `target_effect_ranges/expected` (3838.66 pinned; before 3837.11,
  mesh 3837.13).  Historical records quoting old numbers left as they are: `tests/upgrade/bench/T2.tsv` (10899.8),
  `docs/upgrade/format3_tree.md` l.108 (10899.8), `MERGE_NOTES.md` (109.250542525, 30.9094).

## 2026-11-12 (agent, Phase 5 branch 7 part 2, shared sites) — the points of a data set share the sites of the mesh (`feature/upgrade-5-shared-sites`)

**Done.**
- `include/SlotGroup.h`: the sharing group and its key (segment and term, c.m. angle bits, Legendre order of a
  coefficient row, target effect, stopping factor, frame shift, screening U_e and parameter), one block of sites per
  group in node order, each point a range of its block (ranges overlap).  Off the mesh a point is a group of its own.
- `WorkingCopy` (`src/WorkingCopyUnits.cpp`): `slotBlock_`, `slotWriter_`; `WritesSlot(point, i)` (a slot is written by
  its first reader in data order), `AddSlotGroup`, `LayOutUnits` moves whole blocks, `PlaceSubPointGrids` takes blocks,
  `SeparateSlotsByAngle` (a leaf point whose converted angles differ from its group's gets its own copies),
  `SetSlotRange` (tests), `CheckSlotGroups` (run by `AZURE_NODE_CHECK` / `AZURE_PLAN_CHECK` with the node check).
- Set-up: `AddMeshSlotGroups`, `GroupNodes`, `MeshGroupKey` (`src/MeshSlots.cpp`); in `SetUpTargetEffects` the loop
  of `AddGridSlots` over the queued mesh points (6 lines) is one call.  Hand-over: `ReplaceSubPointGrids` re-places a
  group at once on the union of its points' new ranges (`replaceGroup`); the site construction is one lambda for
  groups and own grids (`laySites`).
- Writers only on written slots: `Unit` set-up recursions (Legendre, channel rows, keys, Coulomb, EC, Recalc),
  `PrepareSubPoints`, `CalculationPlan::TaskUnits` (the energy pass), `ForEachResultUnit` (one record per site), the
  EC read-back of `intEC.dat` (the file keeps every slot of every point), the stopping scale, the shift step.  A slot
  whose solve throws fails every point that reads it (`EnergyPass.cpp`).
- Folds and gradient never write a slot's record: E1/E2 components, A_y, normalised coefficients
  (`src/PointPass.cpp`), the fd weights, A_y and coefficient adjoints and the THM fold weights (`src/AZUREGrad.cpp`)
  fold arrays with `FoldValues`.  A site's cotangent is the sum of its readers' cotangents: each `UnitCotangent`
  names its reader's plan task, and a T-matrix group sums them in reader order (`GroupCotangents`); the energy tasks
  are ranked by reads (`EnergyTask::reads`, `WorkingCopy::SiteReaders`) instead of units.  Both keep the order of the
  gradient's sums of one slot per point, so gradients and fits stay bit for bit (without them the T2 MIGRAD fit of
  200 calls moved in the last digits from call 81 on: its analytic gradient summed in another order).

**Audit of `Site`'s members** (depends on: N node energy, S segment and term, A angle, P the point; read at a slot by):
cm/lab/compound/nominal energy N (all formulas); `settings_`, `segment_key_` S (observables); `cm_angle_` A (Legendre,
Coulomb; key); `lab_angle_` A, a base slot holds the c.m. angle, a leaf's slot its point's converted lab angle (frames
of leaves; separated by angle); `legendre_order_` P, coefficient rows (GenMatrixFunc; key); `row_` P, never read on a
slot (`Row()` is -1); `geofactor_`, `sfactorconv_`, `edep_`, `eDependent*`, `screen_k_` N+S; `screen_ue_`,
`screen_param_` P's target effect (key); `stoppingPower_` N + effect + factor (fold gather; key); `targetThickness_`,
`targetLabEnergy_` P, never set on a slot (the fold reads the point's); `legendreP_`, Coulomb amplitudes A+N (+Q
coefficients of the effect, key); EC amplitudes, `ec_shifted_` N+S (energy task of the writer only); `thm_ps_` THM, off
the mesh; record (`slot_`, angular/natural counts) the site's, written by the energy pass only.

**Verified how** (build `work/b-sites`, scripts in `work/sites/`).  Build 0 warnings; consistency test; `ctest -R
"upgrade_plan|upgrade_history|upgrade_node|upgrade_mesh|kernels|grad|consistency"` 23/23; `check.sh` of 13N,
component_mapping, layer_yield, stopping_yield, rectangular_resolution, fit_modes, ratio_to_reference pass; T1, T2, T3,
13N, natural_neon, thick_target_resonance, layer_yield/models, stopping_yield/models, rectangular_resolution/models,
shift_sqrt_subpoints/models, 18O_p_a_thm, hybrid_potential, beam_profile_kernel (Brune): 193 files byte-identical to
`tab-mesh2` (run.log without timing lines).  Fits against `AZURE2.mesh2`, 1 thread, every output file identical: T2 and
13N `--use-gradient` (30 calls), T2 `--use-lm` (30), T2 MIGRAD and T2 `--use-gradient` (200 calls).  At 3 threads the
200-call fits give the same files as at 1 thread except the call count of `fitStatus.out` (116-118, the base varies
the same way: 117-118).  `upgrade_node_check` extended; the slot-group check shown failing with one point's range
shifted by one site (natural_neon: 17 slots with a wrong writer, 1 range off its mesh plan), then clean.
`upgrade_plan` ranks by reads (shown failing with the old cost by units).  `check.sh --quick` (206 tests) passed on
the work before the gradient order and again on commit 2ad5adfb (206 tests, 0 failed).

**Measured** (1 thread, 3 alternated runs, medians; `logs/meas-final.log`; base `AZURE2.mesh2`):

| project | slots (sites) base -> here | units of the pass | E-dependent values | evaluation | cross sections | wall | peak |
|---|---|---|---|---|---|---|---|
| thick_target_resonance | 74,634 -> 28,924 | 74,634 -> 28,924 | 558 -> 212 ms | 253 -> 133 ms | 163 -> 58 ms | 1.14 -> 0.52 s | 548 -> 272 MB |
| natural_neon | 373,449 -> 47,649 | 373,449 -> 50,813 | 2,179 -> 329 ms | 1,237 -> 353 ms | 913 -> 105 ms | 6.02 -> 1.15 s | 1.68 GB -> 340 MB |
| T2 | 32,032 -> 21,451 | 32,103 -> 21,522 | 56 -> 39 ms | 79 -> 64 ms | 10 -> 5 ms | 0.28 -> 0.23 s | 150 -> 134 MB |
| 13N | | 46,164 -> 3,825 | 120 -> 10 ms | 108 -> 18 ms | 58 -> 6 ms | 7.06 -> 6.60 s | 225 -> 122 MB |
| stopping_yield/models | | 202,888 -> 6,455 | 135 -> 4 ms | 122 -> 7 ms | 68 -> 1 ms | 0.61 -> 0.11 s | 338 -> 105 MB |
| layer_yield/models | | 91,877 -> 14,134 | 60 -> 10 ms | 52 -> 8 ms | 26 -> 2 ms | 0.28 -> 0.10 s | 205 -> 113 MB |
| 18O_p_a_thm (off the mesh) | | 14,688 -> 14,688 | 13 -> 13 ms | 13 -> 13 ms | | 0.10 -> 0.10 s | 116 -> 116 MB |
| T1, T3 (unfolded) | | unchanged | unchanged | unchanged | | 2.63 / 0.15 s both | +0.2 MB |

T2 fits of 200 calls: MIGRAD 31.5 -> 29.9 s (-5 %), `--use-gradient` 32.0 -> 29.9 s (-6 %).  Not shared: T2's lab-frame
differential data (a c.m. angle per point), shifted segments and THM (grids of their own).

**Open.** `NodeTableCache.h` kept (re-placed blocks adopt rows other copies computed).  Lab-frame differential data
remain groups of one point (each point its own c.m. angle) until the angle per node (decision D5).

## 2026-11-12 (agent, Phase 5 branch 7 part 2, after milestone 2) — the stops of the shared mesh resolved, a sharper rule

**The moved unfolded numbers were the comparison, not the engine.** Mechanism, by experiment:
1. d82f3518 run with the Coulomb table gives `param.par`, `parameters.out`, `scatteringLengths.out` of `7Li_p_a`,
   `6Li_d`, `target_effect_ranges` and `ratio_to_reference` `fit_ratio`/`fit_sigma` byte-identical to `tab-trial30`
   (also the director's `tab-mesh2`).
2. The moved files of the earlier quick run (tag `q1`, written at 02:41, its `run.log`s 70 s later by another run) are
   bit for bit the outputs of the same engine with `AZURE2_COULOMB_TABLE` unset: the engine then opens
   `~/.azure2/coulomb_table.bin`, refuses it (format 1) and uses the exact Coulomb functions (charged channels move at
   1e-7).
3. The 53 "unfolded" fallback numbers of `target_effect_ranges` `aa=5` and the 24 of `ratio_to_reference` come from the
   same no-table run.
4. `natural_neon` `aa=2` line 723 is a folded number that k = 4 leaves at today's value to the last digit; against
   `tab-fine8` it is folded and conform.
5. Guard: `tools/run_one.sh` and `tools/runall.sh` refuse to run without `AZURE2_COULOMB_TABLE` naming a file (CLAUDE.md
   recipe says so).

**`hybrid_potential`.** `quadrature points` has no effect on either engine (10, 20, 40, 80 give the same bits): every
interval is set by the resonance tails (step 0.05 |E − E_r| at `per_width 50`), never by the smooth step, so the node
phase of the smooth lattice is not the cause.  In k (`AZURE_REFINE_QUADRATURE`, both engines, distance to today's engine
at k = 64, rms over the 76 folded values): today 5.14e-2, 2.00e-2, 5.14e-3, 1.35e-3, 3.43e-4 at k = 1, 2, 4, 8, 16; mesh
6.17e-2, 1.71e-2, 5.13e-3, 1.35e-3, 3.52e-4.  Both converge at second order to the same value; the mesh is 1.20 x
further at k = 1, 0.85 x at k = 2, equal from k = 4. The error is one sign and nearly one size at every point (rms =
max): a part of the integral common to all thick-target points.

**Rule.** `tools/classify_quadrature.py`: per file also rms of the relative distances <= today's (no factor), both
measures over every folded number of the file; the fine reference per file from a list (`tab-fine8,tab-fine4`), the tag
printed.  Test `tests/consistency/classify_quadrature_test.py` (ctest `tools_classify_quadrature`).  `energy_grid.md`
section 7.1.

**Verified** (build `work/b-mesh3`, scripts `work/mesh3/`): build 0 warnings; consistency test; quick set plus the five
projects, `hybrid_potential`, `17O` (`qset.sh m3-q1`): 138 unfolded files byte-identical to `tab-trial30`, 197 files:
163 identical, 29 quadrature, 5 failure (13N twice, `layer_yield` `aa=3`, `hybrid_potential` twice).  Injected error
(every other node dropped in the lower half of each range, `work/mesh3/inj/AZURE2`): 23 failures; the rms catches 13N,
`beam_profile_kernel`, `ratio_to_reference`, `T2_true` that the maximum passes; `T2.azr` and `thick_target_resonance`
pass both.  1 against 2 threads on T2 and 7Li_p_a identical.  Not run: the whole tree (the director's `tab-mesh2`), 3
threads (budget 2).

**Open.** The rms rule without a factor fails 13N by 0.2 % (rms) and 0.6 % (chi2), `layer_yield` `aa=3` (a fold exact to
3e-8: 1.18 x, 1.44 x) and `hybrid_potential` (1.21 x): a decision of the director.  `layer_yield` `aa=3` (a layer, 1/E
stopping, a nearly constant sigma): today's 10 intervals fit the point's energy loss exactly with both ends as nodes;
the
mesh's nodes are the nucleus's ladder (step <= today's) with the two edge intervals clipped at the support, whose
interpolation error does not cancel as a full trapezoid's does (from the construction, not measured apart).  No engine
file changed on this step.

## 2026-11-12 (agent, Phase 5 branch 7 part 2, milestone 2) — the engine on the shared mesh (`feature/upgrade-5-shared-mesh`)

**Done.**
- One placing function per nucleus: `WorkingCopy::SetUpTargetEffects` queues every folded point of an unshifted, non-THM
  segment with an adaptive grid (`src/MeshSlots.cpp`, `QueueMeshPoint`: support in E_x, settings, a hard lower edge at
  a floor), plans each nucleus once (`PlanMesh`, anchors `MeshAnchorsOf` from the input widths, quantized in E_x, every
  entrance pair of the nucleus's folded points) and gives each point the nodes of its range (`MeshSlotEnergies`).  Both
  `GenerateGrid` branches of the set-up go through it; `AddGridSlots` is the one loop that makes slots.
- The quadrature is `GridQuadrature` on the mesh nodes, clipped at the support: `KernelNodes::supportLow/High`
  (`Low()`, `High()`), set by `PointFold.cpp` for a mesh point; the Gaussian, target and beam-profile kernels read the
  support from there instead of the end nodes.
- Planner: `MeshSupport::lowHard/highHard` (a hard edge is a node and the range ends at it; test `HardEdges`, shown
  failing with the floor node removed: 4 failures); `finest` renamed `otherCoarsest`.
- Hand-over: `ReplaceSubPointGrids` re-plans a nucleus when its quantized anchors change (`WorkingCopy::meshAnchors_`)
  and re-places only the points meeting `MeshChangedRegion`; `AZURE_PLAN_CHECK=1` throws if a point outside would move.
- `tools/classify_quadrature.py`: the per-file rule (1.2 x today's largest relative distance to the fine reference),
  printed beside the strict per-number rule; `energy_grid.md` section 7.1 holds the rule as the director's decision.

**Verified how.** Quick set (`work/mesh2/qset.sh`) and the folded projects of the tree, Brune basis, against
`tab-trial30` / `tab-fine4`; fits against b-trial30; gradient checks; 1-thread measurements (see the branch report).

**Open.** Unfolded numbers moved in 6Li_d, 7Li_p_a, target_effect_ranges, ratio_to_reference fit_ratio/fit_sigma
(param.par, parameters.out, scatteringLengths.out at 1e-7): the moved values equal the base engine's without the
Coulomb table; mechanism not found.  hybrid_potential exceeds the per-file factor (1.22 x).  A point's slots are all
mesh nodes in its range: units of the pass grow 30-55 % (T2, natural_neon, thick_target_resonance) while distinct
energies fall; calculation time up to +43 %.  `beam_profile_kernel` has no fine reference (its tab-fine4 run stopped at
900 s).  NodeTableCache kept (shifted own grids).  Expectations of folded projects not regenerated.

## 2026-11-12 (director, 03:00) — merged: the refusal switch, both readers strides, the extension drill, the mesh planner; the yardstick of the shared mesh

**Merged tonight, each with the tree of a rehearsal that passed the whole verification** (all projects × 3 bases,
`ctest -LE slow`, slow ctest, `tests/run_tests.sh`, spec equivalence, the converted-tree gate, 1 against 3 threads,
eight fits against the build before): the refusal switch; the readers' first stride (every slot on a node table,
energy tasks grouped by node); the extension drill; the mesh planner; the readers' last stride (the sub-point units
are gone).  Numbers of the last rehearsal (`rehearsal/trial30`, logs `work/director/v4_trial30.*`, `fits30.log`): all
projects 823 / 796 / 823 files identical to the drill rehearsal; ctest 206/206, slow 10/10; run_tests 157 passed, 0
failed; spec equivalence 184 of 184 converted equal; gate 483 (456 formal) identical, 38 classified, 0 differ, 0
missing; eight fits identical.  The new reference of the whole-tree comparison is `tab-trial30`.

**Review of the readers' last stride** (agent's entry below).  The design is what the brief asked for in one cut
instead of two: `Site` (what depends on a point at one energy) and `Unit` (a data point: a `Site` plus row, fold,
task, record), the nodes of a folded point as `Site`s in one vector of the working copy.  The +3 % on the T2 MIGRAD
fit accepted with the first stride is gone (agent: 36.2 s against 36.4 s before the node tables, 1 thread, 6
alternated runs; my side-by-side fits: 37 s against 38 s for the build before).  One merge conflict was mine to
resolve: the drill's `SFactorConversionValue(…, shape_perMeV)` and reference-table leaf were written on `Unit`, the
readers' on `Site`; both are on `Site` now (`include/CoulombReference.h`, `src/CoulombReference.cpp`,
`src/ReferenceTableLeaf.cpp`).  Not done by the agent and still open: `include/NodeTableCache.h` with its mutexes
(needed while re-placed nodes share rows across copies; it goes with the shared mesh), the names that still say
sub-point.

**Review of the extension drill** (agent's entry below).  No defect found in the three features; each `check.sh`
fails with its injected error and passes.  What I add to the agent's open list: (1) the definitions of `rsl`, `2l2`,
`4pi` follow `docs/upgrade/data_kinds.md`; the LEXFOR page could not be fetched tonight, so they are **not checked
against the source**; (2) a fold or a shift that leaves the reference table gives a not-a-number point instead of a
refusal at set-up (the shared mesh knows the supports at set-up and can refuse); (3) exp(gE) is written twice
(`SFactorConversionValue` and `SigmaAndS::Conversion` in `src/DataOutputFiles.cpp`); (4) pyazr's `calculate_sfactor`
and the API's `conv` know only the plain S factor.  The cost table and the verdict are in `FOR_THE_OWNER.md` ("How to
add a feature now").

**Review of the mesh planner** (agent's entry below; header, implementation and `energy_grid.md` §2.8 read).
Accepted.  Director's decisions for the engine step, written into its brief (`work/briefs/BA_mesh_engine.md`): a hard
edge at a physical floor (a threshold, E = 0) is an exact node, flagged per edge of a support; anchors are quantized
in E_x; the power-of-two ladder from E_x = 0 stays (its cost is 1.15–1.80 × the smooth nodes, and smooth nodes are
the cheap ones); the width of `MeshChangedRegion` stays; the local variable `finest` of `PlanMesh` is the coarsest
step of the other lattice on the piece and is to be renamed.  One thing the measurement shows and a reader should
know: point by point the mesh interval can be longer than today's own interval on 1.7–3 % of a support's length
(where today's grid is finer than its own declared step); property 1 is stated against the declared step.

**The yardstick of the shared mesh (M0 of `work/briefs/AX_shared_mesh.md`), measured by the director on the
unchanged engine.**  `tab-fine4` = every project with `AZURE_REFINE_QUADRATURE=4` (build of the drill rehearsal);
largest relative distance of today's default quadrature to it, Brune basis, per project (`work/cls-tab-fine4.tsv`):

| project | largest relative distance k = 1 → k = 4 |
|---|---|
| hybrid_potential (a data set with `quadrature points 10`) | 4.4e-2 (a uniform offset) |
| ratio_to_reference, folded | 1.1e-2 (χ² 3.9e-2) |
| T2 | 5.1e-3 (χ² 1.1e-2) |
| T2_true | 4.6e-3 |
| stopping_yield, fit | 3.9e-3 |
| thick_target_resonance | 3.1e-3 (χ² 2.4e-2) |
| layer_yield, fit | 2.0e-3 |
| shift_sqrt, no shift | 1.2e-3 |
| layer_yield | 8e-4 |
| natural_neon | 7.4e-4 |
| 18O_p_a_thm | 4e-4 |
| target_effect_ranges | 3.5e-4 |
| stopping_yield | 2.4e-4 |
| rectangular_resolution, fit | 1.3e-4 |
| 13N | 1.3e-4 |
| 17O | 1e-4 |
| 7Li_p_a | 6.7e-5 |

These are the quadrature errors a user has **today** with the default settings; the shared mesh is accepted where it
is at least as close to the fine reference.  `tab-fine8` (k = 8) is being run to see which projects still move
between 4 and 8; the column is added to `energy_grid.md` §7 with the shared-mesh merge.  The acceptance rule the
shared-mesh agent was given beside the strict per-number rule of `tools/classify_quadrature.py` — per file, the
largest relative distance to the fine reference not above 1.2 × today's largest for that file, and no unfolded
number moved — is a director's decision **to be shown to the PI**.

**In work.**  The engine on the shared mesh (`feature/upgrade-5-shared-mesh`, M2 of the brief): one placing function
from the planner, quadrature on the mesh nodes, the hand-over by decision D1.  After it: the deletions of D3
(automatic tolerance, fold windows, edge blending, the σ(Ē) fallback; the readers' entry below lists the places),
then classification, measurement and documents.

**Small things still open** (none blocks): check 34's five pending files still name `<targetInt>` /
`<segmentsData>` in comments and in error texts (`include/ParameterKind.h`, `src/TargetSetup.cpp`,
`src/DataParameters.cpp`, `include/Unit.h`, `experiment/src/Schema.cpp`; the error texts are user-visible, so they
change with a check of the expectations, after the shared mesh has stopped editing `TargetSetup.cpp`); the "calls"
line of `fitStatus.out` varies by a few between threaded runs (Minuit2's own counter; Minuit2 is not patched): to be
said in the statement on fit reproducibility; `natural_target/models/bad_level` does not set up in `node_count`.

## 2026-11-12 (agent, Phase 5 branch 7 part 3) — the sub-point units are gone: a slot of a folded point is a `Site` (`feature/upgrade-5-grid-readers2`)

**Done.**
- `include/Site.h`: one calculation site, what depends on a point at one energy (c.m./lab/compound energies, angle and
  Legendre rows, channel rows, Coulomb amplitude, screening, stopping, EC amplitudes, its result record).  The formulas
  (channel rows, Legendre rows, Coulomb amplitude, external capture) are `Site` functions writing that site only; the
  `Unit` functions of the same names run them for the point and then each slot (the recursion that lived inside each
  formula; `OwnSite()` for the point alone).  `Unit` (a data point) is a `Site` plus row, fold, plan task, model
  record, fold settings.  `GenMatrixFunc`, `AMatrixFunc` (`PointTBar` too), `RMatrixFunc`, `THMMatrixFunc`, the
  amplitude matrix, `SolveAtOneEnergy`, the Coulomb reference's callback (`at(const Site *)`) read a `Site`.
- The slots are `Site`s in one vector of the working copy (`WorkingCopy::SlotAt`, `Slots`; a point's range
  `SubPointFirst`/`NumSubPoints`), laid out with the points (`LayOutUnits`).  One site index for plan tasks, energy
  tasks, node slots and cotangents: `SiteAt(s)`, the points then the slots at `NumUnits() + k` (`SlotSite`,
  `SiteIndex`): the same numbers as the old one unit vector, so no task, node or record moved.  `SubPointAccess.h`
  deleted; `UnitIndex` deleted.
- Overflow table placed without its hash map (a slot shares the node placed just before at the same key, else a new
  node; the energy tasks group by key).  `AZURE_NODE_CHECK` stays as the invariant check (every slot on a node of a
  table, stage nodes in key order, the node's key equal to the site's key with E_x computed from its c.m. energy bit
  for bit, rows the site's) and now also runs under `AZURE_PLAN_CHECK`.

**Verified** (build `work/b-readers2`, scripts `work/readers2/qset.sh`, `bq.sh`, `speed.sh`; at every commit): build
0 warnings; consistency test; `ctest -R "upgrade_plan|upgrade_history|upgrade_node|kernels|grad|consistency"` 22/22;
`check.sh` of 13N, component_mapping, layer_yield, stopping_yield, rectangular_resolution, fit_modes all pass; T1, T2,
T3, 13N, natural_neon, shift_sqrt_subpoints/models, 18O_p_a_thm, thick_target_resonance, layer_yield/models (Brune)
159 files identical to `tab-trial26` (run.log without `%` lines); T2 and 13N `--use-gradient` (30 calls, 1 thread)
identical to `b-trial27`.  `node_check.sh` shown failing with a slot on the previous overflow node (sqrt: 47165 on a
wrong node; T2 fit: 6175), then clean.  Not run: the whole tree x 3 bases, full ctest, `run_tests.sh`, the gate.
Speed, 1 thread, 6 alternated runs each, medians (machine shared, load about 9; logs `work/readers2/logs/speed1.txt`,
`speed2.txt`): T2 MIGRAD 200 calls, final head 36.24 s against `b-trial25` 36.38 s (-0.4 %) and `b-trial27` 36.59 s
(-1.0 %); after the slot commit (overflow still hashed) 35.80 / 35.91 / 36.55 s.  natural_neon calculation 6.05 s
against 6.23 / 6.21 s.  Peak RSS: natural_neon calculation 1,424,604 kB against 1,460,656 kB (`b-trial27`, -2.5 %);
T2 fit 326 MB against 345 MB.

**Open.** `NodeTableCache.h` (re-placed slots share rows across copies through it; it goes when the shared mesh stops
nodes from moving with the levels) and the slots' copy-on-write rows are not deleted.  The names still say
"sub-point" (`GetSubPoint`, `NumSubPoints`, `AddSubPoint`, ...): a rename to "slot" is mechanical.  For the shared
mesh: slots are placed in four places (`WorkingCopy::SetUpTargetEffects`, two branches, `src/TargetSetup.cpp` around
the two `GenerateGrid` calls; `AddTargetFallbackNodes`, same file, the sigma(E_bar) fallback; `WorkingCopy::
ReplaceSubPointGrids`, same file) and weighted in one (`FoldWeights` / `FoldValues` / `FoldPoint`, `src/PointFold.cpp`,
through the kernels).  D3 lines deletable with a shared mesh: the automatic tolerance and the bare probe
(`Unit::FoldSubPoints`, `src/PointPass.cpp`), the edge blend (`TargetBlendWeight`: `TargetSetup.cpp`, `PointPass.cpp`,
`CoulombReference.cpp`, `AZUREGrad.cpp`), the sigma(E_bar) fallback (`AddTargetFallbackNodes`, `TargetFallbackScope`,
`Unit::FoldThinTarget`), the fold windows of the set-up (`SetUpTargetEffects`).  A shifted folded point still moves
its slots rigidly (`Segment::UpdatePointEnergyWithShift`) onto the overflow table until the shift-kernel branch.

## 2026-11-12 (agent, Phase 5 branch 7 part 2, first milestone) — the planner of the shared mesh, alone (`feature/upgrade-5-mesh-planner`)

**Done.** `include/MeshPlanner.h`, `src/MeshPlanner.cpp`: `PlanMesh(anchors, supports, tails)`, a pure function (no
engine object) that plans one mesh per compound nucleus from quantized anchors in E_x, one support per folded point with
its data set's settings, and the tail constants; `MeshChangedRegion` gives the stretches whose nodes can move when
anchors change (for the hand-over of decision D1). Rules (header comment, `energy_grid.md` §2.8): smooth nodes n·L on a
nested power-of-two ladder from E_x = 0 (L = today's smooth step quantized down); core and tail nodes by today's
formulas
from each anchor; an interval of a lattice kept unless one other lattice, active over all of it, is finer over all of it
(than half its length between a ladder and an anchor's lattice, as today's grid keeps both in that band). Not wired into
the engine; no engine file changed.
Tools: `tests/upgrade/node_count --plan [--detail]` (the planner on every folded project's supports and anchors as the
set-up has them), `tests/upgrade/mesh_fold_check [--with-edges]` (a Breit–Wigner folded on today's grids and on the
mesh against a converged reference).

**Verified** (build `work/b-mesh`, logs `work/mesh/`): build 0 warnings; `upgrade_mesh_planner` 700,675 checks, 0
failures (never coarser than today's step against `AdaptiveIntegrationGrid` on named and 200 random cases, today's
lattice nodes kept bit for bit, shared, order-independent, local), shown failing with the ladder one notch too coarse;
`ctest -R "upgrade_plan|upgrade_node|consistency|mesh"` passed; T1, T2, 13N in three bases: every output file identical
to `tab-trial27` (run.log differs only in the progress bars). Measured: no folded project's mesh above 1.04 × today's
distinct energies (T2 0.54, `natural_neon` 0.37); the fold error on the mesh not larger than today's in 34 of 36 cases,
the two others (+0.3 %, +1 % of an error of 3 × 10⁻⁵) from clipping the edge intervals (equal with the edges as nodes).

**Open.** The support's edges as nodes where a support ends at a physical floor (thick target at 1 keV, THM below
threshold); anchors quantized in E_x instead of each pair's c.m. (numbers move by the quantization); the ladder's cost
(1.15–1.80 × the smooth nodes); `MeshChangedRegion` is wide when a nucleus has a coarse smooth step. Milestone 2 (the
engine on the mesh) waits for the node-storage branch.

## 2026-11-12 (agent, extension drill, `feature/upgrade-5-extension-drill`) — three new kinds of data by the recipes, and what each cost

Done (one commit per lane, then the recipe corrections):
- `quantity modified_sfactor` (S* = E sigma exp(2 pi eta + g E), field `shape <g> 1/MeV`, new dimension
  `inverse_energy`): the S factor's `Scaled` tree with g carried by the S-factor `FactorRef`; physics
  `modified_sfactor_conversion` with its unit test and binding; the AZUREOut columns hold sigma and S*.
  `tests/modified_sfactor`.
- `quantity legendre` with `expansion rsl | 2l2 | 4pi` (data and extrapolations) and `absolute` for extrapolations:
  one table `experiment/LegendreExpansion.h`; on data a `FactorRef::Kind::Legendre` of the point's order; on
  extrapolations `Segment::LegendreCoefficients` (the `.acoeff` writer and the API). `tests/legendre_expansions`.
- `quantity ratio_to_reference` (sigma / a tabulated sigma_ref, field `reference <file> [unit <u>]`): tree kind
  `ReferenceRatio`, `LeafPart::Reference`, `src/ReferenceTableLeaf.cpp`; the Coulomb leaf's fold-and-blend factored
  into `FoldReferenceAtPoint`; `E_REFERENCE_TABLE`. Folded points work through the existing seam.
  `tests/ratio_to_reference`.
- Recipes of `CLAUDE.md` corrected where the drill found them silent or wrong; consistency check 35
  (`IsReferenceLeaf`, not `IsCoulombReference`, where the model pass skips or drops a reference leaf).

Verified how: per lane, build 0 warnings, consistency, the experiment/convert/predictions ctest subset, the lane's
`check.sh` (each shown failing with an injected error), a handful of projects against `tab-trial25c`; at the end the
full comparison, `ctest -LE slow`, `tests/run_tests.sh`, the converted-tree gate (the report gives the numbers).

Open: a fold's sub-point or a shifted energy outside a reference table gives a NaN point (not refused at set-up: the
sub-point grid is made in files another branch is rebuilding); differential data relative to a reference are not
built; pyazr's `calculate_sfactor` and the API's per-point `conv` know only the plain S factor (a modified_sfactor
segment's model is already S*); fitStatus.out's MIGRAD call count varies between 2-thread runs of one fit (also on
the base engine; defect outside the drill).

## 2026-11-12 (agent, Phase 5 branch 7 part 2, stopped early) — every slot on a node table; energy tasks by node (`feature/upgrade-5-grid-readers`)

**Done.**
- Every slot is on a table (`include/NodeTable.h`, `src/NodeTable.cpp`): per plan task one slot per unit it calculates
  (its point, or its sub-points), naming a node of the stage table or of the working copy's overflow table
  (`NodeOverflow`). The hand-over places the overflow table from nothing for the tasks whose units moved (re-placed
  grids: `ReplaceSubPointGrids` now returns their tasks; shifted points) and the tasks already on it
  (`WorkingCopy::PlaceOverflowNodes`, end of `FillEnergyShiftsFromParamsImpl`); its rows are released before units are
  re-energised (`ReleaseOverflowRows`), so they do not trigger the units' copy-on-write. The copy's node state is one
  member (`NodeState`); before, a copy did not copy the node table at all.
- `AZURE_NODE_CHECK=1` checks every unit of every task at every hand-over (also inside Minuit's parallel gradient,
  where there is no point pass) and every point pass: 0 off the tables, 0 on a wrong node, 0 rows differ. Test
  `upgrade_node_check` (`tests/upgrade/node_check.sh`), shown failing with an injected slot on the next node.
- The seam for the shift-kernel branch: the comment of `NodeOverflow` lists what branch 8 deletes.
- First reader moved: `CalculationPlan::EnergyTasks` takes nucleus, route and E_x from the slots' nodes
  (`EnergyTask::node`); its unused compound argument is dropped (also `WorkingCopy::EnergyTasksOf`).

**Verified** (build `work/b-readers`, logs `work/readers/logs/`): all projects x 3 bases against `tab-trial25c` (tag
`tab-readers-full`, after the energy-task commit): 786 / 759 / 786 files, all identical; the Brune tree with
`AZURE_NODE_CHECK=1` (`tab-readers-chk`): 186 checks, 3,255,304 units checked, 653,260 of them on overflow tables, 0
off, 0 on a wrong node, 0 rows differ, the 593 output files other than `run.log` identical; `check.sh --quick` at the
first commit: 204/204; `tests/run_tests.sh` 153 passed, 0 failed; converted-tree gate 483 (456) identical, 38
classified, 0 differ, 0 missing; fits (T2 `--use-gradient`, `--use-lm`, MIGRAD 200 calls, 13N `--use-gradient`,
`layer_yield` and `rectangular_resolution` fits) identical at 1 and 3 threads and to `b-trial25` at 1 thread.
Speed (1 thread, 6 alternated runs with `b-trial25`, median): calculations T2, thick target, THM, `natural_neon`
within -4 to -1 %, 13N +0.1 % (repeat; +3.8 % in a run beside `run_tests.sh`); 13N `--use-gradient` fit -3 %;
`natural_neon` against `b-trial23` (before the node table) +0.6 %; **T2 MIGRAD 200 calls +2.5 to +3.6 %** (three
series): the overflow table placed from nothing in every hand-over that moves a unit (T2 shifts a folded segment at
every evaluation). Peak RSS unchanged within 1.5 %.

**Open (readers not moved, objects not deleted).** The T2 fit's +3 %: place only the moved tasks' overflow nodes
instead of all off-stage tasks (needs a way to drop nodes no slot names).  The observables the energy pass writes
(unit records), the
kernels' gather (`PointFold.cpp`), the result records of sub-points, the gradient's unit cotangents, the
external-capture amplitudes and `intEC.dat`, THM `thm_ps_`; `SubPointAccess.h`, `NodeTableCache.h`, the copy-on-write
tables, sub-point `Unit`s. What a slot must carry for them, found here: the c.m. energy of the entrance pair (not
recoverable bit for bit from E_x: E_x is computed from it), the stopping (moved per copy by a fitted stopping scale:
the per-copy slot data must be copy-on-write), the Legendre rows and Coulomb amplitude (angle and energy), the
screening k, the EC amplitudes (last bits per point), and a record per slot. The `GenMatrixFunc` / `THMMatrixFunc`
observables read all of these through `Unit` getters, so moving the energy pass means giving them a slot view.

## 2026-11-11 (agent, refusal switch) — the engine reads format 3 only; old files refused, the old reading path deleted (`feature/upgrade-4-refusal-switch`)

**Done** (`AU_refusal_switch.md`):
- **The refusal** (`project_format3.md` §8.4, built): `Config::ReadConfigFile` refuses a text that is not a format-3
  project before reading anything of it, with one message, `io::NotFormat3Message` (moved from `convert/` to
  `io/ProjectFormat.h` so that the engine library links nothing of the converter): the format and the first old
  section with its line, `AZURE2 --convert OLD OLD_f3.azr` and `pyazr.convert_azr(...)`, "the converter keeps OLD
  unchanged". Own lines for a `format=3` header over an old section and for a text without `<config>` (empty, not a
  project). The CLI exits 3 with the message as the run's last line; `--validate`, `--describe`,
  `--print-predictions` end with it and exit 3; the API session raises it; pyazr's `require_format3` uses it.
  `--convert OLD NEW` writes format 3 without `--to` (`--to 3` accepted; the format-2 output and
  `convert_azr(to=2)` are gone; `ConvertToFormat2` is the converter's first step, `convert/src/Format2.cpp`).
- **Parameter files**: `AZURE2 --convert-params OLD.sav NEW.sav --project OLD.azr`, `pyazr.convert_params`
  (`convert::ConvertParameterFile`): the map of the OLD project's conversion (`--project` is the old project: a
  format-3 project does not determine the map, its records may have been renamed); `tests/convert_params`.
- **The old path deleted** (map in `format3_tree.md` §5, made with a coverage build): the old producers of specs
  (`ProjectSetup.cpp`), Config's format-1 conversion and `<potential>` reader, the `<scatteringLengths>` reader and
  its chi2 line, the `<thm>` member lists (`ReadFormat3ThmSection` takes members and weight tables from the records;
  `Format3ThmLines`, `ParseThmBlockLines`, `CheckThmBlockLines` gone), the bindings without a Python caller.
  `Config::format3` is `Config::project`. `AZURE2Core` links no `convert/` target (0 undefined `azure2::convert`
  symbols in `libAZURE2Core.a`, 5 before; link lines: `src/CMakeLists.txt`).
- **Small defects** (item 6): density warning names `<target>.areal_density` and `grid_max`; a refused THM edit
  names its experiment once; `Project.set/rename/remove` of an unknown record raise `FormatError`; a `<config>`
  directory gets its trailing '/', `write_output_files` creates the output directory; `sav_fixed_params`,
  `inactive_lines` take `AZURE2_BIN` (their CLI cross-checks now run); `W_PARAMETER_SETTING_UNKNOWN` (a
  `<parameterSettings>` row naming no parameter; a warning: the GUI leaves such rows behind).
- **Tests of old files**: `angular_coefficients` (g) `bad_norm.azr` deleted (refusal on the record,
  `E_MODE_NOT_ALLOWED`), `bad_test.azr` to the converter's fixtures (the converter's refusal);
  `physics_pathway_table` and `upgrade_capture_ay_requests` leave `tests/convert/fixtures/` out (they ran the old
  originals through the engine); `thm_optical` without `segments=`; consistency check 34 (no old section named in
  the engine; five comment-only files of the two parallel branches listed as pending).
- `convert_spec_equivalence` kept unchanged as a converter test: its old side is the converter's old-line reader
  (`convert/SpecsFromOldLines`, which the converter needs for prior-centre rows), no longer an engine producer; the
  guarantee that the engine's results did not move is the converted-tree gate against `tab-trial8`.

**Verified how**: per commit, build 0 warnings, the consistency test, `ctest -R
"convert|cli|consistency|project|validate|...|thm|model_data"` (65), the `check.sh` of `cli_route_exit`,
`convert_params`, `13N`, `thm_fit`, `natural_target`, `angular_coefficients`, `target_density`, `nuisance_kinds`,
`thm_options`, `thm_experiment`, `model_data`, `scattering_length`, `18O_p_a_thm`; T1-T3, 13N, 17O, 18O_p_a_thm,
scattering_length menu 1 byte-identical to `tab-trial24c`. Every build 0 warnings. At the end (engine of d93192f6; the
later commits change tests and documents only):
all projects x 3 bases against `tab-trial24c` (tag `tab-refuse1`): 780/753/780 files, 778/751/778 identical, the
two others per basis the `run.log` of the removed `angular_coefficients/models/bad_norm.azr`, `bad_test.azr`
(only in the reference); `ctest -LE slow` 202/203 (`physics_pathway_table` ran the fixtures through the engine:
fixed in 2756d7ee, it and `upgrade_capture_ay_requests` rerun: pass); `ctest -L slow` 10/10;
`tests/run_tests.sh` 152 passed, 1 failed (`format3_engine`: its old-file run, rewritten in 8483d5b9 and rerun:
pass); converted-tree gate (`grefuse1`) 483 (456) identical, 38 classified, 0 differ, 0 missing; T1-T3 and 13N at 1
against 3 threads, 3 bases: 12 identical; fits at 1 thread, base build (`b-trial24`) against this one, MIGRAD 200
calls: T2, 13N, 13N `--use-gradient`, 13N `--use-lm`: every output file identical; the T1 and T3 generators
(through the converter) with the base build and with this one: every file identical. The refusal for each case
(CLI run, `--validate`, `--describe`, `--print-predictions`, pyazr `azure2`, `Project.open`, the engine session):
the message and exit 3 (`cli_route_exit`, `convert_azr_test`). `convert_params` fails with a wrong name injected
in the map (chi2 1.36e43 against 294.173). Edit audit (`tools/edit_audit.sh`, once, alone, 3 threads, 32 min): 129
checks, 1,682 edit places, 36 hit nothing
(all by design, the audit's 36), 1,506 runs, no vacuous edit (`format3_tree.md` §7); `stopping_yield` fails under the
audit only (the tool's wrapper hides `azure2-stopping`).

**Open**: the five files of the parallel branches still name old sections in comments (`include/Unit.h`,
`include/ParameterKind.h`, `src/TargetSetup.cpp`, `src/DataParameters.cpp` (also its `EffectParameterName`
fallback to the old `density_<n>` names, now unreachable), `experiment/src/Schema.cpp`): check 34 lists them as
pending; after those merges, remove the words and the list. `TargetEffect`'s refusals of the old block texts
(unreachable, kept: the spec carries text members both producers wrote). `tools/convert_tree.sh`'s old-engine pass
needs a pre-switch binary (the gate does not use it). `tests/pyazr/model_data_test.py` edits `<segmentsData>` lines
of a format-3 text (an edit that hits nothing: its gradient case runs on all data sets instead of segment 5; the
check still holds, the data part cancels). User-guide pages of `format3_tree.md` §6 other than THM still show old
lines. Sphinx not installed here: the docs were not built.

**Line counts** (acb0cc3e..the branch head): engine `src/` +352 -657, `include/` +66 -118, `api/` +48 -304, `io/`
+136 -91 (the refusal moved in, ConvertToFormat2 and its test moved out), `experiment/` +56 -4; `convert/` +97 -29;
`pyazr/` +88 -66; tests +332 -165; docs and skills +375 -205. Files of the engine that still name `segmentsData`,
`targetInt`, `format 2` or "old format": `io/src/ProjectFormat.cpp` (the message) and the five pending files above.

## 2026-11-11 (director, 23:10) — merged: one fit-parameter kind for apparatus numbers, the point mapping deleted with the node table beside the sub-points; external capture was missing for a capture term that is only a component

Batch rehearsal `rehearsal/trial25` (worktree `wt-trial25`, build `work/b-trial25`) = the head of 22:00 +
`feature/upgrade-5-apparatus-parameters` + `feature/upgrade-5-grid-nodes` + one fix by the director. Merged into
`feature/upgrade-4-data` and `feature/upgrade`. Nothing pushed.

### The generic apparatus kind (`feature/upgrade-5-apparatus-parameters`) — director's review

- What it does (the agent's entry below has the detail): `ParamKind::Apparatus` replaces `StoppingScale`; a
  parameter of an apparatus effect is a row of `ApparatusRows()` (`include/ApparatusParameters.h`: declare, write,
  read, follow) and its fit plumbing (limits, prior penalty, label, χ² line, hand-over, `parameters.out`) is
  written once in `src/ApparatusParameters.cpp`. The proof is the drill: a triangular resolution kernel with a
  fitted width cost 68 lines in `experiment/src/KernelTriangular.cpp`, 2 lines in `experiment/src/Kernels.cpp`, 1
  CMake line, tests and regenerated documents, and **no file under `src/`, `include/`, `api/`, `pyazr/`** (the
  stopping scale had cost 20).
- Engine diff read by the director. No defect found. What a maintainer should know, not said loudly enough in the
  agent's text:
  - **The limits assume a positive, scale-like number**: [0.1 × nominal, the row's layout]. A parameter that can be
    0 or negative (an energy offset of a detector) needs another limit rule; the place is `ApparatusLayout`.
  - **The derivative is a central difference of the whole χ²** for every apparatus parameter. A kernel's analytic
    derivative with respect to its width exists and is tested (`kernels_test`), but is not wired into the gradient.
    Cost: two full evaluations per apparatus parameter per gradient.
  - **The band does not see apparatus parameters** (as it does not see densities or norms: entry of 22:00).
  - One "warn once" flag is shared by all apparatus parameters: the second parameter that leaves its range does
    not warn.
  - Density, Doppler temperature and layer thickness keep their own kinds (moving them would reorder the
    parameters of every project that fits one, and with them the fits' paths). Decided: not now; after the grid.
- The agent asked whether `Schema::kVersion` should be raised because `width 20 keV` now needs `fixed` or `fit`:
  no — format 3 has not been released, there is nobody to tell apart by a version.
- Left over, small: a build refusal of a detector shared by two effects is printed twice;
  `TargetEffect::GetStoppingScaleMinimum` is unused; `WorkingCopy::SetUpStoppingScaleParams` is an inline forward.

### The point mapping deleted, the node table beside the sub-points (`feature/upgrade-5-grid-nodes`) — director's review

- **This is half of what its brief asked.** Done: the point mapping is gone (`MapData`, `MapLeaves`,
  `SharesCalculation`, the three `AZURE_NO_*_MAPPING` switches, the mapped clauses of 15 engine files: 1,018 lines
  out, 528 in, of which 263 are the node table); the tools the shared mesh needs; a node table built at the end of
  the set-up and compared with the units by `AZURE_NODE_CHECK=1`. **Not done: no reader uses the table.** The
  sub-point units still carry the calculation, so nothing is saved yet and the set-up does one pass more. The
  second half (readers onto the table, sub-point units deleted, the overflow table for re-placed and shifted nodes)
  is brief `AV_grid_part1b_readers.md`, started at once on the rehearsal head (`wt-readers`,
  `feature/upgrade-5-grid-readers`).
- Engine diff read by the director (26 files). No defect found. Outside its brief the agent touched
  `src/AZURECalc.cpp` and `src/AZUREGrad.cpp`: the deleted lines are the mapping's (`CrossSegmentPartners`,
  `MappingConfined`, the second loop over mapped points in `ComputeTHMRows`); accepted.
- `AZURE_REFINE_QUADRATURE=k` is read in `Config::Reset` and ends the process on a bad value. Fine for the command
  line; in pyazr a bad environment variable would end the Python process. To become a flag when `src/AZURE2.cpp`
  has room (it is at its size limit), or an exception.
- Expected output changes, seen by the director: only the "Coulomb table: N lookups" line of `run.log` in
  `angle_mapping` and `sum_lab_frame/models/plain_data` (a work counter).

### External capture was missing for a capture term that is only a component (fixed; older than the upgrade)

- Reported by the grid agent as "a SUM extrapolation of different pairs may be wrong" (elastic + 1000 × capture at
  90° c.m. on 13N: the capture term 0.54 to 2.9 times the separate capture extrapolation, energy dependent). It
  did not investigate.
- Cause, found by the director: the command line switches external capture on when a level with external capture
  exits into a pair that some active data set or extrapolation uses, and it listed **only the leading term** of
  each reaction (`ProjectReactionPairs`, `ProjectExitPairKeys` in `src/ProjectSetup.cpp`). A project whose capture
  pair appears only as the second term of a SUM or RATIO ran **without external capture, without a word**. With
  any other data set or extrapolation that has the capture as leading term the same SUM was right — which is why
  the old check never saw it (its model had both), and why it appeared when the agent's rewritten check changed
  the model. The decision is the same in the code at the branch point with `dev` (`readSegmentFile` reads the
  leading pair of each segment line). The API (pyazr) always switches external capture on and was not affected.
- Fix: every term's pair is listed (`ForEachTermPair`). Test: `tests/component_mapping` (d) — the SUM as the only
  extrapolation equals elastic + 1000 × capture of the separate runs to the digits of the file (7.6e-11 of the
  sum) and writes `intEC.extrap`; on a build without the fix it fails (7.58e-4 of the sum, no `intEC.extrap`).
- No output of the tree changes with the fix (all projects × 3 bases identical before and after): no project of
  the tree had a capture pair only as a component. A user's project could.

### Verification by the director (build `work/b-trial25` at the rehearsal head, 0 warnings)

- All projects × 3 bases against the head of 22:00 (`tab-trial24c`): 778 / 751 / 778 files identical; different are
  the two `run.log` work-counter lines named above (per basis) and 6 files that are new
  (`rectangular_resolution/models/fit.azr`, the width fit). With the external-capture fix against the run without it
  (`tab-trial25b`): 786 / 759 / 786 identical.
- `ctest -LE slow` 203 of 203; `ctest -L slow` 10 of 10; `tests/run_tests.sh` 153 passed, 0 failed; spec equivalence
  184 of 184 converted files equal; the converted-tree gate against `tab-trial8` 0 differ, 0 missing (38 classified
  per basis, as before).
- Calculations at 1 against 3 threads: T1, T2, T3, 13N × 3 bases and `T2_true` with model data identical.
- Fits at 1 thread against the head of 22:00 (build `work/b-trial24`), every output file identical: T2 with
  `--use-lm`, `--use-gradient`, MIGRAD at 200 calls; 13N with `--use-gradient`, `--use-lm`;
  `layer_yield/models/fit.azr` with `--use-gradient`; `stopping_yield/models/fit.azr` (two areal densities and the
  stopping scale, now of the generic kind) with `--use-gradient`, `--use-lm`. At 2 threads against 1: T2
  (`--use-gradient`, `--use-lm`) and the stopping fit identical in every final output (the intermediate snapshot
  `param.fit` of the stopping fit differs: known book-keeping, see FOR_THE_OWNER "Fits and the number of threads").
- `natural_neon`, calculation at 1 thread, 6 alternated runs, median: 6.40 s on the build of 19:25
  (`work/b-trial23`), 6.43 s on the head of 22:00, 6.51 s now (+1.6 % and +1.2 %: the node table's set-up pass).
  Accepted until the readers are on the table; brief AV asks for the time of `work/b-trial23` back.
- `tests/component_mapping` (d) fails on the grid agent's own binary (no fix) and passes on the merged head.
- Under the 12-thread limit of this night: 6 threads for the director's runs, 3 for each agent.

### Next

- Running: the readers onto the node table (brief AV, 3 threads); the engine refusing old files (brief AU, 3
  threads). After AV: the shared mesh (decision D2) with decision D3's deletions, on the node planner alone.
- For the shared mesh's brief: the quadrature gate's "fine" reference should be k = 4 (the agent's classifier
  failed a k = 4 run on 5 of 92 files per basis with k = 2 as reference: `natural_neon`, `T2_true`,
  `layer_yield` const are not asymptotic at k = 2).
- `tests/fit_modes` has a fifth case since 23:30, "apparatus" (the rectangular width alone, four routes, band skipped
  with its message): 37 checks pass.

## 2026-11-11 (agent, Phase 5 branch 7 part 1) — point mapping deleted, tools for the shared mesh, the node table beside the sub-points (`feature/upgrade-5-grid-nodes`)

**Done.**
- The point mapping is deleted (`MapData`, `MapLeaves`, `SharesCalculation`, `AZURE_NO_*_MAPPING`, the mapped clauses
  in 15 engine files; one plan task per point, `WorkingCopy::BuildPlan`). Checks that compared with/without the
  switches now test what is still true (`component_mapping`, `angle_mapping`, `shift_fn_groups`, `sum_lab_frame`,
  `composite_gradient_test`, `failed_point_test`, `plan_test`); consistency check 29 forbids the deleted names.
- Tools for part 2: `AZURE_REFINE_QUADRATURE=k` (every sub-point grid k-fold finer; an environment variable, not the
  flag of `energy_grid.md` §7.2, because `src/AZURE2.cpp` is at its size limit), `tools/classify_quadrature.py`,
  `tests/upgrade/node_count` (node counts per folded project: table in `energy_grid.md` §8, largest mesh/today 1.07),
  `tests/layer_yield` (K') the thick layer against an independent integral (engine +1.54e-4, recorded as the bound;
  second order in the node spacing).
- The node table beside the sub-point units (`include/NodeTable.h`, `AZURE_NODE_CHECK=1`): clean on the whole tree.

**Verified** (logs in `work/nodes/logs/`, binary `work/nodes/AZURE2-m2` = 2feb5ce0): tree x 3 bases after the
deletion (`tab-nodes-m0`) identical but the Coulomb-lookup line of `run.log` in `angle_mapping` and
`sum_lab_frame/models/plain_data`; the tree with `AZURE_NODE_CHECK=1` (`tab-nodes-m2chk`): every output file identical,
0 rows and 0 slots differ in 271 passes; `ctest -LE slow` 203/203; `tests/run_tests.sh` 152 passed, 0 failed; the
converted-tree gate 483 (456) identical, 38 classified, 0 differ, 0 missing; fits (T2 `--use-gradient`, `--use-lm`, 13N,
`layer_yield/fit.azr`) identical at 1 and 2 threads and to `b-trial23` at 1 thread; calculations of eight folded
projects identical at 1 and 3 threads. Speed (menu 1, 1 thread, 6 alternated runs, load from other agents): T2, 13N,
thick target, THM within ±2 %; `natural_neon` median +3.0 %, paired mean +1.2 % (the node table's set-up pass; not
resolved from the noise). Peak RSS unchanged within 0.5 %. Not run: Sphinx (no sphinx-build in the environment),
ThreadSanitizer, fits at 3 threads (budget 2 until late).

**Open (part 2).** The node table's readers are not moved: the energy pass, the kernels (`PointFold.cpp` gathers from
sub-points), the gradient's cotangents and `Results`' sub-point records still use sub-point `Unit`s; no overflow table
(re-placed and shifted nodes are "off the table"); `SubPointAccess.h`, `NodeTableCache.h`, copy-on-write tables not
deleted. The planner to replace: the two grid loops of `WorkingCopy::SetUpTargetEffects` and
`ReplaceSubPointGrids` (`src/TargetSetup.cpp`). `classify_quadrature.py` with k = 2 as "fine" fails a k = 4 run on 5
of 92 files per basis (non-asymptotic values): part 2's gate may need k = 4 as its reference.

## 2026-11-11 (agent, apparatus parameters) — one fit-parameter kind for the numbers of apparatus records (`feature/upgrade-5-apparatus-parameters`)

**Done.**
- Map of the path a fitted apparatus parameter takes (density, Doppler, thickness, stopping scale), the generic kind,
  what the other three would need, the measured cost: `docs/upgrade/apparatus_parameters.md`.
- Kind 12 is `apparatus` (it was the stopping scale, one day old): rows in `ApparatusRows()`
(`include/ApparatusParameters.h`),
  `stopping_scale` (`src/StoppingScale.cpp`, now its row and its follow-up) and `kernel` (any table kernel's parameter
  with `fittable`); one engine file `src/ApparatusParameters.cpp` (groups by sharing key, vector slots, hand-over,
  prior and its `chiSquared.out` line, labels, index map); every switch over the kinds has one case `Apparatus`.
  Names: the owner's values joined by `.` (`stopping_scale.TiN.1H`, `det_fit.energy_resolution.width`). Same code and
  place in the vector (last); API key `apparatus` instead of `stopping_scale`.
- Format 3: a fittable kernel parameter is generated with the nuisance grammar (`width 14 keV fit prior 14 2 grid_max
30`);
  a parameter only when `fit`, a prior or grid_max is written; the grid spans the layout value
  (`KernelParameters::layoutValues`, read by `EffectSupport`). `tests/rectangular_resolution/models/rect.azr` gains
  `fixed` (the grammar has no default mode); its outputs are unchanged.
- Proofs: `tests/rectangular_resolution/check.sh` (c)-(g): the width fitted from 1 %-noise pseudo-data at 20 keV
  (20.026 ± 0.013 keV), noise-free (19.999998 keV), with a prior, one parameter for two data sets/effects, six refusals;
  the stopping scale's fit byte-identical to the base build.
- Drill: a triangular resolution with a fitted base width (`experiment/src/KernelTriangular.cpp`), 0 engine files,
  about 6 minutes (commit bea56a6f).

**Verified** (binary `work/b-appar`, head bea56a6f; scratch `work/appar/`): all projects × 3 bases against
`tab-trial23g`
(tag `tab-appar-final`, 2 jobs): 780 / 753 / 780 identical, the new `fit.azr`'s 6 files per basis only-new;
`ctest -LE slow` 203/203; `tests/run_tests.sh` 152 passed, 0 failed; fits at 1 thread identical to the base build's
(director's `repro23`: T2 --use-gradient, MIGRAD, --use-lm, 13N, stopping_yield fit, layer_yield fit; fitStatus.out
and param.fit excluded); the same fits at 1 and 2 threads identical (3 threads not run: budget 2 at the time);
converted-tree gate (`REF=tab-trial8`, tag `gate-appar`, 3 jobs) 483 (456) identical, 38 classified, 0 differ, 0
missing. Injected errors shown to fail: the fold not reading the fitted width, one
parameter per effect, the stopping scale's change not followed.

**Open.** The kernels' analytic `derivative` is not used by the gradient (central differences of the whole χ²; the
cotangent list `'a'`, `'d'`, `'t'` of `GradEnergyTasks.cpp` is untouched); density, Doppler, thickness keep their own
kinds; a refusal of a shared record is printed once per effect (twice for two effects); `TargetEffect::
GetStoppingScaleMinimum` is now unused; Sphinx not built. For the merge with the node-layout branch: one line in
`src/PointFold.cpp`; `WorkingCopy::SetUpStoppingScaleParams` forwards to `SetUpApparatusParams` (its call is in
`src/TargetSetup.cpp`).

## 2026-11-11 (director, 22:00) — a check of every fit route on every kind of free-parameter set; a second crash found by it, older than the upgrade

Branch `feature/upgrade-5-fit-modes-check` (the director's own, from the head of 21:00), merged into
`feature/upgrade-4-data` and `feature/upgrade`. Nothing pushed.

- **`tests/fit_modes/check.sh`** (new; run by `tests/run_tests.sh`; about 95 s at 2 threads). The four minimizer
  routes (MIGRAD, `--use-gradient`, `--use-lm`, `--use-gsl-lm`), each with `--covariance-band`, on four projects
  that differ in what is free: none (`thick_target_resonance`), level parameters and norms (`upgrade/T1`), level
  parameters and an areal density (`layer_yield/models/fit.azr`), two areal densities and a stopping scale
  (`stopping_yield/models/fit.azr`; `--use-gradient` and `--use-lm` only). Checked: exit 0 and a total χ²; with no
  free parameter the χ² is the calculation's digit for digit; otherwise the routes end at the same minimum (1e-5
  relative; they agree to 1e-7 today: T1 462.27351–462.27354, layer yield 39.30580, stopping yield 72.03660–72.03661);
  the band is written where a level parameter is free and skipped with a message that says why where none is.
- **It fails where it should**: on the build of the head of 19:25 (b8c970b0) in all four routes of "none" (10
  failed checks); on the build of the head of 21:00 in "none" with MIGRAD, plain and `--use-gradient` (4 failed
  checks).
- **The second crash, older than the upgrade.** A MIGRAD fit with `--covariance-band` of a project without free
  parameters ended in a segmentation fault inside Minuit2 (`MnHesse` on a state with no variable parameter). The
  code is the same at the branch point with `dev` (aab1de0d). `src/AZUREMain.cpp` now skips HESSE when the
  minimum has no variable parameter (Minuit2 is not patched). The crash fixed at 21:00 was the least-squares
  route's; the director's reproduction that evening ran MIGRAD without the band and therefore did not see this one.
- **The band's message.** A fit that leaves no level parameter free (only target parameters, or nothing) got "this
  fit produced no usable covariance", which is wrong: the covariance is fine, the band does not use it. It now
  says: "Uncertainty band skipped: the fit has no free level parameter. The band propagates the covariance of level
  energies, widths and the THM coherent background; the parameters of targets and data sets are not in it."
- **A design gap, reported, not changed**: `BandPackedColumns` (`src/CovarianceBand.cpp`) keeps the kinds
  `LevelEnergy`, `Gamma`, `ThmCoherent`. A fitted areal density, stopping scale, Doppler temperature or norm does
  not widen the band of a yield, and a fit of such parameters alone gets no band. It is also one more place that
  names parameter kinds (for the review of the generic apparatus kind). A question for the owner.
- **19 more pyazr tests passed without running** when the engine module could not be imported (their skip path
  was `return 0`; the 33 fixed at 21:00 were the ones with `sys.exit(0)`). All now return 77, which ctest shows as
  skipped; `tests/CLAUDE.md` states the convention with the number.
- Verification by the director: build 0 warnings; the check on three builds as above; the consistency test;
  `ctest -R "pyazr|consistency|upgrade|cli"` on the new build: 125 of 125 pass, none skipped; no committed
  expectation contains the old message.

## 2026-11-11 (director, 21:00) — merged: pyazr reads format 3 only, the audit of the checks with six refusals, the crash without free parameters fixed

### State

- `feature/upgrade-4-data` holds three merge commits (`feature/upgrade-4-pyazr-format3`,
  `feature/upgrade-4-check-audit`, the last two commits of `feature/upgrade-5-energy-tasks-gradient`) with the
  director's fixes of the test harness inside the first, and this documents commit. Its tree before the documents
  commit is the rehearsal's (`rehearsal/trial24`, 58db64bb). `feature/upgrade` is merged with it. Nothing pushed.
- Phase 4: what is left is the refusal switch (the engine refuses old files, its old reading path and the bindings
  pyazr no longer calls are deleted; brief `work/briefs/AU_refusal_switch.md`, `feature/upgrade-4-refusal-switch`,
  started from this head) and the rest of `CNuc`.
- Phase 5, running since the afternoon, both from the head before this one (b8c970b0): branch 7 part 1 (the node
  table, `feature/upgrade-5-grid-nodes`) and the generic fit-parameter kind for apparatus parameters
  (`feature/upgrade-5-apparatus-parameters`). Both were told of the crash and not to fix it themselves; their merges
  bring this head's fix.
- Machine: 12 threads in total this night (PI). Agents 3 each; the director's verification is
  `work/director/full_verify4.sh` (4 threads).

### Director's verification

The rehearsal with the pyazr branch, the harness fixes and the audit branch (36f3d85e, build `work/b-trial24`):

| check | result |
|---|---|
| all projects, 3 bases, against the last verified batch (`tab-trial24` vs `tab-trial23g`) | 780/753/780 files; 779/752/779 identical; the one that differs per basis is the intended `run.log` of `relativistic/models/li7pn_edge_rel` (the set-up refusal is its last line) |
| calculations at 1 against 3 threads (T1, T2, T3, 13N; 12 files, χ²) | identical |
| `ctest -LE slow` | 202/202, none skipped |
| `ctest -L slow` | 10/10 |
| `tests/run_tests.sh` | 152 passed, 0 failed |
| old against format 3 (`spec_equivalence`) | 206 files: 184 converted and equal, 0 differ, 22 not converted |
| converted-tree gate against `tab-trial8` | 483 (456 formal) identical, 38 classified, 0 differ, 0 missing |
| build | 0 warnings |

Then the crash fix on top (58db64bb, the same build directory rebuilt, 0 warnings):

| check | result |
|---|---|
| all projects, 3 bases, against the table above (`tab-trial24c` vs `tab-trial24`) | 780/753/780 files, all identical |
| `ctest -R "upgrade\|grad\|consistency"` | 20/20, with the new `upgrade_energy_task_gradient` |
| `thick_target_resonance` (no free parameter) as a fit with `--use-lm`, `--use-gradient`, plain | exit 0, χ² 30.909423272 in the three (the head before: exit 139 with `--use-lm`) |
| T2 and 13N fits with `--use-lm` and `--use-gradient`, run to convergence, 1 thread, against the build of the head before (`work/b-trial23`) | every output file identical (4 fits) |

Not repeated after the fix: the full `ctest`, `run_tests.sh`, the gate (the fix changes one function that only a
gradient or a row pass reaches; the documents commit after it was checked with the consistency test).
Logs: `work/director/v4_trial24.*.log`, `gate_trial24.log`, `cls_trial24c.log`.

### Open, for the owner or the PI

- The crash (next section) is the second time in this batch that "all tests pass" said less than it seemed (the
  first: tests that passed without running). Proposed and not started: a small matrix of the tree's fit modes
  (plain, `--use-gradient`, `--use-lm`, band) on one project of each kind, the no-free-parameter one included, as a
  ctest of its own.
- New refusals a user can meet: `E_THM_OLD_LIST`, `E_LAB_JACOBIAN_AT_LIMIT`; a refused format-3 set-up ends with
  its own ERROR line. Listed in `FOR_THE_OWNER.md` ("What users notice").
- pyazr 3.0.0.dev0 refuses an old file; the engine still reads it until the refusal switch is merged.
- Found by the director's own slip: `--max-calls 200` (a space instead of `=`) prints "WARNING: Unknown option
  --max-calls." and the fit runs without a limit. An option the engine does not know should end the run; not changed.

### Fits and the number of threads, and a crash the director's last verification missed

- **A crash was in the head merged at 19:25** (data branch b8c970b0): a least-squares fit (`--use-lm`) or a
  covariance band of a project with **no free parameter** ended in a segmentation fault (`EnergyTaskRows` sized its
  buffer by a column count of 0). Reproduced by the director with `thick_target_resonance --use-lm` on the merged
  build; the build before the gradient branch runs it. The gradient agent found it in its own speed series and
  fixed it (`src/GradEnergyTasks.cpp`: with no packed column the row pass returns the models only) with a test
  (`upgrade_energy_task_gradient`: the scalar gradient, the residual Jacobian and the model gradients give the same
  bits at 1 and 3 threads on T2; the no-free-parameter case runs; it segfaults with the fix removed). The fix is in
  this batch. **Why the verification missed it:** no test of the tree ran a least-squares fit without free
  parameters. The director's verification compares what the tests cover; it does not find what no test covers.
- The speed of the gradient branch, the agent's series (1 thread, alternated with the build before it, machine load
  10–13, read by the director from the raw timing files, not repeated): 13N `--use-gradient` 50.7 → 26.9 s, 13N
  `--use-lm` 19.5 → 12.8 s, `natural_neon` `--use-gradient` (30 calls) 219.6 → 105.7 s, `--use-lm` 291.1 → 139.3 s;
  T2 (plain MIGRAD, 200 calls) 39.9 → 37.2 s at 1 thread and 20.7 → 19.1 s at 3. The `thick_target_resonance`
  numbers of that series (13.8 → 5.2 s, 7.9 → 4.6 s) are **not fits**: the project has no free parameter. The
  5-thread T2 timing is still owed.
- **Memory and time since the start of the upgrade, measured by the director** (one `natural_neon` calculation, 1
  thread, both runs one after the other under the same load): the engine at the branch point (aab1de0d, build
  `b-base`) peaks at 4.52 GB and takes 41.2 s; this batch's build peaks at 1.46 GB and takes 8.0 s; χ²/N of the
  first segment 10.9941 in both, total χ² 3342.46 against 3342.51.

### pyazr reads format 3 only (`feature/upgrade-4-pyazr-format3`, Phase 4) — director's review

- What it does: `pyazr/azrfile.py` (3,756 lines, the old file class), `pyazr/datasets.py` and
  `tests/pyazr/old_tree.py` are deleted; version `3.0.0.dev0`. A file that is not format 3 raises `FormatError`
  with the engine's message (`pyazr.project.require_format3`); `pyazr.convert_azr` (default `to=3`) is the only
  reader of an old file. What the old class could do and scripts need exists on records: `masses.mass_check` /
  `fill_masses`, `set_pair(light_excitation=, heavy_excitation=)`, `snapshot.place_level_values`,
  `outputs.fitted_norms`, `load_sav` with the engine's check of the parameter names, and `merge_projects` on format 3.
- Engine diff read by the director: one new piece, `experiment/LevelEdits.h/.cpp` (`MergeLevelRows`: the level rows
  of several projects in one section, pair keys and level IDs shifted; binding `merge_level_rows`). No defect found.
- Tests: of the 54 test files of the old class 46 were moved onto records, 7 are kept as converter tests, 1
  (round trip of the old writer) is deleted with its subject.
- **Coverage lost, listed by the agent** (cases only the old class could build): two different targets tied to one
  density; a thick-target yield of a sum; an inactive scattering-length constraint line; a follower's own ignored
  norm; the `prior_centre` rows of `save_fit_norms`; old byte layouts; some THM refusals at session start; an
  extrapolation without a frame refused; `light_excitation 0` written as "0". None is a format-3 feature without
  a test elsewhere, by the agent's list; the director did not check each.
- **Defects the agent found and did not fix** (engine side: given to the refusal-switch branch, brief AU):
  `Project.set_config(output_dir=)` neither appends "/" nor creates the directory; a refused THM edit repeats its
  prefix; `Project.set` on an unknown record raises `KeyError`; the out-of-range warning of a fitted density uses
  old names; the CLI cross-checks of `sav_fixed_params` and `inactive_lines` search `build*/` and skip in a
  separate build directory; the skills `azure2-setup` and `azure2-targets` describe the old grammar. Not given to
  anyone yet: `check_format3_thm` accepts combinations the session refuses.
- Bindings left without a Python caller (`check_thm_block`, `read_thm_block`, `thm_option`, `thm_spin_text`,
  `thm_global_optical_for`, `read_thm_weight_table`, `parse_segment_line`, `parse_extrapolation_line`,
  `parse_target_effect_line`, `LevelGrouper`): deleted with the engine's old reading path (brief AU).
- Director's decision: the edit of `tests/mass_separation_check/check.sh` (it calls `masses.fill_masses`) is
  accepted; run by the director with and without the engine module: all checks pass, part (e) runs with it.
- **Two defects of the test harness found by the director in the rehearsal and fixed there** (tests that passed
  without running): 33 `tests/pyazr/*_test.py` ended with exit 0 when the engine module could not be imported
  (now 77, which ctest shows as skipped, the convention the other tests state); `pyazr_layer_yield` was registered
  without `PYAZR_ENGINE_DIR` and had been skipped in every build outside the source tree (now it gets the
  directory, runs and passes).

### The audit of the checks and six refusals (`feature/upgrade-4-check-audit`, Phase 4) — director's review

- The question was the director's, from the THM conversion: a check that edits a project's text with `sed` can go
  on passing after a conversion because its pattern no longer matches anything and it compares a run with itself.
  The agent answered it mechanically (`tools/edit_audit.sh`: every `check.sh` run with logging stand-ins for
  `sed`/`awk`/`perl`/`python3` and for the engine; `tools/edit_audit_report.py`): 126 checks, 1,658 edit places,
  40 that change nothing, of which 36 by design and **4 vacuous** (`target_effect_shared`, `thm_beam_profile`:
  both edited the old `<targetInt>` block); one more check tested nothing for another reason (`nuisance_kinds`:
  its `<parameterSettings>` row named a parameter that no longer exists, and the engine ignores such a row
  silently). The three are rewritten on records and each fails with an injected engine error. The 25 pyazr tests
  of the old class were out of its scope (the pyazr branch moved them; not audited since).
- The guard that keeps it so: `tests/lib/edit_guard.sh` (`edited`, `must_differ`: an edit that changes nothing
  fails the check), rule 33 of the consistency test (no unguarded in-place edit in a `check.sh`), the `CLAUDE.md`
  recipe "A new regression check".
- Six refusals, engine diff read by the director (338 lines): `E_THM_OLD_LIST` (a format-3 `<thm>` with
  `segments=` or `weight[k]=` is refused at reading); an experiment no data set names is ignored with one warning;
  one reading of an experiment's name (`experiment/ThmSectionLines.h`); a set-up refusal of a format-3 project
  ends with its own ERROR line (`kFillRefused`; "Could not fill data object from file." no longer follows it: the
  one `run.log` of the tree that changes, `relativistic/models/li7pn_edge_rel`); the converter keeps the sample
  target on an inactive line or refuses with the reason; `E_LAB_JACOBIAN_AT_LIMIT` (a lab-frame differential point
  where dΩ_lab/dΩ_cm ≤ 1e-6 is refused at set-up; without it such a point gave χ² 1.06e17). No defect found in the
  reading. With the stopping codes of the data head the table has 93 diagnostic codes.
- Open, given to brief AU: the silent `<parameterSettings>` row becomes a warning; the THM pages of `docs/source`
  still show `segments=`. Noted: three checks use `\|` in `sed` basic expressions (GNU only); the audit tool's
  key for an engine run drops lines with '/', so refusals that differ only in a file name group together; a full
  audit takes about 2.5 h at one thread.
- For the merge (director): the conflicts of this branch with the data head were lists two branches extended
  (the stopping refusals and this branch's codes in `Diagnostics.h/.cpp`, the reference page, the schema, the
  CMake source list): both sides kept (`work/director/union_conflict.py`), the schema and reference-page tests
  pass on the result.

## 2026-11-11 (agent, check audit) — do the checks still test what they say? Six refusals (`feature/upgrade-4-check-audit`)

**Done** (`AR_check_audit.md`; the audit in `docs/upgrade/format3_tree.md` section 7):
- **The audit**, mechanical: `tools/edit_audit.sh` runs every `tests/*/check.sh` with `sed`/`awk`/`perl`/`python3`
  and the binary replaced by logging stand-ins (`tools/edit_audit/`), `tools/edit_audit_report.py` prints the table.
  Run once over the head `7a8e4d66` with the base engine (`b-trial22`), one thread: 126 checks, all pass; 103 edit a
  project, data or parameter file; 1,658 edit places, 40 hit nothing (36 by design: `\r` stripping, the A-matrix
  flag set to its value, base variants; **4 vacuous**: `target_effect_shared` and `thm_beam_profile` edited the old
  `<targetInt>` block); 1,470 engine runs, 48 checks with runs of different inputs and equal outputs, all asserted
  identities but one (**`nuisance_kinds`**: its `<parameterSettings>` row named `screening_1`, which no parameter is
  called any more: ignored silently). In-scope `tests/pyazr/*.py` text edits checked by hand: all hit.
- **Repaired** on records: `target_effect_shared`, `thm_beam_profile`, `nuisance_kinds`; each passes and fails with an
  injected engine error (section 7). Expectations unchanged; no repaired check shows a number different from its
  committed expectation.
- **Guard**: `tests/lib/edit_guard.sh` (`edited`, `must_differ`); consistency check 33 (in-place edits in a check.sh);
  `CLAUDE.md` recipe "A new regression check".
- **Item 4**, one commit: `E_THM_OLD_LIST` (segments= / weight[k]= in a format-3 `<thm>`, at reading); an experiment
  no data set names is ignored with the warning (was "segments= is required"); one experiment-name grammar
  (`experiment/ThmSectionLines.h`), the reading-time `E_THM_EXPERIMENT_UNKNOWN` the user's message, the run-time test an
  internal guard a project file cannot reach; a set-up refusal of a format-3 project ends with its own ERROR line
  (`kFillRefused` from Fill/MakePoints); the converter writes the sample target on an inactive sample line, refuses one
  whose sample has no target of its own; `E_LAB_JACOBIAN_AT_LIMIT` (dOmega_lab/dOmega_cm <= 1e-6 at set-up; without it
  a point at theta_max gave chi2 1.06e17). 88 diagnostic codes.

**Verified how**: build 0 warnings; unit tests of the change (experiment_frames, _diagnostics, _build, _schema,
_standalone, convert_constructs, _fixtures, _spec_equivalence, pyazr_schema, repo_consistency); the six checks touched
pass, and fail with their injected errors (one scratch rebuild of `work/b-audit` with all injections, reverted and
rebuilt; not a separate `b-audit-inject`: one thread). At the end (build of `3ec01e1e`'s engine, identical to
`2d4c085e`'s): all projects x 3 bases against `tab-trial22` (tag `tab-audit1`): 765/766, 738/739, 765/766 identical,
the one file per basis `relativistic/models/li7pn_edge_rel.azr/run.log` ("Could not fill data object from file." gone
after the set-up refusal, item 4); `ctest -LE slow` 200/201 before the last commit (upgrade_plan counted a refused
project as set up: `Fill == -1`; fixed with five other tests in `3ec01e1e`, the six pass); converted-tree gate
(`gaudit1`) 483 (456) identical, 38 classified, 0 differ, 0 missing; `tests/run_tests.sh` 148 passed, 2 failed, both
pass alone with this build (`cli_help`, known under load; `shared_point_tables` (d): the calculation's peak RSS read
39 MB under load against 184 MB alone, base and new alike).

**Open**: a `<parameterSettings>` row that names no parameter is ignored without a word (found by the audit; engine,
not changed); the audit's stdout key drops lines with '/' (paths), so refusals with different file names group
together (read the check); pyazr tests of the old file class not audited (brief AQ); docs/source pages on THM still
show `segments=` (format3_tree.md section 6).

## 2026-11-11 (pyazr off the old file class, branch `feature/upgrade-4-pyazr-format3`) — pyazr reads format 3 only

**Done.** `pyazr/azrfile.py` (`AzrModel`, 3,756 lines), `pyazr/datasets.py` (`SegmentSet`), `tests/pyazr/old_tree.py`
and the old generators of `doppler_fgm` / `transmission_yield` (under `tests/convert/fixtures/`) are deleted; pyazr
is 3.0.0.dev0 (the version lives in `pyazr/__init__.py`, which `pyproject.toml` reads; the consistency test now
accepts a `.devN` suffix). An old file handed to a session, `Project.open` or any function that opens a project raises
one `FormatError` with the engine's message naming `AZURE2 --convert OLD NEW --to 3` (`project.require_format3`,
binding `not_format3_message`). Every pyazr test that used the old class is moved to the tree's format-3 project
(same committed numbers), kept as a converter test, or (only `roundtrip`) deleted with its counterpart named: the
table and the coverage lost are in `pyazr_on_records.md`. New: `experiment::MergeLevelRows` (binding
`merge_level_rows`) and `merge_projects` on format 3; `SetPair` sets the excitation energies
(`set_pair(light_excitation=,
heavy_excitation=)`); `masses.mass_check` / `fill_masses` on a Project; `snapshot.place_level_values`;
`set_thm_experiment(datasets=)`; `load_sav` matches names by `experiment::MatchParameterFileNames` (binding
`match_parameter_file_names`); `convert_azr` in `pyazr/convert.py` (default `to=3`); `fitted_norms` in
`pyazr/outputs.py`; the examples, `scripts/thm_model_average.py` and the three resonance-search scripts use Project.
Documents: `pyazr_on_records.md` (done), `format3_tree.md` §5 (what is left: the engine and the GUI), the user guide,
the skills, `CLAUDE.md`, `FOR_THE_OWNER.md`.

**Verified** (own runs, 2 threads in all): all projects x 3 bases against `tab-trial22`: 766 / 739 / 766 identical
(`tab-pyoff1`); the converted-tree gate (REF tab-trial8): 483 (456 formal) identical, 38 classified, 0 differ, 0
missing; `ctest -LE slow` 200/200 (one fix on the way: `pyazr_park_dsde` composed an old model at run time, now
converted); `ctest -L slow` 10/10; pyazr tests 112 (113 before, `roundtrip` deleted), all pass; `tests/run_tests.sh`
150/150; consistency test passes. Injected errors (each failed, then reverted): natural_target (no pair-key offset in
the merge), thm_experiment, norm_group, prior_covariance, screening. Not run: 1 against several threads (no engine
computation changed), Sphinx, the GUI, `pyazr/examples/exfor_fetch.py` and `save_fit_to_azr.py` (network / a
parameter file; the other six examples were run once).

**Open.** The engine's refusal of old files (next branch); bindings with no Python caller left (listed in
`pyazr_on_records.md`); defects found: `Project.set_config(output_dir=)` neither appends "/" nor creates the
directory; a refused THM edit repeats its prefix; a new THM experiment without data sets is refused with old wording
("segments= is required"); `Project.set` of an unknown record raises KeyError, not FormatError; the engine's
out-of-range warning of a fitted density still says `density_1` and "densityfit"; `sav_fixed_params` and
`inactive_lines` look for the CLI under `build*/` and skip their cross-check in a separate build directory.

## 2026-11-11 (director, 19:25) — merged: kernels, the stopping wiring, the gradient per compound energy (Phase 5 branches 6, 5, 2)

### State

- `feature/upgrade-4-data` holds the three branches (merge commits of `feature/upgrade-5-kernels`,
  `feature/upgrade-5-stopping-wire`, `feature/upgrade-5-energy-tasks-gradient`) and two director's commits (the pyazr
  test-registration block last again; the documents on the kernel count's limits and the new refusals). Its tree is the
  rehearsal's (`rehearsal/trial23`), which is what was verified. `feature/upgrade` is merged with it. Nothing pushed.
- Phase 5 (`energy_grid.md` §8): branches 1 (one solve per energy), 2 (gradient), 3 (lab frame), 4 (layer yield),
  5 (stopping wiring), 6 (kernels) are merged. Running: branch 7 part 1 (the node table; sub-points stop being
  `Unit`s, the point mapping deleted, bit-identical; brief `work/briefs/AT_grid_part1_nodes.md`,
  `feature/upgrade-5-grid-nodes`) and one generic fit-parameter kind for apparatus parameters (brief
  `work/briefs/AS_apparatus_parameters.md`, `feature/upgrade-5-apparatus-parameters`). After them: branch 7 part 2
  (the shared mesh of decision D2, the deletions of D3, kernels that compose), 8 (shift kernel), 9 (observables, D5),
  10 (THM on the grid).
- Phase 4 close, running: pyazr off its old file class (`feature/upgrade-4-pyazr-off-azrfile`), the audit of the
  regression checks (`feature/upgrade-4-check-audit`). Then the refusal switch (old formats refused, the old reading
  path deleted).
- Machine: 12 threads in total for all agents and the director this night (PI, 2026-11-11); `work/briefs/COMMON.md`
  says so. `work/director/full_verify3.sh` is the 6-thread verification; `full_verify2.sh` (12 threads) is not to be
  used until the PI lifts the limit.

### Director's verification (the tree of the last merge, build `work/b-trial23`)

| check | result |
|---|---|
| all projects, 3 bases, against the last verified batch (`tab-trial23g` vs `tab-trial22`) | 780/753/780 files; 766/739/766 identical, 14 only-new per basis (the two new test directories), 0 differ |
| calculations at 1 against 3 threads (T1, T2, T3, 13N; 12 files, χ²) | identical |
| `ctest -LE slow` | 203/203 |
| `ctest -L slow` | 10/10 |
| `tests/run_tests.sh` | 152 passed, 0 failed |
| old against format 3 (`spec_equivalence`) | 206 files: 184 converted and equal, 0 differ, 22 not converted |
| converted-tree gate against `tab-trial8` | 483 (456 formal) identical, 38 classified, 0 differ, 0 missing |
| build | 0 warnings |

Logs: `work/director/v3_trial23.*.log`, `gate_trial23.log`.

### Fits and the number of threads: what holds now, measured by the director

`work/director/repro23/` (script `run.sh`, table `result.txt`): six fits (T2 with the analytic gradient, T2 MIGRAD,
T2 `--use-lm`, 13N with the gradient, `tests/stopping_yield/models/fit.azr`, `tests/layer_yield/models/fit.azr`),
each at 1 thread and twice at 3 threads.

- **Every result file is identical at 1 and 3 threads and from run to run** (parameters, covariance, χ², the data
  output). Before the gradient branch two runs of one build at 5 threads differed in 8 files.
- **Not identical:** the count `calls` in `fitStatus.out` (T2: 205, 204, 203), the order of the "Iteration" lines of
  the log, and the intermediate snapshot `param.fit`. Cause, read in the code: Minuit2's numerical-gradient loop
  (`external/minuit2/src/Numerical2PGradientCalculator.cxx`) calls the objective from an OpenMP loop and its call
  counter (`MnUserFcn.cxx`, `fNumCall++`) is not atomic; the engine's iteration snapshot is written in the order the
  calls arrive. The numbers of the fit do not depend on it. A fit that ends at `--max-calls` could in principle end
  one iteration apart. The PI decided not to patch Minuit2 (`external/`). A remedy on the engine's side (count the
  calls in the engine's own wrapper, write the snapshot of the lowest χ² instead of the last arrival) is not started.
- For agents and checks: compare fits on the result files, not on `fitStatus.out` or `param.fit`.

### Open, for the owner or the PI

- SRIM's terms (free for non-commercial use, citation) against a GPL repository that ships SRIM's two data files:
  the owner's decision (`stopping/data/srim/README.md`).
- A new kind of fit parameter costs 20 existing files; a new apparatus kernel costs none but does not compose and
  cannot be fitted yet. Both are what the two running branches are for; the count will be repeated on them.
- The accuracy of a yield over a narrow resonance on a thick layer is 1.5e-4 (the present quadrature); the grid
  branch gets the director's integral as a committed test.
- The Minuit2 call counter (above).
- ThreadSanitizer has not been run: there is no clang on this machine.

### The gradient per compound energy (`feature/upgrade-5-energy-tasks-gradient`, Phase 5 branch 2) — director's review

- What it does: the analytic gradient, the least-squares rows and the model gradients run over the same energy
  tasks as the values (`src/GradEnergyTasks.cpp`, 532 lines, read by the director): a forward pass, the cotangent of
  every unit, one fill and factorisation of the level matrix per task, sums over 64 interleaved chunks added in a
  fixed order, row blocks added in task order. That fixed order is what makes the fit's numbers independent of the
  thread count. Deleted: `BuildGradItems` with its grouped and ungrouped routes, `GradOnePoint`, `GroupRows`,
  `ForEachRow`, `TBarAdjointMany`, `ScopedPointMappingOff`, the switch `AZURE_NO_GROUPED_GRADIENT`.
- The 5-thread slowdown of branch 1 (director's finding, 7 % on the T2 fit): the agent measured its cause (the
  tasks rebuilt serially at 156 of 196 passes, 1.15 s of 15.6 s) and removed it (0.37 s). The director has not
  repeated the 5-thread timing: the machine is limited to 12 threads this night and other runs share it. Owed.
- Against the base at 1 thread the fits' end points move by rounding along the path (T2 2.5e-9 in a parameter, 13N
  6.8e-11); `natural_neon` stopped after 30 calls far from a minimum differs by 5.9e-5 (the paths agree to 1e-15 at
  the first call). Accepted: a fit's path is not a regression quantity; the outputs of calculations are identical.
- No defect found in the reading. One thing for extensibility: the parameter cotangents name the kinds of data
  parameter a third time, as characters (`ParamCotangent::kind` 'a', 'd', 't'); brief AS is told to enter new
  apparatus parameters there by the parameter table, so that the next one adds no case.
- The point mapping (`MapData`, `MapLeaves`, `SharesCalculation`) is still there and no longer saves anything: with
  it switched off every output is identical except a look-up count in two `run.log` files; time equal, memory +3.5 %
  on T2. Its deletion is the first milestone of the node-table branch.
- For the merge (director): the record sizes of `tests/consistency/file_sizes.tsv` for files two branches changed
  were set by hand (`--tighten-sizes` only lowers). The agent's speed series and a last documentation commit had not
  arrived at the time of the merge; they are merged when they do.

### Apparatus effects as kernel functions (`feature/upgrade-5-kernels`, Phase 5 branch 6) — director's review

- What it does: the folds of a point (target integration with its four target kinds, Gaussian resolution, beam
  profile, UDR, Doppler broadening) left `Unit` and are plain functions in `experiment/src/Kernel*.cpp`, listed in one
  table (`experiment/src/Kernels.cpp`, `Kernels()`); the engine looks the effect's kernel up and calls it
  (`src/PointFold.cpp`, 112 lines). `src/TargetFold.cpp` (502 lines) is deleted. Value and gradient weights come from
  the same template (on `double` and on `LinearForm`). No number moved (agent's run of every project; the director's
  batch verification repeats it).
- **The extension count the PI asked for.** The agent added a rectangular energy resolution
  (`energy_resolution rectangular width 20 keV` on a detector) by its own recipe: 1 new kernel file (69 lines), 1 line
  in the table, 1 line in `experiment/CMakeLists.txt`, tests, schema and documentation; **0 files under `src/` and
  `include/`** (director read the commit's file list: 14 files, +422/−13), 4 minutes to a passing test. The finite-layer
  yield a week earlier had touched 11 engine files.
- **The limits of that count (director's reading of the code; say them to the owner):**
  1. It holds for a *stand-alone resolution kernel with scalar parameters*, selected by `detector.energy_resolution
     <form>`. A kernel that needs a new datum of the point (as the Doppler kernel needs the temperature, the photon
     kernels Q and mass) or a non-scalar parameter edits `KernelPoint`/`KernelParameters` in `Kernels.h` and the
     engine function that fills them (`TargetEffect::KernelParametersFor`). `KernelParameters` still has named
     members for every older kernel; only the new table kernels use the generic `values`.
  2. **Kernels do not compose yet.** A table kernel on a data set that also has a target, a beam spread, a beam
     profile or a UDR is refused (`src/TargetEffect.cpp`, `SelectKernels`: "the <kernel> kernel cannot be combined
     with a target, a beam spread, a beam profile or a UDR"). So the rectangular resolution cannot be used on a
     target-integrated yield today. The target kernel takes the beam spread and straggling as its own parameters;
     the general chain (target ∘ resolution ∘ broadening on one mesh) is the grid branch (branch 7).
  3. The kernel's parameter row says `fittable = true`, but a kernel parameter is not a fit parameter yet: there is
     no generic "kernel parameter" kind that reads the table. Follow-up after the stopping-wire merge (it adds the
     stopping-scale kind in the same files): one kind, read from the table rows, so that a new kernel's parameters
     are fittable with no further code. Until then the row's flag is a declaration, not a feature.
  4. Of the derivative column of `energy_grid.md` §3.3 only the Doppler temperature and the rectangular width are
     kernel functions; the others (areal density by central differences, etc.) are still engine code.
- New refusal users can meet: two resolution kernels on one data set ("one resolution kernel per data set: X and Y
  cannot be combined"); before, one of them was silently not applied. Accepted; listed under "What users notice".
- Director's own runs (`work/b-kernels`, 1 thread): `ctest -R kernels` passes; `tests/rectangular_resolution/check.sh`
  passes (folded points within 8.6e-5 of the window average of the bare output; the tolerance of the check is 1e-3).
- For the merge: 14 lines of `src/AZUREGrad.cpp` (`FoldPoint`/`FoldWeights`/`FoldTemperatureDerivative` in
  `GradTargetEffectAdjoint` and `ThmPointModelGradient`) meet the gradient branch's rewrite; the stopping callback in
  `KernelParametersFor` must call the stopping-wire branch's accessor.

### The engine takes its stopping from `stopping/` (`feature/upgrade-5-stopping-wire`, Phase 5 branch 5) — director's review

- What it adds: a target says `stopping evaluation` (SRIM-2013 for its `composition`, `phase`, `compound_correction`
  and the projectile of the data set's entrance pair) or `stopping table "file" <units>`; `src/TargetStopping.cpp` is
  the one place a stopping is computed, per active nucleus (the material's stopping per atom over the active
  element's atom fraction), and `TargetEffect::StoppingAt` applies the stopping scale once. The scale is a fit
  parameter (`stopping_scale.<material>.<projectile>`, one per material and projectile, fixed at 1 by default, prior
  the record's or SRIM's stated accuracy when fitted: decision D6). The run's log, `--validate` and `--describe`
  print the provenance and SRIM's citation. `pyazr.stopping` binds the library. ERYA is removed (decision D7:
  `external/erya` only; the director checked that no other path under `external/` is touched).
- Outputs of the present projects: identical (agent's run; the director's batch verification repeats it).
- **Director's own check on a material the agent did not choose** (`work/director/stopcheck/`: Si₃N₄, active ¹⁵N,
  atom fraction 4/7; the narrow resonance of the test; own script with LunaSRIM 1.0.0 for the stopping at the
  resonance and the exit energy, Breit–Wigner integral over the layer): thick layer (45.6 keV), four energies:
  engine/own − 1 between −8.0e-5 and +8.6e-5 (the line-shape and quadrature level of `layer_yield`); a 5.63 keV
  layer scanned every 10 eV over both edges, 127 points: largest |engine − own| 2.0e-4 of the plateau, at the back
  edge, which is 0.03 eV in the exit energy (6e-6 of the layer's loss; the two integrators agree to 5e-6).
- **What a new kind of fit parameter costs (the count the PI asked for): 20 existing files** (parameter table, data
  set, working copy, hand-over, objective, gradient map, finite-difference Jacobian, API, binding, labels, outputs,
  documents), most by 1–20 lines, and one new file. This is not "a few named places". It is the argument for the
  follow-up noted in the kernels review: one generic kind for kernel and target parameters read from the table.
- Limits to know: the scale's derivative is by central differences on the whole χ² (as the density); a fitted scale
  with fold windows or the automatic tolerance is refused; TRIM's straggling is not offered for an evaluated target;
  a scale on an extrapolation's data set is untested; `--validate` cannot check the energy range (refused at set-up);
  the old integrators do not call the library's integrals (brief step 6, not attempted);
  `TargetStopping::LibraryPerActiveNucleus` returns NaN when the data are missing (the set-up refuses first;
  tested with an empty `AZURE2_SRIM_DATA_DIR`).
- Found by the director: `pyazr_stopping` (and five older pyazr tests) were registered after the block of
  `tests/pyazr/CMakeLists.txt` that hands the scratch build's engine module to the tests, so under ctest they were
  skipped unless the caller exported `PYAZR_ENGINE_DIR`. The block is last again (rehearsal commit).
- Accepted: the one line in `src/AZUREGrad.cpp` (`AppendStoppingScaleParams` in `BuildParamIndexMap`); ten printed
  digits of `azure2-stopping`.
- For the merge with the kernels branch (director): the kernels' stopping callback of `KernelParametersFor` (the
  support of a thick target with Bohr straggling) calls `StoppingAt(E, config, nominal scale)`, as the stopping
  branch's `GridExtent` did; `GridExtent`, `GaussianKernel`, the straggling functions stay deleted (they are kernel
  functions). With the pyazr branch: pyazr's version (2.33.0 here, 3.0.0.dev there).
- **For the owner:** SRIM's terms (free for non-commercial use) against a GPL repository that ships its two data
  files: unchanged, the owner's decision.


## 2026-11-11 (agent, Phase 5 branch 2) — the gradient over the energy tasks (`feature/upgrade-5-energy-tasks-gradient`)

**Done.**
- The 5-thread regression of branch 1, measured: a profile stage "energy tasks: set-up" and a counter "energy tasks
  rebuilt" (AZURE_PROFILE) show a T2 fit rebuilding the tasks at 156 of 196 passes (every shift change rebuilds the
  plan), serially, 1.15 s of 15.6 s wall at 5 threads. The tasks are now ranges of flat arrays (`EnergyTaskList`), keys
  and sort in parallel: 0.37 s (5 threads), 0.50 s over 234 passes (3 threads).
- The analytic gradient, the least-squares rows and the model gradients run over the energy tasks
  (`src/GradEnergyTasks.cpp`): the point pass, the cotangents of every unit (`PointCotangents`;
`GradTargetEffectAdjoint`
  hands out sub-point cotangents), one fill and factorisation per task, T-matrix cotangents summed per group, bilinear
  terms merged per task (rows: per row or by the task's basis). The forward of the adjoint is recomputed per task.
  Sums in a thread-independent order (64 interleaved chunks; rows task by task; fitBarFn serial in data order).
- Deleted: `BuildGradItems` and the grouped/ungrouped routes, `GradOnePoint`, `GroupRows`, `ForEachRow`,
  `AZURE_NO_GROUPED_GRADIENT`, `ScopedPointMappingOff` (its two other users run the point pass), `TBarAdjointMany`,
  `AMatrixFunc::PointAdjoint`. A crash found by the speed series fixed: the row pass with no free parameter
  (thick_target_resonance `--use-lm`) divided by the column count; `tests/upgrade/energy_task_gradient_test.cpp`
  (ctest `upgrade_energy_task_gradient`, shown failing with a segfault without the fix) also checks that the scalar
  gradient, the residual Jacobian and the model gradients of T2 are the same bits at 1 and 3 threads.
  `tests/angle_mapping/check.sh` (d) now requires the gradient at 1 thread to equal the one at the check's thread count
  bit for bit (the base fails it: 3, 4, 2 of 21 parameters).

**Verified.** On 2d20a50b (`work/egrad/AZURE2-final2`) unless said:
- all projects x 3 bases identical to `tab-trial21a` (766/739/766; a69f6854, 6 jobs, before the 3-thread budget);
- `AZURE_PLAN_CHECK=1` on the tree: 265 evaluations, 3,653,648 units, 0 differ;
- `ctest -LE slow` 201/201 (-j1, OMP 2); `tests/run_tests.sh` 148 passed, 0 failed (OMP 2);
- converted-tree gate (`REF=tab-trial8`, JOBS=3, tag `gegrad`): 483 (456) identical, 38 classified, 0 differ, 0 missing;
- calculations at 1 and 3 threads identical: T1, T2, T3, 13N, natural_neon, thick target, energy_shift_sqrt, shift_fn,
  18O_p_a_thm, target_effect_ranges, angle_mapping;
- fits at 1, 2, 3 threads and two repeated 3-thread runs identical (param.par, chiSquared.out, param.sav) for T2
  `--use-gradient`, T2 `--use-lm`, 13N `--use-gradient`; the base's two 3-thread T2 `--use-gradient` runs differ
  (chiSquared.out, param.sav); its T2 `--use-lm` and 13N runs agree with each other;
- the gradient at the start (AZURE_GRAD_CHECK, 1 thread, `work/egrad/gstart`): 13 fit projects with checked
  components (18O_p_a_thm and target_density print none; park_formalism, doppler_fgm, transmission_yield and
  natural_target were not found by the script: their .azr lie in subfolders), largest change against the base 3.3e-13
  of the largest component (7Li_p_ay), agreement with finite differences unchanged;
- the gradient at the same parameter points, base against this branch (pyazr, 1 thread, every component,
  `work/egrad/twopt`): natural_neon at the start 3.9e-14 of the largest component (per component 6.2e-13), at the
  base's 30-call end point 1.9e-14 (1.9e-12); rows 2 J^T r <= 1.3e-16; the screening fit 7.4e-14 at its end point;
  T2 <= 1.3e-14. So the fit end points that move (natural_neon 30 calls 5.9e-5 in a parameter, EDM 3.4e3; the
  screening fit 200 calls 0.13, EDM 27) are rounding amplified on unconverged paths (director's class);
- fits at 1 thread against the base: T1, T3, T2 `--use-lm`, angle_mapping, energy_shift_sqrt, sum_free, 18O THM
  identical; T2 2.5e-9 in a parameter (chi2 1.5e-13); 13N 6.8e-11 (2.6e-14);
- speed: `energy_grid.md` section 8 (13N gradient fit 51 -> 27 s, natural_neon 220 -> 106 s; T2 at 3 threads not slower
  than b-trial19). After the fix commit: build, quick set, the new test, thick-target fits.

**Not done / open.**
- The 5-thread acceptance time (budget 3 threads that night): the director times it.
- ThreadSanitizer: not possible here (no clang; `scripts/sanitize_upgrade.sh` needs clang with LLVM's libomp and
  libarcher, a GCC/libgomp build reports every OpenMP barrier as a race).
- The point mapping stays (director's decision): its deletion is a follow-up after branches 5 and 6 merge, written out
  file by file in `energy_grid.md` section 8 with the evidence (tag `tab-egrad-nomap`).
- THM locality (`Unit::GetLoElement` cold rows) not attempted; most of a T2 gradient is the finite-difference part for
  shifts and THM.
- Not edited (other branches' or historical files) that still name deleted functions: `experiment/README.md`,
  `experiment/tests/expression_test.cpp` (GradOnePoint, as history), `MERGE_NOTES.md`, `docs/upgrade/edata_split.md`.

## 2026-11-11 (agent, Phase 5 branch 5) — the engine uses `stopping/` (`feature/upgrade-5-stopping-wire`)

**Done** (`energy_grid.md` §6 and §8 "Branch 5, done", `stopping.md` §14).
- A target says `stopping evaluation` (SRIM-2013 of its `composition`, `phase`, `compound_correction`, for the light
  particle of the entrance pair of each data set naming it) or `stopping table "file" <energy unit> <stopping unit>`
  beside `stopping formula`.  `TargetStopping` (`include/TargetStopping.h`, `src/TargetStopping.cpp`) is the one place
  the stopping is computed: formula | evaluation | table, per active nucleus ε_atom/f_active (the convention is stated
  in `SetupSpecs.h` `StoppingSpec`).  `TargetEffect::StoppingAt(E, scale)` multiplies the scale once and is what every
  reader asks (the sub-points' stopping in the set-up loops and `ReplaceSubPointGrids`, `TargetThickness`,
  `MeanEnergyLoss`, `GridExtent`, the layer refusal).  Build refusals `E_STOPPING_ACTIVE`, `E_STOPPING_PROJECTILE`,
  `E_STOPPING_SCALE_CONFLICT`; set-up refusals `E_STOPPING_DATA`, `E_STOPPING_RANGE`.  The log, `--validate` and
  `--describe` print the library's provenance (SRIM cited) once per target and projectile.
- The stopping scale (D6): parameter kind 12 `stopping_scale.<material>.<projectile>`, one per (material, projectile),
  last in the vector; `fit` without a prior takes SRIM's stated accuracy (printed); its prior is a `TargetPriorTerm`;
  the hand-over (`src/StoppingScale.cpp`) re-sets the sub-points' stopping and the points' thickness; derivative by
  central differences of the whole χ², as the density.
- `tests/stopping_yield` (README): plateau of a thick TiN layer against LunaSRIM's ε (8.3e-5, the line shape), the
  width of a layer's step against LunaSRIM's energy loss (2.4e-7), table of the library's values against the
  evaluation (2.1e-10), a fit of two layers' N and their one scale (pulls −0.66, +0.61, +0.46), the refusals;
  ctest `stopping_convention` (Ta₂O₅ with ¹⁸O, TiN with ¹⁴N by hand).
- `pyazr.stopping` (`api/src/pybind_stopping.cpp`, `pyazr/stopping.py`), ctest `pyazr_stopping` against
  `srim_reference.csv`; pyazr 2.33.0.
- ERYA removed: `external/erya`, `USE_ERYA`, `SRIM2013.xml`, `tests/reference/compound_stopping_reference_test.cpp`
  (its check is `stopping_convention` now).

**Verified** (binary `work/stopwire/AZURE2-final`, HEAD a7b9cac2): all projects × 3 bases against `tab-trial21a`
(tag `tab-stopwire1`, 6 jobs): 766 / 739 / 766 identical, the 12 files per basis of `tests/stopping_yield` only-new;
`ctest -LE slow` 202/202 (`-j3`, OMP 2); `tests/run_tests.sh` (OMP 2) 149 passed, 0 failed; converted-tree gate
(`REF=tab-trial8`, tag `gstopwire1`, 3 jobs) 483 (456) identical, 38 classified, 0 differ, 0 missing; calculations at 1
and 3 threads identical (stopping_yield's two models, layer_yield's two); a fresh build of the tree without
`stopping/data/srim/` (0 warnings): with `AZURE2_SRIM_DATA_DIR` the outputs equal the compiled-in build's, without it
`E_STOPPING_DATA` and exit 1.  Scratch: `work/stopwire/`.

**Open.** §6.4 step 2 (the old integrators on the module's integrals) not done; TRIM straggling for an evaluated
target (the `straggling` field still gives none, Bohr or a coefficient); the scale's analytic derivative (kernels
branch); a scale in an extrapolation's data set is not tested.

## 2026-11-11 (agent, Phase 5 branch 6) — apparatus folds as kernel functions in one table (`feature/upgrade-5-kernels`)

**Done.**
- 14c0cc26: the folds of `src/TargetFold.cpp` (502 lines, deleted) are kernels of `experiment/`
  (`azure2/experiment/Kernels.h`, table `Kernels()` in `experiment/src/Kernels.cpp`): `KernelTarget.cpp`,
  `KernelGaussian.cpp`, `KernelBeamProfile.cpp`, `KernelUdr.cpp`, `KernelDoppler.cpp` (broadening), quadrature in
  `KernelQuadrature.h`; `LinearForm.h` (was `include/TELinearForm.h`), `LayerExit.h`, `UdrKernel.h` moved to
  `experiment/`. Engine fold: `src/PointFold.cpp` (`FoldPoint`, `FoldValues`, `FoldWeights`,
  `FoldTemperatureDerivative`, `FoldBroadeningWeights`), one indirect call per point. `TargetEffect` keeps the
  parameters, chooses its kernels once (`Kernels()`, `KernelParametersFor`), refuses two resolution kernels on one
effect
  (`HasKernelError`). `TargetSetup.cpp`: both copies of the extent code are one `EffectSupport` call. Unit lost 10 fold
  members (717 -> 671 lines), including the hidden-state `DopplerAppliedScope` (now `Broadening::kAlreadyApplied`).
  `AZUREGrad.cpp`: 13 call lines + 1 include line changed (for the merge with the gradient branch).
- d06be05f: docs; CLAUDE.md recipe "A new apparatus kernel".
- 1d3a35ac: rectangular resolution (`KernelRectangular.cpp`) by the recipe alone, no engine file; test project
  `tests/rectangular_resolution/`; f7dab67c adds its lines to `tests/upgrade/predictions.txt` (step added to the
recipe).

**Verified** (scratch `work/kernels/`): all projects x 3 bases vs `tab-trial21a` (tag `tab-kernels-final`, JOBS 6):
766 / 739 / 766 identical, only-new = the new project's run.log; ctest -LE slow 201/202 at -j1, the failure
(`upgrade_predictions`) fixed in f7dab67c and re-run alone; `tests/run_tests.sh` (OMP 2) 149 passed, 0 failed; 1 vs 4
threads identical on 7 folded projects x 3 bases;
fits at 1 thread byte-identical to `b-trial21` (T2 --use-gradient, --use-lm, 13N, layer_yield/fit.azr); gate `gkernels`
483 (456) identical, 38 classified, 0 differ, 0 missing; user-space instructions new/base 13N 0.9999, thick target
0.9989, T2 0.9876 (cycles too noisy on the shared machine). Note: the quick comparison needs
`AZURE2_COULOMB_TABLE` set, else three projects differ from the reference in the 10th digit (base and new alike).

**Open.** Kernel parameters are not fit parameters (no generic parameter kind); the derivative columns of
`energy_grid.md` §3.3 are not yet kernel functions (except Doppler's d/dT); comments in `AZUREGrad.cpp` still name
`Unit::TargetFallbackImpl`/`TargetEffectWeights` (left for the gradient branch's merge).

## 2026-11-11 (director, 16:39) — merged: seven branches (finite-layer yield, pathway table 2, one solve per energy, SRIM stopping, lab frame, THM projects in format 3, format-3 gaps)

**State.** `feature/upgrade-4-data` holds the seven branches, merged in this order; the tree is the one the director
verified as a rehearsal (`wt-trial21`, then `wt-trial22`). All 183 convertible projects are format 3. In work (three
agents, from the head before the last two merges): the engine uses `stopping/` (`feature/upgrade-5-stopping-wire`),
the gradient per compound energy (`feature/upgrade-5-energy-tasks-gradient`), the folds as kernel functions
(`feature/upgrade-5-kernels`). Launched after this entry: pyazr off the old file class (brief AQ), the audit of the
checks' text edits (brief AR).

**Director's verification** (own runs, `work/director/full_verify2.sh`; logs `work/director/v2_trial21a.*`,
`v2_trial22.*`):
- First five branches against the outputs before them (`tab-trial19`): every file identical except the rounding
  class of the lab-frame branch and the two files of `sum_lab_frame`; new projects `layer_yield`, `ec_exit_position`.
- THM and gaps on top (`tab-trial22` against `tab-trial21a`), 766 / 739 / 766 files (Brune / formal / Park):
  721 / 696 / 721 identical. What differs, all expected: `param.par` of the THM projects (parameter names by record:
  `o17_guardo.shift.offset` for `segment_1_energy_shift`), 17O's `normalizations.out` (a THM data set's angle range
  printed 0–180 instead of 0–0), and 39 `run.log` per basis (the refusal line of a project whose data files a check
  writes; `M_thick` no longer announced as an old file; the THM projects' names).
- `ctest -LE slow` 201/201, `ctest -L slow` 10/10, `tests/run_tests.sh` 149 passed and 1 failed, repaired below;
  1 against 6 threads identical; spec equivalence 202 files, 184 converted and equal, 18 not converted (12 written in
  format 3 from the start, 6 old on purpose); the converted-tree gate 483 (456 formal) identical, 38 classified, 0
  differ, 0 missing.
- **The one failure was a meeting of two branches**: both refuse a THM experiment that no `<thm>` line defines, the
  gaps branch when the project is read (`E_THM_EXPERIMENT_UNKNOWN`, with record and field), the THM branch at run
  time. The reading-time one now comes first, so `tests/thm_experiment` (undefined) expected a message that can no
  longer appear. The director changed the check to expect the reading-time diagnostic (it passes); uniting the two
  refusals and the twice-written experiment-name grammar is in brief AR.
- Not run: Sphinx, ThreadSanitizer (owed by the gradient branch), the GUI.

**Open, for the owner or the PI:** SRIM's terms against the GPL repository (below); fits at several threads are not
reproducible run to run (1e-10 in χ² terms, also before these branches; cause not shown); the sub-point quadrature
over a narrow resonance on a thick layer is good to 1.5e-4 (the grid branch's test case).

### Finite-layer yield (`feature/upgrade-5-layer-yield`) — director's review

- What it adds: a data set with `quantity yield` on a target of `kind layer` is Y(E₀) = 10⁻²⁴ ∫ ⟨σ⟩/ε dE from the
  exit energy to the beam energy, the exit energy being where the layer's areal density is used up
  (`include/LayerExit.h`), on the point's own sub-points. Refusals at set-up: `E_LAYER_STOPPING`,
  `E_LAYER_EXIT_BELOW_MODEL`.
- **What it cost (the extensibility count the PI asked for):** 11 existing engine files, about 100 lines, and 5 new
  files. The yield is a branch inside the target fold (`Unit::IntegrateKernelsImpl`), not yet a kernel of its own:
  that is the kernels branch. A new quantity is therefore not yet "a few named places"; say so to the owner.
- **Director's independent check** (`work/director/layercheck/`, own trapezoid integral of σ/ε on a fine grid of
  the bare cross section written by the engine, 38,197 points, analytic exit energy):
  thin layer 1–3e-7; scan over a narrow resonance, largest |ΔY| 1.2e-5 of the plateau (4e-4 relative on the rising
  edge, median 3e-6); **thick layer (N = 1.1e19 atoms/cm², 60 quadrature points) the engine is 1.5e-4 above the
  integral**, and it converges slowly with the number of points (120: 1.5e-4, 240: 1.4e-4, 480: 1.07e-4). This is
  the accuracy of the present sub-point quadrature over a narrow resonance; the infinite-target yield has it too.
  The agent's 3.2e-6 for the thick layer compares the engine with its own infinite-target yield: a consistency
  number, not an accuracy. The grid branch (decision D2) gets this case as a test with the director's integral as
  the reference.
- Limits to know: ∂Y/∂N by central differences; a layer shared by a cross-section data set and a yield data set with
  free screening or temperature is refused (`E_FOLD_CONFLICT`); `convert/tests/spec_equivalence_test.cpp` does not
  compare `layerYield`; a stale comment in `tests/units_scales/models/yield_story.azr`.

### Pathway table, second branch (`feature/upgrade-4-pathway-table-2`) — director's review

- Defect 1 of `cnuc_split.md` §11.9 is fixed: external capture finds the exit to the photon pair by the pair's
  number (`include/ExternalCaptureExit.h`). The base read past the end of the exit list when a pair without a
  pathway from the entrance was numbered below the photon pair: an empty `intEC.dat` and a segmentation fault
  (reproduced by the director with `work/b-trial19`). New test directory `tests/ec_exit_position`.
- Director's own run (`work/director/pathcheck/`): the two models of that directory with an energy shift added to
  two capture data sets (the shifted amplitude takes the same look-up), with the energy tables and with
  `--no-ec-table`: every output file equal in the two orders of the pairs.
- The memo of the shifted amplitude keeps (energy, amplitude) only: the final level's energy has one writer at
  set-up (`CNuc::SetUpECFinalBlocks`), so it cannot depend on the history of a copy; `history_test` project `ec`
  holds that.
- The pathway trees are deleted (`Decay`, `KGroup`, `MGroup`, `ECMGroup`, `KLGroup`, `Interference`: 12 files, 971
  lines; `SortPathways`, `CalcAngularDists`). `CNuc.cpp` 1,774 → 1,548 lines.
- Left: a stale comment in `src/ECIntegral.cpp:438` that names `CNuc::SortPathways` (outside the branch's scope);
  the set-up amplitudes could be shared between copies instead of copied (would move fit bits: not done).
- `cli_help` of `tests/run_tests.sh` failed once under load for two agents and passes alone (cause not shown).

### One solve per compound energy (`feature/upgrade-5-energy-tasks`, Phase 5 branch 1) — director's review

- What it does: the point pass groups the units (points and sub-points) by (nucleus, route, compound energy bits)
  and fills, factorises and solves the level matrix once per group (`src/EnergyPass.cpp`: `SolveAtOneEnergy`);
  units whose external-capture amplitudes differ get their own T-matrix inside the group. `AZURE_PLAN_CHECK=1`
  recomputes every unit the old way and compares: 3,727,625 units over the tree, 0 differ. Outputs of every project
  are byte-identical (director's run of the whole tree on the merged batch).
- Measured by the agent (point pass per evaluation, 1 thread): T2 106 → 45 ms, 13N 317 → 132 ms, thick-target
  project 1,102 → 264 ms, `natural_neon` 3,207 → 1,718 ms, THM unchanged.
- **Director's own timing, the caveat:** the 200-call T2 fit (MIGRAD) at 5 threads, alternated three times, is about
  7 % slower in wall time (base 15.7 / 14.3 / 15.6 s, new 16.2 / 16.3 / 16.4 s; user time equal, 58 s); at 1 thread
  equal (44.4 against 43.6 s). The analytic gradient still solves per host point, so a gradient fit gains nothing
  yet, and a serial section grew. Accepted for the merge; removing it is the first item and an acceptance test of
  the gradient branch (brief AO).
- **Kept on purpose:** the point mapping (`MapData`, `MapLeaves`, `SharesCalculation`, `ScopedPointMappingOff`) and
  `BuildGradItems`; they go with the gradient branch. `tests/component_mapping/check.sh` part (b) asserted that the
  mapping saves level-matrix solutions; that is no longer true (2,769 with and without it), so the director changed
  (b) to assert the two counts are equal.
- **Found on the way:** a fit at several threads is not reproducible run to run, also in the base build (two runs
  of `b-trial19` at 5 threads differ in 8 output files, χ² terms at 1e-10). Fits are compared at 1 thread only.
  Whether the cause is the order of an OpenMP reduction is not shown; the gradient branch owes a ThreadSanitizer run.
- Not run by the agent: ThreadSanitizer, Sphinx, `ctest -L slow`.

### SRIM-2013 stopping through LunaSRIM (`feature/upgrade-5-stopping-2`) — director's review

- PI's decision (2026-11-11): the stopping powers are SRIM's as implemented in the PI's LunaSRIM, SRIM cited.
  `stopping/` now holds a C++ port (`Electronic.cpp`, `Nuclear.cpp`), SRIM's two coefficient files in
  `stopping/data/srim/` (compiled in; `AZURE2_SRIM_DATA_DIR` reads them from a directory instead), the phase of a
  material (solid | gas), the compound correction, TRIM's straggling beside Bohr's. Nothing in the engine calls it
  yet (branch `stopping-wire`).
- **Director's own comparison** (`work/director/srimcheck/own.py`, the installed lunasrim 1.0.0 against the
  `azure2-stopping` tool, 9 ion–material cases the agent did not choose): worst relative difference 5.8e-6, which is
  the tool's printing precision; layers: p 0.3 MeV in Ta₂O₅ 1e18 atoms/cm² exit 283.05601 against 283.05604 keV,
  straggling 2.687 against 2.6874 keV; ¹⁴N 1 MeV in Au 3e17 atoms/cm² 901.22377 against 901.22407 keV, 16.88 against
  16.881 keV. Regenerating the committed reference table with the committed generator reproduces it byte for byte.
  The two data files have the SHA-256 of LunaSRIM's.
- Decision accepted by the director: the compound correction is applied once (as SRIM's table mode, TRIM and
  LunaSRIM's `energy_loss`); LunaSRIM's `compound_stopping` applies it a second time below 10 MeV (SR module list
  mode), documented in LunaSRIM's own physics note.
- **To know:** the uncertainties (H 4.0 %, He 3.9 %, Li 4.8 %, heavier 5.8 %) are SRIM's stated mean deviations from
  the measurements, used as one standard deviation: an interpretation, not a statement by SRIM; nuclear 10 % is
  provisional. Below about 10 keV/u SRIM and PSTAR/ASTAR differ by tens of per cent (p in He −84 % at 2 keV, α in H
  +60 % at 1 keV): `stopping.md` §9.2.
- **For the owner, not decided here:** SRIM's terms, as LunaSRIM's README states them, allow free copying for
  personal, non-commercial use and exclude commercial use without a licence from SRIM.org. AZURE2's public repository
  is GPL. Whether the public repository ships `SCOEF03.DAT` and `SNUC03.dat` is the owner's decision; the library
  works without them when `AZURE2_SRIM_DATA_DIR` points to a copy (`stopping/data/srim/README.md`).

### Data compared in their own frame (`feature/upgrade-5-lab-frame`, Phase 5 branch 3) — director's review

- A laboratory-frame differential data point is no longer converted to the centre of mass at set-up; the model is
  divided by the Jacobian at comparison (`DataResidual`). `FrameRow` loses the stored c.m. value and error.
- Outputs: identical except for last-digit rounding in 23 projects (class `lab-frame-rounding`: 1e-14 relative or
  one unit of the last printed digit) and the `run.log` of `sum_lab_frame`, whose χ² is itself rounding noise
  (4e-17; it moves by 7e-23). The class was not loosened for it: the two files are a listed, explained difference.
  Fit end points move by at most 9.1e-6 relative (`angle_mapping`, one width).
- Open, not blocking: a Jacobian of exactly 0 (the kinematic limit) would give an infinite pull where the old code
  dropped the point: should be a set-up refusal; the pyazr docstrings of `residuals` and `residual_jacobian` still
  say c.m. error; THM lab-differential data are untested.

### The THM projects in format 3 (`feature/upgrade-4-tree-format3-5`, Phase 4 group 7) — director's review

- All 183 convertible projects of the tree are format 3. A THM data set names its experiment on its record
  (`thm <name>`, `thm_weight`); `<thm>` keeps the options and the `experiment[...]` lines and is built empty when
  records name experiments and the file has no section (`src/ConfigThm.cpp`: `Format3ThmSection`).
- Engine diff read by the director: `ConfigThm.cpp`, `ThmBlockEdits.cpp/.h` (`ThmExperimentsWithoutMembers`; idle
  experiment lines skipped), `convert/src/OldProject.cpp` (five refusals of old member/weight forms). 20 checks
  rewritten on records (`tests/lib/thm_records.sh`, `tests/pyazr/thm_records.py`).
- **Known defects, given to the audit branch (brief AR):** a format-3 `<thm>` still accepts the old `segments=`
  list; an experiment no data set names is refused in one place with "segments= is required" and warned about in
  another; **checks whose text edits no longer hit anything in a converted project can pass without testing**
  (`thm_fit` was one, found and repaired by the agent; the others are not audited yet).
- The 13 pyazr tests of the old file class's THM API still read the originals (brief AQ moves them).

### What format 3 could not yet say (`feature/upgrade-4-format3-gaps`) — director's review

- A refused format-3 project ends with exit 1 and its refusal as the last line (was "Could not open …");
  `tests/cli_route_exit`. New diagnostics `E_THM_EXPERIMENT_UNKNOWN`, `E_PARAMETER_FILE_NO_MATCH`,
  `W_PARAMETER_FILE_NAMES` (86 codes in all).
- New in the format: `sigma_total` as a sum of elastic terms on a target with isotopes; `total_optical` on records;
  a detector's `energy_resolution` on β-delayed spectra; with pair keys, the keys decide elasticity. Converter:
  a Doppler-only `<targetInt>` on neutron samples becomes the sample's `doppler` field; `M_thick.azr` converted.
- 24 old-format generators deleted; each README says how the models were made and which commit holds the script.
- `AZURE_PROFILE` prints three reports (set-up, calculation, final evaluation); `tests/profile_report`.
- **Merge with the THM branch** (director): `experiment/src/BuildEffects.cpp` keeps both flags of an effect
  (`layerYield`, `decay`) and both suffixes of its signature. The undefined-experiment refusal now exists twice
  (reading time in `ProjectReader.cpp`, run time in `ThmBlockEdits.cpp`), the experiment-name grammar twice: brief
  AR unites them.
- Left (brief AR): a set-up refusal other than a missing data file still ends with "Could not fill data object from
  file."; the converter writes an inactive natural-target line without its sample target. Left (brief AQ):
  `load_sav` does not check parameter names. Noted: the adjoint gradient's rounding depends on the order of the
  fields of a record (1e-14).
- Not run by the agent: Sphinx, the full `scripts/check.sh`; `tests/run_tests.sh` was not repeated after the last
  fix (the director's batch verification does).

## 2026-11-11 (agent, format-3 gaps) — what format 3 could not say, and the reading path's defects (`feature/upgrade-4-format3-gaps`)

**Done** (`AL_format3_gaps.md`; one commit per item):
- A refused format-3 project exits 1 with its diagnostics and the refusal line last (`Config::kConfigRefused`,
  `ReportFormat3Refused`), also when a data run cannot read a data file (`E_DATA_FILE`, `WorkingCopy::kFillRefused`);
  a missing project file keeps "Could not open", exit 255. `cli_route_exit`.
- A `thm <name>` that names no `experiment[<name>]` of `<thm>` is refused when the project is read
  (`E_THM_EXPERIMENT_UNKNOWN`, `ProjectReader.cpp`).
- A parameter file given as a run's start is checked by name (`experiment/ParameterFileNames.h`): no name of the run
  refuses it (`E_PARAMETER_FILE_NO_MATCH`, exit 1), a partial match warns with both sides (`W_PARAMETER_FILE_NAMES`);
  the project's other run's names and fixed data-set/extrapolation parameters are not reported. `param_file_names`.
- `sigma_total` takes a sum of elastic terms on a target with isotopes (the natural total); `total_optical` runs on
  records (outputs equal to the old models' but param.par names and sample headers).
- A detector's `energy_resolution` on a beta-delayed spectrum is the Gaussian fold of the measured energy;
  `beta_delayed_alpha` (e) on records; the converter writes that detector for the old line.
- With `entrance_pair` / `exit_pair` the keys decide whether a reaction is elastic (`PairsDecideElastic`), not the
label.
- Converter: a Doppler-only `<targetInt>` line on neutron-sample lines becomes the samples' `doppler` field; a target
  shared by lines of a sample and of none is copied for the sample (M_thick's four questions); `M_thick.azr` is format
3.
- The 24 generators that wrote old files are deleted; each README says how the models were made.
- `AZURE_PROFILE` reports the set-up apart and labels a calculation's report.

**Verified how** (on the last code commit; tag `tab-f3gaps`, the new reference): all projects x 3 bases against
`tab-trial19`: the same files, every output file identical; run.log differs in 33 projects per basis (data files not
in the folder: the refusal line of item 1 replaces "Could not fill data object from file.") and in M_thick (its
old-format note gone). Converted-tree gate (`REF=tab-trial8`, `gf3gaps`): 483 (456) identical, 38 classified, 0
differ, 0 missing. Spec equivalence: 198 files, 184 converted (M_thick added), 184 equal, 0 differ, 14 not converted.
`ctest -LE slow` 186/186, `OMP_NUM_THREADS=2 ctest -L slow` 10/10, `tests/run_tests.sh` 147 passed, 1 failed
(model_data: the thm check refused inactive data sets; fixed in the last commit, model_data rerun: passed). 1 vs 4
threads: T1, T1_start, T2, T3, 13N, natural_neon menu 1 byte-identical.

**Open**: pyazr's `load_sav` does not check parameter-file names (pyazr/*.py not touched); the converter writes an
inactive old natural line without its sample target (it would change meaning when switched on); the adjoint gradient
differs from the old path's at 1e-14 relative and moves with the order of a data set's target/beam fields (FD and
Jacobian path identical); `18O_p_a_thm/make_inputs.py` (group 7); `energy_grid.md` section 11 not edited (Phase 5's
page).

## 2026-11-11 (agent, tree to format 3, group 7) — the THM projects in format 3; all 183 convertible files converted (`feature/upgrade-4-tree-format3-5`)

**Done** (`AK_tree_thm.md`; plan and state: `docs/upgrade/format3_tree.md`, group 7 row, §3a, §4, §5). The seven THM
files are format 3 (`18O_p_a_thm`, `17O`, `7Li_p_a`, `6Li_d`, `thm_coherent_background/toy`, examples
`c12c12_tumino2018`,
`o17_guardo2017_fit`): the whole test tree is converted. The checks edit THM records through `tests/lib/thm_records.sh`
(a data set joins an experiment by `thm <experiment>`, weights are `thm_weight`, folds are beam records) and
`tests/pyazr/thm_records.py`; the old member and weight lists go through the converter on the originals. Engine fixes
(`src/ThmBlockEdits.cpp`, `src/ConfigThm.cpp`): an experiment of inactive data sets only is ignored with the old warning
(was refused with "'segments=' is not key=value"); a project without a `<thm>` section no longer drops its data sets'
`thm_weight` tables silently, and an active data set naming an experiment no section defines is refused (it ran as
global).
The converter refuses five more THM member/weight forms the old engine refused. pyazr tests of AzrModel's THM API read
the originals (`old_tree.py`; list in the page's §5).

**Verified how**: page §3a (gate with the base counts; `tab-tree3e2`, the new reference, against `tab-trial19`:
only names, one listed class and run.log in the converted projects; `ctest -LE slow` 186/186, `-L slow` 10/10; 1 vs 4
threads identical; the 18O THM fit equal to the base's fit of the original; every rewritten check number for number
the base's; `tests/run_tests.sh` 145/1, the failure, `model_data`, fixed in the last commit and rerun).

**Open**: the refusal switch (what still reads old files: page §5, now without unconverted groups); the THM pyazr tests
of the old file class to move to `Project`; the engine accepts an old `segments=` list in a format-3 `<thm>` section,
and an experiment no data set names is refused with the old key's message "segments= is required" (both listed, not
fixed: `src/` beyond the brief).

## 2026-11-11 (agent, Phase 5 branch 3) — lab data compared in the lab frame (`feature/upgrade-5-lab-frame`)

**Done** (decision M3 (b), `energy_grid.md` §4 and its new §4.2a with the measured sizes).
- Stage 1 (1ba9772d): `FrameRow` keeps kinematics only (E, θ_c.m., J); `yCM`, `dyCM` and the value arguments of
  `ConvertRow`/`UnconvertedRow` are gone; `Segment::CMValue`/`CMError` are members computing J·y, J·δy for display
  (writers, API, diagnostics).  Byte-identical: all projects × 3 bases against `tab-trial18` (tag `tab-labframe1b`,
  728 / 701 / 728 identical).  One copy after Initialize: T2 11,872,460 -> 11,866,380 bytes, natural_neon
  142,705,164 -> 142,604,828, angle_mapping 501,376 -> 487,712 (probe `work/labframe/cb`).
- Stage 2 (ec234914): one residual, `DataResidual(f, J, y, δy, n)` = (f/J − n y)/(n δy), in the objective, the
  profiled norm, the least-squares rows and norm column, both analytic gradients (model cotangent × 1/J), the
  covariance pulls, the shift columns, the THM profiles (model and columns divided by J; `THMRows::m, Jm` stay c.m.
  for the band).  J = 1 gives the old bits.  `tools/classify.py --class=lab-frame-rounding` (cls.sh forwards it) with
  the probed list of projects with lab-differential rows; `tests/consistency/classify_test.py` and
  `tests/upgrade/lab_frame_residual_test.cpp` (T2; shown failing with f·J).
- Measured (scratch `work/labframe/`): outputs identical as printed except the χ² of `sum_lab_frame/models/sum_data`,
  `sum_free` (3.99e-17, data equal to the model; 1–2e-6 of itself, 7e-23 absolute); model values bit-identical; pulls
  within 1.3 ε of their natural scale (relative to the pull up to 2.7e-12, 1.2e-2 for the 1e-9 pulls of
  sum_lab_frame); χ² ≤ 2e-15 at the start and ≤ 2.5e-15 at the reference fits' end points (sum_free 2e-9);
  fit end points up to 9e-6 in a parameter (angle_mapping), χ² up to 4.3e-9 (spin_orbit_potential).

**Verified on ec234914** (binary `work/labframe/AZURE2-final`, scripts `work/labframe/final.sh`, `bench.sh`): all
projects × 3 bases (tag `tab-labframe2`) identical to `tab-labframe2a`; against `tab-trial18` with the class: 726 / 697
/ 724 identical, 0 / 2 / 2 lab-frame rounding, the 2 `run.log` of `sum_lab_frame` per basis differing; `ctest -LE
slow` 184/184 (check.sh --quick); `tests/run_tests.sh` 146 passed, 0 failed; calculations at 1 and 4 threads
identical (11 projects × 3 bases, 33/33); fits at 1 thread of stage 1 byte-identical to `b-trial18` (T2 lm, grad,
migrad, angle_mapping); converted-tree gate (`glab`) 483 (456) identical, 38 classified, 0 differ, 0 missing.
Cycles (perf stat, 1 thread, min of 6, alternated with `b-trial18`, load average about 18): T2 calculation 649.0 M
against 646.0 M (+0.5 %), 200-call T2 fit 50.32 G against 49.82 G (+1.0 %); a first series gave +1.0 % and +5.3 %
(spread of single runs 50–80 G); stage 2 against stage 1 on the fit, min of 4: −0.7 %.

**Open.**  The class does not accept the `run.log` χ² of `sum_lab_frame` (12 printed digits of a rounding-noise χ²): a
decision for the director (an absolute floor for a χ² that is zero to rounding, or the change as it is).  Reference
tag from this tree: `tab-labframe2` (= `tab-labframe2a`).  pyazr's docstrings (`azure2.py` `residuals`,
`residual_jacobian`) still write the residual with the c.m. error (pyazr is another branch's).

## 2026-11-11 (agent, SRIM stopping) — `stopping/` computes SRIM-2013's stopping powers, ported from LunaSRIM (`feature/upgrade-5-stopping-2`)

**Done.** The first branch's transcribed coefficients (Geant4's ICRU 49 rows, `elements.tsv`, Barkas table, joins,
Pierce–Blann scaling, ZBL universal fit, Chu/Yang forms, the ERYA comparison tool) are gone. `stopping/` now ports
LunaSRIM 1.0.0's `stopping.py` and `straggling.py` to C++ with SRIM's `SCOEF03.DAT` and `SNUC03.dat`
(`stopping/data/srim/`, unchanged, README with SHA-256, terms, citation; compiled in, or read from
`AZURE2_SRIM_DATA_DIR`): electronic stopping of any ion in any element, solid or gas with SRIM's phase factor;
pair-specific nuclear stopping; Bragg's rule with SRIM's compound correction (applied once, as SRIM's tables, TRIM and
LunaSRIM's energy loss do; LunaSRIM's `compound_stopping` applies it twice, following the SR module's list mode);
break points of the integrals at the branch boundaries and at the nuclear stopping's jump (reduced energy 30);
straggling TRIM (default) or Bohr; SRIM's published accuracy as the uncertainty (H 4.0 %, He 3.9 %, Li 4.8 %,
heavier 5.8 %; nuclear 10 % provisional); SRIM cited in the provenance text and in every `azure2-stopping` output.
Documents: `stopping.md` rewritten, `ARCHITECTURE.md` §7, `energy_grid.md`, `ROADMAP.md`, `FOR_THE_OWNER.md`,
`CLAUDE.md`, the summary page.

**Verified how.** Build 0 warnings; `ctest -R "stopping|consistency"`; tests against LunaSRIM's own values
(`stopping/tools/make_srim_reference.py`, 8889 rows): electronic (8 projectiles, 25 elements, both phases for gases,
every branch boundary) and compounds to 1e-10 (largest 1.6e-11), nuclear to 1e-10 (largest 1.1e-12, `log1p`),
eight layers' energy loss and straggling to 1e-5 (largest 4.6e-6, LunaSRIM's integrator); a shifted coefficient
column fails 1001 rows. `ctest -LE slow` (see the branch report). One ε 0.8–1.1 µs (before 3.1–3.5), one E_out
0.6–0.9 ms (before 1.6–2.3).

**Open.** Whether the public repository ships SRIM's two files (owner; LunaSRIM's README reads SRIM's terms as
non-commercial use only); the compound correction by name from `Compound.dat` (not done, the user gives a number);
a source for the nuclear uncertainty; the Python binding and the wiring (next branch).

## 2026-11-11 (agent, Phase 5 branch 1) — one level-matrix solve per compound energy (`feature/upgrade-5-energy-tasks`)

**Done** (commits 3b757db4, a7eef25c, 43c1c5bc and the documents).
- The point pass (`WorkingCopy::CalculatePointsInParallel`, now in `src/EnergyPass.cpp`) runs the plan's energy tasks
  (`CalculationPlan::EnergyTasks`): every unit -- a point without a fold or a sub-point -- with the same nucleus, route
  (A-matrix, THM, R-matrix) and bits of E_x (`AMatrixFunc::InEnergy`, now also a public static) is one task: one fill
  and factorisation, one T-matrix per group (entrance, exit, `ECAmplitudeSource`), one observable per unit
  (`azure2::SolveAtOneEnergy`); then the folds in a second parallel loop (`Unit::FoldSubPoints`, `FoldThinTarget`),
  which read the values already there. Tasks are kept between passes and rebuilt when a unit, its E_x or an
  amplitude source changes; largest first (1 + T-matrices + units), dynamic schedule.
- `Unit::Calculate` (the serial paths: model pass with `LeafValues::Calculate`, gradient, MCMC, reaction rates) has
  no sub-point recursion: a host and its mapped points once per sub-point index through the same `SolveAtOneEnergy`,
  then the fold. The external-capture re-energising of sub-points moved to `Unit::PrepareSubPoints` (before the pass).
- `AZURE_PLAN_CHECK=1`: the write check over energy tasks and folds (batched above 2,000 steps), then every unit solved
  alone and compared bit for bit with its shared solve (`CheckEnergyTasks`: routes and g^T A g of every J-group,
  T-matrix from the unit's provider, record). `AZURE_PROFILE` counts energy tasks and their units.
- The T-matrix key: the pair is not enough. A capture T-matrix carries the unit's external-capture amplitudes, and the
  set-up's amplitudes of two units at one energy can differ in the last bits (`tests/electron_screening/grad.azr`:
  1,602 of 10,125 units; outputs equal at their printed digits, the gradient fit diverged from the base at call 140).
  Groups are now (entrance, exit, source): computed amplitudes together, held ones by a hash of their bits (a shared
  hash with other bits is refused), by segment without the EC table.

**Verified** (on 43c1c5bc unless said; base `work/b-trial18`, reference `tab-trial18`; load 16–26 on 24 cores).
- All projects x 3 bases against `tab-trial18`: brune 728, no-brune 701, use-park 728 files, all identical (tag
  `tab-etasks2`; the same on a7eef25c, `tab-etasks`).
- `AZURE_PLAN_CHECK=1` on the whole tree, three bases: 262 checked evaluations, 3,727,625 units, 0 differ; the screening
  gradient fit under the check: 544 evaluations, 0 differ.
- `ctest -LE slow` 183/183; `upgrade_plan` (151 projects, 116,454 energy tasks for 1,014,989 units),
`upgrade_plan_runtime`
  (now also thick target, `target_effect_ranges`, the screening model), `upgrade_history` pass.
- Calculations at 1 and 5 threads identical: T1, T2, T3, 13N, natural_neon, thick target, energy_shift_sqrt, shift_fn,
  18O_p_a_thm, target_effect_ranges.
- Fits at 1 thread byte-identical to the base: T2 `--use-lm` and `--use-gradient` (200 calls), 13N (200 calls),
  energy_shift_sqrt, the screening fit (`electron_screening/grad.azr` on 13N's data, `--use-gradient`), 18O_p_a_thm.
- `tests/run_tests.sh` (OMP 2): 145 passed, 1 failed -- `component_mapping`, its check (b), see Open.
- Converted-tree gate (`REF=tab-trial8`, tag `getasks`): 483 identical (456 formal), 38 classified, 0 differ, 0 missing,
  exit 0.
- Speed: `energy_grid.md` §8 (per evaluation T2 106 → 45 ms, 13N 317 → 132, thick target 1,102 → 264, natural_neon
  3,207 → 1,718; T1, T3 not slower; THM unchanged). A 200-call T2 fit at 1 thread is not faster as a whole (cycles
  0.96x, within the ±12 % scatter): its 206 evaluations take 24.0 s instead of 28.9, its 19 gradients 29.0 s
  instead of 33.7 (their forward solves are still per host), and the hand-over and Legendre tables are unchanged.
  At 5 threads (2 runs each, load 16–18) the same fit took 17.0 s instead of 14.8 (cycles 0.89x, scatter ±10 %):
  not explained; candidates are the serial rebuild of the energy tasks (each pooled copy and each segment-local
  shift pass rebuilds them: 1.7 % of the 1-thread run) and the units of a task scattered in memory.
  `target_effect_ranges`: solves 5,862 → 2,208 per calculation (2,164 tasks and 44 probes), the point pass not slower
  although every sub-point is now calculated.

**Kept, on purpose.** The mapping (`MapData`, `MapLeaves`, `SharesCalculation`): it still groups the folds and their
automatic-tolerance decision, the energy-dependent values and external-capture amplitudes of mapped points (and the
layout of `intEC.dat`), the serial paths and the analytic gradient (`BuildGradItems`), which must keep its grouping
until the per-energy gradient. The probe copies: the bare value of blended and automatic-tolerance points.

**Open.** `tests/component_mapping/check.sh` (b) fails (`run_tests.sh`): it asserts that the component mapping saves
level-matrix solutions, and the energy tasks now share them with or without it (2,769 = 2,769; extrapolation 271 =
271); its outputs (a) and gradients (c) pass. The check belongs to `tests/` (not edited here): it needs the director's
decision (compare against the base count, or drop (b)). The mapping goes with the per-energy gradient branch; THM
gains nothing (cache misses on the units' rows);
in a fit with THM-only finite differences the energy tasks are rebuilt each time the work list alternates (2 % of a
THM fit); ThreadSanitizer not run.

## 2026-11-11 (agent, brief Y2) — external capture reads the pathway table; the trees are deleted (`feature/upgrade-4-pathway-table-2`)

**Done** (`docs/upgrade/cnuc_split.md` §11.12 has the detail).
- Defect 1 of §11.9 fixed first: external capture found the exit to the photon pair at the photon pair's *number* in
  the list of exits; with a beta pair or a pair of another nucleus numbered below it, it read another exit or past the
  end.  Now `ExternalCaptureExit` (`include/ExternalCaptureExit.h`) looks the exit up by pair and refuses an entrance
  without a pathway, naming the data segment.  Test project `tests/ec_exit_position` (the base: empty `intEC.dat`,
  SIGSEGV) and ctest `upgrade_external_capture_exit`.
- The memo of the shifted capture amplitude cannot depend on the history of a copy (the final state's capture energy
  is written once, in `CNuc::Initialize`); it now keeps (energy, amplitude) only.  `history_test` gained a project with
  external capture and a fitted shift (`--external-capture`) and the sequence B -> C on every project.
- External capture reads the table (step 3): `src/ExternalCaptureSetup.cpp` (walk, signature, amplitude loops of the
  data segments and of the leaves, the latter moved from `WorkingCopy.cpp` as `LeafECAmplitudes`),
  `src/UnitExternalCapture.cpp`.  The intEC signature is unchanged; base and new builds read each other's files.
- The trees are deleted (step 8), after a commit that ran every project with them unbuilt and their accessors
  aborting: `SortPathways`, `CalcAngularDists`, `PPair::decays_`, the classes `Decay`, `KGroup`, `MGroup`, `ECMGroup`,
  `KLGroup`, `Interference` (12 files, 971 lines), `tests/upgrade/pathways_from_cnuc_test.cpp`.  `CNuc.cpp` 1,774 ->
  1,548 lines.
- The double final state is refused since the capture-A_y branch, and `tests/ec_final_state_ambiguous` (d), (e) test
  exactly that case: nothing added.

**Verified** (build `work/b-path2`, 0 warnings): all projects x 3 bases (tag `tab-path2-final`) against `tab-trial18`:
728 / 701 / 728
identical (every `intEC.dat` among them), 22 files per basis new (the two projects of `tests/ec_exit_position`); the
same on the commit with the trees unbuilt and their accessors aborting (`tab-path2-abort`, no abort); `ctest -LE
slow` 182/183 (the failure: `upgrade_predictions`, the two new projects missing from `predictions.txt`; regenerated,
12 lines added, passes); `ctest -L slow` 10/10; `tests/run_tests.sh` 146 passed, 1 failed (`cli_help`, under load;
passes alone, as with the base); the converted-tree gate 483 (456) identical, 38 classified, 0 differ, 0 missing;
calculations at 1 against 5 threads byte-identical (T1-T3, 13N, `beam_profile_kernel`, `energy_shift_sqrt`, the
shifted external-capture project of `history_test`); fits at 1 thread byte-identical to the base (T2 `--use-lm
--max-calls=200`, 13N and `energy_shift_sqrt` `--max-calls=200`); `--no-ec-table` on 13N and the shifted project
identical to the base; 13N's `intEC.dat` written by each build read by the other without recomputation.  Cycles
(perf stat, 1 thread, `--no-coulomb-table`, min of 6-8 alternated with `b-trial18`, a loaded machine): T2 calculation
3.294 G against 3.304 G, 13N 25.28 G against 25.39 G, set-up (13N with one point) 0.987 G against 1.026 G;
instructions -0.05 %, -0.08 %, -0.36 %; peak RSS 57.4 / 149.9 / 12.1 MB against 57.6 / 150.3 / 12.6 MB.

**Open.** The external-capture amplitudes of the set-up still live in the units of the set-up's working copy; a copy
drops them and recomputes each at its first evaluation.  A comment in `src/ECIntegral.cpp` (out of scope) still names
`CNuc::SortPathways`.  Steps left of `CNuc`: the `GetPair`/`GetJGroup` readers of channel data and step 7.

## 2026-11-11 (agent, Phase 5 branch 4) — the yield of a layer of finite thickness (`feature/upgrade-5-layer-yield`)

**Done.** `quantity yield` on a `target` of `kind layer`: Y(E₀) = 10⁻²⁴ ∫_{E_out}^{E₀} ⟨σ⟩/ε dE, ∫ dE/ε = N (lab ε per
active nucleus, the user's stopping formula), per particle or per µC. Today's target kernel (`src/TargetFold.cpp`)
without the division by N; E_out solved on the point's own sub-points with `GridQuadrature`'s rule
(`include/LayerExit.h`), so a constant σ gives 10⁻²⁴ σ N to 1e-15; the grid spans 2 % beyond the layer. A layer yield
is an effect of its own (`TargetEffectSpec::layerYield`, set in `experiment/src/BuildEffects.cpp`); fold windows on it
are refused; the set-up refuses a non-positive stopping and a layer the beam does not leave (`src/TargetLayer.cpp`,
codes `E_LAYER_STOPPING`, `E_LAYER_EXIT_BELOW_MODEL`). ∂Y/∂N stays central differences. Test project
`tests/layer_yield` (README with the derivations; check.py computes every closed form); `layer_kernel` and
`pyazr_layer_yield` in ctest. Docs: `data_formats.rst`, `format3.rst` (generated), `data_kinds.md` §4.6,
`project_format3.md` §7.2/§9.1, `energy_grid.md` §8, ROADMAP, FOR_THE_OWNER, the targets skill.

**Verified.** All projects × 3 bases against `tab-trial18`: 728 / 701 / 728 identical, the new project's files only-new;
`ctest -LE slow` 185 (after `predictions.txt` got the two new projects: 17 lines added); `tests/run_tests.sh` 146
passed,
1 failed (`cli_help` under load; passes alone on this build and on the base); spec equivalence 200 files, 185
converted, 183 equal, 2 recorded, 0 differ, 15 not converted; converted-tree gate (REF tab-trial8) 483 / 456 / 483
identical, 38 classified, 0 differ, 0 missing; menu 1 at 1 and 4 threads identical on the new project and T2 (4, not 6:
the thread budget); the fit of `tests/layer_yield` identical at 1 and 2 threads but for Minuit's call count.

**Open.** The stopping module (scale nuisance, stoichiometric active fraction) plugs in where the sub-points' stopping
is set; an analytic ∂Y/∂N = 10⁻²⁴ σ_lin(E_out) in `AZUREGrad.cpp` if wanted; a target with a free screening or
temperature used by a layer yield and a σ data set is refused (two effects, `E_FOLD_CONFLICT`); the comment in
`tests/units_scales/models/yield_story.azr` still calls the layer yield a Phase 5 kernel (an existing test's .azr, not
edited); `convert/tests/spec_equivalence_test.cpp` does not compare `layerYield` (old lines never set it; convert/ is
out of scope).

## 2026-11-11 (director, 14:00) — merged: the stopping library and the DWBA transfer library (leaves); stopping powers will be SRIM's

Merged on the data branch (`--no-ff`), not pushed: `feature/upgrade-5-stopping-1` (`stopping/`) and
`feature/upgrade-7-transfer-1` (`physics/transfer`). Both are leaf libraries: no engine file includes them, so no
project output can move; the check was the build (0 warnings) and `ctest -LE slow` on the two together (197/197).

What the director ran himself, beyond the agents' reports:
- **Stopping against SRIM-2013** (LunaSRIM 1.0.0, electronic part, p and α in H, He, C, N, O, Ne, Al, Ti, Ta, Au,
  0.02–10 MeV, gases as gases): most values within 3 %; p in Ta +6 to +7 % at 0.1–0.3 MeV; p in C −11 % and α in Au
  −17 % at 0.02 MeV; α in C +7 to +11 % below 0.5 MeV; α above 5 MeV −3 to −5 % in several elements.
- **FRESCO** on the committed (d,n) finite-range case: the cross sections equal `reference.txt` digit for digit.

**Decision of the PI (2026-11-11), which replaces the licence question of the previous entry:** "I would use the
SRIM stopping powers from my LunaSRIM repository. They are free to use as long as SRIM is cited." So:
- the electronic coefficient table of the first branch (transcribed from Geant4's `G4BraggModel.cc` /
  `G4BraggIonModel.cc`) and its Barkas table go; the stopping powers become SRIM-2013's (electronic for any ion,
  gas/solid phase factor, compound correction, pair-specific nuclear stopping), ported to C++ from LunaSRIM
  (github.com/skowrons94/LunaSRIM, MIT) with SRIM's data files `SCOEF03.DAT` and `SNUC03.dat`;
- SRIM is cited in the data directory, in `stopping.md` and in the provenance text the library prints;
- `ARCHITECTURE.md` §7 ("no SRIM, no ERYA, no LunaSRIM") and `stopping.md` §4, §5, §9, §12 are rewritten by that
  branch (`feature/upgrade-5-stopping-2`, launched from this head). Until it is merged, `stopping/` computes with the
  first branch's table: do not wire it into the engine before.
- **For the owner:** LunaSRIM's own LICENSE says the SRIM data files are not under its MIT licence and stay under the
  terms of SRIM's authors. They are copied into this repository on the PI's statement that they are free to use when
  SRIM is cited; whether AZURE2's public repository ships them is the owner's decision (the alternative is a download
  step, with the library reading the two files from a directory the user names).
- The three coefficient discrepancies the first branch reported are moot (the rows go). The Chu table and the
  Yang–O'Connor–Wang constants are not needed for now: the straggling becomes SRIM's (Bohr with TRIM's relativistic
  factor and low-energy limit), propagated as before.

Judgements on what the agents asked:
- DWBA, post interaction: the default stays the whole V of the transferred particle with the core, and the choice
  (whole / nuclear part) becomes an explicit setting when the library is wired to the THM observable.
- Assigned uncertainties of the stopping power: provisional priors until the SRIM branch replaces them by SRIM's
  published accuracy; the stopping scale per (material, projectile) (D6) takes its prior from there.

Defects recorded, not fixed here:
- `physics/coulomb` COULCC returns F and G with the wrong sign at isolated points beyond the turning point (l = 1,
  η = 0.493, ρ = 8.55). The engine does not call COULCC there (it uses COULSM); only `physics/src/transfer/
  RadialNumerics.cpp` and tests do. To be fixed or the routine dropped when the continuum part of the DWBA is written.

Running: pathway table 2; one solve per compound energy; finite-layer yield (resumed to write its report); lab-frame
comparison; THM projects to format 3; the gaps of format 3; SRIM stopping (new).
Still owed to the PI: the count of engine files each new quantity touches (the finite-layer yield is the first), and
the extension drill before the meeting with the owner.

## 2026-11-11 (agent, brief AG) — `physics/transfer`, first part: potentials, bound states, distorted waves, DWBA to a bound state against FRESCO (`feature/upgrade-7-transfer-1`)

**Done.** A leaf library in `physics/` (target `azure2_physics_transfer`, namespace `azure2::physics::transfer`;
nothing of the engine uses it): `OpticalPotential` (Woods-Saxon volume, derivative surface, Thomas spin-orbit in
FRESCO's 2 l.s convention, uniform-sphere Coulomb; three explicit radius conventions; `TermsAt` lists the terms as
FRESCO's fort.34), `BoundState` (depth search by Sturm count and Casoratian, ANC, rms radius, D0; Whittaker and Yukawa
tails), `DistortedWaves` (spin 0, 1/2, 1 with spin-orbit; S, sigma/sigma_R, reaction cross section),
`Transfer` (post-form zero range with recoil and finite range, no remnant, interaction the whole or the nuclear part of
V_sx; explicit spin sum), `ThmVertex` (M_l and the full form factor in the plane-wave limit), `RadialNumerics`.
Formulas, conventions, validation table, limits, the continuum plan and the wiring: `docs/upgrade/transfer.md`.
FRESCO benchmarks with inputs and outputs in `physics/tests/data/fresco/` (8 cases), regenerated by
`physics/tests/data/make_fresco_benchmarks.py --fresco <binary>`; FRESCO receives our constants through its
`unitmass` and `finec`.

**Verified.** Tests `physics_transfer_bound_state`, `_distorted_waves`, `_thm_vertex`, `_fresco` pass; each failed
once with an injected error (listed in the report of the branch). Against FRESCO: bound-state depth, rms, D0, ANC to
the printed digits; S to 8e-6, sigma/sigma_R to 5e-5, sigma_R to 2e-6; transfer dsigma/dOmega (0-90 deg, > 1 % of the
peak) 15N(d,n)16O finite range 3.3e-4, zero range 7e-6, 12C(6Li,d)16O finite range 1.8e-3 (FRESCO's own step
dependence ~1e-3), zero range 4e-5. The engine's THM weight M(E) of `thm_distortion_test` case `optical` (E = 1 MeV)
recomputed with the library: 2.7054435e-4 against the engine's 2.7054438e-4 (-9e-8 at h = 0.01 fm).
`git diff --stat 2363abac -- src include api pyazr experiment convert` is empty.

**Open.** Spin-orbit in the distorted waves of a transfer; remnant and prior form; level functions; the continuum
(FRESCO bins); global optical sets; 18O(3He,d)19F and 12C(14N,d)24Mg benchmarks (not done); l > 40 (AngularCoupling's
limit); the binding of the library in pyazr. Found outside the branch: COULCC (`physics/coulomb`) returns F and G
with the wrong sign at isolated points beyond the turning point (transfer.md section 10); the library uses COULFG there.
With the whole V_alpha-d (Coulomb tail) as post interaction FRESCO does not converge in its non-local range
(transfer.md section 9).

## 2026-11-11 (agent, stopping module, first branch) — `stopping/`: AZURE2's own stopping and straggling as a leaf library (`feature/upgrade-5-stopping-1`)

**Done.** `stopping/` (design and numbers: `docs/upgrade/stopping.md`): materials from formulas or fractions with
phase and density; electronic stopping of p and α from the ICRU 49 coefficient sets (Andersen–Ziegler 1977 and
Ziegler 1977 forms, transcribed from Geant4 and checked against PSTAR/ASTAR) joined by log-log interpolation
between 1 and 3 MeV/u to Bethe with Barkas–Berger shell, Ashley–Ritchie–Brandt Barkas and Bloch corrections; d, t,
³He by velocity, heavier ions by Pierce–Blann effective charge; ZBL nuclear stopping; Bragg additivity with a named
place for the compound correction; user tables and formulas with provenance; an assigned uncertainty on every value;
N(E_in, E_out), E_out, range and their stopping-scale derivatives; Bohr straggling, Chu's factor and Yang's
charge-state term as forms taking data (no data shipped), Tschalär propagation with the beam spread. Two data tables
in `stopping/data/` with a source per row. The CLI `azure2-stopping`; the ERYA comparison
`stopping/tools/compare_erya.sh`. Consistency check 32 (stopping/ includes only stopping/, physics/ and the standard
library); `stopping/tests` registered; `stopping` in the size ratchet.

**Verified how.** Build 0 warnings; the six unit tests and the CLI test pass, each unit test shown failing on
injected errors (list in the report of this branch); consistency test, also shown failing on an injected include and
an unregistered test; `scripts/check.sh --quick` (consistency, build, `ctest -LE slow` 190/190); T1–T3 and
`tests/target_density` x 3 bases against `b-trial18`: 71/71/71 files identical.

**Open.** Coefficients that do not reproduce ICRU 49 (p in C below 40 keV, α in Au below 0.5 MeV; stopping.md §9.1);
no Chu table and no Yang constants; assigned (not fitted) uncertainties; no Python binding (two lines in
`api/src/pybind_api.cpp`, stopping.md §11); the Geant4 notice question for the transcribed numbers (§12). Next:
stopping.md §13.

## 2026-11-11 (director, afternoon) — merged: 176 of 183 project files in format 3; the engine edits the nuclear part of a project; Phase 5 decided and started

**Merged on `feature/upgrade-4-data`, then on `feature/upgrade`** (all `--no-ff`, not pushed):
- `feature/upgrade-4-tree-format3-4` (a573364b): group 8 of `format3_tree.md` (T1–T3 and the files left old in
  converted directories).
- `feature/upgrade-4-pyazr-records` (4d6e12de): the engine reads and edits `<levels>` rows, `<config>` entries and
  `<thm>` lines of a format-3 project; `pyazr.Project` offers the edits; `save_fit`, bands, angular distributions,
  tabulation, model averaging and the level transform work on a format-3 project without `AzrModel`.
- `feature/upgrade-4-tree-format3-3` (05a6898c): group 6 (natural targets). 176 of the 183 convertible files are
  format 3; left: group 7 (THM).
- The PI's decisions on Phase 5 (D2, D3, D5, D6, D7: "yes to all", `energy_grid.md` §12).

**Verified by the director** on the tree of the three merges together (rehearsal `trial19`, clean build, 0 warnings;
`work/director/verify_trial19.out`):
- All projects × 3 bases against the tree with group 8 and pyazr only (`tab-trial19a`, itself identical to the group-8
  agent's tag, 732 / 705 / 732 files): 713 / 686 / 713 identical; every number of every output equal (largest relative
  difference 0). What differs is text: `param.par` of `natural_neon` (parameter names), three sample headers of
  `natural_neon` ("sample natural" for "sample 1"), 15 `run.log` of the converted natural-target and reaction-rate
  models (the old-format note gone; a missing data file is an error in format 3 where the old path warned).
- `ctest -LE slow` 186/186; `ctest -L slow` (2 threads) 10/10; `tests/run_tests.sh` 146 passed, 0 failed; 12
  calculations identical at 1 and 6 threads; model-data χ² identical in three bases; spec equivalence 198 files, 183
  converted, 183 equal, 0 recorded, 0 differ, 15 not converted; gate 483 (456 formal) identical, 38 classified, 0
  differ, 0 missing. Reference tag from here: `tab-trial19` (base build `work/b-trial19`).

**Director's judgements at the merge.**
- Conflicts were resolved by hand in `src/CMakeLists.txt`, `tests/consistency/repo_consistency_test.py`,
  `docs/upgrade/format3_tree.md`, `experiment/src/BuildEffects.cpp`, `tests/pyazr/predictions_test.py`; comments in
  `convert/` that still named `EData` were renamed (check 31).
- Accepted for now, to be tightened by the gaps branch: the total tolerance of `scattering_length` 1e-6 → 1e-4
  (`modelData.out` prints five digits); `cli_route_exit` accepts any non-zero exit.
- Accepted as right: format 3 refuses a Doppler temperature ≤ 0 and an areal density of 0; the converter refuses what
  the old engine refused; `add_level` of pyazr does not renumber the other levels.
- To fix (in the gaps branch): a refused format-3 project exits 255 with "Could not open"; a parameter file whose
  names do not match is ignored without a word (`natural_neon`: χ² 36346 instead of 25647 before the `.sav` was
  renamed); `sigma_total` on a natural target; a detector resolution on beta-delayed data; a data set's `thm <name>`
  that names no experiment is not refused by validation.
- Found by the pyazr agent on the *old* path and not fixed, because that path is to be deleted: `tabulate` on an old
  project with `<targetInt>` applies the target effect to the tabulation (13N: 2.5e-12 for 6.5e-7 b);
  `angular_distribution` on the old 13N returns NaN.

**Running** (2026-11-11 afternoon):
- `feature/upgrade-4-tree-format3-5`: the THM projects to format 3 (group 7).
- `feature/upgrade-4-format3-gaps`: the defects and gaps listed above, the converter's open questions, the model
  generators that still write old files, the set-up profile.
- `feature/upgrade-4-pathway-table-2`: external capture on the pathway table, the old trees deleted.
- `feature/upgrade-5-energy-tasks`: one level-matrix solve per distinct compound energy (bit-identical).
- `feature/upgrade-5-layer-yield`: the yield of a finite layer as a quantity.
- `feature/upgrade-5-lab-frame`: data compared in the laboratory frame; the c.m. copies of the data deleted.
- Finished and being merged as a second batch: `feature/upgrade-5-stopping-1` (`stopping/`) and
  `feature/upgrade-7-transfer-1` (`physics/transfer`), both leaf libraries nothing links.

**Open for the PI and the owner.**
- `stopping/data/electronic_coefficients.tsv` and the Barkas table were transcribed from Geant4's source (ICRU 49's
  numbers; Geant4 Software License). Whether they may be shipped under AZURE2's GPL v3, and with which notice, is the
  maintainers' decision; taking the numbers from the reports themselves or from a fit to the IAEA data base removes
  the question.
- Extensibility is counted from now on: every branch that adds a quantity or an apparatus effect reports the engine
  files it had to touch; a drill with a fresh agent and only the recipes is run before the meeting with the owner.

## 2026-11-11 (agent, pyazr records) — the engine edits the nuclear part; pyazr's modules off `AzrModel` (`feature/upgrade-4-pyazr-records`)

**Done** (`AC_pyazr_off_azrfile.md`; design and state: `docs/upgrade/pyazr_on_records.md`, sections "Edits of the
nuclear part", "Remaining uses", "Tests"). The engine reads and writes `<levels>` rows (`io/LevelRowText.h`,
`experiment/LevelEdits.h`), `<config>` entries (`io/ConfigSection.h`, now the reader of `Config::ReadConfigFile`) and
`<thm>` lines (`include/ThmBlockEdits.h`); bindings in `api/src/pybind_nuclear_part.cpp`. pyazr 2.32.0: `Project`
levels / pairs / config / thm reads and edits (`pyazr/nuclear_part.py`); `save_fit` (now `pyazr/snapshot.py`),
`scattering_length_chi2`, the potential records, `bands`, `angular`, `tabulate`, `modelavg`, `transform` work on a
format-3 project without `AzrModel` (old files: one format test, the old branch). Nine pyazr tests moved to the
format-3 files; three new tests (`io_level_row_text`, `io_config_section`, `pyazr_project_nuclear_part`).

**Verified how**: see the branch report (all projects x 3 bases against `tab-tree3b`, the converted-tree gate,
`ctest -LE slow`, `ctest -L slow`, `tests/run_tests.sh`, 1 vs 4 threads; `save_fit` on 13N and on the THM toy with
a coherent background: the snapshot chi2 equals the old path's to every digit).

**Open**: `nuclei.merge_projects` (waits for the natural-target records of group 6); 13 tests still read originals
through `old_tree.py` and 12 more use `AzrModel` on format-3 files (the page's table says which and why); the skills
and user pages that teach `AzrModel` (`.claude/skills/azure2-*`, `particle_pairs.rst`, ...).

## 2026-11-11 (agent, tree to format 3, third part) — the natural-target projects (group 6); THM not started (`feature/upgrade-4-tree-format3-3`)

**Done** (`AA3_tree_natural_thm.md`; plan and state: `docs/upgrade/format3_tree.md`). Group 6 is format 3: natural_neon
(tests and example, the example's `.sav` renamed to the format-3 parameter names), the natural_target models but
`M_thick`, `reaction_rate/models/n12C.azr`, `p12C.azr` (140 converted files in all). Rewritten on records:
natural_neon, natural_target, shared_effect_pairs, component_subpoint_angle, doppler_fgm and transmission_yield with
their generators (format 3; the old generators kept in the fixtures for the pyazr tests of the old reader);
natural_target/make_models.py deleted (README). The converter refuses six more kinds of old line the old engine
refused and the conversion used to run (component naming a missing pair key, a sample `<abundances>` lacks, doppler
fit forms with vary 2 / sigma < 0 / Tmax <= T, transmission and yield on one line, a thicknessgroup member without its
own block); format 3 refuses a Doppler temperature <= 0 and a sample areal density of 0 (`E_OUT_OF_RANGE`, both ran).

**Verified how**: all projects x 3 bases against `tab-tree3b` (tag `tab-tree3c-1`, the new reference): identical but
natural_neon's param.par names and sample-file headers; gate (`gtree3c`) 483 (456) identical, 38 classified, 0 differ,
0 missing; spec equivalence 198 files, 183 converted, 183 equal, 0 recorded, 0 differ, 15 not converted (bad_key,
bad_sample now refused); ctest, run_tests.sh and the thread comparison: see the branch report.

**Open**: group 7 (THM, not started; the page's group-7 row lists the readers of 17O measured on a trial conversion),
group 8; `total_optical` (a natural sum of totals has no format-3 form: three options in the page's §4, for the
director); `M_thick.azr`; the converter could answer its "two areal densities" question itself for a doppler-only
`<targetInt>` line on a transmission line (one sample target with the doppler field gives byte-identical outputs).

## 2026-11-11 (agent, tree to format 3, group 8) — T1–T3 and the rest of the tree (`feature/upgrade-4-tree-format3-4`)

**Done** (`AA4_tree_rest_group8.md`; plan and state: `docs/upgrade/format3_tree.md`). The 36 old files of group 8 are
format 3 (`upgrade/T1–T3`, `relativistic`, `scattering_length`, `channel_space`, `beta_delayed_alpha`,
`identical_entrance_reaction`, `identical_pp_res`, `hybrid_potential`, `sum_lab_frame`, `ratio_frame` models): 162
of the 185 convertible files. Every check, unit test and tool that reads them edits records (list in the page, §3
row 8); the golden comparisons of T1–T3 compare the old golden files as the format-3 projects write them (the
converter's parameter names, the classes of `tools/format3_gate.py`); the generators of T1–T3 write format 3 through
the converter and give the committed files byte for byte; `light_excitation/make_models.py` (run by its check)
writes records. `tools/convert_tree.sh` drops a committed side file (`prior_cov.f3.txt`) before converting an
original, so the gate converts T2 again.

**Verified how** (page §3a): converted-tree gate `gtree3d` 483 (456) identical, 38 classified, 0 differ, 0 missing;
all projects x 3 bases `tab-tree3d` (the new reference) against `tab-tree3b`: only the converted projects differ, by
listed classes; spec equivalence 198/185/183/2/0/13; `ctest -LE slow` 183/183; `ctest -L slow` (OMP 2) 10/10;
`tests/run_tests.sh` 145 passed, 1 failed (spin_orbit_potential, fixed, rerun passed); 1 vs 4 threads identical on
T1-T3, 13N, identical_pp_res, two relativistic models; T2 benchmark at 1 thread 10899.8 / 691.608.

**Open**: groups 6 and 7 (another branch); 23 developer generators that still write old files (page §5);
`beta_delayed_alpha` (e), a Gaussian resolution on a beta-delayed spectrum, which format 3 cannot state (director);
defects found (page §4): a refused format-3 project exits 255 with "Could not open ... Check that file exists."
(`src/AZURE2.cpp`), a THM experiment whose data sets are all inactive is refused (`ConfigThm.cpp` writes
`segments=`), `sigma` of a label that reads elastic is refused even with `exit_pair` naming another pair.

## 2026-11-11 (agent, Phase 5 design) — the energy grid, kernels, lab-frame comparison, finite layers (`feature/upgrade-5-design`)

**Done.** The design note `docs/upgrade/energy_grid.md`, measured on 2363abac (no code changed). Measured with an
instrumented scratch build (`work/b-p5design`, patch `work/p5design/p5_dump.patch`, never committed) and the Release
build `work/b-trial18`; scripts and runs in `work/p5design/`.
- An evaluation solves the level matrix once per unit; the units lie on far fewer distinct compound energies, bit for
  bit: 13N 45,872 solves on 2,728 energies, `thick_target_resonance` 48,239 on 2,149, `18O_p_a_thm` 14,689 on 704,
  `natural_neon` 291,847 on 56,703, T2 24,300 on 5,442; T1, T3 and the shift projects 1.05-1.5x. A shared mesh at
  today's node density adds at most 2.2x (T2, natural_neon). Level matrix + T-matrix are 40-80 % of a fit's point
  pass (13N 166 of 258 ms, thick target 643 of 795 ms).
- M3 (b) (lab-frame comparison of shifted data) is an identity for the chi^2: today value and error are both
  multiplied by the same J(nu), so the pull is frame-invariant; the decision moves no fit beyond rounding.
- ERYA is not called by the engine (only by `tests/reference/compound_stopping_reference_test.cpp`).
- Ten branches proposed (section 8): bit-identical energy tasks first (one solve per distinct energy), lab frame and
  layer yield at once, stopping wiring after `stopping-1`, kernels out of the unit, then the shared mesh (gate against
  today's own twofold-refined grids), shift kernel, observables, THM. Eight decisions for the PI (section 9).

**Open.** Everything in the note is a proposal for the director and the PI; the speed of the branches is inferred.

## 2026-11-11 (director) — merged: the data classes replaced; two thirds of the tree in format 3; what runs now

**Merged on `feature/upgrade-4-data`, then on `feature/upgrade`** (both `--no-ff`, not pushed):
- `feature/upgrade-4-tree-format3-2` (92b046a8): groups 4b and 5 of `format3_tree.md` (shift, target, screening,
  covariance and apparatus projects), 126 of 185 project files in format 3; the converter refuses what the old engine
  refused and keeps an automatic tolerance without windows; a fitted normalization with a prior centred at 0 is
  refused with a diagnostic code.
- `feature/upgrade-4-data-structure-3` (2363abac): `DataSet` (immutable, `shared_ptr<const DataSet>`), `WorkingCopy`,
  `Segment`, `Unit`; `EData`, `ESegment`, `EPoint` and their files deleted (consistency check 31).

**Verified by the director** (clean builds of the rehearsal trees, 0 warnings).
- The tree branch alone (rehearsal `trial17`, against `tab-trial16`): 685 / 658 / 685 files identical; the others are
  `param.par` of 9 converted projects (parameter names only, values equal) and 36 `run.log`; `ctest -LE slow`
  183/183; `tests/run_tests.sh` 146 passed, 0 failed; 12 calculations identical at 1 and 6 threads; model-data χ²
  identical in three bases; spec equivalence 198 files, 185 converted, 183 equal, 2 recorded, 0 differ; gate 483
  (456 formal) identical, 38 classified, 0 differ, 0 missing; `ctest -L slow` 10/10.
- Both together (rehearsal `trial18` = 2363abac's tree, against `tab-trial17`): 728 / 701 / 728 files, all
  identical; `ctest -LE slow` 183/183; `run_tests.sh` 146 / 0; 12 identical at 1 and 6 threads; χ² identical in three
  bases; spec equivalence and gate with the same counts; `ctest -L slow` 10/10. Reference tag from here: `tab-trial18`.

**Director's judgements at the merge.**
- Two output texts still name the old classes on purpose (`src/WorkingCopy.cpp` "see EData::ECSignature)",
  `convert/src/BuildData.cpp` "(the same index, EData::MakePoints)"): they are bytes of committed outputs; check 31
  lists them as the only exceptions. They go with the outputs' next regeneration.
- The agent's "pyazr behaviour change" (a target's `density` after an evaluation) is not visible to a user: the only
  binding that returns a density builds a fresh `TargetEffect` from the line (`api/src/pybind_api.cpp`,
  `parse_target_effect_line`), and the hand-over reads `EffectValues(n).density`. Nothing to fix.
- `Unit` is 480 bytes but not yet a plain record: about 130 member functions, the apparatus fold among them
  (`src/TargetFold.cpp`). That is Phase 5's work (kernels leave the unit), not a debt of this branch.

**Procedure changed (the PI, 2026-11-11: "Can't we go faster?").** A merge is committed on the data branch as soon as
all outputs are identical, and agents start from it; `feature/upgrade` takes it after the whole verification. One
verification per batch of merges, its parts side by side: `work/director/full_verify2.sh WT TAG NT [REFTAG]`. Work that
depends on nothing in progress starts at once (leaf libraries, the next phase's design).

**Running** (branches from ec04bd55 or 2363abac):
- `feature/upgrade-4-tree-format3-3`: natural targets and THM projects to format 3.
- `feature/upgrade-4-tree-format3-4`: the last group (files left old in converted directories, T1–T3) and the model
  generators.
- `feature/upgrade-4-pyazr-records`: pyazr's own modules and tests off `azrfile.py`.
- `feature/upgrade-4-pathway-table-2`: the decay-position defect of external capture, external capture on the
  pathway table, the trees deleted.
- `feature/upgrade-5-stopping-1`: the stopping-power library `stopping/` (leaf; links `physics/` only).
- `feature/upgrade-7-transfer-1`: potentials, bound states, distorted waves and DWBA transfer in `physics/`, checked
  against FRESCO 3.5 (installed outside the repository: `/data0/skowrons/AZURE2-upgrade/opt/fresco/README.md`).
- `feature/upgrade-5-design` (finished, under review): `docs/upgrade/energy_grid.md`, the design of Phase 5.

**Open** (unchanged unless said): the second validation round waits for the deletion of the old reading path; the
UDR fold through a point has no end-to-end test; comments in `pyazr/azrfile.py` and `pyazr/datasets.py` name the old
classes (the files go with the old path); the 6-thread fit's peak memory (+13 % at milestone 2) is to be measured
again now that a copy is smaller; ThreadSanitizer and the Sphinx build have not been run on this state.

## 2026-11-11 (agent, brief X3, milestone 3) — `DataSet` and `WorkingCopy` apart; the old data classes gone (`feature/upgrade-4-data-structure-3`)

**Done** (commits 075af6db, b1d3daaf, the file moves, the class renames, these docs).
- `DataSet` (`include/DataSet.h`, `src/DataSet.cpp`): what the set-up leaves fixed -- the segments (`Segment`, settings
  and start values, priors, measurement, terms with their settings and frame at load), the target effects, the leaves,
  the `PointSetUp` records of the points (adaptive-grid layout, UDR kernel), the covariance groups, the THM
  experiments, the parameter layout, the prior covariance, the nuclei.  A working copy holds it as
  `shared_ptr<const DataSet>`: the set-up writes it through `WorkingCopy::Draft()`, `Initialize` freezes it, a copy
  shares it; a setting changed later (pyazr's prior covariance, `ReplaceDiagonalPriors`, a parameter layout laid out
  again by `FillMnParams`) makes a new one (`ReviseDataSet`).  Shared and written after the set-up, named in the
  headers: a covariance group's frozen relvector model (its mutex; refrozen only between minimizations) and the
  node-table cache (its shards' mutexes).
- `WorkingCopy` (was `EData`): per segment a `SegmentState` (`include/SegmentState.h`: norm, shifts with their last
  applied values, the shift function `Equation` with its parameters, sample thickness; per term the unit range, the
  present frame, the shift tables), `TargetEffectValues` (U_e, density, Doppler T: the fold and the penalties read
  them there), the THM coherent values, the calculation plan (rebuilt when a shift moves), the unit vector.
  `Segment`'s functions that touch state take the working copy or the state and are const.
- The unit (`Unit`, was `EPoint`) 608 -> 480 bytes: its grid layout and UDR kernel are data-set records named by an
  index, and the angle-limit message of a failed lab -> c.m. map is returned by `Segment::ConvertRowToCM`.  The data
  iterator is gone.  The external-capture amplitudes stay in the unit (`Clone` drops the set-up's, as before).
- Files and classes renamed (`data_structure.md` §3 lists them); consistency check 31 keeps the old names out.
- `history_test`: a copy shares its source's data set; evaluations at A and B and the copy histories leave a
  fingerprint of it unchanged (shown failing with a write of a segment's nominal norm in the hand-over).

**Verified** (clean builds, 0 warnings).  On 075af6db (the data set apart; binary `work/ds3/AZURE2-A`): every
project x 3 bases (tag `tab-ds3A`) against `tab-trial14` 737 / 710 / 737 identical; `ctest -LE slow` 182/182;
`tests/run_tests.sh` 145 / 0; calculations at 1 and 3 threads identical (T1-T3, 13N, natural_neon, energy_shift_sqrt,
shift_sqrt_subpoints/dense, shift_fn, 18O_p_a_thm; three bases, 27/27); fits at 1 thread byte-identical to
`b-trial14` (T2 `--use-lm` and `--use-gradient` 200 calls, 13N, energy_shift_sqrt, screening_6li/fit, 18O_p_a_thm);
the converted-tree gate 483 (456) identical, 38 classified, 0 differ, 0 missing.  The same set on 2bde09d3 (renames,
docs; tag `tab-ds3F`): the same results, but `convert_fixtures` failed (the converter's comment text, fixed in
ef4aebb6); on ef4aebb6: `ctest -LE slow` 182/182 and the gate (`g3ds3G`) as above.  Cycles (perf stat, 1 thread, min
of 6, alternated with `b-trial14`): T2 calculation 645.3 M against 647.8 M, 200-call T2 fit 49.41 G against 50.00 G;
peak RSS T2 fit 305 / 427 MB at 1 / 3 threads against 324 / 452, natural_neon calculation 1465 / 1468 against
1499 / 1506 MB.  One copy after Initialize: T2 11.87 MB in 95 allocations, natural_neon 142.7 MB in 508 (15.06 MB /
236 and 180.6 MB / 779 at the milestone-2 head).

**Open.** pyazr's Python comments (`pyazr/azrfile.py`, `datasets.py`) and two files of the validation branch
(`experiment/src/Build.cpp`, `BuildEffects.cpp`, excepted in check 31) still name the old classes; the
`intEC.dat.sig` header keeps "see EData::ECSignature" (byte identity; change it with the next signature version).
The UDR fold through a point has no project test (only the kernel's reference test).  The external-capture
amplitudes of the set-up move out of the unit with the second pathway branch.

## 2026-11-11 (agent, tree to format 3, second part) — shifts and apparatus projects, the director's decisions (`feature/upgrade-4-tree-format3-2`)

**Done** (`AA2_tree_to_format3_rest.md`; plan and state: `docs/upgrade/format3_tree.md`). 126 of the 185 convertible
project files are format 3 (34 more): the shift projects (`shift_fn`, `shift_groups`, `energy_shift_sqrt`,
`pyazr/data/shift_fn_groups`) and the target, screening, covariance and apparatus projects (`electron_screening`,
`screening_6li`, `target_density`, `data_covariance`, `angular_coefficients`, `beam_profile_kernel`, `hard_sphere`,
`convolution_window`, `convolution_unify/models/bare`, `secondary_gamma/models/grad`). 16 checks rewritten to write or
edit format-3 records (`energy_scale` records for shifts, the target's `screening` and `areal_density`), one generator
(`clone_composite/make_models.py`) writes format 3, `beam_profile_kernel` reads committed extrapolation records.
The director's decisions: the converter refuses an old line with exit key -1 and a type that is not total capture
(old engine's reason, non-zero exit); an old automatic tolerance without windows becomes one window over every energy
with that tolerance (outputs equal, compared); `norm ... fit prior 0 15%` (a centre at 0 or below) is refused by
validation and the run with the new code `E_PRIOR_CENTRE` (`experiment/src/Validate.cpp`). The converter also
refuses four more files the old engine refused and it used to convert into running files (a shiftfn on a follower of a
shift group or of a shift-function group, a line in both kinds of group, densityfit / densitygroup without a thin
target or with a beam profile / UDR). Three slow pyazr tests (`inactive_lines`, `shared_state`, `rmatrix_pole`),
broken by
the first half, edit format-3 records.

**Verified how**: all projects x 3 bases against `tab-trial16` (tag `tab-tree3b`, the new reference): identical but
`param.par` names of 9 converted projects; converted-tree gate (`gtree3b`): 483 (456) identical, 38 classified, 0
differ, 0 missing; spec equivalence 198 files, 185 converted, 183 equal, 2 recorded, 0 differ, 13 not converted (one
more refusal); `ctest -LE slow` 183/183, `ctest -L slow` 10/10; `tests/run_tests.sh` 146 passed, 0 failed; T1-T3, 13N,
natural_neon at 1 and 4 threads byte-identical.

**Open** (the page, §3 and §5): groups 6 (natural targets, with the generators of `doppler_fgm` and
`transmission_yield` and `reaction_rate/models/n12C.azr`, `p12C.azr`), 7 (THM), 8 (T1-T3, relativistic, scattering
lengths, ...); the 28 old-writing generators; the refusal switch (blocked: the unconverted groups, 20 pyazr tests of
the old reader on the originals, the generators).

## 2026-11-11 (director, morning) — the eighth wave, second part merged: the flat working copy; validation that predicts the run's refusals; half of the test tree in format 3

**Merged into `feature/upgrade-4-data`** (`--no-ff`; each rehearsed in a scratch worktree, built and verified once;
the committed trees are those trees): the second milestone of the data structure
(`feature/upgrade-4-data-structure`: a term's settings as one shared record, one flat vector of calculation units,
the plan by index, a copy of the working copy a plain copy), then the validation and pyazr follow-ups
(`feature/upgrade-4-validation`), then the first half of the conversion of the test tree
(`feature/upgrade-4-tree-format3`: 92 of the 185 convertible project files are format 3 under their own names, the
originals in `tests/convert/fixtures/`, 22 `check.sh` scripts edit format-3 records; plan and hand-over page
`docs/upgrade/format3_tree.md`). Their entries are below.

**The merged head as a whole** (build `work/b-trial16`, 0 warnings; reference tag of its outputs `tab-trial16`):
all projects x 3 bases against the previous head (`tab-trial14`, 737 / 710 / 737 files): with the validation alone
(`tab-trial15`) every file identical; with the half-converted tree on top (`tab-trial16`) the 597 / 573 / 597 files
of the projects that are still old identical, and the files that differ belong to the 85 converted projects of
`tests/` and to no other, by the classes of `tools/format3_gate.py` (the records' names in run logs and parameter
files, the sample files of `angle_mapping`, `shifts.out` of the `secondary_gamma` models, the template of
`shift_subpoint_angle` refused without its data files); the outputs of each converted project equal those of the
gate's conversion of its original. The gate (originals laid back, everything converted, against `tab-trial8`):
483 (456 formal) identical, 38 classified, 0 differ, 0 missing. `ctest -LE slow` 183/183; `tests/run_tests.sh` 146
passed, 0 failed; calculations at 1 against 6 threads 12/12 identical; T2 with model data identical in the three
bases; spec equivalence on the originals 198 files, 186 converted, 184 equal, 2 recorded, 0 differ. Milestone 2 by
itself (`tab-trial14` against `tab-trial13`): 724 / 697 / 724 of 726 / 699 / 726 identical, the two run logs of
`sum_lab_frame` differ as text, 11 files new. **Not run on this head:** 1 against 8 threads, ThreadSanitizer,
`ctest -L slow`, Sphinx; Phase 4 as a whole is not yet benchmarked against the Phase 3 head.
The comparison scripts of the director now use the tree's own `tools/runall.sh` and `tools/cls.sh` (they skip
`tests/convert/`); the reference tag for the next branches is `tab-trial16`.

**What the director changed at the merges:** `tests/upgrade/capture_ay_requests_test.cpp` line 153 follows the new
accessor of a term's points (`s.GetPoint(*l.data, k)`); `tests/consistency/file_sizes.tsv` tightened (15 files);
`tests/upgrade/predictions.txt` regenerated on the tree that holds both the validation and the converted projects
(both branches had changed it: the block of `ratio_frame/models/fixed_angle.azr` is new, the SUM data set of
`sum_data` and `sum_free` is described).
**Sent back before the merge:** the rename of an inactive data set left its rows in `<parameterSettings>`, `<mcmc>`
and the prior covariance with the old name (the engine named the parameters of active records only); it now asks the
engine for the names with every record taken as active (`BuildOptions::allRecordsActive`, never set for a run) and
warns, with the lines, of rows it could not rename.
**Checked by hand** (the director, build `work/b-trial15`): a data file with a non-numeric cell in row 7 and a row of
two columns gives two `E_DATA_ROW` diagnostics with file, line, column and the text found, exit 2;
`--describe FILE --record NAME` prints the one record, an unknown name exits 2 with a message; the comparison of the
flat working copy with the previous head flags exactly `sum_data`, `sum_free` (run logs, text) and the new files.

**Decisions of the director.**
- A copy of the working copy keeps one explicit step after the plain copy: its term points drop the set-up's
  external-capture amplitudes. They belong to the external-capture set-up, which the second pathway branch moves to
  the immutable tables; until then the step stays and is named in `Clone`.
- The final state of a photon pair with external capture is decided by one function
  (`experiment::ExternalCaptureFinalStates`), which the engine's set-up (`CNuc::ParseExternalCapture`) and the
  validation both call. The agent edited `src/CNuc.cpp` for it beyond its brief's limit: accepted, the messages and
  every output are unchanged.
- A measured uncertainty ≤ 0 is a **warning** of the validation (`W_DATA_ERROR`), not an error: the run accepts such
  rows (χ² skips a zero and squares a negative, `chi_squared.rst`). Whether the run should refuse a negative
  uncertainty is a question for the maintainers (three committed test data files have such rows: `13N`
  `baumann_1618.dat` negative, `16N` `bda.dat` and `18O_p_a_thm` `lc723_band_mid.dat` zero).
- On the cases the conversion found that format 3 cannot write (`format3_tree.md` §4): (1) an old line with exit
  key −1 and a data type that is not total capture is **refused by the converter**, as the old engine refused it (it
  writes `[total]` today: a conversion must never turn a refused file into one that runs); (2) automatic windows
  together with explicit ranges get **no grammar**; the converter drops the ranges with a note; (3) an automatic
  tolerance without windows must not be dropped silently: the converter writes the one window over every point that
  gives the old outputs, or a note; (4) a fitted normalization with a prior centred at 0 is **refused** by the
  validation and the run (the old engine refused a centre of 0); (5) a fixed normalization group whose follower's
  own line said "free" has no counterpart and needs none. (1), (3), (4) are in the brief of the tree's second half.
- The two expected texts `tests/ratio_to_rutherford/expected/story_describe.txt` and
  `tests/units_scales/expected/story_describe.txt` changed with the wording of an open energy window: accepted.

**Defects and debts found, not yet fixed** (owner in brackets):
- A data row that does not parse makes the run drop the whole data set with a WARNING ("Could Not Fill Segment") and
  go on; the validation refuses it. The run should refuse it too [the deletion of the old reading path].
- Refusals the validation does not predict yet (survey in the validation entry: 2 files of the converted tree, 26
  generated models): the `<levels>`-row checks of `CNuc::Fill`, kinematic limits of a lab angle, combinations of
  target effects, a covariance matrix that is not positive definite, two fitted rows of one covariance file from
  overlapping windows [a second validation round, after the old path is deleted: the rules then move to `experiment/`
  once instead of being written twice].
- `--describe` words the angle window "theta_lab (data converted to c.m.)" also for a file whose angle column is
  `theta_cm` (`experiment/src/Describe.cpp`) [second validation round].
- The names of the target-effect parameters (`.screening`, `.areal_density`, `.doppler`) are written in the data
  classes and again in `experiment::ParameterNamesOf`; `ESegment::Fill` compares the covariance variance with 0
  itself [milestone 3 of the data structure or the round after it].
- The warning of `Project.rename` for rows it could not rename matches the old name as plain text and does not
  cover the names line of the prior-covariance file [second validation round].
- Possible dependence of the external-capture memo on the history of a copy (`ShiftedECAmplitude::levelEnergy`;
  inferred by the agent from the code, not tested) [second pathway branch: a test with a fitted level energy of the
  final state].
- The covariance groups' frozen model and `NodeTableCache` are shared and mutable between copies by design
  [milestone 3: named with the reason they are safe, or reported].
- Still open: the yield on a finite layer (Phase 5), the decay-position defect of the external-capture code (second
  pathway branch), `ScatteringLength::ZeroEnergyShift` for a hybrid pair, an exception inside the external-capture
  path under OpenMP (Phases 5–6).

**For the PI to see** (behaviour a user can meet; additions):
- `AZURE2 --validate` now reports, before any calculation: every data row that does not parse or is too short (file,
  line, column, the text found), a missing or unreadable covariance file or one that does not fit its data set, a
  photon pair whose final state is not a level or is shared with another pair; uncertainties ≤ 0 as warnings.
- `AZURE2 --describe FILE --record NAME`; an open energy window reads "all energies".
- Python: `session.datasets`, `active_datasets`, `dataset_chi2` work on a format-3 project (keyed by the record
  names); `pyazr.records_from_csv` declares data sets from a table; pyazr is 2.31.0.
- The file of the natural decomposition (`AZUREOut_sample_seg=N`) written from a copy of the data (MCMC, fit
  snapshots) now has the columns of the nuclei; before, a copy wrote the total only.
- An S-factor data set's output columns hold σ and S (one formula for the conversion).
- Memory of a fit: T2 at 1 thread 328 MB (was 355), `natural_neon` 1.50 GB (was 1.55); one copy of the working copy
  15 MB for T2 (was 18.6).

**Checks run after the merge** (head c285a312, build `work/b-trial16`; `work/director/extra_trial16.out`):
- Calculations at 1 against 8 threads: 12 of 12 identical.
- T2 benchmark against the head after Phase 3 (28a33359), three alternated repetitions, median (`work/bench-p4w8.tsv`):
  calculation 1.18 s against 1.16 s at 1 thread, 0.39 s against 0.41 s at 6; fit of 200 calls 70.9 s against 72.9 s
  at 1 thread, 38.7 s against 38.5 s at 6; chi-squared 691.608 in every run. Peak memory: calculation 57 MB against
  63 MB; fit at 1 thread 240 MB against 276 MB; fit at 6 threads 570 MB against 505 MB (13 % more: cause not looked
  for; to be measured again after milestone 3, where the copies share the data set).
- `ctest -L slow`: 7 of 10 pass. `pyazr_inactive_lines`, `pyazr_shared_state` and `pyazr_rmatrix_pole` **fail on
  this head**: they edit or load project files that the first tree branch converted to format 3 (the agent and the
  director ran `-LE slow` only). Given to the agent that converts the rest of the tree; `ctest -L slow` is now
  part of the director's verification (`work/director/full_verify.sh`).

**Running now** (Opus agents; briefs in `/data0/skowrons/AZURE2-upgrade/work/briefs/`):
- Milestone 3 of the data structure (`X3_dataset_workingcopy.md`; `feature/upgrade-4-data-structure-3`): an
  immutable `DataSet` and a small `WorkingCopy` apart, the unit slimmed, the names `EData`, `ESegment`, `EPoint` and
  their files deleted.
- The rest of the test tree and of the examples converted to format 3 (`AA2_tree_to_format3_rest.md`;
  `feature/upgrade-4-tree-format3-2`; from `docs/upgrade/format3_tree.md`: the shift tests, targets and apparatus,
  natural targets, THM, T1–T3 and the rest; the 29 model generators; the converter's refusals decided above; the
  switch that makes the engine refuse an old file, last and only if nothing loads one).

**Then.** The second pathway branch (external-capture amplitudes at set-up with the decay-position fix, the trees
deleted); the engine refuses old files, the old reading path and `pyazr/azrfile.py` deleted (inventory:
`docs/upgrade/pyazr_on_records.md`); a second validation round; components as leaves, amplitude filters, `predict`;
the recipe "a new kind of data"; the rest of `CNuc`; Phase 5 (the yield of a finite layer first).

## 2026-11-11 (agent, tree to format 3) — the test tree converted, first half (`feature/upgrade-4-tree-format3`)

**Done** (`AA_tree_to_format3.md`; plan and state: `docs/upgrade/format3_tree.md`). 92 of the 185 convertible
project files of `tests/` and `examples/` are format 3, in place under their own names; the originals are in
`tests/convert/fixtures/<path>` and stay the converter's fixtures (the gate, `convert_spec_equivalence` and
`convert_fixtures` read them there; `tools/convert_tree.sh` lays them back before converting). Groups done: 1 (81
self-contained projects and 7 examples), 2 (13N family and 7Li_p_ay), 3 (thick_target_resonance), 4a (norm_groups,
15N_p_a). 22 checks rewritten to edit format-3 records (or write them); `segline_formats` and parts of
`total_capture_code`, `prior_centre`, `ec_final_state_ambiguous` test old lines through `AZURE2 --convert`. pyazr tests
of the old file class (14) read the originals (`tests/pyazr/old_tree.py`). Tools: `tools/format3_tree_convert.sh`
(convert files in place, originals to the fixtures).

**Verified how**: per group, the touched directories' regression projects and checks, the rewritten checks against
the base checks on the base build (same output, differences listed in the page), `scripts/check.sh --quick`
(ctest -LE slow 182/182). At the last commit: converted-tree gate (`gtree3a`) 483 (456) identical, 38 classified, 0
differ, 0 missing; spec equivalence 197 files, 185 converted, 183 equal, 2 recorded, 0 differ; `convert_fixtures` 86
checks; T1–T3 and 13N at 1 and 4 threads identical; reference tag of this tree: `tab-tree3`
(page §3a); `tests/run_tests.sh` whole: see the report.

**Open** (the page, §3–§5): groups 4b (shifts), 5 (targets and apparatus), 6 (natural targets), 7 (THM, the largest:
about 20 scripts and 12 pyazr tests), 8 (T1–T3, relativistic, scattering lengths, ...); the refusal switch (blocked:
pyazr tests load old files through `pyazr/azrfile.py`; the old `make_models.py` generators). For the director: the
converter writes an old exit −1 / data-type mismatch as `[total]` (the old engine refused it); format 3 cannot state
automatic windows plus explicit ranges, nor an automatic tolerance without windows (the converter drops the latter
without a note); `norm … prior 0 15%` is accepted (a zero-width prior, not applied) where the old engine refused a
centre 0.

## 2026-11-11 (agent, validation follow-ups) — validation reaches the run's refusals; pyazr's session on format 3 (`feature/upgrade-4-validation`)

**Done.**
- `--validate` / `pyazr.validate` / every `Project` edit (a build that uses the data; the engine's own build in
  `Config` does not, its set-up refuses at run time as before) now report what a run refuses: every row of a data file
  that does not parse (`E_DATA_ROW`, line, column, text found; `experiment::CheckDataRows` beside
  `ReadMeasurementRows`), a row with fewer columns than the roles (`E_DATA_FILE`, every row not the first), an
  angular-coefficient row without a Legendre order (`E_DATA_ROW`), a non-positive uncertainty in the window
  (`W_DATA_ERROR`, a warning: the chi^2 skips a 0); the covariance file (`E_COVARIANCE_FILE`: missing, unreadable,
  row counts, data file not listed, a variance <= 0 in the window; its one reader moved to
  `experiment/CovarianceFile`, the engine's `DataCovariance::Read` calls it); the final state of external capture
  (`E_FINAL_STATE_NOT_LEVEL`, `E_FINAL_STATE_SHARED`: one rule, `experiment::ExternalCaptureFinalStates`, called by
  `CNuc::ParseExternalCapture` with its messages unchanged and by `Schema::Build` under `BuildOptions::setUpChecks`).
- `AZURE2 --describe FILE --record NAME`; an open energy window reads "all energies (no energy_window)" (the expected
  describe texts of `ratio_to_rutherford` and `units_scales` follow).
- One function names an abundance parameter (`experiment::AbundanceParameterName`; `CNuc::AbundanceParamName`'s body
  calls it); `experiment::ParameterNamesOf` / `_azure2.parameter_names`; `Project.rename` maps the engine's names.
- pyazr 2.31.0: `azure2.datasets`, `active_datasets`, `extrapolations`, `dataset_chi2` on a format-3 session from the
  engine's description (`pyazr/described.py`); `Project.set` keeps comments between repeats;
  `pyazr.records_from_csv`; `pyazr.describe(path, record=)`; the inventory of `azrfile.py` uses:
  `docs/upgrade/pyazr_on_records.md`.
- Tests: `tests/validate_run_refusals` (format-3 fixtures `*.f3`), `pyazr_session_records`, build_test cases of the
  five new codes, `convert_spec_equivalence` reads with the engine's options.

**Verified** (build `work/b-valid`, 0 warnings, 5 threads): all projects x 3 bases against `tab-trial13` (tag
`tab-valid1`): 726 / 699 / 726 identical; converted-tree gate (`REF=tab-trial8`, tag `gvalid1`): 483 (456) identical,
38 classified, 0 differ, 0 missing, no converted file's bytes moved; `tests/run_tests.sh` 146 passed, 0 failed (the new
`validate_run_refusals`); consistency test: 28 flags, 80 diagnostic codes; `--validate` of T2's conversion 4.6 ms
against 4.0 ms on `b-trial13` (minimum of 20); `ctest -LE slow` 183/183; 13N, T1, T2 and the final-state models at 4
threads: 47 output files byte-identical to the 1-thread run.

**Corrections after the director's review.** `Project.rename` asks the engine for the names with every record taken as
active (`BuildOptions::allRecordsActive`, `parameter_names(text, all_active)`), so the rows of an inactive data set
follow and apply once it is switched on; rows that still contain the old name (no parameter the engine gives) are
listed in a warning.  `Context::Warn` sets `hasFound`: the JSON of `W_DATA_FILE_NOT_READ` carries the path (no
committed expected file held it).

**Open.** Refusals validation does not predict (survey of the converted tree and of the check.sh models): `<levels>`
rows whose level ID has channels of two nuclei and the other `CNuc::Fill` row checks; a lab angle beyond the
kinematic maximum (`relativistic/models/li7pn_edge_rel.azr`); two fitted rows of one covariance file from overlapping
windows; a covariance matrix not positive definite. `ESegment::Fill` still compares a covariance variance with 0 by
itself (the validation reads the variance through the shared reader); EData's `EffectParameterName` should give its
names through experiment (`ParameterNamesOf` copies its three suffixes). An unparseable data row: the run drops the
data set with a WARNING ("Could Not Fill Segment"), validation refuses it.

## 2026-11-11 (agent, brief X2, milestone 2, rest) — the flat working copy (`feature/upgrade-4-data-structure`)

**Done** (four commits after a724d6ba: cfe7a91d, 044928ed, b0759ed8, 60bcb381).
- A point has no pointer to its data set (`EPoint::parentData_` gone).  The data set is an argument of
  `EPoint::Calculate`, `CalculateThinTarget`, `MarkFailed`, `PlanMapped`, `Initialize`, the folds, `ESegment`'s shift
  step and `HasSubPointIntegration`, `azure2::EvaluateModelPoint` / `EvaluateModelSegment` / `ModelPointBackward` /
  `PointResidual`, `CoulombReferenceValue`, the gradient helpers of `AZUREGrad.cpp`; `CalculatePointsInParallel` is a
  member.  Null for a point outside any data set (reaction-rate points).
- One flat vector of calculation units in the working copy (`EData::Units`, `include/UnitRange.h`,
  `src/WorkingCopyUnits.cpp`): a term's points and a point's sub-points are ranges (`ESegment` TermUnits first/count,
  `EPoint` sub_first_/sub_count_); `GetPoints`, `GetPoint`, `TermPoints`, `GetTermPoint`, `GetSubPoint`,
  `GetSubPoints` keep their names and are views taking the data set.  Sub-points given at set-up
  (`EData::AddSubPoint`, held per point) are placed by `CommitSubPoints` at the end of `SetUpTargetEffects`, which
  lays the vector out (term points in data order, then the sub-points in the same order; each unit moved once into
  a vector of the exact size); a re-placed grid (`ReplaceSubPointGrids`) is built in the parallel loop beside the
  units and placed single-threaded (`PlaceSubPointGrids`): in the place of the old one when it has as many nodes,
  else by a new layout.  The points keep their indices after the set-up.  A unit names its term's settings by a plain
  `const TermSettings *` (the term owns the shared pointer; a point outside a data set owns its own).
- The plan names units by index (`PointRef` is an int; `EData::PlanUnits` gives the data order).
- `EData::Clone` is the copy constructor; then the copy's term points drop the set-up's external-capture amplitudes,
  as before (clearing the sub-points' too moved the 13N fit in the 11th digit: kept as it was).  The objective
  counter is `EvaluationCount` (copyable).  THM group records, the counter and the inline sub-point accessors moved to
  `ThmGroup.h`, `EvaluationCount.h`, `SubPointAccess.h` (no data-class file is longer than at a724d6ba).
- `tests/upgrade/history_test.cpp`: a copy of a copy shares no mutable state with its source (shown failing with
  copy-on-write switched off in `SharedTable`).

**Verified** on 60bcb381 (binary `work/dstruct/x2/AZURE2-s5`, build 0 warnings): every project x 3 bases (tag
`tab-x2s5`) against `tab-trial12`: 724 / 697 / 724 identical of 737 / 710 / 737, the rest the run logs of
`sum_data` / `sum_free` and the new files of those two and of `fixed_angle` (as at a724d6ba); against `tab-dsm2a`:
737 / 710 / 737 identical; `ctest -LE slow` 177/177 (the history, poison and `AZURE_PLAN_CHECK=1` tests among them);
`tests/run_tests.sh` 145 passed, 0 failed; calculations at 1 and 3 threads identical (T1-T3, 13N, natural_neon,
energy_shift_sqrt, shift_sqrt_subpoints/dense, shift_fn, 18O_p_a_thm, three bases: 27/27); the converted-tree gate
(`gx2s5`): 483 (456) identical, 38 classified, 0 differ, 0 missing; fits at 1 thread byte-identical to
`b-trial12` (T2 `--use-lm` 200 calls, T2 `--use-gradient` 200 calls, 13N, energy_shift_sqrt, screening_6li/fit with
adaptive grids, 18O_p_a_thm).  Cycles (perf stat, 1 thread, min of 6, alternated with `b-trial12`): T2 calculation
669.4 M against 674.0 M (-0.7 %), 200-call T2 fit 50.33 G against 51.15 G (-1.6 %); peak RSS of the T2 fit 328 / 445
MB at 1 / 3 threads against 355 / 447, natural_neon calculation 1499 / 1504 MB against 1553 / 1557 (a first layout,
which placed the set-up's sub-points through a list and a compaction, was 2.7 % slower on the T2 calculation and
up to 23 % higher in RSS: 60bcb381 places each unit once).  Bytes allocated by one copy
(`work/dstruct/x2/copybytes`, after Initialize): T2 18.61 MB in 549 allocations before, 15.06 MB in 236 after;
natural_neon 223.2 MB in 3901 before, 180.6 MB in 779 after; the copy is now nearly all the unit vector (24,628
units of 608 bytes for T2, 295,849 for natural_neon).

**Open.** Milestone 3: `DataSet` and `WorkingCopy` apart (an `ESegment` still holds both its settings and its
working state: norm, shifts, frames, shift tables), the `EData &` arguments become the working copy, the old names
and files go.  The set-up's external-capture amplitudes still live in the points (a copy drops those of its term
points by an explicit loop in `Clone`): they belong to the external-capture set-up of the second pathway branch.
The unit is 608 bytes; slimming it is what makes a copy cheaper now.

## 2026-11-11 (agent, brief X, milestone 2, first part) — settings out of the points, three debts (`feature/upgrade-4-data-structure`)

**Done.**
- The settings of a term (pairs, observable flags, J/l/s, UPOS, cross-section component, THM attachments) are one
  record, `TermSettings` (`include/TermSettings.h`), built by the segment at set-up and shared, immutable, by the
  term's points, their sub-points and every copy; `EPoint` lost about twenty copied members and their setters.  The raw
  rows pointer of a point is gone: a point keeps its row index, Legendre order, nominal energy and angle; the
  measured value, error and record come from the measurement (`ESegment::LabValue`, `LabError`, `FileRecord`).
- The S factor has one formula (`physics::sfactor_conversion`, E sigma for a neutral pair); the AZUREOut columns of an
  S-factor data set or extrapolation hold sigma (b) and S (MeV b) (`SigmaAndS` in `EDataOutputFiles.cpp`);
  `tests/units_scales/check.sh` reads them so.
- `tests/ratio_frame` (g): a term at a fixed angle, c.m. (not converted) and lab (converted; the c.m. angle by the
  closed form in `make_models.py`).  `tests/sum_lab_frame/models/data/sum.dat` committed, so the SUM data segment
  fills in a run of the project tree (`check.sh` copies models/data with `cp -n`).

**Verified** (build `work/b-dstruct`, 0 warnings; binary `work/dstruct/AZURE2-m2a`): every project x 3 bases against
`tab-trial12` (tag `tab-dsm2a`): all 726 / 699 / 726 reference files identical but the run logs of `sum_data` and
`sum_free` (their segment fills now), plus the new files of those two and of `fixed_angle`; `ctest -LE slow` 176/177,
the one failure `upgrade_predictions` (the committed table lacked the new and now-filled projects; regenerated, nine
added lines, then passed); `tests/run_tests.sh` 145 / 0; the gate (`gdsm2a`): 483 (456) identical, 38 classified, 0
differ, 0 missing; calculations at 1 and 3 threads identical (T1-T3, 13N, natural_neon, energy_shift_sqrt, three
bases, 18/18); fits at 1 thread byte-identical to `b-trial12` (T2 `--use-lm` 200 calls, T2 `--use-gradient` 200
calls, 13N, energy_shift_sqrt).  Speed (min of 6, 1 thread, alternated with `b-trial12`): T2 calculation 659.9 M
cycles against 663.7 M, 200-call T2 fit 49.56 G against 50.48 G; peak RSS of the T2 fit 331 / 412 MB at 1 / 3 threads
against 355 / 447.

**Open (rest of milestone 2).** The flat unit vector (sub-points as ranges), the plan by unit index, the
parent-data pointer of a point (folds, model pass and gradient then take the data set), `Clone` as a plain copy.  A
term's settings are named by a shared pointer, not by an index: the matrix functions ask the point for them and
their bodies may not change while the pathway work runs.

## 2026-11-11 (director, early morning) — the eighth wave, first part merged: the project as records in Python; capture analyzing power at set-up

**Merged into `feature/upgrade-4-data`** (`--no-ff`, in this order; rehearsed together in one scratch worktree, built
and verified once; the committed tree is that tree): `pyazr.Project`, a format-3 project as named records in Python
(`feature/upgrade-4-pyazr-project`); the capture analyzing-power terms built at set-up, with the refusal of two
photon pairs that have one final state (`feature/upgrade-4-capture-ay`). Their entries are below.

**The merged head as a whole** (build `work/b-trial13`, 0 warnings; reference tag of its outputs `tab-trial13`): all
projects × 3 bases identical to the seventh wave's head (726 + 699 + 726 output files, none new or missing);
`ctest -LE slow` 182/182; `tests/run_tests.sh` 145 passed, 0 failed; calculations at 1 and 6 threads 12/12 identical;
T2 with model data identical in the three bases; spec equivalence: 197 files, 185 converted, 183 equal, 2 recorded,
0 differ; the converted tree (exact grids, menu 1) against the fourth-wave reference: 483 (456 formal) identical, 38
classified, 0 differ, 0 missing; consistency test: 28 flags, 75 diagnostic codes. Not run on this head: 1 against 8
threads, ThreadSanitizer, `ctest -L slow`, the Sphinx build (not installed here).

**What the director changed at the merge:** `project_format3.md` §4.1 names the call `add_normalization`.
**Checked by hand** (the director, on the agent's build): one `set(norm=Nuisance(...))` on a converted project changes
exactly one line of the file; a wrong unit, an unknown reference and a value outside a list are refused with the
engine's codes (`E_UNIT_NOT_ALLOWED`, `E_UNKNOWN_REFERENCE`, `E_NOT_IN_LIST`); `remove` of a target that two data sets
use is refused naming both; `rename` updates the references; the generated help of `add_dataset` lists the fields.

**Decisions of the director.**
- pyazr: the calls are named after the record kinds (`add_normalization`; no `add_norm` alias). A calculation is run
  through `Project.session()`, which saves and returns the existing calculation object; no second calculation
  interface.
- A project in which two photon pairs resolve to the same final level is **refused at set-up** (it was read without
  a message and gave a wrong capture cross section, a factor 2.2 at the first point of the 13N test with its photon
  pair repeated). The other case, two final states on one photon pair, was refused already.
- A unit of the working copy keeps a direct reference to its term's immutable settings (`TermSettings`) as the end
  state of the data structure; the term owns them, the unit does not.

**Defects and debts found, not yet fixed** (owner in brackets):
- `--validate` / `pyazr.validate` pass a project whose covariance file is missing, whose data rows do not parse
  beyond the first, or that has the duplicated photon pair above; the run then refuses it [validation follow-ups,
  launched now].
- `azure2.datasets`, `active_datasets`, `dataset_chi2` do not work on a format-3 session (they read the old data
  block through `pyazr/azrfile.py`) [validation follow-ups].
- The parameter name `<target>.isotope.<isotope>` is built in two places; `Project.rename` matches parameter names by
  a pattern in Python; `Project.set` drops whole-line comments between repeats of a replaced field; `--describe`
  prints `1.79769e+308` for an open window; `--describe --record NAME` and `records_from_csv` are not built
  [validation follow-ups].
- The capture-A_y table is built at the configured maximum order (before, a non-differential capture-A_y point
  evaluated first would have left the table empty; no project does that) [none: recorded].
- Still open from the seventh wave: the yield on a finite layer (Phase 5), the decay-position defect of the
  external-capture code in the data classes (second pathway branch), `ScatteringLength::ZeroEnergyShift` for a
  hybrid pair, an exception inside the external-capture path under OpenMP (Phases 5–6).

**For the PI to see** (behaviour a user can meet; additions):
- Python: `p = pyazr.Project.open("x.azr")`; `p.add_dataset(...)`, `p["name"].set(...)`, `p.rename`, `p.remove`,
  `p.validate()`, `p.describe("name")`, `p.save()`, `p.session().calculate()`. Every refusal carries the engine's
  diagnostic code and the allowed values. Recipe: `docs/source/user_guide/project_from_python.rst`.
- A project with the same final state on two photon pairs no longer runs: the message names both pairs and the
  level. Old results from such a project were wrong for the first pair's capture.
- `pyazr.save_fit(norms="fitted")` works on a format-3 project; on an old-format project with a changed potential it
  is refused with the conversion command.

**Running now** (Opus agents; briefs in `/data0/skowrons/AZURE2-upgrade/work/briefs/`):
- The flat working copy of the data, second half of milestone 2 (`X2_flat_working_copy.md`,
  `docs/upgrade/data_structure.md`; `feature/upgrade-4-data-structure`). Its first half is on the branch, verified by
  its agent, not merged: a term's settings as one shared immutable record, one S-factor formula, the output columns
  of an S-factor data set holding σ and S, a test of a term's fixed-angle frame.
- The test tree and the examples converted to format 3, `check.sh` scripts included (`AA_tree_to_format3.md`;
  `feature/upgrade-4-tree-format3`): the first half of "the engine reads only format 3". The refusal of old files is
  switched on only when nothing loads one any more.
- Validation that catches what a run would refuse, and pyazr's session on a format-3 project without the old file
  class (`AB_validation_pyazr_followups.md`; `feature/upgrade-4-validation`).

**Then.** Milestone 3 of the data structure (the old class names and files deleted); the second pathway branch
(external-capture amplitudes at set-up, the trees deleted); the old reading path and `pyazr/azrfile.py` deleted;
`predict`, the recipe "a new kind of data"; the rest of `CNuc`; Phase 5 (targets: the yield on a finite layer
first; the lab-frame comparison of shifted data); Phase 6; the GUI on the new format; `FOR_THE_OWNER.md` and
`ROADMAP.md` for the meeting with the owner; `origin/dev` merged only when the reconstruction is complete.

## 2026-11-11 (agent, pathway table) — capture A_y at set-up (`feature/upgrade-4-capture-ay`)

**Done** (`cnuc_split.md` §11.6 step 6, §11.12).
- `RequestCaptureAyTerms` (`include/CaptureAyRequests.h`, `src/CaptureAyRequests.cpp`) names every term (own and
  leaf) of an analyzing-power data set or extrapolation with a photon exit, after the data are filled and before
  `CNuc::Initialize`, in `AZUREMain`, `runMCMC` and both set-up paths of `AZUREAPI.cpp`.  `BuildPathwayTable` builds
  the requested terms in the table's order (the table does not depend on the order of the requests), none for a beta
  entrance.
- `GenMatrixFunc::CalculateCaptureAnalyzingPower` and the capture-A_y block of `PointTBar` read the table; a point
  whose terms were not requested is refused with a message naming the data set.  `CNuc::CalcCaptureAnalyzingPower`,
  its static mutex and the trees' capture-A_y members are deleted.  The matrix functions read no tree any more.
- Test `upgrade_capture_ay_requests` (every project in data and extrapolation mode: terms present for the capture-A_y
  data sets `13N_capture_ay` and `capture_ay_coupling/models/twolevel.azr`, none for beta entrances, order independence,
  refusal without the requests).

**Verified**: (build `work/b-pathtab`, 0 warnings, 5 threads): 183 x 3 against `tab-trial12` (tag `tab-cayfinal`): 726
/ 699 /
726 identical; `ctest -LE slow` 178/178; `tests/run_tests.sh` 145 passed, 0 failed; converted-tree gate (`gcayfinal`):
exit 0, 483 (456) identical, 38 classified, 0 differ, 0 missing; 1 against 3 threads in three bases on T1-T3, 13N,
13N_capture_ay, beta_delayed_alpha, plus capture_ay_coupling and upos_pair_index through their check.sh: identical;
fits at 1 thread against `b-trial12`, every output file identical: T2 `--use-lm` and `--use-gradient` (200 calls),
13N_capture_ay MIGRAD, `--use-gradient` and `--use-lm`; capture_ay_coupling's twolevel.azr extrapolated in three
bases identical (deuteron.azr is refused, as before); `AZURE_GRAD_CHECK` with the adjoint on 13N_capture_ay: 14
components, largest relative difference to central differences 1.2e-9, the same numbers as the base build.  Cycles
(1 thread, `--no-coulomb-table`, 5 runs alternated with `b-trial12`, minimum): T2 calculation 3.310e9 against 3.303e9
(+0.2 %), 13N_capture_ay calculation 6.765e9 against 6.762e9 (+0.05 %); instructions +0.002 % and -0.004 %.

**Two photon pairs with one final state are refused** (the director's decision): `ParseExternalCapture` stops with an
ERROR naming both pair keys and the level (`FinalStateOfAnotherPair`, `src/CNucExternalCapture.cpp`); before, the
second pair took the level and capture into the first ran without its external part and without a message
(1.95e-9 against 4.32e-9).  Format 3 goes through the same function.  `tests/ec_final_state_ambiguous/check.sh` (d) the
two-pair model, (e) its format-3 conversion (converted, refused when run); shown failing with the check disabled.
Verified on that commit: 183 x 3 against `tab-trial12` (tag `tab-twopair`) 726 / 699 / 726 identical, no project
newly refused; `ctest -LE slow` 178/178; `tests/run_tests.sh` 145 passed, 0 failed (the cases are in the existing
project, so the count stays 145); gate (`gtwopair`) 483 (456) identical, 38 classified, 0 differ, 0 missing.  Two
final states cannot share one photon pair
(a second level in the 1 keV window is refused already).

**Open.** Steps 3 and 8 (after the data classes are replaced).  Follow-up: the format-3 validation could refuse two
photon pairs with one final state earlier, with a diagnostic code (not added: the schema and validation of
`experiment/` are being edited by another agent); today the converter writes such a file and the run refuses it.

## 2026-11-11 (agent, brief Z) — `pyazr.Project`: a format-3 project as records in Python (`feature/upgrade-4-pyazr-project`)

**Done.**
- `pyazr.Project` (`pyazr/project.py`): `open` (a format-1/2 file refused with the conversion command), `new(path,
  nuclear_part=file)` (another format-3 file's record-free sections: `<config>`, `<levels>`, `<thm>`, ...),
  `records()`, `p[name]` (a `Record`: fields as attributes and as a dict; its own kind is `record_kind`, since `kind`
  is a target's field), `add(kind, name, **fields)` and one `add_<kind>` per kind of the schema, generated at import
  from `_azure2.schema()` (signature and help: every field with title, grammar, condition, units, doc), `set`,
  `rename` (header, every reference the validator resolves, parameter names in `<parameterSettings>`/`<mcmc>` and the
  prior-covariance file), `remove` (refused while referenced, naming the referrers), `validate`, `describe(name)`,
  `save`, `session(**options)` (saves, returns `pyazr.azure2` on the file).  The project is the engine's `Document`;
  an edit replaces only its fields' lines.  Every edit is validated by the engine's one validator on the whole project
  as it would stand; an edit that adds an error raises `pyazr.FormatError` (`.diagnostics`, `.codes`: the engine's)
  and leaves the project unchanged.  `pyazr.Nuisance` and `(value, unit)` tuples serialise to the file text.
- Engine: `experiment/ValueText.h` (the one nuisance writer, the converter's too; the edit writing a fitted value into
  its nuisance in its unit), `BuiltProject::parameterSources`, `ValidationResult::references`, describe of one record,
  `io::WriteTokens`; bindings in `api/src/pybind_project.cpp` (the `Document` binding moved there verbatim first).
- Debts fixed: `save_fit(norms="fitted")` of a format-3 project writes free norms, shift terms, shift-function
  parameters and isotope abundances into their records (it raised; and it added an `<abundances>` block to a format-3
  file); `save_fit` of a format-1/2 project after `set_nuclear_potential` is refused with the conversion command;
  `set_nuclear_potential(matching_radius=, step=)`; `io::WriteRecords` aligns a record's values past its longest key
  (fixtures T2.f3, T3.f3, beam_profile_kernel.f3 moved by whitespace only); `tools/classify.py` never counts an
  unreadable file as identical.
- Recipe: `docs/source/user_guide/project_from_python.rst`, `.claude/skills/azure2-setup/references/
  project-from-python.md`, code `pyazr/examples/project_from_python.py` run by `pyazr_project_recipe`.  pyazr 2.30.0.

**Verified** (build `work/b-pyproj`, 0 warnings): 183 × 3 against `tab-trial12` (tag `tab-pyproj1`): 726 / 699 / 726
identical; converted-tree gate (`REF=tab-trial8`, tag `gpyproj1`): 483 (456) identical, 38 classified, 0 differ, 0
missing; `ctest -LE slow` 180/180; `tests/run_tests.sh` 145 passed, 0 failed; T1–T3, 13N at 1 and 4 threads in three
bases: every output file identical (run.log differs in Coulomb-table lookup counts and progress bars only); one edit
on T2 costs 2.7 ms (one validation with build 0.8 ms).  New tests (each shown failing once with an injected error):
`pyazr_project` (37 checks), `pyazr_save_fit_format3`, `pyazr_project_recipe`, `tools_classify`, a case of
`io_record_reader`.

**Open.** `--describe --record NAME` on the command line (the engine function takes the name; only pyazr uses it);
`records_from_csv` (§6.5); `add_norm` of §4.1 is `add_normalization` (the kind's name); the parameter-name form
`<target>.isotope.<isotope>` is now made in `CNuc::AbundanceParamName` and in `Build` (its source), a second copy;
whole-line comments between repeats of a field are dropped when `set` replaces that field; Sphinx not run (not
installed here).  Found, not fixed: `--validate` passes a project whose covariance file is missing or whose data rows
do not parse (it checks the first row's column count only); `azure2.dataset_chi2`/`datasets` read `<segmentsData>`
and do not work on a format-3 session.

## 2026-11-11 (director, night) — the seventh wave merged: units, the measurement held once, spin–orbit, the pathway table

**Merged into `feature/upgrade-4-data`** (`--no-ff`, in this order; rehearsed together in one scratch worktree, built
and verified once; the committed tree is that tree): units of every number, the S factor as a quantity and yields per
charge (`feature/upgrade-4-units`); the first milestone of the new data structure, what was measured exists once
(`feature/upgrade-4-data-structure`); the spin–orbit term of the hybrid potential (`feature/upgrade-4-spin-orbit`);
the reaction pathways as one immutable table, first branch (`feature/upgrade-4-pathway-table-1`). Their entries are
below.

**The merged head as a whole** (build `work/b-trial12`, 0 warnings): every project × 3 bases: the 704 + 677 + 704
output files of the reference tag `tab-trial8` identical, and only the 22 files per basis of the four new project
files (the units story, the spin–orbit project and its two models) are new;
`ctest -LE slow` 177/177; `tests/run_tests.sh` 144 passed / 1 failed in the verification: the units project compares
its `--describe` text, which gained the "rows read, kept" line of the data structure — the expected text is corrected
in the data-structure merge and the check passes alone (the whole script is being run again; its result is in the
next entry if it is not 145 / 0); calculations at 1 and 6 threads 12/12 identical; T2 with a `<modelData>` section
identical in the three bases; spec equivalence: 197 files, 185 converted, 183 equal, 2 recorded, 0 differ (12 not
converted: the four new files are format 3 already); the converted tree: 483 (456 formal) identical, 38 classified,
0 differ, 0 missing.

**What the director changed at the merges** (each in its merge commit's message):
- Units × data structure: the units branch converts the columns of a data row to the engine's units where the row is
  read; the data-structure branch moved that reading into `experiment/Measurement.cpp`. `ReadMeasurementRows` now takes
  the data set's `ColumnConversion` and applies it to the line before the order check, the row window and the
  covariance file; `ESegment::Fill` and `Describe.cpp` pass it. `src/ModelPass.cpp` holds the three factor kinds (term
  frame, ratio frame, S factor).
- pyazr is **2.29.0** (units 2.27.0, the data structure's `describe` fields 2.28.0, spin–orbit 2.29.0; the pathway
  table changed no interface).
- Checked by hand: S = 1 keV b of p + ¹²C at 100 keV gives 1.37383423e-10 b by the formula and by pyazr; the sign of
  the spin–orbit term (with V_so > 0 the j = l + ½ partner is the attractive one, FRESCO's convention).

**Decisions of the director.**
- Spin–orbit: the agent's extra refusal code `E_SO_NO_TERM` and the rule that `E_SO_THM` / `E_SO_RELATIVISTIC` take the
  place of the general codes when the record has the term are accepted; pyazr refuses a term on the default pair.
- Pathway table: a point whose exit pair no pathway of its entrance reaches is refused with a message (before, the
  last decay was read); accepted. The correction of `azure2.ec_integrals` in the formal basis is accepted as a fix.
- Data structure: milestone 1 keeps a raw pointer to the shared rows in a point; it goes with the flat working copy
  (milestone 2).

**Defects and debts found, not yet fixed** (owner in brackets):
- The AZUREOut file of an S-factor data set prints S in the cross-section columns; a second S-factor formula lives in
  the point code [data structure, milestone 2: one formula, `physics::sfactor_conversion`].
- A yield on a finite layer is refused: not built [first item of the target work of Phase 5].
- No project exercises a term's fixed-angle frame angle; `tests/sum_lab_frame/models/sum_data.azr` does not fill its
  segment [data structure, milestone 2].
- The capture-A_y terms are still built lazily under a static mutex, and also for a β entrance; two final states
  that share one photon pair would mismatch rows and final states in `GetECIntegrals` and
  `EPoint::CalculateECAmplitudes` (not known whether such a project can load) [pathway table, capture A_y at set-up].
- The decay-position defect of the external-capture code of the data classes [second pathway branch, after the data
  classes are replaced, with a test project].
- `ScatteringLength::ZeroEnergyShift` takes the point-Coulomb S(0) for a pair with a hybrid potential; an exception
  inside the external-capture path under OpenMP ends the process [Phase 5 / Phase 6].
- pyazr and tools debts of the sixth wave (`save_fit` with `norms="fitted"`, a changed potential in an old file,
  matching radius and step, key padding, `classify.py`) [pyazr on records, running].
- Not run on this head: 1 against 8 threads, ThreadSanitizer, `ctest -L slow`, the Sphinx build (not installed here).

**For the PI to see** (behaviour a user can meet; additions):
- Format 3: every number has a unit; a data file may be in keV, degrees or radians, mb, μb, an E_cm energy column; a
  data set can hold S(E) (the engine converts with the masses of the pair) or a thick-target yield per charge or per
  particle.
- `--describe` / `pyazr.describe`: rows read, rows kept and the kept energy range of every data set.
- A potential record can hold `spin_orbit woods_saxon depth radius diffuseness on light|heavy` for a pair with one
  spinning particle; other cases are refused with a code that says why.
- `azure2.ec_integrals` in the formal basis gives other (correct) numbers than before.
- A data line that does not parse is reported before an angle beyond the kinematic maximum of an earlier row.

**Running now** (Opus agents; briefs in `/data0/skowrons/AZURE2-upgrade/work/briefs/`):
- The data structure, milestones 2 and 3 (`X_data_structure.md`, `docs/upgrade/data_structure.md`): the flat working
  copy addressed by the plan, then the deletion of `EData`, `ESegment`, `EPoint`.
- The pathway table, capture A_y at set-up (`cnuc_split.md` §11.6 step 6): the requests from the data-set specs, the
  two readers onto the table, `CNuc::CalcCaptureAnalyzingPower` and its static mutex deleted.
- `pyazr.Project`: the project as records in Python (`Z_pyazr_project.md`; `project_format3.md` §4, §6).

**Then.** The second pathway branch (external-capture amplitudes at set-up, the trees deleted) after the data classes
are replaced; one set-up path for both formats, `predict`, the recipe "a new kind of data"; the rest of `CNuc`;
Phase 5 (targets: the yield on a finite layer first; the lab-frame comparison of shifted data); Phase 6;
the GUI on the new format; `FOR_THE_OWNER.md` and `ROADMAP.md` for the meeting with the owner; `origin/dev` merged only
when the reconstruction is complete.

## 2026-11-10 (agent, brief Y) — the pathway table, read by everything outside the data classes (`feature/upgrade-4-pathway-table-1`)

**Done** (`cnuc_split.md` §11.6 steps 1, 2, 4, 5 and 7 without capture A_y; §11.11 decisions 1, 4, 5, 6, 7).
- `include/PathwayTable.h`, `src/PathwayTable.cpp`: the immutable `PathwayTable` (per entrance pair: exits, spin
  groups, internal and external rows, flat interferences with `z_upos`, interference pathways, diagonal pathways,
  capture-A_y terms), built by `BuildPathwayTable` from `reaction_pathways` and the coupling coefficients.
  `CNuc::Initialize` builds it (`CNuc::BuildPathways`, `src/CNucPathways.cpp`) and holds a
  `shared_ptr<const PathwayTable>` (`CNuc::Pathways()`), shared by every clone.  Capture-A_y terms are built for the
  (entrance, exit) pairs given to `CNuc::RequestCaptureAy` before `Initialize`; the engine requests none yet.
- Readers moved: the pathways and angular-distributions check files, `SetUpFormalCaptureMasks`,
  `AZUREAPI::GetECIntegrals`, `RMatrixFunc::CalculateTMatrix`, `AMatrixFunc::CalculateTMatrix`,
  `GenMatrixFunc::CalculateBareCrossSection` / `CalculateAmplitudeMatrix` / `BuildCoulombInterference`, and the
  gradient (`PointTBar`, `TBarAdjoint`, `CollectBilinearTerms`, `TBarProducts`).  `rt_decay_` and `tbar_decay_` are
  `const ExitRow *`.  A point whose exit pair no pathway reaches stops with a message (it read the last decay).
  The `TablesBuilt() == false` fall-backs of `GenMatrixFunc` are deleted (an abort run of every project first).
- Fix: `GetECIntegrals` (pyazr `ec_integrals`) uses the final state's capture energy `GetECEnergy()`, not the formal
  energy (test `pyazr_ec_integrals_formal`: formal against Brune 2.995 before, 1.3e-13 after).
- Tests: `physics_pathways_from_cnuc` compares table and trees field by field on every `.azr` (Brune, formal);
  `physics_pathway_table` checks properties without the builder (it stays when the trees go).

**Verified** (final code, build `work/b-pathtab-s7`, 0 warnings): 183 x 3 against `tab-trial8` (tag `tab-ptfinal`):
704 / 677 / 704 identical; `ctest -LE slow` 173/173; `tests/run_tests.sh` 142 passed, 0 failed; 1 against 3 threads on
T1-T3,
13N, beta_delayed_alpha, ec_e1m1_phase (6 models), identical_entrance_reaction (3) in three bases: 43 x 3 runs
identical, upos_pair_index and light_excitation (their check.sh at 1 and 3 threads, work trees identical); R-matrix
route (T1, T3, 13N, beam_profile_kernel with `<config>` false and `--no-ignore-externals`, Brune and formal) and
`--use-rmc` (T1, 13N, Brune and formal), base against branch: 12 cases, every output file identical; the four check
files of T1, T3, 13N, beam_profile_kernel: 16 identical; converted-tree gate (tag `gptfinal`): exit 0, 483 (456)
identical, 38 classified, 0 differ, 0 missing; fits at 1 thread (T2 `--use-lm` 200 calls, T2 MIGRAD with the
analytic gradient, 13N MIGRAD): every output file identical.  Before the fall-backs were deleted, a run of every
project with an abort in them (tag `tab-ptabort`) reached neither.  Cost (1 thread, `--no-coulomb-table`, base
`b-trial11` alternated with the branch): T2 calculation cycles -0.6 %, instructions -0.13 %; T2 200-call fit cycles
-1.3 %, instructions -0.49 %; peak RSS unchanged (64 MB, 275 MB).  Largest table (beam_profile_kernel, every
particle and beta pair an entrance): 396,483 rows, 24.1 MB, 0.59 s to build; it is built beside the trees until
they are deleted.

**Open (second branch, after the data classes are replaced)**: step 3 (external-capture amplitudes at set-up in
`EDataExternalCapture.cpp`, `EPointExternalCapture.cpp`, with the decay-position defect, decision 3), step 6
(capture-A_y requests from the data sets, `CalculateCaptureAnalyzingPower` and the capture-A_y block of `PointTBar`
on the table, the static mutex deleted), step 8 (the trees deleted).  `GenMatrixFunc.cpp` keeps 4 tree accesses and
`AMatrixFunc.cpp` 4, all in the capture-A_y code.

## 2026-11-10 (agent, brief W) — the spin-orbit term of the hybrid potential (`feature/upgrade-4-spin-orbit`)

**Done** (`spin_orbit_potential.md` §5 steps 4 and 5, option (a): pairs with one spinning particle).
- The key: channel functions of a pair whose potential has a spin-orbit term depend on (l, J).  `ChannelTable::Key`
  is (pair, l, two_j, E_c); `CoulFunc`, `ShftFunc`, `WhitFunc`, `ChannelFunc`, `NFIntegral` take (l, two_j); the
  radial scratch and one-deep memo of `CoulFunc`, the shift-derivative memo of `CNuc::CalcShiftFunctions`, the THM
  shift memo and the rate grid of `ECIntegral` (now also keyed by the potential) carry it.  One function forms it:
  `PPair::SpinOrbitTwoJ(J)` (`kNoSpinOrbit` = -1 for a pair without the term); `PPair::SpinOrbitFactor(l, two_j)`
  gives <2 l.i>.  The THM Coulomb term passes `kNoSpinOrbit` (the THM vertex refuses any potential).
- The term: `physics/HybridPotential.h` `SpinOrbitTerm` (optional field of `HybridPotential`, compared by value,
  appended to `PotentialKey` only when present), `TwoLDotI`, `SpinOrbitRadialTerm` (FRESCO: V_so 2.000 fm^2 (1/r)
  df/dr); `RadialPotential`, `IntegrateInward`, `DecayingSolution`, `HybridBoundShift` take <2 l.i>, added last and
  not at all when 0.  Bound states, closed channels, S', S'' and the EC final state use it.
- Format 3: `spin_orbit woods_saxon depth radius diffuseness on light|heavy operator 2l.s` (schema row; --describe;
  the record writer writes it for pyazr's save_fit).  Codes (`experiment/SpinOrbitChecks.h`, shared by the build,
  pyazr and the API): `E_SO_TWO_SPINS`, `E_SO_ZERO_SPIN`, `E_SO_RELATIVISTIC`, `E_SO_THM`, `E_SO_NEEDS_J`,
  `E_SO_NO_TERM`.  `PPair::SetHybridPotential` refuses an impossible term as a last guard.
- pyazr 2.24.0: `set_nuclear_potential(..., spin_orbit={depth, radius, diffuseness, on})`,
`NuclearPotential.spin_orbit`,
  `coulomb_functions(..., J=None)`.  `AZUREAPI::GetCoulombFunctions` moved verbatim to
  `api/src/AZUREAPIChannelFunctions.cpp` first (size limit).
- Tests: `physics_spin_orbit` (m-scheme angular factor, radial term, delta-shell closed form, the DOP853 table
  `physics/tests/data/spin_orbit.csv`, bit-identity at <2 l.i> = 0), `tests/spin_orbit_potential` (χ² pinned; check.sh:
  phase shifts against the table, unitarity, 1 vs 3 threads, gradient, capture), `pyazr_spin_orbit_potential`
  (coulomb_functions with J, API = record, refusals, EC integrals per J against `ec_reference.py`),
  `channel_table_unit` with both j.  Docs: `theory/hybrid_potential.rst`, `examples/spin_orbit_potential.rst`,
  `pyazr.rst`, format3.rst regenerated.

**Verified** (build `work/b-so`, clean, 0 warnings): 183 × 3 against `tab-trial8` (tag `tab-so3`): 704 / 677 / 704
identical, plus 15 / 15 / 15 new files (the new project and its models); channel tables of T1–T3 byte-identical
to the base build (a scratch dumper linked against each); T2 cycles (perf stat, 1 thread, base and branch alternated,
6 repetitions, minimum): calculation 3.326e9 against 3.326e9, 200-call fit 2.122e11 against 2.124e11 (+0.1 %);
instructions +0.38 % and +0.08 %; the remaining checks are in the report of the branch.  Every new test was shown to
fail once with an injected error (l.s for 2 l.s, the sign of df/dr, <2 l.i> halved, the point's row read with
J = l + 1/2, the EC bound state without its J).

**Open.** Web UI / Qt forms for the word (they need: the spin_orbit structured value of the schema row, pyazr's
dict, the J argument of coulomb_functions).  `gaussian` as a spin-orbit form (one schema alternative, one solver
case, one test).  Option (b), pairs with two spins (refused).  `ScatteringLength::ZeroEnergyShift` (the exact
threshold of a scattering-length point) uses the point-Coulomb S(0) for a hybrid pair, not the distorted one (older
than this branch).  Sphinx not run (not installed on this host).

## 2026-11-10 (agent, brief X, milestone 1) — what was measured exists once (`feature/upgrade-4-data-structure`)

**Design.** `docs/upgrade/data_structure.md`: `Measurement` (immutable, shared), the data set, a flat working copy;
where every member of `EPoint`/`ESegment`/`EData` goes; the three milestones.

**Done (milestone 1).**
- `experiment::Measurement` (`experiment/include/azure2/experiment/Measurement.h`): the rows of a data set inside its
  window as read (E_lab, angle, value, error after a covariance file, bins, file record, the Legendre-order flag) and
  their provenance (file, rows read, rows kept, kept E_lab range), held as `shared_ptr<const>` by the segment and its
  copies; an extrapolation's grid is a measurement too.  `ReadMeasurementRows` is the one reader (`DataLine` moved
  to `experiment/`); the engine and `--describe` use it.
- The points read their datum by row: `lab_crosssection_`, `lab_dcrosssection_`, `bin_low_lab_`, `bin_high_lab_`,
  `fileRecord_`, `coefficient_order_`, `coefficient_absolute_`, `original_energy_` (of a row) are gone; a sub-point
  takes its row's Legendre order (`EPoint::BelongToRowOf`).
- One lab -> c.m. conversion, `experiment::ConvertRow`, writing a `MeasurementFrame` row; the load, the shift step and
  the leaves call it (`ESegment::ConvertRowToCM`); a shifted row starts from the frame at load (H4).  y_c.m., dy_c.m.
  and the Jacobian live only in the frame (`ESegment::CMValue`, `CMError`, `Jacobian`); `PointResidual` takes the
  segment and the row.
- Components are leaves of their segment (`ESegment::AddLeaf`, term k >= 1; the segment's own points are term 0):
  `EData::componentSegments_`, `CreateComponentSegment` and the pointer re-mapping of `Clone` are gone; the shift
  step, grid re-placement and target set-up take (segment, term) (`TermRef`); `EData::LeafRefs` keeps creation order
  (plan, layout, intEC file).
- `--describe [--json]` / `pyazr.describe` (pyazr 2.25.0): rows read and kept and the kept E_lab range per data set.

**Verified** (build `work/b-dstruct`, 0 warnings): 183 x 3 against `tab-trial8` (tag `tab-dsm1a`): 704 / 677 / 704
identical; `ctest -LE slow` 172/172 (upgrade_history, upgrade_plan_runtime, upgrade_poison_results among them);
`tests/run_tests.sh` 141 passed, 1 failed (the story `--describe` of `tests/ratio_to_rutherford`: the new rows line;
expected file updated, its check then passed); converted-tree gate (`gdsm1a`): 483 (456) identical, 38 classified,
0 differ, 0 missing; `tests/format3_engine/check.sh` passed; calculations at 1 and 3 threads identical on T1-T3, 13N,
natural_neon (three bases); fits at 1 thread identical to the base build (T1, T2, 13N), at 3 threads they differ from
1 thread on the base build as well (T2: edm, chi2 in the last digits; T1, T3, 13N: the call count) -- not from this
work.  Speed (min of 6 alternated runs, 1 thread, Coulomb table): T2 calculation 668.8 M cycles against 672.3 M,
200-call T2 fit 53.47 G against 53.39 G; peak RSS T2 fit 349 / 427 MB at 1 / 3 threads against 360 / 441,
natural_neon 1532 MB against 1551.

**Open (milestones 2 and 3).** The working copy is not yet flat: units are still `EPoint`s inside segments (term 0 and
leaves), each with its settings; the plan addresses (list, segment, point); `Clone` copies segments.  `EPoint`,
`ESegment`, `EData` still exist.

## 2026-11-10 (agent, brief J) — units in format 3, the S factor as a quantity, yields per charge (`feature/upgrade-4-units`)

**Done.**
- Units: every number of a format-3 project may be written in any unit of its dimension.  One function applies a
  factor (`experiment/Units.h` `Rescale`, through `ToBase`, `ToUnit`, `Context::InBase`): a power-of-ten factor
  divides by the exact integer, so a value exact in the written unit is the same double as the value written in the
  engine unit.  Record fields (windows, shifts incl. the sqrt term, grids, areal density in atoms/b for a layer and a
  sample, U_e in keV, a straggling coefficient in any sqrt-energy unit, stopping and spread formulas in another unit,
  UDR times, model data) are converted in `Build`.  The columns of a data file: `DatasetSpec::columns`
  (`ColumnConversion`: energy, angle, value factors, an `E_cm` column with the entrance pair's masses), applied by
  `experiment/DataRows.h` `RowToEngine` through `DataLine::ToEngineUnits`, called in `ESegment::Fill`; the energy
  window is in the energy column's frame and unit.  An extrapolation in another unit is computed in it.
- S factor: `quantity sfactor` builds; `ReactionSpec::scale` (`ValueScale`) with `sFactor`; tree kind `Scaled`
  (one row of `TreeKinds()`), `FactorRef::Kind::SFactor`, its value `SFactorConversionValue` (src/CoulombReference.cpp)
  from `physics/AstrophysicalFactor.h` `sfactor_conversion` (E exp(2 pi eta), non-relativistic eta, <levels> masses);
  binding `_azure2.sfactor_conversion`.  Adaptive extrapolation grids of `sfactor` refine on S (the model value).
- Yields: `1/uC|mC|C` with `charge_state` (x 1e-6 C/(q e) ...), `counts/...` with a fixed `efficiency`, `arbitrary`;
  a layer target is refused for a yield (Phase 5 kernel), a fitted efficiency is refused.
- pyazr 2.27.0: `sfactor_to_cross_section(S, E_cm_kev, z1, z2, m1_u, m2_u)` via the engine; `to_azr` refuses S-factor
  data; `ExforData.to_sfactor_dataset` (pyazr/exfor_sfactor.py) writes `quantity sfactor`.
- Defects fixed on the way: `unit 1/C` / `counts/C` (factor 1) were built with no scaling; a yield on a layer built as
  the target-averaged cross section; a sqrt shift, straggling coefficient, UDR times read without their unit; percent
  priors of U_e and sample thickness read as absolute.
- Tests: `experiment_build` (units, sfactor, yields), `experiment_prediction` (Scaled), `pyazr_units_scales`,
  `nds_local_exfor` section 6, `tests/units_scales/check.sh` (the section 9.1 example on T2's levels, adaptive S grid).

**Verified** (build `work/b-units`, 0 warnings): 183 + 1 projects x 3 bases against `tab-trial8`: 704 / 677 / 704
identical, the 7 other files per basis those of the new `tests/units_scales/models/yield_story.azr` (tag
`tab-units-final`); converted-tree gate (`gunits2`): 483 (456) identical, 38 classified, 0 differ, 0 missing; spec
equivalence 194 files, 185 converted, 183 equal, 2 recorded, 0 differ; `ctest -LE slow` 172/172 (after adding the new
project's three lines to `tests/upgrade/predictions.txt`); calculations at 1 and 2 threads on T1-T3, 13N in three
bases identical; every new test shown failing once with an injected error (multiplying instead of dividing in
`Rescale`, the S factor at the lab energy, charge state 1, the S-factor node dropped, pi = 3.1416, mass numbers in
`sfactor_to_cross_section`).  `tests/run_tests.sh` (OMP 2) 143 passed, 0 failed (with `units_scales`).

**Open.** AZUREOut files of an S-factor or per-charge data set print the model in the cross-section column and that
times E exp(2 pi eta) in the S column (the output writer is in the data-class files); `calculate_sfactor_rwa` likewise.
`EPoint`'s own S-factor conversion (output column) is a second formula of the same factor (bit-different association).

## 2026-11-10 (director, night) — the sixth wave merged; a faster regime; the data structure is replaced, not carved

**Merged into `feature/upgrade-4-data`** (`--no-ff`; each rehearsed, built and verified in a scratch worktree, the
committed tree is that tree): the tools follow-up (the converted-tree gate as `tools/format3_gate.sh`), the
pathway-table design (`cnuc_split.md` §11, with the director's decisions §11.11), the potential as a record of the
project (`potential` in `<potentials>`; pyazr 2.23.0), and the cut of `EData`, `EPoint`, `ESegment` into files with
one job (`edata_split.md` §5.8). Their entries are below.

**The merged head as a whole** (build `work/b-trial11`, 0 warnings): 183 projects × 3 bases identical to the
reference tag `tab-trial8` (704 + 677 + 704 files); `ctest -LE slow` 171/171; `ctest -L slow` 10/10 (run once, on the
head before the potential record); `tests/run_tests.sh` 142 passed / 0 failed; calculations at 1 and 6 threads 12/12
identical; T2 with a `<modelData>` section identical in the three bases; spec equivalence: 193 files, 185 converted,
183 equal, 2 recorded, 0 differ; the converted tree: 483 (456 formal) identical, 38 classified, 0 differ, 0 missing.
All of it was run on the merged tree of the last merge (`work/director/verify_trial11.out`). Not run on this head:
1 against 8 threads, ThreadSanitizer.

**The PI, 2026-11-10 (evening): go faster; the physics must stay the same.** What changes:
- An agent checks each commit quickly (build, consistency test, the unit tests of the change, a handful of
  projects) and runs the whole comparison once per branch or milestone; speed is measured once, at the end. The
  director verifies once per merge and merges several branches together (`work/director/full_verify.sh`).
- The data classes are not carved further. `EData`, `ESegment`, `EPoint` are **replaced**: what was measured is one
  immutable shared object per data set, the working copy is a flat vector of calculation units addressed by the
  plan, the functions that compute move onto these unchanged, and the three classes are deleted
  (`work/briefs/X_data_structure.md`; its design page will be `docs/upgrade/data_structure.md`). This takes the place
  of the steps S7 and S12 of `data_model.md` §6 and of the data-path part of Phase 6; `edata_split.md` §5 from F1 on
  is superseded.

**Decisions of the director.**
- The pathway table (`cnuc_split.md` §11) is built in two branches: first everything outside the data classes (the
  table, the check files, the set-up readers, the R-matrix route, the forward hot path), then, after the data
  structure is replaced, the external-capture set-up, capture A_y at set-up and the deletion of the trees.
- The spin–orbit term follows the PI's four answers (one spinning particle; FRESCO's convention; the same potential
  for bound states and closed channels; e² = α ħc): `work/briefs/W_spin_orbit.md`.

**Defects and debts found, not yet fixed** (owner in brackets):
- An exit pair number is used as a decay position in the external-capture code of the data classes and in
  `AZUREAPI::GetECIntegrals`; no project of `tests/` triggers it [second pathway branch, with a test project].
- The not-found check of the exit in `GenMatrixFunc`/`AMatrixFunc` cannot be reached; the capture-A_y terms are built
  lazily under a static mutex; `GetECIntegrals` uses the formal energy in every basis [pathway branches].
- pyazr: `save_fit` of a format-3 project fails with `norms="fitted"`; `save_fit` of an old-format project does not
  write a changed potential; `set_nuclear_potential` takes no matching radius or step [pyazr on records].
- `io::WriteRecords` pads keys to 15 characters (`matching_radius` is misaligned); `classify.py` treats two
  unreadable files as equal [tools].

**Then.** Units, S factor and yields per charge (`J_units_scales.md`); pyazr's model on records, one set-up path for
both formats, components as leaves in `predict`, the recipe "a new kind of data"; the rest of `CNuc`; Phase 5 (with
the lab-frame comparison of shifted data); Phase 6 (what is left of the mutable statics); `FOR_THE_OWNER.md` and
`ROADMAP.md` for the meeting with the owner; `origin/dev` merged only when the reconstruction is complete.

**Running now** (Opus agents; briefs in `/data0/skowrons/AZURE2-upgrade/work/briefs/`, the PI allowed all 24 threads):
the replacement of the data structure (`X_data_structure.md`, `feature/upgrade-4-data-structure`), the pathway table
outside the data classes (`Y_pathway_table_1.md`, `feature/upgrade-4-pathway-table-1`), the spin–orbit term
(`W_spin_orbit.md`, `feature/upgrade-4-spin-orbit`), units, S factor and yields per charge (`J_units_scales.md`,
`feature/upgrade-4-units`).

## 2026-11-10 (agent, brief S) — `EData`, `EPoint`, `ESegment` cut into files with one job (`feature/upgrade-4-edata-split-3`)

**What file holds what now.** Eleven new source files and two new headers, each made by a verbatim move (bodies
unchanged, one commit per file); the file map of all three classes is `edata_split.md` §5.8.
- `EData`: `src/EDataTargetEffects.cpp` (target effects, sub-point grids, windows, grid re-placement; the
  component-segment grid block stays inside `SetUpTargetEffects`, which S12 rewrites), `src/EDataExternalCapture.cpp`
  (EC amplitudes of the data, the `intEC` file, its count and signature; `ECProgressRedraw` became a private static
  member of `EData` because `InitializeComponentSegments` in `EData.cpp` uses it too), `src/EDataParameters.cpp`
  (parameter-block offsets, names and start values, shift/norm/thickness groups, density parameters, prior centres),
  `src/EDataHandOver.cpp` (a parameter vector into segments, effects and points; the four target-effect prior totals,
  with four check-11 allowance rows).
- `EPoint`: `src/EPointCalculate.cpp` (`Calculate`, thin target, observable-wise integration, `MarkFailed`,
  `ScopedPointMappingOff` with its thread-local flag), `src/EPointFolds.cpp` with `include/TELinearForm.h` (the
  apparatus folds), `src/EPointChannelTables.cpp` with `include/EPointTables.h` (channel tables, sharing, Coulomb
  amplitude, getters; `EDependentTables`, `EDependentShare`, `ThmPsTable` are namespace-level types that `EPoint`
  names by aliases; check 25 reads the new header), `src/EPointLegendre.cpp` (Legendre polynomials, Q-coefficient
  average), `src/EPointExternalCapture.cpp` (the point's EC amplitudes, also at a shifted energy).
- `ESegment`: `src/ESegmentParameters.cpp` (shift values, shift and norm groups, shift function, shift prior),
  `src/ESegmentShifts.cpp` (re-energising after a shift, channel tables of the last shift states).
Sizes: `EData.cpp` 4,454 → 1,950, `EPoint.cpp` 2,920 → 951, `ESegment.cpp` 1,366 → 881, `EPoint.h` 896 → 766,
`EData.h` 560 → 565; the largest new file is `EDataTargetEffects.cpp` (964).

**Not moved, and who takes it.** What is left in `EData.cpp` is what later steps delete or rewrite: the component
segments (S12), point mapping, result layout, `Clone` and the parallel pass (the rest of S5), `Fill`/`MakePoints` and
the load-time refusals (after 3d).  `EPoint.cpp` keeps the datum and frame code (S7), `ESegment.cpp` the frame code
(S7) and the composition helpers (S10–S12).

**One code change besides the moves.** Once `EPoint.cpp` had shrunk, GCC 11 (-O3) inlined the file-local
`LegendreSeries` into `EPoint::CalcLegendreP`, and the T2 calculation ran 3.6 % more instructions and 3.8 % more
cycles than the base (callgrind: all of it there).  A separate commit keeps the function out of line
(`__attribute__((noinline))`; the CI builds GCC, Apple clang and MinGW GCC, not MSVC).  Cycles, minimum of 6 runs
alternated with the base build at 1 thread: calculation −2.6 %, without the Coulomb table −0.6 %, 200-call fit −2.9 %;
179 × 3 outputs identical before and after.

**Verified.** Every move: 0 warnings; the diff of the removed lines against the new file's body and of the rest against
the old file less those lines, both empty (shown exceptions: the namespace and signature lines of `ECProgressRedraw`,
the alias lines of `EPoint.h`, the indentation of the moved table types); consistency test after `--tighten-sizes`;
179 × 3 identical to `tab-p3f` (689 / 662 / 689; tags `tab-esplit3-c1` … `-c12`); `scripts/check.sh --quick`
(ctest 165/165).  Also `tests/progress_bars` after the EC move and `tests/shift_fn`, `tests/energy_shift_sqrt` after
the shift moves.  End of work (`TMAX=3 NT=3 work/verify-p4c.sh wt-esplit3 esplit3`): 179 × 3 identical, `ctest -LE
slow` 165/165,
`tests/run_tests.sh` 141 passed / 0 failed, 1 against 3 threads on T1–T3 and 13N in three bases 12 identical;
`tests/format3_engine/check.sh` passed.

**Open.** The branch starts from a320e818; the director rebases it onto the moved integration branch (the hybrid
potential and the extrapolation grids edit `src/CMakeLists.txt`, `file_sizes.tsv` and possibly `EData::MakePoints`).

## 2026-11-10 (agent, brief U) — the hybrid potential as a format-3 record (`feature/upgrade-4-potential-record`)

**Done.**
- A format-3 project states a pair's hybrid potential as a `potential` record of a new section `<potentials>`:
  `pair` (the pair label; `pair_key` only where two pairs have the label, the data sets' rule), `coulomb_radius`
  (required), `central woods_saxon depth radius diffuseness | gaussian depth width`, `matching_radius` and `step`
  (defaults 15 fm, 0.0275 fm, noted in the run log when used). No spin–orbit term. A pair without a record has no
  potential; there is no default record. Every word is a row of the schema table (`--schema`, `pyazr.schema()`,
  `format3.rst`); `--describe` prints each pair's potential with units.
- Matching radius and step were already fields of the potential value and in every table key (`PotentialKey`), so
  they are fields of the record. The Coulomb radius was already a field of the value (3.5 fm by default); the old
  block keeps giving 3.5 fm.
- Refusals with codes of the table: `E_POTENTIAL_TWICE`, `E_POTENTIAL_PAIR_KIND` (photon or beta pair),
  `E_POTENTIAL_RELATIVISTIC`, `E_THM_HYBRID_POTENTIAL` (a THM data set or extrapolation whose entrance pair has a
  record; the engine's own check stays for old projects), `E_OLD_POTENTIAL_BLOCK` (`<potential>` in a format-3
  file, naming the converter); the existing `E_UNKNOWN_PAIR`, `E_AMBIGUOUS_PAIR`, `E_PAIR_KEY`, `E_FIELD_REQUIRED`,
  `E_OUT_OF_RANGE`, `E_UNIT_NOT_ALLOWED`, `E_PART_MISSING` for the record's fields.
- Files: `experiment/` (schema row, `BuildPotentials.cpp`, `HybridPotentialSettings` moved here from `include/` and
  `src/`, `build::ResolvePair` shared with the data sets, `Describe.cpp`, `ProjectReader.cpp`); the old block's parser
  moved verbatim to `convert/include/azure2/convert/PotentialBlock.h` (engine and converter read it with one parser);
  `convert/PotentialRecords.h` writes records from settings (the converter and pyazr's writer).
- Converter: one record per particle pair the old block gives a potential (default expanded, `coulomb_radius 3.5 fm`,
  r_m and step written out), notes for sections that write nothing; the format-3 file has no `<potential>`.
- pyazr 2.22.0: `set_nuclear_potential(..., coulomb_radius=)`, `NuclearPotential.coulomb_radius`,
  `azure2.potential_records()`, `save_fit` of a format-3 project writes the session's potentials as records when they
  changed. The potential methods live in `pyazr/potentials.py` (moved verbatim first).

**Verified** (build `work/b-potrec`, 0 warnings, 3 threads): 183 × 3 against `tab-trial8` (tag `tab-pr1`): 704 / 677 /
704 identical; `ctest -LE slow` 170/170; `tests/run_tests.sh` 142 passed / 0 failed; spec equivalence 193 files, 185
converted, 183 equal, 2 recorded, 0 differ (now comparing every particle pair's potential by bits; injection 4 fails
hybrid_potential, T3, T3_start); converted-tree gate (`REF=tab-trial8`, tag `prg1`): 472 (445) identical, 38
classified, 0 differ, 11 missing, as before; converted `hybrid_potential`, T3, T3_start at 1 and 3 threads in the
three bases: 87 / 87 files identical; `pyazr_potential_record`: R_C = 5 fm against DOP853 P 1.7e-8, S 7.2e-8. Every new
check shown failing once with its check removed or its value moved (the report lists the injections).

**Open.** pyazr's `save_fit` of a format-3 project works only with `norms="nominal"` (`"fitted"` reads
`<segmentsData>`); format-2 `save_fit` still does not write a potential changed in the session. Web UI / Qt forms
for the record not made. The spin–orbit term and `two_j` (note §5 steps 4–5).

## 2026-11-10 (agent, brief V) — design of step 6b of the `CNuc` split: the readers of the pathway trees (`docs/upgrade-pathways-6b`)

**Done** (documents only; no code changed).
- `cnuc_split.md` §11, "Step 6b: design": the inventory of every function that reads the pathway trees (320 accesses
  outside `CNuc.cpp`, 86 in it, per function with the fields used, hot or cold, OpenMP or not), the one writer after
  set-up (`CNuc::CalcCaptureAnalyzingPower`, lazily in the point pass under a function-level static mutex), what the
  nodes hold beyond `reaction_pathways`' result (only structural data: three derived lists per decay), the table
  (`PathwayTable`, flat rows per entrance in tree order, engine numbers beside System indices, built once by
  `CNuc::Initialize`, held as `shared_ptr<const>` and shared by every working copy), the hot loop before and after,
  the UPOS, beta-delayed, channel-capture and capture-A_y shapes, eight bit-identical steps with their counts and
  tests, the measurement and stop rule, and six defects met while reading. "Remaining steps" points to §11.

**Verified how.** Counts by `work/pathnote/inventory.py` and `grep -o` at cb8f5912; the decay numbering of 14 of the 21
external-capture projects with the pathways check file (`work/pathnote/ecprobe.sh`, base build `work/b-trial9`, one
thread): every decay position equals its exit pair number; `perf stat` (five runs) and `perf record` of the T2
calculation and of a 30-call T2 fit on the base build, one thread. `python3 tests/consistency/repo_consistency_test.py`
before the commit.

**Open.** The director's decisions listed at the end of the report (where the table's handle lives after step 7, the
capture-A_y request list against building for every photon exit, the defect of an exit pair number used as a decay
position, the fall-back branches). Not run: the cost of building capture-A_y terms for every photon exit; whether a
sub-point of a capture-A_y point can be the first caller with no Legendre order.

## 2026-11-10 (agent, brief T) — the comparison tools in the repository, and small follow-ups (`feature/upgrade-4-tools`)

**Done** (no physics change; every project output identical to the start).
- The comparison of two engines is in the repository: `tools/runall.sh`, `tools/run_one.sh`, `tools/cls.sh`,
  `tools/classify.py`, copied from the scripts outside the repository with their behaviour kept; `JOBS` defaults to 1,
  usage errors exit 2, each has a header (what, environment, output).  `format3_gate.sh` takes `runall.sh` from beside
  itself.  Described in `CLAUDE.md` (recipes) and `ARCHITECTURE.md` §10 "Comparing the outputs of two engines".
- No path of this machine in a tracked file (outside this log): the docs name `tools/...`; consistency check 30 fails
  on a path under /data0 or /home (three allowed lines, of the skills checker and its test).
- The format-3 gate keeps the projects that are format 3 already (the four of `tests/ratio_to_rutherford/models/`) and
  runs them unchanged; a project left out is listed in `tree.txt` with the converter's reason.
- `scripts/check.sh`: one `JOBS` (default 4) for the build and ctest.
- `tests/hybrid_potential/check_per_pair.py`: its reference comes from `expected/chiSquared.out` and `tolerance` (the
  old pin 2462154.2 was stale; the engine gives 2466142.3); it runs in a temporary copy; ctest `pyazr_hybrid_per_pair`
  (16 s).
- The data check file (`checks/data.chk`) of an extrapolation run with an adaptive grid is written after the
  refinement and lists the refined grid (378 points, not the 21 of the start grid; `src/AZUREMain.cpp` only).
- `--print-predictions` and `pyazr.predictions` label a format-3 data set and extrapolation by the record's name
  (`src/Predictions.cpp`); old projects keep their numbers.  Re-pinned: `tests/upgrade/predictions.txt` (8 lines of the
  four format-3 projects), `tests/ratio_to_rutherford/expected/story_predictions.txt`.

**Verified** (build `work/b-tools`, 0 warnings, 3 threads): 183 × 3 against `tab-trial8`: 704 / 677 / 704 identical,
once with the committed tools (`tab-tl1`) and once with the outside ones (`tab-tl2`), the two classification tables
byte-identical; the gate (`REF=tab-trial8`, tag `tlg1`): 521 (494) files per basis, 483 (456) identical, 38
classified, 0 differ, 0 missing; `JOBS=1 scripts/check.sh --quick`: ctest -LE slow 170/170; consistency on the tree of
every commit; `tests/extrapolation_grids/check.sh` (12 ok) and `tests/ratio_to_rutherford/check.sh` (7 ok).  Not run:
the whole `tests/run_tests.sh`, Sphinx, 1 against 8 threads.

**Open.** pyazr's version is not bumped for the new labels of `pyazr.predictions` (format-3 projects only).

## 2026-11-10 (director, evening) — Phase 4: the fifth wave merged (dead code out of the data classes, ratio to Rutherford, the hybrid potential as a value, extrapolation grids)

**Merged into `feature/upgrade-4-data`** (`--no-ff`, in this order; each rehearsed, built and verified in a scratch
worktree, the committed tree is that tree; the four agents' entries are below):
1. ed860f74, `feature/upgrade-4-edata-split-2`: unused code out of `EData`, `EPoint`, `ESegment`; the data covariance
   in `src/EDataCovariance.cpp`.
2. 6392bb62, `feature/upgrade-4-rutherford`: data given as a ratio to the Rutherford or the Mott cross section
   (`quantity ratio_to_rutherford`, `ratio_to_mott`).
3. 8a143945, `feature/upgrade-4-hybrid-potential`: the hybrid potential is a value in `physics/` (the manager
   singleton is gone), e² = α ħc in its Coulomb part, bound states and closed channels of a hybrid pair use the
   solution distorted by the potential.
4. 08a6165e, `feature/upgrade-4-extrap-grids`: extrapolation grids by step, by number of points (linear or
   logarithmic), by list, from a file, in the laboratory or the c.m. frame, and adaptive; the converted-tree gate as
   `tools/format3_gate.sh`.

**The merged head as a whole** (08a6165e; build `work/b-trial8`, 0 warnings):
- 183 projects (179 and the four of `tests/ratio_to_rutherford/models`) × 3 bases: 704 + 677 + 704 output files.
  Of the 689 + 662 + 689 that the Phase 3 head has, 661 + 634 + 661 are identical; the 28 others per basis are
  `hybrid_potential`, `upgrade/T3/T3` and `T3_start`, moved on purpose by merge 3 (cause and size in its entry and
  its merge message). **The reference for bit-identity from here on is the tag `tab-trial8`** (every file of the head).
- `ctest -LE slow` 169/169; `tests/run_tests.sh` 142 passed / 0 failed.
- Calculations at 1 against 3 and at 1 against 6 threads, T1–T3 and 13N in the three bases: 12/12 identical.
- T2 with a `<modelData>` section: identical at 1 and 3 threads in the three bases.
- Spec equivalence: 193 files, 185 converted, 183 equal, 2 recorded, 0 differ, 8 not converted (format 3 already,
  or test inputs that are not projects).
- Converted tree against the same engine's run of the unconverted tree (`tools/format3_gate.sh`, `REF=tab-trial8`):
  472 (445 formal) identical, 38 classified, 0 differ, 11 missing — the two projects written in format 3 with data,
  which the gate does not carry yet (in work: `feature/upgrade-4-tools`).
- T2 benchmark against the Phase 3 head, interleaved, three repetitions, the machine shared with two agents:
  no change within the scatter. Calculation 0.97–1.13 s (Phase 3 head) against 0.98–1.14 s at 1 thread, 0.26–0.32 s
  against 0.27–0.34 s at 6; 200 fit calls 61.7–70.6 s against 63.3–70.4 s at 1 thread, 30.1–33.6 s against
  30.4–33.6 s at 6; peak memory of the fit 276 → 269 MB at 1 thread, 506–513 → 497–509 MB at 6; the same χ²
  (10899.8; 691.608 after 200 calls) in every run (`work/bench-p4w5.tsv`).

**What the director changed at the merges** (each in its merge commit's message):
- `src/AdaptiveExtrapolation.cpp` (merge 4) read a point's model value through `EPoint::ModelSlot()`, which merge 1
  had removed as unused; it uses `EPoint::GetModelValue(results)`. A branch that adds a caller of a name another
  branch deletes builds alone and fails merged: the rehearsal build is what caught it.
- pyazr is **2.21.0** (2.19.0 with the ratio data, 2.20.0 with the potential, 2.21.0 with the grids).
- `tests/consistency/file_sizes.tsv` tightened to the merged counts at each merge.

**Decisions of the director.**
- Ratio to Rutherford: `ExforData.to_azr(rutherford=True)` no longer multiplies the data by σ_R (the engine
  compares the ratio; `to_dataset` writes `quantity ratio_to_rutherford`). The four refusals stay (a ratio with
  Doppler broadening, with electron screening, with attenuation, with `auto_tolerance`): each needs a decision on
  what the reference cross section is under that effect. Nothing more for identical particles (Mott is its own
  quantity).
- Hybrid potential: `GetCoulombFunctions` / `pyazr.coulomb_functions` return the shift function at E ≤ 0 for every
  pair (it returned 0); a THM segment whose entrance pair has a hybrid potential is refused
  (`E_THM_HYBRID_POTENTIAL`); the T1 goldens stay as recorded.
- Extrapolation grids: for `quantity sfactor` the refinement will test the S factor (with the units and S-factor
  step); adaptive grids stay refused in a pyazr session (`pyazr.refine_grid` refines on a Python callable).
- The comparison tools (`runall.sh`, `run_one.sh`, `cls.sh`, `classify.py`) go into the repository's `tools/`; no
  path of one machine in a tracked file (a consistency check).

**Defects and debts found, not yet fixed** (owner in brackets):
- `tests/hybrid_potential/check_per_pair.py` (not in ctest) has a stale χ² pin; the data check file of an adaptive
  extrapolation lists the start grid; `scripts/check.sh` builds with as many jobs as there are cores;
  `--print-predictions` names a format-3 data set by its index [all: `feature/upgrade-4-tools`, running].
- Adaptive refinement finds what the start grid and the level seeds bracket: a narrow structure that is not at a
  level energy (an interference minimum between coarse start points whose ends and midpoint happen to lie on a
  line) is not seen. The director added that limit to the description of adaptive grids (`data_formats.rst`): the
  start step must resolve the structures that are not at a level. A second criterion (curvature from four points)
  is a possible later addition [not scheduled].
- Still not run by anyone: ThreadSanitizer (no runtime on this host), `ctest -L slow`, 1 against 8 threads on this
  head (1 against 6 was run).

**For the PI to see** (behaviour a user can meet; additions):
- Ratio data under a beam spread or a target: the model is ⟨σ⟩/⟨σ_R⟩ (each folded, then divided), which is what an
  experiment that normalises counts to a Rutherford calculation for the mean energy does not exactly give; the
  documentation states the convention.
- pyazr: `to_azr(rutherford=True)` is gone (use `to_dataset`); `coulomb_functions` below threshold returns S.
- Hybrid potential (`hybrid_potential` project): the internal reduced width amplitude of the ¹³N ground state is
  1.060263 MeV^½ (was 1.271943) for the same ANC, because the bound state is now distorted by the potential (the
  PI's answer to question 4); χ² of its segment 3 2238.82 → 2259.87. The e² change moves hybrid outputs by ≤ 5e-5.
- Extrapolations: `range A B points N log`, `list`, `file`, `frame cm`, `adaptive on tolerance T max_points N`; a
  budget that stops the refinement is reported with the largest deviation left.

**Running now** (Opus agents, 3 threads each; briefs in `/data0/skowrons/AZURE2-upgrade/work/briefs/`):
- `feature/upgrade-4-edata-split-3` (`wt-esplit3`, `S_edata_split_3.md`): the remaining cuts of `EData`, `EPoint`,
  `ESegment` into the files of `edata_split.md`, as verbatim moves (started from the tree of merge 2: rebase onto
  6392bb62, then merge over 3 and 4).
- `feature/upgrade-4-tools` (`wt-tools`, `T_tools_followup.md`): the comparison tools into `tools/`, the gate
  carries format-3 projects, the small defects above (started from the rehearsal commit af654a5d of merge 4: rebase
  `--onto 08a6165e af654a5d`).

**Then.** Units and the quantities that are scales (S factor with the `nds.py` mass fix, yields per charge;
`J_units_scales.md`), after the cuts are merged (it edits where the points are made); `Measurement` shared between
working copies (`R_measurement.md`), after the cuts; pyazr's model on records, one set-up path for both formats,
components as leaves, `predict`, the recipe "a new kind of data"; the spin–orbit term (note §5 steps 3–5); the
rest of the `EData` and `CNuc` splits; Phase 5 (targets and kernels, shifted lab data compared in the lab frame);
Phase 6 (mutable statics, the task as unit of work). Merge `origin/dev` only when the reconstruction is complete.

## 2026-11-10 (agent, brief K) — extrapolation grids: list, file, points/log, c.m., adaptive; the gate as a tool (`feature/upgrade-4-extrap-grids`)

**Done.**
- bad7a299: the converted-tree gate is a committed tool, `tools/format3_gate.sh ENGINE SOURCE_TREE TAG` (`WORK`, `REF`,
  `JOBS`, `RUNALL` from the environment) with `tools/format3_gate.py` (`tree`, `compare`; the accepted differences are
  the explicit table `CLASSES`: class, file, projects, why); described in `project_format3.md` §8.5.  runall.sh is not
  in the repository: the gate takes it from `RUNALL` (default: beside the script).
- 05e8f975: one reader of the grid forms, `experiment/Grids.h` (`ReadGridValues`): `range … step` (F10), `range …
  points n [log]`, `list`, `file "path"` (relative to the project, `#` comments), the start grid of `adaptive`; units
  through the unit table; a `cm` energy grid taken to lab by `Frames`' new `BeamLabFromCM` with the `<levels>` masses
  (`PairInfo::lightMass/heavyMass`).  `--describe` prints a grid file's path and count.
- 011880d0: `experiment/AdaptiveGrid.h` (`RefineGrid`, a pure function of a callable model); the extrapolation run
  refines adaptive grids before the data set is initialized (`src/AdaptiveExtrapolation.cpp`, three lines in
  `AZUREMain`; no change to `EData`/`EPoint`/`ESegment`): level seeds, batches through the point pass on a sweep data
  set made by `EData::MakePoints`, union grid over angles, `max_points` warning, then the run's data set on the refined
  grids.  pyazr 2.18.0 unchanged: `pyazr.extrapolation_grid`, `pyazr.refine_grid` (`pyazr/grids.py`,
  `pyazr/__init__.py`); a pyazr session refuses an adaptive grid.  Docs: `data_formats.rst` "Extrapolation Grids",
  `calculate.rst`, `pyazr.rst`, `project_format3.md` §3.8.

**Verified** (build `work/b-grids`, 0 warnings, 3 threads): 179 × 3 identical to `tab-p3f` (689/662/689, `tab-grids1`);
`scripts/check.sh --quick` (ctest -LE slow 165/165); `tests/run_tests.sh` 141 passed, 0 failed (with
`extrapolation_grids`); not run: Sphinx (not installed on this host), 1 against 8 threads (budget 3); the
converted-tree gate
(`format3_gate.sh`, tag `gridsg1`): 472 (445) identical, 38 classified, 0 differ, 0 missing, as on the start head
(`gridsg0`); menu 3 of converted T1, T3 (plain conversion: `range … step`) and T2 with a `range … step` extrapolation:
outputs identical to the start engine `work/b-trial5`; `tests/extrapolation_grids/check.sh` (narrow 1 eV level: 378
points, peak 3.5e-8 from the closed form, integral 4.4e-4 from a dense grid, 1 and 3 threads identical with
`AZURE_PLAN_CHECK=1`).

**Open.** An adaptive grid in a pyazr session (refused now); refinement on the S factor rather than the cross
section for `quantity sfactor` (the model value is refined; the floor 1e-3 max|f| leaves sub-barrier tails on the
start grid); `CHECK_DATA` prints the start grid's points (before the refinement); `runall.sh`/`run_one.sh` into the
repository (`tools/`) so that the gate needs no `RUNALL`.

## 2026-11-10 — agent, brief O: the hybrid potential as a value, one e², distorted bound states (`feature/upgrade-4-hybrid-potential`)

**Part 1 — the potential as a value (bit-identical).**
- `physics/HybridPotential.h`: `FormFactor`, `HybridPotential` (`coulomb_radius_fm` 3.5, `central_depth_MeV`, `central`,
  `matching_radius_fm` 15, `step_fm` 0.0275), compared by value; `PotentialKey` (the numbers of a potential for keys);
  `IntegrateInward` (the Numerov solver, a free function; the loop stops at node n_a - 2, the last one read).
- The settings are values in the run's or session's `Config` (`HybridPotentialSettings`, `include/` +
  `src/HybridPotentialSettings.cpp`); `CNuc::Fill` gives each `PPair` its `std::optional<HybridPotential>`; `CoulFunc`
  reads it from the pair.  Chosen over the System because `CoulFunc` is built from a `PPair` at every point and the
  System is rebuilt from `CNuc` on demand; the pair is set-up state owned by the caller's `CNuc`.
- Deleted: `NuclearPotentialManager` (singleton, tags, counter, `exportState`/`importState`), `NuclearPotential.h`
  (class hierarchy), pyazr's `ConfigScope` potential swap, `CoulFunc`'s dead setters (`setMatchingRadius`,
  `setNumerovGridStep`, `setNuclearPotential`, `setUseHybridMethod`).  The channel table's pair signature, the EC term
  key, the radial scratch key, the shift-derivative memo, the THM line-shape memo and the THM vertex memo compare the
  potential's values; the EC signature hashes `PotentialKey` (hybrid EC caches are recomputed once).
- THM: a THM segment whose entrance pair has a potential is refused, `ERROR (E_THM_HYBRID_POTENTIAL)`
  (`ThmVertexAllowsPair`, `src/ThmVertexBoundary.cpp`); no committed project had one.  Test: `tests/thm_options` (j).
- Words: "outside the channel radius" (`CNuc.cpp`, `particle_pairs.rst`), "inward" (`tests/hybrid_potential/tolerance`).
- Verified: 179 × 3 identical to `tab-p3f` (tag `tab-hyb1`, 689 / 662 / 689); `ctest -LE slow` 163/163; pyazr
  `coulomb_functions` (pairs 1, 2, L = 0–3, six energies) and `calculate` of `tests/hybrid_potential` byte-identical to
  `work/b-trial3` (65 arrays); `potential_sessions_test.py` passes and fails with a function-static shared settings
  object injected into `CNuc::Fill`; `physics_hybrid_potential` fails with the step left out of `operator==`;
  `thm_options` (j) fails on `work/b-trial3` (the run goes through); `tests/run_tests.sh` 138 / 0;
  calculations at 1 and 2 threads identical (T1, T2, T3, 13N, hybrid_potential × 3 bases, 15/15); converted-tree gate
  (`convert_spec_equivalence_test`) 183 equal, 3 recorded, 0 differ.

**Part 2 — e² = αħc (outputs of hybrid projects move).** `physics/src/HybridPotential.cpp`: the Coulomb term's
e² = 1.4399764 MeV fm → `kFineStructure * kHbarC_MeV_fm` (8.3 × 10⁻⁶ relative).  179 × 3 against `tab-p3f` (tag
`tab-hyb2`): only `hybrid_potential` and `upgrade/T3` (T3, T3_start) move, in all three bases; everything else
identical.

| project | file | largest relative change | cause |
|---|---|---|---|
| hybrid_potential | AZUREOut_aa=1_R=1.out (p + ¹²C elastic) | 1.3e-5 | e² |
| hybrid_potential | chiSquared.out (segment 3: 2238.93 → 2238.82) | 4.9e-5 | e² |
| hybrid_potential | param.par, parameters.out, scatteringLengths.out | 2.5e-5 | e² (P, S at a_c) |
| T3 | AZUREOut_aa=1_R=1.out (α + α) | 8.8e-2 on δ = 2.26e-3 → 2.48e-3 at E_lab 0.875 MeV, a phase near zero (absolute 2.2e-4); 5.8e-3 next point; others ≤ 1e-5 | e² |
| T3 | AZUREOut_aa=4_R=1.out | 1.6e-5 | e² |
| T3 | AZUREOut_aa=2/3_* (p + ⁷Li, n + ⁷Be: no potential, through the shared levels) | ≤ 7e-7 | e² |
| T3 | chiSquared.out, param.par, parameters.out | 7.5e-5, 2.3e-6, 2.7e-6 | e² |
| T3_start | chiSquared.out / AZUREOut_aa=1_R=1.out | 2.8e-4 / 2.2e-3 | e² |

Independent check: with V_N = 0 the integration now reproduces the Coulomb P and S to ≤ 9.3e-8 (p + ¹²C, l = 0–2,
0.5–3 MeV; `physics_hybrid_potential`, tolerance 2e-7), against 2.3e-5 with the old e² (the test fails then); the
note's DOP853 integration (`work/so-note/radial.py` B) gives 6e-9 / 2e-8.  References regenerated: T3 goldens
(`golden.py --write`, T1 left as recorded; `ENGINE` says so), `tests/hybrid_potential/expected/chiSquared.out`
(segment 3 2239.69 → 2238.82: 2239.69 was a pin of 2026-09-28; the start head gives 2238.93, so 3.4e-4 of the step is
older drift inside the 2e-2 tolerance and 4.9e-5 is e²).  `pyazr.rst` and `hybrid_derivative_test.py` no longer call
1e-5 the integration's accuracy; the latter now asks for 1e-6.

**Part 3 — distorted bound states and closed channels (outputs of `hybrid_potential` move).**
- `physics/HybridPotential.h`: `DecayingSolution` (W̃_l(r) at E = −B ≤ 0: the Whittaker function beyond r_m, the
  zero-energy decaying solution √r K_{2l+1}(√(8cr)) at B = 0, inward by `IntegrateInward`, which now rescales a
  solution past 1e250 into `log_scale`) and `HybridBoundShift` (S = aW̃′/W̃; closer to threshold than η_b = 1000 the
  linear interpolation to S(0), as `coulomb_bound_shift`).  The Coulomb path's reflection S(E) = 2S(0) − S(−E) in the
  band η > 100 is kept and now reflects the distorted S, so S is continuous through threshold (|S(E) − S(0)| ≤ 5e-4
  for |E| ≤ 1 keV, p + ¹²C); no change of the Coulomb path for pairs without a potential.
- One implementation, reached through the two engine classes every user already calls: `ShftFunc` (S, S′, S″ below
  and at threshold → `ChannelFunc`, the channel table's closed rows, Brune/Park transformations) and `WhitFunc`
  (W̃ for the ANC ↔ γ conversion of `CNucLevelParameters`, the EC conversion factor and NF integral of `CNuc`, the
  final-state tail of `ECIntegral`) dispatch on `PPair::GetHybridPotential()`; a pair without one takes the old lines.
- `coulomb_functions` (API, pyazr 2.17.1) returns S at E ≤ 0 and in the band η > 100 (channel radius only; before: 0,
  and in the band the Coulomb routine's value).
- Tests: `physics_hybrid_potential` parts 3–4 against `physics/tests/data/hybrid_bound.csv`
  (`make_hybrid_bound_reference.py`: scipy DOP853 + mpmath Whittaker/Bessel/Coulomb): S on −2 … +0.5 MeV, l = 0–2, to
  6.6e-8 (tolerance 1e-6), V_N = 0 rows (Whittaker/Coulomb S) to 1.1e-7, ln W̃(a) to 4.5e-10, ∫_a^∞ W̃² dr to 2.6e-8;
  fails with E = +B instead of −B in `DecayingSolution`.  `pyazr_hybrid_threshold`: continuity (1e-3) and S at the
  ¹³N ground state against the table (≤ 6.6e-8); fails with the undistorted S injected in `ShftFunc` (steps 1.06,
  0.92, 0.75 for l = 0, 1, 2, the note's 1.07, 0.92, 0.75).
- 179 × 3 against `tab-hyb2` (Part 2): only `hybrid_potential` moves (T3 has no state below the α + α threshold):

| file | largest relative change | what moves physically | independent number |
|---|---|---|---|
| param.par, parameters.out | 0.166 (Park 0.153) | the ¹³N ground state's g_int (ANC 1.62 fm^-1/2 fixed): 1.271943 → 1.060263 MeV^1/2 | γ_obs = C W̃(a)(ħ²/2μa)^1/2 and g = γ_obs(1 − γ_obs² S′)^-1/2 from the DOP853 script: 1.060263 (old W, S′: 1.271943), to all printed digits |
| AZUREOut_aa=1_R=1.out | 1.2e-2 | elastic p + ¹²C through the transformed bound level | follows from the above |
| chiSquared.out | 9.3e-3 (segment 3 2238.82 → 2259.87) | same | — |
| run.log | 0.69 | the χ² lines, and the Coulomb-table lookup count 66 → 20 (the hybrid pair's bound S no longer reads the table) | — |

  Reference regenerated: `tests/hybrid_potential/expected/chiSquared.out`.  Docs: `pyazr.rst` (the "does not touch"
  paragraph replaced), `particle_pairs.rst`, the hybrid README.
- Verified (head of the branch, 3 threads): `ctest -LE slow` 164/164; `tests/run_tests.sh` see the commit message;
  1 and 2 threads identical (T1–T3, 13N, hybrid_potential × 3 bases, 15/15); converted-tree gate 183 equal / 3
  recorded / 0 differ; consistency at every commit.

**Open / next.** The `<potentials>` record of format 3 (note §3.3) and the spin–orbit term (steps 3–5).
`tests/hybrid_potential/check_per_pair.py` (not in ctest) pins χ² = 2462154.2 and fails on the start head's engine
too (2466121.4 there; 2466142.3 now): a stale pin, not touched.  `E_THM_HYBRID_POTENTIAL` is an engine message, not
yet a row of the format-3 diagnostic table (the potential is not a format-3 record yet).

## 2026-11-10 (night) — S14, brief L: ratio to Rutherford and to Mott, the `rutherford` leaf (`feature/upgrade-4-rutherford`)

**Done** (brief `work/briefs/L_rutherford_ratio.md`).
- `physics/CoulombScattering.h`: `rutherford_cross_section`, `mott_cross_section` (b/sr, c.m., with d/dE) from the
  engine's η and k (`ChannelFunctions.h`); Mott with the interference weight (−1)^{2s}/(2s+1); binding
  `_azure2.rutherford_cross_section`, `mott_cross_section`; test `physics_coulomb_scattering` (mpmath table).
- Format 3: `ratio_to_rutherford`, `ratio_to_mott` are built (`Build.cpp` `ReadObservable`: dσ/dΩ at lab or c.m.
  angles with `ReactionSpec::coulombReference`); refused from the pairs (one elastic reaction, charged, identical as
  the quantity says) and for folds the reference does not follow (`BuildEffects.cpp` `CoulombRatioFolds`: Doppler,
  screening, detector attenuation, automatic tolerance; `E_NOT_BUILT`).  Extrapolations alike.  `--describe` says
  what the value is divided by (JSON `coulomb_reference`).
- Tree: `TreeKind::CoulombRatio`, `ratio(dsdo_cm(r), rutherford(r) | mott(r))`; `LeafPart::Rutherford`, `Mott`.  Engine:
  `src/CoulombReference.cpp` evaluates the leaf at the point (and, under a sub-point effect, at the point's own
  sub-points through `EPoint::IntegrateSubPointValues`, with the edge blend): the model is ⟨dσ/dΩ⟩/⟨σ_C⟩.  Energy
  shift: the rule of every leaf (the point is re-energised, d/dshift by finite differences).  The adjoint drops the
  leaf.  `ESegment::IsDimensionlessRatio`: a ratio's value is not converted at load (one rule for RATIO and these).
- α + α checked: the engine's elastic dσ/dΩ of identical particles is already the symmetrised one.
- `pyazr.nds`: `ExforData.format3_quantity()`, `to_dataset()` (format-3 record, ratios unconverted); `to_azr` refuses a
  ratio (no multiplication by σ_R any more).
- Tests: `tests/ratio_to_rutherford/` (check.sh, truth.py, the §9.2 example with `--describe` committed),
  `pyazr_ratio_to_rutherford`, `pyazr_nds_rutherford` (rewritten), cases in `experiment_build`,
  `experiment_prediction`; io_convert leaves a format-3 project alone.  Docs: `examples/ratio_to_rutherford.rst`,
  `data_formats.rst`, `pyazr.rst`, the setup skill's EXFOR reference, the recipe "A new kind of data set" (CLAUDE.md).

**Verified**: build `work/b-ruth` (3 jobs, 0 warnings).  179 × 3 against `tab-p3f`: 689 / 662 / 689 identical (tag
`tab-ruth`;
the only other files are the four new projects of tests/ratio_to_rutherford, absent from the reference); `ctest -LE
slow` 165/165 (2 skipped as on the base: pyazr_predictions, pyazr_composite_gradient), on the tree before the nds and
docs commits, then the nds tests and `pyazr_ratio_to_rutherford` again after them; `tests/run_tests.sh` (OMP 3) 140
passed, 1 failed: cli_help, which passes three times alone (load-dependent, see Open); `tests/format3_engine/check.sh`
PASS; spec equivalence 193 files, 185 converted, 183 equal, 2 recorded, 0 differ (the 4 new format-3 files "already
format 3"); T1, T2, T3, 13N × 3 bases at 1 and 2 threads: 12/12 identical; consistency checks pass at every commit.
Physics: hard sphere α+¹²C and α+α against mpmath 3.5e-10 / 5.1e-10 (11-digit output); χ² ratio = χ² multiplied
9.3e-15; 5 keV shift: both forms = truth.py to 6e-10, difference 0.39977 on 12.5; beam fold: ratio of folds to
1.8e-6, fold of the ratio off by up to 7e-3; gradient = central differences 5e-7 (plain), 8e-8 (beam); CLI grad
check 4.8e-8.  Every new test shown failing once (injections in the report, work/ruth/).

**Open.** Lab-frame Legendre ratios (`LEG/RTH`) and the electron-scattering Mott cross section are out of scope.
A detector-angle fold (finite angular acceptance) does not exist in the engine; the attenuation coefficients are
refused for these quantities.  Sphinx not run (not installed here).  `tests/cli_help` failed once inside run_tests
under machine load (flags reported "compiled out", then the doc table check), passed alone three times: not caused
by this branch (no CLI or doc table change), not investigated.

## 2026-11-10 (agent, brief Q) — `EData` split W1 and W2 (branch `feature/upgrade-4-edata-split-2`)

**Done** (`edata_split.md` §5.2; one commit per class for W1, one for W2; each item re-proved dead by `git grep -nw`
over `src/ include/ api/ experiment/ io/ physics/ convert/ pyazr/ tests/ tools/` on this tree).
- W1 `EData` (530d5a7f): `Iterations`, `ResetIterations`,
`Get{Screening,Density,Doppler,Thickness,ThmCoherent}ParamOffset`,
  `DensityParamOfEffect`, `ResultsSize`, the write-only `componentError_`; the unreachable `data == nullptr` branch of
  `FillEnergyShiftsFromParamsImpl` with the `data` parameter of `FillEnergyShiftsFromParams` (24 callers passed the
  data set itself; the kept branch reads `this`) and its write-only counter `k`.  §7 items 1, 2, 10 (EData parts).
- W1 `EPoint` (685c2697): `StoreSubpointOffsets`, `ApplySubpointShift` (never defined), `IsTargetFallbackOnly` (the flag
  stays; comments name `SetTargetFallbackOnly`), `ModelSlot`, `ModelAngularSlots`, `SetExcitationEnergy`,
  `angleKinFactor_` with its getter, setter and writes; not in the note: `GetNumAngularDists`, `GetAngularDist`.
  §7 item 6: the three self-assignments through `SetCrossSectionKinFactor` and with them the setter.
- W1 `ESegment` (2518a118): `GetNominalEnergyShiftSqrt` / `energyShiftSqrtNominal_`, `Has{Shift,Norm}GroupFollowers`
  with their members and setters (and the two calls in `EData`), `SetVaryNorm`, `AddComponent`,
  `componentCalculationMutex_`; §7 items 3, 10.
- W2 (3837d59c): `CovarianceFileFor` … `PrepareCovarianceSeeds` (256 lines) verbatim to `src/EDataCovariance.cpp`;
  one check-11 row.
Sizes: `EData.cpp` 4,793 → 4,455, `EPoint.cpp` 2,983 → 2,920, `ESegment.cpp` 1,406 → 1,366, `EData.h` 576 → 560,
`EPoint.h` 911 → 896, `ESegment.h` 590 → 574; new `EDataCovariance.cpp` 275.

**Verified** (build `work/b-esplit2`, 3 threads): every commit 0 warnings, consistency (after `--tighten-sizes`),
`scripts/check.sh --quick` (ctest -LE slow 163/163), 179 × 3 identical to `tab-p3f` (689/662/689; tags
`tab-esplit2-w1a`, `-w1b`, `-w1c`, `-w2`); after the EData commit `tests/shift_fn` and `tests/energy_shift_sqrt` pass;
W2: the move `diff`s empty, and the covariance projects (`data_covariance` diag, dense, rel; T2 Brune and formal;
fit.azr with MIGRAD-gradient, LM, GSL-LM fits) byte for byte against the base build `work/b-trial5`, all output
files (`work/esplit2/cov_cmp.sh`).  At the end (`work/esplit2/final.sh`): `tests/run_tests.sh` (OMP 3) 140 passed,
0 failed; 1 vs 2 threads on T1, T2, T3, 13N in the three bases: 12/12 output directories identical.

**Open.** History comments outside §7 item 10's list remain (`EData::CheckComponentPairs` "as was done",
`EPoint::ConvertLabAngle` "Until 2026-10-29", `ConvertCrossSectionGammas` "Before 2026-10-23", two "Before 2026-11-05"
comments in `EData.cpp`
(the grid re-placement, the gamma pair's final level).  `#include <mutex>` stays in `ESegment.h` (other files may rely
on it).

## 2026-11-10 (director, afternoon) — Phase 4: the fourth wave merged (expression tree, plan, `EData` split N1–N4, the engine runs format 3)

**Merged into `feature/upgrade-4-data`** (`--no-ff`, in this order; each rehearsed, built and verified in a scratch
worktree, the committed tree is that tree): H `feature/upgrade-4-expr` (55fa32ed: S10 the model pass evaluates the
expression tree, S11 the gradient through `backward`, `--print-predictions`); I `feature/upgrade-4-plan` (8a36a61a: the
model-pass follow-ups and S5, the calculation plan); P `feature/upgrade-4-edata-split-1` (9dff49fe: `EData` split
N1–N4); G3 `feature/upgrade-4-format3-reader` (28dbb50a: one `<levels>` table, parameter names from the specs, the
engine
runs a format-3 project, `--describe`). Their entries are below.

**The merged head as a whole** (28dbb50a; build `work/b-trial5`, 0 warnings): 179 × 3 identical to the Phase 3
head (689 + 662 + 689 files); `ctest -LE slow` 163/163; `tests/run_tests.sh` 140 passed / 0 failed; calculations at 1
and
3 threads 12/12 identical; T2 with a `<modelData>` section identical at 1 and 3 threads in the three bases; spec
equivalence: 189 files, 185 converted, 183 equal, 2 recorded, 0 differ; the converted tree (exact grids, menu 1): of
the 510 (483 formal) output files per basis 472 (445) identical, 38 classified, 0 differ, 0 missing. Not yet run on
this head: 1 against 8 threads (the machine is shared with three agents).

**What the director changed at the merges** (each in its merge commit's message):
- pyazr is **2.18.0** (H and G3 both took 2.17.0): `predictions` (H) and `describe` (G3) are both exported
  (`pyazr/__init__.py`, `pyazr/schema.py`).
- `src/ProjectCommands.cpp` holds both commands (`--print-predictions`, `--describe`); the consistency test's list of
  flags handled in `main` has both; `command_line.rst` has both rows.
- `tests/consistency/file_sizes.tsv` at the merged counts: `AZUREAPI.cpp` 1837, `AZUREGrad.cpp` 1850, `CNuc.cpp` 1894,
  `EData.cpp` 4793, `EPoint.cpp` 2983, `ESegment.cpp` 1406; `EDataThm.cpp` (773) is off the record.
- Checked on the combined tree, because no branch could: `--print-predictions` of the converted T2 gives the 13
  lines of the old file; the converted T2 runs with `AZURE_PLAN_CHECK=1`; its outputs are the old file's but for the
  two classified files (`param.par`, `normalizations.out`).

**Decisions of the director.**
- I: the plan's task becomes the unit of work (a task evaluates its host and writes its mapped points without going
  through the data set) in **Phase 6**, with the mutable statics, not now.
- P: `EData::ThmGroupEnergyRange` as a public static member is accepted (its callers are free functions of the
  set-up); it goes into the THM set-up's own file with step M of `edata_split.md`.
- G3, on its report (the agent implemented them; commits in its entry):
  1. The converter keeps every abundance row of a sample that a line uses, also a row no line reaches (the isotope
     from the `<levels>` pairs with the sample's light particle; none or several is a question). The engine's rule
     that a nucleus a sample does not list counts 1.0 in the sum-to-one prior is unchanged and **open for the targets
     step**.
  2. A run that uses no data (menu 3, extrapolation-only pyazr calls) does not need the data files:
     `BuildOptions::dataFiles` (`Required | NotUsed`), `W_DATA_FILE_NOT_READ`. `--validate` and `--describe` stay
     strict.
  3. An unreadable data file refuses a data run (`E_DATA_FILE`); the old engine skipped the segment. Documented.
  4. The start step of a fixed zero shift differs between the formats (`param.par` only): accepted as classified.
- The converted-tree gate lives only in `work/fmt3c` (`gate.sh`, `mktree.sh`, `gatecmp.py`, `classify3.py`); it
  becomes a committed tool with brief K (the director ran it with `work/director/gate.sh`).

**Defects and debts found, not yet fixed** (owner in brackets):
- An `AZURE_PLAN_CHECK=1` breach throws `std::logic_error` that nothing catches (abort without a message naming the
  task) [Phase 6, with the task as unit of work].
- Peak memory did not drop with S5 (the memory is in the tables, not in the points) [S7 measures again; Phase 6].
- `--print-predictions` names a data set by its index (`data 3:`) also in a format-3 project, where it has a name
  [S13/S14, with `predict`].
- `ProjectSetup` re-reads the old sections once per question; `Format3ThmSection` writes old `<thm>` grammar for
  the parser; level parameters are `energy_<n>` / `width_<n>_<m>` in format 3 [3d; S7/3d; the `<levels>` step].
- The two T2 benchmarks of H and I differ by a factor 4 in absolute time (0.29 s against 1.1 s for a calculation):
  each compared with its own start head and found no change, but the conditions differed (not recorded). The
  director measures the merged head against the Phase 3 head when the machine is free, with 1 against 8 threads
  [director, end of wave].
- Not run by anyone: ThreadSanitizer (no runtime on this host), `ctest -L slow`.
- Run-to-run variation of the Minuit call count at 2 threads (`energy_shift_sqrt` 40/41, `shift_fn` 36/37 calls,
  same fitted outputs) is cause B of `fit_threads.md` (the PI: Minuit2 is not patched).

**For the PI to see** (behaviour a user can meet; additions):
- The API's `calculated_segments` of an E1-only / E2-only segment is its E1 / E2 part (confirmed earlier).
- Format 3: a data run with an unreadable data file is refused (the old engine skipped the segment with a warning);
  an extrapolation run does not need the data files; a sample keeps all its isotopes.
- `AZURE2 --describe FILE [--json]` / `pyazr.describe`: how the engine understood a project (pairs, data files as
  opened, effects, parameter names). `AZURE2 --print-predictions FILE` / `pyazr.predictions`: the expression each data
  set is compared with.
- The gradient of an E1-only / E2-only SUM with a component was wrong before S11 (0.29 and 1.0 of the largest
  element against central differences); it is right now (1e-7). Fits with `--use-gradient` of such segments change.

**Running now** (Opus agents, 3 threads each; briefs in `/data0/skowrons/AZURE2-upgrade/work/briefs/`):
- O, `feature/upgrade-4-hybrid-potential` (`wt-hybrid`), `O_hybrid_potential.md`: the potential as an immutable value
  in `physics/` (no manager singleton), e² = α ħc, the distorted decaying solution (the PI's Q4, Q5; moved outputs
  classified).
- Q, `feature/upgrade-4-edata-split-2` (`wt-esplit2`), `Q_edata_split_2.md`: W1 (dead code out of `EData`, `EPoint`,
  `ESegment`), W2 (the data covariance into `src/EDataCovariance.cpp`).
- L, `feature/upgrade-4-rutherford` (`wt-ruth`), `L_rutherford_ratio.md`: S14, the `rutherford` / `mott` leaves and the
  quantities `ratio_to_rutherford`, `ratio_to_mott`.

**The PI on the split, 2026-11-10:** he expected `EData` to be split directly. The director agrees and changes the
order of `edata_split.md` §5: the cuts scheduled "after S7" (M1–M3) and under D1 are made now, as verbatim moves,
right after Q (brief `S_edata_split_3.md`); S7 and Phases 5–6 then edit the small files. F1/F2 still wait for 3d.

**Then.** S (the remaining mechanical cuts of `EData`, `EPoint`, `ESegment`, after Q), K (extrapolation grids: list,
file, points/log, c.m., adaptive; commits the gate), J (units; S factor with the `nds.py` mass fix; yields per
charge), R (S7 `Measurement`, after S), 3c, 3d, S12, S13, S16; the spin-orbit term
(note §5 steps 3–5, after O); `EData` split F1/F2, M, D1; CNuc 3b (iv), 6b, 7; Phase 5 (with M3 (b)); Phase 6 (mutable
statics, the task as unit of work). Merge `origin/dev` only when the reconstruction is complete.

## 2026-11-10 (evening) — format 3, step 3b-1, G-3: the engine runs a format-3 project (`feature/upgrade-4-format3-reader`)

**Done** (brief G-3; report `work/fmt3c/REPORT.md`).
- **A** (b2bb1452): one table of the pairs and levels of `<levels>` (`experiment/LevelTables.h`: `PairInfo`,
`LevelInfo`,
  `PairLabel`, `PairsNamed`, one wording of `E_UNKNOWN_PAIR` / `E_AMBIGUOUS_PAIR`), made once from the rows
  (`LevelTablesOfRows`, reader library) and read by `Build`, `BuildModelData` and the engine's model data (which
checks a
  datum's level indices against the compound nucleus); `PairRow`, the engine's `PairInfos/LevelInfos` and the
converter's
  label copy are gone; `io::IdenticalParticles` is PPair's rule.  `Schema::Build` resolves `<modelData>`
  (`BuiltProject::modelData`); the spec-equivalence test compares the model data.
- **B** (ba0cd97d): parameter names come from the spec (`ParameterNames`; adaptor: `segment_<k>_...`; Build:
`<ds>.norm`,
  `<ds>.shift.offset|sqrt`, `D.value`, `G.offset`, `G.a<j>`, `<target>.areal_density`, a follower's
`<record>.<field>[<ds>]`;
  effects and abundances by their target).  No fixed name buffer (`EData`, `CNuc`).  The converter's map names every
  parameter the engine makes and is injective (checked over every converted file).  Converter: one path for data and
test
  lines with a sample (M_thick asks 4 questions now; no longer a recorded difference); a missing sample of a test line
is
  noted.
- **C**: `Config::ReadConfigFile` reads a format-3 file (records → Validate → `<levels>` tables → Build;
`Config::format3`);
  `src/ProjectSetup.cpp` (was `OldLineSetup.cpp`) gives the set-up the specs of either producer, the `<levels>` rows
read
  once and the exit pairs of the active data sets (`ExitPairKeysInFile` deleted; `AZURE2.cpp`'s pair and THM scans read
  specs too); `<thm>` membership and weights from the records; `AZURE2 --describe FILE [--json]`, `pyazr.describe`
  (pyazr 2.17.0); `tests/format3_engine/check.sh`; docs (`data_formats.rst` format-3 section, `command_line.rst`, the
  note's status, the recipe "A new kind of data set").  Two converter/Build defects found by the gate and fixed: a
  `quadrature` without a straggling flag lost its width multiplier (TargetEffect gated it on the old line's flag), and
the
  converter wrote `width_multiplier 20` where the old line gave none (a pure convolution's core then changed).

**Verified** (build `work/b-fmt3`, 0 warnings, 3 threads): 179 × 3 on the old tree identical to `tab-p3f` after A, B and
with the final engine (689/662/689); `ctest -LE slow` 155/155 (final tree); `tests/run_tests.sh` 138 passed, 0 failed
(with `format3_engine`); consistency checks pass; spec equivalence 189 files, 185 converted, 183 equal, 2 recorded,
0 differ (injections =2, =3 fail). Gate (tree converted into `work/fmt3c/tree3x`, exact grids; `tree3p` plain): of the
510
output files per base of `tab-p3f`, 471 identical, 30 classified (param.par steps of fixed zero shifts, follower rows,
angle window of angle-integrated data, sample header names, scattering lengths as model data), 2 differ (angle_mapping,
the abundance defect), 7 missing (shift_subpoint_angle/template: refused for missing data files); plain grids give the
same in menu 1 and differ in menu 3 in extrapolation outputs only. 1 vs 2 threads on converted T1, T2, T3, 13N
identical. The five refusals of §8.5 hold (three checked with dummy data files).

**Follow-up (the director's decisions on the G-3 report).** 1 (951fd09f): the converter keeps a sample's row of a
compound nucleus no line reaches (its isotope from the <levels> pairs with the lines' light particle; none or several
is a question), Build gives it its abundance row: tests/angle_mapping is now identical but for name classes.
2 (4e6e7f59): `BuildOptions::dataFiles` (`DataFiles::Required | NotUsed`); the engine reads a format-3 project with
NotUsed and `EData::Fill` refuses an unread data file (`E_DATA_FILE`), `EData::MakePoints` warns
(`W_DATA_FILE_NOT_READ`); `--validate` / `--describe` stay Required. 3 (9719363f): an unreadable data file refuses a
data run; documented with the other differences from the old engine. 5 (597c7aad): §3.11's T2 quadratures.
Verified after them (engine `work/fmt3c/binC/AZURE2` = `work/b-fmt3`): 179 × 3 identical to `tab-p3f` (`tab-fmt3cH`);
`ctest -LE slow` 155/155; `tests/run_tests.sh` 138/0; spec equivalence 189/185/183/2/0; gate menu 1, both grids: of 510
(483 no-brune) files 472 (445) identical, 38 classified, 0 differ, 0 missing; menu 3 exact: 302 identical, 9
classified, 8 differ (bad_sample runs; modelData.out added); 1 vs 2 threads on converted T1, T2, T3, 13N identical.

**Open.**
- The engine counts a nucleus a sample does not list at abundance 1.0 in the sum-to-one prior (old behaviour,
unchanged).
  It showed in tests/angle_mapping before decision 1 (format 3: 7_13 at 1.0, Total-Abundance chi^2 25 instead of 0);
  the targets step decides whether an unlisted nucleus counts.
- A member of a varied shared normalization has a fixed parameter of its own also where its old line did not vary (S7
  removes it).

**Debt** (each removed by the step named):
- `src/ProjectSetup.cpp` reads the old sections through the adaptor once per question (`ProjectExitPairKeys`,
  `ProjectReactionPairs`, `ProjectHasThmData`, `EData::Fill`, `MakePoints`) — step 3d deletes the old branch.
- `Format3ThmSection` (src/ConfigThm.cpp) writes old-grammar `segments=` / `weight[k]=` text for the <thm> parser — goes
  when `ThmOptions` takes `ThmMembershipSpec` directly (the THM part of S7/3d).
- Level parameters are still named `energy_<n>` / `width_<n>_<m>` in a format-3 run: `project_format3.md` is silent
  (`<levels>` keeps its grammar, §3.2), so the names of level parameters are open for the `<levels>` step.

## 2026-11-10 (agent, brief P) — `EData` split N1–N4 (branch `feature/upgrade-4-edata-split-1`)

**Done** (`edata_split.md` §5.1, one commit per step):
- N1: the output writer (`WriteNaturalDecomposition`, `WriteBandLine`, `WriteOutputFiles`) to
`src/EDataOutputFiles.cpp`;
  four check-11 allowance rows for the penalty totals of `chiSquared.out` (D2).
- N2: the four data check writers (`PrintData`, `PrintLegendreP`, `PrintEDependentValues`, `PrintCoulombAmplitude`) to
  `src/EDataCheckFiles.cpp`; new ctest `check_files_data` (`tests/check_files_data/check.sh`: T1 with the four
  checks on, files and line counts).
- N3: `WriteThmExperiments` and the three THM tables to `src/EDataThmReports.cpp`; `ThmGroupEnergyRange` is a public
  static member of `EData` (its set-up callers are free functions, so it cannot be private as the note says).
- N4: `WriteOutputFiles` calls `WriteSegmentFiles`, `WriteChiSquaredFile`, the THM writer, `WriteNormalizationsFile`,
  `WriteShiftsFile`; the per-segment χ² rows have their own loop over the same segments.
Sizes: `EData.cpp` 5,571 → 4,841, `EDataThm.cpp` 1,233 → 773; new files 522, 251, 487 lines.

**Verified how.** Each move: the moved lines of the old file `diff` empty against the new file's body. Each step:
build 0 warnings, 179 × 3 identical to `tab-p3f`, `scripts/check.sh --quick`. Outputs the 179 runs do not reach,
base (`work/b-trial2`) against branch, byte for byte: `.band`, lab-frame and `shifts.out` files of seven `check.sh`
(N1, N4); the four check files as `file` and `screen` on T1, T3, 13N, `natural_neon` and a line-shape THM project
(N2); the THM tables through pyazr with `repr` on 38 projects (N3). At the end: `tests/run_tests.sh` 138 passed
(OMP 3), T1–T3 and 13N × 3 bases identical at 1 and 2 threads. Scripts in `work/esplit/`.

**Open.** N3's helper is public, not private (above). The base engine itself varies from run to run at 2 threads in
the Minuit call count of two fits (`energy_shift_sqrt` 40/41, `shift_fn` 36/37 calls, `fitStatus.out`), with the same
fitted outputs; not investigated.

## 2026-11-10 (evening) — brief I: the follow-ups of the model pass, and S5, the calculation plan (`feature/upgrade-4-plan`)

**Done** (brief `work/briefs/I_plan.md`; full report `work/plan/REPORT.md`).
- **Fixed energy shift back to 0** (934c11a3): a segment whose points carry a shift of an earlier evaluation stays
  "shift active", so a fixed shift set back to 0 re-energises its points; `upgrade_history` gains vector C (every fixed
  shift and sqrt(E) coefficient moved by 1 keV): pooled C → A = fresh A, C → A → C = fresh C (8 failures before).
- **API value of E1-only / E2-only segments** (37aac384): `UpdateSegments*` report the model value of every segment
  (one loop, `AZUREAPI::ReportSegments`); pyazr's `_segment_values` no longer swaps the part in itself. Test
  `pyazr_api_model_value` (8 failures on the start head's engine).
- **MCMC snapshot** (3af76e1d): `tests/mcmc_snapshot_components/check.sh` — total, E1-only and E2-only capture blocks of
  the snapshot equal a calculation at `param.mcmc` (8.8e-8); fails with the leaf value injected.
- **Covariance seeds** (87ef450f): `Results::ObjectiveValues::covarianceSeed` (indexed as `Results::model`), not
  `ESegment`; `tests/data_covariance` (d) and `pyazr_data_covariance` fail with the neighbour's seed injected.
- **S5, the calculation plan** (9df9fc92): `include/CalculationPlan.h`, `src/CalculationPlan.cpp`; a point holds its
task
  index, never a pointer to another point; `local_mapped_points_` and `Clone`'s pointer re-mapping are gone; the plan is
  built and checked (one writer per point) in `MapData`; `AZURE_PLAN_CHECK=1` checks each task's slot writes at run
  time; tests `upgrade_plan` (179 projects' set-up), `upgrade_plan_runtime`; check 29.

**Still written into the data set by an evaluation** (after the seeds moved): the hand-over of the parameters (norms,
shifts with their re-energised points, grids and the plan, densities, screening, Doppler T, thicknesses); in the point
pass only the external-capture memo `EPoint::ec_shifted_` (H6, keyed by exact energy) and the scoped flags
`DopplerAppliedScope` / `TargetFallbackScope` (set and reset inside one call); the relvector covariance model frozen
at the first `EnsureFactored` (H9, a shared group, once).  Outside the data set: `CNuc::ObservedOfFitValues`.

**The two projects without data** (`tests/sum_lab_frame/models/sum_data.azr`, `tests/ratio_frame/models/ratio_data.azr`,
and their `*_free.azr`): their data are made at test time — `sum.dat` by `tests/sum_lab_frame/check.sh` (b) from the
plain segments' fitted lab values, `ratio.dat` by `tests/ratio_frame/check.sh` (e); `tests/failed_point/check.sh` and
`tests/pyazr/composite_terms_test.py` write five points of their own.  Nothing to commit; in the 179 × 3 sweep they run
without data in both trees (READMEs say so).

**For the PI to see** (behaviour change): the API's `calculated_segments` of an E1-only / E2-only segment is its E1 /
E2 part (the value its χ² uses), no longer the total; `calculate_sfactor*` of such a segment changes with it.

**Verified** (3 threads, `work/b-plan`, 0 warnings, head 9df9fc92; `work/plan/verify.sh`): 179 × 3 identical to
`tab-p3f` (689 / 662 / 689; tag `tab-plan2`); `ctest -LE slow` 158/158; `tests/run_tests.sh` 138 passed / 0 failed;
calculations at 1 and 3 threads identical on T1, T2, T3, 13N in the three bases (12/12); T1 and T3 short fits
(MIGRAD and `--use-gradient`, 150 calls, 1 thread) every output file identical to the start head (logs differ only in
progress bars). Benchmark T2 (3 runs each, this host): calc 1.08–1.16 s → 1.13–1.14 s (1 thread), 0.45–0.47 →
0.44–0.47 s (3); fit200 68.9–71.0 → 69.9–71.4 s (1), 47.7–49.6 → 48.4–49.5 s (3); peak RSS of the fit 272–273 →
270–274 MB (1), 348–350 → 349–353 MB (3): unchanged, not the drop expected (`sizeof(EPoint)` 840 → 824 B is
5 kB on T2's 330 points; the memory is in the tables).  Not run: ThreadSanitizer (no TSan runtime on this host:
`libtsan.so` missing, no clang), 1 vs 8 threads, `ctest -L slow`.

**Open / next.** S7 (`Measurement`); the plan's tasks still call `EPoint::Calculate` with the mapped points resolved
through the host's data set (Phase 6 makes the task the unit of work); the dead `data == nullptr` branch of
`EData::FillEnergyShiftsFromParamsImpl` (no caller) could go.

## 2026-11-10 (evening) — Phase 4 S10 and S11: the model pass evaluates the expression tree; the gradient through `backward` (`feature/upgrade-4-expr`)

**Done** (brief H; `data_model.md` S10, S11, M7).
- **S10.** `experiment/Prediction.h`: one builder per kind of tree (observable, SUM — also a total capture —, RATIO,
  transmission, yield) in the table `TreeKinds()` (`TreeKindOf` decides), from a `PredictionSetup` (the reaction, the
  terms the set-up made: component segments with scaling, fixed angle, compound nucleus; sample; thickness); the tree's
  leaves, factors and parameters listed beside it; `PredictionText`.  `Expression.h`: `Inputs`/`Cotangents` are views
  of the caller's arrays, the workspace holds a tree of up to 32 nodes at one point in fixed arrays (no allocation).
  Engine: `EData::BuildPredictions` (after Fill / MakePoints, `src/Predictions.cpp`) gives every segment an immutable
  tree shared by the working copies; `ComposeModelValue` (src/ModelPass.cpp) fills the inputs at the point and evaluates
  it; the hand-written SUM / RATIO / sample / transform arithmetic is gone.  `AZURE2 --print-predictions FILE`,
  `pyazr.predictions(path)` (2.17.0); ctest `upgrade_predictions` (189 projects of tests/ and examples/ against
  `tests/upgrade/predictions.txt`; regenerate with `python3 tests/upgrade/predictions_test.py --binary <AZURE2>
--write`),
  `pyazr_predictions`, `experiment_prediction`.  Check 28 extended (the pieces of a composition are read only where the
  tree's inputs or setup are made).
- **S11.** `azure2::ModelPointBackward` (cotangents of leaves, sub-points and parameters from the tree's `Backward`);
  every route of `src/AZUREGrad.cpp` uses it (plain, shared-component, per-term, transform, grouped).
  `GradTransformedPoint`, the SUM/RATIO coefficients, the abundance weights, `ESegment::TransformTerms`, `Transmission`,
  `YieldY0` deleted; AZUREGrad.cpp 2036 → 1850, ESegment.cpp 1446 → 1406.

**Verified** (work/b-expr, 0 warnings, 3 threads, `AZURE2_COULOMB_TABLE` set): 179 × 3 identical to `tab-p3f` after S10
(tag `tab-expr-s10`) and after S11 (`tab-expr-s11`), 689 / 662 / 689; `ctest -LE slow` 158/158 after each part
(poison and history tests pass); `tests/run_tests.sh` (final tree, OMP 3) 137 passed / 0 failed; calculations at 1 and 2
threads on T1–T3 and 13N in three bases 12/12 identical; S11 gate against the start head (pyazr, 1 thread):
`objective_and_grad`, `chi2_and_grad`, `residual_jacobian`, `model_gradients` of T1, T2, T3, natural_neon, two
transmission/yield models (with and without a spread) and a SUM project, 75 arrays all bit for bit; T1 and T3 fits
`--use-gradient --max-calls=150` every output identical; `AZURE_GRAD_CHECK` of tests/transmission_yield unchanged
(3.05e-11, 3.54e-11 of the largest); T2 bench (1 thread, 3 runs each) calc 0.28–0.30 s vs 0.29–0.30 s, fit200
17.4–18.0 s vs 17.0–18.2 s, peak RSS unchanged; malloc calls of a T2 calculation +97 (the 13 trees' set-up, nothing
per point).  Scripts: `work/expr/final.sh`, `grad_dump.py`, `grad_cmp.py`, `fits.sh`.

**Follow-up (director's request): the three S11 changes off the test set, shown.** (1) An E1-only / E2-only SUM with
a component: the forward reads the leading term's E1 / E2 part and each component's full sigma; ctest
`pyazr_composite_gradient` (13N levels, synthetic data) compares the gradient with central differences by both
component routes: 1.1e-7 / 8.9e-8 of the largest with this branch, 0.29 (E1) and 1.0 (E2) with the start head, whose
adjoint took the component's E1 / E2 part.  (2) A RATIO with a zero numerator: `experiment_expression` checks
`RatioBackward` against the closed form (J_den/J_num)/den; no project shows it in a gradient, because a numerator that
is exactly 0 (closed exit channel, zero width) has zero derivative there.  (3) A non-finite SUM term: the forward
propagates it (NaN sum) and the tree's adjoint does not skip it (unit test in `experiment_expression`); the old adjoint
skipped it while the forward did not.  The head of `feature/upgrade-4-data` has moved (docs, `include/Constants.h` named
masses); this branch is not rebased on it (the director merges).

**Open.** The tree holds a node per term (T2 data 8: 4 nodes); the transform path still allocates its sub-point
vectors and the kernel object per point, as before.  Not run: 8 threads, TSan, `ctest -L slow`, sphinx.  Next: S12
(components as leaves), S13 (filters), S14 (`predict`, `rutherford`).

## 2026-11-10 — plan for splitting `EData`, `EPoint`, `ESegment` (branch `docs/upgrade-edata-split`, documents only)

**Done.** `docs/upgrade/edata_split.md`: inventory at 070756a1 of every member function (`EData` 139, `EPoint` 226,
`ESegment` 187 definitions, overloads apart) and data member (43, 93, 98) with lines, responsibility, callers and
destination; a partition of the four source files that sums to their sizes; target files with sizes; what waits for
S5, S7, S10–S12, 3d and Phases 5–6; ordered steps; risks; defects met. A paragraph in `ROADMAP.md` names the track.
The first three steps touch code no running branch edits (checked against `wt-expr`, `wt-plan`, `wt-fmt3` at 10:42):
N1 the output writer to `src/EDataOutputFiles.cpp` (490 lines; `EData.cpp` 5,571 → 5,081; check 11 needs four
allowance rows), N2 the check files to `src/EDataCheckFiles.cpp` (240; → 4,841; no test project writes them: a
by-hand comparison and a new ctest), N3 the THM reports and tables to `src/EDataThmReports.cpp` (461; `EDataThm.cpp`
1,233 → 772).

**Verified how.** All numbers counted by the scripts in `/data0/skowrons/AZURE2-upgrade/work/edata-note/` (brace
matching, name search over `src/`, `include/`, `api/`, `tests/`); the partitions sum back to 5,571, 1,233, 3,006 and
1,446 lines. `python3 tests/consistency/repo_consistency_test.py` passes. Nothing was built or run: the steps are a
plan.

**Open.** Decisions D1–D5 of `edata_split.md` §8 (D1: whether to cut the parts Phases 5 and 6 rewrite into files
before those phases; without it `EData.cpp` stays at about 1,980 and `EPoint.cpp` at about 2,720 lines through Phase 5).
Defects found, not fixed: §7 of the note (twelve items) and the duplicates of §2.6 (eleven).

## 2026-11-10 (agent, brief N) — design note: a spin–orbit term in the hybrid potentials

**Done** (branch `docs/upgrade-spin-orbit-note`, documents only): `docs/upgrade/spin_orbit_potential.md`. It covers
today's hybrid model with file:line (§1) and the physics of a spin–orbit term in AZURE2's channel-spin channels (§2).
With one spinning particle the term is diagonal. With two, the channel spins of one (α, l, J) are coupled by the 6-j
recoupling of `physics/src/Recoupling.cpp`, the channel functions become symmetric blocks, and U keeps the form of
Eq. conv-u with blocks. The note also gives the design of option (a) (§3), tests (§4), steps (§5) and eight questions
for the PI (§6). One line is added to `cnuc_split.md` at step 3b (iii).

**Verified how**: scripts in `work/so-note/`:
- the angular m-scheme checks, 308 cases to 5 × 10⁻¹⁵;
- the engine's Numerov method against DOP853 with a spin–orbit term, to 10⁻⁷;
- a δ-shell closed form;
- a coupled 2 × 2 integration and the unitarity of the block U;
- the reference engine `work/b-trial2` against the Numerov transcription, to 10⁻⁸.

No engine code was written or built.

**Found** (not fixed):
- S of a hybrid pair jumps at threshold: it is distorted above threshold and the plain Whittaker form below. The jump
  is about 1 for `tests/hybrid_potential`.
- The hybrid Coulomb potential uses e² = 1.4399764 MeV fm, not α ħc. This alone accounts for the "10⁻⁵ accuracy of
  the integration" that `pyazr.rst` and `hybrid_derivative_test.py` quote.
- Two wordings are wrong: "inside the channel radius" (`CNuc.cpp:203`, `particle_pairs.rst:439`) and "outward"
  (`CoulFunc.h:100`, `tests/hybrid_potential/tolerance:4`).

**Open**: the PI's answers to §6 of the note; then steps 1–5 of its §5.

## 2026-11-10 (director, day) — Phase 4: the third wave merged (specs, `Schema::Build`, model data, one model pass)

**Merged into `feature/upgrade-4-data`** (`--no-ff`, in this order): G1 `feature/upgrade-4-format3-reader` at 7342c790
(d1711801: the set-up consumes specs, S9); E `feature/upgrade-4-model-data` (5d96d069: `<modelData>`, `ModelDataTerm`,
S15); F `feature/upgrade-4-model-pass` (f7463bd8: S1, the rest of S2, S4); G2 `feature/upgrade-4-format3-reader` at
2f57af31 (fe181256: `Schema::Build`). Their entries are below. **Procedure since the second wave:** every merge is
rehearsed in a scratch worktree, the merged tree is built and verified there, and the tree committed is that tree
(`git diff --cached` against the rehearsal is empty), so no merge commit is unbuilt.

**The merged head as a whole** (fe181256; build `work/b-trial2`, 0 warnings): 179 × 3 identical to the Phase 3 head
(689 + 662 + 689 files); `ctest -LE slow` 155/155; `tests/run_tests.sh` 137 passed / 0 failed; calculations at 1 and 8
threads
12/12 identical; T2 with a `<modelData>` section identical at 1 and 8 threads in the three bases and identical to the
tree
before the model pass; spec-equivalence: 189 files, 186 converted, 183 equal, 3 recorded, 0 differ (one value moved
by one ulp makes 119 files fail).

**What the director changed at the merges** (each in its merge commit's message):
- Consistency checks: **27** is the model data's (every datum quantity is evaluated and tested), **28** the model
  pass's (one composition, one residual; it was 27 on its branch).
- `ModelDataTerm` writes its totals to `Results::objective` (`scatteringLength`; `modelData`, `hasModelData` are new
  members); `chiSquared.out` reads them there. Nothing the objective computes is stored in the data set any more.
- The diagnostic codes of the model data and of `Build` are one list (63 codes); `E_UNKNOWN_PAIR` and
  `E_AMBIGUOUS_PAIR` exist once and have a case in both tests.
- `ESegment` lost both the stored χ² (F) and the grid steps (G2: an extrapolation's grid is a list).
- `tests/schema/schema_v3.json` and `docs/source/reference/format3.rst` regenerated from the merged build;
  `tests/consistency/file_sizes.tsv` set to the merged counts at each merge (lowered: `AZUREAPI.cpp` 2033,
  `AZURECalc.cpp` 2075, `AZURECalcMCMC.cpp` 1434, `EData.cpp` 5571, `ESegment.cpp` 1446).

**Decisions of the director** (the briefs of the next wave carry them):
- G1: a composite block whose operation type is neither sum nor ratio is refused (the old engine ran it
  inconsistently). Shared normalization and energy scale stay group labels with the engine's leader rule until S7.
- E: the gradient of the level quantities (ωγ, widths, energies, branching, ANC) is a central difference of the
  observed values (about 70 µs per evaluation on T2 against 77 ms for the data); the scattering lengths keep their
  exact adjoint. `resonance_strength` and `branching_ratio` refuse a `[total]` exit.
- F: the API's `UpdateSegments` will return the model value for E1-only / E2-only segments (today the leaf σ, while
  χ² uses σ_E1 / σ_E2) — a behaviour change, with a test, in brief I. The MCMC snapshot writes the model value for
  such segments (accepted; test in brief I). Covariance seeds move into `Results` with S5.
- G2: a varied follower of a norm group keeps its own fixed parameter in the engine until S7 and is named
  `<record>.value[<data set>]`; the converter's name map must be injective (test); `varname[50]` becomes a string.
  `bad_key.azr` and `bad_sample.azr` stay recorded differences; `M_thick.azr` is a converter defect (brief G-3).
  Units other than the engine's are refused by `Build` (`E_NOT_BUILT`) until the units step. Data files are read
  relative to the working directory; `--describe` will print the resolved paths.

**Defects found and not yet fixed** (owner in brackets):
- A fixed energy shift changed between evaluations leaves a pooled working copy different from a fresh one
  (`natural_neon`: χ² 5589.035 pooled, 3342.510 fresh, already at 041edabf) [brief I, part 1].
- `tests/sum_lab_frame/sum_data.azr` and `tests/ratio_frame/ratio_data.azr` have no committed data [brief I].
- `EPoint::DopplerAppliedScope` sets a flag on the base point during an evaluation; `src/ModelPass.cpp` has a
  warn-once `static std::atomic<bool>` [Phase 6, mutable statics].
- `CNuc::ObservedOfFitValues` writes the observed slots of the levels during an evaluation [Phase 6 / S16].
- Pair labels are built in three places; `Build` has its own pair table (`PairRow`) beside the model data's
  `PairInfo` [brief G-3, part A].

**For the PI to see** (behaviour a user can meet; additions to the list of the entry of this morning):
- New refusals: a composite block with an operation type other than sum or ratio; β-delayed segments with target
  integration, a beam profile, screening or exit −1.
- The API value of E1-only / E2-only segments (above).
- Format 3: the follower parameter name; units other than the engine's refused for now; data-file paths relative to
  the working directory.

**Decisions of the PI, 2026-11-10** (they replace the list "waiting for the PI" of the earlier entries):
- **Minuit2 is not patched** for now. Cause B of thread-dependent fits (`fit_threads.md`: the non-atomic call counter
  of the vendored Minuit2) stays documented; the bounded fits of the poison test stay at one thread. Do not ask again
  unless a user meets it.
- **The two masses keep their values** (no result moves) and are named in one place, so that each is one line to
  change: `electronMassOfFermiFunction` and `protonMassOfM1Capture` in `include/Constants.h`
  (merged, 46031fa6; 179 × 3 identical).
- **A spin-orbit term may be added to the hybrid potentials**: step 3b (iii) of `cnuc_split.md` is wanted (channel
  functions keyed by (α, l, j)). The design note is written: `spin_orbit_potential.md` (brief N). Director's review:
  the recoupling example (p + ⁷Li, l = 1, J = 1) has the right trace and determinant for λ_j = −1, ½; `angular.py`
  re-run (PASS); e² = 1.4399764 against α ħc = 1.4399645 and the Whittaker S below threshold confirmed in the
  code. Director's decisions on its §6: no channel-spin operator (Q3); R_C a required field, 3.5 fm written by the
  converter (Q6); potential parameters are not fit parameters (Q7); matching radius and step are fields with logged
  defaults (Q8). **Decided by the PI, 2026-11-10 (yes to all four):** Q1 only pairs
  with one spinning particle at first, the others refused until a project needs the coupled treatment; Q2 FRESCO's
  convention (2.000 fm², 2 l·s); Q4 bound states and closed channels are distorted with the potential (moves
  `hybrid_potential`, removes a step of about 1 in S at threshold); Q5 e² = α ħc in the hybrid Coulomb potential
  (moves hybrid outputs by about 10⁻⁵). Both moves are classified with cause and size when they are made.
  The note's steps 1 (potential as a value,
  no singleton), 3 (format-3 record) and 4 (key with `two_j` = −1) move no output and need no decision.
- **M3 (b) is the aim** (`data_model.md` §5): lab data with a fitted shift are compared in the frame of the data —
  the model is converted to the lab, the measured values and errors stay as measured. Phase 5 does it; every fit
  it moves (shifted lab-differential data) is listed with its size. Until then (a) stays, bit-identical.
- **Confirmed by the PI:** the new refusals (a `prior_centre` row naming no segment line; an `<mcmc>` parameter row
  that names no parameter or is short; a composite block whose operation type is neither sum nor ratio; a
  β-delayed segment with target integration, a beam profile, screening or exit −1), the API returning the model
  value for E1-only / E2-only segments, and the format-3 words (`normalization` and `energy_scale` records, the
  follower parameter name `<record>.value[<data set>]`, `quadrature points <n>` required, `sigma_nuclear`,
  `entrance_pair` / `exit_pair`, the target-state grammar, units other than the engine's refused until the units
  step, data files relative to the working directory). Nothing waits for the PI now. No ASan runtime on this host.

**Running now** (Opus agents from fe181256, 3 threads each; briefs in `/data0/skowrons/AZURE2-upgrade/work/briefs/`):
- H, `feature/upgrade-4-expr` (`wt-expr`), `H_expression_pass.md`: S10 (the model pass evaluates the expression
  tree; `--print-predictions`) and S11 (the gradient through `backward`).
- I, `feature/upgrade-4-plan` (`wt-plan`), `I_plan.md`: the follow-ups of the model pass (fixed-shift history, API
  value, snapshot test, covariance seeds), then S5 (the calculation plan; one writer per slot).
- G3, `feature/upgrade-4-format3-reader` (`wt-fmt3`), `G3_engine_runs_format3.md`: the model data joined to the
  reader, the parameter-name stems, the engine runs a format-3 file, `--describe`, the converted-tree gate.

**Then.** S7 `Measurement` (after S5; removes the follower parameters); 3c (pyazr's model on records); 3d (the tree
converted, the old path deleted); units in `Build` with the new quantities — ratio to Rutherford/Mott (S14), S
factor (with the `nds.py` mass fix), yields per charge; S12, S13, S16; extrapolation grids (list, file, adaptive);
CNuc 6b and 7. Merge `origin/dev` only when the reconstruction is complete.

## 2026-11-10 (late morning) — format 3, step 3b-1, G2: records produce the engine's specs (`feature/upgrade-4-format3-reader`)

**Done.**
- `Schema::Build(validated records, pairs of <levels>, options) -> BuiltProject`
(`experiment/include/azure2/experiment/Build.h`,
  `src/Build.cpp`, `src/BuildEffects.cpp`): data sets, extrapolations, shared `normalization`/`energy_scale` (the
record's
  values written into every member; a scale's function only on its first active member, as
`EData::ResolveShiftFunctions`
  requires), effects (one `TargetEffectSpec` per distinct target/beam/detector/quadrature/fold_windows, ordered by the
first
  of its records), sample thicknesses, abundances (one row per isotope and compound nucleus its data sets reach), THM
  membership (`ThmMembershipSpec`).  Reactions resolved against the pairs (states to the digits written; §3.5 pair-key
rule).
- Build diagnostics with table codes and a test each (`experiment/tests/build_test.cpp`): `E_UNKNOWN_PAIR`,
`E_AMBIGUOUS_PAIR`
  (names and wording of the model-data branch), `E_PAIR_KEY`, `E_DATA_FILE`, `E_NOT_BUILT` (what the table takes and
this
  engine cannot set up yet: units other than the engine's, c.m. or file/adaptive grids, column layouts, epochs, ...),
  `E_THM_RELATIVISTIC`, `E_THM_COVARIANCE`, `E_FOLD_CONFLICT`; `E_FIELD_NOT_APPLICABLE` for a phase shift's (J, l, s).
- The reading path `ReadFormat3Text/Project` (`experiment/.../ProjectReader.h`, library `azure2_experiment_reader`,
links io/)
  replaces convert/'s `CheckFormat3Text`; `AZURE2 --validate` and `pyazr.validate` build too (data files relative to the
  working directory, as the engine opens them); the converter's self-check stays validation only.
- `ExtrapolationSpec` holds its grid as lists (`energies`, `angles`): the old adaptor fills them by the old
accumulating loops,
  `Build` by decision F10; `EData::MakePoints` iterates the lists (`ESegment::GetEStep/GetAStep` deleted).
- Spec-equivalence ctest `convert_spec_equivalence` (all 189 files: 186 convert, 183 equal, 3 recorded, 0 differ;
rules in its
  header; `AZURE2_SPEC_INJECT=1` shows it failing).
- Converter: targets with a screening U_e or a density/temperature fit form are never merged (each old line had its own
  parameter; tests/screening_6li/bare.azr had three merged); `PairLabel` writes a light particle's excitation; an
existing
  `<modelData>` section receives the converted `<scatteringLengths>` rows (one section); the `pair` value accepts
`n[0.3]+20Ne`.

**Verified** (build `work/b-fmt3`, 0 warnings, 3 threads): 179 × 3 against `tab-p3f` 689 / 662 / 689 identical (tag
`tab-fmt3b`); `ctest -LE slow` 149/152 then the 3 failures fixed (a test's code list, a note's text, the data-file
directory of `pyazr_validate`) and those plus `repo_consistency` rerun: pass; `tests/run_tests.sh` 135 passed, 0 failed;
consistency checks pass.

**Not done (next):** G2 item 4 (parameter-name stems; a norm-group follower has its own fixed parameter in the engine,
which
the converter's map sends to the same new name as the leader's: a decision); `<modelData>`/`<scatteringLengths>` are not
built (the director joins the model-data branch); G3 entirely (`Config::ReadConfigFile` on format 3, reading order and
`ExitPairKeysInFile`, the five §8.5 refusals, `--describe`, the converted-tree gate, docs).

## 2026-11-10 (day) — Phase 4 S1, S2 remainder and S4: one model pass (`feature/upgrade-4-model-pass`)

**Done** (`data_model.md` S1, S2, S4; brief F).
- **S2** (61dfec8c): a point an energy shift moves beyond the kinematic range keeps its load-time c.m. angle and factors
  (`EPoint::KeepLoadFrame` / `RestoreLoadFrame`, H4); one component-key counter per `EData`, the two
  `CreateComponentSegment` overloads one (H10); a beta-delayed spectrum is refused with exit -1 and with a target
  integration, beam profile or screening on its line (`tests/beta_delayed_alpha/check.sh`; no tree project uses them).
- **S1** (3ca0ff6b): `upgrade_history` (`tests/upgrade/history_test.cpp`) — fresh vs pooled copy, per-segment chi^2
  and model values bit for bit, on T2, natural_neon, shift_sqrt (sqrt, shiftfn), 18O_p_a_thm and a generated 7Li(p,n)
  point that a shift moves out of range; component keys of a second load.  On 041edabf it fails for H4 and H10 only.
- **S4** (18e5ce8b): `Results::model` and `Results::objective`; `include/ModelPass.h` (`EvaluateModel*`, the one
  combination `ComposeModelValue` in `src/ModelPass.cpp`, `TransformedValue` on a local vector, `DataResidual`);
  every former `CalculateTheoreticalCrossSection` loop (20 calls) and residual (13 expressions) goes through it;
  the objective's stored values (profiled norms, THM profiles, segment chi^2, totals) are in the buffer, the writer
  and the API read them there (`EData::NormUsed`); stale reads of S3 (b) fixed (MarkFailed, failed sample terms,
  automatic-tolerance E1/E2); check 28 (numbered 27 on the branch; renumbered at the merge); poison test widened
(natural_neon, transmission + yield under a spread,
  sample transmission).

**Verified** (3 threads, work/b-model): 179 x 3 identical to tab-p3f after S2, S4 and the final commit
(689/662/689); ctest -LE slow 150/150 (final tree); run_tests 136/136 (S4); calculations 1 vs 2 threads 12/12;
T1/T3 fits (MIGRAD and gradient, 150 calls, 1 thread) every output file identical to 041edabf; T2 bench calc ~1.0-1.1 s,
fit200 62-71 s on both, peak RSS 268 -> 271 MB (sizeof(EPoint) 816 -> 840 B).  Not run: 8 threads, TSan, slow ctests.

**Still written into the data set by an evaluation** (not moved here): the covariance seeds (M^-1 r)_i
(`ESegment::SetCovarianceSeed`, read by the gradient after `PrepareCovarianceSeeds`), and what the hand-over sets
(norms, shifts, thicknesses, grids).  **Open:** a fixed energy-shift parameter changed in the vector (pyazr can) moves
the
points, and a copy that evaluated a nonzero value and then 0 does not give a fresh copy's result (history test with
fixed shifts moved, on 041edabf; inferred cause: an inactive shift is not re-applied in
`FillEnergyShiftsFromParamsImpl`) — not fixed.  The API's regular `UpdateSegments` branch returns the leaf sigma of an
unassembled E1/E2-only segment, not its model value (unchanged; a decision).  Next: S5 (the plan).

## 2026-11-10 (morning) — Phase 4 S15: `<modelData>`, `ModelDataTerm`, ωγ first; scattering lengths as its rows (`feature/upgrade-4-model-data`)

**Done** (`data_model.md` §4, M10; `project_format3.md` §3.10; `objective.md` §5).
- `physics/ResonanceQuantities.h`: ω, ωγ = ωΓ_aΓ_b/Γ, Γ_b/Γ, Γ = Σ|Γ_c| with partial derivatives (test
  `physics_resonance_quantities`: T2 level 5 by hand, 12C+12C ω, derivatives against central differences); binding
  `_azure2.statistical_factor`, `resonance_strength`, `branching_ratio` on arrays.
- `experiment/ModelDatum.h`: `ModelQuantity`, `ModelQuantities()`, `BuildModelData(validated, LevelInfo[], PairInfo[])`
  (no engine: `<levels>` as plain tables); codes `E_UNKNOWN_LEVEL`, `E_LEVEL_JPI`, `E_UNKNOWN_PAIR`, `E_AMBIGUOUS_PAIR`,
  `E_NOT_A_CHANNEL`, `W_BROAD_LEVEL`, `W_NEARBY_LEVEL` (cases in `experiment_diagnostics`); table: strength and
  branching ratio refuse a `[total]` exit (schema JSON and `format3.rst` regenerated, `schema_version` 3.0 kept).
- Engine: `ALevel::IDs()` (the `<levels>` IDs of a level, from `CNuc::Fill`); `src/ModelData.cpp` reads
  `<scatteringLengths>` (unchanged grammar) and `<modelData>` (`io::ReadRecords` + `Schema::Validate` +
  `BuildModelData`), evaluates the level quantities from `CNuc::ObservedOfFitValues` (`src/CNucObserved.cpp`:
  TransformOut's two steps without messages, N_f^−½ put back), gradient and rows by central differences through
  `CNuc::SetFitLevelParameter`; scattering lengths exact as before. `ModelDataTerm` replaces `ScatteringLengthTerm` at
  T9; the six hand-written sites (`AZURECalc` gradient, LM rows, GSL rows; `AZUREAPI` gradient, rows; MCMC ×2) call
  `Objective::ModelDataGradient` / `AppendModelDataRows` / `ModelDataResiduals` / `ModelDataChi2`. Check 11 extended.
- Outputs: `modelData.out` and the run log table; `Total-Model-Data-Chi-Squared` in `chiSquared.out` (with the
  section); a project with model data and no segments runs. pyazr 2.16.0: `pyazr.model_data(model)` → `ModelDatum`
  records (from `objective_terms_rwa()["model_data"]`); `penalties()["scattering_length"]` is the whole term
  ("model data").
- Recipe "A model-only datum" in `CLAUDE.md`, check 27 (every datum quantity has its `ModelQuantities()` row and a
  `test_<quantity>`); docs `data_formats.rst`, `output_files.rst`, `chi_squared.rst`, `objective.md` §5.

**Verified** (2 threads, `work/b-mdata`, 0 warnings): 179 × 3 against `tab-p3f` 689 / 662 / 689 identical;
`scripts/check.sh --quick` (ctest −LE slow) 151/151; `tests/run_tests.sh` 136 passed, 0 failed (OMP 2, 2230 s; with
`model_data`). Cost of the term on T2 (model data alone, 1 thread): about 70 µs per evaluation (Brune, Park), 360 µs
(formal: the Brune iteration), against about 77 ms for a T2 evaluation; T2 fit200 (`--use-lm --max-calls=200`, 1
thread, table on) 16.3–17.0 s without, 16.5 s with one strength. Not run: 1 vs 8 threads (budget), sphinx, `ctest -L
slow`.

**Open.** The converter writes `pair` labels without the light particle's excitation (`Builder::PairLabel`), which
`BuildModelData` then matches to the ground-state pair or refuses; converting a format-2 file that already has a
`<modelData>` section writes a second one (the reader takes the first only). Full report: `work/mdata/REPORT.md`.

## 2026-11-10 (morning) — format 3, step 3b-1, G1: the engine's set-up consumes specs (S9) (`feature/upgrade-4-format3-reader`)

**Done.**
- `experiment/include/azure2/experiment/SetupSpecs.h` (plain, standard library only): `DatasetSpec`,
  `ExtrapolationSpec` (sharing `ReactionSpec`: pair keys, an `Observable` instead of the isDiff code, THM, phase-shift
  J/l/s, SUM/RATIO with `ComponentSpec` values, sample, total, transmission/yield), `TargetEffectSpec` (the data sets
  it applies to as keys), `AbundanceRowSpec`/`AbundancesSpec`, and `SpecList<Spec>` (items in order, each with the
  reader's warning/note, and the reason the reader stopped). The prior centres are fields of `DatasetSpec` (F14).
- The engine takes specs: `ESegment(const DatasetSpec&)`, `ESegment(const ExtrapolationSpec&)` (one `SetReaction` for
  the flags), `TargetEffect(const TargetEffectSpec&, const Config&)`, `EData::Fill(cfg, cnuc, const DataSpecs&)`,
  `EData::MakePoints(cfg, cnuc, const ExtrapolationSpecs&)`, `EData::SetUpTargetEffects` (was `ReadTargetEffectsFile`),
  `EData::ApplyPriorCentres` (was `ReadPriorCentres`), `CNuc::SetAbundances`. No engine constructor takes an old line;
  the component list is a vector of values (the "Entrance: X, Exit: Y" string round trip is gone).
- The old-line adaptor, temporary: `convert/`'s `SpecsFromOldLines.h/.cpp` (library `azure2_convert_specs`; the old
  reading loops with their messages, words and order), called by `src/OldLineSetup.cpp` (the old-signature
  `EData::Fill/MakePoints(cfg, cnuc)` and `CNuc::ReadAbundances`). Step 3d deletes both.
- Converter: a `prior_centre` row the engine refuses (a name that matches no `<segmentsData>` line, e.g.
  `segment_2_norm` with one line or `segment_01_norm`) is refused with the engine's reason (it was dropped silently);
  a composite block with an operation type other than 0 / 1 is refused by both (new: the old engine ran it, as a
  RATIO in `Fill` but not a ratio in the c.m. conversion).
- Test `old_line_specs` (every isDiff code of data and test lines against the flags of the old constructors, written
  by hand; components bit for bit; the reading refusals with their words; prior-centre names equal
  `EData::SegmentParameterName`; abundances; `<targetInt>` keys); shown failing with an injected flag error. Three
  `convert_constructs` cases.

**Verified** (build `work/b-fmt3`, 3 threads, 0 warnings): 179 × 3 against `tab-p3f` 689 / 662 / 689 identical
(run.log included; tag `tab-fmt3`); `ctest -LE slow` 150/150 (OMP 2, -j 1); `tests/run_tests.sh` 135 passed, 0 failed
(OMP 3, 1832 s); consistency checks pass; size ratchet lowered
(`EData.cpp` 5883 → 5590, `ESegment.cpp` 1816 → 1743, `CNuc.cpp` 1954 → 1931; `pybind_api.cpp` unchanged at 1299).

**Not done in G1 (left for G2/G3):** a THM spec (the `<thm>` block is still read by `ThmExperiment`'s parser into
`ThmExperiment`, which is already a plain struct; `ThmLines.h` stays included); `<scatteringLengths>` (brief E);
shared norm / energy scale as records of their own (the spec carries the group labels, the engine still takes a
group's values from its first member); parameter-name stems in the spec (G2 item 4); `src/AZURE2.cpp`'s pair check
still reads `SegLine`/`ExtrapLine` (G3 reading order).

**Behaviour that changed only for edge cases** (no test output moved): two `prior_centre` rows for the same parameter of
an inactive line warn once (last row) instead of twice; a `<targetInt>` line shorter than its convolution
coefficient count no longer reads past the vector (missing coefficients are 0).

## 2026-11-10 (director, morning) — Phase 4: the second wave merged (S3 b, parameter kinds, CNuc split, converter)

**Merged into `feature/upgrade-4-data`** (`--no-ff`, in this order; each branch verified by the director's own runs
before its merge — build with 0 warnings, 179 × 3 identical to the Phase 3 head, `ctest -LE slow`,
`tests/run_tests.sh`): `feature/upgrade-4-results-b` (S3 b), `feature/upgrade-cnuc-split-3` (078dc233),
`feature/upgrade-4-param-kinds` (06013f09), `feature/upgrade-4-convert` (041edabf). Their entries are below.
**The merged head as a whole** (041edabf, `work/b-p4r`, `work/verify-p4r.log`): build 0 warnings; 179 × 3 identical to
the Phase 3 head
(689 + 662 + 689 files); `ctest -LE slow` 149/149; `tests/run_tests.sh` 135/135; calculations at 1 and 8 threads
12/12 identical; `upgrade_poison_results`, `upgrade_results_owner` and the two converter tests repeated six
times, all passed; `tools/convert_tree.sh`: 186 of 189 files convert and validate with zero errors, the 3 refused are
refused by the old engine too, 5 files the old engine refuses are converted (the list of `project_format3.md` §8.5).

**What the director changed at the merges** (each in the merge commit's message):
- Check numbers: the size ratchet is check **26** (25 is "no `Results` member in `EPoint.h`/`EData.h`").
- `ALevel::PadChannels` pushes the per-channel vectors once; `ALevel::AddGamma` deleted (8b9ab7cb).
- `<abundances>` has one parser, `convert/`'s `ReadAbundanceLine` (the converter needs it line by line);
  `io/AbundanceLines.{h,cpp}` and the transfer test `io_abundance_lines` are deleted. `convert/` includes no engine
  header (`NucLine.h` is gone; `io::LevelRow`); check 9c no longer allows one.
- `<parameterSettings>` `prior_centre` rows: the engine keeps the reader of the parameter-kinds branch
  (`EData::ReadPriorCentres`: the whole name against `SegmentParameterName`); `convert/`'s `PriorCentreRow.h` serves
  the converter only. **Difference to close in 3b:** the converter drops a row whose N names no `<segmentsData>` line;
  the engine refuses it.
- `tests/thm_distortion/check.sh` drops the `prior_centre` rows of the segment lines its reduced project removes
  (dda649f7): since the parameter-kinds branch a row naming no line is an error. **User-visible:** a hand-reduced
  copy of a project that keeps such rows is refused with a message naming the accepted names; the shipped examples
  are not affected.
- `tests/upgrade/results_owner_test.cpp` called `FillMnParams` with the old argument after the parameter-kinds merge:
  **06013f09 does not compile that test** (the engine builds); repaired in 041edabf. The cause was merging on the
  branch's own verification without building the merged tree first; from now on the merged tree is built before the
  merge is committed.
- `tests/consistency/file_sizes.tsv` was re-baselined at each merge (the record was taken before the branches grew
  files in parallel). Raised in total against the first record: `pybind_api.cpp` 1299, `azrfile.py` 3756,
  `CNuc.cpp` 1954 (from 1942: `PadChannels`, the abundance loop); lowered: `AZURE2.cpp` 1721, `ThmExperiment.cpp`
  1026; `TargetEffect.cpp` (932) left the list. 18 files remain above 1,000 lines.
- pyazr 2.15.0 (`validate`, `convert_azr(..., to=3)`).

**Accepted from the branches' reports, for the PI to see** (behaviour a user can meet):
- F1 `prior_centre` rows as above. F2: a `<mcmc>` `<parameters>` row that names no parameter, or has fewer than five
  tokens, stops the run (before: ignored). F3: limits of nuisance parameters are converted eV → reduced only for
  widths.
- Format decisions that fix what users write (`project_format3.md` §10 and the step 3a entries): record kind
  `normalization`; `quadrature points <n>` required on a folded data set; quantity `sigma_nuclear`; fields
  `entrance_pair` / `exit_pair`; the target state `16O[2.143](n,n)16O[2.143]`; the temporary flag `--to 3`.

**Waiting for the PI** (unchanged): the two atomic increments in the vendored Minuit2 (cause B of thread-dependent
fits; `external/` may not be touched without asking; bounded fits of the poison test run at one thread until then);
M3 (b); the electron mass in `EPoint.cpp` and the proton mass in `ECIntegral.cpp`; whether a spin-orbit hybrid
potential is wanted (CNuc 3b iii). The host has no ASan runtime for this compiler (the sanitizer jobs run in CI
only).

**Running now** (Opus agents from 041edabf; briefs in `/data0/skowrons/AZURE2-upgrade/work/briefs/`):
- F, `feature/upgrade-4-model-pass` (`wt-model`): S1 history test, S2 remainder, S4 one model pass and one residual.
- E, `feature/upgrade-4-model-data` (`wt-mdata`): S15, `<modelData>` with `datum` records, `ModelDataTerm`, ωγ; the
  scattering lengths as rows of the same term.
- G, `feature/upgrade-4-format3-reader` (`wt-fmt3`): 3b-1 with S9 — set-up from specs, `Schema::Build`, the engine
  runs format 3, `--describe`.

**Then.** 3c (pyazr's model on records), 3d (the tree converted, the old path and `convert/`'s engine use deleted,
`--convert-params`); S5, S7, S10–S13 (the tree pass replaces the composition and the hand-written adjoint), S14, S16;
ratio to Rutherford/Mott, S factor (with the `nds.py` mass fix), yields per charge; extrapolation grids; CNuc 6b and
7. Merge `origin/dev` only when the reconstruction is complete.

## 2026-11-10 (night, later) — step 3a follow-up: the director's D1–D5 and review R1–R3 (`feature/upgrade-4-convert`)

**Done.** D1: quantity `sigma_nuclear` (the angle-integrated nuclear elastic integral of a charged pair; warning
`W_NOT_MEASURABLE` in a data set); the converter writes it with a `# CONVERTED:` comment. D2: fields `entrance_pair` /
`exit_pair` (pair keys of `<levels>`, only where a label names two pairs; not on a sum or ratio) — the THM-column rule
was not taken (`project_format3.md` §3.5 says why). D3: a target state in the reaction grammar,
`16O[2.143](n,n)16O[2.143]`. D4: the converter converts what is written; `tools/convert_tree.sh` runs the old engine
on every original and lists the files it refuses that the converter converts (§8.5 records them). R1: the project-file
commands (`--convert`, `--validate`, `--schema`) moved to `src/ProjectCommands.cpp`; `cli_help` reads its literals
too. R2: every `TargetEffect` member has an in-class initialiser (test `target_effect_line`).

**Checked on the work tree before the last commits** (2 threads, `work/b-conv`, 0 warnings): `tools/convert_tree.sh`
186 of 189 convert and validate with no question, the 3 refused are refused by the old engine too, 5 files the old
engine refuses are converted (§8.5); `convert_constructs` 155/155, `convert_fixtures` 86/86, `experiment_*` pass,
`pyazr_validate` 35/35, `cli_help` PASS, `repo_consistency` passes; `target_effect_line` passes, valgrind 0 errors
(with the old header: uninitialised reads). The 179 × 3, `ctest -LE slow` and `run_tests.sh` of the final tree are the
director's.

**For step 3b (recorded here, not done):**
- `--convert-params` for `param.sav` / `param.par` files (the old -> new parameter map is in the converter's report);
- `ExtrapLine` never reads a UPOS block on test lines (`isUPOS_` stays 0): decide whether test lines can carry one;
- `convert/include/azure2/convert/ThmLines.h` duplicates two small generic helpers (split at a separator, whole
  integer) that `src/ThmExperiment.cpp` keeps for its other readers; one copy should remain when `<thm>` moves;
- `CheckFormat3Text` (read the format-3 sections with `io::ReadRecords`, validate with the table) belongs to the
  format-3 reader: move it from `convert/` there, and the CLI's `--validate` and pyazr's `validate` with it;
- the pair-key rule of §3.5 (`entrance_pair`/`exit_pair` accepted only where the label is ambiguous, required there)
  and the target state are checked against `<levels>` by the reader (the validator cannot see `<levels>`);
- every file of §8.5's list (the old engine refuses it, the converter converts it) must be refused by the format-3
  engine.

## 2026-11-10 (night) — format 3, step 3a: named table rows, `--validate`, the converter old → format 3 (`feature/upgrade-4-convert`)

**Done** (`project_format3.md` §§6.2, 8, corrections marked "step 3a").
- **Named rows.** `Quantities.cpp`: one function per quantity filling `Quantity` by member; `Schema.cpp`: one function
  per record kind, field helpers that say where a field applies and whether it is required (no positional booleans),
  structured values from `Keyword`/`Unkeyed` alternatives; `Units.cpp` rows without the empty alias braces. Schema JSON
  byte-identical (cmp). Recipes in `CLAUDE.md` and `experiment/README.md` show a complete row.
- **Step-2 loose ends.** Column role `E_x` removed (no quantity takes it; `frame E_x` is a field value); §3.6 example
  with step 2's corrections; §9.3 at 13.2652 MeV. Five constructs the old engine reads and the table refused, found by
  converting the test tree, now taken (each with a case that fails when reverted): `norm profiled` on `thm global`,
  legendre with apparatus, window columns on `dsigma_domega`, fixed norm on `beta_delayed_spectrum`, charged-pair
  `scattering_length`. `schema_version` stays 3.0 (nothing accepted before changes meaning); JSON and
  `format3.rst` regenerated.
- **`AZURE2 --validate FILE [--json]`**, `pyazr.validate(path)`: every diagnostic, the number of checks left to 3b
  (the quantities' "build" refusals) as a note, exit 0 / 2 / 3 (format 1/2 refused with §8.4's message).
  `io::FormatOf` knows format 3; this engine refuses a format-3 project with a message (it runs format 2 until 3b).
- **`convert/`**: the old line parsers (`SegLine.h`, `ExtrapLine.h`, `ComponentList.h` moved; `TargetIntLine.h`,
  `AbundanceLine.h`, `ScatteringLengthLine.h`, `PriorCentreRow.h`, `ThmLines.h` taken out of the engine's routines
  unchanged) and the converter `azure2_convert` (`ConvertToFormat3`, `CheckFormat3Text`), used by
  `AZURE2 --convert OLD NEW --to 3 [--json-report R] [--answers A] [--exact-grids] [--strict] [--allow-questions]` and
  `pyazr.convert_azr(old, new, to=3, ...)` (pyazr 2.14.0). `tools/convert_tree.sh` converts a copy of the test tree and
  prints the table of §8.1.

**Temporary dependencies of the engine on `convert/` (3b removes exactly these):** `AZURE2Core` links
`azure2_convert_lines`; `src/EData.cpp`, `src/ESegment.cpp`, `src/AZURE2.cpp`, `api/src/pybind_api.cpp` include
`SegLine.h`/`ExtrapLine.h`; `src/TargetEffect.cpp` `TargetIntLine.h`; `src/CNuc.cpp` `AbundanceLine.h`;
`src/ScatteringLength.cpp` `ScatteringLengthLine.h`; `src/EData.cpp` `PriorCentreRow.h`; `src/ThmExperiment.cpp`
`ThmLines.h`. Not temporary in the same sense: the CLI (`AZURE2`) and `_azure2` link `azure2_convert` for `--convert`
and `--validate`; `CheckFormat3Text` (read the format-3 sections with `io::ReadRecords`, validate with the table) is the
driver the 3b reader should own — move it there. `convert/` reads `<levels>` with the engine's `NucLine.h` (check 9c).

**Verified** (2 threads; build `work/b-conv`, 0 warnings): 179 × 3, menu 1, the engine code of these commits against
`tab-p3f` (tag `tab-conv3`): 689 / 662 / 689 files identical (Brune, formal, Park), run.log included; `ctest -LE slow`
145/145 (142 + `convert_constructs`, `convert_fixtures`, `pyazr_validate`); `tests/run_tests.sh` 134 passed, 0 failed
(OMP 2, 2261 s);
`tools/convert_tree.sh`: 168 of 189 files convert and validate with zero errors, the 3 the old engine refuses are
refused (`angular_coefficients/models/bad_norm`, `bad_test`, `electron_screening/bad`), 18 stop at questions of the
three kinds of §8.3. Not run: calculations at 8 threads (budget), sphinx.

**Open (decisions for the director):** a quantity for the angle-integrated nuclear elastic integral of a charged pair
(16 test files); pair names or another way to tell two identical pairs apart (`c12c12_tumino2018`); an excited target in
the reaction grammar (`thm_coherent_background/toy`). `--convert-params` (param.sav) is not written; the parameter map
is in the report.

## 2026-11-10 — parameter kinds: what a parameter is comes from its record, not its name (`feature/upgrade-4-param-kinds`)

Record and recipe: `parameter_kinds.md` (format 3 check 23, decision F3).
- **Done.** `include/ParameterKind.h` (plain, standard library only): `ParamKind` (codes 0–11 = the API's type codes;
  moved from `AZUREGrad.h`, the MCMC's `PARAM_*` and pyazr's `_KINDS` deleted), `ParamKindTable()` (key, unit, owner
  fields, `ownPrior`, `segmentField`, `varyFromProject`, doc) and `ParameterRecord` (kind, ordinal, owner fields).
  `AZUREParams::Add(name, value, step, record)` is the only way a parameter enters the vector (`Record(i)` throws if a
  parameter lacks one); `CNuc::FillMnParams` / `EData::FillMnParams` take `AZUREParams&`. Every name test is gone:
  `AZUREParams` (param files), `ParameterLimitsManager` (GUI numbering, setting names and index fallback, nuisance
  exclusion, width conversion; six dead functions deleted, `FindParameterIndex` that refilled the vector per call
  replaced by the index), `AZURECalc`, `AZUREMain`, `AZURECalcMCMC`, the CLI MCMC (`PrettyParameterName`,
  `FindMcmcPriorParameter`), `EData::ReadPriorCentres` (whole names from `EData::SegmentParameterName`), the API's
  norm / shift index lists, pyazr (`load_sav`, `save_fit`'s cbkg, an abundance's sample). `ParameterLabel.h`'s
  `IsNormName` & co. and the dead `ApplyThmCoherentValues` deleted. Binding: `Session.parameter_records()`,
  `_azure2.parameter_kinds()`; pyazr 2.14.0 (`parameters.parameter_kinds()`, `Parameter.thm_*`, `level_number`).
  Check 23 in `repo_consistency_test.py` (shown failing on the start head: 70 problems). Test
  `pyazr_parameter_kinds`.
- **Second commit:** one element table: `ParameterLabel.cpp` (stopped at Kr) and `pyazr.nuclei.SYMBOLS` take
  `physics/Elements.h` (binding `_azure2.element_symbols()`); check 10b.
- **Review follow-ups (director, 2026-11-10).** F1: a `prior_centre` row is applied to an active segment, reported
  and ignored for a `<segmentsData>` line that is not active, and an error for any other name (names built with
  `SegmentParameterName` over all line keys; the old error for a malformed name, now also for a key beyond the last
  line or a typo). F2: an `<mcmc>` prior row that names no parameter stops the run with a message listing the
  accepted forms. F3: `ParameterLimitsManager` converts nuisance settings for widths only; a nuisance row on a
  screening or cbkg parameter wrote past the end of `CNuc::GetTransformParams()` (an `_GLIBCXX_ASSERTIONS` build
  aborts; the cbkg case also crashed the Release build of the start head). Tests: `tests/prior_centre` (unknown and
  typo rows refused, inactive line warned), `tests/mcmc_classification` (pp_short, pp_badrow), `tests/nuisance_kinds`.
- **Known defects, recorded for later.** (1) The CLI MCMC table labels an energy shift "Segment k Energy Shift
  (keV)" (`PrettyParameterName`, `src/AZURE2.cpp`), but the shift is in MeV; the string is also the GUI's matching
  key, so fix it when that key goes (Phase 9, web UI). (2) pyazr parses the `prior_centre` rows with the same regex
  twice (`pyazr/datasets.py` `_PRIOR_CENTRE`, `pyazr/azrfile.py` `_PRIOR_ROW`); one reader, or the engine's, in the
  format-3 pyazr step.
- **Open.** The layout is still walked a second time by position in `BuildParamIndexMap`, `GetParameterInfo` (now
  checked against the records), `AZURELabel::Parameter`; building them from the records is the next step. pyazr's
  document layer still parses the `prior_centre` rows by regex (no session; listed in check 23).

## 2026-11-10 — `CNuc` split continued, and a file-size ratchet (branch `feature/upgrade-cnuc-split-3`)

**Done** (record: `cnuc_split.md` §2.1 inventory, §6 "Steps (d), (c)", "Remaining steps").
- Check 26 (`tests/consistency/file_sizes.tsv`; numbered 25 on the branch): every hand-written source above 1000 lines
is listed and may not
  grow; an unlisted file may not pass 1000 lines; `--tighten-sizes` lowers the record (never raises). 19 files.
- `CNuc.cpp` 2907 → 1942 lines, without touching the parameter vector: output and check files (`CNucReport.cpp`),
  load-time checks (`CNucChecks.cpp`), transformation messages (`CNucTransformMessages.cpp`), all through const access;
  the `<levels>` and `<abundances>` grammars in `io/` (`LevelRows.h`, `AbundanceLines.h`; `NucLine.h` deleted);
  the unreachable branches of `ParseExternalCapture` and the write-only max-l members deleted.

**Verified** (build `work/b-cnuc3`, GCC 11.5, Release, 0 warnings, 2 threads): after each engine commit 179 × 3
against `tab-p3f` 689 / 662 / 689 identical and ctest −LE slow all passed (142, then 144); check files on and
malformed inputs compared with the base build; `io_level_rows` (10554 rows) and `io_abundance_lines`; 1 vs 2 threads
identical (T1, T2, T3, 13N × 3 bases); `run_tests.sh` 134 passed, 0 failed (OMP 2, final head). One-thread
calculation, mean of three,
base → branch: 13N 7.92 → 7.92 s, T2 0.28 → 0.29 s (Coulomb table on; machine shared).

**Open.** Step 6b and 3b (iii) need a design first; `ExitPairKeysInFile` needs a decision; step 7 waits for the
parameter-vector work.  Found, not fixed: `ALevel::PadChannels` appends `ec_gammas_`, `observed_values_` and
`observed_kinds_` twice per padded channel (a merge leftover; the extra entries are never read).

## 2026-11-10 — Phase 4 S3 (b): the results buffer belongs to the caller (`feature/upgrade-4-results-b`)

**Done** (`data_model.md` S3, decision M2 (c)). b1 (79c89e50): every accessor and engine function that reads or writes a
result takes the buffer (`point->GetFitCrossSection(results)`, `Calculate(cnuc, cfg, results)`,
`CalculateTheoreticalCrossSection(i, cnuc, cfg, data, results)`, `ObjectiveContext::results`, the gradient engine, the
writer, the THM profile, the covariance chi^2, the band, the MCMC, the API); the form is stated at the top of
`include/Results.h`.  b2 (3897fec9): `EPoint` holds only its offset and capacities; `EData` holds the layout
(`LayoutResults`, `ResultsSize`) and no buffer; `FillEnergyShiftsFromParams(p, results, ...)` ends with
`BeginEvaluation(results)` (re-layout carrying values in that buffer, sizing, poisoning).  Owners: `AZUREMain` (passes
its
buffer to `AZURECalc`), the pools of `AZURECalc` and `AZURECalcMCMC` (`WorkingData` = copy + buffer,
`include/WorkingData.h`), the iteration and MCMC snapshots, the trace / pool-check fresh copies, probe copies and
reaction-rate points (`EPoint::UsePrivateRecord`), the Python session (`AZUREAPI::results()`).  The gradient and the
MCMC
took ownership without S4: nothing is left for a "b3".  K1 stays the one listed exception (not fixed: owning the buffer
does not decide what E1/E2 of a differential line should be).

**Tests:** `upgrade_results_owner` (T2 copies, one copied from another, at two parameter vectors whose layouts differ,
sequential and concurrent, bit for bit against fresh copies; fails with a shared buffer injected into `EPoint::Rec`);
`upgrade_poison_results` now also runs a reaction rate; consistency check 25 (no `Results` member in `EPoint.h` /
`EData.h`; fails on an injected member).

**Verified** (b-resb2, this host, 3 threads at most): 179 × 3 identical to tab-p3f (689 / 662 / 689); ctest −LE slow
143/143 (b1: 142/142); `run_tests.sh` 134 passed 0 failed; calculations at 1 and 3 threads 12/12 identical (8 threads
not run: thread budget); T1 fit (menu 2, T1_start) and T3 fit (`--use-gradient`, T3_start) at 1 thread: every output
file identical to a16f2a69.  Side by side with a16f2a69 at 1 thread, three interleaved runs (host load ~8): T2 calc
1.04–1.13 s vs 1.06–1.11 s; fit200 66.4–68.4 s vs 66.0–69.8 s; fit peak RSS 272–276 MB vs 268–271 MB.
`sizeof(EPoint)` 824 → 816 bytes.

**Behaviour that moved only for reads that were already stale** (no output of the tests moves): a clone's first
evaluation no longer starts from the values of the data set it was copied from (its buffer starts at 0); the API
session's buffer is zeroed after `Initialize` (before, a re-initialised data set kept the values of points whose
records did not move).  Both matter only to a read before write.

**Stale reads still open (S4)**, at b2's lines: (1) `EPoint::Calculate`, automatic tolerance, `allPass` branch
(`src/EPoint.cpp` 2095): σ_E1/σ_E2 of an E1/E2 line with a sub-point fold kept from the copy's previous evaluation (now
from its own buffer; a pooled copy's first evaluation reads 0 instead of the master's values); (2)
`ESegment::TransformedValue` / `CalculateTheoreticalCrossSection`: a failed term returns NaN before `SetNaturalTerms`
(`src/ESegment.cpp` 1579 vs 1597, 1677 vs 1706/1737); (3) `MarkFailed` (`src/EPoint.cpp` 1927) leaves NaN σ_E1/σ_E2;
(4) `AMatrixFunc::PointTBar` reads the forward value (`src/AMatrixFunc.cpp` 726) — written by the forward pass of the
same evaluation, so not stale as long as the gradient follows an evaluation into the same buffer.

**Found, not fixed:** `upgrade_poison_results` can fail by chance at 2 threads: `fitStatus.out` of the T1 fit holds
Minuit's call count, which is not atomic (cause B of the fit-thread entry below); a16f2a69 gave 58 instead of 59 calls
in
1 of 12 runs, b2 in 0 of 12.  Comparing that file without its `calls` line, or the Minuit patch, would remove it.

## 2026-11-10 (director) — Phase 4: four branches merged; decisions; what runs next

**Merged into `feature/upgrade-4-data`, and that into `feature/upgrade`** (verified by the director's own runs,
`/data0/skowrons/AZURE2-upgrade/work/verify-p4i.log` and `verify-p4j.log`: 179 × 3 identical to the Phase 3 head —
689 + 662 + 689 files; `ctest -LE slow` 142/142; `tests/run_tests.sh` 134/134; calculations at 1 and 8 threads, 12/12
identical):
- `fix/upgrade-fit-threads` — cause A of the thread-dependent fits (entry below; `fit_threads.md`).
- `feature/upgrade-4-results` — S3 (a), the `Results` buffer behind the old accessors (entry below).
- `feature/upgrade-4-experiment` — **S8**: the expression nodes `Leaf`, `Sum`, `Scale`, `Weighted`, `Ratio`,
  `Transform` with `value` and `backward` (`experiment/README.md` §1; not wired into the engine). **S6**: the
  lab → c.m. functions in `experiment/` (`Frames.h`), `EPoint`'s conversions call them, one conversion sequence
  `ESegment::ConvertPointToCM` and one E_c.m./E_lab factor `ESegment::CMEnergyFactor` (README §2).
- `feature/upgrade-4-format` — **format 3, step 2** (`project_format3.md` §10): the generic record reader
  `io::ReadRecords`, the table (8 record kinds, 18 quantities, 43 units, 49 diagnostic codes) and `Schema::Validate`
  in `experiment/` (README §3), `AZURE2 --schema`, `pyazr.schema()`, the generated page
  `docs/source/reference/format3.rst`, checks 18, 21, 24. Not wired: the engine still reads format 2.

**Decisions (director, with the PI's standing instructions).**
- Data model M1–M13 (`data_model.md` §5) and format F1–F20 (`project_format3.md` §10) as recorded there. PI,
  2026-11-09: no backward compatibility in the engine, a converter instead; the format must read naturally in the
  GUI and in Python.
- F11: `experiment/` may include `azure2/io/Record.h` (header-only target `azure2_io_record`) and nothing else of
  `io/`; check 9b enforces it.
- Step 2's corrections to the note are accepted: the shared-normalization record kind is `normalization` (record
  kinds and field keys are disjoint within a section); keys match `[A-Za-z][A-Za-z0-9_]*`; `profiled` may stand
  alone; a formula states the unit of its result; a folded data set requires `quadrature points <n>`.
- All knowledge of the old formats goes into one directory, `convert/` (target `azure2_convert`), used only by
  `AZURE2 --convert` and `pyazr.convert_azr`; when the format-3 reader is wired the old parsers leave the engine and
  the rule "Old inputs keep working" of `CLAUDE.md` is replaced by "the engine reads one format; the converter reads
  the old ones".
- State and extension points (PI approved 2026-11-10; `ARCHITECTURE.md` §9, `CLAUDE.md`): the only singletons are
  constant tables behind functions returning `const&`; what an evaluation writes is owned by its caller; an
  extension is one explicit table row and a test; no self-registering classes, no abstract factories. The ~50
  mutable function-level statics, `NuclearPotentialManager` and the `thread_local` caches of `ECIntegral.cpp` go in
  Phase 6, with a consistency check.
- The quantity rows of `experiment/src/Quantities.cpp` are positional (twenty values, bare booleans): not accepted
  as the final form; step 3a rewrites them with named members, schema JSON unchanged.
- S6 changed one behaviour: a β-delayed entrance with no exit pair now gets the beam conversion where the old code
  dereferenced a null pair. Not accepted as it is: it must be a refusal with a message (with the S2 remainder).

**Waiting for the PI.**
- Cause B of the thread-dependent fits: Minuit2's call counter is not atomic (`external/minuit2/src/MnUserFcn.cxx:21`,
  `MnFcn.cxx:26`; `mutable int fNumCall`, `inc/Minuit2/MnFcn.h:53`). Two atomic increments fix it (shown on a scratch
  build). `CLAUDE.md` requires asking before `external/` is touched.
- M3 (b) of `data_model.md`; the electron mass in `EPoint.cpp` and the proton mass in `ECIntegral.cpp`.

**Running now (Opus agents; briefs in `/data0/skowrons/AZURE2-upgrade/work/briefs/`).**
- A, `feature/upgrade-4-results-b`: S3 (b) — the buffer passed explicitly, owned by the caller; `EPoint` and `EData`
  hold none.
- B, `feature/upgrade-4-param-kinds`: parameter kinds as an attribute of the parameter vector; no code reads a
  kind from a name (check 23); prerequisite of user-chosen names.
- C, `feature/upgrade-cnuc-split-3`: file-size ratchet; `CNuc` inventory and the moves that do not touch the
  parameter vector (output, reading, N_f^½, 6b, 3b (iii)).
- D, `feature/upgrade-4-convert`: named table rows, `AZURE2 --validate`, the converter old → format 3 in `convert/`
  (`--convert … --to 3` until the reader lands).

**Then, in this order.** Format 3b: the reader (validated records → what the engine consumes), `--describe`, the two
gate runs of F10, the test tree converted, the old parsers out of the engine. S2 remainder (load-time angle and factor
beyond the kinematic range, one component-key counter per `EData`, the β refusal) and S1's history test on more
projects. S4 (one model pass) with the stale reads of the S3 entry, each with a test that reaches the path first; then
S5, S7, S9–S13. S15 (`<modelData>`, ωγ). New quantities in the order of `data_kinds.md` §5.3. `CNuc` step 7 after B.

**Open findings, not yet assigned.**
- `param.fit` (written at evaluation 100, 200, …) holds a schedule-dependent vector under the parallel gradient.
- `EPoint::GetECAmplitudeWithShift` may keep the level energy of a copy's first shifted evaluation (inferred;
  no test reaches it).
- From S8/S6, for S9–S11: a `Transform` holds its own terms (not `Transform(n; Sum)`) to keep the bits; SUM components
  and the RATIO denominator use the total σ whatever the E1/E2 filter of the base (forward and gradient disagree for
  a SUM with an E1/E2 base); a component's fixed angle is not always c.m.; two Jacobian fallbacks differ per call
  site; `backward` differs from `AZUREGrad` by ≤ 2 ulp for a weighted and scaled SUM term; with relativistic=1 the
  angle "kept at the last valid value" is the non-relativistic angle, and the sub-point shift factor is
  non-relativistic; `angleKinFactor_` is dead state.
- pyazr `nds.py` converts S factors with mass numbers (σ +2.9 % for ³He+α at 93 keV, +6.7 % for p+¹²C at 100 keV);
  fixed when the S factor becomes an engine quantity.
- `src/ParameterLabel.cpp` and `pyazr/nuclei.py` have their own element lists beside `physics/…/Elements.h` (in B's
  brief as a second commit).
- Idea of the PI, outside AZURE2: an experiment logbook (yields, target scans) as a separate small Python package on
  `experiment/` and the stopping module; not before Phase 5. AZURE2 stays an R-matrix code.

## 2026-11-09 (night) — fits that depend on the thread count: two causes (branch `fix/upgrade-fit-threads`)

Record: `fit_threads.md`. Diagnostic `AZURE_TRACE_FCN=1|fresh` (every χ² with its parameter vector in hex floats;
`fresh` also on a fresh model copy).
- **Cause A, fixed:** `EData::ReplaceSubPointGrids` kept a sub-point grid whose nodes lay within 10⁻¹² of the new
  ones, so a pooled copy's grid, and its χ² in the last bits, depended on the vectors it had evaluated before. Shown:
  pooled vs fresh copy differ at 60 of 79 T2 evaluations (none with `--no-grid-replacement` or exact comparison);
  1 vs 4 threads, 34 of 59 common vectors with different χ². Fix: exact node comparison. After it, T2 at 1, 2, 4
  threads evaluates the same 216 vectors with the same bits; 8 threads repeatable (3 runs). The 1-thread T2 fit moves
  in its last digits (prior-covariance term 0.306025731622686 → 0.306025733674053; total 691.608). Cost within noise
  (fit200 63.4 → 63.8 s at 1 thread). Test `tests/fit_threads/check.sh` (ctest `fit_threads`), fails with the
  tolerance restored.
- **Cause B, not fixed (needs `external/`):** Minuit2's call counter (`MnUserFcn.cxx:21`) is not atomic; at 4–8
  threads it loses 13–28 counts in 200, so a `--max-calls` fit runs one more MIGRAD iteration, timing-dependent. This
  is the whole 691.608 / 595.877 difference. An atomic increment (scratch build, not committed) makes T2 at 8 threads
  and T1 at 4 threads end where 1 thread ends. Decision for the maintainers: patch the vendored Minuit2 (two lines).
- **Open:** `param.fit` (snapshot at evaluation 100, 200, …) holds a schedule-dependent vector under the parallel
  gradient; the shifted-capture memo `EPoint::GetECAmplitudeWithShift` keeps the level energy of a copy's first
  shifted evaluation (inferred, not covered by T1–T3).

## 2026-11-09 — Phase 4 S3 (a): results in a buffer, behind the old accessors (`feature/upgrade-4-results`)

**Done** (`data_model.md` S3 (a), decision M2 (c)). `include/Results.h` (plain, no engine headers): one record per
point and sub-point — σ, σ_E1, σ_E2, A_y, A_0, then (angular units) a count and max_ang_dist_order + 1 coefficients,
then (points of a sample line) a count and one term per compound nucleus. `EPoint`'s seven result members are gone
(`fitcrosssection_`, `fitE1/E2crosssection_`, `analyzing_power_`, `ang_dist_a0_`, `angularDists_`,
`natural_terms_`); the accessors read the record. The buffer is the `EData`'s (one per working copy):
`LayoutResults` at the end of `Initialize`, in `Clone` (the copy takes the source's values), and `EnsureResultsLayout`
at the start of every evaluation (`BeginEvaluation`, last step of `FillEnergyShiftsFromParams`), which re-lays only
when a point is not where the layout puts it (a re-placed sub-point grid, a copied segment). Probe copies in
`EPoint::Calculate` and the reaction-rate points take a private one-record buffer. `ec_shifted_` stays: a memo
keyed by exact energy (H6), not a result. `sizeof(EPoint)` 888 → 824 bytes; a unit is 24 bytes smaller overall.

**`AZURE_POISON_RESULTS=1`** (documented in `command_line.rst`): the buffer is quiet NaN on entry to each
evaluation. ctest `upgrade_poison_results`: T1–T3 calc / extrap / (T1) fit outputs and the API's arrays equal with
and without it. One known case, listed by the test: K1, the API's E1/E2 arrays of lines that are not
angle-integrated (never calculated there; 0 of a new record).

**Verified:** 179 × 3 (menu 1, 1 thread) identical to the reference build (2040 files, run.log included), and
identical again with `AZURE_POISON_RESULTS=1`; a 19-folder subset at 2 threads identical to 1 thread (only the
Coulomb-table lookup counter in run.log differs); ctest −LE slow 131/131; T1 fit (menu 2) and T3 fit
(`--use-gradient`) outputs identical to the reference, T3 also poisoned. Full T2 fit, 1 thread (bench invocation
without a call limit): χ² 294.173 reference and branch, also poisoned; 292.7 s / 293.5 s, peak RSS 281.3 / 276.7 MB
(run side by side); calc 1.06 / 1.07 s warm.

**Stale reads found by reading (not reached by the test projects):** (1) `EPoint::Calculate`, automatic tolerance,
`allPass` branch (`src/EPoint.cpp`, `this->SetFitCrossSection(bareValue)`): an E1/E2 data line (isDiff 5/6) with a
sub-point fold keeps σ_E1/σ_E2 of the last evaluation that folded fully — a value carried between evaluations, so
pool history (and the thread count) can enter a fit. (2) `ESegment::CalculateTheoreticalCrossSection` /
`TransformedValue`: a failed term returns NaN before `SetNaturalTerms`, so the sample decomposition file prints the
previous evaluation's terms for that point. (3) `MarkFailed` leaves NaN in σ_E1/σ_E2 of a non-capture point for
good (only the API shows them).

**Next:** S3 (b), split as proposed in the report (b1 engine calls take `Results&`; b2 owner moves to the caller;
b3 gradient and MCMC); the stale reads above to S4.

## Phase 3 merged; guardrails (written 2026-11-09, after the plan entry below)

**Merged into `feature/upgrade`:** `feature/upgrade-3-project`, **Phase 3** (design and record in `project.md`;
decisions P1–P5). 31 commits; pyazr −460 lines net (−980 / +517), `io/` +937 lines.
- `io/` (`azure2_io`, standalone like `physics/`): the lossless `io::Document`, the format conversion
  (`io::ConvertToFormat2`, format 2 in `<config>`), the engine's level grouping. Round trip byte for byte of every
  `.azr` under `tests/` and `examples/` in 945 variants (P1: the 48,922-file archive is on another machine).
- The engine reads every section from the document (`Config::SectionStream`); a format-1 file is converted in memory
  with a note; `AZURE2 --convert OLD NEW` and `pyazr.convert_azr` write it (P5: old files keep working).
- pyazr reads through the binding: document, segment / extrapolation / target-effect lines, level grouping, `<thm>`
  (`read_thm_block`; the 650-line Python port is deleted). The places where the two readers disagreed are decided
  as P2/P3 (`project.md` §5). pyazr needs the compiled module to read a project; 2.12.0.
- Step 4: `AzrModel` writes an unedited file back byte for byte in any layout or line ending; an edit changes only
  its lines (`azr_layout_test.py`).
- Step 5: one mass table. The engine's THM nuclei are generated from `pyazr/data/mass_1.mas20` (they were up to
  1 × 10⁻⁸ u off it; no output of the 179 projects moved) and checked by `pyazr_thm_nuclide_masses`; u, ħc, α, m_e
  have one definition (`physics/.../Constants.h`).
- Phase 2 leftover closed: the API's Park / abundance / scattering-length terms come from `Objective::ModelChi2`.

**Guardrails (`CLAUDE.md` "Recipes", `repo_consistency_test.py` 10–12):** a recipe per kind of addition (segment
keyword, `<config>` key or section, objective term, physics function, CLI flag), each with the check that fails when
it is bypassed: a constant copied with another value (10), a χ² term added outside `Objective.cpp` (11), the project
file opened outside `Config.cpp` (12). Each was shown to fail on an injected bypass. The fast pre-merge check is
`scripts/check.sh --quick` (consistency first, under a second).

**Verified** (branch head bba7f965 for the engine; the later commits are API, tests and docs): 179 × 3 — every output
identical except run.log, whose one changed line is the format-1 note; ctest 130/130 and the slow tests 10/10 on the
final head; `run_tests.sh` 133/133.

**Open, for the maintainers (they move results):** `EPoint.cpp`'s m_e is CODATA 2002 (0.510998903 MeV) in the
Fermi-function argument; `ECIntegral.cpp`'s proton mass 1.00727638 u (CODATA 2018: 1.007276466621). The first is a
listed exception of check 10 (which covers u, ħc, α, m_e).

**Open, found by the benchmark (not caused by Phase 3): fits depend on the thread count.** Same binary, T2, MIGRAD
with 200 calls: the final χ² is 0.306025731622686 / …729007276 / …733669139 (prior-covariance term) at 1 / 2 / 4
threads, each reproducible; the writer's total ends at 691.608 (1 thread) and 595.877 (8 threads; one of three
runs of the merged head ended at 691.608, so at 8 threads not even repeatable), because the fit
is unconverged and the 9th-digit differences grow over 200 steps. Present on the Phase 1 build (b-pathways) and the
Phase 2 head too; the 2026-11-07 entry deferred it to Phase 2 ("fixed-order reductions in the numerical gradient"),
which did not address it. Calculations (menu 1) stay identical at 1 and 8 threads. Suspect: the pooled working
copies (`AZURECalc::GetPooledCNuc/EData`) carry state between evaluations, and Minuit2's parallel numerical gradient
hands each copy a thread-dependent subset of calls. Next step: evaluate one parameter vector on a fresh copy and on a
pool copy after a different history, and diff. Belongs with Phase 6 (no pools) unless a smaller cause is found first.

**Benchmark** (this host, three interleaved runs each, rows in `tests/upgrade/bench/T2.tsv`; before = 93f8cbc5, the
start of Phase 1): calc 0.93 → 0.90 s (1 thread), 0.42 → 0.24 s (8); fit200 63.2 → 56.2 s (1), 36.4 → 30.7 s (8);
fit peak RSS 241 → 276 MB (1), 488 → 643 MB (8; the interim shift tables). The pre-upgrade `dev` engine (aab1de0d)
cannot run T2 (χ² 1.2 × 10⁸: keywords it does not know), so it is not a baseline. The 93f8cbc5 build also ends the
8-thread fit at 595.877: the thread dependence above predates Phase 1.

**Document for the owner:** `docs/upgrade/FOR_THE_OWNER.md` (two pages).

Thread check on the merged head (calculations, menu 1): 12/12 identical at 1 and 8 threads.

**Next.** The fit thread dependence above (diagnose before new features); then the first feature through the
recipes, extrapolation segments on a user or adaptive grid.

## 2026-11-10 — plan to the PI's meeting with the owner (next week)

**PI (maintainer):** the goal is to hand the owner a code that is efficient, easy to maintain and easy to extend, while
the owner's AI keeps adding features. Agreed:
- finish Phase 3 (`project.md`; P1–P5 decided);
- the `CLAUDE.md` guardrails: recipes for adding a feature, consistency checks that fail when one is bypassed, a fast
  pre-merge check, a speed benchmark;
- a short document for the owner.

**PI, 2026-11-10: no merge with `origin/dev` during the reconstruction;** it is planned only when everything is
complete. The assessment below is kept for then.

**Requested feature (PI, 2026-11-10), after Phase 3, as the first feature through the new recipes:** extrapolation
segments on a non-uniform grid. Proposed design (format 2, a `<segmentsTest>` keyword):
- `energies <file>`, or `energies E1 E2 ...` inline, replaces min/max/step;
- `adaptive <tol>` bisects, from the coarse uniform grid, every interval where the curve departs from linear
  interpolation by more than `tol` (relative), up to a point budget, seeded with the level energies and widths;
- reuse `AdaptiveIntegrationGrid` where it fits.

**`origin/dev` since the upgrade branched (2026-11-08, 18 commits, measured 2026-11-10):**
- documentation (skills), 10 commits;
- portable as they are: target-effect auto windows (442610f3, 3a104e06), UDR time-of-flight kernels at low energy
  (a68474f0), test and CI fixes (5a08c048);
- to redo in the new structure:
  - the heavy-ion Coulomb threshold fallback (56512d6a): it patches `external/coul`, now `physics/coulomb`; to be
    checked against mpmath;
  - the coincident-level analytic limits (6d3be80a), now in `physics/level_parameters`;
  - pyazr NUTS sampling (bf09666a, 7de653af), which should use `objective_and_grad()`;
- obsolete: the Qt logo (b4be0500; decision 3) and the `libcoul` GSL link (63f0c093; library gone).

## 2026-11-09 (late) — Phase 2 (one objective) merged; band and bound channel functions merged

**Merged into `feature/upgrade`:**
- `shift-band`: ShftFunc interpolates only beyond η_b = 1000; one number moved, 1.2 × 10⁻¹¹ (`coulomb.md` §1.4).
- `bound-functions`: S(0) from its own Bessel K (2 × 10⁻¹⁵ against mpmath); GSL's central derivative ported bit for
  bit; outputs identical.
- `upgrade-2-objective`: **Phase 2**, design and record in `objective.md`. One ordered list of terms (`Objective`)
  for the FCN, the finite-difference and LM cost, the MCMC likelihood (data part) and the API.
  - One enumeration of the diagonal priors serves the gradients, G2, the least-squares rows and the MCMC; no
    parameter is matched by name any more.
  - pyazr computes no prior (`objective()`, `penalties()` and the new `objective_and_grad()` call the engine;
    pyazr 2.11.0).
  - The PI accepted **D1–D7** and decided that **the MCMC shares the fit's configuration**. As built:
    - D1: the MCMC profiles THM norms;
    - D2, D7: Park's wall and the abundance sum-to-one are least-squares rows;
    - D3: the prior covariance replaces its parameters' own priors;
    - D4: no norm prior for a profiled THM norm;
    - D5: a fixed parameter pays no prior;
    - D6: nuisance priors in the API, pyazr, the writer and the MCMC.
  - Every changed output and test is classified in `objective.md` §4. Of note, D1 corrects a diverging THM chain
    (level energy −3.2 × 10⁷ MeV before, 6.29 MeV after).

**Merged head 8f015075 verified as a whole:** 179 × 3 identical to the objective branch's head; ctest 132/132 with the
slow tests; `run_tests.sh` 133/133; 1 vs 8 threads identical.

**Mistake and fix this session:** a new pyazr test was written over the existing `tests/pyazr/objective_test.py`
(86ac2970); it is restored byte-identical (e8468cdc), and the new test is `objective_terms_test.py`.

**Next.** (1) Phase 2 leftovers (`objective.md` §4): the writer from the terms instead of stored totals; `inference/` as
its own library once Phase 4 makes the data values. (2) Phase 3 (`io/Project`: one `.azr` reader/writer). (3) Step 6b
with the observables. (4) Lecture 2.

## 2026-11-09 (night) — six branches merged; step 6a merged; band step under verification

**PI (2026-11-09): "go on with the steps and you can merge yourself, no need for my review."** Merged into
`feature/upgrade` (not pushed): `external-capture` (step 4), `zero-width-fixed`, `whittaker`, `brune-iteration`,
`whittaker-engine`, `coupling` (step 5), in that order. Conflicts: test lists in `tests/pyazr/CMakeLists.txt` (all
kept) and pyazr's version (2.10.0). **After this merge everyone rebuilds the Coulomb table (format 2).** The merged head
db972c84 was verified as a whole:
- 179 × 3 against the `whittaker-engine` run (format-2 table): the only changes are those the branches recorded
  (step 5's two rounding files, `screening_6li` formal, the `brune-iteration` formal run.log and 10⁻¹¹ changes);
- ctest 127/127 including the slow tests; `run_tests.sh` 133/133; 1 vs 8 threads identical.

**Step 6a** (`feature/upgrade-1-pathways`, merged 87721b0c; `cnuc_split.md` §10 and §6). Step 6 as planned would
rewrite about 730 accessor calls in 23 files, which is Phase 2's work, so it is split. 6a moves the pathway
enumeration (`SortPathways`) into `physics/Pathways` (`reaction_pathways` of the System and the external-capture
final states); the engine's trees are filled from it. Evidence:
- the trees match a verbatim copy of the old code on all 179 `.azr` in two bases (51524 nodes);
- 179 × 3 bit-identical; ctest 119/119, run_tests 133/133, threads identical.

`CNuc.cpp` is down to 2913 lines. 6b (the consumers of the trees) goes with Phase 2's observables.

**ShftFunc band** (`feature/upgrade-1-shift-band`, not merged; design in `coulomb.md` §1.4, "The interpolation band
near threshold"). Measured: `whittaker_scaled` is fine near threshold but scatters at 10⁻⁹ beyond η_b ≈ 10⁴. So the
η_b ≤ 1000 limit stays, and the GSL-only condition (η_b > 50 with z < 0.05) is dropped. B_min for p+¹²C goes from
3.3 × 10⁻⁴ to 8.2 × 10⁻⁷ MeV. Verification is running (179 × 3 with and without a table, ctest, threads, run_tests).

**Next.** (1) Band: classify, merge. (2) The rest of `CNuc`'s physics: `CalcShiftFunctions` (Park's numerical S′ and
S″, `g_shiftDerivMemo`) onto `channel_functions`. (3) Phase 2 (`inference/` objective) and 6b.

## 2026-11-09 (evening) — split `CNuc`, step 5: coupling coefficients of the observables

Branch `feature/upgrade-1-coupling` from 9be284f0 (step 4), not merged. Commits: design (`cnuc_split.md` §9,
da3703d2, alone), physics (96917f36), engine (563a7172), record (9143239d). `CalcAngularDists`,
`CalcCaptureAnalyzingPower` and `ThmAngular`'s Z̄ take their coefficients from `physics/ObservableCoupling` (no
GSL; new `wigner_9j`, `racah_w`); `AngCoeff::Wigner9j` deleted; `CNuc.cpp` 3135 → 3061 lines. `SortPathways`
moves to step 6 with the pathway trees.

**Evidence.** The physics tests check the 9j and, independently of the code under test, the particle and photon
coefficients against m-scheme angular distributions (237k checks, ≤ 1.4 × 10⁻¹⁴). The engine-side test checks
every coefficient against the old GSL expression (≤ 4 × 10⁻¹⁵ relative, same terms kept). 179 × 3: 2063/2069
identical; the six files left (`ec_e1m1_phase` s_1, `gamma_lab_frame` cm, in three bases) differ by at most
3.2 × 10⁻¹⁵ of the column maximum. The sign flips in `ec_e1m1_phase` are 10⁻¹⁹ residues of a cross section that
is exactly zero at 0° and 180°. ctest 113/113, run_tests 133/133, threads identical, timing unchanged. **Ready
for PI review/merge**; independent of the four branches waiting.

**Next.** (1) PI: merge order of the five branches (`zero-width-fixed`, `whittaker`, `brune-iteration`,
`whittaker-engine`, `coupling`). (2) Step 6: data accessors and the pathway trees (`SortPathways`) on the System.
(3) `ShftFunc` band; lecture 2.

## 2026-11-09 (day) — PI decisions; the Brune blocker fixed; step 3b (ii)

**PI decisions.** (1) `brune-iteration`: option (i), the iteration must converge where the root exists.
(2) afd8ee6a (`apply_fit_test.py` to §9's rule) confirmed. (3) The Whittaker engine branch waits for step 3b (ii)
— now done, see below.

| branch | state |
|---|---|
| `feature/upgrade-1-brune-iteration` | **verified, blocker fixed** (`level_parameters.md` §10, design a3bb203a + d6afeaeb + ff484603, code 5d8dc1e1, verified e4ca5d9d). E − e(E) is strictly increasing (Weyl, S′ ≥ 0), so each level has one Brune energy; the fixed-point map repels where e′ < −1 (−2.88 on the blocker: 2-cycle −0.973 / 4.462 MeV around the root 2.36820 MeV). A level that ends unconverged is now the root by bracket + bisection (status `Bracketed`). ctest 111/111, run_tests 133/133, threads identical; 179 × 3: only `screening_6li` formal changes (three background levels now transformed, equal to the Brune input digit for digit) |
| `feature/upgrade-1-whittaker-engine` | step 3b (ii) done (design cf5ddcce, code 82ca0985, verified b70ffdf9): final-state W memoized on its exact arguments, bit-identical (2040/2040), 99.8 % hits; T1 5.14 → 2.77 s, 13N 9.71 → 7.08 s (GSL build: 2.89 / 7.60 s). **Ready for PI review/merge**; after the merge everyone rebuilds the Coulomb table (format 2) |

**Next.** (1) PI: merge order (`zero-width-fixed`, `whittaker`, `brune-iteration`, then `whittaker-engine`).
(2) Narrowing `ShftFunc`'s band near threshold (own comparison). (3) Lecture 2.

## 2026-11-09 (night) — the three branches built and verified; the Whittaker engine part on a new branch; lecture 1

Work moved to a 24-core node (local builds, no Condor); at most 12 threads. Environment: conda env with Python 3.12,
numpy, scipy, mpmath, pybind11; build `cmake ... -DBUILD_GUI=OFF -DUSE_API=ON -Dpybind11_DIR=...` (34 s on 12 jobs).
The 179 × 3 comparison takes 51 s with a Coulomb table and ~8 min without; `run_tests.sh` ~20 min on 11 threads.

| branch | state |
|---|---|
| `feature/upgrade-1-external-capture` | `run_tests.sh` rerun: 133/133 (pending item closed); `ctest -LE slow` 111/111 |
| `feature/upgrade-1-zero-width-fixed` | **verified** (`external_capture.md` §9, "Verified"): 2040/2040 identical, ctest 112/112, run_tests 133/133, threads identical. Two test commits: check (c) of the new test called `apply_fit` in a way it refuses (fixed); `apply_fit_test.py` "fixed nonzero, no line: still refused" changed to §9's rule — **PI to review** (afd8ee6a) |
| `feature/upgrade-1-whittaker` | **verified** (`coulomb.md` §1.4, "Part 1 verified"): 2040/2040 identical, ctest 114/114 with the binding test |
| `feature/upgrade-1-brune-iteration` | verified except one test — **blocker for the PI** (`level_parameters.md` §9, "Blocker"): `pyazr_formal_capture_mask` (c) is NaN because `brune_from_formal` 2-cycles (E = −0.973 / 4.462 MeV) at ±10⁻⁶ MeV on `energy_3`; a18fe89e hid the same non-convergence by keeping formal values. Otherwise 179 × 3: Brune, Park identical, formal 657/662 (1.7 × 10⁻¹¹); run_tests 133/133 |
| `feature/upgrade-1-whittaker-engine` (new, from `-whittaker`) | step 3b (i) part 2: design addendum (a2593c96), engine switch (48fb454f), EC signature v6 with a stale-file case (5ed85322). **Changes outputs; PI to review before merging.** Numbers below |
| `docs/lectures-1` (new, from 7704e48f) | lecture 1 text written (`docs/lectures/01-external-region/lecture.md`), figures redone (no label boxes over curves, ranges fitted) |

**Whittaker engine part (`feature/upgrade-1-whittaker-engine`).** `WhitFunc`, `CoulombTable::ExactBound` and
`ThmBoundStateTail` call `whittaker` / `whittaker_scaled`; no engine file includes `gsl_sf_hyperg.h`. The Coulomb
table format is 2 (its bound sheet changed): a format-1 file is refused with a warning, so **everyone must rebuild
their table after the merge** (`AZURE2 --build-coulomb-table`, 112 s on 12 threads). Comparison without a table
against the `-whittaker` head: Brune 595/689, formal 559/662, Park 605/689 identical; every observable within
2.2 × 10⁻⁹ of its column's largest value; parameter files at most 9 × 10⁻⁸ relative (last printed digits), the
apparent large differences being values at rounding (10⁻¹⁸–10⁻³⁷ formal zero-input widths; 10⁻²² cross sections
at 0°/180° that vanish), one 0.99999999 eV width printed as 1000 meV, and −0.000000 signs; the 20 `intEC.dat`
caches change (hence v6). With tables built by each build: `c12c12` identical in the three bases. Note: the
format-1 table in `~/.azure2` on this node (26 Oct) moves `c12c12` by 1.9 × 10⁻⁵ against a format-1 table
built by today's reference build — the table can go stale within one format number. ctest 114/114 (with a format-2
table; `pyazr_park_dsde` needs a usable table, it fails when the only table is refused); threads identical;
run_tests 133/133. **Cost:** one-thread calc 13N 7.10 → 9.14 s, T1 2.62 → 4.73 s, all of it in
set-up (the external-capture table, 1176 integrals at 1320 nodes in T1, recomputed every run; evaluation unchanged,
~1 ms) — the 13× slower W per call. Step 3b (ii), one W table per integral, is the remedy; do it before merging.

**Next.** (1) PI: the brune-iteration blocker (safeguarded iteration, or move the test model); afd8ee6a; whether the
engine switch may merge with its cost or waits for 3b (ii). (2) Step 3b (ii): the per-integral W table in
`ECIntegral`'s data path, then the timing again. (3) Narrowing `ShftFunc`'s band near threshold (own comparison).
(4) Lecture 2 (channels, level and collision matrices, the three bases).

## 2026-11-08, end of session — lectures added to the plan; work stopped here

The PI asked for R-matrix lectures (Notion Markdown, figures, later slides): plan and lecture 1's figures in
`docs/lectures/` (`PLAN.md`). Lecture 1's text is not written. Everything else is as in the entry below: four
branches, the three on top of step 4 BUILD-PENDING (Condor schedd on crcfe02 down).

## 2026-11-08, later — three branches on top of `feature/upgrade-1-external-capture`; builds pending

The PI accepted step 4's fixes and its two proposed defaults (2026-11-08). Each is a branch cut from
`feature/upgrade-1-external-capture` a18fe89e, design note first, then code; they do not depend on each other.
Merge order: this branch, then the three in any order. **Not built yet**: the login node's Condor schedd has been
down since the node rebooted (~07:40); the code is syntax-checked (`g++ -fsyntax-only` with the build's flags, no
warnings in the touched lines) and each code commit says "BUILD AND TESTS PENDING".

- `feature/upgrade-1-brune-iteration` (worktree `azure2-up-ec`): `brune_from_formal` inside an evaluation converges
  at max(10⁻¹² MeV, 10⁻¹²|E|), takes its last iterate at the rounding floor (≤ 10⁻⁹ MeV, status `RoundingLimited`),
  and gives NaN instead of the formal values otherwise; setup failures stop the run (`level_parameters.md` §9).
  Tests: `level_parameters_unit_test` (new cases), `pyazr_formal_brune_iteration` (slow; χ² of
  `beam_profile_kernel` smooth in `width_9_2`, which jumped 40.44 → 42.93).
- `feature/upgrade-1-zero-width-fixed` (worktree `azure2-up-zero`): a width given as 0 is fixed in every basis
  (`ZeroInputWidth`, not the transformed value); `apply_fit` skips fixed unlisted channels; pyazr 2.9.1
  (`external_capture.md` §9). Test `pyazr_zero_width_fixed`. Expected: formal free set = Brune free set (134 → 67 on
  `beam_profile_kernel`), `--no-ignore-externals` input-0 photon widths fixed.
- `feature/upgrade-1-whittaker` (worktree `azure2-up-whit`), step 3b (i) part 1: `physics/coulomb` `whittaker`,
  `whittaker_scaled` without GSL (integral representation of U, double-exponential quadrature; `coulomb.md` §1.4).
  Against mpmath on 1162 points: W 2.1 × 10⁻¹², ρW′/W 3.3 × 10⁻¹³ (GSL as the engine calls it: 6.9 × 10⁻⁹, 10 points
  without a value); COULCC at imaginary argument dropped (0/0 at η_b = l+1, l+2, …). 4.8 µs per call against GSL's
  0.36 µs. `physics_whittaker` and `whittaker_gsl` were run standalone (pass; four mutations each fail it); the
  binding `_azure2.whittaker` and `pyazr.whittaker_wave` (pyazr 2.10.0) wait for the build.

**To do when Condor is back** (per branch): build (`upgrade-work/build.sh SRC BUILD OFF`), `ctest -LE slow`, the new
slow test, `run_tests.sh`, 1 vs 8 threads, the 179-project comparison against a18fe89e (`ec-work/runall.sh`,
`classify.py`), each change classified; then a follow-up commit with the numbers. Also still pending:
`run_tests.sh` on this branch's 709ed040 build (`ec-work/verify.sh`).

**Next (step 3b (i) part 2, design in `coulomb.md` §1.4):** `WhitFunc` over `whittaker`/`whittaker_scaled`,
`ShftFunc`'s bound branch on ρW′/W (whether the interpolated band near threshold narrows is a separate comparison),
`CoulombTable`, `ThmNumerics`, `ThmFunc`; measure `ECIntegral`'s data path (13× per W call) before and after, and
route `pyazr.whittaker` (transform.py) through the engine. Bound-state outputs move in their last digits; classify.

## 2026-11-08 (`feature/upgrade-1-external-capture`) — split `CNuc`, step 4: `physics/external_capture`; the mask on the input

Branch `feature/upgrade-1-external-capture` from `feature/upgrade` 5055bec0, not merged. Commits: the design
(`external_capture.md`, alone), the physics part, the engine through the adapter, the mask on the input, the
lowest-level final state with a free energy (step 2 open item 2); then the addendum (section 7, alone), one adjoint
test for every gradient path with decoupled masked levels on P, `tests/pyazr/formal_capture_mask_test.py`, these docs.

**State.** `physics/external_capture`: N_f, the multipole rule, γ^ext (Eq. 4, radial integral injected), the
projector Π = T⁻¹PT; the adapter `src/CNucExternalCapture.cpp` replaces `CalcExternalWidth` and six copies of
N_f^−½ (`CNuc.cpp` 3349 → 3132 lines). `--ignore-externals` masks a photon channel whose input width is exactly zero
(`ALevel::ZeroInputWidth`), the same in every basis and continuous in a fit; the formal basis carries the Brune mask
with Π. The adjoint does not differentiate Π or the step-2 final-state amplitudes: `FormalPartsOutsideAdjoint`
(`src/AZUREGrad.cpp`) sends the χ² gradient, the residual Jacobian and the model gradients to finite differences
when a slot has a projector or a masked level has a free particle amplitude, or a formal final-state block exists;
a masked level with no particle width is decoupled and uses P (exact). How the parts stay apart: the physics part
gets arrays and returns arrays; the adapter only reads `CNuc` and calls it; the gradient asks `CNuc` one question.

**Evidence** (`external_capture.md` §6). 179 `.azr` × 3 bases against 5055bec0: Brune 689/689, Park 689/689
identical; formal 650/662, every change classified (`upgrade/T1` capture and χ², `beam_profile_kernel` capture,
χ² and three input-0 E2 widths), each changed formal output now equal to the Brune one. The guard commit:
2040/2040 files byte-identical to 3d08b7f0. χ² gradient vs central differences, formal basis, 13N with a masked
coupled 1/2⁺ level: 1.7 × 10⁻³ → 1.4 × 10⁻⁷ (new test; fails on e205038e; other checks fail on 5055bec0). No
warnings; check.sh PASS; `ctest -LE slow` 111/111 (`pyazr_zero_width_level` passes again); `run_tests.sh` 133/133 on
3d08b7f0 (on the final build cut short by the 2026-11-08 node reboot, rerun pending; its
menu-1 outputs are byte-identical to 3d08b7f0's);
1 vs 8 threads byte-identical (T1–T3, 13N, three bases); T2 bench (Brune; mean of three, same host, 5055bec0 →
3d08b7f0; the
guard is a no-op in Brune) fit200 61.58 → 62.16 s (1 thread, +0.9 %), 26.61 → 26.65 s (8); calc 1.03 → 0.98 s,
0.21 → 0.20 s.

**Open (proposed default in brackets).**
- **Blocker for merging the formal basis as "continuous":** `brune_from_formal` at 10⁻¹² MeV is at the rounding
  limit (§8). On `beam_profile_kernel` a 10⁻⁶ change of one width makes the 20 MeV 1⁻ level miss the tolerance
  (last step 7.7 × 10⁻¹² MeV after 1001 iterations); the formal values are then kept silently, its row of T becomes
  a unit vector and χ² jumps 40.44 → 42.93. The merged step-2 final-state amplitudes use the same call. [A
  tolerance relative to the rounding of E, e.g. max(10⁻¹² MeV, 64 ε |E|), in both callers; a level that still does
  not converge is an error of the evaluation, not formal values. Needs the PI: it changes merged formal-basis
  numbers at ~10⁻¹¹ and the no-silent-fallback rule decides the second half.]
- Width parameters with value 0 are fixed at setup (`CNuc.cpp:2259`), tested on the transformed value, so in the
  formal basis rounding decides: `beam_profile_kernel` frees 35 parameters with |value| < 10⁻¹² on 5055bec0 and 48
  here (10 new on the lowest 1⁻, 2⁺, 3⁻ levels, from 8df95d22's reordered sums; 3 masked photon channels now exactly
  γ^int = 0 plus rounding). [Fix on the input, as the mask: a width given as 0 is fixed in every basis. Changes which
  parameters a formal-basis fit varies; PI decision.]
- Formal basis with a projector or a final-state block: finite-difference gradients. [Accept; the derivative of T
  is a Phase 2 item.]
- From step 2, unchanged: `TransformOut`'s 10⁻⁶ MeV report iteration; `ThmLevelWidth`'s own total width.

**Next.** PI decisions on the two items above; then merge; then step 3b (i) (Whittaker functions without GSL).

## 2026-11-08 (`feature/upgrade-1-cnuc-split-2`) — split `CNuc`, step 2: observed quantities; formal-basis external capture

Branch `feature/upgrade-1-cnuc-split-2` from `feature/upgrade` 9d50667f, not merged. Commits: the design
(`cnuc_split.md` §8, alone), `observed_levels` and its callers, the formal-basis external-capture fix with its test,
the restriction to final states above the lowest level, these docs.

**State.** `level_parameters::observed_levels` is the one definition of the observed quantities from total amplitudes
γ^int + γ^ext: W = s 2P|γ^tot|²/N and the `<levels>` value (Γ eV, C fm^−½, μ, Q, B, or the Brune amplitude of a
`gammaIsRWA` column) with its kind. `TransformOut` = `TransformToBrune` + `ObservedWidths`; `GetTransformParams`,
`PrintTransformParams`, `CheckRadiativeWidths` read the stored values (the last now with N and the external part);
`ConvertAmplitudeBasis` calls `park_from_brune`/`brune_from_park`. The five γ → Γ copies are deleted. External capture
takes the final state's Brune amplitude and energy in every basis (`ALevel::ECGamma`, `ECEnergy`; formal basis: a
read-only setup record of the blocks with a final state above their lowest level, `brune_from_formal` to 10⁻¹² MeV in
`FillCompoundFromParams`; γ^ext subtracted before any group is made formal); the intEC signature hashes the energy
used. How the parts stay apart: `level_parameters` gets arrays (E, γ, γ^ext, flags, channel functions) and returns
arrays; the engine only gathers and stores; external capture still supplies γ^ext from `CNuc` (step 4).

**Evidence** (`cnuc_split.md` §6, step 2). 179 `.azr` × 3 bases against 9d50667f: Brune 689/689, Park 689/689
identical; formal 649/662, every change classified (8 `parameters.out` μ/Q instead of G = 0, 3 `screening_6li` logs
lose a warning built on formal amplitudes, `beam_profile_kernel` `param.par` and `parameters.out` whose photon widths
to the 6.049 MeV state now reproduce the input, e.g. 0.000053 → 0.000330 meV for 0.33). No re-baseline. Formal ↔ Brune
capture on `beam_profile_kernel` (111 energies, `--no-ignore-externals`): 6.049 MeV state 0.808 → 1.5 × 10⁻¹¹; the other
four final states 0 before and after. New `tests/beam_profile_kernel/check.sh` (fails on 9d50667f, 0.448). No warnings;
`check.sh` PASS; `ctest -LE slow` 109/109; `run_tests.sh` 133/133; 1 vs 8 threads byte-identical (T1–T3, three bases);
T2 bench (mean of three, same host) fit200 133.1 → 131.8 s (1 thread), 55.8 → 56.2 s (8, +0.7 %); calc 2.45 → 2.22 s,
0.42 → 0.44 s.

**Open (proposed default in brackets).**
- `--ignore-externals` (the default) drops the external resonant capture of levels whose internal photon amplitude is
  < 10⁻⁸ in the basis in use, so formal and Brune still differ by default (9.7 × 10⁻² on the `beam_profile_kernel`
  data). [Step 4: define the mask on the Brune amplitudes in every basis.]
- A formal-basis fit that moves the energy of a lowest-level final state: its Brune amplitude then differs from the
  formal one at second order and is not recomputed. [Accept until step 4; or record every block with a final state
  whose energy is free.]
- Formal basis with a final state above the lowest level: the analytic gradient falls back to finite differences.
  [Accept; the derivative of `brune_from_formal` is a Phase 2 item.]
- The reports (`parameters.out`) of the formal basis keep `TransformOut`'s 10⁻⁶ MeV iteration; deeply bound channels
  amplify it (an ANC input 0 reported as 38.65 fm^−½ in `beam_profile_kernel`, as before). [Tighten to 10⁻¹⁰ MeV in a
  later step, with its own comparison.]
- `ThmLevelWidth` (THM line shape) has its own total-width formula. [Phase 5 with THM.]

**Next.** Merge; then step 4 (`physics/external_capture`) or step 3b (i).

## 2026-11-07 (director) — step 3b merged (4a8ecf02)

Per-stage channel table and the Coulomb memo's deletion are on `feature/upgrade`. Accepted as interim, with the
reason: 13N +1–9 % per pair (inside the base's spread; revisit in step 4, `physics/external_capture`); T2 fit peak
memory +30 % at 8 threads from the three shift-state tables per segment (gone with Phase 2's analytic shift
derivative); the χ² of a fit differs in the 8th digit between 1 and 8 threads, already so on the base (Phase 2:
fixed-order reductions in the numerical gradient). Next: step 2 (observed widths from the total amplitude), step 4
(external capture), then the Whittaker functions without GSL and `--coulomb-method`.

## 2026-11-07 (`feature/upgrade-1-channel-table`) — step 3b: the per-stage channel table; memo deleted

Branch `feature/upgrade-1-channel-table` from `feature/upgrade` 815d9b5f, not merged. Commits: the design
(`channel_table.md`, alone, with the `perf` measurement of where `-nomemo` spent its time), the table, the radial
scratch, the `-nomemo` commits re-cut (c3dc22e5, 37a3637b), the entrance phase through the scratch, keys unique per
point and the table released after setup, the first docs; then (director's answers) the free-shift tables, the
capture amplitude from the table, distinct recorded rows, these docs.

**State.** `ChannelTable` (src/): exact keys (pair α, l, E_c), rows P, S, e^{−iφ}, e^{i(σ_l − σ_0)} from one function
(`ChannelTable::Compute`), built by `EData` before the point loops (keys from `EPoint::ChannelKeys`), one thread per
row, read-only, released after setup; at setup the capture amplitude reads √P and the entrance phase at a_c from it.
Free-shift refills: each segment keeps the tables of its last three shift states (key: the shift values bitwise;
`BeginShiftRefill`/`CommitShiftTable` single-threaded around the parallel loop; replaced, never edited) — interim,
deleted when the objective has the analytic shift derivative (Phase 2). `CoulFunc::SetRadialScratch`: thread-private
exact-argument scratch for the capture integrals (accepted under RULES §2.3, `channel_table.md` §3).
`CoulFuncCache` deleted. Stay: the per-slot rows of `EDependentTables` (Phase 4), `g_shiftDerivMemo` (until
`level_parameters` computes S, S′, S″ at E_λ).

**Evidence** (`channel_table.md` §6). 179 `.azr` × 3 bases against the memo base: observables ≤ 7.0 × 10⁻¹¹ (one
exact zero's residue), `intEC.dat` ≤ 4.0 × 10⁻¹¹; against `-nomemo` 670/689 byte-identical, 19 `intEC.dat`
≤ 2.8 × 10⁻¹³ (one ulp of E_c). T2 fit200 (1 thread) byte-identical to `-nomemo`. No re-baseline. No warnings;
`channel_table_unit`; check.sh PASS; `ctest -LE slow` 109/109; `run_tests.sh` 132/132; THM
byte checks 5/5; 1 vs 8 threads byte-identical on T1, T2, T3, 13N (menu 1); TSan (upgrade section) 0 reports (first
build).

**Speed** (qa-xp-003, alternating, three repeats; base → branch): T2 calc 2.29 → 2.30 s (1 thread), 0.68 → 0.37 s
(8); fit200 140.1 → 121.7 s (1), 68.8 → 53.3 s (8); 13N 15.52 → 15.76 s; natural_neon 21.48 → 18.66 s; `tgt_fixed`
57.65 → 57.41 s. 13N is slower in each of five pairs (+1 to +9 %) but inside the base's spread. Peak RSS of a T2 fit
+13 % (1 thread), +30 % (8): the shift tables.

**Open (proposed default in brackets).**
- 13N +1–9 % per pair. [Accept; the rest goes with step 4 (`physics/external_capture`, unit amplitudes per grid
  energy).]
- Fit peak RSS +30 % at 8 threads (three shift tables per segment and working copy). [Accept until Phase 2 removes
  them; or keep two states, which loses the return to ν₀ after ν₀ ± h.]
- A fit at 1 and 8 threads differs in χ²'s 8th digit — on the memo base and `-nomemo` too. [Separate item for the
  Phase 2 objective (fixed-order reductions in the numerical gradient).]

**Next.** Merge, then step 3b (i), (iii), (iv) or step 2.

## 2026-11-07 (`feature/upgrade-1-coulomb`) — split `CNuc`, step 3: Coulomb functions from pycoulomb (memo kept)

Branch `feature/upgrade-1-coulomb` from `feature/upgrade` f549fa23, not merged. Commits: the design (`coulomb.md`,
alone), `physics/coulomb`, the GSL cross-check, the engine routed and `external/coul` deleted, `channel_functions`
in physics/, the engine calling it, these docs. The memo's deletion (c3dc22e5, with its docs) is on
`feature/upgrade-1-coulomb-nomemo`, to be re-cut on the per-stage table (director's decision; `cnuc_split.md` §3,
step 3b ii′).

**State.** `physics/coulomb`: `coulomb(l, η, ρ, method)` → F, F′, G, G′, σ_l, pure, NaN with a status on failure;
methods coulsm (default; AZURE2's former `external/coul`), coul, coulfg, coulcc, one adaptor file each, sources
byte for byte from pycoulomb azure2-prereqs 7adef98 and coulcc-cpp 37f0f9b (`physics/coulomb/README.md`).
`physics/ChannelFunctions`: P, S, e^{−iφ}, S′ (Richardson), σ_l of an open channel and η, ρ (non-relativistic and
relativistic) with the constants; the engine's `Kinematics`, `Constants.h` and `CoulFunc` call them, and every
Coulomb evaluation goes through `CoulombWavesAt` → `coulomb(…, coulsm)`, behind the unchanged `CoulFuncCache`.
Deleted: `external/coul` (8 files, 4053 lines). Not done (step 3b): Whittaker functions without GSL (`WhitFunc`,
`ShftFunc`, `ChannelFunc`'s bound branch stay), the per-stage table, then the memo's deletion; hybrid channels,
`CoulombTable`, `--gsl-coul`.

**Evidence** (`cnuc_split.md` §6, step 3). 179 `.azr` × 3 bases. Routing: identical except the hybrid-potential
projects (cwfcomp's 10⁻¹⁵ imaginary parts no longer enter the hybrid H⁺): observables ≤ 9.9 × 10⁻¹¹, one formal
width 5.0 × 10⁻⁸. `channel_functions`: bit-identical (2040 files). No re-baseline. physics_coulomb (mpmath, 1296
points): coulsm/coul ≤ 8.3 × 10⁻¹⁰ where the functions oscillate, coulfg 6.3 × 10⁻¹⁴, coulcc 1.5 × 10⁻¹²
everywhere; identity 8.3 × 10⁻¹⁰; σ_l 7.1 × 10⁻¹⁵; 8 threads bit-identical. GSL cross-check 4.6 × 10⁻¹⁰ (448
points). physics_channel_functions: neutral closed forms 6.0 × 10⁻¹⁰, P/S/φ against mpmath ≤ 7.5 × 10⁻⁹, S′ against
coulcc 1.8 × 10⁻⁷.
Checks (GCC 11.5, GUI off, Condor, build of 08516cfb): no compiler warnings; `tests/upgrade/check.sh` PASS (goldens
unchanged); `ctest -LE slow` 110/110 — `pyazr_thm_profiled` and `pyazr_thm_fixed_angle` included, and both also pass
on f549fa23 in 2 of 2 repeats on the same host, so their failure in the director's run is intermittent (the memo is
per thread and call-order dependent); `tests/run_tests.sh` 131/132, the failure the known `shared_point_tables` (d)
peak-RSS ratio (316 MB > 3 × 105 MB). 1 vs 8 threads: T2 and T3 outputs byte-identical; T1 and 13N identical except
`intEC.dat` (≤ 3.4 × 10⁻¹⁴ and ≤ 5.1 × 10⁻¹² of the line's largest value): the memo's per-thread caches, gone on
`-nomemo`, where all four were byte-identical.

**Speed** (T2.tsv rows; `coulomb.md` §7). host d32cepyc086, mean of two, before (f549fa23) → this branch: T2 calc 1.56
→ 1.36 s (1 thread; one
f549fa23 run at 1.89 s), 0.44 → 0.42 s (8); T2 fit200 80.5 → 81.9 s (1 thread), 36.5 → 37.0 s (8); menu-1
calculation, one thread, `13N` 10.05 → 10.16 s, `natural_neon` 12.35 → 12.45 s; `convolution_unify` `tgt_fixed`
(2 threads) 29.9 → 30.6 s. The routing costs at most 2 % (fit and `tgt_fixed`), within the run-to-run spread of the
calculations; peak RSS within 2 %.

**Measured on `-nomemo`** (`coulomb.md` §7): outputs ≤ 7.0 × 10⁻¹¹ (one exact zero's residue), the two THM checks
pass (their failure is the memo's call-order dependence), but 3–17× slower (T2 calc 1.37 → 6.16 s, `natural_neon`
12.3 → 105.6 s) and `convolution_unify` `tgt_fixed` 84 → 1414 s, past its 600 s limit.

**Open (proposed default in brackets).**
- The THM identity checks and the 1-vs-8-thread `intEC.dat` differences are intermittent memo effects. [Accept until
  `-nomemo` is re-cut on the table; there both are deterministic by construction.]
- Whittaker functions without GSL. [COULCC at imaginary argument in `physics/coulomb`, against mpmath and GSL; bound
  outputs will move in their last digits and are classified like this step.]
- `--gsl-coul`, `CoulombTable`. [Add `--coulomb-method=coulsm|coul|coulfg|coulcc` and retire `--gsl-coul` with a
  message; keep `CoulombTable` until the per-stage table is measured against it.]
- coulfg's BSD-3 text (the LLNL COPYRIGHT file) is not in pycoulomb. [Jakub adds it in pycoulomb; copied here then.]

**Next.** Step 3b (ii), the per-stage table, then `-nomemo` re-cut on it; or step 2, which does not depend on it.
## 2026-11-07 (director, end of session) — where to continue

**Branches (local, bundled in `azure2-upgrade.bundle`; pycoulomb branches in `pycoulomb.bundle`).**
- `feature/upgrade`: Phases 0 and 1 so far — test projects T1–T3, goldens, sanitizers, `physics/` (System, JpiBlock,
  Convention, level_parameters), CNuc split steps 1 and 1b, Qt removed, decisions 1–12 recorded below.
- `feature/upgrade-1-coulomb` (not merged): pycoulomb's C++ in `physics/coulomb/` (coulsm default, coul, coulfg,
  coulcc), `channel_functions`, engine routed through them, `external/coul` deleted; output-identical except the
  hybrid-potential rounding (≤ 9.9 × 10⁻¹¹). Its docs/timing commit may be missing if the session ended first: check
  `git log`; if the head is 08516cfb, add the docs (coulomb.md describes the memo as still present) and time T2,
  13N, natural_neon, convolution_unify against f549fa23 before merging.
- `feature/upgrade-1-coulomb-nomemo` (hold): deletes `CoulFuncCache`; correct (≤ 7 × 10⁻¹¹, fixes the two pyazr THM
  byte checks) but 3–17× slower until the per-stage table exists. Re-cut on top of step 3b.
- pycoulomb (Jakub's repo): `azure2-prereqs` (libcoulomb, tests, LICENSING_NOTES.md) and `coulcc-cpp` (COULCC in
  C++, bit-identical to the Fortran). Jakub: add the LLNL BSD-3 COPYRIGHT for coulfg, choose the licence, push, tag.

**Next steps, in order.**
1. ~~Merge `feature/upgrade-1-coulomb`~~ Done (a9cd42c3): routing costs ≤ 2 % (within spread); outputs identical
   except hybrid rounding; `ctest -LE slow` 110/110. The two pyazr THM byte checks are intermittent (the per-thread
   memo's call order), not a merge defect; they become deterministic with step 3b.
2. Step 3b (`cnuc_split.md` §3, `coulomb.md` §6): the per-stage channel-function table — P, S, S′, φ, σ_l per
   (α, l, channel energy) over the evaluation energies (decision 11, exact energies), read-only during an evaluation,
   replaced when a free shift moves; replaces the `EDependentTables` Coulomb rows, `g_shiftDerivMemo` and the memo.
   Gate: outputs ≤ 10⁻⁹ or explained; T2 and `convolution_unify` faster than with the memo. Then re-cut and merge
   `-nomemo`.
3. Whittaker functions without GSL (COULCC at imaginary argument); `--coulomb-method`; retire `--gsl-coul`.
4. Remaining CNuc split steps (`cnuc_split.md` §3), then Phase 2.

**Known open items.** `shared_point_tables` (d) peak-RSS limit has no margin; `light_excitation` (a) external-capture
integrals depend on thread scheduling at 5 × 10⁻¹¹; formal basis with external capture disagrees with Brune (step-1b
entry, item 2); Sphinx never built on CRC (`sphinx-build -W docs/source` before pushing). Agent brief:
`~/AI-R-scratch/upgrade-work/AGENT_BRIEF.md` on CRC (rules are also in CLAUDE.md).

## 2026-11-07 (director) — state after the step-1b and Qt-removal merges: two pyazr THM checks fail

`feature/upgrade` 160439ba, fresh build, GCC 11.5, no Qt (Condor 16205, 16206): no warnings; `tests/upgrade/check.sh`
ok; `tests/run_tests.sh` 131/132, the failure the known `shared_point_tables` (d) peak-RSS ratio (305 MB > 3 × 90 MB;
the limit has no margin). `ctest -LE slow` 105/107: `pyazr_thm_profiled` and `pyazr_thm_fixed_angle` each fail one
check, "write_output_files: AZUREOut_aa=5_R=4.out identical to the CLI's" (χ² agrees to 10⁻⁹; the other output
files are identical). Both passed on each parent of the merges. Under diagnosis on `fix/thm-pyazr-cli-identity`;
the Coulomb-function step waits for it. Decisions 11 (exact energies) and 12 (data as model expressions) added.

## 2026-11-07 (`feature/upgrade-1-cnuc-split-1b`) — split `CNuc`, step 1b: decision 10 on, one solver, one channel energy

Branch `feature/upgrade-1-cnuc-split-1b` from `feature/upgrade` 4bf7aa06, not merged. Commits: the design
(`cnuc_split.md` §7, alone), decision 10 with the options deleted, the channel energy and adapter, the sign rule
corrected, a strongly mixed test, the docs, the benchmark rows, this entry.

**State.** The formal basis has two rules and no options: B_c = S_c at the lowest level of each J^π group, and each
formal level oriented like the Brune level it comes from (`level_parameters.md` §3, §4; `conventions.rst`).
`level_parameters` uses its own eigensolvers only; `BoundaryRule`, `FormalSignRule`, the solver options and
`GslGeneralizedEigen`/`GslSymmetricEigen` are deleted. The adapter forms E_c = (E_x − S) − e, as at a data point, and
returns S, P, S′ at every call (`ChannelEnergyOrder`, `ChannelQuantities` deleted). In the formal basis the engine
prints a note for each group whose first listed level is not its lowest.

**Evidence** (`cnuc_split.md` §6, step 1b). 179 `.azr` × 3 bases, 4bf7aa06 vs 7f4af9d7: Brune 670/689 files
identical, Park 685/689, formal 595/662. Class (a), formal basis: 21 `param.par` (sign, B_c in `natural_neon` and
`be7p_fit`), 25 `parameters.out`, observables ≤ 8.5 × 10⁻¹¹. Class (b), last digits from the E_c order: Brune
observables ≤ 9.8 × 10⁻¹¹ printed (3.6 × 10⁻¹¹ full precision), Park ≤ 6.7 × 10⁻¹¹; Brune amplitudes up to 6.3 × 10⁻⁹ in
the ⁷Li + p models' 0⁺ background level (negative Brune denominator). No class (c). Goldens and `expected/` unchanged
(T2 fit within its own run-to-run spread); no re-baseline. GUI off, Condor: no warnings; `check.sh` PASS;
`ctest -LE slow` 107/107; `run_tests.sh` 131/132 (`shared_point_tables` (d), the flaky peak-RSS ratio of the step-1
entry: the calculation's RSS varies 67–105 MB between runs; fit RSS 301–308 MB here vs 310–312 MB before; 3 repeats:
1 failure here, 0 before). T2 bench, same host, mean of two: fit200 78.54 vs 78.53 s
(1 thread), 34.30 vs 34.36 s (8); calc 1.16 vs 1.21 s, 0.39 vs 0.39 s.

**Found.** (1) Decision 10's "largest |γ| positive" changes observables: it reversed the 6.917 MeV 2⁺ external-capture
final state of `beam_profile_kernel` (largest formal amplitude a photon one) and the formal-basis capture to it changed
by up to 4.51 relative. A final state's sign is shared with the photon amplitudes that end on it and with the
external-capture amplitude, so it is not free. Built instead: the orientation rule, also solver-independent. (2) The
formal and Brune bases disagree in capture with external capture, before and after: 1.5 × 10⁻² to 9.1 relative on
`beam_profile_kernel`'s five final states (the cross section takes the final state's formal amplitude, Γ_γ its Brune
one; with the group order of `TransformIn`, `level_parameters.md` §6); reported formal-basis photon widths with external
capture do not reproduce the input there. (3) The engine's S (about 10⁻¹⁰) and S′ (about 10⁻⁸) are not smooth in their
last digits, so a one-ulp change of E_c moves Brune amplitudes by 10⁻¹¹, more where a level's Brune denominator is
near 0 or negative; the 10⁻¹² expectation for Brune/Park was not met. (4) The step-1 harness ran `X` and `X/models` in
one work tree concurrently; one run of `shift_subpoint_angle` was lost here (fixed in the step-1b harness).

**Open (proposed default in brackets).**
- Sign rule: orientation of the Brune level instead of decision 10's "largest |γ| positive". [Keep the orientation
  rule; PI to confirm. The alternative needs the sign of every final state carried to the photon channels of all
  levels that decay to it, across blocks.]
- Formal-basis external capture differs from the Brune basis (Found 2). [Decide in step 2 or 4: compute γ^ext and the
  capture normalization from the final state's Brune (observed) amplitude in every basis.]
- Formal input written by earlier builds (`param.sav`, `--no-transform`) for out-of-order groups now means another
  R-matrix; the engine notes it. [Accept; no converter. A project that needs explicit B_c values gets a keyword with
  `io/Project`.]
- Brune/Park last digits depend on E_c order and the memo at 10⁻¹¹ (10⁻⁹ in ill-conditioned levels). [Accept until
  `channel_functions` (step 3) gives smooth S, S′.]

**Next.** Step 2 (observed widths from total amplitudes), with the external-capture question above.
## 2026-11-07 (`feature/upgrade-remove-qt`) — the Qt GUI removed (decision 3); parity checklist

Branch `feature/upgrade-remove-qt` from `feature/upgrade` 4bf7aa06, not merged. Commits: the parity checklist
(alone, before the removal), the build removal, the GUI-only engine entry points, packaging/CI, the Sphinx docs, the
skills, CLAUDE.md files and README, this entry with ROADMAP.

**State.** `docs/upgrade/qt_parity.md` lists what the Qt setup utility did at 4bf7aa06 (201 items in 15 sections:
.azr blocks and columns, keyword whitelist, round-trip guarantees, each tab, the THM workspace, network fetches,
environment hooks; 32 lines, one per `tests/gui` file); the Phase 9 web UI is checked against it. Removed: `gui/`
(124 files) and `tests/gui/` (33), 47,651 lines; `BUILD_GUI`, `USE_QWT`, `BUILD_MACOS_BUNDLE`, every Qt/Qwt
find_package, `FindQwt.cmake`, `GUI_BUILD`/`start_gui`/`gui_stub.cpp`, the Qt DLL copies, Qt/Qwt in CI (apt, brew,
MSYS2, windeployqt, linuxdeploy's Qt plugin), both Dockerfiles, `build_linux.sh`, `build_windows.sh`,
`build_crc.csh`; `build_macos.sh` (it built the Qt .app/.dmg) and `examples/run_gui.sh` deleted. Engine entry points
only the GUI called, each with no caller left: `AZURECalcMCMC::RequestStop/ClearStop/SetGUIResultsCallback/
SetGUIMessageCallback` (and the pointers only they fed), `NuclearPotentialManager` getters (7 functions,
`currentType_`), `AdaptiveIntegrationGrid::GetConfig`. Docs: new `user_guide/graphical_interface.rst`; `file_menu`,
`plot`, `gui_structure` pages deleted; tab-framed user-guide pages keep their content (it documents the file) under a
note. Skills, CLAUDE.md (two .azr parsers, not three; no GUI step), tests/CLAUDE.md, pyazr/CLAUDE.md, README,
repo_consistency (no tests/gui, no SegmentKeywords.h) follow.

**Kept, and why.** `--no-gui` is accepted and ignored (scripts, `tests/upgrade` and the goldens' runners pass it).
`ThmDistortion::GlobalEnds` (unused now; the THM Diagnostics page drew it): a physics function the Phase 9 Physics/THM
tab needs through the binding. The Windows icon and the AppImage icon moved to `packaging/icons/`. The rolling
pre-release job (dev only) is unchanged but for the Windows BUILD-INFO text; macOS CI now uploads a tar.gz of the bare
executables and the job keeps the previous .dmg. `docs/upgrade` history (`cnuc_split.md`, `skills_migration.md`,
`THM.md`, `MERGE_NOTES.md`) and the "files written by the Qt GUI" notes in skills are left as written.

**Evidence.** Condor (cluster 16179), modules python + gsl only (no qt, no qwt), fresh configure of 5f09a5ce: the
`--trace-expand` log has 0 find_package/find_library calls for Qt or Qwt, CMakeCache 0 Qt/Qwt entries, `ldd AZURE2`
0 Qt libraries; the 7 CMake warnings are all in `external/` and `numcmc/`, which this branch does not touch; build: 0
compiler warnings.
`tests/upgrade/check.sh` PASS; `ctest -LE slow` 107/107 (check_skills failed in that run because the skills were not
yet edited; 3/3 of check_skills and the consistency tests pass at a90474d4); `tests/run_tests.sh` 131/132, the one
failure `light_excitation` (a), see Found; `tools/check_skills.py` 0 errors. Docs: every :doc:, toctree entry, :ref:
and image of the 37 pages resolves (script; Sphinx is not installed). The engine's sources changed only by deletions
and the `--no-gui` help text, so no output is expected to move.

**Found.** (1) `tests/light_excitation` (a) is flaky on `feature/upgrade` itself: `intEC.dat` of two runs differs at
5 × 10⁻¹¹ relative (tolerance 10⁻¹²). Same host, 4bf7aa06 build vs this branch: 8 threads 3/3 vs 3/3 failures,
2 threads 4/10 vs 4/10. The EC integrals depend on thread scheduling at that level (determinism rule). (2) pyazr's
`AzrModel` cannot create a particle pair (skill azr-levels §1); without the GUI a new compound nucleus is set up by
writing `<levels>` lines. (3) From the parity review: the Pairs table hid the light-E_x column (column 15) and Ctrl+D
help indices predate two tabs; recorded in qt_parity.md, not fixed (Qt stays on dev).

**Open (proposed default in brackets).**
- A self-contained macOS package of the CLI. [Leave the tar.gz until Phase 9's wheels; no .dmg.]
- If this merges into `dev` before the upgrade is released, the dev pre-release stops refreshing the macOS .dmg and
  the GUI downloads. [The upgrade does not merge into dev before Phase 9; the director confirms.]
- `light_excitation` flake. [A fix/ branch: fixed reduction order in the EC integrals or a tolerance with a written
  reason, after finding which loop is scheduled.]
- `ThmDistortion::GlobalEnds` in the binding. [Phase 9, with the Physics tab.]

**Next.** The director merges; Phase 9 checks the web UI against `qt_parity.md`.
## 2026-11-07 (director) — COULCC in C++, bit-identical to the Fortran (pycoulomb, local branch `coulcc-cpp`)

COULCC (Thompson and Barnett, Comput. Phys. Commun. 36 (1985) 363) translated to structured C++17 in Jakub's
pycoulomb clone, branch `coulcc-cpp` (from 9fee957; e378a33 … 37f0f9b; not pushed). One function per Fortran routine,
COMMON blocks as members of a per-call object, no mutable global state. Against gfortran (-O2/-O3, x86-64, GCC 11.5)
on 12 814 cases (real, complex ρ/η/λ, Bessel kinds, all MODE1 × KFN): 618 988 values, every IFAIL, RERR and
diagnostic counter bit-identical. Bit identity needed gfortran's complex semantics (mixed real/complex promotion,
Smith division) and REAL*16 as `__float128`; elsewhere a double-double, equal to rounding. The literal goto port it
replaces had 39 IFAIL mismatches, sign errors of J_ν, Y_ν at negative real ρ, and a data race.

**Algorithm limits found (same in the Fortran; to report to I. J. Thompson).** (1) 2λ non-integer, ρ < ½: G from
the reflection formula divides by sin χ, which cancels (−2.2 × 10⁻¹⁴ at η = 5, λ = 2.7, ρ = 0.1), G 1.3 % wrong with
IFAIL 0; the Wronskian cannot see it. (2) 2λ integer, ρ < ½: the logarithmic origin series loses accuracy to
8 × 10⁻¹¹, not reported in RERR. (3) Re λ_min ≤ −1 (outside the stated domain): Wronskian wrong by up to 54 % with
IFAIL 0. Real-argument accuracy against mpmath otherwise ≤ 7 × 10⁻¹² (CF2), Bessel ≤ 3.4 × 10⁻¹⁴.

**Next.** Merge `coulcc-cpp` with `azure2-prereqs` and add COULCC as a libcoulomb method (one adaptor; complex
arguments through `coulomb_complex()`); then it comes into `physics/coulomb/` with the rest. Not used as AZURE2's
default (`coulsm` stays). Pushing both branches is Jakub's.

## 2026-11-07 (`feature/upgrade-1-cnuc-split`) — split `CNuc`, step 1: the engine calls `level_parameters`

Branch `feature/upgrade-1-cnuc-split` from `feature/upgrade` 93f8cbc5, not merged. Commits: the plan
(`cnuc_split.md`, alone), the physics diagnostics and eigensolver option, the adapters, then one commit per moved
function (`CalcBoundaryConditions`, `TransformIn`, `TransformOut`'s formal iteration), the docs, this entry.

**State.** `cnuc_split.md` classifies all CNuc members (nucleus data, physics functions, parsing, reporting,
parameter vector), states the interface of each target part, and orders the moves (steps 1, 1b, 2–7). Step 1 is done:
B_c, observed → Brune/Park amplitudes, the external-capture subtraction, Brune → formal (GSL's eigensolver passed as
an option) and the formal → Brune iteration of `parameters.out` are `level_parameters`' code; the old bodies are
deleted (`CNuc.cpp` 3628 → 3311 lines). `CNuc` keeps the wording of its warnings, Park's J on the level and N_f^½.

**Evidence.** Engine before (93f8cbc5) vs after, all 179 `.azr` under `tests/` (T1–T3 included) × Brune, `--no-brune`,
`--use-park`: 1503 output files byte-identical; logs identical except `screening_6li` (`--no-brune`), whose warning now
gives the last energy step (118.45, 40.9187, 85.4476 MeV) instead of 0. With the own eigensolver the formal basis
would change 53 files (698 sign flips in `param.par`, observables ≤ 2.2e-11), so it stays on GSL until step 1b.
GUI off, Condor: no warnings; `check.sh` PASS; `ctest -LE slow` 106/107 (`adaptive_grid_reference`, decision 8);
`run_tests.sh` 132/132 (GUI on earlier: 131/132, `shared_point_tables` (d), see Found). T2 bench same host: fit200
79.32 vs 80.10 s (1 thread), 36.16 vs 36.28 s (8); calc 1.17
vs 1.19 s, 0.42 vs 0.41 s; no slowdown beyond one 10 ms tick.

**Found.** (1) The Coulomb memo matches energies within 1e-12 MeV, so the engine's numbers depend on the order of
its Coulomb-function calls at 1e-13; with the three orders of E_c (E − e − S, E − (S + e), E − S − e) this decided
bit-identity (`cnuc_split.md` §5). The adapter reproduces each caller's order and call sequence. (2) A group whose first
listed level is not an R-matrix level had B_c uninitialized; it now takes the first R-matrix level (none in tests/).
(3) `tests/shared_point_tables` (d) peak-RSS ratio sits at 2.9–3.0 against its limit 3.0 for both engines: it failed
once in `run_tests.sh` (GUI build) and once in 3 repeats of the new GUI-off engine (316 > 315 MB), the old engine
passed 3/3 at 305–311 MB; same RSS ranges, so a flaky margin, not this change.

**Open (proposed default in brackets).**
- The adapter options `ChannelEnergyOrder` and `ChannelQuantities` exist only to reproduce the memo history. [Keep
  until step 3; or drop them in step 1b, whose formal re-baseline could absorb a 1e-13 change of the other bases too —
  director's choice.]
- The FACTOR 3.0 of `shared_point_tables` (d) has no margin. [Measure the ratio on 10 runs and set the limit with a
  written reason, separate commit.]
- `CNuc::TransformOut`'s observed widths, `GetTransformParams`, `PrintTransformParams`, `CheckRadiativeWidths` still
  each convert γ → Γ. [Step 2: one observed-width function with γ^ext.]

**Next.** Step 1b (decision 10 on, own eigensolvers, GSL adapters deleted, formal goldens re-baselined), then step 2.
## 2026-11-07 (director) — pycoulomb prerequisites done (in Jakub's repo, local branch)

Decision 8 (a) merged (409fef62): `ctest -LE slow` 107/107 with GUI off. pycoulomb, local branch `azure2-prereqs` of
`~/AI-R-scratch/pycoulomb-src/pycoulomb` (from 9fee957; commits 7e3c834 … 7adef98; not pushed, not tagged):
build products removed; `DESIGN.md`; `libcoulomb/`, one CMake target `pycoulomb::coulomb` without Python, the
vendorable set; pure, thread-safe API `coulomb(l, η, ρ, Method)` → {F, F′, G, G′, σ_l, status}, `coulomb_range`,
`coulomb_phase`, `whittaker`, `whittaker_scaled`; methods `coul` (AZURE2 before 2025-09), `coulsm` (AZURE2 now,
default), `coulfg`, `gsl`. Parity with AZURE2 on 7150 points: 0 bit differences (`coulsm` vs `external/coul`, gsl vs
`CoulFunc`, Whittaker vs `WhitFunc`/`CoulombTable`). Threads: 0 bit mismatches at 8 threads; TSan clean.

**Accuracy against mpmath (40 digits; 1296 points, l 0–8, η 0–20, ρ 10⁻³–50).** `coul`/`coulsm`: ≤ 10⁻⁹ for η > 0
(their built-in precision 10⁻¹⁰). `coulsm` at η = 0: F_l from unstable upward recurrence is wrong for l ≥ 2 at
small ρ (G exact), but P, S, φ depend on F only through F/G, absolute error ≤ 1.4 × 10⁻¹⁶ on the grid; this is the
open neutral-F item of the T1 session. `coulfg`: JWKB inside the turning point, wrong where G > 300. `gsl`: errors to
2.5 × 10⁻² inside the turning point for η ≥ 2. All four agree to 8.3 × 10⁻¹⁰ in the oscillating region.

**For Jakub.** Add the LLNL BSD-3 COPYRIGHT text next to `coulfg2.c` (its header refers to a missing file); choose
the repository licence; push and tag. Whether to fix `coulsm` at η = 0 [keep: AZURE2 stays bit-identical; P, S, φ
unaffected]. **Next.** After the CNuc split: take `libcoulomb/`'s sources into `physics/coulomb/` (not `external/`:
 Jakub, 2026-11-07, the Coulomb functions are physics) and build
`channel_functions` on it (replacing `external/coul`, `CoulFunc`, `CoulFuncCache`, `WhitFunc`), bit-identical with
`coulsm`; COULCC added as a method when `coulcc-cpp` is done. The host must keep GSL's error handler off (`GslQuiet`).

## 2026-11-07 — `fix/adaptive-grid-criterion`: decision 8 (a) in `adaptive_grid_reference`

Branch `fix/adaptive-grid-criterion` from `feature/upgrade` e5dd8ae7, not merged. Commits: the test (47c119e4), this
entry and the notes in `ROADMAP.md` and `sanitizer_findings.md` §4.6. No engine change; the grid is thm's, unchanged.

**Criterion.** The two cases "eV/keV level crossing the window edge" check the total variation of the trapezoid
error over the scan minus |net change| (relative to the area πΓ/2) instead of its swing; the accuracy check (2 × 10⁻⁴)
and the five other cases are unchanged. A one-way drift, here the tail bias (+1.69 × 10⁻⁴) leaving the window with
the level, costs nothing; each back-and-forth excursion of size d costs 2d. Tolerance, both cases:
√3/(2π ppw²) + 10⁻⁴ = 2.10 × 10⁻⁴ at 50 points per width.
- √3/(2π ppw²) = 1.10 × 10⁻⁴: a lattice of pitch h = Γ/ppw cut by a fixed edge has the Euler–Maclaurin end term
  (h²/12) f′(edge); as the peak crosses, f′(edge) runs 0 → −max → +max → 0 (max |f′| = 3√3/(4Γ) at Γ/(2√3) from
  E_R), a variation of 4 max|f′| with no net change. Measured 1.102 × 10⁻⁴ (eV level, scan of 0.00125 Γ steps across
  the edge); in the keV case the cells from 20 Γ below E_R upward give 1.102 × 10⁻⁴. It changes the slope of the
  in-window area in E_R by at most 2/(3 ppw²) = 3 × 10⁻⁴ of itself, so it cannot reverse a gradient.
- 10⁻⁴: the smoothness tolerance of the in-window cases, here for the sum of all other reversals. keV case:
  7.3 × 10⁻⁵ from the cells below the 20-width core (graded tail meeting the fixed 5 keV coarse lattice), a ripple of
  6.8 × 10⁻⁶ peak to peak; the in-window keV case has the same ripple (swing 2.4 × 10⁻⁵).
A sampled scan shows at most the excess of the continuous curve, so the bound holds for any scan step. Values
(2000 / 200 000 steps): eV 1.96 × 10⁻⁵ / 1.14 × 10⁻⁴ (its 1 eV step, 1.25 Γ, samples the edge term coarsely; net
change −1.685 × 10⁻⁴), keV 1.79 × 10⁻⁴ / 1.83 × 10⁻⁴ (net 4 × 10⁻⁷).

**Detection.** The test built against the grid before the graded tails (dev 07581f03, through a scratch shim that
feeds it the resonance list) reproduces the 17O+α failure: in-window eV case max error 7.85 × 10⁻², swing
3.53 × 10⁻². The new crossing check fails on it: variation beyond the net change 2.87 × 10⁻¹ against 2.10 × 10⁻⁴
(net −7.49 × 10⁻², swing 7.87 × 10⁻²); 7 checks fail there, none on the current grid.

**Verified.** Condor 16150/16152, GCC 11.5, `-DBUILD_GUI=OFF`, pyazr on: no warnings in the
rebuild; `adaptive_grid_reference` all checks pass; `ctest -LE slow -j8` 107/107 pass; `tests/upgrade/check.sh` PASS.

## 2026-11-07 (`feature/upgrade-1-level-parameters`) — `Convention` record and `level_parameters` (Phase 1, code part)

Branch `feature/upgrade-1-level-parameters` from `feature/upgrade` e54e6fcc, not merged. Commits: the record
(fbda8bf5), 3j/6j and eigensolvers (becc8bfe), `level_parameters` with its physics tests (45443582), the `CNuc`
adapter and T1–T3 test (7bd0e1e9), the docs (89061935), this entry.

**State.** `physics/` has the `Convention` record (`Convention.h`; the System carries one, `basis()` reads it) and
`level_parameters` (`LevelParameters.h`, `level_parameters.md`): observed ↔ Brune/Park, B_c (first listed or lowest
level), external-capture subtraction, Brune ↔ formal (own Cholesky + Jacobi; the engine's fixed-point inverse),
Brune ↔ Park, and the convention matrix T on the γ columns (LS orders, jj by 6j, iˡ as (−1)^⌊l/2⌋, E/M signs), with
`recoupling_gaps` reporting missing channels. Channel functions are an input (`ChannelFunctions`); the adapter
`src/CNucLevelParameters.cpp` supplies them from `ChannelFunc`/`WhitFunc`. Nothing in the engine calls the new code.
Checked on Condor (GCC 11.5, GUI on, pyazr on): no warnings; `tests/upgrade/check.sh` PASS (goldens unchanged);
`ctest -LE slow -j8` 137/138, the one failure `adaptive_grid_reference` (decision 8). Against `CNuc` on T1–T3: Brune
and Park amplitudes, Park's J and B_c bit-identical; formal E ≤ 6.6e-16, formal γ ≤ 1.0e-14 (up to the level's sign,
which differs from GSL's on 3 of 19 levels); formal → Brune vs `TransformOut` ≤ 6.9e-14; six conventions give the
same formal parameters to 1.5e-15 (parameter level; the engine reads no convention). Physics tests: 10 512 recoupling
coefficients vs m-scheme to 8.9e-16; formal vs Brune level matrix identity 4.4e-15. Six mutations each fail at least
one test (`level_parameters.md` §7).

**Found.** T3's 2⁺ block lists only s = 2 of the two channel spins J = 2 allows for p + ⁷Li and n + ⁷Be at l = 1, so
T3 has no jj form for that block (reported, as intended). In the formal basis `CNuc::TransformIn` evaluates external
capture after earlier groups were made formal, so a final state in an earlier multi-level group would enter with
formal E, γ; on T1–T3 this changes nothing (T1's 14 non-zero external amplitudes go to the single-level ground-state
group). It does not explain the 54 % T1 formal-basis capture difference of `tests/upgrade/README.md`.

**Open (proposed default in brackets).**
- The two proposals (B_c from the lowest level; largest |γ| positive) are options, off. [PI decides; keep off until
  the engine switches to `level_parameters`, then turn both on together with a golden re-baseline of formal-basis
  outputs only.]
- iˡ is represented by the real signs (−1)^⌊l/2⌋ (both i^l and i^−l); no code with that convention was compared.
  [Compare with one external code's amplitudes (SAMMY or FRESCO output for one level) before a project uses it.]
- jj means j = l + j_x with x = light by default (`spin_order` picks the particle); the final coupling J = j + j_y is
  fixed. [Enough for nucleon + target; add a field only when a source needs (j_y, j).]
- External-capture order in the formal basis (above). [Decide whether to compute every γ^ext from Brune final states
  when the engine switches to `level_parameters`; not changed here.]
- E2 static moments scale the input as v × (100 (αħc)^½), the engine as (v × 100) × (αħc)^½: possibly 1 ulp apart.
  [None in T1–T3; accept.]
- No `.azr` keyword, binding or GUI for the record. [Add with `io/Project` and the binding of `level_parameters`
  once `channel_functions` exists.]

**Next.** Split `CNuc`: the engine's TransformIn/TransformOut/CalcBoundaryConditions call `level_parameters`
(bit-identical for Brune/Park shown here; formal to 1e-14, so formal-basis goldens move at that level).

## 2026-11-07 (director) — Phase 1 merges; ignored headers recovered

**Merged into `feature/upgrade`.** `feature/upgrade-1-conventions` (655c5634), `feature/upgrade-1-physics-system`
(f31d2027), `docs/upgrade-1-fix-contradictions` (7e3e5643); ARCHITECTURE §3 amended (B_c per Jπ block, 1af53f98).

**Defect found and fixed (9da89268).** The `.gitignore` pattern `azure2` also matched `physics/include/azure2/`, so
`System.h` and `JpiBlock.h` were never committed and a fresh checkout of f31d2027 did not build. The pattern is now
`/azure2`; the headers were restored unchanged from the session that wrote them. Rule for every branch: run
`git status --short --ignored` before hand-back.

**Verification.** `feature/upgrade` 7e3e5643, GCC 11.5 with GUI, Condor 16128: no compiler warnings;
`tests/upgrade/check.sh` ok; `ctest -LE slow` 133/134, the one failure `adaptive_grid_reference` (decision 8).
Sphinx not built (not installed on CRC): run `sphinx-build -W docs/source` before pushing; the cross-page `:eq:`
references from `level_matrix.rst` into `conventions.rst` are untested.

**Next (Phase 1, without pycoulomb).** `Convention` record and `level_parameters` (formal ↔ Brune/Park, B_c from
the lowest level, sign rule for formal levels), checked by the convention-invariance test on T1–T3; then the `CNuc`
split. `channel_functions` waits for the pycoulomb prerequisites (decision 6).

## 2026-11-07 (`feature/upgrade-1-physics-system`) — `physics/` target, `System` and `JpiBlock` (Phase 1, first step)

Branch `feature/upgrade-1-physics-system` from `feature/upgrade` f8b7bf22, not merged. Commits: the physics target,
the CNuc adapter with its test, the docs.

**State.** `physics/` is its own CMake target `azure2_physics` (headers `physics/include/azure2/physics/`, sources
`physics/src/`, tests `physics/tests/`); it links nothing and is added before the global GSL/OpenMP/Minuit2 include
directories. `System` (nuclei, partitions α, channels c = (α, l, s) with a_c, Jπ blocks, the level basis) and
`JpiBlock` (2J, π, nucleus, channel index list, B_c, E[λ], γ[λ n_c + k]) are immutable values with 0-based indices
and doubled spins; the choices and their reasons are in the header comments and `physics_system.md`. The adapter
`BuildSystemFromCNuc` (`include/CNucSystem.h`, `src/CNucSystem.cpp`, in AZURE2Core) builds a System from a CNuc
after `Initialize`; nothing in the engine calls it. Checked on Condor (GCC 11.5, GUI off, pyazr on):
`physics_system_from_cnuc` compares 2540 numbers of T1, T2, T3 in the Brune, Park and formal bases with the CNuc
accessors, all equal, and the structure with the `<levels>` lines (5 mutations of the adapter each fail it);
`physics_system_unit` and `physics_standalone` pass (an `#include "CNuc.h"` in physics/ fails the latter);
`repo_consistency` now also refuses non-standard includes in physics/ (a `<gsl/gsl_math.h>` fails it);
`tests/upgrade/check.sh` PASS (golden T1, T3 largest difference 0.0, T2 calc_start/calc_truth ok, recovery ok);
`ctest -LE slow -j8` 102/103, the one failure `adaptive_grid_reference`, as before.

**Open (proposed default in brackets).**
- B_c is stored per block because the engine sets it per J-group from the first level; ARCHITECTURE §3 lists it with
  the channel. [Keep it in the block; `level_parameters` may later want B_c = S_c(E_λ) per level, which is a block
  quantity too.]
- The blocks hold the parameters of the project's basis (Brune by default), not formal ones; `System::basis()` says
  which. [Keep until `level_parameters` exists; it then fills the blocks in the internal convention and the basis
  tag moves to the `Convention` record.]
- External-capture final states that are not R-matrix levels, hybrid (α, l, j) channels and the THM binding energy
  are not in the System. [Add with the items that need them: external capture with `channel_functions`, hybrid
  channels with pycoulomb, B_THM with the THM vertex.]
- The adapter reads `ALevel::GetE/GetGamma` (the loaded parameters), not the moving fit values. [Right for setup
  state; per-evaluation E, γ come from `level_parameters`.]

**Next.** `channel_functions(α, l, E)` on the System's partitions; then split `CNuc` (ROADMAP Phase 1).
## 2026-11-07 (`feature/upgrade-1-conventions`) — conventions page, documentation part (ROADMAP Phase 1)

Branch `feature/upgrade-1-conventions` from `feature/upgrade` f8b7bf22, merged 655c5634. Commits: 754ac7d7 (test
`tests/coulomb_nuclear_sign`), 8faa430d (`docs/source/theory/conventions.rst` + Theory toctree), and this
HANDOFF/ROADMAP commit. No C++ changed.

**State.** The page states the present engine's conventions, each with file:line: s = j1 + j2 light first (LT
III.2a), J = s + l in the spin-amplitude route, Condon–Shortley CG/Racah; no i^l on particle channels (only the
E external-capture amplitude carries −i i^(l+L−lf), M1 +i); φ = arctan(F/G), ω = σ_l − σ_0, U = ΩΩ[1 + 2iP^½γAγP^½],
T = e^{2iω} − U; B_c = S_c(E of the first listed level) in the formal basis; sign of γ = sign of the input Γ/ANC/
amplitude (Brune, Park; photon total width with EC; static-moment sign flip for half-integer J); photon widths
never in the level matrix except --use-rmc; Blatt–Biedenharn z = (−1)^(s'−s) Z Z / 4 (photon: (−1)^(1+s−jf));
elastic Coulomb term with C_α = +(4π)^−½ η csc² e^{−2iη ln sin} and Re[i C T*]; identical-particle factors; the
planned Convention record (marked planned). `tests/coulomb_nuclear_sign/check.sh`: p+12C hard sphere vs an
mpmath partial-wave sum, 30 points, A- and R-matrix, Brune and --no-brune, max rel diff 8.4e-11 (Condor 16122);
the interference reversed moves the truth by 7–74 %. Sphinx not built (not installed on CRC; labels and :doc:
targets checked by script).

**Existing docs that contradict the code (not changed here).**
- `docs/source/developer/level_matrix.rst`: "L_c = L_0(c) is the shift function" (code: L⁰ = S − B + iP); U
  written without Ω_cΩ_c'; T_cc = e^{2iω} − e^{iω}e^{iφ}(1+U) and T_cc' = −e^{iω_c}e^{iφ_c'}U (code: Ω_cΩ_c' with
  e^{−iφ} on both channels); "ω_c = σ_l" (code: σ_l − σ_0); [1−RL] = δ − iRP (code: 1 − RL⁰ in the standard
  formalism); interference "(i/√π) g f_C* T P_l + c.c." (code: Re[(i/√π)(π/k²) g C_α T* P_l], opposite sign and
  half the size if f_C means C_α).
- `docs/source/reference/data_formats.rst` (phase-shift data): "−90° < θ_lab < 90°" names the wrong quantity; the
  phase shift is in (−90°, 90°], and in [0°, 180°) for identical pairs (`GenMatrixFunc.cpp:555-564`).

**Open (proposed default in brackets).**
- Sign of internal relative to external capture: fixed by the code's formulas, checked by no independent
  calculation. [An m-scheme or potential-model truth for one resonance + EC, as `tests/ec_e1m1_phase`.]
- The orbital–spin order is pinned only relative to the Blatt–Biedenharn route (polarization Gate 1). [Accept;
  the invariance test of the record will cover it.]
- LT and DB were not re-read (no PDF reader on CRC); the page says which comparisons rest on the code's citations.
  [Re-read LT VII–VIII and DB §2–3 on a machine with a PDF reader and add equation numbers.]
- Formal basis: B_c = S_c(first listed level) makes the formal parameters depend on level order, and GSL's
  eigenvector sign makes a formal level's overall sign arbitrary. Observables are unaffected. [Fix both in
  `level_parameters`: B_c from the lowest level, sign so that the largest |γ| of a level is positive.]
- Fix the two contradicting pages. [Separate `docs/` commit, by the director's choice.]
  Done by 57bcd03d on `docs/upgrade-1-fix-contradictions` (not merged): all listed contradictions confirmed
  and fixed; also the Reich–Moore sign on the A⁻¹ diagonal (`+= iγ²`, code `−= iγ²`).

**Next.** Convention record and `level_parameters` (code part of the ROADMAP item).

## 2026-11-07 (`fix/edep-ownership-tsan`) — copy-on-write ownership TSan can see (sanitizer_findings.md §4.2)

Branch `fix/edep-ownership-tsan` from `feature/upgrade` 5b021003, not merged. Commits: c556f052 (the fix), and the
findings/HANDOFF commit.

**State.** `EPoint::edep_` and `legendreP_` are `SharedTable<T>` (`include/SharedTable.h`): a shared_ptr plus an
atomic holder count (relaxed increment, release decrement, acquire load in `SoleOwner`); the use_count() + fence test
is gone. Checked on Condor: GCC 11.5 build with GUI and pyazr; `tests/upgrade/check.sh` PASS at 8 threads (golden
outputs unchanged); `tests/run_tests.sh` 131/131; `ctest -j8` 135/140, the five failures inherited from origin/thm.
clang TSan on a temporary branch with this fix and `fix/fit-iteration-counter` (`feature/upgrade` + fecf4924 +
c556f052): 54 steps, no report file (before: 25 files, 128 reports, 8 AZURE2 sites). Control: the acquire load made
relaxed brings back 75 reports at `EDependentTables::Clear` (checks section).

**Open (proposed default in brackets).**
- Under TSan, `shared_point_tables` (d) now fails without a report: peak RSS 1162/1164 MB against the limit 1149 MB
  (1133 MB before), from TSan's per-atomic synchronization records; plain builds are unaffected and the ratio is
  near the limit in plain builds too (one failure in six runs, on the branch that does not touch the tables).
  [Run that check under TSan with `FACTOR=3.5` from `sanitize_upgrade.sh`, with the reason in the script; the
  plain-build limit stays 3.0. Not done here: it changes a test limit, which is the director's call.]
- `pyazr_upgrade_thread_independence` failed once in a full `ctest -j8` on this branch (T3 fit `param.fit`, 2.3e-5
  relative, 1 vs 4 threads; parameters equal); four isolated repeats on this branch, on `fix/fit-iteration-counter`
  and on `feature/upgrade` all pass. [Cause not found; not reproduced. Rerun the test
  a few times under load before the next merge.]

**Next.** Merge with `fix/fit-iteration-counter` (both add a HANDOFF entry here: keep both).
## 2026-11-07 (`fix/fit-iteration-counter`) — fit iteration counter atomic (sanitizer_findings.md §4.1)

Branch `fix/fit-iteration-counter` from `feature/upgrade` 5b021003, not merged. Commits: fecf4924 (the fix), and the
findings/HANDOFF commit.

**State.** `EData::iterations_` is `std::atomic<int>`; `EData::Iterate()` returns the number of this evaluation
(`fetch_add`), `AZURECalc::operator()` uses it, and the progress line is written under a mutex. Checked on Condor:
GCC 11.5 build with GUI and pyazr; `tests/upgrade/check.sh` PASS at 8 threads (golden outputs unchanged);
`tests/run_tests.sh` 131/131; `ctest -j8` 135/140, the five failures inherited from origin/thm (fourth session,
below). A T1 fit with Minuit2's numerical gradient prints the same 73 progress lines as `feature/upgrade` at one
thread, and the numbers 10, 20, …, 730 once each at eight. clang TSan on a temporary branch with this fix and
`fix/edep-ownership-tsan`: 54 steps, no report (before: 25 report files; 48 reports at this finding's three sites).

**Next.** Merge with `fix/edep-ownership-tsan` (both add a HANDOFF entry here: keep both).
## 2026-11-07 — `fix/adaptive-grid-edge-swing`: adaptive_grid_reference, analysis only (decision needed)

**What fails.** "eV level crossing the window edge" (Γ = 0.8 eV, 10 keV window, 50 points, level moved from
0.9 to 1.1 of the window in 2000 steps): swing of the error 1.685 × 10⁻⁴ of the area πΓ/2, tolerance 10⁻⁴. The
check is from dev (d216517); the grid is thm's (02ae97b), kept at the merge of dev 2d97b1f0, where MERGE_NOTES.md
already lists this failure as a decision.

**Cause.** thm's grid has geometric tails whose step is r times the distance from the level,
r = TailFactor() = 2.5/50 × √(20/2) = 0.158 for a 20-width core. A chord over a 1/x² wing over-counts by about
r²/2, so the area comes out high by r²Γ/(2π · core) ≈ 2 × 10⁻⁴; measured +1.69 × 10⁻⁴ inside the window, constant
to 2.2 × 10⁻⁷ as the level moves (the "eV level, 10 keV window" case). In the crossing scan the error starts at
that bias and ends at ~0 with the level outside, so the swing equals the net change: swing 1.685 × 10⁻⁴, net
change −1.685 × 10⁻⁴, total variation beyond |net| 2.0 × 10⁻⁵, largest jump between neighbouring positions
(1 eV) 9.4 × 10⁻⁶. For any grid whose bias b passes the case's accuracy check (|b| < 2 × 10⁻⁴), this swing is at
least |b|: here the swing check is a 10⁻⁴ accuracy bound, not the smoothness check it is described as. It is not
an oscillation (the 17O+a failure it was written for had 3.5 × 10⁻² with a 10 eV period). dev's grid passes
because its bias is 7.2 × 10⁻⁵.

**Change.** a test commit that only prints the net change and the variation beyond it for every case; no check,
tolerance or case changed. The test still fails.

**Options.** (a) Keep thm's grid; for the crossing cases check the variation beyond the net change (2.0 × 10⁻⁵
eV case, 1.8 × 10⁻⁴ keV case at present, so its tolerance needs its own argument) and leave the bias to the
accuracy check. (b) Lower tailRatio 2.5 → ~1.9 (b ∝ r², ~30 % more tail nodes): every fold moves by ~10⁻⁴, so the
regression pins and the upgrade goldens (1e-6) move; not done. (c) dev's grid: moves THM outputs. [Proposed
default: (a), decided by the PI, since it changes a reference test's criterion.]
## 2026-11-07 — `fix/c12c12-kernel-pin`: pair_identity and the THM examples' chi-squared

**Cause.** `pair_identity` pinned the engine total of `examples/c12c12_tumino2018` at 112.0376 (data 111.613 +
priors 0.4248), the value of the merge of dev 64ed575. The merge of dev 2d97b1f0 brought dev 7d24c88: Gaussian
resolution kernels to 5σ instead of 3σ. MERGE_NOTES.md measured the move (data 111.613 → 111.661, direct segments
5-8 only; convolutionRange = 3 restores thm's files byte for byte) but the GUI test was not built and the
documented totals were not updated.

**Why the new value is right.** A Gaussian cut at ±kσ and renormalised has variance
σ²[1 − 2kφ(k)/(2Φ(k) − 1)]: 0.9733 σ² at k = 3 (the fold used a resolution 1.3 % narrower than stated),
1 − 1.5 × 10⁻⁵ at k = 5.

**Change.** f0981b53: the pin is 112.0856 (tolerance 0.01 unchanged). bb1bd4f0: examples.rst,
thm_implementation.rst and the azure2-thm calibration table carry the measured totals (build-upgrade, 4 threads):
c12c12 112.086 (111.661; THM segments 1-4 61.02, unchanged), li7_tumino2006 448.620 (439.149), li6_pizzone2011
158.967 (157.041). The other five THM examples give their documented values.
## 2026-11-07 — `fix/mcmc-shift-prior-centre`: mcmc_priors and the MCMC tab's shift-prior centre

**Cause.** `mcmc_priors` (thm, 2026-10-09) expected the energy-shift prior of segment 1 at its start value
(0.004 ± 0.002 MeV) and the sqrt(E) coefficient's at 0.003 ± 0.001. dev 0dd28e5e (2026-10-25), merged later,
centres both on no shift: the fit's penalty (`ESegment::AddShiftPenalty`) and the sampler's automatic prior
(`AZURECalcMCMC::BuildAutoPriors`). The tab follows it; the test was out of date. Reading the code showed a second
difference: the engine centres the shift prior on `GetNominalEnergyShift()`, which a
`segment_N_energy_shift prior_centre c` row of `<parameterSettings>` moves; the tab used 0 regardless.

**Change.** dc4cbaa6: the test expects 0 ± 0.002 MeV and 0 ± 0.001. 2a137b4c: the tab takes
`GetNominalEnergyShift()` (as it takes `GetNominalNorm()`); test section 2b, with the row 0.004, expects 0.004 ±
0.002 MeV for the shift and 0 ± 0.001 for the sqrt(E) term. The old line gives 0 there, so 2b fails on it.

**Overlap.** None with the `EData`/`AZURECalc` counter or `EPoint` ownership work: only `gui/src/MCMCTab.cpp`.

**Open.** `ESegment.h` describes `GetNominalEnergyShiftSqrt()` as "the penalty's reference", but the penalty uses
0. [Correct the comment on the next branch that touches `ESegment`.]
## 2026-11-07 — `fix/gui-pair-columns`: thm_opt_in and light_excitation_gui after the merge of dev 3264896d

**Cause.** The merge of dev 3264896d into thm (made without Qt; MERGE_NOTES.md) added dev's light-particle
excitation e1 to `PairsData` and kept thm's THM columns: pairs model column 15 is e1, 16 the THM binding energy
B (`PairsModel::BindingEnergyColumn`, 15 before), 17 the compound nucleus; the GUI always writes B and the width
flag (columns 32, 33) and e1 as column 34 when it is nonzero; NucLine tells 1, 2 and 3 trailing fields apart. Two
GUI tests were never built against that layout. `thm_opt_in` wrote B into column 15, so it set e1 = 2.2246 MeV and
the project rightly opened with THM off. `light_excitation_gui` (from dev) expected 16 fields and 31/32 columns.
The code is right in both cases; the tests were out of date.

**Change.** `thm_opt_in` uses `BindingEnergyColumn`. `light_excitation_gui` checks SIZE 17, the three named
columns, 33 columns written for a pair without e1 and 34 with it, a 31-column dev file written with 33, and (new)
the 33/34-column output read back with e1 on pair 2 only. `reference/data_formats.rst` now states the
trailing-field rule (it described dev's 32nd column only).

**Verified.** The four `fix/` branches of this session merged on 5b021003 (scratch, detached), GUI on, GCC 11.5,
Condor 16112/16113: `ctest -LE slow` 130 of 131 pass, only `adaptive_grid_reference` fails; `tests/upgrade/check.sh`
PASS (no engine code changed). The failing tests of this branch failed on origin/thm c6ddb741 (plain GUI build).

**Open.** The GUI turns a 31-column classic file into 33 columns on save; the engine and pyazr read both. Whether
the THM columns move behind e1 when thm goes into dev is the decision MERGE_NOTES.md already lists. [Keep the
present layout until that merge.]

- With `build.sh` the pyazr module lands in `BUILD/pyazr-out/` (`PYAZR_OUTPUT_DIR`), but ctest imports
  `pyazr` from the source tree: 31 pyazr tests print `skip` and `pyazr_angular_coefficients` fails at import (it
  has no skip guard). The counts above are with the module copied into the scratch tree's `pyazr/`. [Have ctest
  put `PYAZR_OUTPUT_DIR` on `PYTHONPATH`, and give `angular_coefficients_test.py` the skip guard; own `fix/`
  branch.]

## 2026-11-07 (fourth session) — Phase 0 merged; failures inherited from origin/thm

**State.** `feature/upgrade` carries all Phase 0 branches: rules and skills (09379840), T1/T3 (a80352e6), T2
(356c6770), sanitizers (73290ee2). `tests/upgrade/check.sh` runs T1, T3 and the two T2 calculations against the
golden outputs (e4cadf51); Condor 16086 at 8 threads: PASS. Both race fixes are on `feature/upgrade` (ROADMAP item
ticked).

**Five ctest failures come from origin/thm, not from the upgrade.** A plain build of origin/thm c6ddb741 (GCC 11.5,
Qt GUI on) fails the same five tests as `thm-review` 1174fdb2:

| test | failing check |
|---|---|
| `thm_opt_in` | "binding energy: THM on" |
| `pair_identity` | 112.0856 against the expected 112.0376 |
| `mcmc_priors` | energy shift 0.004 ± 0.002 MeV and the sqrt(E) coefficient prior 0.003 ± 0.001 |
| `light_excitation_gui` | `PairsData::SIZE` 16 and the 31/32-column pair lines |
| `adaptive_grid_reference` | "eV level crossing the window edge": swing 1.685 × 10⁻⁴ of the area, tolerance 10⁻⁴ |

The log is in the scratch build directory of the analysis machine; the commands are `ctest -j8` on a GUI build.
None of the five touches the THM form-factor GUI added on `thm-review` (7c543690). Each needs a `fix/` branch.

**Next.** `fix/fit-iteration-counter` and `fix/edep-ownership-tsan` (TSan green); `fix/` branches for the five
inherited failures; then Phase 1.

**Later the same session.** Four of the five inherited failures fixed and merged (05e08e06, 1d35e1d4, 22214e14;
analysis of `adaptive_grid_reference` in c51f973e, decision 8). `fix/pyazr-tests-engine-dir`: a build with
`PYAZR_OUTPUT_DIR` outside `pyazr/` (as made by the scratch build helper and the sanitizer builds) skipped 31 pyazr
tests; ctest now sets `PYAZR_ENGINE_DIR`, which the package adds to its search path. `ctest -R pyazr -LE slow`:
74/74 passed, none skipped (Condor 16114). The four tests outside `tests/pyazr` that import pyazr
(`beta_delayed_alpha`, `consistency`, `electron_screening`, `hybrid_potential`) do not get the variable yet.

**TSan fixes merged.** `fix/fit-iteration-counter` (a2392e40) and `fix/edep-ownership-tsan` (f8b7bf22): clang TSan on
both, 54 steps at 8 threads, 0 report files (before: 25 files, 128 reports, 8 sites). Merged `feature/upgrade`
f8b7bf22, GCC 11.5 with GUI, Condor 16120: `tests/upgrade/check.sh` PASS; `ctest -LE slow` 130/131, the one failure
`adaptive_grid_reference` (decision 8).

**Open (proposed default in brackets).**
- `shared_point_tables` check (d), a 4-thread fit's peak memory below 3.0 × the calculation's, fails under TSan
  (1163 MB against 1149 MB) with no report, and failed once in six plain runs. [Factor 3.5 under
  `sanitize_upgrade.sh` only, with the reason in the script; PI to confirm, since it is a test limit.]
- `pyazr_upgrade_thread_independence` failed once: the T3 fit's χ² differed by 2.3 × 10⁻⁵ relative between 1 and 4
  threads, fitted parameters equal; twelve reruns on three branches passed. Cause not found. [Watch; investigate if
  it recurs.]

## 2026-11-06 (sanitizer session) — Phase 0: TSan and ASan at 8 threads, script, CI, findings

Branch `feature/upgrade-0-sanitizers` (from `feature/upgrade` 356c6770), local, not pushed. Commits: c692bb8b
(cherry-pick of 189f3023 from `ai-r/dev`: `PYAZR_OUTPUT_DIR`), 557bb0ba (`AZURE2_SANITIZE`), cfd281ab
(`scripts/sanitize_upgrade.sh`, `scripts/sanitizers/`), 16309cf8 (TSan ignorelist), 8572b21c (workflow), 635b0a1a
(`sanitizer_findings.md`, installation option row), and the roadmap/HANDOFF commit.

**State.**
- `-DAZURE2_SANITIZE=thread|address|address,undefined` instruments every target, Minuit2 and coul included.
- `scripts/sanitize_upgrade.sh thread|address`: 54 steps at 8 threads (T1/T3 calc, extrap, rate, gradient fit; T2
  calc and a 30-call fit; 26 mapping/shift/target `check.sh`; 3 regression projects; 4 integrator reference tests;
  12 pyazr tests in-process). Exit 1 on a report, 2 on a failure without one, 3 on setup.
- `.github/workflows/upgrade-sanitizers.yml`: clang 18 + libomp/libarcher on Ubuntu 24.04, nightly, dispatch, and on
  a push touching the workflow/script/suppressions.
- Runs on Condor (`sanitizer_findings.md` §3): ASan clean with clang 23.1.4 and GCC 15.2. TSan, after the two
  third-party entries (libomp suspend/resume, suppressed; Minuit2 call counter, ignorelisted): 25 report files,
  all at two AZURE2 sites.

**Findings for `fix/` branches (not fixed here; details, stacks and proposed fixes in `sanitizer_findings.md` §4).**
1. §4.1 `fix/fit-iteration-counter` — real race: `EData::iterations_` read and incremented, and the progress line
   written to `std::cout`, from every thread of Minuit2's numerical gradient (`AZURECalc.cpp:565-713`). No fit value
   depends on it. Fix: atomic counter with `fetch_add` returning the iteration, progress line under a mutex.
2. §4.2 `fix/edep-ownership-tsan` — the copy-on-write ownership test `SoleOwner` (relaxed `use_count()` + acquire
   fence, 3fba1192) is ordered under the C++ memory model; TSan does not model standalone fences and reports
   `EDependentTables::Clear` against `ClearEDependentRows` (`shared_tables_race`, `shift_fn*`, pyazr `shift_fn`).
   Fix: an ownership count TSan can see (atomic holders count, acquire load), or wait for Phase 6, which removes the
   tables. The TSan job stays red until one of the two lands.
3. Third party, recorded only: Minuit2's `fNumCall++` from the parallel gradient is a real race (call count and
   `--max-calls` approximate under threads).
4. Not a sanitizer finding: `adaptive_grid_reference` fails on every build of 356c6770, plain GCC 11.5 included
   ("eV level crossing the window edge", swing 1.685e-4 > 1e-4). Needs its own branch.

**Not done / not verified.** The workflow has never run on GitHub (no push from here): the apt layout (runtime names,
`libarcher.so` and `omp.h` locations) is handled by lookups, not tested. The cluster's clang module carries no
compiler-rt and no libomp; the runs used both built from the LLVM 23.1.4 sources in scratch (outside the repository).
GCC TSan was not run: with libgomp it cannot be correct (§2 of the findings).

**Open (proposed default in brackets).**
- GitHub fires the schedule only from the default branch (`qt5`, human-only). [Leave the file on `feature/upgrade`;
  rely on the push trigger and manual dispatch until the upgrade merges, or the PI adds the file to the default
  branch.]
- Should the TSan job be green before §4.2 is fixed? [No suppression of AZURE2 code: take `fix/edep-ownership-tsan`
  next, a small change in `EPoint`; it also removes the last fence-based argument from the copy-on-write code.]
- The ASan job warns, not fails, on `adaptive_grid_reference` (exit 2 is a warning in the workflow). [Keep: build.yml
  gates test failures.]

## 2026-11-07 (T2 session) — test project T2, its golden outputs and the benchmark

Branch `feature/upgrade-0-testprojects-t2` (from `feature/upgrade` ff75fc28), commits 3184b1ee … 8870504f, not
merged.

**State.**
- `tests/upgrade/T2/`: p + ¹⁵N → ¹⁶O through an apparatus, 13 segments, 330 points, 19 free parameters (6 level,
  13 nuisance): thin, thick, implanted (fitted density) and gas (Doppler, E-dependent resolution) targets with stopping
  and straggling; beam profile; Q attenuation; constant, √E and f(E) shifts with groups; a norm group; norm,
  shift, density, abundance and covariance priors; a data covariance; SUM and RATIO through α₁; a natural-target
  sample; THM HOES (two runs, one profiled scale, 20 keV c.m. fold) and OES points. `make_t2.py` holds the truth and
  makes the data with the engine (seed 20261106); README states coverage, truth and TODOs.
- Recovery (`check_recovery.py`, HTCondor 8 threads): noisy χ² 294.17/330, max |pull| 2.15; Asimov χ² 0.17,
  max |pull| 0.38. Golden outputs in `tests/upgrade/golden/T2` (1e-5 calculations, 1e-3 fit;
  `regen_compare.py`). Benchmark `scripts/bench_upgrade.sh`, first row in `tests/upgrade/bench/T2.tsv`:
  calc 2.67 s / 1.77 s, 200-call fit 170 s / 96.7 s, peak RSS 238 / 475 MB at 1 / 8 threads.

**Findings for later phases.**
- With THM segments the CLI forces the analytic-gradient MIGRAD; on T2 it spends > 3 min per gradient call and
  did not converge in 10 min. THM has no analytic residual Jacobian, so `--use-lm` falls back to numerical MIGRAD
  (1082 calls, 585 s). The T2 scripts use `--use-lm`. Phase 2/5: a THM Jacobian.
- Exact Coulomb functions at re-placed sub-points are about 2/3 of a fit evaluation; a Coulomb table makes the T2
  fit 8× faster (72 s) with the same parameters to 8 digits.
- Fit speed-up from 1 to 8 threads is 1.76 and peak memory doubles (the 200-call benchmark).
- In the Brune basis the fit fixes γ, but the `.azr` enters Γ at the level energy: moving a level's energy or a
  free width in the file changes the γ of its fixed channels. A synthetic project must write the fixed Γ that give
  the true γ at the start values (`make_t2.py` `match_fixed`); T1 and T3 need the same care.

**Not done.** S(E) fold of the OES points (no engine support; README TODO); a fitted Doppler temperature (¹⁶O has
no resonance narrow enough); a two-nucleus natural target (¹⁴N(p,α) closed); TSan/ASan runs; T2 is not wired into
CTest or `run_tests.sh` (its fit takes 10 min).

**Next.** Merge the branch. Open: T2 has its own scripts (`T2/check_recovery.py`, `golden/T2/regen_compare.py`),
written in parallel with T1/T3's `recover.py` and `golden/golden.py`; folding T2 into their `PROJECTS` lists is a
small follow-up (proposed default: keep the T2 scripts, and let `tests/upgrade/check.sh` call
`golden/T2/regen_compare.py --only calc_start --only calc_truth`, 3 s, but not the 10-min fit). Decide whether the
T2 recovery runs in nightly CI (10 min at 8 threads, or 72 s with a Coulomb table built in the job).

## 2026-11-07 (third session) — Phase 0: rules into CLAUDE.md; ten-skill set; skill checker

**Branch** `feature/upgrade-0-rules-skills` (from `feature/upgrade` ff75fc28), local, not pushed.

**Done.**
- `CLAUDE.md` "Skills" section replaced by `RULES.md` §3 (use / behaviour / change; sizes; no history except
  expiring version notes; one owner; descriptions ≤ 600; no personal paths; the checker). The "date an entry" rule
  is gone; the frames-and-units pointer now names `azure2-check`.
- `.claude/skills/` is the ten-skill set plus `rmfit-campaign` (2,280 → 358 lines + three references). The ten old
  skills are deleted. `docs/upgrade/skills_migration.md` maps every section of every old skill (by line range in the
  ff75fc28 files) to its new location or a reason for dropping it (498 rows; coverage checked by script).
- `tools/check_skills.py` (R4, R5, R6, R7, R8, docs/reference/script links), its unit test
  `tests/consistency/check_skills_test.py` (12 cases), ctest `check_skills_unit` and `check_skills`, a local
  `.pre-commit-config.yaml` hook and `.github/workflows/skills.yml` (runs on push/PR touching skills or docs). The
  checker exits 0 on the branch, also with `--strict`.
- Developer content from the skills moved to `tests/CLAUDE.md` (CI portability, regression vs reference, m-scheme
  traps) and `gui/CLAUDE.md` (headless dialog checks). `reference/command_line.rst` gained the gradient-diagnostic
  and `AZURE_DUMP_NONFINITE` variables.

**Not checked.** `sphinx-build -W` (Sphinx is not installed on the analysis machine; the two RST edits were read by
eye). The CI workflow has not run (nothing is pushed).

**Open.**
- Version notes in the skills expire on 2027-02-01; the checker then fails until each is removed or renewed.
- The user-level `~/.claude/skills/azure2` must narrow its description so it does not load next to the repo skills
  (`RULES.md` §3, outside the repository).
- `pyazr/azrfile.py` (a docstring near line 3255) still says "see the azure2-eval skill"; left alone because this
  branch changes no engine or pyazr code. The owner is now `azure2-calculate`.
- Corrections made while moving content (reviewer please confirm): masses in the old project builder were called
  nuclear (they are atomic); `nuclear_mass()` and `zero_all_widths()` cited by old skills do not exist in pyazr; the
  CLI band prompt is asked in modes 2, 3 and 4; the level-merge rule follows `levels_channels.rst`; in
  `azure2-resonance-search` the gain is "control minus add" (positive = improvement).
## 2026-11-06 (evening, EST) — Phase 0: T1, T3, golden outputs, identity checks (branch `feature/upgrade-0-testprojects-t1t3`)

**State.** Branch `feature/upgrade-0-testprojects-t1t3`, cut from `feature/upgrade` c1cc92c7 (feature/upgrade ff75fc28
merged in for STYLE.md). Ready for review and merge; nothing pushed.
- 63c08471 `tests/upgrade/T1` (p + ¹²C → ¹³N, 442 points, 13 free parameters) and `tests/upgrade/T3` (⁸Be, 494
  points, 8 free parameters): generators with the true parameters, synthetic data (seeds 13 and 8), `.azr`, start
  models, READMEs with truth, data tables and coverage; `tests/upgrade/lib/upgrade_project.py`;
`tests/upgrade/recover.py`.
- 3b48a101 golden outputs of engine c1cc92c7 (`tests/upgrade/golden/`, calc/extrap/rate, `--no-coulomb-table`,
  2 threads) with `golden.py` (compare at rtol 1e-6, `--write` to regenerate); `tests/upgrade/check.sh` (golden +
  recovery, run by `tests/run_tests.sh`, ~50 s).
- fbad9b6a, bf5f1713 seven identity tests `tests/pyazr/upgrade_*_test.py` (ctest `pyazr_upgrade_*`; the thread test is
  `LABELS slow`); `tests/upgrade/README.md` lists every §10 identity with result, and the pending ones.

**Verified (feature/upgrade build c1cc92c7).** Recovery: T1 χ² 462.3 at the fit vs 482.3 at the truth (N = 442), largest
pull 3.64 (5/2⁺ energy), noise-free fit ≤ 0.007 σ; T3 χ² 498.7 vs 509.7 (N = 494), largest pull 1.82, noise-free
≤ 0.036 σ; eight extra noise seeds for T1: pull means within ±0.7, standard deviations 0.7–1.5. Golden: two runs agree
bit for bit; an edited golden value fails. Identities: Coulomb F, G vs mpmath 4.5e-10, Wronskian 4e-16; unitarity to
print precision; optical theorem 2.6e-11; Breit–Wigner area 1.1e-7; Park = Brune exactly, formal = Brune 2.2e-11
(see open 1); gradient vs Richardson differences 1.2e-8 / 2.7e-6 of max|∇χ²|; threads 1 vs 4 equal to 3.6e-11 (open 2).
A build of this branch (CLI + pyazr, no GUI) ran `ctest -R "pyazr_upgrade|consistency"`: 8/8 passed (none
skipped), and `tests/upgrade/check.sh` passed. The full `ctest -LE slow` and `tests/run_tests.sh` were not run: the
branch adds tests and changes no engine, pyazr or GUI code.

**Open (proposed default in brackets).**
1. T1 capture **with external capture** in the formal basis (`--no-brune`) differs from Brune/Park by up to 54 %
   (E_c.m. = 1.568 MeV); without external capture it agrees exactly. The EC amplitude depends on B_c and the
transformation
   of observed photon widths to the formal basis does not carry it. [Treat Brune/Park as the reference; resolve in
Phase 1
   with `level_parameters`; meanwhile a warning for `--no-brune` with an EC mask would protect users — a `fix/` branch
   on dev, after James's group agrees.]
2. Thread dependence at 3.6e-11 relative in T1 capture to the excited (3/2⁻) final state, in the calculation and after a
   fit, also between two 4-thread runs; fitted parameters identical. [Accept until Phase 6; its gate tightens
   `upgrade_thread_independence` to bit for bit. Worth a TSan run on T1 in the Phase 0 CI box.]
3. Neutral channel F_l at small ρ and l ≥ 3 is computed with cancellation (7.9e-6 relative at l = 4, ρ = 0.18); P, S,
   δ_hs unaffected. [Fixed by the pycoulomb functions of Phase 1; the test then tightens to 1e-8.]
4. Photon widths are not in the level denominator (Γ = Γ_particle); `upgrade_breit_wigner` pins it. [Keep; state it in
   `docs/source/theory/conventions.rst` in Phase 1.]
5. `tests/beta_delayed_alpha` labels its energy column "E_cm(12C+a)", but the engine reads a β-delayed data energy as
the
   lab energy of the light particle from the decay at rest (E_c.m. = E·16/12 there). The data are consistent with the
   engine, only the label is wrong. [Fix the label on dev.]
6. Golden files are generated outputs outside `tests/*/expected/`, which the root CLAUDE.md "Commits" rule allows only
   there. [Name `tests/upgrade/golden/` as the second allowed place in the Upgrade section.]
7. T1's convention round trip is a TODO for Phase 1 (no `Convention` record yet).

**Next.** Merge this branch; add T2's golden files and recovery to `golden.py` / `recover.py` (`PROJECTS`, the project
list) when T2 lands; the Phase 0 CI box should run `tests/upgrade/check.sh` and the `pyazr_upgrade_*` tests.

## 2026-11-07 (second session) — branches imported; THM merged in; Phase 0 started

**How to continue on any machine.** Pull `feature/upgrade` and any `feature/upgrade-<phase>-<topic>` branch listed
below, read this file, then `ROADMAP.md`. Every phase branch is cut from `feature/upgrade` and merged back with
`--no-ff` when its gate passes. Nothing in the repository depends on the machine it was written on: scripts take the
engine binary and build directories from arguments or environment variables (`AZURE2_BINARY`, `PYTHONPATH` for the
pyazr module) and document their defaults. Commits are local until Jakub pushes them; an agent never pushes.

**State.**
- `feature/upgrade` now contains: the plan (f99661cf); `origin/thm` c6ddb741 merged in (d20718f8), so the rebuild
  starts from the THM line (which already carries `dev` up to 69826014; the last five `dev` commits, aab1de0d and
  its four parents, come through the plan commit's base); `fix/edep-cow-acquire` 3fba1192 merged (37978592); an
  introduction in `README.md` and the `CLAUDE.md` Upgrade section (c1cc92c7).
- `thm-review` (4882ad04): the uncompiled `vertexModel=ff` GUI code from the analysis machine is applied in the
  working tree and being built and tested before it is committed; then `origin/thm` is merged into it.
- Phase 0 work in progress on `feature/upgrade-0-testprojects-t1t3` (T1, T3, identity checks) and
  `feature/upgrade-0-testprojects-t2` (T2, `scripts/bench_upgrade.sh`). Next: `feature/upgrade-0-rules-skills`.

**Open.** The decisions listed in the entry below are unchanged. The race-fix status of `thm` (a6bc6b9) is resolved on
this branch by the merge; `origin/thm` itself still lacks it.

## 2026-11-07 — plan written; branch created; nothing implemented yet

**State.**
- Branch `feature/upgrade`, cut from `origin/dev` aab1de0d. Contains only this plan (`docs/upgrade/`) and the "Upgrade"
  section of `CLAUDE.md`. No code changed.
- The plan: `ARCHITECTURE.md` (target design: three lines, modules, one-way rule, stopping, Coulomb via pycoulomb,
  conventions, physics explorer, tests, front ends), `THM.md` (simple valid THM with DWBA vertices), `ROADMAP.md`
  (phases 0–9 with checklists and gates), `RULES.md` (rules and the ten-skill set), `figures/structure.html` (open in
  a browser; generated by `figures/make_figures.py`).
- The plan came from six read-only surveys of `dev` (data/point model, experimental effects, fit loop and threading,
  fix history since 2026-08-01, classes vs physics with pyazr, skills), done on 2026-11-06.

**Work that exists only on the analysis machine (almalinux), shipped as patches in the handoff package.**
- `thm-review` (base `origin/thm` 30d3bd9, local, never pushed): 073fe48 `backgroundPositive=on` (engine, pyazr, docs,
  tests), 08ccea3 its GUI box, 4882ad0 `vertexModel=ff` full form factor (engine, pyazr, docs, tests; validation:
  Green identity 2.6e-10, R values vs formula sheet ≤ 1e-6, `ffVolume=off` byte-identical to `perlevel`, light suite
  69/69, THM ctests 30/30). Plus uncommitted GUI code for `ff` in 6 files (`gui/include/ThmExperimentsPage.h`,
  `ThmSettings.h`, `gui/src/ThmExperimentsPage.cpp`, `ThmSettings.cpp`, `ThmWorkspace.cpp`,
  `tests/gui/thm_workspace_test.cpp`): written, **not compiled** (the machine had < 2.7 GB free). Build `AZURE2` and
  `thm_workspace_test` in the GUI build, run the test, then commit in the style of 08ccea3. Note: `origin/thm` has
  since merged dev (c6ddb741); rebase or merge accordingly.
- `fix/edep-cow-acquire` (base `origin/dev` e8727574): 3fba1192, acquire fence on the copy-on-write ownership test of
  the shared point tables (`SoleOwner()` in `src/EPoint.cpp`) and the "a published table is never written again"
  invariant in `include/NodeTableCache.h`. Checked: race test passes at 4 threads; ctest 57/59 and physics 90/96 with
  exactly the same failures on unpatched dev (memory limit of that machine, plus `component_mapping` (c),
  `total_optical` (a), `pyazr_failed_point`, `pyazr_composite_terms`, `shared_point_tables` (d) — worth checking on a
  bigger machine; they may be real dev regressions). Not checked: clang TSan, macOS.
- `origin/thm` lacks the use-after-free fix a6bc6b9 that is on `dev`; merge dev into thm.

**Not in any repository.** `~/THM/` on the analysis machine (THM notes, evaluations, Marco's review in
`~/THM/review/marco_2026-11/`: `MARCO_REVIEW.md`, `IMPLEMENTATION_PLAN.md`, the formula sheet, meeting notes; the
15N D3 refit in `~/THM/evaluations/15N_p_a0/review2026/`). Copy it separately if needed (private 24Mg data must stay
out of any repository).

**Decisions open (Jakub, then James deBoer's group).**
1. Agree the structure and the phase order; who owns which phase.
2. The Phase 0 freeze on perf work in the point model and on new Qt features.
3. Retiring Qt for the web UI after one release with both. **Decided 2026-11-07 (Jakub): no Qt GUI in the upgrade.**
   Builds and checks run with `-DBUILD_GUI=OFF`; `gui/` and `tests/gui/` are removed on their own branch; the web UI
   (Phase 9) is the only graphical front end. Qt stays on `dev`/`qt5` for released versions.
4. The quadrature tolerance rule for Phase 5.
5. Whether pyazr keeps a pure-Python `.azr` reader after Phase 3 (plan: no).
6. pycoulomb licence and the origin of `coulsm` (Jakub, owner). 2026-11-07: Jakub says to take pycoulomb from
   github.com/skowrons94/pycoulomb (his); cloned at 9fee957. Still open: its LICENSE and the origin and terms of
   `coulsm` (SAMMY) and `coulcc` (CPC program library), before either enters this repository. **Decided 2026-11-07
   (Jakub): pycoulomb's C++ sources may be taken directly into AZURE2; permission for them obtained (incl. from Carl).
   Their purpose here is to switch between Coulomb algorithms.** The Fortran `coulcc` is not taken (C++ only).
   Checked by diff: pycoulomb `coulsm/` is AZURE2's present `external/coul` (the SAMMY-side patch, `cwfcomp_accurate`,
   in AZURE2 since 17fd9ded, 2025-09-05) up to a namespace and comments; pycoulomb `coul/` is AZURE2's coul before
   that patch. The default algorithm is therefore `coulsm`, and it must reproduce today's numbers bit for bit.
7. Redistribution terms of any stopping table before shipping one.
8. (Added 2026-11-07.) The criterion of `adaptive_grid_reference` for a level crossing the window edge: thm's grid
   overcounts a level's area by 1.69 × 10⁻⁴ by design, so the swing equals that offset and exceeds the 10⁻⁴
   tolerance. Options (a) test the variation beyond the net change, (b) finer tail steps (moves every fold by about
   10⁻⁴, regenerates the goldens), (c) take dev's grid. Analysis in the `fix/adaptive-grid-edge-swing` entry.
   [Proposed: (a).] **Decided 2026-11-07 (Jakub): (a).** **Implemented by 47c119e4** (`fix/adaptive-grid-criterion`):
variation
   beyond the net change, tolerance √3/(2π ppw²) + 10⁻⁴.
9. (Added 2026-11-07.) The nightly sanitizer schedule fires only from the default branch (`qt5`). [Keep the file on
   `feature/upgrade`; push trigger and manual dispatch until the PI decides.]
10. (Added 2026-11-07.) Formal basis: B_c from the lowest level instead of the first listed, and each formal
   level's largest |γ| positive. Both exist as `Convention` options, off (`level_parameters.md`). [Turn both on when
   the engine switches to `level_parameters`, with a formal-basis golden re-baseline.] **Decided 2026-11-07 (Jakub):
   as proposed.** Implemented 2026-11-07 (`feature/upgrade-1-cnuc-split-1b`) with B_c from the lowest level; the sign
   rule as built is the orientation of the Brune level, because "largest |γ| positive" changed external-capture
   observables (step-1b entry, Found 1). **Confirmed 2026-11-07 (Jakub): keep the Brune orientation.**
11. (Added and decided 2026-11-07, Jakub.) Exact energies by default: observables are calculated at the data
   energies, shifted energies and fold nodes, as AZURE2 does today; interpolation on a fixed mesh only as an opt-in,
   checked mode (`ARCHITECTURE.md` §1, §3, §4).
12. (Added and decided 2026-11-07, Jakub.) Data must be easy to adapt (sums, ratios, yields, E1/E2 parts, …): a
   measurement's prediction is a model expression of observable selectors and combining nodes (`ARCHITECTURE.md`
   §3 experiment/, §9). Designed in Phase 4 with `io/Project` reading today's segment types into it.

**Next actions.**
1. Push `feature/upgrade` (and, if wanted, `thm-review` and `fix/edep-cow-acquire`) from a machine with GitHub access.
2. Phase 0, first tasks: T1–T3 generators and projects, golden outputs, identity tests, benchmark script.
3. pycoulomb prerequisites (Phase 1, first box) can start in parallel in the pycoulomb repository.
4. Build and commit the `ff` GUI code on `thm-review`; the THM re-analyses (D1–D6 of the Marco review plan) continue
   on `thm-review` until Phase 5.
