# The AZURE2 upgrade, for the owner of AZURE2

Every command in section 2 was run on a build of this branch; every other number names the test or document that
holds it. The page as it grew over two weeks is `history/FOR_THE_OWNER_2026-11.md`; the record of merges `HANDOFF.md`.

## 1. What this branch is

- **One engine in libraries with stated dependencies.** `physics/` (formulas, links nothing), `io/` (the project
  file), `stopping/`, `experiment/` (the project format's table, units, apparatus kernels), `convert/` (old formats),
  the engine (`src/`, `include/`), the binding (`api/`). A wrong dependency does not compile (section 5).
- **Project format 3 with a converter.** A data set is a named record with units; the engine reads format 3 only and
  `AZURE2 --convert` turns any older file into it, reporting every change.
- **Commands instead of the menu**: `AZURE2 calculate|fit|extrapolate|minos|rate|mcmc [options] FILE`, no input.
- **pyazr 3.0** on the engine's own session: `Project`, `Session`, `Objective`, `Parameters`, `Evaluation`; no
  physics formula and no project grammar is written in Python.
- **A web page instead of the Qt program** (`python3 -m pyazr ui FILE`), generated from the same tables.
- **ENDF-6 in and out** (`pyazr.endf`: read, check with the IAEA codes when present, start a project, write).
- What it is not: a new reaction model. It stays an R-matrix code; where its numbers move is section 3. New data
  kinds (finite-layer yields, S factors, ratios to Rutherford, literature values as χ² terms) and its own SRIM-2013
  stopping are in the user guide.

**What it costs you**, measured against `dev` at `3e502625`: the C++ outside `external/`, `gui/` and the tests is
611 files and 129,400 lines (dev: 143 and 55,700), the Python 119 files (dev: 37); about 1,700 commits that nobody
can review line by line, so the tests are what there is to trust (section 5); `CNuc`, `EData`, `EPoint` and pyazr 2
do not exist here, so nothing written against them merges (section 7); every older project needs the converter.

## 2. Try it

Run here with `build/src` on `PATH`, `PYTHONPATH` the checkout, `AZURE2_COULOMB_TABLE` naming a Coulomb table
(without one the exact functions are used), `OMP_NUM_THREADS=2`; copied from the terminal, `[...]` marks a cut.

```
$ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -Dpybind11_DIR="$(python3 -m pybind11 --cmakedir)"
$ cmake --build build --parallel 2          # 5 min 43 s here, 0 warnings; the engine module goes into pyazr/
```

An old project (`tests/convert/fixtures/tests/13N/13N.azr`, as `13N_old.azr` beside `tests/13N/data`):

```
$ AZURE2 calculate 13N_old.azr
ERROR: 13N_old.azr is a format-1 project (<segmentsData> at line 66). This AZURE2 reads format 3 only.
       Convert it:  AZURE2 --convert 13N_old.azr 13N_old_f3.azr     (or pyazr.convert_azr("13N_old.azr", [...]
       The converter keeps 13N_old.azr unchanged and reports every change it makes.
$ AZURE2 --convert 13N_old.azr 13N_new.azr
  note: artemov: a thick target (density 0): the old engine compared the numbers as reactions per incident [...]
  note: seagrave: a thick target (density 0): the old engine compared the numbers as reactions per incident [...]
  note: 16 row(s) for parameters format 3 does not have (a fixed zero shift, an inactive or dropped line) removed
Wrote 13N_new.azr (format 3): 3 note(s), 0 question(s).
$ cmp 13N_new.azr 13N.azr && echo IDENTICAL          # 13N.azr: the project of tests/13N
IDENTICAL
```

The synthetic project T1 (`tests/upgrade/T1`: p + ¹²C, capture to two final states, elastic scattering, A_y):

```
$ mkdir output checks                                # the directories its <config> names
$ AZURE2 --validate T1.azr
T1.azr: 0 diagnostic(s), valid.
$ AZURE2 calculate T1.azr                            # 2.0 s, 105 MB
Total Chi-Squared: 482.275239617
$ AZURE2 fit T1_start.azr                            # 2.6 s, 107 MB; T1_start: the same model away from the truth
MIGRAD converged: valid minimum, EDM 4.243e-06 (goal 2.000e-04), 733 function calls, call limit none.
Total Chi-Squared: 462.273515393
$ AZURE2 rate --rate-entrance=1 --rate-exit=2 --temperature-file=temperatures.dat T1.azr     # 9.5 s
$ cat output/reactionrates.out
                  T9                Rate
                0.01          1.9123e-22
                0.03         9.17279e-14
[...]
$ AZURE2 mcmc --walkers=8 --steps=100 --sampler-threads=2 --basis=rwa --overwrite-samples T1.azr    # 3.2 s
Acceptance fraction: 0.361
Total samples generated: 808
$ python3 -m pyazr ui T1.azr --port 28750 --token owner-page-test
http://127.0.0.1:28750/?token=owner-page-test
```

Its `GET /api/evaluation` answered the objective 482.27523961726723. Five lines of pyazr, in the same directory:

```
$ cat session5.py
import pyazr
s = pyazr.Session("T1.azr")
x = s.parameters.x0
print(s.evaluate(x).chi2, len(x))
print(s.objective().value_and_grad(x)[1][:3])
$ python3 session5.py                                # 0.32 s
Coulomb table: [...]
482.27523961726723 13
[ -344.7239501    -93.58136105 52116.49415612]
Coulomb table: 678 lookups, 0 outside the table (exact routine used)
```

## 3. Where the numbers differ from `dev`, and why

`dev` here is the commit where this branch left it (`aab1de0d`) unless a row says otherwise. "Decides": who has to
agree before `dev` follows; "PI" means the PI took the decision and you should know it.

| # | What differs | Size | Since | Decides |
|---|---|---|---|---|
| 1 | Every reaction cross section out of an identical entrance pair has 1 + δ₁₂ | factor 2 | 2026-09 | owner |
| 2 | Column 32 of `<levels>` (a pair's THM binding energy) | E_x moved by B | found 2026-11-15 | owner |
| 3 | B_c of the standard basis at the lowest level of a J^π group | 0 in σ | 2026-11-07 | owner |
| 4 | Two levels of one J^π closer than 1e-6 MeV (dev: 1e-9) | ≤ 1e-5 of χ² | 2026-11-15 | owner |
| 5 | Coulomb functions just above a heavy-ion threshold, without GSL | see line | 2026-11-15 | owner |
| 6 | The centre of a prior is the number written, for every kind | prior kept | 2026-11-13 | PI |
| 7 | The sampler: −½ × the fit's objective, zero outside every limit | see line | 2026-11-13 | PI |
| 8 | Folded data on one shared energy mesh per compound nucleus | ≤ 5e-3 | 2026-11-12 | PI, owner |
| 9 | `useAdaptiveGrid=0` removed (uniform grid of a point's own) | up to 57 % | 2026-11-12 | PI |
| 10 | Fold windows removed: every folded point gets the whole kernel | 1.5 σ, 1 point | 2026-11-12 | PI |
| 11 | A shifted folded point: stopping and thickness at E + shift | 0.9 % | 2026-11-13 | PI |
| 12 | A point given at a lab angle is folded at that lab angle | ≤ 2.2 % | 2026-11-12 | PI |
| 13 | A capture term that is only a later term of a sum or ratio | 0.54–2.9 | 2026-11-11 | owner |
| 14 | One final state on two photon pairs: refused | ×2.2 | 2026-11 | owner |
| 15 | Analytic gradient of an E1-only or E2-only sum with a component | gradient | 2026-11 | owner |
| 16 | A MIGRAD fit stopped by its call limit ends at its lowest objective | see line | 2026-11-14 | owner |
| 17 | Radiative widths are not in a level's total width (as on `dev`) | none | — | — |

1. σ(12→34) = (π/k²) g_I (1 + δ₁₂) Σ …: a ¹²C + ¹²C reaction is twice `dev`'s, which applies the 2 to elastic
   scattering only. `tests/identical_entrance_reaction`; `docs/source/theory/conventions.rst`, Eq. conv-sigma.
2. `dev` (the branch point and `3e502625` alike) reads that column as the light particle's excitation e1 and moves
   the level's E_x by it. Found on a 12C + 12C model of the PI, not in the repository (`HANDOFF.md`, 2026-11-15).
3. Observables do not move while the transformation is exact; formal parameters of a group listed out of energy
   order (`param.sav`, `--no-transform` input) belong to another B_c; the run says so. `conventions.rst`.
4. Below the gap the limit takes the mean of the two levels' S′ (and S″), the chord's slope to second order; just
   above 1e-9 MeV the divided difference loses 1e-5 of χ² to round-off, at 1e-6 it loses 4e-10 |S|. One header,
   `physics/include/azure2/physics/LevelPairShift.h`, for the six places that divided by the gap.
   `pyazr_coincident_levels`; `physics_level_parameters_unit`, section 8.
5. Your repair with this branch's own e^z K_n(z) where the continued fractions give up (η ≳ 230): ¹²C + ¹²C, l = 0,
   1 eV–2.5 keV gives S = −9.594 where the old series gave +4.6 to +11.2. At the switch (3.6 keV) S steps by
   2.3e-3 (2.4e-4 relative). Not compared number by number with `dev`'s GSL route. `physics_coulomb_threshold`.
6. `isotope 15N 0.0037 fit prior 0.0036 0.0002` pays its prior from the first evaluation; `dev` centres an
   abundance prior on the value of its line, so fitted values copied into the file move the prior with them.
   `parameter_kinds.md` section 4; `tests/nuisance_prior_centre`.
7. A chain that never reaches a limit is unchanged bit for bit; `screening_6li` (priors at −300 eV): medians +131,
   −217, −205 eV before, 0.86, 0.60, 0.24 eV now. `objective.md` section 6; `sampler_objective`.
8. A folded number moves by the quadrature error it had: of 577 files 542 identical, 33 as close to an eightfold
   refined reference or closer, 2 of `hybrid_potential` (too coarse settings). Accepted when the largest distance is
   ≤ 1.2 times the old engine's and the r.m.s. ≤ 1.03 times (`energy_grid.md` section 7.1).
9. `examples/c12c12_tumino2018`: χ² 111.661 → 117.516 (117.516 again with 8 times the points); a 9.7 eV level in a
   thick layer: mesh 2.18022e-6, uniform grid 2.06e-8 to 5.09e-7 (64 times the points). `tests/shift_sqrt_subpoints`.
10. `beam_spread_ay`: χ² 3838.66 → 3760.95, one A_y point by 1.5 of its error; equal, bit for bit, to the old engine
    with the window widened to all energies. A window restricted the physics; it was not an approximation.
11. A yield on a resonance was off by (d ln ε / d ln E) × shift / E: 0.9 % for 10 keV at 270 keV. T2's truth
    file: χ² 330.43 → 330.84. `tests/shift_true_energy`.
12. In a resolution window the c.m. angle and dΩ_lab/dΩ_c.m. follow the energy: T2 moves by 2.6e-5; a 300 keV
    resolution 0.4 MeV above the ⁷Li(p,n) threshold at 30° lab was 2.2 % off. `tests/lab_frame_fold`.
13. `dev` decides external capture from the first term of each reaction only, so such a term runs without its
    external part and without a message (0.54 to 2.9 times the right term). `tests/component_mapping`, case (d).
14. `dev` reads it without a message and one pair's capture loses its external part (×2.2 at the first point).
    Refused here, naming both pairs. `tests/ec_final_state_ambiguous`, case (d).
15. Wrong on `dev` by 0.29 and 1.0 of the largest element against central differences; such gradient fits change.
16. Minuit2's last accepted state can be the start: T2 with `--max-calls=40` ends at 10910.1338428 before,
    7087.92418676 now, with a notice; `fitStatus.out` has one writer and new keys (`fitter`, `fitter_status`).
17. The resonance strength stays ω Γ_γ Γ_p / Σ Γ_particle, as on `dev`. New is only the record `gamma_decay` (Γ_γ,
    branchings, mixing ratios; `user_guide/gamma_decay.rst`), which gives the same partial widths as `<levels>`.

Rows 9, 15 and the refusals that stop a project which ran on `dev` (two resolution kernels on one data set, a lab
point at the kinematic limit, ...) are in the history page, "What users notice".

## 4. Speed and memory

Elapsed time and peak resident memory of one command, the same project for both programs, the old one through its
menu; a machine shared with other work, so read the ratios. "Branch point" is `aab1de0d`: the 65 commits `dev` has
added since (to `3e502625`) are not in the comparison. Sources: (a) `HANDOFF.md`, 2026-11-15, the director's entry
of 04:35 and the entries it merged, and the history page's "Speed"; (b) the director's entry of 2026-11-15 09:25;
(c) the director's measurement of 2026-11-11 (history page, "Phase 5 so far").

| Run | Threads | Branch point | This branch | Source |
|---|---|---|---|---|
| 13N, one calculation | 1 | 9.3 s, 216 MB | 7.3 s, 122 MB | a, b |
| 13N, one calculation | 4 | 7.3 s, 220 MB | 5.5 s, 123 MB | a |
| 13N, fit with the analytic gradient | 1 | 92 s, 604 MB | 17.6 s, 136 MB | a |
| 13N, fit with the analytic gradient | 4 | 30.6 s, 610 MB | 8.3 s, 140 MB | a |
| 13N, Levenberg-Marquardt fit | 1 | 28.6 s, 295 MB | 10.1 s, 136 MB | a |
| 13N, plain MIGRAD stopped at 300 calls | 1 | 135 to 141 s, 597 MB | 12.6 s | a, b |
| 13N, plain MIGRAD stopped at 300 calls | 4 | 44.9 s, 602 MB | 7.4 s | a, b |
| 13N, sampler, 16 walkers, 100 steps | 4 | 322.6 s, 3.96 GB | 21.8 s, 136 MB | a |
| `natural_neon`, one calculation | 1 | 41.2 s, 4.52 GB | 1.00 s, 289 MB | c, b |
| a 12C + 12C model of the PI, extrapolation | 1 | 86.9 s, 4.47 GB | 41.5 s, 1.64 GB | a, b |
| the same | 4 | not run | 14.3 s, 1.64 GB | b |
| T2, one calculation | 1 | not run | 0.19 s, 116 MB | b |
| T2, sampler, 16 walkers, 100 steps | 4 | not run | 32.0 s, 260 MB | b |

13N's χ² is the same (336480 printed by the old program, 336479.961 now), so are the converged fits' (28100.2,
28100.172); the stopped fits and the samplers did the same amount of work, not the same steps. The 12C + 12C model
(eight excitation functions folded on one mesh, four of them THM; the old program folds each point on a grid of its
own, ten energies) took 10.99 GB here until a compound energy's
channel rows were held once (`data_model.md` section 9). The gain: one level-matrix solve per compound energy
instead of one per point or sub-point, one calculation site per mesh node, no copy of the compound per thread.

## 5. The structure on one page

The libraries and what each owns: `docs/source/developer/architecture.rst` and its figure `figures/libraries.svg`;
the classes: `classes.rst` and the four figures beside it; who owns what: `state_and_ownership.rst`. The life of
one calculation, in the program and in pyazr alike:

1. **Project → records.** The file is read once (`Config::ReadConfigFile`); `ReadFormat3Project` validates every
   record against the table (`experiment/src/Schema.cpp`), `Schema::Build` makes the specs. Constant afterwards.
2. **Compound system and data set-up.** `SetUpCompoundSystem` (`include/CompoundSystemSetUp.h`): the pairs, J^π
   groups, channels and levels from `<levels>`, then the caller's data set-up (segments, kernels, the energy mesh of
   each compound nucleus, `PlanMesh`), then the last part (B_c, the input transformation, the pathway table, the
   external-capture set-up). It returns a `const CompoundSystem` shared by every thread, and the owner's
   `LevelState` (the levels at a parameter vector) and `LevelInputs` (the input values).
3. **Working copy on the shared mesh.** The data as set up are one immutable `DataSet`; what an evaluation changes
   (norms, shifts, apparatus values) is a `WorkingCopy`. Each mesh node is one calculation site (`Site`).
4. **Model pass.** One level-matrix solve per compound energy (`src/EnergyPass.cpp`) into the caller's `Results`,
   then each folded point over its sites (`src/PointFold.cpp`); each data set's value is an expression tree (sum,
   ratio, scale, yield) with one forward and one adjoint implementation (`src/ModelPass.cpp`).
5. **Objective.** One ordered list of terms (`src/Objective.cpp`): data χ², nuisance priors, Park's wall, model
   data, prior covariance. The fitters (`src/Fitters.cpp`, one row each), the sampler and pyazr's `Objective` all
   evaluate this list.

**No hidden state.** Only constant tables are process-wide; what an evaluation writes belongs to its caller.
Check 38 of `tests/consistency/repo_consistency_test.py` refuses any other `static`, `thread_local`, global or
`mutable`, except the 42 rows of `tests/consistency/allowed_state.tsv`: 14 permanent (GSL's error handler, the
Coulomb tables, the cache of input files, output-stream locks, the level-state version counter) and 28 still to
move (11 members of the fitter's adapter `AZURECalc`, 16 of the sampler's `AZURECalcMCMC`, 1 of `Site`). No source
file may grow: the 12 hand-written files above 1,000 lines have a ceiling in `tests/consistency/file_sizes.tsv`.

## 6. How to add to it

Each kind of extension is a row of one table and a test (`docs/source/developer/extending.rst`; with more detail in
`CLAUDE.md`, "Recipes"). A forgotten place fails the named test or check of `tests/consistency/`:

| Kind | The named places | The test that fails |
|---|---|---|
| A data quantity | `experiment/src/Quantities.cpp` row; its tree in `Schema::Build` | `experiment_quantities` |
| An observable at a mesh energy | `SiteObservable` row; `CalculateBareCrossSection`; `PointTBar` | gradient check |
| An apparatus kernel | one `experiment/src/Kernel<Name>.cpp`; a line in `Kernels()` | `experiment_kernels` |
| A parameter kind | `ParamKindTable()` row; a `FillMnParams`; `BuildParamIndexMap` | `pyazr_parameter_records` |
| A level parameter | a `LevelParameterLayout` slot; `LevelState`, `LevelInputs` | `upgrade_level_parameter_layout` |
| A set-up fact of the compound | a `CompoundSystem` member written by `Fill` or `Initialize` | its own test |
| A term of the objective | an `ObjectiveTerm` in `src/Objective.cpp`; gradient in `src/AZUREGrad.cpp` | check 11 |
| A record kind, key or section | `experiment/src/Schema.cpp`: `FieldDef`, `Kinds()`, `Sections` | `experiment_schema` |
| An option of the calculation | `io::ConfigKeys()` row; `OptionStorages()`; the `Session` keyword | check 1 |
| A command or command option | `MakeQuestions()` row in `src/RunRequest.cpp`; `ApplyAnswers` | `tests/cli_help` |
| A fitter | `src/FitX.cpp`; its row in `Fitters()` | check 46 |
| A physics function for pyazr | the library function; one row in `api/src/pybind_physics.cpp` | `pyazr_physics_table` |
| A route or tab of the web page | `pyazr/ui/routes.py` row; the tab's script | `pyazr_ui_routes_consistency` |

**Measured** (history page, "How to add a feature now"). An assistant given only `CLAUDE.md` added three kinds of
data (2026-11-12): a modified S factor (19 min; 5 engine files, +27 −16 lines), the Legendre forms `rsl`, `2l2`,
`4pi` (18 min; 4 engine files), a cross section relative to a tabulated standard (17 min; 7 engine files, +68 −7);
the recipes were silent 11 times and wrong twice, each now a sentence of a recipe or a check. A triangular
resolution with a fitted width: one kernel file of 68 lines, no engine file, about 6 minutes.

**Your three coming features, traced through this structure by reading the code** (no code was written for them).

**(a) rmfit beside pyazr.** rmfit's code is not on `origin/dev`: only `.claude/skills/rmfit-campaign/MOVED.md`
(the skill moved to rmfit's repository); this branch still carries a ported copy of the skill. Dev's last copy of
the skill names the calls rmfit makes. Their pyazr 3.0 names (`docs/source/user_guide/pyazr.rst`, "Earlier names"):

| rmfit on pyazr 2 | pyazr 3.0 |
|---|---|
| `pyazr.azure2(file, **flags)` | `pyazr.Session(file, **options)` (the same keyword names, also `<config>` keys) |
| `residual_jacobian(x)`, chi2 and gradient | `s.objective().residuals_and_jacobian(x)`, `.value_and_grad(x)` |
| per-segment χ², penalties | `s.evaluate(x).data[i].chi2`, `.terms` by the engine's names |
| `AzrModel.add_data_segment`, a baked `.azr` | `pyazr.Project.open(f).add_dataset(...)`, `save()`, `Session.save_fit` |
| `refreeze_covariance()` | `Session.refreeze_covariance()` |
| `pyazr.nds.fetch_levels` | `pyazr.fetch_levels` (module `pyazr.sources.nds`) |
| its own fold of a target-integrated segment | none needed: the engine's kernels on the shared mesh |

`Objective` also gives `bounds` and `log_prob`; a fit can be the engine's own (`s.fit_with_engine("lm")`). rmfit's
rule "never two sessions in one process" no longer holds (`pyazr_concurrent_sessions`). It would be a package beside
`pyazr/` importing pyazr 3.0 only. Not known without trying: whether its shards, sign rounds and level scans need
anything the 3.0 surface lacks.

**(b) Capture into an unbound final state.** Where the bound-state assumption lives today:
`physics/include/azure2/physics/ExternalCapture.h` (N_f with the Whittaker integral I_fc of each closed channel),
`src/NFIntegral.cpp` (that integral, binding energy |E_f − S|), `src/ECIntegral.cpp` (the radial integral with the
final-state Whittaker function at |E_f − S|, and a memo of it for the rate's grid), `src/ExternalCaptureFinalStates.cpp`
(`ParseExternalCapture`: the ANC-to-amplitude factor from the Whittaker function at the channel radius),
`experiment/src/FinalStates.cpp` (the final state is the level within 1 keV of the photon pair's E_x), and the
photon channel's energy from that one E_x (`physics/src/LevelParameters.cpp`, `src/EngineChannelFunctions.cpp`).
Narrow case: the final level's particle channel is open, so the tail beyond the radius is an outgoing continuum
function at E_f instead of a Whittaker function, and N_f has no closed-channel integral; the places are the three
files above that take |E_f − S|, the observable unchanged. Broad case: E_γ varies across the final state's width,
so the observable is folded over E_f, inside the observable (`src/GenMatrixFunc.cpp`) or as a new tree kind of the
prediction. Today a final state above its threshold is not refused (section 9). Not known without trying: whether
the pathway table and the capture masks need a change for an open final channel.

**(c) Secondary γ-ray transitions as fitted data.** What exists: the unobserved-primary, observed-secondary path
(a data set's field `secondary_gamma` with L, I_c and δ; `experiment/src/Build.cpp`, the R_k weighting in
`GenMatrixFunc::LegendreCoefficientsOf`; `tests/secondary_gamma`), and the `gamma_decay` record (Γ_γ, branchings
and mixing ratios of a level as parameters of kind `RadiativeDecay`; `tests/gamma_decay`). δ of `secondary_gamma` is
a fixed number of each data set today. Fitting it is a parameter kind (one row, its slot in
`WorkingCopy::FillMnParams`) with its derivative through R_k; a δ shared by several data sets would be a record
named by them, as an apparatus record is. A γ line fed through several final states (a capture cascade) is a sum
over final states times branchings: the `sum` tree adds reactions, and whether a term can carry a fitted branching
of a bound level is not known without trying.

## 7. What `dev` has that this branch has not, and the way in

From the plan of the join (outside the repository; sections 1, 4–6), against `origin/dev` at `3e502625`: 65 commits
beyond the branch point, 51 of them not merges.

- **Done here:** section 3, rows 4 and 5; dev's new blocks (`<channelRadii>`, `<isospin>`, spin-operator data types)
  refused by the converter and the reader, never read past; the UDR flight path in metres; the scattering length;
  capture with two compound nuclei (right here by the plan's run of dev's models; the test not yet ported).
- **In work:** per-channel radii; general spin observables (as one quantity with an operator field); the
  time-of-flight resolution at low energy, with dev's external-capture tests; the EXFOR user agent.
- **Not started:** isospin-labelled levels across compound nuclei (the largest: 30–60 h by the plan); dev's
  gradient-based sampler `pyazr.sampling` (waits for the PI; a NUTS script on `Objective` is in the uncertainty
  skill); dev's skill notes.

A plain merge of `3e502625`, measured: 93 conflicted files (45 with 125 hunks, 48 modified on one side and deleted
on the other), 21 merged silently, at least eight of them wrongly; most of dev's new engine code is in `CNuc`,
`EData`, `EPoint` and pyazr 2, which this branch replaced. The way in, in three stages: each feature re-written in
its named places, then a merge commit of a commit of dev's line records dev as merged.

**The one request.** Please write the next features (rmfit's merge, capture into unbound states, secondary γ rays)
against this branch and its recipes, or at least against pyazr 3.0's names (`Session`, `Objective`, `Project`,
`Parameters`) and format-3 records. Each week of separation adds 40–80 hours of porting at dev's present pace.

## 8. What is asked of the owner (each with this branch's proposal)

- **The Coulomb repair's home**: a listed local change in `physics/coulomb/coulsm` (as now) or pycoulomb.
- **The coincidence gap**: 1e-6 MeV with the mean S′ (section 3, row 4); proposed for `dev` too.
- **The covariance written to ENDF MF=32**: HESSE in a formal-basis session (as built); the alternative, the fit's
  covariance carried through the transformation, differs on T1 by up to 3.9 σ_i σ_j.
- **The ten switches the history page lists as used by no test**: seven now have an option and a case in
  `tests/switch_options`, `AZR_DEBUG_GRID` is read by `tests/gamma_decay`; proposed: drop the two print-only
  switches no test reads (`AZURE_GRAD_DEBUG`, `AZURE_CHANNEL_DEBUG`) unless you use them.
- **`rmfit-campaign`**: follow `dev` (the skill lives with rmfit) and drop the ported copy here.
- **The logo**: the original on the web page and in the packages (as now), or `dev`'s new one.
- **The format of dev's new blocks**: record kinds `channel_radius` and `isospin_*`, the quantity `spin_observable`
  with an operator field; the converter reads the old blocks.
- **The acceptance rule of a quadrature change** (section 3, row 8): 1.2 times the largest and 1.03 times the r.m.s.
  distance of the old engine to a refined reference.
- **What an energy shift means for three kernels**: the beam-profile window and its detailed-balance norm, a user's
  resolution table and the 5 σ support of an energy-dependent Gaussian still use the unshifted energy.
- **Fitted target parameters in the band**: today an areal density or a stopping scale does not widen it.
- **Minuit2's call counter**: two atomic increments in `external/` (not applied; AZURE2 counts its own calls).
- **Two constants kept so that no result moves**: the electron mass of the Fermi function (CODATA 2002) and the
  proton mass of the M1 capture integral (`include/Constants.h`).
- **Publishing and CI**: the Coulomb table's address and hash (`pyazr/coulomb_table.py`, `TABLE_URL`, empty); one
  run of `.github/workflows/wheels.yml` (never run); a Windows wheel (MinGW or MSVC); SRIM's two data files under
  their authors' terms (`stopping/data/srim/README.md`); `.github/workflows/build.yml` still passes `USE_READLINE`.

## 9. Known defects and limits

Each checked against the start commit of this page, by a run (R) or by reading the code (C):

- A capture final state above its particle threshold is not refused (R: T1 with its 3/2⁻ final state moved to
  2.5 MeV, S_p = 1.943 MeV: `--validate` says valid; points below E_f give NaN and χ² 1e100 without a message;
  above it the external part takes |E_f − S| as a binding energy, C).
- A project file cannot hold two levels of one J^π at the identical energy: their lines are read as one level, with
  a warning (R). Levels within 1e-6 MeV, and a fit that brings two together, work.
- A repeated top-level section: the run reads the first and warns; `--validate` passes a repeated empty section (R).
- `--max-calls` acts on Minuit2's own counter, approximate by a few calls under threads; the `param.fit` snapshot
  depends on the thread schedule, every other output does not (C; `fit_threads.md` sections 3 and 4).
- Levenberg-Marquardt and GSL's trust region ignore `--refreeze-covariance` without a word (C).
- A `<parameterSettings>` row names a level by its running number, which moves when a level is added (C).
- An options file is applied after the command line: its calculation options win over the command line's (C).
- `pyazr.Session` has no threads argument: a session uses the process's OpenMP threads (C).
- The web page was tested in a headless Chrome and used by the PI; Firefox and Safari were not run.
