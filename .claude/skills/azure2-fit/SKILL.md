---
name: azure2-fit
description: Fit an AZURE2 R-matrix model to data - choose free and fixed parameters, seed widths, stage fits, set normalization and energy-shift nuisances and priors, run CLI MIGRAD or a pyazr least-squares fit, judge convergence, add/remove/deactivate levels and background poles, reuse or convert param.sav files, drop a data set, snapshot a fit with save_fit. Use for any fit, refit or model edit. Use azure2-check to judge a result, azure2-uncertainty for errors, azure2-calculate for observables, rmfit-campaign for campaigns.
---

# Fitting an AZURE2 model

A fit adjusts the level parameters (energies E_λ in MeV, reduced-width amplitudes γ_λc in MeV^½, or the observed widths
and ANCs they map to), the data normalizations n and the energy shifts ΔE of a project so that the objective, the data
χ² plus the penalties of the nuisance parameters, is minimal. The result is the `.azr` together with its parameter file
(`param.sav`): the `.azr` alone is never the fit. This skill gives the procedure (free and fixed parameters, starting
values, stages, nuisances and priors, running and convergence) and the rules for editing levels, background poles and
data sets without corrupting the parameter files. Behaviour reference: `docs/source/user_guide/fitting.rst`,
`docs/source/user_guide/chi_squared.rst`, `docs/source/user_guide/pyazr.rst`.

## When to use, and the siblings

- Use this skill to set up, run, restart or extend a fit, to edit the level scheme or the data-set selection of a
  fitted model, and to carry parameter files across such edits.
- A new project, new particle pairs, the first level scheme, EXFOR data: `azure2-setup`.
- Observables, extrapolation, S(0), rates, decomposition, CLI menus and the `*_rwa` conventions: `azure2-calculate`.
- Whether a fit can be trusted (pulls, residual structure, χ²/ν, Wigner limits, radius and background-pole
  sensitivity, saturated poles): `azure2-check`.
- Parameter uncertainties, bands, MINOS, MCMC, data covariance: `azure2-uncertainty`.
- Target, resolution and Doppler nuisances (density, thickness, temperature): `azure2-targets`.
- Missing or misassigned levels (gated level-add tests): `azure2-resonance-search`.
- Bases (formal, Brune, Park) and their conversions: `azure2-exchange`. THM fits: `azure2-thm`.
- Campaigns of polishes, sign rounds and level rounds on a cluster: `rmfit-campaign`.

## 1. Before the fit

1. **Know which files make the model.** The `.azr` holds the nominal norms and shifts; the parameter file holds the
   fitted values (`references/normalizations-and-shifts.md`). Old-format parameter files (`j=<J>_la=<L>_energy`
   names) load nothing but their `segment_*` lines; convert them first (`references/legacy-params.md`).
2. **Load parameters by name.** In pyazr, `x = m.load_sav(path)`; never index a parameter file by row
   (`references/parameters-and-sav.md`).
3. **Audit the nominals.** `dataNorm` should be the data's own scale (1.0, or the published factor), not a previous
   fit's value; a stale `energyShift` becomes the start of a freed shift.
4. **Check the start.** A mode-1 calculation from the exact files and flags the fit will use must reproduce the
   expected χ² (and the first χ² the fit prints must equal it). A start far above the baseline is a problem even when
   it seems explained (`azure2-check`).
5. **For a new reaction**, build starting parameters level by level before any global fit
   (`references/first-fits.md`): fix the scale with one data set, seed the lowest-l channel of every pair on every
   level, put every channel that opens inside the fit range into the model, polish from low energy up, then add levels.

## 2. Free and fixed

- **A width of exactly 0 is fixed, whatever its fix flag says.** `CNuc.cpp` adds every channel with a Minuit step of
  0.1 × value and fixes it when the value is 0, so a channel written as `0` with `chanFix=0` stays at zero for the whole
  fit, and a token seed of 10⁻¹² has a step of 10⁻¹³ and never moves either. Freeing a width is not enough: give it a
  seed of plausible size. The symptom is a free-parameter count that does not change when a channel is freed (freeing
  six E_x = 40 MeV capture background poles left a model at 366 free parameters; seeding them at 0.01 eV made it 372).
  This holds for particle and photon channels and for `width=0` channels passed to `Project.add_level`.
- **An amplitude at exactly zero never moves under a gradient method** even when it is free and nonzero-seeded
  elsewhere: χ² depends on it through γ². Seed in units of the Wigner amplitude (0.01–0.05 γ_W), never as a physical
  width, and audit closed channels just above the levels (`references/pyazr-fitting.md`, "Traps").
- **Fixed level parameters are fixed in reduced-width space**: their observed widths and ANCs move when the free widths
  of the same level or its energy move (`references/editing-levels.md`).
- **A norm that is not varied comes from the `.azr`; a listed energy shift comes from the parameter file, varied or
  not.** To freeze norms at a fit's values, write them into `dataNorm` and clear `varyNorm`; when a fixed shift changes
  in the `.azr`, change or delete it in every parameter file the run reads (`references/normalizations-and-shifts.md`).
- Background poles: their energies are fixed; for a pole feeding one reaction through one entrance and one exit
  channel, fix one of the two widths (`references/background-poles.md`).

## 3. The objective, nuisances and priors

- **The objective is χ² + penalties.** Per segment with `normErr` > 0 (a percentage):
  ((n − n₀)/(n₀ · normErr/100))², measured from the `.azr` nominal n₀. Per free energy shift: (ΔE/σ_ΔE)², measured
  from zero (also for `sqrtshift` b and each `shiftfn` term). In pyazr this is `m.objective()(x)`; it equals
  `Total-Chi-Squared` + `Total-Norm-Chi-Squared` of `chiSquared.out`. Minimizing the bare data χ²
  (`m.evaluate(x).chi2`) lets the normalizations absorb every discrepancy.
- **Normalizations.** Leave a set's norm free with its quoted systematic error as the penalty width. A large,
  angle-dependent normalization is a symptom, not a result (`azure2-check`). With a free, unpenalized norm a fit can
  effectively drop a data set while it still adds leverage: anchor norms that carry an absolute scale.
- **Energy calibration.** Scan the frozen curve first, choose the form (constant, √E, linear) by the apparatus, centre a
  `shiftfn` inside the data, seed freed shifts at the scan optimum, and share one calibration per experiment with
  `shiftgroup`/`shiftfngroup` (`references/normalizations-and-shifts.md`).
- **Nuisance priors on level parameters** (`<parameterSettings>`, nuisance 1) add ((p − p₀)/σ)² in physical units; on
  builds before 2026-10-27 (expires: 2027-02-01) they hit the wrong parameter (`references/running-fits.md`).
  `<parameterSettings>` written by the Qt GUI of released builds use the engine's level numbering (`references/parameters-and-sav.md`).
- **Target and resolution nuisances** (fitted density, thickness, Doppler temperature, screening): `azure2-targets`.
- **Data covariance** (correlated errors, `covariance <file>`): `azure2-uncertainty`.
- **Alternative objectives** (per-data-set reduced χ², Sivia's outlier-tolerant loss) are campaign options of rmfit
  (`references/first-fits.md`); results under different objectives are not comparable.

## 4. Staged fits and multiple starts

1. **Stage the free set by what the data constrain**, warm-starting each stage: (1) the particle widths of the pair that
   carries the dominant reaction; (2) the other particle pairs; (3) the γ widths at fixed particle widths; (4) the
   background poles; (5) everything, including level energies and normalizations. Run the sequence two or three
   times.
2. **Local refits.** Free a few level parameters (and the norms they interact with) with everything else held; freeze
   held norms by writing them into the `.azr` (section 2).
3. **Multiple starts for interfering amplitudes.** Same-J^π levels, background poles against resonances, internal
   against external capture: the objective has several branches. Test both relative signs; seed from a converged fit of
   the same or a similar reaction, signs included; when attempts keep reverting toward the same large value, seed there
   instead of tightening the guard (`references/interference.md`).
4. **Guard or runaway.** A fit stopped at a bound that terminated on `gtol`/`xtol`/`ftol` found a real minimum past the
   bound; one that exhausted its evaluations or sits at many multiples of the bound ran away
   (`references/interference.md`).
5. **Compare like with like.** A structure test and its control use the same minimizer (preferably `--use-gradient`),
   the same budget and the same objective.

## 5. Running the fit

**CLI (reference implementation).** Mode 2 fits; mode 4 fits and extrapolates. The flags (`--use-brune`,
`--use-park`, `--ignore-externals`, `--use-gradient`, ...) are CLI-only: copy them from the reference run's job
script. The menu answers and the order of the uncertainty-band and parameter-file prompts are in `azure2-calculate`;
confirm `Reading User Parameter File...` in the log (`Creating New param.par File...` means the parameter file was
ignored).

```bash
printf '2\nn\nparam.sav\n\n7\n' | "$AZURE2_BINARY" --no-gui --no-readline --use-brune --use-gradient model.azr
```

**pyazr (your minimizer).** pyazr gives `f = m.objective()`: `f(x)`, `f.value_and_grad`, and `f.residuals_and_jacobian` (the whole
Jacobian for about two forward evaluations; the rows carry the priors). A hand-written fit must minimize `f`, respect the Wigner limits as bounds if they matter, drop zero Jacobian columns, and
catch exceptions from `f.jacobian` (`references/pyazr-fitting.md`).

**Convergence.** Read `output/fitStatus.out` (`status converged|call_limit|invalid`, EDM, calls) or the final
`MIGRAD converged` / `WARNING: MIGRAD did NOT converge` line. `chiSquared.out` and `param.sav` are written either way.
A plateau of a few hundred evaluations on a model with hundreds of free parameters is not a minimum
(`references/running-fits.md`).

**Long fits.** Reserve memory per slot (each of max(4, threads) + 1 model copies holds every sub-point), restart from
`param.fit` after a kill, and verify the snapshot with mode 1 before trusting a restart (`references/running-fits.md`).

**Threads.** Results do not depend on the thread count. With several worker processes, keep processes × threads
within the slots requested.

## 6. Editing levels and background poles

- **Write a new `.azr` and open a fresh session**; `pyazr.Project` edits `<levels>` and the records through the
  engine and keeps every other byte (`references/editing-levels.md`).
- **All levels of one J^π share one channel set.** `add_level` clones an existing group's channels.
- **Every structural edit renames parameters.** Removing, deactivating or adding a level renames every later
  `energy_<N>` / `width_<N>_<c>`; a copied `param.sav` puts values on the wrong levels without an error. Bake the fit
  first (`save_fit`) and reseed only the `segment_*` lines, or remap by physical identity and verify every pair
  (`references/parameters-and-sav.md`).
- **Restrict a background pole to the channels that need it**: widths only in the data set's channels, every other
  channel at exactly zero (`references/background-poles.md`).
- **Frozen probes** (`without_level`, `only_level`) change no parameter vector; a persisted "level off" is
  `deactivate_level`.
- **Same-J^π levels closer than 1 keV** stay two levels only with different level IDs
  (`references/parameters-and-sav.md`).

## 7. Data sets in and out

- **Deactivate, never delete** a `<segmentsData>` line: norms, shifts, `<targetInt>` key lists and every output table
  are keyed by line position. Make the drop reach every `.azr` variant that is run, and verify that the total χ² falls
  by exactly that segment's χ² (`references/drop-dataset.md`).
- Adding or removing data segments, or changing the test grids, invalidates `output/intEC.dat` and
  `output/intEC.extrap`: delete them or use a new output directory (`references/editing-levels.md`).

## 8. After the fit: snapshot and hand-over

- `m.save_fit(path, x)` writes the `.azr` and its `param.sav` and verifies the round trip; it leaves the nominal norms
  and shifts in `<segmentsData>` (`references/editing-levels.md`).
- Check a snapshot with mode 1 and its `.sav` (must reproduce `chiSquared.out` including `Total-Norm-Chi-Squared` to
  about 10⁻⁴), and a bake of `<levels>` with a blank parameter file. Near a destructive-interference minimum, check the
  per-segment χ² of the reopened file, not only the verify.
- Say which home (`.azr` nominal or `.sav` fitted) every norm and shift you hand over came from, and record in the
  readme: files, flags, basis, objective, convergence status, data sets off.
- Judge the result before reporting it (`azure2-check`).

## Pitfalls, short

- Two fits are comparable only with the same data sets, objective, basis, flags, grid placement and build.
- A fit's χ² on a build with different target-grid placement is a different number; re-evaluate the end point in a
  current session before comparing (`references/editing-levels.md`).
- Do not use a pyazr session on a running fit's output directory: it can rewrite the fit's caches
  (`azure2-calculate`).
- `m.datasets[i]` is not the i-th result entry when inactive segments exist: use `m.active_datasets`.
- A `Parameter`'s `pair` is the engine's pair number, not the file's pair key.
- Fitting a derived or extrapolated point together with the data it came from double-counts them (`azure2-check`).

- **A narrow level more than about 1 keV from its peak has no gradient pull** and never comes back; it drifts when
  the objective is dominated by another data set (one 0.4 keV-wide level slid 180 keV while an unnormalised set
  dominated). Give narrow levels an energy window of a few times max(Γ, resolution) and normalise the data first.
  MINUIT's initial steps (0.1 E for a level energy, 0.1 γ for a width) do not prevent fitting such a level: from
  every start within ±4 keV of the peak it converged to the same E and Γ.
- **A free level energy must stay a few keV from every particle threshold.** The fit objective is finite at a
  threshold but the Brune transform is not (log Γ of a zero-energy channel), so the failure appears only when the
  physical parameters are written. Keep ~3 keV clear. Near a threshold the Brune transform also caps a width: at
  E_x = 2.37 MeV in ¹³N an l = 0 width cannot exceed about 45 keV whatever amplitude is given, so put a broad test
  level where the penetrability allows it.

## Pointers

- References: `references/first-fits.md`, `references/parameters-and-sav.md`, `references/legacy-params.md`,
  `references/normalizations-and-shifts.md`, `references/pyazr-fitting.md`, `references/running-fits.md`,
  `references/editing-levels.md`, `references/background-poles.md`, `references/interference.md`,
  `references/formalism.md`, `references/drop-dataset.md`.
- Script: `scripts/convert_sav_names.py` (old-format parameter names).
- Examples: `pyazr/examples/edit_model.py`, `pyazr/examples/save_fit_to_azr.py`, `pyazr/examples/fit_emcee.py`,
  `pyazr/examples/per_dataset_chi2.py`, `pyazr/examples/deactivate_level.py`.
- Docs: `docs/source/user_guide/fitting.rst`, `docs/source/user_guide/chi_squared.rst`,
  `docs/source/user_guide/levels_channels.rst`, `docs/source/user_guide/pyazr.rst`,
  `docs/source/reference/output_files.rst`, `docs/source/reference/command_line.rst`.
