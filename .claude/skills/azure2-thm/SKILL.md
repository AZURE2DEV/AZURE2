---
name: azure2-thm
description: Step-by-step procedure for a complete Trojan Horse Method (THM) evaluation with AZURE2, from published THM spectra and direct data to model-averaged results with a stated model dependence. Covers vetting what published THM points are (HOES yields vs penetrability-corrected on-shell-equivalent values, carrier, binding energy, kinematic-factor convention, resolution, subtracted background, digitised error bars), building the project (THM opt-in, THM segments, experiment lines, binding energy, prior centres), choosing the baseline model (vertex, kinematics, line shape, distortion R(E), DW vertex, spectator window and acceptance, fixed angle, backgrounds, Brune vs Park), staged joint THM + direct fits, deriving strengths through the anchor formula, the model-averaging driver scripts/thm_model_average.py with Akaike/flat weights and the robustness rule, cross-checks, and reporting. Use for any THM fit, HOES data, modified R-matrix, Trojan Horse analysis, THM-derived strengths, S factors or rates, or the model (theory) dependence of THM results. Generic AZURE2 mechanics are in azure2-eval; new compound nuclei in azr-project-builder.
---

# THM evaluation with AZURE2

The procedure for turning published THM material into an evaluated, model-averaged
result. Mechanics of segments and files: **`azure2-eval`** ("Trojan Horse (THM)
segments"). Equations and option reference:
`docs/source/theory/thm_implementation.rst` (cited below as *impl*, "Section").
Level schemes and formalism: **`r-matrix-analysis`**. Data fetching: **`nds-explorer`**.
New compound nucleus: **`azr-project-builder`**. Detailed checklists:
`references/` in this skill folder.

What THM data can and cannot give (the lesson of every cross-check so far):

- A THM spectrum with a free scale fixes **resonance energies, line shapes, ratios
  γ_x²|M_l|² between levels, and ratios of exit widths** (several channels, one scale).
- It never fixes **absolute widths**. The scale comes from overlapping direct data or
  from one directly measured strength (the anchor).
- When direct data overlap, reaction-model choices move results by ≤ 3 %. When
  strengths come from peak areas against an anchor of different l or far away in ρ,
  they move by 20 % to factors of 100. Quote that spread with every number.

## When to use
- A fit that includes THM points (HOES or penetrability-corrected), alone or with direct data.
- A strength, S(0), rate or ANC derived from a THM spectrum.
- Reproducing or checking a published THM analysis.
- Any question of how model-dependent a THM result is.

Prerequisites:

1. AZURE2 built from the `thm` branch, with **both** targets current: the CLI
   (`build/src/AZURE2`) and the pyazr module `_azure2`. They are separate CMake
   targets; `make AZURE2` does not rebuild `_azure2`, and a stale module reproduces
   engine bugs that are already fixed. Build with `make -j2` at most on small machines.
2. Python with numpy, scipy and mpmath (pyazr's core dependencies); `import pyazr`
   from the repository root picks up `_azure2` from `pyazr/`.
3. Memory: a realistic THM session takes 0.2–2 GB. On a machine with ~3 GB and no swap,
   run **one engine process at a time** (section "Memory and process discipline").
4. Keep per-evaluation notes, scripts and fit outputs **outside the repository**
   (a sibling working directory). Only generic features, tests and examples go in the repo.

## The evaluation path

Eight stages. Do not start a stage before the previous checkpoint holds.

### Stage 1. Gather and vet the data

The reading of the published points can matter more than any reaction-model choice.
Fill in the per-dataset record of `references/provenance-checklist.md` first.

1. **Decide HOES or on-shell-equivalent.** This decides the observable of every later fit.
   - *HOES yield*: d³σ/(KF|φ|²) in arbitrary units, no penetrability. It stays finite
     or rises toward threshold, can have points below threshold, and narrow
     high-l levels show peaks comparable with broad ones.
   - *On-shell-equivalent*: penetrability restored (often per l, each l with its own
     normalisation to direct data). It falls like e^{−2πη} toward threshold and usually
     carries absolute units.
   - *Already fitted*: an S factor that went through an R-matrix or polynomial fit.
   - The answer is in the data-reduction section of the primary paper, not in the file
     or in EXFOR.
2. **Record the reaction.** Note the carrier a = x + s, the spectator, B_xs, the beam
   and its energy, normal or inverse kinematics, the spectator-momentum cut, the
   angular acceptance, the resolution (σ or FWHM, which frame), the energy frame of
   the tabulated points, and the background the authors subtracted (its formula).
3. **Record what was divided out.** Possible conventions: the full three-body KF,
   λ₃/λ₂, |φ|² only, or a fit with La Cognata's working formula. This sets
   `kinematics=` (Stage 3).
4. **Check digitised points and errors** against the figure:
   - Prefer PDF vector paths. For a tick-calibrated raster, check every major tick to ≤ 1 px.
   - Compare digitised bars with the visible ones. A bar hidden by the marker is
     bounded by the marker radius: cap it there, do not invent it.
   - Use published tables whenever they exist, even if the paper says "not public".
   - Record the digitisation error.
5. **Check the scatter.** Fit a smooth model (for example one Gaussian per peak) to the
   THM points with the published bars. If χ²/ν ≫ 1, the bars are smaller than the
   scatter (statistical only). Note the factor; Stage 4 uses it.
6. **Collect the direct data.** Use EXFOR (`nds-explorer`, `pyazr.nds`) and the
   published tables. Record units, lab or c.m. frame, systematic (norm) errors, and
   whether per-point errors are given or assigned. Put the source in the data file's
   `#` header.

Checkpoint: every THM dataset has a written reading (HOES / on-shell /
already fitted), carrier, B, frame, KF convention, resolution and background, each with
its source.

Pitfalls met:

- **7Li(p,α) Tumino 2006 and 6Li(d,α) Pizzone 2011 are on-shell-equivalent**
  (penetrability restored, normalised to direct data), with a ³He carrier
  (B = 5.4935 MeV). They were once fitted as HOES with B = 2.2246 / 1.4735 MeV. The
  HOES fit pulled the direct norms to 0.4–0.6 and halved the 7Li S(0) (25–37 against
  63.5 keV b), with −10 to −14σ residuals below 0.4 MeV. For 6Li the χ² even
  preferred the wrong reading: the paper text decides, not χ².
- 15N(p,α₀) La Cognata 2007 Table 3 is on-shell (fitted scale 0.98). 17O(n,α)
  Guardo 2017 is genuine HOES: it has points below threshold, a visible f-wave 5⁻ peak,
  and "arb. units". 18O, 19F and 12C+12C spectra are HOES.
- Paper misprints: Tumino 2006 quotes the deuteron's Hulthén parameters for ³He and
  "5.85 MeV" for 5.49 MeV. A momentum-window limit derived with the wrong
  (inverse) kinematics removed 48 of 66 points that the real kinematics reach.
- Frames: the 17O energies are lab (χ² 13.0 against 17.0 as c.m.; the source figure
  confirms it). The 12C+12C files have E_lab = 2 E_cm. A table may mix frames
  (Guardo Table II). An EXFOR EN-CM column labelled MeV held keV (15N).
- Digitisation: one digitised JUNA point was 17 % low against the published Table I.
  The 19F hidden bars were 30–45 % optimistic in the valleys. Bars ≈ cap size need
  ~5 % added in quadrature.
- Error bars smaller than the scatter: the 19F published bars give χ²/ν ≈ 9 for the
  authors' own Gaussian model. The quoted χ²/ν = 2.0 needs ~15 % errors.
- Background provenance: an example THM file once carried a background tuned to one
  joint fit (plane wave, 5.1 fm), not the authors' line. It biased every other variant
  (6.1 fm: THM χ² 2866 against 540 with the authors' background). Digitise the authors'
  background, subtract it, and let `background=linear` absorb the rest.
- Identical exit particles (α+α): data that count both α need ×½. A free norm far from 1
  (Cassagnou ~1.9) can be a disagreement between datasets, not a double count.

### Stage 2. Build the project

Start with **`azr-project-builder`** (pairs, levels, channels) or an existing `.azr`.
Then:

1. **Enable THM.** In the GUI, tick *Configure > Runtime Options > Use Trojan Horse
   Method (THM)*. It is off for a new project and turns on by itself for a project
   with THM content; it is not stored in the file. All THM editing is in
   *Configure > THM Workspace...* (pages Model, Experiments, Channels, Diagnostics).
   pyazr and the CLI need no switch.
2. **Segments.**
   - HOES data become a THM segment: observable code + 10 (10 angle-integrated,
     11 differential). In the GUI, tick *THM* in the segment dialog; in pyazr,
     `AzrModel.add_data_segment(..., thm=True)`. Free norm, no norm error: AZURE2
     profiles it analytically, n* = S_mm/S_md, with no parameter and no penalty.
   - On-shell-equivalent points become an ordinary segment (σ or dσ/dΩ, converted
     from S if needed) with a free norm. None of the THM options acts on them.
3. **Binding energy.** B_xs of the carrier goes in field 32 of **every** channel line
   of the THM entrance pair (THM Workspace > Channels; pyazr `AzrChannel.binding_energy`).
   The Channels page flags a B that differs from the experiment's masses by > 1 keV.
4. **Experiment lines.** Group the segments of one measurement (several exit channels,
   angular bins, runs) in `<thm>`:
   `experiment[E1] segments=1-4 beam=14N target=12C spectator=d Ebeam=30`
   - One shared profiled norm makes the relative channel heights carry exit-width
     ratios. For 12C+12C, four free scales came out within ±5 % of their mean.
   - `beam/target/spectator/Ebeam` (all four or none) are needed by `lineshape`,
     `ps`, `distortion=coulomb|optical` and `vertexModel=dw`.
   - `background=const|linear|quadratic` is a smooth incoherent term, profiled linearly
     with the norm. Results are in `output/thm_experiments.out` and
     `session.thm_background(name)`.
   - pyazr: `AzrModel.set_thm_experiment(name, segments, ...)`.
5. **Resolution.** Use one `<targetInt>` Gaussian per THM segment. Its σ is a
   **lab** standard deviation: convert c.m. σ by ×(m₁+m₂)/m₂, and convert FWHM to σ
   (÷2.355). Example: 17 keV c.m. gives ≈ 0.01795 MeV for p+18O. Without a stated value,
   start at σ ≈ 30 keV c.m. and fit or scan it. pyazr: `add_target_effect(segments,
   gaussian_sigma=...)`.
6. **Direct-data norms.** Give each direct segment a norm error equal to its quoted
   systematic. Centre the prior on the nominal scale with
   `segment_N_norm prior_centre 1` in `<parameterSettings>` (GUI: Fitting tab,
   *Prior Centre* column; pyazr `AzrModel.set_prior_centres`). Without the row, the
   norm column is both start value and prior centre. A saved fit then re-centres every
   later refit on its own result (15N S(0) 70.7 → 74).
7. **Level scheme.** Take it from ENSDF and direct R-matrix analyses, not from the THM
   paper alone. Include every level the THM spectrum shows, narrow ones too, and check
   that no data window cuts a level.
   - Example: the 18O 3/2⁻ 597.6 keV level (Γ ≈ 2 keV) shows as a step at
     0.58–0.62 MeV in the 2010 spectrum. With LUNA widths and the DW vertex it lowers
     the joint χ² from 481 to 330 with no new parameter.
   - Identical 0⁺ bosons: only even-J natural-parity levels fuse. THM-only levels
     (odd J, unnatural parity) need the workaround in `references/model-choices.md`.

Checkpoint:
- The project runs on the CLI (`printf '1\n\n\n7\n' | build/src/AZURE2 --no-gui
  --no-readline x.azr`).
- `<thm>` is accepted: AZURE2 exits with `ERROR: <thm> ...` otherwise.
- B, σ and the norm priors are as recorded in Stage 1.
- THM Workspace > Diagnostics shows |M_l|² and the HOES/on-shell shapes as expected.

### Stage 3. Baseline model: choices and why

Each choice gets a default, the condition to change it, and its place on the variant
grid of Stage 6. Details and effect sizes: `references/model-choices.md`.

| choice | key | default / when to change |
|---|---|---|
| entrance l coherence | `entranceL=` | `incoherent`: exact for 4π, spin-summed data. For a restricted angular range use `theta=` (below), not `coherent` |
| vertex boundary | `vertex=` | `constant` (B_c = S_c(E₁) of the lowest level of each Jπ): representation-invariant. `onshell` (S+iP) as a variant. `perlevel` only to reproduce pre-Sep-2026 numbers (18O band 23 % vs 5 % peak rms) |
| kinematic factor | `kinematics=` | Match what the data divided by: full three-body KF → `kf3body`; λ₃/λ₂ → `lambda32`; \|φ\|² only → `triple`; La Cognata's formula or unknown → `lacognata`. λ₃/λ₂ leaves 1/k_i (×7 over 0.02–1 MeV) |
| entrance vertex model | `vertexModel=` | `pw` default. Make `dw` (with `distortion=coulomb\|optical`, `opticalAA=`, `opticalSF=`) a variant whenever peak areas of different l, or levels far apart in ρ, become strengths |
| Coulomb integral C_l | `coulombIntegral=` | `0` (no published analysis uses it). `1` only with a plane a+A wave (refused with `dw` or a distorted R(E)); ×3–4 run time |
| Coulomb line shape N_C | `lineshape=on` | Charged spectator only (N_C ≡ 1 for a neutron). Off in published fits; make it a variant (17O, p spectator: +10–17 % on strengths). Refused with `cbackground=` |
| distortion R(E) | `distortion=` | For a charged spectator below its s+F barrier (heavy systems). It multiplies the pw model. Not applied with `dw` (the DW vertex contains it). `distortionRatio=dw` refused with a `ps` window |
| spectator window | `ps=hulthen:0-40`, `gauss:...`, `table:<file>` | Average over the accepted \|p_s\| when the vertex is not smooth over it (near nodes). Excludes `spectatorEnergy` for that pair |
| acceptance | `spectatorAngles=[cm:]a-b`, `table:<file>` | Accepted spectator directions; with `ps` one acceptance (19F: `spectatorAngles=cm:135-180`) |
| fixed angle | `theta=50-70` | Data in a restricted c.m. angular range of the exit pair (7Li θ_cm 50–70°: shape 14 % rms) |
| backgrounds | `background=`, background pole, `cbackground=` | Physics also in direct data → background pole (broad level, same Jπ). Non-quasi-free residue → `background=` (incoherent). `cbackground=` (interfering, THM-only) is a diagnostic, not a model |
| parametrization | CLI `--use-park` / `--no-brune`, pyazr `use_park=True` | Brune (default). Park gives the same THM model (`tests/thm_park`, ≤ 1e-8) and is preferable when single partial widths are fixed, bounded or given priors |

Rules that the engine enforces (*impl*, "Coulomb effects: what each option contains"):
`coulombIntegral=1` is refused with `vertexModel=dw` and with R(E) on a distorted a+A
wave; `lineshape=on` is allowed with R(E) and with `dw`; `ps` is allowed with `dw`
(averaged per node). Do not work around a refusal: it marks double counting.

Checkpoint: write the adopted baseline and the reason for each row. Plot |M_l(E)|² at
the adopted radius and at ±1 fm with the level energies (THM Workspace > Diagnostics,
or `session.thm_vertex(name, E, strict=True)`). A level next to a vertex node warns
that its strength will be model-dependent.

### Stage 4. Fit

AZURE2 computes χ², residuals and Jacobians. pyazr does not fit: the fit loop is
yours (scipy `least_squares` on `residuals()` / `residual_jacobian()`), or the CLI
minimizers. THM Jacobian rows are central differences (the adjoint does not cover
HOES), with the profiled scale included exactly.

1. **Stage the fit.**
   1. Fit direct data only, to learn what they fix.
   2. Fit THM only, to see what shape it wants (absolute widths will be wrong).
   3. Fit jointly, freeing the few parameters both constrain.
   4. Free the rest step by step. Carry the parameter vector between steps.
2. **Joint fit rules.**
   - The THM scale is profiled; the absolute scale comes only from the direct data
     (norm priors) or from penalty rows of measured strengths.
   - Fixed widths in `<levels>` are fixed amplitudes. Under Brune the physical width
     drifts with the other amplitudes; under Park it stays fixed.
   - Large models (12C+12C, 200+ parameters) need priors: Gaussian in ln Γ_c
     (σ = 1), E (σ ≈ 30 keV), and a penalty above the Wigner limit. A width whose
     posterior error is not below half its prior width is not measured; mark it.
3. **Convergence.**
   - Stop on a Δχ² rule (for example < 0.1 over 10 evaluations), not only on an
     evaluation cap. A fit stopped at the cap gives only an approximate covariance.
   - Restart from perturbed points; the same χ² should come back to < 0.01.
   - With penalty rows whose column norms differ by orders of magnitude, use
     `x_scale=1` (driver `--x-scale 1`). `jac` rejected every step on the 19F model.
   - Converge the resolution fold: raise the `<targetInt>` sub-point tokens (for
     example `0.04 5 50` against `0.04 10 200`) until χ² stops moving. Dense grids
     cost memory (a 30-level model: 2.2 GB against 360 MB).
4. **Assess the fit.**
   - Check χ²/N per dataset, not only the total.
   - Direct norms should sit near 1. Norms of 1.25–1.5 signal a wrong THM observable
     or a wrong vertex (18O 597 keV with pw: 1.33).
   - Check θ² against the Wigner limit (`pyazr.widths`), parameters at bounds,
     super-Wigner high-l entrance amplitudes used as HOES shape knobs, and the residual
     pattern near threshold.
5. **Error model.**
   - If the THM bars are smaller than the scatter (Stage 1), scale the statistical
     errors by √(χ²/ν) from the THM points alone and say so.
   - Test alternative error models (10 %, 15 %, c√y): the strengths should not move
     (19F: ≤ 4 %).
   - Bars that are generous (12C+12C THM χ²/N ≈ 0.24) make Akaike weights meaningless
     later.

Checkpoint: a converged adopted fit saved with `save_fit(path, x)`. It writes the
fitted norms, adds `prior_centre` rows and verifies the snapshot. As the last call
of a fit, use `close_session=True`. Running the snapshot alone on the CLI reproduces
the fitted χ² (data + priors).

### Stage 5. Derive the quantities

1. **Widths and strengths.** Take physical widths from the Brune transform
   (`transform_rwa`, `parameters.out`). Keep signs. Never convert a `gammaIsRWA`
   channel by hand. ωγ = g_J Γ_in Γ_out / Γ, × (1+δ₁₂) for identical entrance nuclei.
2. **Anchor formula.** It turns peak areas into strengths (narrow, isolated levels):

   ωγ_i/ωγ_ref = (Y_i/Y_ref) · [C_{l_i}(E_i)/C_{l_ref}(E_ref)] · [K(E_ref)/K(E_i)] · (b_ref/b_i),
   with C_l = P_l/|M_l|².

   Only the vertex ratio |M_{l_ref}(E_ref)/M_{l_i}(E_i)|² carries reaction-model
   dependence.
   - Choose anchors with the **same l and small |Δρ|/d₀**, where ρ = p_xA·a and d₀ is
     the distance to the nearest zero of M_l.
   - 18O: ωγ(20)/ωγ(90) (both l = 2, |Δρ|/d₀ = 0.009) is stable to 1.5 %. ωγ(20) via
     the 144 keV s-wave anchor (near the M₀ zero) varies ×25 over ±1 fm.
   - A peak displaced from E_i is measured at the peak and moved to E_i with P_l. The
     20 keV level moves ×4.6 per keV, so the E_r uncertainty can dominate.
3. **S(E), extrapolations, rates.**
   - Add an extrapolation grid with `add_extrapolation(..., frame="cm")` for c.m.
     energies (the default is lab).
   - Integrate rates numerically on a resonance-resolving grid, or as a narrow-resonance
     sum checked against the integral.
   - State what sets the low-T rate. For 19F below T₉ ≈ 0.05 it is the unconstrained
     Γ_α2(11), not the THM strength.
4. **ANCs.** Sub-threshold levels give ANCs from the fitted widths. A resonance ANC of
   a narrow level adds nothing beyond Γ_p.

Checkpoint: each quantity has a value, a propagated statistical error from the fit
covariance, and its anchor (direct data, a named strength, or a ratio) written next to it.

### Stage 6. Model dependence: refit, average, classify

Every final number is a model average over **full refits** (a re-profiled χ² at fixed
parameters measures goodness of fit only, never parameter shifts).

1. **Choose the axes.** Use only those that apply to the data:
   - always: the entrance channel radius (adopted ± 1 fm);
   - HOES peak-area strengths: `vertexModel` pw/dw (dw with each reasonable pair of
     optical potentials), and `ps` delta / window;
   - charged spectator: `lineshape` on/off; below its barrier also `distortion`
     none/coulomb;
   - data treatment as separate grids: reading (HOES/on-shell), background, anchor.
   - On-shell-equivalent data: only the radius axis acts.
2. **Run the driver.** It fits one variant per fresh subprocess, so no engine state
   leaks between variants. Do not use `--in-process` across radii.

   ```bash
   OMP_NUM_THREADS=1 prlimit --as=1600000000 -- \
   python3 scripts/thm_model_average.py examples/f19_pag_thm/f19_pag_thm.azr \
       --out ../f19_modelavg --radius-pairs 1 4 5 --radii 4.136 5.136 6.136 \
       --vertex-model pw dw --optical ancai06/kd03:extrapolate \
       --ps delta hulthen:0-50 --strength 213=2-@13.057 --strength-in 1 --strength-out 6 \
       --penalty-hook strengths.py:rows --x-scale 1 --max-nfev 40 --rescale best
   ```

   - Run `--dry-run` first: it lists the grid and the known refusals without the engine.
   - `--strength NAME=Jπ@E` takes the **excitation** energy in MeV.
   - `--penalty-hook file.py:func` adds signed rows (direct strengths, priors) to the fit,
     the χ² and the weights. Other quantities go through `--derived-hook`.
   - Other axes: `--vertex constant perlevel onshell`, `--lineshape on off`,
     `--distortion none coulomb optical`, `--pw-distortion`, or a JSON `--spec`.
   - Large prior-held models (12C+12C) need your own fit loop with priors. Average
     those fits with `pyazr.modelavg` (`Variant`, `model_average`, `write_averaged_azr`).
3. **Weights.**
   - Compute both Akaike (`--weights aic`, default) and flat weights.
   - Rescale χ² by s_χ = max(1, χ²/ν of the best variant) (`--rescale best`) whenever
     the best χ²/ν is well above 1. Otherwise Δχ² is inflated and one variant takes all
     the weight for the wrong reason.
   - **ν rule:** ν = N − k. N counts data points (THM + direct) only; penalty and
     prior rows add to χ², not to N. k counts free parameters including profiled
     scales, backgrounds and direct norms. The driver counts this way.
   - 12C+12C counted its 211 prior rows as data once (N = 458). Corrected
     (N = 247, k = 212, ν = 35, s_χ = 2.59), the AIC weights went from 0.87/0.12 to
     0.54/0.25/0.10/0.07.
   - Where AIC weights are not meaningful (generous bars with χ²/N ≪ 1, priors centred
     on one variant), quote flat weights as primary and say so. Treat a variant that
     takes all the weight with suspicion, not as a selection.
4. **Spreads.** Compute m̄ = Σwᵢmᵢ, σ_stat² = Σwᵢσᵢ² and σ_mod² = Σwᵢ(mᵢ − m̄)².
   Quote m̄ ± σ_stat ± σ_mod. Average quantities spanning decades in log₁₀ (dex).
5. **Robustness rule.** Evaluate (i) σ_mod ≤ σ_stat and (ii) σ_mod ≤ 0.2|m̄|
   (≤ 0.08 dex in log₁₀), with Akaike and with flat weights.
   - **robust**: both conditions hold with both weightings.
   - **model-dependent**: neither condition holds with either weighting.
   - **marginal**: everything else.
   - Two further labels: **undetermined** (depends on an input no data constrain) and
     **held by a direct value** (a penalty row fixes it, so it is not a THM result).
6. **Leading ingredient.** Name the axis along which the mean moves most with the other
   axes fixed. Report σ_mod/σ_stat, σ_mod/|m̄|, the range over variants, Δl and
   |Δρ|/d₀ to the anchor, and the spectator charge.

Checkpoint: per quantity, m̄ ± stat ± model with both weightings, the class, the
leading ingredient, and the variant table (χ², N, k, weights). Calibration values from
the cross-checks are in `references/calibration.md`.

### Stage 7. Cross-checks

1. **Engine sanity.** Run `tests/18O_p_a_thm/check.sh build/src/AZURE2`. It reproduces
   the La Cognata et al. ApJ 723 (2010) band: `vertex=constant` ≤ 7 % rms in the peak
   region, `perlevel` ≥ 15 %. Rerun the examples' "file alone" χ² (*impl*, "Examples")
   after any rebuild.
2. **Reproduce the published analysis** where its inputs are public, before replacing
   it. Use one free normalisation over a stated range, the same formula, and an
   independent Python evaluation of that formula. Protocol and verdict wording:
   `references/reproduction.md`.
3. **Pulls against direct data**: (ours − lit)/σ with σ = (stat ⊕ model) ⊕ σ_lit. Use
   independent direct data not in the fit where possible (17O: Koehler & Graff
   0.2–20 keV / model = 0.88 ± 0.07 with the DW vertex).
4. **Consistency of spreads.**
   - σ_stat should not exceed the physical range (a linearised error above any value
     means the data do not fix that width).
   - σ_mod should come from acceptable fits. Flat spreads driven by variants rejected
     at Δχ² ≥ 200 are not model uncertainty.
   - Check that every variant ended converged, not at its evaluation cap.
5. **Data-treatment checks** as separate grids (reading, background, digitisation,
   error model). Keep them out of the reaction-model average unless they are real
   alternatives.

### Stage 8. Report

1. **Per quantity, quote** the value ± stat ± model, the weighting (AIC rescaled or
   flat, and why), the class, the leading ingredient, the anchor, the radii, the l
   values, and the options of the baseline.
2. **State the model dependence plainly.** For example: "model-dependent: the 8 variants
   span 0.16–16 × 10⁻¹⁹ eV, set by the radius through the zero of M₀ near the
   anchor". A model-dependent result is a result, not a failure.
3. **Compare with published analyses constructively.** Locate each difference in a
   step of the analysis: the reading of the points, the area-to-strength conversion,
   the normalisation window, the extrapolation shape, or a table transposition. Phrase
   it as an observation plus what information would resolve it (the list in
   `references/provenance-checklist.md`, "What to ask authors to publish"). Never
   phrase it as an error in the measurement.
4. **Figures.** Data and our fits only (a reproduction figure is the exception); a
   forest plot (ratio to literature, log axis, stat bar, grey model band); |M_l|² vs E
   with nodes at ±1 fm. PRC style, residual panels, the reaction as a bold label inside
   the axes, PDF + PNG; captions state the fold (σ, frame) and χ² per dataset.

## Pitfalls checklist

- [ ] Penetrability-corrected points fitted as HOES (7Li/6Li: S(0) halved through direct norms 0.4–0.6).
- [ ] Wrong carrier or B in field 32 (³He 5.4935, not d 2.2246 or 6Li 1.4735).
- [ ] Wrong `kinematics=`: the KF convention changes relative peak heights, and for several exit channels their α/p weight (×5 for 12C+12C).
- [ ] Lab vs c.m. energy column (17O lab; 12C+12C E_lab = 2E_cm); `<targetInt>` σ is lab and a standard deviation, not FWHM.
- [ ] Digitised bars not checked against visible ones; hidden bars invented rather than capped.
- [ ] Published bars smaller than the scatter, taken at face value (19F χ²/ν ≈ 9).
- [ ] A subtracted background tuned to one model, reused for all variants.
- [ ] A narrow level missing or cut by a data window (18O 597 keV).
- [ ] Norm column used as prior centre after a saved fit: add `prior_centre 1` rows.
- [ ] Fixed Brune amplitudes treated as fixed widths.
- [ ] Fits stopped at the evaluation cap and reported as converged.
- [ ] Resolution fold not converged on eV–keV levels.
- [ ] Strength from a peak area against an anchor of other l, or near a vertex node, quoted without a radius and vertex scan.
- [ ] `vertex=perlevel` used for new work (representation-dependent).
- [ ] `coherent` used as a fixed-angle observable (use `theta=`).
- [ ] Prior or penalty rows counted in N for ν; AIC weights without `--rescale best` at χ²/ν ≫ 1.
- [ ] Several radius variants fitted in one process (`--in-process`). An engine memo once leaked B_c between sessions (fixed in 49b6c30; `tests/pyazr/thm_vertex_memo_sessions_test.py`).
- [ ] Stale `_azure2` after an engine change; stale `output/intEC.dat`/`intEC.extrap` before an extrapolation with external capture (delete them).
- [ ] Old snapshot files run without checking that the "file alone" χ² reproduces the fit.
- [ ] `thm_vertex` read at an energy the window cannot reach (use `strict=True`).
- [ ] Identical exit pair without ×½ in data that count both particles; odd-J levels of an identical-boson entrance pair.
- [ ] Interfering background (`cbackground=`) read as evidence of a mechanism.
- [ ] Rates for heavy systems quoted without the normalisation window, R(E) choice and ≥ 2 radii (12C+12C: plain fit ×14 above the N_C×R fit at 0.5 GK).

## Pointers

| topic | where |
|---|---|
| HOES observable, `<thm>` keys, KF table, HOES vs on-shell | `docs/source/theory/thm_implementation.rst`, "Options (``<thm>`` block)" |
| experiment lines, shared norm, background | same, "THM experiments"; `tests/thm_experiment`, `tests/pyazr/thm_experiment_test.py` |
| interfering background | same, "Coherent background"; `tests/thm_coherent_background` |
| line shape N_C | same, "Coulomb line shape"; `tests/thm_lineshape` |
| spectator window, acceptance | same, "Spectator-momentum window", "Experimental acceptance"; `tests/thm_spectator_window`, `tests/thm_spectator_angles` |
| R(E), DW vertex, optical potentials | same, "Distortion factor R(E)", "Distorted-wave entrance vertex", "Global optical potentials"; `tests/thm_distortion`, `tests/thm_dw_vertex` |
| fixed angle | same, "Fixed-angle observable"; `tests/thm_fixed_angle` |
| Coulomb combination rules | same, "Coulomb effects: what each option contains"; `tests/thm_coulomb_consistency` |
| profiled norm, Jacobian, band | same, "Normalization, gradients and uncertainty bands"; `tests/thm_band` |
| Brune and Park | same, "Brune and Park"; `tests/thm_park` |
| example pins ("file alone") | same, "Examples" |
| model averaging | same, "Model averaging"; `docs/source/user_guide/pyazr.rst`, "Model averaging"; `scripts/thm_model_average.py --help`; `pyazr/modelavg.py`; `tests/pyazr/thm_model_average_test.py` |
| GUI opt-in, THM Workspace | `docs/source/user_guide/configure_menu.rst`, "THM Workspace" |
| prior centres | `docs/source/user_guide/chi_squared.rst` (nominal-norm); `docs/source/user_guide/fitting.rst`; `tests/prior_centre` |
| pyazr THM calls | `docs/source/user_guide/pyazr.rst` (`set_thm_experiment`, `thm_vertex`, `thm_distortion`, `thm_lineshape`, `thm_background`, `save_fit`) |
| option regression, reproduction | `tests/thm_options`, `tests/18O_p_a_thm` (ApJ 723 band) |
| folding, threshold | `tests/thm_narrow_fold`, `tests/thm_threshold_grid`, `tests/thm_energy_shift` |
| identical nuclei | `tests/identical_entrance_reaction` |
| worked examples | `examples/o18_lacognata2008` (narrow, anchor), `examples/o18_lacognata2010` (doublet, DW, authors' background), `examples/o17_guardo2017_fit` (n entrance, p spectator), `examples/f19_pag_thm` (HOES + direct, linear background, direct strengths), `examples/c12c12_tumino2018` (four channels, one experiment, identical nuclei), `examples/li7_tumino2006`, `examples/li6_pizzone2011`, `examples/n15_lacognata2007` (on-shell-equivalent) |
| classic-project differences of the branch | `MERGE_NOTES.md` |
| detailed checklists | `references/provenance-checklist.md`, `references/model-choices.md`, `references/calibration.md`, `references/reproduction.md` |

## Memory and process discipline (small machines)

- Run **one engine process at a time**: the CLI, a pyazr session, the GUI's
  Diagnostics and builds all count. Check `ps` for AZURE2, make or python first.
  Parallel refits (1–2 GB each) have triggered the OOM killer.
- Run every engine job under a cap: `prlimit --as=1600000000 -- <command>`.
- Set `OMP_NUM_THREADS=1` before the first numpy import; the driver does this itself.
  Without it, evaluations under load stalled for tens of minutes.
- Model averaging runs one fresh subprocess per variant (the driver's default), one at
  a time. The driver holds no session itself.
- Close sessions before opening another (`save_fit(..., close_session=True)`); closed
  memory is reused, not returned to the OS.
- Keep `<targetInt>` sub-point density and `ps` node counts as low as convergence
  allows; check one evaluation's memory before a long fit. Run long fits with `nohup`
  and checkpoint the parameter vector so a capped or killed fit can resume.
