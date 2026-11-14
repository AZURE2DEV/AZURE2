---
name: azure2-resonance-search
description: Gated search for missing or misassigned levels in an AZURE2 R-matrix fit. Rank the regions where the fit misses the data (coherent residual runs, data-only peak finder, blurred broad-structure map), list the hypotheses, screen them, test level adds, removals, moves and sign flips refitted against the base (level_test.py: Delta chi2, AIC, BIC, F test, held-out data), and judge on all the data. Use for "is a level missing", residual maps, level adds or removals. Campaigns: rmfit-campaign; level edits: azure2-fit; trust: azure2-check.
---

# Searching for missing or misassigned levels

This skill is the procedure for one question: which structure does the model miss, and is it a level? The input is a
finished fit (an `.azr` with its parameter file and a verified objective χ²) and its data; the output is, per region,
a verdict (accept, reject, or one of the classes in "Verdict") with the evidence the evaluator needs to adopt it.
Running campaigns, polishes and sign rounds is `rmfit-campaign`; adding, removing and freeing levels by hand is
`azure2-fit`; whether an existing level is needed (frozen Δχ²) is `azure2-check`. What a level, channel or
interference pattern means: `azure2-setup` (`references/` there) and `azure2-fit`.

The scripts in `scripts/` use pyazr only (`$AZURE2_SRC` is the repository checkout that holds `pyazr/`, with the
engine built); the campaign commands of §4 and §5 are rmfit's (`rmfit-campaign`). `scripts/selftest.py` runs every
script on `tests/upgrade/T3` with its 2⁻ level at 19.0006 MeV removed and finds it again (14 s at one thread).

```bash
export PYTHONPATH=$AZURE2_SRC OMP_NUM_THREADS=1
S=$AZURE2_SRC/.claude/skills/azure2-resonance-search/scripts
```

## 0. Ground rules

1. **One model, the whole calculation.** Search inside the full fit. To concentrate on an energy region, fix the
   parameters outside it (levelFix / channelFix, or a local-release mask). Do not build a window or partial-data
   model and reconcile it afterwards. A data-only window fit is allowed as a screen for J^π, E_x and Γ (§4c), never
   as a start vector.
2. **Judge on all the data.** The verdict compares the verified objective of the full model with that of a
   same-budget control polish of the incumbent. A candidate that improves a data set other than the one it was aimed
   at still counts. The per-dataset split (ρ, share of the top data set) is reported, not used as a veto.
3. **Broad first.** Fit structure with Γ ≳ 10 keV over the whole range before aligning narrow levels. A narrow spike
   or a line-shape mismatch in one finely stepped data set is end-stage work. Narrow levels added before the broad
   ones grow into the missing broad strength.
4. **Data definition outranks structure.** A wrong sum of exit channels, an unsummed total yield, a bad detector or
   one outlier point has each been worth more than any level add. Rule them out (§2d) before searching.
5. **Only finished fits are starting points** (an export with a verified objective). The evaluator adopts: report,
   do not adopt.

## 1. Locate: rank the regions where the fit misses

Start from the incumbent's point-by-point residuals. The dump uses r = (fit − n·data)/(n·σ), the residual that
enters χ² with the fitted normalization n; "model low" is r < 0.

```bash
# one evaluation of the fitted project (or --sav fit.sav: values loaded by name)
python3 $S/dump_residuals.py fit.azr resid.npz
# seconds; rerun freely
python3 $S/rank_map.py resid.npz regions.txt fit_output/parameters.out fit.azr
```

Check first that Σ r² of the dump equals the fit's data χ² (the script prints both; they differ only where a data
covariance joins points, and then the ranking ignores the correlations). `parameters.out` is the program's output of
a calculation of the fit (`printf '1\n\n\n7\n' | $AZURE2_BINARY --no-readline fit.azr`). On T3 with the 2⁻ level
removed, the first region peaks at 18.962 MeV (truth 19.0006).

How the ranking works:

- Each detector of each segment is one excitation function; a fixed-energy angular distribution is one curve in
  angle. In each curve the maximal runs of n ≥ 3 consecutive same-sign residuals are found.
- A run is scored by its coherent χ², (Σ r)²/n: the part a smooth change of the curve could remove. Scatter,
  including scatter from underestimated error bars (s² = χ²/N of 10–20 is common), alternates in sign and scores
  about 1 per run. Raw χ² would rank the noisiest data set first.
- Each run's score is spread over its E_x extent (at least 10 keV) on a 2 keV grid; the map is smoothed (6 keV) and
  cut into regions at its local minima, ranked by integral. Do not merge overlapping run intervals instead: with
  dense multi-angle data the intervals chain into one region several MeV wide. Data-set credit is by overlap
  fraction. Runs wider than about 0.6 MeV are listed separately as trends (normalization, background, an
  energy-dependent systematic), not levels.
- Each region lists the model's levels inside it (J^π, E_x, total Γ from `parameters.out`, the fixed flag from the
  `.azr`). Look first for a strong region with no level in it (a candidate new level) and for a region whose levels
  are all fixed (the freeze, not the physics, may be the limit). A level several MeV wide inside a region is
  background, not a resonance.
- Read the list also sorted by density (coherent χ² per MeV, floor 20 keV): a sharp 12-point feature in a
  secondary-γ set never outranks a 700-point neutron-total run on totals alone.
- **Curve grouping.** A fixed lab detector's c.m. angle drifts by several degrees with energy. All scripts group
  points into detectors by gap clustering at 2° within a segment (`detector_labels`), never by rounded c.m. angle,
  which fragments every multi-angle segment.

Inside a campaign the equivalent is the candidate step of `rmfit levelscan` (`levelscan.candidate_energies`: runs
with |r| > 3 clustered in 30 keV bins, `--rank density`, `--seg-match`, `--ex-range`). Use it when the next step is
a screen anyway; use the scripts here when the evaluator wants the list first.

Hand the evaluator the list before spending node time: rank, E_x range, centre, coherent χ², data sets, sign (model
low or high), width. The evaluator usually knows several of the regions; ask which and merge the lists.

### 1a. Starting from nothing: the data-only peak finder

With no R-matrix calculation yet, or a model too poor for meaningful residuals, take the first candidate list from
the data alone:

```bash
# per-point table: dump_residuals.py on any .azr that has the pairs and the data segments (only key, name, E_x,
# angle, data and error are used; AZURE2 converts lab -> c.m. -> E_x)
python3 $S/peakfind.py table.npz peaks.txt [--ex-range lo,hi] [--levels parameters.out]
python3 $S/eval_peaks.py ...    # scores a candidate list blind against a level set, with the chance rate
```

Per excitation function (≥ 8 points): a running-median baseline at 0.05, 0.15 and 0.5 MeV; peaks and dips of
z = (y − b)/σ with prominence ≥ 4σ and relative size ≥ 5 % (small error bars alone cannot make a candidate); FWHM at
half prominence; an adjacent peak and dip form one dispersion candidate. All detections enter an E_x density map
(each a unit-area Gaussian of its own width, weight min(z, 10)) whose peaks are the candidates, ranked by the number
of distinct data sets that see them. "(res.)" marks a width at the data's energy step. Recall measured blind against
a chance baseline: `references/locators.md`.

Use it as a recall tool for an overfit-then-eliminate start: put a level at every multi-set candidate, fit, remove
what the data do not need. It finds narrow and intermediate levels (Γ ≲ 50 keV) well above chance; it does not locate
broad levels in a dense spectrum. Build the broad component separately first (integrated and total cross sections,
ENSDF levels, background poles above the range). Once a first fit exists, switch to the residual ranking of §1.
Agreement of the two methods on a level-free candidate is evidence.

**A recovery number means nothing without the chance rate** at the same tolerance: score the same candidates shifted
by 50–300 keV and quote both.

### 1b. Broad structure: blur the data and the fit the same way

Broad levels are hard to see under many narrow ones. Averaging over an energy interval D keeps what is wider than D:
a narrow level (Γ ≪ D) is reduced to its area spread over D, many narrow levels to a smooth average. For a
Lorentzian kernel and the total cross section the average equals the calculation at the complex energy E + iD. The
narrow levels' area does not vanish, so blurred data are compared with the equally blurred fit, never with a curve
of broad levels only.

```bash
# 1. where blurring is valid for this data set: switch off the narrow levels and compare blurred models
python3 $S/make_broad_only.py fit.azr output/parameters.out broad_only.azr --gmax 50   # Γ < 50 keV zeroed, fixed
python3 $S/dump_residuals.py broad_only.azr broad.npz
python3 $S/blur_compare.py full.npz broad.npz blur_compare.txt
# 2. broad-structure map: blurred data (x fitted norm) minus blurred fit, per detector
python3 $S/blur_map.py full.npz blur_map.txt --scales 0.1,0.2 --levels output/parameters.out --azr fit.azr
# 3. figure
python3 $S/blur_figure.py full.npz broad.npz output/parameters.out fig.png 0.1 --panel "set:angle:lo:hi:title"
```

- `blur_compare.py` prints, per data set and D, the rms relative difference full − broad (how much the narrow levels
  matter), blur(full) − blur(broad) (what of them survives the blur: the validity number) and blur(full) − broad.
  Where the validity number is large (≳ 25 %), the narrow levels carry the strength and a blurred deviation can be
  missing narrow area, not a missing broad level. Use only data sets that pass.
- `blur_map.py`: per detector, δ = (blur(n·d) − blur(fit))/|blur(fit)|; a deviation is a one-sign stretch with
  |δ| ≥ 10 % over ≥ D, scored (extent/D) × mean(δ²/(0.05² + σ_b²)), a χ²-like score with a 5 % floor per
  independent blur width so that point count and a fitted normalization offset cannot win. D = 100 keV resolved
  regions where D = 200 keV merged a 6 MeV range into one.
- **Check the global sign first.** Tabulate Σ fit / Σ(n·data) per data set. Values of 0.86–0.90 with the median point
  only 1 % off mean that the fit undershoots the peaks: fitting with data-proportional errors pulls a curve below the
  data (Peelle's pertinent puzzle; G. D'Agostini, Nucl. Instrum. Methods A 346 (1994) 306), and unresolved line shapes
  add to it. Report that table with the map; inside it, rank by the size of the deviation.
- From nothing: fit the broad levels to blurred data with a blurred model, using only the data sets the validity
  table passes; add narrow levels afterwards.

## 2. Characterise each region (minutes, no node)

a. **Model low or model high?** A residual-driven level search only adds strength. Where the model is above the
   data, attribute instead: zero each level's amplitudes in the offending pair one at a time, frozen, and tabulate
   the per-dataset change; the cheapest removal identifies the level to bound.
b. **Width.** The χ²-weighted median E_x extent of the runs in differential segments approximates the FWHM of a peak
   or of one lobe of a dispersion shape (`levelscan.cluster_width`). Skip runs in single-energy angular
   distributions (no extent) and runs covering > 50 % of their segment (normalization). Two adjacent regions of
   opposite sign form a dispersion pair: the level sits between them.
c. **Who sees it.** Channels: an (α,n₀) deficit needs α and n₀ width, a neutron-total deficit only neutron width.
   Angles: backward-only means interference with a specific partial wave. Per-angle mean pulls in 10 keV bands
   separate a broad deficit (same-sign pull ≥ 2–3σ at every angle over ≥ 50 keV) from a narrow spike (two bands,
   pull changes sign: a line-shape or resolution question, not a level search). Legendre coefficients of the
   residual bound J: with s = 1/2 in both channels a₆ needs J ≥ 7/2 and a₈ needs J ≥ 9/2; odd terms need
   interference of opposite parities.
d. **Data first.** Exit-channel sums of totals; points within a few keV of a threshold; 0° detectors; a
   normalization that trends with angle; an isolated outlier (one point with χ² in the thousands); unconvolved
   resolution under a narrow level; duplicate data. Each of these has looked like a missing level.
e. **What is already there.** Every level within ±3Γ of the region (J^π, E_x, partial widths, θ², fixed or free);
   every threshold in the window; the Wigner-limit width Γ_W per channel for each candidate J^π (a channel with
   Γ_W < 1 eV cannot carry a level); the ENSDF list for the window (`pyazr.nds.fetch_levels`, see azure2-setup).
f. Write the anatomy of each region into the work log before fitting anything.

## 3. Hypotheses: list them all, test the cheap ones first

| hypothesis | test | cost |
|---|---|---|
| existing level misplaced, too narrow, wrong partition | `relocate`; width floor on the neighbour (`references/width-floor-recipe.md`); free its fixed energy | 1 job |
| interference sign | local sign round over the J^π groups in the window (`round --local-release W`) | 1 job |
| missing level | §4 screen, then §5 level add with control | 1 screen + 1 job per J^π set |
| wrong J^π of an existing level | remove it, screen all J^π at its own E_x and Γ, re-add best against own | ~30 min per level + 2 jobs |
| missing background in a channel | pole, or seeded amplitudes in the dead channel; frozen check first | 1 job |
| data issue | frozen check with the fix (minutes), then refit | 1 job |

## 4. Screen (frozen, or one Gauss-Newton step)

a. **Template screen.** `rmfit levelscan --ex E ... --gamma auto --scales 0.5,1,2 --rank density --ex-range lo,hi
   --refit 0`: one dummy level per J^π; every width partition × sign pattern × energy offset evaluated frozen in the
   full model; the best per J^π stepped once by Gauss-Newton. `--gamma auto` estimates the width from the residual
   runs. Read the global and the local ranking (differential segments within ±0.25 MeV) together: a candidate that
   leads globally but is absent from the local top five is suspect. Templates use only channels with Γ_W ≥ 1 eV and
   are clipped at the θ² cap; uncapped templates in sub-eV channels manufacture levels. The two-channel (pair)
   family is needed for a level that lives in the product of two channels.
b. **Before believing a null** at a resolved feature, put a level with plausible parameters into the screen's own
   variant by hand and evaluate it frozen. Nulls have come from amplitudes ≈ 0 after a failed Brune transform, from a
   width grid ending below the level, and from a missing two-channel template. `levelscan.grid_coverage` answers in
   seconds whether the grid can represent a given level.
c. **Data-only window multi-start** when the screen is null but the structure is clear: only the complaining data
   set, ±0.5 MeV, normalizations fixed, random sign patterns, all J^π. Read off J^π, E_x and Γ. Never transplant its
   parameter vector into the full model (such transplants came back worse than the control by 3×10³ to 8×10⁴ in χ²).
d. Far from the optimum (after a large frozen removal) the global objective cannot separate J^π (a wrong J^π
   compensates in totals); near it, it can. Screen on a well-polished incumbent.

## 5. Test: level add in the full model, against a control

**Without rmfit** (a model of tens of free parameters): `scripts/level_test.py` makes each variant a `Project` edit
of the fitted project, refits it (least squares on the objective's rows and Jacobian, all free parameters) and
compares it with the base refitted the same way:

```bash
python3 $S/level_test.py fit.azr --variant "scan:2-:18.5,19.0,19.5:2/0/2/100000,3/0/2/50000" --out scan.json
python3 $S/level_test.py fit.azr --variant "remove:2- 19.0006" --variant "off:1+ 17.64" \
    --variant "move:1+ 17.64:17.7" --variant "sign:2- 19.0006:3/0/2" --hold-out n7be_total --out tests.json
```

Variants: `add` / `scan` (J^π, energies, channels pair/l/s/width with widths in eV, free), `remove`, `off`
(inactive), `move`, `sign` (flip one width's sign before the refit). Per variant it prints the objective, the number
k of free parameters, Δ = objective(base) − objective(variant), AIC = objective + 2k, BIC = objective + k ln N, and
for nested models the F test, F = (Δ/Δk) / (objective_larger/(N − k_larger)), with its p value; `--hold-out SET`
refits every model without SET and scores SET at that fit (cross-validation on data the fit did not see). On T3
with the 2⁻ level removed (objective 89 939, k = 5), the scan's 19.0 MeV start returns to 19.00063 MeV with
Δ = 89 359 for 3 parameters (p < 10⁻³⁰⁰); removing the level from the full fit raises the held-out n7be_total χ²
from 32.8 to 34 877. A `sign` variant that refits to the base's objective shows that the sign is free to return;
compare signs with the flipped width held (`azure2-fit`, interference). Accept by the rules of `azure2-evaluation`
(F-test p < 0.01 and BIC lower, held-out data not worse, the level physical: `azure2-check`).

**With rmfit** (hundreds of parameters, a cluster):

```bash
rmfit structure add --candidate i [--jpi ...] --local-release 0.5 --local-nfev 80 --pinned-nfev 60 --nfev 150
```

on a campaign initialised from the incumbent (`candidates.json`: E_x, J^π or null, channels with `gamma` in eV with
sign).

- Stages: pinned (the new level's widths only) → local release (every level within ±0.5 MeV, background,
  normalizations) → global polish with the new level's energy free (±max(50 keV, 3Γ), threshold guard). The control
  gets the same local and global budget with the level off. Check that the control's start objective equals the
  incumbent's: a control that starts elsewhere is broken.
- **A broad hypothesis needs a width floor**, otherwise a released broad seed collapses to Γ ~ eV and is scored
  "nothing there" while the neighbours absorb the misfit. Put the floor on the total width
  (`total_width_floor_windows`); the per-channel `width_floor_windows` forces the width into one named channel.
  Full recipe: `references/width-floor-recipe.md`. A narrow hypothesis takes `width_bound_windows` at about twice the
  feature width and `new_level_theta2_cap`.
- **One campaign directory per alternative** (per J^π when several are tested). An accepted add is adopted as the
  campaign's base, so two adds in one directory are sequential, not alternatives.
- `structure add` exits with status 0 on a failed scan: read the printed outcome, not the exit code.
- Run several J^π per node job (three per 16-slot job is typical), all from the same seed.
- Check that an added level's energy is free in the written `.azr`; a fixed energy at a round number means the add
  did not fit it.

## 6. Verdict

- Gain = verified objective of the control minus that of the add. Threshold τ = max(s² k ln N, 0.002·χ²), with
  s² = χ²/N and k the number of parameters the level adds.
- Classes (`structure.classify`): accept; marginal; worse; collapse (θ² < 10⁻⁴ and Γ < 1 eV); narrow (below the
  campaign's width floor); broad (above `max_level_width_ev`); splits (an unresolved same-J^π neighbour gave up
  ≥ 50 % of the new width: one line shape described by two levels); relocated (energy at its window edge); unphysical
  (channels on the θ² cap, failing transform); redistributes (ρ > 0.5: a flag, not a veto). A `nan` in the screen
  table means the evaluation failed, not "no gain".
- Print for the evaluator: E_x, partial widths, θ² per channel, whether each parameter could move; the nearest
  same-J^π level before and after; per-dataset gains and losses; the share of the gain in the top data set (near
  100 % in a 1 keV-stepped set means line shape); what sits at a bound; and a figure of data, incumbent and
  candidate around E_x, per angle (azure2-plot).
- Cross-check an accept against the alternative hypotheses of §3 at the same budget (widen the neighbour, the
  next-best J^π). A level is believable when it beats the alternatives, not only the control.
- Judge a flip or a new level after its neighbours have moved (local release), never frozen alone.

## 7. After an adoption

1. Full free refit (no floor, no window, 400 evaluations): the level must hold on its own.
2. Sign round over its J^π group (local release), then a continuation polish.
3. Re-rank the residuals on the new incumbent (§1): neighbouring regions often change or vanish.
4. Propagate to every model that freezes this block (`references/width-floor-recipe.md`, step 5).
5. Next region. Expect most J^π to be rejected quickly: a 100 keV floor is above the Wigner limit of the high-spin
   groups.

## 8. Pitfalls

- **Level merge.** AZURE2 merges same-J^π levels closer than 1 keV on read. Keep ≥ 1.5 keV from fitted energies.
- **Parameter files by name.** A positional `.sav` is invalid after any level add or remove; read it by name. Never
  carry a fit between models through `<levels>`; move vectors by key (`rmfit.keys.take`).
- **Spawned workers.** Any script that opens an rmfit evaluator needs `if __name__ == "__main__": main()`.
- **Rebuilding narrow on a stripped model fails**: every candidate grows into the missing broad strength.
- **A level pinned at a threshold guard**: test both sides of the threshold at the same budget.
- **Zero amplitudes never move** (zero gradient). A channel the model has never used must be seeded, in units of
  Γ_W (e.g. 0.01 Γ_W), never in keV.
- **Polishes stopped at `max_nfev` are still descending**: compare against a control of the same budget, never
  against the incumbent's number.
- **Grid gaps look like nulls**: see §4b.

## 9. Costs (13C+α model: N ≈ 45 000 points, ~1 200 free parameters, 16 shards)

| step | cost |
|---|---|
| residual dump and ranking | ~10 min, login node |
| template screen, one energy, all J^π | 25–50 min on a node |
| level add with local release and control, 3 J^π | one 16-slot job, 6–8 h |
| polish, 400 evaluations | ~5 h |

## Pointers

- Scripts: `dump_residuals.py`, `rank_map.py`, `peakfind.py`, `eval_peaks.py`, `blur_map.py`, `make_broad_only.py`,
  `blur_compare.py`, `blur_figure.py`, `level_test.py`; `selftest.py`.
- `references/locators.md`: measured recall of the peak finder, blur validity numbers, open improvements.
- `references/width-floor-recipe.md`: the recipe for a structure that is visible but that no search finds.
- Campaign commands (`levelscan`, `structure add`, `round`, `relocate`): rmfit-campaign.
- Level edits by hand, parameter files by name, the 1 keV merge: azure2-fit.
- χ², pulls and residual conventions: `docs/source/user_guide/chi_squared.rst`.
