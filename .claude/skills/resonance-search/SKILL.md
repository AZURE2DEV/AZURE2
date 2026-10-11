---
name: resonance-search
description: Find and test missing or mis-described resonances in an AZURE2 R-matrix fit -- rank the energy regions where the fit misses the data (coherent residual runs), characterise each one (width, which channels and angles see it, model low vs model high, data issue vs physics), list the hypotheses (new level, existing level wrong, interference sign, J^pi, background, data), screen them cheaply, test the survivors as level adds inside the FULL model against a same-budget control, and judge the result on all the data. Use whenever the task is "where does the fit fail", "is there a missing level", "search for resonances", a residual or band map, a level add / levelscan / width-floor test, or the evaluator points at a structure in the data the fit does not reproduce. Reaction-independent; the tools are rmfit's (rmfit-campaign skill) and AZURE2's (azure2-eval); the theory is in r-matrix-analysis.
---

# Searching for resonances

This skill is the PROCEDURE for one question: *what structure is the model missing, and is it a
level?*  The tools it calls are documented elsewhere and are not repeated here:

| need | where |
|---|---|
| running a campaign, polish, sign rounds, `levelscan`, `structure add`, ledger | `rmfit-campaign` skill (§3 loop, §3c width-floor recipe, §9 log) |
| AZURE2 CLI / pyazr, frames, output files, extrapolation traps | `azure2-eval` skill |
| what a level, channel, partial width, interference pattern means | `r-matrix-analysis` skill |
| ENSDF / TUNL level lists, EXFOR | `nds-explorer` skill (`pyazr.nds.fetch_levels`) |

Worked history (13C+alpha, 2026-09-09 .. 10-11): `13C+a/9-11-26_levelscan_gate`, `9-18-26_levelscan_aa1`,
`9-20-26_levelscan_resid2`, `9-21-26_jpi_study`, `10-5-26_ex94_floor` (the model case), `10-11-26_residual_regions`.

## 0. Ground rules (DeBoer's, each learned the hard way)

1. **One model, the whole calculation.**  Search inside the full fit.  To concentrate on a region,
   FIX the parameters outside it (levelFix / channelFix, or a local-release mask) -- never build a
   window or partial-data model and try to reconcile it afterwards ("we don't waste time trying to
   resolve inconsistencies between the low and high energy versions", 2026-10-11).  A data-only window
   fit is allowed only as a SCREEN for J^pi / energy / width (§4c), never as a start vector.
2. **Judge on all the data.**  The verdict is the verified objective of the full model versus a
   same-budget control polish of the incumbent.  A candidate that fixes a different data set than the
   one it was aimed at still counts; the per-dataset split (rho, top-dataset share) is reported, not a
   veto (2026-10-07).
3. **Broad first.**  Fit Gamma >~ 10 keV structure over the whole range before narrow-level alignment;
   a narrow spike or a lineshape mismatch in one finely stepped data set is end-stage work
   (2026-09-14).
4. **Data-definition questions outrank structure.**  A wrong sum of exit channels, an unsummed total
   yield, a bad detector, a single outlier point have each been worth more than any level add.  Rule
   them out (§2d) before searching.
5. **Only finished fits** are starting points (an export with a verified objective), and the evaluator
   adopts -- report, don't adopt.

## 1. Locate: rank the regions where the fit misses

Start from the incumbent's point-by-point residuals, r = (fit - n*data)/(n*sigma) as AZURE2 puts it into
chi^2 (or its sign-flipped twin; be consistent).  Scripts in this skill directory:

```bash
# one single-session evaluation of the full model, parameters from the .sav BY NAME (~10 min)
PYTHONPATH=<R-matrix>:<AZURE2> OMP_NUM_THREADS=1 python3 scripts/dump_residuals.py fit.azr fit.sav work resid.npz
python3 scripts/rank_map.py resid.npz regions.txt fit_output/parameters.out fit.azr   # seconds; re-run freely
```

Check first that sum r^2 of the dump equals the fit's chi^2 (it must, to the digit).  Residual
convention in the dump: r = (fit - data)/sigma, so "model low" = r < 0.

How the ranking works and why:
- Each angle of each segment is one excitation function (a fixed-energy angular distribution is one
  curve in angle).  In each curve, maximal runs of consecutive SAME-SIGN residuals (n >= 3) are found.
- A run is scored by its **coherent chi^2 (sum r)^2 / n** -- the part a smooth change of the curve could
  remove.  Pure scatter (including scatter from underestimated error bars, typical s^2 = chi^2/N of
  10-20) alternates in sign and scores ~1 per run; raw chi^2 would rank the noisiest data set first.
- Each run's score is spread over its Ex extent (min 10 keV) on a 2 keV grid; the map is smoothed
  (6 keV) and cut into regions at its local minima, ranked by integral.  Do NOT merge overlapping run
  intervals instead: with dense multi-angle data the intervals chain into one 6 MeV "region" (first
  try on 13C+alpha, 2026-10-11).  Data-set credit is by overlap fraction.  Runs wider than ~0.6 MeV are
  listed separately as TRENDS (normalization, background, an energy-dependent systematic -- not a level).
- Each region lists the model's levels inside it (J^pi, E, total width from `parameters.out`, FIXED
  flag from the .azr).  The two shapes to look for first: a strong region with **no level** in it
  (candidate new level), and a region whose levels are all **FIXED** (the freeze, not the physics, may
  be the limit).  A "level" of several MeV inside a region is background, not a resonance.
- Also look at the region list sorted by DENSITY (coherent chi^2 per MeV, floor 20 keV): a sharp
  12-point feature in a secondary-gamma set never outranks a 700-point n-total run on totals alone
  (2026-09-24).

The rmfit equivalent inside a campaign is `levelscan`'s candidate step
(`levelscan.candidate_energies`: runs with |r| > 3, clustered in 30 keV bins, `--rank density`,
`--seg-match`, `--ex-range`) -- use it when the next step is a screen anyway; use the scripts here when
the evaluator wants the list first.

Hand the evaluator the list before spending node time: rank, Ex range, centre, coherent chi^2, data
sets, sign (model low or high), width.  The evaluator usually knows several of the regions already --
ask which ones they see and merge the lists.

## 1a. Starting from nothing: data-only peak finder

With no R-matrix calculation yet (a new reaction, or a model too poor to have meaningful residuals), make the
first candidate list from the data alone:

```bash
# the per-point table: dump_residuals.py on ANY .azr that has the pairs and the data segments (its levels do
# not matter -- only key, name, ex, ang, dat, err are used); AZURE2 itself does the lab -> c.m. -> Ex conversion
python3 scripts/peakfind.py table.npz peaks.txt [--ex-range lo,hi] [--levels parameters.out]
```

Per excitation function (each angle of each segment, >= 8 points): a running-median baseline at 0.05 / 0.15 /
0.5 MeV; peaks and dips of z = (y - b)/sigma with prominence >= 4 sigma AND a relative size >= 5 % (small error bars
alone cannot make a candidate); FWHM at half prominence; an adjacent peak + dip = one dispersion candidate.  All
detections go into an Ex density map (each a unit-area Gaussian of its own width, weight min(z, 10)) whose peaks are
the candidates, ranked by the number of DISTINCT data sets that see them.  "(res.)" marks a width at the data's
energy step.  `--levels` scores the list blind against a level set.

Measured blind on 13C+alpha (2026-10-11; 51 data sets, Ex 4.1-13.4 MeV; truth = the merged incumbent's 143 levels;
chance = the same candidates shifted by 50-300 keV):

| level width | within | found (any) | chance | found (>= 2 data sets) | chance |
|---|---|---|---|---|---|
| 1-10 keV | 5 keV | 11/14 | 4.7 | 10/14 | 2.9 |
| 10-50 keV | 10 keV | 12/16 | 8.1 | 9/16 | 6.1 |
| 50-150 keV | 20 keV | 21/25 | 18.9 | 19/23 | 15.2 |
| >= 150 keV | 40-80 keV | all | ~all | all | ~all |

So: it FINDS narrow and intermediate levels (well above chance below ~50 keV), and it does not locate broad ones
in a dense spectrum -- a 300 keV level overlapping three others is not a peak, and the chance rate is already ~100 %.
Precision cannot be judged against a fitted level set (the model is not the truth); 293 candidates for 143 levels,
many on the shoulders and interference lobes of real ones.  Use it accordingly:
- it is a RECALL tool for an overfit-then-eliminate start (ATARI style: put a level at every multi-set candidate,
  fit, remove what the data do not need -- false candidates are cheap, missed ones are not);
- build the broad component separately first: integrated / total cross sections (n-total, (alpha,n) total) and
  ENSDF/TUNL for the broad levels, background poles above the range (ground rule 3, and the 09-24 lesson that
  narrow levels added before the broad ones grow into the missing broad strength);
- once a first fit exists, switch to §1 (residuals) -- it is the better locator from then on.
- Agreement of the two methods is evidence: on 13C+alpha the peak finder's level-free multi-set candidates at
  9.257 and 11.219 MeV coincide with level-free residual regions (9.24-9.32, 11.18-11.23).

Open improvements: a broad-structure estimator (multi-scale curvature or local Breit-Wigner + background fits on
the integrated data), Legendre-order J lower bounds per candidate, ENSDF cross-match (`pyazr.nds.fetch_levels`),
and a test on a sparse-spectrum reaction where precision can be measured.

## 2. Characterise each region (login node, minutes)

a. **Model low or model high?**  A residual-driven level search only ADDS strength.  Where the model is
   ABOVE the data, do a frozen attribution instead: zero each level's amplitudes in the offending pair
   one at a time and tabulate the per-dataset change; the cheapest removal identifies the level to
   bound (2026-10-05, `13C+a/10-5-26_n2_excess/attribute.py`).
b. **Width.**  The chi^2-weighted median extent of the runs in DIFFERENTIAL segments ~ FWHM of a peak or
   one lobe of a dispersion shape (`levelscan.cluster_width`); skip runs in single-energy angular
   distributions (no extent) and runs covering > 50 % of their segment (normalization).  Two adjacent
   regions of opposite sign = a dispersion pair: the level sits between them.
c. **Who sees it.**  Which channels (an (alpha,n0) deficit needs alpha AND n0 width; an n-total deficit
   only neutron width), which angles (backward-only = interference with a specific partial wave), which
   detectors.  Per-angle mean pulls in 10 keV bands separate a broad deficit (same-sign pull >= 2-3 sigma
   at every angle over >= 50 keV) from a narrow spike (two bands, pull changes sign: a lineshape /
   resolution question).  Legendre coefficients of the residual give J^pi hints (with s = 1/2 in both
   channels, a_6 needs J >= 7/2, a_8 needs J >= 9/2; odd terms need opposite-parity interference)
   (`13C+a/9-22-26_an0_legendre`).
d. **Data first.**  Before blaming the model: exit-channel sums of totals, threshold points (a point
   within a few keV of a threshold), 0-degree detectors, a normalization that trends with angle, an
   isolated outlier (one point with chi^2 in the thousands), unconvolved resolution under a narrow
   level, duplicate data.  Each of these has masqueraded as a missing level here.
e. **What is already there.**  List every level within +-3 Gamma of the region with J^pi, energy, partial
   widths, theta^2, and whether it is fixed; every threshold in the window; Gamma_W per channel for the
   candidate J^pi (a channel with Gamma_W < 1 eV cannot carry a level).  The ENSDF/TUNL list for the
   window (`pyazr.nds.fetch_levels`).
f. Write the anatomy of each region into the readme before fitting anything.

## 3. Hypotheses -- list them all, then test the cheap ones first

| hypothesis | test | cost |
|---|---|---|
| an existing level is mis-placed / too narrow / wrong partition | `relocate`; width floor on the NEIGHBOUR (§3c of rmfit-campaign); free its fixed energy | 1 job |
| interference sign | local sign round over the J^pi groups in the window (`round --local-release W`) | 1 job |
| missing level | §4 screen, §5 level add with control | 1 screen + 1 job per J^pi set |
| wrong J^pi of an existing level | remove it, screen all J^pi at its own energy and width, re-add best vs own (the J^pi study) | ~30 min / level + 2 jobs |
| missing background in a channel | pole or seeded amplitudes in the dead channel, frozen check first | 1 job |
| data issue | frozen check with the fix (minutes), then refit | 1 job |

## 4. Screen (cheap, frozen or one GN step)

a. **Template screen** (`rmfit levelscan --ex E ... --gamma auto --scales 0.5,1,2 --rank density
   --ex-range lo,hi --refit 0`): one dummy level per J^pi, every width partition x sign pattern x
   energy offset evaluated frozen in the full model, best per J^pi GN-stepped.  Read the GLOBAL and the
   LOCAL (differential segments within +-0.25 MeV) rankings together: a candidate that leads globally but
   is absent from the local top five is suspect.  Templates use only channels with Gamma_W >= 1 eV and are
   clipped at the theta^2 cap (2026-09-20: uncapped templates in sub-eV channels manufactured a +1,693
   "level").  The pair family (two-channel shapes) is needed for levels living in the product of two
   channels (2026-09-21).
b. **Before believing a NULL** at a resolved feature: drop a level with plausible parameters into the
   screen's own variant by hand and evaluate it frozen.  Three nulls on 13C+alpha were bugs or grid gaps
   (amplitudes ~0 from a failed Brune transform, width grid topping out below the level, no two-channel
   template).  `levelscan.grid_coverage` answers "could the grid represent this level" in seconds.
c. **Data-only window multi-start** (only the complaining data set, +-0.5 MeV, norms fixed, random sign
   patterns, all J^pi) when the screen is null but the structure is clear: read off J^pi, energy, width.
   Never transplant its vector (2026-10-06: +3k..+83k vs the control).
d. Far from the optimum (a large frozen removal) the global objective cannot tell J^pi apart
   (wrong-J^pi compensation in totals); near it, it can.  Prefer screens on a well-polished incumbent.

## 5. Test: level add in the full model, against a control

`rmfit structure add --candidate i [--jpi ...] --local-release 0.5 --local-nfev 80 --pinned-nfev 60 --nfev 150`
on a campaign initialised from the incumbent (candidates.json: energy, jpi or null, channels with
`gamma` in eV with sign):
- pinned stage (new level's widths only) -> local release (every level within +-0.5 MeV, background,
  norms) -> global polish with the new level's energy free (+-max(50 keV, 3 Gamma), threshold guard);
  the control gets the same local + global budget with the level OFF.  Check the control's START equals
  the incumbent -- a control that starts elsewhere is broken.
- **A broad hypothesis needs a width floor** (`width_floor_windows`, +-0.5 keV on the seed energy, on
  the channel the anatomy says carries it).  Without one, a released broad seed collapses to eV and is
  scored "nothing there" while the neighbours re-absorb the misfit (13C+alpha Ex 9.4: the accepted
  198 keV 1/2- had collapsed in the first test).  Narrow hypothesis: `width_bound_windows` at ~2x the
  feature width, `new_level_theta2_cap`.
- **One campaign directory per alternative** (per J^pi when several are tested): an accepted add is
  ADOPTED as the campaign's base, so two adds in one directory are sequential, not alternatives.
- `structure add` exits 0 on a failed scan: read the printed outcome, never the exit code.
- Several J^pi per node job (three per 16-slot job on 13C+alpha), all from the same seed.

## 6. Verdict

- Gain = verified objective of the add minus that of its control.  tau = max(s^2 k ln N, 0.2 % obj)
  with s^2 = obj/N, k = parameters the level adds.
- Outcomes (`structure.classify`): accept / marginal / worse / collapse (vanished: theta^2 < 1e-4 AND
  Gamma < 1 eV) / narrow (below the campaign's width floor) / broad (above `max_level_width_ev`) /
  splits (an unresolved same-J^pi neighbour gave up >= 50 % of the new width -- one lineshape described
  by two levels) / relocated (energy at its window edge) / unphysical (channels on the theta^2 cap,
  failing transform) / redistributes (rho > 0.5: a flag, not a veto).
- Always print, for the evaluator: the level's energy, partial widths, theta^2 per channel, whether each
  parameter could have moved; the nearest same-J^pi level before/after; the per-dataset gains and
  losses; the share of the gain in the top data set (near 100 % in a 1-keV-stepped set = lineshape);
  what sits at a bound.  And a figure: data, incumbent, candidate around the energy, per angle.
- Cross-check an accept with the alternative hypotheses of §3 at the same budget (widen the neighbour,
  the next-best J^pi): a level is believable when it beats the alternatives, not only the control.

## 7. After an adoption

1. Full free refit (no floor, no window, 400 evaluations): the level must hold on its own.
2. Sign round over its J^pi group (local release); then a continuation polish.
3. Re-rank the residuals on the new incumbent (§1): neighbouring regions often change or vanish.
4. Next region.  Expect most J^pi to be rejected quickly; a 100 keV floor is above the high-spin groups'
   Wigner limits.

## 8. Traps (each one cost at least a day)

- AZURE2 merges same-J^pi levels closer than 1 keV on read; keep >= 1.5 keV from FITTED energies.
- A positional .sav is invalid after any level add/remove; read by name.
- Never carry a fit between models through `<levels>`; move vectors by ParamKey (`rmfit.keys.take`).
- Any script that opens an rmfit evaluator needs `if __name__ == "__main__": main()` (spawn workers).
- Added levels must have FREE energies (pre-09-19 adds were silently fixed at round numbers).
- Stripping broad levels to "rebuild narrow" fails: every candidate grows into the missing broad strength.
- A level pinned at a threshold guard: test both sides at the same budget (13C+alpha 3/2- at the n2
  threshold, -969 below).
- Zero amplitudes never move (zero gradient): a channel the model has never used must be SEEDED, in
  units of Gamma_W, never in keV.
- Polishes that stop at max_nfev are still descending: compare against a control of the same budget,
  never against the incumbent's number.

## 9. Costs (13C+alpha merged model, N ~45k, ~1,200 free, 16 shards)

| step | cost |
|---|---|
| residual dump + ranking | ~10 min login node |
| template screen, one energy, all J^pi | 25-50 min on a node |
| level add with local release + control, 3 J^pi | one 16-slot job, 6-8 h |
| polish 400 | ~5 h |

## 10. Refinement log (append dated entries)

- 2026-10-11 — skill created (DeBoer: "This may be a good time to work on a skill to search for
  resonances").  Consolidates rmfit-campaign §3 step 4, §3c and the 09-12 .. 10-07 log entries; adds the
  coherent-run region ranking (`scripts/`), first used on the merged 13C+alpha incumbent
  (`13C+a/10-11-26_residual_regions`).
- 2026-10-11 (evening) — §1a added (DeBoer: "if you are starting from nothing ... do you have some peak finder type
  routine"): `scripts/peakfind.py` + `scripts/eval_peaks.py`, measured blind on 13C+alpha against a chance baseline.
  Lesson, twice in one day: interval clustering of dense multi-angle data chains everything into one region -- cluster
  by a density map.  And a recovery number means nothing without the chance rate at the same tolerance.
