---
name: azure2-plot
description: Make figures of AZURE2 results and show them to a person -- excitation functions with pulls, angular distributions, S factors with uncertainty bands and level/interference decomposition, rate ratios -- from a freshly recalculated fit ("show me the fit", "plot the fit", "make the S-factor figure"). Use prc-figures for figure style; azure2-calculate for the numbers; azure2-uncertainty for the bands themselves; azure2-check to judge what the figure shows.
---

# Figures and visual review of AZURE2 results

A figure of an R-matrix fit is a comparison of AZURE2's c.m. observable with data that AZURE2 has normalised and
shifted. It is correct only if the curve is the curve AZURE2 computed (the fit column of `AZUREOut_*.out`, or the same
pyazr call), the data are drawn as AZURE2 compares them (multiplied by the fitted norm n, at the fitted energy shift),
and the grid resolves every structure. This skill covers publication figures from pyazr and the output files, and
showing a result to a person. This version of AZURE2 has no graphical front end
(`docs/source/user_guide/graphical_interface.rst`): a result is shown as figures made from fresh output.

## When to use, and when a sibling instead

- Numbers behind the figure (extrapolation, S factor, decomposition, rates): azure2-calculate.
- Bands (covariance, samples): azure2-uncertainty computes them; this skill draws them.
- Whether the fit is acceptable (pulls, χ²/ν, residual structure): azure2-check.
- Style (closed frames, inward ticks, STIX fonts, colour roles, overlap checks): the user-level prc-figures skill.

## Rules for every figure

- **Check every curve against the fit column of `AZUREOut_*.out`.** Interpolate the plotted curve to the data
  energies and compare with column 4 of the output file for the same segment. This catches three failures:
  - observables an extrapolation cannot compute: a `<segmentsTest>` segment is never UPOS, so an extrapolated curve
    at a γ detector's angle is the primary particle's angular distribution, not the secondary-γ one the fit uses (up
    to a factor 2 off on one (α,n₂γ) case). Compute UPOS curves in data mode on a pseudo-data grid;
  - grids too coarse for narrow resonances: a 5 keV grid misdrew peaks by up to a factor 11; keV-wide levels need
    0.25–0.5 keV steps (or `pyazr.tabulate`, azure2-calculate);
  - single-point spikes at channel thresholds: a grid point on a threshold can return nonsense (10²⁸ b). Drop a point
    that departs by more than ~5 % from the mean of its neighbours within a few keV of a threshold, testing the total;
    a partial cross section's own opening is steep but genuine. The same failure in data points is an azure2-check
    matter.
- **Draw the data as AZURE2 compares them**: multiplied by the fitted norm n and at the fitted energy shift. A fixed
  non-unity norm lives only in the `.azr`, a fitted one only in the `.sav`; pyazr's `m.segment_norms(x)` returns the
  norm each segment carries at `x`.
- **State the frame and unit on every axis**: E_c.m. or E_lab (MeV or keV), E_x for plots across entrance channels
  (the only common axis), θ_c.m. or θ_lab (deg), σ (b), dσ/dΩ_c.m. (b/sr), S (keV b or MeV b), N_A⟨σv⟩
  (cm³ mol⁻¹ s⁻¹).
- **Take S from AZURE2** (`calculate_sfactor_rwa`, `m.sfactor`, the S columns); never re-derive the Sommerfeld factor.
- **Bare or folded, say which.** A test-segment curve is folded by any active `<targetInt>` line listing its key; a
  bare curve needs those lines inactive (azure2-calculate, Extrapolation).
- **Evaluate at the fit's parameters** (`m.load_sav("output/param.sav")` or `pyazr.best_fit_params`), not at the
  `.azr`'s `<levels>`, unless the figure is about the `.azr` itself.

## Recipes

Code for each figure type is in `references/recipes.md`:

| Figure | pyazr calls |
|---|---|
| excitation function with pulls | `calculate_rwa`, `calculate_energies`, `m.cross`, `m.cross_err`, `segment_norms` |
| angular distribution | `calculate_angles` per segment; `pyazr.angular_distribution` for Legendre coefficients |
| S factor with band and decomposition | `extrapolation_bands(..., quantity="sfactor")`, `m.sfactor`, `without_level`/`only_level` |
| rate ratio | `reactionrates.out` or `pyazr/examples/reaction_rate.py`, divided by a reference rate on the same T9 grid |

## Showing a result to a person

Figures made from `output/` show the last run, not the `.azr`. After any edit, recalculate before drawing, every time.

1. Finish every `.azr` edit.
2. Run a mode-1 calculate with the project's flags and the right parameter file.
3. Confirm fresh output: `grep Total-Chi output/chiSquared.out` and the `AZUREOut_*.out` timestamps.
4. Draw the figures (Recipes) and check each curve against the fit column (Rules).

```bash
cd <fit directory>
# 2. calculate: menu 1, external parameter file, blank EC file, exit
$AZURE2_BINARY <project flags> --no-readline <name>.azr <<EOF
1
output/param.sav

7
EOF
# 3. confirm
grep Total-Chi output/chiSquared.out
ls -la output/AZUREOut_*.out | head            # timestamps from just now
```

### Things that make the picture wrong without a warning

- **The project's own CLI flags.** The CLI does not read Runtime Options saved with a project. Take `--use-brune`,
  `--ignore-externals`, `--no-long-wavelength` and the rest from the project's job script; a calculate under other
  flags draws different physics.
- **The right parameter file.** A blank parameter file uses `<levels>` with the nominal `dataNorm`/`energyShift`
  fields of `<segmentsData>` (not n = 1); `output/param.sav` uses the fitted values, including norms and shifts. The
  two can differ by a large factor in χ² (127,588 against 96,811 on one model). Blank is right only to verify what
  `<levels>` contains. Never write fitted norms into the nominal fields (azure2-fit).
- **Back up `output/` first** when its contents are a cited record (a published baseline, numbers in a manuscript):
  the calculate overwrites `chiSquared.out`, `normalizations.out`, `parameters.out` and every `AZUREOut_*.out`.
  Copy them to a sibling directory (~1 MB). A data-mode calculate leaves `.extrap` and `.acoeff` files alone.

### Budget the time

A mode-1 calculate on a large multichannel model with target integrations is slow. On one model (107 segments,
~4,800 points, `<targetInt>` entries up to 1,000 integration points): ~33 min in a fresh directory with empty
`output/`, ~5 min in a live directory. Tell the person the calculate is running, and send it to the batch system
rather than a login node: such runs use several GB and effectively one core regardless of `OMP_NUM_THREADS`.

### Before handing the figures over

Say what the person is looking at:

- which `.azr`, and what changed in it since the last review;
- which binary (date) and which parameter file the curves come from;
- the total χ² of the calculate just run;
- which segments are worth checking for this change.

If `<segmentsData>` changed, name the segment keys touched and the data file each holds. A key is the 1-based position
among all `<segmentsData>` lines, inactive ones included: the key `<targetInt>` binds to and `chiSquared.out` and
`normalizations.out` print. Quote those files, not a count by eye.

## Pointers

- No graphical front end: `docs/source/user_guide/graphical_interface.rst`. Output columns:
  `docs/source/reference/output_files.rst`.
- `references/recipes.md`: figure code. Style: prc-figures (user-level skill).
