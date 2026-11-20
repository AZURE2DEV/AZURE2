# Data audit of the active data files of `$azr`

$date.  Requested by rmfit for the campaign `$campaign_dir` ($reaction).  Rules and report format at the end.

Data directory `$data_dir` (files are lab-frame energies of the projectile, see the segment table; columns are
energy, angle, value, error for differential data, energy, value, error for angle-integrated; a `.dat` may carry a
header of `#` lines).  Which rows a segment uses is set by its energy and angle window.  Read the segment table below
first.  Tools: `python3 -I` scripts in `$out_dir/`; EXFOR through `pyazr.nds` (`search_exfor`, `fetch_exfor`; how-to
`~/AZURE2/.claude/skills/nds-explorer/SKILL.md`, fallback dataexplorer REST as in
an archived campaign) with `PYTHONPATH=$pythonpath`.  Worked example of this audit:
an archived campaign and its `audit_*.py` scripts (reuse them).

## What to check (every file of an active segment)

1. **Integrity**: non-numeric rows, zero or negative values or errors, error > value, duplicated rows, the last row
   repeated (export artefact), column count.
2. **Units and scale**: ratios between overlapping sets of the same reaction, band by band (median and 16-84 %); an
   absolute anchor compared with EVERY independent measurement of the same quantity, not only the one in the fit;
   mb vs b, per sr vs integrated, S-factor vs cross section.
3. **Frames**: lab energy of the right projectile (test against the thresholds in the table); lab vs c.m. angles and
   cross sections against the segment type; the first row of a set against the reaction threshold.
4. **Uncertainties**: statistical-only where the source gives a systematic one; an advertised error floor that does
   not hold; implausibly small, constant, or integer-percent errors (an EXFOR FLAG column read as the error); a
   single row whose weight exceeds hundreds of typical points.
5. **Double counting**: the same measurement under two names, or in two overlapping windows of two segments.
6. **Energy scale**: for a set with narrow structure, note whether a shift scan has been done (readme) and whether
   the file's resonance energies agree with the level scheme (ENSDF via `pyazr.nds.fetch_levels`).
7. **Against the source**: re-fetch the EXFOR entry where one exists and compare values, units and error columns
   row by row; list files with no identifiable source.
8. **Experimental effects, from the reference itself** (the evaluator 2026-11-19: "we need to do a careful job for each data
   set that we enter to check the reference for experimental resolution, energy scale uncertainties and observable
   systematic uncertainties and then use AZURE2's experimental effects capabilities to model them").  Read the
   publication (local library and Box first: the `rclone-box-drive` skill; EXFOR alone is not enough) and record, per
   set, with the page, table or equation it comes from:
   - **detector geometry** (size, distance, angles) -> finite-geometry attenuation, Q_k coefficients;
   - **target thickness and beam energy spread** -> target integration or beam-energy convolution (`<targetInt>`);
     whether the quoted energies are beam or effective (mid-target) energies (the `azure2-target-energies` skill);
   - **beam-energy calibration uncertainty** (and separate accelerators or runs) -> an energy shift with a penalty;
   - **source of the absolute scale and its uncertainty** (an adopted reference cross section, a normalization to
     another measurement) -> `dataNormError` in percent, and which sets share one scale;
   - **whether the published values are already corrected** for each of these.
   Write "not given in the reference" where it is not, never a silent default.  Report it as one table, a row per
   set: effect | value from the reference (citation) | how the `.azr` treats it now | what modelling it would take.
   Whether to add the effect to the fit stays the evaluator's decision.  In one archived capture fit all three
   (Q-coefficients, a 9 keV energy error in the database copy, the absolute-scale anchor) were found only after
   fitting.

## The model's segment table (from pyazr; key = `<segmentsData>` line counted from 1)

$segment_table

## Active data files

$data_files

## Already known (do not re-report)

$known

## The evaluator's framing (campaign readme, first lines)

$readme_head

$common
