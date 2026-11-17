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

## The model's segment table (from pyazr; key = `<segmentsData>` line counted from 1)

$segment_table

## Active data files

$data_files

## Already known (do not re-report)

$known

## The evaluator's framing (campaign readme, first lines)

$readme_head

$common
