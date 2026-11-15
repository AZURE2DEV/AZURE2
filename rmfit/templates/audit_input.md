# Input-file audit of `$azr`

$date.  Requested by rmfit for the campaign `$campaign_dir` ($reaction).  Rules and report format at the end.

Read the model file `$azr` (AZURE2 project file; anatomy and token meanings: `~/AZURE2/.claude/skills/azure2-eval/SKILL.md`,
section ".azr file anatomy"; the parsers are `~/AZURE2/include/NucLine.h` and `include/SegLine.h`).  pyazr reads it
in Python: `from pyazr import AzrModel; m = AzrModel.from_file(path)` (levels, channels, `mass_check()`), and
`from pyazr.datasets import SegmentSet; SegmentSet.from_file(path)` (segments; `.table()`), with
`PYTHONPATH=$pythonpath`.  Use `python3 -I` for scripts that read data files.

## What to check

1. **Particle pairs**: mass (AME2020, atomic), charge, spin, parity, excitation energy, separation energy, channel
   radius identical on every level line of a pair and right against the nuclear data (`pyazr.masses`); thresholds
   must reproduce to ~1 keV.
2. **Channel sets**: every (pair, L, S) allowed by angular momentum and parity present, none missing at low L, the
   same set on every level of a J^pi group.
3. **Levels**: same-J^pi pairs within 1.5 keV; free amplitudes that are exactly 0 (never move; already counted by the
   machine audit, so report WHICH open low-L channels are dead, by pair, not the count); levels within 3 keV of a
   threshold; fixed/free energies and widths against the campaign's rules (spec below); round-number fixed energies
   (the old levelscan signature); a level sitting at a bound of the campaign's policy windows.
4. **Segments**: file exists; entrance/exit keys match what the file holds; the energy window against the file's
   range and the exit threshold (rows silently outside the window); isDiff against the angle column and its frame;
   "sum" tails complete for every channel open in the window; composite SUM / RATIO / UPOS tails listed.
5. **Normalizations and shifts**: every `varyNorm 1` with `normError 0` (free, unpenalized); fitted norms outside
   0.7-1.4 in the companion `.sav`; every norm/shift setting against the evaluator's recorded decisions; norm groups
   and covariance files consistent (a set with a covariance file and a free norm would double count its systematic).
6. **`<targetInt>` lines**: keys bind to the intended segments; in sharded runs only angle-integrated segments may be
   convolved; a line naming a `<segmentsTest>` key folds the extrapolation (Golden rule).
7. **Companion files**: `.sav`/`.par` next to the model with legacy `j=..._la=...` names (ignored by AZURE2);
   `output/intEC.*` older than the model.

## The model's segment table (from pyazr; key = `<segmentsData>` line counted from 1)

$segment_table

## Already known (do not re-report)

$known

## The evaluator's framing (campaign readme, first lines)

$readme_head

$common
