# segline_formats

Segment-line layouts the reader must accept (audit 2026-10-25, 06_parameters_outputs_api C7 and C8,
checklist item 29).

A `<segmentsData>` line's three shift fields (`shift shiftError vary`) are optional: older files
carry one shift value or none before the data file. `SegLine.h` decided the layout by whether a
`>> double` read succeeded, so a file whose name starts with digits in the no-shift layout
(`... normErr 12C_p.dat`) was read as a 12 MeV shift with the file `C_p.dat`, and the one-shift
layout (`... normErr 0.001 data/foo.dat`) left the stream failed and the file empty; in both cases
`EData` printed "Could Not Fill Segment" and the fit went on without the data. The layout is now
decided by whole tokens (`WholeNumber`: a token is a number only if `strtod` consumes all of it).

A `<segmentsTest>` line's keyword reader (`ExtrapLine.h`) knew `sample`, `total`, `covariance`
(and now `spin`) and stopped at anything else, so a hand-written `sqrtshift 0 0 0 sample 1` lost
its `sample`. The data-only blocks (`sqrtshift`, `shiftgroup`, `shiftfngroup`, `normgroup`,
`absolute`, `shiftfn`) are now skipped with their values, since a test segment has no shift, norm
or shift function.

## Checks (`check.sh path/to/AZURE2`)

- (a) `tests/13N` with segment 3 in the no-shift layout and its file `84deg_meyer.dat` (no directory, the name itself starts with digits):
  the total equals the expected one;
- (b) segment 3 in the one-shift layout `0.001 data/meyer_84.dat` equals the three-field layout and
  differs from the unshifted total (the 1 keV shift is applied);
- (c) `tests/natural_neon` in mode 3 with an n + 20Ne test segment written
  `sqrtshift 0 0 0 shiftgroup G sample 1` equals the one written `sample 1` (weighted by the 20Ne
  abundance, which itself moves the bare line by more than 1 %).

dev e7d72bd: (a) and (b) drop the segment, (c) loses the sample weighting.
