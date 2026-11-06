# normalizations_out

`normalizations.out` holds one complete row per segment whose normalization is fitted.

`EData::WriteOutputFiles` wrote the first four fields of a row (segment key, data file, angle
bounds) only for segments with a varied normalization, but the remaining six (frames, energy
bounds, normalization, shift) for every segment: a project with fixed and varied normalizations
got a six-field row without key, file or angles for each fixed one (audit 2026-10-25, item 18;
seen in 11B+a/10-25-26_wang_an_global). The docs and the GUI help also claimed the file was read
back when `param.sav` is loaded; nothing reads it, the fitted normalizations travel in `param.sav`.

Now a segment is skipped unless its normalization is varied (a norm-group member keeps its own
flag), every row has ten fields, and the docs describe the file as the report it is.

## Checks (`check.sh path/to/AZURE2`)

- (a) `tests/13N` with the normalization of segment 3 varied: the file is the header plus exactly
  one row of ten fields, starting with the segment's key, file and angle bounds;
- (b) `tests/norm_groups` with the two `normgroup` members' normalizations varied and the other
  segments fixed: the rows are exactly those two, ten fields each.

dev f18b36e: (a) 16 rows, 15 of six fields; (b) passes there too, since both of that model's segments are varied.
