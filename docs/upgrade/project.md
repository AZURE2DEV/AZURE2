# Phase 3: one project model — design

ROADMAP Phase 3: `io/Project`, one `.azr` reader/writer through which the engine and pyazr read; old files keep loading;
byte-identical round trips of the archive of 48,922 `.azr` files. Gate: byte-identical round trips; one mass table.

## 1. What exists (inventory of 2026-11-09, `feature/upgrade` 02bde443)

There are two parsers: the engine's (read only) and pyazr's (read and write). The engine has no `.azr` writer since the
Qt GUI's `AZURESetup::writeProject` was removed (decision 3). `qt_parity.md` §1 is the fullest description of the
format.

The engine re-opens the file once per reader, each with its own scan:
- `Config::ReadConfigFile`, with `ReadPotentialBlock` and `ReadThmBlock`;
- `CNuc::Fill` (`NucLine`) and `ReadAbundances`;
- `EData::Fill` (`SegLine`, `ComponentList`), `MakePoints` (`ExtrapLine`), `ReadTargetEffectsFile` (`TargetEffect`)
  and `ReadPriorCentres`;
- `ScatteringLength::ReadConstraints` and `ParameterLimitsManager::ReadParameterSettings`;
- two CLI pre-scans in `AZURE2.cpp`.

Only `CNuc` and `EData` use the `InputFiles` cache. Some readers match tags exactly, some trim.

pyazr reads every section in Python:
- `AzrModel` in `azrfile.py`, about 4000 lines; `SegmentSet` and `TestSegmentSet` in `datasets.py`;
- `<thm>` twice: a Python port, plus the engine's own parser through the binding (`check_thm_block`).

pyazr also writes every section, by editing the text in place: `AzrModel.to_text` and about 40 setters, `save_fit`,
`modelavg` and `merge_projects`. `tests/pyazr/roundtrip_test.py` checks byte-identical round trips of the 189
project files in `tests/` and `examples/`.

**Where the two readers disagree** (each is a file whose meaning depends on the reader):
1. Tags: pyazr finds `<levels>` as a substring, even inside a comment; it needs `<segmentsData>` with LF line endings;
   it looks for `<scatteringLengths>`, `<abundances>` and `<targetInt>` only after `</levels>`. The engine searches the
   whole file, and some of its readers trim.
2. Level grouping: pyazr groups consecutive lines with the same ID. The engine merges lines with the same ID within
   1 keV, or with the same energy, anywhere in the J^π group.
3. Segment lines: pyazr breaks on quoted data paths with spaces. pyazr reads the one-shift layout as shift 0 in
   `datasets.py`, but as the shift in `azrfile.py`. Keywords are case-insensitive in pyazr and case-sensitive in the
   engine. `isDiff % 10 == 2` against `== 2` for the phase-shift columns.
4. `<config>`: pyazr accepts `relativistic=` and the like anywhere before the block ends, and finds the output
   directory by its comment text. The engine reads by position and only between the tags.
5. Leniency: pyazr skips short `<abundances>` lines (the engine refuses them), accepts nan/inf prior centres (the
   engine refuses them) and skips malformed prior rows silently. The engine silently ignores `<parameterSettings>`
   rows that are not 9, 12 or 14 tokens long, and an unknown segment keyword ends its segment line without a message.
6. Read by one side only: `<lastRun>` and `<mcmc>` (pyazr, kept verbatim); `<potential>`, most `<targetInt>` columns
   and the 12/14-token settings rows (engine only).

## 2. The design

**`io/Project`** (a C++ library, `azure2_io`, depending on nothing else in AZURE2; bound to Python):
- **The document.** The file as a list of sections, each section as its lines with their original bytes. Text outside
  sections, comments and blank lines are kept. `Write()` of an unedited document reproduces the input byte for byte.
- **Typed views.** Per section a parser into plain values: config, potential, thm, levels (`NucLine`), abundances,
  segments (`SegLine`), extrapolations (`ExtrapLine`), targetInt, scatteringLengths, parameterSettings and the mcmc
  priors. Each returns its values, or a list of diagnostics (file, line, column, message, severity). These are the
  engine's present line parsers, moved and made pure (no `Config`, no output stream).
- **Edits.** A typed edit replaces only the bytes of the field it changes (pyazr's rule today); everything else is
  untouched.

**Order of the work.**
1. The document and the round trip: `io/Project` reads and writes the file losslessly. Check: every `.azr` under
   `tests/` and `examples/`, and the archive (§3, P1), byte for byte.
2. The engine reads through it: each of the engine's ten readers takes its section from the `Project`, one reader per
   commit, and the file is opened once. Check: 179 × 3 bit-identical outputs, the parse-regression tests, the same
   diagnostics.
3. pyazr reads through it: `AzrModel`'s parsing calls the binding (as `<thm>` already does), so pyazr and the engine
   read every file the same way. The places where they disagree today (§1) need decision P2.
4. pyazr writes through it: the setters become typed edits of the `Project`; `to_text` is `Write()`. Check:
   `roundtrip_test.py` and the edit tests.
5. One mass table (the gate's second half): `fill_masses` and the engine's masses read the same table.

## 3. Decisions for the PI

- **P1, the archive.** The ROADMAP's 48,922 `.azr` files are mentioned nowhere else in the repository or its history.
  Where are they, and may a copy be used on this node, read only?
- **P2, one reading of every file.** Where the two readers disagree (§1), the engine's reading is what produced every
  published result, so the proposal is that it becomes the reading, with these exceptions, each a diagnostic rather
  than a silent choice:
  - tags: matched as whole lines after trimming, with any line ending (CRLF accepted everywhere);
  - keywords: case-sensitive, as the engine (pyazr's case-insensitive matching goes);
  - quoted paths with spaces: allowed everywhere (the engine's `ReadPathToken`);
  - level grouping: the engine's rule (same ID within 1 keV, or same energy), so that a fit written by pyazr and read
    by the engine has the same levels;
  - the silent cases become messages: a `<parameterSettings>` row of another length (WARNING, with the line), an
    unknown segment keyword (ERROR, with the keyword), a malformed `sample` or `shiftgroup` (ERROR).
  A file that both readers read the same way today reads the same way after the change.
- **P3, strictness.** pyazr is lenient where the engine refuses (short abundance lines, non-finite prior centres).
  Proposal: the engine's strictness everywhere, since a file pyazr accepts and the engine refuses cannot be fitted.
- **P4, unknown sections** (`<lastRun>`, `<mcmc>`'s run settings, anything else): kept verbatim and never interpreted,
  as pyazr does now.
