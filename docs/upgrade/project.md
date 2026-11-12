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

(State when this page was written.  Since pyazr 3.0, 2026-11-11, `azrfile.py`, `datasets.py` and
`roundtrip_test.py` are deleted: pyazr reads format 3 only, through `pyazr.Project` on the engine's document and
grammars; `tests/pyazr/project_test.py` checks the round trip; `docs/upgrade/pyazr_on_records.md`.)

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
   `roundtrip_test.py` and the edit tests (since pyazr 3.0: `pyazr.Project`, checked by `project_test.py`,
   `azr_layout_test.py` and `project_nuclear_part_test.py`).
5. One mass table (the gate's second half): `fill_masses` and the engine's masses read the same table.

## 3. Decisions

**PI, 2026-11-09:** P2 and P3 accepted as proposed. P4 replaced: "remove sections that are not needed". New: "a
converter from the old to the updated file format so people can use the old azr files" (P5, §4). P1 (the archive): it is on another machine; "in here you have to test on what you have or can make". The
gate is therefore the 189 project files under `tests/` and `examples/`, and variants made from them: every file
with CRLF, CR and no final newline; with blank lines, trailing spaces and comments added between and around the
sections; and with each section removed in turn and each single field edited in turn.

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

## 4. The updated format (format 2) and the converter

**What changes** (nothing else; every other byte of a file means what it meant):
- **`format=2`** in `<config>`'s `key=value` lines marks a file of the updated format. A file without it is format 1.
  Older builds skip unknown keys, so a format-2 file still opens in them, apart from the changes below.
- **Removed sections**, which no part of AZURE2 reads:
  - `<lastRun>` (the Qt GUI's runtime options; in 92 of the 189 test files);
  - the run settings of `<mcmc>` (`nwalkers`, `nsteps`, `chainspread`, `energyspread`, `nthreads`,
    `usereducedwidths`), which the CLI asks for. `<mcmc>` keeps its `<parameters>` priors, which are read, and goes
    when it has none.
- **`useAdaptiveGrid` moves from `<potential>` to `<config>`.** It selects the target-integration and rate grids, not
  the potential (`Config.cpp` says so), and in 44 of the 47 files with a `<potential>` block it is that block's only
  reason to exist. `<potential>` remains only where a hybrid potential is used (`useHybridPotential=1`). (Format 3
  has no such key: the converter drops it, `project_format3.md` §3.2.)
- **The reading rules of P2 and P3.** One reading of every file: whole-line tags with any line ending; case-sensitive
  keywords; quoted paths; the engine's level grouping; and messages naming the line where a row or keyword used to be
  skipped silently.

**The converter** (P5): `io::ConvertToFormat2(Project&)`, a function of the library, and three entry points to it:
`AZURE2 --convert old.azr new.azr`, `pyazr.convert_azr(old, new)`, and the engine itself. A format-1 file still
loads: the engine converts it in memory, notes that it did, and reports what changed (which sections it removed or
moved, and a list of the lines it reads differently under P2/P3, if any). A conversion is deterministic, and its
output round-trips byte for byte.

**Checks.**
- Every format-1 file under `tests/` and `examples/` (and the archive, P1): conversion succeeds, and the
  outputs of the original and the converted file are bit-identical (179 × 3).
- Converting a format-2 file changes nothing. The converter's report lists exactly the removed and moved sections.
- A file whose meaning changes under P2/P3 is reported line by line, with both readings.

## 5. pyazr on the engine's document (step 3, as built so far)

- `_azure2.Document` binds `io::Document`. **pyazr now needs the compiled engine module to read a project:** a
  pure-Python fallback would be the second parser that Phase 3 removes.
- `AzrModel.from_file` takes `<levels>` from the document (whole-line tags), and every section search in `azrfile`,
  `datasets`, `azure2`, `bands` and `nuclei` uses the document's tag rule (`_tag_index`: a whole trimmed line). Band
  directories keep their names (the same bytes are hashed).
- Segment, extrapolation and `<targetInt>` lines are read by the engine's own parsers (`SegLine`, `ExtrapLine`,
  `TargetEffect`), bound as `_azure2.parse_segment_line`, `parse_extrapolation_line` and `parse_target_effect_line`.
  Before each switch the Python reader and the engine's were compared on every project file:
  - 644 segment lines and 232 extrapolation lines agree in every field;
  - of 567 screening, Doppler and density views, 561 agree; the 6 others are the two deliberately malformed test files,
    which the engine refuses and the regular expressions read silently (P3).
- Levels: the engine's rule, `io::LevelGrouper`, shared by `CNuc::Fill` and `AzrModel`.
- `<config>` keys only inside `<config>`, the last line winning; the output directory by position.
- Refused as the engine refuses them (P3): malformed `<abundances>` lines and prior centres that are not finite, or a
  norm centre ≤ 0.
- `<thm>` through the engine (after step 4): `_azure2.read_thm_block` (api/src/pybind_thm.cpp, the engine's
  `ParseThmBlockLines`) returns the options and the experiments with each key's value as written; `thm_option`,
  `parse_thm_cbackground`, `thm_nuclides`, `thm_global_opticals` and `read_thm_weight_table` replace the Python copies
  of the grammar, the nuclide mass table, the optical-potential table and the weight-file reader (about 650 lines).
  pyazr keeps the writers (the GUI's compose) and, before the data are read, three checks the engine makes only once
  they are: a global optical potential's validity range, the vertexModel=dw global rules, and the segments of an
  experiment (THM, free norm, exit pairs of cbackground).
- Remaining second reading: `<scatteringLengths>`, whose tokens pyazr checks; the pairs and spins are checked when the
  engine loads the model.
