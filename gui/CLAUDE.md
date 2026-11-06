# gui/ — the Qt setup utility (AZURESetup)

Read the root `CLAUDE.md` first. Most AZURE2 users only ever see the GUI, so a
feature that pyazr can use but the GUI silently drops or corrupts is a bug.

## How it fits together

- `AZURESetup.cpp` reads and writes the whole `.azr` (`readFile` / `writeFile`).
  Each tab owns one block through a Qt model: `PairsModel`, `LevelsModel`,
  `SegmentsDataModel`/`SegmentsTestModel`, `TargetIntModel`,
  `ScatteringLengthsModel`, `AbundancesModel`.
- Runs go through `AZUREMainThread` / `AZUREMCMCThread`, which call the same
  engine (`AZUREMain`) as the CLI. Put no physics in the GUI. If the GUI needs
  a number, the engine computes it.
- Segment-line keyword tails are kept verbatim in `extraTokens`, but only for
  keywords listed in `include/SegmentKeywords.h`. A new engine keyword must be
  added there. `tests/consistency` checks the list against the engine and
  pyazr.
- `InTabDocs.cpp` holds the help text shown in each tab. Update it with the
  Sphinx user guide when a field's meaning changes.
- The Plot tab (`PlotTab`, `AZUREPlot`, Qwt) plots the files in `output/`, not
  the `.azr`. See the `azure2-gui-review` skill before showing a result to
  someone.

## Rules

- **Round trip is sacred.** Loading and saving an `.azr` without edits must
  give back the same file, including inactive lines, unknown tails, comments
  and full numeric precision (level lines are written at full precision since
  2026-10-17). Every new field gets a case in
  `tests/gui/segments_roundtrip_test.cpp` or the relevant `tests/gui/*_test.cpp`.
- Labels say the unit and frame, for example "Energy (MeV, lab)". A field
  label must not contradict the engine: "Mass" is an atomic mass in u, not a
  mass number.
- Qt5 only, for now. Homebrew disables `qt@5` and `qwt-qt5` on 2027-05-19, so
  new code avoids what Qt6 removed or moved (`QRegExp` → `QRegularExpression`,
  `qrand` → `QRandomGenerator`, `QString::SkipEmptyParts` →
  `Qt::SkipEmptyParts`; 65 uses remain). New connections use the
  pointer-to-member `connect`, which the compiler checks, not the
  `SIGNAL`/`SLOT` strings, which fail only at run time. Guard version-specific
  code with `QT_VERSION_CHECK`, as `PairKeyCombo.cpp` does.
- Code that needs Qwt is guarded by `#ifdef USE_QWT`. A `-DUSE_QWT=OFF` build
  must still compile, and so must `-DBUILD_GUI=OFF`.
- Never block the UI thread with a calculation. Long work goes through the run
  threads, which report back with signals.

## Tests

`tests/gui/<name>_test.cpp`, registered in `tests/gui/CMakeLists.txt` with
`QT_QPA_PLATFORM=offscreen`, so they run headless in CI. To check by eye, see
the `azure2-gui-review` skill. A visual check does not replace a round-trip
test.
