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
  the `.azr`. See the `azure2-plot` skill before showing a result to
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
the `azure2-plot` skill. A visual check does not replace a round-trip
test.

## Changing a dialog or a tab: check it headless before showing it

All in `tests/gui/`, all with `QT_QPA_PLATFORM=offscreen`:

- A ctest for the behaviour (`level_dialog_test.cpp` is the pattern: build the tab with its models, drive the
  widgets, assert). With `AZURE2_LEVEL_DIALOG_PNG=<stem>` it also saves the dialog as PNGs under `xvfb-run -a`, so
  the layout is looked at before a person sees it. `tests/gui/pair_key_combo_test` does the same for the pair-key
  combo (`AZURE2_PAIR_COMBO_PNG=<stem>` saves the data-segment dialog).
- `levels_compat_check a.azr b.azr ...` (not a ctest) opens real projects through the GUI reader and checks
  write-back, that the Add/Edit Level list offers every level in the file, and that regenerating channels from the
  pairs reproduces the file's set. Run it on one `.azr` from each of five reaction systems (skip `variants/` and
  `shard*`) and quote the result. Two write-back differences are expected and not caused by a change: older files
  hold `1 0 0` in the three placeholder columns (the GUI writes `0 0 0.0`), and pyazr-written files carry more than
  the 6 significant digits the GUI saves.
- `AZURE2_SCREENSHOT_DIR=<dir> xvfb-run -a -s "-screen 0 1400x900x24" AZURE2 x.azr < /dev/null` opens the project,
  saves every tab as `<dir>/tab<n>_<name>.png` and exits: the way to see the whole GUI without a display, for layout
  reviews and the docs' figures. `AZURE2_THEME=system` turns the azure theme (`gui/src/AzureTheme.cpp`) off for
  comparison.
- Only a GUI build is launched as the GUI (`ldd <binary> | grep -c Qt` > 0). A CLI-only build opens a text menu; with
  an open stdin under `nohup` it re-prompts forever and its log grows without bound (once to hundreds of GB, filling
  a group quota). Launch with `< /dev/null 2>&1 | head -c 200000 > log`; 5 s later check the log is small and
  `ps -eo pid,cmd -u $USER | grep [A]ZURE2` shows one GUI.
- Stop a GUI by pid and verify (`kill $pid; sleep 1; ps -p $pid` empty) before deleting its log. `pkill -f` and
  `pgrep -f` take regular expressions: the `+` in a name like `11B+a.azr` matches nothing. Never `rm` a log whose
  writer may be alive; `: > log` truncates safely.

J^π labels in the Levels and Pairs tables and the level dialog use the superscript signs U+207A/U+207B as plain text
(QComboBox items do not render rich text); the LaTeX export writes `^{+}` itself. The level dialog lists only the J^π
that give channels (`LevelsTab::spinParityOptions`, built by trying each candidate through `calculateChannels`), so
the list and the channels cannot disagree.

Pair keys are never typed in the GUI: the segment dialogs (main and component entrance/exit) and the Run tab's
reaction-rate box use `PairKeyCombo` (`gui/include/PairKeyCombo.h`), which lists the pairs by name ("2: ¹⁴N+n
[0.000 MeV]"; exit items as the reaction with the chosen entrance) and hands back the 1-based key the `.azr` stores.
It rebuilds its items from the pairs model whenever it drops down. A key beyond the pairs list is kept and marked
"(no such pair)".
