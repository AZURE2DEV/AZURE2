# The Qt setup program on the upgrade engine: one boundary

The Qt setup program (`gui/`, AZURESetup) edits a project, runs the program's commands on it and shows what the run
wrote. On `dev` it met the engine in four places: its own reader and writer of the `.azr` text, a shared `Config`
filled by the tabs, the engine's `main` started in a thread with a plain `bool` as its stop flag, and the output files
read back by name. This note puts one narrow boundary there instead: a facade of the engine, `azure2::front`
(`api/front/include/azure2/front/Engine.h`), plain C++, the one engine header `gui/` includes. It states what the
facade does, what it calls, what it never does, and the test that holds it to the program.

## 1. The rule

1. `gui/` includes Qt and `azure2/front/Engine.h`, nothing else of the engine (no `Config.h`, no `AZUREMain.h`, no
   `io/` or `experiment/` header). The engine includes no Qt header. `Engine.h` holds no engine type: strings,
   vectors, numbers.
2. The facade calls the functions the command line and pyazr call, never a copy of them: the lossless document
   `io::Document` and the grammars of `<config>` (`io::ConfigEdits`), of `<levels>` (`io::LevelRow`,
   `io::EditLevelRow`) and of the record sections (`io::ReadRecords`, `io::WriteRecords`), which pyazr's `Project`
   edits through the binding; `ValidateProjectBytes` (what `AZURE2 --validate` prints); `convert::ConvertToFormat3`
   (`--convert`); the schema's JSON (`--schema`); and `RunProgram`, the program's run.
3. No state outside the objects a caller owns: a `front::Project` per open project, a `RunRequest`, an observer and a
   `std::atomic<bool>` stop request per run. No engine global is added, no `mutable`, no `const_cast`.

## 2. The project

`front::Project` is the file's `io::Document`. `Open` reads a format-3 file and refuses any other with the engine's
refusal (`io::NotFormat3Message`, which names the conversion); `ConvertProject` is `--convert` (the same converter,
nothing overwritten). Edits replace only the lines they change: `SetConfig` the entry's line, `SetLevels` the rows that
differ (an edited row keeps its layout, `io::EditLevelRow`), `SetRecord` the field lines whose value changed, or the
record written by `io::WriteRecords` when keys were added or removed. Every row is read back by the engine's own
`io::LevelRow` before it is accepted. A project opened and saved without an edit is written back byte for byte;
`Save` writes a temporary file and renames it. `Validate` returns `--validate`'s diagnostics as data (severity,
code, line, record, field, message): the program's text is not parsed.

The tabs keep their layout and names. The Pairs and Levels tabs read and write `<levels>` rows through `Levels` and
`SetLevels`; the tabs of records (Segments From Data, Segments Without Data, Experimental Effects, and the format-3
sections the upgrade added: model data, gamma decays, channel radii, isospin, potentials) build their forms from the
schema (`SchemaJson`: kinds, fields, titles, closed choices, units), as the web page does; `<config>` from
`ConfigTableJson` (the positional entries and keys of `io::ConfigKeys()` with the old program's titles); the sections
read by their own grammar (`<thm>`, `<parameterSettings>`, `<mcmc>`) are edited as text and checked by `Validate`.

## 3. A run

`front::Run(request, observer, stop)` is the program's run of one command, `AZURE2 <command> [options] FILE`, in the
calling thread: the CLI's `main` and the facade both call `RunProgram` (`src/ProgramMain.cpp`, the part of `main`
after the process's start-up), with a `Config` whose stream is the caller's. The facade's stream hands each line to
`RunObserver::Message` (a line ending in `\r`, the console's progress line, is marked transient). The run reads the
saved file, as the program does: the GUI saves before it runs.

The engine learns of the caller through one pointer, `Config::observer` (`include/ProgramObserver.h`, null for the
command line): its `stop` request, read by the fit (`FitControl::stop`, the objective's `EvaluationObserver::stop`), by
the set-up's channel rows and external-capture integrals (`Config::StopRequested`) and by the sampler's step hook; and
three sinks, each called from the run's threads: every `FitProgress` record, the final evaluation's chi^2 per data set
and the total. A stopped fit ends at its best point (`status stopped` in `fitStatus.out`), runs neither HESSE, MINOS nor
the band, and writes its outputs there. `RunResult` returns the exit code, whether the stop was set, the chi^2 (total
and per data set, the numbers the program prints) and the files the run wrote in the output directory.

The GUI runs `front::Run` in a worker `QThread`. The observer posts lines and progress to the window with queued
signals; the window disables the run's controls until the thread ends; the stop button sets the atomic. Nothing else
is shared between the threads: the worker holds copies of the request's strings, the project is saved before the
thread starts and not edited by the run. One run at a time (the OpenMP thread count is the calling thread's, set by
`--threads` inside `Run`).

## 4. What is dropped, and why

The controls of what the PI removed do not come back: fold windows and `useAdaptiveGrid` (every fold is on the shared
mesh), `--gsl-coul` and libcoul (one Coulomb library), the GSL fitter and the `numcmc/` sampler (being removed on
another branch; the fit tab lists the fitters of `Fitters()`, the commands of `RunCommands()`, so a fitter or
command gone from the engine is gone from the GUI). dev's Plot tab drew with Qwt and read the output
files itself: it is replaced by a plot drawn with QPainter (section 6), and `USE_QWT` is gone.

## 5. The Plot tab: the run's output files as data

`front::ReadCurves(project)` returns the curves of the last run: for every active data set and extrapolation, the
file the engine wrote it to and its block there, found by the engine's own build of the project (each record's pair
keys and observable, `experiment::ReadFormat3Text`) and by the writer's rule (`src/DataOutputFiles.cpp`: one file per
pair keys, one block per segment in segment order, `.band` siblings with the same blocks); a data set the run left out
is skipped through chiSquared.out's keys, and a block whose rows differ in number from chiSquared.out's N is refused
with a note, as is a file with fewer blocks than the project's records (the file is older than the project). Each
point carries E, E_x, theta, the calculated value and S column, the data and errors, and the band's 1 sigma. `gui/`
parses nothing: `PlotTab` turns the curves into series and `PlotWidget` draws them with QPainter (linear or log axes,
error bars, a single point as a marker, the band shaded, a legend, the residuals (data - calc)/error in a lower
panel, zoom by a dragged rectangle) and writes PNG, SVG (QtSvg) or PDF (QPdfWriter). No plotting library.

## 6. Fit parameters, the parameter settings and the sampler's priors

`front::ReadParameters(file)` sets up a session of the saved project as pyazr's `Session` does (`Config` read from
the file, the matrix route resolved, `AZUREAPI::Initialize`; `api/front` links `AZUREAPI` privately) and returns the
parameter vector of `param.sav`: name, kind and unit (`ParamKindTable`), start value in fit units, fixed flag, the
fit's limits (`GetParameterBounds`, `<parameterSettings>` applied), level number and channel.  It writes only what a
session writes (an `intEC.dat` it computes).  `front::ReadParameterFile` reads a `param.sav`/`param.par` as the
engine does ("name value error", the `#parametrization` line skipped).  `front::ParameterSettings` and
`SetParameterSetting` read and write the rows of `<parameterSettings>` in the grammar of
`ParameterLimitsManager::ReadParameterSettings` (9 tokens, or the GUI's 12- and 14-token "Level N ..." forms; a
row is replaced by name, appended, or removed); `front::SamplerPriors` and `SetSamplerPriors` the `<parameters>`
block of `<mcmc>` in the grammar of `src/SamplerRun.cpp`.  "Take the fitted values as the new start" is the engine's
`--parameter-file`: the GUI copies `param.sav` beside the project and passes it to the next run; it does not rewrite
the project (pyazr's `save_fit` does that, in Python).

The sampler's progress: `ProgramObserver::samplerStep` (`include/ProgramObserver.h`), called from the program's step
hook (`src/ProgramMain.cpp`, null on the command line), reaches the observer as a `RunProgress` with stage `mcmc`
(step, steps asked, mean log-posterior and log-likelihood).  `front::HasCommand("mcmc")` decides whether the MCMC
tab exists: a command removed from the engine's table removes the tab, without an edit of `gui/`.

## 7. `<thm>` as forms

`front::ThmKinds()` gives the engine's tables of the keys of `<thm>` (`ThmOptionKeys`, `ThmExperimentKeys`: the
tables its parser reads) as two kinds in the schema's shape; a key that applies only when another is given says so
(`FieldInfo::appliesWhen`, e.g. the `ff*` keys need `vertexModel`).  `front::ThmEntries` lists the key=value
entries as the reader takes them (`ThmBlockEntries`); `SetThmOption` and `SetThmExperimentKey` are the engine's line
edits (`ThmOptionEdits`, `ThmExperimentEdits`, `ApplyThmLineEdits`), the ones pyazr's `Project` calls: a value the
reader would refuse is refused with its message.  The THM tab draws its forms from these and keeps a text view.

## 8. The Tools menu: pyazr as a child process

ENDF-6 export and import, EXFOR and ENSDF imports, the Coulomb table and the web page live in pyazr (Python); the
GUI does not re-implement them.  `PythonTools` starts `<python> -m pyazr <command>` with `QProcess`, its output in
the log pane, the interpreter a setting of the computer (`QSettings`, Configure > Directories), every item disabled
with the reason when `import pyazr` fails.  The web page is `python -m pyazr ui FILE --json` (127.0.0.1; the URL
from its JSON line), stopped with the window.  The samplers (NUTS, emcee) and the band jobs are no commands of the
command line: a disabled item points to the web page and Python.

## 9. What stage 2 did not do

No Results page of the sampler (statistics of `samples.mcmc`); no per-curve style dialog, pan, tracker, print or
level markers in the plot; the parameter table edits limits and priors, not values or fixed flags; dev's EXFOR
search dialog, Coulomb-table builder dialog and THM diagnostics plots are not ported (Tools runs the pyazr commands).

## 10. Tests

- `api_front` (`api/front/tests/front_test.cpp`, no Qt, every build): T1, T2, T3 and 13N opened and saved byte for
  byte, also after every record, row and entry is set to what it holds; an edit changes one line and is undone; a
  broken field is named by its code; an old file refused and converted; `calculate` through the facade gives the same
  output files, byte for byte, as `AZURE2 calculate` on T1 and 13N, and the chi^2 the program prints; a fit stopped
  through the facade ends with `stopped` and its outputs.
- `gui_*` (label `gui`, registered only with `-DBUILD_GUI=ON`, `QT_QPA_PLATFORM=offscreen`): the window opens T1, T2,
  T3 and 13N and saves them byte for byte, runs `calculate` and `fit` through the worker thread, and grabs each tab
  into `AZURE2_QT_SHOTS`.  Stage 2 (2026-11-19) adds: the Plot tab after calculate on T1 (data points = rows of the
  output files, a grab not blank, zoom, PNG/SVG/PDF); extrapolate, fit and minos with the band, rate through the pair
  drop-downs and mcmc from its tab, each equal byte for byte to the program's run of the same command line; the
  parameter table, a `<parameterSettings>` row, the fitted start (= `--parameter-file`); the THM form; records added,
  duplicated, renamed and removed, a pair and levels edited, saved, reopened and valid; the Configure dialogs, the
  formalism, File > New, Open Recent and Convert through its dialogs; the Tools menu and the ENDF-6 export of T1.
- `api_front` also checks `ReadCurves` (points = rows of the files = chiSquared.out's N on T1 and 13N),
  `ReadParameters` and the `<parameterSettings>` and `<mcmc>` rows.

## 11. Stage 3, the shell: a modern, friendly window (agent shell, 2026-11-19)

The PI (2026-11-19): "the GUI could be much more modern. Make some CI for GUI compilation ... Try to get everything
green on CI."  This section is the window around the tabs: its look, its toolbar, status bar and start page, the
Calculate tab, the Plot tab and the CI job.  The tabs of records, Pairs and Levels belong to the agents "forms" and
"levels"; the look reaches their widgets through the application-wide style only.

**CI** (`.github/workflows/gui.yml`, a workflow of its own: the GUI needs a Qt the engine's legs do not install, and
its red or its time must not hide the engine's).  ubuntu-24.04 with Qt 6.4.2 from apt (`qt6-base-dev`,
`qt6-svg-dev`, `qt6-qpa-plugins` for the offscreen platform): configure with `-DBUILD_GUI=ON`, build everything (the
tests need `AZURE2` and pyazr's module), count the warnings of `gui/` and `api/front/` (a warning annotation, as
build.yml, which does not fail on a warning), run `api_front` and `ctest -L gui` offscreen, keep the screenshots as an
artifact.  A failure's first lines become `::error` annotations (`scripts/ci_first_errors.sh`).  macOS (Homebrew
`qt`) and Windows (MSYS2 `qt6-base`, `qt6-svg`) are written but untested, marked `continue-on-error`.  `gui/` keeps
to the API of Qt 6.4: a newer call (the colour scheme of Qt 6.5) is guarded with `QT_VERSION_CHECK`.

**Look** (`AzureTheme`).  Fusion on every platform with two palettes, light and dark, built from one table of
roles; the choice is a View menu item (System, Light, Dark) kept in `QSettings` ("System" follows the platform's
colour scheme where Qt >= 6.5 tells it, else the lightness of the platform's window colour).  One restrained style
sheet over the palette: spacing in em, no font size in pixels, readable table headers, alternating rows in every
table, a visible focus frame, the selection in the accent colour, primary buttons (Run, Save and Run) drawn in the
accent (`property primary=true`), dangerous ones (Stop) in red when enabled.  Icons for the main actions are SVGs drawn
for this program (`gui/images/actions/`, strokes only, one colour that reads on both palettes), scaled by Qt for any
device pixel ratio.

**Window** (`AZURESetup`).
- A tool bar: New, Open, Save, Validate, Calculate, Fit, Stop, Plot, each with its shortcut in the tool tip
  (Ctrl+N, Ctrl+O, Ctrl+S, F7, F5, F6, Esc, Ctrl+P).  Calculate and Fit set the Calculate tab's command and run it;
  Stop sets the run's stop request; Plot shows the Plot tab.  Disabled while they cannot act (no project, a run).
- A status bar: the project file and a modified mark; the validity ("valid", "2 errors, 1 warning"), which opens the
  Diagnostics when clicked; the running job (its command, a progress bar that is busy for a run without a known end
  and counts steps for the sampler, the last chi^2) and a Stop button.
- Diagnostics in a dock that is hidden while the project is valid with no warning and opens by itself when a
  validation finds a problem; a double click goes to the tab and record (`showRecord`).
- A start page while no project is open: New, Open, the recent projects, the example projects found beside the
  program or in a checkout (`tests/upgrade/T1`-`T3`, `13N`), and one line per tab saying what it is for.  Opening or
  creating a project shows the tabs; closing it (File > Close) the start page again.
- Unsaved changes are asked about before New, Open, Close and Quit; the window's geometry, its dock and tool bar,
  the last tab and the recent files are kept in `QSettings`; an `.azr` dropped on the window is opened; an old project
  offers the converter (`ConvertProject`) rather than refusing.

**Calculate** (`RunTab`).  Every option an editor of its kind, made from the facade's row (`front::Options`): a flag a
check box, a whole number or a number a spin box with the row's range, default and unit (`Option::unit`, new: the
unit the row's help names), a word a drop-down, a file a path with a chooser; under each its help in one line; an
option asked only with another (`Option::when`, new: "with --covariance-band") enabled only then.  Run and Stop as
the two primary buttons; while a fit runs a progress strip shows the evaluation, chi^2 and the best so far and a small
curve of the best chi^2 against the evaluation (from the observer's `RunProgress`).  The log highlights warnings and
errors (a line the program starts with "Warning", "WARNING", "Error", "ERROR" or "refused").  At the end a summary:
the total chi^2 and, per data set, its record's name, N, chi^2 and chi^2/N in a sortable table, with "Open Plot".

**Plot** (`PlotTab`, `PlotWidget`).  Wheel zoom about the cursor (Ctrl: x only, Shift: y only), drag a rectangle with
the left button, pan with the right or middle button (or Shift + left), double click or Reset to the whole range,
the cursor's coordinates in the corner, log/linear per axis, legend on/off, the residual panel; the data sets and
extrapolations as check-box lists with All and None; export PNG, SVG, PDF and CSV (the points drawn, one block per
series).  The plot is drawn in the window's palette on screen and on white in a file.

**Sampler results** (`MCMCTab`, a Results page): the chain of the last run read by the engine's one reader of
`samples.mcmc` (`ReadSampleFile`, moved out of an anonymous namespace of `api/src/AZUREAPISampler.cpp` and declared in
`AZUREAPISampler.h`, so that pyazr and the GUI read the file with the same code), through the facade's
`ReadSamplerChain`; the statistics (mean, standard deviation, 16/50/84 % quantiles) after a burn-in chosen there, the
acceptance fraction (moves that changed the position), and a trace, histogram or mean-log-posterior plot drawn by
`PlotWidget`.

**Facade additions** (`api/front/`, each with its case in `front_test.cpp`): `Option::unit`, `when`, `needs`,
`instead`, `lowOpen`; `SamplerChain` and `ReadSamplerChain`.

**Section 9's list, now**: the sampler's Results page, a per-curve colour and marker (the plot lists' context menu),
pan and the cursor's read-out are done.  Still not: print, level markers in the plot, a corner plot of two parameters,
the THM diagnostics plots, dev's EXFOR search and Coulomb-table dialogs; the macOS and Windows GUI jobs have not run.

**Tests** `tests/gui/gui_shell_test.cpp` (label `gui`): the theme's palettes and the View menu choice; the tool bar's
actions and their state; the status bar's file, modified mark and validity, the dock hidden on a valid project and
opened by an error, its double click; the start page and its recent and example entries; the options' editors made
from the facade's rows (a value set in an editor reaches the command line); the summary after a calculation of T1
(rows = data sets, total = the engine's); the plot's wheel zoom, pan, cursor read-out, check-box lists and the CSV
export (rows = points of the files).


## 12. Forms of values: what is chosen is chosen, not typed (upgrade 7, 2026-11-19)

The record tabs of stage 2 showed a record as one line of text and its fields as text lines in which the grammar was
typed by hand (`2e+17 atoms/cm2 fit prior 2e+17 10% grid_max 3e+17`).  The PI's measure is dev's Segments tab:
everything that can be chosen is chosen.  This section is the design of the forms and the tables that replace them.

### 11.1 Where the choices come from

`gui/` spells no keyword, unit, kind, quantity or grid word.  The facade adds, from `experiment/`'s table:

- `FieldInfo::value`, a `ValueInfo`: the value type and, recursively, the schema's `FormSpec` (alternatives led by a
  keyword, positional and keyword parts, each with its own `ValueInfo`: type, dimension, its units, enum words, the
  kinds a reference may name, positive / minimum, epoch, profiled, none).
- `ParseFieldValue(kind, key, text, quantity)`: the value split into its pieces (`ValueParts`: a nuisance's number,
  unit, fit flag, prior centre and width (and whether in per cent), grid_max; a form's alternative and its parts by
  name; a formula's expression and parameters; a grid's kind, span, unit, frame, step or points; the column roles
  and their units ...), by the engine's own parsers (`experiment/Parsers.h`), the ones `--validate` reads the field
  with.  A text the parser refuses gives its message, without the "dataset 'x', field 'y':" head, for the form to
  show beside the field.
- `FieldValueText(kind, key, parts)` / `ValueText(info, parts)`: the text from the pieces, written as the record
  writer writes tokens (`io::WriteTokens`) and a nuisance as `experiment::NuisanceText` writes it.  The words of
  the grid grammar (`GridKinds()`, `Frames()`) and the column roles (`ColumnRoles()`, `experiment::AxisRoles` and
  `PlainRoles`) are lists of the facade; front_test writes every grid kind and reads it back with the engine.
- `Quantities(kind)` (name, title, doc, units), `Pairs(project)` and `ReactionLabel(project, in, out)` (the
  converter's `experiment::ReactionLabel`), `ReactionPairs(project, reaction)` (the pair keys a one-term reaction
  names, by the engine's `PairsNamed`).
- `KindInfo::summaryFields` (the columns of a table: the required fields first, then the references to other
  records, the windows and grids, and, for a kind with at most eight fields set often, all of these; never free text
  or the column roles), `KindInfo::includeField` (an optional yes/no word whose default is its first: a data set's and
  an extrapolation's `active`), `FieldInfo::group` (the field this one belongs with: the one its condition names,
  e.g. a stopping's `phase` with `stopping`, a profile's `truncation` with `profile`).
- `Project::MoveRecord` (a record and its comment lines change places with a neighbour; no other line moves).

### 11.2 Which widget for which value

| Value type | Widget | Choices from |
|---|---|---|
| enum | drop-down (empty item for an optional field) | `ValueInfo::values` |
| reference | drop-down of the names of the records of those kinds in the project, "(none)" when optional; epoch drop-down for a target with epochs | `ValueInfo::values` (kinds), `KindInfo::sections`, the records |
| number, integer | line with a numeric validator (> 0, minimum) | `positive`, `min` |
| unit | drop-down | `units`, or `QuantityUnits` of the record's quantity |
| number with unit, range | number field(s) and a unit drop-down | `units` |
| nuisance | number, unit drop-down, check box "Fit", an expandable "Prior" group (centre, width, "%" check box, grid_max); "(record)" drop-down when it may name a shared record; "profiled" check box where allowed | `refKinds`, `profiled` |
| structured (form) | drop-down of the alternatives (with their doc), then one row per part, each its own widget; an optional part has a check box | `alternatives` |
| formula | expression line and a small table of parameters | |
| path | line and "Browse..." (relative to the project) | |
| grid | drop-down of the grid kinds, then lo, hi, unit, frame, step or points (+ log), or the list, or the file | `GridKinds()`, `Frames()` |
| columns | a row of columns, each a role drop-down and a unit drop-down for an axis role; "+"/"-"; "from the file" check box | `ColumnRoles()` |
| reaction | entrance and exit pair drop-downs (the pairs of `<levels>`; the exit pairs of the same compound nucleus, and total capture) | `Pairs`, `ReactionLabel`, `ReactionPairs` |
| quantity | drop-down of the quantities with their titles | `Quantities(kind)` |
| number list, composition, spin, nucleus, pair, text | a line | |

Each field has an "Edit as text" toggle: the text line with the grammar as placeholder.  A value the parser refuses or
the facade cannot split opens as text, with the parser's message under it.

### 11.3 The round-trip rule

A field's editor keeps the text it was given.  Until a widget of the field is changed, the editor returns that text
unchanged; after a change, the text is written from the pieces, and every piece not changed keeps its spelling (a
number keeps the text it was written with, `6.0e-20` stays `6.0e-20`).  Apply writes the record through
`Project::SetRecord`, which rewrites only the lines whose value changed.  So opening and applying an unedited record
leaves the project byte for byte the same; front_test checks that every field of every record of T1, T2, T3, 13N and
`section_examples.azr` splits and writes back to the same text, and gui_forms_test that Apply of every record
through the new forms moves no byte.

### 11.4 Checking as the user types

Each editor parses its text with `ParseFieldValue` on every change and shows the message under the field in plain
words (red).  Apply validates the project as it would be after the edit (`Validate` of a copy): an error of this
record is shown beside its field, and the user is asked before the record is written anyway (a '?' placeholder of a
new record is such an error until it is filled).

### 11.5 The tables

A section's records are a table (`RecordTableModel`): a check box column for the include field, the name, the kind
(when the section holds several), then the union of the summary fields of its kinds, each cell the field's short
text (`ShortValueText`: a nuisance "1 (fit)", a form its alternative and numbers), the whole record as the row's
tooltip.  Sortable, a filter line, multi-selection; Add (a dialog with the kind, the name and the required fields as
forms, pre-filled: the first reaction of the project's pairs, the first quantity and unit, norm 1 fixed), Remove,
Duplicate, Move up / down, include / exclude of the selected records; a context menu with the same, Ins, Del, Ctrl+D,
Ctrl+Up/Down.  A cell of an enum or reference field is a drop-down in place; a nuisance's "fit" a check box in place.
The Experimental Effects tab shows targets, beams and detectors as three tables, each with its form below; the form
groups a field with the fields of its `group` in a collapsible box, open when the record sets one of them, with the
field's doc as one line and its unit.  "Used by" lists the data sets that name the record, with "Assign to selected
data sets".

Code: `api/front/` (facade and test), `src/ProgramMain.cpp`, `include/ProgramObserver.h`, `gui/`.

## 13. The Levels and Channels tab and the Particle Pairs tab as easy as dev's, and quicker (U7, agent levels)

The PI (2026-11-19): "the level and channels tab should be easy as in the old version, with some QoL improvements
for easier navigations, removal, or addition on the go", and "too many things must be inserted by hand and are not
selectable". The measure is dev's program (`wt-devgui/gui`); the tabs below keep its layout (levels | channels |
configuration and details, the records of the level scheme below) and its every action (the parity table, 11.8), and
add what makes a level scheme of tens of levels quick to work through.  (Implemented on branch
`feature/upgrade-7-gui-levels`, 2026-11-19; parity in 11.8.)

### 11.1 What stays as it is

The tab holds the pairs, levels and channels in its three models (`PairsModel`, `LevelsModel`, `ChannelsModel`) and
writes them into `<levels>` only on a save or a run (`AZURESetup::commitTabs` -> `Project::SetLevels`: an unchanged
row keeps its bytes, an edited one its layout; section 2).  A project opened and saved without an edit is written
back byte for byte.  The channels of a level are never typed: they are the ones its J^pi, the pairs and the limits
(maximum l, maximum multipolarity, multipolarities per decay) give (`LevelsTab::calculateChannels`, the rule the
engine's `experiment::AllowedChannels` states as the Qt program's); a level whose J^pi changes gets the channels of
its new J^pi, the widths of the channels it keeps kept.

### 11.2 On the go: selection, editing in place, adding and removing

- **Selection.**  A project that opens selects its first level (in the order shown): its channels are on screen at
  once.  The selection is kept across every edit, a filter, a sort and an undo, by the level's identity (nucleus,
  J^pi, E), not by its row; a removed level hands the selection to its neighbour.
- **In place.**  Energy (MeV, a number; a level of the same J^pi and E is refused with its name), J^pi (an editor
  with J in a spin box of halves and the parity as a +/- toggle; a J^pi that gives the level no channel under the
  pairs and limits is refused, naming why), the include and fix boxes; in the channel table the value and its fix
  box.  The value column says what the number is and its unit, from the facade (11.4): Gamma (eV) of an open
  particle channel or a photon channel, ANC (fm^-1/2) of a closed one, the reduced-width amplitude (MeV^1/2) of a
  row with `gammaIsRWA`, mu (mu_N) and Q (b) of a static moment, B of a beta channel.  A refused value stays out
  of the model and the reason is shown at the cell (tool tip) and in the tab's status line.  Enter commits and moves
  down a row in the same column; Tab moves to the next editable column.
- **Adding.**  A bar under the levels: J (spin box of halves), parity (+/- toggle), E (MeV), the compound nucleus
  when there are several, and Add (Enter in any of its fields adds).  The J^pi offered start at the selected
  level's, so a series of levels of one J^pi is typed as energies only.  The Insert key puts the cursor in the bar's
  energy.  The dialog of dev (`AddLevelDialog`, "+") stays for those who prefer it.  Ctrl+D duplicates the selected
  levels: same J^pi, nucleus, fix flags, widths and channel fix flags, E + 10 keV (the engine joins two levels of one
  J^pi closer than 1 keV), the copy selected and its energy opened for editing.
- **Removing.**  Del (or "-", or the context menu) removes every selected level after one confirmation that names
  them (J^pi and E, and how many channels go with them).
- **Undo / redo.**  Ctrl+Z, Ctrl+Shift+Z (and Ctrl+Y) undo and redo every change of this tab and of the Particle Pairs
  tab: an add, a removal, a duplicate, an edit in place, a bulk fix, a change of the limits, a pair added, edited or
  removed, a level import.  An entry of the stack is the tab's state before and after (the three models and the
  limits, `LevelsUndo`); the project text is not touched until a save, as for any edit of the tab.  The stack is
  emptied when a project is opened.  Undone back to the opened state with no save in between, the project is
  saved byte for byte as opened; after a save, the rows the save rewrote come back with their content but not
  always their blanks (`SetLevels` writes a row over the line at its place).

### 11.3 Navigation

- **Filter bar** over the levels: a search line (matches J^pi as written, "3/2-", and the energy's digits), a J^pi
  drop-down that lists every J^pi of the levels with its count ("3/2- (4)"), "Only included", and the nucleus
  drop-down of a project with several compound nuclei.  The count of levels shown and in all stands beside them.
- **Sort**: by a click on a header (energy, J^pi); J^pi sorts by J, then parity, then energy, so the levels of one
  J^pi stand together.  The selected level's J^pi group is marked (its other levels' J^pi in bold).
- **Keyboard**: arrows move in the levels; Ctrl+Right (or Tab out of the last column) goes to the channels of the
  selected level, Ctrl+Left back; Space toggles the box under the cursor; F2 or typing edits.
- **Context menu** (levels): add a level of this J^pi, duplicate, remove, edit in the dialog, include / leave out,
  fix / free the energies, fix / free the widths (of every channel of the selected levels), show only this J^pi, add
  a background pole of this J^pi.  Channels: fix / free, set to 0, copy the value.  The buttons Fix All Widths and
  Fix All Energies of dev stay, and act on the selected levels when more than one is selected.
- **Background poles and missing J^pi.**  The red text of stage 2 named the J^pi without a level and said nothing
  of what to do.  It now says: "No active level carries these J^pi, and 'Hard-sphere scattering in every J^pi' is
  off, so they are left out of the calculation (no hard-sphere phase): switch the option on, or add a level of that
  J^pi (a click on one fills the add bar)."  Each J^pi is a link that fills the add bar.  A background pole is an
  ordinary level far from the data, usually with its energy fixed; "Add background pole of this J^pi" puts one 10 MeV
  above the highest level of the J^pi, energy fixed, widths 0, its energy opened for editing (the 10 MeV is a starting
  value to be edited, not a rule).

### 11.4 The channels and the details pane

The channel table: fix box, pair in nuclear notation (`p + 15N`, `γ + 16O(6.13)`), s, l (or the multipole, E1, M1,
GT, F), the value with its meaning and unit (11.2), and, for a particle channel, the reduced width as a fraction of
the Wigner limit, theta^2 = gamma^2 / gamma_W^2.  The details pane: the pair's particles, threshold and channel
radius (the channel's own radius from `<channelRadii>` when a record names it, "own; pair a" beside it, as dev did
with its radii table), E_c = E_x - S - E_x(heavy) - E_x(light), and at the level energy the penetrability P (or the
bound-state factor below threshold), the shift function S and dS/dE, the Wigner limit gamma_W^2 = 3 hbar^2 / (2 mu
a^2), the single-particle width Gamma_sp = 2 P gamma_W^2, the observed reduced width gamma^2 = Gamma / (2 P) (or C^2
/ (2 P) for an ANC) and theta^2.  Every number is the facade's (11.5), from the engine's own functions; the GUI only
shows them.  gamma^2 = Gamma / (2 P) is the observed reduced width (the thin-level value); Brune's gamma~^2 differs
from it by the factor 1 + sum gamma^2 dS/dE (level_parameters.md section 2), which needs all the level's widths and
is the engine's report (`parameters.out`), not this pane's.

### 11.5 Facade additions (`api/front/`, file `src/FrontLevels.cpp`, cases in `front_test.cpp`)

- `LevelColumnNames()`, `LevelColumnIndex(name)`: the columns of a `<levels>` row (`io::LevelColumnNames`), so that
  `gui/` names a column by its name and types no position.
- `ChannelAtLevel(row)`: what the value of a `<levels>` row is (kind, symbol, unit, in words) and the channel's
  quantities at the row's level energy: E_c, P, S, dS/dE (the engine's `ParticleChannelValues` on the pair the
  engine builds from the row, `PPair(LevelRow)`, as the binding's level doors do), gamma_W^2
  (`physics::wigner_limit`), Gamma_sp (`physics::single_particle_width`), gamma^2 and theta^2.
- `ChannelRadius(project, key, l, 2s, 2J, pairRadius)`: the channel's radius from the project's `<channelRadii>`
  (`experiment::ChannelRadiusFor` on the engine's build of the records), with the record that gives it.
- `PairUses(project, key)`: the data sets and extrapolations whose reaction names the pair (the engine's build:
  entrance and exit keys), for the confirmation of a pair's removal.

### 11.6 The Particle Pairs tab

Columns with their units in the header (J^pi of each particle in one column, masses in u, energies in MeV, radius
in fm), the pair in nuclear notation and its type as words (particle, gamma, beta).  Edit in place: separation
energy, excitation energies, channel radius (numbers, refused in the cell with the reason).  Duplicate (Ctrl+D: the
dialog opens on a copy, as a new pair must differ), remove with a confirmation that names what goes with it: the
levels' channels of the pair (the tab's own channels), the data sets and extrapolations that name the pair
(`PairUses`; the engine refuses them at validation once the pair is gone, so the confirmation says so).  Reordering
is not offered: the pair key is the row, and records name pairs by key.  The add/edit dialog: each particle by
element (drop-down of symbols) and mass number (spin box), Z, mass and spin filled from the mass table, J by a spin
box of halves, parity by a drop-down, the pair type by a drop-down.

### 11.7 Tests (`tests/gui/gui_levels_test.cpp`, label gui)

Each action driven through the widgets (the add bar's fields and Enter, the Del key, Ctrl+D, Ctrl+Z, a cell's
editor), then the project text read back: open T1, 13N and T3 and save without an edit (byte for byte); the first
level selected on open and its channels shown; a level added from the bar (its rows in the saved file); a J^pi and an
energy edited in place (rows read back; an invalid energy refused with the model unchanged); a width edited in place;
two levels removed with Del after one confirmation; duplicate; undo of each back to the opened text byte for byte and
redo; the filter, the search and "only included" (rows shown); bulk fix of the selection; the details pane's numbers
equal to the facade's; the channel value's unit by kind; a pair's cell edited in place; a pair's removal message
naming its data sets.

### 11.8 Parity with dev's two tabs, item by item

"dev" is `wt-devgui/gui` (the PI's measure); "stage 2" the program before this section; "now" after it.  Open rows
are marked **open**.

| dev's action or display | Stage 2 | Now |
|---|---|---|
| Pairs: table of light/heavy particle, spins, excitation, separation energy, radius, light excitation | same columns | same; header units (11.6) |
| Pairs: "+" add dialog (light and heavy particle by isotope text "12C", J, parity, Z, M, excitation energies, separation energy, radius, pair type, E1/E2 external capture) | same dialog | same dialog; the isotope fields complete from the mass table (p, d, t, a and every AME2020 nuclide: chosen, not typed). **Open**: J by a spin box of halves and the isotope as element drop-down + A spin box (11.6) are not done; J is still a validated text field (0, 0.5, 1, ...) |
| Pairs: duplicate a pair (11.6) | none | **open** (not done) |
| Pairs: "-" remove; first pair must be a particle pair; no photon entrance | same | same, plus a confirmation naming the channels and the data sets that use the pair (`PairUses`) |
| Pairs: double click edits the pair in the dialog | same | same; simple cells also in place |
| Pairs: "Masses from AME2020..." (compare and fill every pair, MassCheckDialog) | same | same |
| Pairs: "Fill masses and S" in the dialog from the mass table by Z and A | same | same |
| Pairs: compound-nucleus column when several nuclei | same | same |
| Pairs: per-channel radii table (every particle channel of the active levels, own radius bold, "Clear own radius") | the generic records panel of `<channelRadii>` (forms agent's `RecordsPanel`) | records panel unchanged; the channel's own radius shown in the Levels tab's details ("own; pair a"), from `ChannelRadius` |
| Pairs: isotopic abundances of samples (natural targets) | moved to the target records' `isotope` field (Experimental Effects tab) | unchanged (not this tab in format 3) |
| Levels: table Include?/Fix?/J^pi/Energy, compound nucleus column, nucleus filter | same | same, plus search, J^pi filter with counts, only included, count, J^pi group in bold |
| Levels: "+" dialog (J^pi from a list of those with channels, energy, nucleus) | same | same, plus the add bar, Insert, Ctrl+D, background pole |
| Levels: "-" remove one level | same | Del / "-" / menu for several, one confirmation |
| Levels: double click edits J^pi and E in the dialog | same | in place (J spin box + parity toggle, E); the dialog in the context menu |
| Levels: "NDS" import of ENSDF levels from the IAEA (IAEALevelDialog: J^pi choices for ambiguous assignments) | same | same (an undo entry) |
| Levels: Fix All Widths, Fix All Energies | same | same; on the selection when several are selected; fix/free both ways in the menu |
| Levels: Export LaTeX (clipboard or file) | same | same |
| Levels: header click on Include?/Fix? does not sort (LevelsHeaderView) | same | same |
| Channels: Fix?, pair, s, l (E1/M1/GT/F) of the selected level, bold when the width is not 0 | same; empty until a level is clicked | the first level selected on opening; value, meaning and unit, theta^2 columns; value in place |
| Configuration: maximum l, multipolarity, multipolarities per decay | same, plus hard-sphere option and missing-J^pi text | the text says what to do, links fill the add bar |
| Details: level, pair, particles' spins, Z, M, g, excitation, separation energy, channel radius (own from the radii table) | same, pair radius only | same, own radius from `<channelRadii>`; E_c, P, S, dS/dE, gamma_W^2, Gamma_sp, gamma^2, theta^2 from the facade |
| Details: reduced width / ANC / mu / Q / B box with its unit | same (kind decided in `gui/`) | same box; the kind from the facade (`ChannelAtLevel`) |
| Details: "Wigner Limit" button (Gamma_sp of the channel, eV) | same, `SingleParticleWidth` | shown at once (the facade's Gamma_sp at the level energy) |
| Undo | none | Ctrl+Z / Ctrl+Shift+Z over both tabs |
| Pairs: the pair type as a column | hidden | shown in words |
| Levels: the Insert/Del/Ctrl+D keys, the context menus | none | as 11.2, 11.3 |
## 14. The Physics tab's facade (2026-11-19)

The calls a Physics tab needs are in `Engine.h`, section "Physics": the rows of the physics table and one evaluated
(`PhysicsFunctions`, `EvaluatePhysics`: the rows `pyazr.physics` reads, now plain C++ in `api/physics_table/`), and a
`PhysicsSession` of a saved project (pairs, J^pi groups, the levels' physical energies and widths, the session's
channel functions, one J^pi group's U, T, phase shifts, |A| and det A^-1 on an E_x grid, the THM vertex, line shape
and distortion), each answer a `PhysicsCurves` of series over one axis.  The design, the engine function behind
every curve, the shapes and units, and which calls are slow: `docs/upgrade/physics_tab.md` (section 7 is the
interface).  Test: `api_front`, section (f).

## 15. The Overview page: the evaluation and its fit at a glance (U8, agent overview, 2026-11-20)

The PI (2026-11-20): "a first 'welcome' page with a brief overview of the current fit/status of the evaluation could
be nice to have."  The start page (section 11) is shown while no project is open; this page is the first tab of an
open project, selected when a project opens.  It answers, without typing, "what is this evaluation and where does
the fit stand?".  It reads; it never edits the project and never runs the engine on the GUI thread.

**What it shows** (cards in a scroll area; three columns from about 1300 px of width, two from 850, else one):

1. *What to look at next* (the first card, the whole width): plain sentences made from the other cards' numbers
   only, no physics judgement: the project does not validate (errors, warnings), there is no run yet, the outputs
   are older than the project file, the fit ended with a status other than `converged` (the engine's word), N data
   sets have chi^2/N above 3, N normalizations are more than 2 sigma from their prior, N parameters ended at a limit,
   N free parameters have no error after a fit, N active data sets were left out by the run (no point in their
   window).  Each is a link to the tab and record concerned (Diagnostics, Calculate, Segments and the record,
   Fitting Settings).  Quick actions: Validate, Calculate, Fit, Open Plot (the window's own actions).
2. *The project*: the file; the compound nuclei; the pairs in nuclear notation with their thresholds (E_x of the
   compound nucleus, MeV); the levels per J^pi (a row of bars); data sets active / all and extrapolations active /
   all; the records of the experimental effects (`<targets>`: targets, beams, detectors) and how many of them the
   active data sets name; the free parameters by kind (energy, width, norm, shift, apparatus, ...); the basis
   (Brune, formal (Lane-Thomas), Park) and A- or R-matrix; external capture on (which pairs) or off; every option of
   the calculation (`io::ConfigKeys()`' option rows) the project sets away from its default; the fitter chosen in
   Fitting Settings; the validity (a click opens the Diagnostics).
3. *The last run*: the command (fit or calculate, inferred from the files, below), when, a stale mark when the project
   file is newer than the run's outputs ("the project changed since this run"), the fit's status in the engine's
   words (`fitStatus.out`: status, minimizer, calls), and four numbers: total chi^2 of the data (`chiSquared.out`'s
   `Total-Chi-Squared`), the number of points, the free parameters, chi^2 per degree of freedom chi^2/(N - p); the
   objective's other totals of the file (norms, priors, ...) in one line.  A sortable table per data set: name,
   reaction, N, chi^2, chi^2/N with a horizontal bar (the bar's length chi^2/N over the largest), the norm the run used
   and its pull when the norm has a prior; the three worst rows and every row above 3 highlighted.  A click on a row
   opens the data set in the Plot tab (`PlotTab::selectRecords`).  Parameters at a limit or with no error, by name.
   Without a run: "No run yet" and the two buttons Calculate and Fit.
4. *Charts*: chi^2/N per data set as bars (a line at 1), and, after a fit run from this window, the best chi^2
   against the evaluation (the run's progress records, `RunProgress`, drawn by the Calculate tab's `ChiSquaredTrace`):
   the engine writes no file of a fit's progress, so a fit of an earlier session has no curve (the card says so).

**The facade** (`api/front/src/FrontOverview.cpp`, declared in `Engine.h` under "The Overview page"; plain C++):

- `ProjectSummary SummarizeProject(const Project &)`: from the text as it stands and the engine's build of it
  (`experiment::ReadFormat3Text`, the reader `--validate` and `ReadCurves` use; no set-up of the levels): the pairs
  (`Pairs`), the compound nuclei and the levels per J^pi group (`experiment::LevelTablesOfRows`: the engine's own
  grouping), the data sets and extrapolations (active from the build's specs), the effect records by kind and how many
  of them an active data set names (a field whose value is the record's name), the basis and the matrix from `<config>`
  (`use_brune`, `use_park`, the positional `amatrix`), the pairs with external capture on an active row
  (`LevelTables::externalCapturePairs`), and the option rows (`io::ConfigKeys()`) with the value the project gives or
  the default.  Milliseconds: called on the GUI thread.
- `RunSummary ReadLastRun(const Project &)`: the files of the project's output directory as they stand.
  `chiSquared.out` (the engine's writer, `WorkingCopy::WriteChiSquaredFile`: one line per data set the run used with
  key, chi^2, N, the norm used and the norm's prior chi^2, then the total line of `label: value` fields, every field
  kept in its order), matched to the records by key through the build (as `ReadCurves`); `fitStatus.out` (key / value
  lines, `WriteFitStatus`); the modification times of the output files and of the project file, in seconds since 1970.
  The pull of a norm with a prior is sign(norm - centre) sqrt(norm chi^2), the centre being the build's `normCentre`
  (or the start value): the number the file gives, signed.  Milliseconds.
- `bool ReadRunParameters(projectFile, outputDirectory, RunParameters &, why)`: the parameters of the saved project
  as the engine sets them up (`ReadParameters`: kinds, fixed flags, limits), the free ones counted by kind, and,
  read against the run's parameter file (`param.sav` after a fit, else `param.par`, `ReadParameterFile`), those that
  ended at a limit (|value - limit| <= 1e-6 max(1, |limit|), in the fit's units) or, after a fit, with error 0.
  It sets up a session: up to seconds for a project with external capture.  The page runs it in a worker thread
  (`QThread::create`), not while a run is in progress, the newest request only.
- `double ChiSquaredPerDegreeOfFreedom(const RunSummary &, const RunParameters &)`: chi^2/(N - p), NaN when N <= p.
- `std::vector<InputFile> InputFiles(const Project &)`: the files the active records read (a data set's data and
  covariance files, an extrapolation's grid file), resolved as the engine opens them; `ReadLastRun` names the
  records whose file is newer than the run (`changedInputs`, stale too; `projectNewer` for the project file).

**Which command made the files.**  The program does not write its command into the output directory, and the
facade's run must keep writing the program's files byte for byte (`api_front` (d)).  The rule, from the order in
which a run writes: a fit removes `fitStatus.out` at its start and writes it after the fitter ends, before the final
evaluation writes `chiSquared.out`; every run without `--parameter-file` writes `param.par` at its set-up, first.
So the files are a fit's when `fitStatus.out` is there, not newer than `chiSquared.out`, and not older than
`param.par` (or there is no `param.par`); else a calculation's; with only `AZUREOut_*.extrap` files, an
extrapolation's.  A calculation with `--parameter-file` after a fit writes no `param.par` and is taken for the fit
(stated on the page as "inferred from the files").

**Stale.**  The outputs are stale when the project file's modification time is later than `chiSquared.out`'s (or,
for an extrapolation, the newest `.extrap`'s).  An edit that is not saved does not make them stale: the run reads the
saved file, and the page says "modified, not saved" beside it.

**When it refreshes** (the one block of `AZURESetup.cpp` that adds the tab): a project opened (and the tab selected),
the end of a run, a save (the Save and Save As actions), an edit of a section (`ProjectHolder::sectionEdited`, after
0.2 s of quiet), a validation (`DiagnosticsPanel::validated`), and the page shown (after `commitTabs`, so that an edit
of the Levels tab is counted).  A refresh reads files only; the parameters' worker starts when the page is shown.

**Not done**: the points of an inactive data set (the data file would have to be read through the set-up); a
fit's convergence curve from an earlier session (no file); a per-data-set chi^2 of an extrapolation-only run (none
exists).

**Tests** `tests/gui/gui_overview_test.cpp` (label `gui`): T1 and T2 opened: the cards' numbers equal the facade's
(pairs, levels per J^pi, data sets, options); a calculation run through the window: the page shows its total chi^2 and
one row per data set with the chi^2 of the run's `RunResult`; an edit and a save: the stale mark appears; a project
with no outputs: "No run yet" and the two buttons.  `api_front`: `OverviewChecks` (a calculation and a fit run in the
test, and the program's own files of a run of `AZURE2`).

## 16. Changes on disk: the window follows its files (U8, agent overview, 2026-11-20)

The PI (2026-11-20, passed on by the director): "I just want the GUI to auto update if something in the files
changes."  A project is edited by more than this window: a text editor, pyazr, a script, the command line's runs.
The window follows its files; it never overwrites either side without a click.

**What is watched** (`ProjectWatcher`, `gui/src/ProjectWatcher.cpp`, a `QFileSystemWatcher`): the project file and its
directory, the input files of the active records (the facade's `InputFiles`), the parameter file of the last run
started (`--parameter-file=`), the output directory and in it the files the window shows: `chiSquared.out`,
`fitStatus.out`, `AZUREOut_*`, `param.sav`, `samples.mcmc`.  A file written to a temporary name and renamed (as
`Project::Save` and pyazr write) loses its watch: every path is added again after each report.  A change is reported
once every watched file has kept its size and time for a quiet time (0.4 s), so a file being written is never read
half-way; an output file the facade still cannot match (a block count against `chiSquared.out`) only gives a note,
never an error, and is read again at the next change.

**The project file** is compared by content (SHA-1), not by time: bytes equal to the window's own text (its own save,
detected at the notification, before another edit can come) or to the bytes last seen are no change.  A real change:
- no unsaved edit in the window: the project is read again at once; the tab shown and the record selected in each
  records panel are kept where the record still exists; the status bar says "Project reloaded: changed on disk";
- unsaved edits: neither side is overwritten.  A bar under the tool bar, not modal: "The project changed on disk,
  and this window has unsaved changes", with **Reload (discard mine)**, **Keep mine** (the bytes on disk are taken as
  seen: not asked again until they change again; the next Save writes the window's version) and **Show differences**
  (a plain line diff, `ProjectWatcher::lineDiff`, of the disk's text against the window's, in a dialog that is not
  modal either).

**Outputs** of a run outside the window: the Overview, the Plot tab (its check boxes kept, and its zoom when the user
had zoomed: the view is restored after the new curves are read) and the Calculate tab's list of data set names read
again; the status bar says so.  While the window's own run is in progress nothing is reported (its files are being
written); its end refreshes the window as before, and the watcher takes the outputs as they then stand as seen.
**A data file** rewritten: the Overview's stale mark names the data set (`changedInputs`) and its step links to the
record; the plot reads its files again (the data drawn are the last run's copy until the next run).

**View > Follow Changes on Disk**: on by default, kept in `QSettings` (`view/followDisk`); off, nothing is watched.

**Hooks**: the same block of `AZURESetup.cpp` as the Overview's (the watcher follows the project on open, edit and
save; the run's start and end).  `ProjectHolder` is not changed.

**Not done**: the Levels tab's selected level and the scroll positions of tables are not restored after a reload (the
records panels' selections are); the Calculate tab's chi^2 summary is that of the window's last run (the Overview
shows an outside run's); a change of a file the records name but the engine's build does not list (a stopping table, a
beam profile file) is not watched.

**Tests** (`tests/gui/gui_overview_test.cpp`, part (f)): the window's own Save causes no reload; the project edited
outside with no unsaved edit is re-read (the Segments tab still shown, the new norm in its table, no click); with
unsaved edits the bar appears, the window keeps its value and the disk its own, Show differences lists both lines,
Keep mine hides the bar and changes nothing, a Save then writes the window's version without a reload; a run of the
program outside changes the Overview's chi^2 without a click; a data file rewritten makes the stale mark name its data
set; following off, an outside edit is not read and the choice is kept.
