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

**Tests** `tests/gui/gui_shell_test.cpp` (label `gui`): the theme's palettes and the View menu choice; the tool bar's
actions and their state; the status bar's file, modified mark and validity, the dock hidden on a valid project and
opened by an error, its double click; the start page and its recent and example entries; the options' editors made
from the facade's rows (a value set in an editor reaches the command line); the summary after a calculation of T1
(rows = data sets, total = the engine's); the plot's wheel zoom, pan, cursor read-out, check-box lists and the CSV
export (rows = points of the files).

Code: `api/front/` (facade and test), `src/ProgramMain.cpp`, `include/ProgramObserver.h`, `gui/`.
