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

Code: `api/front/` (facade and test), `src/ProgramMain.cpp`, `include/ProgramObserver.h`, `gui/`.

## 11. The Levels and Channels tab and the Particle Pairs tab as easy as dev's, and quicker (U7, agent levels)

The PI (2026-11-19): "the level and channels tab should be easy as in the old version, with some QoL improvements
for easier navigations, removal, or addition on the go", and "too many things must be inserted by hand and are not
selectable". The measure is dev's program (`wt-devgui/gui`); the tabs below keep its layout (levels | channels |
configuration and details, the records of the level scheme below) and its every action (the parity table, 11.6), and
add what makes a level scheme of tens of levels quick to work through.

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
  emptied when a project is opened.

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
