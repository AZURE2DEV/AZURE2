# The web UI of AZURE2: design

This note designs the graphical front end that replaces the Qt GUI (out of the tree since 2026-11-07): `pyazr ui
project.azr`, a FastAPI server on the user's machine and one static page in the browser, as fixed in outline by
`ARCHITECTURE.md` §8 and §11 and `ROADMAP.md` Phase 9. It is built on the objects of pyazr 3.0 (`pyazr3.md`, whose
decisions the PI took on 2026-11-13) and says what can be built on today's pyazr. Every number of section 2 comes from
a script of `work/webui/` (outside the repository; list in 2.8), run on the build `work/b-gecm` at one thread.
"pyazr step N" below is step N of the order of work of `pyazr3.md` section 4; "pyazr decision N" is decision N
of its section 5.

## 0. In one page

- **What it is.** One server process holds the project text and its history (undo = the previous text); the engine
  runs in worker processes: one session worker per open project for evaluations, one process per long job (fit,
  MINOS, sampler, bands, rate). The page is plain JavaScript modules with plotly.js vendored; no build step.
- **The project is the state.** Every edit is an edit of the project through `Project`'s grammar, validated by the
  engine; the session and the evaluation are derived and rebuilt when the text changes. Only the values of free
  parameters move without a rebuild (the fit vector `x`, decision 10 of `pyazr3.md`).
- **One-to-one API.** A route table maps each route to one pyazr call; a consistency test holds the table and pyazr's
  public objects together. The same JSON API serves the page and an AI agent (`/api/describe`, OpenAPI).
- **Generated, not written.** Forms of records come from `pyazr.schema()`, data kinds and parameter kinds from
  `pyazr.kinds()`, physics functions from the physics table. A new record kind, data kind, parameter kind or physics
  function appears in the UI with no edit of the UI (section 10).
- **Measured.** Opening T2 or 13N through the API takes 0.9 s and 1.9 s (13.1–14.4 s the first time on 13N, when the
  external-capture integrals are computed); an evaluation through the API costs 5–13 ms more than in the engine;
  everything the Results tab needs is 27 kB of JSON on T2 and 48 kB on 13N; a stop request ends a fit within one
  objective call (0.01–2.1 s), a kill within 0.1 s. Today the engine's own fit cannot be stopped short of killing it,
  and a killed program loses everything since its last 100-evaluation snapshot.
- **What pyazr 3.0's design lacks for the UI**: twelve items (section 12); the largest are a session from the project
  text without saving the user's file, progress and stop of the engine's fitters, writing a parameter vector into the
  project text, adding and removing pairs with their channels, and the reaction rate.
- **Order.** Six steps, each usable; the first three run on today's pyazr behind one adapter file (section 11). A
  first version for the owner's meeting: open, edit, validate, evaluate, undo, a recipe fit with progress and stop,
  and the χ² per data set, with plots if plotly.js may be vendored.

## 1. What the UI must do

### 1.1 The Qt parity list by tab

Each of the 201 items of `qt_parity.md` is assigned to a tab and a state (`work/webui/parity_classify.py`, output
`parity_class.tsv` with the call or reason per item; the counts below are its output):

- **E** exists in today's pyazr (a call serves it, possibly through the generic record API);
- **D** designed in `pyazr3.md`, not built yet;
- **M** missing from pyazr 3.0's design: a finding (section 12);
- **P** the page alone, no pyazr call needed;
- **X** dropped: format-1/2 grammar the converter handles, Qt-only behaviour, or replaced by a format-3 property
  (named records, a lossless document); the reason is in `parity_class.tsv` and section 13, decision 5.

| tab | E | D | M | P | X | all |
|---|---|---|---|---|---|---|
| Project | 16 | 1 | 2 | 1 | 4 | 24 |
| Nucleus | 30 | 0 | 11 | 1 | 2 | 44 |
| Data | 49 | 0 | 1 | 1 | 8 | 59 |
| Fit | 5 | 12 | 4 | 5 | 5 | 31 |
| Results | 1 | 2 | 1 | 7 | 1 | 12 |
| THM | 14 | 0 | 1 | 3 | 0 | 18 |
| Shell (title, help, theme) | 0 | 2 | 0 | 11 | 0 | 13 |
| all | 115 | 17 | 20 | 29 | 20 | 201 |

The groups and the calls that serve them (group, its parity sections → tab: call; state):

- `<config>`, checks, directories, runtime options (1.1, 3) → Project: `Project.config`, `set_config`; session options
  from `kinds()["options"]`. E (options D).
- Round trip, comments, number format (1.3) → Project: the lossless `Document` under `Project`; edits change only
  their lines. E.
- New project (2) → Project: `Project.new(path, nuclear_part=...)` needs a source project. M (an empty skeleton).
- Save before a run, the snapshot file (1.3, 14) → Project: `Project.session()` saves the user's file first. M (a
  session from the text).
- Pairs: table, masses, AME2020 checks (4) → Nucleus: `Project.pairs`, `set_pair`; `sources.masses` (today
  `masses`). E. As built (11.14): the table in the old program's columns, particles by name (`kinds()["pairs"]`);
  in the first view only those columns (11.16, `pyazr/ui/first_view.py`).
- Pairs: add, delete, channel regeneration (4) → Nucleus: no call. M.
- Levels and channels: table, edit, fix, Wigner limit (5) → Nucleus: `levels`, `add_level`, `remove_level`,
  `set_level`, `set_channel`; `Parameters.units`, `physics.wigner_limit`. E.
- J^π that give channels, channel generation, missing J^π (5) → Nucleus: no call (`add_level` takes the channels as
  given). M.
- IAEA levels (5) → Nucleus: `sources.ensdf` (today `nds.fetch_levels`). E.
- Hybrid potential (8) → Nucleus: `potential` records (`add_potential`, `set`, `remove`). E.
- Data sets, groups, samples, UPOS, sums, scattering lengths (6.1) → Data: `dataset`, `normalization`,
  `energy_scale`, `datum`, `target` records; `validate`; `describe`. E. As built (11.14): "Segments From Data" as a
  list of data sets in words (reaction, file, range, points, normalization, fitted), the data file chosen among the
  folder's files with the engine's reading previewed (the old "Choose…").
- Segments without data (6.2) → Data: `extrapolation` records. E.
- EXFOR (6.3) → Data: `sources.exfor` (today `nds.search_exfor`, `fetch_exfor`, `to_dataset`). E.
- Experimental effects (7) → Data: `target`, `beam`, `detector` records; `pyazr.stopping`. E (the screening
  potential's adiabatic limit M).
- Fitting settings table (9) → Fit: `Parameters` (`of_kind`, `names`, `units`, `bounds`) and a fit's errors. D.
- Edits of fitted values back into the records (9) → Fit: `save_fit` writes files only. M.
- Calculation types, minimisers, MINOS, band (10) → Fit: `Session.evaluate`, `extrapolations`, `fit_with_engine`
  (pyazr decision 8), `pyazr.bands`. E (the engine's fitters and MINOS as jobs, section 11.6; the bands, 11.7).
- Stop and the run log (10) → Fit: `fit_with_engine(progress=, stop=)`, the fitter's lines in the job's log. E
  (section 11.6).
- Reaction rate (10) → Results: no call. M.
- MCMC (11) → Fit: the emcee recipe on `Objective.log_prob` (`pyazr3.md` 2.6), as a job. E (section 11.7).
- Plots (12) → Results: `Evaluation.data[i]`; `physics.sfactor_conversion`; `pyazr.plot` for PDF. E, D, P. As built
  (11.16): the model of each data set on a fine grid over its range (`Session.model_curve`), the points as points,
  the residuals at the points (the old program plotted the model on its own grid as well).
- THM workspace (13) → THM: `Project.thm`, `set_thm`, `set_thm_experiment`, `remove_thm_experiment`;
  `thm_experiments`, `thm_background`, `thm_vertex`, `thm_lineshape`, `thm_distortion`. E (the THM tab, section 11.8).
- THM binding energy and amplitude flag per channel (13) → THM: `set_pair` has neither. M (shown read-only from the
  session, section 11.8).
- Title, About, help, theme, dialogs, screenshots (2, 14) → Shell: `pyazr.describe()`; the page. D, P (the icon,
  About with the citation and the versions: section 7; the first view in the old program's words: section 11.11).

### 1.2 What the engine does that Qt never showed

Each gets a place (capability → where; generated from), and all but the plots come from the engine's tables:

- Format-3 apparatus records: `target` (layers, epochs, isotopes, stopping evaluation), `beam`, `detector` → Data ›
  Apparatus; schema record kinds.
- Yields, transmission, neutron yields, β-delayed spectra, ratios to Rutherford, Mott or a reference → Data and
  Results; the schema's 23 quantities.
- Model data: ωγ, E_r, Γ, Γ_c, branching, ANC, scattering length (`datum`) → Data › Model data and a Results table;
  schema `datum` and quantities.
- `gamma_decay` records (partial widths, branchings, mixing, signs) → Nucleus › γ decays; schema record kind.
- Nuisance parameters with priors and grid maxima, profiled norms → every record form (the nuisance widget) and Fit ›
  Parameters; schema value type `nuisance`, the 14 parameter kinds.
- Data covariance files and the prior covariance → Data (a file field), Project › Options; schema `path`, `<config>`.
- The reaction rate → Results › Rate, a job (its call is missing, section 12).
- THM experiments, vertex, line shape, distortion → THM; `<thm>` keys and the session's THM calls.
- The transfer library (bound states, distorted waves, W_l(p)) → Physics; the physics table (pyazr step 4).
- Stopping and straggling (SRIM-2013) → Physics, and ΔE of a layer in the target form; `pyazr.stopping`.
- The evaluation plan: the sub-point grid E_n and the kernel rows K_i(E_n) of a data set → Results › Data set ›
  Grid; `DataSet.grid()`, `kernel()` (designed).
- Diagnostics: 101 codes with severity and meaning → inline in every form, Project › Problems; schema
  `diagnostic_codes`, `Project.validate()`.
- The engine's description of a record in words → beside every form; `Project.describe(name)`.

### 1.3 Forms generated from the schema

The schema (`pyazr.schema()`, 96 kB of JSON) has 10 record kinds with 114 fields (`work/webui/schema_inventory.log`).
Every field has `title`, `type`, `grammar`, `doc`; 75 say `required`, 38 depend on the quantity (`by_quantity`), 21
apply only under a condition (`applies_when`: field present, field equal to a value, the reaction is single, the
record is folded), 20 carry a `dimension` and 14 a unit list, 13 an enum's `values`, 6 a reference (`refers_to`), 8
are `repeatable`, 27 `structured` fields list their `alternatives` as keyword-led sequences of typed `parts`.

Field types and widgets (21 of the 22 value types occur at the top level of a field; `formula` occurs inside parts;
the count of fields in brackets):

- `structured` (27): a choice of alternatives (radio when few, drop-down otherwise); each alternative's parts drawn
  recursively by their own types.
- `enum` (13): drop-down from `values`; for `quantity`, the options the quantity table refuses for this reaction
  are shown greyed with the refusal.
- `integer` (12), `number` (2), `number_list` (1): number boxes with `min` / `positive`.
- `nuisance` (9): value, unit (when the field has a dimension), mode fixed / fit / profiled (where `profiled`),
  prior centre and width (absolute or %), grid max; or the name of a shared record (`or_refers_to`).
- `number_with_unit` (9), `range` (2), `unit` (2): number boxes and a unit drop-down of the field's dimension.
- `text` (9), `name` (2): a text box; names checked by the engine on the edit.
- `reference` (6): drop-down of the records of the kind in `refers_to`, with "+ new" opening that kind's form;
  `@epoch` for targets.
- `path` (4): the server's file browser, rooted at the project directory.
- `spin` (4), `spin_parity` (2): drop-down of the J (and parity) the pairs allow; free text accepted.
- `reaction_expr` (3), `pair` (2): entrance and exit drop-downs of the project's pair labels; "+ term" for a sum,
  "÷" for a ratio.
- `grid` (2): kind (range / list / file / adaptive / none), bounds, unit, frame, step or points, log.
- `columns` (1): the first rows of the file with a role drop-down above each column, the quantity's roles only.
- `nucleus` (1), `composition` (1): text with the engine's parse shown beside it.

What the form does not decide: whether a value is valid. The page sends the edit, the engine validates it as the
project would stand, and a refused edit comes back as diagnostics with `record` and `field` (for example
`E_UNKNOWN_REFERENCE: artemov.norm: ... no normalization record is named 'abc'`), shown under that field; warnings of
`Project.validate()` (50 on 13N, all `W_DATA_ERROR` of negative uncertainties) are listed in Project › Problems and
beside their record. `applies_when` and `by_quantity` only hide fields; the page never refuses on its own.

What cannot be generated: a plot particular to a data kind that is not "value against energy or angle" (section 10),
and the layout of the Nucleus tab (levels × channels is a table, not a record form).

## 2. The spike (measured)

A FastAPI app on today's pyazr (`work/webui/spike/app.py`, 348 lines; `engine_worker.py`, 168 lines) serving one
project: open, describe, validate, evaluate (χ² and one data set's model values), the whole Results payload, a
parameter-vector edit, a project edit (a channel width) with undo and redo, save inside the project directory, a fit
as a job in a worker process with progress, cooperative stop and kill, and process clean-up. It ran over real HTTP
(uvicorn on 127.0.0.1, port 28431, its PID recorded and checked gone; `measure.py`) and with FastAPI's test client
(`test_spike.py`, 5 tests passed). The fit is scipy's L-BFGS-B on `objective_and_grad`: pyazr cannot run the engine's
fitters today. Projects were copies of T1, T2 and 13N under `work/webui/projects/`.

### 2.1 Opening a project

| | T1 | T2 | 13N |
|---|---|---|---|
| `Project.open` through the API (text, history, working file) | 1.2 ms | 1.1 ms | 2.0 ms |
| `Project.validate` / `describe`, direct | 3 ms / 10 ms | 4 ms / 6 ms | 4 ms / 13 ms |
| through the API: worker spawn, import, session build | 0.76 s | 0.89 s | 1.85 s |
| of which the session build | 0.34 s | 0.42 s | 1.34 s |
| session build, first ever (no `output/intEC.dat`) | 4.9 s | 0.81 s | 14.4 s |
| session build, warm (direct) | — | 0.32 s | 0.89 s |
| rebuild after a width edit (API) | 0.07 s | 0.26 s | 0.83 s |
| rebuild after a channel-radius edit (direct) | — | 0.30 s | 14.8 s |

The first build of 13N computes the external-capture integrals and writes them (17.7 MB) to `output/intEC.dat`; a
later build reads them when the signature matches. A radius change invalidates them, and the file keeps one signature
only, so going back to the old radius costs the 13–15 s again (12.8 s measured). A build holds the GIL for its whole
length (12.8 s on 13N: the longest gap of a thread waking every 10 ms, `measure_open.py`); an evaluation releases it
(longest gap 13–16 ms). A server that built sessions in its own process would freeze for that long: hence the
session worker.

### 2.2 Evaluations through the API

`POST /api/evaluate` (the objective and the χ² per data set), median of 20 over HTTP; the engine time is measured
inside the worker; the direct time is the same two calls in one process (`bench_direct.py`):

| | T1 | T2 | 13N |
|---|---|---|---|
| through the API | 9.8 ms | 213 ms | 111 ms |
| engine time inside the worker | 5.3 ms | 203 ms | 98 ms |
| overhead of HTTP, JSON and the pipe | 4.5 ms | 10.4 ms | 13.0 ms |
| the same calls direct, one process | 3.3 ms | 154 ms | 85 ms |
| objective and gradient, direct | 9 ms | 2.77 s | 0.21 s |

A slider on a free parameter can thus redraw at about 100 Hz on T1, 5 Hz on T2 and 9 Hz on 13N; the evaluation, not
the API, is the bound. Two clients asking at once are served one after the other (T2: 0.24 s and 0.40 s, wall
0.40 s), on one history. The machine is shared: the same measurement run three times on T2 gave 192–213 ms per
evaluation; the tables quote the last run, whose logs are kept.

**Against the program.** T2, menu 1 of the program at one thread (`cli_check.sh`): `Total-Chi-Squared: 10905.2`,
`Total-Prior-Covariance-Chi-Squared: 4.90254073`; the API's objective 10910.133842791602 and data χ² 10905.2313 agree.
Per data set the numbers agree except where a data covariance joins two sets: s02 and s03 are 100.8814 and 413.0152
through `segment_chi2` (what the spike serves), 99.1682 and 414.728 in the program's `chiSquared.out` (the sum is the
same, 513.896). The session's own `write_output_files` writes the program's split (`segment_split.py`). The docstring
of `segment_chi2` documents the difference (whitened rows against Σ r_i (M⁻¹ r)_i). Not a wrong number, two
conventions; the UI must show one, the program's (finding 5).

### 2.3 Sizes of the JSON

| | T1 | T2 | 13N |
|---|---|---|---|
| one data set (E, angle, value, error, model) | 2.0–6.9 kB | 0.8–3.1 kB | 1.2–7.1 kB |
| all the Results tab needs (all sets, χ² terms) | 38.7 kB, 442 points | 26.9 kB, 330 points | 48.2 kB, 593 points |
| the same, gzip | 13.6 kB | 11.6 kB | 17.6 kB |
| bytes per point | 88 | 81 | 81 |
| `Session` description (parameters, data sets) | 5.7 kB | 7.5 kB | 9.0 kB |
| `Project.describe()` | 10.4 kB | 12.4 kB | 14.0 kB |

A number costs 23 bytes of JSON against 8 in binary (`json_vs_binary.py`): 10⁴ numbers encode in 8 ms and decode in
5 ms, 10⁵ in 81 ms and 37 ms (2.3 MB), 10⁶ in 0.9 s and 0.4 s (23 MB). **JSON is enough** for every view of a single
evaluation, a dense extrapolation grid (10⁴ points) and a level-matrix scan. Not for raw sampler chains or band draws
(10⁵–10⁶ numbers): the server sends their quantiles, and the chain stays a file in the output directory (decision 10).

### 2.4 Stopping a fit

**The engine's own fit, today.** pyazr cannot start it. The program has no signal handler: SIGINT after 15 s of a
MIGRAD fit of T2 at one thread (`engine_fit_interrupt.sh`) ended the process after about 95 evaluations, χ² from
10910 to 420, and left only `param.par` (the start values) in `output/`: the progress was lost. The program writes a
snapshot (`param.sav`, the output files, `parameters.out`) every 100 evaluations (`AZURECalc::kOutputInterval`); the
same fit run to its end took 38.7 s. What exists in the engine:

- `Config::stopFlag` (per `Config`, so per session): set by nobody since the Qt GUI left. Read by
  `AZURECalc::operator()`, which then returns χ² = 0 to the minimiser during a fit, and by the point loops of
  `WorkingCopy` and `ExternalCaptureSetup`, which skip their remaining points. Not bound. Returning 0 makes MIGRAD
  end on a false minimum; what it writes then was not run here (no way to set the flag).
- `AZURECalc::CountCallsInto(std::atomic<long>*)`: a call counter the program uses for `--max-calls`. Not bound.
- The per-iteration line `Iteration: N Chi-Squared: X` every 10 calls, on standard output with `\r`.

What the binding must offer (finding 2): `fit_with_engine` (decision 8) running without the GIL, with a stop request
callable from another thread that ends the fit cleanly with status `stopped` and the best point (not χ² = 0), and
progress as a polled record (calls, best χ², current χ², EDM, the best x) or a callback every N calls.

**A recipe fit in a job worker (the spike).** The objective wrapper reads a stop event before each call and sends
(call, χ², best χ², best x) after each:

| | T1 | T2 | 13N |
|---|---|---|---|
| cooperative stop: request to "stopped" | 9 ms | 2.10 s | 0.21 s |
| kill (SIGKILL) to "killed" | 0.08 s | 0.05 s | 0.08 s |
| lost by a kill | the running call | the running call | the running call |
| restart: a new job's session build | 0.25 s | 0.45 s | 1.21 s |
| peak memory of a job worker | 182 MB | 222 MB | 194 MB |

The cooperative stop costs at most one objective with gradient. With the best x in every progress message a kill
loses only the call in flight; with the engine's 100-call snapshots it loses up to 99 calls.

### 2.5 Processes, clients, memory

- Memory: the server 77 MB; a session worker 140 MB (T1), 165 MB (T2), 148 MB (13N); a job worker 182–222 MB at
  its peak.
- Two clients on one project share one state: an evaluation of one waits for the other's; an edit of one is the
  other's next view (same history). An evaluation asked while a job runs was answered 409 with the job's id and how to
  stop it (one compute at a time, the spike's budget).
- **The server killed during a fit** (`server_dies.py`): the session worker ended (its pipe closed), the job worker
  did not: it kept computing as an orphan (4.6 s of CPU used when looked at, 2 s after the server died) until
  killed by its PID. The design makes
  workers die with the server (3.4).
- The engine prints its notes (covariance files, groups, THM kinematics, mesh re-planning) on the standard output of
  whichever process runs it: they reach the server's terminal, not the user's page (finding 4).

### 2.6 The thread budget of the spike

Every process ran with `OMP_NUM_THREADS=1`, and no two computed at once: the server refuses engine requests while a
job runs. Not measured for that reason: a job and evaluations side by side, a session at several threads, an emcee
pool.

### 2.7 Not run

The engine's MIGRAD, MINOS and the sampler through Python (not reachable); `Config::stopFlag`'s effect on a running
fit (not bound); a browser (none installed, none to install: no page was written); plotly.js (not on this machine);
the T3 project; Windows and macOS.

### 2.8 The scripts (`work/webui/`, outside the repository)

`env.sh` (environment, one thread); `probe_api.py <T1|T2|13N>` (first inventory); `schema_inventory.py`
(`schema_inventory.log`); `parity_classify.py` (`parity_items.tsv` → `parity_class.tsv`, `parity_counts.md`);
`spike/app.py`, `spike/engine_worker.py` (the spike); `spike/measure.py <project> <port> --fit` (over HTTP; output
`logs_measure_<project>.out` and `spike/measure_<project>.json`); `spike/measure_open.py <project>` (cold, warm and
edited builds; the GIL; `spike/measure_open_<project>.json`); `spike/bench_direct.py` (`bench_direct.log`);
`spike/json_vs_binary.py` (`json_vs_binary.log`); `spike/cli_check.sh`; `spike/segment_split.py`
(`runs/segment_split_T2.log`); `spike/engine_fit_interrupt.sh <seconds>` (`runs/engine_fit_interrupt.log`);
`spike/server_dies.py <port>` (`server_dies.log`); `spike/test_spike.py` (pytest, `logs_pytest.log`).

## 3. The process model

```
 browser tab(s) ──HTTP/SSE── server (uvicorn, one event loop; no engine session)
                               │  Project text + history + revision; working x; job table; budget
                               ├── session worker (one per open project): Session, evaluations, physics
                               └── job workers (one per job): own Session from the same text; progress, stop
```

### 3.1 What lives where

- **Server**: the `Project` (the engine's lossless document, fast: 2–13 ms per open, validate or describe), the
  history of project texts with a revision number, the working parameter vector `x` (by parameter name), the job
  table, the compute budget, the event stream. It never builds a session (a build holds the GIL: 12.8 s on 13N).
- **Working file**: `.<name>.ui.azr` beside the project, rewritten on every edit. The sessions read it (relative paths
  resolve as for the project itself); the user's file is written only by Save. It is also the crash copy.
- **Session worker**: one process per open project, holding the session built from the current text. It builds
  lazily (on the first request that needs it after a text change) and answers evaluations, data-set arrays,
  transformations and per-session physics. Its build is a visible state (`building`, with the build's own progress
  once the engine reports it).
- **Job workers**: one process per job, its own session from the text of the revision the job was started on; a
  progress queue, a stop event, a checkpoint file `output/.ui-jobs/<id>.json` (state, progress, best x) written with
  each progress message.

### 3.2 The compute budget

`pyazr ui --threads N` (default 1 on a shared machine, the core count less one on a laptop: decision 9). Each worker
gets a thread count; the sum never exceeds N. With N = 1 one process computes at a time: an evaluation asked while a
job runs is answered 409 with the job and its stop route, and the Results tab shows the last evaluation and the job's
best point. With N ≥ 2 the session worker keeps one thread for interactive work and jobs share the rest. Two jobs
asked for at once: the second waits in a queue (its state `queued`), shown in the Fit tab with a cancel button.

### 3.3 Reload, two tabs, agents

A page reload loses nothing but a form not yet sent: the page reads `/api/project`, `/api/jobs` and subscribes to
`/api/events` again. Two tabs or two clients (a person and an agent) share one state: every edit carries the revision
it was made on (`If-Match`); an edit on a stale revision is refused with 412 and the current revision, and the page
reloads the affected view. All clients get revision changes and job progress through one server-sent-events stream
(plain `EventSource`, no library; the spike polled).

### 3.4 When the server dies

Workers die with it: on Linux by `prctl(PR_SET_PDEATHSIG)`, elsewhere by a watchdog thread in the worker that exits
when its pipe to the server closes (the spike's session worker did this; its job worker did not). On the next
`pyazr ui project.azr`: a working file newer than the project is offered back ("unsaved edits from 14:02: restore /
discard"); jobs found in `output/.ui-jobs/` are listed as `interrupted` with their last progress and best x, which can
be evaluated or accepted into the project like a finished fit.

## 4. The `/api`

### 4.1 Routes

Request and response bodies are JSON (one route returns text). The calls are those of `pyazr3.md` section 2;
today's equivalent is in brackets where it differs.

Route → pyazr call: request → response.

*Tables of the engine*

- `GET /api/describe` → `pyazr.describe()` and the route table: → versions, build options, routes with their docs,
  the budget.
- `GET /api/kinds`, `GET /api/schema` → `pyazr.kinds()`, `pyazr.schema()`: → the engine's tables.

*The project (the state; every edit takes the revision `rev` it was made on and returns the new one)*

- `POST /api/project/open`, `/new`, `/save` → `pyazr.open`, `Project.new`, `Project.save`: a path inside the root →
  the project's state.
- `GET /api/project` → path, revision, unsaved, can undo / redo.
- `GET /api/project/text`, `PUT /api/project/text` → `Project.text`; a whole-text replace (agents): text, rev → rev.
- `GET /api/project/validate`, `/describe[/{name}]` → `Project.validate`, `Project.describe`: → diagnostics,
  descriptions.
- `GET`/`POST /api/project/records`, `GET`/`PATCH`/`DELETE /api/project/records/{name}`, `POST .../{name}/rename` →
  `Project.records`, `add`, `set`, `remove`, `rename`: kind, name, fields, rev → rev or diagnostics.
- `GET`/`POST /api/project/levels`, `PATCH`/`DELETE .../levels/{label}` → `levels`, `add_level`, `set_level`,
  `remove_level`.
- `PATCH /api/project/levels/{label}/channels/{pair}/{l}/{s}` → `set_channel`: width, fixed, rev → rev.
- `GET /api/project/pairs`, `PATCH .../pairs/{key}` → `pairs`, `set_pair`.
- `GET`/`PATCH /api/project/config` → `config`, `set_config`.
- `GET`/`PATCH /api/project/thm`, `PUT`/`DELETE .../thm/experiments/{name}` → `thm`, `set_thm`,
  `set_thm_experiment`, `remove_thm_experiment`.  `GET /api/thm/keys` → the keys of `<thm>` as record kinds of a
  schema: `pyazr.kinds()["thm"]`, the engine's table (`src/ThmKeys.cpp`) as it is.
- `POST /api/project/undo`, `/redo`; `GET /api/project/history` → the server's history: → rev; (rev, time, what).
- `POST /api/project/apply` → missing (finding 8): the values of x written into the text: x by name, rev → rev.

*The session and the evaluation (derived)*

- `GET /api/session` → `Session.describe()`: → pairs, levels, data sets, parameters, options; the build state.
- `GET /api/parameters` → `Parameters.table()` / `to_records()`: → rows (name, kind, unit, value, free, bounds, owner).
- `GET`/`PATCH /api/x` → the working vector, `Parameters.put`: {name: value} → x.
- `POST /api/evaluate` → `Session.evaluate(x)` (today `objective` and `segment_chi2`): x or none → objective, χ²,
  terms, χ² per data set.
- `GET /api/evaluation/data/{name}` → `Evaluation.data[i]` (today `calculate_rwa` and the `data_*` bindings): →
  energies, angles, values, errors, model, residuals, extra columns.
- `GET /api/evaluation` → the whole `Evaluation`: → every data set (the Results payload).
- `POST /api/evaluate/grid` → `Session.evaluate` on energies and angles given (`pyazr3.md` 2.3): record, grid → arrays.
- `POST /api/transform` → `Session.transform(x)`: → physical values per level and channel.
- `GET /api/system` → `Session.system`: → pairs, levels, channels, Wigner limits.
- `GET /api/data/{name}/grid`, `/kernel` → `DataSet.grid()`, `kernel()`: → E_n; the rows K_i(E_n).
- `POST /api/outputs` → `Session.write_outputs(x)`: → the files written.

*Physics, libraries, files*

- `GET /api/physics`, `POST /api/physics/{fn}` → `kinds()["physics"]`, `pyazr.physics.<fn>`: arguments → arrays and
  a status per element.  **Served**: the body `{"arguments": {name: number or list}}` (broadcast), always with
  `errors="mask"`; the answer `returns` (fields as lists), `units`, `status` per element and `status_words`, `shape`,
  `axis`; refusals (unknown function 404; unknown or missing argument, wrong type, shapes that do not broadcast,
  more than 20000 elements 422) carry what / allowed / fix built from the row.  An argument may also be an axis
  `{"from", "to", "points", "spacing"}` (linear, log, integer), made by the server with numpy, and a choice a list of
  names; the answer adds `arguments` (each as used, the axis expanded) and `python` (the text of the pyazr call that
  makes the same arrays, section 6.6).  It runs in the server's process
  (`pyazr/ui/physics.py`): the functions need no session, and the C++ loop releases the GIL; the bound keeps one
  request near 1 s (coulsm, 45 µs per element; `shift_derivative` about five times that).
- `POST /api/session/channel_functions`, `/level_matrix` → `Session.channel_functions`, `level_matrix`.
  `channel_functions` is **served** in the session worker (adapter `Engine.channel_functions`, the budget as the
  other session routes), with the same bound of 20000 energies; `level_matrix` stays 501 (no engine call).
- `GET /api/session/functions`, `POST /api/session/functions/{fn}` → the session's functions as rows of the shape of
  `/api/physics`: pyazr's rows (`Session.describe()["functions"]`) with the page's own part (`PAGE` of
  `pyazr/ui/physics.py` and `thm.py`), and a call of one (`Engine.session_function`, by the row's name) answered as
  `/api/physics/{fn}` answers (status -2 where the channel is closed, the Python text opening a session on the
  project file).
  **Served** (the Physics tab; the THM diagnostics below are rows of it too, drawn by the THM tab).
- `POST /api/session/thm_experiments`, `/thm_background`, `/thm_lineshape`, `/thm_vertex`, `/thm_distortion` →
  `Session.thm_experiments`, `thm_background`, `thm_lineshape`, `thm_vertex`, `thm_distortion`: the session's
  answer as plain data, at the working x.  **Served** (adapter `Engine.thm_*`); the same calls as rows of the shape
  of `/api/physics` through `/api/session/functions/{fn}` (sections 11.8, 11.9).
- `GET /api/stopping`, `POST /api/stopping/{fn}` → `kinds()["stopping"]`, `pyazr.stopping.evaluate`: the rows of
  the binding's stopping table (the shape of the physics rows, a text argument `formula` first, SRIM's `citation`)
  and a call of one, answered as `/api/physics/{fn}` answers, with the citation. **Served** (the same generic
  handler, `pyazr/ui/physics.py` `TABLES`).
- `GET /api/sources/exfor`, `/exfor/{id}`, `/ensdf/{nuclide}`, `/masses/{nuclide}` → `pyazr.sources.*` (network
  for EXFOR and ENSDF).
- `GET /api/files[/{path}]` → the server's listing inside the root: → names; the first rows of a data file.

*Jobs and events*

- `GET`/`POST /api/jobs`, `GET /api/jobs/{id}`, `POST .../stop`, `POST .../accept` → the job table (4.4): kind,
  settings, x0 → id; → state, progress, result.
- `GET /api/events` → server-sent events: revision changes, job progress, build state.

### 4.2 The rule that keeps it one-to-one

The routes are rows of one table, `pyazr/ui/routes.py`: (method, path, the pyazr callable by its dotted name, how the
arguments map, which object answers `describe()` for the response). The handlers are made from the rows by four
generic functions (project edit, session read, job, static), as `Project.add_<kind>` is made from the schema; there
is no handler written per route. A consistency test fails when a row names a callable that does not exist, and when a
public method of `Project`, `Session`, `Parameters`, `Evaluation` or a row of the physics table has no route and is
not on the table's list of calls not served, each with its reason (`Project.session`: the server owns sessions;
`Session.close`: the worker's; `Objective`: jobs only).

### 4.3 Errors

pyazr's classes (`pyazr3.md` 2.5) map to statuses; the body keeps the three parts of the message:

| class | status | body |
|---|---|---|
| `FormatError` | 422 | `error`, `message`, `diagnostics` (code, severity, record, field, line, found, text) |
| `ParameterError`, `DomainError` | 422 | `error`, `what`, `allowed`, `fix` |
| `EngineError` | 422 (set-up or evaluation refused); 500 (a worker crashed) | `error`, `message`, the record or level |
| `SessionClosed` | 409 | the call; that the session is being rebuilt |
| a job running, budget full | 409 | the job's id, its stop route |
| stale revision | 412 | the current revision |
| an unknown record, level, job | 404 | the names that exist (as pyazr's "did you mean") |
| a path outside the root | 403 | the root |

The page shows `what` as the message (`message` where `what` is null), `allowed` and `fix` beneath it.  The three
parts are `e.what`, `e.allowed`, `e.fix` of pyazr's error (carried from a worker with the class and the message);
they are null where the message is the engine's own text, and the server takes no message apart.

### 4.4 Jobs

A job is a row of the job table (`pyazr/ui/jobs.py`): its kind, the pyazr call it runs, its settings (from
`kinds()` where they exist), and what it hands back. Kinds: `fit` (engine: MIGRAD, MIGRAD with gradient, LM, GSL-LM
through `fit_with_engine`; recipes: scipy L-BFGS-B, scipy least squares, iminuit MIGRAD + HESSE), `minos` (the
engine's through `fit_with_engine(minos=True)`, or iminuit's), `sample` (the emcee recipe), `bands`
(`pyazr.bands.uncertainty_bands`, `extrapolation_bands`), `rate` (missing, finding 6), `coulomb_table` (decision 7 of
pyazr). State: `queued`, `building`, `running`, `stopping`, `stopped`, `converged`, `ended` (with the minimiser's
message), `failed` (with the error), `killed`, `interrupted`. Progress: calls, χ², best χ², elapsed time, best x,
and per kind more (EDM; acceptance fraction and τ of the sampler). `accept` writes the result's x into the project
as one history entry (finding 8).

### 4.5 For an AI agent

The page is not needed: `GET /api/describe` lists every route with the docstring of its pyazr call and the shape of
its request; `/openapi.json` is FastAPI's (the spike's lists its routes); every response is JSON; every refusal
carries the engine's code. `pyazr ui project.azr --no-browser` prints the URL with its token on one line and writes
it to `output/.ui-url` (mode 0600) for an agent on the same account. An agent that does not need to share the user's
live session uses pyazr directly; the API is for working beside the user, on what the user sees.

## 5. The state rule

- **The project text is the state.** Every change of the model is an edit through `Project` (records, levels,
  channels, pairs, config, THM), validated by the engine on the project as it would stand; a refused edit leaves the
  project unchanged and comes back as diagnostics. JavaScript never validates.
- **Undo and redo** step through the history of project texts (the spike's history: undo restored the objective
  exactly on T1, T2, 13N, and the text byte for byte by construction: the history keeps the texts). A history
  entry is one user action (a slider's release, not its ticks).
- **Derived state**: the session (rebuilt when the text differs from the text it was built from), the evaluation
  (recomputed when the session or `x` changes), the working vector `x` (reset to `x0` on a rebuild, except the values
  of parameters of the same name, which the server carries over).
- **Session options** (Brune, Park, transform, Coulomb table, ...) are keys of `<config>` (decision 7, built
  2026-11-13: `project_format3.md` §3.12): setting one is a project edit through the config route (a history entry,
  a new session); the server opens its sessions without keywords, so the project's key or the default decides.

**Which edits rebuild, and what that costs** (2.1, 2.2):

| edit | needs | T2 | 13N |
|---|---|---|---|
| value of a free parameter (energy, width, norm, shift, density) | a new x only | 0.21 s | 0.11 s |
| a fixed value, a fix/free flag, a record field, a level or channel added | a new session | 0.26–0.42 s | 0.8–1.3 s |
| a radius, a pair, what the external-capture integrals depend on | a session and the integrals | 0.30 s | 14.8 s |

**A slider on a level width.** On a free width the slider moves `x` and evaluates on each tick, the newest request
replacing any not yet started (5 Hz on T2, 9 Hz on 13N); on release the value is written into the project text (one
history entry; finding 8). The slider shows the physical width (`Session.transform(x)`) but moves the fit basis;
setting an exact physical value needs the inverse transformation, which pyazr 3.0 does not have (finding 7): until
then a typed value is a project edit (a rebuild). On a fixed width or a radius the slider is a project edit, applied
on release only, with the rebuild's progress shown (15 s for a 13N radius).

## 6. The tabs

The page has a header (project name, unsaved mark, revision, undo/redo, Save, the budget and the running job), seven
tabs and a status line (the last message of the engine or a job). Keys: `1`–`9` and `0` switch to the first ten
tabs shown, `[` and `]` to the previous and the next (built: 11.16), `Ctrl+Z`/`Ctrl+Shift+Z`
undo and redo, `Ctrl+S` save, `/` searches the current table, `Esc` closes a form.

### 6.1 Project

Shows the file, the `<config>` entries, session options, checks and directories, the problems of `validate()`, the
text itself (read-only, with a "replace text" for experts and agents), the history. Routes: `/api/project*`,
`/api/kinds` (options).

```
┌ Project ─ T2.azr ───────────────────────────────── rev 14 · unsaved ─ [Undo] [Redo] [Save] [Save as] ┐
│ File      ~/T2/T2.azr                    Working copy .T2.ui.azr (14:02)                             │
│ Formalism (•) A-matrix ( ) R-matrix        Output dir output/   Checks dir checks/                   │
│ Options   [x] Brune  [ ] Park  [x] transform  [x] Coulomb table (v2, hash 3f…)  [ ] GSL Coulomb      │
│           [ ] relativistic   prior covariance [prior_cov.f3.txt     ][Browse]                        │
│           hard-sphere in every Jπ [x]   max l [2]   max γ multipolarity [2]                          │
│ Checks    compound nucleus [none ▾] boundary [none ▾] data [none ▾] ... (8)                          │
│ ─ Problems (validate) ───────────────────────────────────────────────────────────────────────────── │
│ ⚠ W_DATA_ERROR  baumann_1618.file  line 2: the uncertainty is negative ...          [go to record]   │
│ ⚠ W_BROAD_LEVEL level 5 ...                                                          [go to level]   │
│ ─ History ───────────────────────────────────────────────────────────────────────────────────────── │
│ 14  14:02  set_channel 1- 12.4258 (1, 0, 1/2) width 547.9 → 602.7 eV                              │
│ 13  13:58  set s02_diff090.norm "1 fit prior 1 8%"                                                │
│ 12  13:41  add_dataset s13_thm_oes                                                                │
│ ─ Text ─────────────────────────────────────────────────────────── [Show] [Replace text…] ──────── │
│ ─ Engine ── AZURE2 <version> (<commit>) · build: MCMC off, NLopt on · threads 1 of 1 · SRIM-2013 ── │
└──────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### 6.2 Nucleus

Pairs, levels with their channels, γ decays (`gamma_decay` records), hybrid potentials (`potential` records). Routes:
`/api/project/pairs`, `/levels`, `/records?kind=gamma_decay|potential`, `/api/system`, `/api/sources/masses`,
`/api/sources/ensdf`.

```
┌ Nucleus ───────────────────────────────────────────────────────────────────────────────────────────┐
│ Pairs                                                                     [+ pair]¹ [AME2020 check] │
│ key label          J1π  J2π   M1 (u)    M2 (u)     S (MeV)  a (fm)  E_x2   check                   │
│ 1   p+15N          1/2+ 1/2-  1.007825  15.000109  12.1274  4.85    0      S ok (0.2 keV)          │
│ 2   a+12C          0+   0+    4.002603  12.000000   7.1619  5.43    0      S ok                    │
│ 3   a+12C[4.43891] 0+   2+    4.002603  12.000000   7.1619  5.43    4.4389 S ok                    │
│ Levels  filter [      ] nucleus [16O ▾]                       [+ level] [IAEA…] [Fix all] [LaTeX]  │
│ on fix Jπ   E_x (MeV)   channels                                  Γ / ANC (row unit)   Wigner      │
│ [x] [ ] 1-  12.425800   (1,0,1/2) (1,2,1/2) (2,1,0) (3,1,2)       547.9 eV …          [Γ_W]        │
│ [x] [x] 1-  13.105000   (1,0,1/2) (1,2,1/2) (2,1,0) (3,1,2)       40000 eV …                       │
│ [x] [x] 3-  13.073400   ...                                                                         │
│ ─ Level 1- 12.425800 ──────────────────────────────────────────────────────────────────────────── │
│ ch  pair   l  s    kind        width        unit       fixed  slider                      Γ_W      │
│ 1   p+15N  0  1/2  particle    547.90014    eV         [ ]    ├────●──────┤ (free: x)      2.1 MeV  │
│ 2   p+15N  2  1/2  particle    0            eV         [x]                                  0.3 MeV │
│ 3   a+12C  1  0    particle    75000        eV         [x]    (fixed: project edit)        …       │
│ Engine: "level 1- 12.4258 MeV: Γ_p = 547.9 eV, ..." (describe)                                    │
│ γ decays [+]   Potentials [+]   (record forms, from the schema)                                    │
│ ¹ add_pair, remove_pair (its channels too as a second action), set_pair with B (section 11.10)     │
└────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### 6.3 Data

Every record of the sections `data`, `targets`, `extrapolations`, `modelData`, grouped by section as the schema lists
them; one generated form per record with the engine's description beside it and a preview plot of the data with the
current model. Routes: `/api/project/records*`, `/api/project/describe/{name}`, `/api/schema`, `/api/files`,
`/api/sources/exfor*`, `/api/evaluation/data/{name}`.

```
┌ Data ─────────────────────────────────────────────────────────────────────────────────────────────┐
│ [Data sets 13] [Norms & scales 2] [Apparatus 9] [Extrapolations 3] [Model data 1]   [+ new ▾][EXFOR]│
│ on name           reaction      quantity        unit  target        norm   χ²/N   rows             │
│ [x] s01_int_thin  15N(p,a)12C   sigma           b     tgt_s01…      fit    17.4   40/40            │
│ [x] s02_diff090   15N(p,a)12C   dsigma_domega   b/sr  —             D      3.3    30/30            │
│ [x] s05_implanted 15N(p,a)12C   yield           1/p   tgt_s05…      fit    116.6  26/26  ⚠1        │
│ ─ s02_diff090 (dataset) ──────────────────────────────── [Duplicate] [Rename] [Delete] ──────────── │
│ reaction   [p+15N ▾] → [a+12C ▾]  [+ term] [÷]       quantity [dsigma_domega ▾]  unit [b/sr ▾]     │
│ file       [data/s02_diff090.dat ][Browse]   columns  E_lab[MeV]▾ theta_lab[deg]▾ value▾ error▾    │
│            0.25  90  0.0213  0.0011     (first rows of the file)                                   │
│ windows    E [0.249999] – [1.350001] [MeV ▾]    θ [90] – [90] [deg ▾]                              │
│ norm       ( ) value [1] mode [fit ▾] prior [1] ± [8 %]   (•) shared: [D ▾]                        │
│ shift      (•) shared: [G ▾]                                                                        │
│ covariance [cov_diff.txt][Browse]   beam [beam_targetInt2 ▾][+]   detector [det_targetInt2 ▾][+]   │
│ quadrature points [20]                                                                              │
│ ⚠ E_...: messages of the engine under the field they concern                                       │
│ Engine: "differential cross section dσ/dΩ(lab) of 15N(p,α)12C at 90°, ... normalization D (fit,    │
│          prior 1 ± 8 %), 30 rows read, 0 excluded"                                                  │
│ ┌ preview: data (points) and model (line) against E_lab ───────────────────────────────────────┐   │
│ └───────────────────────────────────────────────────────────────────────────────────────────────┘   │
└────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### 6.4 Fit

The parameter table (every parameter, filtered by kind; free ones editable as `x`, fixed ones as project edits), the
job form, the running job, and the result with errors. Routes: `/api/parameters`, `/api/x`, `/api/jobs*`,
`/api/events`, `/api/project/apply`.

```
┌ Fit ──────────────────────────────────────────────────────────────────────────────────────────────┐
│ Parameters  kind [all ▾] [free only ✓]  /filter                         62 parameters, 19 free    │
│ name                       kind      unit     value        bounds         error    prior          │
│ energy_1                   energy    MeV      12.4258      —              ±0.0011  —              │
│ width_1_1                  width     MeV^1/2  0.4710       —              ±0.012   —              │
│ s01_int_thin.norm          norm      1        1.000        —              ±0.021   1 ± 5 %        │
│ tgt_s05_implanted.areal…   density   at/cm²   1.2e18       [0, 5e18]      ±4e16    —              │
│ ─ Job ──────────────────────────────────────────────────────────────────────────────────────────  │
│ kind [fit ▾]  method [engine MIGRAD (gradient) ▾ | iminuit MIGRAD+HESSE | scipy L-BFGS-B | ...]     │
│ start [current x ▾]  max calls [ ]  threads [1 of 1]   [ ] MINOS after   [ ] bands after   [Run]   │
│ ─ Running: fit 36fea7be  engine MIGRAD  ──────────────────────────── [Stop] [Kill] ───────────── │
│ calls 412   χ² 689.26  best 689.26   EDM 3e-4   12 s   ▁▂▃▅▇ (χ² against call, log)               │
│ log: Iteration 410 Chi-Squared 689.26 · Note: the energy shifts moved 130 folded points ...        │
│ ─ Result ───────────────────────────────────────────────────────────────────────────────────────── │
│ converged · χ² 689.25 (N 330) · covariance from HESSE · 25 iterations, 31 calls, 23 s              │
│ [Accept into project]  [Evaluate]  [Save .sav]  [Covariance ▸]  [MINOS ▸]  [Compare with start]    │
│ ─ Sampler ──────────────────────────────────────────────────────────────────────────────────────── │
│ emcee  walkers [50] steps [10000] spread [5 %] thin [10]   progress 3400/10000  acc 0.31  τ ≈ 120  │
│ [Corner ▸] [Traces ▸] [Statistics: mean, sd, 2.5 / 50 / 97.5 %] [Bands from draws ▸]              │
└────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

**Engine fitters or recipes.** Both, as rows of one job table (decision 2). The engine's MIGRAD and MINOS are the
reference behind the program, so the default fit is the engine's (`fit_with_engine`, decision 8 of pyazr); the
recipes (iminuit, scipy, emcee) are the objective-based alternatives and run today. MINOS: the engine's when the fit
was the engine's, iminuit's after an iminuit fit; both give asymmetric intervals in the same table.

**The sampler: emcee, in a job worker.** `numcmc` goes (pyazr decision 3 (b)); pyazr's sampler is the emcee recipe on
`Objective.log_prob` (`pyazr3.md` 2.6); its progress is natural (one `sample` step at a time), it stops between
steps, and its pool is one session per worker process, counted in the budget. The chain is written to
`output/chain_<id>.npy` with a JSON summary; the page receives quantiles, not draws (2.3).

### 6.5 Results

Plots of every evaluated data set and extrapolation, the χ² table, model data, outputs, the reaction rate. Routes:
`/api/evaluation*`, `/api/evaluate/grid`, `/api/transform`, `/api/jobs` (bands, rate), `/api/outputs`.

```
┌ Results ─ at: [current x ▾ | job 36fea7be best | start x0]  χ² 689.25 + priors 4.90 = 694.15 ─────┐
│ data sets [x] s01 [x] s02 [ ] s03 ...  extrapolations [ ] sfactor_low   [all] [none]               │
│ x [E_cm ▾ | E_x | E_lab | θ_cm]  [log x]   y [as measured ▾ | σ | S-factor | ratio]  [log y]       │
│ show [x] model [x] residuals [x] levels [ ] 1σ band (job ▾) [ ] components (E1, E2)                │
│ ┌──────────────────────────────────────────────────────────────────────────────────────────────┐   │
│ │  data with error bars (one colour per set) · model line · level markers at E_x               │   │
│ │                                                                                              │   │
│ └──────────────────────────────────────────────────────────────────────────────────────────────┘   │
│ ┌ residuals (data − model)/error ──────────────────────────────────────────────────────────────┐   │
│ └──────────────────────────────────────────────────────────────────────────────────────────────┘   │
│ [Export data (CSV)] [SVG] [PNG] [PDF (matplotlib)]                                                 │
│ ─ χ² per data set (the program's split) ───────────────────────────────────────────────────────── │
│ s01_int_thin 697.5/40   s02_diff090 99.2/30   s03_diff150 414.7/30   ...   norm terms ...          │
│ ─ Angular distributions ─ E [1.02 MeV ▾]  sets at that energy, model against θ ──────────────────  │
│ ─ Model data ─ ωγ(5) 0.4800 eV (data 0.25 ± 0.02 eV), ...                                          │
│ ─ Reaction rate ─ entrance [p+15N ▾] exit [a+12C ▾] T [0.01]–[10] GK [Run rate job] ¹             │
│ ─ Outputs ─ [Write the program's output files] → output/ (AZUREOut_*, chiSquared.out, ...)         │
│ ¹ waits for a pyazr call (finding 6)                                                               │
└────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

Plots: data with model (points with error bars, line); residuals; σ ↔ S-factor (the engine's conversion per point)
↔ the measured quantity; angular distributions at a chosen energy; analyzing power and phase shifts on linear axes;
Legendre coefficients per order; bands of a sampler or of the covariance (`pyazr.bands`) as shaded regions;
extrapolations on dense grids (`/api/evaluate/grid`); the reaction rate against T₉ with its ratio to a reference.

### 6.6 Physics

One generated form and one plot per row of the physics table (`pyazr3.md` 2.7), plus the per-session functions.
Routes: `/api/physics*`, `/api/session/channel_functions`, `/level_matrix`, `/api/stopping*`, `/api/data/{name}/grid`
and `/kernel`.

```
┌ Physics ──────────────────────────────────────────────────────────────────────────────────────────┐
│ search [        ]          │ coulomb_waves(l, eta, rho, method)                                     │
│ kinematics                 │ F, G, F′, G′ and σ_l at (l, η, ρ), by every method of pycoulomb        │
│   sommerfeld_eta …         │ plot against [rho ▾]                                                   │
│ Coulomb functions          │ rho   from [0.01] to [20]   points [400] (≤ 5000 for 4 lines) [log ▾]  │
│ ▸ coulomb_waves            │ l     [0 1 2 3]   (several values: one line each)                      │
│ channel functions …        │ eta   [0.5]                                                            │
│ stopping …                 │ method [coulsm ▾]                                    [Plot] [Defaults] │
│ this project               │ ┌ F, Fp, G, Gp against rho, one line per l ──────────────────────────┐ │
│   channel_functions        │ └────────────────────────────────────────────────────────────────────┘ │
│                            │ 1580 of 1600 elements computed; rho outside its domain (finite > 0): 20│
│                            │ fields [x] F [x] Fp [x] G [x] Gp [ ] sigma   [log x] [log y]            │
│                            │ [Export (CSV)] [SVG] [PNG] [Copy as Python]                            │
│                            │ import numpy as np … pyazr.physics.coulomb_waves(l=…, rho=np.geomspace…)│
└────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

**The list.** The rows of `GET /api/physics` (the binding's physics table), of the stopping table (`GET
/api/stopping`, the second part of this step) and of the session's functions (`GET /api/session/functions`), grouped
by the row's `group` in the table's order, the session's under "this project"; a search box filters by name,
summary and group. A row of any table appears in the list and gets its form with no line of JavaScript.

**What the rows must say** (the smallest addition to the table, made in `api/src/pybind_physics.cpp`, binding code
only, read by `kinds()` unchanged): per row a `group` (the comment lines of the table become a field); per argument
an `example` (a typical value: the form's start value where the argument has no `default`) and, for an argument a
plot may run over, a `sweep` (`lower`, `upper`, `spacing` linear or log: the range of a first plot). The row's `axis`
has a sweep; every argument without a default has an example or is the axis (`physics_table_test.py` checks both).
A list of defaults in JavaScript was the alternative and is not taken: the row would not be complete for an agent.

**How a row becomes a form.** A form of the existing generator (`forms.js` `describeForm`, drawn by `formview.js`):
the tab writes the row as a record kind of a schema of its own, one field per argument, with these value types (one
table, `ARGUMENT_KINDS` in `physicsform.js`, maps an argument kind to its field):

- the axis (the row's `axis`; the selector offers every argument with a `sweep`): a `range` field (from, to, the
  unit shown), `points` (an `integer`; the bound shown is the route's 20000 elements divided by the number of
  lines) and `spacing` (`enum` linear / log, the sweep's by default); an integer axis takes every whole number from
  `from` to `to` and has no points;
- a real argument: a `number_list` (one value, or several separated by blanks: one line each);
- an integer: an `integer` box (also several values: `0 1 2 3` gives four lines);
- a choice (`method`): an `enum` of its choices and "all" (one line per choice);
- a flag: an `enum` of 0, 1 and "both";
- a pair of the project (the session's functions): an `enum` of the project's pairs (`GET /api/project/pairs`).

Start values: the argument's `default`, else its `example`; the axis from its `sweep` with 200 points. Every row
plots without typing. Several values of several arguments give the product of their lines (at most 24); each list
is a broadcast dimension of its own, in the order of the row, the axis the last.

**The request and the answer.** `POST /api/physics/{fn}` with `{"arguments": {...}}`: a list argument as a nested list
of the shape that broadcasts (`l: [[0], [1], [2], [3]]`), a choice list by names, the axis as `{"from", "to",
"points", "spacing"}`, which the server turns into the array with numpy (`linspace`, `geomspace`, `arange` for an
integer axis): the page makes no number. The answer adds `arguments` (each argument as used, the axis expanded) and
`python`, the text of the call that makes the same arrays.

**The plot.** One plot per unit of the returned fields shown (checkboxes, all shown at first); one trace per field
and line, named by the line's values (`F, l = 2`); an element whose status is not 0 is a gap (null in the trace,
never a zero) and is counted under the plot. The status line: "N of M elements computed", then for each status that
occurs its words (`status_words` of the answer) and its count. Log x and log y on request (log x at first where the
spacing is log; what a log axis drops is counted). Axis titles: the argument and field names with the row's units.

**Export.** CSV of exactly the arrays returned: one line per element, the arguments at that element (a choice by its
name), every returned field to 17 significant digits ("nan" where not computed) and the status; units in the header.
SVG and PNG of the figure through plotly's export (as the other tabs).

**Copy as Python.** The answer's `python` (made by the server from the request, `pyazr/ui/physics.py`): `import numpy
as np`, `import pyazr`, the axis as `np.linspace` / `np.geomspace` / `np.arange` of the numbers sent, a list as
`np.array([...]).reshape(-1, 1, ...)`, and `result, status = pyazr.physics.<fn>(..., errors="mask")`; shown under the
plot and copied by the button. The test runs the text and compares the arrays with the route's bit for bit.

**The session's functions.** `channel_functions` (pair, `e_cm` the axis, l, radius, J) is pyazr's row
(`Session.describe()["functions"]`), served by `POST /api/session/functions/channel_functions` in the session worker
(section 11.9). When the project does not build, the form is drawn and
the plot area shows the engine's refusal (what / allowed / fix) with "the functions above need no project"; the
free functions keep working. Its Python text opens a session on the project file.

**Stopping (the second part).** The stopping library's functions become rows of the same shape (arguments with unit,
kind, domain, example and sweep; returns with units; `group` "stopping"; and the `citation` of SRIM, which the PI
requires beside every use of the SRIM data), made by the binding with the same row types as the physics table
(shared, not copied) and read as `kinds()["stopping"]`; a text argument (the material's formula, its phase) is a new
argument kind, one value per call. `POST /api/stopping/{fn}` is the same generic handler as the physics route, on the
stopping rows. Every stopping plot carries the citation under it, and the CSV and the Python text carry it too.

**State (section 5).** The tab holds no model state. The form's texts per function, the axis, the fields shown and
the log choices are page state, kept in the module while the page lives and lost on reload, as the Results tab keeps
its view; hash parameters set them, so a view can be linked: `#/physics/coulomb_waves?axis=rho&from=-2&to=20&v.l=0 1`
(`v.<argument>` an argument's text; `axis`, `from`, `to`, `points`, `spacing`, `fields`, `logx`, `logy`, `plot=table`).

### 6.7 THM (shown when the project has THM content)

Model options, experiments, channels (binding energies, amplitude flags), diagnostics. Routes: `/api/project/thm*`,
`/api/project/pairs`, `/api/session` (THM calls), `/api/physics` (vertex, distortion).

```
┌ THM ──────────────────────────────────────────────────────────────────────────────────────────────┐
│ Model   entrance waves [incoherent ▾]  vertex [perlevel ▾]  kinematics [lacognata ▾]               │
│         external Coulomb [x]   spectator energy [  ] MeV   per pair [table]   weights [table]       │
│ Experiments                                                                    [+ experiment]      │
│ name   data sets                 background  reaction                     Ebeam   acceptance       │
│ LC09   s11_thm_hoes_a, s12_…     none        15N + d (d = x + n)          60 MeV  Hulthén          │
│ ─ LC09 ─ beam [d] target [15N] spectator [n]  B(x+s) 2.22457 MeV  E_qf 1.55188 MeV               │
│ distortion [none ▾]  bound state [Whittaker ▾]  optical a+A [global ▾] s+F [plane ▾]              │
│ coherent background  Jπ  exit  channels  form  Re c0  Im c0  c1                                    │
│ Channels  pair  B (MeV)  amplitude?   (as the session reads them; edited in the Nucleus tab)       │
│ ─ Diagnostics ─ data set [s11 ▾]  Jπ [1- ▾] [Compute]                                              │
│ ┌ |M_l|² (vertex) · HOES vs on-shell · line shape |N_C|² and ζ(E) · w(E) · R(E) ────────────────┐ │
│ └──────────────────────────────────────────────────────────────────────────────────────────────┘   │
└────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

## 7. The page

**Files** (`pyazr/ui/static/`, served by the server; no build step, no package manager): `index.html` (60 lines),
`app.js` (tabs, routes `#/data/s02_diff090`, keys; 300), `api.js` (fetch, revision, errors, events; 200), `forms.js`
(schema → form, one widget per value type; 700), `plots.js` (plot kinds, CSV export, palette; 400), `tables.js`
(sortable, filterable, keyboard; 250), `tabs/project.js`, `nucleus.js`, `data.js`, `fit.js`, `results.js`,
`physics.js`, `thm.js` (200–400 each), `style.css` (400), `vendor/`. About 4500 lines; JavaScript modules loaded by the
browser as they are.

**plotly.js, vendored** (described, not fetched: decision 11). The partial bundle `plotly-cartesian.min.js` (scatter
with error bars and filled bands, bar, heatmap, histogram, contour: all the plots above), version fixed at the time
it is vendored (the 3.x line, the current one; the 2.x line if a plot needs a trace 3.x removed), licence MIT, at
`pyazr/ui/static/vendor/plotly-cartesian-<version>.min.js` with its `LICENSE` beside it and a line in the wheel's
licence list. Its size is of the order of 1–1.5 MB (the full bundle about 4–5 MB); to be measured when vendored. The
page works without it: plots are replaced by their data table and CSV export, so the tests need no plotly.

**The look.** A scientific tool: dense tables in a monospace figure font, little chrome, one accent colour; every
action reachable by keyboard and listed in a `?` overlay; light and dark themes following the system, switchable;
numbers right-aligned with units in the column head; the engine's own words for every description and message.

**The original's face.** The header and the browser tab show the original program's icon
(`pyazr/ui/static/azure2_icon.png`, byte-identical to `packaging/icons/azure2_icon.png`: consistency check 57).
The colour tokens of `style.css` are built on the icon's measured colours: its blue `#004293` is the accent, its red
`#990000` the error colour, its green `#04800d` (darkened to `#04720c` for 4.5:1 on every panel) the success colour;
the dark theme uses lighter tints of the same hues (`#79aef0`, `#ff8f8a`, `#6fd27a`). A state is never said by its
colour alone: the error, warning and success boxes and the status line carry ✕, ⚠ or ✓ beside it (CSS `::before`).
`tests/pyazr/ui/static_test.py` computes the WCAG contrast of every text token on every background it is drawn on,
in both themes, from the CSS (at least 4.5:1). **About** (header) gives the program's and pyazr's versions, the
reference to cite (as `docs/source/index.rst` gives it; the test compares the two) and the documentation's address.

**Plots and accessibility.** A colour-blind-safe palette (Okabe–Ito for up to 8 data sets, then markers vary as
well as colours); every plot has "Export data (CSV)" with the exact numbers drawn and their units, SVG and PNG
export from plotly, PDF through `pyazr.plot` (matplotlib, server side) for a paper figure; every plot's title and
axes name the quantity, frame and unit; a plot is never the only place a number is.

## 8. Security

- Binds to 127.0.0.1 only; `--host` other than a loopback address is refused (remote use goes through an SSH tunnel:
  decision 4).
- **A token, as Jupyter does**: a random token in the URL the command prints, exchanged for a cookie on the first
  load; every `/api` request needs it. Reason: on a shared machine (this one, a cluster's login node) every local user
  can reach 127.0.0.1 ports; without a token any of them could edit and save a user's project. `--token` sets a
  fixed one for agents and tests. The cookie is `SameSite=Lax` and named for the port (`azure2_ui_token_<port>`);
  a browser that asks for the page without the token gets a page with a field for it, a program the JSON refusal
  (as built: 11.1, item 2).
- The `Host` header must be `127.0.0.1` or `localhost` (DNS rebinding); no CORS headers; JSON bodies only.
- File access is limited to a root: the project's directory by default, `--root DIR` to open several; paths are
  resolved (symlinks included) and refused outside the root (403; the spike refused `../../escape.azr`).
- No shell, no code evaluated from requests: formulas go to the engine's expression parser; network only for the
  EXFOR and ENSDF imports, on the user's request.

## 9. Tests

- **pytest on `/api`** with T1, T2, T3 (`tests/pyazr/ui/`), through FastAPI's test client, workers at one thread:
  the objective and every data set's model values equal a direct session's bit for bit (the spike's test); the totals
  equal the program's `chiSquared.out` on the same text; an edit, undo and redo restore the text byte for byte and the
  χ² exactly; every class of 4.3 gives its status and its fields (a table of wrong calls, those of `pyazr3.md` 1.5
  among them); a stale revision gives 412; a path outside the root 403; a job reports increasing call numbers, stops
  within one objective call, a kill leaves no process (PIDs checked), `interrupted` jobs are found after a restart;
  the JSON of each route validates against its OpenAPI schema.
- **The route-table consistency test** (4.2), and a check that `forms.js` has a widget for every value type of the
  schema (by name, in its widget table) and `plots.js` a rule for every quantity's axes.
- **Browser smoke tests**: Playwright for Python (needs the `playwright` package and a Chromium of about 150 MB
  downloaded by `playwright install`; not installed here, and not to be). Open T1, visit each tab, assert no console
  error, one plot drawn where plotly is present, a screenshot per tab (this replaces Qt's `AZURE2_SCREENSHOT_DIR`).
  Run in a workflow of their own, not in `ctest`.
- **The parity list as a checklist**: `qt_parity.md` gets a test name or a reason per item: E and D items a test of
  `/api` (about 130), P items a smoke-test step (29), X items their reason (20), M items wait for pyazr (20).

## 10. Extensibility (the PI's measure)

What must know of a new thing under this design:

- **A data kind** (a quantity row of the engine): none. Forms, the Data table, χ² and the default plot come from the
  quantity row: `columns.energy` and `angle_axis` give the x axis, `engine_unit` the y unit, `columns.order` one
  trace per order, the row's extra outputs further traces. A plot unlike "value against energy or angle" adds one row
  to the plot table of `plots.js`.
- **A record kind**: none. Its section and form come from the schema; a new section becomes a new group of the Data
  tab.
- **A fit-parameter kind**: none. The Fit table's kind filter, units and owner columns come from
  `kinds()["parameters"]`.
- **A physics function**: none. Its form and plot come from its row, provided the row describes its outputs (names,
  units, which argument is the natural axis: finding 12).
- **A value type of the schema**: one widget in `forms.js` and its test (the consistency check fails until it
  exists).
- **A kind of plot**: one function in `plots.js`, one row of the plot table, one smoke-test step.
- **A kind of job**: one row of `jobs.py` (its pyazr call, settings, progress fields) and its test.
- **A session option**: none. The Project tab's options come from `kinds()["options"]`; an option of the calculation
  is a `<config>` key, edited with the other entries: its row of `io::ConfigKeys()` (with its title, group, absence
  and rarity) gives its control in the run settings.

## 11. The order of work

Each step leaves something a user can run. Sizes are estimates; "today" means it runs on today's pyazr behind one
adapter, `pyazr/ui/engine.py` (about 250 lines), the only file that changes when pyazr's `Session` replaces
`azure2` (pyazr step 5).

1. **Server and process model.** `pyazr/ui/` with `server.py`, `state.py` (text, history, revision, working file),
   `workers.py` (session worker, budget, death with the server), `jobs.py` (the table; recipe fits with progress,
   stop, kill, checkpoints), `routes.py` (the table, generic handlers), `errors.py`, `engine.py`; the routes of
   Project, records, levels, config, THM, evaluate, data. About 1500 lines; tests about 600 (the assertions of section
   9 on T1, T2). **Today.** Needs no plotly. **Done 2026-11-13** (section 11.1).
2. **The page shell, Project, Data, Nucleus.** `index.html`, `app.js`, `api.js`, `forms.js`, `tables.js`,
   `tabs/project.js`, `data.js`, `nucleus.js`, `style.css`. About 2200 lines; the widget consistency check; first
   smoke tests written (run where Playwright exists). **Today** (pair add/remove and channel generation wait).
   **Done 2026-11-13** with plotly.js vendored and a preview plot, smoke tests in headless Chrome (section 11.2).
3. **Results.** `plots.js`, `tabs/results.js`, the Results payload, grids, CSV/SVG/PNG export. About 900 lines; tests
   of the payload against the program's output files. **Today**, with plotly vendored (decision 11); `Evaluation`
   (pyazr step 3) replaces the adapter's calls. **Done 2026-11-13** (section 11.4).
4. **Fit.** `tabs/fit.js`, job kinds: recipe fits and MINOS (iminuit), the emcee sampler, bands; then the engine's
   fitters once pyazr has `fit_with_engine` with progress and stop (decision 8 and finding 2). About 1000 lines.
   Recipes need `Objective` (pyazr step 3; today `objective_and_grad` serves); engine fits need pyazr's call. **First
   part done 2026-11-13** (section 11.3); **the engine's fitters and MINOS done 2026-11-14** (section 11.6); **the
   sampler and the bands done 2026-11-14** (section 11.7).
5. **Physics and THM.** `tabs/physics.js` generated from the physics table (pyazr step 4; before it, stopping only),
   `tabs/thm.js` (today: `set_thm`, the session's THM calls). About 800 lines. **Physics tab done 2026-11-13**
   (section 11.5, stopping included); **THM tab done 2026-11-14** (section 11.8).
6. **Close-out.** `pyazr ui` as a command (pyazr step 7), static files and plotly in the wheel (pyazr step 8), the
   parity checklist ticked or dropped item by item, the smoke tests in their workflow, `numcmc` and the menus removed
   (Phase 9). About 400 lines and the checklist.

Steps 1–3 need nothing of pyazr 3.0; 4 needs pyazr step 3 for the recipes in their final form and pyazr's engine-fit
call; 5 needs pyazr step 4; 6 needs pyazr steps 7 and 8. Step 2 and 3 can run side by side once 1 is merged.

**For the owner's meeting next week**: step 1 and the parts of 2 and 3 that show the idea: `pyazr ui T2.azr`; the
Project tab with options and problems; the Data tab's list with one generated form (a data set) and the engine's
description beside it; the Nucleus tab's levels with a width edit, the χ² changing, undo; a recipe fit job with its
progress and Stop; the χ² per data set and one plot per data set. About 2500 lines, one agent, a week; the plots only
if plotly.js may be vendored by then (otherwise tables and CSV). Not in it: pairs added or removed, engine fits,
the sampler, Physics, THM.

### 11.1 Step 1 as built (2026-11-13)

`pyazr/ui/` (about 2000 lines in nine files) on today's pyazr: `python -m pyazr.ui project.azr [--port] [--threads]
[--token] [--root]`. Tests in `tests/pyazr/ui/` (four files, 33 tests; the rule each guards and the error injected
to show it fails are in its `README.md`); the user guide's page `docs/source/user_guide/web_ui_api.rst`.

**Served** (55 routes of 70; `GET /api/describe` lists them all with their status): describe (with
`pyazr.describe()` and the session options), server (PIDs and memory), schema, kinds, the physics functions' list;
the project's state, open, new, save, restore, text (get, replace), validate, describe, records
(list, get, add, edit, remove, rename), levels (list, get, add, edit, remove), channels, pairs (list, edit), the
reaction label, config, THM (get, options, experiments set and removed), undo, redo, history, apply; the session,
parameters, x (get, set), evaluate, the whole evaluation, one data set's arrays, transform; jobs (kinds, list, start,
get, stop, kill, accept). **501, a pyazr call missing**: adding and removing pairs (`Project.add_pair`,
`remove_pair`, finding 10), the model on a grid (`Session.evaluate`), a data set's grid and kernel (`DataSet.grid`,
`kernel`), the level matrix; jobs: the
engine's MIGRAD and LM and MINOS (`fit_with_engine`), the sampler (`Objective.log_prob`), bands, the rate
(`Session.reaction_rate`). **501, a later step** (pyazr has the call): system, channel functions, outputs,
EXFOR, ENSDF, files, events.  (Served since: a physics function on arrays, the channel functions, the
session's functions as rows and their calls, the stopping rows and their calls.)

**Measured on the product** (a server process over HTTP at one thread, `work/uiserver/measure.py`): opening T2 takes
1.2–1.5 s to the first answer and 1.2–1.8 s more for the first evaluation (worker start, import, session build
0.9–1.3 s); 13N 0.9–1.4 s and then 13.3 s the first time (the external-capture integrals into the UI's own output
directory), 1.9 s once they are cached. An evaluation through the API: 229 ms on T2 against 218 ms for the same
engine calls in one process, 140 ms against 124 ms on 13N; the numbers are equal bit for bit. Memory (resident):
server 78–80 MB, session worker 149–166 MB, job worker 187–204 MB (peak 192–211 MB). A stop request to "stopped":
0.07 s on T1, 0.3 s on 13N, 2.7–3.4 s on T2 (one objective with gradient); a kill 0.07 s. The workers of a server
killed with SIGKILL during a fit end within 0.05 s (`PR_SET_PDEATHSIG`).

**Where the product departs from this note**, each in one named place:

1. *The sessions' files lie beside the project* (`.<name>.ui-session.azr`, `.<name>.ui-job-<id>.azr`), not in the
   UI's directory: the engine resolves the prior covariance file against the directory of the file it reads, not
   against the session's base directory (`azure2(file, cwd=...)`), so a session file elsewhere does not open T2
   (`pyazr/ui/state.py`, `State.session_file`; a defect of the engine, reported).
2. *The sessions write into the UI's own output directory* `.<name>.ui/engine/`: the session's text is the project's
   with `<config>`'s output directory pointed there (`State.session_text`), so neither an evaluation nor a job
   writes into the user's `output/`. The χ² per data set is read from the `chiSquared.out` the session writes there
   (`pyazr/ui/engine.py`): the program's split, at the file's six significant digits; no pyazr call gives it at full
   precision (finding 5 stands). The data χ² and the objective are full doubles.
3. *Writing x into the text* (apply, accept) goes through `save_fit` (a snapshot file in the UI's directory, read
   back and checked by pyazr, then merged by lines into the project text: the lines where the session's text differs
   from the project's stay the project's; `State.merge_snapshot`), not an in-memory call (finding 8 stands).
4. *Checkpoints, logs and the URL file* are in `.<name>.ui/` (`jobs/<id>.json`, `logs/`, `url`), not under
   `output/`. The engine's notes go to the workers' log files there (finding 4 stands).
5. *A second job while the budget is full is refused* (409 with the running job's id), not queued; no `queued`
   state (`pyazr/ui/workers.py`, `Budget`). With N ≥ 2 threads each job gets one thread and the session worker one;
   that case was not run (budget 1).
6. *Routes added*: `GET /api/server`, `POST /api/project/restore` (the crash copy), `GET /api/jobs/kinds`,
   `POST /api/jobs/{id}/kill`, `GET /api/project/levels/{label}`, `GET /api/project/reaction`. A level is addressed
   by its label in the path (`{label:path}`, it may hold "/"); a channel spin as a decimal (0.5).
7. *Revisions*: every change (edit, undo, redo, restore) gives a new number; an edit without one answers 428.
8. *Errors*: a status is the first row of `errors.STATUS` that names the exception's class or any base (a failed
   lookup 404 though also a `ValueError`, a closed session 409, a missing file 404, any `pyazr.Error` 422); a worker
   sends the names of the class and its bases. `what` is the whole message; `allowed` and `fix` come from the first
   error diagnostic of a FormatError and are empty otherwise: pyazr's messages hold the three parts as one text, not
   as attributes (finding 9 stands).
9. *Session options*: the engine's table (`pyazr.kinds()["options"]`: keyword, key, type, default, flag, effect,
   doc) with the value in use and its `source` ("default", "project"), in `GET /api/describe` (from the project's
   `<config>`) and `GET /api/session` (from the session). Since decision 7 the options of the calculation are
   `<config>` keys: edited with `PATCH /api/project/config` like any entry (a contradiction is refused with the
   engine's code); the provisional `SESSION_OPTIONS` is gone.
10. *Workers* are `python -m pyazr.ui.workers` processes speaking JSON lines (not `multiprocessing`: no resource
    tracker process, a thread count per worker through its environment); every worker is started from one
    long-lived thread, since `PR_SET_PDEATHSIG` follows the thread that started the process.
11. *Engine notices*: the warnings of a request in the session worker (`pyazr.EngineWarning`, the engine's notices,
    among them) come back in the answer's `warnings` list (category, message); a job's go to its log.
12. *Origin*: a request with an `Origin` header other than the server's own (scheme, host and port) is refused (403):
    every port of localhost is one site for the browser, so the cookie alone would let a page of another local server
    send edits. A Host header `[::1]` with or without a port is a loopback name.
13. *Not built*: the event stream (polling serves the tests; the Fit tab long-polls a job, 11.3), the file browser, the
    interrupted job's evaluation beside "accept" (accept works on it).

### 11.2 Step 2 as built (2026-11-13)

`python -m pyazr.ui project.azr` prints the URL; the page is `pyazr/ui/static/` (about 2400 lines in 16 files and
the vendored plotly.js): `index.html`, `app.js` (header, the tab table, hash routes `#/data/<name>`,
`#/nucleus/level/<label>`, keys and the `?` overlay, themes), `api.js`, `tables.js`, `forms.js` (the widget table),
`formview.js` (drawing a form), `layout.js` (layout hints), `plots.js` (the plot table, palette, CSV), `ui.js`,
`tabs/project.js`, `tabs/nucleus.js`, `tabs/data.js`, `tabs/records.js` (a record's panel, shared by Data and
Nucleus), `selftest.js`, `style.css`, `vendor/plotly-cartesian-4.1.2.min.js` (1.5 MB, MIT, with its licence).
Tests: `tests/pyazr/ui/static_test.py`, `node_test.py` (node's own runner on `js/*.test.mjs`), `smoke_test.py`
(label `browser`; headless Chrome).

**What a user can do.** Project: the file, the validator's problems with links to the record, `<config>` edited as
text through the config route, the session options read-only (from `/api/describe`), the history with Undo and Redo,
the text with "replace text", the engine's line. Nucleus: pairs (edited through `set_pair`), levels with their
channels (energy, fixed, active; a channel's width and fixed flag; each edit one history entry; the status line
gives the χ² before and after), the `gamma_decay` and `potential` records through the generated forms. Data: every
section of the schema except those of the Nucleus tab as a group, the program's χ² per data set, one generated form
per record with the engine's description beside it and the validator's messages under their fields, add,
duplicate, rename, delete, the first rows of a data file, the preview plot (data with errors and the model at the
working x against c.m. energy or angle; plotly.js or, without it, a table; CSV export). Fit, Results, Physics and THM
are tabs that say which step brings them.

**Measured** (a server process over HTTP at one thread, `work/uipage/measure.py` and `walk.py`, outside the
repository): the requests of the first paint of the Data list (page, modules, plotly.js, describe, schema, project,
records, validate) 0.21 s on T2 and 0.63 s on the largest project of `tests/` (`reaction_rate/models/
moved_azr_truth.azr`, 900 kB); the first evaluation (worker start and session build) 0.96 s on T2 after it; a width
edit with the evaluation that follows 0.43–0.54 s on T2, its undo 0.49 s. The χ² after the width edit on T2
(7465.8873676603125) equals a direct pyazr session's on the edited text bit for bit, and Undo gives back
10905.2313020616 exactly. Headless Chrome, its start included, dumps the finished Data list in 2.0 s.

**Where the built page departs from this note**:

1. *Static files* are served by `STATIC_ROUTES` and one generic handler in `routes.py`, apart from `ROUTES` (whose
   rows are pyazr calls or server methods, and whose tests count them); the page revalidates every file on each load.
2. *The token*: `/?token=` answers with a redirect to `/` that sets the cookie, so the token leaves the address bar.
   The cookie is `SameSite=Lax` (since 2026-11-14; it was `Strict`): the printed URL is usually followed from
   somewhere else (a link in a chat or a mail page), and a browser keeps a strict cookie back on a navigation that
   started on another site and on the redirect that follows it, so the server answered 401 to its own link (met by
   the PI; reproduced in headless Chrome, `smoke_test.py`). A lax cookie goes with a navigation (GET) from elsewhere
   and with nothing else from there. This is safe while **every GET route only reads** (the rule for a new route:
   a request that changes something is not a GET; it then carries an `Origin` header, checked by item 12). The
   cookie's name carries the port the browser named (`azure2_ui_token_28400`): a browser has one cookie jar per
   host whatever the port, so two servers on one machine overwrote each other's token. A navigation without the
   token or with a wrong one (a GET that accepts HTML, outside `/api`) gets a small page that says so and takes
   the token in a field (its form is `GET /?token=`, the printed URL's door; its own policy allows its inline
   style and that form, no script); every other request keeps the JSON refusal with its three parts.
   *Content-Security-Policy* on every answer: scripts from the server only (no inline script, no eval); inline styles
   are allowed because plotly.js sets style attributes and adds style elements; images also from `data:` and
   `blob:` (plotly's image export).
3. *A route added*: `GET /api/files/{path}` (the first rows of a file inside the root, `Server.file_head`); the
   listing and file browser (`GET /api/files`) stay a later step.
4. *No event stream*: the page polls `GET /api/project` every 3 s while visible; a newer revision redraws.
5. *Self-test and smoke test in headless Chrome*, not Playwright: `#selftest` visits every view and writes its result
   into the document; `#selftest?edit=1` also types a width and presses Undo. Screenshots use Chrome's old headless
   mode, whose viewport is the window size.
6. *Forms*: a widget splits a field's text into boxes and joins it again; a text it cannot split is one text box. The
   engine validates; conditions `field_present` and `field_is` hide fields, the other conditions (`fold`,
   `single_reaction`, ...) show the field with their text. `<config>` is edited as key and text (the schema has no
   table of its keys: finding 13).
7. *Plots*: the model is drawn at the data points (the engine's values there), not on a dense grid (step 3); a data
   set is drawn against energy, or against angle when its energy is one value; points of a few distinct angles are
   one trace each. Legendre orders are not separated (finding 14).
8. *Pairs* are shown with their spins as 2J; the channel width is shown "as written" (finding 15).

**Extensibility as built** (section 10): a new value type: one row of `WIDGETS` in `forms.js` (the static test and
the node test fail until it exists); a new record kind: nothing (its form is generic; a line in `layout.js` only to
put fields on one line; `gamma_decay`, `beam`, `normalization` have no hints and get full forms); a new section:
nothing (a new group of the Data tab), unless the Nucleus tab should show it (`NUCLEUS_SECTIONS` in
`tabs/records.js`); a new data kind: nothing for the list, the form and the preview plot when it is a value against
energy or angle (axis rule from the quantity row's `angle_axis`, units from the dimensions); another kind of plot: a
row of `PLOTS` in `plots.js` and its node test; a new tab: a module `tabs/<name>.js` with `render()` and a row of
`TABS` in `app.js` (the self-test visits it; add it to `VIEWS` of `smoke_test.py` for a screenshot).

### 11.3 Step 4, first part, as built (2026-11-13): the Fit tab

`tabs/fit.js` (the tab, the parameter table, following a job), `tabs/fitjob.js` (the job form, the job's panel and
result, the jobs list, the widget table of setting types), `fitting.js` (no document: rows, choices, the progress
and result views), `selftestfit.js` (the self-test's walk through the tab); a row of `PLOTS` (the objective against
the call, log y). Server: `jobs.py` (settings as typed rows, iminuit's MINOS as a job kind, covariance and
correlation, the trace, the log tail, the long poll), `server.py` (the parameter table and x readable while a job holds
the budget; the accepted revision kept in the job's record). Tests: `jobs_test.py` (5 new), `js/fitting.test.mjs`
(8), the smoke test's walk through the tab.

**What a user can do.** Choose a job from the server's job table (the kinds and methods as the server lists them;
the rows that wait for pyazr greyed with the call they wait for), set its settings (drawn from their types and
defaults), its start point (the working x, the text's values, a job's best point) and Run; follow it (calls, the
objective at the last call and the best, time, the session build, the objective against the call on a log scale, the
last lines of its log), Stop it between calls or Kill it; read the result in the recipe's own words, the objective
at the start and at the best point, the χ² and N at the best point (evaluated on the job's revision); **Accept into
project** (one history entry; Undo takes it back), after which the χ² of the project as written stands beside the
fit's with the known causes of a difference; start iminuit's MINOS from an iminuit fit. The table of every parameter:
filters by kind (the kinds present, in the order and with the docs of `pyazr.kinds()`), free only, text; the working
x of a free parameter edited in place (the header's χ² follows); **Reset x**; the job's best point, HESSE errors,
MINOS intervals and, on demand, its start point. The covariance and correlation as tables and CSV.

**Measured** (`work/uifit/walk.py` on T2, a server process over HTTP at one thread): a job's GET 2.5 ms (2.1 kB with its
trace); a stop request answered by the waiting GET 0.88 s later (the objective call in flight, with its gradient, about
2 s on T2); an iminuit fit from the stopped point to its end 416 calls in 779 s, objective 7087.924 → 298.917; χ² at the
best point 295.0675134652589, of the project as written after Accept 295.06757397066707 (T2 fits a target density:
likely the open point above; on T1, which does not, the two differ by 3 × 10⁻¹⁵), Undo back to 10905.2313020616 exactly.

**Where it departs from this note:**

1. *No event stream.* A running job is followed by a long poll: `GET /api/jobs/{id}?wait=s` answers when the job's
   state changes or after s seconds (at most 10); the page keeps one such request in flight with s = 1, so it shows
   the progress once a second and a Stop, a kill or the end at once. Server-sent events were not built: they need a
   streaming handler beside the generic JSON handlers of the route table, the test client and agents with curl read a
   long poll as any other answer, and headless Chrome's virtual time (the smoke test) is paused while a request is in
   flight, which a never-ending stream would make endless. The revision is still polled every 3 s.
2. *The job's trace* (call, objective, best so far; at most 500 points, thinned by two past that) and the last eight
   lines of its log are in the job's record; `GET /api/jobs` leaves them out.
3. *A job's record names the revision of its text `text_rev`*: an answer's `rev` is the project's current revision,
   which the page re-reads on (a job's older `rev` made it re-read without end).
4. *Settings are typed rows* (`Setting`: name, type, default, doc; types `integer`, `boolean`, `names` in
   `SETTING_TYPES` with their conversion); a request's value is converted or refused (422) with what the type allows.
5. *MINOS*: iminuit's, as the job kind `minos/iminuit` from a fit's best point (MIGRAD there, HESSE, MINOS on the
   parameters named); `follows` in a row makes the page offer it beside the result of the row it follows. The
   engine's MINOS stays a row that waits for `fit_with_engine`.
6. *Errors and covariance* come from iminuit's HESSE only: scipy's L-BFGS-B gives an approximate inverse Hessian,
   not shown.
7. *A fixed parameter* is not edited from the table: its row does not name the record or level that holds it
   (finding 17); the page links to the record where the parameter's name begins with a record's name.
8. *The "physical value"* is the session's at the text's start values; its unit is not sent (finding 18). The prior
   column shows the parameter row's own prior (density, Doppler temperature, abundance, ...); normalizations and
   shifts carry their priors in their records only (finding 19).
9. *The objective* (χ² plus every prior term) is what the recipes minimise and what the progress shows; the data χ²
   and N at the best point come from an evaluation there.
10. *Threads*: one per job (the budget's); the form shows it, it is not a setting.
11. *Screenshots of a running job* need Chrome's real-time `--timeout`, not its virtual-time budget (the long poll is
    always in flight).

**Extensibility as built.** A new parameter kind: nothing (the kind filter, its doc and the units come from
`/api/parameters` and `/api/kinds`). A new fit recipe: a function and a `JobRow` in `jobs.py` (and its test); the
form, the table of jobs and the result follow from the row. A new job kind: the same; its result shows the keys the
page knows (`errors`, `covariance`, `correlation`, `minos`), other keys need a line in `fitting.js`'s `resultView`. A
new option of a job: a `Setting` of its row; a new type of option: a row of `SETTING_TYPES` (`jobs.py`) and of
`SETTING_WIDGETS` (`tabs/fitjob.js`).

### 11.4 Step 3 as built (2026-11-13)

The Results tab (`pyazr/ui/static/tabs/results.js`), its plot rows in `plots.js` (`results_against_x`,
`angular_distribution`, with `X_CHOICES` and `Y_CHOICES`), the Results payload in the adapter
(`Engine.evaluation`), the output-file routes. Tests: `tests/pyazr/ui/results_test.py`, `js/results.test.mjs`; the
self-test visits `#/results`, its table view and `x=e_x&levels=1&logy=1`. The director's corrections to step 2
(bounded file preview, tables that wrap, repeated groups as small tables, the full χ² in the status line, paths
broken at `/`, the tests registered above the last block of `tests/pyazr/CMakeLists.txt`) came first.

**The payload** (`GET /api/evaluation`, at the working x; 69 kB on T2, 92 kB on T1, 101 kB on 13N): the program's
χ² split (the rows and the Total- line of the `chiSquared.out` the session writes), the objective and its terms;
per data set from one forward pass the c.m. energy, excitation energy and c.m. angle of each point after energy
shifts, the datum and error as read and times the norm in use (`segment_norms`), the model at the point, the
engine's standardized residual (m − n·d)/(n·σ) (`residuals`) and the sum of their squares (`segment_chi2`), the S
factor (model, normalized datum and error times `calculated_conv`) where the engine unit is a cross section and the
entrance pair charged; per active extrapolation the same on its grid (a second session made with `extrapolations=True`); the model data
(`pyazr.model_data`); the residual's definition. Every column equals the program's `AZUREOut_*` (menu 1) and
`.extrap` (menu 3) to their eleven digits on T1 and T2, and the χ² table equals `chiSquared.out` (`results_test`).
The squares of the residuals add up to each data set's χ² of `chiSquared.out` (to its six digits), except for the
two data sets of T2's data covariance, whose residuals are whitened rows: their squares add up to the group's χ²
(513.896) but split it 100.881 / 413.015 against the program's 99.1682 / 414.728 (finding 5); the table shows both
and the page names the covariance file. Load (a fresh server, the session's build included; then again): T2 1.15 s
/ 0.54 s, T1 7.0 s / 0.07 s (the extrapolations' external capture the first time), 13N 10.1 s / 0.18 s.

**Where it departs from section 6.5**:

1. *The model line of a data set* joins the engine's values at the data points; a resonance narrower than their
   spacing is not drawn, and the page says so. Extrapolations are drawn on their grids, an `adaptive` grid refined
   by the engine. A dense grid between data points waits for `Session.evaluate` (finding 21).
2. *x*: c.m. energy, excitation energy, c.m. angle. The lab energy is absent (finding 22). *y*: as measured (the
   datum times the norm, which the model is compared with) or the S factor; "cross section" of a data set measured
   as an S factor is absent (finding 23).
3. *One plot per series* by default; **together** draws one plot per y unit. Every trace is in the legend.
4. *Level markers* on the excitation-energy axis only (the project text's levels: E_x needs no conversion).
5. *"At"* is the working x; the job / x0 selector of the Fit tab is not on this tab yet (a sentence in its place,
   no control).
6. *More than 5000 points*: drawn as they are (SVG, no decimation), and the plot says so; no project of `tests/`
   comes near (13N: 593 points).
7. *Routes*: `POST /api/outputs` (was a later step) writes the program's output files at the working x into the
   session's directory and copies them into the project's `<config>` output directory, the one session request
   that writes into the user's directory; `GET /api/outputs` lists that directory (added). `GET /api/evaluation`
   and `/api/evaluation/data/{name}` carry the payload above (energies now after the energy shifts).
8. *Legendre orders* are one series (finding 14 stands, said on the page). The reaction rate and the bands are a
   panel saying what they wait for.

**Extensibility as built**: a new data kind: nothing (series, axes from the dimensions, the S factor from the
engine unit's dimension and the pair's charges, χ² and model data from the engine); a new kind of plot: one row of
`PLOTS` (or a builder like `resultsPlots`) and its node test; a new x or y choice: one row of `X_CHOICES` /
`Y_CHOICES` in `plots.js` and, if it needs arrays the payload lacks, their lines in `Engine._dataset_arrays`; a new
section whose records are evaluated: a group of the selector by itself; its arrays need one method of the adapter
beside `_extrapolation_arrays` and a key of the payload read by `seriesOf` in `tabs/results.js`.

### 11.5 Step 5, first part, as built (2026-11-13): the Physics tab

`tabs/physics.js` (the tab: the list, the form, the answer drawn), `physicsform.js` (no document: the rules of
section 6.6, the table `ARGUMENT_KINDS`), two options of shared modules (`renderForm(..., {whole})`: every field's
text to the caller, "Defaults" instead of "Revert"; `drawPlot(..., {csv})`: the CSV text from the caller), the routes'
side in `pyazr/ui/physics.py` (the axis, choices by name, the arguments as used, the Python text, the session rows),
two rows of `routes.py`, and the rows' new fields in the binding (`group`; `example` and `sweep` per argument).
Tests: `tests/pyazr/ui/physics_tab_test.py` (4), `js/physics.test.mjs` (9), the self-test's four Physics views and
the smoke test's checks of them, a check of the rows in `tests/pyazr/physics_table_test.py`.

**What a user can do.** Pick a function in the list (grouped as the table groups its rows, with a search box), see
its signature, summary, formula and library function, choose the argument to plot against (any with a sweep), edit
the generated form (the axis as from, to, points, spacing; several values of any other argument give one line each,
"all" of a choice, "both" of a flag) and press Plot; every function plots at once with its start values. The
answer: the status line ("400 of 444 elements computed; rho outside its domain (finite > 0): 44"), the fields to
show, log x and log y, one plot per unit with a trace per field and line and gaps where an element is not computed,
CSV of the arrays as returned (one line per element with the arguments at it and the status), SVG and PNG, and the
pyazr call that makes the same arrays (shown, copied, or saved as a .py file). "this project": the session's channel
functions of a pair of the project (its pairs offered by key and label), refused with the engine's words when the
project does not build. A view is a link: `#/physics/coulomb_waves?from=-2&to=20&spacing=linear&v.l=0 1 2 3`.

**Where it departs from section 6.6**: none in the rules. The radius of `channel_functions` defaults to 0 (the
pair's channel radius) and `J` is not offered (a pair with a spin-orbit potential is refused with the engine's
message). Its Python text opens the saved file: the session is built from the working text, which equals it when
nothing is unsaved (finding 1).

**Extensibility as built** (the PI's measure). A new physics function: one row of `BuildTable()` in
`api/src/pybind_physics.cpp` with its `group`, an `example` for each argument without a default and a `sweep` for its
axis (the check of `physics_table_test.py` fails without them), and its `row_test`; nothing in the page or the server
(the node test draws a form for every row and for one injected there). A new kind of argument (a choice among named
potentials): its kind in the binding's `ArgKind` and `KindName`, how the row's adapter converts it, and one row of
`ARGUMENT_KINDS` in `physicsform.js` (the field of the form, how its text becomes values; the node test fails until
the row exists); a choice whose names come from the project (as the session's `pair`) also says where the page reads
them (`fctx` in `tabs/physics.js`). A new table of functions (the stopping rows): one row of `SOURCES` in
`tabs/physics.js` and a route that answers rows of this shape.

**Stopping, as built (the second part).** The row types of the physics table (`ArgSpec`, `OutSpec`, the domain
check, their description as data) moved into `api/src/pybind_rows.h`, shared by `pybind_physics.cpp` and
`pybind_stopping.cpp`; `ArgKind` gained `text` (one string per call, not broadcast). `pybind_stopping.cpp` has a
table of five rows (`stopping_power`, `range`, `energy_after`, `areal_density_between`, `energy_straggling`): the
material (`formula` text, `phase` choice, `compound_correction`) and projectile (`z1`, `mass_u`) arguments first,
then the row's own; each row carries `group` "stopping", SRIM's `citation` and `implementation`. Each adapter calls
the library function `pyazr.stopping`'s functions call, with the material and projectile made as they make them, so
the values are theirs bit for bit (tested). A library refusal of one element (a layer at least the range) is status
-2; a formula the library refuses is a refusal of the request (422). Python: `pyazr.stopping.table()` and
`evaluate(name, formula, **arguments)`; `kinds()["stopping"]` is the table (it was a signature and a doc per binding
function). The page lists the rows under "stopping" (one row of `SOURCES`, one row `text` of `ARGUMENT_KINDS`) and
shows the citation under every stopping plot, in the CSV's first line and in the Python text. Not done: the
provenance per element (`pyazr.stopping.provenance`) is not shown; user stopping tables of a project are not
reachable (the rows evaluate SRIM's stopping only).

*Extensibility of stopping*: a new stopping function is one row of `BuildStoppingTable()` in
`api/src/pybind_stopping.cpp` (its own arguments with unit, domain, example, and a sweep for its axis; returns with
units; an adapter calling the library); `physics_tab_test.py` plots every row without typing and runs its Python text,
the node test draws its form; nothing in pyazr, the server or the page.

### 11.6 Step 4, second part, as built (2026-11-14): the engine's fitters and MINOS as jobs

`jobs.py` (the rows of the engine's fitters made from its table, `engine_rows`; the job function `fit_engine`; the
whole table `job_rows()`; a thread of the job's worker that reads the stop request, `listen_for_stop`; progress
numbers beyond the four common ones passed on by name; non-finite numbers sent as null), `engine.py` (`open_for_job`,
`objective_of`, `fit_with_engine`: the adapter's three functions for a job's worker), `fitting.js` (`factsText`, the
progress view's further numbers, where a job's errors come from), `tabs/fitjob.js` (the facts under the method and
as a column of the job table, the further progress numbers). Tests: `tests/pyazr/ui/engine_jobs_test.py` (6), the
generated rows against `pyazr.kinds()["fitters"]` in `routes_consistency_test.py`, `js/fitting.test.mjs` (2 new).

**What a user can do.** Choose an engine fitter in the Fit tab's job form: the rows are the engine's fitter table
(`pyazr.kinds()["fitters"]`), the program's default fit (MIGRAD) first and the form's default, then the others in the
table's order (`engine_gsl-lm`, `engine_lm`, `engine_migrad-gradient`), then the recipes. Under the method: the
engine's doc of the fitter and its facts (the program's switch, whether it gives a covariance, whether MINOS can
follow, gradient, residual rows); the job table has them as a column. Settings: `max_calls` (the program's
`--max-calls`; MIGRAD's and NLopt's call limit, the least-squares fitters stop at their own iteration limit),
`minos` on the rows whose fitter allows MINOS, the start point as for the recipes. Run follows the job as a recipe is
followed: calls (the fit's, then HESSE's and MINOS's counted on), the objective at the last record and the best,
the stage (`fit`, `hesse`, `minos`) and the fitter's measure (LM's damping) by name, the objective against the call,
the log. **Stop** sets the engine's stop request: the fit ends at the best point it has seen, state `stopped`, the
best point kept and acceptable. The result: the engine's status in its words ("the engine's migrad: converged",
"... call_limit"; the state `converged`, `stopped`, or `ended` for any other status), the objective at the start and
at the end point, the χ² and N at the end point (an evaluation, as for the recipes), the errors (as the program
writes them into `param.sav`), covariance and correlation (HESSE after MIGRAD; (JᵀJ + priors)⁻¹ after LM and GSL-LM),
the MINOS intervals when asked; **Accept into project** and Undo as for the recipes. After an engine MIGRAD fit (with
or without the gradient) the result offers **minos (engine_migrad)…**: the engine's MINOS as a job from the fit's
end point (MIGRAD from there, HESSE, MINOS on every free parameter).

**Measured** (`work/uifit2/walk_engine.py`: a server process over HTTP at one thread, every engine fit row from the
text's start; then `Session.fit_with_engine` in another process on a copy of the project, and the program with the
row's switch, menu 2):

| project, row | state | job's objective | `fit_with_engine` | program (12 digits) | calls | wall time |
|---|---|---|---|---|---|---|
| T3 `engine_migrad` | converged | 509.3125016649391 | equal | 509.312501665, equal | 276 | 0.7 s |
| T3 `engine_gsl-lm` | converged | 509.31249471170713 | equal | 509.312494712, equal | 31 | 0.5 s |
| T3 `engine_lm` | converged | 509.3124947112076 | equal | 509.312494711, equal | 9 | 0.4 s |
| T3 `engine_migrad-gradient` | converged | 509.3124963010905 | equal | 509.312496301, equal | 180 | 0.7 s |
| T2 `engine_migrad` | converged | 297.8603087301522 | equal | (no switch for it, item 7) | 1109 | 48.3 s |
| T2 `engine_gsl-lm` | converged | 297.8603087301522 | equal | 297.86030873, equal | 861 | 37.7 s |
| T2 `engine_lm` | converged | 297.8603087301522 | equal | 297.86030873, equal | 861 | 37.8 s |
| T2 `engine_migrad-gradient` | converged | 297.8603085507331 | equal | 297.860308551, equal | 421 | 31.3 s |

"equal" for `fit_with_engine` is equality of the binary numbers (the same code from the same start in a fresh
session), for the program equality of the printed 12 digits. The calls include HESSE's after MIGRAD. On T2 LM and
GSL-LM run MIGRAD in their place (a data set without an analytic Jacobian; the result says so), as the program does.
The session builds in 0.2 s on both (the Coulomb table). Stop on T2 (`engine_migrad`, after 20 calls, 0.059 s per
call): the waiting GET answered `stopped` 0.205 s after the request, one progress record after it (the evaluation
running then), the engine's status `stopped`, best objective 7087.924186759944 (from 10910.13 at the start). Accept
of that stopped job: χ² at its best point 7083.021646029943, of the project as written 7083.021646029921 (3 × 10⁻¹⁵
relative: the written digits); Undo: 10905.2313020616, the start χ², exactly. On T3 after `engine_migrad`: χ² at the
best point 509.3125016649391, as written 509.31250166493913; Undo exact. `engine_jobs_test.py`: the four fitters on
T3 equal `fit_with_engine` bit for bit (end point, objective, covariance, errors); the stop of MIGRAD and of LM on T2
answered 0.16-0.22 s after the request with one record after it.

**Where it departs from this note:**

1. *The job's best point is the fitter's end point*, with the engine's objective there; while it runs, the best is
   the lowest objective reported (HESSE and MINOS evaluate near the minimum and can report a lower one).
2. *The start objective* is one evaluation at x0 after the fit: an evaluation before it would change the state the
   engine's pooled working copies carry into the fit (applied energy shifts, sub-point grids; `fit.md` section 4), and
   the job would no longer end where `fit_with_engine` ends. The same x evaluated in the fit and after it may differ
   in the last digit (3 × 10⁻¹⁶ relative on T2).
3. *Calls* of the job count the fit's records, then HESSE's and MINOS's on from there (each stage of the engine
   counts from 1); `fit_calls` in the result is the engine's count of the fit alone. LM and GSL-LM report per
   iteration, so their trace has a point per iteration.
4. *MINOS validity* is not given by the engine: the interval's `valid` is null.
5. *The stop request* is read by a thread of the job's worker (a blocking read of its input), not polled before each
   objective call: it reaches the engine while an evaluation runs (the engine releases the GIL), and a recipe sees
   it before its next call as before. `reports_after_stop` in the result counts the progress records after the
   request (the evaluation, or LM's iteration, running when it came).
6. *`max_calls` is offered on every engine row*: the fitter table has no column saying which fitters honour a call
   limit (the doc of `fit_with_engine` says MIGRAD and NLopt); a column would let the row say so.  The engine prints
   no notice when a least-squares fitter ignores it.  Where LM or GSL-LM run MIGRAD in their place, the result has
   `fitter` (`asked`, `ran`) and `notices` (the engine's line of the fitter's log), shown as a notice.
7. *Equality with the program* is measured against the program's printed total (12 digits) with the fitter's switch;
   on T2 the program fits THM data with `migrad-gradient` when no switch is given, so the plain MIGRAD row of T2 is
   compared with `fit_with_engine` only.

**Extensibility as built.** A new fitter of the engine: a row of its table (`src/Fitters.cpp`) and its binding; the
page shows its row, settings, facts and (where its row allows MINOS) a MINOS job after it, with no line of
`pyazr/ui` changed (`routes_consistency_test.py` builds the rows of a table with an extra fitter). A new recipe: a
function `(tracker, x0, settings, names)` and a `JobRow` of `JOBS`; a job that calls the engine otherwise uses
`tracker.session` (the adapter) and reports with `tracker.count(calls, objective, x, **numbers)`. A new setting type:
a row of `SETTING_TYPES` (`jobs.py`) and of `SETTING_WIDGETS` (`tabs/fitjob.js`).

### 11.7 Step 4, second part, as built (2026-11-14): the sampler and the bands

`jobs.py` (the job functions `sample_emcee`, `bands_covariance`, `bands_draws`; the chain's summary
`chain_summary` with the split R-hat; the start walkers inside the fit's limits; the setting type `number`; a row's
`takes`, and `follows` with any method of a kind; the state `finished`; the tracker's `source`, `product` and
progress without an objective), `engine.py` (`band_points`, `models_at`, `linear_sigma`: the model of every data set
and extrapolation, at one point and at many, and the band's sigma through `pyazr.bands.uncertainty_bands`),
`server.py` (a job started from a job gets that job's record as its source; the project's output directory for the
files a job writes; the list of jobs leaves a band's arrays out), `fitting.js` (`follows`, `chainRows`, the state),
`tabs/fit.js` (the chain's summary, a band job's result), `tabs/fitjob.js` (the widget of a number), `plots.js` (a
trace of role `band`: a shaded region in plotly, two columns in the CSV), `tabs/results.js` (the newest band job's
band on the series of the same name and points, a flag "band", the band's source in the head). Tests:
`tests/pyazr/ui/sampler_bands_test.py` (4), `js/fitting.test.mjs` (1 new), `js/results.test.mjs` (1 new).

**What a user can do.** After any fit, **sample (emcee)…** fills the job form for the sampler from the fit's best
point: walkers (empty: twice the free parameters plus two), steps (stored), spread (relative, around the start
point), thin, seed, discard (the burn-in left out of the summary; empty: the first half). The job follows as a fit
does: the objective at the best walker, the best so far, the step, the mean acceptance and a running autocorrelation
time τ; **Stop** ends it between stored steps with the chain so far. The result: the chain in
`output/chain_<id>.npy` of the project (stored steps × walkers × parameters) and its summary beside it as JSON; on the
page the mean, standard deviation, the 2.5, 16, 50, 84 and 97.5 % quantiles, τ, the effective draws and the split
R-hat of each parameter, the draws' covariance and correlation (the parameter table's error column is the chain's
sd). A sampler started from a sampler job (its best point as the start) continues that chain: the new file holds the
old steps, unchanged, and the new ones. After a fit with a covariance, **bands (covariance)…** runs the linear band
at the fit's best point (σ² = gᵀCg, `pyazr.bands`, over the level parameters); after a sampler, **bands (draws)…**
evaluates the model at draws of its chain after the burn-in and takes the 16th, 50th and 84th percentiles at each
point. Every data set and extrapolation gets its band, on its own points, in the engine's unit and as an S factor
where the Results tab offers it. The Results tab draws the band of the newest band job that finished as a shaded
region under each series (flag **band**), says which job and text it belongs to, and exports its edges in the
CSV.

**Measured** (`work/uifit2/walk_sampler.py`: a server process over HTTP at one thread on a copy of T3 with two free
parameters, energy_5 and width_5_1 of the level 2⁻ 19.0 MeV; every other level and width fixed):

- The engine's MIGRAD from the text: objective 516.240237935057; HESSE errors (Laplace) energy_5 4.36245e-4 MeV,
  width_5_1 1.05673e-3 MeV^1/2. The posterior's standard deviations on a 61 × 61 grid of ±5 σ (exp(−objective/2)
  summed, outside the server): 4.36251e-4 and 1.05676e-3, the Laplace values to 3 × 10⁻⁵ (the posterior is Gaussian
  here; on the level 2⁺ 3.03 MeV instead, a Park wall 1 σ below width_1_1 makes it one-sided: chain and grid agree
  there to 3 %, and the grid's lie 7 % and 17 % under the Laplace values).
- emcee from the fit's best point, 8 walkers, 6000 stored steps, spread 10⁻⁴: 56.1 s (48008 log-probability calls),
  acceptance 0.713, τ 31.8 and 30.3 steps, split R̂ 1.012 and 1.021, about 760 and 790 effective draws after the
  burn-in of 3000 steps. Chain sd over Laplace sd: **0.9989** (energy_5) and **1.0011** (width_5_1).
- The band from 1000 draws (3.3 s) against the linear band from the fit's covariance, the ratio of their half widths
  at all 333 points of the 9 data sets and 9 extrapolations with a band: median **0.994**, 5-95 % 0.965-1.014, every
  point within 10 %.
- `sampler_bands_test.py` (8 walkers, 2500 steps, 400 draws; 31 s): ratios 1.0002 and 1.038; bands median 0.999,
  96 % of the points within 10 %; a stop answered 0.06 s after the request with one step after it; the chain
  continued with its 20 new steps after the old ones, unchanged.
- Screenshots: `work/uifit2/tmp189154ig/fit_sample.png` (the sampler's panel and result), `results_band.png` (the
  band of the draws around the model after Accept).

**Where it departs from this note:**

1. *A band job takes its source from its start point*: `x0: {"job": id}` names the fit (with a covariance) or the
   sampler (with a chain), and the server hands that job's record to the new job (`takes`: the key its result must
   have; refused with 422 otherwise). No new setting type names a job, and the page offers the band beside the
   result that has what it takes.
2. *The linear band spans the level parameters only*: `pyazr.bands.uncertainty_bands` propagates the covariance's
   block over the columns a band spans (`rmatrix_columns`: energies, widths, the THM coherent background); the band
   from draws moves every free parameter of the chain, nuisances included.  Each band job's result says which
   (`varied`, `free`), and the Fit and Results tabs show them. The two agree where no nuisance is free
   (the evidence above); on a project with free shifts or densities they differ by the nuisances' share.
3. *The band's centre*: the linear band is the model at the fit's best point ± σ; the band from draws is the median
   of the draws. The Results tab draws the band around the model of the band job's start point and the model line at
   the working x: they coincide once that point is accepted into the project.  The legend names the band's job and
   the job whose best point it is drawn around; where the text is an accepted best point and the band is around
   another, the Results tab offers a band job from the accepted job.
4. *A sampler that ends* (all its steps) is in the state `finished`, not `converged`: a chain does not converge by
   a criterion of the job. Its burn-in is a setting; the summary does not test the chain's length against τ (the
   page states τ and the effective draws; the user guide's rule of 50 τ applies).
5. *One process*: the sampler runs without a pool (the budget gives a job one thread). The split R-hat treats each
   half of a walker's series as a chain; walkers of an ensemble sampler are not independent, so a value near 1 is
   necessary, not sufficient.
6. *A killed sampler* keeps the chain written at its last save (every 50 stored steps); its record has no result,
   and a sampler started from it continues that file (`Tracker.source_product`), as after a restart of the server
   (an `interrupted` job); the page offers **sample (emcee)…** beside any sampler job (the row follows a fit and a
   sampler).
7. *A start outside the fit's limits* is refused by the job (state `failed`, `ParameterError` naming each parameter,
   its value and its limits), not by the request: the limits are the session's (`Objective.bounds`), which only the
   job's worker has.

**Extensibility as built.** A new sampler: a function `(tracker, x0, settings, names)` that calls the objective it
needs through `tracker.session` (`objective_of`), reports each step with `tracker.count(calls, objective, x,
**numbers)` (the numbers appear on the page by name), stops when `tracker.stop_requested()`, writes its files at
`tracker.product(name, suffix)` and returns the summary keys the page knows (`quantiles`, `errors`, `covariance`,
`chain`); one `JobRow` with `follows` and `takes`. A new band: the same, returning `bands` in the shape of
`_band_entries`; the Results tab draws it with no change. A new setting type: a row of `SETTING_TYPES` and of
`SETTING_WIDGETS` (`number` was added so).

### 11.8 Step 5, second part, as built (2026-11-14): the THM tab

`tabs/thm.js` (the tab: options, experiments, channels, diagnostics), `functionview.js` (a row computed and its
answer drawn, moved out of `tabs/physics.js` and shared by both tabs: status line, facts, tables, plots, CSV, Python
text), `physicsform.js` (an argument kind `experiment`; a row without an axis answered as tables; a field with lines
of its own; a row's own bound of elements and the fields left off the first plot), `app.js` (a row of `TABS` may say
when it is `shown`), `pyazr/ui/thm.py` (the `<thm>` keys as record kinds, `KEYS`; the five THM calls as rows,
`THM_ROWS`, and how each answer is shaped, `ANSWERS`), `engine.py` (`Engine.thm_*`, the session's pairs in its
description), `physics.py` (the session's rows from two tables, text kinds, status words per row), `routes.py` (five
session rows and `GET /api/thm/keys`).  Tests: `tests/pyazr/ui/thm_routes_test.py` (8), `js/thm.test.mjs` (4),
`js/physics.test.mjs` (3 new), the self-test's THM views and the smoke test's checks of them.

**What a user can do.** On a project with THM content (a `<thm>` section, or a data set with a `thm` field) the tab
**THM** is listed; on another its link is hidden and `#/thm` says why (setting an option there adds the section).
*Model options*: a generated form of the global keys (vertex, kinematics, entranceL, coulombIntegral,
spectatorEnergy; a key the table lacks, as `spectatorEnergy[2]`, is a text field of its own, and "other key" adds
one); Apply is one `PATCH /api/project/thm`.  *Experiments*: their table (name, data sets that name it, keys as
written); a row opens its generated form (the twenty keys of an experiment line, with choices where the grammar has
them and its grammar where it is free text; `psNodes` and `spectatorAngleNodes` only with their window) with "data
sets that join" (set in the same edit) and **Remove**; **+ experiment** adds one by name.  Each edit is one history
entry, Undo takes it back, a refusal shows the engine's words.  *Channels*: the session's THM pairs (binding energy
B, radius) and their channels (level, l, s, whether the width is an amplitude), read-only with the reason (finding
11).  *Diagnostics*: the five THM calls of the session as links; each has the generated form of its row (the
experiment a choice of the project's, the energies an axis) and its answer drawn by the Physics tab's code: the
status line, the call's facts with units (B, E_aA, the window, the distortion's settings, ...), its tables (the
experiments' profiled norms and backgrounds with sigma, the coherent background's parameters, a background's
covariance), one plot per unit with a trace per channel and level (vertex), per exit pair and level (line shape) or
per field (distortion), the CSV of the arrays, Copy as Python.  A view is a link: `#/thm/LC09?fn=thm_vertex`.

**Measured** (`work/uithm/walk_thm.py`: a server process over HTTP at one thread, pinned to two CPUs with the
browser; the direct `Session` calls on the same file at x0; 200 energies from 0.05 to 2 MeV):

| project | kind | build | experiments, background, vertex | line shape | distortion | option edited |
|---|---|---|---|---|---|---|
| T2 (LC09) | plane wave | 0.5 s | equal, 0.04-0.06 s | 422 | 422 | kinematics=triple |
| f19_pag_thm | plane wave, linear bkg | 0.4 s | equal, 0.02-0.06 s | 422 | 422 | kinematics=lambda32 |
| 18O copy, cbkg | coherent background | 0.3 s | equal, 0.01 s | 422 | 422 | vertex=perlevel |
| 18O copy, dist | Coulomb distortion | 2.7 s | equal, 0.02 s | equal, 0.02 s | equal, 3.2 s | vertex=onshell |
| o18_lacognata2010 | DW vertex, optical | 6.4 s | equal, 0.01 s (98 unreached) | 422 | 422 | kinematics=triple |

"equal" is equality of every number of the session's answer (the raw route) and of every number the tab shows
(fields, facts, tables), as binary doubles; 422 is the engine's refusal of a call the experiment has no key for
("has no line shape (lineshape=on)").  After each option edit the per-data-set rows of the server's χ² split equal
the program's `chiSquared.out` on the edited text (six digits; the totals: T2 10896.788794139226 against the
program's 10896.7822, f19 1310.5367 against 1310.53998).  Undo gave the start χ² and objective exactly on all five.
The DW vertex is built with the session (6.3 s at every rebuild), so a vertex request costs no more than a plane-wave
one; the distortion factor costs about 17 ms per energy (Coulomb and optical alike: 3.2 s and 3.4 s for 200
energies), so its row has a bound of 1000 energies per request (about 20 s) instead of 20000.  No THM call needs a job.

**Where it departs from section 6.7 and the brief:**

1. *The keys' forms come from a table in `pyazr/ui/thm.py`*, not from the engine's schema: the engine has none for
   `<thm>` (finding 27).  The table's experiment keys are checked against the list the engine names when it refuses
   an unknown key; the global options are not (the engine's refusal of an unknown option names none).
2. *The diagnostics are rows of `/api/session/functions`* with a `tab` field; the Physics tab leaves out a row that
   names another tab.  The raw calls have routes of their own (`/api/session/thm_*`), one-to-one.
3. *The THM calls' answers are reshaped on the server* (fields with lines, facts, tables; `ANSWERS`), one function
   per call: the session returns nested dictionaries (finding 28).  The page has no code per function.
4. *The binding energy and the amplitude flag are read-only* (finding 11).  Data sets leave an experiment through
   their own `thm` field (the Data tab); the THM tab only adds members.
5. *The diagnostics run at the working x*; the Python text opens a session on the saved file and says that
   `params=x` gives another point.

**Extensibility as built** (before section 11.9, which removed `KEYS`, `THM_ROWS` and `ANSWERS`). A new key of
`<thm>`: one entry in `KEYS` of `pyazr/ui/thm.py` (value type of the
schema, choices, doc; the test fails for an experiment key the engine reads and the table lacks); until then the
page shows it as a text field.  A new THM diagnostic of the session: a row in `THM_ROWS` (arguments of the kinds
`ARGUMENT_KINDS` knows, returns with units and `lines`, facts) and its function in `ANSWERS`, an adapter operation
`Engine.<name>` and a route row; nothing in the page.  A new tab: a module `tabs/<name>.js` with `render` (and
`shown` if it is about part of a project) and one row of `TABS`; a session function drawn by it names the tab in
its row's `tab`.

### 11.9 The page on pyazr's tables (2026-11-14): THM keys, session functions, error parts, lists

The page's own copies of tables pyazr serves are gone.

- **The `<thm>` keys**: `pyazr/ui/thm.py` `keys()` is `pyazr.kinds()["thm"]` (the engine's table, `src/ThmKeys.cpp`,
  the one its parser reads); the hand-written table and the test that compared it with the engine's refusal are
  gone (`tests/pyazr/ui/thm_routes_test.py` checks that the route serves the table as it is and that a key added to
  it is served).
- **The session's functions**: the rows of `GET /api/session/functions` are pyazr's
  (`Session.describe()["functions"]`, `pyazr/session_functions.py`) in the shape the page draws (`physics.page_row`):
  every argument but the fit vector (kind `vector`: the server gives the working x), each marked `required` from the
  method's signature (`engine.required_arguments`: a row's `default` None does not tell an optional `J` from a
  required `pair`); the fields of shape `axis` but the axis's echo and the status field (`returns`); the top-level
  scalars and texts (`facts`).  A call goes to `Engine.session_function(fn, arguments, x)`, which calls the session's
  method by the row's name (with `params=x` where the row has `at_x`); `physics.answer_parts` walks the answer by the
  row: a field with `lines` gives one list per element of `lines.over`, labelled by the fields of `lines.label` with
  their values and units (`pair 2, jgroup 1, level 1, J 0.5, pi 1, E_level 0.609 MeV, Gamma 0.19047 MeV`; the page
  formats no J^pi).  What stays the page's (`PAGE` of `physics.py` and of `thm.py`): the tab, the group, the bound of
  one request (1000 energies for the distortion), the field that echoes the axis (`E` of the THM calls), the status
  field and its code (`open` → -2, `reached` → -3; else -1 where a field is not finite) with its words, the fields
  left off the first plot, and the two answers made into tables (`thm_experiments`, `thm_background`).
  `channel_functions` gains its argument `J` (a pair with a spin-orbit term) and `thm_vertex` its `strict`.
- **Errors**: `errors.parts_of` reads `e.what`, `e.allowed`, `e.fix`; a worker sends them with the class names and
  the message (`WorkerError.parts`); the body carries them, null where pyazr has none (the engine's own texts, as
  "no THM experiment 'none' in use"), and the page shows `message` then.  The messages are unchanged.
- **Lists**: the THM key forms offer the nuclide names of `kinds()["nuclides"]` in the beam, target and spectator
  boxes and the global sets of `kinds()["optical_potentials"]` in the optical-potential boxes (a free-text field
  whose grammar names `<nucleus>` or `<global set>`: `forms.js` `suggestionsOf`, `tabs/thm.js` `choicesOf`).  The
  masses of that nuclide list are not used for pairs (they are the THM code's nuclear masses).

**Measured** (`work/uiport/capture.py`, `compare.py`: the JSON of every session-function route and of every raw THM
route, before and after, on T1, T2 (LC09), two THM copies of `tests/18O_p_a_thm` (Coulomb distortion with the line
shape; coherent background), `examples/f19_pag_thm` and `examples/o18_lacognata2010` (DW vertex); 96 answers): the
raw routes are identical (36 of 36); in the 60 answers of the function routes every returned array, status, table
and fact value is identical, and the refusals have the same status and message.  What differs: the line labels
(above), the order of the vertex's facts (window before model), `norm` added to the background's facts, `J: null`
in the arguments of `channel_functions` and `strict: false` in those of `thm_vertex`, the order of
`channel_functions`'s keywords in the Python text, and the `units` of the two table answers (empty: the tables carry
their columns).

**Extensibility as built.** A new key of `<thm>`: a row of `src/ThmKeys.cpp`, nothing in the page.  A new
diagnostic of the session: its method and its row in `pyazr/session_functions.py`; the Physics tab lists and draws it
with no line of the page (an entry in `PAGE` only to put it on another tab or give it status words);
`routes_consistency_test.py` counts the rows of that table as served.

### 11.10 Pairs and THM columns in the Nucleus tab; what a fit returns in the Fit tab (2026-11-14)

**Nucleus tab** (`tabs/nucleus.js`): every edit is a `Project` call through a project route (revision checked, one
history entry, Undo by the server's history, the engine's words for a refusal).
- *Pairs*: a row opens `#/nucleus/pair/<key>`: its numbers as fields of `PATCH /api/project/pairs/{key}`, the fields
  being that route's arguments read from `GET /api/describe` (`routeArguments`), so `binding_energy` (B) is among
  them and a field the route gains appears with no line of the page (only its unit, from the page's `UNITS`, would
  be missing).  **Remove pair** sends `DELETE` (`remove_pair`); its refusal is shown as the engine gave it, and only
  then **Remove pair … and its channels** appears (`remove_channels: true`): the second, explicit action.  **+ pair**
  (`#/nucleus/pair/+new`): kind (particle, photon, beta), the partners as `add_pair`'s dicts `{J, parity, Z, mass,
  g}` (none for a photon's light partner), the route's other arguments as text boxes; numbers are sent as numbers,
  any other text as typed (pyazr refuses it in its words).  The THM nuclide list of `kinds()` is not offered here:
  its masses are nuclear masses, the pairs' are atomic.
- *Channels*: the amplitude flag is a check box per channel row (`PATCH .../channels/{pair}/{l}/{s}`,
  `width_is_amplitude`; the width is not converted, said in the column's doc); a photon channel's refusal is the
  engine's.
- *New level* (`#/nucleus/level/+new`): J (a number), parity, E, energy fixed, `l max`; **Channels** asks `GET
  /api/project/allowed_channels` (a new route, `Project.allowed_channels`; its row left `PROJECT_NOT_SERVED`) and
  lists them checked; **Add level** posts the checked ones (width 0, fixed), or `"allowed"` when every one is
  checked so that a J^π the project has takes its group's channels (`add_level`'s rule).
- The THM tab's Channels section stays read-only and links each row to the pair or level in the Nucleus tab.

**Fit tab**: an engine fit job's result carries `outcome`, the engine's `notices` (no longer made by the server from
the log) and, after MINOS, `sides` per parameter (`minos_valid` of pyazr: `{valid, at_limit, call_limit,
new_minimum}` per side) with `valid` both sides; the fallback fitter's name is the fitter table's `fallback` column
(the server's constant `FALLBACK_FITTER` is gone).  The page shows the outcome in a box (`ok-box` converged,
`warn-box` otherwise, with "not converged: the best point is where the fitter stopped"), a warning when a MINOS
interval is not valid, and the column **MINOS validity** (`fitting.js` `minosValidity`: the flags set per side, one
word for iminuit's).  The fitter rows' facts add "honours a call limit" and "falls back to" from
`kinds()["fitters"]`, shown under the method and in the job table.

**Measured.**  T2, engine MIGRAD over HTTP (`work/uiport/fitlimit.py`): with `max_calls` 40 the job ends `ended`,
outcome `call_limit`, best 7087.924186759944, the engine's notice that the end point is the lowest objective
evaluated; without a limit `converged` at 297.8603087301522; headless Chrome draws the first in a `warn-box` and
the second in an `ok-box`.  Routes against the direct `Project` calls (`tests/pyazr/ui/nucleus_routes_test.py`):
equal texts for add, remove, remove with channels, B, the amplitude flag and a new level; an unused pair leaves T1's
objective exactly; T2 with B 2.224566 → 2.374566 MeV equals a session on the direct edit (11039.423103526951).  In
the self-test (`#selftest?edit=1`) a pair is typed into the form, added, removed and both undone with the buttons on
T1 and T2: the χ² unchanged throughout.

**Extensibility.**  A new argument of `set_pair` or `add_pair`: its route row's `kwargs` (one word) and nothing in
the page (a unit in `UNITS` if it has one).  A new column of the fitter table: a fact in `engine_rows` (one entry);
a new outcome word or notice: nothing.

### 11.11 The page for the median user (2026-11-14): the original's face, two depths, Calculate, run settings

The PI after using the page: a user who knows R-matrix and the old program, but not the project's grammar, must not
be overwhelmed, should set pairs, levels and data easily, and check an evaluation easily; the experts work through
pyazr. Measured first (the typed places, the six tasks in steps: the HANDOFF entry of the branch
`feature/upgrade-9-webui-first-view`).

- **The original's face** (section 7): the icon, colours from it, symbols beside every state colour, About.
- **Two depths, one page.** The first view (`layout.js` `FIRST_VIEW`, the default; `ctx.expert` false) shows the
  old program's path and words: Project, Particle Pairs, Levels and Channels, Segments, Experimental Effects, Fit,
  Plot. Each is a view of one of `TABS` (the Nucleus tab's `pair` part or the rest; the Data tab's `targets` group
  or every other group); a route keeps its meaning in both depths and marks the view it belongs to (`tabOfRoute`,
  after the tab has drawn: the Data tab says which group it drew). Physics and THM are in the expert's depth. In the
  first view the Nucleus tab's γ decays and potentials, the Project tab's session options, history, text and engine,
  and the Results tab's totals and band are under a `details.more`; the expert's depth (`expert` in the header, key
  `e`, kept in `localStorage`, `?depth=` for one view) shows what the page showed before.
- **Rare fields: a column of the engine's table.** `FieldDef::rare` (`experiment/src/Schema.cpp`, `Rare(...)`),
  written as `"rare": true` by the schema: 35 fields of 6 kinds (provenance texts, pair keys, THM, reference, UPOS,
  multipole, shape constant; a target's phase, compound correction, isotope sum, Doppler, screening; a potential's
  pair key, matching radius, Numerov step; a γ decay's mixing and sign). A first-view form (`foldRare`) puts a rare
  field under "More fields" unless the record sets it, needs it or has a message about it (`forms.js` `folded`);
  an edit of a folded field is found as any other (`changedFields` reads both parts). The field's doc is written
  under its title in the first view (`showMeaning`), on hover in the expert's.
- **Calculate** (header, every tab, key `c`): `POST /api/evaluate` of the current text, remembered as a snapshot
  (revision, χ² and objective, χ² and N per data set: `calculate.js`, page state lost on reload), then the Plot tab
  (`#/results?calculated=1`) with "Calculate: what changed": Δχ² in total and of the objective, one row per data
  set (before, now, Δ; the largest change first), each row scrolling to its plot, the plots in that order. The
  baseline is the newest calculation that is not the current one (same revision and χ²); the status line gives the
  change since the last Calculate after every edit. Clicks from an edited width to the plot of the data set it moves:
  one (Calculate), the set that moved most first; before, one (Results) and the data set found by eye.
- **Run settings** (Project tab, both depths; `tabs/settings.js`): the rows of `kinds()["config"]`, a record kind
  made by `pyazr/described.py` from the engine's table (`config_table()` of the binding: `io::ConfigPositionals()`
  and `io::ConfigKeys()`, which gained the columns `title` (the old program's words where it had the setting),
  `group`, `absent` (what the engine does without the key) and `rare`; the positional entries their `values`, now
  also the edit's check). A switch for a 0/1 key, a choice for a positional entry, a number box with its bound, a
  text box for a path or a directory; each change is one `PATCH /api/project/config` (one history entry); `×`
  removes a key; a search over key, title, meaning and group; rare rows under "More settings" unless set; the free
  key and value table under "Advanced". A new key of the table appears with no line of the page.
- **Next step** (Project tab, first view): the first missing thing in the old program's order (`nextStep`: pairs,
  levels, data (records of a kind with a quantity), the engine's refusal, else Calculate) with its button. Empty
  pairs, levels and Data groups say what belongs there and offer the first one.

**Tests.** `tests/pyazr/ui/js/firstview.test.mjs` (the views and their marking, the groups of the Data views, the
folding on every kind with rare fields of the engine's schema, an edit of a folded field, the next step, the
comparison of two calculations); the self-test visits every view in the expert's depth, then each first-view tab,
a data set's form (one "More fields") and Calculate twice; with `edit=1` it presses Calculate after the width edit
and `smoke_test.py` checks that the comparison is the engine's two χ² (before, after) and that the Plot tab is
marked.

### 11.12 The covariance builder and the reaction rate (2026-11-14)

Both are arithmetic of pyazr with a screen of the page that calls it; the page computes no matrix, rate, ratio or
window.

- **Covariance builder** (`pyazr/covariance.py`, user guide `covariance_builder.rst`).  The rows of a data set's
  file come from the engine's reader (`_azure2.data_rows`, `api/src/pybind_data_rows.cpp`: every data row in the
  engine's unit, its record number, whether the window keeps it, the normalization's fit flag and prior).  Parts
  (`PART_KINDS`: own, scale, range, length, correlation) become the engine's terms (diagonal, vector, relvector,
  block); the checks are the engine's (own variance positive definite on the fitted rows, sqrt(C_ii) against the
  column at 1 %), plus a refusal of a part that is not positive semi-definite.  The parts are kept on the file's
  `#pyazr.covariance` line.  Routes `GET /api/covariance/{name}`, `POST /api/covariance/build`, `POST
  /api/covariance/write` (server methods over `pyazr/ui/utilities.py`); the write is `State.edit_with_files`: the
  file and the field in one history entry whose undo restores the file as it was (removed when new) and whose redo
  writes it again.  The screen is `tabs/covariance.js` (a tab without a header link, `#/covariance/<name>`), opened
  by "Covariance…" in `tabs/data.js` for a record whose kind has a `covariance` field in the schema.
- **Reaction rate** (`pyazr/rate.py`, user guide `reaction_rate.rst`).  A rate session is the binding's Session made
  with `rate_entrance=`, `rate_exit=` (set up as the program's `rate` command; the API now checks the
  pairs with `CNuc::RatePairsProblem`, so exit -1 works), and `Session.reaction_rate(params, min, max, step, file)`
  runs `ReactionRate` as `RunReactionRate` does (`ReactionRate::Rates()` added in `include/ReactionRate.h` to read the
  numbers; `src/ReactionRate.cpp` unchanged).  Named grids `TEMPERATURE_GRIDS` (a new grid: one row), the band
  (`rate_band`: linear in ln rate with `bands.step_sizes`, or 16/50/84 % of draws), `gamow_window` (from
  `physics.sommerfeld_eta`), CSV with header, a comparison table.  The page: the job `rate`/`engine` (`utilities.py
  rate_engine`, a row of `jobs.py`, setting type `text` added with its widget in `fitjob.js`), the tab
  `tabs/rate.js` (`#/rate`), route
  `GET /api/rate/grids`.

**Measured** (`tests/pyazr/covariance_builder_test.py`, `reaction_rate_test.py`, `ui/utilities_test.py`): the
engine's chi-squared of p12c_60 with a built file (own, 5 % scale, 3 % over 0.2 MeV) 9.38942514115 equals v^T C^-1
v; pyazr's rates equal the program's `reactionrates.out` on a grid, a list, with a parameter file, and T1's and T3's
golden files; the bands equal direct loops; the routes' answers equal the direct calls; undo removes the file.

**Rate libraries** (`pyazr/rate_formats.py`, a second commit): the rate job also writes the rate as a STARLIB
entry (on `grid="starlib"` with a band only) and as REACLIB sets (fitted on 0.01-10 GK within a tolerance), each
the file's text or the refusal's three parts in `result.formats`; the Rate tab has the labels and flags as controls
and the two downloads beside the CSV.

**Extensibility.**  A new kind of covariance part: a row of `PART_KINDS` and a branch in `_terms` returning its
terms (the page's form comes from the row).  A new named temperature grid: a row of `TEMPERATURE_GRIDS` (the page
lists it from `GET /api/rate/grids`).

### 11.13 The ENDF tab: an evaluation over the data, the check of a file (2026-11-14)

The PI: an easy way of checking the evaluations, and of checking that an ENDF file is correctly formatted with the
ENDF tools.  The tab is named **ENDF** (id `endf`): the word an evaluator looks for, and both screens start from an
ENDF-6 file.  It is a row of `TABS` (expert's depth); a first-view entry is one row of `FIRST_VIEW`
(`{ title: "ENDF", tab: "endf", hash: "#/endf" }`), left to the merge.

- **Over the data** (`#/endf?files=a.endf,d/b.endf`): the files of `GET /api/endf/files` with a box and a colour
  each (and a typed path inside the root), what each holds (`POST /api/endf/summary` = `read(path).summary()`),
  and per data set the plot of the data times its norm, the model and each file's values at the points, the ratio
  plot below (file / data with the data's error, file / model), the sections and coverage under it; where a file
  has nothing for a data set, its refusal (what, allowed, fix) in that place (`POST /api/endf/curves`, a session
  operation `Engine.endf_curves`: `Material.curves(data, session=s)` for every data set at the working x).
- **Where the shown files live** (section 5): in the route (`?files=`, `?file=`), not in the project text and not
  in the server; the tab remembers the last list per project for the page's life (`ctx.endfShown`).
- **Check and correct** (`#/endf/check?file=a.endf[&job=id]`): the codes of `GET /api/endf/checkers` (found or not,
  what each does and does to numbers, how they are found); **Check** is the job `endf_check`/`codes`
  (`pyazr.endf.check`): the own findings (line, MF, MT, message, page) apart from each code's (MF, MT, record, text,
  `meaning`), the counts, the raw report in a `details`; **Correct** is the job `endf_standardize`/`program`
  (`pyazr.endf.standardize`), the text of `PROGRAMS[program][3]` said before the run, `Candidate.diff()` as a table,
  then **Keep as a new file** / **Replace** (asks a second time) through `POST /api/endf/accept`.
- **Jobs without a session**: `JobRow.session = False` (the worker builds no session; x0 and names are empty) and
  `JobRow.files` (settings naming a file: the server refuses a path outside the root, 403, before the job starts).
- **Files**: `pyazr/ui/evaluations.py` (the listing `endf_files`, the only place that says what an ENDF-6 file is:
  its first two lines read by `pyazr.endf.records.split_line`; the listing is a filter over the page's one file
  listing, `Server.files` of `GET /api/files`, since the merge; its `job_rows()`, the two job rows, joined to the
  table by `jobs.job_rows`), `tabs/endf.js`, five rows of `routes.py`, `Engine.endf_curves`, server methods `endf_*`; lines in shared files: `app.js` (one import, one `TABS` row), `style.css` (a block at the end),
  `selftest.js` (the ENDF views and Check pressed), `smoke_test.py` (two views, the ENDF files and the stand-in codes
  in the server's environment).
- **Extensibility**: a new checking code is a row of `PROGRAMS` (with its numbers text) and its deck in
  `checkers.py`; a new quantity of `curves` is a branch of `pyazr.endf.curves`; neither changes a file of the page.

**Tests**: `tests/pyazr/ui/endf_test.py` (the routes and jobs against the direct calls on T1 with the two test
files and the stand-in codes; refusals; accept without and with overwrite), `js/endf.test.mjs` (the plots and rows
from an answer), the self-test and two screenshots (`expert-endf`, `expert-endf-check`).

### 11.14 The lists in words, the Reaction Rate in the first view, nothing with a structure typed (2026-11-14)

**The PI's decisions of 2026-11-14** (answers to the director's questions): the Reaction Rate tab is in the first
view as its last tab; the Fit tab has two depths, the sampler, the bands and MINOS one click away in the first view;
the first view is the default also for a returning user (the page remembers the expert's depth only when the user
switched to it with `expert` or the key `e`, `localStorage` `azure2-depth`; nothing else stores it); when the project
has THM content the Physics and THM tabs join the first view (after Plot, before Reaction Rate).  The REACLIB
defaults (tolerance 5 %, set label azr2) stay as they are.

- **Lists in words** (`lists.js`).  A group of the Data tab shows one list per record kind under the kind's title
  (the schema's), the data sets first (`listKinds`: kinds with a `quantity` field), the records they name after.  A
  kind the engine describes shows the entries of its description (`GET /api/project/describe`) that pyazr's table
  marks for the list: `pyazr.project.DESCRIPTION_FIELDS`, served as `kinds()["description"]` (reaction, observable,
  data file, energy range (MeV, lab), angle window, points, active), the normalization in use with fitted or fixed
  and the shared record (from the description's parameter of kind `norm`) and χ²; any other kind shows the fields
  its records set (not those the schema marks rare), headed by the schema's `title`, the key in the hint.  "Used
  by" is the engine's references (`Project.references`, `GET /api/project/references`).  The pairs' list is the old
  program's columns from `ProjectPair.FIELDS` (`kinds()["pairs"]`: title, unit, doc, `list`, `rare`; the particles
  by name split from the engine's pair label, J^π and the kind in words, new fields of `ProjectPair`).  The levels'
  list shows J^π and the energy once.  A new record kind costs nothing: its list comes from the schema.
- **The model as a curve** (not built; measured).  The Plot tab draws the model at the data points.  What exists:
  `bands.trimmed_model(azr, grids=[extrapolation fields])` writes a `.azr` with extrapolation records only, and a
  `Session(..., extrapolations=True)` on it evaluates them (`_extrapolation_vector` carries x); on T2 one excitation
  function of 601 points costs 0.02 s (trim) + 0.01 s (session) + 0.006 s (evaluate).  It is not one call: the
  extrapolation record of a data set (reaction or sum or ratio, quantity, engine unit, energy range from the
  description, the angle of a differential set, J, l, s of a phase shift, the target, beam and detector of a folded
  set) must be made from the data set, the temporary file lives in the project's folder, and the session must be
  kept per text.  Needed: a pyazr call `Session.model_curve(name, points=, params=)` (or `pyazr.curves`) returning
  the model on a fine grid over a data set that is a function of one variable, with the rule for the others
  (several angles, a ratio of two pairs), and its route `GET /api/evaluation/curve/{name}`.  Not built: the brief
  stops this item where a new pyazr call is needed.  (Built since: section 11.16.)
- **First view**: `FIRST_VIEW` ends with Physics and THM (`when: "thm"`, shown by app.js `CONDITIONS.thm`, the THM
  tab's own rule `thm.js hasThm`: the server's `<thm>` or a data set naming an experiment) and Reaction Rate (the
  pairs, the temperatures and Calculate; the parameters, the band, the comparison and the library labels under one
  `details.more`).  Keys `1`–`9` follow the tab bar as shown.
- **Files** (`GET /api/files`, `Server.files`: the project directory to 4 levels, at most 2000 files, not
  `output/`, `checks/` or hidden ones, never outside the root; kind at a glance: project, table with its columns,
  text, binary).  Every field of type `path` offers them (`forms.js` row `path`, `ctx.files`), as do the covariance
  builder's file and the Rate tab's comparison (`ui.js withFiles`).  A file chosen in a field is previewed under it:
  its first rows and what the engine makes of the record with it on a copy of the project (`GET
  /api/project/records/{name}/preview?field=&value=`, `Server.field_preview`: `Project.set` then
  `Project.describe`; a refusal with its three parts).  `?field=file&preview=<path>` links to that state.
- **References**: a reference's choice lists the existing records and "new <kind> …" (`forms.js NEW_RECORD`), which
  opens the form of a new record of that kind (`#/data/+new?kind=`); never a value sent.
- **Lists of numbers** (`number_list`, a detector's attenuation): a box per number and two empty ones; an emptied
  box removes its number.  **A job's free parameters** (MINOS's `parameters`, setting type `names`): a check box per
  free parameter of the working x (`ctx.freeNames` from `GET /api/x`), none ticked meaning every one.
- **The pair form** from `kinds()["pairs"]`: `partner_fields` (J, parity as a choice of its `values`, Z, M (u),
  g-factor), the route's other arguments with the words and units of `fields` or `arguments`; the page's lists
  `UNITS`, `PARTNER`, `PAIR_KINDS` are gone.
- **The Fit tab's two depths**: the first view's job form offers the kinds `jobs.py FIRST_VIEW_KINDS` marks (a fit;
  the row's `first_view`), "More: minos, sample, bands, rate …" offers the others; the job table under the form is in
  the expert's depth.  `#/fit?job=<kind>&method=<method>` opens a job's form.
- **Still typed**: a data set's `columns`, a `reaction`, a `grid` (extrapolation energies and angles), a target's
  `composition` and `nucleus`, the structured fields (stopping, straggling, profile, isotopes) and a THM key's
  name; the covariance builder's list of data sets.  A list of numbers grows by two boxes per Apply (repaired:
  section 11.16, a box follows as the last one is typed in).

**Tests.** `tests/pyazr/ui/js/lists.test.mjs` (heads per kind, the description's columns, the norm's words, used
by, the pairs' columns, "new …", a list of numbers), `nucleus.test.mjs` (the pair form's table),
`firstview.test.mjs` (ten views, `when`), `nucleus_routes_test.py` (pairs' names, references, a target's
description refused), `static_test.py` (the files route, the preview against `Project.set` + `describe`), the
smoke test (the tabs shown on T1 and T2, the preview and "new …" pressed by the self-test, the Fit tab's kinds in
both depths; views first-rate, first-physics, first-dataset, first-file-preview, first-new-detector, first-fit,
fit-minos, first-new-pair).

### 11.15 The prior builder (2026-11-14)

The prior on the fit parameters (`priorCovariance`) from named parts, as the covariance builder is for a data set:
pyazr's arithmetic with a screen that calls it (user guide `prior_builder.rst`).

- **pyazr** (`pyazr/priors.py`): `build(parameters, parts, directory=)` over a session's parameter table (names,
  kinds, units, fixed marks, limits, the priors of their own, the basis come from `Session.parameters`; nothing is
  listed twice); parts `PART_KINDS`: `fit` (an earlier fit's `param.sav` and `covariance.dat`, the numbers of
  `prior_from_fit`, the tag checked against the run's basis), `width`, `correlation`, `chain` (a sampler job's
  `chain_<id>.npy` with its JSON: `read_chain`; the JSON's `parametrization`, written by the sampler job since the
  PI's decision of 2026-11-14 with the values of `BASES`, is read by `chain_basis` and checked as a fit's tag is:
  another parametrization refused, none recorded taken as the run's with a note; the page asked nothing before
  and showed the assumption as a note, it now shows the recorded parametrization in the same place).  The `Prior` record prints
  itself; its checks are the engine's (symmetric, the engine's Cholesky rule) with the eigenvalues and the condition
  number; its notes say what the engine does: a fixed parameter marginalized, a prior of its own replaced (measured:
  `ReplaceDiagonalPriors`), a centre outside the limits.  `read_prior` reads the grammar as `PriorCovariance::Read`
  does (a table of 35 files, `tests/pyazr/prior_builder_test.py`).
- **Server** (`pyazr/ui/prior.py`, three methods of `server.py`, three rows of `routes.py`): `GET /api/prior` (the
  key, the parts from the file's `#pyazr.prior` line, the proposed name, `kinds_table()`, the project's files: the
  page's one listing, `Server.files(outputs=True)`, which adds the program's `output`, where a fit's `param.sav` and
  `covariance.dat` and a sampler's chain are, since the merge), `POST /api/prior/build` (the session's
  parameter rows made back into a `pyazr.Parameters` with `session.basis_of(options)`, then `build`), `POST
  /api/prior/write` (`State.edit_with_files`: the file and `set_config(priorCovariance=)` as one entry).  Removing is
  the config route.  `State.edit_with_files` now keeps an entry when only a file changed (a rewrite under the name
  the text already holds); before, such a write was undone at once and the answer said it was written (the
  covariance builder had the same defect).
- **Page** (`tabs/prior.js`, `#/prior`, no header link): the parameter rows of `fitting.js` with check boxes; a part
  per row with a control per field type (`path:<suffixes>` a select of the project's files, `names` the checked
  parameters, `number`, `integer`, `pairs` two selects and a coefficient, `matrix` shown only); Build; the heat map
  of the server's correlation; Write and use; Remove the prior.  Shared files, one line each: `app.js` (the import,
  the `TABS` row), `tabs/fit.js` (the link "Prior on the parameters…"), `smoke_test.py` (`VIEWS`).

**Tests.** `tests/pyazr/prior_builder_test.py` (pyazr: the readers' table, the term through a session = d^T C^-1 d,
the fit part = `prior_from_fit`, the chain part, every refusal and note), `tests/pyazr/ui/prior_builder_test.py`
(the routes against the direct calls on T1, the history), the smoke test's view `prior`.

**Extensibility.** A new kind of part: a row of `PART_KINDS` and a function of `_BLOCKS`; a new type of field: its
control in `fieldInput` of `tabs/prior.js`.

### 11.16 The model as a curve; a project from an evaluation and back; the first view's rough edges (2026-11-14)

- **The model as a curve** (`pyazr/model_curve.py`; `Session.model_curve(name, x=None, points=200)` forwards to
  it).  The argument is `x` because `Session.evaluate(x)`, `transform(x)` and `write_outputs(x)` take the fit vector
  so; `points` is the number of grid values per branch.  For each data set an extrapolation record is made from its
  record (the fields the engine's extrapolation kind takes, `kinds()["records"]["extrapolation"]["fields"]`: the
  reaction, quantity, unit, pair keys, target, beam, detector, quadrature, J, l, s, ...) with a grid of `points`
  values between its first and last point and its own points; one record per branch (a fixed angle of an
  excitation function, in the frame of the file's angle column where the angle is one value in both frames; a fixed
  lab energy of an angular distribution, the angles c.m.).  Energies are the points' lab energies after the energy
  shift at x (`extra["energy_lab"]`).  `bands._trimmed_project` makes the project of these records alone; a session
  of the extrapolations (`CurveSession`) evaluates them at the data session's x carried by name
  (`carry_model_parameters`).  `model_curves(session, names)` makes several in one session: their folds are then
  planned on one mesh, as the data session's.  Measured (`tests/pyazr/model_curve_test.py`): at the data points the
  curve equals `Session.evaluate` to rounding for every unfolded data set of T1, T2 and an alpha + 12C ratio to
  Rutherford, to 4e-9 for T2's folded data sets without a shift, to 8e-5 for those with a shift and THM (the data
  session folds them on grids of their own); one data set alone (`Session.model_curve`) on its own mesh: 2.4e-3 on
  T2.  Between the points: the curve's grid read as a data file gives the curve to rounding (four data sets).
  Cost: T2 0.6 s for all thirteen curves, T1 3.9 s, 13N 14.9 s (its external-capture integrals), then one
  evaluation per x (13N: 0.015 s).  Refused with the three parts: a column that gives each point its own value
  (T2's beam-profile windows: `CURVE_COLUMNS`), `columns from_file`, points not one function of one variable, more
  than 64 branches or 20000 points, a record the engine refuses on a grid (an analyzing power at a lab angle: the
  engine's words; the other data sets are made).
- **Route and page**: `GET /api/evaluation/curve/{name}` (`Engine.curve`: the curves of every data set in one
  `CurveSession`, kept while its records are the same, evaluated once per x).  The Plot tab draws the points at
  once and the curves when they come (one line per branch, the residuals at the points; a refused data set keeps
  its line through the points with the refusal under the plots).  The ENDF tab draws the model and each file on
  the curve's grid (`POST /api/endf/curves` adds `curve` and per file `on_curve`: `Material.curves` takes a
  `ModelCurve` as it takes a data set); the ratios stay at the points.
- **A new project from an ENDF-6 file** (the ENDF tab, screen `#/endf/new?file=`, and "New project…" beside each
  file): what MF=2 holds (`POST /api/endf/resonances`: `Material.resonances(isotope, energy_range)`, nothing
  written, or the refusal); Convert, the job `endf_import`/`to_project` (`pyazr.endf.to_project` into a directory
  of the job's own: the separation energy and its source, the masses against AME2020, the pairs, the levels with
  their Brune iteration, what was not written, the notes, the project's text); "Write the new project"
  (`POST /api/endf/keep_project`: a new `.azr` of the folder, never an existing file); then "Open … here"
  (`POST /api/project/open`; disabled while the open project has unsaved edits) or `pyazr ui <file>` for a page of
  its own.  The ENDF tab, because both ways start or end with an ENDF-6 file and the tab lists the folder's files.
- **Export** (screen `#/endf/export`): before the run, `GET /api/endf/export` (`pyazr.endf.from_session(s, x)`: the
  losses and the photon widths, or the refusal; the data's energy range; a file name not in use); the job
  `endf_write`/`writer` (`pyazr.endf.write` on a session of the job's text under the project's name, so that the
  MF=1 text names the project and not the server's working file; settings: the MF=1 fields, the range, the
  tolerance, the number form); the result: the sections, each MF=3 grid, the MF=2 check, the losses, what is not
  written, and "Check…" on the file, which the tab's file list now shows.
- **The first view's rough edges.** (1) A finished job's follow-up buttons follow the job form's rule
  (`fitting.js followersShown`: the kinds of `jobs.py FIRST_VIEW_KINDS`, the others behind "More: …", which opens
  them all as the form's "More" does).  (2) A tab lists the jobs of its own rows (`fitting.js jobsOfTab`, by
  `JobRow.tab`); the Fit tab says how many jobs other tabs started and lists every job with "Every job".  (3) The
  first view's lists of data sets and pairs hold the old program's columns: one table,
  `pyazr/ui/first_view.py FIRST_VIEW_COLUMNS`, served as `first_view` of `GET /api/describe`; the others stay in
  the expert's depth and the record's form (the PI's decision of 2026-11-14).  (4) A finding of no section reads "whole file".  (5) A list of numbers
  shows one empty box, and the form adds the next as soon as it is typed in (`grows`).  (6) Keys `0` (the tenth tab)
  and `[` `]` (the previous and the next) reach every tab.  (7) The smoke test's header lines are within 120 columns.

**Tests.** `tests/pyazr/model_curve_test.py` (the curve against `evaluate` and against the same points read as
data; shifts; the grid; the refusals; an injected frame error, the energies written as c.m., failed 29 checks),
`tests/pyazr/ui/curves_test.py` (the route against `model_curves` bit for bit on T1 and T2, at a moved x; the ENDF
file on the curve's grid against `Material.curves`), `tests/pyazr/ui/endf_project_test.py` (resonances, the import
job, keep and its refusals, opening the kept project; the export's preview against `from_session`; the written file
byte for byte the direct `write`), node tests (`results`, `endf`, `fitting`, `firstview`, `lists`), the self-test
(Convert pressed) and the smoke test's views `expert-endf-new`, `first-endf-export`.

## 12. What pyazr 3.0's design lacks for the UI (findings)

None of these changes the decided objects; each is a call or a field to add:

1. **A session from the project text without saving the user's file.** `Project.session()` saves the project first
   (`project.py`), so every evaluation of an edit would overwrite the user's file. The spike writes a working file
   beside the project and opens a session on it (today's `azure2(file, cwd=...)` allows the working file to live
   elsewhere). Wanted: `Project.session(**options)` on the text as it stands, the project directory as base, the
   file untouched.
2. **Progress and stop of long engine calls**: `fit_with_engine` (decision 8), MINOS, the band, a session build that
   computes external-capture integrals (13–15 s), the Coulomb table build. A stop request from another thread ending
   cleanly with the best point (today `Config::stopFlag` makes the objective return 0; not bound), and progress as a
   polled record or callback (the call counter `CountCallsInto` exists, unbound). The build should release the GIL.
   *Delivered for the fit* (`fit.md`): `Session.fit_with_engine(method, minos, x0, max_calls, progress, stop)` runs
   without the GIL; `stop` is a `pyazr.session_fit.StopRequest` (or a `threading.Event`) set from any thread, and the
   fit ends at the best point seen with status `stopped`, the session usable; `progress` is a callback per record.
   MINOS (inside the call) obeys the same request. Not delivered: progress and stop of the band, of a session build
   and of the Coulomb table build.
   *Closed for the page* (section 11.6): the Fit tab runs the engine's fitters and MINOS as jobs with progress and
   Stop.
3. **A record of the engine's progress lines** (`Iteration: N Chi-Squared: X`) as data, not text on standard output.
   *Delivered*: `FitProgress` (fitter, kind evaluation or iteration, evaluation number, iteration, calls, objective,
   best, the fitter's measure, x); the program prints its lines from the same records.
4. **The engine's notes as a per-session log**, not standard output (2.5): `pyazr3.md` 2.5 turns notices into
   `EngineWarning`; the many informational notes (covariance read, groups, THM kinematics, mesh re-planning) need a
   log the session hands to its caller. *Partly delivered (pyazr's errors and self-description)*: the notices a
   session holds back are `EngineWarning`s, and the server passes them to the client; the informational notes still
   go to standard output (the workers' logs).
5. **χ² per data set defined as the program's**: `Evaluation.data[i].chi2` should equal the row of `chiSquared.out`
   (s02/s03 of T2 differ under today's `segment_chi2`, 2.2).
6. **The reaction rate**: no pyazr call is designed (menu 5 of the program only; `pyazr3.md` deletes the Python rate
   integral of `examples/reaction_rate.py`). Wanted: `Session.reaction_rate(entrance, exit, temperatures, x)`.
7. **The inverse transformation**: physical values come from `Session.transform(x)`; setting a physical value of a
   free parameter needs the way back (a basis-aware `Parameters.put`), otherwise it is a project edit and a rebuild.
8. **A parameter vector written into the project text in memory**: `save_fit(path, x)` writes files; the UI needs
   the text with the values of `x` (accept a fit, release a slider) as an edit of `Project`, giving a history entry.
9. **The three parts of an error as fields** (`what`, `allowed`, `fix`) of the exception, not only its text.
   *Partly delivered*: pyazr's error classes (`pyazr.errors`, each also the builtin it was) and messages that say the
   three parts; the parts are not attributes, so the server keeps the whole text in `what`.
   *Closed* (section 11.9): the error body's `what`, `allowed`, `fix` are these attributes, carried from the worker.
   *Delivered in pyazr*: `e.what`, `e.allowed`, `e.fix` on every `pyazr.Error` (None where the message is the
   engine's own text, as an `AZURE2Error` from C++); set where the message is built from the parts (the session's
   argument checks, a closed session, a missing working directory, the physics functions' domains, unknown records,
   parameters, levels, kinds, THM experiments) and, for a `FormatError`, from its first diagnostic.  The server's
   error handler (`pyazr/ui/`) can read them instead of putting the whole text in `what`.
10. **Pairs added and removed, with channel generation** (*closed* for the page: the Nucleus tab, section 11.10):
    `add_pair`, `remove_pair`, and the channels a J^π allows from
    the pairs and `l_max` (Qt's `LevelsTab` rules; `add_level` takes channels as given).
    *Delivered in pyazr*: `Project.add_pair(kind=, light=, heavy=, separation_energy=, radius=, ...)` (returns the key;
    a pair no level names is held by a pair row, inactive with level ID 0), `Project.remove_pair(key,
    remove_channels=False)` (refused while a level has a channel of it, naming them), `Project.allowed_channels(J,
    parity, l_max=2, multipolarity_max=2, multipolarities_max=2, pairs=None)` (the rule of `LevelsTab::calculateChannels`,
    in the engine: `experiment::AllowedChannels`), `add_level(..., "allowed")`.  The Nucleus tab can call them.
11. **THM binding energy B and the width-is-amplitude flag** per channel (*closed* for the page: the Nucleus tab,
    section 11.10): `set_pair` and `set_channel` have neither.
    *The page shows them read-only* (section 11.8): the session's pair records carry `binding_energy` and its
    parameter rows `width_is_amplitude`.  Wanted: `set_pair(key, binding_energy=...)` and `set_channel(...,
    width_is_amplitude=...)` (or the record fields they are), so that the THM tab edits them as project edits.
    *Delivered in pyazr*: `Project.set_pair(key, binding_energy=B)` (a particle pair, B >= 0, every row of the pair)
    and `Project.set_channel(level, pair, l, s, width=, width_is_amplitude=)` (a particle channel); without a session,
    `Project.pairs()[k].binding_energy` and `Project.level(label).channel(pair, l, s).amplitude` read them.
12. **Outputs of a physics function described in its row** (names, units, which argument is the natural axis), so
    the Physics tab can draw it without a per-function plot.  **Closed**: every row of `/api/physics` has its
    arguments (unit, kind, domain, bounds, default, choices), its returned fields with units and `axis`; the two
    routes are served.  The Physics tab still lacks: the level matrix (no engine call), the transfer functions on
    potentials, structured rows for stopping (`kinds()["stopping"]` gives signature and doc only), the reaction rate,
    and the page code itself (a later step).

Delivered by pyazr's errors and self-description besides: `pyazr.kinds()` (served at `/api/kinds`; its rows of
session options are the server's options table, its physics rows the list at `/api/physics`), `pyazr.describe()`
(in `/api/describe`), `azure2.describe()` (the adapter's parameter rows, with units). Still missing for the server:
findings 1–3, 5–8, 10–12 as written, the three parts as attributes (9), the informational notes as a log (4).

Found while building the page (step 2):

13. **A schema of `<config>`**: its keys with type, allowed values, default and doc (as the record kinds have), so the
    Project tab can draw a generated form instead of key and text.
14. **The order of each point** of a Legendre data set in the data arrays (`Evaluation.data[i]`), so the plot draws
    one trace per order as the quantity row says (`columns.order`).
15. **The convention of each channel's width** in `Project.levels()` (partial width in eV, ANC in fm^-1/2, amplitude
    in MeV^1/2), so the Nucleus tab can put the unit beside each width; and the channel's Wigner limit and physical
    value next to it without a session (today only `/api/parameters` after a build).
16. **`Project.describe(name)` for every record kind**: it refuses the kinds that compute nothing (`gamma_decay`:
    "the engine describes the records that compute something"); the page shows that sentence where a description
    would be.

Found while building the Fit tab:

17. **The owner of each parameter as a record name or level label and field** (or `Project.set_parameter(name,
    value)`), so the table edits a fixed value as a project edit; the parameter rows carry the engine's owner
    numbers (J group, level index, segment key) only.
18. **The unit of each parameter's physical value** (`physical_value` of a width is in eV, fm^-1/2 or MeV^1/2 by the
    channel's convention; the row's `unit` is that of x).
19. **The prior of every nuisance in its parameter row** (normalizations, energy shifts, shift-function terms have
    `prior_nominal` and `prior_sigma` only in their records).
20. **Progress, stop and the result's covariance of the engine's fitters** (finding 2 again): the Fit tab's rows for
    the engine's MIGRAD, LM and MINOS wait for `fit_with_engine`; the sampler and bands for `Objective`.
    *Delivered in pyazr* (`fit.md`): the rows are `pyazr.kinds()["fitters"]`; the result has `x`, `objective`,
    `status`, `converged`, `edm`, `calls`, `covariance`, `errors`, `minos`. What the server should call (not done
    here, `pyazr/ui/` is not touched): a job `fit_engine` in `pyazr/ui/jobs.py` that runs
    `session.fit_with_engine(method, minos=..., x0=x, progress=<post the record to the job's progress>,
    stop=<the job's threading.Event>)` in the worker, offers the methods of `kinds()["fitters"]` (with `minos`
    only where the row allows it), and on success hands `r.x` to the history as a new revision and `r.covariance`
    to the bands; `AZURE2_NOT_SERVED` loses `fit_with_engine` then.
    **Closed** (section 11.6): the rows are made from `kinds()["fitters"]`, the job runs `fit_with_engine` with the
    job's progress and stop, Accept writes its end point; `SESSION_NOT_SERVED` no longer lists `fit_with_engine`.
    The sampler (emcee on `Objective.log_prob`) and the bands (`pyazr.bands`, and percentiles of draws) are jobs too
    (section 11.7).

Found while building the Results tab (step 3):

21. **The model on given energies and angles** (`Session.evaluate` on a grid, `POST /api/evaluate/grid`), or a data
    set's sub-point grid (`DataSet.grid`), so the model line passes between the data points.
22. **The lab energy (and lab angle) of each point** in the data arrays: the program writes them (`kinflag.dat`),
    pyazr gives c.m. values only.
23. **The program's σ and S columns** (`SigmaAndS` of `DataOutputFiles.cpp`: a data set of the S factor, the
    modified S factor's shape) and the normalized datum as arrays of the evaluation: the adapter forms datum × norm
    and × conversion itself for cross sections only, as the program does.

**Delivered by pyazr's `Objective` and `Evaluation`** (`pyazr3.md` section 7; 2026-11-13).  Closed: finding 5
(`Evaluation.data[i].chi2` is the row of `chiSquared.out`, data covariance by its per-set partition), finding 22
(`data[i].extra["energy_lab"]`, `["angle_lab"]`), the sampler's part of 20 (`Objective.log_prob`, the engine
sampler's rule, and a pool recipe in `pyazr/examples/fit_emcee.py`).  Partly: 23 (`extra["sfactor_conversion"]`
per point; the σ and S columns as the writer forms them, `SigmaAndS`, are not served yet), 14 (`extra["legendre"]`
gives each point's coefficients, not the order of a datum).  Open: 21 (the dense model line), bands, the rest.
The calls the server should make instead of today's (`pyazr/ui/engine.py`, not changed by that step):
- `Engine.evaluate` (`objective_terms_rwa`, `calculate_chi2_rwa`, `write_output_files` and a parse of
  `chiSquared.out`) → `ev = s.evaluate(x)`: `ev.objective`, `ev.chi2`, `ev.terms`, per data set `d.chi2`, `d.n`,
  `d.norm`, `d.norm_chi2` (the files are then needed only for the Outputs route);
- `Engine._dataset_arrays` (`segment_norms`, `segment_chi2`, `residuals`, `calculate_rwa`, `calculated_*`, `cross`)
  → the same `ev.data[i]`: `energies`, `excitation_energies`, `angles`, `model`, `values`, `errors`, `norm`,
  `residuals` (plain pulls, not whitened: the `RESIDUAL` text changes; their squares add up to `d.chi2` except where
  a data covariance joins sets), `extra`; the S factor from `extra["sfactor_conversion"]` as today;
- `model_data(self.s, x)` → `ev.model_data`;
- `open_objective` (`objective_and_grad`, `params_rwa`, `_free_names`) → `f = s.fit_objective()`:
  `f.value_and_grad`, `f.x0`, `f.names`, `f.bounds`; a sampler job on `f.log_prob` (done: `open_objective` is
  now `open_for_job` and `objective_of`, the sampler a job, section 11.7);
- `AZURE2_NOT_SERVED` (or a call): `fit_objective` is a public name of the session now, so
  `tests/pyazr/ui/routes_consistency_test.py` fails until the server calls it or lists it (one entry).

Smaller, in the parity list: an empty format-3 project (`Project.new` needs a source); the screening potential's
adiabatic limit from electron binding energies; Gaussian priors on level parameters for the sampler (records carry
priors of nuisances only).

Found while building the Physics tab (step 5):

24. **The session's physics functions as rows** (*closed*: section 11.9): `channel_functions` (and later
    `level_matrix`) have no row in the
    binding's table; the server writes the row of `channel_functions` by hand (`pyazr/ui/physics.py`
    `SESSION_ROWS`, from its docstring). Wanted: the session's functions described as the physics table describes
    its rows (arguments with unit, domain, example and sweep; returns with units), `J` among them.
    *Delivered in pyazr* (with 28): `Session.describe()["functions"]["channel_functions"]` (`pyazr/session_functions.py`):
    its arguments (`J` among them) and returns in the physics rows' shape.  `SESSION_ROWS` of `pyazr/ui/physics.py`
    can be made from it (the page's own fields `group`, `op`, `broadcast` added there).
25. **A status per element from the session's functions**: `channel_functions` returns `open` and 0 where the
    channel is closed; the tab turns `open` into status -2. Finding 12's remaining items (the level matrix, the
    transfer functions on potentials, structured stopping rows, the rate) stand; its "page code" is built (11.5).

26. **Stopping with the project's own targets**: the stopping rows take a formula; a target of the project (its
    layers, a user stopping table, a fitted stopping scale) has no call that evaluates its stopping on arrays.

Found while building the THM tab (step 5, second part):

27. **A schema of `<thm>`** (*closed*: `GET /api/thm/keys` serves `pyazr.kinds()["thm"]`, section 11.9): the engine
    reads its keys with a grammar of if-chains (`ParseThmOptionLine`,
    `ParseThmExperimentLine` in `src/ThmExperiment.cpp`) and describes none of them as data.  The page's forms come
    from a table written in `pyazr/ui/thm.py` (`KEYS`: value type of the schema, choices, unit, doc), checked against
    the list of keys the engine names when it refuses an unknown one.  Wanted: the keys of the options and of an
    experiment line as record kinds of `pyazr.schema()` (or `kinds()["thm"]`), read by the parser itself.
    *Delivered in pyazr*: `pyazr.kinds()["thm"]`, the rows of `src/ThmKeys.cpp` the parser reads (each key's
    description and reader in one row), in the shape of `KEYS`, equal to it today (`tests/pyazr/thm_keys_test.py`).
    The page replaces `thm.keys()` by `copy.deepcopy(pyazr.kinds()["thm"])` (or `KEYS = pyazr.kinds()["thm"]`); the
    reference page `docs/source/reference/thm_keys.rst` is generated from the same rows.
28. **The THM calls described as rows** (*closed*: the page's rows and answers are made from
    `Session.describe()["functions"]`, section 11.9): `thm_experiments`, `thm_background`, `thm_lineshape`,
    `thm_vertex`,
    `thm_distortion` return nested dictionaries; `pyazr/ui/thm.py` (`THM_ROWS`, `ANSWERS`) writes their rows from the
    docstrings and turns the answer into fields on the energies with lines (per channel and level, per exit pair),
    facts with units and tables.  Wanted: the rows in the session's description (as finding 24 asks for
    `channel_functions`) and the arrays per line, so that the server reshapes nothing.
    *Delivered in pyazr*: `Session.describe()["functions"]` has a row per THM call (arguments with unit and kind;
    returned fields with unit, `path` in the answer, `shape`, and `lines` -- the list a field runs over and the fields
    that label one line, e.g. `{"over": "channels[].levels[]", "label": ["jgroup", "channel", "J", "pi", "l", "s",
    "level"]}`).  The answers are unchanged; a test checks on four THM set-ups that every key of every answer is a path
    of its row.  `THM_ROWS` of `pyazr/ui/thm.py` can drop its argument and return lists (from the rows' `arguments`
    and `returns`) and keep only its page fields (`tab`, `op`, `max_elements`, `status_words`, `shown`); `ANSWERS`
    can drop `_lineshape` and `_vertex`'s loops by walking `lines.over` and labelling a line from `lines.label`, and
    `_distortion` (every return is a top-level key of shape `axis`); `_experiments` and `_background` stay tables.

## 13. Decisions for the PI

Each with the options, the recommendation and what changing it later costs.

1. **One page or several.** (a) one page with tabs and hash routes; (b) one page per tab. *(a)*: one state, one event
   stream, one set of keys. Later: small (the tabs are modules either way).
2. **Which fitters the UI runs.** (a) the engine's only; (b) the engine's and the recipes (iminuit, scipy, emcee) as
   rows of one job table. *(b)*: the recipes run today, the engine's need pyazr's call; the engine's stays the
   default because it is the program's reference. Later: small (a row).
3. **Where the sampler lives.** (a) the emcee recipe in a job worker; (b) the engine's `numcmc` behind a binding.
   *(a)*, as pyazr decided `numcmc` goes with Phase 9. Later: medium.
4. **Remote use** (server on a cluster, browser on a laptop through an SSH tunnel). (a) a goal now: loopback binding
   with a token, no absolute client paths, plotly vendored, small payloads; (b) later. *(a)*: it costs nothing beyond
   the token, and the measured payloads (27–48 kB per Results view) suit a tunnel. Binding to other addresses stays
   refused. Later: dear without the token (every client changes), small with it.
5. **Parity items dropped** (20, `parity_class.tsv`): old `<parameterSettings>` rows, `<externalCapture>`, segment
   keyword tails, their whitelist and block rules (the converter's), format-1/2 write-back artefacts (no rewrite of
   untouched lines now), segment renumbering and the status message about it (records are named), up/down order of
   segments, `<lastRun>` (a job carries its settings), `<mcmc>` (numcmc goes), the GSL warning on open (no run mask),
   the integrals file choice (cached by signature), the save-before-run rule (sessions read the working text), the
   AZUREOut file reader of the plot (arrays come from the API), the shifting of potential settings on pair deletion
   (potentials name their pair). *Recommendation: drop these 20*; the 20 missing items wait for pyazr, not dropped.
   Later: any can come back as an item.
6. **The command's name.** (a) `pyazr ui project.azr` (as `ARCHITECTURE.md` §11 and `pyazr3.md` 2.8 say); (b)
   `azure2 ui`, the program starting Python; (c) a separate `azure2-ui`. *(a)*: the UI is pyazr's, the program stays
   free of Python. Later: small (an alias).
7. **Where session options live.** (a) keys of `<config>` (formalism, transform, Coulomb table use), so the project
   is the whole state and the program reads them too; (b) a side file of the UI; (c) the page only. *(a)*; it is an
   engine change (new keys) and a choice of defaults for flags. Later: medium (projects written meanwhile lack them).
   **Taken by the PI on 2026-11-13: (a).** Built the same day (`project_format3.md` §3.12): eleven keys under the
   pyazr keywords; flag or keyword > key > default; the defaults unchanged.
8. **A token.** (a) always, Jupyter-like; (b) none on loopback. *(a)*, for shared machines (section 8). Later: dear
   to add once agents use the API without one.
9. **The default compute budget.** (a) 1 thread: one computation at a time, evaluations refused during a job;
   (b) the core count less one, with one thread kept for interactive evaluations. *(b) on a personal machine, with
   `--threads` to lower it; (a) wherever a budget is imposed, as here.* Later: none (a flag).
10. **Arrays: JSON only, or binary too.** (a) JSON, with sampler chains as files and quantiles in JSON; (b) a binary
    route for large arrays now. *(a)*: 2.3 shows JSON enough below 10⁵ numbers per view. Later: small (one route).
11. **plotly.js: which bundle and when to fetch it.** It needs a download (not allowed so far): (a) the cartesian
    partial bundle (MIT, about 1–1.5 MB), fetched once with the PI's yes before step 3; (b) the full bundle; (c) no
    plotly: SVG drawn by our own code. *(a)*. Later: small (a file).
    **Taken by the PI on 2026-11-13: (a).** Version 4.1.2 of `plotly.js-cartesian-dist-min`, fetched and checked
    against the registry's hash (`HANDOFF.md`, entry of that day); vendored with its licence by the step that
    draws the first plot.

## 14. Not looked at

A browser of any kind (none installed; the spike has no page); plotly.js and its sizes (not on the machine);
Playwright; T3 and the 183 projects (the spike used T1, T2, 13N); the engine's MIGRAD, MINOS and sampler through
Python, and the effect of `Config::stopFlag` on a running fit (neither reachable); several threads per worker and
jobs beside evaluations (the budget was one thread); Windows and macOS (process start-up by spawn was used, as there;
`PR_SET_PDEATHSIG` is Linux-only, the watchdog was not written); the LunaSRIM web app whose stack this follows (not
read); the logbook tool (a separate package); accessibility checked with a screen reader.
