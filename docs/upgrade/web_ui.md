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
  `masses`). E.
- Pairs: add, delete, channel regeneration (4) → Nucleus: no call. M.
- Levels and channels: table, edit, fix, Wigner limit (5) → Nucleus: `levels`, `add_level`, `remove_level`,
  `set_level`, `set_channel`; `Parameters.units`, `physics.wigner_limit`. E.
- J^π that give channels, channel generation, missing J^π (5) → Nucleus: no call (`add_level` takes the channels as
  given). M.
- IAEA levels (5) → Nucleus: `sources.ensdf` (today `nds.fetch_levels`). E.
- Hybrid potential (8) → Nucleus: `potential` records (`add_potential`, `set`, `remove`). E.
- Data sets, groups, samples, UPOS, sums, scattering lengths (6.1) → Data: `dataset`, `normalization`,
  `energy_scale`, `datum`, `target` records; `validate`; `describe`. E.
- Segments without data (6.2) → Data: `extrapolation` records. E.
- EXFOR (6.3) → Data: `sources.exfor` (today `nds.search_exfor`, `fetch_exfor`, `to_dataset`). E.
- Experimental effects (7) → Data: `target`, `beam`, `detector` records; `pyazr.stopping`. E (the screening
  potential's adiabatic limit M).
- Fitting settings table (9) → Fit: `Parameters` (`of_kind`, `names`, `units`, `bounds`) and a fit's errors. D.
- Edits of fitted values back into the records (9) → Fit: `save_fit` writes files only. M.
- Calculation types, minimisers, MINOS, band (10) → Fit: `Session.evaluate`, `extrapolations`, `fit_with_engine`
  (pyazr decision 8), `pyazr.bands`. D.
- Stop and the run log (10) → Fit: no call; the stop flag is not bound, the engine prints on standard output. M.
- Reaction rate (10) → Results: no call. M.
- MCMC (11) → Fit: the emcee recipe on `Objective.log_prob` (`pyazr3.md` 2.6), as a job. D, P.
- Plots (12) → Results: `Evaluation.data[i]`; `physics.sfactor_conversion`; `pyazr.plot` for PDF. E, D, P.
- THM workspace (13) → THM: `Project.thm`, `set_thm`, `set_thm_experiment`; `thm_vertex`, `thm_lineshape`,
  `thm_distortion`. E.
- THM binding energy and amplitude flag per channel (13) → THM: `set_pair` has neither. M.
- Title, About, help, theme, dialogs, screenshots (2, 14) → Shell: `pyazr.describe()`; the page. D, P.

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
  `set_thm_experiment`, `remove_thm_experiment`.
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
  more than 20000 elements 422) carry what / allowed / fix built from the row.  It runs in the server's process
  (`pyazr/ui/physics.py`): the functions need no session, and the C++ loop releases the GIL; the bound keeps one
  request near 1 s (coulsm, 45 µs per element; `shift_derivative` about five times that).
- `POST /api/session/channel_functions`, `/level_matrix` → `Session.channel_functions`, `level_matrix`.
  `channel_functions` is **served** in the session worker (adapter `Engine.channel_functions`, the budget as the
  other session routes), with the same bound of 20000 energies; `level_matrix` stays 501 (no engine call).
- `POST /api/stopping/{fn}` → `pyazr.stopping.*`: material, ion, energies → arrays.
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

The page shows `what` as the message, `allowed` and `fix` beneath it. pyazr 3.0 designs the three parts as text of
one message; the API needs them as fields of the exception (finding 9).

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
- **Session options** (Brune, Park, transform, Coulomb table, ...) are not in the project file today: they are
  program flags and session arguments. For "the project is the state" they need a home (decision 7).

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
tabs and a status line (the last message of the engine or a job). Keys: `1`–`7` switch tabs, `Ctrl+Z`/`Ctrl+Shift+Z`
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
│ ¹ adding and removing pairs with channel regeneration waits for pyazr (finding 10)                 │
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
Routes: `/api/physics*`, `/api/session/channel_functions`, `/level_matrix`, `/api/stopping/*`, `/api/data/{name}/grid`
and `/kernel`.

```
┌ Physics ──────────────────────────────────────────────────────────────────────────────────────────┐
│ function [coulomb_waves ▾]   search [        ]                                                     │
│   kinematics · Coulomb · channel functions · widths · resonance quantities · transfer/THM ·       │
│   stopping · this project: channel functions, level matrix, data-set grid and kernel              │
│ ─ coulomb_waves(l, eta, rho, method) ─ "F, G, F′, G′ on ρ, by every method of pycoulomb" ───────  │
│ l    [0] … [3]  (integer ≥ 0)                                                                      │
│ eta  [0.5]                       (≥ 0)                                                             │
│ rho  range [0.01] – [20] points [400] [log]   ← the plotted axis                                   │
│ method [all ▾]                                                                                     │
│ ┌ F_l, G_l against ρ, one line per l and method ───────────────────────────────────────────────┐   │
│ └──────────────────────────────────────────────────────────────────────────────────────────────┘   │
│ status: 400 of 400 elements ok                         [Export (CSV)] [SVG] [Copy as Python]       │
└────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

"Copy as Python" gives the pyazr call that makes the same arrays: the page teaches the API.

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
│ Channels  pair  B (MeV)  amplitude?   (needs pyazr: finding 11)                                    │
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
  fixed one for agents and tests.
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
- **A session option**: none. The Project tab's options come from `kinds()["options"]`.

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
   (pyazr step 3) replaces the adapter's calls.
4. **Fit.** `tabs/fit.js`, job kinds: recipe fits and MINOS (iminuit), the emcee sampler, bands; then the engine's
   fitters once pyazr has `fit_with_engine` with progress and stop (decision 8 and finding 2). About 1000 lines.
   Recipes need `Objective` (pyazr step 3; today `objective_and_grad` serves); engine fits need pyazr's call.
5. **Physics and THM.** `tabs/physics.js` generated from the physics table (pyazr step 4; before it, stopping only),
   `tabs/thm.js` (today: `set_thm`, the session's THM calls). About 800 lines.
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
stopping, EXFOR, ENSDF, files, events.  (Served since: a physics function on arrays, the channel functions.)

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
9. *Session options*: the engine's table (`pyazr.kinds()["options"]`: keyword, type, default, flag, doc) with the
   value in use, in `GET /api/describe` and `GET /api/session`; the server uses the defaults. `SESSION_OPTIONS` in
   `state.py` is the provisional place of other values (decision 7 not taken); not editable through the API.
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
3. **A record of the engine's progress lines** (`Iteration: N Chi-Squared: X`) as data, not text on standard output.
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
10. **Pairs added and removed, with channel generation**: `add_pair`, `remove_pair`, and the channels a J^π allows from
    the pairs and `l_max` (Qt's `LevelsTab` rules; `add_level` takes channels as given).
11. **THM binding energy B and the width-is-amplitude flag** per channel: `set_pair` and `set_channel` have neither.
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

Smaller, in the parity list: an empty format-3 project (`Project.new` needs a source); the screening potential's
adiabatic limit from electron binding energies; Gaussian priors on level parameters for the sampler (records carry
priors of nuisances only).

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
