# pyazr/ — the Python API

Read the root `CLAUDE.md` first. pyazr is how most users script AZURE2, so its
public names are a contract.

## Using pyazr as an agent: ask the package first

- `help(pyazr)`: the objects of today's API and a ten-line session (tested by `tests/pyazr/describe_test.py`).
- `pyazr.describe()`: version, engine build, Coulomb table in use (path, content identity or "exact"), stopping and
  mass data, threads.  `pyazr.describe(path)` is still the engine's description of a project.
- `pyazr.kinds()`: record kinds and fields, data kinds (units, frames), parameter kinds (unit, owner fields), session
  options (type, default, CLI flag, doc), bound physics functions, error classes.  Read from the engine's tables
  (`described.py`); never add a hand-typed list of kinds in Python.
- `Session(...).describe()`: pairs, levels, data sets, parameters (their records: name, kind, unit, x, free, bounds,
  owner fields), options in effect and their sources.  The options of the calculation are `<config>` keys too
  (`session_options.py`: a keyword given wins over the key, the key over the default; the signature's defaults are
  unchanged, a keyword left out is not given).
- `Session(...).parameters` (a `Parameters`, `parameters.py`): the table of the fit parameters, one row per
  parameter, built once per session and layout from the engine's records (`Session.parameter_records()`; the pairs
  from `Session.pair_records()`).  The field names, units and kinds are one table of the binding
  (`api/src/pybind_records.cpp`, `pyazr.kinds()["parameter_fields"]`, `["pair_fields"]`): never type a field name, a
  column position or a numeric kind code in Python; a new field is a row there.  The fit vector `x` is the only
  vector of every call (`x0`, `put`, `take`, `read`, `write`; a selection is `of_kind`, `of_level`, `free`);
  `print(parameters.table())` shows what each row is.  Physical values are `Session.transform(x)` (one per row, the
  engine's transformation); `Session.describe()["parameters"]` is the same records as plain data.  A parameter file
  with its fixed values is taken when the session is made (`Session(file, parameter_file=path)`,
  `parameter_file.py`); no call sets a vector inside a live session.
- Errors (`errors.py`): every deliberate error is a `pyazr.Error`; `FormatError`, `ParameterError` (with
  `ParameterKeyError`, `ParameterIndexError`, `ParameterVectorError`), `EngineError` (= `AZURE2Error`, one class,
  which keeps that name), `DomainError` (`DomainZeroDivision`), `SessionClosed`, `MissingFile`; engine notices are
  `EngineWarning`.  A new error derives from the class above and from the builtin the call raised before; its
  message says what was wrong, what is allowed and which call shows or fixes it (`tests/pyazr/errors_test.py`
  holds the table of wrong calls); give the three parts as keywords too, `Error(message, what=, allowed=, fix=)`, each
  a substring of the message (the test checks it), so that a caller (the web page) reads them as `e.what`,
  `e.allowed`, `e.fix` (a `FormatError` takes them from its first diagnostic).
  The argument checks of the session's calls live in `session_checks.py`, not in `session.py`.
- An engine notice becomes a warning when its line is a row of `api/src/pybind_notices.cpp` (the session's
  `Config::outStream` passes through that buffer); every other line still goes to standard output.
- `pyazr.physics.<name>`: the engine's physics functions on arrays (`help()` gives arguments, units, domains);
  `pyazr.kinds()["physics"]` the same rows as data; `errors="mask"` gives NaN and a status per element instead of
  `DomainError`.  The per-session ones (the session's Coulomb table or potential) are `Session(...).channel_functions`.

## Recipe: add a physics function

1. The library computes it (`physics/`, `experiment/`, `stopping/`); if it does not, a thin free function there with
   its unit test (as `physics/include/azure2/physics/Widths.h`).  New physics is a design note first, not a row.
2. One row in `BuildTable()` of `api/src/pybind_physics.cpp`: name, summary, doc, library function, axis, arguments
   (`Real`/`Integer`/`Flag`/`Choice` with unit and domain: `.gt`, `.ge`, `.lt`, `.le`, `.otherwise` for a default),
   returns (name, unit, doc), and an adapter that converts the arguments and calls the library (no formula there);
   the adapter returns 0, `kJointly` (-2) for arguments outside the domain together, `kFailed` (-1), or the 1-based
   index of an argument.
3. A `row_test("<name>", ...)` in `tests/pyazr/physics_table_test.py` against mpmath, a published table or the
   engine's own use on a project (consistency check 41 fails without one).
4. `python3 tools/gen_physics_docs.py` (the user-guide page; the test fails until it is run).
Nothing in Python: `pyazr/physics.py` builds the function, its signature and docstring from the row.

## Recipe: add a Project edit

An edit of the nuclear part (`<levels>`, `<config>`, `<thm>`) is the engine's: a function in `experiment/` (the
`<levels>` rows: `experiment/src/LevelEdits.cpp`, as `SetPair`, `AddPair`, `RemovePair`) that reads the document,
refuses with one sentence (what was wrong, what is allowed, which call fixes it) and returns line replacements; a
binding in `api/src/pybind_nuclear_part.cpp` (`*_edits`, keyword arguments, `ValueError` with the reason); one method of
`NuclearPart` (`pyazr/nuclear_part.py`) that names things (a level by label, a channel by pair, l, s) and hands the
replacements to `Project._apply`, which validates the whole project as it would stand.  No column, position or grammar
in Python.  A new column of a pair is a field of `PairChange` / `NewPair` and of `PairInfo` (`LevelTables.h`, read by
`PairInfoOfRow`), its key in `level_table`'s pair dict and its field of `ProjectPair`.  Test: the edited file gives the
engine's result of the same edit made by hand (`tests/pyazr/thm_pair_edits_test.py`), add then remove gives the text
back, each refusal leaves the project unchanged (`tests/pyazr/pair_edits_test.py`).

## Recipe: describe a session call

A diagnostic of the session (as `channel_functions`, the THM calls of `session_thm.py`) is a method where it
computes and one row of `SESSION_FUNCTIONS` (`pyazr/session_functions.py`): arguments (kind, unit, domain, default,
example, sweep), returns with unit, `path` in the answer, `shape` and, under a list, `lines` (what one line is, the
fields that label it).  `Session.describe()["functions"]` serves the rows; the answer of the call is not reshaped for
them.  Its case in `tests/pyazr/session_functions_test.py` runs the call and fails while a key of the answer is not a
path of the row, or a path of the row never appears.

## Recipe: use the objective and the evaluation

- A fit or a sampler: `f = s.objective()` (`objective.py`): `f(x)`, `f.grad`, `f.value_and_grad`,
  `f.residuals`/`f.jacobian(x, skip=None)` (the objective as rows, `sum(r**2) == f(x)`), `f.log_prob` (the engine's sampler's
  rule: -f/2, -inf outside its support), `f.x0`, `f.names`, `f.bounds`, `f.errordef`.  Recipes:
  `examples/fit_scipy.py`, `fit_iminuit.py`, `fit_emcee.py` (tested by `tests/pyazr/fit_recipes_test.py`).
- The program's fitters on a session: `r = s.fit_with_engine(method, minos=, max_calls=, progress=, stop=)`
  (`session_fit.py`, rows of `kinds()["fitters"]` with their columns `call_limit` and `fallback`).  Every fitter's
  `EngineFitResult` has the same fields: `x`/`objective` the lowest objective reached (at a call limit the lowest
  evaluated), `status` the fitter's word, `outcome` the shared word (`fitStatus.out`'s status), `notices` (fallback,
  a setting ignored, the best point kept), `minos` with `minos_valid` per side.  A new field is a member of
  `FitResult` (`include/FitControl.h`), a key of `FitWithEngineDict` and a field here.
- The numbers at one x: `ev = s.evaluate(x)` (`evaluate.py`), one engine evaluation: `ev.objective`, `ev.chi2`
  (`Total-Chi-Squared`), `ev.terms` (by the engine's term names), `ev.data[i]` (model, residuals, chi2 as the row of
  `chiSquared.out`, `rows` its rows of `f.residuals(x)`, extra columns), `ev.table()`; of a session made with
  `Session(file, extrapolations=True)` the extrapolations, with every chi-squared field `None`.  These two are the
  only ways to evaluate: a number of `chiSquared.out` or `AZUREOut_*` is read from them, never rebuilt.  A session
  never switches between data and extrapolations: the other is a second session
  (`parameters.carry_model_parameters` carries a data session's vector to it).
- **A new term of the objective**: a class in `src/Objective.cpp` (its `Name()`), nothing in `pyazr/`; it appears in
  `Evaluation.terms`; if it is a sum of squares, its rows go into `AZUREAPI::ObjectiveRowsRWA`
  (`api/src/AZUREAPIObjective.cpp`; `objective_object_test.py` checks `sum(r**2) == f(x)`).
- **A new number a point carries** (an extra column of a data kind): one row of `PointOutputs()` in
  `api/src/pybind_evaluation.cpp` (name, unit, when, doc, the reader); it appears in `Evaluation.data[i].extra` and
  `pyazr.kinds()["data_outputs"]`.
- **A new field of `Evaluation`**: an entry of the dict `EvaluationDict` builds (`pybind_evaluation.cpp`), the
  attribute in `evaluate.py` and its line in the class's `FIELDS` (`kinds()["objects"]`).
- Calls of both objects take the session's lock (`pyazr.locking.session_lock`), as every method of the session that
  enters the engine does.
- `Session._engine` is the binding's session object (`_azure2.Session`): private, no lock, no promise of stability.
  A module of `pyazr/` reads it only inside an `engine_call` method or a `with session_lock(s):` block
  (`tests/pyazr/session_lock_test.py` checks every module by its syntax tree); a test may use it only to test the
  binding itself.

## How it is built

- `session.py` drives the engine **in-process** through the pybind11 module
  `_azure2` (`api/src/pybind_api.cpp` → `api/src/AZUREAPI.cpp`). It runs no
  subprocess and opens no socket. Rebuilding with CMake (`-DUSE_API=ON` and
  pybind11 found) writes `pyazr/_azure2*.so` here. That file is git-ignored;
  never commit it.
- pyazr reads project format 3 only.  `project.py` (`Project`, records) and
  `nuclear_part.py` (`<levels>`, `<config>`, `<thm>`) edit the file on the
  engine's reading of it: the lossless document (`_azure2.Document`), the
  record grammar and validator, the level edits (`experiment/LevelEdits.h`),
  the `<config>` and `<thm>` edits. They need the compiled module. Never add a
  Python parser or a column position for a section: bind the engine's. A
  project must round-trip **byte for byte** (any layout, any line ending), and
  an edit changes only the lines it edits (`tests/pyazr/project_test.py`,
  `azr_layout_test.py`). A format-1/2 file is refused in one place,
  `project.require_format3` (the engine's message, io::NotFormat3Message, names
  `AZURE2 --convert OLD NEW`); `convert.py` is the converter's binding
  (`convert_azr`, `convert_params`).
- No physics formula in Python: physical values are `Session.transform(x)` (the engine's `TransformOut`), the
  functions the rows of `pyazr.physics` or `Session.channel_functions`.  The one formula of pyazr's own is the
  Weisskopf estimate (`widths.weisskopf_width`, tested against its closed form by `tests/pyazr/weisskopf_test.py`:
  the engine has none).  Never add another: a missing function is a row of the physics table.
- `sources/nds.py` (EXFOR/ENSDF client) is the only EXFOR reader on this branch (the
  Qt GUI's `gui/src/ExforData.cpp` exists on `dev` only).

## The wheel (`pyproject.toml`)

- `pip install .` / `pip wheel .` build with scikit-build-core: the same CMake files with `BUILD_TESTING`
  off (`USE_MCMC` on: the binding needs the sampler), targets `_azure2` and `AZURE2`, the install rules of component `python` (module into the
  package, program into the scripts), the package directory copied as it is less `wheel.exclude`.  The in-tree build
  (`work/build.sh`, `PYAZR_OUTPUT_DIR`, `PYAZR_ENGINE_DIR`) is what every test uses; the wheel is tested by
  `pyazr_wheel` (slow: builds, unpacks, never installs) and `tests/pyazr/wheel_smoke.py` (the workflow's test).
- **A run-time data file**: put it under `pyazr/`, read it from next to its module (`__file__`), and add its row to
  `DATA_FILES` in `datafiles.py` (or to `NOT_SHIPPED` and `wheel.exclude` if the wheel must not carry it); check 48
  and `pyazr_wheel` fail otherwise.  Never an environment variable for a file of the package.
- The Coulomb table: `coulomb_table.py` (`where`, `fetch`, `build`; `pyazr coulomb-table ...`).  `TABLE_URL` and
  `TABLE_SHA256` are set together by the owner when he publishes the file; never download anything in a test (the
  test serves a small table from `http.server` on 127.0.0.1).

## Rules

- **Public API = `__all__` in `__init__.py`.** Adding a name is a minor version
  bump of `__version__` in the same file. Renaming or removing one within a
  major version needs a deprecation period first (3.0 removed the old-format
  API at once: the converter handles old files): keep the old name, have it emit
  `DeprecationWarning`, and note it in the docstring for at least one minor
  release. `tests/consistency` checks that every exported name exists.
- **Units and frames are in the name or the docstring** of every argument and
  return value, for example "lab energy, MeV" or "centre-of-mass, barns/sr".
  Return NumPy arrays, not lists, for numeric columns.
- Every public function has a docstring that says what it computes, its
  units, and what it raises. Match the style of the module you are editing;
  numpydoc sections (Parameters / Returns) are used in `session.py` and
  `bands.py`. A new public name also gets an entry in
  `docs/source/user_guide/pyazr.rst`. More than half of `__all__` is not
  documented there yet, so document what you touch.
- Errors: raise a specific exception with a message that names the file, line
  or parameter. Never return `None` or `nan` to signal failure, and never
  `print` a warning where `warnings.warn` belongs.
- Python ≥ 3.10 (`pyproject.toml`). The core dependency is numpy only (`import pyazr` imports neither scipy nor
  mpmath). Anything else is an optional extra (scipy: `web`, `examples`) and is imported inside the function that
  needs it.
- Engine semantics (χ², penalties, normalizations, shifts, parameter order) are
  defined by the C++ code. pyazr must reproduce them exactly. Tests compare
  against the engine's `chiSquared.out` and parameter files, not against
  hand-computed expectations.
- `pyazr/examples/*.py` are user documentation. When an API changes, run the
  examples it affects.

## Tests

`tests/pyazr/<name>_test.py`, a plain script registered in
`tests/pyazr/CMakeLists.txt`, with `OMP_NUM_THREADS=2` and `LABELS slow` if it
takes more than a minute. Without the engine module, a test prints `skip:` and
exits 0. Read the output and confirm it says `ok`, not `skip`.

## The command line (`python -m pyazr`, `cli.py`)

`python -m pyazr COMMAND ... [--json]` (`__main__.py` calls `cli.main`): describe, validate, convert, schema, kinds,
stopping, coulomb-table check, import exfor, import ensdf, ui, commands.  With `--json` one document on standard output
(`{"command", "ok", "exit", "result"}`, and `"error": {"class", "what", "allowed", "fix"}`); exit 0 done, 1 failed,
2 a wrong call or a refused input, 3 an old project (the program's project commands).  `cli.py` imports nothing of
the engine at import time (check 42 of the consistency test reads its table without it).

- **A command of the command line**: one `Command` row of `COMMANDS` in `cli.py` (name, `Arg` rows with a type of
  `ARG_TYPES`, the pyazr function it `calls`, the help line) and its `run` function, which only calls that function
  (a `text` printer and a `status` when the result can be a failure, as validate's errors); its case in `JSON_CASES`
  of `tests/pyazr/cli_test.py` (run with `--json` on T1 or T2, the JSON read); then regenerate the reference page:
  `python3 -c "from pyazr import cli; cli.write_reference_page()"` (from the repository's root, with the engine module)
  (the consistency test, check 42, fails for a row without a case or a stale page).  A command that would need more
  than a call and a printer shows a gap of pyazr: add the function to pyazr first, or give the row `missing=(what,
  allowed, fix)` (it then exits 2 with them).
- A new argument type is one row of `ARG_TYPES` (description, converter raising `CallError` with its three parts,
  the value's name in a usage line).

## The web UI's server (`pyazr/ui`)

A subpackage `import pyazr` does not load (it needs fastapi and uvicorn: the extra `web`). Design and status:
`docs/upgrade/web_ui.md` (section 11.1 is what is built). An agent that works beside the user on the user's live
project starts or reuses the server and asks it what it serves:

```bash
python -m pyazr.ui project.azr --port 28400 --token "$TOKEN"   # prints http://127.0.0.1:28400/?token=...
curl -s -H "Authorization: Bearer $TOKEN" http://127.0.0.1:28400/api/describe   # every route, its call, status
```

The URL of a running server is in `.<name>.ui/url` beside the project (mode 0600). Every edit carries the revision
it was made on (`"rev"`); refusals are JSON with `error`, `message` and the engine's `diagnostics`. An agent that
does not need the user's live state uses pyazr directly.

- `engine.py` is the only file that names the class `Session`. Every public method of `Session` is either called
  there or listed in `SESSION_NOT_SERVED` with a reason.
- **A diagnostic of the session on a tab** (as the THM tab's): its method and its row in
  `pyazr/session_functions.py` (arguments of a kind `ARGUMENT_KINDS` of `static/physicsform.js` knows; returns with
  `path`, `shape` and `lines`) are all the page needs: `Engine.session_function` calls it by name, the server walks
  the answer by the row (`physics.answer_parts`) and the Physics tab draws it.  To draw it on another tab, or to give
  it status words, a bound or fields left off the first plot, add its entry to `PAGE` (`physics.py`; the THM tab's
  in `thm.py`, with the two answers made into tables there).  Test it against the direct call in
  `tests/pyazr/ui/thm_routes_test.py`.  A key of `<thm>` is a row of `src/ThmKeys.cpp` (its description and its
  reader), served as `pyazr.kinds()["thm"]` and by `GET /api/thm/keys` as it is.
- A route is one row of `routes.ROUTES`; a job kind one row of `jobs.JOBS`.
- **A job kind or a recipe** (a new fitter, MINOS of another library, a sampler): a function `(tracker, x0, settings,
  names) -> (status, message, extra)` in `jobs.py` that calls `tracker.value` / `tracker.grad` (or `tracker(x)` for
  both: the stop request and the progress come with it) and one `JobRow` of `JOBS` naming it, the package it imports
  (`needs`) and its `Setting` rows (name, type, default, doc); `follows=(kind, method)` makes the Fit tab offer it
  beside that job's result. `extra` may carry `errors` {name: value}, `covariance` and `correlation` ({names, matrix})
  and `minos` {name: {lower, upper, valid}}: the Fit tab shows these keys. A row the engine cannot run yet names the
  pyazr call it waits for in `missing`. A new setting type is a row of `SETTING_TYPES` and a widget of
  `SETTING_WIDGETS` in `static/tabs/fitjob.js` (the node test fails until it exists). Test it in
  `tests/pyazr/ui/jobs_test.py` on T1 at one thread (a fit ends within 20 s there). No file of `pyazr/ui` names a
  data, record or parameter kind of the engine (`tests/pyazr/ui/routes_consistency_test.py` checks all of this).
- The server never builds a session (a build holds the GIL); workers are separate processes and die with it.
- Tests: `tests/pyazr/ui/*_test.py` (pytest files run as scripts, `ctest -R pyazr_ui`).

## The web UI's page (`pyazr/ui/static`)

Plain JavaScript modules served as they are (no build step, no package manager, nothing from the network); design
`docs/upgrade/web_ui.md` (section 11.2 is what is built). The page holds no state of the project: it shows what the
server answers, sends every edit with its revision and re-reads; it never computes a physical number.

- `app.js` (header, the tab table `TABS`, hash routes, keys `KEYS`), `api.js` (fetch, revision, errors),
  `tables.js`, `forms.js` (schema → form: the widget table `WIDGETS`), `formview.js` (draws a form), `layout.js`
  (the only per-kind table: which fields share a line), `plots.js` (the plot table `PLOTS`, palette, CSV), `ui.js`,
  `tabs/*.js`, `selftest.js`, `vendor/` (plotly.js 4.1.2, MIT; its hash is checked by `static_test.py`).
- **A widget** (a new value type of the schema): one row of `WIDGETS` in `forms.js` turning (field, text, ctx) into
  controls and back into text; `static_test.py` and `js/forms.test.mjs` fail until it exists.
- **A tab**: a module `tabs/<name>.js` exporting `render(ctx, element, args)` and one row of `TABS` in `app.js`
  (replacing its `later` text); it appears in the self-test by itself; add its route to `VIEWS` of `smoke_test.py`.
- **A plot row**: one row of `PLOTS` in `plots.js` (`applies`, `build` → axes with quantity and unit, traces with a
  role: `data` with errors, `model` (dash for a grid), `residual` on the axis under the plot, `levels`) and its case
  in `js/plots.test.mjs` or `js/results.test.mjs`; the figure, the CSV (17 significant digits) and the table without
  plotly follow from the description.  A new x or y choice of the Results tab is one row of `X_CHOICES` or
  `Y_CHOICES` (an array of the payload; a converted y needs its arrays from the server, `engine.py`).
- **The Results payload** is `Engine.evaluation` (`GET /api/evaluation`): every number is the engine's or pyazr's;
  the only products formed are datum × norm and × S-factor conversion, as the program forms them for `AZUREOut_*`
  (`tests/pyazr/ui/results_test.py` compares them with the program's files).
- No module but `layout.js` names a record kind, quantity or parameter kind (`static_test.py` checks it); the
  sections the Nucleus tab shows are `NUCLEUS_SECTIONS` in `tabs/records.js`, every other section is a group of Data.
- Tests: `ctest -R pyazr_ui` (with `node` on PATH or `NODE=...` for `pyazr_ui_node`). The smoke test:
  `AZURE2_UI_SCREENSHOTS=<dir> [AZURE2_UI_SHOTS=all] python3 tests/pyazr/ui/smoke_test.py` (or `ctest -L browser`):
  a server on T1 and T2, headless Chrome (`CHROME=...` to choose) runs `#selftest?edit=1` and takes one screenshot
  per tab into `<dir>` (a temporary directory otherwise); more projects with `AZURE2_UI_SMOKE_PROJECTS=a.azr:b.azr`.
