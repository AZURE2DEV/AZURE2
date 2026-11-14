# Tests of the web UI's server (`pyazr/ui`)

pytest files, each also a script registered in `tests/pyazr/CMakeLists.txt` (`ctest -R pyazr_ui`); without fastapi,
httpx, uvicorn, pytest or the engine module a file prints `skip: <reason>` and exits 77. The thread budget of the
apps is `OMP_NUM_THREADS` of the environment (ctest gives 1); no test sets it. `AZURE2_BIN` (ctest gives it) is the
program whose `chiSquared.out` the API's numbers are compared with.

- `api_test.py`: the API through FastAPI's test client on copies of T1 and T2.
- `jobs_test.py`: recipe fits as jobs on T1: to the end, stop, kill, resume, budget, accept, iminuit; the settings
  table with types; Accept writing what pyazr's `save_fit` writes, bit for bit, and Undo; a stopped job's best point
  accepted; the trace, the log tail, the long poll (`?wait=`), the parameter table during a job; iminuit's MINOS.
- `server_test.py`: a real server process on a free port above 20000 over HTTP; killed and stopped, its workers
  checked gone by PID; the user guide's session (`docs/source/user_guide/web_ui_api.rst`) run against it.
- `routes_consistency_test.py`: the route table, the app, the OpenAPI document, pyazr and the adapter agree.
- `static_test.py`: the page's files as served: the token exchanged for a cookie by a redirect, the
  Content-Security-Policy, every file the page loads present and tracked, plotly.js's hash and licence, a widget per
  value type of the schema, no engine kind named in the page's modules, lines of at most 120 columns (the vendored
  directory excepted: minified third-party code), the file-head route inside the root.
- `node_test.py`: runs `js/*.test.mjs` with node's own test runner (`node --test`; node from PATH or `NODE`; skipped
  without it) on the page's modules `forms.js`, `tables.js`, `api.js`, `plots.js`, `fitting.js` and the setting
  widgets of `tabs/fitjob.js`, the Physics tab's rules of `physicsform.js` (a form for every row of the engine's
  tables on every axis it allows, and for a row injected by the test; lines, gaps, the status line, the CSV, the
  bound, hash parameters), with the engine's schema, the server's job table and the function tables.
- `results_test.py`: the Results payload (`GET /api/evaluation`) against the program's own files: every data set's
  points, model, normalized datum and error and S factor against `AZUREOut_*` (menu 1) and every extrapolation
  against `AZUREOut_*.extrap` (menu 3) to eleven digits on T1 and T2; the χ² table against `chiSquared.out`; the
  residuals' squares against each data set's χ² (and the covariance group's sum on T2); model data on T3; the output
  files written into the project's output directory (`POST`/`GET /api/outputs`).
- `engine_jobs_test.py`: the engine's fitters as jobs: every fit row of `pyazr.kinds()["fitters"]` on T3 ends at
  the objective, end point, covariance and errors of `Session.fit_with_engine` called directly from the same start,
  bit for bit; the engine's MINOS as a job after an engine fit and as the fit's setting; a call limit said in the
  engine's words; Accept and Undo; on T2 a stop of MIGRAD and of LM answered within one progress record.
- `sampler_bands_test.py`: emcee and the bands as jobs on a copy of T3 with two free parameters: the chain's
  standard deviations against the fit's (Laplace) within 10 %, the chain's file and summary, no draws sent; a stop
  between steps and a chain continued; the linear band against the band from draws (median half-width ratio within
  10 %), a band for every data set and extrapolation on the Results payload's points; refusals of a band without
  its source; the start walkers inside the fit's limits.
- `physics_routes_test.py`: the physics routes: the rows are `pyazr.kinds()["physics"]`, the values a direct call's,
  NaN with the status outside the domain, the refusals, the bound, the session's channel functions.
- `physics_tab_test.py`: what the Physics tab asks of the server: every row plots with its examples and its axis's
  sweep, every element computed; the answer's Python text executed gives the route's arrays bit for bit (every row,
  a masked axis, several values as lines, a list of choices by name, an integer axis); the axis refused with what /
  allowed / fix; the session's functions as rows, their answer equal to the session's call, their Python text too;
  the stopping rows: SRIM's citation in every row, answer and Python text, the values bit for bit those of
  `pyazr.stopping`'s functions, a layer beyond the range status -2, a refused formula 422.
- `smoke_test.py` (label `browser`): a real server on T1 and T2, headless Chrome runs the page's `#selftest?edit=1`
  (every tab, a level, a data set's form and plot by plotly.js and as a table, a kind without layout hints, a new
  record's form, a width typed in and undone with the button; in the Fit tab a fit run, followed, stopped, accepted
  and undone with the buttons; in the Physics tab the first function, a Coulomb function with elements outside its
  domain, its table view and the session's channel functions, each with its status line and Python text) and takes a
  screenshot per tab (and of the masked Coulomb function)
  (`AZURE2_UI_SCREENSHOTS`, `AZURE2_UI_SHOTS=all`, `AZURE2_UI_SMOKE_PROJECTS`; see its docstring).
- `uitest.py`: shared helpers.

## Each rule and the error injected to show its test fails

Each line: the rule; the change made (then removed); the test that failed. Run on 2026-11-13 with
`work/uiserver/inject/inject.py` (outside the repository), at one thread; every injection gave `1 failed`.

- the user's file is written only by save: `State._push` also calls `project.save()`; fails
  `api_test.py::test_user_file_written_only_by_save`.
- a refused edit leaves text and revision: `State.edit` raises the revision before the edit; fails
  `api_test.py::test_refused_edit_leaves_state`.
- a whole text is validated before it enters the history: `State.replace_text` ignores the new errors; fails
  `api_test.py::test_refused_edit_leaves_state`.
- a stale revision is refused: `State.check_rev` never refuses; fails `api_test.py::test_stale_revision_refused`.
- a request over the budget is refused with the job's id: `Budget.require_session` never refuses; fails
  `jobs_test.py::test_budget_refuses_with_the_job`.
- binds to a loopback address only: `main` accepts any host; fails `server_test.py::test_non_loopback_host_refused`.
- the token is required on every route: the guard checks the token on `/api/project` only; fails
  `api_test.py::test_token_required_on_every_route`.
- the Host header names the loopback: the guard skips the Host check; fails
  `api_test.py::test_host_header_must_be_loopback`.
- no path outside the root is read or written: `State.inside_root` never refuses; fails
  `api_test.py::test_paths_outside_the_root_refused`.
- the route table and the handlers agree: a row names `pyazr.Project.levelz`; fails
  `routes_consistency_test.py::test_every_row_names_an_existing_callable`.
- the OpenAPI document lists every route: `/api/x` left out of the schema; fails
  `routes_consistency_test.py::test_app_serves_the_rows_and_openapi_lists_them`.
- workers die with the server: `die_with_parent` returns at once; fails
  `server_test.py::test_server_process_workers_die_with_it`.
- the chi2 per data set is the program's: the rows' chi2 taken from `evaluate().data[i].chi2`; fails
  `api_test.py::test_chi2_per_dataset_is_the_programs[T2]`.
- a stop ends a job before its next call: the tracker ignores the stop request; fails
  `jobs_test.py::test_stop_kill_and_resume_from_the_checkpoint`.
- no engine kind is spelled in the UI: the adapter drops data sets whose quantity is `"sigma"`; fails
  `routes_consistency_test.py::test_no_engine_kind_is_spelled_in_the_ui`.

Added after the merge of pyazr's errors and self-description (each shown failing the same way, with
`work/uiserver/inject/inject2.py`):

- every public exception class of `pyazr.errors` has a status other than 500, by the order of `errors.STATUS` (a
  failed lookup 404 though also a `ValueError`): the `KeyError` row moved below `ParameterError`; fails
  `api_test.py::test_every_pyazr_error_class_has_a_status`.
- an exception raised in the session worker maps by its class and bases: the worker sends its class name only;
  fails `api_test.py::test_error_statuses_and_fields`.
- a request with an `Origin` other than the server's own (scheme, host, port) is refused: the Origin check
  removed; fails `api_test.py::test_origin_must_be_the_servers_own`.
- a Host header `[::1]` (with or without a port) is a loopback name: the old parse (`rsplit(":", 1)`); fails
  `api_test.py::test_host_header_must_be_loopback`.
- an evaluation carries the engine's notices to the client: the worker drops its warnings; fails
  `api_test.py::test_evaluation_carries_the_engines_notices`.

Added with the page (each shown failing on 2026-11-13 with `work/uipage/inject.sh`, outside the repository; the
file restored from git after each):

- forms.js has a widget for every value type: the `grid` row renamed; fails `static_test.py::
  test_widget_table_has_every_value_type_of_the_schema` and `node_test.py` (forms.test.mjs).
- every file the page loads exists: `tabs/data.js` imports `../uix.js`; fails
  `static_test.py::test_every_file_the_page_loads_exists_and_is_tracked`.
- a Content-Security-Policy on every answer: the header line removed from the guard; the token leaves the address
  bar: the redirect of the page's row disabled; each fails `static_test.py::test_page_served_with_token_cookie_and_csp`.
- no engine kind named in the page: `tabs/data.js` spells `"dataset"`; fails
  `static_test.py::test_no_engine_kind_is_spelled_in_the_page`.
- the vendored plotly.js is the checked bundle: one word of its banner changed; fails
  `static_test.py::test_vendored_plotly_is_the_checked_bundle_with_its_licence`.
- an edit carries the revision: `api.js` does not add `rev`; empty cells sort last both ways: the guard removed from
  `sortRows`; the y unit of a plot is the data set's: the x unit used; each fails `node_test.py`.
- a script error in a tab is caught: `tabs/project.js` dereferences null; the Undo button takes the edit back: it
  sends redo; each fails `smoke_test.py` on T1 and T2.

Added with the Fit tab (each shown failing on 2026-11-13 with `work/uifit/inject.py`, outside the repository; the
file restored from git after each):

- a request's setting is converted by its type: `merged_settings` keeps the value as sent; fails
  `jobs_test.py::test_settings_listed_with_types_and_converted`.
- Accept writes the job's best point as pyazr's `save_fit` does: `job_accept` applies the working x instead; fails
  `jobs_test.py::test_accept_writes_what_pyazr_writes_back_and_undo_restores`.
- a stopped job keeps its best point: the stop drops it; fails
  `jobs_test.py::test_stopped_job_leaves_its_best_point_which_is_accepted`.
- the trace is thinned past `TRACE_POINTS`; a GET with `wait` answers at a state change; the parameter table is
  readable while a job holds the budget; a job's record has no `rev` (only `text_rev`): each in turn broken (no
  thinning; the wait ignores the state; the budget refuses the parameter table; the field named `rev`); each fails
  `jobs_test.py::test_trace_log_wait_and_parameters_during_a_job`.
- MINOS's intervals are in the result: the key renamed; fails `jobs_test.py::test_minos_after_an_iminuit_fit`.
- Stop is not offered while a job is stopping; a widget exists for every setting type; the trace's y axis is
  logarithmic; the job's columns appear only with a job: each broken in turn; each fails `node_test.py`
  (`js/fitting.test.mjs`).
- the Stop button of the job panel is found and pressed: its `data-action` removed; fails `smoke_test.py` on T1 and
  T2.

Added with the Results tab (each shown failing once, the file restored from git after each):

- the first rows of a file are bounded (1 to 200 rows, refused outside with the range; a line read to 4096 bytes and
  marked cut): the line read whole (`readline()` without its bound); the rows check removed; each fails
  `static_test.py::test_file_head_is_bounded`.
- the Results payload is the program's: the datum not multiplied by the norm (fails `results_test.py::
  test_payload_is_the_programs_output[T2]`, s11 with its THM norm); the S factor offered without the charge rule
  (fails `test_model_data_and_neutral_entrance`); the residuals scaled by 1.001 (fails `...[T1]`).
- the CSV keeps 17 significant digits (`toPrecision(15)` fails `js/plots.test.mjs`); residuals on the second y axis
  (drawn on the first fails `js/results.test.mjs`); level markers only against excitation energy (drawn on every x
  fails `js/results.test.mjs`).

The smoke test's "Failed to fetch" (2026-11-13): requests of the page cancelled in Chrome's renderer
(`Network.loadingFailed`, net::ERR_ABORTED, canceled; absent from Chrome's network log and from the server's request
log, which saw every other request answered), on T2 in 1 of 10 runs of this branch and 4 of 10 of the page's merge
commit, in real time as with the virtual-time budget. Not caused by the page's navigation (hash or history API: none
of about 2200 fetches in flight cancelled in a direct test), duplicate requests, keep-alive or the server. Caused by
Chrome's work in each run's fresh profile: 0 of 25 runs with `--disable-extensions`, 0 of 25 with
`--disable-component-update --disable-background-networking`, 6 of 25 with Safe Browsing off. The test's Chrome now
runs with those three flags (`QUIET_PROFILE` of `smoke_test.py`); a test artefact (a browser in use has loaded its
extensions before the page opens), not hidden by a retry. Scripts: `work/uiresults/diag/` (outside the repository).

Added with the Physics tab (each shown failing on 2026-11-13 with `work/uiphys/inject/inject.py`, outside the
repository, the file restored after each; every one failed):

- the copied Python makes the route's arrays: a log axis written as `np.linspace`; fails `physics_tab_test.py::
  test_every_row_plots_without_typing_and_its_python_is_the_route`.
- a log axis with from <= 0 is refused: the check removed; fails `physics_tab_test.py::test_axis_refusals`.
- the answer's arguments are those used: the axis echoed as sent; fails `test_masked_elements_lines_and_choices`.
- a session function's status comes from `open`: `open` ignored; fails `test_session_functions`.
- an element not computed is a gap: drawn with its value (a 0 of a closed channel); the status line lists each other
  status's words: words left out; the CSV keeps 17 digits: the shortest text; the CSV's arguments are broadcast to
  each element: the first value on every line; a kind of argument has its row: the `flag` row renamed; several
  values are dimensions of their own: all on the last axis; an integer axis has integer spacing and no points: the
  sweep's spacing; the page's bound is the server's: 10000; each fails `node_test.py` (`js/physics.test.mjs`).
- the page's status line counts what the answer computed: every element counted; fails `smoke_test.py` on T1.
- every row of the binding has a group, its axis a sweep, every argument a default or an example (a first plot
  without typing): the radius's example removed and the module rebuilt; fails `tests/pyazr/physics_table_test.py`.

Added with the stopping rows (`work/uiphys/inject/inject2.py`, the same way; every one failed):

- the Python text of a stopping row carries SRIM's citation: left out; a stopping choice by name reaches the library
  as that choice: always the first; each fails `physics_tab_test.py::test_stopping_rows_values_and_citation`.
- a text argument has its row in `ARGUMENT_KINDS`: the `text` row renamed; fails `node_test.py`.
- the citation is drawn under a stopping plot: not drawn; fails `smoke_test.py` on T1.

Added with the engine's fitters as jobs (each shown failing on 2026-11-14, the file restored after each):

- the rows of the engine's fitters are its table: every engine row given the `minos` setting; fails
  `routes_consistency_test.py::test_engine_rows_are_the_engine_fitter_table`.
- a stop request reaches a running engine fit: the worker's stop reader ignores the stop message; fails
  `engine_jobs_test.py::test_stop_ends_a_running_engine_fit_within_one_call` (the fit runs to its end).
- a row's facts are said as words, and progress numbers beyond the common four are shown by name: `factsText` says
  "yes" for false; `progressView` drops the further numbers; each fails `node_test.py` (`js/fitting.test.mjs`).

Added with the sampler and the bands (each shown failing on 2026-11-14, the file restored after each):

- a start outside the fit's limits is refused: the check removed; fails `sampler_bands_test.py::
  test_start_walkers_inside_the_limits_and_a_start_outside_refused`.
- a sampler from a sampler job continues its chain: the previous chain dropped after taking its last walkers; fails
  `sampler_bands_test.py::test_stop_between_steps_and_resume`.
- the linear band is the model ± sigma of `pyazr.bands`: sigma doubled; fails `sampler_bands_test.py::
  test_linear_band_and_band_from_draws_agree`.
- the band is a closed region (the lower edge drawn back); a band of other points is not drawn: each broken in
  turn; each fails `node_test.py` (`js/results.test.mjs`).
