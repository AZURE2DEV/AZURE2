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
  widgets of `tabs/fitjob.js`, with the engine's schema and the server's job table.
- `smoke_test.py` (label `browser`): a real server on T1 and T2, headless Chrome runs the page's `#selftest?edit=1`
  (every tab, a level, a data set's form and plot by plotly.js and as a table, a kind without layout hints, a new
  record's form, a width typed in and undone with the button; in the Fit tab a fit run, followed, stopped, accepted
  and undone with the buttons) and takes a screenshot per tab
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
- the chi2 per data set is the program's: the rows' chi2 taken from `segment_chi2`; fails
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

INJECTIONS_PLACEHOLDER
