# Tests of the web UI's server (`pyazr/ui`)

pytest files, each also a script registered in `tests/pyazr/CMakeLists.txt` (`ctest -R pyazr_ui`); without fastapi,
httpx, uvicorn, pytest or the engine module a file prints `skip: <reason>` and exits 77. The thread budget of the
apps is `OMP_NUM_THREADS` of the environment (ctest gives 1); no test sets it. `AZURE2_BIN` (ctest gives it) is the
program whose `chiSquared.out` the API's numbers are compared with.

- `api_test.py`: the API through FastAPI's test client on copies of T1 and T2.
- `jobs_test.py`: recipe fits as jobs on T1: to the end, stop, kill, resume, budget, accept, iminuit.
- `server_test.py`: a real server process on a free port above 20000 over HTTP; killed and stopped, its workers
  checked gone by PID; the user guide's session (`docs/source/user_guide/web_ui_api.rst`) run against it.
- `routes_consistency_test.py`: the route table, the app, the OpenAPI document, pyazr and the adapter agree.
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
