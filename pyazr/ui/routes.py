"""The route table of the web UI's API (docs/upgrade/web_ui.md sections 4.1 and 4.2) and its generic handlers.

Every route is one row of :data:`ROUTES`; there is no handler written per route.  A row names the callable it serves
and how a request's path parameters, query and JSON body map onto the callable's arguments; one of five generic
handlers runs it:

- ``project``: a method of :class:`pyazr.Project` on the open project.  A GET reads; any other method is an edit
  made on the revision the request names (``rev`` in the body or query, or an ``If-Match`` header), validated by the
  engine before it enters the history (:meth:`pyazr.ui.state.State.edit`), the text values it carries checked to
  name no path outside the root;
- ``session``: an operation of the adapter :class:`pyazr.ui.engine.Engine` in the session worker
  (:meth:`pyazr.ui.server.Server.session_request`);
- ``server``: a method of :class:`pyazr.ui.server.Server` (the history, the working x, jobs, the server's own
  description); ``serves`` names the pyazr calls it makes;
- ``missing``: a call pyazr does not have yet (501 naming it, web_ui.md section 12);
- ``later``: a route of a later step of web_ui.md section 11 (501 naming the step and today's call).

An argument spec is ``name`` or ``name:int`` / ``name:float`` (converted); positional ``args`` are looked up in the
path, the query, then the body; ``kwargs`` in the body (or query) and passed only when given; ``spread`` names the
body key whose object is passed as keyword arguments (a record's fields); ``fixed`` adds constant keywords.
The consistency test (tests/pyazr/ui/routes_consistency_test.py) holds this table, the app and pyazr together.
"""
import dataclasses
import json
import math
import warnings
from dataclasses import dataclass, field

from .errors import Missing, to_http


@dataclass(frozen=True)
class Route:
    method: str
    path: str
    handler: str        # project | session | server | missing | later
    call: str           # pyazr.Project.<m>, engine.Engine.<op>, server.Server.<m>, or the missing pyazr call
    doc: str
    args: tuple = ()
    kwargs: tuple = ()
    spread: str = None
    fixed: dict = field(default_factory=dict)
    serves: tuple = ()  # the pyazr calls a server row makes
    text: bool = False  # the response is plain text


P, S, V, M, L = "project", "session", "server", "missing", "later"

ROUTES = (
    # -- tables of the engine and the server
    Route("GET", "/api/describe", V, "server.Server.describe",
          "The server: pyazr.describe() (versions, build, Coulomb table, data in use), the session options, the "
          "compute budget, the open project, the job table and this route table.", serves=("pyazr.describe",)),
    Route("GET", "/api/server", V, "server.Server.processes",
          "PIDs and memory (kB) of the server and its worker processes."),
    Route("GET", "/api/schema", V, "server.Server.schema", "The engine's schema of the project format.",
          serves=("pyazr.schema",)),
    Route("GET", "/api/kinds", V, "server.Server.kinds",
          "Record, data, parameter kinds, session options, physics functions and error classes of the engine.",
          serves=("pyazr.kinds",)),
    # -- the project: the state
    Route("GET", "/api/project", V, "server.Server.project_state",
          "Path, root, revision, unsaved, can undo / redo, the crash copy found at start-up."),
    Route("POST", "/api/project/open", V, "server.Server.open", "Open another project inside the root.",
          args=("path",), serves=("pyazr.Project.open",)),
    Route("POST", "/api/project/new", V, "server.Server.new",
          "A new project at path from the nuclear part of another (unsaved until saved).",
          args=("path", "nuclear_part"), serves=("pyazr.Project.new",)),
    Route("POST", "/api/project/save", V, "server.Server.save",
          "Write the project to its file, or to path inside the root: the only call that writes the user's file.",
          kwargs=("rev", "path"), serves=("pyazr.Project.save",)),
    Route("POST", "/api/project/restore", V, "server.Server.restore",
          "Take the crash copy (the working file of a server that ended with unsaved edits) as an edit.",
          kwargs=("rev",)),
    Route("GET", "/api/project/text", V, "server.Server.text", "The project's text.",
          serves=("pyazr.Project.text", "pyazr.Project.to_bytes"), text=True),
    Route("PUT", "/api/project/text", V, "server.Server.put_text",
          "Replace the whole text (agents): validated by the engine, refused with its diagnostics.",
          args=("text",), kwargs=("rev",), serves=("pyazr.Project.validate",)),
    Route("GET", "/api/project/validate", P, "pyazr.Project.validate",
          "Every problem of the project as the engine's validator reports it."),
    Route("GET", "/api/project/describe", P, "pyazr.Project.describe",
          "What the engine makes of every record.", fixed={"as_json": True}),
    Route("GET", "/api/project/describe/{name}", P, "pyazr.Project.describe",
          "What the engine makes of one record.", args=("name",), fixed={"as_json": True}),
    Route("GET", "/api/project/records", P, "pyazr.Project.records", "Every record: kind, name, line, fields."),
    Route("POST", "/api/project/records", P, "pyazr.Project.add",
          "Add a record: kind, name, fields ({key: text or number or list}).", args=("kind", "name"),
          spread="fields"),
    Route("GET", "/api/project/records/{name}", P, "pyazr.Project.__getitem__", "One record.", args=("name",)),
    Route("PATCH", "/api/project/records/{name}", P, "pyazr.Project.set",
          "Edit fields of a record (null removes a field).", args=("name",), spread="fields"),
    Route("DELETE", "/api/project/records/{name}", P, "pyazr.Project.remove",
          "Remove a record (refused while another refers to it).", args=("name",)),
    Route("POST", "/api/project/records/{name}/rename", P, "pyazr.Project.rename",
          "Rename a record with every reference to it.", args=("name", "new")),
    Route("GET", "/api/project/levels", P, "pyazr.Project.levels", "Every level of <levels> with its channels."),
    Route("POST", "/api/project/levels", P, "pyazr.Project.add_level",
          "Add a level: J, parity (+1/-1), energy (MeV), channels [{pair, l, s, width, fixed}].",
          args=("J", "parity:int", "energy:float", "channels"), kwargs=("fixed",)),
    Route("PATCH", "/api/project/levels/{label:path}/channels/{pair}/{l}/{s}", P, "pyazr.Project.set_channel",
          "Set a channel's width (the row's convention) or whether it is fixed; s as a decimal (0.5).",
          args=("label", "pair:int", "l:int", "s:float"), kwargs=("width", "fixed")),
    Route("GET", "/api/project/levels/{label:path}", P, "pyazr.Project.level", "One level by its label.",
          args=("label",)),
    Route("PATCH", "/api/project/levels/{label:path}", P, "pyazr.Project.set_level",
          "Set a level's energy (MeV), fixed, active.", args=("label",), kwargs=("energy", "fixed", "active")),
    Route("DELETE", "/api/project/levels/{label:path}", P, "pyazr.Project.remove_level", "Remove a level.",
          args=("label",)),
    Route("GET", "/api/project/pairs", P, "pyazr.Project.pairs", "The particle pairs of <levels>."),
    Route("PATCH", "/api/project/pairs/{key}", P, "pyazr.Project.set_pair",
          "Set a pair's radius (fm), masses (u), separation energy and excitations (MeV).", args=("key:int",),
          kwargs=("radius", "light_mass", "heavy_mass", "separation_energy", "light_excitation",
                  "heavy_excitation")),
    Route("POST", "/api/project/pairs", M, "pyazr.Project.add_pair",
          "Add a pair with its channels (web_ui.md finding 10)."),
    Route("DELETE", "/api/project/pairs/{key}", M, "pyazr.Project.remove_pair",
          "Remove a pair with its channels (web_ui.md finding 10)."),
    Route("GET", "/api/project/reaction", P, "pyazr.Project.reaction",
          "The reaction label from pair entrance to pair exit (query).", args=("entrance:int", "exit:int")),
    Route("GET", "/api/project/config", P, "pyazr.Project.config", "The <config> entries."),
    Route("PATCH", "/api/project/config", P, "pyazr.Project.set_config",
          "Set <config> entries by name (entries: {name: value}; null removes a key).", spread="entries"),
    Route("GET", "/api/project/thm", P, "pyazr.Project.thm", "The <thm> section (options, experiments)."),
    Route("PATCH", "/api/project/thm", P, "pyazr.Project.set_thm", "Set global <thm> options (options: {...}).",
          args=("options",)),
    Route("PUT", "/api/project/thm/experiments/{name}", P, "pyazr.Project.set_thm_experiment",
          "Set keys of a THM experiment (keys: {...}; datasets: names that join it).", args=("name",),
          kwargs=("datasets",), spread="keys"),
    Route("DELETE", "/api/project/thm/experiments/{name}", P, "pyazr.Project.remove_thm_experiment",
          "Remove a THM experiment.", args=("name",)),
    Route("POST", "/api/project/undo", V, "server.Server.undo", "Back to the previous text.", kwargs=("rev",)),
    Route("POST", "/api/project/redo", V, "server.Server.redo", "Forward to the next text.", kwargs=("rev",)),
    Route("GET", "/api/project/history", V, "server.Server.history", "The history: revision, time, what."),
    Route("POST", "/api/project/apply", V, "server.Server.apply",
          "Write the working x (entries x given over it) into the project text; one history entry.",
          kwargs=("rev", "x"), serves=("engine.Engine.snapshot",)),
    # -- the session and the evaluation (derived from the text)
    Route("GET", "/api/session", V, "server.Server.session_describe",
          "The session of the current text: options, data sets, parameters; built when the text changed."),
    Route("GET", "/api/parameters", V, "server.Server.parameters",
          "Every parameter: name, kind, free, value as read, x0 and the working x (fit basis)."),
    Route("GET", "/api/x", V, "server.Server.get_x", "The working vector x by name, and x0."),
    Route("PATCH", "/api/x", V, "server.Server.patch_x",
          "Set entries of x (values: {name: value}, fit basis; reset: back to x0 first); no rebuild.",
          kwargs=("values", "reset")),
    Route("POST", "/api/evaluate", S, "engine.Engine.evaluate",
          "Objective, its terms, data chi2, and the program's chi2 per data set, at x (entries given over the "
          "working x).", kwargs=("x",)),
    Route("GET", "/api/evaluation", S, "engine.Engine.evaluation", "Every data set's arrays at the working x."),
    Route("GET", "/api/evaluation/data/{name}", S, "engine.Engine.data",
          "One data set's arrays at the working x: energies, angles, values, errors, model, residuals.",
          args=("name",)),
    Route("POST", "/api/transform", S, "engine.Engine.transform",
          "Physical values of every parameter at x (entries given over the working x).", kwargs=("x",)),
    Route("POST", "/api/evaluate/grid", M, "pyazr.Session.evaluate",
          "The model on energies and angles given (pyazr3.md 2.3)."),
    Route("GET", "/api/data/{name}/grid", M, "pyazr.DataSet.grid", "The sub-point grid E_n of a data set."),
    Route("GET", "/api/data/{name}/kernel", M, "pyazr.DataSet.kernel", "The kernel rows K_i(E_n) of a data set."),
    Route("POST", "/api/session/level_matrix", M, "pyazr.Session.level_matrix",
          "R and U matrix elements, eigenphases, poles at energies."),
    Route("GET", "/api/physics", V, "server.Server.physics",
          "The physics functions the engine binds: signature and doc (the rows of pyazr.kinds()).",
          serves=("pyazr.kinds",)),
    Route("POST", "/api/physics/{fn}", M, "pyazr.physics.functions", "A physics function on arrays."),
    Route("GET", "/api/system", L, "step 2 (Nucleus tab); today the session's pairs and physical_levels",
          "The built system: pairs, levels, channels, Wigner limits."),
    Route("POST", "/api/session/channel_functions", L, "step 5 (Physics tab); today the session's coulomb_functions",
          "Channel functions of a pair at energies."),
    Route("POST", "/api/outputs", L, "step 3 (Results tab); today the session's write_output_files",
          "Write the program's output files into the project's output directory."),
    Route("POST", "/api/stopping/{fn}", L, "step 5 (Physics tab); today pyazr.stopping", "Stopping powers."),
    Route("GET", "/api/sources/exfor", L, "step 2 (Data tab); today pyazr.search_exfor", "EXFOR search."),
    Route("GET", "/api/sources/ensdf/{nuclide}", L, "step 2 (Nucleus tab); today pyazr.fetch_levels", "ENSDF levels."),
    Route("GET", "/api/files", L, "step 2 (the file browser); no pyazr call", "Files inside the root."),
    Route("GET", "/api/events", L, "step 2 (the page); no pyazr call",
          "Server-sent events: revisions, job progress, build state."),
    # -- jobs
    Route("GET", "/api/jobs/kinds", V, "server.Server.job_kinds",
          "The job table: kind, method, settings, the package or pyazr call each needs."),
    Route("GET", "/api/jobs", V, "server.Server.jobs_list", "Every job of the project (interrupted ones too)."),
    Route("POST", "/api/jobs", V, "server.Server.job_start",
          "Start a job: kind, method, settings, x0 (omitted: the working x; \"text\"; {\"job\": id}; {name: value}).",
          args=("kind", "method"), kwargs=("settings", "x0", "rev")),
    Route("GET", "/api/jobs/{job_id}", V, "server.Server.job_get", "A job: state, progress, best x, result.",
          args=("job_id",)),
    Route("POST", "/api/jobs/{job_id}/stop", V, "server.Server.job_stop",
          "Cooperative stop: ends before the next objective call with the best point.", args=("job_id",)),
    Route("POST", "/api/jobs/{job_id}/kill", V, "server.Server.job_kill", "Kill the job's worker.",
          args=("job_id",)),
    Route("POST", "/api/jobs/{job_id}/accept", V, "server.Server.job_accept",
          "Write the job's best point into the project text (one history entry).", args=("job_id",),
          kwargs=("rev",), serves=("engine.Engine.snapshot",)),
)

# Public calls of pyazr.Project the table does not serve, with the reason (the consistency test needs every public
# name of Project in a row or here).
PROJECT_NOT_SERVED = {
    "session": "the server owns sessions (the session worker); a session of the text, not of the saved file",
    "directory": "a property: the project's directory is in GET /api/project (root)",
    "unsaved": "a property: GET /api/project",
    "add_<kind>": "POST /api/project/records with the kind (Project.add)",
}


def describe_routes():
    """The route table as data (GET /api/describe): method, path, kind, the call, its doc, status."""
    out = []
    for r in ROUTES:
        out.append({"method": r.method, "path": r.path, "kind": r.handler, "call": r.call, "doc": r.doc,
                    "status": {M: "missing", L: "later"}.get(r.handler, "served"),
                    "arguments": list(r.args) + list(r.kwargs) + ([r.spread] if r.spread else [])})
    return out


# -- generic handlers -------------------------------------------------------------------------------------------

def _convert(spec, value):
    name, _, kind = spec.partition(":")
    if value is None:
        return None
    return {"int": int, "float": float}.get(kind, lambda v: v)(value)


def _arguments(route, params, query, body):
    """(positional, keyword) arguments of ``route`` from the request."""
    def find(name):
        for source in (params, query, body):
            if name in source:
                return source[name]
        raise KeyError(f"{route.method} {route.path}: the argument {name!r} is missing")

    args = [_convert(spec, find(spec.split(":")[0])) for spec in route.args]
    kwargs = dict(route.fixed)
    for spec in route.kwargs:
        name = spec.split(":")[0]
        if name in body or name in query:
            kwargs[name] = _convert(spec, body.get(name, query.get(name)))
    if route.spread:
        extra = body.get(route.spread) or {}
        if not isinstance(extra, dict):
            raise TypeError(f"{route.spread}: an object of keyword arguments, not {type(extra).__name__}")
        kwargs.update(extra)
    return args, kwargs


def _rev(query, body, headers):
    for value in (body.get("rev"), query.get("rev"), headers.get("if-match")):
        if value is not None:
            return int(str(value).strip('"'))
    return None


def run(route, server, params, query, body, headers):
    """Serve one request of ``route``: the generic handler of its kind."""
    if route.handler in (M, L):
        what = "the pyazr call" if route.handler == M else "a later step"
        raise Missing(f"{route.method} {route.path} is not served: it needs {what}: {route.call}",
                      missing=route.call, status_kind=route.handler)
    name = route.call.rsplit(".", 1)[1]
    if route.handler == V:
        if "rev" in route.kwargs and "rev" not in body:
            body = dict(body, rev=_rev(query, body, headers))
        args, kwargs = _arguments(route, params, query, body)
        return getattr(server, name)(*args, **kwargs)
    if route.handler == S:
        args, kwargs = _arguments(route, params, query, body)
        keys = [s.split(":")[0] for s in route.args]
        return server.session_request(name, **dict(zip(keys, args)), **kwargs)
    # a method of pyazr.Project
    from pyazr import Project
    method = getattr(Project, name)
    args, kwargs = _arguments(route, params, query, body)
    state = server.state
    if route.method == "GET":
        with state.lock:
            return method(state.project, *args, **kwargs)
    state.check_paths([args, kwargs])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = state.edit(_rev(query, body, headers), f"{name} {' '.join(map(str, args))}".strip(),
                            lambda project: method(project, *args, **kwargs))
    return {"rev": state.rev, "result": result, "warnings": [str(w.message) for w in caught]}


def plain(value):
    """JSON data of a result: records, levels and pairs as objects, numpy as lists, a non-finite float as the
    strings "inf", "-inf", "nan" (JSON has no such numbers)."""
    from pyazr import Record
    if isinstance(value, Record):
        return {"name": value.record_name, "kind": value.record_kind, "line": value.line, "fields": value.fields}
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return plain(dataclasses.asdict(value))
    if isinstance(value, dict):
        if value and all(isinstance(v, Record) for v in value.values()):
            return [plain(v) for v in value.values()]
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if hasattr(value, "tolist"):
        return plain(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return "nan" if math.isnan(value) else ("inf" if value > 0 else "-inf")
    return value


def install(app, server):
    """Add every row of :data:`ROUTES` to ``app``, each served by :func:`run`."""
    from fastapi import Request
    from fastapi.responses import JSONResponse, PlainTextResponse
    from starlette.concurrency import run_in_threadpool

    def endpoint_of(route):
        async def endpoint(request: Request):
            try:
                raw = await request.body()
                body = json.loads(raw) if raw else {}
                if not isinstance(body, dict):
                    raise TypeError("the request body is a JSON object")
                result = await run_in_threadpool(run, route, server, dict(request.path_params),
                                                 dict(request.query_params), body, request.headers)
            except Exception as e:  # noqa: BLE001  -- every exception is a JSON refusal (errors.to_http)
                status, payload = to_http(e)
                return JSONResponse(plain(payload), status_code=status)
            if route.text:
                return PlainTextResponse(result)
            return JSONResponse(plain(result))

        endpoint.__doc__ = f"{route.doc}\n\nCall: {route.call}."
        return endpoint

    for r in ROUTES:
        app.add_api_route(r.path, endpoint_of(r), methods=[r.method], summary=r.doc.split(".")[0],
                          description=f"{r.doc}\n\nCall: `{r.call}` ({r.handler}).",
                          name=f"{r.method} {r.path}", operation_id=f"{r.method.lower()}_{_slug(r.path)}")


def _slug(path):
    return "".join(c if c.isalnum() else "_" for c in path.replace("/api/", "")).strip("_")
