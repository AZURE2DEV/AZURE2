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

The page's files are not API rows: :data:`STATIC_ROUTES` lists them, served by one generic static handler.

An argument spec is ``name`` or ``name:int`` / ``name:float`` (converted); positional ``args`` are looked up in the
path, the query, then the body; ``kwargs`` in the body (or query) and passed only when given; ``spread`` names the
body key whose object is passed as keyword arguments (a record's fields); ``fixed`` adds constant keywords.
The consistency test (tests/pyazr/ui/routes_consistency_test.py) holds this table, the app and pyazr together.
"""
import dataclasses
import json
import math
import os
import warnings
from dataclasses import dataclass, field

from .errors import Missing, to_http
from .physics import MAX_ELEMENTS


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
    Route("GET", "/api/evaluation", S, "engine.Engine.evaluation",
          "The Results payload at the working x: the program's chi-squared split (as /api/evaluate), every data "
          "set's arrays with its chi-squared row, every extrapolation on its grid, the model data, and the residual's "
          "definition (residual)."),
    Route("GET", "/api/evaluation/data/{name}", S, "engine.Engine.data",
          "One data set's arrays at the working x: c.m. and excitation energies, c.m. angles (after energy shifts), "
          "the datum and its error as read and times the norm in use, the model, the engine's standardized "
          "residuals, the S factor where offered.", args=("name",)),
    Route("POST", "/api/transform", S, "engine.Engine.transform",
          "Physical values of every parameter at x (entries given over the working x).", kwargs=("x",)),
    Route("POST", "/api/evaluate/grid", M, "pyazr.Session.evaluate.energies",
          "The model on energies and angles given (pyazr3.md 2.3)."),
    Route("GET", "/api/data/{name}/grid", M, "pyazr.DataSet.grid", "The sub-point grid E_n of a data set."),
    Route("GET", "/api/data/{name}/kernel", M, "pyazr.DataSet.kernel", "The kernel rows K_i(E_n) of a data set."),
    Route("POST", "/api/session/level_matrix", M, "pyazr.Session.level_matrix",
          "R and U matrix elements, eigenphases, poles at energies."),
    Route("GET", "/api/physics", V, "server.Server.physics",
          "The physics functions, one row each (pyazr.kinds()['physics']): arguments with unit, kind and domain, "
          "returned fields with units, the plot axis.", serves=("pyazr.kinds",)),
    Route("POST", "/api/physics/{fn}", V, "server.Server.physics_call",
          "pyazr.physics.<fn> on {\"arguments\": {name: number, list or axis {from, to, points, spacing}}} "
          "(broadcast), errors='mask': the returned fields as lists, the status per element and its words, the "
          f"arguments as used, the Python text of the call; at most {MAX_ELEMENTS} elements per request; runs in the "
          "server's process (no session).", args=("fn",), spread="arguments", serves=("pyazr.physics.functions",)),
    Route("GET", "/api/system", L, "step 2 (Nucleus tab); today the session's pairs and its levels",
          "The built system: pairs, levels, channels, Wigner limits."),
    Route("POST", "/api/session/channel_functions", S, "engine.Engine.channel_functions",
          "The session's channel functions of a pair (key) at c.m. energies e_cm (MeV): F, Fp, G, Gp, P, S, delta_hs, "
          f"open; l, radius (fm, 0: the channel radius), J; at most {MAX_ELEMENTS} energies per request.",
          args=("pair:int", "e_cm"), kwargs=("l:int", "radius:float", "J:float")),
    Route("POST", "/api/session/thm_experiments", S, "engine.Engine.thm_experiments",
          "Every THM experiment as a chi-squared evaluation at x profiles it (entries given over the working x): by "
          "name, its data sets, background, points, chi2, status, norm, background coefficients, covariance, sigma "
          "and the coherent background's parameters.", kwargs=("x",)),
    Route("POST", "/api/session/thm_background", S, "engine.Engine.thm_background",
          "The profiled norm and background of one THM experiment at x, with their covariance.",
          args=("experiment",), kwargs=("x",)),
    Route("POST", "/api/session/thm_lineshape", S, "engine.Engine.thm_lineshape",
          "The Coulomb line shape of a THM experiment (lineshape=on) at c.m. energies (MeV) of x + A, the level "
          f"poles of x: per exit pair zeta, eta_sb and per level |N_C|^2; at most {MAX_ELEMENTS} energies.",
          args=("experiment", "energies"), kwargs=("x",)),
    Route("POST", "/api/session/thm_vertex", S, "engine.Engine.thm_vertex",
          "The THM entrance vertex of an experiment at c.m. energies (MeV) of x + A: per channel and level "
          f"<|M_l|^2> over the spectator window and at p_s = 0, the window's nodes, reached; at most {MAX_ELEMENTS} "
          "energies.", args=("experiment", "energies"), kwargs=("x",)),
    Route("POST", "/api/session/thm_distortion", S, "engine.Engine.thm_distortion",
          "The distortion factor R(E) of a THM experiment at c.m. energies (MeV) of x + A, with |M|^2, |M_PW|^2 and "
          f"the spectator-F channel; at most {MAX_ELEMENTS} energies.", args=("experiment", "energies")),
    Route("GET", "/api/thm/keys", V, "server.Server.thm_keys",
          "The keys of <thm> (the options, an experiment's) as record kinds of a schema: value type, choices, unit, "
          "doc (pyazr/ui/thm.py KEYS: the engine describes them by its grammar only)."),
    Route("GET", "/api/session/functions", V, "server.Server.session_functions",
          "The session's functions as rows of the shape of /api/physics (pyazr/ui/physics.py SESSION_ROWS and the THM "
          "diagnostics of pyazr/ui/thm.py: the session's methods have no table in the binding); a row with `tab` is "
          "drawn by that tab."),
    Route("POST", "/api/session/functions/{fn}", V, "server.Server.session_function_call",
          "A session function on {\"arguments\": {...}} (one value each, the axis a list or {from, to, points, "
          "spacing}) in the session worker, answered as /api/physics/{fn} answers (returns, units, status, its words, "
          "the arguments as used, the Python text; a THM call's facts and tables).", args=("fn",), spread="arguments",
          serves=("engine.Engine.channel_functions", "engine.Engine.thm_experiments", "engine.Engine.thm_background",
                  "engine.Engine.thm_lineshape", "engine.Engine.thm_vertex", "engine.Engine.thm_distortion")),
    Route("POST", "/api/outputs", V, "server.Server.outputs",
          "Write the program's output files at the working x (pyazr's write_output_files) and copy them into the "
          "project's output directory (<config> output_dir): the files written.", kwargs=("rev",),
          serves=("engine.Engine.outputs",)),
    Route("GET", "/api/outputs", V, "server.Server.outputs_list",
          "The files of the project's output directory: name, size (bytes), time written (s since the epoch)."),
    Route("GET", "/api/stopping", V, "server.Server.stopping",
          "The stopping functions, one row each (pyazr.kinds()['stopping']): the shape of /api/physics's rows, a text "
          "argument formula first, and SRIM's citation, owed by every use.", serves=("pyazr.kinds",)),
    Route("POST", "/api/stopping/{fn}", V, "server.Server.stopping_call",
          "pyazr.stopping.evaluate(fn, formula, ...) on {\"arguments\": {formula: text, name: number, list or axis}}, "
          "answered as /api/physics/{fn} answers, with the citation; at most "
          f"{MAX_ELEMENTS} elements per request; in the server's process.", args=("fn",), spread="arguments",
          serves=("pyazr.stopping.evaluate",)),
    Route("GET", "/api/sources/exfor", L, "step 2 (Data tab); today pyazr.search_exfor", "EXFOR search."),
    Route("GET", "/api/sources/ensdf/{nuclide}", L, "step 2 (Nucleus tab); today pyazr.fetch_levels", "ENSDF levels."),
    Route("GET", "/api/files", L, "step 2 (the file browser); no pyazr call", "Files inside the root."),
    Route("GET", "/api/files/{path:path}", V, "server.Server.file_head",
          "The first rows (query: rows, 1 to 200, default 20) of a text file inside the root, as written; a line is "
          "read to 4096 bytes (cut lines listed in cut).", args=("path",),
          kwargs=("rows:int",)),
    Route("GET", "/api/events", L, "a later step; the page polls the revision and long-polls a job "
          "(GET /api/jobs/{job_id}?wait=)", "Server-sent events: revisions, job progress, build state."),
    # -- jobs
    Route("GET", "/api/jobs/kinds", V, "server.Server.job_kinds",
          "The job table: kind, method, settings, the package or pyazr call each needs."),
    Route("GET", "/api/jobs", V, "server.Server.jobs_list", "Every job of the project (interrupted ones too)."),
    Route("POST", "/api/jobs", V, "server.Server.job_start",
          "Start a job: kind, method, settings, x0 (omitted: the working x; \"text\"; {\"job\": id}; {name: value}).",
          args=("kind", "method"), kwargs=("settings", "x0", "rev")),
    Route("GET", "/api/jobs/{job_id}", V, "server.Server.job_get",
          "A job: state, progress, best x, result, its trace (objective against call) and the last lines of its log; "
          "with wait (s, at most 10), once its state changes or after that long.", args=("job_id",),
          kwargs=("wait:float",)),
    Route("POST", "/api/jobs/{job_id}/stop", V, "server.Server.job_stop",
          "Cooperative stop: ends before the next objective call with the best point.", args=("job_id",)),
    Route("POST", "/api/jobs/{job_id}/kill", V, "server.Server.job_kill", "Kill the job's worker.",
          args=("job_id",)),
    Route("POST", "/api/jobs/{job_id}/accept", V, "server.Server.job_accept",
          "Write the job's best point into the project text (one history entry).", args=("job_id",),
          kwargs=("rev",), serves=("engine.Engine.snapshot",)),
)

# The page's static files (docs/upgrade/web_ui.md section 7), served by the generic static handler from
# pyazr/ui/static: (path, the file under static/ or None for the path's own, doc).  Kept apart from ROUTES, which holds
# the API (every row of ROUTES is a pyazr call or a server method; these are files).
STATIC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")
STATIC_ROUTES = (
    ("/", "index.html", "The page; a request with ?token= is answered by a redirect to / that sets the cookie."),
    ("/static/{path:path}", None, "The page's modules, style sheet and the vendored plotly.js."),
)
STATIC_TYPES = {".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8",
                ".css": "text/css; charset=utf-8", ".svg": "image/svg+xml", ".txt": "text/plain; charset=utf-8"}

# The Content-Security-Policy of every response.  Scripts only from the server ('self': the page's modules and the
# vendored plotly.js; no inline script, no eval); styles from the server and inline, because plotly.js sets style
# attributes and adds <style> elements while it draws; images also as data: and blob: URLs (plotly's PNG and SVG
# export); requests only to the server itself; no frames, forms, objects or base URL.
CONTENT_SECURITY_POLICY = ("default-src 'none'; script-src 'self'; style-src 'self' 'unsafe-inline'; "
                           "img-src 'self' data: blob:; connect-src 'self'; font-src 'self'; base-uri 'none'; "
                           "form-action 'none'; frame-ancestors 'none'; object-src 'none'")


def static_file(path):
    """(the absolute file, its media type) of ``path`` under :data:`STATIC_DIR`; NotFound outside it or absent."""
    from .errors import NotFound
    full = os.path.realpath(os.path.join(STATIC_DIR, path))
    if not full.startswith(os.path.realpath(STATIC_DIR) + os.sep) or not os.path.isfile(full):
        raise NotFound(f"no static file {path!r}")
    media = STATIC_TYPES.get(os.path.splitext(full)[1])
    if media is None:
        raise NotFound(f"{path!r}: not a file type the page serves ({', '.join(STATIC_TYPES)})")
    return full, media


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
    for path, fixed, doc in STATIC_ROUTES:
        app.add_api_route(path, static_endpoint_of(fixed), methods=["GET"], summary=doc.split(";")[0],
                          include_in_schema=False, name=f"static {path}")


def static_endpoint_of(fixed):
    """The generic static handler: the file ``fixed`` (or the request's ``path``) of :data:`STATIC_DIR`, revalidated
    on every load (the page and the server change together); ``?token=`` on the page is answered by a redirect to the
    same path without it (the server's guard sets the cookie on that answer, and the token leaves the address bar)."""
    from fastapi import Request
    from fastapi.responses import FileResponse, JSONResponse, RedirectResponse

    async def endpoint(request: Request):
        if fixed is not None and "token" in request.query_params:
            return RedirectResponse(request.url.path, status_code=303)
        try:
            full, media = static_file(fixed or request.path_params["path"])
        except Exception as e:  # noqa: BLE001  -- a refusal is JSON, as for the API
            status, payload = to_http(e)
            return JSONResponse(plain(payload), status_code=status)
        return FileResponse(full, media_type=media, headers={"Cache-Control": "no-cache"})

    return endpoint


def _slug(path):
    return "".join(c if c.isalnum() else "_" for c in path.replace("/api/", "")).strip("_")
