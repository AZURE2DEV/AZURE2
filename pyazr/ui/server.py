"""The web UI's server (docs/upgrade/web_ui.md sections 3 and 8): :class:`Server` (the open project, its session
worker, its jobs, the working parameter vector, the budget), :func:`create_app` (the FastAPI app of the route table
:mod:`pyazr.ui.routes`) and :func:`main` (``python -m pyazr.ui project.azr``).

The server builds no engine session (a build holds the GIL); it holds the project text and its history
(:mod:`pyazr.ui.state`), and asks a worker (:mod:`pyazr.ui.workers`) for every number of the engine.  It binds to a
loopback address only, and every request needs the token (``Authorization: Bearer <token>``, ``?token=`` or the
cookie set by the first ``?token=`` request) and a ``Host`` header naming the loopback (DNS rebinding).
"""
import hmac
import inspect
import os
import secrets
import shutil
import socket
import sys
import time

from .errors import Conflict, NotFound, Refused
from .jobs import JobTable, job_table
from .routes import CONTENT_SECURITY_POLICY
from .state import SESSION_OPTIONS, State, session_options
from .workers import Budget, SessionWorker, Spawner

LOOPBACK = ("127.0.0.1", "localhost", "::1")
TOKEN_COOKIE = "azure2_ui_token"
# The bounds of the first rows of a file (GET /api/files/{path}): rows per request, bytes read of one line.
FILE_HEAD_ROWS = 200
FILE_HEAD_LINE_BYTES = 4096


class Server:
    """Everything one server holds for its open project."""

    def __init__(self, project, *, threads=1, root=None):
        self.budget = Budget(threads)
        self.spawner = Spawner()
        self.started = time.time()
        self._open_state(State(project, root=root))

    def _open_state(self, state):
        self.state = state
        self.worker = SessionWorker(self.spawner, log_path=state.ui_path("logs", "session.log"), cwd=state.directory)
        self.jobs = JobTable(state.ui_path("jobs", "x").rsplit(os.sep, 1)[0], self.budget, self.spawner)
        self.session_text_built = None
        self.free_names, self.x0, self.x, self.build_s, self.session = [], {}, {}, None, None

    def close(self):
        """Kill the jobs' workers, stop the session worker; the working file stays only while edits are unsaved."""
        self.jobs.close()
        self.worker.close()
        if os.path.exists(self.state.session_file()):
            os.remove(self.state.session_file())
        if not self.state.unsaved and self.state.crash_text is None and os.path.exists(self.state.working_file):
            os.remove(self.state.working_file)

    # -- the project ------------------------------------------------------------------------------------------------

    def project_state(self):
        return self.state.describe()

    def text(self):
        return self.state.project.text

    def put_text(self, text, rev=None):
        self.state.replace_text(rev, text)
        return {"rev": self.state.rev}

    def history(self):
        return self.state.history_view()

    def undo(self, rev=None):
        return {"rev": self.state.undo(rev)}

    def redo(self, rev=None):
        return {"rev": self.state.redo(rev)}

    def save(self, rev=None, path=None):
        return {"path": self.state.save(rev, path), "rev": self.state.rev}

    def restore(self, rev=None):
        return {"rev": self.state.restore(rev)}

    def file_head(self, path, rows=20):
        """The first ``rows`` lines (1 to :data:`FILE_HEAD_ROWS`) of a text file inside the root (relative paths from
        the project's directory), as written: the page shows the head of a data file beside its record.  At most
        :data:`FILE_HEAD_LINE_BYTES` bytes of a line are read (a binary file without line ends is not read whole);
        ``cut`` lists the numbers (from 1) of the lines cut there."""
        rows = int(rows)
        if not 1 <= rows <= FILE_HEAD_ROWS:
            raise Refused(f"rows={rows}: the first rows of a file are 1 to {FILE_HEAD_ROWS}",
                          allowed=f"rows from 1 to {FILE_HEAD_ROWS}", fix=f"GET /api/files/{path}?rows=20")
        full = self.state.inside_root(path)
        if not os.path.isfile(full):
            raise NotFound(f"no file {path!r} (relative to {self.state.directory})")
        lines, cut = [], []
        with open(full, "rb") as f:
            while len(lines) < rows:
                raw = f.readline(FILE_HEAD_LINE_BYTES + 1)
                if not raw:
                    return {"path": path, "lines": lines, "more": False, "cut": cut}
                if len(raw) > FILE_HEAD_LINE_BYTES and not raw.endswith(b"\n"):
                    raw = raw[:FILE_HEAD_LINE_BYTES]
                    cut.append(len(lines) + 1)
                    while True:  # skip the rest of the line, a block at a time
                        rest = f.readline(65536)
                        if not rest or rest.endswith(b"\n"):
                            break
                lines.append(raw.decode("utf-8", "replace").rstrip("\r\n"))
            more = bool(f.read(1))
        return {"path": path, "lines": lines, "more": more, "cut": cut}

    def outputs(self, rev=None):
        """Write the program's output files at the working x and copy them into the project's output directory
        (``<config>`` output_dir, inside the root): the user asked for them, so this is the one session request that
        writes into the user's directory.  Returns the directory (relative to the project's) and the files."""
        if rev is not None:
            self.state.check_rev(rev)
        target = self._output_directory()
        with self.worker.lock:
            self.budget.require_session()
            rev = self._ensure_built()
            written = self.worker.call("outputs", x=self._vector())
        os.makedirs(target, exist_ok=True)
        for name in written["files"]:
            if not name or os.path.basename(name) != name or name in (".", ".."):
                raise Refused(f"the worker named {name!r} among the output files: not a plain file name")
            shutil.copy2(os.path.join(written["directory"], name), os.path.join(target, name))
        return {"rev": rev, "directory": os.path.relpath(target, self.state.directory), "files": written["files"],
                "listing": self.outputs_list()["files"]}

    def outputs_list(self):
        """The files of the project's output directory: name, size (bytes), time written (s since the epoch)."""
        target = self._output_directory()
        files = []
        if os.path.isdir(target):
            for name in sorted(os.listdir(target)):
                full = os.path.join(target, name)
                if os.path.isfile(full):
                    st = os.stat(full)
                    files.append({"name": name, "size": st.st_size, "mtime": st.st_mtime})
        return {"directory": os.path.relpath(target, self.state.directory), "files": files}

    def _output_directory(self):
        """The project's output directory (``<config>`` output_dir, default output/), resolved inside the root."""
        with self.state.lock:
            entries = self.state.project.config()
        value = entries.get("output_dir") if isinstance(entries, dict) else None
        return self.state.inside_root(str(value or "output/").strip().strip('"'))

    def open(self, path):
        """Open another project inside the root (no job may run)."""
        self._replace_state(lambda root: State(path, root=root))
        return self.state.describe()

    def new(self, path, nuclear_part):
        """A new project at ``path`` from the nuclear part of the project ``nuclear_part`` (both inside the root);
        unsaved until saved."""
        from pyazr import Project
        source, target = self.state.inside_root(nuclear_part), self.state.inside_root(path)
        self._replace_state(lambda root: State(target, root=root, project=Project.new(target, nuclear_part=source)))
        return self.state.describe()

    def _replace_state(self, make):
        if self.jobs.running():
            job = self.jobs.running()[0].record["id"]
            raise Conflict(f"job {job} is running; stop it before opening another project", job=job)
        state = make(self.state.root)
        self.close()
        self._open_state(state)

    @staticmethod
    def schema():
        import pyazr
        return pyazr.schema()

    @staticmethod
    def kinds():
        import pyazr
        return pyazr.kinds()

    @staticmethod
    def physics():
        import pyazr
        return pyazr.kinds()["physics"]

    @staticmethod
    def physics_call(fn, **arguments):
        from .physics import call  # in the server's process: no session needed (pyazr/ui/physics.py)
        return call(fn, arguments)

    @staticmethod
    def session_functions():
        from .physics import session_rows
        return session_rows()

    def session_function_call(self, fn, **arguments):
        """A session function (``physics.SESSION_ROWS``) on the session of the current text, answered in the shape of
        a physics function's answer, with the Python text that opens a session on the project's file."""
        from .engine import session_text
        from .physics import session_answer, session_request
        row, kwargs, used, lines, texts = session_request(fn, arguments)
        out = self.session_request(row["op"], **kwargs)
        answer = session_answer(fn, row, out, used, lines, texts, session_text(self.state.path))
        return dict(answer, rev=out.get("rev"), warnings=out.get("warnings", []))

    # -- the session and the working vector x ---------------------------------------------------------------------

    def session_request(self, op, **args):
        """``op`` of the adapter (:class:`pyazr.ui.engine.Engine`) on the session of the current text, built first
        when the text changed; an ``x`` argument is the working vector with the entries given ({name: value})."""
        from .engine import Engine
        with self.worker.lock:
            self.budget.require_session()
            notices = []
            rev = self._ensure_built(notices)
            if "x" in inspect.signature(getattr(Engine, op)).parameters:
                args["x"] = self._vector(args.get("x"))
            result = self.worker.call(op, **args)
            notices += self.worker.warnings
        if isinstance(result, dict):
            result = dict(result, rev=rev, warnings=notices)
        return result

    def _ensure_built(self, notices=None):
        """(Re)build the worker's session when the project text differs from the one it was built from.  The working
        x keeps the value of a free parameter whose start value the edit did not change.  The build's warnings (the
        engine's notices) are added to ``notices``."""
        with self.state.lock:
            text, rev = self.state.project.text, self.state.rev
        if self.worker.built_from == text and self.worker.pid() is not None:
            return rev
        session_text = self.state.session_text(text)
        path = self.state.session_file()
        with open(path, "wb") as f:
            f.write(session_text.encode("utf-8", "surrogateescape"))
        t0 = time.perf_counter()
        description = self.worker.call("open", file=path, cwd=self.state.directory, **SESSION_OPTIONS)
        self.build_s = time.perf_counter() - t0
        if notices is not None:
            notices += self.worker.warnings
        self.worker.built_from, self.session_text_built = text, session_text
        names = [p["name"] for p in description["parameters"] if p["free"]]
        x0 = {p["name"]: p["x0"] for p in description["parameters"] if p["free"]}
        self.x = {n: (self.x[n] if n in self.x and self.x0.get(n) == x0[n] else x0[n]) for n in names}
        self.free_names, self.x0, self.session = names, x0, description
        return rev

    def _vector(self, values=None):
        values = values or {}
        self._check_names(values)
        return [float(values.get(n, self.x[n])) for n in self.free_names]

    def _check_names(self, values):
        unknown = sorted(set(values) - set(self.free_names))
        if unknown:
            fixed = [n for n in unknown if any(p["name"] == n for p in self.session["parameters"])]
            if fixed:
                raise Refused(f"{', '.join(fixed)}: fixed in the project; a fixed value is a project edit, not a "
                              "value of x", fix="PATCH the record or level that holds it")
            raise NotFound(f"no free parameter {', '.join(unknown)}", known=self.free_names)

    def _require_session_unless_built(self):
        """Refuse over the budget only when the worker must build: a session of the current text already built is
        read without computing (the parameter table and x stay readable while a job holds the budget)."""
        with self.state.lock:
            current = self.worker.built_from == self.state.project.text and self.worker.pid() is not None
        if not current:
            self.budget.require_session()

    def session_describe(self):
        with self.worker.lock:
            self._require_session_unless_built()
            rev = self._ensure_built()
            return {"rev": rev, "pid": self.worker.pid(), "build_s": self.build_s, "options": self.session["options"],
                    "n_free": self.session["n_free"], "datasets": self.session["datasets"],
                    "parameters": self.session["parameters"], "session_options": session_options()}

    def parameters(self):
        rows = self.session_describe()["parameters"]
        for r in rows:
            r["x"] = self.x.get(r["name"]) if r["free"] else None
        return rows

    def get_x(self):
        with self.worker.lock:
            self._require_session_unless_built()
            rev = self._ensure_built()
            return {"rev": rev, "x": dict(self.x), "x0": dict(self.x0)}

    def patch_x(self, values=None, reset=False):
        """Set entries of the working x ({name: value}, fit basis); ``reset``: back to the text's x0 first."""
        with self.worker.lock:
            self._require_session_unless_built()
            rev = self._ensure_built()
            values = values or {}
            self._check_names(values)
            if reset:
                self.x = dict(self.x0)
            self.x.update({k: float(v) for k, v in values.items()})
            return {"rev": rev, "x": dict(self.x)}

    def apply(self, rev=None, x=None, what="apply x"):
        """Write the working x (with the entries ``x`` given) into the project text: one history entry."""
        self.state.check_rev(rev)
        with self.worker.lock:
            self.budget.require_session()
            self._ensure_built()
            vector = self._vector(x)
            path = self.state.ui_path("snapshots", f"rev{self.state.rev}.azr")
            snapshot = self.worker.call("snapshot", path=path, x=vector)
            os.remove(path)
            merged = self.state.merge_snapshot(self.worker.built_from, self.session_text_built, snapshot)
        self.state.replace_text(rev, merged, what)
        return {"rev": self.state.rev}

    # -- jobs -------------------------------------------------------------------------------------------------------

    @staticmethod
    def job_kinds():
        return job_table()

    def jobs_list(self):
        bulky = ("best", "trace", "log")  # GET /api/jobs/{id} serves them
        return [{k: v for k, v in j.view().items() if k not in bulky} for j in self.jobs.jobs.values()]

    def job_start(self, kind, method, settings=None, x0=None, rev=None):
        """Start a job on the current text from ``x0``: omitted, the working x; ``"text"``, the text's x0;
        ``{"job": id}``, that job's best point (its checkpoint); a {name: value} dict, those entries over the
        working x."""
        if rev is not None:
            self.state.check_rev(rev)
        with self.worker.lock:
            with self.state.lock:
                text, rev = self.state.project.text, self.state.rev
            if isinstance(x0, dict) and set(x0) == {"job"}:
                x0 = self.jobs.best_x(x0["job"])
            elif x0 is None and self.worker.built_from == text:
                x0 = dict(self.x)
            elif x0 == "text":
                x0 = None
            return self.jobs.start(kind, method, settings, x0, rev, self.state.session_text(text),
                                   self.state.directory, self.state.session_file)

    def job_get(self, job_id, wait=None):
        """A job's record, its trace and the last lines of its log; with ``wait`` (s), once its state changes or
        after that long (a long poll: the page follows a running job with one request in flight)."""
        return (self.jobs.wait(job_id, wait) if wait else self.jobs.get(job_id)).view()

    def job_stop(self, job_id):
        return self.jobs.stop(job_id)

    def job_kill(self, job_id):
        return self.jobs.kill(job_id)

    def job_accept(self, job_id, rev=None):
        """Write the job's best point into the project text (one history entry); the job's record keeps the revision
        that made and the chi-squared, objective and N at the best point on the job's own text (``accepted``), so a
        client compares them with the project as written, also after a reload."""
        job = self.jobs.get(job_id)
        best = self.jobs.best_x(job_id)
        at_best = None
        if job.record.get("text_rev") == self.state.rev:  # the fit's own numbers, on the text it ran on
            e = self.session_request("evaluate", x=best)
            at_best = {"chi2": e["chi2"], "objective": e["objective"], "n": e["totals"].get("Total-N")}
        answer = self.apply(rev, best, what=f"accept job {job_id}")
        with job.lock:
            job.record["accepted"] = {"rev": answer["rev"], "time": time.time(), "at_best": at_best}
        return answer

    # -- the server itself ------------------------------------------------------------------------------------------

    def describe(self):
        import pyazr
        from .routes import describe_routes
        return {"server": "pyazr.ui", "package": pyazr.describe(), "python": sys.version.split()[0],
                "budget": self.budget.describe(), "project": self.state.describe(), "jobs": job_table(),
                "session_options": session_options(), "routes": describe_routes()}

    def processes(self):
        """PIDs and memory (kB, from /proc: resident and its peak) of the server and its workers."""
        out = [dict(role="server", pid=os.getpid(), **_memory(os.getpid()))]
        if self.worker.pid():
            out.append(dict(role="session", pid=self.worker.pid(), **_memory(self.worker.pid())))
        for j in self.jobs.jobs.values():
            if j.process is not None and j.process.alive():
                out.append(dict(role=f"job {j.record['id']}", pid=j.process.pid, **_memory(j.process.pid)))
        return out


def _memory(pid):
    try:
        with open(f"/proc/{pid}/status") as f:
            fields = dict(line.split(":", 1) for line in f if ":" in line)
        return {"rss_kb": int(fields["VmRSS"].split()[0]), "peak_kb": int(fields["VmHWM"].split()[0])}
    except (OSError, KeyError, ValueError):
        return {"rss_kb": None, "peak_kb": None}


# -- the app ------------------------------------------------------------------------------------------------------

def create_app(project, *, threads=1, token, root=None):
    """The FastAPI app serving ``project`` (a format-3 .azr) with the compute budget ``threads`` and the access
    ``token`` (required: every request must carry it).  The workers stop with the app (its lifespan)."""
    from contextlib import asynccontextmanager

    from fastapi import FastAPI, Request
    from fastapi.responses import JSONResponse

    from .routes import install

    if not token:
        raise ValueError("a token is required (main() makes one when none is given)")
    server = Server(project, threads=threads, root=root)

    @asynccontextmanager
    async def lifespan(app):
        yield
        server.close()

    app = FastAPI(title="AZURE2 web UI API", version="0.1", lifespan=lifespan,
                  description="The API of the web UI of AZURE2 on pyazr: GET /api/describe lists every route.")
    app.state.server = server

    @app.middleware("http")
    async def guard(request: Request, call_next):
        host_header = request.headers.get("host") or ""
        host = host_name(host_header)
        if host not in LOOPBACK:
            return JSONResponse({"error": "Forbidden", "message": f"Host {host!r} is not a loopback name",
                                 "allowed": list(LOOPBACK)}, status_code=403)
        origin = request.headers.get("origin")
        if origin is not None and origin != f"{request.url.scheme}://{host_header}":
            # every port of localhost is one site for the browser: a page of another local server would carry the
            # cookie; only the server's own page (same scheme, host and port) may send a request with an Origin
            return JSONResponse({"error": "Forbidden", "message": f"Origin {origin!r} is not this server's",
                                 "allowed": f"{request.url.scheme}://{host_header}"}, status_code=403)
        given = _token_of(request)
        if given is None or not hmac.compare_digest(given, token):
            return JSONResponse({"error": "Unauthorized", "message": "the token is missing or wrong",
                                 "fix": "Authorization: Bearer <token>, or ?token=<token> (the URL the server "
                                        "printed)"}, status_code=401)
        response = await call_next(request)
        if request.query_params.get("token"):
            response.set_cookie(TOKEN_COOKIE, token, httponly=True, samesite="strict")
        response.headers["Content-Security-Policy"] = CONTENT_SECURITY_POLICY
        return response

    install(app, server)
    return app


def host_name(header):
    """The host of a Host header without its port: ``[::1]:8000`` and ``[::1]`` give ``::1``, ``localhost:80`` gives
    ``localhost``."""
    header = header.strip()
    if header.startswith("["):
        return header[1:header.index("]")] if "]" in header else header
    return header.rsplit(":", 1)[0] if header.count(":") == 1 else header


def _token_of(request):
    auth = request.headers.get("authorization", "")
    if auth.lower().startswith("bearer "):
        return auth[7:].strip()
    return request.query_params.get("token") or request.cookies.get(TOKEN_COOKIE)


def free_port():
    """A free TCP port of 127.0.0.1 above 20000."""
    for _ in range(100):
        port = secrets.randbelow(40000) + 20001
        with socket.socket() as s:
            try:
                s.bind(("127.0.0.1", port))
                return port
            except OSError:
                continue
    raise RuntimeError("no free port found above 20000")


def main(project, port=0, threads=1, token=None, *, host="127.0.0.1", root=None):
    """Serve ``project`` on ``host:port`` (loopback only; port 0: a free one above 20000) with the compute budget
    ``threads`` and ``token`` (a random one when None).  Prints the URL with its token on one line and writes it to
    ``<ui dir>/url`` (mode 0600); runs until interrupted, then stops every worker."""
    if host not in LOOPBACK:
        raise SystemExit(f"pyazr.ui: refusing to bind to {host!r}: the server binds to a loopback address only "
                         f"({', '.join(LOOPBACK)}); use an SSH tunnel for remote use")
    import uvicorn
    token = token or secrets.token_urlsafe(24)
    port = int(port) or free_port()
    app = create_app(project, threads=threads, token=token, root=root)
    url = f"http://{'[::1]' if host == '::1' else host}:{port}/?token={token}"
    url_file = app.state.server.state.ui_path("url")
    fd = os.open(url_file, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as f:
        f.write(url + "\n")
    print(url, flush=True)
    uvicorn.Server(uvicorn.Config(app, host=host, port=port, log_level="warning")).run()
