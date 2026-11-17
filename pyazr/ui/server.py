"""The web UI's server (docs/upgrade/web_ui.md sections 3 and 8): :class:`Server` (the open project, its session
worker, its jobs, the working parameter vector, the budget), :func:`create_app` (the FastAPI app of the route table
:mod:`pyazr.ui.routes`) and :func:`main` (``python -m pyazr.ui project.azr``).

The server builds no engine session (a build holds the GIL); it holds the project text and its history
(:mod:`pyazr.ui.state`), and asks a worker (:mod:`pyazr.ui.workers`) for every number of the engine.  It binds to a
loopback address only, and every request needs the token (``Authorization: Bearer <token>``, ``?token=`` or the
cookie set by the first ``?token=`` request; a browser without it gets a page that asks for it) and a ``Host``
header naming the loopback (DNS rebinding).
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

FILES_MAX = 2000  # the files GET /api/files lists at most
FILES_DEPTH = 4  # how deep below the project's directory it looks
GLANCE_BYTES = 4096  # what it reads of a file to say what it holds


def _inside(path, root):
    full = os.path.realpath(path)
    return os.path.commonpath([full, root]) == root


def _glance(full):
    """What a file holds at a glance (its first bytes): ``{"kind": project | table | text | binary}``, a table with
    the number of numbers on its first line that is not a comment (``#``, ``!``) or empty."""
    if full.endswith(".azr"):
        return {"kind": "project"}
    with open(full, "rb") as f:
        head = f.read(GLANCE_BYTES)
    if b"\0" in head:
        return {"kind": "binary"}
    for line in head.decode("utf-8", "replace").splitlines()[:-1] or head.decode("utf-8", "replace").splitlines():
        words = line.replace(",", " ").split()
        if not words or words[0][0] in "#!":
            continue
        try:
            [float(w) for w in words]
        except ValueError:
            return {"kind": "text"}
        return {"kind": "table", "columns": len(words)}
    return {"kind": "text"}
from .jobs import JobTable, job_table
from .routes import CONTENT_SECURITY_POLICY
from .state import State, session_options
from . import evaluations, prior, utilities
from .workers import Budget, SessionWorker, Spawner

LOOPBACK = ("127.0.0.1", "localhost", "::1")
# The cookie's name is TOKEN_COOKIE_<port> (cookie_name): a browser keeps one cookie jar per host whatever the port,
# so two servers on one machine would otherwise overwrite each other's token.
TOKEN_COOKIE = "azure2_ui_token"
# The policy of the page that asks for the token (token_page): its own inline style and its form, nothing else.
TOKEN_PAGE_SECURITY_POLICY = ("default-src 'none'; style-src 'unsafe-inline'; form-action 'self'; base-uri 'none'; "
                              "frame-ancestors 'none'")
# The bounds of the first rows of a file (GET /api/files/{path}): rows per request, bytes read of one line.
FILE_HEAD_ROWS = 200
FILE_HEAD_LINE_BYTES = 4096


class Server:
    """Everything one server holds for its open project."""

    def __init__(self, project, *, threads=1, root=None):
        self.budget = Budget(threads)
        self.spawner = Spawner()
        self.started = time.time()
        self._open_state(self._first_state(project, root))

    @staticmethod
    def _first_state(project, root):
        """The open project: the file at ``project``, or, where it does not exist and its folder does, a new empty
        project there (``pyazr.Project.new(path)``; written when saved) with the output and checks directories the
        engine writes into."""
        path = os.path.abspath(project)
        if os.path.exists(path) or not path.endswith(".azr") or not os.path.isdir(os.path.dirname(path)):
            return State(project, root=root)
        from pyazr import Project
        for directory in ("output", "checks"):
            os.makedirs(os.path.join(os.path.dirname(path), directory), exist_ok=True)
        return State(path, root=root, project=Project.new(path))

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

    def files(self, outputs=False):
        """The files of the project's directory (and of its subdirectories, to :data:`FILES_DEPTH` levels; not the
        program's ``output`` and ``checks``, nor hidden ones), at most :data:`FILES_MAX`: each with its path relative
        to the project (as a record's field names it), its size in bytes and what it holds at a glance (``kind``:
        ``project`` for a ``.azr``, ``table`` with ``columns`` when its first line that is not a comment is numbers,
        ``text``, ``binary``).  Whether a data set can read it is the engine's answer: :meth:`field_preview`.
        ``outputs`` lists the program's ``output`` too (a fit's ``param.sav`` and ``covariance.dat``, a sampler's
        chain: what the parts of the prior builder name)."""
        with self.state.lock:
            top = self.state.directory
            hidden = {"checks"} if outputs else {"output", "checks"}
            out, more = [], False
            for here, dirs, names in os.walk(top):
                depth = os.path.relpath(here, top).count(os.sep) + (here != top)
                dirs[:] = sorted(d for d in dirs if not d.startswith(".") and not (here == top and d in hidden)
                                 and depth < FILES_DEPTH and _inside(os.path.join(here, d), self.state.root))
                for name in sorted(names):
                    full = os.path.join(here, name)
                    if name.startswith(".") or not os.path.isfile(full) or not _inside(full, self.state.root):
                        continue
                    if len(out) == FILES_MAX:
                        more = True
                        break
                    out.append(dict(path=os.path.relpath(full, top).replace(os.sep, "/"),
                                    size=os.path.getsize(full), **_glance(full)))
                if more:
                    break
            return {"directory": top, "files": out, "more": more}

    def field_preview(self, name, field, value):
        """What the engine makes of the record ``name`` with its field ``field`` set to ``value`` (a file, say), on a
        copy of the project (the project is not changed): the record's description (rows read and kept, their energy
        range, the reader's error) or the engine's refusal with its three parts."""
        from pyazr import Project
        with self.state.lock:
            copy = Project(self.state.project.to_bytes(), self.state.path)
        self.state.check_paths([value])
        copy.set(name, **{field: value})
        try:
            return {"name": name, "field": field, "value": value, "description": copy.describe(name, as_json=True)}
        except ValueError as e:  # a record the engine does not describe: the edit is accepted, nothing to read
            return {"name": name, "field": field, "value": value, "description": None, "note": str(e)}

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
        """Open another project inside the root (no job may run); a relative ``path`` is relative to the open
        project's folder (as every path of the API), not to the server's working directory."""
        full = self.state.inside_root(path)
        self._replace_state(lambda root: State(full, root=root))
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
    def masses(heavy, light=None):
        """A pair's charges, atomic masses (u) and separation energy (MeV) from the AME2020 table
        (``pyazr.sources.masses.pair_masses``: ``Z1``, ``M1``, ``Z2``, ``M2``, ``S``, ``compound``); without ``light``
        (a photon pair) the heavy partner's ``Z2`` and ``M2`` alone.  Names as ``"p"``, ``"12C"``, ``"4He"``."""
        from pyazr.sources import masses
        if not light:
            entry = masses.lookup(heavy)
            return {"Z2": entry.Z, "M2": entry.atomic_mass, "table": masses.AME_VERSION}
        return dict(masses.pair_masses(light, heavy)._asdict(), table=masses.AME_VERSION)

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
    def stopping():
        import pyazr
        return pyazr.kinds()["stopping"]

    @staticmethod
    def stopping_call(fn, **arguments):
        from .physics import call  # in the server's process, as the physics functions
        return call(fn, arguments, table="stopping")

    @staticmethod
    def thm_keys():
        from .thm import keys
        return keys()

    @staticmethod
    def session_functions():
        from .physics import session_rows
        return session_rows()

    def session_function_call(self, fn, **arguments):
        """A session function (pyazr's rows, ``Session.describe()["functions"]``) on the session of the current text
        at the working x, answered in the shape of a physics function's answer, with the Python text that opens a
        session on the project's file."""
        from .engine import session_text
        from .physics import session_answer, session_request
        row, kwargs, used, lines, texts = session_request(fn, arguments)
        out = self.session_request("session_function", fn=fn, arguments=kwargs)
        own = {k: v for k, v in out.items() if k not in ("rev", "warnings")}  # what session_request adds
        answer = session_answer(fn, row, own, used, lines, texts, session_text(self.state.path))
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
        description = self.worker.call("open", file=path, cwd=self.state.directory)
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
                    "pairs": self.session["pairs"],
                    "parameters": self.session["parameters"],
                    "session_options": session_options(session=self.session)}

    def parameters(self):
        rows = [dict(r) for r in self.session_describe()["parameters"]]  # the session's description stays as built
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

    # -- utilities: the covariance builder and the rate's tables (pyazr/ui/utilities.py) ------------------------------

    def covariance_open(self, name):
        return utilities.covariance_open(self.state, name)

    def covariance_build(self, datasets, parts):
        return utilities.covariance_build(self.state, datasets, parts)

    def covariance_write(self, datasets, parts, path, rev=None):
        return utilities.covariance_write(self.state, datasets, parts, path, rev)

    @staticmethod
    def rate_grids():
        return utilities.rate_grids()

    # -- ENDF-6 evaluations: the files, a summary, the curves over the data, the codes (pyazr/ui/evaluations.py) ----

    def endf_files(self):
        """The ENDF-6 files of the project folder (relative paths), the folder, and how a file is recognised: the
        files of :meth:`files` whose first two lines read as ENDF-6 lines."""
        listing = self.files()
        return {"directory": listing["directory"],
                "files": evaluations.endf_files(listing["directory"], listing["files"]),
                "depth": FILES_DEPTH, "more": listing["more"],
                "how": "files whose first two lines read as ENDF-6 lines (80 columns with MAT, MF, MT), in the "
                       f"project's directory and {FILES_DEPTH} levels below; a path inside the root "
                       "may also be typed"}

    def endf_summary(self, path):
        """What an ENDF-6 file inside the root holds (``pyazr.endf.read(path).summary()``), in the server's process."""
        return evaluations.summary(self.state.inside_root(path), path)

    def endf_curves(self, files):
        """The files' values at every data set's points of the current text at the working x, with the ratios to the
        datum times its normalization and to the model (the session worker: ``Material.curves(data, session=s)``)."""
        names = [files] if isinstance(files, str) else list(files or [])
        if not names:
            raise Refused("name at least one ENDF-6 file", allowed="files of GET /api/endf/files, or paths inside "
                          "the root", fix='POST /api/endf/curves {"files": ["a.endf"]}')
        resolved = [[n, self.state.inside_root(n)] for n in names]
        return self.session_request("endf_curves", files=resolved)

    @staticmethod
    def endf_checkers():
        return evaluations.checkers()

    def endf_accept(self, job_id, to=None, overwrite=False):
        """Keep the candidate of a finished rewrite job: at ``to`` (inside the root; none: beside the original), an
        existing file replaced only with ``overwrite``.  No history entry: an ENDF file is not the project text."""
        job = self.jobs.get(job_id)
        answer = evaluations.accept(self.state, job.view(), to, overwrite)
        with job.lock:
            job.record.setdefault("kept", []).append(dict(answer, time=time.time()))
        return answer

    def endf_resonances(self, path, isotope=0, energy_range=0):
        """What MF=2 of an ENDF-6 file inside the root holds for a new project (``Material.resonances``), nothing
        written; in the server's process."""
        return evaluations.resonances(self.state.inside_root(path), path, isotope, energy_range)

    def endf_keep_project(self, job_id, to=None):
        """Put the new project of a finished ``endf_import`` job into the project folder (a new file only); the open
        project is not changed."""
        job = self.jobs.get(job_id)
        answer = evaluations.keep_project(self.state, job.view(), to)
        with job.lock:
            job.record.setdefault("kept", []).append(dict(answer, time=time.time()))
        return answer

    def endf_export(self):
        """What an ENDF-6 file of the current text at the working x would lose (the session worker:
        ``pyazr.endf.from_session``), the writer's default energy range, and a file name not in use for it."""
        answer = self.session_request("endf_export")
        stem = os.path.splitext(os.path.basename(self.state.path))[0]
        name, k = f"{stem}.endf", 1
        while os.path.exists(os.path.join(self.state.directory, name)):
            k += 1
            name = f"{stem}_{k}.endf"
        return dict(answer, name=name)

    # -- the prior builder (pyazr/ui/prior.py) ----------------------------------------------------------------------

    def prior_open(self):
        return prior.prior_open(self.state, [f["path"] for f in self.files(outputs=True)["files"]])

    def prior_build(self, parts):
        return prior.prior_build(self.state, self.session_describe(), parts)

    def prior_write(self, parts, path, rev=None):
        return prior.prior_write(self.state, self.session_describe(), parts, path, rev)

    # -- jobs -------------------------------------------------------------------------------------------------------

    @staticmethod
    def job_kinds():
        return job_table()

    def jobs_list(self):
        bulky = ("best", "trace", "log")  # GET /api/jobs/{id} serves them, and a result's bands
        out = []
        for j in self.jobs.jobs.values():
            record = {k: v for k, v in j.view().items() if k not in bulky}
            if (record.get("result") or {}).get("bands"):
                record["result"] = dict(record["result"], bands="GET /api/jobs/{id}")
            out.append(record)
        return out

    def job_start(self, kind, method, settings=None, x0=None, rev=None):
        """Start a job on the current text from ``x0``: omitted, the working x; ``"text"``, the text's x0;
        ``{"job": id}``, that job's best point (its checkpoint), and that job's record handed to the new job as its
        source (a covariance to propagate, a chain to continue); a {name: value} dict, those entries over the working
        x.  Files a job writes (a chain) go to the project's output directory."""
        if rev is not None:
            self.state.check_rev(rev)
        from .jobs import job_row
        for name in job_row(kind, method).files:  # a file a job reads: inside the root (403 otherwise)
            if (settings or {}).get(name):
                self.state.inside_root(settings[name])
        source = None
        with self.worker.lock:
            with self.state.lock:
                text, rev = self.state.project.text, self.state.rev
                products = os.path.join(self.state.directory, self.state.project.config().get("output_dir") or "")
            if isinstance(x0, dict) and set(x0) == {"job"}:
                source = {k: v for k, v in self.jobs.get(x0["job"]).view().items() if k not in ("trace", "log")}
                x0 = self.jobs.best_x(x0["job"])
            elif x0 is None and self.worker.built_from == text:
                x0 = dict(self.x)
            elif x0 == "text":
                x0 = None
            return self.jobs.start(kind, method, settings, x0, rev, self.state.session_text(text),
                                   self.state.directory, self.state.session_file, source=source, products=products)

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
        from .first_view import first_view_columns
        return {"server": "pyazr.ui", "package": pyazr.describe(), "python": sys.version.split()[0],
                "budget": self.budget.describe(), "project": self.state.describe(), "jobs": job_table(),
                "session_options": self._project_options(), "routes": describe_routes(),
                "first_view": first_view_columns()}

    def _project_options(self):
        with self.state.lock:
            return session_options(project=self.state.project)

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
    from fastapi.responses import HTMLResponse, JSONResponse

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
            if asks_for_a_page(request):
                return HTMLResponse(token_page(wrong=given is not None), status_code=401,
                                    headers={"Content-Security-Policy": TOKEN_PAGE_SECURITY_POLICY,
                                             "Cache-Control": "no-store"})
            return JSONResponse({"error": "Unauthorized", "message": "the token is missing or wrong",
                                 "fix": "Authorization: Bearer <token>, or ?token=<token> (the URL the server "
                                        "printed)"}, status_code=401)
        response = await call_next(request)
        if request.query_params.get("token"):
            # "lax", not "strict": the printed URL is often followed from another site or application (a link in a
            # chat or a mail page), and the browser keeps a strict cookie back on that navigation and on the redirect
            # that follows it, so the page answered 401 to its own link.  A lax cookie goes with a navigation (GET)
            # from elsewhere and with nothing else from there; every GET route only reads, and a request that
            # changes something carries an Origin header, checked above.
            response.set_cookie(cookie_name(host_header), token, httponly=True, samesite="lax")
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


def cookie_name(host_header):
    """The token cookie's name for the port the browser names in its Host header (the port it reached the server on,
    a tunnel's local port included): ``azure2_ui_token_8000``; without a port, the scheme's default 80."""
    header = host_header.strip()
    port = header.rsplit(":", 1)[1] if ":" in header and not header.endswith("]") else "80"
    return f"{TOKEN_COOKIE}_{port}" if port.isdigit() else f"{TOKEN_COOKIE}_80"


def asks_for_a_page(request):
    """A browser's navigation (a GET that accepts HTML, outside /api): its refusal is a page, not JSON."""
    return (request.method == "GET" and not request.url.path.startswith("/api/")
            and "text/html" in request.headers.get("accept", ""))


def token_page(*, wrong):
    """The page a browser gets instead of the JSON refusal when the token is missing or wrong: what happened, where
    the token is, and a field that takes it (a GET of ``/?token=...``: the same door as the printed URL)."""
    what = "The token in this address is not this server's." if wrong else "This address was opened without its token."
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light dark"><title>AZURE2 — token needed</title>
<style>body {{ font: 16px/1.5 system-ui, sans-serif; max-width: 40em; margin: 3em auto; padding: 0 1em; }}
input {{ font: inherit; width: 22em; max-width: 100%; }} button {{ font: inherit; }}</style></head>
<body><h1>AZURE2: the token is needed</h1>
<p id="token-refusal">{what}</p>
<p>Open the address the server printed when it started, with its <code>?token=...</code> part (it is also in the
file <code>.&lt;project&gt;.ui/url</code> beside the project), or paste the token here:</p>
<form method="get" action="/"><input name="token" type="password" autocomplete="off" autofocus
aria-label="token"> <button type="submit">Open the page</button></form>
</body></html>
"""


def _token_of(request):
    auth = request.headers.get("authorization", "")
    if auth.lower().startswith("bearer "):
        return auth[7:].strip()
    return (request.query_params.get("token")
            or request.cookies.get(cookie_name(request.headers.get("host") or "")))


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


def main(project, port=0, threads=1, token=None, *, host="127.0.0.1", root=None, log_requests=False):
    """Serve ``project`` on ``host:port`` (loopback only; port 0: a free one above 20000) with the compute budget
    ``threads`` and ``token`` (a random one when None).  Prints the URL with its token on one line and writes it to
    ``<ui dir>/url`` (mode 0600); runs until interrupted, then stops every worker.  ``log_requests``: every request
    on standard error (uvicorn's access log; the browser test keeps it beside a failure)."""
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
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write(url + "\n")
    print(url, flush=True)
    config = uvicorn.Config(app, host=host, port=port, log_level="info" if log_requests else "warning",
                            access_log=log_requests, log_config=None if log_requests else uvicorn.config.LOGGING_CONFIG)
    if log_requests:  # the access log on standard error with the time of each line (standard output holds the URL)
        import logging
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter("%(asctime)s.%(msecs)03d %(name)s %(message)s", "%H:%M:%S"))
        logger = logging.getLogger("uvicorn")  # uvicorn.access and uvicorn.error reach it
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
    uvicorn.Server(config).run()
