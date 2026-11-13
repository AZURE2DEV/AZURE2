"""Jobs of the web UI (docs/upgrade/web_ui.md section 4.4): the job table, the recipes it runs, and the job records.

A job kind is one row of :data:`JOBS`: its kind and method, the function that runs it (or the pyazr call it waits
for: the route answers 501 naming it), the package it needs, its settings with their defaults, and its doc.  A job
runs in a worker process of its own (``python -m pyazr.ui.workers job``) on its own session of the text of the
revision it was started on; it reports progress after every objective call (calls, objective, best objective, elapsed
time, best x by name), stops cooperatively before the next call when asked, and writes a checkpoint
``<ui dir>/jobs/<id>.json`` (state, progress, best x) on every improvement and at its end, from which a new job can
start (``x0: {"job": id}``) and which a restarted server lists as ``interrupted``.

A row's settings are rows of :class:`Setting` (name, type, default, doc): the page draws the job form from them and the
server converts a request's values by the type's row of :data:`SETTING_TYPES`.  The server keeps a decimated trace of
the objective against the call number (at most :data:`TRACE_POINTS` points) and serves the last lines of the job's
log beside its record; ``GET /api/jobs/{id}?wait=s`` answers when the job's state changes or after ``s`` seconds.
"""
import json
import os
import threading
import time
import uuid
from dataclasses import dataclass
from functools import lru_cache, partial

import numpy as np

from .errors import Conflict, Missing, NotFound, Refused, UIError, describe_exception

FINAL = ("stopped", "converged", "ended", "failed", "killed", "interrupted")
TRACE_POINTS = 500  # the most points of a job's trace (objective against call); older points are thinned by two
LOG_LINES = 8  # the last lines of a job's log served with its record
WAIT_MAX = 10.0  # s, the longest a GET of a job waits for its state to change


class StopRequested(Exception):
    """Raised inside the objective when a stop was asked: the recipe ends with the best point seen."""


class Tracker:
    """The job's handle for its function: wraps the session's objective-and-gradient for a recipe (reads the stop
    request before each call, counts it, records the best point and reports it with ``report(message)``; ``None`` in
    a direct call), and counts the progress of calls the job makes itself (:meth:`count`: the engine's fitters, the
    sampler).  ``session`` is the job's session adapter (:class:`pyazr.ui.engine.Engine`), for a job that calls the
    engine other than through the objective."""

    def __init__(self, fg, names, *, report=None, stop=None, session=None):
        self.fg, self.names, self.report, self.session = fg, names, report, session
        self.stop = stop if stop is not None else threading.Event()  # set: the job ends with its best point
        self.calls, self.best, self.best_x, self.t0 = 0, np.inf, None, time.perf_counter()
        self.first = None  # the objective at the start point
        self.after_stop = 0  # progress reports after the stop request was seen (the calls running then)
        self._last = None  # (x bytes, f, g) of the last call, for recipes that ask f and g separately

    def __call__(self, x):
        """(f, g) at x."""
        x = np.asarray(x, float)
        key = x.tobytes()
        if self._last is not None and self._last[0] == key:
            return self._last[1], self._last[2]
        if self.stop_requested():
            raise StopRequested()
        f, g = self.fg(x)
        f = float(f)
        self._last = (key, f, np.asarray(g, float))
        self.count(self.calls + 1, f, x)
        return f, self._last[2]

    def count(self, calls, objective, x, **more):
        """Progress: ``calls`` so far (counting up), the objective at the free vector ``x`` of the last of them, and
        further numbers by name (``more``: a fitter's stage, a sampler's acceptance; None where there is none),
        reported with the best point when it improved."""
        objective = float(objective)
        if self.stop.is_set():
            self.after_stop += 1
        self.calls = int(calls)
        if self.first is None:
            self.first = objective
        improved = objective < self.best
        if improved:
            self.best, self.best_x = objective, np.array(x, float)
        if self.report is not None:
            self.report({"event": "call", "calls": self.calls, "objective": objective, "best": self.best,
                         "elapsed": time.perf_counter() - self.t0,
                         "best_x": dict(zip(self.names, self.best_x.tolist())) if improved else None, **more})

    def stop_requested(self):
        return self.stop.is_set()

    def value(self, x):
        return self(x)[0]

    def grad(self, x):
        return self(x)[1]


# -- the recipes: functions of (tracker, x0, settings, names) -> (status, message, extra) -------------------------

def fit_scipy_lbfgsb(tracker, x0, settings, names):
    """scipy.optimize.minimize, L-BFGS-B, on the engine's objective and gradient (pyazr3.md 2.6)."""
    import scipy.optimize as so
    r = so.minimize(tracker, x0, jac=True, method="L-BFGS-B", options={"maxiter": int(settings["max_iterations"])})
    return ("converged" if r.success else "ended"), str(r.message), {"iterations": int(r.nit)}


def fit_iminuit_migrad(tracker, x0, settings, names):
    """iminuit's MIGRAD with the engine's gradient, then HESSE (errors, covariance and correlation by name), errordef 1
    (a chi-squared)."""
    m = _minuit(tracker, x0, names)
    m.migrad(ncall=settings["max_calls"])
    extra = {"edm": float(m.fmin.edm), "valid": bool(m.valid)}
    if not m.valid:
        return "ended", "MIGRAD: invalid minimum", extra
    if not settings["hesse"]:
        return "converged", "MIGRAD: valid minimum", extra
    m.hesse()  # a valid minimum can turn invalid here (the EDM recomputed with HESSE's matrix)
    extra.update(_errors_of(m, names), edm=float(m.fmin.edm), valid=bool(m.valid))
    extra["errors_from"] = "HESSE (iminuit)"
    if m.valid:
        return "converged", "MIGRAD: valid minimum; HESSE: still valid", extra
    return "ended", f"MIGRAD: valid minimum; after HESSE invalid (EDM {m.fmin.edm:.3g})", extra


def minos_iminuit(tracker, x0, settings, names):
    """iminuit's MINOS from a fit's best point: MIGRAD (a few calls at a minimum), HESSE, then MINOS on the parameters
    named (all free ones when none): asymmetric intervals {name: {lower, upper, valid}} in the fit basis."""
    m = _minuit(tracker, x0, names)
    m.migrad(ncall=settings["max_calls"])
    extra = {"edm": float(m.fmin.edm), "valid": bool(m.valid)}
    if not m.valid:
        return "ended", "MIGRAD from the start point: invalid minimum; MINOS not run", extra
    m.hesse()
    extra.update(_errors_of(m, names), errors_from="HESSE (iminuit)")
    which = settings["parameters"] or list(names)
    m.minos(*which, ncall=settings["max_calls"])
    extra["minos"] = {n: {"lower": float(m.merrors[n].lower), "upper": float(m.merrors[n].upper),
                          "valid": bool(m.merrors[n].is_valid)} for n in which}
    bad = [n for n in which if not m.merrors[n].is_valid]
    message = "MINOS: every interval valid" if not bad else f"MINOS: invalid intervals of {', '.join(bad)}"
    return ("converged" if not bad else "ended"), message, extra


def _minuit(tracker, x0, names):
    import iminuit
    m = iminuit.Minuit(tracker.value, np.asarray(x0, float), grad=tracker.grad, name=names)
    m.errordef = iminuit.Minuit.LEAST_SQUARES
    return m


def _errors_of(m, names):
    """HESSE's errors by name, and its covariance and correlation matrices in the order of ``names``."""
    return {"errors": dict(zip(names, [float(e) for e in m.errors])),
            "covariance": {"names": list(names), "matrix": np.asarray(m.covariance, float).tolist()},
            "correlation": {"names": list(names), "matrix": np.asarray(m.covariance.correlation(), float).tolist()}}


def _covariance_by_name(names, matrix):
    """The covariance and the correlation ({names, matrix}) of a covariance matrix in the order of ``names``."""
    c = np.asarray(matrix, float)
    sd = np.sqrt(np.clip(np.diag(c), 0.0, None))
    scale = np.where(sd > 0, sd, 1.0)
    return {"covariance": {"names": list(names), "matrix": c.tolist()},
            "correlation": {"names": list(names), "matrix": (c / np.outer(scale, scale)).tolist()}}


# The job's state after an engine fit, by the engine's status (any other status: "ended", with the status said).
ENGINE_STATES = {"converged": "converged", "stopped": "stopped"}


def fit_engine(tracker, x0, settings, names, *, method, minos=None):
    """The engine's fitter ``method`` (``Session.fit_with_engine``, the program's fit) from x0: each progress record
    counted into the tracker (the calls of the fit, then of HESSE and MINOS, counted on), the tracker's stop request
    passed as the fit's; MINOS when ``minos`` (or the setting ``minos``).  The job's best point is the fitter's end
    point and its objective the engine's there; the start objective is one evaluation at x0 after the fit (an
    evaluation before it would change the state the engine's pooled working copies carry into the fit).  The
    fitter's lines go to the job's log."""
    from .engine import fit_with_engine, objective_of
    stage = {"name": None, "offset": 0, "last": 0}

    def progress(record):
        if record["stage"] != stage["name"]:  # HESSE and MINOS count their calls from 1: the job counts on
            stage.update(name=record["stage"], offset=stage["offset"] + stage["last"])
        stage["last"] = record["calls"]
        measure = record["measure"]
        tracker.count(stage["offset"] + record["calls"], record["objective"], record["x"], stage=record["stage"],
                      measure=float(measure) if np.isfinite(measure) else None)

    run_minos = settings.get("minos", False) if minos is None else minos
    r = fit_with_engine(tracker.session, method, run_minos, x0, settings["max_calls"], progress, tracker.stop)
    print(r["log"], flush=True)
    tracker.first = float(objective_of(tracker.session)(x0))
    tracker.best, tracker.best_x = r["objective"], np.array([r["x"][n] for n in names], float)
    words = f"the engine's {r['method']}: {r['status']}"
    if r["fell_back_to_migrad"]:
        words += " (no analytic Jacobian for a data set: it ran MIGRAD, without HESSE or MINOS)"
    extra = {"engine_status": r["status"], "valid": r["converged"], "edm": r["edm"], "fit_calls": r["calls"]}
    if r["covariance"] is not None:
        extra.update(_covariance_by_name(r["names"], r["covariance"]), errors=r["errors"],
                     errors_from="HESSE after the engine's MIGRAD" if r["covariance_from_hesse"]
                     else "the engine's (J^T J + priors)^-1")
    if r["minos"]:
        extra["minos"] = {n: {"lower": lo, "upper": hi, "valid": None} for n, (lo, hi) in r["minos"].items()}
    return ENGINE_STATES.get(r["status"], "ended"), words, extra


# -- settings: their types and the conversion of a request's value ------------------------------------------------

def _as_integer(value):
    if value is None or value == "":
        return None
    number = float(value)
    if not number.is_integer():
        raise ValueError(f"{value!r} is not an integer")
    return int(number)


def _as_boolean(value):
    if isinstance(value, bool):
        return value
    if str(value).lower() in ("true", "1", "yes"):
        return True
    if str(value).lower() in ("false", "0", "no"):
        return False
    raise ValueError(f"{value!r} is not true or false")


def _as_names(value):
    if value is None or value == "":
        return []
    words = value.replace(",", " ").split() if isinstance(value, str) else list(value)
    return [str(w) for w in words]


# The value types of a job's settings: type -> (conversion of a request's value, what the type allows).  The page has
# one widget per type (tabs/fit.js, SETTING_WIDGETS); a new type is a row here and a widget there.
SETTING_TYPES = {
    "integer": (_as_integer, "an integer, or empty for the default"),
    "boolean": (_as_boolean, "true or false"),
    "names": (_as_names, "names of free parameters, separated by blanks or commas; empty: every free parameter"),
}


@dataclass(frozen=True)
class Setting:
    """One setting of a job kind: its name, type (a key of :data:`SETTING_TYPES`), default and doc."""
    name: str
    type: str
    default: object
    doc: str


@dataclass(frozen=True)
class JobRow:
    """One kind of job: ``kind`` and ``method`` name it; ``run`` its function ``(tracker, x0, settings, names) ->
    (state, message, extra)``, or ``missing`` the pyazr call it waits for; ``needs`` the package it imports;
    ``settings`` its :class:`Setting` rows; ``follows`` the (kind, method) of the job whose best point it starts from
    (the page offers it beside that job's result); ``facts`` what the page shows of it beside its doc (name -> value:
    the program's switch, whether it gives a covariance, whether MINOS can follow)."""
    kind: str
    method: str
    doc: str
    run: object = None
    needs: str = None
    settings: tuple = ()
    missing: str = None
    follows: tuple = None
    facts: tuple = ()  # (name, value) pairs

    def defaults(self):
        return {s.name: s.default for s in self.settings}


MAX_CALLS_ENGINE = Setting("max_calls", "integer", None, "the most objective calls of the fit (the program's "
                           "--max-calls; empty: no limit)")

# The rows written here: the recipes on the objective, and the jobs that wait for a pyazr call.  The engine's fitters
# are rows made from its table (engine_rows); job_rows() is the whole table.
JOBS = (
    JobRow("fit", "iminuit_migrad", "iminuit's MIGRAD with the gradient, then HESSE (errors, covariance)",
           run=fit_iminuit_migrad, needs="iminuit",
           settings=(Setting("max_calls", "integer", None, "the most objective calls of MIGRAD (empty: iminuit's "
                             "default)"),
                     Setting("hesse", "boolean", True, "HESSE after a valid minimum: errors and covariance")),
           facts=(("covariance", True), ("MINOS can follow", True))),
    JobRow("fit", "scipy_lbfgsb", "scipy's L-BFGS-B on the objective and its gradient (no errors)",
           run=fit_scipy_lbfgsb, needs="scipy",
           settings=(Setting("max_iterations", "integer", 1000, "the most iterations of L-BFGS-B"),),
           facts=(("covariance", False), ("MINOS can follow", False))),
    JobRow("minos", "iminuit", "iminuit's MINOS intervals from the best point of an iminuit fit",
           run=minos_iminuit, needs="iminuit", follows=("fit", "iminuit_migrad"),
           settings=(Setting("parameters", "names", [], "the parameters MINOS scans (empty: every free one)"),
                     Setting("max_calls", "integer", None, "the most objective calls of each MIGRAD and MINOS scan "
                             "(empty: iminuit's default)"))),
    JobRow("sample", "emcee", "the emcee sampler on the log-probability", missing="Objective.log_prob (pyazr3.md 2.6)"),
    JobRow("bands", "covariance", "uncertainty bands of a fit",
           missing="pyazr.bands on a Session and an Evaluation (pyazr3.md 2.2)"),
    JobRow("rate", "engine", "the reaction rate", missing="Session.reaction_rate (web_ui.md finding 6)"),
)


def engine_rows(fitters):
    """The job rows of the engine's fitter table (``pyazr.kinds()["fitters"]``: name -> row): a fit job per fitter
    (the program's default first, then the table's order), and a MINOS job after each fitter whose row says MINOS can
    follow it.  A new fitter of the engine is a new row here with no line of the UI changed."""
    order = sorted(fitters, key=lambda name: fitters[name]["option"] != "")
    fits, minos = [], []
    for name in order:
        f = fitters[name]
        method = f"engine_{name}"
        settings = (MAX_CALLS_ENGINE,)
        if f["minos"]:
            settings += (Setting("minos", "boolean", False, "MINOS after the fit and its HESSE: asymmetric "
                                 "intervals of every free parameter (Up = the session's chiVariance)"),)
        facts = (("program switch", f["option"] or "none (the program's default)"), ("covariance", f["covariance"]),
                 ("MINOS can follow", f["minos"]), ("gradient", f["gradient"]), ("residual rows", f["residual_rows"]))
        fits.append(JobRow("fit", method, f"the engine's {name}: {f['doc']}", settings=settings, facts=facts,
                           run=partial(fit_engine, method=name)))
        if f["minos"]:
            minos.append(JobRow("minos", method, f"the engine's MINOS after its {name} (the fit from the start point, "
                                "HESSE, MINOS on every free parameter)", settings=(MAX_CALLS_ENGINE,),
                                follows=("fit", method), run=partial(fit_engine, method=name, minos=True)))
    return tuple(fits), tuple(minos)


@lru_cache(maxsize=1)
def job_rows():
    """The whole job table: the engine's fit rows, the rows of :data:`JOBS`, the engine's MINOS rows."""
    import pyazr
    fits, minos = engine_rows(pyazr.kinds()["fitters"])
    return fits + JOBS + minos


def job_row(kind, method):
    rows = job_rows()
    for row in rows:
        if row.kind == kind and row.method == method:
            return row
    known = [f"{r.kind}/{r.method}" for r in rows]
    raise NotFound(f"no job {kind}/{method} (the job table: {', '.join(known)})", known=known)


def job_table():
    """The job table as data (GET /api/jobs/kinds)."""
    out = []
    for r in job_rows():
        available = r.missing is None and _importable(r.needs)
        settings = [{"name": x.name, "type": x.type, "default": x.default, "doc": x.doc,
                     "allowed": SETTING_TYPES[x.type][1]} for x in r.settings]
        out.append({"kind": r.kind, "method": r.method, "doc": r.doc, "settings": settings, "needs": r.needs,
                    "available": available, "missing": r.missing, "follows": list(r.follows) if r.follows else None,
                    "facts": dict(r.facts)})
    return out


def merged_settings(row, given):
    """The row's defaults with the request's values over them, each converted by its type; Refused on an unknown
    setting or a value its type does not allow."""
    by_name = {x.name: x for x in row.settings}
    unknown = sorted(set(given or {}) - set(by_name))
    if unknown:
        raise Refused(f"job {row.kind}/{row.method} has no setting {', '.join(unknown)}", allowed=sorted(by_name))
    out = row.defaults()
    for name, value in (given or {}).items():
        convert, allowed = SETTING_TYPES[by_name[name].type]
        try:
            out[name] = convert(value)
        except (TypeError, ValueError) as e:
            raise Refused(f"job {row.kind}/{row.method}, setting {name}: {e}", allowed=allowed) from None
    return out


def _importable(name):
    if name is None:
        return True
    import importlib.util
    return importlib.util.find_spec(name) is not None


# -- worker side ------------------------------------------------------------------------------------------------------

def write_checkpoint(path, record):
    """Write ``record`` as JSON to ``path`` atomically (a reader never sees half a file)."""
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(record, f)
    os.replace(tmp, path)


def plain_json(value):
    """``value`` with every non-finite float as None (JSON has no NaN or infinity; a browser refuses them)."""
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, dict):
        return {k: plain_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain_json(v) for v in value]
    return value


def listen_for_stop(channel, stop):
    """Set ``stop`` when the server asks the job to stop, or closes its input (a thread of the job's worker: the
    request is seen while the engine computes, which releases the GIL)."""
    def listen():
        while True:
            message = channel.receive()
            if message is None or message.get("op") == "stop":
                stop.set()
                return

    threading.Thread(target=listen, name="pyazr-ui-job-stop", daemon=True).start()


def run_job(channel):
    """A job worker: read the job's spec, build the session, run the row's function, report, write the checkpoint."""
    from .engine import objective_of, open_for_job
    spec = channel.receive()
    record = {"id": spec["id"], "kind": spec["kind"], "method": spec["method"], "text_rev": spec["rev"],
              "settings": spec["settings"], "state": "building", "progress": None, "best": None, "result": None}
    checkpoint = spec["checkpoint"]
    stop = threading.Event()
    listen_for_stop(channel, stop)

    def report(message):
        message = plain_json(message)
        channel.send(message)
        record["progress"] = {k: v for k, v in message.items() if k not in ("event", "best_x")}
        if message["best_x"] is not None:
            record["best"] = {"objective": message["best"], "x": message["best_x"]}
            write_checkpoint(checkpoint, record)

    try:
        write_checkpoint(checkpoint, record)
        t0 = time.perf_counter()
        session, x0, names = open_for_job(spec["session_file"], spec["cwd"])
        channel.send({"event": "built", "build_s": time.perf_counter() - t0, "names": names})
        if spec.get("x0"):
            x0 = vector_from_names(spec["x0"], names, x0)
        record["state"] = "running"
        write_checkpoint(checkpoint, record)
        row = job_row(spec["kind"], spec["method"])
        tracker = Tracker(objective_of(session).value_and_grad, names, report=report, stop=stop, session=session)
        try:
            status, message, extra = row.run(tracker, x0, spec["settings"], names)
        except StopRequested:
            status, message, extra = "stopped", "stopped on request", {}
        best_x = None if tracker.best_x is None else dict(zip(names, tracker.best_x.tolist()))
        record.update(state=status, best={"objective": tracker.best, "x": best_x} if best_x else None,
                      result={"status": status, "message": message, "calls": tracker.calls,
                              "start_objective": tracker.first, "start_x": dict(zip(names, np.asarray(x0, float)
                                                                                    .tolist())),
                              "elapsed": time.perf_counter() - tracker.t0, **extra})
        if stop.is_set():
            record["result"]["reports_after_stop"] = tracker.after_stop
    except Exception as e:  # noqa: BLE001  -- reported as the job's failure
        classes, text, _ = describe_exception(e)
        cls = classes[0]
        record.update(state="failed", result={"status": "failed", "error": cls, "message": text})
    record = plain_json(record)
    write_checkpoint(checkpoint, record)
    channel.send({"event": "done", "record": record})


def vector_from_names(values, names, x0):
    """x with the entries named in ``values`` ({name: value}); every name must be a free parameter of ``names``."""
    unknown = sorted(set(values) - set(names))
    if unknown:
        raise KeyError(f"not free parameters of this session: {', '.join(unknown[:8])} (free: {', '.join(names)})")
    x = np.array(x0, float)
    for i, n in enumerate(names):
        if n in values:
            x[i] = float(values[n])
    return x


# -- server side ------------------------------------------------------------------------------------------------------

class Job:
    """A job as the server sees it: its record (the checkpoint's fields) and its worker."""

    def __init__(self, record, process=None, log_path=None):
        self.record, self.process, self.session_file, self.log_path = record, process, None, log_path
        self.lock = threading.Lock()
        self.stop_asked_at = None
        self.trace_stride = 1  # a call enters the trace when its number is a multiple of this

    def add_to_trace(self, calls, objective, best):
        """Keep (call, objective, best) in the record's trace: every call until TRACE_POINTS, then every second one
        of those kept and every second call after (the trace stays a picture of the whole run).  Under the lock."""
        trace = self.record.setdefault("trace", [])
        if calls % self.trace_stride:
            return
        trace.append([calls, objective, best])
        if len(trace) > TRACE_POINTS:
            self.trace_stride *= 2
            trace[:] = [p for p in trace if p[0] % self.trace_stride == 0]

    def view(self):
        with self.lock:
            out = dict(self.record)
            out["pid"] = self.process.pid if self.process is not None else None
        out["log"] = log_tail(self.log_path)
        return out


def log_tail(path, lines=LOG_LINES):
    """The last ``lines`` non-empty lines of a job's log (the engine's notes; a carriage return ends a line too)."""
    if not path or not os.path.exists(path):
        return []
    with open(path, "rb") as f:
        f.seek(max(0, os.path.getsize(path) - 8192))
        text = f.read().decode("utf-8", "replace")
    return [line for line in text.replace("\r", "\n").split("\n") if line.strip()][-lines:]


class JobTable:
    """The jobs of the open project: started, followed by a reader thread each, stopped, killed, listed; those
    found in the jobs directory at start-up (a server that died) are ``interrupted``."""

    def __init__(self, directory, budget, spawner):
        self.directory, self.budget, self.spawner = directory, budget, spawner
        self.jobs = {}
        os.makedirs(directory, exist_ok=True)
        for name in sorted(os.listdir(directory)):
            if name.endswith(".json"):
                with open(os.path.join(directory, name)) as f:
                    record = json.load(f)
                if record.get("state") not in FINAL:
                    record["state"] = "interrupted"
                self.jobs[record["id"]] = Job(record, log_path=os.path.join(directory, f"{record['id']}.log"))

    def start(self, kind, method, settings, x0, rev, session_text, cwd, file_of):
        """Start a job of the table's row (kind, method) on ``session_text`` (the session file's text of revision
        ``rev``, written to ``file_of(job id)``) from ``x0`` ({name: value} or None for the text's x0)."""
        from .workers import WorkerProcess
        row = job_row(kind, method)
        if row.missing:
            raise Missing(f"job {kind}/{method} needs a pyazr call that does not exist yet: {row.missing}",
                          missing=row.missing)
        if not _importable(row.needs):
            raise Missing(f"job {kind}/{method} needs the package {row.needs}, not installed", missing=row.needs)
        merged = merged_settings(row, settings)
        job_id = uuid.uuid4().hex[:8]
        threads = self.budget.take_job(job_id)
        try:
            session_file = file_of(job_id)
            with open(session_file, "wb") as f:
                f.write(session_text.encode("utf-8", "surrogateescape"))
            spec = {"id": job_id, "kind": kind, "method": method, "settings": merged, "x0": x0, "rev": rev,
                    "session_file": session_file, "cwd": cwd, "checkpoint": os.path.join(self.directory,
                                                                                         f"{job_id}.json")}
            process = WorkerProcess(self.spawner, "job", threads=threads,
                                    log_path=os.path.join(self.directory, f"{job_id}.log"), cwd=cwd)
            process.send(spec)
        except Exception:
            self.budget.give_back(job_id)
            raise
        job = Job({"id": job_id, "kind": kind, "method": method, "text_rev": rev, "settings": merged,
                   "state": "building", "progress": None, "best": None, "result": None,
                   "started": time.time(), "trace": []}, process, log_path=process.log_path)
        job.session_file = session_file
        self.jobs[job_id] = job
        threading.Thread(target=self._follow, args=(job,), name=f"pyazr-ui-job-{job_id}", daemon=True).start()
        return job.view()

    def _follow(self, job):
        """Read the job's messages until it ends; then end its worker, give its thread back to the budget, and only
        then show its final state (a client that sees the state final finds the budget free)."""
        final = None
        try:
            while True:
                message = job.process.receive()
                if message["event"] == "done":
                    final = dict(message["record"])
                    break
                with job.lock:
                    if message["event"] == "built":
                        job.record["state"] = "stopping" if job.stop_asked_at else "running"
                        job.record["build_s"] = message["build_s"]
                    elif message["event"] == "call":
                        job.record["progress"] = {k: v for k, v in message.items() if k not in ("event", "best_x")}
                        job.add_to_trace(message["calls"], message["objective"], message["best"])
                        if message["best_x"] is not None:
                            job.record["best"] = {"objective": message["best"], "x": message["best_x"]}
        except UIError:  # the worker ended without "done": killed, or crashed
            with job.lock:
                final = dict(job.record)
            final["state"] = "killed" if final.get("kill_asked") else "failed"
            if final["state"] == "failed":
                final["result"] = {"status": "failed", "message": "the worker ended without a result",
                                   "log": job.process.log_path}
            write_checkpoint(os.path.join(self.directory, f"{final['id']}.json"), final)
        finally:
            job.process.kill()
            self.budget.give_back(job.record["id"])
            if job.session_file and os.path.exists(job.session_file):
                os.remove(job.session_file)
            with job.lock:
                job.record.update(final or {"state": "failed"})
                job.record["ended"] = time.time()

    def wait(self, job_id, seconds):
        """The job once its state differs from the one it has now, or after ``seconds`` (at most WAIT_MAX): a client
        follows a running job with one request in flight and sees a stop, a kill or the end at once."""
        job = self.get(job_id)
        state, end = job.record["state"], time.monotonic() + min(max(float(seconds), 0.0), WAIT_MAX)
        if state in FINAL:
            return job
        while job.record["state"] == state and time.monotonic() < end:
            time.sleep(0.02)
        return job

    def get(self, job_id):
        if job_id not in self.jobs:
            raise NotFound(f"no job {job_id!r} (the jobs: {', '.join(self.jobs) or 'none'})", known=list(self.jobs))
        return self.jobs[job_id]

    def stop(self, job_id):
        """Cooperative stop: the job ends before its next objective call, with the best point seen."""
        job = self.get(job_id)
        with job.lock:
            if job.record["state"] in FINAL:
                return dict(job.record)
            job.stop_asked_at = time.time()
            job.record["state"] = "stopping"
        try:
            job.process.send({"op": "stop"})
        except UIError:
            pass
        return job.view()

    def kill(self, job_id):
        """SIGKILL the job's worker; its state becomes ``killed`` (the best point is the last reported)."""
        job = self.get(job_id)
        with job.lock:
            if job.record["state"] in FINAL:
                return dict(job.record)
            job.record["kill_asked"] = True
        if job.process is not None:
            job.process.proc.kill()
        return job.view()

    def running(self):
        return [j for j in self.jobs.values() if j.record["state"] not in FINAL]

    def best_x(self, job_id):
        """{name: value} of the job's best point (its checkpoint when the server restarted)."""
        job = self.get(job_id)
        with job.lock:
            best = job.record.get("best")
        if not best or not best.get("x"):
            raise Conflict(f"job {job_id} has no best point yet", fix=f"GET /api/jobs/{job_id}")
        return best["x"]

    def close(self):
        """Kill every running job's worker (server shutdown)."""
        for job in list(self.jobs.values()):
            if job.process is not None and job.process.alive():
                job.record["kill_asked"] = True
                job.process.kill()
