"""Jobs of the web UI (docs/upgrade/web_ui.md section 4.4): the job table, the recipes it runs, and the job records.

A job kind is one row of :data:`JOBS`: its kind and method, the function that runs it (or the pyazr call it waits
for: the route answers 501 naming it), the package it needs, its settings with their defaults, and its doc.  A job
runs in a worker process of its own (``python -m pyazr.ui.workers job``) on its own session of the text of the
revision it was started on; it reports progress after every objective call (calls, objective, best objective, elapsed
time, best x by name), stops cooperatively before the next call when asked, and writes a checkpoint
``<ui dir>/jobs/<id>.json`` (state, progress, best x) on every improvement and at its end, from which a new job can
start (``x0: {"job": id}``) and which a restarted server lists as ``interrupted``.
"""
import json
import os
import threading
import time
import uuid
from dataclasses import dataclass, field

import numpy as np

from .errors import Conflict, Missing, NotFound, Refused, UIError, describe_exception

FINAL = ("stopped", "converged", "ended", "failed", "killed", "interrupted")


class StopRequested(Exception):
    """Raised inside the objective when a stop was asked: the recipe ends with the best point seen."""


class Tracker:
    """Wraps an objective-and-gradient for a recipe: reads the stop request before each call, records the best
    point, and reports each call (``report(message)``; ``None`` in a direct call)."""

    def __init__(self, fg, names, *, report=None, stop_requested=lambda: False):
        self.fg, self.names, self.report, self.stop_requested = fg, names, report, stop_requested
        self.calls, self.best, self.best_x, self.t0 = 0, np.inf, None, time.perf_counter()
        self.first = None  # the objective at the start point
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
        self.calls += 1
        if self.first is None:
            self.first = f
        improved = f < self.best
        if improved:
            self.best, self.best_x = f, x.copy()
        self._last = (key, f, np.asarray(g, float))
        if self.report is not None:
            self.report({"event": "call", "calls": self.calls, "objective": f, "best": self.best,
                         "elapsed": time.perf_counter() - self.t0,
                         "best_x": dict(zip(self.names, self.best_x.tolist())) if improved else None})
        return f, self._last[2]

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
    """iminuit's MIGRAD with the engine's gradient, then HESSE (errors by name), errordef 1 (a chi-squared)."""
    import iminuit
    m = iminuit.Minuit(tracker.value, np.asarray(x0, float), grad=tracker.grad, name=names)
    m.errordef = iminuit.Minuit.LEAST_SQUARES
    m.migrad(ncall=settings["max_calls"])
    extra = {"edm": float(m.fmin.edm), "valid": bool(m.valid)}
    if settings["hesse"] and m.valid:
        m.hesse()
        extra["errors"] = dict(zip(names, [float(e) for e in m.errors]))
    return ("converged", "MIGRAD: valid minimum", extra) if m.valid else ("ended", "MIGRAD: invalid minimum", extra)


@dataclass(frozen=True)
class JobRow:
    """One kind of job: ``kind`` and ``method`` name it; ``run`` the recipe, or ``missing`` the pyazr call it
    waits for; ``needs`` the package it imports; ``settings`` the defaults a request may override."""
    kind: str
    method: str
    doc: str
    run: object = None
    needs: str = None
    settings: dict = field(default_factory=dict)
    missing: str = None


JOBS = (
    JobRow("fit", "scipy_lbfgsb", "scipy's L-BFGS-B on the objective and its gradient", run=fit_scipy_lbfgsb,
           needs="scipy", settings={"max_iterations": 1000}),
    JobRow("fit", "iminuit_migrad", "iminuit's MIGRAD with the gradient, then HESSE", run=fit_iminuit_migrad,
           needs="iminuit", settings={"max_calls": None, "hesse": True}),
    JobRow("fit", "engine_migrad", "the engine's MIGRAD (the program's reference fit)",
           missing="Session.fit_with_engine(method='migrad') with progress and stop (pyazr3.md decision 8)"),
    JobRow("fit", "engine_lm", "the engine's Levenberg-Marquardt", missing="Session.fit_with_engine(method='lm')"),
    JobRow("minos", "engine", "MINOS intervals of the engine's fit",
           missing="Session.fit_with_engine(minos=True)"),
    JobRow("sample", "emcee", "the emcee sampler on the log-probability", missing="Objective.log_prob (pyazr3.md 2.6)"),
    JobRow("bands", "covariance", "uncertainty bands of a fit",
           missing="pyazr.bands on a Session and an Evaluation (pyazr3.md 2.2)"),
    JobRow("rate", "engine", "the reaction rate", missing="Session.reaction_rate (web_ui.md finding 6)"),
)


def job_row(kind, method):
    for row in JOBS:
        if row.kind == kind and row.method == method:
            return row
    known = ", ".join(f"{r.kind}/{r.method}" for r in JOBS)
    raise NotFound(f"no job {kind}/{method} (the job table: {known})", known=[f"{r.kind}/{r.method}" for r in JOBS])


def job_table():
    """The job table as data (GET /api/jobs/kinds)."""
    out = []
    for r in JOBS:
        available = r.missing is None and _importable(r.needs)
        out.append({"kind": r.kind, "method": r.method, "doc": r.doc, "settings": dict(r.settings),
                    "needs": r.needs, "available": available, "missing": r.missing})
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


def run_job(channel):
    """A job worker: read the job's spec, build the session, run the recipe, report, write the checkpoint."""
    from .engine import open_objective
    spec = channel.receive()
    record = {"id": spec["id"], "kind": spec["kind"], "method": spec["method"], "rev": spec["rev"],
              "settings": spec["settings"], "state": "building", "progress": None, "best": None, "result": None}
    checkpoint = spec["checkpoint"]
    stop = {"asked": False}

    def stop_requested():
        if not stop["asked"] and channel.poll():
            message = channel.receive()
            stop["asked"] = message is None or message.get("op") == "stop"
        return stop["asked"]

    def report(message):
        channel.send(message)
        record["progress"] = {k: message[k] for k in ("calls", "objective", "best", "elapsed")}
        if message["best_x"] is not None:
            record["best"] = {"objective": message["best"], "x": message["best_x"]}
            write_checkpoint(checkpoint, record)

    try:
        write_checkpoint(checkpoint, record)
        t0 = time.perf_counter()
        fg, x0, names = open_objective(spec["session_file"], spec["cwd"])
        channel.send({"event": "built", "build_s": time.perf_counter() - t0, "names": names})
        if spec.get("x0"):
            x0 = vector_from_names(spec["x0"], names, x0)
        record["state"] = "running"
        write_checkpoint(checkpoint, record)
        row = job_row(spec["kind"], spec["method"])
        tracker = Tracker(fg, names, report=report, stop_requested=stop_requested)
        try:
            status, message, extra = row.run(tracker, x0, spec["settings"], names)
        except StopRequested:
            status, message, extra = "stopped", "stopped on request", {}
        best_x = None if tracker.best_x is None else dict(zip(names, tracker.best_x.tolist()))
        record.update(state=status, best={"objective": tracker.best, "x": best_x} if best_x else None,
                      result={"status": status, "message": message, "calls": tracker.calls,
                              "start_objective": tracker.first,
                              "elapsed": time.perf_counter() - tracker.t0, **extra})
    except Exception as e:  # noqa: BLE001  -- reported as the job's failure
        cls, text, _ = describe_exception(e)
        record.update(state="failed", result={"status": "failed", "error": cls, "message": text})
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

    def __init__(self, record, process=None):
        self.record, self.process, self.session_file = record, process, None
        self.lock = threading.Lock()
        self.stop_asked_at = None

    def view(self):
        with self.lock:
            out = dict(self.record)
            out["pid"] = self.process.pid if self.process is not None else None
            return out


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
                self.jobs[record["id"]] = Job(record)

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
        unknown = sorted(set(settings or {}) - set(row.settings))
        if unknown:
            raise Refused(f"job {kind}/{method} has no setting {', '.join(unknown)}", allowed=sorted(row.settings))
        merged = dict(row.settings, **(settings or {}))
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
        job = Job({"id": job_id, "kind": kind, "method": method, "rev": rev, "settings": merged,
                   "state": "building", "progress": None, "best": None, "result": None,
                   "started": time.time()}, process)
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
                        job.record["progress"] = {k: message[k] for k in ("calls", "objective", "best", "elapsed")}
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
