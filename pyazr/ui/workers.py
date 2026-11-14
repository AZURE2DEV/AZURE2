"""Worker processes of the web UI (docs/upgrade/web_ui.md sections 3 and 3.4) and the compute budget (3.2).

The server never builds an engine session.  A **session worker** (one per open project) holds the session of the
current project text and answers requests one at a time; a **job worker** (one per job, :mod:`pyazr.ui.jobs`) builds
its own session and runs a recipe.  Both are ``python -m pyazr.ui.workers session|job``: requests arrive as JSON
lines on standard input, answers and progress leave as JSON lines on a pipe of their own (``--out-fd``), and the
engine's notes, which it prints on standard output, go to the worker's log file instead.

**Death with the server.**  A worker sets ``PR_SET_PDEATHSIG`` (Linux; SIGKILL when the thread that started it ends:
every worker is started from one thread that lives as long as the server, :class:`Spawner`) and runs a watchdog
thread that exits when its parent changes (any POSIX system; it needs the GIL, so it acts between engine calls).
"""
import json
import os
import signal
import subprocess
import sys
import threading
import time
import warnings

from .errors import Busy, WorkerCrashed, WorkerError, describe_exception


# -- the compute budget -------------------------------------------------------------------------------------------

class Budget:
    """The thread budget of the server: ``threads`` in all (``main(threads=...)``, default 1).

    With one thread, one process computes at a time: a session request while a job runs, or a second job, is refused
    (:class:`pyazr.ui.errors.Busy`, with the job's id and its stop route).  With N >= 2 threads the session worker
    keeps one for interactive work and the jobs share the other N - 1, one thread each.  Each worker is started with
    ``OMP_NUM_THREADS`` set to its share."""

    def __init__(self, threads):
        if int(threads) < 1:
            raise ValueError(f"the thread budget is at least 1, not {threads}")
        self.threads = int(threads)
        self.jobs = {}  # job id -> threads held
        self.lock = threading.Lock()

    def job_slots(self):
        return 1 if self.threads == 1 else self.threads - 1

    def require_session(self):
        """Refuse a session request while jobs hold every thread."""
        with self.lock:
            if self.threads == 1 and self.jobs:
                self._refuse("an evaluation")

    def take_job(self, job_id):
        """Reserve one thread for job ``job_id``, or refuse; returns the threads given."""
        with self.lock:
            if len(self.jobs) >= self.job_slots():
                self._refuse("another job")
            self.jobs[job_id] = 1
            return 1

    def give_back(self, job_id):
        with self.lock:
            self.jobs.pop(job_id, None)

    def _refuse(self, what):
        job = next(iter(self.jobs))
        raise Busy(f"the compute budget ({self.threads} thread(s)) is held by job {job}: {what} waits until it ends",
                   allowed=f"{self.threads} thread(s) in all", fix=f"POST /api/jobs/{job}/stop, or wait",
                   job=job, stop=f"/api/jobs/{job}/stop")

    def describe(self):
        with self.lock:
            return {"threads": self.threads, "job_slots": self.job_slots(), "jobs": dict(self.jobs)}


# -- server side: starting and talking to workers ---------------------------------------------------------------

class Spawner:
    """Starts every worker from one long-lived thread: ``PR_SET_PDEATHSIG`` fires when the *thread* that started a
    process ends, and the request handlers run on pool threads that come and go."""

    def __init__(self):
        self._lock = threading.Lock()
        self._requests = []
        self._wake = threading.Condition(self._lock)
        self._thread = threading.Thread(target=self._loop, name="pyazr-ui-spawner", daemon=True)
        self._thread.start()

    def _loop(self):
        while True:
            with self._wake:
                while not self._requests:
                    self._wake.wait()
                args, kwargs, box, done = self._requests.pop(0)
            try:
                box.append(subprocess.Popen(*args, **kwargs))
            except Exception as e:  # noqa: BLE001  (handed to the caller)
                box.append(e)
            done.set()

    def popen(self, *args, **kwargs):
        box, done = [], threading.Event()
        with self._wake:
            self._requests.append((args, kwargs, box, done))
            self._wake.notify()
        done.wait()
        if isinstance(box[0], Exception):
            raise box[0]
        return box[0]


class WorkerProcess:
    """A worker process seen from the server: ``send`` a JSON message on its standard input, ``receive`` the next
    JSON line of its output pipe; ``kill`` it; ``pid``."""

    def __init__(self, spawner, role, *, threads, log_path, cwd):
        read_fd, write_fd = os.pipe()
        env = dict(os.environ, OMP_NUM_THREADS=str(threads), PYTHONUNBUFFERED="1")
        env["PYTHONPATH"] = os.pathsep.join([_package_root()] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep)
                                                                 if p])
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        self.log_path = log_path
        with open(log_path, "ab") as log:
            self.proc = spawner.popen(
                [sys.executable, "-m", "pyazr.ui.workers", role, "--out-fd", str(write_fd), "--ppid", str(os.getpid())],
                stdin=subprocess.PIPE, stdout=log, stderr=log, pass_fds=(write_fd,), env=env, cwd=cwd)
        os.close(write_fd)
        self.out = os.fdopen(read_fd, "rb")
        self.pid = self.proc.pid

    def send(self, message):
        try:
            self.proc.stdin.write((json.dumps(message) + "\n").encode())
            self.proc.stdin.flush()
        except (BrokenPipeError, ValueError, OSError):
            raise self._crashed() from None

    def receive(self):
        """The next message, or WorkerCrashed when the worker ended without one."""
        line = self.out.readline()
        if not line:
            raise self._crashed()
        return json.loads(line)

    def _crashed(self):
        code = self.proc.poll()
        if code is None:
            try:
                code = self.proc.wait(5)
            except subprocess.TimeoutExpired:
                code = None
        return WorkerCrashed(f"worker {self.pid} ended without an answer (exit status {code}); its log: "
                             f"{self.log_path}", returncode=code, log=self.log_path)

    def alive(self):
        return self.proc.poll() is None

    def kill(self):
        if self.alive():
            self.proc.kill()
        self.proc.wait()
        for f in (self.proc.stdin, self.out):
            try:
                f.close()
            except OSError:
                pass

    def stop(self, timeout=5.0):
        """Close its input (a session worker then returns) and wait; kill it when it does not end."""
        try:
            self.proc.stdin.close()
        except OSError:
            pass
        try:
            self.proc.wait(timeout)
        except subprocess.TimeoutExpired:
            pass
        self.kill()


class SessionWorker:
    """The session worker of the open project: one request at a time (a lock), started on the first call."""

    def __init__(self, spawner, *, log_path, cwd):
        self.spawner, self.log_path, self.cwd = spawner, log_path, cwd
        self.process = None
        self.lock = threading.Lock()
        self.built_from = None  # the session file text the worker's session was built from
        self.warnings = []  # the warnings of the last call (the engine's notices among them)

    def call(self, op, **args):
        """``op`` of :class:`pyazr.ui.engine.Engine` with ``args``; its answer, or its exception as WorkerError.
        The caller holds :attr:`lock`."""
        if self.process is None or not self.process.alive():
            self.process = WorkerProcess(self.spawner, "session", threads=1, log_path=self.log_path, cwd=self.cwd)
            self.built_from = None
        self.process.send({"op": op, "args": args})
        try:
            answer = self.process.receive()
        except WorkerCrashed:
            self.process, self.built_from = None, None
            raise
        self.warnings = answer.get("warnings", [])
        if not answer["ok"]:
            raise WorkerError(answer["classes"], answer["message"], answer.get("diagnostics"))
        return answer["result"]

    def pid(self):
        return self.process.pid if self.process is not None and self.process.alive() else None

    def close(self):
        with self.lock:
            if self.process is not None:
                self.process.stop()
                self.process = None
                self.built_from = None


def _package_root():
    """The directory that holds the package ``pyazr`` (the worker imports the same pyazr as the server)."""
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# -- worker side ----------------------------------------------------------------------------------------------------

def die_with_parent(ppid):
    """Make this process end when the server ``ppid`` ends: PR_SET_PDEATHSIG on Linux, a watchdog thread elsewhere
    (and on Linux too, for a parent that died before the call)."""
    if sys.platform.startswith("linux"):
        try:
            import ctypes
            libc = ctypes.CDLL(None, use_errno=True)
            PR_SET_PDEATHSIG = 1
            libc.prctl(PR_SET_PDEATHSIG, signal.SIGKILL, 0, 0, 0)
        except (OSError, AttributeError):
            pass
    if os.getppid() != ppid:
        os._exit(3)

    def watch():
        while True:
            time.sleep(0.2)
            if os.getppid() != ppid:
                os._exit(3)

    threading.Thread(target=watch, name="pyazr-ui-watchdog", daemon=True).start()


class Channel:
    """The worker's end: JSON lines in on standard input, out on ``out_fd``."""

    def __init__(self, out_fd):
        self.out = os.fdopen(out_fd, "wb")
        self.inp = sys.stdin.buffer

    def send(self, message):
        self.out.write((json.dumps(message) + "\n").encode())
        self.out.flush()

    def receive(self):
        line = self.inp.readline()
        return json.loads(line) if line else None


def serve_session(channel):
    """The session worker's loop: one request, one answer, until standard input closes."""
    from .engine import Engine
    engine = Engine()
    while True:
        message = channel.receive()
        if message is None:
            return
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                result = getattr(engine, message["op"])(**message.get("args", {}))
                answer = {"ok": True, "result": result}
            except Exception as e:  # noqa: BLE001  -- the server maps the classes to an HTTP status
                classes, text, diagnostics = describe_exception(e)
                answer = {"ok": False, "classes": classes, "message": text, "diagnostics": diagnostics}
        # the engine's notices (pyazr.EngineWarning) and pyazr's other warnings of this request, for the client
        answer["warnings"] = [{"category": w.category.__name__, "message": str(w.message)} for w in caught]
        channel.send(answer)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(prog="python -m pyazr.ui.workers", description="a worker of the web UI's server")
    ap.add_argument("role", choices=["session", "job"])
    ap.add_argument("--out-fd", type=int, required=True)
    ap.add_argument("--ppid", type=int, required=True)
    a = ap.parse_args(argv)
    die_with_parent(a.ppid)
    channel = Channel(a.out_fd)
    if a.role == "session":
        serve_session(channel)
    else:
        from .jobs import run_job
        run_job(channel)


if __name__ == "__main__":
    main()
