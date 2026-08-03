"""Spawn and supervise an AZURE2 ``--use-api`` subprocess.

Each :class:`server` owns exactly one AZURE2 process bound to one port.  Unlike
the original implementation -- which launched the binary through a shell on a
detached thread and never killed it -- this uses :class:`subprocess.Popen` so
the process can be terminated deterministically and probed for liveness.
"""

import atexit
import os
import shutil
import signal
import subprocess
import sys
import threading
import weakref

# Printed by AZURESocket.cpp once the server is bound and listening; the trailing
# token is the actual (possibly OS-assigned) port.
_LISTENING_MARKER = "AZURE2_API_LISTENING"


def _default_binary():
    """Locate the AZURE2 executable.

    Resolution order:
      1. ``$AZURE2_BINARY`` environment variable.
      2. ``build/src/AZURE2`` relative to the repository root (two levels up
         from this file).
      3. ``AZURE2`` on ``$PATH``.
    """
    env = os.environ.get("AZURE2_BINARY")
    if env:
        return env

    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(here)
    candidate = os.path.join(repo_root, "build", "src", "AZURE2")
    if os.path.isfile(candidate):
        return candidate

    found = shutil.which("AZURE2")
    if found:
        return found

    # Fall back to the conventional location; the error surfaces clearly on
    # spawn if it is missing.
    return candidate


class ServerError(RuntimeError):
    """Raised when the AZURE2 subprocess cannot be started."""


# Every live server, so interpreter exit can stop them all. __del__ alone is not
# enough: it may not run at shutdown, and it certainly does not run when a
# notebook kernel is restarted or the interpreter is killed outright, which is
# how stray AZURE2 processes accumulate. A WeakSet keeps this from holding
# servers alive past their natural lifetime.
_LIVE_SERVERS = weakref.WeakSet()


def _stop_all_servers():
    """Stop every server this interpreter started. Registered with atexit."""
    for srv in list(_LIVE_SERVERS):
        try:
            srv.stop()
        except Exception:      # never let cleanup raise during shutdown
            pass


atexit.register(_stop_all_servers)


class server:

    def __init__(self, port, file, binary=None, verbose=False, extra_args=None,
                 cwd=None):
        """Start an AZURE2 API server.

        Parameters
        ----------
        port : TCP port to request.  Pass ``0`` to let the OS assign a
            guaranteed-unique free port; the actual port is read back from the
            server and exposed as :attr:`port` (see :meth:`wait_until_listening`).
        file : the ``.azr`` configuration file to load.
        binary : path to the AZURE2 executable (auto-detected if ``None``).
        verbose : if ``False`` (default) the subprocess' output is discarded; if
            ``True`` it is echoed to this process' stdout.
        extra_args : optional list of additional command-line arguments.
        cwd : working directory for the subprocess.  A ``.azr`` file names its
            output and checks directories *relative to the process' cwd* (e.g.
            ``output/``), so a caller running a model that lives elsewhere must
            set this to the model's directory or AZURE2 will scatter its results
            into the caller's cwd.  Defaults to the ``.azr`` file's own
            directory, which is what that relative convention means.
        """
        self.port = port                      # updated to the real port on bind
        self.file = file
        self.binary = binary or _default_binary()
        self.verbose = verbose
        self.extra_args = list(extra_args) if extra_args else []
        self.cwd = cwd if cwd is not None else (
            os.path.dirname(os.path.abspath(file)) or None)

        # Popen resolves a *relative* executable path against the child's cwd,
        # not ours -- so once cwd is set, a caller's relative binary (the
        # documented ``binary='../../build/src/AZURE2'`` form) would silently
        # resolve against the model directory and fail to launch.  Anchor it to
        # the caller's cwd now, while that is still what it means.  A bare name
        # with no separator is left alone so $PATH lookup keeps working.
        if os.sep in self.binary and not os.path.isabs(self.binary):
            self.binary = os.path.abspath(self.binary)
        self.process = None
        self._listening = threading.Event()   # set once the port is known
        self._reader = None
        self.start()
        _LIVE_SERVERS.add(self)

    def __del__(self):
        self.stop()

    def start(self):
        if not os.path.isfile(self.binary):
            raise ServerError(
                f"AZURE2 binary not found at {self.binary!r}. Set the "
                "AZURE2_BINARY environment variable or pass binary=..."
            )

        # The file is passed absolute: we hand the subprocess a cwd of its own
        # (the model's directory), so a relative path would resolve against the
        # wrong root.
        cmd = [self.binary, "--no-gui", "--use-api",
               str(self.port), os.path.abspath(self.file), *self.extra_args]

        # start_new_session lets us signal the whole group on POSIX, ensuring
        # no orphaned children survive.
        kwargs = {}
        if os.name == "posix":
            kwargs["start_new_session"] = True
        if self.cwd is not None:
            kwargs["cwd"] = self.cwd

        try:
            # Capture stdout (with stderr merged in) so we can read the
            # "AZURE2_API_LISTENING <port>" line; a reader thread drains the pipe
            # continuously so a chatty subprocess can never fill it and deadlock.
            self.process = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, **kwargs
            )
        except OSError as err:
            raise ServerError(f"Failed to launch AZURE2: {err}") from err

        self._reader = threading.Thread(target=self._drain_output, daemon=True)
        self._reader.start()

    def _drain_output(self):
        """Read the subprocess' output, capturing the listening port."""
        try:
            for line in self.process.stdout:
                if not self._listening.is_set() and _LISTENING_MARKER in line:
                    try:
                        self.port = int(line.split()[-1])
                    except (ValueError, IndexError):
                        pass
                    self._listening.set()
                if self.verbose:
                    sys.stdout.write(line)
        finally:
            # On EOF (process exited) unblock any waiter so it can detect death.
            self._listening.set()

    def wait_until_listening(self, timeout=60.0):
        """Block until the server reports its bound port; return that port."""
        if not self._listening.wait(timeout):
            raise ServerError(
                f"AZURE2 did not report a listening port within {timeout}s.")
        if self.process.poll() is not None:
            raise ServerError(
                f"AZURE2 exited (code {self.process.returncode}) before "
                "listening; check the .azr file and binary.")
        return self.port

    def is_alive(self):
        return self.process is not None and self.process.poll() is None

    def stop(self, timeout=10.0):
        """Terminate the subprocess, escalating to SIGKILL if needed."""
        if self.process is None:
            return
        if self.process.poll() is None:
            try:
                if os.name == "posix":
                    os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                else:
                    self.process.terminate()
            except (ProcessLookupError, OSError):
                pass
            try:
                self.process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                try:
                    if os.name == "posix":
                        os.killpg(os.getpgid(self.process.pid), signal.SIGKILL)
                    else:
                        self.process.kill()
                except (ProcessLookupError, OSError):
                    pass
                try:
                    self.process.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    pass
        # Closing stdout ends the reader thread (it's a daemon, but tidy up).
        if self.process.stdout is not None:
            try:
                self.process.stdout.close()
            except OSError:
                pass
        # Join the reader so no background thread is left alive -- important if
        # the caller later forks (a live thread at fork can deadlock the child).
        if self._reader is not None:
            self._reader.join(timeout=2.0)
            self._reader = None
        self.process = None
