"""Spawn and supervise an AZURE2 ``--use-api`` subprocess.

Each :class:`server` owns exactly one AZURE2 process bound to one port.  Unlike
the original implementation -- which launched the binary through a shell on a
detached thread and never killed it -- this uses :class:`subprocess.Popen` so
the process can be terminated deterministically and probed for liveness.
"""

import os
import shutil
import signal
import subprocess


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


class server:

    def __init__(self, port, file, binary=None, verbose=False, extra_args=None):
        """Start an AZURE2 API server.

        Parameters
        ----------
        port : TCP port the server listens on.
        file : the ``.azr`` configuration file to load.
        binary : path to the AZURE2 executable (auto-detected if ``None``).
        verbose : if ``False`` (default) the subprocess' stdout/stderr are
            discarded; if ``True`` they are inherited from this process.
        extra_args : optional list of additional command-line arguments.
        """
        self.port = port
        self.file = file
        self.binary = binary or _default_binary()
        self.verbose = verbose
        self.extra_args = list(extra_args) if extra_args else []
        self.process = None
        self.start()

    def __del__(self):
        self.stop()

    def start(self):
        if not os.path.isfile(self.binary):
            raise ServerError(
                f"AZURE2 binary not found at {self.binary!r}. Set the "
                "AZURE2_BINARY environment variable or pass binary=..."
            )

        cmd = [self.binary, "--no-gui", "--use-api",
               str(self.port), self.file, *self.extra_args]

        out = None if self.verbose else subprocess.DEVNULL

        # start_new_session lets us signal the whole group on POSIX, ensuring
        # no orphaned children survive.
        kwargs = {}
        if os.name == "posix":
            kwargs["start_new_session"] = True

        try:
            self.process = subprocess.Popen(
                cmd, stdout=out, stderr=out, **kwargs
            )
        except OSError as err:
            raise ServerError(f"Failed to launch AZURE2: {err}") from err

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
        self.process = None
