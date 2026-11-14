"""The engine's text during a library call: kept off standard output (:func:`held`).

The engine prints what it does on the process's standard output (file descriptor 1): progress bars of a reaction rate,
the external-capture integrals it recomputes, the ERROR lines of a project it refuses, a level its transformation
cannot reach.  A library call is quiet, so the calls that enter the engine -- every call of a :class:`pyazr.Session`
(``pyazr.locking.engine_call``) and its set-up, a :class:`pyazr.rate.RateSession` and its rate -- run inside
:func:`held`: descriptor 1 points at a temporary file while the engine runs, and what the engine wrote is
then handed to the logger ``pyazr.engine`` (Python's :mod:`logging`) line by line, the frames of a progress bar left
out: the ERROR and WARNING lines at level WARNING (shown on standard error unless the caller configures logging
otherwise), every other line at INFO (not shown unless asked: ``logging.basicConfig(level=logging.INFO)``).  The
lines are also returned to the caller, so that a refusal can quote them (``Session``'s set-up errors do); the lines
of a call that raised are logged at INFO only, since its error carries them.

Python's own standard output is not held: while descriptor 1 points at the file, ``sys.stdout`` (when it is the
stream of descriptor 1) is a stream on a copy of the real descriptor 1, so that ``print`` -- in a fit's ``progress``
callback, or in another thread -- is shown as before.  Descriptor 1 belongs to the process, not to a thread: what
C code of any thread writes there while a call is inside goes to the file.  Calls from several threads share one
redirection (the first in makes it, the last out restores descriptor 1 and ``sys.stdout``).
"""
import contextlib
import logging
import os
import sys
import tempfile
import threading

LOGGER = logging.getLogger("pyazr.engine")

_lock = threading.Lock()
_depth = 0       # how many calls are inside held() now
_saved = None    # a duplicate of descriptor 1 as it was before the first of them
_sink = None     # the temporary file descriptor 1 points at; the last call out logs all of it
_python = None   # (sys.stdout as it was, the stream on the copy of descriptor 1 that replaces it), or None


try:  # the C library, loaded once (a session closed while the interpreter shuts down can no longer import)
    import ctypes
    _LIBC = ctypes.CDLL(None)
except (ImportError, OSError):  # pragma: no cover
    _LIBC = None


def _flush_c_stdout():
    """Flush the C library's buffer of standard output (the engine's C++ stream writes through it)."""
    try:
        if _LIBC is not None:
            _LIBC.fflush(None)
    except (OSError, AttributeError, ValueError):
        pass


def _flush_python_stdout():
    try:
        sys.stdout.flush()
    except (AttributeError, ValueError, OSError):
        pass


def _python_stdout_to(saved):
    """Point ``sys.stdout`` at a stream on a copy of ``saved`` (the real descriptor 1) when it is the stream of
    descriptor 1; returns (the old stream, the new one), or None when it is another stream (a notebook's, a test
    runner's), which the redirection of descriptor 1 does not touch."""
    old = sys.stdout
    try:
        if old is None or old.fileno() != 1:
            return None
    except (AttributeError, ValueError, OSError):
        return None
    new = os.fdopen(os.dup(saved), "w", buffering=1, encoding=getattr(old, "encoding", None) or "utf-8",
                    errors=getattr(old, "errors", None) or "strict")
    sys.stdout = new
    return old, new


def _python_stdout_back(old, new):
    try:
        new.flush()
    finally:
        if sys.stdout is new:
            sys.stdout = old
        new.close()


def _engine_lines(text):
    """The lines of the engine's text, without the frames of its progress bars (a carriage return starts a new
    frame of the same line: only the last frame is kept) and without empty lines."""
    out = []
    for line in text.splitlines():
        line = line.rsplit("\r", 1)[-1].rstrip()
        if line.strip():
            out.append(line)
    return out


def _log(lines, refused):
    """Log the engine's lines; when the call ``refused`` (raised), every line at INFO: the error quotes them."""
    warn = False
    for line in lines:
        head = line.lstrip("* \t")
        if head.startswith(("ERROR", "WARNING", "Warning", "Error")):
            warn = True
        elif not line[:1].isspace():
            warn = False  # an indented line continues the message above it
        LOGGER.log(logging.WARNING if warn and not refused else logging.INFO, "%s", line)


@contextlib.contextmanager
def held():
    """Run the block with descriptor 1 pointed at a temporary file; yields a list that holds, after the block, the
    engine's lines written while this call was inside (with what calls of other threads wrote in that time).  Where
    descriptor 1 is not open (a daemon), the block runs as it is and the list stays empty."""
    global _depth, _saved, _sink, _python
    mine = []
    if sys.is_finalizing():  # a session closed while the interpreter shuts down: nothing to hold
        yield mine
        return
    with _lock:
        if _depth == 0:
            _flush_python_stdout()
            _flush_c_stdout()
            try:
                saved = os.dup(1)
            except OSError:
                saved = None
            if saved is not None:
                _sink = tempfile.TemporaryFile()
                os.dup2(_sink.fileno(), 1)
                _python = _python_stdout_to(saved)
            _saved = saved
        else:
            _flush_c_stdout()
        start = os.lseek(_sink.fileno(), 0, os.SEEK_CUR) if _sink is not None else 0
        _depth += 1
    refused = False
    try:
        yield mine
    except BaseException:
        refused = True
        raise
    finally:
        everything = None
        with _lock:
            _depth -= 1
            _flush_python_stdout()
            _flush_c_stdout()
            if _sink is not None:
                fd = _sink.fileno()
                end = os.lseek(fd, 0, os.SEEK_CUR)
                mine.extend(_engine_lines(os.pread(fd, end - start, start).decode("utf-8", "replace")))
                if _depth == 0:
                    everything = _engine_lines(os.pread(fd, end, 0).decode("utf-8", "replace"))
            if _depth == 0:
                if _python is not None:
                    _python_stdout_back(*_python)
                if _saved is not None:
                    os.dup2(_saved, 1)
                    os.close(_saved)
                    _sink.close()
                _saved = _sink = _python = None
        if everything:
            _log(everything, refused)
