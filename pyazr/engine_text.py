"""The engine's text during a library call: kept off standard output (:func:`held`).

The engine prints what it does on the process's standard output (file descriptor 1): progress bars of a reaction rate,
the external-capture integrals it recomputes, the ERROR lines of a project it refuses.  A library call is quiet, so
the calls that print -- the set-up of a :class:`pyazr.Session`, a :class:`pyazr.rate.RateSession` and its rate --
run inside :func:`held`: descriptor 1 points at a temporary file while the engine runs, and what the engine wrote is
then handed to the logger ``pyazr.engine`` (Python's :mod:`logging`) line by line, the frames of a progress bar left
out: the ERROR and WARNING lines at level WARNING (shown on standard error unless the caller configures logging
otherwise), every other line at INFO (not shown unless asked: ``logging.basicConfig(level=logging.INFO)``).  The
lines are also returned to the caller, so that a refusal can quote them (``Session``'s set-up errors do); the lines
of a call that raised are logged at INFO only, since its error carries them.

Descriptor 1 belongs to the process, not to a thread: while any thread is inside :func:`held`, whatever any thread
writes to it (the engine, or Python's ``print`` from another thread) goes to the file and then to the logger.  Calls
from several threads share one redirection (the first in makes it, the last out restores descriptor 1).
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


def _flush_c_stdout():
    """Flush the C library's buffer of standard output (the engine's C++ stream writes through it)."""
    try:
        import ctypes
        ctypes.CDLL(None).fflush(None)
    except (OSError, AttributeError):
        pass


def _flush_python_stdout():
    try:
        sys.stdout.flush()
    except (AttributeError, ValueError, OSError):
        pass


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
    global _depth, _saved, _sink
    mine = []
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
                if _saved is not None:
                    os.dup2(_saved, 1)
                    os.close(_saved)
                    _sink.close()
                _saved = _sink = None
        if everything:
            _log(everything, refused)
