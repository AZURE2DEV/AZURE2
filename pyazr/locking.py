"""The lock of a session: every call that enters the engine of a :class:`pyazr.Session` holds it.

The engine evaluates on the session's one compound nucleus and data, so two calls at once on one session would race
(two threads calling one ``Objective`` crash the process without it).  Each session owns one ``threading.RLock``
(:func:`session_lock`); a method of the session that enters the engine is wrapped by :func:`engine_call`, and the
objects built on a session (:class:`pyazr.Objective`, :func:`pyazr.evaluate.evaluate`) take the same lock.  It is
re-entrant, so a locked method may call another.  Two sessions have two locks and run side by side.
"""
import functools
import threading


def session_lock(session):
    """The lock of ``session`` (made on first use; a session makes it when it opens)."""
    lock = session.__dict__.get("_engine_lock")
    if lock is None:
        lock = session.__dict__.setdefault("_engine_lock", threading.RLock())
    return lock


def engine_call(method):
    """Wrap a method of the session so that it holds the session's lock while it runs (``method.takes_lock`` is
    True on the wrapper: ``tests/pyazr/session_lock_test.py`` checks that every public member that enters the engine
    has it)."""
    @functools.wraps(method)
    def locked(self, *args, **kwargs):
        with session_lock(self):
            return method(self, *args, **kwargs)
    locked.takes_lock = True
    return locked
