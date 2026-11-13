"""Errors of the web UI's API: pyazr's exceptions and the server's own refusals as HTTP statuses (docs/upgrade/web_ui.md
section 4.3).

Every refusal is a JSON body with ``error`` (the class name), ``message`` (the whole text) and the three parts of
pyazr 3.0's messages as fields: ``what`` (what was wrong), ``allowed`` (what is allowed), ``fix`` (which call shows or
fixes it); a :class:`pyazr.FormatError` adds the engine's ``diagnostics``.  Today's pyazr exceptions carry one text,
not the three parts (web_ui.md finding 9): ``what`` is that text and ``allowed`` / ``fix`` stay empty unless the
diagnostics give them.  :func:`status_of` is the one place that maps a class to a status, by class name, so that an
exception raised in a worker process (which reaches the server as its class name and message) maps as one raised in
the server.
"""

# pyazr's classes (today's and those of pyazr3.md section 2.5) and Python's builtins -> HTTP status.  A class not
# listed falls back along its bases (FormatError is a ValueError, so a subclass of it maps as it does).
STATUS = {
    "FormatError": 422, "ParameterError": 422, "DomainError": 422, "ValueError": 422, "TypeError": 422,
    "EngineError": 422, "AZURE2Error": 422,
    "SessionClosed": 409,
    "KeyError": 404, "IndexError": 404, "FileNotFoundError": 404,
    "WorkerCrashed": 500, "RuntimeError": 500,
}


class UIError(Exception):
    """A refusal of the server itself (not of pyazr): ``status``, the message and extra fields of the body."""
    status = 500

    def __init__(self, message, *, what=None, allowed=None, fix=None, **fields):
        super().__init__(message)
        self.what, self.allowed, self.fix, self.fields = what or message, allowed, fix, fields


class Refused(UIError):
    """A request the server refuses as it stands (a fixed parameter given as a value of x, a setting a job lacks)."""
    status = 422


class StaleRevision(UIError):
    """An edit made on a revision that is no longer the current one (another client edited first)."""
    status = 412


class RevisionRequired(UIError):
    """An edit without the revision it was made on."""
    status = 428


class Busy(UIError):
    """The compute budget is full: a job holds it (``job`` and its ``stop`` route in the body)."""
    status = 409


class Conflict(UIError):
    """A request that cannot be served in the server's present state (nothing to undo, the project changed)."""
    status = 409


class NotFound(UIError):
    """An unknown record, level, job or data set; ``known`` lists the names that exist."""
    status = 404


class OutsideRoot(UIError):
    """A path that resolves outside the root directory the server may read and write."""
    status = 403


class Missing(UIError):
    """A route whose pyazr call does not exist yet (``missing``: the call), or one planned for a later step."""
    status = 501


class WorkerCrashed(UIError):
    """A worker process ended without answering (``returncode``)."""
    status = 500


class WorkerError(Exception):
    """An exception raised inside a worker process, re-raised in the server: the class name, message and
    diagnostics travel, the class does not."""

    def __init__(self, cls, message, diagnostics=None):
        super().__init__(message)
        self.cls, self.diagnostics = cls, diagnostics


def status_of(names):
    """The HTTP status of an exception whose class and bases are ``names`` (most derived first)."""
    for name in names:
        if name in STATUS:
            return STATUS[name]
    return 500


def describe_exception(e):
    """(class name, message, diagnostics or None) of an exception, for the pipe of a worker."""
    return type(e).__name__, str(e), _diagnostics(e)


def _diagnostics(e):
    d = getattr(e, "diagnostics", None)
    return [dict(x) for x in d] if d else None


def to_http(e):
    """(status, body) of any exception raised while serving a request."""
    if isinstance(e, UIError):
        body = {"error": type(e).__name__, "message": str(e), "what": e.what, "allowed": e.allowed, "fix": e.fix}
        body.update(e.fields)
        return e.status, body
    if isinstance(e, WorkerError):
        cls, message, diagnostics = e.cls, str(e), e.diagnostics
        status = status_of([cls])
    else:
        cls, message, diagnostics = describe_exception(e)
        status = status_of([c.__name__ for c in type(e).__mro__])
    if isinstance(e, KeyError) and e.args:
        message = str(e.args[0])  # KeyError's str() quotes its argument
    body = {"error": cls, "message": message, "what": message, "allowed": None, "fix": None}
    if diagnostics:
        body["diagnostics"] = diagnostics
        firsts = [d for d in diagnostics if d.get("severity") == "error"] or diagnostics
        body["allowed"] = firsts[0].get("allowed")
        body["fix"] = firsts[0].get("fix")
    return status, body
