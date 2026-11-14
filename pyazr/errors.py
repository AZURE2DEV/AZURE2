"""The errors and warnings of pyazr (docs/upgrade/pyazr3.md section 2.5).

Every error pyazr raises on purpose is a :class:`Error`.  Each class also derives from the builtin exception a user
would catch at that place, so that ``except ValueError`` (or ``KeyError``, ``RuntimeError``, ...) written against an
older pyazr keeps catching it:

- :class:`FormatError` (a ``ValueError``): the engine refuses a project text or a record edit.
- :class:`ParameterError` (``ValueError``): a parameter vector, name, index or value that the call cannot take.
  Where a lookup raised a bare builtin before, the class derives from it too: :class:`ParameterKeyError` (also a
  ``KeyError``: an unknown parameter, level, pair, data set or THM experiment), :class:`ParameterIndexError` (also an
  ``IndexError``: an unknown pair number), :class:`ParameterVectorError` (also an :class:`EngineError`: a
  free-parameter vector of the wrong length or with a NaN, refused by the engine).
- :class:`EngineError` (``RuntimeError``): the engine refuses a set-up or an evaluation; the same class as
  ``AZURE2Error`` (``pyazr.AZURE2Error``, ``pyazr._azure2.AZURE2Error``).
- :class:`DomainError` (``ValueError``): an argument of a physics function outside its domain
  (:class:`DomainZeroDivision`, also a ``ZeroDivisionError``, where the formula divided by zero before).
- :class:`SessionClosed` (``RuntimeError``): a call on a closed session.
- :class:`MissingFile` (``FileNotFoundError``): a project, working directory or parameter file that does not exist.

A message says (1) what was wrong with the value given, (2) what is allowed, (3) which call shows or fixes it, and
names objects by their label.  Where the message is built from these three parts, they are attributes too:
``e.what``, ``e.allowed``, ``e.fix`` (texts; ``None`` where a message is not built from parts, as the engine's own
messages); a :class:`FormatError` takes them from its first diagnostic of severity error (the validator's ``message``,
``allowed`` and ``fix``).  An engine notice that a session meets (the engine printed it on standard output before)
is a :class:`EngineWarning` with the file concerned and the consequence (``warnings.warn``; turn it into an error with
``warnings.simplefilter("error", pyazr.EngineWarning)``).

This module imports nothing of the engine: the compiled module ``_azure2`` takes :class:`EngineError` and
:class:`ParameterVectorError` from here when it is imported (``api/src/pybind_api.cpp``), so that the class a C++
refusal raises and the class a Python check raises are the same.
"""

import difflib

__all__ = ["Error", "FormatError", "ParameterError", "ParameterKeyError", "ParameterIndexError",
           "ParameterVectorError", "EngineError", "AZURE2Error", "DomainError", "DomainZeroDivision", "SessionClosed",
           "MissingFile", "EngineWarning"]


class Error(Exception):
    """The base of every error pyazr raises on purpose (``except pyazr.Error`` catches them all).

    ``what`` (what was wrong with the value given), ``allowed`` (what is allowed) and ``fix`` (which call shows or
    fixes it): the three parts of the message where it is built from them, else ``None``.  They are given as keyword
    arguments, ``Error(message, what=..., allowed=..., fix=...)``; the message is the text as written."""

    what = None
    allowed = None
    fix = None

    def __init__(self, *args, what=None, allowed=None, fix=None):
        super().__init__(*args)
        self.what, self.allowed, self.fix = what, allowed, fix


class FormatError(Error, ValueError):
    """A project or an edit the engine refuses.

    ``diagnostics`` is the engine's list (``AZURE2 --validate --json``, docs/upgrade/project_format3.md section 6.2):
    dicts with ``severity``, ``code``, ``line``, ``record``, ``field``, ``found``, ``allowed``, ``message`` and a
    syntax-only ``fix``; ``codes`` the codes in order.  The message is the engine's text, one line per diagnostic.
    """

    def __init__(self, message, diagnostics=(), **parts):
        self.diagnostics = list(diagnostics)
        lines = [message] + ["  " + _diagnostic_text(d) for d in self.diagnostics]
        first = next((d for d in self.diagnostics if d.get("severity", "error") == "error"), None)
        if first is not None and not parts:
            parts = {"what": first.get("message") or None, "allowed": first.get("allowed") or None,
                     "fix": first.get("fix") or None}
        super().__init__("\n".join(lines), **parts)

    @property
    def codes(self):
        return [d.get("code") for d in self.diagnostics]


def _diagnostic_text(d):
    where = ".".join(x for x in (d.get("record") or "", d.get("field") or "") if x)
    return f"{d.get('code')}: {where + ': ' if where else ''}{d.get('message')}"


class ParameterError(Error, ValueError):
    """A parameter vector, a parameter name or index, or a parameter value that the call cannot take."""


class ParameterKeyError(ParameterError, KeyError):
    """An unknown name or key: a parameter, level, pair, data set or THM experiment that the session does not have.

    A ``KeyError`` too (these lookups raised a bare ``KeyError`` before); unlike ``KeyError`` its text is the message,
    not the message's repr."""

    def __str__(self):
        return str(self.args[0]) if self.args else ""


class ParameterIndexError(ParameterError, IndexError):
    """An index out of range (a pair number that the session does not have); an ``IndexError`` too."""


class AZURE2Error(Error, RuntimeError):
    """The engine refuses a set-up or an evaluation; the message is the engine's.

    Its name in pyazr 3.0 is :data:`EngineError` (one class, two names): ``pyazr.EngineError``,
    ``pyazr.AZURE2Error`` and ``pyazr._azure2.AZURE2Error`` are this class.  The class keeps the name ``AZURE2Error``
    (``type(e).__name__``) that scripts and tests compare against."""


#: The name of :class:`AZURE2Error` in pyazr 3.0 (one class, two names).
EngineError = AZURE2Error


class ParameterVectorError(ParameterError, AZURE2Error):
    """A free-parameter vector of the wrong length, or with a NaN or an infinity, refused by the engine before it
    computes anything.  A ``ValueError`` (as a :class:`ParameterError`) and an :class:`EngineError` (``AZURE2Error``,
    which the engine raised for a wrong length before)."""


class DomainError(Error, ValueError):
    """An argument of a physics function outside its domain: the message names the function, the element's index and
    value, and the domain."""


class DomainZeroDivision(DomainError, ZeroDivisionError):
    """A :class:`DomainError` where the formula divided by zero before (L = -1, an energy of exactly 0 in
    ``pyazr.penetrability``, ``shift``, ``shift_derivative``, ``whittaker``); a ``ZeroDivisionError`` too."""


class SessionClosed(Error, RuntimeError):
    """A call on a session that has been closed (``close()``, or the end of its ``with`` block); a new session is
    needed."""


class MissingFile(Error, FileNotFoundError):
    """A file or directory that a call needs does not exist; a ``FileNotFoundError`` too, with ``filename`` set."""

    def __init__(self, path, message, **parts):
        super().__init__(2, message, str(path), **parts)  # errno ENOENT, as FileNotFoundError

    def __str__(self):
        return self.strerror


class EngineWarning(UserWarning):
    """A notice of the engine that a session met (it was printed on standard output before): the file concerned and
    the consequence."""


def close_names(name, names, n=3):
    """Up to ``n`` entries of ``names`` close to ``name`` (difflib), as a text for a message: "" when none is close."""
    close = difflib.get_close_matches(str(name), [str(x) for x in names], n=n)
    return f" (did you mean {', '.join(repr(c) for c in close)}?)" if close else ""


def listed(names, limit=12):
    """``names`` as a short text for a message: every one up to ``limit``, then how many more."""
    names = [str(x) for x in names]
    text = ", ".join(names[:limit])
    return text + (f", ... ({len(names) - limit} more)" if len(names) > limit else "") if names else "none"
