"""The argument checks and engine notices of a calculation session (:class:`pyazr.Session`).

Each check raises one of pyazr's errors (:mod:`pyazr.errors`) with a message in three parts: what was wrong with
the value given, what is allowed, which call shows or fixes it.  :func:`warn_engine_notices` turns the notices the
binding held back (``_azure2.Session.take_notices``, ``api/src/pybind_notices.cpp``) into warnings.
"""
import os
import warnings

import numpy as np

from .errors import (AZURE2Error, DomainError, EngineWarning, MissingFile, ParameterIndexError, ParameterVectorError,
                     SessionClosed, listed)


def session_closed(project_file, caller):
    """The error of a call on a closed session; ``caller`` is the name of the session's method that was called."""
    call = f"Session.{caller}" if caller not in ("<module>", "<lambda>") else "this call"
    return SessionClosed(
        f"{call}: this session of {os.path.basename(project_file)} is closed (close() or the end of its with "
        f"block released the engine); a closed session cannot be reopened.  Open a new one: "
        f"pyazr.Session({project_file!r}); is_alive() tells whether a session is open.",
        what=f"this session of {os.path.basename(project_file)} is closed", allowed="a closed session cannot be reopened",
        fix=f"Open a new one: pyazr.Session({project_file!r}); is_alive() tells whether a session is open")


def missing_directory(cwd, project_file):
    """The error of a working directory (``cwd=``) that does not exist."""
    return MissingFile(cwd, f"working directory {cwd!r} (cwd=) does not exist; the project's relative paths are read "
                            f"and written there.  Give an existing directory, or leave cwd=None for the project "
                            f"file's own ({os.path.dirname(project_file)}).",
                       what=f"working directory {cwd!r} (cwd=) does not exist",
                       allowed="an existing directory",
                       fix=f"Give an existing directory, or leave cwd=None for the project file's own "
                           f"({os.path.dirname(project_file)})")


def require_extrapolations(project_file, cwd):
    """Refuse a session for the extrapolations (``extrapolations=True``) of a project without an active extrapolation
    record, before the engine is asked (it would refuse a run without an active segment)."""
    from .described import session_records
    described = session_records("extrapolations", project_file, cwd)
    if not any(e.active for e in described):
        what = (f"its {len(described)} extrapolation record(s) are all inactive" if len(described)
                else "it has no extrapolation record")
        fix = (f"Add one (pyazr.Project.open({project_file!r}).add_extrapolation(...), then save), or make a session "
               f"for the data sets: pyazr.Session({project_file!r})")
        raise AZURE2Error(
            f"Session(extrapolations=True): {os.path.basename(project_file)} has nothing to extrapolate ({what}); a "
            f"session for the extrapolations needs at least one active extrapolation record.  Add one "
            f"(pyazr.Project.open({project_file!r}).add_extrapolation(...), then save), or make a session for the "
            f"data sets: pyazr.Session({project_file!r}).",
            what=f"{os.path.basename(project_file)} has nothing to extrapolate ({what})",
            allowed="a session for the extrapolations needs at least one active extrapolation record", fix=fix)


def setup_refused(err, said, project_file):
    """The error of a set-up the engine refused (``err``, its ``AZURE2Error``), with the engine's lines of the set-up
    (``said``, :func:`pyazr.engine_text.held`) in its message and as its three parts: what the engine refused (its
    ERROR lines, or its message), what a set-up needs, and where to look."""
    errors = [line.strip() for line in said if line.lstrip("* \t").startswith(("ERROR", "Error"))]
    text = str(err).replace("; check the output above.", ".").replace(" (see the ERROR lines above; ", " (")
    text = text.replace(" (see the ERROR above).", ".")
    name = os.path.basename(project_file)
    what = f"the engine did not set up {name}: " + ("; ".join(errors) if errors else text)
    allowed = "a project whose every set-up check passes (the engine's lines name the one that failed)"
    fix = (f"correct what they name (pyazr.Project.open({project_file!r}) and its edits); Project.validate() lists "
           f"what the validator finds, Project.describe() what the engine makes of each record")
    shown = said[:40] + ([f"... ({len(said) - 40} more lines)"] if len(said) > 40 else [])
    message = (f"{text}\n{what}\nallowed: {allowed}\n{fix}"
               + ("\nThe engine's lines of the set-up:\n  " + "\n  ".join(shown) if shown else ""))
    return type(err)(message, what=what, allowed=allowed, fix=fix)


VECTOR_HINT = ("Session.parameters.x0 is such a vector (the start values); Session.parameters.free names its entries "
               "in order")


def fit_vector(x, n, call, names=None):
    """``x`` as a flat float array, checked to be a free parameter vector of the session: ``n`` entries, every one
    finite.  ``names``: a callable giving the free parameters' names (asked for only to name a bad entry).  The
    refusal is a :class:`~pyazr.ParameterVectorError` with its three parts (the engine makes the same checks, with a
    message only, when a vector reaches it unchecked)."""
    x = np.asarray(x, float).ravel()
    if x.size != n:
        what = f"{call}: a vector of {x.size} values"
        allowed = f"expected the free parameter vector ({n} values, in the order of the fit vector)"
        raise ParameterVectorError(f"{what}; {allowed}, got {x.size}.  {VECTOR_HINT}.", what=what, allowed=allowed,
                                   fix=VECTOR_HINT)
    finite = np.isfinite(x)
    if not finite.all():
        bad = np.flatnonzero(~finite)
        k = int(bad[0])
        name = list(names())[k] if names is not None else f"free parameter {k}"
        more = f" ({bad.size} entries are not finite)" if bad.size > 1 else ""
        what = f"{call}: entry {k} of the free parameter vector ({name}) is {x[k]}{more}"
        allowed = "every entry must be a finite number"
        raise ParameterVectorError(f"{what}; {allowed}.  {VECTOR_HINT}.", what=what, allowed=allowed,
                                   fix=VECTOR_HINT)
    return x


def snapshot_mismatch(what, got, want, rtol):
    """save_fit's refusal of a snapshot whose ``what`` ("chi2", the data chi-squared, or "objective") differs from
    the fit's."""
    cause = ("the objective holds the model terms (abundance priors, Park's wall, measured scattering lengths): a "
             "free target abundance whose prior the engine centres on its start value gives such a difference "
             "(docs/upgrade/pyazr3.md section 1.6).  " if what == "objective" else "")
    found = f"save_fit: the snapshot run on its own gives {what} = {got:.10g}, the fit {want:.10g}"
    allowed = f"relative {rtol:g}"
    fix = "Session.evaluate(x).terms shows the terms; save_fit(..., verify=False) writes the files without this check"
    return AZURE2Error(f"{found} (allowed: {allowed}); it has been removed.  {cause}{fix}.", what=found,
                       allowed=allowed, fix=fix)


def channel_request(pairs, pair, L):
    """``(pair key, L)`` of a request of channel functions (:meth:`pyazr.Session.channel_functions`), checked:
    :class:`~pyazr.ParameterIndexError` for a pair key the model does not have (an ``IndexError`` came from numpy
    before), :class:`~pyazr.DomainError` for a photon pair or an L that is not a whole number >= 0 (zeros came back
    before)."""
    known = {p.key: p for p in pairs}
    try:
        key = int(pair)
    except (TypeError, ValueError):
        key = None
    if key is None or key != pair or key not in known:
        pairs = listed([f'{k} ({q.label})' for k, q in known.items()])
        raise ParameterIndexError(
            f"channel_functions: pair {pair!r} is not a pair key of this model; its pairs are "
            f"{pairs}; Session.pairs lists them (Pair.key).",
            what=f"pair {pair!r} is not a pair key of this model", allowed=f"its pairs are {pairs}",
            fix="Session.pairs lists them (Pair.key)")
    if known[key].is_photon:
        particles = [f"{k} ({q.label})" for k, q in known.items() if not q.is_photon]
        raise DomainError(
            f"channel_functions: pair {key} ({known[key].label}) is a photon pair, which has no Coulomb functions "
            f"(they were returned as zeros before); the particle pairs are {listed(particles)} "
            f"(Session.pairs.particles).",
            what=f"pair {key} ({known[key].label}) is a photon pair, which has no Coulomb functions",
            allowed=f"the particle pairs are {listed(particles)}", fix="Session.pairs.particles")
    try:
        whole = float(L) == int(L) and int(L) >= 0
    except (TypeError, ValueError):
        whole = False
    if not whole:
        raise DomainError(
            f"channel_functions: L = {L!r} is not an orbital angular momentum; L must be a whole number >= 0.  The "
            f"channels of a pair: Session.parameters.of_kind('width') (Parameter.L).",
            what=f"L = {L!r} is not an orbital angular momentum", allowed="L must be a whole number >= 0",
            fix="Session.parameters.of_kind('width') (Parameter.L)")
    return key, int(L)


def warn_engine_notices(engine, stacklevel=4):
    """Warn (:class:`~pyazr.EngineWarning`) with each notice the engine ``engine`` (an ``_azure2.Session``) held
    back since the last call: the file it concerns and the consequence."""
    if engine is None or not hasattr(engine, "take_notices"):
        return
    for n in engine.take_notices():
        warnings.warn(f"{n['file']}: the engine could not write it (\"{n['text']}\").  {n['consequence']}",
                      EngineWarning, stacklevel=stacklevel)
