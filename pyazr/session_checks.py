"""The argument checks and engine notices of a calculation session (:class:`pyazr.azure2`).

Each check raises one of pyazr's errors (:mod:`pyazr.errors`) with a message in three parts: what was wrong with
the value given, what is allowed, which call shows or fixes it.  :func:`warn_engine_notices` turns the notices the
binding held back (``_azure2.Session.take_notices``, ``api/src/pybind_notices.cpp``) into warnings.
"""
import os
import warnings

from .errors import (AZURE2Error, DomainError, EngineWarning, MissingFile, ParameterIndexError, SessionClosed,
                     listed)


def session_closed(project_file, caller):
    """The error of a call on a closed session; ``caller`` is the name of the session's method that was called."""
    call = f"azure2.{caller}" if caller not in ("<module>", "<lambda>") else "this call"
    return SessionClosed(
        f"{call}: this azure2 session of {os.path.basename(project_file)} is closed (close() or the end of its with "
        f"block released the engine); a closed session cannot be reopened.  Open a new one: "
        f"pyazr.azure2({project_file!r}); is_alive() tells whether a session is open.")


def missing_directory(cwd, project_file):
    """The error of a working directory (``cwd=``) that does not exist."""
    return MissingFile(cwd, f"working directory {cwd!r} (cwd=) does not exist; the project's relative paths are read "
                            f"and written there.  Give an existing directory, or leave cwd=None for the project "
                            f"file's own ({os.path.dirname(project_file)}).")


def missing_parameter_file(path, call):
    """The error of a parameter file (``param.sav``) that a call updates and that does not exist."""
    return MissingFile(path, f"{call}: {path} does not exist; this call updates the rows of an existing param.sav "
                             f"(the program's fit writes one).  azure2.save_fit(path) writes a project and its .sav "
                             f"from a free vector without one.")


def snapshot_mismatch(what, got, want, rtol):
    """save_fit's refusal of a snapshot whose ``what`` ("chi2" or "priors") differs from the fit's."""
    cause = ("calculate_chi2_rwa includes the model terms (abundance priors, Park's wall, measured scattering "
             "lengths): a free target abundance whose prior the engine centres on its start value gives such a "
             "difference (docs/upgrade/pyazr3.md section 1.6).  " if what == "chi2" else "")
    return AZURE2Error(f"save_fit: the snapshot run on its own gives {what} = {got:.10g}, the fit {want:.10g} "
                       f"(allowed: relative {rtol:g}); it has been removed.  {cause}azure2.penalties(x) shows the "
                       f"terms; save_fit(..., verify=False) writes the files without this check.")


def channel_request(pairs, pair, L):
    """``(pair key, L)`` of a request of channel functions (:meth:`pyazr.azure2.coulomb_functions`), checked:
    :class:`~pyazr.ParameterIndexError` for a pair key the model does not have (an ``IndexError`` came from numpy
    before), :class:`~pyazr.DomainError` for a photon pair or an L that is not a whole number >= 0 (zeros came back
    before)."""
    known = {p.key: p for p in pairs}
    try:
        key = int(pair)
    except (TypeError, ValueError):
        key = None
    if key is None or key != pair or key not in known:
        raise ParameterIndexError(
            f"coulomb_functions: pair {pair!r} is not a pair key of this model; its pairs are "
            f"{listed([f'{k} ({q.label})' for k, q in known.items()])}; azure2.pairs lists them (Pair.key).")
    if known[key].is_photon:
        particles = [f"{k} ({q.label})" for k, q in known.items() if not q.is_photon]
        raise DomainError(
            f"coulomb_functions: pair {key} ({known[key].label}) is a photon pair, which has no Coulomb functions "
            f"(they were returned as zeros before); the particle pairs are {listed(particles)} "
            f"(azure2.pairs.particles).")
    try:
        whole = float(L) == int(L) and int(L) >= 0
    except (TypeError, ValueError):
        whole = False
    if not whole:
        raise DomainError(
            f"coulomb_functions: L = {L!r} is not an orbital angular momentum; L must be a whole number >= 0 (the "
            f"result was all zeros before).  The channels of a pair: azure2.parameters.widths (Parameter.L).")
    return key, int(L)


def warn_engine_notices(engine, stacklevel=4):
    """Warn (:class:`~pyazr.EngineWarning`) with each notice the engine ``engine`` (an ``_azure2.Session``) held
    back since the last call: the file it concerns and the consequence."""
    if engine is None or not hasattr(engine, "take_notices"):
        return
    for n in engine.take_notices():
        warnings.warn(f"{n['file']}: the engine could not write it (\"{n['text']}\").  {n['consequence']}",
                      EngineWarning, stacklevel=stacklevel)
