"""The engine's own fitters on a session (docs/upgrade/fit.md): :meth:`SessionFit.fit_with_engine` runs the program's
fit -- MIGRAD with or without the analytic gradient, Levenberg-Marquardt, the GSL trust region, NLopt where built --
on the session's objective, through the same engine functions as ``AZURE2`` (the table ``pyazr.kinds()["fitters"]``),
with progress records and a stop request another thread can set.  Mixed into the session class
(``class Session(..., SessionFit)``).

Example (T2 of ``tests/upgrade``)::

    s = pyazr.Session("T2.azr")
    r = s.fit_with_engine("migrad")            # the program's MIGRAD, then HESSE
    r.x, r.objective, r.status, r.errors       # free vector, objective (as the program prints it), status, errors
    stop = pyazr.session_fit.StopRequest()     # or a threading.Event
    threading.Timer(10.0, stop.set).start()
    r = s.fit_with_engine("lm", progress=print, stop=stop)
"""
import threading
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import _azure2
from .errors import EngineError
from .locking import engine_call
from .objective import require_data_mode

#: A stop request for :meth:`SessionFit.fit_with_engine`: ``set()`` from any thread ends the fit at its best point.
StopRequest = _azure2.StopRequest


@dataclass
class EngineFitResult:
    """What :meth:`SessionFit.fit_with_engine` returns.

    ``x`` is the free vector at the end (the session's fit basis, the order of ``names``); ``objective`` the objective
    there from one engine evaluation, the number the program prints as ``Total Chi-Squared``; ``status`` how the fitter
    ended (``converged``, ``call_limit``, ``invalid``, ``stopped``, ``no_decrease``, ``iteration_limit``,
    ``failed_<reason>``; NLopt: ``converged_<reason>``); ``converged`` whether it reached a minimum by its criterion;
    ``outcome`` the same in the words every fitter shares (``converged``, ``call_limit``, ``iteration_limit``,
    ``stopped``, ``invalid``, ``failed``; ``fitStatus.out``'s ``status``); ``x`` and ``objective`` are the lowest
    objective the fitter reached: a MIGRAD that stops at its call limit ends at the lowest objective it evaluated
    (a notice says so when that is not Minuit2's last state);
    ``edm`` MIGRAD's estimated distance to the minimum (NaN for the other fitters); ``calls`` the evaluations of the
    objective (LM, GSL-LM: residual passes); ``covariance`` of the free parameters (HESSE after MIGRAD,
    (J^T J + priors)^-1 after LM and GSL-LM; None without one, e.g. NLopt or a stopped fit), ``errors`` the errors
    the program writes into param.sav (the square roots of the diagonal after LM and GSL-LM; Minuit's parabolic errors
    after MIGRAD, which differ from them for a parameter with limits); ``minos`` name -> (lower, upper) when MINOS ran;
    ``minos_valid`` name -> {"lower": side, "upper": side} beside it, each side Minuit2's verdict as a dict: ``valid``
    (the crossing was found), ``at_limit`` (the scan reached the parameter's limit), ``call_limit`` (it ran out of
    calls), ``new_minimum`` (it found a lower point: the fit was not at its minimum); ``notices`` what the caller
    should know, one sentence each (a fallback to another fitter, a setting the fitter ignored such as ``max_calls``
    on LM, the best point kept at a call limit); ``log`` the fitter's lines.
    """
    method: str
    x: np.ndarray
    objective: float
    status: str
    converged: bool
    edm: float
    calls: int
    covariance: Optional[np.ndarray]
    errors: Optional[np.ndarray]
    minos: Dict[str, Tuple[float, float]]
    names: List[str]
    x_full: np.ndarray = field(repr=False, default=None)
    covariance_from_hesse: bool = False
    fell_back_to_migrad: bool = False
    outcome: str = ""
    minos_valid: Dict[str, Dict[str, Dict[str, bool]]] = field(default_factory=dict)
    notices: List[str] = field(default_factory=list)
    log: str = field(repr=False, default="")


def _stop_request(stop, progress):
    """The engine's StopRequest for ``stop`` (None, a StopRequest, or anything with ``is_set()`` such as a
    threading.Event), the progress callable to hand to the engine, and the function that ends the watching.

    The engine reads its own StopRequest only.  Anything else is looked at after every progress record, in the
    record's thread, so that the evaluations running when it was set are the last ones; a helper thread looks too
    (every 20 ms), for the time in which no record comes.  The helper alone let further evaluations through where
    one takes less than its period (seen on the CI's macOS machine: two records after the request)."""
    if stop is None or isinstance(stop, StopRequest):
        return stop, progress, lambda: None
    if not callable(getattr(stop, "is_set", None)):
        raise TypeError("stop must be a pyazr.session_fit.StopRequest or a threading.Event (has is_set()).")
    request = StopRequest()
    done = threading.Event()

    def record(rec):
        try:
            if progress is not None:
                progress(rec)
        finally:
            if stop.is_set():
                request.set()

    def watch():
        while not done.is_set():
            if stop.is_set():
                request.set()
                return
            done.wait(0.02)

    watcher = threading.Thread(target=watch, daemon=True)
    watcher.start()

    def finish():
        done.set()
        watcher.join()
    return request, record, finish


class SessionFit:
    """Mixed into :class:`pyazr.Session`: uses its ``_engine``, ``mode`` and ``parameters``."""

    @engine_call
    def fit_with_engine(self, method="migrad", minos=False, x0=None, max_calls=None, progress=None, stop=None):
        """The program's fit on this session's objective, through the engine's own fitters.

        Parameters
        ----------
        method : a fitter of ``pyazr.kinds()["fitters"]`` (the engine's table: each row gives the program's switch
            ``option``, what the fitter yields and its ``doc``; a build without NLopt has no ``"nlopt"`` row).  The
            default ``"migrad"`` is the row the program runs without a switch.  An unknown name is refused with the list.  The program
            runs ``migrad-gradient`` for a project with THM data when no switch is given; here the method is the one
            named.
        minos : run MINOS after the fit (Up = the session's ``chiVariance``, 1 by default); for a fitter whose row
            says ``minos`` only.
        x0 : the start fit vector (default ``parameters.x0``).
        max_calls : the call limit (``--max-calls``) of a fitter that honours one (``call_limit`` of its row in
            ``pyazr.kinds()["fitters"]``); the others ignore it and say so in ``notices``.  None: no
            limit.
        progress : called with a dict per progress record, from the engine's threads with the GIL taken: ``fitter``,
            ``stage`` ("fit", then "hesse" and "minos" for the error analysis after a MIGRAD row; each stage counts its
            calls and its best from the start), ``kind`` ("evaluation" or "iteration"), ``evaluation``, ``iteration``, ``calls``, ``objective``, ``best``,
            ``measure`` (LM: the damping; NaN where none), ``error_analysis``, ``x`` (the free vector of the record).
            It must not call the session.  An exception in it stops the fit and is raised afterwards.
        stop : a :data:`StopRequest` or a ``threading.Event``; set from another thread, the fit ends at the best point
            it has seen with status ``"stopped"`` (the evaluations running then finish; every later call returns at
            once) and no covariance.

        Returns an :class:`EngineFitResult`.  Writes no file (``save_fit`` and ``write_outputs`` do) and does
        not change the session's parameter values; the session's evaluation buffer holds the end point's model
        values.  The fit applies ``<parameterSettings>`` (limits, nuisance priors) as the program does.  After a
        MIGRAD row HESSE refines the covariance (the program runs it with ``--covariance-band``).  A MIGRAD row on a
        data covariance with relvector terms re-freezes its model at the fitted point, as the program does.
        """
        require_data_mode(self, "fit_with_engine()")
        rows = {row["name"]: row for row in _azure2.fitters()}
        if method not in rows:
            raise ValueError(f"unknown fit method {method!r}; the engine's fitters: {', '.join(rows)}.")
        if minos and not rows[method]["minos"]:
            raise ValueError(f"MINOS follows a Minuit2 minimum only; {method!r} has none "
                             f"(use one of {', '.join(n for n, r in rows.items() if r['minos'])}).")
        free = self.parameters.free
        names = list(free.names)
        start = [] if x0 is None else [float(v) for v in np.asarray(x0, float).ravel()]
        if x0 is not None and len(start) != len(names):
            raise ValueError(f"x0 has {len(start)} values; the session has {len(names)} free parameters.")
        calls = 0 if max_calls is None else int(max_calls)
        if calls < 0:
            raise ValueError("max_calls must be a whole number >= 0 (None: no limit).")
        request, record, finish = _stop_request(stop, progress)
        try:
            out = self._engine.fit_with_engine(method, bool(minos), start, calls, record, request)
            fixed = np.asarray(self._engine.params_fixed(), bool)
        except ValueError as e:  # the engine's refusal (std::invalid_argument)
            raise EngineError(str(e)) from e
        finally:
            finish()
        x_full = np.asarray(out["x_full"], float)
        idx = np.asarray(out["free_indices"], int)
        x = x_full[~fixed] if x_full.size == fixed.size else x_full
        covariance = errors = None
        cov = np.asarray(out["covariance"], float)
        if cov.size and idx.size and cov.size == idx.size * idx.size:
            full_cov = cov.reshape(idx.size, idx.size)
            # rows in the order of the session's free vector
            position = {int(f): k for k, f in enumerate(np.flatnonzero(~fixed))}
            order = [position.get(int(f), -1) for f in idx]
            if all(o >= 0 for o in order) and len(order) == len(names):
                covariance = np.zeros((len(names), len(names)))
                for a, oa in enumerate(order):
                    for b, ob in enumerate(order):
                        covariance[oa, ob] = full_cov[a, b]
        errors_full = np.asarray(out["errors_full"], float)
        if covariance is not None and errors_full.size == fixed.size:
            errors = errors_full[~fixed]
        intervals, valid = {}, {}
        if out["minos"]:
            free_positions = np.flatnonzero(~fixed)
            for k, f in enumerate(free_positions):
                if f < len(out["minos"]):
                    intervals[names[k]] = tuple(float(v) for v in out["minos"][f])
                if f < len(out["minos_validity"]):
                    valid[names[k]] = {side: dict(flags) for side, flags in out["minos_validity"][f].items()}
        return EngineFitResult(method=out["method"], x=np.asarray(x, float), objective=float(out["objective"]),
                               status=out["status"], converged=bool(out["converged"]), edm=float(out["edm"]),
                               calls=int(out["calls"]), covariance=covariance, errors=errors, minos=intervals,
                               names=names, x_full=x_full,
                               covariance_from_hesse=bool(out["covariance_from_hesse"]),
                               fell_back_to_migrad=bool(out["fell_back_to_migrad"]), outcome=out["outcome"],
                               minos_valid=valid, notices=list(out["notices"]), log=out["log"])
