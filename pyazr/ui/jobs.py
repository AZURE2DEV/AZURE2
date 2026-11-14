"""Jobs of the web UI (docs/upgrade/web_ui.md section 4.4): the job table, the recipes it runs, and the job records.

A job kind is one row of :data:`JOBS`: its kind and method, the function that runs it (or the pyazr call it waits
for: the route answers 501 naming it), the package it needs, its settings with their defaults, and its doc.  A job
runs in a worker process of its own (``python -m pyazr.ui.workers job``) on its own session of the text of the
revision it was started on; it reports progress after every objective call (calls, objective, best objective, elapsed
time, best x by name), stops cooperatively before the next call when asked, and writes a checkpoint
``<ui dir>/jobs/<id>.json`` (state, progress, best x) on every improvement and at its end, from which a new job can
start (``x0: {"job": id}``) and which a restarted server lists as ``interrupted``.

A row's settings are rows of :class:`Setting` (name, type, default, doc): the page draws the job form from them and the
server converts a request's values by the type's row of :data:`SETTING_TYPES`.  The server keeps a decimated trace of
the objective against the call number (at most :data:`TRACE_POINTS` points) and serves the last lines of the job's
log beside its record; ``GET /api/jobs/{id}?wait=s`` answers when the job's state changes or after ``s`` seconds.
"""
import json
import os
import threading
import time
import uuid
from dataclasses import dataclass
from functools import lru_cache, partial

import numpy as np

from .errors import Conflict, Missing, NotFound, Refused, UIError, describe_exception
from .utilities import rate_engine

FINAL = ("stopped", "converged", "ended", "finished", "failed", "killed", "interrupted")
TRACE_POINTS = 500  # the most points of a job's trace (objective against call); older points are thinned by two
LOG_LINES = 8  # the last lines of a job's log served with its record
WAIT_MAX = 10.0  # s, the longest a GET of a job waits for its state to change


class StopRequested(Exception):
    """Raised inside the objective when a stop was asked: the recipe ends with the best point seen."""


class Tracker:
    """The job's handle for its function: wraps the session's objective-and-gradient for a recipe (reads the stop
    request before each call, counts it, records the best point and reports it with ``report(message)``; ``None`` in
    a direct call), and counts the progress of calls the job makes itself (:meth:`count`: the engine's fitters, the
    sampler).  ``session`` is the job's session adapter (:class:`pyazr.ui.engine.Engine`), for a job that calls the
    engine other than through the objective; ``source`` the record of the job whose best point the job starts from
    (its result: a covariance, a chain), or None; :meth:`product` the path of a file the job writes."""

    def __init__(self, fg, names, *, report=None, stop=None, session=None, source=None, job_id="direct",
                 products="."):
        self.fg, self.names, self.report, self.session = fg, names, report, session
        self.source, self.job_id, self.products = source, job_id, products
        self.stop = stop if stop is not None else threading.Event()  # set: the job ends with its best point
        self.calls, self.best, self.best_x, self.t0 = 0, np.inf, None, time.perf_counter()
        self.first = None  # the objective at the start point
        self.after_stop = 0  # progress reports after the stop request was seen (the calls running then)
        self._last = None  # (x bytes, f, g) of the last call, for recipes that ask f and g separately

    def __call__(self, x):
        """(f, g) at x."""
        x = np.asarray(x, float)
        key = x.tobytes()
        if self._last is not None and self._last[0] == key:
            return self._last[1], self._last[2]
        if self.stop_requested():
            raise StopRequested()
        f, g = self.fg(x)
        f = float(f)
        self._last = (key, f, np.asarray(g, float))
        self.count(self.calls + 1, f, x)
        return f, self._last[2]

    def count(self, calls, objective, x, **more):
        """Progress: ``calls`` so far (counting up), the objective at the free vector ``x`` of the last of them (None
        where the call computes none: a model on an extrapolation grid), and further numbers by name (``more``: a
        fitter's stage, a sampler's acceptance; None where there is none), reported with the best point when it
        improved."""
        objective = None if objective is None else float(objective)
        if self.stop.is_set():
            self.after_stop += 1
        self.calls = int(calls)
        if self.first is None and objective is not None:
            self.first = objective
        improved = objective is not None and objective < self.best
        if improved:
            self.best, self.best_x = objective, np.array(x, float)
        if self.report is not None:
            self.report({"event": "call", "calls": self.calls, "objective": objective, "best": self.best,
                         "elapsed": time.perf_counter() - self.t0,
                         "best_x": dict(zip(self.names, self.best_x.tolist())) if improved else None, **more})

    def stop_requested(self):
        return self.stop.is_set()

    def product(self, name, suffix):
        """The path of a file the job writes: ``<products>/<name>_<job id><suffix>``."""
        os.makedirs(self.products, exist_ok=True)
        return os.path.join(self.products, f"{name}_{self.job_id}{suffix}")

    def source_product(self, name, suffix):
        """The path of the file ``name`` the source job wrote (its :meth:`product`), or None where it wrote none: a
        job killed or interrupted (a server restart) has no result, but its files as last saved."""
        if not self.source:
            return None
        path = os.path.join(self.products, f"{name}_{self.source['id']}{suffix}")
        return path if os.path.isfile(path) else None

    def value(self, x):
        return self(x)[0]

    def grad(self, x):
        return self(x)[1]


# -- the recipes: functions of (tracker, x0, settings, names) -> (status, message, extra) -------------------------

def fit_scipy_lbfgsb(tracker, x0, settings, names):
    """scipy.optimize.minimize, L-BFGS-B, on the engine's objective and gradient (pyazr3.md 2.6)."""
    import scipy.optimize as so
    r = so.minimize(tracker, x0, jac=True, method="L-BFGS-B", options={"maxiter": int(settings["max_iterations"])})
    return ("converged" if r.success else "ended"), str(r.message), {"iterations": int(r.nit)}


def fit_iminuit_migrad(tracker, x0, settings, names):
    """iminuit's MIGRAD with the engine's gradient, then HESSE (errors, covariance and correlation by name), errordef 1
    (a chi-squared)."""
    m = _minuit(tracker, x0, names)
    m.migrad(ncall=settings["max_calls"])
    extra = {"edm": float(m.fmin.edm), "valid": bool(m.valid)}
    if not m.valid:
        return "ended", "MIGRAD: invalid minimum", extra
    if not settings["hesse"]:
        return "converged", "MIGRAD: valid minimum", extra
    m.hesse()  # a valid minimum can turn invalid here (the EDM recomputed with HESSE's matrix)
    extra.update(_errors_of(m, names), edm=float(m.fmin.edm), valid=bool(m.valid))
    extra["errors_from"] = "HESSE (iminuit)"
    if m.valid:
        return "converged", "MIGRAD: valid minimum; HESSE: still valid", extra
    return "ended", f"MIGRAD: valid minimum; after HESSE invalid (EDM {m.fmin.edm:.3g})", extra


def minos_iminuit(tracker, x0, settings, names):
    """iminuit's MINOS from a fit's best point: MIGRAD (a few calls at a minimum), HESSE, then MINOS on the parameters
    named (all free ones when none): asymmetric intervals {name: {lower, upper, valid}} in the fit basis."""
    m = _minuit(tracker, x0, names)
    m.migrad(ncall=settings["max_calls"])
    extra = {"edm": float(m.fmin.edm), "valid": bool(m.valid)}
    if not m.valid:
        return "ended", "MIGRAD from the start point: invalid minimum; MINOS not run", extra
    m.hesse()
    extra.update(_errors_of(m, names), errors_from="HESSE (iminuit)")
    which = settings["parameters"] or list(names)
    m.minos(*which, ncall=settings["max_calls"])
    extra["minos"] = {n: {"lower": float(m.merrors[n].lower), "upper": float(m.merrors[n].upper),
                          "valid": bool(m.merrors[n].is_valid)} for n in which}
    bad = [n for n in which if not m.merrors[n].is_valid]
    message = "MINOS: every interval valid" if not bad else f"MINOS: invalid intervals of {', '.join(bad)}"
    return ("converged" if not bad else "ended"), message, extra


def _minuit(tracker, x0, names):
    import iminuit
    m = iminuit.Minuit(tracker.value, np.asarray(x0, float), grad=tracker.grad, name=names)
    m.errordef = iminuit.Minuit.LEAST_SQUARES
    return m


def _errors_of(m, names):
    """HESSE's errors by name, and its covariance and correlation matrices in the order of ``names``."""
    return {"errors": dict(zip(names, [float(e) for e in m.errors])),
            "covariance": {"names": list(names), "matrix": np.asarray(m.covariance, float).tolist()},
            "correlation": {"names": list(names), "matrix": np.asarray(m.covariance.correlation(), float).tolist()}}


def _covariance_by_name(names, matrix):
    """The covariance and the correlation ({names, matrix}) of a covariance matrix in the order of ``names``."""
    c = np.asarray(matrix, float)
    sd = np.sqrt(np.clip(np.diag(c), 0.0, None))
    scale = np.where(sd > 0, sd, 1.0)
    return {"covariance": {"names": list(names), "matrix": c.tolist()},
            "correlation": {"names": list(names), "matrix": (c / np.outer(scale, scale)).tolist()}}


# The job's state after an engine fit, by the engine's status (any other status: "ended", with the status said).
ENGINE_STATES = {"converged": "converged", "stopped": "stopped"}


def _fitters():
    import pyazr
    return pyazr.kinds()["fitters"]


def fit_engine(tracker, x0, settings, names, *, method, minos=None):
    """The engine's fitter ``method`` (``Session.fit_with_engine``, the program's fit) from x0: each progress record
    counted into the tracker (the calls of the fit, then of HESSE and MINOS, counted on), the tracker's stop request
    passed as the fit's; MINOS when ``minos`` (or the setting ``minos``).  The job's best point is the fitter's end
    point and its objective the engine's there; the start objective is one evaluation at x0 after the fit (an
    evaluation before it would change the state the engine's pooled working copies carry into the fit).  The
    fitter's lines go to the job's log."""
    from .engine import fit_with_engine, objective_of
    stage = {"name": None, "offset": 0, "last": 0}

    def progress(record):
        if record["stage"] != stage["name"]:  # HESSE and MINOS count their calls from 1: the job counts on
            stage.update(name=record["stage"], offset=stage["offset"] + stage["last"])
        stage["last"] = record["calls"]
        measure = record["measure"]
        tracker.count(stage["offset"] + record["calls"], record["objective"], record["x"], stage=record["stage"],
                      measure=float(measure) if np.isfinite(measure) else None)

    run_minos = settings.get("minos", False) if minos is None else minos
    r = fit_with_engine(tracker.session, method, run_minos, x0, settings["max_calls"], progress, tracker.stop)
    print(r["log"], flush=True)
    tracker.first = float(objective_of(tracker.session)(x0))
    tracker.best, tracker.best_x = r["objective"], np.array([r["x"][n] for n in names], float)
    # the fitter that ran in place of the one asked for where a data set has no analytic Jacobian: the fitter table's
    # "fallback" of the row asked for
    fallback = _fitters()[r["method"]]["fallback"]
    ran = fallback if r["fell_back_to_migrad"] else r["method"]
    words = f"the engine's {ran}: {r['status']}"
    # the fit's notices (a fallback to another fitter, a call limit the fitter ignored) and its outcome in the words
    # every fitter shares (converged, call_limit, ...), as the engine gives them
    extra = {"engine_status": r["status"], "outcome": r["outcome"], "valid": r["converged"], "edm": r["edm"],
             "fit_calls": r["calls"], "fitter": {"asked": r["method"], "ran": ran}, "notices": list(r["notices"])}
    if r["fell_back_to_migrad"]:
        words += (f" ({fallback} in place of {r['method']}: no analytic Jacobian for a data set; no HESSE, "
                  "no MINOS)")
    if r["covariance"] is not None:
        extra.update(_covariance_by_name(r["names"], r["covariance"]), errors=r["errors"],
                     errors_from="HESSE after the engine's MIGRAD" if r["covariance_from_hesse"]
                     else f"the covariance of the {fallback} that ran in place of {r['method']}"
                     if r["fell_back_to_migrad"] else "the engine's (J^T J + priors)^-1")
    if r["minos"]:
        # each side's verdict as the engine gives it ({valid, at_limit, call_limit, new_minimum}); valid: both sides
        sides = r["minos_valid"]
        extra["minos"] = {n: {"lower": lo, "upper": hi,
                              "valid": all(v.get("valid") for v in sides[n].values()) if n in sides else None,
                              "sides": sides.get(n)} for n, (lo, hi) in r["minos"].items()}
    return ENGINE_STATES.get(r["status"], "ended"), words, extra


# -- the sampler and the bands -------------------------------------------------------------------------------------

QUANTILES = (2.5, 16.0, 50.0, 84.0, 97.5)  # %, of each parameter's marginal in a sampler's summary
TAU_EVERY = 10  # stored steps between two running estimates of the autocorrelation time
SAVE_EVERY = 50  # stored steps between two writes of the chain while it runs


def _start_walkers(log_prob, x0, bounds, names, walkers, spread, rng):
    """Walkers around x0, each ``x0 + spread * |x0| * N(0, 1)`` (``spread`` itself where x0 is 0), drawn again until
    it lies inside the fit's limits with a finite log-probability.  A start outside the limits is refused with the
    parameters named."""
    from pyazr.errors import ParameterError
    lower, upper = np.array(bounds, float).T
    outside = [f"{n} = {v:g} (limits [{lo:g}, {hi:g}])" for n, v, lo, hi in zip(names, x0, lower, upper)
               if not lo <= v <= hi]
    if outside:
        raise ParameterError(f"the sampler's start point is outside the fit's limits: {'; '.join(outside)}; start "
                             "from a point inside them (a fit's best point, or the working x moved inside)")
    width = spread * np.where(x0 != 0.0, np.abs(x0), 1.0)
    p0 = np.empty((walkers, len(x0)))
    for k in range(walkers):
        for _ in range(1000):
            p0[k] = x0 + width * rng.standard_normal(len(x0))
            if np.all(p0[k] >= lower) and np.all(p0[k] <= upper) and np.isfinite(log_prob(p0[k])):
                break
        else:
            raise ParameterError(f"no start inside the fit's limits near x0 for walker {k} after 1000 draws: make "
                                 "the spread smaller")
    return p0


def _autocorrelation_time(chain):
    """emcee's integrated autocorrelation time of each parameter (in stored steps) of ``chain`` (steps, walkers,
    parameters), without its test of the chain's length (a running estimate: the result says the length)."""
    import emcee
    return np.asarray(emcee.autocorr.integrated_time(chain, tol=0, quiet=True), float)


def _split_rhat(chain):
    """The split R-hat of each parameter (Gelman et al., Bayesian Data Analysis, 3rd ed., eq. 11.4): every walker's
    series cut in two halves taken as chains; sqrt(((m-1)/m W + B/m) / W) with W the mean of the halves' variances and
    B/m the variance of their means.  The walkers of an ensemble sampler are not independent chains: a value near 1 is
    necessary, not sufficient."""
    m = chain.shape[0] // 2
    halves = np.concatenate([chain[:m], chain[m:2 * m]], axis=1)  # (m, 2 walkers, parameters)
    within = halves.var(axis=0, ddof=1).mean(axis=0)
    between = halves.mean(axis=0).var(axis=0, ddof=1)
    return np.sqrt(((m - 1) / m * within + between) / within)


def chain_summary(chain, names, discard):
    """What the page receives of a chain (steps, walkers, parameters) after ``discard`` stored steps: per parameter
    the mean, standard deviation, quantiles (:data:`QUANTILES`), autocorrelation time, effective draws and split
    R-hat; the covariance and correlation of the draws by name.  No draw."""
    kept = chain[discard:]
    flat = kept.reshape(-1, kept.shape[-1])
    tau = _autocorrelation_time(kept)
    rhat = _split_rhat(kept) if kept.shape[0] >= 4 else np.full(len(names), np.nan)
    q = np.percentile(flat, QUANTILES, axis=0)
    out = {"draws_used": int(flat.shape[0]), "discard": int(discard), "steps_stored": int(chain.shape[0]),
           "walkers": int(chain.shape[1]),
           "mean": dict(zip(names, flat.mean(axis=0).tolist())),
           "errors": dict(zip(names, flat.std(axis=0, ddof=1).tolist())),
           "quantiles": {n: {f"{p:g}": float(q[i, j]) for i, p in enumerate(QUANTILES)} for j, n in enumerate(names)},
           "tau": dict(zip(names, tau.tolist())),
           "effective_draws": dict(zip(names, (flat.shape[0] / np.where(tau > 0, tau, np.nan)).tolist())),
           "rhat": dict(zip(names, rhat.tolist())),
           "errors_from": "the chain's standard deviations after its burn-in"}
    out.update(_covariance_by_name(names, np.atleast_2d(np.cov(flat, rowvar=False))))
    return out


def sample_emcee(tracker, x0, settings, names):
    """emcee's ensemble sampler on the objective's log-probability (``Objective.log_prob``: -0.5 objective, -inf
    outside the fit's limits), one process, stopping between stored steps.  Walkers start around x0 inside the
    limits, or, from a sampler job, where its chain ended (the chain continues: the new file holds the old steps and
    the new).  The chain (stored steps, walkers, parameters) goes to ``chain_<id>.npy`` in the project's output
    directory every :data:`SAVE_EVERY` stored steps and at the end, the summary (:func:`chain_summary`) beside it as
    JSON; the page receives the summary.  Progress per stored step: the objective at the best walker, the step, the
    mean acceptance and a running autocorrelation time."""
    import emcee
    from pyazr.errors import ParameterError
    from .engine import objective_of
    f = objective_of(tracker.session)
    ndim = len(names)
    calls = {"n": 0}

    def log_prob(x):
        calls["n"] += 1
        return f.log_prob(x)

    previous = None
    # the source's chain: its result's, or the file a killed or interrupted sampler saved last (every SAVE_EVERY)
    chain_of_source = ((tracker.source or {}).get("result") or {}).get("chain") or tracker.source_product("chain",
                                                                                                          ".npy")
    if chain_of_source:
        previous = np.load(chain_of_source)
        if previous.shape[2] != ndim:
            raise ParameterError(f"the chain of job {tracker.source['id']} has {previous.shape[2]} parameters, this "
                                 f"session {ndim}: a chain continues on the text it was sampled on")
        p0 = previous[-1]
    else:
        walkers = settings["walkers"] or 2 * ndim + 2
        if walkers < 2 * ndim:
            raise ParameterError(f"walkers {walkers}: emcee needs at least twice the free parameters ({2 * ndim})")
        p0 = _start_walkers(f.log_prob, np.asarray(x0, float), f.bounds, names, walkers, settings["spread"],
                            np.random.default_rng(settings["seed"]))
    sampler = emcee.EnsembleSampler(len(p0), ndim, log_prob)
    sampler.random_state = np.random.RandomState(settings["seed"]).get_state()
    chain_path, summary_path = tracker.product("chain", ".npy"), tracker.product("chain", ".json")

    def whole():
        new = sampler.get_chain()
        return new if previous is None else np.concatenate([previous, new])

    start_steps = 0 if previous is None else previous.shape[0]
    tau, stopped = None, False
    for k, state in enumerate(sampler.sample(p0, iterations=settings["steps"], thin_by=settings["thin"]), 1):
        if k % TAU_EVERY == 0 or k == settings["steps"]:
            tau = float(np.max(_autocorrelation_time(whole())))
        if k % SAVE_EVERY == 0:
            np.save(chain_path, whole())
        best = int(np.argmax(state.log_prob))
        tracker.count(calls["n"], -2.0 * state.log_prob[best], state.coords[best], step=start_steps + k,
                      steps=start_steps + settings["steps"], acceptance=float(np.mean(sampler.acceptance_fraction)),
                      tau=tau)
        if tracker.stop_requested():
            stopped = True
            break
    chain = whole()
    np.save(chain_path, chain)
    discard = chain.shape[0] // 2 if settings["discard"] is None else min(settings["discard"], chain.shape[0] - 1)
    summary = chain_summary(chain, names, discard)
    summary.update(chain=chain_path, acceptance=float(np.mean(sampler.acceptance_fraction)), thin=settings["thin"],
                   continues=None if previous is None else tracker.source["id"],
                   continued_from=None if previous is None else chain_of_source)
    write_checkpoint(summary_path, plain_json(dict(summary, names=list(names))))
    summary["summary"] = summary_path
    taus = [t for t in summary["tau"].values() if np.isfinite(t) and t > 0]
    tau_words = f"up to {max(taus):.3g} steps" if taus else "not estimated (too few steps)"
    words = (f"{chain.shape[0]} stored steps of {chain.shape[1]} walkers ({start_steps} of them continued), "
             f"acceptance {summary['acceptance']:.2f}, autocorrelation time {tau_words} after a burn-in of {discard}")
    return ("stopped" if stopped else "finished"), words, summary


def _band_entries(points, values, lower, upper, how):
    """The band of each data set or extrapolation: the points' coordinates, the central value, lower and upper
    edges in the engine's unit, and the same times the S-factor conversion where it is offered."""
    out = []
    for p, v, lo, hi in zip(points, values, lower, upper):
        entry = {k: p[k] for k in ("name", "e_cm", "e_x", "angle_cm")}
        entry.update(value=np.asarray(v).tolist(), lower=np.asarray(lo).tolist(), upper=np.asarray(hi).tolist())
        if p["conv"] is not None:
            conv = np.asarray(p["conv"])
            entry["s_factor"] = {"value": (v * conv).tolist(), "lower": (lo * conv).tolist(),
                                 "upper": (hi * conv).tolist()}
        out.append(entry)
    return {"how": how, "sets": out}


def bands_covariance(tracker, x0, settings, names):
    """The 1-sigma band of the model of every data set and extrapolation at a fit's best point from its covariance
    (linear: model +- sigma, sigma^2 = g^T C g, ``pyazr.bands``)."""
    from .engine import band_points, linear_sigma
    cov = tracker.source["result"]["covariance"]
    points = band_points(tracker.session, x0)
    sigma_data, sigma_grids, varied = linear_sigma(tracker.session, x0, cov["matrix"], cov["names"])
    every = points["datasets"] + points["extrapolations"]
    values = [np.asarray(p["model"]) for p in every]
    sigma = sigma_data + sigma_grids
    tracker.count(1, None, x0)
    how = (f"linear, from the covariance of job {tracker.source['id']} at its best point: the model +- sigma, "
           "sigma^2 = g^T C g over the level parameters (pyazr.bands)")
    bands = _band_entries(every, values, [v - s_ for v, s_ in zip(values, sigma)],
                          [v + s_ for v, s_ in zip(values, sigma)], how)
    return "finished", how, {"bands": bands, "source": tracker.source["id"], "varied": varied,
                             "free": list(names)}


def bands_draws(tracker, x0, settings, names):
    """The band of the model of every data set and extrapolation from draws of a sampler's chain (after its burn-in):
    the 16th, 50th and 84th percentiles of the model over the draws at each point.  Stops between draws."""
    from pyazr.errors import ParameterError
    from .engine import band_points, models_at
    result = tracker.source["result"]
    chain = np.load(result["chain"])
    if chain.shape[2] != len(names):
        raise ParameterError(f"the chain of job {tracker.source['id']} has {chain.shape[2]} parameters, this session "
                             f"{len(names)}: a band from draws needs the text the chain was sampled on")
    flat = chain[result["discard"]:].reshape(-1, chain.shape[2])
    rng = np.random.default_rng(settings["seed"])
    xs = flat[rng.choice(len(flat), size=min(settings["draws"], len(flat)), replace=False)]

    def progress(k, objective, x):
        if tracker.stop_requested():
            raise StopRequested()
        tracker.count(k, objective, x)

    points = band_points(tracker.session, x0)
    data, grids = models_at(tracker.session, xs, progress)
    every = points["datasets"] + points["extrapolations"]
    models = data + grids
    lower, middle, upper = zip(*[np.percentile(m, [16.0, 50.0, 84.0], axis=0) for m in models])
    how = (f"from {len(xs)} draws of the chain of job {tracker.source['id']} after its burn-in of "
           f"{result['discard']} steps: the 16th, 50th and 84th percentiles of the model at each point")
    return "finished", how, {"bands": _band_entries(every, middle, lower, upper, how), "source": tracker.source["id"],
                             "draws": len(xs), "varied": list(names), "free": list(names)}


# -- settings: their types and the conversion of a request's value ------------------------------------------------

def _as_integer(value):
    if value is None or value == "":
        return None
    number = float(value)
    if not number.is_integer():
        raise ValueError(f"{value!r} is not an integer")
    return int(number)


def _as_boolean(value):
    if isinstance(value, bool):
        return value
    if str(value).lower() in ("true", "1", "yes"):
        return True
    if str(value).lower() in ("false", "0", "no"):
        return False
    raise ValueError(f"{value!r} is not true or false")


def _as_number(value):
    if value is None or value == "":
        return None
    number = float(value)
    if not np.isfinite(number):
        raise ValueError(f"{value!r} is not a finite number")
    return number


def _as_text(value):
    return "" if value is None else str(value).strip()


def _as_names(value):
    if value is None or value == "":
        return []
    words = value.replace(",", " ").split() if isinstance(value, str) else list(value)
    return [str(w) for w in words]


# The value types of a job's settings: type -> (conversion of a request's value, what the type allows).  The page has
# one widget per type (tabs/fit.js, SETTING_WIDGETS); a new type is a row here and a widget there.
SETTING_TYPES = {
    "integer": (_as_integer, "an integer, or empty for the default"),
    "boolean": (_as_boolean, "true or false"),
    "names": (_as_names, "names of free parameters, separated by blanks or commas; empty: every free parameter"),
    "number": (_as_number, "a finite number, or empty for the default"),
    "text": (_as_text, "a word or a path, as the setting's doc says"),
}


@dataclass(frozen=True)
class Setting:
    """One setting of a job kind: its name, type (a key of :data:`SETTING_TYPES`), default and doc."""
    name: str
    type: str
    default: object
    doc: str


@dataclass(frozen=True)
class JobRow:
    """One kind of job: ``kind`` and ``method`` name it; ``run`` its function ``(tracker, x0, settings, names) ->
    (state, message, extra)``, or ``missing`` the pyazr call it waits for; ``needs`` the package it imports;
    ``settings`` its :class:`Setting` rows; ``follows`` the (kind, method) of the job whose best point it starts from
    (the page offers it beside that job's result; method None: any job of that kind), or a tuple of such pairs;
    ``takes`` the key of that job's result it needs (a covariance, a chain: the job is refused without one, and offered
    only beside a result that has it); ``facts`` what the page shows of it beside its doc (name -> value: the program's
    switch, whether it gives a covariance, whether MINOS can follow); ``files`` the names of its settings that name a
    file (the server refuses one outside its root before the job starts); ``session`` False for a job that needs no
    session of the project (its worker builds none; x0 and the free names are empty); ``tab`` the page's tab whose
    form starts it (the Fit tab's job form offers the rows of "fit" only; every row is in the expert's job table).
    """
    kind: str
    method: str
    doc: str
    run: object = None
    needs: str = None
    settings: tuple = ()
    missing: str = None
    follows: tuple = None
    takes: str = None
    facts: tuple = ()  # (name, value) pairs
    files: tuple = ()
    session: bool = True
    tab: str = "fit"

    def defaults(self):
        return {s.name: s.default for s in self.settings}


MAX_CALLS_ENGINE = Setting("max_calls", "integer", None, "the most objective calls of the fit (the program's "
                           "--max-calls; empty: no limit)")

# The rows written here: the recipes on the objective, and the jobs that wait for a pyazr call.  The engine's fitters
# are rows made from its table (engine_rows); job_rows() is the whole table.
JOBS = (
    JobRow("fit", "iminuit_migrad", "iminuit's MIGRAD with the gradient, then HESSE (errors, covariance)",
           run=fit_iminuit_migrad, needs="iminuit",
           settings=(Setting("max_calls", "integer", None, "the most objective calls of MIGRAD (empty: iminuit's "
                             "default)"),
                     Setting("hesse", "boolean", True, "HESSE after a valid minimum: errors and covariance")),
           facts=(("covariance", True), ("MINOS can follow", True))),
    JobRow("fit", "scipy_lbfgsb", "scipy's L-BFGS-B on the objective and its gradient (no errors)",
           run=fit_scipy_lbfgsb, needs="scipy",
           settings=(Setting("max_iterations", "integer", 1000, "the most iterations of L-BFGS-B"),),
           facts=(("covariance", False), ("MINOS can follow", False))),
    JobRow("minos", "iminuit", "iminuit's MINOS intervals from the best point of an iminuit fit",
           run=minos_iminuit, needs="iminuit", follows=("fit", "iminuit_migrad"),
           settings=(Setting("parameters", "names", [], "the parameters MINOS scans (empty: every free one)"),
                     Setting("max_calls", "integer", None, "the most objective calls of each MIGRAD and MINOS scan "
                             "(empty: iminuit's default)"))),
    JobRow("sample", "emcee", "emcee's ensemble sampler on the log-probability (-0.5 objective, zero outside the fit's "
           "limits): walkers around the start point, or a sampler job's chain continued; one process",
           run=sample_emcee, needs="emcee", follows=(("fit", None), ("sample", "emcee")),
           settings=(Setting("walkers", "integer", None, "the walkers (empty: twice the free parameters plus two)"),
                     Setting("steps", "integer", 1000, "the stored steps of the chain"),
                     Setting("spread", "number", 1e-3, "the start walkers' spread around the start point, relative to "
                             "each value (absolute where it is 0)"),
                     Setting("thin", "integer", 1, "the sampler's steps between two stored steps"),
                     Setting("seed", "integer", 1, "the seed of the start walkers and of the sampler's draws"),
                     Setting("discard", "integer", None, "the stored steps left out of the summary as burn-in "
                             "(empty: the first half)")),
           facts=(("covariance", True), ("draws", "the chain in output/chain_<id>.npy"))),
    JobRow("bands", "covariance", "the 1-sigma band of the model of every data set and extrapolation at a fit's best "
           "point, from its covariance (linear, pyazr.bands)", run=bands_covariance, follows=("fit", None),
           takes="covariance", facts=(("from", "a fit job with a covariance (its best point as the start)"),)),
    JobRow("bands", "draws", "the band of the model of every data set and extrapolation from draws of a sampler's "
           "chain (16th, 50th, 84th percentiles at each point)", run=bands_draws, follows=("sample", "emcee"),
           takes="chain", settings=(Setting("draws", "integer", 200, "the draws of the chain (after its burn-in) "
                                            "the model is evaluated at"),
                                    Setting("seed", "integer", 1, "the seed of the choice of draws")),
           facts=(("from", "a sampler job (its best point as the start)"),)),
    JobRow("rate", "engine", "the reaction rate N_A<sigma v> of a pair of pairs at the start point, computed by the "
           "engine as the program's rate command (pyazr.rate), with the Gamow window, an optional band and an "
           "optional comparison table", run=rate_engine, tab="rate",
           follows=(("fit", None), ("sample", "emcee")),
           settings=(Setting("entrance", "integer", 1, "the entrance pair's key"),
                     Setting("exit", "integer", 2, "the exit pair's key; -1 the total radiative capture"),
                     Setting("grid", "text", "table", "table or log20 (pyazr.rate.TEMPERATURE_GRIDS), or custom: "
                             "t_min to t_max in steps of t_step (GK, the program's grid)"),
                     Setting("t_min", "number", None, "custom grid: the lowest temperature (GK)"),
                     Setting("t_max", "number", None, "custom grid: the highest temperature (GK)"),
                     Setting("t_step", "number", None, "custom grid: the step (GK)"),
                     Setting("band", "text", "none", "none; covariance (the start job's, linear in ln rate); draws "
                             "(the start job's chain: 16th, 50th, 84th percentiles)"),
                     Setting("draws", "integer", 100, "band from draws: the vectors of the chain the rate is "
                             "computed at"),
                     Setting("seed", "integer", 1, "band from draws: the seed of the choice of vectors"),
                     Setting("comparison", "text", "", "a rate table of the project folder to compare with: T9 (GK), "
                             "N_A<sigma v>, optionally its uncertainty"),
                     Setting("starlib_label", "text", "azr2", "STARLIB: the reference label (one to four characters)"),
                     Setting("starlib_flag", "text", "", "STARLIB: the interaction label: empty, w, g or v"),
                     Setting("reaclib_label", "text", "azr2", "REACLIB: the set label (four characters)"),
                     Setting("reaclib_flag", "text", "n", "REACLIB: the resonance flag: n, r or w"),
                     Setting("reaclib_tolerance", "number", 0.05, "REACLIB: the largest relative deviation of the "
                             "sets from the computed rate on 0.01-10 GK"),
                     Setting("reaclib_sets", "integer", 3, "REACLIB: the most sets tried")),
           facts=(("program command", "AZURE2 rate"), ("table", "output/rate_<id>.csv"))),
)


def engine_rows(fitters):
    """The job rows of the engine's fitter table (``pyazr.kinds()["fitters"]``: name -> row): a fit job per fitter
    (the program's default first, then the table's order), and a MINOS job after each fitter whose row says MINOS can
    follow it.  A new fitter of the engine is a new row here with no line of the UI changed."""
    order = sorted(fitters, key=lambda name: fitters[name]["option"] != "")
    fits, minos = [], []
    for name in order:
        f = fitters[name]
        method = f"engine_{name}"
        settings = (MAX_CALLS_ENGINE,)
        if f["minos"]:
            settings += (Setting("minos", "boolean", False, "MINOS after the fit and its HESSE: asymmetric "
                                 "intervals of every free parameter (Up = the session's chiVariance)"),)
        facts = (("program switch", f["option"] or "none (the program's default)"), ("covariance", f["covariance"]),
                 ("MINOS can follow", f["minos"]), ("gradient", f["gradient"]), ("residual rows", f["residual_rows"]),
                 ("honours a call limit", f["call_limit"]), ("falls back to", f["fallback"] or "none"))
        fits.append(JobRow("fit", method, f"the engine's {name}: {f['doc']}", settings=settings, facts=facts,
                           run=partial(fit_engine, method=name)))
        if f["minos"]:
            minos.append(JobRow("minos", method, f"the engine's MINOS after its {name} (the fit from the start point, "
                                "HESSE, MINOS on every free parameter)", settings=(MAX_CALLS_ENGINE,),
                                follows=("fit", method), run=partial(fit_engine, method=name, minos=True)))
    return tuple(fits), tuple(minos)


@lru_cache(maxsize=1)
def job_rows():
    """The whole job table: the engine's fit rows, the rows of :data:`JOBS`, the ENDF rows
    (:func:`pyazr.ui.evaluations.job_rows`), the engine's MINOS rows."""
    import pyazr
    from .evaluations import job_rows as endf_rows
    fits, minos = engine_rows(pyazr.kinds()["fitters"])
    return fits + JOBS + endf_rows() + minos


def job_row(kind, method):
    rows = job_rows()
    for row in rows:
        if row.kind == kind and row.method == method:
            return row
    known = [f"{r.kind}/{r.method}" for r in rows]
    raise NotFound(f"no job {kind}/{method} (the job table: {', '.join(known)})", known=known)


# The kinds of job the Fit tab's first view offers (the PI, 2026-11-14: fit, stop, accept, undo and the result against
# the data first); the others (MINOS, the sampler, the bands, ...) are one click away.  A row's ``first_view``.
FIRST_VIEW_KINDS = ("fit",)


def job_table():
    """The job table as data (GET /api/jobs/kinds)."""
    out = []
    for r in job_rows():
        available = r.missing is None and _importable(r.needs)
        settings = [{"name": x.name, "type": x.type, "default": x.default, "doc": x.doc,
                     "allowed": SETTING_TYPES[x.type][1]} for x in r.settings]
        out.append({"kind": r.kind, "method": r.method, "doc": r.doc, "settings": settings, "needs": r.needs,
                    "available": available, "missing": r.missing, "follows": _follows(r.follows),
                    "takes": r.takes, "facts": dict(r.facts), "files": list(r.files), "session": r.session,
                    "tab": r.tab, "first_view": r.kind in FIRST_VIEW_KINDS})
    return out


def _follows(follows):
    """A row's ``follows`` as data: [kind, method], or a list of such pairs (None where it follows no job)."""
    if not follows:
        return None
    return [list(p) for p in follows] if isinstance(follows[0], tuple) else list(follows)


def merged_settings(row, given):
    """The row's defaults with the request's values over them, each converted by its type; Refused on an unknown
    setting or a value its type does not allow."""
    by_name = {x.name: x for x in row.settings}
    unknown = sorted(set(given or {}) - set(by_name))
    if unknown:
        raise Refused(f"job {row.kind}/{row.method} has no setting {', '.join(unknown)}", allowed=sorted(by_name))
    out = row.defaults()
    for name, value in (given or {}).items():
        convert, allowed = SETTING_TYPES[by_name[name].type]
        try:
            out[name] = convert(value)
        except (TypeError, ValueError) as e:
            raise Refused(f"job {row.kind}/{row.method}, setting {name}: {e}", allowed=allowed) from None
    return out


def _importable(name):
    if name is None:
        return True
    import importlib.util
    return importlib.util.find_spec(name) is not None


# -- worker side ------------------------------------------------------------------------------------------------------

def write_checkpoint(path, record):
    """Write ``record`` as JSON to ``path`` atomically (a reader never sees half a file)."""
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(record, f)
    os.replace(tmp, path)


def plain_json(value):
    """``value`` with every non-finite float as None (JSON has no NaN or infinity; a browser refuses them)."""
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, dict):
        return {k: plain_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain_json(v) for v in value]
    return value


def listen_for_stop(channel, stop):
    """Set ``stop`` when the server asks the job to stop, or closes its input (a thread of the job's worker: the
    request is seen while the engine computes, which releases the GIL)."""
    def listen():
        while True:
            message = channel.receive()
            if message is None or message.get("op") == "stop":
                stop.set()
                return

    threading.Thread(target=listen, name="pyazr-ui-job-stop", daemon=True).start()


def run_job(channel):
    """A job worker: read the job's spec, build the session (unless the row needs none), run the row's function,
    report, write the checkpoint."""
    from .engine import objective_of, open_for_job
    spec = channel.receive()
    record = {"id": spec["id"], "kind": spec["kind"], "method": spec["method"], "text_rev": spec["rev"],
              "settings": spec["settings"], "state": "building", "progress": None, "best": None, "result": None}
    checkpoint = spec["checkpoint"]
    stop = threading.Event()
    listen_for_stop(channel, stop)

    def report(message):
        message = plain_json(message)
        channel.send(message)
        record["progress"] = {k: v for k, v in message.items() if k not in ("event", "best_x")}
        if message["best_x"] is not None:
            record["best"] = {"objective": message["best"], "x": message["best_x"]}
            write_checkpoint(checkpoint, record)

    try:
        write_checkpoint(checkpoint, record)
        row = job_row(spec["kind"], spec["method"])
        t0 = time.perf_counter()
        session, x0, names, fg = None, [], [], None
        if row.session:
            session, x0, names = open_for_job(spec["session_file"], spec["cwd"])
            fg = objective_of(session).value_and_grad
        channel.send({"event": "built", "build_s": time.perf_counter() - t0, "names": names})
        if spec.get("x0") and row.session:
            x0 = vector_from_names(spec["x0"], names, x0)
        record["state"] = "running"
        write_checkpoint(checkpoint, record)
        tracker = Tracker(fg, names, report=report, stop=stop, session=session, source=spec.get("source"),
                          job_id=spec["id"], products=spec["products"])
        try:
            status, message, extra = row.run(tracker, x0, spec["settings"], names)
        except StopRequested:
            status, message, extra = "stopped", "stopped on request", {}
        best_x = None if tracker.best_x is None else dict(zip(names, tracker.best_x.tolist()))
        record.update(state=status, best={"objective": tracker.best, "x": best_x} if best_x else None,
                      result={"status": status, "message": message, "calls": tracker.calls,
                              "start_objective": tracker.first, "start_x": dict(zip(names, np.asarray(x0, float)
                                                                                    .tolist())),
                              "elapsed": time.perf_counter() - tracker.t0, **extra})
        if stop.is_set():
            record["result"]["reports_after_stop"] = tracker.after_stop
    except Exception as e:  # noqa: BLE001  -- reported as the job's failure
        classes, text, _ = describe_exception(e)
        cls = classes[0]
        record.update(state="failed", result={"status": "failed", "error": cls, "message": text})
    record = plain_json(record)
    write_checkpoint(checkpoint, record)
    channel.send({"event": "done", "record": record})


def vector_from_names(values, names, x0):
    """x with the entries named in ``values`` ({name: value}); every name must be a free parameter of ``names``."""
    unknown = sorted(set(values) - set(names))
    if unknown:
        raise KeyError(f"not free parameters of this session: {', '.join(unknown[:8])} (free: {', '.join(names)})")
    x = np.array(x0, float)
    for i, n in enumerate(names):
        if n in values:
            x[i] = float(values[n])
    return x


# -- server side ------------------------------------------------------------------------------------------------------

class Job:
    """A job as the server sees it: its record (the checkpoint's fields) and its worker."""

    def __init__(self, record, process=None, log_path=None):
        self.record, self.process, self.session_file, self.log_path = record, process, None, log_path
        self.lock = threading.Lock()
        self.stop_asked_at = None
        self.trace_stride = 1  # a call enters the trace when its number is a multiple of this

    def add_to_trace(self, calls, objective, best):
        """Keep (call, objective, best) in the record's trace: every call until TRACE_POINTS, then every second one
        of those kept and every second call after (the trace stays a picture of the whole run).  Under the lock."""
        trace = self.record.setdefault("trace", [])
        if calls % self.trace_stride:
            return
        trace.append([calls, objective, best])
        if len(trace) > TRACE_POINTS:
            self.trace_stride *= 2
            trace[:] = [p for p in trace if p[0] % self.trace_stride == 0]

    def view(self):
        with self.lock:
            out = dict(self.record)
            out["pid"] = self.process.pid if self.process is not None else None
        out["log"] = log_tail(self.log_path)
        return out


def log_tail(path, lines=LOG_LINES):
    """The last ``lines`` non-empty lines of a job's log (the engine's notes; a carriage return ends a line too)."""
    if not path or not os.path.exists(path):
        return []
    with open(path, "rb") as f:
        f.seek(max(0, os.path.getsize(path) - 8192))
        text = f.read().decode("utf-8", "replace")
    return [line for line in text.replace("\r", "\n").split("\n") if line.strip()][-lines:]


class JobTable:
    """The jobs of the open project: started, followed by a reader thread each, stopped, killed, listed; those
    found in the jobs directory at start-up (a server that died) are ``interrupted``."""

    def __init__(self, directory, budget, spawner):
        self.directory, self.budget, self.spawner = directory, budget, spawner
        self.jobs = {}
        os.makedirs(directory, exist_ok=True)
        for name in sorted(os.listdir(directory)):
            if name.endswith(".json"):
                with open(os.path.join(directory, name)) as f:
                    record = json.load(f)
                if record.get("state") not in FINAL:
                    record["state"] = "interrupted"
                self.jobs[record["id"]] = Job(record, log_path=os.path.join(directory, f"{record['id']}.log"))

    def start(self, kind, method, settings, x0, rev, session_text, cwd, file_of, source=None, products=None):
        """Start a job of the table's row (kind, method) on ``session_text`` (the session file's text of revision
        ``rev``, written to ``file_of(job id)``) from ``x0`` ({name: value} or None for the text's x0); ``source``:
        the record of the job x0 was taken from (a row that ``takes`` a key of its result needs it); ``products``:
        the directory of the files the job writes (the project's output directory)."""
        from .workers import WorkerProcess
        row = job_row(kind, method)
        if row.missing:
            raise Missing(f"job {kind}/{method} needs a pyazr call that does not exist yet: {row.missing}",
                          missing=row.missing)
        if not _importable(row.needs):
            raise Missing(f"job {kind}/{method} needs the package {row.needs}, not installed", missing=row.needs)
        if row.takes and not ((source or {}).get("result") or {}).get(row.takes):
            what = f"job {source['id']} ({source['kind']} {source['method']}) has no {row.takes}" if source \
                else "no job is its start point"
            raise Refused(f"job {kind}/{method} starts from a job whose result has a {row.takes}: {what}",
                          allowed=f"x0: {{\"job\": id}} of a job with a {row.takes}",
                          fix="start it from such a job (the Fit tab offers it beside the job's result)")
        merged = merged_settings(row, settings)
        job_id = uuid.uuid4().hex[:8]
        threads = self.budget.take_job(job_id)
        try:
            session_file = file_of(job_id)
            with open(session_file, "wb") as f:
                f.write(session_text.encode("utf-8", "surrogateescape"))
            spec = {"id": job_id, "kind": kind, "method": method, "settings": merged, "x0": x0, "rev": rev,
                    "session_file": session_file, "cwd": cwd, "source": source, "products": products or cwd,
                    "checkpoint": os.path.join(self.directory, f"{job_id}.json")}
            process = WorkerProcess(self.spawner, "job", threads=threads,
                                    log_path=os.path.join(self.directory, f"{job_id}.log"), cwd=cwd)
            process.send(spec)
        except Exception:
            self.budget.give_back(job_id)
            raise
        job = Job({"id": job_id, "kind": kind, "method": method, "text_rev": rev, "settings": merged,
                   "state": "building", "progress": None, "best": None, "result": None,
                   "started": time.time(), "trace": []}, process, log_path=process.log_path)
        job.session_file = session_file
        self.jobs[job_id] = job
        threading.Thread(target=self._follow, args=(job,), name=f"pyazr-ui-job-{job_id}", daemon=True).start()
        return job.view()

    def _follow(self, job):
        """Read the job's messages until it ends; then end its worker, give its thread back to the budget, and only
        then show its final state (a client that sees the state final finds the budget free)."""
        final = None
        try:
            while True:
                message = job.process.receive()
                if message["event"] == "done":
                    final = dict(message["record"])
                    break
                with job.lock:
                    if message["event"] == "built":
                        job.record["state"] = "stopping" if job.stop_asked_at else "running"
                        job.record["build_s"] = message["build_s"]
                    elif message["event"] == "call":
                        job.record["progress"] = {k: v for k, v in message.items() if k not in ("event", "best_x")}
                        job.add_to_trace(message["calls"], message["objective"], message["best"])
                        if message["best_x"] is not None:
                            job.record["best"] = {"objective": message["best"], "x": message["best_x"]}
        except UIError:  # the worker ended without "done": killed, or crashed
            with job.lock:
                final = dict(job.record)
            final["state"] = "killed" if final.get("kill_asked") else "failed"
            if final["state"] == "failed":
                final["result"] = {"status": "failed", "message": "the worker ended without a result",
                                   "log": job.process.log_path}
            write_checkpoint(os.path.join(self.directory, f"{final['id']}.json"), final)
        finally:
            job.process.kill()
            self.budget.give_back(job.record["id"])
            if job.session_file and os.path.exists(job.session_file):
                os.remove(job.session_file)
            with job.lock:
                job.record.update(final or {"state": "failed"})
                job.record["ended"] = time.time()

    def wait(self, job_id, seconds):
        """The job once its state differs from the one it has now, or after ``seconds`` (at most WAIT_MAX): a client
        follows a running job with one request in flight and sees a stop, a kill or the end at once."""
        job = self.get(job_id)
        state, end = job.record["state"], time.monotonic() + min(max(float(seconds), 0.0), WAIT_MAX)
        if state in FINAL:
            return job
        while job.record["state"] == state and time.monotonic() < end:
            time.sleep(0.02)
        return job

    def get(self, job_id):
        if job_id not in self.jobs:
            raise NotFound(f"no job {job_id!r} (the jobs: {', '.join(self.jobs) or 'none'})", known=list(self.jobs))
        return self.jobs[job_id]

    def stop(self, job_id):
        """Cooperative stop: the job ends before its next objective call, with the best point seen."""
        job = self.get(job_id)
        with job.lock:
            if job.record["state"] in FINAL:
                return dict(job.record)
            job.stop_asked_at = time.time()
            job.record["state"] = "stopping"
        try:
            job.process.send({"op": "stop"})
        except UIError:
            pass
        return job.view()

    def kill(self, job_id):
        """SIGKILL the job's worker; its state becomes ``killed`` (the best point is the last reported)."""
        job = self.get(job_id)
        with job.lock:
            if job.record["state"] in FINAL:
                return dict(job.record)
            job.record["kill_asked"] = True
        if job.process is not None:
            job.process.proc.kill()
        return job.view()

    def running(self):
        return [j for j in self.jobs.values() if j.record["state"] not in FINAL]

    def best_x(self, job_id):
        """{name: value} of the job's best point (its checkpoint when the server restarted)."""
        job = self.get(job_id)
        with job.lock:
            best = job.record.get("best")
        if not best or not best.get("x"):
            raise Conflict(f"job {job_id} has no best point yet", fix=f"GET /api/jobs/{job_id}")
        return best["x"]

    def close(self):
        """Kill every running job's worker (server shutdown)."""
        for job in list(self.jobs.values()):
            if job.process is not None and job.process.alive():
                job.record["kill_asked"] = True
                job.process.kill()
