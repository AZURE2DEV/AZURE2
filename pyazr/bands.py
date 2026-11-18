"""Covariance-propagated uncertainty bands for AZURE2 observables.

A fit leaves behind two things: the best-fit parameters (``output/param.sav``)
and their covariance (``output/covariance.dat``).  The 1-sigma band on any
observable the model can compute follows from linear propagation,

.. math::  \\sigma_y^2 = g^{\\mathsf T} C\\, g, \\qquad g_j = \\partial y/\\partial p_j

which is the same expression AZURE2 uses internally to write its ``.band``
files (SAMMY Eq. IV E4.2).  ``covariance.dat`` spans exactly the **free
R-matrix parameters** -- level energies and reduced-width amplitudes, in their
fit-vector order -- because an observable is insensitive to the dataset
normalizations and energy shifts that make up the rest of the fit vector.

The sensitivities :math:`g` come from AZURE2 itself -- one reverse-mode adjoint
per point, so the whole matrix costs about two forward evaluations however many
parameters the model has.  Finite differences are available as a cross-check
(``method='central'``), at ``2 * ncols`` forward passes.

Two entry points:

``uncertainty_bands(model, segments, ...)``
    band for selected segments of a *live* :class:`~pyazr.Session`
    session, in whichever mode it is in (data or extrapolation).  The session
    is the in-process engine itself -- one ``_azure2.Session``, no subprocess
    or instance pool to talk to.

``extrapolation_bands(azr_file, keys=..., grids=..., ...)``
    band on an extrapolation grid, computed in a **dedicated session whose
    ``.azr`` carries only the requested segments**.  Each finite-difference
    adjoint pass covers every active segment of the model, so trimming the
    extrapolations to what you actually want to plot still pays:
    on the full 7Be model one pass over all 110 grids is ~31 s, against
    ~0.1 s for a handful of them.

Examples
--------
>>> from pyazr import Session, extrapolation_bands
>>> bands = extrapolation_bands("7Be.azr", keys=[113], quantity="sfactor")
>>> b = bands[113]
>>> b.value[0], b.sigma[0]                          # S(E->0) +- 1 sigma
"""

import os
import shutil
import sys
import tempfile
import warnings
from dataclasses import dataclass, field

import numpy as np

from .errors import MissingFile
from .locking import session_lock

__all__ = ["Band", "load_covariance", "rmatrix_columns", "best_fit_params",
           "live_parameters", "step_sizes", "sensitivities", "active_indices",
           "trimmed_model", "uncertainty_bands", "extrapolation_bands"]


# -- the band -----------------------------------------------------------------

@dataclass
class Band:
    """One segment's best-fit curve and its covariance uncertainty.

    ``value`` and ``sigma`` are in the units selected by ``quantity``: the
    cross section as AZURE2 reports it (b, or b/sr for a differential segment)
    or the S factor (MeV b).  ``energy``/``angle`` are AZURE2's own calculated
    center-of-mass grids for the segment.
    """

    segment: int                     # calculated-segment index (0-based)
    key: int                         # 1-based segment number in the .azr file
    energy: np.ndarray
    angle: np.ndarray
    excitation: np.ndarray
    value: np.ndarray
    sigma: np.ndarray
    quantity: str = "cross"
    sensitivity: np.ndarray = field(default=None, repr=False)

    def interval(self, nsigma=1.0):
        """``(lower, upper)`` at ``nsigma`` standard deviations."""
        d = nsigma * self.sigma
        return self.value - d, self.value + d

    @property
    def relative(self):
        """Fractional uncertainty, ``sigma / |value|`` (0 where the value is)."""
        out = np.zeros_like(self.sigma)
        nz = self.value != 0
        out[nz] = self.sigma[nz] / np.abs(self.value[nz])
        return out

    def __len__(self):
        return len(self.value)

    def __repr__(self):
        rel = self.relative
        return (f"Band(segment={self.segment}, key={self.key}, "
                f"{len(self.value)} pts, quantity={self.quantity!r}, "
                f"E={self.energy.min():g}-{self.energy.max():g}, "
                f"rel. unc. {100 * rel.min():.2f}-{100 * rel.max():.2f}%)")


# -- inputs -------------------------------------------------------------------

def _model_dir(model):
    """The directory a session's ``output/`` lives in.

    ``Session`` switches the process into the ``.azr``'s directory for the
    session's lifetime, and AZURE2 names ``output/`` relative to that -- so
    that is where ``param.sav`` and ``covariance.dat`` are.
    """
    return os.path.abspath(model.cwd or os.path.dirname(model.file) or ".")


def load_covariance(path="output/covariance.dat", ncols=None, names=None):
    """Read AZURE2's saved parameter covariance (a bare ``N x N`` matrix).

    ``covariance.dat`` is written by a fit run with the uncertainty band
    enabled and holds the free parameters a band spans (``pyazr.kinds()["parameters"]``' ``band``),
    in packed parameter order, after a ``#columns`` header that names them.
    ``ncols``, when given, is the number of such parameters the current model
    has; a mismatch means the covariance belongs to a different model and is
    refused rather than silently misaligned -- unless ``names`` (the current
    model's band columns) are given and the header names each of them: then
    the matrix is taken on those columns (a session of the extrapolations
    after a fit with a free target parameter, which that session has not).
    """
    C = np.loadtxt(path, ndmin=2)
    if C.ndim != 2 or C.shape[0] != C.shape[1]:
        raise ValueError(f"{path}: expected a square matrix, got {C.shape}.")
    if ncols is not None and C.shape[0] != ncols and names is not None and len(names) == ncols:
        with open(path) as f:
            header = next((line.split()[1:] for line in f if line.startswith("#columns")), [])
        if len(header) == C.shape[0] and all(n in header for n in names):
            pick = [header.index(n) for n in names]
            C = C[np.ix_(pick, pick)]
    if ncols is not None and C.shape[0] != ncols:
        raise ValueError(
            f"{path} is {C.shape[0]}x{C.shape[0]} but the model has {ncols} "
            f"free R-matrix parameters.  The covariance was saved for a "
            f"different model (or a different set of fixed parameters); re-run "
            f"the fit with the uncertainty band enabled.")
    return C


_RMATRIX_KINDS = ("energy", "width")


def _band_kinds():
    """The kinds of fit parameter a band spans (the engine's ``IsBandParameter``, ``pyazr.kinds()["parameters"]``'s
    ``band``): the R-matrix parameters, the THM coherent background (``cbkg``) and, since 2026-11-18, the parameters of
    targets, samples and apparatus that change the model (abundance, screening, density, doppler, thickness,
    apparatus).  Normalizations and energy shifts move the data, not the model, and are not in it."""
    return tuple(k["key"] for k in _engine_module().parameter_kinds() if k["band"])


def _engine_module():
    from . import _azure2
    return _azure2


def live_parameters(model):
    """``(free_kinds, free_identities)`` queried from the running session.

    ``model.parameters`` is built once, from whatever mode the session was in
    at the time -- but the parameter vector is **not** the same in the two
    modes.  Extrapolation segments carry no normalizations, and AZURE2 drops
    the fixed R-matrix parameters its grids cannot reach, so a model with 483
    parameters (197 free) in a session of the data can report 234 (72 free) in
    a session of the extrapolations.  The free R-matrix block survives, in order, and it is
    that block the covariance describes.  This reads the engine's *current*
    parameter records (``Session.parameters``), rebuilt on every
    ``initialize()``, so it is right for the mode the session is in now.

    Returns the ``kind`` of each free parameter (keys of ``pyazr.kinds()["parameters"]``) and its identity
    ``(kind, jgroup, level, channel)`` (``None`` where a field does not apply) -- the key AZURE2 uses to match
    covariance columns across runs.
    """
    with session_lock(model):
        free = [r for r in model._engine.parameter_records() if not r["fixed"]]
    return ([r["kind"] for r in free],
            [(r["kind"], r.get("jgroup"), r.get("level"), r.get("channel")) for r in free])


def rmatrix_columns(model):
    """Packed indices of the free parameters a band spans (``pyazr.kinds()["parameters"]``' ``band``), in the order of the fit vector.

    These are the columns ``covariance.dat`` spans: level energies and reduced
    width amplitudes, the THM coherent background (``cbackground=``, kind
    ``cbkg``) when the project has one, and the free target and apparatus
    parameters (a fitted areal density, a stopping scale, ...; since 2026-11-18).
    Normalizations and energy shifts are excluded -- they move the data, not the
    model.  Queried live, so the answer is right for the mode the session is in now.
    """
    kinds, _ = live_parameters(model)
    band = _band_kinds()
    cols = [j for j, k in enumerate(kinds) if k in band]
    if not cols:
        raise RuntimeError("the model has no free R-matrix parameters.")
    return cols


def best_fit_params(model, path=None):
    """The fit vector of the parameter file ``path`` (default ``output/param.sav`` beside the model) for the session
    ``model``: the file's free values, by name, in the order of ``model.parameters``.

    A model's ``.azr`` file need not carry the fitted values -- ``param.sav``
    is what a fit writes -- so a band must be evaluated at the parameters the
    covariance belongs to.  On the 7Be model the difference is not subtle: the
    ``.azr``'s own parameters give chi2 = 173006, ``param.sav`` gives 4130.

    The file holds the fixed level parameters too, and a fit vector cannot carry them: the session must hold the
    file's fixed values.  A session made with ``Session(file, parameter_file=path)`` does, and its ``parameters.x0``
    is returned.  Otherwise the file is read by a second session of the same project and options made with
    ``parameter_file=path`` (``model`` is not changed); its vector is returned when its fixed parameters equal
    ``model``'s, and :class:`pyazr.ParameterError` is raised when the file sets a fixed parameter to another value
    (then make the session with ``parameter_file=path``).
    """
    from .errors import ParameterError
    from .session import Session

    if path is None:
        path = os.path.join(_model_dir(model), "output", "param.sav")
    path = os.path.abspath(os.fspath(path))
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found: the band must be evaluated at the fitted "
            f"parameters the covariance belongs to.  Pass params= explicitly "
            f"to override.")
    if model.options.get("parameter_file") == path:
        return model.parameters.x0
    with Session(model.file, cwd=model.cwd, **dict(model.options, parameter_file=path)) as other:
        held = {p.name: p.start for p in model.parameters if p.fixed}
        given = {p.name: p.start for p in other.parameters if p.fixed}
        moved = sorted(n for n in held if given.get(n) != held[n])
        if moved:
            raise ParameterError(
                f"best_fit_params: {path} sets {len(moved)} fixed parameter(s) to other values than the session "
                f"holds ({', '.join(moved[:5])}{', ...' if len(moved) > 5 else ''}); a fit vector cannot carry "
                f"fixed values: make the session with Session(file, parameter_file={path!r}) and use its "
                f"parameters.x0.")
        return other.parameters.x0


# -- sensitivities ------------------------------------------------------------

def _fetch(sess, segments, quantity):
    """Observable values for ``segments`` from the last UPDATE_SEGMENTS call."""
    out = []
    for i in segments:
        v = sess.calculated_segments(i)
        if quantity == "sfactor":
            v = v * sess.calculated_conv(i)
        out.append(v)
    return out


def _evaluate(model, params, segments, quantity):
    with session_lock(model):  # the forward pass and the arrays read after it, one after the other
        s = model._engine
        n = int(s.update_segments_rwa(np.asarray(params, float)))
        bad = [i for i in segments if i < 0 or i >= n]
        if bad:
            raise IndexError(f"segment index/indices {bad} out of range: the model "
                             f"calculates {n} (active) segments.")
        return _fetch(s, segments, quantity)


# -- finite differences across processes ---------------------------------------
#
# The engine is one non-reentrant C++ object per process, so threads cannot
# spread the columns the way the socket API's instance pool did -- each worker
# needs an engine of its own, which means a process of its own.  The pieces
# below are module-level because a "spawn" pool has to import them by name.

_WORKER = {}


def _worker_setup(spec):
    """Build this worker's own session, with an output directory of its own.

    Workers must not share ``output/``: AZURE2 writes its external-capture
    cache there, and several processes deciding at once that the cache is
    stale would write over each other.  Each gets a private directory seeded
    from the parent's caches, so the integrals are reused rather than rebuilt
    -- the same trick :func:`trimmed_model` uses.
    """
    from .session import Session
    from .project import Project

    cwd, workdir = spec["cwd"], os.path.join(spec["workdir"], str(os.getpid()))
    os.makedirs(workdir, exist_ok=True)
    # Best effort: a missing seed only means this worker rebuilds its own
    # integrals, which is correct, just slower.
    # The .sig sidecar goes with its file: without it AZURE2 cannot tell the
    # copy is valid and recomputes it.
    for cache in ("intEC.dat", "intEC.dat.sig", "intEC.extrap", "intEC.extrap.sig"):
        seed = os.path.join(spec["outputs"], cache)
        if os.path.exists(seed):
            shutil.copyfile(seed, os.path.join(workdir, cache))

    # The edited .azr lives in the private directory too, so the parent removing
    # that directory takes the temporary file with it.
    # the engine reads a format-3 project's data files relative to cwd
    edited = Project.open(spec["file"])
    edited.set_config(output_dir=workdir.rstrip("/") + "/")
    fd, path = tempfile.mkstemp(suffix=".azr", dir=workdir)
    os.close(fd)
    edited.save(path)

    model = Session(path, cwd=cwd, **spec["options"])  # the options carry extrapolations=
    _WORKER.update(model=model, base=None)


def _worker_column(task):
    """One finite-difference column, evaluated in this worker's own engine."""
    k, j, h_k, p0, segments, quantity, method = task
    model = _WORKER["model"]
    p0 = np.asarray(p0, float)

    pp = p0.copy(); pp[j] += h_k
    plus = _evaluate(model, pp, segments, quantity)
    if method == "central":
        pm = p0.copy(); pm[j] -= h_k
        minus = _evaluate(model, pm, segments, quantity)
        return k, [(a - b) / (2.0 * h_k) for a, b in zip(plus, minus)]

    if _WORKER["base"] is None:                    # forward: one shared baseline
        _WORKER["base"] = _evaluate(model, p0, segments, quantity)
    return k, [(a - b) / h_k for a, b in zip(plus, _WORKER["base"])]


def _parallel_columns(model, cols, h, p0, segments, quantity, method,
                      nprocs, progress):
    """Finite-difference columns spread over ``nprocs`` worker processes."""
    import multiprocessing

    spec = {
        "file": model.file,
        "cwd": model.cwd,
        "options": dict(getattr(model, "options", {})),
        "outputs": model.output_dir,
        # each worker's output directory under a temporary one, removed at the end: nothing beside the project
        "workdir": tempfile.mkdtemp(prefix="pyazr_sensitivities_"),
    }
    tasks = [(k, cols[k], float(h[k]), p0, segments, quantity, method)
             for k in range(len(cols))]

    columns = [None] * len(cols)
    # "spawn", not "fork": forking a live engine would duplicate a C++ object
    # graph and its OpenMP runtime into a child that never asked for it.
    ctx = multiprocessing.get_context("spawn")
    try:
        with ctx.Pool(nprocs, initializer=_worker_setup, initargs=(spec,)) as pool:
            for done, (k, column) in enumerate(
                    pool.imap_unordered(_worker_column, tasks), start=1):
                columns[k] = column
                if progress:
                    print(f"\r  sensitivities: {done}/{len(cols)} columns "
                          f"({nprocs} processes)", end="", flush=True)
        if progress:
            print()
    finally:
        shutil.rmtree(spec["workdir"], ignore_errors=True)
    return columns


def step_sizes(p0, cols, step=1e-4):
    """Finite-difference steps: ``step`` times each parameter's own magnitude.

    The step has to stay inside the model's linear regime.  Scaling it by the
    parameter's *uncertainty* is tempting -- one number then means the same
    thing for every column -- but a fit typically has a few nearly
    unconstrained directions (this 7Be fit has reduced widths with sigma ~ 380),
    and a step of even a few percent of such a sigma throws the model far from
    the minimum: the difference quotient then measures a chord across a
    resonance, not the derivative at the fit.  Scaling by the value keeps every
    perturbation small in the sense that matters.

    A parameter at (or near) zero gets the median magnitude as its scale, so
    its column is neither zero nor absurdly large.
    """
    p = np.abs(np.asarray(p0, float)[cols])
    typical = float(np.median(p[p > 0])) if np.any(p > 0) else 1.0
    return step * np.maximum(p, typical)


def sensitivities(model, segments, params=None, quantity="cross", step=1e-4,
                  method="analytic", cols=None, progress=False, nprocs=1):
    """``d(observable)/d(parameter)`` for selected segments.

    Parameters
    ----------
    model : a live :class:`~pyazr.Session` session.
    segments : calculated-segment indices (0-based, counting **active**
        segments only -- see :func:`active_indices`).
    params : the free parameter vector to differentiate about; defaults to
        ``model.parameters.x0``.
    quantity : ``'cross'`` or ``'sfactor'``.
    step : relative perturbation, or a per-column array of absolute steps, for
        the finite-difference methods only.  See :func:`step_sizes`.
    method : ``'analytic'`` (the default) asks AZURE2 for the exact
        sensitivities -- one reverse-mode adjoint per point, about two forward
        evaluations for the whole matrix.  ``'central'`` (2 evaluations per
        column) and ``'forward'`` (1 per column) finite-difference instead, and
        cost ``2 * ncols`` / ``ncols`` forward passes: on the 7Be capture grid
        that is 40 s against 0.8 s.  Use them only to cross-check, or with an
        AZURE2 that predates the analytic command.
    cols : packed parameter indices to differentiate with respect to;
        defaults to :func:`rmatrix_columns`.  Ignored by ``'analytic'``, which
        always returns the free R-matrix columns.
    progress : print a one-line counter as columns are evaluated.
    nprocs : worker processes to spread the finite-difference columns over
        (ignored by ``'analytic'``, which is one call).  The engine is one
        non-reentrant object per process, so each worker builds its own
        session: a fixed cost of one model initialization -- and, in
        a session of the extrapolations, one grid rebuild -- before it evaluates anything.
        The columns have to be worth more than that.  As a rule of thumb it
        needs ``len(cols)`` in the hundreds: on the 14-parameter ``tests/13N``
        model it is a net loss at any ``nprocs``, in a session of the data (0.8 s serial)
        and on an 1845-point extrapolation grid (53 s serial) alike.  Measure
        on your own model rather than assuming.

    Returns
    -------
    dict
        ``{segment: array of shape (npoints, len(cols))}``.
    """
    p0 = np.asarray(model.parameters.x0 if params is None else params, float)

    if method == "analytic":
        G = model.model_gradients(p0)
        bad = [i for i in segments if i < 0 or i >= len(G)]
        if bad:
            raise IndexError(f"segment index/indices {bad} out of range: the "
                             f"model calculates {len(G)} (active) segments.")
        if quantity == "sfactor":
            # S = conv * sigma with conv fixed by the kinematics, so the
            # sensitivities scale by the same factor, row by row.
            with session_lock(model):
                _evaluate(model, p0, segments, quantity)
                return {seg: G[seg] * model._engine.calculated_conv(seg)[:, None]
                        for seg in segments}
        return {seg: G[seg] for seg in segments}

    cols = list(rmatrix_columns(model) if cols is None else cols)
    h = np.atleast_1d(np.asarray(step, float))
    h = step_sizes(p0, cols, float(h[0])) if h.size == 1 else h
    if h.size != len(cols):
        raise ValueError(f"step has {h.size} entries, expected {len(cols)}.")

    if nprocs > 1 and len(cols) > 1:
        columns = _parallel_columns(model, cols, h, p0, segments, quantity,
                                    method, min(nprocs, len(cols)), progress)
    else:
        base = None if method == "central" else _evaluate(
            model, p0, segments, quantity)

        columns = []
        for k in range(len(cols)):
            j = cols[k]
            pp = p0.copy(); pp[j] += h[k]
            plus = _evaluate(model, pp, segments, quantity)
            if method == "central":
                pm = p0.copy(); pm[j] -= h[k]
                minus = _evaluate(model, pm, segments, quantity)
                columns.append([(a - b) / (2.0 * h[k]) for a, b in zip(plus, minus)])
            else:
                columns.append([(a - b) / h[k] for a, b in zip(plus, base)])
            if progress:
                print(f"\r  sensitivities: {k + 1}/{len(cols)} columns",
                      end="", flush=True)
        if progress:
            print()

    return {seg: np.column_stack([columns[k][s] for k in range(len(cols))])
            for s, seg in enumerate(segments)}


def active_indices(segments):
    """Map file order to calculated order for the session's data sets or extrapolations
    (:class:`pyazr.DescribedDatasets`, :class:`pyazr.DescribedExtrapolations`).

    AZURE2 only calculates the segments flagged active, and the API indexes
    *those*, so a model with inactive segments has calculated index != file
    index.  Returns ``{1-based key: calculated index}``.
    """
    out, k = {}, 0
    for s in segments:
        if getattr(s, "active", True):
            out[s.key] = k
            k += 1
    return out


# -- bands --------------------------------------------------------------------

def uncertainty_bands(model, segments=None, keys=None, params=None,
                      covariance=None, quantity="cross",
                      step=1e-4, method="analytic", progress=False, nprocs=1):
    """Covariance uncertainty bands for selected segments of a live session.

    Parameters
    ----------
    model : a live :class:`~pyazr.Session`, in data or extrapolation
        mode; the bands are computed for whichever segments that mode exposes.
    segments : calculated-segment indices (0-based).  Mutually exclusive with
        ``keys``.
    keys : 1-based segment numbers as they appear in the ``.azr`` file, which
        is what :attr:`~pyazr.Session.extrapolations` /
        :attr:`~pyazr.Session.datasets` report.  Translated to calculated
        indices for you (inactive segments are not calculated).
    params : parameters to evaluate at; defaults to ``output/param.sav``, since
        that -- not the ``.azr`` -- is what a fit writes and what the
        covariance describes.  Pass ``model.parameters.x0`` to use the file's own
        values.
    covariance : path to ``covariance.dat``, or the matrix itself.  Defaults to
        ``output/covariance.dat`` beside the model.
    quantity : ``'cross'`` or ``'sfactor'``.
    step : finite-difference step, relative to each parameter's own magnitude
        (see :func:`step_sizes`); ignored by the analytic method.
    method : ``'analytic'``, ``'central'`` or ``'forward'`` -- see
        :func:`sensitivities`.  ``'analytic'`` falls back to ``'central'``,
        with a warning, on an AZURE2 that lacks the sensitivity command.
    progress : print a progress counter (finite differences only).
    nprocs : worker processes for the finite-difference columns -- see
        :func:`sensitivities`.  Worth setting when the analytic method is
        unavailable and the model has many free R-matrix parameters.

    Returns
    -------
    dict
        ``{1-based segment key: Band}``.
    """
    kinds, _ = live_parameters(model)
    band = _band_kinds()
    cols = [j for j, k in enumerate(kinds) if k in band]
    if not cols:
        raise RuntimeError("the model has no free R-matrix parameters.")

    C = covariance
    if C is None:
        C = os.path.join(_model_dir(model), "output", "covariance.dat")
    if isinstance(C, (str, os.PathLike)):
        free_names = [p.name for p in model.parameters if not p.fixed]
        C = load_covariance(C, ncols=len(cols),
                            names=[free_names[j] for j in cols] if len(free_names) == len(kinds) else None)
    C = np.asarray(C, float)
    if C.shape[0] != len(cols):
        raise ValueError(f"the covariance is {C.shape[0]}x{C.shape[0]} but the "
                         f"model has {len(cols)} free R-matrix parameters.")

    if params is None:
        params = best_fit_params(model)
    p0 = np.asarray(params, float)
    if p0.size != len(kinds):
        raise ValueError(
            f"params has {p0.size} entries but the running instance expects "
            f"{len(kinds)} free parameters.  In a session of the extrapolations the vector "
            f"is the R-matrix block alone -- see pyazr.bands.live_parameters.")

    catalog = (model.extrapolations if getattr(model, "mode", "data") == "extrap"
               else model.datasets)
    calc_index = active_indices(catalog)
    key_of = {v: k for k, v in calc_index.items()}

    if (segments is None) == (keys is None):
        raise ValueError("give exactly one of segments= or keys=.")
    if keys is not None:
        missing = [k for k in keys if k not in calc_index]
        if missing:
            raise KeyError(f"segment key(s) {missing} are not active in "
                           f"{model.file} (or do not exist).")
        segments = [calc_index[k] for k in keys]
    segments = [int(s) for s in segments]

    try:
        G = sensitivities(model, segments, params=p0, quantity=quantity,
                          step=step, method=method, cols=cols,
                          progress=progress, nprocs=nprocs)
    except RuntimeError as err:
        if method != "analytic" or "MODEL_GRADIENTS" not in str(err):
            raise
        # Loud, and on stderr as well as through warnings: the fallback is
        # ~len(cols) times slower, so silently taking it looks like a hang.
        note = (f"{err}\n  Falling back to central differences: "
                f"{2 * len(cols)} forward passes instead of ~2.  Install an "
                f"AZURE2 built with the sensitivity command to avoid this.")
        print(f"pyazr.bands: {note}", file=sys.stderr, flush=True)
        warnings.warn(note, RuntimeWarning, stacklevel=2)
        G = sensitivities(model, segments, params=p0, quantity=quantity,
                          step=step, method="central", cols=cols,
                          progress=progress, nprocs=nprocs)

    # Leave the session holding the best fit again, and take the central
    # values from that same pass.
    with session_lock(model):  # the forward pass and the arrays read after it, one after the other
        values = _evaluate(model, p0, segments, quantity)
        s = model._engine

        out = {}
        for s_, seg in enumerate(segments):
            g = G[seg]
            if g.shape[1] != C.shape[0]:
                raise RuntimeError(
                    f"segment {seg}: {g.shape[1]} sensitivity columns against a "
                    f"{C.shape[0]}x{C.shape[0]} covariance.")
            if g.shape[0] != len(values[s_]):
                raise RuntimeError(
                    f"segment {seg}: {g.shape[0]} sensitivity rows for "
                    f"{len(values[s_])} calculated points.")
            var = np.einsum("ij,jk,ik->i", g, C, g)
            out[key_of[seg]] = Band(
                segment=seg, key=key_of[seg],
                energy=s.calculated_energies(seg),
                angle=s.calculated_angles(seg),
                excitation=s.calculated_excitation_energies(seg),
                value=values[s_], sigma=np.sqrt(np.clip(var, 0.0, None)),
                quantity=quantity, sensitivity=g)
    return out


# -- the convenient path: a dedicated, trimmed session ------------------------

def trimmed_model(azr_file, keys=None, grids=None, cwd=None, workdir=None):
    """A ``.azr`` carrying only the requested extrapolation grids.

    Returns ``(path, workdir)``: the trimmed project (its extrapolations are just ``keys`` / ``grids``), written into
    ``workdir``, which is also its output directory.  ``workdir`` defaults to a new temporary directory; the caller
    owns it and removes it (``shutil.rmtree(workdir)``).  Nothing is written beside the project.

    A session of the trimmed project takes ``cwd`` (default: the project's directory) as its working directory, so
    that it reads the project's data and tables there: ``Session(path, cwd=cwd, extrapolations=True)``.

    The private output directory keeps the project's own ``output/`` as it is: the external-capture cache of the
    trimmed grids (``intEC.extrap``) is written there, not over the project's.  The engine checks a cache against
    the grid it was computed for (its ``.sig`` file) and recomputes it when the grid has changed, so a ``workdir``
    the caller keeps between calls reuses its integrals while the grid is the same.  The project's data-segment
    cache (``output/intEC.dat`` under ``cwd``, valid since the data sets are not trimmed) is copied in when there is
    one, which saves recomputing it.
    """
    if (keys is None) == (grids is None):
        raise ValueError("give exactly one of keys= or grids=.")
    if cwd is None:
        cwd = os.path.dirname(os.path.abspath(azr_file)) or "."
    cwd = os.path.abspath(cwd)

    model_file, _ = _trimmed_project(azr_file, keys, grids)
    workdir = os.path.abspath(workdir) if workdir is not None else tempfile.mkdtemp(prefix="pyazr_bands_")
    os.makedirs(workdir, exist_ok=True)
    model_file.set_config(output_dir=workdir.rstrip("/") + "/")

    for cache in ("intEC.dat", "intEC.dat.sig"):
        seed = os.path.join(cwd, "output", cache)
        target = os.path.join(workdir, cache)
        if os.path.exists(seed) and not os.path.exists(target):
            shutil.copyfile(seed, target)

    fd, path = tempfile.mkstemp(suffix=".azr", dir=workdir)
    os.close(fd)
    model_file.save(path)
    return path, workdir


def _trimmed_project(azr_file, keys, grids):
    """A format-3 project with only the extrapolations ``keys`` (1-based, in file order; kept in the order given)
    or the extrapolation records ``grids`` (dicts of :meth:`pyazr.Project.add_extrapolation` fields, with an
    optional ``name``), and the text the band directory is named from (the kept records)."""
    from .project import Project
    model = Project.open(azr_file)
    names = [n for n, r in model.records().items() if r.record_kind == "extrapolation"]
    if keys is not None:
        bad = [k for k in keys if not 1 <= k <= len(names)]
        if bad:
            raise KeyError(f"{azr_file} declares {len(names)} extrapolations; no {bad}.")
        kept = [(names[k - 1], model[names[k - 1]].fields) for k in keys]
    else:
        kept = []
        for i, g in enumerate(grids, start=1):
            g = dict(g)
            kept.append((g.pop("name", f"band_{i}"), g))
    for n in names:
        model.remove(n)
    for name, fields in kept:
        model.add("extrapolation", name, **fields)
    block = "\n".join(model[name].text for name, _ in kept)
    return model, block


def extrapolation_bands(azr_file, keys=None, grids=None, params=None,
                        covariance=None, quantity="cross",
                        step=1e-4, method="analytic", cwd=None,
                        progress=False, keep_tempfile=False, nprocs=1, workdir=None):
    """Uncertainty bands on extrapolation grids, in a dedicated session.

    Writes a temporary ``.azr`` whose extrapolations are **only**
    the requested grids (see :func:`trimmed_model`), opens its own AZURE2
    session on it, and computes the band there.  Every finite-difference step
    re-evaluates all active segments, so this is what keeps the cost
    proportional to what you asked for rather than to the whole model.

    Parameters
    ----------
    azr_file : the model.
    keys : 1-based extrapolation numbers of the file to keep, in file order (as reported
        by :attr:`~pyazr.Session.extrapolations`).  The returned dict is
        keyed by these same numbers.
    grids : alternatively, an iterable of grids the file does not contain, keyed 1, 2, ... in the order given:
        dicts of :meth:`pyazr.Project.add_extrapolation` fields (``reaction``, ``quantity``, ``unit``,
        ``energies``, ...; an optional ``name``).
    params : defaults to ``output/param.sav`` beside the model.
    covariance, quantity, step, method, progress, nprocs : see
        :func:`uncertainty_bands`.
    cwd : working directory for the session (its data and tables are read there); defaults to the ``.azr``'s.
    workdir : the directory the trimmed project and its outputs are written to (:func:`trimmed_model`); default a
        temporary directory, removed at the end.  A directory of the caller's is kept (only the trimmed ``.azr`` is
        removed), so that a second call on the same grid reuses its external-capture integrals.
    keep_tempfile : leave the trimmed ``.azr`` and the temporary directory on disk (for debugging).

    Notes
    -----
    Nothing is written beside the project: the trimmed model writes into ``workdir``, so the project's ``output/``
    (with its external-capture caches) is left alone.  Building the integrals of a capture grid is the bulk of the
    start-up cost; a ``workdir`` kept between calls on the same grid saves it.
    """
    from .session import Session

    if cwd is None:
        cwd = os.path.dirname(os.path.abspath(azr_file)) or "."
    out_keys = list(keys) if keys is not None else list(
        range(1, len(list(grids)) + 1))
    own_workdir = workdir is None
    tmp, workdir = trimmed_model(azr_file, keys=keys, grids=grids, cwd=cwd, workdir=workdir)

    try:
        # Without params, the session takes the fit's parameter file, fixed values included (best_fit_params).
        sav = None if params is not None else os.path.join(cwd, "output", "param.sav")
        if sav is not None and not os.path.exists(sav):
            raise MissingFile(sav, f"{sav} not found: the band must be evaluated at the fitted parameters the "
                                   f"covariance belongs to.  Pass params= explicitly to override.",
                              what=f"{sav} not found",
                              allowed="the band must be evaluated at the fitted parameters the covariance belongs to",
                              fix="Pass params= explicitly to override")
        with Session(tmp, cwd=cwd, extrapolations=True, parameter_file=sav) as m:
            if params is None:
                params = best_fit_params(m, sav)
            if covariance is None:
                covariance = os.path.join(cwd, "output", "covariance.dat")
            bands = uncertainty_bands(
                m, keys=list(range(1, len(out_keys) + 1)), params=params,
                covariance=covariance, quantity=quantity, step=step,
                method=method, progress=progress, nprocs=nprocs)
    finally:
        if not keep_tempfile:
            if own_workdir:
                shutil.rmtree(workdir, ignore_errors=True)
            elif os.path.exists(tmp):
                os.remove(tmp)

    # Re-key from the trimmed file's 1..N numbering to the caller's.
    return {out_keys[b.segment]: _rekey(b, out_keys[b.segment])
            for b in bands.values()}


def _rekey(band, key):
    band.key = key
    return band
