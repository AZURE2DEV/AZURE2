"""Covariance-propagated uncertainty bands for AZURE2 observables.

A fit leaves behind two things: the best-fit parameters (``output/param.sav``)
and their covariance (``output/covariance.dat``).  The 1-sigma band on any
observable the model can compute follows from linear propagation,

.. math::  \\sigma_y^2 = g^{\\mathsf T} C\\, g, \\qquad g_j = \\partial y/\\partial p_j

which is the same expression AZURE2 uses internally to write its ``.band``
files (SAMMY Eq. IV E4.2).  ``covariance.dat`` spans exactly the **free
R-matrix parameters** -- level energies and reduced-width amplitudes, in their
``params_rwa`` order -- because an observable is insensitive to the dataset
normalizations and energy shifts that make up the rest of the fit vector.

The sensitivities :math:`g` come from AZURE2 itself -- one reverse-mode adjoint
per point, so the whole matrix costs about two forward evaluations however many
parameters the model has.  Finite differences are available as a cross-check
(``method='central'``), at ``2 * ncols`` forward passes.

Two entry points:

``uncertainty_bands(model, segments, ...)``
    band for selected segments of a *live* :class:`~pyazr.azure2.azure2`
    session, in whichever mode it is in (data or extrapolation).  The session
    is the in-process engine itself -- one ``_azure2.Session``, no subprocess
    or instance pool to talk to.

``extrapolation_bands(azr_file, keys=..., grids=..., ...)``
    band on an extrapolation grid, computed in a **dedicated session whose
    ``.azr`` carries only the requested segments**.  Each finite-difference
    adjoint pass covers every active segment of the model, so trimming the
    ``<segmentsTest>`` block to what you actually want to plot still pays:
    on the full 7Be model one pass over all 110 grids is ~31 s, against
    ~0.1 s for a handful of them.

Examples
--------
>>> from pyazr import azure2, extrapolation_bands
>>> bands = extrapolation_bands("7Be.azr", keys=[113], quantity="sfactor")
>>> b = bands[113]
>>> b.value[0], b.sigma[0]                          # S(E->0) +- 1 sigma
"""

import os
import shutil
import sys
import warnings
from dataclasses import dataclass, field

import numpy as np

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

    ``azure2`` switches the process into the ``.azr``'s directory for the
    session's lifetime, and AZURE2 names ``output/`` relative to that -- so
    that is where ``param.sav`` and ``covariance.dat`` are.
    """
    return os.path.abspath(model.cwd or os.path.dirname(model.file) or ".")


def load_covariance(path="output/covariance.dat", ncols=None):
    """Read AZURE2's saved parameter covariance (a bare ``N x N`` matrix).

    ``covariance.dat`` is written by a fit run with the uncertainty band
    enabled and holds the free **R-matrix** block only, in packed parameter
    order.  ``ncols``, when given, is the number of free R-matrix parameters
    the current model has; a mismatch means the covariance belongs to a
    different model and is refused rather than silently misaligned.
    """
    C = np.loadtxt(path, ndmin=2)
    if C.ndim != 2 or C.shape[0] != C.shape[1]:
        raise ValueError(f"{path}: expected a square matrix, got {C.shape}.")
    if ncols is not None and C.shape[0] != ncols:
        raise ValueError(
            f"{path} is {C.shape[0]}x{C.shape[0]} but the model has {ncols} "
            f"free R-matrix parameters.  The covariance was saved for a "
            f"different model (or a different set of fixed parameters); re-run "
            f"the fit with the uncertainty band enabled.")
    return C


_RMATRIX_KINDS = (0, 1)      # GET_PARAMS_INFO type codes: 0 energy, 1 width
_NFIELDS = 15                # doubles per parameter record


def live_parameters(model):
    """``(free_kinds, free_identities)`` queried from the running session.

    ``model.parameters`` is built once, from whatever mode the session was in
    at the time -- but the parameter vector is **not** the same in the two
    modes.  Extrapolation segments carry no normalizations, and AZURE2 drops
    the fixed R-matrix parameters its grids cannot reach, so a model with 483
    parameters (197 free) in data mode can report 234 (72 free) in
    extrapolation mode.  The free R-matrix block survives, in order, and it is
    that block the covariance describes.  This queries the engine's *current*
    parameter bookkeeping (``params_fixed`` / ``parameter_info``), which is
    rebuilt on every ``initialize()``, so it is right for the mode the session
    is in now.

    Returns the per-free-parameter ``kind`` codes and their identity tuples
    ``(kind, jgroup, level, channel)`` -- the same key AZURE2 uses to match
    covariance columns across runs.
    """
    s = model.sess
    fixed = np.asarray(s.params_fixed(), float).round().astype(bool)
    info = np.asarray(s.parameter_info(), float).reshape(len(fixed), _NFIELDS)
    kinds, ids = [], []
    for i in range(len(fixed)):
        if fixed[i]:
            continue
        k = int(round(info[i, 0]))
        kinds.append(k)
        ids.append((k, int(round(info[i, 1])), int(round(info[i, 4])),
                    int(round(info[i, 6]))))
    return kinds, ids


def rmatrix_columns(model):
    """Packed indices of the free R-matrix parameters, in ``params_rwa`` order.

    These are the columns ``covariance.dat`` spans: level energies and reduced
    width amplitudes.  Normalizations and energy shifts are excluded -- no
    calculated observable depends on them.  Queried live, so the answer is
    right for the mode the session is in now.
    """
    kinds, _ = live_parameters(model)
    cols = [j for j, k in enumerate(kinds) if k in _RMATRIX_KINDS]
    if not cols:
        raise RuntimeError("the model has no free R-matrix parameters.")
    return cols


def best_fit_params(model, path=None):
    """The best-fit free parameter vector, read from ``output/param.sav``.

    A model's ``.azr`` file need not carry the fitted values -- ``param.sav``
    is what a fit writes -- so a band must be evaluated at the parameters the
    covariance belongs to.  On the 7Be model the difference is not subtle: the
    ``.azr``'s own parameters give chi2 = 173006, ``param.sav`` gives 4130.

    ``param.sav`` is written in the *data*-mode layout (one line per parameter,
    fixed ones included).  In extrapolation mode the session holds a shorter
    vector -- see :func:`live_parameters` -- so the leading R-matrix block is
    taken, after checking that it is what the running session expects.
    """
    if path is None:
        path = os.path.join(_model_dir(model), "output", "param.sav")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found: the band must be evaluated at the fitted "
            f"parameters the covariance belongs to.  Pass params= explicitly "
            f"to override.")
    full = np.loadtxt(path, usecols=(1,))
    fixed = np.asarray(model.fixed_params, float).round().astype(bool)
    if full.size != fixed.size:
        raise ValueError(f"{path} has {full.size} parameters but the model has "
                         f"{fixed.size}.")
    free = full[~fixed]

    kinds, _ = live_parameters(model)
    if len(kinds) == free.size:
        return free
    # Fewer parameters live than param.sav describes: extrapolation mode.  The
    # R-matrix block leads the vector in both layouts, so the leading slice is
    # the same parameters -- but only if nothing but R-matrix parameters is
    # left, which is exactly the extrapolation-mode case.
    if len(kinds) < free.size and all(k in _RMATRIX_KINDS for k in kinds):
        return free[:len(kinds)]
    raise ValueError(
        f"{path} gives {free.size} free parameters but the running session "
        f"expects {len(kinds)}, and they are not the leading R-matrix block.  "
        f"Pass params= explicitly.")


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
    s = model.sess
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
    from .azrfile import AzrModel
    from .azure2 import azure2

    cwd, workdir = spec["cwd"], os.path.join(spec["workdir"], str(os.getpid()))
    os.makedirs(workdir, exist_ok=True)
    # Best effort: a missing seed only means this worker rebuilds its own
    # integrals, which is correct, just slower.
    for cache in ("intEC.dat", "intEC.extrap"):
        seed = os.path.join(cwd, "output", cache)
        if os.path.exists(seed):
            shutil.copyfile(seed, os.path.join(workdir, cache))

    # The edited .azr lives in the private directory too, so the parent removing
    # that directory takes the temporary file with it.
    edited = AzrModel.from_file(spec["file"])
    edited.set_output_dir(os.path.relpath(workdir, cwd))
    path = edited.to_tempfile(dir=workdir)

    model = azure2(path, cwd=cwd, **spec["options"])
    if spec["mode"] == "extrap":
        model.extrap_mode()
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
        "mode": getattr(model, "mode", "data"),
        "workdir": os.path.join(model.cwd, "output", "sensitivities"),
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
    model : a live :class:`~pyazr.azure2.azure2` session.
    segments : calculated-segment indices (0-based, counting **active**
        segments only -- see :func:`active_indices`).
    params : the free parameter vector to differentiate about; defaults to
        ``model.params_rwa``.
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
        extrapolation mode, one grid rebuild -- before it evaluates anything.
        The columns have to be worth more than that.  As a rule of thumb it
        needs ``len(cols)`` in the hundreds: on the 14-parameter ``tests/13N``
        model it is a net loss at any ``nprocs``, in data mode (0.8 s serial)
        and on an 1845-point extrapolation grid (53 s serial) alike.  Measure
        on your own model rather than assuming.

    Returns
    -------
    dict
        ``{segment: array of shape (npoints, len(cols))}``.
    """
    p0 = np.asarray(model.params_rwa if params is None else params, float)

    if method == "analytic":
        G = model.model_gradients(p0)
        bad = [i for i in segments if i < 0 or i >= len(G)]
        if bad:
            raise IndexError(f"segment index/indices {bad} out of range: the "
                             f"model calculates {len(G)} (active) segments.")
        if quantity == "sfactor":
            # S = conv * sigma with conv fixed by the kinematics, so the
            # sensitivities scale by the same factor, row by row.
            _evaluate(model, p0, segments, quantity)
            return {seg: G[seg] * model.sess.calculated_conv(seg)[:, None]
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
    """Map file order to calculated order for a :class:`TestSegmentSet`/
    :class:`SegmentSet`.

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
    model : a live :class:`~pyazr.azure2.azure2`, in data or extrapolation
        mode; the bands are computed for whichever segments that mode exposes.
    segments : calculated-segment indices (0-based).  Mutually exclusive with
        ``keys``.
    keys : 1-based segment numbers as they appear in the ``.azr`` file, which
        is what :attr:`~pyazr.azure2.azure2.extrapolations` /
        :attr:`~pyazr.azure2.azure2.datasets` report.  Translated to calculated
        indices for you (inactive segments are not calculated).
    params : parameters to evaluate at; defaults to ``output/param.sav``, since
        that -- not the ``.azr`` -- is what a fit writes and what the
        covariance describes.  Pass ``model.params_rwa`` to use the file's own
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
    cols = [j for j, k in enumerate(kinds) if k in _RMATRIX_KINDS]
    if not cols:
        raise RuntimeError("the model has no free R-matrix parameters.")

    C = covariance
    if C is None:
        C = os.path.join(_model_dir(model), "output", "covariance.dat")
    if isinstance(C, (str, os.PathLike)):
        C = load_covariance(C, ncols=len(cols))
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
            f"{len(kinds)} free parameters.  In extrapolation mode the vector "
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
    values = _evaluate(model, p0, segments, quantity)
    s = model.sess

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

    Returns ``(path, workdir)``: a temporary ``.azr`` (the caller owns it) whose
    ``<segmentsTest>`` block is just ``keys`` / ``grids``, and the private
    output directory it was pointed at.

    The private output directory is what makes this safe *and* fast.  AZURE2
    caches its external-capture integrals in ``<output>/intEC.extrap`` and reads
    them back even when the grid has changed -- silently mixing integrals from
    another set of energies into the calculation.  Giving each grid its own
    directory (named from a hash of the trimmed block) means the main run's
    cache is never touched, and a repeat call on the same grid reuses its own
    valid cache instead of rebuilding it.
    """
    from .azrfile import AzrModel
    from .datasets import TestSegmentSet
    import hashlib

    if (keys is None) == (grids is None):
        raise ValueError("give exactly one of keys= or grids=.")
    if cwd is None:
        cwd = os.path.dirname(os.path.abspath(azr_file)) or "."

    model_file = AzrModel.from_file(azr_file)
    if keys is not None:
        n = len(TestSegmentSet.from_file(azr_file))
        bad = [k for k in keys if not 1 <= k <= n]
        if bad:
            raise KeyError(f"{azr_file} declares {n} <segmentsTest> segments; "
                           f"no {bad}.")
        model_file.keep_extrapolations(list(keys))
    else:
        model_file.set_extrapolations(list(grids))

    if workdir is None:
        block = model_file._suffix.split("<segmentsTest>")[1].split(
            "</segmentsTest>")[0]
        tag = hashlib.sha1(block.encode()).hexdigest()[:10]
        workdir = os.path.join(cwd, "output", "bands", tag)
    model_file.set_output_dir(os.path.relpath(workdir, cwd))

    # Trimming touches <segmentsTest> only, so the data-segment external-capture
    # integrals are unchanged and the main run's cache is valid here.  Seeding
    # it saves rebuilding them -- the bulk of a capture model's startup.
    seed = os.path.join(cwd, "output", "intEC.dat")
    target = os.path.join(workdir, "intEC.dat")
    if os.path.exists(seed) and not os.path.exists(target):
        shutil.copyfile(seed, target)

    return model_file.to_tempfile(dir=cwd), workdir


def extrapolation_bands(azr_file, keys=None, grids=None, params=None,
                        covariance=None, quantity="cross",
                        step=1e-4, method="analytic", cwd=None,
                        progress=False, keep_tempfile=False, nprocs=1):
    """Uncertainty bands on extrapolation grids, in a dedicated session.

    Writes a temporary ``.azr`` whose ``<segmentsTest>`` block holds **only**
    the requested grids (see :func:`trimmed_model`), opens its own AZURE2
    session on it, and computes the band there.  Every finite-difference step
    re-evaluates all active segments, so this is what keeps the cost
    proportional to what you asked for rather than to the whole model.

    Parameters
    ----------
    azr_file : the model.
    keys : 1-based ``<segmentsTest>`` numbers of the file to keep (as reported
        by :attr:`~pyazr.azure2.azure2.extrapolations`).  The returned dict is
        keyed by these same numbers.
    grids : alternatively, an iterable of keyword bundles for
        :meth:`~pyazr.AzrModel.add_extrapolation` -- a grid the file does not
        contain.  Keyed 1, 2, ... in the order given.
    params : defaults to ``output/param.sav`` beside the model.
    covariance, quantity, step, method, progress, nprocs : see
        :func:`uncertainty_bands`.
    cwd : working directory for the session; defaults to the ``.azr``'s.
    keep_tempfile : leave the trimmed ``.azr`` on disk (for debugging).

    Notes
    -----
    The trimmed model writes to its own ``output/bands/<hash>/`` directory, so
    the main run's ``output/`` -- including the external-capture cache
    ``intEC.extrap``, which AZURE2 would otherwise reuse on the changed grid --
    is left alone.  Building those integrals is the bulk of the startup cost
    for a capture grid, so the second call on the same grid is much faster.
    """
    from .azure2 import azure2

    if cwd is None:
        cwd = os.path.dirname(os.path.abspath(azr_file)) or "."
    out_keys = list(keys) if keys is not None else list(
        range(1, len(list(grids)) + 1))
    tmp, workdir = trimmed_model(azr_file, keys=keys, grids=grids, cwd=cwd)

    try:
        with azure2(tmp, cwd=cwd) as m:
            m.extrap_mode()
            if params is None:
                params = best_fit_params(
                    m, os.path.join(cwd, "output", "param.sav"))
            if covariance is None:
                covariance = os.path.join(cwd, "output", "covariance.dat")
            bands = uncertainty_bands(
                m, keys=list(range(1, len(out_keys) + 1)), params=params,
                covariance=covariance, quantity=quantity, step=step,
                method=method, progress=progress, nprocs=nprocs)
    finally:
        if not keep_tempfile and os.path.exists(tmp):
            os.remove(tmp)

    # Re-key from the trimmed file's 1..N numbering to the caller's.
    return {out_keys[b.segment]: _rekey(b, out_keys[b.segment])
            for b in bands.values()}


def _rekey(band, key):
    band.key = key
    return band
