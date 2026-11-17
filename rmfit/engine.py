"""Evaluators: the objective, residuals and Jacobian of a model at a free vector.

Two implementations share one interface:

- :class:`SingleSession` -- one engine session on the full model.  Slow
  (23 s per evaluation on a large model) but authoritative: every
  number that is recorded as *verified* comes from here.
- :class:`Sharded` -- the data split across N worker processes (see
  :mod:`rmfit.model`), each evaluating its share; residuals are concatenated,
  Jacobian columns of split segments summed, and AZURE2's normalization /
  energy-shift penalties (and optional prior rows) appended once, here, from
  the base model's segment table rather than from pyazr.

Both return residual vectors whose sum of squares *is* the objective, so a
least-squares solver sees the penalties' Jacobian too.
"""
import hashlib
import json
import multiprocessing as mp
import os
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import AZURE2_ROOT  # noqa: F401
from .keys import ParamKey, index_map, build_keys
from .model import ShardSpec, build_shards, base_segments, single_variant, plan_balanced, shard_shift_count
from .covariance import Covariance, NormGroups, RowBlock, has_keyword, load_reference
from .azr3 import SegmentTable, has_field, without_fields


class EvalError(RuntimeError):
    """The engine could not evaluate this vector (non-finite, crash, ...)."""


#: order of a segment's shift-like keys in Sharded's key tail (norms first, then by segment)
_TAIL_ORDER = {"shift": 0, "shift_sqrt": 1, "shift_fn": 2}


def _rows_key(rows) -> bytes:
    """A content key of a data-row vector (Evaluator._cov_memo)."""
    return hashlib.blake2b(np.ascontiguousarray(rows, float).tobytes(), digest_size=16).digest()


@dataclass
class Penalties:
    """Nominal values and sigmas of the norm / shift penalty rows.

    AZURE2 adds ``((p - centre) / sigma)^2`` for every free nuisance with a prior (a normalization, an energy-shift
    term, ...); centre and sigma are the engine's own, read from the session's parameter rows
    (:meth:`from_nuisances`).
    """
    cols: np.ndarray            # free-vector positions
    nominal: np.ndarray
    sigma: np.ndarray

    @classmethod
    def from_nuisances(cls, keys: Sequence[ParamKey], rows, key_of=None):
        """The engine's Gaussian priors of the free parameters as rows: one ``(x - centre)/sigma`` per free
        parameter whose row of the session's table has ``prior_sigma`` > 0 (``rows``: :func:`rmfit.azr3.nuisance_rows`
        of a session whose fit vector has the order of ``keys``; ``key_of``: for the rows of another session, a
        function row -> ParamKey or None).  Normalizations, energy-shift terms, shift-function parameters, target
        densities, Doppler temperatures, sample thicknesses, screening potentials, abundances and the priors of
        ``<parameterSettings>`` all come this way, with the centre the engine uses.  The engine's terms that couple
        several parameters (the abundances' sum, a prior covariance) are not rows here."""
        where = {k: j for j, k in enumerate(keys)}
        got: Dict[int, Tuple[float, float]] = {}
        for r in rows:
            sigma = r.get("prior_sigma") or 0.0
            if r["fixed"] or r.get("free_index") is None or not sigma > 0:
                continue
            j = int(r["free_index"]) if key_of is None else where.get(key_of(r))
            if j is None:
                continue
            got[int(j)] = (float(r["prior_centre"]), float(sigma))
        cols = sorted(got)
        return cls(np.array(cols, int), np.array([got[j][0] for j in cols], float),
                   np.array([got[j][1] for j in cols], float))

    def with_param_priors(self, keys: Sequence[ParamKey], priors) -> "Penalties":
        """Append Gaussian rows ``((x - mean) / sigma)^2`` on named free parameters.

        ``priors``: ``[{"key": "g[1/2+#2 p2 L1 S0.5]", "mean": 0.537, "sigma": 0.062}, ...]`` with
        ``key`` the parameter's ``str(ParamKey)``.  The row acts on the free vector itself, i.e. on
        the reduced-width amplitude for a width (Brune: the observed amplitude, the space AZURE2's own
        nuisance priors use) and on MeV for an energy; give ``mean`` with the amplitude's sign.  Like
        the norm rows the rows are never weighted, and they count in ``penalty``, not in ``chi2``.
        A prior on a key that is not free is an error, so a frozen parameter cannot hide a prior.
        """
        if not priors:
            return self
        where = {str(k): j for j, k in enumerate(keys)}
        cols, nom, sig = list(self.cols), list(self.nominal), list(self.sigma)
        for pr in priors:
            key = str(pr["key"])
            if key not in where:
                raise KeyError(f"param_priors: {key} is not a free parameter of this model")
            s = float(pr["sigma"])
            if not s > 0:
                raise ValueError(f"param_priors: {key} needs sigma > 0")
            cols.append(where[key]); nom.append(float(pr["mean"])); sig.append(s)
        return Penalties(np.array(cols, int), np.array(nom, float), np.array(sig, float))

    def residuals(self, x):
        if self.cols.size == 0:
            return np.zeros(0)
        return (np.asarray(x, float)[self.cols] - self.nominal) / self.sigma

    def jacobian(self, nfree):
        J = np.zeros((self.cols.size, nfree))
        for r, (c, s) in enumerate(zip(self.cols, self.sigma)):
            J[r, c] = 1.0 / s
        return J

    def value(self, x):
        return float(np.sum(self.residuals(x) ** 2))


@dataclass
class Prior:
    """Optional Gaussian rows on width amplitudes: ``x_w / (scale * gamma_W)``.

    A regularizer for exploration only (it damps the flat directions that let
    an unbounded fit run away); switch it off for every reported number.
    """
    cols: np.ndarray
    sigma: np.ndarray

    def residuals(self, x):
        return np.asarray(x, float)[self.cols] / self.sigma

    def jacobian(self, nfree):
        J = np.zeros((self.cols.size, nfree))
        for r, (c, s) in enumerate(zip(self.cols, self.sigma)):
            J[r, c] = 1.0 / s
        return J

    def value(self, x):
        return float(np.sum(self.residuals(x) ** 2))


@dataclass
class ParkWall:
    """Park's J > 0 wall as residual rows (2026-10-24, the evaluator: "add the J rows").

    In Park's parametrization (flags ``use_park``) every level has the overlap
    ``J = 1 - sum_c gamma_c^2 dS_c/dE`` (pyazr ``park_norms``); a parameter set with J <= 0 has
    no R-matrix meaning.  The engine adds ``sum (J/1e-3)^2`` over levels with J < 0 to its own
    scalar chi-squared, but rmfit builds its objective from residual rows, so the wall is one
    row per level here, ``r_l = min(J_l, 0) / 1e-3`` -- zero while J > 0, the engine's term
    beyond it.  J of a level depends on that level's own energy and amplitudes only; its
    Jacobian rows come from the worker's ``park_jac`` (central differences of park_norms).
    The rows go last (after the penalty and prior rows), are never weighted, and count in
    ``penalty``.
    """
    labels: List[str]                 # level label per J ("5/2+#3"), park_norms order
    cols: List[List[int]]             # free R-matrix columns of each level
    scale: float = 1e-3

    @classmethod
    def build(cls, keys, n_rmatrix, labels):
        if not labels:
            raise RuntimeError("use_park: the worker reported no levels for the J > 0 wall")
        by = {}
        for j, k in enumerate(keys[:n_rmatrix]):
            by.setdefault(k.level, []).append(j)
        return cls(list(labels), [by.get(lab, []) for lab in labels])


@dataclass
class ModelRows:
    """The rows of the engine's objective that belong to the model as a whole, not to a data set and not to one
    nuisance's prior: model data (measured scattering lengths, widths, branching ratios, ... of ``<modelData>``), a
    prior covariance, Park's J > 0 wall, the priors of ``gamma_decay`` records and the abundance sums.

    pyazr 3.0's ``Objective.residuals`` gives them after the data rows and before the nuisance priors; the worker
    cuts that block out (``model_rows``) and the evaluator appends it once, from one worker that holds the whole
    R-matrix block (a shard keeps the model data and the decay records, so any shard has it).  The rows are never
    weighted and count in ``penalty``.  pyazr 2 had no such rows: rmfit rebuilt the scattering lengths by finite
    differences (ScatLenRows) and left the other terms out."""
    n: int

    @property
    def parts(self):
        return range(self.n)


@dataclass
class ObjectiveSpec:
    """How data residuals enter the objective (campaign.json ``objective``).

    ``dataset_weight``: ``"none"`` (plain chi-squared), ``"scaled"`` (2026-09-28:
    every data set's error bars inflated so that its chi2/N is 1 at the
    reference model -- row weight ``1 / max(scales[i], 1)`` with ``scales`` the
    baseline chi2/N per data set that ``Campaign.init`` records; the objective
    then has the magnitude of N and the level threshold tau = s^2 k ln N runs
    with s^2 ~ 1, set by every data set alike instead of the worst-fitted ones)
    or ``"reduced"`` -- the
    reduced-chi-squared method of deBoer et al. (2017, Eq. for the EDA-style
    fit): every data set's chi-squared is divided by its number of points so
    that sets with many points do not dominate.  The row weights are
    ``(N / K) / max(N_i - nu, 1)`` for the ``N_i`` points of data set ``i``
    (``K`` data sets, ``N`` points in all), i.e. the paper's ``1/(N_i - nu)``
    rescaled so the objective keeps the magnitude of a chi-squared (thresholds
    such as ``tau_flip`` stay meaningful).  A data set is a data file
    (``dataset_key="file"``) or a segment (``"segment"``).  Normalization and
    shift penalty rows are never weighted.  ``weight_power`` tempers the
    weighting: ``w_i = ((N/K)/N_i)^p`` renormalized to a total of ``N``; ``p = 1``
    is the paper's method, ``p = 0.5`` a compromise (on a large model the literal method
    spans per-point weights over three decades across many files and trades the three
    10,000-point excitation functions for the small sets: plain chi2 up
    2.7x in one polish).

    ``loss``: ``"chi2"`` or ``"sivia"`` -- the broader probability density of
    Sivia (2006) used in deBoer et al. (2017): instead of chi-squared,
    ``-2 log[(1 - exp(-R^2/2)) / R^2] - 2 log 2`` is summed over the data
    rows (zero at R = 0, ~R^2/2 for small residuals, ~2 log R^2 for outliers:
    an outlier is down-weighted like a point with an inflated error bar).
    Applied to the weighted residuals; penalty rows stay quadratic.

    ``loss_scale`` (sivia only): the residuals of data set ``i`` are
    standardized by its reference reduced chi-squared ``s_i^2`` before the
    loss, ``c rho(z / c)`` with ``c = w_i s_i^2`` per row, so a residual
    counts as an outlier only relative to the *typical* residual of its own
    data set.  Without this (``"none"``), on a model with chi2/N ~ 20-300 every
    point sits in the logarithmic tail, the data lose almost all weight
    against the quadratic norm penalties and the fit drifts (chi2
    up by 70 % in one polish, a data set abandoned).  ``"global"`` uses one scale, the
    overall chi2/N.  ``scales`` holds the reference values (``Campaign.init``
    computes them from the baseline's plain per-dataset chi-squared).
    """
    dataset_weight: str = "none"
    dataset_key: str = "file"
    loss: str = "chi2"
    nu: int = 0
    loss_scale: str = "dataset"     # sivia: "dataset" (each set's reference chi2/N), "global", or "none"
    scales: dict = field(default_factory=dict)   # dataset -> reference chi2/N (filled by Campaign.init)
    weight_power: float = 1.0       # reduced: w_i = ((N/K)/N_i)^p, renormalized to sum N; 1 = the paper, 0.5 = tempered
    point_weights: list = field(default_factory=list)   # [[ex_lo, ex_hi, w], ...] (MeV): data rows whose excitation
                                    # energy lies in a window have their chi2 term multiplied by w (windows multiply
                                    # where they overlap).  the evaluator 2026-09-28: up-weight the resonance peaks while a
                                    # narrow level is fitted so it cannot run off into the smooth regions.  Penalty
                                    # rows are never weighted; per-segment chi2 tables stay plain.
    param_priors: list = field(default_factory=list)    # [{"key": str(ParamKey), "mean": m, "sigma": s}, ...]:
                                    # Gaussian penalty rows on free parameters (Penalties.with_param_priors), e.g. an
                                    # ANC prior on a near-threshold alpha amplitude (the evaluator 2026-10-20).
    cov_reference: str = ""         # .npz of covariance.save_reference: the data covariance's model-scaled (relvector)
                                    # terms are scaled by this FIXED model at every evaluation instead of the model of
                                    # the evaluator's first evaluation (the evaluator 2026-10-25, after the refreeze iteration
                                    # oscillated on a large model).  Campaign.init writes the seed's (cov_reference.npz).
    cov_free_norm_ok: list = field(default_factory=list)   # covariance file basenames whose relvector rows may carry
                                    # a FREE data normalization.  Campaign.init refuses any other (2026-10-26, the evaluator:
                                    # "build the t0 version"): a free norm on relvector rows is biased whichever way
                                    # the reference follows it (a norm ran up fivefold data-unit, down threefold
                                    # norm-following); the t0 method fixes the norm and puts its prior in C as a
                                    # relvector term (NNPDF, JHEP 05 (2010) 075; an archived campaign).

    @property
    def plain(self) -> bool:
        return self.dataset_weight == "none" and self.loss == "chi2" and not self.point_weights

    @property
    def label(self) -> str:
        return ("chi2" if self.dataset_weight == "none" else f"chi2 with errors scaled per {self.dataset_key} (baseline chi2/N = 1)"
                if self.dataset_weight == "scaled" else f"reduced-chi2 per {self.dataset_key}"
                + ("" if self.weight_power == 1.0 else f" (power {self.weight_power:g})")) + \
               ("" if self.loss == "chi2" else f" with {self.loss} loss" +
                ("" if self.loss_scale == "none" else f" ({self.loss_scale}-scaled)")) + \
               ("" if not self.point_weights else " with " + ", ".join(
                   f"x{float(w[2]):g} on Ex {float(w[0]):.4g}-{float(w[1]):.4g}" for w in self.point_weights)) + \
               ("" if not self.param_priors else " + priors " + ", ".join(
                   f"{p['key']}={float(p['mean']):g}({float(p['sigma']):g})" for p in self.param_priors)) + \
               ("" if not self.cov_reference else f" (relvector model fixed: {os.path.basename(self.cov_reference)})")


_LOG2 = float(np.log(2.0))


def rho_sivia(z):
    """Sivia's conservative loss on ``z = R^2``: ``(rho, rho', rho'')`` with
    ``rho = -2 log[(1 - e^{-z/2}) / z] - 2 log 2``."""
    z = np.asarray(z, float)
    rho = np.empty_like(z)
    d1 = np.empty_like(z)
    d2 = np.empty_like(z)
    small = z < 1e-3
    zs = z[small]
    rho[small] = zs / 2 - zs * zs / 48
    d1[small] = 0.5 - zs / 24
    d2[small] = -1.0 / 24
    zb = np.clip(z[~small], 1e-3, 1400.0)
    em = np.expm1(zb / 2)
    rho[~small] = 2 * np.log(zb) - 2 * np.log1p(-np.exp(-zb / 2)) - 2 * _LOG2
    d1[~small] = 2 / zb - 1 / em
    d2[~small] = -2 / zb ** 2 + 0.5 * np.exp(zb / 2) / em ** 2
    big = ~small & (z > 1400.0)
    if np.any(big):                       # beyond the clip: pure 2 log z tail
        rho[big] = 2 * np.log(z[big]) - 2 * _LOG2
        d1[big] = 2 / z[big]
        d2[big] = -2 / z[big] ** 2
    return rho, d1, d2


def rho_sivia_scaled(z, c):
    """``c rho(z / c)`` and its derivatives: the Sivia loss with the residual
    standardized by the scale ``c`` (per row)."""
    z = np.asarray(z, float)
    c = np.asarray(c, float)
    r, d1, d2 = rho_sivia(z / c)
    return c * r, d1, d2 / c


def make_ls_loss(npts: int, scale=None):
    """A ``scipy.optimize.least_squares`` loss callable: (scaled) Sivia on the
    first ``npts`` rows (the data), quadratic on the rest (penalties, priors)."""
    def loss(z):
        z = np.asarray(z, float)
        out = np.empty((3, z.size))
        out[0] = z
        out[1] = 1.0
        out[2] = 0.0
        if npts > 0:
            if scale is None:
                r, d1, d2 = rho_sivia(z[:npts])
            else:
                r, d1, d2 = rho_sivia_scaled(z[:npts], scale)
            out[0, :npts] = r
            out[1, :npts] = d1
            out[2, :npts] = d2
        return out
    return loss


@dataclass
class Score:
    objective: float
    chi2: float
    penalty: float
    prior: float
    npts: int
    segment_chi2: Dict[int, float] = field(default_factory=dict)
    wall: float = 0.0
    roundtrip: Optional[dict] = None    # set by bake(): rwa->physical->rwa noise
    objective_kind: str = "chi2"       # ObjectiveSpec.label of the evaluator that scored it
    chi2_self: Optional[float] = None  # with a fixed covariance reference: the chi2 with the relvector terms frozen at
                                       # this point itself (= a reference-free fresh session's), for comparison only


class _Worker:
    """A restartable worker process holding one engine session."""

    def __init__(self, azr, cwd, flags, seg_map, name="worker", omp_threads=1):
        self.azr, self.cwd, self.flags, self.seg_map, self.name = azr, cwd, flags, seg_map, name
        self.omp_threads = omp_threads
        self.proc = None
        self.conn = None
        self.info = None
        self.restarts = 0

    def start(self):
        from ._worker import worker_main
        # the OpenMP runtime reads this when the engine library loads, which a
        # spawned child does while importing its main module -- before
        # worker_main runs -- so it has to be in the environment it inherits
        os.environ["OMP_NUM_THREADS"] = str(self.omp_threads)
        ctx = mp.get_context("spawn")
        parent, child = ctx.Pipe()
        self.proc = ctx.Process(target=worker_main,
                                args=(child, self.azr, self.cwd, self.flags, self.seg_map),
                                name=self.name, daemon=True)
        self.proc.start()
        child.close()
        self.conn = parent
        self.warnings = []
        while True:
            status, payload = self.conn.recv()
            if status == "warn":
                self.warnings.append(payload)
                continue
            break
        if status != "ready":
            raise RuntimeError(f"{self.name} failed to start on {self.azr}:\n{payload}")
        self.info = payload
        return self

    def stop(self):
        try:
            if self.conn is not None:
                self.conn.send(("quit", None))
                self.conn.recv()
        except Exception:
            pass
        if self.proc is not None:
            self.proc.join(timeout=10)
            if self.proc.is_alive():
                self.proc.kill()
        self.proc = self.conn = None

    def send(self, cmd, payload):
        self.conn.send((cmd, payload))

    def recv(self):
        try:
            status, payload = self.conn.recv()
        except (EOFError, ConnectionResetError, BrokenPipeError) as err:
            raise EvalError(f"{self.name} died ({err.__class__.__name__})")
        if status == "ok":
            return payload
        raise EvalError(f"{self.name}: {payload}")

    def call(self, cmd, payload):
        self.send(cmd, payload)
        return self.recv()

    def restart(self):
        self.restarts += 1
        try:
            if self.proc is not None and self.proc.is_alive():
                self.proc.kill()
        except Exception:
            pass
        self.proc = self.conn = None
        return self.start()


class Evaluator:
    """Common interface (see the module docstring)."""
    keys: List[ParamKey]
    x0: np.ndarray
    npts: int
    n_rmatrix: int
    penalties: Penalties
    prior: Optional[Prior] = None

    spec: "ObjectiveSpec"
    row_sw: Optional[np.ndarray] = None   # sqrt of the data-row weights (None = all one)
    cov: Optional[Covariance] = None      # data covariance groups (None = no `covariance` segments): the
                                          # data rows are whitened here (rmfit.covariance), weights before

    def objective(self, x) -> float:
        raise NotImplementedError

    def residuals(self, x, jac=False, cols=None) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """Residual rows at ``x`` and, with ``jac``, their Jacobian over the free vector.

        ``cols`` (indices or a bool mask over the free vector): the Jacobian columns the caller reads.  The workers'
        finite-differenced columns (energy shifts, sqrt(E) shifts, shift functions, target densities: two engine
        evaluations each) outside it are then not computed and come back zero (2026-10-30; needs a pyazr with
        ``residual_jacobian(skip=)``, AZURE2 dev >= the jac-skip commit -- older builds compute every column).
        ``None``: every column."""
        raise NotImplementedError

    # -- the objective definition (ObjectiveSpec) ------------------------------------------

    def dataset_of_key(self) -> Dict[int, object]:
        """Segment key -> data set id (data file basename, or the key itself)."""
        out = {}
        for d in self.segments:
            out[int(d.key)] = (os.path.basename(d.data_file) if self.spec.dataset_key == "file" else int(d.key))
        return out

    def _init_objective(self, spec):
        self.spec = spec if spec is not None else ObjectiveSpec()
        if self.spec.cov_reference and self.cov is not None and self.cov.relative:
            self.cov.set_reference(load_reference(self.spec.cov_reference))
        if self.cov is not None and self.spec.loss != "chi2":
            raise NotImplementedError(f"loss {self.spec.loss!r} on whitened (covariance) rows is not defined; "
                                      f"use loss 'chi2' with a data covariance")
        if self.spec.param_priors:
            if not hasattr(self.penalties, "with_param_priors"):
                raise NotImplementedError("param_priors are not implemented for this evaluator (rmfit.mixture)")
            self.penalties = self.penalties.with_param_priors(self.keys, self.spec.param_priors)
        self.row_sw = None
        self.row_c = None
        if self.spec.dataset_weight not in ("none", "reduced", "scaled"):
            raise ValueError(f"unknown dataset_weight {self.spec.dataset_weight!r}")
        seg_of_key = self.dataset_of_key()
        n_of = {}
        for key, (a, b) in self.row_ranges:
            ds = seg_of_key.get(int(key), int(key))
            n_of[ds] = n_of.get(ds, 0) + (b - a)
        self.dataset_points = n_of
        K = len(n_of)
        w = np.ones(self.npts)
        if self.spec.dataset_weight == "reduced":
            for key, (a, b) in self.row_ranges:
                ds = seg_of_key.get(int(key), int(key))
                w[a:b] = ((self.npts / K) / max(n_of[ds] - self.spec.nu, 1)) ** float(self.spec.weight_power)
            w *= self.npts / float(np.sum(w))
            self.row_sw = np.sqrt(w)
        if self.spec.dataset_weight == "scaled":
            if not self.spec.scales:
                raise ValueError("ObjectiveSpec: dataset_weight 'scaled' needs reference scales (Campaign.init fills them)")
            for key, (a, b) in self.row_ranges:
                ds = seg_of_key.get(int(key), int(key))
                s2 = self.spec.scales.get(str(ds), self.spec.scales.get("*"))
                if s2 is None:
                    raise KeyError(f"no reference scale for data set {ds!r} (dataset_weight 'scaled')")
                w[a:b] = 1.0 / max(float(s2), 1.0)
            self.row_sw = np.sqrt(w)
        if self.spec.point_weights:
            ex = np.asarray(self.row_excitation(), float)
            if ex.size != self.npts:
                raise RuntimeError(f"point_weights: {ex.size} excitation energies for {self.npts} data rows")
            for win in self.spec.point_weights:
                lo, hi, wt = float(win[0]), float(win[1]), float(win[2])
                if wt <= 0:
                    raise ValueError(f"point_weights: weight must be positive, got {wt}")
                w[(ex >= lo) & (ex <= hi)] *= wt
            self.row_sw = np.sqrt(w)
        if self.spec.loss == "sivia" and self.spec.loss_scale != "none":
            if not self.spec.scales:
                raise ValueError("ObjectiveSpec: loss_scale needs reference scales (Campaign.init fills them)")
            c = np.ones(self.npts)
            for key, (a, b) in self.row_ranges:
                ds = seg_of_key.get(int(key), int(key))
                if self.spec.loss_scale == "global":
                    s2 = self.spec.scales.get("*")
                else:
                    s2 = self.spec.scales.get(str(ds), self.spec.scales.get("*"))
                if s2 is None:
                    raise KeyError(f"no reference scale for data set {ds!r} (loss_scale {self.spec.loss_scale})")
                c[a:b] = w[a:b] * float(s2)
            self.row_c = c

    def row_excitation(self) -> np.ndarray:
        """Excitation energy (MeV) of every data row, in row order (NaN where unknown)."""
        d = self.data()
        out = np.full(self.npts, np.nan)
        exc = d.get("excitation", [])
        for i, (key, (a, b)) in enumerate(self.row_ranges):
            if i < len(exc) and len(exc[i]) == b - a:
                out[a:b] = np.asarray(exc[i], float)
        return out

    def _weight_rows(self, r, J=None):
        """Apply the data-row weights to freshly assembled (plain) rows."""
        if self.row_sw is not None:
            r = r.copy()
            r[:self.npts] *= self.row_sw
            if J is not None:
                J = J.copy()
                J[:self.npts] *= self.row_sw[:, None]
        return r, J

    park: Optional[ParkWall] = None

    def _park_handle(self):
        """``(worker, map)``: a worker holding the whole R-matrix block, and the positions of its
        free vector in the global one."""
        raise NotImplementedError

    def park_norms(self, x) -> np.ndarray:
        """Park's J of every level at ``x`` (``ParkWall.labels`` order)."""
        w, m_ = self._park_handle()
        return np.asarray(w.call("park_norms", np.asarray(x, float)[m_]), float)

    def _park_rows(self, x, jac=False):
        """The J > 0 wall rows (and their Jacobian rows) at ``x``."""
        x = np.asarray(x, float)
        pw = self.park
        J = self.park_norms(x)
        if J.size != len(pw.labels):
            raise RuntimeError(f"park_norms gave {J.size} values for {len(pw.labels)} levels")
        r = np.minimum(J, 0.0) / pw.scale
        if not jac:
            return r, None
        Jm = np.zeros((J.size, x.size))
        neg = np.where(J < 0)[0]
        if neg.size:
            w, m_ = self._park_handle()
            inv = {int(g): s for s, g in enumerate(m_)}
            gcols = sorted({int(c) for l in neg for c in pw.cols[l]})
            if gcols:
                D = np.asarray(w.call("park_jac", (x[m_], [inv[c] for c in gcols])), float)
                where = {c: i for i, c in enumerate(gcols)}
                for l in neg:
                    for c in pw.cols[l]:
                        Jm[l, c] = D[where[int(c)], l] / pw.scale
        return r, Jm

    def park_value(self, r) -> float:
        """The wall's share of the objective, from a residual vector as ``residuals`` returns it."""
        if self.park is None:
            return 0.0
        return float(np.sum(np.asarray(r, float)[-len(self.park.labels):] ** 2))

    scat: Optional["ModelRows"] = None

    def _init_scat(self, worker, m_, azr):
        """The engine's own rows of the model (:class:`ModelRows`), None when it has none (or only Park's wall,
        which :class:`ParkWall` carries)."""
        n = int(worker.call("model_rows", (np.asarray(self.x0, float)[m_], False))[0].size)
        if n == 0:
            self.scat = None
            return
        if self.park is not None:
            if n == len(self.park.labels):
                self.scat = None              # the block is the wall alone: ParkWall's rows
                return
            # the engine's block holds the wall among its other rows: take the whole block from the engine
            print(f"rmfit: Park's J > 0 wall is taken from the engine's rows (with its other model rows)", flush=True)
            self.park = None
        self.scat = ModelRows(n)
        print(f"rmfit: {n} model row(s) from the engine's objective (model data such as measured scattering "
              f"lengths, a prior covariance, decay priors, abundance sums"
              + (", Park's wall" if self.flags.get("use_park") else "") + ")", flush=True)

    def _scat_rows(self, x, jac=False):
        """The engine's model rows (and their Jacobian rows) at ``x``."""
        x = np.asarray(x, float)
        w, m_ = self._park_handle()          # a worker holding the whole R-matrix block
        r, Jl = w.call("model_rows", (x[m_], bool(jac)))
        r = np.asarray(r, float)
        if r.size != self.scat.n:
            raise RuntimeError(f"the engine gave {r.size} model rows, {self.scat.n} at start-up")
        if not jac:
            return r, None
        J = np.zeros((r.size, x.size))
        J[:, np.asarray(m_, int)] = np.asarray(Jl, float)
        return r, J

    def scat_value(self, r) -> float:
        """The model rows' share of the objective (they sit just before the Park wall rows)."""
        if self.scat is None:
            return 0.0
        r = np.asarray(r, float)
        end = r.size - (len(self.park.labels) if self.park is not None else 0)
        return float(np.sum(r[end - len(self.scat.parts):end] ** 2))

    def extra_penalty_value(self, r) -> float:
        """Rows the engine keeps inside its scalar chi-squared but rmfit counts in ``penalty``: Park's J > 0 wall
        and the measured scattering lengths."""
        return self.park_value(r) + self.scat_value(r)

    def rows_objective(self, r) -> float:
        """The objective from a (weighted) residual vector as ``residuals`` returns it."""
        r = np.asarray(r, float)
        data, rest = r[:self.npts], r[self.npts:]
        if self.spec.loss == "chi2":
            v = float(np.sum(data ** 2))
        elif self.spec.loss == "sivia":
            v = float(np.sum(rho_sivia(data ** 2)[0])) if self.row_c is None else \
                float(np.sum(rho_sivia_scaled(data ** 2, self.row_c)[0]))
        else:
            raise ValueError(f"unknown loss {self.spec.loss!r}")
        return v + float(np.sum(rest ** 2))

    def rows_chi2(self, r) -> float:
        """The plain data chi-squared from a (weighted) residual vector."""
        if self.cov is not None and self.row_sw is not None:
            return float(np.sum(self._cov_memo_get(r)[0] ** 2))
        return float(np.sum(self._plain_data_rows(r) ** 2))

    def point_chi2(self, r) -> np.ndarray:
        """Per data row, its share of the plain chi-squared (they add up to it): ``R_i^2``, or on
        the rows of a data covariance ``r_i (M^-1 r)_i`` (the engine's partition, can be negative)."""
        if self.cov is None:
            return self._plain_data_rows(r) ** 2
        return self._cov_memo_get(r)[1]

    def point_residuals(self, r) -> np.ndarray:
        """Per data row, the plain normalized residual (fit - n d)/(n sigma) -- before any weight and,
        on the rows of a data covariance, before whitening (the correlations ignored)."""
        if self.cov is None:
            return self._plain_data_rows(r)
        return self._cov_memo_get(r)[2]

    # -- data covariance (rmfit.covariance) --------------------------------------------------

    def _cov_norm_rows(self, x):
        """``(norms, cols)``: every data row's data normalization at ``x`` and the free column it is (-1: fixed,
        the .azr value) -- a normgroup member uses its leader's (rmfit.covariance: a fixed reference's relvector
        terms follow the normalization)."""
        if getattr(self, "_norm_rows", None) is None:
            groups = getattr(self, "norm_groups", None)
            if groups is None:
                groups = NormGroups.from_azr(self.azr)
            pos = {k: j for j, k in enumerate(self.keys)}
            col = np.full(self.npts, -1, int)
            nom = np.ones(self.npts)
            for key, (a, b) in self.row_ranges:
                lead = groups.owner(int(key))
                d = self.segments.by_key(lead)
                j = pos.get(ParamKey("norm", seg=lead)) if d.vary_norm else None
                col[a:b] = -1 if j is None else j
                nom[a:b] = float(d.norm)
            self._norm_rows = (col, nom)
        col, nom = self._norm_rows
        n = nom.copy()
        m = col >= 0
        n[m] = np.asarray(x, float)[col[m]]
        return n, col

    def free_norms_on_relvector(self) -> Dict[str, List[int]]:
        """``{covariance file basename: [segment keys]}``: the segments with relvector-covariance rows whose data
        normalization is a FREE parameter (their own or their normgroup leader's) -- the configuration the t0 method
        excludes (ObjectiveSpec.cov_free_norm_ok)."""
        if self.cov is None or not self.cov.relative:
            return {}
        _, cols = self._cov_norm_rows(np.zeros(len(self.keys)))
        out = {}
        for g in self.cov.groups:
            if not g.relative:
                continue
            mine = np.zeros(self.npts, bool)
            mine[g.pos] = True
            segs = sorted({int(k) for k, (a, b) in self.row_ranges if np.any(mine[a:b] & (cols[a:b] >= 0))})
            if segs:
                out[os.path.basename(g.path)] = segs
        return out

    def _cov_rows(self, r, J, x):
        """Plain data rows -> the rows the objective sees: row weights, then whitening."""
        norms, cols = self._cov_norm_rows(x)
        w, Jw, w0, contrib = self.cov.apply(r, J, self.row_sw, norms=norms, norm_cols=cols)
        memo = self.__dict__.setdefault("_cov_memo", OrderedDict())
        memo[_rows_key(w)] = (w0, contrib, np.array(r, float, copy=True))
        while len(memo) > 16:
            memo.popitem(last=False)
        return w, Jw

    def _cov_memo_get(self, r):
        memo = self.__dict__.get("_cov_memo") or {}
        got = memo.get(_rows_key(np.asarray(r, float)[:self.npts]))
        if got is None:
            raise RuntimeError("covariance: these residual rows did not come from this evaluator's residuals()")
        return got

    def refreeze_covariance(self) -> int:
        """Forget the model frozen for the ``relvector`` terms; the next evaluation freezes it again
        (the engine's ``refreeze_covariance``).  Returns the number of groups with relvector terms."""
        if self.cov is None:
            return 0
        self.cov.refreeze()
        return sum(1 for g in self.cov.groups if g.relative)

    def _plain_data_rows(self, r):
        data = np.asarray(r, float)[:self.npts]
        return data / self.row_sw if self.row_sw is not None else data

    @property
    def ls_loss(self):
        """``loss`` argument for ``scipy.optimize.least_squares`` (None = linear)."""
        return make_ls_loss(self.npts, self.row_c) if self.spec.loss == "sivia" else None

    def score(self, x) -> Score:
        t0 = time.time()
        r, _ = self.residuals(x, jac=False)
        if self.cov is not None:
            chi2 = self.rows_chi2(r)
            seg = self._segment_chi2_from_contrib(self.point_chi2(r))
        else:
            data = self._plain_data_rows(r)
            chi2 = float(np.sum(data ** 2))
            seg = self._segment_chi2_from_rows(data)
        pen = self.penalties.value(x) + self.extra_penalty_value(r)
        pri = self.prior.value(x) if self.prior is not None else 0.0
        obj = self.rows_objective(r) if not self.spec.plain else chi2 + pen + pri
        c2s = (self.cov.chi2_self(self._cov_memo_get(r)[2]) if self.cov is not None and self.cov.has_reference
               else None)
        return Score(obj, chi2, pen, pri, self.npts, seg, time.time() - t0, objective_kind=self.spec.label,
                     chi2_self=c2s)

    def _segment_chi2_from_rows(self, data_rows) -> Dict[int, float]:
        out: Dict[int, float] = {}
        for key, (a, b) in self.row_ranges:
            out[key] = out.get(key, 0.0) + float(np.sum(data_rows[a:b] ** 2))
        return out

    def _segment_chi2_from_contrib(self, contrib) -> Dict[int, float]:
        out: Dict[int, float] = {}
        for key, (a, b) in self.row_ranges:
            out[key] = out.get(key, 0.0) + float(np.sum(contrib[a:b]))
        return out

    def wigner_amplitudes(self) -> Dict[int, float]:
        """``{free position: gamma_W}`` in MeV^1/2 (the bound on |rwa|)."""
        return {c: float(np.sqrt(w)) for c, w in self.wigner_limits.items()}

    def wigner_widths(self, x) -> Dict[int, float]:
        """``{free position: Gamma_W}`` in eV -- the physical width of each free
        channel at the Wigner limit, from the engine (R-matrix block only)."""
        raise NotImplementedError

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()



def _skip_columns(cols, local_to_global, nlocal):
    """Local free-column indices a worker need not finite-difference: those whose global column is not in ``cols``
    (indices or a bool mask over the global free vector).  ``None`` when every column is wanted."""
    if cols is None:
        return None
    c = np.asarray(cols)
    want = set(np.where(c)[0].tolist()) if c.dtype == bool else set(int(v) for v in c.ravel())
    g = range(nlocal) if local_to_global is None else local_to_global
    return [j for j, gj in enumerate(g) if int(gj) not in want]

class SingleSession(Evaluator):
    """The full model in one worker process (authoritative)."""

    def __init__(self, azr, cwd=None, flags=None, name="single", objective=None):
        self.azr = os.path.abspath(azr)
        self.cwd = cwd or os.path.dirname(self.azr)
        self.flags = dict(flags or {})
        self.name = name
        # a model with data covariances runs WITHOUT the keyword (plain rows, whitened here: see
        # rmfit.covariance); normalization groups stay with the engine
        run_azr = self.azr
        if has_field(self.azr, "covariance"):
            run_azr = without_fields(self.azr, os.path.splitext(self.azr)[0] + ".nocov.azr", ("covariance",))
        self.run_azr = run_azr
        self.worker = _Worker(run_azr, self.cwd, self.flags, None, name=name).start()
        info = self.worker.info
        self.keys = info["keys"]
        self.x0 = np.asarray(info["x0"], float)
        self.n_rmatrix = info["n_rmatrix"]
        counts = info["npts"]
        edges = np.cumsum([0] + counts)
        self.row_ranges = [(k, (int(edges[i]), int(edges[i + 1])))
                           for i, k in enumerate(info["seg_keys"])]
        self.npts = int(edges[-1])
        self.nuisances = list(info.get("nuisances", []))
        self.segments = base_segments(self.azr).with_nuisances(self.nuisances)
        self.penalties = Penalties.from_nuisances(self.keys, self.nuisances)
        # the engine's rows of the fitted nuisances by free column (prior, limits), for the bounds
        # (local.default_bounds)
        self.nuisance_of_col = {int(r["free_index"]): r for r in self.nuisances
                                if not r["fixed"] and r.get("free_index") is not None}
        self.wigner_limits = dict(info["wigner"])
        self.level_energy = info["level_energy"]
        self.fixed_amplitudes = list(info.get("fixed_amplitudes", []))
        self.prior = None
        self.park = (ParkWall.build(self.keys, self.n_rmatrix, info.get("park_levels"))
                     if self.flags.get("use_park") else None)
        self._init_scat(self.worker, np.arange(len(self.keys)), self.azr)
        self.cov = None
        self._cov_check = False
        if run_azr != self.azr:
            blocks = [RowBlock(a, b, int(key), -np.inf, np.inf) for key, (a, b) in self.row_ranges]
            self.cov = Covariance.build(self.azr, self.cwd, self.segments, blocks, self.npts)
            if self.cov is not None:
                d = self.data()
                with np.errstate(divide="ignore", invalid="ignore"):
                    q = np.concatenate([np.asarray(c, float) / np.asarray(e, float)
                                        for c, e in zip(d["cross"], d["cross_err"])])
                self.cov.check_ratios(q)
                # at the first evaluation, compare with an engine session on the model as written
                self._cov_check = os.environ.get("RMFIT_COV_CHECK", "1") != "0"
        self._init_objective(objective)

    def _check_covariance(self, x, r):
        """The engine's chi-squared of the model as written (covariance included, the relvector
        model frozen at this same first evaluation) must equal the whitened rows' (rtol 1e-8)."""
        from .keys import take
        # the engine freezes the relvector model at this evaluation; with a fixed reference, rmfit's objective uses
        # the reference instead, so compare the engine with the rows frozen at this point (chi2_self)
        chi2 = (self.cov.chi2_self(self._cov_memo_get(r)[2]) if self.cov.has_reference
                else self.rows_chi2(r)) + self.extra_penalty_value(r)   # the engine's chi2 holds the J > 0 wall
                                                                        # and the scattering lengths
        chk = _Worker(self.azr, self.cwd, self.flags, None, name=f"{self.name}_covcheck").start()
        try:
            ck = chk.info["keys"]
            xe = take(np.asarray(x, float), self.keys, ck, np.asarray(chk.info["x0"], float))
            eng = float(chk.call("eval", xe))
        finally:
            chk.stop()
        if not np.isclose(eng, chi2, rtol=1e-8, atol=1e-10):
            raise RuntimeError(f"covariance: the engine's chi2 of {os.path.basename(self.azr)} is {eng:.10g}, rmfit's "
                               f"whitened rows give {chi2:.10g} (rel {abs(eng - chi2) / max(abs(eng), 1e-300):.2e}). "
                               f"A pyazr build without data covariances (before AZURE2 4b71177) ignores the keyword; "
                               f"otherwise rmfit's row mapping or whitening differs from the engine's.")
        print(f"rmfit covariance: {self.cov.summary()}; chi2 {chi2:.6f} = engine {eng:.6f} at the first evaluation"
              + (f" (frozen there; the objective uses the fixed reference {os.path.basename(self.spec.cov_reference)})"
                 if self.cov.has_reference else ""), flush=True)

    def _call(self, cmd, payload):
        try:
            return self.worker.call(cmd, payload)
        except EvalError:
            self.worker.restart()
            raise

    def objective(self, x):
        x = np.asarray(x, float)
        if not self.spec.plain or self.cov is not None or self.park is not None or self.scat is not None:
            r, _ = self.residuals(x, jac=False)
            return self.rows_objective(r)
        chi2 = self._call("eval", x)
        if not np.isfinite(chi2):
            raise EvalError("non-finite chi-squared")
        v = chi2 + self.penalties.value(x)
        if self.prior is not None:
            v += self.prior.value(x)
        return float(v)

    def residuals(self, x, jac=False, cols=None):
        x = np.asarray(x, float)
        skip = _skip_columns(cols, None, x.size) if jac else None
        r, J = self._call("resid", (x, True, skip) if skip else (x, bool(jac)))
        if self.cov is not None:
            if not np.all(np.isfinite(r)):
                raise EvalError("non-finite residuals")
            r, J = self._cov_rows(r, J, x)
        parts = [r, self.penalties.residuals(x)]
        Jparts = [J, self.penalties.jacobian(x.size)] if jac else None
        if self.prior is not None:
            parts.append(self.prior.residuals(x))
            if jac:
                Jparts.append(self.prior.jacobian(x.size))
        if self.scat is not None:
            rs, Js = self._scat_rows(x, jac)
            parts.append(rs)
            if jac:
                Jparts.append(Js)
        if self.park is not None:
            rp, Jp = self._park_rows(x, jac)
            parts.append(rp)
            if jac:
                Jparts.append(Jp)
        r = np.concatenate(parts)
        if not np.all(np.isfinite(r)):
            raise EvalError("non-finite residuals")
        if self.cov is not None:
            if self._cov_check:
                self._cov_check = False
                self._check_covariance(x, r)
            return r, (np.vstack(Jparts) if jac else None)
        return self._weight_rows(r, np.vstack(Jparts) if jac else None)

    def _park_handle(self):
        return self.worker, np.arange(len(self.keys))

    def transform(self, x):
        return self._call("transform", np.asarray(x, float))

    def wigner_widths(self, x):
        out = self._call("wigner_widths", np.asarray(x, float))
        return {int(k): float(v) for k, v in out.items() if int(k) < self.n_rmatrix}

    def write_output_files(self, x):
        return self._call("write_output", np.asarray(x, float))

    def data(self):
        return self._call("data", None)

    def cross(self, x):
        return self._call("cross", np.asarray(x, float))

    def close(self):
        self.worker.stop()


class Sharded(Evaluator):
    """The data split across ``nshards`` worker processes."""

    def __init__(self, base_azr, workdir, nshards=24, flags=None, weights=None,
                 no_split=(), reuse=True, objective=None, balance=None):
        self.base_azr = os.path.abspath(base_azr)
        self.workdir = os.path.abspath(workdir)
        self.flags = dict(flags or {})
        import pyazr
        if pyazr.Project.open(self.base_azr).config().get("priorCovariance"):
            # a prior covariance couples parameters of several data sets; a shard holds a subset of them and its
            # file drops the key, so the term would be missing here and present in a single session
            raise NotImplementedError(f"{os.path.basename(self.base_azr)}: a prior covariance (priorCovariance of "
                                      f"<config>) is supported in SingleSession only")
        self.segments = base_segments(self.base_azr)
        # shard files carry neither `covariance` nor `normgroup` (rmfit.model.write_shard): a follower's
        # pieces run at its own fixed norm and are rescaled here, covariance rows are whitened here
        self.norm_groups = NormGroups.from_azr(self.base_azr)
        # 2026-10-26 Jacobian-cost balancing (the evaluator: "go ahead with the Jacobian-cost shard balancing"): measure
        # each shard's residual and resid+Jacobian time once, fit R (alpha + beta k) (k = finite-differenced shift
        # columns), repack, cache the model in <workdir>/jac_balance.json for later evaluators of the same segments.
        # RMFIT_JAC_BALANCE=0 (or balance=False) keeps the objective-cost packing.
        if balance is None:
            balance = os.environ.get("RMFIT_JAC_BALANCE", "1") != "0"
        self._plan_args = (nshards, weights, tuple(no_split))
        self.balance_info = None
        cached = self._balance_cache_load() if balance else None
        self.shards: List[ShardSpec] = build_shards(self.base_azr, nshards, self.workdir,
                                                    weights=weights, no_split=no_split,
                                                    balance=cached["model"] if cached else None)
        self._start_workers()
        if balance and cached is None and len(self.shards) > 1:
            self._rebalance()
        elif cached is not None:
            self.balance_info = dict(cached, source="cache")
        self.prior = None
        self.park = (ParkWall.build(self.keys, self.n_rmatrix, self.workers[0].info.get("park_levels"))
                     if self.flags.get("use_park") else None)
        # the shards carry no <scatteringLengths> (model.write_shard); the rows come from the base model's lines
        self._init_scat(self.workers[0], self.maps[0], self.base_azr)
        self.timings: List[float] = []
        self._init_covariance()
        self._init_objective(objective)

    def _start_workers(self):
        self.workers: List[_Worker] = []
        for s in self.shards:
            self.workers.append(_Worker(s.azr, self.workdir, self.flags, s.seg_map, name=s.name))
        for w in self.workers:
            w.start()
        self._assemble()

    # -- Jacobian-cost balancing -------------------------------------------------

    def _balance_signature(self) -> str:
        # every record of the project that is not a level: the data sets and what they name (targets, beams, ...)
        import pyazr
        records = [r.text for r in pyazr.Project.open(self.base_azr).records().values()]
        nshards, weights, no_split = self._plan_args
        env = {k: os.environ.get(k, "") for k in ("RMFIT_SPLIT_SUM", "RMFIT_CONV_PPS", "RMFIT_BALANCE_PIECE")}
        blob = json.dumps(dict(records=records,
                               nshards=nshards, weights=weights, no_split=list(no_split), env=env), sort_keys=True)
        return hashlib.md5(blob.encode()).hexdigest()

    def _balance_cache_path(self):
        return os.path.join(self.workdir, "jac_balance.json")

    def _balance_cache_load(self):
        try:
            d = json.load(open(self._balance_cache_path()))
        except (OSError, ValueError):
            return None
        return d if d.get("signature") == self._balance_signature() and d.get("model") else None

    def _measure_shards(self, x):
        segs = self.segments
        rows = []
        for w, m_, s in zip(self.workers, self.maps, self.shards):
            t0 = time.time(); w.call("resid", (x[m_], False)); tr = time.time() - t0
            t0 = time.time(); w.call("resid", (x[m_], True)); tj = time.time() - t0
            rows.append(dict(name=s.name, plain=tr, jac=tj, cost=float(sum(p.cost for p in s.pieces)),
                             k=shard_shift_count(self.base_azr, s.pieces, segs),
                             pieces=[(int(p.base_key), float(p.cost)) for p in s.pieces]))
        return rows

    def _rebalance(self):
        """Measure, fit the cost model, repack if it predicts a >= 15 % shorter Jacobian; up to
        RMFIT_BALANCE_ROUNDS (2) rounds, each correcting the per-segment costs by the measured / predicted residual
        time of the shards holding them; cache the model."""
        x = np.asarray(self.x0, float)
        t_start = time.time()
        nshards, weights, no_split = self._plan_args
        gain = float(os.environ.get("RMFIT_BALANCE_MIN_GAIN", "0.15"))   # repack only for a predicted gain >= 15 %
        nround = max(1, int(os.environ.get("RMFIT_BALANCE_ROUNDS", "2")))
        from scipy.optimize import nnls
        allrows, rounds = [], []
        factors: Dict[str, float] = {}
        model = None
        msgs = []
        for rnd in range(nround):
            try:
                rows = self._measure_shards(x)
            except EvalError as err:
                msgs.append(f"round {rnd + 1}: measurement failed at x0 ({err})")
                break
            allrows += rows
            Jmax = max(r["jac"] for r in rows)
            if model is not None:
                rounds[-1]["measured_max"] = Jmax          # what the previous round's plan gave
            R = np.array([r["plain"] for r in allrows])
            Jt = np.array([r["jac"] for r in allrows])
            k = np.array([r["k"] for r in allrows], float)
            if np.any(k > 0):
                (alpha, beta), _ = nnls(np.column_stack([R, k * R]), Jt)
            else:
                alpha, beta = float(np.sum(Jt * R) / max(np.sum(R * R), 1e-300)), 0.0
            alpha = max(float(alpha), 1.0)
            # per-segment correction: the cost-weighted mean of measured / predicted over the shards holding it
            num, den = {}, {}
            for r in rows:
                q = r["plain"] / max(r["cost"], 1e-300)
                for key, c in r["pieces"]:
                    num[key] = num.get(key, 0.0) + c * q
                    den[key] = den.get(key, 0.0) + c
            for key in num:
                if den[key] > 0:
                    factors[str(key)] = factors.get(str(key), 1.0) * num[key] / den[key]
            model = dict(alpha=float(alpha), beta=float(beta), factors=dict(factors),
                         default_factor=float(np.median(list(factors.values()))) if factors else 1.0)
            _, loads = plan_balanced(self.base_azr, nshards, weights, no_split, True, model)
            pred = float(np.max(loads))
            rec = dict(round=rnd + 1, measured=rows, measured_max_before=Jmax, alpha=float(alpha), beta=float(beta),
                       predicted=[float(t) for t in loads], predicted_max=pred)
            rounds.append(rec)
            if not pred < (1.0 - gain) * Jmax:
                rec["repacked"] = False
                msgs.append(f"round {rnd + 1}: slowest shard {Jmax:.3g} s, plan predicts {pred:.3g} s -- kept")
                break
            self.close()
            self.shards = build_shards(self.base_azr, nshards, self.workdir, weights=weights, no_split=no_split,
                                       balance=model)
            self._start_workers()
            rec["repacked"] = True
            msgs.append(f"round {rnd + 1}: slowest shard {Jmax:.3g} s (alpha {alpha:.2f}, beta {beta:.2f}), "
                        f"repacked for {pred:.3g} s")
        info = dict(signature=self._balance_signature(), model=model, rounds=rounds,
                    repacked=any(r.get("repacked") for r in rounds), created=time.strftime("%Y-%m-%d %H:%M"))
        if model is not None and info["repacked"]:
            t0 = time.time()
            self._broadcast("resid", lambda m_: (x[m_], True))
            info["after_parallel_jac"] = time.time() - t0
            msgs.append(f"a parallel resid+Jacobian now takes {info['after_parallel_jac']:.3g} s")
        info["overhead_s"] = time.time() - t_start
        self.balance_info = info
        if model is not None:
            try:
                tmp = self._balance_cache_path() + ".tmp"
                json.dump(info, open(tmp, "w"), indent=1)
                os.replace(tmp, self._balance_cache_path())
            except OSError:
                pass
        print("rmfit shards (Jacobian balancing): " + "; ".join(msgs) + f" ({info['overhead_s']:.0f} s)", flush=True)

    def _norm_owner(self, base_key):
        """``(column of the norm segment base_key uses or None, its nominal value)``: the
        group leader's for a normgroup follower."""
        lead = self.norm_groups.owner(int(base_key))
        d = self.segments.by_key(lead)
        pos = {k: j for j, k in enumerate(self.keys)}
        return (pos.get(ParamKey("norm", seg=lead)) if d.vary_norm else None), float(d.norm)

    def _init_covariance(self):
        """Follower rescaling data (``self.follow``) and the covariance groups (``self.cov``)."""
        self.follow = {}
        need_data = bool(self.norm_groups) or has_keyword(self.base_azr, "covariance")
        parts = self._broadcast("data", lambda m_: None) if need_data else None
        q_all = np.full(self.npts, np.nan)
        for si, entries in enumerate(self.layout):
            for ei, (p, (sa, sb), (ga, gb), ncol, nominal) in enumerate(entries):
                if parts is None:
                    continue
                d = parts[si]
                with np.errstate(divide="ignore", invalid="ignore"):
                    q = np.asarray(d["cross"][ei], float) / np.asarray(d["cross_err"][ei], float)
                q_all[ga:gb] = q
                if p.base_key in self.norm_groups.leader:
                    col, nom = self._norm_owner(p.base_key)
                    n0 = float(p.norm0 if p.norm0 is not None else 1.0)
                    self.follow[(si, ei)] = (col, nom, n0, np.where(np.isfinite(q), q, 0.0),
                                             np.asarray(d["cross_err"][ei], float) != 0)
        self.cov = None
        if has_keyword(self.base_azr, "covariance"):
            blocks = []
            for entries in self.layout:
                for p, _, (ga, gb), _, _ in entries:
                    blocks.append(RowBlock(ga, gb, int(p.base_key), float(p.e_lo), float(p.e_hi), "file"))
            self.cov = Covariance.build(self.base_azr, os.path.dirname(self.base_azr), self.segments,
                                        blocks, self.npts)
            if self.cov is not None:
                self.cov.check_ratios(q_all)

    # -- layout -----------------------------------------------------------------

    def _assemble(self):
        """Global keys = the (identical) R-matrix block of every shard, then the
        union of norm/shift keys ordered by base segment key (plus the norms of
        convolved segments, which no shard carries)."""
        infos = [w.info for w in self.workers]
        rblocks = [i["keys"][:i["n_rmatrix"]] for i in infos]
        for rb in rblocks[1:]:
            if rb != rblocks[0]:
                raise RuntimeError("shards disagree on the R-matrix parameter block")
        self.n_rmatrix = len(rblocks[0])
        tailset = {k for i in infos for k in i["keys"][i["n_rmatrix"]:]}
        tail = sorted(tailset, key=lambda k: (k.kind != "norm", k.seg, _TAIL_ORDER.get(k.kind, 0), k.term or 0))
        self.keys = list(rblocks[0]) + tail
        nfree = len(self.keys)
        self.maps = [index_map(i["keys"], self.keys) for i in infos]
        for m_, i in zip(self.maps, infos):
            if np.any(m_ < 0):
                raise RuntimeError("a shard key is missing from the global keys")
        # x0: R-block from shard 0, norms/shifts from whichever shard has them
        x0 = np.zeros(nfree)
        for m_, i in zip(self.maps, infos):
            x0[m_] = i["x0"]
        self.x0 = x0
        # rows: shard order, piece order inside a shard.  self.layout[s] is a
        # list of (piece, shard_slice, global_slice, norm_col or None, nominal)
        pos = {k: j for j, k in enumerate(self.keys)}
        self.row_ranges = []
        self.shard_rows = []
        self.layout = []
        off = 0
        for s, i in zip(self.shards, infos):
            a = off
            soff = 0
            entries = []
            if len(i["npts"]) != len(s.pieces):
                raise RuntimeError(f"{s.name}: engine reports {len(i['npts'])} segments for {len(s.pieces)} pieces")
            for p, n in zip(s.pieces, i["npts"]):
                if n != p.npts:
                    raise RuntimeError(f"{s.name}: piece {p.base_key} [{p.e_lo},{p.e_hi}] loaded {n} points, planned {p.npts}")
                entries.append((p, (soff, soff + n), (off, off + n), None, None))
                self.row_ranges.append((int(p.base_key), (off, off + n)))
                off += n
                soff += n
            self.layout.append(entries)
            self.shard_rows.append((a, off))
        self.npts = off
        # the nuisances of every shard on the global columns and the base model's data-set keys: the priors
        # (Penalties) and the values of the base table's nuisance fields come from the shards' parameter rows
        self.nuisances = []
        for s, m_, i in zip(self.shards, self.maps, infos):
            smap = s.seg_map
            for r in i.get("nuisances", []):
                r = dict(r)
                if r.get("free_index") is not None:
                    r["free_index"] = int(m_[int(r["free_index"])])
                if r.get("segment_key") is not None:
                    r["segment_key"] = int(smap[int(r["segment_key"])])
                self.nuisances.append(r)
        # a normalization-group follower is written into its shard with the fixed norm 1 (model.write_shard): its
        # row says nothing about the base model, where it follows the group's first member
        followers = set(self.norm_groups.leader)
        self.segments.with_nuisances([r for r in self.nuisances
                                      if not (r["kind"] == "norm" and r.get("segment_key") in followers)])
        self.penalties = Penalties.from_nuisances(self.keys, self.nuisances)
        self.nuisance_of_col = {int(r["free_index"]): r for r in self.nuisances
                                if not r["fixed"] and r.get("free_index") is not None}
        self.wigner_limits = dict(infos[0]["wigner"])
        self.level_energy = infos[0]["level_energy"]
        self.fixed_amplitudes = list(infos[0].get("fixed_amplitudes", []))

    # -- evaluation -------------------------------------------------------------

    def _broadcast(self, cmd, payload_for):
        """Send to every worker, collect in order, restarting a dead worker
        once and retrying its task."""
        for w, m_ in zip(self.workers, self.maps):
            w.send(cmd, payload_for(m_))
        out = []
        for w, m_ in zip(self.workers, self.maps):
            try:
                out.append(w.recv())
            except EvalError as err:
                w.restart()
                try:
                    out.append(w.call(cmd, payload_for(m_)))
                except EvalError as err2:
                    raise EvalError(f"{err}; retry after restart also failed: {err2}")
        return out

    def objective(self, x):
        x = np.asarray(x, float)
        t0 = time.time()
        if (not self.spec.plain or self.follow or self.cov is not None or self.park is not None
                or self.scat is not None):
            r, _ = self.residuals(x, jac=False)
            self.timings.append(time.time() - t0)
            return self.rows_objective(r)
        parts = self._broadcast("eval", lambda m_: x[m_])
        self.timings.append(time.time() - t0)
        chi2 = float(np.sum(parts))
        if not np.isfinite(chi2):
            raise EvalError("non-finite chi-squared")
        v = chi2 + self.penalties.value(x)
        if self.prior is not None:
            v += self.prior.value(x)
        return v

    def shard_times(self, x, repeat=1):
        """Per-shard wall time of one evaluation (for weight calibration)."""
        x = np.asarray(x, float)
        times = np.zeros(len(self.workers))
        for _ in range(repeat):
            for i, (w, m_) in enumerate(zip(self.workers, self.maps)):
                t0 = time.time()
                w.call("eval", x[m_])
                times[i] += time.time() - t0
        return times / repeat

    def residuals(self, x, jac=False, cols=None):
        x = np.asarray(x, float)
        nfree = x.size

        def payload(m_):
            skip = _skip_columns(cols, m_, len(m_)) if jac else None
            return (x[m_], True, skip) if skip else (x[m_], bool(jac))
        parts = self._broadcast("resid", payload)
        r = np.empty(self.npts)
        J = np.zeros((self.npts, nfree)) if jac else None
        for si, (entries, m_, (rs, Js)) in enumerate(zip(self.layout, self.maps, parts)):
            for ei, (p, (sa, sb), (ga, gb), ncol, nominal) in enumerate(entries):
                if p.kind == "data":
                    fw = self.follow.get((si, ei)) if self.follow else None
                    if fw is not None:
                        # a normgroup follower, evaluated at its own fixed norm n0: rescale to the
                        # leader's norm n, r = (n0/n)(r0 + q) - q with q = d/sigma per row
                        col, nom, n0, q, ok = fw
                        n = float(x[col]) if col is not None else nom
                        a = rs[sa:sb] + q
                        r[ga:gb] = np.where(ok, (n0 / n) * a - q, 0.0)
                        if jac:
                            for j, gj in enumerate(m_):
                                J[ga:gb, gj] += (n0 / n) * Js[sa:sb, j]
                            if col is not None:
                                J[ga:gb, col] += np.where(ok, -(n0 / (n * n)) * a, 0.0)
                        continue
                    r[ga:gb] = rs[sa:sb]
                    if jac:
                        for j, gj in enumerate(m_):
                            J[ga:gb, gj] += Js[sa:sb, j]
                    continue
        if self.cov is not None:
            if not np.all(np.isfinite(r)):
                raise EvalError("non-finite residuals")
            r, J = self._cov_rows(r, J, x)
        rows = [r, self.penalties.residuals(x)]
        Jrows = [J, self.penalties.jacobian(nfree)] if jac else None
        if self.prior is not None:
            rows.append(self.prior.residuals(x))
            if jac:
                Jrows.append(self.prior.jacobian(nfree))
        if self.scat is not None:
            rs, Js = self._scat_rows(x, jac)
            rows.append(rs)
            if jac:
                Jrows.append(Js)
        if self.park is not None:
            rp, Jp = self._park_rows(x, jac)
            rows.append(rp)
            if jac:
                Jrows.append(Jp)
        r = np.concatenate(rows)
        if not np.all(np.isfinite(r)):
            raise EvalError("non-finite residuals")
        if self.cov is not None:
            return r, (np.vstack(Jrows) if jac else None)
        return self._weight_rows(r, np.vstack(Jrows) if jac else None)

    def transform(self, x):
        x = np.asarray(x, float)
        # the R-matrix block is what transform reports; any shard can do it
        w, m_ = self.workers[0], self.maps[0]
        return w.call("transform", x[m_])

    def row_excitation(self) -> np.ndarray:
        out = np.full(self.npts, np.nan)
        parts = self._broadcast("data", lambda m_: None)
        for entries, d in zip(self.layout, parts):
            exc = d.get("excitation", [])
            for q, (p, (sa, sb), (ga, gb), ncol, nom) in enumerate(entries):
                if q < len(exc) and len(exc[q]) == gb - ga:
                    out[ga:gb] = np.asarray(exc[q], float)
        return out

    def _park_handle(self):
        return self.workers[0], self.maps[0]

    def wigner_widths(self, x):
        x = np.asarray(x, float)
        w, m_ = self.workers[0], self.maps[0]
        out = w.call("wigner_widths", x[m_])
        return {int(k): float(v) for k, v in out.items() if int(k) < self.n_rmatrix}

    def close(self):
        for w in self.workers:
            w.stop()


def nudge_zero_widths(x, keys, eps=1e-6):
    """A copy of ``x`` with every free width amplitude that is EXACTLY 0 set to ``eps`` (2026-09-30: AZURE2
    writes no parameter for a zero width, so a baked model re-read by the engine has one key fewer per zero
    and the bake check refuses it; 1e-6 MeV^1/2 is negligible).  Returns ``(x, [nudged key names])``."""
    x = np.array(x, float, copy=True)
    nudged = []
    for j, k in enumerate(keys):
        if k.kind == "width" and j < x.size and x[j] == 0.0:
            x[j] = eps
            nudged.append(str(k))
    return x, nudged


def bake(base_azr, x, keys, path, flags, workdir=None, output_dir=None,
         verify_rtol=1e-4, name="bake", objective=None, separate=True, separate_max_cost=10.0):
    """Write ``x`` (given in ``keys`` order) into a copy of ``base_azr`` at
    ``path`` plus a companion ``<path stem>.sav`` carrying every parameter
    (norms and shifts included), and verify the pair by re-evaluating it in a
    fresh session: the objective must agree to ``verify_rtol``.

    The files are written by ``pyazr.Session.save_fit`` (in a worker process); the check is rmfit's own, with
    rmfit's objective, every session in a process of its own.

    Returns ``(azr_path, sav_path, Score)``.
    """
    from .keys import take
    workdir = workdir or os.path.join(os.path.dirname(os.path.abspath(base_azr)), "variants")
    x, nudged = nudge_zero_widths(x, keys)
    if nudged:
        print(f"bake: {len(nudged)} free width(s) at exactly 0 nudged to 1e-6 so the written model keeps them as parameters: "
              + ", ".join(nudged[:6]) + (" ..." if len(nudged) > 6 else ""))
    src = single_variant(base_azr, f"{name}_src", workdir)
    with SingleSession(src, cwd=workdir, flags=flags, name=f"{name}_src", objective=objective) as ev:
        xs = take(np.asarray(x, float), keys, ev.keys, ev.x0)
        if len(set(keys) & set(ev.keys)) != len(ev.keys):
            missing = set(ev.keys) - set(keys)
            raise ValueError(f"bake: {len(missing)} free parameters of the base model are not in keys: "
                             + ", ".join(str(k) for k in list(missing)[:6]))
        want = ev.score(xs)
        from .enginepin import keeps_close_levels
        if separate and not keeps_close_levels():
            # 2026-10-10: on a build before dev e239a12 same-J^pi levels < 1 keV apart merge when the written file is read
            # back; move the broader member of any pair < 1.5 keV away to 2.5 keV first (local.separate_close_levels), if
            # that costs less than ``separate_max_cost`` in the objective.  Narrow pairs are refused there and fail below as
            # before.  2026-11-08 (the evaluator): a newer engine reads such a pair back as two levels, so the fit is written as it is.
            from .local import separate_close_levels
            xn, moved, refused = separate_close_levels(ev, xs)
            for r in refused:
                print(f"bake: same-J^pi pair NOT separated: {r}")
            if moved:
                trial = ev.score(xn); cost = trial.objective - want.objective
                desc = "; ".join(f"{m['level']} {m['e_from']:.6f} -> {m['e_to']:.6f} MeV ({m['shift_kev']:+.2f} keV, Gamma "
                                 f"{m['gamma_kev']:.4g} keV, was {m['gap_kev']:.3f} keV from {m['other']})" for m in moved)
                if cost > separate_max_cost:
                    raise ValueError(f"bake: separating near-degenerate levels would cost {cost:+.3f} > {separate_max_cost:g} "
                                     f"in the objective: {desc}")
                print(f"bake: separated near-degenerate same-J^pi level(s) before writing, objective {want.objective:.3f} -> "
                      f"{trial.objective:.3f} ({cost:+.3f}): {desc}")
                xs, want = xn, trial
        rmask = np.array([k.kind in ("energy", "width") for k in ev.keys])
        want_levels = ev.score(np.where(rmask, xs, ev.x0))   # norms/shifts at nominal
        src_keys = list(ev.keys)
    # The snapshot is Session.save_fit, in a process of its own on the base model as written (its covariance
    # fields included; the evaluator above ran a .nocov copy): every level's physical values into the rows they
    # come from, fixed parameters included (a fixed parameter is fixed in the fit's basis, so its observed value
    # moves with the other amplitudes of its level), the nuisances left at the file's values (norms="nominal": the
    # prior centres stay where they are), and the companion parameter file with every parameter and its basis tag.
    # pyazr's own re-read check is off: the check below does it in fresh processes with rmfit's objective.
    path = os.path.abspath(path)
    # save_fit rewrites a snapshot's relative paths for its directory but cannot rewrite those inside a data
    # covariance file, and refuses such a project in another directory: so the session that writes the snapshot
    # opens a copy of the base model in the snapshot's own directory (its output in the work directory)
    from .azr3 import write_variant
    from .model import link_data_dirs
    dest = os.path.dirname(path)
    os.makedirs(dest, exist_ok=True)
    link_data_dirs(base_azr, dest)
    save_out = os.path.join(workdir, f"output_{name}_save")
    os.makedirs(save_out, exist_ok=True)
    save_src, _ = write_variant(base_azr, os.path.join(dest, f".{name}_save_src.azr"), save_out)
    w = _Worker(save_src, dest, flags, None, name=f"{name}_save").start()
    try:
        xw = take(xs, src_keys, w.info["keys"], np.asarray(w.info["x0"], float))
        azr_written, sav_written = w.call("save_fit", (path, xw, dict(param_sav=True, verify=False, norms="nominal")))
    finally:
        w.stop()
        if os.path.exists(save_src):
            os.remove(save_src)
    sav = os.path.splitext(path)[0] + ".sav"
    if sav_written and os.path.abspath(sav_written) != sav:
        os.replace(sav_written, sav)
    # the snapshot carries the source variant's output directory: give it the caller's, else the base model's
    import pyazr
    out = pyazr.Project.open(path)
    out.set_config(output_dir=(output_dir if output_dir is not None
                               else pyazr.Project.open(os.path.abspath(base_azr)).config().get("output_dir", "output/")))
    out.save()
    chk = single_variant(path, f"{name}_chk", workdir)
    with SingleSession(chk, cwd=workdir, flags=flags, name=f"{name}_chk", objective=objective) as ev2:
        x2 = take(xs, ev.keys, ev2.keys, ev2.x0)
        if set(ev2.keys) != set(ev.keys):
            # 2026-09-16: this fires when the written model is not the same model.  The
            # case seen: the polish had drifted two same-J^pi levels to within 0.66 keV
            # (under a keV apart), and the engine reading the written file treats them as
            # ONE level, so the group's labels shift (9/2+#6 -> #5) and 8 keys disappear.
            # Delete the files as the rtol branch does -- leaving them on disk left a
            # directory holding a `<name>.azr` that looks like an export and is not.
            missing = sorted(map(str, set(ev.keys) - set(ev2.keys)))
            extra = sorted(map(str, set(ev2.keys) - set(ev.keys)))
            for f in (path, sav):
                if os.path.exists(f):
                    os.remove(f)
            raise RuntimeError(
                "the baked model has a different free-parameter set than its source (files removed): "
                f"{len(missing)} key(s) lost, e.g. {missing[:4]}; {len(extra)} gained, e.g. {extra[:4]}. "
                "Two levels of one J^pi read as one are the usual cause: at the identical energy, with one level ID, "
                "or (builds before dev e239a12) within 1 keV.")
        got = ev2.score(x2)
        # the baked <levels> block on its own (norms/shifts nominal) must
        # reproduce the source evaluated the same way: that is what a CLI run
        # with a blank external parameter file sees
        levels_only = ev2.score(ev2.x0)
    d_obj = abs(got.objective - want.objective) / abs(want.objective)
    d_lev = abs(levels_only.chi2 - want_levels.chi2) / abs(want_levels.chi2)
    if d_obj > verify_rtol or d_lev > verify_rtol:
        for f in (path, sav):
            if os.path.exists(f):
                os.remove(f)
        raise RuntimeError(
            f"bake verification failed (rtol {verify_rtol:g}): source objective {want.objective:.6f}, "
            f"baked {got.objective:.6f} (rel {d_obj:.2e}); <levels>-only chi2 {levels_only.chi2:.6f} vs "
            f"{want_levels.chi2:.6f} (rel {d_lev:.2e}); files removed")
    got.roundtrip = dict(objective_rel=d_obj, levels_rel=d_lev, levels_only_chi2=levels_only.chi2)
    return path, sav, got


def verify(base_azr, x, flags, workdir=None, name="verify", objective=None) -> Tuple[Score, List[ParamKey]]:
    """Authoritative ``(score, keys)`` of ``x`` from a fresh single session on
    the full model, in its own output directory.  ``x`` must be in the base
    model's own free-vector order (use :func:`rmfit.keys.take` to map a
    sharded vector onto ``keys`` first if the orders differ)."""
    workdir = workdir or os.path.join(os.path.dirname(os.path.abspath(base_azr)), "variants")
    azr = single_variant(base_azr, name, workdir)
    with SingleSession(azr, cwd=workdir, flags=flags, name=name, objective=objective) as ev:
        return ev.score(x), ev.keys
