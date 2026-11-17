"""Joint fit of several AZURE2 projects to data they share (stage 2 of the
natural-abundance / compound-target work, the evaluator 2026-10-12).

A natural target's cross section is an abundance-weighted sum over isotopes,
each with its own compound nucleus and therefore its own ``.azr``.  AZURE2
holds one compound nucleus per process, so here every isotope ("member") runs
in its own worker process (:class:`rmfit.engine._Worker`) and the combination
is done in Python:

    F(E) = sum_k w_k f_k(E)                  (shared segment)
    r    = (F - y n) / (e n)                 (the engine's residual convention)

with ``f_k`` the member's bare model on the shared data points, ``w_k`` its
weight (abundance, or atoms per formula unit), ``y``/``e`` the data and ``n``
the ONE normalization of the shared data set.  The Jacobian follows:
``dr/dtheta_k = w_k G_k / (e n)`` with ``G_k = df_k/dtheta_k`` from the member's
own residual Jacobian, ``dr/dn = -F / (e n^2)``, ``dr/dw_k = f_k / (e n)``.
Segments a member has to itself ("private") pass through unchanged.

The global free vector is every member's R-matrix block (member-tagged keys,
level labels ``"22Ne:1/2+#3"``), then every private normalization / shift,
then one normalization per shared data set, then the free weights.  AZURE2's
normalization and shift penalties are added here as rows, as in
:class:`rmfit.engine.Sharded`; free weights get a Gaussian prior and,
optionally, a sum-to-one row (without one of those a free weight and the
shared normalization are degenerate).

Every member's ``.azr`` carries the shared data segment as an ordinary line
(same data file, same energy window, same observable; its own normalization
is pinned to 1 during the joint fit, its energy shift fixed at the file's
value).  :func:`freeze_and_exchange` writes every member back as a standalone
stage-1 file -- the fitted parameters, the shared normalization, and the other
members' contributions as ``mixture`` tables -- and checks that each
reproduces the joint fit on its own.

Driver::

    python3 -m rmfit.mixture <dir> init      # read mixture.json, start, check, score
    python3 -m rmfit.mixture <dir> polish    # local least-squares from state.json
    python3 -m rmfit.mixture <dir> freeze    # standalone files + tables + check
    python3 -m rmfit.mixture <dir> status
"""
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import AZURE2_ROOT, ARCHIVE_FLAGS
from .backend import require as _require
# pyazr 2 only (the evaluator 2026-11-16: "pyazr-2 only"): superseded since 2026-10-13 by several compound nuclei in one project,
# and the pyazr calls it needs do not exist on pyazr 3
_require(2, "rmfit.mixture (the joint fit of several single-nucleus projects; use one project with several compound nuclei)")
from .engine import _Worker, Evaluator, EvalError, ObjectiveSpec, Penalties
from .keys import ParamKey
from .model import single_variant, replace_block, set_output_dir_text, block_lines
from pyazr.datasets import SegmentSet

SHARED_GID_BASE = 900000       # global segment ids of shared data sets
MEMBER_GID_STRIDE = 10000      # private segment gid = (member index + 1) * stride + key


# ---------------------------------------------------------------------------
# the specification

@dataclass
class MemberSpec:
    name: str
    azr: str
    flags: dict = field(default_factory=lambda: dict(ARCHIVE_FLAGS))
    weight: float = 1.0
    weight_free: bool = False
    weight_sigma: float = 0.0          # Gaussian prior on a free weight (0 = none)


@dataclass
class SharedSpec:
    name: str
    segments: Dict[str, int]           # member name -> that member's <segmentsData> key


@dataclass
class MixtureSpec:
    members: List[MemberSpec]
    shared: List[SharedSpec]
    sum_to_one: bool = True
    sum_sigma: float = 1e-3
    objective: dict = field(default_factory=dict)
    path: str = ""

    @classmethod
    def load(cls, path):
        path = os.path.abspath(path)
        if os.path.isdir(path):
            path = os.path.join(path, "mixture.json")
        base = os.path.dirname(path)
        with open(path) as f:
            cfg = json.load(f)
        members = []
        for m in cfg["members"]:
            azr = m["azr"]
            if not os.path.isabs(azr):
                azr = os.path.normpath(os.path.join(base, azr))
            members.append(MemberSpec(name=str(m["name"]), azr=azr,
                                      flags=dict(m.get("flags", ARCHIVE_FLAGS)),
                                      weight=float(m.get("weight", 1.0)),
                                      weight_free=bool(m.get("weight_free", False)),
                                      weight_sigma=float(m.get("weight_sigma", 0.0))))
        shared = [SharedSpec(name=str(s["name"]), segments={str(k): int(v) for k, v in s["segments"].items()})
                  for s in cfg.get("shared", [])]
        w = cfg.get("weights", {})
        spec = cls(members, shared, sum_to_one=bool(w.get("sum_to_one", True)),
                   sum_sigma=float(w.get("sum_sigma", 1e-3)), objective=dict(cfg.get("objective", {})),
                   path=path)
        spec.validate()
        return spec

    @property
    def directory(self):
        return os.path.dirname(self.path) if self.path else os.getcwd()

    def member(self, name) -> MemberSpec:
        for m in self.members:
            if m.name == name:
                return m
        raise KeyError(name)

    def validate(self):
        names = [m.name for m in self.members]
        if len(set(names)) != len(names):
            raise ValueError(f"mixture: duplicate member names {names}")
        if not self.members:
            raise ValueError("mixture: no members")
        for m in self.members:
            if not os.path.exists(m.azr):
                raise FileNotFoundError(f"mixture: member {m.name}: {m.azr}")
            if m.weight < 0:
                raise ValueError(f"mixture: member {m.name}: negative weight")
        seen = set()
        for s in self.shared:
            if s.name in seen:
                raise ValueError(f"mixture: shared data set {s.name!r} listed twice")
            seen.add(s.name)
            if not s.segments:
                raise ValueError(f"mixture: shared data set {s.name!r} names no segment")
            for mname, key in s.segments.items():
                if mname not in names:
                    raise ValueError(f"mixture: shared {s.name!r} names unknown member {mname!r}")
                segs = SegmentSet.from_file(self.member(mname).azr)
                d = segs.by_key(key)
                if d is None:
                    raise ValueError(f"mixture: shared {s.name!r}: member {mname} has no segment {key}")
                if not d.active:
                    raise ValueError(f"mixture: shared {s.name!r}: member {mname} segment {key} is inactive")
                if d.operation == "ratio":
                    raise ValueError(f"mixture: shared {s.name!r}: member {mname} segment {key} is a ratio")
                if d.mixture_self is not None:
                    raise ValueError(f"mixture: shared {s.name!r}: member {mname} segment {key} already carries a "
                                     "mixture block; the joint fit replaces it -- remove it")
            # the same data file, window and observable everywhere
            ref = None
            for mname, key in s.segments.items():
                d = SegmentSet.from_file(self.member(mname).azr).by_key(key)
                sig = (os.path.realpath(os.path.join(os.path.dirname(self.member(mname).azr), d.data_file)),
                       round(d.energy_min, 9), round(d.energy_max, 9), round(d.angle_min, 6), round(d.angle_max, 6),
                       d.observable)
                if ref is None:
                    ref = (mname, sig)
                elif sig != ref[1]:
                    raise ValueError(f"mixture: shared {s.name!r}: members {ref[0]} and {mname} disagree on the "
                                     f"segment ({ref[1]} vs {sig})")
        free_w = [m for m in self.members if m.weight_free]
        if free_w and not self.sum_to_one and any(m.weight_sigma <= 0 for m in free_w):
            raise ValueError("mixture: a free weight without a prior and without sum_to_one is degenerate "
                             "with the shared normalization")


class _SegList:
    """What :class:`Evaluator` expects of ``segments``: iteration and ``by_key``."""

    def __init__(self, segs):
        self._segs = list(segs)
        self._by = {int(s.key): s for s in self._segs}

    def __iter__(self):
        return iter(self._segs)

    def __len__(self):
        return len(self._segs)

    def by_key(self, key):
        return self._by.get(int(key))


class MixturePenalties:
    """AZURE2's norm/shift penalty rows, the weight priors, and a sum-to-one row."""

    def __init__(self, base: Penalties, sum_cols=None, sum_const=0.0, sum_target=1.0, sum_sigma=1e-3):
        self.base = base
        self.sum_cols = None if sum_cols is None or len(sum_cols) == 0 else np.asarray(sum_cols, int)
        self.sum_const, self.sum_target, self.sum_sigma = float(sum_const), float(sum_target), float(sum_sigma)

    @property
    def cols(self):
        return self.base.cols

    @property
    def nominal(self):
        return self.base.nominal

    @property
    def sigma(self):
        return self.base.sigma

    def residuals(self, x):
        r = self.base.residuals(x)
        if self.sum_cols is not None:
            s = self.sum_const + float(np.sum(np.asarray(x, float)[self.sum_cols]))
            r = np.concatenate([r, [(s - self.sum_target) / self.sum_sigma]])
        return r

    def jacobian(self, nfree):
        J = self.base.jacobian(nfree)
        if self.sum_cols is not None:
            row = np.zeros((1, nfree))
            row[0, self.sum_cols] = 1.0 / self.sum_sigma
            J = np.vstack([J, row])
        return J

    def value(self, x):
        return float(np.sum(self.residuals(x) ** 2))


def _tag(k: ParamKey, member: str, seg_gid: Optional[int] = None) -> ParamKey:
    """A member's key as a global one: member-prefixed level label, member field, global segment id."""
    level = None if k.level is None else f"{member}:{k.level}"
    return replace(k, level=level, member=member, seg=(k.seg if seg_gid is None else seg_gid))


def member_gid(member_index: int, seg_key: int) -> int:
    return (member_index + 1) * MEMBER_GID_STRIDE + int(seg_key)


def shared_gid(group_index: int) -> int:
    return SHARED_GID_BASE + group_index


# ---------------------------------------------------------------------------
# the evaluator

class Mixture(Evaluator):
    """The joint objective of a :class:`MixtureSpec` (see the module docstring)."""

    def __init__(self, spec: MixtureSpec, workdir: Optional[str] = None, objective=None,
                 omp_threads: int = 1, check: bool = True, log=None):
        self.mix = spec if isinstance(spec, MixtureSpec) else MixtureSpec.load(spec)
        self.log = log or (lambda *a: None)
        self.workdir = os.path.abspath(workdir or os.path.join(self.mix.directory, "variants"))
        os.makedirs(self.workdir, exist_ok=True)
        self.members = list(self.mix.members)
        self.workers: List[_Worker] = []
        self.member_azr: Dict[str, str] = {}
        for m in self.members:
            vdir = os.path.join(self.workdir, m.name)
            azr = single_variant(m.azr, m.name, vdir)
            text = relocate_data_paths(open(azr).read(), os.path.dirname(m.azr), vdir)
            with open(azr, "w") as fh:
                fh.write(text)
            self.member_azr[m.name] = azr
            self.workers.append(_Worker(azr, vdir, m.flags, None, name=m.name, omp_threads=omp_threads))
        for w in self.workers:
            w.start()
        self.infos = [w.info for w in self.workers]
        self.member_data = [w.call("data", None) for w in self.workers]
        self.member_segments = [SegmentSet.from_file(m.azr) for m in self.members]
        self.prior = None
        self.timings: List[float] = []
        self._assemble()
        self._init_objective(objective if objective is not None else
                             (ObjectiveSpec(**self.mix.objective) if self.mix.objective else None))
        if check:
            self._self_check()

    # -- layout -----------------------------------------------------------------

    def _assemble(self):
        spec = self.mix
        K = len(self.members)
        names = [m.name for m in self.members]
        # shared segment key per (member index) for quick lookup
        self.shared_of_member: List[Dict[int, int]] = [{} for _ in range(K)]      # member seg key -> group index
        for g, s in enumerate(spec.shared):
            for mname, key in s.segments.items():
                self.shared_of_member[names.index(mname)][int(key)] = g

        keys: List[ParamKey] = []
        x0: List[float] = []
        self.maps: List[np.ndarray] = []
        self.pinned: List[np.ndarray] = []
        self.nR: List[int] = []
        self.level_energy: Dict[str, float] = {}
        self.wigner_limits: Dict[int, float] = {}
        self.fixed_amplitudes: List[ParamKey] = []
        # 1. R-matrix blocks
        rstart = []
        for k, (m, info) in enumerate(zip(self.members, self.infos)):
            nR = int(info["n_rmatrix"])
            self.nR.append(nR)
            rstart.append(len(keys))
            for j in range(nR):
                keys.append(_tag(info["keys"][j], m.name))
                x0.append(float(info["x0"][j]))
            for fi, w in info["wigner"]:
                if int(fi) < nR:
                    self.wigner_limits[rstart[k] + int(fi)] = float(w)
            for lvl, e in info["level_energy"].items():
                self.level_energy[f"{m.name}:{lvl}"] = float(e)
            for fk in info.get("fixed_amplitudes", []):
                self.fixed_amplitudes.append(_tag(fk, m.name))
        self.n_rmatrix = len(keys)
        # 2. private norms / shifts
        pos_private: List[Dict[int, int]] = []
        for k, (m, info) in enumerate(zip(self.members, self.infos)):
            pp = {}
            for j in range(self.nR[k], len(info["keys"])):
                key = info["keys"][j]
                if int(key.seg) in self.shared_of_member[k]:
                    continue
                keys.append(_tag(key, m.name, member_gid(k, key.seg)))
                x0.append(float(info["x0"][j]))
                pp[j] = len(keys) - 1
            pos_private.append(pp)
        # 3. one normalization per shared data set
        self.shared_norm_col: List[Optional[int]] = []
        self.shared_norm_nominal: List[float] = []
        self.shared_norm_sigma: List[float] = []
        for g, s in enumerate(spec.shared):
            vary, nominal, sigma = False, 1.0, 0.0
            first = True
            for mname, key in s.segments.items():
                d = self.member_segments[names.index(mname)].by_key(key)
                if first:
                    nominal = float(d.norm)
                    first = False
                if d.vary_norm:
                    vary = True
                    sigma = max(sigma, float(d.norm) / 100.0 * float(d.norm_error))
            if vary:
                keys.append(ParamKey("norm", seg=shared_gid(g)))
                x0.append(nominal)
                self.shared_norm_col.append(len(keys) - 1)
            else:
                self.shared_norm_col.append(None)
            self.shared_norm_nominal.append(nominal)
            self.shared_norm_sigma.append(sigma)
        # 4. free weights
        self.weight_col: List[Optional[int]] = []
        for m in self.members:
            if m.weight_free:
                keys.append(ParamKey("weight", level=m.name))
                x0.append(float(m.weight))
                self.weight_col.append(len(keys) - 1)
            else:
                self.weight_col.append(None)
        self.keys = keys
        self.x0 = np.asarray(x0, float)
        nfree = len(keys)
        # member free vectors: global position or pinned value
        for k, (m, info) in enumerate(zip(self.members, self.infos)):
            n_k = len(info["keys"])
            mp = np.full(n_k, -1, int)
            pin = np.zeros(n_k)
            for j in range(self.nR[k]):
                mp[j] = rstart[k] + j
            for j in range(self.nR[k], n_k):
                if j in pos_private[k]:
                    mp[j] = pos_private[k][j]
                else:
                    key = info["keys"][j]
                    pin[j] = 1.0 if key.kind == "norm" else float(info["x0"][j])
            self.maps.append(mp)
            self.pinned.append(pin)
        # rows: private segments per member, then the shared data sets
        self.row_ranges: List[Tuple[int, Tuple[int, int]]] = []
        self.member_private: List[List[Tuple[int, int, int, int]]] = []   # (seg_index, a_k, b_k, g_a)
        off = 0
        segs = []
        for k, (m, info) in enumerate(zip(self.members, self.infos)):
            edges = np.cumsum([0] + list(info["npts"]))
            entries = []
            for i, skey in enumerate(info["seg_keys"]):
                a, b = int(edges[i]), int(edges[i + 1])
                if int(skey) in self.shared_of_member[k]:
                    continue
                gid = member_gid(k, skey)
                entries.append((i, a, b, off))
                self.row_ranges.append((gid, (off, off + (b - a))))
                d = self.member_segments[k].by_key(skey)
                segs.append(replace(d, key=gid))
                off += b - a
            self.member_private.append(entries)
        self.shared_rows: List[Tuple[int, int]] = []
        self.shared_members: List[List[Tuple[int, int, int]]] = []       # per group: (member index, a_k, b_k)
        self.shared_y: List[np.ndarray] = []
        self.shared_e: List[np.ndarray] = []
        for g, s in enumerate(spec.shared):
            mems = []
            y = e = en = None
            ref_name = None
            for mname, key in s.segments.items():
                k = names.index(mname)
                info = self.infos[k]
                i = list(map(int, info["seg_keys"])).index(int(key))
                edges = np.cumsum([0] + list(info["npts"]))
                a, b = int(edges[i]), int(edges[i + 1])
                yk = np.asarray(self.member_data[k]["cross"][i], float)
                ek = np.asarray(self.member_data[k]["cross_err"][i], float)
                enk = np.asarray(self.member_data[k]["energies"][i], float)
                if y is None:
                    y, e, en, ref_name = yk, ek, enk, mname
                elif yk.shape != y.shape or not (np.allclose(yk, y) and np.allclose(ek, e)):
                    raise RuntimeError(f"mixture: shared {s.name!r}: members {ref_name} and {mname} loaded "
                                       f"different data ({yk.size} vs {y.size} points)")
                mems.append((k, a, b))
            n = y.size
            gid = shared_gid(g)
            self.shared_members.append(mems)
            self.shared_y.append(y)
            self.shared_e.append(e)
            self.shared_rows.append((off, off + n))
            self.row_ranges.append((gid, (off, off + n)))
            k0 = mems[0][0]
            d = self.member_segments[k0].by_key(s.segments[names[k0]])
            segs.append(replace(d, key=gid, norm=self.shared_norm_nominal[g]))
            off += n
        self.npts = off
        self.segments = _SegList(segs)
        # penalties: private norms/shifts, shared norms, weight priors, sum-to-one
        cols, nom, sig = [], [], []
        for j, key in enumerate(self.keys):
            if key.kind in ("norm", "shift") and key.member:
                d = self.segments.by_key(key.seg)
                if key.kind == "norm":
                    s_ = d.norm / 100.0 * d.norm_error
                    if s_ > 0:
                        cols.append(j); nom.append(d.norm); sig.append(s_)
                elif d.vary_shift and d.energy_shift_error:
                    cols.append(j); nom.append(d.energy_shift); sig.append(d.energy_shift_error)
            elif key.kind == "norm":
                g = key.seg - SHARED_GID_BASE
                if self.shared_norm_sigma[g] > 0:
                    cols.append(j); nom.append(self.shared_norm_nominal[g]); sig.append(self.shared_norm_sigma[g])
            elif key.kind == "weight":
                m = self.mix.member(key.level)
                if m.weight_sigma > 0:
                    cols.append(j); nom.append(m.weight); sig.append(m.weight_sigma)
        base = Penalties(np.array(cols, int), np.array(nom, float), np.array(sig, float))
        sum_cols = [c for c in self.weight_col if c is not None]
        sum_const = sum(m.weight for m in self.members if not m.weight_free)
        self.penalties = MixturePenalties(base, sum_cols if (self.mix.sum_to_one and sum_cols) else None,
                                          sum_const=sum_const, sum_sigma=self.mix.sum_sigma)
        self.nfree = nfree

    # -- member vectors ---------------------------------------------------------

    def member_x(self, x, k: int, shared_norms: bool = False) -> np.ndarray:
        """Member ``k``'s free vector at the global ``x``.  The member's own
        normalization of a shared segment is pinned to 1 (the shared
        normalization is applied in Python); with ``shared_norms`` it takes the
        fitted shared value instead, which is what a standalone file needs."""
        x = np.asarray(x, float)
        mp, pin = self.maps[k], self.pinned[k]
        out = np.where(mp >= 0, x[np.clip(mp, 0, None)], pin)
        if shared_norms:
            info = self.infos[k]
            for j in range(self.nR[k], len(info["keys"])):
                key = info["keys"][j]
                if key.kind == "norm" and int(key.seg) in self.shared_of_member[k]:
                    out[j] = self.shared_norm(x, self.shared_of_member[k][int(key.seg)])
        return out

    def shared_norm(self, x, g: int) -> float:
        c = self.shared_norm_col[g]
        return float(np.asarray(x, float)[c]) if c is not None else self.shared_norm_nominal[g]

    def weight(self, x, k: int) -> float:
        c = self.weight_col[k]
        return float(np.asarray(x, float)[c]) if c is not None else float(self.members[k].weight)

    def weights(self, x) -> Dict[str, float]:
        return {m.name: self.weight(x, k) for k, m in enumerate(self.members)}

    # -- evaluation -------------------------------------------------------------

    def _broadcast(self, cmd, payload_for):
        for k, w in enumerate(self.workers):
            w.send(cmd, payload_for(k))
        out = []
        for k, w in enumerate(self.workers):
            try:
                out.append(w.recv())
            except EvalError as err:
                w.restart()
                try:
                    out.append(w.call(cmd, payload_for(k)))
                except EvalError as err2:
                    raise EvalError(f"{err}; retry after restart also failed: {err2}")
        return out

    def residuals(self, x, jac=False):
        x = np.asarray(x, float)
        t0 = time.time()
        parts = self._broadcast("resid", lambda k: (self.member_x(x, k), bool(jac)))
        nfree = x.size
        r = np.zeros(self.npts)
        J = np.zeros((self.npts, nfree)) if jac else None
        for k, (rk, Jk) in enumerate(parts):
            rk = np.asarray(rk, float)
            for (i, a, b, ga) in self.member_private[k]:
                n = b - a
                r[ga:ga + n] = rk[a:b]
                if jac:
                    mp = self.maps[k]
                    good = mp >= 0
                    J[ga:ga + n][:, mp[good]] = np.asarray(Jk, float)[a:b][:, good]
        for g, mems in enumerate(self.shared_members):
            y, e = self.shared_y[g], self.shared_e[g]
            n = self.shared_norm(x, g)
            ga, gb = self.shared_rows[g]
            F = np.zeros(y.size)
            fks = []
            for (k, a, b) in mems:
                rk = np.asarray(parts[k][0], float)[a:b]
                fk = rk * e + y                       # the member's norm is pinned to 1
                fks.append((k, fk))
                F += self.weight(x, k) * fk
            r[ga:gb] = (F - y * n) / (e * n)
            if jac:
                for (k, a, b), (_, fk) in zip(mems, fks):
                    Gk = np.asarray(parts[k][1], float)[a:b, :self.nR[k]] * e[:, None]
                    cols = self.maps[k][:self.nR[k]]
                    J[ga:gb][:, cols] = self.weight(x, k) * Gk / (e * n)[:, None]
                    wc = self.weight_col[k]
                    if wc is not None:
                        J[ga:gb, wc] = fk / (e * n)
                nc = self.shared_norm_col[g]
                if nc is not None:
                    J[ga:gb, nc] = -F / (e * n * n)
        pen = self.penalties.residuals(x)
        r = np.concatenate([r, pen])
        if jac:
            J = np.vstack([J, self.penalties.jacobian(nfree)])
        if not np.all(np.isfinite(r)):
            raise EvalError("non-finite residuals")
        self.timings.append(time.time() - t0)
        return self._weight_rows(r, J)

    def objective(self, x):
        r, _ = self.residuals(x, jac=False)
        return self.rows_objective(r)

    def transform(self, x):
        x = np.asarray(x, float)
        parts = self._broadcast("transform", lambda k: self.member_x(x, k))
        return np.concatenate([np.asarray(p, float)[:self.nR[k]] for k, p in enumerate(parts)])

    def wigner_widths(self, x):
        x = np.asarray(x, float)
        parts = self._broadcast("wigner_widths", lambda k: self.member_x(x, k))
        out = {}
        for k, p in enumerate(parts):
            for fi, v in p.items():
                if int(fi) < self.nR[k]:
                    out[int(self.maps[k][int(fi)])] = float(v)
        return out

    def data(self):
        out = dict(energies=[], excitation=[], cross=[], cross_err=[], angles=[])
        for k in range(len(self.members)):
            d = self.member_data[k]
            for (i, a, b, ga) in self.member_private[k]:
                for key in out:
                    arr = d.get(key, [])
                    out[key].append(np.asarray(arr[i], float) if i < len(arr) else np.zeros(0))
        for g, mems in enumerate(self.shared_members):
            k, a, b = mems[0]
            i = list(map(int, self.infos[k]["seg_keys"])).index(int(self.mix.shared[g].segments[self.members[k].name]))
            d = self.member_data[k]
            for key in out:
                arr = d.get(key, [])
                out[key].append(np.asarray(arr[i], float) if i < len(arr) else np.zeros(0))
        return out

    def cross(self, x):
        x = np.asarray(x, float)
        parts = self._broadcast("cross", lambda k: self.member_x(x, k))
        return {m.name: p for m, p in zip(self.members, parts)}

    def shared_model(self, x) -> List[dict]:
        """Per shared data set: the energies, data, the total F and each member's w_k f_k at ``x``."""
        x = np.asarray(x, float)
        parts = self._broadcast("resid", lambda k: (self.member_x(x, k), False))
        out = []
        for g, mems in enumerate(self.shared_members):
            y, e = self.shared_y[g], self.shared_e[g]
            k0, a0, b0 = mems[0]
            i = list(map(int, self.infos[k0]["seg_keys"])).index(int(self.mix.shared[g].segments[self.members[k0].name]))
            en = np.asarray(self.member_data[k0]["energies"][i], float)
            terms = {}
            F = np.zeros(y.size)
            for (k, a, b) in mems:
                fk = np.asarray(parts[k][0], float)[a:b] * e + y
                terms[self.members[k].name] = self.weight(x, k) * fk
                F += terms[self.members[k].name]
            out.append(dict(name=self.mix.shared[g].name, energies=en, data=y, error=e, total=F, terms=terms,
                            norm=self.shared_norm(x, g)))
        return out

    def write_output_files(self, x):
        x = np.asarray(x, float)
        self._broadcast("write_output", lambda k: self.member_x(x, k, shared_norms=True))

    def close(self):
        for w in self.workers:
            w.stop()

    # -- checks -----------------------------------------------------------------

    def _self_check(self):
        """The shared-segment model reconstructed from the members' residuals
        equals their direct calculation, and the members agree on the data."""
        x = self.x0
        parts = self._broadcast("resid", lambda k: (self.member_x(x, k), False))
        crosses = self._broadcast("cross", lambda k: self.member_x(x, k))
        for g, mems in enumerate(self.shared_members):
            y, e = self.shared_y[g], self.shared_e[g]
            for (k, a, b) in mems:
                fk = np.asarray(parts[k][0], float)[a:b] * e + y
                i = list(map(int, self.infos[k]["seg_keys"])).index(
                    int(self.mix.shared[g].segments[self.members[k].name]))
                direct = np.asarray(crosses[k][i], float)
                if direct.shape != fk.shape or not np.allclose(direct, fk, rtol=1e-8, atol=1e-12 * max(1.0, np.max(np.abs(direct)))):
                    raise RuntimeError(f"mixture self-check: member {self.members[k].name}'s model on shared "
                                       f"{self.mix.shared[g].name!r} does not reconstruct from its residuals "
                                       f"(max rel diff {np.max(np.abs(direct - fk) / np.maximum(np.abs(direct), 1e-30)):.2e}); "
                                       "is its normalization of that segment free and nominal 1?")
        self.log(f"mixture: {len(self.members)} members, {len(self.mix.shared)} shared data set(s), "
                 f"{self.npts} data rows, {self.nfree} free parameters ({self.n_rmatrix} R-matrix)")

    def describe(self) -> str:
        lines = []
        for k, m in enumerate(self.members):
            nprivate = sum(b - a for (_, a, b, _) in self.member_private[k])
            lines.append(f"  {m.name}: {self.nR[k]} R-matrix parameters, {nprivate} private points, weight "
                         f"{m.weight:g}{' (free)' if m.weight_free else ''}")
        for g, s in enumerate(self.mix.shared):
            a, b = self.shared_rows[g]
            lines.append(f"  shared {s.name}: {b - a} points, members {', '.join(s.segments)}, "
                         f"norm {'free' if self.shared_norm_col[g] is not None else 'fixed'} "
                         f"(nominal {self.shared_norm_nominal[g]:g})")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# freeze and exchange: every member as a standalone stage-1 file

def _run_pyazr(code: str, cwd: str, timeout=3600) -> dict:
    """Run ``code`` in a fresh python process (one engine session per process)
    with ``result`` (a dict) returned as JSON."""
    prog = ("import os, sys, json\n"
            f"sys.path.insert(0, {AZURE2_ROOT!r})\n"
            f"sys.path.insert(0, {os.path.dirname(os.path.dirname(os.path.abspath(__file__)))!r})\n"
            "import numpy as np\n"
            "from pyazr import azure2, AzrModel\n"
            "result = {}\n" + code + "\nprint('@@RESULT@@' + json.dumps(result))\n")
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "1")
    p = subprocess.run([sys.executable, "-c", prog], cwd=cwd, capture_output=True, text=True,
                       timeout=timeout, env=env)
    tag = "@@RESULT@@"
    for line in p.stdout.splitlines():
        if line.startswith(tag):
            return json.loads(line[len(tag):])
    raise RuntimeError(f"pyazr subprocess failed in {cwd}:\n{p.stdout[-2000:]}\n{p.stderr[-4000:]}")


def table_grid(energies: np.ndarray, step: Optional[float] = None, margin: float = 0.02,
               max_points: int = 20000) -> np.ndarray:
    """Lab energies for a frozen table: the data energies plus a uniform grid
    over their span (default step: span/4000, floored so the table stays
    under ``max_points``) widened by ``margin`` (relative) at both ends, so a
    kernel window at the edge is covered."""
    e = np.unique(np.asarray(energies, float))
    lo, hi = e.min(), e.max()
    span = max(hi - lo, 1e-6)
    lo, hi = lo - margin * span, hi + margin * span
    lo = max(lo, 1e-6)
    if step is None:
        step = span / 4000.0
    step = max(step, (hi - lo) / max_points)
    grid = np.arange(lo, hi + step / 2, step)
    # rounded to 1 eV: two points closer than that would print alike and the
    # engine refuses a table whose energies do not increase strictly
    return np.unique(np.round(np.concatenate([grid, e]), 9))


def freeze_and_exchange(ev: Mixture, x, outdir: str, step: Optional[float] = None,
                        margin: float = 0.02, verify: bool = True, log=print) -> dict:
    """Write every member as a standalone stage-1 file reproducing the joint fit.

    For each member: its fitted ``.azr`` (``save_fit``, every level parameter
    at its fitted value, the shared normalization in the ``.sav`` and on the
    line), then, per shared data set, a ``mixture`` block on that line with
    the member's own weight as ``w_self`` and the other members' fitted
    contributions as tables (bare cross section on :func:`table_grid`, in lab
    energy).  With ``verify`` each standalone file is re-evaluated and its
    shared-segment chi-squared compared with the joint fit's.
    """
    x = np.asarray(x, float)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    report = dict(members={}, tables={}, verify={})
    joint = ev.score(x)
    members = ev.members
    # 1. standalone fitted files
    fitted = {}
    for k, m in enumerate(members):
        mdir = os.path.join(outdir, m.name)
        os.makedirs(mdir, exist_ok=True)
        from .model import link_data_dirs
        link_data_dirs(m.azr, mdir)
        os.makedirs(os.path.join(mdir, "output"), exist_ok=True)
        xk = ev.member_x(x, k, shared_norms=True)
        out_azr = os.path.join(mdir, f"{m.name}.azr")
        res = _run_pyazr(f"""
xk = np.array({xk.tolist()!r}, float)
with azure2({ev.member_azr[m.name]!r}, cwd={os.path.dirname(ev.member_azr[m.name])!r}, **{m.flags!r}) as s:
    azr, sav = s.save_fit({out_azr!r}, x=xk, param_sav=True, verify=True)
    result['azr'] = azr; result['sav'] = sav
    result['chi2'] = [float(v) for v in s.segment_chi2(xk)]
    result['seg_keys'] = [int(d.key) for d in s.active_datasets]
""", cwd=os.path.dirname(ev.member_azr[m.name]))
        text = relocate_data_paths(open(res["azr"]).read(), os.path.dirname(ev.member_azr[m.name]), mdir)
        with open(res["azr"], "w") as fh:
            fh.write(text)
        fitted[m.name] = res
        report["members"][m.name] = dict(azr=res["azr"], sav=res["sav"])
        log(f"freeze: {m.name}: {os.path.relpath(res['azr'], outdir)} written and verified")
    # 2. tables of every member's contribution to every shared data set
    tdir = os.path.join(outdir, "tables")
    os.makedirs(tdir, exist_ok=True)
    tables: Dict[Tuple[int, str], str] = {}
    for g, s in enumerate(ev.mix.shared):
        k0, a0, b0 = ev.shared_members[g][0]
        i0 = list(map(int, ev.infos[k0]["seg_keys"])).index(int(s.segments[members[k0].name]))
        en = np.asarray(ev.member_data[k0]["energies"][i0], float)
        d0 = ev.member_segments[k0].by_key(s.segments[members[k0].name])
        # the engine's energies are c.m.; the data file and the table are lab
        cm_per_lab = _cm_per_lab(ev.member_azr[members[k0].name], d0.entrance_key)
        grid = table_grid(en / cm_per_lab, step=step, margin=margin)
        angle = float(d0.angle_min)
        gridfile = os.path.join(tdir, f"grid_{s.name}.dat")
        with open(gridfile, "w") as f:
            for e in grid:
                f.write(f"{e:.9f} {angle:g} 1.0 1.0\n")
        for (k, a, b) in ev.shared_members[g]:
            m = members[k]
            seg_key = int(s.segments[m.name])
            fit_azr = fitted[m.name]["azr"]
            text = open(fit_azr).read()
            lines = block_lines(text, "segmentsData")
            keep = None
            for ln in lines:
                t = ln.split()
                if not t:
                    continue
                # the key is the position among non-empty lines
                if lines_key(lines, ln) == seg_key:
                    keep = ln
                    break
            if keep is None:
                raise RuntimeError(f"freeze: member {m.name}: segment {seg_key} not found in {fit_azr}")
            keep = _replace_data_file(keep, gridfile)
            keep = _force_plain_norm(keep, grid[0] * 0.999, grid[-1] * 1.001)
            var = replace_block(text, "segmentsData", [keep])
            var = replace_block(var, "segmentsTest", [])
            var = replace_block(var, "targetInt", [])
            vdir = os.path.join(tdir, f"var_{s.name}_{m.name}")
            os.makedirs(os.path.join(vdir, "output"), exist_ok=True)
            vazr = os.path.join(vdir, f"{m.name}.azr")
            var = set_output_dir_text(var, os.path.join(vdir, "output"))
            with open(vazr, "w") as f:
                f.write(var)
            res = _run_pyazr(f"""
with azure2({vazr!r}, cwd={vdir!r}, **{m.flags!r}) as s:
    vals = np.asarray(s.calculate_rwa(s.params_rwa)[0], float)
    result['values'] = vals.tolist()
""", cwd=vdir)
            vals = np.asarray(res["values"], float)
            if vals.size != grid.size:
                raise RuntimeError(f"freeze: table of {m.name} on {s.name}: {vals.size} values for {grid.size} grid points")
            tpath = os.path.join(tdir, f"{s.name}_{m.name}.dat")
            with open(tpath, "w") as f:
                f.write(f"# {m.name} contribution to {s.name}: bare model of {os.path.relpath(fit_azr, outdir)} "
                        f"at the joint fit, lab energy [MeV], cross section as the segment measures it\n")
                for e, v in zip(grid, vals):
                    f.write(f"{e:.9f} {v:.8e}\n")
            tables[(g, m.name)] = tpath
            report["tables"][f"{s.name}/{m.name}"] = tpath
            log(f"freeze: table {os.path.relpath(tpath, outdir)} ({grid.size} points)")
    # 3. the mixture block on every member's shared line
    for g, s in enumerate(ev.mix.shared):
        for (k, a, b) in ev.shared_members[g]:
            m = members[k]
            fit_azr = fitted[m.name]["azr"]
            others = [(os.path.relpath(tables[(g, o.name)], os.path.dirname(fit_azr)), ev.weight(x, ko))
                      for ko, o in enumerate(members) if o.name != m.name and (g, o.name) in tables]
            mdl = AzrModel_from(fit_azr)
            mdl.set_mixture(int(s.segments[m.name]), ev.weight(x, k), others)
            mdl.write(fit_azr)
    # 4. verify: each standalone reproduces the joint fit's shared chi-squared
    if verify:
        for g, s in enumerate(ev.mix.shared):
            gid = shared_gid(g)
            joint_chi2 = joint.segment_chi2.get(gid, float("nan"))
            for (k, a, b) in ev.shared_members[g]:
                m = members[k]
                fit_azr = fitted[m.name]["azr"]
                sav = fitted[m.name]["sav"]
                res = _run_pyazr(f"""
with azure2({fit_azr!r}, cwd={os.path.dirname(fit_azr)!r}, **{m.flags!r}) as s:
    xk = s.load_sav({sav!r}) if {sav!r} else np.asarray(s.params_rwa, float)
    xk = np.asarray(xk, float)
    result['chi2'] = [float(v) for v in s.segment_chi2(xk)]
    result['seg_keys'] = [int(d.key) for d in s.active_datasets]
""", cwd=os.path.dirname(fit_azr))
                i = res["seg_keys"].index(int(s.segments[m.name]))
                own = res["chi2"][i]
                ok = abs(own - joint_chi2) <= 1e-3 * max(1.0, abs(joint_chi2))
                report["verify"][f"{s.name}/{m.name}"] = dict(standalone=own, joint=joint_chi2, ok=ok)
                log(f"freeze: verify {s.name}/{m.name}: standalone chi2 {own:.4f} vs joint {joint_chi2:.4f} "
                    f"{'ok' if ok else 'MISMATCH'}")
    with open(os.path.join(outdir, "freeze_report.json"), "w") as f:
        json.dump(report, f, indent=1)
    return report


def AzrModel_from(path):
    from pyazr.azrfile import AzrModel
    return AzrModel.from_file(path)


def lines_key(lines: Sequence[str], line: str) -> int:
    """1-based key of a <segmentsData> line = its position among the non-empty lines."""
    k = 0
    for ln in lines:
        if ln.strip():
            k += 1
            if ln is line:
                return k
    return -1


def relocate_data_paths(text: str, old_dir: str, new_dir: str) -> str:
    """Rewrite every relative data-file path of the <segmentsData> block so a
    copy of the file in ``new_dir`` finds the files the original in ``old_dir``
    did.  The engine drops a segment whose file is missing without a word, and
    a project that keeps a shared file outside its own data/ (as the natural
    neon members do) would lose it on every copy."""
    lines = block_lines(text, "segmentsData")
    out = []
    for ln in lines:
        t = ln.split()
        if not t:
            out.append(ln)
            continue
        for i, tok in enumerate(t):
            try:
                float(tok)
            except ValueError:
                if not os.path.isabs(tok):
                    t[i] = os.path.relpath(os.path.normpath(os.path.join(old_dir, tok)), new_dir)
                break
        out.append(" ".join(t))
    return replace_block(text, "segmentsData", out)


def _replace_data_file(line: str, new_path: str) -> str:
    t = line.split()
    for i, tok in enumerate(t):
        try:
            float(tok)
        except ValueError:
            t[i] = new_path
            return " ".join(t)
    raise ValueError(f"no data file token in segment line: {line}")


def _force_plain_norm(line: str, e_lo: Optional[float] = None, e_hi: Optional[float] = None) -> str:
    """Normalization 1, not varied (the table is the bare model), and an energy
    window that covers the whole table grid (the data's window would drop the
    grid's margin points)."""
    t = line.split()
    if e_lo is not None:
        t[3] = f"{e_lo:.7g}"
    if e_hi is not None:
        t[4] = f"{e_hi:.7g}"
    is_diff = int(t[7])
    j = 8 + (2 if is_diff == 2 else 0)
    t[j], t[j + 1], t[j + 2] = "1", "0", "0"
    return " ".join(t)


def _cm_per_lab(azr: str, entrance_key: int) -> float:
    from .model import cm_conversion
    return cm_conversion(open(azr).read(), entrance_key)


# ---------------------------------------------------------------------------
# driver

def _state_path(d):
    return os.path.join(d, "state.json")


def _load_state(d):
    p = _state_path(d)
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return None


def _save_state(d, ev: Mixture, x, score, note=""):
    st = _load_state(d) or dict(history=[])
    st.update(x=[float(v) for v in x], keys=[str(k) for k in ev.keys],
              objective=score.objective, chi2=score.chi2, penalty=score.penalty, npts=score.npts,
              weights=ev.weights(x), norms={ev.mix.shared[g].name: ev.shared_norm(x, g) for g in range(len(ev.mix.shared))},
              segment_chi2={str(k): v for k, v in score.segment_chi2.items()},
              updated=time.strftime("%Y-%m-%d %H:%M"))
    st.setdefault("history", []).append(dict(time=st["updated"], objective=score.objective, chi2=score.chi2, note=note))
    with open(_state_path(d), "w") as f:
        json.dump(st, f, indent=1)


def main(argv=None):
    import argparse
    from .local import polish, default_bounds, BoundsPolicy
    ap = argparse.ArgumentParser(prog="python3 -m rmfit.mixture")
    ap.add_argument("directory")
    ap.add_argument("command", choices=["init", "polish", "freeze", "status", "score"])
    ap.add_argument("--max-nfev", type=int, default=100)
    ap.add_argument("--threads", type=int, default=1, help="OpenMP threads per member worker")
    ap.add_argument("--step", type=float, default=None, help="table grid step (MeV, lab) for freeze")
    ap.add_argument("--margin", type=float, default=0.02,
                    help="table grid margin beyond the data span, relative (cover the resolution kernel)")
    ap.add_argument("--out", default=None, help="freeze output directory (default <dir>/frozen)")
    ap.add_argument("--no-verify", action="store_true")
    args = ap.parse_args(argv)
    d = os.path.abspath(args.directory)
    spec = MixtureSpec.load(d)
    if args.command == "status":
        st = _load_state(d)
        print(f"mixture in {d}: {len(spec.members)} members, {len(spec.shared)} shared data sets")
        if st:
            print(f"  objective {st['objective']:.4f} (chi2 {st['chi2']:.4f}, penalty {st['penalty']:.4f}), "
                  f"weights {st['weights']}, norms {st['norms']}, updated {st['updated']}")
            for h in st.get("history", [])[-5:]:
                print(f"    {h['time']}  {h['objective']:.4f}  {h.get('note', '')}")
        return 0
    with Mixture(spec, omp_threads=args.threads, log=print) as ev:
        print(ev.describe())
        st = _load_state(d)
        x = np.asarray(st["x"], float) if (st and args.command != "init" and len(st["x"]) == ev.nfree) else ev.x0.copy()
        if args.command == "init":
            sc = ev.score(x)
            print(f"init: objective {sc.objective:.4f}, chi2 {sc.chi2:.4f} on {sc.npts} points, penalty {sc.penalty:.4f}")
            for g, s in enumerate(spec.shared):
                print(f"  shared {s.name}: chi2 {sc.segment_chi2.get(shared_gid(g), float('nan')):.4f}")
            _save_state(d, ev, x, sc, note="init")
        elif args.command == "score":
            sc = ev.score(x)
            print(f"score: objective {sc.objective:.4f}, chi2 {sc.chi2:.4f}, penalty {sc.penalty:.4f}")
        elif args.command == "polish":
            policy = BoundsPolicy()
            bounds = default_bounds(ev, x, policy)
            res = polish(ev, x, bounds=bounds, policy=policy, max_nfev=args.max_nfev, verbose=1)
            sc = ev.score(res.x)
            print(f"polish: {res.start_objective:.4f} -> {sc.objective:.4f} ({res.status}); weights {ev.weights(res.x)}")
            _save_state(d, ev, res.x, sc, note=f"polish {res.status}")
        elif args.command == "freeze":
            out = args.out or os.path.join(d, "frozen")
            freeze_and_exchange(ev, x, out, step=args.step, margin=args.margin, verify=not args.no_verify, log=print)
    return 0


if __name__ == "__main__":
    sys.exit(main())
