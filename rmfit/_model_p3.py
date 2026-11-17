"""rmfit.model on pyazr 3: project format 3 through :mod:`rmfit.azr3` (data sets as records, shards and variants
as ``Project`` edits; the engine folds, so there are no convolution pieces).  Imported by :mod:`rmfit.model`;
import from there."""
from __future__ import annotations

import os
import re
import shutil
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple
import numpy as np
from . import AZURE2_ROOT  # noqa: F401  (puts pyazr on sys.path)
from .azr3 import SegmentTable as SegmentSet, Segment, data_paths, read_points, write_variant

from .model import (ShardSpec, free_shift, segment_class, pack, pack_units, shift_columns, shard_shift_count, pack_jac, _build_shards_balanced, calibrate_weights)  # noqa: F401  (shared part of rmfit.model)

__all__ = ['SegmentSet', 'Segment', 'rows_in_segment', 'points_in_segment', 'Piece', 'DEFAULT_WEIGHTS', 'plan_pieces', 'data_dirs', 'link_data_dirs', 'write_shard', 'build_shards', 'plan_balanced', 'single_variant', 'base_segments']


def rows_in_segment(seg: Segment, base_dir: str) -> Optional[np.ndarray]:
    """The data rows AZURE2 loads for this data set, sorted by energy (columns: lab energy in MeV, angle), or
    None when rmfit cannot place the rows of its file (:func:`rmfit.azr3.read_points`)."""
    pts = read_points(seg, base_dir)
    if pts is None:
        return None
    if pts.size == 0:
        return pts.reshape(0, 2)
    e, a = pts[:, 0], pts[:, 1]
    sel = (e >= seg.energy_min) & (e <= seg.energy_max)
    if seg.differential:
        tol = 1e-6
        sel &= (a >= seg.angle_min - tol) & (a <= seg.angle_max + tol)
    rows = pts[sel]
    return rows[np.argsort(rows[:, 0], kind="stable")]


def points_in_segment(seg: Segment, base_dir: str) -> Optional[np.ndarray]:
    """Lab energies (MeV) of the data rows AZURE2 loads for this data set (sorted), or None
    (:func:`rows_in_segment`)."""
    rows = rows_in_segment(seg, base_dir)
    return None if rows is None else rows[:, 0]


# -- folded data sets ---------------------------------------------------------------
#
# pyazr 2: a data set whose <targetInt> entry was a pure Gaussian convolution became a "conv" piece: the shard
# evaluated the bare cross section on a grid and rmfit folded it.  On pyazr 3.0 the engine folds every data set
# with its own kernels on the mesh of its session, so there is no conv piece: a folded data set is a plain piece,
# never split.


@dataclass
class Piece:
    base_key: int
    record: str                 # the name of the base model's dataset record
    npts: int                   # rows this piece contributes to the assembled residual
    cost: float
    cls: str
    e_lo: Optional[float] = None
    e_hi: Optional[float] = None
    kind: str = "data"
    split: bool = False         # one of several energy pieces of its data set
    name: str = ""              # the record's name in the shard file (write_shard)
    norm0: Optional[float] = None   # a normalization-group follower's fixed norm in the shard (write_shard)


DEFAULT_WEIGHTS = {"2->2:dsigma_domega": 8.5, "2->1:dsigma_domega": 2.3}


def plan_pieces(base_azr: str, weights: Optional[Dict[str, float]] = None,
                target_cost: Optional[float] = None,
                no_split: Sequence[int] = (), conv_grid: bool = True,
                split_composite: Optional[bool] = None,
                seg_factors: Optional[Dict[int, float]] = None, default_factor: float = 1.0) -> List[Piece]:
    """One :class:`Piece` per active data set, big ones cut into energy pieces.

    A data set is never split if it carries a free energy shift, is folded by the engine (its record names a
    target, a beam or a detector), is a ratio, has a file whose rows rmfit cannot place
    (:func:`rmfit.azr3.read_points`) -- or if its key is in ``no_split``.  A sum of reactions is split like a
    plain one (``split_composite``, default on; RMFIT_SPLIT_SUM=0 turns it off) and costs one calculation per
    term.  ``conv_grid`` is accepted and ignored (pyazr 2's conv pieces; see above).
    """
    base_dir = os.path.dirname(os.path.abspath(base_azr))
    weights = dict(DEFAULT_WEIGHTS, **(weights or {}))
    segs = SegmentSet.from_file(base_azr)
    if split_composite is None:
        split_composite = os.environ.get("RMFIT_SPLIT_SUM", "1") != "0"
    pieces: List[Piece] = []
    for seg in segs:
        if not seg.active:
            continue
        energies = points_in_segment(seg, base_dir)
        n = int(seg.rows_kept if energies is None else energies.size)
        if energies is not None and n != seg.rows_kept:
            # rmfit's reading of the file disagrees with the engine's: trust the engine and do not split
            energies, n = None, int(seg.rows_kept)
        cls = segment_class(seg)
        w = weights.get(cls, 1.0)
        if seg_factors is not None:            # measured cost correction (Jacobian balancing): seconds per cost unit
            w *= float(seg_factors.get(int(seg.key), default_factor))
        if split_composite and seg.operation == "sum":
            w *= max(1, len(seg.components))       # one calculation per term
        cost = n * w
        splittable = (energies is not None and seg.key not in no_split and not free_shift(seg)
                      and not seg.folded and not seg.shift_group
                      and (seg.operation is None or (split_composite and seg.operation == "sum"))
                      and n >= 40)
        if target_cost is None or cost <= target_cost or not splittable:
            pieces.append(Piece(seg.key, seg.name, n, cost, cls, seg.energy_min, seg.energy_max))
            continue
        k = int(np.ceil(cost / target_cost))
        # cut between distinct energies so no point straddles two pieces
        uniq = np.unique(energies)
        k = min(k, uniq.size)
        cuts = [uniq[int(round(j * uniq.size / k))] for j in range(1, k)]
        edges = [seg.energy_min]
        for c in cuts:
            prev = uniq[uniq < c]
            if prev.size == 0:
                continue
            edges.append(0.5 * (prev[-1] + c))
        edges.append(seg.energy_max)
        for lo, hi in zip(edges[:-1], edges[1:]):
            if hi <= lo:
                continue
            m = int(np.sum((energies >= lo) & (energies <= hi)))
            if m == 0:
                continue
            pieces.append(Piece(seg.key, seg.name, m, m * w, cls, float(lo), float(hi), split=True))
    return pieces


def data_dirs(azr_path: str) -> List[str]:
    """Directories, relative to the model, that hold the files its records read (data files, covariance files,
    grid files, the prior covariance): the first path component of each such name that is an existing directory
    beside the model, with ``data`` always first."""
    base_dir = os.path.dirname(os.path.abspath(azr_path))
    out = ["data"]
    for name in data_paths(azr_path):
        first = os.path.normpath(name).split(os.sep)[0]
        if first in ("", ".", "..") or first in out:
            continue
        if os.path.isdir(os.path.join(base_dir, first)):
            out.append(first)
    return out


def link_data_dirs(azr_path: str, workdir: str, src_dir: Optional[str] = None) -> List[str]:
    """Link every directory of :func:`data_dirs` from ``src_dir`` (default: the model's own directory) into
    ``workdir``, and every file the model names beside itself (a covariance file, a prior covariance), so a copy
    of the model written there finds all its files.  Existing entries are left alone.  Returns the names linked
    or already present."""
    src_dir = os.path.abspath(src_dir or os.path.dirname(os.path.abspath(azr_path)))
    os.makedirs(workdir, exist_ok=True)
    done = []
    for sub in data_dirs(azr_path):
        src, dst = os.path.join(src_dir, sub), os.path.join(workdir, sub)
        if os.path.isdir(src) and not os.path.lexists(dst):
            os.symlink(src, dst)
        if os.path.lexists(dst):
            done.append(sub)
    for name in data_paths(azr_path):
        if os.sep in os.path.normpath(name):
            continue                              # inside a linked directory (or above the model: left alone)
        src, dst = os.path.join(src_dir, name), os.path.join(workdir, name)
        if os.path.isfile(src) and not os.path.lexists(dst):
            os.symlink(src, dst)
    return done


def write_shard(base_azr: str, name: str, pieces: List[Piece], workdir: str) -> ShardSpec:
    """Write ``<workdir>/<name>.azr`` holding only ``pieces`` as its data sets (record edits of the base project,
    :func:`rmfit.azr3.write_variant`).

    Every piece is its base record with the piece's energy window; no ``covariance`` (the evaluator whitens the
    assembled rows).  A normalization shared by several data sets (a ``normalization`` record) is written out per
    data set, because a group's members may sit in different shards: the group's first member carries the
    record's nuisance as its own ``norm``; a follower gets the fixed norm 1 (``Piece.norm0``) and the evaluator
    rescales its rows to the leader's (rmfit.covariance).  A shared energy scale stays as it is: its members share
    a shard.  Extrapolations, model data and a prior covariance belong to the whole model and are left out."""
    from .covariance import NormGroups
    import pyazr
    segs = SegmentSet.from_file(base_azr)
    groups = NormGroups.from_azr(base_azr, segments=segs)
    project = pyazr.Project.open(os.path.abspath(base_azr))
    os.makedirs(workdir, exist_ok=True)
    # data files are referenced relative to the project: link every directory they live in (data/ and any other)
    link_data_dirs(base_azr, workdir)
    seen: Dict[int, int] = {}
    datasets = []
    for p in pieces:
        seg = segs.by_key(p.base_key)
        seen[p.base_key] = seen.get(p.base_key, 0) + 1
        p.name = seg.name if seen[p.base_key] == 1 else f"{seg.name}__p{seen[p.base_key]}"
        fields: Dict[str, object] = {"covariance": None}
        if p.split:
            fields["energy_window"] = (repr(float(p.e_lo)), repr(float(p.e_hi)), "MeV")
        if seg.norm_group:
            if p.base_key in groups.leader:
                p.norm0 = 1.0
                fields["norm"] = pyazr.Nuisance(1.0, "fixed")
            else:
                fields["norm"] = project[seg.norm_group].get("value")
        datasets.append((seg.name, p.name, fields))
    outdir = os.path.join(workdir, f"output_{name}")
    os.makedirs(outdir, exist_ok=True)
    for f in ("intEC.dat", "intEC.extrap"):
        fp = os.path.join(outdir, f)
        if os.path.exists(fp):
            os.remove(fp)
    path, keys = write_variant(base_azr, os.path.join(workdir, f"{name}.azr"), outdir, datasets=datasets)
    ordered = sorted(pieces, key=lambda p: keys[p.name])
    return ShardSpec(name=name, azr=path, output_dir=outdir, pieces=ordered)


def build_shards(base_azr: str, nshards: int, workdir: str,
                 weights: Optional[Dict[str, float]] = None,
                 no_split: Sequence[int] = (), conv_grid: bool = True,
                 shift_solo: Optional[bool] = None, balance: Optional[dict] = None) -> List[ShardSpec]:
    """Plan, pack and write ``nshards`` shard files of ``base_azr``.

    ``balance`` (2026-10-26): a measured Jacobian cost model ``{"alpha", "beta", "factors": {segment key: s per
    cost unit}, "default_factor"}`` (Sharded measures it, ``engine.Sharded._rebalance``); the pieces are then
    costed in seconds and packed by :func:`pack_jac` (shift_solo is not used).

    ``shift_solo`` (default on; RMFIT_SHIFT_SOLO=0 turns it off): every
    segment with a free energy shift gets a shard of its own, the rest is
    packed into the others.  A shift is finite-differenced in the Jacobian
    and each step re-evaluates its whole shard (2026-10-10, a large merged
    model, 16 shards: shards holding the 6 shifted segments took 28-92 s of
    a 90 s Jacobian, the others 10-17 s; solo, the Jacobian is ~30 s).
    Used only while the shifted segments are fewer than half the shards."""
    if balance is not None:
        return _build_shards_balanced(base_azr, nshards, workdir, weights, no_split, conv_grid, balance)
    if shift_solo is None:
        shift_solo = os.environ.get("RMFIT_SHIFT_SOLO", "1") != "0"
    pieces = plan_pieces(base_azr, weights=weights, conv_grid=conv_grid)
    total = sum(p.cost for p in pieces)
    target = total / nshards
    from .covariance import shift_groups, shift_fn_groups
    sg = shift_groups(base_azr).units()
    # 2026-10-24: a shift-function group's members use the function of its first active member, which declares it;
    # a split member would declare it twice in one shard (an engine error) or not at all, so members are never split
    fg = {}     # format 3: one energy_scale record holds the offset and the shift function (sg has its members)
    both = sorted(set(sg) & set(fg))
    if both:
        raise ValueError(f"segments {both} carry both shiftgroup and shiftfngroup (the engine refuses that)")
    pieces = plan_pieces(base_azr, weights=weights, target_cost=target * 0.5,
                         no_split=tuple(no_split) + tuple(fg), conv_grid=conv_grid)
    unit_of = {k: ("g", l) for k, l in sg.items()}
    unit_of.update({k: ("f", l) for k, l in fg.items()})
    if unit_of:
        # 2026-10-22: a shift group's members follow the constant shift of its first member (the engine, within
        # one file), so all their pieces go to ONE shard, which keeps the keyword: the group is exact there, and a
        # finite-difference step of the group's shift moves every member.  A group with a free shift is one solo unit.
        # 2026-10-24: the same for a shift-function group (shiftfngroup); bins keep file order, so the declaring
        # member stays the first of its group in the shard file.
        units: Dict[object, List[Piece]] = {}
        for p in pieces:
            units.setdefault(unit_of.get(p.base_key, ("p", id(p))), []).append(p)
        units_l = list(units.values())
        solo_u = []
        if shift_solo:
            shift_keys = {s.key for s in SegmentSet.from_file(base_azr) if s.active and free_shift(s)}
            solo_u = [u for u in units_l if any(p.base_key in shift_keys for p in u)]
        if solo_u and len(solo_u) <= nshards // 2:
            rest = [u for u in units_l if all(u is not v for v in solo_u)]
            bins = [sorted(u, key=lambda p: (p.base_key, p.e_lo if p.e_lo is not None else 0.0)) for u in solo_u] + \
                pack_units(rest, nshards - len(solo_u))
        else:
            bins = pack_units(units_l, nshards)
    else:
        solo = []
        if shift_solo:
            shift_keys = {s.key for s in SegmentSet.from_file(base_azr) if s.active and free_shift(s)}
            solo = [p for p in pieces if p.base_key in shift_keys]
        if solo and len(solo) <= nshards // 2:
            rest = [p for p in pieces if p.base_key not in shift_keys]
            bins = [[p] for p in solo] + pack(rest, nshards - len(solo))
        else:
            bins = pack(pieces, nshards)
    shards = []
    for j, b in enumerate(bins):
        if not b:
            continue
        shards.append(write_shard(base_azr, f"shard{j:02d}", b, workdir))
    return shards


def plan_balanced(base_azr: str, nshards: int, weights, no_split, conv_grid, balance: dict):
    """The bins and predicted times of a Jacobian-balanced plan (no files written)."""
    from .covariance import shift_groups, shift_fn_groups
    f = {int(k): float(v) for k, v in (balance.get("factors") or {}).items()}
    dflt = float(balance.get("default_factor", 1.0))
    alpha, beta = float(balance["alpha"]), float(balance["beta"])
    sg = shift_groups(base_azr).units()
    fg = {}     # format 3: one energy_scale record holds the offset and the shift function (sg has its members)
    p0 = plan_pieces(base_azr, weights=weights, conv_grid=conv_grid, seg_factors=f, default_factor=dflt)
    total = sum(p.cost for p in p0)
    frac = float(os.environ.get("RMFIT_BALANCE_PIECE", "0.35"))      # split target as a fraction of a shard's share
    pieces = plan_pieces(base_azr, weights=weights, target_cost=frac * total / nshards,
                         no_split=tuple(no_split) + tuple(fg), conv_grid=conv_grid, seg_factors=f,
                         default_factor=dflt)
    segs = SegmentSet.from_file(base_azr)
    shifted = {s.key for s in segs if s.active and free_shift(s)}
    unit_of = {k: ("g", l) for k, l in sg.items()}
    unit_of.update({k: ("f", l) for k, l in fg.items()})
    units: Dict[object, List[Piece]] = {}
    for p in pieces:
        key = unit_of.get(p.base_key) or (("s", p.base_key) if p.base_key in shifted else ("p", id(p)))
        units.setdefault(key, []).append(p)
    kc = lambda u: shard_shift_count(base_azr, u, segs)
    return pack_jac(list(units.values()), nshards, alpha, beta, kc)


def single_variant(base_azr: str, name: str, workdir: str) -> str:
    """A copy of the base model with its own output directory (for an
    authoritative single session that must not touch the campaign's output/)."""
    os.makedirs(workdir, exist_ok=True)
    link_data_dirs(base_azr, workdir)
    outdir = os.path.join(workdir, f"output_{name}")
    os.makedirs(outdir, exist_ok=True)
    for f in ("intEC.dat", "intEC.extrap"):
        fp = os.path.join(outdir, f)
        if os.path.exists(fp):
            os.remove(fp)
    path, _ = write_variant(base_azr, os.path.join(workdir, f"{name}.azr"), outdir)
    return path


def base_segments(base_azr: str) -> SegmentSet:
    return SegmentSet.from_file(base_azr)
