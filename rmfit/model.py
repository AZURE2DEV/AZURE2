"""Model variants written as new ``.azr`` files.

The engine reads its model from a file, so every variant -- a shard holding a
subset of the data segments, an energy window, a level switched off -- is a
new ``.azr`` written next to the campaign and opened in a fresh session.  The
original is never modified.

Shards are the workhorse: chi-squared is exactly additive over data segments
and the engine's cost is per data point (and does not parallelize with OpenMP),
so splitting the segments across N single-thread sessions divides the cost of
one objective or Jacobian by ~N.  Big segments are cut into energy pieces so
the shards balance; pieces of one segment share its normalization, which the
:mod:`rmfit.engine` assembler sums back into one column.

Two backends (:mod:`rmfit.backend`): the shared planning and packing are here; the file work (reading data rows,
writing shards and variants, the data-set table) is ``_model_p2`` (pyazr 2, ``.azr`` text) or ``_model_p3`` (pyazr 3,
project format 3 through :mod:`rmfit.azr3`), imported at the end of this module.
"""
from __future__ import annotations

import os
import re
import shutil
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple
import numpy as np



@dataclass
class ShardSpec:
    name: str
    azr: str
    output_dir: str
    pieces: List[Piece] = field(default_factory=list)
    predicted: Optional[float] = None   # Jacobian-balanced plans: the predicted resid+Jacobian time (s)

    @property
    def seg_map(self) -> Dict[int, int]:
        """shard segment key (1-based, file order) -> base segment key"""
        return {i + 1: p.base_key for i, p in enumerate(self.pieces)}

    @property
    def npts(self):
        return sum(p.npts for p in self.pieces)

    @property
    def cost(self):
        return sum(p.cost for p in self.pieces)


def free_shift(seg: Segment) -> bool:
    """Any free energy-shift term: the constant, the sqrt(E) coefficient or a shift-function parameter.
    Such a segment is finite-differenced in the Jacobian, so it is never split and gets a shift-solo shard."""
    return bool(seg.vary_shift or getattr(seg, "vary_shift_sqrt", False)
                or any(v for _, _, v in getattr(seg, "shift_fn_params", ())))


def segment_class(seg: Segment) -> str:
    return f"{seg.entrance_key}->{seg.exit_key}:{seg.observable}"


#: Relative cost per data point by segment class, calibrated on a large
#: model (rmfit.model.calibrate_weights): fixed-angle differential segments
#: cost more per point than angle-integrated ones.


def pack(pieces: List[Piece], nshards: int) -> List[List[Piece]]:
    """Longest-processing-time greedy bin packing by cost."""
    bins: List[List[Piece]] = [[] for _ in range(nshards)]
    load = np.zeros(nshards)
    for p in sorted(pieces, key=lambda p: -p.cost):
        j = int(np.argmin(load))
        bins[j].append(p)
        load[j] += p.cost
    # keep file order inside a shard (stable row order for bookkeeping)
    return [sorted(b, key=lambda p: (p.base_key, p.e_lo if p.e_lo is not None else 0.0))
            for b in bins]


def pack_units(units: List[List[Piece]], nshards: int) -> List[List[Piece]]:
    """:func:`pack` for units that must share a shard (a shift group's pieces): LPT by unit cost."""
    bins: List[List[Piece]] = [[] for _ in range(nshards)]
    load = np.zeros(nshards)
    for u in sorted(units, key=lambda u: -sum(p.cost for p in u)):
        j = int(np.argmin(load))
        bins[j].extend(u)
        load[j] += sum(p.cost for p in u)
    return [sorted(b, key=lambda p: (p.base_key, p.e_lo if p.e_lo is not None else 0.0))
            for b in bins]



# -- data directories ----------------------------------------------------------


def shift_columns(seg: Segment) -> int:
    """Number of finite-differenced (shift-type) Jacobian columns a segment carries: the constant shift, the
    sqrt(E) coefficient and each free shift-function parameter."""
    return (int(bool(seg.vary_shift)) + int(bool(getattr(seg, "vary_shift_sqrt", False)))
            + sum(1 for _, _, v in getattr(seg, "shift_fn_params", ()) if v))


def shard_shift_count(base_azr: str, pieces: Sequence[Piece], segs=None) -> int:
    """Shift-type columns a shard holding ``pieces`` finite-differences (each segment counted once)."""
    segs = segs if segs is not None else SegmentSet.from_file(base_azr)
    return sum(shift_columns(segs.by_key(k)) for k in sorted({int(p.base_key) for p in pieces}))


def pack_jac(units: List[List[Piece]], nshards: int, alpha: float, beta: float,
             kcount) -> Tuple[List[List[Piece]], List[float]]:
    """Greedy packing for the Jacobian (2026-10-26, the evaluator: "go ahead with the Jacobian-cost shard balancing").

    A shard's resid+Jacobian time is modelled as ``R (alpha + beta k)``: R its data cost (s, one residual
    evaluation), k its finite-differenced shift columns (each re-evaluates the WHOLE shard, ~beta R).  Units with
    shifts are placed first (largest first), then the rest; each unit goes where the resulting time is smallest, so
    the shifts spread out and the shards holding them stay light.  Returns the bins (file order) and their predicted
    times."""
    bins: List[List[Piece]] = [[] for _ in range(nshards)]
    R = np.zeros(nshards)
    K = np.zeros(nshards)
    cost = lambda u: float(sum(p.cost for p in u))
    order = sorted(units, key=lambda u: (kcount(u) == 0, -cost(u)))
    for u in order:
        r, k = cost(u), kcount(u)
        t = (R + r) * (alpha + beta * (K + k))
        j = int(np.argmin(t))
        bins[j].extend(u)
        R[j] += r
        K[j] += k
    loads = list(R * (alpha + beta * K))
    return ([sorted(b, key=lambda p: (p.base_key, p.e_lo if p.e_lo is not None else 0.0)) for b in bins], loads)


def _build_shards_balanced(base_azr, nshards, workdir, weights, no_split, conv_grid, balance):
    bins, loads = plan_balanced(base_azr, nshards, weights, no_split, conv_grid, balance)
    shards = []
    for j, (b, t) in enumerate(zip(bins, loads)):
        if not b:
            continue
        spec = write_shard(base_azr, f"shard{j:02d}", b, workdir)
        spec.predicted = float(t)
        shards.append(spec)
    return shards


def calibrate_weights(shards: List[ShardSpec], times: Sequence[float]) -> Dict[str, float]:
    """Per-class cost per point from measured shard evaluation times
    (non-negative least squares), normalized so the cheapest class is 1."""
    from scipy.optimize import nnls
    classes = sorted({p.cls for s in shards for p in s.pieces})
    A = np.zeros((len(shards), len(classes)))
    for i, s in enumerate(shards):
        for p in s.pieces:
            A[i, classes.index(p.cls)] += p.npts
    c, _ = nnls(A, np.asarray(times, float))
    c = np.where(c > 0, c, np.nan)
    base = np.nanmin(c) if np.isfinite(c).any() else 1.0
    return {cls: (float(v / base) if np.isfinite(v) else 1.0) for cls, v in zip(classes, c)}


# -- other variants ------------------------------------------------------------



# -- the backend's file work (rmfit.backend) ----------------------------------------------------------------------
from .backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    from ._model_p3 import *  # noqa: E402,F401,F403
    from ._model_p3 import __all__ as _BACKEND_NAMES  # noqa: E402
else:
    from ._model_p2 import *  # noqa: E402,F401,F403
    from ._model_p2 import __all__ as _BACKEND_NAMES  # noqa: E402
