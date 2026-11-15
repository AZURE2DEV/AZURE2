"""Tables and scans that turn a score into something a modeler acts on.

- per-dataset chi-squared (an experiment usually owns several segments);
- the change between two scores, per dataset, with the redistribution ratio
  ``rho = sum(losses) / sum(gains)`` that flags an "improvement" that merely
  moved chi-squared from one experiment to another;
- residual runs: consecutive points more than ``z`` sigma on the same side of
  the fit -- the archive's way of spotting a missing level;
- theta^2 of every free amplitude against the Wigner limit.
"""
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .engine import Evaluator, Score


def dataset_name(seg) -> str:
    return os.path.splitext(os.path.basename(seg.data_file))[0]


def dataset_table(ev: Evaluator, score: Score) -> List[dict]:
    """``[{name, chi2, n, per_pt, segments}]`` sorted by chi2 per point, worst first."""
    by: Dict[str, dict] = {}
    for key, (a, b) in ev.row_ranges:
        seg = ev.segments.by_key(key)
        name = dataset_name(seg)
        d = by.setdefault(name, dict(name=name, chi2=0.0, n=0, segments=set()))
        # a Sharded evaluator lists one row range per PIECE of a split segment, while
        # segment_chi2 is per base segment: count the chi2 once per key, the points per piece
        # (until 2026-09-18 a segment split in two was counted twice)
        if key not in d["segments"]:
            d["chi2"] += float(score.segment_chi2.get(key, 0.0))
        d["n"] += b - a
        d["segments"].add(key)
    out = []
    for d in by.values():
        d["segments"] = sorted(d["segments"])
        d["per_pt"] = d["chi2"] / max(d["n"], 1)
        out.append(d)
    out.sort(key=lambda d: -d["per_pt"])
    return out


def compare(ev: Evaluator, before: Score, after: Score) -> dict:
    """Per-dataset change and the redistribution ratio."""
    tb = {d["name"]: d for d in dataset_table(ev, before)}
    ta = {d["name"]: d for d in dataset_table(ev, after)}
    rows = []
    for name in tb:
        rows.append(dict(name=name, before=tb[name]["chi2"], after=ta[name]["chi2"],
                         delta=ta[name]["chi2"] - tb[name]["chi2"], n=tb[name]["n"]))
    rows.sort(key=lambda r: r["delta"])
    gains = -sum(r["delta"] for r in rows if r["delta"] < 0)
    losses = sum(r["delta"] for r in rows if r["delta"] > 0)
    rho = losses / gains if gains > 0 else (np.inf if losses > 0 else 0.0)
    return dict(rows=rows, gains=gains, losses=losses, rho=rho,
                total_delta=after.objective - before.objective)


@dataclass
class Run:
    key: int
    e_lo: float
    e_hi: float
    n: int
    sign: int
    max_z: float


def residual_runs(energies: np.ndarray, resid: np.ndarray, key: int, z: float = 3.0,
                  min_run: int = 3) -> List[Run]:
    """Consecutive points (in energy order) with |residual| > z on one side."""
    order = np.argsort(energies)
    e, r = np.asarray(energies)[order], np.asarray(resid)[order]
    runs, i = [], 0
    while i < e.size:
        if abs(r[i]) > z:
            s = 1 if r[i] > 0 else -1
            j = i
            while j + 1 < e.size and abs(r[j + 1]) > z and (r[j + 1] > 0) == (s > 0):
                j += 1
            if j - i + 1 >= min_run:
                runs.append(Run(key, float(e[i]), float(e[j]), j - i + 1, s, float(np.max(np.abs(r[i:j + 1])))))
            i = j + 1
        else:
            i += 1
    return runs


def theta2_table(ev: Evaluator, x) -> List[dict]:
    gW = ev.wigner_amplitudes()
    out = []
    for j, k in enumerate(ev.keys):
        if k.kind == "width" and j in gW and gW[j] > 0:
            out.append(dict(key=str(k), level=k.level, pair=k.pair, L=k.L, S=k.S,
                            rwa=float(x[j]), theta2=float((x[j] / gW[j]) ** 2)))
    out.sort(key=lambda d: -d["theta2"])
    return out


def format_dataset_table(rows: Sequence[dict], top: Optional[int] = None) -> str:
    lines = ["| dataset | chi2/N | chi2 | N | segments |", "|---|---:|---:|---:|---|"]
    for d in (rows if top is None else rows[:top]):
        segs = ",".join(str(s) for s in d["segments"][:6]) + ("..." if len(d["segments"]) > 6 else "")
        lines.append(f"| {d['name']} | {d['per_pt']:.2f} | {d['chi2']:.1f} | {d['n']} | {segs} |")
    return "\n".join(lines)


def format_compare(cmp: dict, top: int = 12) -> str:
    lines = [f"total delta {cmp['total_delta']:+.1f}; gains {cmp['gains']:.1f}, losses {cmp['losses']:.1f}, "
             f"rho = {cmp['rho']:.2f}", "| dataset | before | after | delta |", "|---|---:|---:|---:|"]
    rows = cmp["rows"]
    shown = rows[:top // 2] + rows[-(top - top // 2):] if len(rows) > top else rows
    for r in shown:
        lines.append(f"| {r['name']} | {r['before']:.1f} | {r['after']:.1f} | {r['delta']:+.1f} |")
    return "\n".join(lines)


def window_delta(ev: Evaluator, x_with, x_without, e_r: float, gamma_mev: float,
                 multiples: Sequence[float] = (1, 2, 3), min_half_mev: float = 0.0) -> dict:
    """chi2 with and without a resonance inside E_R +- m Gamma (the evaluator 2026-10-12: "give the difference between
    the chi2 with and without the new resonance based on something like +/- 3 times the width of the resonance
    to either side").  Plain data chi2 (no penalty / prior rows, no robust loss), both vectors on the SAME
    evaluator, so a build change between the two fits cannot enter.  Rows are placed by their excitation energy
    (``ev.row_excitation``); rows without one are counted outside every window.

    Returns ``dict(e_r, gamma_mev, total_without, total_with, total_delta, windows=[...], datasets=[...])`` with,
    per window, ``m, half, lo, hi, n, without, with, delta, outside`` and, for the LARGEST window, the per-data-set
    ``name, n, without, with, delta`` sorted by delta (largest gain first)."""
    ex = np.asarray(ev.row_excitation(), float)
    rw, _ = ev.residuals(np.asarray(x_with, float), jac=False)
    ro, _ = ev.residuals(np.asarray(x_without, float), jac=False)
    # point_chi2 (2026-10-22) is absent from an Evaluator loaded before it existed (lazy import after an install)
    _pc = ev.point_chi2 if hasattr(ev, "point_chi2") else (lambda r: ev._plain_data_rows(r) ** 2)
    cw = _pc(rw)
    co = _pc(ro)
    name = np.empty(ev.npts, dtype=object)
    for key, (a, b) in ev.row_ranges:
        name[a:b] = dataset_name(ev.segments.by_key(key))
    dist = np.abs(ex - float(e_r))
    ok = np.isfinite(dist)
    out = dict(e_r=float(e_r), gamma_mev=float(gamma_mev), total_without=float(co.sum()),
               total_with=float(cw.sum()), total_delta=float(cw.sum() - co.sum()), windows=[], datasets=[])
    sel = None
    for m in multiples:
        half = max(float(m) * float(gamma_mev), float(min_half_mev))
        sel = ok & (dist <= half)
        d_in = float(cw[sel].sum() - co[sel].sum())
        out["windows"].append(dict(m=float(m), half=half, lo=float(e_r) - half, hi=float(e_r) + half,
                                   n=int(sel.sum()), without=float(co[sel].sum()), **{"with": float(cw[sel].sum())},
                                   delta=d_in, outside=out["total_delta"] - d_in))
    if sel is not None:
        by: Dict[str, list] = {}
        for i in np.flatnonzero(sel):
            r = by.setdefault(name[i], [0, 0.0, 0.0])
            r[0] += 1; r[1] += co[i]; r[2] += cw[i]
        out["datasets"] = sorted((dict(name=k, n=v[0], without=v[1], **{"with": v[2]}, delta=v[2] - v[1])
                                  for k, v in by.items()), key=lambda r: r["delta"])
    return out


def format_window_delta(res: dict, top: int = 6) -> List[str]:
    """Log lines for :func:`window_delta`: one per window, then the largest gains and losses inside the widest."""
    lines = [f"chi2 with / without the level inside E_R +- m Gamma (E_R {res['e_r']:.4f} MeV, Gamma "
             f"{1e3 * res['gamma_mev']:.1f} keV; plain data chi2, same evaluator):"]
    for w in res["windows"]:
        lines.append(f"  +-{w['m']:g} Gamma  Ex {w['lo']:.3f}-{w['hi']:.3f}  N {w['n']:6d}  without {w['without']:11.1f}  "
                     f"with {w['with']:11.1f}  diff {w['delta']:+9.1f}   outside {w['outside']:+9.1f}")
    ds = res["datasets"]
    if ds:
        m = res["windows"][-1]["m"]
        gains = [r for r in ds if r["delta"] < 0][:top]
        losses = [r for r in reversed(ds) if r["delta"] > 0][:top]
        lines.append(f"  inside +-{m:g} Gamma, by data set (N, without -> with, diff):")
        for r in gains + losses:
            lines.append(f"    {r['name'][:44]:44s} {r['n']:5d} {r['without']:10.1f} -> {r['with']:10.1f}  {r['delta']:+8.1f}")
    return lines
