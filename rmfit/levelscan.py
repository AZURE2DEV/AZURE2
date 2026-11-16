"""Residual-driven level search (plan Phase 7, the evaluator's proposal).

1. ``candidate_energies``: residual runs of the incumbent on every active
   segment, converted to excitation energy and clustered across segments --
   where the fit misses something, ranked by the chi-squared in the runs.
2. ``template_screen``: at one candidate energy, a dummy level of every J^pi
   is appended to the model (all allowed channels, energies fixed).  Each
   trial pattern -- a J^pi, a partition of an estimated total width over the
   group's channels, and a sign pattern -- is written into the amplitude vector
   (rwa from the engine's own physical-width conversion, rescaled by
   sqrt(Gamma_c / Gamma_ref)) and evaluated frozen in the full model, so the
   interference with the existing levels of that J^pi is part of the template;
   the best patterns then get one damped Gauss-Newton step on the level's own
   amplitudes and the norms.  This is the evaluator's "one level of every J^pi at the
   estimated width, compare with the residual" done inside the current model.
   The response to a new level is quadratic in its amplitudes and saturates at
   resonance, so a finite trial width is essential (a linear matched filter is
   identically zero).
3. The best (J^pi, pattern) candidates are handed to the existing level-add
   refit (``Campaign.run_level_add``) with pinned energy.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence
import itertools
import os
import time

import numpy as np

from .diagnose import residual_runs
from .engine import EvalError, Evaluator
from .keys import ParamKey
from .local import BoundsPolicy, default_bounds, gn_step
from .structure import JPI_GROUPS, LevelCandidate, jpi_parse, lowest_L, resolve_jpis  # noqa: F401 (JPI_GROUPS: explicit use only)


@dataclass
class ExCandidate:
    ex: float                     # excitation energy (MeV), chi2-weighted centre of the cluster
    score: float                  # chi2 carried by the residual runs of the cluster
    npts: int
    segments: List[int]
    runs: List[dict] = field(default_factory=list)
    width: float = 0.0            # MeV: width of the structure in the DIFFERENTIAL runs (see cluster_width)
    span: float = 0.0             # MeV: full extent of the cluster, lowest e_lo to highest e_hi
    width_total: float = 0.0      # MeV: the same statistic over the total/integrated-segment runs
    density: float = 0.0          # chi2 per MeV of structure (cluster_density), the rank="density" key

    def to_json(self):
        return dict(ex=self.ex, score=self.score, npts=self.npts, segments=self.segments, runs=self.runs,
                    width=self.width, span=self.span, width_total=self.width_total, density=self.density)


def _weighted_median(values, weights) -> float:
    v = np.asarray(values, float); w = np.asarray(weights, float)
    if v.size == 0:
        return 0.0
    o = np.argsort(v); v, w = v[o], w[o]
    c = np.cumsum(w)
    return float(v[int(np.searchsorted(c, 0.5 * c[-1]))])


def cluster_width(runs: Sequence[dict], min_energies: int = 3, max_cover: float = 0.5) -> dict:
    """The width the residual structure itself shows, before any template is
    tried (2026-09-24: a high-energy window turned out to hold
    600 keV - 2 MeV levels while the screen sampled 1.5-50 keV and sat at the
    top of its grid 95 times out of 120).

    ``width``: chi2-weighted MEDIAN of the energy extents (e_hi - e_lo) of the
    runs in differential segments that cover at least ``min_energies``
    distinct energies -- a proxy for the FWHM of a peak or one lobe of a
    dispersion shape.  A median, not a mean: on the first try the mean was the
    Cierjacks n-total alone, one 791-point run over 560 keV (a broadband
    deficit, not a resonance), and angular distributions at a single energy
    give runs of zero extent.  A run covering more than ``max_cover`` of its
    segment's energy range (``cover``) is a normalization or broadband deficit
    -- the whole-range Ge secondary-gamma runs at 12.5-12.9 -- and does not
    count either.  Runs without ``diff``/``ne``/``cover`` fields (older
    ledgers) count as differential, over enough energies, not broadband.
    Falls back to every eligible run when no differential one qualifies.
    ``width_total``: the same over the non-differential (total, integrated) runs.
    ``span``: full extent of the cluster.  All MeV; 0 when nothing qualifies."""
    runs = list(runs or [])
    if not runs:
        return dict(width=0.0, span=0.0, width_total=0.0, n_runs=0)

    def ext(r):
        return max(float(r["e_hi"]) - float(r["e_lo"]), 0.0)

    def wt(r):
        return max(float(r.get("chi2", 1.0)), 1e-30)

    ok = [r for r in runs if int(r.get("ne", min_energies)) >= min_energies and ext(r) > 0.0
          and float(r.get("cover", 0.0)) <= max_cover]
    diff = [r for r in ok if bool(r.get("diff", True))]
    tot = [r for r in ok if not bool(r.get("diff", True))]
    width = _weighted_median([ext(r) for r in diff], [wt(r) for r in diff]) if diff else \
        _weighted_median([ext(r) for r in ok], [wt(r) for r in ok])
    width_total = _weighted_median([ext(r) for r in tot], [wt(r) for r in tot]) if tot else 0.0
    span = float(max(r["e_hi"] for r in runs) - min(r["e_lo"] for r in runs))
    return dict(width=width, span=max(span, 0.0), width_total=width_total, n_runs=len(runs))


def auto_gamma_kev(cand: dict, bounds_kev: Sequence[float] = (10.0, 1000.0), default_kev: float = 10.0) -> float:
    """Template total width (keV) for a screen at this cluster: the cluster's
    run width clamped into ``bounds_kev``; ``default_kev`` when the candidate
    carries no run width (an ``--ex`` list, or single-point runs)."""
    w = float(cand.get("width", 0.0) or 0.0) * 1e3
    if not np.isfinite(w) or w <= 0.0:
        return float(default_kev)
    lo, hi = float(bounds_kev[0]), float(bounds_kev[1])
    return float(min(max(w, lo), hi))


def candidate_energies(ev: Evaluator, x, z: float = 3.0, min_run: int = 3, bin_width: float = 0.03,
                       top: int = 20, exclude_segments: Sequence[int] = (), pair_width: float = 0.15,
                       include_match: Optional[str] = None, rank: str = "chi2",
                       density_floor_kev: float = 20.0) -> List[ExCandidate]:
    """Residual peaks of ``x`` clustered in excitation energy.  ``ev`` must be a
    single session (its ``data()`` gives the excitation energy of every point
    in row order).

    2026-09-24 (the evaluator: "there are many narrow structures, especially in the
    secondary gamma ray cross sections, but you failed to fit any of them"):
    ``include_match`` -- a regex on the segment's data-file name; only matching
    segments contribute runs (e.g. ``NoFeed`` = the secondary-gamma sets, 12-125
    points each, whose runs never outrank a 791-point n-total run).
    ``rank="density"`` -- order clusters by chi2 per MeV of structure,
    score / max(width or span, ``density_floor_kev``), so a sharp feature
    outranks a smooth deficit of the same chi2; ``"chi2"`` is the old order."""
    r, _ = ev.residuals(x, jac=False)
    # point_residuals (2026-10-22) is absent from an Evaluator loaded before it existed: a running process
    # that imports this module lazily after an install must still work
    plain = ev.point_residuals(r) if hasattr(ev, "point_residuals") else ev._plain_data_rows(r)
    d = ev.data()
    exc = d["excitation"]
    runs = []
    import re as _re
    pat = _re.compile(include_match) if include_match else None
    for i, (key, (a, b)) in enumerate(ev.row_ranges):
        if key in exclude_segments or i >= len(exc):
            continue
        if pat is not None:
            try:
                fname = os.path.basename(str(getattr(ev.segments.by_key(key), "data_file", "") or ""))
            except Exception:
                fname = ""
            if not pat.search(fname):
                continue
        Ex = np.asarray(exc[i], float)
        res = plain[a:b]
        if Ex.size != res.size:
            continue
        try:   # differential (excitation function / angular distribution) or total, for cluster_width
            is_diff = str(getattr(ev.segments.by_key(key), "observable", "")).lower() in DIFFERENTIAL_OBSERVABLES
        except Exception:
            is_diff = True
        for run in residual_runs(Ex, res, key, z=z, min_run=min_run):
            sel = (Ex >= run.e_lo - 1e-12) & (Ex <= run.e_hi + 1e-12)
            chi = float(np.sum(res[sel] ** 2))
            centre = float(np.sum(Ex[sel] * res[sel] ** 2) / max(chi, 1e-30))
            runs.append(dict(segment=int(key), e_lo=run.e_lo, e_hi=run.e_hi, n=int(sel.sum()),
                             sign=int(run.sign), zmax=float(run.max_z), chi2=chi, centre=centre,
                             ne=int(np.unique(np.round(Ex[sel], 6)).size), diff=bool(is_diff),
                             cover=float((run.e_hi - run.e_lo) / max(float(Ex.max() - Ex.min()), 1e-12))))
    runs.sort(key=lambda t: t["centre"])
    clusters: List[dict] = []
    for run in runs:
        if clusters and abs(run["centre"] - clusters[-1]["centre"]) <= bin_width:
            c = clusters[-1]
            w = c["chi2"] + run["chi2"]
            c["centre"] = (c["centre"] * c["chi2"] + run["centre"] * run["chi2"]) / max(w, 1e-30)
            c["chi2"] = w; c["npts"] += run["n"]; c["runs"].append(run)
        else:
            clusters.append(dict(centre=run["centre"], chi2=run["chi2"], npts=run["n"], runs=[run]))
    out = [ExCandidate(ex=round(c["centre"], 6), score=c["chi2"], npts=c["npts"],
                       segments=sorted({r["segment"] for r in c["runs"]}), runs=c["runs"],
                       **{k: v for k, v in cluster_width(c["runs"]).items() if k != "n_runs"}) for c in clusters]
    # a level interfering with its neighbours leaves a dispersion-shaped
    # residual: two adjacent clusters of opposite sign with the level between
    # them -- add a candidate at their chi2-weighted boundary
    def dominant_sign(c):
        return 1 if sum(r["chi2"] * r["sign"] for r in c["runs"]) > 0 else -1
    for a, b in zip(clusters, clusters[1:]):
        if dominant_sign(a) != dominant_sign(b) and 0 < b["centre"] - a["centre"] <= pair_width:
            hi_a = max(r["e_hi"] for r in a["runs"]); lo_b = min(r["e_lo"] for r in b["runs"])
            out.append(ExCandidate(ex=round(0.5 * (hi_a + lo_b), 6), score=a["chi2"] + b["chi2"], npts=a["npts"] + b["npts"],
                                   segments=sorted({r["segment"] for r in a["runs"] + b["runs"]}),
                                   runs=a["runs"] + b["runs"],
                                   **{k: v for k, v in cluster_width(a["runs"] + b["runs"]).items() if k != "n_runs"}))
    for c in out:
        c.density = cluster_density(c.score, c.width, c.span, density_floor_kev)
    if rank == "density":
        out.sort(key=lambda c: -c.density)
    elif rank == "chi2":
        out.sort(key=lambda c: -c.score)
    else:
        raise ValueError(f"rank must be 'chi2' or 'density', not {rank!r}")
    return out[:top]


def cluster_density(score: float, width: float, span: float, floor_kev: float = 20.0) -> float:
    """chi2 per MeV of structure: ``score`` / max(``width`` (or ``span`` when the
    width is 0), ``floor_kev``).  The floor keeps a zero-extent or single-energy
    cluster from ranking infinitely high."""
    ext = width if width and width > 0 else (span or 0.0)
    return float(score) / max(float(ext), float(floor_kev) * 1e-3)


def template_variant(base_azr: str, path: str, ex: float, jpis: Optional[Sequence[str]] = None,
                     gamma_ref_ev: float = 1.0, gamma_ref_photon_ev: float = 1e-3,
                     output_dir: Optional[str] = None) -> dict:
    """``path`` = ``base_azr`` plus one dummy level per J^pi at ``ex`` (energy
    free, so the screen can offset it) carrying every channel the group already uses (or the lowest waves
    into pairs 1 and 2 for a new group), each at ``gamma_ref_ev`` so that the
    engine's initial vector gives the rwa of that physical width.
    ``jpis`` None/empty = the model's own groups; a J^pi no particle pair can couple to
    raises (structure.resolve_jpis, 2026-11-04 -- the old default was the half-integer
    list (17O-like), unphysical on integer-spin systems)."""
    from pyazr import AzrModel
    # AZURE2's level-line reader mis-parses energies written with ~16 significant
    # digits (13N: a 17-character energy turned the model into garbage); keep 1 eV
    ex = round(float(ex), 6)
    mdl = AzrModel.from_file(base_azr)
    jpis = resolve_jpis(mdl, jpis)
    levels = {}
    for jpi in jpis:
        J, par = jpi_parse(jpi)
        existing = mdl.find(jpi=jpi)
        photon = []
        if existing:
            seen, chans = set(), []
            for lv in existing:
                for c in lv.channels:
                    key = (c.pair, c.L, round(float(c.S), 3))
                    if key in seen:
                        continue
                    seen.add(key)
                    isph = bool(getattr(c, "is_photon", False))
                    chans.append(dict(pair=c.pair, L=c.L, S=c.S, gamma=(gamma_ref_photon_ev if isph else gamma_ref_ev),
                                      fixed=False))
                    photon.append(isph)
        else:
            chans = []
            for pair in (1, 2, 3):
                try:
                    c = mdl._pair_template(pair)
                except ValueError:
                    continue
                if getattr(c, "is_photon", False):
                    continue                              # a new group gets particle channels only
                for L, S in lowest_L(J, par, c.J1, c.parity1, c.J2, c.parity2)[:1]:
                    chans.append(dict(pair=pair, L=L, S=S, gamma=gamma_ref_ev, fixed=False))
                    photon.append(False)
        if not chans or not any(not ph for ph in photon):
            continue
        # the partition of each EXISTING level of this group, on the dummy's own
        # particle-channel order: the only real prior on how the group spreads width
        part_order = [(c["pair"], c["L"], round(float(c["S"]), 3)) for c, ph in zip(chans, photon) if not ph]
        priors = []
        for lv in existing:
            gm = {(c.pair, c.L, round(float(c.S), 3)): float(c.gamma) for c in lv.channels}
            v = np.array([gm.get(k, 0.0) for k in part_order], float)
            if np.abs(v).sum() <= 0:
                continue
            fr, sg = partition_fractions(v)
            priors.append(dict(label=f"{jpi}@{lv.energy:.3f}", fractions=[float(t) for t in fr],
                               signs=[int(t) for t in sg]))
        mdl.add_level(J=J, parity=par, energy=ex, channels=chans, level_fixed=False)   # energy free: the screen moves it
        levels[jpi] = dict(label=f"{jpi}#{len(existing) + 1}", channels=chans, photon=photon,
                           priors=priors)   # label: file order; the engine's is resolved by dummy_labels
    if output_dir is not None:
        mdl.set_output_dir(output_dir)
    # the variant's <segmentsData> paths are relative ("data/..."), so the directory
    # we write it into needs the same data/ link the base model has beside it --
    # same idiom as structure.window_screen and windows.ladder_step
    wdir = os.path.dirname(os.path.abspath(path))
    os.makedirs(wdir, exist_ok=True)
    from .model import link_data_dirs
    link_data_dirs(base_azr, wdir)
    mdl.write(path)
    return dict(path=path, ex=ex, levels=levels, gamma_ref_ev=gamma_ref_ev, gamma_ref_photon_ev=gamma_ref_photon_ev)


def dummy_labels(ev: Evaluator, variant: dict) -> Dict[str, str]:
    """``{jpi: level label}`` of the dummy levels as the engine numbers them
    (levels are numbered inside a J^pi group in the engine's order, which is
    not the file's): the level of that J^pi whose seed energy is ``ex``."""
    out = {}
    with_widths = {k.level for k in ev.keys if k.kind == "width"}
    for jpi in variant["levels"]:
        hits = [lv for lv, e in ev.level_energy.items() if lv.startswith(jpi + "#") and abs(float(e) - variant["ex"]) < 1e-6]
        # a level of the model that was "removed" by zeroing its widths may sit at
        # the same energy: the dummy is the one whose widths are free parameters
        hits = [lv for lv in hits if lv in with_widths] or hits
        if len(hits) == 1:
            out[jpi] = hits[0]
    return out


def incumbent_on_variant(ev: Evaluator, x, keys, variant: dict) -> np.ndarray:
    """The incumbent mapped onto the template variant's keys: the variant's own
    initial values for the new keys (dummy energies at ``ex``), then the dummy
    amplitudes set to zero."""
    from .keys import take
    labels = set(dummy_labels(ev, variant).values())
    xv = take(x, keys, ev.keys, ev.x0)
    for j, k in enumerate(ev.keys):
        if k.kind == "width" and k.level in labels:
            xv[j] = 0.0
    return xv


def partition_fractions(widths):
    """``(fractions, signs)`` of a real level's channel widths: the share of
    the total |Gamma| each channel carries, and its sign."""
    w = np.asarray(widths, float)
    tot = float(np.abs(w).sum())
    if tot <= 0:
        return None, None
    return np.abs(w) / tot, np.where(w >= 0, 1, -1).astype(int)


def pattern_distance(f_true, s_true, f_try, s_try, floor: float = 0.02) -> float:
    """How far a template pattern is from a real level's partition: ``0`` means
    the same shape and the same relative signs, ``1`` orthogonal shapes.  A
    global sign flip is an exact symmetry of a level, so only the RELATIVE
    signs of the channels that carry strength are compared."""
    ft, fy = np.asarray(f_true, float), np.asarray(f_try, float)
    nt, ny = np.linalg.norm(ft), np.linalg.norm(fy)
    shape = 1.0 - float(ft @ fy) / max(nt * ny, 1e-30)
    st, sy = np.asarray(s_true, int), np.asarray(s_try, int)
    m = (ft > floor) & (np.abs(fy) > floor)
    if not m.any():
        return shape + 1.0
    agree = float(np.mean(st[m] == sy[m]))
    return shape + (1.0 - max(agree, 1.0 - agree))


def grid_coverage(widths_true, patterns, totals):
    """Can the template grid represent this level AT ALL?  Pure geometry, no
    evaluation: returns the closest grid point to ``widths_true`` and its
    distance, plus the total-width ratio.  A screen cannot find a level whose
    closest grid point is far away, however many patterns it evaluates -- run
    this before blaming the search."""
    f_true, s_true = partition_fractions(widths_true)
    if f_true is None:
        raise ValueError("the level carries no width")
    gt_true = float(np.abs(np.asarray(widths_true, float)).sum())
    best = None
    for fr, sg, desc in patterns:
        d = pattern_distance(f_true, s_true, fr, sg)
        for gt in totals:
            ratio = gt / gt_true
            row = dict(distance=d, desc=desc, fractions=list(fr), signs=list(sg),
                       gamma_total_ev=float(gt), width_ratio=ratio,
                       penalty=d + abs(np.log(max(ratio, 1e-12))))
            if best is None or row["penalty"] < best["penalty"]:
                best = row
    best.update(gamma_true_ev=gt_true, fractions_true=[float(v) for v in f_true],
                signs_true=[int(v) for v in s_true])
    return best


def segment_ex_ranges(ev) -> Dict[int, tuple]:
    """``{segment key: (Ex_lo, Ex_hi)}`` in MeV, from a session that can report
    the excitation energy of every point (a single session).  Used to pick the
    segments a candidate energy can actually affect."""
    d = ev.data()
    exc = d["excitation"]
    out = {}
    for i, (key, _) in enumerate(ev.row_ranges):
        if i >= len(exc):
            continue
        E = np.asarray(exc[i], float)
        if E.size:
            out[int(key)] = (float(E.min()), float(E.max()))
    return out


DIFFERENTIAL_OBSERVABLES = ("differential", "differential-cm")


def scoring_segments(ev, ex: float, half_width: float = 0.25, differential_only: bool = True,
                     ex_ranges: Optional[Dict[int, tuple]] = None) -> List[int]:
    """The segments a level at ``ex`` should be judged on: those whose energy
    range comes within ``half_width`` MeV of it and -- by default -- only the
    DIFFERENTIAL ones, because angle-integrated and total cross sections barely
    distinguish J^pi, so a wrong-J^pi level that adds strength at the right
    energy scores as well as the right one on them.  Judging on the angular
    data is what the method was supposed to do."""
    ranges = ex_ranges if ex_ranges is not None else segment_ex_ranges(ev)
    keys = []
    for key, (lo, hi) in ranges.items():
        if hi < ex - half_width or lo > ex + half_width:
            continue
        if differential_only:
            seg = ev.segments.by_key(key)
            # "differential-cm" (isDiff 4: secondary gamma rays, CM-frame data) is differential
            # data too; until 2026-09-18 it was silently left out of the local ranking
            if str(getattr(seg, "observable", "")).lower() not in DIFFERENTIAL_OBSERVABLES:
                continue
        keys.append(int(key))
    return sorted(keys)


def _patterns(n: int, max_signs: int, seed: int = 0, priors: Sequence[dict] = (),
              dominant: Sequence[float] = (0.5, 0.8), dominant_top: int = 3,
              weights: Sequence[float] = (), pair_top: int = 4, pair_splits: Sequence[float] = (0.5, 0.7)):
    """``(fractions, signs, description)`` templates for a level with ``n``
    particle channels.

    Four families, in order of how much they know:

    * ``prior``  -- the partition of an EXISTING level of the same J^pi in this
      model, which is the only real information available about how that group
      distributes strength.  This is the family that matters: a single-channel
      or equal-split guess cannot represent a level spread unevenly over four
      or five channels, and that is exactly why an accepted
      level (five channels, 40:1 spread, mixed signs) scored
      nothing in the first gate run.
    * ``dominant`` -- one channel takes a fraction ``f``, the rest share the
      remainder; covers a level whose strength sits mostly in one wave.
    * ``channel c only`` -- the single-channel limit.
    * ``equal split`` -- the uninformative default.
    * ``pair`` (2026-09-21) -- TWO channels share everything (50/50, 70/30, 30/70, both relative
      signs), over the pairs among the ``pair_top`` channels with the largest ``weights`` (the
      Wigner widths of the carrying channels, i.e. the channels that can physically carry
      width) plus any channel a prior puts >= 10% on.  In one J^pi study a level carried
      44% n0 g-wave + 56% n2 p-wave and nothing else -- no single-channel, dominant-plus-tail
      or prior template came close (the nearest prior put 23% into alpha0 L5, which the elastic
      data reject), the screen scored the level's own J^pi at +2,724 of a 33,979 removal while
      its own partition dropped in recovers 99%.  A level that lives in the PRODUCT of two
      channels needs both in one template.
    """
    pats = []
    for pr in priors:
        f = np.abs(np.asarray(pr["fractions"], float))[:n]
        if f.size < n:
            f = np.concatenate([f, np.zeros(n - f.size)])
        tot = float(f.sum())
        if tot <= 0:
            continue
        sg = np.asarray(pr["signs"], int)[:n]
        if sg.size < n:
            sg = np.concatenate([sg, np.ones(n - sg.size, int)])
        pats.append((tuple(f / tot), tuple(int(v) for v in sg), f"prior {pr.get('label', '?')}"))
        pats.append((tuple(f / tot), tuple(-int(v) for v in sg), f"prior {pr.get('label', '?')} flipped"))
    if n >= 2:
        # Only the channels the group's existing levels actually use strongly are
        # worth a dominant-channel template, and at most `dominant_top` of them:
        # one J^pi group of a large model has 19 channels, and 4n patterns per scale per group is how
        # a screen ends up at 77 min an energy.
        if priors:
            strength = np.zeros(n)
            for pr in priors:
                v = np.abs(np.asarray(pr["fractions"], float))[:n]
                strength[:v.size] += v
            cand = list(np.argsort(-strength)[:dominant_top])
        else:
            cand = list(range(min(n, dominant_top)))
        for c in cand:
            for fdom in dominant:
                f = [(1.0 - fdom) / (n - 1)] * n
                f[c] = fdom
                base = [1] * n
                for sflip in (None, c):
                    sg = list(base)
                    if sflip is not None:
                        sg[sflip] = -1
                    pats.append((tuple(f), tuple(sg), f"channel {c} dominant {fdom:.0%}"))
    for c in range(n):
        f = [0.0] * n; f[c] = 1.0
        for s in (1, -1):
            sg = [0] * n; sg[c] = s
            pats.append((tuple(f), tuple(sg), f"channel {c} only"))
    if n >= 2 and pair_top >= 2:
        w = np.asarray(list(weights) + [0.0] * max(0, n - len(weights)), float)[:n]
        cand = list(np.argsort(-w)[:pair_top]) if np.any(w > 0) else list(range(min(n, pair_top)))
        for pr in priors:
            v = np.abs(np.asarray(pr["fractions"], float))[:n]
            cand += [int(i) for i in np.flatnonzero(v >= 0.1) if int(i) not in cand]
        cand = sorted(set(int(c) for c in cand))[:pair_top + 1]
        for a, b in itertools.combinations(cand, 2):
            fracs = set()
            for fs in pair_splits:
                fracs.add((round(float(fs), 3), round(1.0 - float(fs), 3))); fracs.add((round(1.0 - float(fs), 3), round(float(fs), 3)))
            for fa, fb in sorted(fracs, reverse=True):
                f = [0.0] * n; f[a] = fa; f[b] = fb
                for sb in (1, -1):
                    sg = [0] * n; sg[a] = 1; sg[b] = sb
                    pats.append((tuple(f), tuple(sg), f"pair {a}+{b} {fa:.0%}/{fb:.0%}"))
    if n >= 2:
        # 2026-11-06: never materialise all 2^n sign vectors -- at n = 30-31 channels the
        # list is 300-650 GB and the job is OOM-killed.  Draw the
        # same indices as before (rng.choice on the integer 2^n) and decode index i into itertools.product((1, -1),
        # repeat=n) order (bit n-1-k of i set -> sign -1 at position k): identical templates, O(max_signs) memory.
        nall = 2 ** n
        if nall > max_signs:
            rng = np.random.default_rng(seed)
            idx = [int(i) for i in sorted(rng.choice(nall, max_signs, replace=False))]
        else:
            idx = list(range(nall))
        allsigns = [tuple(-1 if (i >> (n - 1 - k)) & 1 else 1 for k in range(n)) for i in idx]
        for sg in allsigns:
            pats.append((tuple([1.0 / n] * n), tuple(sg), "equal split"))
    return pats


def align_priors(priors: Sequence[dict], chan_order: Sequence[tuple], carrying: Sequence[tuple]) -> List[dict]:
    """Re-express each prior partition on exactly the channels a template uses.

    ``priors[i]["fractions"/"signs"]`` run over ``chan_order`` (the dummy's particle channels in
    file order, as ``template_variant`` wrote them); ``carrying`` is the (pair, L, S) of the channels
    that enter the templates, in template order.  ``_patterns`` used to take ``fractions[:n]``,
    which puts a prior's strength on the WRONG channels whenever a closed or dead channel sits
    before an open one (2026-09-20).  Priors with no strength left on the carrying channels are
    dropped."""
    out = []
    for pr in priors:
        m = {tuple(k): (abs(float(f)), int(sg)) for k, f, sg in zip(chan_order, pr.get("fractions", ()), pr.get("signs", ()))}
        fr = [m.get(tuple(k), (0.0, 1))[0] for k in carrying]
        if sum(fr) <= 0:
            continue
        out.append(dict(label=pr.get("label", "?"), fractions=fr, signs=[m.get(tuple(k), (0.0, 1))[1] for k in carrying]))
    return out


def carrying_channels(gamma_w_ev: Sequence[float], gamma_w_amp: Sequence[float], dead_channel_ev: float = 1.0):
    """Which particle channels of a dummy level may enter a template.

    ``gamma_w_ev[i]`` / ``gamma_w_amp[i]``: the channel's Wigner-limit width (eV, at the level
    energy) and amplitude; 0 for a closed channel.  Returns ``(open_idx, dead_idx, closed_idx)``:
    *open* = Gamma_W >= ``dead_channel_ev`` (the rule ``run_level_add`` uses for its dead-channel
    zeroing), *dead* = open in principle but unable to hold ``dead_channel_ev`` even at the Wigner
    limit (high L, barely above threshold), *closed* = no Wigner width at all.

    2026-09-20: until now every channel with Gamma_W > 0 entered the templates.  On one model an
    "equal split" of 10 keV over 15 channels put 667 eV into five channels with Gamma_W of 0.4 eV
    .. 1e-15 eV (theta^2 1e3 .. 1e17); the screen evaluates templates frozen and uncapped, and
    reported a 7/2+ worth +1,693 that the honest add could not find (+28)."""
    open_idx, dead_idx, closed_idx = [], [], []
    for i, (gw, ga) in enumerate(zip(gamma_w_ev, gamma_w_amp)):
        if not (gw > 0 and ga > 0):
            closed_idx.append(i)
        elif gw < dead_channel_ev:
            dead_idx.append(i)
        else:
            open_idx.append(i)
    return open_idx, dead_idx, closed_idx


def template_amplitude(width_ev: float, gamma_w_ev: float, gamma_w_amp: float, theta2_cap: float):
    """``(|amplitude|, |width_ev|)`` of one template channel asked to carry ``width_ev``:
    gamma = gamma_W sqrt(Gamma / Gamma_W), CLIPPED at theta^2 = ``theta2_cap`` (a template must be
    a level the bounded refit could actually reach)."""
    width_ev = abs(float(width_ev))
    if not (gamma_w_ev > 0 and gamma_w_amp > 0) or width_ev == 0.0:
        return 0.0, 0.0
    theta2 = min(width_ev / gamma_w_ev, float(theta2_cap))
    return float(gamma_w_amp * np.sqrt(theta2)), float(theta2 * gamma_w_ev)


def template_screen(ev: Evaluator, x_inc, variant: dict, gamma_total_ev: float,
                    scales: Sequence[float] = (0.3, 0.6, 1.2, 2.5, 5.0, 10.0),
                    max_signs: int = 32, gn_top: int = 20, gn_per_jpi: int = 2,
                    policy: BoundsPolicy = BoundsPolicy(),
                    jpis: Optional[Sequence[str]] = None, ex_offsets: Sequence[float] = (-0.04, 0.0, 0.04),
                    gamma_photon_ev: float = 1.0, score_keys: Optional[Sequence[int]] = None,
                    stage_scales: Optional[Sequence[float]] = (1.2, 5.0), refine_top: int = 3,
                    dead_channel_ev: float = 1.0, log=print) -> List[dict]:
    """Frozen evaluation of every (J^pi, width partition, signs, scale)
    pattern of the dummy levels in ``variant`` on the sharded evaluator of the
    variant, then a damped GN step on the level's amplitudes plus the norms for
    the ``gn_top`` best.  ``x_inc`` is the incumbent mapped onto ``ev.keys``
    (dummy-level amplitudes zero).  Returns rows sorted by the GN objective.

    Two stages (2026-09-14, after one candidate cost tens of thousands of patterns and
    4 h): every partition x sign pattern is first evaluated at ``stage_scales``
    only and at the candidate energy itself; the ``refine_top`` best SHAPES per
    J^pi are then evaluated at every scale in ``scales`` and every energy
    offset.  Shape ranking is nearly independent of the trial width and of a
    40 keV energy offset, and the width/offset refinement is cheap once the
    shapes are few.  ``stage_scales=None`` (or ``refine_top=0``) restores the
    single exhaustive stage."""
    x_inc = np.asarray(x_inc, float)
    pos = {k: i for i, k in enumerate(ev.keys)}
    norm_cols = [i for i, k in enumerate(ev.keys) if k.kind == "norm"]
    bounds = default_bounds(ev, x_inc, policy)
    t0 = time.time()
    obj0 = ev.objective(x_inc)
    log(f"  template screen at Ex {variant['ex']:.4f} MeV: incumbent {obj0:.1f}; Gamma_total {gamma_total_ev:.0f} eV"
        f" x {list(scales)}")
    rows = []
    labels = dummy_labels(ev, variant)
    gW_amp = ev.wigner_amplitudes()                       # {col: gamma_W}  (rwa at the Wigner limit)
    gW_wid = ev.wigner_widths(x_inc)                      # {col: Gamma_W eV} at the dummies' energies; open channels only
    two_stage = bool(stage_scales) and refine_top > 0 and (len(scales) > 1 or len(ex_offsets) > 1)
    if two_stage:
        s1 = [sc for sc in scales if any(abs(sc - t) < 1e-9 for t in stage_scales)]
        if not s1:  # none of the requested stage scales is in the grid: take the two nearest 1x
            s1 = sorted(scales, key=lambda sc: abs(np.log(sc)))[:2]
        log(f"  two-stage screen: every shape at scales {s1} and the candidate energy, then the best {refine_top} "
            f"shapes per J^pi at all {len(scales)} scales x {len(ex_offsets)} offsets")
    else:
        s1 = list(scales)
    gref_p = float(variant.get("gamma_ref_ev", 1.0)); gref_ph = float(variant.get("gamma_ref_photon_ev", 1e-3))
    for jpi, info in variant["levels"].items():
        if jpis and jpi not in jpis:
            continue
        label = labels.get(jpi)
        if label is None:
            log(f"    {jpi}: dummy level not found in the engine's key table; skipped"); continue
        amp_keys = [k for k in ev.keys if k.kind == "width" and k.level == label]
        e_keys = [k for k in ev.keys if k.kind == "energy" and k.level == label]
        if not amp_keys:
            continue
        cols = np.array([pos[k] for k in amp_keys])
        ecol = pos[e_keys[0]] if e_keys else None
        # which of the level's channels are photon channels (by the variant's channel order)
        ch_specs = info["channels"]; ph_flags = info.get("photon", [False] * len(ch_specs))
        def is_photon_key(k):
            for spec, ph in zip(ch_specs, ph_flags):
                if int(spec["pair"]) == k.pair and int(spec["L"]) == k.L and abs(float(spec["S"]) - k.S) < 1e-6:
                    return ph
            return False
        photon_idx = [i for i, k in enumerate(amp_keys) if is_photon_key(k)]
        # Amplitude for a physical width Gamma_c.  Until 2026-09-14 this was
        # |x0| sqrt(Gamma_c / Gamma_ref): the rwa the engine derived from the
        # variant's 1 eV reference width.  In a dense model that derivation goes
        # through a Brune transform that FAILS for a 1 eV level inside a group of
        # keV-wide neighbours ("Denominator less than zero", once per dummy) and
        # returns an rwa ~2,400x too small -- every template of every group was
        # evaluated at an amplitude ~0, which is why two screens said null at
        # energies where `structure add` then found +2,266 (1/2+ 9.862) and the
        # accepted 5/2+ 10.55 scored nothing.  Now: Gamma_W (eV, at the dummy's
        # energy) and gamma_W (the rwa at the Wigner limit) from the engine,
        # gamma = gamma_W sqrt(Gamma_c / Gamma_W), no transform of the dummy
        # involved; closed channels (no Gamma_W: an ANC, not a width) are left out.
        ref = np.zeros(len(amp_keys))
        gw_ev = [0.0 if i in photon_idx else float(gW_wid.get(int(cols[i]), 0.0)) for i in range(len(amp_keys))]
        gw_am = [0.0 if i in photon_idx else float(gW_amp.get(int(cols[i]), 0.0)) for i in range(len(amp_keys))]
        for i in photon_idx:
            ref[i] = abs(float(ev.x0[cols[i]]))              # photon: linear, the file's value stands
        open_idx, dead_idx, _closed_idx = carrying_channels(gw_ev, gw_am, dead_channel_ev)
        open_idx = [i for i in open_idx if i not in photon_idx]
        dead_idx = [i for i in dead_idx if i not in photon_idx]
        for i in open_idx:
            ref[i] = gw_am[i] / np.sqrt(gw_ev[i])
        t2cap = float(policy.cap_for(variant["ex"]))
        part_idx = open_idx
        n = len(part_idx)
        closed = [k for i, k in enumerate(amp_keys) if i not in photon_idx and i not in open_idx and i not in dead_idx]
        if dead_idx:
            log(f"    {jpi}: {len(dead_idx)} channel(s) with Gamma_W < {dead_channel_ev:g} eV left out of the templates: "
                + ", ".join(f"{str(amp_keys[i])} ({gw_ev[i]:.2g} eV)" for i in dead_idx))
        x0_ratio = [abs(float(ev.x0[cols[i]])) / np.sqrt(gref_p) / ref[i] for i in part_idx if ref[i] > 0]
        log(f"    {jpi}: {n} open particle channels" + (f", {len(closed)} closed left out" if closed else "")
            + (f"; x0-based / Wigner-based reference amplitude {min(x0_ratio):.2g}..{max(x0_ratio):.2g}"
               + (" -- the file-derived rwa is unusable here (Brune transform of the dummy)" if min(x0_ratio) < 0.3 or max(x0_ratio) > 3 else "")
               if x0_ratio else ""))
        if n == 0:
            continue
        ph_sign_sets = [()] if not photon_idx else (list(itertools.product((1, -1), repeat=len(photon_idx)))
                                                   if len(photon_idx) <= 2 else [(1,) * len(photon_idx), (-1,) * len(photon_idx)])
        _order = [(c["pair"], c["L"], round(float(c["S"]), 3)) for c, ph in zip(info.get("channels", ()), info.get("photon", ())) if not ph]
        _carry = [(amp_keys[i].pair, amp_keys[i].L, round(float(amp_keys[i].S), 3)) for i in part_idx]
        pats = _patterns(n, max_signs, priors=align_priors(info.get('priors', ()), _order, _carry),
                         weights=[gw_ev[i] for i in part_idx])
        shapes = [(fr, sg, desc, phs) for fr, sg, desc in pats for phs in ph_sign_sets]
        offsets = list(ex_offsets) if ecol is not None else [0.0]
        done = set()
        group_rows = []

        def evaluate(fr, sg, desc, phs, scale, de):
            key = (desc, tuple(sg), tuple(phs), float(scale), float(de))
            if key in done:
                return
            done.add(key)
            gtot = gamma_total_ev * scale
            gph = gamma_photon_ev * scale
            amps = np.zeros(len(amp_keys)); widths = [0.0] * len(amp_keys)
            for f, sgn, i in zip(fr, sg, part_idx):
                a_, w_ = template_amplitude(f * gtot, gw_ev[i], gw_am[i], t2cap)      # clipped at the theta^2 cap
                amps[i] = sgn * a_; widths[i] = float(sgn * w_)
            for sgn, i in zip(phs, photon_idx):
                amps[i] = sgn * ref[i] * np.sqrt(gph / gref_ph); widths[i] = float(sgn * gph)
            full_desc = desc + (f", photon signs {list(phs)}" if phs else "")
            x = x_inc.copy(); x[cols] = amps
            if ecol is not None:
                x[ecol] = variant["ex"] + de
            x = bounds.clip(x)
            try:
                of = ev.objective(x)
            except EvalError as err:
                log(f"    {jpi} {full_desc} {sg}: engine failure ({err})"); return
            group_rows.append(dict(jpi=jpi, level=label, desc=full_desc, fractions=list(fr), signs=list(sg),
                                   gamma_total_ev=gtot, ex=variant["ex"] + de, frozen=float(of), gn=float("nan"),
                                   x=x, cols=(np.append(cols, ecol) if ecol is not None else cols),
                                   widths_ev=widths, shape=(desc, tuple(sg), tuple(phs)),
                                   channels=[dict(pair=k.pair, L=k.L, S=k.S) for k in amp_keys]))

        tg = time.time()
        for scale in s1:
            for fr, sg, desc, phs in shapes:
                evaluate(fr, sg, desc, phs, scale, 0.0)
        n1 = len(group_rows)
        if two_stage:
            # the best shapes of THIS group at the stage-1 widths, then every width and offset for them
            best_shapes, seen_sh = [], set()
            for r in sorted(group_rows, key=lambda r: r["frozen"]):
                if r["shape"] not in seen_sh:
                    seen_sh.add(r["shape"]); best_shapes.append(r["shape"])
                if len(best_shapes) >= refine_top:
                    break
            for desc, sg, phs in best_shapes:
                fr = next(f for f, g, d in pats if d == desc and tuple(g) == sg)
                for scale in scales:
                    for de in offsets:
                        evaluate(fr, sg, desc, phs, scale, de)
        bestg = min(group_rows, key=lambda r: r["frozen"]) if group_rows else None
        log(f"    {jpi}: {n} channels, {len(shapes)} shapes; {n1} stage-1 + {len(group_rows) - n1} refinement evaluations "
            f"in {time.time() - tg:.0f}s" + (f"; best frozen {bestg['desc']} signs {bestg['signs']} x{bestg['gamma_total_ev']/gamma_total_ev:g} "
                                            f"@{bestg['ex']:.3f} {bestg['frozen'] - obj0:+.1f}" if bestg else ""))
        rows.extend(group_rows)
    rows.sort(key=lambda r: r["frozen"])
    n_impr = sum(1 for r in rows if r["frozen"] < obj0)
    log(f"  {len(rows)} patterns frozen-evaluated in {time.time()-t0:.0f}s; {n_impr} below the incumbent; best frozen: "
        + ", ".join(f"{r['jpi']} {r['desc']} @{r['ex']:.3f} {r['frozen']:.1f}" for r in rows[:3]))
    t1 = time.time()
    # allocate the GN budget PER J^pi GROUP first.  Ranking globally let one group
    # take every slot on a large model, so five groups reported a gain of ~0 that had never
    # been GN-stepped at their own best pattern.
    picked, per = [], {}
    for r in rows:
        if len(picked) >= gn_top:
            break
        if per.get(r["jpi"], 0) < gn_per_jpi:
            per[r["jpi"]] = per.get(r["jpi"], 0) + 1
            picked.append(r)
    for r in rows:
        if len(picked) >= gn_top:
            break
        if r not in picked:
            picked.append(r)
    for r in picked:
        cols = np.concatenate([r["cols"], norm_cols]).astype(int)
        try:
            xg, og, _ = gn_step(ev, r["x"], cols, damping=1e-2, bounds=bounds)
            if og > r["frozen"]:
                xg, og, _ = gn_step(ev, r["x"], cols, damping=1.0, bounds=bounds)
            r["gn"] = float(min(og, r["frozen"])); r["x_gn"] = xg if og <= r["frozen"] else r["x"]
        except EvalError as err:
            log(f"    GN failure for {r['jpi']} {r['desc']}: {err}"); r["gn"] = r["frozen"]; r["x_gn"] = r["x"]
    rows.sort(key=lambda r: (r["gn"] if np.isfinite(r["gn"]) else r["frozen"]))
    log(f"  GN screen of the best {min(gn_top, len(picked))} in {time.time()-t1:.0f}s; best after GN: "
        + ", ".join(f"{r['jpi']} {r['desc']} signs {r['signs']} @{r['ex']:.3f} "
                    + (f"{r['gn']:.1f} (gain {obj0 - r['gn']:+.1f})" if np.isfinite(r["gn"])
                       else f"{r['frozen']:.1f} (frozen, not GN-stepped: {obj0 - r['frozen']:+.1f})") for r in rows[:5]))
    for r in rows:
        r["gain"] = obj0 - (r["gn"] if np.isfinite(r["gn"]) else r["frozen"])
        r["obj0"] = obj0
    if score_keys:
        # The total objective barely distinguishes J^pi: angle-integrated and total
        # cross sections only care that strength sits at the right energy, so a
        # wrong-J^pi level scores as well as the right one (a new
        # level won by a wide margin on a model far from its optimum).  Re-rank the GN'd
        # candidates on the DIFFERENTIAL segments near the candidate energy, which
        # is the shape information that actually carries J^pi.
        keys = [int(k) for k in score_keys]
        def local(x):
            sc = ev.score(x)
            return float(sum(sc.segment_chi2.get(k, 0.0) for k in keys))
        base_local = local(x_inc)
        scored = [r for r in rows if "x_gn" in r or np.isfinite(r["gn"])]
        for r in scored:
            try:
                r["local"] = local(r.get("x_gn", r["x"]))
                r["local_gain"] = base_local - r["local"]
            except EvalError as err:
                log(f"    local score failed for {r['jpi']} {r['desc']}: {err}")
        ranked = [r for r in scored if "local_gain" in r]
        ranked.sort(key=lambda r: -r["local_gain"])
        rest = [r for r in rows if r not in ranked]
        rows = ranked + rest
        log(f"  local re-rank on {len(keys)} differential segments near Ex {variant['ex']:.3f} "
            f"(reference chi2 {base_local:.1f}): "
            + ", ".join(f"{r['jpi']} {r['desc']} local {r['local_gain']:+.1f} (global {r['gain']:+.1f})"
                        for r in ranked[:5]))
    return rows


def best_per_jpi(rows: List[dict]) -> List[dict]:
    seen, out = set(), []
    for r in rows:
        if r["jpi"] not in seen:
            seen.add(r["jpi"]); out.append(r)
    return out


def to_level_candidate(row: dict, ex: Optional[float] = None, source: str = "levelscan") -> LevelCandidate:
    chans = [dict(pair=c["pair"], L=c["L"], S=c["S"], gamma=w) for c, w in zip(row["channels"], row["widths_ev"]) if w != 0]
    ex = float(row.get("ex", ex) if ex is None else ex)
    return LevelCandidate(energy=round(float(ex), 6), jpi=row["jpi"], channels=chans, source=source,
                          note=f"{row['desc']} signs {row['signs']} Gamma_total {row['gamma_total_ev']:.0f} eV; "
                               f"screen gain {row['gain']:+.1f}", energy_fixed=True)


def relocate(ev: Evaluator, x_inc, level: str, ex_new: float, policy: BoundsPolicy, half_window_kev: float = 5.0,
             max_nfev: int = 100, log=print):
    """Move an existing level to ``ex_new`` (MeV) and re-polish everything, the
    moved level's energy confined to +-``half_window_kev`` around ``ex_new`` so
    the polish cannot carry a narrow level away again.  Returns the FitResult."""
    from .moves import energy_kick, apply
    from .local import polish
    x = apply(x_inc, ev.keys, energy_kick(ev.keys, level, float(ex_new)))
    b = default_bounds(ev, x, policy)
    for j, k in enumerate(ev.keys):
        if k.kind == "energy" and k.level == level:
            b.lb[j] = ex_new - half_window_kev * 1e-3
            b.ub[j] = ex_new + half_window_kev * 1e-3
            b.fixed[j] = False
    log(f"  relocate {level} -> {ex_new:.4f} MeV (window +-{half_window_kev:g} keV), polish {max_nfev} evaluations")
    return polish(ev, x, bounds=b, policy=policy, max_nfev=max_nfev, log=log)


def review_figure(ev: Evaluator, x_before, x_after, ex: float, path: str, half_width: float = 0.06,
                  title: str = "", label_before: str = "before", label_after: str = "after") -> dict:
    """The human check of a level change: for every active segment with data
    inside ``ex`` +- ``half_width`` MeV, the data and the two calculations,
    excitation functions versus Ex for angle-integrated / fixed-angle segments
    and the angular distributions at the energy nearest ``ex`` for segments
    that hold several angles at one energy.  ``ev`` must be a single session
    (``data()`` and ``cross()``).  Writes a PNG (and PDF) to ``path``."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = ev.data()
    c0 = ev.cross(np.asarray(x_before, float)); c1 = ev.cross(np.asarray(x_after, float))
    panels = []
    for i, (key, (a, b)) in enumerate(ev.row_ranges):
        if i >= len(d["excitation"]):
            continue
        Ex = np.asarray(d["excitation"][i], float); y = np.asarray(d["cross"][i], float); e = np.abs(np.asarray(d["cross_err"][i], float))
        sel = np.abs(Ex - ex) <= half_width
        if sel.sum() < 3:
            continue
        seg = ev.segments.by_key(key)
        name = os.path.basename(seg.data_file) if seg is not None else str(key)
        xaxis, xlabel = Ex[sel], "Ex (MeV)"
        if Ex[sel].max() - Ex[sel].min() < 2e-3 and seg is not None:
            # one energy, many angles (an angular distribution / analyzing power):
            # take the angles from the data file rows of this segment, in order
            try:
                from .model import read_points
                base_dir = os.path.dirname(os.path.abspath(ev.azr)) if hasattr(ev, "azr") else os.getcwd()
                pts = read_points(os.path.join(base_dir, seg.data_file))
                lab = np.asarray(d["energies"][i], float)
                rows = [r for r in pts if seg.energy_min - 1e-9 <= r[0] <= seg.energy_max + 1e-9]
                if len(rows) == lab.size:
                    xaxis, xlabel = np.array([r[1] for r in rows])[sel], "angle (deg)"
            except Exception:
                pass
        panels.append(dict(key=key, name=name, x=xaxis, xlabel=xlabel, y=y[sel], e=e[sel], f0=np.asarray(c0[i], float)[sel],
                           f1=np.asarray(c1[i], float)[sel], angle=(seg.angle_min if seg is not None else None),
                           ex_here=float(np.mean(Ex[sel])) if xlabel != "Ex (MeV)" else None))
    if not panels:
        return dict(path=None, panels=0)
    n = len(panels); cols = 2 if n > 1 else 1; rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.2 * cols, 3.4 * rows), squeeze=False)
    for ax, pnl in zip(axes.ravel(), panels):
        o = np.argsort(pnl["x"])
        ax.errorbar(pnl["x"][o], pnl["y"][o], pnl["e"][o], fmt="o", ms=2.5, color="k", lw=0.7, label="data")
        ax.plot(pnl["x"][o], pnl["f0"][o], "-", color="tab:red", lw=1.2, label=label_before)
        ax.plot(pnl["x"][o], pnl["f1"][o], "-", color="tab:blue", lw=1.2, label=label_after)
        if pnl["xlabel"] == "Ex (MeV)":
            ax.axvline(ex, color="gray", ls=":", lw=0.8)
            ang = f" {pnl['angle']:.0f} deg" if pnl["angle"] not in (None, 0.0) else ""
            ax.set_title(f"seg {pnl['key']} {pnl['name']}{ang}", fontsize=8)
        else:
            ax.set_title(f"seg {pnl['key']} {pnl['name']} at Ex {pnl['ex_here']:.3f} MeV", fontsize=8)
        ax.set_xlabel(pnl["xlabel"], fontsize=8)
        ax.tick_params(labelsize=7); ax.legend(fontsize=6.5)
    for ax in axes.ravel()[n:]:
        ax.set_visible(False)
    if title:
        fig.suptitle(title, fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=130)
    if path.endswith(".png"):
        fig.savefig(path[:-4] + ".pdf")
    plt.close(fig)
    return dict(path=path, panels=n)
