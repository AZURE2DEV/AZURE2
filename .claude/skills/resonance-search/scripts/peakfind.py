#!/usr/bin/env python3
"""Data-only resonance candidates: peaks, dips and dispersion shapes in every excitation function, clustered in Ex.

For a reaction with NO R-matrix calculation yet (resonance-search skill, section 1a).  Input = the per-point table written by
dump_residuals.py (only key, name, ex, ang, dat, err are used; the fit columns are ignored, so any .azr with the pairs and the data
segments will do -- its levels do not matter).  For each excitation function (one angle of one segment, >= MIN_PTS points):
  1. baseline b(E) = running median over a window of W MeV, for each W in SCALES (a structure is seen at the scales wider than it);
  2. significance z = (y - b) / sigma and relative size rel = (y - b) / |b|;
  3. scipy.signal.find_peaks on +z (peaks) and -z (dips): prominence >= Z_MIN and |rel| >= REL_MIN at the extremum (underestimated
     error bars alone cannot make a candidate); FWHM from peak_widths at half prominence, in MeV; detections of the same extremum at
     several scales are merged (the narrowest scale that sees it gives the width);
  4. an adjacent peak and dip of the same curve closer than the larger of their widths form one DISPERSION candidate at the midpoint.
Detections from all curves go into an Ex density map (each a unit-area Gaussian of its own FWHM, min 3 keV, weight min(z, ZCAP)); the
map's peaks (prominence >= MAP_PROM x its maximum) are the candidates, each collecting the detections inside its half-maximum width.
(Interval clustering chains dense multi-angle data into one multi-MeV cluster -- first try, 2026-10-11.)  Each candidate lists the number of
distinct data sets, channels (exit pairs from the segment name) and curves that see it, the median width, the summed capped
significance (min(z, ZCAP) per curve, so one 20-angle set cannot outvote everything) and the peak / dip / dispersion mix.
Ranked by distinct data sets, then by score.  Optional: --levels parameters.out to score the list blind against a known level set.
usage: peakfind.py table.npz out.txt [--levels parameters.out] [--ex-range lo,hi]"""
import sys, re, argparse
import numpy as np
from collections import defaultdict
from scipy.signal import find_peaks, peak_widths
from scipy.ndimage import median_filter


def detector_labels(ang, gap=2.0):
    """One label per detector: the c.m. angle of a fixed lab detector drifts with energy (several degrees over a long run), so
    rounding c.m. angles splits one detector into many fragments (13C+alpha ND 2021: 558 'curves' for 146 detectors).  Sorted
    angles are cut where consecutive values differ by > gap degrees; each label is its cluster's mean angle (0.1 deg)."""
    a = np.nan_to_num(np.asarray(ang, float), nan=-1.0)
    order = np.argsort(a); s = a[order]
    cid = np.concatenate([[0], np.cumsum(np.diff(s) > gap)])
    lab = np.empty_like(a)
    for c in np.unique(cid):
        lab[order[cid == c]] = round(float(s[cid == c].mean()), 1)
    return lab

SCALES = (0.05, 0.15, 0.5)      # MeV baseline windows
MIN_PTS, Z_MIN, REL_MIN, ZCAP, MAP_PROM = 8, 4.0, 0.05, 10.0, 0.002


def running_median(e, y, w):
    """Median of y over |E - E_i| <= w/2 (irregular grids: by index window sized from the local step)."""
    out = np.empty_like(y)
    lo = np.searchsorted(e, e - w / 2); hi = np.searchsorted(e, e + w / 2, side="right")
    for i in range(len(e)):
        out[i] = np.median(y[lo[i]:hi[i]])
    return out


def detect(e, y, s):
    """Peaks (+1) and dips (-1) of one excitation function, merged over the baseline scales."""
    found = []
    step = np.median(np.diff(e)) if len(e) > 1 else 0.0
    for w in SCALES:
        if w < 4 * step or w > 0.7 * (e[-1] - e[0]):
            continue
        b = running_median(e, y, w)
        z = (y - b) / np.where(s > 0, s, np.inf)
        rel = (y - b) / np.maximum(np.abs(b), 1e-300)
        for sign in (1, -1):
            idx, props = find_peaks(sign * z, prominence=Z_MIN)
            if not len(idx):
                continue
            fw = peak_widths(sign * z, idx, rel_height=0.5)
            for k, i in enumerate(idx):
                if abs(rel[i]) < REL_MIN:
                    continue
                el = np.interp(fw[2][k], np.arange(len(e)), e); er = np.interp(fw[3][k], np.arange(len(e)), e)
                found.append(dict(ex=float(e[i]), sign=sign, z=float(props["prominences"][k]), rel=float(rel[i]),
                                  width=float(max(er - el, step)), step=float(step), scale=w))
    merged = []                       # same extremum at several scales -> keep the narrowest-scale width, the largest z
    for d in sorted(found, key=lambda d: d["scale"]):
        m = [q for q in merged if q["sign"] == d["sign"] and abs(q["ex"] - d["ex"]) <= max(2 * d["step"], 0.002)]
        if m:
            m[0]["z"] = max(m[0]["z"], d["z"])
        else:
            merged.append(dict(d))
    merged.sort(key=lambda d: d["ex"])
    out, used = [], set()             # dispersion: adjacent peak + dip
    for i in range(len(merged) - 1):
        a, c = merged[i], merged[i + 1]
        if i in used or a["sign"] == c["sign"]:
            continue
        if c["ex"] - a["ex"] <= max(a["width"], c["width"]):
            out.append(dict(ex=0.5 * (a["ex"] + c["ex"]), sign=0, z=max(a["z"], c["z"]), rel=max(abs(a["rel"]), abs(c["rel"])),
                            width=c["ex"] - a["ex"] + 0.5 * (a["width"] + c["width"]), step=a["step"], scale=a["scale"]))
            used |= {i, i + 1}
    out += [d for i, d in enumerate(merged) if i not in used]
    return out


def channel_of(name):
    m = re.search(r"(nn\d|na\d|an\d|aa\d|n\dg|a\dg|ntot|N_TOT|nng|AN_TOT|aa_|_an_|13Can|13Caa)", name)
    return m.group(1) if m else name


def read_levels(path):
    lev, cur = [], None
    for line in open(path):
        m = re.match(r"J =\s*(\S+)\s+E_level =\s*(\S+)", line)
        if m:
            cur = dict(jpi=m.group(1), e=float(m.group(2)), g=0.0); lev.append(cur); continue
        m = re.search(r"G\s+=\s+(\S+) keV", line)
        if m and cur is not None:
            cur["g"] += abs(float(m.group(1)))
    return lev


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("table"); ap.add_argument("out")
    ap.add_argument("--levels"); ap.add_argument("--ex-range")
    a = ap.parse_args()
    z = np.load(a.table)
    key, name, ex, ang, dat, err = z["key"], z["name"], z["ex"], z["ang"], z["dat"], z["err"]
    dets = []
    for k in np.unique(key):
        idx = np.where(key == k)[0]; an = detector_labels(ang[idx])
        for av in np.unique(an):
            s = idx[an == av]; s = s[np.argsort(ex[s])]
            e, y, sg = ex[s], dat[s], err[s]
            keep = np.concatenate([[True], np.diff(e) > 1e-7])      # duplicate energies: keep the first
            e, y, sg = e[keep], y[keep], sg[keep]
            if len(e) < MIN_PTS:
                continue
            for d in detect(e, y, sg):
                d.update(key=int(k), name=str(name[idx[0]]), ang=float(av), chan=channel_of(str(name[idx[0]])))
                dets.append(d)
    if a.ex_range:
        lo, hi = map(float, a.ex_range.split(",")); dets = [d for d in dets if lo <= d["ex"] <= hi]
    dets.sort(key=lambda d: d["ex"])
    # density map: each detection = a Gaussian of its own FWHM (min 3 keV), weight min(z, ZCAP), normalised to unit area
    grid = np.arange(min(d["ex"] for d in dets) - 0.2, max(d["ex"] for d in dets) + 0.2, 0.001)
    Mp = np.zeros_like(grid)
    for d in dets:
        sig = max(d["width"], 0.003) / 2.355
        g = np.exp(-0.5 * ((grid - d["ex"]) / sig) ** 2); Mp += min(d["z"], ZCAP) * g / g.sum()
    pk, pr = find_peaks(Mp, prominence=MAP_PROM * Mp.max())
    hw = peak_widths(Mp, pk, rel_height=0.5)
    clusters = []
    for j, p in enumerate(pk):
        lo_e, hi_e = np.interp(hw[2][j], np.arange(len(grid)), grid), np.interp(hw[3][j], np.arange(len(grid)), grid)
        D = [d for d in dets if lo_e - 0.002 <= d["ex"] <= hi_e + 0.002]
        if not D:
            continue
        clusters.append(dict(d=D, ex_map=float(grid[p]), lo=float(lo_e), hi=float(hi_e)))
    for c in clusters:
        D = c["d"]; w = np.array([min(q["z"], ZCAP) for q in D])
        c["ex"] = c["ex_map"]
        c["width"] = float(np.median([q["width"] for q in D])); c["step"] = float(np.median([q["step"] for q in D]))
        c["sets"] = sorted({q["name"] for q in D}); c["chans"] = sorted({q["chan"] for q in D})
        c["curves"] = len({(q["key"], q["ang"]) for q in D}); c["score"] = float(w.sum())
        c["mix"] = "/".join(f"{sum(1 for q in D if q['sign'] == s)}{t}" for s, t in ((1, "p"), (-1, "d"), (0, "x")))
    clusters.sort(key=lambda c: (-len(c["sets"]), -c["score"]))
    out = open(a.out, "w")
    print(f"{len(dets)} detections in {len({(d['key'], d['ang']) for d in dets})} curves -> {len(clusters)} clusters "
          f"(prominence >= {Z_MIN} sigma, |rel| >= {REL_MIN:.0%}, scales {SCALES} MeV)\n", file=out)
    print(f"{'rank':>4} {'Ex':>8} {'width':>8} {'sets':>4} {'curves':>6} {'score':>6} {'p/d/x':>9}  channels | data sets", file=out)
    for i, c in enumerate(clusters, 1):
        res = " (res.)" if c["width"] <= 2.5 * c["step"] else ""
        print(f"{i:4d} {c['ex']:8.4f} {1000*c['width']:6.1f}k{res:6s} {len(c['sets']):4d} {c['curves']:6d} {c['score']:6.0f} {c['mix']:>9}  "
              f"{','.join(c['chans'])} | {'; '.join(c['sets'][:5])}", file=out)
    if a.levels:
        lev = [L for L in read_levels(a.levels) if ex.min() <= L["e"] <= ex.max() and L["g"] >= 1.0]
        if a.ex_range:
            lev = [L for L in lev if lo <= L["e"] <= hi]
        print("\nBLIND CHECK against the level set (levels with Gamma >= 1 keV inside the data's Ex range):", file=out)
        bins = [(1, 10), (10, 100), (100, 500), (500, 1e9)]
        for b0, b1 in bins:
            Ls = [L for L in lev if b0 <= L["g"] < b1]
            hit = [L for L in Ls if any(abs(c["ex"] - L["e"]) <= max(0.010, 0.5e-3 * L["g"]) for c in clusters)]
            print(f"  Gamma {b0:g}-{b1:g} keV: {len(hit)}/{len(Ls)} found within max(10 keV, Gamma/2)", file=out)
        multi = [c for c in clusters if len(c["sets"]) >= 2]
        near = [c for c in multi if any(abs(c["ex"] - L["e"]) <= max(0.010, 0.5e-3 * L["g"]) for L in lev)]
        print(f"  clusters seen by >= 2 data sets: {len(multi)}, of which {len(near)} lie at a level; the other "
              f"{len(multi) - len(near)} are false alarms OR structure the level set lacks", file=out)
        for L in sorted(lev, key=lambda L: -L["g"]):
            if 10 <= L["g"] < 500 and not any(abs(c["ex"] - L["e"]) <= max(0.010, 0.5e-3 * L["g"]) for c in clusters):
                print(f"    missed: {L['jpi']} {L['e']:.4f} MeV, Gamma {L['g']:.1f} keV", file=out)
    print(f"wrote {a.out}: {len(clusters)} clusters")


if __name__ == "__main__":
    main()
