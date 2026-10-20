#!/usr/bin/env python3
"""Ex density map of coherent residual runs -> ranked regions.  Run finding as rank_regions.py (each angle / energy of each segment
one curve; maximal same-sign runs, n >= 3, coherent chi2 (sum r)^2/n >= 25).  Each run's coherent chi2 is spread uniformly over its Ex
extent (min 10 keV) on a 2 keV grid; runs wider than MAX_SPAN are trends (listed, not mapped).  The map is smoothed (sigma 6 keV) and
cut into regions at its local minima; a region = one basin, ranked by integral.  Residual convention r = (fit - data)/sigma:
"model low" = runs with r < 0.  usage: rank_map.py resid.npz out.txt [parameters.out fit.azr]   (the last two list the model's levels in each region)"""
import sys
import numpy as np
from collections import defaultdict
from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelmin

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

MIN_RUN, MIN_SCORE, MAX_SPAN, MINW, DE, SIG = 3, 25.0, 0.6, 0.010, 0.002, 0.006
z = np.load(sys.argv[1]); out = open(sys.argv[2], "w") if len(sys.argv) > 2 else sys.stdout
key, name, ex, ang, res = z["key"], z["name"], z["ex"], z["ang"], z["res"]
curves = {}
for k in np.unique(key):
    idx = np.where(key == k)[0]; a = detector_labels(ang[idx]); e = np.round(ex[idx], 5)
    if len(np.unique(e)) >= len(np.unique(a)):
        for av in np.unique(a): s = idx[a == av]; curves[(k, "ang", av)] = s[np.argsort(ex[s])]
    else:
        for ev in np.unique(e): s = idx[e == ev]; curves[(k, "E", ev)] = s[np.argsort(ang[s])]
runs = []
for (k, kind, v), idx in curves.items():
    r = res[idx]; sg = np.sign(r); i = 0
    while i < len(r):
        j = i
        while j + 1 < len(r) and sg[j + 1] == sg[i]: j += 1
        n = j - i + 1
        if n >= MIN_RUN and r[i:j + 1].sum() ** 2 / n >= MIN_SCORE:
            sel = idx[i:j + 1]
            runs.append(dict(name=str(name[idx[0]]), key=int(k), kind=kind, v=float(v), lo=float(ex[sel].min()), hi=float(ex[sel].max()),
                             n=n, sign=int(sg[i]), score=float(r[i:j + 1].sum() ** 2 / n)))
        i = j + 1
trends = [r for r in runs if r["hi"] - r["lo"] > MAX_SPAN]; loc = [r for r in runs if r["hi"] - r["lo"] <= MAX_SPAN]
grid = np.arange(ex.min() - 0.05, ex.max() + 0.05, DE); M = np.zeros_like(grid); Mlow = np.zeros_like(grid)
for r in loc:
    lo, hi = r["lo"], r["hi"]
    if hi - lo < MINW: c = 0.5 * (lo + hi); lo, hi = c - MINW / 2, c + MINW / 2
    m = (grid >= lo) & (grid <= hi)
    if not m.any(): m[np.argmin(abs(grid - 0.5 * (lo + hi)))] = True
    M[m] += r["score"] / m.sum()
    if r["sign"] < 0: Mlow[m] += r["score"] / m.sum()
S = gaussian_filter1d(M, SIG / DE)
cuts = [0] + list(argrelmin(S, order=3)[0]) + [len(grid) - 1]
regions = []
for a, b in zip(cuts[:-1], cuts[1:]):
    tot = M[a:b + 1].sum()
    if tot < 100: continue
    pk = a + np.argmax(S[a:b + 1]); half = S[pk] / 2
    l = pk
    while l > a and S[l] > half: l -= 1
    h = pk
    while h < b and S[h] > half: h += 1
    lo_e, hi_e = grid[a], grid[b]
    by = defaultdict(float)
    for r in loc:
        lo_r, hi_r = r["lo"], r["hi"]
        if hi_r - lo_r < MINW: c = 0.5 * (lo_r + hi_r); lo_r, hi_r = c - MINW / 2, c + MINW / 2
        ov = max(0.0, min(hi_r, hi_e) - max(lo_r, lo_e)) / (hi_r - lo_r)
        if ov > 0: by[r["name"]] += r["score"] * ov
    regions.append(dict(lo=lo_e, hi=hi_e, peak=grid[pk], fwhm=grid[h] - grid[l], tot=tot, low=Mlow[a:b + 1].sum() / tot,
                        sets=sorted(by.items(), key=lambda t: -t[1])))
regions.sort(key=lambda g: -g["tot"])
levels = []
if len(sys.argv) > 4:
    import re
    sys.path.insert(0, "/users/rdeboer1/AZURE2"); from pyazr import AzrModel
    fixed = [(l.jpi, l.energy, l.fixed) for l in AzrModel.from_file(sys.argv[4]).levels if l.active]
    cur = None
    for line in open(sys.argv[3]):
        m = re.match(r"J =\s*(\S+)\s+E_level =\s*(\S+)", line)
        if m:
            J = m.group(1); jpi = J.replace("0.5", "1/2").replace("1.5", "3/2").replace("2.5", "5/2").replace("3.5", "7/2").replace("4.5", "9/2").replace("5.5", "11/2").replace("6.5", "13/2")
            cur = dict(jpi=jpi, e=float(m.group(2)), g=0.0); levels.append(cur); continue
        m = re.search(r"G\s+=\s+(\S+) keV", line)
        if m and cur is not None: cur["g"] += float(m.group(1))
    for L in levels:
        c = [f for f in fixed if f[0] == L["jpi"] and abs(f[1] - L["e"]) < 0.003]
        L["fixed"] = bool(c and c[0][2])
np.savez(sys.argv[1].replace(".npz", "_map.npz"), grid=grid, M=M, Mlow=Mlow, S=S)
print(f"{len(runs)} runs (n >= {MIN_RUN}, coherent chi2 >= {MIN_SCORE:g}); coherent chi2 mapped {M.sum():.0f}; "
      f"{len(regions)} regions >= 100; {len(trends)} trends wider than {MAX_SPAN} MeV\n", file=out)
print(f"{'rank':>4} {'Ex range (MeV)':>15} {'peak':>7} {'FWHM':>6} {'coh chi2':>9} {'model low':>9}  data sets (coherent chi2)", file=out)
for i, g in enumerate(regions[:45], 1):
    s = "; ".join(f"{n} {v:.0f}" for n, v in g["sets"][:4])
    print(f"{i:4d} {g['lo']:7.3f}-{g['hi']:<7.3f} {g['peak']:7.3f} {1000*g['fwhm']:5.0f}k {g['tot']:9.0f} {100*g['low']:8.0f}%  {s}", file=out)
    inl = [L for L in levels if g["lo"] - 0.01 <= L["e"] <= g["hi"] + 0.01]
    if levels: print("      levels: " + ("; ".join(f"{L['jpi']} {L['e']:.4f} {L['g']:.3g} keV{' FIXED' if L['fixed'] else ''}" for L in inl) or "none"), file=out)
print(f"\ntrends (single runs wider than {MAX_SPAN} MeV):", file=out)
for r in sorted(trends, key=lambda r: -r["score"])[:12]:
    print(f"  {r['lo']:7.3f}-{r['hi']:<7.3f} {r['name']} key {r['key']} {r['kind']}={r['v']:g} n {r['n']} {'model low' if r['sign']<0 else 'model high'} coh {r['score']:.0f}", file=out)
