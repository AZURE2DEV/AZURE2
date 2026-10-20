#!/usr/bin/env python3
"""Broad-structure locator: where does the BLURRED data differ from the equally blurred fit?  (resonance-search skill, section 1b;
DeBoer 2026-10-11: "blur the image of the cross section so only the broad levels are left visible")
Input: the per-point table of a fit (dump_residuals.py).  Per detector curve (one lab angle of one segment, >= MIN_PTS energies) and
blur width D (Gaussian kernel regression, sigma = D, points within 1.5 D of the curve ends dropped):
    data  n*d = fit * d / (d + r*e)    (r = (fit - n d)/(n e) from the table, so the data carry their fitted normalization)
    delta(E) = (blur(n*d) - blur(fit)) / |blur(fit)|
A DEVIATION is a maximal stretch of one sign with |delta| >= REL_MIN and extent >= D; its score is
    (extent / D) * mean( delta^2 / (FLOOR^2 + sigma_b^2) )      sigma_b = statistical error of the blurred data, relative
i.e. roughly chi^2 with a FLOOR (5 %) systematic per independent blur width -- so a dense data set cannot win on point count and a
fitted normalization offset cannot win at all (it is absorbed by the norm).  Deviations of all curves are spread over their extent on
a 2 keV Ex map, smoothed, cut at the map's minima; regions ranked by integral, with sign (data above / below the fit), data sets and
the model's levels inside (optional parameters.out + .azr).
READ WITH the validity table of blur_compare.py: where narrow levels carry the strength (13C+alpha ND 2021 (a,n0) below 10 MeV), a
blurred deviation can be MISSING NARROW STRENGTH (area / D), not a missing broad level.
usage: blur_map.py table.npz out.txt [--scales 0.1,0.2] [--levels parameters.out --azr fit.azr]"""
import sys, re, argparse
import numpy as np
from collections import defaultdict
from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelmin

MIN_PTS, REL_MIN, FLOOR, DE = 12, 0.10, 0.05, 0.002


def detector_labels(ang, gap=2.0):
    """One label per lab detector (the c.m. angle drifts with energy; cut sorted angles at gaps > gap deg)."""
    a = np.nan_to_num(np.asarray(ang, float), nan=-1.0)
    order = np.argsort(a); s = a[order]
    cid = np.concatenate([[0], np.cumsum(np.diff(s) > gap)])
    lab = np.empty_like(a)
    for c in np.unique(cid):
        lab[order[cid == c]] = round(float(s[cid == c].mean()), 1)
    return lab


def blur(e, y, s, w=None):
    """Gaussian kernel regression at the points e (sorted), window +-4 sigma; optional per-point weights w (returns value, sum w^2 / (sum w)^2)."""
    lo = np.searchsorted(e, e - 4 * s); hi = np.searchsorted(e, e + 4 * s, side="right")
    out = np.empty_like(y); var = np.empty_like(y)
    for i in range(len(e)):
        k = np.exp(-0.5 * ((e[lo[i]:hi[i]] - e[i]) / s) ** 2)
        out[i] = k @ y[lo[i]:hi[i]] / k.sum()
        if w is not None:
            var[i] = (k ** 2) @ (w[lo[i]:hi[i]] ** 2) / k.sum() ** 2
    return out, var


def curves(z):
    out = []
    for k in np.unique(z["key"]):
        idx = np.where(z["key"] == k)[0]; a = detector_labels(z["ang"][idx]); e = np.round(z["ex"][idx], 6)
        if len(np.unique(e)) < len(np.unique(a)):
            continue                      # single-energy angular distributions: no excitation function to blur
        for av in np.unique(a):
            s = idx[a == av]; s = s[np.argsort(z["ex"][s])]
            keep = np.concatenate([[True], np.diff(z["ex"][s]) > 1e-7]); s = s[keep]
            if len(s) >= MIN_PTS:
                out.append((int(k), float(av), s))
    return out


def deviations(z, D):
    devs = []
    for k, av, s in curves(z):
        e = z["ex"][s]
        if e[-1] - e[0] < 4 * D:
            continue
        d, er, r, f = z["dat"][s], z["err"][s], z["res"][s], z["fit"][s]
        den = d + r * er
        ok = np.abs(den) > 0
        nd = np.where(ok, f * d / np.where(ok, den, 1), d); ne = np.where(ok, np.abs(f * er / np.where(ok, den, 1)), er)
        bd, vd = blur(e, nd, D, ne); bf, _ = blur(e, f, D)
        m = (e >= e[0] + 1.5 * D) & (e <= e[-1] - 1.5 * D)
        if m.sum() < 3:
            continue
        E, dl = e[m], (bd[m] - bf[m]) / np.maximum(np.abs(bf[m]), 1e-300)
        sb2 = vd[m] / np.maximum(bf[m] ** 2, 1e-300)
        sg = np.sign(dl) * (np.abs(dl) >= REL_MIN)
        i = 0
        while i < len(E):
            if sg[i] == 0:
                i += 1; continue
            j = i
            while j + 1 < len(E) and sg[j + 1] == sg[i]:
                j += 1
            ext = E[j] - E[i]
            if ext >= D:
                sc = (ext / D) * float(np.mean(dl[i:j + 1] ** 2 / (FLOOR ** 2 + sb2[i:j + 1])))
                devs.append(dict(key=k, ang=av, name=str(z["name"][s[0]]), lo=float(E[i]), hi=float(E[j]), sign=int(sg[i]),
                                 peak=float(E[i + np.argmax(np.abs(dl[i:j + 1]))]), dmax=float(dl[i + np.argmax(np.abs(dl[i:j + 1]))]),
                                 score=sc))
            i = j + 1
    return devs


def read_levels(par, azr):
    lev, cur = [], None
    for line in open(par):
        m = re.match(r"J =\s*(\S+)\s+E_level =\s*(\S+)", line)
        if m:
            J = m.group(1); jpi = J.replace("0.5", "1/2").replace("1.5", "3/2").replace("2.5", "5/2").replace("3.5", "7/2") \
                .replace("4.5", "9/2").replace("5.5", "11/2").replace("6.5", "13/2")
            cur = dict(jpi=jpi, e=float(m.group(2)), g=0.0, fixed=False); lev.append(cur); continue
        m = re.search(r"G\s+=\s+(\S+) keV", line)
        if m and cur is not None:
            cur["g"] += abs(float(m.group(1)))
    if azr:
        sys.path.insert(0, "/users/rdeboer1/AZURE2")
        from pyazr import AzrModel
        fx = [(l.jpi, l.energy, l.fixed) for l in AzrModel.from_file(azr).levels if l.active]
        for L in lev:
            c = [f for f in fx if f[0] == L["jpi"] and abs(f[1] - L["e"]) < 0.003]
            L["fixed"] = bool(c and c[0][2])
    return lev


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("table"); ap.add_argument("out")
    ap.add_argument("--scales", default="0.1,0.2"); ap.add_argument("--levels"); ap.add_argument("--azr")
    a = ap.parse_args()
    z = np.load(a.table); out = open(a.out, "w")
    lev = read_levels(a.levels, a.azr) if a.levels else []
    for D in [float(x) for x in a.scales.split(",")]:
        devs = deviations(z, D)
        grid = np.arange(z["ex"].min() - 0.05, z["ex"].max() + 0.05, DE)
        M, Mlow = np.zeros_like(grid), np.zeros_like(grid)
        for d in devs:
            m = (grid >= d["lo"]) & (grid <= d["hi"])
            if not m.any():
                continue
            M[m] += d["score"] / m.sum()
            if d["sign"] < 0:
                Mlow[m] += d["score"] / m.sum()
        S = gaussian_filter1d(M, max(1.0, 0.25 * D / DE))
        cuts = [0] + list(argrelmin(S, order=3)[0]) + [len(grid) - 1]
        regions = []
        for i0, i1 in zip(cuts[:-1], cuts[1:]):
            tot = M[i0:i1 + 1].sum()
            if tot < 1.0:
                continue
            lo_e, hi_e = grid[i0], grid[i1]; pk = grid[i0 + np.argmax(S[i0:i1 + 1])]
            by = defaultdict(float)
            for d in devs:
                ov = max(0.0, min(d["hi"], hi_e) - max(d["lo"], lo_e)) / max(d["hi"] - d["lo"], 1e-9)
                if ov > 0:
                    by[f"{d['name']} {d['ang']:g}°"] += d["score"] * ov
            regions.append(dict(lo=lo_e, hi=hi_e, peak=pk, tot=tot, below=Mlow[i0:i1 + 1].sum() / tot,
                                cur=sorted(by.items(), key=lambda t: -t[1])))
        regions.sort(key=lambda g: -g["tot"])
        print(f"=== blur sigma D = {1000*D:.0f} keV: {len(devs)} deviations (|delta| >= {REL_MIN:.0%} over >= D, "
              f"{FLOOR:.0%} floor) in {len({(d['key'], d['ang']) for d in devs})} curves; {len(regions)} regions\n", file=out)
        print(f"{'rank':>4} {'Ex range (MeV)':>15} {'peak':>7} {'score':>7} {'data<fit':>8}  curves (score)", file=out)
        for i, g in enumerate(regions[:30], 1):
            s = "; ".join(f"{n} {v:.0f}" for n, v in g["cur"][:4])
            print(f"{i:4d} {g['lo']:7.3f}-{g['hi']:<7.3f} {g['peak']:7.3f} {g['tot']:7.0f} {100*g['below']:7.0f}%  {s}", file=out)
            if lev:
                inl = [L for L in lev if g["lo"] - 0.01 <= L["e"] <= g["hi"] + 0.01 and L["g"] >= 10]
                print("      levels (Gamma >= 10 keV): " + ("; ".join(f"{L['jpi']} {L['e']:.3f} {L['g']:.3g} keV" + (" FIXED" if L["fixed"] else "")
                                                             for L in inl) or "none"), file=out)
        print("", file=out)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
