#!/usr/bin/env python3
"""Does blurring the cross section remove the narrow levels and leave the broad ones?  (DeBoer 2026-10-11: "blur the image of the
cross section so only the broad levels are left visible")
Inputs: the full fit's and the broad-only model's per-point tables (dump_residuals.py; same rows).  Per excitation function (one angle
of one segment, >= 12 energies), Gaussian kernel regression with sigma = D at the data energies, points within 1.5 D of either end
dropped.  Per data set and D, point-weighted rms of relative differences (normalised by |blur(full)|):
  A = full - broad              (unblurred: how much the narrow levels matter)
  B = blur(full) - blur(broad)  (what of the narrow levels survives the blur)
  C = blur(full) - broad        (is the blurred full model the broad component itself?)
usage: blur_compare.py full.npz broad.npz out.txt"""
import sys
import numpy as np
from collections import defaultdict

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

SC = (0.02, 0.05, 0.1, 0.2, 0.4)


def blur(e, y, s):
    """Gaussian kernel regression at the points e (sorted), window +-4 sigma."""
    lo = np.searchsorted(e, e - 4 * s); hi = np.searchsorted(e, e + 4 * s, side="right"); out = np.empty_like(y)
    for i in range(len(e)):
        w = np.exp(-0.5 * ((e[lo[i]:hi[i]] - e[i]) / s) ** 2); out[i] = w @ y[lo[i]:hi[i]] / w.sum()
    return out


def curves(z):
    out = []
    for k in np.unique(z["key"]):
        idx = np.where(z["key"] == k)[0]; a = detector_labels(z["ang"][idx])
        for av in np.unique(a):
            s = idx[a == av]; s = s[np.argsort(z["ex"][s])]
            if len(np.unique(np.round(z["ex"][s], 6))) >= 12: out.append(s)
    return out


def main():
    F, Bz, out = np.load(sys.argv[1]), np.load(sys.argv[2]), open(sys.argv[3], "w")
    assert np.allclose(F["ex"], Bz["ex"]) and (F["key"] == Bz["key"]).all()
    full, broad, ex, name = F["fit"], Bz["fit"], F["ex"], F["name"]
    acc = defaultdict(lambda: defaultdict(lambda: np.zeros(4)))   # [set][scale] -> sumA2, sumB2, sumC2, n
    for s in curves(F):
        e = ex[s]; yf, yb = full[s], broad[s]
        for D in SC:
            if e[-1] - e[0] < 6 * D: continue
            bf, bb = blur(e, yf, D), blur(e, yb, D)
            m = (e >= e[0] + 1.5 * D) & (e <= e[-1] - 1.5 * D)
            if m.sum() < 3: continue
            nrm = np.maximum(np.abs(bf[m]), 1e-30)
            v = acc[str(name[s[0]])][D]
            v += [np.sum(((yf - yb)[m] / nrm) ** 2), np.sum(((bf - bb)[m] / nrm) ** 2), np.sum(((bf - yb)[m] / nrm) ** 2), m.sum()]
    print("rms relative difference (%), per data set and blur sigma D:  A = full-broad unblurred | B = blur(full)-blur(broad) | "
          "C = blur(full)-broad\n", file=out)
    print(f"{'data set':40s} " + " ".join(f"{'D='+str(int(1000*D))+'k A/B/C':>20s}" for D in SC), file=out)
    tot = defaultdict(lambda: np.zeros(4))
    for nm in sorted(acc, key=lambda n: -sum(v[3] for v in acc[n].values())):
        row = []
        for D in SC:
            v = acc[nm].get(D)
            if v is None or v[3] == 0: row.append(f"{'-':>20s}"); continue
            tot[D] += v
            row.append(f"{100*np.sqrt(v[0]/v[3]):6.1f}/{100*np.sqrt(v[1]/v[3]):5.1f}/{100*np.sqrt(v[2]/v[3]):5.1f}".rjust(20))
        print(f"{nm[:40]:40s} " + " ".join(row), file=out)
    print(f"{'ALL (point-weighted)':40s} " + " ".join(f"{100*np.sqrt(tot[D][0]/tot[D][3]):6.1f}/{100*np.sqrt(tot[D][1]/tot[D][3]):5.1f}/"
                                                  f"{100*np.sqrt(tot[D][2]/tot[D][3]):5.1f}".rjust(20) for D in SC), file=out)


if __name__ == "__main__":
    main()
