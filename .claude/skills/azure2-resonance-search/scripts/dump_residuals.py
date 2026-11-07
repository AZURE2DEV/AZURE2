#!/usr/bin/env python3
"""Point-by-point residuals of a merged-fit export (single session, parameters from its .sav BY NAME).
usage: dump_residuals.py <fit.azr> <fit.sav> <workdir> <out.npz>
Saves, per data row: segment key, dataset name, E_cm, Ex, c.m. angle, data, data error, fit, normalized residual
(the residual AZURE2/rmfit put into chi^2, i.e. with the fitted norm), so the ranking needs no re-evaluation."""
import sys, os, json
import numpy as np
from rmfit.engine import SingleSession
from rmfit.model import single_variant
from rmfit.diagnose import dataset_name
FLAGS = dict(use_long_wavelength=False, ignore_externals=True, use_brune=True)
def main():
    azr, sav, wd, out = sys.argv[1:5]
    os.makedirs(wd, exist_ok=True)
    if not os.path.exists(os.path.join(wd, "data")):
        os.symlink(os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(azr)), "data")), os.path.join(wd, "data"))
    vals = {}
    for line in open(sav):
        t = line.split()
        if len(t) >= 2:
            try: vals[t[0]] = float(t[1])
            except ValueError: pass
    with SingleSession(single_variant(azr, "res", wd), cwd=wd, flags=FLAGS, name="res") as ev:
        x = np.array(ev.x0, float)
        for p in ev.worker.call("snapshot", None)["parameters"]:
            if not p.fixed and p.free_index is not None and p.name in vals: x[p.free_index] = vals[p.name]
        s = ev.score(x); print(f"objective {s.objective:.3f} chi2 {s.chi2:.3f} N {s.npts}", flush=True)
        r, _ = ev.residuals(x, jac=False); res = ev._plain_data_rows(r)
        d = ev.data(); fit = ev.cross(x)
        cols = dict(key=[], name=[], ecm=[], ex=[], ang=[], dat=[], err=[], fit=[], res=[])
        for i, (key, (a, b)) in enumerate(ev.row_ranges):
            n = b - a; seg = ev.segments.by_key(key)
            ang = d["angles"][i] if i < len(d["angles"]) and len(d["angles"][i]) == n else np.full(n, np.nan)
            assert len(d["excitation"][i]) == n and len(fit[i]) == n, (key, n, len(d["excitation"][i]), len(fit[i]))
            cols["key"] += [key] * n; cols["name"] += [dataset_name(seg)] * n
            cols["ecm"] += list(d["energies"][i]); cols["ex"] += list(d["excitation"][i]); cols["ang"] += list(ang)
            cols["dat"] += list(d["cross"][i]); cols["err"] += list(d["cross_err"][i]); cols["fit"] += list(fit[i])
            cols["res"] += list(res[a:b])
        np.savez(out, **{k: np.asarray(v) for k, v in cols.items()})
        print(f"rows {len(cols['res'])}, sum res^2 {np.sum(np.square(cols['res'])):.3f}", flush=True)

if __name__ == "__main__":
    main()
