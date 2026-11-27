#!/usr/bin/env python3
"""theta^2 cap per J^pi and pair (BoundsPolicy.theta2_cap_jpi_pairs, 2026-11-26:
"Cap theta2 at 1" for one J^pi's channels except the elastic one).  Pure checks of BoundsPolicy.channel_cap, then default_bounds and sanity on
the 13N fixture.  Backend-agnostic (pyazr 2 and 3)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
import numpy as np
from rmfit.local import BoundsPolicy, default_bounds, sanity
from rmfit.keys import ParamKey
ok = True
def check(name, cond, detail=""):
    global ok
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)
def main():
    p = BoundsPolicy(theta2_cap=3.0, theta2_cap_jpi_pairs=[["7/2+", [2, 3], 1.0]])
    k = lambda jpi, pair, kind="width": ParamKey(kind, level=f"{jpi}#1", pair=pair, L=3, S=0.5)
    check("matching J^pi and pair -> min(cap, rule)", p.channel_cap(k("7/2+", 2), 3.0) == 1.0 and p.channel_cap(k("7/2+", 3), 3.0) == 1.0)
    check("rule never raises a cap", p.channel_cap(k("7/2+", 2), 0.5) == 0.5)
    check("other pair of that J^pi unchanged", p.channel_cap(k("7/2+", 1), 3.0) == 3.0)
    check("other J^pi unchanged", p.channel_cap(k("7/2-", 2), 3.0) == 3.0 and p.channel_cap(k("5/2+", 3), 3.0) == 3.0)
    check("energy / norm keys unchanged", p.channel_cap(ParamKey("energy", level="7/2+#1"), 3.0) == 3.0 and p.channel_cap(ParamKey("norm", seg=1), 3.0) == 3.0)
    check("default (no rules) unchanged", BoundsPolicy().channel_cap(k("7/2+", 2), 3.0) == 3.0)
    from rmfit.engine import SingleSession
    from rmfit import AZURE2_ROOT
    azr = os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr")
    with SingleSession(azr, cwd=os.path.dirname(azr), flags=dict(use_long_wavelength=False, ignore_externals=True, use_brune=True), name="capjpi") as ev:
        x = np.array(ev.x0, float); gW = ev.wigner_amplitudes()
        W = [(j, kk) for j, kk in enumerate(ev.keys) if kk.kind == "width" and j in gW and gW[j] > 0]
        jpi, pair = W[0][1].jpi, W[0][1].pair
        base = dict(theta2_cap=3.0, bg_theta2_cap=3.0, bg_energy=99.0)
        b0 = default_bounds(ev, x, BoundsPolicy(**base))
        b1 = default_bounds(ev, x, BoundsPolicy(**base, theta2_cap_jpi_pairs=[[jpi, [pair], 1.0]]))
        hit = [j for j, kk in W if kk.jpi == jpi and kk.pair == pair]; miss = [j for j, kk in W if j not in hit]
        check(f"13N: {len(hit)} widths of {jpi} pair {pair} get 1/sqrt(3) x the bound",
              len(hit) > 0 and all(abs(b1.ub[j] / b0.ub[j] - 1 / np.sqrt(3)) < 1e-9 for j in hit), [(b1.ub[j], b0.ub[j]) for j in hit])
        check(f"13N: the other {len(miss)} width bounds unchanged", all(abs(b1.ub[j] - b0.ub[j]) < 1e-12 for j in miss))
        j0 = hit[0]; xs = x.copy(); xs[j0] = np.sqrt(2.0) * gW[j0]          # theta^2 = 2: inside 3, above the rule's 1
        r0 = sanity(ev, xs, BoundsPolicy(**base)); r1 = sanity(ev, xs, BoundsPolicy(**base, theta2_cap_jpi_pairs=[[jpi, [pair], 1.0]]))
        check("13N: sanity flags theta^2 = 2 only under the rule",
              not any(str(ev.keys[j0]) in v for v in r0.theta2_violations) and any(str(ev.keys[j0]) in v for v in r1.theta2_violations),
              (r0.theta2_violations, r1.theta2_violations))
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1
if __name__ == "__main__":
    sys.exit(main())
