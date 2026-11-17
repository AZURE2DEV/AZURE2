#!/usr/bin/env python3
"""width_bound_pins (2026-09-28): amplitudes on a WIDTH bound (below the theta^2-cap bound) are the bound
doing its job, not absorbers; only cap pins count as unphysical.  Run: python3 rmfit/tests/test_width_bound_pins.py"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.keys import ParamKey  # noqa: E402
from rmfit.local import width_bound_pins  # noqa: E402
failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok: failures.append(name)
keys = [ParamKey("width", level="3/2+#9", pair=1, L=2, S=0.5), ParamKey("width", level="3/2+#9", pair=2, L=1, S=0.5),
        ParamKey("energy", level="3/2+#9"), ParamKey("width", level="5/2-#2", pair=2, L=3, S=0.5)]
gW = {0: 2.0, 1: 0.3, 3: 2.0}
cap = 0.5
ub = np.array([np.sqrt(0.04) * 2.0, np.sqrt(0.5) * 0.3, np.inf, np.sqrt(0.5) * 2.0])   # key 0 bounded by a width bound, key 1 by the cap
pins = width_bound_pins(keys, [0, 1, 2], ub, gW, cap)
check("key on the width bound is listed", pins == ["g[3/2+#9 p1 L2 S0.5]"], pins)
check("key on the cap bound is not", "g[3/2+#9 p2 L1 S0.5]" not in pins)
check("other level ignored", width_bound_pins(keys, [3], ub, gW, cap) == [])
check("no width bound -> empty", width_bound_pins(keys, [0, 1], np.array([np.sqrt(0.5) * 2.0, np.sqrt(0.5) * 0.3, np.inf, 1.0]), gW, cap) == [])
print(); 
if failures: print(f"FAILED: {failures}"); sys.exit(1)
print("all width-bound-pin checks passed")
