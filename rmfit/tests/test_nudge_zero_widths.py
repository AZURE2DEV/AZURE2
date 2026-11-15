#!/usr/bin/env python3
"""bake's zero-width guard (2026-09-30): a free width at exactly 0 drops out of a re-read model; nudge to 1e-6.
Run: python3 rmfit/tests/test_nudge_zero_widths.py"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.keys import ParamKey  # noqa: E402
from rmfit.engine import nudge_zero_widths  # noqa: E402
keys = [ParamKey("width", level="1/2+#3", pair=1, L=0, S=0.5), ParamKey("width", level="1/2+#3", pair=2, L=1, S=0.5),
        ParamKey("energy", level="1/2+#3"), ParamKey("norm", seg=4)]
x = np.array([0.0, -0.3, 12.5, 0.0])
y, nudged = nudge_zero_widths(x, keys)
ok = (y[0] == 1e-6 and y[1] == -0.3 and y[2] == 12.5 and y[3] == 0.0 and nudged == ["g[1/2+#3 p1 L0 S0.5]"] and x[0] == 0.0)
print(("  ok    " if ok else "  FAIL  ") + "only the zero WIDTH is nudged, the input untouched", y, nudged)
y2, n2 = nudge_zero_widths(np.array([0.2, -0.3, 12.5, 1.0]), keys)
print(("  ok    " if not n2 else "  FAIL  ") + "nothing to nudge -> empty list")
sys.exit(0 if ok and not n2 else 1)
