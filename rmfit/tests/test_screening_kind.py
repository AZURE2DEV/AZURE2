#!/usr/bin/env python3
"""rmfit with fitted electron-screening potentials (screening_<n>, 'screening Ue vary sigma' on a <targetInt> line).
Until 2026-11-05 build_keys stopped at this kind ("unknown parameter kind 'screening'").  Model: AZURE2's
tests/electron_screening/grad.azr (12C(p,g), LUNA segments, U_e = 1000 eV free with a 200 eV prior; data from
tests/13N/data).
  (a) build_keys gives screening[t1];
  (b) at U_e = 1100 eV SingleSession's penalty includes ((1100 - 1000)/200)^2 = 0.25 and its objective equals
      pyazr's objective() (same flags) to 1e-9;
  (c) default_bounds: U_e in [0, inf);
  (d) a shard map is refused (SingleSession only).
    python3 rmfit/tests/test_screening_kind.py      (installed pyazr; about 20 s)"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit import AZURE2_ROOT  # noqa: E402
from rmfit.engine import SingleSession  # noqa: E402
from rmfit.keys import build_keys  # noqa: E402
from rmfit.local import BoundsPolicy, default_bounds  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    from pyazr import azure2
    src = os.path.join(AZURE2_ROOT, "tests", "electron_screening", "grad.azr")
    if not os.path.exists(src):
        print(f"SKIP: no {src}")
        return 2
    tmp = tempfile.mkdtemp(prefix="rmfit_screening_")
    try:
        for d in ("output", "checks"):
            os.makedirs(os.path.join(tmp, d))
        shutil.copy(src, tmp)
        os.symlink(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
        azr = os.path.join(tmp, "grad.azr")
        with SingleSession(azr, cwd=tmp, name="scr") as ev:
            labels = [str(k) for k in ev.keys]
            check("(a) key screening[t1]", "screening[t1]" in labels, labels[-4:])
            j = labels.index("screening[t1]")
            x = ev.x0.copy()
            x[j] = 1100.0
            r0 = ev.penalties.value(ev.x0.copy())
            pen = float(ev.penalties.value(x))
            check("(b) penalty includes ((1100 - 1000)/200)^2", abs((pen - r0) - 0.25) < 1e-9, f"{pen} - {r0}")
            obj = float(ev.objective(x))
            with azure2(azr, cwd=tmp, **ev.flags) as s:
                ref = float(s.objective(x))
            check("(b) objective = pyazr's objective()", abs(obj - ref) <= 1e-9 * abs(ref), f"{obj} vs {ref}")
            b = default_bounds(ev, ev.x0.copy(), BoundsPolicy())
            check("(c) U_e bounds [0, inf)", b.lb[j] == 0.0 and np.isinf(b.ub[j]) and not b.fixed[j], (b.lb[j], b.ub[j], b.fixed[j]))
            with azure2(azr, cwd=tmp, **ev.flags) as s:
                try:
                    build_keys(s, seg_map={k: k for k in range(1, 20)})
                    check("(d) a shard map is refused", False, "no exception")
                except NotImplementedError as e:
                    check("(d) a shard map is refused", "screening" in str(e) and "SingleSession only" in str(e), str(e))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("PASS" if not failures else f"FAIL: {failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
