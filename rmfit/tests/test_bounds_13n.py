#!/usr/bin/env python3
"""The narrow-level energy-window guard (BoundsPolicy.narrow_window) on 13N."""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.engine import SingleSession
from rmfit.local import BoundsPolicy, default_bounds
from rmfit.model import single_variant
from pyazr import AzrModel
failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok: failures.append(name)
def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_bounds_"); work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True); os.makedirs(os.path.join(work, "output"))
    base = AzrModel.from_file(os.path.join(work, "13N.azr")).write(os.path.join(work, "base.azr"))
    with SingleSession(single_variant(base, "b", tmp), cwd=tmp, name="b") as ev:
        x = ev.x0.copy()
        e_keys = [(j, k) for j, k in enumerate(ev.keys) if k.kind == "energy"]
        b0 = default_bounds(ev, x, BoundsPolicy(bg_energy=15.0, e_window=0.5, threshold_margin_mev=0.0))
        bT = default_bounds(ev, x, BoundsPolicy(bg_energy=15.0, e_window=0.5))   # default 3 keV threshold margin
        b1 = default_bounds(ev, x, BoundsPolicy(bg_energy=15.0, e_window=0.5, narrow_window=dict(gamma_kev=100.0, factor=3.0, min_kev=1.0)))
        phys = np.asarray(ev.transform(x), float)
        widths = {}
        for j, k in enumerate(ev.keys[:ev.n_rmatrix]):
            if k.kind == "width": widths[k.level] = widths.get(k.level, 0.0) + abs(phys[j]) * 1e-3
        print("  level total widths (keV):", {k: round(v, 1) for k, v in widths.items()})
        narrow = [(j, k) for j, k in e_keys if not b1.fixed[j] and widths.get(k.level, 1e9) < 100.0]
        wide = [(j, k) for j, k in e_keys if not b1.fixed[j] and widths.get(k.level, 1e9) >= 100.0]
        T = 1.94351   # 12C + p threshold of the 13N fixture
        check("threshold guard: every free energy window stays >= 3 keV from the p threshold, on its own side",
              all((bT.lb[j] >= T + 0.003 - 1e-9 if x[j] >= T else bT.ub[j] <= T - 0.003 + 1e-9) for j, k in e_keys if not bT.fixed[j]))
        check("threshold guard: windows that do not touch the threshold are unchanged",
              all(abs(bT.ub[j] - b0.ub[j]) < 1e-12 and abs(bT.lb[j] - b0.lb[j]) < 1e-12 for j, k in e_keys
                  if not bT.fixed[j] and not (b0.lb[j] < T + 0.003 and b0.ub[j] > T - 0.003)))
        check("threshold guard: at least one 13N window is actually cut (the test has teeth)",
              any(abs(bT.ub[j] - b0.ub[j]) > 1e-9 or abs(bT.lb[j] - b0.lb[j]) > 1e-9 for j, k in e_keys if not bT.fixed[j]))
        check("default policy unchanged (+-0.5 MeV) for free energies (margin off)", all(abs((b0.ub[j] - b0.lb[j]) - 1.0) < 1e-9 for j, k in e_keys if not b0.fixed[j]))
        check("narrow levels get +-max(1 keV, 3 Gamma)", narrow and all(abs((b1.ub[j] - b1.lb[j]) - 2 * max(1.0, 3 * widths[k.level]) * 1e-3) < 1e-9 for j, k in narrow),
              str([(k.level, b1.ub[j] - b1.lb[j]) for j, k in narrow]))
        check("wide levels keep the full window", all(abs((b1.ub[j] - b1.lb[j]) - 1.0) < 1e-9 for j, k in wide))
        check("fixed flags identical", np.array_equal(b0.fixed, b1.fixed))
    shutil.rmtree(tmp, ignore_errors=True)
    print("FAILURES: " + ", ".join(failures) if failures else "ALL OK"); return 1 if failures else 0
if __name__ == "__main__":
    sys.exit(main())
