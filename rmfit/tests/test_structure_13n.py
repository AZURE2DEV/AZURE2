#!/usr/bin/env python3
"""Window variants, level add/deactivate variants and the pinned J^pi scan on
the 13N model.  Run:  python3 rmfit/tests/test_structure_13n.py
"""
import os, shutil, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.engine import SingleSession
from rmfit.keys import take
from rmfit.windows import window_variant, ladder_variants
from rmfit.structure import (LevelCandidate, add_level_variant, deactivate_level_variant,
                             scan_jpi_on_window, frozen_removal_screen)
from pyazr import AzrModel

failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok: failures.append(name)

def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_struct_")
    work = os.path.join(tmp, "13N"); shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True); os.makedirs(os.path.join(work, "output"))
    base = os.path.join(work, "13N.azr")
    with SingleSession(base, cwd=work, flags={}) as ev:
        keys0, x0, n0 = ev.keys, ev.x0.copy(), ev.npts
        o0 = ev.objective(x0)
    print("1. window variant")
    win = os.path.join(work, "win.azr")
    info = window_variant(base, win, 1.0, 3.0, level_margin=0.3, output_dir=os.path.join(work, "output_win"))
    print("  ", {k: v for k, v in info.items() if k != "path"})
    with SingleSession(win, cwd=work, flags={}) as ev:
        check("window keys are a subset of the base keys", set(ev.keys) <= set(keys0), str(set(ev.keys) - set(keys0)))
        check("window has fewer points", 0 < ev.npts < n0, f"{ev.npts} vs {n0}")
        xw = take(x0, keys0, ev.keys, ev.x0)
        ow = ev.objective(xw)
        check("window objective is finite and below the full one", np.isfinite(ow) and ow < o0, f"{ow} vs {o0}")
    print("2. ladder variants")
    lad = ladder_variants(base, os.path.join(tmp, "ladder"), cutoffs=[2.0, 3.0, 5.0], ex_lo=0.0)
    check("three ladder steps", len(lad) == 3 and all(os.path.exists(l["path"]) for l in lad))
    print("  ", [(l["cutoff"], l["segments_kept"], len(l["levels_fixed"])) for l in lad])
    print("3. add / deactivate level variants")
    cand = LevelCandidate(energy=3.0, channels=[dict(pair=1, L=0, S=0.5, gamma=20000.0)], source="test")
    added = add_level_variant(base, os.path.join(work, "add.azr"), cand, "1/2+", output_dir=os.path.join(work, "output_add"))
    with SingleSession(added["path"], cwd=work, flags={}) as ev:
        new = [k for k in ev.keys if k.level == added["level"]]
        check("the added level has free width keys", len(new) >= 1, str(added))
        check("all base keys survive the addition", set(keys0) <= set(ev.keys))
        xa = take(x0, keys0, ev.keys, ev.x0)
        oa = ev.objective(xa)
        print(f"   objective with the seeded level: {oa:.3f} vs {o0:.3f}")
    deact = deactivate_level_variant(base, os.path.join(work, "deact.azr"), "5/2+", 3.5453, output_dir=os.path.join(work, "output_deact"))
    with SingleSession(deact["path"], cwd=work, flags={}) as ev:
        gone = [k for k in keys0 if k.level == "5/2+#1"]
        check("deactivated level's widths leave the free vector", not any(k in ev.keys for k in gone), str(gone))
    print("4. frozen removal screen (base)")
    with SingleSession(base, cwd=work, flags={}) as ev:
        scr = frozen_removal_screen(ev, x0, log=print)
        check("screen ranks every level with a nonzero width", len(scr) >= 3)
    print("5. pinned J^pi scan on a window (3 groups)")
    t0 = time.time()
    res = scan_jpi_on_window(base, os.path.join(tmp, "scan"), cand, {}, jpis=["1/2+", "3/2-", "5/2+"], nshards=3,
                             half_width=0.6, x_base=x0, base_keys=keys0, max_nfev=25, log=print)
    print(f"   {time.time()-t0:.0f}s; " + "; ".join(f"{r.jpi}:{r.signs}:{r.objective:.1f}({r.status})" for r in res))
    check("scan returns one result per J^pi", len(res) == 3)
    check("at least one J^pi polished", any(r.status != "crashed" for r in res))
    shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all structure checks passed")

if __name__ == "__main__":
    main()
