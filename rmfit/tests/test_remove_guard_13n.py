#!/usr/bin/env python3
"""run_remove's "last level of a J^pi" guard counts every active level of the group, fixed ones included (2026-11-20).

On a windowed model the levels above the window are fixed, so a level added inside the window was the only one of its
J^pi with FREE widths, and run_remove kept it as "the only level of 6+ (it carries the group's hard-sphere term)" --
a spurious narrow 6+ added in a benchmark (frozen removal cost well below the threshold).  Here: 13N
with the 20 MeV 5/2+ background level fixed, so the 5/2+ 3.5453 MeV level is the group's only free one; run_remove must
not call it "the only level of 5/2+".  The 3/2+ group has a single level: that one must still be kept.
Run:  python3 rmfit/tests/test_remove_guard_13n.py"""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3")
    sys.exit(77)
from rmfit import AZURE2_ROOT
from rmfit.campaign import Campaign
from pyazr import AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_remove_guard_")
    try:
        shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
        os.makedirs(os.path.join(tmp, "output"))
        mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
        for lv in mdl.levels:
            if lv.jpi == "5/2+" and float(lv.energy) > 10:
                lv.set_fixed(True)
                for c in lv.channels:
                    c.channel_fixed = True
        mdl.set_output_dir(os.path.join(tmp, "output"))
        base = mdl.write(os.path.join(tmp, "13N.azr"))
        c = Campaign.init(tmp, base, flags={}, nshards=2, policy=dict(bg_energy=15.0, e_window=0.5), note="remove guard test")
        lines = []
        c.log = lambda s, *a, **k: lines.append(str(s))
        c.run_remove(top=12, nfev=10)
        text = "\n".join(lines)
        only = [l for l in lines if "the only level of" in l]
        check("5/2+ 3.5453 (only FREE 5/2+; the 20 MeV 5/2+ is fixed) is not called the only level of 5/2+",
              not any("of 5/2+" in l for l in only), "; ".join(only))
        check("the single 3/2+ level is still kept as the only level of its group", any("of 3/2+" in l for l in only)
              or "3/2+" not in text, "; ".join(only))
        c.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all remove-guard checks passed")


if __name__ == "__main__":
    main()
