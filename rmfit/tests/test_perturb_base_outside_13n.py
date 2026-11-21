#!/usr/bin/env python3
"""perturb on a campaign whose base model sits OUTSIDE the campaign directory (2026-11-20).

perturb linked data/ into the new campaign from the campaign directory only; when the base .azr lives elsewhere (e.g.
rungs/model.azr with rungs/data, campaign camp/) nothing was linked and AZURE2 stopped with
"Could Not Fill Segment #1 from file".  Here: 13N with the model in model/ and the campaign in camp/; perturb must
build the sibling campaign.  Run:  python3 rmfit/tests/test_perturb_base_outside_13n.py"""
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
    tmp = tempfile.mkdtemp(prefix="rmfit_perturb_outside_")
    try:
        mdir, cdir = os.path.join(tmp, "model"), os.path.join(tmp, "camp")
        shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(mdir, "data"))
        os.makedirs(os.path.join(mdir, "output"))
        mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
        mdl.set_output_dir(os.path.join(mdir, "output"))
        base = mdl.write(os.path.join(mdir, "13N.azr"))
        c = Campaign.init(cdir, base, flags={}, nshards=2, policy=dict(bg_energy=15.0, e_window=0.5), note="perturb test")
        check("campaign directory has no data/ of its own (the case under test)", not os.path.lexists(os.path.join(cdir, "data")))
        out = os.path.join(tmp, "pert")
        try:
            st = c.perturb(out, n_flips=2, seed=1, jitter=0.0, e_kick=0.0, polish_nfev=3)
            ok, detail = True, ""
        except Exception as err:  # noqa: BLE001
            ok, detail = False, f"{type(err).__name__}: {str(err)[:200]}"
        check("perturb builds the sibling campaign", ok, detail)
        if ok:
            p = Campaign(out, enforce_engine=False)
            obj = p.incumbent()[2]
            check("the sibling has a scored incumbent", obj == obj and obj > 0, str(obj))
            p.close()
        c.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all perturb checks passed")


if __name__ == "__main__":
    main()
