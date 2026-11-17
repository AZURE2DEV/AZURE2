#!/usr/bin/env python3
"""Relocation of a narrow level on 13N: shift the 3/2- 3.5032 MeV level by +60 keV
in the file, polish (it must not find its way back inside a narrow window), then
relocate it to 3.5032 and check the loss is recovered."""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT); os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3 (its pyazr 3 version is test_p3_relocate_13n.py)")
    sys.exit(77)
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.engine import Sharded
from rmfit.local import BoundsPolicy, default_bounds, polish
from rmfit.levelscan import relocate
from pyazr import AzrModel
FLAGS = {"use_brune": True}
failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok: failures.append(name)
def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_reloc_"); work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True); os.makedirs(os.path.join(work, "output"))
    mdl = AzrModel.from_file(os.path.join(work, "13N.azr")); mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    lv = mdl.find(jpi="3/2-", energy=3.5032)[0]; lv.set_energy(3.5632)
    base = mdl.write(os.path.join(work, "shifted.azr"))
    policy = BoundsPolicy(bg_energy=15.0, e_window=0.5, narrow_window=dict(gamma_kev=100.0, factor=0.2, min_kev=5.0))
    sh = Sharded(base, os.path.join(tmp, "shards"), nshards=3, flags=FLAGS)
    label = [k.level for k in sh.keys if k.kind == "energy" and k.jpi == "3/2-" and abs(sh.level_energy[k.level] - 3.5632) < 1e-6][0]
    obj0 = sh.objective(sh.x0)
    f1 = polish(sh, sh.x0, bounds=default_bounds(sh, sh.x0, policy), policy=policy, max_nfev=60)
    e1 = f1.x[[j for j, k in enumerate(sh.keys) if k.kind == "energy" and k.level == label][0]]
    print(f"  shifted model: {obj0:.1f} -> polished {f1.objective:.1f}; {label} energy after the polish {e1:.4f} MeV")
    check("the narrow window kept the shifted level near 3.563 MeV", abs(e1 - 3.5632) <= 0.012, str(e1))
    f2 = relocate(sh, f1.x, label, 3.5032, policy, half_window_kev=5.0, max_nfev=60, log=print)
    e2 = f2.x[[j for j, k in enumerate(sh.keys) if k.kind == "energy" and k.level == label][0]]
    print(f"  relocated: {f1.objective:.1f} -> {f2.objective:.1f}; energy {e2:.4f} MeV ({f2.status})")
    check("relocation improves the objective a lot", f2.objective < 0.5 * f1.objective, f"{f2.objective} vs {f1.objective}")
    check("relocated energy within the window of 3.5032", abs(e2 - 3.5032) <= 0.0051, str(e2))
    sh.close(); shutil.rmtree(tmp, ignore_errors=True)
    print("FAILURES: " + ", ".join(failures) if failures else "ALL OK"); return 1 if failures else 0


def figure_check():
    """The review figure on the relocation case (data, shifted fit, relocated fit)."""
    import tempfile, shutil, os
    from rmfit.engine import SingleSession
    from rmfit.model import single_variant
    from rmfit.levelscan import review_figure
    from rmfit.keys import take
    tmp = tempfile.mkdtemp(prefix="rmfit_fig_"); work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True); os.makedirs(os.path.join(work, "output"))
    mdl = AzrModel.from_file(os.path.join(work, "13N.azr")); mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    lv = mdl.find(jpi="3/2-", energy=3.5032)[0]; lv.set_energy(3.5632)
    base = mdl.write(os.path.join(work, "shifted.azr"))
    with SingleSession(single_variant(base, "fig", tmp), cwd=tmp, name="fig", flags=FLAGS) as ev:
        x0 = ev.x0.copy(); x1 = x0.copy()
        j = [j for j, k in enumerate(ev.keys) if k.kind == "energy" and k.jpi == "3/2-" and abs(ev.level_energy[k.level] - 3.5632) < 1e-6][0]
        x1[j] = 3.5032
        out = review_figure(ev, x0, x1, 3.5032, os.path.join(tmp, "review.png"), half_width=0.08,
                            title="13N test: 3/2- shifted vs relocated", label_before="shifted", label_after="relocated")
    ok = out["path"] is not None and os.path.exists(out["path"]) and os.path.getsize(out["path"]) > 5000 and out["panels"] >= 2
    print(f"  {'ok  ' if ok else 'FAIL'}  review figure written ({out['panels']} panels)")
    shutil.rmtree(tmp, ignore_errors=True)
    return ok


if __name__ == "__main__":
    if os.environ.get("RELOCATE_FIGURE") == "1":
        sys.exit(0 if figure_check() else 1)
    sys.exit(main())
