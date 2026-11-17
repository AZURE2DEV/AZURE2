#!/usr/bin/env python3
"""bake writes the FIXED parameters' current physical values too (2026-10-24), with the real engine on 13N.
A fixed parameter is fixed in reduced-width space; its observed value moves with the other amplitudes of its level (Brune's
1 + sum gamma^2 dS/dE).  Until 2026-10-24 bake wrote only the free parameters, so a fixed width kept the seed file's observed
value and the written file was another model (one archived candidate: 14 fixed widths of frozen
levels whose free ANCs had moved, up to +13 %, export refused).
Model: 13N with the E1 photon width of the 3/2- 3.5032 MeV level FIXED (-0.49 eV), its proton width free; the vector scales the
proton amplitude by 1.5.
  - bake passes its fresh-session check (objective and <levels>-only chi2);
  - the written fixed photon width is the engine's current observed value, not the seed's -0.49 eV;
  - the written file's rwa for that width equals the source's (reads back the same amplitude);
  - the old bake (free parameters only) fails the same check (only if engine.py.bak-2026-10-24 is beside this package).
Run: python3 rmfit/tests/test_bake_fixed_13n.py   (a few minutes)"""
import os, sys, shutil, tempfile, io, contextlib, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
PKG_PARENT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PKG_PARENT)
sys.path.insert(0, os.environ.get("AZURE2_ROOT", os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))))
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if not _PYAZR3:
    print("SKIP: a pyazr 3 test; this process runs pyazr 2 (the pyazr 2 version is test_bake_fixed_13n.py)")
    sys.exit(77)
import numpy as np
from rmfit.azr3 import Model as AzrModel
from rmfit.engine import SingleSession, bake
from rmfit.model import single_variant
AZ = os.environ.get("AZURE2_ROOT", os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
failures = []
FLAGS = dict(use_brune=True)


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def build(tmp):
    # fix the 3/2- 3.5032 E1 photon channel (pair 2, l = 2, -0.49 eV)
    mdl = AzrModel.from_file(os.path.join(AZ, "tests", "13N", "13N.azr"))
    lv = mdl.find(jpi="3/2-", energy=3.5032)[0]
    c = next(c for c in lv.channels if c.pair == 2 and abs(c.gamma + 0.49) < 1e-12)
    c.channel_fixed = True
    mdl.set_output_dir(os.path.join(tmp, "output"))
    return mdl.write(os.path.join(tmp, "13N_fixed.azr"))


def fixed_value(azr):
    for l in AzrModel.from_file(azr).levels:
        if l.jpi == "3/2-" and abs(l.energy - 3.5032) < 1e-3:
            for c in l.channels:                    # the channel build() fixed (pair 2, the only fixed nonzero one)
                if c.pair == 2 and c.channel_fixed and c.gamma != 0.0:
                    return float(c.gamma)
    raise AssertionError("3/2- 3.5032 photon line not found in " + azr)


def fixed_rwa(azr, tmp, tag):
    from pyazr import Session
    with Session(azr, cwd=tmp, **FLAGS) as m:
        for p in m.parameters:
            if p.kind == "width" and p.fixed and p.start != 0.0 and abs((p.level_energy or 0) - 3.5032) < 1e-2 \
                    and p.pair_key == 2:
                return float(p.start)
    return None


def vector(base, tmp):
    with SingleSession(single_variant(base, "v", tmp), cwd=tmp, flags=FLAGS, name="v") as ev:
        x = np.array(ev.x0, float); keys = list(ev.keys)
    lvl = [k.level for i, k in enumerate(keys) if k.kind == "energy" and abs(x[i] - 3.5032) < 1e-6]
    assert len(lvl) == 1, lvl
    j = [i for i, k in enumerate(keys) if k.kind == "width" and k.level == lvl[0] and k.pair == 1]
    assert len(j) == 1, [str(k) for k in keys]
    x[j[0]] *= 1.5
    return x, keys


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_bakefix_")
    shutil.copytree(os.path.join(AZ, "tests", "13N", "data"), os.path.join(tmp, "data")); os.makedirs(os.path.join(tmp, "output"))
    base = build(tmp)
    check("model: the photon width is fixed at -0.49 eV", abs(fixed_value(base) + 0.49) < 1e-12, str(fixed_value(base)))
    x, keys = vector(base, tmp)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            azr, sav, sc = bake(base, x, keys, os.path.join(tmp, "baked.azr"), FLAGS, workdir=tmp)
        rt = sc.roundtrip
        check("bake passes its fresh-session check", rt["objective_rel"] < 1e-6 and rt["levels_rel"] < 1e-6, str(rt))
        v = fixed_value(azr)
        check("the written fixed width moved from the seed's -0.49 eV", abs(v + 0.49) > 1e-3 * 0.49, str(v))
        r0, r1 = fixed_rwa(base, tmp, "src"), fixed_rwa(azr, tmp, "baked")
        check("the written file reads back the same fixed amplitude", r0 is not None and r1 is not None
              and abs(r1 - r0) <= 1e-7 * abs(r0), f"{r0} vs {r1}")
        print(f"   fixed photon width: seed -0.49 eV -> written {v:.6g} eV; rwa {r0:.9g} / {r1:.9g}; check {rt}")
    except Exception as e:
        check("bake passes its fresh-session check", False, f"{e}\n{buf.getvalue()}")
    old = os.path.join(HERE, "..", "engine.py.bak-2026-10-24")
    if os.path.exists(old):
        import importlib.machinery
        loader = importlib.machinery.SourceFileLoader("rmfit._engine_old", old)    # not a .py name: an explicit loader
        spec = importlib.util.spec_from_loader("rmfit._engine_old", loader)
        mod = importlib.util.module_from_spec(spec); mod.__package__ = "rmfit"; loader.exec_module(mod)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                mod.bake(base, x, keys, os.path.join(tmp, "baked_old.azr"), FLAGS, workdir=tmp)
            check("control: the old bake (free parameters only) fails", False, "it passed: the test does not discriminate")
        except RuntimeError as e:
            check("control: the old bake (free parameters only) fails", "verification failed" in str(e), str(e))
            print("   old bake: " + str(e).splitlines()[0])
    shutil.rmtree(tmp, ignore_errors=True)
    print("ALL PASSED" if not failures else f"FAILED: {failures}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
