#!/usr/bin/env python3
"""rmfit with the parameter kinds of 2026-11-05: a fitted free-gas Doppler temperature (doppler_<n>, 'doppler T vary
sigma [Tmax]' on a <targetInt> line) and fitted sample thicknesses (segment_<key>_thickness, 'transmission|yield n
vary sigma').  Model: AZURE2's tests/doppler_fgm dopT (n + 12C, a transmission and a yield segment under a
Doppler-only line), made the fit forms 'doppler 300 1 20 900', 'transmission 0.15 1 0.01', 'yield 0.15 1 0.01'.
  (a) build_keys gives doppler[t1], thickness[seg1], thickness[seg2];
  (b) at a moved vector, SingleSession's penalty = ((T-300)/20)^2 + sum ((n-0.15)/0.01)^2, and its objective equals
      pyazr's objective() (same flags) to 1e-9;
  (c) default_bounds: T in [75, 900] K, n in [1.5e-7, inf);
  (d) a shard map is refused for both kinds (SingleSession only).
Needs a pyazr that knows kinds 9 and 10:
    AZURE2_ROOT=<build> python3 rmfit/tests/test_doppler_thickness_kinds.py      (about 20 s)"""
import importlib.util
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
    from pyazr import parameters as P
    if "thickness" not in getattr(P, "_KINDS", {}).values():
        print(f"SKIP: the pyazr at AZURE2_ROOT={AZURE2_ROOT} knows no 'thickness' kind (needs AZURE2 >= 2026-11-05)")
        return 2
    tmp = tempfile.mkdtemp(prefix="rmfit_kinds_")
    try:
        spec = importlib.util.spec_from_file_location("dfgm", os.path.join(AZURE2_ROOT, "tests", "doppler_fgm", "make_models.py"))
        mm = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mm)
        mm.main(os.path.join(AZURE2_ROOT, "tests", "reaction_rate", "models", "n12C.azr"), tmp)
        text = open(os.path.join(tmp, "dopT.azr")).read()
        for old, new in (("doppler 300\n", "doppler 300 1 20 900\n"), ("transmission 0.15\n", "transmission 0.15 1 0.01\n"),
                         ("yield 0.15\n", "yield 0.15 1 0.01\n")):
            assert text.count(old) == 1, old
            text = text.replace(old, new)
        azr = os.path.join(tmp, "kinds.azr")
        open(azr, "w").write(text)
        os.makedirs(os.path.join(tmp, "output"), exist_ok=True)
        os.makedirs(os.path.join(tmp, "checks"), exist_ok=True)

        with SingleSession(azr, cwd=tmp, name="kinds") as ev:
            labels = [str(k) for k in ev.keys]
            check("(a) keys doppler[t1], thickness[seg1], thickness[seg2]",
                  {"doppler[t1]", "thickness[seg1]", "thickness[seg2]"} <= set(labels), labels[-6:])
            j = {s: labels.index(s) for s in ("doppler[t1]", "thickness[seg1]", "thickness[seg2]")}
            x = ev.x0.copy()
            x[j["doppler[t1]"]] = 330.0
            x[j["thickness[seg1]"]] = 0.16
            x[j["thickness[seg2]"]] = 0.17
            want = ((330.0 - 300.0) / 20.0) ** 2 + ((0.16 - 0.15) / 0.01) ** 2 + ((0.17 - 0.15) / 0.01) ** 2
            pen = float(ev.penalties.value(x))
            check("(b) penalty = the three priors", abs(pen - want) < 1e-9 * want, f"{pen} vs {want}")
            obj = float(ev.objective(x))
            with azure2(azr, cwd=tmp, **ev.flags) as s:
                ref = float(s.objective(x))
            check("(b) objective = pyazr's objective()", abs(obj - ref) <= 1e-9 * abs(ref), f"{obj} vs {ref}")
            b = default_bounds(ev, ev.x0.copy(), BoundsPolicy())
            jd, j1 = j["doppler[t1]"], j["thickness[seg1]"]
            check("(c) doppler bounds [75, 900] K", abs(b.lb[jd] - 75.0) < 1e-9 and abs(b.ub[jd] - 900.0) < 1e-9 and not b.fixed[jd],
                  (b.lb[jd], b.ub[jd], b.fixed[jd]))
            check("(c) thickness bounds [1.5e-7, inf)", abs(b.lb[j1] - 1.5e-7) < 1e-15 and np.isinf(b.ub[j1]) and not b.fixed[j1],
                  (b.lb[j1], b.ub[j1], b.fixed[j1]))
        # (d) one model per kind, so each refusal is reached: the Doppler fit form alone, then the thicknesses alone
        only = {"doppler": text.replace("transmission 0.15 1 0.01\n", "transmission 0.15\n").replace("yield 0.15 1 0.01\n", "yield 0.15\n"),
                "thickness": text.replace("doppler 300 1 20 900\n", "doppler 300\n")}
        for kind, t in only.items():
            path = os.path.join(tmp, f"only_{kind}.azr")
            open(path, "w").write(t)
            with azure2(path, cwd=tmp) as s:
                kinds = {p.kind for p in s.parameters if not p.fixed}
                try:
                    build_keys(s, seg_map={1: 1, 2: 2})
                    check(f"(d) a shard map is refused ({kind})", False, "no exception")
                except NotImplementedError as e:
                    check(f"(d) a shard map is refused ({kind})", kind in kinds and "SingleSession only" in str(e)
                          and ("Doppler" in str(e)) == (kind == "doppler"), str(e))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("PASS" if not failures else f"FAIL: {failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
