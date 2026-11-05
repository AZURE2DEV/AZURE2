#!/usr/bin/env python3
"""The light particle's excitation energy, column 32 of a <levels> line (2026-11-05), through pyazr:

  * AzrModel reads a 32-column line, writes it back byte for byte, and a 31-column file stays 31 columns;
    AzrChannel.light_excitation / heavy_excitation / excitation (= e1 + e2, the threshold quantity); setting
    light_excitation appends the column, setting it to 0 drops it again;
  * Pair (the GET_PAIRS_INFO record) has light_excitation, heavy_excitation, threshold, label, and still reads
    a 19-field record of an older engine;
  * through the engine (skipped without the compiled module): the pair table reports e1 and the total, and the
    excitation moved from the heavy to the light partner of alpha' + 12C(2+) gives the same cross sections
    at the lab-frame data points.

Models: tests/light_excitation/make_models.py.

Run from anywhere:  python3 tests/pyazr/light_excitation_test.py
"""
import importlib.util
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SUITE = os.path.join(ROOT, "tests", "light_excitation")
sys.path.insert(0, ROOT)

try:
    from pyazr.azrfile import AzrModel  # noqa: E402
    from pyazr.parameters import Pair  # noqa: E402
except ImportError as e:
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(0)

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({detail})"))
    if not ok:
        failures.append(name)


E2P = 4.43891
tmp = tempfile.mkdtemp(prefix="pyazr_light_excitation_")
try:
    spec = importlib.util.spec_from_file_location("le_models", os.path.join(SUITE, "make_models.py"))
    mm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mm)
    mm.main(tmp)

    # --- AzrModel (no engine) ---------------------------------------------------------------------
    print("1. AzrModel and the 32nd column")
    for name in ("aC_light", "aC_heavy", "aC_heavy0", "cc_new"):
        path = os.path.join(tmp, name + ".azr")
        check(f"{name}.azr round-trips byte for byte", AzrModel.from_file(path).to_text() == open(path).read())
    light = AzrModel.from_file(os.path.join(tmp, "aC_light.azr"))
    p2 = [c for lv in light.levels for c in lv.channels if c.pair == 2]
    p1 = [c for lv in light.levels for c in lv.channels if c.pair == 1]
    check("pair 2 lines have 32 fields, pair 1 lines 31",
          all(len(c.to_line().split()) == 32 for c in p2) and all(len(c.to_line().split()) == 31 for c in p1))
    c = p2[0]
    check("light_excitation = e1, heavy_excitation = e2 = 0, excitation = e1 + e2",
          (c.light_excitation, c.heavy_excitation, c.excitation) == (E2P, 0.0, E2P),
          (c.light_excitation, c.heavy_excitation, c.excitation))
    check("a 31-field line has light_excitation 0", p1[0].light_excitation == 0.0)
    cc = AzrModel.from_file(os.path.join(tmp, "cc_new.azr")).levels[0].channels[0]
    check("cc_new: excitation = e1 + e2 = 2 E*", abs(cc.excitation - 2 * E2P) < 1e-12, cc.excitation)

    heavy = AzrModel.from_file(os.path.join(tmp, "aC_heavy.azr"))
    ch = heavy.levels[0].channels[0]
    before = ch.to_line()
    ch.light_excitation = 1.25
    after = ch.to_line()
    check("setting light_excitation appends the column and changes nothing else",
          after == before + " 1.25" and AzrModel.from_file(os.path.join(tmp, "aC_heavy.azr")).levels[0].channels[0].to_line() == before,
          repr(after))
    ch.light_excitation = 0
    check("setting it back to 0 drops the column (written only when nonzero)", ch.to_line() == before, repr(ch.to_line()))
    ch.light_excitation = 0
    check("setting 0 on a 31-field line leaves it alone", ch.to_line() == before and len(ch.tokens) == 31)
    ch2 = AzrModel.from_file(os.path.join(tmp, "aC_light.azr")).levels[0].channels[1]   # pair 2, 32 columns
    ch2.gamma = 12345.0
    line = ch2.to_line()
    check("editing another field of a 32-column line keeps column 32",
          line.split()[-1] == "4.43891" and line.split()[11] == "12345", line)
    ch2.light_excitation = 0
    check("dropping column 32 after an edit keeps the edit", len(ch2.to_line().split()) == 31 and
          ch2.to_line().split()[11] == "12345", ch2.to_line())
    try:
        from pyazr.azrfile import AzrChannel
        AzrChannel.from_line(" ".join(["0"] * 33))
        check("a 33-field line is refused", False)
    except ValueError:
        check("a 33-field line is refused", True)

    # clone / write / reread keeps e1
    m = AzrModel.from_file(os.path.join(tmp, "aC_light.azr"))
    out = os.path.join(tmp, "rewritten.azr")
    m.add_level(J=2.0, parity=1, energy=15.0, channels=[dict(pair=2, L=0, S=2.0, gamma=10.0)])
    m.write(out)
    added = [c for lv in AzrModel.from_file(out).levels if abs(lv.energy - 15.0) < 1e-9 for c in lv.channels]
    check("add_level: the new level's pair-2 lines carry the light excitation, its pair-1 line none",
          any(a.pair == 2 for a in added) and all(a.light_excitation == (E2P if a.pair == 2 else 0.0) for a in added),
          [a.to_line() for a in added])

    # --- Pair (no engine) -------------------------------------------------------------------------
    print("2. Pair records")
    rec = [2, 2, 0, 0, 2.0, 1, 2, 4.0026, 0.0, 1, 6, 12.0, 7.16191, E2P, 5.43, 0.2, 8, 16, 1.0, E2P]
    p = Pair.from_record(rec)
    check("a 20-field record: light_excitation, heavy_excitation, threshold",
          p.light_excitation == E2P and p.heavy_excitation == 0.0 and abs(p.threshold - (7.16191 + E2P)) < 1e-12,
          (p.light_excitation, p.heavy_excitation, p.threshold))
    check("label marks the excited light particle", p.label == "2+*(4.439 MeV) + 0+", p.label)
    check("repr shows light_excitation", "light_excitation=4.43891" in repr(p), repr(p))
    old = Pair.from_record(rec[:19])
    check("a 19-field record (older engine) reads, light_excitation = 0, heavy = total",
          old.light_excitation == 0.0 and old.heavy_excitation == E2P, (old.light_excitation, old.heavy_excitation))
    check("Pair._NFIELDS is 20", Pair._NFIELDS == 20)

    # --- engine ----------------------------------------------------------------------------------
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import numpy as np
        from pyazr import azure2
        from pyazr import _azure2  # noqa: F401
    except Exception as err:
        print(f"skip engine part: {type(err).__name__}: {err}")
    else:
        print("3. through the engine")

        def session(name):
            work = os.path.join(tmp, "run_" + name)
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))
            shutil.copy(os.path.join(tmp, name + ".azr"), work)
            os.symlink(os.path.join(tmp, "data"), os.path.join(work, "data"))
            return azure2(os.path.join(work, name + ".azr"), cwd=work)

        res = {}
        for name in ("aC_heavy", "aC_light"):
            with session(name) as s:
                q = s.pairs[1]
                res[name] = ([np.asarray(v, float) for v in s.calculate_rwa(np.asarray(s.params_rwa, float))],
                             q.light_excitation, q.excitation)
        check("aC_light pair 2: light_excitation = E*, aC_heavy pair 2: 0; both excitation = E*",
              res["aC_light"][1] == E2P and res["aC_heavy"][1] == 0.0 and res["aC_light"][2] == res["aC_heavy"][2] == E2P,
              (res["aC_light"][1:], res["aC_heavy"][1:]))
        a, b = res["aC_heavy"][0], res["aC_light"][0]
        worst = max(float(np.max(np.abs(x - y) / np.maximum(np.abs(x), 1e-300))) for x, y in zip(a, b))
        check("calculate_rwa: excitation on the light vs the heavy partner, same values (1e-10)",
              len(a) == len(b) == 2 and worst <= 1e-10, worst)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: pyazr light_excitation")
