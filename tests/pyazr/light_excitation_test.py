#!/usr/bin/env python3
"""The light particle's excitation energy, column 32 of a <levels> line (2026-11-05), through pyazr:

  * pyazr.Project reads a 32-column row and gives the file back byte for byte; ProjectPair.light_excitation /
    heavy_excitation (e1, e2); Project.set_pair(light_excitation=) writes the column on every row of the pair
    (appended where a row has none, nothing else of the row changed), 0 on a row without it changes nothing; a new
    level (add_level) takes the light excitation of its pair's rows;
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
    from pyazr import Project  # noqa: E402
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

    # --- Project ------------------------------------------------------------------------------------
    print("1. Project and the 32nd column")
    for name in ("aC_light", "aC_heavy", "aC_heavy0", "cc_new"):
        path = os.path.join(tmp, name + ".azr")
        check(f"{name}.azr round-trips byte for byte", Project.open(path).to_bytes() == open(path, "rb").read())
    light = Project.open(os.path.join(tmp, "aC_light.azr"))
    lines = light.text.split("\n")
    rows = lambda p, key: [lines[c.line - 1] for lv in p.levels() for c in lv.channels if c.pair == key]
    check("pair 2 rows have 32 fields, pair 1 rows 31",
          all(len(r.split()) == 32 for r in rows(light, 2)) and all(len(r.split()) == 31 for r in rows(light, 1)))
    q = {p.key: p for p in light.pairs()}
    check("light_excitation = e1, heavy_excitation = e2 = 0",
          (q[2].light_excitation, q[2].heavy_excitation) == (E2P, 0.0), (q[2].light_excitation, q[2].heavy_excitation))
    check("a 31-field row has light_excitation 0", q[1].light_excitation == 0.0)
    cc = Project.open(os.path.join(tmp, "cc_new.azr")).pairs()[0]
    check("cc_new: e1 + e2 = 2 E*", abs(cc.light_excitation + cc.heavy_excitation - 2 * E2P) < 1e-12,
          (cc.light_excitation, cc.heavy_excitation))

    # the nuclear part alone (no data set names a pair, so a pair whose state changes stays valid)
    bare_heavy = Project.new(os.path.join(tmp, "aC_heavy_bare.azr"), nuclear_part=os.path.join(tmp, "aC_heavy.azr"))
    bare_heavy.save()
    heavy = Project.open(os.path.join(tmp, "aC_heavy_bare.azr"))
    key = heavy.levels()[0].channels[0].pair
    before = rows(heavy, key)
    heavy.set_pair(key, light_excitation=1.25)
    lines = heavy.text.split("\n")
    after = rows(heavy, key)
    check("set_pair(light_excitation=) appends the column on every row of the pair and changes nothing else",
          after == [r + " 1.25" for r in before] and {p.key: p for p in heavy.pairs()}[key].light_excitation == 1.25, repr(after[:1]))
    heavy.set_pair(key, light_excitation=0)
    lines = heavy.text.split("\n")
    check("setting it back to 0 writes 0 in the column", all(r.split()[-1] == "0" and len(r.split()) == 32
                                                          for r in rows(heavy, key)), repr(rows(heavy, key)[:1]))
    plain = Project.open(os.path.join(tmp, "aC_heavy_bare.azr"))
    plain.set_pair(key, light_excitation=0)
    check("setting 0 on 31-field rows leaves them alone",
          plain.to_bytes() == open(os.path.join(tmp, "aC_heavy_bare.azr"), "rb").read())
    try:
        Project.open(os.path.join(tmp, "aC_heavy_bare.azr")).set_pair(key, light_excitation=-1.0)
        check("a negative excitation is refused", False)
    except ValueError:
        check("a negative excitation is refused", True)

    # a new level of pair 2 carries the light excitation, its pair-1 row none
    m = Project.new(os.path.join(tmp, "aC_light_bare.azr"), nuclear_part=os.path.join(tmp, "aC_light.azr"))
    m.add_level(2.0, +1, 15.0, [dict(pair=2, l=0, s=2.0, width=10.0)])
    out = m.save(os.path.join(tmp, "rewritten.azr"))
    back = Project.open(out)
    lines = back.text.split("\n")
    added = [(c.pair, lines[c.line - 1]) for c in back.level(("2+", 15.0)).channels]
    check("add_level: the new level's pair-2 rows carry the light excitation, its pair-1 row none",
          any(k == 2 for k, _ in added)
          and all((r.split()[-1] == "4.43891" and len(r.split()) == 32) if k == 2 else len(r.split()) == 31
                  for k, r in added), added)

    # --- Pair (no engine) -------------------------------------------------------------------------
    print("2. Pair records")
    rec = [2, 2, 0, 0, 2.0, 1, 2, 4.0026, 0.0, 1, 6, 12.0, 7.16191, E2P, 5.43, 0.2, 8, 16, 1.0, E2P]  # dev's 20
    p = Pair.from_record(rec)
    check("a 20-field record: light_excitation, heavy_excitation, threshold",
          p.light_excitation == E2P and p.heavy_excitation == 0.0 and abs(p.threshold - (7.16191 + E2P)) < 1e-12,
          (p.light_excitation, p.heavy_excitation, p.threshold))
    check("label marks the excited light particle", p.label == "2+*(4.439 MeV) + 0+", p.label)
    check("repr shows light_excitation", "light_excitation=4.43891" in repr(p), repr(p))
    old = Pair.from_record(rec[:19])
    check("a 19-field record (older engine) reads, light_excitation = 0, heavy = total",
          old.light_excitation == 0.0 and old.heavy_excitation == E2P, (old.light_excitation, old.heavy_excitation))
    check("Pair._NFIELDS is 21 (thm: binding_energy last)", Pair._NFIELDS == 21)
    thm_rec = Pair.from_record(rec + [2.2246])
    check("a 21-field record: light_excitation and binding_energy",
          thm_rec.light_excitation == E2P and thm_rec.binding_energy == 2.2246)

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
