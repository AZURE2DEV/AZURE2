#!/usr/bin/env python3
"""Channel order inside a level is presentation, not physics.

A J-group's channels can be written in any order -- the group owns the set, and
every level of the group carries a width for each of them. Listing them
differently must not change a result, and must certainly not crash.

It did both. CNuc::TransformIn and CNuc::TransformOut each padded their
per-channel arrays for a non-particle channel by copying entry [0], to keep the
array indexed by channel number. When channel 1 was *itself* a photon or beta
channel that array was still empty, so the copy read off the end of an empty
vector: a segfault on a file whose <levels> simply happens to list its capture
channel first. Every reader of those entries is guarded by RadType == 'P', so
the padding value is never used and zero does as well.

TransformIn is on the default path. TransformOut's copy is reached only with
the Brune formalism off, which is why neither had been seen.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/channel_order_test.py
"""
import importlib.util
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2
    from pyazr.azrfile import AzrModel
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)

# 13N is the fixture with a photon channel to move: reversing each level's
# channels puts it first, which is the condition that used to crash.
PROJECT = os.path.join(ROOT, "tests", "13N")

with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "13N")
    shutil.copytree(PROJECT, work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))

    plain = os.path.join(work, "13N.azr")
    model = AzrModel.from_file(plain)
    photon_first = 0
    for lv in model.levels:
        lv.channels = list(reversed(lv.channels))
        if lv.channels[0].is_photon:
            photon_first += 1
    reversed_path = os.path.join(work, "reversed.azr")
    model.write(reversed_path)
    check("the reversed file leads with a photon channel", photon_first > 0,
          "nothing to test: no level starts with one")

    # Every variant gets its own output directory. AZURE2 caches its
    # external-capture integrals in output/intEC.dat and reuses them when the
    # file is there, so runs sharing one directory are not independent -- on
    # this model that alone moves chi-squared by 0.6, which would swamp what
    # this test is trying to see.
    results = {}
    for name, source in (("as written", plain), ("reversed", reversed_path)):
        for brune in (True, False):
            tag = f"{name.replace(' ', '_')}_{brune}"
            outdir = os.path.join(work, "out_" + tag)
            os.makedirs(outdir, exist_ok=True)
            variant = AzrModel.from_file(source)
            variant.set_output_dir(outdir)
            path = os.path.join(work, tag + ".azr")
            variant.write(path)
            try:
                with azure2(path, cwd=work, use_brune=brune) as m:
                    x = np.asarray(m.params_rwa, float)
                    results[(name, brune)] = float(np.sum(m.calculate_chi2_rwa(x)))
            except Exception as err:
                results[(name, brune)] = f"{type(err).__name__}: {err}"

    for key, value in sorted(results.items(), key=lambda kv: str(kv[0])):
        print(f"        {key[0]:11} brune={str(key[1]):5} -> {value}")

    for brune in (True, False):
        a, b = results[("as written", brune)], results[("reversed", brune)]
        check(f"reversed channels load (brune={brune})", isinstance(b, float), str(b))
        if isinstance(a, float) and isinstance(b, float):
            check(f"and give the same chi-squared (brune={brune})",
                  abs(a - b) <= 1e-9 * max(abs(a), 1.0), f"{a} vs {b}")

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all channel-order checks passed")
