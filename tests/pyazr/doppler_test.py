#!/usr/bin/env python3
"""Free-gas Doppler broadening through pyazr: AzrModel.add_target_effect(doppler=T) writes the ``doppler`` keyword
and refuses T <= 0 (pure Python); through the engine, a Doppler-only line changes sigma_t at the narrowed n + 12C
resonance of tests/doppler_fgm by tens of percent and a transmission line gives exp(-N sigma_D) (the API path
assembles transformed segments through ESegment::NeedsAssembly, which must include the Doppler step).

The engine part needs the compiled module and skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/doppler_test.py
"""
import importlib.util
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SUITE = os.path.join(ROOT, "tests", "doppler_fgm")
sys.path.insert(0, ROOT)

try:
    from pyazr.azrfile import AzrModel  # noqa: E402
except ImportError as e:
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(0)

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({detail})"))
    if not ok:
        failures.append(name)


tmp = tempfile.mkdtemp(prefix="pyazr_doppler_")
try:
    spec = importlib.util.spec_from_file_location("doppler_models", os.path.join(SUITE, "make_models.py"))
    mm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mm)
    mm.main(os.path.join(ROOT, "tests", "reaction_rate", "models", "n12C.azr"), tmp)

    # --- AzrModel (no engine) ---------------------------------------------------------------------
    b = AzrModel.from_file(os.path.join(tmp, "bare.azr"))
    b.add_target_effect("1,2", n_points=400, doppler=293.6)
    path = os.path.join(tmp, "added.azr")
    b.write(path)
    lines = open(path).read().split("<targetInt>")[1].split("</targetInt>")[0].strip().splitlines()
    check("add_target_effect(doppler=293.6) writes one line ending in 'doppler 293.6'",
          len(lines) == 1 and lines[0].split()[-2:] == ["doppler", "293.6"], lines)
    for bad in (0, -5):
        try:
            AzrModel.from_file(os.path.join(tmp, "bare.azr")).add_target_effect("1", doppler=bad)
            check(f"add_target_effect(doppler={bad}) refused", False)
        except ValueError:
            check(f"add_target_effect(doppler={bad}) refused", True)

    # --- engine ----------------------------------------------------------------------------------
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import numpy as np
        from pyazr import azure2
        from pyazr import _azure2  # noqa: F401
    except Exception as err:
        print(f"skip engine part: {type(err).__name__}: {err}")
    else:
        def evaluate(name):
            work = os.path.join(tmp, "run_" + name)
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))
            shutil.copy(os.path.join(tmp, name + ".azr"), work)
            os.symlink(os.path.join(tmp, "data"), os.path.join(work, "data"))
            with azure2(os.path.join(work, name + ".azr"), cwd=work) as s:
                return [np.asarray(v, float) for v in s.calculate_rwa(np.asarray(s.params_rwa, float))]
        bare, dop, dopT = evaluate("bare"), evaluate("dop"), evaluate("dopT")
        rel = float(np.max(np.abs(dop[0] / bare[0] - 1.0)))
        check("a Doppler-only line changes sigma_t by more than 5 % somewhere", rel > 0.05, rel)
        want = np.exp(-mm.N * dop[0])
        err = float(np.max(np.abs(dopT[0] - want) / want))
        check("transmission with doppler alone = exp(-N sigma_D) of the Doppler-only sigma_t (1e-9)", err < 1e-9, err)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: pyazr doppler")
