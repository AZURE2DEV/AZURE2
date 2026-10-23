#!/usr/bin/env python3
"""pyazr and relativistic kinematics (relativistic=1 in <config>).

AzrModel.relativistic reads and writes the key; an engine session on tests/relativistic/models/na_hs_rel.azr
(n + 4He s-wave hard sphere) computes sigma = 4 pi sin^2(k a)/k^2 with k from the invariant mass at the c.m.
energies it reports, and the non-relativistic model with k^2 = 2 mu T.

Run from anywhere:  python3 tests/pyazr/relativistic_test.py
"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    import pyazr  # noqa: E402
except ImportError as e:  # as the other pyazr tests: a missing package dependency (e.g. mpmath) -> nothing to test
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(0)

U, HBARC = 931.4940880, 197.32696310
MN, MA = 1.00866491595, 4.00260325413


def k2(T, rel):
    a, b = MN * U, MA * U
    m = a + b
    if not rel:
        return 2 * a * b / m * T / HBARC ** 2
    return T * (T + 2 * m) * (T + 2 * a) * (T + 2 * b) / (4 * (m + T) ** 2) / HBARC ** 2


def main():
    models = os.path.join(ROOT, "tests", "relativistic", "models")
    work = tempfile.mkdtemp(prefix="pyazr_rel_")
    fails = 0
    try:
        m = pyazr.AzrModel.from_file(os.path.join(models, "na_hs.azr"))
        assert not m.relativistic
        m.relativistic = True
        path = os.path.join(work, "flip.azr")
        m.write(path)
        assert pyazr.AzrModel.from_file(path).relativistic, "relativistic=1 not written"
        m.relativistic = False
        m.write(path)
        assert "relativistic" not in open(path).read(), "key not removed"
        print("  ok    AzrModel.relativistic round trip")
        try:
            from pyazr import _azure2  # noqa: F401  the compiled engine
        except Exception as err:  # engine not built
            print(f"  skip  engine part: {type(err).__name__}: {err}")
            return
        for rel in (False, True):
            name = "na_hs" + ("_rel" if rel else "")
            d = os.path.join(work, name)
            os.makedirs(os.path.join(d, "output"))
            os.makedirs(os.path.join(d, "checks"))
            shutil.copy(os.path.join(models, name + ".azr"), os.path.join(d, "t.azr"))
            # the model has extrapolation segments only, so the session starts in extrapolation mode
            with pyazr.azure2(os.path.join(d, "t.azr"), cwd=d, data_mode=False, use_brune=False) as s:
                sig = np.asarray(s.calculate(s.params)[0], float)
                T = np.asarray(s.calculate_energies(s.params)[0], float)
            k = np.sqrt(k2(T, rel))
            ref = 4 * np.pi * np.sin(3.0 * k) ** 2 / k ** 2 / 100.0
            r = float(np.max(np.abs(sig - ref) / ref))
            ok = r < 1e-8
            fails += 0 if ok else 1
            print(f"  {'ok  ' if ok else 'FAIL'}  {name}: sigma = 4 pi sin^2(ka)/k^2 (k {'invariant' if rel else 'NR'}) to {r:.1e}")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print(f"{fails} failure(s)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
