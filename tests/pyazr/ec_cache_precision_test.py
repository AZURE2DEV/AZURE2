#!/usr/bin/env python3
"""A session that reads the external-capture integrals back from intEC.dat equals the session that computed them.

Up to 2026-10-30 output/intEC.dat was written with the stream's default 6 significant digits, e.g.
"(6.35126e-10,1.50813e-05)", so every run after the first in a directory computed its capture cross sections from
rounded amplitudes: tests/13N chi2 336479.7107 when the integrals were computed, 336479.0869 when they were read
back (-0.62, ~2e-6 relative), and up to ~1e-5 relative on single capture points.  The file now holds 17
significant digits (signature v5, so older files are recomputed once).

Checks on tests/13N (capture), each in a fresh process: (a) the first session computes the integrals and writes
output/intEC.dat with its signature; (b) a second session in the same directory reads them back (no new file is
written) and its chi2 and every calculated segment equal the first's to 1e-12 relative.

Run from anywhere:  python3 tests/pyazr/ec_cache_precision_test.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)


def worker(model_dir, out):
    import numpy as np
    from pyazr import azure2
    a = azure2(os.path.join(model_dir, "13N.azr"), cwd=model_dir)
    x = np.array(a.params_rwa, float)
    res = {"chi2": [float(np.sum(a.calculate_chi2_rwa(x)))],
           "segments": [np.asarray(c, float).ravel().tolist() for c in a.calculate_rwa(x)]}
    json.dump(res, open(out, "w"))


def main():
    try:
        import numpy as np
        from pyazr import azure2  # noqa: F401
    except Exception as err:  # same skip as the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 0
    failures = []

    def check(name, ok, detail=""):
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
        if not ok:
            failures.append(name)

    work = tempfile.mkdtemp(prefix="pyazr_ec_cache_")
    try:
        d = os.path.join(work, "13N")
        shutil.copytree(os.path.join(ROOT, "tests", "13N"), d)
        for sub in ("output", "checks"):
            shutil.rmtree(os.path.join(d, sub), ignore_errors=True)
            os.makedirs(os.path.join(d, sub))
        ec = os.path.join(d, "output", "intEC.dat")
        res = []
        for k in range(2):
            out = os.path.join(work, f"run{k}.json")
            p = subprocess.run([sys.executable, os.path.abspath(__file__), "worker", d, out], capture_output=True,
                               text=True, timeout=900)
            if p.returncode != 0:
                print(p.stderr[-2000:])
                check(f"session {k + 1} ran", False, f"exit {p.returncode}")
                return 1
            res.append(json.load(open(out)))
            if k == 0:
                check("(a) the first session writes intEC.dat and its signature",
                      os.path.exists(ec) and os.path.exists(ec + ".sig"))
                stamp = os.stat(ec).st_mtime_ns
        check("(b) the second session reads intEC.dat back (not rewritten)", os.stat(ec).st_mtime_ns == stamp)
        c0, c1 = res[0]["chi2"][0], res[1]["chi2"][0]
        rel = abs(c1 - c0) / abs(c0)
        check(f"(b) chi2 read back = computed to 1e-12 ({c0:.10f} vs {c1:.10f}, rel {rel:.1e})", rel <= 1e-12)
        worst = 0.0
        for u, v in zip(res[0]["segments"], res[1]["segments"]):
            u, v = np.asarray(u, float), np.asarray(v, float)
            ok = ~np.isnan(u) & ~np.isnan(v)
            if ok.any():
                worst = max(worst, float(np.max(np.abs(u[ok] - v[ok]) / np.maximum(np.abs(u[ok]), 1e-300))))
        check(f"(b) every segment read back = computed to 1e-12 (largest {worst:.1e})", worst <= 1e-12)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(sys.argv[2], sys.argv[3])
        sys.exit(0)
    sys.exit(main())
