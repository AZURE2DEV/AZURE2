#!/usr/bin/env python3
"""The API calculates the points of an evaluation in one parallel pass, with the serial path's results.

Up to 2026-10-30 the CLI's chi-squared calculated every point in one OpenMP pass (AZURECalc::PrecalculatePoints)
but the API calls behind pyazr (calculate_chi2_rwa, calculate_chi2, calculate_rwa, calculate, calculate_all_rwa)
walked the points one at a time on one thread: a 17O+a evaluation ran at ~1.1 cpu/wall with OMP_NUM_THREADS=4
(2,490 s on 2,719 points; R-matrix/AZURE2_code_audit/10-25-26/readme).  They now share the CLI's pass
(EData::CalculatePointsInParallel); after any point failure they repeat the calculation point by point, so a failed
point keeps the serial handling.

Checks, each model in fresh processes with OMP_NUM_THREADS=4:
(a) with AZURE_PROFILE=1 every one of those calls reports a "points" stage (the parallel pass ran).  The old
    API never times that stage, so this fails on a build before the change;
(b) every result -- the chi2 of both bases, every calculated segment array -- equals the serial path's
    (AZURE_SERIAL_CHI2=1) to 1e-12 relative, at the file's parameters and at moved ones.  Not bit for bit: a
    point's last bits depend on what its thread calculated before it (thread-reused workspaces), so the parallel
    pass, the serial path and the same pass on one thread (only the point order differs) disagree by an ulp or
    two at a few natural_neon points (2026-10-30: up to ~2e-15 relative, run to run); tests/13N agreed bit for bit.
    Each mode starts from its own fresh copy of the model: a session that reads the external-capture integrals
    back from output/intEC.dat (written with 6 significant digits) differs from one that computed them by
    ~1e-6 relative (13N chi2 by 0.62), which is not what this test is about.

Models: tests/13N (capture, angular distributions, mapped points) and tests/natural_neon (SUM components,
natural-sample segments).

Run from anywhere:  python3 tests/pyazr/parallel_points_test.py
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

MODELS = [("13N", "13N.azr"), ("natural_neon", "natural_neon.azr")]
CALLS = ["calculate_chi2_rwa", "calculate_chi2", "calculate_rwa", "calculate", "calculate_all_rwa"]


def worker(model_dir, azr, out):
    """One process: every call at x0 and at a moved x; results and the profile labels to `out`."""
    import numpy as np
    from pyazr import azure2
    a = azure2(os.path.join(model_dir, azr), cwd=model_dir)
    res = {}
    x0 = np.array(a.params_rwa, float)
    p0 = np.array(a.params, float)
    full0 = np.array(a.sess.params_all_rwa(), float)
    for tag, s in [("x0", 1.0), ("moved", 1.03)]:
        x, p, full = x0 * s, p0 * s, full0 * s
        for call in CALLS:
            sys.stderr.write(f"@@ {tag} {call}\n")
            sys.stderr.flush()
            arg = {"calculate_chi2_rwa": x, "calculate_rwa": x, "calculate_chi2": p, "calculate": p,
                   "calculate_all_rwa": full}[call]
            r = getattr(a, call)(arg)
            res[f"{tag} {call}"] = [np.asarray(c, float).ravel().tolist() for c in r]
    json.dump(res, open(out, "w"))


def run(model_dir, azr, env_extra):
    out = tempfile.mktemp(suffix=".json", dir=model_dir)
    env = {k: v for k, v in os.environ.items() if k not in ("AZURE_SERIAL_CHI2", "AZURE_PROFILE")}
    env.update(OMP_NUM_THREADS="4", **env_extra)
    p = subprocess.run([sys.executable, os.path.abspath(__file__), "worker", model_dir, azr, out], env=env,
                       capture_output=True, text=True, timeout=1800)
    if p.returncode != 0:
        raise RuntimeError(f"worker failed ({p.returncode}): {p.stderr[-2000:]}")
    return json.load(open(out)), p.stderr


def close(a, b, rel=1e-12):
    """Equal lists of arrays to `rel` relative (NaN equal to NaN); returns (ok, largest relative difference)."""
    import numpy as np
    if len(a) != len(b):
        return False, float("inf")
    worst = 0.0
    for u, v in zip(a, b):
        u, v = np.asarray(u, float), np.asarray(v, float)
        if u.shape != v.shape or not np.array_equal(np.isnan(u), np.isnan(v)):
            return False, float("inf")
        ok = ~np.isnan(u)
        d = np.abs(u[ok] - v[ok]) / np.maximum(np.abs(v[ok]), 1e-300)
        if d.size:
            worst = max(worst, float(d.max()))
    return worst <= rel, worst


def main():
    try:
        import numpy  # noqa: F401
        from pyazr import azure2  # noqa: F401
    except Exception as err:  # same skip as the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 0
    failures = []

    def check(name, ok, detail=""):
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
        if not ok:
            failures.append(name)

    for sub, azr in MODELS:
        print(f"\n{sub}")
        work = tempfile.mkdtemp(prefix="pyazr_parallel_points_")
        try:
            dirs = {}
            for mode in ("par", "ser"):        # a fresh copy each: no intEC.dat from the other run
                d = dirs[mode] = os.path.join(work, mode, sub)
                shutil.copytree(os.path.join(ROOT, "tests", sub), d)
                os.makedirs(os.path.join(d, "output"), exist_ok=True)
                os.makedirs(os.path.join(d, "checks"), exist_ok=True)
            par, log = run(dirs["par"], azr, {"AZURE_PROFILE": "1"})
            ser, _ = run(dirs["ser"], azr, {"AZURE_SERIAL_CHI2": "1"})
            # (a) the parallel pass ran in every call: the "points" stage between this call's marker and the next
            blocks = log.split("@@ ")[1:]
            for b in blocks:
                head = b.split("\n", 1)[0]
                check(f"(a) {head}: parallel point pass", "[profile]   points" in b,
                      "no 'points' stage in the profile report")
            check("(a) every call reported", len(blocks) == 2 * len(CALLS), f"{len(blocks)} markers")
            # (b) identical to the serial path
            for k in sorted(par):
                ok, worst = close(par[k], ser.get(k, []))
                check(f"(b) {k}: serial path to 1e-12 (largest {worst:.1e})", ok)
        finally:
            shutil.rmtree(work, ignore_errors=True)
    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s): {failures}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(*sys.argv[2:5])
        sys.exit(0)
    sys.exit(main())
