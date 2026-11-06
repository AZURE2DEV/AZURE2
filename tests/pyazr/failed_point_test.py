#!/usr/bin/env python3
"""A point whose calculation fails is NaN, also in a composite (SUM / RATIO) segment, and the run says so.

ESegment::CalculateTheoreticalCrossSection evaluates a segment's base point and the points of its components.  Before
2026-11-01 it caught every exception there: a failed base point gave 0, a failed SUM component was left out of the sum
(as was a component whose value was NaN or inf), a failed RATIO denominator gave 0, all without a word; the parallel
pre-pass (EData::CalculatePointsInParallel) also set a failed point to 0.  Seen with tests/natural_neon on the R-matrix
route (Brune), where a singular matrix made 84 points of the natural SUM segments 0.  Now such a point is NaN and a
warning naming the segment, energy and error is printed once.

The failure is injected with the engine's test hook AZURE_TEST_FAIL_E=<lab energy> (AZURE_TEST_FAIL_EXIT=<exit pair
key> to fail one term only; AZURE_NO_COMPONENT_MAPPING=1 so that the component is a point of its own).  Models: tests/sum_lab_frame/models/sum_data.azr (2->3 plus 2->5) and
tests/ratio_frame/models/ratio_data.azr (2->3 over 2->2 at 60 deg), with data files written here.  For the base term
and for the component/denominator alone: the injected point is NaN, every other point is unchanged, and the warning is
printed.  Each case runs in its own process (the hook is read once).

Run from anywhere:  python3 tests/pyazr/failed_point_test.py
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
ENERGIES = ["5.95", "6.05", "6.15", "6.25", "6.35"]
FAIL = ENERGIES[2]
# suite, model, data file, base exit key, the other term's exit key and its role in the warning, the segment's
# energy shift (MeV; the hook compares the shifted lab energy)
MODELS = [("sum_lab_frame", "sum_data", "sum.dat", 3, 5, "SUM component", 0.05),
          ("ratio_frame", "ratio_data", "ratio.dat", 3, 2, "RATIO denominator", 0.0)]


def worker(work, model, out):
    import numpy as np
    from pyazr import azure2
    with azure2(os.path.join(work, model + ".azr"), cwd=work) as a:
        v = np.asarray(a.calculate_rwa(np.array(a.params_rwa, float))[0], float)
    json.dump([None if not np.isfinite(x) else float(x) for x in v], open(out, "w"))


def main():
    try:
        import numpy as np  # noqa: F401
        from pyazr import azure2  # noqa: F401
    except Exception as err:  # same skip as the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 0
    failures = []

    def check(name, ok, detail=""):
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
        if not ok:
            failures.append(name)

    for suite, model, datafile, base_exit, other_exit, role, shift in MODELS:
        work = tempfile.mkdtemp(prefix="pyazr_failed_point_")
        try:
            src = os.path.join(ROOT, "tests", suite, "models")
            shutil.copy(os.path.join(src, model + ".azr"), work)
            if os.path.isdir(os.path.join(src, "data")):
                shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
            else:
                os.makedirs(os.path.join(work, "data"))
            with open(os.path.join(work, "data", datafile), "w") as f:
                for e in ENERGIES:
                    f.write(f"{e}  90  1.0  0.01\n")
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))

            def run(env_extra):
                out = os.path.join(work, "res.json")
                if os.path.exists(out):
                    os.remove(out)
                env = dict(os.environ, OMP_NUM_THREADS="2", **env_extra)
                p = subprocess.run([sys.executable, os.path.abspath(__file__), "worker", work, model, out],
                                   env=env, capture_output=True, text=True, timeout=600)
                vals = json.load(open(out)) if p.returncode == 0 and os.path.exists(out) else None
                return vals, p.stdout + p.stderr

            base, log0 = run({})
            if base is None or len(base) != len(ENERGIES) or any(x is None for x in base):
                check(f"{model}: reference run", False, log0[-800:])
                continue
            k = ENERGIES.index(FAIL)
            # The other term alone: with AZURE_NO_COMPONENT_MAPPING=1 it is a point of its own (mapped onto the base
            # point, it is solved inside the base point's calculation, whose failure is the first case).
            for label, exit_key, what, extra in (("base term", base_exit, "base term", {}),
                                                 (role, other_exit, role, {"AZURE_NO_COMPONENT_MAPPING": "1"})):
                vals, log = run(dict({"AZURE_TEST_FAIL_E": repr(float(FAIL) + shift),
                                      "AZURE_TEST_FAIL_EXIT": str(exit_key)}, **extra))
                if vals is None:
                    check(f"{model}, failed {label}: ran", False, log[-800:])
                    continue
                check(f"{model}, failed {label}: that point is NaN", vals[k] is None,
                      f"value {vals[k]!r} (reference {base[k]!r})")
                others = [i for i in range(len(base)) if i != k]
                same = all(vals[i] is not None and abs(vals[i] - base[i]) <= 1e-12 * abs(base[i]) for i in others)
                check(f"{model}, failed {label}: the other points unchanged", same,
                      f"{[vals[i] for i in others]} vs {[base[i] for i in others]}")
                # Reported by whichever path meets it first (the parallel pre-pass: "point").
                check(f"{model}, failed {label}: warning printed", "WARNING: the calculation failed at segment" in log
                      and "injected failure" in log, "no warning in the output")
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
        worker(*sys.argv[2:5])
        sys.exit(0)
    sys.exit(main())
