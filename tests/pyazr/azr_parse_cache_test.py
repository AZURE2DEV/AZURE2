#!/usr/bin/env python3
"""A session parses its .azr once, not on every evaluation.

Before 2026-11-03 every azure2.calculate_rwa / calculate / calculate_all_rwa call re-read and
re-parsed the whole .azr once PER SEGMENT (_segment_values looked up the segment's observable
through active_datasets -> datasets -> a Python reader of the file), and every penalties()/objective()
call once more.  On the 3H+d model from /groups (NFS) that was 83 ms per calculate_rwa against
8 ms from local disk, and it is what made forked pyazr workers there slow (fork_bench, 2026-11-02).

Checks, on tests/channel_space/channel_space.azr plus an E1 copy of its capture data set (a record with
`multipole E1`, so the E1-component path of calculate_rwa is exercised):
  1. after the first evaluation, calculate_rwa / calculate / calculate_all_rwa / penalties /
     objective / segment_chi2 open the .azr ZERO times (the old code: once per segment per call);
  2. the values are unchanged: each segment is the engine's total, the E1 segment the engine's
     E1 component, and objective = chi2 + penalties as before;
  3. after a change to the file on disk the penalties stay the engine's, on the model it loaded (since pyazr
     2.11), and configure() re-reads the file once, the next call not at all;
  4. the public datasets / active_datasets / extrapolations are copies and cannot be edited (the engine's
     description, frozen): the session's penalties do not move.
The old code fails 1.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/azr_parse_cache_test.py
"""
import builtins
import os
import shutil
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


class OpenCounter:
    """Counts Python-level opens of one file (the engine's C++ reads are not counted)."""

    def __init__(self, path):
        self.path = os.path.realpath(path)
        self.n = 0
        self._real = builtins.open

    def __enter__(self):
        real = self._real

        def counting(f, *a, **k):
            try:
                if os.path.realpath(os.fspath(f)) == self.path:
                    self.n += 1
            except TypeError:
                pass
            return real(f, *a, **k)
        builtins.open = counting
        return self

    def __exit__(self, *exc):
        builtins.open = self._real


def main():
    sys.path.insert(0, ROOT)
    try:
        import numpy as np
        from pyazr import azure2, Project
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"SKIP: engine not importable ({e})")
        return 0

    work = tempfile.mkdtemp(prefix="pyazr_parse_cache_")
    try:
        src = os.path.join(ROOT, "tests", "channel_space")
        shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
        for d in ("output", "checks"):
            os.makedirs(os.path.join(work, d), exist_ok=True)
        azr = os.path.join(work, "cache.azr")
        p = Project.open(os.path.join(src, "channel_space.azr"))
        cap = p["capture"]
        p.add("dataset", "capture_e1", **{k: v for k, v in cap.fields.items()}, multipole="E1")
        p.save(azr)
        text = open(azr).read()

        with azure2(azr, cwd=work) as m:
            x = np.asarray(m.params_rwa, float)
            n = m.nsegments
            names = [d.name for d in m.active_datasets]
            check("the model has the E1 data set", "capture_e1" in names, f"{names}")
            k_e1 = names.index("capture_e1") if "capture_e1" in names else n - 1

            m.calculate_rwa(x)
            m.objective(x)                      # warm whatever is cached

            print("1. evaluations do not re-read the .azr")
            calls = {
                "calculate_rwa": lambda: m.calculate_rwa(x),
                "calculate": lambda: m.calculate(np.asarray(m.params, float)),
                "calculate_all_rwa": lambda: m.calculate_all_rwa(np.asarray(m.sess.params_all_rwa(), float)),
                "penalties": lambda: m.penalties(x),
                "objective": lambda: m.objective(x),
                "segment_chi2": lambda: m.segment_chi2(x),
            }
            for name, fn in calls.items():
                with OpenCounter(azr) as c:
                    for _ in range(5):
                        fn()
                check(f"{name}: 0 opens of the .azr in 5 calls", c.n == 0, f"{c.n} opens ({n} segments)")

            print("\n2. values unchanged")
            vals = m.calculate_rwa(x)
            ok = True
            for i in range(n):
                ref = (m.sess.calculated_segments_e1(i) if i == k_e1 else m.sess.calculated_segments(i))
                ok &= np.array_equal(np.asarray(vals[i]), np.asarray(ref))
            check("each segment = the engine's total, the E1 segment its E1 component", ok)
            tot = np.asarray(m.sess.calculated_segments(k_e1))
            e1 = np.asarray(m.sess.calculated_segments_e1(k_e1))
            print(f"     (E1 segment: max |total - E1| / total = {np.max(np.abs(tot - e1) / np.abs(tot)):.3g})")
            pen = m.penalties(x)
            obj = m.objective(x)
            ref_obj = (float(np.sum(m.calculate_chi2_rwa(x))) + float(np.sum(pen["norm"]))
                       + float(np.sum(pen["shift"])) + float(np.sum(pen["shift_sqrt"]))
                       + float(np.sum(pen["shift_fn"])) + pen["screening"] + pen["density"])
            # The objective is the engine's (docs/upgrade/objective.md step 4): its terms summed in the engine's
            # order, so the Python sum may differ in the last bit.
            check("objective = chi2 + penalties", abs(obj - ref_obj) <= 1e-13 * abs(ref_obj), f"{obj} vs {ref_obj}")

            print("\n3. an edit of the file on disk is still seen")
            # data set 1: nominal norm 1 -> 1.1 (5 % error, norm parameter still 1.0 in the engine)
            norm_before = float(m.penalties(x)["norm"][0])
            old = "    norm           1 fit prior 1 5%\n"
            new = "    norm           1.1 fit prior 1.1 5%\n"
            check("the edit target is in the file", old in text)
            time.sleep(0.01)
            with open(azr, "w") as f:
                f.write(text.replace(old, new, 1))
            # Since pyazr 2.11 the penalties are the engine's (objective.md step 4), computed on the model the
            # engine loaded -- the one its chi-squared uses; a session does not re-read the .azr (configure()
            # re-reads metadata from the engine).  Before, penalties() parsed the file in Python and followed an
            # edit at once, so they could describe another model than the chi-squared.  The Python parse cache
            # still follows the edit (its own views of the file).
            norm_unconfigured = float(m.penalties(x)["norm"][0])
            want = ((1.0 - 1.1) / (1.1 / 100.0 * 5.0)) ** 2
            check("norm penalty before the edit is 0", norm_before == 0.0, f"{norm_before}")
            check("the engine's penalty keeps the loaded model", norm_unconfigured == 0.0, f"{norm_unconfigured}")
            with OpenCounter(azr) as c:
                m.configure()
                m.active_datasets
                norm_after = float(m.penalties(x)["norm"][0])
            with OpenCounter(azr) as c2:
                m.active_datasets
                m.penalties(x)
            check("... also after configure()", norm_after == 0.0, f"{norm_after}")
            del want
            check("configure() re-reads, the next call does not", c.n >= 1 and c2.n == 0,
                  f"{c.n} then {c2.n} opens")
            with open(azr, "w") as f:           # put it back
                f.write(text)

            print("\n4. public views are copies")
            p0 = m.penalties(x)["norm"].copy()
            import dataclasses
            refused = 0
            for view in (m.active_datasets, m.datasets, m.extrapolations):
                if len(view):
                    try:
                        view[0].active = not view[0].active
                    except dataclasses.FrozenInstanceError:
                        refused += 1
            p1 = m.penalties(x)["norm"]
            check("the views cannot be edited", refused == 2 + (1 if len(m.extrapolations) else 0), str(refused))
            check("the penalties do not move", np.array_equal(p0, p1), f"{p0} vs {p1}")
            check("a fresh view is another object", m.active_datasets is not m.active_datasets)

            t0 = time.perf_counter()
            for _ in range(50):
                m.calculate_rwa(x)
            print(f"\n  ({(time.perf_counter() - t0) / 50 * 1e3:.2f} ms per calculate_rwa, {n} segments, {work})")
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
