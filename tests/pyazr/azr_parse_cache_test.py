#!/usr/bin/env python3
"""A session parses its .azr once, not on every evaluation.

Before 2026-11-03 every azure2.calculate_rwa / calculate / calculate_all_rwa call re-read and
re-parsed the whole .azr once PER SEGMENT (_segment_values looked up the segment's observable
through active_datasets -> datasets -> SegmentSet.from_file), and every penalties()/objective()
call once more.  On the 3H+d model from /groups (NFS) that was 83 ms per calculate_rwa against
8 ms from local disk, and it is what made forked pyazr workers there slow (fork_bench, 2026-11-02).

Checks, on tests/channel_space/channel_space.azr plus an angle-integrated-E1 copy of its capture
segment (isDiff 5, so the E1-component path of calculate_rwa is exercised):
  1. after the first evaluation, calculate_rwa / calculate / calculate_all_rwa / penalties /
     objective / segment_chi2 open the .azr ZERO times (the old code: once per segment per call);
  2. the values are unchanged: each segment is the engine's total, the E1 segment the engine's
     E1 component, and objective = chi2 + penalties as before;
  3. a change to the file on disk is still seen (the penalty follows an edited nominal norm),
     and configure() re-reads;
  4. the public datasets / active_datasets / extrapolations are copies: editing them does not
     change the session's penalties.
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
        from pyazr import azure2
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
        text = open(os.path.join(src, "channel_space.azr")).read()
        e1_line = "1 1 3 0.5 7 0 180 5 1 0 0 0 0 0 data/capture.dat 0 0\n"
        text = text.replace("</segmentsData>\n", e1_line + "</segmentsData>\n", 1)
        with open(azr, "w") as f:
            f.write(text)

        with azure2(azr, cwd=work) as m:
            x = np.asarray(m.params_rwa, float)
            n = m.nsegments
            obs = [d.observable for d in m.active_datasets]
            check("the model has an angle-integrated-E1 segment", "angle-integrated-E1" in obs, f"{obs}")
            k_e1 = obs.index("angle-integrated-E1") if "angle-integrated-E1" in obs else n - 1

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
            check("objective = chi2 + penalties", obj == ref_obj, f"{obj} vs {ref_obj}")

            print("\n3. an edit of the file on disk is still seen")
            # segment 1: nominal norm 1 -> 1.1 (5 % error, norm parameter still 1.0 in the engine)
            norm_before = float(m.penalties(x)["norm"][0])
            old = "1 1 1 0.5 7 0 180 1 1 1 5 0 0 0 data/elastic_90.dat 0 0\n"
            new = "1 1 1 0.5 7 0 180 1 1.1 1 5 0 0 0 data/elastic_90.dat 0 0\n"
            check("the edit target is in the file", old in text)
            time.sleep(0.01)
            with open(azr, "w") as f:
                f.write(text.replace(old, new, 1))
            norm_after = float(m.penalties(x)["norm"][0])
            want = ((1.0 - 1.1) / (1.1 / 100.0 * 5.0)) ** 2
            check("norm penalty before the edit is 0", norm_before == 0.0, f"{norm_before}")
            check("norm penalty follows the edited nominal", abs(norm_after - want) <= 1e-12 * want,
                  f"{norm_after} vs {want}")
            with OpenCounter(azr) as c:
                m.configure()
                m.penalties(x)
            with OpenCounter(azr) as c2:
                m.penalties(x)
            # two parses after configure(): the segment list and the AzrModel (scattering lengths)
            check("configure() re-reads, the next call does not", c.n >= 1 and c2.n == 0,
                  f"{c.n} then {c2.n} opens")
            with open(azr, "w") as f:           # put it back
                f.write(text)

            print("\n4. public views are copies")
            p0 = m.penalties(x)["norm"].copy()
            ds = m.active_datasets
            ds[0].norm = 5.0
            alld = m.datasets
            alld[0].norm = 5.0
            ex = m.extrapolations
            if len(ex):
                ex[0].active = not ex[0].active
            p1 = m.penalties(x)["norm"]
            check("editing active_datasets/datasets does not move the penalties", np.array_equal(p0, p1),
                  f"{p0} vs {p1}")
            check("a fresh active_datasets is unedited", m.active_datasets[0].norm == 1.0,
                  f"{m.active_datasets[0].norm}")

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
