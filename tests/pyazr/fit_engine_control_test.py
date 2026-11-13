#!/usr/bin/env python3
"""The controls of azure2.fit_with_engine: the fitter table, progress records and the stop request (fit.md).

(a) pyazr.kinds()["fitters"] is the engine's table (_azure2.fitters()): the program's switches by name; an unknown
    method is refused with the list of fitters; MINOS after a least-squares row is refused.
(b) Progress on T2 (MIGRAD, call limit 40): records arrive (count > 0), their best objective never increases, the calls
    count up, every record's best is the lowest objective seen so far; the fit ends with status call_limit.
(c) A stop from the progress callable at the 15th record (MIGRAD and LM): status "stopped", no record after it (one
    thread: the evaluation running when the request came is the last one computed), the end point is the best point
    seen and its objective the best record's; evaluating the start afterwards gives the numbers of before, bit for bit.
(d) A stop from a second thread (a StopRequest, and a threading.Event) ends the call: status "stopped", at most one
    more record after the request (the evaluation that was running), well before the fit would end on its own.
(e) An exception in the progress callable stops the fit and is raised; the session evaluates as before.

Run: python3 tests/pyazr/fit_engine_control_test.py
"""
import os
import shutil
import sys
import tempfile
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
TESTS = os.path.join(ROOT, "tests")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    import numpy as np
    import pyazr
    from pyazr import _azure2
    from pyazr.session_fit import StopRequest
    pyazr.kinds()
except Exception as err:  # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def t2_copy(tmp, name):
    dst = os.path.join(tmp, name)
    shutil.copytree(os.path.join(TESTS, "upgrade", "T2"), dst, ignore=shutil.ignore_patterns("output", "checks"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(dst, d), exist_ok=True)
    return os.path.join(dst, "T2.azr")


def main():
    # (a) the table
    rows = _azure2.fitters()
    kinds = pyazr.kinds()["fitters"]
    check("kinds()['fitters'] lists the engine's rows in order", list(kinds) == [r["name"] for r in rows],
          f"{list(kinds)}")
    options = {n: k["option"] for n, k in kinds.items()}
    check("the program's switches by fitter", options.get("migrad") == "" and options.get("lm") == "--use-lm" and
          options.get("gsl-lm") == "--use-gsl-lm" and options.get("migrad-gradient") == "--use-gradient", f"{options}")
    with tempfile.TemporaryDirectory() as tmp, pyazr.azure2(t2_copy(tmp, "a")) as s:
        try:
            s.fit_with_engine("simplex")
            check("an unknown method is refused", False, "no error")
        except ValueError as e:
            check("an unknown method is refused with the list", "migrad" in str(e) and "lm" in str(e), str(e))
        try:
            s.fit_with_engine("lm", minos=True)
            check("MINOS after lm is refused", False, "no error")
        except ValueError as e:
            check("MINOS after lm is refused", "MINOS" in str(e), str(e))

        x0 = np.array(s.params_rwa, float)
        before = s.evaluate(x0)

        # (b) progress records
        everything = []
        r = s.fit_with_engine("migrad", max_calls=40, progress=everything.append)
        records = [rec for rec in everything if rec["stage"] == "fit"]
        check("progress: the fit's records, then HESSE's", {rec["stage"] for rec in everything} == {"fit", "hesse"}
              and [rec["stage"] for rec in everything] == sorted((rec["stage"] for rec in everything),
                                                                 key=["fit", "hesse"].index))
        bests = [rec["best"] for rec in records]
        check(f"progress: {len(records)} records, status {r.status}", len(records) > 0 and r.status == "call_limit")
        check("progress: the best objective never increases", all(b2 <= b1 for b1, b2 in zip(bests, bests[1:])))
        check("progress: each record's best is the lowest objective so far",
              all(rec["best"] == min(q["objective"] for q in records[:k + 1]) for k, rec in enumerate(records)))
        check("progress: the calls count up", [rec["calls"] for rec in records] == list(range(1, len(records) + 1)))
        check("progress: records carry the free vector", records[0]["x"] is not None and len(records[0]["x"]) == len(x0))

        # (c) a stop at the 15th record, MIGRAD and LM
        for method in ("migrad", "lm"):
            seen, after = [], []
            stop = StopRequest()

            def progress(rec, seen=seen, after=after, stop=stop):
                (after if stop.is_set() else seen).append(rec)
                if len(seen) == (15 if method == "migrad" else 2) and not stop.is_set():
                    stop.set()
            r = s.fit_with_engine(method, progress=progress, stop=stop)
            best = min(rec["objective"] for rec in seen)
            check(f"{method} stopped from its progress: status {r.status}, {len(seen)} records, {len(after)} after",
                  r.status == "stopped" and not r.converged and len(after) == 0)
            check(f"{method} stopped: the end point is the best point seen ({r.objective:.12g} = {best:.12g})",
                  abs(r.objective - best) <= 1e-12 * abs(best))
            check(f"{method} stopped: no covariance", r.covariance is None)
            again = s.evaluate(x0)
            check(f"after the {method} stop the start evaluates as before (bit for bit)",
                  again.objective == before.objective and again.chi2 == before.chi2,
                  f"{again.objective!r} against {before.objective!r}")

        # (d) a stop from a second thread
        for kind in ("StopRequest", "threading.Event"):
            stop = StopRequest() if kind == "StopRequest" else threading.Event()
            started, records, at_stop = threading.Event(), [], []

            def progress(rec, records=records, started=started):
                records.append(rec)
                if len(records) == 5:
                    started.set()

            def stopper(stop=stop, started=started, records=records, at_stop=at_stop):
                started.wait(600)
                at_stop.append(len(records))
                stop.set()
            t = threading.Thread(target=stopper)
            t.start()
            r = s.fit_with_engine("migrad", progress=progress, stop=stop)
            t.join()
            more = len(records) - at_stop[0] if at_stop else None
            check(f"a {kind} set from another thread stops the fit: status {r.status}, {more} record(s) after it",
                  r.status == "stopped" and more is not None and more <= 1)

        # (e) an exception in the progress callable
        def broken(rec):
            if rec["calls"] == 3:
                raise RuntimeError("progress failed")
        try:
            s.fit_with_engine("migrad", progress=broken)
            check("an exception in progress is raised", False, "no error")
        except RuntimeError as e:
            check("an exception in progress stops the fit and is raised", "progress failed" in str(e), str(e))
        again = s.evaluate(x0)
        check("after the exception the start evaluates as before", again.objective == before.objective)
    print("FAIL" if failures else "OK", f"({len(failures)} failure(s))" if failures else "")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
