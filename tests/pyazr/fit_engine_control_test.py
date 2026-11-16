#!/usr/bin/env python3
"""The controls of Session.fit_with_engine: the fitter table, progress records and the stop request (fit.md).

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
(f) What a fit returns, the same fields for every fitter: the table's columns call_limit and fallback; a MIGRAD that
    stops at its call limit on T2 (Minuit2 returns its start there) ends at the lowest objective of its records, with
    outcome call_limit and a notice; LM on T2 runs MIGRAD (no analytic Jacobian), which honours the limit, and says so
    in its notices; LM on T3 ignores max_calls and says so in its notices.
(g) MINOS on tests/identical_pp_res: every free parameter has an interval and the validity of each side
    (minos_valid: valid, at_limit, call_limit, new_minimum); MINOS's progress records are marked error_analysis.

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
    columns = {n: (k["call_limit"], k["fallback"]) for n, k in kinds.items()}
    check("the columns call_limit and fallback: MIGRAD honours a limit, LM and GSL-LM ignore it, fall back to migrad",
          columns.get("migrad") == (True, "") and columns.get("migrad-gradient") == (True, "") and
          columns.get("lm") == (False, "migrad") and columns.get("gsl-lm") == (False, "migrad"), f"{columns}")
    options = {n: k["option"] for n, k in kinds.items()}
    check("the program's switches by fitter", options.get("migrad") == "" and options.get("lm") == "--use-lm" and
          options.get("gsl-lm") == "--use-gsl-lm" and options.get("migrad-gradient") == "--use-gradient", f"{options}")
    with tempfile.TemporaryDirectory() as tmp, pyazr.Session(t2_copy(tmp, "a")) as s:
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

        x0 = s.parameters.x0
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
        lowest = min(rec["objective"] for rec in records)
        check(f"call limit: the result is the lowest objective of the records ({r.objective:.10g} = {lowest:.10g}), "
              "below the start", r.objective == lowest and lowest < records[0]["objective"])
        check(f"call limit: outcome {r.outcome!r}, a notice of the best point kept",
              r.outcome == "call_limit" and any("lowest objective it evaluated" in n for n in r.notices),
              f"{r.notices}")
        r = s.fit_with_engine("lm", max_calls=40)
        check(f"lm on T2 runs MIGRAD, which honours the limit: {r.status}, notices say the fallback",
              r.fell_back_to_migrad and r.status == "call_limit" and "falling back to MIGRAD" in r.notices[0] and
              not any("does not honour" in n for n in r.notices), f"{r.notices}")

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
                stop.set()  # first the request, then the count: a record between the two is not one after it
                at_stop.append(len(records))
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

    # (f) LM where it runs (T3) ignores max_calls and says so; (g) MINOS's validity per side
    with tempfile.TemporaryDirectory() as tmp:
        t3 = os.path.join(tmp, "T3")
        shutil.copytree(os.path.join(TESTS, "upgrade", "T3"), t3, ignore=shutil.ignore_patterns("output", "checks"))
        for d in ("output", "checks"):
            os.makedirs(os.path.join(t3, d), exist_ok=True)
        with pyazr.Session(os.path.join(t3, "T3.azr")) as s:
            r = s.fit_with_engine("lm", max_calls=3)
            check(f"lm on T3 with max_calls 3: ends by its own criterion ({r.calls} passes, {r.outcome}) and says the "
                  "limit was ignored", r.calls > 3 and not r.fell_back_to_migrad and r.outcome == "converged" and
                  any("lm does not honour a call limit" in n for n in r.notices), f"{r.notices}")
        pp = os.path.join(tmp, "pp")
        shutil.copytree(os.path.join(TESTS, "identical_pp_res"), pp, ignore=shutil.ignore_patterns("output", "checks",
                                                                                                    "expected"))
        for d in ("output", "checks"):
            os.makedirs(os.path.join(pp, d), exist_ok=True)
        with pyazr.Session(os.path.join(pp, "identical_pp_res.azr")) as s:
            records = []
            r = s.fit_with_engine("migrad", minos=True, progress=records.append)
            minos_records = [rec for rec in records if rec["stage"] == "minos"]
            sides = [r.minos_valid.get(n, {}).get(side, {}) for n in r.names for side in ("lower", "upper")]
            check(f"MINOS: an interval and the validity of both sides for each of the {len(r.names)} free parameters",
                  set(r.minos) == set(r.names) and set(r.minos_valid) == set(r.names) and
                  all(set(v) == {"valid", "at_limit", "call_limit", "new_minimum"} for v in sides), f"{r.minos_valid}")
            check("MINOS: every side valid here (a converged fit, no limits reached)",
                  all(v["valid"] and not v["call_limit"] and not v["new_minimum"] for v in sides), f"{r.minos_valid}")
            check(f"MINOS: its {len(minos_records)} progress records are marked error_analysis",
                  minos_records and all(rec["error_analysis"] for rec in minos_records) and
                  not any(rec["error_analysis"] for rec in records if rec["stage"] != "minos"))
    print("FAIL" if failures else "OK", f"({len(failures)} failure(s))" if failures else "")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
