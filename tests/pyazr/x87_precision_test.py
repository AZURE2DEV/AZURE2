#!/usr/bin/env python3
"""A pyazr call computes the same bits in any thread, whatever the thread's x87 precision
(support/include/azure2/support/FloatingPointPrecision.h; api/src/pybind_engine_call.h).

On x86 `long double` is the x87 format, and the mantissa an operation keeps is the precision control of the thread's
control word: 64 bits (0x37F) or 53 (0x27F).  On Windows a Python thread starts with 53 bits and the main thread with
64 (MinGW's start-up), and the same call of T1 gave other bits in a threading.Thread than in the main thread
(session_lock_test.py (b), objective_object_test.py (d); 2026-11-16).  The engine now takes the 64-bit precision for
the length of every call from the binding and in every thread of its OpenMP regions.

(a) What the threads hold (printed for the CI's logs: the precision a Python thread and the OpenMP runtime's threads
    start with on this system): the main thread's control word, a threading.Thread's, and those of the threads of
    an OpenMP region of each, as found and inside the engine's scope.  Checked: inside every scope the precision is
    64 bits, and the calling thread's word is as it was after the call.
(b) T1 at a point: the objective's value and gradient, its rows and Jacobian, and evaluate(x)'s chi^2 and model
    values, computed in the main thread, equal bit for bit those computed
      - in a threading.Thread as it starts,
      - in a threading.Thread whose precision is lowered to 53 bits before the calls (as on Windows; on Linux the
        OpenMP threads that thread starts inherit the lowered word, so the regions' scopes are tested too),
      - in the main thread with its precision lowered to 53 bits (its OpenMP threads, made earlier, keep 64),
    and the lowered word is still in force after the calls.
Without an x87 (arm64) the control words are 0 and (b) compares the threads as they start only.

Can fail: with the build before the scope (trial 101), lowering a thread's precision with a ctypes helper moves 27 of
T1's 442 residuals by up to 1.4e-12 relative, its value by 8e-16 and 7 of 13 gradient components (2026-11-16).

Run: python3 tests/pyazr/x87_precision_test.py
"""
import os
import shutil
import sys
import tempfile
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
LOWERED = 0x27F  # precision control 10: a 53-bit mantissa, all exceptions masked, round to nearest
PC_BITS, PC_64 = 0x300, 0x300

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import Session
    from pyazr import _azure2
    _azure2._x87_control_words
except Exception as err:  # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def in_thread(fn):
    out = {}

    def run():
        try:
            out["value"] = fn()
        except BaseException as e:  # reported by the caller
            out["error"] = e
    t = threading.Thread(target=run)
    t.start()
    t.join()
    if "error" in out:
        raise out["error"]
    return out["value"]


def words_text(w):
    return (f"calling {w['calling']:#06x}, in the call {w['in_call']:#06x}, after {w['after_call']:#06x}; OpenMP "
            f"threads found {[hex(v) for v in w['omp_found']]}, in their scope {[hex(v) for v in w['omp_in_region']]}")


print("(a) the x87 control words")
main_words = _azure2._x87_control_words()
worker_words = in_thread(_azure2._x87_control_words)
x87 = main_words["calling"] != 0
print(f"  main thread: {words_text(main_words)}")
print(f"  a threading.Thread: {words_text(worker_words)}")
if x87:
    for label, w in (("main thread", main_words), ("a threading.Thread", worker_words)):
        inside = [w["in_call"]] + list(w["omp_in_region"])
        check(f"{label}: 64-bit precision in the call and in every thread of the region",
              all(v & PC_BITS == PC_64 for v in inside), words_text(w))
        check(f"{label}: its control word after the call is the word before it", w["after_call"] == w["calling"],
              words_text(w))
    lowered_words = in_thread(lambda: (_azure2._set_x87_control_word(LOWERED), _azure2._x87_control_words())[1])
    print(f"  a threading.Thread lowered to 53 bits: {words_text(lowered_words)}")
    check("a thread lowered to 53 bits: 64 bits in the call and in every thread of the region, 53 after the call",
          all(v & PC_BITS == PC_64 for v in [lowered_words["in_call"]] + list(lowered_words["omp_in_region"]))
          and lowered_words["calling"] == LOWERED and lowered_words["after_call"] == LOWERED, words_text(lowered_words))
else:
    print("  no x87 on this machine: the words are 0")

print("(b) T1 in the main thread and in other threads, bit for bit")
tmp = tempfile.mkdtemp(prefix="x87_precision_")
try:
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), os.path.join(tmp, "T1"),
                    ignore=shutil.ignore_patterns("output", "checks"))
    os.makedirs(os.path.join(tmp, "T1", "output"), exist_ok=True)
    with Session(os.path.join(tmp, "T1", "T1.azr")) as s:
        f = s.objective()
        x0 = s.parameters.x0
        x = x0 * (1 + 1e-3 * np.cos(np.arange(x0.size)))

        def compute():
            v, g = f.value_and_grad(x)
            r, J = f.residuals_and_jacobian(x)
            ev = s.evaluate(x)
            return {"value": np.array([v]), "gradient": np.asarray(g), "rows": np.asarray(r), "Jacobian": np.asarray(J),
                    "evaluate chi2": np.array([ev.chi2]),
                    "evaluate model": np.concatenate([np.ravel(d.model) for d in ev.data])}

        def lowered():
            _azure2._set_x87_control_word(LOWERED)
            out = compute()
            out["word after"] = _azure2._x87_control_words()["calling"]
            return out

        def differences(a, b):
            """'' when every entry of a equals b's bit for bit, else what differs and by how much."""
            text = []
            for k in b:
                u, v = np.ravel(a[k]).astype(float), np.ravel(b[k]).astype(float)
                if u.shape != v.shape or u.tobytes() != v.tobytes():
                    n = int(np.sum(u != v)) if u.shape == v.shape else -1
                    rel = float(np.max(np.abs(u - v) / np.maximum(np.abs(v), 1e-300))) if n > 0 else 0.0
                    text.append(f"{k}: {n} of {v.size} differ, largest relative {rel:.2e}")
            return "; ".join(text)

        main = compute()
        check("T1 has rows and a Jacobian to compare", main["rows"].size > 100 and main["Jacobian"].size > 1000,
              {k: v.size for k, v in main.items()})
        d = differences(in_thread(compute), main)
        check("a threading.Thread as it starts = the main thread", not d, d)
        if x87:
            w = in_thread(lowered)
            d = differences(w, main)
            check("a threading.Thread lowered to 53 bits = the main thread", not d, d)
            check("  and its word is still the lowered one after the calls", w["word after"] == LOWERED,
                  hex(w["word after"]))
            found = _azure2._x87_control_words()["calling"]
            try:
                m = lowered()
            finally:
                _azure2._set_x87_control_word(found)
            d = differences(m, main)
            check("the main thread lowered to 53 bits = the main thread at 64", not d, d)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok: the same bits in every thread")
