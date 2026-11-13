#!/usr/bin/env python3
"""The lock of a pyazr.Session: every call that enters the engine holds it (pyazr/locking.py).

(a) Every public method and property of pyazr.Session (its own and its mixins') that reaches the engine
    (its source names ``self.sess`` or ``self._sess``, or it calls another member that does) is wrapped by
    ``engine_call`` (``takes_lock``); the members that do not reach the engine are listed here by name.
(b) Two threads on one session of T1 at once -- one calling ``evaluate(x)`` and the old ``calculate_rwa(x)``, the
    other ``fit_objective()``'s value and gradient and ``segment_chi2(x)`` -- each at its own point, get what each gets
    alone, bit for bit.  (Without the lock two threads on one engine crash the process or mix their points.)

Run: python3 tests/pyazr/session_lock_test.py
"""
import inspect
import os
import re
import shutil
import sys
import tempfile
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    import pyazr
    from pyazr import Session
    pyazr.kinds()
except Exception as err:  # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)

# Public members that never reach the engine: the raw handle itself, values read from the options or the file.
NO_ENGINE = {"sess", "output_dir", "is_alive", "basis", "read_sav", "read_sav_basis"}


def members():
    """{name: function} of the public methods and property getters of Session (its MRO, without object)."""
    out = {}
    for klass in reversed(Session.__mro__[:-1]):
        for name, value in vars(klass).items():
            if name.startswith("_"):
                continue
            if isinstance(value, property):
                out[name] = value.fget
            elif isinstance(value, staticmethod):
                out[name] = value.__func__
            elif inspect.isfunction(value):
                out[name] = value
    return out


def source(fn):
    return inspect.getsource(inspect.unwrap(fn))


print("(a) every public member that reaches the engine takes the lock")
found = members()
reaches = {n for n, fn in found.items() if re.search(r"self\.(_?sess)\b", source(fn))}
for _ in range(5):  # a member that calls one that reaches the engine reaches it too
    reaches |= {n for n, fn in found.items() if any(re.search(rf"self\.{m}\b", source(fn)) for m in reaches)}
unlocked = sorted(n for n in reaches - NO_ENGINE if not getattr(found[n], "takes_lock", False))
check("every public member that reaches the engine is an engine_call", not unlocked, unlocked)
check("the members listed as not reaching the engine exist", NO_ENGINE <= set(found), sorted(NO_ENGINE - set(found)))
check("close() takes the lock (it frees the engine a running call uses)", getattr(found["close"], "takes_lock", False))

print("(b) two threads on one session")
tmp = tempfile.mkdtemp(prefix="session_lock_")
try:
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), os.path.join(tmp, "T1"),
                    ignore=shutil.ignore_patterns("output", "checks"))
    os.makedirs(os.path.join(tmp, "T1", "output"))
    with Session(os.path.join(tmp, "T1", "T1.azr")) as s:
        f = s.fit_objective()
        x0 = s.parameters.x0
        points = [x0, x0 * (1 + 1e-3 * np.cos(np.arange(x0.size)))]

        def first(x):
            ev = s.evaluate(x)
            return [ev.objective, ev.chi2] + [d.model for d in ev.data] + list(s.calculate_rwa(x))

        def second(x):
            v, g = f.value_and_grad(x)
            return [v, g, s.segment_chi2(x)]

        work = [first, second]
        alone = [work[i](points[i]) for i in range(2)]
        got = {}

        def run(i):
            got[i] = [work[i](points[i]) for _ in range(4)]

        threads = [threading.Thread(target=run, args=(i,)) for i in range(2)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        def same(a, b):
            return len(a) == len(b) and all(np.array_equal(np.asarray(u), np.asarray(v)) for u, v in zip(a, b))

        check("evaluate and calculate_rwa beside objective and segment_chi2 at once: each thread's values alone, "
              "bit for bit", all(same(r, alone[i]) for i in range(2) for r in got.get(i, [[]])),
              {i: len(got.get(i, [])) for i in range(2)})
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok: the session's lock")
