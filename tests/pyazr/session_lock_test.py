#!/usr/bin/env python3
"""The lock of a pyazr.Session: every call that enters the engine holds it (pyazr/locking.py).

(a) Every public method and property of pyazr.Session (its own and its mixins') that reaches the engine
    (its source names ``self._engine`` or ``self._sess``, or it calls another member that does) is wrapped by
    ``engine_call`` (``takes_lock``); the members that do not reach the engine are listed here by name.
(a') Every use of the raw handle ``Session._engine`` in the other modules of ``pyazr/`` is inside a function wrapped
    by ``engine_call`` or a ``with session_lock(...)`` (or an ``Objective``'s ``with self._lock``) block (by their
    syntax trees).
(b) Two threads on one session of T1 at once -- one calling ``evaluate(x)`` twice, the other ``objective()``'s
    value and gradient, rows, and rows with their Jacobian -- each at its own point, get what each gets alone, bit
    for bit.  (Without the lock two threads on one engine crash the process or mix their points.)

Run: python3 tests/pyazr/session_lock_test.py
"""
import ast
import glob
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
NO_ENGINE = {"output_dir", "is_alive", "basis"}


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
reaches = {n for n, fn in found.items() if re.search(r"self\.(_sess|_engine)\b", source(fn))}
for _ in range(5):  # a member that calls one that reaches the engine reaches it too
    reaches |= {n for n, fn in found.items() if any(re.search(rf"self\.{m}\b", source(fn)) for m in reaches)}
unlocked = sorted(n for n in reaches - NO_ENGINE if not getattr(found[n], "takes_lock", False))
check("every public member that reaches the engine is an engine_call", not unlocked, unlocked)
check("the members listed as not reaching the engine exist", NO_ENGINE <= set(found), sorted(NO_ENGINE - set(found)))
check("close() takes the lock (it frees the engine a running call uses)", getattr(found["close"], "takes_lock", False))
check("the raw handle is private (Session._engine; no public sess)", "sess" not in found and "_engine" in dir(Session))

print("(a') every use of the raw handle in pyazr's other modules holds the lock")
SESSION_CLASSES = {"Session", "_HybridPotentials", "SessionPhysics", "SessionTHM"}  # check (a) covers them


def unlocked_uses(path):
    """The lines of `path` that read `<x>._engine` outside a class of the session, outside a function wrapped by
    engine_call and outside a `with session_lock(...)` / `with self._lock` block."""
    tree = ast.parse(open(path, encoding="utf-8").read())
    bad = []

    def locked_with(node):
        for item in node.items:
            e = item.context_expr
            if isinstance(e, ast.Call) and getattr(e.func, "id", getattr(e.func, "attr", "")) == "session_lock":
                return True
            if isinstance(e, ast.Attribute) and e.attr == "_lock":
                return True
        return False

    def visit(node, inside):
        if isinstance(node, ast.ClassDef) and node.name in SESSION_CLASSES:
            return
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if any(getattr(d, "id", getattr(d, "attr", "")) == "engine_call" for d in node.decorator_list):
                inside = True
        if isinstance(node, (ast.With, ast.AsyncWith)) and locked_with(node):
            inside = True
        if isinstance(node, ast.Attribute) and node.attr == "_engine" and not inside:
            bad.append(f"{os.path.relpath(path, ROOT)}:{node.lineno}")
        for child in ast.iter_child_nodes(node):
            visit(child, inside)

    visit(tree, False)
    return bad


modules = sorted(glob.glob(os.path.join(ROOT, "pyazr", "**", "*.py"), recursive=True))
bad = [u for path in modules for u in unlocked_uses(path)]
check(f"no use of Session._engine outside the lock in {len(modules)} modules of pyazr/", not bad, bad)

print("(b) two threads on one session")
tmp = tempfile.mkdtemp(prefix="session_lock_")
try:
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), os.path.join(tmp, "T1"),
                    ignore=shutil.ignore_patterns("output", "checks"))
    os.makedirs(os.path.join(tmp, "T1", "output"))
    with Session(os.path.join(tmp, "T1", "T1.azr")) as s:
        f = s.objective()
        x0 = s.parameters.x0
        points = [x0, x0 * (1 + 1e-3 * np.cos(np.arange(x0.size)))]

        def first(x):
            ev = s.evaluate(x)
            return [ev.objective, ev.chi2] + [d.model for d in ev.data] + [d.chi2 for d in s.evaluate(x).data]

        def second(x):
            v, g = f.value_and_grad(x)
            r, J = f.residuals_and_jacobian(x)
            return [v, g, f.residuals(x), r, J]

        work = [first, second]
        alone = [work[i](points[i]) for i in range(2)]

        def same(a, b):
            return len(a) == len(b) and all(np.array_equal(np.asarray(u), np.asarray(v)) for u, v in zip(a, b))

        def differing(a, b):
            """Which of the values differ, and the largest relative difference of each (the failure's detail)."""
            out = {}
            for k, (u, v) in enumerate(zip(a, b)):
                u, v = np.ravel(np.asarray(u, float)), np.ravel(np.asarray(v, float))
                if u.shape != v.shape:
                    out[k] = "shape"
                elif not np.array_equal(u, v):
                    rel = float(np.max(np.abs(u - v) / np.maximum(np.abs(v), 1e-300)))
                    out[k] = f"{int(np.sum(u != v))}/{v.size}, {rel:.1e}"
            return out

        # each thread's work alone in a threading.Thread first: a difference here is the thread's, not the two at
        # once (on Windows a Python thread's x87 precision was 53 bits: azure2/support/FloatingPointPrecision.h)
        one = {}
        for i in range(2):
            t = threading.Thread(target=lambda i=i: one.__setitem__(i, work[i](points[i])))
            t.start()
            t.join()
        check("each thread's values alone in a threading.Thread: the main thread's, bit for bit",
              all(same(one.get(i, []), alone[i]) for i in range(2)),
              {i: differing(one[i], alone[i]) for i in one})
        got = {}

        def run(i):
            got[i] = [work[i](points[i]) for _ in range(4)]

        threads = [threading.Thread(target=run, args=(i,)) for i in range(2)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        check("evaluate twice beside the objective's value, gradient, rows and Jacobian at once: each thread's values alone, "
              "bit for bit", all(same(r, alone[i]) for i in range(2) for r in got.get(i, [[]])),
              {i: [differing(r, alone[i]) for r in got.get(i, [])] for i in range(2)})
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok: the session's lock")
