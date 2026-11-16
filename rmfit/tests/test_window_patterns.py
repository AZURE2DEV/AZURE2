"""structure.scan_jpi_on_window builds only the 16 sign patterns it screens (2026-11-13: building all 2^(n-1) first ran out of
memory for n ~ 40 in one model).  Pure python: the patched expression equals the first 16 of the old full list, and n = 60 is
instant.  Run: python3 rmfit/tests/test_window_patterns.py"""
import time
def old(n):
    return [[1] * n] if n == 1 else [[1] + [1 if (m >> i) & 1 else -1 for i in range(n - 1)] for m in range(2 ** (n - 1))]
def new(n):
    return [[1] * n] if n == 1 else [[1] + [1 if (m >> i) & 1 else -1 for i in range(n - 1)] for m in range(min(2 ** (n - 1), 16))]
ok = all(new(n) == old(n)[:16] for n in range(1, 19))
t = time.time(); p = new(60); dt = time.time() - t
ok = ok and len(p) == 16 and len(p[0]) == 60 and dt < 0.1
import inspect, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
src = inspect.getsource(__import__("rmfit.structure", fromlist=["x"]).scan_jpi_on_window)
ok = ok and "min(2 ** (n - 1), 16)" in src
print("test_window_patterns:", "ok" if ok else "FAILED", f"(n=60 in {dt*1e3:.2f} ms)")
sys.exit(0 if ok else 1)
