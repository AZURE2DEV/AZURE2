#!/usr/bin/env python3
"""levelscan._patterns equal-split signs without materialising all 2^n sign vectors (2026-11-06).

Before: list(itertools.product((1, -1), repeat=n)) then a sample of --max-signs: ~75 GB at n=28, 300-650 GB at n=30-31,
so every levelscan with a >=30-channel J^pi was OOM-killed.
Checks: (1) the templates are IDENTICAL to the old enumerate-then-sample code for small n (same rng indices, same order);
(2) n = 31 and n = 40 finish at once in bounded memory, with --max-signs distinct sign vectors."""
import itertools, os, sys, time
try:
    import resource   # POSIX only: without it (Windows) the time bound is checked and the memory bound is not
except ImportError:
    resource = None
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
import numpy as np
from rmfit.levelscan import _patterns
ok = True
def check(name, cond, detail=""):
    global ok
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)

def old_equal_split(n, max_signs, seed=0):
    """the pre-2026-11-06 code, verbatim (only safe for small n)"""
    allsigns = list(itertools.product((1, -1), repeat=n))
    if len(allsigns) > max_signs:
        rng = np.random.default_rng(seed)
        allsigns = [allsigns[i] for i in sorted(rng.choice(len(allsigns), max_signs, replace=False))]
    return [(tuple([1.0 / n] * n), tuple(sg), "equal split") for sg in allsigns]

def main():
    same = True
    for n in (2, 3, 5, 8, 12, 16):
        for ms in (8, 32, 4096):
            new = [p for p in _patterns(n, ms) if p[2] == "equal split"]
            same &= new == old_equal_split(n, ms)
    check("equal-split templates identical to the old code for n = 2..16, max_signs 8/32/4096", same)
    rss = (lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) if resource else (lambda: 0)
    t0, m0 = time.time(), rss()
    p31 = _patterns(31, 32); p40 = _patterns(40, 32)
    dt, dm = time.time() - t0, (rss() - m0) / 1024
    eq31 = {s for _, s, d in p31 if d == "equal split"}
    check("n = 31 and 40 build in < 5 s and < 200 MB extra", dt < 5 and dm < 200, f"{dt:.1f} s, {dm:.0f} MB")
    check("n = 31: 32 distinct equal-split sign vectors of length 31", len(eq31) == 32 and all(len(s) == 31 for s in eq31), len(eq31))
    check("n = 40: 32 distinct equal-split sign vectors", len({s for _, s, d in p40 if d == "equal split"}) == 32)
    # 2026-11-26: n >= 63 overflowed rng.choice (a window with about 80 channels)
    for n in (62, 63, 64, 82, 120):
        pn = _patterns(n, 32)
        eq = [s for _, s, d in pn if d == "equal split"]
        check(f"n = {n}: 32 distinct equal-split sign vectors of length {n}, entries +-1",
              len(set(eq)) == 32 and len(eq) == 32 and all(len(s) == n and set(s) <= {1, -1} for s in eq), len(set(eq)))
    check("n = 82 is reproducible (same seed -> same templates)", _patterns(82, 32) == _patterns(82, 32))

if __name__ == "__main__":
    main()
    print("ALL OK" if ok else "FAILURES")
    sys.exit(0 if ok else 1)
