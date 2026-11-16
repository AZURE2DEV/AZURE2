#!/usr/bin/env python3
"""structure.effective_k -- which of a new level's parameters tau charges for (2026-09-20).
Pure function.  Run:  python3 rmfit/tests/test_effective_k.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.structure import effective_k  # noqa: E402
failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok: failures.append(name)
W = lambda pair, L, S, t2, g: dict(kind="width", pair=pair, L=L, S=S, theta2=t2, gamma_ev=g)
def main():
    E = dict(kind="energy")
    # nothing seeded: identical to the old k_new
    p = [E, W(1, 2, 0.5, 0.02, 2000.0), W(2, 3, 0.5, 0.001, 50.0)]
    k = effective_k(p, [])
    check("no seeding: k_effective == k_all == old k_new", k["k_pattern"] == k["k_effective"] == k["k_all"] == 3, str(k))
    # a one-channel pattern + 20 seeded channels that stayed at their seeds (1 eV on a 3 keV level)
    idle = [W(4 + i % 6, i % 4, 0.5 + (i % 3), 1e-4, 1.0) for i in range(20)]
    seeds = [(q["pair"], q["L"], q["S"]) for q in idle]
    k = effective_k([E, W(1, 2, 0.5, 0.03, 3000.0)] + idle, seeds)
    check("idle seeded channels are not charged", (k["k_pattern"], k["k_effective"], k["k_all"]) == (2, 2, 22), str(k))
    # one seeded channel grew to 6 keV on a channel with a huge Wigner limit: theta2 tiny, but it IS used
    grown = [W(4, 0, 0.5, 1e-3, 6000.0)] + idle[1:]
    k = effective_k([E, W(1, 2, 0.5, 0.03, 3000.0)] + grown, [(4, 0, 0.5)] + seeds[1:])
    check("a seeded channel with a real share of the width is charged (theta2 floor alone would miss it)",
          k["k_effective"] == 3 and k["used"] == [(4, 0, 0.5)], str(k))
    # theta2 trigger on a weak channel
    k = effective_k([E, W(1, 2, 0.5, 0.03, 3000.0), W(8, 4, 2.5, 0.2, 5.0)], [(8, 4, 2.5)])
    check("theta2 >= 0.05 also counts as used", k["k_effective"] == 3, str(k))
    # pattern channels count even when tiny
    k = effective_k([E, W(1, 2, 0.5, 1e-9, 1e-6)], [])
    check("pattern channels always count", k["k_effective"] == 2, str(k))
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all effective-k checks passed")
if __name__ == "__main__":
    main()
