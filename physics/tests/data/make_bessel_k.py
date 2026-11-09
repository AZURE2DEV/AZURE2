#!/usr/bin/env python3
"""Reference values of e^x K_n(x) from mpmath (40 digits), for tests/bound_functions_test.cpp.

Grid: n = 0 ... 21 (S(0) of l <= 10 needs K_2l and K_2l+1); x = 1e-3 ... 500 (17 points).
Usage: make_bessel_k.py OUTFILE
"""
import sys

import mpmath as mp

mp.mp.dps = 40
XS = ["1e-3", "3e-3", "0.01", "0.03", "0.1", "0.3", "0.75", "1", "2", "3", "5", "10", "20", "30", "50", "100", "500"]
with open(sys.argv[1], "w") as out:
    out.write("# n x exp(x)*K_n(x)  (mpmath, 40 digits)\n")
    for n in range(0, 22):
        for xs in XS:
            x = mp.mpf(xs)
            out.write(f"{n} {xs} {mp.nstr(mp.besselk(n, x) * mp.exp(x), 20, min_fixed=0, max_fixed=0)}\n")
