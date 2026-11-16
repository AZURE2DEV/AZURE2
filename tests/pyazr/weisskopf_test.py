#!/usr/bin/env python3
"""pyazr.weisskopf_width: the engine's Weisskopf estimate (physics/Widths.h, the row ``weisskopf_width`` of
pyazr.physics) for one transition.

(1) Against the closed form it tabulates, computed here from the constants:
Gamma_W = hbar lambda with lambda(XL) = 8 pi (L + 1) / (hbar L [(2L+1)!!]^2) (E/hbar c)^(2L+1) B_W(XL), the
Weisskopf reduced strengths B_W(EL) = (1/4 pi) (3/(L+3))^2 R^(2L) e^2 and B_W(ML) = (10/pi) (3/(L+3))^2 R^(2L-2)
mu_N^2, R = 1.2 A^(1/3) fm (Bohr and Mottelson, vol. 1, eqs. 3C-38 and 3C-39).  The table's coefficients have two
significant digits: E1-E3 and M1-M3 agree to 4 %, E4 and M4 to 8 % (the table rounds them furthest).
(2) Against the Python formula pyazr had before the estimate became a function of physics/ (its table of
coefficients, kept here): the same numbers, bit for bit, on 648 transitions.
(3) A multipolarity outside E1-E4, M1-M4 or a transition energy <= 0 gives None; the top-level name and the row
give one number.

Run from anywhere:  python3 tests/pyazr/weisskopf_test.py
"""
import math
import os
import sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, ROOT)
try:
    import pyazr
    from pyazr.widths import weisskopf_width
except ImportError as e:  # pyazr needs its engine module to import
    print(f"skip: {e}")
    sys.exit(0)

ALPHA, HBARC, MP = 1 / 137.035999084, 197.3269804, 938.27208816   # -, MeV fm, MeV (CODATA 2018)


def closed_form(kind, L, e_gamma, A):
    """The Weisskopf width in eV from the formula above."""
    R = 1.2 * A ** (1 / 3)
    rate = (8 * math.pi * (L + 1) / (L * math.prod(range(2 * L + 1, 0, -2)) ** 2) * (e_gamma / HBARC) ** (2 * L + 1)
            * (3 / (L + 3)) ** 2)
    if kind == "E":
        width = rate * ALPHA * HBARC * R ** (2 * L) / (4 * math.pi)
    else:
        width = rate * (10 / math.pi) * ALPHA * HBARC * (HBARC / (2 * MP)) ** 2 * R ** (2 * L - 2)
    return width * 1e6


failures = []
for kind in "EM":
    for L in range(1, 5):
        tol = 0.04 if L < 4 else 0.08
        worst = max(abs(weisskopf_width(kind, L, e, A) / closed_form(kind, L, e, A) - 1)
                    for A in (4, 13, 27, 100) for e in (0.5, 2.0, 8.0))
        ok = worst <= tol
        print(f"  {'ok  ' if ok else 'FAIL'}  {kind}{L}: largest relative difference {worst:.3f} (<= {tol})")
        if not ok:
            failures.append(f"{kind}{L}")
# (2) the formula pyazr had (pyazr/widths.py before pyazr 3.0.0.dev5): eV, E in MeV
BEFORE = {
    ("E", 1): lambda A, e: 6.8e-2 * A ** (2 / 3) * e ** 3,
    ("E", 2): lambda A, e: 4.9e-8 * A ** (4 / 3) * e ** 5,
    ("E", 3): lambda A, e: 2.3e-14 * A ** 2 * e ** 7,
    ("E", 4): lambda A, e: 6.8e-21 * A ** (8 / 3) * e ** 9,
    ("M", 1): lambda A, e: 2.1e-2 * e ** 3,
    ("M", 2): lambda A, e: 1.5e-8 * A ** (2 / 3) * e ** 5,
    ("M", 3): lambda A, e: 6.8e-15 * A ** (4 / 3) * e ** 7,
    ("M", 4): lambda A, e: 2.0e-21 * A ** 2 * e ** 9,
}
cases = [(k, L, e, A) for (k, L) in BEFORE for A in (4, 7, 13, 14, 16, 27, 100, 208, 238)
         for e in (0.01, 0.137, 0.5, 1.0, 2.3654, 4.439, 8.0, 15.1, 20.0)]
worst = max(abs(weisskopf_width(k, L, e, A) / BEFORE[(k, L)](float(A), e) - 1) for k, L, e, A in cases)
ok = worst <= 1e-15
print(f"  {'ok  ' if ok else 'FAIL'}  the numbers of the Python formula it replaces, {len(cases)} transitions: largest "
      f"relative difference {worst:.1e} (<= 1e-15)")
if not ok:
    failures.append("the earlier formula")
ok = (weisskopf_width("M", 2, 4.439, 13) == pyazr.physics.weisskopf_width(4.439, 2, 13, "M")
      and pyazr.weisskopf_width is weisskopf_width)
print(f"  {'ok  ' if ok else 'FAIL'}  pyazr.weisskopf_width('M', 2, 4.439, 13) is the row's number")
if not ok:
    failures.append("the row")
for args in (("E", 5, 1.0, 13), ("M", 0, 1.0, 13), ("E", 1, 0.0, 13), ("E", 1, None, 13), ("X", 1, 1.0, 13)):
    ok = weisskopf_width(*args) is None
    print(f"  {'ok  ' if ok else 'FAIL'}  weisskopf_width{args} is None")
    if not ok:
        failures.append(str(args))
if failures:
    print(f"{len(failures)} check(s) FAILED: {failures}")
    sys.exit(1)
print("ok: all checks passed")
