#!/usr/bin/env python3
"""local.separate_close_levels / sanity's separable pairs (2026-10-10, the evaluator: give the broader of two near-degenerate
same-J^pi levels a small energy adjustment instead of fixing the energies).  Stub evaluator: transform returns the
vector with the width entries x 1e6, so a width key's value is its physical width in MeV (amplitude-sized, below
sanity's runaway limit) and the transform returns eV.  Run: python3 rmfit/tests/test_separate_close_levels.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from types import SimpleNamespace as NS
import numpy as np
from rmfit.local import separate_close_levels, sanity, BoundsPolicy


class K(NS):
    def __str__(self): return f"{self.kind}[{self.level}]"


def make(e_a, e_b, g_a, g_b, free_a=True, free_b=True, extra=()):
    """two 1/2+ levels A (#1) and B (#2); widths in eV; optional extra fixed levels in level_energy"""
    keys, x = [], []
    if free_a: keys.append(K(kind="energy", level="1/2+#1", pair=None)); x.append(e_a)
    if free_b: keys.append(K(kind="energy", level="1/2+#2", pair=None)); x.append(e_b)
    keys += [K(kind="width", level="1/2+#1", pair=1), K(kind="width", level="1/2+#2", pair=1)]; x += [g_a / 1e6, g_b / 1e6]
    le = {"1/2+#1": e_a, "1/2+#2": e_b, "3/2-#1": e_a}
    le.update(dict(extra))
    scale = np.array([1e6 if k.kind == "width" else 1.0 for k in keys])
    ev = NS(keys=keys, level_energy=le, wigner_amplitudes=lambda: {}, transform=lambda v: np.asarray(v, float) * scale)
    return ev, np.array(x, float)


def main():
    ok = True
    def check(name, cond, detail=""):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)

    # 1. a ~6 MeV level 0.5 keV BELOW a ~0.6 MeV level -> the broad one moves down to 2.5 keV
    ev, x = make(11.4760, 11.4765, 6.0e6, 6.0e5)
    y, moved, refused = separate_close_levels(ev, x)
    check("one level moved, nothing refused", len(moved) == 1 and not refused, f"{moved} {refused}")
    check("the BROADER level moved", moved and moved[0]["level"] == "1/2+#1", str(moved))
    check("it moved AWAY (down) to 2.5 keV separation", abs((11.4765 - y[0]) * 1e3 - 2.5) < 1e-9 and y[1] == 11.4765, str(y))
    check("the input vector is untouched", x[0] == 11.4760, str(x))
    check("other-J^pi levels at the same energy are ignored", all(m["other"] != "3/2-#1" for m in moved))

    # 2. broader level ABOVE: moves up
    ev, x = make(11.4760, 11.4765, 6.0e5, 6.0e6)
    y, moved, _ = separate_close_levels(ev, x)
    check("broader level above moves up", moved and moved[0]["level"] == "1/2+#2" and abs((y[1] - 11.4760) * 1e3 - 2.5) < 1e-9, f"{moved} {y}")

    # 3. a narrow pair (the 9/2+ doublet kind: keV widths) is refused, nothing moves
    ev, x = make(9.8620, 9.86267, 4400.0, 330.0)
    y, moved, refused = separate_close_levels(ev, x)
    check("narrow pair refused, vector unchanged", not moved and len(refused) == 1 and np.array_equal(x, y), f"{moved} {refused}")

    # 4. broader level's energy fixed -> the other moves if IT is broad enough, else refused
    ev, x = make(11.4760, 11.4765, 6.0e6, 6.0e5, free_a=False)
    y, moved, refused = separate_close_levels(ev, x)
    check("broad member fixed, narrower (600 keV) free -> narrower moves up", moved and moved[0]["level"] == "1/2+#2", f"{moved} {refused}")
    ev, x = make(11.4760, 11.4765, 6.0e6, 6.0e4, free_a=False)
    y, moved, refused = separate_close_levels(ev, x)
    check("broad member fixed, other too narrow (60 keV) -> refused", not moved and len(refused) == 1, f"{moved} {refused}")
    ev, x = make(11.475997, 11.477468, 6.0e6, 6.0e5, free_a=False, free_b=False)
    y, moved, refused = separate_close_levels(ev, x)
    check("both energies fixed (today's C0/A/B case, 1.47 keV) -> nothing moves, refusal reported", not moved and len(refused) == 1, f"{refused}")

    # 5. pairs >= 1.5 keV apart are left alone
    ev, x = make(11.4760, 11.4776, 6.0e6, 6.0e5)
    y, moved, refused = separate_close_levels(ev, x)
    check("1.6 keV apart: untouched", not moved and not refused and np.array_equal(x, y))

    # 6. sanity: a broad pair < 1 keV is separable (ok), a narrow one still fails
    pol = BoundsPolicy(e_window=0.5)
    ev, x = make(11.4760, 11.4765, 6.0e6, 6.0e5)
    r = sanity(ev, x, pol, transform=True, x_ref=x)
    check("sanity: broad 0.5 keV pair is ok and listed as separable", r.ok and len(r.separable_pairs) == 1 and not r.merge_pairs, r.summary())
    ev, x = make(9.8620, 9.86267, 4400.0, 330.0)
    r = sanity(ev, x, pol, transform=True, x_ref=x)
    check("sanity: narrow 0.67 keV pair still FAILS", (not r.ok) and len(r.merge_pairs) == 1, r.summary())

    # 7. a chain of three within 1.5 keV ends with every gap >= 2.5 keV or refused
    keys = [K(kind="energy", level=f"1/2+#{i}", pair=None) for i in (1, 2, 3)] + [K(kind="width", level=f"1/2+#{i}", pair=1) for i in (1, 2, 3)]
    x = np.array([11.4760, 11.4765, 11.4770, 5.0, 4.0, 3.0])
    sc3 = np.array([1.0] * 3 + [1e6] * 3)
    ev = NS(keys=keys, level_energy={f"1/2+#{i}": e for i, e in zip((1, 2, 3), x[:3])}, wigner_amplitudes=lambda: {}, transform=lambda v: np.asarray(v, float) * sc3)
    y, moved, refused = separate_close_levels(ev, x)
    es = np.sort(y[:3]); gaps = np.diff(es) * 1e3
    check("chain of three: all gaps >= 1.5 keV afterwards", np.all(gaps >= 1.5 - 1e-9), f"{gaps} {moved} {refused}")
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
