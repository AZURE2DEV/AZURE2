#!/usr/bin/env python3
"""local.sanity must flag same-J^pi levels that AZURE2 would merge on read (< 1 keV apart) and
warn about pairs within 1.5 keV.  Regression for 2026-09-19: a fit carried two same-J^pi levels
0.67 keV apart; when the engine was rebuilt with the stock 1 keV merge tolerance the incumbent
read back as a different model (301,445 -> 327,514) and an export failed.  Stub evaluator."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if not _PYAZR3:
    print("SKIP: a pyazr 3 test; this process runs pyazr 2 (the pyazr 2 version is test_sanity_same_jpi.py)")
    sys.exit(77)
from types import SimpleNamespace as NS
import numpy as np
import rmfit.enginepin as _enginepin
# these checks are the rule for an engine that merges same-J^pi levels within 1 keV on read (before dev e239a12);
# the engine of this tree keeps them apart (enginepin.keeps_close_levels), where a close pair is only a warning
_enginepin.keeps_close_levels = lambda *a, **k: False
from rmfit.local import sanity, BoundsPolicy


class K(NS):
    def __str__(self): return f"{self.kind}[{self.level}]"


def make(e_free):
    keys = [K(kind="energy", level="9/2+#2", pair=None), K(kind="energy", level="1/2+#1", pair=None)]
    ev = NS(keys=keys, wigner_amplitudes=lambda: {},
            level_energy={"9/2+#1": 9.8620, "9/2+#2": 9.8700, "9/2+#3": 10.89, "1/2+#1": 9.4906, "1/2+#2": 9.4920, "7/2-#1": 9.8622})
    return ev, np.array([e_free, 9.4906])


def main():
    ok = True
    def check(name, cond, detail=""):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)
    pol = BoundsPolicy(e_window=0.5)
    ev, x = make(9.862671)                      # the free 9/2+ sits 0.671 keV from the fixed one
    r = sanity(ev, x, pol, transform=False, x_ref=x)
    check("a 0.67 keV same-J^pi pair fails the check", (not r.ok) and len(r.merge_pairs) == 1, f"{r.merge_pairs} {r.summary()}")
    check("the summary names the pair", "9/2+#1" in r.summary() and "9/2+#2" in r.summary(), r.summary())
    check("the FREE energy is taken from x, not from the file", "9.862671" in " ".join(r.merge_pairs), str(r.merge_pairs))
    check("a 1.4 keV pair is a warning, not a failure", len(r.close_pairs) == 1 and "1/2+#1" in r.close_pairs[0], str(r.close_pairs))
    check("different J^pi at the same energy is not flagged", not any("7/2-" in m for m in r.merge_pairs + r.close_pairs))
    ev, x = make(9.8700)                        # moved 8 keV away
    r = sanity(ev, x, pol, transform=False, x_ref=x)
    check("with the pair 8 keV apart the check passes", r.ok and not r.merge_pairs, r.summary())
    check("and the 1.4 keV warning still shows in the ok summary", "WARNING" in r.summary(), r.summary())
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
