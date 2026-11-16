"""Ledger.population / duplicate_of compare sign strings recomputed from the saved vectors
(2026-11-12), so a ledger whose signs column holds the pre-fix canonical form, or a mix,
still finds gauge-equivalent candidates.  Pure python.  Run: python3 rmfit/tests/test_ledger_signs.py"""
import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
from rmfit.keys import ParamKey, flip
from rmfit.ledger import Ledger

keys = [ParamKey("width", f"1/2+#{l}", p, 0, 0.5) for l in (1, 2, 3) for p in (1, 2, 3)]
x = np.array([1., -2, 3, 4, 5, -6, 7, 8, 9])
gauge = flip(x, keys, [k for k in keys if k.level == "1/2+#2"] + [k for k in keys if k.pair == 3])   # level + pair flip
other = flip(x, keys, [keys[4]])                                                                      # a genuine single flip
with tempfile.TemporaryDirectory() as d:
    L = Ledger(os.path.join(d, "ledger.sqlite"))
    sid = L.add_structure("s", os.path.join(d, "m.azr"), keys)
    class S:                                    # a minimal Score
        def __init__(self, o): self.objective = o
    c1 = L.add_candidate(sid, x, keys); L.record(c1, S(100.0), "verified")
    c2 = L.add_candidate(sid, other, keys); L.record(c2, S(100.01), "verified")
    L.db.execute("UPDATE candidates SET signs='stale-old-form' WHERE id=?", (c1,)); L.db.commit()
    ok = []
    ok.append(L.duplicate_of(sid, gauge, keys, 100.0) == c1)
    ok.append(L.duplicate_of(sid, other, keys, 100.0) == c2)
    ok.append(L.duplicate_of(sid, gauge, keys, 120.0) is None)
    c3 = L.add_candidate(sid, gauge, keys); L.record(c3, S(99.0), "verified")
    pop = L.population(sid, k=4)
    ok.append([p[0] for p in pop] == [c3, c2])          # c1 is c3's pattern: deduplicated
    L.close()
print("test_ledger_signs:", "ok" if all(ok) else f"FAILED {ok}")
sys.exit(0 if all(ok) else 1)
