"""Unit check of the flip_floor filter (search.py, 2026-09-16): a flip whose FROZEN response
is numerical noise must not take a staged-polish slot.  Numbers measured on a real
sign round: the six excluded flips are amplitudes in particle pairs that no active
data segment observes (pairs 4, 6 here) or high-L channels of the new level -- one of them,
g[9/2+#6 p4 L4 S0.5], has theta^2 = 0.98 and still moves the objective by only 0.24, which is
why an amplitude-size threshold (sign_eps) cannot catch them.  No fit, no engine: runs in ms."""
import sys, types
import os; sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from rmfit.search import SearchConfig, tau_flip
obj_inc, npts = 349629.49, 34951
cfg = SearchConfig()
tau = tau_flip(obj_inc, npts, cfg)
# measured frozen responses of the 9/2+ round (flip - incumbent)
cases = [("g[9/2+#5 p1 L4 S0.5]", 1866.04), ("g[9/2+#5 p2 L5 S0.5]", 1865.79),
         ("g[9/2+#6 p1 L4 S0.5]", 442.32), ("g[9/2+#6 p2 L5 S0.5]", 591.05),
         ("g[9/2+#6 p5 L1 S3.5]", 461.46),
         ("g[9/2+#5 p4 L4 S0.5]", 0.2436), ("g[9/2+#6 p3 L4 S0.5]", 0.1339),
         ("g[9/2+#6 p4 L4 S0.5]", 0.2436), ("g[9/2+#6 p5 L3 S2.5]", 0.0096),
         ("g[9/2+#6 p5 L3 S3.5]", 0.0164), ("g[9/2+#6 p6 L3 S1.5]", 0.0000)]
scr = [types.SimpleNamespace(move=types.SimpleNamespace(description=d, kind="flip"), obj_flip=obj_inc + dv, obj_gn=obj_inc + dv) for d, dv in cases]
noise = [s for s in scr if abs(s.obj_flip - obj_inc) < cfg.flip_floor * tau]
kept = [s for s in scr if s not in noise]
print(f"tau {tau:.1f}; floor {cfg.flip_floor:g} x tau = {cfg.flip_floor*tau:.2f}")
print("excluded:", [s.move.description for s in noise])
print("kept    :", [s.move.description for s in kept])
ok = (len(noise) == 6 and all(abs(s.obj_flip-obj_inc) < 1.0 for s in noise)
      and all(abs(s.obj_flip-obj_inc) > 100 for s in kept))
print("ALL OK" if ok else "FAILURES: flip_floor selection")
sys.exit(0 if ok else 1)
