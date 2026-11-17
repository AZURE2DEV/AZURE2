#!/usr/bin/env python3
"""A level add that takes its width from an unresolved same-J^pi neighbour is a SPLIT (2026-09-21).
Numbers are from real archive cases (widths read from the exported .azr files)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.structure import split_evidence, classify

ok = True
def check(name, cond, detail=""):
    global ok
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)

def main():
    # 1/2- added at 9.8637 (11.06 keV); the 1/2- 9.8693 went 43.93 -> 36.27 keV
    nb = [dict(label="1/2-#7", energy=9.86998, gamma_before_ev=43927.0, gamma_after_ev=36271.0),
          dict(label="1/2-#8", energy=10.1556, gamma_before_ev=2.0e5, gamma_after_ev=2.0e5)]
    s = split_evidence(9.86371, 11060.0, nb)
    check("the 9.8637 1/2- is a split of 1/2-#7", s is not None and s["label"] == "1/2-#7", s)
    import json, numpy as np
    sj = split_evidence(np.float64(9.86371), np.float64(11060.0), [dict(nb[0], gamma_before_ev=np.float64(43927.0))])
    check("the result is JSON-serialisable even from numpy inputs", bool(json.dumps(sj)))
    check("separation and lost width are reported", s and abs(s["sep_kev"] - 6.27) < 0.05 and abs(s["lost_ev"] - 7656) < 1 and 0.68 < s["lost_frac"] < 0.70, s)
    # 5/2- added at 10.886 (363.6 keV); the 5/2- 10.940 GAINED width 182.5 -> 485.9 keV: overlap, no transfer
    check("overlapping broad levels without width transfer are NOT a split",
          split_evidence(10.886, 363600.0, [dict(label="5/2-#9", energy=10.940, gamma_before_ev=182500.0, gamma_after_ev=485900.0)]) is None)
    check("a neighbour that lost width but is RESOLVED is not a split",
          split_evidence(9.80, 11060.0, [dict(label="x", energy=9.87, gamma_before_ev=43927.0, gamma_after_ev=30000.0)]) is None)
    check("a small loss (< half the new width) is not a split",
          split_evidence(9.86371, 11060.0, [dict(label="x", energy=9.86998, gamma_before_ev=43927.0, gamma_after_ev=40000.0)]) is None)
    check("transfer_frac is adjustable",
          split_evidence(9.86371, 11060.0, [dict(label="x", energy=9.86998, gamma_before_ev=43927.0, gamma_after_ev=40000.0)], transfer_frac=0.3) is not None)
    check("the narrow-doublet shape: 330 eV level 1.2 keV from a 3.8 keV level that gave up 300 eV",
          split_evidence(9.8620, 330.0, [dict(label="9/2+#2", energy=9.8632, gamma_before_ev=3800.0, gamma_after_ev=3500.0)]) is not None)
    pole = dict(label="1/2-#10", energy=12.0605, gamma_before_ev=2.4e6, gamma_after_ev=2.4e6 - 16093.0)
    check("a far, MeV-wide background level that shed 0.7% of itself is not the partner",
          split_evidence(9.86371, 11060.0, [pole]) is None)
    s2 = split_evidence(9.86371, 11060.0, [pole] + nb)
    check("... and the real partner is still found with it present", s2 is not None and s2["label"] == "1/2-#7", s2)
    # 2026-10-16 recall test (10-15-26_recall): the re-added 3/2+ 8.073 (82.5 keV) was called a split of the 3/2+ pole at
    # 50 MeV ("41927.2 keV away ... gave up 7.721e+06 eV"); the 9/2+ 10.33 add likewise named the 9/2+ 50 MeV pole
    pole50 = dict(label="3/2+#10", energy=50.0, gamma_before_ev=9.0e7, gamma_after_ev=9.0e7 - 7.721e6)
    check("a 50 MeV background pole is never a split partner (recall 3/2+ 8.073)", split_evidence(8.0728, 82540.0, [pole50]) is None)
    check("... 9/2+ 10.33 case", split_evidence(10.2906, 17780.0, [dict(label="9/2+#4", energy=50.0, gamma_before_ev=1.0e9,
          gamma_after_ev=1.0e9 - 7.677e7)]) is None)
    far = dict(label="x", energy=9.5, gamma_before_ev=2.0e6, gamma_after_ev=1.8e6)
    check("a broad level 1.4 MeV away that shed 10% of itself is not a partner (1 MeV cap)", split_evidence(8.1, 80000.0, [far]) is None)
    check("... unless the cap is raised", split_evidence(8.1, 80000.0, [far], max_sep_kev=2000.0) is not None)
    s3 = split_evidence(9.86371, 11060.0, [pole50] + nb)
    check("the real 9.8637 partner is still found with a 50 MeV pole present", s3 is not None and s3["label"] == "1/2-#7", s3)
    check("no width known -> no verdict", split_evidence(9.86, 0.0, nb) is None and split_evidence(9.86, None, nb) is None)
    check("no neighbours -> None", split_evidence(9.86, 1e4, []) is None)
    tiny = {"n0 L1": 8.7e-4}
    kw = dict(cap=1.0, tau=610.0, rho=0.39, level_gamma_ev=11060.0, narrow_width_ev=1e4)
    check("classify: accept without the flag (old behaviour)", classify(305047.8, 303818.1, tiny, 0.0003, **kw) == "accept")
    check("classify: 'splits' with it", classify(305047.8, 303818.1, tiny, 0.0003, splits=True, **kw) == "splits")
    check("classify: splits never masks marginal", classify(305047.8, 304900.0, tiny, 0.0003, splits=True, **kw) == "marginal")
    check("classify: splits never masks redistributes", classify(305047.8, 303818.1, tiny, 0.0003, splits=True, **dict(kw, rho=0.7)) == "redistributes")
    check("classify: splits comes before narrow", classify(305047.8, 303818.1, tiny, 0.0003, splits=True, **dict(kw, level_gamma_ev=500.0)) == "splits")
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
