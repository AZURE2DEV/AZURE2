#!/usr/bin/env python3
"""Writes mcn.azr: two ODD-A compound nuclei with the same J^pi levels (2026-10-30, rmfit multi-CN level labels).

Built with the helpers of AZURE2 tests/natural_target/make_models.py (the three-level p+12C skeleton):
  pair 1  p + 12C -> 13N (S_p 1.9435), levels 1/2+, 3/2-, 5/2+ at the skeleton's c.m. energies
  pair 2  p + 14C -> 15N (S_p 10.2074), the SAME J^pi set at other c.m. energies
so the engine has J-groups 1/2+ (13N), 1/2+ (15N), ... each numbered from 1: "1/2+#1" names two levels.
Energies and widths free; one data segment per nucleus (angle-integrated, data_13n.dat / data_15n.dat: a flat 1 b
+- 1 %, only the keys and level table matter).  No <abundances>, no natural keyword.
    AZURE2_ROOT=<tree> python3 make_fixture.py
"""
import importlib.util, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("AZURE2_ROOT", os.path.expanduser("~/AZURE2"))
spec = importlib.util.spec_from_file_location("mm", os.path.join(ROOT, "tests", "natural_target", "make_models.py"))
mm = importlib.util.module_from_spec(spec); spec.loader.exec_module(mm)
PAIR_C = dict(key=2, j1=0.5, pi1=1, j2=0.0, pi2=1, e2=0.0, m1=1.00782, m2=14.003242, z1=1, z2=6, sep=10.2074, rad=5.0)
LEVELS_C = [(0.5, 1, 0.80, 0.40), (1.5, -1, 1.30, 0.30), (2.5, 1, 1.70, 0.25)]
lv, _ = mm.levels([(mm.PAIR_A, mm.LEVELS_A), (PAIR_C, LEVELS_C)])
seg = ("1 1 1 1.0 3.0 0 0 0 1 0 0 0 0 0 data_13n.dat\n"
       "1 2 2 1.0 3.0 0 0 0 1 0 0 0 0 0 data_15n.dat")
text = mm.model(lv, "", segments_data=seg)
open(os.path.join(HERE, "mcn.azr"), "w").write(text)
for f in ("data_13n.dat", "data_15n.dat"):
    with open(os.path.join(HERE, f), "w") as fh:
        for i in range(41):
            e = 1.0 + 0.05 * i
            fh.write(f"{e:.4f} 0 1.0 0.01\n")
print("wrote mcn.azr, data_13n.dat, data_15n.dat")
