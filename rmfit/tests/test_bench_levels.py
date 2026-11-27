#!/usr/bin/env python3
"""bench.level_difference (2026-11-25): a run that beats its reference with a different level set is EXTRA-LEVELS, not
REFERENCE-SUSPECT.  Engine-free: two synthetic REVIEW.md level tables."""
import os, sys, tempfile
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
from rmfit.bench import bench as bb

HEAD = "## Levels inside the data range\n\n| level | Ex | Gamma | data sets | points | local Ex spacing | target-int segments |\n|---|---|---|---|---|---|---|\n"
REF = HEAD + "| 2+ | 12.9487 | 357 keV | a | 1 | | |\n| 3- | 13.1374 | 105 keV | a | 1 | | |\n\n## Figures\n| 9- | 1.0 | 1 keV | x | 1 | | |\n"
RUN = (HEAD + "| 2+ | 12.9482 | 358 keV | a | 1 | | |\n| 3- | 12.9897 | 233 eV | a | 1 | | |\n| 3- | 13.1371 | 104 keV | a | 1 | | |\n"
       "| 6+ | 11.4389 | ? | a | 1 | | |\n")
fails = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + ("" if ok else f"  -- {info}"))
    if not ok:
        fails.append(name)


with tempfile.TemporaryDirectory() as d:
    pr, pu = os.path.join(d, "ref.md"), os.path.join(d, "run.md")
    open(pr, "w").write(REF)
    open(pu, "w").write(RUN)
    lv = bb.review_levels(pr)
    check("only the level table is read (the row after ## Figures is not)", [x[0] for x in lv] == ["2+", "3-"], str(lv))
    extra, missing = bb.level_difference(pu, pr)
    check("extra = the narrow 3- only (matched within 30 keV, '?' = switched off)", extra == ["3- 12.9897 (233 eV)"], str(extra))
    check("nothing missing", missing == [], str(missing))
    check("same file: no difference", bb.level_difference(pr, pr) == ([], []))
    check("no REVIEW.md: (None, None)", bb.level_difference(os.path.join(d, "x.md"), pr) == (None, None))
print(f"FAILED {len(fails)}: " + ", ".join(fails) if fails else "all ok")
sys.exit(1 if fails else 0)
