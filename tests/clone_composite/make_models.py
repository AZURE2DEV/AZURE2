#!/usr/bin/env python3
"""Models for tests/clone_composite, built at run time (check.sh) into a work directory.

From tests/target_density (12C(p,p) at 84.3 deg lab, segment 3, every level fixed, a thin target with a fitted areal
density), with segment 3 made a SUM of its own elastic term and one identical elastic component (2 sigma), so the
fitted density acts on a segment that holds components:

  true.azr   N = 1.1e18, the SUM segment alone (check.sh makes the fit data from its output)
  fit.azr    N free from 1e18 ("densityfit 1 0"), data from true.azr

    python3 make_models.py <tests/target_density> <workdir>
"""
import sys

src, work = sys.argv[1], sys.argv[2]
SUM_TAIL = " 1 0 1 1 1 -999 0"


def sum_segment3(text, data_dir):
    lines = text.split("\n")
    i0, i1 = lines.index("<segmentsData>"), lines.index("</segmentsData>")
    n3 = 0
    for k, i in enumerate(range(i0 + 1, i1), start=1):
        f = lines[i].split()
        if not f:
            continue
        if k == 3:
            f[0] = "1"
            f[14] = f[14].replace("fitdata/", "data/").replace("data/", data_dir + "/")
            if f[15:17] != ["0", "0"]:
                sys.exit("make_models: segment 3 is not a plain line")
            lines[i] = " ".join(f[:15]) + SUM_TAIL
            n3 += 1
        else:
            f[0] = "0"
            lines[i] = " ".join(f)
    if n3 != 1:
        sys.exit("make_models: segment 3 not found")
    j0, j1 = lines.index("<targetInt>"), lines.index("</targetInt>")
    keep = [l for l in lines[j0 + 1:j1] if l.split()[:2] == ["1", '"3"']]
    if len(keep) != 1:
        sys.exit("make_models: no <targetInt> line for segment 3")
    lines[j0 + 1:j1] = keep
    return "\n".join(lines)


open(work + "/true.azr", "w").write(sum_segment3(open(src + "/true.azr").read(), "data"))
open(work + "/fit.azr", "w").write(sum_segment3(open(src + "/fit.azr").read(), "fitdata"))
