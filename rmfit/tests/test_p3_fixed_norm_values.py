#!/usr/bin/env python3
"""Fixed normalizations on pyazr 3 (2026-11-26, an engine migration): a FIXED norm has no parameter row in a
pyazr 3 session, so SegmentTable.with_nuisances never set it and rmfit read 1.0 for a fixed shared record
(``normalization SI`` / ``value 1.9537733 fixed``) and for a data set's own ``norm 0.8 fixed``; a fixed relvector reference
(model x norm) was then divided by the wrong norm (about 2x, and a fixed-reference objective below the old engine's).
Fixture fixtures_p3/relvector_p12c/ngrel.azr, edited: record B -> ``value 1.9537733 fixed``, data set p12c_90 -> ``norm 0.8
fixed``.  Checks: (a) SegmentTable.from_file gives those values (no session); (b) a SingleSession's segments and
_cov_norm_rows give them on the rows of those data sets, and the free record A keeps its session value."""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if not _PYAZR3:
    print("SKIP: a pyazr 3 test (pyazr 2 reads <segmentsData> norms directly)")
    sys.exit(77)
from _fixtures import fixture  # noqa: E402
import numpy as np
from rmfit.azr3 import SegmentTable
from rmfit.engine import SingleSession
failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)
def main():
    src = fixture("relvector_p12c")
    tmp = tempfile.mkdtemp(prefix="p3fixnorm_")
    try:
        os.symlink(os.path.join(src, "data"), os.path.join(tmp, "data"))
        for d in ("output", "checks"):
            os.makedirs(os.path.join(tmp, d))
        t = open(os.path.join(src, "ngrel.azr")).read()
        a = "normalization B\n    value          0.95 fit prior 0.95 3%\n"
        b = "dataset p12c_90\n"
        assert t.count(a) == 1 and t.count(b) == 1
        t = t.replace(a, "normalization B\n    value          1.9537733 fixed\n")
        i = t.index(b); j = t.index("    norm           1 fixed\n", i)
        t = t[:j] + "    norm           0.8 fixed\n" + t[j + len("    norm           1 fixed\n"):]
        azr = os.path.join(tmp, "fixnorm.azr"); open(azr, "w").write(t)
        segs = SegmentTable.from_file(azr)
        B = [s for s in segs if s.norm_group == "B"]; s90 = segs.by_name("p12c_90")
        check("(a) from_file: fixed record B -> 1.9537733 on its data sets", B and all(abs(s.norm - 1.9537733) < 1e-12 for s in B), [s.norm for s in B])
        check("(a) from_file: own 'norm 0.8 fixed' -> 0.8", abs(s90.norm - 0.8) < 1e-12, s90.norm)
        with SingleSession(azr, cwd=tmp, flags=dict(use_brune=True), name="fixnorm") as ev:
            seg = {s.key: s for s in ev.segments}
            check("(b) session segments: B members 1.9537733, p12c_90 0.8",
                  all(abs(seg[s.key].norm - 1.9537733) < 1e-12 for s in B) and abs(seg[s90.key].norm - 0.8) < 1e-12,
                  [(s.key, seg[s.key].norm) for s in B] + [seg[s90.key].norm])
            n, col = ev._cov_norm_rows(ev.x0)
            rr = dict(ev.row_ranges)
            okB = all(np.allclose(n[slice(*rr[s.key])], 1.9537733) and np.all(col[slice(*rr[s.key])] == -1) for s in B if s.key in rr)
            ok90 = s90.key in rr and np.allclose(n[slice(*rr[s90.key])], 0.8)
            check("(b) _cov_norm_rows: fixed rows carry the written values", okB and ok90)
            A = [s for s in segs if s.norm_group == "A" and s.key in rr]
            check("(b) free record A follows its free column", A and all(np.all(col[slice(*rr[s.key])] >= 0) for s in A))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("ALL PASSED" if not failures else f"FAILED: {failures}")
    return 1 if failures else 0
if __name__ == "__main__":
    sys.exit(main())
