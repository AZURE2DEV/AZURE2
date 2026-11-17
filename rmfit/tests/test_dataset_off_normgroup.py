#!/usr/bin/env python3
"""Dataset-off with a normgroup LEADER switched off.

The engine takes a normgroup's norm, vary flag and norm error from the group's first ACTIVE segment, and rmfit keys
the group's fitted norm by that segment (ParamKey("norm", seg=leader)).  Switching the leader off used to hand the
group to the next member with that member's own .azr fields (another penalty centre and width) and to start its norm
from that field instead of the fitted value (in one archived fit the members' fields were 6 % apart and the fitted
norm 30 % below both).

Fixture fixtures/relvector_p12c/ngrel.azr with its `covariance` keywords stripped (so the chi2 is a plain sum over
segments): group A = segments 1 (norm 1, free, 5 %), 3 and 4 (norm 1.1, free, error 0).  Switch p12c_60 (segment 1)
off.  Checks, on Sharded evaluators:
  (a) dataset_off_variant reports moved = {1: 3} and segment 3's norm triple is segment 1's (1, 1, 5);
  (b) the group norm starts at the incumbent's fitted value (0.8 here, distinct from both .azr fields);
  (c) at that start every remaining segment's chi2 and the total penalty equal the base model's at the incumbent
      (removing a dataset changes nothing else when nothing is refitted);
  (d) teeth: the old path (no field copy, take() only) changes both the penalty and segments 3-4's chi2.
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test (in format 3 a normgroup's norm is one shared record; nothing to copy)")
    sys.exit(77)
import numpy as np  # noqa: E402
from _fixtures import fixture  # noqa: E402
from rmfit.engine import Sharded  # noqa: E402
from rmfit.keys import ParamKey, take  # noqa: E402
from rmfit.covariance import strip_text, segment_lines, norm_field  # noqa: E402
from rmfit.structure import dataset_off_variant, carry_moved_norms  # noqa: E402
from pyazr import AzrModel  # noqa: E402

FLAGS = dict(use_brune=True)
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def triple(azr, key):
    t = segment_lines(azr)[key - 1]
    f = norm_field(t)
    return [float(v) for v in t[f:f + 3]]


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_offng_")
    src = fixture("relvector_p12c")
    shutil.copytree(os.path.join(src, "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    base = os.path.join(tmp, "ng.azr")
    open(base, "w").write(strip_text(open(os.path.join(src, "ngrel.azr")).read(), ("covariance",)))
    mdl = AzrModel.from_file(base)
    mdl.set_output_dir(os.path.join(tmp, "output"))
    mdl.write(base)
    check("fixture: group A led by segment 1 (1, 1, 5); segment 3 carries (1.1, 1, 0)",
          triple(base, 1) == [1, 1, 5] and triple(base, 3) == [1.1, 1, 0], f"{triple(base, 1)} {triple(base, 3)}")

    off = os.path.join(tmp, "off.azr")
    info = dataset_off_variant(base, off, "p12c_60", output_dir=os.path.join(tmp, "output_off"))
    os.makedirs(os.path.join(tmp, "output_off"), exist_ok=True)
    check("(a) moved = {1: 3}", info.get("moved") == {1: 3}, f"{info}")
    check("(a) segment 3 now carries the group's norm triple (1, 1, 5)", triple(off, 3) == [1, 1, 5], f"{triple(off, 3)}")
    check("(a) segment 4 (still a follower) is unchanged", triple(off, 4) == [1.1, 1, 0], f"{triple(off, 4)}")

    old = os.path.join(tmp, "off_old.azr")                  # the pre-fix variant: active flag only
    m2 = AzrModel.from_file(base)
    m2.set_segment_active("p12c_60", False)
    m2.set_output_dir(os.path.join(tmp, "output_off"))
    m2.write(old)

    ev = Sharded(base, os.path.join(tmp, "sh_base"), nshards=2, flags=FLAGS)
    ev2 = Sharded(off, os.path.join(tmp, "sh_off"), nshards=2, flags=FLAGS)
    ev3 = Sharded(old, os.path.join(tmp, "sh_old"), nshards=2, flags=FLAGS)
    try:
        x = np.array(ev.x0, float)
        jA = ev.keys.index(ParamKey("norm", seg=1))
        x[jA] = 0.8
        s = ev.score(x)
        x2 = carry_moved_norms(x, ev.keys, take(x, ev.keys, ev2.keys, ev2.x0), ev2.keys, info["moved"])
        j2 = ev2.keys.index(ParamKey("norm", seg=3))
        check("(b) the group norm starts at the fitted 0.8", abs(x2[j2] - 0.8) < 1e-15, f"{x2[j2]}")
        s2 = ev2.score(x2)
        rest = {k: v for k, v in s.segment_chi2.items() if k != 1}
        dseg = max(abs(s2.segment_chi2.get(k, np.nan) - v) / max(abs(v), 1.0) for k, v in rest.items())
        print(f"  base: chi2 {s.chi2:.6f} (seg 1 {s.segment_chi2.get(1, 0):.6f}), penalty {s.penalty:.6f}; "
              f"off: chi2 {s2.chi2:.6f}, penalty {s2.penalty:.6f}; max rel seg diff {dseg:.2e}")
        check("(c) every remaining segment's chi2 is unchanged", dseg < 1e-9, f"{dseg}")
        check("(c) the penalty is unchanged", abs(s2.penalty - s.penalty) < 1e-9 * max(1.0, s.penalty),
              f"{s.penalty} vs {s2.penalty}")
        check("(c) total: off objective = base objective - segment 1's chi2",
              abs(s2.objective - (s.objective - s.segment_chi2[1])) < 1e-8 * s.objective,
              f"{s2.objective} vs {s.objective - s.segment_chi2[1]}")
        x3 = take(x, ev.keys, ev3.keys, ev3.x0)
        s3 = ev3.score(x3)
        j3 = ev3.keys.index(ParamKey("norm", seg=3))
        d3 = max(abs(s3.segment_chi2[k] - s.segment_chi2[k]) for k in (3, 4))
        print(f"  old path: group norm starts at {x3[j3]:.4f}, penalty {s3.penalty:.6f}, max seg 3-4 chi2 change {d3:.3g}")
        check("(d) teeth: the old path starts the group at segment 3's own field, not 0.8", abs(x3[j3] - 0.8) > 0.1,
              f"{x3[j3]}")
        check("(d) teeth: the old path changes the penalty", abs(s3.penalty - s.penalty) > 1e-3, f"{s3.penalty}")
        check("(d) teeth: the old path changes segments 3-4's chi2", d3 > 1e-3, f"{d3}")
    finally:
        for e in (ev, ev2, ev3):
            e.close()
    shutil.rmtree(tmp, ignore_errors=True)
    print("FAILURES: " + ", ".join(failures) if failures else "ALL OK")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
