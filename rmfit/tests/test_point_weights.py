#!/usr/bin/env python3
"""ObjectiveSpec.point_weights (2026-09-28, the evaluator: weight the chi2 to the resonance peaks while a narrow
level is fitted).  Rows whose excitation energy lies in a window get their chi2 term multiplied by w;
penalty rows never; overlapping windows multiply; combines with the reduced-chi2 weighting; the spec
is no longer "plain".  Stubbed evaluator, no engine.  Run:  python3 rmfit/tests/test_point_weights.py"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.engine import Evaluator, ObjectiveSpec  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


class Seg:
    def __init__(self, key, f):
        self.key, self.data_file = key, f


class Stub(Evaluator):
    def __init__(self, spec):
        self.row_ranges = [(1, (0, 4)), (2, (4, 8))]
        self.npts = 8
        self.segments = [Seg(1, "a.dat"), Seg(2, "b.dat")]
        self.ex = np.array([11.0, 11.1, 11.2, 11.3, 11.15, 11.25, 11.35, 11.45])
        self._init_objective(spec)

    def row_excitation(self):
        return self.ex


def main():
    ev = Stub(ObjectiveSpec())
    check("no windows: rows unweighted", ev.row_sw is None)
    check("no windows: plain", ev.spec.plain)

    ev = Stub(ObjectiveSpec(point_weights=[[11.1, 11.3, 10.0]]))
    w = ev.row_sw ** 2
    check("rows inside the window x10", np.allclose(w, [1, 10, 10, 10, 10, 10, 1, 1]), w)
    check("weighted spec is not plain", not ev.spec.plain)
    check("label names the window", "x10 on Ex 11.1-11.3" in ev.spec.label, ev.spec.label)

    ev = Stub(ObjectiveSpec(point_weights=[[11.1, 11.3, 10.0], [11.25, 11.5, 2.0]]))
    w = ev.row_sw ** 2
    check("overlapping windows multiply", np.allclose(w, [1, 10, 10, 20, 10, 20, 2, 2]), w)

    ev = Stub(ObjectiveSpec(dataset_weight="reduced", point_weights=[[11.0, 11.0, 4.0]]))
    w = ev.row_sw ** 2
    check("combines with reduced weighting (equal sets -> 1 each, then x4 on one row)", np.allclose(w, [4, 1, 1, 1, 1, 1, 1, 1]), w)

    try:
        Stub(ObjectiveSpec(point_weights=[[11.0, 11.5, 0.0]])); bad = False
    except ValueError:
        bad = True
    check("non-positive weight rejected", bad)

    # the weights enter residual rows, but plain per-row chi2 is recovered by _plain_data_rows
    ev = Stub(ObjectiveSpec(point_weights=[[11.1, 11.3, 9.0]]))
    r = np.ones(8)
    rw = r * ev.row_sw
    check("plain rows recovered from weighted rows", np.allclose(ev._plain_data_rows(rw), r))

    # dataset_weight "scaled" (2026-09-28): rows of a data set weighted 1/max(chi2/N, 1) at the baseline
    ev = Stub(ObjectiveSpec(dataset_weight="scaled", scales={"a.dat": 25.0, "b.dat": 0.4}))
    w = ev.row_sw ** 2
    check("scaled: set a (chi2/N 25) weighted 1/25, set b (0.4) floored at 1", np.allclose(w, [0.04] * 4 + [1.0] * 4), w)
    check("scaled: not plain", not ev.spec.plain)
    check("scaled: label", "errors scaled" in ev.spec.label, ev.spec.label)
    try:
        Stub(ObjectiveSpec(dataset_weight="scaled")); bad = False
    except ValueError:
        bad = True
    check("scaled without scales rejected", bad)
    ev = Stub(ObjectiveSpec(dataset_weight="scaled", scales={"a.dat": 4.0, "b.dat": 1.0}, point_weights=[[11.1, 11.1, 10.0]]))
    w = ev.row_sw ** 2
    check("scaled combines with point weights", np.allclose(w, [0.25, 2.5, 0.25, 0.25, 1, 1, 1, 1]), w)

    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all point-weight checks passed")


if __name__ == "__main__":
    main()
