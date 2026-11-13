#!/usr/bin/env python3
"""A rebuild drops what the session derived from the old model (AZUREAPI::ReplaceModel).

T2 with a <parameterSettings> nuisance prior on a free width (Level 1 Channel 3, 80000 +- 5000 eV, 1 at the start;
the prior is converted to the reduced amplitude at the session's channel radius).  A session evaluates the objective
(which builds the nuisance priors' manager), then set_channel_radius rebuilds the model; its objective, gradient and
nuisance term must equal a fresh session's given the same radius before its first evaluation.  Until 2026-11-13 the
manager survived the rebuild, holding the freed model and the old radius's conversion: nuisance 1.217 instead of 1
at 5 fm.

Needs the compiled module.  Run from anywhere:  python3 tests/pyazr/rebuild_nuisance_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
T2 = os.path.join(ROOT, "tests", "upgrade", "T2")
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
from pyazr import Session  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def main():
    work = tempfile.mkdtemp(prefix="rebuild_nuisance.")
    try:
        d = os.path.join(work, "T2")
        shutil.copytree(T2, d)
        os.makedirs(os.path.join(d, "output"), exist_ok=True)
        p = os.path.join(d, "T2.azr")
        text = open(p).read().replace("<parameterSettings>\n",
                                      "<parameterSettings>\nLevel 1 Channel 3 Width (eV) 80000 0 0 5000 0 1 level 0\n")
        open(p, "w").write(text)
        with Session(p, cwd=d) as m:
            x = np.array(m.parameters.x0, float)
            m.objective().value_and_grad(x)
            check("the nuisance prior is 1 at the start", abs(m.evaluate(x).terms["nuisance priors"] - 1.0) < 1e-9,
                  m.evaluate(x).terms["nuisance priors"])
            m.set_channel_radius(1, 5.0)
            x = np.array(m.parameters.x0, float)
            o_a, g_a = m.objective().value_and_grad(x)
            n_a = m.evaluate(x).terms["nuisance priors"]
        with Session(p, cwd=d) as m:
            m.set_channel_radius(1, 5.0)
            x = np.array(m.parameters.x0, float)
            o_b, g_b = m.objective().value_and_grad(x)
            n_b = m.evaluate(x).terms["nuisance priors"]
        check("after the rebuild: the nuisance term of a fresh session", n_a == n_b, (n_a, n_b))
        check("after the rebuild: the objective of a fresh session", o_a == o_b, (o_a, o_b))
        check("after the rebuild: the gradient of a fresh session", bool(np.array_equal(g_a, g_b)),
              float(np.max(np.abs(g_a - g_b))))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print(f"{len(failures)} failure(s)" if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
