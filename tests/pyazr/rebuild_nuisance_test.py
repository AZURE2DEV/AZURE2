#!/usr/bin/env python3
"""The nuisance priors of a session belong to its model: a session of a Project edit made after another session has
evaluated its objective in the same process gives the numbers of a session made first.

T2 with a <parameterSettings> nuisance prior on a free width (Level 1 Channel 3, 80000 +- 5000 eV, 1 at the start;
the prior is converted to the reduced amplitude at the session's channel radius).  A session at the file radius
evaluates the objective (which builds the nuisance priors' manager); then a session of the project with pair 1 at
5.0 fm (Project.set_pair) gives its objective, gradient and nuisance term, and a second session of that edited project
must give the same, bit for bit, with the nuisance term 1 at the start (a manager that kept the other model's
conversion gave 1.217 at 5 fm).

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
from pyazr import Project, Session  # noqa: E402

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
        edited = Project.open(p)
        edited.set_pair(1, radius=5.0)
        p5 = edited.save(os.path.join(d, "T2_5fm.azr"))
        with Session(p, cwd=d) as m:
            x = np.array(m.parameters.x0, float)
            m.objective().value_and_grad(x)
            check("the nuisance prior is 1 at the start", abs(m.evaluate(x).terms["nuisance priors"] - 1.0) < 1e-9,
                  m.evaluate(x).terms["nuisance priors"])
        results = []
        for _ in range(2):
            with Session(p5, cwd=d) as m:
                x = np.array(m.parameters.x0, float)
                o, g = m.objective().value_and_grad(x)
                results.append((o, g, m.evaluate(x).terms["nuisance priors"]))
        (o_a, g_a, n_a), (o_b, g_b, n_b) = results
        check("at 5.0 fm the nuisance prior is 1 at the start", abs(n_a - 1.0) < 1e-9, n_a)
        check("the edited project's nuisance term, twice", n_a == n_b, (n_a, n_b))
        check("its objective, twice", o_a == o_b, (o_a, o_b))
        check("its gradient, twice", bool(np.array_equal(g_a, g_b)),
              float(np.max(np.abs(g_a - g_b))))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print(f"{len(failures)} failure(s)" if failures else "PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
