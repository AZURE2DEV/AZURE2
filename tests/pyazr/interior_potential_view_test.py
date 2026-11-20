#!/usr/bin/env python3
"""pyazr's view of a potential that also acts inside the channel radius (docs/upgrade/interior_potential.md): the
record words `interior_boundary`, `l` and a signed depth reach Session.nuclear_potential.

On tests/upgrade/interior_bound (alpha + 12C, one record with interior_boundary 0):
  (1) nuclear_potential(1) carries interior_boundary 0.0 and V0 45.891; with l=1 the same record (the pair has no
      record of single l);
  (2) Project.set(a12C, interior_boundary="-1") gives -1.0 in a new session;
  (3) a second record with `l 1`, a negative depth and no interior_boundary: nuclear_potential(1, l=1) is that record
      (V0 < 0, interior_boundary None), l=0 and the pair's own view stay the first record;
  (4) Project.set(a12C, interior_charge_e1="0.012") gives interior_charge {"E1": 0.012, "E2": 0.0} in a new session
      (0 and 0 before).
Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/interior_potential_view_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "upgrade", "interior_bound", "interior_bound.azr")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    from pyazr import Project, Session
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

work = tempfile.mkdtemp(prefix="interior_view_")
try:
    for d in ("output", "checks"):
        os.makedirs(os.path.join(work, d))
    base = os.path.join(work, "interior_bound.azr")
    shutil.copy(SRC, base)
    p = Session(base, cwd=work, extrapolations=True).nuclear_potential(1)
    check("(1) the pair's view: interior_boundary 0.0, V0 45.891", p.interior_boundary == 0.0 and p.V0 == 45.891 and p.l is None,
          repr(p))
    q = Session(base, cwd=work, extrapolations=True).nuclear_potential(1, l=1)
    check("(1) l=1 without a record of single l: the pair's record", q.interior_boundary == 0.0 and q.V0 == 45.891 and q.l == 1,
          repr(q))

    edited = Project.open(base)
    edited.set("a12C", interior_boundary="-1")
    f2 = os.path.join(work, "edited.azr")
    edited.save(f2)
    p2 = Session(f2, cwd=work, extrapolations=True).nuclear_potential(1)
    check("(2) Project.set(a12C, interior_boundary=-1), a new session: -1.0", p2.interior_boundary == -1.0, repr(p2))

    text = open(base).read()
    record = ("potential a12C_p\n    pair            a+12C\n    l               1\n    coulomb_radius  3.55 fm\n"
              "    central         woods_saxon depth -5 MeV radius 3.55 fm diffuseness 0.6 fm\n")
    assert "</potentials>" in text
    f3 = os.path.join(work, "by_l.azr")
    with open(f3, "w") as f:
        f.write(text.replace("</potentials>", record + "</potentials>"))
    S = Session(f3, cwd=work, extrapolations=True)
    l1, l0, own = S.nuclear_potential(1, l=1), S.nuclear_potential(1, l=0), S.nuclear_potential(1)
    check("(3) l=1: the record of l = 1, a repulsive potential without interior_boundary",
          l1.V0 == -5.0 and l1.interior_boundary is None and l1.own, repr(l1))
    check("(3) l=0 and the pair's view: the first record", l0.V0 == 45.891 and l0.interior_boundary == 0.0 and own.V0 == 45.891,
          f"{l0!r} {own!r}")

    charged = Project.open(base)
    charged.set("a12C", interior_charge_e1="0.012")
    f4 = os.path.join(work, "charged.azr")
    charged.save(f4)
    p4 = Session(f4, cwd=work, extrapolations=True).nuclear_potential(1)
    check("(4) Project.set(a12C, interior_charge_e1=0.012), a new session: E1 0.012, E2 0; 0 and 0 before",
          p4.interior_charge == {"E1": 0.012, "E2": 0.0} and p.interior_charge == {"E1": 0.0, "E2": 0.0},
          f"{p4.interior_charge!r} {p.interior_charge!r}")
finally:
    shutil.rmtree(work, ignore_errors=True)

print("FAILED: " + ", ".join(failures) if failures else "all checks passed")
sys.exit(1 if failures else 0)
