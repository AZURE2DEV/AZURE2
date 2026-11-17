#!/usr/bin/env python3
"""A channel's own radius is a `channel_radius` record of the Project and a new session, as a pair's radius is a
Project.set_pair and a new session (radius_override_test.py): the session has no method that changes a radius.

Checks, on tests/channel_space/channel_space.azr (p+12C and p+12C*, both 5.0 fm, plus a gamma pair):
  1. Project.add("channel_radius", ...) for every l of pair 1 at 4.2 fm gives the session of Project.set_pair(1,
     radius=4.2): the same reduced widths (rel 1e-12) and chi2 (rel 1e-10), while the pair still reads 5.0 fm;
  2. a record for l = 1 alone gives a model that is neither the file's nor the pair edit's (chi2 differs by more
     than 1e-6 from both);
  3. channel_functions(pair, E, l) with no radius is evaluated at the channel's own radius for the l of the record
     and at the pair's for another l;
  4. records with s that give the channels of one l of pair 2 two radii: channel_functions with no radius is refused, with a
     radius it answers;
  5. Project.set(record, radius=) and Project.remove(record) are the edit and the way back: after the removal the
     session is the file's (rel 1e-12, 1e-10);
  6. a second record of the same (pair, l) is refused by the edit (E_CHANNEL_RADIUS_TWICE).

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/channel_radius_record_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def main():
    sys.path.insert(0, ROOT)
    try:
        import numpy as np
        from pyazr import Session, Project, FormatError
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"SKIP: engine not importable ({e})")
        return 77

    R = 4.2
    E = np.array([0.5, 1.0, 2.0])
    work = tempfile.mkdtemp(prefix="pyazr_channel_radius_record_")
    try:
        src = os.path.join(ROOT, "tests", "channel_space")
        shutil.copy(os.path.join(src, "channel_space.azr"), work)
        shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
        for d in ("output", "checks"):
            os.makedirs(os.path.join(work, d), exist_ok=True)
        base = os.path.join(work, "channel_space.azr")

        def result(path):
            with Session(path, cwd=work) as s:
                x = np.asarray(s.parameters.x0, float)
                return x, s.evaluate(x).chi2, [p.channel_radius for p in s.pairs][:2]

        def same(a, b):
            return a[0].shape == b[0].shape and np.allclose(a[0], b[0], rtol=1e-12, atol=0.0) and rel(a[1], b[1]) <= 1e-10

        file_result = result(base)
        with Session(base, cwd=work) as s:
            ls = sorted({q.L for q in s.parameters if q.kind == "width" and q.pair == 1})
        p = Project.open(base)
        p.set_pair(1, radius=R)
        pair_edit = result(p.save(os.path.join(work, "pair.azr")))

        def with_records(name, records):
            q = Project.open(base)
            for rec, fields in records:
                q.add("channel_radius", rec, **fields)
            return q, q.save(os.path.join(work, name))

        def entry(l, radius, **more):
            return dict(pair="p+12C", l=l, radius=f"{radius} fm", **more)

        print(f"1. a record for every l of pair 1 ({ls}) at {R} fm == Project.set_pair(1, radius={R})")
        _, path = with_records("every_l.azr", [(f"r{l}", entry(l, R)) for l in ls])
        every = result(path)
        check("the pair's radius is still the file's", every[2] == file_result[2], f"{every[2]}")
        check("reduced widths and chi2 == the pair edit's", same(every, pair_edit), f"{every[1]} vs {pair_edit[1]}")
        check("the radius matters (chi2 differs from the file's)", rel(every[1], file_result[1]) > 1e-3,
              f"{every[1]} vs {file_result[1]}")

        print("\n2. a record for l = 1 alone")
        project, one_path = with_records("l1.azr", [("r1", entry(1, R))])
        one = result(one_path)
        check("not the file's model", rel(one[1], file_result[1]) > 1e-6, f"{one[1]} vs {file_result[1]}")
        check("not the pair edit's model", rel(one[1], pair_edit[1]) > 1e-6, f"{one[1]} vs {pair_edit[1]}")

        print("\n3. channel_functions at the channel's own radius")
        with Session(one_path, cwd=work) as s:
            for l, want in ((1, R), (0, 5.0)):
                got, at = s.channel_functions(1, E, l=l), s.channel_functions(1, E, l=l, radius=want)
                check(f"l = {l}: no radius == radius {want} fm",
                      all(np.array_equal(got[k], at[k]) for k in ("P", "S", "delta_hs")))
            other = s.channel_functions(1, E, l=1, radius=5.0)
            check("l = 1: the pair's radius gives another penetrability",
                  not np.allclose(other["P"], s.channel_functions(1, E, l=1)["P"], rtol=1e-6))

        print("\n4. two radii for one l (pair 2, l = 1: records with s = 3/2 and s = 5/2)")
        q = Project.open(base)
        q.add("channel_radius", "ra", pair="p+12C[0.8]", l=1, radius=f"{R} fm", s="3/2")
        q.add("channel_radius", "rb", pair="p+12C[0.8]", l=1, radius="4.6 fm", s="5/2")
        with Session(q.save(os.path.join(work, "two.azr")), cwd=work) as s:
            try:
                s.channel_functions(2, E, l=1)
                check("no radius is refused", False, "no exception")
            except ValueError as e:
                check("no radius is refused", "more than one radius" in str(e), str(e))
            check("a radius answers", bool(np.all(s.channel_functions(2, E, l=1, radius=R)["P"] > 0)))
            got, at = s.channel_functions(2, E, l=0), s.channel_functions(2, E, l=0, radius=5.0)
            check("another l keeps the pair's radius", np.array_equal(got["P"], at["P"]))

        print("\n5. the edit and the way back")
        project.set("r1", radius="4.6 fm")
        moved = result(project.save(os.path.join(work, "moved.azr")))
        check("Project.set(radius=): another model", rel(moved[1], one[1]) > 1e-6, f"{moved[1]} vs {one[1]}")
        project.remove("r1")
        back = result(project.save(os.path.join(work, "back.azr")))
        check("Project.remove: the file's reduced widths and chi2", same(back, file_result), f"{back[1]} vs {file_result[1]}")

        print("\n6. one record per (pair, l, s, J)")
        q, _ = with_records("twice.azr", [("r1", entry(1, R))])
        try:
            q.add("channel_radius", "r1b", **entry(1, 4.6))
            check("a second record of the selection is refused", False, "no exception")
        except FormatError as e:
            check("a second record of the selection is refused", "E_CHANNEL_RADIUS_TWICE" in str(e), str(e))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
