#!/usr/bin/env python3
"""A channel radius is an edit of the Project and a new session: the radii written by Project.set_pair are the ones the
session computes with, and a scan over radii in one process gives every radius its own model.

Checks, on tests/channel_space/channel_space.azr (p+12C and p+12C*, both 5.0 fm, plus a gamma pair):
  1. Project.set_pair(1, radius=4.2) and (2, 6.3) in one edit: the session reads both radii, and the radii matter
     (chi2 differs from the file's by more than 1e-3);
  2. the same edit made in the other order, and as two saves of one project, gives the same reduced widths
     (rel 1e-12) and chi2 (rel 1e-10) as (1);
  3. a scan in one process (4.2, 6.3 on pair 2, then the file radius back): each session has the radius it was
     saved with, and the session at the file radius gives the file's reduced widths and chi2 (rel 1e-12, 1e-10);
  4. the session has no method that changes a channel radius (set_channel_radius is gone).

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/radius_override_test.py
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
        from pyazr import Session, Project
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"SKIP: engine not importable ({e})")
        return 77

    R1, R2 = 4.2, 6.3
    work = tempfile.mkdtemp(prefix="pyazr_radius_override_")
    try:
        src = os.path.join(ROOT, "tests", "channel_space")
        shutil.copy(os.path.join(src, "channel_space.azr"), work)
        shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
        for d in ("output", "checks"):
            os.makedirs(os.path.join(work, d), exist_ok=True)
        base = os.path.join(work, "channel_space.azr")

        # reference: both radii written into the file, fresh session, own output directory
        ref_dir = os.path.join(work, "ref")
        for d in ("output", "checks"):
            os.makedirs(os.path.join(ref_dir, d), exist_ok=True)
        shutil.copytree(os.path.join(src, "data"), os.path.join(ref_dir, "data"))
        mdl = Project.open(base)
        mdl.set_pair(1, radius=R1)
        mdl.set_pair(2, radius=R2)
        ref_azr = os.path.join(ref_dir, "both.azr")
        mdl.save(ref_azr)
        with Session(ref_azr, cwd=ref_dir) as r:
            check("reference file has both radii", [p.channel_radius for p in r.pairs][:2] == [R1, R2],
                  f"{[p.channel_radius for p in r.pairs]}")
            x_ref = np.asarray(r.parameters.x0, float)
            chi_ref = r.evaluate(x_ref).chi2
        with Session(base, cwd=work) as f:
            x_file = np.asarray(f.parameters.x0, float)
            chi_file = f.evaluate(x_file).chi2
        check("the radii matter (reference chi2 differs from the file's)", rel(chi_ref, chi_file) > 1e-3,
              f"{chi_ref} vs {chi_file}")

        # The .azr holds physical widths, so the session's reduced widths are computed from the radii in force:
        # compare them with the reference file's, and chi2 at them.
        def session_of(edits, name):
            p = Project.open(base)
            for pair, radius in edits:
                p.set_pair(pair, radius=radius)
            return Session(p.save(os.path.join(work, name)), cwd=work)

        def same_vector(x, ref):
            return x.shape == ref.shape and np.allclose(x, ref, rtol=1e-12, atol=0.0)

        print("1-2. both radii through Project edits, in either order")
        for edits, name in (([(1, R1), (2, R2)], "a.azr"), ([(2, R2), (1, R1)], "b.azr")):
            with session_of(edits, name) as m:
                x = np.asarray(m.parameters.x0, float)
                check(f"{name}: both radii in force", [p.channel_radius for p in m.pairs][:2] == [R1, R2],
                      f"{[p.channel_radius for p in m.pairs]}")
                check(f"{name}: reduced widths == the reference file's", same_vector(x, x_ref))
                c = m.evaluate(x).chi2
                check(f"{name}: chi2 == the reference file's", rel(c, chi_ref) <= 1e-10, f"{c} vs {chi_ref}")
        p = Project.open(base)
        p.set_pair(1, radius=R1)
        p.save(os.path.join(work, "c.azr"))
        p.set_pair(2, radius=R2)
        with Session(p.save(os.path.join(work, "c.azr")), cwd=work) as m:
            check("two saves of one project: both radii in force",
                  [q.channel_radius for q in m.pairs][:2] == [R1, R2], f"{[q.channel_radius for q in m.pairs]}")

        print("\n3. a scan in one process")
        scan = (([(1, R1)], [R1, 5.0]), ([(2, R2)], [5.0, R2]), ([(1, 5.0), (2, 5.0)], [5.0, 5.0]))
        for i, (edits, want) in enumerate(scan):
            with session_of(edits, f"scan{i}.azr") as m:
                got = [q.channel_radius for q in m.pairs][:2]
                check(f"radii {want}", got == want, f"{got}")
                if want == [5.0, 5.0]:
                    x = np.asarray(m.parameters.x0, float)
                    check("the file radius back: the file's reduced widths", same_vector(x, x_file))
                    c = m.evaluate(x).chi2
                    check("the file radius back: the file's chi2", rel(c, chi_file) <= 1e-10, f"{c} vs {chi_file}")

        print("\n4. no live radius edit")
        with Session(base, cwd=work) as m:
            check("Session has no set_channel_radius", not hasattr(m, "set_channel_radius"))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
