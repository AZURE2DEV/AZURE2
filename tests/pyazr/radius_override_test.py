#!/usr/bin/env python3
"""A channel radius set with azure2.set_channel_radius must survive every later rebuild of the session.

Before 2026-11-02 the API's SetRadius rebuilt the model from the .azr with only its own call's
override, and Rebuild()/Initialize() with none: setting pair 2 put pair 1 back to its file radius,
and so did a nuclear-potential change or a data/extrapolation mode switch -- silently, while the
caller believed the earlier radius was in force (3H+d radius scan, 2026-11-02).  The wrapper
returned the read-back radius without comparing it with the request.

Checks, on tests/channel_space/channel_space.azr (p+12C and p+12C*, both 5.0 fm, plus a gamma pair):
  1. set_channel_radius(1, 4.2) then (2, 6.3): both radii are in force, and chi2 equals a fresh
     session on an .azr with both radii written by AzrModel.set_channel_radius;
  2. a nuclear-potential change (clear_nuclear_potential, which rebuilds) keeps both;
  3. a mode switch (data_mode, which re-initializes) keeps both;
  4. setting a pair back to its file radius is an override like any other (chi2 as the file);
  5. the wrapper raises RuntimeError when the rebuilt model does not have the radius asked for.
The old engine fails 1-3 (pair 1 back at 5.0 fm), the old wrapper 5.

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
        from pyazr import azure2, AzrModel
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"SKIP: engine not importable ({e})")
        return 0

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
        mdl = AzrModel.from_file(base)
        mdl.set_channel_radius(1, R1)
        mdl.set_channel_radius(2, R2)
        ref_azr = mdl.write(os.path.join(ref_dir, "both.azr"))
        with azure2(ref_azr, cwd=ref_dir) as r:
            check("reference file has both radii", [p.channel_radius for p in r.pairs][:2] == [R1, R2],
                  f"{[p.channel_radius for p in r.pairs]}")
            x_ref = np.asarray(r.params_rwa, float)
            chi_ref = float(np.sum(r.calculate_chi2_rwa(x_ref)))
        with azure2(base, cwd=work) as f:
            x_file = np.asarray(f.params_rwa, float)
            chi_file = float(np.sum(f.calculate_chi2_rwa(x_file)))
        check("the radii matter (reference chi2 differs from the file's)", rel(chi_ref, chi_file) > 1e-3,
              f"{chi_ref} vs {chi_file}")

        # The .azr holds physical widths, so the session's reduced widths are recomputed at each rebuild
        # from the radii in force: compare them with the reference file's, and chi2 at them.
        with azure2(base, cwd=work) as m:

            def radii():
                return [p.channel_radius for p in m.pairs][:2]

            def chi():
                return float(np.sum(m.calculate_chi2_rwa(np.asarray(m.params_rwa, float))))

            def same_vector(ref):
                x = np.asarray(m.params_rwa, float)
                return x.shape == ref.shape and np.allclose(x, ref, rtol=1e-12, atol=0.0)

            print("1. two set_channel_radius calls in one session")
            got1 = m.set_channel_radius(1, R1)
            check("pair 1 reads back its radius", got1 == R1, f"{got1}")
            m.set_channel_radius(2, R2)
            check("pair 1 kept after setting pair 2", radii() == [R1, R2], f"{radii()}")
            check("reduced widths == the reference file's", same_vector(x_ref))
            c = chi()
            check("chi2 == fresh session on the file with both radii", rel(c, chi_ref) <= 1e-10,
                  f"{c} vs {chi_ref}")

            print("\n2. a nuclear-potential change (rebuild) keeps them")
            m.clear_nuclear_potential(1)
            check("radii after clear_nuclear_potential", radii() == [R1, R2], f"{radii()}")
            c = chi()
            check("chi2 after clear_nuclear_potential", rel(c, chi_ref) <= 1e-10, f"{c} vs {chi_ref}")

            print("\n3. a mode switch (re-initialize) keeps them")
            m.data_mode()
            m.configure()
            check("radii after data_mode", radii() == [R1, R2], f"{radii()}")
            c = chi()
            check("chi2 after data_mode", rel(c, chi_ref) <= 1e-10, f"{c} vs {chi_ref}")

            print("\n4. back to the file radii")
            m.set_channel_radius(1, 5.0)
            m.set_channel_radius(2, 5.0)
            check("radii back at 5.0", radii() == [5.0, 5.0], f"{radii()}")
            check("reduced widths == the file's", same_vector(x_file))
            c = chi()
            check("chi2 == the file's", rel(c, chi_file) <= 1e-10, f"{c} vs {chi_file}")

            print("\n5. a rebuild that does not take the radius is refused")

            class NoRadius:
                """The session, except that set_radius reports success and changes nothing."""

                def __init__(self, sess):
                    self._sess = sess

                def set_radius(self, pair, radius):
                    return True

                def __getattr__(self, name):
                    return getattr(self._sess, name)

            real = m._sess
            m._sess = NoRadius(real)
            try:
                m.set_channel_radius(1, 3.3)
                check("RuntimeError on a read-back mismatch", False, "no exception")
            except RuntimeError as e:
                check("RuntimeError on a read-back mismatch", "3.3" in str(e), str(e))
            finally:
                m._sess = real
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
