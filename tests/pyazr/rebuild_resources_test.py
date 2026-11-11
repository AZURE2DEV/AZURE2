#!/usr/bin/env python3
"""A pyazr rebuild (set_channel_radius, a potential change) must not accumulate memory, re-read unchanged
input files, or print the setup Notes again.

Before 2026-11-04, sampling a channel radius in a likelihood (3H+d free-radius MCMC, 2026-11-03):
  - grew the Coulomb-function memo by one set of keys per new radius, never released (the key holds the
    radius): ~96 kB per call, 239 GB over 16 workers;
  - re-read the .azr and every data file at every rebuild, which on /groups stalled 16 forked workers;
  - printed every setup "Note:" line (shift function, shift group, ...) again at every rebuild.

Checks, on tests/angle_mapping/angle_mapping.azr (p+12C and p+13C, a shift group and a shift function),
without the Coulomb table (which covers these channels and would leave the memo unused):
  1. 40 rebuilds at random radii of pair 1: the memo holds no more keys than after one calculation at one
     radius, and chi2 at a radius equals a fresh session on an .azr with that radius written in (rel 1e-12);
  2. the Notes of the session's first build are not printed again by rebuilds, and a new session prints them;
  3. rebuilds read no input file from disk; an edited data file is read again (only it) and its new values
     are used (chi2 equals a fresh session on the edited file).
The old engine fails 1, 2 and 3.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/rebuild_resources_test.py
"""
import contextlib
import os
import time
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


@contextlib.contextmanager
def captured_stdout(sink):
    """Collects what the engine (C++ std::cout) writes to file descriptor 1."""
    sys.stdout.flush()
    saved = os.dup(1)
    with tempfile.TemporaryFile(mode="w+b") as f:
        os.dup2(f.fileno(), 1)
        try:
            yield
        finally:
            sys.stdout.flush()
            os.dup2(saved, 1)
            os.close(saved)
            f.seek(0)
            sink.append(f.read().decode(errors="replace"))


def notes(text, what):
    return sum(1 for line in text.splitlines() if line.startswith("Note:") and what in line)


def main():
    sys.path.insert(0, ROOT)
    try:
        import numpy as np
        from pyazr import azure2, AzrModel
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"SKIP: engine not importable ({e})")
        return 0

    work = tempfile.mkdtemp(prefix="pyazr_rebuild_")
    try:
        src = os.path.join(ROOT, "tests", "angle_mapping")

        def project(name, radius=None, edit=None):
            d = os.path.join(work, name)
            shutil.copytree(os.path.join(src, "data"), os.path.join(d, "data"))
            for sub in ("output", "checks"):
                os.makedirs(os.path.join(d, sub), exist_ok=True)
            mdl = AzrModel.from_file(os.path.join(src, "angle_mapping.azr"))
            if radius is not None:
                mdl.set_channel_radius(1, radius)
            azr = mdl.write(os.path.join(d, "m.azr"))
            if edit:
                edit(d)
            return d, azr

        def chi2(m):
            return float(np.sum(m.calculate_chi2_rwa(np.asarray(m.params_rwa, float))))

        R_STAR = 4.37
        d_ref, azr_ref = project("ref", radius=R_STAR)
        with azure2(azr_ref, cwd=d_ref, use_coulomb_table=False) as r:
            chi_ref = chi2(r)

        d, azr = project("main")
        # The copies were written just now; InputFiles re-reads a file modified within 3 s of being read (a later
        # same-size write in that window could keep its stamp, 2026-11-05), so date them back an hour, as the data
        # files of a real project are, for the "read once" checks of section 3.
        old = time.time() - 3600.0
        for root, _, files in os.walk(d):
            for name in files:
                os.utime(os.path.join(root, name), (old, old))
        out = []
        with captured_stdout(out):
            m = azure2(azr, cwd=d, use_coulomb_table=False)
        first = out[-1]
        try:
            what_fn, what_grp = "shift function on segment", "shift group 'G'"
            check("the first build prints the shift-function and shift-group Notes",
                  notes(first, what_fn) == 1 and notes(first, what_grp) == 1,
                  f"{notes(first, what_fn)} / {notes(first, what_grp)}")

            print("1. Coulomb memo across rebuilds")
            chi2(m)
            k0 = m.cache_stats()["keys"]
            rng = np.random.default_rng(20261104)
            out = []
            with captured_stdout(out):
                for _ in range(40):
                    m.set_channel_radius(1, 4.0 + 1.5 * rng.random())
                    chi2(m)
            k40 = m.cache_stats()["keys"]
            check("keys after 40 random radii <= keys after one radius", 0 < k40 <= k0, f"{k0} -> {k40}")
            m.set_channel_radius(1, R_STAR)
            c = chi2(m)
            check("chi2 at a radius == fresh session with it in the file", rel(c, chi_ref) <= 1e-12,
                  f"{c!r} vs {chi_ref!r} (rel {rel(c, chi_ref):.2e})")

            print("\n2. Notes once per session")
            rebuilds = out[-1]
            check("40 rebuilds print no shift-function Note", notes(rebuilds, what_fn) == 0,
                  f"{notes(rebuilds, what_fn)} lines")
            check("40 rebuilds print no shift-group Note", notes(rebuilds, what_grp) == 0,
                  f"{notes(rebuilds, what_grp)} lines")
            out = []
            with captured_stdout(out):
                with azure2(azr, cwd=d) as m2:
                    chi2(m2)
            check("a new session prints them again", notes(out[-1], what_fn) == 1 and notes(out[-1], what_grp) == 1,
                  f"{notes(out[-1], what_fn)} / {notes(out[-1], what_grp)}")

            print("\n3. input files read once")
            st = m.cache_stats()
            check("cache_stats reports the input files", "input_files_read" in st, f"{sorted(st)}")
            if "input_files_read" in st:
                for k in range(5):
                    m.set_channel_radius(1, 4.5 + 0.1 * k)
                st2 = m.cache_stats()
                check("5 rebuilds read nothing from disk", st2["input_files_read"] == st["input_files_read"],
                      f"{st['input_files_read']} -> {st2['input_files_read']}")
                check("... and are served from memory", st2["input_files_served"] > st["input_files_served"],
                      f"{st['input_files_served']} -> {st2['input_files_served']}")

                def scale_165(dd):
                    p = os.path.join(dd, "data", "p12c_165.dat")
                    rows = []
                    for line in open(p):
                        f = line.split()
                        if len(f) >= 4 and not line.lstrip().startswith("#"):
                            f[2] = repr(float(f[2]) * 1.07)
                            line = " ".join(f) + "\n"
                        rows.append(line)
                    with open(p, "w") as fh:
                        fh.writelines(rows)

                d_edit, azr_edit = project("edited", radius=R_STAR, edit=scale_165)
                with azure2(azr_edit, cwd=d_edit, use_coulomb_table=False) as e:
                    chi_edit = chi2(e)
                check("the edit changes chi2", rel(chi_edit, chi_ref) > 1e-6, f"{chi_edit} vs {chi_ref}")
                before = m.cache_stats()["input_files_read"]
                scale_165(d)
                m.set_channel_radius(1, R_STAR)
                after = m.cache_stats()["input_files_read"]
                check("the edited file (only) is read again", after - before == 1, f"{before} -> {after}")
                c = chi2(m)
                check("its new values are used", rel(c, chi_edit) <= 1e-12, f"{c!r} vs {chi_edit!r}")
        finally:
            m.close() if hasattr(m, "close") else m.__exit__(None, None, None)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
