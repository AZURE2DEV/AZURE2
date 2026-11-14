#!/usr/bin/env python3
"""A scan over a channel radius is a Project edit and a new session per radius: the sessions of one process must not
accumulate memory, re-read unchanged input files, or print one another's setup Notes.

Sampling a channel radius in a likelihood (3H+d free-radius MCMC, 2026-11-03) once grew a Coulomb-function memo by one
set of keys per radius, re-read the .azr and every data file per radius (which on /groups stalled 16 forked workers),
and printed every setup "Note:" line again per radius.

Checks, on tests/angle_mapping/angle_mapping.azr (p+12C and p+13C, a shift group and a shift function),
without the Coulomb table (which covers these channels and would leave the memo unused):
  1. 8 sessions at random radii of pair 1, one after the other: no Coulomb memo holds anything (keys 0 after each),
     and a session at a radius after them gives the chi2 of the first session at that radius (rel 1e-12);
  2. every session prints the shift-function and shift-group Notes once;
  3. the input files are cached per process: a session at a new radius reads its new .azr from disk and serves the
     unchanged data files from memory; an edited data file is read again (only it, with the new .azr) and its new
     values are used (chi2 equals a session on the edited file in its own directory).

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
        from pyazr import Session, Project
        import pyazr._azure2  # noqa: F401
    except Exception as e:  # pragma: no cover
        print(f"SKIP: engine not importable ({e})")
        return 77

    work = tempfile.mkdtemp(prefix="pyazr_rebuild_")
    try:
        src = os.path.join(ROOT, "tests", "angle_mapping")

        def project(name, radius=None, edit=None):
            d = os.path.join(work, name)
            shutil.copytree(os.path.join(src, "data"), os.path.join(d, "data"))
            for sub in ("output", "checks"):
                os.makedirs(os.path.join(d, sub), exist_ok=True)
            mdl = Project.open(os.path.join(src, "angle_mapping.azr"))
            if radius is not None:
                mdl.set_pair(1, radius=radius)
            azr = mdl.save(os.path.join(d, "m.azr"))
            if edit:
                edit(d)
            return d, azr

        def chi2(m):
            return m.evaluate(np.asarray(m.parameters.x0, float)).chi2

        R_STAR = 4.37
        d_ref, azr_ref = project("ref", radius=R_STAR)
        with Session(azr_ref, cwd=d_ref, use_coulomb_table=False) as r:
            chi_ref = chi2(r)

        d, azr = project("main")
        # The copies were written just now; InputFiles re-reads a file modified within 3 s of being read (a later
        # same-size write in that window could keep its stamp, 2026-11-05), so date them back an hour, as the data
        # files of a real project are, for the "read once" checks of section 3.
        old = time.time() - 3600.0
        for root, _, files in os.walk(d):
            for name in files:
                os.utime(os.path.join(root, name), (old, old))
        what_fn, what_grp = "shift function on segment", "shift group 'G'"

        def at_radius(radius, name):
            """A session of the main copy with pair 1 at ``radius`` (fm), its .azr saved as ``name``: (chi2, stdout)."""
            mdl = Project.open(azr)
            mdl.set_pair(1, radius=radius)
            path = mdl.save(os.path.join(d, name))
            os.utime(path, (old, old))  # dated back as the rest (a file written within 3 s is read at every open)
            out = []
            with captured_stdout(out):
                with Session(path, cwd=d, use_coulomb_table=False) as m:
                    c = chi2(m)
                    keys = m.cache_stats()["keys"]
            return c, out[-1], keys

        print("1. No Coulomb memo across sessions")
        rng = np.random.default_rng(20261104)
        runs = [at_radius(4.0 + 1.5 * rng.random(), f"r{i}.azr") for i in range(8)]
        check("no Coulomb memo: keys 0 after each of 8 sessions at random radii", all(k == 0 for _, _, k in runs),
              f"{[k for _, _, k in runs]}")
        c, text, _ = at_radius(R_STAR, "rstar.azr")
        check("chi2 at a radius after them == the first session at it", rel(c, chi_ref) <= 1e-12,
              f"{c!r} vs {chi_ref!r} (rel {rel(c, chi_ref):.2e})")

        print("\n2. Notes once per session")
        counts = [(notes(t, what_fn), notes(t, what_grp)) for _, t, _ in runs] + [(notes(text, what_fn),
                                                                                 notes(text, what_grp))]
        check("every session prints the shift-function and shift-group Notes once", all(n == (1, 1) for n in counts),
              f"{counts}")

        print("\n3. input files read once per process")
        with Session(azr, cwd=d, use_coulomb_table=False) as m0:
            st = m0.cache_stats()
        check("cache_stats reports the input files", "input_files_read" in st, f"{sorted(st)}")
        if "input_files_read" in st:
            for k in range(3):
                at_radius(4.5 + 0.1 * k, f"s{k}.azr")
            with Session(azr, cwd=d, use_coulomb_table=False) as m0:
                st2 = m0.cache_stats()
            check("3 sessions at new radii read only their 3 new .azr files from disk",
                  st2["input_files_read"] - st["input_files_read"] == 3,
                  f"{st['input_files_read']} -> {st2['input_files_read']}")
            check("... and serve the rest from memory", st2["input_files_served"] > st["input_files_served"],
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
            with Session(azr_edit, cwd=d_edit, use_coulomb_table=False) as e:
                chi_edit = chi2(e)
            check("the edit changes chi2", rel(chi_edit, chi_ref) > 1e-6, f"{chi_edit} vs {chi_ref}")
            with Session(azr, cwd=d, use_coulomb_table=False) as m0:
                before = m0.cache_stats()["input_files_read"]
            scale_165(d)
            os.utime(os.path.join(d, "data", "p12c_165.dat"), (old, old))
            c, _, _ = at_radius(R_STAR, "edited.azr")
            with Session(azr, cwd=d, use_coulomb_table=False) as m0:
                after = m0.cache_stats()["input_files_read"]
            check("the edited file (only) is read again, with the new .azr", after - before == 2,
                  f"{before} -> {after}")
            check("its new values are used", rel(c, chi_edit) <= 1e-12, f"{c!r} vs {chi_edit!r}")
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
