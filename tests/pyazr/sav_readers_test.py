#!/usr/bin/env python3
"""pyazr reads parameter files by NAME, never by row position.

AZURE2 matches `param.sav` rows to parameters by name (AZUREParams::
ReadUserParameters), appends new parameters as the code grows (the
energy-shift sqrt rows, `<data set>.shift.sqrt` in format 3) and now writes a
`#parametrization` tag line first (2026-10-06/07).  A positional reader is
off by one or more rows on any file from another layout, with no error --
an 8Be+alpha rate came out exactly 2x that way.  Pinned here on
tests/identical_pp_res:

  * an old-layout file (no .shift.sqrt rows, no tag, one width nudged): the nudged
    width is read (Parameters.read, Session(parameter_file=)), every missing row
    keeps the session's value, and bands.best_fit_params returns the file's vector of the free-parameter length
    and leaves the session as it was (a session made with the file gives its own x0); on tests/upgrade/T2, a file
    that moves a fixed width is refused by best_fit_params with the call that takes it (Session(file,
    parameter_file=));
  * Parameters.write writes the tag line first and every parameter by name;
  * a new-layout file (tag first) reads every value; a `parametrization` row
    without the '#' (files of 2026-10-06) is skipped the same way;
  * a file with no name in common with the model is refused.

Needs numpy and the compiled engine; exits 77 (skipped) without them.
Run from anywhere:  python3 tests/pyazr/sav_readers_test.py

This test reads the binding itself (`_azure2.Session`), so it uses the private handle `Session._engine`.
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


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    import pyazr
    from pyazr import Session
    from pyazr.bands import best_fit_params
    from pyazr.parameters import read_parameter_file
except Exception as err:
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)

with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "identical_pp_res")
    shutil.copytree(os.path.join(ROOT, "tests", "identical_pp_res"), work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
        os.makedirs(os.path.join(work, junk))
    azr = os.path.join(work, "identical_pp_res.azr")
    sav = os.path.join(work, "output", "param.sav")

    with Session(azr) as m:
        P = m.parameters
        names = P.names
        allrwa = np.asarray(m._engine.params_all_rwa(), float)
        x0 = P.x0
        nfree = x0.size
        widths = P.of_kind("width").free
        target = widths[0]
        nudged = allrwa[target.index] * 1.05

        # --- old layout: no .shift.sqrt rows, no tag, one width nudged
        with open(sav, "w") as fh:
            for p in P:
                if p.name.endswith(".shift.sqrt"):
                    continue
                v = nudged if p.index == target.index else allrwa[p.index]
                fh.write(f"{p.name:>28s} {v: .7e} {0.0: .7e}\n")
        nsqrt = sum(1 for n in names if n.endswith(".shift.sqrt"))
        check("the old-layout file is shorter than the model", nsqrt > 0, f"{nsqrt} sqrt rows")
        x = P.read(sav)
        check("nudged width read from the old-layout file", np.isclose(x[target.free_index], nudged))
        others = [i for i in range(nfree) if i != target.free_index]
        check("every other value unchanged (missing rows keep the start)", np.allclose(x[others], x0[others]))
        with Session(azr, parameter_file=sav) as f:
            full = np.asarray(f._engine.params_all_rwa(), float)
            rest = [i for i in range(allrwa.size) if i != target.index]
            check("a session made with the file holds it (missing rows keep the project's)",
                  np.isclose(full[target.index], nudged) and np.allclose(full[rest], allrwa[rest]))
            check("its x0 is the file's free vector",
                  f.parameters.x0.size == nfree and np.isclose(f.parameters.x0[target.free_index], nudged))
        bf = best_fit_params(m, sav)
        check("bands.best_fit_params accepts the old layout",
              bf.size == nfree and np.isclose(bf[target.free_index], nudged))
        check("... and leaves the session as it was", np.array_equal(m.parameters.x0, x0)
              and np.array_equal(np.asarray(m._engine.params_all_rwa(), float), allrwa))
        with Session(azr, parameter_file=sav) as f:
            check("... a session made with the file: its own x0",
                  np.array_equal(best_fit_params(f, sav), f.parameters.x0) and np.array_equal(f.parameters.x0, bf))
        # --- the file written by name: the tag first, every parameter, the values of x
        best = P.x0
        best[target.free_index] = allrwa[target.index] * 1.10
        new = os.path.join(work, "output", "param.sav.new")
        P.write(new, best)
        lines = [l.split() for l in open(new) if l.split()]
        check("the written file has the tag line first", lines[0][0] == "#parametrization")
        written = read_parameter_file(new)[1]
        check("the written file lists every parameter", all(n in written for n in names))
        check("the written file holds x's value of the nudged width",
              np.isclose(written[target.name], allrwa[target.index] * 1.10))
        check("numpy.loadtxt sees only parameter rows", np.loadtxt(new, usecols=(1,)).size == len(lines) - 1)

        # --- new layout with the tag, and the 10-06 spelling without the '#'
        for tagname in ("#parametrization", "parametrization"):
            with open(sav, "w") as fh:
                fh.write(f"{tagname:>28s} {1.0: .7e} {0.0: .7e}\n")
                for p in P:
                    fh.write(f"{p.name:>28s} {allrwa[p.index]: .7e} {0.0: .7e}\n")
            basis, values = read_parameter_file(sav)
            check(f"a file tagged '{tagname}' reads every value",
                  basis == 1 and np.allclose([values[n] for n in names], allrwa))

        # --- a foreign file is refused
        with open(sav, "w") as fh:
            fh.write("foo_1 1.0 0.0\nbar_2 2.0 0.0\n")
        try:
            P.read(sav)
            check("a file with no matching names is refused (Parameters.read)", False, "no error raised")
        except ValueError:
            check("a file with no matching names is refused (Parameters.read)", True)
        try:
            Session(azr, parameter_file=sav).close()
            check("a file with no matching names is refused (parameter_file=)", False, "no error raised")
        except ValueError:
            check("a file with no matching names is refused (parameter_file=)", True)

# --- a file that moves a fixed level parameter (T2 has fixed widths): best_fit_params refuses it
with tempfile.TemporaryDirectory() as tmp:
    t2 = os.path.join(tmp, "T2")
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T2"), t2)
    os.makedirs(os.path.join(t2, "output"), exist_ok=True)
    with Session(os.path.join(t2, "T2.azr")) as m:
        allrwa = np.asarray(m._engine.params_all_rwa(), float)
        fixed = next(p for p in m.parameters if p.fixed and p.kind == "width" and allrwa[p.index] != 0.0)
        moved = os.path.join(t2, "output", "moved_fixed.sav")
        with open(moved, "w") as fh:
            for p in m.parameters:
                v = allrwa[p.index] * (1.5 if p.index == fixed.index else 1.0)
                fh.write(f"{p.name:>28s} {v: .7e} {0.0: .7e}\n")
        try:
            best_fit_params(m, moved)
            check("a file that moves a fixed width is refused (best_fit_params)", False, "no error raised")
        except pyazr.ParameterError as e:
            check("a file that moves a fixed width is refused (best_fit_params)",
                  fixed.name in str(e) and "parameter_file=" in str(e), str(e))

print()
if failures:
    print(f"FAILED ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
print("PASS")
