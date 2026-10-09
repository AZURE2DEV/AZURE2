#!/usr/bin/env python3
"""pyazr reads parameter files by NAME, never by row position.

AZURE2 matches `param.sav` rows to parameters by name (AZUREParams::
ReadUserParameters), appends new parameters as the code grows (the
`segment_<key>_energy_shift_sqrt` rows, 2026-10-05) and now writes a
`#parametrization` tag line first (2026-10-06/07).  A positional reader is
off by one or more rows on any file from another layout, with no error --
an 8Be+alpha rate came out exactly 2x that way.  Pinned here on
tests/identical_pp_res:

  * an old-layout file (no _sqrt rows, no tag, one width nudged): the nudged
    width is read, every missing row keeps the session's value, and
    bands.best_fit_params returns a vector of the free-parameter length;
  * update_sav_from_rwa_params writes rows back by name, appends the free
    parameters the file lacks and keeps the tag line;
  * a new-layout file (tag first) reads every value; a `parametrization` row
    without the '#' (files of 2026-10-06) is skipped the same way;
  * a file with no name in common with the model is refused.

Needs numpy and the compiled engine; exits 77 (skipped) without them.
Run from anywhere:  python3 tests/pyazr/sav_readers_test.py
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
    from pyazr import azure2
    from pyazr.bands import best_fit_params
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

    with azure2(azr) as m:
        names = [p.name for p in m.parameters]
        allrwa = np.asarray(m.sess.params_all_rwa(), float)
        fixed = np.asarray(m.fixed_params, float).round().astype(bool)
        widths = [p for p in m.parameters if p.kind == "width" and not p.fixed]
        target = widths[0]
        nudged = allrwa[target.index] * 1.05

        # --- old layout: no _sqrt rows, no tag, one width nudged
        with open(sav, "w") as fh:
            for p in m.parameters:
                if p.name.endswith("_energy_shift_sqrt"):
                    continue
                v = nudged if p.index == target.index else allrwa[p.index]
                fh.write(f"{p.name:>28s} {v: .7e} {0.0: .7e}\n")
        nsqrt = sum(1 for n in names if n.endswith("_energy_shift_sqrt"))
        check("the old-layout file is shorter than the model", nsqrt > 0, f"{nsqrt} sqrt rows")
        full = m.full_rwa_from_sav(sav)
        check("nudged width read from the old-layout file", np.isclose(full[target.index], nudged))
        others = [i for i in range(allrwa.size) if i != target.index]
        check("every other value unchanged (missing rows keep the session's)",
              np.allclose(full[others], allrwa[others]))
        m.update_rwa_params_from_sav()
        check("update_rwa_params_from_sav gives the free vector",
              len(m.params_rwa) == int((~fixed).sum()) and np.isclose(
                  m.params_rwa[target.free_index], nudged))
        bf = best_fit_params(m, sav)
        check("bands.best_fit_params accepts the old layout",
              bf.size == int((~fixed).sum()) and np.isclose(bf[target.free_index], nudged))

        # --- write-back by name appends the missing free rows and keeps a tag
        with open(sav) as fh:
            body = fh.read()
        with open(sav, "w") as fh:
            fh.write(f"{'#parametrization':>28s} {1.0: .7e} {0.0: .7e}\n" + body)
        best = np.asarray(m.params_rwa, float)
        best[target.free_index] = allrwa[target.index] * 1.10
        m.update_sav_from_rwa_params(best)
        new = os.path.join(work, "output", "param.sav.new")
        lines = [l.split() for l in open(new) if l.split()]
        check("write-back keeps the tag line first", lines[0][0] == "#parametrization")
        written = m.read_sav(new)
        free_names = [p.name for p in m.parameters if not p.fixed]
        check("write-back lists every free parameter (missing ones appended)",
              all(n in written for n in free_names))
        check("write-back updated the nudged width by name",
              np.isclose(written[target.name], allrwa[target.index] * 1.10))
        check("numpy.loadtxt sees only parameter rows", np.loadtxt(new, usecols=(1,)).size == len(lines) - 1)

        # --- new layout with the tag, and the 10-06 spelling without the '#'
        for tagname in ("#parametrization", "parametrization"):
            with open(sav, "w") as fh:
                fh.write(f"{tagname:>28s} {1.0: .7e} {0.0: .7e}\n")
                for p in m.parameters:
                    fh.write(f"{p.name:>28s} {allrwa[p.index]: .7e} {0.0: .7e}\n")
            full = m.full_rwa_from_sav(sav)
            check(f"a file tagged '{tagname}' reads every value", np.allclose(full, allrwa))

        # --- a foreign file is refused
        with open(sav, "w") as fh:
            fh.write("foo_1 1.0 0.0\nbar_2 2.0 0.0\n")
        try:
            m.full_rwa_from_sav(sav)
            check("a file with no matching names is refused", False, "no error raised")
        except ValueError:
            check("a file with no matching names is refused", True)

print()
if failures:
    print(f"FAILED ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
print("PASS")
