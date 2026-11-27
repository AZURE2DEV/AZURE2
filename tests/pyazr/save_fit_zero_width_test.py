#!/usr/bin/env python3
"""save_fit keeps a width stated as zero at zero, also when external capture gives its channel a physical width.

A width stated as 0 in <levels> is fixed whatever its flag.  With external capture a photon channel whose width is
stated as 0 still has a physical width: the external part has an imaginary portion no internal amplitude cancels, the
engine holds the internal amplitude that brings the total closest to zero and warns ("Imaginary portion of the
external width is greater than the total width").  save_fit wrote that physical width into the row.  It read back as
the same amplitude, so the chi-squared agreed and save_fit's own check passed, but the row was no longer zero and the
channel was free: the snapshot had more free parameters than the fit it was taken from (found by rmfit's bake check,
2026-11-27).

Model: tests/13N.  Its 3/2+ and 5/2- background levels at 20 MeV state 0 for their E1 and E2 photon channels.

  1. external capture on: the session warns about the two channels, and their physical widths are not zero (the case
     is there);
  2. the snapshot states 0 for every width the project states 0 for, and reopens with the same parameters, the same
     free ones, the same amplitudes and the same objective;
  3. external capture ignored: the same, and no width stated as 0 has a physical width.

Run from anywhere:  python3 tests/pyazr/save_fit_zero_width_test.py
"""
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "13N")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import Project, Session
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def zero_rows(path):
    """(level label, pair, l, 2s) of every active row that states a width of 0."""
    return {(level.label, c.pair, c.l, int(round(2 * c.s)))
            for level in Project.open(path).levels() if level.active for c in level.channels if c.width == 0.0}


def state(session):
    p = session.parameters
    x = np.asarray(p.x0, float)
    return ([q.name for q in p], [q.name for q in p if not q.fixed], np.asarray(session._all_rwa(x), float),
            np.asarray(session.transform(x), float), float(session.evaluate(x).objective))


with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "13N")
    shutil.copytree(SOURCE, work)
    for d in ("output", "checks"):
        os.makedirs(os.path.join(work, d), exist_ok=True)
    source = os.path.join(work, "13N.azr")
    stated = zero_rows(source)
    for part, ignore in (("1-2", False), ("3", True)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s = Session(source, ignore_externals=ignore, use_brune=True, use_long_wavelength=True)
            names, free, rwa, physical, objective = state(s)
            zero_stated = [q for q in s.parameters if q.kind == "width" and q.fixed and rwa[q.index] != 0.0
                           and physical[q.index] != 0.0 and
                           any(k[1:] == (next(pr.key for pr in s.pairs if pr.number == q.pair), q.L,
                                         int(round(2 * (q.S or 0.0)))) for k in stated)]
            with_width = [q.name for q in zero_stated
                          if Project.open(source).level((q.jgroup, q.level)).channel(
                              next(pr.key for pr in s.pairs if pr.number == q.pair), q.L, q.S or 0.0).width == 0.0]
            if not ignore:
                check("1. external capture on: two widths stated as 0 have a physical width",
                      sorted(with_width) == ["width_7_2", "width_8_2"], f"{with_width}")
            else:
                check("3. external capture ignored: no width stated as 0 has a physical width", not with_width,
                      f"{with_width}")
            out = os.path.join(work, f"saved_{int(ignore)}.azr")
            s.save_fit(out, s.parameters.x0, param_sav=False)
            back = Session(out, ignore_externals=ignore, use_brune=True, use_long_wavelength=True)
            names2, free2, rwa2, physical2, objective2 = state(back)
        label = "2" if not ignore else "3"
        check(f"{label}. the snapshot states 0 for the widths the project states 0 for", zero_rows(out) == stated,
              f"{sorted(stated - zero_rows(out))} are no longer 0")
        check(f"{label}. the snapshot has the fit's parameters and the same free ones",
              names2 == names and free2 == free, f"{len(names)} -> {len(names2)} parameters, {len(free)} -> {len(free2)} free")
        check(f"{label}. the amplitudes read back (1e-9)",
              rwa2.shape == rwa.shape and bool(np.all(np.abs(rwa2 - rwa) <= 1e-9 * np.maximum(1.0, np.abs(rwa)))),
              "the amplitudes differ")
        check(f"{label}. the objective reads back (1e-9)", abs(objective2 - objective) <= 1e-9 * abs(objective),
              f"{objective2} vs {objective}")

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: save_fit keeps a width stated as zero at zero")
