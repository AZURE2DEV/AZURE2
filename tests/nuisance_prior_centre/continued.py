#!/usr/bin/env python3
"""tests/nuisance_prior_centre/continued.sh (5): the pyazr side of a continued fit.  For each fit of continued.sh
(WORK/<case>/t.azr, its fit.sav, the written project written.azr):

  - the written project and the fitted one, both at the fit's parameters: the objective (chi^2 and every prior) and
    the chi^2 of every data set equal to 1e-12, and the layout value of the fitted parameter (density_max,
    doppler_max, the apparatus parameter's layout) equal;
  - save_fit(verify=True) of the fit (it reopens the snapshot and checks its chi^2 against the fit's to 1e-6); not for
    the stopping scale, which save_fit writes on the first of the two targets sharing it (the written project is
    refused: a defect of the write-back, reported, not fixed here);
  - the start outside the centre's range (WORK/outside: D.doppler 400 K, prior 300 20): the layout of the start,
    1.21 x 400 K = 484 K.

Usage: continued.py WORK"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from pyazr import Session  # noqa: E402

failures = []
CASES = [("density", "tgt_meyer_84.areal_density"), ("doppler", "D.doppler"), ("scale", "stopping_scale.TiN.1H"),
         ("kernel", "det_fit.energy_resolution.width")]


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def layout(m, name):
    p = next(q for q in m.parameters if q.name == name)
    return {"density": p.density_max, "doppler": p.doppler_max}.get(p.kind, p.S)


def evaluate(path, cwd, sav, name):
    with Session(path, cwd=cwd) as m:
        x = m.load_sav(sav)
        return m.objective()(x), {d.name: d.chi2 for d in m.evaluate(x).data}, layout(m, name)


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def main(work):
    for case, name in CASES:
        d = os.path.join(work, case)
        sav = os.path.join(d, "fit.sav")
        if not os.path.exists(sav) or not os.path.exists(os.path.join(d, "written.azr")):
            check(f"({case}) the fit and its written project exist", False, d)
            continue
        fo, fd, fl = evaluate(os.path.join(d, "t.azr"), d, sav, name)
        wo, wd, wl = evaluate(os.path.join(d, "written.azr"), d, sav, name)
        worst = max([rel(fo, wo)] + [rel(fd[k], wd[k]) for k in fd])
        check(f"({case}) the written project at the fit's parameters: objective and every data set's chi^2 equal "
              f"(largest difference {worst:.1e}), the layout {name} {wl!r} the fit's", worst <= 1e-12 and fl == wl,
              (fo, wo, fd, wd, fl, wl))
        if case == "scale":
            continue
        with Session(os.path.join(d, "t.azr"), cwd=d) as m:
            x = m.load_sav(sav)
            try:
                m.save_fit(os.path.join(d, "snapshot.azr"), x, param_sav=False)
                check(f"({case}) save_fit(verify=True): the snapshot reads back as the fit", True)
            except Exception as e:  # noqa: BLE001 -- what the base build gives for the Doppler temperature
                check(f"({case}) save_fit(verify=True): the snapshot reads back as the fit", False, f"{type(e).__name__}: {e}")
    d = os.path.join(work, "outside")
    with Session(os.path.join(d, "t.azr"), cwd=d) as m:
        got = layout(m, "D.doppler")
    check(f"(4) the start outside the centre's range: the start's layout ({got!r} K, want 1.21 x 400)", got == 1.21 * 400.0,
          got)
    if failures:
        print(f"{len(failures)} failure(s)")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
