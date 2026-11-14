#!/usr/bin/env python3
"""pyazr follows the engine's parameter layout of a data session and of an extrapolation session.

The engine lists one normalization (when varied) and two energy-shift
parameters per segment of the current segment set, so extrapolation mode
(the <segmentsTest> grids) has a different parameter count from data mode.
pyazr once cached the fixed flags in data mode only: after the mode switch of that time the
parameter table raised "N parameter records but M fixed flags" and ``_all_rwa`` indexed past the cached flags (seen on a 12C+12C
model: 220 against 212 with eight test segments).

On tests/15N_p_a (three data segments with varied norms) with five test
segments added:

  1. the engine's counts differ between the modes (else the test is vacuous);
  2. in each session, parameter_records, params_fixed (engine and cached) and ``parameters`` have the same
     length;
  3. the free R-matrix block (names, count) is the same in both sessions;
  4. the data session's free vector carried to the extrapolation session
     (``pyazr.parameters.carry_model_parameters``): ``_all_rwa`` and ``evaluate().data[i].model`` give what
     the extrapolation session's own vector gives;
  5. a data session made after it has the first one's layout, x0 and chi2 (1e-12).

Run from anywhere:  python3 tests/pyazr/extrap_layout_test.py

This test reads the binding itself (`_azure2.Session`), so it uses the private handle `Session._engine`.
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "15N_p_a")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import Session, Project
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def counts(m, table=True):
    s = m._engine
    rows = len(s.parameter_records())
    return rows, len(s.params_fixed()), len(m._fixed), (len(m.parameters) if table else rows)


def rmatrix_names(m):
    return [p.name for p in m.parameters if p.kind in ("energy", "width") and not p.fixed]


with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "15N")
    shutil.copytree(os.path.join(SOURCE, "data"), os.path.join(work, "data"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(work, d))
    azr = [f for f in os.listdir(SOURCE) if f.endswith(".azr")][0]
    mdl = Project.open(os.path.join(SOURCE, azr))
    for name, r in mdl.records().items():
        if r.record_kind == "extrapolation":
            mdl.remove(name)
    for k in range(5):  # the angle-integrated cross section of pair 1 -> 2, 0.1 + 0.1 k to 0.2 + 0.1 k MeV lab
        mdl.add_extrapolation(f"x{k + 1}", reaction=mdl.reaction(1, 2), quantity="sigma", unit="b",
                              energies=f"range {0.1 + 0.1 * k:.1f} {0.2 + 0.1 * k:.1f} MeV lab step 0.05")
    path = mdl.save(os.path.join(work, "layout.azr"))

    from pyazr.parameters import carry_model_parameters
    with Session(path, cwd=work) as m:
        x0 = np.asarray(m.parameters.x0, float).copy()
        chi0 = m.evaluate(x0).chi2
        n_data = counts(m, table=False)       # leave pyazr's parameter table unbuilt
        check("data session: all counts agree", len(set(counts(m))) == 1, str(counts(m)))
        names_data, nrm_data = rmatrix_names(m), len(m.parameters.free.of_kind("energy", "width"))
        P_data = m.parameters

    # 2, 3, 4: a session made for the extrapolations (its table first built there)
    with Session(path, cwd=work, extrapolations=True) as m:
        try:
            n_ext = counts(m)
            built = True
        except RuntimeError as err:
            n_ext, built = None, False
            check("parameter table builds in an extrapolation session", False, str(err))
        if built:
            check("engine layouts differ between the two sessions", n_ext[0] != n_data[0], f"{n_data} {n_ext}")
            check("extrapolation session: all counts agree", len(set(n_ext)) == 1, str(n_ext))
            check("free R-matrix block identical in both sessions",
                  rmatrix_names(m) == names_data and len(m.parameters.free.of_kind("energy", "width")) == nrm_data)
            x_carried = carry_model_parameters(x0, P_data, m.parameters)
            y_own = [d.model for d in m.evaluate(m.parameters.x0).data]
            y_x0 = [d.model for d in m.evaluate(x_carried).data]
            same = all(np.allclose(a, b, rtol=1e-12, atol=0) for a, b in zip(y_own, y_x0))
            check("extrapolation session: the data session's vector carried over gives the engine's curves", same)
            a = m._all_rwa(x_carried)
            ok = a.size == n_ext[0] and np.array_equal(a[:nrm_data], m._all_rwa(m.parameters.x0)[:nrm_data])
            check("extrapolation session: _all_rwa of the carried vector", ok, f"{a.size} vs {n_ext[0]}")

    # 5: a data session made after it has the first data session's layout, x0 and chi2
    with Session(path, cwd=work) as m:
        n_back = counts(m)
        check("data session again: the first data session's layout", len(set(n_back)) == 1 and n_back == n_data,
              f"{n_data} {n_back}")
        check("data session again: x0 is the first x0", np.array_equal(m.parameters.x0, x0))
        chi1 = m.evaluate(x0).chi2
        check("data session again: chi2 unchanged", abs(chi1 / chi0 - 1) < 1e-12, f"{chi0!r} {chi1!r}")

print("FAILED: " + ", ".join(failures) if failures else "all checks passed")
sys.exit(1 if failures else 0)
