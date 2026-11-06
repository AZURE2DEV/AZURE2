#!/usr/bin/env python3
"""pyazr follows the engine's parameter layout across data and extrapolation mode.

The engine lists one normalization (when varied) and two energy-shift
parameters per segment of the current segment set, so extrapolation mode
(the <segmentsTest> grids) has a different parameter count from data mode.
pyazr cached ``fixed_params`` in data mode only: after ``extrap_mode()`` the
parameter table raised "parameter_info returned N parameters but params_fixed
reported M" and ``_all_rwa`` indexed past the cached flags (seen on a 12C+12C
model: 220 against 212 with eight test segments).

On tests/15N_p_a (three data segments with varied norms) with five test
segments added:

  1. the engine's counts differ between the modes (else the test is vacuous);
  2. in every mode, parameter_info, params_fixed (engine and cached) and
     ``parameters`` have the same length, also when the table is first built
     in extrapolation mode;
  3. the free R-matrix block (names, n_rmatrix) is the same in both modes;
  4. in extrapolation mode a data-mode free vector still works:
     ``_all_rwa`` and ``calculate_rwa`` give what the engine's own vector gives;
  5. a free vector set by the caller survives extrap -> data, and chi2 back in
     data mode equals the chi2 before the switch (1e-12);
  6. a session started in extrapolation mode switches to data mode with the
     data layout and the chi2 of a fresh data-mode session.

Run from anywhere:  python3 tests/pyazr/extrap_layout_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
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
    from pyazr import azure2, AzrModel
    from pyazr.parameters import Parameter
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def counts(m, table=True):
    s = m.sess
    rows = np.asarray(s.parameter_info(), float).reshape(-1, Parameter._NFIELDS).shape[0]
    return rows, len(s.params_fixed()), len(m.fixed_params), (len(m.parameters) if table else rows)


def rmatrix_names(m):
    return [p.name for p in m.parameters if p.kind in ("energy", "width") and not p.fixed]


with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "15N")
    shutil.copytree(os.path.join(SOURCE, "data"), os.path.join(work, "data"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(work, d))
    azr = [f for f in os.listdir(SOURCE) if f.endswith(".azr")][0]
    mdl = AzrModel.from_file(os.path.join(SOURCE, azr))
    mdl.clear_extrapolations()
    for k in range(5):
        mdl.add_extrapolation(1, 2, 0.1 + 0.1 * k, 0.2 + 0.1 * k, 0.05)
    path = os.path.join(work, "layout.azr")
    mdl.write(path)

    with azure2(path, cwd=work) as m:
        x0 = np.asarray(m.params_rwa, float).copy()
        chi0 = m.calculate_chi2_rwa(x0)[0]
        n_data = counts(m, table=False)       # leave pyazr's parameter table unbuilt
        # 2 (extrapolation first): the table is first built in extrapolation mode
        m.extrap_mode()
        try:
            n_ext = counts(m)
            built = True
        except RuntimeError as err:
            n_ext, built = None, False
            check("parameter table builds in extrapolation mode", False, str(err))
        if built:
            check("engine layouts differ between the modes", n_ext[0] != n_data[0], f"{n_data} {n_ext}")
            check("extrapolation mode: all counts agree", len(set(n_ext)) == 1, str(n_ext))
            names_ext, nrm_ext = rmatrix_names(m), m.n_rmatrix
            y_own = m.calculate_rwa(m.params_rwa)
            y_x0 = m.calculate_rwa(x0)
            same = all(np.allclose(a, b, rtol=1e-12, atol=0) for a, b in zip(y_own, y_x0))
            check("extrapolation mode: data-mode vector gives the engine's curves", same)
            try:
                a = m._all_rwa(x0)
                ok = a.size == n_ext[0] and np.allclose(a[:nrm_ext], m._all_rwa(m.params_rwa)[:nrm_ext], rtol=0, atol=0)
                check("extrapolation mode: _all_rwa of a data-mode vector", ok, f"{a.size} vs {n_ext[0]}")
            except Exception as err:
                check("extrapolation mode: _all_rwa of a data-mode vector", False, repr(err))
        # 5: a caller's vector survives the round trip
        xs = np.asarray(m.params_rwa, float).copy()
        xs[0] *= 1.01
        m.params_rwa = xs
        m.data_mode()
        n_back = counts(m)
        check("data mode again: all counts agree with the first data-mode table", len(set(n_back)) == 1 and n_back == n_data,
              f"{n_data} {n_back}")
        if built:
            check("free R-matrix block identical in both modes", rmatrix_names(m) == names_ext and m.n_rmatrix == nrm_ext)
        xb = np.asarray(m.params_rwa, float)
        check("caller's R-matrix value kept through extrap -> data", xb.size == x0.size and abs(xb[0] / (x0[0] * 1.01) - 1) < 1e-14,
              f"{xb[0]} vs {x0[0] * 1.01}")
        check("caller's norm/shift tail kept as the engine's", np.array_equal(xb[m.n_rmatrix:], x0[m.n_rmatrix:]))
        chi1 = m.calculate_chi2_rwa(x0)[0]
        check("chi2 unchanged after the round trip", abs(chi1 / chi0 - 1) < 1e-12, f"{chi0!r} {chi1!r}")

    # 6: a session started in extrapolation mode
    with azure2(path, cwd=work, data_mode=False) as m:
        n0 = counts(m)
        check("extrapolation start: counts agree", len(set(n0)) == 1, str(n0))
        m.data_mode()
        n1 = counts(m)
        check("extrapolation start -> data mode: data layout", n1 == n_data, f"{n1} vs {n_data}")
        chi2 = m.calculate_chi2_rwa(np.asarray(m.params_rwa, float))[0]
        check("extrapolation start -> data mode: chi2 of a data session", abs(chi2 / chi0 - 1) < 1e-12, f"{chi0!r} {chi2!r}")

print("FAILED: " + ", ".join(failures) if failures else "all checks passed")
sys.exit(1 if failures else 0)
