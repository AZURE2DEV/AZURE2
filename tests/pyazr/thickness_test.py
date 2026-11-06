#!/usr/bin/env python3
"""The fitted sample thickness of transmission / self-shielded-yield segments through pyazr.

Pure Python: AzrModel.sample_thicknesses / set_sample_thickness / add_data_segment(transmission=..., yield_=...,
thickness_group=...) round trips and refusals; the datasets Segment fields (thickness_vary, thickness_sigma,
thickness_group) of the fit form.

Through the engine, on tests/transmission_yield's grad_fit model (eight transmission / yield lines with the fit form
"n vary sigma": 1, 2 without a target effect, 3, 4 under a Gaussian spread with a 0.01 atoms/b prior, 5, 6 a
thicknessgroup under doppler alone, 7, 8 under doppler + the Gaussian; every width free): the type-10 parameter
records (segment 6 the fixed follower of 5), penalties()["thickness"] and objective() at moved thicknesses, and the
segment_<k>_thickness entries of chi2_and_grad and the residual_jacobian columns against central differences.

The engine part needs the compiled module and skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/thickness_test.py
"""
import importlib.util
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SUITE = os.path.join(ROOT, "tests", "transmission_yield")
sys.path.insert(0, ROOT)

try:
    from pyazr.azrfile import AzrModel  # noqa: E402
    from pyazr.datasets import SegmentSet  # noqa: E402
except ImportError as e:
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(0)

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({detail})"))
    if not ok:
        failures.append(name)


def data_lines(path):
    text = open(path).read()
    return text.split("<segmentsData>\n")[1].split("</segmentsData>")[0].strip().splitlines()


tmp = tempfile.mkdtemp(prefix="pyazr_thickness_")
try:
    spec = importlib.util.spec_from_file_location("ty_models", os.path.join(SUITE, "make_models.py"))
    mm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mm)
    mm.n12c_models(os.path.join(ROOT, "tests", "reaction_rate", "models", "n12C.azr"), tmp)

    # --- AzrModel and datasets (no engine) ------------------------------------------------------
    m = AzrModel.from_file(os.path.join(tmp, "thin.azr"))
    before = data_lines(os.path.join(tmp, "thin.azr"))
    t = m.sample_thicknesses
    check("sample_thicknesses reads the plain blocks",
          t == {1: dict(kind="transmission", n=0.15, fit=False, vary=False, sigma=0.0, group=None, active=True),
                2: dict(kind="yield", n=0.15, fit=False, vary=False, sigma=0.0, group=None, active=True)}, t)
    path = os.path.join(tmp, "rt.azr")
    m.set_sample_thickness(1, vary=True, sigma=0.01, group="s")
    m.set_sample_thickness(2, n=0.12, vary=False, sigma=0, group="s")
    m.write(path)
    after = data_lines(path)
    check("set_sample_thickness writes the fit form and the group, the rest of the line kept",
          after[0] == before[0].replace("transmission 0.15", "transmission 0.15 1 0.01") + "  thicknessgroup s" and
          after[1] == before[1].replace("yield 0.15", "yield 0.12 0 0") + "  thicknessgroup s", after)
    r = AzrModel.from_file(path).sample_thicknesses
    check("the file reads back (fit, vary, sigma, group)",
          r[1] == dict(kind="transmission", n=0.15, fit=True, vary=True, sigma=0.01, group="s", active=True) and
          r[2] == dict(kind="yield", n=0.12, fit=True, vary=False, sigma=0.0, group="s", active=True), r)
    segs = SegmentSet.from_file(path)
    check("datasets Segment: thickness, thickness_vary, thickness_sigma, thickness_group",
          (segs[0].thickness, segs[0].thickness_vary, segs[0].thickness_sigma, segs[0].thickness_group,
           segs[0].thickness_fit, segs[0].total) == (0.15, True, 0.01, "s", True, True) and
          (segs[1].thickness, segs[1].thickness_vary, segs[1].thickness_group) == (0.12, False, "s"),
          (segs[0], segs[1]))
    plain = SegmentSet.from_file(os.path.join(tmp, "thin.azr"))[0]
    check("a plain block: thickness_vary None, not thickness_fit",
          plain.thickness_vary is None and not plain.thickness_fit and plain.thickness == 0.15, plain)
    m2 = AzrModel.from_file(path)
    m2.set_sample_thickness(1, n=0.2, group="")
    m2.write(path)
    check("set_sample_thickness(n, group='') keeps the fit form and drops the group",
          data_lines(path)[0] == before[0].replace("transmission 0.15", "transmission 0.2 1 0.01"), data_lines(path)[0])
    m3 = AzrModel.from_file(os.path.join(tmp, "bare.azr"))
    m3.set_sample_thickness(2, kind="yield", n=0.1)
    m3.write(path)
    check("set_sample_thickness(kind, n) adds a plain block to a line without one",
          data_lines(path)[1].split()[-2:] == ["yield", "0.1"], data_lines(path)[1])
    for what, call in (("a negative sigma", lambda: AzrModel.from_file(os.path.join(tmp, "thin.azr"))
                        .set_sample_thickness(1, vary=True, sigma=-1)),
                       ("n <= 0", lambda: AzrModel.from_file(os.path.join(tmp, "thin.azr")).set_sample_thickness(1, n=0)),
                       ("a line without a block and no kind", lambda: AzrModel.from_file(os.path.join(tmp, "bare.azr"))
                        .set_sample_thickness(1, n=0.1)),
                       ("a group label with a blank", lambda: AzrModel.from_file(os.path.join(tmp, "thin.azr"))
                        .set_sample_thickness(1, group="a b"))):
        try:
            call()
            check(f"set_sample_thickness refuses {what}", False)
        except ValueError:
            check(f"set_sample_thickness refuses {what}", True)
    a = AzrModel.from_file(os.path.join(tmp, "bare.azr"))
    a.add_data_segment("data/d.dat", 1, 1, energy_min=0.29, energy_max=0.66, transmission=(0.15, True, 0.02),
                       thickness_group="g")
    a.add_data_segment("data/g.dat", 1, 2, energy_min=0.29, energy_max=0.66, yield_=0.15, thickness_group="g")
    a.write(path)
    added = data_lines(path)[-2:]
    check("add_data_segment(transmission=(n, vary, sigma), yield_=n, thickness_group=...)",
          added[0].split()[-6:] == ["transmission", "0.15", "1", "0.02", "thicknessgroup", "g"] and
          added[1].split()[-4:] == ["yield", "0.15", "thicknessgroup", "g"], added)
    for what, kw in (("both", dict(transmission=0.1, yield_=0.1)), ("a group without a block", dict(thickness_group="x")),
                     ("a bad tuple", dict(transmission=(0.1, 1)))):
        try:
            AzrModel.from_file(os.path.join(tmp, "bare.azr")).add_data_segment("data/d.dat", 1, 1, **kw)
            check(f"add_data_segment refuses {what}", False)
        except ValueError:
            check(f"add_data_segment refuses {what}", True)

    # --- engine ----------------------------------------------------------------------------------
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import numpy as np
        from pyazr import azure2
        from pyazr import _azure2  # noqa: F401
    except Exception as err:
        print(f"skip engine part: {type(err).__name__}: {err}")
    else:
        work = os.path.join(tmp, "run_grad_fit")
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        shutil.copy(os.path.join(tmp, "grad_fit.azr"), work)
        os.symlink(os.path.join(tmp, "data"), os.path.join(work, "data"))
        with azure2(os.path.join(work, "grad_fit.azr"), cwd=work) as s:
            tp = s.parameters.thicknesses
            ok = (len(tp) == 8 and [p.segment_key for p in tp] == list(range(1, 9)) and
                  [p.name for p in tp] == ["segment_%d_thickness" % k for k in range(1, 9)] and
                  [p.fixed for p in tp] == [False] * 5 + [True] + [False] * 2 and
                  [p.prior_sigma for p in tp] == [0, 0, 0.01, 0.01, 0, 0, 0, 0] and
                  all(p.prior_nominal == 0.15 for p in tp) and
                  [p.sample_kind for p in tp] == ["transmission", "yield"] * 4 and
                  [p.thickness_leader for p in tp] == [None] * 5 + [5] + [None] * 2)
            check("eight type-10 records: priors on 3, 4; segment 6 the fixed follower of 5", ok, list(tp))
            free = [p for p in tp if not p.fixed]
            x = np.asarray(s.params_rwa, float)
            pen0 = s.penalties(x)["thickness"]
            xm = x.copy()
            for p in free:
                xm[p.free_index] = {3: 0.16, 4: 0.17}.get(p.segment_key, 0.15 + 0.004 * p.segment_key)
            pen = s.penalties(xm)["thickness"]
            obj = s.objective(xm) - float(np.sum(s.calculate_chi2_rwa(xm)))
            check("penalties()['thickness'] = 1 + 4 at segment_3 = 0.16, segment_4 = 0.17; 0 at the start",
                  abs(pen - 5.0) < 1e-9 and pen0 == 0.0, (pen, pen0))
            check("objective() carries it", abs(obj - 5.0) < 1e-6, obj)
            c0 = float(np.sum(s.calculate_chi2_rwa(x)))
            c1 = float(np.sum(s.calculate_chi2_rwa(xm)))
            check("moving the thicknesses moves the data chi2", abs(c1 - c0) > 1e-6 * abs(c0), (c0, c1))
            for name, xx in (("start", x), ("moved", xm)):
                c, g = s.chi2_and_grad(xx)
                g = np.asarray(g, float)
                worst = 0.0
                for p in free:
                    h = 1e-6
                    xp, xn = xx.copy(), xx.copy()
                    xp[p.free_index] += h
                    xn[p.free_index] -= h
                    fd = (float(np.sum(s.calculate_chi2_rwa(xp))) - float(np.sum(s.calculate_chi2_rwa(xn)))) / (2 * h)
                    worst = max(worst, abs(g[p.free_index] - fd) / max(abs(fd), 1e-300))
                check(f"chi2_and_grad d chi2/d segment_k_thickness = central differences ({name}, 7 of them, rel 1e-6)",
                      worst < 1e-6, worst)
            r0, J = s.residual_jacobian(xm)
            J = np.asarray(J, float)
            worst = 0.0
            for p in free:
                h = 1e-6
                xp, xn = xm.copy(), xm.copy()
                xp[p.free_index] += h
                xn[p.free_index] -= h
                rp, _ = s.residual_jacobian(xp)
                rn, _ = s.residual_jacobian(xn)
                fd = (np.asarray(rp, float) - np.asarray(rn, float)) / (2 * h)
                worst = max(worst, float(np.max(np.abs(J[:, p.free_index] - fd)) / np.max(np.abs(fd))))
            check("residual_jacobian columns segment_k_thickness = central differences (1e-6 of the largest)",
                  worst < 1e-6, worst)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: pyazr sample thickness")
