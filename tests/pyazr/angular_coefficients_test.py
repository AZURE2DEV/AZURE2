#!/usr/bin/env python3
"""Angular-coefficient data segments (isDiff 9) through pyazr (2026-10-26).

  1. AzrModel.add_data_segment(observable="angular-coefficients") writes an isDiff 9 line, with "absolute 1"
     when asked, and refuses absolute= on another observable; the datasets reader gives the observable and
     the flag back.
  2. With the compiled engine (skipped without it): a model whose coefficient data are the model's own
     Legendre coefficients (from an isDiff 3 extrapolation of the same model) evaluates to chi2 ~ 0, and
     calculate_rwa returns the coefficients in the data slots.

Run from anywhere:  python3 tests/pyazr/angular_coefficients_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
from pyazr import AzrModel  # noqa: E402
from pyazr.datasets import SegmentSet  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def main():
    import numpy as np
    src = os.path.join(ROOT, "tests", "angular_coefficients", "models", "extrap.azr")
    work = tempfile.mkdtemp(prefix="pyazr_acoef_")
    try:
        print("1. the segment line")
        model = AzrModel.from_file(src)
        os.makedirs(os.path.join(work, "data"), exist_ok=True)
        open(os.path.join(work, "data", "ak.dat"), "w").write("6.0  2  -0.5  0.01\n")
        model.add_data_segment("data/ak.dat", 2, 5, observable="angular-coefficients", energy_min=5.9, energy_max=6.4)
        model.add_data_segment("data/ak.dat", 2, 5, observable="angular-coefficients", energy_min=5.9, energy_max=6.4,
                               absolute=True)
        try:
            model.add_data_segment("data/ak.dat", 2, 5, observable="differential", absolute=True)
            check("absolute= on another observable is refused", False)
        except ValueError:
            check("absolute= on another observable is refused", True)
        azr = os.path.join(work, "m.azr")
        model.write(azr)
        lines = [l for l in open(azr).read().split("<segmentsData>")[1].split("</segmentsData>")[0].splitlines() if l.strip()]
        check("two isDiff 9 lines, the second with 'absolute 1'",
              len(lines) == 2 and lines[0].split()[7] == "9" and "absolute" not in lines[0] and lines[1].split()[-2:] == ["absolute", "1"], lines)
        segs = SegmentSet.from_file(azr)
        check("datasets: observable angular-coefficients, absolute flag read back",
              [s.observable for s in segs] == ["angular-coefficients"] * 2 and [s.absolute for s in segs] == [False, True],
              [(s.observable, s.absolute) for s in segs])

        print("\n2. through the engine")
        try:
            from pyazr import azure2
            import pyazr._azure2  # noqa: F401
        except Exception as e:
            print(f"  skip  engine not importable ({e})")
            return
        # the model's own coefficients: an isDiff 3 extrapolation at the data energies
        ext = AzrModel.from_file(src)
        ext.set_output_dir("output/")
        os.makedirs(os.path.join(work, "ext", "output"), exist_ok=True)
        os.makedirs(os.path.join(work, "ext", "checks"), exist_ok=True)
        ext.write(os.path.join(work, "ext", "m.azr"))
        with azure2(os.path.join(work, "ext", "m.azr"), cwd=os.path.join(work, "ext"), data_mode=False) as m:
            coef = m.calculate_angular_dists_rwa(m.params_rwa)[0]   # test segment 1: one array per energy
        n_e = len(coef)
        check("isDiff 3 extrapolation gives coefficient arrays", n_e == 11 and all(len(c) >= 5 for c in coef), (n_e, [len(c) for c in coef][:3]))
        # data at lab 6.0-6.3 (grid 5.9 + 0.05 i, i = 2..5): rows E_lab, k, a_k, err
        with open(os.path.join(work, "data", "ak.dat"), "w") as f:
            for i in range(2, 6):
                for k in (2, 4):
                    f.write(f"{5.9 + 0.05 * i:.6f}  {k}  {coef[i][k]:.12e}  0.01\n")
        dm = AzrModel.from_file(src)
        dm.set_output_dir("output/")
        dm.add_data_segment("data/ak.dat", 2, 5, observable="angular-coefficients", energy_min=5.9, energy_max=6.4)
        os.makedirs(os.path.join(work, "output"), exist_ok=True)
        os.makedirs(os.path.join(work, "checks"), exist_ok=True)
        dm.write(azr)
        with azure2(azr, cwd=work) as m:
            chi2 = float(m.calculate_chi2_rwa(m.params_rwa)[0])
            vals = np.asarray(m.calculate_rwa(m.params_rwa)[0], float)
            data = np.asarray(m.cross[0], float)
            ks = np.asarray(m.angles[0], float)
        check("the coefficient data segment evaluates to chi2 ~ 0", chi2 < 1e-8, chi2)
        check("calculate_rwa returns the coefficients (k in the angle slot)",
              np.allclose(vals, data, rtol=1e-8, atol=0) and sorted(set(ks.astype(int))) == [2, 4], (vals[:3], data[:3], ks[:3]))
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
    print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
    sys.exit(1 if failures else 0)
