#!/usr/bin/env python3
"""Angular-coefficient data sets (``quantity legendre``) through pyazr (2026-10-26).

  1. Project.add writes a ``legendre`` data set with ``expansion a_over_a0`` (unit 1) or ``absolute`` (b/sr), and
     refuses an expansion on another quantity; the records give the expansion back.
  2. With the compiled engine (skipped without it): a model whose coefficient data are the model's own
     Legendre coefficients (from the ``legendre`` extrapolation of the same model) evaluates to chi2 ~ 0, and
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
try:
    from pyazr import Project, FormatError  # noqa: E402
except Exception as err:                                   # engine not built: the package does not import
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


REACTION = "13C(a,n)16O[6.12989]"


def coefficient_set(model, name, absolute=False):
    """A legendre data set on data/ak.dat at lab 5.9-6.4 MeV."""
    model.add("dataset", name, reaction=REACTION, quantity="legendre", unit="b/sr" if absolute else "1",
              expansion="absolute" if absolute else "a_over_a0", file="data/ak.dat",
              columns="E_lab[MeV] order value error", energy_window="5.9 6.4 MeV", norm="1 fixed")


def main():
    import numpy as np
    src = os.path.join(ROOT, "tests", "angular_coefficients", "models", "extrap.azr")
    work = tempfile.mkdtemp(prefix="pyazr_acoef_")
    try:
        print("1. the data set records")
        os.makedirs(os.path.join(work, "data"), exist_ok=True)
        open(os.path.join(work, "data", "ak.dat"), "w").write("6.0  2  -0.5  0.01\n")
        azr = os.path.join(work, "m.azr")
        shutil.copy(src, azr)
        model = Project.open(azr)
        coefficient_set(model, "ak")
        coefficient_set(model, "ak_abs", absolute=True)
        try:
            model.add("dataset", "bad", reaction=REACTION, quantity="dsigma_domega", unit="b/sr", expansion="absolute",
                      file="data/ak.dat", columns="E_lab[MeV] theta_cm[deg] value error", norm="1 fixed")
            check("an expansion on another quantity is refused", False)
        except FormatError:
            check("an expansion on another quantity is refused", "bad" not in model)
        model.save()
        back = Project.open(azr)
        sets = [r for r in back.records().values() if r.record_kind == "dataset"]
        check("two legendre data sets, the second absolute (b/sr)",
              [(r["quantity"], r["expansion"], r["unit"]) for r in sets]
              == [("legendre", "a_over_a0", "1"), ("legendre", "absolute", "b/sr")],
              [r.fields for r in sets])

        print("\n2. through the engine")
        try:
            from pyazr import azure2
            import pyazr._azure2  # noqa: F401
        except Exception as e:
            print(f"  skip  engine not importable ({e})")
            return
        # the model's own coefficients: its legendre extrapolation at the data energies
        os.makedirs(os.path.join(work, "ext", "output"), exist_ok=True)
        os.makedirs(os.path.join(work, "ext", "checks"), exist_ok=True)
        shutil.copy(src, os.path.join(work, "ext", "m.azr"))
        with azure2(os.path.join(work, "ext", "m.azr"), cwd=os.path.join(work, "ext"), data_mode=False) as m:
            coef = m.calculate_angular_dists_rwa(m.params_rwa)[0]   # test segment 1: one array per energy
        n_e = len(coef)
        check("the legendre extrapolation gives coefficient arrays", n_e == 11 and all(len(c) >= 5 for c in coef), (n_e, [len(c) for c in coef][:3]))
        # data at lab 6.0-6.3 (grid 5.9 + 0.05 i, i = 2..5): rows E_lab, k, a_k, err
        with open(os.path.join(work, "data", "ak.dat"), "w") as f:
            for i in range(2, 6):
                for k in (2, 4):
                    f.write(f"{5.9 + 0.05 * i:.6f}  {k}  {coef[i][k]:.12e}  0.01\n")
        shutil.copy(src, azr)
        dm = Project.open(azr)
        coefficient_set(dm, "ak")
        os.makedirs(os.path.join(work, "output"), exist_ok=True)
        os.makedirs(os.path.join(work, "checks"), exist_ok=True)
        dm.save()
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
