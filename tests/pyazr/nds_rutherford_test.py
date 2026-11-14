#!/usr/bin/env python3
"""Ratio-to-Rutherford EXFOR data (,DA,,RTH; ,DA,,MOT) become the format-3 quantities ratio_to_rutherford and
ratio_to_mott (docs/upgrade/data_kinds.md sections 3.2, 4.3; step S14): ExforData.format3_quantity maps the code,
ExforData.to_dataset writes the ratios unconverted with the record, and to_azr refuses a ratio instead of
multiplying it by a Rutherford cross section (which froze sigma_R under an energy shift, did not fold it and used
mass numbers; audit item 8 had fixed its frame on 2026-10-26).

Fixture: tests/pyazr/data/exfor_master T0001004 (12C(p,el), EN/MEV, ANG lab, DATA NO-DIM ,,RTH) and
T0001005 (the same with ANG-CM).  Checks:
  1. elastic_lab_to_cm against the engine's own formulae (Unit::ConvertLabAngle / ConvertCrossSection);
  2. format3_quantity: RTH -> ratio_to_rutherford, MOT of identical particles -> ratio_to_mott; refused (ValueError,
     with the reason): RTH of identical particles (K2), MOT of different ones (point-nucleus Mott), an inelastic group;
  3. to_azr refuses a ratio; a no-dim differential that is not one is an analyzing power, written through;
  4. to_dataset: the ratios and errors unconverted, lab or c.m. angles as the heading says, the record's fields; no
     ERR-SYS and no norm= is refused;
  5. with the compiled engine (skipped without it): the code table equals the exfor column of the quantity table, and
     the written record runs in a format-3 p + 12C project: the data are the ratios, the angles the c.m. ones.

Run from anywhere:  python3 tests/pyazr/nds_rutherford_test.py
"""
import dataclasses
import importlib.util
import math
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
os.environ["AZURE2_EXFOR_MASTER"] = os.path.join(HERE, "data", "exfor_master")
spec = importlib.util.spec_from_file_location("pyazr_nds_standalone",
                                              os.path.join(ROOT, "pyazr", "sources", "nds.py"))
N = importlib.util.module_from_spec(spec)
spec.loader.exec_module(N)
import numpy as np  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def close(a, b, tol=1e-12):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return a.shape == b.shape and np.all(np.abs(a - b) <= tol * np.maximum(1.0, np.abs(b)))


def refused(f, *words):
    try:
        f()
    except ValueError as e:
        return all(w in str(e) for w in words)
    return False


MP, MT = 1.00782, 12.0          # the masses of the p + 12C test model
Z1, Z2 = 1, 6

print("1. elastic lab -> c.m. kinematics against the engine's formulae")
lab = np.array([0.0, 10.0, 30.0, 60.0, 90.0, 120.0, 150.0, 170.0, 180.0])
tc, jac = N.elastic_lab_to_cm(lab, MT, MP)
g = MP / MT
t = np.radians(lab)
tc_engine = np.degrees(t + np.arcsin(g * np.sin(t)))                      # Unit::ConvertLabAngle(pPair)
check("theta_cm = theta_lab + asin(gamma sin theta_lab)", close(tc, tc_engine), (tc, tc_engine))
inner = (lab > 0) & (lab < 180)
tcr = np.radians(tc[inner])
f_engine = np.sin(t[inner]) ** 2 / np.sin(tcr) ** 2 * np.cos(tcr - t[inner])  # Unit::ConvertCrossSection
check("J = sigma_cm/sigma_lab equals the engine's sin^2 cos form away from 0/180", close(jac[inner], f_engine, 1e-12),
      (jac[inner], f_engine))
check("J finite at 0 and 180 degrees, (1+gamma)^-2 and (1-gamma)^-2",
      close(jac[[0, -1]], [1 / (1 + g) ** 2, 1 / (1 - g) ** 2]), jac[[0, -1]])
check("theta_lab = atan2(sin theta_cm, gamma + cos theta_cm)",
      close(np.degrees(np.arctan2(np.sin(tcr), g + np.cos(tcr))), lab[inner], 1e-12))

print("\n2. the quantity code")
cm = N.fetch_exfor("T0001005", source="local")
lb = N.fetch_exfor("T0001004", source="local")
check("12C(p,el) ,,DA,,RTH -> ratio_to_rutherford", lb.format3_quantity() == "ratio_to_rutherford")
check("PAR,DA,,RTH -> ratio_to_rutherford",
      dataclasses.replace(lb, reaction="6-C-12(P,EL)6-C-12,PAR,DA,,RTH").format3_quantity() == "ratio_to_rutherford")
aa_mot = dataclasses.replace(lb, reaction="2-HE-4(A,EL)2-HE-4,,DA,,MOT")
check("4He(a,el) ,,DA,,MOT -> ratio_to_mott", aa_mot.format3_quantity() == "ratio_to_mott")
check("RTH of identical particles: refused, Rutherford or Mott is asked (K2)",
      refused(dataclasses.replace(lb, reaction="2-HE-4(A,EL)2-HE-4,,DA,,RTH").format3_quantity, "K2", "ratio_to_mott"))
check("MOT of different particles: refused (relativistic point-nucleus Mott)",
      refused(dataclasses.replace(lb, reaction="6-C-12(P,EL)6-C-12,,DA,,MOT").format3_quantity, "point-nucleus"))
check("an inelastic group: refused", refused(dataclasses.replace(lb, reaction="6-C-12(A,INL)6-C-12,,DA,,RTH").format3_quantity,
                                               "not elastic"))
check("a plain dsigma/dOmega is no Coulomb ratio", dataclasses.replace(lb, reaction="6-C-12(P,EL)6-C-12,,DA").format3_quantity() is None)

print("\n3. to_azr")
work = tempfile.mkdtemp(prefix="pyazr_nds_rth_")
try:
    check("to_azr refuses a ratio to Rutherford, naming to_dataset",
          refused(lambda: lb.to_azr(os.path.join(work, "data"), 1, 1, target_mass=MT, projectile_mass=MP), "to_dataset"))
    ay = dataclasses.replace(lb, reaction="6-C-12(P,EL)6-C-12,,POL/DA,,ANA")
    kw = ay.to_azr(os.path.join(work, "data"), 1, 1, file_name="ay")
    check("an analyzing power: analyzing-power, values through", kw["observable"] == "analyzing-power" and
          close(np.loadtxt(kw["data_file"])[:, 2], [0.9, 1.2, 1.5]), kw["observable"])
    check("the observable guess of the capture fixture: total-capture",
          N.fetch_exfor("T0001002", source="local")._guess_observable() == "total-capture")

    print("\n4. to_dataset")
    check("no ERR-SYS and no norm=: refused", refused(lambda: lb.to_dataset(os.path.join(work, "data"), "rth"), "norm="))
    rec_lb = lb.to_dataset(os.path.join(work, "data"), "rth_lab", norm="1 fit prior 1 3%", target_mass=MT,
                           projectile_mass=MP, file_name="lab")
    rec_cm = cm.to_dataset(os.path.join(work, "data"), "rth_cm", norm="1 fixed", target_mass=MT, projectile_mass=MP,
                           file_name="cm")
    rows_lb = np.loadtxt(os.path.join(work, "data", "lab.dat"))
    rows_cm = np.loadtxt(os.path.join(work, "data", "cm.dat"))
    ratio, err, angles = [0.9, 1.2, 1.5], [0.02, 0.03, 0.04], [30.0, 90.0, 150.0]
    check("lab and c.m. files: E_lab 2 MeV, the angles and the ratios and errors as given",
          close(rows_lb, np.column_stack([[2.0] * 3, angles, ratio, err])) and close(rows_cm, rows_lb), (rows_lb, rows_cm))
    want = ["dataset rth_lab", "reaction       12C(p,p)12C", "quantity       ratio_to_rutherford", "unit           1",
            "columns        E_lab[MeV] theta_lab[deg] value error", "norm           1 fit prior 1 3%", "exfor          T0001004"]
    check("the record of the lab set", all(w in rec_lb for w in want), rec_lb)
    check("the record of the c.m. set: theta_cm", "columns        E_lab[MeV] theta_cm[deg] value error" in rec_cm, rec_cm)
    check("identical particles with quantity='ratio_to_mott' given: 4He(a,a)4He",
          "reaction       4He(a,a)4He" in dataclasses.replace(lb, reaction="2-HE-4(A,EL)2-HE-4,,DA,,RTH").to_dataset(
              os.path.join(work, "data"), "aa", quantity="ratio_to_mott", norm="1 fixed", file_name="aa"))

    print("\n5. with the engine")
    sys.path.insert(0, ROOT)
    try:
        import pyazr
        from pyazr import Session
    except Exception as e:
        print(f"  skip  engine not importable ({e})")
    else:
        q = pyazr.schema()["quantities"]
        table = {code: name for name in ("ratio_to_rutherford", "ratio_to_mott") for code in q[name]["exfor"]}
        check("the code table = the exfor column of the two quantities",
              table == N._exfor_coulomb().COULOMB_RATIO_CODES, table)
        model = open(os.path.join(ROOT, "tests", "convolution_unify", "models", "fixed.azr")).read()
        head = model[:model.index("</config>")] + "format=3\n" + model[model.index("</config>"):model.index("</levels>") + 10]
        os.makedirs(os.path.join(work, "output"), exist_ok=True)
        os.makedirs(os.path.join(work, "checks"), exist_ok=True)
        azr = os.path.join(work, "rth.azr")
        with open(azr, "w") as f:
            f.write(head + "\n<data>\n" + rec_lb + "</data>\n<extrapolations>\n</extrapolations>\n<targets>\n</targets>\n")
        with Session(azr, cwd=work) as m:
            ang = np.asarray(m.angles[0], float)
            xs = np.asarray(m.cross[0], float)
        tcm, _ = N.elastic_lab_to_cm(np.array(angles), MT, MP)
        check("engine: lab angles converted to theta_cm", close(ang, tcm, 1e-8), (ang, tcm))
        check("engine: the data are the ratios, not converted", close(xs, ratio, 1e-12), xs)
finally:
    shutil.rmtree(work, ignore_errors=True)

print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
sys.exit(1 if failures else 0)
