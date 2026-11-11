#!/usr/bin/env python3
"""Ratio-to-Rutherford EXFOR data with a LAB angle column (audit item 8, fixed 2026-10-26).

The ratio sigma/sigma_R is the same in both frames, but to_azr writes a cross section in the frame of
the angle column.  For an ``ANG`` (lab) column it used to multiply by the c.m. Rutherford cross section
at the lab angle: 1.13x off for p + 16O, 1.8-2.2x for alpha + 12C.  Now the lab angle goes to the c.m.
frame with the elastic kinematics (``elastic_lab_to_cm``), and the c.m. value ratio x sigma_R(theta_cm)
is divided by the lab-to-c.m. Jacobian AZURE2 applies when it reads a lab file, so the engine sees
exactly ratio x sigma_R,cm(theta_cm).  ``_guess_observable`` also called every angle column c.m.; a
lab ``ANG`` column is ``differential`` now, and a ratio to Rutherford is no analyzing power.

Fixture: tests/pyazr/data/exfor_master T0001004 (12C(p,el), EN/MEV, ANG lab, DATA NO-DIM ,,RTH) and
T0001005 (the same with ANG-CM).  Checks:
  1. elastic_lab_to_cm against the engine's own formulae (EPoint::ConvertLabAngle / ConvertCrossSection);
  2. ANG-CM: values = ratio x sigma_R,cm(E_cm, theta);  ANG: values = ratio x sigma_R,cm(theta_cm) / J;
  3. the observable guess: lab -> differential, c.m. -> differential-cm, no-dim without RTH -> analyzing
     power, the capture fixture -> total-capture;
  4. with the compiled engine (skipped without it): the lab file through a p + 12C model, the engine's
     converted c.m. angles and data = theta_cm and ratio x sigma_R,cm(theta_cm) to 1e-8.

Run from anywhere:  python3 tests/pyazr/nds_rutherford_test.py
"""
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
spec = importlib.util.spec_from_file_location("pyazr_nds_standalone", os.path.join(ROOT, "pyazr", "nds.py"))
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


MP, MT = 1.00782, 12.0          # the masses of the p + 12C test model
Z1, Z2 = 1, 6

print("1. elastic lab -> c.m. kinematics against the engine's formulae")
lab = np.array([0.0, 10.0, 30.0, 60.0, 90.0, 120.0, 150.0, 170.0, 180.0])
tc, jac = N.elastic_lab_to_cm(lab, MT, MP)
g = MP / MT
t = np.radians(lab)
tc_engine = np.degrees(t + np.arcsin(g * np.sin(t)))                      # EPoint::ConvertLabAngle(pPair)
check("theta_cm = theta_lab + asin(gamma sin theta_lab)", close(tc, tc_engine), (tc, tc_engine))
inner = (lab > 0) & (lab < 180)
tcr = np.radians(tc[inner])
f_engine = np.sin(t[inner]) ** 2 / np.sin(tcr) ** 2 * np.cos(tcr - t[inner])  # EPoint::ConvertCrossSection
check("J = sigma_cm/sigma_lab equals the engine's sin^2 cos form away from 0/180", close(jac[inner], f_engine, 1e-12),
      (jac[inner], f_engine))
check("J finite at 0 and 180 degrees, (1+gamma)^-2 and (1-gamma)^-2",
      close(jac[[0, -1]], [1 / (1 + g) ** 2, 1 / (1 - g) ** 2]), jac[[0, -1]])
check("theta_lab = atan2(sin theta_cm, gamma + cos theta_cm)",
      close(np.degrees(np.arctan2(np.sin(tcr), g + np.cos(tcr))), lab[inner], 1e-12))

print("\n2. the converted values")
cm = N.fetch_exfor("T0001005", source="local")
lb = N.fetch_exfor("T0001004", source="local")
ratio = np.array([0.9, 1.2, 1.5])
err = np.array([0.02, 0.03, 0.04])
angles = np.array([30.0, 90.0, 150.0])
e_cm = 2.0 * MT / (MT + MP)
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    rows_cm = N._azr_columns(cm, target_mass=MT, projectile_mass=MP)
    rows_lb = N._azr_columns(lb, target_mass=MT, projectile_mass=MP)
check("no warning for a ,,RTH reaction code", not w, [str(x.message) for x in w])
r_cm = N._rutherford_barn_per_sr(Z1, Z2, e_cm, angles)
check("ANG-CM: ratio x sigma_R,cm(theta)", close(rows_cm[:, 2], ratio * r_cm) and close(rows_cm[:, 3], err * r_cm),
      rows_cm[:, 2:])
check("ANG-CM: the angle column is written through", close(rows_cm[:, 1], angles))
tcm, j = N.elastic_lab_to_cm(angles, MT, MP)
want = ratio * N._rutherford_barn_per_sr(Z1, Z2, e_cm, tcm) / j
check("ANG (lab): ratio x sigma_R,cm(theta_cm) / J", close(rows_lb[:, 2], want) and close(rows_lb[:, 3], err / ratio * want),
      (rows_lb[:, 2], want))
check("ANG (lab): the lab angle column is written through", close(rows_lb[:, 1], angles))
old = ratio * r_cm                                           # what was written before the fix
check("the old value (c.m. Rutherford at the lab angle) is 1.05-1.3x off here",
      np.all(old / want > 1.04) and np.all(old / want < 1.3), old / want)
# the Rutherford helper itself: (Z1 Z2 alpha hbar c / 4 E)^2 / sin^4(theta/2) in barn
amp = Z1 * Z2 * N._FSTRUC * N._HBARC / (4.0 * e_cm)
check("sigma_R,cm formula", close(r_cm, amp ** 2 / np.sin(np.radians(angles) / 2) ** 4 * 1e-2))

print("\n3. the observable guess")
check("lab ANG ratio -> differential", lb._guess_observable(rutherford=True) == "differential")
check("ANG-CM ratio -> differential-cm", cm._guess_observable(rutherford=True) == "differential-cm")
check("a no-dim differential that is not a ratio -> analyzing-power",
      lb._guess_observable() == "analyzing-power" and cm._guess_observable() == "analyzing-power")
cap = N.fetch_exfor("T0001002", source="local")
check("the capture fixture -> total-capture", cap._guess_observable() == "total-capture")
work = tempfile.mkdtemp(prefix="pyazr_nds_rth_")
try:
    kw_lb = lb.to_azr(os.path.join(work, "data"), 1, 1, target_mass=MT, projectile_mass=MP)
    kw_cm = cm.to_azr(os.path.join(work, "data"), 1, 1, target_mass=MT, projectile_mass=MP, file_name="cm")
    check("to_azr: lab set -> observable differential, c.m. set -> differential-cm",
          kw_lb["observable"] == "differential" and kw_cm["observable"] == "differential-cm",
          (kw_lb["observable"], kw_cm["observable"]))
    check("to_azr: angle range from the lab angles", (kw_lb["angle_min"], kw_lb["angle_max"]) == (30.0, 150.0))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        kw_forced = lb.to_azr(os.path.join(work, "data"), 1, 1, rutherford=False, file_name="ay")
    check("rutherford=False: values through, analyzing-power", kw_forced["observable"] == "analyzing-power"
          and close(np.loadtxt(kw_forced["data_file"])[:, 2], ratio), kw_forced["observable"])

    print("\n4. the lab file through the engine")
    sys.path.insert(0, ROOT)
    try:
        from pyazr import azure2, AzrModel
        import pyazr._azure2  # noqa: F401
    except Exception as e:
        print(f"  skip  engine not importable ({e})")
    else:
        model = AzrModel.from_file(os.path.join(ROOT, "tests", "convolution_unify", "models", "fixed.azr"))
        model.add_data_segment(**kw_lb)
        azr = os.path.join(work, "rth.azr")
        os.makedirs(os.path.join(work, "output"), exist_ok=True)
        os.makedirs(os.path.join(work, "checks"), exist_ok=True)
        model.write(azr)
        with azure2(azr, cwd=work) as m:
            ang = np.asarray(m.angles[0], float)
            xs = np.asarray(m.cross[0], float)
            en = np.asarray(m.energies[0], float)
        check("engine: lab angles converted to theta_cm", close(ang, tcm, 1e-8), (ang, tcm))
        check("engine: c.m. energy", close(en, [e_cm] * 3, 1e-8), en)
        want_cm = ratio * N._rutherford_barn_per_sr(Z1, Z2, e_cm, tcm)
        # the data file carries 8 significant digits
        check("engine: converted data = ratio x sigma_R,cm(theta_cm)", close(xs, want_cm, 1e-7), (xs, want_cm))
finally:
    shutil.rmtree(work, ignore_errors=True)

print("\n" + ("ALL PASSED" if not failures else f"FAILED: {failures}"))
sys.exit(1 if failures else 0)
