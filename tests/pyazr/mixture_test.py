#!/usr/bin/env python3
"""Mixture segments (natural-abundance targets) as pyazr sees them.

A segment line may end with ``mixture w_self n {"file" w}...`` (MixtureTable.h):
the model becomes w_self * (this project) + the folded tables.  pyazr must read
and write the block, build tables with extrap_to_lab_table, and -- with the
engine -- return the mixture values and derivatives scaled by w_self.

Run from anywhere:  python3 tests/pyazr/mixture_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:                                           # pyazr itself needs numpy
    import numpy as np                         # noqa: E402
    from pyazr.datasets import SegmentSet, TestSegmentSet      # noqa: E402
    from pyazr.azrfile import AzrModel         # noqa: E402
    from pyazr.mixture import (pair_masses, extrap_to_lab_table, write_table,
                               lab_from_cm)    # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(0)

MODELS = os.path.join(ROOT, "tests", "mixture", "models")
FIXED = os.path.join(ROOT, "tests", "convolution_unify", "models", "fixed.azr")   # plain parent, for the parsing checks
GAUSS = os.path.join(MODELS, "gauss.azr")   # the same model with broad resonances and the 0.35-1.2 MeV grid, for the engine checks

print("1. the block is read from data and test lines")
m = AzrModel.from_file(os.path.join(MODELS, "mix_data.azr"))
s = SegmentSet.from_file(os.path.join(MODELS, "mix_data.azr"))
check("data line: w_self 0.5, one table", s[0].mixture_self == 0.5 and s[0].mixture_tables == (("table.dat", 0.5),),
      (s[0].mixture_self, s[0].mixture_tables))
check("data file still parsed", s[0].data_file == "data.dat", s[0].data_file)
t = TestSegmentSet.from_file(os.path.join(MODELS, "mix_two.azr"))
check("test line: two tables", t[0].mixture_self == 0.5 and len(t[0].mixture_tables) == 2
      and t[0].mixture_tables[1] == ("table.dat", 0.25), t[0].mixture_tables)
check("a plain model has no block", SegmentSet.from_file(FIXED) == [] and
      TestSegmentSet.from_file(FIXED)[0].mixture_self is None)
check("byte-identical re-write", m.to_text() == open(os.path.join(MODELS, "mix_data.azr")).read())

print("\n2. editing")
m2 = AzrModel.from_file(FIXED)
m2.set_test_mixture(1, 0.9048, [("tables/22 Ne.dat", 0.0925), ("tables/21Ne.dat", 0.0027)])
line = [l for l in m2.to_text().splitlines() if "mixture" in l]
check("set_test_mixture writes the block last, paths quoted", len(line) == 1 and
      line[0].endswith('mixture 0.9048 2 "tables/22 Ne.dat" 0.0925 "tables/21Ne.dat" 0.0027'), line)
d = tempfile.mkdtemp(prefix="pyazr_mixture_")
try:
    p = os.path.join(d, "e.azr"); m2.write(p)
    ts = TestSegmentSet.from_file(p)
    check("a quoted path with a space round-trips", ts[0].mixture_tables[0] == ("tables/22 Ne.dat", 0.0925), ts[0].mixture_tables)
    m2.add_data_segment("data/nat.dat", 1, 1, observable="angle-integrated", energy_min=0.1, energy_max=2,
                        energy_shift_sqrt=0.001, mixture=(0.9, [("t.dat", 0.1)]))
    m2.write(p); ss = SegmentSet.from_file(p)
    check("add_data_segment: sqrtshift and mixture together", ss[0].energy_shift_sqrt == 0.001 and ss[0].mixture_self == 0.9
          and ss[0].mixture_tables == (("t.dat", 0.1),), (ss[0].energy_shift_sqrt, ss[0].mixture_self, ss[0].mixture_tables))
    m2.set_mixture("nat.dat", 0.8, [("u.dat", 0.2)]); m2.write(p); ss = SegmentSet.from_file(p)
    check("set_mixture replaces the block and keeps sqrtshift", ss[0].mixture_tables == (("u.dat", 0.2),) and ss[0].energy_shift_sqrt == 0.001)
    m2.remove_mixture(1); m2.write(p); ss = SegmentSet.from_file(p)
    check("remove_mixture", ss[0].mixture_self is None and ss[0].energy_shift_sqrt == 0.001)
    try:
        m2.set_mixture("nowhere.dat", 1.0); check("unknown segment raises", False)
    except KeyError:
        check("unknown segment raises", True)
finally:
    shutil.rmtree(d, ignore_errors=True)

print("\n2b. Breit-Wigner reconstruction of an ENDF resonance range")
from pyazr.mixture import breit_wigner_total, _penetrability_shift_phase
# one isolated s-wave resonance with Gamma_n = Gamma on a spin-0 target: at E_r the total is
# 4 pi/k^2 [sin^2(phi) + g cos(2 phi)] with g = 1 (ENDF-102 D.1.1), and far away the hard sphere alone
awri, ap, er, gn = 20.0, 0.50, 1.0e5, 2.0e3
blocks = [(awri, 0.0, 0, 0, np.array([[er, 0.5, gn, gn, 0.0, 0.0]]))]
k = 2.196771e-3 * awri / (awri + 1.0) * np.sqrt(er)
phi = k * ap
sig = breit_wigner_total(np.array([er, 2.0e6]), 0.0, ap, blocks, naps=1)
want = 4 * np.pi / k ** 2 * (np.sin(phi) ** 2 + np.cos(2 * phi))
check("isolated s-wave resonance: peak total = 4pi/k^2 [sin^2 phi + g cos 2phi]", abs(sig[0] / want - 1) < 1e-9, (sig[0], want))
k2 = 2.196771e-3 * awri / (awri + 1.0) * np.sqrt(2.0e6)
check("far from it: hard sphere + a small tail", abs(sig[1] - 4 * np.pi / k2 ** 2 * np.sin(k2 * ap) ** 2) < 0.05 * sig[1], sig[1])
P, S, ph = _penetrability_shift_phase(2, 1.3)
check("d-wave penetrability / shift / phase (ENDF-102 Table D.1)",
      abs(P - 1.3 ** 5 / (9 + 3 * 1.69 + 1.69 ** 2)) < 1e-12 and abs(S + (18 + 3 * 1.69) / (9 + 3 * 1.69 + 1.69 ** 2)) < 1e-12
      and abs(ph - (1.3 - np.arctan(3 * 1.3 / (3 - 1.69)))) < 1e-12)

print("\n3. with the engine: values and derivatives")
try:
    from pyazr import azure2
    import pyazr._azure2  # noqa: F401
except Exception as err:
    print(f"skip (engine not built: {type(err).__name__}: {err})")
    sys.exit(1 if failures else 0)

d = tempfile.mkdtemp(prefix="pyazr_mixture_engine_")
try:
    def place(name, model):
        sub = os.path.join(d, name); os.makedirs(os.path.join(sub, "output")); os.makedirs(os.path.join(sub, "checks"))
        path = os.path.join(sub, "m.azr"); model.write(path); return path

    base = AzrModel.from_file(GAUSS)
    m1, m2mass = pair_masses(base, 1)
    check("pair masses of p+12C", abs(m1 - 1.00782) < 1e-9 and abs(m2mass - 12.0) < 1e-9, (m1, m2mass))

    # the table: this model's bare cross section, c.m. differential at 120 deg c.m., 0.2-1.8 MeV lab
    bare = AzrModel.from_file(GAUSS).clear_target_effects().clear_extrapolations()
    bare.add_extrapolation(1, 1, 0.2, 1.8, 0.0001, observable="differential-cm", angle=120)
    with azure2(place("bare", bare), data_mode=False, use_brune=False, ignore_externals=False) as mb:
        x0 = np.array(mb.params_rwa)
        e_cm = np.array(mb.calculate_energies(x0)[0]); sig = np.array(mb.calculate_rwa(x0)[0])
    table = os.path.join(d, "table.dat")
    write_table(lab_from_cm(e_cm, m1, m2mass), sig, table, header="bare p+12C 120 deg c.m.")
    check("table written", os.path.getsize(table) > 0 and len(e_cm) > 15000, len(e_cm))

    # the reference: the folded model on the data grid (same kernel as fixed.azr), c.m. differential
    ref = AzrModel.from_file(GAUSS).clear_extrapolations()
    ref.add_extrapolation(1, 1, 0.35, 1.2, 0.002, observable="differential-cm", angle=120)
    with azure2(place("ref", ref), data_mode=False, use_brune=False, ignore_externals=False) as mr:
        e_ref = np.array(mr.calculate_energies(x0)[0]); f_ref = np.array(mr.calculate_rwa(x0)[0])
    data = os.path.join(d, "data.dat")
    np.savetxt(data, np.column_stack([lab_from_cm(e_ref, m1, m2mass), np.full_like(e_ref, 120.0), f_ref, 0.01 * f_ref]))

    def data_model(mixture):
        mm = AzrModel.from_file(GAUSS).clear_extrapolations()
        mm.add_data_segment(data, 1, 1, observable="differential-cm", energy_min=0.3, energy_max=1.3,
                            angle_min=0, angle_max=180, norm_error=1.0, mixture=mixture)
        return mm

    with azure2(place("plain", data_model(None)), use_brune=False, ignore_externals=False) as mp, \
         azure2(place("mix", data_model((0.5, [(table, 0.5)]))), use_brune=False, ignore_externals=False) as mx, \
         azure2(place("half", data_model((0.5, []))), use_brune=False, ignore_externals=False) as mh:
        xp = np.array(mp.params_rwa)
        check("same free parameters in every variant", len(xp) == len(mx.params_rwa) == len(mh.params_rwa))
        fp = np.array(mp.calculate_rwa(xp)[0]); fx = np.array(mx.calculate_rwa(xp)[0]); fh = np.array(mh.calculate_rwa(xp)[0])
        check("plain data model reproduces its data", np.max(np.abs(fp / f_ref - 1)) < 1e-9, np.max(np.abs(fp / f_ref - 1)))
        check("half the model + half the table == the folded model (1e-4)", np.max(np.abs(fx / fp - 1)) < 1e-4, np.max(np.abs(fx / fp - 1)))
        check("mixture 0.5 with no table == half the model (1e-12)", np.max(np.abs(fh / (0.5 * fp) - 1)) < 1e-12, np.max(np.abs(fh / (0.5 * fp) - 1)))
        chi2 = float(np.sum(np.asarray(mx.calculate_chi2_rwa(xp), dtype=float)))
        check("chi2 of the mixture ~ 0 (table interpolation only)", chi2 < 1e-2, chi2)
        rp, Jp = mp.residual_jacobian(xp); rx, Jx = mx.residual_jacobian(xp); rh, Jh = mh.residual_jacobian(xp)
        n = mp.n_rmatrix
        check("Jacobian R-matrix columns of 'half' == 0.5 x plain (1e-10)",
              np.max(np.abs(Jh[:, :n] - 0.5 * Jp[:, :n])) < 1e-10 * max(1.0, np.max(np.abs(Jp[:, :n]))),
              np.max(np.abs(Jh[:, :n] - 0.5 * Jp[:, :n])))
        check("Jacobian R-matrix columns of the mixture == 0.5 x plain (1e-10)",
              np.max(np.abs(Jx[:, :n] - 0.5 * Jp[:, :n])) < 1e-10 * max(1.0, np.max(np.abs(Jp[:, :n]))),
              np.max(np.abs(Jx[:, :n] - 0.5 * Jp[:, :n])))
        err = 0.01 * f_ref
        check("residuals of 'half' = (0.5 f - y)/e", np.max(np.abs(rh - (0.5 * fp - f_ref) / err)) < 1e-8, np.max(np.abs(rh - (0.5 * fp - f_ref) / err)))
        # central finite differences on three R-matrix columns: the mixture model must agree with
        # its finite differences as well as the plain model does (the FD truncation error on a
        # resonant model is the same for both; what is tested is that the table adds nothing and
        # the w_self scaling is right)
        def fd_mismatch(m, J):
            worst = 0.0
            for col in sorted(set([0, n // 2, n - 1])):
                h = 1e-4 * max(1.0, abs(xp[col]))
                xpl = xp.copy(); xpl[col] += h; xmi = xp.copy(); xmi[col] -= h
                fd = (m.residual_jacobian(xpl)[0] - m.residual_jacobian(xmi)[0]) / (2 * h)
                worst = max(worst, np.max(np.abs(fd - J[:, col])) / max(np.max(np.abs(J[:, col])), 1e-12))
            return worst
        wp, wx = fd_mismatch(mp, Jp), fd_mismatch(mx, Jx)
        check("analytic Jacobian vs finite differences: mixture no worse than the plain model",
              wx < max(1e-5, 1.5 * wp), (wx, wp))
finally:
    shutil.rmtree(d, ignore_errors=True)

print("\n" + ("FAILED: " + ", ".join(failures) if failures else "PASSED"))
sys.exit(1 if failures else 0)
