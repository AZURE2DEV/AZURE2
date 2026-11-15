#!/usr/bin/env python3
"""The engine's door of the nuclear amplitude expansion (ENDF-102, 2023 edition, LAW=5 LTP=1, eqs. 6.9-6.14,
pp. 145-147; docs/upgrade/endf.md section 13.12): the point output ``nuclear_amplitudes`` of an elastic data set.

(a) For T1 (p + 12C, spin 1/2 on 0), identical_pp_res (p + p, identical fermions) and identical_dd_mott (d + d,
    identical bosons of spin 1): the elastic dsigma/dOmega rebuilt from the coefficients by the manual's formula,
    written out below from the page (not the module's), equals the engine's dsigma/dOmega at 20 lab energies x 12
    c.m. angles to 1e-9 (the constants of eqs. 6.11-6.12 here are CODATA 2018 and differ from the engine's in the
    last digits: measured 1e-11).
(b) The list has the manual's length: NW = 4 NL + 3 (distinguishable), 3 NL + 3 (identical).
(c) The module's reader of LTP=1 (pyazr.endf.distributions.amplitude_expansion) gives the same sum as (a).
(d) A data set with a target-effect fold or of another reaction carries no coefficients (an empty list per point).
Exits 77 (skipped) without the compiled engine module.  Run from anywhere:  python3 tests/pyazr/endf_amplitudes_test.py
"""
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    import numpy as np
    from numpy.polynomial import legendre
    import pyazr
    from pyazr import _azure2  # noqa: F401
    from pyazr.endf.distributions import amplitude_expansion
    from pyazr.endf.writer import _pair_keys
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)

warnings.simplefilter("ignore")
ALPHA = 7.2973525693e-3      # CODATA 2018
HBARC = 197.3269804          # MeV fm
UC2 = 931.49410242           # MeV
ANGLES = np.array([10, 20, 30, 45, 60, 75, 90, 105, 120, 135, 150, 170.0])
TMP = tempfile.mkdtemp(prefix="endf_amp_")
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:600]}"))
    if not ok:
        failures.append(name)


def P(l, mu):
    return legendre.legval(mu, [0] * l + [1])


def manual(values, e_lab, m1, m2, z1, z2, mu, identical, two_s):
    """dsigma/dOmega (b/sr) of one LTP=1 record by eqs. 6.9-6.14 (p. 145), E the lab energy (MeV), masses in u."""
    A = m2 / m1
    k = A / (1 + A) * np.sqrt(2 * UC2 * m1 * e_lab) / HBARC * 10.0       # eq. 6.11, b^-1/2
    eta = z1 * z2 * ALPHA * np.sqrt(UC2 * m1 / (2 * e_lab))              # eq. 6.12
    v = np.asarray(values, float)
    nl = (len(v) - 3) // (3 if identical else 4)
    b = v[:nl + 1] if identical else v[:2 * nl + 1]
    a = v[len(b):][0::2] + 1j * v[len(b):][1::2]
    if not identical:                                                     # eqs. 6.9, 6.13
        sc = eta ** 2 / (k ** 2 * (1 - mu) ** 2)
        s = sum((2 * l + 1) / 2 * a[l] * P(l, mu) for l in range(nl + 1))
        it = -2 * eta / (1 - mu) * np.real(np.exp(1j * eta * np.log((1 - mu) / 2)) * s)
        return sc + it + sum((2 * l + 1) / 2 * b[l] * P(l, mu) for l in range(2 * nl + 1)), sc
    s2 = two_s / 2                                                        # eqs. 6.10, 6.14
    sc = 2 * eta ** 2 / (k ** 2 * (1 - mu ** 2)) * ((1 + mu ** 2) / (1 - mu ** 2) + (-1) ** (2 * s2) / (2 * s2 + 1)
                                                     * np.cos(eta * np.log((1 + mu) / (1 - mu))))
    s = sum(((1 + mu) * np.exp(1j * eta * np.log((1 - mu) / 2)) + (-1) ** l * (1 - mu)
             * np.exp(1j * eta * np.log((1 + mu) / 2))) * (2 * l + 1) / 2 * a[l] * P(l, mu) for l in range(nl + 1))
    nuc = sum((4 * l + 1) / 2 * b[l] * P(2 * l, mu) for l in range(nl + 1))
    return sc - 2 * eta / (1 - mu ** 2) * np.real(s) + nuc, sc


def probe(project_dir, azr, n_energies=20):
    """The engine's dsigma/dOmega and coefficients of the project's elastic channel on a data set of placeholder
    values at n_energies lab energies over the data's range x ANGLES."""
    d = os.path.join(TMP, os.path.basename(project_dir))
    shutil.copytree(project_dir, d, ignore=shutil.ignore_patterns("output*", "checks*"))
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)
    s0 = pyazr.Session(os.path.join(d, azr))
    ent = next(p.to_dict() for p in s0.pairs if p.is_entrance)
    lab = [x for dd in s0.evaluate(s0.parameters.x0).data for x in np.asarray(dd.extra.get("energy_lab", []), float)]
    energies = np.linspace(min(lab), max(lab), n_energies)
    proj = pyazr.Project.open(os.path.join(d, azr))
    p = pyazr.Project.new(os.path.join(d, "probe.azr"), nuclear_part=s0.file)
    np.savetxt(os.path.join(d, "probe.dat"), [(e, a, 1.0, 1.0) for e in energies for a in ANGLES], fmt="%.17g")
    key = ent["key"]
    p.add_dataset("el", reaction=proj.reaction(key, key), quantity="dsigma_domega", unit="b/sr",
                  file=os.path.join(d, "probe.dat"), columns="E_lab[MeV] theta_cm[deg] value error",
                  norm=pyazr.Nuisance(1, "fixed"), **_pair_keys(proj, key, key))
    p.save()
    s0.close()
    s = pyazr.Session(p.path)
    try:
        dd = s.evaluate(s.parameters.x0).data[0]
        return ent, dd
    finally:
        s.close()


CASES = [("upgrade/T1", "T1.azr"), ("identical_pp_res", "identical_pp_res.azr"),
         ("identical_dd_mott", "identical_dd_mott.azr")]
for folder, azr in CASES:
    print(azr)
    ent, d = probe(os.path.join(ROOT, "tests", folder), azr)
    identical = ent["Z1"] == ent["Z2"] and round(ent["M1"]) == round(ent["M2"])
    two_s = int(round(2 * ent["J1"]))
    model = np.asarray(d.model, float)
    elab = np.asarray(d.extra["energy_lab"], float)
    mu = np.cos(np.radians(np.asarray(d.angles, float)))
    coeffs = d.extra["nuclear_amplitudes"]
    got = [manual(c, e, ent["M1"], ent["M2"], ent["Z1"], ent["Z2"], m, identical, two_s) for c, e, m in
           zip(coeffs, elab, mu)]
    rebuilt = np.array([g[0] for g in got])
    rel = np.abs(rebuilt - model) / np.abs(model)
    print(f"      {len(model)} points: rebuilt by eqs. 6.9-6.14 against the engine, largest relative {rel.max():.2e}, "
          f"median {np.median(rel):.1e}")
    check(f"{azr}: dsigma/dOmega from the coefficients (eqs. 6.9-6.14) = the engine's to 1e-9", rel.max() < 1e-9,
          rel.max())
    lengths = {len(c) for c in coeffs}
    nl = [(n - 3) / (3 if identical else 4) for n in lengths]
    check(f"{azr}: every list has the manual's length NW = {'3' if identical else '4'} NL + 3 (lengths {lengths})",
          all(x == int(x) and x >= 0 for x in nl), lengths)
    # (c) the module's reader of LTP=1 with the same eta
    eta = ent["Z1"] * ent["Z2"] * ALPHA * np.sqrt(UC2 * ent["M1"] / (2 * elab))
    module = np.array([amplitude_expansion(c, int((len(c) - 3) // (3 if identical else 4)), identical, h, m)
                       for c, h, m in zip(coeffs, eta, mu)])
    own = rebuilt - np.array([g[1] for g in got])
    dev = np.max(np.abs(module - own)) / np.max(np.abs(model))
    check(f"{azr}: the module's LTP=1 sum (amplitude_expansion) = the test's own lines to 1e-13 of the cross "
          "section", dev < 1e-13, dev)

# (d) a data set of another reaction or with a target-effect fold carries no coefficients
print("other data sets")
d = os.path.join(TMP, "T1_own")
shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), d, ignore=shutil.ignore_patterns("output*", "checks*"))
for sub in ("output", "checks"):
    os.makedirs(os.path.join(d, sub), exist_ok=True)
s = pyazr.Session(os.path.join(d, "T1.azr"))
ev = s.evaluate(s.parameters.x0)
kinds = {dd.name: (dd.quantity, dd.extra.get("nuclear_amplitudes")) for dd in ev.data}
elastic = [n for n, (q, v) in kinds.items() if q == "dsigma_domega" and n.startswith("el_")]
check("T1: its elastic data sets carry coefficients at every point, its capture and A_y sets none or empty lists",
      elastic and all(len(kinds[n][1]) == len(ev.data[[x.name for x in ev.data].index(n)].model) and
                      all(len(v) > 0 for v in kinds[n][1]) for n in elastic) and
      all(v is None or all(len(x) == 0 for x in v) for n, (q, v) in kinds.items() if n not in elastic),
      {n: (q, None if v is None else len(v)) for n, (q, v) in kinds.items()})
s.close()

shutil.rmtree(TMP, ignore_errors=True)
print(f"{len(failures)} failure(s)" + (": " + ", ".join(failures) if failures else ""))
sys.exit(1 if failures else 0)
