#!/usr/bin/env python3
"""Model data through pyazr (<modelData>, src/ModelData.cpp): pyazr.model_data's records on T2's generator
parameters (tests/upgrade/T2/T2_true.azr) against the hand values of tests/model_data/README.md; the term in
penalties(); the term's share of chi2_and_grad against central differences; the binding's resonance_strength,
statistical_factor and branching_ratio against the hand values.

Needs the compiled module; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/model_data_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
T2 = os.path.join(ROOT, "tests", "upgrade", "T2")
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def close(a, b, tol):
    return abs(a - b) <= tol * max(abs(b), 1e-300)


try:
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    import pyazr
    from pyazr import azure2, _azure2  # noqa: F401
except Exception as err:
    print(f"skip: the engine is not importable ({type(err).__name__}: {err})")
    sys.exit(0)

MODEL_DATA = """<modelData>
datum wg_L5
    quantity   resonance_strength
    level      5
    J          3-
    reaction   15N(p,a)12C
    value      2.7318028 keV
    error      0.1 keV
    source     "hand value + 1 sigma"
datum G_L5
    quantity   total_width
    level      5
    J          3-
    value      15.622 keV
    error      1 keV
datum a_p15N
    quantity   scattering_length
    pair       p+15N
    s          1
    value      0 fm
    error      100 fm
</modelData>
"""

# The binding's physics functions (physics/ResonanceQuantities.h), on arrays.
w = _azure2.statistical_factor(np.array([6, 4]), np.array([1, 0]), np.array([1, 0]), np.array([False, True]))
check("statistical_factor: 7/4 (p + 15N, 3-), 10 (12C + 12C, 2+)", list(w) == [1.75, 10.0], w)
s = _azure2.resonance_strength(np.array([1.75]), np.array([1.702]), np.array([12.92]), np.array([14.622]))
check("resonance_strength: 2.6318028 keV", close(s["value"][0], 2.6318028, 5e-8), s)
check("resonance_strength: d/dGamma = -value/Gamma", close(s["d_gamma"][0], -s["value"][0] / 14.622, 1e-14), s)
b = _azure2.branching_ratio(np.array([12.92]), np.array([14.622]))
check("branching_ratio: 0.88360005", close(b["value"][0], 0.88360005, 1e-8), b)

tmp = tempfile.mkdtemp(prefix="pyazr_model_data_")
try:
    work = os.path.join(tmp, "t2")
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    shutil.copytree(os.path.join(T2, "data"), os.path.join(work, "data"))
    with open(os.path.join(T2, "T2_true.azr")) as f:
        text = f.read()
    path = os.path.join(work, "m.azr")
    with open(path, "w") as f:
        f.write(text + MODEL_DATA)
    m = azure2(path)
    x = np.asarray(m.params_rwa, float)
    records = pyazr.model_data(m)
    check("three records, in file order", [r.name for r in records] == ["wg_L5", "G_L5", "a_p15N"], records)
    by = {r.name: r for r in records}
    check("records are ModelDatum", all(isinstance(r, pyazr.ModelDatum) for r in records))
    check("omega gamma of level 5 = 2.6318028 keV, in keV", close(by["wg_L5"].model, 2.6318028, 5e-7) and by["wg_L5"].unit == "keV",
          by["wg_L5"])
    check("its residual = -1", close(by["wg_L5"].residual, -1.0, 5e-6), by["wg_L5"])
    check("Gamma of level 5 = 14.622 keV", close(by["G_L5"].model, 14.622, 1e-9), by["G_L5"])
    a0 = {(v.pair, v.jpi): v for v in m.scattering_lengths()}[(1, "1-")].a0.real
    check("scattering length = scattering_lengths()", close(by["a_p15N"].model, a0, 1e-12), (by["a_p15N"].model, a0))
    check("source kept", by["wg_L5"].source == "hand value + 1 sigma", by["wg_L5"].source)
    term = m.penalties(x)["scattering_length"]  # the whole model-data term
    check("penalties: the term = sum of the records' chi2", close(term, sum(r.chi2 for r in records), 1e-12),
          (term, [r.chi2 for r in records]))

    # The gradient with the term: a free alpha0 width and energy of level 5 (central differences of the objective).
    # Segment 5 alone (the implanted target across level 5), as tests/model_data/check.sh.
    fit_text, in_segments, n_seg = [], False, 0
    for line in text.splitlines():
        t = line.split()
        if line.startswith("<segmentsData>"):
            in_segments = True
        elif line.startswith("</segmentsData>"):
            in_segments = False
        elif in_segments and t:
            n_seg += 1
            if n_seg != 5:
                t[0] = "0"
                line = " ".join(t)
        if line.startswith("3 -1 13.269200 "):
            t[3] = "0"
            if t[5] == "2":
                t[10] = "0"
            line = " ".join(t)
        fit_text.append(line)
    # The term's share of the gradient: with the data minus without them (the data part's own finite-difference
    # error, about 5e-5 relative in the energy here, cancels).
    shares = {}
    for name, extra in (("with", MODEL_DATA), ("without", "")):
        fpath = os.path.join(work, f"fit_{name}.azr")
        with open(fpath, "w") as f:
            f.write("\n".join(fit_text) + "\n" + extra)
        fit = azure2(fpath)
        xf = np.asarray(fit.params_rwa, float)
        value, grad = fit.chi2_and_grad(xf)
        h = 1e-6 * np.maximum(1.0, np.abs(xf))
        obj = lambda v: float(np.sum(fit.calculate_chi2_rwa(v)))
        fd = np.array([(obj(xf + h[i] * e) - obj(xf - h[i] * e)) / (2 * h[i]) for i, e in enumerate(np.eye(xf.size))])
        shares[name] = (grad, fd)
    an = shares["with"][0] - shares["without"][0]
    fd = shares["with"][1] - shares["without"][1]
    worst = float(np.max(np.abs(an - fd) / np.maximum(1.0, np.abs(fd))))
    check(f"chi2_and_grad: the term's share = central differences ({worst:.2e})", worst < 1e-6, (an, fd))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL pyazr_model_data ({len(failures)})")
    sys.exit(1)
print("PASS pyazr_model_data")
