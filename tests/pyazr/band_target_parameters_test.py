#!/usr/bin/env python3
"""A fitted target parameter in the uncertainty band (docs/source/user_guide/fitting.rst, "Uncertainty bands"; the
owner's answer of 2026-11-18, FOR_THE_OWNER.md section 8).

Project: tests/layer_yield/models/fit.azr, a thick-target yield whose areal density is fitted beside the level
parameters (fit.azr: ``areal_density 1e18 atoms/cm2 fit``).  The program (AZURE2_BIN) fits it with
``--covariance-band``; then
(a) covariance.dat names the density among its columns (its "#columns" header) beside the level parameters;
(b) the band equals a direct propagation: sigma_i^2 = g_i^T C g_i with C = covariance.dat and g_i the central
    differences of the model of point i in every free parameter of C (pyazr's Session.evaluate at the fitted
    vector, steps 1e-5 (|x| + 1)), to 1e-4 of sigma_i at every point of the band file (the band's sensitivities of
    the levels are the adjoint's, its density column central differences of step 1e-6 (|x| + 1); the fitted vector
    comes back from param.sav's 8 digits, so the values agree to 1e-5); without the
    density in the band (the engine before 2026-11-18) the band is the level block's alone and misses this by the
    density's share;
(c) the same project with the density fixed at its fitted value (8 digits), fitted from the first fit's param.sav (the
    same curve to 1e-4): its band is
    nowhere wider than the first's (to 1e-3 of sigma, the noise of two numerical HESSE matrices) and narrower by more
    than 1 % somewhere: fixing a parameter removes its share sigma_full^2 - sigma_fixed^2 = (g . C e_d)^2 / C_dd >= 0.
(d) an extrapolation run of the project after the fit (no density parameter in that run) takes covariance.dat on its
    own column by name: its band equals, byte for byte, the band from a covariance.dat that holds that column alone.
Exits 77 (skipped) without the engine module or AZURE2_BIN.  About a minute at one thread.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    import numpy as np
    import pyazr
    from pyazr import _azure2  # noqa: F401
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)
BIN = os.environ.get("AZURE2_BIN")
if not BIN or not os.path.exists(BIN):
    print("skip: AZURE2_BIN not set")
    sys.exit(77)

SRC = os.path.join(ROOT, "tests", "layer_yield", "models")
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:600]}"))
    if not ok:
        failures.append(name)


def copy(dst, text=None):
    shutil.copytree(SRC, dst, ignore=shutil.ignore_patterns("output*", "checks*"))
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(dst, sub), exist_ok=True)
    if text is not None:
        with open(os.path.join(dst, "fit.azr"), "w") as f:
            f.write(text)
    return dst


def fit(d, *extra):
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "1"))
    p = subprocess.run([BIN, "fit", "--no-gui", "--covariance-band", *extra, "fit.azr"], cwd=d, env=env,
                       capture_output=True, text=True, timeout=1800)
    if p.returncode != 0:
        print(p.stdout[-2000:], p.stderr[-2000:])
    return p.returncode


def band_files(d):
    """{file: array of rows (E, Ex, angle, value, sigma, ...)} of the .band files of a run."""
    out = {}
    for name in sorted(os.listdir(os.path.join(d, "output"))):
        if name.endswith(".out.band"):
            rows = [[float(v) for v in line.split()] for line in open(os.path.join(d, "output", name)) if line.strip()]
            out[name] = np.array(rows)
    return out


def covariance(d):
    lines = open(os.path.join(d, "output", "covariance.dat")).read().splitlines()
    names = next((l.split()[1:] for l in lines if l.startswith("#columns")), [])
    rows = [[float(v) for v in l.split()] for l in lines if l.strip() and not l.startswith("#")]
    return names, np.array(rows)


tmp = tempfile.mkdtemp(prefix="band_target_")
try:
    free_dir = copy(os.path.join(tmp, "free"))
    check("the fit with the density free ends (exit 0)", fit(free_dir) == 0)
    names, C = covariance(free_dir)
    density = [n for n in names if n.endswith("areal_density")]
    levels = [n for n in names if re.match(r"^(energy|width)_", n)]
    print(f"      covariance.dat: {len(names)} columns: {names}")
    # (a)
    check("(a) covariance.dat names the fitted areal density beside the level parameters",
          len(density) == 1 and len(levels) >= 1 and C.shape == (len(names), len(names)), names)

    # (b) direct propagation with pyazr: the model of every data point at x +- h e_k
    os.chdir(free_dir)
    with pyazr.Session("fit.azr", parameter_file=os.path.join(free_dir, "output", "param.sav")) as s:
        x = np.asarray(s.parameters.x0, float)
        free_names = [p.name for p in s.parameters if not p.fixed]
        idx = [free_names.index(n) for n in names]

        def models(v):
            return np.concatenate([np.asarray(d_.model, float) for d_ in s.evaluate(v).data])

        m0 = models(x)
        G = np.zeros((m0.size, len(names)))
        for c, k in enumerate(idx):
            h = 1e-5 * (abs(x[k]) + 1.0)
            xp, xm = x.copy(), x.copy()
            xp[k] += h
            xm[k] -= h
            G[:, c] = (models(xp) - models(xm)) / (2 * h)
        sigma_direct = np.sqrt(np.einsum("ia,ab,ib->i", G, C, G))
        d_col = names.index(density[0]) if density else 0
        keep = [c for c in range(len(names)) if c != d_col]
        sigma_levels = np.sqrt(np.einsum("ia,ab,ib->i", G[:, keep], C[np.ix_(keep, keep)], G[:, keep]))
    os.chdir(HERE)
    bands = band_files(free_dir)
    band = np.concatenate([b for b in bands.values()]) if bands else np.zeros((0, 5))
    n = min(len(band), m0.size)
    ok_shape = len(band) == m0.size
    value_rel = np.max(np.abs(band[:n, 3] - m0[:n]) / np.abs(m0[:n])) if n else 1.0
    rel = np.abs(band[:n, 4] - sigma_direct[:n]) / sigma_direct[:n] if n else np.array([1.0])
    share = np.max(np.abs(sigma_direct - sigma_levels) / sigma_direct) if n and density else 0.0
    print(f"      {n} points; band values against the session's model to {value_rel:.1e}; the band's sigma against "
          f"the direct propagation to {rel.max():.2e} (largest); the density's share of sigma up to {share:.1%}")
    check("(b) the band file's points are the data points, its values the model's (1e-5: param.sav keeps 8 digits)",
          ok_shape and value_rel < 1e-5, (len(band), m0.size, value_rel))
    check("(b) the band = sqrt(g^T C g) of the direct propagation to 1e-4 of sigma at every point, and the density's "
          "share is larger than that tolerance (a band without it fails)", rel.max() < 1e-4 and share > 1e-2,
          (rel.max(), share))

    # (c) the density fixed at its fitted value
    sav = os.path.join(free_dir, "output", "param.sav")
    dname = density[0] if density else "TF.areal_density"
    fitted = next(float(l.split()[1]) for l in open(sav) if l.split() and l.split()[0] == dname)
    text = open(os.path.join(SRC, "fit.azr")).read()
    fixed_text = re.sub(r"(areal_density\s+)1e18 atoms/cm2 fit", lambda m_: f"{m_.group(1)}{fitted!r} atoms/cm2 fixed",
                        text, count=1)
    check("(c) the variant with the density fixed differs from fit.azr", fixed_text != text)
    fixed_dir = copy(os.path.join(tmp, "fixed"), fixed_text)
    check("(c) the fit with the density fixed ends (exit 0)", fit(fixed_dir, f"--parameter-file={sav}") == 0)
    names_f, _ = covariance(fixed_dir)
    bands_f = band_files(fixed_dir)
    band_f = np.concatenate([b for b in bands_f.values()]) if bands_f else np.zeros((0, 5))
    same_points = band_f.shape == band.shape and np.allclose(band_f[:, :4], band[:, :4], rtol=1e-4, atol=0)
    ratio = band[:, 4] / band_f[:, 4] if same_points else np.array([0.0])
    if not same_points and band_f.shape == band.shape:
        print(f"      fixed against free, columns E..value: largest relative difference "
              f"{np.max(np.abs(band_f[:, :4] - band[:, :4]) / np.maximum(np.abs(band[:, :4]), 1e-300)):.2e}")
    print(f"      fixed: covariance.dat {len(names_f)} columns; sigma(free) / sigma(fixed) from {ratio.min():.4f} to "
          f"{ratio.max():.4f}")
    check("(c) with the density fixed: no density column, the same points and values (1e-4: the density written to 8 "
          "digits, the width refitted)",
          density and density[0] not in names_f and same_points, (names_f, band_f.shape, band.shape))
    check("(c) the band with the density free is nowhere narrower (1e-3) and somewhere wider by more than 1 %",
          ratio.min() > 1 - 1e-3 and ratio.max() > 1.01, (ratio.min(), ratio.max()))

    # (d) an extrapolation after the fit: its run has no density parameter; covariance.dat is taken on its own column
    # by name, exactly as a covariance.dat that holds that column alone
    extra = ("<extrapolations>\nextrapolation x1\n    reaction       15N(p,a)12C\n    quantity       sigma\n"
             "    unit           b\n    energies       list 0.9585 0.96 0.962 0.965 0.97 MeV lab\n</extrapolations>\n")
    with open(os.path.join(free_dir, "fit_x.azr"), "w") as f:
        f.write(text.rstrip("\n") + "\n" + extra)
    env = dict(os.environ, OMP_NUM_THREADS="1")

    def extrapolate(d):
        p = subprocess.run([BIN, "extrapolate", "--no-gui", "--covariance-band", f"--parameter-file={sav}", "fit_x.azr"],
                           cwd=d, env=env, capture_output=True, text=True, timeout=900)
        files = [n for n in os.listdir(os.path.join(d, "output")) if n.endswith(".extrap.band")]
        body = open(os.path.join(d, "output", files[0])).read() if files else ""
        return p, body

    p_named, band_named = extrapolate(free_dir)
    one_dir = os.path.join(tmp, "one")
    shutil.copytree(free_dir, one_dir)
    w = names.index("width_1_1") if "width_1_1" in names else 0
    with open(os.path.join(one_dir, "output", "covariance.dat"), "w") as f:
        f.write(f"#columns width_1_1\n{C[w, w]:.17e}\n")
    p_one, band_one = extrapolate(one_dir)
    sig = [float(l.split()[4]) for l in band_named.splitlines() if l.strip()]
    print(f"      extrapolation: exit {p_named.returncode}, {len(sig)} band lines, sigma {sig[:2]}...")
    check("(d) the extrapolation takes covariance.dat (2 columns) on its one column by name: a band, no size mismatch, "
          "the same file as from a covariance.dat of that column alone",
          p_named.returncode == 0 and len(sig) == 5 and all(v > 0 for v in sig) and "mismatch" not in p_named.stdout
          and band_named == band_one, (p_named.returncode, sig, p_named.stdout[-600:]))
finally:
    os.chdir(HERE)
    shutil.rmtree(tmp, ignore_errors=True)
print(f"{len(failures)} failure(s)" + (": " + ", ".join(failures) if failures else ""))
sys.exit(1 if failures else 0)
