#!/usr/bin/env python3
"""A level whose widths are all zero changes nothing, also at a data point exactly at its energy.

Such a level decouples: its row and column of the level matrix are zero apart from the diagonal E_lambda - E, so at
E = E_lambda the level matrix is singular although the level contributes nothing.  Before 2026-10-31 the A-matrix's
level route then dropped the whole J-group at that point (the factorization failed and every bilinear form of the
group came back zero), and the R-matrix route divided 0 by 0.  Seen in 12C+12C (a zero-width 4+ placed at
E_cm 3.55 MeV, a Tan2024 data energy, in a 16-level group: chi2 4353 -> 4505).

Model (tests/13N): the 84.3 deg elastic and the capture segment (A-matrix; the R-matrix variants the elastic only), plus a 1/2+ level with every width 0 at
E_x = E_cm + S_p of one elastic point (formed as the engine forms its energy, so the hit is exact).  Each variant
runs in its own process: A-matrix Brune and standard, each with the default route and with the level route forced
(AZURE_LEVEL_SPACE=1); R-matrix Brune and standard.  Every calculated point must equal the model without the level
(relative 1e-10), and the residual Jacobian must be finite.

Run from anywhere:  python3 tests/pyazr/zero_width_level_test.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

SRC = os.path.join(ROOT, "tests", "13N", "13N.azr")
# <segmentsData> lines kept: capture (0-5 MeV) and the 84.3 deg elastic, whose middle point is hit.  The R-matrix
# variants keep the elastic only: their capture is not the A-matrix's (no external capture, and NaN above the
# 1/2+ resonance in the standard formalism, with or without the extra level -- a separate matter).
KEEP = {False: (0, 2), True: (2,)}


def write_model(path, rmatrix, level_energy=None):
    from pyazr import AzrModel
    m = AzrModel.from_file(SRC)
    if level_energy is not None:
        m.add_level(J=0.5, parity=+1, energy=level_energy, channels=[], level_fixed=True)
    m.write(path)
    lines = open(path).read().split("\n")
    i0, i1 = lines.index("<segmentsData>"), lines.index("</segmentsData>")
    seg = [l for l in lines[i0 + 1:i1] if l.split()]
    for k, l in enumerate(seg):
        if k not in KEEP[rmatrix]:
            seg[k] = "0" + l.lstrip()[1:]
    lines[i0 + 1:i1] = seg
    c = lines.index("<config>") + 1
    rest = (lines[c].split(None, 1) + [""])[1]
    lines[c] = ("false" if rmatrix else "true") + "  " + rest                  # Perform A-Matrix Calculation
    open(path, "w").write("\n".join(lines))


def worker(work, rmatrix, brune, out):
    import numpy as np
    from pyazr import azure2
    rmatrix, brune = rmatrix == "1", brune == "1"
    base = os.path.join(work, f"base_{int(rmatrix)}.azr")
    write_model(base, rmatrix)
    with azure2(base, cwd=work, use_brune=brune) as a:
        x = np.array(a.params_rwa, float)
        vb = [np.asarray(v, float).tolist() for v in a.calculate_rwa(x)]
        en = [np.asarray(e, float) for e in a.energies]
        pair = [p for p in a.pairs if p.number == 1][0]
        el = len(en) - 1                     # the elastic segment
        k = len(en[el]) // 2
        e_hit = float(en[el][k] + pair.sep_energy + pair.excitation)
    zw = os.path.join(work, f"zw_{int(rmatrix)}.azr")
    write_model(zw, rmatrix, e_hit)
    with azure2(zw, cwd=work, use_brune=brune) as a:
        x = np.array(a.params_rwa, float)
        vz = [np.asarray(v, float).tolist() for v in a.calculate_rwa(x)]
        r, J = (np.asarray(v, float) for v in a.residual_jacobian(x))
        en2 = [np.asarray(e, float) for e in a.energies]
        hit = float(en2[el][k] + pair.sep_energy + pair.excitation) == e_hit
    json.dump({"vb": vb, "vz": vz, "hit": hit, "k": k, "e_hit": e_hit,
               "jac_finite": bool(np.all(np.isfinite(J)) and np.all(np.isfinite(r)))}, open(out, "w"))


def main():
    try:
        import numpy as np
        from pyazr import azure2  # noqa: F401
    except Exception as err:  # same skip as the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 0
    failures = []

    def check(name, ok, detail=""):
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
        if not ok:
            failures.append(name)

    variants = [("A-matrix Brune", "0", "1", {}),
                ("A-matrix Brune, level route", "0", "1", {"AZURE_LEVEL_SPACE": "1"}),
                ("A-matrix standard", "0", "0", {}),
                ("A-matrix standard, level route", "0", "0", {"AZURE_LEVEL_SPACE": "1"}),
                ("R-matrix Brune", "1", "1", {}),
                ("R-matrix standard", "1", "0", {})]
    for label, rmatrix, brune, extra in variants:
        work = tempfile.mkdtemp(prefix="pyazr_zero_width_")
        try:
            shutil.copytree(os.path.join(ROOT, "tests", "13N", "data"), os.path.join(work, "data"))
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))
            out = os.path.join(work, "res.json")
            env = dict(os.environ, **extra)
            p = subprocess.run([sys.executable, os.path.abspath(__file__), "worker", work, rmatrix, brune, out],
                               env=env, capture_output=True, text=True, timeout=900)
            if p.returncode != 0 or not os.path.exists(out):
                print(p.stderr[-1500:])
                check(f"{label}: ran", False, f"exit {p.returncode}")
                continue
            res = json.load(open(out))
            check(f"{label}: the level sits exactly on elastic point {res['k']} (E_x {res['e_hit']!r})", res["hit"])
            worst, where = 0.0, ""
            for s, (a, b) in enumerate(zip(res["vz"], res["vb"])):
                a, b = np.asarray(a), np.asarray(b)
                d = np.abs(a - b) / np.maximum(np.abs(b), 1e-300)
                d[~np.isfinite(a)] = np.inf
                if d.size and d.max() > worst:
                    worst, where = float(d.max()), f"segment {s} point {int(d.argmax())}: {b[d.argmax()]:.6e} -> " \
                                                   f"{a[d.argmax()]:.6e}"
            check(f"{label}: every point as without the level (worst {worst:.1e})", worst <= 1e-10, where)
            check(f"{label}: residual Jacobian finite", res["jac_finite"])
        finally:
            shutil.rmtree(work, ignore_errors=True)
    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(*sys.argv[2:6])
        sys.exit(0)
    sys.exit(main())
