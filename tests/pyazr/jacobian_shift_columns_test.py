#!/usr/bin/env python3
"""residual_jacobian differences its energy-shift and density columns forward only, and a self-contained
segment's shift on that segment's rows alone.

The adjoint does not cover a shift or a fitted target density, so those columns are central differences.  Up to
2026-10-30 pyazr formed each side by re-solving the WHOLE residual Jacobian, adjoints included, and kept only its
residuals: a 17O+a rmfit polish spent ~16 of its ~18 evaluations per Jacobian (41 min each) there.  The CLI's
least-squares path had long differenced a self-contained segment's shift on its own rows; both now share that code
(FillShiftColumns), and the other columns take one forward point pass per side.

Models (tests/13N, each in a fresh process with AZURE_PROFILE=1):
  L  segments 3 and 4 (Meyer, 84.3 and 114.5 deg), a free constant shift each       -> two segment-local columns
  F  segment 3 at 0.3-0.9 MeV with target integration and ``densityfit``, segment 4 with a free shift
                                                                                    -> one full, one segment-local
Checks: (a) the "api jacobian" profile report counts the columns of each kind (a build before the change has no
such stages and fails here); (b) every differenced column equals the central difference of the residuals that
residual_jacobian itself returns at x +- h (h as the engine's), to 1e-6 of the column's largest entry; (c) the
residuals and every other column are bit-identical to a call that skips the differenced columns, and a repeated
call gives the same result (the model is left at x).

Run from anywhere:  python3 tests/pyazr/jacobian_shift_columns_test.py
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

SRC = os.path.join(ROOT, "tests", "13N", "13N.azr")
KINDS = ("shift", "shift_sqrt", "shift_fn", "density")


def make_model(work, density):
    lines = open(SRC).read().split("\n")
    i0, i1 = lines.index("<segmentsData>"), lines.index("</segmentsData>")
    seg = [l.split() for l in lines[i0 + 1:i1] if l.split()]
    for k, f in enumerate(seg):
        if k == 3:                             # segment 4: a free shift
            f[11], f[13] = "0.01", "1"
        elif k == 2:                           # segment 3: a free shift, or a fitted density at 0.3-0.9 MeV
            if density:
                f[3], f[4] = "0.3", "0.9"
            else:
                f[11], f[13] = "0.01", "1"
        else:
            f[0] = "0"
    lines[i0 + 1:i1] = [" ".join(f) for f in seg]
    if density:
        j0, j1 = lines.index("<targetInt>"), lines.index("</targetInt>")
        lines[j0 + 1:j1] = ['1 "3" 50 0 0 1 1e+18 "a0" 1 1e-20 0 0 0 "" 0 0 densityfit 1 5%']
    path = os.path.join(work, "model.azr")
    open(path, "w").write("\n".join(lines))
    return path


def worker(work, density, out):
    import numpy as np
    from pyazr import azure2
    a = azure2(make_model(work, density == "1"), cwd=work)
    x = np.array(a.params_rwa, float)
    free = [p for p in a.parameters if not p.fixed]
    cols = [i for i, p in enumerate(free) if p.kind in KINDS]
    sys.stderr.write("@@ full\n"); sys.stderr.flush()
    r, J = (np.asarray(v, float) for v in a.residual_jacobian(x))
    sys.stderr.write("@@ rest\n"); sys.stderr.flush()
    r2, J2 = (np.asarray(v, float) for v in a.residual_jacobian(x))
    try:
        rs, Js = (np.asarray(v, float) for v in a.residual_jacobian(x, skip=cols))
        kw = {"skip": cols}
    except TypeError:                          # a build before residual_jacobian(skip=): reach check (a) anyway
        rs, Js = r2, J2
        kw = {}
    fd = {}
    for c in cols:
        h = 1.0e-6 * (abs(x[c]) + 1.0)
        side = []
        for s in (+1, -1):
            xs = x.copy(); xs[c] += s * h
            side.append(np.asarray(a.residual_jacobian(xs, **kw)[0], float))
        fd[c] = ((side[0] - side[1]) / (2 * h)).tolist()
    json.dump({"kinds": [free[c].kind for c in cols], "cols": cols, "r": r.tolist(), "J": J.tolist(),
               "r2": r2.tolist(), "J2": J2.tolist(), "rs": rs.tolist(), "Js": Js.tolist(),
               "fd": {str(k): v for k, v in fd.items()}}, open(out, "w"))


def profile_calls(log, stage):
    """Calls of `stage` in the first "api jacobian" report of `log` (the full call), 0 when absent."""
    block = log.split("@@ full", 1)[-1].split("@@ rest", 1)[0]
    m = re.search(r"\[profile\]\s+" + re.escape(stage) + r"\s+[\d.]+\s+(\d+)", block)
    return int(m.group(1)) if m else 0


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

    for label, density, want in (("L", "0", {"shift cols (segment)": 2, "shift cols (full)": 0}),
                                 ("F", "1", {"shift cols (segment)": 1, "shift cols (full)": 1})):
        print(f"\nmodel {label}")
        work = tempfile.mkdtemp(prefix="pyazr_jac_shiftcols_")
        try:
            shutil.copytree(os.path.join(ROOT, "tests", "13N", "data"), os.path.join(work, "data"))
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))
            out = os.path.join(work, "res.json")
            env = dict(os.environ, AZURE_PROFILE="1")
            p = subprocess.run([sys.executable, os.path.abspath(__file__), "worker", work, density, out], env=env,
                               capture_output=True, text=True, timeout=1500)
            if p.returncode != 0:
                print(p.stderr[-2000:])
                check(f"model {label} ran", False, f"exit {p.returncode}")
                continue
            res = json.load(open(out))
            check(f"model {label}: differenced columns {res['kinds']}",
                  len(res["cols"]) == 2 and ("density" in res["kinds"]) == (density == "1"))
            # (a)
            for stage, n in want.items():
                got = profile_calls(p.stderr, stage)
                check(f"(a) profile '{stage}': {n} column(s)", got == n, f"{got}")
            # (b)
            J = np.asarray(res["J"])
            for c in res["cols"]:
                col, ref = J[:, c], np.asarray(res["fd"][str(c)])
                scale = max(float(np.max(np.abs(ref))), 1e-300)
                dev = float(np.max(np.abs(col - ref))) / scale
                check(f"(b) column {c}: central difference of the returned residuals (dev {dev:.1e} of max "
                      f"{scale:.3g})", dev <= 1e-6 and np.any(col != 0))
            # (c)
            other = [k for k in range(J.shape[1]) if k not in res["cols"]]
            Js = np.asarray(res["Js"])
            check("(c) residuals bit-identical to the call skipping the differenced columns",
                  res["r"] == res["rs"])
            check("(c) every other column bit-identical", np.array_equal(J[:, other], Js[:, other]))
            check("(c) a repeated call gives the same result", res["r"] == res["r2"] and res["J"] == res["J2"])
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
        worker(*sys.argv[2:5])
        sys.exit(0)
    sys.exit(main())
