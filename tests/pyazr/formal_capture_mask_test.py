#!/usr/bin/env python3
"""The --ignore-externals mask: defined on the input, the same in every basis, continuous, and exact in the gradient.

Model (tests/13N): the capture segment (artemov) and the 84.3 deg elastic, plus a 1/2+ level at E_x 4.5 MeV with a
proton width of 2e4 eV and the photon width under test.  With 0 that level is masked (no external part) while its
proton width couples it to the 2.368 MeV resonance of the same block -- in the formal basis the mask is then the
projector of docs/upgrade/external_capture.md section 3, not a diagonal.

Checks:
  (a) continuity, default flags: the widths 1e-20 eV and 1e-30 eV give the same capture (relative 1e-9), and the
      width 0 (masked: no external part) differs from them; builds before 2026-11-08 masked |gamma| < 1e-8 MeV^1/2,
      so 1e-20 eV lost its external part like 0 (--no-ignore-externals is no reference here: it also un-masks the
      other levels' zero photon widths);
  (b) the formal basis gives the Brune capture with the default flags (relative 1e-9);
  (c) the chi2 gradient equals central differences in the formal and Brune bases (1e-5 per component, scaled by
      the larger of its size and 1e-3 of the largest component) -- the formal-basis chi2 gradient once missed the
      dependence of the projector on the parameters (1.7e-3 here; 1.4e-7 with the fallback, 1.1e-7 in Brune;
      external_capture.md section 7).

Margins (2026-11-08): (a) the masked and unmasked captures differ by 4.5e-4 (5e-13 on builds that masked 1e-20 eV);
(b) 6e-15 against 5.7e-4 before.

Run from anywhere:  python3 tests/pyazr/formal_capture_mask_test.py
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
KEEP = (0, 2)          # <segmentsData> lines kept: capture (artemov) and the 84.3 deg elastic
CAPTURE = 0            # its index among the kept segments


def write_model(path, photon_width_eV):
    from pyazr import AzrModel
    m = AzrModel.from_file(SRC)
    # A 1/2+ level at E_x 4.5 MeV with a proton width (fixed) and the photon width under test (fixed): with 0 it is
    # masked, and its proton width couples it to the 2.368 MeV resonance of the block.
    m.add_level(J=0.5, parity=+1, energy=4.5, level_fixed=True,
                channels=[{"pair": 1, "L": 0, "S": 0.5, "gamma": 2.0e4, "fixed": True},
                          {"pair": 2, "L": 1, "S": 0.5, "gamma": photon_width_eV, "fixed": True}])
    m.write(path)
    lines = open(path).read().split("\n")
    i0, i1 = lines.index("<segmentsData>"), lines.index("</segmentsData>")
    seg = [l for l in lines[i0 + 1:i1] if l.split()]
    for k, l in enumerate(seg):
        if k not in KEEP:
            seg[k] = "0" + l.lstrip()[1:]
    lines[i0 + 1:i1] = seg
    open(path, "w").write("\n".join(lines))


def worker(work, width, brune, ignore, grad, out):
    import numpy as np
    from pyazr import azure2
    azr = os.path.join(work, "m.azr")
    write_model(azr, float(width))
    res = {}
    with azure2(azr, cwd=work, use_brune=brune == "1", ignore_externals=ignore == "1") as a:
        x = np.array(a.params_rwa, float)
        res["capture"] = np.asarray(a.calculate_rwa(x)[CAPTURE], float).tolist()
        if grad == "1":
            c0, g = a.chi2_and_grad(x)
            g = np.asarray(g, float)
            fd = np.zeros_like(g)
            energy = {p.free_index for p in a.parameters
                      if p.kind == "energy" and not p.fixed and p.free_index is not None}
            for i in range(len(x)):
                # chi2 ~ 3e5: a level energy steps by 1e-6 MeV (resonances of tens of keV: 1e-4 relative, 240 eV,
                # leaves a truncation error of 1e-3); a width by 1e-4 relative (1e-6 loses 1e-3 of the small
                # components to cancellation in the difference)
                h = 1e-6 if i in energy else 1e-4 * max(abs(x[i]), 1e-3)
                xp, xm = x.copy(), x.copy()
                xp[i] += h
                xm[i] -= h
                fd[i] = (a.calculate_chi2_rwa(xp)[0] - a.calculate_chi2_rwa(xm)[0]) / (2 * h)
            # each component against its own size, or 1e-3 of the largest one, whichever is larger
            scale = np.maximum(np.abs(fd), 1e-3 * np.max(np.abs(fd)))
            res["grad_worst"] = float(np.max(np.abs(g - fd) / scale))
    json.dump(res, open(out, "w"))


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

    def run(width, brune, ignore, grad=False):
        work = tempfile.mkdtemp(prefix="pyazr_capture_mask_")
        try:
            shutil.copytree(os.path.join(ROOT, "tests", "13N", "data"), os.path.join(work, "data"))
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))
            out = os.path.join(work, "res.json")
            p = subprocess.run([sys.executable, os.path.abspath(__file__), "worker", work, repr(width),
                                str(int(brune)), str(int(ignore)), str(int(grad)), out],
                               capture_output=True, text=True, timeout=900)
            if p.returncode != 0 or not os.path.exists(out):
                print(p.stderr[-1500:])
                return None
            return json.load(open(out))
        finally:
            shutil.rmtree(work, ignore_errors=True)

    def rel(a, b):
        a, b = np.asarray(a), np.asarray(b)
        return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))

    tiny = run(1e-20, True, True)
    tinier = run(1e-30, True, True)
    zero_default = run(0.0, True, True, grad=True)
    zero_formal = run(0.0, False, True, grad=True)
    if None in (tiny, tinier, zero_default, zero_formal):
        check("every configuration ran", False, "see the traceback above")
    else:
        d = rel(tiny["capture"], tinier["capture"])
        check(f"(a) continuous: 1e-20 eV and 1e-30 eV give the same capture ({d:.1e})", d <= 1e-9)
        d = rel(tiny["capture"], zero_default["capture"])
        check(f"(a) 1e-20 eV keeps its external part, 0 is masked: they differ ({d:.1e})", d > 1e-5)
        d = rel(zero_formal["capture"], zero_default["capture"])
        check(f"(b) formal basis = Brune basis, default flags ({d:.1e})", d <= 1e-9)
        for label, r in (("Brune", zero_default), ("formal", zero_formal)):
            check(f"(c) {label} chi2 gradient = central differences ({r['grad_worst']:.1e})", r["grad_worst"] <= 1e-5)
    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(*sys.argv[2:8])
        sys.exit(0)
    sys.exit(main())
