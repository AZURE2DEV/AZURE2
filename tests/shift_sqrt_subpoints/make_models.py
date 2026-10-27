#!/usr/bin/env python3
"""Models for the energy-dependent-shift sub-point test (audit item 11, 2026-10-26).

All from tests/thick_target_resonance (14N(p,g) across the 0.97 keV-wide 259.7 keV resonance, 3.2e19 at/cm2
target, no external capture), the x20 target-integration line, four yield points at lab 272-320 keV:

  sqrt    : one segment, sqrtshift b = 0.02 MeV^1/2 (fixed), constant shift 0
  const   : four one-point segments, each with the CONSTANT shift b sqrt(E_k) of its own point -- the same
            calibration point by point, so the yields must be identical if sub-points move rigidly with
            their data point
  shiftfn : one segment, shiftfn "a0*sqrt(x)" with a0 = b (the same shift through the other code path)
  dense   : the sqrt model on a uniform 2000-point sub-point grid (useAdaptiveGrid=0), a grid-free check
  noshift : the sqrt model with b = 0 (to show the shift matters)
  narrow_sqrt, narrow_const : sqrt and const with the resonance's proton width 100x smaller (9.7 eV), where
            the fine grid's position matters: the old engine's grid missed the shifted resonance

Run from this directory: python3 make_models.py  (writes models/ and models/data/).
"""
import os
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "thick_target_resonance")
B = 0.02


def main():
    base = open(os.path.join(SRC, "thick_target_resonance.azr")).read()
    os.makedirs(os.path.join(HERE, "models", "data"), exist_ok=True)
    shutil.copy(os.path.join(SRC, "data", "yield_scan.dat"), os.path.join(HERE, "models", "data"))
    rows = [l.split() for l in open(os.path.join(SRC, "data", "yield_scan.dat")) if l.strip()]
    energies = [float(r[0]) for r in rows]
    for k, l in enumerate(open(os.path.join(SRC, "data", "yield_scan.dat"))):
        if l.strip():
            open(os.path.join(HERE, "models", "data", f"yield_{k + 1}.dat"), "w").write(l)

    head = base[: base.index("<segmentsData>")]
    tail = base[base.index("<segmentsTest>"):]
    # the x20 target-integration line of the source model, with the segment list replaced
    tint = [l for l in base[base.index("<targetInt>"):].splitlines() if l.startswith('1  "2"')][0]
    assert tint.endswith("20  50"), tint

    def seg(datafile, shift, extra=""):
        return (f"1  1  -1  0  6  0  180  3  360000  0  0  {shift:.10g}  0.0005  0  data/{datafile}  0  0{extra}\n")

    def model(name, segments, seglist, points=20, adaptive=True, narrow=False):
        h = head if adaptive else head.replace("useAdaptiveGrid=1", "useAdaptiveGrid=0")
        if narrow:  # the 259.7 keV resonance's proton width 974.937 -> 9.74937 eV
            assert "974.937" in h
            h = h.replace("974.937", "9.74937", 1)
        assert adaptive or "useAdaptiveGrid=0" in h, "no useAdaptiveGrid key in the source model"
        t = tint.replace('"2"', f'"{seglist}"').replace("  20  0  0.0001", f"  {points}  0  0.0001", 1)
        txt = h + "<segmentsData>\n" + "".join(segments) + "</segmentsData>\n" + tail
        txt = txt[: txt.index("<targetInt>")] + "<targetInt>\n" + t + "\n</targetInt>\n"
        open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)

    model("sqrt", [seg("yield_scan.dat", 0.0, f"  sqrtshift {B} 0 0")], "1")
    # the same with the 259.7 keV resonance 100x narrower (proton width 974.937 -> 9.74937 eV): the fine grid
    # (+-20 widths = +-0.2 keV) must sit on the shifted resonance to 0.2 keV, while the shift's slope times the
    # resonance's depth in the target is up to 0.019 x 38 keV = 0.7 keV
    model("narrow_sqrt", [seg("yield_scan.dat", 0.0, f"  sqrtshift {B} 0 0")], "1", narrow=True)
    model("narrow_const", [seg(f"yield_{k + 1}.dat", B * e ** 0.5) for k, e in enumerate(energies)], "1,2,3,4", narrow=True)
    model("const", [seg(f"yield_{k + 1}.dat", B * e ** 0.5) for k, e in enumerate(energies)], "1,2,3,4")
    model("shiftfn", [seg("yield_scan.dat", 0.0, f'  shiftfn "a0*sqrt(x)" 1 {B} 0 0')], "1")
    model("dense", [seg("yield_scan.dat", 0.0, f"  sqrtshift {B} 0 0")], "1", points=2000, adaptive=False)
    model("noshift", [seg("yield_scan.dat", 0.0)], "1")


if __name__ == "__main__":
    main()
