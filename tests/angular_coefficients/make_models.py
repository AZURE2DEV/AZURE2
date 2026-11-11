#!/usr/bin/env python3
"""Models for the angular-coefficient data segments (isDiff 9), 2026-10-26.

13C + alpha in 17O (the two-level model of tests/sum_lab_frame: 5/2+ and 3/2+, exits n + 16O(0+) key 3 and
n + 16O(3-) key 5, no external capture), reaction alpha + 13C -> n + 16O(3-), whose angular distribution has
Legendre orders up to 4 from the 5/2+ / 3/2+ interference.

  extrap    : test segments at lab 5.9-6.4 MeV, step 0.05: T1 the Legendre coefficients (isDiff 3, order 6),
              T2 the c.m. differential cross section 0-180 deg step 5 (isDiff 5), T3 the angle-integrated
              cross section (isDiff 0); verify.py projects T2 on P_k and compares with T1 and T3
  data      : two angular-coefficient DATA segments on the four energies 6.0-6.3: normalised a_k (k = 1-4)
              and absolute A_k ("absolute 1"), data written by verify.py from the extrapolation
  shift     : the normalised segment with a fixed 0.05 MeV energy shift (the model moves to E + 0.05)
  conv      : the normalised and the absolute segments under a 30 keV Gaussian beam-energy spread
              (<targetInt> convolution, 60 sub-points), and the same spread on the T2 test segment, so the
              folded coefficients can be compared with the projection of the folded distribution
  fit       : the normalised, absolute and convolved segments with the seven widths free and the data a few
              per cent off the model, for the AZURE_GRAD_CHECK gradient test
  bad_norm  : the normalised segment with a free normalization (refused)
  bad_elastic: a coefficient segment on the entrance pair (refused)
  bad_test  : isDiff 9 in <segmentsTest> (refused)

Run from this directory: python3 make_models.py  (writes models/; data/ is written by verify.py at check time).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "sum_lab_frame", "models", "extrap.azr")


def main():
    src = open(SRC).read()
    head = src[: src.index("<segmentsData>")]           # config + levels (widths fixed)
    # the level lines are "J pi E fixE aa key 2s 2l id active fixWidth width ...": fixWidth 0 for the free model
    lines = []
    for l in head.splitlines():
        t = l.split()
        if len(t) > 12 and t[0] in ("2.5", "1.5") and t[1] == "1":
            t[10] = "0"
            lines.append("  ".join(t))
        else:
            lines.append(l)
    head_free = "\n".join(lines) + "\n"
    os.makedirs(os.path.join(HERE, "models"), exist_ok=True)

    def write(name, data_lines, test_lines, tint_lines=(), free=False):
        txt = ((head_free if free else head) + "<segmentsData>\n" + "".join(s + "\n" for s in data_lines)
               + "</segmentsData>\n<segmentsTest>\n" + "".join(s + "\n" for s in test_lines)
               + "</segmentsTest>\n<targetInt>\n" + "".join(s + "\n" for s in tint_lines) + "</targetInt>\n")
        open(os.path.join(HERE, "models", name + ".azr"), "w").write(txt)

    grid = "5.9  6.4  0.05"
    tests = [f"1  2  5  {grid}  0  0  0  3  6",       # T1 coefficients, order 6
             f"1  2  5  {grid}  0  180  5  5",        # T2 c.m. differential
             f"1  2  5  {grid}  0  0  0  0"]          # T3 angle-integrated
    write("extrap", [], tests)

    def seg(datafile, shift=0.0, vary_norm=0, extra=""):
        return f"1  2  5  5.9  6.4  0  180  9  1.0  {vary_norm}  0  {shift:g}  0  0  data/{datafile}  0  0{extra}"

    write("data", [seg("coef_norm.dat"), seg("coef_abs.dat", extra="  absolute 1")], [])
    write("shift", [seg("coef_norm.dat", shift=0.05)], [])
    # a Gaussian beam-energy spread of 30 keV (lab), no target thickness: targetInt
    #   active "segs" points isConvolution sigma isTargetIntegration density stoppingEq nparams ...
    conv = '1  "{segs}"  60  1  0.03  0  0  ""  0  0  0  0  ""  0  0  0.04  5  50'
    write("conv", [seg("coef_norm.dat"), seg("coef_abs.dat", extra="  absolute 1")], [], [conv.format(segs="1,2")])
    write("conv_diff", [], [f"1  2  5  {grid}  0  180  5  5"], [conv.format(segs="1")])
    write("fit", [seg("coef_norm_off.dat"), seg("coef_abs_off.dat", extra="  absolute 1"), seg("coef_norm_off.dat")],
          [], [conv.format(segs="3")], free=True)
    write("bad_norm", [seg("coef_norm.dat", vary_norm=1)], [])
    write("bad_elastic", ["1  2  2  5.9  6.4  0  180  9  1.0  0  0  0  0  0  data/coef_norm.dat  0  0"], [])
    write("bad_test", [], [f"1  2  5  {grid}  0  0  0  9"])


if __name__ == "__main__":
    main()
