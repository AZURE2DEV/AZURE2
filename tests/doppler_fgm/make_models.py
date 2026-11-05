#!/usr/bin/env python3
"""Models for tests/doppler_fgm, built at run time (check.sh) into a work directory.

From tests/reaction_rate/models/n12C.azr (n + 12C -> 13C, elastic and capture, external capture off), with the
s-wave 1/2+ resonance at E_cm = 0.4247 MeV narrowed from Gamma_n = 34 keV to 200 eV, so that free-gas Doppler
broadening at 300 K (Doppler width about 60 eV at E_lab = 0.46 MeV) changes the cross section by tens of percent:

  bare     sigma_t ("total 1", 1 -> 1) and sigma_gamma (1 -> 2) on the data grid, no effect
  dop      the same with "doppler T" alone on the <targetInt> line (sigma_D at the point)
  dopg     doppler T and a Gaussian resolution SIGMA_LAB (Doppler first, then the Gaussian)
  gauss    the Gaussian alone (to show the Doppler step matters under the Gaussian too)
  dopT     "transmission N" (on sigma_t) and "yield N" (on sigma_gamma) with doppler T alone
  dopgT    the same with doppler T and the Gaussian: Doppler, then the transform, then the Gaussian
  truth    test segments sigma_t and sigma_gamma on a 1 eV lab grid, no effect (mode 3)
  grad     dop (segments 1, 2), dopgT (3, 4) and dopT (5, 6) in one project, for the gradient check
  refuse_* <targetInt> lines the engine must refuse

    python3 make_models.py <n12C.azr> <workdir>
"""
import os
import re
import sys

N = 0.15            # sample thickness, atoms/barn
T_K = 300.0         # effective temperature, K
SIGMA_LAB = 1.0e-4  # Gaussian resolution, MeV lab (100 eV)
E_DATA = [round(0.4589 + 0.00005 * k, 7) for k in range(63)]   # 0.45890 .. 0.46200 MeV lab
GAMMA_N = "200"     # eV, the narrowed 1/2+ neutron width


def blocks(text):
    m = re.search(r"(.*<segmentsData>\n)(.*?)(</segmentsData>\n<segmentsTest>\n)(.*?)(</segmentsTest>\n.*?<targetInt>\n)(.*?)"
                  r"(</targetInt>.*)", text, re.S)
    if not m:
        sys.exit("make_models: unexpected .azr layout")
    return list(m.groups())


def narrow(head):
    def level(line):
        f = line.split()
        if len(f) > 12 and f[0] == "0.5" and f[1] == "1" and f[11] == "34000":
            f[11] = GAMMA_N
            return "  ".join(f) + "\n"
        return line
    out = re.sub(r"(<levels>\n)(.*?)(</levels>)",
                 lambda m: m.group(1) + "".join(level(l) for l in m.group(2).splitlines(True)) + m.group(3),
                 head, flags=re.S)
    if out.count(" " + GAMMA_N + " ") < 1:
        sys.exit("make_models: the 1/2+ level was not found")
    return out


def write(path, parts):
    with open(path, "w") as f:
        f.write("".join(parts))


def main(src, work):
    head, _, mid1, _, mid2, _, tail = blocks(open(src).read())
    head = narrow(head)
    os.makedirs(os.path.join(work, "data"), exist_ok=True)
    with open(os.path.join(work, "data", "t.dat"), "w") as f:
        for e in E_DATA:
            f.write("%.7f 0 0.5 0.05\n" % e)
    with open(os.path.join(work, "data", "g.dat"), "w") as f:
        for e in E_DATA:
            f.write("%.7f 0 1.0e-6 1.0e-7\n" % e)

    def data_line(entr, exit_, dat, kw):
        return "1 %d %d 0.4585 0.4625 0 180 0 1 0 0 0 0 0 data/%s 0 0%s\n" % (entr, exit_, dat, kw)

    # 1000 sub-points under the Gaussian: its window is wider than the Doppler one, and the 200 eV resonance must
    # still be resolved by the grid (400 gives 4e-4 at the dip of sigma_t, 1600 gives 3e-5).
    def effect(segs, gauss, doppler=True, n=None):
        n = n if n is not None else (1000 if gauss else 400)
        g = "1 %.6g" % SIGMA_LAB if gauss else "0 0"
        return '1 "%s" %d %s 0 0 "" 0 0 0 0 "" 0%s\n' % (segs, n, g, " doppler %g" % T_K if doppler else "")

    bare = data_line(1, 1, "t.dat", " total 1") + data_line(1, 2, "g.dat", "")
    thin = data_line(1, 1, "t.dat", " total 1 transmission %g" % N) + data_line(1, 2, "g.dat", " yield %g" % N)
    write(os.path.join(work, "bare.azr"), [head, bare, mid1, "", mid2, "", tail])
    write(os.path.join(work, "dop.azr"), [head, bare, mid1, "", mid2, effect("1,2", False), tail])
    write(os.path.join(work, "dopg.azr"), [head, bare, mid1, "", mid2, effect("1,2", True), tail])
    write(os.path.join(work, "gauss.azr"), [head, bare, mid1, "", mid2, effect("1,2", True, False), tail])
    write(os.path.join(work, "dopT.azr"), [head, thin, mid1, "", mid2, effect("1,2", False), tail])
    write(os.path.join(work, "dopgT.azr"), [head, thin, mid1, "", mid2, effect("1,2", True), tail])
    test = "1 1 1 0.4560 0.4650 0.000001 0 0 0 0 total 1\n1 1 2 0.4560 0.4650 0.000001 0 0 0 0\n"
    write(os.path.join(work, "truth.azr"), [head, "", mid1, test, mid2, "", tail])
    grad = bare + thin + thin
    write(os.path.join(work, "grad.azr"),
          [head, grad, mid1, "", mid2, effect("1,2,5,6", False, n=120) + effect("3,4", True, n=120), tail])
    refuse = {
        "refuse_zero": '1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler 0\n',
        "refuse_missing": '1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler\n',
    }
    for name, line in refuse.items():
        write(os.path.join(work, name + ".azr"), [head, bare, mid1, "", mid2, line, tail])


if __name__ == "__main__":
    os.makedirs(sys.argv[2], exist_ok=True)
    main(sys.argv[1], sys.argv[2])
