#!/usr/bin/env python3
"""Models for tests/transmission_yield, built at run time (check.sh) into a work directory.

From tests/reaction_rate/models/n12C.azr (n + 12C -> 13C, elastic and capture, external capture off; the s-wave 1/2+
resonance at E_cm = 0.425 MeV with Gamma_n = 34 keV):

  bare     segments sigma_t ("total 1", 1 -> 1) and sigma_gamma (1 -> 2) on the data grid
  thin     the same lines with "transmission N" and "yield N": T = exp(-N sigma_t), Y0 = (sg/st)(1 - exp(-N st))
  conv     thin with a Gaussian beam spread of SIGMA_LAB on both segments (the transform before the kernel)
  bconv    bare with the same spread (broadening sigma first: what the transform must NOT do)
  truth    test segments sigma_t and sigma_gamma on a 0.1 keV lab grid, no spread (mode 3)
  xbare, xthin   test segments (mode 3): sigma_t, sigma_gamma and their transmission and yield
  grad     thin + conv in one project, every width free, for the gradient check
  refuse_* lines the engine must refuse

From tests/natural_neon/natural_neon.azr (three compound nuclei, natural sums of totals):

  nat      only segment 1 (Vaughn, natural sum, no target effect) active
  natT     the same line with "transmission NNAT": T = exp(-NNAT sigma_nat)
  ybare, yelast, yreac   n + 20Ne with 1 keV reaction widths: sigma_t ("total 1") and sigma(n,alpha); the elastic;
                         a yield to alpha + 17O (its sigma_t must be the total, which differs from the elastic here)
  refuse_natsub   segment 7 (Junghans, natural sum under an energy-dependent spread) with "transmission"
  refuse_natyield a natural reaction line (n+20Ne -> alpha+17O) with "yield"

    python3 make_models.py <n12C.azr> <natural_neon.azr> <workdir>
"""
import os
import re
import sys

N = 0.15            # sample thickness, atoms/barn (n sigma_t ~ 1 at the 1/2+ peak)
NNAT = 0.05
SIGMA_LAB = 0.015   # Gaussian beam spread, MeV lab
E_DATA = [round(0.300 + 0.005 * k, 6) for k in range(71)]   # 0.300 .. 0.650 MeV lab


def blocks(text):
    """(head, segmentsData body, between, segmentsTest body, between, targetInt body, tail)"""
    m = re.search(r"(.*<segmentsData>\n)(.*?)(</segmentsData>\n<segmentsTest>\n)(.*?)(</segmentsTest>\n.*?<targetInt>\n)(.*?)"
                  r"(</targetInt>.*)", text, re.S)
    if not m:
        sys.exit("make_models: unexpected .azr layout")
    return list(m.groups())


def write(path, parts):
    with open(path, "w") as f:
        f.write("".join(parts))


def n12c_models(src, work):
    text = open(src).read()
    # free every width that is fixed only by the fix flag (column 11) of a level line, for the gradient check
    head, _, mid1, _, mid2, _, tail = blocks(text)
    os.makedirs(os.path.join(work, "data"), exist_ok=True)
    with open(os.path.join(work, "data", "d.dat"), "w") as f:
        for e in E_DATA:
            f.write("%.6f 0 0.5 0.05\n" % e)
    with open(os.path.join(work, "data", "g.dat"), "w") as f:
        for e in E_DATA:
            f.write("%.6f 0 1.0e-7 1.0e-8\n" % e)

    def data_line(entr, exit_, dat, kw):
        return "1 %d %d 0.29 0.66 0 180 0 1 0 0 0 0 0 data/%s 0 0%s\n" % (entr, exit_, dat, kw)

    conv = '1 "1,2" 200 1 %.6f 0 0 "" 0 0 0 0 "" 0\n' % SIGMA_LAB
    off = '0 "1,2" 200 1 %.6f 0 0 "" 0 0 0 0 "" 0\n' % SIGMA_LAB
    bare = data_line(1, 1, "d.dat", " total 1") + data_line(1, 2, "g.dat", "")
    thin = data_line(1, 1, "d.dat", " total 1 transmission %g" % N) + data_line(1, 2, "g.dat", " yield %g" % N)
    write(os.path.join(work, "bare.azr"), [head, bare, mid1, "", mid2, off, tail])
    write(os.path.join(work, "thin.azr"), [head, thin, mid1, "", mid2, off, tail])
    write(os.path.join(work, "conv.azr"), [head, thin, mid1, "", mid2, conv, tail])
    write(os.path.join(work, "bconv.azr"), [head, bare, mid1, "", mid2, conv, tail])
    test = "1 1 1 0.20 0.75 0.0001 0 0 0 0 total 1\n1 1 2 0.20 0.75 0.0001 0 0 0 0\n"
    write(os.path.join(work, "truth.azr"), [head, "", mid1, test, mid2, off, tail])
    # test segments (mode 3): bare and transformed curves
    xb = "1 1 1 0.30 0.65 0.001 0 0 0 0 total 1\n1 1 2 0.30 0.65 0.001 0 0 0 0\n"
    xt = "1 1 1 0.30 0.65 0.001 0 0 0 0 total 1 transmission %g\n1 1 2 0.30 0.65 0.001 0 0 0 0 yield %g\n" % (N, N)
    write(os.path.join(work, "xbare.azr"), [head, "", mid1, xb, mid2, off, tail])
    write(os.path.join(work, "xthin.azr"), [head, "", mid1, xt, mid2, off, tail])
    # gradient: thin (points without sub-points) and conv (with) in one project, segments 3,4 under the spread
    grad_conv = '1 "3,4" 60 1 %.6f 0 0 "" 0 0 0 0 "" 0\n' % SIGMA_LAB
    g = re.sub(r"(<levels>\n)(.*?)(</levels>)",
               lambda m: m.group(1) + "".join(free_width(l) for l in m.group(2).splitlines(True)) + m.group(3),
               head, flags=re.S)
    write(os.path.join(work, "grad.azr"), [g, thin + thin, mid1, "", mid2, grad_conv, tail])
    # refusals (one data line each)
    refuse = {
        "refuse_both": data_line(1, 2, "g.dat", " yield %g transmission %g" % (N, N)),
        "refuse_yield_total": data_line(1, 1, "d.dat", " total 1 yield %g" % N),
        "refuse_diff": data_line(1, 1, "d.dat", " transmission %g" % N).replace(" 0 180 0 1 ", " 0 180 1 1 ", 1),
    }
    for name, line in refuse.items():
        write(os.path.join(work, name + ".azr"), [head, line, mid1, "", mid2, off, tail])
    # charged entrance: the same line on the p + 12C model next to it
    p12c = os.path.join(os.path.dirname(src), "p12C.azr")
    ph, _, pm1, _, pm2, _, pt = blocks(open(p12c).read())
    write(os.path.join(work, "refuse_charged.azr"),
          [ph, data_line(1, 1, "d.dat", " transmission %g" % N), pm1, "", pm2, off, pt])


def free_width(line):
    f = line.split()
    # a level line: column 10 (1-based 11) is the width's fix flag; free the open channels with a nonzero width
    if len(f) > 30 and f[10] == "1" and float(f[11]) != 0.0 and float(f[2]) < 10.0:
        f[10] = "0"
        return "  ".join(f) + "\n"
    return line


def neon_models(src, work):
    text = open(src).read()
    head, data, mid1, _, mid2, target, tail = blocks(text)
    lines = data.splitlines(True)

    def only(k, line=None):
        out = []
        for i, l in enumerate(lines):
            if i + 1 == k:
                out.append(line if line is not None else l)
            else:
                out.append("0" + l[1:] if l.startswith("1 ") else l)
        return "".join(out)

    seg1 = lines[0].rstrip("\n")
    write(os.path.join(work, "nat.azr"), [head, only(1), mid1, "", mid2, target, tail])
    write(os.path.join(work, "natT.azr"), [head, only(1, seg1 + " transmission %g\n" % NNAT), mid1, "", mid2, target, tail])
    # the yield's sigma_t is the TOTAL, not the elastic: n + 20Ne with a free 1 keV width on every open reaction
    # channel (as tests/total_optical), a yield to alpha + 17O (pair 1) next to sigma_t, sigma(n,alpha) and elastic
    k = [0]

    def react(l):
        f = l.split()
        if len(f) > 30 and f[5] in ("1", "3", "5", "7") and float(f[2]) > float(f[22]) + float(f[25]) + 0.05:
            f[10] = "0"
            f[11] = "-1000" if k[0] % 2 else "1000"
            k[0] += 1
            return " ".join(f) + "\n"
        return l
    rhead = re.sub(r"(<levels>\n)(.*?)(</levels>)",
                   lambda m: m.group(1) + "".join(react(l) for l in m.group(2).splitlines(True)) + m.group(3),
                   head, flags=re.S)
    off = "".join(("0" + l[1:]) if l.startswith("1 ") else l for l in target.splitlines(True))
    line = "1 2 %d 0 2.33 0 180 0 1 0 0 0 0 0 data/Vaughn_ntot.dat 0 0%s\n"
    write(os.path.join(work, "ybare.azr"), [rhead, line % (2, " total 1") + line % (1, ""), mid1, "", mid2, off, tail])
    write(os.path.join(work, "yelast.azr"), [rhead, line % (2, ""), mid1, "", mid2, off, tail])
    write(os.path.join(work, "yreac.azr"), [rhead, line % (1, " yield %g" % NNAT), mid1, "", mid2, off, tail])
    seg7 = lines[6].rstrip("\n")
    write(os.path.join(work, "refuse_natsub.azr"),
          [head, only(7, seg7 + " transmission %g\n" % NNAT), mid1, "", mid2, target, tail])
    write(os.path.join(work, "refuse_natyield.azr"),
          [head, only(1, "1 2 3 0 2.33 0 180 0 1 0 0 0 0 0 data/Vaughn_ntot.dat 0 0 natural 1 yield %g\n" % NNAT),
           mid1, "", mid2, target, tail])


if __name__ == "__main__":
    n12c, neon, work = sys.argv[1:4]
    os.makedirs(work, exist_ok=True)
    n12c_models(n12c, work)
    neon_models(neon, work)
