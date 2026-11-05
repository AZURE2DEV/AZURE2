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
  grad_fit the same with fitted temperatures: "doppler 300 1 0 900" on the line of 1, 2, 5, 6 and
           "doppler 300 1 20" on the line of 3, 4 (doppler_1, doppler_2 free with the widths)
  dop450   dop at T_TRUE = 450 K: its sigma_t and sigma_gamma become the data of the fits below
  refuse_* <targetInt> lines the engine must refuse

    python3 make_models.py <n12C.azr> <workdir>

and, once dop450 has been calculated (mode 1),

    python3 make_models.py fitdata <n12C.azr> <dop450 dir> <workdir>

writes data/t450.dat, data/g450.dat (the dop450 values, 1 % errors), data/t450w.dat, data/g450w.dat (the same
values, 20 % errors) and, with every level parameter fixed,

  fit       the 1 % data, the line "doppler 300 1 0 900": start at 300 K, no prior, grids laid out for 900 K
  fit_prior the 20 % data (T alone known to about 60 K), the line "doppler 300 1 20": a 20 K prior about 300 K
            (grids for max(300 + 60, 1.21 * 300) = 363 K), so the fit lands between 300 K and 363 K
"""
import os
import re
import sys

N = 0.15            # sample thickness, atoms/barn
T_K = 300.0         # effective temperature, K
SIGMA_LAB = 1.0e-4  # Gaussian resolution, MeV lab (100 eV)
E_DATA = [round(0.4589 + 0.00005 * k, 7) for k in range(63)]   # 0.45890 .. 0.46200 MeV lab
GAMMA_N = "200"     # eV, the narrowed 1/2+ neutron width
T_TRUE = 450.0      # K, the temperature of the fit data


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


def fix_all(head):
    """Every level energy and width fixed (fields 4 and 11 of a <levels> line set to 1)."""
    def level(line):
        f = line.split()
        if len(f) > 12:
            f[3] = "1"
            f[10] = "1"
            return "  ".join(f) + "\n"
        return line
    return re.sub(r"(<levels>\n)(.*?)(</levels>)",
                  lambda m: m.group(1) + "".join(level(l) for l in m.group(2).splitlines(True)) + m.group(3),
                  head, flags=re.S)


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
    def effect(segs, gauss, doppler=True, n=None, fit=""):
        n = n if n is not None else (1000 if gauss else 400)
        g = "1 %.6g" % SIGMA_LAB if gauss else "0 0"
        return '1 "%s" %d %s 0 0 "" 0 0 0 0 "" 0%s\n' % (segs, n, g, (" doppler %g" % T_K + fit) if doppler else "")

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
    write(os.path.join(work, "grad_fit.azr"),
          [head, grad, mid1, "", mid2,
           effect("1,2,5,6", False, n=120, fit=" 1 0 900") + effect("3,4", True, n=120, fit=" 1 20"), tail])
    dop450 = effect("1,2", False).replace(" doppler %g" % T_K, " doppler %g" % T_TRUE)
    write(os.path.join(work, "dop450.azr"), [head, bare, mid1, "", mid2, dop450, tail])
    refuse = {
        "refuse_zero": '1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler 0\n',
        "refuse_missing": '1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler\n',
        # the fit form: vary not 0/1, sigma missing or negative, Tmax not above T, T <= 0 with a fit form
        "refuse_vary": '1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler 300 2 0\n',
        "refuse_nosigma": '1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler 300 1\n',
        "refuse_sigma": '1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler 300 1 -5\n',
        "refuse_tmax": '1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler 300 1 0 200\n',
        "refuse_negT": '1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler -300 1 0 900\n',
    }
    for name, line in refuse.items():
        write(os.path.join(work, name + ".azr"), [head, bare, mid1, "", mid2, line, tail])


def fitdata(src, dop450, work):
    """The fit models on the dop450 values (sigma_t, sigma_gamma at T_TRUE), 1 % errors, every level fixed."""
    head, _, mid1, _, mid2, _, tail = blocks(open(src).read())
    head = fix_all(narrow(head))
    for name, exit_, dat in (("t450.dat", 1, "1_R=1"), ("g450.dat", 2, "1_R=2")):
        rows = []
        for line in open(os.path.join(dop450, "output", "AZUREOut_aa=%s.out" % dat)):
            f = line.split()
            if len(f) >= 4:
                rows.append(float(f[3]))
        if len(rows) != len(E_DATA):
            sys.exit("make_models fitdata: %d values in %s, want %d" % (len(rows), dat, len(E_DATA)))
        for suffix, err in (("", 0.01), ("w", 0.20)):
            with open(os.path.join(work, "data", name.replace(".dat", suffix + ".dat")), "w") as f:
                for e, v in zip(E_DATA, rows):
                    f.write("%.7f 0 %.10e %.10e\n" % (e, v, err * v))

    def data_line(exit_, dat, kw):
        return "1 1 %d 0.4585 0.4625 0 180 0 1 0 0 0 0 0 data/%s 0 0%s\n" % (exit_, dat, kw)

    segs = data_line(1, "t450.dat", " total 1") + data_line(2, "g450.dat", "")
    wide = data_line(1, "t450w.dat", " total 1") + data_line(2, "g450w.dat", "")
    line = '1 "1,2" 400 0 0 0 0 "" 0 0 0 0 "" 0 doppler %g%s\n'
    write(os.path.join(work, "fit.azr"), [head, segs, mid1, "", mid2, line % (T_K, " 1 0 900"), tail])
    write(os.path.join(work, "fit_prior.azr"), [head, wide, mid1, "", mid2, line % (T_K, " 1 20"), tail])


if __name__ == "__main__":
    if sys.argv[1] == "fitdata":
        fitdata(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        os.makedirs(sys.argv[2], exist_ok=True)
        main(sys.argv[1], sys.argv[2])
