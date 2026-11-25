#!/usr/bin/env python3
"""Old projects as each era's setup program wrote them (the writers read from the repository's history; README.md).
usage (from the repository root): python3 tests/convert/fixtures/history/make_fixtures.py
Writes tests/convert/fixtures/history/<form>/<form>.azr and its data/ (two files of tests/13N/data).
One folder per era with p.azr and data/; the same model in each: 13N (p+12C, g+13N), two data lines, a SUM and a
RATIO composite line where the era has them, two test lines (one composite where the era has them), one target line."""
import os, re, shutil, sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", "..", "..", ".."))
src = os.path.join(root, "tests/convert/fixtures/tests/13N/13N.azr")
data = os.path.join(root, "tests/13N/data")
out = here
text = open(src).read()
levels = re.search(r"<levels>\n(.*?)</levels>", text, re.S).group(1)


def g6(x):
    """QTextStream's default for a double: %g with 6 significant digits."""
    s = "%g" % x
    return s


def f15(x):
    """qSetFieldWidth(15) on a left-aligned stream (the config writer sets AlignLeft first)."""
    s = g6(x) if isinstance(x, float) else str(x)
    return s.ljust(15)


def config(extra=""):
    rows = ["true", "output/", "checks/", "none", "none", "none", "none", "none", "none", "none", "none"]
    notes = ["#Perform A-Matrix Calculation", "#Full Path to Output Directory", "#Full Path to Checks Directory",
             "#Compond Nucleus Check", "#Boundary Condition Check", "#Data Check", "#Lo-Matrix and Penetrability Check",
             "#Legendre Polynomial Check", "#Coulomb Amplitudes Check", "#Reaction Pathway Check",
             "#Angular Distributions Check"]
    return "<config>\n" + "".join(r.ljust(100) + n + "\n" for r, n in zip(rows, notes)) + extra + "</config>\n"


def last_run():
    return '<lastRun>\n4294967295\n""\n""\n0 0\n0 ""\n-1. -1. -1.\n</lastRun>\n'


# ---- data lines: (active, entr, exit, minE, maxE, minA, maxA, isDiff, norm, vary, err, shift, shiftErr, varyShift, file)
D1 = (1, 1, 2, 0.0, 5.0, 0.0, 180.0, 0, 1.0, 0, 0.0, 0.0, 0.0, 0, "data/artemov.dat")
D2 = (1, 1, 1, 0.0, 5.0, 84.3, 84.3, 1, 1.0, 1, 5.0, 0.0, 0.005, 1, "data/meyer_84.dat")
T1 = (1, 1, 2, 0.1, 1.0, 0.01, 0.0, 0.0, 0.0, 0)          # active entr exit minE maxE eStep minA maxA aStep isDiff
T2 = (1, 1, 1, 0.1, 1.0, 0.01, 90.0, 90.0, 0.0, 1)


def data_line(d, era, adv=None, upos=None):
    head = "".join(f15(v) for v in d[:8])
    if era == "E0":  # 2020: norm vary err file (the rest of the line)
        return head + f15(d[8]) + f15(d[9]) + f15(d[10]) + d[14]
    body = head + "".join(f15(v) for v in d[8:14])
    s = body + d[14]
    if era == "E1":  # 2025-08-01: + energy shift, its error and vary flag before the file
        return s
    # 2025-08-30 on: the advanced block (or " 0") after the file; 2026-03-09 on: the UPOS field after it
    s += adv if adv is not None else " 0"
    if era in ("E4", "E5"):
        s += upos if upos is not None else " 0"
    return s


def test_line(t, era, adv=None):
    s = "".join(f15(v) for v in t[:9]) + str(t[9])
    if era in ("E0", "E1"):
        return s
    return s + (adv if adv is not None else " 0")


def comps(era, op, pairs, angles=None, scalings=None):
    """The component block as the era's writer wrote it."""
    if era == "E2":  # 2025-08-30: entrance exit
        body = "".join(" %d %d" % p for p in pairs)
        return " 1 %d %d%s" % (op, len(pairs), body)
    if era in ("E3", "E4"):  # 2025-12-22: entrance exit angle (-999 = the segment's own)
        angles = angles or [-999] * len(pairs)
        body = "".join(" %d %d %s" % (p[0], p[1], g6(float(a))) for p, a in zip(pairs, angles))
        return " 1 %d %d%s" % (op, len(pairs), body)
    if era == "E5":  # 2026-05-21: -1 marker, entrance exit angle scaling when any scaling
        angles = angles or [-999] * len(pairs)
        scalings = scalings or [-999] * len(pairs)
        body = "".join(" %d %d %s %s" % (p[0], p[1], g6(float(a)), g6(float(s))) for p, a, s in zip(pairs, angles, scalings))
        return " 1 %d -1 %d%s" % (op, len(pairs), body)


def target_line(era):
    head = f15(1) + f15('"1,2"') + f15(10) + f15(0) + f15(0.0) + f15(1) + f15(0.0)
    eq = ' "((a0*(x*1e3)^a1)/1e+21)" ' + "2 " + "2.40289 0.00491497 "
    q = "              0              0 "
    if era in ("E0", "E1-early"):
        return head + eq + q
    conv_noeq = "              0              0 "
    conv = "              0               \"\" 0 "
    if era == "T1":  # 2025-01-26 .. 01-30: the convolution coefficients without an equation token
        return head + eq + q + conv_noeq
    if era == "T1c":  # the same with an energy-dependent sigma: c0 + c1 E (lab MeV)
        return head.replace(f15(0) + f15(0.0) + f15(1), f15(1) + f15(0.0) + f15(1), 1) + eq + q + "              1              2 0.001 0.0005 "
    s = head + eq + q + conv
    if era in ("E1", "E2", "E3"):
        return s
    s += "              0 0.04"  # 2026-01-22: straggling flag and coefficient
    if era == "E4":
        return s + " 5 50"  # 2026-03-28: adaptive grid
    return s + " 5 50"


def project(era):
    lv = levels
    ext = ""
    if era in ("E-1", "E-1n"):  # before 2020: 30 columns (no ecMultMask), external capture in its own section
        rows = []
        for l in lv.splitlines():
            t = l.split()
            if len(t) >= 31:
                t = t[:30]
            rows.append(" ".join(t))
        lv = "\n".join(rows) + "\n"
        ext = "<externalCapture>\n1 2 0.5 0.5 5\n</externalCapture>\n" if era == "E-1" else ""
    e = "E0" if era in ("E-1", "E-1n", "T1", "T1c") else era
    dl = [data_line(D1, e), data_line(D2, e)]
    tl = [test_line(T1, e), test_line(T2, e)]
    if e in ("E2", "E3", "E4", "E5"):
        dl.append(data_line(D2, e, adv=comps(e, 0, [(1, 1), (1, 2)])))
        dl.append(data_line(D2, e, adv=comps(e, 1, [(1, 2)])))
        tl.append(test_line(T1, e, adv=comps(e, 0, [(1, 2), (1, 2)])))
    if e in ("E3", "E4"):
        dl.append(data_line(D2, e, adv=comps(e, 0, [(1, 1), (1, 2)], angles=[30, 60])))
    if e == "E5":
        dl.append(data_line(D2, e, adv=comps(e, 0, [(1, 1), (1, 2)], angles=[-999, 60], scalings=[0.5, 2])))
        dl.append(data_line(D2, e, adv=comps(e, 0, [(1, 1)])))  # a line without scaling in a file that has some
    if e in ("E4", "E5"):
        dl.append(data_line(D1, e, upos=" 1 2 1.5 0.1"))
    tg = target_line(era if era in ("T1", "T1c") else e)
    s = config() + "<levels>\n" + lv + "</levels>\n"
    s += "<segmentsData>\n" + "\n".join(dl) + "\n</segmentsData>\n"
    s += "<segmentsTest>\n" + "\n".join(tl) + "\n</segmentsTest>\n"
    s += "<targetInt>\n" + tg + "\n</targetInt>\n"
    if e not in ("E0",):
        s += "<parameterSettings>\n# Fitting parameter settings (only non-fixed parameters shown)\n" \
             "# Format: name value lower_limit upper_limit error fit_error nuisance category minuit_index\n" \
             "</parameterSettings>\n"
    s += last_run()
    if e in ("E2", "E3", "E4", "E5"):
        s += "<mcmc>\nnwalkers 32\nnsteps 1000\nchainspread 0.01\nnthreads 4\nusereducedwidths 0\n<parameters>\n</parameters>\n</mcmc>\n"
    return s + ext


NAMES = {"E-1": "pre2020_external_capture", "E-1n": "pre2020_30_columns", "E0": "gui_2020", "T1": "target_2025_01_no_equation",
         "T1c": "target_2025_01_coefficients", "E1": "energy_shift_2025_08", "E2": "components_2025_08",
         "E3": "components_angle_2025_12", "E4": "secondary_decay_2026_03", "E5": "components_scaling_2026_05"}
for era, name in NAMES.items():
    d = os.path.join(out, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "data"))
    for f in ("artemov.dat", "meyer_84.dat"):
        shutil.copy(os.path.join(data, f), os.path.join(d, "data", f))
    open(os.path.join(d, name + ".azr"), "w").write(project(era))
