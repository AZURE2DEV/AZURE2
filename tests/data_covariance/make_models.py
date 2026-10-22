#!/usr/bin/env python3
"""Write the models and covariance files of tests/data_covariance (run once; outputs are committed).

    python3 make_models.py

Base: tests/angle_mapping/angle_mapping.azr (p+12C and p+13C, thirteen segments on one energy
grid).  Changes, in every model:
  * segment 3 (data/p12c_120_150.dat, two angles in one file) is split into two segments by
    angle window (115-125 and 145-155 lab degrees): one data file, two segments;
  * segment 2 (90 deg) is cut at 2.0 MeV, so rows of the covariance file drop out.
The covariance (models with "covariance <file>") covers segments 1 (free norm), 2, 3a, 3b,
4 (A_y), 6 and 7 (shift group with a free shift) and 8 (free sqrt shift); its header also lists
data/p12c_165.dat, which no covariance segment uses (its rows drop out as well).

  plain.azr   no covariance
  diag.azr    covariance cov_diag.txt  = diagonal(e_i^2): chi2 identical to plain.azr
  dense.azr   covariance cov_dense.txt = 0.5 diag(e^2) + one block per beam energy (all
              rows at that energy, non-contiguous: 0.4 e_i e_j off the diagonal) + per file
              an absolute 2 % vector + one smooth global vector
  rel.azr     covariance cov_rel.txt   = cov_dense.txt with the per-file 2 % vectors as
              relvector 0.02 (model-scaled) instead
  fit.azr     dense.azr with every level energy and width fixed: only the norms, shifts and
              abundances are free, a problem with one minimum, for the fitter comparison
              (with the levels free, MIGRAD, LM and GSL reach different minima of this
              synthetic model, depending on the platform)
"""
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "..", "angle_mapping")

COV_FILES = ["data/p12c_60.dat", "data/p12c_90.dat", "data/p12c_120_150.dat", "data/p12c_ay100.dat",
             "data/p12c_75.dat", "data/p12c_135.dat", "data/p12c_105.dat", "data/p12c_165.dat"]
COV_SEGMENT_FILES = set(COV_FILES[:-1])  # every listed file but p12c_165 carries the keyword


def records(path):
    rows = []
    with open(os.path.join(BASE, path)) as f:
        for line in f:
            t = line.split()
            if not t or t[0].startswith("#"):
                continue
            rows.append([float(x) for x in t[:4]])
    return rows


def segment_lines(azr):
    text = open(azr).read()
    a = text.index("<segmentsData>") + len("<segmentsData>\n")
    b = text.index("</segmentsData>")
    return text, a, b, [l for l in text[a:b].splitlines() if l.strip()]


def model(cov):
    text, a, b, lines = segment_lines(os.path.join(BASE, "angle_mapping.azr"))
    out = []
    for k, line in enumerate(lines):
        t = line.split()
        if k == 1:
            t[4] = "2.0"  # segment 2: energy cut
        if k == 2:
            for lo, hi in (("115", "125"), ("145", "155")):
                u = list(t)
                u[5], u[6] = lo, hi
                out.append(u)
            continue
        out.append(t)
    lines = []
    for t in out:
        path = t[14]
        if cov and path in COV_SEGMENT_FILES:
            t = t + ["covariance", cov]
        lines.append(" ".join(t))
    return text[:a] + "\n".join(lines) + "\n" + text[b:]


def header(f):
    for path in COV_FILES:
        f.write("file %s %d\n" % (path, len(records(path))))


def numbers(f, xs, per=6):
    for i in range(0, len(xs), per):
        f.write(" ".join("%.15e" % x for x in xs[i:i + per]) + "\n")


def main():
    rows, energies = [], []
    for path in COV_FILES:
        for r in records(path):
            rows.append((path, r))
    N = len(rows)
    err = [r[3] for _, r in rows]
    val = [r[2] for _, r in rows]
    ener = [round(r[0], 9) for _, r in rows]

    for name in ("plain", "diag", "dense", "rel"):
        cov = None if name == "plain" else "cov_%s.txt" % name
        with open(os.path.join(HERE, name + ".azr"), "w") as f:
            f.write(model(cov))

    text = open(os.path.join(HERE, "dense.azr")).read()
    a, b = text.index("<levels>\n"), text.index("</levels>")
    lines = []
    for line in text[a + 9:b].splitlines():
        t = line.split()
        if t:
            t[3], t[10] = "1", "1"  # energy fixed, width fixed
            line = " ".join(t)
        lines.append(line)
    with open(os.path.join(HERE, "fit.azr"), "w") as f:
        f.write(text[:a + 9] + "\n".join(lines) + "\n" + text[b:])

    with open(os.path.join(HERE, "cov_diag.txt"), "w") as f:
        f.write("# tests/data_covariance: diagonal(e_i^2) -- must give the chi2 of no covariance\n")
        header(f)
        f.write("diagonal\n")
        numbers(f, [e * e for e in err])

    def dense(f, relative):
        header(f)
        f.write("diagonal\n")
        numbers(f, [0.5 * e * e for e in err])
        for E in sorted(set(ener)):
            idx = [i for i in range(N) if ener[i] == E]
            f.write("block %d\n" % len(idx))
            f.write(" ".join(str(i + 1) for i in idx) + "\n")
            vals = []
            for i in idx:
                for j in idx:
                    vals.append(0.0 if i == j else 0.4 * err[i] * err[j])
            numbers(f, vals)
        for path in COV_FILES:
            v = [0.0] * N
            for i in range(N):
                if rows[i][0] == path:
                    v[i] = 0.02 if relative else 0.02 * val[i]
            f.write("relvector\n" if relative else "vector\n")
            numbers(f, v)
        f.write("vector  # a smooth global term\n")
        numbers(f, [0.01 * val[i] * math.sin(3.0 * ener[i]) for i in range(N)])

    with open(os.path.join(HERE, "cov_dense.txt"), "w") as f:
        f.write("# tests/data_covariance: dense absolute covariance (see make_models.py)\n")
        dense(f, False)
    with open(os.path.join(HERE, "cov_rel.txt"), "w") as f:
        f.write("# tests/data_covariance: as cov_dense.txt with the per-file 2 % terms model-scaled (relvector)\n")
        dense(f, True)


if __name__ == "__main__":
    main()
