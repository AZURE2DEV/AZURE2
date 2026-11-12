#!/usr/bin/env python3
"""A cross section relative to a tabulated reference, sigma(E) / sigma_ref(E) (quantity ratio_to_reference), from
the engine's own sigma output and the table read here (mpmath; no AZURE2 code).

sigma_ref is models/data/ref_cm.dat (E_cm in MeV, sigma in mb), interpolated linearly in ln sigma against E_cm
between the two rows around E (the engine interpolates in the lab energy, a linear function of E_cm: the same
interpolant).

  truth.py check <AZUREOut .extrap of models/extrap.azr>
      blocks: 1 sigma, 2 the ratio to ref_cm.dat, 3 the ratio to a constant 1 b.  Prints the largest relative
      difference of block 2's value from sigma / sigma_ref and whether block 3 equals block 1 character for
      character; exit 1 beyond 1e-10 or when the text differs.
  truth.py data
      writes models/data/ratio.dat: the rows of sigma.dat (E_cm, value, error) divided by sigma_ref(E_cm) in b.
"""
import os
import sys

import mpmath as mp

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))


def table():
    rows = []
    for line in open(os.path.join(HERE, "models", "data", "ref_cm.dat")):
        line = line.split("#")[0].split()
        if len(line) >= 2:
            rows.append((mp.mpf(line[0]), mp.mpf(line[1]) * mp.mpf("1e-3")))  # mb -> b
    return rows


def reference(rows, e):
    e = mp.mpf(e)
    for (e0, s0), (e1, s1) in zip(rows, rows[1:]):
        if e0 <= e <= e1:
            t = (e - e0) / (e1 - e0)
            return mp.exp((1 - t) * mp.log(s0) + t * mp.log(s1))
    raise ValueError("E = %s outside the table" % e)


def blocks(path):
    out, cur = [], []
    for line in open(path):
        if line.strip():
            cur.append(line.rstrip("\n"))
        elif cur:
            out.append(cur)
            cur = []
    if cur:
        out.append(cur)
    return out


def check(path):
    b, rows = blocks(path), table()
    if len(b) != 3 or len({len(x) for x in b}) != 1:
        print("expected 3 blocks of equal length, got", [len(x) for x in b])
        return 1
    worst = mp.mpf(0)
    for s_row, r_row in zip(b[0], b[1]):
        s, r = s_row.split(), r_row.split()
        if s[0] != r[0]:
            print("grids differ:", s[0], r[0])
            return 1
        want = mp.mpf(s[3]) / reference(rows, s[0])
        worst = max(worst, abs(mp.mpf(r[3]) / want - 1))
    same = b[2] == b[0]
    print("the ratio against sigma / sigma_ref: %.2e; the ratio to 1 b equals sigma: %s"
          % (worst, "yes" if same else "no"))
    return 0 if worst <= mp.mpf("1e-10") and same else 1


def data():
    rows = table()
    with open(os.path.join(HERE, "models", "data", "ratio.dat"), "w") as f:
        f.write("# E_cm [MeV]  (unused)  sigma / sigma_ref  error  (truth.py data: the rows of sigma.dat over ref_cm.dat)\n")
        for line in open(os.path.join(HERE, "models", "data", "sigma.dat")):
            if line.startswith("#") or not line.strip():
                continue
            e, _, v, dv = line.split()
            ref = reference(rows, e)
            f.write("%s  0  %s  %s\n" % (e, mp.nstr(mp.mpf(v) / ref, 17, min_fixed=0, max_fixed=0),
                                          mp.nstr(mp.mpf(dv) / ref, 17, min_fixed=0, max_fixed=0)))
    return 0


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "data":
        sys.exit(data())
    if len(sys.argv) != 3 or sys.argv[1] != "check":
        print(__doc__)
        sys.exit(2)
    sys.exit(check(sys.argv[2]))
