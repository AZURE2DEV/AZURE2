#!/usr/bin/env python3
"""Numerical comparisons for tests/light_excitation/check.sh.  Every command prints
one line starting with a number (or exits non-zero), so an empty or crashed
analysis is a FAIL in check.sh, never a pass.

  analyze.py compare <dirA> <dirB> <glob>...
      max relative difference of every numeric token of the matching files of
      two output directories (same files, same lines, same token counts or exit 1);
      prints "<maxrel> <files> <numbers>"
  analyze.py symmetry <extrap file>
      c.m. differential rows (E_cm, E_x, theta, sigma, ...): max |s(t) - s(180-t)| / s
      over every energy and angle pair; prints "<maxrel> <pairs>"
  analyze.py channels <run.log> <identical 0|1> <j> <maxL> <levelJ> <levelPi>
      the "Hard-sphere group added: J = ... with N channels" lines against an
      independent enumeration of the (l, s) of a pair of two spin-j particles
      of positive parity (l + s even when identical); prints "<mismatches> <groups> <channels>"
"""
import glob
import math
import os
import re
import sys


def tokens(path):
    rows = []
    with open(path) as f:
        for ln in f:
            rows.append(ln.split())
    return rows


def isnum(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def compare(a, b, patterns):
    files = sorted({os.path.basename(p) for pat in patterns for p in glob.glob(os.path.join(a, pat))})
    if not files:
        sys.exit("no files")
    worst, n = 0.0, 0
    for name in files:
        pb = os.path.join(b, name)
        if not os.path.exists(pb):
            sys.exit(f"{name} missing from {b}")
        ra, rb = tokens(os.path.join(a, name)), tokens(pb)
        if len(ra) != len(rb):
            sys.exit(f"{name}: {len(ra)} vs {len(rb)} lines")
        for la, lb in zip(ra, rb):
            if len(la) != len(lb):
                sys.exit(f"{name}: line {la} vs {lb}")
            for x, y in zip(la, lb):
                if isnum(x) and isnum(y):
                    x, y = float(x), float(y)
                    if not (math.isfinite(x) and math.isfinite(y)):
                        sys.exit(f"{name}: non-finite value")
                    d = abs(x - y) / max(abs(x), abs(y), 1e-300) if (x or y) else 0.0
                    worst = max(worst, d)
                    n += 1
                elif x != y:
                    sys.exit(f"{name}: token {x!r} vs {y!r}")
    if n == 0:
        sys.exit("no numbers")
    print(f"{worst:.3e} {len(files)} {n}")


def symmetry(path):
    by_e = {}
    for row in tokens(path):
        if len(row) < 4 or not all(isnum(t) for t in row[:4]):
            continue
        e, th, s = float(row[0]), float(row[2]), float(row[3])
        by_e.setdefault(round(e, 9), {})[round(th, 6)] = s
    worst, pairs = 0.0, 0
    for e, d in by_e.items():
        for th, s in d.items():
            if th >= 90:
                continue
            m = d.get(round(180 - th, 6))
            if m is None:
                continue
            if not (s > 0 and m > 0 and math.isfinite(s) and math.isfinite(m)):
                sys.exit(f"bad value at E = {e}, theta = {th}")
            worst = max(worst, abs(s - m) / max(s, m))
            pairs += 1
    if pairs == 0:
        sys.exit("no angle pairs")
    print(f"{worst:.3e} {pairs}")


def channels(log, identical, j, maxl, level_j, level_pi):
    j, maxl = float(j), int(maxl)
    identical = identical == "1"
    want = {}
    s_values = [abs(j - j) + k for k in range(int(round(2 * j)) + 1)]
    for s in s_values:
        for l in range(maxl + 1):
            pi = 1 if l % 2 == 0 else -1
            J = abs(s - l)
            while J <= s + l + 1e-6:
                want.setdefault((J, pi), set())
                J += 1.0
    for (J, pi) in list(want):
        for s in s_values:
            for l in range(maxl + 1):
                if (1 if l % 2 == 0 else -1) != pi or not (abs(s - l) <= J + 1e-6 <= s + l + 2e-6):
                    continue
                if identical and (l + int(round(s))) % 2:
                    continue
                want[(J, pi)].add((l, s))
    want.pop((float(level_j), int(level_pi)), None)      # the level's own group: not a hard-sphere group
    got = {}
    pat = re.compile(r"Hard-sphere group added: J = ([0-9.]+)([+-]) with (\d+) channels")
    for ln in open(log):
        m = pat.search(ln)
        if m:
            got[(float(m.group(1)), 1 if m.group(2) == "+" else -1)] = int(m.group(3))
    if not got:
        sys.exit("no hard-sphere groups in the log")
    bad = 0
    for key in set(want) | set(got):
        w = len(want.get(key, ()))
        if w == 0 and key not in got:
            continue
        if got.get(key, -1) != w:
            print(f"# J^pi {key}: {got.get(key)} channels, expected {w}", file=sys.stderr)
            bad += 1
    print(f"{bad} {len(got)} {sum(got.values())}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "compare" and len(sys.argv) >= 5:
        compare(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif cmd == "symmetry" and len(sys.argv) == 3:
        symmetry(sys.argv[2])
    elif cmd == "channels" and len(sys.argv) == 8:
        channels(*sys.argv[2:])
    else:
        sys.exit(__doc__)
