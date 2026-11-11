#!/usr/bin/env python3
"""Comparisons for tests/transmission_yield/check.sh.  Each prints one line "<max> <detail>"; the shell decides.

  formula  T|Y  <bare dir> <thin dir> N      max relative |model - formula(bare)| over the data points
  xformula T|Y  <bare dir> <thin dir> N      the same for test segments (mode 3, .extrap)
  conv     T|Y  <truth dir> <conv dir> N SIGMA_LAB M1 M2
           max relative |model - integral| with the integral of the transformed dense, unbroadened cross sections
           against the Gaussian kernel (+-5 sigma, normalised over the window, c.m. energies) at every data point
  wrong    T|Y  <truth dir> <bconv dir> N SIGMA_LAB M1 M2
           the same integral against the transform of the BROADENED cross sections (what must differ)
  natural  <nat dir> <natT dir> N                  max relative |T - exp(-N sigma_nat)|
"""
import math
import sys


def table(path):
    """Numeric rows of an output file; none if the file is missing (the caller prints "nan rows")."""
    rows = []
    try:
        lines = open(path).readlines()
    except OSError:
        return rows
    for line in lines:
        f = line.split()
        if len(f) >= 4:
            try:
                rows.append([float(x) for x in f])
            except ValueError:
                pass
    return rows


def out(d, ent, ex, ext="out"):
    return table("%s/output/AZUREOut_aa=%d_R=%d.%s" % (d, ent, ex, ext))


def y0(n, s, st):
    u = n * st
    return s * (n if abs(u) < 1e-12 else -math.expm1(-u) / st)


def transform(kind, n, s, st):
    return math.exp(-n * st) if kind == "T" else y0(n, s, st)


def maxrel(pairs):
    w, at = 0.0, None
    for e, got, want in pairs:
        d = abs(got - want) / abs(want)
        if d > w:
            w, at = d, (e, got, want)
    return w, at


def report(w, at, n):
    if at is None:
        print("%.3e n=%d" % (w, n))
    else:
        print("%.3e n=%d at E_cm %.5f got %.9g want %.9g" % (w, n, at[0], at[1], at[2]))


def formula(kind, bare, thin, n, ext="out"):
    st = out(bare, 1, 1, ext)
    sg = out(bare, 1, 2, ext)
    model = out(thin, 1, 1, ext) if kind == "T" else out(thin, 1, 2, ext)
    if not (len(st) == len(sg) == len(model) > 0):
        print("nan rows %d %d %d" % (len(st), len(sg), len(model)))
        return
    pairs = [(m[0], m[3], transform(kind, n, g[3], t[3])) for m, t, g in zip(model, st, sg)]
    report(*maxrel(pairs), len(pairs))


def kernel_avg(grid, values, e0, sigma):
    lo, hi = e0 - 5 * sigma, e0 + 5 * sigma
    num = den = 0.0
    prev = None
    for e, v in zip(grid, values):
        if e < lo or e > hi:
            prev = None
            continue
        g = math.exp(-0.5 * ((e - e0) / sigma) ** 2)
        if prev is not None:
            h = e - prev[0]
            num += 0.5 * h * (g * v + prev[1] * prev[2])
            den += 0.5 * h * (g + prev[1])
        prev = (e, g, v)
    return num / den


def conv(kind, truth, model_dir, n, sigma_lab, m1, m2, wrong=False):
    st = out(truth, 1, 1, "extrap")
    sg = out(truth, 1, 2, "extrap")
    grid = [r[0] for r in st]
    sigma = sigma_lab * m2 / (m1 + m2)
    model = out(model_dir, 1, 1) if kind == "T" else out(model_dir, 1, 2)
    if not model or len(st) != len(sg) or len(st) < 1000:
        print("nan rows %d %d %d" % (len(model), len(st), len(sg)))
        return
    if not wrong:
        f = [transform(kind, n, g[3], t[3]) for t, g in zip(st, sg)]
        pairs = [(m[0], m[3], kernel_avg(grid, f, m[0], sigma)) for m in model]
    else:
        # model_dir holds the broadened sigma_t (1,1) and sigma_gamma (1,2); the truth is the integral of the transform
        bt, bg = out(model_dir, 1, 1), out(model_dir, 1, 2)
        f = [transform(kind, n, g[3], t[3]) for t, g in zip(st, sg)]
        pairs = [(t[0], transform(kind, n, g[3], t[3]), kernel_avg(grid, f, t[0], sigma)) for t, g in zip(bt, bg)]
    report(*maxrel(pairs), len(pairs))


def reaction_yield(bare, elast, thin, n):
    """Y0 of n + 20Ne -> alpha + 17O against the formula with sigma_t = "total 1", and how far the elastic is from it."""
    st, sx, el, y = out(bare, 2, 2), out(bare, 2, 1), out(elast, 2, 2), out(thin, 2, 1)
    if not (len(st) == len(sx) == len(el) == len(y) > 0):
        print("nan rows %d %d %d %d" % (len(st), len(sx), len(el), len(y)))
        return
    pairs, zero_bad, open_ = [], 0, 0
    for t, x, m in zip(st, sx, y):
        want = y0(n, x[3], t[3])
        if want == 0.0:
            zero_bad += (m[3] != 0.0)
            continue
        open_ += 1
        pairs.append((m[0], m[3], want))
    w, at = maxrel(pairs)
    sep = max(abs(e[3] - t[3]) / t[3] for e, t in zip(el, st))
    if zero_bad:
        w = float("inf")
    print("%.3e n=%d elastic_vs_total %.3e" % (w, open_, sep))


def natural(nat, natT, n):
    s = out(nat, 2, 2)
    t = out(natT, 2, 2)
    if not (len(s) == len(t) > 0):
        print("nan rows %d %d" % (len(s), len(t)))
        return
    report(*maxrel([(a[0], b[3], math.exp(-n * a[3])) for a, b in zip(s, t)]), len(s))


if __name__ == "__main__":
    what = sys.argv[1]
    a = sys.argv[2:]
    if what == "formula":
        formula(a[0], a[1], a[2], float(a[3]))
    elif what == "xformula":
        formula(a[0], a[1], a[2], float(a[3]), "extrap")
    elif what in ("conv", "wrong"):
        conv(a[0], a[1], a[2], float(a[3]), float(a[4]), float(a[5]), float(a[6]), wrong=(what == "wrong"))
    elif what == "reaction":
        reaction_yield(a[0], a[1], a[2], float(a[3]))
    elif what == "natural":
        natural(a[0], a[1], float(a[2]))
