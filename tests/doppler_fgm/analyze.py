#!/usr/bin/env python3
"""Comparisons for tests/doppler_fgm/check.sh.  Each prints one line "<max> <detail>"; the shell decides.

The independent free-gas integral (SAMMY manual, "Free-Gas Model of Doppler Broadening"), on the 1 eV truth grid:
  sigma_D(E) = 1/(u sqrt(pi) v^2) int [exp(-(v-w)^2/u^2) - exp(-(v+w)^2/u^2)] w^2 sigma(w^2) dw,
v = sqrt(E_L), w = sqrt(E'_L) in lab-equivalent energies E_L = E_cm (M1+M2)/M2, u^2 = kT M1/M2; trapezoid in w.

  doppler  K T|G <truth> <dop dir> T_K M1 M2           sigma_D at the data points
  dopgauss K T|G <truth> <dir> T_K M1 M2 SIGMA_LAB     Gaussian average (c.m., +-5 sigma, normalised) of sigma_D
  dopT     K T|Y <truth> <dir> T_K M1 M2 N [SIGMA_LAB]   transform of sigma_D at the point, or the Gaussian average
                                                       of the transform of sigma_D (Doppler, transform, kernel)
  effect   <bare dir> <dop dir>                        max |sigma_D/sigma - 1| over sigma_t (how much the step does)
  order    <truth> <dopg dir> T_K M1 M2 N SIGMA_LAB    max relative |exp(-N <sigma_D>_G) - <exp(-N sigma_D)>_G|
"""
import bisect
import math
import sys

KB = 8.617333262e-11  # MeV/K


def table(path):
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


def report(w, at, n):
    if at is None:
        print("%.3e n=%d" % (w, n))
    else:
        print("%.3e n=%d at E_cm %.7f got %.9g want %.9g" % (w, n, at[0], at[1], at[2]))


def maxrel(pairs):
    w, at = float("nan"), None
    for e, got, want in pairs:
        d = abs(got - want) / abs(want)
        if at is None or d > w:
            w, at = d, (e, got, want)
    return w, at


class Truth:
    def __init__(self, d, t_k, m1, m2):
        st, sg = out(d, 1, 1, "extrap"), out(d, 1, 2, "extrap")
        self.ok = len(st) > 1000 and len(st) == len(sg)
        self.c = m2 / (m1 + m2)
        self.u = math.sqrt(KB * t_k * m1 / m2)
        self.ecm = [r[0] for r in st]
        self.st = [r[3] for r in st]
        self.sg = [r[3] for r in sg]

    def broaden(self, sig, ecm):
        """The free-gas integral at c.m. energy ecm, over truth points within 6 u of v (the kernel is < e^-36
        beyond).  Returns nan if that window is not inside the truth grid."""
        v = math.sqrt(ecm / self.c)
        u = self.u
        if not hasattr(self, "w"):
            self.w = [math.sqrt(e / self.c) for e in self.ecm]
        i0 = bisect.bisect_left(self.w, v - 6 * u)
        i1 = bisect.bisect_right(self.w, v + 6 * u)
        if i0 == 0 or i1 >= len(self.w):
            return float("nan")
        num, prev = 0.0, None
        for i in range(i0 - 1, i1 + 1):
            w, s = self.w[i], sig[i]
            f = (math.exp(-((v - w) / u) ** 2) - math.exp(-((v + w) / u) ** 2)) * w * w * s
            if prev is not None:
                num += 0.5 * (w - prev[0]) * (f + prev[1])
            prev = (w, f)
        return num / (u * math.sqrt(math.pi) * v * v)


def gauss_avg(grid, values, e0, sigma):
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


def y0(n, s, st):
    u = n * st
    return s * (n if abs(u) < 1e-12 else -math.expm1(-u) / st)


def main():
    what, a = sys.argv[1], sys.argv[2:]
    if what == "effect":
        b, d = out(a[0], 1, 1), out(a[1], 1, 1)
        if not b or len(b) != len(d):
            print("nan rows %d %d" % (len(b), len(d)))
            return
        report(*maxrel([(x[0], y[3], x[3]) for x, y in zip(b, d)]), len(b))
        return
    if what == "order":
        tr = Truth(a[0], float(a[2]), float(a[3]), float(a[4]))
        n, sig_lab = float(a[5]), float(a[6])
        model = out(a[1], 1, 1)
        if not tr.ok or not model:
            print("nan rows")
            return
        sd = [tr.broaden(tr.st, e) for e in tr.ecm[::5]]
        grid = tr.ecm[::5]
        sigma = sig_lab * tr.c
        pairs = []
        for m in model:
            right = gauss_avg(grid, [math.exp(-n * s) for s in sd], m[0], sigma)
            wrong = math.exp(-n * gauss_avg(grid, sd, m[0], sigma))
            pairs.append((m[0], wrong, right))
        report(*maxrel(pairs), len(pairs))
        return
    kind = a[0]
    tr = Truth(a[1], float(a[3]), float(a[4]), float(a[5]))
    d = a[2]
    if not tr.ok:
        print("nan truth rows %d" % len(tr.ecm))
        return
    model = out(d, 1, 1) if kind in ("T",) else out(d, 1, 2)
    if not model:
        print("nan rows 0")
        return
    if what == "doppler":
        sig = tr.st if kind == "T" else tr.sg
        report(*maxrel([(m[0], m[3], tr.broaden(sig, m[0])) for m in model]), len(model))
    elif what == "dopgauss":
        sig = tr.st if kind == "T" else tr.sg
        sigma = float(a[6]) * tr.c
        lo, hi = model[0][0] - 6 * sigma, model[-1][0] + 6 * sigma
        idx = [i for i, e in enumerate(tr.ecm) if lo <= e <= hi][::2]
        grid = [tr.ecm[i] for i in idx]
        sd = [tr.broaden(sig, e) for e in grid]
        report(*maxrel([(m[0], m[3], gauss_avg(grid, sd, m[0], sigma)) for m in model]), len(model))
    elif what == "dopT":
        n = float(a[6])
        f = (lambda st, sg: math.exp(-n * st)) if kind == "T" else (lambda st, sg: y0(n, sg, st))
        if len(a) > 7:
            sigma = float(a[7]) * tr.c
            lo, hi = model[0][0] - 6 * sigma, model[-1][0] + 6 * sigma
            idx = [i for i, e in enumerate(tr.ecm) if lo <= e <= hi][::2]
            grid = [tr.ecm[i] for i in idx]
            vals = [f(tr.broaden(tr.st, e), tr.broaden(tr.sg, e)) for e in grid]
            report(*maxrel([(m[0], m[3], gauss_avg(grid, vals, m[0], sigma)) for m in model]), len(model))
        else:
            report(*maxrel([(m[0], m[3], f(tr.broaden(tr.st, m[0]), tr.broaden(tr.sg, m[0]))) for m in model]),
                   len(model))


if __name__ == "__main__":
    main()
