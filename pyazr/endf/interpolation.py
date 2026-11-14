"""The interpolation laws of ENDF-6 (ENDF-102, 2023, Table 16 and eqs. 1-6, pp. 25-26).

INT=1 y constant in x (the value at the lower end of the interval), 2 linear-linear, 3 y linear in ln x, 4 ln y
linear in x, 5 ln y linear in ln x, 6 the charged-particle law between two points,
sigma(E) = (sigma_2 E_2)^alpha (sigma_1 E_1)^(1 - alpha) / E with alpha = (x(E) - x_1)/(x_2 - x_1),
x(E) = 1/sqrt(E - T), T = 0 for Q > 0 and the threshold for Q <= 0 (eqs. 2, 5, 6).  The manual leaves where T
comes from to the section: a caller gives it (``threshold``).  The two-dimensional schemes (11-15 corresponding
points, 21-25 unit base, pp. 26-28) are interpolations between distributions, not values: refused here.

A value outside the table is NaN (never extrapolated, never zero).  Where x repeats (a discontinuity, p. 25) the
value at that x is the one of the later point (the interval that starts there).  INT=4 and 5 with a zero end
of the interval give 0 inside it (the law's limit), with a change of sign NaN (the logarithm has no value).
"""

import numpy as np

from ..errors import FormatRefused

LAWS = {1: "constant (histogram)", 2: "linear-linear", 3: "linear in ln(x)", 4: "ln(y) linear in x",
        5: "ln(y) linear in ln(x)", 6: "charged-particle law (eq. 6)"}


def _law_of(interval, nbt, laws):
    """The law of the interval [x(n), x(n+1)] (n 0-based): the range m whose NBT(m) is the first above n+1."""
    for bound, law in zip(nbt, laws):
        if interval + 1 < bound:
            return law
    return laws[-1]


def _one(x1, x2, y1, y2, x, law, threshold):
    if law == 1:
        return y1
    if x2 == x1:
        return y2
    if law == 2:
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    if law == 3:
        return y1 + (y2 - y1) * np.log(x / x1) / np.log(x2 / x1)
    if law in (4, 5):
        # ln y is undefined at y = 0: the law's limit inside the interval is 0; a change of sign has no value
        if y1 == 0.0 or y2 == 0.0:
            return 0.0
        if y1 * y2 < 0.0:
            return np.nan
        t = (x - x1) / (x2 - x1) if law == 4 else np.log(x / x1) / np.log(x2 / x1)
        return y1 * np.exp(np.log(y2 / y1) * t)
    if law == 6:
        if y1 == 0.0 or y2 == 0.0:  # p. 26: "the result is clearly zero"
            return 0.0
        t = 0.0 if threshold is None else threshold
        if x1 <= t:
            return np.nan
        u1, u2, u = (1.0 / np.sqrt(v - t) for v in (x1, x2, x))
        alpha = (u - u1) / (u2 - u1)
        return (y2 * x2) ** alpha * (y1 * x1) ** (1.0 - alpha) / x
    raise FormatRefused(f"interpolation law INT={law}: not an interpolation of values",
                        what=f"interpolation law INT={law}", allowed="INT 1-6 (Table 16, p. 26)",
                        fix="laws 11-15 and 21-25 interpolate between distributions; read the record that holds them")


def interpolate(x, xs, ys, nbt, laws, threshold=None):
    """y at ``x`` (a number or an array) from a table: points ``xs``, ``ys`` (xs non-decreasing), interpolation
    ranges ``nbt`` (1-based last point of each range) and ``laws`` (INT of each range, p. 25).  ``threshold``: T of
    INT=6 (eq. 2), in the units of x.  NaN outside [xs[0], xs[-1]]."""
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    for law in laws:
        if law not in LAWS:
            raise FormatRefused(f"interpolation law INT={law} is not one of 1-6",
                                what=f"interpolation law INT={law}", allowed="INT 1-6 (Table 16, p. 26)",
                                fix="laws 11-15 and 21-25 interpolate between distributions, not values")
    scalar = np.ndim(x) == 0
    out = []
    for v in np.atleast_1d(np.asarray(x, float)):
        if not xs.size or not (xs[0] <= v <= xs[-1]) or np.isnan(v):
            out.append(np.nan)
            continue
        n = int(np.searchsorted(xs, v, side="right")) - 1  # the last point at or below v
        if n >= xs.size - 1:
            out.append(float(ys[-1]))
            continue
        if xs[n] == v:
            out.append(float(ys[n]))
            continue
        out.append(float(_one(xs[n], xs[n + 1], ys[n], ys[n + 1], v, _law_of(n, nbt, laws), threshold)))
    return out[0] if scalar else np.asarray(out)
