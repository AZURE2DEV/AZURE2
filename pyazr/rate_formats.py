"""A rate of :mod:`pyazr.rate` in the formats of the rate libraries: STARLIB's table and REACLIB's parameter sets.

Sources of the layouts (the two documents were read for this module; nothing here is from memory):

- **REACLIB**: "REACLIB File Format", the help page of the JINA REACLIB Database (reaclib.jinaweb.org, help topic
  ``reaclib_format``, last updated 2010-11-29): the chapters (its Table 1), the REACLIB 1 file (a header entry of
  three lines, the chapter number on the first and two blank lines, then the set entries), every line 74
  characters padded with blanks; a set entry is three lines, ``format(i1,4x,6a5,8x,a4,a1,a1,3x,1pe12.5)`` (five
  blanks, six nuclide fields of five characters, eight blanks, the four-character set label, the resonance flag
  (blank or ``n`` non-resonant, ``r`` resonant, ``w`` weak), the reverse flag (``v``), three blanks, the Q value
  in twelve characters, ten blanks), ``format(4e13.6)`` (a0 to a3, 22 blanks), ``format(3e13.6)`` (a4 to a6,
  35 blanks).  The rate of a set is ``exp(a0 + sum_{i=1..5} a_i T9^((2i-5)/3) + a6 ln T9)``; a rate is the sum of
  its sets (Cyburt et al. 2010, ApJS 189, 240, eq. 1; the help page shows it as an image).  The numbers are written
  as ``%13.6e`` (`` 1.234567e+00``), which ``e13.6`` reads.
- **STARLIB**: Sallaska et al. 2013, ApJS 207, 18 (arXiv:1304.7811), Appendix A "STARLIB Format": per interaction
  one header line of 64 fields: the REACLIB chapter left-aligned in fields 1-2, the species in fields 6-35 (five
  fields each, right-aligned, entrance channel first, lightest to heaviest, element symbol then mass number, ``n``,
  ``p``, ``t`` for the neutron, proton and triton, ``d`` for the deuteron as in the appendix's example), the
  reference label right-aligned in fields 44-47, the interaction label in field 48, the Q value (MeV) in fields
  53-64; then three columns: T9 (GK), the recommended rate (cm^3 mol^-1 s^-1) and its factor uncertainty f.u. =
  exp(sigma) of a lognormal rate (x_low = e^(mu - sigma), x_med = e^mu, x_high = e^(mu + sigma), 68 %).  The grid
  has 60 temperatures from 1 MK to 10 GK; the paper does not print them: :data:`STARLIB_T9` is 0.001 to 0.009 GK
  in steps of 0.001 and the 51 temperatures of the paper's rate tables (its Appendix C, 0.010 to 10 GK), which
  make 60.  The paper gives no number format for the columns: they are written ``%.3e %.6e %.4e``.

The factor uncertainty comes from a band of :func:`pyazr.rate.rate_band` and nothing else: from a covariance band,
the recommended rate is the rate at the parameters and f.u. = exp(sigma) of ln N_A<sigma v> (``high / best``); from
draws, the recommended rate is the median (50th percentile) and f.u. = (x84 / x16)^(1/2).  Without a band the
table is not written.  The species, the chapter and the Q value come from the project's pairs (``Q`` = the entrance
pair's threshold minus the exit pair's, the thresholds being the pairs' separation energies relative to the
compound nucleus plus their excitation energies; the total capture, exit -1, ends in the compound nucleus's ground
state).  The labels a user chooses: a four-character reference / set label and the flags, among the allowed values.
"""
import math

import numpy as np

__all__ = ["STARLIB_T9", "REACLIB_RANGE", "species_of", "write_starlib", "read_starlib", "fit_reaclib",
           "write_reaclib", "read_reaclib", "reaclib_rate", "ReaclibFit"]

#: STARLIB's 60 temperatures (GK); see the module's text for where they come from.
STARLIB_T9 = [round(0.001 * k, 3) for k in range(1, 10)] + [
    0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.018, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
    0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

#: The temperatures (GK) a REACLIB set must reproduce: the computed table must cover them.
REACLIB_RANGE = (0.01, 10.0)

#: The allowed values of the labels a user chooses: REACLIB's resonance flag, STARLIB's interaction label.
REACLIB_FLAGS = ("n", "r", "w")
STARLIB_FLAGS = ("", "w", "g", "v")


def _refusal(what, allowed, fix):
    from .errors import ParameterError
    return ParameterError(f"{what}; allowed: {allowed}; {fix}", what=what, allowed=allowed, fix=fix)


def _name(Z, A):
    from .system.nuclei import SYMBOLS
    special = {(0, 1): "n", (1, 1): "p", (1, 2): "d", (1, 3): "t"}
    if (Z, A) in special:
        return special[(Z, A)]
    return f"{SYMBOLS[Z].lower()}{A}"


def species_of(file, entrance, exit):
    """``{"chapter", "reactants", "products", "Q"}`` of the reaction from pair ``entrance`` to pair ``exit`` of the
    project ``file``: the libraries' names (lightest first), the REACLIB chapter (4: two reactants, one product; 5:
    two and two), Q (MeV) from the pairs' thresholds.  A beta pair, a photon entrance or an unknown pair is refused."""
    from .project import Project
    pairs = {q.key: q for q in Project.open(file).pairs()}
    if entrance not in pairs or pairs[entrance].kind != 0:
        raise _refusal(f"entrance {entrance} is not a particle pair of the project",
                       ", ".join(str(k) for k, q in pairs.items() if q.kind == 0), "choose a particle pair")
    a = pairs[entrance]
    reactants = sorted([(a.light, a.light_mass), (a.heavy, a.heavy_mass)], key=lambda x: x[1])
    if exit == -1:
        Z, A = a.light[0] + a.heavy[0], a.light[1] + a.heavy[1]
        products, q = [(Z, A)], a.threshold
    else:
        if exit not in pairs or pairs[exit].kind not in (0, 10):
            raise _refusal(f"exit {exit} is not a particle or photon pair of the project",
                           ", ".join(str(k) for k, q in pairs.items() if q.kind in (0, 10)) + ", -1",
                           "choose a particle or photon pair, or -1 for the total capture")
        b = pairs[exit]
        q = a.threshold - b.threshold
        if b.kind == 10:
            products = [b.heavy]
        else:
            products = [x for x, _ in sorted([(b.light, b.light_mass), (b.heavy, b.heavy_mass)], key=lambda x: x[1])]
    names_in = [_name(*z) for z, _ in reactants]
    names_out = [_name(*z) for z in products]
    return {"chapter": 4 if len(names_out) == 1 else 5, "reactants": names_in, "products": names_out, "Q": float(q)}


# -- STARLIB ------------------------------------------------------------------------------------------------------

def _starlib_columns(rate, band):
    if band is None:
        raise _refusal("no band: STARLIB's third column is the factor uncertainty of the rate",
                       "a band of pyazr.rate.rate_band (from a fit's covariance or a sampler's draws)",
                       "compute the rate with a band; nothing is written without one")
    if not np.array_equal(np.asarray(rate.T9, float), np.asarray(STARLIB_T9, float)):
        raise _refusal("the rate is not on STARLIB's 60 temperatures", "grid='starlib' (pyazr.rate)",
                       "compute the rate on grid='starlib'")
    if band.method == "covariance":
        return np.asarray(band.best, float), np.asarray(band.high, float) / np.asarray(band.best, float)
    return np.asarray(band.median, float), np.sqrt(np.asarray(band.high, float) / np.asarray(band.low, float))


def write_starlib(path, rate, band, species, label="azr2", flag=""):
    """Write ``rate`` (on :data:`STARLIB_T9`) with the factor uncertainty of ``band`` as one STARLIB entry:
    ``species`` from :func:`species_of`, ``label`` the reference label (at most four characters), ``flag`` the
    interaction label (one of :data:`STARLIB_FLAGS`).  Returns ``path``."""
    if len(label) > 4 or not label.strip():
        raise _refusal(f"reference label {label!r}", "one to four characters", "give a shorter label")
    if flag not in STARLIB_FLAGS:
        raise _refusal(f"interaction label {flag!r}", ", ".join(repr(f) for f in STARLIB_FLAGS),
                       "w weak, g gamma-ray transition, v from detailed balance, empty otherwise")
    recommended, fu = _starlib_columns(rate, band)
    names = species["reactants"] + species["products"]
    if len(names) > 6:
        raise _refusal(f"{len(names)} species", "at most six (fields 6-35)", "this reaction has too many")
    head = [" "] * 64
    head[0:2] = f"{species['chapter']:<2d}"
    for k, n in enumerate(names):
        head[5 + 5 * k:10 + 5 * k] = f"{n:>5}"[:5]
    head[43:47] = f"{label:>4}"
    head[47] = flag or " "
    head[52:64] = f"{species['Q']:12.5e}"
    lines = ["".join(head)]
    lines += [f"{t:.3e} {r:.6e} {f:.4e}" for t, r, f in zip(rate.T9, recommended, fu)]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def read_starlib(path, entry=0):
    """The entry ``entry`` (0-based) of a STARLIB file: ``{"header", "chapter", "species", "Q", "label", "T9",
    "rate", "fu"}``.  An entry is a header line and the three-column lines after it."""
    entries, current = [], None
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            t = line.split()
            try:
                row = [float(v) for v in t[:3]]
            except ValueError:
                row = None
            if row is None or len(t) != 3:
                current = {"header": line.rstrip("\n"), "rows": []}
                entries.append(current)
            elif current is not None:
                current["rows"].append(row)
    if not entries or entry >= len(entries):
        raise ValueError(f"{path}: no STARLIB entry {entry} ({len(entries)} found)")
    e = entries[entry]
    h = e["header"].ljust(64)
    rows = np.array(e["rows"], float)
    return {"header": e["header"], "chapter": int(h[0:2]), "species": h[5:35].split(), "label": h[43:47].strip(),
            "flag": h[47].strip(), "Q": float(h[52:64]), "T9": rows[:, 0], "rate": rows[:, 1], "fu": rows[:, 2],
            "path": path}


# -- REACLIB ------------------------------------------------------------------------------------------------------

def _basis(T9):
    t = np.asarray(T9, float)
    return np.stack([np.ones_like(t), 1 / t, t ** (-1 / 3), t ** (1 / 3), t, t ** (5 / 3), np.log(t)], axis=1)


def reaclib_rate(sets, T9):
    """N_A<sigma v> of REACLIB ``sets`` (each seven numbers a0..a6) at ``T9`` (GK): the sum over the sets."""
    B = _basis(T9)
    return sum(np.exp(B @ np.asarray(a, float)) for a in sets)


class ReaclibFit:
    """REACLIB sets fitted to a computed rate: ``sets`` (each a0..a6), ``T9`` and ``rate`` (the table fitted),
    ``deviation`` (fit / table - 1 at each temperature), ``max_deviation``, ``tolerance``, ``range`` (GK)."""

    def __init__(self, sets, T9, rate, tolerance):
        self.sets = [np.asarray(a, float) for a in sets]
        self.T9, self.rate, self.tolerance = np.asarray(T9, float), np.asarray(rate, float), tolerance
        self.deviation = reaclib_rate(self.sets, self.T9) / self.rate - 1.0
        self.max_deviation = float(np.max(np.abs(self.deviation)))
        self.range = (float(self.T9.min()), float(self.T9.max()))

    @property
    def n_sets(self):
        return len(self.sets)

    def summary(self):
        worst = int(np.argmax(np.abs(self.deviation)))
        return (f"{self.n_sets} set(s) fitted to {self.T9.size} temperatures from {self.range[0]:g} to "
                f"{self.range[1]:g} GK: largest deviation {100 * self.max_deviation:.3g} % at T9 = "
                f"{self.T9[worst]:g} (tolerance {100 * self.tolerance:g} %)")


def _fit_sets(lnR, B, n, t):
    """``n`` sets fitted to ln R: one by linear least squares; several by scipy's least squares in ln of their sum,
    started from linear fits of ``n`` contiguous pieces of the temperatures (each piece's set fitted to its own
    temperatures only)."""
    if n == 1:
        return [np.linalg.lstsq(B, lnR, rcond=None)[0]]
    from scipy.optimize import least_squares
    pieces = np.array_split(np.arange(t.size), n)
    start = []
    for idx in pieces:
        sub = B[idx][:, [0, 1, 6]] if idx.size < 7 else B[idx]
        a = np.linalg.lstsq(sub, lnR[idx], rcond=None)[0]
        full = np.zeros(7)
        full[[0, 1, 6] if idx.size < 7 else slice(None)] = a
        start.append(full)
    x0 = np.concatenate(start)

    def resid(x):
        terms = B @ x.reshape(n, 7).T
        m = terms.max(axis=1)
        return m + np.log(np.exp(terms - m[:, None]).sum(axis=1)) - lnR

    r = least_squares(resid, x0, method="lm", max_nfev=20000 * n, xtol=1e-14, ftol=1e-14)
    return list(r.x.reshape(n, 7))


def fit_reaclib(rate, tolerance=0.05, max_sets=3):
    """REACLIB sets reproducing ``rate`` (a :class:`pyazr.rate.Rate`) at every computed temperature in
    :data:`REACLIB_RANGE` to within ``tolerance`` (relative): one set, then two, ... up to ``max_sets``; the first
    number of sets that reaches the tolerance is returned (:class:`ReaclibFit`).  Refused (ParameterError with the
    three parts) when the computed temperatures do not cover the range (no extrapolation), when a rate there is not
    positive, or when ``max_sets`` sets do not reach the tolerance (the best fit's deviation is said)."""
    t = np.asarray(rate.T9, float)
    lo, hi = REACLIB_RANGE
    if t.size == 0 or t.min() > lo or t.max() < hi:
        raise _refusal(f"the rate is computed from {t.min():g} to {t.max():g} GK; a REACLIB set stands for "
                       f"{lo:g} to {hi:g} GK", f"a rate computed over {lo:g} to {hi:g} GK",
                       "compute it on grid='table' or 'starlib' (pyazr.rate); nothing is extrapolated")
    keep = (t >= lo) & (t <= hi)
    t, R = t[keep], np.asarray(rate.rate, float)[keep]
    if np.any(~(R > 0)):
        raise _refusal(f"the rate is not positive at T9 = {t[~(R > 0)][0]:g} GK", "a rate > 0 at every temperature",
                       "a rate that underflows cannot be fitted in ln; restrict the reaction or the temperatures")
    lnR, B = np.log(R), _basis(t)
    best = None
    for n in range(1, max_sets + 1):
        try:
            fit = ReaclibFit(_fit_sets(lnR, B, n, t), t, R, tolerance)
        except (ImportError, np.linalg.LinAlgError, ValueError, FloatingPointError):
            continue
        if not math.isfinite(fit.max_deviation):
            continue
        if best is None or fit.max_deviation < best.max_deviation:
            best = fit
        if fit.max_deviation <= tolerance:
            return fit
    worst = "none could be fitted" if best is None else best.summary()
    raise _refusal(f"{max_sets} REACLIB set(s) do not reproduce the rate within {100 * tolerance:g} %: {worst}",
                   f"a largest deviation <= {100 * tolerance:g} %",
                   "allow more sets (max_sets) or a larger tolerance, and say so where the sets are used")


def write_reaclib(path, fit, species, label="azr2", flag="n"):
    """Write the sets of ``fit`` (:class:`ReaclibFit`) as a REACLIB 1 file: the chapter's header entry, then one set
    entry per set; ``species`` from :func:`species_of`, ``label`` the set label (four characters), ``flag`` the
    resonance flag (one of :data:`REACLIB_FLAGS`).  Returns ``path``."""
    if len(label) != 4:
        raise _refusal(f"set label {label!r}", "four characters", "give a four-character label")
    if flag not in REACLIB_FLAGS:
        raise _refusal(f"resonance flag {flag!r}", ", ".join(REACLIB_FLAGS), "n non-resonant, r resonant, w weak")
    names = species["reactants"] + species["products"]
    blank = " " * 74
    lines = [f"{species['chapter']:<74d}", blank, blank]
    for a in fit.sets:
        first = " " * 5 + "".join(f"{n:>5}" for n in names).ljust(30) + " " * 8 + label + flag + " " + " " * 3 \
            + f"{species['Q']:12.5e}" + " " * 10
        lines += [first, "".join(f"{v:13.6e}" for v in a[:4]) + " " * 22,
                  "".join(f"{v:13.6e}" for v in a[4:]) + " " * 35]
    if any(len(ln) != 74 for ln in lines):
        raise RuntimeError("a REACLIB line is not 74 characters")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def read_reaclib(path):
    """The sets of a REACLIB 1 file: ``[{"chapter", "species", "label", "flag", "reverse", "Q", "a"}]``, read by the
    columns of the format (a set's first line: nuclides in columns 6-35, label 44-47, flags 48-49, Q 53-64)."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    out, chapter, i = [], None, 0
    while i < len(lines):
        ln = lines[i].ljust(74)
        if ln[:5].strip().isdigit() and not ln[5:].strip():
            chapter, i = int(ln[:5]), i + 3
            continue
        if not ln.strip():
            i += 1
            continue
        if i + 2 >= len(lines):
            break
        l2, l3 = lines[i + 1].ljust(74), lines[i + 2].ljust(74)
        a = [float(l2[13 * k:13 * k + 13]) for k in range(4)] + [float(l3[13 * k:13 * k + 13]) for k in range(3)]
        out.append({"chapter": chapter, "species": [ln[5 + 5 * k:10 + 5 * k].strip() for k in range(6)
                                                    if ln[5 + 5 * k:10 + 5 * k].strip()],
                    "label": ln[43:47], "flag": ln[47].strip(), "reverse": ln[48].strip(), "Q": float(ln[52:64]),
                    "a": np.array(a)})
        i += 3
    return out
