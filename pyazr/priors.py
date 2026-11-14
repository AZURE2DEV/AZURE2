"""Prior parameter covariance files: a multivariate Gaussian prior on fit parameters, its reader, writer and builder.

The ``<config>`` key ``priorCovariance=<path>`` (path relative to the ``.azr``'s directory) makes AZURE2 add

    (p - mu)^T C^-1 (p - mu)

to chi-squared over the parameters the file lists (``Total-Prior-Covariance-Chi-Squared`` in ``chiSquared.out``;
``evaluate(x).terms["prior covariance"]`` in a pyazr session).  It is SAMMY's "PUP" / sequential Bayesian update:
fit A's parameters and covariance become the prior of fit B, which then adds new data.

The file is plain text with ``#`` comments::

    #parametrization 1                      optional: param.sav's basis tag (0 standard, 1 Brune, 2 Park)
    energy_2 width_2_1 segment_1_norm       the parameter NAMES, as in param.sav
    0.8 0.73 1.0                            their MEANS
    1e-6 2e-7 0                             the n rows of the symmetric n x n COVARIANCE, in the same order
    2e-7 4e-4 0
    0 0 1e-4

The values are in the parameter vector's own space, the one ``param.sav`` holds: level energies (MeV), reduced-width
amplitudes of the run's parametrization, normalizations, shifts, ...  A file whose tag differs from the run's
parametrization is refused (the amplitudes would be in another basis; nothing is converted).

The engine's reader (``src/PriorCovariance.cpp``, ``PriorCovariance::Read``) defines the grammar; :func:`read_prior`
reads as it does and refuses what it refuses: a ``#parametrization`` tag counts only before the names line and must be
0, 1 or 2 (a tag without a value is refused); a number is a whole C ``strtod`` token and finite (``nan``, ``inf``,
``1_0``, ``1d0`` are refused, a hexadecimal ``0x1p-3`` is read); tokens are separated by ASCII white space only; a
name listed twice, a matrix not symmetric to 1e-10 of ``sqrt(C_ii C_jj)``, or one whose Cholesky factorization meets
a pivot not above 1e-13 of its diagonal entry is refused; the matrix returned is the symmetrized one.  What only a run
can check (a name the model does not have, a tag of another basis) :func:`build` checks against a session's parameter
table, and the engine when the run starts.

- :func:`read_prior`, :func:`write_prior`: the file;  :func:`read_covariance_dat`: a fit's ``covariance.dat``;
  :func:`prior_from_fit`: the file from a fit's ``param.sav`` and ``covariance.dat`` without a session (the numbers of
  the builder's part ``fit``); :func:`read_chain`: a sampler job's chain.
- :func:`build`: the prior made from named parts over a session's parameter table (:data:`PART_KINDS`: ``fit``,
  ``width``, ``correlation``, ``chain``), checked, with notes on what the engine will do with it: a :class:`Prior`,
  which prints itself, writes the file (:meth:`Prior.write`) and names it in a project (:meth:`Prior.use`).  The parts
  are kept on the file's ``#pyazr.prior`` comment line (:func:`read_parts`), so a builder reopens with them.  A new
  kind of part is one row of :data:`PART_KINDS` and one function in :data:`_BLOCKS` that returns its block.
"""

import json
import os
import re
from dataclasses import dataclass, field

import numpy as np

from .errors import FormatError, ParameterError, ParameterKeyError, close_names, listed

__all__ = ["prior_from_fit", "read_prior", "write_prior", "read_covariance_dat", "read_chain", "build", "Prior",
           "PART_KINDS", "kinds_table", "read_parts", "MARKER", "BASES"]

#: The comment line that carries the parts in a file :meth:`Prior.write` wrote.
MARKER = "#pyazr.prior"

#: The values of the ``#parametrization`` tag (param.sav's) and their names.
BASES = {0: "0 (standard)", 1: "1 (Brune)", 2: "2 (Park)"}


# -- the grammar: the engine's reading ------------------------------------------------------------------------------

_DECIMAL = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?\Z")
_HEX = re.compile(r"[+-]?0[xX]([0-9a-fA-F]+\.?[0-9a-fA-F]*|\.[0-9a-fA-F]+)([pP][+-]?\d+)?\Z")
_SPACE = re.compile(r"[ \t\n\v\f\r]+")   # isspace() of the "C" locale: what std::istringstream splits on


def _number(token):
    """The value of a whole C ``strtod`` token when it is finite, else None (the engine's ``Number``)."""
    if _DECIMAL.match(token):
        v = float(token)
    elif _HEX.match(token):
        v = (-1.0 if token[0] == "-" else 1.0) * float.fromhex(token.lstrip("+-"))
    else:
        return None
    return v if np.isfinite(v) else None


def _tokens(text):
    return [t for t in _SPACE.split(text) if t]


def _cholesky_pivot(cov):
    """The row at which the engine's Cholesky factorization fails (a pivot not above 1e-13 of its diagonal entry),
    or None when it succeeds: ``Cholesky`` of ``src/PriorCovariance.cpp``, its sums in the same order."""
    A = np.array(cov, float)
    n = A.shape[0]
    for j in range(n):
        d = A[j, j]
        diag = d
        for k in range(j):
            d -= A[j, k] * A[j, k]
        if not d > 0.0 or d <= 1.0e-13 * abs(diag):
            return j
        ljj = np.sqrt(d)
        A[j, j] = ljj
        for i in range(j + 1, n):
            s = A[i, j]
            for k in range(j):
                s -= A[i, k] * A[j, k]
            A[i, j] = s / ljj
    return None


def _asymmetry(cov):
    """``(i, j)`` of the first pair the engine finds not symmetric (|a - b| above 1e-10 of max(sqrt|C_ii C_jj|, |a|,
    |b|)), or None."""
    n = cov.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = cov[i, j], cov[j, i]
            scale = max(np.sqrt(abs(cov[i, i] * cov[j, j])), abs(a), abs(b))
            if abs(a - b) > 1.0e-10 * scale:
                return i, j
    return None


def _matrix_refusal(names, cov, where):
    """The engine's refusal of a matrix as its three parts (what, allowed, fix), or None."""
    pair = _asymmetry(cov)
    if pair is not None:
        i, j = pair
        return {"what": f"{where}: the covariance is not symmetric: C({names[i]}, {names[j]}) = {cov[i, j]:.10g} but "
                        f"C({names[j]}, {names[i]}) = {cov[j, i]:.10g}",
                "allowed": "a symmetric matrix (to 1e-10 of sqrt(C_ii C_jj))",
                "fix": "write both halves alike (a matrix built by pyazr.priors.build is symmetric)"}
    bad = _cholesky_pivot(cov)
    if bad is not None:
        return {"what": f"{where}: the covariance is not positive definite (the Cholesky factorization fails at "
                        f"{names[bad]})",
                "allowed": "a positive definite matrix (every Cholesky pivot above 1e-13 of its diagonal entry)",
                "fix": "correlations of magnitude below 1 that fit together, a sigma above 0 on every parameter"}
    return None


def _refuse(what, allowed, fix, error=ParameterError):
    """An error of class ``error`` built from its three parts."""
    return error(f"{what}; allowed: {allowed}; {fix}", what=what, allowed=allowed, fix=fix)


def read_prior(path):
    """``{"names", "mean", "cov", "basis"}`` of a prior covariance file, read as the engine reads it (``basis`` None
    without a tag; ``cov`` symmetrized as the engine symmetrizes it).  A file the engine's reader refuses is refused
    (:class:`pyazr.FormatError` with what, allowed and fix); the module's docstring states the grammar."""
    with open(path, "rb") as fh:
        raw = fh.read().decode("utf-8", "surrogateescape")
    basis, rows, line_of = None, [], []
    for number, line in enumerate(raw.split("\n"), 1):   # std::getline: lines end at '\n' only
        s = line.strip(" \t\r\n")
        if not s:
            continue
        if s[0] == "#":
            t = _tokens(s[1:])
            if t and t[0] == "parametrization" and not rows:
                v = _number(t[1]) if len(t) > 1 else None
                if v not in (0.0, 1.0, 2.0):
                    raise _refuse(f"{path} line {number}: the #parametrization tag must be 0, 1 or 2",
                                  "#parametrization 0, 1 or 2 (param.sav's tag: standard, Brune, Park)",
                                  "write the run's tag, or leave the line out (the file is then taken as it is)",
                                  FormatError)
                basis = int(v)
            continue
        t = _tokens(s.split("#", 1)[0])
        if t:
            rows.append(t)
            line_of.append(number)
    if len(rows) < 2:
        raise _refuse(f"{path}: expected a line of parameter names, a line of means and the covariance rows",
                      "names, means and n rows of n numbers", "pyazr.priors.write_prior writes them", FormatError)
    names = rows[0]
    n = len(names)
    twice = [x for i, x in enumerate(names) if x in names[:i]]
    if twice:
        raise _refuse(f"{path}: parameter {twice[0]} is listed twice", "every name once", "remove one of them",
                      FormatError)
    if len(rows[1]) != n:
        raise _refuse(f"{path}: size mismatch, {n} names but {len(rows[1])} means (line {line_of[1]})",
                      "one mean per name", "write one mean per name, in the names' order", FormatError)
    if len(rows) - 2 != n:
        raise _refuse(f"{path}: size mismatch, {n} parameters but {len(rows) - 2} covariance rows",
                      "n rows of the n x n covariance", "write one row per name", FormatError)

    def numbers(tokens, line):
        out = []
        for tok in tokens:
            v = _number(tok)
            if v is None:
                raise _refuse(f"{path} line {line}: '{tok}' is not a number",
                              "a finite number that C's strtod reads whole (1, -2.5, 1e-3, 0x1p-3)",
                              "write it in decimal (no nan, inf, '_' or Fortran 'd')", FormatError)
            out.append(v)
        return out

    mean = np.array(numbers(rows[1], line_of[1]))
    cov = np.zeros((n, n))
    for i in range(n):
        row = rows[i + 2]
        if len(row) != n:
            raise _refuse(f"{path} line {line_of[i + 2]}: size mismatch, covariance row {i + 1} has {len(row)} "
                          f"entries, expected {n}", "n entries per row", "write one entry per name", FormatError)
        cov[i] = numbers(row, line_of[i + 2])
    refusal = _matrix_refusal(names, cov, path)
    if refusal is not None:
        raise _refuse(refusal["what"], refusal["allowed"], refusal["fix"], FormatError)
    return {"names": names, "mean": mean, "cov": 0.5 * (cov + cov.T), "basis": basis}


def _file_text(names, mean, cov, basis=None, comment=None, parts=None):
    lines = [] if basis is None else [f"#parametrization {int(basis)}"]
    lines += [f"# {c}" for c in (comment.splitlines() if comment else [])]
    if parts is not None:
        lines.append(MARKER + " " + json.dumps({"parts": parts}, separators=(",", ":")))
    lines.append(" ".join(names))
    lines.append(" ".join(f"{float(v):.17g}" for v in mean))
    lines += [" ".join(f"{float(v):.17g}" for v in row) for row in cov]
    return "\n".join(lines) + "\n"


def write_prior(path, names, mean, cov, basis=None, comment=None):
    """Write a prior covariance file (names, means, covariance; optional ``#parametrization`` tag).

    Refused (:class:`pyazr.ParameterError`) where the engine's reader would refuse the file: sizes, a name twice or
    holding white space or ``#``, a number that is not finite, a matrix not symmetric or not positive definite."""
    names = [str(n) for n in names]
    mean = np.asarray(mean, float).ravel()
    cov = np.asarray(cov, float)
    n = len(names)
    if mean.size != n or cov.shape != (n, n):
        raise _refuse(f"size mismatch: {n} names, {mean.size} means, covariance {cov.shape}",
                      "one mean per name and an n x n covariance", "give them in the names' order")
    if len(set(names)) != n:
        raise _refuse("a parameter is listed twice", "every name once", "remove one of them")
    odd = [x for x in names if not x or "#" in x or _SPACE.search(x)]
    if odd:
        raise _refuse(f"names with white space or '#': {listed(odd)}", "the names of param.sav",
                      "Session.parameters.names")
    if basis is not None and int(basis) not in BASES:
        raise _refuse(f"basis {basis}", "0, 1 or 2", "Session.basis")
    if not (np.all(np.isfinite(mean)) and np.all(np.isfinite(cov))):
        raise _refuse("a mean or covariance entry is not finite", "finite numbers", "check where they come from")
    refusal = _matrix_refusal(names, cov, "write_prior")
    if refusal is not None:
        raise _refuse(refusal["what"], refusal["allowed"], refusal["fix"])
    with open(path, "w") as fh:
        fh.write(_file_text(names, mean, cov, basis, comment))
    return path


# -- a fit's files and a chain --------------------------------------------------------------------------------------

def _read_sav(path):
    """``({name: value}, basis or None)`` from a param.sav / param.par / param.fit."""
    vals, basis = {}, None
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if len(t) < 2:
                continue
            if t[0] in ("#parametrization", "parametrization"):
                try:
                    basis = int(round(float(t[1])))
                except ValueError:
                    pass
                continue
            if t[0].startswith("#"):
                continue
            try:
                vals[t[0]] = float(t[1])
            except ValueError:
                continue
    return vals, basis


def read_covariance_dat(path):
    """``(C, columns, scale)`` from an AZURE2 ``output/covariance.dat``.

    The file holds the covariance of the fit's FREE level energies and reduced-width amplitudes (no normalization,
    shift or other nuisance parameter: the cross-section band does not depend on them), in parameter-vector order --
    the order of ``param.sav`` with the fixed parameters and every non-R-matrix parameter left out -- in the
    param.sav space (Minuit's external parameters).  A fit run with ``--scale-covariance`` (or the "scale to reduced
    chi-squared = 1" answer) multiplied it by the reduced chi-squared when that was above 1.

    Builds since 2026-11-05 head it with ``#columns <names>`` and, when scaled, ``#scaled <factor>``; for an older
    file ``columns`` is None and ``scale`` 1.0 (unknown).
    """
    columns, scale = None, 1.0
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if not t or not t[0].startswith("#"):
                continue
            if t[0] == "#columns":
                columns = t[1:]
            elif t[0] == "#scaled" and len(t) > 1:
                scale = float(t[1])
    C = np.loadtxt(path, ndmin=2)
    if C.ndim != 2 or C.shape[0] != C.shape[1]:
        raise ValueError(f"{path}: expected a square matrix, got shape {C.shape}")
    if columns is not None and len(columns) != C.shape[0]:
        raise ValueError(f"{path}: {len(columns)} column names for a {C.shape[0]}x{C.shape[0]} matrix")
    return C, columns, scale


def _fit_block(sav_path, covariance_path, names=None, columns=None):
    """``(names, mean, cov, basis, comment, scale)`` of an earlier fit's ``param.sav`` and ``covariance.dat``: the
    means of ``names`` (every column when None) and their marginal block.  The one reading of a fit's files, for
    :func:`prior_from_fit` and the builder's part ``fit``."""
    vals, basis = _read_sav(sav_path)
    C, cols, scale = read_covariance_dat(covariance_path)
    if cols is None:
        if not columns:
            raise _refuse(f"{covariance_path} has no '#columns' header (written before 2026-11-05)",
                          f"the {C.shape[0]} column names given (columns=[...])",
                          "pass columns=[...]: the free energy_*/width_* names of that fit in param.sav order")
        cols = list(columns)
        if len(cols) != C.shape[0]:
            raise _refuse(f"{len(cols)} column names for a {C.shape[0]}x{C.shape[0]} covariance",
                          f"{C.shape[0]} names", "one name per column of covariance.dat")
    missing = [c for c in cols if c not in vals]
    if missing:
        raise _refuse(f"{sav_path} does not hold {listed(missing)} (the covariance belongs to another fit?)",
                      "the param.sav written by the fit that wrote covariance.dat",
                      "take both files from one fit's output directory")
    if not names:
        idx = list(range(len(cols)))
    else:
        pos = {c: i for i, c in enumerate(cols)}
        bad = [n for n in names if n not in pos]
        if bad:
            raise _refuse(f"{listed(bad)} not in {covariance_path}",
                          f"its columns, the fit's free level energies and widths: {listed(cols)}",
                          "choose among the columns, or give the others a part of their own (width)")
        idx = [pos[n] for n in names]
    sel = [cols[i] for i in idx]
    mean = np.array([vals[c] for c in sel])
    cov = C[np.ix_(idx, idx)]
    comment = (f"prior from {os.path.abspath(sav_path)}\n"
               f"and {os.path.abspath(covariance_path)}" + (f" (scaled by {scale:g})" if scale != 1.0 else ""))
    return sel, mean, cov, basis, comment, scale


def prior_from_fit(sav_path, covariance_path, out_path, names=None, columns=None):
    """Write a prior covariance file from a previous fit: means from its ``param.sav``, covariance from its
    ``covariance.dat``.  No session is needed; with one, :func:`build` with a part ``{"kind": "fit", ...}`` makes the
    same names, means and covariance and checks them against the run.

    ``covariance.dat`` (written by a fit with the uncertainty band on: ``--covariance-band`` or the band prompt's
    "y") covers the fit's FREE level energies and reduced-width amplitudes only, in param.sav order and space (see
    :func:`read_covariance_dat`); fixed parameters and normalizations, shifts and the other nuisance parameters are
    not in it, so they are not in the prior either.  A covariance that was scaled to reduced chi-squared = 1 is used
    as it is (the ``#scaled`` factor is copied into a comment).

    ``names``: optional subset (in any order) of the covariance's parameters; the prior is the corresponding
    marginal (the sub-block), which is what a Gaussian prior on fewer parameters is.  ``columns``: the column names
    of a ``covariance.dat`` written before 2026-11-05 (no ``#columns`` header) -- the free ``energy_*`` /
    ``width_*`` names of that fit in param.sav order; without the header and without ``columns`` the file is refused,
    since param.sav does not say which parameters were free.

    The file carries param.sav's ``#parametrization`` tag, so a run in another parametrization refuses it.  Returns
    ``out_path``.
    """
    sel, mean, cov, basis, comment, _ = _fit_block(sav_path, covariance_path, names, columns)
    return write_prior(out_path, sel, mean, cov, basis=basis, comment=comment)


def read_chain(path, discard=None):
    """``(names, draws, discard)`` of a sampler's chain as pyazr writes it: ``chain_<id>.npy`` (stored steps, walkers,
    parameters), written by the sampler job of the web page (``pyazr/ui/jobs.py``), with ``chain_<id>.json`` beside it
    (its ``names``, the free parameters in the order of x, and its burn-in ``discard``).  ``draws``: the stored steps
    after ``discard`` (the JSON's when None), every walker, as (draws, parameters)."""
    stem = os.path.splitext(path)[0]
    try:
        chain = np.load(path)
    except (OSError, ValueError) as e:
        raise _refuse(f"cannot read the chain {path} ({e})", "a chain_<id>.npy a sampler job wrote",
                      "choose it in the project's output directory") from None
    try:
        with open(stem + ".json") as fh:
            meta = json.load(fh)
    except (OSError, ValueError):
        raise _refuse(f"{os.path.basename(stem)}.json, beside the chain, is missing or unreadable (it names the "
                      f"chain's parameters)", "a chain with the summary its sampler job wrote beside it",
                      "keep chain_<id>.json beside chain_<id>.npy") from None
    names = [str(n) for n in meta.get("names") or []]
    if chain.ndim != 3 or chain.shape[2] != len(names):
        raise _refuse(f"{path}: a chain of shape {chain.shape} for {len(names)} names",
                      f"(steps, walkers, {len(names)})", "the chain and the summary of one sampler job")
    if discard is None or discard == "":
        discard = meta.get("discard", chain.shape[0] // 2)
    discard = int(discard)
    if not 0 <= discard < chain.shape[0]:
        raise _refuse(f"a burn-in of {discard} for a chain of {chain.shape[0]} stored steps",
                      f"0 to {chain.shape[0] - 1}", "a burn-in shorter than the chain")
    return names, chain[discard:].reshape(-1, chain.shape[2]), discard


# -- the builder ----------------------------------------------------------------------------------------------------

#: The kinds of part: their fields (name, type, default, unit, doc) and what they are.  A field with default None is
#: required.  Types: ``path:<suffixes>`` (a file relative to the project's directory; the suffixes a page offers),
#: ``names`` (parameter names), ``number`` (unit ``unit``: the parameter's own), ``integer``, ``pairs``
#: (``[name, name, rho]`` rows), ``matrix`` (rows of numbers).
PART_KINDS = {
    "fit": {"doc": "an earlier fit: the means of its param.sav and the covariance of its covariance.dat (its free "
                   "level energies and widths), for the parameters chosen",
            "fields": (("sav", "path:.sav,.par,.fit", None, "", "the fit's param.sav"),
                       ("covariance", "path:.dat", None, "", "the fit's covariance.dat (a fit with --covariance-band)"),
                       ("names", "names", "", "", "the parameters taken (empty: every column of covariance.dat)"),
                       ("columns", "names", "", "",
                        "only for a covariance.dat without its #columns line: its column names in order"))},
    "width": {"doc": "a standard deviation per parameter, uncorrelated, around the parameter's start value or a "
                     "given centre",
              "fields": (("names", "names", None, "", "the parameters"),
                         ("sigma", "number", "", "unit", "the standard deviation in the parameter's unit"),
                         ("percent", "number", "", "%", "or the standard deviation in percent of the centre"),
                         ("centre", "number", "", "unit",
                          "the centre of every parameter named (empty: each parameter's start value)"))},
    "correlation": {"doc": "correlation coefficients between parameters of the other parts (their sigmas stay)",
                    "fields": (("pairs", "pairs", "", "", "[name, name, rho] rows, -1 < rho < 1"),
                               ("names", "names", "", "", "with matrix: the parameters of its rows"),
                               ("matrix", "matrix", "", "", "a correlation matrix over names (diagonal 1)"))},
    "chain": {"doc": "the sample mean and covariance of a sampler's chain after its burn-in, for the parameters "
                     "chosen",
              "fields": (("file", "path:.npy", None, "", "chain_<id>.npy of a sampler job, chain_<id>.json beside it"),
                         ("names", "names", "", "", "the parameters taken (empty: every parameter of the chain)"),
                         ("discard", "integer", "", "",
                          "stored steps dropped at the start (empty: the chain's own burn-in)"))},
}


def kinds_table():
    """:data:`PART_KINDS` as plain data (a page's forms): ``{kind: {"doc", "fields": [{name, type, default, unit,
    doc}]}}``."""
    return {k: {"doc": v["doc"], "fields": [dict(zip(("name", "type", "default", "unit", "doc"), f))
                                            for f in v["fields"]]} for k, v in PART_KINDS.items()}


def _normal_part(part):
    """The part with every field of its kind (defaults filled, an absent optional field ""), numbers as floats,
    names as lists."""
    if not isinstance(part, dict) or "kind" not in part:
        raise _refuse(f"a part {part!r}", f"a dict with a 'kind': {', '.join(PART_KINDS)}",
                      "pyazr.priors.PART_KINDS lists the kinds and their fields")
    kind = part["kind"]
    if kind not in PART_KINDS:
        raise _refuse(f"no kind of part {kind!r}{close_names(kind, PART_KINDS)}", ", ".join(PART_KINDS),
                      "pyazr.priors.PART_KINDS")
    fields = PART_KINDS[kind]["fields"]
    known = [f[0] for f in fields]
    extra = sorted(set(part) - set(known) - {"kind"})
    if extra:
        raise _refuse(f"part {kind!r} has no field {', '.join(extra)}", f"its fields: {', '.join(known)}",
                      "pyazr.priors.PART_KINDS")
    out = {"kind": kind}
    for name, typ, default, unit, doc in fields:
        value = part.get(name, default)
        if value is None or (isinstance(value, (str, list, tuple)) and len(value) == 0):
            if default is None:
                raise _refuse(f"part {kind!r} needs {name}", doc, f"give {name}")
            out[name] = ""
            continue
        if typ == "names":
            value = [value] if isinstance(value, str) else [str(v) for v in value]
        elif typ in ("number", "integer"):
            try:
                value = float(value)
            except (TypeError, ValueError):
                raise _refuse(f"part {kind!r}: {name} {value!r}", "a number" + (f" ({unit})" if unit else ""),
                              doc) from None
            if not np.isfinite(value):
                raise _refuse(f"part {kind!r}: {name} {value}", "a finite number", doc)
            if typ == "integer":
                if value != int(value):
                    raise _refuse(f"part {kind!r}: {name} {value:g}", "a whole number", doc)
                value = int(value)
        elif typ == "pairs":
            rows = []
            for row in value:
                if not isinstance(row, (list, tuple)) or len(row) != 3:
                    raise _refuse(f"part {kind!r}: a pair {row!r}", "[name, name, rho]", doc)
                try:
                    rho = float(row[2])
                except (TypeError, ValueError):
                    raise _refuse(f"part {kind!r}: rho {row[2]!r}", "a number between -1 and 1", doc) from None
                rows.append([str(row[0]), str(row[1]), rho])
            value = rows
        elif typ == "matrix":
            try:
                value = np.asarray(value, float).tolist()
            except (TypeError, ValueError):
                raise _refuse(f"part {kind!r}: matrix {value!r}", "rows of numbers", doc) from None
        else:
            value = str(value)
        out[name] = value
    return out


@dataclass
class _Block:
    """What a part gives: its parameters, their centres and covariance, its notes."""
    names: list
    mean: np.ndarray
    cov: np.ndarray
    notes: list = field(default_factory=list)


def _path(directory, p):
    return p if os.path.isabs(p) else os.path.join(directory, p)


def _fit_part(part, index, table, directory, basis):
    names, mean, cov, tag, _, scale = _fit_block(_path(directory, part["sav"]), _path(directory, part["covariance"]),
                                                 part["names"] or None, part["columns"] or None)
    if tag is not None and basis is not None and tag != basis:
        raise _refuse(f"part {index + 1} (fit): {part['sav']} holds parametrization {BASES.get(tag, tag)} amplitudes "
                      f"and this run uses parametrization {BASES.get(basis, basis)}",
                      "a fit in the run's parametrization (the engine converts nothing)",
                      "rerun the fit the prior comes from in this run's parametrization, or run this one in the fit's")
    notes = []
    if tag is None:
        notes.append(f"part {index + 1} (fit): {part['sav']} has no #parametrization tag: its amplitudes are taken "
                     f"as this run's, basis {BASES.get(basis, basis)}.")
    if scale != 1.0:
        notes.append(f"part {index + 1} (fit): its covariance.dat was scaled by {scale:g} (to a reduced chi-squared "
                     f"of 1); it is used as it is.")
    return _Block(names, mean, cov, notes)


def _width_part(part, index, table, directory, basis):
    has_sigma, has_percent = part["sigma"] != "", part["percent"] != ""
    if has_sigma == has_percent:
        raise _refuse(f"part {index + 1} (width): " + ("both sigma and percent" if has_sigma else
                                                       "neither sigma nor percent"),
                      "one of them: sigma in the parameter's unit, or percent of the centre", "give sigma or percent")
    centres = np.array([part["centre"] if part["centre"] != "" else table[n].start for n in part["names"]], float)
    sigma = np.full(len(centres), part["sigma"]) if has_sigma else part["percent"] / 100.0 * np.abs(centres)
    bad = [n for n, s in zip(part["names"], sigma) if not s > 0]
    if bad:
        raise _refuse(f"part {index + 1} (width): a standard deviation of 0 or less on {listed(bad)}"
                      + (" (percent of a centre 0)" if has_percent else ""), "a sigma above 0",
                      "give sigma in the parameter's unit" if has_percent else "a positive sigma")
    return _Block(list(part["names"]), centres, np.diag(sigma ** 2))


def _chain_part(part, index, table, directory, basis):
    names, draws, discard = read_chain(_path(directory, part["file"]), part["discard"])
    chosen = part["names"] or names
    pos = {n: i for i, n in enumerate(names)}
    bad = [n for n in chosen if n not in pos]
    if bad:
        raise _refuse(f"part {index + 1} (chain): {listed(bad)} not in the chain", f"its parameters: {listed(names)}",
                      "choose among them, or give the others a part of their own (width)")
    idx = [pos[n] for n in chosen]
    if draws.shape[0] <= len(idx):
        raise _refuse(f"part {index + 1} (chain): {draws.shape[0]} draws after a burn-in of {discard} for "
                      f"{len(idx)} parameters", "more draws than parameters", "a longer chain or a shorter burn-in")
    sub = draws[:, idx]
    notes = [f"part {index + 1} (chain): the mean and covariance of {sub.shape[0]} draws after a burn-in of "
             f"{discard} stored steps; a chain does not record its basis: it is taken as this run's, "
             f"{BASES.get(basis, basis)}."]
    return _Block(list(chosen), sub.mean(axis=0), np.atleast_2d(np.cov(sub, rowvar=False)), notes)


#: The function that makes the block of each kind of part (a part ``correlation`` acts on the others' blocks).
_BLOCKS = {"fit": _fit_part, "width": _width_part, "chain": _chain_part}


def _correlations(part, index, owner):
    """``{(a, b): rho}`` of a part ``correlation``; every name a parameter of another part."""
    pairs = list(part["pairs"] or [])
    if part["matrix"] != "" or part["names"] != "":
        names, m = part["names"] or [], np.asarray(part["matrix"] or [], float)
        if m.shape != (len(names), len(names)):
            raise _refuse(f"part {index + 1} (correlation): a matrix of shape {m.shape} over {len(names)} names",
                          "an n x n matrix over the n names", "give names and matrix together")
        if np.max(np.abs(np.diag(m) - 1.0)) > 1e-12 or np.max(np.abs(m - m.T)) > 1e-12:
            raise _refuse(f"part {index + 1} (correlation): not a correlation matrix", "symmetric with diagonal 1",
                          "write rho_ij = rho_ji and 1 on the diagonal")
        pairs += [[names[i], names[j], m[i, j]] for i in range(len(names)) for j in range(i + 1, len(names))]
    if not pairs:
        raise _refuse(f"part {index + 1} (correlation): no pair", "pairs [name, name, rho], or names and a matrix",
                      "give at least one pair")
    out = {}
    for a, b, rho in pairs:
        for n in (a, b):
            if n not in owner:
                raise _refuse(f"part {index + 1} (correlation): {n} is in no other part",
                              f"parameters of the other parts: {listed(list(owner))}",
                              "give it a part of its own first (width), then correlate it")
        if a == b:
            raise _refuse(f"part {index + 1} (correlation): {a} with itself", "two different parameters",
                          "choose two parameters")
        if not -1.0 < rho < 1.0:
            raise _refuse(f"part {index + 1} (correlation): rho {rho:g} between {a} and {b}", "-1 < rho < 1",
                          "a coefficient of magnitude below 1")
        if owner[a] == owner[b] and owner[a][1] != "width":
            raise _refuse(f"part {index + 1} (correlation): {a} and {b} are already correlated by part "
                          f"{owner[a][0] + 1} ({owner[a][1]})", "pairs of different parts, or of one part width",
                          "leave that pair to its part")
        key = tuple(sorted((a, b)))
        if key in out:
            raise _refuse(f"part {index + 1} (correlation): the pair {a}, {b} twice", "each pair once",
                          "remove one of them")
        out[key] = rho
    return out


@dataclass
class Prior:
    """A prior on fit parameters built from parts (:func:`build`); ``print(prior)`` shows it.

    - ``names``, ``mean`` (the centres), ``cov`` (the covariance), in the order of the parts; ``sigma`` and
      ``correlation`` from them; ``units``, ``kinds``, ``fixed`` and ``bounds`` per name from the parameter table;
      ``basis`` the run's (the file's ``#parametrization`` tag); ``parts`` normalized, ``part_of`` each name's part
      (0-based).
    - ``checks``: ``symmetric``, ``positive_definite`` (the engine's Cholesky), ``smallest_eigenvalue``,
      ``largest_eigenvalue``, ``condition_number``; ``refusal``: why the engine would refuse the file (``what``,
      ``allowed``, ``fix``; None when it would not), ``refusal_parts`` the parts concerned (0-based).
    - ``notes``: what the engine does with it that the user should know (sentences).
    """
    names: list
    mean: np.ndarray
    cov: np.ndarray
    basis: int
    parts: list
    part_of: list
    units: list
    kinds: list
    fixed: list
    bounds: list
    checks: dict = field(default_factory=dict)
    refusal: dict = None
    refusal_parts: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    @property
    def sigma(self):
        return np.sqrt(np.abs(np.diag(self.cov)))

    @property
    def correlation(self):
        s = self.sigma
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.cov / np.outer(s, s)

    def __str__(self):
        lines = [f"prior on {len(self.names)} parameter(s), parametrization {BASES.get(self.basis, self.basis)}, "
                 f"from {len(self.parts)} part(s)"]
        rows = [("name", "part", "unit", "centre", "sigma", "")]
        for i, n in enumerate(self.names):
            rows.append((n, f"{self.part_of[i] + 1} {self.parts[self.part_of[i]]['kind']}", self.units[i],
                         f"{self.mean[i]:.8g}", f"{self.sigma[i]:.4g}", "fixed" if self.fixed[i] else ""))
        widths = [max(len(r[c]) for r in rows) for c in range(len(rows[0]))]
        lines += ["  " + "  ".join(c.ljust(widths[k]) for k, c in enumerate(r)).rstrip() for r in rows]
        if len(self.names) > 1:
            lines.append("correlation:")
            w = max(len(n) for n in self.names)
            lines += ["  " + n.ljust(w) + " " + " ".join(f"{v:6.3f}" for v in self.correlation[i])
                      for i, n in enumerate(self.names)]
        c = self.checks
        if self.refusal:
            lines.append(f"refused: {self.refusal['what']}; allowed: {self.refusal['allowed']}; {self.refusal['fix']}")
        else:
            lines.append(f"positive definite: eigenvalues {c['smallest_eigenvalue']:.4g} to "
                         f"{c['largest_eigenvalue']:.4g}, condition number {c['condition_number']:.4g}")
        lines += [f"note: {n}" for n in self.notes]
        return "\n".join(lines)

    def text(self, comment=None):
        """The file's text: the run's ``#parametrization`` tag, a comment, the parts' line (``#pyazr.prior``), the
        names, the centres and the covariance (17 significant digits)."""
        words = f"prior covariance built from {len(self.parts)} part(s) by pyazr.priors.build"
        return _file_text(self.names, self.mean, self.cov, self.basis, words + ("\n" + comment if comment else ""),
                          self.parts)

    def write(self, path, comment=None):
        """Write the file the engine reads (refused while :attr:`refusal` is set); returns ``path``."""
        if self.refusal:
            raise _refuse(f"the prior is not written: {self.refusal['what']}", self.refusal["allowed"],
                          self.refusal["fix"])
        with open(path, "w") as fh:
            fh.write(self.text(comment))
        return path

    def use(self, project, path):
        """Write the file at ``path`` (relative to the directory of the :class:`pyazr.Project` ``project``) and set
        the project's ``<config>`` key ``priorCovariance`` to it (``Project.set_config``; saving the project is the
        caller's).  Returns ``path``."""
        if os.path.isabs(str(path)):
            raise _refuse(f"{path}: an absolute path", "a path relative to the project's directory",
                          "give the file's name inside the project's directory")
        self.write(os.path.join(project.directory, path))
        project.set_config(priorCovariance=str(path))
        return path

    def summary(self):
        """Plain data for a caller that draws it (the web page): the record's fields, the checks and the file's
        text (None while refused)."""
        return {"names": list(self.names), "mean": self.mean.tolist(), "sigma": self.sigma.tolist(),
                "correlation": self.correlation.tolist(), "covariance": self.cov.tolist(), "basis": self.basis,
                "basis_name": BASES.get(self.basis, str(self.basis)), "parts": self.parts, "part_of": self.part_of,
                "units": self.units, "kinds": self.kinds, "fixed": self.fixed,
                "bounds": [[float(v) if np.isfinite(v) else None for v in b] for b in self.bounds],
                "checks": self.checks, "refusal": self.refusal, "refusal_parts": self.refusal_parts,
                "notes": self.notes, "text": None if self.refusal else self.text()}


def build(parameters, parts, *, directory="."):
    """The prior made from ``parts`` (dicts, :data:`PART_KINDS`) over the parameter table ``parameters`` (a
    session's :attr:`pyazr.Session.parameters`: names, kinds, units, fixed marks, limits, the priors of their own and
    the basis), the files the parts name read relative to ``directory`` (the project's).  Returns a :class:`Prior`.

    Refused (:class:`pyazr.ParameterError`, with what, allowed and fix): a part that cannot be made, a name the table
    does not have (:class:`pyazr.ParameterKeyError`, the close names said), a parameter in two parts, a fit of another
    basis.  A matrix the engine would refuse comes back with ``refusal`` set.  Notes: a parameter the run fixes (the
    engine marginalizes it), a parameter with a prior of its own (the engine replaces that prior by this one), a
    centre outside the parameter's limits, a fit or a chain whose basis is not recorded."""
    if not parts:
        raise _refuse("no part", "at least one part fit, width or chain", "add a part")
    parts = [_normal_part(p) for p in parts]
    table = {p.name: p for p in parameters}
    basis = getattr(parameters, "basis", None)

    def known(index, kind, names):
        for n in names:
            if n not in table:
                raise _refuse(f"part {index + 1} ({kind}): {n} is not a parameter of this project"
                              f"{close_names(n, table)}", "the names of the session's parameter table",
                              "print(session.parameters.table()) lists them", ParameterKeyError)

    blocks, owner, correlation_parts = [], {}, []
    for i, part in enumerate(parts):
        if part["kind"] == "correlation":
            correlation_parts.append(i)
            continue
        if part["kind"] == "width":
            known(i, "width", part["names"])
        block = _BLOCKS[part["kind"]](part, i, table, directory, basis)
        known(i, part["kind"], block.names)
        for n in block.names:
            if n in owner:
                raise _refuse(f"{n} is in part {owner[n][0] + 1} ({owner[n][1]}) and in part {i + 1} "
                              f"({part['kind']})", "one part per parameter", "remove it from one of them")
            owner[n] = (i, part["kind"])
        blocks.append(block)
    if not blocks:
        raise _refuse("only correlation parts", "a part fit, width or chain that gives the parameters",
                      "add their part first")
    names = [n for b in blocks for n in b.names]
    mean = np.concatenate([b.mean for b in blocks])
    cov = np.zeros((len(names), len(names)))
    k = 0
    for b in blocks:
        cov[k:k + len(b.names), k:k + len(b.names)] = b.cov
        k += len(b.names)
    pos = {n: i for i, n in enumerate(names)}
    for i in correlation_parts:
        for (a, b), rho in _correlations(parts[i], i, owner).items():
            ia, ib = pos[a], pos[b]
            cov[ia, ib] = cov[ib, ia] = rho * np.sqrt(cov[ia, ia] * cov[ib, ib])
    rows = [table[n] for n in names]
    prior = Prior(names=names, mean=mean, cov=cov, basis=basis, parts=parts, part_of=[owner[n][0] for n in names],
                  units=[p.unit for p in rows], kinds=[p.kind for p in rows], fixed=[bool(p.fixed) for p in rows],
                  bounds=[tuple(float(v) for v in p.bounds) for p in rows])
    for b in blocks:
        prior.notes += b.notes
    _check(prior, correlation_parts)
    _notes(prior, rows)
    return prior


def _check(prior, correlation_parts):
    cov = prior.cov
    eig = np.linalg.eigvalsh(0.5 * (cov + cov.T))
    bad = _cholesky_pivot(cov)
    prior.checks.update(symmetric=_asymmetry(cov) is None, positive_definite=bad is None,
                        smallest_eigenvalue=float(eig[0]), largest_eigenvalue=float(eig[-1]),
                        condition_number=float(eig[-1] / eig[0]) if eig[0] > 0 else float("inf"))
    refusal = _matrix_refusal(prior.names, cov, "the prior")
    if refusal is not None:
        prior.refusal = refusal
        concerned = set(correlation_parts)
        if bad is not None:
            concerned.add(prior.part_of[bad])
        prior.refusal_parts = sorted(concerned)


def _notes(prior, rows):
    notes = prior.notes
    fixed = [p.name for p in rows if p.fixed]
    if fixed:
        one = len(fixed) == 1
        notes.append(f"{listed(fixed)} {'is' if one else 'are'} fixed in this run: the engine marginalizes "
                     f"{'it' if one else 'them'} (row and column dropped, a Note in the run); the prior acts on the "
                     f"other {len(prior.names) - len(fixed)}.")
    for p in rows:
        if (p.prior_sigma or 0) > 0:
            notes.append(f"{p.name} has a prior of its own in the project (centre {p.prior_centre:.6g}, sigma "
                         f"{p.prior_sigma:.4g} {p.unit}): the engine replaces it by this multivariate prior (its "
                         f"one-parameter term is not added; a Note in the run says so).")
    for i, p in enumerate(rows):
        lo, hi = p.bounds
        if not lo <= prior.mean[i] <= hi:
            notes.append(f"{p.name}: the centre {prior.mean[i]:.6g} is outside the parameter's limits "
                         f"[{lo:.6g}, {hi:.6g}] {p.unit}: the fit cannot reach it.")


def read_parts(path):
    """``{"parts": [...]}`` of a file :meth:`Prior.write` wrote (its ``#pyazr.prior`` line); None for a file without
    the line (written by hand, or by :func:`prior_from_fit`)."""
    with open(path) as fh:
        for line in fh:
            if line.startswith(MARKER + " "):
                return json.loads(line[len(MARKER) + 1:])
    return None
