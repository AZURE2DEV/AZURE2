"""The data covariance builder: a data set's covariance file made from named parts.

A data set whose record names ``covariance <file>`` has its points' chi-squared computed as ``v^T C^-1 v``, ``v`` the
model minus the data (``docs/source/user_guide/chi_squared.rst``; the engine's reader is
``experiment/src/CovarianceFile.cpp``, its use ``include/DataCovariance.h``).  What the engine does with the file, and
so what this module writes:

- The file lists the data files it covers (``file <path> <rows>``, every data row of each file, as the engine counts
  them: lines that are neither blank nor comments), then terms that add up to ``C``: ``diagonal`` (N numbers),
  ``vector`` (adds ``v v^T``), ``relvector`` (fractions ``f``: adds ``(f m)(f m)^T`` with ``m`` the model, frozen at the
  run's first evaluation), ``block`` (a symmetric block on chosen rows) and ``matrix`` (N x N).
- The numbers are compared with the rows in the engine's units of the quantity (the row's value and error after the
  data set's column units, ``experiment/src/Measurement.cpp``); this module writes them in those units.
- Only the rows a data set fits (inside its window) enter; the others drop out.  The ``diagonal``, ``block`` and
  ``matrix`` terms on the fitted rows must be positive definite by themselves: every fitted row needs a variance of its
  own (a fully correlated part alone is refused).
- The error column of a fitted row is replaced by ``sqrt(C_ii)`` (the data value standing in for ``m``); the run warns
  when the column differs by more than 1 %.  The chi-squared does not depend on the column.
- A data set's normalization ``n`` scales the data and the covariance together (``M = D^-1 C D^-1`` on residuals
  ``(fit - n d)/(n sigma)``), as it scales the error column.  A fitted normalization with a prior of e % and a fully
  correlated scale part of e % are the same uncertainty counted twice: keep one (:func:`build` says so in its notes).

The parts (one dict each, ``{"kind": ..., <fields>}``; :data:`PART_KINDS` is the table of kinds and their fields):

- ``own``: the points' own, uncorrelated errors (a ``diagonal`` term): the file's error column (``source="column"``),
  a percentage of the value (``source="percent"``, ``percent``), or what the error column holds beyond the other parts
  (``source="remainder"``: the column is the total error; refused where the other parts exceed it).
- ``scale``: a common scale uncertainty of ``percent`` %, fully correlated over every row (``basis="model"``: a
  ``relvector`` term, the model's scale, free of Peelle's pertinent puzzle; ``basis="data"``: a ``vector`` term of the
  data values).
- ``range``: as ``scale``, on the rows with ``e_min <= E_lab <= e_max`` (MeV) only.
- ``length``: ``percent`` % of the value, correlated as ``exp(-|E_i - E_j| / length)`` (``shape="exponential"``) or
  ``exp(-(E_i - E_j)^2 / (2 length^2))`` (``shape="gaussian"``), ``length`` in MeV (lab); optionally on ``e_min`` to
  ``e_max`` only.  A ``block`` term on the data values.
- ``correlation``: a given correlation matrix (``file``, a text file of N x N numbers over every row the covariance
  covers, ``#`` comments) with the error column (``percent`` absent) or ``percent`` % of the value as the standard
  deviation.  A ``block`` term.

:func:`build` returns a :class:`Covariance` (the matrix on the fitted rows, the total error per point beside the
file's, the checks, the notes); :meth:`Covariance.write` writes the file with the parts in a comment line
(``#pyazr.covariance {...}``) that :func:`read_parts` reads back, so a builder reopens with them.  A new kind of part
is one row of :data:`PART_KINDS` and one function that returns its terms.
"""
import json
import os
from dataclasses import dataclass, field

import numpy as np

from .paths import is_absolute

__all__ = ["PART_KINDS", "Covariance", "build", "read_parts", "rows_of", "MARKER"]

#: The comment line that carries the parts in a written file.
MARKER = "#pyazr.covariance"

#: Rows above which the dense checks are not made (the matrix on the fitted rows is formed for them).
DENSE_LIMIT = 8000

#: The kinds of part: their fields (name, type, default, unit, doc) and what they are.  A field with default None is
#: required unless the doc says otherwise.
PART_KINDS = {
    "own": {"doc": "the points' own, uncorrelated errors (diagonal)",
            "fields": (("source", "choice:column,percent,remainder", "column", "",
                        "column: the error column; percent: a percentage of the value; remainder: the column is the "
                        "total error, this part what the other parts leave of it"),
                       ("percent", "number", "", "%", "with source=percent: the error as a percentage of the value"))},
    "scale": {"doc": "a common scale uncertainty, fully correlated over every row",
              "fields": (("percent", "number", None, "%", "the scale's standard deviation"),
                         ("basis", "choice:model,data", "model",  "",
                          "model: relvector (the model's scale, frozen at the first evaluation); data: vector of "
                          "the data values"))},
    "range": {"doc": "a scale uncertainty fully correlated over an energy range",
              "fields": (("percent", "number", None, "%", "the standard deviation in the range"),
                         ("e_min", "number", None, "MeV", "lab energy, lower end (inclusive)"),
                         ("e_max", "number", None, "MeV", "lab energy, upper end (inclusive)"),
                         ("basis", "choice:model,data", "model", "", "as for scale"))},
    "length": {"doc": "a part correlated over a correlation length in energy",
               "fields": (("percent", "number", None, "%", "standard deviation, of the data value"),
                          ("length", "number", None, "MeV", "correlation length (lab energy)"),
                          ("shape", "choice:exponential,gaussian", "exponential", "",
                           "exp(-|dE|/length) or exp(-dE^2/(2 length^2))"),
                          ("e_min", "number", "", "MeV", "optional: only rows from this lab energy"),
                          ("e_max", "number", "", "MeV", "optional: only rows up to this lab energy"))},
    "correlation": {"doc": "a given correlation matrix",
                    "fields": (("file", "path", None, "", "N x N correlation matrix over every row covered "
                                                          "(path relative to the project's directory)"),
                               ("percent", "number", "", "%",
                                "the standard deviation as a percentage of the value; empty: the error column"))},
}


def _field_default(kind, name):
    for f in PART_KINDS[kind]["fields"]:
        if f[0] == name:
            return f[2]
    raise KeyError(name)


def _normal_part(part):
    """The part with every field of its kind (defaults filled), numbers as floats; ValueError naming the problem."""
    if not isinstance(part, dict) or "kind" not in part:
        raise ValueError(f"a part is a dict with a 'kind' ({', '.join(PART_KINDS)}), not {part!r}")
    kind = part["kind"]
    if kind not in PART_KINDS:
        raise ValueError(f"unknown kind of part {kind!r}: one of {', '.join(PART_KINDS)}")
    names = [f[0] for f in PART_KINDS[kind]["fields"]]
    extra = sorted(set(part) - set(names) - {"kind"})
    if extra:
        raise ValueError(f"part {kind!r} has no field {', '.join(extra)} (its fields: {', '.join(names)})")
    out = {"kind": kind}
    for name, typ, default, unit, doc in PART_KINDS[kind]["fields"]:
        value = part.get(name, default)
        if value is None or (value == "" and default is None):
            raise ValueError(f"part {kind!r} needs {name} ({doc})")
        if typ == "number" and value != "":
            try:
                value = float(value)
            except (TypeError, ValueError):
                raise ValueError(f"part {kind!r}: {name} is a number ({unit}), not {value!r}") from None
            if not np.isfinite(value):
                raise ValueError(f"part {kind!r}: {name} must be finite")
        if typ.startswith("choice:") and value not in typ[7:].split(","):
            raise ValueError(f"part {kind!r}: {name} is one of {typ[7:]}, not {value!r}")
        out[name] = value
    for name in ("percent", "length"):
        if isinstance(out.get(name), float) and out[name] <= 0:
            raise ValueError(f"part {kind!r}: {name} must be positive")
    if kind == "own" and out["source"] == "percent" and out["percent"] == "":
        raise ValueError("part 'own' with source percent needs percent (the error as a percentage of the value)")
    if kind == "range" and not out["e_min"] <= out["e_max"]:
        raise ValueError("part 'range': e_min must not exceed e_max")
    return out


# -- the rows ------------------------------------------------------------------------------------------------------

@dataclass
class Rows:
    """Every row a covariance covers, file after file (the order of its ``file`` lines), as the engine reads them."""
    files: list           # [(data file as written, number of data rows)]
    dataset: np.ndarray   # the data set whose window keeps the row ("" when none)
    record: np.ndarray    # 1-based data row in its file
    energy: np.ndarray    # MeV, lab
    angle: np.ndarray     # deg
    value: np.ndarray     # engine unit
    error: np.ndarray     # engine unit, the error column
    fitted: np.ndarray    # bool: a data set fits the row
    unit: str             # the engine unit of value and error
    norms: dict           # data set -> (fitted, prior %, group)

    @property
    def n(self):
        return len(self.value)


def rows_of(project, datasets):
    """The rows a covariance of the data sets ``datasets`` (names of ``project``, a :class:`pyazr.Project`) covers:
    every data row of their data files, read by the engine (``_azure2.data_rows``), in the order of the files' first
    use.  A row two data sets fit (overlapping windows on one file) is refused, as the engine refuses it."""
    from .described import _engine
    if isinstance(datasets, str):
        datasets = [datasets]
    if not datasets:
        raise ValueError("rows_of: name at least one data set")
    # The rows are read with the data sets' covariance fields removed, so a missing or broken file can be rebuilt.
    from .project import Project
    copy = Project(project.to_bytes(), project.path)
    for name in datasets:
        if copy[name].get("covariance") is not None:
            copy.set(name, covariance=None)
    by_file, order, units, norms = {}, [], set(), {}
    for name in datasets:
        r = _engine().data_rows(copy.to_bytes(), copy.directory, name)
        units.add(r["engine_unit"])
        norms[name] = (bool(r["norm_fit"]), float(r["norm_prior_percent"]), r["norm_group"])
        key = os.path.normpath(r["path"])
        if key not in by_file:
            by_file[key] = (r, [])
            order.append(key)
        by_file[key][1].append((name, np.asarray(r["in_window"], bool)))
    if len(units) > 1:
        raise ValueError(f"rows_of: the data sets have different units ({', '.join(sorted(units))}); one covariance "
                         f"holds one quantity's numbers")
    parts = {k: [] for k in ("dataset", "record", "energy", "angle", "value", "error", "fitted")}
    files = []
    for key in order:
        r, users = by_file[key]
        n = len(r["value"])
        owner = np.array([""] * n, dtype=object)
        for name, kept in users:
            twice = kept & (owner != "")
            if np.any(twice):
                raise ValueError(f"rows_of: data rows {list(np.flatnonzero(twice)[:5] + 1)} of {r['data_file']} "
                                 f"are fitted by {owner[twice][0]} and {name}; the engine refuses a row fitted twice")
            owner[kept] = name
        files.append((r["data_file"], n))
        parts["dataset"].append(owner)
        parts["record"].append(np.asarray(r["record"], int))
        for k in ("energy", "angle", "value", "error"):
            parts[k].append(np.asarray(r[k], float))
        parts["fitted"].append(owner != "")
    cat = {k: np.concatenate(v) for k, v in parts.items()}
    return Rows(files=files, unit=units.pop(), norms=norms, **cat)


# -- the parts' terms ----------------------------------------------------------------------------------------------

@dataclass
class Term:
    """One term of the file: kind (diagonal, vector, relvector, block), the part it comes from (0-based), the rows
    (0-based, for a block), and its numbers (N for diagonal / vector / relvector; the block's n x n matrix)."""
    kind: str
    part: int
    values: np.ndarray
    rows: np.ndarray = None

    def dense(self, rows):
        """The term as an N x N matrix, the data values standing in for the model (as the engine's error column)."""
        n = rows.n
        if self.kind == "diagonal":
            return np.diag(self.values)
        if self.kind == "vector":
            return np.outer(self.values, self.values)
        if self.kind == "relvector":
            u = self.values * rows.value
            return np.outer(u, u)
        out = np.zeros((n, n))
        out[np.ix_(self.rows, self.rows)] = self.values
        return out

    def diagonal(self, rows):
        """C_ii of the term, the data standing in for the model."""
        if self.kind == "diagonal":
            return self.values.copy()
        if self.kind == "vector":
            return self.values ** 2
        if self.kind == "relvector":
            return (self.values * rows.value) ** 2
        out = np.zeros(rows.n)
        out[self.rows] = np.diag(self.values)
        return out


def _in_range(rows, part):
    lo = -np.inf if part.get("e_min", "") == "" else part["e_min"]
    hi = np.inf if part.get("e_max", "") == "" else part["e_max"]
    return (rows.energy >= lo) & (rows.energy <= hi)


def _scale_terms(rows, part, index, mask=None):
    f = np.full(rows.n, part["percent"] / 100.0)
    if mask is not None:
        f[~mask] = 0.0
    if part["basis"] == "model":
        return [Term("relvector", index, f)]
    return [Term("vector", index, f * rows.value)]


def _own_terms(rows, part, index, others_diagonal):
    if part["source"] == "column":
        return [Term("diagonal", index, rows.error ** 2)]
    if part["source"] == "percent":
        return [Term("diagonal", index, (part["percent"] / 100.0 * rows.value) ** 2)]
    left = rows.error ** 2 - others_diagonal
    bad = np.flatnonzero(rows.fitted & (left <= 0))
    if bad.size:
        raise ValueError(f"part 'own' (remainder): the other parts exceed the error column on {bad.size} fitted "
                         f"row(s), first record {rows.record[bad[0]]} (E_lab {rows.energy[bad[0]]:g} MeV: error "
                         f"{rows.error[bad[0]]:g}, other parts {np.sqrt(others_diagonal[bad[0]]):g})")
    return [Term("diagonal", index, np.where(left > 0, left, 0.0))]


def _length_terms(rows, part, index):
    idx = np.flatnonzero(_in_range(rows, part))
    if idx.size == 0:
        return []
    s = part["percent"] / 100.0 * np.abs(rows.value[idx])
    d = np.abs(rows.energy[idx, None] - rows.energy[None, idx]) / part["length"]
    rho = np.exp(-d) if part["shape"] == "exponential" else np.exp(-0.5 * d * d)
    return [Term("block", index, s[:, None] * rho * s[None, :], idx)]


def _correlation_terms(rows, part, index, directory):
    path = part["file"] if is_absolute(part["file"]) else os.path.join(directory, part["file"])
    try:
        rho = np.loadtxt(path, comments="#", ndmin=2)
    except OSError as e:
        raise ValueError(f"part 'correlation': cannot read {part['file']} ({e})") from None
    except ValueError as e:
        raise ValueError(f"part 'correlation': {part['file']} is not a table of numbers ({e})") from None
    if rho.shape != (rows.n, rows.n):
        raise ValueError(f"part 'correlation': {part['file']} holds a {rho.shape[0]} x {rho.shape[1]} matrix; the "
                         f"covariance covers {rows.n} rows (every data row of {', '.join(f for f, _ in rows.files)})")
    if np.max(np.abs(rho - rho.T)) > 1e-10:
        raise ValueError(f"part 'correlation': {part['file']} is not symmetric")
    if np.max(np.abs(np.diag(rho) - 1.0)) > 1e-10 or np.max(np.abs(rho)) > 1.0 + 1e-12:
        raise ValueError(f"part 'correlation': {part['file']} is not a correlation matrix (diagonal 1, |entries| <= 1)")
    s = rows.error if part["percent"] == "" else part["percent"] / 100.0 * np.abs(rows.value)
    return [Term("block", index, s[:, None] * rho * s[None, :], np.arange(rows.n))]


def _terms(rows, parts, directory):
    """Every term, in the order of the parts; an ``own`` part with source remainder is made last (it needs the
    others' diagonal)."""
    terms, later = [], []
    for i, part in enumerate(parts):
        k = part["kind"]
        if k == "own" and part["source"] == "remainder":
            later.append(i)
        elif k == "own":
            terms += _own_terms(rows, part, i, None)
        elif k == "scale":
            terms += _scale_terms(rows, part, i)
        elif k == "range":
            terms += _scale_terms(rows, part, i, _in_range(rows, part))
        elif k == "length":
            terms += _length_terms(rows, part, i)
        elif k == "correlation":
            terms += _correlation_terms(rows, part, i, directory)
    if len(later) > 1:
        raise ValueError("only one part 'own' with source remainder: it is what the column leaves of the others")
    for i in later:
        others = sum((t.diagonal(rows) for t in terms), np.zeros(rows.n))
        terms += _own_terms(rows, parts[i], i, others)
    return terms


# -- the result ----------------------------------------------------------------------------------------------------

@dataclass
class Covariance:
    """A covariance built from parts (:func:`build`).

    - ``rows``: the rows covered (:class:`Rows`); ``parts``: the parts, normalized; ``terms``: the file's terms.
    - ``variance``: C_ii of every row (engine unit squared), the data standing in for the model;
      ``total_error`` its square root, ``error_column`` the file's, ``ratio`` their quotient.
    - ``fitted``: the 0-based rows the data sets fit; ``matrix`` / ``correlation``: C and its correlation on them
      (None above :data:`DENSE_LIMIT` rows).
    - ``checks``: ``symmetric``, ``positive_definite`` (C on the fitted rows), ``own_variance`` (the diagonal and
      block terms alone, which the engine needs positive definite), ``smallest_eigenvalue``, ``sigma_mismatch``
      (fitted rows whose error column differs from sqrt(C_ii) by more than 1 %, as the engine counts them) and
      ``worst_mismatch``; ``refusal``: why the engine would refuse it (None when it would not), with ``refusal_parts``
      the parts concerned.
    - ``notes``: what the engine does with it that the user should know (sentences).
    """
    rows: Rows
    parts: list
    terms: list
    variance: np.ndarray
    fitted: np.ndarray
    matrix: np.ndarray = None
    checks: dict = field(default_factory=dict)
    refusal: str = None
    refusal_parts: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    datasets: list = field(default_factory=list)

    @property
    def total_error(self):
        return np.sqrt(self.variance)

    @property
    def error_column(self):
        return self.rows.error

    @property
    def ratio(self):
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.total_error / self.rows.error

    @property
    def correlation(self):
        if self.matrix is None:
            return None
        s = np.sqrt(np.diag(self.matrix))
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.matrix / np.outer(s, s)

    def text(self, comment=None):
        """The file's text: a comment, the parts' line, the ``file`` lines and the terms (17 significant digits)."""
        num = lambda a: " ".join(f"{float(v):.17g}" for v in a)
        lines = [f"# data covariance of {', '.join(self.datasets)}, built from {len(self.parts)} part(s) by "
                 f"pyazr.covariance; numbers in {self.rows.unit or 'the engine unit'} (squared)"]
        lines += [f"# {c}" for c in (comment or "").splitlines()]
        lines.append(MARKER + " " + json.dumps({"datasets": self.datasets, "parts": self.parts},
                                               separators=(",", ":")))
        lines += [f"file {f} {n}" for f, n in self.rows.files]
        for t in self.terms:
            lines.append(f"# part {t.part + 1} ({self.parts[t.part]['kind']})")
            if t.kind == "block":
                lines.append(f"block {len(t.rows)}")
                lines.append(" ".join(str(int(i) + 1) for i in t.rows))
                lines += [num(r) for r in t.values]
            else:
                lines.append(t.kind)
                lines.append(num(t.values))
        return "\n".join(lines) + "\n"

    def write(self, path, comment=None):
        """Write the file the engine reads (refused while :attr:`refusal` is set); returns ``path``."""
        if self.refusal:
            raise ValueError(f"the covariance is not written: {self.refusal}")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(self.text(comment))
        return path

    def summary(self):
        """Plain data for a caller that draws it (the web page): rows, totals, checks, notes, the matrix's
        correlation on the fitted rows."""
        r = self.rows
        corr = self.correlation
        return {"unit": r.unit, "files": [list(f) for f in r.files], "datasets": self.datasets, "parts": self.parts,
                "kinds": {k: {"doc": v["doc"], "fields": [dict(zip(("name", "type", "default", "unit", "doc"), f))
                                                          for f in v["fields"]]} for k, v in PART_KINDS.items()},
                "rows": {"dataset": [str(x) for x in r.dataset], "record": r.record.tolist(),
                         "energy": r.energy.tolist(), "value": r.value.tolist(), "error": r.error.tolist(),
                         "total_error": self.total_error.tolist(), "fitted": r.fitted.tolist()},
                "terms": [{"kind": t.kind, "part": t.part} for t in self.terms],
                "checks": self.checks, "refusal": self.refusal, "refusal_parts": self.refusal_parts,
                "notes": self.notes,
                "correlation": None if corr is None else {"rows": self.fitted.tolist(), "matrix": corr.tolist()}}


def build(project, datasets, parts):
    """The covariance of the data sets ``datasets`` (names of the :class:`pyazr.Project` ``project``) made from
    ``parts`` (dicts, :data:`PART_KINDS`): its terms, C_ii, the checks and the notes (:class:`Covariance`).  A part
    that cannot be made raises ValueError; a matrix the engine would refuse comes back with ``refusal`` set."""
    if isinstance(datasets, str):
        datasets = [datasets]
    parts = [_normal_part(p) for p in parts]
    if not parts:
        raise ValueError("build: give at least one part (an 'own' part: the engine needs every row's own variance)")
    rows = rows_of(project, datasets)
    terms = _terms(rows, parts, project.directory)
    variance = sum((t.diagonal(rows) for t in terms), np.zeros(rows.n))
    fitted = np.flatnonzero(rows.fitted)
    cov = Covariance(rows=rows, parts=parts, terms=terms, variance=variance, fitted=fitted, datasets=list(datasets))
    _check(cov)
    _notes(cov)
    return cov


def _check(cov):
    rows, f = cov.rows, cov.fitted
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(rows.error[f] - np.sqrt(cov.variance[f])) / np.sqrt(cov.variance[f])
    cov.checks["sigma_mismatch"] = int(np.sum(rel > 0.01))
    cov.checks["worst_mismatch"] = float(np.max(rel)) if rel.size else 0.0
    cov.checks["fitted_rows"] = int(f.size)
    zero = f[cov.variance[f] <= 0]
    if zero.size:
        cov.refusal = (f"{zero.size} fitted row(s) have no variance (first: record {rows.record[zero[0]]} of "
                       f"{rows.dataset[zero[0]]}); the engine refuses a variance <= 0")
        cov.refusal_parts = []
        return
    if f.size > DENSE_LIMIT:
        cov.checks["dense"] = f"not made: {f.size} fitted rows (more than {DENSE_LIMIT})"
        return
    own_terms = [t for t in cov.terms if t.kind in ("diagonal", "block")]
    B = sum((t.dense(rows)[np.ix_(f, f)] for t in own_terms), np.zeros((f.size, f.size)))
    C = sum((t.dense(rows)[np.ix_(f, f)] for t in cov.terms), np.zeros((f.size, f.size)))
    cov.matrix = C
    cov.checks["symmetric"] = bool(np.allclose(C, C.T, rtol=0, atol=1e-12 * np.max(np.abs(C))))
    eig = np.linalg.eigvalsh(C)
    cov.checks["smallest_eigenvalue"] = float(eig[0])
    cov.checks["largest_eigenvalue"] = float(eig[-1])
    cov.checks["positive_definite"] = _cholesky_ok(C)
    cov.checks["own_variance"] = _cholesky_ok(B)
    # A part is a covariance itself: a block that is not positive semi-definite is refused even where the others
    # hide it (the engine would take the sum).
    concerned = []
    for t in own_terms:
        if t.kind == "block":
            m = t.dense(rows)[np.ix_(f, f)]
            if np.linalg.eigvalsh(m)[0] < -1e-12 * max(np.max(np.abs(m)), 1e-300):
                concerned.append(t.part)
    if concerned:
        cov.refusal = ("not positive semi-definite on the fitted rows (a part must be a covariance itself): part(s) "
                       + ", ".join(f"{p + 1} ({cov.parts[p]['kind']})" for p in sorted(set(concerned))))
    elif not cov.checks["own_variance"]:
        concerned = [p for p, part in enumerate(cov.parts) if part["kind"] in ("scale", "range")]
        cov.refusal = ("the diagonal and block parts alone are not positive definite on the fitted rows: the engine "
                       "needs every fitted row's own variance (add an 'own' part); fully correlated parts "
                       "(scale, range) cannot stand alone")
    elif not cov.checks["positive_definite"]:
        cov.refusal = "the matrix is not positive definite on the fitted rows"
    cov.refusal_parts = sorted(set(concerned))


def _cholesky_ok(m):
    if m.size == 0:
        return True
    try:
        np.linalg.cholesky(m)
        return True
    except np.linalg.LinAlgError:
        return False


def _notes(cov):
    notes = cov.notes
    notes.append("The engine replaces the error column of every fitted row by sqrt(C_ii) (the data value standing in "
                 "for the model) and warns when they differ by more than 1 %"
                 + (f": here {cov.checks['sigma_mismatch']} row(s), worst {100 * cov.checks['worst_mismatch']:.3g} %."
                    if cov.checks.get("sigma_mismatch") else ": here none does."))
    if any(t.kind == "relvector" for t in cov.terms):
        notes.append("relvector terms scale with the model, frozen at the first evaluation of a run (written to "
                     "covarianceModel.out); refreeze_covariance() or a new run from the fitted parameters updates it.")
    scales = [p for p in cov.parts if p["kind"] == "scale"]
    for name, (fit, prior, group) in cov.rows.norms.items():
        if scales and fit and prior > 0:
            sizes = ", ".join(f"{p['percent']:g}" for p in scales)
            notes.append(f"{name}: its normalization is fitted with a {prior:g} % prior and a scale part of "
                         f"{sizes} % is in the covariance: the same "
                         f"uncertainty counted twice; keep one (the prior, or the part with the normalization fixed).")
    unfitted = cov.rows.n - cov.fitted.size
    if unfitted:
        notes.append(f"{unfitted} of {cov.rows.n} rows are outside the data sets' windows: they are in the file "
                     f"(the file lines count every row) and drop out of the chi-squared.")


def read_parts(path):
    """``{"datasets": [...], "parts": [...]}`` of a file :meth:`Covariance.write` wrote (its ``#pyazr.covariance``
    line); None for a file without the line (written by hand or by another program)."""
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(MARKER + " "):
                return json.loads(line[len(MARKER) + 1:])
    return None
