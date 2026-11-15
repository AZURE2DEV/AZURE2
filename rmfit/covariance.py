"""Data covariance matrices (``covariance <file>``) and normalization groups
(``normgroup <label>``) for the rmfit evaluators.

AZURE2 (dev 4b71177, 2026-10-22) whitens the residual rows of the segments that
name one covariance file with ``W = (I + V T V^T) L_B^-1`` (``W^T W = M^-1``,
``M = D^-1 C D^-1``; ``src/DataCovariance.cpp``).  rmfit cannot use the engine's
whitened rows directly: a shard sees only its own rows (the marginal, with the
correlations to the other shards dropped), and the engine offers no way to whiten
an arbitrary vector, so the worker's one-evaluation residual path could not
follow.  Instead every worker evaluates a model with the ``covariance`` keywords
removed (plain normalized residuals, the error column of the data file), and the
evaluator whitens the assembled rows here, with the engine's own algorithm:

- the same rows (data-file records counted as ``ESegment::Fill`` counts them,
  the points actually fitted, in file-row order),
- the same ``D`` (``sigma_i = sqrt(C_ii)`` with the data value in the
  ``relvector`` terms; the plain rows are rescaled by ``sigma_file / sigma_i``),
- the same ``B`` components and Cholesky factors, the same ``V`` and ``T``,
- the same model freeze for ``relvector`` terms: ``m = d + sigma r`` at the
  FIRST evaluation of the evaluator (``refreeze()`` forgets it),

unless a fixed REFERENCE model is set (``Covariance.set_reference``, 2026-10-25, the evaluator: "use a fixed
reference model"): the ``relvector`` terms are then scaled by that model at every evaluation and
``refreeze()`` returns to it.  rmfit campaigns save the seed's frozen model at ``init``
(``cov_reference.npz``, ObjectiveSpec.cov_reference) so a polish, its verify and its export all see one
fixed objective; the self-consistent chi-squared (frozen at the evaluated point itself, what the engine
or a reference-free session reports) is ``Covariance.chi2_self``.

2026-10-26 (the evaluator: "build the fix and test it"): the reference is a CROSS SECTION (fit units, the frozen
data-unit model times the data normalization there), and the relvector terms use ``ref / n`` with ``n``
each row's CURRENT normalization.  With a data-unit reference, raising a free norm shrank every residual's
shape part as 1/n while the overall miss lay along the reference relvector directions at fixed variance:
one large model ran a norm from 0.63 to 3.08 (stopped only by its prior), its objective
a quarter below its self-consistent chi2.  Now the variances scale like the residuals; the factors are
rebuilt when a norm moves, and the norm columns of the Jacobian carry ``d(W)/dn y``
(``CovGroup.dwhiten_dnorm``, exact).  The first-evaluation freeze without a reference is unchanged (engine).

so the whitened rows equal the engine's ``residual_jacobian`` rows to rounding,
and the chi-squared equals the engine's (``SingleSession`` checks the latter
against an engine session on the unstripped model at its first evaluation).

Normalization groups: the engine gives every member the norm of the group's
first ACTIVE segment in file order (``EData::ResolveNormGroups``).  A single
session keeps the keyword (the engine does it); shards drop it, evaluate a
follower's pieces at the follower's own fixed norm ``n0`` and the evaluator
rescales them to the leader's norm ``n``: ``r = (n0/n)(r0 + q) - q`` with
``q = d/sigma`` per row.
"""
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# -- segment lines -----------------------------------------------------------------

#: The keyword blocks a segment line may end with and their token counts
#: (pyazr.datasets._KEYWORDS of AZURE2 dev; shiftfn takes a variable number).
_KEYWORDS = {"sqrtshift": 3, "sample": 1, "natural": 1, "shiftgroup": 1, "normgroup": 1,
             "covariance": 1, "shiftfn": None,
             "shiftfngroup": 1,   # feature/shift-fn-groups (2026-10-22), known ahead so it cannot hide what follows
             "total": 1}          # optical-theorem total cross section, AZURE2 dev d37d4e6 (2026-10-24)


def split_keyword_tail(tokens: Sequence[str]) -> Tuple[List[str], List[Tuple[str, List[str]]]]:
    """``(positional tokens, [(keyword, its tokens), ...])`` of a segment line, the
    blocks in their order (the parser of pyazr.datasets, kept here so rmfit reads
    the keywords with any installed pyazr)."""
    tokens = list(tokens)
    first = next((k for k, tok in enumerate(tokens) if tok.lower() in _KEYWORDS), None)
    if first is None:
        return tokens, []
    blocks, k = [], first
    while k < len(tokens):
        kw = tokens[k].lower()
        if kw not in _KEYWORDS:
            break
        if kw == "shiftfn":
            j = k + 1
            if j < len(tokens) and tokens[j].startswith('"'):
                while not (len(tokens[j]) > 1 and tokens[j].endswith('"')) and j + 1 < len(tokens):
                    j += 1
            try:
                n = int(tokens[j + 1])
            except (IndexError, ValueError):
                raise ValueError("malformed shiftfn block: " + " ".join(tokens))
            end = j + 2 + 3 * n
            blocks.append((kw, tokens[k + 1:end]))
            k = end
            continue
        n = _KEYWORDS[kw]
        if k + 1 + n > len(tokens):
            raise ValueError(f"{kw} block without its value: " + " ".join(tokens))
        blocks.append((kw, tokens[k + 1:k + 1 + n]))
        k += 1 + n
    return tokens[:first], blocks + ([("", tokens[k:])] if k < len(tokens) else [])


def line_keyword(tokens: Sequence[str], kw: str) -> Optional[str]:
    """The (first) value of keyword block ``kw`` on a segment line, or None."""
    for name, vals in split_keyword_tail(tokens)[1]:
        if name == kw and vals:
            return vals[0].strip('"')
    return None


def strip_keywords(tokens: Sequence[str], drop: Sequence[str] = ("covariance", "normgroup")) -> List[str]:
    """The segment line without the keyword blocks in ``drop`` (the others verbatim)."""
    pos, blocks = split_keyword_tail(tokens)
    out = list(pos)
    for name, vals in blocks:
        if name in drop:
            continue
        out += ([name] if name else []) + list(vals)
    return out


def norm_field(tokens: Sequence[str]) -> int:
    """Index of the norm value on a segment line (norm, vary flag and norm error follow)."""
    return 8 + (2 if int(float(tokens[7])) == 2 else 0)


def segment_lines(azr: str) -> List[List[str]]:
    """The tokens of every non-blank <segmentsData> line; line ``i`` is segment key ``i + 1``."""
    text = open(azr).read()
    a = text.find("<segmentsData>")
    b = text.find("</segmentsData>")
    if a < 0 or b < 0:
        raise ValueError(f"{azr}: no <segmentsData> block")
    body = text[text.index("\n", a) + 1:b]
    return [l.split() for l in body.splitlines() if l.strip()]


def has_keyword(azr: str, kw: str) -> bool:
    return any(line_keyword(t, kw) is not None for t in segment_lines(azr))


def strip_text(text: str, drop: Sequence[str] = ("covariance",)) -> str:
    """``text`` (a whole .azr) with the keyword blocks in ``drop`` removed from its data segment lines."""
    a = text.find("<segmentsData>")
    b = text.find("</segmentsData>")
    if a < 0 or b < 0:
        return text
    start = text.index("\n", a) + 1
    lines = []
    for l in text[start:b].splitlines():
        t = l.split()
        lines.append("  ".join(strip_keywords(t, drop)) if t and any(line_keyword(t, k) for k in drop) else l)
    return text[:start] + "".join(l + "\n" for l in lines) + text[b:]


@dataclass
class NormGroups:
    """``normgroup`` (or, with ``kw="shiftgroup"``, shift-group) membership of a model:
    ``leader[key]`` is the key whose norm (constant shift) a follower uses -- the group's
    first ACTIVE segment in file order (EData::ResolveNormGroups / ResolveShiftGroups)."""
    leader: Dict[int, int] = field(default_factory=dict)        # follower key -> leader key
    label: Dict[int, str] = field(default_factory=dict)         # member key -> label

    @classmethod
    def from_azr(cls, azr: str, kw: str = "normgroup") -> "NormGroups":
        out = cls()
        first: Dict[str, int] = {}
        for i, t in enumerate(segment_lines(azr), start=1):
            if not t or int(float(t[0])) != 1:
                continue
            lab = line_keyword(t, kw)
            if lab is None:
                continue
            out.label[i] = lab
            if lab in first:
                out.leader[i] = first[lab]
            else:
                first[lab] = i
        return out

    def owner(self, key: int) -> int:
        """The key whose norm segment ``key`` uses."""
        return self.leader.get(int(key), int(key))

    def __bool__(self):
        return bool(self.leader)

    def units(self) -> Dict[int, int]:
        """``{member key: leader key}`` for every member of a group with followers (leaders included)."""
        out = {f: l for f, l in self.leader.items()}
        out.update({l: l for l in self.leader.values()})
        return out


def shift_groups(azr: str) -> NormGroups:
    """``shiftgroup`` membership (same leader rule as normgroup)."""
    return NormGroups.from_azr(azr, "shiftgroup")


def shift_fn_groups(azr: str) -> NormGroups:
    """``shiftfngroup`` membership (AZURE2 dev 62392fc, 2026-10-24): every member uses the shift function of the
    group's first ACTIVE member, which must declare it (EData::ResolveShiftFunctions); same leader rule."""
    return NormGroups.from_azr(azr, "shiftfngroup")


# -- data files and covariance files -------------------------------------------------

def read_records(path: str) -> np.ndarray:
    """(energy, angle, value, error) of every data RECORD of a data file, counted as
    ``ESegment::Fill``/``DataCovariance::CountRecords`` count them: every line (CR also
    ends a line) that is neither blank nor starts with '#'."""
    rows = []
    with open(path, newline="") as f:
        for raw in f.read().replace("\r\n", "\n").replace("\r", "\n").split("\n"):
            s = raw.strip(" \t")
            if not s or s[0] == "#":
                continue
            t = s.split()
            try:
                rows.append([float(v) for v in t[:4]] + [np.nan] * (4 - len(t[:4])))
            except ValueError:
                raise ValueError(f"{path}: cannot parse data line {raw!r}")
    return np.array(rows, float).reshape(-1, 4)


def _canon(path: str, cwd: str) -> str:
    return os.path.realpath(os.path.join(cwd, path))


def _tokens(path: str):
    """Whitespace tokens with '#' comments; "quoted" tokens (DataCovariance::Tokens)."""
    with open(path) as f:
        for line in f:
            i, n = 0, len(line)
            while i < n:
                if line[i].isspace():
                    i += 1
                    continue
                if line[i] == "#":
                    break
                if line[i] == '"':
                    j = line.find('"', i + 1)
                    j = n if j < 0 else j
                    yield line[i + 1:j]
                    i = j + 1
                    continue
                j = i
                while j < n and not line[j].isspace() and line[j] != "#":
                    j += 1
                yield line[i:j]
                i = j


@dataclass
class CovFile:
    """A parsed covariance file (DataCovariance::File)."""
    path: str
    files: List[str]                 # canonical data-file paths
    counts: List[int]
    offsets: List[int]
    N: int
    diag: np.ndarray
    blocks: List[Tuple[np.ndarray, np.ndarray]]      # (0-based rows, m x m values)
    vectors: List[np.ndarray]
    rel: List[np.ndarray]

    def row_of(self, data_file: str, record: int) -> int:
        """0-based matrix row of 1-based ``record`` of canonical ``data_file`` (-1: none)."""
        if data_file not in self.files:
            return -1
        k = self.files.index(data_file)
        return self.offsets[k] + record - 1 if 1 <= record <= self.counts[k] else -1

    def diag_b(self) -> np.ndarray:
        d = self.diag.copy()
        for rows, vals in self.blocks:
            d[rows] += np.diag(vals)
        return d

    def variance(self, rows: np.ndarray, model: np.ndarray) -> np.ndarray:
        """C_ii with ``model`` standing in for the model in the relvector terms."""
        c = self.diag_b()[rows]
        for v in self.vectors:
            c = c + v[rows] ** 2
        for f in self.rel:
            c = c + (f[rows] * model) ** 2
        return c


def read_covariance(path: str, cwd: str) -> CovFile:
    """Parse a covariance file as DataCovariance::Read does (paths relative to ``cwd``)."""
    full = os.path.join(cwd, path)
    if not os.path.exists(full):
        raise FileNotFoundError(f"covariance file {path} not found from {cwd} (rmfit's shard and variant "
                                f"directories link only data/: keep covariance files under data/)")
    tok = _tokens(full)
    files, counts, offsets, N = [], [], [], 0
    diag = None
    blocks, vectors, rel = [], [], []

    def num():
        t = next(tok, None)
        if t is None:
            raise ValueError(f"{path}: the file ended inside a term")
        return float(t)

    def nums(n):
        return np.array([num() for _ in range(n)], float)

    def sym(a, what):
        s = np.sqrt(np.abs(np.outer(np.diag(a), np.diag(a))))
        s = np.maximum(s, np.maximum(np.abs(a), np.abs(a.T)))
        if np.any(np.abs(a - a.T) > 1e-6 * s):
            raise ValueError(f"{path}: {what} is not symmetric")
        return 0.5 * (a + a.T)

    for t in tok:
        if t == "file":
            if diag is not None:
                raise ValueError(f"{path}: every file line must come before the matrix terms")
            name = next(tok)
            c = int(num())
            canon = _canon(name, cwd)
            if canon in files:
                raise ValueError(f"{path}: data file {name} is listed twice")
            actual = read_records(os.path.join(cwd, name)).shape[0]
            if actual != c:
                raise ValueError(f"{path}: data file {name} has {actual} data rows, the header says {c}")
            files.append(canon); counts.append(c); offsets.append(N); N += c
        elif t == "size":
            if int(num()) != N:
                raise ValueError(f"{path}: size does not match the {N} rows of the file lines")
        elif t in ("block", "diagonal", "vector", "relvector", "matrix"):
            if N == 0:
                raise ValueError(f"{path}: the matrix terms must follow the file lines")
            if diag is None:
                diag = np.zeros(N)
            if t == "block":
                m = int(num())
                rows = nums(m).astype(int) - 1
                if np.any(rows < 0) or np.any(rows >= N) or len(set(rows.tolist())) != m:
                    raise ValueError(f"{path}: bad block row indices")
                blocks.append((rows, sym(nums(m * m).reshape(m, m), "block")))
            elif t == "diagonal":
                diag += nums(N)
            elif t == "vector":
                vectors.append(nums(N))
            elif t == "relvector":
                rel.append(nums(N))
            else:
                blocks.append((np.arange(N), sym(nums(N * N).reshape(N, N), "matrix")))
        else:
            raise ValueError(f"{path}: unknown term {t!r}")
    if N == 0 or diag is None:
        raise ValueError(f"{path}: no file lines or no matrix terms")
    return CovFile(path, files, counts, offsets, N, diag, blocks, vectors, rel)


# -- one covariance group ----------------------------------------------------------------

class CovGroup:
    """The fitted rows of the segments naming one covariance file, in file-row order,
    and the whitening of DataCovariance::Group."""

    def __init__(self, cf: CovFile, rows: np.ndarray, pos: np.ndarray, data: np.ndarray,
                 sigma_file: np.ndarray, label: str):
        order = np.argsort(rows, kind="stable")
        self.cf, self.path, self.label = cf, cf.path, label
        self.rows, self.pos = rows[order], pos[order]
        self.data, self.sigma_file = data[order], sigma_file[order]
        self.n = int(self.rows.size)
        if np.unique(self.rows).size != self.n:
            raise ValueError(f"{cf.path}: a data row is fitted by two segments (overlapping windows on one file)")
        if not np.all(self.sigma_file > 0):
            bad = int(np.argmax(~(self.sigma_file > 0)))
            raise ValueError(f"{cf.path}: row {self.rows[bad] + 1} has a zero or missing error column; rmfit "
                             f"evaluates covariance segments without the keyword and needs it nonzero")
        var = cf.variance(self.rows, self.data)
        if not np.all(var > 0):
            raise ValueError(f"{cf.path}: a fitted row has a variance <= 0")
        self.sigma = np.sqrt(var)                 # D of the engine
        self.scale = self.sigma_file / self.sigma  # plain rows -> the engine's normalized residuals
        self._setup()
        self.model = None                          # the relvector model of the current factors (data units)
        self.model_norm = None                     # the data normalizations at the first-evaluation freeze
        self.ref = None                            # a fixed reference CROSS SECTION (set_reference: model x norm,
                                                   # fit units) or None; the factors use ref / current norm
        self.fnorm = None                          # the normalizations the factors were built for (reference mode)
        self.V = self.T = None

    def _setup(self):
        n, cf = self.n, self.cf
        where = np.full(cf.N, -1)
        where[self.rows] = np.arange(n)
        parent = np.arange(n)

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        for brows, _ in cf.blocks:
            mem = where[brows]
            mem = mem[mem >= 0]
            for a in mem[1:]:
                x, y = find(mem[0]), find(a)
                if x != y:
                    parent[max(x, y)] = min(x, y)
        roots = np.array([find(a) for a in range(n)])
        comps: Dict[int, List[int]] = {}
        for a in range(n):
            comps.setdefault(int(roots[a]), []).append(a)
        s = self.sigma
        single, multi = [], []
        for idx in comps.values():
            (single if len(idx) == 1 else multi).append(np.array(idx))
        # B' = D^-1 B D^-1 per component
        self.single = np.array([i[0] for i in single], int)
        # a lone member's B'_aa: the diagonal plus every block term on its row (a block can hold
        # no other fitted member, or the two would share a component)
        self.single_l = (np.sqrt(cf.diag_b()[self.rows[self.single]] / s[self.single] ** 2)
                         if self.single.size else np.zeros(0))
        if np.any(~(self.single_l > 0)):
            raise ValueError(f"{cf.path}: B is not positive definite (a row without a variance of its own)")
        loc = {}
        for c, idx in enumerate(multi):
            for q, a in enumerate(idx):
                loc[int(a)] = (c, q)
        A = [np.diag(cf.diag[self.rows[idx]] / s[idx] ** 2) for idx in multi]
        for brows, vals in cf.blocks:
            mem = where[brows]
            ok = np.nonzero(mem >= 0)[0]
            if ok.size == 0 or int(mem[ok[0]]) not in loc:
                continue                    # no fitted row, or a lone member (in single_l)
            a = mem[ok]
            c = loc[int(a[0])][0]           # one block's fitted rows share one component
            q = np.array([loc[int(e)][1] for e in a])
            A[c][np.ix_(q, q)] += vals[np.ix_(ok, ok)] / np.outer(s[a], s[a])
        self.multi = multi
        self.multi_l = []
        for idx, Ac in zip(multi, A):
            try:
                self.multi_l.append(np.linalg.cholesky(Ac))
            except np.linalg.LinAlgError:
                raise ValueError(f"{cf.path}: B is not positive definite on the fitted rows {self.rows[idx][:5] + 1}...")
        self.absolute = [v[self.rows] for v in cf.vectors if np.any(v[self.rows] != 0)]
        self.rel = [f[self.rows] for f in cf.rel if np.any(f[self.rows] != 0)]

    # L_B^-1 / L_B^-T on a vector or a matrix (rows = members)
    def _solve_b(self, y, trans=False):
        from scipy.linalg import solve_triangular
        y = np.array(y, float, copy=True)
        if self.single.size:
            l = self.single_l if y.ndim == 1 else self.single_l[:, None]
            y[self.single] = y[self.single] / l
        for idx, L in zip(self.multi, self.multi_l):
            y[idx] = solve_triangular(L, y[idx], lower=True, trans=1 if trans else 0, check_finite=False)
        return y

    def _low_rank(self, y):
        if self.V is None or self.V.shape[1] == 0:
            return y
        return y + self.V @ (self.T @ (self.V.T @ y))

    @property
    def relative(self) -> bool:
        return bool(self.rel)

    @staticmethod
    def _cfun(s):
        """``c(s) = (1/sqrt(1+s) - 1)/s`` (T = P c(S) P^T) and its derivative, stable at small s."""
        s = np.asarray(s, float)
        small = s < 1e-8
        with np.errstate(divide="ignore", invalid="ignore"):
            g = 1.0 / np.sqrt(1.0 + s) - 1.0
            c = np.where(small, -0.5 + 0.375 * s, g / s)
            dc = np.where(small, 0.375 - 0.625 * s, (-0.5 * (1.0 + s) ** -1.5 * s - g) / s ** 2)
        return c, dc

    def _factors(self, model: np.ndarray, full: bool = False):
        U = [u for u in self.absolute] + [f * model for f in self.rel]
        if not U:
            V, T = np.zeros((self.n, 0)), np.zeros((0, 0))
            return (V, T, None) if full else (V, T)
        V = self._solve_b(np.column_stack(U) / self.sigma[:, None])
        s, P = np.linalg.eigh(V.T @ V)
        s = np.maximum(s, 0.0)
        c, _ = self._cfun(s)
        T = (P * c) @ P.T
        return (V, T, (s, P, c)) if full else (V, T)

    def factor(self, model: np.ndarray, norm: Optional[np.ndarray] = None):
        """Build V and T for the model ``m`` (data units, as the data file holds it; used by the relvector terms);
        ``norm``: the data normalizations it goes with (recorded for the reference)."""
        self.V, self.T, self._eig = self._factors(model, full=True)
        self.model = np.asarray(model, float).copy()
        self.fnorm = None if norm is None else np.asarray(norm, float).copy()

    def ensure_factored(self, y, norm: Optional[np.ndarray] = None):
        """Factors for the rows now (y: engine-scaled rows; norm: each row's data normalization now, None = 1).
        Reference mode (2026-10-26): the relvector model is the fixed reference CROSS SECTION divided by the
        current normalization, ``ref / n`` -- refactored whenever a normalization moves, so raising a norm
        shrinks these variances together with the residuals (no free direction).  Without a reference: the
        engine's first-evaluation freeze of the data-unit model ``m = d + sigma y`` (the normalizations there
        are recorded, :meth:`Covariance.frozen_models`)."""
        n = np.ones(self.n) if norm is None else np.asarray(norm, float)
        if self.ref is not None:
            if self.V is None or self.fnorm is None or not np.array_equal(n, self.fnorm):
                self.factor(self.ref / n, n)
        elif self.V is None:
            self.factor(self.data + self.sigma * y, n)
            self.model_norm = n.copy()

    def set_reference(self, xs: Optional[np.ndarray]):
        """Scale the relvector terms by this fixed CROSS SECTION (fit units: the reference model in data units
        times the data normalization there; this group's rows) divided by the current normalization, from now
        on; None returns to the first-evaluation freeze."""
        if xs is None:
            self.ref = None
            self.refreeze()
            return
        xs = np.asarray(xs, float)
        if xs.shape != (self.n,) or not np.all(np.isfinite(xs)):
            raise ValueError(f"covariance {self.path}: reference has shape {xs.shape} or non-finite "
                             f"values; the group has {self.n} rows")
        self.ref = xs.copy()
        self.V = self.T = None
        self.fnorm = None

    def refreeze(self):
        if self.relative:
            if self.ref is not None:
                return                                 # a fixed reference never follows the model
            self.V = self.T = None
            self.model = None
            self.model_norm = None

    def _components_of(self, mask: np.ndarray):
        """``(single positions, multi indices, Q)``: the B components holding a row of ``mask``, and every
        member row of them (Q, a bool mask)."""
        if not hasattr(self, "_comp"):
            comp = np.full(self.n, -1, int)
            comp[self.single] = -2 - np.arange(self.single.size)          # singles: -2, -3, ...
            for c, idx in enumerate(self.multi):
                comp[idx] = c
            self._comp = comp
        hit = np.unique(self._comp[mask])
        multi = [int(c) for c in hit if c >= 0]
        single = np.array([-2 - int(c) for c in hit if c <= -2], int)
        Q = np.isin(self._comp, hit)
        return single, multi, Q

    def _solve_b_on(self, Y, single, multi):
        """L_B^-1 on the rows of the listed components only (other rows returned unchanged)."""
        from scipy.linalg import solve_triangular
        Y = np.array(Y, float, copy=True)
        if single.size:
            rows = self.single[single]
            l = self.single_l[single] if Y.ndim == 1 else self.single_l[single][:, None]
            Y[rows] = Y[rows] / l
        for c in multi:
            idx = self.multi[c]
            Y[idx] = solve_triangular(self.multi_l[c], Y[idx], lower=True, check_finite=False)
        return Y

    def dwhiten_dnorm(self, y, mask: np.ndarray, n_c: float) -> np.ndarray:
        """d(W y)/dn for the normalization ``n`` shared by the rows of ``mask`` (reference mode: the relvector
        model there is ``ref / n``); the factors must be current (:meth:`ensure_factored`).  Exact: dV from
        dU = -U_rel / n on those rows, dT by the Daleckii-Krein formula in the eigenbasis of V^T V."""
        if self.ref is None or self.V is None or self.V.shape[1] == 0 or not np.any(mask):
            return np.zeros(self.n)
        s, P, c = self._eig
        V = self.V
        k_abs = len(self.absolute)
        single, multi, Q = self._components_of(mask)
        dU = np.zeros((self.n, V.shape[1]))
        mrow = np.where(mask, -self.model / n_c, 0.0)                 # d(ref/n)/dn on the rows of mask
        for j, f in enumerate(self.rel):
            dU[:, k_abs + j] = f * mrow
        dV = self._solve_b_on(dU / self.sigma[:, None], single, multi)
        dV[~Q] = 0.0
        VQ, dVQ = V[Q], dV[Q]
        dG = dVQ.T @ VQ
        dG = dG + dG.T
        Gt = P.T @ dG @ P
        _, dc = self._cfun(s)
        ds = s[:, None] - s[None, :]
        close = np.abs(ds) < 1e-6 * (1.0 + np.abs(s[:, None]) + np.abs(s[None, :]))
        with np.errstate(divide="ignore", invalid="ignore"):
            Gam = np.where(close, self._cfun(0.5 * (s[:, None] + s[None, :]))[1], (c[:, None] - c[None, :]) / ds)
        z = self._solve_b(y)
        a = V.T @ z
        da = dVQ.T @ z[Q]
        Ta = self.T @ a
        dTa = P @ ((Gt * Gam) @ (P.T @ a))            # dT a without forming dT (2026-10-26: saves 2 k^3 per column)
        return dV @ Ta + V @ (dTa + self.T @ da)

    def chi2_self(self, y) -> float:
        """Chi-squared of these rows (engine-scaled) with the relvector terms frozen at their OWN model
        ``d + sigma y`` (what a fresh reference-free session would report); the state is not touched."""
        V, T = self._factors(self.data + self.sigma * y)
        z = self._solve_b(y)
        if V.shape[1]:
            z = z + V @ (T @ (V.T @ z))
        return float(np.sum(z ** 2))

    def whiten(self, y):
        return self._low_rank(self._solve_b(y))

    def whiten_t(self, w):
        return self._solve_b(self._low_rank(w), trans=True)


# -- the evaluator-side transform ------------------------------------------------------------

@dataclass
class RowBlock:
    """Evaluator rows ``[a, b)`` holding the points of base segment ``key`` with lab
    energy in ``[e_lo, e_hi]``, in file order (``order="file"``, engine pieces) or
    sorted by energy (``"energy"``: conv pieces, ``rmfit.model.rows_in_segment``)."""
    a: int
    b: int
    key: int
    e_lo: float
    e_hi: float
    order: str = "file"


def segment_records(seg, path: str, e_lo: float, e_hi: float, order: str = "file") -> np.ndarray:
    """1-based records of ``path`` that the engine loads for segment ``seg`` restricted
    to ``[e_lo, e_hi]`` (the filter of rmfit.model.rows_in_segment), in point order."""
    from .model import DIFFERENTIAL
    rec = read_records(path)
    e, ang = rec[:, 0], rec[:, 1]
    sel = (e >= e_lo) & (e <= e_hi) & (e >= seg.energy_min) & (e <= seg.energy_max)
    if seg.observable in DIFFERENTIAL:
        sel &= (ang >= seg.angle_min - 1e-6) & (ang <= seg.angle_max + 1e-6)
    idx = np.nonzero(sel)[0]
    if order == "energy":
        idx = idx[np.argsort(e[idx], kind="stable")]
    return idx + 1


class Covariance:
    """Every covariance group of a model, on an evaluator's row layout."""

    def __init__(self, groups: List[CovGroup], npts: int):
        self.groups = groups
        self.npts = npts
        self.rows = np.concatenate([g.pos for g in groups]) if groups else np.zeros(0, int)

    @classmethod
    def build(cls, azr: str, cwd: str, segments, blocks: Sequence[RowBlock], npts: int) -> Optional["Covariance"]:
        """The covariance groups of model ``azr`` (paths relative to ``cwd``, where the engine
        would run) on the evaluator rows ``blocks``; None if no active segment names one."""
        lines = segment_lines(azr)
        cov_of = {}
        for i, t in enumerate(lines, start=1):
            c = line_keyword(t, "covariance")
            if c is not None and int(float(t[0])) == 1:
                cov_of[i] = c
        if not cov_of:
            return None
        files: Dict[str, CovFile] = {}
        parts: Dict[str, List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]] = {}
        for blk in blocks:
            path = cov_of.get(int(blk.key))
            if path is None:
                continue
            canon = _canon(path, cwd)
            if canon not in files:
                files[canon] = read_covariance(path, cwd)
            cf = files[canon]
            seg = segments.by_key(int(blk.key))
            dfile = os.path.join(cwd, seg.data_file)
            rec = segment_records(seg, dfile, blk.e_lo, blk.e_hi, blk.order)
            if rec.size != blk.b - blk.a:
                raise RuntimeError(f"covariance: segment {blk.key} [{blk.e_lo}, {blk.e_hi}] has {rec.size} records "
                                   f"by rmfit's filter, the engine loaded {blk.b - blk.a}")
            rows = np.array([cf.row_of(_canon(seg.data_file, cwd), int(r)) for r in rec], int)
            if np.any(rows < 0):
                raise ValueError(f"covariance {path}: data file {seg.data_file} of segment {blk.key} is not among "
                                 f"its file lines (or a record is past its count)")
            vals = read_records(dfile)[rec - 1]
            parts.setdefault(canon, []).append((rows, np.arange(blk.a, blk.b), vals[:, 2], vals[:, 3]))
        groups = []
        for canon, ps in parts.items():
            groups.append(CovGroup(files[canon], np.concatenate([p[0] for p in ps]), np.concatenate([p[1] for p in ps]),
                                   np.concatenate([p[2] for p in ps]), np.concatenate([p[3] for p in ps]),
                                   label=files[canon].path))
        return cls(groups, npts) if groups else None

    def check_ratios(self, q_engine: np.ndarray, rtol=1e-7):
        """Row-mapping check: ``q_engine`` (value / error of every evaluator row as the
        engine holds it, c.m. frame) must equal data / error of the mapped records."""
        for g in self.groups:
            qe = q_engine[g.pos]
            qf = g.data / g.sigma_file
            ok = np.isfinite(qe)
            bad = ok & ~np.isclose(qe, qf, rtol=rtol, atol=1e-300)
            if np.any(bad):
                i = int(np.argmax(bad))
                raise RuntimeError(f"covariance {g.path}: evaluator row {g.pos[i]} does not hold file row "
                                   f"{g.rows[i] + 1} (value/error {qe[i]!r} vs {qf[i]!r}); rmfit's row mapping "
                                   f"disagrees with the engine")

    def refreeze(self):
        for g in self.groups:
            g.refreeze()

    @property
    def relative(self) -> bool:
        return any(g.relative for g in self.groups)

    @property
    def has_reference(self) -> bool:
        return any(g.ref is not None for g in self.groups if g.relative)

    def frozen_models(self) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """``{covariance file basename: (file rows, frozen model, data normalizations)}`` of every relvector group
        (each must have been evaluated, i.e. frozen, already): the model in data units and the normalizations it
        was frozen with, so the reference cross section is ``model * norm`` (:func:`save_reference`)."""
        out = {}
        for g in self.groups:
            if not g.relative:
                continue
            if g.model is None:
                raise RuntimeError(f"covariance {g.path}: no frozen model yet (evaluate once first)")
            name = os.path.basename(g.path)
            if name in out:
                raise ValueError(f"two covariance groups share the file name {name}")
            norm = g.fnorm if g.ref is not None else g.model_norm
            if norm is None:
                raise RuntimeError(f"covariance {g.path}: the normalizations of the frozen model are unknown")
            out[name] = (g.rows.copy(), g.model.copy(), np.asarray(norm, float).copy())
        return out

    def set_reference(self, ref: Dict[str, Tuple[np.ndarray, ...]]):
        """Fix every relvector group's reference CROSS SECTION (model x normalization, see :meth:`frozen_models`),
        matched by covariance file name and file row.  Every fitted row must be in the reference.  A reference
        without normalizations (written before 2026-10-26, data-unit model only) is refused: with it a free
        normalization could shrink the residuals at no cost (a norm once ran up fivefold)."""
        for g in self.groups:
            if not g.relative:
                continue
            name = os.path.basename(g.path)
            if name not in ref:
                raise ValueError(f"covariance reference has no model for {name} (it has {sorted(ref)})")
            item = ref[name]
            if len(item) < 3 or item[2] is None:
                raise ValueError(f"covariance reference for {name} has no data normalizations (an rmfit reference "
                                 f"written before 2026-10-26, a data-unit model only): it cannot be scaled with the "
                                 f"normalizations; re-initialize the campaign (rmfit init) to write a new one")
            rows, model, norm = item[0], item[1], item[2]
            at = {int(r): i for i, r in enumerate(np.asarray(rows, int))}
            miss = [int(r) for r in g.rows if int(r) not in at]
            if miss:
                raise ValueError(f"covariance reference for {name} lacks {len(miss)} fitted row(s), e.g. file row "
                                 f"{miss[0] + 1}: it was saved from a fit with other data windows")
            sel = [at[int(r)] for r in g.rows]
            g.set_reference(np.asarray(model, float)[sel] * np.asarray(norm, float)[sel])

    def chi2_self(self, y_plain: np.ndarray) -> float:
        """Plain data chi-squared with every relvector term frozen at the evaluated model itself (rows of
        :meth:`apply`'s ``y_plain``, unweighted); equals a reference-free fresh session's."""
        y_plain = np.asarray(y_plain, float)
        mask = np.ones(y_plain.size, bool)
        tot = 0.0
        for g in self.groups:
            mask[g.pos] = False
            tot += g.chi2_self(y_plain[g.pos] * g.scale)
        return tot + float(np.sum(y_plain[mask] ** 2))


    def apply(self, y_plain: np.ndarray, J: Optional[np.ndarray] = None, sw: Optional[np.ndarray] = None,
              norms: Optional[np.ndarray] = None, norm_cols: Optional[np.ndarray] = None):
        """Whiten the data rows.

        ``y_plain``: the plain normalized data rows (error column of the data files);
        ``J`` their Jacobian or None; ``sw``: sqrt of the row weights (applied before
        whitening) or None; ``norms``: every evaluator row's data normalization now and ``norm_cols``
        the column of the free parameter it is (-1: fixed) -- needed with a fixed reference, whose
        relvector terms follow the normalization (and give the norm columns of ``Jw`` a second term,
        :meth:`CovGroup.dwhiten_dnorm`).  Returns ``(w, Jw, w_unweighted, contrib)``: the rows the
        objective sees, their Jacobian, the unweighted whitened rows (chi2 = their sum of
        squares) and per-row contributions (``r_i (M^-1 r)_i`` on covariance rows,
        ``r_i^2`` elsewhere; they add up to chi2 exactly)."""
        y_plain = np.asarray(y_plain, float)
        w0 = y_plain.copy()
        contrib = y_plain ** 2
        w = None if sw is None else y_plain * sw
        Jw = None
        if J is not None:
            Jw = J.copy() if sw is None else J * sw[:, None]
        for g in self.groups:
            y = y_plain[g.pos] * g.scale
            if g.ref is not None and norms is None:
                raise ValueError(f"covariance {g.path}: a fixed reference needs the rows' normalizations")
            g.ensure_factored(y, None if norms is None else np.asarray(norms, float)[g.pos])
            wg = g.whiten(y)
            w0[g.pos] = wg
            contrib[g.pos] = y * g.whiten_t(wg)
            if w is not None:
                w[g.pos] = g.whiten(y * sw[g.pos])
            if Jw is not None:
                Jg = Jw[g.pos] * g.scale[:, None]
                nz = np.flatnonzero(np.any(Jg != 0, axis=0))     # W of a zero column is zero (the engine skips them too)
                Jg[:, nz] = g.whiten(Jg[:, nz])
                Jw[g.pos] = Jg
                if g.ref is not None and norm_cols is not None:
                    # the relvector terms follow the normalization (ref / n): d(W y)/dn on its column
                    cols = np.asarray(norm_cols, int)[g.pos]
                    nrow = np.asarray(norms, float)[g.pos]
                    yd = y if sw is None else y * sw[g.pos]
                    for col in np.unique(cols[cols >= 0]):
                        m = cols == col
                        Jw[g.pos, col] += g.dwhiten_dnorm(yd, m, float(nrow[m][0]))
        return (w0 if w is None else w), Jw, w0, contrib

    def summary(self) -> str:
        out = []
        for g in self.groups:
            out.append(f"{os.path.basename(g.path)}: {g.n} rows, {len(g.multi)} block component(s) "
                       f"(largest {max([i.size for i in g.multi] + [1])}), {len(g.absolute)} vector + "
                       f"{len(g.rel)} relvector term(s)")
        return "; ".join(out)


def save_reference(path: str, ref: Dict[str, Tuple[np.ndarray, ...]], note: str = "") -> str:
    """Write a covariance reference (``Covariance.frozen_models``: rows, data-unit model, normalizations) as .npz."""
    arrays = {"names": np.array(sorted(ref)), "note": np.array(note), "format": np.array("xs-norm 2026-10-26")}
    for i, name in enumerate(sorted(ref)):
        item = ref[name]
        arrays[f"rows_{i}"], arrays[f"model_{i}"] = item[0], item[1]
        if len(item) > 2 and item[2] is not None:
            arrays[f"norm_{i}"] = item[2]
    np.savez(path, **arrays)
    return path


def load_reference(path: str) -> Dict[str, Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]]:
    """``{name: (rows, model, norm or None)}``; norm is None in a reference written before 2026-10-26."""
    with np.load(path) as z:
        return {str(n): (z[f"rows_{i}"], z[f"model_{i}"], z[f"norm_{i}"] if f"norm_{i}" in z.files else None)
                for i, n in enumerate(z["names"])}
