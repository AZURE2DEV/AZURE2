#!/usr/bin/env python3
"""Write, or check, a data covariance file of AZURE2 (the key `covariance <file>` of a data set; format in
docs/source/reference/data_formats.rst, "Data Covariance Files").

Write mode builds C, in the units of the data files' value column (lab frame, as in the files), from components:
  --stat               the files' error column squared, as a `diagonal` (on by default; --no-stat to leave out)
  --extra-diag F       a second `diagonal` (F * value)^2: scatter beyond the quoted errors
  --common F           one term of fraction F on every row of every file (a scale shared by the files)
  --common-per-file F  one term of fraction F per file
  --shape PATH         one term per file PATH: N fractions (one per row of all files, whitespace-separated)
  --matrix PATH        a full N x N matrix (.npy or text), absolute units, written as a `matrix` term
The fraction terms are model-scaled `relvector`s; --absolute writes them as data-scaled `vector`s instead (biased
low: Peelle's puzzle; use it only for a term that is absolute by nature).  An overall scale that a free
normalization with a prior already carries must not be in C as well.

Check mode (--check COV) reads a covariance file, forms C with the data values standing for the model in the
`relvector` terms (as the program does for the error column), and reports: the rows per file against the data
files, whether the block-diagonal part B (diagonal, block, matrix) and C are positive definite (B must be: the
program whitens with it), the smallest and largest eigenvalue, and the largest relative difference between
sqrt(C_ii) and the error column (the program warns above 1 %).

Reads:  the data files (paths as on the data set records, relative to --base, the project's directory).
Writes: OUT (write mode); a report on standard output.

    python3 write_covariance.py cov.txt --base . --file data/a.dat --file data/b.dat --common-per-file 0.02
    python3 write_covariance.py --check cov.txt --base .

Exit status: 0 written / checked and positive definite, 1 a check failed (not positive definite, row count),
2 bad arguments.
"""
import argparse
import os
import sys

import numpy as np


def data_rows(path, value_column=3, error_column=4):
    """(values, errors) of a data file, its rows counted as the program counts them: lines that are neither blank
    nor comments (#).  Columns are 1-based."""
    values, errors = [], []
    with open(path) as fh:
        for line in fh:
            t = line.split("#", 1)[0].split()
            if not t:
                continue
            values.append(float(t[value_column - 1]))
            errors.append(float(t[error_column - 1]) if len(t) >= error_column else np.nan)
    return np.array(values), np.array(errors)


def build_terms(files, base=".", stat=True, extra_diag=0.0, common=0.0, common_per_file=0.0, shapes=(),
                matrix=None, absolute=False, value_column=3, error_column=4):
    """The terms [(kind, array)] and (values, errors, counts) of the data files."""
    vals, errs, counts = [], [], []
    for f in files:
        v, e = data_rows(os.path.join(base, f), value_column, error_column)
        vals.append(v)
        errs.append(e)
        counts.append(len(v))
    v, e = np.concatenate(vals), np.concatenate(errs)
    n = len(v)
    terms = []
    frac_kind = "vector" if absolute else "relvector"

    def frac_term(fr):
        return (frac_kind, fr * v if absolute else fr)

    if stat:
        if not np.all(np.isfinite(e)) or np.any(e <= 0):
            raise ValueError("a data row has no positive error in the error column; give --no-stat and a matrix")
        terms.append(("diagonal", e ** 2))
    if extra_diag:
        terms.append(("diagonal", (extra_diag * v) ** 2))
    if common:
        terms.append(frac_term(np.full(n, common)))
    if common_per_file:
        start = 0
        for c in counts:
            fr = np.zeros(n)
            fr[start:start + c] = common_per_file
            terms.append(frac_term(fr))
            start += c
    for s in shapes:
        fr = np.loadtxt(s).ravel()
        if len(fr) != n:
            raise ValueError(f"{s} has {len(fr)} fractions; the files have {n} rows")
        terms.append(frac_term(fr))
    if matrix is not None:
        M = np.load(matrix) if str(matrix).endswith(".npy") else np.loadtxt(matrix)
        if M.shape != (n, n):
            raise ValueError(f"{matrix} is {M.shape}; the files have {n} rows")
        terms.append(("matrix", 0.5 * (M + M.T)))
    return terms, (v, e, counts)


def write(path, files, counts, terms, comment="written by write_covariance.py"):
    with open(path, "w") as fh:
        fh.write(f"# {comment}\n")
        for f, c in zip(files, counts):
            fh.write(f"file {f} {c}\n")
        fh.write(f"size {sum(counts)}\n")
        for kind, a in terms:
            fh.write(f"{kind}\n")
            a = np.atleast_2d(a) if kind == "matrix" else a[None, :]
            for row in a:
                for k in range(0, len(row), 6):
                    fh.write(" ".join(f"{x:.16e}" for x in row[k:k + 6]) + "\n")


def read(path):
    """(files [(name, rows)], terms [(kind, rows or None, array)]) of a covariance file."""
    tokens = []
    with open(path) as fh:
        for line in fh:
            tokens += line.split("#", 1)[0].split()
    files, terms, i = [], [], 0
    while i < len(tokens) and tokens[i] in ("file", "size"):
        if tokens[i] == "file":
            files.append((tokens[i + 1], int(tokens[i + 2])))
            i += 3
        else:
            i += 2
    n = sum(c for _, c in files)
    while i < len(tokens):
        kind = tokens[i]
        i += 1
        if kind == "block":
            k = int(tokens[i])
            rows = [int(t) - 1 for t in tokens[i + 1:i + 1 + k]]
            vals = np.array(tokens[i + 1 + k:i + 1 + k + k * k], float).reshape(k, k)
            terms.append((kind, rows, vals))
            i += 1 + k + k * k
        elif kind in ("diagonal", "vector", "relvector"):
            terms.append((kind, None, np.array(tokens[i:i + n], float)))
            i += n
        elif kind == "matrix":
            terms.append((kind, None, np.array(tokens[i:i + n * n], float).reshape(n, n)))
            i += n * n
        else:
            raise ValueError(f"{path}: unknown term {kind!r}")
    return files, terms


def dense(terms, n, model):
    """(B, C): the block-diagonal part and the whole matrix, `model` standing for m in the relvector terms."""
    B = np.zeros((n, n))
    for kind, rows, a in terms:
        if kind == "diagonal":
            B[np.diag_indices(n)] += a
        elif kind == "block":
            B[np.ix_(rows, rows)] += a
        elif kind == "matrix":
            B += a
    C = B.copy()
    for kind, rows, a in terms:
        if kind in ("vector", "relvector"):
            u = a * model if kind == "relvector" else a
            C += np.outer(u, u)
    return B, C


def positive_definite(M):
    try:
        np.linalg.cholesky(M)
        return True
    except np.linalg.LinAlgError:
        return False


def check(path, base=".", value_column=3, error_column=4, out=print):
    files, terms = read(path)
    ok = True
    vals, errs = [], []
    for f, c in files:
        v, e = data_rows(os.path.join(base, f), value_column, error_column)
        if len(v) != c:
            out(f"FAIL {f}: {len(v)} data rows, the covariance file says {c}")
            ok = False
        vals.append(v)
        errs.append(e)
    if not ok:
        return False
    v, e = np.concatenate(vals), np.concatenate(errs)
    B, C = dense(terms, len(v), v)
    lam = np.linalg.eigvalsh(C)
    pd_b, pd_c = positive_definite(B), positive_definite(C)
    rel = np.abs(np.sqrt(np.diag(C)) - e) / np.abs(e)
    out(f"{path}: {len(files)} files, {len(v)} rows, {len(terms)} terms")
    out(f"  B (diagonal, block, matrix) positive definite: {pd_b}; C positive definite: {pd_c}")
    out(f"  eigenvalues of C: {lam.min():.4g} .. {lam.max():.4g}")
    out(f"  sqrt(C_ii) against the error column: largest relative difference {np.nanmax(rel):.3g}"
        + (" (the program warns above 0.01: data and covariance files from different tables?)"
           if np.nanmax(rel) > 0.01 else ""))
    return pd_b and pd_c


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("out", nargs="?", help="the covariance file to write")
    ap.add_argument("--check", metavar="COV", help="check an existing covariance file instead")
    ap.add_argument("--base", default=".", help="directory the data file paths are relative to (the project's)")
    ap.add_argument("--file", action="append", default=[], help="data file, as on its data set record (repeatable)")
    ap.add_argument("--no-stat", dest="stat", action="store_false", help="leave out the error column's diagonal")
    ap.add_argument("--extra-diag", type=float, default=0.0)
    ap.add_argument("--common", type=float, default=0.0)
    ap.add_argument("--common-per-file", type=float, default=0.0)
    ap.add_argument("--shape", action="append", default=[])
    ap.add_argument("--matrix")
    ap.add_argument("--absolute", action="store_true", help="fraction terms as data-scaled vectors")
    ap.add_argument("--value-column", type=int, default=3)
    ap.add_argument("--error-column", type=int, default=4)
    a = ap.parse_args(argv)
    if a.check:
        return 0 if check(a.check, a.base, a.value_column, a.error_column) else 1
    if not a.out or not a.file:
        ap.error("give OUT and at least one --file (or --check COV)")
    try:
        terms, (v, e, counts) = build_terms(a.file, a.base, a.stat, a.extra_diag, a.common, a.common_per_file,
                                            a.shape, a.matrix, a.absolute, a.value_column, a.error_column)
    except (ValueError, OSError) as err:
        print(f"write_covariance: {err}", file=sys.stderr)
        return 2
    write(a.out, a.file, counts, terms)
    return 0 if check(a.out, a.base, a.value_column, a.error_column) else 1


if __name__ == "__main__":
    sys.exit(main())
