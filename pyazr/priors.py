"""Prior parameter covariance files: a multivariate Gaussian prior on fit parameters.

The ``<config>`` key ``priorCovariance=<path>`` (path relative to the ``.azr``'s directory) makes AZURE2 add

    (p - mu)^T C^-1 (p - mu)

to chi-squared over the parameters the file lists (``Total-Prior-Covariance-Chi-Squared`` in ``chiSquared.out``;
``penalties()["prior_covariance"]`` in a pyazr session).  It is SAMMY's "PUP" / sequential Bayesian update: fit A's
parameters and covariance become the prior of fit B, which then adds new data.

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

This module needs only numpy (no engine): :func:`prior_from_fit` writes the file from a previous fit's
``output/param.sav`` and ``output/covariance.dat``; :func:`write_prior` / :func:`read_prior` write and read any prior.
"""

import os

import numpy as np

__all__ = ["prior_from_fit", "read_prior", "write_prior", "read_covariance_dat"]


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


def write_prior(path, names, mean, cov, basis=None, comment=None):
    """Write a prior covariance file (names, means, covariance; optional ``#parametrization`` tag).

    ``cov`` must be square, symmetric and positive definite -- checked here as AZURE2 checks it, so a bad matrix is
    refused before a fit is started with it.
    """
    names = [str(n) for n in names]
    mean = np.asarray(mean, float).ravel()
    cov = np.asarray(cov, float)
    n = len(names)
    if mean.size != n or cov.shape != (n, n):
        raise ValueError(f"size mismatch: {n} names, {mean.size} means, covariance {cov.shape}")
    if len(set(names)) != n:
        raise ValueError("a parameter is listed twice")
    scale = np.sqrt(np.abs(np.outer(np.diag(cov), np.diag(cov))))
    if np.any(np.abs(cov - cov.T) > 1e-10 * np.maximum(scale, np.maximum(np.abs(cov), np.abs(cov.T)))):
        raise ValueError("the covariance is not symmetric")
    try:
        np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        raise ValueError("the covariance is not positive definite") from None
    with open(path, "w") as fh:
        if basis is not None:
            fh.write(f"#parametrization {int(basis)}\n")
        for c in (comment.splitlines() if comment else []):
            fh.write(f"# {c}\n")
        fh.write(" ".join(names) + "\n")
        fh.write(" ".join(f"{v:.17g}" for v in mean) + "\n")
        for row in cov:
            fh.write(" ".join(f"{v:.17g}" for v in row) + "\n")
    return path


def read_prior(path):
    """``{"names", "mean", "cov", "basis"}`` of a prior covariance file (basis None without a tag)."""
    basis, rows = None, []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                t = s[1:].split()
                if len(t) >= 2 and t[0] == "parametrization" and not rows:
                    basis = int(round(float(t[1])))
                continue
            t = s.split("#", 1)[0].split()
            if t:
                rows.append(t)
    if len(rows) < 2:
        raise ValueError(f"{path}: expected names, means and covariance rows")
    names = rows[0]
    n = len(names)
    mean = np.array([float(v) for v in rows[1]])
    cov = np.array([[float(v) for v in r] for r in rows[2:]])
    if mean.size != n or cov.shape != (n, n):
        raise ValueError(f"{path}: size mismatch ({n} names, {mean.size} means, covariance {cov.shape})")
    return {"names": names, "mean": mean, "cov": cov, "basis": basis}


def prior_from_fit(sav_path, covariance_path, out_path, names=None, columns=None):
    """Write a prior covariance file from a previous fit: means from its ``param.sav``, covariance from its
    ``covariance.dat``.

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
    vals, basis = _read_sav(sav_path)
    C, cols, scale = read_covariance_dat(covariance_path)
    if cols is None:
        if columns is None:
            raise ValueError(
                f"{covariance_path} has no '#columns' header (written before 2026-11-05): pass columns=[...], the "
                f"free energy_*/width_* names of that fit in param.sav order ({C.shape[0]} of them)")
        cols = list(columns)
        if len(cols) != C.shape[0]:
            raise ValueError(f"{len(cols)} column names for a {C.shape[0]}x{C.shape[0]} covariance")
    missing = [c for c in cols if c not in vals]
    if missing:
        raise ValueError(f"{sav_path} does not hold {', '.join(missing)} (the covariance belongs to another fit?)")
    if names is None:
        idx = list(range(len(cols)))
    else:
        pos = {c: i for i, c in enumerate(cols)}
        bad = [n for n in names if n not in pos]
        if bad:
            raise ValueError(f"{', '.join(bad)} not in {covariance_path} (free level energies and widths only: "
                             f"{' '.join(cols)})")
        idx = [pos[n] for n in names]
    sel = [cols[i] for i in idx]
    mean = np.array([vals[c] for c in sel])
    cov = C[np.ix_(idx, idx)]
    comment = (f"prior from {os.path.abspath(sav_path)}\n"
               f"and {os.path.abspath(covariance_path)}" + (f" (scaled by {scale:g})" if scale != 1.0 else ""))
    return write_prior(out_path, sel, mean, cov, basis=basis, comment=comment)
