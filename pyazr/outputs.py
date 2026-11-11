"""Readers of the engine's output files."""
import numpy as np

__all__ = ["fitted_norms"]


def fitted_norms(path="output/chiSquared.out", nsegments=None):
    """Per-segment normalizations from a fit, indexed by 0-based segment.

    ``chiSquared.out`` is one comma-separated line per segment --
    ``key, chi2, N, norm,`` -- followed by a totals line that is *not* CSV, so
    it has to be parsed by segment key rather than by row position.  (Reading
    it with ``skiprows=2`` drops segment 1 and shifts every normalization onto
    the wrong segment.)

    AZURE2's residual is ``(fit - data*n)/(err*n)``: the normalization scales
    the *data*, and its own ``AZUREOut_*.out`` files report the scaled values.
    Multiply measured points by ``fitted_norms()[i]`` to put them on the same
    footing as the model curve.
    """
    norms = {}
    with open(path) as f:
        for line in f:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 4:
                continue                       # header, blank, or totals line
            try:
                key = int(parts[0])
                norms[key] = float(parts[3])
            except ValueError:
                continue
    if not norms:
        raise ValueError(f"{path}: no per-segment lines found.")
    n = nsegments if nsegments is not None else max(norms)
    missing = [k for k in range(1, n + 1) if k not in norms]
    if missing:
        raise ValueError(f"{path} has no normalization for segment(s) "
                         f"{missing[:5]}{'...' if len(missing) > 5 else ''}.")
    return np.array([norms[k] for k in range(1, n + 1)])
