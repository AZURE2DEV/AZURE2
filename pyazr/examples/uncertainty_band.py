"""Uncertainty bands from a saved fit covariance.

A fit run with the uncertainty band enabled leaves two files behind:

    output/param.sav        the best-fit parameters
    output/covariance.dat   their covariance (the free R-matrix block)

Given those, the 1-sigma band on any observable the model can calculate is
linear propagation,

    sigma_y^2 = g^T C g,        g_j = d y / d p_j

where g is the per-point sensitivity matrix.  AZURE2 computes g analytically --
one reverse-mode adjoint per point, so the whole matrix costs about two forward
evaluations no matter how many parameters the fit has -- and pyazr asks for it
over the API (``azure2.model_gradients``).  ``uncertainty_bands`` puts the two
together and reproduces the ``.band`` files AZURE2 writes itself to six
decimals, on whatever grid you ask for, without re-running a fit.

Three things to know:

* **Evaluate at ``param.sav``, not at the ``.azr``.**  A ``.azr`` need not carry
  the fitted parameters.  On this model the difference is chi2 = 173006 (file)
  against 4130 (fit).  ``params=`` defaults to ``param.sav`` for that reason.

* **Ask for the segments you want, and no more.**  The adjoint pass covers every
  active segment, so a band over the full 7Be ``<segmentsTest>`` block (110
  grids, ~390k points) is far more work than the handful you plot.
  ``extrapolation_bands`` trims the block in a private copy of the model.

* **The first call on a new set of grids is the slow one.**  AZURE2 builds the
  external-capture integrals for the grid at startup -- minutes on a fine
  capture grid.  They are cached per grid under ``output/bands/<hash>/``, so
  every later call on the same grids skips it.  (That directory is private to
  the trimmed model, which also keeps a stale ``intEC.extrap`` from ever being
  read back into the main run.)

Run from the directory holding the ``.azr`` (as a module, so ``pyazr`` beside
it is importable):

    python -m pyazr.examples.uncertainty_band
"""

import os
import time

os.environ.setdefault("OMP_NUM_THREADS", "4")   # set BEFORE importing numpy

import numpy as np

from pyazr import azure2, extrapolation_bands, trimmed_model, uncertainty_bands

AZR = "7Be.azr"


def main():
    # ------------------------------------------------------------------
    # 1.  The one-liner: bands on segments the .azr already declares.
    #     Keys are the 1-based <segmentsTest> numbers -- print
    #     azure2(...).extrapolations.table() to see them.
    # ------------------------------------------------------------------
    t = time.time()
    b = extrapolation_bands(AZR, keys=[113], quantity="sfactor")[113]
    print(f"{b}   [{time.time() - t:.1f} s]")
    print(f"  S(E -> 0)  = {b.value[0] * 1e3:.4f} +- {b.sigma[0] * 1e3:.4f} keV b")
    print(f"  band width = {100 * b.relative.min():.2f}% to "
          f"{100 * b.relative.max():.2f}%")

    lo, hi = b.interval(2)                       # 2-sigma envelope for a plot
    print(f"  2 sigma at threshold: [{lo[0] * 1e3:.4f}, {hi[0] * 1e3:.4f}] keV b")

    # ------------------------------------------------------------------
    # 2.  A grid the file does not contain: give the spec instead of a key.
    #     Energies are LAB, as everywhere in <segmentsTest>.
    # ------------------------------------------------------------------
    grid = dict(entrance=1, exit=-1, e_min=0.01, e_max=1.0, e_step=0.01,
                observable="total-capture")
    print(extrapolation_bands(AZR, grids=[grid], quantity="sfactor")[1])

    # ------------------------------------------------------------------
    # 3.  Data and band from a single session.  Opening two azure2 sessions in
    #     one interpreter is asking for trouble, so when a plot needs both the
    #     measured points and the band, trim the model yourself, read the data
    #     in data mode, then switch and call uncertainty_bands.  This is what
    #     plot-capture.py and plot-angular.py do.
    # ------------------------------------------------------------------
    tmp, workdir = trimmed_model(AZR, keys=[113])
    try:
        with azure2(tmp, cwd=os.path.dirname(os.path.abspath(AZR)) or ".") as m:
            capture = [i for i, s in enumerate(m.datasets)
                       if s.observable == "total-capture" and s.entrance_key == 1]
            print(f"\n{len(capture)} 3He+4He capture datasets, "
                  f"{sum(len(m.energies[i]) for i in capture)} points")

            m.extrap_mode()
            # The trimmed file renumbers the grids it kept 1..N.
            print(uncertainty_bands(m, keys=[1], quantity="sfactor")[1])
    finally:
        os.remove(tmp)

    # ------------------------------------------------------------------
    # 4.  Cross-check the analytic sensitivities against finite differences.
    #     Same answer, but 2*ncols forward passes instead of ~2 -- use it to
    #     validate, not to produce.
    # ------------------------------------------------------------------
    t = time.time()
    fd = extrapolation_bands(AZR, keys=[113], quantity="sfactor",
                             method="central")[113]
    print(f"\ncentral differences: {time.time() - t:.1f} s, "
          f"sigma agrees to {np.abs(fd.sigma / b.sigma - 1).max():.1e}")


if __name__ == "__main__":
    main()
