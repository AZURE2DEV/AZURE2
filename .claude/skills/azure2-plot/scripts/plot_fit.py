#!/usr/bin/env python3
"""Figures of an AZURE2 fit: every data set with the model and the pulls, drawn as the engine compares them.

Reads a format-3 project and, optionally, a parameter file (by name); evaluates once and writes one figure per data
set (``<outdir>/<data set>.png``, or ``.pdf`` with ``--format pdf``): top panel the datum times the data set's
normalization n with its error times n, and the model, against the c.m. energy (or the c.m. angle when the data set has
one energy); bottom panel the pulls (model - datum n)/(error n) of the engine. With ``--sfactor`` cross sections are
drawn as S factors (the engine's conversion, MeV b; never a hand-made Sommerfeld factor). The model is the engine's
value at the data points (folded as the data set is); a smooth bare curve on a fine grid comes from azure2-calculate's
``extrapolate.py``.

It prints, per data set, the largest difference between the pulls it drew and the engine's residuals (a figure that
does not show what the fit compares is wrong), and writes ``<outdir>/figures.txt`` with the data chi2 and the
parameter file used.

Style is plain matplotlib; for publication style use the user-level prc-figures skill on the same arrays.
Exit status: 0, 1 when a drawn pull differs from the engine's residual by more than 1e-6, 2 on a usage error.
"""
import argparse
import sys
from pathlib import Path

import numpy as np


def option_value(text):
    low = text.lower()
    if low in ("true", "yes", "on", "1"):
        return True
    if low in ("false", "no", "off", "0"):
        return False
    try:
        return float(text)
    except ValueError:
        return text


def arrays(d, sfactor):
    """(x, x label, model, datum*n, error*n, pulls, y label) of one evaluated data set."""
    e, th = np.asarray(d.energies, float), np.asarray(d.angles, float)
    by_angle = np.ptp(e) == 0 and np.ptp(th) > 0
    x, xlab = (th, r"$\theta_{\rm c.m.}$ (deg)") if by_angle else (e, r"$E_{\rm c.m.}$ (MeV)")
    n = float(d.norm)
    model, y, dy = np.asarray(d.model, float), np.asarray(d.values, float) * n, np.asarray(d.errors, float) * n
    pulls = (model - y) / dy
    ylab = f"{d.quantity} ({d.unit})"
    conv = d.extra.get("sfactor_conversion")
    if sfactor and conv is not None and d.quantity in ("sigma", "sigma_total", "sigma_nuclear"):
        c = np.asarray(conv, float)
        model, y, dy, ylab = model * c, y * c, dy * c, "S (MeV b)"
    return x, xlab, model, y, dy, pulls, ylab


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr")
    ap.add_argument("--params", help="parameter file (param.sav), by name; default: the project's values")
    ap.add_argument("--option", action="append", default=[], metavar="KEY=VALUE", help="a session option")
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--format", default="png", choices=["png", "pdf"])
    ap.add_argument("--sfactor", action="store_true", help="draw cross sections as S factors")
    ap.add_argument("--dataset", action="append", help="only these data sets (record names)")
    a = ap.parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pyazr
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    opts = {k: option_value(v) for k, v in (o.split("=", 1) for o in a.option)}
    status = 0
    with pyazr.Session(a.azr, **opts) as s:
        x = s.parameters.read(a.params) if a.params else s.parameters.x0
        ev = s.evaluate(x)
        for d in ev.data:
            if a.dataset and d.name not in a.dataset:
                continue
            xs, xlab, model, y, dy, pulls, ylab = arrays(d, a.sfactor)
            diff = float(np.max(np.abs(pulls - np.asarray(d.residuals, float)))) if d.n else 0.0
            print(f"{d.name}: {d.n} points, chi2 {d.chi2:.4f}, norm {d.norm:.5g}, drawn pulls vs engine residuals "
                  f"{diff:.1e}")
            if diff > 1e-6:
                status = 1
            order = np.argsort(xs, kind="stable")
            fig, (top, bot) = plt.subplots(2, 1, sharex=True, figsize=(6, 5), gridspec_kw=dict(height_ratios=[3, 1]))
            top.errorbar(xs, y, dy, fmt="o", ms=3, color="k", label=f"data x n (n = {d.norm:.4g})")
            top.plot(xs[order], model[order], "-", color="C3", label="model at the data points")
            top.set_ylabel(ylab)
            top.legend(frameon=False, fontsize=8)
            top.set_title(f"{d.name}: chi2 = {d.chi2:.2f}, N = {d.n}", fontsize=9)
            bot.axhline(0, color="0.5", lw=0.8)
            bot.plot(xs, pulls, "o", ms=3, color="k")
            bot.set_ylabel("pull")
            bot.set_xlabel(xlab)
            fig.tight_layout()
            fig.savefig(out / f"{d.name}.{a.format}", dpi=150)
            plt.close(fig)
        (out / "figures.txt").write_text(f"project {a.azr}\nparameters {a.params or 'the project values'}\n"
                                         f"data chi2 {ev.chi2:.10g}\n")
    return status


if __name__ == "__main__":
    sys.exit(main())
