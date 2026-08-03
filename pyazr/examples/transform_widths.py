"""Turn reduced-width amplitudes into physical partial widths, outside AZURE2.

AZURE2 fits reduced-width amplitudes but reports partial widths, and the step
between them is not ``Gamma = 2 gamma^2 P``: it carries the level-shift
normalisation ``1 + sum_c gamma_c^2 dS_c/dE``.  Anyone reading an MCMC chain, a
covariance sample, or a hand-made parameter vector has to redo that step, and
for a strongly coupled level the denominator is large -- 34 for the 3/2+ state
in d+t, so the naive formula overstates the width 34-fold.

Print the width table for the parameters in ``output/param.par``::

    python transform_widths.py 3H+d.azr --brune

Map a whole MCMC chain onto partial widths, one level at a time::

    python transform_widths.py 3H+d.azr --brune \
        --chain output/samples.mcmc --level 1 \
        --energy-col param0 --gamma-cols param1,param2 --channels 2,4

``--gamma-cols`` names the chain columns holding the reduced widths and
``--channels`` says which channels of the level they are (1-based, in the order
the .azr lists them); channels not named keep their .azr value, normally zero.

Warning signs the script prints for you
---------------------------------------
* ``1 + sum g^2 dS/dE`` far from 1 -- the naive formula is wrong by that factor.
* a reduced width whose partial width sits near ``2 P_c / (dS_c/dE)``.  That is
  the ceiling the transformation saturates at, and a chain pressed against it is
  telling you the data do not constrain the coupling strength: gamma can grow
  without bound while Gamma barely moves.
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

from pyazr.transform import (levels_from_azr, transform_out,   # noqa: E402
                             LevelWidthTransformer)


def read_param_par(path):
    """``{level_index: [gamma per channel]}`` from an AZURE2 parameter file."""
    per_level = {}
    with open(path) as fh:
        for line in fh:
            tok = line.split()
            if len(tok) < 2:
                continue
            m = re.match(r"width_(\d+)_(\d+)$", tok[0])
            if m:
                per_level.setdefault(int(m.group(1)), {})[int(m.group(2))] = float(tok[1])
    return {k: [v[j] for j in sorted(v)] for k, v in per_level.items()}


def saturation_limit(channel, energy):
    """``2 P_c / (dS_c/dE)`` -- the largest width this channel can reach, eV."""
    d = channel.shift_derivative(energy)
    if d <= 0:
        return None
    return 2e6 * channel.penetrability(energy) / d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", help="the .azr file the parameters belong to")
    ap.add_argument("--params", help="parameter file (default <dir>/output/param.par)")
    ap.add_argument("--brune", action="store_true",
                    help="the fit used --use-brune / the GUI's Brune option")
    ap.add_argument("--chain", help="MCMC samples (CSV with a header row)")
    ap.add_argument("--level", type=int, default=1,
                    help="1-based level index for --chain (default 1)")
    ap.add_argument("--energy-col", default="param0")
    ap.add_argument("--gamma-cols", default="",
                    help="comma-separated chain columns holding reduced widths")
    ap.add_argument("--channels", default="",
                    help="comma-separated 1-based channel numbers for those columns")
    ap.add_argument("--burn-in", type=int, default=0, help="steps to discard")
    args = ap.parse_args()

    par = args.params or os.path.join(os.path.dirname(args.azr) or ".",
                                      "output", "param.par")
    gammas = read_param_par(par) if os.path.exists(par) else {}
    if not gammas:
        print(f"no widths found in {par}; using the .azr values (zeros)",
              file=sys.stderr)

    levels = levels_from_azr(args.azr)
    for i, lv in enumerate(levels, 1):
        if i in gammas:
            if len(gammas[i]) != len(lv.channels):
                sys.exit(f"level {i}: {par} has {len(gammas[i])} widths but the "
                         f".azr has {len(lv.channels)} channels -- are they the "
                         f"same model?")
            for c, g in zip(lv.channels, gammas[i]):
                c.gamma = g

    if not args.chain:
        for i, out in enumerate(transform_out(levels, brune=args.brune), 1):
            print(f"\n[level {i}]")
            print(out.table())
            if out.normalization > 1.05:
                print(f"  note: the level-shift normalisation is "
                      f"{out.normalization:.4g}; Gamma = 2 gamma^2 P would be "
                      f"{out.normalization:.4g}x too large here")
            for tc in out.channels:
                lim = saturation_limit(tc.channel, out.energy)
                w = tc.width_eV
                if lim and w and w > 0.9 * lim:
                    print(f"  note: {tc.label} is at {100 * w / lim:.1f}% of its "
                          f"saturation ceiling 2P/(dS/dE) = {lim:.6g} eV -- the "
                          f"width barely responds to gamma above this point")
        return

    # -- chain mode ---------------------------------------------------------
    import numpy as np
    import pandas as pd

    if not args.brune:
        sys.exit("--chain needs --brune: without it the levels of a J-group mix "
                 "and the transformation is not a per-sample closed form.")

    df = pd.read_csv(args.chain)
    if args.burn_in and "step" in df:
        df = df[df["step"] >= args.burn_in]
    level = levels[args.level - 1]

    cols = [c.strip() for c in args.gamma_cols.split(",") if c.strip()]
    chans = [int(c) for c in args.channels.split(",") if c.strip()]
    if len(cols) != len(chans):
        sys.exit("--gamma-cols and --channels must have the same length")

    energy = df[args.energy_col].to_numpy(float)
    g = np.tile([c.gamma for c in level.channels], (len(df), 1))
    for col, ch in zip(cols, chans):
        g[:, ch - 1] = df[col].to_numpy(float)

    tr = LevelWidthTransformer(level, energies=energy)
    widths = tr.widths(energy, g, unit="keV")
    norm = tr.normalization(energy, g)

    print(f"{len(df)} samples, level {args.level} "
          f"(E in [{energy.min():.6g}, {energy.max():.6g}] MeV)\n")
    print(f"{'channel':<22} {'16th':>12} {'median':>12} {'84th':>12}   unit")
    for k, c in enumerate(level.channels):
        col = widths[:, k]
        if not np.isfinite(col).any() or np.allclose(col[np.isfinite(col)], 0):
            continue
        q = np.percentile(col[np.isfinite(col)], [16, 50, 84])
        lab = f"pair{c.pair} l={c.L}" + ("" if c.S is None else f" s={c.S:g}")
        print(f"{lab:<22} {q[0]:12.6g} {q[1]:12.6g} {q[2]:12.6g}   keV")

    print(f"\n1 + sum g^2 dS/dE : median {np.median(norm):.6g}, "
          f"range [{norm.min():.6g}, {norm.max():.6g}]")
    if np.median(norm) > 1.05:
        print("  the naive Gamma = 2 gamma^2 P would be wrong by that factor")

    for k, c in enumerate(level.channels):
        lim = saturation_limit(c, float(np.median(energy)))
        col = widths[:, k]
        if lim and np.isfinite(col).any():
            frac = np.nanmedian(col) * 1e3 / lim
            if frac > 0.9:
                print(f"  pair{c.pair} l={c.L}: median width is {100 * frac:.1f}% "
                      f"of the ceiling 2P/(dS/dE) = {lim * 1e-3:.6g} keV -- the "
                      f"chain is on the flat part of the transformation, so its "
                      f"width uncertainty is set by the prior, not the data")


if __name__ == "__main__":
    main()
