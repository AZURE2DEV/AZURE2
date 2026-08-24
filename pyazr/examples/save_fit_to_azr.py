"""Save a fit back into a .azr you can reopen in the GUI.

``azure2.save_fit`` does the whole thing: it converts the fit to the physical
values a ``<levels>`` line holds, writes a companion ``param.sav`` for the
normalizations and energy shifts that block cannot carry, and verifies the
result before handing it back.  This example is the wrapper around it -- load a
fit from ``param.sav`` or a numpy file, optionally re-radius the model, and
write the snapshot:

    python save_fit_to_azr.py 7Be.azr --params output/param.sav --out 7Be_fit.azr

Three things ``save_fit`` handles that used to be the caller's problem, and
which are worth knowing about even now:

1. **The ``<levels>`` ``gamma`` field is not a reduced-width amplitude.**  It
   holds the physical value AZURE2 prints in ``parameters.out``: a partial width
   in **eV** for an open particle channel, an **ANC in fm^-1/2** for a closed
   (sub-threshold) one.  Those differ from the rwa by factors of 10^2 to 10^7,
   and a file written with the rwa loads without complaint and is wrong.

2. **A ``Parameter``'s ``pair`` is not the file's pair key.**  It is the
   engine's number, which counts pairs in the order ``<levels>`` first mentions
   them.  On the 8Be model engine pair 1 is file key 2, so matching one against
   the other put every width on the wrong channel.

3. **The snapshot is not a fit until it reads back as one.**  ``save_fit``
   reopens what it wrote and compares; if anything disagrees it removes the
   files and raises, rather than leaving a plausible-looking blend of two fits.

You can snapshot a *different* model at the same time -- e.g. the same fit at
another channel radius, as a starting point for a refit:

    python save_fit_to_azr.py 7Be.azr --params output/param.sav \\
        --radius 1=4.70 --radius 2=4.40 --out 7Be_r4.7.azr

Remember that the reduced widths then mean something different, so that file is
a *start*, not a fit -- refit before believing anything in it.
"""

import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np                                            # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

from pyazr import azure2, AzrModel                            # noqa: E402


def load_params(model, path):
    """A free-parameter vector from ``param.sav`` / ``param.par``, or a numpy file."""
    if path.endswith((".sav", ".par")):
        sav = {}
        for line in open(path):
            f = line.split()
            if len(f) >= 2:
                sav[f[0]] = float(f[1])
        free = sorted((p for p in model.parameters if not p.fixed),
                      key=lambda p: p.free_index)
        missing = [p.name for p in free if p.name not in sav]
        if missing:
            raise SystemExit(
                f"{path} has no entry for {len(missing)} free parameter(s): "
                f"{missing[:5]}{' ...' if len(missing) > 5 else ''}\n"
                "It was written for a different .azr.")
        extra = sorted(set(sav) - {p.name for p in model.parameters})
        if extra:
            print(f"note: {path} carries {len(extra)} parameter(s) this .azr does "
                  f"not have ({', '.join(extra[:4])}"
                  f"{' ...' if len(extra) > 4 else ''}); they are ignored.",
                  file=sys.stderr)
        return np.array([sav[p.name] for p in free], float)
    src = np.load(path)
    return np.asarray(src["x"] if hasattr(src, "files") else src, float)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("azr", help="the .azr the fit was made against")
    ap.add_argument("--params", required=True, help="output/param.sav or an .npz/.npy")
    ap.add_argument("--out", required=True, help="the .azr to write")
    ap.add_argument("--radius", action="append", default=[], metavar="PAIR=FM",
                    help="also change a channel radius, e.g. --radius 1=4.70 "
                         "(repeatable; the result is a starting point, not a fit)")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(args.azr)) or "."
    mdl = AzrModel.from_file(os.path.abspath(args.azr))

    for spec in args.radius:
        pair, _, fm = spec.partition("=")
        mdl.set_channel_radius(int(pair), float(fm))
    if args.radius:
        print(f"channel radii: {mdl.channel_radii()}")
        # AZURE2 caches its external-capture integrals in output/intEC*, and
        # those belong to the radii they were built with -- it will silently
        # reuse them on the new ones.  They are only a cache; drop them.
        for cache in ("intEC.dat", "intEC.extrap"):
            f = os.path.join(here, "output", cache)
            if os.path.exists(f):
                os.remove(f)
                print(f"removed output/{cache} (belongs to the old radii)")

    # The instance must be built from the model we are about to write into --
    # same radii, same channels -- or the parameter metadata will not line up.
    src = mdl.to_tempfile() if args.radius else os.path.abspath(args.azr)

    with azure2(src, cwd=here) as m:
        x = load_params(m, args.params)
        chi2 = float(np.sum(m.calculate_chi2_rwa(x)))
        print(f"loaded {len(x)} free parameters, chi2(data) = {chi2:.4f}")

        # >>> the point of the example <<<
        # Applies the transform, writes the companion param.sav, and verifies
        # the result -- raising, and removing both files, if it does not read
        # back as the fit it was made from.
        out, sav = m.save_fit(os.path.abspath(args.out), x)
        print(f"wrote {out}")
        print(f"wrote {sav}  (normalizations too -- <levels> cannot hold them)")

    # The .azr alone leaves the normalizations at their file defaults, so a
    # calculate run on it sits above the fit by whatever they were absorbing.
    # Feed the companion .sav back to reproduce the fitted chi-squared.
    with azure2(out, cwd=here) as m2:
        free = sorted((p for p in m2.parameters if not p.fixed),
                      key=lambda p: p.free_index)
        vals = {}
        for line in open(sav):
            f = line.split()
            if len(f) >= 2:
                vals[f[0]] = float(f[1])
        xv = np.array([vals[p.name] for p in free], float)
        print(f"chi2(data) from the written file = "
              f"{float(np.sum(m2.calculate_chi2_rwa(xv))):.4f}")


if __name__ == "__main__":
    main()
