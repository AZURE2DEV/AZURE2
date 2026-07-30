"""Save a fit back into a .azr you can reopen in the GUI -- and prove it round-trips.

A fit result lives in a numpy vector or an ``output/param.sav``; neither is
something you can open, plot or hand to a colleague.  ``AzrModel.apply_fit``
turns one into a real ``.azr`` snapshot.  Two things make that less trivial than
it looks, and this example handles both:

1. **The ``<levels>`` ``gamma`` field is not a reduced-width amplitude.**  It
   holds the physical value AZURE2 prints in ``parameters.out``: a partial width
   in **eV** for an open particle channel, an **ANC in fm^-1/2** for a closed
   (sub-threshold) one, and a partial width in eV for a photon channel.  Those
   differ from the rwa by factors of 10^2 to 10^7, and a file written with the
   rwa loads without complaint and is wrong.  So the conversion must be
   explicit: pass ``transform=m.transform_rwa`` (or convert yourself and pass
   ``physical=True``).  ``apply_fit`` raises if given neither.

2. **``apply_fit`` writes ``<levels>`` only.**  Normalizations are not in that
   block, so a fit that moved them is only half saved by the ``.azr``.  This
   writes a companion ``param.sav`` carrying every parameter, which is what you
   hand to AZURE2 as the external parameter file.

Run it on a saved fit:

    python save_fit_to_azr.py 7Be.azr --params output/param.sav --out 7Be_fit.azr

and it will also verify the result: reload the written file, transform its own
parameters back to physical values, and check every R-matrix entry against the
fit.  A file that fails that check is not a snapshot of anything.

You can snapshot a *different* model at the same time -- e.g. the same fit at
another channel radius, as a starting point for a refit:

    python save_fit_to_azr.py 7Be.azr --params output/param.sav \\
        --radius 1=4.70 --radius 2=4.40 --radius 3=4.40 --out 7Be_r4.7.azr

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


def write_param_sav(model, x, path, errors_from=None):
    """Write every parameter (free and fixed) in ``param.sav`` format.

    This is the half ``apply_fit`` cannot cover: normalizations and energy
    shifts are not in ``<levels>``.  Hand this file to AZURE2 as the external
    parameter file and the model is fully restored.  Minuit errors are copied
    from ``errors_from`` (the original ``param.sav``) where available.
    """
    errs = {}
    if errors_from and os.path.exists(errors_from):
        for line in open(errors_from):
            f = line.split()
            if len(f) >= 3:
                errs[f[0]] = float(f[2])
    free = {p.name: p.free_index for p in model.parameters if not p.fixed}
    with open(path, "w") as fh:
        for p in model.parameters:
            v = float(x[free[p.name]]) if p.name in free else float(p.value)
            fh.write(f"{p.name:>28s} {v: .7e} {errs.get(p.name, 0.0): .7e}\n")


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

    with azure2(src, nprocs=1, cwd=here) as m:
        x = load_params(m, args.params)
        chi2 = float(np.sum(m.calculate_chi2_rwa(x)))
        print(f"loaded {len(x)} free parameters, chi2(data) = {chi2:.4f}")

        # >>> the point of the example <<<
        mdl.apply_fit(m.parameters, x, transform=m.transform_rwa)
        out = mdl.write(os.path.abspath(args.out))
        print(f"wrote {out}  ({mdl.applied} level values)")

        sav = os.path.splitext(out)[0] + ".sav"
        write_param_sav(m, x, sav, errors_from=args.params
                        if args.params.endswith(".sav") else None)
        print(f"wrote {sav}  (normalizations too -- apply_fit covers <levels> only)")

        physical = np.asarray(m.transform_rwa(x), float)

    # Verify: a file that does not reproduce the fit is not a snapshot of it.
    with azure2(out, nprocs=1, cwd=here) as m2:
        free = sorted((p for p in m2.parameters if not p.fixed),
                      key=lambda p: p.free_index)
        back = np.asarray(m2.transform_rwa(np.asarray(m2.params_rwa, float)), float)
        n = min(len(physical), len(back))
        bad = [(free[i].name, physical[i], back[i]) for i in range(n)
               if free[i].kind in ("energy", "width")
               and not np.isclose(physical[i], back[i], rtol=2e-5, atol=1e-9)]
        print(f"\nround-trip: {n - len(bad)}/{n} R-matrix parameters match the fit")
        for name, a, b in bad[:5]:
            print(f"   MISMATCH {name}: fit {a:.6g}  file {b:.6g}")

        # the .azr alone leaves norms at their file defaults, so patch them from
        # the companion .sav to reproduce the fitted chi-squared
        vals = {}
        for line in open(sav):
            f = line.split()
            if len(f) >= 2:
                vals[f[0]] = float(f[1])
        xv = np.array([vals[p.name] for p in free], float)
        print(f"chi2(data) from the written file = "
              f"{float(np.sum(m2.calculate_chi2_rwa(xv))):.4f}")
        if bad:
            raise SystemExit("the written file does not reproduce the fit.")


if __name__ == "__main__":
    main()
