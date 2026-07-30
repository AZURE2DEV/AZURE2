"""Dimensionless widths of a fit: theta^2 for particle channels, W.u. for gammas.

A partial width in eV is hard to judge on its own.  This example turns every
fitted width of a model into the number you can actually argue about:

  * particle channels -> theta^2 = Gamma_c / Gamma_W, with the Wigner width
    Gamma_W = 2 P_l(E_r) gamma_W^2 that the AZURE2 GUI shows next to the
    channel.  theta^2 <= 1 is the single-particle sum-rule expectation; a fit
    that needs theta^2 = 3 is telling you something is wrong.
  * photon channels -> the strength in Weisskopf units.  A few W.u. is normal
    for E1/E2; hundreds means the level is soaking up strength it should not.

Run it on the model as it stands in the .azr:

    python dimensionless_widths.py 7Be.azr

or on a saved best fit -- either AZURE2's own ``output/param.sav`` or an
``.npz``/``.npy`` from a Python fit:

    python dimensionless_widths.py 7Be.azr --params output/param.sav

Note that ``param.sav`` is matched to the model *by parameter name*, so it must
come from a run of the same .azr; if a normalization was free then and is fixed
now (or vice versa) the script says so instead of silently misaligning.
"""

import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np                                            # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

from pyazr import azure2                                      # noqa: E402


def load_params(model, path):
    """A free-parameter vector from ``param.sav`` or a numpy file."""
    if path.endswith(".sav") or path.endswith(".par"):
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
                "It was written for a different .azr -- refit or fix the model.")
        # The reverse mismatch is quieter and nastier: the fit had a parameter
        # this .azr no longer does (a normalization since fixed, a level since
        # removed).  The widths still load, but the chi-squared will not
        # reproduce the fit, so say so rather than let it pass as agreement.
        extra = sorted(set(sav) - {p.name for p in model.parameters})
        if extra:
            print(f"warning: {path} carries {len(extra)} parameter(s) this .azr "
                  f"does not have ({', '.join(extra[:4])}"
                  f"{' ...' if len(extra) > 4 else ''}).\n"
                  f"         They are ignored -- expect the chi-squared below to "
                  f"differ from the fit's.\n", file=sys.stderr)
        return np.array([sav[p.name] for p in free], float)
    src = np.load(path)
    return np.asarray(src["x"] if hasattr(src, "files") else src, float)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("azr", help="the .azr model file")
    ap.add_argument("--params", help="output/param.sav, or an .npz/.npy fit result")
    ap.add_argument("--all", action="store_true",
                    help="also list the inert (zero-width) channels")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(args.azr)) or "."
    with azure2(os.path.abspath(args.azr), nprocs=1, cwd=here) as m:
        x = load_params(m, args.params) if args.params else np.asarray(m.params_rwa, float)
        print(f"compound A = {m.mass_number}, "
              f"chi2 = {np.sum(m.calculate_chi2_rwa(x)):.2f}\n")

        widths = m.dimensionless_widths(x)

        print("=== particle channels: theta^2 = Gamma / Gamma_W ===")
        print(widths.particles.table(nonzero_only=not args.all))

        print("\n=== photon channels: Weisskopf units ===")
        print(widths.photons.table(nonzero_only=not args.all))

        # what a referee looks at first: anything past the sum-rule limit
        over = [c for c in widths.particles if c.theta2 and c.theta2 > 1.0]
        if over:
            print(f"\n{len(over)} channel(s) above the Wigner limit:")
            for c in sorted(over, key=lambda c: -c.theta2):
                print(f"  {c.jpi:>5s} {c.level_energy:8.4f} MeV  {c.channel:<22s}"
                      f"  theta^2 = {c.theta2:.2f}")
        else:
            print("\nno particle channel exceeds the Wigner limit.")

        strong = [c for c in widths.photons if c.wu and c.wu > 10.0]
        if strong:
            print(f"\n{len(strong)} gamma channel(s) above 10 W.u.:")
            for c in sorted(strong, key=lambda c: -c.wu):
                print(f"  {c.jpi:>5s} {c.level_energy:8.4f} MeV  {c.channel:<22s}"
                      f"  {c.wu:.1f} W.u.  (E_gamma = {c.e_gamma:.3f} MeV)")

        # theta^2_formal uses the API's energy-independent limit instead, and
        # differs by the level-shift derivative -- worth seeing for broad levels
        print("\n=== where the two conventions disagree most ===")
        both = [c for c in widths.particles
                if c.theta2 and c.theta2_formal]
        for c in sorted(both, key=lambda c: -abs(c.theta2_formal / c.theta2 - 1))[:5]:
            print(f"  {c.jpi:>5s} {c.level_energy:8.4f} MeV  {c.channel:<22s}"
                  f"  theta^2 = {c.theta2:.3f}  vs  theta^2_formal = "
                  f"{c.theta2_formal:.3f}  ({c.theta2_formal/c.theta2:.2f}x)")


if __name__ == "__main__":
    main()
