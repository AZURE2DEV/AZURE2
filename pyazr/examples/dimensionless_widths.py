"""Dimensionless widths of a fit: theta^2 for particle channels, W.u. for gammas.

A partial width in eV is hard to judge on its own.  This example turns every
fitted width of a model into the number you can actually argue about:

  * particle channels -> theta^2 = Gamma_c / Gamma_W (Gamma in eV), with the Wigner width
    Gamma_W = 2 P_l(E_r) gamma_W^2 that the AZURE2 GUI shows next to the
    channel.  theta^2 <= 1 is the single-particle sum-rule expectation; a fit
    that needs theta^2 = 3 is telling you something is wrong.
  * photon channels -> the strength in Weisskopf units.  A few W.u. is normal
    for E1/E2; hundreds means the level is soaking up strength it should not.

Run it on the project as it stands, or at a fit's parameter file (``param.sav``, matched to the project by the
parameters' names; the session takes it when it is made and warns about a name it does not have):

Usage (from the repository's root):

    python3 pyazr/examples/dimensionless_widths.py tests/upgrade/T1/T1.azr [--params output/param.sav] [--all]
"""

import argparse
import os

os.environ.setdefault("OMP_NUM_THREADS", "2")

from pyazr import Session                                      # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("azr", help="the format-3 project (.azr)")
    ap.add_argument("--params", help="a parameter file of a fit (param.sav, param.par)")
    ap.add_argument("--all", action="store_true",
                    help="also list the inert (zero-width) channels")
    args = ap.parse_args()

    with Session(args.azr, parameter_file=args.params) as m:
        x = m.parameters.x0
        print(f"compound A = {m.mass_number}, "
              f"chi2 = {m.evaluate(x).chi2:.2f}\n")

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
