#!/usr/bin/env python3
"""Read a project's level scheme and its data sets: what the engine made of the file, before any calculation.

``Session.level_scheme`` is the engine's view of ``<levels>`` (pairs, J-groups, levels, channels, physical widths in
eV and the Wigner limit), ``Session.datasets.table()`` the data sets (file, reaction, quantity, rows kept), and
``level_scheme.resonances(energy_max=)`` a selection by excitation energy.

Usage (from the repository's root):

    python3 pyazr/examples/print_scheme.py tests/upgrade/T1/T1.azr [--emax 5]
"""
import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "2")

from pyazr import Session  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr)")
    ap.add_argument("--emax", type=float, default=5.0, help="list the levels below this excitation energy (MeV)")
    args = ap.parse_args()

    with Session(args.project) as s:
        scheme = s.level_scheme
        print(scheme)                   # pairs, J-groups, levels, channels, widths, Wigner limits
        print()
        print(s.datasets.table())       # data file, reaction, quantity, rows kept of rows read

        print(f"\nlevels below E_x = {args.emax:g} MeV:")
        for lv in scheme.resonances(energy_max=args.emax):
            print(f"  {lv.jpi:>4}  E_x = {lv.energy:.4f} MeV  ({len(lv.channels)} channels, {lv.n_free} free)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
