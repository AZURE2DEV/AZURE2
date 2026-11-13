#!/usr/bin/env python3
"""A target layer from SRIM's stopping (AZURE2's stopping library): energy loss, straggling, and the record to write.

Reads nothing but its arguments. For a material (chemical formula), a projectile and a lab beam energy E0 (MeV) and
either an areal density of the material (atoms/cm^2, all atoms) or an energy loss quoted by a paper (keV), prints:

- the stopping cross section eps(E0) per atom of the material and per atom of the active element (eV/1e15 atoms/cm^2;
  the engine's ``stopping evaluation`` uses the latter, eps / atom fraction);
- the energy loss Delta = E0 - E_out from the range integral int_{E_out}^{E0} dE/eps = N (what the engine integrates),
  and the first-order eps(E0) N beside it;
- the straggling at the exit (TRIM's, and Bohr's), 1 sigma in keV;
- the fields of a format-3 ``target`` record of ``kind layer`` with this stopping: ``active``, ``composition``,
  ``areal_density`` of the active nucleus (atoms/cm^2 = N x atom fraction of its element), ``stopping evaluation``,
  ``phase``, ``straggling``.

Every use of these stopping powers cites SRIM (printed). Energies: MeV, lab, of the projectile.
Exit status: 0, 2 on a usage error or an energy outside SRIM's range.
"""
import argparse
import re
import sys

import numpy as np


def layer(formula, projectile, e0, areal_density=None, loss_kev=None, active=None, phase=None):
    """The numbers of the layer, as a dict (see the module's docstring)."""
    from pyazr import stopping
    from pyazr.masses import parse_nuclide
    m = stopping.Material(formula, phase=phase)
    if areal_density is None:
        areal_density = float(stopping.areal_density_between(m, projectile, e0, e0 - loss_kev * 1e-3))
    eps = float(stopping.stopping(m, projectile, e0)["total"])
    out = float(stopping.energy_after(m, projectile, e0, areal_density))
    trim = float(stopping.straggling(m, projectile, e0, areal_density)["straggling"])
    bohr = float(stopping.straggling(m, projectile, e0, areal_density, model="bohr")["straggling"])
    counts = re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", formula)
    fraction = m.atom_fraction(parse_nuclide(active)[0]) if active else None
    return dict(material=m, eps=eps, eps_active=eps / fraction if fraction else None, fraction=fraction,
                areal_density=areal_density, loss=e0 - out, loss_first_order=eps * areal_density * 1e-21,
                straggling_trim=trim, straggling_bohr=bohr,
                composition=" ".join(f"{s} {n or 1}" for s, n in counts), citation=stopping.CITATION)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--material", required=True, help='chemical formula: "Ta2O5", "TiN", "C", "CH4"')
    ap.add_argument("--projectile", required=True, help='"p", "d", "t", "3He", "a"')
    ap.add_argument("--energy", type=float, required=True, help="beam energy E0 at the layer's entrance, MeV lab")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--areal-density", type=float, help="atoms/cm^2 of the material (all atoms)")
    g.add_argument("--loss-kev", type=float, help="energy loss in the layer at E0, keV (lab)")
    ap.add_argument("--active", help='the active nucleus, "15N" (for the record and the per-active-atom stopping)')
    ap.add_argument("--phase", choices=["solid", "gas"], help="SRIM's phase (default: SRIM's choice)")
    a = ap.parse_args(argv)
    try:
        r = layer(a.material, a.projectile, a.energy, a.areal_density, a.loss_kev, a.active, a.phase)
    except (ValueError, RuntimeError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    print(f"{r['material']}, {a.projectile} at {a.energy:g} MeV (lab)")
    print(f"stopping eps(E0) = {r['eps']:.5g} eV/(1e15 atoms/cm2) per atom of {a.material}"
          + (f"; {r['eps_active']:.5g} per atom of the active element (fraction {r['fraction']:.4g})"
             if r["eps_active"] else ""))
    print(f"areal density N = {r['areal_density']:.6g} atoms/cm2 of {a.material}")
    print(f"energy loss (range integral) {r['loss'] * 1e3:.4f} keV; first order eps(E0) N "
          f"{r['loss_first_order'] * 1e3:.4f}"
          f" keV ({(r['loss_first_order'] / r['loss'] - 1) * 100:+.2f} %)")
    print(f"straggling at the exit, 1 sigma: TRIM {r['straggling_trim'] * 1e3:.4f} keV, Bohr "
          f"{r['straggling_bohr'] * 1e3:.4f} keV")
    if a.active and r["fraction"]:
        print("\nformat-3 target record (fields; the areal density is of the active nucleus):")
        print(f"    kind           layer\n    active         {a.active}\n    composition    {r['composition']}")
        print(f"    areal_density  {r['areal_density'] * r['fraction']:.6g} atoms/cm2 fixed")
        print("    stopping       evaluation" + (f"\n    phase          {a.phase}" if a.phase else ""))
        print("    straggling     bohr")
    print(f"\nStopping: {r['citation']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
