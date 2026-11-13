#!/usr/bin/env python3
"""Check the engine's thin-target integration of one data set against an independent integral.

Reads a format-3 project, the name of a data set of ``quantity sigma`` that names a ``target`` of ``kind layer`` with
``stopping evaluation`` (SRIM through the stopping library), ``straggling none`` and no beam or detector, and,
optionally, a parameter file. Computes, outside the engine's kernel,

    Y(E_0) = (1/N) int_{E_0 - Delta}^{E_0} sigma(E') / eps(E') dE',   int_{E_0-Delta}^{E_0} dE'/eps = N,

with sigma the bare model on a fine lab grid (an extrapolation record of the same reaction in a copy of the project,
folded by nothing), eps and Delta from ``pyazr.stopping`` for the target's composition (per atom of the material, N
of the material = N of the active nucleus / its element's atom fraction), the trapezoidal rule on the grid (1/N as
1 / int dE'/eps on the same nodes), and compares Y
with the engine's folded model of the data set at every point. A difference above ``--rtol`` (default 2e-3) points
at a modelling error of the record (wrong N, wrong composition) or at the engine; equal numbers show that the record
does what was meant.

Prints per point E_lab, the engine's value, the reference and their ratio; writes nothing but a temporary copy.
Exit status: 0 when every point agrees within --rtol, 1 when one does not, 2 on a usage error (a data set the check
does not cover).
"""
import argparse
import re
import shutil
import sys
from pathlib import Path

import numpy as np


def target_fields(project, name):
    rec = project[name]
    return {k: str(rec[k]) for k in ("kind", "active", "composition", "areal_density", "stopping", "straggling")
            if rec.get(k) is not None}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--params", help="parameter file (param.sav), by name")
    ap.add_argument("--points", type=int, default=4001, help="points of the fine bare grid (default 4001)")
    ap.add_argument("--rtol", type=float, default=2e-3)
    a = ap.parse_args(argv)
    import pyazr
    from pyazr import stopping
    from pyazr.masses import parse_nuclide
    azr = a.azr.resolve()
    p = pyazr.Project.open(azr)
    ds = p[a.dataset]
    if str(ds["quantity"]) != "sigma" or ds.get("target") is None or ds.get("beam") or ds.get("detector"):
        print("error: the check covers a sigma data set on a target, without beam or detector", file=sys.stderr)
        return 2
    tg = target_fields(p, str(ds["target"]).split("@")[0])
    if tg.get("kind") != "layer" or not tg.get("stopping", "").startswith("evaluation") \
            or tg.get("straggling", "none") != "none":
        print(f"error: the target must be a layer with 'stopping evaluation' and 'straggling none' ({tg})",
              file=sys.stderr)
        return 2
    formula = "".join(f"{s}{n}" for s, n in re.findall(r"([A-Z][a-z]?)\s+(\S+)", tg["composition"]))
    material = stopping.Material(formula)
    fraction = material.atom_fraction(parse_nuclide(tg["active"])[0])
    n_material = float(tg["areal_density"].split()[0]) / fraction
    reaction = str(ds["reaction"])
    projectile = re.match(r"\s*\S+?\((\w+),", reaction).group(1)
    proj = {"p": "p", "d": "d", "t": "t", "a": "a", "3He": "3He"}.get(projectile, projectile)

    with pyazr.Session(str(azr)) as s:
        x = s.parameters.read(a.params) if a.params else s.parameters.x0
        d = next(d for d in s.evaluate(x).data if d.name == a.dataset)
        e_lab = np.asarray(d.extra["energy_lab"], float)
        folded = np.asarray(d.model, float)
    e_out = np.asarray(stopping.energy_after(material, proj, e_lab, np.full(e_lab.size, n_material)), float)
    lo, hi = float(e_out.min()) * 0.999, float(e_lab.max()) * 1.0001
    grid = np.linspace(lo, hi, a.points)
    q = pyazr.Project.open(azr)
    q.add("extrapolation", "target_check_bare", reaction=reaction, quantity="sigma", unit="b",
          energies="list " + " ".join(repr(float(e)) for e in grid) + " MeV lab")
    out = azr.parent / "output_target_check"  # its own output directory: no capture cache of the project is touched
    out.mkdir(exist_ok=True)
    q.set_config(output_dir=f"{out.name}/")
    copy = q.save(azr.with_name(f"{azr.stem}_target_check.azr"))
    try:
        with pyazr.Session(str(copy), data_mode=False) as s:
            xb = s.parameters.read(a.params) if a.params else s.parameters.x0
            sigma = np.asarray(next(d for d in s.evaluate(xb).data if d.name == "target_check_bare").model, float)
    finally:
        Path(copy).unlink(missing_ok=True)
        shutil.rmtree(out, ignore_errors=True)
    inv_eps = 1e21 / np.asarray(stopping.stopping(material, proj, grid)["total"], float)  # cm^2/MeV per atom

    status = 0
    print(f"# {a.dataset}: {reaction} on {formula}, N(material) = {n_material:.6g} atoms/cm2")
    print(f"# {'E_lab_MeV':>12s} {'engine':>14s} {'reference':>14s} {'ratio':>10s}")
    for e0, e1, y in zip(e_lab, e_out, folded):
        m = (grid >= e1) & (grid <= e0)
        xs = np.concatenate([[e1], grid[m], [e0]])
        f = np.interp(xs, grid, sigma * inv_eps)
        w = np.interp(xs, grid, inv_eps)
        ref = np.trapezoid(f, xs) / np.trapezoid(w, xs)  # the depth average: int sigma/eps over int 1/eps (= N)
        ratio = y / ref if ref else float("nan")
        if not abs(ratio - 1) <= a.rtol:
            status = 1
        print(f"  {e0:12.6f} {y:14.6e} {ref:14.6e} {ratio:10.6f}")
    print("agree" if status == 0 else f"DIFFER beyond rtol {a.rtol:g}")
    return status


if __name__ == "__main__":
    sys.exit(main())
