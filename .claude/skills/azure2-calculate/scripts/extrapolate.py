#!/usr/bin/env python3
"""A cross section, S factor, differential cross section or phase shift of an AZURE2 model on any grid.

Reads a format-3 project and, optionally, a parameter file (by name). Writes a copy of the project beside it
(``<stem>_extrap.azr``, so that relative paths hold) with one ``extrapolation`` record added -- bare: it names no
target, beam or detector, so nothing folds it -- opens a session in extrapolation mode on the copy, evaluates at the
parameters and writes the record's points as text columns:

    E_cm_MeV  E_x_MeV  theta_cm_deg  value  S_MeV_b

``value`` in the record's unit (b, b/sr, deg, 1); ``S_MeV_b`` is ``model * extra["sfactor_conversion"]`` of the
engine (0 where it gives none, e.g. a phase shift). The grid is the record's ``energies`` field (format 3:
``"range 0.01 0.6 MeV cm points 120 log"``, ``"list 0.1 0.2 MeV lab"``, ``"range 0.1 2 MeV lab step 0.01"``), its frame
stated in it; the output energies are c.m. The copy is removed unless ``--keep``.

S(0): give a grid that starts well below the data (a log grid from about 1 keV) and read S at its lowest energies;
check that halving the lowest energy does not move it.

Exit status: 0, 1 when the engine refuses the record or the evaluation, 2 on a usage error.
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


def extrapolate(azr, fields, params=None, options=None, name="skill_grid", keep=False):
    """{'energies', 'excitation_energies', 'angles', 'model', 'sfactor'} of one bare extrapolation record."""
    import pyazr
    azr = Path(azr).resolve()
    p = pyazr.Project.open(azr)
    p.add("extrapolation", name, **fields)
    copy = p.save(azr.with_name(f"{azr.stem}_extrap.azr"))
    try:
        with pyazr.Session(str(copy), data_mode=False, **(options or {})) as s:
            x = s.parameters.read(params) if params else s.parameters.x0
            d = next(d for d in s.evaluate(x).data if d.name == name)
            conv = d.extra.get("sfactor_conversion")
            model = np.asarray(d.model, float)
            return dict(energies=np.asarray(d.energies, float), excitation_energies=np.asarray(
                d.excitation_energies, float), angles=np.asarray(d.angles, float), model=model,
                sfactor=model * np.asarray(conv, float) if conv is not None else np.zeros_like(model))
    finally:
        if not keep:
            Path(copy).unlink(missing_ok=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr")
    ap.add_argument("--params", help="parameter file (param.sav), by name; default: the project's values")
    ap.add_argument("--option", action="append", default=[], metavar="KEY=VALUE", help="a session option")
    ap.add_argument("--reaction", required=True, help='as in a record: "12C(p,g)13N", "12C(p,p)12C", '
                    '"12C(p,g)13N[1.0]" (capture to the 1.0 MeV state)')
    ap.add_argument("--quantity", default="sigma", help="sigma, dsigma_domega, phase_shift, analyzing_power, ...")
    ap.add_argument("--unit", default=None, help="b, b/sr, deg, 1 (default: from the quantity)")
    ap.add_argument("--energies", required=True, help='the grid: "range 0.01 0.6 MeV cm points 120 log"')
    ap.add_argument("--angles", help='for a differential quantity: "list 90 deg cm"')
    ap.add_argument("--field", action="append", default=[], metavar="KEY=VALUE",
                    help="another field of the record (J=1/2, l=0, s=1/2 for a phase shift)")
    ap.add_argument("--out", help="write the table here (default: standard output)")
    ap.add_argument("--keep", action="store_true", help="keep the copy of the project")
    a = ap.parse_args(argv)
    unit = a.unit or {"sigma": "b", "dsigma_domega": "b/sr", "phase_shift": "deg", "analyzing_power": "1"}.get(
        a.quantity)
    if unit is None:
        print("error: give --unit for this quantity", file=sys.stderr)
        return 2
    fields = dict(reaction=a.reaction, quantity=a.quantity, unit=unit, energies=a.energies)
    if a.angles:
        fields["angles"] = a.angles
    fields.update(f.split("=", 1) for f in a.field)
    opts = {k: option_value(v) for k, v in (o.split("=", 1) for o in a.option)}
    try:
        r = extrapolate(a.azr, fields, a.params, opts, keep=a.keep)
    except Exception as err:
        print(f"error: {type(err).__name__}: {err}", file=sys.stderr)
        return 1
    lines = [f"# {a.reaction} {a.quantity} [{unit}] at {a.energies}" + (f", {a.angles}" if a.angles else ""),
             f"# {'E_cm_MeV':>18s} {'E_x_MeV':>18s} {'theta_cm_deg':>12s} {'value':>20s} {'S_MeV_b':>20s}"]
    lines += [f"  {e:18.12g} {ex:18.12g} {th:12.6g} {v:20.12e} {sf:20.12e}" for e, ex, th, v, sf in
              zip(r["energies"], r["excitation_energies"], r["angles"], r["model"], r["sfactor"])]
    text = "\n".join(lines) + "\n"
    if a.out:
        Path(a.out).write_text(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
