#!/usr/bin/env python3
"""Wigner and Weisskopf limits of every channel of an AZURE2 project, and the level norm J of every level.

Reads a format-3 project and, optionally, a parameter file (by name); prints one row per channel from
``Session.dimensionless_widths(x)``: the level (J^pi, E_x in MeV), the pair, l, s, the physical value (Gamma in eV for
an open particle channel, the ANC in fm^-1/2 for a closed one, Gamma_gamma in eV for a photon), and

- particle channels: theta^2 = Gamma/Gamma_W (Gamma_W = 2 P_l(E_lambda) gamma^2_W, in eV) and
  theta^2_formal = gamma^2/gamma^2_W, with gamma^2_W = 3 hbar^2/(2 mu a_c^2) (Teichmann-Wigner, MeV); a fixed channel
  (not in the fit vector) gets Gamma_W from ``pyazr.physics.single_particle_width`` at the level energy, and no
  theta^2_formal;
- photon channels: the Weisskopf estimate (eV) and the strength in Weisskopf units (W.u.).

In the Brune basis the ratio theta^2/theta^2_formal of an open channel is the level's norm
J = 1/(1 + sum_c gamma_c^2 dS_c/dE) (Park, Phys. Rev. C 104, 064612 (2021)), the same for every channel of the level;
it is printed per level. J well below 0.5 marks a saturated level (its widths are not constrained by the data).
In the Park basis ``Session.park_norms(x)`` gives J directly; it is printed too.

Flags: theta^2 > --max-theta2 (default 1) on an open particle channel, a photon above --max-wu W.u. (default 100),
J < --min-norm (default 0.3; background poles above --background-above MeV are exempt).

Exit status: 0 without a flag, 1 with one, 2 on a usage or set-up error. ``--json PATH`` writes the rows.
"""
import argparse
import json
import sys

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


def widths_report(session, x):
    """(channel rows, level rows) at the fit vector x."""
    from pyazr import physics
    table = session.dimensionless_widths(x)
    channels = []
    levels = {}
    for c in table:
        jpi = c.jpi
        theta2, wigner_width = getattr(c, "theta2", None), getattr(c, "wigner_width", None)
        if not c.is_photon and c.is_open and wigner_width is None and c.value is not None:
            # a fixed channel: the session gives no Gamma_W; the engine's single-particle width at the level energy
            q = session.pairs.by_number(c.pair)
            wigner_width = 1e6 * float(physics.single_particle_width(c.L, q.Z1, q.Z2, q.M1, q.M2, q.channel_radius,
                                                                      c.level_energy - c.threshold))
            theta2 = abs(c.value) / wigner_width
        row = dict(name=c.name, jpi=jpi, level_energy=c.level_energy, pair=c.pair, L=c.L, S=c.S,
                   photon=bool(c.is_photon), open=bool(c.is_open), value=c.value, free=not c.fixed,
                   theta2=theta2, theta2_formal=getattr(c, "theta2_formal", None),
                   wigner_width=wigner_width, wigner_gamma2=getattr(c, "wigner_gamma2", None),
                   weisskopf=getattr(c, "weisskopf", None), wu=getattr(c, "wu", None))
        channels.append(row)
        key = (c.jgroup, c.level)
        lv = levels.setdefault(key, dict(jpi=jpi, level_energy=c.level_energy, norms=[]))
        if (not c.is_photon and c.is_open and row["theta2"] and row["theta2_formal"]
                and row["theta2_formal"] > 1e-12):
            lv["norms"].append(row["theta2"] / row["theta2_formal"])
    level_rows = []
    for key, lv in sorted(levels.items()):
        n = lv["norms"]
        level_rows.append(dict(jgroup=key[0], level=key[1], jpi=lv["jpi"], level_energy=lv["level_energy"],
                               norm=float(np.median(n)) if n else None,
                               norm_spread=float(np.ptp(n)) if len(n) > 1 else 0.0))
    if session.basis == 2:
        for lv, j in zip(level_rows, session.park_norms(x)):
            lv["park_norm"] = float(j)
    return channels, level_rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr")
    ap.add_argument("--params", dest="parameter_file",
                    help="parameter file (param.sav), by name; default: the project's values")
    ap.add_argument("--option", action="append", default=[], metavar="KEY=VALUE", help="a session option")
    ap.add_argument("--max-theta2", type=float, default=1.0)
    ap.add_argument("--max-wu", type=float, default=100.0)
    ap.add_argument("--min-norm", type=float, default=0.3)
    ap.add_argument("--background-above", type=float, default=float("inf"),
                    help="levels above this E_x (MeV) are background poles: no norm flag")
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    import pyazr
    try:
        s = pyazr.Session(a.azr, **{k: option_value(v) for k, v in (o.split("=", 1) for o in a.option)})
    except Exception as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    with s:
        x = s.parameters.read(a.parameter_file) if a.parameter_file else s.parameters.x0
        channels, levels = widths_report(s, x)
    flags = []
    fmt = lambda v, f: "-" if v is None else format(v, f)
    print(f"{'J^pi':6s} {'E_x MeV':>9s} {'pair':>4s} {'l':>2s} {'s':>4s} {'value':>12s} {'theta2':>10s} "
          f"{'theta2_f':>10s} {'Gamma_W eV':>11s} {'W.u.':>9s}")
    for c in channels:
        print(f"{c['jpi']:6s} {fmt(c['level_energy'], '9.4f')} {c['pair']:4d} {c['L']:2d} {c['S']:4.1f} "
              f"{fmt(c['value'], '12.5g')} {fmt(c['theta2'], '10.4g')} {fmt(c['theta2_formal'], '10.4g')} "
              f"{fmt(c['wigner_width'], '11.4g')} {fmt(c['wu'], '9.3g')}" + ("" if c["open"] else "  closed"))
        if not c["photon"] and c["open"] and c["theta2"] is not None and c["theta2"] > a.max_theta2:
            flags.append(f"{c['name']}: theta^2 = {c['theta2']:.3g} > {a.max_theta2:g}")
        if c["photon"] and c["wu"] is not None and c["wu"] > a.max_wu:
            flags.append(f"{c['name']}: {c['wu']:.3g} W.u. > {a.max_wu:g}")
    print(f"\n{'J^pi':6s} {'E_x MeV':>9s} {'J (norm)':>9s}")
    for lv in levels:
        j = lv.get("park_norm", lv["norm"])
        print(f"{lv['jpi']:6s} {fmt(lv['level_energy'], '9.4f')} {fmt(j, '9.4f')}")
        if j is not None and j < a.min_norm and (lv["level_energy"] or 0.0) <= a.background_above:
            flags.append(f"level {lv['jpi']} {lv['level_energy']:.4f} MeV: J = {j:.3g} < {a.min_norm:g} (saturated)")
    for f in flags:
        print("FLAG", f)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(dict(channels=channels, levels=levels, flags=flags), fh, indent=1)
    return 1 if flags else 0


if __name__ == "__main__":
    sys.exit(main())
