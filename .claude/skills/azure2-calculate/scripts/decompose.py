#!/usr/bin/env python3
"""Decompose an AZURE2 cross section into what one level does, its bare resonance and its interference.

At fixed parameters (no refit), with the free widths of chosen levels set to zero in the fit vector
(``Parameters.of_level(label)``, ``Parameters.put``), four evaluations of one data set or extrapolation:

    full   every level as fitted
    off    the level's free widths zeroed: the model without it
    only   every other level's free widths zeroed (bound final states kept, see below)
    bare   every level's free widths zeroed (bound final states kept): hard sphere and direct capture

and writes per point: E_cm, E_x, full, off, only, bare, full - off (what the level does: resonance plus its
interference), only - bare (the bare resonance) and their difference (the interference). With ``--with LEVEL2``
also the pairwise term only{T,L} - only{T} - only{L} + bare = 2 Re(A_T A_L*), which is zero for an
angle-integrated cross section when the two levels have different J^pi (levels of different J^pi do not interfere
there) -- a check of any decomposition.

Bound final states: an external-capture amplitude scales with the ANC of its final state, which is the width of a
closed channel of a bound level. ``only`` and ``bare`` keep the levels below every particle threshold of their
pairs (the final states) unless ``--zero-bound``. Fixed widths are not entries of the fit vector and stay as they
are; the script lists them (zero them in a copy of the project with ``Project.set_channel`` if they must go).

Levels are named by ``Parameter.level_label`` (``"1/2+ 2.365"``, as ``Project.levels()`` gives them).
Run on the data sets (default) or on the project's extrapolations (``--extrapolations``).
Exit status: 0, 1 when a level or data set is unknown or the engine refuses, 2 on a usage error.
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


def free_widths(P, labels):
    """The free width rows of the levels ``labels``."""
    rows = []
    for lab in labels:
        rows += [p for p in P.of_level(lab).of_kind("width") if not p.fixed]
    return rows


def bound_levels(session):
    """Labels of levels below every particle threshold of the pairs of their channels (capture final states)."""
    thresholds = {q.number: q.threshold for q in session.pairs if not q.is_photon}
    out = set()
    for lab in dict.fromkeys(p.level_label for p in session.parameters if p.level_label):
        rows = session.parameters.of_level(lab)
        e = next((p.level_energy for p in rows if p.level_energy is not None), 0.0)
        pairs = [p.pair for p in rows.of_kind("width") if p.pair in thresholds]
        if pairs and all(e < thresholds[q] for q in pairs):
            out.add(lab)
    return out


def decompose(session, x, level, dataset, other=None, zero_bound=False):
    P = session.parameters
    labels = list(dict.fromkeys(p.level_label for p in P if p.level_label))
    for lab in [level] + ([other] if other else []):
        if lab not in labels:
            raise KeyError(f"no level {lab!r}; the levels: {labels}")
    keep = set() if zero_bound else bound_levels(session)

    def model(zeroed):
        rows = free_widths(P, zeroed)
        xx = P.put(x, rows, 0.0) if rows else np.array(x, float)
        ev = session.evaluate(xx)
        d = next((d for d in ev.data if d.name == dataset), None) if isinstance(dataset, str) else ev.data[dataset]
        if d is None:
            raise KeyError(f"no data set or extrapolation {dataset!r}; the names: {[d.name for d in ev.data]}")
        return d, np.asarray(d.model, float)

    others = [lab for lab in labels if lab not in keep]
    d, full = model([])
    out = dict(energies=np.asarray(d.energies, float), excitation_energies=np.asarray(d.excitation_energies, float),
               full=full, off=model([level])[1], only=model([lab for lab in others if lab != level])[1],
               bare=model(others)[1], kept=sorted(keep),
               fixed=[p.name for lab in labels for p in P.of_level(lab).of_kind("width") if p.fixed])
    if other:
        out["only_both"] = model([lab for lab in others if lab not in (level, other)])[1]
        out["only_other"] = model([lab for lab in others if lab != other])[1]
        out["pair_term"] = out["only_both"] - out["only"] - out["only_other"] + out["bare"]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr")
    ap.add_argument("--level", required=True, help='the level, "1/2+ 2.365"')
    ap.add_argument("--with", dest="other", help="a second level for the pairwise interference term")
    ap.add_argument("--dataset", required=True, help="the data set's (or extrapolation's) record name")
    ap.add_argument("--extrapolations", action="store_true", help="evaluate the project's extrapolations")
    ap.add_argument("--params", dest="parameter_file", help="parameter file (param.sav), by name")
    ap.add_argument("--option", action="append", default=[], metavar="KEY=VALUE", help="a session option")
    ap.add_argument("--zero-bound", action="store_true", help="also zero the bound final states' ANCs")
    ap.add_argument("--out", help="write the table here (default: standard output)")
    a = ap.parse_args(argv)
    import pyazr
    opts = {k: option_value(v) for k, v in (o.split("=", 1) for o in a.option)}
    try:
        with pyazr.Session(a.azr, extrapolations=a.extrapolations, **opts) as s:
            x = s.parameters.read(a.parameter_file) if a.parameter_file else s.parameters.x0
            r = decompose(s, x, a.level, a.dataset, a.other, a.zero_bound)
    except (KeyError, pyazr.Error) as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    cols = ["full", "off", "only", "bare"] + (["pair_term"] if a.other else [])
    lines = [f"# level {a.level}" + (f" with {a.other}" if a.other else "") + f", {a.dataset}; kept (final states): "
             f"{r['kept'] or 'none'}; fixed widths left as they are: {r['fixed'] or 'none'}",
             "# " + " ".join(f"{c:>14s}" for c in ["E_cm_MeV", "E_x_MeV"] + cols + ["level", "resonance",
                                                                                   "interference"])]
    for i in range(r["full"].size):
        level, res = r["full"][i] - r["off"][i], r["only"][i] - r["bare"][i]
        vals = [r["energies"][i], r["excitation_energies"][i]] + [r[c][i] for c in cols] + [level, res, level - res]
        lines.append("  " + " ".join(f"{v:14.7e}" for v in vals))
    text = "\n".join(lines) + "\n"
    if a.out:
        Path(a.out).write_text(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
