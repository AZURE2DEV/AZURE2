#!/usr/bin/env python3
"""Edit a project from Python: levels and data sets added, removed and switched off; the edited copy calculated.

``pyazr.Project`` edits what a project is made of -- the level scheme (``<levels>``: a level by its label, a channel
by pair, l and s) and the records (``dataset``, ``extrapolation``, ...) -- with the engine's own grammars: every edit
is checked as it is made and a refused one raises ``pyazr.FormatError`` with the project unchanged.  The edits here:

  1. a background pole removed (``remove_level``), a 3/2+ resonance added with the channels its J^pi allows
     (``allowed_channels``, ``add_level``),
  2. a resonance switched off in the file (``set_level(label, active=False)``: its rows stay, the engine reads none
     of them; ``active=True`` brings it back),
  3. a data set removed (``remove``) and one added on an existing data file with another energy window
     (``add_dataset``),
  4. a prior on a level's energy and limits on one of its widths (``set_parameter_setting``: rows of
     ``<parameterSettings>``, read back by ``parameter_settings()``); the edited session's objective holds the
     prior as its term "nuisance priors", and the parameter's bounds hold the limits.

The edited project is written to ``--out`` (default: a new temporary directory, printed); the original is not
changed.  Its session reads the data files from the original project's directory (``Session(..., cwd=)``).  The
external-capture integrals need no care: the engine checks the signature of a saved ``intEC.dat`` and recomputes
what does not match.  The result is not a physical model; the point is the mechanics of the edits.

Usage (from the repository's root):

    python3 pyazr/examples/edit_model.py tests/upgrade/T1/T1.azr [--out DIR]
"""
import argparse
import os
import sys
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "2")

from pyazr import Project, Session  # noqa: E402


def summary(path, cwd):
    """Data sets, points, free parameters and the data chi-squared of a project's session."""
    with Session(path, cwd=cwd) as s:
        ev = s.evaluate()
        return len(ev.data), sum(d.n for d in ev.data), len(s.parameters.x0), ev.chi2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr); tested on tests/upgrade/T1/T1.azr")
    ap.add_argument("--out", default=None, help="the directory of the edited copy (default: a new temporary one)")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(args.project))
    out_dir = args.out or tempfile.mkdtemp(prefix="edit_model.")
    os.makedirs(out_dir, exist_ok=True)

    p = Project.open(args.project)
    print(f"{os.path.basename(args.project)}: levels {[lv.label for lv in p.levels()]}")

    # -- 1. levels: remove a background pole, add a resonance -------------------------------------------------------
    pole = max(p.levels(), key=lambda lv: lv.energy)                    # the highest level: a background pole
    p.remove_level(pole.label)
    print(f"\nremoved the level {pole.label} (E_x in MeV)")
    channels = p.allowed_channels(1.5, +1)                               # what a 3/2+ level can couple to
    print(f"a 3/2+ level allows {[(c['pair'], c['l'], c['s']) for c in channels]} (pair, l, s)")
    particle = [c for c in channels if c["pair"] == 1]
    p.add_level(1.5, +1, 4.1, channels=[dict(c, width=1000.0) for c in particle])   # width in eV
    print(f"added the level {p.levels()[-1].label}: Gamma = 1000 eV in {[(c['l'], c['s']) for c in particle]}")

    # -- 2. a resonance switched off in the file -----------------------------------------------------------------------
    off = [lv for lv in p.levels() if not lv.fixed][-1]                   # a resonance whose energy is fitted
    p.set_level(off.label, active=False)
    print(f"switched off the level {off.label} (its rows stay in the file)")

    # -- 3. data sets ----------------------------------------------------------------------------------------------
    records = p.records()
    data = [name for name, r in records.items() if r.record_kind == "dataset"]
    p.remove(data[-1])
    print(f"\nremoved the data set {data[-1]}")
    first = dict(records[data[0]].items())                               # its fields, as the file has them
    p.add_dataset(f"{data[0]}_low", **dict(first, energy_window="0.3 1.0 MeV"))
    print(f"added the data set {data[0]}_low: {first['file']} between 0.3 and 1.0 MeV (lab)")

    # -- 4. a prior and limits on level parameters --------------------------------------------------------------------
    level = next(lv for lv in p.levels() if lv.active and not lv.fixed)  # a level whose energy is fitted
    channel = next(c for c in level.channels if not c.fixed)              # and one of its fitted widths
    p.set_parameter_setting(level.label, prior=(level.energy - 0.002, 0.01))             # MeV
    p.set_parameter_setting(level.label, (channel.pair, channel.l, channel.s),
                            limits=(0.5 * channel.width, 1.25 * channel.width))          # the row's unit (eV here)
    print()
    for row in p.parameter_settings():
        print(f"<parameterSettings> {row.label} -> level {row.level}, channel {row.channel}: prior {row.prior}, "
              f"limits {row.limits}")

    # -- 5. the edited copy, beside the original's data -------------------------------------------------------------
    path = p.save(os.path.join(out_dir, os.path.basename(args.project)))
    print(f"\nwrote {path}")
    n0, m0, f0, c0 = summary(args.project, here)
    n1, m1, f1, c1 = summary(path, here)
    print(f"original: {n0} data sets, {m0} points, {f0} free parameters, chi2 = {c0:.4f}")
    print(f"edited  : {n1} data sets, {m1} points, {f1} free parameters, chi2 = {c1:.4f}")
    with Session(path, cwd=here) as s:
        prior = s.evaluate().terms.get("nuisance priors", 0.0)
        limited = [(q.name, q.bounds) for q in s.parameters.free if all(abs(b) != float("inf") for b in q.bounds)]
    print(f"edited  : the prior's term of the objective, nuisance priors = {prior:.4f}; limited: "
          + ", ".join(f"{n} in [{lo:.4g}, {hi:.4g}] (fit units)" for n, (lo, hi) in limited))
    return 0


if __name__ == "__main__":
    sys.exit(main())
