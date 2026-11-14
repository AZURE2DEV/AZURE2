#!/usr/bin/env python3
"""A prior or a limit on a level parameter as a Project edit (Project.set_parameter_setting, <parameterSettings>).

On a copy of tests/upgrade/T1 (no <parameterSettings>):

  (a) a prior on the energy of 3/2- 3.502 (free): the row "Level 3 Energy (MeV)" is written in a new section; a
      session of the saved project has the objective term "nuisance priors" = ((3.502 - 3.5)/0.01)^2 and the data
      chi-squared of before;
  (b) limits on the width of its channel (pair 1, l 1, s 1/2): the row "Level 3 Channel 1 Width (eV)"; the session's
      parameter width_3_1 has finite bounds (the limits converted to its amplitude) that hold its start, and the
      other parameters' bounds are unchanged; the engine's numbers of the rows are those of its parameter names
      (energy_3, width_3_1);
  (c) Project.parameter_settings() reads the rows back as records, with the level and channel they name;
      kinds()["parameter_settings"] lists the fields with units;
  (d) the same parameter again replaces its row; both None removes it; the section stays valid;
  (e) refusals with their three parts, the project unchanged: a fixed parameter (the bound level 1/2- 0.0's energy),
      a channel the level has not, a sigma <= 0, limits not in order, an inactive level.

Run from anywhere:  python3 tests/pyazr/parameter_settings_test.py
"""
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import numpy as np  # noqa: E402
    import pyazr  # noqa: E402
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

warnings.simplefilter("ignore", pyazr.EngineWarning)
LEVEL = "3/2- 3.502"
tmp = tempfile.mkdtemp(prefix="parameter_settings_")
try:
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), os.path.join(tmp, "T1"),
                    ignore=shutil.ignore_patterns("output", "checks"))
    T1 = os.path.join(tmp, "T1", "T1.azr")
    with pyazr.Session(T1) as s:
        before = s.evaluate()
        bounds_before = {p.name: tuple(p.bounds) for p in s.parameters}

    print("(a), (b) a prior on an energy, limits on a width")
    p = pyazr.Project.open(T1)
    check("T1 has no <parameterSettings>", "<parameterSettings>" not in p.text)
    r1 = p.set_parameter_setting(LEVEL, prior=(3.5, 0.01))
    r2 = p.set_parameter_setting(LEVEL, (1, 1, 0.5), limits=(1e3, 1e5))
    check("the rows' labels", r1.label == "Level 3 Energy (MeV)" and r2.label == "Level 3 Channel 1 Width (eV)",
          (r1, r2))
    path = p.save(os.path.join(tmp, "T1", "settings.azr"))
    with pyazr.Session(path) as s:
        ev = s.evaluate()
        names = [q.name for q in s.parameters]
        energy, width = s.parameters.by_name("energy_3"), s.parameters.by_name("width_3_1")
        check("the rows name the engine's parameters energy_3 and width_3_1",
              energy.level_label == LEVEL and width.level_label == LEVEL and width.channel == 1
              and energy.level_number == 3 and width.level_number == 3, (energy, width))
        check("the objective has the prior: nuisance priors = ((3.502 - 3.5)/0.01)^2",
              abs(ev.terms.get("nuisance priors", 0.0) - 0.04) < 1e-9, ev.terms)
        check("the data chi-squared is that of before", ev.chi2 == before.chi2, (ev.chi2, before.chi2))
        lo, hi = width.bounds
        check("width_3_1 is limited, its start inside", np.isfinite(lo) and np.isfinite(hi) and lo < width.start < hi,
              (width.bounds, width.start))
        others = {q.name: tuple(q.bounds) for q in s.parameters if q.name != "width_3_1"}
        check("no other parameter's bounds changed", all(others[n] == bounds_before[n] for n in others),
              [n for n in others if others[n] != bounds_before[n]])
        f = s.objective()
        check("the Objective's bounds are the parameters'", dict(zip(f.names, f.bounds))["width_3_1"] == (lo, hi), "")

    print("(c) the records")
    rows = pyazr.Project.open(path).parameter_settings()
    check("two rows read back", [(r.label, r.level, r.channel, r.prior, r.limits) for r in rows] ==
          [("Level 3 Energy (MeV)", LEVEL, None, (3.5, 0.01), None),
           ("Level 3 Channel 1 Width (eV)", LEVEL, (1, 1, 0.5), None, (1e3, 1e5))], rows)
    k = pyazr.kinds()["parameter_settings"]
    check("kinds()['parameter_settings']: fields with units, the engine's columns",
          {"value", "error", "lower", "upper"} <= set(k["fields"]) and all("unit" in v for v in k["fields"].values())
          and k["columns"][:3] == ["value", "lower", "upper"], k)

    print("(d) replace and remove")
    p = pyazr.Project.open(path)
    p.set_parameter_setting(LEVEL, prior=(3.51, 0.02), limits=(3.4, 3.6))
    rows = p.parameter_settings()
    check("the energy's row replaced, not added", len(rows) == 2 and rows[0].prior == (3.51, 0.02)
          and rows[0].limits == (3.4, 3.6), rows)
    check("both None removes the row", p.set_parameter_setting(LEVEL, (1, 1, 0.5)) is None
          and [r.label for r in p.parameter_settings()] == ["Level 3 Energy (MeV)"], p.parameter_settings())
    check("the project validates", not [d for d in p.validate() if d["severity"] == "error"], "")

    print("(e) refusals")
    p = pyazr.Project.open(T1)
    text = p.text
    for label, call, words in (
            ("a fixed parameter", lambda: p.set_parameter_setting("1/2- 0.0", prior=(0.0, 0.1)),
             ["is fixed", "set_level('1/2- 0.0', fixed=False)"]),
            ("a channel the level has not", lambda: p.set_parameter_setting(LEVEL, (1, 3, 0.5), limits=(0, 1)),
             ["no such channel", "(1, 1, 0.5)"]),
            ("sigma <= 0", lambda: p.set_parameter_setting(LEVEL, prior=(3.5, 0.0)), ["sigma > 0"]),
            ("limits not in order", lambda: p.set_parameter_setting(LEVEL, limits=(4.0, 3.0)), ["lower < upper"]),
            ("an unknown level", lambda: p.set_parameter_setting("7/2+ 9.9", prior=(1, 1)), ["7/2+"])):
        try:
            call()
            check(f"{label}: refused", False, "accepted")
        except pyazr.FormatError as e:
            parts = (e.what, e.allowed, e.fix)
            check(f"{label}: refused with its parts in the message, saying {words}",
                  (all(parts) and all(x in str(e) for x in parts) or label == "an unknown level")
                  and all(w in str(e) for w in words), (parts, str(e)))
    check("the project is unchanged", p.text == text, "")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p.set_level("5/2+ 3.547", active=False)
    try:
        p.set_parameter_setting("5/2+ 3.547", prior=(3.5, 0.1))
        check("an inactive level: refused", False, "accepted")
    except pyazr.FormatError as e:
        check("an inactive level: refused, saying how to activate it", "inactive" in e.what and "active=True" in e.fix,
              str(e))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"{'FAILED: ' + str(len(failures)) if failures else 'all ok'}")
sys.exit(1 if failures else 0)
