#!/usr/bin/env python3
"""The words and units of pyazr's answers, as the example scripts met them.

  (a) Session.level_scheme names each pair by its key in <levels> and its kind's word (particle, photon, beta), as
      Project.pairs() does: T3's beta pair (key 4, the engine's pair 2) is "4: beta", not "2: photon"; a channel's
      label carries the key and the kind; T1's photon pairs are "photon" with their keys; a channel's value is named
      by what the row holds: T1's bound 1/2- level holds an ANC (fm^-1/2), T3's beta channel a feeding amplitude.
  (b) Session.ec_integrals: every pathway carries its unit, fm^(L+1) with L its multipolarity.
  (c) evaluate(x).data[i].unit in a session of the extrapolations is the engine's unit of the extrapolation's quantity
      (T1: b for sigma, the unit of Project.describe's engine_unit), as in a session of the data.
  (d) a gamma_decay parameter's unit is its own part's (total_width eV, branching and mixing 1), not the kind's text.

Run from anywhere:  python3 tests/pyazr/library_answers_test.py
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
    import pyazr  # noqa: E402
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

warnings.simplefilter("ignore", pyazr.EngineWarning)
tmp = tempfile.mkdtemp(prefix="library_answers_")
try:
    for name in ("T1", "T3"):
        shutil.copytree(os.path.join(ROOT, "tests", "upgrade", name), os.path.join(tmp, name))
    shutil.copytree(os.path.join(ROOT, "tests", "gamma_decay", "models"), os.path.join(tmp, "gamma_decay"))
    T1 = os.path.join(tmp, "T1", "T1.azr")
    T3 = os.path.join(tmp, "T3", "T3.azr")

    # (a) the level scheme's pairs: the project's keys and kinds
    print("(a) Session.level_scheme")
    for path in (T1, T3):
        project_pairs = {p.key: p.kind_name for p in pyazr.Project.open(path).pairs()}
        with pyazr.Session(path) as s:
            scheme = s.level_scheme
            text = scheme.format()
        listed = {}
        block = text.split("Pairs (by their key in <levels>):", 1)[1].split("\n\n", 1)[0]
        for line in block.strip().splitlines():
            key, rest = line.strip().split(":", 1)
            listed[int(key)] = rest.split()[0]
        check(f"{os.path.basename(path)}: the pairs listed are the project's keys and kinds", listed == project_pairs,
              (listed, project_pairs))
        for lv in scheme.levels:
            for c in lv.channels:
                ok = c.pair_kind == project_pairs[c.pair_key] and (
                    f"pair {c.pair_key} " in c.label if c.pair_kind == "particle"
                    else c.label.endswith(f"({c.pair_kind} pair {c.pair_key})"))
                if not ok:
                    check(f"{os.path.basename(path)}: channel label {c.label!r}", False, (c.pair_key, c.pair_kind))
        if path == T1:
            bound = scheme.levels[0].channels[0]
            check("T1's bound 1/2- level: its pair-1 channel holds the ANC", bound.quantity == "anc"
                  and "C     =         1.85 fm^-1/2" in text, (bound, text[:800]))
    check("T3's beta pair is key 4 and called beta", "  4: beta" in text and "photon" not in text, text[:600])
    check("T3's beta channel: the feeding amplitude B, no unit", "G1 (beta pair 4)        B     =            1 "
          in text, text[:1200])

    # (b) the external-capture integrals' unit
    print("(b) Session.ec_integrals")
    with pyazr.Session(T1) as s:
        paths = s.ec_integrals(1, [0.5, 1.0])
    check("T1: pathways returned", len(paths) > 0, len(paths))
    check("every pathway carries unit fm^(L+1)",
          all(p.get("unit") == f"fm^{p['multipolarity'] + 1}" for p in paths), [p.get("unit") for p in paths])

    # (c) the unit of an extrapolation's model values
    print("(c) evaluate(x).data[i].unit, extrapolations")
    described = {e["name"]: e["engine_unit"] for e in pyazr.Project.open(T1).describe(as_json=True)["extrapolations"]}
    with pyazr.Session(T1, extrapolations=True) as s:
        data = s.evaluate().data
    check("T1: every extrapolation has the engine unit of Project.describe",
          data and all(d.unit == described[d.name] and d.unit for d in data), [(d.name, d.unit) for d in data])
    check("T1: a sigma extrapolation is in b", any(d.quantity == "sigma" and d.unit == "b" for d in data),
          [(d.quantity, d.unit) for d in data])

    # (d) a gamma_decay parameter's unit
    print("(d) the unit of a gamma_decay (radiative) parameter")
    with pyazr.Session(os.path.join(tmp, "gamma_decay", "record.azr")) as s:
        radiative = [p for p in s.parameters if p.kind == "radiative"]
    want = {"total_width": "eV", "branching": "1", "mixing": "1"}
    check("record.azr has radiative parameters of the three parts",
          {p.decay_part for p in radiative} == set(want), {p.decay_part for p in radiative})
    check("each has its part's unit", all(p.unit == want[p.decay_part] for p in radiative),
          [(p.name, p.unit) for p in radiative])
    check("the table's units are the records'", all(u == p.unit for u, p in zip(
        pyazr.Parameters(radiative).units, radiative)), "")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"{'FAILED: ' + str(len(failures)) if failures else 'all ok'}")
sys.exit(1 if failures else 0)
