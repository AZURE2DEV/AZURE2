#!/usr/bin/env python3
"""The words and units of pyazr's answers, as the example scripts met them.

  (a) Session.level_scheme names each pair by its key in <levels> and its kind's word (particle, photon, beta), as
      Project.pairs() does: T3's beta pair (key 4, the engine's pair 2) is "4: beta", not "2: photon"; a channel's
      label carries the key and the kind; T1's photon pairs are "photon" with their keys; a channel's value is named
      by what the row holds: T1's bound 1/2- level holds an ANC (fm^-1/2), T3's beta channel a feeding amplitude.

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
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"{'FAILED: ' + str(len(failures)) if failures else 'all ok'}")
sys.exit(1 if failures else 0)
