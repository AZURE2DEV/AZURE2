#!/usr/bin/env python3
"""Pairs added and removed as Project edits, and the channels a J^pi allows (Project.add_pair, remove_pair,
allowed_channels, add_level(channels="allowed")).

  (a) T1's nuclear part built from nothing: an empty <levels>, its three pairs by add_pair (each held by a pair row
      until a level names it), its seven levels by add_level (explicit channels, allowed_channels with a limit, and
      channels="allowed" for a later level of a J^pi group) and set_channel; no pair row is left; that <levels> in
      place of T1's gives T1's chi-squared and parameter count from the engine;
  (b) add_pair then remove_pair gives the starting text back, byte for byte (T1, and the empty section);
  (c) allowed_channels by the Qt rule on T1's pairs: particle channels by parity and l_max, photon multipolarities
      by the triangle rule and the two limits, a beta pair by its parent's spin;
  (d) refusals, the project unchanged: remove_pair of a pair a level names (the levels and channels named), of the
      only channel of a level even with remove_channels, of a pair a data set names (the validator's), an unknown
      key; add_pair with a key in use, a radius that is not > 0, a photon pair given a light partner;
      remove_pair(remove_channels=True) of a photon pair no data set names removes its channel rows only.

Run from anywhere:  python3 tests/pyazr/pair_edits_test.py
"""
import os
import re
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
    from pyazr import Project, FormatError  # noqa: E402
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

T1 = os.path.join(ROOT, "tests", "upgrade", "T1", "T1.azr")
LEVELS = re.compile(rb"<levels>\r?\n.*?</levels>", re.S)


def evaluation(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pyazr.Session(path) as s:
            return s.evaluate(s.parameters.x0).chi2, len(s.parameters.x0)


tmp = tempfile.mkdtemp(prefix="pair_edits_")
try:
    print("(a) T1's nuclear part from nothing")
    t1 = Project.open(T1)
    p = Project.new(os.path.join(tmp, "built.azr"), nuclear_part=T1)
    for lv in reversed(p.levels()):
        p.remove_level(lv.label)
    empty = p.to_bytes()
    check("an empty <levels>", p.levels() == [] and p.pairs() == [], str(p.pairs()))
    proton = dict(J=0.5, parity=1, Z=1, mass=1.00782503223, g=5.5856946893)
    c12 = dict(J=0, parity=1, Z=6, mass=12.0)
    k1 = p.add_pair(light=proton, heavy=c12, separation_energy=1.9434900024034978, radius=3.6)
    n13 = dict(J=0.5, parity=-1, Z=7, mass=13.00573861)
    k2 = p.add_pair(kind="photon", heavy=n13, separation_energy=0.0, external_capture=7)
    k3 = p.add_pair(kind="photon", heavy=dict(n13, J=1.5), heavy_excitation=1.0, separation_energy=0.0,
                    external_capture=7)
    check("keys 1, 2, 3", (k1, k2, k3) == (1, 2, 3), str((k1, k2, k3)))
    check("the pairs read as T1's (columns of the first row)",
          [(q.key, q.kind, q.light, q.heavy, q.radius, q.separation_energy, q.heavy_excitation) for q in p.pairs()]
          == [(q.key, q.kind, q.light, q.heavy, q.radius, q.separation_energy, q.heavy_excitation)
              for q in t1.pairs()], str(p.pairs()))
    check("pair rows are no levels", p.levels() == [])

    def channels_of(level):
        return [dict(pair=c.pair, l=c.l, s=c.s, width=c.width, fixed=c.fixed) for c in level.channels]

    allowed_seen = {}
    for lv in t1.levels():
        here = {(x.J, x.parity) for x in p.levels()}
        if (lv.J, lv.parity) in here:
            p.add_level(lv.J, lv.parity, lv.energy, "allowed", fixed=lv.fixed)  # the group's channels, width 0
            for c in lv.channels:
                p.set_channel(f"{lv.jpi} {lv.energy_text}", c.pair, c.l, c.s, width=c.width, fixed=c.fixed)
            continue
        allowed = p.allowed_channels(lv.J, lv.parity, multipolarity_max=1)
        wanted = [dict(pair=c.pair, l=c.l, s=c.s) for c in lv.channels]
        allowed_seen[lv.label] = allowed == wanted
        p.add_level(lv.J, lv.parity, lv.energy, channels_of(lv), fixed=lv.fixed)
    check("allowed_channels(multipolarity_max=1) gives T1's channel sets of 1/2+ and 5/2+",
          allowed_seen.get("1/2+ 2.365") and allowed_seen.get("5/2+ 3.547"), str(allowed_seen))
    check("no pair row is left once levels name the pairs", b" 0 0 1 0 " not in LEVELS.search(p.to_bytes()).group(0)
          and len(p.levels()) == len(t1.levels()))
    check("the levels read as T1's (J^pi, E_x, fixed, channels with width and fixed flag)",
          [(x.jpi, x.energy, x.fixed, [(c.pair, c.l, c.s, c.width, c.fixed) for c in x.channels]) for x in p.levels()]
          == [(x.jpi, x.energy, x.fixed, [(c.pair, c.l, c.s, c.width, c.fixed) for c in x.channels])
              for x in t1.levels()])
    d = os.path.join(tmp, "T1")
    shutil.copytree(os.path.dirname(T1), d, ignore=shutil.ignore_patterns("output*", "checks*"))
    original = open(T1, "rb").read()
    built = LEVELS.sub(lambda m: LEVELS.search(p.to_bytes()).group(0), original)
    check("the built <levels> differs in text from T1's (it is not a copy)", built != original)
    with open(os.path.join(d, "built.azr"), "wb") as f:
        f.write(built)
    reference, built_eval = evaluation(os.path.join(d, "T1.azr")), evaluation(os.path.join(d, "built.azr"))
    check(f"the engine on T1 with the built <levels>: chi2 {built_eval[0]!r} and {built_eval[1]} parameters as T1's",
          built_eval == reference, f"{built_eval} vs {reference}")

    print("(b) add then remove")
    for name, start in (("T1", None), ("the empty section", empty)):
        q = Project.open(T1) if start is None else Project(start, os.path.join(tmp, "e.azr"))
        before = q.to_bytes()
        k = q.add_pair(light=dict(J=0, parity=1, Z=2, mass=4.00260325413), heavy=dict(J=1.5, parity=-1, Z=4,
                       mass=9.0121831), separation_energy=10.6, radius=5.0, binding_energy=2.2)
        check(f"{name}: the new pair {k} is listed", any(x.key == k and x.binding_energy == 2.2 for x in q.pairs()))
        q.remove_pair(k)
        check(f"{name}: add_pair then remove_pair gives the text back", q.to_bytes() == before)

    print("(c) allowed_channels")
    q = Project.open(T1)
    got = q.allowed_channels(1.5, -1)
    check("3/2-: p + 12C l = 1; M1 + E2 to 1/2-; M1 + E2 to 3/2- (multipolarity_max 2)",
          [(c["pair"], c["l"], c["s"]) for c in got] == [(1, 1, 0.5), (2, 1, 0.5), (2, 2, 0.5), (3, 1, 1.5),
                                                         (3, 2, 1.5)], str(got))
    got = q.allowed_channels(2.5, 1, l_max=1, multipolarities_max=1)
    check("5/2+ with l_max 1: no proton channel (l = 2 needed); E2 to 1/2-; one multipolarity (M1) to 3/2-",
          [(c["pair"], c["l"], c["s"]) for c in got] == [(2, 2, 0.5), (3, 1, 1.5)], str(got))
    got = q.allowed_channels(0.5, 1, pairs=[1])
    check("1/2+ of pair 1 only: l = 0", [(c["pair"], c["l"], c["s"]) for c in got] == [(1, 0, 0.5)], str(got))
    b = Project(empty, os.path.join(tmp, "b.azr"))
    kb = b.add_pair(kind="beta", light=dict(J=0.5, parity=1, Z=-1, mass=0.000548579909), heavy=dict(J=1, parity=1,
                    Z=6, mass=12.0), separation_energy=0.0)
    check("beta, parent 1+: 1+ gives Gamow-Teller and Fermi, 0+ and 2+ Gamow-Teller, 1- none",
          [[(c["l"], c["s"]) for c in b.allowed_channels(J, par)] for J, par in ((1, 1), (0, 1), (2, 1), (1, -1))]
          == [[(1, 0.0), (0, 0.0)], [(1, 1.0)], [(1, 1.0)], []] and kb == 1)

    print("(d) refusals and remove_channels")

    def refused(what, call, *words, exc=FormatError):
        q = Project.open(T1)
        before = q.to_bytes()
        try:
            call(q)
            check(f"{what} refused", False, "no error")
        except exc as e:
            check(f"{what} refused, the project unchanged",
                  q.to_bytes() == before and all(w in str(e) for w in words), str(e)[:400])

    refused("remove_pair of a pair levels name", lambda q: q.remove_pair(2), "'3/2- 1.0' (l = 1, s = 1/2; l = 2",
            "remove_channels=True")
    refused("remove_pair of the only channel of a level", lambda q: q.remove_pair(1, remove_channels=True),
            "only channel", "'1/2- 0.0'", "remove_level")
    refused("remove_pair of a pair a data set names", lambda q: q.remove_pair(3, remove_channels=True),
            "E_UNKNOWN_PAIR", "cap_g1.reaction")
    refused("remove_pair of an unknown key", lambda q: q.remove_pair(9), "no pair has the key 9", "the pairs: 1")
    refused("add_pair with a key in use", lambda q: q.add_pair(light=proton, heavy=c12, separation_energy=1.94,
                                                               radius=3.6, key=1), "key 1 is taken", "next free")
    refused("add_pair with a radius that is not > 0", lambda q: q.add_pair(light=proton, heavy=c12,
                                                                           separation_energy=1.94, radius=-1),
            "-1", "> 0")
    refused("a photon pair given a light partner", lambda q: q.add_pair(kind="photon", light=proton, heavy=n13,
                                                                        separation_energy=0.0),
            "photon is implied", exc=ValueError)
    q = Project.new(os.path.join(tmp, "nodata.azr"), nuclear_part=T1)  # no data set names pair 3
    q.remove_pair(3, remove_channels=True)
    check("remove_pair(3, remove_channels=True) without data: pair 3's channel rows go, every level stays",
          [x.label for x in q.levels()] == [x.label for x in t1.levels()]
          and all(c.pair != 3 for x in q.levels() for c in x.channels)
          and [x.key for x in q.pairs()] == [1, 2], str([x.key for x in q.pairs()]))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok")
