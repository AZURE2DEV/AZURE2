#!/usr/bin/env python3
"""pyazr.Project on the nuclear part of a format-3 project: <levels>, <config> and <thm> read and edited by the engine
(docs/upgrade/pyazr_on_records.md, "Edits of the nuclear part").

On tests/13N (format 3):
  (a) levels() names every level as the engine numbers it: the energy of each level parameter of a session, by its
      (jgroup, level), is the energy of the Level with that (group, level); every row once; pairs() as <levels> says;
  (b) an edit changes only the rows it names, and only their tokens: set_level energy / fixed, set_channel width
      (right-aligned in its column), set_pair radius on every row of the pair; a no-op edit keeps every byte;
  (c) add_level: the rows of the J^pi group (channels not given: width 0, fixed), or a new J^pi from a row of each
      pair; the next free level ID; remove_level of it gives back the file byte for byte;
  (d) refusals leave the project unchanged, each with the engine's reason: an unknown or ambiguous level, a channel
      outside the group's set, a width that is not finite, an unknown pair; an unknown <config> entry, a value the
      engine would not read, the format;
  (e) <config>: config() as the engine reads it; set_config of a positional entry keeps its comment column, a key is
      added, replaced and removed;
  (f) save_fit writes the fitted level values into the rows they come from: the saved file run again gives the fit's
      chi2 (rel 1e-6, save_fit's own check in the Brune basis: the physical values read back through the transform),
      a fixed level stays fixed, and a level whose widths are all zero keeps them.
On examples/o18_lacognata2010 (format 3, a THM experiment):
  (g) thm() reads the options and the experiment; set_thm / set_thm_experiment write them, the block checked as the
      engine reads a format-3 one; a misspelt option, a refused experiment key, segments= and the removal of an
      experiment a data set names are refused.

Run from anywhere:  python3 tests/pyazr/project_nuclear_part_test.py
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
    from pyazr import Project, FormatError  # noqa: E402
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

SRC = os.path.join(ROOT, "tests", "13N", "13N.azr")
THM = os.path.join(ROOT, "examples", "o18_lacognata2010", "o18_lacognata2010.azr")


def changed_lines(a, b):
    a, b = a.split(b"\n"), b.split(b"\n")
    if len(a) != len(b):
        return None
    return [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]


def refused(what, call, *codes, prepare=None):
    p = Project.open(SRC)
    if prepare:
        prepare(p)
    before = p.to_bytes()
    try:
        call(p)
        check(f"(d) {what} refused", False, "no error")
    except FormatError as e:
        check(f"(d) {what} refused, the project unchanged", p.to_bytes() == before and all(c in str(e) for c in codes),
              str(e)[:200])


print("(a) levels and pairs")
p = Project.open(SRC)
levels = p.levels()
check("(a) twelve levels, every row once", len(levels) == 12 and sum(len(lv.channels) for lv in levels) == 30,
      [len(levels), sum(len(lv.channels) for lv in levels)])
check("(a) the labels", [lv.label for lv in levels[:3]] == ["1/2- 0", "1/2- 20", "1/2+ 2.3682"],
      [lv.label for lv in levels[:3]])
check("(a) pairs", [(q.key, q.label, q.radius) for q in p.pairs()] == [(1, "p+12C", 3.6), (2, "g+13N", 0.0)],
      p.pairs())
check("(a) a level by J^pi and E_x within 1 keV", p.level(("3/2-", 3.5035)).label == "3/2- 3.5032")
check("(a) the reaction from pair keys", p.reaction(1, 2) == "12C(p,g)13N" and p.reaction(1, -1) == "12C(p,g)13N[total]",
      (p.reaction(1, 2), p.reaction(1, -1)))
work = tempfile.mkdtemp(prefix="pyazr_nuclear_part_")
try:
    shutil.copytree(os.path.join(ROOT, "tests", "13N"), os.path.join(work, "13N"),
                    ignore=shutil.ignore_patterns("output", "checks"))
    os.makedirs(os.path.join(work, "13N", "output"))
    os.makedirs(os.path.join(work, "13N", "checks"))
    model = os.path.join(work, "13N", "13N.azr")
    with pyazr.Session(model, cwd=os.path.dirname(model)) as a:
        by_index = {(lv.group, lv.level): lv for lv in levels}
        energies = [(prm.jgroup, prm.level, prm.level_energy) for prm in a.parameters if prm.kind == "energy"]
        check("(a) every level parameter's (jgroup, level) names the level of that energy",
              all((g, n) in by_index and abs(by_index[(g, n)].energy - (e or 0.0)) < 1e-9 for g, n, e in energies),
              energies)

    print("(b) edits keep the bytes they do not change")
    before = p.to_bytes()
    p.set_level("3/2- 3.5032", energy=3.5032)
    check("(b) a no-op edit keeps every byte", p.to_bytes() == before)
    rows = [c.line for c in p.level("3/2- 3.5032").channels]
    p.set_level("3/2- 3.5032", energy=3.51, fixed=False)
    check("(b) set_level changes the level's rows only", changed_lines(before, p.to_bytes()) == rows,
          (changed_lines(before, p.to_bytes()), rows))
    lv = p.level("3/2- 3.51")
    check("(b) the energy and the fixed flag read back", lv.energy == 3.51 and not lv.fixed, lv)
    line = p.text.split("\n")[lv.channel(1, 1, 0.5).line - 1]
    p.set_channel(lv, 1, 1, 0.5, width=6e4)
    new = p.text.split("\n")[lv.channel(1, 1, 0.5).line - 1]
    check("(b) set_channel: the width right-aligned in its column, the line's length kept",
          len(new) == len(line) and new.replace("60000", "55200") == line, (line, new))
    after = p.to_bytes()
    p.set_pair(1, radius=4.0)
    pair1 = sorted(c.line for lv in p.levels() for c in lv.channels if c.pair == 1)
    check("(b) set_pair: every row of the pair, no other", changed_lines(after, p.to_bytes()) == pair1)
    check("(b) the radius reads back", [q.radius for q in p.pairs()] == [4.0, 0.0])

    print("(c) adding and removing a level")
    q = Project.open(SRC)
    before = q.to_bytes()
    q.add_level(J=1.5, parity=-1, energy=7.0, channels=[dict(pair=1, l=1, s=0.5, width=100.0)])
    new = q.level("3/2- 7")
    check("(c) the group's three channels, those not given width 0 and fixed",
          [(c.pair, c.l, c.s, c.width, c.fixed) for c in new.channels]
          == [(1, 1, 0.5, 100.0, False), (2, 1, 0.5, 0.0, True), (2, 2, 0.5, 0.0, True)], new.channels)
    check("(c) the next free level ID, the engine's next level of the group", new.ids == (13,) and new.level == 3,
          (new.ids, new.level))
    q.remove_level("3/2- 7")
    check("(c) removing it gives back the file byte for byte", q.to_bytes() == before)
    q.add_level(J=4.5, parity=1, energy=9.0, channels=[dict(pair=1, l=4, s=0.5, width=10.0), dict(pair=2, l=2, s=0.5)],
                fixed=False)
    new = q.level("9/2+ 9")
    check("(c) a new J^pi: one row per channel, as given",
          [(c.pair, c.l, c.s, c.width, c.fixed) for c in new.channels] == [(1, 4, 0.5, 10.0, False), (2, 2, 0.5, 0.0, False)]
          and not new.fixed and new.group == 9 and new.level == 1, new)

    print("(d) refusals")
    refused("an unknown level", lambda x: x.set_level("3/2- 99", energy=1.0), "no 3/2- level within 1 keV")
    refused("an ambiguous J^pi and E_x", lambda x: x.set_level(("3/2-", 20.0), energy=21.0), "names 2 levels",
            prepare=lambda x: x.add_level(J=1.5, parity=-1, energy=20.0005, channels=[]))
    refused("a channel outside the group's set", lambda x: x.add_level(J=0.5, parity=1, energy=5.0,
            channels=[dict(pair=1, l=2, s=0.5)]), "share one channel set")
    refused("a width that is not finite", lambda x: x.set_channel("1/2- 0", 1, 1, 0.5, width=float("nan")), "finite")
    refused("an unknown pair", lambda x: x.set_pair(7, radius=5.0), "no pair has the key 7")
    refused("a row the level does not have", lambda x: x.set_channel("1/2- 0", 1, 2, 0.5, width=1.0), "has no row")
    refused("an unknown <config> entry", lambda x: x.set_config(relativstic=True), "not an entry of <config>")
    refused("the removed key useAdaptiveGrid", lambda x: x.set_config(useAdaptiveGrid=True), "not an entry of <config>")
    refused("a value the engine would not read", lambda x: x.set_config(check_data="print"), "none, screen or file")
    refused("the format", lambda x: x.set_config(format="2"), "not edited")

    print("(e) <config>")
    q = Project.open(SRC)
    c = q.config()
    check("(e) config() as the engine reads it", c["output_dir"] == "output/" and c["amatrix"] == "true"
          and "useAdaptiveGrid" not in c and c["format"] == "3", c)
    old_line = q.text.split("\n")[2]
    q.set_config(output_dir="run2/", relativistic=True)
    new_line = q.text.split("\n")[2]
    check("(e) the output directory, its comment in its column",
          new_line.index("#") == old_line.index("#") and new_line.startswith("run2/"), new_line)
    check("(e) a key added", q.config()["relativistic"] == "1")
    q.set_config(relativistic=None)
    check("(e) a key removed", "relativistic" not in q.config())
    # useAdaptiveGrid is no key of format 3: a file that still carries the line is refused by the engine's validator
    with open(SRC, "rb") as f:
        text = f.read().replace(b"</config>", b"useAdaptiveGrid=1\n</config>", 1)
    old_key = Project(text, SRC)
    codes = [d.get("code") for d in old_key.validate()]
    check("(e) useAdaptiveGrid in <config>: E_CONFIG_KEY_REMOVED", codes == ["E_CONFIG_KEY_REMOVED"], codes)

    print("(f) save_fit on the format-3 project")
    with pyazr.Session(model, cwd=os.path.dirname(model)) as a:
        x = np.asarray(a.parameters.x0, float) * 1.001
        chi = a.evaluate(x).chi2
        fixed_levels = {(prm.jgroup, prm.level) for prm in a.parameters if prm.kind == "energy" and prm.fixed}
        snap, _ = a.save_fit(os.path.join(work, "13N", "snap.azr"), x=x, param_sav=False)
    s = Project.open(snap)
    with pyazr.Session(snap, cwd=os.path.dirname(snap)) as b:
        chi_back = b.evaluate(b.parameters.x0).chi2
    check("(f) the saved file gives the fit's chi2", abs(chi_back - chi) <= 1e-6 * abs(chi), (chi, chi_back))
    check("(f) a fixed level stays fixed", all(lv.fixed for lv in s.levels() if (lv.group, lv.level) in fixed_levels)
          and len(fixed_levels) > 0)
    zero = [lv for lv in Project.open(SRC).levels() if all(c.width == 0.0 for c in lv.channels)]
    check("(f) a level whose widths are all zero keeps them",
          zero and all(all(c.width == 0.0 for c in s.level(lv.label).channels) for lv in zero), [z.label for z in zero])
finally:
    shutil.rmtree(work, ignore_errors=True)

print("(g) <thm>")
t = Project.open(THM)
th = t.thm()
check("(g) thm() reads the options and the experiment",
      th["options"]["kinematics"] == "lacognata" and th["experiments"]["LC10"]["Ebeam"] == "54", th)
before = t.to_bytes()
t.set_thm(kinematics="lacognata")
check("(g) a no-op option keeps every byte", t.to_bytes() == before)
t.set_thm(kinematics="triple", spectatorEnergy=0.1)
t.set_thm_experiment("LC10", Ebeam=55)
th = t.thm()
check("(g) the option and the experiment key written",
      th["options"]["kinematics"] == "triple" and th["options"]["spectatorEnergy"] == "0.1"
      and th["experiments"]["LC10"]["Ebeam"] == "55" and th["experiments"]["LC10"]["beam"] == "18O", th)
t.set_thm(spectatorEnergy=None)
check("(g) an option removed", "spectatorEnergy" not in t.thm()["options"])
for what, call, text in (
        ("a misspelt option", lambda: t.set_thm(kinemtics="triple"), "not understood"),
        ("a refused experiment key", lambda: t.set_thm_experiment("LC10", lineshape="maybe"), "lineshape"),
        ("segments=", lambda: t.set_thm_experiment("LC10", segments="1"), "segments="),
        ("removing an experiment a data set names", lambda: t.remove_thm_experiment("LC10"), "LC10")):
    before = t.to_bytes()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            call()
        check(f"(g) {what} refused", False, "no error")
    except FormatError as e:
        check(f"(g) {what} refused, the project unchanged", t.to_bytes() == before and text in str(e), str(e)[:200])
        if "experiment key" in what or what == "segments=":
            check(f"(g) {what}: the experiment named once ('{str(e)[:40]}...')",
                  str(e).startswith("<thm> experiment[LC10]: ") and str(e).count("experiment[LC10]") == 1, str(e)[:200])

print()
if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("all checks passed")
