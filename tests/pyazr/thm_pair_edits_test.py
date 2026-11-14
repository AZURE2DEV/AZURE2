#!/usr/bin/env python3
"""The THM columns of <levels> as Project edits: Project.set_pair(key, binding_energy=...) (column 32, B of the
Trojan horse, on every row of the pair) and Project.set_channel(level, pair, l, s, width_is_amplitude=...) (column 33,
gammaIsRWA, on the channel's row).

On two THM projects (tests/upgrade/T2, examples/o18_lacognata2008):
  (a) the edited project read by the engine gives the chi-squared of the same edit made by hand in the text (the
      columns written on the rows of the pair, or of the channel), exactly, and a chi-squared other than the
      unedited project's (the edit reaches the calculation); pairs() reads the new B back;
  (b) B = 0 on rows that do not write the THM columns leaves them as they are (no byte changes), and setting the
      value a row has is no change;
  (c) refusals, the project unchanged: B on a photon pair, a negative or non-finite B, the amplitude flag on a
      photon channel, a flag that is not a bool.

Run from anywhere:  python3 tests/pyazr/thm_pair_edits_test.py
"""
import math
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

PROJECTS = [os.path.join(ROOT, "tests", "upgrade", "T2", "T2.azr"),
            os.path.join(ROOT, "examples", "o18_lacognata2008", "o18_lacognata2008.azr")]
T1 = os.path.join(ROOT, "tests", "upgrade", "T1", "T1.azr")


def chi2(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pyazr.Session(path) as s:
            return s.evaluate(s.parameters.x0).chi2


def hand_edit(text, rows, column, token):
    """The <levels> rows chosen by ``rows(tokens)`` with ``column`` (0-based: 31 B, 32 gammaIsRWA) set to ``token``,
    the THM pair of columns appended ("0" for the other one) where a row has 31 or 32 tokens."""
    head, body, tail = re.split(r"(?<=<levels>\n)|(?=</levels>)", text, maxsplit=2)
    out = []
    for line in body.split("\n"):
        t = line.split()
        if len(t) >= 31 and rows(t):
            e1 = t[31:] if len(t) in (32, 34) else []
            e1 = e1[-1:] if e1 else []
            thm = t[31:33] if len(t) >= 33 else ["0", "0"]
            thm[column - 31] = token
            line = " ".join(t[:31] + thm + e1)
        out.append(line)
    return head + "\n".join(out) + tail


def save(project_text, directory, name):
    path = os.path.join(directory, name)
    with open(path, "w") as f:
        f.write(project_text)
    return path


for src in PROJECTS:
    tag = os.path.basename(src)
    print(f"(a) {tag}")
    base = Project.open(src)
    entrance = next(p for p in base.pairs() if p.binding_energy > 0)
    new_b = round(entrance.binding_energy + 0.15, 6)
    # A particle channel of the THM entrance pair with a width, for the amplitude flag.
    level, channel = next((lv, c) for lv in base.levels() if lv.active for c in lv.channels
                          if c.pair == entrance.key and c.width != 0 and not c.amplitude)
    with tempfile.TemporaryDirectory() as tmp:
        # A copy of the project's directory: data files are named relative to the project.
        d = os.path.join(tmp, "project")
        shutil.copytree(os.path.dirname(src), d, ignore=shutil.ignore_patterns("output*", "checks*"))
        src = os.path.join(d, tag)
        if True:
            p = Project.open(src)
            p.set_pair(entrance.key, binding_energy=new_b)
            check(f"{tag}: pairs() reads B back", next(q for q in p.pairs() if q.key == entrance.key).binding_energy
                  == new_b)
            by_call = p.save(os.path.join(d, "by_call.azr"))
            text = base.text
            by_hand = save(hand_edit(text, lambda t: int(t[5]) == entrance.key, 31, repr(new_b)), d, "by_hand.azr")
            c0, c1, c2 = chi2(src), chi2(by_call), chi2(by_hand)
            check(f"{tag}: B {entrance.binding_energy} -> {new_b} MeV: chi2 by the call == by hand ({c1!r})",
                  c1 == c2, f"{c1!r} {c2!r}")
            check(f"{tag}: the edit reaches the calculation (chi2 {c0:.6g} -> {c1:.6g})", c1 != c0)

            # The amplitude flag on one channel row, with a width in MeV^1/2 given in the same call.
            p = Project.open(src)
            amp = 0.05
            p.set_channel(level, channel.pair, channel.l, channel.s, width=amp, width_is_amplitude=True)
            row = p.level(level.label).channel(channel.pair, channel.l, channel.s)
            check(f"{tag}: the row reads as an amplitude of {amp}", row.amplitude and row.width == amp)
            by_call = p.save(by_call)
            line = channel.line
            lines = text.split("\n")
            toks = lines[line - 1].split()
            toks[11] = repr(amp)
            if len(toks) in (31, 32):
                toks = toks[:31] + ["0", "1"] + toks[31:]
            else:
                toks[32] = "1"
            lines[line - 1] = " ".join(toks)
            save("\n".join(lines), d, "by_hand.azr")
            c1, c2 = chi2(by_call), chi2(by_hand)
            check(f"{tag}: level {level.label}, channel ({channel.pair}, {channel.l}, {channel.s}) as an amplitude: "
                  f"chi2 by the call == by hand ({c1!r})", c1 == c2, f"{c1!r} {c2!r}")
            check(f"{tag}: the flag reaches the calculation", c1 != c0)

print("(b) no-op edits")
p = Project.open(T1)
before = p.to_bytes()
p.set_pair(1, binding_energy=0.0)
p.set_channel("1/2- 0.0", 1, 1, 0.5, width_is_amplitude=False)
check("B = 0 and the flag off on 31-column rows: no byte changes", p.to_bytes() == before)
q = Project.open(PROJECTS[0])
b = next(x for x in q.pairs() if x.binding_energy > 0)
before = q.to_bytes()
q.set_pair(b.key, binding_energy=b.binding_energy)
check("the B a pair has: no byte changes", q.to_bytes() == before)

print("(c) refusals")


def refused(what, call, *words, exc=FormatError):
    p = Project.open(T1)
    before = p.to_bytes()
    try:
        call(p)
        check(f"{what} refused", False, "no error")
    except exc as e:
        check(f"{what} refused, the project unchanged", p.to_bytes() == before and all(w in str(e) for w in words),
              str(e)[:300])


refused("B on a photon pair", lambda p: p.set_pair(2, binding_energy=2.2), "photon or beta pair", "particle pair")
refused("a negative B", lambda p: p.set_pair(1, binding_energy=-1.0), "-1", "out of range", ">= 0")
refused("a B that is not finite", lambda p: p.set_pair(1, binding_energy=math.nan), "out of range")
refused("the amplitude flag on a photon channel",
        lambda p: p.set_channel("3/2- 1.0", 2, 1, 0.5, width_is_amplitude=True), "photon or beta channel",
        "particle pair")
refused("a flag that is not a bool", lambda p: p.set_channel("1/2- 0.0", 1, 1, 0.5, width_is_amplitude="yes"),
        "True or False", exc=TypeError)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok")
