#!/usr/bin/env python3
"""AzrModel must return a .azr byte for byte unless something was edited.

A project file is hand-maintained and column-aligned, and people diff it. If
reading and writing it back reflows every line of <levels>, every edit made
through pyazr shows up as a whole-block change and the real one-field change is
invisible -- so the round-trip has to be exact, not merely equivalent.

Imports pyazr.azrfile directly rather than the package, so this runs without the
compiled engine module: AzrModel is pure Python and parses the file itself.

Run from anywhere:  python3 tests/pyazr/roundtrip_test.py
"""
import glob
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

spec = importlib.util.spec_from_file_location(
    "azrfile", os.path.join(ROOT, "pyazr", "azrfile.py"))
azrfile = importlib.util.module_from_spec(spec)
spec.loader.exec_module(azrfile)
AzrModel = azrfile.AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


print("1. every project file in tests/ round-trips unchanged")
files = sorted(glob.glob(os.path.join(ROOT, "tests", "*", "*.azr")))
if not files:
    check("found project files", False, "no .azr under tests/")
for f in files:
    src = open(f).read()
    check(os.path.relpath(f, ROOT), AzrModel.from_file(f).to_text() == src)

print("\n2. an edit changes that field and nothing else")
model = AzrModel.from_file(files[0])
before = [c.to_line() for lv in model.levels for c in lv.channels]
channel = model.levels[0].channels[0]
original_gamma = channel.gamma
channel.gamma = 1.75
after = [c.to_line() for lv in model.levels for c in lv.channels]
changed = [i for i, (a, b) in enumerate(zip(before, after)) if a != b]
check("exactly one line changed", changed == [0], f"changed {changed}")
check("only the gamma token differs",
      [t for i, t in enumerate(before[0].split()) if t != after[0].split()[i]]
      == [azrfile._fmt(original_gamma)])
check("value reads back", channel.gamma == 1.75)

print("\n3. a value that fits keeps the column width")
check("line length unchanged", len(after[0]) == len(before[0]),
      f"{len(before[0])} -> {len(after[0])}")

print("\n4. setting a field to what it already holds is not an edit")
model2 = AzrModel.from_file(files[0])
level = model2.levels[0]
level.set_energy(level.energy)
level.set_fixed(level.fixed)
check("still byte-identical", model2.to_text() == open(files[0]).read())

print("\n5. renumbering on every write does not disturb the lines")
model3 = AzrModel.from_file(files[0])
model3.to_text()
check("second write identical to the first", model3.to_text() == open(files[0]).read())

print("\n6. a channel built from tokens, with no raw line, still emits")
tokens = model.levels[0].channels[0].tokens
check("emits a full line", len(azrfile.AzrChannel(tokens).to_line().split())
      == azrfile._NFIELDS)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all round-trip checks passed")
