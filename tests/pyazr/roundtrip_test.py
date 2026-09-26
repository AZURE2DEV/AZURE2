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
n = len(azrfile.AzrChannel(tokens).to_line().split())
# 31 fields for a conventional line, up to 33 with the optional THM columns.
check("emits a full line",
      azrfile._NFIELDS <= n <= azrfile._NFIELDS_MAX, f"{n} fields")

print("\n7. a <thm> options block survives edits wherever it sits")
# The engine reads <thm> from anywhere in the file (Config::ReadThmBlock);
# AzrModel interprets only <levels>, so the block has to come back verbatim.
thm = "<thm>\nkinematics=kf3body   # data / full three-body KF\nspectatorEnergy[5]=0.5\n</thm>\n"
src7 = open(files[0]).read()
at_end = src7 + "\n" + thm
before_levels = src7.replace("<levels>", thm + "<levels>", 1)
import tempfile
for label, text in (("at end of file", at_end), ("before <levels>", before_levels)):
    fd, path = tempfile.mkstemp(suffix=".azr")
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    try:
        m = AzrModel.from_file(path)
        check(f"{label}: byte-identical", m.to_text() == text)
        m.levels[0].channels[0].gamma = 1.75
        m.clear_data_segments()
        m.clear_target_effects()
        check(f"{label}: block intact after edits", thm in m.to_text())
    finally:
        os.remove(path)

print("\n8. every edit keeps the file's final newline -- or its absence")
# "\n".join(text.splitlines()) drops it, and every splicer of the blocks after
# <levels> did exactly that.
src8 = open(files[0]).read().rstrip("\n")
edits = [
    ("clear_data_segments", lambda m: m.clear_data_segments()),
    ("clear_target_effects", lambda m: m.clear_target_effects()),
    ("clear_extrapolations", lambda m: m.clear_extrapolations()),
    ("add_extrapolation", lambda m: m.add_extrapolation(1, 2, 0.1, 1.0, 0.1)),
    ("set_thm_option", lambda m: m.set_thm_option("vertex", "onshell")),
]
for ending in ("\n", ""):
    for name, edit in edits:
        fd, path = tempfile.mkstemp(suffix=".azr")
        with os.fdopen(fd, "w") as fh:
            fh.write(src8 + ending)
        try:
            m = AzrModel.from_file(path)
            edit(m)
            text = m.to_text()
            ok = text.endswith("\n") == (ending == "\n") and not text.endswith("\n\n")
            check(f"{name}, file {'with' if ending else 'without'} final newline", ok,
                  repr(text[-30:]))
        finally:
            os.remove(path)

print("\n9. the <thm> block: read, set, clear -- the GUI's rules")
# tests/7Li_p_a has a THM segment (line 1 of <segmentsData>) and no block.
thm_src = os.path.join(ROOT, "tests", "7Li_p_a", "7Li_p_a.azr")
tmpdir = tempfile.mkdtemp()


def project(text):
    path = os.path.join(tmpdir, "p.azr")
    with open(path, "w") as fh:
        fh.write(text)
    return AzrModel.from_file(path)


def refuses(what, fn):
    try:
        fn()
    except ValueError:
        return True
    print(f"        ({what}: no ValueError)")
    return False


base = open(thm_src).read()
m = project(base)
check("no block: defaults", m.thm_options() == {} and m.thm_options(defaults=True) == {
    "entranceL": "incoherent", "vertex": "constant", "kinematics": "lacognata",
    "coulombIntegral": False, "spectatorEnergy": 0.0})
m.set_thm_option("vertex", "constant")
check("setting a default adds no block", m.to_text() == base)
m.set_thm_option("kinematics", "kf3body")
m.set_thm_option("spectatorEnergy[5]", 0.5)
m.set_thm_option("coulombIntegral", True)
text = m.to_text()
check("new block after </targetInt>",
      "</targetInt>\n<thm>\nkinematics=kf3body\nspectatorEnergy[5]=0.5\ncoulombIntegral=1\n</thm>\n"
      in text, text[-200:])
check("reads back", m.thm_options() == {"kinematics": "kf3body", "coulombIntegral": True,
                                         "spectatorEnergy[5]": 0.5})
m2 = project(text)
check("written block parses the same", m2.thm_options() == m.thm_options())
for key in ("kinematics", "spectatorEnergy[5]", "coulombIntegral"):
    m2.clear_thm_option(key)
check("clearing every option removes the block", m2.to_text() == base)

# A hand-written block: comments, indentation, inline comments, aliases.
hand = ("<thm>\n# THM options of the refit\n  vertex = real   # alias of perlevel\n"
        "kinematics=triple\n\nspectatorEnergy[5]=0.50\n</thm>\n")
m = project(base + "\n" + hand)
check("aliases and spacing read as the engine does",
      m.thm_options() == {"vertex": "perlevel", "kinematics": "triple",
                          "spectatorEnergy[5]": 0.5})
m.set_thm_option("spectatorEnergy[5]", 0.5)
check("an option set to its value keeps the block verbatim", m.to_text() == base + "\n" + hand)
m.set_thm_option("vertex", "onshell")
m.set_thm_option("kinematics", "lacognata")
m.set_thm_option("entranceL", "coherent")
want = ("<thm>\n# THM options of the refit\n  vertex=onshell   # alias of perlevel\n"
        "\nspectatorEnergy[5]=0.50\nentranceL=coherent\n</thm>\n")
check("changed line rewritten in place, default removed, new key appended, "
      "comments kept", m.to_text() == base + "\n" + want, m.to_text()[-160:])
for key in ("vertex", "spectatorEnergy[5]", "entranceL"):
    m.clear_thm_option(key)
check("all default: block removed, comments included", m.to_text() == base + "\n")

# Refusals, each leaving the model untouched.
m = project(base)
bad = [
    ("unknown key", lambda: m.set_thm_option("vertexx", "onshell")),
    ("unknown vertex", lambda: m.set_thm_option("vertex", "offshell")),
    ("unknown kinematics", lambda: m.set_thm_option("kinematics", "kf")),
    ("coulombIntegral=maybe", lambda: m.set_thm_option("coulombIntegral", "maybe")),
    ("negative spectator energy", lambda: m.set_thm_option("spectatorEnergy", -0.1)),
    ("spectator energy not a number", lambda: m.set_thm_option("spectatorEnergy", "x")),
    ("bad pair key", lambda: m.set_thm_option("spectatorEnergy[p]", 0.1)),
    ("weight on a missing file", lambda: m.set_thm_weight(1, "no_such_file.dat")),
    ("weight on a segment beyond the block", lambda: m.set_thm_weight(2, "w.dat")),
    ("weight with '#' in the path", lambda: m.set_thm_option("weight[1]", "a#b.dat")),
    ("weight[0]", lambda: m.set_thm_option("weight[0]", "w.dat")),
]
for what, fn in bad:
    check(f"refused: {what}", refuses(what, fn))
check("refusals leave the model unchanged", m.to_text() == base)
check("clear of an unknown key refused", refuses("clear", lambda: m.clear_thm_option("foo")))


def table(name, rows):
    with open(os.path.join(tmpdir, name), "w") as fh:
        fh.write(rows)
    return name


for what, rows in (("one row", "1 1\n"), ("w = 0", "0 1\n1 0\n"),
                   ("E not increasing", "1 1\n1 2\n"), ("three columns", "0 1 2\n1 1 1\n"),
                   ("not a number", "0 1\nx 2\n")):
    name = table("bad.dat", rows)
    check(f"weight table refused: {what}", refuses(what, lambda: m.set_thm_weight(1, name)))
name = table("w.dat", "# E w\n-1.0 0.8\n10.0 3   # ramp w = 1 + E/5\n")
m.set_thm_weight(1, name)
check("weight set (relative to the .azr)", m.thm_options() == {"weight[1]": "w.dat"})
check("weight line written", "<thm>\nweight[1]=w.dat\n</thm>" in m.to_text())
m.set_thm_weight(1, os.path.join(tmpdir, name))
check("absolute weight path", m.thm_options()["weight[1]"] == os.path.join(tmpdir, name))
m.clear_thm_weight(1)
check("weight cleared", m.to_text() == base)

# A non-THM segment: the same data as an ordinary segment (isDiff 0).
line1 = [ln for ln in base.splitlines() if ln.startswith("1  5  4  0  8.2")][0]
nonthm = base.replace(line1, line1 + "\n" + line1.replace("  180  10  ", "  180  0  "), 1)
m = project(nonthm)
check("weight on a non-THM segment refused", refuses("non-THM", lambda: m.set_thm_weight(2, name)))
m.set_thm_weight(1, name)
check("weight on the THM segment accepted", m.thm_options() == {"weight[1]": "w.dat"})

# A block the engine refuses is refused here too.
m = project(base + "\n<thm>\nvertex=offshell\n</thm>\n")
check("a block with a line the engine refuses: thm_options raises",
      refuses("read", lambda: m.thm_options()))
check("... and set_thm_option raises", refuses("set", lambda: m.set_thm_option("vertex", "onshell")))
m = project(base + "\n<thm>\nvertex=onshell\n")
check("an unterminated block raises", refuses("unterminated", lambda: m.thm_options()))

# The engine, if built, accepts what was written and refuses what was refused.
import shutil
import subprocess
binaries = [b for b in glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2")) if os.access(b, os.X_OK)]
if binaries:
    binary = max(binaries, key=os.path.getmtime)
    shutil.copytree(os.path.join(ROOT, "tests", "7Li_p_a", "data"), os.path.join(tmpdir, "data"))
    os.makedirs(os.path.join(tmpdir, "output"), exist_ok=True)
    os.makedirs(os.path.join(tmpdir, "checks"), exist_ok=True)

    def engine(model):
        path = model.write(os.path.join(tmpdir, "run.azr"))
        r = subprocess.run([binary, "--no-gui", "--no-readline", "run.azr"], cwd=tmpdir,
                           input="1\n\n\n7\n", capture_output=True, text=True, timeout=600,
                           env=dict(os.environ, OMP_NUM_THREADS="1"))
        chi = [ln for ln in r.stdout.splitlines() if ln.startswith("Total Chi-Squared:")]
        return r.returncode, (float(chi[0].split()[-1]) if chi else None), r.stdout

    m = project(base)
    m.set_thm_option("vertex", "onshell")
    rc, chi, _ = engine(m)
    check(f"engine runs the written block (vertex=onshell, pin 2070.84)",
          rc == 0 and chi is not None and abs(chi - 2070.84) < 0.01, f"rc {rc} chi2 {chi}")
    m.set_thm_weight(1, name)
    rc, chi2w, _ = engine(m)
    check("engine accepts the written weight (w = 1 + E/5 moves chi2)",
          rc == 0 and chi2w is not None and abs(chi2w - chi) > 1e-3, f"rc {rc} chi2 {chi2w}")
    rc, _, out = engine(project(base + "\n<thm>\nvertex=offshell\n</thm>\n"))
    check("engine refuses what pyazr refuses", rc != 0 and "ERROR: <thm>" in out, f"rc {rc}")
else:
    print("  skip  engine checks (no build*/src/AZURE2)")
shutil.rmtree(tmpdir)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all round-trip checks passed")
