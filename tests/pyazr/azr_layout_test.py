#!/usr/bin/env python3
"""pyazr.Project keeps a file's own layout (the lossless document, io::Document; the engine's row edits).

On every format-3 project file under tests/ and on variants of tests/13N, tests/17O and tests/upgrade/T2 -- CRLF
line endings, no blank lines in <levels>, the <levels> rows interleaved (every other row first) -- :
  1. an unedited project gives the file back byte for byte;
  2. a field edit (one level's energy, Project.set_level) changes exactly the rows of that level, wherever they stand;
  3. a structural edit (the last level removed, Project.remove_level) keeps the line ending and drops the level.
Needs the compiled engine; skips cleanly without it.  Run from anywhere:  python3 tests/pyazr/azr_layout_test.py
"""
import glob
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
failures = []


def check(name, ok, detail=""):
    if not ok:
        print(f"  FAIL  {name}  -- {detail}")
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    from pyazr import Project, FormatError
    from pyazr import _azure2  # noqa: F401
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

work = tempfile.mkdtemp(prefix="azr_layout.")


def model_of(text, name):
    p = os.path.join(work, name)
    with open(p, "w", newline="") as f:
        f.write(text)
    return Project.open(p)


# 1. every project file, unedited
n = 0
for p in sorted(glob.glob(os.path.join(ROOT, "tests", "**", "*.azr"), recursive=True)):
    if "/output/" in p:
        continue
    raw = open(p, newline="").read()
    try:
        m = Project.open(p)
    except (FormatError, ValueError):
        continue  # a file the engine refuses (tests of the refusal; the originals kept for the converter)
    n += 1
    check(f"{p}: unedited round trip", m.to_bytes().decode() == raw)
print(f"  {n} project files: unedited round trip checked")

# 2, 3. layout variants
for rel in ("13N/13N.azr", "17O/17O.azr", "upgrade/T2/T2.azr"):
    src = open(os.path.join(ROOT, "tests", rel), newline="").read()
    lines = src.split("\n")
    a = next(i for i, l in enumerate(lines) if l.strip() == "<levels>")
    b = next(i for i, l in enumerate(lines) if l.strip() == "</levels>")
    body = [l for l in lines[a + 1:b] if l.strip()]
    variants = {
        "CRLF": src.replace("\n", "\r\n"),
        "no blank lines": "\n".join(lines[:a + 1] + body + lines[b:]),
        "interleaved": "\n".join(lines[:a + 1] + body[::2] + body[1::2] + lines[b:]),
    }
    for vname, text in variants.items():
        tag = f"{rel} ({vname})"
        m = model_of(text, "v.azr")
        check(f"{tag}: unedited round trip", m.to_bytes().decode() == text)
        lv = m.levels()[0]
        rows = {c.line - 1 for c in lv.channels}
        m.set_level(lv, energy=lv.energy + 0.001)
        new_text = m.to_bytes().decode()
        nl = "\r\n" if vname == "CRLF" else "\n"
        before, after = text.split(nl), new_text.split(nl)
        changed = {i for i, (x, y) in enumerate(zip(before, after)) if x != y}
        check(f"{tag}: a field edit changes only that level's {len(lv.channels)} rows",
              len(before) == len(after) and changed == rows, f"{sorted(changed)} vs {sorted(rows)}")
        m2 = model_of(text, "w.azr")
        m2.remove_level(m2.levels()[-1])
        out = m2.to_bytes().decode()
        check(f"{tag}: a structural edit keeps the line ending", (nl == "\r\n") == ("\r\n" in out))
        check(f"{tag}: ... and drops the level", len(Project.open(m2.save(os.path.join(work, "x.azr"))).levels())
              == len(m.levels()) - 1)

print("PASS" if not failures else f"FAIL: {len(failures)}")
sys.exit(1 if failures else 0)
