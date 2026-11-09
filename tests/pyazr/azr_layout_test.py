#!/usr/bin/env python3
"""AzrModel keeps a file's own layout (docs/upgrade/project.md section 5, step 4).

On every project file under tests/ and on variants of tests/13N, tests/17O and tests/upgrade/T2 -- CRLF line
endings, no blank lines in <levels>, the <levels> lines interleaved (every other line first) -- :
  1. an unedited model writes the file back byte for byte (the former writer re-emitted <levels> in its canonical
     layout, so only files already in that layout round-tripped);
  2. a field edit (one level's energy) changes exactly the lines of that level's channels, wherever they stand;
  3. a structural edit (a level removed) re-emits <levels> in the canonical layout, and the line ending is kept.
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
    from pyazr import AzrModel
    from pyazr import _azure2  # noqa: F401
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(0)

work = tempfile.mkdtemp(prefix="azr_layout.")


def model_of(text, name):
    p = os.path.join(work, name)
    with open(p, "w", newline="") as f:
        f.write(text)
    return AzrModel.from_file(p)


# 1. every project file, unedited
n = 0
for p in sorted(glob.glob(os.path.join(ROOT, "tests", "**", "*.azr"), recursive=True)):
    if "/output/" in p:
        continue
    raw = open(p, newline="").read()
    try:
        m = AzrModel.from_file(p)
    except ValueError:
        continue  # a file AZURE2 refuses (tests of the refusal)
    n += 1
    check(f"{p}: unedited round trip", m.to_text() == raw)
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
        check(f"{tag}: unedited round trip", m.to_text() == text)
        lv = m.levels[0]
        old_lines = {c.to_line() for c in lv.channels}
        lv.set_energy(lv.energy + 0.001)
        new_text = m.to_text()
        nl = "\r\n" if vname == "CRLF" else "\n"
        before, after = text.split(nl), new_text.split(nl)
        changed = [i for i, (x, y) in enumerate(zip(before, after)) if x != y]
        check(f"{tag}: a field edit changes only that level's {len(lv.channels)} lines",
              len(before) == len(after) and len(changed) == len(lv.channels)
              and all(before[i] in old_lines or before[i].rstrip() in old_lines for i in changed),
              f"{len(changed)} lines changed")
        m2 = model_of(text, "w.azr")
        m2.levels.pop()
        out = m2.to_text()
        check(f"{tag}: a structural edit keeps the line ending", (nl == "\r\n") == ("\r\n" in out))
        check(f"{tag}: ... and drops the level", len(AzrModel.from_file(m2.write(os.path.join(work, "x.azr"))).levels)
              == len(m.levels) - 1)

print("PASS" if not failures else f"FAIL: {len(failures)}")
sys.exit(1 if failures else 0)
