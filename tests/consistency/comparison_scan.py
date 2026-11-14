#!/usr/bin/env python3
"""Comparisons in the shell tests whose failure cannot reach the caller (check 56 of repo_consistency_test.py).

A test that passes when its inputs differ is worse than none.  Two forms can be found by reading the awk programs
of tests/*/check.sh and tests/lib/*.sh:

  mark   an awk program prints a mark from its body ("x", "grid", ...: the two files' grids differ) and then, from
         END, a number as well, so the caller reads both lines and never sees the mark alone.  END must test the flag
         the body set together with the mark before it prints anything.
  bound  an upper bound (<, <=) on a value handed in with -v that the program never tests to be a number
         (`v ~ /^[-+]?[0-9]/...`): awk compares an empty value or a mark as text ("" < "1e-06") or, after `v + 0`,
         as 0, so a missing number or a mismatch mark passes.

    python3 tests/consistency/comparison_scan.py          # lists every place on the tree (path:line: kind: text)

Exits 1 when it lists anything.  Pure standard library.
"""

import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# awk [-F x] [-v name=value]... 'program': the options and the single-quoted program.
_AWK = re.compile(r"""\b(?:g|m)?awk((?:\s+-F\s*(?:'[^']*'|"[^"]*"|\S+)|\s+-v\s*[A-Za-z_]\w*=(?:"(?:\\.|[^"\\])*"|'[^']*'|\S+))*)"""
                  r"""\s+'([^']*)'""")
_VAR = re.compile(r"-v\s*([A-Za-z_]\w*)=")
_END = re.compile(r"\bEND\s*\{")
_MARK = re.compile(r"""print\s+"[A-Za-z_][\w ]*"\s*;\s*([A-Za-z_]\w*)\s*=\s*1\b""")
_PRINT = re.compile(r"\bprintf?\b")
# an upper bound whose left side is a name (or an element of an array), bare or as `name + 0`
_UPPER = re.compile(r"(?<![\w\]\)$.])([A-Za-z_]\w*)(\[\w+\])?\s*(?:\+\s*0\s*)?<=?(?!=)")


def awk_programs(text):
    """(line number, -v names, program) of every awk call of a shell script with a single-quoted program."""
    for m in _AWK.finditer(text):
        yield text.count("\n", 0, m.start()) + 1, set(_VAR.findall(m.group(1))), m.group(2)


def _split_names(program, names):
    """The arrays split() fills from a -v value: `split(s, x, ...)` makes x[1] a -v value too."""
    return {a for s, a in re.findall(r"\bsplit\(\s*([A-Za-z_]\w*)\s*,\s*([A-Za-z_]\w*)", program) if s in names}


def problems_of_program(names, program):
    """(kind, detail) of each form above in one awk program with its -v names."""
    found = []
    end = _END.search(program)
    if end:
        body, tail = program[:end.start()], program[end.end():]
        marks = _MARK.findall(body)
        unflagged = re.search(r"""print\s+"[A-Za-z_][\w ]*"\s*;\s*(?:next|exit)\b""", body)
        first_print = _PRINT.search(tail)
        if first_print and (marks or unflagged):
            before = tail[:first_print.start()]
            if not marks or not any(re.search(r"\b%s\b" % f, before) for f in marks):
                found.append(("mark", "the body prints a mark and END then prints as well"))
    if names:
        arrays = _split_names(program, names)
        for m in _UPPER.finditer(program):
            name, index = m.group(1), m.group(2) or ""
            if not (name in names or (index and name in arrays)):
                continue
            ref = re.escape(name + index)
            rest = program[m.end():]
            if re.match(r"\s*0\s*\?", rest) or re.match(r"\s*0\s*\)\s*%s\s*=\s*-" % ref, rest):
                continue  # an absolute value (`v < 0 ? -v : v`, `if (v < 0) v = -v`), not a bound
            if not re.search(ref + r"\s*!?~", program):
                found.append(("bound", f"{name + index} has an upper bound and is never tested to be a number"))
    return found


def scan_text(text):
    """(line number, kind, detail) of every problem in a shell script."""
    return [(line, kind, detail) for line, names, program in awk_programs(text)
            for kind, detail in problems_of_program(names, program)]


def scanned_files(root=ROOT):
    """The shell scripts this scan reads: tests/*/check.sh and tests/lib/*.sh (relative paths)."""
    tests = os.path.join(root, "tests")
    out = []
    for name in sorted(os.listdir(tests)):
        if os.path.isfile(os.path.join(tests, name, "check.sh")):
            out.append(f"tests/{name}/check.sh")
    lib = os.path.join(tests, "lib")
    out += [f"tests/lib/{n}" for n in sorted(os.listdir(lib)) if n.endswith(".sh")]
    return out


def scan_tree(root=ROOT):
    """(path, line, kind, detail) over every scanned file."""
    out = []
    for rel in scanned_files(root):
        with open(os.path.join(root, rel), encoding="utf-8") as f:
            out += [(rel, line, kind, detail) for line, kind, detail in scan_text(f.read())]
    return out


# Forms the scan must find (one each) and safe forms it must leave alone; check 56 runs them, so a change of the
# patterns that blinds the scan fails the consistency test.
FINDS = [
    """awk '{ if ($1 != $3) { print "grid"; e = 1; next } w = $2 } END { printf "%.3e\\n", w }'""",
    """awk '{ if (e) next; if ($1 != $3) { print "ENERGY"; e = 1; next } } END { printf "%.2e %d", m, n; exit e }'""",
    """awk '{ if ($1 != $3) { print "x"; next } } END { print w }'""",
    """awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t + 0) }'""",
    """awk -v d="$d" 'BEGIN { exit !(d < 1e-12) }'""",
    """awk -v s="$s" 'BEGIN { split(s, x, " "); exit !(x[2] > 0 && x[1] + 0 <= 1e-9) }'""",
    """awk -v d="$d" 'BEGIN { exit !(d != "mismatch" && d + 0 <= 1e-10) }'""",
]
SAFE = [
    """awk '{ if (e) next; if ($1 != $3) { print "grid"; e = 1; next } w = $2 } END { if (e) exit; printf "%.3e\\n", w }'""",
    """awk '{ if ($1 != $3) bad = 1 } END { print bad ? "x" : w }'""",
    """awk -v r="$1" -v t="$2" 'BEGIN { exit !(r ~ /^[0-9.]+$/ && r + 0 <= t + 0) }'""",
    """awk -v a="$1" -v b="$2" -v re="$num" 'BEGIN { exit !(a ~ re && a + 0 <= b + 0) }'""",
    """awk -v n="$n" 'BEGIN { exit !(n + 0 > 0) }'""",
    """awk '{ if ($1 + 0 < 3) print }'""",
    """awk -v a="$1" 'BEGIN { if (a !~ /^[0-9]/) exit 1; m = a < 0 ? -a : a; if (a < 0) a = -a }'""",
    """awk -v s="$s" 'BEGIN { split(s, x, " "); exit !(x[1] ~ /^[0-9]/ && x[1] + 0 <= 1e-9) }'""",
]


def self_test():
    """The forms of FINDS that the scan misses and the forms of SAFE it flags (both empty when it works)."""
    return [s for s in FINDS if len(scan_text(s)) != 1], [s for s in SAFE if scan_text(s)]


if __name__ == "__main__":
    missed, flagged = self_test()
    if missed or flagged:
        print("the scan misses", missed, "and flags", flagged)
        sys.exit(2)
    found = scan_tree()
    for rel, line, kind, detail in found:
        print(f"{rel}:{line}: {kind}: {detail}")
    print(f"{len(found)} place(s) in {len(scanned_files())} scripts")
    sys.exit(1 if found else 0)
