"""rmfit's sources stay valid on Python < 3.12 (pyazr declares >= 3.8 on dev, >= 3.10 on feature/upgrade).

Python 3.12 (PEP 701) accepts an f-string whose replacement field reuses the string's own quote or holds a backslash,
f"{d['key']}" inside f'...' for instance; 3.11 and older refuse the whole module with a SyntaxError.  Such a line
passed every test on the 3.12 cluster and broke `import rmfit.campaign` and `rmfit.review` on CI's Python 3.10
(2026-11-17).  On Python >= 3.12 this test reads the f-string tokens; on an older Python, compiling the files is the
test.
"""
import io
import os
import sys
import tokenize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def sources():
    for d, subs, files in os.walk(ROOT):
        subs[:] = [s for s in subs if not s.startswith((".", "_queue", "__pycache__")) and s != "attic"]
        for f in files:
            if f.endswith(".py"):
                yield os.path.join(d, f)


def offending(path):
    """Lines of `path` with a replacement field that reuses its f-string's quote or holds a backslash."""
    src = open(path, encoding="utf-8").read()
    if not hasattr(tokenize, "FSTRING_START"):
        try:
            compile(src, path, "exec", dont_inherit=True)
        except SyntaxError as err:
            return [err.lineno]
        return []
    out, stack = [], []
    for t in tokenize.generate_tokens(io.StringIO(src).readline):
        if t.type == tokenize.FSTRING_START:
            if stack and stack[-1][1] > 0 and stack[-1][0] in t.string:
                out.append(t.start[0])
            stack.append([t.string.lstrip("fFrRbB")[0], 0])
        elif t.type == tokenize.FSTRING_END:
            stack.pop()
        elif stack and t.type == tokenize.OP and t.string in "{}":
            stack[-1][1] += 1 if t.string == "{" else -1
        elif stack and stack[-1][1] > 0 and ((t.type == tokenize.STRING and stack[-1][0] in t.string) or "\\" in t.string):
            out.append(t.start[0])
    return sorted(set(out))


def main():
    bad = []
    for path in sorted(sources()):
        bad += [f"{os.path.relpath(path, ROOT)}:{n}" for n in offending(path)]
    # the detector itself: a nested same-quote field is caught, a different-quote one is not
    probe = os.path.join(HERE, "_py_compat_probe.tmp")
    try:
        with open(probe, "w") as fh:
            fh.write("d = {'a': 1}\nx = f\"{d['a']}\"\n")
        clean = offending(probe)
        with open(probe, "w") as fh:
            fh.write("d = {'a': 1}\nx = f'{d['a']}'\n")
        caught = offending(probe) if hasattr(tokenize, "FSTRING_START") else [2]
    finally:
        os.remove(probe)
    ok = not bad and clean == [] and caught == [2]
    for b in bad:
        print("  FAIL ", b, "-- f-string field reuses its quote or holds a backslash (Python < 3.12 refuses the module)")
    print(f"  {'ok  ' if clean == [] and caught == [2] else 'FAIL'}  detector: clean probe {clean}, nested probe {caught}")
    print("ALL PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
