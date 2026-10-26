#!/usr/bin/env python3
"""Repository consistency: the drift that no physics test sees.

Every check here is a way the repository has gone wrong without any test
failing -- a flag the parser accepts that --help and the documentation never
mention, a test file nobody registered, a check.sh without its executable bit
(tests/run_tests.sh skips those silently).  Pure standard library, reads the
source tree only, runs in well under a second.

    python3 tests/consistency/repo_consistency_test.py

Exits 0 when everything agrees, 1 with one line per problem otherwise.
"""

import os
import re
import stat
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
problems = []


def read(rel):
    with open(os.path.join(ROOT, rel), encoding="utf-8") as f:
        return f.read()


def function_body(source, signature):
    """The text of the C++ function whose definition starts with `signature`."""
    start = source.index(signature)
    brace = source.index("{", start)
    depth = 0
    for i in range(brace, len(source)):
        if source[i] == "{":
            depth += 1
        elif source[i] == "}":
            depth -= 1
            if depth == 0:
                return source[brace:i]
    raise ValueError("unbalanced braces after " + signature)


# 1. Every command-line flag the parser accepts is in --help and in the
#    command-line reference.  A flag that is only in the source is a flag only
#    the person who wrote it knows about.
main_cpp = read("src/AZURE2.cpp")
parser = function_body(main_cpp, "bool parseOptions(")
help_text = function_body(main_cpp, "void printHelp(")
cli_doc = read("docs/source/reference/command_line.rst")

flags = set(re.findall(r'"(--[a-z0-9-]+)=?"', parser))
# Handled in main() before parseOptions, so accepted although parseOptions
# never sees them.
handled_in_main = {"--help", "--build-coulomb-table"}
flags -= handled_in_main
# Removed: rejected with an error, documented as removed, not offered by --help.
removed = {"--use-api"}

for flag in sorted(flags):
    if flag not in cli_doc:
        problems.append(f"{flag}: parsed in src/AZURE2.cpp, missing from docs/source/reference/command_line.rst")
    if flag not in removed and flag not in help_text:
        problems.append(f"{flag}: parsed in src/AZURE2.cpp, missing from printHelp()")

for flag in sorted(set(re.findall(r"\\t(--[a-z0-9-]+)", help_text))):
    if flag not in flags | handled_in_main:
        problems.append(f"{flag}: offered by printHelp() but not accepted by parseOptions()")


# 2. Every test source is registered with CTest.  An unregistered test is never
#    run, by CI or anyone else, and looks exactly like a passing one.
def check_registered(directory, pattern):
    cmake = read(os.path.join(directory, "CMakeLists.txt"))
    for name in sorted(os.listdir(os.path.join(ROOT, directory))):
        if re.fullmatch(pattern, name) and name not in cmake:
            problems.append(f"{directory}/{name}: not registered in {directory}/CMakeLists.txt")


check_registered("tests/pyazr", r".+_test\.py")
check_registered("tests/reference", r".+_test\.cpp")
check_registered("tests/gui", r".+_test\.cpp")
check_registered("tests/consistency", r".+_test\.py")


# 3. Every regression directory is actually exercised by tests/run_tests.sh:
#    an .azr needs expected/chiSquared.out or a check.sh, and a check.sh must
#    be executable -- the runner tests -x and skips the rest without a word.
#    The mode is taken from git's index where available, since that is what a
#    fresh clone (and CI) gets.
def git_mode(rel):
    try:
        out = subprocess.run(["git", "ls-files", "-s", "--", rel], cwd=ROOT,
                             capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.split()[0] if out else None


for name in sorted(os.listdir(os.path.join(ROOT, "tests"))):
    d = os.path.join(ROOT, "tests", name)
    if not os.path.isdir(d):
        continue
    has_azr = any(f.endswith(".azr") for f in os.listdir(d))
    has_expected = os.path.isfile(os.path.join(d, "expected", "chiSquared.out"))
    check = os.path.join(d, "check.sh")
    if has_azr and not has_expected and not os.path.isfile(check):
        problems.append(f"tests/{name}: has an .azr but neither expected/chiSquared.out nor check.sh, so it never runs")
    if os.path.isfile(check):
        mode = git_mode(f"tests/{name}/check.sh")
        executable = (mode == "100755") if mode else bool(os.stat(check).st_mode & stat.S_IXUSR)
        if not executable:
            problems.append(f"tests/{name}/check.sh: not executable, so tests/run_tests.sh skips it (git update-index --chmod=+x)")


# 4. pyazr's version is stated once, as MAJOR.MINOR.PATCH (pyproject.toml
#    reads it from there).
init = read("pyazr/__init__.py")
versions = re.findall(r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", init, re.M)
if len(versions) != 1 or not re.fullmatch(r"\d+\.\d+\.\d+", versions[0]):
    problems.append(f"pyazr/__init__.py: expected one MAJOR.MINOR.PATCH __version__, found {versions}")


# 5. Every name pyazr exports exists.  A stale __all__ entry breaks
#    `from pyazr import *` only, so nothing else notices.
all_match = re.search(r"__all__\s*=\s*\[(.*?)\]", init, re.S)
if all_match:
    exported = re.findall(r"['\"]([A-Za-z_][A-Za-z0-9_]*)['\"]", all_match.group(1))
    imported = set(re.findall(r"\b([A-Za-z_][A-Za-z0-9_]*)\b", init[:all_match.start()]))
    for name in exported:
        if name not in imported:
            problems.append(f"pyazr.__all__: '{name}' is exported but never imported in pyazr/__init__.py")
    duplicates = sorted({n for n in exported if exported.count(n) > 1})
    if duplicates:
        problems.append(f"pyazr.__all__: listed twice: {', '.join(duplicates)}")


# 6. The segment-line keyword blocks are the same set in all three .azr
#    readers.  The GUI keeps a tail verbatim only for keywords it knows
#    (SegmentKeywords.h), and pyazr splits the tail by its own table, so a
#    keyword the engine learns and either front end does not is mangled by a
#    load and save there.
engine_kw = set()
for header in ("include/SegLine.h", "include/ExtrapLine.h"):
    engine_kw |= set(re.findall(r'keyword == "([a-z]+)"', read(header)))
gui_src = read("gui/include/SegmentKeywords.h")
gui_kw = set(re.findall(r'token\.compare\("([a-z]+)"', function_body(gui_src, "inline bool isKeyword(")))
gui_kw -= {"mixture"}   # retired block, still carried through an old file's load/save
py_match = re.search(r"_KEYWORDS\s*=\s*\{(.*?)\}", read("pyazr/datasets.py"), re.S)
py_kw = set(re.findall(r'"([a-z]+)"\s*:', py_match.group(1))) if py_match else set()
for name, kw in (("gui/include/SegmentKeywords.h", gui_kw), ("pyazr/datasets.py _KEYWORDS", py_kw)):
    for k in sorted(engine_kw - kw):
        problems.append(f"segment keyword '{k}': read by the engine (SegLine.h/ExtrapLine.h), unknown to {name}")
    for k in sorted(kw - engine_kw):
        problems.append(f"segment keyword '{k}': known to {name}, not read by the engine")
if not engine_kw:
    problems.append("could not find the segment keywords in include/SegLine.h; update this test")


# 7. The repository carries its licence.
if os.path.getsize(os.path.join(ROOT, "LICENSE")) < 1000:
    problems.append("LICENSE: empty or truncated (README and pyproject.toml declare GPL-3.0)")


if problems:
    for p in problems:
        print("FAIL", p)
    print(f"{len(problems)} consistency problem(s)")
    sys.exit(1)
print(f"repository consistent: {len(flags)} CLI flags, segment keywords, test registration, check.sh modes, pyazr version and exports")
