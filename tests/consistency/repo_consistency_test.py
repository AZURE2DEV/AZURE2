#!/usr/bin/env python3
"""Repository consistency: the drift that no physics test sees.

Every check here is a way the repository has gone wrong without any test
failing -- a flag the parser accepts that --help and the documentation never
mention, a test file nobody registered, a check.sh without its executable bit
(tests/run_tests.sh skips those silently).  Pure standard library, reads the
source tree only, runs in well under a second.

    python3 tests/consistency/repo_consistency_test.py
    python3 tests/consistency/repo_consistency_test.py --tighten-sizes   # lower the size record (check 26), only

Exits 0 when everything agrees, 1 with one line per problem otherwise.
"""

import glob
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


# 26. The size ratchet (CLAUDE.md "Code structure"): a file that can be read whole is one a person or an agent can
#     maintain.  Every hand-written source file (.cpp, .h, .py and the other C/C++ suffixes) of the engine, the
#     modules, the API, pyazr and tools/ that is longer than SIZE_LIMIT lines is listed in SIZE_RECORD with its line
#     count, and may not be longer than recorded; an unlisted file above the limit fails, so a new file stays below
#     it (entries are added by the director only).  Tests, external/ and generated files are not counted.
#     `--tighten-sizes` rewrites the record to the present counts where they are lower and drops the files that fell
#     to the limit or below; it never raises a number and never adds a file.  Exact equality is not required: other
#     branches shorten files too, and the record is tightened at each merge.  It runs first, because --tighten-sizes
#     runs it alone.
SIZE_LIMIT = 1000
SIZE_DIRS = ("src", "include", "physics", "io", "experiment", "inference", "api/src", "pyazr", "tools", "stopping")
SIZE_SUFFIXES = (".cpp", ".h", ".py", ".hpp", ".cc", ".c", ".H")
SIZE_SKIP_DIRS = {"tests", "__pycache__"}  # test sources, bytecode
SIZE_GENERATED = set()                     # generated files in SIZE_DIRS (none today), relative paths
SIZE_RECORD = "tests/consistency/file_sizes.tsv"
SIZE_HEADER = f"""# The size ratchet (tests/consistency/repo_consistency_test.py, check 26): every hand-written source file above
# {SIZE_LIMIT} lines, with the most lines it may have.  A listed file may not grow; an unlisted file may not pass {SIZE_LIMIT}
# lines.  Lower the numbers (never raise them) with
#   python3 tests/consistency/repo_consistency_test.py --tighten-sizes
# Entries are added by the director only.
"""


def line_count(rel):
    """Lines of a file as wc -l counts them, plus a last line without its newline."""
    with open(os.path.join(ROOT, rel), "rb") as f:
        data = f.read()
    return data.count(b"\n") + (1 if data and not data.endswith(b"\n") else 0)


def sized_files():
    """{relative path: lines} of every file the ratchet counts."""
    sizes = {}
    for top in SIZE_DIRS:
        for dirpath, dirnames, files in os.walk(os.path.join(ROOT, top)):
            dirnames[:] = sorted(d for d in dirnames if d not in SIZE_SKIP_DIRS and not d.startswith((".", "build")))
            for name in sorted(files):
                rel = os.path.relpath(os.path.join(dirpath, name), ROOT).replace(os.sep, "/")
                if name.endswith(SIZE_SUFFIXES) and rel not in SIZE_GENERATED:
                    sizes[rel] = line_count(rel)
    return sizes


def read_size_record():
    """{relative path: recorded lines}, in file order."""
    record = {}
    for number, line in enumerate(read(SIZE_RECORD).splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 2 or not fields[1].strip().isdigit():
            raise ValueError(f"{SIZE_RECORD}:{number}: not '<path><TAB><lines>': {line!r}")
        record[fields[0].strip()] = int(fields[1])
    return record


def check_sizes():
    try:
        record = read_size_record()
    except (OSError, ValueError) as e:
        problems.append(f"check 26: {e}")
        return
    for rel, lines in sorted(sized_files().items()):
        if rel in record and lines > record[rel]:
            problems.append(f"{rel}: {lines} lines, more than the {record[rel]} recorded in {SIZE_RECORD}; move code "
                            "out of it into a file of its own (CLAUDE.md, \"Code structure\")")
        elif rel not in record and lines > SIZE_LIMIT:
            problems.append(f"{rel}: {lines} lines, above the limit of {SIZE_LIMIT} and not in {SIZE_RECORD}; split "
                            "it (a new file stays below the limit; only the director adds an entry)")


def tighten_sizes():
    record = read_size_record()
    sizes = sized_files()
    kept = {}
    for rel, recorded in record.items():
        lines = sizes.get(rel, 0)
        if lines <= SIZE_LIMIT:
            print(f"dropped {rel}: {lines} lines, recorded {recorded}")
        elif lines < recorded:
            print(f"lowered {rel}: {recorded} -> {lines}")
            kept[rel] = lines
        else:
            kept[rel] = recorded
    with open(os.path.join(ROOT, SIZE_RECORD), "w", encoding="utf-8", newline="\n") as f:
        f.write(SIZE_HEADER)
        for rel in sorted(kept):
            f.write(f"{rel}\t{kept[rel]}\n")
    print(f"{SIZE_RECORD}: {len(kept)} files")


if sys.argv[1:] == ["--tighten-sizes"]:
    tighten_sizes()
    sys.exit(0)
if sys.argv[1:]:
    print(f"usage: {sys.argv[0]} [--tighten-sizes]")
    sys.exit(2)
check_sizes()


# 1. Every command-line flag the parser accepts is in --help and in the
#    command-line reference.  A flag that is only in the source is a flag only
#    the person who wrote it knows about.
main_cpp = read("src/AZURE2.cpp")
parser = function_body(main_cpp, "void parseOptions(")
help_text = function_body(main_cpp, "void printHelp(")
cli_doc = read("docs/source/reference/command_line.rst")

flags = set(re.findall(r'"(--[a-z0-9-]+)=?"', parser))
# The flags of the options of the calculation: their rows of io::ConfigKeys() (Option(name, flag, negation, ...)),
# parsed from the table (ApplyOptionFlag) and printed from it (printOptionFlags).
option_flags = {f for row in re.findall(r'Option\("[a-z_]+", "(--[a-z0-9-]*)?", "(--[a-z0-9-]*)?"',
                                        read("io/src/ConfigSection.cpp")) for f in row if f}
if not option_flags or "printOptionFlags();" not in help_text:
    problems.append("io/src/ConfigSection.cpp: no option flags found, or printHelp() does not call printOptionFlags()")
flags |= option_flags
# The fitters' switches and --max-calls: printed from the fitter table (FitterSwitchesHelp, src/Fitters.cpp; check 46).
fitter_help = function_body(read("src/Fitters.cpp"), "std::string FitterSwitchesHelp(") \
    if "FitterSwitchesHelp()" in help_text else ""
# Handled in main() before parseOptions, so accepted although parseOptions
# never sees them.
handled_in_main = {"--help", "--build-coulomb-table", "--convert", "--convert-params", "--schema", "--schema-out", "--validate",
                   "--describe", "--print-predictions",
                   # option of --convert-params
                   "--project",
                   # options of --convert and of --validate, read in main() with their command
                   "--to", "--json-report", "--answers", "--exact-grids", "--strict", "--allow-questions", "--json",
                   # option of --describe
                   "--record"}
flags -= handled_in_main
# Removed: rejected with an error, documented as removed, not offered by --help.
removed = {"--use-api", "--no-readline"}

for flag in sorted(flags):
    if flag not in cli_doc:
        problems.append(f"{flag}: parsed in src/AZURE2.cpp, missing from docs/source/reference/command_line.rst")
    if flag not in removed and flag not in option_flags and flag not in help_text and flag not in fitter_help:
        problems.append(f"{flag}: parsed in src/AZURE2.cpp, missing from printHelp()")

for flag in sorted(set(re.findall(r"\\t(--[a-z0-9-]+)", help_text))):
    if flag not in flags | handled_in_main:
        problems.append(f"{flag}: offered by printHelp() but not accepted by parseOptions()")

# 36. A value-taking flag reads its value through include/CommandLineValues.h (OptionWholeNumber, OptionText), which
#     refuses a missing or unparseable value with the usage line and exit 2; a value parsed by hand in parseOptions
#     (atoi, strtol, substr past the "--") was ignored with a warning or read as 0 (--max-calls N,
#     --refreeze-covariance=x).
VALUE_BY_HAND = re.compile(r"\b(atoi|atol|strtol|strtod|stoi|stol|stod|substr)\((?!0, 2\))")
for bad in VALUE_BY_HAND.findall(parser):
    problems.append(f"src/AZURE2.cpp parseOptions: a flag's value read with {bad}; read it with OptionWholeNumber or "
                    "OptionText (include/CommandLineValues.h, check 36)")
if not VALUE_BY_HAND.search("std::atoi(it->substr(22).c_str())") or VALUE_BY_HAND.search('arg.substr(0, 2) == "--"'):
    problems.append("check 36: the pattern no longer matches the form it guards; update this check")

# 2. Every test source is registered with CTest.  An unregistered test is never
#    run, by CI or anyone else, and looks exactly like a passing one.
def check_registered(directory, pattern, cmake_dir=None):
    cmake_dir = cmake_dir or directory
    cmake = read(os.path.join(cmake_dir, "CMakeLists.txt"))
    for name in sorted(os.listdir(os.path.join(ROOT, directory))):
        if re.fullmatch(pattern, name) and name not in cmake:
            problems.append(f"{directory}/{name}: not registered in {cmake_dir}/CMakeLists.txt")


check_registered("tests/pyazr", r".+_test\.py")
check_registered("tests/reference", r".+_test\.cpp")
check_registered("tests/consistency", r".+_test\.py")
check_registered("tests/upgrade", r".+_test\.cpp")
check_registered("physics/tests", r".+_test\.cpp", "physics")
check_registered("io/tests", r".+_test\.cpp", "io")
check_registered("experiment/tests", r".+_test\.cpp", "experiment")
check_registered("stopping/tests", r".+_test\.cpp", "stopping")
check_registered("inference/tests", r".+_test\.cpp", "inference")


# 3. Every regression directory is actually exercised by tests/run_tests.sh:
#    an .azr needs expected/chiSquared.out or a check.sh, and a check.sh must
#    be executable -- the runner tests -x and skips the rest without a word.
#    The mode is taken from git's index where available, since that is what a
#    fresh clone (and CI) gets.
def git_mode(rel):
    try:
        out = subprocess.run(["git", "-c", "safe.directory=*", "ls-files", "-s", "--", rel], cwd=ROOT,
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
        if mode:
            executable = mode == "100755"
        elif os.name == "nt" or sys.platform.startswith(("win", "msys", "cygwin")):
            executable = True   # no executable bit on Windows file systems and no git answer: cannot judge here
        else:
            executable = bool(os.stat(check).st_mode & stat.S_IXUSR)
        if not executable:
            problems.append(f"tests/{name}/check.sh: not executable, so tests/run_tests.sh skips it (git update-index --chmod=+x)")


# 4. pyazr's version is stated once, in pyazr/__init__.py (pyproject.toml reads it from there), as
#    MAJOR.MINOR.PATCH with an optional PEP 440 development suffix (.devN) while a major version is in work.
init = read("pyazr/__init__.py")
versions = re.findall(r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", init, re.M)
if len(versions) != 1 or not re.fullmatch(r"\d+\.\d+\.\d+(\.dev\d+)?", versions[0]):
    problems.append(f"pyazr/__init__.py: expected one MAJOR.MINOR.PATCH[.devN] __version__, found {versions}")


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
    # ... and the user guide lists exactly these names (its section "The names of the package").
    guide = read("docs/source/user_guide/pyazr.rst")
    block = re.search(r"\nThe names of the package\n-+\n.*?\.\. code-block:: text\n\n((?:   [^\n]*\n)+)", guide, re.S)
    if not block:
        problems.append("docs/source/user_guide/pyazr.rst: no section 'The names of the package' with a text block")
    else:
        listed_names = block.group(1).split()
        for name, where in ((n, "the user guide") for n in sorted(set(exported) - set(listed_names))):
            problems.append(f"pyazr.__all__: '{name}' is exported but not listed in {where} (The names of the package)")
        for name in sorted(set(listed_names) - set(exported)):
            problems.append(f"docs/source/user_guide/pyazr.rst: '{name}' is listed (The names of the package) but "
                            "not in pyazr.__all__")


# 7. The repository carries its licence.
if os.path.getsize(os.path.join(ROOT, "LICENSE")) < 1000:
    problems.append("LICENSE: empty or truncated (README and pyproject.toml declare GPL-3.0)")


# 8. physics/ uses nothing (docs/upgrade/ARCHITECTURE.md section 2): its sources include only its own headers
#    ("azure2/physics/...") and the C++ standard library (<vector>, <cmath>: lower-case names, no path, no ".h").
#    The ctest physics_standalone builds physics/ alone and catches an AZURE2 header or library; a header that
#    happens to be on the system include path (GSL, Qt, omp.h) compiles there, so it is caught here.
COULOMB_HEADERS = {n for _, _, fs in os.walk(os.path.join(ROOT, "physics", "coulomb")) for n in fs if n.endswith((".h", ".H"))}
for dirpath, _, files in os.walk(os.path.join(ROOT, "physics")):
    for name in sorted(files):
        if not name.endswith((".h", ".hpp", ".cpp")):
            continue
        rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
        for inc in re.findall(r'^\s*#\s*include\s*([<"][^>"]+[>"])', read(rel), re.M):
            if inc.startswith('"') and inc[1:-1].startswith("azure2/physics/"):
                continue
            # physics/coulomb/: the algorithms taken from pycoulomb include their own headers by bare name
            # (docs/upgrade/coulomb.md section 1.3); a name that is a file of physics/coulomb/ is allowed there.
            if rel.startswith(os.path.join("physics", "coulomb") + os.sep) and inc.startswith('"') and inc[1:-1] in COULOMB_HEADERS:
                continue
            # physics/tests/: a test may include a helper header of physics/tests/ itself by bare name (shared
            # readers of reference data, such as transfer_fresco_common.h).
            if (rel.startswith(os.path.join("physics", "tests") + os.sep) and inc.startswith('"') and "/" not in inc
                    and os.path.isfile(os.path.join(ROOT, "physics", "tests", inc[1:-1]))):
                continue
            if inc.startswith("<") and re.fullmatch(r"[a-z_]+", inc[1:-1]):
                continue
            problems.append(f"{rel}: includes {inc}; physics/ may include only azure2/physics/ headers and the C++ standard library")


# 8b. inference/ uses nothing (docs/upgrade/inference.md section 1): its own headers ("azure2/inference/...") and the
#     C++ standard library only; a test may include the helper header of inference/tests/ by bare name.  The ctest
#     inference_standalone builds the directory alone; a header on the system include path (Minuit2, GSL) compiles
#     there, so it is caught here.  An engine header (Objective.h, WorkingCopy.h) would make the library a renamed
#     directory of the engine.
def inference_include_problems(rel, text):
    found = []
    for inc in re.findall(r'^\s*#\s*include\s*([<"][^>"]+[>"])', text, re.M):
        if inc.startswith('"') and inc[1:-1].startswith("azure2/inference/"):
            continue
        if (rel.startswith("inference/tests/") and inc.startswith('"') and "/" not in inc
                and os.path.isfile(os.path.join(ROOT, "inference", "tests", inc[1:-1]))):
            continue
        if inc.startswith("<") and re.fullmatch(r"[a-z_]+", inc[1:-1]):
            continue
        found.append(f"{rel}: includes {inc}; inference/ may include only azure2/inference/ headers and the C++ "
                     "standard library (check 8b)")
    return found


for dirpath, _, files in os.walk(os.path.join(ROOT, "inference")):
    for name in sorted(files):
        if name.endswith((".h", ".hpp", ".cpp")):
            rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
            problems += inference_include_problems(rel, read(rel))
if not inference_include_problems("inference/src/X.cpp", '#include "Objective.h"\n') or \
        not inference_include_problems("inference/src/X.cpp", '#include "Minuit2/MnUserParameters.h"\n') or \
        not inference_include_problems("inference/src/X.cpp", '#include "azure2/physics/Constants.h"\n') or \
        inference_include_problems("inference/src/X.cpp", '#include "azure2/inference/ParameterKind.h"\n#include <vector>\n'):
    problems.append("check 8b: a planted engine, Minuit2 or physics/ include in inference/ is not caught; update this check")


# 9. io/ uses nothing either (docs/upgrade/project.md section 2): its own headers ("azure2/io/...") and the C++
#    standard library only.
for dirpath, _, files in os.walk(os.path.join(ROOT, "io")):
    for name in sorted(files):
        if not name.endswith((".h", ".hpp", ".cpp")):
            continue
        rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
        for inc in re.findall(r'^\s*#\s*include\s*([<"][^>"]+[>"])', read(rel), re.M):
            if inc.startswith('"') and inc[1:-1].startswith("azure2/io/"):
                continue
            if inc.startswith("<") and re.fullmatch(r"[a-z_]+", inc[1:-1]):
                continue
            problems.append(f"{rel}: includes {inc}; io/ may include only azure2/io/ headers and the C++ standard library")


# 32. stopping/ is a leaf on physics/ (docs/upgrade/stopping.md section 1): its headers ("azure2/stopping/..."),
#     physics/'s ("azure2/physics/..."), the environment record of support/ ("azure2/support/Environment.h", where
#     AZURE2_SRIM_DATA_DIR is read), the C++ standard library, and the bare names of its own files where they live
#     (src/'s generated-table header, the tests' check.h).
STOPPING_OWN = {"src": {"TableText.h"}, "tests": {"check.h"}}
for dirpath, _, files in os.walk(os.path.join(ROOT, "stopping")):
    for name in sorted(files):
        if not name.endswith((".h", ".hpp", ".cpp", ".in")):
            continue
        rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
        part = rel.split(os.sep)[1] if len(rel.split(os.sep)) > 2 else ""
        for inc in re.findall(r'^\s*#\s*include\s*([<"][^>"]+[>"])', read(rel), re.M):
            if inc.startswith('"') and (inc[1:-1].startswith(("azure2/stopping/", "azure2/physics/")) or
                                        inc[1:-1] == "azure2/support/Environment.h"):
                continue
            if inc.startswith('"') and inc[1:-1] in STOPPING_OWN.get(part, set()):
                continue
            if inc.startswith("<") and re.fullmatch(r"[a-z_]+", inc[1:-1]):
                continue
            problems.append(f"{rel}: includes {inc}; stopping/ may include only azure2/stopping/ and azure2/physics/ "
                            "headers, azure2/support/Environment.h, its own files and the C++ standard library "
                            "(check 32)")



# 9b. experiment/ uses physics/ and io/'s plain record structs only (docs/upgrade/data_model.md section 2.5, decision
#     M1; docs/upgrade/project_format3.md section 7.3, the director's decision of 2026-11-09): its own headers
#     ("azure2/experiment/..."), physics/'s ("azure2/physics/..."), the one header "azure2/io/Record.h" (header-only
#     structs) and the C++ standard library.  Its tests may also include io/'s reader ("azure2/io/..."): they read
#     their fixtures as the engine will; so may the files of the reading path, src/ProjectReader.cpp,
#     src/LevelTablesOfRows.cpp and src/LevelEdits.cpp (library azure2_experiment_reader, which links io/: records ->
#     Validate -> <levels> tables -> Build, step 3b; the <levels> rows and their edits), which the table library
#     azure2_experiment does not contain.  An engine header (Unit.h, Config.h), another io/ header in
#     the table library or a system library makes the module depend on what it must not.
for dirpath, _, files in os.walk(os.path.join(ROOT, "experiment")):
    for name in sorted(files):
        if not name.endswith((".h", ".hpp", ".cpp")):
            continue
        rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
        in_tests = rel.startswith(os.path.join("experiment", "tests") + os.sep) or \
            rel in (os.path.join("experiment", "src", "ProjectReader.cpp"),
                    os.path.join("experiment", "src", "LevelTablesOfRows.cpp"),
                    os.path.join("experiment", "src", "LevelEdits.cpp"))
        for inc in re.findall(r'^\s*#\s*include\s*([<"][^>"]+[>"])', read(rel), re.M):
            if inc.startswith('"') and inc[1:-1].startswith(("azure2/experiment/", "azure2/physics/")):
                continue
            if inc == '"azure2/io/Record.h"' or (in_tests and inc.startswith('"azure2/io/')):
                continue
            if inc.startswith("<") and re.fullmatch(r"[a-z_]+", inc[1:-1]):
                continue
            problems.append(f"{rel}: includes {inc}; experiment/ may include only azure2/experiment/, azure2/physics/ "
                            "and azure2/io/Record.h headers (its tests also azure2/io/) and the C++ standard library")


# 9c. convert/ (the old formats and the converter, docs/upgrade/project_format3.md section 8) uses io/, experiment/ and
#     physics/ and nothing else of AZURE2 (the <levels> line, whose grammar format 3 keeps, is io/'s LevelRow).  Its
#     own headers are "azure2/convert/..." and, in src/, the plain names of src/'s headers.
for dirpath, _, files in os.walk(os.path.join(ROOT, "convert")):
    own = {n for n in os.listdir(os.path.join(ROOT, "convert", "src"))} if os.path.isdir(os.path.join(ROOT, "convert", "src")) else set()
    for name in sorted(files):
        if not name.endswith((".h", ".hpp", ".cpp")):
            continue
        rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
        for inc in re.findall(r'^\s*#\s*include\s*([<"][^>"]+[>"])', read(rel), re.M):
            if inc.startswith('"') and inc[1:-1].startswith(("azure2/convert/", "azure2/io/", "azure2/experiment/", "azure2/physics/")):
                continue
            if inc.startswith('"') and inc[1:-1] in own and rel.startswith(os.path.join("convert", "src") + os.sep):
                continue
            if inc.startswith("<") and re.fullmatch(r"[a-z_]+", inc[1:-1]):
                continue
            problems.append(f"{rel}: includes {inc}; convert/ may include only azure2/convert/, azure2/io/, azure2/experiment/, "
                            "azure2/physics/ headers, its own src/ headers and the C++ standard library")


# 10. One set of physical constants (docs/upgrade/project.md section 5, step 5): physics/'s Constants.h holds u, hbar c,
#     alpha and m_e; any other copy in the engine, the API, io/ or pyazr must carry the same value (comments aside),
#     so that a constant changed in one place cannot leave a stale copy that only some results use.  KNOWN lists the
#     copies that differ on purpose, each with its reason; each is a decision for the maintainers, not drift.
consts_src = read("physics/include/azure2/physics/Constants.h")
CONSTANTS = {}
for name, prefix in (("kHbarC_MeV_fm", "197.32"), ("kAtomicMassUnit_MeV", "931.49"), ("kElectronMass_MeV", "0.51099"),
                     ("kFineStructure", "137.03")):
    m = re.search(name + r"\s*=\s*(?:1\.0*\s*/\s*)?([0-9.]+)", consts_src)
    if not m:
        problems.append(f"physics Constants.h: {name} not found; update this test")
        continue
    CONSTANTS[prefix] = (name, float(m.group(1)))
KNOWN = {
    # m_e of CODATA 2002 in the beta-decay Fermi-function argument, named and kept on purpose (maintainers' decision
    # 2026-11-10: no result moves; one line to change); CODATA 2018 moves beta-delayed outputs (~1e-7).
    ("include/Constants.h", "0.510998903"),
}
for top in ("src", "include", "api", "io", "pyazr"):
    for dirpath, _, files in os.walk(os.path.join(ROOT, top)):
        for name in sorted(files):
            if not name.endswith((".h", ".cpp", ".py")):
                continue
            rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
            in_block = False
            for line in read(rel).splitlines():
                code = line
                if name.endswith(".py"):
                    code = code.split("#", 1)[0]
                else:
                    if in_block:
                        if "*/" not in code:
                            continue
                        code, in_block = code.split("*/", 1)[1], False
                    code = re.sub(r"/\*.*?\*/", "", code)
                    if "/*" in code:
                        code, in_block = code.split("/*", 1)[0], True
                    code = code.split("//", 1)[0]
                for literal in re.findall(r"(?<![0-9.])[0-9]+\.[0-9]+", code):
                    for prefix, (cname, value) in CONSTANTS.items():
                        if literal.startswith(prefix) and float(literal) != value and (rel, literal) not in KNOWN:
                            problems.append(f"{rel}: {literal} is a copy of {cname} = {value!r} (physics Constants.h) "
                                            "with another value; use the constant")
for rel, literal in sorted(KNOWN):
    if literal not in read(rel):
        problems.append(f"{rel}: {literal} is gone; remove it from KNOWN in check 10")



# 10b. One table of element symbols: physics/'s Elements.h (ElementSymbol, ElementZ; the binding's element_symbols(),
#      which pyazr.system.nuclei.SYMBOLS is).  A second list in the engine, the API, io/, experiment/ or pyazr -- three
#      consecutive symbols written out -- is a copy that drifts (src/ParameterLabel.cpp's stopped at Kr).
ELEMENT_RUN = re.compile(r'\bHe\b["\',\s]{1,6}\bLi\b["\',\s]{1,6}\bBe\b')
for top in ("src", "include", "api", "io", "experiment", "pyazr"):
    for dirpath, _, files in os.walk(os.path.join(ROOT, top)):
        for fname in sorted(files):
            rel = os.path.relpath(os.path.join(dirpath, fname), ROOT)
            if fname.endswith((".h", ".cpp", ".py")) and "/tests/" not in "/" + rel and ELEMENT_RUN.search(read(rel)):
                problems.append(f"{rel}: a list of element symbols; use physics/'s Elements.h "
                                "(pyazr: _azure2.element_symbols())")
if not ELEMENT_RUN.search(read("physics/src/Elements.cpp")):
    problems.append("check 10b: physics/src/Elements.cpp does not match the element pattern; update this check")

# 11. One objective (docs/upgrade/objective.md): the chi^2 terms beyond the data points -- Park's wall, abundances,
#     the model data (scattering lengths, <modelData>), target-effect priors, the data covariance -- are added in src/Objective.cpp, which every
#     path (CLI fit, gradient, API, writer) asks.  A new path that adds one itself drifts from the rest the first time
#     a term changes.  Allowed elsewhere: the files that define the terms, and the listed uses with their reason.
PENALTIES = ("ScatteringLength::Chi2",
             "ScatteringLength::Residuals", "ScatteringLength::Chi2Gradient", "ModelData::Chi2", "ModelData::Residuals",
             "ModelData::Gradient", "ModelData::AppendRows",
             "AppendTargetPriors", "CovarianceChi2Correction")
DEFINERS = {"src/Objective.cpp", "src/ScatteringLength.cpp", "src/ModelData.cpp", "src/WorkingCopy.cpp",
            "src/Segment.cpp", "src/TargetEffect.cpp", "src/DiagonalPriors.cpp"}
ALLOWED_PENALTY_USES = {
    # The API's analytic gradient: T2 for its seeds and rows; pyazr_objective checks it against the CLI's objective.
    # (The API's data chi^2, EvaluateFilledChi2, asks Objective::DataChi2 since 2026-11-16.)
    ("api/src/AZUREAPIGradient.cpp", "CovarianceChi2Correction"),
    # The output writer (docs/upgrade/objective_arrays.md section 4): chiSquared.out prints the target-effect and
    # apparatus prior totals from the term's own gather at the working copy's present values (a stored total could be
    # stale after pyazr's write_outputs).
    ("src/DataOutputFiles.cpp", "AppendTargetPriors"),
    ("src/ApparatusParameters.cpp", "AppendTargetPriors"),
    # The data covariance, moved out of src/WorkingCopy.cpp (docs/upgrade/edata_split.md W2): the file defines the
    # covariance term, and PrepareCovarianceSeeds evaluates it for the seeds the gradient reads.
    ("src/DataCovarianceGroups.cpp", "CovarianceChi2Correction"),
}
for top in ("src", "api/src"):
    for name in sorted(os.listdir(os.path.join(ROOT, top))):
        rel = f"{top}/{name}"
        if not name.endswith(".cpp") or rel in DEFINERS:
            continue
        text = read(rel)
        for fn in PENALTIES:
            for m in re.finditer(r"(?<![\w:])" + re.escape(fn) + r"\(", text):
                line = text[text.rfind("\n", 0, m.start()) + 1:text.find("\n", m.start())]
                if line.lstrip().startswith(("//", "*")) or (rel, fn) in ALLOWED_PENALTY_USES:
                    continue
                problems.append(f"{rel}: calls {fn}() -- a chi^2 term outside src/Objective.cpp; ask the objective "
                                "(Objective::Engine().Chi2 / ModelChi2 / DataChi2) or add a term there "
                                "(docs/upgrade/objective.md)")


# 12. One reading of the project file (docs/upgrade/project.md): the engine opens the .azr only in src/Config.cpp
#     (InputFiles + io::Document); every reader takes its section from Config::SectionStream.  Opening the file
#     anywhere else is a second parser.
ALLOWED_OPENS = {"src/ParameterLimitsManager.cpp"}  # only whether the file exists (false = no file, true = no block)
for top in ("src", "api/src"):
    for name in sorted(os.listdir(os.path.join(ROOT, top))):
        rel = f"{top}/{name}"
        if not name.endswith(".cpp") or rel == "src/Config.cpp" or rel in ALLOWED_OPENS:
            continue
        for line in read(rel).splitlines():
            if re.search(r"(ifstream\s*\w*\s*[({]|\.open\s*\(|fopen\s*\(|InputFiles::Open\s*\()[^;]*configfile", line):
                problems.append(f"{rel}: opens the project file ({line.strip()}); take the section from "
                                "Config::SectionStream")

# 18. The format-3 reference page is the table's, not prose that can drift: docs/source/reference/format3.rst equals
#     what tools/gen_schema_docs.py makes of tests/schema/schema_v3.json (docs/upgrade/project_format3.md section 7.4),
#     which in turn equals AZURE2 --schema (check 21: the ctests experiment_schema and pyazr_schema).
try:
    import importlib.util
    _spec = importlib.util.spec_from_file_location("gen_schema_docs", os.path.join(ROOT, "tools", "gen_schema_docs.py"))
    _gen = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_gen)
    import json
    if _gen.generate(json.loads(read("tests/schema/schema_v3.json"))) != read("docs/source/reference/format3.rst"):
        problems.append("docs/source/reference/format3.rst: not the output of tools/gen_schema_docs.py; run it "
                        "(after AZURE2 --schema-out=tests/schema/schema_v3.json if the table changed)")
except (OSError, ValueError) as e:
    problems.append(f"check 18: {e}")


# 27. Every model-data quantity is evaluated and tested (CLAUDE.md recipe "A model-only datum"; data_model.md M10): each
#     quantity of the table for record kind datum (experiment/src/Quantities.cpp, q.recordKinds = kDatum) has its row
#     in ModelQuantities() (experiment/src/ModelDatum.cpp: the engine's evaluation, src/ModelData.cpp, switches on it)
#     and a test named after it, a shell function `test_<quantity>()` or a Python `def test_<quantity>(` under tests/.
quantities_text = read("experiment/src/Quantities.cpp")
datum_quantities = []
for body in re.split(r"\nQuantity \w+\(\) \{", quantities_text)[1:]:
    body = body.split("\n}\n", 1)[0]
    m = re.search(r'q\.name = "(\w+)";', body)
    if m and re.search(r"q\.recordKinds = kDatum;", body):
        datum_quantities.append(m.group(1))
if not datum_quantities:
    problems.append("check 27: no datum quantities found in experiment/src/Quantities.cpp; update this test")
evaluated = set(re.findall(r'(?m)^\s*\{"(\w+)", ModelQuantity::\w+\}', read("experiment/src/ModelDatum.cpp")))
test_names = set()
for dirpath, _, files in os.walk(os.path.join(ROOT, "tests")):
    for name in files:
        if name.endswith((".sh", ".py")):
            text = open(os.path.join(dirpath, name), encoding="utf-8", errors="replace").read()
            test_names.update(re.findall(r"(?m)^\s*(?:def )?test_(\w+)\s*\(", text))
for q in datum_quantities:
    if q not in evaluated:
        problems.append(f"{q}: a datum quantity of the table without a row in ModelQuantities() "
                        "(experiment/src/ModelDatum.cpp) and its evaluation in src/ModelData.cpp")
    if q not in test_names:
        problems.append(f"{q}: a datum quantity without a test named after it (test_{q} under tests/)")


# 24. Every diagnostic code is declared, listed in the schema and tested (docs/upgrade/project_format3.md section 7.4):
#     the codes are the constants of experiment/'s Diagnostics.h (namespace codes) and io/'s Record.h (record_codes);
#     each is in DiagnosticCodes() (experiment/src/Diagnostics.cpp, which the schema exports); each appears quoted in
#     a test of experiment/tests/ or io/tests/; and every "E_..."/"W_..." literal in experiment/src/ and io/src/ is a
#     declared code (a code made up at the call site would bypass the list).
declared = {}
for rel in ("experiment/include/azure2/experiment/Diagnostics.h", "io/include/azure2/io/Record.h"):
    for const, code in re.findall(r'constexpr const char \*(k\w+) = "([EWA][A-Z_]*)"', read(rel)):
        declared[const] = code
registry = read("experiment/src/Diagnostics.cpp")
tests_text = ""
for top in ("experiment/tests", "io/tests"):
    for name in sorted(os.listdir(os.path.join(ROOT, top))):
        if name.endswith(".cpp"):
            tests_text += read(f"{top}/{name}")
if not declared:
    problems.append("check 24: no diagnostic codes found; update this test")
for const, code in sorted(declared.items()):
    if not re.search(r"\{\s*(?:rc::)?" + const + r"\s*,", registry):
        problems.append(f"{code}: declared but not in DiagnosticCodes() (experiment/src/Diagnostics.cpp)")
    if f'"{code}"' not in tests_text:
        problems.append(f"{code}: no test in experiment/tests/ or io/tests/ names it (each code needs one that triggers it)")
for top in ("experiment/src", "io/src"):
    for name in sorted(os.listdir(os.path.join(ROOT, top))):
        if not name.endswith(".cpp"):
            continue
        for code in sorted(set(re.findall(r'"([EW]_[A-Z][A-Z_]+)"', read(f"{top}/{name}")))):
            if code not in declared.values():
                problems.append(f"{top}/{name}: diagnostic code {code} is not declared in Diagnostics.h or Record.h")


# 23. No code recognises a parameter's kind by its name (docs/upgrade/parameter_kinds.md; project_format3.md section
#     7.4): what a parameter is comes from its record (inference/include/azure2/inference/ParameterKind.h), written
#     where it is created; a name
#     is a label, and where a file names a parameter the name is looked up whole (AZUREParams::IndexOfName).  Caught:
#     in src/, include/, api/ a find / rfind / compare / starts_with / ends_with / EndsWith / StartsWith whose string
#     is a stem of a parameter name ("segment_", "_norm", "width", ...), and the deleted helpers by name; in pyazr/ a
#     startswith / endswith on such a stem or on a parameter's .name against a computed prefix, and a regular
#     expression containing a stem.  Not caught: building a name (the creation sites), command-line flags and other
#     strings.  ALLOWED_NAME_TESTS lists the file grammars read without a parameter vector, each with its reason.
STEM = (r'(?:energy_|width_|segment_|abundance_|screening_|density_|doppler_|cbkg_|_norm|_energy_shift|_thickness|'
        r'_sqrt|_fn|segment"|width"|Width"|norm"|shift")')
CPP_NAME_TEST = re.compile(r'(?:\.\s*(?:r?find|compare|starts_with|ends_with)\s*\(\s*(?:0\s*,\s*\d+\s*,\s*)?'
                           r'|\b(?:EndsWith|StartsWith)\s*\([^,()]*,\s*)"' + STEM)
CPP_HELPERS = re.compile(r'\b(?:IsNormName|IsEnergyShiftName|IsEnergyShiftFnName|IsThicknessName)\b')
PY_NAME_TEST = re.compile(r'\.(?:startswith|endswith)\(\s*[rbf]*["\']' + STEM.replace('"', '["\']')
                          + r'|\.name\.(?:startswith|endswith)\(\s*[A-Za-z_]')
PY_NAME_REGEX = re.compile(r're\.(?:match|fullmatch|search|compile)\(\s*[rbf]*["\'][^"\']*'
                           r'(?:energy_|width_|segment_|abundance_|screening_|density_|doppler_|cbkg_)')
ALLOWED_NAME_TESTS = set()  # (file, text of the line): none today
name_test_hits = set()
for top in ("src", "include", "api", "pyazr", "inference"):
    for dirpath, _, files in os.walk(os.path.join(ROOT, top)):
        for fname in sorted(files):
            if not fname.endswith((".h", ".cpp", ".py")):
                continue
            rel = os.path.relpath(os.path.join(dirpath, fname), ROOT)
            in_block = False
            for number, line in enumerate(read(rel).splitlines(), 1):
                code = line
                if fname.endswith(".py"):
                    code = re.sub(r'#(?=[^"\']*$).*', "", code)
                    hit = PY_NAME_TEST.search(code) or PY_NAME_REGEX.search(code)
                else:
                    if in_block:
                        if "*/" not in code:
                            continue
                        code, in_block = code.split("*/", 1)[1], False
                    code = re.sub(r"/\*.*?\*/", "", code)
                    if "/*" in code:
                        code, in_block = code.split("/*", 1)[0], True
                    code = code.split("//", 1)[0]
                    hit = CPP_NAME_TEST.search(code) or CPP_HELPERS.search(code)
                if not hit:
                    continue
                allowed = [a for a in ALLOWED_NAME_TESTS if a[0] == rel and a[1] in line]
                if allowed:
                    name_test_hits.update(allowed)
                    continue
                problems.append(f"{rel}:{number}: tests a parameter's name ({line.strip()}); ask its record "
                                "(AZUREParams::Record / Kind, pyazr Parameter.kind and owner fields; "
                                "docs/upgrade/parameter_kinds.md)")
for rel, text in sorted(ALLOWED_NAME_TESTS - name_test_hits):
    problems.append(f"{rel}: '{text}' is gone; remove it from ALLOWED_NAME_TESTS in check 23")
# The kinds: one enum, one table row per kind in the order of the codes, each documented for maintainers and users.
KINDS_H = "inference/include/azure2/inference/ParameterKind.h"
kinds_h = read(KINDS_H)
enum_body = kinds_h[kinds_h.index("enum class ParamKind"):]
enum_body = enum_body[:enum_body.index("};")]
enum_kinds = re.findall(r"^\s*(\w+) = (\d+),?", enum_body, re.M)
row_kinds = re.findall(r'\{ParamKind::(\w+), "(\w+)"', kinds_h)
if [k for k, _ in enum_kinds] != [k for k, _ in row_kinds] or \
        [int(c) for _, c in enum_kinds] != list(range(len(enum_kinds))):
    problems.append(KINDS_H + ": ParamKindTable() must have one row per ParamKind, in the order of the "
                    "codes 0, 1, 2, ...")
for rel in ("docs/upgrade/parameter_kinds.md", "docs/source/reference/output_files.rst"):
    text = read(rel)
    for kind, key in row_kinds:
        if f"``{key}``" not in text and f"`{key}`" not in text:
            problems.append(f"{rel}: the parameter kind {kind} ('{key}') is not documented")
for other in ("include/AZUREGrad.h", "include/AZURECalcMCMC.h", "pyazr/parameters.py"):
    if re.search(r"enum (?:class )?ParamKind\b|PARAM_ENERGY|^_KINDS\s*=", read(other), re.M):
        problems.append(f"{other}: a second enumeration of parameter kinds; use ParameterKind.h "
                        "(pyazr: parameter_kinds())")

# 25. The results of an evaluation belong to its caller (docs/upgrade/data_model.md S3, include/Results.h): a point and a
#     data set hold where a record lies, never the values, so a copy of either cannot write into another owner's
#     buffer.  Site.h, Unit.h (with UnitTables.h, its tables) and WorkingCopy.h declare no member of a type that names
#     Results -- a value, pointer, reference, container or smart pointer.  Function declarations may name it (the
#     buffer is passed in).
for rel in ("include/Site.h", "include/Unit.h", "include/UnitTables.h", "include/WorkingCopy.h"):
    text = re.sub(r"//[^\n]*|/\*.*?\*/|^\s*#[^\n]*", "", read(rel), flags=re.S | re.M)
    for stmt in re.split(r"[;{}]", text):
        if "(" in stmt or not re.search(r"\bResults\b", stmt):
            continue  # a function (its parameters may name the buffer), or nothing about results
        problems.append(f"{rel}: declares a results member ({' '.join(stmt.split())}); the buffer belongs to the "
                        "caller of the evaluation and is passed in (include/Results.h)")

# 28. One model pass and one residual (docs/upgrade/data_model.md S4, include/ModelPass.h): what a data set predicts
#     is combined in src/ModelPass.cpp only (ComposeModelValue, TransformedValue, reached through
#     azure2::EvaluateModelPoint / EvaluateModelSegment / EvaluateModel), and the residual (f - y n)/(dy n) is written
#     once (azure2::DataResidual, PointResidual).  A second combination loop or residual expression drifts from the
#     first the next time either changes (before S4: 19 loops and 9 residuals).
COMPOSITION = re.compile(r"(?<![\w:])(?:ComposeModelValue|TransformedValue|CalculateTheoreticalCrossSection)\s*\(")
RESIDUAL = re.compile(r"[\w)\]]\s*-(?![>-])\s*[\w\->.()]*?(?:GetCMCrossSection\(\)|\bdataval\b|\bd\[\w+\])\s*\*\s*"
                      r"[\w\[\].>()-]*?\b(?:n|norm|dataNorm|GetNorm)\b")
for top in ("src", "api/src", "include", "tests/upgrade", "tests/reference"):
    for name in sorted(os.listdir(os.path.join(ROOT, top))):
        rel = f"{top}/{name}"
        if not name.endswith((".cpp", ".h")) or rel in ("src/ModelPass.cpp", "include/ModelPass.h"):
            continue
        text = re.sub(r"//[^\n]*|/\*.*?\*/", "", read(rel), flags=re.S)
        if COMPOSITION.search(text):
            problems.append(f"{rel}: combines a segment's model value itself; call azure2::EvaluateModelPoint "
                            "(include/ModelPass.h), the one model pass")
        for m in RESIDUAL.finditer(text):
            problems.append(f"{rel}: writes a data residual ({m.group(0).strip()}); use azure2::DataResidual or "
                            "PointResidual (include/ModelPass.h)")
if not RESIDUAL.search("double r = th - pt->GetCMCrossSection() * norm;"):
    problems.append("check 28: the residual pattern no longer matches the form it guards; update this check")
# ... and the composition is the expression tree's (step S10, include/Predictions.h): the pieces it combines -- a
# term's scaling and frame factor, a sample's abundance weight, the transmission / yield transform, the transform's
# term list -- are read only where the tree's inputs are filled (src/ModelPass.cpp) or its setup is made
# (src/Predictions.cpp), and defined in Segment; experiment/ holds the arithmetic.  The adjoint of the composition is
# the tree's backward too (step S11, azure2::ModelPointBackward), so the gradient reads none of them.
PIECES = re.compile(r"(?<![\w:])(?:Segment::)?(?:TermFrameFactor|RatioTermFrameFactor|TransformTerms|Transmission|YieldY0)"
                    r"\s*\(|(?:->|\.)GetComponentScaling\s*\(|(?:->|\.)GetAbundance\s*\(\s*\w*[Ss]ample")
PIECES_ALLOWED = {"src/ModelPass.cpp", "src/Predictions.cpp", "src/Segment.cpp", "include/Segment.h"}
for top in ("src", "api/src", "include"):
    for name in sorted(os.listdir(os.path.join(ROOT, top))):
        rel = f"{top}/{name}"
        if not name.endswith((".cpp", ".h")) or rel in PIECES_ALLOWED:
            continue
        text = re.sub(r"//[^\n]*|/\*.*?\*/", "", read(rel), flags=re.S)
        for m in PIECES.finditer(text):
            problems.append(f"{rel}: composes a model value from its pieces ({m.group(0).strip()}); a data set's "
                            "prediction is its expression tree (include/Predictions.h, evaluated by "
                            "azure2::EvaluateModelPoint)")
if not PIECES.search("w = compound->GetAbundance(sample, nucleus) * Segment::TermFrameFactor(p, q);"):
    problems.append("check 28: the composition pattern no longer matches the form it guards; update this check")


# 29. The calculation plan (docs/upgrade/data_model.md S5, include/CalculationPlan.h): which point calculates which is
#     the plan's, held by position in the data set (PointRef); a point holds the index of its task, never a pointer
#     to another point.  A pointer from one point into another made Clone re-map every list, and a missed one had a
#     working copy write into the original's points.  So: (a) no data member of Unit is a pointer to an Unit or a
#     container of such pointers; (b) the point mapping (a point calculated by another point at its energy: the old
#     mapped-point lists, MapData, MapLeaves, SharesCalculation, the AZURE_NO_*_MAPPING switches) is gone -- shared
#     energies are the energy tasks'; (c) a plan position is turned into a point (WorkingCopy::PlanPoint) only by the
#     plan (src/CalculationPlan.cpp) and the plan's tests -- everything else asks the plan.
POINT_PTR = re.compile(r"^\s*(?:mutable\s+)?(?:std::\w+<\s*)?(?:const\s+)?Unit\s*\*(?:\s*>)?\s*\w+_\s*(?:=[^;]*)?;", re.M)
eph = re.sub(r"//[^\n]*|/\*.*?\*/", "", read("include/Site.h") + read("include/Unit.h"), flags=re.S)
# the RAII scopes nested in Unit (DopplerAppliedScope, ...) hold the point they set a flag on: not a member of Unit
eph = re.sub(r"\bclass\s+\w+Scope\b[^{;]*\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\};", "", eph)
if "plan_task_" not in eph:
    problems.append("check 29: Unit's members are no longer found in include/Unit.h; update this check")
for m in POINT_PTR.finditer(eph):
    problems.append(f"include/Unit.h: a point holds a pointer to a point ({m.group(0).strip()}); the calculation "
                    "plan (include/CalculationPlan.h) says who calculates whom")
if not POINT_PTR.search("  std::vector<Unit *> mapped_;") or not POINT_PTR.search("  Unit *host_ = nullptr;"):
    problems.append("check 29: the member pattern no longer matches the forms it guards; update this check")
OLD_MAPPING = re.compile(r"\b(?:\w*LocalMappedPoints?|local_mapped_points_|GetMappedPoints|EPointMapIterator|PlanMapped|"
                         r"MappedPoints|MapData|MapLeaves|SharesCalculation|IsMapped|AZURE_NO_\w*MAPPING)\b")
PLAN_POINT = re.compile(r"\bPlanPoint\s*\(")
PLAN_POINT_FILES = {"src/CalculationPlan.cpp", "include/WorkingCopy.h", "tests/upgrade/plan_test.cpp"}
for top in ("src", "api/src", "include", "tests/upgrade", "tests/reference"):
    for name in sorted(os.listdir(os.path.join(ROOT, top))):
        rel = f"{top}/{name}"
        if not name.endswith((".cpp", ".h")):
            continue
        text = re.sub(r"//[^\n]*|/\*.*?\*/", "", read(rel), flags=re.S)
        if OLD_MAPPING.search(text):
            problems.append(f"{rel}: uses the deleted point mapping ({OLD_MAPPING.search(text).group(0)}); a point is "
                            "calculated by its own plan task, a shared energy by the energy tasks (CalculationPlan.h)")
        if rel not in PLAN_POINT_FILES and PLAN_POINT.search(text):
            problems.append(f"{rel}: resolves a plan position itself (WorkingCopy::PlanPoint); ask the calculation plan "
                            "(src/CalculationPlan.cpp)")


# 30. No path of one machine in the repository: a tracked file that names a directory under /data0 or /home works on
#     the machine it was written on only, and tells a reader of a clone to run something they do not have.  Tools are
#     found beside each other (tools/), scratch is an environment variable ($WORK), data files are relative to the
#     project.  docs/upgrade/HANDOFF.md is a log of what was run where and is left out; ALLOWED_MACHINE_PATHS lists
#     the lines that name such a path on purpose (a rule against it, or its test), each with its reason.  Binary
#     files are skipped.
MACHINE_PATH = re.compile(r"/(?:data0|home)/")
MACHINE_PATH_FREE_FILES = {"docs/upgrade/HANDOFF.md"}
HOME = "/" + "home/"  # spelled in two parts so that this file passes its own check
ALLOWED_MACHINE_PATHS = {
    # The skills checker's rule R8 refuses personal paths in skills: its description and its pattern name them.
    ("tools/check_skills.py", "R8  no personal absolute paths"),
    ("tools/check_skills.py", r'r"(?:(?<![\w.$})])(?:/groups/|/Users/|/users/|' + HOME),
    # ... and its test feeds such a path to the rule.
    ("tests/consistency/check_skills_test.py", f'"See {HOME}me/x."'),
}


def tracked_files():
    """The files git tracks (a fresh clone's), or, without git, every file of the tree outside build and output
    folders."""
    try:
        out = subprocess.run(["git", "-c", "safe.directory=*", "ls-files", "-z"], cwd=ROOT, capture_output=True,
                             check=True).stdout
        return sorted(f for f in out.decode("utf-8", "replace").split("\0") if f)
    except (OSError, subprocess.CalledProcessError):
        found = []
        for dirpath, dirnames, files in os.walk(ROOT):
            dirnames[:] = [d for d in dirnames if d != ".git" and not d.startswith("build")
                           and d not in ("output", "checks", "__pycache__", "_build")]
            found += [os.path.relpath(os.path.join(dirpath, f), ROOT) for f in files]
        return sorted(found)


machine_path_hits = set()
for rel in tracked_files():
    if rel in MACHINE_PATH_FREE_FILES:
        continue
    try:
        with open(os.path.join(ROOT, rel), "rb") as f:
            raw = f.read()
    except OSError:
        continue  # a tracked file deleted in the work tree
    if b"\0" in raw:
        continue
    for number, line in enumerate(raw.decode("utf-8", "replace").splitlines(), 1):
        if not MACHINE_PATH.search(line):
            continue
        allowed = [a for a in ALLOWED_MACHINE_PATHS if a[0] == rel and a[1] in line]
        if allowed:
            machine_path_hits.update(allowed)
            continue
        problems.append(f"{rel}:{number}: names a path of one machine ({MACHINE_PATH.search(line).group(0)}...); "
                        "use a path relative to the repository or the project, or an environment variable")
for rel, text in sorted(ALLOWED_MACHINE_PATHS - machine_path_hits):
    problems.append(f"{rel}: '{text}' is gone; remove it from ALLOWED_MACHINE_PATHS in check 30")
if not MACHINE_PATH.search("tools/runall.sh is in /" + "data0/x/tools") or not MACHINE_PATH.search("see " + HOME):
    problems.append("check 30: the path pattern no longer matches the forms it guards; update this check")

# 31. The data classes are named after what they hold (docs/upgrade/data_structure.md): the data set (DataSet), the
#     working copy (WorkingCopy), its units (Unit), a data set's segment (Segment).  The old names EData, ESegment and
#     EPoint, which mixed settings and state, do not come back: no tracked file of the engine is named after them and no
#     code names them.
OLD_DATA_NAME = re.compile(r"\b(?:EData|ESegment|EPoint)\b")
#     Two texts output files carry, kept as they were written so that those files stay byte-identical: the header of
#     intEC.dat.sig (changes with the next signature version) and a comment the converter writes into a format-3
#     project (changes with the converter's next fixture regeneration).
OLD_DATA_NAME_OUTPUT = {("src/WorkingCopy.cpp", "see EData::ECSignature)"),
                        ("convert/src/BuildData.cpp", "(the same index, EData::MakePoints)")}
for rel in tracked_files():
    if not rel.startswith(("src/", "include/", "api/", "tests/", "experiment/", "convert/", "io/", "physics/", "inference/")):
        continue
    if re.match(r"(?:EData|ESegment|EPoint)", os.path.basename(rel)):
        problems.append(f"{rel}: a file named after an old data class; see check 31")
    if not rel.endswith((".cpp", ".h", ".hpp")):
        continue
    try:
        text = read(rel)
    except (OSError, UnicodeDecodeError):
        continue
    for number, line in enumerate(text.splitlines(), 1):
        if OLD_DATA_NAME.search(line) and not any(r == rel and t in line for r, t in OLD_DATA_NAME_OUTPUT):
            problems.append(f"{rel}:{number}: names an old data class ({OLD_DATA_NAME.search(line).group(0)}); the data "
                            "set is DataSet, the working copy WorkingCopy, a unit Unit, a segment Segment (check 31)")
if not OLD_DATA_NAME.search("EData *ld = data();") or OLD_DATA_NAME.search("EDataLike ECAmplitudeCache"):
    problems.append("check 31: the old-name pattern no longer matches the forms it guards; update this check")

# 33. A check's text edits can be shown to bite (docs/upgrade/format3_tree.md section 7): an edit in place (sed -i,
#     perl -i, awk -i inplace) leaves no copy to compare with, so a pattern that matches nothing goes unnoticed; a
#     tests/*/check.sh edits through tests/lib/edit_guard.sh (`edited SRC DEST awk ...`) instead.  The checks listed
#     in IN_PLACE_AUDITED had in-place edits when the rule came in, and tools/edit_audit.sh showed that each of them
#     changes its file; a new in-place edit elsewhere is listed.
IN_PLACE = re.compile(r"\bsed\s+(?:-[a-zA-Z]+\s+)*-i|\bperl\s+-[a-zA-Z0-9]*i|\bawk\s+-i\s+inplace")
IN_PLACE_AUDITED = {"data_covariance", "format3_engine", "extrapolation_grids", "validate_run_refusals"}
for name in sorted(os.listdir(os.path.join(ROOT, "tests"))):
    rel = f"tests/{name}/check.sh"
    if name in IN_PLACE_AUDITED or not os.path.isfile(os.path.join(ROOT, rel)):
        continue
    for number, line in enumerate(read(rel).splitlines(), 1):
        if not line.lstrip().startswith("#") and IN_PLACE.search(line):
            problems.append(f"{rel}:{number}: edits a file in place; edit through tests/lib/edit_guard.sh (`edited SRC "
                            "DEST awk ...`), which fails when the edit changes nothing (check 33)")
if not IN_PLACE.search("sed -i.bak 's/a/b/' t.azr") or not IN_PLACE.search("perl -pi -e 's/a/b/' t.azr") or \
        IN_PLACE.search("sed 's/a/b/' t.azr > u.azr"):
    problems.append("check 33: the in-place pattern no longer matches the forms it guards; update this check")

# 34. The engine reads format 3 only (docs/upgrade/project_format3.md section 8.4): what knows the old sections and the
#     old formats is the converter (convert/), and the one refusal of an old file names them (io/src/
#     ProjectFormat.cpp).  A file of the engine (src/, include/, api/, io/, experiment/) that names <segmentsData>,
#     <targetInt>, format 2 or "old format" is listed: a new reader of an old section would be a second grammar.
#     OLD_WORDS_PENDING: files of the two branches that work in parallel (node layout, fit-parameter plumbing) whose
#     comments still name them; they leave this list when those branches are merged and the words are gone.
OLD_WORDS = re.compile(r"segmentsData|targetInt\b|format[ -]2\b|old format", re.I)
OLD_WORDS_ALLOWED = {"io/src/ProjectFormat.cpp"}
OLD_WORDS_PENDING = {"include/Unit.h", "inference/include/azure2/inference/ParameterKind.h", "src/TargetSetup.cpp",
                     "src/DataParameters.cpp",
                     "experiment/src/Schema.cpp"}
for top in ("src", "include", "api", "io", "experiment", "inference"):
    for dirpath, _, files in os.walk(os.path.join(ROOT, top)):
        if "/tests" in dirpath[len(ROOT):]:
            continue
        for f in files:
            if not f.endswith((".cpp", ".h")):
                continue
            rel = os.path.relpath(os.path.join(dirpath, f), ROOT)
            if rel in OLD_WORDS_ALLOWED or rel in OLD_WORDS_PENDING:
                continue
            for number, line in enumerate(read(rel).splitlines(), 1):
                if OLD_WORDS.search(line):
                    problems.append(f"{rel}:{number}: names an old section or format ({OLD_WORDS.search(line).group(0)}); "
                                    "the engine reads format 3 only, the old grammar is the converter's (check 34)")
if not OLD_WORDS.search("<segmentsData> line 3") or OLD_WORDS.search("bool targetIntegration = false;"):
    problems.append("check 34: the old-word pattern no longer matches the forms it guards; update this check")

# 35. A leaf the R-matrix does not compute (a Coulomb reference, a tabulated reference) is asked for with
#     IsReferenceLeaf (experiment/Prediction.h) wherever the model pass skips or drops it (ComposeModelValue,
#     ModelPointBackward): a test with the narrower IsCoulombReference would miss a new reference leaf.  In src/ it is
#     named only where the Coulomb value itself is chosen: once in src/ModelPass.cpp, and in src/CoulombReference.cpp.
NARROW_LEAF = re.compile(r"\bIsCoulombReference\(")
for rel in sorted(f"src/{n}" for n in os.listdir(os.path.join(ROOT, "src")) if n.endswith(".cpp")):
    if rel == "src/CoulombReference.cpp":
        continue
    hits = [n for n, line in enumerate(read(rel).splitlines(), 1) if NARROW_LEAF.search(line)]
    if len(hits) > (1 if rel == "src/ModelPass.cpp" else 0):
        problems.append(f"{rel}:{hits[-1]}: IsCoulombReference where a reference leaf is skipped or dropped; ask "
                        "IsReferenceLeaf (check 35)")
if not NARROW_LEAF.search("if (azure2::experiment::IsCoulombReference(part))") or NARROW_LEAF.search("IsReferenceLeaf(p)"):
    problems.append("check 35: the pattern no longer matches the form it guards; update this check")

# 37. The sharing key of the mesh compares the term settings field by field (SameObservable, include/TermSettings.h):
#     a field added to TermSettings and forgotten there makes two different observables share their sites silently.
#     Every field is either compared in SameObservable (a.f == b.f, or by value: SameX(a.f, b.f)) or named in
#     kReaderOnlyTermSettings beside it,
#     never both, and both lists name fields that exist.
def term_settings_problems(text):
    """The problems of check 37 in the text of TermSettings.h (a function, so that the check can be shown to bite)."""
    found = []
    text = re.sub(r"//[^\n]*|/\*.*?\*/", "", text, flags=re.S)
    struct = re.search(r"\bstruct\s+TermSettings\s*\{(.*?)\n\};", text, re.S)
    compare = re.search(r"\bbool\s+SameObservable\s*\([^)]*\)\s*\{(.*?)\n\}", text, re.S)
    listed = re.search(r"\bkReaderOnlyTermSettings\[\]\s*=\s*\{([^}]*)\}", text)
    if not (struct and compare and listed):
        return ["check 37: TermSettings, SameObservable or kReaderOnlyTermSettings is no longer found in "
                "include/TermSettings.h; update this check"]
    fields = []
    for statement in struct.group(1).split(";"):
        statement = re.sub(r"=[^,]*", "", statement.strip())  # initialisers (none holds a comma)
        for declarator in statement.split(","):
            name = re.search(r"(\w+)\s*$", declarator.strip())
            if name:
                fields.append(name.group(1))
    # a.f == b.f, or a comparison by value of an attachment, SameX(a.f, b.f)
    compared = set(re.findall(r"\ba\.(\w+)\s*==\s*b\.\1\b", compare.group(1))) | \
        set(re.findall(r"\b\w+\(\s*a\.(\w+)\s*,\s*b\.\1\s*\)", compare.group(1)))
    reader_only = set(re.findall(r'"(\w+)"', listed.group(1)))
    for f in fields:
        if f not in compared and f not in reader_only:
            found.append(f"include/TermSettings.h: the field {f} of TermSettings is neither compared in SameObservable "
                         "nor named in kReaderOnlyTermSettings (check 37): does the energy pass read it?")
        if f in compared and f in reader_only:
            found.append(f"include/TermSettings.h: the field {f} is both compared in SameObservable and named "
                         "reader-only (check 37)")
    for f in sorted((compared | reader_only) - set(fields)):
        found.append(f"include/TermSettings.h: {f} is compared or named reader-only but is no field of TermSettings "
                     "(check 37)")
    return found


problems += term_settings_problems(read("include/TermSettings.h"))
_probe = read("include/TermSettings.h").replace("int maxAngDistOrder = 0;", "int maxAngDistOrder = 0;\n  int probe = 0;")
if not any(" probe " in p for p in term_settings_problems(_probe)):
    problems.append("check 37: a field added to TermSettings is not caught; update this check")

# 38. No hidden state (docs/upgrade/evaluation.md; ARCHITECTURE.md section 9): no thread_local, no non-const static or
#     namespace-scope variable, no mutable member in the engine's sources (tests/consistency/state_scan.py) outside the
#     allow-list tests/consistency/allowed_state.tsv (file, symbol, kind, reason); every allowed row still exists.  New
#     state is a table built at set-up or scratch owned by the caller (CLAUDE.md, recipe "Add a cache").
sys.path.insert(0, os.path.join(ROOT, "tests", "consistency"))
import state_scan  # noqa: E402

_allowed = state_scan.read_allowed(os.path.join(ROOT, "tests", "consistency", "allowed_state.tsv"))
_found = state_scan.scan(ROOT)
for (_f, _sym, _kind), _lines in sorted(_found.items()):
    if (_f, _sym, _kind) not in _allowed:
        problems.append(f"{_f}:{_lines[0]}: {_kind} {_sym} is state outside the rule (check 38): make it a table built "
                        "at set-up or the caller's scratch, or a const (docs/upgrade/evaluation.md)")
for (_f, _sym, _kind), _reason in sorted(_allowed.items()):
    if _f == "?" or not _reason.strip():
        problems.append(f"tests/consistency/allowed_state.tsv: a row without file, symbol, kind and reason: {_sym} (check 38)")
    elif (_f, _sym, _kind) not in _found:
        problems.append(f"tests/consistency/allowed_state.tsv: {_f} {_sym} ({_kind}) no longer exists; remove the row "
                        "(check 38)")
_planted = "namespace {\nint g_planted = 0;\n}\nvoid F() {\n  static double memo = 0.0;\n  thread_local int t;\n}\n" \
           "struct S {\n  mutable int m_;\n  static const char *const kOk;\n  std::once_flag once;\n};\n"
if sorted(k for _, k, _ in state_scan.scan_text(_planted)) != ["mutable", "namespace", "once", "static", "thread_local"]:
    problems.append("check 38: a planted namespace variable, static, thread_local, mutable member or once-written "
                    "member is not caught; update tests/consistency/state_scan.py")

# 39. One reader of the environment (docs/upgrade/evaluation.md section 8; support/include/azure2/support/Environment.h):
#     no getenv in the engine, its libraries or the binding outside support/src/Environment.cpp.  Test programs (any
#     tests/ directory) are not the engine.  A switch is a row of that file's table, read through the record.
ENV_READER = os.path.join("support", "src", "Environment.cpp")
GETENV_DIRS = ("src", "include", "api", "physics", "experiment", "stopping", "io", "convert", "support", "inference")


def getenv_problems(rel, text):
    found = []
    code = re.sub(r"//[^\n]*|/\*.*?\*/", "", text, flags=re.S)
    code = re.sub(r'"(?:\\.|[^"\\\n])*"', '""', code)
    if re.search(r"\b(?:std::)?(?:getenv|secure_getenv|_dupenv_s|GetEnvironmentVariable[AW]?)\s*\(", code):
        found.append(f"{rel}: reads the environment (getenv); read the switch from the record "
                     "(azure2/support/Environment.h) or the run's Config::environment (check 39)")
    return found


_getenv_files = 0
for _d in GETENV_DIRS:
    for dirpath, dirnames, files in os.walk(os.path.join(ROOT, _d)):
        dirnames[:] = [x for x in dirnames if x != "tests"]
        for name in sorted(files):
            if not name.endswith((".h", ".hpp", ".cpp", ".cc", ".c")):
                continue
            rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
            if rel == ENV_READER:
                continue
            _getenv_files += 1
            problems += getenv_problems(rel, read(rel))
if not getenv_problems("probe.cpp", 'int f() { return std::getenv("AZURE_X") != nullptr; }'):
    problems.append("check 39: a getenv call is not caught; update this check")
if getenv_problems("probe.cpp", '// getenv("AZURE_X") in a comment\nconst char *t = "getenv(";'):
    problems.append("check 39: a getenv in a comment or a string is reported; update this check")
if "getenv" not in read(ENV_READER):
    problems.append(f"{ENV_READER}: the one reader of the environment no longer calls getenv (check 39)")

# 40. The table of switches and the manual's list are one list (docs/source/reference/command_line.rst, "Environment
#     switches"): every row of support/src/Environment.cpp's table is a row there with its kind, and every row there is
#     in the table.


def switch_table(text):
    return re.findall(r'\{Switch::k\w+, "([A-Z0-9_]+)", T::k\w+, "[^"]*", K::k(\w+),', text)


def manual_switches(text):
    part = text.split("Environment switches\n--------------------\n", 1)
    if len(part) < 2:
        return None
    return [(n, k) for n, _o, k in manual_switch_rows(part[1])]


def manual_switch_rows(part):
    """(variable, option or "none", kind) of each row of the manual's list of switches."""
    return re.findall(r"^   \* - ``([A-Z0-9_]+)``\n     - (?:``(--[a-z0-9-]+(?:=[A-Z]+)?)``|none)\n     - (\w+)$",
                      part.split("\nExamples\n", 1)[0], re.M)


_table = [(n, k.lower()) for n, k in switch_table(read(ENV_READER))]
_manual = manual_switches(read("docs/source/reference/command_line.rst"))
if len(_table) < 30:
    problems.append(f"{ENV_READER}: {len(_table)} switch rows read; the pattern of check 40 no longer matches the table")
elif _manual is None:
    problems.append("docs/source/reference/command_line.rst: no section 'Environment switches' (check 40)")
else:
    for _n, _k in sorted(set(_table) - set(_manual)):
        problems.append(f"docs/source/reference/command_line.rst: the switch {_n} ({_k}) of the table is not listed "
                        "with that kind under 'Environment switches' (check 40)")
    for _n, _k in sorted(set(_manual) - set(_table)):
        problems.append(f"docs/source/reference/command_line.rst: {_n} ({_k}) is listed but is no row of {ENV_READER} "
                        "(check 40)")

# 43. Every switch that can move a number is an option of the program (docs/upgrade/evaluation.md section 9): a row
#     of kind numbers has its option, no row of another kind has one, the options are unique, and the manual's list
#     names each row's option (and none that is no row's).


def switch_options(text):
    return re.findall(r'\{Switch::k\w+, "([A-Z0-9_]+)", T::k\w+, "[^"]*", K::k(\w+),\s*"[^"]*",\s*kManual,\s*'
                      r'"([^"]*)", (-?\d+), (-?\d+)\}', text, re.S)


_options = switch_options(read(ENV_READER))
if len(_options) != len(_table):
    problems.append(f"{ENV_READER}: {len(_options)} rows with an option column read of {len(_table)} (check 43: the "
                    "pattern no longer matches the table)")
for _n, _k, _o, _lo, _hi in _options:
    if _k == "Numbers" and not _o:
        problems.append(f"{ENV_READER}: the switch {_n} can move a number and has no option of the program (check 43)")
    if _k != "Numbers" and _o:
        problems.append(f"{ENV_READER}: the switch {_n} ({_k.lower()}) has an option; only a switch of kind numbers "
                        "is one (check 43)")
_names = [_o.split("=")[0] for _n, _k, _o, _lo, _hi in _options if _o]
for _dup in sorted({o for o in _names if _names.count(o) > 1}):
    problems.append(f"{ENV_READER}: the option {_dup} is the option of two switches (check 43)")
_doc_part = read("docs/source/reference/command_line.rst").split("Environment switches\n--------------------\n", 1)
_doc_options = {n: o for n, o, k in manual_switch_rows(_doc_part[1])} if len(_doc_part) == 2 else {}
for _n, _k, _o, _lo, _hi in _options:
    if _doc_options.get(_n, "") != _o:
        problems.append(f"docs/source/reference/command_line.rst: the switch {_n} is listed with the option "
                        f"'{_doc_options.get(_n, '')}', the table's is '{_o}' (check 43)")

# 41. Every physics function has a test (docs/upgrade/pyazr3.md section 2.7): each row of the binding's physics table
#     (api/src/pybind_physics.cpp, `t.push_back({"<name>", ...`) is named by a test as `row_test("<name>"` in a file
#     under tests/; a row without one, or a row_test of a name the table does not have, fails.
_rows = re.findall(r'^\s*t\.push_back\(\{"(\w+)"', read("api/src/pybind_physics.cpp"), re.M)
_row_tests = set()
for _dir, _, _files in os.walk(os.path.join(ROOT, "tests")):
    for _name in _files:
        if _name.endswith(".py"):
            with open(os.path.join(_dir, _name), errors="replace") as _fh:
                _row_tests |= set(re.findall(r'\brow_test\("(\w+)"', _fh.read()))
if not _rows:
    problems.append("api/src/pybind_physics.cpp: no row of the physics table found (check 41: the pattern changed?)")
for _r in _rows:
    if _r not in _row_tests:
        problems.append(f"api/src/pybind_physics.cpp: the physics function {_r} has no test (check 41): add "
                        f"row_test(\"{_r}\", ...) to tests/pyazr/physics_table_test.py")
for _r in sorted(_row_tests - set(_rows)):
    problems.append(f"tests: row_test(\"{_r}\") names no row of api/src/pybind_physics.cpp (check 41)")

# 42. The pyazr command line is one table (pyazr/cli.py, COMMANDS; docs/upgrade/pyazr3.md section 2.8): every row has
#     its case in tests/pyazr/cli_test.py (JSON_CASES, keyed by the command's name) and no case is left without a row;
#     the reference page docs/source/user_guide/pyazr_command_line.rst is what cli.reference_page() makes of the table.
#     pyazr/cli.py imports nothing of the engine at import time, so it is read here without it.
try:
    import ast
    import importlib.util
    _spec = importlib.util.spec_from_file_location("pyazr_cli_table", os.path.join(ROOT, "pyazr", "cli.py"))
    _cli = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_cli)
    _cmd_rows = [c.name for c in _cli.COMMANDS]
    _cases = None
    for _node in ast.walk(ast.parse(read("tests/pyazr/cli_test.py"))):
        if isinstance(_node, ast.Assign) and any(getattr(t, "id", "") == "JSON_CASES" for t in _node.targets):
            _cases = [k.value for k in _node.value.keys]
    if _cases is None:
        problems.append("tests/pyazr/cli_test.py: no JSON_CASES dict (check 42)")
    else:
        for _name in sorted(set(_cmd_rows) - set(_cases)):
            problems.append(f"pyazr/cli.py: the command '{_name}' has no case in JSON_CASES of tests/pyazr/cli_test.py "
                            "(check 42)")
        for _name in sorted(set(_cases) - set(_cmd_rows)):
            problems.append(f"tests/pyazr/cli_test.py: the case '{_name}' is no row of COMMANDS in pyazr/cli.py "
                            "(check 42)")
    for _c in _cli.COMMANDS:
        for _a in _c.args:
            if _a.type not in _cli.ARG_TYPES:
                problems.append(f"pyazr/cli.py: {_c.name} {_a.name}: type '{_a.type}' is not in ARG_TYPES (check 42)")
    if _cli.reference_page() + "\n" != read("docs/source/user_guide/pyazr_command_line.rst"):
        problems.append("docs/source/user_guide/pyazr_command_line.rst: not what pyazr.cli.reference_page() makes of "
                        "the table; regenerate it (pyazr/CLAUDE.md, recipe 'A command of the command line') (check 42)")
except (OSError, ValueError, AttributeError) as e:
    problems.append(f"check 42: {e}")

# 44. The class scheme of the Developer Guide names types that exist (docs/source/developer/classes.rst, drawn from
#     docs/source/developer/figures/classes.tsv by make_figures.py): every row's type is declared in the file the row
#     names -- a class or struct by "class X" / "struct X" (X the last component of a qualified name), a CMake target
#     by add_library / add_executable / pybind11_add_module of X, a package by its file -- so that a renamed or
#     deleted type fails here until the scheme is updated.  Every arrow of arrows.tsv joins two boxes of its figure.
_scheme_rows = []
try:
    _scheme_lines = [l for l in read("docs/source/developer/figures/classes.tsv").splitlines() if not l.startswith("#")]
    _scheme_head = _scheme_lines[0].split("\t")
    _scheme_rows = [dict(zip(_scheme_head, l.split("\t"))) for l in _scheme_lines[1:] if l.strip()]
    for _r in _scheme_rows:
        _where = f"docs/source/developer/figures/classes.tsv: {_r['type']} ({_r['figure']})"
        if not os.path.isfile(os.path.join(ROOT, _r["file"])):
            problems.append(f"{_where}: {_r['file']} does not exist (check 44)")
            continue
        _short = re.escape(_r["type"].split("::")[-1])
        if _r["kind"] == "class":
            _decl = re.compile(r"^\s*(?:template\s*<[^>]*>\s*)?(?:class|struct)\s+(?:\w+\s+)?" + _short +
                               r"\b(?!\s*;)", re.M)
        elif _r["kind"] == "target":
            _decl = re.compile(r"(?:add_library|add_executable|pybind11_add_module)\s*\(\s*" + _short + r"\s")
        elif _r["kind"] == "package":
            continue
        else:
            problems.append(f"{_where}: unknown kind '{_r['kind']}' (check 44)")
            continue
        if not _decl.search(read(_r["file"])):
            problems.append(f"{_where}: no declaration of {_r['type']} in {_r['file']}; update the scheme "
                            "(classes.tsv, arrows.tsv, make_figures.py) (check 44)")
    _scheme_boxes = {(_r["figure"], _r["type"]) for _r in _scheme_rows}
    _arrow_lines = [l for l in read("docs/source/developer/figures/arrows.tsv").splitlines() if not l.startswith("#")]
    _arrow_head = _arrow_lines[0].split("\t")
    for _l in _arrow_lines[1:]:
        _a = dict(zip(_arrow_head, _l.split("\t")))
        for _end in (_a.get("from", ""), _a.get("to", "")):
            if (_a.get("figure", ""), _end) not in _scheme_boxes:
                problems.append(f"docs/source/developer/figures/arrows.tsv: {_a.get('from')} -> {_a.get('to')}: "
                                f"{_end} is no box of the figure {_a.get('figure')} in classes.tsv (check 44)")
except (OSError, IndexError, KeyError) as e:
    problems.append(f"check 44: {e}")

# 46. The engine's fitters are one table (src/Fitters.cpp, Fitters(); docs/upgrade/fit.md): the program parses the
#     fitters' switches from it (ApplyFitterSwitch) and prints their usage lines from it (FitterSwitchesHelp), and
#     src/AZURE2.cpp names no row's switch itself (no second list); every row's switch, and the switch of a fitter the
#     build may leave out (AbsentFitters()), is a row of docs/source/reference/command_line.rst; every row's name is
#     in the fitter table of docs/source/user_guide/fitting.rst; a fitter is chosen only through the table: include/
#     Config.h declares no fitter bit (no enumerator USE_*_MINIMIZER or USE_ANALYTIC_GRADIENT: the choice is the row's
#     name, Config::fitter, docs/upgrade/fit.md section 5) and Config::fitter is assigned in src/Fitters.cpp
#     (ChooseFitter) and Config::Reset only.  The prose that lists the fitters beside that table names every row too: the
#     section "The engine's own fitters" of user_guide/pyazr.rst (the method), the page user_guide/web_interface.rst
#     (the job method ``engine_<name>``) and the section "Minimizer Selection" of user_guide/calculate.rst (the switch).
_fitter_rows = []
try:
    _fitter_rows = [(_n, _o, "row") for _n, _o in
                    re.findall(r'\{"([\w-]+)", "((?:--[\w-]+)?)", (?:true|false),', read("src/Fitters.cpp"))]
    if not _fitter_rows:
        problems.append("src/Fitters.cpp: no row found in Fitters() (check 46)")
    _fitters_src = read("src/Fitters.cpp")
    _fitter_rows += [(_n, _o, "absent") for _n, _o in
                     re.findall(r'\{"([\w-]+)", "(--[\w-]+)", "[^"]*"\}', _fitters_src)]
    _parser = read("src/AZURE2.cpp")
    for _call in ("ApplyFitterSwitch(configure, *it)", "FitterSwitchesHelp()"):
        if _call not in _parser:
            problems.append(f"src/AZURE2.cpp: {_call} is not called: the fitters' switches and their usage lines come "
                            "from the fitter table (check 46)")
    _cli_page = read("docs/source/reference/command_line.rst")
    _fit_page = read("docs/source/user_guide/fitting.rst")
    _fit_table = _fit_page.split("The Engine's Fitters", 1)[1] if "The Engine's Fitters" in _fit_page else ""
    if not _fit_table:
        problems.append("docs/source/user_guide/fitting.rst: no section \"The Engine's Fitters\" (check 46)")
    for _name, _option, _bit in _fitter_rows:
        if _option and f'"{_option}"' in _parser:
            problems.append(f"src/AZURE2.cpp: names the switch {_option} of the fitter '{_name}' itself; the parser "
                            "and the usage text take the fitters' switches from src/Fitters.cpp (check 46)")
        if _option and f"``{_option}``" not in _cli_page:
            problems.append(f"docs/source/reference/command_line.rst: the switch {_option} of the fitter '{_name}' is "
                            "not listed (check 46)")
        if _bit != "absent" and _fit_table and f"``{_name}``" not in _fit_table:
            problems.append(f"docs/source/user_guide/fitting.rst: the fitter '{_name}' is not in the table of "
                            "\"The Engine's Fitters\" (check 46)")
    def _prose_section(_text, _title):
        """The text from the heading `_title` to the next heading (a line followed by an underline of =-~^)."""
        _heading = re.compile(r"^(\S[^\n]*)\n([=\-~^]{4,})$", re.M)
        _found = [m for m in _heading.finditer(_text)]
        for _k, _m in enumerate(_found):
            if _m.group(1).strip() == _title:
                return _text[_m.end():_found[_k + 1].start() if _k + 1 < len(_found) else len(_text)]
        return ""
    _prose = [("docs/source/user_guide/pyazr.rst", "The engine's own fitters", lambda n, o: f"``{n}``"),
              ("docs/source/user_guide/web_interface.rst", None, lambda n, o: f"``engine_{n}``"),
              ("docs/source/user_guide/calculate.rst", "Minimizer Selection", lambda n, o: f"``{o}``" if o else "")]
    for _page, _title, _form in _prose:
        _text = read(_page)
        _scope = _prose_section(_text, _title) if _title else _text
        if not _scope:
            problems.append(f"{_page}: no section \"{_title}\" (check 46)")
            continue
        for _name, _option in {_n: _o for _n, _o, _b in _fitter_rows}.items():  # a row once, absent or not
            _word = _form(_name, _option)
            if _word and _word not in _scope:
                problems.append(f"{_page}: the fitter '{_name}' of Fitters() is not named ({_word}) in "
                                f"{'the section \"' + _title + '\"' if _title else 'the page'} (check 46)")
    for _bit in re.findall(r"(?m)^\s*(USE_\w+_MINIMIZER|USE_ANALYTIC_GRADIENT)\s*=", read("include/Config.h")):
        problems.append(f"include/Config.h: the bit {_bit} selects a fitter; the chosen fitter is the row's name, "
                        "Config::fitter, set by ChooseFitter (docs/upgrade/fit.md section 5) (check 46)")
    for _rel in tracked_files():
        if not _rel.endswith((".cpp", ".h")) or _rel.startswith("external/") or _rel in ("src/Fitters.cpp",
                                                                                         "src/Config.cpp"):
            continue
        if re.search(r"\bconfig\w*(?:\(\))?\s*(?:\.|->)\s*fitter\s*(?:=(?!=)|\+=|\.(?:assign|clear|swap)\()",
                     read(_rel)):  # a Config by its usual names: configure, configure(), config_, config
            problems.append(f"{_rel}: sets Config::fitter itself; a fitter is chosen by ChooseFitter "
                            "(src/Fitters.cpp), which keeps the table's precedence (check 46)")
except (OSError, IndexError) as e:
    problems.append(f"check 46: {e}")

# 48. The wheel carries every file pyazr reads at run time (docs/upgrade/pyazr3.md section 2.9; pyazr/datafiles.py): every
#     tracked file under pyazr/ that is not Python is a data file the wheel carries (DATA_FILES) or a file it leaves
#     out (NOT_SHIPPED, with the reason), never both; every DATA_FILES entry exists; the patterns of NOT_SHIPPED are
#     the wheel.exclude list of pyproject.toml, in the same order.  tests/pyazr/wheel_test.py compares a built wheel
#     with the list these make.
_wheel_data = 0
try:
    _spec = importlib.util.spec_from_file_location("pyazr_datafiles", os.path.join(ROOT, "pyazr", "datafiles.py"))
    _datafiles = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_datafiles)
    for _rel in tracked_files():
        if not _rel.startswith("pyazr/") or _rel.endswith(".py"):
            continue
        _rel = _rel[len("pyazr/"):]
        if _rel.startswith("_azure2"):
            continue  # the engine module of an in-tree build (no file of git; seen where git cannot list: Windows CI)
        _in, _out = _datafiles.shipped(_rel), _datafiles.excluded(_rel)
        if _in and _out:
            problems.append(f"pyazr/{_rel}: both in DATA_FILES and in NOT_SHIPPED of pyazr/datafiles.py (check 48)")
        elif not (_in or _out):
            problems.append(f"pyazr/{_rel}: not Python, and neither a data file of the wheel (DATA_FILES) nor left "
                            "out with a reason (NOT_SHIPPED) in pyazr/datafiles.py (check 48)")
        _wheel_data += _in
    for _name, _ in _datafiles.DATA_FILES:
        _full = _datafiles.path(_name)
        if not (os.path.isfile(_full) or (_name.endswith("/") and os.path.isdir(_full) and os.listdir(_full))):
            problems.append(f"pyazr/datafiles.py: DATA_FILES names pyazr/{_name}, which does not exist (check 48)")
    _exclude = re.search(r"^wheel\.exclude\s*=\s*\[([^\]]*)\]", read("pyproject.toml"), re.M)
    _patterns = re.findall(r'"([^"]*)"', _exclude.group(1)) if _exclude else None
    if _patterns != [_p for _p, _ in _datafiles.NOT_SHIPPED]:
        problems.append("pyproject.toml: wheel.exclude is not the patterns of NOT_SHIPPED in pyazr/datafiles.py "
                        f"({_patterns} against {[_p for _p, _ in _datafiles.NOT_SHIPPED]}) (check 48)")
except (OSError, AttributeError, ImportError) as e:
    problems.append(f"check 48: {e}")

# 50. No pipe into a reader that leaves early, in a shell script of tests/ or tools/ that sets pipefail
#     (tests/CLAUDE.md, "Writing one"): `printf '%s\n' "$list" | grep -q X` or `cmd | head -1` passes or fails at
#     random on a loaded machine.  The reader (grep -q, -m, -l, -L; head; sed with a q command; awk with an exit
#     outside END) leaves at its first match, the writer's next write gets SIGPIPE, and under pipefail the pipeline's
#     status is then 141: a match reads as "not found", and under `!` a miss as a pass.  A variable is searched with
#     a here-string (`grep -q X <<< "$list"`), a command's output is captured first or read to its end (`sed -n 1p`,
#     a flag in awk).  head -c is a cap, not a race: it reads its whole input up to the cap.  tests/lib/ is sourced
#     by scripts that set pipefail and counts as setting it.
SHELL_TOKEN = re.compile(r"""'[^']*'|"(?:\\.|[^"\\])*"|[;|&)]|[^\s;|&)]+""")  # a word, or what ends a command


def shell_logical_lines(text):
    """(first line number, the line with its continuations and open quotes) of a shell script, comments dropped."""
    lines, start, buf, quote = [], None, "", None
    for number, raw in enumerate(text.split("\n"), 1):
        cut, escaped = len(raw), False
        for i, c in enumerate(raw):
            if escaped:
                escaped = False
            elif c == "\\" and quote != "'":
                escaped = True
            elif quote:
                quote = None if c == quote else quote
            elif c in "'\"":
                quote = c
            elif c == "#" and (i == 0 or raw[i - 1] in " \t;("):
                cut = i
                break
        line = raw[:cut].rstrip()
        start = number if start is None else start
        if quote:  # a quoted text (an awk program) goes on on the next line
            buf += raw[:cut] + "\n"
            continue
        if line.endswith("\\") or re.search(r"(?:\|\||&&|\|)$", line):
            buf += line.rstrip("\\") + " "
            continue
        lines.append((start, buf + line))
        start, buf = None, ""
    return lines


def shell_pipes(line):
    """The text after each pipe '|' of a line (not '||', not '>|'), also inside $( ) within double quotes."""
    stack, i, after = [["", 0]], 0, []  # frames: quote ("", "'", '"') and open parentheses
    while i < len(line):
        c, frame = line[i], stack[-1]
        if frame[0] == "'":
            if c == "'":
                stack.pop()
        elif frame[0] == '"':
            if c == "\\":
                i += 1
            elif c == '"':
                stack.pop()
            elif line.startswith("$(", i):
                stack.append(["", 0])
                i += 1
        elif c == "\\":
            i += 1
        elif c in "'\"":
            stack.append([c, 0])
        elif c == "(":
            frame[1] += 1
        elif c == ")":
            if frame[1]:
                frame[1] -= 1
            elif len(stack) > 1:
                stack.pop()
        elif c == "|" and line.startswith("||", i):
            i += 1
        elif c == "|" and (i == 0 or line[i - 1] != ">"):
            after.append(line[i + 1:])
        i += 1
    return after


def early_reader(command):
    """The name of the reader that leaves before the end of its input, or None."""
    words = SHELL_TOKEN.findall(command)
    words = words[:next((k for k, w in enumerate(words) if w in ";|&)"), len(words))]
    while words and re.match(r"[A-Za-z_]\w*=", words[0]):
        words = words[1:]
    if not words:
        return None
    name, args = os.path.basename(words[0]), words[1:]
    if name in ("grep", "egrep", "fgrep"):
        skip = False
        for a in args:
            if skip:
                skip = False
                continue
            if a == "--" or not a.startswith("-"):
                return None
            if a in ("--quiet", "--silent", "--files-with-matches", "--files-without-match") or \
                    a.startswith("--max-count"):
                return f"{name} {a}"
            if a.startswith("--"):
                continue
            for k, ch in enumerate(a[1:], 1):
                if ch in "qmlL":
                    return f"{name} -{ch}"
                if ch in "ef":  # the rest of the word, or the next word, is a pattern or its file
                    skip = k == len(a) - 1
                    break
        return None
    if name == "head":
        return None if any(a.startswith("-c") for a in args) else "head"
    if name == "sed" and re.search(r"(?:^|['\";{}/0-9$\s])[qQ]\s*\d*\s*(?:['\";}]|$)", " ".join(args)):
        return "sed q"
    if name in ("awk", "gawk", "mawk") and re.search(r"\bexit\b", re.sub(r"\bEND\s*\{.*", "", " ".join(args), flags=re.S)):
        return "awk exit"
    return None


def early_pipes(text):
    """(line number, reader) of every pipe into an early reader."""
    return [(number, reader) for number, line in shell_logical_lines(text) for after in shell_pipes(line)
            for reader in [early_reader(after)] if reader]


PIPEFAIL = re.compile(r"^\s*set\s+-\w*o\s+pipefail\b", re.M)
_pipefail_scripts = 0
for _top in ("tests", "tools"):
    for _dirpath, _dirnames, _files in os.walk(os.path.join(ROOT, _top)):
        _dirnames[:] = sorted(d for d in _dirnames if not d.startswith((".", "build")))
        for _name in sorted(_files):
            if not _name.endswith(".sh"):
                continue
            _rel = os.path.relpath(os.path.join(_dirpath, _name), ROOT).replace(os.sep, "/")
            _text = read(_rel)
            if not (PIPEFAIL.search(_text) or _rel.startswith("tests/lib/")):
                continue
            _pipefail_scripts += 1
            for _number, _reader in early_pipes(_text):
                problems.append(f"{_rel}:{_number}: a pipe into '{_reader}', which leaves before the end of its input; "
                                "under pipefail the writer's SIGPIPE (141) decides the test at random: search a "
                                "variable with a here-string (grep -q X <<< \"$v\"), capture a command's output first "
                                "(check 50)")
_EARLY = ["printf '%s\\n' \"$list\" | grep -q x", "cmd | grep -Fqx 2", "ok \"$(echo \"$m\" | grep -m1 E)\"",
          "! cmd |\n  grep -E -q x", "x=\"$(cmd | head -1)\"", "cmd | head -n 2", "cmd | sed -n '/a/{p;q;}'",
          "cmd | awk '{ print $2; exit }'", "cmd | grep --quiet x", "cmd | LC_ALL=C grep -qF x"]
_LATE = ["grep -q x <<< \"$list\"", "cmd 2>&1 | head -c 1000 > log", "cmd | awk '{ n++ } END { exit !n }'",
         "a || grep -q x f", "cmd | grep -c x", "cmd | sed -n 1p", "awk '{ print \"a|grep -q\" }' f | sort",
         "cmd | grep -e -q", "# cmd | grep -q x", "cmd | sed 's/q/Q/'"]
if [len(early_pipes(s)) for s in _EARLY] != [1] * len(_EARLY) or any(early_pipes(s) for s in _LATE):
    problems.append("check 50: the early-reader scan no longer finds the forms it guards (or finds a safe one): "
                    f"{[s for s in _EARLY if len(early_pipes(s)) != 1]} {[s for s in _LATE if early_pipes(s)]}")

# 56. A comparison of a shell test can fail (tests/CLAUDE.md, "Writing one"; tests/consistency/comparison_scan.py):
#     in tests/*/check.sh and tests/lib/*.sh no awk program prints a mismatch mark from its body and then a number
#     from END (the caller reads both lines and never sees the mark), and no awk program puts an upper bound on a -v
#     value it never tests to be a number (an empty value or a mark compares as text, or after `v + 0` as 0, and
#     passes).  The scan's own forms (FINDS, SAFE) are run first, so a change that blinds it fails here.
import comparison_scan  # noqa: E402

_missed, _flagged = comparison_scan.self_test()
if _missed or _flagged:
    problems.append(f"check 56: the comparison scan misses {_missed} or flags the safe forms {_flagged}; update "
                    "tests/consistency/comparison_scan.py")
for _rel, _line, _kind, _detail in comparison_scan.scan_tree(ROOT):
    problems.append(f"{_rel}:{_line}: {_detail} ({_kind}): a comparison whose failure cannot reach the caller; "
                    "print the mark alone (END: if (flag) exit) or test the value (v ~ /^ *[-+]?\\.?[0-9]/ && ...) "
                    "(check 56)")
_comparison_scripts = len(comparison_scan.scanned_files(ROOT))

# 57. The web page's icon is the original program's (docs/upgrade/web_ui.md section 7): pyazr/ui/static/azure2_icon.png,
#     which the wheel carries, is byte-identical to packaging/icons/azure2_icon.png (no second drawing of it).
_page_icon, _original_icon = "pyazr/ui/static/azure2_icon.png", "packaging/icons/azure2_icon.png"
try:
    with open(os.path.join(ROOT, _page_icon), "rb") as _a, open(os.path.join(ROOT, _original_icon), "rb") as _b:
        if _a.read() != _b.read():
            problems.append(f"{_page_icon} differs from {_original_icon}: copy the original's icon (check 57)")
except OSError as e:
    problems.append(f"check 57: {e}")

# 58. Every example is run and listed (docs/source/user_guide/pyazr_examples.rst): each script of pyazr/examples has
#     a row in ROWS of tests/pyazr/examples_test.py (the table that runs it) and a row of the documentation's list,
#     and neither names a script that is not there; no page of the repository's skills names an example that is not
#     there.  A new example is the file and its two rows.
_examples_page, _examples_test = "docs/source/user_guide/pyazr_examples.rst", "tests/pyazr/examples_test.py"
try:
    import ast
    _example_files = {f[:-3] for f in os.listdir(os.path.join(ROOT, "pyazr", "examples"))
                      if f.endswith(".py") and not f.startswith("_")}
    _example_rows = set()
    for _node in ast.walk(ast.parse(read(_examples_test))):
        if isinstance(_node, ast.Call) and getattr(_node.func, "id", "") == "Row" and _node.args \
                and isinstance(_node.args[0], ast.Constant):
            _example_rows.add(_node.args[0].value)
    _page = read(_examples_page)
    _example_listed = set(re.findall(r"^ +\* - ``([A-Za-z0-9_]+)\.py``", _page, re.M))
    for _name in sorted(_example_files - _example_rows):
        problems.append(f"pyazr/examples/{_name}.py has no row in ROWS of {_examples_test} (check 58)")
    for _name in sorted(_example_rows - _example_files):
        problems.append(f"{_examples_test}: the row '{_name}' names no script of pyazr/examples (check 58)")
    for _name in sorted(_example_files - _example_listed):
        problems.append(f"pyazr/examples/{_name}.py has no row in the list of {_examples_page} (check 58)")
    for _name in sorted(_example_listed - _example_files):
        problems.append(f"{_examples_page}: the list names {_name}.py, which is not in pyazr/examples (check 58)")
    # The repository's skills name examples too: in a paragraph, list item or table row of a skill's page that names
    # pyazr/examples/ (or says "Examples:"), and in the table of azure2-calculate/references/examples.md, every
    # script named (`pyazr/examples/x.py` or a bare `x.py`) is a file of pyazr/examples.
    for _md in sorted(glob.glob(os.path.join(ROOT, ".claude", "skills", "**", "*.md"), recursive=True)):
        _rel = os.path.relpath(_md, ROOT)
        _whole = _rel.endswith(os.path.join("azure2-calculate", "references", "examples.md"))
        for _item in re.split(r"\n\s*\n|\n(?=\s*- )|\n(?=\|)", read(_rel)):  # paragraphs, list items, table rows
            if not (_whole or "pyazr/examples/" in _item or "Examples:" in _item):
                continue
            for _token in re.findall(r"`((?:pyazr/examples/)?[A-Za-z0-9_]+\.py)`", _item):
                _name = _token.rsplit("/", 1)[-1][:-3]
                if _name not in _example_files:
                    problems.append(f"{_rel}: names the example {_token}, which is not in pyazr/examples (check 58)")
except (OSError, ValueError, SyntaxError) as e:
    problems.append(f"check 58: {e}")

# 59. The compound nuclei as set up are a CompoundSystem (include/CompoundSystem.h; docs/upgrade/cnuc_split.md section
#     12): one constant object shared by every evaluation, never copied.  The old name CNuc (in any identifier or file
#     name), a Clone() of the compound and the pools of compounds (cnuc_pool_) do not come back in the engine, the
#     binding or their headers.
OLD_COMPOUND_NAME = re.compile(r"CNuc|cnuc_pool_")
COMPOUND_CLONE = re.compile(r"\bClone\s*\(")
_compound_files = 0
for rel in tracked_files():
    if not rel.startswith(("src/", "include/", "api/")):
        continue
    if "CNuc" in os.path.basename(rel):
        problems.append(f"{rel}: a file named after CNuc; the compound system is CompoundSystem (check 59)")
    if not rel.endswith((".cpp", ".h", ".hpp")):
        continue
    try:
        text = read(rel)
    except (OSError, UnicodeDecodeError):
        continue
    _compound_files += 1
    for number, line in enumerate(text.splitlines(), 1):
        if OLD_COMPOUND_NAME.search(line):
            problems.append(f"{rel}:{number}: names {OLD_COMPOUND_NAME.search(line).group(0)}; the compound system is "
                            "CompoundSystem, shared constant, with a LevelState per owner (check 59)")
    if rel == "include/CompoundSystem.h" and COMPOUND_CLONE.search(text):
        problems.append(f"{rel}: CompoundSystem has a Clone(); it is shared constant, never copied (check 59)")
if not OLD_COMPOUND_NAME.search("const CNuc *c") or not OLD_COMPOUND_NAME.search("std::vector<int> theCNucList") \
        or not OLD_COMPOUND_NAME.search("cnuc_pool_.top()") or OLD_COMPOUND_NAME.search("const CompoundSystem &c") \
        or not COMPOUND_CLONE.search("std::unique_ptr<CompoundSystem> Clone() const;"):
    problems.append("check 59: the old-name pattern no longer matches the forms it guards; update this check")

if problems:
    for p in problems:
        print("FAIL", p)
    print(f"{len(problems)} consistency problem(s)")
    sys.exit(1)
print(f"repository consistent: {len(flags)} CLI flags, test registration, check.sh modes, pyazr version and exports, "
      f"format-3 reference page, {len(declared)} diagnostic codes, file sizes, one model pass and residual, one calculation plan, "
      f"no machine paths, no old data class names, no old section in the engine, reference leaves by one predicate, "
      f"the sharing key's term settings, {len(_allowed)} allowed state rows, one reader of the environment "
      f"({_getenv_files} files), {len(_table)} environment switches in the manual, {len(_names)} switch options, "
      f"{len(_rows)} physics rows tested, {len(_cli.COMMANDS)} pyazr commands, {len(_scheme_rows)} boxes of the class scheme, "
      f"{len(_fitter_rows)} fitters, {_wheel_data} data files of the wheel, "
      f"{_pipefail_scripts} shell scripts under pipefail without a pipe into an early reader, "
      f"{_comparison_scripts} test scripts whose comparisons can fail, the page's icon the original's, "
      f"{len(_example_files)} examples run and listed, no CNuc in {_compound_files} engine files")
