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
SIZE_DIRS = ("src", "include", "physics", "io", "experiment", "api/src", "pyazr", "tools", "stopping")
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
parser = function_body(main_cpp, "bool parseOptions(")
help_text = function_body(main_cpp, "void printHelp(")
cli_doc = read("docs/source/reference/command_line.rst")

flags = set(re.findall(r'"(--[a-z0-9-]+)=?"', parser))
# Handled in main() before parseOptions, so accepted although parseOptions
# never sees them.
handled_in_main = {"--help", "--build-coulomb-table", "--convert", "--schema", "--validate", "--describe",
                   "--print-predictions",
                   # options of --convert ... --to 3 and of --validate, read in main() with their command
                   "--to", "--json-report", "--answers", "--exact-grids", "--strict", "--allow-questions", "--json",
                   # option of --describe
                   "--record"}
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
            if inc.startswith("<") and re.fullmatch(r"[a-z_]+", inc[1:-1]):
                continue
            problems.append(f"{rel}: includes {inc}; physics/ may include only azure2/physics/ headers and the C++ standard library")


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
#     physics/'s ("azure2/physics/..."), the C++ standard library, and the bare names of its own files where they live
#     (src/'s generated-table header, the tests' check.h).
STOPPING_OWN = {"src": {"TableText.h"}, "tests": {"check.h"}}
for dirpath, _, files in os.walk(os.path.join(ROOT, "stopping")):
    for name in sorted(files):
        if not name.endswith((".h", ".hpp", ".cpp", ".in")):
            continue
        rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
        part = rel.split(os.sep)[1] if len(rel.split(os.sep)) > 2 else ""
        for inc in re.findall(r'^\s*#\s*include\s*([<"][^>"]+[>"])', read(rel), re.M):
            if inc.startswith('"') and inc[1:-1].startswith(("azure2/stopping/", "azure2/physics/")):
                continue
            if inc.startswith('"') and inc[1:-1] in STOPPING_OWN.get(part, set()):
                continue
            if inc.startswith("<") and re.fullmatch(r"[a-z_]+", inc[1:-1]):
                continue
            problems.append(f"{rel}: includes {inc}; stopping/ may include only azure2/stopping/ and azure2/physics/ "
                            "headers, its own files and the C++ standard library (check 32)")



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
#      which pyazr.nuclei.SYMBOLS is).  A second list in the engine, the API, io/, experiment/ or pyazr -- three
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
PENALTIES = ("ParkNormPenalty", "AbundancePenalty", "AbundanceSumPenalty", "ScatteringLength::Chi2",
             "ScatteringLength::Residuals", "ScatteringLength::Chi2Gradient", "ModelData::Chi2", "ModelData::Residuals",
             "ModelData::Gradient", "ModelData::AppendRows",
             "ScreeningPenalty", "DensityPenalty", "DopplerPenalty", "ThicknessPenalty", "CovarianceChi2Correction")
DEFINERS = {"src/Objective.cpp", "src/CNuc.cpp", "src/ScatteringLength.cpp", "src/ModelData.cpp", "src/WorkingCopy.cpp",
            "src/Segment.cpp", "src/TargetEffect.cpp"}
ALLOWED_PENALTY_USES = {
    # The MCMC splits the posterior: per-nucleus abundance priors are log-priors (CalculateLogPrior), so the
    # likelihood carries only the sum-to-one term (and the model data, Objective::ModelDataChi2); Park's wall rejects
    # the step outright.
    ("src/AZURECalcMCMC.cpp", "AbundanceSumPenalty"),
    ("src/AZURECalcMCMC.cpp", "ParkNormPenalty"),
    # The API's data chi^2 in residual form (EvaluateFilledChi2, the analytic gradient): T2 for its rows;
    # pyazr_objective checks it against the CLI's objective.
    ("api/src/AZUREAPI.cpp", "CovarianceChi2Correction"),
    # The output writer, moved out of src/WorkingCopy.cpp (docs/upgrade/edata_split.md N1): chiSquared.out prints the
    # target-effect prior totals, which Results::objective does not hold yet (decision D2 of the note).
    ("src/DataOutputFiles.cpp", "ScreeningPenalty"),
    ("src/DataOutputFiles.cpp", "DensityPenalty"),
    ("src/DataOutputFiles.cpp", "DopplerPenalty"),
    ("src/DataOutputFiles.cpp", "ThicknessPenalty"),
    # The data covariance, moved out of src/WorkingCopy.cpp (docs/upgrade/edata_split.md W2): the file defines the
    # covariance term, and PrepareCovarianceSeeds evaluates it for the seeds the gradient reads.
    ("src/DataCovarianceGroups.cpp", "CovarianceChi2Correction"),
    # The hand-over of a parameter vector, moved out of src/WorkingCopy.cpp (docs/upgrade/edata_split.md M2): the file
    # defines the target-effect prior totals (WorkingCopy::ScreeningPenalty ...), each the sum of its effects' terms.
    ("src/HandOver.cpp", "ScreeningPenalty"),
    ("src/HandOver.cpp", "DensityPenalty"),
    ("src/HandOver.cpp", "DopplerPenalty"),
    ("src/HandOver.cpp", "ThicknessPenalty"),
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
                        "(after AZURE2 --schema > tests/schema/schema_v3.json if the table changed)")
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
#     7.4): what a parameter is comes from its record (include/ParameterKind.h), written where it is created; a name
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
ALLOWED_NAME_TESTS = {
    # An example script that reads param.par without the engine (pyazr.transform only).
    ("pyazr/examples/transform_widths.py", 'm = re.match(r"width_'),
}
name_test_hits = set()
for top in ("src", "include", "api", "pyazr"):
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
kinds_h = read("include/ParameterKind.h")
enum_body = kinds_h[kinds_h.index("enum class ParamKind"):]
enum_body = enum_body[:enum_body.index("};")]
enum_kinds = re.findall(r"^\s*(\w+) = (\d+),?", enum_body, re.M)
row_kinds = re.findall(r'\{ParamKind::(\w+), "(\w+)"', kinds_h)
if [k for k, _ in enum_kinds] != [k for k, _ in row_kinds] or \
        [int(c) for _, c in enum_kinds] != list(range(len(enum_kinds))):
    problems.append("include/ParameterKind.h: ParamKindTable() must have one row per ParamKind, in the order of the "
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
#     buffer.  Unit.h (with UnitTables.h, its tables) and WorkingCopy.h declare no member of a type that names
#     Results -- a value, pointer, reference, container or smart pointer.  Function declarations may name it (the
#     buffer is passed in).
for rel in ("include/Unit.h", "include/UnitTables.h", "include/WorkingCopy.h"):
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
eph = re.sub(r"//[^\n]*|/\*.*?\*/", "", read("include/Unit.h"), flags=re.S)
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
    if not rel.startswith(("src/", "include/", "api/", "tests/", "experiment/", "convert/", "io/", "physics/")):
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

if problems:
    for p in problems:
        print("FAIL", p)
    print(f"{len(problems)} consistency problem(s)")
    sys.exit(1)
print(f"repository consistent: {len(flags)} CLI flags, test registration, check.sh modes, pyazr version and exports, "
      f"format-3 reference page, {len(declared)} diagnostic codes, file sizes, one model pass and residual, one calculation plan, "
      f"no machine paths, no old data class names")
