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
SIZE_DIRS = ("src", "include", "physics", "io", "experiment", "api/src", "pyazr", "tools")
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
handled_in_main = {"--help", "--build-coulomb-table", "--convert", "--schema", "--validate",
                   # options of --convert ... --to 3 and of --validate, read in main() with their command
                   "--to", "--json-report", "--answers", "--exact-grids", "--strict", "--allow-questions", "--json"}
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


# 6. The segment-line keyword blocks are the same set in both .azr readers.
#    pyazr splits the tail by its own table, so a keyword the engine learns and
#    pyazr does not is mangled by a load and save there.
engine_kw = set()
for header in ("convert/include/azure2/convert/SegLine.h", "convert/include/azure2/convert/ExtrapLine.h"):
    engine_kw |= set(re.findall(r'keyword == "([a-z]+)"', read(header)))
py_match = re.search(r"_KEYWORDS\s*=\s*\{(.*?)\}", read("pyazr/datasets.py"), re.S)
py_kw = set(re.findall(r'"([a-z]+)"\s*:', py_match.group(1))) if py_match else set()
for name, kw in (("pyazr/datasets.py _KEYWORDS", py_kw),):
    for k in sorted(engine_kw - kw):
        problems.append(f"segment keyword '{k}': read by the engine (convert/'s SegLine.h/ExtrapLine.h), unknown to {name}")
    for k in sorted(kw - engine_kw):
        problems.append(f"segment keyword '{k}': known to {name}, not read by the engine")
if not engine_kw:
    problems.append("could not find the segment keywords in convert/include/azure2/convert/SegLine.h; update this test")


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



# 9b. experiment/ uses physics/ and io/'s plain record structs only (docs/upgrade/data_model.md section 2.5, decision
#     M1; docs/upgrade/project_format3.md section 7.3, the director's decision of 2026-11-09): its own headers
#     ("azure2/experiment/..."), physics/'s ("azure2/physics/..."), the one header "azure2/io/Record.h" (header-only
#     structs) and the C++ standard library.  Its tests may also include io/'s reader ("azure2/io/..."): they read
#     their fixtures as the engine will.  An engine header (EPoint.h, Config.h), another io/ header in the library
#     or a system library makes the module depend on what it must not.
for dirpath, _, files in os.walk(os.path.join(ROOT, "experiment")):
    for name in sorted(files):
        if not name.endswith((".h", ".hpp", ".cpp")):
            continue
        rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
        in_tests = rel.startswith(os.path.join("experiment", "tests") + os.sep)
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
    # m_e of CODATA 2002 in the beta-decay Fermi-function argument; CODATA 2018 moves beta-delayed outputs (~1e-7).
    ("src/EPoint.cpp", "0.510998903"),
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
DEFINERS = {"src/Objective.cpp", "src/CNuc.cpp", "src/ScatteringLength.cpp", "src/ModelData.cpp", "src/EData.cpp",
            "src/ESegment.cpp", "src/TargetEffect.cpp"}
ALLOWED_PENALTY_USES = {
    # The MCMC splits the posterior: per-nucleus abundance priors are log-priors (CalculateLogPrior), so the
    # likelihood carries only the sum-to-one term (and the model data, Objective::ModelDataChi2); Park's wall rejects
    # the step outright.
    ("src/AZURECalcMCMC.cpp", "AbundanceSumPenalty"),
    ("src/AZURECalcMCMC.cpp", "ParkNormPenalty"),
    # The API's data chi^2 in residual form (EvaluateFilledChi2, the analytic gradient): T2 for its rows;
    # pyazr_objective checks it against the CLI's objective.
    ("api/src/AZUREAPI.cpp", "CovarianceChi2Correction"),
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
    # The <parameterSettings> "segment_<key>_norm|energy_shift prior_centre c" row, read and rewritten by pyazr's
    # document layer (AzrModel.prior_centres / set_prior_centres, SegmentSet) without an engine session, so without a
    # parameter vector to look the name up in.  The engine matches the row's name whole (EData::ReadPriorCentres).
    ("pyazr/datasets.py", "_PRIOR_CENTRE = re.compile("),
    ("pyazr/azrfile.py", "_PRIOR_ROW = re.compile("),
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
#     buffer.  EPoint.h and EData.h declare no member of a type that names Results -- a value, pointer, reference,
#     container or smart pointer.  Function declarations may name it (the buffer is passed in).
for rel in ("include/EPoint.h", "include/EData.h"):
    text = re.sub(r"//[^\n]*|/\*.*?\*/|^\s*#[^\n]*", "", read(rel), flags=re.S | re.M)
    for stmt in re.split(r"[;{}]", text):
        if "(" in stmt or not re.search(r"\bResults\b", stmt):
            continue  # a function (its parameters may name the buffer), or nothing about results
        problems.append(f"{rel}: declares a results member ({' '.join(stmt.split())}); the buffer belongs to the "
                        "caller of the evaluation and is passed in (include/Results.h)")


if problems:
    for p in problems:
        print("FAIL", p)
    print(f"{len(problems)} consistency problem(s)")
    sys.exit(1)
print(f"repository consistent: {len(flags)} CLI flags, segment keywords, test registration, check.sh modes, pyazr version and exports, "
      f"format-3 reference page, {len(declared)} diagnostic codes, file sizes")
