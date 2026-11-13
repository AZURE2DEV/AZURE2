#!/usr/bin/env python3
"""The pyazr command line (pyazr/cli.py, ``python -m pyazr``; docs/upgrade/pyazr3.md section 2.8).

  1. every row of the table of commands runs with --json on T1 or T2 (or a file of the repository) and its JSON is
     read: JSON_CASES has one case per row (tests/consistency/repo_consistency_test.py check 41 fails for a row
     without one); standard output holds exactly one JSON document;
  2. describe and validate give what the program gives (AZURE2 --describe/--validate --json on T2 and on a project
     with errors, with AZURE2_BIN); schema equals AZURE2 --schema;
  3. wrong calls (a missing file, an unknown command or option, a missing or bad argument) exit 2 and an old project
     3, each with its three parts (what, allowed, fix) in both modes, nothing on standard output without --json;
  4. --help of every command names every argument of its row.
Nothing here uses the network: EXFOR is read from tests/pyazr/data/exfor_master, ENSDF from a saved answer of the
NNDC API (tests/pyazr/data/ensdf_api_14N_levels.json) in this process, and the commands run in a subprocess have an
unreachable proxy.  Exits 77 (skipped) without the compiled engine module.
Run from anywhere:  python3 tests/pyazr/cli_test.py
"""
import json
import os
import select
import shutil
import signal
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
T1 = os.path.join(ROOT, "tests", "upgrade", "T1", "T1.azr")
T2 = os.path.join(ROOT, "tests", "upgrade", "T2", "T2.azr")
OLD_T2 = os.path.join(ROOT, "tests", "convert", "fixtures", "tests", "upgrade", "T2", "T2.azr")
EXFOR_MASTER = os.path.join(HERE, "data", "exfor_master")
ENSDF_PAGE = os.path.join(HERE, "data", "ensdf_api_14N_levels.json")
BIN = os.environ.get("AZURE2_BIN")

sys.path.insert(0, ROOT)
try:
    import pyazr
    from pyazr import _azure2  # noqa: F401
    from pyazr import cli
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)

failures = []
OFFLINE = {"http_proxy": "http://127.0.0.1:9", "https_proxy": "http://127.0.0.1:9",
           "HTTP_PROXY": "http://127.0.0.1:9", "HTTPS_PROXY": "http://127.0.0.1:9", "no_proxy": ""}


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:600]}"))
    if not ok:
        failures.append(name)


def env(**extra):
    e = dict(os.environ, OMP_NUM_THREADS="1", **OFFLINE, **extra)
    e["PYTHONPATH"] = ROOT + (os.pathsep + e["PYTHONPATH"] if e.get("PYTHONPATH") else "")
    return e


def run(*args, cwd=None, **extra):
    """``python -m pyazr ARGS``: (exit code, the JSON document when --json was given else the text, stderr)."""
    p = subprocess.run([sys.executable, "-m", "pyazr", *args], capture_output=True, text=True, cwd=cwd or TMP,
                       env=env(**extra), timeout=300)
    if "--json" not in args:
        return p.returncode, p.stdout, p.stderr
    try:
        return p.returncode, json.loads(p.stdout), p.stderr  # one document and nothing else
    except ValueError as e:
        check(f"{' '.join(args)}: standard output is one JSON document", False, f"{e}: {p.stdout[:300]!r}")
        return p.returncode, None, p.stderr


def program(*args):
    p = subprocess.run([BIN, *args], capture_output=True, text=True, env=env(), cwd=TMP)
    return p.returncode, p.stdout


def three_parts(doc):
    err = (doc or {}).get("error") or {}
    return all(isinstance(err.get(k), str) and err.get(k).strip() for k in ("what", "allowed", "fix"))


def ok_doc(name, code, doc, command):
    check(f"{name}: exit 0, ok, command named", code == 0 and doc and doc["ok"] and doc["command"] == command
          and "result" in doc, doc)
    return (doc or {}).get("result")


# ---------------------------------------------------------------------------------------------------------------------
# One case per row of cli.COMMANDS
# ---------------------------------------------------------------------------------------------------------------------

def case_commands():
    code, doc, _ = run("commands", "--json")
    result = ok_doc("commands", code, doc, "commands") or []
    check("commands: every row, in order", [r["name"] for r in result] == [c.name for c in cli.COMMANDS], result)
    check("commands: every argument with its type and help",
          all(len(r["arguments"]) == len(c.args)
              and all(a["type"] in cli.ARG_TYPES and a["help"] for a in r["arguments"])
              for r, c in zip(result, cli.COMMANDS)))


def case_describe():
    code, doc, _ = run("describe", T2, "--json")
    result = ok_doc("describe T2", code, doc, "describe")
    check("describe T2: the data sets", result and len(result["datasets"]) > 0, result)
    if BIN:
        pcode, out = program("--describe", T2, "--json")
        check("describe T2 --json equals AZURE2 --describe T2 --json", pcode == 0 and json.loads(out) == result)
    code, doc, _ = run("describe", T2, "--record", "s01_int_thin", "--json")
    result = ok_doc("describe --record", code, doc, "describe")
    check("describe --record: that data set only",
          result and [d["name"] for d in result["datasets"]] == ["s01_int_thin"])
    code, doc, _ = run("describe", "--json")
    result = ok_doc("describe (the package)", code, doc, "describe")
    check("describe (the package): version and engine", result and result["pyazr"] == pyazr.__version__
          and "engine" in result, result)
    code, text, err = run("describe", T2)
    check("describe T2 (text): the engine's text", code == 0 and text.startswith("dataset s01_int_thin"), text[:200])


BROKEN = """<config>
format=3
</config>
<data>
dataset d1
    reaction       15N(p,a)12C
    quantity       sigma
    unit           mb/sr
    file           d1.dat
    columns        E_lab[MeV] value error
    norm           1 fixed
    taget          t1
</data>
"""


def case_validate():
    code, doc, _ = run("validate", T2, "--json")
    check("validate T2: exit 0, no diagnostic", ok_doc("validate T2", code, doc, "validate") == [], doc)
    if BIN:
        pcode, out = program("--validate", T2, "--json")
        check("validate T2 --json equals AZURE2 --validate T2 --json",
              pcode == 0 and doc and json.loads(out) == doc["result"])
    broken = os.path.join(TMP, "broken.azr")
    with open(broken, "w") as f:
        f.write(BROKEN)
    code, doc, _ = run("validate", broken, "--json")
    check("validate a project with errors: exit 2, the diagnostics, three parts",
          code == 2 and doc and not doc["ok"] and len(doc["result"]) >= 2 and three_parts(doc), doc)
    if BIN:
        pcode, out = program("--validate", broken, "--json")
        check("validate broken: the program's exit and diagnostics",
              pcode == code and doc and json.loads(out) == doc["result"], out)
    code, text, err = run("validate", broken)
    check("validate broken (text): file:line: record.field lines, the summary, error on stderr",
          code == 2 and f"{broken}:12: d1.taget: E_UNKNOWN_FIELD" in text and "not valid." in text and "fix:" in err,
          text + err)


def case_convert():
    folder = shutil.copytree(os.path.dirname(T2), os.path.join(TMP, "T2old"))
    os.remove(os.path.join(folder, "prior_cov.f3.txt"))  # the conversion writes its own beside the result
    old, new = os.path.join(folder, "T2.azr"), os.path.join(folder, "T2_f3.azr")
    shutil.copy(OLD_T2, old)
    code, doc, _ = run("convert", old, new, "--json")
    result = ok_doc("convert the old T2", code, doc, "convert")
    check("convert: the report, the parameters renamed", result and result["ok"] and result["parameters"], result)
    code, doc, _ = run("validate", new, "--json")
    check("convert: the result validates", code == 0 and doc["result"] == [], doc)
    code, doc, _ = run("convert", old, new, "--json")
    check("convert onto an existing file: exit 2, three parts, --overwrite named",
          code == 2 and three_parts(doc) and "--overwrite" in doc["error"]["fix"], doc)


def case_schema():
    code, doc, _ = run("schema", "--json")
    result = ok_doc("schema", code, doc, "schema")
    check("schema equals pyazr.schema()", result == json.loads(json.dumps(pyazr.schema())))
    if BIN:
        pcode, out = program("--schema")
        check("schema equals AZURE2 --schema", pcode == 0 and json.loads(out) == result)
    code, doc, _ = run("schema", "dataset", "--json")
    result = ok_doc("schema dataset", code, doc, "schema")
    check("schema dataset: the record kind", result and result["table"] == "records" and "fields" in result["entry"])
    code, doc, _ = run("schema", "datset", "--json")
    check("schema datset: exit 2, the close name", code == 2 and "'dataset'" in doc["error"]["what"], doc)


def case_kinds():
    code, doc, _ = run("kinds", "--json")
    result = ok_doc("kinds", code, doc, "kinds")
    check("kinds: the tables of pyazr.kinds()", result and sorted(result) == sorted(pyazr.kinds()),
          result and list(result))
    code, doc, _ = run("kinds", "options", "--json")
    result = ok_doc("kinds options", code, doc, "kinds")
    check("kinds options: the session options", result and "use_brune" in result, result)


def case_stopping():
    code, doc, _ = run("stopping", "Ta2O5", "p", "0.5", "1.0", "--json")
    result = ok_doc("stopping", code, doc, "stopping")
    ref = pyazr.stopping.stopping("Ta2O5", "p", [0.5, 1.0])
    check("stopping: pyazr.stopping.stopping's numbers, the citation",
          result and result["total"] == ref["total"].tolist() and result["energy_MeV"] == [0.5, 1.0]
          and result["citation"] == pyazr.stopping.CITATION, result)
    code, doc, _ = run("stopping", "Ar", "a", "1", "--phase", "solid", "--json")
    check("stopping --phase solid", code == 0 and doc["result"]["phase"] == "solid", doc)


def case_coulomb_check():
    code, doc, _ = run("coulomb-table", "check", T1, "--json")
    check("coulomb-table check on a file that is no table: exit 1, not loaded, three parts",
          code == 1 and doc["result"]["loaded"] is False and doc["result"]["in_use"] == "exact" and three_parts(doc),
          doc)
    table = os.environ.get("AZURE2_COULOMB_TABLE")
    if table and os.path.isfile(table):
        code, doc, _ = run("coulomb-table", "check", "--json")
        result = ok_doc("coulomb-table check $AZURE2_COULOMB_TABLE", code, doc, "coulomb-table check")
        check("coulomb-table check: the path and the content identity",
              result and result["loaded"] and result["path"] == table and result["identity"], result)
    else:
        print("  skip  coulomb-table check of a real table: AZURE2_COULOMB_TABLE is not set to a file")


def _not_available(words):
    code, doc, _ = run(*words, "--json")
    check(f"{' '.join(words)}: exit 2, not available, three parts",
          code == 2 and doc["error"]["class"] == "NotAvailable" and three_parts(doc), doc)


def case_ui():
    try:
        import fastapi  # noqa: F401
        import uvicorn  # noqa: F401
    except ImportError as e:
        print(f"  skip  ui: {e}")
        return
    folder = shutil.copytree(os.path.dirname(T1), os.path.join(TMP, "T1ui"))
    p = subprocess.Popen([sys.executable, "-m", "pyazr", "ui", os.path.join(folder, "T1.azr"), "--json"],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=TMP, env=env())
    try:
        ready, _, _ = select.select([p.stdout], [], [], 60)
        line = p.stdout.readline() if ready else ""
        doc = json.loads(line) if line else None
        check("ui --json: the document with the URL as soon as it listens",
              doc and doc["ok"] and doc["command"] == "ui" and doc["result"]["url"].startswith("http://127.0.0.1:"),
              line)
    finally:
        p.send_signal(signal.SIGTERM)
        try:
            p.wait(30)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
    check("ui: nothing more on standard output", p.stdout.read() == "")


def case_import_exfor():
    code, doc, _ = run("import", "exfor", "T0001002", "--source", "local", "--json", AZURE2_EXFOR_MASTER=EXFOR_MASTER)
    result = ok_doc("import exfor", code, doc, "import exfor")
    check("import exfor: reaction, author, the values",
          result and result["reaction"] == "6-C-12(A,G)8-O-16,,SIG,,SFC" and result["author"] == "A.Author"
          and result["arrays"]["DATA"][:2] == [0.000123, 0.00025], result)


def in_process(argv):
    """cli.main(argv) in this process: (exit code, standard output as written to descriptor 1)."""
    with tempfile.TemporaryFile("w+") as f:
        sys.stdout.flush()
        saved = os.dup(1)
        os.dup2(f.fileno(), 1)
        try:
            code = cli.main(argv)
        finally:
            sys.stdout.flush()
            os.dup2(saved, 1)
            os.close(saved)
        f.seek(0)
        return code, f.read()


def case_import_ensdf():
    import socket
    from pyazr import nds
    page = json.load(open(ENSDF_PAGE))

    def web(url, params=None, timeout=60):
        assert url == f"{nds.ENSDF_API_URL}/levels", url
        return json.dumps(page)

    def no_network(*a, **k):
        raise AssertionError("the network was used")

    saved_get, saved_connect = nds._http_get, socket.socket.connect
    nds._http_get, socket.socket.connect = web, no_network
    try:
        code, out = in_process(["import", "ensdf", "14n", "--source", "nndc", "--json"])
        reference = cli.plain(pyazr.fetch_levels("14n", source="nndc"))
    finally:
        nds._http_get, socket.socket.connect = saved_get, saved_connect
    doc = json.loads(out)
    result = ok_doc("import ensdf", code, doc, "import ensdf")
    check("import ensdf: pyazr.fetch_levels's levels", result == reference and len(result) == 3, result)


JSON_CASES = {
    "commands": case_commands,
    "describe": case_describe,
    "validate": case_validate,
    "convert": case_convert,
    "schema": case_schema,
    "kinds": case_kinds,
    "stopping": case_stopping,
    "coulomb-table check": case_coulomb_check,
    "coulomb-table build": lambda: _not_available(["coulomb-table", "build"]),
    "coulomb-table fetch": lambda: _not_available(["coulomb-table", "fetch"]),
    "ui": case_ui,
    "import exfor": case_import_exfor,
    "import ensdf": case_import_ensdf,
}

# (arguments, exit code, error class): wrong calls; each is run with --json and without
WRONG = [
    (["describe", "missing.azr"], 2, "MissingFile"),
    (["frobnicate"], 2, "UsageError"),
    (["coulomb-table"], 2, "UsageError"),
    (["validate"], 2, "UsageError"),
    (["describe", T2, "--recrd", "x"], 2, "UsageError"),
    (["stopping", "C", "p", "x"], 2, "BadArgument"),
    (["stopping", "C", "q", "1"], 2, "BadArgument"),
    (["stopping", "Xx2", "p", "1"], 2, "ValueError"),
    (["import", "exfor", "T0001002", "--plus", "one"], 2, "BadArgument"),
    (["validate", OLD_T2], 3, "FormatError"),
]


def main():
    check("one case per row of the table", sorted(JSON_CASES) == sorted(c.name for c in cli.COMMANDS),
          sorted(set(JSON_CASES) ^ {c.name for c in cli.COMMANDS}))
    for name, case in JSON_CASES.items():
        print(f"{name}:")
        case()

    print("wrong calls:")
    for args, expected, kind in WRONG:
        code, doc, _ = run(*args, "--json")
        check(f"{' '.join(args)} --json: exit {expected}, {kind}, three parts",
              code == expected and doc and not doc["ok"] and doc["exit"] == expected and doc["error"]["class"] == kind
              and three_parts(doc), doc)
        code, out, err = run(*args)
        check(f"{' '.join(args)}: exit {expected}, three parts on stderr, nothing on stdout",
              code == expected and out == "" and "error: " in err and "\n  allowed: " in err and "\n  fix: " in err,
              out + err)

    print("--help:")
    code, out, _ = run("--help")
    check("--help names every command", code == 0 and all(c.name.split()[0] in out for c in cli.COMMANDS), out)
    for c in cli.COMMANDS:
        code, out, _ = run(*c.name.split(), "--help")
        missing = [a.name for a in c.args if a.name not in out] + ([] if "--json" in out else ["--json"])
        check(f"{c.name} --help names every argument", code == 0 and not missing and c.help.split(",")[0] in
              " ".join(out.split()), missing or out)

    print("reference page:")
    page = os.path.join(ROOT, "docs", "source", "user_guide", "pyazr_command_line.rst")
    check("the reference page is the table's", open(page).read() == cli.reference_page() + "\n")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as TMP:
        t0 = time.time()
        main()
        print(f"({time.time() - t0:.0f} s)")
    if failures:
        print(f"FAIL: {len(failures)} check(s): {failures}")
        sys.exit(1)
    print("ok")
