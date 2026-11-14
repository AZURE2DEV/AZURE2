#!/usr/bin/env python3
"""The names the published pages show exist: every program option and every pyazr name of docs/source.

Reads the inline literals (``...``) and the literal blocks (after ``::`` or ``.. code-block::``) of every page of
docs/source/ and lists

  option  an option of the program (a literal that starts with --, or a line that runs AZURE2) that the program's
          --help does not list;
  pyazr   a name pyazr.<a>.<b>... that does not resolve after `import pyazr`, or an attribute of Session, Project,
          Parameters, Objective or Evaluation (written Class.attr, or s.attr / p.attr for a session / project) that
          the class lacks.

The table "Earlier names" of docs/source/user_guide/pyazr.rst lists removed names on purpose and is skipped.  The
program and pyazr come from the environment the script runs in (AZURE2_BIN, PYTHONPATH and PYAZR_ENGINE_DIR):

    AZURE2_BIN=<build>/src/AZURE2 PYTHONPATH=<worktree> PYAZR_ENGINE_DIR=<build>/pyazr-out \\
        python3 tools/check_doc_names.py

Exits 1 when it lists anything.
"""

import importlib
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PAGES = os.path.join(ROOT, "docs", "source")
# a literal or a line that names another program than AZURE2: its options are not the program's
OTHER_PROGRAMS = re.compile(r"\b(pyazr|azure2-coulomb-table|cmake|ctest|pip|git|python3?|sphinx|make|conda|"
                            r"run_tests|check\.sh|runall|cls\.sh|pytest|npm|node|gh|uv|brew|apt|ssh|curl|"
                            r"chrome|scikit|pyproject)\b", re.I)
OPTION = re.compile(r"(?<![\w-])--[a-z][a-z0-9-]*")
PYAZR_NAME = re.compile(r"\bpyazr(?:\.[A-Za-z_]\w*)+")
CLASS_ATTR = re.compile(r"\b(Session|Project|Parameters|Objective|Evaluation)\.([A-Za-z_]\w*)")
INSTANCE_ATTR = re.compile(r"(?<![\w.+/-])(s|p)\.([A-Za-z_]\w*)\b(?!\.)")
FILE_SUFFIXES = ("cpp", "h", "py", "rst", "md", "json", "azr", "out", "sh", "txt")
# names a page shows on purpose although they do not exist, with the reason
SHOWN_ON_PURPOSE = {
    "--use-gradiant": "reference/command_line.rst: the example of a mistyped option",
    "pyazr.cleanup": "user_guide/pyazr.rst: 'no pyazr.cleanup' (what a session does not need)",
}
# pages whose bare ``--option`` literals are the program's (others document pyazr's command line, the Coulomb
# table tool, the web page's server or scripts)
PROGRAM_PAGES = re.compile(r"docs/source/(reference/(command_line|format3|output_files|data_formats)|"
                           r"user_guide/(?!pyazr|web_|graphical)[\w]+)\.rst$")


def literals(text):
    """(line number, text) of every inline literal and literal-block line of an .rst page."""
    out = []
    lines = text.split("\n")
    block_indent = None
    for number, line in enumerate(lines, 1):
        stripped = line.strip()
        if block_indent is not None:
            indent = len(line) - len(line.lstrip())
            if stripped and indent <= block_indent:
                block_indent = None
            elif stripped:
                out.append((number, stripped))
                continue
            else:
                continue
        if stripped.endswith("::") or stripped.startswith(".. code-block::"):
            block_indent = len(line) - len(line.lstrip())
            if stripped.startswith(".. code-block::") or stripped == "::":
                pass
        for m in re.finditer(r"``(.+?)``", line):
            out.append((number, m.group(1)))
    return out


def skip_earlier_names(rel, text):
    """The text with the section "Earlier names" of the pyazr page blanked (line numbers kept)."""
    if not rel.endswith("user_guide/pyazr.rst"):
        return text
    lines = text.split("\n")
    out, skipping = [], False
    for i, line in enumerate(lines):
        nxt = lines[i + 1] if i + 1 < len(lines) else ""
        if re.match(r"^[-=~^]{3,}$", nxt):
            skipping = line.strip() == "Earlier names"
        out.append("" if skipping else line)
    return "\n".join(out)


def program_options(binary):
    help_text = subprocess.run([binary, "--help"], capture_output=True, text=True).stdout
    return set(OPTION.findall(help_text))


def resolve(dotted):
    """True if pyazr.<...> resolves to a module or an attribute."""
    parts = dotted.split(".")
    for k in range(len(parts), 0, -1):
        try:
            obj = importlib.import_module(".".join(parts[:k]))
        except ImportError:
            continue
        for p in parts[k:]:
            if not hasattr(obj, p):
                return False
            obj = getattr(obj, p)
        return True
    return False


def has(cls, attr):
    """The class has the attribute, or its module sets it on an instance (self.<attr> = ...)."""
    if hasattr(cls, attr):
        return True
    import inspect
    try:
        source = inspect.getsource(inspect.getmodule(cls))
    except (OSError, TypeError):
        return False
    return re.search(r"\bself\.%s\b(?:\s*,[^=\n]*)?\s*=[^=]" % re.escape(attr), source) is not None


def main():
    binary = os.environ.get("AZURE2_BIN", "")
    options = program_options(binary) if binary else None
    import pyazr  # noqa: F401
    classes = {"Session": pyazr.Session, "Project": pyazr.Project, "Parameters": pyazr.Parameters,
               "Objective": pyazr.Objective, "Evaluation": pyazr.Evaluation}
    found = []
    checked = {"option": 0, "pyazr": 0}
    for dirpath, _, files in os.walk(PAGES):
        for name in sorted(files):
            if not name.endswith(".rst"):
                continue
            rel = os.path.relpath(os.path.join(dirpath, name), ROOT)
            with open(os.path.join(dirpath, name), encoding="utf-8") as f:
                text = skip_earlier_names(rel, f.read())
            for number, lit in literals(text):
                runs_program = re.match(r"^(\$ )?(\S*/)?AZURE2(\.exe)?\s", lit) is not None
                if options is not None and ((lit.startswith("--") and PROGRAM_PAGES.search(rel)) or runs_program) \
                        and not OTHER_PROGRAMS.search(lit):
                    for opt in OPTION.findall(lit):
                        checked["option"] += 1
                        if opt not in options and opt not in SHOWN_ON_PURPOSE:
                            found.append((rel, number, "option", opt))
                for dotted in PYAZR_NAME.findall(lit):
                    checked["pyazr"] += 1
                    if dotted.endswith((".py", ".rst", ".json")) or "/" in lit and dotted in lit.split("/")[-1]:
                        continue
                    if not resolve(dotted) and dotted not in SHOWN_ON_PURPOSE:
                        found.append((rel, number, "pyazr", dotted))
                for cls, attr in CLASS_ATTR.findall(lit):
                    checked["pyazr"] += 1
                    if attr not in FILE_SUFFIXES and not has(classes[cls], attr):
                        found.append((rel, number, "pyazr", f"{cls}.{attr}"))
                for var, attr in INSTANCE_ATTR.findall(lit):
                    if re.search(r"\bfor\s+%s\s+in\b" % var, lit):
                        continue  # a loop variable (a Parameter, a record), not a session or a project
                    checked["pyazr"] += 1
                    cls = classes["Session" if var == "s" else "Project"]
                    if attr not in FILE_SUFFIXES and not has(cls, attr):
                        found.append((rel, number, "pyazr", f"{var}.{attr} ({cls.__name__})"))
    for rel, number, kind, what in sorted(set(found)):
        print(f"{rel}:{number}: {kind}: {what}")
    print(f"{len(set(found))} name(s) the pages show that do not exist ({checked['option']} options and "
          f"{checked['pyazr']} pyazr names read)")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
