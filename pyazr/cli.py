"""The ``pyazr`` command line: ``python -m pyazr COMMAND ...`` (docs/upgrade/pyazr3.md section 2.8).

What pyazr can say and do about a project, as commands for a person at a shell and for an agent without a Python
session.  Each command is one row of :data:`COMMANDS` (its name, its arguments with their types, the pyazr function it
calls, its help line); the parser, every ``--help`` text, ``pyazr commands --json`` and the reference page
(:func:`reference_page`, docs/source/user_guide/pyazr_command_line.rst) are made from the rows.  A command is a call
of pyazr and a printer, nothing more: what a command would need beyond that is a gap of pyazr, named in the row's
``missing`` until pyazr has it.

Answers.  Without ``--json``: text a person reads on standard output.  With ``--json``: one JSON document on standard
output and nothing else there (what the engine prints meanwhile goes to standard error)::

    {"command": "validate", "ok": true, "exit": 0, "result": ...}
    {"command": "validate", "ok": false, "exit": 2, "result": ..., "error": {"class": ..., "what": ...,
                                                                           "allowed": ..., "fix": ...}}

``result`` is absent when the command refused its call; numbers that are not finite are ``null``.

Exit codes, those of the program's project commands (``AZURE2 --validate``, ``--describe``, ``--convert``;
docs/source/reference/command_line.rst, "Exit status"): 0 done; 1 the command ran and failed (a converter refusal, a
Coulomb table that cannot be read, a network error, an optional package missing); 2 a wrong call or a refused input
(an unknown command or option, a value that does not parse, a missing file, a project with errors, a command pyazr
cannot run yet); 3 a file that is not a format-3 project (``pyazr convert`` converts it).  An error is printed with
its three parts (what was wrong, what is allowed, the call that shows or fixes it): on standard error without
``--json``, in the document's ``error`` with it.

This module imports nothing of the engine at import time (the table is read by the consistency check without it);
each command imports what it calls.
"""

import argparse
import contextlib
import dataclasses
import difflib
import json
import math
import os
import sys
import textwrap
import traceback

PROG = "python -m pyazr"


class CallError(Exception):
    """A refusal with its exit code and three parts: what was wrong, what is allowed, the call that fixes it."""

    def __init__(self, exit_code, what, allowed, fix, kind="BadArgument"):
        super().__init__(what)
        self.exit_code, self.what, self.allowed, self.fix, self.kind = exit_code, what, allowed, fix, kind


# ---------------------------------------------------------------------------------------------------------------------
# Argument types: each a description (for --help, the reference page and `pyazr commands`) and a converter that turns
# the text given into the value the call takes, or raises CallError (exit 2, or 3 for an old project).
# ---------------------------------------------------------------------------------------------------------------------

def _existing_file(text, arg):
    if not os.path.isfile(text):
        folder = os.path.dirname(os.path.abspath(text))
        what = "is a directory, not a file" if os.path.isdir(text) else "does not exist"
        raise CallError(2, f"{arg.name} {text}: {what} (resolved: {os.path.abspath(text)})",
                        f"{arg.name}: an existing file", f"ls {folder}", "MissingFile")
    return text


def _project(text, arg):
    from . import errors, project
    try:
        project.require_format3(text)
    except errors.MissingFile as e:
        raise CallError(2, str(e), f"{arg.name}: an existing format-3 project file (.azr)",
                        f"ls {os.path.dirname(os.path.abspath(text))}", "MissingFile")
    except errors.FormatError as e:
        raise CallError(3, str(e), f"{arg.name}: a project in format 3 (this pyazr reads format 3 only)",
                        f"{PROG} convert {text} NEW.azr", "FormatError")
    return text


def _output(text, arg):
    folder = os.path.dirname(os.path.abspath(text))
    if not os.path.isdir(folder):
        raise CallError(2, f"{arg.name} {text}: the folder {folder} does not exist",
                        f"{arg.name}: a path in an existing folder", f"mkdir -p {folder}", "MissingFile")
    return text


def _number(kind, description):
    def convert(text, arg):
        try:
            return kind(text)
        except ValueError:
            raise CallError(2, f"{arg.name} {text!r}: not {description}", f"{arg.name}: {description}",
                            f"{PROG} {arg.command} --help")
    return convert


def _choice(text, arg):
    if text not in arg.choices:
        close = difflib.get_close_matches(text, arg.choices, n=1)
        raise CallError(2, f"{arg.name} {text!r}: not one of the choices" + (f" (did you mean {close[0]!r}?)" if close
                                                                               else ""),
                        f"{arg.name}: one of {', '.join(arg.choices)}", f"{PROG} {arg.command} --help")
    return text


def _json_file(text, arg):
    _existing_file(text, arg)
    try:
        with open(text, encoding="utf-8") as f:
            value = json.load(f)
    except ValueError as e:
        raise CallError(2, f"{arg.name} {text}: not JSON ({e})", f"{arg.name}: a JSON file holding one object",
                        f"python -m json.tool {text}")
    if not isinstance(value, dict):
        raise CallError(2, f"{arg.name} {text}: holds a {type(value).__name__}, not an object",
                        f"{arg.name}: a JSON file holding one object", f"python -m json.tool {text}")
    return value


def _project_or_new(text, arg):
    """An existing format-3 project, or a new ``.azr`` in an existing folder: the web UI starts an empty project
    there (written when saved)."""
    if not os.path.exists(text) and text.endswith(".azr") and os.path.isdir(os.path.dirname(os.path.abspath(text))):
        return text
    return _project(text, arg)


#: name -> (description, converter, the value's name in a usage line); a converter of None keeps the text.
ARG_TYPES = {
    "project": ("an existing format-3 project file (.azr)", _project, "FILE"),
    "project-or-new": ("an existing format-3 project file (.azr), or a new one in an existing folder",
                       _project_or_new, "FILE"),
    "file": ("an existing file", _existing_file, "FILE"),
    "output": ("a file to write, in an existing folder", _output, "FILE"),
    "number": ("a number", _number(float, "a number"), "X"),
    "int": ("a whole number", _number(int, "a whole number"), "N"),
    "text": ("a text", None, "TEXT"),
    "choice": ("one of the choices", _choice, "CHOICE"),
    "json-file": ("a JSON file holding one object", _json_file, "FILE"),
    "flag": ("an option without a value", None, ""),
}


@dataclasses.dataclass(frozen=True)
class Arg:
    """One argument of a command: ``name`` "FILE" (positional) or "--record" (an option), its type (a key of
    :data:`ARG_TYPES`), its help line; ``optional`` for a positional that may be left out, ``many`` for one or more
    values, ``default`` of an option, ``choices`` of a "choice"."""
    name: str
    type: str
    help: str
    optional: bool = False
    many: bool = False
    default: object = None
    choices: tuple = ()
    command: str = ""  # filled from the row: the command's name, for messages

    @property
    def dest(self):
        return self.name.lstrip("-").lower().replace("-", "_")

    @property
    def is_option(self):
        return self.name.startswith("--")

    @property
    def type_text(self):
        """What the argument takes, in words: its type, its choices, its default."""
        text = f"one of {', '.join(self.choices)}" if self.choices else ARG_TYPES[self.type][0]
        return text + (f"; default {self.default}" if self.default is not None else "")


@dataclasses.dataclass(frozen=True)
class Command:
    """One row of the table of commands.

    ``name``: the words of the command ("describe", "coulomb-table check"); ``args``: its :class:`Arg` rows;
    ``calls``: the pyazr function it calls (dotted, as a user imports it); ``help``: its help line; ``run``: a
    function ``(values, ctx) -> result`` that calls ``calls`` with the converted values (``ctx.json``: the --json
    mode); ``text``: the printer of the result for a person (None: the generic one); ``status``: a function
    ``result -> None or CallError`` for a result that is a failure (a project with errors); ``refused_exit``: the exit
    code of a ``ValueError`` the call raises (2, a refused input; 1 where the program says 1); ``fix``: the call that
    helps when the input is refused; ``missing``: (what, allowed, fix) when pyazr cannot run the command yet (it then
    exits 2 with these three parts)."""
    name: str
    args: tuple
    calls: str
    help: str
    run: object = None
    text: object = None
    status: object = None
    refused_exit: int = 2
    fix: str = ""
    missing: tuple = ()

    def arguments(self):
        return tuple(dataclasses.replace(a, command=self.name) for a in self.args)


# ---------------------------------------------------------------------------------------------------------------------
# The calls (each a call of pyazr; the result is plain data or the engine's text) and their printers.
# ---------------------------------------------------------------------------------------------------------------------

def _run_commands(v, ctx):
    return [_row_record(c) for c in COMMANDS]


def _text_commands(result, v):
    return "\n".join(f"{_usage(c)}\n    {c.help}" + (f" Calls {c.calls}." if c.calls else "") for c in COMMANDS)


def _run_describe(v, ctx):
    from . import describe
    if v.file is None:
        return describe()
    return describe(v.file, as_json=ctx.json, record=v.record)


def _run_validate(v, ctx):
    from . import validate
    return validate(v.file)


def _validate_status(result, v):
    errors = [d for d in result if d.get("severity") == "error"]
    if errors:
        return CallError(2, f"{v.file}: {len(errors)} error(s)", "a project without errors",
                         "the 'fix' of each diagnostic; then " + f"{PROG} describe {v.file}", "FormatError")
    return None


def _text_validate(result, v):
    lines = [f"{d.get('file')}:{d.get('line')}: "
             + ".".join(x for x in (d.get("record") or "", d.get("field") or "") if x)
             + f": {d.get('code')}: {d.get('message')}" for d in result]
    ok = not any(d.get("severity") == "error" for d in result)
    return "\n".join(lines + [f"{v.file}: {len(result)} diagnostic(s), {'valid' if ok else 'not valid'}."])


def _run_convert(v, ctx):
    from . import convert_azr
    return convert_azr(v.old, v.new, overwrite=v.overwrite, report=v.report, answers=v.answers,
                       exact_grids=v.exact_grids, strict=v.strict, allow_questions=v.allow_questions)


def _text_convert(result, v):
    lines = [f"{result.get('input')} -> {result.get('output')}: converted"]
    lines += [f"  changed: {c}" for c in result.get("format2_changes", [])]
    lines += [f"  note: {n.get('record') or n.get('construct')}: {n.get('text')}" for n in result.get("notes", [])]
    lines += [f"  wrote: {f}" for f in result.get("side_files", [])]
    return "\n".join(lines)


def _run_schema(v, ctx):
    from . import schema
    table = schema()
    if v.kind is None:
        return table
    for part in ("records", "quantities"):
        if v.kind in table[part]:
            return {"table": part, "name": v.kind, "entry": table[part][v.kind]}
    names = sorted(table["records"]) + sorted(table["quantities"])
    raise CallError(2, f"KIND {v.kind!r}: neither a record kind nor a quantity of the schema"
                    + _close(v.kind, names), "KIND: one of " + ", ".join(names), f"{PROG} schema --json")


def _run_kinds(v, ctx):
    from . import kinds
    table = kinds()
    if v.table is None:
        return table
    if v.table not in table:
        raise CallError(2, f"TABLE {v.table!r}: not a table of pyazr.kinds()" + _close(v.table, table),
                        "TABLE: one of " + ", ".join(table), f"{PROG} kinds")
    return table[v.table]


def _text_kinds(result, v):
    if v.table is None:
        return "\n".join(f"{k}: {', '.join(t) if isinstance(t, dict) else t}" for k, t in result.items())
    if not isinstance(result, dict):
        return _text(result)
    return "\n".join(f"{name}: {_first(entry)}" for name, entry in result.items())


def _first(entry):
    if isinstance(entry, dict):
        for key in ("title", "doc", "signature", "unit"):
            if entry.get(key):
                return str(entry[key]).splitlines()[0]
        return ""
    return str(entry)


def _run_stopping(v, ctx):
    from . import stopping
    material = stopping.Material(v.material, phase=v.phase, compound_correction=v.compound_correction)
    values = stopping.stopping(material, v.ion, v.energy)
    return {"material": v.material, "phase": material.phase_used, "ion": v.ion, "energy_MeV": v.energy,
            **values, "unit": "eV/(1e15 atoms/cm^2) per atom of the material; rel_uncertainty relative, 1 sigma",
            "citation": stopping.CITATION}


def _text_stopping(result, v):
    columns = ["energy_MeV", "electronic", "nuclear", "total", "rel_uncertainty"]
    rows = [[f"{result[c][i]:.6g}" for c in columns] for i in range(len(result["energy_MeV"]))]
    head = f"{result['ion']} in {result['material']} ({result['phase']}); {result['unit']}"
    return "\n".join([head, _table(columns, rows), f"Cite: {result['citation']}"])


def _run_coulomb_check(v, ctx):
    from .described import coulomb_table
    return coulomb_table(coulomb_table=v.file or "")


def _coulomb_status(result, v):
    if result.get("loaded"):
        return None
    return CallError(1, f"the Coulomb table {result.get('path')} cannot be read (the engine uses its exact routines)",
                     "a table built by pyazr (or AZURE2 --build-coulomb-table[=FILE]); FILE, else "
                     "$AZURE2_COULOMB_TABLE, else ~/.azure2/coulomb_table.bin",
                     f"{PROG} coulomb-table build {result.get('path')}", "CoulombTableMissing")


def _run_coulomb_where(v, ctx):
    from . import coulomb_table
    return coulomb_table.where(v.file or "")


def _run_coulomb_fetch(v, ctx):
    from . import coulomb_table
    return coulomb_table.fetch(v.file or "", url=v.url, sha256=v.sha256, overwrite=v.overwrite)


def _run_coulomb_build(v, ctx):
    from . import coulomb_table
    return coulomb_table.build(v.file or "", threads=v.threads, l_max=v.l_max, per_decade=v.per_decade)


def _run_ui(v, ctx):
    try:
        from .ui.server import main
    except ImportError as e:
        raise CallError(1, f"the web UI cannot start here: {e}", "pyazr with its extra 'web' (fastapi, uvicorn)",
                        "pip install pyazr[web]  (not done by this command)", "ImportError")
    with ctx.first_line(lambda line: {"url": line.strip()}):
        try:
            main(v.file, v.port, v.threads, v.token)
        except KeyboardInterrupt:  # the end the command names: "until interrupted"
            pass
    return None


def _run_import_exfor(v, ctx):
    from . import fetch_exfor
    return fetch_exfor(v.id, plus=v.plus, source=v.source)


def _text_exfor(result, v):
    head = f"{result['dataset_id']} {result['reaction']} ({result['author']} {result['year']})"
    names = [c["name"] for c in result["columns"] if c["name"] in result["arrays"]]
    columns = [f"{c['name']} ({c['unit']})" for c in result["columns"] if c["name"] in result["arrays"]]
    n = len(result["arrays"][names[0]]) if names else 0
    rows = [[_cell(result["arrays"][name][i]) for name in names] for i in range(n)]
    return "\n".join([head, _table(columns, rows)])


def _run_import_ensdf(v, ctx):
    from . import fetch_levels
    return fetch_levels(v.nuclide, source=v.source)


def _text_ensdf(result, v):
    columns = ["index", "energy", "energy_unc", "jp", "half_life", "width_kev"]
    return "\n".join([f"{v.nuclide}: ENSDF adopted levels (energies in keV)",
                      _table(columns, [[_cell(level[c]) for c in columns] for level in result])])


def _run_endf_summary(v, ctx):
    from . import endf
    summary = endf.read(v.file, mat=v.mat).summary()
    return summary.to_dict() if ctx.json else str(summary)


def _run_endf_curves(v, ctx):
    from . import Session, endf
    from .errors import FormatRefused
    material = endf.read(v.file, mat=v.mat)
    out = []
    with Session(v.project) as s:
        for d in s.evaluate(s.parameters.x0).data:
            if v.record and d.name != v.record:
                continue
            try:
                c = material.curves(d, session=s)
                row = {"data set": d.name, "quantity": c.quantity, "unit": c.unit, "source": c.source,
                       "covered": int(c.coverage.sum()), "points": len(c.values), "notes": c.notes}
                if ctx.json:
                    row.update(energies=c.energies, angles=c.angles, values=c.values, model=d.model)
            except FormatRefused as e:
                row = {"data set": d.name, "quantity": d.quantity, "unit": "", "source": "", "covered": 0,
                       "points": d.n, "notes": [f"refused: {e.what}; allowed: {e.allowed}"]}
            out.append(row)
    if v.record and not out:
        raise CallError(2, f"--record {v.record}: no data set of that name in {v.project}",
                        "--record: the name of a data set of the project", f"{PROG} describe {v.project}")
    return out


def _text_endf_curves(result, v):
    rows = [[r["data set"], r["quantity"], r["source"], f"{r['covered']}/{r['points']}", "; ".join(r["notes"])]
            for r in result]
    return _table(["data set", "quantity", "from", "covered", "notes"], rows)


def _run_endf_to_project(v, ctx):
    from . import endf
    material = endf.read(v.file, mat=v.mat)
    sep = "ame2020" if v.separation_energy is None else v.separation_energy
    _, report = endf.to_project(material.resonances(isotope=v.isotope, energy_range=v.range), v.new,
                                separation_energy=sep, material=material, overwrite=v.overwrite)
    if not ctx.json:
        return str(report)
    return {"project": report.path, "separation_energy_MeV": report.separation_MeV,
            "separation_source": report.separation_source,
            "mass_checks": [{"nuclide": n, "file_u": f, "ame2020_u": a} for n, f, a in report.mass_checks],
            "pairs": list(report.pairs), "levels": [dict(vars(lv)) for lv in report.levels],
            "not_written": list(report.skipped), "notes": list(report.notes)}


def _run_endf_write(v, ctx):
    from . import Session, endf
    header = v.header  # the JSON object, read by the argument's type
    rng = (v.low, v.high) if v.low is not None and v.high is not None else None
    with Session(v.project) as s:
        written = endf.write(s, v.new, header=header, energy_range_MeV=rng, tolerance=v.tolerance,
                             number_form=v.number_form, overwrite=v.overwrite, elastic_form=v.elastic_form or endf.writer.ELASTIC_FORM,
                             covariance="hesse" if v.covariance else None)
    if not ctx.json:
        return str(written)
    return {"file": written.path, "sections": [{"mf": a, "mt": b, "lines": c} for a, b, c in written.sections],
            "grids": {str(k): {"points": g[0], "low_eV": g[1], "high_eV": g[2]} for k, g in written.grids.items()},
            "mf2_check": {str(k): v for k, v in written.mf2_check.items()}, "losses": list(written.losses),
            "not_written": list(written.not_written)}


def _run_endf_check(v, ctx):
    from . import endf
    result = endf.check(v.file, codes=[v.program] if v.program else None, timeout=v.timeout,
                        resonances=v.resonances)
    return result.to_dict() if ctx.json else str(result)


def _run_endf_standardize(v, ctx):
    from . import endf
    candidate = endf.standardize(v.file, program=v.program, timeout=v.timeout, mf2=v.mf2)
    written = candidate.accept(to=v.write, overwrite=v.overwrite) if v.write else None
    if ctx.json:
        return {"original": v.file, "candidate": candidate.path, "program": candidate.program, "written": written,
                "changes": [c for c in candidate.diff() if c.kind != "unchanged"]}
    return candidate.summary() + (f"\nwritten: {written}" if written else
                                  "\nnothing written (--write NEW keeps the candidate)")


_COULOMB_FILE = Arg("FILE", "text", "the table file (default: $AZURE2_COULOMB_TABLE, else ~/.azure2/coulomb_table.bin)",
                    optional=True)

#: The table of commands.  A new command is one row here (with its call, a few lines above) and its case in
#: tests/pyazr/cli_test.py; the parser, --help, `pyazr commands` and the reference page follow from the row.
COMMANDS = (
    Command("commands", (), "pyazr.cli.COMMANDS",
            "List the commands: name, arguments with their types, the pyazr function each calls, help.",
            run=_run_commands, text=_text_commands),
    Command("describe", (Arg("FILE", "project", "the project; without it, the package and its engine", optional=True),
                         Arg("--record", "text", "only the data set, extrapolation, datum or potential of this name")),
            "pyazr.describe",
            "What the engine makes of a project (as AZURE2 --describe), or the package and its engine.",
            run=_run_describe, fix="python -m pyazr validate FILE lists every diagnostic"),
    Command("validate", (Arg("FILE", "project", "the project"),), "pyazr.validate",
            "Every problem of a project, as AZURE2 --validate (exit 2 when one is an error).",
            run=_run_validate, text=_text_validate, status=_validate_status),
    Command("convert", (Arg("OLD", "file", "the project in format 1 or 2"),
                        Arg("NEW", "output", "the project in format 3 to write"),
                        Arg("--overwrite", "flag", "replace NEW, the report and the side files if they exist"),
                        Arg("--report", "output", "write the converter's report (JSON) to this file"),
                        Arg("--answers", "json-file", "answers to the converter's questions, "
                                                      "{\"<record>.<field>\": value}"),
                        Arg("--exact-grids", "flag", "keep the old energy grids exactly"),
                        Arg("--strict", "flag", "refuse instead of noting an implicit meaning"),
                        Arg("--allow-questions", "flag", "write '?' for an unanswered question instead of refusing")),
            "pyazr.convert_azr",
            "Write an old project (format 1 or 2) in format 3, as AZURE2 --convert; OLD is kept.",
            run=_run_convert, text=_text_convert, refused_exit=1,
            fix="the converter's report (--report FILE) names the construct it refused"),
    Command("schema", (Arg("KIND", "text", "a record kind or a quantity: only its entry", optional=True),),
            "pyazr.schema", "The table of project format 3 (as AZURE2 --schema), or one record kind or quantity of it.",
            run=_run_schema),
    Command("kinds", (Arg("TABLE", "text", "records, data, parameters, options, physics, units, dimensions or errors",
                          optional=True),),
            "pyazr.kinds", "The kinds of thing the engine and pyazr know (record, data, parameter kinds, options, "
                           "physics functions, errors).", run=_run_kinds, text=_text_kinds),
    Command("stopping", (Arg("MATERIAL", "text", "the target's chemical formula (Ta2O5, C, CH4)"),
                         Arg("ION", "choice", "the projectile", choices=("p", "d", "t", "3He", "a")),
                         Arg("ENERGY", "number", "lab kinetic energies of the projectile, MeV", many=True),
                         Arg("--phase", "choice", "solid or gas (default: SRIM's for the material)",
                             choices=("solid", "gas")),
                         Arg("--compound-correction", "number", "SRIM's factor of Bragg's rule", default=1.0)),
            "pyazr.stopping.stopping",
            "Stopping cross sections of a material for an ion (SRIM-2013 as ported from LunaSRIM; cite SRIM).",
            run=_run_stopping, text=_text_stopping, fix="python -m pyazr stopping --help"),
    Command("coulomb-table check", (_COULOMB_FILE,), "pyazr.described.coulomb_table",
            "Load the Coulomb table as a session does and give its path and content identity (exit 1: not readable).",
            run=_run_coulomb_check, status=_coulomb_status),
    Command("coulomb-table where", (_COULOMB_FILE,), "pyazr.coulomb_table.where",
            "Where the Coulomb table is looked for and whether a file is there (without loading it).",
            run=_run_coulomb_where),
    Command("coulomb-table fetch", (_COULOMB_FILE,
                                    Arg("--url", "text", "the copy to download (required; no table is published)"),
                                    Arg("--sha256", "text", "the SHA-256 the copy must have (required, with --url)"),
                                    Arg("--overwrite", "flag", "replace an existing file")),
            "pyazr.coulomb_table.fetch",
            "Copy a Coulomb table from a URL you name and check its SHA-256 (exit 2 without --url and --sha256: "
            "build the table here).",
            run=_run_coulomb_fetch, fix="python -m pyazr coulomb-table build"),
    Command("coulomb-table build", (_COULOMB_FILE,
                                    Arg("--threads", "int", "OpenMP threads of the build (default: the engine's)"),
                                    Arg("--l-max", "int", "the largest l of the table (default: the program's, 10)"),
                                    Arg("--per-decade", "number", "nodes per decade of eta and rho (default 100)")),
            "pyazr.coulomb_table.build",
            "Build the Coulomb table here from the exact routines (about 5 min on 4 threads, 88 MB).",
            run=_run_coulomb_build),
    Command("ui", (Arg("FILE", "project-or-new", "the project to serve; a new .azr starts an empty project"),
                   Arg("--port", "int", "the TCP port (0: a free one above 20000)", default=0),
                   Arg("--threads", "int", "the compute budget in threads", default=1),
                   Arg("--token", "text", "a fixed access token (default: a random one)")),
            "pyazr.ui.server.main",
            "Serve a project to the web UI on 127.0.0.1 until interrupted (--json: the URL as soon as it listens).",
            run=_run_ui),
    Command("import exfor", (Arg("ID", "text", "the EXFOR data set (O2599004; a pointer as O2599004.1)"),
                             Arg("--source", "choice", "the Web-API, the local EXFOR master files "
                                                       "($AZURE2_EXFOR_MASTER), or the first that answers",
                                 choices=("auto", "web", "local"), default="auto"),
                             Arg("--plus", "int", "0 original units, 1 computational units, 2 the universal grid",
                                 default=0)),
            "pyazr.fetch_exfor", "One EXFOR data set: reaction, author, year, columns with units, values.",
            run=_run_import_exfor, text=_text_exfor, fix="python -m pyazr import exfor --help"),
    Command("import ensdf", (Arg("NUCLIDE", "text", "the nuclide, as 14n or 20Ne"),
                             Arg("--source", "choice", "IAEA LiveChart, the NNDC ENSDF API, or the first that answers",
                                 choices=("auto", "livechart", "nndc"), default="auto")),
            "pyazr.fetch_levels", "The ENSDF adopted levels of a nuclide (energies in keV, J^pi, half-lives).",
            run=_run_import_ensdf, text=_text_ensdf),
    Command("endf summary", (Arg("FILE", "file", "an ENDF-6 file"),
                             Arg("--mat", "int", "the material (MAT) of a file that holds several")),
            "pyazr.endf.read",
            "What an ENDF-6 file holds: material, projectile, library, sections in words, the header text.",
            run=_run_endf_summary, fix="python -m pyazr endf summary --help"),
    Command("endf curves", (Arg("FILE", "file", "an ENDF-6 file (the evaluation)"),
                            Arg("PROJECT", "project", "the project whose data sets are compared"),
                            Arg("--record", "text", "only the data set of this name"),
                            Arg("--mat", "int", "the material (MAT) of a file that holds several")),
            "pyazr.endf.Material.curves",
            "The evaluation at each data set's points in its quantity and unit (--json: the values beside the "
            "model at the project's start values).",
            run=_run_endf_curves, text=_text_endf_curves, fix="python -m pyazr endf summary FILE lists what it holds"),
    Command("endf to-project", (Arg("FILE", "file", "an ENDF-6 file with MF=2/MT=151 in LRF=7 (R-Matrix Limited)"),
                                Arg("NEW", "output", "the format-3 project to write (Brune basis, every level fixed)"),
                                Arg("--separation-energy", "number", "S of the incident pair, MeV (default: AME2020)"),
                                Arg("--isotope", "int", "the isotope of MF=2 (0-based)", default=0),
                                Arg("--range", "int", "the energy range of the isotope (0-based)", default=0),
                                Arg("--overwrite", "flag", "replace NEW if it exists"),
                                Arg("--mat", "int", "the material (MAT) of a file that holds several")),
            "pyazr.endf.to_project",
            "An evaluation's resonance parameters as a new project: pairs, levels and channels, the conversions and "
            "what was not written.", run=_run_endf_to_project,
            fix="python -m pyazr endf summary FILE shows its MF=2"),
    Command("endf write", (Arg("PROJECT", "project", "the project whose model is written (at its start values)"),
                           Arg("NEW", "output", "the ENDF-6 file to write"),
                           Arg("--header", "json-file", "the MF=1 fields: {MAT, NLIB, NVER, LREL, NMOD, ALAB, EDATE, "
                                                       "AUTH, REF, DDATE, RDATE, LIBRARY, TEXT}"),
                           Arg("--low", "number", "lowest lab energy of the projectile, MeV (default: the data's)"),
                           Arg("--high", "number", "highest lab energy of the projectile, MeV (default: the data's)"),
                           Arg("--tolerance", "number", "relative tolerance of the engine's adaptive grid",
                               default=1e-3),
                           Arg("--number-form", "choice", "how numbers are written", default="widest",
                               choices=("widest", "classic", "e")),
                           Arg("--elastic-form", "choice", "charged-particle elastic scattering as the nuclear "
                               "amplitude expansion (LTP=1) or tabulated (LTP=12) (default: pyazr.endf.writer.ELASTIC_FORM)",
                               choices=("amplitudes", "tabulated")),
                           Arg("--covariance", "flag", "write MF=32: the fit's covariance (HESSE at x) carried through the transformation"),
                           Arg("--overwrite", "flag", "replace NEW if it exists")),
            "pyazr.endf.write",
            "The model as an ENDF-6 material: MF=1, MF=2 (LRF=7, formal basis), MF=3 and MF=6 from the engine "
            "(LRP=2); what the file does not carry is listed.", run=_run_endf_write,
            fix="python -m pyazr endf write --help"),
    Command("endf check", (Arg("FILE", "file", "an ENDF-6 file"),
                           Arg("--program", "choice", "the one checking code to run (default: CHECKR, FIZCON and "
                               "PSYCHE, those found)", choices=("checkr", "fizcon", "psyche", "inter")),
                           Arg("--timeout", "number", "seconds each code may take", default=300.0),
                           Arg("--resonances", "flag", "also check the resonance parameters: MF=2 against MF=3 by "
                               "the engine, the masses against AME2020, the text against the records")),
            "pyazr.endf.check",
            "The module's own format check, then the ENDF checking codes found ($AZURE2_ENDF_UTILITIES), their "
            "reports read into findings.", run=_run_endf_check,
            fix="python -m pyazr endf check --help; pyazr.endf.checkers() lists the codes found"),
    Command("endf standardize", (Arg("FILE", "file", "an ENDF-6 file (never changed)"),
                                 Arg("--program", "choice", "the code that rewrites a copy", default="stanef",
                                     choices=("stanef", "fixup", "dictin", "endf2c")),
                                 Arg("--mf2", "flag", "STANEF: rewrite the numbers of MF=2/MT=151 too"),
                                 Arg("--write", "output", "keep the candidate as this file"),
                                 Arg("--overwrite", "flag", "with --write: replace an existing file"),
                                 Arg("--timeout", "number", "seconds the code may take", default=300.0)),
            "pyazr.endf.standardize",
            "Run a code that rewrites a copy of the file and show what changed per section (numbers, digits lost); "
            "nothing is written without --write.", run=_run_endf_standardize,
            fix="python -m pyazr endf standardize --help; pyazr.endf.checkers() lists the codes found"),
)


# ---------------------------------------------------------------------------------------------------------------------
# The printer: plain data for --json, text for a person.
# ---------------------------------------------------------------------------------------------------------------------

def plain(value):
    """``value`` as JSON data: arrays and tuples as lists, dataclasses as dicts, numbers that are not finite as None."""
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {f.name: plain(getattr(value, f.name)) for f in dataclasses.fields(value)}
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if hasattr(value, "tolist"):  # numpy arrays and scalars
        return plain(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _cell(x):
    if x is None:
        return ""
    return f"{x:.6g}" if isinstance(x, float) else str(x)


def _table(columns, rows):
    widths = [max([len(c)] + [len(r[i]) for r in rows]) for i, c in enumerate(columns)]
    return "\n".join("  ".join(s.ljust(w) for s, w in zip(r, widths)).rstrip() for r in [columns] + rows)


def _text(value):
    """The generic printer: a text as it is; a dict as 'key: value' lines; a list of dicts as a table, of texts or
    numbers one per line."""
    if isinstance(value, str):
        return value.rstrip("\n")
    if isinstance(value, list) and all(isinstance(r, (str, int, float)) for r in value):
        return "\n".join(_cell(r) for r in value)
    if isinstance(value, list) and value and all(isinstance(r, dict) for r in value):
        columns = list(value[0])
        return _table(columns, [[_cell(r.get(c)) if not isinstance(r.get(c), (dict, list))
                                 else json.dumps(r.get(c)) for c in columns] for r in value])
    if isinstance(value, dict):
        return "\n".join(f"{k}: {json.dumps(v) if isinstance(v, (dict, list)) else _cell(v)}" for k, v in value.items())
    return json.dumps(value, indent=1)


def _close(name, names):
    close = difflib.get_close_matches(str(name), [str(n) for n in names], n=2)
    return f" (did you mean {' or '.join(repr(c) for c in close)}?)" if close else ""


def _row_record(c):
    """A row as plain data (``pyazr commands --json``)."""
    return {"name": c.name, "help": c.help, "calls": c.calls or None,
            "missing": dict(zip(("what", "allowed", "fix"), c.missing)) if c.missing else None,
            "arguments": [{"name": a.name, "type": a.type, "type_doc": a.type_text, "help": a.help,
                           "option": a.is_option, "required": not (a.optional or a.is_option),
                           "many": a.many, "default": a.default, "choices": list(a.choices) or None}
                          for a in c.args]}


def _usage(c):
    words = []
    for a in c.args:
        if a.is_option:
            words.append(f"[{a.name}]" if a.type == "flag" else f"[{a.name} {ARG_TYPES[a.type][2]}]")
        else:
            words.append(f"[{a.name}]" if a.optional else f"{a.name} [{a.name} ...]" if a.many else a.name)
    return " ".join([PROG, c.name] + words + ["[--json]"])


def reference_page():
    """The reference page of the command line (docs/source/user_guide/pyazr_command_line.rst), made from the table."""
    title = "The pyazr command line"
    out = [title, "=" * len(title), "",
           "``python -m pyazr COMMAND ...`` runs what pyazr can say and do about a project from a shell: each command",
           "is a call of one pyazr function.  ``--json`` gives one JSON document on standard output",
           "(``{\"command\", \"ok\", \"exit\", \"result\"}``, and ``\"error\"`` with its parts ``what``, ``allowed``,",
           "``fix`` when it fails); without it the answer is text.  Exit codes: 0 done, 1 the command ran and failed,",
           "2 a wrong call or a refused input, 3 a file that is not a format-3 project (as the program's",
           ":doc:`../reference/command_line`).  ``python -m pyazr commands --json`` gives this table as data.",
           "",
           "This page is generated from the table of commands (``pyazr/cli.py``, ``reference_page()``);",
           "``tests/consistency/repo_consistency_test.py`` checks that it is.", ""]
    for c in COMMANDS:
        head = f"``{c.name}``"
        out += [head, "-" * len(head), "", c.help, "", ".. code-block:: text", "", "   " + _usage(c), ""]
        if c.calls:
            out += [f"Calls ``{c.calls}``.", ""]
        if c.missing:
            out += textwrap.wrap(f"Not available yet: {c.missing[0]}.  Instead: {c.missing[1]}.", 116) + [""]
        for a in c.args:
            out += textwrap.wrap(f"- ``{a.name}``: {a.help} ({a.type_text}).", 116, subsequent_indent="  ")
        out += ["- ``--json``: one JSON document on standard output.", ""]
    return "\n".join(out)


REFERENCE_PAGE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "source", "user_guide",
                              "pyazr_command_line.rst")


def write_reference_page(path=REFERENCE_PAGE):
    """Write :func:`reference_page` to ``path`` (the repository's page by default)."""
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(reference_page() + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# Parser and main
# ---------------------------------------------------------------------------------------------------------------------

class _Parser(argparse.ArgumentParser):
    """argparse with its errors raised as CallError (exit 2, three parts) instead of printed."""

    def error(self, message):
        raise CallError(2, message, self.format_usage().strip(), f"{self.prog} --help", "UsageError")


def build_parser():
    """The parser, made from :data:`COMMANDS`."""
    top = _Parser(prog=PROG, description="The pyazr command line: what pyazr can say and do about a project. "
                  f"'{PROG} commands' lists the commands; every command takes --json.")
    groups = {(): top.add_subparsers(dest="word0", metavar="COMMAND", parser_class=_Parser)}
    parsers = {}
    for c in COMMANDS:
        words = c.name.split()
        for depth in range(len(words)):
            key = tuple(words[:depth + 1])
            if key in parsers:
                continue
            last = depth == len(words) - 1
            p = groups[key[:-1]].add_parser(words[depth], help=c.help if last else f"{words[0]} commands",
                                            description=c.help if last else None)
            parsers[key] = p
            if not last:
                groups[key] = p.add_subparsers(dest=f"word{depth + 1}", metavar="COMMAND", parser_class=_Parser)
        p = parsers[tuple(words)]
        p.set_defaults(command=c)
        for a in c.args:
            text = f"{a.help} ({a.type_text})"
            if a.is_option and a.type == "flag":
                p.add_argument(a.name, action="store_true", help=text)
            elif a.is_option:
                p.add_argument(a.name, default=a.default, metavar=ARG_TYPES[a.type][2], help=text)
            else:
                p.add_argument(a.dest, metavar=a.name, nargs="+" if a.many else "?" if a.optional else None,
                               default=a.default, help=text)
        p.add_argument("--json", action="store_true", help="one JSON document on standard output, nothing else")
    return top


def _convert(c, ns):
    for a in c.arguments():
        convert = ARG_TYPES[a.type][1]
        value = getattr(ns, a.dest)
        if convert is None or not isinstance(value, (str, list)):  # a default given as its value, or None
            continue
        setattr(ns, a.dest, [convert(x, a) for x in value] if a.many else convert(value, a))


def _failure(e, c):
    """An exception of a call as a CallError with its exit code and three parts."""
    if isinstance(e, CallError):
        return e
    name = type(e).__name__
    if getattr(e, "what", None) and getattr(e, "allowed", None) and getattr(e, "fix", None):  # a pyazr.Error's parts
        return CallError(c.refused_exit if isinstance(e, ValueError) else 1, str(e), e.allowed, e.fix, name)
    allowed = "; ".join(f"{a.name}: {a.type_text}" for a in c.args) or "the command without arguments"
    fix = c.fix or f"{PROG} {c.name} --help"
    if isinstance(e, FileExistsError):
        return CallError(2, str(e), "a path that does not exist yet", f"{PROG} {c.name} ... --overwrite", name)
    if isinstance(e, FileNotFoundError):
        return CallError(2, str(e), allowed, fix, name)
    if isinstance(e, KeyError) and e.args and isinstance(e.args[0], str):
        return CallError(c.refused_exit, e.args[0], allowed, fix, name)  # str(KeyError) is the repr of its text
    if isinstance(e, (ValueError, KeyError)):
        return CallError(c.refused_exit, str(e), allowed, fix, name)
    if isinstance(e, OSError):
        return CallError(1, str(e), "a reachable service or a readable file", fix, name)
    return CallError(1, f"{name}: {e}", allowed, fix, name)


class _Context:
    """What a call may need of the printer: the mode, and the first line of a call that prints and then blocks."""

    def __init__(self, as_json, out):
        self.json, self.out, self.emitted, self.command = as_json, out, False, None

    def emit(self, result):
        self.out.write(json.dumps(_document(self.command, 0, result), allow_nan=False) + "\n")
        self.out.flush()
        self.emitted = True

    @contextlib.contextmanager
    def first_line(self, make):
        """In --json mode, the call's first printed line becomes the document's result (``make(line)``), sent at
        once; whatever it prints later goes to standard error."""
        if not self.json:
            yield
            return
        ctx, saved = self, sys.stdout

        class Catcher:
            buffer = ""

            def write(self, text):
                if ctx.emitted:
                    return sys.stderr.write(text)
                self.buffer += text
                if "\n" in self.buffer:
                    line, rest = self.buffer.split("\n", 1)
                    ctx.emit(make(line))
                    sys.stderr.write(rest)
                return len(text)

            def flush(self):
                sys.stderr.flush()

            def isatty(self):  # asked by logging set-ups (uvicorn's formatter) of the stream they write to
                return False

        sys.stdout = Catcher()
        try:
            yield
        finally:
            sys.stdout = saved


def _document(command, code, result=None, error=None):
    doc = {"command": command.name if command else None, "ok": code == 0, "exit": code}
    if result is not None or error is None:
        doc["result"] = plain(result)
    if error is not None:
        doc["error"] = {"class": error.kind, "what": error.what, "allowed": error.allowed, "fix": error.fix}
    return doc


@contextlib.contextmanager
def _stdout_to_stderr():
    """Standard output (Python's and the engine's, file descriptor 1) sent to standard error; yields the real one."""
    sys.stdout.flush()
    saved = os.dup(1)
    os.dup2(2, 1)
    real, sys.stdout = sys.stdout, sys.stderr
    out = os.fdopen(os.dup(saved), "w", encoding="utf-8")
    try:
        yield out
    finally:
        from . import engine_text
        engine_text.flush_c_stdout()  # the engine's buffered output, before descriptor 1 is restored
        sys.stdout = real
        os.dup2(saved, 1)
        os.close(saved)
        out.flush()


def main(argv=None):
    """Run one command; returns the exit code (0, 1, 2 or 3; see the module's docstring)."""
    argv = list(sys.argv[1:] if argv is None else argv)
    as_json = "--json" in argv
    command, result, error = None, None, None
    with (_stdout_to_stderr() if as_json else contextlib.nullcontext(sys.stdout)) as out:
        ctx = _Context(as_json, out)
        try:
            ns = build_parser().parse_args(argv)
            command = ctx.command = getattr(ns, "command", None)
            if command is None:
                given = [w for w in (getattr(ns, f"word{i}", None) for i in range(3)) if w]
                names = [c.name for c in COMMANDS if c.name.split()[:len(given)] == given]
                given = " ".join(given)
                raise CallError(2, f"'{given}' needs a command after it" if given else "no command given",
                                "one of: " + ", ".join(names), f"{PROG} commands", "UsageError")
            if command.missing:
                raise CallError(2, *command.missing, kind="NotAvailable")
            _convert(command, ns)
            result = command.run(ns, ctx)
            error = command.status(result, ns) if command.status else None
        except Exception as e:  # every refusal becomes an exit code and its three parts
            if command or not isinstance(e, CallError):
                error = _failure(e, command or COMMANDS[0])
            else:
                error = e
            if not isinstance(e, CallError) and not isinstance(e, (ValueError, KeyError, OSError)):
                traceback.print_exc(file=sys.stderr)
        code = error.exit_code if error else 0
        if ctx.emitted:
            return code
        if as_json:
            out.write(json.dumps(_document(command, code, result, error), allow_nan=False) + "\n")
        else:
            if result is not None:
                printer = command.text if command.text else (lambda r, v: _text(r))
                text = printer(plain(result), ns)
                if text:
                    print(text)
            if error is not None:
                where = f"{PROG} {command.name}" if command else PROG
                print(f"{where}: error: {error.what}\n  allowed: {error.allowed}\n  fix: {error.fix}", file=sys.stderr)
        return code
