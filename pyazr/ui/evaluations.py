"""The page's ENDF screens on :mod:`pyazr.endf`: the ENDF-6 files of the project folder, what a file holds
(``Material.summary()``), the file's values at every data set's points (``Material.curves``, through the session
worker: :meth:`pyazr.ui.engine.Engine.endf_curves`), the checking codes found (``pyazr.endf.checkers()``), and the
two jobs a user starts with a button: the check of a file (``pyazr.endf.check``) and a rewrite by a code
(``pyazr.endf.standardize``) whose candidate replaces nothing until it is kept (``Candidate.accept``).

The format work and every number are pyazr's.  Two ratios are formed here, named in the answer: the file's value
over the datum times its normalization (``ratio_to_data``, with ``ratio_to_data_error``, the ratio times the
datum's relative error) and over the model (``ratio_to_model``), at each point where the file has a value.  The
files a user names are checked to lie inside the server's root by the server before any of these functions sees
them (the state rule, docs/upgrade/web_ui.md section 5).
"""
import dataclasses
import os

import numpy as np

#: How deep below the project's directory the listing of ENDF files looks (0: the directory itself).
SEARCH_DEPTH = 2
#: The most files the listing reads the first lines of (a project folder with a whole library in it is cut there).
MAX_FILES = 2000
#: What the page says of how the checking codes are found (the order of ``pyazr.endf.checkers()``).
HOW_FOUND = ("The checking codes are found in the directories named by the environment variables "
             "AZURE2_ENDF_UTILITIES and AZURE2_PREPRO of the server, then on its PATH (pyazr.endf.checkers()); a "
             "program found is run only by a job you start, on a copy of the file, with a time limit.")


# -- the files ------------------------------------------------------------------------------------------------------

def _looks_like_endf(path):
    """Whether the file's first two lines read as ENDF-6 lines (pyazr.endf's line reader: 80 columns with MAT, MF
    and MT as integers) and the second names a material (MAT > 0: the TPID line may carry 0)."""
    from pyazr.endf.records import split_line
    from pyazr.errors import Error
    try:
        with open(path, "r", encoding="ascii", errors="strict") as f:
            lines = [f.readline(200) for _ in range(2)]
        if len(lines[1].rstrip("\r\n")) < 75:
            return False
        parsed = [split_line(path, k + 1, text) for k, text in enumerate(lines)]
    except (OSError, UnicodeDecodeError, Error):
        return False
    return parsed[1].mat > 0


def endf_files(directory, skip=()):
    """The ENDF-6 files under the project's ``directory`` (to :data:`SEARCH_DEPTH` levels, hidden directories and
    the directories ``skip`` left out), each ``{path (relative), size (bytes), mtime (s since the epoch)}``, sorted by
    path.  The one place the page's list of files is made: when the server has a general file listing, this keeps
    only the filter."""
    out, seen = [], 0
    skip = {os.path.realpath(s) for s in skip}
    for here, dirs, names in os.walk(directory):
        depth = os.path.relpath(here, directory).count(os.sep) + (here != directory)
        dirs[:] = sorted(d for d in dirs if not d.startswith(".") and depth < SEARCH_DEPTH
                         and os.path.realpath(os.path.join(here, d)) not in skip)
        for name in sorted(names):
            full = os.path.join(here, name)
            if name.startswith(".") or not os.path.isfile(full) or name.endswith((".azr", ".py")):
                continue
            seen += 1
            if seen > MAX_FILES:
                return out
            if _looks_like_endf(full):
                st = os.stat(full)
                out.append({"path": os.path.relpath(full, directory), "size": st.st_size, "mtime": st.st_mtime})
    return sorted(out, key=lambda f: f["path"])


def summary(path, name):
    """``pyazr.endf.read(path).summary().to_dict()`` with the name the page gave (``name``)."""
    from pyazr import endf
    return dict(endf.read(path).summary().to_dict(), name=name)


def checkers():
    """The checking codes found (``pyazr.endf.checkers()``: name, path, family, kind, what it does, what it does to
    the numbers), every program pyazr knows (:data:`pyazr.endf.PROGRAMS`), where they were looked for and how."""
    from pyazr import endf
    from pyazr.endf.checkers import CHECK_DEFAULT
    found = [dataclasses.asdict(c) for c in endf.checkers()]
    known = {name: dict(zip(("family", "kind", "does", "numbers"), row)) for name, row in endf.PROGRAMS.items()}
    search = {v: os.environ.get(v) or None for v in ("AZURE2_ENDF_UTILITIES", "AZURE2_PREPRO")}
    return {"found": found, "known": known, "search": search, "how": HOW_FOUND,
            "check_default": list(CHECK_DEFAULT)}


# -- the curves (called by the session worker's adapter) --------------------------------------------------------

def _plain(values):
    return [float(v) if np.isfinite(v) else None for v in np.asarray(values, float)]


def _parts(e):
    return {"error": type(e).__name__, "message": str(e), "what": e.what, "allowed": e.allowed, "fix": e.fix}


def _ratio(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _plain(np.where(b != 0, a / b, np.nan))


def curves_answer(files, session, evaluated, arrays):
    """The answer of the curves route: per file (``files``: [(name, absolute path)]) whether it was read; per data
    set of ``evaluated`` (``Session.evaluate(x).data``, with ``arrays``, the adapter's arrays of the same data sets)
    its points, datum times norm with its error and the model, and per file the curve (``Material.curves``: values in
    the data set's unit, NaN as None, coverage, source, notes) with the two ratios, or the refusal's parts."""
    from pyazr import endf
    from pyazr.errors import Error
    materials, shown = {}, []
    for name, path in files:
        try:
            materials[name] = endf.read(path)
            shown.append({"name": name, "refused": None})
        except Error as e:
            shown.append({"name": name, "refused": _parts(e)})
    sets = []
    for d, a in zip(evaluated, arrays):
        entry = {k: a[k] for k in ("name", "quantity", "engine_unit", "e_cm", "e_x", "angle_cm", "value_norm",
                                   "error_norm", "model")}
        entry["curves"] = {}
        for name, material in materials.items():
            try:
                c = material.curves(d, session=session)
                if c.unit != a["engine_unit"]:
                    raise endf.FormatRefused(what=f"the file's values in {c.unit}, the data set's in "
                                             f"{a['engine_unit']}", allowed="the data set's unit",
                                             fix="the page converts no unit: compare in Python with pyazr.endf")
                entry["curves"][name] = {"values": _plain(c.values), "coverage": [bool(v) for v in c.coverage],
                                         "unit": c.unit, "frame": c.frame, "source": c.source,
                                         "notes": list(c.notes), "covered": int(np.sum(c.coverage)),
                                         "ratio_to_data": _ratio(c.values, a["value_norm"]),
                                         "ratio_to_data_error": _ratio(np.abs(c.values) * np.asarray(
                                             a["error_norm"], float), np.square(a["value_norm"])),
                                         "ratio_to_model": _ratio(c.values, a["model"])}
            except Error as e:
                entry["curves"][name] = {"refused": _parts(e)}
        sets.append(entry)
    return {"files": shown, "datasets": sets,
            "ratios": {"ratio_to_data": "the file's value / (datum x normalization in use)",
                       "ratio_to_data_error": "the ratio to the datum times the datum's relative error",
                       "ratio_to_model": "the file's value / the model at the working x"}}


# -- the jobs (rows of jobs.JOBS; they run in a job's worker, without a session) ----------------------------------

def _counts(findings):
    return {s: sum(f["severity"] == s for f in findings) for s in ("error", "warning", "note")}


def check_job(tracker, x0, settings, names):
    """``pyazr.endf.check(path, codes, timeout)``: the module's own format check, then each checking code chosen
    (default: CHECKR, FIZCON, PSYCHE where found), each on a copy in a directory of its own; the answer's
    ``to_dict()`` with the counts per severity of each part and the counts each report states."""
    from pyazr import endf
    result = endf.check(settings["path"], codes=settings["codes"] or None, timeout=settings["timeout"])
    tracker.count(1, None, x0)
    out = result.to_dict()
    out["file"] = settings["path"]
    out["own"]["counts"] = _counts(out["own"]["findings"])
    for r in out["reports"]:
        r["counts"] = _counts(r["findings"])
    ran = ", ".join(r["program"] for r in out["reports"]) or "no checking code (none found)"
    errors = out["own"]["counts"]["error"] + sum(r["counts"]["error"] for r in out["reports"])
    return "finished", f"the own format check and {ran}: {errors} errors in all", {"check": out}


def standardize_job(tracker, x0, settings, names):
    """``pyazr.endf.standardize(path, program, mf2, timeout)`` on a copy: the candidate's path (outside the project
    folder until it is kept), its difference per section (``Candidate.diff()``), its summary and the program's
    terminal text and listing; what the program does to numbers (its row of :data:`pyazr.endf.PROGRAMS`)."""
    from pyazr import endf
    program = settings["program"] or "stanef"
    candidate = endf.standardize(settings["path"], program=program, mf2=settings["mf2"], timeout=settings["timeout"])
    tracker.count(1, None, x0)
    diff = [dataclasses.asdict(c) for c in candidate.diff()]
    changed = [c for c in diff if c["kind"] != "unchanged"]
    return "finished", f"{len(changed)} of {len(diff)} sections changed", {"candidate": {
        "file": settings["path"], "original": os.path.abspath(settings["path"]), "path": candidate.path,
        "program": program, "mf2": settings["mf2"], "diff": diff, "summary": candidate.summary(),
        "report": candidate.report, "numbers": endf.PROGRAMS[program][3]}}


def accept(state, record, to=None, overwrite=False):
    """Keep the candidate of a finished rewrite job (its ``record``): ``Candidate.accept(to=, overwrite=)`` with
    ``to`` inside the root (relative to the project's directory; none: beside the original as
    ``<name>.<program>.endf``), an existing file replaced only with ``overwrite`` (the page asks twice).  The project
    text is not changed (an ENDF file is not part of it): no history entry."""
    from pyazr import endf
    from .errors import Refused
    result = record.get("result") or {}
    cand = result.get("candidate")
    if record.get("state") != "finished" or not cand:
        raise Refused(f"job {record.get('id')} has no candidate file ({record.get('kind')} {record.get('method')}, "
                      f"{record.get('state')})", allowed="a finished rewrite job (endf_standardize)",
                      fix="run Correct on the file, then keep its candidate")
    original = state.inside_root(record["settings"]["path"])
    if not os.path.isfile(cand["path"]):
        raise Refused(f"the candidate {cand['path']} is gone (a temporary file of the job)",
                      allowed="a candidate that still exists", fix="run Correct again")
    target = state.inside_root(to) if to else None
    if target and not os.path.isdir(os.path.dirname(target)):
        raise Refused(f"{to}: no directory {os.path.dirname(to) or '.'} to write it in", allowed="a file in an "
                      "existing directory inside the root", fix="name a file beside the original")
    written = endf.Candidate(original, cand["path"], cand["program"], cand["report"]).accept(
        to=target, overwrite=bool(overwrite))
    state.inside_root(written)
    return {"path": os.path.relpath(written, state.directory), "replaced": os.path.realpath(written) == original,
            "overwrite": bool(overwrite)}


def job_rows():
    """The two jobs of the ENDF tab as rows of the job table (:class:`pyazr.ui.jobs.JobRow`): they build no session
    of the project, and their ``path`` setting is refused outside the root before they start."""
    from .jobs import JobRow, Setting
    return (
        JobRow("endf_check", "codes", "check an ENDF-6 file: pyazr.endf's own format check, then the checking codes "
               "found (pyazr.endf.check), each on a copy in a directory of its own", run=check_job, files=("path",),
               session=False,
               settings=(Setting("path", "text", "", "the ENDF-6 file, inside the project's root (relative to the "
                                 "project's directory)"),
                         Setting("codes", "names", [], "the checking codes to run (checkr, fizcon, psyche, inter); "
                                 "empty: CHECKR, FIZCON and PSYCHE where found"),
                         Setting("timeout", "integer", 300, "the time limit of each code's run (s)")),
               facts=(("pyazr call", "pyazr.endf.check"), ("writes", "nothing in the project folder"))),
        JobRow("endf_standardize", "program", "rewrite an ENDF-6 file with a code that rewrites (STANEF, FIXUP, "
               "DICTIN, ENDF2C) on a copy: a candidate file and its difference per section, kept only on request "
               "(pyazr.endf.standardize)", run=standardize_job, files=("path",), session=False,
               settings=(Setting("path", "text", "", "the ENDF-6 file, inside the project's root (relative to the "
                                 "project's directory)"),
                         Setting("program", "text", "stanef", "the code that rewrites: stanef, fixup, dictin or "
                                 "endf2c"),
                         Setting("mf2", "boolean", False, "STANEF: also the numbers of MF=2/MT=151 (its options "
                                 "0,Y,Y,Y)"),
                         Setting("timeout", "integer", 300, "the time limit of the run (s)")),
               facts=(("pyazr call", "pyazr.endf.standardize"),
                      ("writes", "a candidate outside the project folder; POST /api/endf/accept keeps it"))),
    )
