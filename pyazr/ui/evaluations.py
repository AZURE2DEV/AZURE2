"""The page's ENDF screens on :mod:`pyazr.endf`: the ENDF-6 files of the project folder, what a file holds
(``Material.summary()``), the file's values at every data set's points and on the grid of its model curve
(``Material.curves``, through the session worker: :meth:`pyazr.ui.engine.Engine.endf_curves`), the checking codes
found (``pyazr.endf.checkers()``), and the two jobs a user starts with a button: the check of a file
(``pyazr.endf.check``) and a rewrite by a code (``pyazr.endf.standardize``) whose candidate replaces nothing until it
is kept (``Candidate.accept``).

The format work and every number are pyazr's.  Two ratios are formed here, named in the answer: the file's value
over the datum times its normalization (``ratio_to_data``, with ``ratio_to_data_error``, the ratio times the
datum's relative error) and over the model (``ratio_to_model``), at each point where the file has a value.  The
files a user names are checked to lie inside the server's root by the server before any of these functions sees
them (the state rule, docs/upgrade/web_ui.md section 5).
"""
import dataclasses
import os

import numpy as np

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


def endf_files(directory, listing):
    """The ENDF-6 files among ``listing`` (the rows of the server's file listing, ``GET /api/files``: paths relative
    to ``directory``), each ``{path (relative), size (bytes), mtime (s since the epoch)}``, sorted by path.  The one
    place that says which files of the folder the ENDF tab offers: a filter over the page's one file listing."""
    out = []
    for row in listing:
        full = os.path.join(directory, row["path"])
        if row.get("kind") in ("project", "binary") or row["path"].endswith(".py") or not _looks_like_endf(full):
            continue
        out.append({"path": row["path"], "size": row["size"], "mtime": os.stat(full).st_mtime})
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


def _curve_entry(material, data, session, unit):
    """One file's values at ``data``'s points (an ``EvaluatedData``, or a ``ModelCurve``: ``Material.curves`` takes
    both), checked to be in the data set's unit."""
    from pyazr import endf
    c = material.curves(data, session=session)
    if c.unit != unit:
        raise endf.FormatRefused(what=f"the file's values in {c.unit}, the data set's in {unit}",
                                 allowed="the data set's unit",
                                 fix="the page converts no unit: compare in Python with pyazr.endf")
    return c


def curves_answer(files, session, evaluated, arrays, model_curves=None):
    """The answer of the curves route: per file (``files``: [(name, absolute path)]) whether it was read; per data
    set of ``evaluated`` (``Session.evaluate(x).data``, with ``arrays``, the adapter's arrays of the same data sets)
    its points, datum times norm with its error and the model, and per file the curve (``Material.curves``: values in
    the data set's unit, NaN as None, coverage, source, notes) with the two ratios, or the refusal's parts.  Where
    ``model_curves`` gives a data set's model curve (``{name: (ModelCurve, its arrays) or the refusal}``), the entry
    has it as ``curve`` (its arrays, or ``refused``), and each file's values on the curve's grid as ``on_curve``
    (values, coverage; or the refusal)."""
    from pyazr import endf
    from pyazr.errors import Error
    materials, shown = {}, []
    for name, path in files:
        try:
            materials[name] = endf.read(path)
            shown.append({"name": name, "refused": None})
        except Error as e:
            shown.append({"name": name, "refused": _parts(e)})
    model_curves = model_curves or {}
    sets = []
    for d, a in zip(evaluated, arrays):
        entry = {k: a[k] for k in ("name", "quantity", "engine_unit", "e_cm", "e_x", "angle_cm", "value_norm",
                                   "error_norm", "model")}
        entry["curves"] = {}
        model_curve = model_curves.get(d.name)
        if isinstance(model_curve, Exception):
            entry["curve"] = {"refused": _parts(model_curve)}
            model_curve = None
        elif model_curve is not None:
            model_curve, entry["curve"] = model_curve
        for name, material in materials.items():
            try:
                c = _curve_entry(material, d, session, a["engine_unit"])
                entry["curves"][name] = {"values": _plain(c.values), "coverage": [bool(v) for v in c.coverage],
                                         "unit": c.unit, "frame": c.frame, "source": c.source,
                                         "notes": list(c.notes), "covered": int(np.sum(c.coverage)),
                                         "ratio_to_data": _ratio(c.values, a["value_norm"]),
                                         "ratio_to_data_error": _ratio(np.abs(c.values) * np.asarray(
                                             a["error_norm"], float), np.square(a["value_norm"])),
                                         "ratio_to_model": _ratio(c.values, a["model"])}
            except Error as e:
                entry["curves"][name] = {"refused": _parts(e)}
                continue
            if model_curve is not None:
                try:
                    g = _curve_entry(material, model_curve, session, a["engine_unit"])
                    entry["curves"][name]["on_curve"] = {"values": _plain(g.values),
                                                         "coverage": [bool(v) for v in g.coverage]}
                except Error as e:
                    entry["curves"][name]["on_curve"] = {"refused": _parts(e)}
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
    ``<name>.<program>.endf``), an existing file replaced only with ``overwrite`` (the page asks twice) and only
    when it is itself an ENDF-6 file (never the project file or a data file).  The project
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
    if target and os.path.exists(target) and not _looks_like_endf(target):
        raise Refused(f"{to} exists and is not an ENDF-6 file: a candidate replaces an ENDF-6 file only",
                      allowed="a new file, or (with overwrite) an existing ENDF-6 file inside the root",
                      fix="name a new file beside the original")
    written = endf.Candidate(original, cand["path"], cand["program"], cand["report"]).accept(
        to=target, overwrite=bool(overwrite))
    state.inside_root(written)
    return {"path": os.path.relpath(written, state.directory), "replaced": os.path.realpath(written) == original,
            "overwrite": bool(overwrite)}


# -- a new project from a file, a file from the project (pyazr.endf.to_project, pyazr.endf.write) -------------------

def _refusal(what, allowed, fix):
    from pyazr.errors import ParameterError
    return ParameterError(f"{what}; allowed: {allowed}; {fix}", what=what, allowed=allowed, fix=fix)


def resonances(path, name, isotope=0, energy_range=0):
    """What MF=2 of a file holds for a new project (``read(path).resonances(isotope, energy_range)``), nothing
    written: the record's text, EL and EH (eV), the incident pair, the pairs (MT, reaction) and the spin groups with
    their counts of levels and channels; or the refusal's three parts (``refused``)."""
    from pyazr import endf
    from pyazr.endf import names
    from pyazr.errors import Error
    try:
        record = endf.read(path).resonances(isotope=int(isotope), energy_range=int(energy_range))
        entrance = record.entrance
    except Error as e:
        return {"name": name, "isotope": int(isotope), "energy_range": int(energy_range), "refused": _parts(e)}
    return {"name": name, "isotope": int(isotope), "energy_range": int(energy_range), "refused": None,
            "text": str(record), "el_eV": record.el, "eh_eV": record.eh, "entrance": entrance,
            "pairs": [{"pair": k, "mt": p.mt, "reaction": names.reaction_name(p.mt), "q_eV": p.q}
                      for k, p in enumerate(record.pairs, start=1)],
            "groups": [{"jpi": g.jpi, "levels": len(g.levels), "channels": len(g.channels)} for g in record.groups]}


def _separation(text):
    """The separation-energy setting: ``"ame2020"`` (or empty) or a number in MeV."""
    text = (text or "ame2020").strip()
    if text.lower() == "ame2020":
        return "ame2020"
    try:
        return float(text)
    except ValueError:
        raise _refusal(f"the separation energy {text!r}", "ame2020, or a number in MeV",
                       "leave it empty for AME2020's") from None


def import_job(tracker, x0, settings, names):
    """``pyazr.endf.to_project`` of the file's MF=2 (``Material.resonances(isotope, energy_range)``) into a new project
    in a directory of the job's own, outside the project folder: the report of every conversion (the separation
    energy and its source, each mass against AME2020, the pairs, the levels with their Brune iteration, what was not
    written, the notes) and the new project's text; :func:`keep_project` puts it into the folder."""
    import tempfile
    from pyazr import endf
    m = endf.read(settings["path"])
    record = m.resonances(isotope=settings["isotope"] or 0, energy_range=settings["energy_range"] or 0)
    name = settings["name"] or os.path.splitext(os.path.basename(settings["path"]))[0] + ".azr"
    if not name.endswith(".azr"):
        raise _refusal(f"the new project's name {name!r}", "a file name ending in .azr", "name it <name>.azr")
    work = tempfile.mkdtemp(prefix="pyazr_endf_project_")
    project, conversion = endf.to_project(record, os.path.join(work, os.path.basename(name)),
                                          separation_energy=_separation(settings["separation_energy"]), material=m)
    tracker.count(1, None, x0)
    report = dataclasses.asdict(conversion)
    return "finished", f"{len(conversion.levels)} levels in {len(conversion.pairs)} pairs", {"project": {
        "file": settings["path"], "name": name, "candidate": project.path, "report": report,
        "report_text": str(conversion), "record_text": str(record), "text": project.text}}


def keep_project(state, record, to=None):
    """Put the new project of a finished ``endf_import`` job (its ``record``) into the project folder: at ``to``
    (relative to the project's directory; none: the job's name), a new file ending in ``.azr`` -- never an existing
    file, never the open project.  The open project is not changed; the page then offers to open the new one."""
    import shutil
    from .errors import Refused
    result = (record.get("result") or {}).get("project")
    if record.get("state") != "finished" or not result:
        raise Refused(f"job {record.get('id')} has no new project ({record.get('kind')} {record.get('method')}, "
                      f"{record.get('state')})", allowed="a finished endf_import job",
                      fix="run New project from this file, then keep its project")
    name = to or result["name"]
    if not str(name).endswith(".azr"):
        raise Refused(f"{name}: a project file ends in .azr", allowed="a new file <name>.azr", fix="name it .azr")
    target = state.inside_root(name)
    if os.path.exists(target):
        raise Refused(f"{name} exists: a new project is a new file (the open project and any other are kept)",
                      allowed="a file name not in use", fix="give another name")
    if not os.path.isdir(os.path.dirname(target)):
        raise Refused(f"{name}: no directory {os.path.dirname(name) or '.'} to write it in",
                      allowed="a file in an existing directory inside the root", fix="name a file beside the project")
    if not os.path.isfile(result["candidate"]):
        raise Refused(f"the job's project {result['candidate']} is gone (a temporary file of the job)",
                      allowed="a project that still exists", fix="run New project from this file again")
    shutil.copyfile(result["candidate"], target)
    # open_path: what POST /api/project/open takes (it resolves a relative path in the server's working directory)
    return {"path": os.path.relpath(target, state.directory), "open_path": target}


def export_preview(session, x):
    """What an ENDF-6 file of the session at x would lose, before it is written: ``pyazr.endf.from_session(s, x)``'s
    losses and photon widths (or its refusal's three parts), the sections known not to be written
    (``pyazr.endf.writer.not_written_known``), and the energy range the writer takes by default (the data's lab
    energies of the projectile, MeV)."""
    from pyazr import endf
    from pyazr.errors import Error
    lab = [e for d in session.evaluate(x).data for e in np.asarray(d.extra.get("energy_lab", []), float)]
    out = {"energy_range_MeV": [min(lab), max(lab)] if lab else None}
    try:
        record, losses = endf.from_session(session, x)
    except Error as e:
        return dict(out, refused=_parts(e), losses=[], photon_widths=[], not_written=[])
    widths = [{"jpi": jpi, "energy_MeV": e, "pair": k, "multipole": f"{r}{l}", "width_eV": w}
              for jpi, e, k, r, l, w in losses.photon_widths]
    # what is known before the run; the job's result lists every section the run then left out (Written.not_written)
    return dict(out, refused=None, losses=list(losses.items), photon_widths=widths,
                not_written=list(endf.writer.not_written_known(record)))


def write_job(tracker, x0, settings, names):
    """``pyazr.endf.write`` of the job's session at its start point into a new file of the project folder (never
    over an existing one): its sections, the grid of each MF=3 section, the check that MF=2 reproduces MF=3, the
    losses and what is not written."""
    from pyazr import endf
    s = tracker.session.s
    name = settings["name"]
    if not name:
        raise _refusal("no file name", "a new file of the project folder, <name>.endf", "name the file")
    lo, hi = settings["e_low"], settings["e_high"]
    if (lo is None) != (hi is None):
        raise _refusal("one end of the energy range", "both ends (lab MeV), or neither for the data's range",
                       "give e_low and e_high")
    header = {key: settings[setting] for key, setting in HEADER_SETTINGS
              if settings.get(setting) not in (None, "")}
    # the writer names the session's file in MF=1: a session of the same text under the project's name
    from .engine import session_named
    from .state import project_file_of
    with session_named(s, project_file_of(s.file) or os.path.basename(s.file)) as named:
        written = endf.write(named, os.path.join(s.cwd, name), np.asarray(x0, float), header=header,
                             energy_range_MeV=None if lo is None else (lo, hi),
                             tolerance=settings["tolerance"] or 1e-3, number_form=settings["number_form"] or "widest")
    tracker.count(1, None, x0)
    path = os.path.relpath(written.path, s.cwd)
    return "finished", f"{path}: {len(written.sections)} sections", {"written": {
        "path": path, "sections": [list(x) for x in written.sections],
        "grids": {str(mt): list(g) for mt, g in written.grids.items()},
        "mf2_check": {str(mt): v for mt, v in written.mf2_check.items()}, "losses": list(written.losses),
        "not_written": list(written.not_written)}}


#: The MF=1 fields of an export (pyazr.endf.write's header keys, p. 42) and the job's settings that give them.
HEADER_SETTINGS = (("MAT", "mat"), ("AUTH", "authors"), ("ALAB", "lab"), ("EDATE", "date"), ("REF", "reference"),
                   ("LIBRARY", "library"))


def job_rows():
    """The jobs of the ENDF tab as rows of the job table (:class:`pyazr.ui.jobs.JobRow`): the check, the rewrite and
    the new project build no session of the project; the export builds the session of the job's text.  A setting
    naming a file (``path``, ``name``) is refused outside the root before the job starts."""
    from .jobs import JobRow, Setting
    return (
        JobRow("endf_check", "codes", "check an ENDF-6 file: pyazr.endf's own format check, then the checking codes "
               "found (pyazr.endf.check), each on a copy in a directory of its own", run=check_job, files=("path",),
               session=False, tab="endf",
               settings=(Setting("path", "text", "", "the ENDF-6 file, inside the project's root (relative to the "
                                 "project's directory)"),
                         Setting("codes", "names", [], "the checking codes to run (checkr, fizcon, psyche, inter); "
                                 "empty: CHECKR, FIZCON and PSYCHE where found"),
                         Setting("timeout", "integer", 300, "the time limit of each code's run (s)")),
               facts=(("pyazr call", "pyazr.endf.check"), ("writes", "nothing in the project folder"))),
        JobRow("endf_standardize", "program", "rewrite an ENDF-6 file with a code that rewrites (STANEF, FIXUP, "
               "DICTIN, ENDF2C) on a copy: a candidate file and its difference per section, kept only on request "
               "(pyazr.endf.standardize)", run=standardize_job, files=("path",), session=False, tab="endf",
               settings=(Setting("path", "text", "", "the ENDF-6 file, inside the project's root (relative to the "
                                 "project's directory)"),
                         Setting("program", "text", "stanef", "the code that rewrites: stanef, fixup, dictin or "
                                 "endf2c"),
                         Setting("mf2", "boolean", False, "STANEF: also the numbers of MF=2/MT=151 (its options "
                                 "0,Y,Y,Y)"),
                         Setting("timeout", "integer", 300, "the time limit of the run (s)")),
               facts=(("pyazr call", "pyazr.endf.standardize"),
                      ("writes", "a candidate outside the project folder; POST /api/endf/accept keeps it"))),
        JobRow("endf_import", "to_project", "a new project from an ENDF-6 file's resonance parameters (MF=2, LRF=7): "
               "pyazr.endf.to_project in a directory of the job's own, the report of every conversion; the "
               "project is put into the folder only when kept", run=import_job, files=("path",), session=False,
               tab="endf",
               settings=(Setting("path", "text", "", "the ENDF-6 file, inside the project's root (relative to the "
                                 "project's directory)"),
                         Setting("isotope", "integer", 0, "the isotope of MF=2 (0: the first)"),
                         Setting("energy_range", "integer", 0, "the energy range of the isotope (0: the first)"),
                         Setting("separation_energy", "text", "ame2020", "the incident pair's separation energy "
                                 "from the compound nucleus: ame2020, or a number in MeV"),
                         Setting("name", "text", "", "the new project's file name (.azr; empty: the ENDF file's "
                                 "name with .azr)")),
               facts=(("pyazr call", "pyazr.endf.to_project"),
                      ("writes", "a project outside the project folder; POST /api/endf/keep_project keeps it"))),
        JobRow("endf_write", "writer", "the project at the job's start point as an ENDF-6 file (pyazr.endf.write: "
               "MF=1, MF=2 LRF=7, MF=3 of every exit and of capture), a new file of the project folder",
               run=write_job, files=("name",), tab="endf",
               settings=(Setting("name", "text", "", "the new file, relative to the project's directory"),
                         Setting("mat", "integer", 9999, "MAT, the material number (p. 42)"),
                         Setting("authors", "text", "", "AUTH, the evaluators"),
                         Setting("lab", "text", "", "ALAB, the laboratory"),
                         Setting("date", "text", "", "EDATE, the evaluation's date"),
                         Setting("reference", "text", "", "REF, the reference"),
                         Setting("library", "text", "AZURE2", "the library's name (HSUB's first line)"),
                         Setting("e_low", "number", None, "the lowest lab energy of the projectile (MeV; empty: "
                                 "the data's)"),
                         Setting("e_high", "number", None, "the highest lab energy of the projectile (MeV; empty: "
                                 "the data's)"),
                         Setting("tolerance", "number", 1e-3, "the grid's tolerance: linear interpolation between "
                                 "its points within this (relative)"),
                         Setting("number_form", "text", "widest", "the number form: widest, classic or e "
                                 "(pyazr.endf.NUMBER_FORMS)")),
               facts=(("pyazr call", "pyazr.endf.write"), ("writes", "a new ENDF-6 file of the project folder"))),
    )
