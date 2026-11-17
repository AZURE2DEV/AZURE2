"""The page's utilities on pyazr's modules: the covariance builder (:mod:`pyazr.covariance`) and the reaction rate's
tables and its job (:mod:`pyazr.rate`; :func:`rate_engine`, a row of ``jobs.JOBS``). The arithmetic is pyazr's; these
functions read the open project, hand it over, and write what the page asks to write as one history entry.

Where a builder's parts are kept: in the covariance file itself, on its ``#pyazr.covariance`` comment line
(:func:`pyazr.covariance.read_parts`), so that the builder reopens with them; nothing is kept in the browser or in
the server beyond the project's history.
"""
import os

import numpy as np

from .errors import Refused
from ..paths import is_absolute


def _with_covariance(project):
    """The records whose kind has a covariance field (the schema's): the data sets a covariance can join."""
    import pyazr
    kinds = {k for k, row in pyazr.schema()["records"].items() if COVARIANCE_FIELD in row["fields"]}
    return [n for n, r in project.records().items() if r.record_kind in kinds]


#: The field of a record that names its covariance file.
COVARIANCE_FIELD = "covariance"


def _datasets(project, datasets):
    names = [datasets] if isinstance(datasets, str) else list(datasets or [])
    if not names:
        raise Refused("name at least one data set", allowed="the data sets of the project")
    known = _with_covariance(project)
    unknown = [n for n in names if n not in known]
    if unknown:
        raise Refused(f"not data sets of the project: {', '.join(unknown)}", allowed=", ".join(known), known=known)
    return names


def covariance_open(state, name):
    """What the builder opens with for the data set ``name``: its covariance field, the parts of the file it names
    (read from the file's parts line; None for a file written by hand), the data sets naming the same file, a file
    name to propose (the field's, else ``cov_<name>.txt`` not yet taken) and the kinds of part."""
    from pyazr import covariance
    with state.lock:
        project = state.project
        _datasets(project, [name])
        field = project[name].get(COVARIANCE_FIELD)
        sharing = [n for n in _with_covariance(project) if field and project[n].get(COVARIANCE_FIELD) == field]
        directory = state.directory
    saved, hand_written = None, False
    if field:
        path = os.path.join(directory, field)
        if os.path.isfile(path):
            saved = covariance.read_parts(path)
            hand_written = saved is None
    proposed = field
    if not proposed:
        stem, k = f"cov_{name}", 1
        proposed = f"{stem}.txt"
        while os.path.exists(os.path.join(directory, proposed)):
            k += 1
            proposed = f"{stem}_{k}.txt"
    parts = saved["parts"] if saved else [{"kind": "own"}]
    datasets = saved["datasets"] if saved else (sharing or [name])
    return {"name": name, "covariance": field, "hand_written": hand_written, "parts": parts,
            "datasets": datasets, "proposed": proposed, "kinds": _kinds(covariance)}


def _kinds(covariance):
    return {k: {"doc": v["doc"], "fields": [dict(zip(("name", "type", "default", "unit", "doc"), f))
                                            for f in v["fields"]]} for k, v in covariance.PART_KINDS.items()}


def covariance_build(state, datasets, parts):
    """The covariance of ``datasets`` from ``parts`` on the current text (:func:`pyazr.covariance.build`): the
    summary the page draws (total error per point beside the file's, correlation on the fitted rows, checks, the
    refusal and its parts, notes)."""
    from pyazr import Project, covariance
    with state.lock:
        project = Project(state.project.to_bytes(), state.path)
    names = _datasets(project, datasets)
    return covariance.build(project, names, parts).summary()


def covariance_write(state, datasets, parts, path, rev=None):
    """Build the covariance, write it at ``path`` (relative to the project's directory, inside the root) and set
    the ``covariance`` field of every data set of ``datasets`` to it: one history entry, which undo takes back with
    the file (removed, or restored as it was).  Refused when the engine would refuse the matrix."""
    from pyazr import Project, covariance
    if not path or is_absolute(str(path)):
        raise Refused(f"{path!r}: a file name relative to the project's directory", allowed="a relative path",
                      fix="give the name the builder proposes, or another inside the project's directory")
    with state.lock:
        project = Project(state.project.to_bytes(), state.path)
        names = _datasets(project, datasets)
        state.inside_root(path)
        cov = covariance.build(project, names, parts)
        if cov.refusal:
            raise Refused(f"the covariance is not written: {cov.refusal}", allowed="a positive definite matrix with "
                          "every fitted row's own variance", fix="change the parts the refusal names")
        text = cov.text()

        def change(p):
            for n in names:
                p.set(n, **{COVARIANCE_FIELD: str(path)})

        state.edit_with_files(rev, f"covariance {path} for {', '.join(names)}", change, {path: text})
        return {"rev": state.rev, "path": path, "datasets": names}


def rate_grids():
    """The named temperature grids of :mod:`pyazr.rate` (name -> doc, temperatures in GK) and the rate's units."""
    from pyazr import rate
    return {"grids": {k: {"doc": doc, "T9": list(t)} for k, (doc, t) in rate.TEMPERATURE_GRIDS.items()},
            "units": rate.UNITS}


# -- the rate job (a row of jobs.JOBS; it runs in the job's worker) ----------------------------------------------

def _refusal(what, allowed, fix):
    """A pyazr.ParameterError with its three parts (each a part of the message)."""
    from pyazr.errors import ParameterError
    return ParameterError(f"{what}; allowed: {allowed}; {fix}", what=what, allowed=allowed, fix=fix)


def rate_engine(tracker, x0, settings, names):
    """The reaction rate N_A<sigma v> of an entrance and an exit pair at the start point, computed by the engine as the
    program's rate command computes it (:mod:`pyazr.rate`), on a named or a custom temperature grid; with the band
    of the source job's covariance (linear in ln rate) or of its chain's draws (16th, 50th, 84th percentiles); the
    Gamow window per temperature; a comparison table of the project folder with the ratio; the table as CSV in the
    output directory.  Stops between the band's rates."""
    from pyazr import rate
    s = tracker.session.s
    grid = settings["grid"] or "table"
    if grid == "custom":
        temps = {"t_min": settings["t_min"], "t_max": settings["t_max"], "t_step": settings["t_step"]}
    elif grid in rate.TEMPERATURE_GRIDS:
        temps = {"grid": grid}
    else:
        raise _refusal(f"no temperature grid {grid!r}", ", ".join(list(rate.TEMPERATURE_GRIDS) + ["custom"]),
                       "choose a named grid, or custom with t_min, t_max and t_step (GK)")
    band_kind = settings["band"] or "none"
    if band_kind not in ("none", "covariance", "draws"):
        raise _refusal(f"no band {band_kind!r}", "none, covariance, draws",
                       "covariance after a fit with a covariance, draws after a sampler")
    result = (tracker.source or {}).get("result") or {}
    if band_kind == "covariance" and not result.get("covariance"):
        raise _refusal("no covariance to propagate: a band from a covariance starts from a fit job",
                       "a fit job's best point with its covariance",
                       "start from a fit job that gave a covariance (iminuit MIGRAD + HESSE, the engine's)")
    if band_kind == "draws" and not result.get("chain"):
        raise _refusal("no chain to draw from: a band from draws starts from a sampler job",
                       "a sampler job's best point with its chain", "start from a sampler job")
    comparison = None
    if settings["comparison"]:
        path = os.path.realpath(os.path.join(s.cwd, settings["comparison"]))
        if os.path.commonpath([path, os.path.realpath(s.cwd)]) != os.path.realpath(s.cwd):
            raise _refusal(f"{settings['comparison']} is outside the project folder", f"a file inside {s.cwd}",
                           "copy the table into the project folder")
        comparison = rate.read_rate_table(path)
    values = dict(zip(names, np.asarray(x0, float).tolist()))
    source = (f"the best point of job {tracker.source['id']}" if tracker.source
              else "the job's start point (the working x of its revision)")
    with rate.RateSession(s.file, settings["entrance"], settings["exit"], cwd=s.cwd) as rs:
        best = rs.rate(values, source=source, **temps)
        tracker.count(1, None, x0)
        band = None
        progress = lambda done, total: tracker.count(1 + done, None, x0)  # noqa: E731
        if band_kind == "covariance":
            cov = result["covariance"]
            band = rate.rate_band(rs, values, covariance=cov["matrix"], names=cov["names"], progress=progress,
                                  stop=tracker.stop_requested, **temps)
        elif band_kind == "draws":
            chain = np.load(result["chain"])
            flat = chain[result["discard"]:].reshape(-1, chain.shape[2])
            band = rate.rate_band(rs, values, draws=flat, names=names, max_draws=settings["draws"],
                                  seed=settings["seed"], progress=progress, stop=tracker.stop_requested, **temps)
    gamow = rate.gamow_window(s.file, settings["entrance"], best.T9)
    path = tracker.product("rate", ".csv")
    rate.write_rate_table(path, best, band=band, comparison=comparison)
    with open(path, encoding="utf-8") as fh:
        table = fh.read()
    out = {"rate": {"T9": best.T9.tolist(), "rate": best.rate.tolist(), "reaction": best.reaction,
                    "entrance": best.entrance, "exit": best.exit, "source": best.source, "units": best.units},
           "gamow": {k: np.asarray(v, float).tolist() for k, v in gamow.items()}, "table": table, "csv": path}
    if band is not None:
        out["band"] = {"low": band.low.tolist(), "median": band.median.tolist(), "high": band.high.tolist(),
                       "how": band.how, "varied": band.varied, "draws": band.draws, "method": band.method}
    if comparison is not None:
        out["comparison"] = {k: (None if v is None else np.asarray(v, float).tolist()) if k != "path" else v
                             for k, v in comparison.items()}
    # the ratios the page draws under the rate: each curve over the median (the band's, else the rate itself)
    median = band.median if band is not None else best.rate
    ratios = {"rate": (best.rate / median).tolist()}
    if band is not None:
        ratios.update(low=(band.low / median).tolist(), high=(band.high / median).tolist())
    if comparison is not None:
        ref = rate.table_at(comparison, best.T9)
        ratios["comparison"] = (ref / median).tolist()
        out["comparison"]["at_T9"] = ref.tolist()
    out["ratios"] = ratios
    out["formats"] = _library_formats(s.file, settings, best, band, tracker)
    how = f"N_A<sigma v> of {best.reaction} at {best.T9.size} temperatures" + (f"; {band.how}" if band else "")
    return "finished", how, out


def _library_formats(file, settings, best, band, tracker):
    """The rate as a STARLIB entry and as REACLIB sets (:mod:`pyazr.rate_formats`), each the file's text with what
    it holds, or the refusal's three parts (STARLIB without a band or off its grid; REACLIB off 0.01-10 GK or above
    the tolerance); the species, chapter and Q from the project's pairs."""
    from pyazr import rate_formats as rf
    from pyazr.errors import Error
    out = {}
    try:
        species = rf.species_of(file, settings["entrance"], settings["exit"])
    except Error as e:
        refusal = {"refused": str(e), "what": e.what, "allowed": e.allowed, "fix": e.fix}
        return {"starlib": refusal, "reaclib": refusal}
    out["species"] = species

    def attempt(name, make):
        path = tracker.product(name, ".dat")
        try:
            facts = make(path)
            with open(path, encoding="utf-8") as fh:
                out[name] = {"text": fh.read(), "path": path, **facts}
        except Error as e:
            out[name] = {"refused": str(e), "what": e.what, "allowed": e.allowed, "fix": e.fix}

    def starlib(path):
        rf.write_starlib(path, best, band, species, label=settings["starlib_label"], flag=settings["starlib_flag"])
        return {"uncertainty": "f.u. = exp(sigma) of the covariance band" if band.method == "covariance"
                else "f.u. = (x84 / x16)^(1/2) of the draws, the median as the rate"}

    def reaclib(path):
        fit = rf.fit_reaclib(best, tolerance=settings["reaclib_tolerance"], max_sets=settings["reaclib_sets"])
        rf.write_reaclib(path, fit, species, label=settings["reaclib_label"], flag=settings["reaclib_flag"])
        return {"summary": fit.summary(), "sets": fit.n_sets, "max_deviation": fit.max_deviation,
                "deviation": fit.deviation.tolist(), "T9": fit.T9.tolist()}

    attempt("starlib", starlib)
    attempt("reaclib", reaclib)
    return out
