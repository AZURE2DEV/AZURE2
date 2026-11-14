"""What the package, the engine and a session say about themselves, as plain data (docs/upgrade/pyazr3.md section 2.5).

- :func:`describe_package` (``pyazr.describe()`` without a path): the version, the engine's build, the Coulomb table,
  the stopping and mass data in use, the thread count.
- :func:`kinds` (``pyazr.kinds()``): the record kinds, data kinds (quantities), parameter kinds, session options,
  bound physics functions and error classes.  Every table is read from the engine's own tables at the first call
  (the schema of ``AZURE2 --schema``, ``ParamKindTable``, the binding's option and physics tables); none is
  written a second time here.
- :func:`describe_session` (``pyazr.Session.describe()``): what one session holds.
- The data sets and extrapolations of a format-3 project as the engine built them (:class:`DescribedDataset`, ...):
  what the engine made of the records (``AZURE2 --describe --json``) -- the quantity, the observable it computes, the
  reaction tree with the ``<levels>`` pair keys, the lab-energy window, the data file and its rows read and kept, the
  effects and the fit parameters; a calculation session hands them out (:attr:`pyazr.Session.datasets`,
  :attr:`~pyazr.Session.active_datasets`, :attr:`~pyazr.Session.extrapolations`).  Nothing here reads a project file.

Every result is made of dicts, lists, strings, numbers, booleans and None: ``json.dumps`` takes it.
"""

import copy
import hashlib
import json
import os
import sys
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .errors import ParameterKeyError, close_names, listed

__all__ = ["DescribedDataset", "DescribedDatasets", "DescribedExtrapolation", "DescribedExtrapolations"]


def _engine():
    from . import _azure2
    return _azure2


def describe_project(path, cwd):
    """The engine's description of the format-3 project ``path`` as a dict (``AZURE2 --describe --json``), its data
    files read relative to ``cwd`` as a session opens them.  RuntimeError when the engine does not describe it."""
    with open(path, "rb") as f:
        text = f.read()
    refusal, diagnostics, description = _engine().describe_format3_text(text, path, True, cwd, "")
    if refusal or not description:
        raise RuntimeError(f"{path}: the engine does not describe the project: {refusal or diagnostics}")
    return json.loads(description)


def session_records(kind, path, cwd):
    """What a calculation session hands out as ``kind`` ("datasets" or "extrapolations") for the format-3 project
    ``path``: the engine's description (:class:`DescribedDatasets`, :class:`DescribedExtrapolations`)."""
    described = DescribedDatasets if kind == "datasets" else DescribedExtrapolations
    return described.from_description(describe_project(path, cwd))


def _first_term(entry):
    terms = entry.get("terms") or []
    return (int(terms[0]["entrance_pair"]), int(terms[0]["exit_pair"])) if terms else (0, 0)


@dataclass(frozen=True)
class DescribedDataset:
    """One ``dataset`` record as the engine built it (one entry of ``describe(as_json=True)["datasets"]``).

    ``key`` is its number among the data sets (1-based, file order), the label of its outputs; ``name`` the record's
    name.  ``entrance_key``/``exit_key`` are the ``<levels>`` pair keys of the reaction's first term (-1: every capture
    pair); ``terms`` every term as the engine resolved it.  ``energy_min``/``energy_max`` the lab-energy window (MeV;
    ``energy_max`` is the largest double when the record has no upper end), ``angle_min``/``angle_max`` the angle window
    (deg) of an angle-dependent observable, else None.  ``rows_read``/``rows_kept`` the rows of the data file and those
    inside the windows (0 for an inactive data set, whose file is not read).  ``parameters``: dicts with ``name``,
    ``kind``, ``status`` (fit, fixed, shared, profiled)."""
    key: int
    name: str
    active: bool
    quantity: str
    observable: str
    tree: str
    entrance_key: int
    exit_key: int
    terms: tuple
    unit: str
    engine_unit: str
    energy_min: float
    energy_max: float
    angle_min: Optional[float]
    angle_max: Optional[float]
    data_file: str
    data_file_resolved: str
    rows_read: int
    rows_kept: int
    effects: tuple
    sample: str
    parameters: tuple

    @classmethod
    def from_json(cls, d):
        entrance, exit_ = _first_term(d)
        lo, hi = d.get("E_lab_window_MeV", [0.0, float("inf")])
        angles = d.get("angle_window_deg")
        return cls(key=int(d["key"]), name=d["name"], active=bool(d["active"]), quantity=d.get("quantity", ""),
                   observable=d.get("observable", ""), tree=d.get("tree", ""), entrance_key=entrance, exit_key=exit_,
                   terms=tuple(dict(t) for t in d.get("terms", [])), unit=d.get("unit", ""),
                   engine_unit=d.get("engine_unit", ""), energy_min=float(lo), energy_max=float(hi),
                   angle_min=float(angles[0]) if angles else None, angle_max=float(angles[1]) if angles else None,
                   data_file=d.get("data_file", ""), data_file_resolved=d.get("data_file_resolved", ""),
                   rows_read=int(d.get("rows_read", 0)), rows_kept=int(d.get("rows_kept", 0)),
                   effects=tuple(d.get("effects", [])), sample=d.get("sample", ""),
                   parameters=tuple(dict(p) for p in d.get("parameters", [])))

    def reaction(self, pairs=None):
        """The reaction as written (its tree); ``pairs`` is accepted and not needed (the tree names the pairs)."""
        return self.tree

    def describe(self):
        """One line: name, quantity, reaction with pair keys, rows kept of rows read."""
        return (f"{self.name} (data set {self.key}{'' if self.active else ', inactive'}): {self.quantity} of "
                f"{self.tree} [pairs {self.entrance_key} -> {self.exit_key}], {self.rows_kept} of {self.rows_read} "
                f"rows of {self.data_file}")


@dataclass(frozen=True)
class DescribedExtrapolation:
    """One ``extrapolation`` record as the engine built it (``describe(as_json=True)["extrapolations"]``): ``key``
    its number among the extrapolations (1-based), the pair keys of its first term, the numbers of grid energies and
    angles, the grid file ("" for a written grid) and whether the grid is adaptive."""
    key: int
    name: str
    active: bool
    quantity: str
    observable: str
    tree: str
    entrance_key: int
    exit_key: int
    terms: tuple
    energies: int
    angles: int
    grid_file: str
    adaptive: bool
    effects: tuple

    @classmethod
    def from_json(cls, d):
        entrance, exit_ = _first_term(d)
        return cls(key=int(d["key"]), name=d["name"], active=bool(d["active"]), quantity=d.get("quantity", ""),
                   observable=d.get("observable", ""), tree=d.get("tree", ""), entrance_key=entrance, exit_key=exit_,
                   terms=tuple(dict(t) for t in d.get("terms", [])), energies=int(d.get("energies", 0)),
                   angles=int(d.get("angles", 0)), grid_file=d.get("grid_file", ""), adaptive=bool(d.get("adaptive")),
                   effects=tuple(d.get("effects", [])))

    def reaction(self, pairs=None):
        """The reaction as written (its tree)."""
        return self.tree


class _Described(list):
    """A list of described records in file order (``key`` 1, 2, ...)."""

    @property
    def active(self):
        """The active ones, in the order the engine evaluates them (results index ``i`` is ``active[i]``)."""
        return type(self)(x for x in self if x.active)

    def by_key(self, key):
        """The record with this key (1, 2, ... in file order); :class:`pyazr.ParameterKeyError` (a ``KeyError``)
        when there is none."""
        for x in self:
            if x.key == key:
                return x
        raise ParameterKeyError(f"no record with key {key!r} (the keys: {[x.key for x in self]}); "
                                f"Session.datasets and Session.extrapolations list them.")

    def by_name(self, name):
        """The record with this name; :class:`pyazr.ParameterKeyError` (a ``KeyError``) when there is none."""
        for x in self:
            if x.name == name:
                return x
        names = [x.name for x in self]
        raise ParameterKeyError(f"no record named {name!r}{close_names(name, names)}; the names: {listed(names)}; "
                                f"Session.datasets and Session.extrapolations list them.")

    def by_reaction(self, entrance_key=None, exit_key=None):
        return type(self)(x for x in self if (entrance_key is None or x.entrance_key == entrance_key)
                          and (exit_key is None or x.exit_key == exit_key))


class DescribedDatasets(_Described):
    """The data sets of a format-3 project as the engine built them (:class:`DescribedDataset` each)."""

    @classmethod
    def from_description(cls, description):
        """From the dict of ``pyazr.describe(path, as_json=True)``."""
        return cls(DescribedDataset.from_json(d) for d in description.get("datasets", []))

    def files(self):
        """The data files, as written, in first-appearance order."""
        return list(dict.fromkeys(d.data_file for d in self))

    def table(self):
        """A text table: key, name, active, quantity, pairs, rows kept / read, file."""
        rows = [f"{'key':>3}  {'name':<20} {'act':<3} {'quantity':<16} {'pairs':<8} {'rows':>11}  file"]
        for d in self:
            rows.append(f"{d.key:>3}  {d.name:<20} {'yes' if d.active else 'no':<3} {d.quantity:<16} "
                        f"{f'{d.entrance_key}->{d.exit_key}':<8} {f'{d.rows_kept}/{d.rows_read}':>11}  {d.data_file}")
        return "\n".join(rows)


class DescribedExtrapolations(_Described):
    """The extrapolations of a format-3 project as the engine built them (:class:`DescribedExtrapolation` each)."""

    @classmethod
    def from_description(cls, description):
        return cls(DescribedExtrapolation.from_json(d) for d in description.get("extrapolations", []))


# -- the package and its engine --------------------------------------------------------------------------------------

def _first_line(text):
    return (text or "").strip().splitlines()[0].strip() if (text or "").strip() else ""


def coulomb_table(use_coulomb_table=True, coulomb_table=""):
    """The Coulomb table a session with these options uses: ``{"in_use": "table" or "exact", "path", "loaded",
    "identity"}``.  ``identity`` is the engine's content hash of the table (``CoulombTable::ContentIdentity``: the
    same text for two files that give the same lookups); "exact" where the option is off or the file cannot be read
    (the engine then uses its exact routines).  Loads the table as a session does (once per process)."""
    if not use_coulomb_table:
        return {"in_use": "exact", "path": "", "loaded": False, "identity": ""}
    info = dict(_engine().coulomb_table_info(str(coulomb_table or "")))
    info["in_use"] = "table" if info["loaded"] else "exact"
    return info


def _file_digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def describe_package():
    """The package and its engine, as plain data (``pyazr.describe()``).

    Keys: ``pyazr`` (version), ``engine`` (version and ``git``, the ``git describe`` of the source at CMake configure
    time; the compiled module's path), ``build`` (mcmc, nlopt, openmp), ``threads`` (the threads an
    evaluation uses now, ``omp_get_max_threads()``, and ``OMP_NUM_THREADS``), ``coulomb_table`` (the default
    session's: :func:`coulomb_table`), ``stopping`` (the stopping library: implementation and citation), ``masses``
    (the mass table: version, file, sha256), ``python`` (version).
    """
    from . import __version__
    from . import stopping
    from .sources import masses
    build = dict(_engine().build_info())
    masses_path = masses._TABLE_PATH
    return {
        "pyazr": __version__,
        "engine": {"version": build.pop("version"), "git": build.pop("git"), "module": _engine().__file__},
        "build": {k: build[k] for k in ("mcmc", "nlopt", "openmp")},
        "threads": {"max_threads": int(build["max_threads"]),
                    "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS")},
        "coulomb_table": coulomb_table(),
        "stopping": {"implementation": stopping.IMPLEMENTATION, "citation": stopping.CITATION},
        "masses": {"version": masses.AME_VERSION, "file": masses_path,
                   "sha256": _file_digest(masses_path) if os.path.isfile(masses_path) else None},
        "python": sys.version.split()[0],
    }


_KINDS = None


def _build_kinds():
    from . import errors
    engine = _engine()
    schema = engine.schema()
    records = {name: {"title": r.get("title", ""), "doc": r.get("doc", ""), "sections": list(r.get("sections", [])),
                      "fields": r.get("fields", {})}
               for name, r in schema["records"].items()}
    data = {}
    for name, q in schema["quantities"].items():
        columns = q.get("columns", {})
        data[name] = {"title": q.get("title", ""), "doc": q.get("doc", ""), "aliases": list(q.get("aliases", [])),
                      "record_kinds": list(q.get("record_kinds", [])), "engine_unit": q.get("engine_unit", ""),
                      "units": list(q.get("units", {})), "frames": list(columns.get("energy", [])),
                      "angle_columns": list(columns.get("angle", [])), "angle_axis": columns.get("angle_axis", ""),
                      "target_kinds": list(q.get("target_kinds", [])),
                      # the numbers a point carries depend on the data set, not on the quantity alone: they are the
                      # rows of "data_outputs" (docs/upgrade/pyazr3.md section 7.1)
                      "extra_outputs": None}
    # the binding's table of point outputs (api/src/pybind_evaluation.cpp): Evaluation.data[i].extra
    data_outputs = {row["name"]: {k: v for k, v in row.items() if k != "name"} for row in engine.point_output_fields()}
    from .evaluate import EvaluatedData, Evaluation
    from .objective import Objective
    objects = {cls.__name__: {"doc": _first_line(cls.__doc__), "built_by": cls.BUILT_BY, "fields": dict(cls.FIELDS)}
               for cls in (Objective, Evaluation, EvaluatedData)}
    parameters = {row["key"]: {"unit": row["unit"], "owner": list(row["owner"]),
                               "own_prior": bool(row["own_prior"]), "segment_field": bool(row["segment_field"]),
                               "vary_from_project": bool(row["vary_from_project"]), "doc": row["doc"]}
                  for row in engine.parameter_kinds()}
    parameter_fields = {f["name"]: {k: (list(v) if isinstance(v, (list, tuple)) else v)
                                    for k, v in f.items() if k != "name"} for f in engine.parameter_record_fields()}
    pair_fields = {f["name"]: {"unit": f["unit"], "doc": f["doc"]} for f in engine.pair_record_fields()}
    options = {}
    for row in engine.session_options():
        options[row["keyword"]] = {"engine_name": row["name"], "type": row["type"], "default": row["default"],
                                   "key": row["key"], "flag": row["flag"], "negation": row["negation"],
                                   "effect": row["effect"], "reason": row["reason"], "doc": row["doc"]}
    # one entry per row of the binding's physics table (api/src/pybind_physics.cpp); pyazr.physics is built from it
    physics = {row["name"]: row for row in engine.physics_functions()}
    # one entry per row of the binding's stopping table (api/src/pybind_stopping.cpp), with SRIM's citation
    stopping = {row["name"]: row for row in engine.stopping.functions()}
    # one entry per row of the engine's fitter table (src/Fitters.cpp), in the program's precedence
    fitters = {row["name"]: {k: v for k, v in row.items() if k != "name"} for row in engine.fitters()}
    error_classes = {}
    for name in errors.__all__:
        cls = getattr(errors, name)
        if isinstance(cls, type) and cls.__name__ == name:
            error_classes[name] = {"bases": [b.__name__ for b in cls.__mro__[1:] if b not in (BaseException, object)],
                                   "doc": _first_line(cls.__doc__)}
    error_classes["EngineError"] = {"same_as": "AZURE2Error", "doc": _first_line(errors.AZURE2Error.__doc__)}
    return {"schema_version": schema.get("schema_version", ""), "records": records, "data": data,
            "parameters": parameters, "parameter_fields": parameter_fields, "pair_fields": pair_fields,
            "options": options, "physics": physics, "stopping": stopping, "fitters": fitters, "units": schema.get("units", {}),
            "dimensions": schema.get("dimensions", {}), "errors": error_classes, "data_outputs": data_outputs,
            "objects": objects, "thm": engine.thm_keys(),
            "nuclides": {n: {"Z": z, "A": a, "mass": m} for n, z, a, m in engine.thm_nuclides()},
            "optical_potentials": {g["name"]: {"projectiles": g["projectiles"], "reference": g["reference"],
                                               "a_range": list(g["a_range"]), "e_range": list(g["e_range"])}
                                   for g in engine.thm_global_opticals()}}


def kinds():
    """What kinds of thing the engine and pyazr know, as plain data (``pyazr.kinds()``).

    Keys:

    - ``records``: the record kinds of project format 3: title, doc, sections and ``fields`` (each with its type,
      grammar, whether required, units or dimension, doc; ``Project.add_<kind>`` takes them);
    - ``data``: the data kinds (quantities of a ``dataset``): title, doc, aliases, record kinds, ``engine_unit``,
      ``units`` accepted, ``frames`` (energy columns), angle columns and axis, target kinds, ``extra_outputs``
      (None: they depend on the data set; see ``data_outputs``);
    - ``data_outputs``: the numbers a point carries beside its model value (``Evaluation.data[i].extra``): unit
      (``"model"``: the data set's engine unit), when (the data sets it applies to), doc, per_point;
    - ``objects``: pyazr's objects of a fit and a calculation (``Objective``, ``Evaluation``, ``EvaluatedData``):
      what builds each and its fields with a line of doc;
    - ``parameters``: the kinds of fit parameter (``Parameter.kind``), in their order in the vector: unit (of the
      value in the fit vector), ``owner`` (the record fields that say what it belongs to), own_prior, segment_field,
      vary_from_project, doc;
    - ``parameter_fields``, ``pair_fields``: the fields of a parameter's and a pair's record (:class:`pyazr.Parameter`,
      :class:`pyazr.system.Pair`): unit (``"unit"``: the parameter's own), doc; a parameter field also the ``kinds`` it
      describes, ``every_kind`` and ``owner`` (shown in ``Parameters.table()`` as what the parameter belongs to);
    - ``options``: the session options by their keyword of :class:`pyazr.Session`: type, default, ``key`` (the
      ``<config>`` key of an option of the calculation, the same word; "" for one that is not project state, with its
      ``reason``), the program's ``flag`` and ``negation``, ``effect`` (physics, convention, numerics), doc;
    - ``physics``: the rows of the binding's physics table, by name (:mod:`pyazr.physics` has one function per row):
      summary, doc, library (the C++ function called), axis (the argument a plot puts on its horizontal axis),
      arguments (name, unit, doc, kind, domain, lower, upper, lower_open, upper_open, choices, default) and returns
      (name, unit, doc);
    - ``physics`` rows also carry ``group`` (their section of the table) and, per argument, ``example`` (a typical
      value) and ``sweep`` (``lower``, ``upper``, ``spacing``: the range of a first plot over it);
    - ``stopping``: the rows of the binding's stopping table, by name, in the shape of the physics rows (a ``text``
      argument ``formula`` first) with ``citation``, SRIM's, owed by every use (:func:`pyazr.stopping.evaluate`);
    - ``fitters``: the engine's fitters by name (``Session.fit_with_engine(method=)``), in the program's precedence:
      option (the program's switch, "" for the default), gradient, residual_rows, bounds, covariance, minos,
      call_limit (honours ``max_calls``; the others ignore it with a notice), fallback (the fitter that runs in its
      place for a model without an analytic Jacobian, "" for none), progress_label, doc;
    - ``errors``: pyazr's error classes with their bases;
    - ``thm``: the keys of the ``<thm>`` section as two record kinds, ``"<thm>"`` (the global options) and
      ``"experiment"`` (the keys of an ``experiment[<name>]`` line): title, doc, ``fields`` by key (type of the schema,
      doc, required, title with the unit, unit, values of an enum, default, grammar of a text, applies_when), read
      from the tables the engine's parser reads (``src/ThmKeys.cpp``); :meth:`pyazr.Project.set_thm` and
      ``set_thm_experiment`` take these keys;
    - ``nuclides``: the built-in nuclei a ``<thm>`` experiment names by name (``beam=18O``), by name: ``Z``, ``A`` and
      the nuclear ``mass`` in u (AME2020); any other nucleus is written ``Z,A,mass``;
    - ``optical_potentials``: the global optical potentials ``opticalAA=`` / ``opticalSF=`` take by name: the
      ``projectiles``, the ``reference``, the target mass numbers ``a_range`` and lab energies ``e_range`` (MeV) of
      the fit (``:extrapolate`` uses one outside them);
    - ``units``, ``dimensions``: the unit table; ``schema_version``.

    Read from the engine's tables at the first call (``AZURE2 --schema``, the parameter-kind table, the binding's
    option table, the binding's physics table, the fitter table); a new row there appears here without an edit of
    Python.  Each call returns a fresh copy.
    """
    global _KINDS
    if _KINDS is None:
        _KINDS = _build_kinds()
    return copy.deepcopy(_KINDS)


def describe_session(session):
    """What one :class:`pyazr.Session` session holds, as plain data (``Session.describe()``); see there."""
    from .session_functions import session_functions
    s = session
    parameters = []
    for p in s.parameters:  # the records, as Session.parameters holds them; an infinite bound is None (JSON)
        entry = p.to_dict()
        entry["bounds"] = [b if np.isfinite(b) else None for b in p.bounds]
        parameters.append(entry)
    pairs = [dict(q.to_dict(), label=q.label, photon=bool(q.is_photon)) for q in s.pairs]
    levels = [{"label": p.level_label, "jgroup": p.jgroup, "J": p.J, "parity": p.parity, "level": p.level,
               "energy": p.level_energy} for p in s.parameters.of_kind("energy")]
    description = describe_project(s.file, s.cwd)
    from .project import Project
    records = Project.open(s.file).records()
    points = {}
    for i, d in enumerate(s.active_datasets if s.mode == "data" else []):
        points[d.name] = int(len(s.energies[i])) if i < len(s.energies) else None
    datasets = []
    for d in description.get("datasets", []):
        entry = {k: d.get(k) for k in ("key", "name", "active", "quantity", "observable", "tree", "unit",
                                         "engine_unit", "data_file", "rows_read", "rows_kept")}
        record = records.get(d.get("name"))
        columns = record.get("columns") if record is not None else None
        entry["columns"] = columns
        entry["frame"] = columns.split()[0].split("[")[0] if isinstance(columns, str) and columns.split() else None
        entry["points"] = points.get(d.get("name"))
        datasets.append(entry)
    extrapolations = [{k: d.get(k) for k in ("key", "name", "active", "quantity", "observable", "tree", "energies",
                                             "angles", "grid_file", "adaptive")}
                      for d in description.get("extrapolations", [])]
    options = dict(s.options)
    options["cwd"] = s.cwd
    return json.loads(json.dumps({
        "project": s.file, "cwd": s.cwd, "output_dir": s.output_dir, "mode": s.mode, "options": options,
        "option_sources": dict(s.option_sources),
        "coulomb_table": coulomb_table(s.options.get("use_coulomb_table", True), s.options.get("coulomb_table", "")),
        "pairs": pairs, "levels": levels, "datasets": datasets, "extrapolations": extrapolations,
        "parameters": parameters, "free_parameters": int(sum(1 for p in parameters if not p["fixed"])),
        "functions": session_functions(),
    }, default=float))
