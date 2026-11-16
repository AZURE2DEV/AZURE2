"""The project file as rmfit reads it: pyazr 3.0, project format 3.

rmfit's other modules name data sets by their key (the 1-based place among the ``dataset`` records of ``<data>``,
the number of the engine's outputs) and read a handful of facts about each.  This module is the one place those
facts are taken from pyazr: the engine's description of the project (``pyazr.describe(path)``: reaction, windows,
quantity, data file, rows) and the record's own fields (``pyazr.Project``: the texts of ``norm``, ``shift``,
``covariance``).  Nothing here parses a line of the file.

The values, priors and limits of normalizations, energy shifts and the other nuisances are not read from the file
at all: they are fields of the session's parameter rows (``prior_centre``, ``prior_sigma``, ``fixed``,
``segment_key``), which the worker hands to the parent (:func:`nuisance_rows`).
"""
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from . import AZURE2_ROOT  # noqa: F401


@dataclass
class Segment:
    """One ``dataset`` record.  ``key``: its number among the data sets of the file (inactive ones counted);
    ``name``: the record's name; ``observable``: the engine's quantity (``sigma``, ``dsigma_domega``, ...);
    ``entrance_key``/``exit_key``: the pair keys of the reaction's first term (-1: every capture pair);
    ``composite``: the reaction is a sum or a ratio (more than one term); ``covariance``: the record's covariance
    file or ""; ``norm_text``/``shift_text``: the fields as written (a nuisance, or the name of a shared
    ``normalization`` / ``energy_scale`` record, which is then ``norm_group`` / ``shift_group``)."""
    key: int
    name: str
    active: bool
    entrance_key: int
    exit_key: int
    energy_min: float
    energy_max: float
    angle_min: Optional[float]
    angle_max: Optional[float]
    observable: str
    unit: str
    tree: str
    terms: tuple
    data_file: str
    data_file_resolved: str
    rows_read: int
    rows_kept: int
    covariance: str = ""
    norm_text: str = ""
    shift_text: str = ""
    norm_group: str = ""
    shift_group: str = ""
    fields: Tuple[Tuple[str, str], ...] = ()
    # the nuisances, from the session's parameter rows (SegmentTable.with_nuisances); the defaults of a table read
    # without a session say "not known": norm 1 fixed, no shift
    norm: float = 1.0                     # the start value (the prior's centre when it has a prior)
    vary_norm: bool = False
    norm_error: float = 0.0               # the prior's sigma in per cent of its centre (0: no prior)
    energy_shift: float = 0.0             # MeV
    vary_shift: bool = False
    energy_shift_error: float = 0.0
    energy_shift_sqrt: float = 0.0        # MeV^1/2
    vary_shift_sqrt: bool = False
    energy_shift_sqrt_error: float = 0.0
    shift_fn_params: tuple = ()           # (a_j start, sigma, vary) per shift-function parameter

    @property
    def differential(self) -> bool:
        """The quantity is measured at an angle (the engine gives such a data set an angle window)."""
        return self.angle_min is not None

    @property
    def composite(self) -> bool:
        return len(self.terms) > 1

    @property
    def operation(self) -> Optional[str]:
        """None for one reaction; ``"ratio"`` or ``"sum"`` for a reaction of several terms (the engine's tree)."""
        if len(self.terms) <= 1:
            return None
        return "ratio" if "/" in self.tree else "sum"

    @property
    def components(self) -> tuple:
        return self.terms if len(self.terms) > 1 else ()

    @property
    def folded(self) -> bool:
        """The record names a target, a beam or a detector: its points are folded by the engine."""
        return any(self.get(f) for f in ("target", "beam", "detector"))

    def get(self, name, default=""):
        """The text of a field of the record (the first, of a repeatable one), or ``default``."""
        for k, v in self.fields:
            if k == name:
                return v
        return default


class SegmentTable(list):
    """The data sets of a project in file order."""

    def by_key(self, key) -> Optional[Segment]:
        for s in self:
            if s.key == int(key):
                return s
        return None

    def by_name(self, name) -> Optional[Segment]:
        for s in self:
            if s.name == name:
                return s
        return None

    @property
    def active(self) -> "SegmentTable":
        return SegmentTable(s for s in self if s.active)

    def with_nuisances(self, rows) -> "SegmentTable":
        """Fill the nuisance fields of every data set from the session's parameter rows (:func:`nuisance_rows`).
        A member of a shared normalization or energy scale other than the first has a fixed row of its own (the
        group's parameter carries the first member's key): its ``vary_*`` is False, as the engine fits it."""
        fn: Dict[int, dict] = {}
        for r in rows:
            d = self.by_key(r["segment_key"]) if r.get("segment_key") is not None else None
            if d is None:
                continue
            centre, sigma = r.get("prior_centre"), r.get("prior_sigma") or 0.0
            free = not r["fixed"]
            if r["kind"] == "norm":
                d.norm = float(centre if (sigma > 0 and centre is not None) else r["start"])
                d.vary_norm = free
                d.norm_error = 100.0 * sigma / abs(d.norm) if (sigma > 0 and d.norm) else 0.0
            elif r["kind"] == "shift":
                d.energy_shift, d.vary_shift, d.energy_shift_error = float(r["start"]), free, float(sigma)
            elif r["kind"] == "shift_sqrt":
                d.energy_shift_sqrt, d.vary_shift_sqrt, d.energy_shift_sqrt_error = float(r["start"]), free, float(sigma)
            elif r["kind"] == "shift_fn":
                fn.setdefault(d.key, {})[int(r["fn_term"] or 0)] = (float(r["start"]), float(sigma), free)
        for key, terms in fn.items():
            self.by_key(key).shift_fn_params = tuple(terms[j] for j in sorted(terms))
        return self

    @classmethod
    def from_file(cls, path, cwd=None) -> "SegmentTable":
        import pyazr
        path = os.path.abspath(path)
        project = pyazr.Project.open(path)
        described = pyazr.DescribedDatasets.from_description(
            pyazr.described.describe_project(path, cwd or os.path.dirname(path)))
        kinds = {name: rec.record_kind for name, rec in project.records().items()}
        out = cls()
        for d in described:
            rec = project[d.name]
            flds = tuple((str(k), str(v)) for k, v in rec.items())
            get = dict(reversed(flds)).get
            norm, shift = get("norm", ""), get("shift", "")
            out.append(Segment(
                key=d.key, name=d.name, active=d.active, entrance_key=d.entrance_key, exit_key=d.exit_key,
                energy_min=d.energy_min, energy_max=d.energy_max, angle_min=d.angle_min, angle_max=d.angle_max,
                observable=d.quantity, unit=d.unit, tree=d.tree, terms=d.terms, data_file=d.data_file,
                data_file_resolved=d.data_file_resolved, rows_read=d.rows_read, rows_kept=d.rows_kept,
                covariance=get("covariance", ""), norm_text=norm, shift_text=shift,
                norm_group=norm if kinds.get(norm) == "normalization" else "",
                shift_group=shift if kinds.get(shift) == "energy_scale" else "", fields=flds))
            # which nuisances are fitted, as the engine describes the record (no session needed): a member of a
            # shared record other than the first reads "fixed (the group's is its first member's)"
            s, fn = out[-1], []
            for q in d.parameters:
                fit = str(q.get("status", "")) == "fit"
                if q.get("kind") == "norm":
                    s.vary_norm = fit
                elif q.get("kind") == "energy shift":
                    s.vary_shift = fit
                elif q.get("kind") == "energy shift sqrt":
                    s.vary_shift_sqrt = fit
                elif q.get("kind") == "energy shift function":
                    fn.append((0.0, 0.0, fit))
            s.shift_fn_params = tuple(fn)
        return out


def segments(path, cwd=None) -> SegmentTable:
    return SegmentTable.from_file(path, cwd)


def record_kinds(path) -> set:
    """The kinds of the records of the project (``dataset``, ``target``, ``datum``, ...)."""
    import pyazr
    return {rec.record_kind for rec in pyazr.Project.open(os.path.abspath(path)).records().values()}


def has_field(path, name) -> bool:
    """True when a ``dataset`` record of the project has the field ``name`` (``covariance``, ...)."""
    import pyazr
    project = pyazr.Project.open(os.path.abspath(path))
    return any(rec.record_kind == "dataset" and name in rec.keys() for rec in project.records().values())


def without_fields(path, out_path, names) -> str:
    """The project written to ``out_path`` with the fields ``names`` removed from every ``dataset`` record."""
    import pyazr
    project = pyazr.Project.open(os.path.abspath(path))
    for rname, rec in project.records().items():
        if rec.record_kind != "dataset":
            continue
        drop = {n: None for n in names if n in rec.keys()}
        if drop:
            project.set(rname, **drop)
    return project.save(out_path)


def nuisance_rows(parameters) -> List[dict]:
    """The session's parameters that are not level parameters, as plain data for the pipe: ``name``, ``kind``,
    ``fixed``, ``free_index``, ``start``, ``segment_key``, ``fn_term``, ``prior_centre``, ``prior_sigma``,
    ``bounds``.  Level parameters with a prior of ``<parameterSettings>`` are included too."""
    out = []
    for p in parameters:
        sigma = getattr(p, "prior_sigma", None)
        if p.kind in ("energy", "width") and not (sigma and sigma > 0):
            continue
        out.append(dict(name=p.name, kind=p.kind, fixed=bool(p.fixed), free_index=p.free_index,
                        start=float(p.start), segment_key=p.segment_key, fn_term=getattr(p, "fn_term", None),
                        prior_centre=getattr(p, "prior_centre", None), prior_sigma=sigma,
                        bounds=tuple(p.bounds) if p.bounds is not None else None))
    return out


#: the factor to MeV of the energy units a data file's energy column may be written in
_ENERGY_UNITS = {"MeV": 1.0, "keV": 1e-3, "eV": 1e-6}


def read_points(seg: Segment, base_dir: str):
    """``(E_lab in MeV, angle in deg or NaN)`` of every data row of the file of ``seg``, in file order, or None
    when rmfit cannot place the rows itself: the record's ``columns`` must name a lab-energy column
    (``E_lab[MeV|keV|eV]``) and, for a quantity at an angle, an angle column (``theta_lab[deg]`` or
    ``theta_cm[deg]``); a file that names its own columns (``columns from_file``) or another frame gives None, and
    such a data set is then never split.  The windows of the record are not applied here."""
    import numpy as np
    roles = seg.get("columns").split()
    ie = next((i for i, r in enumerate(roles) if r.startswith("E_lab[")), None)
    if ie is None or roles[ie][6:-1] not in _ENERGY_UNITS or not seg.get("energy_window", "x MeV").endswith("MeV"):
        return None
    ia = next((i for i, r in enumerate(roles) if r.startswith("theta_") and r.endswith("[deg]")), None)
    if seg.differential and ia is None:
        return None
    scale = _ENERGY_UNITS[roles[ie][6:-1]]
    rows = []
    path = seg.data_file_resolved or os.path.join(base_dir, seg.data_file)
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s[0] in "#!":
                continue
            t = s.split()
            try:
                rows.append((float(t[ie]) * scale, float(t[ia]) if ia is not None else float("nan")))
            except (ValueError, IndexError):
                continue
    return np.array(rows, float).reshape(-1, 2)


def data_paths(path) -> List[str]:
    """Every file a record of the project names relative to the project (data files, covariance files, grid
    files) and the files of ``<config>`` (the prior covariance), as written."""
    import pyazr
    project = pyazr.Project.open(os.path.abspath(path))
    out = []
    for rec in project.records().values():
        for k, v in rec.items():
            if k in ("file", "covariance", "grid_file"):
                out.append(str(v).split()[0].strip('"'))
    v = project.config().get("priorCovariance")
    if v:
        out.append(str(v))
    return [p for p in dict.fromkeys(out) if p and not os.path.isabs(p)]


def write_variant(base_azr: str, out_path: str, output_dir: str, datasets=None, drop_fields=()) -> Tuple[str, Dict[str, int]]:
    """Write a variant of project ``base_azr`` to ``out_path`` with the output directory ``output_dir``.

    ``datasets``: None keeps the data sets as they are; else a list of ``(base record name, new name, {field:
    value})``: the data sets of the variant, each the base record with the given fields replaced (None removes a
    field); a base record listed more than once is copied under the new names (pieces of a split data set).  Data
    sets not listed are removed, and with them (``datasets`` given) the extrapolations, the model data and the
    prior covariance of ``<config>``, which belong to the whole model.  ``drop_fields``: fields removed from every
    data set kept.  Returns ``(path, {record name: its key in the variant})``."""
    import pyazr
    project = pyazr.Project.open(os.path.abspath(base_azr))
    if datasets is not None:
        used = {b for b, _, _ in datasets}
        for name, rec in list(project.records().items()):
            kind = rec.record_kind
            if kind in ("extrapolation", "datum") or (kind == "dataset" and name not in used):
                project.remove(name)
        if project.config().get("priorCovariance"):
            project.set_config(priorCovariance=None)
        first = {}
        for base, new, fields in datasets:
            fields = dict(fields)
            for f in drop_fields:
                fields.setdefault(f, None)
            if base not in first:
                first[base] = new
                original = dict(project[base].items())
                project.set(base, **{k: v for k, v in fields.items() if not (v is None and k not in original)})
                if new != base:
                    project.rename(base, new)
                first[base] = (new, original)
            else:
                full = dict(first[base][1])
                full.update(fields)
                project.add_dataset(new, **{k: v for k, v in full.items() if v is not None})
    elif drop_fields:
        for name, rec in list(project.records().items()):
            if rec.record_kind == "dataset":
                drop = {f: None for f in drop_fields if f in rec.keys()}
                if drop:
                    project.set(name, **drop)
    project.set_config(output_dir=output_dir.rstrip("/") + "/")
    path = project.save(out_path)
    keys, n = {}, 0
    for name, rec in project.records().items():
        if rec.record_kind == "dataset":
            n += 1
            keys[name] = n
    return path, keys
