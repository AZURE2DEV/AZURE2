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
    def composite(self) -> bool:
        return len(self.terms) > 1

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
