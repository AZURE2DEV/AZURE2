"""The data sets and extrapolations of a format-3 project as the engine built them.

A format-3 project's data sets are records (``dataset``, ``extrapolation``; docs/upgrade/project_format3.md); what the
engine made of them -- the quantity, the observable it computes, the reaction tree with the ``<levels>`` pair keys,
the lab-energy window, the data file and its rows read and kept, the effects and the fit parameters -- is the engine's
description of the built project (``AZURE2 --describe --json``, ``pyazr.describe``).  This module turns that
description into the objects a calculation session hands out (:attr:`pyazr.azure2.datasets`,
:attr:`~pyazr.azure2.active_datasets`, :attr:`~pyazr.azure2.extrapolations`); nothing here reads the file.  The old
formats' segment lines are read by :mod:`pyazr.datasets` instead.
"""

import json
from dataclasses import dataclass
from typing import Optional

__all__ = ["DescribedDataset", "DescribedDatasets", "DescribedExtrapolation", "DescribedExtrapolations"]


def _engine():
    from . import _azure2
    return _azure2


def is_format3(path):
    """Is the project file ``path`` in format 3 (the engine's reading of its ``<config>``)?"""
    with open(path, "rb") as f:
        return _engine().project_format(f.read()) == 3


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
    """What a calculation session hands out as ``kind`` ("datasets" or "extrapolations") for the project ``path``:
    of a format-3 project the engine's description (:class:`DescribedDatasets`, :class:`DescribedExtrapolations`), of
    an older one its segment lines (:class:`pyazr.datasets.SegmentSet`, :class:`~pyazr.datasets.TestSegmentSet`)."""
    from .datasets import SegmentSet, TestSegmentSet
    if not is_format3(path):
        return (SegmentSet if kind == "datasets" else TestSegmentSet).from_file(path)
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
        """The reaction as written (its tree); ``pairs`` is accepted as by :meth:`pyazr.datasets.Segment.reaction`."""
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
        for x in self:
            if x.key == key:
                return x
        raise KeyError(f"no record with key {key}")

    def by_name(self, name):
        for x in self:
            if x.name == name:
                return x
        raise KeyError(f"no record named {name!r}")

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
