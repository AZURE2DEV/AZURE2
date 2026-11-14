"""The fit parameters and the particle pairs of a session, as the engine's records.

A session hands each fit parameter over as a record with named fields (``Session.parameter_records()``) and each
particle pair likewise (``Session.pair_records()``).  The fields, their units and the kinds of parameter they describe
are one table of the engine (``api/src/pybind_records.cpp``; ``pyazr.kinds()["parameter_fields"]`` and
``["pair_fields"]``): this module types no field name of its own.

- :class:`Parameter`: one record; ``p.name``, ``p.kind``, ``p.unit``, ``p.start``, ``p.fixed``, ``p.free_index``,
  ``p.bounds`` and the fields of its kind (``p.level_label``, ``p.L``, ``p.segment_key``, ``p.prior_sigma``, ...);
  a field of another kind reads ``None``.
- :class:`Parameters`: the table of a session's parameters, one row per parameter in the engine's order.  The fit
  vector ``x`` is the only vector: ``x0`` (a copy of the start), ``names``, ``kinds``, ``units``, ``free_mask``,
  ``bounds``; ``of_kind``, ``of_level``, ``by_name``, ``free`` select rows (a selection is again a
  :class:`Parameters`); ``put``, ``take`` write and read rows of ``x``; ``read``, ``write`` a parameter file;
  ``table()``, ``to_records()``.
- the particle pairs' records: :class:`pyazr.system.Pair` (:mod:`pyazr.system`).
"""

import functools
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .errors import ParameterError, ParameterKeyError, close_names, listed


def _engine():
    from . import _azure2
    return _azure2


def parameter_kinds():
    """The kinds of fit parameter as the engine defines them (``ParameterKind.h``): one dict per kind, in their order in
    the vector, with ``key`` (:attr:`Parameter.kind`), ``unit``, ``owner`` (the record fields that say what it belongs
    to), ``own_prior``, ``segment_field``, ``vary_from_project`` and ``doc``."""
    return _engine().parameter_kinds()


@functools.lru_cache(maxsize=None)
def _field_table(which):
    """The engine's table of record fields, ``"parameter"`` or ``"pair"``, as a tuple of dicts (read once)."""
    rows = _engine().parameter_record_fields() if which == "parameter" else _engine().pair_record_fields()
    return tuple(rows)


def record_fields():
    """The fields of the records, as the engine's table states them: ``{"parameter": [...], "pair": [...]}``, one dict
    per field with ``name``, ``unit`` (``"unit"``: the parameter's own unit), ``doc``; a parameter field also with
    ``kinds`` (the kinds it describes), ``every_kind`` and ``owner`` (shown as what the parameter belongs to)."""
    return {which: [dict(r) for r in _field_table(which)] for which in ("parameter", "pair")}


@functools.lru_cache(maxsize=None)
def _field_names(which):
    return frozenset(r["name"] for r in _field_table(which))


@functools.lru_cache(maxsize=None)
def _owner_fields(kind):
    """The fields a parameter of ``kind`` shows as its owner, in the table's order."""
    return tuple(r["name"] for r in _field_table("parameter") if r["owner"] and kind in r["kinds"])


class _Record:
    """A read-only record of the engine: its fields are attributes; a field of the table that this record does not
    carry (it describes another kind) reads ``None``."""

    __slots__ = ("_record",)
    _which = None  # "parameter" or "pair": the engine's table of fields

    def __init__(self, record=(), **fields):
        """From the engine's record (a dict), and/or its fields as keywords (a record made by hand, as in a test)."""
        object.__setattr__(self, "_record", dict(record, **fields))

    def __getattr__(self, name):
        record = object.__getattribute__(self, "_record")
        if name in record:
            return record[name]
        if name in _field_names(self._which):
            return None
        raise AttributeError(f"{type(self).__name__} has no field {name!r}{close_names(name, record)}; its fields: "
                             f"{listed(sorted(record), 60)} (pyazr.kinds()['{self._which}_fields'] describes them).")

    def __setattr__(self, name, value):
        raise AttributeError(f"{type(self).__name__} is the engine's record and cannot be changed; a value goes into "
                             f"the fit vector (Parameters.put) or the project (pyazr.Project).")

    def __reduce__(self):
        return (type(self), (self._record,))

    def __eq__(self, other):
        return type(other) is type(self) and self._record == other._record

    __hash__ = None

    def __dir__(self):
        return sorted(set(super().__dir__()) | set(self._record))

    def to_dict(self):
        """The record as a plain dict (a copy)."""
        return dict(self._record)


class Parameter(_Record):
    """One fit parameter: the engine's record of it (``Session.parameter_records()``).

    Every parameter has ``index`` (among all), ``name``, ``kind``, ``ordinal``, ``unit``, ``start`` (in the fit
    vector's basis and ``unit``), ``value`` (as the project states it), ``fixed``, ``free_index`` (its entry of the fit
    vector, ``None`` when fixed) and ``bounds``; the fields of its kind say what it belongs to (a level's
    ``level_label``, a channel's ``L``, ``S``, ``pair``, a data set's ``segment_key``, a nuisance's ``prior_centre``
    and ``prior_sigma``, ...).  ``pyazr.kinds()["parameter_fields"]`` lists every field with its unit and doc.
    """

    __slots__ = ()
    _which = "parameter"

    @property
    def nucleus(self):
        """``(Z, A)`` of an abundance parameter's compound nucleus, else ``None``."""
        return None if self.nucleus_Z is None else (self.nucleus_Z, self.nucleus_A)

    @property
    def jpi(self) -> Optional[str]:
        """The J group's spin-parity as text, e.g. ``"2+"``, ``"0.5-"`` (``None`` for a parameter of no J group)."""
        if self.J is None or self.parity is None:
            return None
        j = int(self.J) if float(self.J).is_integer() else self.J
        return f"{j}{'+' if self.parity > 0 else '-'}"

    def owner_text(self):
        """What the parameter belongs to, as the parameter table shows it: its kind's owner fields, ``field=value``."""
        return " ".join(f"{f}={_text(getattr(self, f))}" for f in _owner_fields(self.kind)
                        if getattr(self, f) is not None)

    def __repr__(self):
        return (f"Parameter(#{self.index}, {self.name}, {self.kind}, {'fixed' if self.fixed else 'free'}, "
                f"start={self.start:.6g} {self.unit}" + (f", {self.owner_text()}" if self.owner_text() else "") + ")")


def _text(value):
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, (list, tuple)):
        return ".".join(str(v[1] if isinstance(v, (list, tuple)) else v) for v in value)
    return str(value)


def read_parameter_file(path):
    """``(basis, {name: value})`` of a ``param.par`` / ``param.sav`` / ``param.fit``: ``basis`` is the
    ``#parametrization`` tag (0 standard, 1 Brune, 2 Park), ``None`` without one.

    Blank lines and ``#`` comments are skipped.  A parameter file is read by NAME, never by row position: AZURE2
    matches names (``AZUREParams::ReadUserParameters``), and a positional read of a file of another layout is silently
    off by one or more rows."""
    values, basis = {}, None
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if len(t) >= 2 and t[0] in ("#parametrization", "parametrization"):
                try:
                    basis = int(round(float(t[1])))
                except ValueError:
                    pass
                continue
            if len(t) < 2 or t[0].startswith("#"):
                continue
            try:
                values[t[0]] = float(t[1])
            except ValueError:
                continue
    return basis, values


class Parameters(list):
    """The table of a session's fit parameters: one :class:`Parameter` per row, in the engine's order
    (``Session.parameters``; built once per session and layout).

    The fit vector ``x`` (the free parameters, in the order of their ``free_index``, the engine's basis) is the only
    vector of the interface.  ``x0`` is a copy of its start; every selection carries the same ``x0``, so
    ``sel.put(sel.x0, ...)`` and ``table.put(table.x0, sel, ...)`` are the same.

    >>> P = s.parameters
    >>> w = P.of_level(P.of_kind("energy")[0].level_label).of_kind("width")
    >>> x1 = P.put(P.x0, w, 2 * P.take(P.x0, w))       # a copy of x0 with the level's widths doubled
    >>> print(P.table())
    """

    def __init__(self, rows=(), *, x0=None, basis=None):
        super().__init__(rows)
        self._x0 = None if x0 is None else np.array(x0, float)
        self._basis = basis

    def _select(self, rows):
        return type(self)(rows, x0=self._x0, basis=self._basis)

    def __getitem__(self, i):
        got = super().__getitem__(i)
        return self._select(got) if isinstance(i, slice) else got

    # -- the columns ----------------------------------------------------------------------------------------------

    @property
    def x0(self):
        """A copy of the start of the fit vector (every free parameter's ``start``, in ``free_index`` order)."""
        if self._x0 is None:
            raise ParameterError("this set of parameters was built without its session, so it has no fit vector; "
                                 "Session.parameters is the session's table.")
        return self._x0.copy()

    @property
    def names(self):
        """The parameters' names, one per row."""
        return [p.name for p in self]

    @property
    def kinds(self):
        """The parameters' kinds (keys of ``pyazr.kinds()["parameters"]``), one per row."""
        return [p.kind for p in self]

    @property
    def units(self):
        """The units of ``start``, ``bounds`` and the entries of x, one per row."""
        return [p.unit for p in self]

    @property
    def free_mask(self):
        """``True`` for a row the fit varies (an entry of x), as a numpy array."""
        return np.array([not p.fixed for p in self], bool)

    @property
    def bounds(self):
        """``(lower, upper)`` of each row, ``-inf`` / ``inf`` where the engine sets none (for an optimiser:
        ``P.free.bounds``, in the order of x)."""
        return [tuple(p.bounds) for p in self]

    # -- selections -------------------------------------------------------------------------------------------------

    @property
    def free(self):
        """The rows the fit varies, in the order of x (a :class:`Parameters`)."""
        return self._select(sorted((p for p in self if not p.fixed), key=lambda p: p.free_index))

    def of_kind(self, *kinds):
        """The rows of these kinds (keys of ``pyazr.kinds()["parameters"]``: ``"width"``, ``"norm"``, ...)."""
        known = [k["key"] for k in parameter_kinds()]
        for kind in kinds:
            if kind not in known:
                raise ParameterKeyError(f"no parameter kind {kind!r}{close_names(kind, known)}; the kinds are "
                                        f"{listed(known, 20)} (pyazr.kinds()['parameters']).",
                                        what=f"no parameter kind {kind!r}", allowed=f"the kinds are {listed(known, 20)}",
                                        fix="pyazr.kinds()['parameters']")
        return self._select(p for p in self if p.kind in kinds)

    def of_level(self, level):
        """The rows of one level (its energy, widths and gamma_decay numbers), named by its label as
        :attr:`Parameter.level_label` and ``Project.levels()`` give it (``"3/2- 3.5032"``), by ``(jgroup, level)``
        or by one of its parameters."""
        if isinstance(level, Parameter):
            key = (level.jgroup, level.level)
        elif isinstance(level, tuple):
            key = tuple(level)
        else:
            key = next(((p.jgroup, p.level) for p in self if p.level_label == level), None)
        rows = [p for p in self if p.level is not None and (p.jgroup, p.level) == key]
        if not rows:
            labels = list(dict.fromkeys(p.level_label for p in self if p.level_label))
            raise ParameterKeyError(f"no level {level!r} among these parameters{close_names(level, labels)}; the "
                                    f"levels: {listed(labels, 10)} (Parameter.level_label, Project.levels()).",
                                    what=f"no level {level!r} among these parameters",
                                    allowed=f"the levels: {listed(labels, 10)}",
                                    fix="Parameter.level_label, Project.levels()")
        return self._select(rows)

    def by_name(self, name):
        """The parameter with this name (:attr:`Parameter.name`).

        Raises :class:`pyazr.ParameterKeyError` (a ``KeyError``) naming the close names when none has it."""
        for p in self:
            if p.name == name:
                return p
        names = self.names
        raise ParameterKeyError(f"no parameter named {name!r} among the {len(names)} of this set"
                                f"{close_names(name, names)}; the names are Parameter.name, e.g. "
                                f"{listed(names, 6)}; print(parameters.table()) lists them all.",
                                what=f"no parameter named {name!r}", allowed=f"the names are Parameter.name, e.g. {listed(names, 6)}",
                                fix="print(parameters.table()) lists them all")

    # -- the fit vector ---------------------------------------------------------------------------------------------

    def _rows(self, rows):
        """``rows`` (a Parameter, a name, or a sequence of them) as a list of Parameters, each free."""
        if isinstance(rows, (Parameter, str)):
            rows = [rows]
        out = [self.by_name(r) if isinstance(r, str) else r for r in rows]
        fixed = [p.name for p in out if p.fixed]
        if fixed:
            raise ParameterError(f"{listed(fixed, 6)} {'is' if len(fixed) == 1 else 'are'} fixed, not an entry of the "
                                 f"fit vector; P.free lists the entries (a fixed value is changed in the project).")
        return out

    def _vector(self, x):
        x = np.array(x, float).ravel()
        if self._x0 is not None and x.size != self._x0.size:
            raise ParameterError(f"x has {x.size} entries; the fit vector has {self._x0.size} (P.x0 is one, P.free "
                                 f"names its entries).")
        return x

    def put(self, x, rows, values):
        """A copy of the fit vector ``x`` with rows set: ``rows`` a :class:`Parameter`, a name or a sequence of them
        (all free: ``P.of_kind("width").free``), ``values`` one number for all or one per row."""
        rows = self._rows(rows)
        out = self._vector(x)
        idx = np.array([p.free_index for p in rows], int)
        values = np.asarray(values, float)
        if values.ndim == 0:
            out[idx] = float(values)
        elif values.size == idx.size:
            out[idx] = values.ravel()
        else:
            raise ParameterError(f"put() got {values.size} value(s) for {idx.size} free parameter(s) in this set; give "
                                 f"one number for all, or one per free parameter in this set's order (set.free lists "
                                 f"them).")
        return out

    def take(self, x, rows):
        """The entries of ``x`` of ``rows`` (a :class:`Parameter`, a name or a sequence of them, all free:
        ``P.of_kind("width").free``), as an array."""
        chosen = self._rows(rows)
        x = self._vector(x)
        return x[np.array([p.free_index for p in chosen], int)] if chosen else np.empty(0)

    def read(self, path):
        """The fit vector with the values a parameter file (``param.sav``, ``param.par``) gives its free rows, matched
        by name; a free row the file lacks keeps its start.  The file's fixed values are not applied: a session opened
        with ``parameter_file=path`` takes every value of the file, fixed ones included, as the program does.  Refused:
        a file of another basis, a file with no name of these rows."""
        basis, values = read_parameter_file(path)
        if basis is not None and self._basis is not None and basis != self._basis:
            raise ParameterError(f"{path} holds basis {basis} amplitudes (0 standard, 1 Brune, 2 Park), this session "
                                 f"basis {self._basis}; open the session in that basis, or convert the file.")
        x = self.x0
        hit = 0
        for p in self:
            if p.name in values:
                hit += 1
                if not p.fixed:
                    x[p.free_index] = values[p.name]
        if not hit:
            raise ParameterError(f"{path}: no parameter name matches these parameters (a file of another project?); "
                                 f"the names: {listed(self.names, 6)}.")
        return x

    def write(self, path, x):
        """Write a parameter file (as ``param.sav``): the basis tag, then every row with its value, a free one from
        ``x``, a fixed one its start.  The table must hold every parameter of the session."""
        if [p.index for p in self] != list(range(len(self))):
            raise ParameterError("write needs the session's whole table (Session.parameters), not a selection: a "
                                 "parameter file names every parameter.")
        x = self._vector(x)
        with open(os.fspath(path), "w") as fh:
            # The engine reads the tag and converts a file of the other alternative basis; '#' keeps it a comment for
            # readers that index rows by position.
            fh.write(f"{'#parametrization':>28s} {float(self._basis or 0): .7e} {0.0: .7e}\n")
            for p in self:
                v = p.start if p.fixed else x[p.free_index]
                fh.write(f"{p.name:>28s} {float(v): .7e} {0.0: .7e}\n")
        return os.fspath(path)

    # -- text and records -------------------------------------------------------------------------------------------

    def to_records(self):
        """Every row as a plain ``dict`` of its record's fields (e.g. for a DataFrame)."""
        return [p.to_dict() for p in self]

    def table(self):
        """The rows as aligned text: index, name, kind, unit, start, free ('*'), bounds and what it belongs to."""
        def bound(v):
            return "" if not np.isfinite(v) else f"{v:.4g}"
        rows = [("#", "name", "kind", "unit", "start", "free", "bounds", "owner")]
        for p in self:
            lo, hi = p.bounds
            rows.append((str(p.index), p.name, p.kind, p.unit, f"{p.start:.6g}", "" if p.fixed else "*",
                         "" if not (np.isfinite(lo) or np.isfinite(hi)) else f"[{bound(lo)}, {bound(hi)}]",
                         p.owner_text()))
        widths = [max(len(r[c]) for r in rows) for c in range(len(rows[0]))]
        return "\n".join("  ".join(cell.ljust(widths[c]) for c, cell in enumerate(r)).rstrip() for r in rows)

    def __repr__(self):
        return f"{type(self).__name__}({len(self)} parameters, {sum(1 for p in self if not p.fixed)} free)"


@dataclass
class NuclearPotential:
    """The hybrid Coulomb model's setting for one particle pair.

    A nuclear potential belongs to a pair: it bends the radial wave functions
    of that channel and no other.  ``pair`` is the 1-based :attr:`Pair.number`,
    or ``0`` for the default that every pair without one of its own falls back
    to.

    ``own`` distinguishes a pair that carries its own setting from one merely
    reporting the default -- the values are the same either way, but clearing a
    pair that has none is a no-op.

    ``coulomb_radius`` (fm) is the radius of the uniformly charged sphere of the
    potential's Coulomb term: the ``coulomb_radius`` field of a format-3
    ``potential`` record, 3.5 fm in an older project's ``<potential>`` block.
    Depths are MeV, ``R``, ``a``, ``r0`` fm.  ``matching_radius`` (fm) is where
    the solution inside the potential is matched to the Coulomb functions and
    ``step`` (fm) the radial step of the integration (the record's
    ``matching_radius`` and ``step``; a session reports the values it uses).

    ``spin_orbit`` is ``None`` or the pair's spin-orbit term as a dict
    ``{"depth": MeV, "radius": fm, "diffuseness": fm, "on": "light"|"heavy"}``:
    v_so(r) <2 l.i> with v_so = V_so (2.000 fm^2) (1/r) df/dr, FRESCO's
    convention (``docs/source/theory/hybrid_potential.rst``).
    """

    pair: int
    enabled: bool
    type: str
    V0: float
    R: float
    a: float
    r0: float
    own: bool
    coulomb_radius: float = 3.5
    spin_orbit: Optional[dict] = None
    matching_radius: Optional[float] = None
    step: Optional[float] = None

    @property
    def is_woods_saxon(self):
        """Is the potential a Woods-Saxon rather than a Gaussian?"""
        return self.type == "WoodsSaxon"

    def __repr__(self):
        who = "default" if self.pair == 0 else f"pair {self.pair}"
        if not self.enabled:
            return f"NuclearPotential({who}, off)"
        if self.is_woods_saxon:
            shape = f"WoodsSaxon V0={self.V0:g} R={self.R:g} a={self.a:g}"
        else:
            shape = f"Gaussian V0={self.V0:g} r0={self.r0:g}"
        so = self.spin_orbit
        if so:
            shape += (f" spin-orbit V_so={so['depth']:g} R_so={so['radius']:g} "
                      f"a_so={so['diffuseness']:g} on {so['on']}")
        return (f"NuclearPotential({who}, {shape} R_C={self.coulomb_radius:g}"
                f"{'' if self.own or self.pair == 0 else ', inherited'})")


#: The kinds of parameter a data session and an extrapolation session of one project share, in the same order at the
#: head of the fit vector (the nuisances that follow are the data sets' or the extrapolations').
MODEL_KINDS = ("energy", "width", "abundance")


def carry_model_parameters(x, source, target):
    """A copy of ``target.x0`` (a :class:`Parameters`) with the free model parameters (``MODEL_KINDS``: level
    energies, widths, abundances) taken from the fit vector ``x`` of ``source``, by name: a data session's fit
    vector carried to the extrapolation session of the same project and options.  Raises
    :class:`pyazr.ParameterError` when ``x`` is not a vector of ``source`` and
    :class:`pyazr.ParameterKeyError` when a free model parameter of ``target`` is not one of ``source``."""
    x = source._vector(x)
    rows = [p for p in target.free if p.kind in MODEL_KINDS]
    return target.put(target.x0, rows, [x[source.by_name(p.name).free_index] for p in rows]) if rows else target.x0
