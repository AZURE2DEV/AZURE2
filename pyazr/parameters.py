"""Structured description of the AZURE2 fit parameters.

The AZURE2 API exposes its parameters only as flat, parallel arrays (names,
values, fixed flags).  This module turns that into a small, self-describing
object model so callers can ask *what* each parameter actually is -- which
level (J^pi, energy) an R-matrix parameter belongs to, and for a width which
channel it is (orbital angular momentum L, channel spin S, particle pair,
radiation type) -- without having to parse the parameter names by hand.

The numeric metadata is produced by ``AZUREAPI::GetParameterInfo`` (command
``GET_PARAMS_INFO``); see ``api/include/AZUREAPI.h`` for the field layout, which
``Parameter.from_record`` mirrors.
"""

from dataclasses import dataclass, fields
from typing import Optional


# Parameter kinds, indexed by the integer ``type`` code in the C++ record.
_KINDS = {0: "energy", 1: "width", 2: "norm", 3: "shift"}


def _opt(value, *, integer=False, sentinel=-1.0):
    """Return ``None`` for the C++ "not applicable" sentinel, else the value."""
    if value == sentinel:
        return None
    return int(round(value)) if integer else float(value)


@dataclass
class Parameter:
    """One AZURE2 fit parameter and everything known about it.

    Attributes that do not apply to a given parameter kind are ``None`` (e.g. a
    normalization has no ``L``; an energy has no ``channel``).
    """

    # Number of doubles per record emitted by AZUREAPI::GetParameterInfo.
    # Must match AZUREAPI::kParamInfoFields.
    _NFIELDS = 16

    index: int                      # position among *all* parameters
    name: str                       # raw AZURE2 parameter name
    kind: str                       # 'energy' | 'width' | 'norm' | 'shift'
    fixed: bool                     # held fixed during the fit?
    value: float                    # current (physical) value
    free_index: Optional[int]       # position among the non-fixed parameters

    # Level / J-group information (R-matrix parameters).
    jgroup: Optional[int] = None
    J: Optional[float] = None
    parity: Optional[int] = None
    level: Optional[int] = None
    level_energy: Optional[float] = None

    # Channel information (width parameters).
    channel: Optional[int] = None
    L: Optional[int] = None
    S: Optional[float] = None
    pair: Optional[int] = None
    radiation_type: Optional[str] = None   # 'P' (particle), 'E', 'M', ...

    # Wigner limit of the channel's reduced width (width parameters only).
    # This is the bound AZURE2 places on the reduced-width amplitude when the
    # Wigner-limit constraint is enabled, in the same units as ``value`` for a
    # width parameter; ``None`` for non-width parameters.
    wigner_limit: Optional[float] = None

    # True if this width's .azr input value was declared as a reduced width
    # amplitude (MeV^(1/2)) and thus AZURE2 should treat it accordingly.
    # ``None`` for non-width parameters.
    input_is_rwa: Optional[bool] = None

    # Data-segment information (norm / shift parameters).
    segment_key: Optional[int] = None

    @property
    def jpi(self) -> Optional[str]:
        """The level's spin-parity as a string, e.g. ``"2+"`` (or ``None``)."""
        if self.J is None or self.parity is None:
            return None
        j = int(self.J) if float(self.J).is_integer() else self.J
        return f"{j}{'+' if self.parity > 0 else '-'}"

    @classmethod
    def from_record(cls, index, name, record, free_index=None):
        """Build a :class:`Parameter` from one numeric record + its name.

        ``record`` is the slice of ``_NFIELDS`` doubles for this parameter as
        returned by ``GET_PARAMS_INFO``.
        """
        (type_code, jgroup, J, parity, level, level_energy, channel, L, S,
         pair, radtype, fixed, value, segment_key, wigner_limit,
         input_is_rwa) = record

        rad = chr(int(round(radtype))) if radtype != -1 else None

        return cls(
            index=index,
            name=name,
            kind=_KINDS.get(int(round(type_code)), "unknown"),
            fixed=bool(round(fixed)),
            value=float(value),
            free_index=free_index,
            jgroup=_opt(jgroup, integer=True),
            J=_opt(J),
            parity=_opt(parity, integer=True, sentinel=0.0),
            level=_opt(level, integer=True),
            level_energy=_opt(level_energy, sentinel=0.0) if type_code in (0, 1) else None,
            channel=_opt(channel, integer=True),
            L=_opt(L, integer=True),
            S=_opt(S),
            pair=_opt(pair, integer=True),
            radiation_type=rad,
            wigner_limit=_opt(wigner_limit),
            input_is_rwa=(bool(round(input_is_rwa))
                          if input_is_rwa != -1 else None),
            segment_key=_opt(segment_key, integer=True),
        )

    def __repr__(self):
        bits = [f"#{self.index}", self.name, self.kind,
                "fixed" if self.fixed else "free", f"value={self.value:.6g}"]
        if self.kind in ("energy", "width"):
            bits.append(f"J^pi={self.jpi}")
            bits.append(f"E={self.level_energy:.4g}MeV")
        if self.kind == "width":
            bits.append(f"L={self.L}")
            bits.append(f"S={self.S}")
            bits.append(f"pair={self.pair}")
            bits.append(f"rad={self.radiation_type}")
            if self.wigner_limit is not None:
                bits.append(f"wigner={self.wigner_limit:.4g}")
            if self.input_is_rwa:
                bits.append("input=RWA")
        if self.kind in ("norm", "shift"):
            bits.append(f"segment={self.segment_key}")
        return "Parameter(" + ", ".join(bits) + ")"


class ParameterSet(list):
    """An ordered collection of :class:`Parameter` with convenient views.

    Behaves like a ``list`` (indexable, iterable) but adds filtered views and a
    couple of lookups so callers can pull out exactly the parameters they care
    about, e.g. ``params.free``, ``params.widths``, ``params.by_level(2, +1)``.
    """

    # -- filtered views -------------------------------------------------------

    @property
    def free(self):
        return ParameterSet(p for p in self if not p.fixed)

    @property
    def fixed(self):
        return ParameterSet(p for p in self if p.fixed)

    @property
    def energies(self):
        return ParameterSet(p for p in self if p.kind == "energy")

    @property
    def widths(self):
        return ParameterSet(p for p in self if p.kind == "width")

    @property
    def norms(self):
        return ParameterSet(p for p in self if p.kind == "norm")

    @property
    def shifts(self):
        return ParameterSet(p for p in self if p.kind == "shift")

    # -- lookups --------------------------------------------------------------

    def by_name(self, name):
        for p in self:
            if p.name == name:
                return p
        raise KeyError(name)

    def by_level(self, jgroup=None, parity=None, level=None):
        """All R-matrix parameters (energy + widths) for a level / J-group."""
        def match(p):
            if p.jgroup is None:
                return False
            if jgroup is not None and p.jgroup != jgroup:
                return False
            if parity is not None and p.parity != parity:
                return False
            if level is not None and p.level != level:
                return False
            return True
        return ParameterSet(p for p in self if match(p))

    # -- conversions ----------------------------------------------------------

    def values(self, free_only=False):
        """The parameter values as a plain list."""
        src = self.free if free_only else self
        return [p.value for p in src]

    def wigner_limits(self, free_only=False):
        """Per-parameter Wigner limits, aligned one-to-one with the set.

        Each entry is the Wigner limit of the corresponding width parameter's
        reduced width (the bound on the reduced-width amplitude), or ``None``
        for parameters that are not widths.  The list lines up index-for-index
        with :meth:`values`, so ``zip(params.values(), params.wigner_limits())``
        pairs every parameter with its limit.
        """
        src = self.free if free_only else self
        return [p.wigner_limit for p in src]

    def to_records(self):
        """Every parameter as a plain ``dict`` (e.g. for a DataFrame)."""
        keys = [f.name for f in fields(Parameter) if not f.name.startswith("_")]
        return [{k: getattr(p, k) for k in keys} for p in self]

    def table(self):
        """A human-readable, aligned text table of the parameters."""
        rows = [("idx", "name", "kind", "free", "J^pi", "E(MeV)",
                 "L", "S", "pair", "rad", "seg", "value", "wigner")]
        for p in self:
            rows.append((
                str(p.index),
                p.name,
                p.kind,
                "" if p.fixed else "*",
                p.jpi or "",
                "" if p.level_energy is None else f"{p.level_energy:.4g}",
                "" if p.L is None else str(p.L),
                "" if p.S is None else f"{p.S:g}",
                "" if p.pair is None else str(p.pair),
                p.radiation_type or "",
                "" if p.segment_key is None else str(p.segment_key),
                f"{p.value:.6g}",
                "" if p.wigner_limit is None else f"{p.wigner_limit:.4g}",
            ))
        widths = [max(len(r[c]) for r in rows) for c in range(len(rows[0]))]
        lines = []
        for r in rows:
            lines.append("  ".join(cell.ljust(widths[c])
                                   for c, cell in enumerate(r)))
        return "\n".join(lines)

    def __repr__(self):
        return f"ParameterSet({len(self)} parameters, {len(self.free)} free)"
