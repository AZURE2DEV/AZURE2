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


@dataclass(frozen=True)
class LevelKey:
    """Identity of one physical R-matrix level (a resonance / pole).

    AZURE2 numbers levels *within* a J-group, so ``level`` alone is ambiguous
    (every J-group has a level 1).  The unique identity is ``(jgroup, level)``;
    ``J``/``parity``/``energy`` are carried for readability.  Hashable, so it
    can key a dict of per-level parameter subsets.
    """
    jgroup: int
    J: float
    parity: int
    level: int
    energy: Optional[float]

    @property
    def jpi(self) -> str:
        """The level's J^pi as text, e.g. "3/2-"."""
        j = int(self.J) if float(self.J).is_integer() else f"{int(round(2*self.J))}/2"
        return f"{j}{'+' if self.parity > 0 else '-'}"

    def __str__(self):
        e = "?" if self.energy is None else f"{self.energy:.3f}"
        return f"{self.jpi}#{self.level}@{e}MeV"


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
    _NFIELDS = 15

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
         pair, radtype, fixed, value, segment_key, wigner_limit) = record

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
        if self.kind in ("norm", "shift"):
            bits.append(f"segment={self.segment_key}")
        return "Parameter(" + ", ".join(bits) + ")"


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
    """

    pair: int
    enabled: bool
    type: str
    V0: float
    R: float
    a: float
    r0: float
    own: bool

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
        return (f"NuclearPotential({who}, {shape}"
                f"{'' if self.own or self.pair == 0 else ', inherited'})")


@dataclass
class Pair:
    """One AZURE2 particle pair (a reaction channel's two constituents).

    Built from the numeric records emitted by ``AZUREAPI::GetPairsInfo``
    (command ``GET_PAIRS_INFO``); see ``api/include/AZUREAPI.h`` for the field
    layout, which :meth:`from_record` mirrors.  A width :class:`Parameter`'s
    ``pair`` attribute is this object's :attr:`number`.
    """

    # Number of doubles per record emitted by AZUREAPI::GetPairsInfo.
    # Must match AZUREAPI::kPairInfoFields.
    _NFIELDS = 16

    number: int             # 1-based pair number (matches Parameter.pair)
    key: int                # user pair key from the .azr file
    ptype: int              # 0 = particle channel, otherwise photon
    is_entrance: bool       # the reaction's entrance pair?
    J1: float               # intrinsic spin of particle 1
    parity1: int
    Z1: int
    M1: float               # mass of particle 1 (amu)
    J2: float               # intrinsic spin of particle 2
    parity2: int
    Z2: int
    M2: float
    sep_energy: float       # separation energy (MeV)
    excitation: float       # pair excitation energy (MeV)
    channel_radius: float   # channel radius (fm)
    i1i2factor: float       # 1 / ((2 J1 + 1)(2 J2 + 1))

    @property
    def is_photon(self) -> bool:
        """True for the capture / radiative (gamma) pair."""
        return self.ptype != 0

    @classmethod
    def from_record(cls, record):
        """Build a :class:`Pair` from one ``GET_PAIRS_INFO`` record."""
        (number, key, ptype, is_entrance, J1, parity1, Z1, M1, J2, parity2,
         Z2, M2, sepE, exE, chRad, i1i2factor) = record
        return cls(
            number=int(round(number)),
            key=int(round(key)),
            ptype=int(round(ptype)),
            is_entrance=bool(round(is_entrance)),
            J1=float(J1), parity1=int(round(parity1)), Z1=int(round(Z1)),
            M1=float(M1),
            J2=float(J2), parity2=int(round(parity2)), Z2=int(round(Z2)),
            M2=float(M2),
            sep_energy=float(sepE), excitation=float(exE),
            channel_radius=float(chRad), i1i2factor=float(i1i2factor),
        )

    def __repr__(self):
        return (f"Pair(#{self.number}"
                f"{' entrance' if self.is_entrance else ''}"
                f"{' photon' if self.is_photon else ''}, "
                f"J1={self.J1:g}, J2={self.J2:g}, "
                f"i1i2factor={self.i1i2factor:g})")


class PairSet(list):
    """An ordered collection of :class:`Pair` (1-based number == index+1)."""

    @property
    def entrance(self):
        """The entrance :class:`Pair` (or ``None`` if none is flagged)."""
        for p in self:
            if p.is_entrance:
                return p
        return None

    @property
    def photons(self):
        """All radiative (capture) pairs."""
        return PairSet(p for p in self if p.is_photon)

    @property
    def particles(self):
        """All particle pairs."""
        return PairSet(p for p in self if not p.is_photon)

    def by_number(self, number):
        """The pair with this 1-based number."""
        for p in self:
            if p.number == number:
                return p
        raise KeyError(number)

    def __repr__(self):
        return f"PairSet({len(self)} pairs)"


class ParameterSet(list):
    """An ordered collection of :class:`Parameter` with convenient views.

    Behaves like a ``list`` (indexable, iterable) but adds filtered views and a
    couple of lookups so callers can pull out exactly the parameters they care
    about, e.g. ``params.free``, ``params.widths``, ``params.by_level(2, +1)``.
    """

    # -- the free vector ------------------------------------------------------
    #
    # A calculation takes a flat vector of the *free* parameters, and a
    # Parameter's free_index is its position in it.  Threading that by hand is
    # the usual source of silent misalignment, so these do it.

    def indices(self):
        """Free-vector positions of these parameters, in this set's order.

        Fixed parameters have no position and are dropped, so
        ``len(indices()) <= len(self)``.
        """
        return [p.free_index for p in self
                if not p.fixed and p.free_index is not None]

    def take(self, x):
        """The entries of the free vector ``x`` belonging to these parameters."""
        import numpy as _np
        x = _np.asarray(x, float)
        return x[_np.asarray(self.indices(), int)] if self.indices() else _np.empty(0)

    def put(self, x, values):
        """A *copy* of the free vector ``x`` with these parameters set.

        ``values`` is either one number for all of them or one per parameter,
        in this set's order:

        >>> x = m.parameters.widths.put(m.params_rwa, 0.0)   # zero every width
        """
        import numpy as _np
        out = _np.array(x, float, copy=True)
        idx = _np.asarray(self.indices(), int)
        values = _np.asarray(values, float)
        if values.ndim == 0:
            out[idx] = float(values)
        elif values.size == idx.size:
            out[idx] = values
        else:
            raise ValueError(
                f"put() got {values.size} value(s) for {idx.size} free "
                f"parameter(s) in this set.")
        return out

    # -- filtered views -------------------------------------------------------

    @property
    def free(self):
        """The parameters the fit varies."""
        return ParameterSet(p for p in self if not p.fixed)

    @property
    def fixed(self):
        """The parameters held fixed."""
        return ParameterSet(p for p in self if p.fixed)

    @property
    def energies(self):
        """The level-energy parameters."""
        return ParameterSet(p for p in self if p.kind == "energy")

    @property
    def widths(self):
        """The reduced-width parameters."""
        return ParameterSet(p for p in self if p.kind == "width")

    @property
    def norms(self):
        """The data-normalization parameters."""
        return ParameterSet(p for p in self if p.kind == "norm")

    @property
    def shifts(self):
        """The energy-shift parameters."""
        return ParameterSet(p for p in self if p.kind == "shift")

    # -- lookups --------------------------------------------------------------

    def by_name(self, name):
        """The parameter with this name."""
        for p in self:
            if p.name == name:
                return p
        raise KeyError(name)

    def by_physical_level(self):
        """Group the R-matrix parameters (energy + widths) into physical levels.

        Returns a ``dict`` mapping each :class:`LevelKey` to the
        :class:`ParameterSet` of that level's parameters, ordered by
        ``(jgroup, level)``.  A "physical level" is one resonance / pole -- the
        unit you turn on or off -- as opposed to AZURE2's per-J-group level
        numbering, in which every J-group restarts at level 1.

        Examples
        --------
        >>> for key, ps in azr.parameters.by_physical_level().items():
        ...     print(key, len(ps.widths), "widths")
        """
        groups = {}
        for p in self:
            if p.jgroup is None or p.kind not in ("energy", "width"):
                continue
            k = LevelKey(p.jgroup, p.J, p.parity, p.level, p.level_energy)
            groups.setdefault(k, ParameterSet()).append(p)
        return dict(sorted(groups.items(),
                           key=lambda kv: (kv[0].jgroup, kv[0].level or 0)))

    def by_level(self, jgroup=None, parity=None, level=None):
        """All R-matrix parameters (energy + widths) for a level / J-group."""
        def match(p):
            """The parameters matching every keyword given, e.g. match(jpi="1-", pair=2)."""
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
