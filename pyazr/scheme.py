"""A structured, printable view of an AZURE2 R-matrix level scheme.

Turns the flat parameter / pair metadata (``azure2.parameters`` and
``azure2.pairs``) into a small object model grouped the way a physicist reads it
-- particle pairs, then J-groups, then levels, then the channels of each level
with their (L, S, radiation type, partial width, fixed flag, Wigner limit).

This is a *read-only* view of a running instance.  To edit the scheme (add or
remove levels) and write it back to a file, see :mod:`pyazr.azrfile`.
"""

from dataclasses import dataclass, field
from typing import List, Optional


def _fmt_spin(x):
    """Format a spin/parity nicely: 0.5 -> '1/2', 2.0 -> '2'."""
    if x is None:
        return "?"
    if float(x).is_integer():
        return str(int(x))
    # half-integers as p/2
    return f"{int(round(2 * x))}/2"


@dataclass
class SchemeChannel:
    """One channel of a level (a reduced-width / partial-width parameter)."""

    param_index: int                # global parameter index of the width
    pair: int                       # particle-pair number
    L: int
    S: float
    radiation_type: str             # 'P' (particle) or 'E'/'M'/... (photon)
    width: float                    # current physical value (partial width, eV)
    fixed: bool
    wigner_limit: Optional[float]
    free_index: Optional[int]

    @property
    def is_photon(self) -> bool:
        """Is this a photon channel rather than a particle one?"""
        return self.radiation_type not in (None, "P")

    @property
    def label(self) -> str:
        """Compact channel tag, e.g. ``'p1 L=1 S=1/2'`` or ``'E1 (pair 2)'``."""
        if self.is_photon:
            return f"{self.radiation_type}{self.L} (pair {self.pair})"
        return f"pair{self.pair} L={self.L} S={_fmt_spin(self.S)}"


@dataclass
class SchemeLevel:
    """One R-matrix level: a spin-parity, an energy, and its channels."""

    jgroup: int
    level: int
    J: float
    parity: int
    energy: float
    energy_fixed: bool
    energy_free_index: Optional[int]
    channels: List[SchemeChannel] = field(default_factory=list)

    @property
    def jpi(self) -> str:
        """J^pi as text."""
        return f"{_fmt_spin(self.J)}{'+' if self.parity > 0 else '-'}"

    @property
    def total_width(self) -> float:
        """Sum of |partial widths| over the level's channels (eV)."""
        return sum(abs(c.width) for c in self.channels)

    @property
    def n_free(self) -> int:
        """How many of this level's parameters the fit varies."""
        return ((0 if self.energy_fixed else 1)
                + sum(0 if c.fixed else 1 for c in self.channels))


class LevelScheme:
    """The full level scheme of an AZURE2 model: pairs + levels + channels.

    Build one from a running instance with :meth:`from_azr` (or the convenience
    property ``azure2.level_scheme``).  ``str(scheme)`` / ``print(scheme)`` gives
    a grouped, human-readable listing.
    """

    def __init__(self, levels, pairs, name=None):
        self.levels: List[SchemeLevel] = list(levels)
        self.pairs = pairs
        self.name = name

    # -- construction ---------------------------------------------------------

    @classmethod
    def from_azr(cls, azr, name=None):
        """Build the scheme by reading a .azr file."""
        params = azr.parameters
        pairs = azr.pairs

        # index levels by (jgroup, level) preserving first-seen order
        order = []
        by_key = {}
        for p in params:
            if p.kind not in ("energy", "width"):
                continue
            key = (p.jgroup, p.level)
            if key not in by_key:
                by_key[key] = SchemeLevel(
                    jgroup=p.jgroup, level=p.level, J=p.J, parity=p.parity,
                    energy=(p.level_energy if p.level_energy is not None else 0.0),
                    energy_fixed=True, energy_free_index=None)
                order.append(key)
            lvl = by_key[key]
            if p.kind == "energy":
                lvl.energy = p.level_energy if p.level_energy is not None else 0.0
                lvl.energy_fixed = p.fixed
                lvl.energy_free_index = p.free_index
            else:  # width
                lvl.channels.append(SchemeChannel(
                    param_index=p.index, pair=p.pair, L=p.L, S=p.S,
                    radiation_type=p.radiation_type, width=p.value,
                    fixed=p.fixed, wigner_limit=p.wigner_limit,
                    free_index=p.free_index))

        levels = [by_key[k] for k in order]
        return cls(levels, pairs, name=name or getattr(azr, "file", None))

    # -- queries --------------------------------------------------------------

    def by_jpi(self, jpi):
        """The levels of one J^pi group."""
        return [lv for lv in self.levels if lv.jpi == jpi]

    def get(self, jgroup, level):
        """The level at this (jgroup, level) key."""
        for lv in self.levels:
            if lv.jgroup == jgroup and lv.level == level:
                return lv
        raise KeyError((jgroup, level))

    def resonances(self, energy_min=None, energy_max=None):
        """Levels with a physical energy in the given window (excludes far-off
        background poles when ``energy_max`` is set)."""
        out = []
        for lv in self.levels:
            if energy_min is not None and lv.energy < energy_min:
                continue
            if energy_max is not None and lv.energy > energy_max:
                continue
            out.append(lv)
        return out

    # -- rendering ------------------------------------------------------------

    def _pairs_block(self):
        lines = ["Particle pairs:"]
        for p in self.pairs:
            kind = "photon  " if p.is_photon else "particle"
            flags = " [entrance]" if p.is_entrance else ""
            lines.append(
                f"  {p.number}: {kind}  "
                f"J={_fmt_spin(p.J1)}{'+' if p.parity1 > 0 else '-'}"
                f" & {_fmt_spin(p.J2)}{'+' if p.parity2 > 0 else '-'}"
                f"   Z={p.Z1}+{p.Z2}"
                f"   Sep={p.sep_energy:.4g} MeV   r={p.channel_radius:g} fm"
                f"{flags}")
        return lines

    def format(self, widths=True):
        """A grouped, aligned text rendering of the whole scheme."""
        head = f"Level scheme"
        if self.name:
            head += f"  [{self.name}]"
        head += (f"   {len(self.pairs)} pairs, {len(self.levels)} levels, "
                 f"{sum(lv.n_free for lv in self.levels)} free R-matrix params")
        out = [head, "=" * len(head), ""]
        out += self._pairs_block()
        out.append("")

        # group by J^pi in first-seen order
        seen = []
        for lv in self.levels:
            if lv.jpi not in seen:
                seen.append(lv.jpi)
        for jpi in seen:
            group = [lv for lv in self.levels if lv.jpi == jpi]
            out.append(f"J^pi = {jpi}   (J-group {group[0].jgroup})")
            for lv in group:
                efix = "fixed" if lv.energy_fixed else "FREE "
                out.append(f"  level {lv.level}:  E = {lv.energy:>9.4g} MeV  "
                           f"[{efix}]")
                if widths:
                    for c in lv.channels:
                        fix = "fixed" if c.fixed else "FREE "
                        wl = ("" if c.wigner_limit is None
                              else f"  wigner={c.wigner_limit:.3g}")
                        out.append(
                            f"      {c.label:<22}  Gamma = {c.width:>12.5g} eV "
                            f" [{fix}]{wl}")
            out.append("")
        return "\n".join(out).rstrip()

    def __str__(self):
        return self.format()

    def __repr__(self):
        return (f"LevelScheme({len(self.levels)} levels, "
                f"{len(self.pairs)} pairs)")
