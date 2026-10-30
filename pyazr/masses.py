"""Atomic masses and separation energies from the AME2020 mass table.

A particle pair carries both its two masses and its separation energy, and the
two are linked: ``S = m1 + m2 - M`` of the compound nucleus.  Filling both from
one table keeps them consistent; AZURE2 warns at load time when a project's
pairs disagree by more than 1 keV (``CNuc::CheckMassSeparationConsistency``).

The table is ``pyazr/data/mass_1.mas20``, the AME2020 file as the Atomic Mass
Data Center distributes it.  The GUI embeds the same file, so both front ends
fill identical numbers.  Masses are **atomic** masses in u (electrons
included), which is the convention of AZURE2 projects; the separation energy is
computed from the mass excesses, ``S = D1 + D2 - Dc``, where the electrons
cancel because ``Z1 + Z2 = Zc``.

    >>> from pyazr import masses
    >>> masses.atomic_mass("12C")
    12.0
    >>> round(masses.separation_energy("p", "12C"), 6)    # S_p(13N), MeV
    1.94349
    >>> masses.pair_masses("a", "12C")
    PairMasses(Z1=2, M1=4.00260325413, Z2=6, M2=12.0, S=7.16191804, compound=(8, 16))

Particles are named by isotope (``"12C"``, ``"C12"``, ``"c-12"``, ``"4He"``) or
by the usual shorthands: ``n``, ``p``, ``d``, ``t``, ``a``/``alpha``, ``helion``.
``"h"`` is not accepted, since it reads as hydrogen as easily as the helion.  ``(Z, A)`` tuples work everywhere a name does.
"""

import os
import re
from typing import NamedTuple, Optional, Tuple

__all__ = ["AME_VERSION", "MassEntry", "PairMasses", "lookup", "parse_nuclide",
           "atomic_mass", "mass_excess", "separation_energy", "pair_masses",
           "identify", "U_MEV"]

#: The evaluation the table comes from; shown wherever a filled value is reported.
AME_VERSION = "AME2020"

#: The atomic mass unit in MeV used by the AME2020 tables (931494.10242 keV).
U_MEV = 931.49410242

_TABLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "mass_1.mas20")

# Shorthands -> (Z, A).  No "h": it reads as hydrogen as easily as the helion
# (write "3He" or "helion").  The photon and the electron have no entry: they are
# not nuclides and their pair columns are fixed conventions, not table values.
_SHORTHAND = {
    "n": (0, 1), "neutron": (0, 1),
    "p": (1, 1), "proton": (1, 1),
    "d": (1, 2), "deuteron": (1, 2),
    "t": (1, 3), "triton": (1, 3),
    "helion": (2, 3),
    "a": (2, 4), "alpha": (2, 4),
}


class MassEntry(NamedTuple):
    """One row of the mass table."""
    Z: int
    A: int
    symbol: str
    mass_excess: float          #: keV
    mass_excess_unc: float      #: keV
    atomic_mass: float          #: u (atomic, electrons included)
    atomic_mass_unc: float      #: u
    estimated: bool             #: the evaluators' extrapolation ('#' in the file), not a measurement

    @property
    def label(self) -> str:
        """``"12C"``; the neutron is ``"n"``."""
        return "n" if self.Z == 0 else f"{self.A}{self.symbol}"


class PairMasses(NamedTuple):
    """What a particle pair's columns should hold, from one table."""
    Z1: int
    M1: float                   #: u
    Z2: int
    M2: float                   #: u
    S: float                    #: MeV, separation energy of the pair from the compound nucleus
    compound: Tuple[int, int]   #: (Z, A) of the compound nucleus


_TABLE = None


def _number(text):
    """A fixed-width AME field: '#' marks an estimated value in place of the decimal point."""
    t = text.strip()
    if not t or t == "*":
        return None, False
    estimated = "#" in t
    return float(t.replace("#", ".")), estimated


def _load():
    """Parse mass_1.mas20 once (format a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,...,i3,1x,f13.6,f12.6)."""
    global _TABLE
    if _TABLE is not None:
        return _TABLE
    table = {}
    with open(_TABLE_PATH, "r", encoding="ascii") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if len(line) < 123:
                continue
            try:
                Z = int(line[9:14])
                A = int(line[14:19])
            except ValueError:
                continue
            symbol = line[20:23].strip()
            if not symbol.isalpha():
                continue
            dm, est1 = _number(line[28:42])
            ddm, _ = _number(line[42:54])
            try:
                whole = int(line[106:109])
            except ValueError:
                continue
            micro, est2 = _number(line[110:123])
            dmicro, _ = _number(line[123:135])
            if dm is None or micro is None:
                continue
            table[(Z, A)] = MassEntry(
                Z=Z, A=A, symbol=symbol.capitalize() if Z > 0 else "n",
                mass_excess=dm, mass_excess_unc=ddm or 0.0,
                atomic_mass=round(whole + micro * 1e-6, 12),   # the file's 1e-12 u resolution, no float noise
                atomic_mass_unc=(dmicro or 0.0) * 1e-6,
                estimated=est1 or est2)
    if len(table) < 3000:
        raise RuntimeError(f"{_TABLE_PATH}: read only {len(table)} nuclides; is the file complete?")
    _TABLE = table
    return table


def _symbols():
    try:
        from .nuclei import SYMBOLS
    except ImportError:                 # loaded standalone by path, no parent package
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "pyazr_nuclei_standalone", os.path.join(os.path.dirname(os.path.abspath(__file__)), "nuclei.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        SYMBOLS = mod.SYMBOLS
    return SYMBOLS


def parse_nuclide(spec) -> Tuple[int, int]:
    """``(Z, A)`` of ``spec``: a shorthand (``"p"``, ``"a"``), an isotope name
    (``"12C"``, ``"C12"``, ``"c-12"``, ``"4He"``), a ``(Z, A)`` pair or a
    :class:`pyazr.nuclei.Nucleus`."""
    if isinstance(spec, (tuple, list)) and len(spec) == 2:
        return int(spec[0]), int(spec[1])
    s = str(spec).strip()
    if s.lower() in _SHORTHAND:
        return _SHORTHAND[s.lower()]
    m = re.fullmatch(r"\s*(\d+)\s*-?\s*([A-Za-z]{1,3})\s*", s) or None
    if m:
        a_text, sym = m.group(1), m.group(2)
    else:
        m = re.fullmatch(r"\s*([A-Za-z]{1,3})\s*-?\s*(\d+)\s*", s)
        if not m:
            raise ValueError(f"cannot read a nuclide from {spec!r}; give e.g. '12C', 'C12', 'p', 'a' or (Z, A)")
        sym, a_text = m.group(1), m.group(2)
    sym = sym[0].upper() + sym[1:].lower()
    symbols = _symbols()
    if sym == "N" and int(a_text) == 1:
        return 0, 1                     # "1n" / "n1" is the neutron, not nitrogen-1
    if sym not in symbols:
        raise ValueError(f"unknown element symbol {sym!r} in {spec!r}")
    return symbols.index(sym), int(a_text)


def lookup(spec) -> MassEntry:
    """The table row of a nuclide (see :func:`parse_nuclide` for ``spec``)."""
    Z, A = parse_nuclide(spec)
    entry = _load().get((Z, A))
    if entry is None:
        raise KeyError(f"{AME_VERSION} has no nuclide Z={Z}, A={A}")
    return entry


def atomic_mass(spec) -> float:
    """Atomic mass in u."""
    return lookup(spec).atomic_mass


def mass_excess(spec) -> float:
    """Mass excess in keV."""
    return lookup(spec).mass_excess


def separation_energy(light, heavy) -> float:
    """Separation energy of ``light`` from the compound nucleus ``light + heavy``, in MeV.

    ``S = D(light) + D(heavy) - D(compound)`` from the mass excesses; the
    compound nucleus is ``(Z1 + Z2, A1 + A2)``.  This is the number a pair's
    separation-energy column holds, for any excitation energy of the residual
    nucleus (the threshold of an excited-state pair is ``S + E_x``).
    """
    z1, a1 = parse_nuclide(light)
    z2, a2 = parse_nuclide(heavy)
    compound = lookup((z1 + z2, a1 + a2))
    # mass excesses are given to 1e-6 keV, so S to 1e-9 MeV: round off the float noise
    return round((lookup((z1, a1)).mass_excess + lookup((z2, a2)).mass_excess - compound.mass_excess) / 1000.0, 9)


def pair_masses(light, heavy) -> PairMasses:
    """Charges, atomic masses and separation energy of a particle pair."""
    l, h = lookup(light), lookup(heavy)
    return PairMasses(Z1=l.Z, M1=l.atomic_mass, Z2=h.Z, M2=h.atomic_mass,
                      S=separation_energy((l.Z, l.A), (h.Z, h.A)),
                      compound=(l.Z + h.Z, l.A + h.A))


def identify(Z, M) -> Optional[MassEntry]:
    """The nuclide a pair column describes, by its charge and the mass rounded
    to a mass number (the engine's rule); ``None`` for a photon, an electron
    or a nuclide the table does not have."""
    A = int(round(float(M)))
    if A <= 0 or int(Z) < 0:
        return None
    return _load().get((int(Z), A))
