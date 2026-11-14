"""The particle pairs of a session, as the engine's records (``Session.pair_records()``; the fields with their units:
``pyazr.kinds()["pair_fields"]``).  :class:`Pair` is one record, :class:`PairSet` the session's pairs in the
engine's order (``Session.pairs``)."""

from ..errors import ParameterKeyError
from ..parameters import _Record


class Pair(_Record):
    """One particle pair (a reaction channel's two constituents): the engine's record of it
    (``Session.pair_records()``; ``pyazr.kinds()["pair_fields"]`` lists the fields with their units).

    ``number`` is the engine's pair number (:attr:`Parameter.pair`), ``key`` the pair's key in ``<levels>``; the two
    differ when ``<levels>`` does not introduce the pairs in key order.
    """

    __slots__ = ()
    _which = "pair"

    @property
    def heavy_excitation(self) -> float:
        """Excitation energy of particle 2 (MeV): :attr:`excitation` less :attr:`light_excitation`."""
        return self.excitation - (self.light_excitation or 0.0)  # a record made by hand may leave it out

    @property
    def threshold(self) -> float:
        """Excitation energy of the compound nucleus at which the pair's channels open (MeV): ``sep_energy +
        excitation``."""
        return self.sep_energy + self.excitation

    @property
    def is_photon(self) -> bool:
        """True for a pair that is not a particle pair (a photon or a beta-decay pair)."""
        return self.ptype != 0

    @property
    def nucleus(self):
        """``(Z, A)`` of the compound nucleus this pair forms."""
        return (self.nucleus_Z, self.nucleus_A)

    @property
    def label(self) -> str:
        """The two particles' J^pi, an excited one marked with its excitation energy: ``'0+ + 2+*(4.439 MeV)'``."""
        def part(J, parity, ex):
            return f"{J:g}{'+' if parity > 0 else '-'}" + (f"*({ex:.4g} MeV)" if ex else "")
        if self.is_photon:
            return f"gamma + {part(self.J2, self.parity2, self.heavy_excitation)}"
        return (f"{part(self.J1, self.parity1, self.light_excitation)} + "
                f"{part(self.J2, self.parity2, self.heavy_excitation)}")

    def __repr__(self):
        # the engine's number AND the key: engine messages and isospin.out name pairs by key, Parameter.pair by number
        return (f"Pair(#{self.number} key={self.key}{' entrance' if self.is_entrance else ''}"
                f"{' photon' if self.is_photon else ''}, J1={self.J1:g}, J2={self.J2:g}, "
                f"{f'light_excitation={self.light_excitation:g}, ' if self.light_excitation else ''}"
                f"{f'heavy_excitation={self.heavy_excitation:g}, ' if self.heavy_excitation else ''}"
                f"i1i2factor={self.i1i2factor:g})")


class PairSet(list):
    """An ordered collection of :class:`Pair` (1-based number == index+1)."""

    @property
    def entrance(self):
        """The entrance :class:`Pair` (or ``None`` if none is flagged)."""
        return next((p for p in self if p.is_entrance), None)

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
        raise ParameterKeyError(f"no pair {number!r} (the pair numbers: {[p.number for p in self]}); "
                                f"Session.pairs lists the pairs.")

    def __repr__(self):
        return f"PairSet({len(self)} pairs)"
