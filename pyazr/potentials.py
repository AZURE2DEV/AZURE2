"""The hybrid nuclear potentials of a session, per particle pair, as the engine resolved them from the project's
``potential`` records (the methods of :class:`pyazr.Session` that read them; ``docs/source/user_guide/pyazr.rst``, "The
hybrid nuclear potential, per particle pair").  A potential is changed by an edit of the :class:`pyazr.Project` (its
``potential`` record) and a new session."""

from .parameters import NuclearPotential
from .locking import engine_call


class _HybridPotentials:
    """Mixed into :class:`pyazr.Session`: uses its ``_engine`` and ``pairs``."""


    @engine_call
    def nuclear_potential(self, pair=0):
        """The hybrid model's setting for one particle pair.

        Returns a :class:`NuclearPotential` -- ``enabled``, ``type``
        (``"WoodsSaxon"`` or ``"Gaussian"``), the shape parameters ``V0``/``R``/
        ``a``/``r0``, ``coulomb_radius`` (fm), and ``own``, which says whether
        this pair carries a setting of its own or is following the default.

        ``pair=0`` addresses that default: the setting every pair falls back to
        when it has none.  A model where no pair is named therefore behaves as
        it did when the potential was a single global object.
        """
        pair = self._check_potential_pair(pair)
        enabled, type_, V0, R, a, r0, own, rc, so, rm, step = self._engine.get_potential(pair)
        spin_orbit = None
        if len(so):
            spin_orbit = {"depth": float(so[0]), "radius": float(so[1]), "diffuseness": float(so[2]),
                          "on": "heavy" if so[3] else "light"}
        return NuclearPotential(pair=pair, enabled=bool(enabled), type=str(type_),
                                V0=float(V0), R=float(R), a=float(a),
                                r0=float(r0), own=bool(own), coulomb_radius=float(rc),
                                spin_orbit=spin_orbit, matching_radius=float(rm), step=float(step))

    @engine_call
    def nuclear_potentials(self):
        """Every pair's resolved setting, plus the default under key ``0``."""
        out = {0: self.nuclear_potential(0)}
        for p in self.pairs:
            out[p.number] = self.nuclear_potential(p.number)
        return out

    def _check_potential_pair(self, pair):
        pair = int(pair)
        if pair and not any(p.number == pair for p in self.pairs):
            raise KeyError(f"no pair {pair} in this model "
                           f"(have {[p.number for p in self.pairs]}); "
                           f"pair=0 is the default.")
        if pair < 0:
            raise ValueError(f"pair must be >= 0, got {pair}.")
        return pair
