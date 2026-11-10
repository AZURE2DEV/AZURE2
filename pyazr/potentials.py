"""The hybrid nuclear potentials of a session, per particle pair (the methods of :class:`pyazr.azure2` that read and
set them; ``docs/source/user_guide/pyazr.rst``, "The hybrid nuclear potential, per particle pair")."""
from .parameters import NuclearPotential


class _HybridPotentials:
    """Mixed into :class:`pyazr.azure2`: uses its ``sess``, ``pairs`` and ``configure``."""


    def nuclear_potential(self, pair=0):
        """The hybrid model's setting for one particle pair.

        Returns a :class:`NuclearPotential` -- ``enabled``, ``type``
        (``"WoodsSaxon"`` or ``"Gaussian"``), the shape parameters ``V0``/``R``/
        ``a``/``r0``, and ``own``, which says whether this pair carries a
        setting of its own or is following the default.

        ``pair=0`` addresses that default: the setting every pair falls back to
        when it has none.  A model where no pair is named therefore behaves as
        it did when the potential was a single global object.
        """
        pair = self._check_potential_pair(pair)
        enabled, type_, V0, R, a, r0, own = self.sess.get_potential(pair)
        return NuclearPotential(pair=pair, enabled=bool(enabled), type=str(type_),
                                V0=float(V0), R=float(R), a=float(a),
                                r0=float(r0), own=bool(own))

    def nuclear_potentials(self):
        """Every pair's resolved setting, plus the default under key ``0``."""
        out = {0: self.nuclear_potential(0)}
        for p in self.pairs:
            out[p.number] = self.nuclear_potential(p.number)
        return out

    def set_nuclear_potential(self, pair=0, type=None, enabled=True,
                              V0=None, R=None, a=None, r0=None,
                              reinitialize=True):
        """Give one particle pair its own hybrid nuclear potential.

        A nuclear potential belongs to a pair -- it bends the radial wave
        functions of that channel and no other -- so each pair carries its own,
        and ``pair=0`` sets the default that unnamed pairs inherit.  Anything
        left at ``None`` keeps the value the pair already resolves to, so
        turning one pair off is just
        ``set_nuclear_potential(2, enabled=False)``.

        ``type`` is ``"WoodsSaxon"`` (``V0``, ``R``, ``a``) or ``"Gaussian"``
        (``V0``, ``r0``); depths are MeV and lengths fm.

        The potential is read by ``CoulFunc`` when it is constructed, deep
        inside a calculation, so a change only reaches the model once the model
        is rebuilt.  That is what ``reinitialize=True`` does -- it costs about
        as much as opening the session did, so pass ``reinitialize=False`` when
        setting several pairs and let the last call do it.

        Like :meth:`set_channel_radius` this changes only the running session,
        and it changes the physics: a fit made without the potential is not a
        fit made with it.

        **Re-read the parameter vector afterwards.**  The potential changes the
        penetrabilities and shift functions, and those are what map physical
        widths to reduced-width amplitudes -- so ``params_rwa`` is re-derived by
        the rebuild and a vector captured beforehand no longer describes the
        same model.  Feeding the old one back gives a quietly wrong chi-squared;
        re-reading ``m.params_rwa`` after the call reproduces exactly what
        loading a ``.azr`` carrying the same potential gives.
        """
        pair = self._check_potential_pair(pair)
        cur = self.nuclear_potential(pair)
        type_ = cur.type if type is None else str(type)
        if type_ not in ("WoodsSaxon", "Gaussian"):
            raise ValueError(f"unknown potential type {type_!r}; expected "
                             f"'WoodsSaxon' or 'Gaussian'.")
        V0 = cur.V0 if V0 is None else float(V0)
        R = cur.R if R is None else float(R)
        a = cur.a if a is None else float(a)
        r0 = cur.r0 if r0 is None else float(r0)
        if type_ == "WoodsSaxon":
            if not R > 0 or not a > 0:
                raise ValueError(f"Woods-Saxon needs R > 0 and a > 0, "
                                 f"got R={R}, a={a}.")
        elif not r0 > 0:
            raise ValueError(f"Gaussian needs r0 > 0, got r0={r0}.")
        self.sess.set_potential(pair, type_, bool(enabled), V0, R, a, r0)
        if reinitialize:
            self._rebuild_after_potential_change()
        return self.nuclear_potential(pair)

    def clear_nuclear_potential(self, pair=0, reinitialize=True):
        """Drop a pair's own setting so it follows the default again.

        ``pair=0`` resets the default *and* every per-pair setting, which is
        how one gets back to a model with no hybrid potential at all.
        """
        pair = self._check_potential_pair(pair)
        self.sess.clear_potential(pair)
        if reinitialize:
            self._rebuild_after_potential_change()
        return self.nuclear_potential(pair)

    def _rebuild_after_potential_change(self):
        """Rebuild so the new potential reaches the model.

        A plain ``initialize()`` is not enough: it decides whether to reuse the
        external-capture integrals by looking for ``output/intEC.dat``, which is
        the right call at startup and the wrong one here -- the integrals on
        disk belong to the Coulomb functions the old potential produced.
        ``rebuild()`` recomputes them, the way a channel-radius change does.
        """
        if not self.sess.rebuild():
            raise RuntimeError(
                "AZURE2 could not rebuild with the new nuclear potential; "
                "the session is no longer usable.")
        self.configure()

    def _check_potential_pair(self, pair):
        pair = int(pair)
        if pair and not any(p.number == pair for p in self.pairs):
            raise KeyError(f"no pair {pair} in this model "
                           f"(have {[p.number for p in self.pairs]}); "
                           f"pair=0 is the default.")
        if pair < 0:
            raise ValueError(f"pair must be >= 0, got {pair}.")
        return pair
