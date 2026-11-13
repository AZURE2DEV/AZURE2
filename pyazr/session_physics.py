"""Physics of one session's level scheme (docs/upgrade/pyazr3.md section 2.7, "Per session"): the methods of
:class:`pyazr.azure2` that need the session's pairs, radii and options, where the free functions of
:mod:`pyazr.physics` need none.  Mixed into the session class (``class azure2(..., SessionPhysics)``)."""

import numpy as np

from . import session_checks


class SessionPhysics:
    """Mixed into :class:`pyazr.azure2`: uses its ``sess`` and ``pairs``."""

    def channel_functions(self, pair, e_cm, l=0, radius=0.0, J=None):
        """Channel functions of one pair of this session at c.m. energies, as the session computes them.

        Parameters
        ----------
        pair : particle-pair key (1-based, as in ``<levels>``).
        e_cm : c.m. energies of the pair (MeV), a number or an array of any shape.
        l : orbital angular momentum (a whole number >= 0).
        radius : evaluation radius (fm); 0 means the pair's channel radius (where P, S and the hard-sphere phase
            belong).
        J : the channel's total J, required (and only allowed) for a pair whose potential has a spin-orbit term.

        Returns
        -------
        dict of arrays of the shape of ``e_cm``: ``F``, ``Fp``, ``G``, ``Gp`` (the Coulomb functions and their
        rho-derivatives), ``P`` (penetrability), ``S`` (shift function), ``delta_hs`` (hard-sphere phase shift, rad,
        = -phi of :func:`pyazr.physics.hard_sphere_phase`), ``open`` (bool: E > 0 and the engine's functions defined
        there) and ``e_cm``.  The numbers follow the session's options: the Coulomb table when one is in use (about
        1e-7 from the exact routine of :mod:`pyazr.physics`), GSL with ``use_gsl_coul``, the hybrid model's Numerov
        solution for a pair with a nuclear potential.  Where ``open`` is False every field but ``S`` is 0, and ``S``
        is the engine's closed-channel shift function at the channel radius (0 at another radius).

        Raises :class:`pyazr.ParameterIndexError` for an unknown pair, :class:`pyazr.DomainError` for a photon pair
        or an l that is not a whole number >= 0, ``ValueError`` for a J a spin-orbit pair needs or cannot make.
        The older name of this call is :meth:`coulomb_functions` (keys ``dF``, ``dG``, flat arrays).
        """
        pair, l = session_checks.channel_request(self.pairs, pair, l)
        energies = np.asarray(e_cm, float)
        e = energies.ravel()
        req = np.concatenate([[pair, l, -1 if J is None else round(2 * J), radius, e.size], e])
        resp = np.asarray(self.sess.coulomb_functions(req), dtype=float)
        block = resp[1:].reshape(int(round(resp[0])), 7)
        out = {k: block[:, i].reshape(energies.shape)
               for i, k in enumerate(("F", "Fp", "G", "Gp", "P", "S", "delta_hs"))}
        out["open"] = (energies > 0) & ((out["F"] != 0) | (out["G"] != 0))
        out["e_cm"] = energies.copy()
        return out
