"""Physics of one session's level scheme (docs/upgrade/pyazr3.md section 2.7, "Per session"): the methods of
:class:`pyazr.Session` that need the session's pairs, radii and options, where the free functions of
:mod:`pyazr.physics` need none.  Mixed into the session class (``class Session(..., SessionPhysics)``)."""

from typing import NamedTuple

import numpy as np

from . import session_checks
from .locking import engine_call


class ScatteringLength(NamedTuple):
    """The s-wave scattering length of one channel spin of a particle pair
    (:meth:`Session.scattering_lengths`).  ``a0`` is complex when a lower channel
    is open at the pair's threshold; ``has_levels`` is False when no level has a
    width in the channel (``a0`` is then the hard-sphere value)."""
    pair: int
    s: float
    parity: int
    a0: complex
    has_levels: bool

    @property
    def jpi(self):
        """J^pi of the channel (J = s for l = 0) as the engine writes it, e.g. ``"1/2+"``."""
        from . import _azure2
        return _azure2.jpi_text(int(round(2 * self.s)), self.parity)


class SessionPhysics:
    """Mixed into :class:`pyazr.Session`: uses its ``_engine``, ``pairs``, ``parameters`` and ``transform``."""

    @engine_call
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
        1e-7 from the exact routine of :mod:`pyazr.physics`), the hybrid model's Numerov
        solution for a pair with a nuclear potential.  Where ``open`` is False every field but ``S`` is 0, and ``S``
        is the engine's closed-channel shift function at the channel radius (0 at another radius).

        Raises :class:`pyazr.ParameterIndexError` for an unknown pair, :class:`pyazr.DomainError` for a photon pair
        or an l that is not a whole number >= 0, ``ValueError`` for a J a spin-orbit pair needs or cannot make.
        """
        pair, l = session_checks.channel_request(self.pairs, pair, l)
        energies = np.asarray(e_cm, float)
        e = energies.ravel()
        req = np.concatenate([[pair, l, -1 if J is None else round(2 * J), radius, e.size], e])
        resp = np.asarray(self._engine.coulomb_functions(req), dtype=float)
        block = resp[1:].reshape(int(round(resp[0])), 7)
        out = {k: block[:, i].reshape(energies.shape)
               for i, k in enumerate(("F", "Fp", "G", "Gp", "P", "S", "delta_hs"))}
        out["open"] = (energies > 0) & ((out["F"] != 0) | (out["G"] != 0))
        out["e_cm"] = energies.copy()
        return out

    # -- dimensionless widths -------------------------------------------------

    @property
    @engine_call
    def mass_number(self):
        """Mass number A of the compound nucleus, from the particle pairs."""
        for p in self.pairs:
            if not p.is_photon:
                return int(round(p.M1 + p.M2))
        raise ValueError("no particle pair to take the compound mass from.")

    @engine_call
    def wigner_widths(self, params=None, eps=None):
        """``{free_index: Gamma_W}`` in eV: the Wigner-limit width ``Gamma_W = 2 P_l(E_r) gamma_W^2`` of every free
        particle channel that is open at its level's energy, at the fit vector ``params`` (default
        ``parameters.x0``).

        - ``gamma_W^2 = 3 hbar^2 / (2 mu a^2)`` (MeV) is the engine's Wigner limit of the channel, the row's
          ``wigner_limit`` (the Teichmann-Wigner sum-rule unit; :func:`pyazr.physics.wigner_limit`).
        - ``P_l(E_r)`` is the penetrability of this session's own channel functions (:meth:`channel_functions`: the
          Coulomb table when one is in use, the hybrid model's solution for a pair with a nuclear potential) at
          the pair's channel radius.
        - ``E_r`` is the pair's c.m. energy at the level's *physical* energy: the level's row of
          :meth:`transform`, the observed energy the program writes into ``<levels>``, minus the pair's separation
          and excitation energy.  In Brune's and Park's parametrizations that is the level's energy in ``params``; in
          the formal one ``params`` holds the formal pole energy, which is not where the observed width is defined.

        It is the number the AZURE2 GUI shows beside a channel, the denominator of ``theta^2 = Gamma_c / Gamma_W``
        (:meth:`dimensionless_widths`), and the same in the three parametrizations of one level scheme.

        No entry, as the key says: a fixed channel (no free index), a photon channel (its scale is the Weisskopf
        estimate), a channel closed at ``E_r`` (``E_r <= 0``: it has an ANC, not a width, and P = 0) and a channel
        without a Wigner limit in its record.

        ``eps`` is accepted and ignored.  Before pyazr 3.0.0.dev5 the number was a probe of the engine's
        transformation with every free amplitude set to ``eps``: equal to this one under Park and for a Brune level
        whose particle widths are all free, smaller by the level's ``1 / (1 + sum gamma^2 dS/dE)`` of its fixed
        widths otherwise, and not a property of the channel in the formal basis (``docs/upgrade/pyazr3.md``,
        section 11, has the table).
        """
        x = np.asarray(self.parameters.x0 if params is None else params, float)
        physical = self.transform(x)
        level_energy = {(p.jgroup, p.level): float(physical[p.index]) for p in self.parameters.of_kind("energy")}
        out = {}
        for p in self.parameters.of_kind("width"):
            if (p.fixed or p.free_index is None or p.free_index >= x.size
                    or p.radiation_type in ("E", "M") or p.wigner_limit is None):
                continue
            ex = level_energy.get((p.jgroup, p.level))
            if ex is None:
                continue
            pair = self.pairs.by_number(p.pair)
            e_r = ex - (pair.sep_energy + pair.excitation)
            if not e_r > 0:
                continue
            potential = self.nuclear_potential(pair.number)
            spin_orbit = potential.enabled and potential.spin_orbit is not None
            functions = self.channel_functions(pair.key, e_r, p.L, J=p.J if spin_orbit else None)
            if not functions["open"]:
                continue
            out[p.free_index] = 2.0 * float(functions["P"]) * p.wigner_limit * 1e6  # MeV -> eV
        return out

    def _level_energies(self, params=None):
        """``{(jgroup, level): Ex}`` -- level energies at ``params`` (MeV)."""
        x = np.asarray(self.parameters.x0 if params is None else params, float)
        energies = {}
        for p in self.parameters:
            if p.kind != "energy":
                continue
            key = (p.jgroup, p.level)
            if not p.fixed and p.free_index is not None and p.free_index < x.size:
                energies[key] = float(x[p.free_index])
            else:
                # Parameter.level_energy carries None for a level at Ex = 0
                # (the API's sentinel); as an energy that is a real 0.0 MeV.
                energies[key] = 0.0 if p.level_energy is None else p.level_energy
        return energies

    def _channel_open(self, param, x, pair, energies=None):
        """Is ``param``'s channel above threshold at the level's energy?"""
        energies = energies or self._level_energies(x)
        ex = energies.get((param.jgroup, param.level), param.level_energy)
        if ex is None:
            return False
        return ex > pair.sep_energy + pair.excitation

    @engine_call
    def dimensionless_widths(self, params=None, eps=None):
        """Every channel's width made dimensionless -- a :class:`WidthTable`.

        Particle channels get the Wigner limit and ``theta^2 = Gamma/Gamma_W``
        (plus ``theta^2_formal = gamma^2/gamma^2_W``), with ``Gamma_W = 2 P_l(E_r)
        gamma_W^2`` from :meth:`wigner_widths` (free open channels; a fixed channel
        has ``wigner_width`` and ``theta2`` ``None``); photon channels get the
        Weisskopf single-particle estimate and the strength in W.u.  ``params``
        is a fit vector -- pass a fit result to report
        the fitted widths, or leave it out for the ``.azr``'s own values.  ``eps``
        is accepted and ignored (:meth:`wigner_widths`).

        >>> t = m.dimensionless_widths(best)
        >>> print(t.photons.nonzero.table())
        >>> max(c.theta2 for c in t.particles if c.theta2)
        """
        from .widths import ChannelWidth, WidthTable, weisskopf_width

        x = np.asarray(self.parameters.x0 if params is None else params, float)
        phys = self.transform(x)[self.parameters.free_mask]
        gamma_w = self.wigner_widths(x, eps=eps)
        energies = self._level_energies(x)
        A = self.mass_number

        table = WidthTable()
        for p in self.parameters.of_kind("width"):
            pair = self.pairs.by_number(p.pair)
            ex = energies.get((p.jgroup, p.level), p.level_energy)
            free = p.free_index if (not p.fixed and p.free_index is not None) else None
            # the rwa of a fixed channel is not exposed by the API; its physical
            # value does not move with the fit, so report that and leave gamma out
            g = float(x[free]) if free is not None and free < x.size else None
            value = (float(phys[free]) if free is not None and free < phys.size
                     else p.value)
            photon = p.radiation_type in ("E", "M")
            row = ChannelWidth(
                name=p.name, index=p.index, free_index=p.free_index,
                fixed=p.fixed, jgroup=p.jgroup, level=p.level,
                jpi=p.jpi,  # the engine's J^pi text, "1/2+"
                level_energy=ex, pair=p.pair, L=p.L, S=p.S,
                radiation_type=p.radiation_type or "P",
                gamma=g, value=value, is_photon=photon,
                threshold=None if photon else pair.sep_energy + pair.excitation,
            )
            if photon:
                row.is_open = True
                row.e_gamma = None if ex is None else ex - pair.excitation
                row.weisskopf = weisskopf_width(p.radiation_type, p.L,
                                                row.e_gamma, A)
                if row.weisskopf and value is not None:
                    row.wu = abs(value) / row.weisskopf
            else:
                row.is_open = self._channel_open(p, x, pair, energies)
                row.wigner_gamma2 = p.wigner_limit
                if p.wigner_limit and g is not None:
                    row.theta2_formal = g ** 2 / p.wigner_limit
                row.wigner_width = gamma_w.get(free)
                if row.wigner_width and value is not None and row.is_open:
                    row.theta2 = abs(value) / row.wigner_width
            table.append(row)
        return table

    @engine_call
    def scattering_lengths(self, params=None):
        """The s-wave scattering lengths of every particle pair at a free RWA
        vector (default ``parameters.x0``): a list of :class:`ScatteringLength`,
        one per allowed channel spin (identical pairs: even s only).

        a0 = -a [M/(x^2 K1^2) - 2 I1/(x^2 K1)] with M = gamma^T A gamma at the
        pair's threshold (Paneru et al., PRC 99, 045807 (2019), Eq. 7), a (1 - M)
        for a neutral pair; the same in every parametrization.  Call it on each
        sample of an MCMC chain for the posterior of a0.  The measured ones of
        the ``scattering_length`` data of ``<modelData>`` are already part of
        the objective (``evaluate(x).terms["model data"]``)."""
        x = np.asarray(self.parameters.x0 if params is None else params, float)
        flat = np.asarray(self._engine.scattering_lengths_rwa(x), float).reshape(-1, 6)
        return [ScatteringLength(int(r[0]), float(r[1]), int(r[2]), complex(r[3], r[4]), bool(r[5])) for r in flat]

    @engine_call
    def park_norms(self, params=None):
        """Park formalism: the overlap J = 1 - sum_c gamma_c^2 dS_c/dE of every
        R-matrix level, in the order of the energy rows of ``parameters`` (``of_kind("energy")``).

        J is the squared norm of the level's basis state (Park 2021, Eq. 28)
        and the factor between the two parametrizations,
        ``gamma_Park = gamma_Brune * sqrt(J)``.  Every real Brune parameter set
        has 0 < J <= 1; a Park parameter set with J <= 0 has no R-matrix
        counterpart (its widths exceed what the channel radii allow) and the
        engine penalizes it (``evaluate(x).terms["Park J > 0"]``).  All ones in Brune mode.
        """
        x = np.asarray(self.parameters.x0 if params is None else params, float)
        return np.asarray(self._engine.park_norms(x), float)

    @engine_call
    def model_gradients(self, params):
        """Analytic ``d(observable)/d(parameter)`` for every calculated point.

        Returns a list with one ``(npoints, ncols)`` array per calculated
        segment, in the order and point ordering of ``evaluate(x).data``.  The
        columns are the free **R-matrix** parameters -- level energies and
        reduced-width amplitudes, in the order of the fit vector -- and the THM
        coherent background (``cbkg``) when there is one, which is exactly
        what ``output/covariance.dat`` spans; normalizations and energy shifts
        are omitted because no calculated observable depends on them.

        This is the sensitivity matrix an uncertainty band needs
        (``sigma^2 = g^T C g``).  Each row comes from one reverse-mode adjoint,
        so the whole thing costs about two forward evaluations no matter how
        many parameters the model has -- as against the ``2 * ncols`` forward
        passes a finite-difference estimate would take.  See
        :func:`pyazr.bands.uncertainty_bands`.

        THM (HOES) points have their own adjoint (see
        ``objective().jacobian``).

        Raises ``RuntimeError`` if the model contains a segment outside the
        supported analytic path.
        """
        resp = np.asarray(
            self._engine.calculate_model_gradients_rwa(
                np.asarray(params, float).ravel()), dtype=float)
        if resp.size == 1 and resp[0] == -1.0:
            raise RuntimeError("model_gradients: an analytically-unsupported "
                               "segment/config is present.")
        nseg = int(round(resp[0]))
        ncols = int(round(resp[1]))
        counts = [int(round(x)) for x in resp[2:2 + nseg]]
        flat = resp[2 + nseg:]
        if flat.size != sum(counts) * ncols:
            raise RuntimeError(
                f"model_gradients: got {flat.size} values, expected "
                f"{sum(counts) * ncols}.")
        out, start = [], 0
        for n in counts:
            out.append(flat[start:start + n * ncols].reshape(n, ncols))
            start += n * ncols
        return out

    @engine_call
    def ec_integrals(self, pair, energies):
        """External-capture radial integrals on an energy grid.

        ``pair`` is the entrance pair's key in ``<levels>``, ``energies`` its c.m. energies (MeV).  Every
        external-capture pathway the compound nucleus generates from this entrance pair is evaluated at every energy.
        Returns a list of dicts, one per pathway, each carrying its quantum numbers (``li``, ``lf``, ``si``, ``sf``,
        ``multipolarity``, ``radiation``), ``energy`` (the c.m. energies, MeV), the complex integral as ``value`` and
        its ``unit``, ``"fm^(L+1)"`` with L the multipolarity (``"fm^2"`` for an E1 or M1 pathway).

        The unit is the one the engine's use of the number gives it: the T-matrix element of the pathway is
        dimensionless and is the final state's ANC (fm^-1/2) times ``k_gamma^(L+1/2)`` (fm^-(L+1/2)) times this
        integral, with phases (``Site::ECPathwayAmplitude``, ``AMatrixFunc``).

        The pathways come in the engine's order (by final level, then spin group); the final state each belongs to
        is not part of the engine's answer, so two pathways to two final states with the same quantum numbers are
        told apart by their position in the list only.

        These integrals are the most expensive part of a capture calculation --
        which is why AZURE2 caches them.  Asking for them twice and watching
        :meth:`cache_stats` is the direct way to see that.
        """
        e = np.asarray(energies, float).ravel()
        req = np.concatenate([[pair, e.size], e])
        resp = np.asarray(self._engine.ec_integrals(req), dtype=float)
        npath = int(round(resp[0]))
        nE = int(round(resp[1]))
        stride = 6 + 2 * nE
        body = resp[2:]
        if body.size != npath * stride:
            raise RuntimeError(
                f"ec_integrals: got {body.size} values, expected {npath * stride}.")
        out = []
        for p in range(npath):
            b = body[p * stride:(p + 1) * stride]
            vals = b[6:].reshape(nE, 2)
            out.append({
                "li": int(round(b[0])),
                "lf": int(round(b[1])),
                "si": b[2] / 2.0,
                "sf": b[3] / 2.0,
                "multipolarity": int(round(b[4])),
                "radiation": "E" if b[5] > 0.5 else "M",
                "energy": e[:nE],
                "value": vals[:, 0] + 1j * vals[:, 1],
                "unit": f"fm^{int(round(b[4])) + 1}",
            })
        return out
