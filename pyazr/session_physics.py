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
        radius : evaluation radius (fm); 0 means the channel radius (where P, S and the hard-sphere phase
            belong): the pair's, or the one a ``channel_radius`` record gives the pair's channels of this ``l``.
            When records with ``s`` or ``J`` give those channels more than one radius, 0 is refused
            (``ValueError``): name the radius.
        J : the channel's total J, required (and only allowed) for a pair whose functions depend on J: a spin-orbit
            term, or a ``wave`` of its potential with its own J.

        Returns
        -------
        dict of arrays of the shape of ``e_cm``: ``F``, ``Fp``, ``G``, ``Gp`` (the Coulomb functions and their
        rho-derivatives), ``P`` (penetrability), ``S`` (shift function), ``delta_hs`` (hard-sphere phase shift, rad,
        = -phi of :func:`pyazr.physics.hard_sphere_phase`), ``open`` (bool: E > 0 and the engine's functions defined
        there) and ``e_cm``.  The numbers follow the session's options: the Coulomb table when one is in use (about
        1e-7 from the exact routine of :mod:`pyazr.physics`), the hybrid model's Numerov
        solution for a pair with a nuclear potential.  In a partial wave with an interior background every field is
        the dressed one (``docs/upgrade/potential_background.md`` Eq. (9)): the functions the levels on top see, so
        that ``delta_hs`` = -phi~ and ``P`` = rho/|H~|^2 belong together; a wave with ``only bound_state`` gives the
        point-Coulomb functions.  Where ``open`` is False every field but ``S`` is 0, and ``S``
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

    @engine_call
    def potential_box_states(self, pair, l=0, J=None, radius=0.0, boundary=0.0, e_min=-60.0, e_max=20.0):
        """The box states of a pair's potential in one partial wave: the poles of its interior R-function.

        The eigenvalues E_p of the radial equation of the potential (the partial wave's own depth and width, as the
        session reads the ``potential`` record) inside ``radius`` with r u'(r)/u(r) = ``boundary`` there, their node
        numbers and residues gamma_p^2: R_pot(E) = sum_p gamma_p^2 / (E_p - E)
        (``docs/upgrade/potential_background.md`` sections 1 and 3; the node number is what ``subtract`` names).

        Parameters
        ----------
        pair : particle-pair key (1-based, as in ``<levels>``) of a pair with a potential.
        l : orbital angular momentum (whole number >= 0).
        J : the channel's total J; required when the potential depends on J (a spin-orbit term or a ``wave`` with its
            own J), ignored when it does not.
        radius : the box radius (fm); 0: the pair's channel radius.
        boundary : the boundary-condition constant B (dimensionless), as ``background interior boundary``.
        e_min, e_max : the energy window (MeV, c.m.).

        Returns
        -------
        dict: ``nodes`` (int array), ``E`` (MeV, c.m.), ``gamma2`` (MeV), one entry per state, ascending in E.

        Raises ``ValueError`` for a pair without a potential, a J needed or not possible, an empty window.
        """
        pair, l = session_checks.channel_request(self.pairs, pair, l)
        try:
            out = self._engine.potential_box_states(pair, l, -1 if J is None else round(2 * J), float(radius),
                                                    float(boundary), float(e_min), float(e_max))
        except Exception as e:
            raise ValueError(str(e)) from None
        return {"nodes": np.asarray(out["nodes"], int), "E": np.asarray(out["E"], float),
                "gamma2": np.asarray(out["gamma2"], float)}

    @engine_call
    def potential_bound_state(self, pair, l=0, J=None, energy=-1.0, window=2.0):
        """The bound state of a pair's potential in one partial wave nearest to an energy, and its ANC.

        The eigenvalue of the radial equation of the potential (0 to infinity: its own exterior) nearest ``energy``
        within ``window``, and the ANC of the normalised state, u = C W(r) beyond the potential: what the rule of
        ``capture interior`` compares the final level with (``docs/upgrade/potential_background.md`` section 5.3).

        Parameters
        ----------
        pair, l, J : as :meth:`potential_box_states`.
        energy : c.m. energy (MeV) below threshold (< 0) near which to look.
        window : half-width of the search (MeV).

        Returns
        -------
        dict: ``found`` (bool), ``E`` (MeV, c.m.; 0 when not found), ``anc`` (fm^-1/2; 0 when not found).

        Raises ``ValueError`` as :meth:`potential_box_states`, and for an energy >= 0.
        """
        pair, l = session_checks.channel_request(self.pairs, pair, l)
        try:
            out = self._engine.potential_bound_state(pair, l, -1 if J is None else round(2 * J), float(energy),
                                                     float(window))
        except Exception as e:
            raise ValueError(str(e)) from None
        return {"found": bool(out["found"]), "E": float(out["E"]), "anc": float(out["anc"])}

    @engine_call
    def jpi_block(self, jgroup, e_x, params=None):
        """One J^pi group of this session on excitation energies: channel functions, U, T, phase shifts, |A|, det A^-1.

        Parameters
        ----------
        jgroup : the J group, 1-based as ``parameters.out`` numbers them (``describe()["levels"]`` gives each level's).
        e_x : excitation energies E_x of the compound nucleus (MeV), a number or a 1-d array.
        params : the fit vector x (default ``parameters.x0``); the fixed values are the session's.

        Returns
        -------
        dict: ``jgroup``, ``J``, ``pi``, ``parametrization`` (standard, brune or park), ``E_x`` (MeV); ``levels`` (per
        level in the R-matrix: ``level``, ``E_level`` MeV, the energy of ``params`` in the session's parametrization);
        ``channels`` (per channel: ``channel``, ``pair``, ``radiation`` P, E or M, ``l``, ``s``, ``radius`` fm,
        ``threshold`` MeV of E_x, ``boundary`` B_c, and over E_x: ``E_c`` MeV the pair's c.m. energy, ``P``, ``S``,
        ``L_re`` = S - B_c and ``L_im`` = P as the level matrix reads them; a photon channel's are 0);
        ``elements`` (per pair c <= c' of particle channels: ``U_re``, ``U_im``, ``U_abs``, ``T_re``, ``T_im`` over
        E_x, NaN where either channel is closed); ``phase_shifts`` (per particle channel: ``delta`` in deg, the nuclear
        phase shift as the engine's ``phase_shift`` data kind forms it, and ``eta`` = |U_cc|); ``level_pairs`` (per
        pair of levels l <= m: ``A_abs`` = |A_lm| in MeV^-1); ``det_re``, ``det_im`` (det A^-1, MeV^n) and
        ``log10_abs_det``, least near each resonance; ``singular`` (bool over E_x).

        The arithmetic is the engine's (``physics/CollisionMatrix.h`` on the level route, the session's channel
        functions); a Coulomb table, when the session uses one, interpolates P and S to about 1e-7.  Raises
        ``ValueError`` for a J group that is not there, :class:`pyazr.ParameterVectorError` for a vector of another
        length.
        """
        x = session_checks.fit_vector(self.parameters.x0 if params is None else params, len(self.parameters.x0),
                                      "jpi_block")
        energies = np.atleast_1d(np.asarray(e_x, float)).ravel()
        return self._engine.jpi_block(self._all_rwa(x), int(jgroup), energies)

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
    def _observed_levels(self, x):
        """The engine's observed-level record at the fit vector ``x`` (``docs/upgrade/data_model.md`` section 10):
        ``{(jgroup, level, channel): row}``, a row ``[jgroup, level, channel, E (MeV, the level's physical energy),
        kind, W (MeV), value, P (the factor W was formed with: an open particle channel's penetrability), mu (u),
        a (fm)]``."""
        x = session_checks.fit_vector(x, len(self.parameters.x0), "observed levels")
        rows = np.asarray(self._engine.observed_levels_all_rwa(self._all_rwa(x)), float).reshape(-1, 10)
        return {(int(r[0]), int(r[1]), int(r[2])): r for r in rows}

    def _wigner_width(self, p, record):
        """Gamma_W = 2 P_l(E_r) gamma_W^2 (eV) of the particle channel of width parameter ``p`` from the record
        (:meth:`_observed_levels`), or None: a photon channel, no Wigner limit, or closed at the level's physical
        energy (E_r <= 0, or P = 0 just above a Coulomb threshold)."""
        if p.radiation_type in ("E", "M") or p.wigner_limit is None:
            return None
        row = record.get((p.jgroup, p.level, p.channel))
        if row is None:
            return None
        pair = self.pairs.by_number(p.pair)
        e_r = float(row[3]) - (pair.sep_energy + pair.excitation)
        P = float(row[7])
        if not (e_r > 0 and P > 0):
            return None
        return 2.0 * P * p.wigner_limit * 1e6  # MeV -> eV

    @engine_call
    def wigner_widths(self, params=None, eps=None):
        """``{free_index: Gamma_W}`` in eV: the Wigner-limit width ``Gamma_W = 2 P_l(E_r) gamma_W^2`` of every free
        particle channel that is open at its level's energy, at the fit vector ``params`` (default
        ``parameters.x0``).

        - ``gamma_W^2 = 3 hbar^2 / (2 mu a^2)`` (MeV) is the engine's Wigner limit of the channel, the row's
          ``wigner_limit`` (the Teichmann-Wigner sum-rule unit; :func:`pyazr.physics.wigner_limit`).
        - ``P_l(E_r)`` is the penetrability the engine formed the channel's observed width with (its observed-level
          record): this session's own channel functions (the Coulomb table when one is in use, the hybrid model's
          solution for a pair with a nuclear potential) at the pair's channel radius, the number
          :meth:`channel_functions` gives at ``E_r``.
        - ``E_r`` is the pair's c.m. energy at the level's *physical* energy: the level's row of
          :meth:`transform`, the observed energy the program writes into ``<levels>``, minus the pair's separation
          and excitation energy.  In Brune's and Park's parametrizations that is the level's energy in ``params``; in
          the formal one ``params`` holds the formal pole energy, which is not where the observed width is defined.

        It is the number the AZURE2 GUI shows beside a channel, the denominator of ``theta^2 = Gamma_c / Gamma_W``
        (:meth:`dimensionless_widths`), and the same in the three parametrizations of one level scheme.

        No entry, as the key says: a fixed channel (no free index; :meth:`dimensionless_widths` gives its
        ``wigner_width``), a photon channel (its scale is the Weisskopf estimate), a channel closed at ``E_r``
        (``E_r <= 0``: it has an ANC, not a width, and P = 0) and a channel without a Wigner limit in its record.

        ``eps`` is accepted and ignored.  Before pyazr 3.0.0.dev5 the number was a probe of the engine's
        transformation with every free amplitude set to ``eps``: equal to this one under Park and for a Brune level
        whose particle widths are all free, smaller by the level's ``1 / (1 + sum gamma^2 dS/dE)`` of its fixed
        widths otherwise, and not a property of the channel in the formal basis (``docs/upgrade/pyazr3.md``,
        section 11, has the table).
        """
        x = np.asarray(self.parameters.x0 if params is None else params, float)
        record = self._observed_levels(x)
        out = {}
        for p in self.parameters.of_kind("width"):
            if p.fixed or p.free_index is None or p.free_index >= x.size:
                continue
            width = self._wigner_width(p, record)
            if width is not None:
                out[p.free_index] = width
        return out

    @engine_call
    def dimensionless_widths(self, params=None, eps=None):
        """Every channel's width made dimensionless -- a :class:`WidthTable`.

        Particle channels get the Wigner limit and ``theta^2 = Gamma/Gamma_W``
        (plus ``theta^2_formal = gamma^2/gamma^2_W`` for a free channel), with
        ``Gamma_W = 2 P_l(E_r) gamma_W^2`` as :meth:`wigner_widths` defines it,
        for free and fixed channels open at the level's energy; photon channels
        get the Weisskopf single-particle estimate and the strength in W.u.
        ``level_energy`` is the level's physical energy (the row of
        :meth:`transform`: in the formal basis not the pole energy of the fit
        vector), and ``is_open`` and a photon's ``e_gamma`` are taken there.
        ``params`` is a fit vector -- pass a fit result to report the fitted
        widths, or leave it out for the ``.azr``'s own values.  ``eps`` is
        accepted and ignored (:meth:`wigner_widths`).

        Every row's ``value`` -- free and fixed -- is the physical value at
        ``params``, the row of :meth:`transform` (``params``): a fixed
        channel holds its amplitude, so under Brune's and the formal
        parametrization its observed width moves with the level's other
        widths (and in every basis with the level's energy, through P).

        Before pyazr 3.0.0.dev6 a fixed channel had ``wigner_width`` and
        ``theta2`` ``None``, its ``value`` (and so ``theta2`` and a photon's
        ``wu``) was the value as the session was made (``Parameter.value``),
        not at ``params``, and in the formal basis ``level_energy``,
        ``is_open`` and ``e_gamma`` were taken at the formal pole energy.

        >>> t = m.dimensionless_widths(best)
        >>> print(t.photons.nonzero.table())
        >>> max(c.theta2 for c in t.particles if c.theta2)
        """
        from .widths import ChannelWidth, WidthTable, weisskopf_width

        x = np.asarray(self.parameters.x0 if params is None else params, float)
        phys = self.transform(x)  # one entry per row of parameters, fixed ones included: the observed values at x
        record = self._observed_levels(x)
        A = self.mass_number

        table = WidthTable()
        for p in self.parameters.of_kind("width"):
            pair = self.pairs.by_number(p.pair)
            ex = float(record[(p.jgroup, p.level, p.channel)][3])  # the level's physical energy (MeV)
            free = p.free_index if (not p.fixed and p.free_index is not None) else None
            g = float(x[free]) if free is not None and free < x.size else None
            # the physical value at x, free and fixed alike: a fixed channel holds its amplitude, so its observed
            # width moves with the level's other parameters (pyazr3.md section 13.5)
            value = float(phys[p.index])
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
                row.e_gamma = ex - pair.excitation
                row.weisskopf = weisskopf_width(p.radiation_type, p.L,
                                                row.e_gamma, A)
                if row.weisskopf and value is not None:
                    row.wu = abs(value) / row.weisskopf
            else:
                row.is_open = ex > pair.sep_energy + pair.excitation
                row.wigner_gamma2 = p.wigner_limit
                if p.wigner_limit and g is not None:
                    row.theta2_formal = g ** 2 / p.wigner_limit
                row.wigner_width = self._wigner_width(p, record)
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
