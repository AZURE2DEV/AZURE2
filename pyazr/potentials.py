"""The hybrid nuclear potentials of a session, per particle pair (the methods of :class:`pyazr.azure2` that read and
set them, and write them as format-3 records; ``docs/source/user_guide/pyazr.rst``, "The hybrid nuclear potential, per
particle pair")."""

from .parameters import NuclearPotential


class _HybridPotentials:
    """Mixed into :class:`pyazr.azure2`: uses its ``sess``, ``pairs`` and ``configure``."""


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
        enabled, type_, V0, R, a, r0, own, rc, so, rm, step = self.sess.get_potential(pair)
        spin_orbit = None
        if len(so):
            spin_orbit = {"depth": float(so[0]), "radius": float(so[1]), "diffuseness": float(so[2]),
                          "on": "heavy" if so[3] else "light"}
        return NuclearPotential(pair=pair, enabled=bool(enabled), type=str(type_),
                                V0=float(V0), R=float(R), a=float(a),
                                r0=float(r0), own=bool(own), coulomb_radius=float(rc),
                                spin_orbit=spin_orbit, matching_radius=float(rm), step=float(step))

    def nuclear_potentials(self):
        """Every pair's resolved setting, plus the default under key ``0``."""
        out = {0: self.nuclear_potential(0)}
        for p in self.pairs:
            out[p.number] = self.nuclear_potential(p.number)
        return out

    def set_nuclear_potential(self, pair=0, type=None, enabled=True,
                              V0=None, R=None, a=None, r0=None,
                              coulomb_radius=None, spin_orbit=None, matching_radius=None,
                              step=None, reinitialize=True):
        """Give one particle pair its own hybrid nuclear potential.

        A nuclear potential belongs to a pair -- it bends the radial wave
        functions of that channel and no other -- so each pair carries its own,
        and ``pair=0`` sets the default that unnamed pairs inherit.  Anything
        left at ``None`` keeps the value the pair already resolves to, so
        turning one pair off is just
        ``set_nuclear_potential(2, enabled=False)``.

        ``type`` is ``"WoodsSaxon"`` (``V0``, ``R``, ``a``) or ``"Gaussian"``
        (``V0``, ``r0``); depths are MeV and lengths fm.  ``coulomb_radius``
        (fm) is the radius of the uniformly charged sphere of the Coulomb term;
        a pair without a setting of its own starts from 3.5 fm (the old
        ``<potential>`` block's value), a format-3 record states it.
        ``matching_radius`` (fm) is the radius where the solution inside the
        potential is matched to the Coulomb functions, ``step`` (fm) the radial
        step of the integration (the record's ``matching_radius`` and ``step``;
        the engine's defaults are 15 fm and 0.0275 fm); every channel table
        is keyed by them, so a change reaches the model like any other.

        ``spin_orbit`` adds a spin-orbit term to a pair with one spinning
        particle: a dict ``{"depth": V_so MeV, "radius": R_so fm,
        "diffuseness": a_so fm, "on": "light"|"heavy"}``, the term
        v_so(r) <2 l.i> with v_so = V_so (2.000 fm^2) (1/r) df/dr (FRESCO's
        convention; a strength quoted for l.s without the 2 is twice this V_so).
        ``None`` keeps the pair's term, ``False`` removes it.  A pair with two
        spinning particles (E_SO_TWO_SPINS), the term on a spin-0 particle
        (E_SO_ZERO_SPIN) and a term on the default (``pair=0``) are refused.
        The channel functions of such a pair depend on J:
        :meth:`coulomb_functions` then needs ``J``.

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
        rc = cur.coulomb_radius if coulomb_radius is None else float(coulomb_radius)
        rm = cur.matching_radius if matching_radius is None else float(matching_radius)
        step_fm = cur.step if step is None else float(step)
        if not rc > 0:
            raise ValueError(f"coulomb_radius must be > 0 fm, got {rc}.")
        if type_ == "WoodsSaxon":
            if not R > 0 or not a > 0:
                raise ValueError(f"Woods-Saxon needs R > 0 and a > 0, "
                                 f"got R={R}, a={a}.")
        elif not r0 > 0:
            raise ValueError(f"Gaussian needs r0 > 0, got r0={r0}.")
        if spin_orbit is None:
            so = cur.spin_orbit
        elif spin_orbit is False:
            so = None
        else:
            so = dict(spin_orbit)
            if set(so) != {"depth", "radius", "diffuseness", "on"}:
                raise ValueError(f"spin_orbit needs exactly the keys depth, radius, diffuseness, on; got {sorted(so)}.")
            if so["on"] not in ("light", "heavy"):
                raise ValueError(f"spin_orbit 'on' is 'light' or 'heavy', got {so['on']!r}.")
        so_list = [] if not so else [float(so["depth"]), float(so["radius"]), float(so["diffuseness"]),
                                     1.0 if so["on"] == "heavy" else 0.0]
        try:
            self.sess.set_potential(pair, type_, bool(enabled), V0, R, a, r0, rc, rm, step_fm, so_list)
        except Exception as e:  # the engine's refusals (SpinOrbitChecks.h) and parameter checks
            raise ValueError(str(e)) from None
        self._potentials_edited = True
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
        self._potentials_edited = True
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

    def potential_records(self):
        """The session's potentials as the ``potential`` records of format 3: ``(changed, lines, notes)``.

        ``lines`` is the body of ``<potentials>`` as the converter writes it (one record per particle pair with a
        potential, a default expanded into the pairs it applies to, every number in fm and MeV); ``changed`` says
        whether a pair's potential differs from the one the project's records give; ``notes`` lists the settings
        that write no record.  Raises ValueError for a project that is not in format 3.
        """
        try:
            changed, lines, notes = self.sess.potential_records()
        except Exception as e:  # the engine refuses a session whose project has no format-3 records
            raise ValueError(str(e)) from None
        return bool(changed), list(lines), list(notes)

    def _write_potential_records(self, model):
        """Into ``model`` (a :class:`pyazr.Project` of this session's format-3 project, or an AzrModel of an older
        one): the session's potentials, when set_nuclear_potential changed them; the section as read otherwise.
        Returns ``model``.  A project of format 1 or 2 whose potential this session changed is refused with the
        conversion command (ValueError): its ``<potential>`` block is not written; the format-3 records are."""
        from .project import Project
        if isinstance(model, Project):
            changed, lines, _ = self.potential_records()
            if changed:
                model._replace_record_section("potentials", lines)
        elif getattr(self, "_potentials_edited", False):
            from .described import _engine as _engine_module
            with open(self.file, "rb") as f:
                refusal = _engine_module().validate_format3_text(f.read(), self.file)[0]
            raise ValueError("save_fit: this session changed a hybrid potential, which is written only to a format-3 "
                             "project (its potential records); convert the project first.\n" + refusal)
        return model
