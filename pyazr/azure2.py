"""High-level driver for one in-process AZURE2 engine instance.

An :class:`azure2` object owns one ``_azure2.Session`` -- the AZUREAPI C++
object bound into Python via pybind11 -- loaded from a single ``.azr`` file.
There are no subprocesses, no sockets and no instance pools to manage: the
engine runs inside the interpreter, so a session is exactly one Python object.
Opening several ``azure2()`` objects in one interpreter is fine, each is
independent; drive them from one thread, as the engine is not reentrant.
"""

import os
from contextlib import contextmanager

import numpy as np

try:
    from . import _azure2
except ImportError as err:                            # pragma: no cover
    raise ImportError(
        "the AZURE2 engine module (pyazr._azure2) is not built.  It is a C++ "
        "extension, not pure Python: configure with CMake (USE_API=ON, the "
        "default) and build -- the module lands in pyazr/."
    ) from err

from .parameters import (NuclearPotential, Pair, PairSet, Parameter,
                         ParameterSet)


@contextmanager
def _in_dir(path):
    """Run the block with the process cwd at ``path``, then put it back."""
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


class _Engine:
    """A ``_azure2.Session`` whose every call is made from the model's directory.

    AZURE2 names ``output/``, ``checks/`` and its data files relative to the
    process' working directory, and a process has only one.  Entering that
    directory per call rather than holding it for the session's lifetime is
    what lets two sessions in different directories coexist -- and it leaves
    the caller's own cwd alone, as the subprocess-based API did.
    """

    def __init__(self, session, cwd):
        self._session = session
        self._cwd = cwd

    def __getattr__(self, name):
        # __getattr__ runs before _session exists (unpickling, a failed
        # __init__), and forwarding then recurses forever.
        if name.startswith("_"):
            raise AttributeError(name)
        attr = getattr(self._session, name)
        if not callable(attr):
            return attr

        def call(*args, **kwargs):
            with _in_dir(self._cwd):
                return attr(*args, **kwargs)
        return call


class azure2:

    def __init__(self, file, cwd=None, data_mode=True, use_brune=True,
                 ignore_externals=True, transform=True,
                 use_long_wavelength=True, use_gsl_coul=False, use_rmc=False):
        """Build the R-matrix engine for ``file`` in this interpreter.

        Parameters
        ----------
        file : the ``.azr`` configuration file.
        cwd : working directory for the run.  A model names its ``output/``,
            ``checks/`` and data paths *relative to the process cwd*, so this
            defaults to the ``.azr`` file's own directory -- the convention the
            file itself is written against.  Each engine call enters it and
            leaves again, so the caller's cwd is never left changed.
        data_mode : evaluate the ``<segmentsData>`` blocks (the startup state);
            ``False`` starts in extrapolation mode (``<segmentsTest>``).
        use_brune, ignore_externals, transform, use_long_wavelength :
            R-matrix formalism options, matching the CLI's ``--use-brune``,
            ``--ignore-externals``, ``--no-transform`` and ``--no-long-wavelength``
            (the ``use_*`` options are ON by default, ``transform`` too).
        use_gsl_coul : use GSL's Coulomb functions instead of AZURE2's own.
        use_rmc : use the R-matrix-with-channels (RMC) formalism instead of
            Brune; mutually exclusive with ``use_brune`` (Brune wins).

        Raises
        ------
        _azure2.AZURE2Error
            if the file is missing or the model will not initialize.
        """
        self._sess = None
        self.file = os.path.abspath(file)
        self.cwd = os.path.abspath(cwd if cwd is not None else (
            os.path.dirname(self.file) or "."))
        if not os.path.isdir(self.cwd):
            raise FileNotFoundError(f"working directory {self.cwd!r} does not exist.")

        # Kept so a session can be reproduced elsewhere -- a worker process
        # rebuilding this same model needs the formalism options, not just the
        # file (see pyazr.bands.sensitivities with nprocs > 1).
        self.options = dict(
            data_mode=bool(data_mode), use_brune=bool(use_brune),
            ignore_externals=bool(ignore_externals), transform=bool(transform),
            use_long_wavelength=bool(use_long_wavelength),
            use_gsl_coul=bool(use_gsl_coul), use_rmc=bool(use_rmc))

        opts = _azure2.RuntimeOptions()
        for name, value in self.options.items():
            setattr(opts, name, value)
        # Construction reads the .azr and its data files, so it runs from the
        # model's directory like every call afterwards.
        with _in_dir(self.cwd):
            self._sess = _Engine(_azure2.Session(self.file, opts), self.cwd)

        self.mode = "data" if data_mode else "extrap"
        self.configure()

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # -- lifecycle ------------------------------------------------------------

    @property
    def sess(self):
        """The underlying engine, or a clear error once it has been released."""
        if self._sess is None:
            raise RuntimeError(
                "this azure2 session is closed; open a new one for more work.")
        return self._sess

    def close(self):
        """Release the engine.

        The C++ model is several MB per session and is freed here rather than
        whenever the garbage collector gets to it -- which matters in a loop
        over model variants, where the collector may hold two at once.
        """
        self._sess = None

    def is_alive(self):
        """Is the engine still open? False once close() has run."""
        return self._sess is not None

    def configure(self):
        """Re-read the parameter, pair and dataset metadata from the engine, after anything that rebuilt the model."""
        s = self.sess
        self.nsegments = int(s.update_data())
        rng = range(self.nsegments)
        self.energies = [s.data_energies(i) for i in rng]
        self.excitation_energies = [s.data_excitation_energies(i) for i in rng]
        self.angles = [s.data_angles(i) for i in rng]
        self.cross = [s.data_segments(i) for i in rng]
        self.cross_err = [s.data_segments_errors(i) for i in rng]
        self.conv = [s.data_conv(i) for i in rng]
        self.sfactor = [self.cross[i] * self.conv[i] for i in rng]
        self.sfactor_err = [self.cross_err[i] * self.conv[i] for i in rng]
        self.params = s.params_values()
        self.params_rwa = s.params_values_rwa()
        self.fixed_params = s.params_fixed()
        self._parameters = None
        self._pairs = None

    # -- parameter metadata ---------------------------------------------------

    @property
    def parameters(self):
        """A :class:`ParameterSet` describing every fit parameter.

        Each entry says what the parameter is: for R-matrix parameters which
        level (J^pi, energy) it belongs to, and for widths which channel (L, S,
        particle pair, radiation type).  Built lazily and cached; call
        :meth:`refresh_parameters` to rebuild after the model changes.

        Examples
        --------
        >>> azr.parameters.widths                 # all width parameters
        >>> azr.parameters.free                    # only the non-fixed ones
        >>> azr.parameters.by_level(jgroup=1)      # one level's energy + widths
        >>> print(azr.parameters.table())          # readable overview
        """
        if self._parameters is None:
            self._parameters = self._build_parameters()
        return self._parameters

    def refresh_parameters(self):
        """Re-fetch and rebuild the cached :attr:`parameters`."""
        self._parameters = self._build_parameters()
        return self._parameters

    @property
    def pairs(self):
        """A :class:`PairSet` describing every particle pair (channel).

        Each :class:`Pair` carries the two constituents' spins and parities and
        a flag for the reaction entrance pair.  A width parameter's ``pair``
        attribute is the :attr:`Pair.number`.  Built lazily and cached.
        """
        if self._pairs is None:
            self._pairs = self._build_pairs()
        return self._pairs

    def refresh_pairs(self):
        """Re-fetch and rebuild the cached :attr:`pairs`."""
        self._pairs = self._build_pairs()
        return self._pairs

    # -- channel radius -------------------------------------------------------

    def set_channel_radius(self, pair, radius):
        """Change one particle pair's channel radius.

        ``pair`` is the 1-based :attr:`Pair.number`, ``radius`` is in fm.  This
        is not a light edit: the radius sets the matching surface, so AZURE2
        rebuilds the compound nucleus and the data from the ``.azr``, redoes the
        penetrabilities, shift functions, boundary conditions and Wigner limits,
        and **recomputes every external-capture integral** rather than reading
        back ``output/intEC*`` (which belong to the old radius).  Expect it to
        cost about as much as starting a fresh instance.

        Everything cached on the Python side is re-read afterwards, so
        :attr:`params_rwa`, :attr:`parameters` (including ``wigner_limit``),
        :attr:`pairs` and the data arrays are consistent with the new radius.

        The reduced widths keep their numerical values, but a reduced width
        means something different at a different radius -- **the model is no
        longer fitted**.  Refit before reading anything off it.

        Note this changes only the running session; the ``.azr`` on disk is
        untouched.  To persist a radius (and to scan several radii, which is
        the safer pattern) use :meth:`pyazr.AzrModel.set_channel_radius` and
        write a new file.
        """
        pair = int(pair)
        if not any(p.number == pair for p in self.pairs):
            raise KeyError(f"no pair {pair} in this model "
                           f"(have {[p.number for p in self.pairs]}).")
        if not radius > 0:
            raise ValueError(f"channel radius must be positive, got {radius}.")
        if not self.sess.set_radius(pair, float(radius)):
            raise RuntimeError(
                f"AZURE2 could not rebuild with pair {pair} at "
                f"{radius} fm; the session is no longer usable.")
        self.configure()
        return self.pairs.by_number(pair).channel_radius

    # -- hybrid nuclear potential, per particle pair --------------------------

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

    def _build_pairs(self):
        flat = np.asarray(self.sess.pairs_info(), dtype=float)
        nfields = Pair._NFIELDS
        records = flat.reshape(-1, nfields)
        return PairSet(Pair.from_record(rec) for rec in records)

    @property
    def level_scheme(self):
        """A structured, printable :class:`~pyazr.scheme.LevelScheme`.

        Groups the model the way it reads physically -- particle pairs, then
        J-groups, then levels and their channels (L, S, radiation type, partial
        width, fixed flag, Wigner limit).  ``print(azr.level_scheme)`` gives a
        human-readable overview.  Read-only; to add/remove levels and write the
        result to a file see :class:`pyazr.AzrModel`.
        """
        from .scheme import LevelScheme
        return LevelScheme.from_azr(self)

    @property
    def datasets(self):
        """Per-segment dataset provenance parsed from the ``.azr`` file.

        A :class:`~pyazr.datasets.SegmentSet`: for each data segment, the data
        file it came from, the reaction channel (entrance/exit pairs), energy /
        angle range, observable type, and normalization systematic error.
        ``print(azr.datasets.table())`` gives an overview; ``azr.datasets
        .sys_errors()`` returns the per-segment systematics the fits use.
        """
        from .datasets import SegmentSet
        return SegmentSet.from_file(self.file)

    @property
    def extrapolations(self):
        """Per-segment extrapolation grids parsed from the ``.azr`` file.

        A :class:`~pyazr.datasets.TestSegmentSet`: for each ``<segmentsTest>``
        entry, the reaction channel, the energy / angle grid, and the observable
        type.  These describe the segments AZURE2 reports in extrapolation mode
        (:meth:`extrap_mode`), in the same order, so segment ``i`` of this set
        corresponds to index ``i`` of ``calculate``/``calculate_energies``.
        """
        from .datasets import TestSegmentSet
        return TestSegmentSet.from_file(self.file)

    def _build_parameters(self):
        s = self.sess
        n = len(self.fixed_params)

        flat = np.asarray(s.parameter_info(), dtype=float)
        nfields = Parameter._NFIELDS
        records = flat.reshape(-1, nfields)
        if records.shape[0] != n:
            raise RuntimeError(
                f"parameter_info returned {records.shape[0]} parameters but "
                f"params_fixed reported {n}."
            )

        params = ParameterSet()
        free_counter = 0
        for i in range(n):
            name = s.params_names(i)
            fixed = bool(self.fixed_params[i])
            free_index = None if fixed else free_counter
            if not fixed:
                free_counter += 1
            params.append(
                Parameter.from_record(i, name, records[i], free_index)
            )
        return params

    # -- physical levels (turn a resonance on/off) ----------------------------

    def physical_levels(self):
        """The model's physical levels (resonances / poles) as
        :class:`~pyazr.parameters.LevelKey` objects, ordered ``(jgroup, level)``.

        Each key is a handle you can pass to :meth:`without_level` /
        :meth:`only_level` to switch that resonance off or isolate it.
        """
        return list(self.parameters.by_physical_level().keys())

    def _match_level(self, level=None, jpi=None, energy=None, tol=1e-2):
        """Resolve a level selector to a single :class:`LevelKey`.

        Accepts a ``LevelKey`` directly (``level=``), or a ``jpi`` string
        (e.g. ``"5/2-"``) and/or an ``energy`` (MeV, matched within ``tol``).
        Raises if the selector is ambiguous or matches nothing.
        """
        keys = self.physical_levels()
        if level is not None:
            if level in keys:
                return level
            raise KeyError(f"{level!r} is not a level of this model.")
        hits = [k for k in keys
                if (jpi is None or k.jpi == jpi)
                and (energy is None or (k.energy is not None
                                        and abs(k.energy - energy) <= tol))]
        if not hits:
            raise KeyError(f"no level matches jpi={jpi} energy={energy}.")
        if len(hits) > 1:
            raise KeyError(f"ambiguous selector jpi={jpi} energy={energy}: "
                           f"matches {[str(h) for h in hits]}.")
        return hits[0]

    def level_free_width_indices(self, level=None, jpi=None, energy=None,
                                 tol=1e-2):
        """``free_index`` list of a level's non-fixed reduced-width parameters.

        These are the positions, within the free-parameter vector
        (:attr:`params_rwa` order), of the reduced-width amplitudes that carry
        the level's coupling to its channels.  Zeroing them removes the level
        from the calculation (its fixed widths, if any, are already inert).
        """
        key = self._match_level(level=level, jpi=jpi, energy=energy, tol=tol)
        ps = self.parameters.by_physical_level()[key]
        return [p.free_index for p in ps
                if p.kind == "width" and not p.fixed and p.free_index is not None]

    def without_level(self, params, level=None, jpi=None, energy=None,
                      tol=1e-2):
        """A copy of the free-parameter vector with one level switched OFF.

        The level's reduced widths are set to zero, decoupling it from every
        channel; all other parameters are untouched.  ``params`` is a free
        (``params_rwa``-order) vector.  Use with :meth:`calculate_rwa` /
        :meth:`residual_jacobian` to see the model *without* that resonance.
        """
        import numpy as _np
        out = _np.array(params, dtype=float)
        out[self.level_free_width_indices(level, jpi, energy, tol)] = 0.0
        return out

    def only_level(self, params, level=None, jpi=None, energy=None, tol=1e-2):
        """A copy of the free-parameter vector with ONLY one level active.

        Every *other* level's reduced widths are zeroed, leaving the target
        level's widths (and all energies / normalizations) in place -- the bare
        contribution of one resonance on top of the non-resonant (Coulomb /
        hard-sphere) background.  Compare against :meth:`without_level` and the
        full model to read off interference.
        """
        import numpy as _np
        keep = set(self.level_free_width_indices(level, jpi, energy, tol))
        out = _np.array(params, dtype=float)
        for p in self.parameters.widths:
            if not p.fixed and p.free_index is not None and p.free_index not in keep:
                out[p.free_index] = 0.0
        return out

    # -- index queries --------------------------------------------------------

    def norm_indices(self):
        """Free-vector positions of the normalization parameters."""
        return self.sess.normalization_indices()

    def shift_indices(self):
        """Free-vector positions of the energy-shift parameters."""
        return self.sess.energy_shift_indices()

    def energy_indices(self):
        """Packed indices of the free R-matrix level-energy parameters.

        Returns the positions of the level-energy parameters within
        ``params_rwa`` (the non-fixed parameter vector), the same convention
        as :meth:`norm_indices` and :meth:`shift_indices`.  Unlike those --
        which the C++ side derives by substring-matching parameter names --
        this uses the structured ``type`` code from ``parameter_info``
        (``kind == "energy"``), so it never depends on how energies are named.
        """
        return [p.free_index for p in self.parameters
                if p.kind == "energy" and not p.fixed]

    # -- param.sav helpers ----------------------------------------------------

    # param.sav lives beside the model, so these read and write from the
    # session's directory rather than wherever the caller happens to be.

    def update_rwa_params_from_sav(self):
        """Reload params_rwa from the run's param.sav."""
        with _in_dir(self.cwd):
            all_rwa_params = np.loadtxt('output/param.sav', usecols=(1,))
        self.params_rwa = []
        for i in range(len(all_rwa_params)):
            if self.fixed_params[i]:
                continue
            else:
                self.params_rwa.append(all_rwa_params[i])

    def update_sav_from_rwa_params(self, best):
        """Write a free RWA vector back out to param.sav."""
        params_full = []
        with _in_dir(self.cwd):
            with open('output/param.sav', 'r') as f:
                for line in f.readlines():
                    l = line.split()
                    params_full.append([l[0], float(l[1]), float(l[2])])

            idx = 0
            for i in range(len(params_full)):
                if self.fixed_params[i]:
                    continue
                else:
                    params_full[i][1] = best[idx]
                    idx += 1

            with open('output/param.sav.new', 'w') as f:
                for param in params_full:
                    f.write(f"{param[0]} {param[1]} {param[2]}\n")

    # -- transforms -----------------------------------------------------------

    def transform_rwa(self, params):
        """Map a reduced-width-amplitude vector to physical parameters: level energies in MeV, partial widths in eV, ANCs in fm^-1/2."""
        return self.sess.transform_rwa(params)

    def transform_all_rwa(self, params):
        """As transform_rwa, over every parameter rather than the free ones."""
        return self.sess.transform_all_rwa(params)

    def transform_physical(self, params):
        """The inverse of transform_rwa: physical parameters back to reduced-width amplitudes."""
        raise NotImplementedError(
            "AZURE2 exposes no physical->RWA transform; only RWA->physical "
            "is available (transform_rwa / transform_all_rwa).")

    # -- dimensionless widths -------------------------------------------------

    @property
    def mass_number(self):
        """Mass number A of the compound nucleus, from the particle pairs."""
        for p in self.pairs:
            if not p.is_photon:
                return int(round(p.M1 + p.M2))
        raise ValueError("no particle pair to take the compound mass from.")

    def wigner_widths(self, params=None, eps=1e-4):
        """``{free_index: Gamma_W}`` in eV -- the GUI's Wigner width per channel.

        ``Gamma_W = 2 P_l(E_r) gamma^2_W`` is what the AZURE2 GUI shows next to
        each particle channel, and ``theta^2 = Gamma_c / Gamma_W``.  It needs
        the penetrability at the level energy, which the API does not expose, so
        this asks AZURE2 for it: every reduced-width amplitude is set to a tiny
        ``eps`` and the model transformed, which makes the level-shift
        denominator 1 and the J-group level mixing vanish to O(eps^2), leaving
        ``Gamma_c(eps) = 2 P_c eps^2``.  One extra transform, evaluated with
        AZURE2's own radii, boundary conditions and Coulomb functions.

        Only *free* channels appear: the transform returns non-fixed parameters
        only.  Closed (sub-threshold) channels are omitted -- there
        ``P_l = 0`` and AZURE2 reports an ANC rather than a width.  Level
        energies are taken from ``params`` (default: the current
        :attr:`params_rwa`), so the limits sit at the fitted resonance energies.
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        probe = x.copy()
        for p in self.parameters.widths:
            if p.fixed or p.free_index is None or p.free_index >= probe.size:
                continue
            # photon channels are linear in the amplitude (external capture),
            # so they must not pollute the probe
            probe[p.free_index] = 0.0 if p.radiation_type in ("E", "M") else eps
        probed = np.asarray(self.transform_rwa(probe), float)

        out = {}
        for p in self.parameters.widths:
            if (p.fixed or p.free_index is None or p.free_index >= probed.size
                    or p.radiation_type in ("E", "M") or p.wigner_limit is None):
                continue
            pair = self.pairs.by_number(p.pair)
            if not self._channel_open(p, x, pair):
                continue
            two_p = abs(float(probed[p.free_index])) / eps ** 2
            out[p.free_index] = two_p * p.wigner_limit
        return out

    def _level_energies(self, params=None):
        """``{(jgroup, level): Ex}`` -- level energies at ``params`` (MeV)."""
        x = np.asarray(self.params_rwa if params is None else params, float)
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

    def dimensionless_widths(self, params=None, eps=1e-4):
        """Every channel's width made dimensionless -- a :class:`WidthTable`.

        Particle channels get the Wigner limit and ``theta^2 = Gamma/Gamma_W``
        (plus ``theta^2_formal = gamma^2/gamma^2_W``); photon channels get the
        Weisskopf single-particle estimate and the strength in W.u.  ``params``
        is a free (``params_rwa``-order) vector -- pass a fit result to report
        the fitted widths, or leave it out for the ``.azr``'s own values.

        >>> t = m.dimensionless_widths(best)
        >>> print(t.photons.nonzero.table())
        >>> max(c.theta2 for c in t.particles if c.theta2)
        """
        from .widths import ChannelWidth, WidthTable, weisskopf_width

        x = np.asarray(self.params_rwa if params is None else params, float)
        phys = np.asarray(self.transform_rwa(x), float)
        gamma_w = self.wigner_widths(x, eps=eps)
        energies = self._level_energies(x)
        A = self.mass_number

        table = WidthTable()
        for p in self.parameters.widths:
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
                fixed=p.fixed, jgroup=p.jgroup, level=p.level, jpi=p.jpi,
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

    # -- chi-squared ----------------------------------------------------------

    def calculate_chi2_rwa(self, params):
        """Total data chi-squared at a free RWA vector, as a one-element list. Excludes the normalization and energy-shift penalties AZURE2's own objective adds."""
        return [float(self.sess.calculate_chi2_rwa(params))]

    def calculate_chi2(self, params):
        """As calculate_chi2_rwa, taking the physical parameter vector instead."""
        return [float(self.sess.calculate_chi2_physical(params))]

    def chi2_and_grad(self, params):
        """Value and analytic gradient of the (data) chi-squared.

        Parameters
        ----------
        params : the non-fixed RWA parameters (energies, reduced widths,
            normalizations), in the same order as ``self.params_rwa``.

        Returns
        -------
        (float, numpy.ndarray)
            ``(chi2, grad)`` with one gradient entry per input parameter.
            Energies / reduced widths / normalizations are analytic; energy
            shifts are finite-differenced.  (For a Gaussian log-likelihood use
            ``lnL = -0.5*(chi2 + const)`` and ``grad_lnL = -0.5*grad``.)
        """
        resp = np.asarray(
            self.sess.calculate_chi2_grad_rwa(np.asarray(params, float).ravel()),
            dtype=float)
        return float(resp[0]), np.asarray(resp[1:], dtype=float)

    def residual_jacobian(self, params):
        """Standardized residuals and their Jacobian.

        ``r_i = (fit_i - data_i*n)/(cmErr_i*n)`` so ``sum(r_i**2) == chi2``.

        Returns ``(r, J)`` with ``r`` shape ``(n_res,)`` and ``J`` shape
        ``(n_res, n_params)``; columns match the non-fixed RWA parameters (the
        input ordering).

        Level-energy, reduced-width and normalization columns come from the
        reverse-mode adjoint, so that block costs ~2 forward evaluations
        regardless of the parameter count -- for Gauss-Newton /
        Levenberg-Marquardt.

        Energy-shift columns are finite-differenced, at two extra residual
        evaluations each.  A shift translates the energy axis of a whole
        segment, so the derivative wanted is d(model)/dE, and AZURE2 applies a
        shift by rebuilding every energy-dependent quantity of the affected
        points; there is no cheap analytic route through the forward code.  A
        model with many free shifts is dominated by those columns.

        (Before pyazr 2.7 these columns came back as zero, which a
        least-squares driver accepts silently by never moving those
        parameters.)
        """
        resp = np.asarray(
            self.sess.calculate_residual_jacobian_rwa(
                np.asarray(params, float).ravel()), dtype=float)
        if resp.size == 1 and resp[0] == -1.0:
            raise RuntimeError("residual_jacobian: an analytically-unsupported "
                               "segment/config is present.")
        n_res = int(round(resp[0]))
        n_cols = int(round(resp[1]))
        r = resp[2:2 + n_res]
        J = resp[2 + n_res:].reshape(n_res, n_cols)
        return r, J

    def model_gradients(self, params):
        """Analytic ``d(observable)/d(parameter)`` for every calculated point.

        Returns a list with one ``(npoints, ncols)`` array per calculated
        segment, in the order and point ordering of :meth:`calculate_rwa`.  The
        columns are the free **R-matrix** parameters -- level energies and
        reduced-width amplitudes, in ``params_rwa`` order -- which is exactly
        what ``output/covariance.dat`` spans; normalizations and energy shifts
        are omitted because no calculated observable depends on them.

        This is the sensitivity matrix an uncertainty band needs
        (``sigma^2 = g^T C g``).  Each row comes from one reverse-mode adjoint,
        so the whole thing costs about two forward evaluations no matter how
        many parameters the model has -- as against the ``2 * ncols`` forward
        passes a finite-difference estimate would take.  See
        :func:`pyazr.bands.uncertainty_bands`.

        Raises ``RuntimeError`` if the model contains a segment outside the
        supported analytic path.
        """
        resp = np.asarray(
            self.sess.calculate_model_gradients_rwa(
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

    # -- the external region, and the caches that make it affordable ----------

    def coulomb_functions(self, pair, energies, L=0, radius=0.0):
        """Coulomb wave functions on an energy grid.

        Parameters
        ----------
        pair : particle-pair key (1-based, as in the .azr).
        energies : centre-of-mass energies in MeV.
        L : orbital angular momentum.
        radius : evaluation radius in fm; 0 means the pair's channel radius,
            which is where penetrabilities and hard-sphere phases are wanted.

        Returns
        -------
        dict of arrays, all the same length as `energies`:
        ``F``, ``dF``, ``G``, ``dG`` (the Coulomb functions and their
        derivatives with respect to rho), ``P`` (penetrability), ``S`` (shift
        function) and ``delta_hs`` (hard-sphere phase shift, radians).

        The values follow the run's own configuration, so the same call returns
        the accurate Coulomb routine's answer, GSL's (``use_gsl_coul=True`` at
        construction), or the Numerov solution through a nuclear potential (the
        hybrid model), and comparing them is how one sees what those options do.
        """
        e = np.asarray(energies, float).ravel()
        req = np.concatenate([[pair, L, radius, e.size], e])
        resp = np.asarray(self.sess.coulomb_functions(req), dtype=float)
        n = int(round(resp[0]))
        block = resp[1:].reshape(n, 7)
        keys = ("F", "dF", "G", "dG", "P", "S", "delta_hs")
        out = {k: block[:, i] for i, k in enumerate(keys)}
        out["energy"] = e[:n]
        return out

    def ec_integrals(self, pair, energies):
        """External-capture radial integrals on an energy grid.

        Every external-capture pathway the compound nucleus generates from this
        entrance pair is evaluated at every energy.  Returns a list of dicts,
        one per pathway, each carrying its quantum numbers (``li``, ``lf``,
        ``si``, ``sf``, ``multipolarity``, ``radiation``) and the complex
        integral as ``value``.

        These integrals are the most expensive part of a capture calculation --
        which is why AZURE2 caches them.  Asking for them twice and watching
        :meth:`cache_stats` is the direct way to see that.
        """
        e = np.asarray(energies, float).ravel()
        req = np.concatenate([[pair, e.size], e])
        resp = np.asarray(self.sess.ec_integrals(req), dtype=float)
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
            })
        return out

    def recalculate_external_capture(self):
        """Force a from-scratch recomputation of every external-capture
        integral.

        AZURE2 caches these in ``output/intEC.dat`` (data mode) /
        ``output/intEC.extrap`` (extrapolation mode) and reads them back on the
        next session start.  The cache is keyed on the *grid*, not on the
        segment selection, so after you add or remove data segments (or change
        the channel radius) the cached integrals belong to other energies and
        must not be reused.  This clears AZURE2's reuse flag and recomputes
        them in the running session.

        For a file you are about to run from scratch, deleting the cache files
        (or :meth:`pyazr.AzrModel.set_output_dir` to a fresh directory) is the
        equivalent and usually simpler.
        """
        if not self.sess.calculate_external_capture():
            raise RuntimeError(
                "AZURE2 could not recompute the external-capture integrals; "
                "check that output/ exists and is writable.")
        return self

    def cache_stats(self):
        """Coulomb-function cache counters, summed over threads.

        Returns a dict with ``queries``, ``hits``, ``hit_rate``, ``entries``,
        ``keys``, ``disabled_keys`` and ``threads``.  ``disabled_keys`` counts
        the keys that gave up on their memo because too few of their entries
        were being asked for twice -- which is what happens when a free energy
        shift moves every point energy at every iteration.
        """
        resp = np.asarray(self.sess.cache_stats(), dtype=float)
        q, h = float(resp[0]), float(resp[1])
        return {
            "queries": int(q),
            "hits": int(h),
            "hit_rate": (h / q) if q else 0.0,
            "entries": int(resp[2]),
            "keys": int(resp[3]),
            "disabled_keys": int(resp[4]),
            "threads": int(resp[5]),
        }

    # -- calculations ---------------------------------------------------------

    def calculate_excitation_energy(self, params):
        """Compound-nucleus excitation energy per segment -- the one axis every entrance pair shares."""
        s = self.sess
        nsegments = int(s.update_segments(params))
        return [s.calculated_excitation_energies(i) for i in range(nsegments)]

    def calculate_angles(self, params):
        """Per-segment angles of the calculated points, one array per segment.

        The companion to :meth:`calculate_energies`: for a differential segment
        the returned angles are AZURE2's own (center-of-mass) values, which
        differ from the lab angles declared in the ``.azr`` file.
        """
        s = self.sess
        nsegments = int(s.update_segments(params))
        return [s.calculated_angles(i) for i in range(nsegments)]

    @staticmethod
    def _unpack_angular_dists(flat):
        """Split the self-describing angular-distribution frame into per-point rows.

        The frame format is a count followed by that many Legendre coefficients,
        repeated once per point, so segments whose points carry different orders
        survive the round trip.
        """
        import numpy as _np
        rows, i, n = [], 0, len(flat)
        while i < n:
            count = int(round(float(flat[i]))); i += 1
            if count < 0 or i + count > n:      # malformed frame; stop rather than guess
                break
            rows.append(_np.asarray(flat[i:i + count], dtype=float))
            i += count
        return rows

    def calculate_analyzing_power(self, params):
        """Vector analyzing power A_y per segment, for analyzing-power segments.

        A_y is returned in the slot a cross section would occupy, because
        AZURE2 treats it as another observable rather than a separate quantity:
        declare a segment with ``observable="analyzing-power"`` and it comes
        back from :meth:`calculate` like any other. This is a named alias for
        that, so the intent is visible at the call site.

        Segments that are not analyzing-power segments return their cross
        section as usual, so check what you asked for.

        Defined for a spin-1/2 projectile on a spin-0 target, in the Madison
        convention with y along k_in x k_out. Bounded by one in magnitude, zero
        in the pure-Coulomb limit, and zero at 0 and 180 degrees.
        """
        return self.calculate(params)

    def calculate_analyzing_power_rwa(self, params):
        """:meth:`calculate_analyzing_power` from reduced-width amplitudes."""
        return self.calculate_rwa(params)

    def calculate_angular_dists(self, params):
        r"""Legendre coefficients of the angular distribution, per segment.

        Returns one entry per segment; each is a list with one array of
        coefficients per calculated point. Points that do not belong to an
        angular-distribution segment give an empty array, so the outer shape
        always matches :meth:`calculate_energies`.

        AZURE2 computes these only for segments declared as angular
        distributions -- ``observable="angular-distribution"`` with an ``order``
        in ``<segmentsTest>``. Everything else returns empty arrays. Use
        :meth:`angular_dist_at` to obtain them at one chosen energy.

        The coefficients are the :math:`a_k` of

        .. math:: W(\theta) = \sum_k a_k P_k(\cos\theta)

        normalised as AZURE2 writes them into ``AZUREOut_*`` files.
        """
        s = self.sess
        nsegments = int(s.update_segments(params))
        return [self._unpack_angular_dists(s.calculated_angular_dists(i))
                for i in range(nsegments)]

    def calculate_angular_dists_rwa(self, params):
        """:meth:`calculate_angular_dists` from reduced-width amplitudes."""
        s = self.sess
        nsegments = int(s.update_segments_rwa(params))
        return [self._unpack_angular_dists(s.calculated_angular_dists(i))
                for i in range(nsegments)]

    def calculate(self, params):
        """Calculated observable per segment, from a physical parameter vector. Cross sections are centre-of-mass, in b or b/sr."""
        s = self.sess
        nsegments = int(s.update_segments(params))
        return [s.calculated_segments(i) for i in range(nsegments)]

    def calculate_rwa(self, params):
        """As calculate, from a free RWA vector. This is the usual entry point."""
        s = self.sess
        nsegments = int(s.update_segments_rwa(params))
        return [s.calculated_segments(i) for i in range(nsegments)]

    def calculate_all_rwa(self, params):
        """As calculate_rwa, taking every parameter rather than the free ones."""
        s = self.sess
        nsegments = int(s.update_segments_all_rwa(params))
        return [s.calculated_segments(i) for i in range(nsegments)]

    def calculate_energies(self, params):
        """Centre-of-mass energies of the calculated points, per segment."""
        s = self.sess
        nsegments = int(s.update_segments(params))
        return [s.calculated_energies(i) for i in range(nsegments)]

    def calculate_sfactor(self, params):
        """Astrophysical S-factor per segment (MeV b), from a physical parameter vector."""
        s = self.sess
        nsegments = int(s.update_segments(params))
        segments = [s.calculated_segments(i) for i in range(nsegments)]
        conv = [s.calculated_conv(i) for i in range(nsegments)]
        for i in range(nsegments):
            segments[i] = segments[i] * conv[i]
        return segments

    def calculate_sfactor_rwa(self, params):
        """As calculate_sfactor, from a free RWA vector."""
        s = self.sess
        nsegments = int(s.update_segments_rwa(params))
        segments = [s.calculated_segments(i) for i in range(nsegments)]
        conv = [s.calculated_conv(i) for i in range(nsegments)]
        for i in range(nsegments):
            segments[i] = segments[i] * conv[i]
        return segments

    # -- modes ----------------------------------------------------------------

    def _set_mode(self, extrap):
        """Switch the engine to data or extrapolation mode and re-initialize."""
        if extrap:
            self.sess.set_extrap()
        else:
            self.sess.set_data()
        self.sess.initialize()
        self.mode = "extrap" if extrap else "data"

    def extrap_mode(self):
        """Switch to the <segmentsTest> grids. Re-initializes the engine."""
        self._set_mode(True)

    def data_mode(self):
        """Switch back to the <segmentsData> segments, which is what chi-squared uses. Re-initializes the engine."""
        self._set_mode(False)
