"""The session: one project set up by the engine in this interpreter (:class:`Session`).

A :class:`Session` owns one ``_azure2.Session`` -- the AZUREAPI C++ object bound into Python via pybind11 -- loaded
from a single ``.azr`` file.  There are no subprocesses, no sockets and no instance pools: the engine runs inside the
interpreter, so a session is one Python object.  Several sessions in one interpreter share nothing and run side by
side; calls on one session from several threads are safe: every call that enters the engine holds the session's lock
(:mod:`pyazr.locking`) and they run one after the other.
"""

import copy
import os
import re
import shutil
import sys
import tempfile
import threading

import numpy as np
from typing import NamedTuple

try:
    from . import _azure2
except ImportError as err:                            # pragma: no cover
    raise ImportError(
        "the AZURE2 engine module (pyazr._azure2) is not built.  It is a C++ "
        "extension, not pure Python: configure with CMake (USE_API=ON, the "
        "default) and build -- the module lands in pyazr/, or in "
        "PYAZR_OUTPUT_DIR when that is set; then export PYAZR_ENGINE_DIR=<it>."
    ) from err

from .errors import ParameterError, ParameterKeyError, close_names, listed
from . import session_checks, session_fit, session_options, session_physics, session_thm
from .locking import engine_call
from .parameters import Pair, PairSet, Parameter, ParameterSet, read_parameter_file
from .project import require_format3
from .potentials import _HybridPotentials


def _config_output_dir(path):
    """The output directory of a ``.azr``, parsed as Config::ReadConfigFile
    does: the line after the first one of ``<config>``, cut at its last
    ``#`` and trimmed ("" if there is none)."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    try:
        i = next(k for k, line in enumerate(lines) if line.strip() == "<config>")  # a whole trimmed line (P2)
    except StopIteration:
        return ""
    if i + 2 >= len(lines):
        return ""
    line = lines[i + 2]
    if "#" in line:
        line = line[:line.rindex("#")]
    return line.strip(" \n\t\r")

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
        """J^pi of the channel (J = s for l = 0), e.g. ``"1/2+"``."""
        j = self.s
        text = f"{int(round(2 * j))}/2" if abs(j - round(j)) > 1e-6 else f"{int(round(j))}"
        return text + ("+" if self.parity > 0 else "-")


class Session(
        _HybridPotentials,
        session_physics.SessionPhysics,
        session_thm.SessionTHM, session_fit.SessionFit,
):
    """One calculation session: a format-3 project set up by the engine in this interpreter (options: see __init__).

    Parameters: :attr:`params_rwa` (the free vector, start values, engine order), :attr:`parameters` (what each entry
    is), :attr:`pairs`, :meth:`physical_levels`.  Data: :attr:`datasets`, :attr:`energies`, :attr:`cross`.  At a
    free vector x: :meth:`evaluate` (every number at x: model values, chi-squared, the objective's terms) and
    :meth:`objective` (the fit's objective as one callable: value, gradient, rows, Jacobian, log-probability).
    Files: :meth:`write_output_files`, :meth:`save_fit`.
    :meth:`describe` (plain data), :meth:`close` (or a ``with`` block).  Errors: :mod:`pyazr.errors`.
    """

    #: Energy-shift penalties are measured from 0, (dE/sigma)^2, not from the
    #: .azr value (2026-10-25).  rmfit reads this to match the engine it runs.
    SHIFT_PENALTY_FROM_ZERO = True

    @session_options.records_given
    def __init__(self, file, *, cwd=None, data_mode=True, use_brune=True,
                 ignore_externals=True, transform=True,
                 use_long_wavelength=True, use_gsl_coul=False, use_rmc=False,
                 use_park=False, use_coulomb_table=True, coulomb_table=None,
                 replace_grids=True, use_ec_table=True, ec_convection_current=True):
        """Build the R-matrix engine for ``file`` in this interpreter.

        Parameters
        ----------
        file : the ``.azr`` configuration file.
        cwd : the directory the model's relative paths (``output/``,
            ``checks/``, data, covariance and table files) are read and written
            in; defaults to the ``.azr`` file's own directory -- the convention
            the file itself is written against.  The engine resolves them
            against it (Config::baseDirectory); the process's working
            directory is neither read nor changed, so sessions in different
            directories can run on different threads at once.
        data_mode : evaluate the ``<segmentsData>`` blocks (the startup state);
            ``False`` starts in extrapolation mode (``<segmentsTest>``).
        use_brune, ignore_externals, transform, use_long_wavelength :
            R-matrix formalism options, matching the CLI's ``--use-brune``,
            ``--ignore-externals``, ``--no-transform`` and ``--no-long-wavelength``
            (the ``use_*`` options are ON by default, ``transform`` too).
            Every option below except ``cwd``, ``data_mode`` and ``coulomb_table`` is also a ``<config>`` key of the
            project under its name (``use_brune=0``; ``Project.set_config(use_brune=False)``): a keyword given here
            wins over the key, the key over the default shown; a keyword left out is not given.
            :attr:`options` holds the values in effect, :attr:`option_sources` where each came from
            (``"default"``, ``"project"``, ``"keyword"``).
        use_gsl_coul : use GSL's Coulomb functions instead of AZURE2's own.
        use_coulomb_table : take P, S and the hard-sphere phase of charged
            channels from the universal Coulomb table (the default; CLI
            ``--no-coulomb-table`` turns it off), falling back to the exact
            routines outside it or when there is no table file.
        coulomb_table : the table file; ``None`` = ``$AZURE2_COULOMB_TABLE``,
            else ``~/.azure2/coulomb_table.bin``.
        replace_grids : lay the adaptive sub-point grids of target effects out
            again on every evaluation, for the current levels and energy
            shifts (the default; CLI ``--no-grid-replacement``).  ``False``
            keeps the grids placed from the file's levels.
        use_ec_table : external-capture integrals from per-term energy tables,
            computed at nodes and interpolated (the default; CLI
            ``--no-ec-table``).  Shifted points are then evaluated at their
            shifted energies.  ``False`` integrates at every energy and
            interpolates shifted amplitudes between the segment's points.
        ec_convection_current : with ``use_long_wavelength=False``, external
            capture carries the convection-current part of the electric
            multipole operator as well as its charge (extended Siegert) part
            (the default; CLI ``--no-ec-convection-current`` turns it off).
            No effect in the long-wavelength approximation.
        use_rmc : Reich-Moore approximation for capture (neutron capture only), CLI ``--use-rmc``.  The program
            then drops ``--use-brune`` and ``--use-park``; a session does not: pass ``use_brune=False`` with it.
        use_park : use Park's level-dependent boundary conditions (Phys. Rev.
            C 104, 064612), the CLI's ``--use-park``.  The rwa vector then
            holds the *observed* reduced width amplitudes,
            ``Gamma = 2 P gamma**2`` with no shift-derivative factor.  Switches
            ``use_brune`` on.  Analytic derivatives are available; the engine
            adds a penalty to chi-squared for any level whose overlap
            ``J = 1 - sum gamma**2 dS/dE`` is not positive (:meth:`park_norms`).

        Raises
        ------
        pyazr.MissingFile
            (a ``FileNotFoundError``) if the file or ``cwd`` does not exist.
        pyazr.FormatError
            if the file is not a format-3 project (the message names the conversion).
        pyazr.AZURE2Error
            (``pyazr.EngineError``) if the model will not initialize.

        Warns
        -----
        pyazr.EngineWarning
            for an engine notice met while setting up (an external-capture amplitude file that cannot be written).
        """
        self._sess = None
        self._engine_lock = threading.RLock()  # every call that enters the engine holds it (pyazr.locking)
        self._azr_cache = {}
        require_format3(file)  # pyazr reads format 3 only: the engine's message names the conversion
        self.file = os.path.abspath(file)
        self.cwd = os.path.abspath(cwd if cwd is not None else (
            os.path.dirname(self.file) or "."))
        if not os.path.isdir(self.cwd):
            raise session_checks.missing_directory(self.cwd, self.file)

        # The keywords this call named win over the project's <config> keys, which win over the defaults (the engine's
        # rule); self.options holds the values in effect, so that a worker process can rebuild this same model
        # (pyazr.bands.sensitivities with nprocs > 1), and self.option_sources where each came from.
        given = {k: v for k, v in locals().items() if k in self._given_options}
        opts = session_options.runtime_options(_azure2, given, data_mode, coulomb_table, self.cwd)
        self._sess = _azure2.Session(self.file, opts)
        self.options, self.option_sources = session_options.in_use(self._sess, data_mode, coulomb_table)

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
    def _engine(self):
        """The binding's _azure2.Session, private: no lock, no promise of stability (used under the session's lock)."""
        if self._sess is None:
            raise session_checks.session_closed(self.file, sys._getframe(1).f_code.co_name)
        return self._sess

    @engine_call
    def close(self):
        """Release the engine.

        The C++ model is several MB per session and is freed here rather than
        whenever the garbage collector gets to it -- which matters in a loop
        over model variants, where the collector may hold two at once.
        """
        self._sess = None

    @engine_call
    def describe(self):
        """What this session holds, as plain data: project, options and option_sources (where each came from),
        Coulomb table, pairs, levels, data sets, extrapolations, parameters (their records: name, kind, unit, x,
        free, bounds, owner fields).
        The keys: :func:`pyazr.described.describe_session`."""
        from .described import describe_session
        return describe_session(self)

    def is_alive(self):
        """Is the engine still open? False once close() has run."""
        return self._sess is not None

    @property
    def output_dir(self):
        """The project's output directory, as the engine resolves it.

        The second line of ``<config>`` (up to its ``#`` comment, trimmed),
        relative to :attr:`cwd` unless absolute -- the directory AZURE2
        writes ``AZUREOut_*``, ``chiSquared.out`` and ``param.sav`` to.  It is
        read from the file each time, so a project whose output lives
        elsewhere than ``output/`` is reported as such.
        """
        outdir = _config_output_dir(self.file)
        return os.path.normpath(os.path.join(self.cwd, outdir)) if outdir else self.cwd

    @engine_call
    def configure(self):
        """Re-read the parameter, pair and dataset metadata from the engine, after anything that rebuilt the model."""
        s = self._engine
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
        self._azr_cache = {}
        session_checks.warn_engine_notices(self._sess)

    # -- the .azr, parsed once per version of the file ------------------------

    def _azr_stamp(self):
        st = os.stat(self.file)
        return (st.st_mtime_ns, st.st_ctime_ns, st.st_size, st.st_ino)

    def _parsed(self, kind):
        """The project's records as the engine describes them, as ``kind`` -- ``"datasets"``
        (:class:`~pyazr.described.DescribedDatasets`), ``"active"`` (its active data sets, checked against the
        engine's count) or ``"extrapolations"`` (:class:`~pyazr.described.DescribedExtrapolations`).

        Parsed once and kept until the file changes on disk (modification time,
        size or inode) or :meth:`configure` runs: an evaluation reads it, and a parse
        of the whole file per call costs 83 ms against 8 ms on a model on NFS.

        The objects are shared: internal read-only use only.  The public
        properties (:attr:`datasets`, :attr:`active_datasets`,
        :attr:`extrapolations`) hand out copies, so a caller who edits what they
        get cannot change the session's view.
        """
        stamp = self._azr_stamp()
        cache = self.__dict__.setdefault("_azr_cache", {})
        if cache.get("stamp") != stamp:
            cache.clear()
            cache["stamp"] = stamp
        if kind not in cache:
            if kind in ("datasets", "extrapolations"):  # the engine's description
                from .described import session_records
                cache[kind] = session_records(kind, self.file, self.cwd)
            elif kind == "active":
                active = self._parsed("datasets").active
                if len(active) != self.nsegments:
                    raise RuntimeError(
                        f"{len(active)} active data sets in {self.file} but "
                        f"the engine reports {self.nsegments} segments.")
                cache[kind] = active
            else:
                raise ValueError(kind)
        return cache[kind]

    # -- parameter metadata ---------------------------------------------------

    @property
    @engine_call
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

    @property
    @engine_call
    def n_rmatrix(self):
        """How many leading entries of the free vector are R-matrix parameters.

        The vector is laid out with the level energies and reduced widths
        first, then the normalizations and energy shifts.  ``transform_rwa``
        accepts either the whole thing or just this leading block, so this is
        the slice point:

        >>> physical = m.transform_rwa(x[:m.n_rmatrix])
        """
        idx = [p.free_index for p in self.parameters
               if p.kind in ("energy", "width") and not p.fixed
               and p.free_index is not None]
        return 1 + max(idx) if idx else 0

    @engine_call
    def refresh_parameters(self):
        """Re-fetch and rebuild the cached :attr:`parameters`."""
        self._parameters = self._build_parameters()
        return self._parameters

    @property
    @engine_call
    def pairs(self):
        """A :class:`PairSet` describing every particle pair (channel).

        Each :class:`Pair` carries the two constituents' spins and parities and
        a flag for the reaction entrance pair.  A width parameter's ``pair``
        attribute is the :attr:`Pair.number`.  Built lazily and cached.
        """
        if self._pairs is None:
            self._pairs = self._build_pairs()
        return self._pairs

    @engine_call
    def refresh_pairs(self):
        """Re-fetch and rebuild the cached :attr:`pairs`."""
        self._pairs = self._build_pairs()
        return self._pairs

    # -- channel radius -------------------------------------------------------

    @engine_call
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
        the safer pattern) use :meth:`pyazr.Project.set_pair` (``radius=``) and
        write a new file.

        The radius stays set for the life of the session, through later
        rebuilds: another pair's ``set_channel_radius``, a nuclear-potential
        change and :meth:`extrap_mode`/:meth:`data_mode` keep it.  (Before
        2026-11-02 each rebuild re-read the ``.azr`` with only its own call's
        radius, so setting pair 2 put pair 1 back to its file value, silently.)
        The radius read back from the rebuilt model is checked against the
        request; a mismatch raises ``RuntimeError``.
        """
        pair = int(pair)
        if not any(p.number == pair for p in self.pairs):
            raise ParameterKeyError(f"set_channel_radius: no pair {pair} in this model (have "
                                    f"{[p.number for p in self.pairs]}); Session.pairs lists the pairs.")
        if not radius > 0:
            raise ParameterError(f"set_channel_radius: the channel radius of pair {pair} must be a positive length "
                                 f"in fm, got {radius}; Session.pairs.by_number({pair}).channel_radius is the "
                                 f"current one.")
        if not self._engine.set_radius(pair, float(radius)):
            raise RuntimeError(
                f"AZURE2 could not rebuild with pair {pair} at "
                f"{radius} fm; the session is no longer usable.")
        self.configure()
        got = self.pairs.by_number(pair).channel_radius
        if abs(got - float(radius)) > 1e-9 * float(radius):
            raise RuntimeError(
                f"pair {pair}: asked for a channel radius of {radius} fm, the "
                f"rebuilt model has {got} fm.")
        return got

    def _build_pairs(self):
        return PairSet(Pair(r) for r in self._engine.pair_records())

    @property
    @engine_call
    def level_scheme(self):
        """A structured, printable :class:`~pyazr.scheme.LevelScheme`.

        Groups the model the way it reads physically -- particle pairs, then
        J-groups, then levels and their channels (L, S, radiation type, partial
        width, fixed flag, Wigner limit).  ``print(azr.level_scheme)`` gives a
        human-readable overview.  Read-only; to add/remove levels and write the
        result to a file see :meth:`pyazr.Project.add_level` and the other level edits of :class:`pyazr.Project`.
        """
        from .scheme import LevelScheme
        return LevelScheme.from_azr(self)

    @property
    @engine_call
    def datasets(self):
        """Per-segment dataset provenance (a :class:`~pyazr.described.DescribedDatasets`, the engine's
        description): for each data set, the data file it came from, the reaction channel (entrance/exit pairs), energy /
        angle range, observable type, and normalization systematic error.
        ``print(azr.datasets.table())`` gives an overview; ``azr.datasets
        .sys_errors()`` returns the per-segment systematics the fits use.
        A copy of the session's parse of the file (see :meth:`_parsed`).
        """
        return copy.deepcopy(self._parsed("datasets"))

    @property
    @engine_call
    def extrapolations(self):
        """Per-segment extrapolation grids (:class:`~pyazr.described.DescribedExtrapolations`).

        For each extrapolation record, the reaction channel, the energy / angle grid, and the observable
        type.  These describe the segments AZURE2 reports in extrapolation mode
        (:meth:`extrap_mode`), in the same order, so segment ``i`` of this set
        corresponds to ``evaluate(x).data[i]`` there.
        A copy of the session's parse of the file (see :meth:`_parsed`).
        """
        return copy.deepcopy(self._parsed("extrapolations"))

    def _build_parameters(self):
        """The table of the engine's parameter records (built once per layout; :meth:`configure` and a mode switch
        drop it)."""
        rows = [Parameter(r) for r in self._engine.parameter_records()]
        return ParameterSet(rows, x0=[p.start for p in rows if not p.fixed], basis=self.basis)

    # -- physical levels (turn a resonance on/off) ----------------------------

    @engine_call
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
        have = listed([str(k) for k in keys])
        if level is not None:
            if level in keys:
                return level
            raise ParameterKeyError(f"{level!r} is not a level of this model (its levels: {have}); "
                                    f"Session.physical_levels() gives the keys.")
        hits = [k for k in keys
                if (jpi is None or k.jpi == jpi)
                and (energy is None or (k.energy is not None
                                        and abs(k.energy - energy) <= tol))]
        if not hits:
            raise ParameterKeyError(f"no level matches jpi={jpi} energy={energy} (MeV, within {tol}); the levels: "
                                    f"{have}; Session.physical_levels() gives the keys.")
        if len(hits) > 1:
            raise ParameterKeyError(f"ambiguous selector jpi={jpi} energy={energy}: "
                                    f"matches {[str(h) for h in hits]}; give energy= (MeV) as well, or level= one "
                                    f"key of Session.physical_levels().")
        return hits[0]

    @engine_call
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

    @engine_call
    def without_level(self, params, level=None, jpi=None, energy=None,
                      tol=1e-2):
        """A copy of the free-parameter vector with one level switched OFF.

        The level's reduced widths are set to zero, decoupling it from every
        channel; all other parameters are untouched.  ``params`` is a free
        (``params_rwa``-order) vector.  Use with :meth:`evaluate` /
        :meth:`objective` to see the model *without* that resonance.
        """
        import numpy as _np
        out = _np.array(params, dtype=float)
        out[self.level_free_width_indices(level, jpi, energy, tol)] = 0.0
        return out

    @engine_call
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

    @engine_call
    def norm_indices(self):
        """Free-vector positions of the normalization parameters."""
        return self._engine.normalization_indices()

    @engine_call
    def shift_indices(self):
        """Free-vector positions of the energy-shift parameters."""
        return self._engine.energy_shift_indices()

    @engine_call
    def energy_indices(self):
        """Packed indices of the free R-matrix level-energy parameters.

        Returns the positions of the level-energy parameters within
        ``params_rwa`` (the non-fixed parameter vector), the same convention
        as :meth:`norm_indices` and :meth:`shift_indices`, from the parameters'
        records (``kind == "energy"``).
        """
        return [p.free_index for p in self.parameters
                if p.kind == "energy" and not p.fixed]

    # -- param.sav helpers ----------------------------------------------------

    # param.sav lives beside the model, so these read and write from the
    # session's directory rather than wherever the caller happens to be.

    @staticmethod
    def read_sav(path):
        """``{name: value}`` from a ``param.par`` / ``param.sav`` / ``param.fit``, by NAME
        (:func:`pyazr.parameters.read_parameter_file`)."""
        return read_parameter_file(path)[1]

    @engine_call
    def full_rwa_from_sav(self, path):
        """The model's full RWA vector (fixed parameters included, ``p.index``
        order) with the values a parameter file names applied.

        A name the file lacks keeps the session's current value, exactly as
        AZURE2 does on read; a file with no name in common with the model is
        refused rather than applied to nothing.

        This only *reads*.  The free-vector calls (``evaluate(x)``,
        ``update_segments_rwa(x)``, ``transform_rwa(x)``, ...) take every FIXED
        parameter from the session, i.e. from the ``.azr``, so passing just the
        free slice of this vector drops the file's fixed level parameters --
        which the CLI does apply.  To evaluate "the ``.azr`` + this ``.sav``"
        as the CLI does, use :meth:`load_sav`.
        """
        vals = self.read_sav(path)
        full = np.asarray(self._engine.params_all_rwa(), float).copy()
        hit = 0
        for p in self.parameters:
            if p.name in vals and 0 <= p.index < full.size:
                full[p.index] = vals[p.name]
                hit += 1
        if hit == 0:
            raise ValueError(f"{path}: no parameter name matches this model "
                             f"(is it a file from another project?)")
        return full

    @staticmethod
    def read_sav_basis(path):
        """The ``#parametrization`` tag of a parameter file (0 standard,
        1 Brune, 2 Park), or None for a file written before the tag existed."""
        return read_parameter_file(path)[0]

    @engine_call
    def load_sav(self, path):
        """Load a parameter file into the session the way the CLI does, and
        return the free RWA vector to evaluate at.

        The CLI's external parameter file (mode 1/2/3 prompt; ``AZUREParams::
        ReadUserParameters``) sets every parameter it names, FIXED level
        energies and widths included.  The free-vector calls of this class
        take the fixed values from the session instead, so
        ``x = full_rwa_from_sav(p)[~fixed]`` followed by
        ``evaluate(x)`` evaluates the ``.azr``'s fixed parameters,
        not the file's.  When they differ the two are different models:
        11B+alpha (10-12-26) had 10 fixed widths off by up to 1.6 %, and pyazr
        gave chi2 9925.03 against the CLI's 9868.04 at the same file, with no
        warning (key 72 alone +58).

        This pushes the whole vector (``full_rwa_from_sav``) into the engine,
        so every later free-vector call -- chi2, segments, Jacobian,
        ``transform_rwa``, ``save_fit`` -- uses the file's fixed values, and
        sets :attr:`params_rwa` to the file's free values, which it returns.

        Fixed segment norms and shifts are carried along but the engine never
        reads them (only varied ones are filled from parameters), so they stay
        at the ``.azr`` value, as in the CLI.  Unlike the CLI, the free/fixed
        split is NOT changed: the CLI frees a parameter fixed at 0 that the
        file gives a nonzero value, and fixes a free one the file sets to 0.
        Values are applied either way, so a calculation agrees; a fit started
        from here keeps the ``.azr``'s split.  A file whose
        ``#parametrization`` tag differs from this session's basis is refused
        (the CLI converts; this does not).

        Switching mode (:meth:`extrap_mode` / :meth:`data_mode`) re-initializes
        the engine from the ``.azr`` and discards what this loaded: call it
        again after the switch.  Names are matched against the engine's live
        parameter list, so it works in either mode, as the CLI's reader does.
        In extrapolation mode the file's ``segment_*`` entries (data-segment
        normalizations and energy shifts) are not applied: they belong to data
        segment k, and until 2026-10-26 they moved test segment k's energy grid.
        """
        tag = self.read_sav_basis(path)
        if tag is not None and tag != self.basis:
            raise ValueError(
                f"{path} holds basis {tag} amplitudes (0 standard, 1 Brune, 2 Park) "
                f"but this session uses basis {self.basis}; open the model with the "
                f"matching use_brune/use_park, or convert the file with the CLI.")
        # The engine's LIVE bookkeeping, matched by name as the CLI does: the
        # cached `parameters` describe the mode the session was first built
        # in, and extrapolation mode has a shorter vector.
        s = self._engine
        fixed = np.asarray(s.params_fixed(), bool)
        full = np.asarray(s.params_all_rwa(), float).copy()
        vals = self.read_sav(path)
        hit = 0
        skipped = 0
        segment_field = {k["key"]: k["segment_field"] for k in _azure2.parameter_kinds()}
        records = s.parameter_records()
        self._match_parameter_file_names(path, list(vals), fixed, records, segment_field)
        for i in range(full.size):
            name = s.params_names(i)
            if name in vals:
                if self.mode == "extrap" and segment_field[records[i]["kind"]]:
                    # A test segment carries no normalization or energy shift of its own; the
                    # file's entries for segment parameters belong to DATA segment k and would move
                    # test segment k's grid (audit 2026-10-25 item 12; the CLI skips them too).
                    skipped += 1
                    continue
                full[i] = vals[name]
                hit += 1
        if skipped:
            import warnings
            warnings.warn(f"{path}: {skipped} segment_* entries (data-segment normalizations and "
                          f"energy shifts) not applied in extrapolation mode; a test segment's grid "
                          f"is the one its line asks for.")
        if hit == 0:
            raise ValueError(f"{path}: no parameter name matches this model "
                             f"(is it a file from another project?)")
        setter = getattr(s, "set_params_all_rwa", None)
        if setter is not None:
            setter(full)
        else:
            # Modules built before set_params_all_rwa: the same effect, at the
            # cost of one forward pass.
            self._engine.update_segments_all_rwa(full)
        self.fixed_params = [bool(f) for f in fixed]
        self.params_rwa = [float(v) for v in full[~fixed]]
        self.params = s.params_values()
        self._parameters = None
        return np.asarray(self.params_rwa, float)

    def _match_parameter_file_names(self, path, in_file, fixed, records, segment_field):
        """The file's names against this run's parameters, by the engine's rule (experiment::MatchParameterFileNames,
        as AZURE2 reads a param.sav): no name of the run raises ValueError with E_PARAMETER_FILE_NO_MATCH, names the
        project lacks or run parameters the file leaves out warn with W_PARAMETER_FILE_NAMES.  Fixed data-set and
        extrapolation nuisances may be left out; names the project gives to another run's parameters (an
        extrapolation's in a data run) are neither applied nor unknown."""
        import json
        import warnings
        s = self._engine
        of_run = [s.params_names(i) for i in range(len(records))]
        optional = [of_run[i] for i, r in enumerate(records) if segment_field[r["kind"]] and fixed[i]]
        with open(self.file, "rb") as f:
            other = list(_azure2.parameter_names(f.read()))
        for text in _azure2.match_parameter_file_names(os.fspath(path), in_file, of_run, optional, other):
            d = json.loads(text)
            message = f"{d.get('code', '')}: {d.get('message', '')}"
            if d.get("severity") == "error":
                raise ValueError(message)
            warnings.warn(message, stacklevel=3)

    @engine_call
    def update_rwa_params_from_sav(self):
        """Reload the run's param.sav (matched by name) into the session --
        fixed level parameters included, see :meth:`load_sav` -- and set
        params_rwa to its free values."""
        self.load_sav(os.path.join(self.cwd, self.output_dir, 'param.sav'))

    @engine_call
    def update_sav_from_rwa_params(self, best):
        """Write a free RWA vector back out to param.sav.new, by name.

        Rows the file already has are updated in place (other columns and the
        tag line kept); a free parameter the file does not list is appended, so
        the result describes the whole model.
        """
        best = np.asarray(best, float).ravel()
        free = {p.name: float(best[p.free_index]) for p in self.parameters
                if not p.fixed and p.free_index is not None and p.free_index < best.size}
        sav = os.path.join(self.output_dir, 'param.sav')
        if not os.path.isfile(sav):
            raise session_checks.missing_parameter_file(sav, "update_sav_from_rwa_params")
        rows, tag, seen = [], None, set()
        with open(sav) as fh:
            for line in fh:
                t = line.split()
                if not t:
                    continue
                if t[0].startswith("#") or t[0] == "parametrization":
                    tag = t
                    continue
                rows.append(t)
        for t in rows:
            if t[0] in free:
                t[1] = repr(free[t[0]])
                seen.add(t[0])
        for name, value in free.items():
            if name not in seen:
                rows.append([name, repr(value), "0.0"])
        with open(sav + '.new', 'w') as f:
            if tag is not None:
                f.write(" ".join(tag) + "\n")
            for t in rows:
                f.write(" ".join(t) + "\n")

    @engine_call
    def save_fit(self, path, x=None, param_sav=True, verify=True, norms="fitted",
                 close_session=False):
        """Snapshot a fit as a ``.azr`` you can reopen, plot, or hand over.

        ``path`` is required and is written to explicitly -- nothing is ever
        saved in place over the model you loaded.

        ``x`` is a free vector in :attr:`params_rwa` order; it defaults to the
        session's current parameters. The conversion to the physical values a
        ``<levels>`` line holds (partial widths in eV, ANCs in fm^-1/2) is done
        here, so a caller cannot forget it.

        **Fixed parameters are written too.** A width marked fixed holds its
        reduced-width amplitude during a fit -- in pyazr as in the CLI's MIGRAD
        and LM -- so its *physical* value moves whenever the level's other
        widths do (the Brune denominator ``1 + sum_c gamma_c^2 dS_c/dE``
        contains them all).  The snapshot carries the physical value the fit
        actually had, so the fixed amplitude reads back unchanged.  (Before,
        the file kept the fixed width as it was read: a fit that tripled a
        level's free width reloaded as a different model, 4 % off in chi2 on
        tests/15N_p_a.)

        Two things a ``<levels>`` block cannot carry, which is why this is not
        just a write of the level values into ``<levels>``:

        - **Normalizations, energy shifts and the other nuisances are not in
          it.** They are values of the data-set and effect records (their
          nuisances: a start value, ``fit`` or ``fixed``, and a prior that
          states its centre).  ``norms`` chooses:

          ``"fitted"`` (default)
              every free parameter that is not a level parameter (norms,
              shift terms, shift-function parameters, abundances, ...) is
              written into the record it comes from, in that record's unit;
              the priors keep the centres they state.  The snapshot run on
              its own gives the fitted data chi-squared *and* the fit's
              priors (both checked by ``verify``).
          ``"nominal"``
              the records keep the loaded file's values (but the abundances,
              written at the fit's values); the fitted norms and shifts are
              in the ``.sav`` only (needs ``param_sav``).  Run with that
              ``.sav`` as the external parameter file (CLI: "use previous
              parameters") to reproduce the fit.

          ``param_sav`` (the default) writes the companion ``<name>.sav``
          beside the file, carrying every parameter, free and fixed.
        - **A written file is not a fit until it reads back as one.** With
          ``verify`` (the default) the result is reopened and checked two
          ways: its physical values transformed back must match, and the
          reduced-width amplitudes it reads for the particle and beta-decay
          channels of every level must be ``x``'s, up to a common sign of
          the level (the standard parametrization's eigenvector convention,
          no change of model). If either disagrees the files are removed and
          this raises, rather than leaving a snapshot that is quietly a
          mixture. The second check exists because the first cannot see a
          value two parameters map onto: the engine returned ``|B|`` for a
          beta-decay amplitude until 2026-10-12, and a snapshot with the
          sign lost passed. See :meth:`_snapshot_mismatch` for why photon
          channels are judged on physical values only.  With
          ``norms="fitted"`` (data mode) the reopened file's data chi-squared
          and its objective (:meth:`evaluate`: ``chi2``, ``objective``) must also equal the fit's
          (rel. 1e-6; 1e-4 in the standard parametrization, whose transform is
          inverted iteratively on reading).

        A ``.azr`` names its data files and its output directory *relative to
        itself*, so a snapshot only runs from the directory the original did.
        Write it beside the model, or move its ``data/`` with it.

        ``close_session`` closes this session once the files are written and
        what the check compares with has been computed, *before* the
        snapshot is reopened: one engine in memory instead of two (a 19F THM
        session is about 220 MB, 390 MB with the verifying one beside it).
        The session is unusable afterwards (:meth:`is_alive` is False); use
        it as the last call of a fit, as ``scripts/thm_model_average.py``
        does.

        Returns ``(azr_path, sav_path_or_None)``.
        """
        if norms not in ("fitted", "nominal"):
            raise ParameterError(f"save_fit: norms must be 'fitted' or 'nominal', not {norms!r} (see "
                                 f"help(pyazr.Session.save_fit)).")
        if norms == "nominal" and not param_sav:
            raise ParameterError("save_fit: norms='nominal' keeps the fitted norms in the .sav "
                                 "only; it needs param_sav=True.")
        x = np.asarray(self.params_rwa if x is None else x, float).ravel()
        free = [i for i, p in enumerate(self.parameters) if not p.fixed]
        if x.size != len(free):
            raise ParameterError(f"save_fit: x has {x.size} entries; the session has {len(free)} free parameters "
                                 f"(Session.params_rwa is such a vector, Session.parameters.free names its entries).")
        # Every parameter, with the free ones from x, and every parameter's
        # physical value, fixed ones included (see above).
        full = self._all_rwa(x)
        want = np.asarray(self.transform_all_rwa(full, include_fixed=True), float)

        path = os.path.abspath(path)
        from .snapshot import write_format3
        write_format3(self, path, x, want, full, norms)

        sav = None
        if param_sav:
            sav = self.parameters.write(os.path.splitext(path)[0] + ".sav", x)

        # What the check compares with, from this session before it may close.
        chi_want = pen_want = None
        if verify and norms == "fitted" and self.mode == "data":
            ev = self.evaluate(x)
            chi_want, pen_want = ev.chi2, ev.objective
        check_options = {k: v for k, v in self.options.items() if k != "data_mode"}
        check_options["data_mode"] = self.mode == "data"
        if close_session:
            self.close()

        if verify:
            # Verify a *copy*, in the session's own directory and with its
            # output redirected. A .azr resolves its data files and output
            # directory relative to itself, so the snapshot only loads where
            # the original did -- and opening it there would otherwise
            # overwrite the real run's output/param.par.
            probe = tmpdir = None
            try:
                tmpdir = tempfile.mkdtemp(prefix=".azr_verify_", dir=self.cwd)
                fd, probe = tempfile.mkstemp(suffix=".azr", dir=self.cwd)
                os.close(fd)
                from .project import Project
                copy = Project.open(path)
                copy.set_config(output_dir=tmpdir.rstrip("/") + "/")
                copy.save(probe)
                with Session(probe, cwd=self.cwd, **check_options) as check:
                    back_full = np.asarray(check._engine.params_all_rwa(), float)
                    got = np.asarray(check.transform_all_rwa(back_full, include_fixed=True), float)
                    back = {p.name: float(back_full[i])
                            for i, p in enumerate(check.parameters)}
                    chi_got = None
                    if chi_want is not None:
                        ev = check.evaluate(check.params_rwa)
                        chi_got, pen_got = ev.chi2, ev.objective
            except Exception as err:
                self._discard(path, sav)
                raise RuntimeError(
                    f"the snapshot was written but would not load back "
                    f"({err}); it has been removed.") from err
            finally:
                if probe and os.path.exists(probe):
                    os.remove(probe)
                if tmpdir and os.path.isdir(tmpdir):
                    shutil.rmtree(tmpdir, ignore_errors=True)
            bad = self._snapshot_mismatch(full, want, got, back)
            if bad:
                self._discard(path, sav)
                raise RuntimeError(
                    "the snapshot does not read back as the fit it was made "
                    f"from ({'; '.join(bad)}); it has been removed. This means "
                    "the .azr and the parameter set do not describe the same "
                    "model.")
            if chi_got is not None:
                # The standard parametrization (basis 0) inverts its transform
                # iteratively on reading (energies to 1e-6 MeV), which leaves
                # ~1e-5 in the chi-squared of a steep spectrum
                # (tests/pyazr/beta_sign_test.py); Brune and Park read back
                # exactly.
                rtol = 1e-4 if self.basis == 0 else 1e-6
                if not np.isclose(chi_got, chi_want, rtol=rtol, atol=1e-9):
                    self._discard(path, sav)
                    raise session_checks.snapshot_mismatch("chi2", chi_got, chi_want, rtol)
                if not np.isclose(pen_got, pen_want, rtol=rtol, atol=1e-9):
                    self._discard(path, sav)
                    raise session_checks.snapshot_mismatch("objective", pen_got, pen_want, rtol)
        return path, sav

    @staticmethod
    def _discard(*paths):
        """Remove a snapshot that failed verification, so nothing bad survives."""
        for f in paths:
            if f and os.path.exists(f):
                os.remove(f)

    def _snapshot_mismatch(self, x, want, got, back):
        """Why a re-read snapshot is not the fit it was made from; ``[]`` if it is.

        ``x``, ``want``, ``got`` and ``back`` cover every parameter, fixed ones
        included (``params_all_rwa`` order; ``back`` by name).  Two
        comparisons, because neither alone is enough:

        - the **physical values** (``transform_rwa`` of both sides), which
          is what the file holds.  A transform that maps two parameters onto
          one value passes this alone: the engine returned ``|B|`` for a
          beta-decay amplitude until 2026-10-12, and a snapshot with the
          sign lost read back as a different model (12C+alpha).
        - the **re-read reduced-width amplitudes** of the particle and
          beta-decay channels, and the level energies, against ``x`` --
          level by level, allowing a common sign of the whole level, which
          the standard parametrization's eigenvector convention may flip and
          which is no change of model.

        Photon channels are left to the first comparison on purpose: their
        amplitude is derived from the written width together with the
        external width of the level, a sum over the level's particle
        amplitudes, so it is only as exact as the written widths are, and the
        physical comparison is the direct one.
        """
        if got.shape != want.shape:
            return [f"{got.size} physical parameters instead of {want.size}"]
        bad = []
        levels = {}
        for i, p in enumerate(self.parameters):
            if p.kind not in ("energy", "width") or i >= want.size:
                continue
            if p.name not in back:
                bad.append(f"{p.name} is not in the re-read file")
                continue
            levels.setdefault((p.jgroup, p.level), []).append(
                (p, float(want[i]), float(got[i]), float(x[i]), back[p.name]))

        def same(pairs, rtol):
            """Equal, or equal with the level's common sign flipped."""
            if not pairs:
                return True
            a = np.array([v[0] for v in pairs])
            b = np.array([v[1] for v in pairs])
            return bool(np.allclose(a, b, rtol=rtol, atol=1e-12)
                        or np.allclose(a, -b, rtol=rtol, atol=1e-12))

        for (jg, lv), items in sorted(levels.items()):
            ok = all(abs(w - g) <= 1e-5 * max(1.0, abs(w))
                     for p, w, g, a, b in items if p.kind == "energy")
            # the written physical widths, every channel (sign carried)
            ok = ok and same([(w, g) for p, w, g, a, b in items
                              if p.kind == "width"], 1e-4)
            # the re-read amplitudes: particle and beta channels only, since
            # a photon channel's amplitude depends on the external-capture
            # integrals the re-read computed afresh
            ok = ok and same([(a, b) for p, w, g, a, b in items
                              if p.kind == "width"
                              and p.radiation_type in ("P", "F", "G")], 1e-3)
            ok = ok and all(abs(a - b) <= 1e-5 * max(1.0, abs(a))
                            for p, w, g, a, b in items if p.kind == "energy")
            if not ok:
                p0 = items[0][0]
                bad.append(f"level {p0.jpi} #{lv} (J-group {jg}, E = {p0.level_energy}) "
                           "re-reads to different parameters")
        return bad

    # -- transforms -----------------------------------------------------------

    @engine_call
    def transform_rwa(self, params):
        """Map a reduced-width-amplitude vector to physical parameters: level energies in MeV, partial widths in eV, ANCs in fm^-1/2."""
        return self._engine.transform_rwa(params)

    @engine_call
    def transform_all_rwa(self, params, include_fixed=False, free_only=None):
        """As transform_rwa, taking every parameter rather than the free ones.

        Returns the physical values of the free parameters, or with
        ``include_fixed`` (equivalently ``free_only=False``) of every
        parameter, in ``params_all_rwa`` order (the level parameters come
        first) -- what a snapshot of a fit writes.  A *fixed* width holds its
        reduced-width amplitude (as in the CLI's fits), so its physical value
        is not constant: it follows the level's other widths through the
        Brune denominator ``1 + sum_c gamma_c^2 dS_c/dE``.  Under Park
        (``use_park``) the amplitude is the observed width and stays put.  A
        channel entered as an amplitude (``gammaIsRWA``) is returned as
        Brune's amplitude in both modes, the convention of the .azr.
        """
        if free_only is not None:
            include_fixed = not free_only
        return self._engine.transform_all_rwa(params, include_fixed=bool(include_fixed))

    def _all_rwa(self, x):
        """The full RWA vector with the free vector ``x`` in its free slots."""
        allx = np.asarray(self._engine.params_all_rwa(), float).copy()
        free = [i for i in range(allx.size) if not self.fixed_params[i]]
        x = np.asarray(x, float).ravel()
        for k, i in enumerate(free[:x.size]):
            allx[i] = x[k]
        return allx

    @engine_call
    def transform_physical(self, params):
        """The inverse of transform_rwa: physical parameters back to reduced-width amplitudes."""
        raise NotImplementedError(
            "AZURE2 exposes no physical->RWA transform; only RWA->physical "
            "is available (transform_rwa / transform_all_rwa).")

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

    @engine_call
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

    # -- output files and norms -----------------------------------------------

    @engine_call
    def write_output_files(self, params=None):
        """Write the run's standard output files, as the CLI does at the end.

        ``AZUREOut_*``, ``chiSquared.out`` and the rest land in the ``.azr``'s
        output directory, created when it does not exist.  Without this a Python
        session could compute everything and still had to re-run the binary to
        get the files a colleague or the GUI expects.

        A forward pass is run first at ``params`` (the current parameters by
        default), so the files describe the parameters you asked about rather
        than whatever was evaluated last.

        Returns the output directory the engine wrote to (:attr:`output_dir`:
        the ``<config>`` one, relative to :attr:`cwd`), not a fixed
        ``cwd/output``.
        """
        x = np.asarray(self.params_rwa if params is None else params, float).ravel()
        os.makedirs(self.output_dir, exist_ok=True)
        if self.mode == "data":
            # chiSquared.out reports the per-segment chi-squared *stored on
            # each segment*, and only the chi-squared evaluation sets it --
            # a bare forward pass leaves it at zero and the file comes out
            # full of zeros.  This populates both that and the point cross
            # sections AZUREOut_* is written from.
            self._engine.calculate_chi2_rwa(x)
        else:
            # No data to compare against; chiSquared.out is meaningless in
            # extrapolation mode, and the forward pass is all AZUREOut_*
            # needs.
            self._engine.update_segments_rwa(x)
        if not self._engine.write_output_files():
            raise RuntimeError("AZURE2 has no data object to write from.")
        session_checks.warn_engine_notices(self._sess, stacklevel=3)
        return self.output_dir

    @property
    @engine_call
    def active_datasets(self):
        """The data segments the engine evaluates, in the order its results
        come back: ``m.energies[i]``, ``m.cross[i]``, ``evaluate(x).data[i]``
        all belong to ``active_datasets[i]``.

        :attr:`datasets` lists every ``<segmentsData>`` line, inactive ones
        included, so ``datasets[i]`` is *not* segment ``i`` of a result on any
        model that has one -- the 13C+alpha archive's model has 25 -- and code
        that indexed it that way dropped the normalization penalty of every
        active segment past the active count.
        """
        return copy.deepcopy(self._parsed("active"))

    # -- the objective, the evaluation and what the engine adds to them ------

    @engine_call
    def prior_covariance_residuals(self, params=None):
        """Whitened residuals r = L^-1 (p - mu) of the prior parameter covariance (``priorCovariance=<file>`` in
        ``<config>``, C = L L^T), one per prior parameter in :attr:`prior_covariance_names` order; empty without
        the key (or in extrapolation mode).  ``sum(r**2)`` is ``evaluate(x).terms["prior covariance"]``, and
        ``objective().residuals(x)`` holds these rows (Jacobian L^-1 on the prior's columns)."""
        x = np.asarray(self.params_rwa if params is None else params, float).ravel()
        if not hasattr(self._engine, "prior_covariance_residuals"):    # an engine older than 2026-11-05
            return np.zeros(0)
        return np.asarray(self._engine.prior_covariance_residuals(x), float)

    @property
    @engine_call
    def prior_covariance_names(self):
        """The prior parameter covariance's parameter names, in its row order, after the engine marginalized the
        listed parameters that are fixed in this model; empty without the ``priorCovariance`` key."""
        if not hasattr(self._engine, "prior_covariance_names"):
            return []
        return list(self._engine.prior_covariance_names())

    @property
    @engine_call
    def nuclei(self):
        """The compound nuclei of the project: ``{(Z, A): PairSet}`` of the
        pairs that form each, in pair-number order (one entry for every project
        written before several nuclei could share a file)."""
        from .parameters import PairSet
        out = {}
        for pair in self.pairs:
            out.setdefault(pair.nucleus, PairSet()).append(pair)
        return out

    @property
    def basis(self):
        """Which amplitudes ``params_rwa`` are: 0 standard (constant boundary
        condition), 1 Brune, 2 Park -- the ``#parametrization`` tag AZURE2 writes
        first in ``param.sav``."""
        if self.options.get("use_park"):
            return 2
        return 1 if self.options.get("use_brune") else 0

    @engine_call
    def scattering_lengths(self, params=None):
        """The s-wave scattering lengths of every particle pair at a free RWA
        vector (default ``params_rwa``): a list of :class:`ScatteringLength`,
        one per allowed channel spin (identical pairs: even s only).

        a0 = -a [M/(x^2 K1^2) - 2 I1/(x^2 K1)] with M = gamma^T A gamma at the
        pair's threshold (Paneru et al., PRC 99, 045807 (2019), Eq. 7), a (1 - M)
        for a neutral pair; the same in every parametrization.  Call it on each
        sample of an MCMC chain for the posterior of a0.  The measured ones of
        the ``scattering_length`` data of ``<modelData>`` are already part of
        the objective (``evaluate(x).terms["model data"]``)."""
        x = np.asarray(self.params_rwa if params is None else params, float)
        flat = np.asarray(self._engine.scattering_lengths_rwa(x), float).reshape(-1, 6)
        return [ScatteringLength(int(r[0]), float(r[1]), int(r[2]), complex(r[3], r[4]), bool(r[5])) for r in flat]

    @engine_call
    def park_norms(self, params=None):
        """Park formalism: the overlap J = 1 - sum_c gamma_c^2 dS_c/dE of every
        R-matrix level, in ``physical_levels()`` order.

        J is the squared norm of the level's basis state (Park 2021, Eq. 28)
        and the factor between the two parametrizations,
        ``gamma_Park = gamma_Brune * sqrt(J)``.  Every real Brune parameter set
        has 0 < J <= 1; a Park parameter set with J <= 0 has no R-matrix
        counterpart (its widths exceed what the channel radii allow) and the
        engine penalizes it (``evaluate(x).terms["Park J > 0"]``).  All ones in Brune mode.
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        return np.asarray(self._engine.park_norms(x), float)

    @engine_call
    def refreeze_covariance(self):
        """Data covariances with ``relvector`` terms: forget the model frozen for them.

        The model m in the multiplicative terms (f m)(f m)^T is frozen at the first
        chi-squared evaluation, so the covariance stays constant through a minimization.
        After this call the next evaluation freezes it again at its parameters: the
        iterative scheme is fit, ``refreeze_covariance()``, evaluate at the best fit,
        fit again.  Returns the number of covariance groups with relvector terms.
        """
        return int(self._engine.refreeze_covariance())

    @engine_call
    def objective(self):
        """The fit's objective as one object for any minimiser or sampler (:class:`pyazr.objective.Objective`):
        ``f(x)`` (the data chi-squared plus every prior and constraint), ``f.value_and_grad``, ``f.grad``,
        ``f.residuals``, ``f.jacobian``, ``f.residuals_and_jacobian``, ``f.log_prob``, ``f.x0``, ``f.names``,
        ``f.bounds``."""
        from .objective import Objective
        return Objective(self)

    @engine_call
    def evaluate(self, x=None):
        """Every number of the session at the free vector ``x`` (default ``params_rwa``) from one engine evaluation:
        a :class:`pyazr.evaluate.Evaluation` (the model of every point, ``data[i].model``, ``.energies``, ``.extra``;
        in data mode also the objective, its terms, the data chi-squared and each data set's).  In extrapolation mode
        the data sets are the extrapolations and every chi-squared field is ``None``."""
        from .evaluate import evaluate
        return evaluate(self, x)

    @engine_call
    def model_gradients(self, params):
        """Analytic ``d(observable)/d(parameter)`` for every calculated point.

        Returns a list with one ``(npoints, ncols)`` array per calculated
        segment, in the order and point ordering of ``evaluate(x).data``.  The
        columns are the free **R-matrix** parameters -- level energies and
        reduced-width amplitudes, in ``params_rwa`` order -- and the THM
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

    # -- the external region, and the caches that make it affordable ----------

    @engine_call
    def coulomb_functions(self, pair, energies, L=0, radius=0.0, J=None):
        """Coulomb wave functions on an energy grid.

        Parameters
        ----------
        pair : particle-pair key (1-based, as in the .azr).
        energies : centre-of-mass energies in MeV.
        L : orbital angular momentum; J : the channel's total J, required (and only allowed) for a spin-orbit pair.
        radius : evaluation radius in fm; 0 means the pair's channel radius,
            which is where penetrabilities and hard-sphere phases are wanted.

        Returns
        -------
        dict of arrays, all the same length as `energies`:
        ``F``, ``dF``, ``G``, ``dG`` (the Coulomb functions and their
        derivatives with respect to rho), ``P`` (penetrability), ``S`` (shift
        function) and ``delta_hs`` (hard-sphere phase shift, radians); for a
        closed channel (E <= 0, eta > 100) S alone, at the channel radius.

        The values follow the run's own configuration, so the same call returns
        the accurate Coulomb routine's answer, GSL's (``use_gsl_coul=True`` at
        construction), or the hybrid model's Numerov solution: compare them.

        Zeros by design: at E <= 0 (a closed channel) all but ``S`` (the shift
        function at the channel radius; 0 at another radius); at an energy
        where the engine's Coulomb routine fails, all.  Raises
        :class:`pyazr.ParameterIndexError` (an ``IndexError``) for an unknown
        pair key, :class:`pyazr.DomainError` (a ``ValueError``) for a photon
        pair or an L that is not a whole number >= 0, ``ValueError`` for a J
        a spin-orbit pair needs or cannot make (E_SO_ codes).
        """
        pair, L = session_checks.channel_request(self.pairs, pair, L)
        e = np.asarray(energies, float).ravel()
        req = np.concatenate([[pair, L, -1 if J is None else round(2 * J), radius, e.size], e])
        resp = np.asarray(self._engine.coulomb_functions(req), dtype=float)
        n = int(round(resp[0]))
        block = resp[1:].reshape(n, 7)
        keys = ("F", "dF", "G", "dG", "P", "S", "delta_hs")
        out = {k: block[:, i] for i, k in enumerate(keys)}
        out["energy"] = e[:n]
        return out

    @engine_call
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
            })
        return out

    @engine_call
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
        (or :meth:`pyazr.Project.set_config` ``output_dir=`` a fresh directory) is the
        equivalent and usually simpler.
        """
        if not self._engine.calculate_external_capture():
            raise RuntimeError(
                "AZURE2 could not recompute the external-capture integrals; "
                "check that output/ exists and is writable.")
        session_checks.warn_engine_notices(self._sess, stacklevel=3)
        return self

    @engine_call
    def cache_stats(self):
        """Input-file counters of the process (and the former Coulomb-memo counters, now always 0).

        Returns a dict with ``queries``, ``hits``, ``hit_rate``, ``entries``,
        ``keys``, ``disabled_keys`` and ``threads``: the counters of the
        process-wide Coulomb-function memo of earlier builds.  The engine has
        no such memo any more (the Coulomb functions are computed at every
        call; docs/upgrade/coulomb.md), so they are 0; the keys are kept so
        that scripts reading them keep working.

        Also ``input_files_read`` (input files read from
        disk), ``input_files_served`` (opens answered from memory) and
        ``input_files_held``, for the whole process: a session reads the .azr and
        its data files once and again only when a file changes on disk, so a
        rebuild normally reads nothing.
        """
        resp = np.asarray(self._engine.cache_stats(), dtype=float)
        q, h = float(resp[0]), float(resp[1])
        out = {
            "queries": int(q),
            "hits": int(h),
            "hit_rate": (h / q) if q else 0.0,
            "entries": int(resp[2]),
            "keys": int(resp[3]),
            "disabled_keys": int(resp[4]),
            "threads": int(resp[5]),
        }
        if len(resp) >= 9:
            out.update(input_files_read=int(resp[6]), input_files_served=int(resp[7]),
                       input_files_held=int(resp[8]))
        return out

    def _mode_free_vector(self, params):
        """A free RWA vector for the engine's current layout.

        In extrapolation mode a data-mode vector (longer: the data segments'
        norms and shifts) is accepted, as before the binding refused vectors
        of the wrong length (dev, 2026-10): its free R-matrix (and abundance)
        block, which comes first and is the same in both modes (see
        :meth:`_sync_layout`), is taken, the tail is the session's.  Any other
        length goes to the engine as given, which refuses it.
        """
        x = np.asarray(params, float).ravel()
        cur = np.asarray(self.params_rwa, float).ravel()
        if self.mode != "extrap" or x.size <= cur.size:
            return params
        head = sum(1 for r in self._engine.parameter_records()
                   if not r["fixed"] and r["kind"] in ("energy", "width", "abundance"))
        y = cur.copy()
        y[:head] = x[:head]
        return y

    # -- modes ----------------------------------------------------------------

    def _set_mode(self, extrap):
        """Switch the engine to data or extrapolation mode and re-initialize."""
        if extrap:
            self._engine.set_extrap()
        else:
            self._engine.set_data()
        self._engine.initialize()
        self.mode = "extrap" if extrap else "data"
        self._sync_layout()
        session_checks.warn_engine_notices(self._sess, stacklevel=3)

    def _sync_layout(self):
        """Follow the engine's parameter layout after a mode switch.

        The engine lists one normalization (if varied) and two energy-shift
        parameters per segment of the *current* segment set, so the
        ``<segmentsTest>`` grids of extrapolation mode give a different
        parameter count from the ``<segmentsData>`` segments (e.g. 220 against
        212 on a 12C+12C model with eight test segments).  The cached
        ``fixed_params`` used to keep the data-mode layout, so the parameter
        table (:attr:`parameters`), :meth:`_all_rwa` and :meth:`save_fit` mixed
        the two (220 parameter records against 212 fixed flags).  Refresh the cached layout from the engine whenever it
        changed, and carry the free R-matrix block of ``params_rwa`` /
        ``params`` over -- energies and widths come first and are the same in
        both modes, so a vector the caller set (``update_rwa_params_from_sav``)
        is kept; the norm/shift tail is the engine's for the new segment set.
        """
        s = self._engine
        fixed = list(s.params_fixed())
        if fixed == list(self.fixed_params):
            return
        nrm = sum(1 for r in s.parameter_records() if not r["fixed"] and r["kind"] in ("energy", "width"))
        rwa = np.asarray(s.params_values_rwa(), float).copy()
        val = np.asarray(s.params_values(), float).copy()
        old_rwa = np.asarray(self.params_rwa, float).ravel()
        old_val = np.asarray(self.params, float).ravel()
        k = min(nrm, rwa.size, old_rwa.size)
        rwa[:k] = old_rwa[:k]
        k = min(nrm, val.size, old_val.size)
        val[:k] = old_val[:k]
        self.params_rwa = rwa
        self.params = val
        self.fixed_params = fixed
        self._parameters = None

    @engine_call
    def extrap_mode(self):
        """Switch to the <segmentsTest> grids. Re-initializes the engine."""
        self._set_mode(True)

    @engine_call
    def data_mode(self):
        """Switch back to the <segmentsData> segments, which is what chi-squared uses. Re-initializes the engine."""
        self._set_mode(False)
