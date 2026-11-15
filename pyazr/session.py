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
import sys
import threading

import numpy as np

try:
    from . import _azure2
except ImportError as err:                            # pragma: no cover
    raise ImportError(
        "the AZURE2 engine module (pyazr._azure2) is not built.  It is a C++ "
        "extension, not pure Python: configure with CMake (USE_API=ON, the "
        "default) and build -- the module lands in pyazr/, or in "
        "PYAZR_OUTPUT_DIR when that is set; then export PYAZR_ENGINE_DIR=<it>."
    ) from err

from .errors import AZURE2Error, close_names, listed
from . import engine_text, session_checks, session_fit, session_options, session_physics, session_thm
from .locking import engine_call
from .parameter_file import apply_parameter_file
from .parameters import Parameter, Parameters
from .system import Pair, PairSet
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


def basis_of(options):
    """The basis of the fit vector of a session with these options (:attr:`Session.options`): 2 Park when
    ``use_park``, else 1 Brune when ``use_brune``, else 0 standard."""
    if options.get("use_park"):
        return 2
    return 1 if options.get("use_brune") else 0


class Session(
        _HybridPotentials,
        session_physics.SessionPhysics,
        session_thm.SessionTHM, session_fit.SessionFit,
):
    """One calculation session: a format-3 project set up by the engine in this interpreter (options: see __init__).

    Parameters: :attr:`parameters` (the table of the fit parameters: what each entry of the fit vector x is, its start
    ``parameters.x0``), :meth:`transform` (the physical values at x), :attr:`pairs`.  Data: :attr:`datasets`,
    :attr:`energies`, :attr:`cross`.  At a fit vector x: :meth:`evaluate` (every number at x: model values,
    chi-squared, the objective's terms) and :meth:`objective` (the fit's objective as one callable: value, gradient,
    rows, Jacobian, log-probability).  Files: :meth:`write_outputs`, :meth:`save_fit`.
    :meth:`describe` (plain data), :meth:`close` (or a ``with`` block).  Errors: :mod:`pyazr.errors`.
    """

    #: Energy-shift penalties are measured from 0, (dE/sigma)^2, not from the
    #: .azr value (2026-10-25).  rmfit reads this to match the engine it runs.
    SHIFT_PENALTY_FROM_ZERO = True

    @session_options.records_given
    def __init__(self, file, *, cwd=None, extrapolations=False, use_brune=True,
                 ignore_externals=True, transform=True,
                 use_long_wavelength=True, use_gsl_coul=False, use_rmc=False,
                 use_park=False, use_coulomb_table=True, coulomb_table=None,
                 replace_grids=True, use_ec_table=True, ec_convection_current=True, parameter_file=None):
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
        extrapolations : make the session for the project's extrapolations (the program's menu choice 3): then
            :meth:`evaluate` gives the model on the extrapolations' grids (``evaluate(x).data[i]`` belongs to
            ``extrapolations[i]``), every chi-squared field ``None``, and the fit vector holds the extrapolations'
            nuisances; :meth:`objective` and the fits refuse.  ``False`` (the default): the data sets.  Chosen when
            the session is made: a session for the other is a second ``Session`` (``**s.options`` with this keyword
            changed rebuilds the same model).  A project without an active extrapolation record is refused.
        use_brune, ignore_externals, transform, use_long_wavelength :
            R-matrix formalism options, matching the CLI's ``--use-brune``,
            ``--ignore-externals``, ``--no-transform`` and ``--no-long-wavelength``
            (the ``use_*`` options are ON by default, ``transform`` too).
            Every option below except ``cwd``, ``extrapolations`` and ``coulomb_table`` is also a ``<config>`` key of
            the project under its name (``use_brune=0``; ``Project.set_config(use_brune=False)``): a keyword given here
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
        parameter_file : a parameter file (``param.sav``, ``param.par``) whose every value, fixed level parameters
            included, the session takes when it is made, as the program takes its external parameter file
            (:mod:`pyazr.parameter_file`); :attr:`parameters` ``.x0`` is then the file's free values.  Not a
            ``<config>`` key.  To read only the free values into a fit vector: ``Session.parameters.read(path)``.

        Raises
        ------
        pyazr.MissingFile
            (a ``FileNotFoundError``) if the file or ``cwd`` does not exist.
        pyazr.FormatError
            if the file is not a format-3 project (the message names the conversion).
        pyazr.AZURE2Error
            (``pyazr.EngineError``) if the model will not initialize, or with ``extrapolations=True`` if the project
            has no active extrapolation record.

        Warns
        -----
        pyazr.EngineWarning
            for an engine notice met while setting up (an external-capture amplitude file that cannot be written).
        """
        self._sess = None
        self._parameter_file = None if parameter_file is None else os.path.abspath(os.fspath(parameter_file))
        self._engine_lock = threading.RLock()  # every call that enters the engine holds it (pyazr.locking)
        self._azr_cache = {}
        require_format3(file)  # pyazr reads format 3 only: the engine's message names the conversion
        self.file = os.path.abspath(file)
        self.cwd = os.path.abspath(cwd if cwd is not None else (
            os.path.dirname(self.file) or "."))
        if not os.path.isdir(self.cwd):
            raise session_checks.missing_directory(self.cwd, self.file)
        if extrapolations:
            session_checks.require_extrapolations(self.file, self.cwd)

        # The keywords this call named win over the project's <config> keys, which win over the defaults (the engine's
        # rule); self.options holds the values in effect, so that a worker process can rebuild this same model
        # (pyazr.bands.sensitivities with nprocs > 1), and self.option_sources where each came from.
        given = {k: v for k, v in locals().items() if k in self._given_options}
        opts = session_options.runtime_options(_azure2, given, not extrapolations, coulomb_table, self.cwd)
        # the engine's text of the set-up goes to the logger pyazr.engine, not to standard output (pyazr.engine_text);
        # a refused set-up quotes it in its error
        said = []
        try:
            with engine_text.held() as said:
                self._sess = _azure2.Session(self.file, opts)
        except AZURE2Error as err:
            raise session_checks.setup_refused(err, said, self.file) from None
        self.options, self.option_sources = session_options.in_use(self._sess, extrapolations, coulomb_table,
                                                                    self._parameter_file)

        self.mode = "extrap" if extrapolations else "data"  # what the session evaluates, fixed when it is made
        with engine_text.held():
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
        free, bounds, owner fields), and ``functions``: the session's diagnostic calls (``channel_functions`` and the
        THM calls) described as rows -- arguments with unit, kind, bounds and default, returned fields with unit,
        the axis, and where each array sits in the answer (:mod:`pyazr.session_functions`).
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
        self._apply_parameter_file()
        self._fixed = [bool(f) for f in s.params_fixed()]
        self._parameters = None
        self._pairs = None
        self._azr_cache = {}
        session_checks.warn_engine_notices(self._sess)

    def _apply_parameter_file(self):
        """Put the session's parameter file (``parameter_file=``) into the engine, after a set-up or a rebuild."""
        if self._parameter_file is not None:
            apply_parameter_file(self._engine, self.file, self._parameter_file, self.basis, self.mode == "extrap")

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
        """The table of the fit parameters (:class:`pyazr.Parameters`): one row per parameter, the engine's record of
        it (name, kind, unit, start, fixed, free_index, bounds and what it belongs to); ``x0`` is a copy of the start
        of the fit vector, the only vector of the session's calls.  Built once per layout.

        >>> P = s.parameters
        >>> P.of_kind("width").free                 # the free widths, in the order of x
        >>> P.of_level("3/2- 3.5032")               # one level's energy and widths
        >>> print(P.table())                        # every row
        """
        if self._parameters is None:
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

    def _build_pairs(self):
        return PairSet(Pair(r) for r in self._engine.pair_records())

    @property
    @engine_call
    def level_scheme(self):
        """A structured, printable :class:`~pyazr.system.LevelScheme`.

        Groups the model the way it reads physically -- particle pairs, then
        J-groups, then levels and their channels (L, S, radiation type, partial
        width, fixed flag, Wigner limit).  ``print(azr.level_scheme)`` gives a
        human-readable overview.  Read-only; to add/remove levels and write the
        result to a file see :meth:`pyazr.Project.add_level` and the other level edits of :class:`pyazr.Project`.
        """
        from .system import LevelScheme
        return LevelScheme.from_azr(self)

    @property
    @engine_call
    def datasets(self):
        """Per-segment dataset provenance (a :class:`~pyazr.described.DescribedDatasets`, the engine's
        description): for each data set, the data file it came from, the reaction channel (entrance/exit pairs), energy
        / angle range, observable type, and normalization systematic error.
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
        type.  These describe the segments a session made with ``extrapolations=True`` evaluates, in the same
        order, so ``evaluate(x).data[i]`` of such a session belongs to ``extrapolations[i]``.
        A copy of the session's parse of the file (see :meth:`_parsed`).
        """
        return copy.deepcopy(self._parsed("extrapolations"))

    def _build_parameters(self):
        """The table of the engine's parameter records (built once per layout; :meth:`configure` and a mode switch
        drop it)."""
        rows = [Parameter(r) for r in self._engine.parameter_records()]
        return Parameters(rows, x0=self._engine.params_values_rwa(), basis=self.basis)

    # -- the physical values ---------------------------------------------------

    @engine_call
    def transform(self, x):
        """The physical value of every parameter at the fit vector ``x``: a numpy array, one entry per row of
        :attr:`parameters` (fixed ones included), in the project's units: a level energy in MeV, a partial width in
        eV (an ANC in fm^-1/2 for a closed channel, a channel entered as an amplitude as Brune's amplitude), every other
        kind as its fit value (``parameters.units``).

        The engine's transformation (``CompoundSystem::TransformOut``), as the program writes ``<levels>`` from a fit: a fixed
        width holds its reduced-width amplitude, so its physical value follows the level's other widths through the
        Brune denominator ``1 + sum_c gamma_c^2 dS_c/dE`` (under Park the amplitude is the observed width).  The
        inverse is not offered by the engine: a fit vector is made from ``parameters.x0`` (``Parameters.put``) or a
        parameter file (``Parameters.read``).

        >>> phys = s.transform(x)
        >>> dict(zip(s.parameters.names, phys))
        """
        x = session_checks.fit_vector(x, len(self.parameters.x0), "transform")
        full = self._all_rwa(x)
        rmatrix = np.asarray(self._engine.transform_all_rwa(full, include_fixed=True), float)
        out = full.copy()
        out[:rmatrix.size] = rmatrix
        return out

    def _all_rwa(self, x):
        """The full RWA vector with the free vector ``x`` in its free slots."""
        allx = np.asarray(self._engine.params_all_rwa(), float).copy()
        free = [i for i in range(allx.size) if not self._fixed[i]]
        x = np.asarray(x, float).ravel()
        for k, i in enumerate(free[:x.size]):
            allx[i] = x[k]
        return allx

    # -- output files and norms -----------------------------------------------

    @engine_call
    def save_fit(self, path, x=None, param_sav=True, verify=True, norms="fitted", close_session=False, *,
                 relative_to=None):
        """Snapshot the fit vector ``x`` (default ``parameters.x0``) as a project ``path`` that reopens as the fit.

        The level parameters are written as the ``<levels>`` rows state them (:meth:`transform`: energies in MeV,
        partial widths in eV, ANCs), fixed ones included (a fixed width holds its amplitude, so its physical value
        moves with the level's other widths); ``norms="fitted"`` (default) writes every free nuisance into the record
        it comes from, ``"nominal"`` keeps the records' values and leaves the fitted ones to the ``.sav``.  With
        ``param_sav`` (default) ``<name>.sav`` is written beside it (:meth:`Parameters.write`).  With ``verify``
        (default) the file is reopened and must give back the physical values, the particle amplitudes level by level
        and (a session of the data, ``"fitted"``) the data chi-squared and the objective; otherwise both files are
        removed and a ``RuntimeError`` is raised.  ``close_session`` closes this session before the check opens the
        copy.  Relative paths (data files, covariance and THM tables) are written relative to the snapshot's directory,
        so it reopens there; ``relative_to=<dir>`` writes them relative to that directory instead (a text to be opened
        with that working directory).  The rules in full: :func:`pyazr.snapshot.save_fit`.  Returns ``(azr_path,
        sav_path_or_None)``.
        """
        from .snapshot import save_fit
        return save_fit(self, path, x, param_sav=param_sav, verify=verify, norms=norms, close_session=close_session,
                        relative_to=relative_to)

    @engine_call
    def write_outputs(self, x=None):
        """Write the files the program writes at the end of a calculation, for the fit vector ``x``.

        ``AZUREOut_*``, ``chiSquared.out`` and the rest land in the ``.azr``'s output directory, created when it does
        not exist.  An evaluation at ``x`` (default ``parameters.x0``) is run first, so the files describe that
        vector rather than whatever was evaluated last.

        Returns the output directory the engine wrote to (:attr:`output_dir`: the ``<config>`` one, relative to
        :attr:`cwd`), not a fixed ``cwd/output``.
        """
        x = np.asarray(self.parameters.x0 if x is None else x, float).ravel()
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
            # a session of the extrapolations, and the forward pass is all AZUREOut_*
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

        :attr:`datasets` lists every data set record, inactive ones
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
        the key (or in a session of the extrapolations).  ``sum(r**2)`` is ``evaluate(x).terms["prior covariance"]``,
        and ``objective().residuals(x)`` holds these rows (Jacobian L^-1 on the prior's columns)."""
        x = np.asarray(self.parameters.x0 if params is None else params, float).ravel()
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
        out = {}
        for pair in self.pairs:
            out.setdefault(pair.nucleus, PairSet()).append(pair)
        return out

    @property
    def basis(self):
        """Which amplitudes the fit vector holds: 0 standard (constant boundary
        condition), 1 Brune, 2 Park -- the ``#parametrization`` tag AZURE2 writes
        first in ``param.sav``."""
        return basis_of(self.options)

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
        """Every number of the session at the fit vector ``x`` (default ``parameters.x0``) from one engine evaluation:
        a :class:`pyazr.evaluate.Evaluation` (the model of every point, ``data[i].model``, ``.energies``, ``.extra``;
        in a session of the data also the objective, its terms, the data chi-squared and each data set's).  In a session
        of the extrapolations the data sets are the extrapolations and every chi-squared field is ``None``."""
        from .evaluate import evaluate
        return evaluate(self, x)

    def model_curve(self, name, x=None, points=200, *, alone=False, grid="refined", tolerance=3e-3):
        """The model of the data set ``name`` on a fine grid over its range at the fit vector ``x`` (default
        ``parameters.x0``): a :class:`pyazr.model_curve.ModelCurve` (c.m. and lab energies, c.m. angles and the model
        in the unit of ``evaluate(x).data[i].model``; one branch per fixed angle of an excitation function, or per
        fixed energy of an angular distribution; the data set's own points among the grid values).  The engine's
        evaluation of extrapolation records made from the data set's record, in sessions of the extrapolations opened
        and closed by the call (:mod:`pyazr.model_curve`).

        ``grid="refined"`` (the default): ``points`` values per branch seeded with the level energies in range and
        refined where the curve bends to the relative ``tolerance`` (the ENDF writer's rule,
        :func:`pyazr.model_curve.refine_plan`: a narrow level is drawn through; one session per sweep, about ten);
        ``grid="even"``: ``points`` values evenly spaced, one session.  By default the curves of every active data
        set are made together, so that a folded data set's curve is folded on the mesh the data session plans from
        every data set and agrees with ``evaluate`` at the points (4e-9 on T2's folded data sets without a shift);
        ``alone=True`` makes this one alone (cheaper; folded on a mesh of its own points: 2.4e-3 on T2's thin target).
        A data set that is not one function of one variable, or whose columns give each point its own window or
        order, raises :class:`pyazr.ParameterError`."""
        from .model_curve import model_curve
        return model_curve(self, name, x, points, alone=alone, grid=grid, tolerance=tolerance)

    # -- the external region, and the caches that make it affordable ----------

    @engine_call
    def recalculate_external_capture(self):
        """Force a from-scratch recomputation of every external-capture
        integral.

        AZURE2 caches these in ``output/intEC.dat`` (a session of the data) /
        ``output/intEC.extrap`` (of the extrapolations) and reads them back on the
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
