"""High-level driver for one in-process AZURE2 engine instance.

An :class:`azure2` object owns one ``_azure2.Session`` -- the AZUREAPI C++
object bound into Python via pybind11 -- loaded from a single ``.azr`` file.
There are no subprocesses, no sockets and no instance pools to manage: the
engine runs inside the interpreter, so a session is exactly one Python object.
Opening several ``azure2()`` objects in one interpreter is fine, each is
independent; drive them from one thread, as the engine is not reentrant.
"""

import os
import re
import shutil
import tempfile
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
# Imported inside the module rather than at package level: azrfile does not
# depend on the compiled engine, and this keeps that direction of the
# dependency one-way.
from .azrfile import AzrModel


@contextmanager
def _in_dir(path):
    """Run the block with the process cwd at ``path``, then put it back."""
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


_THM_BACKGROUND_TERMS = {"none": 0, "const": 1, "linear": 2, "quadratic": 3}


def _config_output_dir(path):
    """The output directory of a ``.azr``, parsed as Config::ReadConfigFile
    does: the line after the first one of ``<config>``, cut at its last
    ``#`` and trimmed ("" if there is none)."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    try:
        i = next(k for k, line in enumerate(lines) if line == "<config>")
    except StopIteration:
        return ""
    if i + 2 >= len(lines):
        return ""
    line = lines[i + 2]
    if "#" in line:
        line = line[:line.rindex("#")]
    return line.strip(" \n\t\r")

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
                 use_long_wavelength=True, use_gsl_coul=False, use_rmc=False,
                 use_park=False):
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
        use_park : use Park's level-dependent boundary conditions (Phys. Rev.
            C 104, 064612), the CLI's ``--use-park``.  The rwa vector then
            holds the *observed* reduced width amplitudes,
            ``Gamma = 2 P gamma**2`` with no shift-derivative factor.  Switches
            ``use_brune`` on.  Analytic derivatives are available; the engine
            adds a penalty to chi-squared for any level whose overlap
            ``J = 1 - sum gamma**2 dS/dE`` is not positive (:meth:`park_norms`).

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
            use_gsl_coul=bool(use_gsl_coul), use_rmc=bool(use_rmc),
            use_park=bool(use_park))

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

    @property
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

    @staticmethod
    def read_sav(path):
        """``{name: value}`` from a ``param.par`` / ``param.sav`` / ``param.fit``.

        Blank lines, ``#`` comments (the ``#parametrization`` tag) and a bare
        ``parametrization`` row (files written 2026-10-06, before the tag became
        a comment) are skipped.  Read parameter files by NAME, never by row
        position: AZURE2 itself matches names (``AZUREParams::
        ReadUserParameters``), new parameters are appended as the code grows
        (the ``segment_<key>_energy_shift_sqrt`` rows of 2026-10-05), and a
        positional read of a file from another layout is silently off by one
        or more rows (an 8Be+alpha rate came out exactly 2x that way).
        """
        vals = {}
        with open(path) as fh:
            for line in fh:
                t = line.split()
                if len(t) < 2 or t[0].startswith("#") or t[0] == "parametrization":
                    continue
                try:
                    vals[t[0]] = float(t[1])
                except ValueError:
                    continue
        return vals

    def full_rwa_from_sav(self, path):
        """The model's full RWA vector (fixed parameters included, ``p.index``
        order) with the values a parameter file names applied.

        A name the file lacks keeps the session's current value, exactly as
        AZURE2 does on read; a file with no name in common with the model is
        refused rather than applied to nothing.
        """
        vals = self.read_sav(path)
        full = np.asarray(self.sess.params_all_rwa(), float).copy()
        hit = 0
        for p in self.parameters:
            if p.name in vals and 0 <= p.index < full.size:
                full[p.index] = vals[p.name]
                hit += 1
        if hit == 0:
            raise ValueError(f"{path}: no parameter name matches this model "
                             f"(is it a file from another project?)")
        return full

    def update_rwa_params_from_sav(self):
        """Reload params_rwa from the run's param.sav (matched by name)."""
        with _in_dir(self.cwd):
            full = self.full_rwa_from_sav(os.path.join(self.output_dir, 'param.sav'))
        self.params_rwa = [float(v) for v, f in zip(full, self.fixed_params) if not f]

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
        with _in_dir(self.cwd):
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
        just :meth:`~pyazr.azrfile.AzrModel.apply_fit`:

        - **Normalizations and energy shifts are not in it.** They live in the
          ``<segmentsData>`` lines, whose norm (and shift) field AZURE2 reads
          as the *start value* of a fit or calculation and, unless the file
          has an explicit ``prior_centre`` row for it in
          ``<parameterSettings>``, also as the *centre of the prior* (the
          penalty ``((n - n0)/(n0 sigma%/100))^2``, likewise for the shift).
          ``norms`` chooses:

          ``"fitted"`` (default)
              the free norms and free shifts of the fit are written into
              their segment lines, and the centre of every prior that thereby
              leaves its field is written as an explicit ``prior_centre`` row
              (the centre the fit used: the loaded field, or the loaded
              row).  The snapshot run on its own -- CLI calculate, GUI, a
              fresh pyazr session -- gives the fitted data chi-squared *and*
              the fit's priors (both checked by ``verify``), and a fit
              started from it is pulled to the same centres, not to the
              fitted values.  Before October 2026 this wrote the norms to the
              ``.sav`` only (a ``_fit.azr`` run alone gave another
              chi-squared: 18O 2010, 7483 against 1872), and in between
              without the centre rows (a refit was centred on the fitted
              norms).  AZURE2 versions before the rows existed skip them and
              centre the prior on the field.
          ``"nominal"``
              the segment lines keep the loaded file's values; the fitted
              norms and shifts are in the ``.sav`` only (needs
              ``param_sav``).  Run with that ``.sav`` as the external
              parameter file (CLI: "use previous parameters"; GUI: Run tab,
              Parameters File) to reproduce the fit.

          ``param_sav`` (the default) writes the companion ``<name>.sav``
          beside the file, carrying every parameter, free and fixed.
        - **A written file is not a fit until it reads back as one.** With
          ``verify`` (the default) the result is reopened and its own
          parameters transformed back; if any R-matrix value disagrees the
          files are removed and this raises, rather than leaving a snapshot
          that is quietly a mixture.  With ``norms="fitted"`` (data mode) the
          reopened file's data chi-squared and its prior terms
          (:meth:`penalties`) must also equal the fit's (rel. 1e-6).

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
            raise ValueError(f"norms must be 'fitted' or 'nominal', not {norms!r}.")
        if norms == "nominal" and not param_sav:
            raise ValueError("norms='nominal' keeps the fitted norms in the .sav "
                             "only; it needs param_sav=True.")
        x = np.asarray(self.params_rwa if x is None else x, float).ravel()
        # Every parameter's physical value, fixed ones included (see above).
        want = np.asarray(self.transform_all_rwa(self._all_rwa(x), include_fixed=True), float)

        path = os.path.abspath(path)
        model = AzrModel.from_file(self.file)
        model.apply_fit(self.parameters, want, physical=True, pairs=self.pairs,
                        include_fixed=True)
        self._save_cbkg(model, np.asarray(self._all_rwa(x), float))
        if norms == "fitted":
            values, centres = {}, {}
            datasets = self.datasets
            for p in self.parameters:
                if p.kind not in ("norm", "shift") or p.fixed or p.free_index is None:
                    continue
                n, sh = values.get(p.segment_key, (None, None))
                v = float(x[p.free_index])
                values[p.segment_key] = (v, sh) if p.kind == "norm" else (n, v)
                # Keep the prior where the fit had it: an explicit centre for
                # every prior whose field now holds another value.
                d = datasets.by_key(p.segment_key)
                cn, cs = centres.get(p.segment_key, (None, None))
                if p.kind == "norm" and d.norm_error and v != d.nominal_norm:
                    cn = d.nominal_norm
                elif (p.kind == "shift" and d.vary_shift and d.energy_shift_error
                      and v != d.nominal_shift):
                    cs = d.nominal_shift
                if cn is not None or cs is not None:
                    centres[p.segment_key] = (cn, cs)
            if values:
                model.set_segment_values(values)
            if centres:
                model.set_prior_centres(centres)
        model.write(path)

        sav = None
        if param_sav:
            sav = os.path.splitext(path)[0] + ".sav"
            allrwa = list(np.asarray(self.sess.params_all_rwa(), float))
            free = iter(range(len(x)))
            with open(sav, "w") as fh:
                # AZURE2 tags every parameter file with the amplitudes it holds
                # (0 standard, 1 Brune, 2 Park) and converts on read when a run
                # in the other alternative basis loads it.  The '#' keeps
                # numpy.loadtxt-based readers, which index rows by position,
                # unaffected.
                fh.write(f"{'#parametrization':>28s} {float(self.basis): .7e} {0.0: .7e}\n")
                for i, p in enumerate(self.parameters):
                    v = x[next(free)] if not self.fixed_params[i] else allrwa[i]
                    fh.write(f"{p.name:>28s} {float(v): .7e} {0.0: .7e}\n")

        # What the check compares with, from this session before it may close.
        chi_want = pen_want = None
        if verify and norms == "fitted" and self.mode == "data":
            chi_want = self.calculate_chi2_rwa(x)[0]
            pen_want = self.penalties(x)
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
                copy = AzrModel.from_file(path)
                copy.set_output_dir(tmpdir)
                fd, probe = tempfile.mkstemp(suffix=".azr", dir=self.cwd)
                os.close(fd)
                copy.write(probe)
                with azure2(probe, cwd=self.cwd, **check_options) as check:
                    got = np.asarray(check.transform_all_rwa(
                        check._all_rwa(check.params_rwa), include_fixed=True), float)
                    chi_got = None
                    if chi_want is not None:
                        chi_got = check.calculate_chi2_rwa(check.params_rwa)[0]
                        pen_got = check.penalties(check.params_rwa)
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
            if got.shape != want.shape or not np.allclose(got, want, rtol=1e-4,
                                                          atol=0.0):
                n = (0 if got.shape != want.shape
                     else int(np.sum(~np.isclose(got, want, rtol=1e-4))))
                self._discard(path, sav)
                raise RuntimeError(
                    "the snapshot does not read back as the fit it was made "
                    f"from ({n or 'a different number of'} parameter(s) "
                    "disagree); it has been removed. This means the .azr and "
                    "the parameter set do not describe the same model.")
            if chi_got is not None:
                if not np.isclose(chi_got, chi_want, rtol=1e-6, atol=1e-9):
                    self._discard(path, sav)
                    raise RuntimeError(
                        f"the snapshot run on its own gives chi2 = {chi_got:.10g}, "
                        f"the fit {chi_want:.10g}; it has been removed.")
                got_p = sum(float(np.sum(v)) for v in pen_got.values())
                want_p = sum(float(np.sum(v)) for v in pen_want.values())
                if not np.isclose(got_p, want_p, rtol=1e-6, atol=1e-9):
                    self._discard(path, sav)
                    raise RuntimeError(
                        f"the snapshot's priors give {got_p:.10g}, the fit's "
                        f"{want_p:.10g}; it has been removed.")
        return path, sav

    _CBKG_NAME = re.compile(r"([0-9/]+[+-])_(\d+)_([^_]+)_(re|im)([01])")

    def _save_cbkg(self, model, values):
        """Write the THM coherent-background parameters (``cbkg_*``) back into
        the ``cbackground=`` values of their experiments, one explicit term
        per (entrance, exit) combination with its values (``f`` for fixed),
        so the snapshot reads back as the fit.  ``values`` is the full
        parameter vector in the session's order (the cbkg entries are not
        transformed)."""
        cb = [(i, p) for i, p in enumerate(self.parameters) if p.kind == "cbkg"]
        if not cb:
            return
        for name in model.thm_experiments():
            prefix = f"cbkg_{name}_"
            combos = {}
            for i, p in cb:
                if not p.name.startswith(prefix):
                    continue
                m = self._CBKG_NAME.fullmatch(p.name[len(prefix):])
                if not m:
                    continue
                key = (m.group(1), m.group(2), m.group(3))
                slot = (0 if m.group(4) == "re" else 1) + 2 * int(m.group(5))
                combos.setdefault(key, {})[slot] = (float(values[i]),
                                                    bool(self.fixed_params[i]))
            if not combos:
                continue
            terms = []
            for (jpi, exit_key, chans), v in combos.items():
                n = 4 if 2 in v or 3 in v else 2
                vals = ",".join(repr(v.get(k, (0.0, False))[0])
                                + ("f" if v.get(k, (0.0, False))[1] else "")
                                for k in range(n))
                terms.append(f"{jpi}:{exit_key}:{chans}"
                             + (":linear" if n == 4 else "") + "=" + vals)
            model.set_thm_cbackground(name, terms)

    @staticmethod
    def _discard(*paths):
        """Remove a snapshot that failed verification, so nothing bad survives."""
        for f in paths:
            if f and os.path.exists(f):
                os.remove(f)

    # -- transforms -----------------------------------------------------------

    def transform_rwa(self, params):
        """Map a reduced-width-amplitude vector to physical parameters: level energies in MeV, partial widths in eV, ANCs in fm^-1/2."""
        return self.sess.transform_rwa(params)

    def transform_all_rwa(self, params, include_fixed=False):
        """As transform_rwa, taking every parameter rather than the free ones.

        Returns the physical values of the free parameters, or with
        ``include_fixed`` of every parameter.  A *fixed* width holds its
        reduced-width amplitude (as in the CLI's fits), so its physical value
        is not constant: it follows the level's other widths through the
        Brune denominator ``1 + sum_c gamma_c^2 dS_c/dE``.
        """
        return self.sess.transform_all_rwa(params, include_fixed=bool(include_fixed))

    def _all_rwa(self, x):
        """The full RWA vector with the free vector ``x`` in its free slots."""
        allx = np.asarray(self.sess.params_all_rwa(), float).copy()
        free = [i for i in range(allx.size) if not self.fixed_params[i]]
        x = np.asarray(x, float).ravel()
        for k, i in enumerate(free[:x.size]):
            allx[i] = x[k]
        return allx

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
        """Total data chi-squared at a free RWA vector, as a one-element list.

        Excludes the normalization and energy-shift penalties AZURE2's own
        objective adds.  A THM segment with a free norm is treated as the CLI
        treats it: its arbitrary scale is profiled out (n* = S_mm/S_md, chi2 =
        S_dd - S_md^2/S_mm, no penalty) and is not an entry of ``params`` --
        see :meth:`segment_norms` for the n* it lands on.
        """
        return [float(self.sess.calculate_chi2_rwa(params))]

    def calculate_chi2(self, params):
        """As calculate_chi2_rwa, taking the physical parameter vector instead."""
        return [float(self.sess.calculate_chi2_physical(params))]

    def segment_norms(self, params=None):
        """The normalization each segment carries at ``params``.

        One entry per segment, in ``segment_chi2`` order.  For most segments
        this is the free-vector value (or the fixed one from the file); for a
        THM segment with a free norm it is the profiled optimum ``n*`` -- which
        is not a fit parameter, so this is the only way to see it.  The
        segments of a THM experiment (``<thm>`` ``experiment[...]``) all carry
        its shared ``n*`` (see :meth:`thm_background`).  Runs a chi-squared
        evaluation.
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        self.sess.calculate_chi2_rwa(x)
        return np.asarray(self.sess.current_norms(), float)

    def thm_experiments(self, params=None):
        """The THM experiments (``<thm>`` ``experiment[<name>]`` lines) as a
        chi-squared evaluation at ``params`` profiles them.

        Returns ``{name: report}``; see :meth:`thm_background` for a report.
        Empty without experiments (or in extrapolation mode).  Runs a
        chi-squared evaluation.
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        if self.mode == "data":
            self.sess.calculate_chi2_rwa(x)
        out = {}
        for r in self.sess.thm_experiments():
            q = 1 + _THM_BACKGROUND_TERMS[r["background"]]
            cov = np.asarray(r["covariance"], float).reshape(4, 4)[:q, :q]
            value = np.asarray(r["value"], float)
            out[r["name"]] = {
                "segments": list(r["segments"]),
                "background": r["background"],
                "points": int(r["points"]),
                "chi2": float(r["chi2"]),
                "status": r["status"],
                "norm": float(value[0]),
                "b": value[1:q].copy(),
                "cov": cov,
                "sigma": np.sqrt(np.maximum(np.diag(cov), 0.0)),
            }
            if "cbkg_names" in r:
                # Coherent background (cbackground=): its fit parameters,
                # name -> value, as they are in the parameter vector.
                out[r["name"]]["cbkg"] = dict(zip(r["cbkg_names"],
                                                  map(float, r["cbkg_values"])))
        return out

    def thm_background(self, experiment, params=None):
        """Shared norm and background of one THM experiment at ``params``.

        The segments of ``experiment[<name>]`` share one profiled norm ``n*``
        (it multiplies the data, as ``segment_norms``) and a background
        ``b(E) = b0 + b1 E + b2 E^2`` added to the folded HOES model (model
        units, E the c.m. energy of the THM entrance pair, MeV); both are
        eliminated by closed-form weighted linear least squares at every
        evaluation.  Returns a dict:

        ``norm``        n*
        ``b``           array of the background coefficients (length 0-3)
        ``cov``         covariance of (norm, b0, ...) from the profile, at fixed
                        R-matrix parameters (not scaled by chi2/nu)
        ``sigma``       sqrt(diag(cov))
        ``chi2``        the experiment's chi-squared, ``points``, ``segments``
        ``background``  none | const | linear | quadratic
        ``status``      "profiled", or what a degenerate profile fell back to

        The model at a point, as the output files show it next to the data
        scaled by n*, is ``calculate_rwa`` plus ``b(E)``.  Runs a chi-squared
        evaluation.  KeyError for an unknown experiment.
        """
        reports = self.thm_experiments(params)
        if experiment not in reports:
            raise KeyError(f"no THM experiment {experiment!r} "
                           f"(have: {', '.join(reports) or 'none'})")
        return reports[experiment]

    def thm_lineshape(self, experiment, energies, params=None):
        """The Coulomb line shape N_C of a THM experiment (``lineshape=on``).

        In A + a(x+s) -> s + F*(x+A) -> s + b + B the spectator s leaves in
        the Coulomb field of F*, and of b and B once F* has decayed; each
        level's exit amplitude is multiplied by
        ``N_C = exp(pi zeta/2) (E_l - E - i G_l/2)^(-i zeta)``, so for an
        isolated level ``|N_C|^2 = exp[2 zeta arctan(2 (E_l - E)/G_l)]``
        (Mukhamedzhanov, Kadyrov & Pang, EPJA 56 (2020) 233, eqs. 55-62), with
        ``zeta = eta_sB - eta_0`` (their case 2: B much heavier than s and b,
        the s-b term eta_sb neglected), ``eta_0 = Z_s Z_F alpha mu_sF/k_sF``,
        ``eta_sB = Z_s Z_B alpha mu_sB/k_sF`` and ``E_sF = E_aA - B - E``.
        See docs/source/theory/thm_implementation.rst, "Coulomb line shape".

        ``energies``: c.m. energies E of x + A (MeV).  The level poles
        (observed energy and total width, Brune) are those of ``params``
        (default: the current parameters; a chi-squared evaluation is run to
        load them in data mode).  Returns a dict:

        ``E``, ``E_sF``, ``eta_0``   arrays on the grid
        ``E_aA``, ``B``, ``Zs``, ``ZF``, ``spectator``
        ``exits``   one dict per exit pair of the experiment's segments:
                    ``pair`` (key), ``Zb``/``ZB``, ``mb``/``mB`` (light and
                    heavy nucleus), ``zeta`` and ``eta_sb`` (the neglected
                    s-b parameter averaged over the b direction, a validity
                    check: it should be << 1) on the grid, and ``levels``: per
                    level of a J group coupling entrance and exit, ``jgroup``,
                    ``level`` (1-based, as parameters.out), ``J``, ``pi``,
                    ``E_level`` (c.m. of x + A), ``Gamma`` (MeV) and ``NC2``
                    (|N_C|^2 on the grid).

        Raises the engine error if the experiment is unknown or has no line
        shape.
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        if self.mode == "data":
            self.sess.calculate_chi2_rwa(x)
        r = self.sess.thm_lineshape(str(experiment), np.asarray(energies, float).ravel())
        return r

    def thm_vertex(self, experiment, energies, params=None, strict=False):
        """The THM entrance vertex of an experiment, averaged over its
        spectator-momentum window (``ps=`` on the experiment line).

        The off-shell x-A momentum is ``p_xA^2/2mu_xA = E + B + p_s^2/2mu_sx``
        (Mukhamedzhanov et al., PRC 96 (2017) 024623, eq. 31), so the vertex
        ``M_l = (B_c - 1) j_l(rho) - rho j_l'(rho) (+ C_l)``, ``rho = p_xA a``,
        changes with the spectator momentum p_s inside the accepted window.
        The engine averages the cross section (not the amplitude) over the
        window with the event weight |phi(p_s)|^2 p_s^2; this reports the
        vertex part of it.  See docs/source/theory/thm_implementation.rst,
        "Spectator-momentum window".

        ``energies``: c.m. energies E of x + A (MeV).  The per-level vertex
        boundaries are those of ``params`` (default: the current parameters;
        a chi-squared evaluation is run to load them in data mode).  Returns a
        dict:

        ``window``    description (``"delta"`` without a window)
        ``p_s``, ``weights``, ``T_s``   the Gauss-Legendre nodes (MeV/c), their
                      normalized weights and T_s = p_s^2/2mu_sx (MeV); a delta
                      experiment has one node, p_s = 0 and T_s the
                      ``spectatorEnergy`` of its pair
        ``mu_sx``, ``B``, ``radius``, ``pair``   (MeV, MeV, fm, entrance key)
        ``E``         the grid; ``rho``: list over E of the rho at each node
        ``channels``  one dict per entrance channel of every J group:
                      ``jgroup``, ``channel`` (1-based), ``J``, ``pi``, ``l``,
                      ``s`` and ``levels``: per level ``level`` (1-based),
                      ``boundary`` (the vertex boundary B_c for that level;
                      NaN for ``vertex=onshell``, where it is S(E) + iP(E)),
                      ``M2`` (<|M_l|^2> over the window) and ``M2_qf``
                      (|M_l|^2 at p_s = 0) on the grid.
        ``model``     ``"pw"``, or ``"dw"`` for an experiment with
                      ``vertexModel=dw``: then the vertex is the surface term of
                      the prior-form DWBA (docs, "Distorted-wave entrance
                      vertex"), ``M2`` is its average over the nodes the engine
                      uses (the reachable part of the ``ps`` window at each E, or
                      the one ``spectatorAngle`` node), ``M2_qf`` its value at
                      the ``spectatorAngle`` direction and ``M2_pw`` the
                      plane-wave |M_l(p)|^2 there, p = |k_aA - alpha k_sF| (what
                      the DW vertex becomes without distortion); ``rho`` is p a
                      at each node, ``dw_q``/``dw_weights`` the nodes (q in MeV/c)
                      and weights per energy, ``dw_q_delta`` (MeV/c) and
                      ``dw_p_delta`` (fm^-1) those of the ``spectatorAngle``
                      node.  ``p_s``, ``weights``, ``T_s`` stay the nominal
                      window of the line.  With ``spectatorAngles=`` the
                      vertex is averaged over the accepted spectator
                      directions (the model is linear in the Gram matrix, so
                      that is one node with the averaged matrix: ``dw_q`` its
                      weighted mean q), and ``angle_theta_cm`` (deg, c.m. to
                      the beam), ``angle_q`` (MeV/c), ``angle_weights`` give
                      the directions at the grid energy nearest to each E.
        ``reached``   bool array over E: False where the window (or, for
                      ``dw``, the vertex grid) does not reach E.  The engine
                      then uses the nodes and vertex of the nearest data point
                      (nearest grid energy for ``dw``) -- right for a folding
                      sub-point next to the data, but the row at such an E is
                      not the vertex *at* E (several rows can come out
                      identical).

        Energies out of reach are not silent: a ``UserWarning`` names them,
        or with ``strict=True`` a ``ValueError`` is raised instead.

        Raises the engine error if the experiment is unknown.
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        if self.mode == "data":
            self.sess.calculate_chi2_rwa(x)
        r = self.sess.thm_vertex(str(experiment), np.asarray(energies, float).ravel())
        reached = np.asarray(r.get("reached", []), bool)
        if reached.size and not reached.all():
            far = np.asarray(r["E"], float)[~reached]
            shown = ", ".join(f"{e:.6g}" for e in far[:8]) + (" ..." if far.size > 8 else "")
            where = ("outside the DW vertex grid or beyond the reach of its spectator "
                     "directions" if r.get("model") == "dw" else
                     "beyond the reach of the experiment's spectator window")
            msg = (f"thm_vertex({experiment!r}): {far.size} of {reached.size} energies "
                   f"are {where} ({r.get('window', '')}); the vertex reported there "
                   f"is that of the nearest data point / grid energy, not the vertex "
                   f"at E (E = {shown} MeV; see the 'reached' entry)")
            if strict:
                raise ValueError(msg)
            import warnings
            warnings.warn(msg, UserWarning, stacklevel=2)
        return r

    def thm_distortion(self, experiment, energies):
        """The distortion factor R(E) of a THM experiment (``distortion=`` on
        its experiment line).

        Zero-range prior-form DWBA (Mukhamedzhanov & Pang, PRC 99 (2019)
        064618, eqs. 20-24; Mukhamedzhanov, arXiv:2609.04498, eqs. 22-30):
        ``M(E) = <chi(-)_{k_sF} | phi_sx | chi(+)_{k_aA}(beta r)>``, beta =
        m_s/m_a, summed over partial waves; ``M_PW`` is its plane-wave limit,
        the Fourier transform of the s-x bound state at q = |k_sF - beta
        k_aA|, which the PWA data reduction divides by.  The model is
        multiplied by ``R(E) = rho(E)/rho(E_ref)`` before the folding, with
        rho = |M|^2/|M_PW|^2 (``distortionRatio=dwpw``, default) or |M|^2
        (``dw``, the published convention); equivalently the PWA-extracted
        S* divided by R.  See docs/source/theory/thm_implementation.rst,
        "Distortion factor R(E)".  R does not depend on the R-matrix
        parameters.

        ``energies``: c.m. energies E of x + A (MeV).  Returns a dict:

        ``kind``, ``description``   coulomb | optical | table, the settings
        ``E``, ``R_model``          the grid and the w(E) the model uses
                                    (interpolated; the table's for
                                    ``table:<file>``, then nothing else)
        ``R``                       R evaluated directly at each E
        ``M2``, ``M2_PW``           |M|^2 and |M_PW|^2 (fm^6 up to a constant)
        ``E_sF``, ``k_sF``, ``eta_sF``   the spectator-F channel (MeV, fm^-1)
        ``theta_cm``, ``x``, ``q``  spectator c.m. angle to the beam (deg),
                                    cos(k_sF, k_aA), |k_sF - beta k_aA|
        ``lmax``                    the highest partial wave summed

        With ``spectatorAngles=`` (a window of spectator directions) ``M2``
        and ``M2_PW`` are the averages over the accepted directions (weight
        d cos theta_cm x acceptance), R is their ratio (``dwpw``), and
        ``theta_cm``, ``x``, ``q`` are weighted means.
        ``E_ref``, ``E_aA``, ``B``, ``k_aA``, ``eta_aA``, ``kappa``,
        ``eta_b``, ``beta``         scalars of the experiment

        Raises the engine error if the experiment is unknown, has no
        distortion, or an energy leaves the spectator no energy.
        """
        return self.sess.thm_distortion(str(experiment), np.asarray(energies, float).ravel())

    def residuals(self, params=None):
        """Standardized residuals ``(fit_i - data_i*n)/(cmErr_i*n)``, from a
        forward pass alone.

        The residual vector of :meth:`residual_jacobian` (same rows, in segment
        order) without its Jacobian, so it is as cheap as one chi-squared and
        never refuses a model.  ``sum(r**2) == calculate_chi2_rwa``: THM
        segments with a free norm use their profiled ``n*``.
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        return np.asarray(self.sess.calculate_residuals_rwa(x.ravel()), float)

    def write_output_files(self, params=None):
        """Write the run's standard output files, as the CLI does at the end.

        ``AZUREOut_*``, ``chiSquared.out`` and the rest land in the ``.azr``'s
        output directory, which must already exist.  Without this a Python
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
        with _in_dir(self.cwd):
            if self.mode == "data":
                # chiSquared.out reports the per-segment chi-squared *stored on
                # each segment*, and only the chi-squared evaluation sets it --
                # a bare forward pass leaves it at zero and the file comes out
                # full of zeros.  This populates both that and the point cross
                # sections AZUREOut_* is written from.
                self.sess.calculate_chi2_rwa(x)
            else:
                # No data to compare against; chiSquared.out is meaningless in
                # extrapolation mode, and the forward pass is all AZUREOut_*
                # needs.
                self.sess.update_segments_rwa(x)
            if not self.sess.write_output_files():
                raise RuntimeError("AZURE2 has no data object to write from.")
        return self.output_dir

    # -- chi-squared, per segment and per dataset -----------------------------

    def segment_chi2(self, params=None):
        """Chi-squared per data segment, in ``<segmentsData>`` order.

        ``calculate_chi2_rwa`` returns only the total, but the standardized
        residuals carry the split: they come back in segment order, so summing
        their squares between the segment boundaries recovers each one.  The
        sum of this equals ``calculate_chi2_rwa`` exactly (THM segments with
        a free norm at their profiled optimum, as in the CLI's chiSquared.out).
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        r = self.residuals(x)
        edges = np.cumsum([0] + [len(self.energies[i]) for i in range(self.nsegments)])
        return np.array([float(np.sum(r[edges[i]:edges[i + 1]] ** 2))
                         for i in range(self.nsegments)])

    def dataset_chi2(self, params=None):
        """Chi-squared per *dataset*, which is how a fit is usually judged.

        An experiment often owns several segments -- one per angle, or one per
        final state -- so the per-segment split is finer than the question
        being asked.  Returns ``{name: (chi2, points, segments)}`` keyed by the
        data file's name, in first-appearance order.
        """
        seg = self.segment_chi2(params)
        out = {}
        for i, d in enumerate(self.active_datasets):
            name = d.name
            chi2, n, k = out.get(name, (0.0, 0, 0))
            out[name] = (chi2 + float(seg[i]), n + len(self.energies[i]), k + 1)
        return out

    @property
    def active_datasets(self):
        """The data segments the engine evaluates, in the order its results
        come back: ``m.energies[i]``, ``m.cross[i]``, ``segment_chi2(x)[i]``
        all belong to ``active_datasets[i]``.

        :attr:`datasets` lists every ``<segmentsData>`` line, inactive ones
        included, so ``datasets[i]`` is *not* segment ``i`` of a result on any
        model that has one -- the 13C+alpha archive's model has 25 -- and code
        that indexed it that way dropped the normalization penalty of every
        active segment past the active count.
        """
        active = self.datasets.active
        if len(active) != self.nsegments:
            raise RuntimeError(
                f"{len(active)} active <segmentsData> lines in {self.file} but "
                f"the engine reports {self.nsegments} segments.")
        return active

    # -- AZURE2's own objective -----------------------------------------------

    def penalties(self, params=None):
        """The prior terms AZURE2 adds to chi-squared, per segment.

        ``AZURECalc::operator()`` does not minimize the data chi-squared alone.
        For every segment it also adds

            ((norm - nominal) / (nominal/100 * norm_error))^2

        and, for a segment whose energy shift is free,

            ((shift - nominal_shift) / shift_error)^2

        and the same again for a free sqrt(E) shift coefficient (the b of
        E' = E + shift + b*sqrt(E), key ``"shift_sqrt"``).

        Minimize the bare chi-squared instead and the normalizations drift to
        absorb every discrepancy -- a "better" number AZURE2 would never have
        found, worth -480 on the 7Be model.  Roll your own Minuit or
        least-squares fit against :meth:`objective`, not
        :meth:`calculate_chi2_rwa`.

        ``nominal`` is the segment's ``nominal_norm`` (``nominal_shift``): the
        norm (shift) field of its ``<segmentsData>`` line, or its explicit
        ``prior_centre`` row in ``<parameterSettings>`` when it has one (the
        field is then only the start value).  A fixed norm away from an
        explicit centre pays its penalty too, as in the engine.
        Note the denominator uses the *nominal* normalization, and that
        ``norm_error`` is a percentage.  A THM segment with a free norm has no
        norm penalty: its scale is profiled out, not a parameter.  Returns
        ``{"norm": array, "shift": array, "shift_sqrt": array, "park": float}``,
        one entry per segment (``park`` is one number).
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        norm = np.zeros(self.nsegments)
        shift = np.zeros(self.nsegments)
        shift_sqrt = np.zeros(self.nsegments)
        current = {p.segment_key: p for p in self.parameters.norms}
        shifting = {p.segment_key: p for p in self.parameters.shifts}
        sqrting = {p.segment_key: p for p in self.parameters.sqrt_shifts}
        for i, d in enumerate(self.active_datasets):
            p = current.get(d.key)
            value = (float(x[p.free_index])
                     if p is not None and not p.fixed and p.free_index is not None
                     and p.free_index < x.size else d.norm)
            sigma = d.nominal_norm / 100.0 * d.norm_error
            if sigma:
                norm[i] = ((value - d.nominal_norm) / sigma) ** 2
            if d.vary_shift and d.energy_shift_error:
                q = shifting.get(d.key)
                sval = (float(x[q.free_index])
                        if q is not None and not q.fixed and q.free_index is not None
                        and q.free_index < x.size else d.energy_shift)
                shift[i] = ((sval - d.nominal_shift) / d.energy_shift_error) ** 2
            if d.vary_shift_sqrt and d.energy_shift_sqrt_error:
                q = sqrting.get(d.key)
                bval = (float(x[q.free_index])
                        if q is not None and not q.fixed and q.free_index is not None
                        and q.free_index < x.size else d.energy_shift_sqrt)
                shift_sqrt[i] = ((bval - d.energy_shift_sqrt)
                                 / d.energy_shift_sqrt_error) ** 2
        # Park formalism: the J > 0 wall, sum (J/1e-3)^2 over levels with J < 0.
        # The engine adds it to calculate_chi2_rwa itself (it is part of what
        # AZURE2 minimizes), so objective() does not add it again; it is
        # reported here so a fit can see when it is active.  Always 0 in Brune
        # mode.  Note this is one number, not one per segment.
        park = 0.0
        if self.options.get("use_park"):
            J = np.asarray(self.sess.park_norms(x), float)
            park = float(np.sum((J[J < 0] / 1.0e-3) ** 2))
        return {"norm": norm, "shift": shift, "shift_sqrt": shift_sqrt, "park": park}

    @property
    def basis(self):
        """Which amplitudes ``params_rwa`` are: 0 standard (constant boundary
        condition), 1 Brune, 2 Park -- the ``#parametrization`` tag AZURE2 writes
        first in ``param.sav``."""
        if self.options.get("use_park"):
            return 2
        return 1 if self.options.get("use_brune") else 0

    def park_norms(self, params=None):
        """Park formalism: the overlap J = 1 - sum_c gamma_c^2 dS_c/dE of every
        R-matrix level, in ``physical_levels()`` order.

        J is the squared norm of the level's basis state (Park 2021, Eq. 28)
        and the factor between the two parametrizations,
        ``gamma_Park = gamma_Brune * sqrt(J)``.  Every real Brune parameter set
        has 0 < J <= 1; a Park parameter set with J <= 0 has no R-matrix
        counterpart (its widths exceed what the channel radii allow) and the
        engine penalizes it -- see :meth:`penalties`.  All ones in Brune mode.
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        return np.asarray(self.sess.park_norms(x), float)

    def objective(self, params=None):
        """What AZURE2's own fit minimizes: chi-squared plus the penalties.

        This is the number to hand a minimizer.  ``chiSquared.out`` reports the
        two halves separately, as ``Total-Chi-Squared`` and
        ``Total-Norm-Chi-Squared``.
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        pen = self.penalties(x)
        # pen["park"] is already inside calculate_chi2_rwa (engine side).
        return (float(np.sum(self.calculate_chi2_rwa(x)))
                + float(np.sum(pen["norm"])) + float(np.sum(pen["shift"]))
                + float(np.sum(pen["shift_sqrt"])))

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
            shifts are finite-differenced.  THM (HOES) points are outside the
            adjoint: their contribution, with the norm profiled, is
            ``2 J^T r`` over the THM rows of :meth:`residual_jacobian` (model
            Jacobian by central differences, two THM forward passes per free
            R-matrix parameter, like the CLI's own THM gradient).  (For a Gaussian log-likelihood use
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

        THM (HOES) segments: the adjoint does not cover the HOES observable,
        so the model Jacobian ``J_m`` of a THM point is taken by central
        differences (step ``1e-6 (|x| + 1)``) in every free level energy and
        reduced width -- two THM forward passes per column.  With a free norm
        the scale is profiled: ``r_i = (s m_i - d_i)/e_i`` with
        ``s = 1/n* = S_md/S_mm`` a function of the parameters, and its
        dependence is added analytically,
        ``dr_i/dp = (s J_mi + m_i ds/dp)/e_i``,
        ``ds/dp = (sum d_k J_mk/e_k^2 - 2 s sum m_k J_mk/e_k^2)/S_mm``.
        The segments of a THM experiment share the scale and a background,
        profiled together by linear least squares; ``J`` then carries the
        exact derivative through that profile (Golub-Pereyra; see
        docs/source/theory/thm_implementation.rst, "THM experiments").
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

        THM (HOES) points are outside the adjoint; their rows are central
        differences of the HOES model (see :meth:`residual_jacobian`).

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
        shift moves every point energy at every iteration.  Such a key keeps
        its values for the exact same energy only (counted in the totals), so
        a later session in this process finds them; it drops even those if
        they fill up without being asked for again.
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

    def _geometry_rwa(self, params):
        """The free RWA vector a geometry query (energies, angles) runs at.

        Point energies and angles do not depend on the R-matrix parameters,
        only on the energy shifts (and nothing on the norms): so the R-matrix
        block is the session's own amplitudes, and the rest -- norms, shifts,
        THM coherent-background terms, identical in the physical and the RWA
        vector -- is taken from ``params``.  Running the RWA path means no
        observed-to-formal transformation, which is what used to print
        "Denominator less than zero while transforming" when a caller handed
        these methods ``params_rwa``: an amplitude read as a width (or a
        sub-threshold ANC) that no reduced width reproduces.  The CLI never
        transforms an RWA vector, so it never printed it.
        """
        rwa = np.array(self.params_rwa, float).ravel()
        if params is None:
            return rwa
        p = np.asarray(params, float).ravel()
        if p.size != rwa.size:
            raise ValueError(f"expected {rwa.size} free parameters, got {p.size}.")
        n = self.n_rmatrix
        rwa[n:] = p[n:]
        return rwa

    def calculate_excitation_energy(self, params=None):
        """Compound-nucleus excitation energy per segment -- the one axis every entrance pair shares.

        ``params`` may be the physical or the RWA free vector (or None): only
        its energy shifts matter (see :meth:`calculate_energies`).
        """
        s = self.sess
        nsegments = int(s.update_segments_rwa(self._geometry_rwa(params)))
        return [s.calculated_excitation_energies(i) for i in range(nsegments)]

    def calculate_angles(self, params=None):
        """Per-segment angles of the calculated points, one array per segment.

        The companion to :meth:`calculate_energies`: for a differential segment
        the returned angles are AZURE2's own (center-of-mass) values, which
        differ from the lab angles declared in the ``.azr`` file.  ``params``
        as for :meth:`calculate_energies`.
        """
        s = self.sess
        nsegments = int(s.update_segments_rwa(self._geometry_rwa(params)))
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

    def calculate_energies(self, params=None):
        """Centre-of-mass energies of the calculated points, per segment.

        The energies do not depend on the R-matrix parameters, only on the
        segments' energy shifts, so ``params`` may be the physical vector (as
        :meth:`calculate` takes), the RWA vector (:attr:`params_rwa`) or None
        (the session's own): its norm/shift tail is used, its R-matrix block
        is not.  The evaluation runs on the session's amplitudes without the
        observed-to-formal transformation -- before October 2026 it took the
        physical path, and ``params_rwa`` passed here was read as widths and
        ANCs (a sub-threshold ANC channel then printed "Denominator less than
        zero while transforming", which the CLI never does).
        """
        s = self.sess
        nsegments = int(s.update_segments_rwa(self._geometry_rwa(params)))
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
