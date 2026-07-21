"""High-level driver for one or more AZURE2 API instances.

An :class:`azure2` object owns a pool of ``nprocs`` AZURE2 subprocesses, each
paired with a client connection.  Methods that take a ``proc`` argument let a
caller dispatch work to a specific instance (e.g. for parallel chi-squared
evaluations across a worker pool).
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np

from .client import client
from .parameters import Pair, PairSet, Parameter, ParameterSet
from .server import server


class azure2:

    def __init__(self, file, nprocs=1, port=20000, binary=None,
                 verbose=False, auto_port=True, cwd=None):
        """Launch ``nprocs`` AZURE2 instances bound to ``file``.

        Parameters
        ----------
        file : the ``.azr`` configuration file.
        nprocs : number of parallel AZURE2 instances to spawn.
        port : base TCP port used only when ``auto_port=False`` (instance ``i``
            then binds ``port + i``).  Ignored when ``auto_port`` is set.
        binary : path to the AZURE2 executable (auto-detected if ``None``).
        verbose : forward subprocess output to the console.
        auto_port : let the OS assign each instance a unique free port (bind
            port 0; the server reports back the port it got).  This is race-free,
            unlike probing, so concurrent instances can never collide.  The
            actual ports are available as :attr:`ports` after construction.
        cwd : working directory for the subprocesses.  Defaults to the ``.azr``
            file's directory, since a model names its ``output/`` and ``checks/``
            directories relative to the running process.
        """
        self.file = file
        self.nprocs = nprocs
        self.binary = binary
        self.verbose = verbose
        self.cwd = cwd

        # Instance-level lists.  (The original code declared `servers` as a
        # *class* attribute, so every azure2 object shared -- and leaked into
        # -- the same list.)
        self.servers = []
        self.clients = []

        # With auto_port we no longer *probe* for free ports (which raced: the
        # probe released the port before the server bound it).  Instead we ask
        # the OS to assign one via port 0; each server reports the port it
        # actually got, so two instances can never collide.  Explicit ports
        # (auto_port=False) are still honored verbatim for callers who pin them.
        if auto_port:
            self.requested_ports = [0] * nprocs
        else:
            self.requested_ports = [port + i for i in range(nprocs)]
        self.ports = list(self.requested_ports)   # filled in with real ports

        self._closed = False
        try:
            self.spawn()
            self.configure()
        except Exception:
            self.close()
            raise

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # -- lifecycle ------------------------------------------------------------

    def spawn(self):
        """Start the subprocesses, connect clients, and initialize them.

        Connecting and (especially) the ``INITIALIZE`` step are run
        concurrently across all instances.  ``INITIALIZE`` triggers the heavy
        R-matrix setup inside each AZURE2 process; done sequentially the total
        cost grows linearly with ``nprocs``.  Because ``communicate`` blocks on
        socket I/O -- which releases the GIL -- dispatching it from a thread per
        instance lets every AZURE2 process initialize in parallel, so the wall
        time stays close to a single instance's setup time.
        """
        # Launch every subprocess first; Popen is non-blocking, so all the
        # AZURE2 processes start booting concurrently.
        for p in self.requested_ports:
            self.servers.append(
                server(p, self.file, binary=self.binary, verbose=self.verbose,
                       cwd=self.cwd)
            )

        if self.nprocs == 1:
            # Avoid thread-pool overhead in the common single-instance case.
            # Wait for the server to report its actual (OS-assigned) port, then
            # connect to exactly that port -- no probing, no race.
            self.ports = [self.servers[0].wait_until_listening()]
            self.clients = [client(port=self.ports[0])]
            self.clients[0].communicate("INITIALIZE", [0])
            return

        with ThreadPoolExecutor(max_workers=self.nprocs) as pool:
            # Wait for every server's bound port in parallel (each blocks on the
            # heavy startup), keeping order so self.ports/clients stay aligned
            # with `proc`.
            self.ports = list(pool.map(lambda s: s.wait_until_listening(),
                                       self.servers))

            # connect() still polls until the port is accepting, run in parallel
            # to overlap any residual latency.
            self.clients = list(pool.map(lambda p: client(port=p), self.ports))

            # Fan out the heavy INITIALIZE so the processes set up in parallel.
            list(pool.map(lambda c: c.communicate("INITIALIZE", [0]),
                          self.clients))

    def close(self):
        """Disconnect all clients and terminate all subprocesses."""
        if getattr(self, "_closed", False):
            return
        self._closed = True
        for c in getattr(self, "clients", []):
            try:
                c.disconnect()
            except Exception:
                pass
        for s in getattr(self, "servers", []):
            try:
                s.stop()
            except Exception:
                pass
        self.clients = []
        self.servers = []

    def is_alive(self):
        return bool(self.servers) and all(s.is_alive() for s in self.servers)

    def configure(self):
        c = self.clients[0]
        self.nsegments = int(c.communicate("UPDATE_DATA", [0])[0])
        rng = range(self.nsegments)
        self.energies = [c.communicate("GET_DATA_ENERGIES", [i]) for i in rng]
        self.excitation_energies = [c.communicate("GET_DATA_EXCITATION_ENERGY", [i]) for i in rng]
        self.angles = [c.communicate("GET_DATA_ANGLES", [i]) for i in rng]
        self.cross = [c.communicate("GET_DATA_SEGMENTS", [i]) for i in rng]
        self.cross_err = [c.communicate("GET_DATA_SEGMENTS_ERRORS", [i]) for i in rng]
        self.conv = [c.communicate("GET_DATA_CONV", [i]) for i in rng]
        self.sfactor = [self.cross[i] * self.conv[i] for i in rng]
        self.sfactor_err = [self.cross_err[i] * self.conv[i] for i in rng]
        self.params = c.communicate("GET_PARAMS", [0])
        self.params_rwa = c.communicate("GET_PARAMS_RWA", [0])
        self.fixed_params = c.communicate("GET_PARAMS_FIXED", [0])
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

    def refresh_parameters(self, proc=0):
        """Re-fetch and rebuild the cached :attr:`parameters`."""
        self._parameters = self._build_parameters(proc=proc)
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

    def refresh_pairs(self, proc=0):
        """Re-fetch and rebuild the cached :attr:`pairs`."""
        self._pairs = self._build_pairs(proc=proc)
        return self._pairs

    def _build_pairs(self, proc=0):
        flat = self.clients[proc].communicate("GET_PAIRS_INFO", [0])
        nfields = Pair._NFIELDS
        records = np.asarray(flat, dtype=float).reshape(-1, nfields)
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

    def _build_parameters(self, proc=0):
        c = self.clients[proc]
        n = len(self.fixed_params)

        flat = c.communicate("GET_PARAMS_INFO", [0])
        nfields = Parameter._NFIELDS
        records = np.asarray(flat, dtype=float).reshape(-1, nfields)
        if records.shape[0] != n:
            raise RuntimeError(
                f"GET_PARAMS_INFO returned {records.shape[0]} parameters but "
                f"GET_PARAMS_FIXED reported {n}."
            )

        params = ParameterSet()
        free_counter = 0
        for i in range(n):
            raw_name = c.communicate("GET_PARAMS_NAMES", [i])
            name = "".join(chr(int(round(x))) for x in raw_name)
            fixed = bool(round(self.fixed_params[i]))
            free_index = None if fixed else free_counter
            if not fixed:
                free_counter += 1
            params.append(
                Parameter.from_record(i, name, records[i], free_index)
            )
        return params

    # -- index queries --------------------------------------------------------

    def norm_indices(self, proc=0):
        return self.clients[proc].communicate("GET_NORM_INDICES", [0])

    def shift_indices(self, proc=0):
        return self.clients[proc].communicate("GET_SHIFT_INDICES", [0])

    def energy_indices(self, proc=0):
        """Packed indices of the free R-matrix level-energy parameters.

        Returns the positions of the level-energy parameters within
        ``params_rwa`` (the non-fixed parameter vector), the same convention
        as :meth:`norm_indices` and :meth:`shift_indices`.  Unlike those --
        which the C++ side derives by substring-matching parameter names --
        this uses the structured ``type`` code from ``GET_PARAMS_INFO``
        (``kind == "energy"``), so it never depends on how energies are named.
        """
        return [p.free_index for p in self.parameters
                if p.kind == "energy" and not p.fixed]

    # -- param.sav helpers ----------------------------------------------------

    def update_rwa_params_from_sav(self):
        all_rwa_params = np.loadtxt('output/param.sav', usecols=(1,))
        self.params_rwa = []
        for i in range(len(all_rwa_params)):
            if self.fixed_params[i]:
                continue
            else:
                self.params_rwa.append(all_rwa_params[i])

    def update_sav_from_rwa_params(self, best):
        params_full = []
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

    def transform_rwa(self, params, proc=0):
        return self.clients[proc].communicate("TRANSFORM_RWA", params)

    def transform_physical(self, params, proc=0):
        return self.clients[proc].communicate("TRANSFORM_PHYSICAL", params)

    def transform_all_rwa(self, params, proc=0):
        return self.clients[proc].communicate("TRANSFORM_ALL_RWA", params)

    # -- chi-squared ----------------------------------------------------------

    def calculate_chi2_rwa(self, params, proc=0):
        return self.clients[proc].communicate("CALCULATE_CHI2_RWA", params).tolist()

    def calculate_chi2(self, params, proc=0):
        return self.clients[proc].communicate("CALCULATE_CHI2", params).tolist()

    def chi2_and_grad(self, params, proc=0):
        """Value and analytic gradient of the (data) chi-squared.

        Parameters
        ----------
        params : the non-fixed RWA parameters (energies, reduced widths,
            normalizations), in the same order as ``self.params_rwa``.
        proc : which AZURE2 instance to dispatch to.

        Returns
        -------
        (float, numpy.ndarray)
            ``(chi2, grad)`` with one gradient entry per input parameter.
            Energies / reduced widths / normalizations are analytic; energy
            shifts are finite-differenced.  (For a Gaussian log-likelihood use
            ``lnL = -0.5*(chi2 + const)`` and ``grad_lnL = -0.5*grad``.)
        """
        resp = self.clients[proc].communicate(
            "CALCULATE_CHI2_GRAD_RWA", np.asarray(params, float).ravel())
        return float(resp[0]), np.asarray(resp[1:], dtype=float)

    def residual_jacobian(self, params, proc=0):
        """Standardized residuals and their analytic Jacobian.

        ``r_i = (fit_i - data_i*n)/(cmErr_i*n)`` so ``sum(r_i**2) == chi2``.

        Returns ``(r, J)`` with ``r`` shape ``(n_res,)`` and ``J`` shape
        ``(n_res, n_params)``; columns match the non-fixed RWA parameters (the
        input ordering).  Built from the reverse-mode adjoint, so the whole
        Jacobian costs ~2 forward evaluations regardless of the parameter count
        -- for Gauss-Newton / Levenberg-Marquardt.  Energy-shift columns are
        returned as zero.
        """
        resp = self.clients[proc].communicate(
            "CALCULATE_RESIDUAL_JACOBIAN_RWA", np.asarray(params, float).ravel())
        if resp.size == 1 and resp[0] == -1.0:
            raise RuntimeError("residual_jacobian: an analytically-unsupported "
                               "segment/config is present.")
        n_res = int(round(resp[0]))
        n_cols = int(round(resp[1]))
        r = resp[2:2 + n_res]
        J = resp[2 + n_res:].reshape(n_res, n_cols)
        return r, J

    # -- calculations ---------------------------------------------------------

    def calculate_excitation_energy(self, params, proc=0):
        c = self.clients[proc]
        nsegments = int(c.communicate("UPDATE_SEGMENTS", params)[0])
        return [c.communicate("GET_EXCITATION_ENERGY", [i]) for i in range(nsegments)]

    def calculate_angles(self, params, proc=0):
        """Per-segment angles of the calculated points, one array per segment.

        The companion to :meth:`calculate_energies`: for a differential segment
        the returned angles are AZURE2's own (center-of-mass) values, which
        differ from the lab angles declared in the ``.azr`` file.

        (The underlying command is spelled ``GET_CALCUALTED_ANGLES`` -- the typo
        is in AZURE2's opcode table, not here.  This previously issued
        ``GET_DATA_ANGLES`` with the *parameter vector* in place of a segment
        index, which returned one arbitrary segment's data angles.)
        """
        c = self.clients[proc]
        nsegments = int(c.communicate("UPDATE_SEGMENTS", params)[0])
        return [c.communicate("GET_CALCUALTED_ANGLES", [i])
                for i in range(nsegments)]

    def calculate(self, params, proc=0):
        c = self.clients[proc]
        nsegments = int(c.communicate("UPDATE_SEGMENTS", params)[0])
        return [c.communicate("GET_CALCULATED_SEGMENT", [i]) for i in range(nsegments)]

    def calculate_rwa(self, params, proc=0):
        c = self.clients[proc]
        nsegments = int(c.communicate("UPDATE_SEGMENTS_RWA", params)[0])
        return [c.communicate("GET_CALCULATED_SEGMENT", [i]) for i in range(nsegments)]

    def calculate_all_rwa(self, params, proc=0):
        c = self.clients[proc]
        nsegments = int(c.communicate("UPDATE_SEGMENTS_ALL_RWA", params)[0])
        return [c.communicate("GET_CALCULATED_SEGMENT", [i]) for i in range(nsegments)]

    def calculate_energies(self, params, proc=0):
        c = self.clients[proc]
        nsegments = int(c.communicate("UPDATE_SEGMENTS", params)[0])
        return [c.communicate("GET_CALCULATED_ENERGIES", [i]) for i in range(nsegments)]

    def calculate_sfactor(self, params, proc=0):
        c = self.clients[proc]
        nsegments = int(c.communicate("UPDATE_SEGMENTS", params)[0])
        segments = [c.communicate("GET_CALCULATED_SEGMENT", [i]) for i in range(nsegments)]
        conv = [c.communicate("GET_CALCULATED_CONV", [i]) for i in range(nsegments)]
        for i in range(nsegments):
            segments[i] = segments[i] * conv[i]
        return segments

    def calculate_sfactor_rwa(self, params, proc=0):
        c = self.clients[proc]
        nsegments = int(c.communicate("UPDATE_SEGMENTS_RWA", params)[0])
        segments = [c.communicate("GET_CALCULATED_SEGMENT", [i]) for i in range(nsegments)]
        conv = [c.communicate("GET_CALCULATED_CONV", [i]) for i in range(nsegments)]
        for i in range(nsegments):
            segments[i] = segments[i] * conv[i]
        return segments

    # -- modes ----------------------------------------------------------------

    def _set_mode(self, mode_cmd):
        """Switch every instance to ``mode_cmd`` and re-INITIALIZE in parallel.

        Like spawn(), the INITIALIZE re-run is the expensive part and is fanned
        out across threads so all instances re-initialize concurrently.
        """
        def switch(c):
            c.communicate(mode_cmd, [0])
            c.communicate("INITIALIZE", [0])

        if self.nprocs == 1:
            switch(self.clients[0])
            return

        with ThreadPoolExecutor(max_workers=self.nprocs) as pool:
            list(pool.map(switch, self.clients))

    def extrap_mode(self):
        self._set_mode("SET_EXTRAP_MODE")

    def data_mode(self):
        self._set_mode("SET_DATA_MODE")
