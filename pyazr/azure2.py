"""High-level driver for one or more AZURE2 API instances.

An :class:`azure2` object owns a pool of ``nprocs`` AZURE2 subprocesses, each
paired with a client connection.  Methods that take a ``proc`` argument let a
caller dispatch work to a specific instance (e.g. for parallel chi-squared
evaluations across a worker pool).
"""

import socket
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from .client import client
from .parameters import Parameter, ParameterSet
from .server import server


def _find_free_port(start, count):
    """Return a list of ``count`` free TCP ports starting near ``start``.

    Probing avoids the silent failures that happen when a requested port is
    already in use by a stale instance.
    """
    ports = []
    candidate = start
    while len(ports) < count:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                s.bind(("", candidate))
                ports.append(candidate)
            except OSError:
                pass
        candidate += 1
        if candidate > start + count + 1000:
            raise RuntimeError("Could not find enough free ports for AZURE2.")
    return ports


class azure2:

    def __init__(self, file, nprocs=1, port=20000, binary=None,
                 verbose=False, auto_port=True):
        """Launch ``nprocs`` AZURE2 instances bound to ``file``.

        Parameters
        ----------
        file : the ``.azr`` configuration file.
        nprocs : number of parallel AZURE2 instances to spawn.
        port : base TCP port; instance ``i`` uses ``port + i`` (or the next
            free ports when ``auto_port`` is set).
        binary : path to the AZURE2 executable (auto-detected if ``None``).
        verbose : forward subprocess output to the console.
        auto_port : probe for free ports instead of assuming they are free.
        """
        self.file = file
        self.nprocs = nprocs
        self.binary = binary
        self.verbose = verbose

        # Instance-level lists.  (The original code declared `servers` as a
        # *class* attribute, so every azure2 object shared -- and leaked into
        # -- the same list.)
        self.servers = []
        self.clients = []

        if auto_port:
            self.ports = _find_free_port(port, nprocs)
        else:
            self.ports = [port + i for i in range(nprocs)]

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
        for p in self.ports:
            self.servers.append(
                server(p, self.file, binary=self.binary, verbose=self.verbose)
            )

        if self.nprocs == 1:
            # Avoid thread-pool overhead in the common single-instance case.
            # client.connect() polls until the port is accepting connections,
            # so no fixed sleep is needed here.
            self.clients = [client(port=p) for p in self.ports]
            for c in self.clients:
                c.communicate("INITIALIZE", [0])
            return

        with ThreadPoolExecutor(max_workers=self.nprocs) as pool:
            # connect() polls until each server's port is accepting
            # connections; running the polls in parallel overlaps the
            # subprocesses' startup latency.  map preserves order, so
            # self.clients stays aligned with self.ports (and thus `proc`).
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
        return self.clients[proc].communicate("GET_DATA_ANGLES", params)

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
