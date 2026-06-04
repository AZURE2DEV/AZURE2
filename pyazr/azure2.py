"""High-level driver for one or more AZURE2 API instances.

An :class:`azure2` object owns a pool of ``nprocs`` AZURE2 subprocesses, each
paired with a client connection.  Methods that take a ``proc`` argument let a
caller dispatch work to a specific instance (e.g. for parallel chi-squared
evaluations across a worker pool).
"""

import socket

import numpy as np

from .client import client
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
        """Start the subprocesses, connect clients, and initialize them."""
        for p in self.ports:
            self.servers.append(
                server(p, self.file, binary=self.binary, verbose=self.verbose)
            )

        # client.connect() polls until the port is accepting connections, so
        # no fixed sleep is needed here.
        for p in self.ports:
            self.clients.append(client(port=p))

        for c in self.clients:
            c.communicate("INITIALIZE", [0])

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

    def extrap_mode(self):
        for c in self.clients:
            c.communicate("SET_EXTRAP_MODE", [0])
            c.communicate("INITIALIZE", [0])

    def data_mode(self):
        for c in self.clients:
            c.communicate("SET_DATA_MODE", [0])
            c.communicate("INITIALIZE", [0])
