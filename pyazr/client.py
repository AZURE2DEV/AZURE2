"""TCP client for the AZURE2 API.

Wire protocol (matches ``api/src/AZURESocket.cpp``)::

    frame = [uint64 count][count * float64]   (native byte order)

A *request* frame carries ``[cmd, arg0, arg1, ...]`` and a *response* frame
carries the raw result values.  Because every frame is length-prefixed, there
is no fixed maximum payload size -- arbitrarily large segments no longer
overflow a buffer the way the old fixed-size protocol did.
"""

import socket
import struct
import time

import numpy as np

# Native-endian unsigned 64-bit length prefix; payload is native float64.
_LEN = struct.Struct("=Q")
_F8 = np.dtype("=f8")


class ClientError(RuntimeError):
    """Raised when the AZURE2 API connection fails or misbehaves."""


class client:

    CMD = {
        'INITIALIZE': 0,
        'UPDATE_SEGMENTS': 1,
        'GET_PARAMS': 2,
        'GET_PARAMS_NAMES': 3,
        'GET_PARAMS_ALL': 4,
        'CALCULATE_EXTERNAL_CAPTURE': 5,
        'GET_DATA_ENERGIES': 6,
        'GET_DATA_SEGMENTS': 7,
        'GET_DATA_SEGMENTS_ERRORS': 8,
        'UPDATE_DATA': 9,
        'SET_DATA_MODE': 10,
        'SET_EXTRAP_MODE': 11,
        'GET_CALCULATED_SEGMENT': 12,
        'GET_CALCULATED_ENERGIES': 13,
        'SET_CHANNEL_RADIUS': 14,
        'GET_NORMS': 15,
        'GET_NORMS_ERRORS': 16,
        'GET_CALCULATED_SEGMENT_E1': 17,
        'GET_CALCULATED_SEGMENT_E2': 18,
        'GET_DATA_CONV': 19,
        'GET_CALCULATED_CONV': 20,
        'GET_DATA_ANGLES': 21,
        'GET_CALCUALTED_ANGLES': 22,
        'GET_PARAMS_RWA': 23,
        'UPDATE_SEGMENTS_RWA': 24,
        'TRANSFORM_RWA': 25,
        'TRANSFORM_ALL_RWA': 26,
        'CALCULATE_CHI2_RWA': 27,
        'CALCULATE_CHI2': 28,
        'UPDATE_SEGMENTS_ALL_RWA': 29,
        'GET_NORM_INDICES': 30,
        'GET_SHIFT_INDICES': 31,
        'GET_PARAMS_FIXED': 32,
        'GET_DATA_EXCITATION_ENERGY': 33,
        'GET_EXCITATION_ENERGY': 34,
        'GET_PARAMS_INFO': 37,
        'CALCULATE_CHI2_GRAD_RWA': 41,
        'CALCULATE_RESIDUAL_JACOBIAN_RWA': 42,
    }

    def __init__(self, server='localhost', port=20000,
                 timeout=300.0, connect_timeout=30.0, connect_retry=0.1):
        """Connect to an AZURE2 API server.

        Parameters
        ----------
        server, port : address of the server.
        timeout : per-operation socket timeout, in seconds, once connected.
        connect_timeout : how long to keep retrying the initial connection
            while the server process is still starting up.
        connect_retry : delay between connection attempts.
        """
        self.server_address = server
        self.server_port = port
        self.timeout = timeout
        self.connect_timeout = connect_timeout
        self.connect_retry = connect_retry
        self.client_socket = None
        self.connect()

    def __del__(self):
        self.disconnect()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.disconnect()

    # -- connection management ------------------------------------------------

    def connect(self):
        """Open the connection, retrying until ``connect_timeout`` elapses.

        The server is a freshly spawned subprocess, so the port may not be
        accepting connections yet; we poll rather than guess with a fixed
        sleep.
        """
        deadline = time.monotonic() + self.connect_timeout
        last_err = None
        while True:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            try:
                sock.settimeout(max(0.5, self.connect_retry * 5))
                sock.connect((self.server_address, self.server_port))
                sock.settimeout(self.timeout)
                sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
                self.client_socket = sock
                return
            except (ConnectionRefusedError, socket.timeout, OSError) as err:
                last_err = err
                sock.close()
                if time.monotonic() >= deadline:
                    raise ClientError(
                        "Could not connect to AZURE2 API at "
                        f"{self.server_address}:{self.server_port}: {err}"
                    ) from err
                time.sleep(self.connect_retry)

    def disconnect(self):
        if self.client_socket is not None:
            try:
                self.client_socket.close()
            except OSError:
                pass
            self.client_socket = None

    # -- low-level framing ----------------------------------------------------

    def _send_all(self, payload):
        try:
            self.client_socket.sendall(payload)
        except OSError as err:
            raise ClientError(f"Failed sending to AZURE2 API: {err}") from err

    def _recv_exactly(self, nbytes):
        """Read exactly ``nbytes`` bytes, looping over short reads."""
        chunks = []
        remaining = nbytes
        while remaining > 0:
            try:
                chunk = self.client_socket.recv(remaining)
            except socket.timeout as err:
                raise ClientError("Timed out waiting for AZURE2 API "
                                  "response.") from err
            except OSError as err:
                raise ClientError(f"Failed reading from AZURE2 API: {err}") from err
            if not chunk:
                raise ClientError("AZURE2 API closed the connection "
                                  "unexpectedly.")
            chunks.append(chunk)
            remaining -= len(chunk)
        return b"".join(chunks)

    # -- public protocol ------------------------------------------------------

    def send(self, cmd, data):
        """Send a request frame ``[cmd, *data]``."""
        payload = np.empty(1 + len(data), dtype=_F8)
        payload[0] = cmd
        if len(data):
            payload[1:] = np.asarray(data, dtype=_F8)
        body = payload.tobytes()
        self._send_all(_LEN.pack(len(payload)) + body)

    def receive(self):
        """Read a response frame and return it as a float64 numpy array."""
        (count,) = _LEN.unpack(self._recv_exactly(_LEN.size))
        if count == 0:
            return np.empty(0, dtype=_F8)
        raw = self._recv_exactly(count * _F8.itemsize)
        # Copy so the returned array owns its memory (frombuffer is read-only).
        return np.frombuffer(raw, dtype=_F8).copy()

    def communicate(self, cmd, data, type=None):
        """Send ``cmd`` with ``data`` and return the response array.

        ``type`` is accepted for backwards compatibility and ignored: the
        protocol is always float64.
        """
        if cmd not in self.CMD:
            raise KeyError(f"Unknown AZURE2 command: {cmd!r}")
        self.send(self.CMD[cmd], data)
        return self.receive()
