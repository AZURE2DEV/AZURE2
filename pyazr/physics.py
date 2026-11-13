"""Physics functions of the engine as free functions on arrays (docs/upgrade/ARCHITECTURE.md section 8).

Each one calls the C++ function of ``physics/`` through ``_azure2`` element-wise; no session or project is needed.
"""
import numpy as np

from .azure2 import _azure2
from .errors import DomainError, EngineError

__all__ = ["whittaker_wave"]

def whittaker_wave(l, eta_b, rho):
    """Whittaker function of a bound channel, W_l(rho) = W_{-eta_b, l+1/2}(2 rho).

    The decaying solution of u'' - [1 + 2 eta_b/rho + l(l+1)/rho^2] u = 0, normalized so that
    W ~ exp(-rho) (2 rho)^(-eta_b) for large rho (unit ANC), computed without GSL
    (docs/upgrade/coulomb.md section 1.4; relative error below 1e-12 against mpmath on l = 0-10,
    eta_b = 0-1000, rho = 1e-4-200).  In dimensionless variables, unlike :func:`pyazr.whittaker`, which takes the
    pair, radius and binding energy (a Python port of the engine's GSL path).

    Parameters
    ----------
    l : orbital angular momentum (integer >= 0).
    eta_b : bound-state Sommerfeld parameter Z1 Z2 e^2 mu / (hbar^2 kappa), >= 0 (dimensionless).
    rho : kappa r > 0 (dimensionless).  Arrays broadcast against each other.

    Returns
    -------
    dict of arrays of the broadcast shape: ``mantissa`` and ``log10_scale`` (W = mantissa * 10**log10_scale, kept
    apart because W leaves the double range just below a charged threshold), ``log10_W``, ``W`` (0 or inf where W
    is outside the double range), ``rho_dlog`` (rho W'/W, primes d/drho: the bound-channel shift function at
    rho = kappa a).

    Outside the double range W is returned as 0 or inf in ``W`` by design (``mantissa`` and ``log10_scale`` keep
    the value); that is not an error.

    Raises
    ------
    pyazr.DomainError
        (a ``ValueError``) if l is not an integer in 0-10000, or an element of eta_b is negative or not finite, or of
        rho not positive and finite (the first such element is named, with the domain).
    pyazr.EngineError
        (a ``RuntimeError``) if the computation fails for an element (not seen on the tested range).
    """
    L, E, R = np.broadcast_arrays(np.asarray(l), np.asarray(eta_b, float), np.asarray(rho, float))
    if not np.all(np.equal(np.mod(L, 1), 0)):
        i = int(np.flatnonzero(~np.equal(np.mod(L, 1), 0).ravel())[0])
        raise DomainError(f"whittaker_wave: l must be an integer in 0-10000, got {float(L.ravel()[i])!r} at "
                          f"element {i} (help(pyazr.whittaker_wave))")
    shape = L.shape
    res = _azure2.whittaker(L.astype(np.int32).ravel(), E.ravel(), R.ravel())
    status = np.asarray(res["status"])
    for code, kind in ((1, DomainError), (2, EngineError)):
        bad = np.flatnonzero(status == code)
        if bad.size:
            i = bad[0]
            domain = ("; the domain: l an integer in 0-10000, eta_b finite >= 0, rho finite > 0 "
                      "(help(pyazr.whittaker_wave); pyazr.kinds()['physics']['whittaker'])" if code == 1 else
                      "; please report it with these arguments")
            raise kind(f"whittaker_wave: {'invalid argument' if code == 1 else 'computation failed'} at element {i} "
                       f"(l = {int(L.ravel()[i])}, eta_b = {float(E.ravel()[i])!r}, rho = {float(R.ravel()[i])!r}); "
                       f"{bad.size} element(s) in all{domain}")
    mantissa = np.asarray(res["mantissa"], float).reshape(shape)
    scale = np.asarray(res["log10_scale"]).reshape(shape)
    log10_W = np.log10(mantissa) + scale
    with np.errstate(over="ignore", under="ignore"):
        W = np.power(10.0, log10_W)
    return {"mantissa": mantissa, "log10_scale": scale, "log10_W": log10_W, "W": W,
            "rho_dlog": np.asarray(res["rho_dlog"], float).reshape(shape)}
