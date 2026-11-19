"""Physics functions of the engine's libraries as free functions on arrays (docs/upgrade/pyazr3.md section 2.7).

One Python function per row of the binding's table (``api/physics_table/src/PhysicsTable.cpp``), built here at import from
``_azure2.physics_functions()``: its signature, docstring, units and domains are the row's, and each call evaluates
the row's library function element by element in C++ (no session, no project, no physics in Python).

Calling convention, the same for every function:

- arguments broadcast against each other (numpy rules); scalars in give scalars out, arrays give arrays of the
  broadcast shape;
- a function with one returned field returns it (a float or an array); one with several returns a dict of them;
- ``errors="raise"`` (the default): an element outside the domain raises :class:`pyazr.DomainError` naming the
  argument, its value, the element and the allowed range; a computation that does not converge raises
  :class:`pyazr.EngineError`;
- ``errors="mask"``: nothing is raised; such elements are NaN and the call returns ``(result, status)``, ``status`` an
  integer array of the broadcast shape: 0 computed, k > 0 argument k (1-based, in the order of the signature) outside
  its domain, -2 the arguments together outside the domain (each in its own: a lab angle beyond the kinematic
  maximum, an uncharged pair), -1 not converged.  :func:`status_text` turns codes into words.

``functions`` maps each name to its function; ``table()`` is the rows as plain data (``pyazr.kinds()["physics"]``
is the same), with the argument a plot puts on its horizontal axis (``axis``) and every returned field's unit.
A new function is a row of the binding's table and a test naming it (``pyazr/CLAUDE.md``, "Add a physics function").

Two functions are not rows, because their arguments or answers are not one number per element:
:func:`collision_matrix`, the level matrix and the collision matrix of one J^pi block at one energy from the arrays
of its levels and channels (``physics/CollisionMatrix.h``), and :func:`optical_scattering`, the partial-wave S-matrix,
the elastic angular distribution and the reaction cross section of a pair in an optical potential
(``physics/transfer``).
"""
import copy
import inspect

import numpy as np

from . import _azure2
from .errors import DomainError, EngineError

STATUS = {0: "computed", -1: "not converged", -2: "the arguments together are outside the domain"}
_ROWS = {row["name"]: row for row in _azure2.physics_functions()}


def table():
    """The rows of the physics table as plain data: a dict by name of ``summary``, ``doc``, ``library`` (the C++
    function called), ``axis``, ``arguments`` (name, unit, doc, kind, domain in words, lower, upper, lower_open,
    upper_open, choices, default) and ``returns`` (name, unit, doc).  A fresh copy at each call."""
    return copy.deepcopy(_ROWS)


def status_text(name, status):
    """Words for status codes of function ``name`` (an int or an array of them): ``"computed"``, ``"e_cm outside
    its domain (finite > 0)"``, ``"not converged"``, ...  Returns a str or a list of str."""
    row = _ROWS[name]

    def one(code):
        code = int(code)
        if code > 0:
            arg = row["arguments"][code - 1]
            return f"{arg['name']} outside its domain ({arg['domain']})"
        return STATUS.get(code, f"unknown status {code}")
    return one(status) if np.ndim(status) == 0 else [one(c) for c in np.ravel(status)]


def _value_text(value, unit):
    return f"{value!r}" + (f" {unit}" if unit else "")


def _raise_for(row, status, columns, shape):
    """The error of the first element whose status is not 0 (columns: the flat arguments)."""
    name = row["name"]
    bad = np.flatnonzero(status != 0)
    i = int(bad[0])
    code = int(status[i])
    where = f"element {i}" + (f" (index {np.unravel_index(i, shape)})" if len(shape) > 1 else "")
    count = f"{bad.size} of {status.size} element(s) are not computed"
    fix = f"help(pyazr.physics.{name}) gives the domain; errors='mask' returns NaN with a status per element instead"
    if code > 0:
        arg = row["arguments"][code - 1]
        value = float(columns[code - 1][i])
        if arg["kind"] == "choice":
            value = arg["choices"][int(value)] if value == int(value) and 0 <= value < len(arg["choices"]) else value
        raise DomainError(f"{name}: {arg['name']} = {_value_text(value, arg['unit'])} at {where} is outside the "
                          f"domain; allowed: {arg['name']} {arg['domain']} ({arg['doc']}); {count}.  {fix}.",
                          what=f"{arg['name']} = {_value_text(value, arg['unit'])} at {where} is outside the domain",
                          allowed=f"{arg['name']} {arg['domain']} ({arg['doc']})", fix=fix)
    given = ", ".join(f"{a['name']} = {_value_text(float(columns[k][i]), a['unit'])}"
                      for k, a in enumerate(row["arguments"]))
    if code == -2:
        raise DomainError(f"{name}: the arguments at {where} ({given}) are outside the domain together, though each "
                          f"is in its own: {row['doc']}  {count}.  {fix}.",
                          what=f"the arguments at {where} ({given}) are outside the domain together",
                          allowed=row["doc"], fix=fix)
    raise EngineError(f"{name}: the computation did not converge at {where} ({given}); {count}.  Please report it "
                      f"with these arguments; errors='mask' gives the other elements.")


def _argument_column(row, k, value):
    """One argument as a float array: a choice by its name, a flag from a bool."""
    arg = row["arguments"][k]
    if arg["kind"] == "choice" and isinstance(value, str):
        if value not in arg["choices"]:
            raise DomainError(f"{row['name']}: {arg['name']} = {value!r} is not known; allowed: {arg['name']} "
                              f"{arg['domain']}.  help(pyazr.physics.{row['name']}).",
                              what=f"{arg['name']} = {value!r} is not known", allowed=f"{arg['name']} {arg['domain']}",
                              fix=f"help(pyazr.physics.{row['name']})")
        return np.asarray(float(arg["choices"].index(value)))
    try:
        return np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{row['name']}: {arg['name']} must be a number or an array of numbers "
                        f"({arg['unit'] or 'no unit'}, {arg['domain']}), got {type(value).__name__}") from None


def _call(row, values, errors):
    if errors not in ("raise", "mask"):
        raise ValueError(f"{row['name']}: errors = {errors!r}; allowed: 'raise' (DomainError, EngineError) or 'mask' "
                         f"(NaN and a status per element)")
    arrays = [_argument_column(row, k, v) for k, v in enumerate(values)]
    scalar = all(a.ndim == 0 for a in arrays)
    broadcast = np.broadcast_arrays(*arrays)
    shape = broadcast[0].shape
    columns = [np.ascontiguousarray(a, dtype=float).ravel() for a in broadcast]
    out, status = _azure2.physics_evaluate(row["name"], columns)
    if errors == "raise" and np.any(status != 0):
        _raise_for(row, status, columns, shape)
    fields = {}
    for j, ret in enumerate(row["returns"]):
        fields[ret["name"]] = float(out[j][0]) if scalar else out[j].reshape(shape)
    result = fields[row["returns"][0]["name"]] if len(fields) == 1 else fields
    if errors == "mask":
        return result, (int(status[0]) if scalar else status.reshape(shape))
    return result


def _docstring(row):
    lines = [row["summary"], "", row["doc"], "", "Parameters", "----------"]
    for a in row["arguments"]:
        unit = f", {a['unit']}" if a["unit"] else ""
        default = f", default {a['default']!r}" if a["default"] is not None else ""
        lines.append(f"{a['name']} : {a['kind']}{unit}{default}")
        lines.append(f"    {a['doc']}; domain: {a['domain']}.")
    lines.append("errors : 'raise' (default) or 'mask'")
    lines.append("    'mask': NaN outside the domain, and the call returns (result, status) (help(pyazr.physics)).")
    lines += ["", "Returns", "-------"]
    if len(row["returns"]) > 1:
        lines.append("dict of arrays of the broadcast shape (floats for scalar arguments):")
    for r in row["returns"]:
        lines.append(f"{r['name']} : {r['unit'] or 'dimensionless'}")
        lines.append(f"    {r['doc']}")
    lines += ["", "Raises", "------",
              "pyazr.DomainError", "    (a ValueError) an element outside the domain, named with its value.",
              "pyazr.EngineError", "    (a RuntimeError) a computation that does not converge.",
              "", f"Library: {row['library']} (element-wise in C++ through _azure2.physics_evaluate)."]
    return "\n".join(lines)


def _make(row):
    params = [inspect.Parameter(a["name"], inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                default=inspect.Parameter.empty if a["default"] is None else a["default"])
              for a in row["arguments"]]
    params.append(inspect.Parameter("errors", inspect.Parameter.KEYWORD_ONLY, default="raise"))
    signature = inspect.Signature(params)

    def function(*args, **kwargs):
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        return _call(row, [bound.arguments[a["name"]] for a in row["arguments"]], bound.arguments["errors"])

    function.__name__ = function.__qualname__ = row["name"]
    function.__module__ = __name__
    function.__doc__ = _docstring(row)
    function.__signature__ = signature
    return function


def collision_matrix(level_energies, reduced_widths, energy, *, penetrability, shift, boundary=None,
                     hard_sphere_phase, coulomb_phase=None, parametrization="formal", shift_at_levels=None,
                     shift_derivative_at_levels_per_MeV=None, park_norms=None, photon_channels=None, use_rmc=False):
    """The level matrix and the collision matrix of one J^pi block at one energy, without a session or a project.

    The block is n levels and nc channels.  With formal (standard R-matrix) levels and L_c = S_c - B_c + i P_c the
    inverse level matrix is ``M[l, m] = (E_l - E) delta_lm - sum_c gamma_lc gamma_mc L_c``, ``A = M^-1``, and

        U[c, c'] = Omega_c Omega_c' (delta_cc' + 2i sqrt(P_c P_c') gamma_c^T A gamma_c'),
        Omega_c = exp(i (omega_c - phi_c)),      T[c, c'] = exp(2i omega_c) delta_cc' - U[c, c'],

    the engine's own arithmetic (``physics/CollisionMatrix.h``: ``collision_matrix``; one factorization of M, one
    solve per channel).  The channel functions at ``energy`` are the caller's: for a pair of the rows of this
    module, ``penetrability``, ``shift``, ``hard_sphere_phase`` and ``coulomb_phase(l, eta) - coulomb_phase(0, eta)``.

    ``parametrization`` says what ``level_energies`` and ``reduced_widths`` are, and with it the form of M (U and T
    follow from A as above in each).  With S_lc = S_c at the energy of level l (``shift_at_levels``), over the
    particle channels:

    - ``"formal"``: the matrix above, with the constants B_c (``boundary``).
    - ``"brune"`` (Brune, Phys. Rev. C 66, 044611 (2002), Eq. (33)): the observed energies and Brune's amplitudes,
      ``M[l, l] = E_l - E - sum_c gamma_lc^2 (S_c - S_lc + i P_c)`` and, for l != m,
      ``M[l, m] = -sum_c gamma_lc gamma_mc (S_c + i P_c - (S_lc (E - E_m) - S_mc (E - E_l)) / (E_l - E_m))``.
    - ``"park"`` (Park, Phys. Rev. C 104, 064612 (2021), Eqs. (21), (22), (28)): the observed energies and the
      amplitudes of the observed widths, Gamma_lc = 2 P_c gamma_lc^2,
      ``M[l, m] = (E_l - E) J_lm - sum_c gamma_lc gamma_mc (S_c - S_lc + i P_c)`` with
      ``J_lm = -sum_c gamma_lc gamma_mc (S_lc - S_mc) / (E_l - E_m)`` for l != m and J_ll = ``park_norms[l]``,
      which is ``1 - sum_c gamma_lc^2 dS_c/dE`` at E_l (:meth:`pyazr.Session.park_norms` for a session's levels).

    Two levels closer than 1e-6 MeV are coincident: the divided differences above are then their Taylor forms with
    the mean of the two derivatives dS_c/dE (``physics/LevelPairShift.h``), and
    ``shift_derivative_at_levels_per_MeV`` is required.

    A photon channel (``photon_channels``) has L_c = 0: its widths do not enter M, and it has no shift function,
    boundary value or phase (each must be 0 there).  Its ``penetrability`` is the caller's factor between gamma^2
    and the width; the engine's is ``(E_gamma / (hbar c))^(2L+1)`` (fm^-(2L+1)) for multipolarity L, with
    Gamma = 2 P gamma^2.  ``use_rmc`` adds the Reich-Moore term ``-i gamma_lc^2`` of every photon channel to
    ``M[l, l]`` (the option ``use_rmc`` of a session, which then takes P = 1 for a photon channel).

    Parameters
    ----------
    level_energies : array (n,), MeV
        E_l, on the energy scale of ``energy`` (the engine's is the excitation energy of the compound nucleus):
        formal energies, or the observed ones (Brune, Park).
    reduced_widths : array (n, nc), MeV^1/2
        The reduced-width amplitudes gamma_lc of the parametrization (a photon channel: in the unit that makes
        ``penetrability * gamma^2`` an energy in MeV).
    energy : float, MeV
        E.
    penetrability, shift, hard_sphere_phase : arrays (nc,)
        P_c (0 for a closed channel), S_c and phi_c (rad) of each channel at ``energy``.
    boundary : array (nc,), ``"formal"`` only
        B_c.  Not an argument of the other two parametrizations, whose matrices do not contain it.
    coulomb_phase : array (nc,), rad, optional
        omega_c = sigma_l - sigma_0 of each channel; zero (neutral pairs) when not given.
    parametrization : ``"formal"`` (default), ``"brune"`` or ``"park"``
    shift_at_levels : array (n, nc), ``"brune"`` and ``"park"``
        S_c at the energy of each level (the channel energy that E_l is, in the centre of mass of that channel's
        pair), dimensionless; 0 for a photon channel.
    shift_derivative_at_levels_per_MeV : array (n, nc), MeV^-1, optional
        dS_c/dE at the energy of each level.  Read only for coincident levels.
    park_norms : array (n,), ``"park"`` only
        J_ll, dimensionless.
    photon_channels : array (nc,) of bool, optional
        True for a photon channel; every channel is a particle channel when not given.
    use_rmc : bool
        The Reich-Moore term (``"formal"`` only, with at least one photon channel).

    Returns
    -------
    dict
        ``A`` (n, n), ``U`` (nc, nc) and ``T`` (nc, nc), complex.  The cross sections follow from T
        (``docs/source/developer/level_matrix.rst``, "Cross Section Calculation").

    Raises
    ------
    ValueError
        A shape that does not fit, a value that is not finite, a negative penetrability, an unknown
        parametrization, an array a parametrization needs and did not get or got and does not read, a shift or
        phase of a photon channel that is not 0.
    pyazr.EngineError
        The inverse level matrix is singular at this energy.
    """
    def array(values):
        return None if values is None else np.ascontiguousarray(values, dtype=float)

    widths = np.ascontiguousarray(reduced_widths, dtype=float)
    channels = widths.shape[1] if widths.ndim == 2 else 0
    omega = np.zeros(channels) if coulomb_phase is None else coulomb_phase
    photons = None if photon_channels is None else np.ascontiguousarray(photon_channels, dtype=bool)
    out = _azure2.collision_matrix(array(level_energies), widths, float(energy), array(penetrability), array(shift),
                                   array(boundary), array(hard_sphere_phase), array(omega), str(parametrization),
                                   array(shift_at_levels), array(shift_derivative_at_levels_per_MeV),
                                   array(park_norms), photons, bool(use_rmc))
    if out.pop("singular"):
        raise EngineError(f"collision_matrix: the inverse level matrix is singular at energy = {float(energy)!r} MeV")
    return out


def optical_scattering(z1, z2, m1_u, m2_u, e_cm, potential, *, angles_deg=(), projectile_spin=0.0, l_max=20,
                       step_fm=0.02, matching_radius_fm=20.0):
    """Elastic scattering of a pair in an optical potential at one energy, without a session or a project.

    The projectile has spin 0, 1/2 or 1; the target's spin is ignored.  The potential at the distance r (fm) is

        U(r) = V_C(r) - V f_V(r) - i W f_W(r) - i 4 W_D f_D(r) (1 - f_D(r))
               + (V_so + i W_so) (2 fm^2) (1/r) df_so/dr <2 l.s>,       f(r) = 1 / (1 + exp((r - R) / a)),

    in MeV (V, W, W_D > 0 attract and absorb; ``physics/transfer/OpticalPotential.h`` Eq. (2.1), FRESCO's forms;
    V_C the potential of a uniformly charged sphere of radius R_C).  The radial equation of each partial wave is
    integrated by Numerov's method from the origin to ``matching_radius_fm`` and matched there to Coulomb functions
    (``physics/transfer/DistortedWaves.h``, checked against FRESCO: ``docs/upgrade/transfer.md``); the kinematics
    are non-relativistic.  These are the waves the engine's THM distortion uses.

    Parameters
    ----------
    z1, z2 : int
        Charge numbers of the two nuclei.
    m1_u, m2_u : float, u
        Their masses.
    e_cm : float, MeV
        Kinetic energy of the pair in its centre of mass, > 0.
    potential : dict
        ``radius_convention``: ``"absolute"`` (R = r0, the default), ``"target_mass"`` (R = r0 A_t^1/3) or
        ``"both_masses"`` (R = r0 (A_p^1/3 + A_t^1/3)); ``target_mass_number`` A_t and
        ``projectile_mass_number`` A_p of that convention; ``coulomb_r0_fm`` (the parameter of R_C, fm; 0: a point
        charge); and the terms ``real_volume`` (V), ``imaginary_volume`` (W), ``imaginary_surface`` (W_D),
        ``real_spin_orbit`` (V_so), ``imaginary_spin_orbit`` (W_so), each ``(depth in MeV, r0 in fm, diffuseness
        a in fm)``.  A member left out is zero.
    angles_deg : array, deg, optional
        Scattering angles in the centre of mass, 0 < theta <= 180.
    projectile_spin : 0, 0.5 or 1
    l_max : int
        The highest orbital angular momentum solved.
    step_fm, matching_radius_fm : float, fm
        The step of the radial grid and its last node.

    Returns
    -------
    dict
        ``k`` (fm^-1) and ``eta`` (the wave number and the Sommerfeld parameter); ``coulomb_phase`` (sigma_l, rad,
        l = 0 ... l_max); ``l``, ``two_j`` and ``S``: one entry per partial wave, S the nuclear S-matrix element
        (the Coulomb phase exp(2i sigma_l) is not in it; one wave per l with ``two_j`` = 2l when the potential has
        no spin-orbit term); ``ratio_to_rutherford`` at ``angles_deg`` (dimensionless: the elastic cross section
        over Rutherford's, averaged over the projectile's spin projections) or, for a neutral pair,
        ``dsigma_domega_mb_sr`` (mb/sr, c.m.); ``reaction_mb`` (the reaction cross section, mb); ``spin_orbit``
        (bool); ``description`` (the potential and its radii in words).

    Raises
    ------
    ValueError
        A member the potential does not have, a term that is not three numbers or has a depth without a positive
        radius and diffuseness, an energy <= 0, a spin other than 0, 1/2, 1, an angle outside (0, 180], a grid
        that cannot hold the matching.
    pyazr.EngineError
        A Coulomb function at the matching radius failed.
    """
    two_s = int(round(2 * float(projectile_spin)))
    if abs(two_s - 2 * float(projectile_spin)) > 1e-12:
        raise ValueError(f"optical_scattering: projectile_spin = {projectile_spin!r} is not 0, 0.5 or 1")
    try:
        return _azure2.optical_scattering(int(z1), int(z2), float(m1_u), float(m2_u), float(e_cm), dict(potential),
                                          np.atleast_1d(np.ascontiguousarray(angles_deg, dtype=float)), two_s,
                                          int(l_max), float(step_fm), float(matching_radius_fm))
    except RuntimeError as err:
        raise EngineError(f"optical_scattering: {err}") from err


functions = {name: _make(row) for name, row in _ROWS.items()}
globals().update(functions)
__all__ = sorted(functions) + ["functions", "table", "status_text", "STATUS", "collision_matrix", "optical_scattering"]
