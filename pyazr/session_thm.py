"""The Trojan-horse experiments of one session (the ``<thm>`` ``experiment[...]`` lines): their profiled norms and
backgrounds, line shapes, vertex functions and distortion factors as the engine evaluates them.  Mixed into
:class:`pyazr.Session` (``class Session(..., SessionTHM)``); uses its ``_engine``, ``mode`` and ``params_rwa``."""

import numpy as np

from .errors import ParameterKeyError, close_names
from .locking import engine_call

_THM_BACKGROUND_TERMS = {"none": 0, "const": 1, "linear": 2, "quadratic": 3}


class SessionTHM:
    """Mixed into :class:`pyazr.Session`: the THM experiments of its project."""

    @engine_call
    def thm_experiments(self, params=None):
        """The THM experiments (``<thm>`` ``experiment[<name>]`` lines) as a
        chi-squared evaluation at ``params`` profiles them.

        Returns ``{name: report}``; see :meth:`thm_background` for a report.
        Empty without experiments (or in extrapolation mode).  Runs a
        chi-squared evaluation.
        """
        x = np.asarray(self.params_rwa if params is None else params, float)
        if self.mode == "data":
            self._engine.calculate_chi2_rwa(x)
        out = {}
        for r in self._engine.thm_experiments():
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

    @engine_call
    def thm_background(self, experiment, params=None):
        """Shared norm and background of one THM experiment at ``params``.

        The segments of ``experiment[<name>]`` share one profiled norm ``n*``
        (it multiplies the data, as ``evaluate(x).data[i].norm``) and a background
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
        scaled by n*, is ``evaluate(x).data[i].model`` plus ``b(E)``.  Runs a chi-squared
        evaluation.  KeyError for an unknown experiment.
        """
        reports = self.thm_experiments(params)
        if experiment not in reports:
            raise ParameterKeyError(f"no THM experiment {experiment!r} (have: {', '.join(reports) or 'none'})"
                                    f"{close_names(experiment, reports)}; Session.thm_experiments() lists them.")
        return reports[experiment]

    @engine_call
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
            self._engine.calculate_chi2_rwa(x)
        r = self._engine.thm_lineshape(str(experiment), np.asarray(energies, float).ravel())
        return r

    @engine_call
    def thm_vertex(self, experiment, energies, params=None, strict=False):
        """The THM entrance vertex of an experiment, averaged over its
        spectator-momentum window (``ps=`` on the experiment line).

        The off-shell x-A momentum is ``p_xA^2/2mu_xA = E + B + p_s^2/2mu_sx``
        (Mukhamedzhanov et al., PRC 96 (2017) 024623, eq. 31), so the vertex
        ``M_l = (B_c - 1) j_l(rho) - rho j_l'(rho) (+ C_l)``, ``rho = p_xA a``,
        changes with the spectator momentum p_s inside the accepted window.
        At fixed E the spectator direction fixes |p_s|, and the engine
        averages the cross section (not the amplitude) over the accepted
        directions with the fixed-E event weight |phi(q)|^2 d cos(theta_cm),
        i.e. |phi(q)|^2 q dq over the q = |p_s| the kinematics reach (not the
        E-integrated |phi|^2 p_s^2 dp_s); this reports the vertex part of it.
        See docs/source/theory/thm_implementation.rst, "Spectator-momentum
        window".

        ``energies``: c.m. energies E of x + A (MeV).  The per-level vertex
        boundaries are those of ``params`` (default: the current parameters;
        a chi-squared evaluation is run to load them in data mode).  Returns a
        dict:

        ``window``    description (``"delta"`` without a window)
        ``p_s``, ``weights``, ``T_s``   lists over E: at each E the nodes
                      (Gauss-Legendre in cos theta_cm on the directions inside
                      the window) as |p_s| (MeV/c), their normalized weights
                      and T_s = p_s^2/2mu_sx (MeV); ``theta_cm`` (deg) their
                      spectator c.m. angles to the beam.  A delta experiment
                      has one node at each E, p_s = 0 and T_s the
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
                      node.  ``p_s``, ``weights``, ``T_s`` are the nodes of
                      the line's ``ps`` window (empty lists without one).
                      With ``spectatorAngles=`` the
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
            self._engine.calculate_chi2_rwa(x)
        r = self._engine.thm_vertex(str(experiment), np.asarray(energies, float).ravel())
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

    @engine_call
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
        return self._engine.thm_distortion(str(experiment), np.asarray(energies, float).ravel())
