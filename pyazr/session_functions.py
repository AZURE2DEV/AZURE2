"""The session's diagnostic calls described as rows (``Session.describe()["functions"]``): :meth:`Session.channel_functions`,
:meth:`Session.potential_box_states`, :meth:`Session.potential_bound_state`,
:meth:`Session.potential_scattering_length`, :meth:`Session.potential_depth_for`, :meth:`Session.jpi_block` and the
five THM calls of :mod:`pyazr.session_thm`.

Each row has the shape of a row of the physics table (``pyazr.kinds()["physics"]``): ``name``, ``summary``, ``doc``,
``axis`` (the argument a plot puts on its horizontal axis, None for none), ``at_x`` (the answer depends on the fit
vector ``params``), ``arguments`` (name, kind, unit, doc, domain, lower, upper, default, example, sweep) and
``returns``.  A returned field has ``name``, ``unit``, ``doc``, ``path`` (where it sits in the answer: ``"E"`` a key,
``"exits[].levels[].NC2"`` a key of each element of a list of a list, ``"{experiment}.chi2"`` a key under each
experiment's name) and ``shape``: ``"axis"`` (an array over the axis argument), ``"axis x nodes"`` (per axis value a
list of node values), ``"scalar"``, ``"text"``, ``"list"`` or ``"matrix"``.  A field under a list has ``lines``: what
one element is and the fields that name it (``{"over": "exits[].levels[]", "label": ["pair", "level", ...]}``), so
that a plot draws one line per element without knowing the call.  The answers themselves are not changed by these
rows; ``tests/pyazr/session_functions_test.py`` checks on THM projects that every path of a row is in the answer and
every key of the answer is a path of its row.

A new diagnostic of the session is a method (where it computes) and one row here, with its case in that test.
"""
import copy

__all__ = ["SESSION_FUNCTIONS", "session_functions"]

_ARG = {"unit": "", "lower": None, "upper": None, "default": None, "example": None, "sweep": None}
_PARAMS = dict(_ARG, name="params", kind="vector", doc="the fit vector x (default: Session.parameters.x0)",
               domain="a vector of the session's length")
_EXPERIMENT = dict(_ARG, name="experiment", kind="experiment", doc="a THM experiment of the project",
                   domain="the name of an experiment[<name>] of <thm>")
_ENERGIES = dict(_ARG, name="energies", kind="real", unit="MeV", doc="c.m. energy E of x + A", domain="finite",
                 example=0.5, sweep={"lower": 0.05, "upper": 2.0, "spacing": "linear"})


def _r(path, unit, doc, shape="axis", lines=None):
    """A returned field: its name is the last part of its path."""
    out = {"name": path.rsplit(".", 1)[-1], "unit": unit, "doc": doc, "path": path, "shape": shape}
    if lines:
        out["lines"] = lines
    return out


def _row(name, summary, doc, arguments, returns, *, axis=None, at_x=True):
    return {"name": name, "summary": summary, "doc": doc, "axis": axis, "at_x": at_x,
            "library": f"pyazr: Session.{name}", "arguments": [dict(a) for a in arguments], "returns": list(returns)}


_EXITS = {"over": "exits[]", "label": ["pair"]}
_EXIT_LEVELS = {"over": "exits[].levels[]", "label": ["pair", "jgroup", "level", "J", "pi", "E_level", "Gamma"]}
_CHANNELS = {"over": "channels[]", "label": ["jgroup", "channel", "J", "pi", "l", "s"]}
_CHANNEL_LEVELS = {"over": "channels[].levels[]", "label": ["jgroup", "channel", "J", "pi", "l", "s", "level"]}
_EXPERIMENTS = {"over": "{experiment}", "label": ["experiment"]}
_BLOCK_LEVELS = {"over": "levels[]", "label": ["level", "E_level"]}
_BLOCK_CHANNELS = {"over": "channels[]", "label": ["channel", "pair", "radiation", "l", "s"]}
_BLOCK_ELEMENTS = {"over": "elements[]", "label": ["channel", "channel_out"]}
_BLOCK_PHASES = {"over": "phase_shifts[]", "label": ["channel", "pair", "l", "s"]}
_BLOCK_LEVEL_PAIRS = {"over": "level_pairs[]", "label": ["level", "level_out"]}


def _profile(prefix, lines):
    return [_r(f"{prefix}segments", "", "the keys of its data sets", "list", lines),
            _r(f"{prefix}background", "", "none, const, linear or quadratic", "text", lines),
            _r(f"{prefix}points", "", "data points", "scalar", lines),
            _r(f"{prefix}chi2", "", "the experiment's chi-squared", "scalar", lines),
            _r(f"{prefix}status", "", "profiled, or what a degenerate profile fell back to", "text", lines),
            _r(f"{prefix}norm", "", "the profiled norm n* (multiplies the data)", "scalar", lines),
            _r(f"{prefix}b", "model, model/MeV, model/MeV^2", "the background coefficients b0, b1, b2 (0 to 3 of them)",
               "list", lines),
            _r(f"{prefix}cov", "", "covariance of (norm, b0, ...) from the profile, at fixed R-matrix parameters",
               "matrix", lines),
            _r(f"{prefix}sigma", "", "sqrt of the covariance's diagonal", "list", lines),
            _r(f"{prefix}cbkg", "", "the coherent background's fit parameters, name -> value (with cbackground=)",
               "list", lines)]


SESSION_FUNCTIONS = {
    "channel_functions": _row(
        "channel_functions", "The session's channel functions of a pair: Coulomb functions, P, S and the hard-sphere "
        "phase.",
        "As the session computes them: the Coulomb table when one is in use, the hybrid model's "
        "Numerov solution for a pair with a nuclear potential.  Where open is False every field but S is 0 and S is "
        "the closed-channel shift function at the channel radius.",
        [dict(_ARG, name="pair", kind="pair", doc="particle pair key of <levels>", domain="a particle pair of the "
              "project"),
         dict(_ARG, name="e_cm", kind="real", unit="MeV", doc="c.m. energy of the pair", domain="finite",
              example=1.0, sweep={"lower": 0.01, "upper": 10.0, "spacing": "log"}),
         dict(_ARG, name="l", kind="integer", doc="orbital angular momentum", domain="integer >= 0", lower=0.0,
              default=0, example=0, sweep={"lower": 0.0, "upper": 10.0, "spacing": "linear"}),
         dict(_ARG, name="radius", kind="real", unit="fm", doc="evaluation radius; 0: the pair's channel radius",
              domain="finite >= 0", lower=0.0, default=0.0, sweep={"lower": 1.0, "upper": 10.0, "spacing": "linear"}),
         dict(_ARG, name="J", kind="real", doc="the channel's total J: required (and only allowed) for a pair whose "
              "potential has a spin-orbit term", domain="a multiple of 1/2 >= 0, or None")],
        [_r("F", "", "regular Coulomb function"), _r("Fp", "", "dF/drho"), _r("G", "", "irregular Coulomb function"),
         _r("Gp", "", "dG/drho"), _r("P", "", "penetrability"), _r("S", "", "shift function"),
         _r("delta_hs", "rad", "hard-sphere phase shift"),
         _r("open", "", "E > 0 and the engine's functions defined there (bool)"), _r("e_cm", "MeV", "the energies")],
        axis="e_cm", at_x=False),
    "potential_box_states": _row(
        "potential_box_states", "The box states of a pair's potential in one partial wave: the poles E_p, gamma_p^2 of "
        "its interior R-function.",
        "The eigenvalues of the radial equation of the potential inside the radius with r u'/u = boundary there, their "
        "node numbers (what `subtract` names) and residues gamma_p^2 (docs/upgrade/potential_background.md sections 1 "
        "and 3).",
        [dict(_ARG, name="pair", kind="pair", doc="particle pair key of <levels> (a pair with a potential)",
              domain="a particle pair of the project"),
         dict(_ARG, name="l", kind="integer", doc="orbital angular momentum", domain="integer >= 0", lower=0.0,
              default=0, example=0),
         dict(_ARG, name="J", kind="real", doc="the channel's total J: required when the potential depends on J, "
              "ignored when it does not", domain="a multiple of 1/2 >= 0, or None", example=0.5),
         dict(_ARG, name="radius", kind="real", unit="fm", doc="box radius; 0: the pair's channel radius",
              domain="finite >= 0", lower=0.0, default=0.0),
         dict(_ARG, name="boundary", kind="real", doc="boundary-condition constant B", domain="finite", default=0.0),
         dict(_ARG, name="e_min", kind="real", unit="MeV", doc="lower end of the energy window (c.m.)",
              domain="finite", default=-60.0),
         dict(_ARG, name="e_max", kind="real", unit="MeV", doc="upper end of the energy window (c.m.)",
              domain="finite > e_min", default=20.0)],
        [_r("nodes", "", "node numbers of the box states", "list"), _r("E", "MeV", "their energies E_p (c.m.)", "list"),
         _r("gamma2", "MeV", "their residues gamma_p^2", "list")],
        at_x=False),
    "potential_bound_state": _row(
        "potential_bound_state", "The bound state of a pair's potential in one partial wave nearest an energy, and its "
        "ANC.",
        "The eigenvalue of the potential's radial equation from 0 to infinity nearest the energy within the window, "
        "and the ANC of the normalised state: what `capture interior` compares a final level with "
        "(docs/upgrade/potential_background.md section 5.3).",
        [dict(_ARG, name="pair", kind="pair", doc="particle pair key of <levels> (a pair with a potential)",
              domain="a particle pair of the project"),
         dict(_ARG, name="l", kind="integer", doc="orbital angular momentum", domain="integer >= 0", lower=0.0,
              default=0, example=1),
         dict(_ARG, name="J", kind="real", doc="the channel's total J: required when the potential depends on J, "
              "ignored when it does not", domain="a multiple of 1/2 >= 0, or None", example=0.5),
         dict(_ARG, name="energy", kind="real", unit="MeV", doc="c.m. energy below threshold near which to look",
              domain="finite < 0", upper=0.0, default=-1.0, example=-1.5),
         dict(_ARG, name="window", kind="real", unit="MeV", doc="half-width of the search", domain="finite > 0",
              lower=0.0, default=2.0)],
        [_r("found", "", "an eigenvalue in the window (bool)", "scalar"),
         _r("E", "MeV", "the eigenvalue (c.m.)", "scalar"), _r("anc", "fm^-1/2", "the ANC of the normalised state",
                                                               "scalar")],
        at_x=False),
    "potential_scattering_length": _row(
        "potential_scattering_length", "The scattering length and effective range of a pair's potential alone in one "
        "partial wave.",
        "The effective-range expansion of the potential's own phase shifts near threshold: Coulomb-modified for a "
        "charged pair (l = 0 only), k^(2l+1) cot delta for a neutral one, a cubic in k^2 fitted on a window set by "
        "the potential's range (docs/upgrade/potential_background.md section 13.1).",
        [dict(_ARG, name="pair", kind="pair", doc="particle pair key of <levels> (a pair with a potential)",
              domain="a particle pair of the project"),
         dict(_ARG, name="l", kind="integer", doc="orbital angular momentum (0 for a charged pair)",
              domain="integer >= 0", lower=0.0, default=0, example=0),
         dict(_ARG, name="J", kind="real", doc="the channel's total J: required when the potential depends on J, "
              "ignored when it does not", domain="a multiple of 1/2 >= 0, or None", example=0.5)],
        [_r("a", "fm^(2l+1)", "scattering length", "scalar"), _r("r", "fm^(1-2l)", "effective range", "scalar"),
         _r("a_spread", "fm^(2l+1)", "difference of a to two other fits", "scalar"),
         _r("r_spread", "fm^(1-2l)", "difference of r to two other fits", "scalar"),
         _r("e_window", "MeV", "the window (low, high), c.m.", "list"),
         _r("points", "", "energies fitted", "scalar"),
         _r("bound_states", "", "bound states of the potential in the wave", "scalar")],
        at_x=False),
    "potential_depth_for": _row(
        "potential_depth_for", "The depth (and Gaussian width) of a pair's potential in one partial wave for a chosen "
        "scattering length (and effective range).",
        "On the branch of its present number of bound states: the depth for a (one per branch), or the depth and width for "
        "a and r by nested searches over widths 0.4-2.5 times the present one; not found, with the reason, when there "
        "is none (docs/upgrade/potential_background.md sections 13.2, 13.4).",
        [dict(_ARG, name="pair", kind="pair", doc="particle pair key of <levels> (a pair with a potential)",
              domain="a particle pair of the project"),
         dict(_ARG, name="a", kind="real", unit="fm^(2l+1)", doc="the scattering length asked", domain="finite",
              example=40.0),
         dict(_ARG, name="l", kind="integer", doc="orbital angular momentum (0 for a charged pair)",
              domain="integer >= 0", lower=0.0, default=0, example=0),
         dict(_ARG, name="J", kind="real", doc="the channel's total J: required when the potential depends on J, "
              "ignored when it does not", domain="a multiple of 1/2 >= 0, or None", example=0.5),
         dict(_ARG, name="r", kind="real", unit="fm^(1-2l)", doc="the effective range asked; None: the depth alone",
              domain="finite, or None", example=0.93)],
        [_r("found", "", "a solution on the branch (bool)", "scalar"),
         _r("why", "", "the reason when none", "text"),
         _r("depth", "MeV", "V0 (None when not found)", "scalar"),
         _r("width", "fm", "the Gaussian r0 (None when not found)", "scalar"),
         _r("a", "fm^(2l+1)", "a at that depth and width", "scalar"),
         _r("r", "fm^(1-2l)", "r at that depth and width", "scalar"),
         _r("bound_states", "", "the branch: bound states of the wave", "scalar"),
         _r("nodes", "", "nodes of the zero-energy solution inside r_m", "scalar"),
         _r("r_range", "fm^(1-2l)", "the r found with a over the widths scanned (low, high)", "list")],
        at_x=False),
    "jpi_block": _row(
        "jpi_block", "One J^pi group on excitation energies: channel functions, U, T, phase shifts, |A| and det A^-1.",
        "The level matrix A, the collision matrix U and T of the group's open particle channels (physics/"
        "CollisionMatrix.h), the nuclear phase shift as the engine's phase_shift data kind forms it, |A_lm| per pair "
        "of levels and det A^-1, whose modulus is least near each resonance (docs/upgrade/physics_tab.md sections 4 "
        "and 5).  The energies are E_x of the compound nucleus; each channel's E_c = E_x - its threshold.  A photon "
        "channel enters A only with R-matrix capture and has no U or T element here.",
        [dict(_ARG, name="jgroup", kind="integer", doc="the J group, 1-based as parameters.out numbers them",
              domain="integer >= 1", lower=1.0, example=1),
         dict(_ARG, name="e_x", kind="real", unit="MeV", doc="excitation energy E_x of the compound nucleus",
              domain="finite", example=2.0, sweep={"lower": 0.5, "upper": 10.0, "spacing": "linear"}),
         _PARAMS],
        [_r("E_x", "MeV", "the energies"), _r("jgroup", "", "the J group", "scalar"),
         _r("J", "", "its spin", "scalar"), _r("pi", "", "its parity", "scalar"),
         _r("parametrization", "", "standard, brune or park", "text"),
         _r("levels[].level", "", "a level in the R-matrix (1-based in the group)", "scalar", _BLOCK_LEVELS),
         _r("levels[].E_level", "MeV", "its energy E_lambda at params (formal in the standard basis)", "scalar",
            _BLOCK_LEVELS),
         _r("channels[].channel", "", "the channel (1-based)", "scalar", _BLOCK_CHANNELS),
         _r("channels[].pair", "", "its pair key", "scalar", _BLOCK_CHANNELS),
         _r("channels[].radiation", "", "P particle, E or M photon", "text", _BLOCK_CHANNELS),
         _r("channels[].l", "", "orbital angular momentum (a photon's multipolarity)", "scalar", _BLOCK_CHANNELS),
         _r("channels[].s", "", "channel spin", "scalar", _BLOCK_CHANNELS),
         _r("channels[].radius", "fm", "channel radius a_c (0 for a photon)", "scalar", _BLOCK_CHANNELS),
         _r("channels[].threshold", "MeV", "E_x of the pair's threshold", "scalar", _BLOCK_CHANNELS),
         _r("channels[].boundary", "", "boundary value B_c as the level matrix reads it", "scalar", _BLOCK_CHANNELS),
         _r("channels[].E_c", "MeV", "c.m. energy of the channel's pair", "axis", _BLOCK_CHANNELS),
         _r("channels[].P", "", "penetrability (0 closed)", "axis", _BLOCK_CHANNELS),
         _r("channels[].S", "", "shift function", "axis", _BLOCK_CHANNELS),
         _r("channels[].L_re", "", "S - B_c, the real part of L_c - B_c", "axis", _BLOCK_CHANNELS),
         _r("channels[].L_im", "", "P, the imaginary part of L_c - B_c", "axis", _BLOCK_CHANNELS),
         _r("elements[].channel", "", "entrance channel c", "scalar", _BLOCK_ELEMENTS),
         _r("elements[].channel_out", "", "exit channel c'", "scalar", _BLOCK_ELEMENTS),
         _r("elements[].U_re", "", "Re U_cc'", "axis", _BLOCK_ELEMENTS),
         _r("elements[].U_im", "", "Im U_cc'", "axis", _BLOCK_ELEMENTS),
         _r("elements[].U_abs", "", "|U_cc'|", "axis", _BLOCK_ELEMENTS),
         _r("elements[].T_re", "", "Re T_cc'", "axis", _BLOCK_ELEMENTS),
         _r("elements[].T_im", "", "Im T_cc'", "axis", _BLOCK_ELEMENTS),
         _r("phase_shifts[].channel", "", "the channel", "scalar", _BLOCK_PHASES),
         _r("phase_shifts[].pair", "", "its pair key", "scalar", _BLOCK_PHASES),
         _r("phase_shifts[].l", "", "orbital angular momentum", "scalar", _BLOCK_PHASES),
         _r("phase_shifts[].s", "", "channel spin", "scalar", _BLOCK_PHASES),
         _r("phase_shifts[].delta", "deg", "nuclear phase shift arg(U)/2", "axis", _BLOCK_PHASES),
         _r("phase_shifts[].eta", "", "inelasticity |U_cc|", "axis", _BLOCK_PHASES),
         _r("level_pairs[].level", "", "level l", "scalar", _BLOCK_LEVEL_PAIRS),
         _r("level_pairs[].level_out", "", "level m", "scalar", _BLOCK_LEVEL_PAIRS),
         _r("level_pairs[].A_abs", "MeV^-1", "|A_lm|", "axis", _BLOCK_LEVEL_PAIRS),
         _r("det_re", "MeV^n", "Re det A^-1"), _r("det_im", "MeV^n", "Im det A^-1"),
         _r("log10_abs_det", "", "log10 |det A^-1|, least near a resonance"),
         _r("singular", "", "A^-1 singular at E_x (bool)")],
        axis="e_x"),
    "thm_experiments": _row(
        "thm_experiments", "Every THM experiment as a chi-squared evaluation at params profiles it.",
        "Per experiment (the answer's keys are their names): its data sets, the profiled norm n* and background b(E) = "
        "b0 + b1 E + b2 E^2 with their covariance from the profile, its chi-squared and points.  Runs a chi-squared "
        "evaluation.", [_PARAMS], _profile("{experiment}.", _EXPERIMENTS)),
    "thm_background": _row(
        "thm_background", "The shared norm and background of one THM experiment at params.",
        "The profiled norm n* and background coefficients with their covariance from the profile (at fixed R-matrix "
        "parameters, not scaled by chi2/nu); E in b(E) is the c.m. energy of the THM entrance pair.",
        [_EXPERIMENT, _PARAMS], _profile("", None)),
    "thm_lineshape": _row(
        "thm_lineshape", "The Coulomb line shape |N_C|^2 of an experiment with lineshape=on.",
        "|N_C|^2 = exp[2 zeta arctan(2 (E_l - E)/G_l)] per level of a J group coupling entrance and exit "
        "(Mukhamedzhanov, Kadyrov & Pang, EPJA 56 (2020) 233), zeta = eta_sB - eta_0 per exit pair, and eta_sb, the "
        "neglected s-b parameter (a validity check: it should be << 1).  The level poles are those of params.",
        [_EXPERIMENT, _ENERGIES, _PARAMS],
        [_r("experiment", "", "the experiment's name", "text"), _r("E", "MeV", "the energies"),
         _r("E_sF", "MeV", "E_aA - B - E"), _r("eta_0", "", "Z_s Z_F alpha mu_sF/k_sF"),
         _r("E_aA", "MeV", "c.m. energy of a + A", "scalar"), _r("B", "MeV", "binding energy of a = x + s", "scalar"),
         _r("Zs", "", "spectator charge", "scalar"), _r("ZF", "", "charge of F", "scalar"),
         _r("spectator", "", "the spectator", "text"),
         _r("exits[].pair", "", "the exit pair's key", "scalar", _EXITS),
         _r("exits[].Zb", "", "charge of the light exit nucleus b", "scalar", _EXITS),
         _r("exits[].ZB", "", "charge of the heavy exit nucleus B", "scalar", _EXITS),
         _r("exits[].mb", "u", "mass of b", "scalar", _EXITS), _r("exits[].mB", "u", "mass of B", "scalar", _EXITS),
         _r("exits[].zeta", "", "eta_sB - eta_0", "axis", _EXITS),
         _r("exits[].eta_sb", "", "the neglected s-b Sommerfeld parameter, averaged over the b direction", "axis",
            _EXITS),
         _r("exits[].levels[].jgroup", "", "the level's J group (1-based)", "scalar", _EXIT_LEVELS),
         _r("exits[].levels[].level", "", "the level (1-based, as parameters.out)", "scalar", _EXIT_LEVELS),
         _r("exits[].levels[].J", "", "its spin", "scalar", _EXIT_LEVELS),
         _r("exits[].levels[].pi", "", "its parity", "scalar", _EXIT_LEVELS),
         _r("exits[].levels[].E_level", "MeV", "its pole energy, c.m. of x + A", "scalar", _EXIT_LEVELS),
         _r("exits[].levels[].Gamma", "MeV", "its total width", "scalar", _EXIT_LEVELS),
         _r("exits[].levels[].NC2", "", "|N_C|^2", "axis", _EXIT_LEVELS)],
        axis="energies"),
    "thm_vertex": _row(
        "thm_vertex", "The THM entrance vertex |M_l|^2 per channel and level, averaged over the spectator window.",
        "<|M_l|^2> over the experiment's spectator-momentum window and |M_l|^2 at p_s = 0 per entrance channel and "
        "level (the vertex boundary of each level from params); with vertexModel=dw the surface term of the prior-form "
        "DWBA and its plane-wave limit.  Where reached is False the engine uses the nearest data point's vertex.",
        [_EXPERIMENT, _ENERGIES, _PARAMS,
         dict(_ARG, name="strict", kind="flag", doc="refuse (ValueError) energies out of reach instead of a warning",
              domain="True or False", default=False)],
        [_r("experiment", "", "the experiment's name", "text"), _r("E", "MeV", "the energies"),
         _r("window", "", "the ps window (delta without one)", "text"),
         _r("model", "", "pw or dw", "text"), _r("pair", "", "the entrance pair's key", "scalar"),
         _r("B", "MeV", "binding energy used", "scalar"), _r("mu_sx", "MeV", "reduced mass of s + x", "scalar"),
         _r("radius", "fm", "channel radius", "scalar"),
         _r("reached", "", "False where the window (dw: the vertex grid) does not reach E (bool)"),
         _r("p_s", "MeV/c", "the window's nodes |p_s| at each E", "axis x nodes"),
         _r("weights", "", "their normalized weights", "axis x nodes"),
         _r("T_s", "MeV", "p_s^2/2mu_sx at each node", "axis x nodes"),
         _r("theta_cm", "deg", "the nodes' spectator c.m. angles to the beam", "axis x nodes"),
         _r("rho", "", "rho = p_xA a at each node (dw: p a)", "axis x nodes"),
         _r("dw_q", "MeV/c", "dw: the nodes q at each E", "axis x nodes"),
         _r("dw_weights", "", "dw: their weights", "axis x nodes"),
         _r("dw_q_delta", "MeV/c", "dw: q of the spectatorAngle node"),
         _r("dw_p_delta", "fm^-1", "dw: p of the spectatorAngle node"),
         _r("angle_theta_cm", "deg", "spectatorAngles=: the accepted directions (c.m., to the beam)", "axis x nodes"),
         _r("angle_q", "MeV/c", "spectatorAngles=: their q", "axis x nodes"),
         _r("angle_weights", "", "spectatorAngles=: their weights", "axis x nodes"),
         _r("channels[].jgroup", "", "the J group (1-based)", "scalar", _CHANNELS),
         _r("channels[].channel", "", "the channel (1-based)", "scalar", _CHANNELS),
         _r("channels[].J", "", "J of the group", "scalar", _CHANNELS),
         _r("channels[].pi", "", "its parity", "scalar", _CHANNELS),
         _r("channels[].l", "", "orbital angular momentum", "scalar", _CHANNELS),
         _r("channels[].s", "", "channel spin", "scalar", _CHANNELS),
         _r("channels[].levels[].level", "", "the level (1-based)", "scalar", _CHANNEL_LEVELS),
         _r("channels[].levels[].boundary", "", "the vertex boundary B_c of the level (NaN for vertex=onshell)",
            "scalar", _CHANNEL_LEVELS),
         _r("channels[].levels[].M2", "", "<|M_l|^2> over the window", "axis", _CHANNEL_LEVELS),
         _r("channels[].levels[].M2_qf", "", "|M_l|^2 at p_s = 0 (dw: at the spectatorAngle direction)", "axis",
            _CHANNEL_LEVELS),
         _r("channels[].levels[].M2_pw", "", "dw: the plane-wave |M_l(p)|^2 at p = |k_aA - alpha k_sF|", "axis",
            _CHANNEL_LEVELS)],
        axis="energies"),
    "thm_distortion": _row(
        "thm_distortion", "The distortion factor R(E) of an experiment with distortion=.",
        "Zero-range prior-form DWBA: R(E) = rho(E)/rho(E_ref), rho = |M|^2/|M_PW|^2 (dwpw) or |M|^2 (dw); R_model is "
        "the R the model uses (interpolated; the table's for distortion=table).  R does not depend on the R-matrix "
        "parameters.",
        [_EXPERIMENT, _ENERGIES],
        [_r("experiment", "", "the experiment's name", "text"), _r("E", "MeV", "the energies"),
         _r("R_model", "", "the R the model uses"), _r("R", "", "R evaluated at E"),
         _r("M2", "fm^6", "|M|^2 (up to a constant)"), _r("M2_PW", "fm^6", "|M_PW|^2 (up to a constant)"),
         _r("E_sF", "MeV", "c.m. energy of s + F"), _r("k_sF", "fm^-1", "its wave number"),
         _r("eta_sF", "", "its Sommerfeld parameter"), _r("theta_cm", "deg", "spectator c.m. angle to the beam"),
         _r("x", "", "cos(k_sF, k_aA)"), _r("q", "fm^-1", "|k_sF - beta k_aA|"),
         _r("lmax", "", "the highest partial wave summed"),
         _r("kind", "", "coulomb, optical or table", "text"), _r("description", "", "the settings", "text"),
         _r("E_ref", "MeV", "reference energy", "scalar"), _r("E_aA", "MeV", "c.m. energy of a + A", "scalar"),
         _r("B", "MeV", "binding energy of a = x + s", "scalar"), _r("k_aA", "fm^-1", "a + A wave number", "scalar"),
         _r("eta_aA", "", "a + A Sommerfeld parameter", "scalar"),
         _r("kappa", "fm^-1", "bound-state wave number", "scalar"),
         _r("eta_b", "", "bound-state Sommerfeld parameter", "scalar"), _r("beta", "", "m_s/m_a", "scalar")],
        axis="energies", at_x=False),
}


def session_functions():
    """The rows of :data:`SESSION_FUNCTIONS`, a fresh copy (plain data)."""
    return copy.deepcopy(SESSION_FUNCTIONS)
