"""The energies and panels of charged-particle elastic scattering as the writer writes it, LAW=5 LTP=12 with
sigma_NI in MF=3/MT=2 (ENDF-102, 2023 edition, pp. 144-147), so that the file holds its tolerance under the
manual's reading.

The manual's reading (p. 146: "The angular distribution p_NI is given in File 6 as a tabulated function of mu, and
sigma_NI(E) in barns is given in File 3"; eqs. 6.19-6.20): a processing code interpolates sigma_NI(E) in MF=3 and
P_NI(mu, E) in MF=6, each under its own law, and multiplies.  Between two tabulated energies where sigma_NI changes
sign, P_NI = (sigma_e - sigma_c)/sigma_NI is large and of opposite sign on the two sides, and the product of the two
interpolations misses sigma_e - sigma_c by about its own size, however close the two energies are.  So:

1. **Brackets.**  Every sign change of sigma_NI between two energies of the writer's grid is closed in by new
   energies (the engine's values there) until the two energies that hold it are at most ``BRACKET_WIDTH`` x E apart:
   a ladder of energies E_0 -+ (w/2) ``LADDER_RATIO``^k around the zero E_0 estimated from the chord, as many
   rounds as needed.  Inside a bracket the file's elastic cross section under the manual's reading is wrong by up to
   the nuclear-plus-interference part; the bracket's width is what the file states (MF=1, ``Written.elastic``).
2. **The energies kept**: from the grid and the ladders, by the rule of :mod:`pyazr.endf.grid`
   (:func:`pyazr.endf.grid.thin_to`) applied to sigma_e on the cosine grid, with the reader's own reconstruction
   between two kept energies (sigma_NI and P_NI each interpolated linearly, INT=2, then multiplied) instead of the
   chord; the two ends of every bracket are kept.
3. **The quarter points checked**: the engine computes the quarter points of every kept interval that holds no
   computed energy (brackets and intervals narrower than ``SMALLEST_INTERVAL`` x E aside), and the choice is made
   again, until every such point passes (the refinement's rule under the reader's reconstruction; quarter points,
   since near a zero of sigma_NI the reader's error is not largest at the midpoint) or the energies computed would
   exceed ``CANDIDATES_PER_POINT`` x ``max_points`` (``checked`` says which).  The value compared at a computed
   energy is sigma_e as its own panel gives it, so that the rule measures the interpolation between energies; each
   panel's cosines hold the rule by their own thinning.
4. **The tolerance stated**: the engine computes the quarter points and the midpoint of every written interval
   (brackets aside) once more, and the largest deviation of the manual's reading there, cosines' error included,
   is part of the tolerance the file states.
"""

import dataclasses

import numpy as np

from ..errors import FormatRefused
from . import grid

#: The largest width of the bracket of a sign change of sigma_NI, relative to its energy (lab): inside it the
#: manual's reading of the file is wrong by up to sigma_e - sigma_c.
BRACKET_WIDTH = 1e-5
#: The ratio of the ladder's steps away from an estimated zero of sigma_NI (the energies the thinning chooses from:
#: the reader's error between two energies whose sigma_NI differ by a ratio r is about (r-1)^2/(4r) of
#: sigma_e - sigma_c, 6e-4 for 1.05).
LADDER_RATIO = 1.05
#: The energies the search may compute, per point of the writer's ``max_points``.
CANDIDATES_PER_POINT = 2
#: An interval narrower than this (relative to its energy) is not split by the midpoint check.
SMALLEST_INTERVAL = 1e-7
#: Rounds of the search (each one engine session): ladders while a sign change is not closed in, then the
#: midpoints of the kept intervals that hold no computed energy, until every such midpoint passes the rule.
REFINE_ROUNDS = 16


@dataclasses.dataclass(frozen=True)
class ElasticSection:
    """The elastic section as written: ``energies`` (lab MeV), ``panels`` [(cosines, sigma_e - sigma_c in b/sr)] per
    energy, ``sigma_ni`` (b) per energy, ``brackets`` [(E_low, E_high) in lab MeV] of the sign changes of sigma_NI,
    ``reached`` (the largest deviation by the rule at an energy the engine computed and the file dropped, under the
    manual's reading), ``candidates`` (energies the engine computed), ``sessions`` (the search's) and ``checked``
    (True when every written interval passed the check of its quarter points; False when the search stopped at
    its budget or its rounds first), ``cosines_reached`` (the largest deviation by the rule at a cosine a panel
    dropped: above the tolerance only where the cosine limit cut a panel)."""
    energies: np.ndarray
    panels: list
    sigma_ni: np.ndarray
    brackets: list
    reached: float
    candidates: int
    sessions: int
    checked: bool = False
    cosines_reached: float = 0.0


def panel_sigma_ni(panel):
    """sigma_NI of a panel (eq. 6.19): the trapezoid integral over its cosines of sigma_e - sigma_c, exact for P_NI
    linear in mu (LTP=12)."""
    trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    mu, v = panel
    return float(trapezoid(v, mu))


def _ladder(e_lo, e_hi, s_lo, s_hi):
    """Energies around the chord's zero of sigma_NI between e_lo and e_hi (lab MeV), inside the interval."""
    zero = e_lo - s_lo * (e_hi - e_lo) / (s_hi - s_lo)
    half = 0.5 * BRACKET_WIDTH * zero
    out = []
    for side in (-1.0, 1.0):
        step = half
        while True:
            e = zero + side * step
            if not e_lo < e < e_hi:
                break
            out.append(e)
            step *= LADDER_RATIO
    if not out:
        out.append(0.5 * (e_lo + e_hi))
    return out


def section(sweeps, reaction, energies_MeV, mu_range, coulomb_minus, tolerance, max_points, limit=None,
            cosine_limit=None):
    """The elastic section from the writer's energy grid ``energies_MeV`` (the module's text); at most ``limit``
    energies (None: as many as the tolerance asks) and ``cosine_limit`` cosines per energy.  Returns an
    :class:`ElasticSection` and the cosine refinement's report (``budget_reached`` when it stopped short).  Refused
    (:class:`pyazr.FormatRefused`): more sign changes of sigma_NI than ``limit`` can bracket."""
    mu, block, mu_report = sweeps.elastic_cosines(reaction, energies_MeV, mu_range, coulomb_minus, tolerance,
                                                  max_points)
    st = {"e": np.asarray(energies_MeV, float), "v": block}
    st["panels"] = [_panel(mu, row, coulomb_minus, ek, tolerance, mu_range, cosine_limit)
                    for ek, row in zip(st["e"], st["v"])]
    sessions = 0

    def add(new):
        new = np.array(sorted(set(new) - set(st["e"].tolist())))
        v_new = sweeps.elastic_values(reaction, new, mu, coulomb_minus)
        p_new = [_panel(mu, row, coulomb_minus, ek, tolerance, mu_range, cosine_limit) for ek, row in zip(new, v_new)]
        e = np.concatenate([st["e"], new])
        order = np.argsort(e)
        st["e"], st["v"] = e[order], np.vstack([st["v"], v_new])[order]
        allp = st["panels"] + p_new
        st["panels"] = [allp[k] for k in order]

    converged = False
    for _ in range(REFINE_ROUNDS):
        e = st["e"]
        s = np.array([panel_sigma_ni(p[:2]) for p in st["panels"]])
        ends = np.nonzero(s[:-1] * s[1:] < 0)[0]
        ladders = [x for i in ends if e[i + 1] - e[i] > BRACKET_WIDTH * e[i]
                   for x in _ladder(e[i], e[i + 1], s[i], s[i + 1])]
        if ladders:  # a sign change not yet closed in
            sessions += 1
            add(ladders)
            continue
        kept, reached = _thin(st, s, ends, mu, coulomb_minus, tolerance, limit)
        if limit is not None and len(kept) >= limit and reached > tolerance:
            break  # the limit, not the tolerance, ends the choice: no midpoint can be kept
        bracket_starts = set(ends.tolist())
        checks = [e[a] + q * (e[b] - e[a]) for a, b in zip(kept[:-1], kept[1:])
                  if b == a + 1 and a not in bracket_starts and e[b] - e[a] > SMALLEST_INTERVAL * e[a]
                  for q in (0.25, 0.75)]
        if not checks:
            converged = True
            break
        if len(e) + len(checks) > CANDIDATES_PER_POINT * max_points:
            break
        sessions += 1
        add(checks)
    e = st["e"]
    s = np.array([panel_sigma_ni(p[:2]) for p in st["panels"]])
    ends = np.nonzero(s[:-1] * s[1:] < 0)[0]
    kept, reached = _thin(st, s, ends, mu, coulomb_minus, tolerance, limit)
    brackets = [(float(e[i]), float(e[i + 1])) for i in ends]
    if limit is not None and 2 * len(ends) + 2 > limit:
        raise FormatRefused(what=f"elastic scattering: sigma_NI changes sign {len(ends)} times; the brackets and the "
                            f"two ends need {2 * len(ends) + 2} energies", allowed=f"at most {limit} energies per "
                            "MF=6 section", fix="write a narrower energy range, or ask for the manual's limits")
    between = _between_written(sweeps, reaction, st, s, kept, ends, mu, coulomb_minus)
    sessions += 1
    reached = max(reached, between)
    return ElasticSection(e[kept], [st["panels"][k][:2] for k in kept], s[kept], brackets, reached, len(e), sessions,
                          converged, max(st["panels"][k][2] for k in kept)), mu_report


def _between_written(sweeps, reaction, st, s, kept, ends, mu, coulomb_minus):
    """The largest deviation by the rule of the manual's reading at the quarter points and midpoint of every written
    interval (brackets aside), against the engine there (one session): what the file gives between its energies,
    with the cosines' own error."""
    e = st["e"]
    starts = set(ends.tolist())
    pairs = [(a, b) for a, b in zip(kept[:-1], kept[1:]) if a not in starts]
    if not pairs:
        return 0.0
    q = np.array([0.25, 0.5, 0.75])
    points = np.concatenate([e[a] + q * (e[b] - e[a]) for a, b in pairs])
    truth = sweeps.elastic_values(reaction, points, mu, coulomb_minus) + coulomb_minus(points, mu)
    p_at = _p_at_mu(mu, st["panels"], s)
    floor = 1e-3 * np.max(np.abs(st["v"] + coulomb_minus(e, mu)), axis=0)
    pred = []
    for a, b in pairs:
        t = q[:, None]
        pred.append(s[a] + (s[b] - s[a]) * t)
        pred[-1] = pred[-1] * (p_at[a] + (p_at[b] - p_at[a]) * t)
    pred = np.vstack(pred) + coulomb_minus(points, mu)
    return float(np.max(grid._deviation(truth, pred, floor)))


def _p_at_mu(mu, panels, s):
    """P_NI of every panel at the cosines ``mu`` (linear in mu, zero outside the panel, as the reader takes it)."""
    return np.array([np.interp(mu, p[0], p[1] / sk if sk != 0.0 else np.zeros(len(p[1])), left=0.0, right=0.0)
                     for p, sk in zip(panels, s)])


def _thin(st, s, ends, mu, coulomb_minus, tolerance, limit):
    """The energies kept of the computed ones (:func:`pyazr.endf.grid.thin_to` with the manual's reading between
    two kept energies; the two ends of every bracket kept) and the tolerance reached."""
    e, v, panels = st["e"], st["v"], [p[:2] for p in st["panels"]]
    coulomb = coulomb_minus(e, mu)
    p_at_mu = _p_at_mu(mu, panels, s)
    # sigma_e at the cosine grid as each energy's own panel gives it (the cosines' error is the panel's, held to the
    # rule by its own thinning): the rule here measures the interpolation between energies alone
    full = coulomb + s[:, None] * p_at_mu

    def manual_reading(i, j, inner):
        """sigma_e at the energies ``inner`` from the panels i and j as the manual's reader gets it."""
        t = ((e[inner] - e[i]) / (e[j] - e[i]))[:, None]
        sigma_ni = s[i] + (s[j] - s[i]) * t
        p_ni = p_at_mu[i] + (p_at_mu[j] - p_at_mu[i]) * t
        return coulomb[inner] + sigma_ni * p_ni

    required = sorted(set(ends.tolist()) | set((ends + 1).tolist()))
    return grid.thin_to(e, full, tolerance, limit=limit, required=required, reconstruct=manual_reading)


def _panel(mu, row, coulomb_minus, energy, tolerance, mu_range, cosine_limit=None):
    """The cosines one energy keeps (:func:`pyazr.endf.grid.thin` on sigma_e - sigma_c plus a constant, the smallest
    |sigma_e| at the start cosines of :meth:`pyazr.endf.engine_values.Sweeps.elastic_cosines`, which are points of
    ``mu``; at most ``cosine_limit`` of them), sigma_e - sigma_c there and the tolerance reached."""
    from .engine_values import ELASTIC_START_COSINES
    start = np.linspace(*mu_range, ELASTIC_START_COSINES)
    full = row + coulomb_minus(np.array([energy]), mu)[0]
    offset = float(np.min(np.abs(np.interp(start, mu, full))))
    k, reached = grid.thin_to(mu, row + offset, tolerance, limit=cosine_limit)
    return mu[k], row[k], reached
