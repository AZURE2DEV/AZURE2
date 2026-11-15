"""The energies of charged-particle elastic scattering as the writer writes it in the nuclear amplitude expansion,
LAW=5 LTP=1 (ENDF-102, 2023 edition, pp. 145-147), with 1.0 in MF=3/MT=2 (p. 145: "Since an integrated cross
section is not defined for this representation, a value of 1.0 is used in File 3").

The coefficients b_l and a_l at an energy are the engine's (:meth:`pyazr.endf.engine_values.Sweeps.
nuclear_amplitudes`, the door ``nuclear_amplitudes``): from them eqs. 6.13 and 6.14 give the engine's elastic
cross section at every angle, to rounding.  Between two energies of the file a reader interpolates sigma_e - sigma_c
at fixed mu, the Coulomb part at the energy asked (:mod:`pyazr.endf.distributions`; the manual is silent on the
quantity interpolated for LTP=1, and this reading is the module's for every LTP but 12-15).  So the file's error
lies between its energies alone, and there is no product of two interpolations and no sign change to bracket:

1. **The candidates**: the coefficients at every energy of the writer's grid (the one grid of MF=3 and MF=6).
2. **The energies kept**: by the rule of :mod:`pyazr.endf.grid` (:func:`pyazr.endf.grid.thin_to`) on the cross
   section the coefficients give at ``RULE_COSINES`` c.m. cosines, with the reader's interpolation between two kept
   energies instead of the chord; at most ``limit`` of them.
3. **The quarter points checked**: the engine gives the coefficients at the quarter points of every kept interval
   that holds no computed energy, and the choice is made again, until every such point passes, the limit binds or
   the energies computed would exceed ``CANDIDATES_PER_POINT`` x ``max_points`` (``checked`` says which).
4. **The tolerance stated**: the quarter points and the midpoint of every written interval are computed once more;
   the largest deviation there is part of the tolerance the file states (``reached``).

The cross section the rule holds is the engine's: the coefficients with the Coulomb terms on the project's masses
(the engine's own eta and k).  The reader computes the Coulomb terms on the file's masses (AWI nuclear, AWR atomic,
p. 24), which differ from the project's by the electron masses; that difference is in the measured deviation of
the file (the writer's test), not in the rule.
"""

import dataclasses

import numpy as np

from .distributions import amplitude_expansion

#: The cosines (c.m.) at which the rule holds the elastic cross section, evenly spaced over the range of LAW=5
#: (from -1, or 0 for identical particles, to the writer's largest cosine).
RULE_COSINES = 17
#: The energies the search may compute, per point of the writer's ``max_points``.
CANDIDATES_PER_POINT = 2
#: An interval narrower than this (relative to its energy) is not split by the quarter-point check.
SMALLEST_INTERVAL = 1e-7
#: Rounds of the quarter-point check (each one engine session).
REFINE_ROUNDS = 12


@dataclasses.dataclass(frozen=True)
class AmplitudeSection:
    """The elastic section as written in LTP=1: ``energies`` (lab MeV), ``records`` (the LIST values per energy,
    b/sr), ``reached`` (the largest deviation by the rule at an energy computed and not written, and at the quarter
    points and midpoint of every written interval), ``candidates`` (energies the engine computed), ``sessions`` (the
    search's) and ``checked`` (True when every written interval passed its quarter points)."""
    energies: np.ndarray
    records: list
    reached: float
    candidates: int
    sessions: int
    checked: bool


def section(sweeps, reaction, energies_MeV, mu_range, coulomb_engine, identical, two_s, tolerance, max_points,
            limit=None):
    """The LTP=1 section from the writer's grid ``energies_MeV`` (the module's text): at most ``limit`` energies
    (None: as many as the tolerance asks).  ``coulomb_engine``: the Coulomb terms on the project's masses
    (:class:`pyazr.endf.curves.CoulombTerms`), the engine's eta and sigma_c; ``two_s`` twice the spin of identical
    particles (the Mott term, eq. 6.10)."""
    mu = np.linspace(mu_range[0], mu_range[1], RULE_COSINES)
    st = {"e": np.zeros(0), "c": []}

    def values_at(e_MeV, records):
        """sigma_c and sigma_e - sigma_c (n_E, n_mu) of the records at their energies."""
        e_eV = np.asarray(e_MeV, float) * 1e6
        eta = coulomb_engine.eta(e_eV)
        part = np.array([amplitude_expansion(r, _nl(r, identical), identical, h, mu) for r, h in zip(records, eta)])
        sc = np.array([coulomb_engine.sigma(np.full(len(mu), x), mu, identical=identical, two_s=two_s) for x in e_eV])
        return sc, part.reshape(len(e_eV), len(mu))

    def add(new):
        new = np.array(sorted(set(np.asarray(new, float).tolist()) - set(st["e"].tolist())))
        if not new.size:
            return
        recs = sweeps.nuclear_amplitudes(reaction, new)
        e = np.concatenate([st["e"], new])
        order = np.argsort(e)
        allc = st["c"] + recs
        st["e"], st["c"] = e[order], [allc[k] for k in order]
        sc, part = values_at(st["e"], st["c"])
        st["sc"], st["part"] = sc, part

    add(energies_MeV)
    sessions = 1
    converged = False
    kept, reached = _thin(st, tolerance, limit)
    for _ in range(REFINE_ROUNDS):
        if limit is not None and len(kept) >= limit and reached > tolerance:
            break  # the limit, not the tolerance, ends the choice
        e = st["e"]
        checks = [e[a] + q * (e[b] - e[a]) for a, b in zip(kept[:-1], kept[1:])
                  if b == a + 1 and e[b] - e[a] > SMALLEST_INTERVAL * e[a] for q in (0.25, 0.75)]
        if not checks:
            converged = True
            break
        if len(e) + len(checks) > CANDIDATES_PER_POINT * max_points:
            break
        sessions += 1
        add(checks)
        kept, reached = _thin(st, tolerance, limit)
    # the tolerance stated: the quarter points and the midpoint of every written interval, once more
    e = st["e"]
    q = np.array([0.25, 0.5, 0.75])
    points = np.concatenate([e[a] + q * (e[b] - e[a]) for a, b in zip(kept[:-1], kept[1:])])
    between = 0.0
    if points.size:
        sessions += 1
        recs = sweeps.nuclear_amplitudes(reaction, points)
        sc, part = values_at(points, recs)
        truth = sc + part
        pred = []
        for a, b in zip(kept[:-1], kept[1:]):
            pred.append(st["part"][a] + (st["part"][b] - st["part"][a]) * q[:, None])
        pred = np.vstack(pred) + sc
        floor = 1e-3 * np.max(np.abs(st["sc"] + st["part"]), axis=0)
        from .grid import _deviation
        between = float(np.max(_deviation(truth, pred, floor)))
    return AmplitudeSection(e[kept], [st["c"][k] for k in kept], max(reached, between), len(e), sessions,
                            converged)


def _nl(values, identical):
    """NL of a LIST record of LTP=1 from its length (p. 147: NW = 4 NL + 3, or 3 NL + 3 for identical particles)."""
    n = len(values)
    return (n - 3) // 3 if identical else (n - 3) // 4


def _thin(st, tolerance, limit):
    """The energies kept of the computed ones (:func:`pyazr.endf.grid.thin_to` with the reader's interpolation of
    sigma_e - sigma_c between two kept energies, the Coulomb part at the energy) and the tolerance reached."""
    from .grid import thin_to
    e, sc, part = st["e"], st["sc"], st["part"]

    def reader(i, j, inner):
        t = ((e[inner] - e[i]) / (e[j] - e[i]))[:, None]
        return sc[inner] + part[i] + (part[j] - part[i]) * t

    return thin_to(e, sc + part, tolerance, limit=limit, reconstruct=reader)
