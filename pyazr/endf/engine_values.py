"""The engine's values of the reactions a written file holds, on the writer's grids (ENDF-102, 2023 edition, cited
by printed page).

Every number is the engine's: a sweep is one session of a project made of the session's nuclear part at the fit
vector x, with one extrapolation per quantity at the sweep's lab energies (and, for elastic scattering, c.m.
angles).  The quantities:

- ``sigma`` of every exit, of capture to every final state and of capture to each final state (b);
- ``legendre`` with the expansion ``rsl``, A_L / ((2L+1) A_0), of every exit and every capture final state: the
  coefficients a_l of LAW=2 (eq. 6.8, p. 143: p(mu) = 1/2 + sum_l (2l+1)/2 a_l P_l(mu)) as the engine gives them;
- ``dsigma_domega`` of elastic scattering at c.m. angles (b/sr), from which the module subtracts the Coulomb term of
  LAW=5 as the file's reader adds it back (eqs. 6.9-6.10, the engine's rows on the file's masses, p. 145).

The energy grid of every section is one (:func:`energy_grid`): the rule of :mod:`pyazr.endf.grid` on each sigma
(MF=3), each yield sigma_f / sigma (MF=6/MT=102), each LAW=2 distribution p(mu) at the cosines ``LAW2_COSINES``
and the elastic cross section sigma_e at ``ELASTIC_COSINES`` cosines of its range (the observable; the Coulomb term
the reader adds is exact at every energy).  The cosine grid of LAW=5 (:func:`elastic_cosines`) is made by the same
rule on sigma_e - sigma_c at every energy of the grid, held against the size of the cross section.
"""

import dataclasses
import os

import numpy as np
from numpy.polynomial import legendre as legendre_poly

from ..errors import FormatRefused
from . import grid

#: The cosines (c.m.) at which a LAW=2 distribution enters the rule of the energy grid.
LAW2_COSINES = np.linspace(-1.0, 1.0, 9)
#: How many cosines of the LAW=5 range enter the rule of the energy grid (evenly spaced, both ends included).
ELASTIC_COSINES = 9
#: Evenly spaced cosines the LAW=5 cosine grid starts from.
ELASTIC_START_COSINES = 21
#: The session options a sweep keeps (the calculation's; the others are the session's own business).
SWEEP_OPTIONS = ("use_brune", "use_park", "use_rmc", "use_gsl_coul", "use_coulomb_table", "ignore_externals",
                 "use_long_wavelength", "ec_convection_current")


@dataclasses.dataclass(frozen=True)
class Reaction:
    """One reaction of the written file: ``mt``; ``role`` -- ``"elastic"`` (MT=2), ``"exit"`` (a particle exit),
    ``"capture"`` (MT=102, every final state), ``"final_state"`` (MT=102, one photon pair); the project's reaction
    ``label`` and the ``pair_keys`` an extrapolation of it names (field, key); the session's ``exit_key`` (-1 for
    every final state); ``excitation`` of the final state (MeV, for a final state)."""
    mt: int
    role: str
    label: str
    pair_keys: tuple
    exit_key: int
    excitation: float = 0.0

    @property
    def key(self):
        """A name of the reaction's columns: the MT, and the pair key for a final state."""
        return (self.mt, self.exit_key) if self.role == "final_state" else (self.mt,)

    @property
    def has_legendre(self):
        return self.role in ("exit", "final_state")


@dataclasses.dataclass
class EnergyGrid:
    """The energy grid (lab MeV, increasing) and the engine's values on it: ``sigma`` {key: b}, ``legendre`` {key:
    (n_E, L+1) array of a_l, a_0 = 1}, ``elastic`` (n_E, n_mu) sigma_e - sigma_c at ``elastic_mu`` (b/sr) or None,
    and the refinement's report."""
    energies: np.ndarray
    sigma: dict
    legendre: dict
    elastic: object
    elastic_mu: object
    report: dict


def law2_probability(coefficients, mu):
    """p(mu) of LAW=2 from a_0 = 1, a_1, ... (eq. 6.8, p. 143): sum_l (2l+1)/2 a_l P_l(mu)."""
    c = np.asarray(coefficients, float)
    return legendre_poly.legval(np.asarray(mu, float), c * (2 * np.arange(len(c)) + 1) / 2.0)


class Sweeps:
    """The sessions of the writer: a project of ``session``'s nuclear part at the fit vector ``x`` in ``workdir``,
    evaluated once per call.  ``legendre_order``: the highest order the engine is asked for (the last one must come
    back zero: :meth:`energy_grid` refuses otherwise)."""

    def __init__(self, session, x, workdir, legendre_order):
        self.session, self.workdir, self.order = session, workdir, int(legendre_order)
        for d in ("output", "checks"):
            os.makedirs(os.path.join(workdir, d), exist_ok=True)
        self.parameter_file = self._parameter_file(x)
        self.options = {k: v for k, v in session.options.items() if k in SWEEP_OPTIONS}
        self.count = 0

    def _parameter_file(self, x):
        """Every parameter of the session at x, full precision, in the session's basis (read by name)."""
        P = self.session.parameters
        x = P.x0 if x is None else np.asarray(x, float)
        path = os.path.join(self.workdir, "x.par")
        with open(path, "w") as fh:
            fh.write(f"#parametrization {self.session.basis}\n")
            for p in P:
                fh.write(f"{p.name} {float(p.start if p.fixed else x[p.free_index])!r}\n")
        return path

    def evaluate(self, records):
        """The engine's ``EvaluatedData`` of a project with the extrapolation ``records`` (dicts of
        ``Project.add_extrapolation``'s fields with a ``name``), in their order."""
        from ..project import Project
        from ..session import Session
        self.count += 1
        path = os.path.join(self.workdir, "sweep.azr")
        if os.path.exists(path):
            os.remove(path)
        p = Project.new(path, nuclear_part=self.session.file)
        for r in records:
            r = dict(r)
            p.add_extrapolation(r.pop("name"), **r)
        p.save()
        s = Session(path, extrapolations=True, parameter_file=self.parameter_file, **self.options)
        try:
            return s.evaluate(s.parameters.x0).data
        finally:
            s.close()

    def _energies(self, name, energies_MeV):
        """A file of lab energies (MeV), absolute: the project is read without its directory too."""
        path = os.path.join(self.workdir, f"{name}.dat")
        np.savetxt(path, np.asarray(energies_MeV, float), fmt="%.17g")
        return f"file {path} MeV lab"

    @staticmethod
    def _angles(mu):
        return "list " + " ".join(repr(float(np.degrees(np.arccos(np.clip(m, -1.0, 1.0))))) for m in mu) + " deg cm"

    def _sigma_and_legendre(self, reactions, energies):
        """The records of sigma (every reaction but elastic) and Legendre coefficients (exits and final states)."""
        records, roles = [], []
        for r in reactions:
            if r.role == "elastic":
                continue
            records.append(dict(name=f"s{len(records)}", reaction=r.label, quantity="sigma", unit="b",
                                energies=energies, **dict(r.pair_keys)))
            roles.append(("sigma", r))
            if r.has_legendre:
                records.append(dict(name=f"s{len(records)}", reaction=r.label, quantity="legendre", unit="1",
                                    expansion="rsl", orders=f"0 {self.order}", energies=energies,
                                    **dict(r.pair_keys)))
                roles.append(("legendre", r))
        return records, roles

    def _elastic_record(self, reaction, energies, mu):
        return dict(name="elastic", reaction=reaction.label, quantity="dsigma_domega", unit="b/sr",
                    energies=energies, angles=self._angles(mu), **dict(reaction.pair_keys))

    @staticmethod
    def _elastic_block(data, n_energies, mu, coulomb_minus, energies_MeV):
        """sigma_e - sigma_c (n_E, n_mu) from an elastic extrapolation evaluated angle by angle (the engine lists
        every energy of an angle together); sigma_e itself when ``coulomb_minus`` is None."""
        model = np.asarray(data.model, float).reshape(len(mu), n_energies)
        angles = np.asarray(data.angles, float).reshape(len(mu), n_energies)
        expected = np.degrees(np.arccos(np.clip(mu, -1.0, 1.0)))
        if not np.allclose(angles[:, 0], expected, atol=1e-9):
            raise FormatRefused(what="the engine's elastic points do not come angle by angle",
                                allowed="the order of an extrapolation with an angle list",
                                fix="report this to the AZURE2 developers")
        return model.T if coulomb_minus is None else model.T - coulomb_minus(energies_MeV, mu)

    def energy_grid(self, reactions, start, tolerance, max_points, elastic_range=None):
        """The one energy grid of every section (the module's text), refined from ``start`` (lab MeV).
        ``elastic_range``: (mu_min, mu_max) of LAW=5, or None when elastic scattering is not written."""
        elastic = next((r for r in reactions if r.role == "elastic"), None) if elastic_range else None
        mu_el = np.linspace(*elastic_range, ELASTIC_COSINES) if elastic is not None else None
        capture = next((r for r in reactions if r.role == "capture"), None)
        cache = {}

        def evaluate_sorted(xs):
            """The engine's values at xs (lab MeV), computed in increasing order and returned in the order of xs."""
            order = np.argsort(xs)
            e = np.asarray(xs, float)[order]
            energies = self._energies("grid", e)
            records, roles = self._sigma_and_legendre(reactions, energies)
            if elastic is not None:
                records.append(self._elastic_record(elastic, energies, mu_el))
            data = self.evaluate(records)
            back = np.empty_like(order)
            back[order] = np.arange(len(order))
            out = {"sigma": {}, "legendre": {}}
            for (what, r), d in zip(roles, data):
                if what == "sigma":
                    out["sigma"][r.key] = np.asarray(d.model, float)[back]
                else:
                    a = np.zeros((len(e), self.order + 1))
                    for k, v in enumerate(d.extra.get("legendre", [])):
                        v = np.nan_to_num(np.asarray(v, float))[:self.order + 1]
                        a[k, :len(v)] = v
                    a[:, 0] = 1.0  # A_0 / A_0; a point without a distribution (sigma = 0) is isotropic
                    out["legendre"][r.key] = a[back]
            if elastic is not None:
                out["elastic"] = self._elastic_block(data[-1], len(e), mu_el, None, e)[back]
            return out

        def model(xs):
            out = evaluate_sorted(xs)
            columns = []  # each (n_E,) or (n_E, n_mu)
            for r in reactions:
                if r.role in ("exit", "capture"):
                    columns.append(out["sigma"][r.key])
                elif r.role == "final_state":
                    total = out["sigma"][capture.key]
                    columns.append(np.divide(out["sigma"][r.key], total, out=np.zeros(len(xs)), where=total != 0))
                if r.has_legendre:
                    a = out["legendre"][r.key]
                    columns.append(np.array([law2_probability(a[k], LAW2_COSINES) for k in range(len(xs))]))
            if elastic is not None:
                columns.append(out["elastic"])
            for k, x in enumerate(xs):
                cache[x] = {"sigma": {key: v[k] for key, v in out["sigma"].items()},
                            "legendre": {key: v[k] for key, v in out["legendre"].items()},
                            "elastic": out["elastic"][k] if elastic is not None else None}
            return [np.concatenate([np.atleast_1d(c[k]) for c in columns]) for k in range(len(xs))]

        refined = grid.refine(start, model, tolerance, max_points)
        e = np.asarray(refined["x"], float)
        sigma = {r.key: np.array([cache[x]["sigma"][r.key] for x in e]) for r in reactions if r.role != "elastic"}
        legendre = {r.key: np.array([cache[x]["legendre"][r.key] for x in e]) for r in reactions if r.has_legendre}
        for key, a in legendre.items():
            if a.size and np.any(np.abs(a[:, -1]) > 1e-12):
                raise FormatRefused(what=f"the engine's Legendre coefficients of MT={key[0]} reach order "
                                    f"{self.order}, the highest asked", allowed="coefficients that end below it",
                                    fix="report the project to the AZURE2 developers")
        el = np.array([cache[x]["elastic"] for x in e]) if elastic is not None else None
        return EnergyGrid(e, sigma, legendre, el, mu_el, refined)

    def elastic_cosines(self, reaction, energies_MeV, mu_range, coulomb_minus, tolerance, max_points):
        """The cosine grid of LAW=5, from ``ELASTIC_START_COSINES`` even cosines over ``mu_range``, by the rule of
        :mod:`pyazr.endf.grid` with one component per energy of the grid: sigma_e - sigma_c (what the file tabulates)
        plus a constant C(E), the smallest sigma_e of the start cosines at that energy.  The constant cancels in the
        rule's chord, so the rule holds the interpolation error of sigma_e - sigma_c below tolerance x
        max(|sigma_e - sigma_c + C|, ...): about tolerance x sigma_e where the nuclear part is small against the
        cross section, and relative to the nuclear part where it is large.  (The rule on sigma_e - sigma_c alone asks
        for thousands of cosines where it is a small, oscillating part of a cross section; the file's error is in the
        cross section.)  The grid serves every energy; each energy keeps the cosines its own values need
        (:func:`pyazr.endf.grid.thin`, the same rule).  Returns ([(cosines, sigma_e - sigma_c in b/sr) per energy],
        the refinement's report)."""
        energies = self._energies("elastic", energies_MeV)
        e = np.asarray(energies_MeV, float)
        offset = {}

        def model(mus):
            order = np.argsort(mus)
            m = np.asarray(mus, float)[order]
            data = self.evaluate([self._elastic_record(reaction, energies, m)])
            full = self._elastic_block(data[0], len(e), m, None, e)
            block = full - coulomb_minus(e, m)
            if "C" not in offset:
                offset["C"] = np.min(np.abs(full), axis=1)
            back = np.empty_like(order)
            back[order] = np.arange(len(order))
            return [block[:, k] + offset["C"] for k in back]

        refined = grid.refine(np.linspace(*mu_range, ELASTIC_START_COSINES), model, tolerance, max_points)
        mu = np.asarray(refined["x"], float)
        shifted = np.asarray(refined["values"], float).reshape(len(mu), len(e)).T
        panels = []
        for row, c in zip(shifted, offset["C"]):  # each energy keeps the cosines it needs (grid.thin, the same rule)
            k = grid.thin(mu, row, tolerance)
            panels.append((mu[k], row[k] - c))
        return panels, refined
