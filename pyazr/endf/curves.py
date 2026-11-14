"""An evaluation's values at a data set's points, in the data set's quantity, frame and units (``Material.curves``).

Which section (ENDF-102, 2023): elastic scattering (the data set's exit pair is its entrance pair) is MF=6/MT=2
LAW=5 (pp. 144-147); a particle exit pair is the discrete-level reaction MT = first + k of the emitted particle's
family (Appendix B, pp. 332, 339-341), k found by the file's own Q values, QM - QI = the residual's excitation
(p. 105); capture to every final state is MF=3/MT=102.  The sums of the format are the module's
(:mod:`pyazr.endf.distributions`); every conversion of physics is the engine's, as rows of :mod:`pyazr.physics`:

- lab energy of the file's projectile: the engine's lab energy of the point (``EvaluatedData.extra["energy_lab"]``)
  when the projectile is the pair's first particle, else the row ``cm_to_lab_energy`` with the projectile as beam;
- the Coulomb terms of LAW=5 (eqs. 6.9-6.12): ``lab_to_cm_energy``, ``sommerfeld_eta``, ``rutherford_cross_section``
  and ``mott_cross_section`` on the file's masses (AWI, AWR x the neutron mass of p. 394), non-relativistic as the
  manual writes them;
- ``ratio_to_rutherford`` and ``ratio_to_mott``: the rows ``rutherford_cross_section`` and ``mott_cross_section``
  on the data set's pair; ``sfactor``: the row ``sfactor_conversion``.

Angles stay c.m. cosines (pp. 143-144); the data set's angles are already c.m. (the engine's ``angles``).
A quantity the file cannot give (an analyzing power, a phase shift, capture to one final state, ...) is refused
with its three parts, never drawn as zero.
"""

import dataclasses
import math

import numpy as np

from ..errors import FormatRefused
from . import names
from .material import M_NEUTRON_U

QUANTITIES = ("sigma", "dsigma_domega", "ratio_to_rutherford", "ratio_to_mott", "sfactor")


@dataclasses.dataclass(frozen=True)
class PairSpec:
    """One particle pair as a curve needs it: charges, masses (u), mass numbers, the second particle's excitation
    (MeV) and spins; ``photon`` for a capture pair.  ``PairSpec.of(pair)`` takes a session's pair record."""
    z1: int
    z2: int
    m1_u: float
    m2_u: float
    a1: int
    a2: int
    excitation: float = 0.0
    light_excitation: float = 0.0
    j1: float = 0.0
    photon: bool = False

    @classmethod
    def of(cls, pair):
        return cls(int(pair.Z1), int(pair.Z2), float(pair.M1), float(pair.M2), int(round(pair.M1)),
                   int(round(pair.M2)), float(pair.heavy_excitation), float(pair.light_excitation or 0.0),
                   float(pair.J1), bool(pair.is_photon))

    @property
    def za1(self):
        return 1000 * self.z1 + self.a1

    @property
    def za2(self):
        return 1000 * self.z2 + self.a2


@dataclasses.dataclass(frozen=True)
class CurveRequest:
    """What a curve is asked for: the quantity (a data kind of the engine), the points (c.m. energies in MeV and
    c.m. angles in deg, as ``EvaluatedData``), the entrance pair and the exit pair (``None``: every capture pair),
    and optionally the engine's lab energies of the beam (the first particle of the entrance pair) per point."""
    quantity: str
    energies: tuple
    angles: tuple
    entrance: PairSpec
    exit: object
    energies_lab: tuple = None
    name: str = ""
    relativistic: bool = False


@dataclasses.dataclass
class Curve:
    """The evaluation at a data set's points: ``energies`` (MeV c.m.) and ``angles`` (deg c.m.) as the data set,
    ``values`` in ``unit`` (NaN where the file has no value), ``coverage`` (True where it has), ``source`` (the
    sections used), ``notes`` (what to know)."""
    name: str
    quantity: str
    unit: str
    frame: str
    energies: np.ndarray
    angles: np.ndarray
    values: np.ndarray
    coverage: np.ndarray
    source: str
    notes: list

    def __str__(self):
        lines = [f"{self.name or 'curve'}: {self.quantity} ({self.unit}, {self.frame}) from {self.source}; "
                 f"{int(self.coverage.sum())} of {len(self.values)} points covered"]
        lines += [f"  note: {n}" for n in self.notes]
        lines.append(f"{'E_cm/MeV':>14} {'angle/deg':>10} {'value':>14}")
        lines += [f"{e:14.7g} {a:10.5g} {v:14.7g}" for e, a, v in zip(self.energies, self.angles, self.values)]
        return "\n".join(lines)


class CoulombTerms:
    """eta(E) and the Coulomb cross section sigma_c(E, mu) of LAW=5 (eqs. 6.9-6.12, p. 145) from the engine's rows,
    on the file's masses: projectile AWI and target AWR in neutron masses (pp. 39, 24), the neutron mass of p. 394;
    E the lab energy (eV), mu the c.m. cosine."""

    def __init__(self, material):
        from .. import physics
        h = material.require_header()
        self._physics = physics
        self.z1 = h.projectile_za // 1000
        self.z2 = int(round(h.za)) // 1000
        self.m1 = h.awi * M_NEUTRON_U
        self.m2 = h.awr * M_NEUTRON_U

    def _ecm(self, energies_ev):
        return self._physics.lab_to_cm_energy(self.m1, self.m2, np.asarray(energies_ev, float) * 1e-6)

    def eta(self, energies_ev):
        return np.atleast_1d(self._physics.sommerfeld_eta(self.z1, self.z2, self.m1, self.m2, self._ecm(energies_ev)))

    def sigma(self, energies_ev, mu, identical=False, two_s=0):
        theta = np.degrees(np.arccos(np.clip(np.asarray(mu, float), -1, 1)))
        e = self._ecm(energies_ev)
        if identical:
            return self._physics.mott_cross_section(self.z1, self.m1, two_s, e, theta)["value"]
        return self._physics.rutherford_cross_section(self.z1, self.z2, self.m1, self.m2, e, theta)["value"]


def _refuse(what, allowed, fix):
    return FormatRefused(what=what, allowed=allowed, fix=fix)


def request_of(session, data):
    """A :class:`CurveRequest` from a session and one of its ``evaluate(x).data[i]`` (the data set's pairs from the
    session's description, its points and the engine's lab energies from ``data``)."""
    described = [d for d in list(session.datasets) + list(session.extrapolations) if d.name == data.name]
    if not described:
        raise _refuse(f"no data set or extrapolation named {data.name!r} in the session",
                      "an entry of Session.evaluate(x).data of this session", "s.evaluate(x).data lists them")
    d = described[0]
    by_key = {p.key: p for p in session.pairs}
    entrance = PairSpec.of(by_key[d.entrance_key])
    exit_ = None if d.exit_key < 0 else PairSpec.of(by_key[d.exit_key])
    lab = data.extra.get("energy_lab") if isinstance(data.extra, dict) else None
    relativistic = False
    try:
        from ..project import Project
        relativistic = Project(session.file).config().get("relativistic", "0").strip() in ("1", "true")
    except Exception:  # noqa: BLE001 - a project that cannot be read again keeps the default
        pass
    return CurveRequest(data.quantity, tuple(np.asarray(data.energies, float)), tuple(np.asarray(data.angles, float)),
                        entrance, exit_, tuple(np.asarray(lab, float)) if lab is not None else None, data.name,
                        relativistic)


def _lab_energies_ev(request, projectile_za):
    """The lab energies (eV) of the file's projectile at the request's points."""
    e = request.entrance
    if e.za1 == projectile_za and request.energies_lab is not None:
        return np.asarray(request.energies_lab, float) * 1e6
    from .. import physics
    beam, target = (e.m1_u, e.m2_u) if e.za1 == projectile_za else (e.m2_u, e.m1_u)
    return np.asarray(physics.cm_to_lab_energy(beam, target, np.asarray(request.energies, float),
                                                request.relativistic), float) * 1e6


def _particle_mt(material, request):
    """The MT of the exit pair's reaction: the family of its first particle, k by QM - QI = the excitation."""
    x = request.exit
    if x.light_excitation:
        raise _refuse("an exit pair with an excited light particle", "a light particle in its ground state "
                      "(the MT families of Appendix B, pp. 339-341)", "compare with a data set of another exit pair")
    first = names.family_of(x.za1)
    if first is None:
        raise _refuse(f"an exit particle of ZA={x.za1}", "n, p, d, t, 3He or alpha (Appendix B)",
                      "the file has no family of reactions for it")
    candidates = [mt for (mf, mt) in material.sections if mf == 3 and first <= mt < first + 49]
    for mt in candidates:
        xs = material.cross_section(mt)
        ex = (xs.qm - xs.qi) * 1e-6
        if abs(ex - x.excitation) <= max(1e-3, 1e-3 * x.excitation):
            return mt
    raise _refuse(f"no MF=3 section of the family MT {first}-{first + 48} with the residual at "
                  f"{x.excitation:g} MeV (QM - QI, p. 105)",
                  "the sections: " + ", ".join(f"MT={mt}" for mt in candidates) if candidates else "none in the file",
                  "Material.summary() lists the reactions the file holds")


def curves(material, request, session=None):
    """The file's values at the request's points (see the module's text).  ``request``: a :class:`CurveRequest`, or
    an ``EvaluatedData`` with ``session`` (the session it came from)."""
    if not isinstance(request, CurveRequest):
        if session is None:
            raise _refuse("a data set without its session", "Material.curves(request) or "
                          "Material.curves(s.evaluate(x).data[i], session=s)", "pass session=")
        request = request_of(session, request)
    h = material.require_header()
    q = request.quantity
    if q not in QUANTITIES:
        raise _refuse(f"the quantity {q!r}", "one of " + ", ".join(QUANTITIES),
                      "an ENDF file holds cross sections and distributions; other quantities are not in it")
    e = request.entrance
    projectile, target = h.projectile_za, int(round(h.za))
    if {e.za1, e.za2} != {projectile, target}:
        raise _refuse(f"the data set's entrance pair (ZA {e.za1} + {e.za2}) is not the file's "
                      f"projectile + target ({projectile} + {target})", f"a data set of ZA {projectile} on {target}",
                      "Material.summary() names the file's projectile and target")
    energies = np.asarray(request.energies, float)
    angles = np.asarray(request.angles, float)
    mu = np.cos(np.radians(angles))
    lab = _lab_energies_ev(request, projectile)
    notes = []
    from .. import physics
    if request.exit is None or request.exit.photon:
        if request.exit is not None:
            raise _refuse("capture to one final state", "capture to every final state (the data set's exit pair "
                          "-1) from MF=3/MT=102", "the yields of each final state (MF=6/MT=102) are not built")
        if q not in ("sigma", "sfactor"):
            raise _refuse(f"{q} of capture", "sigma or sfactor of MF=3/MT=102", "use an integrated data set")
        value, source = material.cross_section(102)(lab), "MF=3/MT=102"
    elif request.exit == e:
        source = "MF=6/MT=2 LAW=5"
        if q in ("sigma", "sfactor"):
            raise _refuse(f"{q} of elastic scattering", "dsigma_domega, ratio_to_rutherford, ratio_to_mott "
                          "(an integrated charged-particle elastic cross section is not defined, p. 145)",
                          "use a differential data set")
        from .distributions import elastic_cross_section
        product = material.distribution(2).product(projectile)
        if product.law != 5:
            raise _refuse(f"MF=6/MT=2 LAW={product.law}", "LAW=5 for charged-particle elastic scattering (p. 144)",
                          "the file's elastic distribution is not of a charged pair")
        sigma_ni = material.cross_section(2) if (3, 2) in material.sections else None
        # the data set's angle is its first particle's; the file's mu the projectile's (the same when it is first)
        mu_projectile = mu if e.za1 == projectile else -mu
        value = elastic_cross_section(product, lab, mu_projectile, CoulombTerms(material), sigma_ni)
    else:
        mt = _particle_mt(material, request)
        xs = material.cross_section(mt)
        source = f"MF=3/MT={mt}"
        if q in ("sigma", "sfactor"):
            value = xs(lab)
        else:
            if (6, mt) not in material.sections:
                raise _refuse(f"{q} of MT={mt} without MF=6/MT={mt}", "an angular distribution in MF=6 (p. 143)",
                              "compare an integrated data set with MF=3")
            from .distributions import two_body_probability
            d = material.distribution(mt)
            if d.lct != 2:
                raise _refuse(f"MF=6/MT={mt} with LCT={d.lct}", "LCT=2 (c.m.) for LAW=2 (p. 143)",
                              "the file breaks the format: pyazr.endf.check_format(path)")
            product = d.product(request.exit.za1)
            if product.law != 2:
                raise _refuse(f"MF=6/MT={mt} product ZAP={product.zap} LAW={product.law}", "LAW=2 (p. 143)",
                              "other laws are not built")
            from .interpolation import interpolate
            t = product.yield_table
            y = interpolate(lab, t.x, t.y, t.nbt, t.interpolation)
            # eq. 6.1 integrated over E': dsigma/dOmega = sigma(E) y(E) p(mu, E) / 2 pi
            value = xs(lab) * y * two_body_probability(product, lab, mu) / (2 * math.pi)
            source += f" x MF=6/MT={mt} LAW=2"
    value = np.asarray(value, float)
    unit = {"sigma": "b", "dsigma_domega": "b/sr", "ratio_to_rutherford": "1", "ratio_to_mott": "1",
            "sfactor": "MeV*b"}[q]
    if q == "sfactor":
        value = value * physics.sfactor_conversion(e.z1, e.z2, e.m1_u, e.m2_u, energies)
    elif q == "ratio_to_rutherford":
        value = value / physics.rutherford_cross_section(e.z1, e.z2, e.m1_u, e.m2_u, energies, angles,
                                                         request.relativistic)["value"]
    elif q == "ratio_to_mott":
        value = value / physics.mott_cross_section(e.z1, e.m1_u, int(round(2 * e.j1)), energies, angles,
                                                   request.relativistic)["value"]
    coverage = np.isfinite(value)
    if not coverage.all():
        notes.append(f"{int((~coverage).sum())} points outside the file's tables (lab energies "
                     f"{lab[~coverage].min() * 1e-6:.6g}-{lab[~coverage].max() * 1e-6:.6g} MeV)")
    return Curve(request.name, q, unit, "c.m.", energies, angles, value, coverage, source, notes)
