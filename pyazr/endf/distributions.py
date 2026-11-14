"""MF=3 cross sections and MF=6 angular distributions evaluated at given energies and angles (ENDF-102, 2023).

- MF=3 (pp. 103-105): HEAD (ZA, AWR), TAB1 (QM, QI, 0, LR; sigma(E)), E the lab energy of the incident particle
  in eV, sigma in barns.  INT=6 needs the threshold T of eq. 2 (p. 26), which the manual names and does not give
  a formula for: 0 for QI > 0, else the lab energy of the projectile at a c.m. energy of -QI, the engine's row
  ``cm_to_lab_energy`` on the file's masses (AWI of MF=1, AWR of the section, x the neutron mass of p. 394).
- MF=6 (pp. 133-147): HEAD (ZA, AWR, JP, LCT, NK); per product a TAB1 (ZAP, AWP, LIP, LAW; yield y(E)) and the
  law's records.  Built: LAW=2 (two-body, p. 143: LANG=0 Legendre, eq. 6.8; LANG=12 p linear in mu; LANG=14 ln p
  linear in mu) and LAW=5 (charged-particle elastic, pp. 144-147: LTP=1 eqs. 6.13-6.14, LTP=2 eqs. 6.15-6.18,
  LTP=12, 14, 15 eqs. 6.19-6.20).  LAW=0, 3, 4 have no records; LAW=1, 6, 7 are read past (their records counted,
  pp. 137-150) and refused when asked for a value.
- Between incident energies (the TAB2 law): the value at fixed mu is interpolated (Cartesian, eq. 7, p. 26); the
  unit-base laws 21-25 are the same as 1-5 for distributions that span the same range of mu (p. 27, note (i)),
  as a Legendre sum and an LAW=5 distribution do; the corresponding-point laws 11-15 are refused.
- LAW=5 with LTP=12, 14, 15 (the manual's reading, the default): sigma_e = sigma_c + sigma_NI(E) P_NI(mu, E), with
  sigma_NI(E) the TAB1 of MF=3/MT=2 under its own interpolation and P_NI(mu, E) the tabulated distribution of MF=6
  under the TAB2's law between incident energies at fixed mu (p. 146: "The angular distribution p_NI is given in
  File 6 as a tabulated function of mu, and sigma_NI(E) in barns is given in File 3"; eqs. 6.19-6.20).  The manual
  states no other reading.  Where sigma_NI changes sign between two incident energies, P_NI = (sigma_e - sigma_c) /
  sigma_NI changes sign too and is large on both sides, and this product of two interpolations differs from
  sigma_e - sigma_c by about its own size there (docs/upgrade/endf.md 13.10).  The other reading, ``"product"``,
  interpolates sigma_NI(E_i) P_NI(mu, E_i) of the panels (the tabulated sigma_e - sigma_c) instead: a named option,
  for files written so that only it holds.
- LAW=5 with LTP=1 and 2: the quantity interpolated between incident energies is sigma_e - sigma_c (the nuclear and
  interference part, b/sr; the manual says nothing of it for these LTP).  The Coulomb parts (eta and the Rutherford or
  Mott cross section) are computed at the energy asked, the engine's: a :class:`pyazr.endf.curves.CoulombTerms`
  given by the caller.

Values outside a table (energy or cosine) are NaN; a tabulated distribution is zero outside its cosines (p. 142).
"""

import dataclasses

import numpy as np
from numpy.polynomial import legendre

from ..errors import FormatRefused
from .interpolation import interpolate
from .records import Tab1, refuse


@dataclasses.dataclass(frozen=True)
class CrossSection:
    """MF=3/MT (pp. 103-104): ZA, AWR, QM and QI (eV), LR, and the table sigma(E) (``table``, a TAB1: lab eV, b);
    ``awi`` the projectile's mass (neutron masses, MF=1) for the threshold of INT=6."""
    mt: int
    za: float
    awr: float
    qm: float
    qi: float
    lr: int
    table: Tab1
    awi: float = 1.0

    @property
    def energy_range(self):
        """(first, last) lab energy of the table, eV."""
        return (self.table.x[0], self.table.x[-1]) if self.table.x else (np.nan, np.nan)

    @property
    def threshold(self):
        """T of INT=6 (eq. 2, p. 26), eV: 0 for QI > 0, else the lab energy at a c.m. energy -QI (the engine's row)."""
        if self.qi > 0:
            return 0.0
        from .. import physics
        from .material import M_NEUTRON_U
        return float(physics.cm_to_lab_energy(self.awi * M_NEUTRON_U, self.awr * M_NEUTRON_U, -self.qi * 1e-6)) * 1e6

    def __call__(self, energies_ev):
        """sigma (b) at lab energies (eV); NaN outside the table."""
        t = self.table
        threshold = self.threshold if 6 in t.interpolation else None
        return interpolate(energies_ev, t.x, t.y, t.nbt, t.interpolation, threshold=threshold)


def read_cross_section(material, mt):
    s = material.section(3, mt)
    c = s.cursor(material.path)
    head = c.cont()
    tab = c.tab1()
    _require_done(material, c, 3, mt)
    h = material.header
    return CrossSection(mt, head.c1, head.c2, tab.c1, tab.c2, tab.l2, tab, h.awi if h is not None else 1.0)


# ---------------------------------------------------------------------------------------------------- MF=6
@dataclasses.dataclass(frozen=True)
class Panel:
    """The record of one incident energy: E (lab eV), the representation flag (LANG of LAW=2, LTP of LAW=5), NL
    and the values A of the LIST (p. 143, p. 147)."""
    energy: float
    flag: int
    nl: int
    values: tuple


@dataclasses.dataclass(frozen=True)
class Product:
    """One product subsection (p. 136): ZAP, AWP, LIP, LAW, the yield y(E) (a TAB1), and for LAW=2 and 5 the TAB2
    (``tab2``: SPI in C1 and LIDP in L1 for LAW=5, p. 147) and the panels; other laws keep no records."""
    zap: int
    awp: float
    lip: int
    law: int
    yield_table: Tab1
    tab2: object = None
    panels: tuple = ()

    @property
    def energies(self):
        return np.array([p.energy for p in self.panels])

    @property
    def lidp(self):
        return self.tab2.l1 if self.tab2 is not None else 0

    @property
    def spin(self):
        return self.tab2.c1 if self.tab2 is not None else 0.0


@dataclasses.dataclass(frozen=True)
class Distribution:
    """MF=6/MT (p. 136): ZA, AWR, JP, LCT (1 lab, 2 c.m., 3, 4) and its products."""
    mt: int
    za: float
    awr: float
    jp: int
    lct: int
    products: tuple

    def product(self, zap=None):
        """The product subsection of ``zap`` (the first with LAW 2 or 5 when None)."""
        for p in self.products:
            if (zap is None and p.law in (2, 5)) or p.zap == zap:
                return p
        raise FormatRefused(what=f"MF=6/MT={self.mt} has no product ZAP={zap}",
                            allowed="its products: " + ", ".join(f"ZAP={p.zap} LAW={p.law}" for p in self.products),
                            fix="Material.summary() lists the products of each section")


def _skip_law(cursor, law, where):
    """Read past the records of LAW=1, 6, 7 (pp. 137-150) so that the next product is found."""
    if law == 1:
        t2 = cursor.tab2()
        for _ in range(t2.nz):
            cursor.list()
    elif law == 6:
        cursor.cont()
    elif law == 7:
        t2 = cursor.tab2()
        for _ in range(t2.nz):
            inner = cursor.tab2()
            for _ in range(inner.nz):
                cursor.tab1()
    elif law not in (0, 3, 4):
        raise FormatRefused(what=f"{where}: LAW={law} is not a law of MF=6",
                            allowed="LAW 0-7 (p. 135-136)", fix="pyazr.endf.check_format(path) reports the record")


def read_distribution(material, mt):
    s = material.section(6, mt)
    c = s.cursor(material.path)
    head = c.cont()
    products = []
    for _ in range(max(head.n1, 0)):
        y = c.tab1()
        zap, awp, lip, law = int(y.c1), y.c2, y.l1, y.l2
        if law in (2, 5):
            t2 = c.tab2()
            panels = []
            for _ in range(t2.nz):
                rec = c.list()
                panels.append(Panel(rec.c2, rec.l1, rec.n2, rec.values))
            products.append(Product(zap, awp, lip, law, y, t2, tuple(panels)))
        else:
            _skip_law(c, law, f"{material.path} MF=6/MT={mt} line {y.line}")
            products.append(Product(zap, awp, lip, law, y))
    _require_done(material, c, 6, mt)
    return Distribution(mt, head.c1, head.c2, head.l1, head.l2, tuple(products))


def _require_done(material, cursor, mf, mt):
    """Lines left after the last record mean a count smaller than the record (pp. 34-35): refused."""
    if not cursor.done:
        line = cursor.lines[cursor.i].number
        raise refuse(material.path, line, None, f"{len(cursor.lines) - cursor.i} lines",
                     f"MF={mf}/MT={mt}: lines after the last record (a count of a record is smaller than its lines)",
                     "counts (NPL, NR, NP, NZ, NK) that match the lines of the section", 34, "ENDF_COUNT")


def _panel_law(product, where):
    laws = []
    for law in product.tab2.interpolation:
        if 21 <= law <= 25:
            law -= 20
        if not 1 <= law <= 5:
            raise FormatRefused(what=f"{where}: interpolation law {law} between incident energies",
                                allowed="INT 1-5, or 21-25 (unit base, the same as 1-5 for these distributions, p. 27)",
                                fix="corresponding-point laws (11-15) are not built")
        laws.append(law)
    return laws


def _between_panels(product, energies, values_at, where):
    """Interpolate ``values_at(panel_index, k)`` (the value of panel i at point k) between incident energies at each
    point's energy (Cartesian, eq. 7 p. 26)."""
    laws = _panel_law(product, where)
    e = product.energies
    out = np.full(len(energies), np.nan)
    for k, energy in enumerate(energies):
        if not e.size or not (e[0] <= energy <= e[-1]):
            continue
        i = int(np.searchsorted(e, energy, side="right")) - 1
        if i >= e.size - 1 or e[i] == energy:
            out[k] = values_at(min(i, e.size - 1), k)
            continue
        ys = [values_at(i, k), values_at(i + 1, k)]
        out[k] = interpolate(energy, e[i:i + 2], ys, (2,), (laws_for(product.tab2.nbt, laws, i),))
    return out


def laws_for(nbt, laws, interval):
    for bound, law in zip(nbt, laws):
        if interval + 1 < bound:
            return law
    return laws[-1]


def _legendre_sum(coefficients, mu, first=0):
    """sum_l (2l+1)/2 c_l P_l(mu), l from ``first``."""
    c = np.zeros(len(coefficients) + first)
    for l, a in enumerate(coefficients, start=first):
        c[l] = (2 * l + 1) / 2.0 * a
    return legendre.legval(mu, c)


def two_body_probability(product, energies_ev, mu, where="MF=6"):
    """p(mu, E) of a LAW=2 product (eq. 6.8, p. 143), normalised so that its integral over mu is 1."""
    def at(i, k):
        p = product.panels[i]
        if p.flag == 0:
            return 0.5 + _legendre_sum(p.values[:p.nl], mu[k], first=1)
        if p.flag in (12, 14):
            m, v = p.values[0::2], p.values[1::2]
            if not m[0] <= mu[k] <= m[-1]:
                return 0.0
            return interpolate(mu[k], m, v, (len(m),), (2 if p.flag == 12 else 4,))
        raise FormatRefused(what=f"{where}: LANG={p.flag} at E={p.energy:g} eV", allowed="LANG 0, 12, 14 (p. 143)",
                            fix="the file breaks the format: pyazr.endf.check_format(path)")
    return _between_panels(product, np.asarray(energies_ev, float), at, where)


#: The readings of a tabulated LAW=5 distribution (LTP=12, 14, 15) between incident energies (the module's text).
ELASTIC_READINGS = ("manual", "product")


def elastic_cross_section(product, energies_ev, mu, coulomb, sigma_ni=None, where="MF=6/MT=2", reading="manual"):
    """dsigma/dOmega (b/sr, c.m.) of LAW=5 at lab energies (eV) and c.m. cosines (pp. 144-147).

    ``coulomb``: the Coulomb terms (``eta(E)``, ``sigma(E, mu)``; distinguishable or identical by LIDP);
    ``sigma_ni``: sigma_NI(E) of MF=3/MT=2 (a callable, eV -> b) for LTP=12-15; ``reading``: for LTP=12-15,
    ``"manual"`` (sigma_NI(E) of MF=3 times P_NI(mu, E) interpolated between the panels, p. 146) or ``"product"``
    (sigma_NI(E_i) P_NI(mu, E_i) of the panels interpolated), the module's text."""
    if reading not in ELASTIC_READINGS:
        raise FormatRefused(what=f"{where}: reading {reading!r} of a tabulated LAW=5 distribution",
                            allowed=" or ".join(repr(r) for r in ELASTIC_READINGS),
                            fix="leave the reading at its default, the manual's (p. 146)")
    energies_ev = np.asarray(energies_ev, float)
    mu = np.asarray(mu, float)
    identical = product.lidp == 1
    two_s = int(round(2 * product.spin))

    def nuclear_and_interference(i, k):
        p = product.panels[i]
        m = mu[k]
        if p.flag == 1:
            eta = float(coulomb.eta(np.array([p.energy]))[0])
            nl = p.nl
            nb = nl + 1 if identical else 2 * nl + 1
            b = p.values[:nb]
            a = np.array(p.values[nb:nb + 2 * (nl + 1)][0::2]) + 1j * np.array(p.values[nb:nb + 2 * (nl + 1)][1::2])
            terms = [(2 * l + 1) / 2.0 * a[l] * legendre.legval(m, [0] * l + [1]) for l in range(nl + 1)]
            if not identical:  # eq. 6.13
                interference = -2 * eta / (1 - m) * np.real(np.exp(1j * eta * np.log((1 - m) / 2)) * sum(terms))
                return interference + _legendre_sum(b, m)
            # eq. 6.14
            s = sum(((1 + m) * np.exp(1j * eta * np.log((1 - m) / 2))
                     + (-1) ** l * (1 - m) * np.exp(1j * eta * np.log((1 + m) / 2))) * terms[l] for l in range(nl + 1))
            nuclear = sum((4 * l + 1) / 2.0 * b[l] * legendre.legval(m, [0] * (2 * l) + [1]) for l in range(len(b)))
            return -2 * eta / (1 - m * m) * np.real(s) + nuclear
        if p.flag == 2:  # eqs. 6.15-6.18
            c = p.values[:p.nl + 1]
            if not identical:
                return _legendre_sum(c, m) / (1 - m)
            return sum((4 * l + 1) / 2.0 * c[l] * legendre.legval(m, [0] * (2 * l) + [1])
                       for l in range(len(c))) / (1 - m * m)
        if p.flag in (12, 14, 15):  # eqs. 6.19-6.20: P_NI(mu, E_i), times sigma_NI(E_i) for the product reading
            ms, v = p.values[0::2], p.values[1::2]
            if not ms[0] <= m <= ms[-1]:
                return 0.0
            p_ni = interpolate(m, ms, v, (len(ms),), (4 if p.flag == 14 else 2,))
            return p_ni if reading == "manual" else float(sigma_ni(np.array([p.energy]))[0]) * p_ni
        raise FormatRefused(what=f"{where}: LTP={p.flag} at E={p.energy:g} eV", allowed="LTP 1, 2, 12, 14, 15 (p. 146)",
                            fix="the file breaks the format: pyazr.endf.check_format(path)")

    flags = {p.flag for p in product.panels}
    if flags & {12, 14, 15} and sigma_ni is None:
        raise FormatRefused(what=f"{where}: LTP={sorted(flags)} needs sigma_NI(E) of MF=3/MT=2",
                            allowed="a material with MF=3/MT=2", fix="the file lacks MF=3/MT=2")
    # interpolated between the panels at fixed mu: sigma_e - sigma_c for LTP 1 and 2 and for the product reading;
    # P_NI for the manual's reading of LTP 12-15, then times sigma_NI(E) of MF=3 under its own interpolation
    part = _between_panels(product, energies_ev, nuclear_and_interference, where)
    if reading == "manual" and flags & {12, 14, 15}:
        if not flags <= {12, 14, 15}:
            raise FormatRefused(what=f"{where}: LTP={sorted(flags)} mixed in one section",
                                allowed="one representation, or LTP 12, 14, 15 together (p. 147)",
                                fix="the file breaks the format: pyazr.endf.check_format(path)")
        part = np.asarray(sigma_ni(energies_ev), float) * part
    return coulomb.sigma(energies_ev, mu, identical=identical, two_s=two_s) + part
