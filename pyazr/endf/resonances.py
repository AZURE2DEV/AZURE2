"""MF=2/MT=151 in the R-Matrix Limited format (LRF=7): read into a plain record, written back (ENDF-102, 2023 edition,
cited by printed page).

The section (pp. 64-65): HEAD ``ZA, AWR, 0, 0, NIS, 0``; per isotope a CONT ``ZAI, ABN, 0, LFW, NER, 0``; per energy
range a CONT ``EL, EH, LRU, LRF, NRO, NAPS`` and its subsection.  The LRF=7 subsection (pp. 71-74): a CONT ``0, 0,
IFG, KRM, NJS, KRL``; a LIST of the NPP particle pairs, 12 numbers each (MA, MB, ZA, ZB, IA, IB, Q, PNT, SHF, MT,
PA, PB); then per spin group a LIST ``AJ, PJ, KBK, KPS, 6*NCH, NCH`` of six numbers per channel (PPI, L, SCH, BND,
APE, APT) and a LIST ``0, 0, 0, NRS, 6*NX, NX`` of the levels, each ER and its NCH values GAM starting a new line
(NX = (NCH/6 + 1) NRS, NX = 1 with no level, p. 74).

Units and frame, as the file holds them: masses in neutron masses (p. 72), Q and ER in eV (pp. 72-73), ER the
laboratory energy of the incident particle (p. 66), GAM the reduced-width amplitude in eV^1/2 for IFG=1 (p. 73) on
the same scale (eqs. D.39, D.44, pp. 355-356), radii in 1e-12 cm (p. 73).  The conversions to AZURE2's frame and
units are :mod:`pyazr.endf.model`'s; nothing here converts.

What AZURE2 cannot read (docs/upgrade/endf.md section 3.2) is refused here, each with its three parts (what,
allowed, fix) as :class:`pyazr.FormatRefused`.
"""

import dataclasses
import math
import os

from ..errors import FormatRefused
from . import names
from .records import RecordWriter

#: A particle whose mass is below this (neutron masses) is a photon: MA = 0 for a gamma pair (p. 73).
PHOTON_MASS = 1e-6


@dataclasses.dataclass(frozen=True)
class ParticlePair:
    """One particle pair of an LRF=7 subsection (pp. 72-73), as the file holds it: masses ``ma``, ``mb`` (neutron
    masses), charges ``za``, ``zb``, spins ``ia``, ``ib`` (the sign of a non-zero spin is the parity), ``q`` (eV),
    ``pnt``, ``shf``, ``mt``, parities ``pa``, ``pb`` (used when the spin is 0)."""
    ma: float
    mb: float
    za: float
    zb: float
    ia: float
    ib: float
    q: float
    pnt: int
    shf: int
    mt: int
    pa: float
    pb: float

    @property
    def is_photon(self):
        return self.ma < PHOTON_MASS or self.mt == 102

    @staticmethod
    def _parity(spin, parity):
        """p. 72: IA "spin (and parity, if non-zero)"; PA "used only in the case where IA is zero and the parity is
        negative"."""
        if spin != 0.0:
            return 1 if spin > 0 else -1
        return -1 if parity < 0 else 1

    @property
    def spin_a(self):
        return abs(self.ia)

    @property
    def spin_b(self):
        return abs(self.ib)

    @property
    def parity_a(self):
        return self._parity(self.ia, self.pa)

    @property
    def parity_b(self):
        return self._parity(self.ib, self.pb)

    def values(self):
        return [self.ma, self.mb, self.za, self.zb, self.ia, self.ib, self.q, self.pnt, self.shf, self.mt, self.pa,
                self.pb]


@dataclasses.dataclass(frozen=True)
class Channel:
    """One channel of a spin group (p. 73): ``pair`` the particle-pair number PPI (1-based), ``l``, channel spin
    ``s``, boundary condition ``boundary`` (BND), radii ``radius_phase`` (APE) and ``radius`` (APT), 1e-12 cm."""
    pair: int
    l: int
    s: float
    boundary: float
    radius_phase: float
    radius: float


@dataclasses.dataclass(frozen=True)
class Level:
    """One level of a spin group (p. 73): ``energy`` ER (eV, lab energy of the incident particle) and its values
    ``gam`` per channel (eV^1/2 for IFG=1, eV for IFG=0)."""
    energy: float
    gam: tuple


@dataclasses.dataclass(frozen=True)
class SpinGroup:
    """One J^pi spin group (p. 74): ``J``, ``parity``, its channels and levels (``kbk``, ``kps`` are 0: other values
    are refused when read)."""
    J: float
    parity: int
    channels: tuple
    levels: tuple
    kbk: int = 0
    kps: int = 0

    @property
    def jpi(self):
        twice = round(2 * self.J)
        j = f"{twice // 2}" if twice % 2 == 0 else f"{twice}/2"
        return j + ("+" if self.parity > 0 else "-")


@dataclasses.dataclass(frozen=True)
class LRF7:
    """An R-Matrix Limited subsection with what surrounds it in MF=2/MT=151 (pp. 64-65, 71-74), in the file's units
    and frame: ``za``, ``awr`` (HEAD), ``zai``, ``abn`` (isotope), ``el``, ``eh`` (eV), ``naps`` (range), ``ifg``,
    ``krm``, ``krl``, the particle ``pairs`` and the spin ``groups``.  ``source`` names where it was read from (empty
    for one built from a session).  ``print(record)`` lists it."""
    za: float
    awr: float
    zai: float
    abn: float
    el: float
    eh: float
    naps: int
    ifg: int
    krm: int
    krl: int
    pairs: tuple
    groups: tuple
    source: str = ""

    @property
    def entrance(self):
        """The 1-based number of the incident pair: MT=2 (p. 71: "MA = 1.0 with Q = 0 or MT = 2")."""
        for k, p in enumerate(self.pairs, start=1):
            if p.mt == 2:
                return k
        for k, p in enumerate(self.pairs, start=1):
            if abs(p.ma - 1.0) < 1e-3 and p.q == 0.0:
                return k
        raise FormatRefused(what=f"{self.source or 'the LRF=7 subsection'} has no incident pair",
                            allowed="a particle pair with MT=2, or MA = 1 with Q = 0 (ENDF-102 p. 71)",
                            fix="pyazr.endf.read(path).summary() shows the pairs")

    def __str__(self):
        out = [f"LRF=7 (R-Matrix Limited, {names.RML_FORMULAS.get(self.krm, f'KRM={self.krm}')}), IFG={self.ifg}, "
               f"KRL={self.krl}, ZA={self.za:g}, EL={self.el:g} eV, EH={self.eh:g} eV"
               + (f"  [{self.source}]" if self.source else "")]
        for k, p in enumerate(self.pairs, start=1):
            out.append(f"  pair {k}: MT={p.mt} {names.reaction_name(p.mt)}, MA={p.ma:.7g} MB={p.mb:.7g} m_n, "
                       f"Z={p.za:g}+{p.zb:g}, spins {p.spin_a:g}{'+' if p.parity_a > 0 else '-'} "
                       f"{p.spin_b:g}{'+' if p.parity_b > 0 else '-'}, Q={p.q:.10g} eV, PNT={p.pnt} SHF={p.shf}")
        for g in self.groups:
            chans = ", ".join(f"(pair {c.pair}, l={c.l}, s={c.s:g}, B={c.boundary:g}, a={c.radius * 10:.6g} fm)"
                              for c in g.channels)
            out.append(f"  {g.jpi}: {len(g.levels)} level(s); channels {chans}")
            for lv in g.levels:
                out.append(f"    ER={lv.energy:.10g} eV  GAM=" + " ".join(f"{v:.8g}" for v in lv.gam))
        return "\n".join(out)

    # -- writing --------------------------------------------------------------------------------------------------

    def rows(self, number_form="widest"):
        """MF=2/MT=151 as rows of 66 columns (HEAD, isotope, range, the LRF=7 subsection; pp. 64-65, 71-74), one
        isotope and one energy range (NIS = NER = 1, LFW = 0, NRO = 0)."""
        w = RecordWriter(number_form)
        w.cont(self.za, self.awr, 0, 0, 1, 0)
        w.cont(self.zai, self.abn, 0, 0, 1, 0)
        w.cont(self.el, self.eh, 1, 7, 0, self.naps)
        w.cont(0.0, 0.0, self.ifg, self.krm, len(self.groups), self.krl)
        values = [v for p in self.pairs for v in p.values()]
        w.list(0.0, 0.0, len(self.pairs), 0, 2 * len(self.pairs), values)
        for g in self.groups:
            cv = [v for c in g.channels for v in (float(c.pair), float(c.l), c.s, c.boundary, c.radius_phase, c.radius)]
            w.list(g.J * (1 if g.parity > 0 else -1) if g.J != 0 else 0.0, float(g.parity) if g.J == 0 else 0.0,
                   g.kbk, g.kps, len(g.channels), cv)
            per_level = 6 * (len(g.channels) // 6 + 1)  # each level starts a new line (p. 74)
            rv = []
            for lv in g.levels:
                row = [lv.energy] + list(lv.gam)
                rv.extend(row + [0.0] * (per_level - len(row)))
            nx = len(rv) // 6 if g.levels else 1
            w.list(0.0, 0.0, 0, len(g.levels), nx, rv if g.levels else [0.0] * 6)
        return w.rows


def _refuse(source, what, allowed, fix):
    return FormatRefused(what=f"{source}: {what}", allowed=allowed, fix=fix)


def read_resonances(material, isotope=0, energy_range=0):
    """MF=2/MT=151 of ``material``: the LRF=7 subsection of isotope ``isotope`` and energy range ``energy_range``
    (0-based) as an :class:`LRF7`.  Refused (:class:`pyazr.FormatRefused`), each with what, allowed and fix: no
    MF=2, a range that is not LRF=7 or has NRO=1, KRM other than 4, KRL=1, IFG other than 0 and 1, a photon or
    fission pair, KBK > 0 or KPS > 0, APE != APT or two radii in one pair, SHF=0 or PNT=-1 on a particle channel.
    An IFG=0 record keeps its widths (eV); :func:`pyazr.endf.to_project` converts them (p. 73)."""
    source = f"{os.path.basename(material.path)} MAT={material.mat} MF=2/MT=151"
    section = material.section(2, 151)
    c = section.cursor(material.path)
    head = c.cont()
    nis = head.n1
    if not 0 <= isotope < max(nis, 1):
        raise _refuse(source, f"has {nis} isotope(s); isotope={isotope} asked", f"isotope = 0 ... {nis - 1}",
                      "Material.summary() lists the isotopes")
    chosen = None
    for i in range(nis):
        iso = c.cont()
        for r in range(iso.n1):
            rng = c.cont()
            if (i, r) == (isotope, energy_range):
                chosen = (iso, rng)
                break
            if rng.l1 == 1 and rng.l2 == 7:
                _skip_lrf7(c)
            else:
                raise _refuse(source, f"isotope {i}, range {r} is LRU={rng.l1} LRF={rng.l2}; the ranges before the "
                              "one asked must be LRF=7 to be passed over", "isotope and range of an LRF=7 file",
                              "ask for the first range (energy_range=0)")
        if chosen:
            break
        if i == isotope:
            raise _refuse(source, f"isotope {isotope} has {iso.n1} range(s); energy_range={energy_range} asked",
                          f"energy_range = 0 ... {iso.n1 - 1}", "Material.summary() lists the ranges")
    iso, rng = chosen
    if rng.l1 != 1 or rng.l2 != 7:
        raise _refuse(source, f"the range is LRU={rng.l1}, LRF={rng.l2} "
                      f"({names.RESOLVED_FORMATS.get(rng.l2, '?') if rng.l1 == 1 else 'not resolved'})",
                      "LRU=1, LRF=7 (R-Matrix Limited, pp. 71-74)",
                      "AZURE2 reads R-matrix parameters only; Material.curves() gives the file's cross sections")
    if rng.n1 != 0:
        raise _refuse(source, "NRO=1: an energy-dependent scattering radius (p. 63)", "NRO=0",
                      "AZURE2 has one radius per pair; give the radius as a constant in the evaluation")
    sub = c.cont()
    ifg, krm, njs, krl = sub.l1, sub.l2, sub.n1, sub.n2
    if krm == 3:
        raise _refuse(source, "KRM=3 (Reich-Moore): the first pair is gamma + compound nucleus, whose widths are "
                      "eliminated into the level energies (pp. 73-74, Appendix D.1.7)", "KRM=4, the full R-matrix "
                      "(p. 74)", "AZURE2's Reich-Moore option would take the eliminated widths, but turning them into "
                      "the Brune basis of a start model is not built; compare the file's cross sections with the "
                      "data instead")
    if krm != 4:
        raise _refuse(source, f"KRM={krm} ({names.RML_FORMULAS.get(krm, 'unknown')})",
                      "KRM=4, the full R-matrix (p. 74)", "a Breit-Wigner file is not read into a start model; "
                      "compare the file's cross sections with the data instead")
    if krl != 0:
        raise _refuse(source, "KRL=1, relativistic kinematics (p. 71: a flag without equations)", "KRL=0",
                      "the convention of KRL=1 is not decided for AZURE2; read a non-relativistic evaluation")
    if ifg not in (0, 1):
        raise _refuse(source, f"IFG={ifg}", "IFG=0 (widths in eV) or IFG=1 (reduced-width amplitudes in eV^1/2), "
                      "p. 73", "correct IFG in the file")
    plist = c.list()
    pairs = []
    for k in range(plist.l1):
        v = plist.values[12 * k:12 * k + 12]
        pairs.append(ParticlePair(v[0], v[1], v[2], v[3], v[4], v[5], v[6], int(round(v[7])), int(round(v[8])),
                                  int(round(v[9])), v[10], v[11]))
    for k, p in enumerate(pairs, start=1):
        if p.is_photon:
            raise _refuse(source, f"pair {k} is a photon pair (MT={p.mt}): with KRM=4 its channels enter the level "
                          "matrix with P = 1 (p. 74, p. 356)", "particle pairs only (AZURE2's photon channels stay "
                          "out of the level matrix: radiative widths are not part of the total width)",
                          "remove the photon pair from MF=2; capture belongs in MF=3/MT=102")
        if p.mt in (18, 19, 20, 21, 38):
            raise _refuse(source, f"pair {k} is a fission pair (MT={p.mt})", "particle pairs of two nuclei",
                          "AZURE2 has no fission channel")
    groups = []
    for _ in range(njs):
        clist = c.list()
        aj, pj, kbk, kps, nch = clist.c1, clist.c2, clist.l1, clist.l2, clist.n2
        J = abs(aj)
        parity = (1 if aj > 0 else -1) if aj != 0 else (1 if pj >= 0 else -1)
        jpi = SpinGroup(J, parity, (), ()).jpi
        if kbk or kps:
            raise _refuse(source, f"spin group {jpi} has KBK={kbk}, KPS={kps}: a background R-matrix or tabulated "
                          "phase shifts (pp. 74-77)", "KBK=0, KPS=0",
                          "represent the background by distant levels in the evaluation")
        channels = []
        for k in range(nch):
            v = clist.values[6 * k:6 * k + 6]
            ch = Channel(int(round(v[0])), int(round(v[1])), v[2], v[3], v[4], v[5])
            if not 1 <= ch.pair <= len(pairs):
                raise _refuse(source, f"spin group {jpi}, channel {k + 1} names pair {ch.pair}",
                              f"a pair 1 ... {len(pairs)}", "correct PPI in the file")
            p = pairs[ch.pair - 1]
            if ch.radius_phase != ch.radius:
                raise _refuse(source, f"spin group {jpi}, channel {k + 1}: APE={ch.radius_phase:g} differs from "
                              f"APT={ch.radius:g}", "APE = APT (AZURE2 uses one radius for P, S and the phase)",
                              "write one radius per channel")
            if p.shf != 1 or p.pnt == -1:
                raise _refuse(source, f"spin group {jpi}, channel {k + 1}: pair {ch.pair} has SHF={p.shf}, "
                              f"PNT={p.pnt}", "SHF=1 and PNT=1 (or 0) on particle channels: the shift function and "
                              "the penetrability computed (p. 72)", "AZURE2 computes both; set the flags so")
            channels.append(ch)
        llist = c.list()
        nrs = llist.l2
        per = 6 * (nch // 6 + 1) if nrs else 0
        levels = tuple(Level(llist.values[per * i], tuple(llist.values[per * i + 1:per * i + 1 + nch]))
                       for i in range(nrs))
        groups.append(SpinGroup(J, parity, tuple(channels), levels, kbk, kps))
    radii = {}
    for g in groups:
        for ch in g.channels:
            if radii.setdefault(ch.pair, ch.radius) != ch.radius:
                raise _refuse(source, f"pair {ch.pair} has channels of radii {radii[ch.pair]:g} and {ch.radius:g} "
                              "(1e-12 cm)", "one radius per pair", "AZURE2 has one channel radius per pair")
    return LRF7(head.c1, head.c2, iso.c1, iso.c2, rng.c1, rng.c2, rng.n2, ifg, krm, krl, tuple(pairs), tuple(groups),
                source)


def _skip_lrf7(c):
    """Reads past one LRF=7 subsection (pp. 73-74) without keeping it."""
    sub = c.cont()
    c.list()
    for _ in range(sub.n1):
        c.list()
        c.list()


def lab_factor(pair):
    """E_lab / E_cm = (m_A + m_B) / m_B of an incident pair (eq. D.39, p. 355; the masses of the file)."""
    return (pair.ma + pair.mb) / pair.mb


def amplitude_scale(pair):
    """GAM (eV^1/2, lab scale) / gamma (MeV^1/2, c.m. scale): gamma^2 / (E_lambda - E) unchanged when both energies
    are lab energies (eq. D.44, p. 356), so gamma_lab = gamma_cm sqrt(E_lab / E_cm), times 1e3 for MeV^1/2 -> eV^1/2."""
    return math.sqrt(lab_factor(pair)) * 1e3
