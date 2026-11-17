"""A session's model written as an ENDF-6 material (ENDF-102, 2023 edition, cited by printed page;
docs/upgrade/endf.md sections 4 and 11).

What is written:
- MF=1/MT=451 (pp. 38-42): the header (LRP=2, NFOR=6, the projectile's AWI, EMAX, NSUB), the three standard text
  lines and HSUB, the text (what wrote the file, the project, the sentence on MF=2 and MF=3, the photon widths, the
  losses), the directory made from the sections (NC without SEND, p. 41).
- MF=2/MT=151 (pp. 64-65, 71-74): the particle channels in the formal basis, LRF=7 KRM=4 IFG=1
  (:func:`pyazr.endf.from_session`, the engine's door).
- MF=3 per reaction (pp. 103-104): every particle exit (MT 600 + k, 650 + k, ..., 800 + k, 50 + k) and the total
  capture (MT=102) from the engine on one grid (:mod:`pyazr.endf.grid`: the engine's refinement of an adaptive
  extrapolation grid, started at the level energies and thresholds, checked at quarter points) to the tolerance of
  linear interpolation, lab energies of the projectile in eV,
  barns, INT=2.  Before the file is written, MF=2 is read back into a project and its
  reaction cross sections compared with MF=3 at every grid point (the PI's decision: LRP=2 with the statement that
  MF=2 reproduces MF=3); a difference above ``MF2_TOLERANCE`` refuses the write.

- MF=6 (pp. 133-147): LAW=2 of every exit with its residual (LAW=4), the primary photons of MT=102, and
  charged-particle elastic scattering as LAW=5 in one of :data:`ELASTIC_FORMS`: the nuclear amplitude expansion
  LTP=1 with 1.0 in MF=3/MT=2 (:mod:`pyazr.endf.amplitudes`, the engine's door ``nuclear_amplitudes``), or LTP=12
  with sigma_NI in MF=3/MT=2 (:mod:`pyazr.endf.elastic`).

What is not written is listed in the MF=1 text and in the report (never dropped silently): every section a reaction
of the project needs and the file lacks, the photon channels (MF=1 text), and what :class:`pyazr.endf.Losses`
names.
"""

import dataclasses
import datetime
import os
import tempfile

import numpy as np

from ..errors import FormatRefused
from . import amplitudes as amplitude_section
from . import covariance as covariance_section
from .covariance import FormalCovariance
from . import elastic as elastic_section
from . import engine_values, grid, names
from . import sections as section_rows
from .material import M_NEUTRON_U
from .model import EV_PER_MEV, from_session, to_project
from .records import RecordWriter, NUMBER_FORMS

#: LRP of MF=1 (pp. 38-39): 2, "parameters given in File 2 but not to be used in calculating cross sections, which
#: are assumed to be represented completely in File 3" (the PI's decision; p. 42 reserves it for derived files, which
#: CHECKR reports).  One named value, so that another choice is one line.
LRP = 2
#: MF=2 must reproduce MF=3 to this (relative, at every grid point where MF=3 is above 1e-12 of its largest value).
MF2_TOLERANCE = 1e-5
#: Why a section a project's reactions need is not written: (MF, role) -> the reason.  The writer lists every section
#: that a reaction of the project needs and the file does not hold (``Written.not_written`` and the MF=1 text).
NOT_WRITTEN_REASONS = {
    (3, "neutral_elastic"): "elastic scattering of a neutral projectile (MF=3/MT=2 with MF=4 or MF=6 LAW=2) is not "
                            "written",
    (6, "neutral_elastic"): "elastic scattering of a neutral projectile (MF=3/MT=2 with MF=4 or MF=6 LAW=2) is not "
                            "written",
    (3, "elastic"): "charged-particle elastic scattering: sigma_NI was not computed",
    (6, "elastic"): "charged-particle elastic scattering: LAW=5 was not computed",
    (3, "exit"): "the exit's cross section was not computed",
    (6, "exit"): "the exit's angular distribution was not computed",
    (3, "capture"): "the capture cross section was not computed",
    (6, "capture"): "the primary photons were not computed"}
#: The limits of a written MF=6 section, (incident energies, cosines of a LAW=5 panel), by name.  ``"checkr"``, the
#: default (the PI's decision): the limits of CHECKR 8.24 (``checkr.f``: NES5MAX = 200 incident energies, File 5's,
#: applied to MF=6; NANGMAX = 201 angular points, File 4's), which a file must meet to be checked by it.
#: ``"manual"``: ENDF-102 Appendix G (p. 391): MF=6 NE 2,000; no limit of MF=6 on the cosines of LAW=5, so the
#: general mesh size NP 10,000.  Where the tolerance asks for more, a section is thinned by the grid's rule to the
#: limit and the tolerance it reaches is stated (MF=1, ``Written.accuracy``).
MF6_LIMITS = {"checkr": (200, 201), "manual": (2000, 10000)}
LIMIT_WORDS = {"checkr": "CHECKR's limit of", "manual": "the manual's limit of"}
#: The largest c.m. cosine of the LAW=5 table (p. 146: "as close to 1.0 as possible"; sigma_e - sigma_c grows as
#: 1/(1 - mu) towards 1): about 11.5 degrees.  Below it the file's elastic cross section is the Coulomb term alone.
ELASTIC_MAX_COSINE = 0.98
#: The forms of charged-particle elastic scattering (MF=6/MT=2 LAW=5, pp. 145-147): ``"amplitudes"``, the nuclear
#: amplitude expansion LTP=1 with the engine's coefficients and 1.0 in MF=3/MT=2 (p. 150: "Whenever possible, the
#: nuclear amplitude expansion should be used with LAW=5"); ``"tabulated"``, P_NI tabulated in mu, LTP=12, with
#: sigma_NI in MF=3/MT=2.
ELASTIC_FORMS = ("amplitudes", "tabulated")
#: The form written when none is asked (one named value: the PI decides it).
ELASTIC_FORM = "amplitudes"


@dataclasses.dataclass(frozen=True)
class ElasticReport:
    """The elastic section's energies (:mod:`pyazr.endf.elastic`): ``energies`` written, ``candidates`` the engine
    computed, ``brackets`` ((E_low, E_high) lab MeV) of the sign changes of sigma_NI, inside which the manual's
    reading of the file is wrong by up to sigma_e - sigma_c, ``bracket_width`` (their largest width relative to E),
    ``reached`` (the largest deviation by the grid's rule under the manual's reading at an energy computed and not
    written), ``checked`` (every written interval checked at its quarter points)."""
    energies: int
    candidates: int
    brackets: tuple
    bracket_width: float
    reached: float
    checked: bool = True
    cosines_reached: float = 0.0
    tolerance: float = 0.0
    limits: tuple = (None, None)
    limits_name: str = ""

    def lines(self):
        """The statement of MF=1 and of the report, in words."""
        limited = self.limits[0] is not None and self.energies >= self.limits[0] and self.reached > self.tolerance
        cut = self.limits[1] is not None and self.cosines_reached > self.tolerance
        words = LIMIT_WORDS.get(self.limits_name, "the limit of")
        out = [f"Elastic scattering (MF=3/MT=2 and MF=6/MT=2) on {self.energies} energies chosen of {self.candidates} "
               f"the engine computed{f' ({words} {self.limits[0]} per MF=6 section)' if limited else ''}, so that "
               f"sigma_NI(E) of MF=3 times P_NI(mu, E) of MF=6, each interpolated "
               f"linearly (the manual's reading, p. 146), gives the elastic cross section to {self.reached:.1e} by the "
               f"grid's rule at every energy the engine computed outside the brackets below, the quarter points "
               f"and the midpoint of every written interval among them"
               + ("." if self.checked else "; the search for energies stopped at its budget or its limit before "
                                            "every interval met the tolerance.")]
        if cut:
            out.append(f"Each energy holds at most {self.limits[1]} cosines ({words} {self.limits[1]}); between them "
                       f"sigma_e - sigma_c is good to {self.cosines_reached:.1e} by the grid's rule (asked: "
                       f"{self.tolerance:g}).")
        if self.brackets:
            out.append(f"sigma_NI changes sign {len(self.brackets)} time(s); each change lies between two written "
                       f"energies at most {self.bracket_width:g} x E apart, inside which the manual's reading is "
                       "wrong by up to sigma_e - sigma_c: " + ", ".join(f"{lo * 1e6:.3f}-{hi * 1e6:.3f} eV"
                                                                       for lo, hi in self.brackets) + ".")
        return out


#: The steps of :func:`write`, in order, as ``progress`` names them.
WRITE_STEPS = ("the levels in the formal basis (MF=2)", "the energy grid of MF=3 and MF=6",
               "elastic scattering (MF=6/MT=2)", "the covariance (MF=32)", "MF=2 against MF=3", "the file")


@dataclasses.dataclass(frozen=True)
class AmplitudeReport:
    """The elastic section in the nuclear amplitude expansion (:mod:`pyazr.endf.amplitudes`): ``energies`` written,
    ``candidates`` the engine computed, ``reached`` (the largest deviation by the grid's rule of the reader's
    interpolation, at the energies computed and the quarter points and midpoint of every written interval),
    ``checked`` (every written interval passed its quarter points), the tolerance asked, the limits.  ``form`` is
    "amplitudes"; ``brackets`` is empty (no sigma_NI, no sign change)."""
    energies: int
    candidates: int
    reached: float
    checked: bool
    tolerance: float
    limits: tuple = (None, None)
    limits_name: str = ""
    form: str = "amplitudes"
    brackets: tuple = ()

    def lines(self):
        """The statement of MF=1 and of the report, in words."""
        limited = self.limits[0] is not None and self.energies >= self.limits[0] and self.reached > self.tolerance
        words = LIMIT_WORDS.get(self.limits_name, "the limit of")
        return [f"Elastic scattering (MF=6/MT=2, LAW=5 LTP=1, the nuclear amplitude expansion of eqs. 6.13-6.14 with "
                f"the engine's coefficients; MF=3/MT=2 holds 1.0) on {self.energies} energies chosen of "
                f"{self.candidates} the engine computed{f' ({words} {self.limits[0]} per MF=6 section)' if limited else ''}"
                f", so that sigma_e - sigma_c interpolated linearly at fixed mu between them, the Coulomb term at the "
                f"energy, gives the elastic cross section to {max(self.reached, 0.0):.1e} by the grid's rule at every "
                f"energy the engine computed, the quarter points and the midpoint of every written interval"
                + ("." if self.checked else "; the search for energies stopped at its budget or its limit before "
                                            "every interval met the tolerance.")]


@dataclasses.dataclass(frozen=True)
class Written:
    """What :func:`write` wrote: the file, its sections (MF, MT, lines), the grid of each MF=3 section (points, lab
    eV range; MF=6 has the same energies), the largest relative difference of MF=2's reconstruction from MF=3 per
    MT, the losses, the sections not written, the cosine grid of LAW=5 per MT (points, first, last), the number
    of engine sessions the grids took, the elastic section's report and ``accuracy``: a line per MF=6 section that
    the limits (``MF6_LIMITS``) or the manual's reading of LAW=5 keep from the tolerance asked, with the tolerance
    reached (the same lines end ``losses``, which the page shows), and ``reached``: {MT: that tolerance} of those
    MF=6 sections.  ``print()`` lists it."""
    path: str
    sections: tuple
    grids: dict
    mf2_check: dict
    losses: tuple
    not_written: tuple
    cosines: dict = dataclasses.field(default_factory=dict)
    sessions: int = 0
    elastic: object = None
    accuracy: tuple = ()
    reached: dict = dataclasses.field(default_factory=dict)

    def __str__(self):
        out = [f"{self.path}"]
        for mf, mt, nc in self.sections:
            g = self.grids.get(mt) if mf == 3 else None
            if mf == 6 and mt in self.cosines:
                c = self.cosines[mt]
                g = None
                nc = f"{nc} lines, LAW=5 on up to {c[0]} cosines {c[1]:g} to {c[2]:g} per energy,"
            extra = f", {g[0]} points {g[1]:.6g}-{g[2]:.6g} eV" if g else ""
            check = f", MF=2 reproduces it to {self.mf2_check[mt]:.1e}" if mf == 3 and mt in self.mf2_check else ""
            out.append(f"  MF={mf}/MT={mt} {names.reaction_name(mt)}: {nc}{'' if isinstance(nc, str) else ' lines'}"
                       f"{extra}{check}")
        out += [f"  lost: {t}" for t in self.losses]
        out += [f"  not written: {t}" for t in self.not_written]
        return "\n".join(out)


def _symbol(z):
    from .. import _azure2
    return _azure2.element_symbols()[int(z)]


def _zsymam(z, a):
    """ZSYMAM (p. 40): Z right-justified in columns 1-3, '-', the symbol in 5-6, '-', A in 8-10."""
    return f"{int(z):3d}-{_symbol(z):<2s}-{int(a):3d}"


def _za(z, mass_neutron):
    """ZA = 1000 Z + A of a particle of the record (charge, mass in neutron masses)."""
    return 1000 * int(round(z)) + int(round(mass_neutron * M_NEUTRON_U))


def _law2_columns(coefficients):
    """p(mu) of a LAW=2 table at the cosines of the energy grid's rule, (n_E, 9)."""
    return np.array([engine_values.law2_probability(a, engine_values.LAW2_COSINES) for a in coefficients])


def _mf6_energies(mt, e_MeV, columns, tolerance, limit, limits_name, accuracy, reached_by_mt):
    """The energies of an MF=6 section of the grid: all of them up to ``limit``; above it the subset of
    :func:`pyazr.endf.grid.thin_to` (the grid's rule on ``columns``, (n_E, k) arrays, to the limit), with a line in
    ``accuracy`` that says the tolerance reached."""
    if len(e_MeV) <= limit:
        return np.arange(len(e_MeV))
    keep, reached = grid.thin_to(e_MeV, np.hstack(columns), tolerance, limit=limit)
    reached_by_mt[mt] = max(reached, tolerance)
    accuracy.append(f"MF=6/MT={mt} ({names.reaction_name(mt)}): {len(keep)} of the {len(e_MeV)} energies of MF=3 "
                    f"({LIMIT_WORDS[limits_name]} {limit}); between them the distribution is good to "
                    f"{max(reached, tolerance):.1e} by the grid's rule (asked: {tolerance:g}), MF=3 to the tolerance.")
    return keep


def not_written_known(record):
    """The sections known not to be written before the engine computes anything (what a preview can say): elastic
    scattering of a neutral projectile.  :func:`write` lists the same lines first in ``Written.not_written``."""
    ent = record.pairs[0]
    if round(ent.za) > 0 and round(ent.zb) > 0:
        return ()
    return tuple(f"MF={mf}/MT=2 ({names.reaction_name(2)}): {NOT_WRITTEN_REASONS[(mf, 'neutral_elastic')]}"
                 for mf in (3, 6))


def _reactions(session, project, record):
    """The reactions the file holds, as :class:`pyazr.endf.engine_values.Reaction`: elastic scattering of a charged
    projectile (MT=2), every particle exit (its MT from the record), capture to every final state (MT=102) and to each
    final state (one per photon pair, for MF=6/MT=102).  Labels and pair keys are the project's."""
    pairs = {p.key: p.to_dict() for p in session.pairs}
    entrance = next(p for p in pairs.values() if p["is_entrance"])
    e = entrance["key"]
    particle = [k for k, p in pairs.items() if p["pair_kind"] == "particle" and k != e]
    mt_of = dict(zip([e] + particle, [p.mt for p in record.pairs]))
    out = []
    if entrance["Z1"] > 0 and entrance["Z2"] > 0:
        out.append(engine_values.Reaction(2, "elastic", project.reaction(e, e),
                                          tuple(_pair_keys(project, e, e).items()), e))
    for k in particle:
        out.append(engine_values.Reaction(mt_of[k], "exit", project.reaction(e, k),
                                          tuple(_pair_keys(project, e, k).items()), k))
    photons = [k for k, p in pairs.items() if p["pair_kind"] == "photon"]
    if photons:
        out.append(engine_values.Reaction(102, "capture", project.reaction(e, -1), (), -1))
        for k in photons:
            out.append(engine_values.Reaction(102, "final_state", project.reaction(e, k),
                                              tuple(_pair_keys(project, e, k).items()), k, pairs[k]["excitation"]))
    return out


def _pair_keys(project, entrance_key, exit_key):
    """The pair keys an extrapolation of entrance -> exit names: a key only where its label names two pairs of the
    same nuclei (the engine's rule, E_PAIR_KEY / E_AMBIGUOUS_PAIR)."""
    pairs = {p.key: (p.light, p.heavy, p.light_excitation, p.heavy_excitation) for p in project.pairs()}
    out = {}
    for field, key in (("entrance_pair", entrance_key), ("exit_pair", exit_key)):
        if key > 0 and sum(v == pairs[key] for v in pairs.values()) > 1:
            out[field] = key
    return out


def _reconstruct(sections, s_entrance, energies_MeV, workdir):
    """MF=2 and MF=3 as written (MF=3 gives the residuals' excitations, QM - QI), read back into a project
    (:func:`pyazr.endf.to_project`) and its reaction cross sections at
    the given lab energies: {MT: values in b}."""
    from . import read
    path = os.path.join(workdir, "mf2.endf")
    _write_material(path, 9999, sections, header_rows=None, tpid="MF=2 check")
    m = read(path)
    project, conversion = to_project(m.resonances(), os.path.join(workdir, "mf2.azr"), separation_energy=s_entrance,
                                     material=m)
    key_of = {p["mt"]: p["key"] for p in conversion.pairs}
    out = {}
    for mt, e in energies_MeV.items():
        np.savetxt(os.path.join(workdir, f"e{mt}.dat"), e, fmt="%.17g")
        exit_key = key_of[mt]
        project.add_extrapolation(f"mt{mt}", reaction=project.reaction(1, exit_key), quantity="sigma", unit="b",
                                  energies=f"file {os.path.join(workdir, f'e{mt}.dat')} MeV lab",
                                  **_pair_keys(project, 1, exit_key))
    project.save()
    s = project.session(extrapolations=True)
    try:
        ev = s.evaluate(s.parameters.x0)
        for mt, d in zip(energies_MeV, ev.data):
            out[mt] = np.asarray(d.model, float)
    finally:
        s.close()
    return out


def _write_material(path, mat, sections, header_rows, tpid):
    """The lines of a tape (p. 32): TPID, the sections in (MF, MT) order with sequence numbers, SEND (NS 99999), FEND,
    MEND, TEND (p. 34).  ``header_rows`` are MF=1/MT=451's rows (None: no MF=1)."""
    allsec = dict(sections)
    if header_rows is not None:
        allsec[(1, 451)] = header_rows
    zero = NUMBER_FORMS["classic"](0.0) * 2 + f"{0:11d}" * 4

    def line(content, m, mf, mt, ns):
        return content.ljust(66) + f"{m:4d}{mf:2d}{mt:3d}{ns:5d}"

    out = [line(tpid[:66], 1, 0, 0, 0)]
    last_mf = None
    for mf, mt in sorted(allsec):
        if last_mf is not None and mf != last_mf:
            out.append(line(zero, mat, 0, 0, 0))
        for k, row in enumerate(allsec[(mf, mt)], start=1):
            out.append(line(row, mat, mf, mt, k))
        out.append(line(zero, mat, mf, 0, 99999))
        last_mf = mf
    out.append(line(zero, mat, 0, 0, 0))
    out.append(line(zero, 0, 0, 0, 0))
    out.append(line(zero, -1, 0, 0, 0))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")


def _header_rows(record, header, sections, text, number_form):
    """MF=1/MT=451 (p. 42): HEAD, three CONT, NWD text lines, NXC directory lines (NC of every section and of this
    one; MOD = NMOD)."""
    ent = record.pairs[0]
    nxc = len(sections) + 1
    nwd = len(text)
    n451 = 4 + nwd + nxc
    w = RecordWriter(number_form)
    za_projectile = 1000 * int(round(ent.za)) + int(round(ent.ma * M_NEUTRON_U))
    # AWI: the projectile's mass in neutron masses, nuclear (p. 24), from AME2020, not the project's mass, which a
    # project may give nuclear or atomic
    awi = section_rows.product_mass(round(ent.za), round(ent.ma * M_NEUTRON_U))
    w.cont(record.za, record.awr, LRP, 0, header["NLIB"], header["NMOD"])
    w.cont(0.0, 0.0, 0, 0, 0, 6)
    w.cont(awi, record.eh, header["LREL"], 0, 10 * za_projectile, header["NVER"])
    w.cont(0.0, 0.0, 0, 0, nwd, nxc)
    for t in text:
        w.text(t)
    for mf, mt in [(1, 451)] + sorted(sections):
        nc = n451 if (mf, mt) == (1, 451) else len(sections[(mf, mt)])
        w.rows.append(" " * 22 + "".join(f"{v:11d}" for v in (mf, mt, nc, header["NMOD"])))
    return w.rows


def _wrap(text, width=66):
    """Text wrapped into TEXT records of at most 66 characters (p. 33)."""
    out, line = [], ""
    for word in text.split():
        if line and len(line) + 1 + len(word) > width:
            out.append(line)
            line = word
        else:
            line = f"{line} {word}" if line else word
    if line:
        out.append(line)
    return out


def write(session, path, x=None, *, header=None, energy_range_MeV=None, tolerance=1e-3, max_points=4000,
          number_form="widest", overwrite=False, elastic_max_cosine=ELASTIC_MAX_COSINE, mf6_limits="checkr",
          project_name=None, elastic_form=ELASTIC_FORM, covariance=None, progress=None):
    """Write ``session``'s model at the fit vector ``x`` (default its start) as an ENDF-6 material at ``path``.

    ``header``: the evaluator's fields of MF=1 (p. 42), a dict with any of ``MAT`` (default 9999), ``NLIB``
    (0), ``NVER`` (1), ``LREL`` (0), ``NMOD`` (1), ``ALAB``, ``EDATE``, ``AUTH``, ``REF``, ``DDATE``, ``RDATE``,
    ``LIBRARY`` (HSUB's first line), ``TEXT`` (lines of free text).  ``energy_range_MeV``: (low, high) lab energy
    of the projectile (default: the data's, as the session has it); EL, EH of MF=2 and EMAX.  ``tolerance`` and
    ``max_points``: of the engine's adaptive grid (relative, for linear interpolation between its points).
    ``number_form``: ``"widest"`` (default), ``"classic"`` or ``"e"`` (:data:`pyazr.endf.NUMBER_FORMS`).
    ``mf6_limits``: ``"checkr"`` (default: at most 200 energies per MF=6 section and 201 cosines per LAW=5 panel)
    or ``"manual"`` (Appendix G's: 2,000 and 10,000), :data:`MF6_LIMITS`.  ``project_name``: the name of the
    project that MF=1's text and the tape's first line give (default: the session's file name).  ``elastic_form``:
    charged-particle elastic scattering as ``"amplitudes"`` (LAW=5 LTP=1) or ``"tabulated"`` (LTP=12 with sigma_NI
    in MF=3/MT=2), :data:`ELASTIC_FORMS`; default :data:`ELASTIC_FORM`.  ``covariance``: None (no MF=32),
    ``"hesse"`` (HESSE in a formal-basis session at ``x``: :func:`pyazr.endf.covariance.formal_covariance`), or a
    :class:`pyazr.endf.covariance.FormalCovariance` of this session at ``x``; written as MF=32/MT=151 LCOMP=1.  A
    covariance of the fit's own basis (``EngineFitResult.covariance``, Brune or Park) or a prior's matrix is refused:
    MF=32 holds the formal parameters of MF=2 (p. 253).  ``progress``: called as ``progress(done, total, step)``
    after each of the :data:`WRITE_STEPS` (``step`` its words; a step the file does not need is passed over and
    counted), from the caller's thread; an exception in it stops the write.

    Returns a :class:`Written`.  Refused (:class:`pyazr.FormatRefused`): an existing ``path`` without ``overwrite``,
    what :func:`pyazr.endf.from_session` refuses, a session without exit reactions, MF=2 that does not reproduce
    MF=3 to ``MF2_TOLERANCE``."""
    path = os.path.abspath(path)
    if os.path.exists(path) and not overwrite:
        raise FormatRefused(what=f"{path} exists", allowed="a new file (an existing file is replaced only when that is "
                            "asked for)", fix="write to a file name that does not exist yet, or ask for the "
                            "existing file to be replaced (overwrite)")
    if mf6_limits not in MF6_LIMITS:
        raise FormatRefused(what=f"mf6_limits {mf6_limits!r}", allowed=", ".join(repr(k) for k in MF6_LIMITS),
                            fix="leave it at \"checkr\" (files every checking code reads) or ask for \"manual\" "
                            "(the format's own limits)")
    energy_limit, cosine_limit = MF6_LIMITS[mf6_limits]
    if not (covariance is None or (isinstance(covariance, str) and covariance == "hesse") or
            isinstance(covariance, FormalCovariance)):
        raise FormatRefused(what=f"covariance of type {type(covariance).__name__}",
                            allowed="None, \"hesse\" or a pyazr.endf.covariance.FormalCovariance",
                            fix="MF=32 holds the covariance of MF=2's formal parameters: ask for \"hesse\" (HESSE in a "
                            "formal-basis session at the fit vector), not the fit's covariance in its own basis")
    if elastic_form not in ELASTIC_FORMS:
        raise FormatRefused(what=f"elastic_form {elastic_form!r}", allowed=", ".join(repr(k) for k in ELASTIC_FORMS),
                            fix="leave it at the default or name one of the two forms")
    shown_name = os.path.basename(session.file) if project_name is None else str(project_name)
    if number_form not in NUMBER_FORMS:
        raise FormatRefused(what=f"number_form {number_form!r}", allowed=", ".join(NUMBER_FORMS),
                            fix="choose the widest form (the default): fixed point where it fits, 8-10 digits")
    h = dict(MAT=9999, NLIB=0, NVER=1, LREL=0, NMOD=1, ALAB="", EDATE="", AUTH="", REF="", DDATE="", RDATE="",
             LIBRARY="AZURE2", TEXT=())
    unknown = set(header or {}) - set(h)
    if unknown:
        raise FormatRefused(what=f"header fields {sorted(unknown)}", allowed=", ".join(h),
                            fix="name the MF=1 fields as listed")
    h.update(header or {})
    if energy_range_MeV is None:
        lab = [e for d in session.evaluate(session.parameters.x0 if x is None else x).data
               for e in np.asarray(d.extra.get("energy_lab", []), float)]
        if not lab:
            raise FormatRefused(what="the session has no data to take an energy range from",
                                allowed="an energy range given with the request (lab MeV of the projectile)",
                                fix="give the lowest and highest lab energy of the file")
        energy_range_MeV = (min(lab), max(lab))
    lo, hi = (float(v) for v in energy_range_MeV)
    record, losses = from_session(session, x, energy_range_eV=(lo * EV_PER_MEV, hi * EV_PER_MEV))

    def step(k):
        if progress is not None:
            progress(k + 1, len(WRITE_STEPS), WRITE_STEPS[k])

    step(0)
    from ..project import Project
    project = Project.open(session.file)
    pairs = {p.key: p.to_dict() for p in session.pairs}
    entrance = next(p for p in pairs.values() if p["is_entrance"])
    s_ent = entrance["sep_energy"] + entrance["excitation"]
    reactions = _reactions(session, project, record)
    ent = record.pairs[0]
    awi = section_rows.product_mass(round(ent.za), round(ent.ma * M_NEUTRON_U))
    identical = round(ent.za) == round(ent.zb) and round(ent.ma * M_NEUTRON_U) == round(ent.mb * M_NEUTRON_U)
    mu_range = (0.0 if identical else -1.0, float(elastic_max_cosine))
    from .curves import CoulombTerms
    coulomb = CoulombTerms.of(round(ent.za), round(ent.zb), awi * M_NEUTRON_U, record.awr * M_NEUTRON_U)

    def coulomb_minus(e_MeV, mu):
        """The Coulomb term of LAW=5 as the file's reader adds it (b/sr), (n_E, n_mu)."""
        return np.array([coulomb.sigma(np.full(len(mu), e * EV_PER_MEV), mu, identical=identical,
                                       two_s=int(round(2 * ent.spin_a))) for e in e_MeV])

    work = tempfile.mkdtemp(prefix="endf_write_")
    order = 2 * max([c.l for g in record.groups for c in g.channels] + [0]) + 4
    sweeps = engine_values.Sweeps(session, x, work, order)
    has_elastic = any(r.role == "elastic" for r in reactions)
    values = sweeps.energy_grid(reactions, grid.start_grid(session, x, (lo, hi)), tolerance, max_points,
                                mu_range if has_elastic else None)
    refined = values.report
    if refined.get("budget_reached"):
        raise FormatRefused(what=f"the grid reached {max_points} points before the tolerance {tolerance:g} "
                            f"(largest deviation {refined.get('largest_deviation'):.2e} at "
                            f"{refined.get('largest_deviation_at')} MeV)", allowed="a grid that meets the tolerance",
                            fix="allow more grid points, a looser tolerance or a narrower energy range")
    e_MeV = values.energies
    e_eV = e_MeV * EV_PER_MEV
    written = {(2, 151): record.rows(number_form)}
    grids, cosines = {}, {}
    za, awr = record.za, record.awr
    capture_total = next((r for r in reactions if r.role == "capture"), None)
    finals = []
    accuracy = []
    reached = {}  # MT of an MF=6 section -> the tolerance it reaches, where the limits thinned it
    for r in reactions:
        if r.role == "exit":
            p = next(q for q in record.pairs if q.mt == r.mt)
            qm = p.q + pairs[r.exit_key]["excitation"] * EV_PER_MEV
            written[(3, r.mt)] = section_rows.cross_section(za, awr, qm, p.q, e_eV, values.sigma[r.key], number_form)
            emitted = (_za(p.za, p.ma), section_rows.product_mass(round(p.za), round(p.ma * M_NEUTRON_U)))
            residual = (_za(p.zb, p.mb), section_rows.product_mass(round(p.zb), round(p.mb * M_NEUTRON_U)))
            a = values.legendre[r.key]
            keep = _mf6_energies(r.mt, e_MeV, [_law2_columns(a)], tolerance, energy_limit, mf6_limits, accuracy,
                                 reached)
            written[(6, r.mt)] = section_rows.two_body(za, awr, emitted, residual, e_eV[keep], a[keep], number_form)
        elif r.role == "capture":
            q = s_ent * EV_PER_MEV
            written[(3, 102)] = section_rows.cross_section(za, awr, q, q, e_eV, values.sigma[r.key], number_form)
        elif r.role == "final_state":
            total = values.sigma[capture_total.key]
            y = np.divide(values.sigma[r.key], total, out=np.zeros(len(total)), where=total != 0)
            finals.append((r.excitation, (s_ent - r.excitation) * EV_PER_MEV, y, values.legendre[r.key]))
        if r.role in ("exit", "capture"):
            grids[r.mt] = (len(e_eV), e_eV[0], e_eV[-1])
    if finals:
        z_c, a_c = round(ent.za + ent.zb), round((ent.ma + ent.mb) * M_NEUTRON_U)
        compound = (1000 * z_c + a_c, section_rows.product_mass(z_c, a_c))
        finals = sorted(finals, key=lambda f: f[0])
        keep = _mf6_energies(102, e_MeV, [np.asarray(f[2])[:, None] for f in finals] +
                             [_law2_columns(f[3]) for f in finals], tolerance, energy_limit, mf6_limits, accuracy,
                             reached)
        written[(6, 102)] = section_rows.capture(za, awr, compound, [(f[1], np.asarray(f[2])[keep], f[3][keep])
                                                                     for f in finals], e_eV[keep], number_form)
    step(1)
    elastic = next((r for r in reactions if r.role == "elastic"), None)
    elastic_report = None
    if elastic is not None and elastic_form == "amplitudes":
        projectile = (_za(ent.za, ent.ma), awi)
        target = (_za(ent.zb, ent.mb), section_rows.product_mass(round(ent.zb), round(ent.mb * M_NEUTRON_U)))
        # the engine's own Coulomb terms (the project's masses): the cross section the rule holds is the engine's
        engine_coulomb = CoulombTerms.of(entrance["Z1"], entrance["Z2"], entrance["M1"], entrance["M2"])
        am = amplitude_section.section(sweeps, elastic, e_MeV, mu_range, engine_coulomb, identical,
                                       int(round(2 * ent.spin_a)), tolerance, max_points, limit=energy_limit)
        el_eV = am.energies * EV_PER_MEV
        written[(3, 2)], written[(6, 2)] = section_rows.elastic_amplitudes(za, awr, projectile, target, identical,
                                                                          ent.spin_a, el_eV, am.records, number_form)
        grids[2] = (len(el_eV), el_eV[0], el_eV[-1])
        elastic_report = AmplitudeReport(len(el_eV), am.candidates, am.reached, am.checked, tolerance,
                                         MF6_LIMITS[mf6_limits], mf6_limits)
        accuracy += elastic_report.lines()
        reached[2] = max(am.reached, tolerance)
    elif elastic is not None:
        el, mu_report = elastic_section.section(sweeps, elastic, e_MeV, mu_range, coulomb_minus, tolerance,
                                                max_points, limit=energy_limit, cosine_limit=cosine_limit)
        if mu_report.get("budget_reached"):
            raise FormatRefused(what=f"the cosine grid of elastic scattering reached {max_points} points before the "
                                f"tolerance {tolerance:g}", allowed="a grid that meets the tolerance",
                                fix="allow more grid points, a looser tolerance, a smaller largest cosine or a "
                                "narrower energy range")
        projectile = (_za(ent.za, ent.ma), awi)
        target = (_za(ent.zb, ent.mb), section_rows.product_mass(round(ent.zb), round(ent.mb * M_NEUTRON_U)))
        el_eV = el.energies * EV_PER_MEV
        panels = el.panels
        written[(3, 2)], written[(6, 2)], _ = section_rows.elastic(za, awr, projectile, target, identical, ent.spin_a,
                                                               el_eV, panels, number_form)
        grids[2] = (len(el_eV), el_eV[0], el_eV[-1])
        elastic_report = ElasticReport(len(el_eV), el.candidates, tuple(el.brackets),
                                       elastic_section.BRACKET_WIDTH, el.reached, el.checked, el.cosines_reached,
                                       tolerance, MF6_LIMITS[mf6_limits], mf6_limits)
        accuracy += elastic_report.lines()
        reached[2] = max(el.reached, el.cosines_reached, tolerance)
        counts = [len(m) for m, _ in panels]
        cosines[2] = (max(counts), float(panels[0][0][0]), float(panels[0][0][-1]))
    step(2)
    covariance_text = []
    if covariance is not None:
        cov = covariance_section.formal_covariance(session, x, record) if isinstance(covariance, str) else covariance
        sizes = {(gi, li): len(g.channels) + 1 for gi, g in enumerate(record.groups) for li in range(len(g.levels))}
        if any((gi, li) not in sizes or col >= sizes[(gi, li)] for gi, li, col in cov.parameters):
            raise FormatRefused(what="the covariance names parameters MF=2 does not have",
                                allowed="a FormalCovariance of this session at this fit vector",
                                fix="ask for covariance=\"hesse\"")
        written[(32, 151)] = covariance_section.rows(record, cov, number_form)
        fixed = [p for p, f in zip(cov.parameters, cov.fixed) if f]
        covariance_text = [
            f"MF=32 (LCOMP=1, one block): the covariance of {len(cov.parameters)} parameters of MF=2 (ER, GAM of "
            f"{sum(len(v) for v in cov.levels.values())} resonances in {len(cov.groups)} spin group(s)): {cov.source} "
            f"(objective {cov.objective_fit:.10g} of the fit, {cov.objective_formal:.10g} of the formal basis); the "
            "norms and other parameters of the fit and the photon channels are marginalised (their rows dropped)"
            + (f"; {len(fixed)} parameter(s) of listed resonances were fixed in the fit and have zero rows: "
               + ", ".join(f"spin group {gi + 1} resonance {li + 1} {'ER' if col == 0 else f'GAM{col}'}"
                           for gi, li, col in fixed) if fixed else "") + "."]
    step(3)
    # MF=2 must reproduce MF=3 for the particle reactions (the photon channels are not in MF=2)
    particle = {r.mt: e_MeV for r in reactions if r.role == "exit"}
    mf2_check = {}
    if particle:
        back = _reconstruct({k: v for k, v in written.items() if k[0] in (2, 3)}, s_ent, particle, work)
        for mt in particle:
            sigma = values.sigma[(mt,)]
            keep = np.abs(sigma) > 1e-12 * np.max(np.abs(sigma)) if np.any(sigma) else np.zeros(len(sigma), bool)
            mf2_check[mt] = float(np.max(np.abs(back[mt][keep] - sigma[keep]) / np.abs(sigma[keep]))) \
                if keep.any() else 0.0
            if mf2_check[mt] > MF2_TOLERANCE:
                raise FormatRefused(what=f"MF=2 as written reproduces MF=3/MT={mt} only to {mf2_check[mt]:.2e}",
                                    allowed=f"agreement to {MF2_TOLERANCE:g} (LRP=2: MF=2 reproduces MF=3)",
                                    fix=f"report the project to the AZURE2 developers (the files are kept in {work})")
    # every section a reaction of the project needs, and why it is not written where it is not
    # (a neutral projectile's elastic scattering is not LAW=5: known before the run, not_written_known)
    needed = [(mf, r.mt, r.role) for r in reactions if r.role != "final_state" for mf in (3, 6)]
    not_written = not_written_known(record) + tuple(
        f"MF={mf}/MT={mt} ({names.reaction_name(mt)}): {NOT_WRITTEN_REASONS[(mf, role)]}"
        for mf, mt, role in needed if (mf, mt) not in written)
    sections = written
    # MF=1
    ent = record.pairs[0]
    target_z, target_a = int(round(ent.zb)), int(round(ent.mb * M_NEUTRON_U))
    projectile = names.SUBLIBRARIES.get(10 * (1000 * int(round(ent.za)) + int(round(ent.ma * M_NEUTRON_U))),
                                        "INCIDENT PARTICLE")
    today = datetime.date.today()
    std = [f"{_zsymam(target_z, target_a):<11s}{h['ALAB'][:11]:<11s}{h['EDATE'][:10]:<10s} {h['AUTH'][:33]}",
           f" {h['REF'][:21]:<21s}{h['DDATE'][:10]:<10s} {h['RDATE'][:10]:<10s}"
           f"{'':12s}{today.strftime('%Y%m%d')}",
           f"---- {h['LIBRARY'][:16]:<16s}  MATERIAL {int(h['MAT']):4d}",
           f"----- {projectile.upper()} DATA",
           "------ ENDF-6"]
    body = [f"Written by pyazr.endf.write (AZURE2) on {today.isoformat()} from the project "
            f"{shown_name} ({'Park' if session.basis == 2 else 'Brune'} fit vector); the "
            "particle channels converted to the formal basis with the engine's boundary values B_c = S_c(E_1) of "
            "each spin group.",
            f"LRP={LRP}: MF=3 holds the cross sections computed by AZURE2 on an adaptive grid (relative tolerance "
            f"{tolerance:g} of linear interpolation, not below 1e-3 of the largest value); MF=2 (LRF=7, KRM=4, "
            "IFG=1) reproduces the particle reactions "
            f"of MF=3 (checked before writing to {MF2_TOLERANCE:g}). MF=2 holds no photon channel: AZURE2 keeps "
            "radiative widths out of the total width; capture is MF=3/MT=102 (internal and external capture).",
            "MF=6 holds the engine's angular distributions in the c.m. (LCT=2), on the energy grid of MF=3 (elastic "
            "scattering on its own, below): LAW=2 "
            "Legendre coefficients of each two-body exit with the residual as LAW=4"
            + ("; charged-particle elastic scattering as LAW=5 LTP=1, the nuclear amplitude expansion (eqs. "
               "6.13-6.14) with the engine's coefficients b_l and a_l at each energy, valid at every angle, and 1.0 "
               "in MF=3/MT=2" if (6, 2) in sections and elastic_form == "amplitudes" else "")
            + (f"; charged-particle elastic scattering as LAW=5 LTP=12, P_NI tabulated on up to {cosines[2][0]} "
               f"cosines from {cosines[2][1]:g} to {cosines[2][2]:g} with sigma_NI in MF=3/MT=2 (below the largest "
               "cosine "
               "the cross section is the Coulomb term alone)" if 2 in cosines else "")
            + ("; for capture (MT=102) the primary photon to each final state (LAW=2, AWP its energy at zero "
               "incident energy, the yield sigma_f/sigma) and its recoil (LAW=4); the cascade below the final "
               "states is not given (AZURE2 does not hold it)" if (6, 102) in sections else "") + ".",
            *covariance_text,
            "Separation energy of the incident pair: "
            f"{s_ent:.9f} MeV (the project's); masses as the project's atomic masses."]
    text = std
    for paragraph in body:
        text = text + [""] + _wrap(paragraph)
    if losses.photon_widths:
        text += ["", "Photon widths of the levels (eV, at the level energy):"]
        for jpi, e, k, r, l, wv in losses.photon_widths:
            text.append(f"  {jpi:5s} E_x {e:11.6f} MeV to pair {k} ({r}{l}): {wv:.6g}")
    if accuracy:
        text += ["", "Accuracy:"]
        for t in accuracy:
            text += _wrap("- " + t)
    text += ["", "Not carried by this file:"]
    for t in losses.items + not_written:
        text += _wrap("- " + t)
    for t in h["TEXT"]:
        text += _wrap(t) or [""]
    step(4)
    header_rows = _header_rows(record, h, sections, text, number_form)
    _write_material(path, int(h["MAT"]), sections, header_rows,
                    tpid=f" pyazr.endf.write, {shown_name}")
    listed_sections = [(1, 451, len(header_rows))] + [(mf, mt, len(rows))
                                                       for (mf, mt), rows in sorted(sections.items())]
    step(5)
    return Written(path, tuple(listed_sections), grids, mf2_check, tuple(losses.items) + tuple(accuracy), not_written,
                   cosines, sweeps.count, elastic_report, tuple(accuracy), reached)
