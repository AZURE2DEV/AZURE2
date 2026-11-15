"""The check of an evaluation's resonance parameters (docs/upgrade/endf.md sections 7.2 and 13.14): what the
checking codes do not do for LRF=7 (they count its records, endf.md 9.1).  Three parts, each a list of findings in
the form of the module's own format check (:class:`pyazr.endf.format_check.Finding`, ``source``
``"pyazr.endf resonances"``), so that a page that shows the own findings shows these too:

- **reconstruction** (code ``RECONSTRUCTION``, ``THRESHOLD``): MF=2 read into a project (:func:`pyazr.endf.to_project`)
  and the engine's cross section of each particle exit at the energies of its MF=3 section, against MF=3: the
  largest and the median relative difference per MT, the points off by more than ``tolerance``; points where MF=3
  exceeds the reconstruction by more than ``THRESHOLD_ORDERS`` orders of magnitude are named as such (they lie
  where the penetrability is small: the first points above a threshold).  With LRP=2 MF=3 is the evaluation and the
  reconstruction a check; with LRP=1 MF=3 holds backgrounds and the comparison is not made (said).
- **masses** (``MASS``, ``Q``, ``Q_MASSES``): each nuclide's mass in MF=2 against AME2020 (atomic, and nuclear for a
  light particle: ENDF-102 p. 72 does not say which MA and MB are), each pair's Q against AME2020's, and against the
  Q the file's own masses give (p. 105: Q from atomic masses).
- **text** (``TEXT``): every flag the MF=1 text names (``LRF=``, ``KRM=``, ``IFG=``, ``LTP=``, ``LRP=``) against the
  records that hold it; a flag whose section is not in the file is said so.
"""

import dataclasses
import os
import re
import tempfile

import numpy as np

from ..errors import FormatRefused
from .format_check import Finding
from .material import M_NEUTRON_U

#: The relative difference of the reconstruction from MF=3 below which a point agrees.
TOLERANCE = 1e-4
#: MF=3 above the reconstruction by more than this many orders of magnitude: a threshold finding.
THRESHOLD_ORDERS = 3.0
#: A mass of the file off AME2020 by more than this (keV) is a finding; a Q value likewise.
MASS_TOLERANCE_KEV = 1.0
UC2_KEV = 931494.10242  # u c^2 in keV (CODATA 2018)
SOURCE = "pyazr.endf resonances"
ELECTRON_MASS_U = 0.000548579909


@dataclasses.dataclass
class ResonanceCheck:
    """The answer of :func:`check_resonances`: the file, ``findings`` (:class:`Finding`, in the order reconstruction,
    masses, text), ``parts`` {part: [rows]} with the numbers (``reconstruction``: per MT the points, the largest and
    median relative difference and where; ``masses``: per nuclide and pair the differences in keV; ``text``: per flag
    the text's value and the record's), and ``notes``.  ``print()`` lists it; ``to_dict()`` is plain data."""
    path: str
    findings: list
    parts: dict
    notes: list

    @property
    def errors(self):
        return [f for f in self.findings if f.severity == "error"]

    def to_dict(self):
        return {"path": self.path, "findings": [dataclasses.asdict(f) for f in self.findings], "parts": self.parts,
                "notes": list(self.notes)}

    def __str__(self):
        counts = {s: sum(f.severity == s for f in self.findings) for s in ("error", "warning", "note")}
        out = [f"{os.path.basename(self.path)}: the resonance parameters checked: {counts['error']} errors, "
               f"{counts['warning']} warnings, {counts['note']} notes"]
        out += [f"  {f}" for f in self.findings]
        out += [f"  note: {n}" for n in self.notes]
        return "\n".join(out)


def check_resonances(material, tolerance=TOLERANCE, separation_energy="ame2020"):
    """Check the resonance parameters of ``material`` (a :class:`pyazr.endf.Material` or a path) against its own
    cross sections, the mass table and its text (the module's text).  Refused (:class:`pyazr.FormatRefused`): a file
    without MF=2 in LRF=7, or one :func:`pyazr.endf.to_project` refuses (the refusal's three parts)."""
    from .material import read
    m = read(material) if isinstance(material, str) else material
    mat = m.mat
    findings, notes = [], []
    parts = {"reconstruction": [], "masses": [], "text": []}
    record = m.resonances()
    work = tempfile.mkdtemp(prefix="endf_check_")
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(work, sub), exist_ok=True)
    from .model import to_project
    project, conversion = to_project(record, os.path.join(work, "evaluation.azr"),
                                     separation_energy=separation_energy, material=m)
    _reconstruction(m, project, conversion, tolerance, findings, parts, notes, work)
    _masses(m, record, conversion, findings, parts)
    _text(m, record, findings, parts)
    findings = [dataclasses.replace(f, mat=mat) for f in findings]
    return ResonanceCheck(m.path, findings, parts, notes)


def _finding(severity, code, mf, mt, message, page):
    return Finding(severity, code, 0, 0, mf, mt, message, page, source=SOURCE)


def _reconstruction(m, project, conversion, tolerance, findings, parts, notes, work):
    """The engine's cross section of each particle exit of MF=2 at MF=3's energies against MF=3 (the module's text)."""
    lrp = m.header.lrp if m.header is not None else None
    if lrp == 1:
        notes.append("LRP=1: MF=3 holds backgrounds to be added to MF=2's cross sections; the reconstruction is not "
                     "compared with it")
        return
    key_of = {p["mt"]: p["key"] for p in conversion.pairs}
    mts = [mt for (mf, mt) in sorted(m.sections) if mf == 3 and mt in key_of and mt != 2]
    if not mts:
        notes.append("no MF=3 section of a particle exit of MF=2: nothing to reconstruct")
        return
    tables = {}
    for mt in mts:
        xs = m.cross_section(mt)
        e = np.asarray(xs.table.x, float) / 1e6
        tables[mt] = (e, np.asarray(xs.table.y, float))
        np.savetxt(os.path.join(work, f"e{mt}.dat"), e, fmt="%.17g")
        project.add_extrapolation(f"mt{mt}", reaction=project.reaction(1, key_of[mt]), quantity="sigma", unit="b",
                                  energies=f"file {os.path.join(work, f'e{mt}.dat')} MeV lab")
    project.save()
    s = project.session(extrapolations=True)
    try:
        data = s.evaluate(s.parameters.x0).data
    finally:
        s.close()
    for mt, d in zip(mts, data):
        e, file_sigma = tables[mt]
        engine = np.asarray(d.model, float)
        ok = (file_sigma != 0) | (engine != 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(ok, np.abs(engine - file_sigma) / np.maximum(np.abs(file_sigma), np.abs(engine)), 0.0)
            orders = np.where((file_sigma > 0) & (engine > 0), np.log10(file_sigma / engine), np.nan)
        off = rel > tolerance
        k = int(np.argmax(rel))
        row = {"mt": mt, "points": int(len(e)), "largest": float(rel[k]), "at_eV": float(e[k] * 1e6),
               "median": float(np.median(rel)), "off": int(off.sum()), "tolerance": tolerance}
        parts["reconstruction"].append(row)
        severity = "warning" if off.any() else "note"
        findings.append(_finding(severity, "RECONSTRUCTION", 3, mt,
                                 f"MF=2 reconstructed by the engine against MF=3/MT={mt} at its {len(e)} energies: "
                                 f"largest relative difference {rel[k]:.2e} at {e[k] * 1e6:.6g} eV, median "
                                 f"{np.median(rel):.1e}; {int(off.sum())} point(s) off by more than {tolerance:g}",
                                 "66, 354-364"))
        high = np.nonzero(orders > THRESHOLD_ORDERS)[0]
        if high.size:
            findings.append(_finding("warning", "THRESHOLD", 3, mt,
                                     f"MF=3/MT={mt} lies {np.nanmax(orders[high]):.0f} orders of magnitude above what "
                                     f"MF=2 gives at {len(high)} point(s), "
                                     + ", ".join(f"{e[i] * 1e6:.6g} eV ({file_sigma[i]:.3g} b against "
                                                 f"{engine[i]:.3g} b)" for i in high[:4])
                                     + ": the file's first points above the threshold do not follow from its MF=2",
                                     "105"))


def _masses(m, record, conversion, findings, parts):
    """The masses and Q values of MF=2 against AME2020 and against each other (the module's text)."""
    from ..sources import masses
    for name, file_u, ame_u in conversion.mass_checks:
        if ame_u is None:
            continue
        d_atomic = (file_u - ame_u) * UC2_KEV
        parts["masses"].append({"nuclide": name, "file_u": file_u, "ame2020_u": ame_u, "difference_keV": d_atomic})
        if abs(d_atomic) > MASS_TOLERANCE_KEV:
            findings.append(_finding("warning", "MASS", 2, 151,
                                     f"the mass of {name} in MF=2 ({file_u:.9f} u) is {d_atomic:+.1f} keV from "
                                     f"AME2020's atomic mass ({ame_u:.9f} u)", "72"))
    ent = record.pairs[0]

    def ame(z, mass_neutron):
        try:
            return masses.atomic_mass((int(round(z)), int(round(mass_neutron * M_NEUTRON_U))))
        except Exception:  # noqa: BLE001 - a nuclide AME2020 does not list: no comparison
            return None

    a_in = [ame(ent.za, ent.ma), ame(ent.zb, ent.mb)]
    for p in record.pairs[1:]:
        if p.mt in (2, 102) or p.mt < 0:
            continue
        a_out = [ame(p.za, p.ma), ame(p.zb, p.mb)]
        q_masses = (ent.ma + ent.mb - p.ma - p.mb) * M_NEUTRON_U * UC2_KEV * 1e3  # eV, the file's masses
        q_ame = None if None in a_in + a_out else (sum(a_in) - sum(a_out)) * UC2_KEV * 1e3
        row = {"mt": p.mt, "q_file_eV": p.q, "q_ame2020_eV": q_ame, "q_file_masses_eV": q_masses}
        parts["masses"].append(row)
        if q_ame is not None:
            d = (p.q - q_ame) / 1e3
            findings.append(_finding("warning" if abs(d) > MASS_TOLERANCE_KEV else "note", "Q", 2, 151,
                                     f"Q of the pair MT={p.mt} in MF=2 ({p.q:.10g} eV) is {d:+.3f} keV from "
                                     f"AME2020's ({q_ame:.10g} eV)", "105"))
        d = (p.q - q_masses) / 1e3
        if abs(d) > MASS_TOLERANCE_KEV:
            findings.append(_finding("warning", "Q_MASSES", 2, 151,
                                     f"Q of the pair MT={p.mt} in MF=2 ({p.q:.10g} eV) and the Q its masses give "
                                     f"({q_masses:.10g} eV) differ by {d:+.1f} keV: the file's masses and Q values "
                                     "do not come from one mass table", "72, 105"))


FLAG = re.compile(r"\b(LRF|KRM|IFG|LTP|LRP)\s*=\s*(-?\d+)")


def _text(m, record, findings, parts):
    """Each flag the MF=1 text names against the record that holds it (the module's text)."""
    text = list(m.header.text) if m.header is not None else []
    held = {"LRF": ({7}, (2, 151)), "KRM": ({record.krm}, (2, 151)), "IFG": ({record.ifg}, (2, 151)),
            "LRP": ({m.header.lrp} if m.header is not None else set(), (1, 451))}
    if (6, 2) in m.sections:
        try:
            prod = m.distribution(2).product(m.header.projectile_za)
            held["LTP"] = ({p.flag for p in prod.panels}, (6, 2))
        except FormatRefused:
            pass
    for number, line in enumerate(text, start=1):
        for flag, value in FLAG.findall(line):
            value = int(value)
            if flag == "LTP" and "LTP" not in held:
                parts["text"].append({"flag": flag, "text": value, "record": None, "line": line.strip()})
                findings.append(_finding("warning", "TEXT", 1, 451,
                                         f"the MF=1 text (line {number}: \"{line.strip()}\") names LTP={value} for "
                                         "elastic scattering; MF=6/MT=2, which would hold it, is not in the file"
                                         + (" (the directory lists it)" if _in_directory(m, 6, 2) else ""),
                                         "146-147"))
                continue
            if flag not in held:
                continue
            values, (mf, mt) = held[flag]
            parts["text"].append({"flag": flag, "text": value, "record": sorted(values), "line": line.strip()})
            if value not in values:
                findings.append(_finding("warning", "TEXT", mf, mt,
                                         f"the MF=1 text (line {number}: \"{line.strip()}\") names {flag}={value}; "
                                         f"the record of MF={mf}/MT={mt} holds {flag}="
                                         + ", ".join(str(v) for v in sorted(values)), "38-42"))


def _in_directory(m, mf, mt):
    h = m.header
    return h is not None and any((d[0], d[1]) == (mf, mt) for d in h.directory)
