"""An ENDF-6 file read into materials and sections (ENDF-102, 2023, cited by printed page).

A tape (p. 32): an optional first TEXT record TPID (MF = MT = 0, p. 33), then materials, each a sequence of files
(MF) of sections (MT), every line carrying MAT, MF, MT; a section ends with SEND (MT = 0), a file with FEND (MF = 0),
a material with MEND (MAT = 0), the tape with TEND (MAT = -1) (p. 34).  MAT, MF, MT increase (p. 32).
MF=1/MT=451 holds the header and the directory (pp. 38-42).

``read(path)`` reads the structure of every line (the control numbers) and groups the sections; the records of a
section are read when a call asks for them.  What breaks the grouping (a control field that is not a number, a line
shorter than 75 columns) raises :class:`pyazr.FormatError`; everything else a file may get wrong is reported by
``pyazr.endf.check_format`` without refusing to read.
"""

import dataclasses
import os

from ..errors import FormatRefused, MissingFile, listed
from . import names
from .records import Cursor, split_line

M_NEUTRON_U = 1.00866491595  # the neutron mass in u, ENDF-102 Appendix H, Table H.2, p. 394


@dataclasses.dataclass(frozen=True)
class Section:
    """One section of a material: MF, MT, its lines without the SEND record, the SEND line (None when missing)."""
    mf: int
    mt: int
    lines: tuple
    send: object

    @property
    def first_line(self):
        return self.lines[0].number if self.lines else (self.send.number if self.send else 0)

    def cursor(self, path):
        return Cursor(path, self.lines)


@dataclasses.dataclass(frozen=True)
class Header:
    """MF=1/MT=451 (p. 42): the material's numbers, its three standard text lines and the rest of its text, and
    the directory, one (MF, MT, NC, MOD) per section (NC: the lines of the section without SEND, p. 41)."""
    za: float
    awr: float
    lrp: int
    lfi: int
    nlib: int
    nmod: int
    elis: float
    sta: float
    lis: int
    liso: int
    nfor: int
    awi: float
    emax: float
    lrel: int
    nsub: int
    nver: int
    temp: float
    ldrv: int
    nwd: int
    nxc: int
    text: tuple
    directory: tuple

    @property
    def projectile_za(self):
        """The incident particle's ZA: NSUB = 10 x IPART + ITYPE (Table 3, p. 7)."""
        return self.nsub // 10

    @property
    def standard_lines(self):
        """The fields of the three standard text lines (p. 40-41), stripped."""
        t = [line.ljust(66) for line in self.text] + [""] * 3
        return {"ZSYMAM": t[0][0:11].strip(), "ALAB": t[0][11:22].strip(), "EDATE": t[0][22:32].strip(),
                "AUTH": t[0][33:66].strip(), "REF": t[1][1:22].strip(), "DDATE": t[1][22:32].strip(),
                "RDATE": t[1][33:43].strip(), "ENDATE": t[1][55:63].strip(), "HSUB": t[2].strip()}


class Material:
    """One material of an ENDF-6 file: ``mat``, its sections by ``(MF, MT)``, its header (``header``, None
    without MF=1/MT=451) and the calls that read them: :meth:`summary`, :meth:`cross_section` (MF=3),
    :meth:`distribution` (MF=6), :meth:`curves` (the file's values at a data set's points)."""

    def __init__(self, path, mat, sections, tpid=None):
        self.path, self.mat, self.sections, self.tpid = path, mat, dict(sections), tpid
        self._header = None

    def __repr__(self):
        h = self.header
        what = f"{h.standard_lines['ZSYMAM']}, {names.SUBLIBRARIES.get(h.nsub, f'NSUB={h.nsub}')}" if h else "no header"
        return f"Material(MAT={self.mat}, {what}, {len(self.sections)} sections, {os.path.basename(self.path)})"

    def section(self, mf, mt):
        """The section (MF, MT); :class:`pyazr.FormatRefused` when the material has none."""
        if (mf, mt) not in self.sections:
            raise FormatRefused(what=f"{os.path.basename(self.path)} MAT={self.mat} has no section MF={mf}/MT={mt}",
                                allowed="its sections: " + listed(f"MF={a}/MT={b}" for a, b in sorted(self.sections)),
                                fix="Material.summary() lists what the file holds")
        return self.sections[(mf, mt)]

    @property
    def header(self):
        """MF=1/MT=451 as a :class:`Header` (None when the material has no such section)."""
        if self._header is None and (1, 451) in self.sections:
            self._header = read_header(self.path, self.sections[(1, 451)])
        return self._header

    def require_header(self):
        if self.header is None:
            raise FormatRefused(what=f"{os.path.basename(self.path)} MAT={self.mat} has no MF=1/MT=451",
                                allowed="a material that starts with its header and directory (p. 38)",
                                fix="pyazr.endf.check_format(path) lists what the file holds")
        return self.header

    def cross_section(self, mt):
        """MF=3/MT as a :class:`pyazr.endf.distributions.CrossSection` (lab energies in eV, barns)."""
        from .distributions import read_cross_section
        return read_cross_section(self, mt)

    def distribution(self, mt):
        """MF=6/MT as a :class:`pyazr.endf.distributions.Distribution` (its products and their laws)."""
        from .distributions import read_distribution
        return read_distribution(self, mt)

    def summary(self):
        """What the file holds, as a :class:`pyazr.endf.summary.Summary` that prints itself."""
        from .summary import summarize
        return summarize(self)

    def curves(self, request, session=None):
        """The file's values at a data set's points, in its quantity, frame and units (:func:`pyazr.endf.curves`)."""
        from .curves import curves
        return curves(self, request, session=session)


def read_header(path, section):
    """MF=1/MT=451 (p. 42): HEAD, three CONT, NWD TEXT lines, NXC directory CONT lines."""
    c = section.cursor(path)
    h, c2, c3, c4 = c.cont(), c.cont(), c.cont(), c.cont()
    text = tuple(c.text() for _ in range(max(c4.n1, 0)))
    directory = []
    for _ in range(max(c4.n2, 0)):
        d = c.cont()
        directory.append((d.l1, d.l2, d.n1, d.n2))
    return Header(za=h.c1, awr=h.c2, lrp=h.l1, lfi=h.l2, nlib=h.n1, nmod=h.n2, elis=c2.c1, sta=c2.c2, lis=c2.l1,
                  liso=c2.l2, nfor=c2.n2, awi=c3.c1, emax=c3.c2, lrel=c3.l1, nsub=c3.n1, nver=c3.n2, temp=c4.c1,
                  ldrv=c4.l1, nwd=c4.n1, nxc=c4.n2, text=text, directory=tuple(directory))


def scan(path):
    """The lines of a file as :class:`pyazr.endf.records.Line` (raises FormatError where a control number cannot be
    read) and the TPID line (or None)."""
    if not os.path.isfile(path):
        raise MissingFile(path, f"{path}: no such file", what=f"{path}: no such file",
                          allowed="an existing ENDF-6 file", fix=f"ls {os.path.dirname(os.path.abspath(path))}")
    with open(path, "r", encoding="ascii", errors="replace", newline="") as f:
        raw = f.read().splitlines()
    lines = [split_line(path, k + 1, text) for k, text in enumerate(raw) if text.strip()]
    return lines, (lines[0] if lines and is_tpid(lines[0]) else None)


def is_tpid(line):
    """The tape's first line is TPID when its MF and MT are 0 (p. 33: MAT holds the tape number NTAPE, which some
    libraries write as 0) and it is not TEND; a TEXT record, not an end record."""
    return line.mf == 0 and line.mt == 0 and line.mat != -1


def group(lines):
    """Materials as {MAT: {(MF, MT): Section}} in file order, from the lines (TPID, FEND, MEND, TEND skipped)."""
    materials = {}
    current = None
    for line in lines:
        if line.mat <= 0 or (line.mf == 0 and line.mt == 0):
            current = None  # MEND, TEND, FEND, TPID
            continue
        sections = materials.setdefault(line.mat, {})
        if line.mt == 0:  # SEND
            if current is not None and current[0] == line.mat and current[1] == line.mf:
                key = (current[1], current[2])
                sections[key] = dataclasses.replace(sections[key], send=line)
            current = None
            continue
        key = (line.mf, line.mt)
        if current != (line.mat, line.mf, line.mt):
            current = (line.mat, line.mf, line.mt)
            # a section given twice: its lines are joined; check_format reports it
            sections.setdefault(key, Section(line.mf, line.mt, (), None))
        s = sections[key]
        sections[key] = Section(s.mf, s.mt, s.lines + (line,), s.send)
    return materials


def read(path, mat=None):
    """The material of an ENDF-6 file (``mat``: its MAT number, needed when the file holds several).

    >>> m = pyazr.endf.read("h-002_He_004.endf")
    >>> print(m.summary())
    """
    lines, tpid = scan(path)
    materials = group(lines)
    if not materials:
        raise FormatRefused(what=f"{path}: no material (no line with MAT > 0)",
                            allowed="an ENDF-6 file: lines of 80 columns with MAT, MF, MT (p. 30)",
                            fix="pyazr.endf.check_format(path) says what the file holds")
    if mat is None:
        if len(materials) > 1:
            raise FormatRefused(what=f"{path} holds {len(materials)} materials", allowed="mat= one of "
                                + listed(materials), fix=f"pyazr.endf.read({path!r}, mat={next(iter(materials))})")
        mat = next(iter(materials))
    if mat not in materials:
        raise FormatRefused(what=f"{path} has no material MAT={mat}", allowed="mat= one of " + listed(materials),
                            fix=f"pyazr.endf.read({path!r}, mat={next(iter(materials))})")
    return Material(path, mat, materials[mat], tpid=tpid)
