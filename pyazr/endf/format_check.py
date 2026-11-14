"""The module's own check of an ENDF-6 file's format (``pyazr.endf.check_format``): no outside program needed.

What it checks, each finding with the page of ENDF-102 (2023) that states the rule:
- lines: the control numbers MAT, MF, MT are integers (p. 30); nothing beyond column 80 (p. 30: undefined, "we
  recommend not using them"); text characters 32-126 only (p. 33);
- order and ends: MAT, MF, MT increase (p. 32); every section ends with SEND, every file with FEND, the material
  with MEND, the tape with TEND, their six fields zero (p. 34); a section given once;
- sequence numbers NS (p. 33): optional and declared obsolete in 2013 (p. 30); where a file writes them, each section
  counts 1, 2, ... and SEND has 99999, FEND, MEND, TEND 0;
- the directory of MF=1/MT=451 (pp. 40-42): NXC lines; each section of the material listed once with NC its lines
  without SEND; nothing listed that is not there; MOD <= NMOD; NWD text lines;
- the records of MF=1/MT=451, MF=3 (HEAD and one TAB1, pp. 103-104) and MF=6 (pp. 136-147) read field by field: every
  number readable (p. 31), counts NPL, NR, NP, NZ that fit the section (pp. 34-35), no line left after the last
  record; the interpolation laws 1-6 in MF=3 (Table 16, p. 26); energies that do not decrease (p. 25); the HEAD's ZA
  and AWR those of MF=1 (p. 34);
- MF=2/MT=151 with LRF=7 (pp. 71-74): the record counts (the R-matrix itself is the engine's business).

Severity: ``error`` breaks a rule a reader relies on; ``warning`` breaks a rule of the manual a reader can live with;
``note`` is information (a file without sequence numbers).
"""

import dataclasses
import os

from ..errors import Error
from .material import Material, group, is_tpid, read_header
from .records import Cursor, Line, read_float, read_int

ALLOWED_TEXT = set(range(32, 127))


@dataclasses.dataclass(frozen=True)
class Finding:
    """One finding: severity (error, warning, note), a short code, the line (1-based, 0 for the whole file), MAT,
    MF, MT, the text and the page of ENDF-102 (2023) with the rule; ``source`` the checker ("pyazr.endf", or the
    checking code's name); for a code, ``record`` its record or sequence number and ``meaning`` what a known message
    means for an evaluator."""
    severity: str
    code: str
    line: int
    mat: int
    mf: int
    mt: int
    message: str
    page: str
    source: str = "pyazr.endf"
    record: int = 0
    meaning: str = ""

    def __str__(self):
        where = f"line {self.line}" if self.line else "file"
        sec = f" MAT={self.mat} MF={self.mf} MT={self.mt}" if self.mat else ""
        page = f" (ENDF-102 p. {self.page})" if self.page else ""
        return f"{self.severity:7s} {self.code:16s} {where}{sec}: {self.message}{page}"


@dataclasses.dataclass
class FormatReport:
    """The findings of :func:`check_format` (a list, file order) and the file; prints itself."""
    path: str
    findings: list

    @property
    def errors(self):
        return [f for f in self.findings if f.severity == "error"]

    @property
    def ok(self):
        """True when no finding is an error."""
        return not self.errors

    def to_dict(self):
        return {"path": self.path, "ok": self.ok, "findings": [dataclasses.asdict(f) for f in self.findings]}

    def __str__(self):
        counts = {s: sum(f.severity == s for f in self.findings) for s in ("error", "warning", "note")}
        head = (f"{os.path.basename(self.path)}: format check by pyazr.endf: {counts['error']} errors, "
                f"{counts['warning']} warnings, {counts['note']} notes")
        return "\n".join([head] + [f"  {f}" for f in self.findings])


def _lines(path, add):
    """The lines that can be split (a finding for each that cannot), and the raw lengths."""
    out = []
    with open(path, "r", encoding="ascii", errors="replace", newline="") as f:
        raw = f.read().splitlines()
    for k, text in enumerate(raw, start=1):
        if not text.strip():
            add("warning", "BLANK_LINE", k, 0, 0, 0, "an empty line (every line carries MAT, MF, MT)", "30")
            continue
        body = text.rstrip("\r\n")
        if len(body.rstrip()) > 80:
            add("warning", "BEYOND_COL_80", k, 0, 0, 0, f"characters beyond column 80 ({len(body.rstrip())} "
                "columns); columns greater than 80 are undefined", "30")
        bad = sorted({c for c in body[:80] if ord(c) not in ALLOWED_TEXT})
        if bad:
            add("error", "CHARACTER", k, 0, 0, 0, f"characters outside ASCII 32-126 ({', '.join(map(repr, bad))})",
                "33")
        padded = body.ljust(80)
        try:
            mat = read_int(padded[66:70], path, k, 67)
            mf = read_int(padded[70:72], path, k, 71)
            mt = read_int(padded[72:75], path, k, 73)
        except Error as e:
            add("error", "CONTROL", k, 0, 0, 0, f"MAT, MF or MT not an integer: {e.what}", "30")
            continue
        if len(body.rstrip()) < 75:
            add("error", "SHORT_LINE", k, mat, mf, mt, "the line ends before columns 67-75 (MAT, MF, MT)", "30")
        out.append(Line(k, padded[:66], mat, mf, mt, padded[75:80].strip(), len(body)))
    return out


def _check_order_and_ends(lines, add):
    body = lines[1:] if is_tpid(lines[0]) else lines
    prev = None
    kind = None
    for ln in body:
        kind = ("TEND" if ln.mat == -1 else "MEND" if ln.mat == 0 else "FEND" if ln.mf == 0
                else "SEND" if ln.mt == 0 else "DATA")
        if kind != "DATA" and any(_nonzero(ln.content[11 * k:11 * k + 11]) for k in range(6)):
            add("warning", "END_FIELDS", ln.number, ln.mat, ln.mf, ln.mt, f"a {kind} record with non-zero fields",
                "34")
        if prev is not None and prev.mat > 0 and kind == "DATA" and ln.mat == prev.mat:
            if (ln.mf, ln.mt) < (prev.mf, prev.mt) and prev.mt != 0 and prev.mf != 0:
                add("error", "ORDER", ln.number, ln.mat, ln.mf, ln.mt, f"MF={ln.mf}/MT={ln.mt} after "
                    f"MF={prev.mf}/MT={prev.mt}: MF and MT must increase", "32")
            if (ln.mf, ln.mt) != (prev.mf, prev.mt) and prev.mt != 0 and prev.mf != 0 and prev.mat > 0:
                add("error", "NO_SEND", prev.number, prev.mat, prev.mf, prev.mt,
                    f"the section MF={prev.mf}/MT={prev.mt} ends without SEND", "34")
            if prev.mt == 0 and prev.mf != 0 and ln.mf != prev.mf and ln.mt != 0:
                add("error", "NO_FEND", prev.number, prev.mat, prev.mf, 0, f"the file MF={prev.mf} ends without FEND",
                    "34")
        if prev is not None and prev.mat > 0 and ln.mat != prev.mat and kind != "MEND":
            add("error", "NO_MEND", prev.number, prev.mat, prev.mf, prev.mt,
                f"the material MAT={prev.mat} ends without MEND", "34")
        if prev is not None and kind == "MEND" and prev.mf != 0:
            add("error", "NO_FEND", prev.number, prev.mat, prev.mf, prev.mt,
                f"the material ends without FEND after MF={prev.mf}", "34")
        prev = ln
    if kind != "TEND":
        add("error", "NO_TEND", body[-1].number if body else 0, 0, 0, 0, "the file does not end with TEND (MAT=-1)",
            "34")


def _nonzero(field):
    if not field.strip():
        return False
    try:
        return read_float(field) != 0.0
    except Error:
        return True


def _check_sequence(lines, add):
    if not any(ln.ns for ln in lines):
        add("note", "UNSEQUENCED", 0, 0, 0, 0, "no sequence numbers NS (optional, declared obsolete in 2013)", "30")
        return
    count = {}
    for ln in lines:
        try:
            ns = int(ln.ns) if ln.ns else None
        except ValueError:
            add("warning", "SEQUENCE", ln.number, ln.mat, ln.mf, ln.mt, f"NS {ln.ns!r} is not a number", "33")
            continue
        if ln.mat <= 0 or ln.mf == 0:
            want = 0
        elif ln.mt == 0:
            want = 99999
        else:
            key = (ln.mat, ln.mf, ln.mt)
            count[key] = count.get(key, 0) + 1
            want = count[key]
        if ns != want:
            add("warning", "SEQUENCE", ln.number, ln.mat, ln.mf, ln.mt, f"NS {ln.ns or 'blank'} where {want} is "
                "expected (NS restarts at 1 in each section; SEND 99999; FEND, MEND, TEND 0)", "33")


def _check_records(material, add):
    path = material.path
    h = None
    if (1, 451) in material.sections:
        try:
            h = read_header(path, material.sections[(1, 451)])
        except Error as e:
            add("error", "RECORD", _line_of(e), material.mat, 1, 451, e.what or str(e), _page_of(e, "42"))
    else:
        add("error", "NO_HEADER", 0, material.mat, 1, 451, "the material has no MF=1/MT=451", "38")
    if h is not None and h.lrp == 2:
        add("warning", "LRP", material.sections[(1, 451)].first_line, material.mat, 1, 451, "LRP=2 (File 3 already "
            "holds the cross sections of File 2) is for derived files: not allowed in primary evaluations", "42")
    if h is not None and h.lrp not in (-1, 0, 1, 2):
        add("error", "LRP", material.sections[(1, 451)].first_line, material.mat, 1, 451, f"LRP={h.lrp} is not "
            "defined (-1, 0, 1, 2)", "38-39")
    for (mf, mt), s in sorted(material.sections.items()):
        c = Cursor(path, s.lines)
        try:
            if (mf, mt) == (1, 451):
                read_header(path, s)
                c.i = len(s.lines) if h is None else 4 + h.nwd + h.nxc
            elif mf == 3:
                head = c.cont()
                tab = c.tab1()
                _check_head(material, h, head, mf, mt, add)
                _check_table(tab, mt, mf, material.mat, add, laws=range(1, 7))
                _limit(add, material, mf, mt, tab.line, "NP", len(tab.x), 50000)
                _limit(add, material, mf, mt, tab.line, "NR", len(tab.nbt), 20)
            elif mf == 6:
                head = c.cont()
                _check_head(material, h, head, mf, mt, add)
                _limit(add, material, mf, mt, head.line, "NK", head.n1, 2000)
                for _ in range(max(head.n1, 0)):
                    y = c.tab1()
                    law = y.l2
                    if law in (1, 2, 5):
                        t2 = c.tab2()
                        _limit(add, material, mf, mt, t2.line, "NE", t2.nz, 2000)
                        for _ in range(t2.nz):
                            rec = c.list()
                            if law in (2, 5) and rec.l1 in (0, 1, 2):
                                _limit(add, material, mf, mt, rec.line, "NL", rec.n2, 64)
                    elif law == 6:
                        c.cont()
                    elif law == 7:
                        t2 = c.tab2()
                        for _ in range(t2.nz):
                            inner = c.tab2()
                            for _ in range(inner.nz):
                                c.tab1()
                    elif law not in (0, 3, 4) and law >= 0:
                        add("error", "LAW", y.line, material.mat, mf, mt, f"LAW={law} is not a law of MF=6", "135")
                        c.i = len(s.lines)
            elif (mf, mt) == (2, 151):
                _read_mf2(c, material, add)
            else:
                continue  # sections whose records are not read here
        except Error as e:
            add("error", "RECORD", _line_of(e), material.mat, mf, mt, e.what or str(e), _page_of(e, "34"))
            continue
        for number, column in c.trailing:
            add("warning", "COUNT", number, material.mat, mf, mt, f"fields from column {column} on are not blank after "
                "the last value of a record and not zero: a count (NPL, NR, NP) smaller than the values written?", "34-35")
        if not c.done:
            add("error", "COUNT", s.lines[c.i].number, material.mat, mf, mt, f"{len(s.lines) - c.i} lines after the "
                "last record of the section (a count of a record is smaller than its lines)", "34")
        if s.send is None:
            add("error", "NO_SEND", s.lines[-1].number if s.lines else 0, material.mat, mf, mt,
                "the section has no SEND record", "34")
    return h


def _read_mf2(c, material, add):
    """MF=2/MT=151: HEAD, per isotope a CONT, per range a CONT; LRF=7 records counted (pp. 63-64, 71-74); other
    formats are not read (the rest of the section is skipped)."""
    head = c.cont()
    for _ in range(max(head.n1, 0)):
        iso = c.cont()
        for _ in range(max(iso.n1, 0)):
            rng = c.cont()
            if rng.l1 == 1 and rng.l2 == 7:
                if rng.n1:
                    c.tab1()
                rml = c.cont()
                c.list()
                for _ in range(max(rml.n1, 0)):
                    group_ = c.list()
                    c.list()
                    if group_.l1:  # KBK > 0: the background R-matrix records (pp. 75-76)
                        add("note", "MF2_NOT_READ", group_.line, material.mat, 2, 151,
                            "a background R-matrix (KBK > 0): its records are not counted here", "75")
                        c.i = len(c.lines)
                        return
                    if group_.l2:
                        add("note", "MF2_NOT_READ", group_.line, material.mat, 2, 151,
                            "tabulated phase shifts (KPS > 0): their records are not counted here", "76")
                        c.i = len(c.lines)
                        return
            else:
                add("note", "MF2_NOT_READ", rng.line, material.mat, 2, 151,
                    f"LRU={rng.l1} LRF={rng.l2}: only LRF=7 records are counted here", "63")
                c.i = len(c.lines)
                return


def _limit(add, material, mf, mt, line, name, value, limit):
    """Appendix G (p. 391): the limits a primary evaluation must observe (p. 35; derived libraries may exceed them)."""
    if value > limit:
        add("warning", "LIMIT", line, material.mat, mf, mt, f"{name} = {value} above the limit {limit} of Appendix G "
            "(strict for primary evaluations, relaxed for derived files, p. 35)", "391")


def _check_head(material, h, head, mf, mt, add):
    if h is not None and (abs(head.c1 - h.za) > 1e-6 * max(abs(h.za), 1) or
                          abs(head.c2 - h.awr) > 1e-5 * max(abs(h.awr), 1)):
        add("warning", "HEAD", head.line, material.mat, mf, mt, f"HEAD ZA, AWR = {head.c1:g}, {head.c2:g} where MF=1 "
            f"has {h.za:g}, {h.awr:g}", "34")


def _check_table(tab, mt, mf, mat, add, laws):
    for law in tab.interpolation:
        if law not in laws:
            add("error", "INT", tab.line, mat, mf, mt, f"interpolation law INT={law} in a one-dimensional table",
                "26")
    for a, b in zip(tab.x, tab.x[1:]):
        if b < a:
            add("error", "ENERGY_ORDER", tab.line, mat, mf, mt, f"energies decrease ({a:g} then {b:g} eV)", "25")
            break


def _check_directory(material, h, add):
    if h is None:
        return
    listed = {}
    for mf, mt, nc, mod in h.directory:
        if (mf, mt) in listed:
            add("error", "DIRECTORY", 0, material.mat, 1, 451, f"MF={mf}/MT={mt} listed twice in the directory", "41")
        listed[(mf, mt)] = nc
        if mod > h.nmod:
            add("warning", "DIRECTORY", 0, material.mat, 1, 451, f"MOD={mod} of MF={mf}/MT={mt} above NMOD={h.nmod}",
                "41")
    for key, nc in listed.items():
        s = material.sections.get(key)
        if s is None:
            add("error", "DIRECTORY", 0, material.mat, 1, 451, f"the directory lists MF={key[0]}/MT={key[1]}, which "
                "the material does not hold", "41")
        elif len(s.lines) != nc:
            add("error", "DIRECTORY", s.first_line, material.mat, key[0], key[1], f"the directory gives NC={nc} lines "
                f"for MF={key[0]}/MT={key[1]}, the section has {len(s.lines)} (without SEND)", "41")
    for key in material.sections:
        if key not in listed:
            add("error", "DIRECTORY", material.sections[key].first_line, material.mat, key[0], key[1],
                f"MF={key[0]}/MT={key[1]} is not in the directory of MF=1/MT=451", "41")


def _line_of(e):
    d = getattr(e, "diagnostics", None)
    return d[0].get("line", 0) if d else 0


def _page_of(e, default):
    text = str(e)
    k = text.rfind("ENDF-102 p. ")
    return text[k + 12:].rstrip(")") if k >= 0 else default


def check_format(path):
    """The format of an ENDF-6 file checked by this module (no outside program): a :class:`FormatReport` whose
    findings each name the line, the section and the page of ENDF-102 (2023) with the rule (see the module's text)."""
    findings = []

    def add(severity, code, line, mat, mf, mt, message, page):
        findings.append(Finding(severity, code, line, mat, mf, mt, message, page))

    if not os.path.isfile(path):
        add("error", "NO_FILE", 0, 0, 0, 0, f"{path}: no such file", "30")
        return FormatReport(path, findings)
    lines = _lines(path, add)
    if not lines:
        add("error", "EMPTY", 0, 0, 0, 0, "no line with MAT, MF, MT", "30")
        return FormatReport(path, findings)
    _check_order_and_ends(lines, add)
    _check_sequence(lines, add)
    materials = group(lines)
    seen = {}
    for ln in lines:
        if ln.mat > 0 and ln.mf > 0 and ln.mt > 0:
            key = (ln.mat, ln.mf, ln.mt)
            if key in seen and seen[key] != ln.number - 1 and _starts_again(lines, seen[key], ln.number):
                add("error", "TWICE", ln.number, ln.mat, ln.mf, ln.mt, "the section appears a second time", "32")
            seen[key] = ln.number
    for mat, sections in materials.items():
        m = Material(path, mat, sections)
        h = _check_records(m, add)
        _check_directory(m, h, add)
    findings.sort(key=lambda f: (f.line, f.code))
    return FormatReport(path, findings)


def _starts_again(lines, previous, current):
    """True when another section's line lies between two lines of one section."""
    between = [ln for ln in lines if previous < ln.number < current]
    return any(ln.mt != 0 for ln in between)
