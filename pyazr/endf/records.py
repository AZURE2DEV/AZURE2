"""ENDF-6 lines, numbers and records (ENDF-102, 2023 edition, cited by printed page).

A line (p. 30): columns 1-66 the content, 67-70 MAT, 71-72 MF, 73-75 MT, 76-80 the optional sequence number NS;
columns beyond 80 are undefined.  The content is six 11-column fields (p. 30), or 66 characters of text in a TEXT
record (p. 33).  Numbers (p. 31): an integer field is an optional sign and digits, blank is zero; a float is an
optional sign, digits with an optional decimal point, and an optional exponent written either as a sign and digits
or as ``e``/``E``, an optional sign and digits; spaces inside a number are ignored; ``d``/``D`` is not allowed.
Records (pp. 32-36): TEXT, CONT (with its special cases HEAD, SEND, FEND, MEND, TEND, DIR), LIST, TAB1, TAB2, INTG.

The three number forms a writer may use, all read here (p. 31): the classic one without ``E`` (``-1.234567+6``, 7
significant digits), the form with ``E`` (``-1.23456E+6``, 6 digits) and the widest one, fixed point where it fits
(``-9546287.6``, up to 10 digits), as the ENDF/B-VIII.1 3He + 4He file has.

A malformed field or record raises :class:`pyazr.FormatError` with one diagnostic: the file, the line (1-based), the
columns, the text found, what is allowed and the manual's page.
"""

import dataclasses
import re

from ..errors import FormatError

#: A float field after its spaces are removed (p. 31): mantissa, then an exponent as [eE][sign]digits or sign digits.
_FLOAT = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)(([eE][+-]?\d+)|([+-]\d+))?$")
_INT = re.compile(r"^[+-]?\d+$")


@dataclasses.dataclass(frozen=True)
class Line:
    """One line of a file: its number in the file (1-based), the content (columns 1-66, padded), MAT, MF, MT, the
    sequence-number field as written (columns 76-80, ``""`` when blank) and the raw length of the line."""
    number: int
    content: str
    mat: int
    mf: int
    mt: int
    ns: str
    length: int


def refuse(path, line, columns, found, message, allowed, page, code="ENDF_FIELD"):
    """A :class:`FormatError` for a field or record: where (file:line:columns), what was found, what is allowed and
    the page of ENDF-102 that says so."""
    where = f"{path}:{line}:{columns[0]}-{columns[1]}" if columns else f"{path}:{line}"
    fix = f"correct line {line} of {path} (ENDF-102 p. {page})"
    diagnostic = {"severity": "error", "code": code, "line": line, "record": "", "field": "",
                  "found": found, "allowed": allowed, "message": f"{where}: {message}", "fix": fix}
    return FormatError(f"{where}: {message} (found {found!r}; allowed: {allowed}; ENDF-102 p. {page})",
                       [diagnostic], what=f"{where}: {message}", allowed=allowed, fix=fix)


def split_line(path, number, text):
    """A :class:`Line` from the text of one line (p. 30).  The control fields MAT (I4), MF (I2), MT (I3) must be
    integers; a line shorter than 75 columns cannot hold them and is refused."""
    raw = text.rstrip("\r\n")
    body = raw.rstrip()
    if len(body) < 75 and body.strip():
        raise refuse(path, number, (1, len(body)), body[-20:], "the line ends before its MAT, MF, MT fields "
                     "(columns 67-75)", "80 columns: content 1-66, MAT 67-70, MF 71-72, MT 73-75, NS 76-80", 30,
                     "ENDF_LINE")
    padded = raw.ljust(80)
    ids = []
    for (a, b) in ((67, 70), (71, 72), (73, 75)):
        field = padded[a - 1:b]
        if not _INT.match(field.replace(" ", "") or "0"):
            raise refuse(path, number, (a, b), field, "a control number (MAT, MF or MT) is not an integer",
                         "an integer field (I4, I2, I3)", 30, "ENDF_CONTROL")
        ids.append(int(field.replace(" ", "") or "0"))
    return Line(number, padded[:66], ids[0], ids[1], ids[2], padded[75:80].strip(), len(raw))


def read_float(field, path="", line=0, column=1):
    """The value of an 11-column float field (p. 31): blank is zero; the three number forms; ``d``/``D`` refused."""
    t = field.replace(" ", "")
    if not t:
        return 0.0
    if not _FLOAT.match(t):
        raise refuse(path, line, (column, column + 10), field, "not an ENDF number",
                     "an optional sign, digits with an optional decimal point, an optional exponent as e/E with "
                     "digits or as a sign with digits (no d or D)", 31, "ENDF_NUMBER")
    if "e" in t or "E" in t:
        return float(t)
    for i in range(len(t) - 1, 0, -1):
        if t[i] in "+-":
            return float(t[:i] + "e" + t[i:])
    return float(t)


def read_int(field, path="", line=0, column=1):
    """The value of an 11-column integer field (p. 31): blank is zero."""
    t = field.replace(" ", "")
    if not t:
        return 0
    if not _INT.match(t):
        raise refuse(path, line, (column, column + 10), field, "not an integer", "an optional sign and digits "
                     "(blank is zero)", 31, "ENDF_INTEGER")
    return int(t)


# ---------------------------------------------------------------------------------------------------- number forms
def form_classic(x):
    """The classic form without ``E`` (``-1.234567+6``): 7 significant digits, 6 with a two-digit exponent."""
    if x == 0.0:
        return " 0.000000+0"
    for nd in (6, 5, 4, 3):
        m, e = f"{x:.{nd}e}".split("e")
        s = f"{m}{'+' if int(e) >= 0 else '-'}{abs(int(e))}"
        if len(s) <= 10 or (len(s) == 11 and s[0] == "-"):
            return s.rjust(11)
    raise ValueError(f"{x!r} does not fit an 11-column ENDF field")


def form_e(x):
    """The form with ``E`` (``-1.23456E+6``): 6 significant digits, 5 with a two-digit exponent."""
    if x == 0.0:
        return " 0.00000E+0"
    for nd in (5, 4, 3, 2):
        m, e = f"{x:.{nd}e}".split("e")
        s = f"{m}E{'+' if int(e) >= 0 else '-'}{abs(int(e))}"
        if len(s) <= 10 or (len(s) == 11 and s[0] == "-"):
            return s.rjust(11)
    raise ValueError(f"{x!r} does not fit an 11-column ENDF field")


def form_widest(x):
    """The most digits 11 columns hold: fixed point where it fits, else the exponent without ``E``; whichever of the
    two is closer to ``x`` (the leading ``+`` is optional, p. 31, so a positive number may use column 1)."""
    if x == 0.0:
        return " 0.00000000"
    candidates = []
    for nd in range(10, -1, -1):
        s = f"{x:.{nd}f}"
        if len(s) <= 11 and float(s) != 0.0:
            candidates.append(s)
            break
    for nd in range(9, 0, -1):
        m, e = f"{x:.{nd}e}".split("e")
        s = f"{m}{'+' if int(e) >= 0 else '-'}{abs(int(e))}"
        if len(s) <= 11:
            candidates.append(s)
            break
    if not candidates:
        raise ValueError(f"{x!r} does not fit an 11-column ENDF field")
    return min(candidates, key=lambda s: abs(read_float(s) - x)).rjust(11)


#: The number forms by name (a writer's choice; the reader reads all three).
NUMBER_FORMS = {"classic": form_classic, "e": form_e, "widest": form_widest}


def significant_digits(field):
    """How many significant digits a float field as written carries (for the digits a rewrite loses)."""
    t = field.replace(" ", "").lstrip("+-")
    mantissa = re.split(r"[eE]|(?<=\d)[+-]", t, maxsplit=1)[0].replace(".", "").lstrip("0")
    return max(len(mantissa), 1)


# ---------------------------------------------------------------------------------------------------- records
@dataclasses.dataclass(frozen=True)
class Cont:
    """A CONT record (pp. 33-34) or its HEAD case (C1 = ZA, C2 = AWR): six fields; ``line`` its line number."""
    c1: float
    c2: float
    l1: int
    l2: int
    n1: int
    n2: int
    line: int


@dataclasses.dataclass(frozen=True)
class ListRecord:
    """A LIST record (p. 34): the CONT fields (N1 = NPL, the count of values) and the NPL values B(n)."""
    c1: float
    c2: float
    l1: int
    l2: int
    n1: int
    n2: int
    values: tuple
    line: int


@dataclasses.dataclass(frozen=True)
class Tab1:
    """A TAB1 record (p. 35): the CONT fields, the interpolation ranges (NBT, INT; NR of them) and NP pairs x, y."""
    c1: float
    c2: float
    l1: int
    l2: int
    nbt: tuple
    interpolation: tuple
    x: tuple
    y: tuple
    line: int


@dataclasses.dataclass(frozen=True)
class Tab2:
    """A TAB2 record (p. 35): the CONT fields, the interpolation ranges between the NZ values of z (NBT, INT)."""
    c1: float
    c2: float
    l1: int
    l2: int
    nbt: tuple
    interpolation: tuple
    nz: int
    line: int


class Cursor:
    """Reads the records of one section in order (the lines of the section without its SEND)."""

    def __init__(self, path, lines):
        self.path, self.lines, self.i = path, list(lines), 0
        #: (line number, first column) of the last line of a LIST, TAB1 or TAB2 whose fields after the record's
        #: last value are neither blank nor zero: a count smaller than what was written, which the reader cannot see
        #: otherwise
        self.trailing = []

    @property
    def done(self):
        return self.i >= len(self.lines)

    def _next(self, what):
        if self.i >= len(self.lines):
            last = self.lines[-1].number if self.lines else 0
            raise refuse(self.path, last, None, "end of section", f"the section ends inside a {what} record (a count "
                         "of the record asks for more lines than the section holds)", "counts (NPL, NR, NP, NZ) that "
                         "match the lines of the record", 34, "ENDF_COUNT")
        line = self.lines[self.i]
        self.i += 1
        return line

    def _floats(self, line, n=6):
        return [read_float(line.content[11 * k:11 * k + 11], self.path, line.number, 11 * k + 1) for k in range(n)]

    def _ints(self, line, first=0, n=6):
        return [read_int(line.content[11 * k:11 * k + 11], self.path, line.number, 11 * k + 1)
                for k in range(first, first + n)]

    def text(self):
        """A TEXT record (p. 33): the 66 characters of content."""
        return self._next("TEXT").content

    def cont(self):
        """A CONT (or HEAD) record (pp. 33-34)."""
        line = self._next("CONT")
        c1, c2 = self._floats(line, 2)
        return Cont(c1, c2, *self._ints(line, 2, 4), line.number)

    def _values(self, n, what):
        out = []
        while len(out) < n:
            line = self._next(what)
            k = min(6, n - len(out))
            out.extend(self._floats(line, k))
            self._note_trailing(line, k)
        return tuple(out)

    def _note_trailing(self, line, used):
        """Fields after the last value hold blanks or zeros (some writers pad with 0); anything else is noted."""
        for k in range(used, 6):
            field = line.content[11 * k:11 * k + 11]
            try:
                zero = read_float(field) == 0.0
            except FormatError:
                zero = False
            if not zero:
                self.trailing.append((line.number, 11 * k + 1))
                return

    def list(self):
        """A LIST record (p. 34): NPL = N1 values, six per line."""
        head = self.cont()
        if head.n1 < 0:
            raise refuse(self.path, head.line, (45, 55), str(head.n1), "a negative count NPL", "NPL >= 0", 34,
                         "ENDF_COUNT")
        return ListRecord(head.c1, head.c2, head.l1, head.l2, head.n1, head.n2, self._values(head.n1, "LIST"),
                          head.line)

    def _ranges(self, nr, what, line):
        nbt, law = [], []
        while len(nbt) < nr:
            row = self._next(what)
            k = min(3, nr - len(nbt))
            ints = self._ints(row, 0, 2 * k)
            self._note_trailing(row, 2 * k)
            nbt.extend(ints[0::2])
            law.extend(ints[1::2])
        return tuple(nbt), tuple(law)

    def tab1(self):
        """A TAB1 record (p. 35): NR interpolation ranges, then NP pairs (x, y)."""
        head = self.cont()
        nr, np_ = head.n1, head.n2
        if nr < 0 or np_ < 0:
            raise refuse(self.path, head.line, (45, 66), f"NR={nr} NP={np_}", "a negative count", "NR, NP >= 0", 35,
                         "ENDF_COUNT")
        nbt, law = self._ranges(nr, "TAB1", head.line)
        _check_ranges(self.path, head.line, nbt, law, np_)
        xy = self._values(2 * np_, "TAB1")
        return Tab1(head.c1, head.c2, head.l1, head.l2, nbt, law, xy[0::2], xy[1::2], head.line)

    def tab2(self):
        """A TAB2 record (p. 35): NR interpolation ranges over NZ values of z (the records that follow)."""
        head = self.cont()
        nbt, law = self._ranges(head.n1, "TAB2", head.line)
        _check_ranges(self.path, head.line, nbt, law, head.n2)
        return Tab2(head.c1, head.c2, head.l1, head.l2, nbt, law, head.n2, head.line)


def _check_ranges(path, line, nbt, law, n):
    """The interpolation table (p. 25): NBT increasing, the last one equal to the number of points."""
    if n and nbt and (list(nbt) != sorted(nbt) or nbt[-1] != n):
        raise refuse(path, line, None, f"NBT={list(nbt)} with {n} points", "the interpolation ranges do not end at "
                     "the last point", "NBT(m) increasing and NBT(NR) equal to the number of points", 25,
                     "ENDF_COUNT")
    if n and not nbt:
        raise refuse(path, line, None, "NR=0", "no interpolation range for a table with points", "NR >= 1", 25,
                     "ENDF_COUNT")


