"""The ENDF checking codes, when they are found: CHECKR, FIZCON, PSYCHE, INTER, STANEF (the ENDF utility codes of
the NNDC and the IAEA) and PREPRO's LINEAR, RECENT, FIXUP, DICTIN, ENDF2C (D.E. Cullen, IAEA-NDS-0241).  Nothing
here installs or fetches a program: ``checkers()`` finds what is installed, ``check()`` runs the checking codes and
reads their reports, ``standardize()`` runs a code that rewrites a file and returns the result as a candidate that
replaces nothing until it is accepted.

Where the manuals and the sources disagree about a deck, the source is taken (and the behaviour of the built
programs, run on this machine with the PI's agreement):
- CHECKR, FIZCON, PSYCHE, INTER read on standard input the input file, the message file, "Y" (check the entire
  file) and, after each file, the next input file; a line DONE ends the run (checkr.f lines 255-270 and 806-861).
  Without DONE the program reaches the end of its input and stops with a Fortran runtime error and exit status 2
  after writing a complete report: the exit status does not say whether there are findings; the report does.
- STANEF reads the input file, the output file, the tape label (blank: keep the old one) and one line of options
  (blank: the standard ones; "0,Y,Y,Y": character output, new directory, numbers in standard form, MF=2 too),
  stanef.f lines 145-180 and 536-612.  Its standard options rewrite the numbers of MF=2 too (measured: the manual's
  "MF=2 on request" does not hold for the built 8.10).
- PREPRO's codes read and write fixed names in the working directory: ENDFB.IN, ENDFB.OUT and <CODE>.INP
  (dictin.f lines 205-256, fixup.f lines 587-615; ENDF2C has no input deck, endf2c.f lines 44-46).

Reports are read into findings (:class:`pyazr.endf.format_check.Finding`, ``source`` the program): CHECKR and
FIZCON by their "ERROR(S) FOUND IN MAT=, MF=, MT=" and "WARNING(S) IN MAT=" blocks (checkr.f 7158, 7214) and
"ERROR(S) - MISSING SECTIONS IN MAT" (fizcon.f 10312), counted against their closing "Encountered N errors, M
warnings"; PSYCHE by its "FILE n / SECTION mt" headings, every comment line a warning ("no warnings detected in
sections without comments").  A crash, a time-out or a report that cannot be read is a finding of the run, never an
exception.  Known messages carry their meaning for an evaluator (``Finding.meaning``).
"""

import dataclasses
import os
import re
import shutil
import subprocess
import tempfile

from ..errors import FormatRefused
from .format_check import Finding, check_format
from .records import read_float, significant_digits

#: The programs known: name -> (family, kind, what it does, what it does to the numbers of a file).  kind: "check"
#: (writes a report), "rewrite" (writes a new ENDF file: a candidate), "derive" (writes a derived file, not a
#: correction).  The last column is what an evaluator is told before a run (docs/upgrade/endf.md section 13.2: the
#: built programs as measured on the ENDF/B-VIII.1 3He + 4He file).
PROGRAMS = {
    "checkr": ("utility", "check", "format: every record against the ENDF-6 rules",
               "reads the file and writes a report: no number of the file is changed"),
    "fizcon": ("utility", "check", "procedures and simple physics: energy order, sums, Q values, missing sections",
               "reads the file and writes a report: no number of the file is changed"),
    "psyche": ("utility", "check", "physics: energy conservation, Wick's limit, resonance statistics",
               "reads the file and writes a report: no number of the file is changed"),
    "inter": ("utility", "check", "integral values: thermal, resonance integrals, 14 MeV (linearised MF=3 only)",
              "reads the file and writes a report: no number of the file is changed"),
    "stanef": ("utility", "rewrite", "bookkeeping: tape label, directory, sequence numbers, numbers in standard form",
               "rewrites every number in the classic standard form, 7 significant digits (measured: 2 digits lost "
               "in MF=2, MF=3 and MF=6 of a file written with more), MF=2 included even with the standard options; "
               "a number written without a decimal point (35070343) is misread and written 1e4 too small "
               "(3.507034+3): the code reads with E11.4 (stanef.f line 1131) where the manual says F11.0 (p. 31)"),
    "fixup": ("prepro", "rewrite", "ZA and AWR made consistent, sums redefined, sequence numbers",
              "rewrites the directory and the summed cross sections; measured: 4 numbers of an MF=3 section lost "
              "1 digit (2.5e-6 relative), the numbers of MF=2 and MF=6 written another way with equal values"),
    "dictin": ("prepro", "rewrite", "the directory of MF=1/MT=451 made from the sections",
               "rewrites the directory of MF=1/MT=451 only: no number of another section is changed"),
    "endf2c": ("prepro", "rewrite", "numbers rewritten in a form C and Fortran read alike",
               "rewrites numbers with the letter E; measured: a number with a two-digit exponent loses 1 digit "
               "(-1.38025-12 to -1.3803E-12; up to 3.6e-5 relative in MF=6)"),
    "linear": ("prepro", "derive", "MF=3 linearised",
               "writes a derived file (MF=3 on a linear grid within a tolerance), not a correction of the original"),
    "recent": ("prepro", "derive", "MF=3 reconstructed from the resonance parameters of MF=2",
               "writes a derived file (MF=3 from MF=2), not a correction of the original"),
}
CHECK_DEFAULT = ("checkr", "fizcon", "psyche")

#: Messages of the codes with what they mean for an evaluator of a charged-particle file: (pattern, meaning).
MEANINGS = (
    (r"LRP=2 ONLY ALLOWED IN DERIVED FILES", "ENDF-102 p. 42 reserves LRP=2 for derived files; the ENDF/B-VIII.1 "
     "charged-particle files with LRF=7 use it for a full MF=3 beside MF=2 (a decision of the evaluation)"),
    (r"NE = +\d+ OUT OF RANGE +1 - +200", "the code's own array of 200 incident energies; Appendix G (p. 391) allows "
     "2000 in MF=6: the section is not checked by the code, not wrong"),
    (r"SECTION CANNOT BE CHECKED", "the code skipped the rest of the section after the previous message"),
    (r"IS MISSING", "the directory of MF=1/MT=451 lists a section the material does not hold (p. 41)"),
    (r"REQUIRES AN EQUIVALENT SECTION IN FILE 6", "a charged-particle reaction of MF=3 needs its distribution in MF=6 "
     "(MF=4 is for neutrons and photons only, p. 112)"),
    (r"AP NOT IN RANGE", "PSYCHE reads the record after the range as SPI, AP of LRF 1-3 (psyche.f lines 1428-1434); "
     "in LRF=7 that record holds 0, 0, IFG, KRM, NJS, KRL (p. 73): the number tested is not a radius"),
    (r"FILE IS UNSEQUENCED", "no sequence numbers; optional since 2013 (p. 30)"),
    (r"INVALID ANGULAR REPRESENTATION LTP=15", "CHECKR 8.24 does not know LTP=15, which the 2023 manual defines "
     "(p. 147: the nuclear part P_NI(mu) tabulated, linear in mu): not an error of the file"),
    (r"NW SHOULD BE SET TO \*+", "CHECKR has no rule for the record's representation and prints its expected "
     "length as asterisks: read it with the message before it (LTP=15: not an error of the file)"),
)


def _meaning(text):
    return next((m for pattern, m in MEANINGS if re.search(pattern, text)), "")


@dataclasses.dataclass(frozen=True)
class Checker:
    """A program found: its name, path, family (utility or prepro), kind (check, rewrite, derive), what it does,
    what it does to the numbers of a file (``numbers``: said before a run), and its version (read from its report
    after a run; None before)."""
    name: str
    path: str
    family: str
    kind: str
    does: str
    numbers: str = ""
    version: str = None


def _search_dirs(search):
    if search is not None:
        return [search] if isinstance(search, str) else list(search)
    dirs = [os.environ[v] for v in ("AZURE2_ENDF_UTILITIES", "AZURE2_PREPRO") if os.environ.get(v)]
    return dirs + os.environ.get("PATH", "").split(os.pathsep)


def checkers(search=None):
    """The checking codes found: ``search`` a directory, a program path or a list of them; else the directories of
    ``$AZURE2_ENDF_UTILITIES`` and ``$AZURE2_PREPRO``, then ``PATH``.  Names as :data:`PROGRAMS` (any case, with or
    without ``.exe``); the first found of each name wins.  No program is run here."""
    found = {}
    for place in _search_dirs(search):
        if os.path.isfile(place):
            candidates = [place]
        elif os.path.isdir(place):
            candidates = [os.path.join(place, f) for f in sorted(os.listdir(place))]
        else:
            continue
        for path in candidates:
            name = os.path.basename(path).lower()
            name = name[:-4] if name.endswith(".exe") else name
            if name in PROGRAMS and name not in found and os.access(path, os.X_OK) and os.path.isfile(path):
                family, kind, does, numbers = PROGRAMS[name]
                found[name] = Checker(name, os.path.abspath(path), family, kind, does, numbers)
    return [found[n] for n in PROGRAMS if n in found]


@dataclasses.dataclass
class CheckReport:
    """One program's run on a file: the program, its version (from the report), the exit status (None after a
    time-out), its findings, the report's text, the terminal text, the working directory (kept for a look)."""
    program: str
    path: str
    version: object
    exit_status: object
    findings: list
    report: str
    terminal: str
    workdir: str
    totals: object = None

    def __str__(self):
        head = (f"{self.program}" + (f" {self.version}" if self.version else "") +
                f": exit {self.exit_status}, {len(self.findings)} findings"
                + (f" (the report counts {self.totals[0]} errors, {self.totals[1]} warnings)" if self.totals else ""))
        lines = [head]
        for f in self.findings:
            lines.append(f"  {f.severity:7s} MAT={f.mat} MF={f.mf} MT={f.mt}: {f.message}"
                         + (f"\n          meaning: {f.meaning}" if f.meaning else ""))
        return "\n".join(lines)


@dataclasses.dataclass
class CheckResult:
    """``check()``'s answer: the module's own format check (``own``) and one :class:`CheckReport` per program."""
    path: str
    own: object
    reports: list

    @property
    def findings(self):
        return list(self.own.findings) + [f for r in self.reports for f in r.findings]

    def to_dict(self):
        return {"path": self.path, "own": self.own.to_dict(),
                "reports": [dict(dataclasses.asdict(r), findings=[dataclasses.asdict(f) for f in r.findings])
                            for r in self.reports]}

    def __str__(self):
        out = [str(self.own)]
        out += [str(r) for r in self.reports]
        if not self.reports:
            out.append("no checking code was run (none found: pyazr.endf.checkers() lists what is installed)")
        return "\n".join(out)


def _workdir(path, program):
    folder = tempfile.mkdtemp(prefix=f"endf_{program}_")
    name = "ENDFB.IN" if PROGRAMS[program][0] == "prepro" else os.path.basename(path)[-60:].replace(" ", "_")
    shutil.copyfile(path, os.path.join(folder, name))
    return folder, name


def _run(checker, folder, deck, timeout):
    """Run a program in ``folder`` with ``deck`` on standard input: (exit status or None, terminal text)."""
    try:
        p = subprocess.run([checker.path], input=deck, capture_output=True, text=True, cwd=folder, timeout=timeout,
                           errors="replace")
        return p.returncode, (p.stdout or "") + (p.stderr or "")
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode(errors="replace") if isinstance(e.stdout, bytes) else (e.stdout or "")
        return None, out
    except OSError as e:
        return -1, f"could not start {checker.path}: {e}"


_BLOCK = re.compile(r"^\s*(ERROR\(S\) FOUND|WARNING\(S\))\s+IN MAT=\s*(-?\d+),\s*MF=\s*(\d+),\s*MT=\s*(\d+)")
_MISSING = re.compile(r"^\s*ERROR\(S\) - MISSING SECTIONS IN MAT\s+(-?\d+)\s+MF\s+(\d+)")
_TOTALS = re.compile(r"Encountered\s+(\d+)\s+errors,\s+(\d+)\s+warnings")
_VERSION = re.compile(r"PROGRAM\s+(\w+)\s+VERSION\s+([\w.]+)")
_RECORD = re.compile(r"(RECORD|SEQUENCE) NUMBER\s+(\d+)")


def parse_report(program, text, header_text=()):
    """Findings of a report of CHECKR, FIZCON or PSYCHE (see the module's text), the version and the totals the
    report states (None when it states none)."""
    findings, current, continued = [], None, False
    version = _VERSION.search(text)
    totals = _TOTALS.search(text)
    echo = {t.strip() for t in header_text if t.strip()}
    if program in ("checkr", "fizcon"):
        for line in text.splitlines():
            m, miss = _BLOCK.match(line), _MISSING.match(line)
            if m:
                current = ("error" if m.group(1).startswith("ERROR") else "warning", int(m.group(2)),
                           int(m.group(3)), int(m.group(4)))
                continue
            if miss:
                current = ("error", int(miss.group(1)), int(miss.group(2)), 0)
                continue
            if "FILE IS UNSEQUENCED" in line:
                findings.append(_finding(program, ("note", 0, 0, 0), line))
                continue
            if _TOTALS.search(line) or line.strip().startswith("Done ") or not line.strip():
                if not line.strip():
                    continue
                current = None
                continue
            if current is not None and line.startswith("    "):
                previous = findings[-1] if findings else None
                # FIZCON wraps a message over lines (seen: "ENERGY RANGE FOR DISTRIBUTIONS IN LIST RECORDS" /
                # "INCONSISTENT WITH TAB1 RECORD SEQUENCE NUMBER 1", one error; "THE MAXIMUM INCIDENT ENERGY OF
                # 3.00000E+06 (EV)" / "SHOULD BE GREATER THAN OR EQUAL TO 2.00000E+07", one error): a line that
                # does not end in a number goes on
                if (program == "fizcon" and previous is not None and continued and
                        (previous.mat, previous.mf, previous.mt) == current[1:]):
                    findings[-1] = _finding(program, current, previous.message + " " + line)
                else:
                    findings.append(_finding(program, current, line))
                continued = not re.search(r"\d$", line.rstrip())
                continue
            continued = False
    elif program == "psyche":
        mf = mt = mat = 0
        for line in text.splitlines():
            s = line.strip()
            if re.match(r"^CHECK MATERIAL\s+(\d+)", s):
                mat = int(s.split()[2])
            elif re.match(r"^FILE\s+(\d+)$", s):
                mf, mt = int(s.split()[1]), 0
            elif re.match(r"^SECTION\s+(\d+)$", s):
                mt = int(s.split()[1])
            elif (mf and mt and s and not s.startswith("---") and s not in echo
                  and s != "CHECK ON ISOTOPE PROPERTIES" and not s.startswith("Done PSYCHE")):
                # PSYCHE prints tables and comments under a section: a comment that says something is out of range
                # or wrong is a warning, every other line a note
                warn = re.search(r"NOT IN RANGE|WARNING|ERROR|EXCEED|VIOLAT|INCONSISTENT|NEGATIVE", s)
                findings.append(_finding(program, ("warning" if warn else "note", mat, mf, mt), line))
    counted = sum(f.severity in ("error", "warning") for f in findings)  # CHECKR counts "SECTION CANNOT BE CHECKED"
    return findings, (version.group(2) if version else None), (
        (int(totals.group(1)), int(totals.group(2))) if totals else None), counted


def _finding(program, where, line):
    severity, mat, mf, mt = where
    text = " ".join(line.split())
    record = _RECORD.search(text)
    return Finding(severity, program.upper(), int(record.group(2)) if record and record.group(1) == "SEQUENCE"
                   else 0, mat, mf, mt, text, "", source=program, record=int(record.group(2)) if record else 0,
                   meaning=_meaning(text))


def _run_check(checker, path, timeout, header_text):
    folder, name = _workdir(path, checker.name)
    out = f"{checker.name}.out"
    deck = f"{name}\n{out}\nY\nDONE\n"
    status, terminal = _run(checker, folder, deck, timeout)
    report_path = os.path.join(folder, out)
    report = open(report_path, errors="replace").read() if os.path.isfile(report_path) else ""
    findings, version, totals, counted = [], None, None, 0
    if checker.name in ("checkr", "fizcon", "psyche") and report:
        findings, version, totals, counted = parse_report(checker.name, report, header_text)
    where = ("error", 0, 0, 0)
    if status is None:
        findings.append(_finding(checker.name, where, f"no answer within {timeout} s: the run was stopped"))
    elif not report:
        findings.append(_finding(checker.name, where, f"no report written (exit {status}); terminal: "
                                 + " ".join(terminal.split()[-30:])))
    elif checker.name in ("checkr", "fizcon", "psyche") and f"Done {checker.name.upper()}" not in report:
        findings.append(_finding(checker.name, where, f"the report ends before 'Done {checker.name.upper()}' "
                                 f"(exit {status}): the program stopped early"))
    elif totals is not None and totals[0] + totals[1] != counted:
        findings.append(_finding(checker.name, ("warning", 0, 0, 0), f"the report counts {sum(totals)} errors and "
                                 f"warnings, {counted} were read: see the report's text"))
    return CheckReport(checker.name, checker.path, version, status, findings, report, terminal, folder, totals)


def check(path, codes=None, search=None, timeout=300):
    """Check an ENDF-6 file: the module's own format check first (``own``, needs no program), then each checking
    code found (``codes``: names of :data:`PROGRAMS` of kind check; default CHECKR, FIZCON, PSYCHE), each in a fresh
    directory of its own with its deck, one process, ``timeout`` seconds.  A code asked for and not found is refused
    (:class:`pyazr.FormatRefused`)."""
    from .material import read
    own = check_format(path)
    found = {c.name: c for c in checkers(search)}
    wanted = list(codes) if codes else [n for n in CHECK_DEFAULT if n in found]
    for name in wanted:
        if name not in PROGRAMS or PROGRAMS[name][1] != "check":
            raise FormatRefused(what=f"{name!r} is not a checking code", allowed="one of " + ", ".join(
                n for n, v in PROGRAMS.items() if v[1] == "check"), fix="pyazr.endf.standardize() runs the codes "
                "that rewrite a file")
        if name not in found:
            raise FormatRefused(what=f"{name} is not found", allowed="programs found: " + (", ".join(found) or
                                "none"), fix="set AZURE2_ENDF_UTILITIES (or AZURE2_PREPRO) to the directory of the "
                                "programs, or pass search=")
    try:
        header_text = read(path).header.text
    except Exception:  # noqa: BLE001 - a file the reader refuses is still given to the codes
        header_text = ()
    return CheckResult(path, own, [_run_check(found[n], path, timeout, header_text) for n in wanted])


# ---------------------------------------------------------------------------------------------------- candidates
@dataclasses.dataclass(frozen=True)
class SectionChange:
    """What a rewrite changed in one section: ``kind`` unchanged, bookkeeping (MF=1/MT=451: text and directory),
    form (numbers written another way, every value equal), numbers (values changed), restructured (other lines or
    fields), added or removed; for numbers, how many
    fields changed, the largest relative change and the most significant digits a field lost."""
    mf: int
    mt: int
    kind: str
    numbers_changed: int = 0
    max_relative_change: float = 0.0
    digits_lost: int = 0
    lines_before: int = 0
    lines_after: int = 0
    example: str = ""


@dataclasses.dataclass
class Candidate:
    """A rewritten file beside the original (``path`` in a scratch directory): nothing is replaced until
    :meth:`accept`.  ``program``, ``report`` (its terminal text and listing)."""
    original: str
    path: str
    program: str
    report: str

    def diff(self):
        """Per section, what changed (:class:`SectionChange`), in file order."""
        return compare(self.original, self.path)

    def accept(self, to=None, overwrite=False):
        """Write the candidate to ``to`` (default: beside the original, ``<name>.<program>.endf``); an existing file
        is replaced only with ``overwrite=True`` (the original too).  Returns the path written."""
        if to is None:
            stem, _ = os.path.splitext(self.original)
            to = f"{stem}.{self.program}.endf"
        if os.path.exists(to) and not overwrite:
            raise FormatRefused(what=f"{to} exists", allowed="a file name that does not exist yet (an existing file is "
                                "replaced only when that is asked for)", fix="keep the candidate as a new file, or "
                                "ask for the existing file to be replaced (overwrite)")
        shutil.copyfile(self.path, to)
        return to

    def summary(self):
        lines = [f"{self.program}: candidate {self.path} for {self.original}"]
        for c in self.diff():
            if c.kind == "unchanged":
                continue
            extra = (f": {c.numbers_changed} numbers, largest relative change {c.max_relative_change:.3g}, up to "
                     f"{c.digits_lost} digits lost" if c.kind == "numbers" else
                     f" ({c.lines_before} -> {c.lines_after} lines)")
            lines.append(f"  MF={c.mf} MT={c.mt}: {c.kind}{extra}" + (f"; e.g. {c.example}" if c.example else ""))
        return "\n".join(lines)


def _sections_of(path):
    """{(MF, MT): [content of each line]} of the first material (MF=1 included), and the order."""
    out, order = {}, []
    with open(path, errors="replace") as f:
        for text in f.read().splitlines():
            t = text.ljust(80)
            try:
                mat, mf, mt = int(t[66:70] or 0), int(t[70:72] or 0), int(t[72:75] or 0)
            except ValueError:
                continue
            if mat <= 0 or mf == 0 or mt == 0:
                continue
            if (mf, mt) not in out:
                out[(mf, mt)] = []
                order.append((mf, mt))
            out[(mf, mt)].append(t[:66])
    return out, order


def _number(field):
    s = field.strip()
    if not s or not re.match(r"^[+-]?[\d.]", s):
        return None
    try:
        return read_float(field)
    except Exception:  # noqa: BLE001
        return None


def compare(original, rewritten):
    """The sections of two ENDF-6 files compared field by field (see :class:`SectionChange`)."""
    a, order_a = _sections_of(original)
    b, order_b = _sections_of(rewritten)
    changes = []
    for key in order_a + [k for k in order_b if k not in a]:
        mf, mt = key
        la, lb = a.get(key), b.get(key)
        if la is None or lb is None:
            changes.append(SectionChange(mf, mt, "added" if la is None else "removed", lines_before=len(la or []),
                                         lines_after=len(lb or [])))
            continue
        if la == lb:
            changes.append(SectionChange(mf, mt, "unchanged", lines_before=len(la), lines_after=len(lb)))
            continue
        if len(la) != len(lb) or key == (1, 451):
            changes.append(SectionChange(mf, mt, "bookkeeping" if key == (1, 451) else "restructured",
                                         lines_before=len(la), lines_after=len(lb)))
            continue
        n, worst, lost, example, other = 0, 0.0, 0, "", False
        for x, y in zip(la, lb):
            for k in range(6):
                fx, fy = x[11 * k:11 * k + 11], y[11 * k:11 * k + 11]
                if fx == fy:
                    continue
                vx, vy = _number(fx), _number(fy)
                if vx is None or vy is None:
                    other = True
                    continue
                if vx == vy:
                    continue
                n += 1
                rel = abs(vy - vx) / abs(vx) if vx else abs(vy)
                if rel > worst:
                    worst, example = rel, f"{fx.strip()} -> {fy.strip()}"
                lost = max(lost, significant_digits(fx) - significant_digits(fy))
        kind = "numbers" if n else ("restructured" if other else "form")
        changes.append(SectionChange(mf, mt, kind, n, worst, lost, len(la), len(lb), example))
    return changes


def _prepro_deck(program):
    if program == "dictin":
        return "ENDFB.IN\nENDFB.OUT\n"
    if program == "fixup":
        # fixup.f lines 593-615: column 1 (ZA and AWR, required) and column 5 = 2 (the built-in summation rules: the
        # built FIXUP 2023-1 refuses to run without option 5, "Cross Section Reconstruction is MANDATORY"); every
        # other test off; the standard file names
        return "10002000000000\nENDFB.IN\nENDFB.OUT\n\n"
    return ""


def standardize(path, program="stanef", search=None, timeout=300, mf2=False):
    """Run a code that rewrites a file (STANEF by default; FIXUP, DICTIN, ENDF2C) on a copy, in a directory of its
    own, and return the result as a :class:`Candidate` (``.diff()``, ``.summary()``, ``.accept(to=...)``).  The
    original is never touched.  ``mf2``: STANEF's option that also rewrites the numbers of MF=2/MT=151."""
    found = {c.name: c for c in checkers(search)}
    if program not in PROGRAMS or PROGRAMS[program][1] != "rewrite":
        raise FormatRefused(what=f"{program!r} does not rewrite a file", allowed="one of " + ", ".join(
            n for n, v in PROGRAMS.items() if v[1] == "rewrite"), fix="pyazr.endf.check() runs the checking codes")
    if program not in found:
        raise FormatRefused(what=f"{program} is not found", allowed="programs found: " + (", ".join(found) or "none"),
                            fix="set AZURE2_ENDF_UTILITIES or AZURE2_PREPRO to the directory of the programs")
    checker = found[program]
    folder, name = _workdir(path, program)
    if checker.family == "prepro":
        deck = _prepro_deck(program)
        if deck:
            with open(os.path.join(folder, f"{program.upper()}.INP"), "w") as f:
                f.write(deck)
        status, terminal = _run(checker, folder, "", timeout)
        out = os.path.join(folder, "ENDFB.OUT")
    else:
        out = os.path.join(folder, "candidate.endf")
        status, terminal = _run(checker, folder, f"{name}\ncandidate.endf\n\n{'0,Y,Y,Y' if mf2 else ''}\nDONE\n",
                                timeout)
    listing = "".join(open(os.path.join(folder, f), errors="replace").read() for f in sorted(os.listdir(folder))
                      if f.upper().endswith(".LST"))
    if status is None or not os.path.isfile(out) or os.path.getsize(out) == 0:
        raise FormatRefused(what=f"{program} wrote no file (exit {status})", allowed="a run that writes its output",
                            fix=f"look at {folder}: " + " ".join((terminal + listing).split()[-40:]))
    return Candidate(path, out, program, terminal + listing)
