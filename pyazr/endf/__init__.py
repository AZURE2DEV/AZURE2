"""ENDF-6 for evaluators: read an evaluation, see what it holds, draw it over a data set, check a file.

Every format statement comes from the ENDF-6 formats manual ENDF-102 (2023 edition, BNL), cited by its printed page
in the module that codes it.  The format work (records, number forms, interpolation laws, the sums of MF=6 LAW=2 and
LAW=5) is this package's; every conversion of physics (lab and c.m. energies, eta, the Rutherford and Mott cross
sections, the S-factor conversion) is a row of :mod:`pyazr.physics`, the engine's.

    m = pyazr.endf.read("h-002_He_004.endf")     # a Material: its sections by (MF, MT), its header
    print(m.summary())                            # MAT, target, projectile, library, sections in words, the text
    xs = m.cross_section(600)                     # MF=3/MT=600; xs(energies_eV) -> b
    c = m.curves(s.evaluate(x).data[i], session=s)  # the evaluation at the data set's points, its quantity and unit
    r = m.resonances()                            # MF=2/MT=151 LRF=7 as a record (the file's units and frame)
    project, report = pyazr.endf.to_project(r, "new.azr", material=m)   # a new project, Brune basis, all fixed
    record, losses = pyazr.endf.from_session(s, x)  # a session's particle channels as LRF=7, formal basis
    written = pyazr.endf.write(s, "new.endf", x)  # MF=1, MF=2 (LRF=7), MF=3 from the engine; LRP=2
    findings = pyazr.endf.check_format(path)      # the module's own check of the format, each finding with its page
    result = pyazr.endf.check(path)               # the own check, then CHECKR, FIZCON, PSYCHE when they are found
    candidate = pyazr.endf.standardize(path)      # STANEF on a copy: candidate.diff(), candidate.accept(to=...)

Errors: a malformed file is a :class:`pyazr.FormatError` (file:line:columns, what was found, what is allowed, the
manual's page); a valid file, or a part of it, that is not read or cannot give what was asked is a
:class:`pyazr.FormatRefused`; both carry ``what``, ``allowed``, ``fix``.
"""

from ..errors import FormatError, FormatRefused
from .records import NUMBER_FORMS, form_classic, form_e, form_widest, read_float, read_int
from .interpolation import interpolate
from .material import Material, Header, Section, read
from .distributions import CrossSection, Distribution, Product
from .summary import Summary
from .curves import CurveRequest, Curve, PairSpec, CoulombTerms, request_of, curves
from .format_check import Finding, FormatReport, check_format
from .resonances import LRF7, ParticlePair, Channel, Level, SpinGroup, read_resonances
from .model import Conversion, ConvertedLevel, Losses, to_project, from_session
from .writer import Written, write
from .checkers import (PROGRAMS, Checker, CheckReport, CheckResult, Candidate, SectionChange, checkers, check,
                       standardize, compare)

__all__ = ["read", "Material", "Header", "Section", "Summary", "CrossSection", "Distribution", "Product",
           "CurveRequest", "Curve", "PairSpec", "CoulombTerms", "request_of", "curves", "interpolate",
           "NUMBER_FORMS", "form_classic", "form_e", "form_widest", "read_float", "read_int", "FormatError",
           "FormatRefused", "Finding", "FormatReport", "check_format", "PROGRAMS", "Checker", "CheckReport",
           "CheckResult", "Candidate", "SectionChange", "checkers", "check", "standardize", "compare", "LRF7",
           "ParticlePair", "Channel", "Level", "SpinGroup", "read_resonances", "Conversion", "ConvertedLevel", "Losses",
           "to_project", "from_session", "Written", "write"]
