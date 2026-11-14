"""ENDF-6 for evaluators: read an evaluation, see what it holds, draw it over a data set, check a file.

Every format statement comes from the ENDF-6 formats manual ENDF-102 (2023 edition, BNL), cited by its printed page
in the module that codes it.  The format work (records, number forms, interpolation laws, the sums of MF=6 LAW=2 and
LAW=5) is this package's; every conversion of physics (lab and c.m. energies, eta, the Rutherford and Mott cross
sections, the S-factor conversion) is a row of :mod:`pyazr.physics`, the engine's.

    m = pyazr.endf.read("h-002_He_004.endf")     # a Material: its sections by (MF, MT), its header
    print(m.summary())                            # MAT, target, projectile, library, sections in words, the text
    xs = m.cross_section(600)                     # MF=3/MT=600; xs(energies_eV) -> b
    c = m.curves(s.evaluate(x).data[i], session=s)  # the evaluation at the data set's points, its quantity and unit
    findings = pyazr.endf.check_format(path)      # the module's own check of the format, each finding with its page

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

__all__ = ["read", "Material", "Header", "Section", "Summary", "CrossSection", "Distribution", "Product",
           "CurveRequest", "Curve", "PairSpec", "CoulombTerms", "request_of", "curves", "interpolate",
           "NUMBER_FORMS", "form_classic", "form_e", "form_widest", "read_float", "read_int", "FormatError",
           "FormatRefused", "Finding", "FormatReport", "check_format"]
