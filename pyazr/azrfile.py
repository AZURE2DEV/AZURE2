"""Editable, file-backed model of an AZURE2 ``.azr`` level scheme.

``AZURE2`` builds its R-matrix model from the ``<levels>`` block of a ``.azr``
file: one whitespace-separated line per *channel*, grouped into levels.  This
module parses that block into structured, mutable objects so a caller can add or
remove levels and channels from Python and write the result to a new file --
**without touching the original**.  Everything outside ``<levels>`` (config,
potential, data segments, target integration, ...) is preserved verbatim.

Because AZURE2 reads its model from a file, an edited scheme is *applied* by
writing a file and launching a fresh instance from it::

    from pyazr import AzrModel, azure2

    model = AzrModel.from_file("13N.azr")
    print(model)                                  # inspect the scheme
    model.remove_level(jpi="1/2+", energy=20)     # drop a background pole
    model.add_level(J=1.5, parity=+1, energy=4.1,
                    channels=[dict(pair=1, L=2, S=0.5, gamma=1000.0),
                              dict(pair=2, L=1, S=0.5, gamma=0.1)])
    path = model.write("13N_edited.azr")          # original file untouched
    azr = azure2(path)                            # run with the edited scheme

The 31 fields of a channel line match ``NucLine`` in the AZURE2 source
(``include/NucLine.h``); note the file stores ``2*S`` and ``2*L`` as integers,
which the ``S`` / ``L`` properties convert.
"""

import math
import os
import re
import tempfile
from typing import List, Optional

# Field order of a <levels> channel line, matching NucLine's read order.
_FIELDS = [
    "levelJ", "levelPi", "levelE", "levelFix", "aa", "ir", "s2", "l2",
    "levelID", "isActive", "channelFix", "gamma", "j1", "pi1", "j2", "pi2",
    "e2", "m1", "m2", "z1", "z2", "entranceSepE", "sepE", "j3", "pi3", "e3",
    "pType", "chRad", "g1", "g2", "ecMultMask",
]
# Two fields were appended later (binding energy, and a flag saying the gamma
# column holds a reduced width amplitude rather than a partial width or ANC).
# Files written before they existed have 31 fields and must keep loading, so
# they are optional and a short line is padded on read rather than rejected.
_OPTIONAL_FIELDS = ["bindingEnergy", "gammaIsRWA"]
_FIELDS_ALL = _FIELDS + _OPTIONAL_FIELDS
_IDX = {name: i for i, name in enumerate(_FIELDS_ALL)}
_NFIELDS = len(_FIELDS)
_NFIELDS_MAX = len(_FIELDS_ALL)
_OPTIONAL_DEFAULTS = ["0", "0"]


def _token_spans(line):
    """(start, end) of every whitespace-separated token in ``line``."""
    return [m.span() for m in re.finditer(r"\S+", line)]


def _isnum(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def _fmt(x):
    """Round-trippable, compact token for a number."""
    if isinstance(x, int):
        return str(x)
    if float(x).is_integer() and abs(x) < 1e15:
        return str(int(x))
    return repr(float(x))


# -- the <thm> block: the engine's (Config::ReadThmBlock) and the GUI's
#    (ThmSettings in gui/src/ThmOptionsDialog.cpp) rules, ported -------------

_THM_GLOBAL_DEFAULTS = {"entranceL": "incoherent", "vertex": "constant",
                        "kinematics": "lacognata", "coulombIntegral": False,
                        "spectatorEnergy": 0.0}
_THM_FLOAT_PREFIX = re.compile(r"[ \t\n\r\f\v]*[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")
_THM_INT_PREFIX = re.compile(r"[ \t\n\r\f\v]*[+-]?\d+")


# An explicit prior centre in <parameterSettings> (EData::ReadPriorCentres).
_PRIOR_ROW = re.compile(r"^\s*segment_(\d+)_(norm|energy_shift)\s+prior_centre\s+(\S+)\s*$")
_PRIOR_COMMENT = "# Prior centres: segment_N_norm|segment_N_energy_shift prior_centre value"


class _Classic:
    """Marker for :meth:`AzrModel.set_prior_centres`: remove the explicit
    centre, so that the segment field is the centre again."""

    def __repr__(self):
        return "CLASSIC"


CLASSIC = _Classic()


def _thm_default_settings():
    s = dict(_THM_GLOBAL_DEFAULTS)
    s.update(spectatorByPair={}, weight={}, weightTest={}, experiments={})
    return s


# -- THM experiments: experiment[<name>] key=value ... (src/ThmExperiment.cpp) --

# AME2020 nuclear masses (u): atomic mass - Z m_e + electron binding (Lunney
# et al. 2003), the engine's table (ThmExperiment.cpp, kNuclides).
_THM_NUCLIDES = {
    "n": (0, 1, 1.0086649159), "p": (1, 1, 1.0072764675), "d": (1, 2, 2.0135532134),
    "t": (1, 3, 3.0155007169), "3He": (2, 3, 3.0149322434), "4He": (2, 4, 4.0015061756),
    "6Li": (3, 6, 6.0134773618), "7Li": (3, 7, 7.0143579087), "9Be": (4, 9, 9.0099891714),
    "10B": (5, 10, 10.0101946866), "11B": (5, 11, 11.0065629921),
    "12C": (6, 12, 11.9967096429), "13C": (6, 13, 13.0000644777),
    "14N": (7, 14, 13.9992355671), "15N": (7, 15, 14.9962704619),
    "16O": (8, 16, 15.9905282131), "17O": (8, 17, 16.9947453497),
    "18O": (8, 18, 17.9947732059), "19F": (9, 19, 18.9934689019),
    "20Ne": (10, 20, 19.9869581831), "23Na": (11, 23, 22.9837396828),
    "24Mg": (12, 24, 23.9784646239),
}
_THM_BACKGROUNDS = ("none", "const", "linear", "quadratic")
_THM_EXPERIMENT_KEYS = ("segments", "background", "beam", "target", "spectator", "Ebeam",
                        "lineshape", "ps", "psNodes", "distortion", "opticalAA", "opticalSF",
                        "spectatorAngle", "distortionRef", "distortionRatio", "boundState",
                        "theta", "vertexModel", "spectatorAngles", "spectatorAngleNodes",
                        "cbackground")
_THM_DISTORTION_DETAIL = ("spectatorAngle", "distortionRef", "distortionRatio", "boundState")
_THM_NAME = re.compile(r"[A-Za-z0-9_.+-]+")
_THM_WHOLE_INT = re.compile(r"[+-]?\d+")


def _thm_whole_double(text):
    """A whole token as a finite number (operator>> consuming everything)."""
    m = _THM_FLOAT_PREFIX.match(text)
    if not m or m.start(0) != 0 or text[:1].isspace() or m.end(0) != len(text):
        return None
    x = float(m.group(0))
    return x if math.isfinite(x) else None


def _thm_whole_int(text):
    return int(text) if _THM_WHOLE_INT.fullmatch(text) else None


def _thm_nuclide(text):
    """ThmNuclide::Parse: (name, Z, A, mass) or ValueError."""
    if text in _THM_NUCLIDES:
        return (text,) + _THM_NUCLIDES[text]
    f = text.split(",")
    if len(f) == 3:
        Z, A, mass = _thm_whole_int(f[0]), _thm_whole_int(f[1]), _thm_whole_double(f[2])
        if (Z is not None and A is not None and mass is not None and Z >= 0
                and A >= 1 and Z <= A and mass > 0.0 and abs(mass - A) < 1.0):
            return (text, Z, A, mass)
        raise ValueError(f"'{text}': expected Z,A,mass with 0 <= Z <= A and the "
                         "nuclear mass in u (within 1 u of A)")
    raise ValueError(f"unknown nuclide '{text}' (known: {' '.join(_THM_NUCLIDES)}; "
                     "or Z,A,mass in u)")


def _thm_segment_list(text):
    out = []
    for item in text.split(","):
        dash = item.find("-", 1)
        if dash < 0:
            lo = hi = _thm_whole_int(item)
            ok = lo is not None
        else:
            lo, hi = _thm_whole_int(item[:dash]), _thm_whole_int(item[dash + 1:])
            ok = lo is not None and hi is not None and hi >= lo
        if not ok or lo < 1:
            raise ValueError(f"segments='{text}': expected segment numbers >= 1 "
                             "like 1,2,4-6")
        if hi > 100000:                       # ThmExperiment.cpp kMaxThmSegment
            raise ValueError(f"segments='{text}': segment numbers go up to 100000")
        seen = set(out)
        for k in range(lo, hi + 1):
            if k in seen:
                raise ValueError(f"segments='{text}': segment {k} is listed twice")
            seen.add(k)
            out.append(k)
    return sorted(out)


def _thm_ps_window(text):
    """"pmin-pmax" (MeV/c, 0 <= pmin <= pmax): the '-' with a number on both
    sides; (lo, hi) or None."""
    for k in range(1, len(text) - 1):
        if text[k] != "-":
            continue
        lo, hi = _thm_whole_double(text[:k]), _thm_whole_double(text[k + 1:])
        if lo is not None and hi is not None:
            return (lo, hi) if 0.0 <= lo <= hi else None
    return None


def _thm_parse_ps(value):
    """ParsePs (src/ThmExperiment.cpp): the value of ps= as a normalized
    string, or ValueError (message without the experiment prefix)."""
    usage = (f"ps='{value}': expected delta, hulthen:pmin-pmax, hulthen:a,b:pmin-pmax "
             "(a, b in fm^-1), gauss:FWHM:pmin-pmax or table:<file> (momenta in MeV/c, "
             "0 <= pmin <= pmax)")
    if value == "delta":
        return value
    f = value.split(":")
    if len(f) >= 2 and f[0] == "table":
        if not value[6:]:
            raise ValueError(usage)
        return value
    if f[0] == "hulthen" and len(f) in (2, 3):
        if len(f) == 3:
            ab = f[1].split(",")
            a = _thm_whole_double(ab[0]) if len(ab) == 2 else None
            b = _thm_whole_double(ab[1]) if len(ab) == 2 else None
            if a is None or b is None or not (a > 0.0 and b > a):
                raise ValueError(f"ps='{value}': Hulthen a,b in fm^-1 with 0 < a < b "
                                 "(deuteron: 0.2317,1.202)")
        if _thm_ps_window(f[-1]) is None:
            raise ValueError(usage)
        return value
    if f[0] == "gauss" and len(f) == 3:
        w = _thm_whole_double(f[1])
        if w is None or not w > 0.0:
            raise ValueError(f"ps='{value}': the FWHM of |phi(p_s)|^2 must be a number "
                             "> 0 (MeV/c)")
        if _thm_ps_window(f[2]) is None:
            raise ValueError(usage)
        return value
    raise ValueError(usage)


# Global optical potentials (src/ThmOptical.cpp): name -> projectiles (Z, A),
# reference, target mass range, projectile lab-energy range (MeV).
_THM_GLOBAL_OPTICALS = {
    "ancai06": (((1, 2),), "An & Cai, PRC 73 (2006) 054605", 12, 238, 0.0, 183.0),
    "daehnick80": (((1, 2),), "Daehnick, Childs & Vrcelj, PRC 21 (1980) 2253 "
                   "(global set)", 27, 238, 11.8, 90.0),
    "kd03": (((0, 1), (1, 1)), "Koning & Delaroche, NPA 713 (2003) 231 (global)", 24, 209,
             0.001, 200.0),
    "bg71": (((1, 3), (2, 3)), "Becchetti & Greenlees (1971)", 40, 208, 1.0, 40.0),
    "liang09": (((2, 3),), "Liang, Li & Cai, J. Phys. G 36 (2009) 085104", 9, 208, 0.0, 270.0),
    "mcfadden66": (((2, 4),), "McFadden & Satchler, NPA 84 (1966) 177", 16, 208, 1.0, 25.0),
    "avrigeanu94": (((2, 4),), "Avrigeanu, Hodgson & Avrigeanu, PRC 49 (1994) 2136", 16, 250,
                    1.0, 73.0),
}
_THM_SPECIES = {(0, 1): "n", (1, 1): "p", (1, 2): "d", (1, 3): "t", (2, 3): "3He", (2, 4): "4He"}


def _thm_global_optical(value):
    """(name, extrapolate) of a global optical potential value, or None."""
    name, colon, option = value.partition(":")
    if name in _THM_GLOBAL_OPTICALS and (not colon or option == "extrapolate"):
        return name, bool(colon)
    return None


def _thm_parse_optical(key, value):
    """ParseOptical (src/ThmExperiment.cpp): plane | coulomb | a global
    optical potential <name>[:extrapolate] | ten numbers V,R,a,W,RW,aW,WD,RD,aD,RC;
    the value as given, or ValueError."""
    if value in ("plane", "coulomb"):
        return value
    if value[:1].isalpha():
        if _thm_global_optical(value) is None:
            raise ValueError(f"{key}='{value}': expected plane, coulomb, a global optical "
                             f"potential ({', '.join(_THM_GLOBAL_OPTICALS)}, optionally "
                             ":extrapolate) or ten numbers V,R,a,W,RW,aW,WD,RD,aD,RC")
        return value
    f = value.split(",")
    p = [_thm_whole_double(t) for t in f] if len(f) == 10 else None
    ok = p is not None and all(v is not None for v in p)
    if ok:
        for t in range(3):
            if p[3 * t] != 0.0 and not (p[3 * t + 1] > 0.0 and p[3 * t + 2] > 0.0):
                ok = False
        ok = ok and p[9] >= 0.0
    if not ok:
        raise ValueError(f"{key}='{value}': expected plane, coulomb or ten numbers "
                         "V,R,a,W,RW,aW,WD,RD,aD,RC (MeV and fm: real volume, imaginary "
                         "volume and imaginary surface Woods-Saxon, depths > 0 "
                         "attractive/absorptive, radii and diffusenesses > 0 where the "
                         "depth is not 0, the Coulomb radius RC >= 0, 0 = point charge)")
    return value


def _thm_parse_distortion_key(key, value):
    """The distortion keys of an experiment line: the value to keep, or
    ValueError (message without the experiment prefix)."""
    if key == "distortion":
        if value in ("none", "coulomb", "optical") or (value.startswith("table:")
                                                        and len(value) > 6):
            return value
        raise ValueError(f"distortion='{value}': expected none, coulomb, optical or "
                         "table:<file>")
    if key in ("opticalAA", "opticalSF"):
        return _thm_parse_optical(key, value)
    if key == "spectatorAngle":
        if value == "qf":
            return value
        a = _thm_whole_double(value[3:]) if value.startswith("cm:") else _thm_whole_double(value)
        if a is None or not 0.0 <= a <= 180.0:
            raise ValueError(f"spectatorAngle='{value}': expected qf, a lab angle in "
                             "degrees (0-180) or cm:<degrees> (0-180)")
        return value
    if key == "distortionRef":
        if _thm_whole_double(value) is None:
            raise ValueError(f"distortionRef='{value}': expected the reference energy "
                             "E_ref in MeV (c.m. of x + A)")
        return value
    if key == "distortionRatio":
        if value not in ("dwpw", "dw"):
            raise ValueError(f"distortionRatio='{value}': expected dwpw or dw")
        return value
    if key == "boundState":
        f = value.split(":")
        ok = len(f) in (1, 2) and f[0] in ("whittaker", "yukawa")
        if ok and len(f) == 2:
            r = _thm_whole_double(f[1])
            ok = r is not None and 0.0 <= r <= 50.0
        if not ok:
            raise ValueError(f"boundState='{value}': expected whittaker or yukawa, "
                             "optionally :rmin in fm (0-50)")
        return value
    raise KeyError(key)


def _thm_spin(text):
    """ReadSpin (ThmExperiment.cpp): "1/2", "0.5", "2" -> a multiple of 1/2, or None."""
    if not text or text[0] in "+-":
        return None
    if "/" in text:
        num, _, den = text.partition("/")
        n = _thm_whole_int(num)
        return n / 2.0 if den == "2" and n is not None and n >= 0 else None
    x = _thm_whole_double(text)
    if x is None or x < 0.0 or abs(2.0 * x - round(2.0 * x)) > 1e-9:
        return None
    return x


def _thm_spin_text(j):
    twice = int(round(2.0 * j))
    return str(twice // 2) if twice % 2 == 0 else f"{twice}/2"


def _thm_parse_cbackground(value):
    """ParseThmCoherentBackground: a list of term dicts, ValueError if AZURE2
    would refuse the value."""
    form = ("expected <J><+|->:<exit pair key>[:<s>,<l>,<s'>,<l'>][:const|:linear]"
            "[=<Re c0>,<Im c0>[,<Re c1>,<Im c1>]] terms separated by ';' (a value "
            "followed by f is fixed)")
    items = value.split(";")
    terms = []
    for item in items:
        bad = f"cbackground: '{item}': "
        if not item:
            raise ValueError(f"cbackground='{value}': an empty term; " + form)
        spec, eq, values = item.partition("=")
        f = spec.split(":")
        if not 2 <= len(f) <= 4:
            raise ValueError(bad + form)
        jpi = f[0]
        J = _thm_spin(jpi[:-1]) if len(jpi) >= 2 and jpi[-1] in "+-" else None
        if J is None:
            raise ValueError(bad + f"J^pi '{jpi}': expected e.g. 1/2+ or 2-")
        key = _thm_whole_int(f[1])
        if key is None or key < 1 or f[1][0] in "+-":
            raise ValueError(bad + f"exit pair key '{f[1]}': expected a positive whole number")
        t = {"J": J, "parity": 1 if jpi[-1] == "+" else -1, "exit": key,
             "channels": None, "form": "const", "values": None}
        have_form = False
        for x in f[2:]:
            if x in ("const", "linear"):
                if have_form:
                    raise ValueError(bad + form)
                have_form = True
                t["form"] = x
                continue
            if t["channels"] is not None or have_form:
                raise ValueError(bad + form)
            c = x.split(",")
            ok = len(c) == 4
            if ok:
                s_in, l_in = _thm_spin(c[0]), _thm_whole_int(c[1])
                s_out, l_out = _thm_spin(c[2]), _thm_whole_int(c[3])
                ok = (None not in (s_in, l_in, s_out, l_out) and l_in >= 0 and l_out >= 0
                      and c[1][0] != "+" and c[3][0] != "+")
            if not ok:
                raise ValueError(bad + f"channels '{x}': expected <s>,<l>,<s'>,<l'> "
                                 "(entrance and exit channel spin and l)")
            t["channels"] = (s_in, l_in, s_out, l_out)
        if eq:
            v = values.split(",")
            n = 2 if t["form"] == "const" else 4
            if len(v) != n:
                raise ValueError(bad + f"expected {n} values after '=' (" +
                                 ("Re c0, Im c0" if n == 2 else "Re c0, Im c0, Re c1, Im c1")
                                 + ")")
            vals = []
            for w in v:
                fixed = w.endswith("f")
                x = _thm_whole_double(w[:-1] if fixed else w)
                if x is None:
                    raise ValueError(bad + f"value '{w}': expected a number, optionally "
                                     "followed by f (fixed)")
                vals.append((x, fixed))
            t["values"] = vals
        terms.append(t)
    return terms


def _thm_format_cbackground(terms):
    """FormatThmCoherentBackground: the canonical value of a term list."""
    out = []
    for t in terms:
        s = f"{_thm_spin_text(t['J'])}{'+' if t['parity'] > 0 else '-'}:{t['exit']}"
        if t["channels"] is not None:
            a, b, c, d = t["channels"]
            s += f":{_thm_spin_text(a)},{b},{_thm_spin_text(c)},{d}"
        if t["form"] == "linear":
            s += ":linear"
        if t["values"] is not None:
            s += "=" + ",".join(_thm_plain_number(x) + ("f" if fx else "")
                                for x, fx in t["values"])
        out.append(s)
    return ";".join(out)


def _thm_parse_experiment(line, experiments):
    """ParseThmExperimentLine: merge one experiment line (comment stripped,
    trimmed) into ``experiments`` (name -> record); ValueError if AZURE2
    would refuse it."""
    close = line.find("]")
    if not line.startswith("experiment[") or close < 0:
        raise ValueError("<thm> expected experiment[<name>] key=value ...")
    name = line[len("experiment["):close]
    if not _THM_NAME.fullmatch(name):
        raise ValueError(f"<thm> experiment[{name}]: a name is letters, digits "
                         "and _ - . + only")
    where = f"<thm> experiment[{name}]: "
    rest = line[close + 1:]
    if rest and rest[0] not in " \t":
        raise ValueError(where + "expected a space after ']'")
    work = dict(experiments.get(name, {"keys": []}))
    work["keys"] = list(work["keys"])
    tokens = rest.split()
    if not tokens:
        raise ValueError(where + "no key=value given")
    for token in tokens:
        eq = token.find("=")
        if eq <= 0 or eq + 1 == len(token):
            raise ValueError(where + f"'{token}' is not key=value")
        key, value = token[:eq], token[eq + 1:]
        if key in work["keys"]:
            raise ValueError(where + f"key '{key}' is given twice")
        if key == "segments":
            try:
                work["segments"] = _thm_segment_list(value)
            except ValueError as err:
                raise ValueError(where + str(err)) from None
        elif key == "background":
            if value not in _THM_BACKGROUNDS:
                raise ValueError(where + f"background='{value}': expected none, "
                                 "const, linear or quadratic")
            work["background"] = value
        elif key in ("beam", "target", "spectator"):
            try:
                work[key] = _thm_nuclide(value)
            except ValueError as err:
                raise ValueError(where + f"{key}: {err}") from None
        elif key == "Ebeam":
            x = _thm_whole_double(value)
            if x is None or not x > 0.0:
                raise ValueError(where + f"Ebeam='{value}': expected the lab beam "
                                 "energy in MeV, > 0")
            work["Ebeam"] = x
        elif key == "lineshape":
            if value not in ("on", "off"):
                raise ValueError(where + f"lineshape='{value}': expected on or off")
            work["lineshape"] = value == "on"
        elif key == "ps":
            try:
                work["ps"] = _thm_parse_ps(value)
            except ValueError as err:
                raise ValueError(where + str(err)) from None
        elif key == "psNodes":
            n = _thm_whole_int(value)
            if n is None or not 1 <= n <= 64:
                raise ValueError(where + f"psNodes='{value}': expected a whole number of "
                                 "Gauss-Legendre nodes, 1 to 64")
            work["psNodes"] = n
        elif key in ("distortion", "opticalAA", "opticalSF") + _THM_DISTORTION_DETAIL:
            try:
                work[key] = _thm_parse_distortion_key(key, value)
            except ValueError as err:
                raise ValueError(where + str(err)) from None
        elif key == "theta":
            window = None if value == "all" else _thm_ps_window(value)
            if value != "all" and (window is None or window[1] > 180.0):
                raise ValueError(where + f"theta='{value}': expected all or thmin-thmax, "
                                 "the c.m. angles of the exit pair relative to p_xA in "
                                 "degrees, 0 <= thmin <= thmax <= 180")
            work["theta"] = value
        elif key == "vertexModel":
            if value not in ("pw", "dw"):
                raise ValueError(where + f"vertexModel='{value}': expected pw or dw")
            work["vertexModel"] = value
        elif key == "spectatorAngles":
            v = value[3:] if value.startswith("cm:") else value
            window = None if v.startswith("table:") else _thm_ps_window(v)
            if not ((v.startswith("table:") and len(v) > 6)
                    or (window is not None and window[1] <= 180.0)):
                raise ValueError(where + f"spectatorAngles='{value}': expected thmin-thmax "
                                 "(lab polar angles of the spectator to the beam, degrees, "
                                 "0 <= thmin <= thmax <= 180), cm:thmin-thmax (c.m.), "
                                 "table:<file> or cm:table:<file> (angle, acceptance)")
            work["spectatorAngles"] = value
        elif key == "spectatorAngleNodes":
            n = _thm_whole_int(value)
            if n is None or not 1 <= n <= 64:
                raise ValueError(where + f"spectatorAngleNodes='{value}': expected a whole "
                                 "number of Gauss-Legendre nodes, 1 to 64")
            work["spectatorAngleNodes"] = n
        elif key == "cbackground":
            try:
                work["cbackground"] = _thm_format_cbackground(_thm_parse_cbackground(value))
            except ValueError as err:
                raise ValueError(where + str(err)) from None
        else:
            raise ValueError(where + f"unknown key '{key}' (keys: segments, "
                             "background, beam, target, spectator, Ebeam, lineshape, ps, "
                             "psNodes, distortion, opticalAA, opticalSF, spectatorAngle, "
                             "distortionRef, distortionRatio, boundState, theta, "
                             "vertexModel, spectatorAngles, spectatorAngleNodes, "
                             "cbackground)")
        work["keys"].append(key)
    experiments[name] = work


def _thm_check_experiments(experiments):
    """CheckThmExperiments: ValueError if the set is refused."""
    owner = {}
    for name, x in experiments.items():
        where = f"<thm> experiment[{name}]: "
        if "segments" not in x["keys"]:
            raise ValueError(where + "segments= is required")
        kin = sum(k in x["keys"] for k in ("beam", "target", "spectator", "Ebeam"))
        if kin not in (0, 4):
            raise ValueError(where + "beam, target, spectator and Ebeam go together "
                             "(all four or none)")
        if x.get("lineshape") and kin != 4:
            raise ValueError(where + "lineshape=on needs the kinematics of the reaction: "
                             "beam, target, spectator and Ebeam")
        if x.get("lineshape") and x.get("cbackground"):
            raise ValueError(where + "lineshape=on and cbackground= cannot be combined: "
                             "the Coulomb line shape N_C multiplies the level amplitudes "
                             "only, and its phase (E_lambda - E - i Gamma/2)^(-i zeta) "
                             "depends on the energy unit unless every term of the "
                             "amplitude carries it")
        if x.get("ps", "delta") != "delta" and kin != 4:
            raise ValueError(where + "a ps window (ps=hulthen|gauss|table) needs the "
                             "kinematics of the reaction: beam, target, spectator and "
                             "Ebeam (mu_sx and the spectator momenta reached at each "
                             "energy)")
        if "psNodes" in x["keys"] and x.get("ps", "delta") == "delta":
            raise ValueError(where + "psNodes= needs a ps window (ps=hulthen|gauss|table)")
        dist = x.get("distortion", "none")
        computed = dist in ("coulomb", "optical")
        if computed and kin != 4:
            raise ValueError(where + f"distortion={dist} needs the kinematics of the "
                             "reaction: beam, target, spectator and Ebeam")
        if ("opticalAA" in x["keys"] or "opticalSF" in x["keys"]) and dist != "optical":
            raise ValueError(where + "opticalAA= and opticalSF= need distortion=optical")
        for key in _THM_DISTORTION_DETAIL:
            if key in x["keys"] and not computed:
                raise ValueError(where + f"{key}= needs distortion=coulomb or "
                                 "distortion=optical")
        if x.get("vertexModel", "pw") == "dw":
            if not computed:
                raise ValueError(where + "vertexModel=dw builds the vertex from the "
                                 "distorted waves of a + A and s + F: it needs "
                                 "distortion=coulomb or distortion=optical (and the "
                                 "kinematics)")
            for key in ("distortionRef", "distortionRatio"):
                if key in x["keys"]:
                    raise ValueError(where + f"{key}= belongs to the distortion factor "
                                     "R(E), which vertexModel=dw replaces (R is not "
                                     "applied; the DW vertex carries the energy "
                                     "dependence)")
            if x.get("theta", "all") != "all":
                raise ValueError(where + "theta= (fixed-angle observable) is not "
                                 "available with vertexModel=dw: the distorted source has "
                                 "every m_l about p_xA, which the fixed-angle sum does "
                                 "not carry")
        # The acceptance (ThmExperiment.cpp CheckThmExperiments): ps= (a |p_s|
        # cut with the momentum distribution) and spectatorAngles= describe one
        # set of directions at fixed E.
        has_ps = x.get("ps", "delta") != "delta"
        if ("spectatorAngles" in x or has_ps) and "spectatorAngle" in x["keys"]:
            if "spectatorAngles" in x:
                raise ValueError(where + "spectatorAngle= (one direction) and "
                                 "spectatorAngles= (a window) exclude each other")
            raise ValueError(where + "spectatorAngle= (one direction) and a ps window "
                             "exclude each other: the window accepts every direction "
                             "whose |p_s| lies in it (write spectatorAngles=cm:t-t for one "
                             "direction with the cut)")
        if "spectatorAngles" in x:
            if not computed and not has_ps:
                raise ValueError(where + "spectatorAngles= averages over the accepted "
                                 "spectator directions the plane-wave vertex (with the "
                                 "momentum distribution of ps=hulthen|gauss|table), or R(E) "
                                 "and the DW vertex (distortion=coulomb|optical); with "
                                 "neither it has nothing to average")
            if "psNodes" in x["keys"]:
                raise ValueError(where + "with spectatorAngles= the nodes are "
                                 "spectatorAngleNodes= (per c.m. interval); psNodes= is the "
                                 "node count of a ps window without spectatorAngles=")
        elif "spectatorAngleNodes" in x["keys"]:
            raise ValueError(where + "spectatorAngleNodes= needs a spectator-direction "
                             "window (spectatorAngles=)")
        if (computed and x.get("vertexModel", "pw") != "dw"
                and x.get("ps", "delta") != "delta"
                and x.get("distortionRatio", "dwpw") == "dw"):
            raise ValueError(where + "distortionRatio=dw multiplies the model by |M|^2, "
                             "which carries the momentum distribution |phi(q)|^2 at the "
                             "spectator direction; a ps window already weights the model "
                             "with |phi(p_s)|^2 of data divided by it, so it would count "
                             "twice. Use distortionRatio=dwpw (the default) with a window")
        if dist == "optical" and kin == 4:
            _thm_check_global_optical(where, x)
        for k in x["segments"]:
            if k in owner:
                raise ValueError(where + f"segment {k} is already in "
                                 f"experiment[{owner[k]}]")
            owner[k] = name


def _thm_nucleus_name(Z, A):
    for name, (z, a, _m) in _THM_NUCLIDES.items():
        if (z, a) == (Z, A):
            return name
    return f"(Z,A)=({Z},{A})"


def _thm_check_global_optical(where, x):
    """The checks of ThmDistortion::Setup that need no data: the global
    potential describes a partner of its channel, and the target mass (and the
    a + A lab energy) lies in its validity range unless :extrapolate.  The
    s + F lab energy follows E and is checked by the engine over the data."""
    beam, target, spec = x["beam"], x["target"], x["spectator"]
    eaa = x["Ebeam"] * target[3] / (beam[3] + target[3])
    sfF = (beam[1] + target[1] - spec[1], beam[2] + target[2] - spec[2],
           beam[3] + target[3] - spec[3])
    channels = (("opticalAA", "a + A", (beam[1], beam[2], beam[3]),
                 (target[1], target[2], target[3]), eaa),
                ("opticalSF", "s + F", (spec[1], spec[2], spec[3]), sfF, None))
    for key, label, one, two, ecm in channels:
        g = _thm_global_optical(x.get(key, "coulomb"))
        if g is None:
            continue
        name, extrapolate = g
        projectiles, ref, amin, amax, emin, emax = _THM_GLOBAL_OPTICALS[name]
        if (one[0], one[1]) in projectiles:
            proj, targ = one, two
        elif (two[0], two[1]) in projectiles:
            proj, targ = two, one
        else:
            species = " ".join(_THM_SPECIES[p] for p in projectiles)
            raise ValueError(where + f"{key}={name} ({ref}, for {species}): the {label} "
                             f"channel is {_thm_nucleus_name(one[0], one[1])} + "
                             f"{_thm_nucleus_name(two[0], two[1])}, which it does not "
                             "describe (heavy-ion and other channels take the ten numbers "
                             "V,R,a,W,RW,aW,WD,RD,aD,RC)")
        problems = []
        if not amin <= targ[1] <= amax:
            problems.append(f"target A = {targ[1]} (valid {amin}-{amax})")
        if ecm is not None:
            elab = ecm * (proj[2] + targ[2]) / targ[2]
            if not emin <= elab <= emax:
                problems.append(f"lab energy of {_thm_nucleus_name(proj[0], proj[1])} "
                                f"{elab:.4g} MeV (valid {emin:g}-{emax:g} MeV)")
        if problems and not extrapolate:
            raise ValueError(where + f"{key}={name} ({ref}) is outside its validity range "
                             f"for {_thm_nucleus_name(proj[0], proj[1])} + "
                             f"{_thm_nucleus_name(targ[0], targ[1])}: " + ", ".join(problems)
                             + f". Write {name}:extrapolate to use it there anyway (warned), "
                             "or give the ten numbers")


def _thm_experiment_record(x):
    """The public form of a parsed experiment."""
    out = {"segments": list(x["segments"]), "background": x.get("background", "none")}
    for key in ("beam", "target", "spectator"):
        if key in x:
            out[key] = x[key][0]
    if "Ebeam" in x:
        out["Ebeam"] = x["Ebeam"]
    if x.get("lineshape"):
        out["lineshape"] = True
    if "ps" in x:
        out["ps"] = x["ps"]
    if "psNodes" in x:
        out["psNodes"] = x["psNodes"]
    for key in ("distortion", "opticalAA", "opticalSF") + _THM_DISTORTION_DETAIL:
        if key in x:
            out[key] = x[key]
    if "theta" in x:
        out["theta"] = x["theta"]
    if "vertexModel" in x:
        out["vertexModel"] = x["vertexModel"]
    for key in ("spectatorAngles", "spectatorAngleNodes", "cbackground"):
        if key in x:
            out[key] = x[key]
    return out


def _thm_experiment_line(name, rec):
    """The canonical line for a record (segments as a plain comma list)."""
    parts = [f"experiment[{name}]", "segments=" + ",".join(str(k) for k in rec["segments"])]
    if rec.get("background", "none") != "none":
        parts.append(f"background={rec['background']}")
    for key in ("beam", "target", "spectator"):
        if key in rec:
            parts.append(f"{key}={rec[key]}")
    if "Ebeam" in rec:
        parts.append(f"Ebeam={_thm_number(float(rec['Ebeam']))}")
    if rec.get("lineshape"):
        parts.append("lineshape=on")
    if "ps" in rec:
        parts.append(f"ps={rec['ps']}")
    if "psNodes" in rec:
        parts.append(f"psNodes={int(rec['psNodes'])}")
    for key in ("distortion", "opticalAA", "opticalSF") + _THM_DISTORTION_DETAIL:
        if key in rec:
            parts.append(f"{key}={rec[key]}")
    if "theta" in rec:
        parts.append(f"theta={rec['theta']}")
    if "vertexModel" in rec:
        parts.append(f"vertexModel={rec['vertexModel']}")
    if "spectatorAngles" in rec:
        parts.append(f"spectatorAngles={rec['spectatorAngles']}")
    if "spectatorAngleNodes" in rec:
        parts.append(f"spectatorAngleNodes={int(rec['spectatorAngleNodes'])}")
    if "cbackground" in rec:
        parts.append(f"cbackground={rec['cbackground']}")
    return " ".join(parts)


def _thm_plain_number(x):
    """Shortest round-trip text of a number, without exponent for whole
    numbers ('100', '4.2', '1e-05')."""
    t = repr(float(x))
    return t[:-2] if t.endswith(".0") else t


def _thm_is_default(s):
    return s == _thm_default_settings()


def _thm_number(x):
    """Shortest text that reads back as the same double, as Qt's
    QString::number(x, 'g', QLocale::FloatingPointShortest) writes it."""
    for p in range(1, 18):
        t = f"{x:.{p}g}"
        if float(t) == x:
            return t
    return repr(x)


def _thm_read_double(text):
    """operator>> on an istringstream: a leading number, the rest ignored."""
    m = _THM_FLOAT_PREFIX.match(text)
    return float(m.group(0)) if m else None


def _thm_read_int(text, strict):
    m = _THM_INT_PREFIX.match(text)
    if not m:
        return None
    if strict and text[m.end():].strip():
        return None
    return int(m.group(0))


def _thm_parse_line(raw, s):
    """One block line into settings ``s``: returns (key, value) normalized
    as the GUI writes them, ("", "") for a blank or comment line; raises
    ValueError for a line the engine would refuse."""
    line = raw.split("#", 1)[0].strip()
    if not line:
        return "", ""
    if line.startswith("experiment["):
        # An experiment record: validated and merged, but kept as a line of
        # its own (like a comment) when the options are rewritten.
        _thm_parse_experiment(line, s["experiments"])
        return "", ""
    if "=" not in line:
        raise ValueError(f"<thm> line not understood: '{raw.strip()}'")
    k, v = line.split("=", 1)
    k = k.rstrip(" \t")
    v = v.lstrip(" \t")
    bad = ValueError(f"<thm> line not understood: '{raw.strip()}'")
    if k == "vertex":
        if v == "real":
            v = "perlevel"
        if v not in ("onshell", "constant", "perlevel"):
            raise bad
        s["vertex"] = v
    elif k == "kinematics":
        if v not in ("lacognata", "triple", "kf3body", "lambda32"):
            raise bad
        s["kinematics"] = v
    elif k == "entranceL":
        if v not in ("coherent", "incoherent"):
            raise bad
        s["entranceL"] = v
    elif k == "coulombIntegral":
        if v in ("1", "true", "on"):
            s["coulombIntegral"] = True
        elif v in ("0", "false", "off"):
            s["coulombIntegral"] = False
        else:
            raise bad
        v = "1" if s["coulombIntegral"] else "0"
    elif k.startswith("spectatorEnergy"):
        x = _thm_read_double(v)
        if x is None or not x >= 0.0:
            raise bad
        v = _thm_number(x)
        if k == "spectatorEnergy":
            s["spectatorEnergy"] = x
        elif len(k) > 17 and k[15] == "[" and k.endswith("]"):
            pair = _thm_read_int(k[16:-1], strict=False)
            if pair is None:
                raise bad
            s["spectatorByPair"][pair] = x
            k = f"spectatorEnergy[{pair}]"
        else:
            raise bad
    elif k.startswith("weight"):
        test = k.startswith("weightTest[")
        open_ = 10 if test else 6
        if not (len(k) > open_ + 2 and k[open_] == "[" and k.endswith("]") and v):
            raise bad
        seg = _thm_read_int(k[open_ + 1:-1], strict=True)
        if seg is None or seg < 1:
            raise bad
        s["weightTest" if test else "weight"][seg] = v
        k = f"{'weightTest' if test else 'weight'}[{seg}]"
    else:
        raise bad
    return k, v


def _thm_key_values(s):
    """ThmSettings::keyValues: the non-default options, in the GUI's order."""
    d = _thm_default_settings()
    kv = []
    for key in ("entranceL", "vertex", "kinematics"):
        if s[key] != d[key]:
            kv.append((key, s[key]))
    if s["coulombIntegral"] != d["coulombIntegral"]:
        kv.append(("coulombIntegral", "1"))
    if s["spectatorEnergy"] != d["spectatorEnergy"]:
        kv.append(("spectatorEnergy", _thm_number(s["spectatorEnergy"])))
    for k in sorted(s["spectatorByPair"]):
        kv.append((f"spectatorEnergy[{k}]", _thm_number(s["spectatorByPair"][k])))
    for k in sorted(s["weight"]):
        kv.append((f"weight[{k}]", s["weight"][k]))
    for k in sorted(s["weightTest"]):
        kv.append((f"weightTest[{k}]", s["weightTest"][k]))
    return kv


def _thm_compose(old_lines, s):
    """ThmSettings::compose: the block body for settings ``s``, keeping
    comments, blank lines and unchanged lines of ``old_lines``."""
    kv = _thm_key_values(s)
    wanted = dict(kv)
    written, out = set(), []
    for line in old_lines:
        try:
            key, value = _thm_parse_line(line, _thm_default_settings())
        except ValueError:
            key = ""
        if not key:
            out.append(line)           # comment or blank line
            continue
        if key not in wanted or key in written:
            continue                   # back to default, or repeated
        written.add(key)
        if value == wanted[key]:
            out.append(line)
            continue
        # Rewrite the value, keeping the indentation and any inline comment.
        hash_ = line.find("#")
        code = line[:hash_] if hash_ >= 0 else line
        comment = line[hash_:] if hash_ >= 0 else ""
        start = len(code) - len(code.lstrip())
        end = len(code.rstrip())
        out.append(code[:start] + key + "=" + wanted[key] + code[end:] + comment)
    for key, value in kv:
        if key not in written:
            out.append(key + "=" + value)
    return out


def _thm_read_weight_table(path):
    """ThmWeightTable::Read: "" if the table is usable, else why not."""
    try:
        with open(path) as f:
            text = f.read()
    except OSError:
        return f"cannot read the weight file '{path}'"
    energies = []
    for n, line in enumerate(text.splitlines(), 1):
        line = line.split("#", 1)[0]
        if not line.strip():
            continue
        where = f"'{path}' line {n}: "
        tok = line.split()
        try:
            if len(tok) != 2 or any("_" in t for t in tok):
                raise ValueError
            e, w = float(tok[0]), float(tok[1])
        except ValueError:
            return where + "expected two numbers, E (MeV) and w"
        if not (math.isfinite(e) and math.isfinite(w) and w > 0.0):
            return where + "the weight must be finite and > 0"
        if energies and not e > energies[-1]:
            return where + "the energies must be strictly increasing"
        energies.append(e)
    if len(energies) < 2:
        return f"'{path}' needs at least two rows (E w)"
    return ""


class AzrChannel:
    """One channel line (31 fields).  Level fields are shared by a level's lines.

    Keeps the line exactly as it was read, so an untouched channel round-trips
    byte for byte; typed accessors parse on demand and mutators rewrite only the
    field they change, in place, keeping the column the file had it in.
    """

    def __init__(self, tokens, raw=None):
        if not (_NFIELDS <= len(tokens) <= _NFIELDS_MAX):
            raise ValueError(
                f"a <levels> line needs {_NFIELDS} to {_NFIELDS_MAX} fields, "
                f"got {len(tokens)}: {tokens}")
        self.tokens = list(tokens)
        #: The line as read, or None for a channel built from tokens alone.
        self._raw = raw
        #: Field names whose token no longer matches ``_raw``.
        self._dirty = set()
        # Remember whether the line carried the optional (THM) fields, so a
        # file that did not have them is written back without them.
        self._n_written = len(tokens)
        while len(self.tokens) < _NFIELDS_MAX:
            self.tokens.append(_OPTIONAL_DEFAULTS[len(self.tokens) - _NFIELDS])

    # -- typed field access ---------------------------------------------------

    def _get(self, name, cast):
        return cast(self.tokens[_IDX[name]])

    def _set(self, name, value):
        formatted = _fmt(value)
        if self.tokens[_IDX[name]] == formatted:
            return                      # unchanged: do not disturb the raw line
        self.tokens[_IDX[name]] = formatted
        self._dirty.add(name)

    # channel identity

    @property
    def pair(self):
        """1-based particle-pair number this channel decays to."""
        return self._get("ir", lambda v: int(float(v)))

    @property
    def entrance_key(self):
        """Entrance-pair key stored on the line."""
        return self._get("aa", lambda v: int(float(v)))

    @property
    def L(self):
        """Orbital angular momentum of the channel."""
        return self._get("l2", lambda v: int(float(v))) // 2

    @L.setter
    def L(self, v):        self._set("l2", int(2 * v))

    @property
    def S(self):
        """Channel spin."""
        return self._get("s2", lambda v: int(float(v))) / 2.0

    @S.setter
    def S(self, v):        self._set("s2", int(round(2 * v)))

    @property
    def gamma(self):
        """The <levels> width field: a partial width in eV for an open channel, an ANC in fm^-1/2 for a closed one -- not a reduced-width amplitude."""
        return self._get("gamma", float)

    @gamma.setter
    def gamma(self, v):    self._set("gamma", float(v))

    @property
    def channel_fixed(self):
        """Is this channel's width held fixed in the fit?"""
        return self._get("channelFix", lambda v: int(float(v))) != 0

    @channel_fixed.setter
    def channel_fixed(self, v): self._set("channelFix", 1 if v else 0)

    @property
    def ptype(self):
        """Particle type code: 0 for a particle channel, nonzero for a photon."""
        return self._get("pType", lambda v: int(float(v)))

    @property
    def is_photon(self):
        """Is this a photon channel?"""
        return self.ptype != 0

    @property
    def channel_radius(self):
        """Channel radius in fm."""
        return self._get("chRad", float)

    @channel_radius.setter
    def channel_radius(self, v): self._set("chRad", float(v))

    @property
    def active(self):
        """Is the level active?"""
        return self._get("isActive", lambda v: int(float(v))) != 0

    # pair physics -- the same quantities GET_PAIRS_INFO reports at runtime, but
    # readable straight from the file, so a caller can identify a model's
    # channels without launching AZURE2.

    @property
    def Z1(self):
        """Charge number of the light particle."""
        return self._get("z1", lambda v: int(float(v)))

    @property
    def Z2(self):
        """Charge number of the heavy particle."""
        return self._get("z2", lambda v: int(float(v)))

    @property
    def M1(self):
        """Mass of the light particle, in u."""
        return self._get("m1", float)

    @property
    def M2(self):
        """Mass of the heavy particle, in u."""
        return self._get("m2", float)

    @property
    def J1(self):
        """Intrinsic spin of the light particle."""
        return self._get("j1", float)

    @property
    def parity1(self):
        """Parity of the light particle."""
        return self._get("pi1", lambda v: int(float(v)))

    @property
    def J2(self):
        """Intrinsic spin of the heavy particle."""
        return self._get("j2", float)

    @property
    def parity2(self):
        """Parity of the heavy particle."""
        return self._get("pi2", lambda v: int(float(v)))

    @property
    def excitation(self):
        """Excitation energy of the pair's residual nucleus (MeV).

        Zero for a ground-state pair; this is what distinguishes the capture
        channels of a multi-transition model (gamma_0, gamma_1, ...).
        """
        return self._get("e2", float)

    @property
    def sep_energy(self):
        """Separation energy of the pair, in MeV."""
        return self._get("sepE", float)

    # level fields (shared across a level's channel lines)

    @property
    def levelJ(self):
        """Total angular momentum of the level."""
        return self._get("levelJ", float)

    @property
    def levelPi(self):
        """Parity of the level."""
        return self._get("levelPi", lambda v: int(float(v)))

    @property
    def levelE(self):
        """Level energy in MeV -- an excitation energy of the compound nucleus."""
        return self._get("levelE", float)

    @property
    def levelID(self):
        """The level's number, shared by all of its channel lines."""
        return self._get("levelID", lambda v: int(float(v)))

    @property
    def level_fixed(self):
        """Is the level energy held fixed in the fit?"""
        return self._get("levelFix", lambda v: int(float(v))) != 0

    def _set_level(self, J, parity, energy, level_fixed, level_id):
        self._set("levelJ", float(J))
        self._set("levelPi", int(parity))
        self._set("levelE", float(energy))
        self._set("levelFix", 1 if level_fixed else 0)
        self._set("levelID", int(level_id))

    def clone(self):
        """A copy of this channel line, sharing nothing with the original.

        The raw line comes with it, so a cloned channel whose couplings are then
        changed keeps the column layout of the one it was cloned from.
        """
        copy = AzrChannel(self.tokens[:self._n_written], self._raw)
        copy._dirty = set(self._dirty)
        return copy

    def to_line(self):
        """The channel as one line of the ``<levels>`` block.

        An untouched line is returned exactly as it was read.  An edited one is
        the same line with the changed fields substituted in place, right
        aligned in the width the file used, so a single edit does not reflow the
        whole block.  Two fallbacks: a channel with no raw line behind it --
        built from tokens, not parsed -- joins with single spaces; and an edit
        to an optional (THM) field a shorter line never carried appends the
        optional columns, since there is no column to substitute into.
        """
        if self._raw is None:
            return " ".join(self.tokens[:self._n_written])
        if not self._dirty:
            return self._raw
        spans = _token_spans(self._raw)
        if len(spans) < self._n_written:  # not the line we parsed; do not guess
            return " ".join(self.tokens[:self._n_written])
        out = self._raw
        appended = self._n_written > len(spans)
        for name in sorted(self._dirty, key=lambda n: -_IDX[n]):
            i = _IDX[name]
            if i >= len(spans):
                appended = True
                continue                  # no column in the raw line: append below
            start, end = spans[i]
            new = self.tokens[i]
            width = end - start
            out = out[:start] + (new.rjust(width) if len(new) <= width else new) + out[end:]
        if appended:
            out = out.rstrip() + "  " + "  ".join(
                self.tokens[len(spans):self._n_written])
        return out

    # -- the optional trailing fields ----------------------------------------
    @property
    def gamma_is_rwa(self):
        """True when the ``gamma`` column is a reduced width amplitude
        (MeV^(1/2)) rather than a partial width in eV or an ANC.

        AZURE2 keeps this convention on output as well: ``CNuc::TransformOut``
        returns the amplitude unchanged for such a channel, so a fit written
        back through :meth:`AzrModel.apply_fit` does not flip the convention.
        """
        return self._get("gammaIsRWA", lambda v: int(float(v))) == 1

    @gamma_is_rwa.setter
    def gamma_is_rwa(self, v):
        self._set("gammaIsRWA", 1 if v else 0)
        self._n_written = _NFIELDS_MAX      # the field must now be emitted

    @property
    def binding_energy(self):
        """THM binding energy B_xs (MeV) of the transferred particle; 0 for a
        conventional pair."""
        return self._get("bindingEnergy", float)

    @binding_energy.setter
    def binding_energy(self, v):
        self._set("bindingEnergy", float(v))
        self._n_written = _NFIELDS_MAX      # the field must now be emitted

    @classmethod
    def from_line(cls, line):
        """Parse one .azr channel line, keeping it verbatim for round-tripping."""
        return cls(line.split(), raw=line)


class AzrLevel:
    """A level: shared (J, parity, energy, fixed) + its channels."""

    def __init__(self, channels: List[AzrChannel]):
        if not channels:
            raise ValueError("a level needs at least one channel.")
        self.channels = list(channels)

    @property
    def J(self):
        """Total angular momentum of the level."""
        return self.channels[0].levelJ

    @property
    def parity(self):
        """Parity of the level."""
        return self.channels[0].levelPi

    @property
    def energy(self):
        """Level energy in MeV."""
        return self.channels[0].levelE

    @property
    def level_id(self):
        """The level's number."""
        return self.channels[0].levelID

    @property
    def fixed(self):
        """Is the level energy held fixed?"""
        return self.channels[0].level_fixed

    @property
    def jpi(self):
        """J^pi as text, e.g. "3/2-"."""
        j = int(self.J) if float(self.J).is_integer() else f"{int(round(2*self.J))}/2"
        return f"{j}{'+' if self.parity > 0 else '-'}"

    @property
    def active(self):
        """Is the level active (the ``isActive`` field of its lines)?

        AZURE2 skips every ``<levels>`` line whose ``isActive`` is 0
        (``CNuc::Fill``), so an inactive level is invisible to the engine: it
        gets no parameters and does not count in the engine's level numbering.
        The flag is per line, but a level whose lines disagree is not something
        the engine can represent consistently, so that raises here.
        """
        flags = {c.active for c in self.channels}
        if len(flags) != 1:
            raise ValueError(
                f"level {self.jpi} at {self.energy} MeV has both active and "
                f"inactive channel lines; AZURE2 reads them line by line, so "
                f"the level would be half present.")
        return flags.pop()

    def set_active(self, active):
        """Set the ``isActive`` flag on every channel line of the level."""
        for c in self.channels:
            c._set("isActive", 1 if active else 0)

    def set_energy(self, energy):
        """Set the level energy on every channel line of the level."""
        for c in self.channels:
            c._set("levelE", float(energy))

    def set_fixed(self, fixed):
        """Fix or free the level energy on every channel line."""
        for c in self.channels:
            c._set("levelFix", 1 if fixed else 0)

    def _renumber(self, level_id):
        for c in self.channels:
            c._set("levelID", int(level_id))

    def __repr__(self):
        return (f"AzrLevel(J^pi={self.jpi}, E={self.energy:g} MeV, "
                f"{len(self.channels)} channels"
                f"{', fixed' if self.fixed else ''})")


class AzrModel:
    """A ``.azr`` file parsed into an editable level scheme.

    Only the ``<levels>`` block is interpreted; the surrounding text is kept
    verbatim and re-emitted unchanged by :meth:`write`.
    """

    def __init__(self, prefix, levels, suffix, source=None):
        self._prefix = prefix       # text up to and including "<levels>\n"
        self.levels: List[AzrLevel] = levels
        self._suffix = suffix       # text from "</levels>" onward
        self.source = source

    # -- parsing / serialization ---------------------------------------------

    @classmethod
    def from_file(cls, path):
        """Parse a .azr file; only <levels> is interpreted, the rest is kept verbatim."""
        with open(path) as f:
            text = f.read()
        start = text.find("<levels>")
        end = text.find("</levels>")
        if start == -1 or end == -1:
            raise ValueError(f"{path}: no <levels> ... </levels> block.")
        body_start = text.index("\n", start) + 1
        prefix = text[:body_start]
        body = text[body_start:end]
        suffix = text[end:]

        # group channel lines into levels by their levelID
        levels, cur, cur_id = [], [], None
        for line in body.splitlines():
            if not line.strip():
                continue
            ch = AzrChannel.from_line(line)
            if cur and ch.levelID != cur_id:
                levels.append(AzrLevel(cur))
                cur = []
            cur.append(ch)
            cur_id = ch.levelID
        if cur:
            levels.append(AzrLevel(cur))
        return cls(prefix, levels, suffix, source=path)

    def _final_newline(self):
        """"\n" if the file ends with a newline, else "" -- what it keeps."""
        return "\n" if self._suffix.endswith("\n") else ""

    def _set_suffix_lines(self, lines):
        """Replace the verbatim tail by ``lines``, keeping whether the file
        ends with a newline.  (``"\n".join(text.splitlines())`` alone drops
        it, so every edit of a block after ``<levels>`` used to leave the
        written file without its final newline.)"""
        self._suffix = "\n".join(lines) + self._final_newline()

    def set_output_dir(self, path):
        """Point the model's ``<config>`` output directory somewhere else.

        AZURE2 writes its results *and* its external-capture integral caches
        (``intEC.dat`` / ``intEC.extrap``) here.  Giving an edited model its own
        directory keeps it from overwriting the main run's outputs -- and, since
        the cache is only valid for the grid it was built on, keeps a variant's
        integrals from being read back into a run with different segments.

        The directory is created if it does not exist; AZURE2 does not create it
        and fails quietly otherwise.
        """
        path = str(path)
        if not path.endswith("/"):
            path += "/"
        os.makedirs(path, exist_ok=True)
        lines = self._prefix.splitlines(keepends=True)
        for i, line in enumerate(lines):
            if line.rstrip().endswith("#Full Path to Output Directory"):
                comment = line[line.index("#"):]
                lines[i] = f"{path:<100}{comment}"
                break
        else:
            raise ValueError("no output-directory line in the <config> block.")
        self._prefix = "".join(lines)
        return self

    def to_text(self):
        """The whole file as text, with <levels> re-emitted and everything else unchanged."""
        self._renumber()
        blocks = ["\n".join(c.to_line() for c in lv.channels)
                  for lv in self.levels]
        body = "\n\n".join(blocks)
        return self._prefix + body + "\n\n" + self._suffix

    def write(self, path):
        """Write the (edited) model to ``path`` and return it.  The original
        file is never modified."""
        with open(path, "w") as f:
            f.write(self.to_text())
        return path

    def to_tempfile(self, suffix=".azr", dir=None):
        """Write to a fresh temp file (for launching an edited instance) and
        return its path.  The caller owns the file."""
        fd, path = tempfile.mkstemp(suffix=suffix, dir=dir)
        os.close(fd)
        return self.write(path)

    # -- queries --------------------------------------------------------------

    def find(self, jpi=None, energy=None, tol=1e-3, index=None):
        """Return the levels matching ``jpi`` and/or ``energy`` (or by index)."""
        if index is not None:
            return [self.levels[index]]
        out = []
        for lv in self.levels:
            if jpi is not None and lv.jpi != jpi:
                continue
            if energy is not None and abs(lv.energy - energy) > tol:
                continue
            out.append(lv)
        return out

    def _pair_template(self, pair):
        """An existing channel that uses ``pair`` (for cloning its pair data)."""
        for lv in self.levels:
            for c in lv.channels:
                if c.pair == pair:
                    return c
        raise ValueError(
            f"pair {pair} is not used by any existing channel; adding a brand-"
            f"new particle pair is not supported -- reference an existing pair.")

    def _renumber(self):
        """Give every level a fresh consecutive levelID (1-based)."""
        for i, lv in enumerate(self.levels, start=1):
            lv._renumber(i)

    # -- editing --------------------------------------------------------------

    def remove_level(self, jpi=None, energy=None, index=None, tol=1e-3):
        """Remove the matching level(s).  Returns the removed :class:`AzrLevel`
        objects.  Raises if the selector matches nothing."""
        victims = self.find(jpi=jpi, energy=energy, index=index, tol=tol)
        if not victims:
            raise KeyError(f"no level matches jpi={jpi} energy={energy} "
                           f"index={index}.")
        self.levels = [lv for lv in self.levels if lv not in victims]
        self._renumber()
        return victims

    def add_level(self, J, parity, energy, channels, level_fixed=True,
                  at=None):
        """Add a level with the given channels.

        ``channels`` is a list of dicts, each ``{pair, L, S, gamma[, fixed]}``,
        referencing an existing particle pair by number; the pair's physical
        data (masses, charges, separation energy, radius) is copied from an
        existing channel that uses that pair.

        **All levels of a given J^pi share one channel set** (an R-matrix
        J-group requirement).  If a level of this ``(J, parity)`` already exists,
        this method clones that group's exact channel structure and only applies
        your ``gamma`` / ``fixed`` values, matched by ``(pair, L, S)`` -- any
        channel you don't mention is added inert (``gamma=0``, fixed).  A spec
        that matches none of the group's channels raises, listing the allowed
        ones.  For a brand-new J^pi the channels you give define the group.

        Returns the new :class:`AzrLevel`.
        """
        new_id = max((lv.level_id for lv in self.levels), default=0) + 1
        existing = [lv for lv in self.levels
                    if lv.J == J and lv.parity == parity]

        if existing:
            # Clone the group's channel structure; start every channel inert.
            template = existing[0]
            chans = [c.clone() for c in template.channels]
            for ch in chans:
                ch._set_level(J, parity, energy, level_fixed, new_id)
                ch._set("isActive", 1)
                ch.gamma = 0.0
                ch.channel_fixed = True
            allowed = [(ch.pair, ch.L, ch.S) for ch in chans]
            for spec in channels:
                match = [ch for ch in chans
                         if ch.pair == spec["pair"] and ch.L == spec["L"]
                         and abs(ch.S - spec["S"]) < 1e-6]
                if not match:
                    raise ValueError(
                        f"J^pi {template.jpi} already exists with a fixed "
                        f"channel set {allowed}; the spec "
                        f"(pair={spec['pair']}, L={spec['L']}, S={spec['S']}) "
                        f"matches none of them.")
                match[0].gamma = spec.get("gamma", 0.0)
                match[0].channel_fixed = spec.get("fixed", False)
        else:
            chans = []
            for spec in channels:
                ch = self._pair_template(spec["pair"]).clone()
                ch._set_level(J, parity, energy, level_fixed, new_id)
                ch._set("isActive", 1)
                ch.L = spec["L"]
                ch.S = spec["S"]
                ch.gamma = spec.get("gamma", 0.0)
                ch.channel_fixed = spec.get("fixed", False)
                chans.append(ch)

        level = AzrLevel(chans)
        if at is None:
            self.levels.append(level)
        else:
            self.levels.insert(at, level)
        self._renumber()
        return level

    def remove_channel(self, jpi, pair, L, S, tol=1e-6):
        """Remove a channel ``(pair, L, S)`` from *every* level of a J^pi group.

        All levels of a J-group share one channel set, so a channel is removed
        from the whole group at once (removing it from a single level would make
        AZURE2 reject the file).  Returns the number of channel lines removed.
        Raises if no level of ``jpi`` has that channel, or if it is a level's
        only channel.
        """
        removed = 0
        for lv in self.levels:
            if lv.jpi != jpi:
                continue
            keep = [c for c in lv.channels
                    if not (c.pair == pair and c.L == L
                            and abs(c.S - S) < tol)]
            if len(keep) == len(lv.channels):
                continue
            if not keep:
                raise ValueError(
                    f"channel (pair={pair}, L={L}, S={S}) is the only channel "
                    f"of a {jpi} level; cannot remove it.")
            removed += len(lv.channels) - len(keep)
            lv.channels = keep
        if removed == 0:
            raise KeyError(f"no {jpi} level has channel "
                           f"(pair={pair}, L={L}, S={S}).")
        return removed

    def add_channel(self, level, pair, L, S, gamma=0.0, fixed=False):
        """Add one channel to an existing level (referencing an existing pair)."""
        template = self._pair_template(pair)
        ch = template.clone()
        ch._set_level(level.J, level.parity, level.energy, level.fixed,
                      level.level_id)
        ch._set("isActive", 1)
        ch.L, ch.S, ch.gamma = L, S, gamma
        ch.channel_fixed = fixed
        level.channels.append(ch)
        return ch

    # -- channel radius -------------------------------------------------------

    def channel_radii(self):
        """``{pair: radius}`` in fm, for every particle pair in the file."""
        out = {}
        for lv in self.levels:
            for c in lv.channels:
                if not c.is_photon:
                    out[c.pair] = c.channel_radius
        return dict(sorted(out.items()))

    def set_channel_radius(self, pair, radius):
        """Set the channel radius (fm) of one particle pair, on every line.

        AZURE2 stores the radius per channel line, so it must be written to all
        of them or the model is inconsistent; this does that and returns the
        number of lines changed.

        The radius is the matching surface between the internal R-matrix region
        and the external Coulomb solutions, so changing it changes the
        penetrabilities, shift functions, boundary conditions, Wigner limits and
        the lower limit of every external-capture integral. A reduced width
        therefore *means* something different afterwards, and the level scheme
        is no longer fitted -- **refit before reading anything off the result**.
        Delete ``output/intEC.dat`` / ``output/intEC.extrap`` (or write the new
        model into its own output directory) so the cached external-capture
        integrals, which belong to the old radius, cannot be reused.

        >>> mdl = AzrModel.from_file("7Be.azr")
        >>> mdl.channel_radii()                 # {1: 4.24151, 2: 3.94396, 3: 3.94396}
        >>> mdl.set_channel_radius(1, 5.0)      # 3He+alpha
        >>> path = mdl.write("a5.0.azr")
        """
        radius = float(radius)
        if not radius > 0:
            raise ValueError(f"channel radius must be positive, got {radius}.")
        pair = int(pair)
        changed = 0
        for lv in self.levels:
            for c in lv.channels:
                if c.pair == pair and not c.is_photon:
                    c.channel_radius = radius
                    changed += 1
        if not changed:
            raise KeyError(
                f"pair {pair} has no particle channel in this model "
                f"(radii: {self.channel_radii()}).")
        return changed

    # -- level activation -----------------------------------------------------

    def deactivate_level(self, jpi=None, energy=None, index=None, tol=1e-3):
        """Switch off the matching level(s) without removing them.

        Every channel of a matched level has its reduced width (``gamma``) set
        to zero and is fixed, so the level stays in the file (and in its
        J-group's shared channel set) but couples to nothing -- it contributes
        exactly nothing to the calculation.  This is the file-level counterpart
        of :meth:`azure2.without_level`; use it when you want to persist a
        "level removed" model and refit it.  Returns the affected levels.
        """
        victims = self.find(jpi=jpi, energy=energy, index=index, tol=tol)
        if not victims:
            raise KeyError(f"no level matches jpi={jpi} energy={energy} "
                           f"index={index}.")
        for lv in victims:
            for c in lv.channels:
                c.gamma = 0.0
                c.channel_fixed = True
        return victims

    # -- extrapolations (edits the <segmentsTest> block) ----------------------

    # observable name -> isDiff code for a <segmentsTest> line
    # (ESegment::ESegment(ExtrapLine); mirrors datasets._EXTRAP_OBSERVABLE).
    _EXTRAP_CODE = {
        "angle-integrated": 0, "differential": 1, "phase-shift": 2,
        "angular-distribution": 3, "total-capture": 4, "differential-cm": 5,
        # Vector analyzing power A_y. Differential in the centre-of-mass frame
        # like code 5, but the quantity is a dimensionless ratio, not a cross
        # section.
        "analyzing-power": 7,
    }

    def _splice_segments_test(self, new_lines):
        """Replace the body of the <segmentsTest> block in the verbatim tail.

        ``new_lines`` is a list of pre-formatted extrapolation lines (no
        trailing newline).  Creates the block if the file has none.
        """
        lines = self._suffix.splitlines()
        try:
            start = lines.index("<segmentsTest>")
            end = lines.index("</segmentsTest>")
        except ValueError:
            # No block yet: append one at the end of the tail.
            self._suffix = (self._suffix.rstrip("\n") + "\n\n<segmentsTest>\n"
                            + "\n".join(new_lines) + "\n</segmentsTest>"
                            + self._final_newline())
            return
        self._set_suffix_lines(lines[:start + 1] + new_lines
                                 + lines[end:])

    def clear_extrapolations(self):
        """Remove every ``<segmentsTest>`` line (leave the block empty)."""
        self._splice_segments_test([])
        return self

    def add_extrapolation(self, entrance, exit, e_min, e_max, e_step,
                          observable="angle-integrated", angle=None,
                          angle_min=0.0, angle_max=0.0, angle_step=0.0,
                          phase_J=None, phase_L=None, order=None, active=True,
                          frame="lab"):
        """Append one extrapolation segment (a grid to evaluate the model on).

        ``entrance`` / ``exit`` are particle-pair keys (``exit=-1`` for a summed
        / total observable such as capture).  The model is evaluated on the
        energy grid ``e_min : e_max : e_step`` in MeV.

        ``frame`` says what those energies are.  ``"lab"`` (the default) is
        what a ``<segmentsTest>`` line holds and AZURE2 reads: the lab energy
        of the light particle of the entrance pair on the heavy one at rest,
        as for data files and in the GUI.  ``"cm"`` gives entrance-channel
        c.m. energies, converted here to lab with the pair's masses,
        E_lab = E_cm (M1 + M2) / M2 -- the inverse of AZURE2's own conversion,
        so the engine evaluates, and :meth:`~pyazr.azure2.azure2.calculate_energies`
        reports, the c.m. grid asked for.  (Until pyazr 2.8.0 this docstring called
        the grid c.m. while the engine read it as lab: a c.m. grid came out
        stretched by (M1 + M2)/M2 in c.m. energy, 6 % for p + 17O.)

        ``observable`` is one of ``angle-integrated``, ``differential``,
        ``differential-cm``, ``total-capture``, ``angular-distribution``,
        ``phase-shift``.  For a single-angle differential give ``angle=`` (sets
        min=max, step 0); for an angular grid give ``angle_min/max/step``.
        ``phase_J``/``phase_L`` are required for a phase-shift, ``order`` for an
        angular distribution.
        """
        if observable not in self._EXTRAP_CODE:
            raise ValueError(f"unknown observable {observable!r}; expected one "
                             f"of {sorted(self._EXTRAP_CODE)}.")
        if frame not in ("lab", "cm"):
            raise ValueError(f"frame must be 'lab' or 'cm', not {frame!r}.")
        if frame == "cm":
            pair = self._pair_template(int(entrance))
            if pair.ptype != 0:
                raise ValueError(
                    f"frame='cm' converts with a particle pair's masses; entrance "
                    f"pair {entrance} has pType {pair.ptype} -- give its energies "
                    f"as AZURE2 reads them (frame='lab').")
            to_lab = (pair.M1 + pair.M2) / pair.M2
            e_min, e_max, e_step = e_min * to_lab, e_max * to_lab, e_step * to_lab
            # AZURE2 steps e_min + e_step + ... while <= e_max; the converted
            # numbers are not exact, so a last point that should sit on e_max
            # could land just above it and be lost.  A billionth of a step of
            # slack keeps it.
            e_max += 1e-9 * abs(e_step)
        isDiff = self._EXTRAP_CODE[observable]
        if angle is not None:
            angle_min = angle_max = angle
            angle_step = 0.0
        toks = [1 if active else 0, int(entrance), int(exit),
                _fmt(e_min), _fmt(e_max), _fmt(e_step),
                _fmt(angle_min), _fmt(angle_max), _fmt(angle_step), isDiff]
        if isDiff == 2:                       # phase shift carries J, L
            if phase_J is None or phase_L is None:
                raise ValueError("a phase-shift extrapolation needs phase_J "
                                 "and phase_L.")
            toks += [_fmt(phase_J), int(phase_L)]
        elif isDiff == 3:                     # angular distribution carries order
            if order is None:
                raise ValueError("an angular-distribution extrapolation needs "
                                 "order.")
            toks += [int(order)]
        toks += [0]                           # trailing isAdvanced flag
        line = "  ".join(_fmt(t) if not isinstance(t, str) else t for t in toks)
        lines = self._suffix.splitlines()
        if "<segmentsTest>" in lines:
            end = lines.index("</segmentsTest>")
            self._set_suffix_lines(lines[:end] + [line] + lines[end:])
        else:
            self._splice_segments_test([line])
        return self

    def set_extrapolations(self, specs):
        """Replace the whole ``<segmentsTest>`` block with ``specs``.

        ``specs`` is an iterable of dicts, each a keyword bundle for
        :meth:`add_extrapolation` (e.g. ``dict(entrance=1, exit=-1, e_min=0.05,
        e_max=5, e_step=0.05, observable="total-capture")``).
        """
        self.clear_extrapolations()
        for spec in specs:
            self.add_extrapolation(**spec)
        return self

    def keep_extrapolations(self, keys):
        """Drop every ``<segmentsTest>`` line except the given 1-based ``keys``,
        in the order they are listed here.

        The surviving lines are copied verbatim, so the grids are exactly the
        ones the file declares.  AZURE2 re-evaluates *every* active segment on
        each forward pass, so trimming the block to the handful you actually
        want is what makes a finite-difference uncertainty band affordable --
        see :func:`pyazr.bands.extrapolation_bands`.

        Note this invalidates ``output/intEC.extrap``: AZURE2 caches the
        external-capture integrals there and reuses them on a changed grid,
        silently corrupting the result.  Delete it before running the edited
        model.
        """
        lines = self._suffix.splitlines()
        try:
            start = lines.index("<segmentsTest>") + 1
            end = lines.index("</segmentsTest>")
        except ValueError:
            raise ValueError("no <segmentsTest> block to trim.")
        body = [ln for ln in lines[start:end] if ln.strip()]
        missing = [k for k in keys if not 1 <= k <= len(body)]
        if missing:
            raise KeyError(f"<segmentsTest> has {len(body)} segments; no "
                           f"{missing}.")
        self._splice_segments_test([body[k - 1] for k in keys])
        return self

    # -- segment normalizations (edits the <segmentsData> block) --------------

    def set_segment_norm(self, file_substr, vary=None, sys_error=None):
        """Edit the normalization of every ``<segmentsData>`` line whose text
        matches ``file_substr`` (e.g. a data-file name).

        ``vary=True/False`` frees/fixes the normalization; ``sys_error`` sets its
        systematic (percent, as stored in the file).  Returns the number of
        segment lines changed.  Operates on the verbatim tail text, so the
        levels editing is unaffected.
        """
        if "<segmentsData>" not in self._suffix:
            raise ValueError("no <segmentsData> block to edit.")
        out, changed, inside = [], 0, False
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
            elif s == "</segmentsData>":
                inside = False
            elif inside and s and file_substr in line:
                t = line.split()
                isDiff = int(float(t[7]))
                i = 8 + (2 if isDiff % 10 == 2 else 0)   # dataNorm index
                if vary is not None:
                    t[i + 1] = "1" if vary else "0"
                if sys_error is not None:
                    t[i + 2] = _fmt(sys_error)
                line = " ".join(t)
                changed += 1
            out.append(line)
        self._set_suffix_lines(out)
        if changed == 0:
            raise KeyError(f"no <segmentsData> line matches {file_substr!r}.")
        return changed

    def segment_values(self):
        """``{key: (norm, shift)}`` of every ``<segmentsData>`` line, ``key``
        the engine's segment key (1-based line number, inactive lines
        counted: the ``segment_<key>_norm`` names of the parameters); shift
        is None on a line without the energy-shift fields."""
        out = {}
        for key, (t, i) in self._segment_lines().items():
            shift = float(t[i + 3]) if len(t) > i + 3 and _isnum(t[i + 3]) else None
            out[key] = (float(t[i]), shift)
        return out

    def _segment_lines(self):
        """``{key: (tokens, norm index)}`` of the ``<segmentsData>`` lines."""
        out, inside, key = {}, False, 0
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
            elif s == "</segmentsData>":
                inside = False
            elif inside and s:
                key += 1
                t = line.split()
                isDiff = int(float(t[7]))
                out[key] = (t, 8 + (2 if isDiff % 10 == 2 else 0))
        return out

    def set_segment_values(self, values):
        """Write normalizations and energy shifts into ``<segmentsData>``.

        ``values`` maps the engine's segment key (see :meth:`segment_values`)
        to ``(norm, shift)``; either may be None to leave it.  Only those two
        fields of those lines change (the rest of each line, spacing
        included, is kept).  **The norm field is also the centre of the
        segment's normalization prior** (and of its percentage error) when a
        fit starts from this file, and the shift field the centre of the
        shift prior: writing fitted values here moves those centres -- see
        :meth:`pyazr.azure2.azure2.save_fit`.  Returns the number of lines
        changed; KeyError for a key with no line.
        """
        lines = self._suffix.splitlines()
        inside, key, changed, seen = False, 0, 0, set()
        for n, line in enumerate(lines):
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
                continue
            if s == "</segmentsData>":
                inside = False
                continue
            if not (inside and s):
                continue
            key += 1
            if key not in values:
                continue
            seen.add(key)
            norm, shift = values[key]
            spans = _token_spans(line)
            t = [line[a:b] for a, b in spans]
            i = 8 + (2 if int(float(t[7])) % 10 == 2 else 0)
            edits = {}
            if norm is not None:
                edits[i] = repr(float(norm))
            if shift is not None:
                if not (len(t) > i + 3 and _isnum(t[i + 3])):
                    raise ValueError(f"segment {key} has no energy-shift field to write.")
                edits[i + 3] = repr(float(shift))
            if not edits:
                continue
            for k in sorted(edits, reverse=True):
                a, b = spans[k]
                line = line[:a] + edits[k] + line[b:]
            lines[n] = line
            changed += 1
        missing = set(values) - seen
        if missing:
            raise KeyError(f"no <segmentsData> line with key(s) {sorted(missing)}.")
        self._set_suffix_lines(lines)
        return changed

    # -- explicit prior centres (rows of the <parameterSettings> block) -------

    def prior_centres(self):
        """The explicit prior centres of the file, ``{key: (norm, shift)}``.

        A row ``segment_<key>_norm prior_centre <c>`` (or
        ``segment_<key>_energy_shift prior_centre <c>``) in
        ``<parameterSettings>`` holds the centre of that segment's
        normalization (energy-shift) prior apart from the norm (shift) field
        of its ``<segmentsData>`` line, which is then only the start value.
        Without a row the field is both, as in every classic file.  Only the
        keys with a row are listed; the other entry of a pair is None.
        """
        out = {}
        for line in self._parameter_settings_lines():
            m = _PRIOR_ROW.match(line)
            if not m:
                continue
            key = int(m.group(1))
            n, sh = out.get(key, (None, None))
            v = float(m.group(3))
            out[key] = (v, sh) if m.group(2) == "norm" else (n, v)
        return out

    def _parameter_settings_lines(self):
        lines, inside = [], False
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<parameterSettings>":
                inside = True
            elif s.startswith("<"):
                inside = False
            elif inside:
                lines.append(line)
        return lines

    def set_prior_centres(self, values):
        """Pin the prior centres of normalizations and energy shifts.

        ``values`` maps a segment key (see :meth:`segment_values`) to
        ``(norm, shift)``: a number writes that centre, None leaves it as it
        is, and :data:`CLASSIC` removes the row (the field is the centre
        again).  The rows are kept together, in key order, at the end of
        ``<parameterSettings>`` -- the order the GUI writes them in; the block
        is created after ``</targetInt>`` if the file has none.  Older
        AZURE2 versions skip these three-token rows and read the field as the
        centre.  Returns the model.
        """
        have = self.prior_centres()
        keys = set(self.segment_values())
        for key, (n, sh) in values.items():
            if key not in keys:
                raise KeyError(f"no <segmentsData> line with key {key}.")
            old = list(have.get(key, (None, None)))
            for i, v in enumerate((n, sh)):
                if v is CLASSIC:
                    old[i] = None
                elif v is not None:
                    v = float(v)
                    if not math.isfinite(v) or (i == 0 and v <= 0.0):
                        raise ValueError(f"segment {key}: prior centre {v!r} "
                                         "(a norm centre must be positive).")
                    old[i] = v
            have[key] = tuple(old)
        rows = []
        for key in sorted(have):
            for name, v in zip(("norm", "energy_shift"), have[key]):
                if v is not None:
                    rows.append(f"segment_{key}_{name} prior_centre {_fmt(v)}")
        lines = self._suffix.splitlines()
        lines = [l for l in lines
                 if not _PRIOR_ROW.match(l) and l.strip() != _PRIOR_COMMENT]
        block = [_PRIOR_COMMENT] + rows if rows else []
        try:
            start = next(i for i, l in enumerate(lines)
                         if l.strip() == "<parameterSettings>")
            end = next(i for i in range(start + 1, len(lines))
                       if lines[i].strip().startswith("<"))
            if lines[end].strip() != "</parameterSettings>":
                raise ValueError("<parameterSettings> is not closed.")
            lines[end:end] = block
        except StopIteration:
            if rows:
                at = next((i + 1 for i, l in enumerate(lines)
                           if l.strip() == "</targetInt>"), len(lines))
                lines[at:at] = (["<parameterSettings>"] + block
                                + ["</parameterSettings>"])
        self._set_suffix_lines(lines)
        return self

    def engine_level_keys(self):
        """``{(jgroup, level): AzrLevel}`` -- the numbering AZURE2 itself uses.

        The engine builds a J-group the first time it meets a ``(J, parity)``
        while reading ``<levels>``, and numbers levels within a group in file
        order.  Reproducing that here gives the same ``(jgroup, level)`` a
        :class:`~pyazr.parameters.Parameter` reports, which is the only reliable
        way to say *which* level a fitted value belongs to: a level at
        ``Ex = 0`` comes back from the API with ``level_energy = None``, so an
        energy-based key cannot identify it.

        Both indices are 1-based, as the API reports them.

        **Inactive levels are skipped**, because the engine skips them
        (``CNuc::Fill`` reads only lines with ``isActive == 1``): an inactive
        line neither opens a J-group nor takes a level number.  Counting them
        here shifted every later level of the group by one, so a fit applied
        through :meth:`apply_fit` landed one level early -- the 13C+alpha
        archive's 5/2+ block was silently corrupted that way at two bakes.
        Use :meth:`purge_inactive_levels` to remove such lines for good.
        """
        order, seen, out = [], {}, {}
        for lv in self.levels:
            if not lv.active:
                continue
            k = (int(round(2 * lv.J)), int(lv.parity))
            if k not in seen:
                order.append(k)
                seen[k] = 0
            seen[k] += 1
            out[(order.index(k) + 1, seen[k])] = lv
        return out

    @property
    def active_levels(self):
        """The levels AZURE2 will actually read (``isActive == 1``), in file order."""
        return [lv for lv in self.levels if lv.active]

    def purge_inactive_levels(self):
        """Delete every inactive level from the file.  Returns the removed
        :class:`AzrLevel` objects.

        An inactive level contributes nothing to the engine, but it is a trap
        for anything that numbers levels from the file text (older versions of
        :meth:`engine_level_keys`, hand edits, the GUI's level table), so a
        model that is going to be edited or baked programmatically is safer
        without them.  Level IDs are renumbered on the next write.
        """
        gone = [lv for lv in self.levels if not lv.active]
        self.levels = [lv for lv in self.levels if lv.active]
        self._renumber()
        return gone

    def apply_fit(self, parameters, x, transform=None, physical=False,
                  pairs=None, strict=True, include_fixed=False):
        """Write a fitted parameter vector into the levels block.

        **The ``gamma`` field of a ``<levels>`` line is not a reduced-width
        amplitude.**  It holds the same *physical* value AZURE2 prints in
        ``parameters.out`` and shows in the GUI: a partial width in eV for an
        open particle channel, an ANC in fm^-1/2 for a closed (sub-threshold)
        one, and a partial width in eV for a photon channel.  Writing an rwa
        there produces a file that loads without complaint and is wrong -- for
        the 7Be model the two differ by factors of 10^2 to 10^7.

        So ``x`` (a free vector in ``params_rwa`` order) must be converted
        first.  Either hand over the transform and let this do it:

        >>> mdl.apply_fit(m.parameters, x_best, transform=m.transform_rwa)

        or convert yourself and say so:

        >>> mdl.apply_fit(m.parameters, m.transform_rwa(x_best), physical=True)

        Passing neither raises, rather than silently writing an rwa.

        Levels are matched on the ``(jgroup, level)`` the engine reports, via
        :meth:`engine_level_keys`, and channels on ``(pair, L, S)`` within the
        matched level.

        **Pass ``pairs=m.pairs``.**  A ``Parameter``'s ``pair`` is the engine's
        number, which counts particle pairs in the order ``<levels>`` first
        mentions them -- not the pair key the file writes.  The two differ
        whenever the levels do not introduce the pairs in key order: on the 8Be
        model engine pair 1 is file key 2, and file key 1 is engine pair 6, so
        matching one against the other writes every width to the wrong channel.
        The :class:`~pyazr.parameters.PairSet` carries both, and is the only
        thing that can translate.

        ``strict`` (the default) raises if any free parameter finds no home,
        rather than skipping it: a partial write produces a file that loads
        cleanly and is a mixture of two fits.

        ``include_fixed``: ``x`` is a vector over *every* parameter (index
        :attr:`Parameter.index`, as ``transform_all_rwa(..., include_fixed=True)``
        returns it) and the fixed energies and widths are written as well.  A
        fixed width holds its reduced-width amplitude in a fit, so its physical
        value moves with the level's other widths; writing it keeps the file
        the model the fit had.

        Note this covers ``<levels>`` only -- normalizations and energy shifts
        are not in that block, so a fit that moved them is only half saved.
        :meth:`pyazr.azure2.azure2.save_fit` writes the companion
        ``param.sav`` and verifies the result; prefer it to calling this
        directly.  The number of values written is left in :attr:`applied`.
        """
        if transform is None and not physical:
            raise ValueError(
                "apply_fit needs physical values, not reduced-width amplitudes: "
                "pass transform=m.transform_rwa, or convert with "
                "m.transform_rwa(x) yourself and pass physical=True.")
        if transform is not None:
            if physical:
                raise ValueError("pass transform= or physical=True, not both.")
            x = transform(x)

        lvlmap = self.engine_level_keys()
        # engine pair number -> the pair key the file writes
        pairkey = {p.number: p.key for p in pairs} if pairs is not None else {}
        written, unplaced = 0, []
        for p in parameters:
            slot = p.index if include_fixed else p.free_index
            if not include_fixed and (p.fixed or p.free_index is None):
                continue
            if p.kind not in ("energy", "width"):
                continue            # norms and shifts do not live in <levels>
            if slot is None or slot >= len(x):
                unplaced.append(f"{p.name} (index {slot} beyond the vector)")
                continue
            lv = lvlmap.get((p.jgroup, p.level))
            if lv is None:
                unplaced.append(f"{p.name} (no level at jgroup {p.jgroup}, level {p.level})")
                continue
            v = float(x[slot])
            if p.kind == "energy":
                lv.set_energy(v)
                written += 1
                continue
            want_pair = pairkey.get(p.pair, p.pair)
            for c in lv.channels:
                if (c.pair == want_pair and c.L == p.L
                        and abs(c.S - (p.S or 0.0)) < 1e-6):
                    c.gamma = v
                    written += 1
                    break
            else:
                unplaced.append(
                    f"{p.name} (no channel pair={want_pair} L={p.L} S={p.S} "
                    f"in {lv.jpi} at {lv.energy} MeV)")

        if unplaced and strict:
            raise ValueError(
                f"apply_fit could not place {len(unplaced)} of the free "
                f"parameters, so the result would be a mixture of two fits:\n  "
                + "\n  ".join(unplaced[:8])
                + (f"\n  ... and {len(unplaced) - 8} more" if len(unplaced) > 8 else "")
                + "\nThe .azr and the parameter set do not describe the same model. "
                  "Pass strict=False to write the rest anyway."
                + ("\nNote apply_fit was given no pairs=, so it assumed the "
                   "engine's pair numbers are the file's pair keys. Pass "
                   "pairs=m.pairs if they are not." if pairs is None else ""))
        self.applied = written
        self.unplaced = unplaced
        return self

    def set_segment_datafile(self, file_substr, new_path):
        """Repoint matching ``<segmentsData>`` lines to ``new_path`` (the data
        file is the first non-numeric token).  Returns the number changed."""
        out, changed, inside = [], 0, False
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
            elif s == "</segmentsData>":
                inside = False
            elif inside and s and file_substr in line:
                t = line.split()
                for i, tok in enumerate(t):
                    if not _isnum(tok):
                        t[i] = new_path
                        break
                line = " ".join(t)
                changed += 1
            out.append(line)
        self._set_suffix_lines(out)
        if changed == 0:
            raise KeyError(f"no <segmentsData> line matches {file_substr!r}.")
        return changed

    def set_segment_active(self, file_substr, active):
        """Activate/deactivate every ``<segmentsData>`` line matching
        ``file_substr`` (sets the leading isActive field).  Deactivated segments
        are ignored by AZURE2.  Returns the number of lines changed."""
        if "<segmentsData>" not in self._suffix:
            raise ValueError("no <segmentsData> block to edit.")
        out, changed, inside = [], 0, False
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
            elif s == "</segmentsData>":
                inside = False
            elif inside and s and file_substr in line:
                t = line.split()
                t[0] = "1" if active else "0"
                line = " ".join(t)
                changed += 1
            out.append(line)
        self._set_suffix_lines(out)
        if changed == 0:
            raise KeyError(f"no <segmentsData> line matches {file_substr!r}.")
        return changed

    def set_segment_energy_range(self, file_substr, energy_min=None,
                                 energy_max=None):
        """Cap or widen the lab-energy window of every ``<segmentsData>`` line
        matching ``file_substr`` (the ``minE`` / ``maxE`` fields, lab MeV).
        A ``None`` leaves that bound as it is.  Returns the number of lines
        changed.

        AZURE2 only loads the data points inside the window, so this is how an
        energy-stepped fit grows its data set -- and, like any change to the
        loaded grid, it invalidates ``output/intEC.dat`` (delete it, or give
        the edited model its own output directory).
        """
        if "<segmentsData>" not in self._suffix:
            raise ValueError("no <segmentsData> block to edit.")
        out, changed, inside = [], 0, False
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
            elif s == "</segmentsData>":
                inside = False
            elif inside and s and file_substr in line:
                t = line.split()
                if energy_min is not None:
                    t[3] = _fmt(float(energy_min))
                if energy_max is not None:
                    t[4] = _fmt(float(energy_max))
                line = " ".join(t)
                changed += 1
            out.append(line)
        self._set_suffix_lines(out)
        if changed == 0:
            raise KeyError(f"no <segmentsData> line matches {file_substr!r}.")
        return changed

    # -- adding / removing whole data segments --------------------------------

    # observable name -> isDiff code for a <segmentsData> line
    # (ESegment::ESegment(SegLine); mirrors datasets._OBSERVABLE).
    _DATA_CODE = {
        "angle-integrated": 0, "differential": 1, "phase-shift": 2,
        "total-capture": 3, "differential-cm": 4, "angle-integrated-E1": 5,
        "angle-integrated-E2": 6, "analyzing-power": 7,
    }

    def add_data_segment(self, data_file, entrance, exit,
                         observable="angle-integrated",
                         energy_min=0.0, energy_max=5.0,
                         angle_min=0.0, angle_max=180.0,
                         norm=1.0, vary_norm=False, norm_error=0.0, thm=False,
                         energy_shift=0.0, energy_shift_error=0.0,
                         vary_shift=False, phase_J=None, phase_L=None,
                         active=True):
        """Append one data segment (a ``<segmentsData>`` line).

        ``data_file`` is the data file path (relative to the run directory,
        which for a model living in its own folder is usually ``data/...``).
        ``entrance`` / ``exit`` are particle-pair keys (``exit=-1`` for a
        summed/total observable).  ``observable`` is one of
        ``angle-integrated``, ``differential``, ``differential-cm``,
        ``total-capture``, ``phase-shift``, ``angle-integrated-E1``,
        ``angle-integrated-E2``.

        ``norm`` is the normalization applied to the data, ``norm_error`` its
        systematic error (percent, as stored in the file), ``energy_shift``
        the beam-energy shift (MeV) with its ``_error``.  ``vary_norm`` /
        ``vary_shift`` free the corresponding parameter.

        Returns ``self`` so calls chain.
        """
        if observable not in self._DATA_CODE:
            raise ValueError(f"unknown observable {observable!r}; expected one "
                             f"of {sorted(self._DATA_CODE)}.")
        isDiff = self._DATA_CODE[observable]
        if thm:
            # A Trojan-Horse (half-off-shell) segment is its ordinary
            # observable with an isDiff offset of +10 -- the flag the engine
            # reads in ESegment::ESegment(SegLine).  The data are the
            # THM-extracted two-body cross section, whose scale is arbitrary:
            # give the segment a free normalization (vary_norm=True) and no
            # penalty (norm_error=0), so the fit sets the scale and nothing
            # pulls it back.
            isDiff += 10
        toks = [1 if active else 0, int(entrance), int(exit),
                _fmt(energy_min), _fmt(energy_max),
                _fmt(angle_min), _fmt(angle_max), isDiff]
        if isDiff == 2:                       # phase shift carries J, L
            if phase_J is None or phase_L is None:
                raise ValueError("a phase-shift data segment needs phase_J "
                                 "and phase_L.")
            toks += [_fmt(phase_J), int(phase_L)]
        toks += [_fmt(norm), 1 if vary_norm else 0, _fmt(norm_error),
                 _fmt(energy_shift), _fmt(energy_shift_error),
                 1 if vary_shift else 0, str(data_file), 0, 0]
        line = "  ".join(t if isinstance(t, str) else _fmt(t) for t in toks)
        lines = self._suffix.splitlines()
        if "<segmentsData>" not in lines:
            raise ValueError("no <segmentsData> block to add to.")
        end = lines.index("</segmentsData>")
        self._set_suffix_lines(lines[:end] + [line] + lines[end:])
        return self

    def remove_data_segments(self, file_substr):
        """Remove every ``<segmentsData>`` line whose text matches
        ``file_substr`` (e.g. a data-file name).  Returns the number of
        segments removed.  Raises if nothing matches.

        Removing data changes which energies AZURE2 evaluates, so the
        external-capture integrals must be recalculated before the next run --
        delete ``output/intEC.dat`` / ``output/intEC.extrap`` (or write the
        model into its own output directory) so the stale cache cannot be
        reused.
        """
        if "<segmentsData>" not in self._suffix:
            raise ValueError("no <segmentsData> block to edit.")
        out, removed, inside = [], 0, False
        for line in self._suffix.splitlines():
            s = line.strip()
            if s == "<segmentsData>":
                inside = True
            elif s == "</segmentsData>":
                inside = False
            elif inside and s and file_substr in line:
                removed += 1
                continue
            out.append(line)
        self._set_suffix_lines(out)
        if removed == 0:
            raise KeyError(f"no <segmentsData> line matches {file_substr!r}.")
        return removed

    def clear_data_segments(self):
        """Remove every ``<segmentsData>`` line (leave the block empty).

        The external-capture caches ``output/intEC.dat`` / ``output/intEC.extrap``
        belong to the removed grids and must be deleted before the next run.
        """
        lines = self._suffix.splitlines()
        try:
            start = lines.index("<segmentsData>")
            end = lines.index("</segmentsData>")
        except ValueError:
            raise ValueError("no <segmentsData> block to clear.")
        self._set_suffix_lines(lines[:start + 1] + lines[end:])
        return self

    # -- experimental effects (edits the <targetInt> block) -------------------

    def target_effects(self):
        """The raw ``<targetInt>`` lines (one experimental effect each)."""
        lines = self._suffix.splitlines()
        try:
            start = lines.index("<targetInt>") + 1
            end = lines.index("</targetInt>")
        except ValueError:
            return []
        return [ln for ln in lines[start:end] if ln.strip()]

    def _splice_target_int(self, new_lines):
        lines = self._suffix.splitlines()
        try:
            start = lines.index("<targetInt>")
            end = lines.index("</targetInt>")
        except ValueError:
            self._suffix = (self._suffix.rstrip("\n") + "\n\n<targetInt>\n"
                            + "\n".join(new_lines) + "\n</targetInt>"
                            + self._final_newline())
            return
        self._set_suffix_lines(lines[:start + 1] + new_lines + lines[end:])

    def clear_target_effects(self):
        """Remove every ``<targetInt>`` line (leave the block empty)."""
        self._splice_target_int([])
        return self

    def add_target_effect(self, segments, n_points=200, gaussian_sigma=None,
                          beam_profile=None, tpc_sigma=0.0, truncation=0.0,
                          photodissociation=False, q_coefficients=None,
                          resonance_width_multiplier=20.0, points_per_width=50.0,
                          active=True):
        """Append one experimental effect (a ``<targetInt>`` line).

        ``segments`` is the segment-key list as AZURE2 writes it (``"3"``,
        ``"3-5"``, ``"3,7-9"`` or an iterable of ints).  Keys count *every*
        ``<segmentsData>`` line, active or not, and a ``<segmentsTest>`` line
        with the same key gets the effect too (see the azure2-eval skill).

        ``gaussian_sigma`` (lab MeV) is the classic beam-energy Gaussian
        convolution.  ``beam_profile`` is the beam-profile kernel: a list of
        ``(xi, omega, alpha, weight)`` skewed-Gaussian components in lab
        entrance-channel energy (MeV), with the detector energy resolution
        ``tpc_sigma`` (lab MeV; the per-point energy window comes from
        columns 5-6 of the data file), an optional ``truncation`` of each
        component at mean +- n standard deviations (0 = none) and, for the
        inverse reaction of a photodissociation measurement,
        ``photodissociation=True`` to weight the average with the
        detailed-balance factor.  ``q_coefficients`` are the finite-geometry
        attenuation coefficients Q_0..Q_n.  Target integration and straggling
        are not exposed here.  Returns ``self``.
        """
        if not isinstance(segments, str):
            segments = ",".join(str(int(k)) for k in segments)
        toks = [1 if active else 0, f'"{segments}"', int(n_points)]
        if gaussian_sigma is not None:
            toks += [1, _fmt(gaussian_sigma)]
        else:
            toks += [0, 0]
        toks += [0, 0, '""', 0]                        # no target integration
        if q_coefficients:
            toks += [1, len(q_coefficients)] + [_fmt(q) for q in q_coefficients]
        else:
            toks += [0, 0]
        toks += [0, '""', 0]                           # no energy-dependent sigma
        toks += [0, 0.04, _fmt(resonance_width_multiplier), _fmt(points_per_width)]
        if beam_profile:
            toks += ["beamprofile", len(beam_profile)]
            for xi, omega, alpha, weight in beam_profile:
                toks += [_fmt(xi), _fmt(omega), _fmt(alpha), _fmt(weight)]
            toks += [_fmt(tpc_sigma), _fmt(truncation), 1 if photodissociation else 0]
        line = "  ".join(t if isinstance(t, str) else _fmt(t) for t in toks)
        lines = self._suffix.splitlines()
        if "<targetInt>" in lines:
            end = lines.index("</targetInt>")
            self._set_suffix_lines(lines[:end] + [line] + lines[end:])
        else:
            self._splice_target_int([line])
        return self

    # -- THM options (the <thm> block) ----------------------------------------
    #
    # The engine reads the block from anywhere in the file (Config::ReadThmBlock)
    # and refuses a line it does not understand; the GUI's THM Options dialog
    # (gui/src/ThmOptionsDialog.cpp) edits it with the same rules.  These
    # methods are the GUI's, line for line: the same keys and values, the same
    # normalization of a value, and the same way of writing the block back --
    # comment and blank lines stay where they were, a line whose value did not
    # change is kept verbatim, a changed one is rewritten keeping its
    # indentation and inline comment, a key that returns to its default is
    # removed, a new key is appended, and options that are all default remove
    # the block (comments included), since an empty block and no block are the
    # same to the engine.  A new block goes after </targetInt>, where the GUI
    # writes it.

    def _thm_locate(self):
        """(attribute, lines, open, close) of the <thm> block, or None.

        ``lines`` is the text of ``self._prefix`` or ``self._suffix`` split with
        line ends kept; ``open``/``close`` index its ``<thm>`` and ``</thm>``
        lines.  Raises ValueError for an unterminated block (the engine refuses
        the file)."""
        for attr in ("_prefix", "_suffix"):
            lines = getattr(self, attr).splitlines(keepends=True)
            for i, line in enumerate(lines):
                if line.strip().startswith("<thm>"):
                    for j in range(i + 1, len(lines)):
                        if lines[j].split("#", 1)[0].strip() == "</thm>":
                            return attr, lines, i, j
                    raise ValueError("the <thm> block is not terminated by "
                                     "</thm> (AZURE2 refuses the file).")
        return None

    def _thm_body(self):
        loc = self._thm_locate()
        if loc is None:
            return []
        _, lines, i, j = loc
        return [ln.rstrip("\r\n") for ln in lines[i + 1:j]]

    def _thm_settings(self):
        """The block parsed into the GUI's ThmSettings (a dict); ValueError on
        a line the engine would refuse."""
        s = _thm_default_settings()
        for line in self._thm_body():
            _thm_parse_line(line, s)
        _thm_check_experiments(s["experiments"])
        return s

    def thm_options(self, defaults=False):
        """The options of the ``<thm>`` block, as a dict.

        Keys are the block's keys (``vertex``, ``kinematics``, ``entranceL``,
        ``coulombIntegral``, ``spectatorEnergy``, ``spectatorEnergy[<pair>]``,
        ``weight[<k>]``, ``weightTest[<k>]``); values are str, bool
        (``coulombIntegral``), float (spectator energies, MeV) and the file as
        written (weights).  Only the options that differ from the engine's
        defaults are listed -- with ``defaults=True`` the five global options
        are always there.  Raises ValueError if the block has a line AZURE2
        would refuse.  See docs/source/theory/thm_implementation.rst.
        """
        s = self._thm_settings()
        out = {}
        d = _thm_default_settings()
        for key in ("entranceL", "vertex", "kinematics", "coulombIntegral",
                    "spectatorEnergy"):
            if defaults or s[key] != d[key]:
                out[key] = s[key]
        for k in sorted(s["spectatorByPair"]):
            out[f"spectatorEnergy[{k}]"] = s["spectatorByPair"][k]
        for k in sorted(s["weight"]):
            out[f"weight[{k}]"] = s["weight"][k]
        for k in sorted(s["weightTest"]):
            out[f"weightTest[{k}]"] = s["weightTest"][k]
        return out

    def set_thm_option(self, key, value):
        """Set one option of the ``<thm>`` block, validated as AZURE2 does.

        ``key`` is a block key (see :meth:`thm_options`); ``value`` a str, a
        bool for ``coulombIntegral``, a number >= 0 for a spectator energy (MeV).
        ``weight[<k>]`` / ``weightTest[<k>]`` take a file name; prefer
        :meth:`set_thm_weight`, which also checks the segment.  Setting an
        option to its default removes its line.  Raises ValueError for an
        unknown key or a value the engine would refuse, and leaves the model
        unchanged then.
        """
        if isinstance(value, bool):
            text = "1" if value else "0"
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            text = _thm_number(float(value))
        else:
            text = str(value)
        if "\n" in text or "\r" in text:
            raise ValueError(f"<thm> {key}: a value must be one line.")
        s = self._thm_settings()
        try:
            canon_key, canon_value = _thm_parse_line(f"{key}={text}", s)
        except ValueError:
            raise ValueError(f"<thm> {key}={text}: not an option AZURE2 "
                             "accepts (see thm_options()).") from None
        if not canon_key:
            raise ValueError(f"<thm> {key!r}: not an option.")
        if canon_key.startswith("weight"):
            self._thm_check_weight_file(canon_key, canon_value)
        if s["entranceL"] == "coherent":
            for name, x in s["experiments"].items():
                if x.get("theta", "all") != "all":
                    raise ValueError(f"<thm> experiment[{name}]: theta= computes the "
                                     "interference of the entrance partial waves exactly; "
                                     "entranceL=coherent cannot be combined with it.")
                if "cbackground" in x:
                    raise ValueError(f"<thm> experiment[{name}]: cbackground= adds an "
                                     "amplitude per entrance bucket (s, l); "
                                     "entranceL=coherent cannot be combined with it.")
        self._thm_check_dw_globals(s)
        self._thm_write(s)
        return self

    def _thm_check_dw_globals(self, s, only=None):
        """The global options an experiment refuses: coulombIntegral=1 with
        vertexModel=dw or with a computed R(E) whose a + A wave is distorted
        (CheckThmCoulombConsistency); entranceL=coherent and a spectator
        energy for its entrance pair with vertexModel=dw (EData::
        BuildThmGroups)."""
        seg_lines = self._block_lines("segmentsData") or []
        for name, x in s["experiments"].items():
            if only is not None and name != only:
                continue
            where = f"<thm> experiment[{name}]: "
            dist = x.get("distortion", "none")
            if (s["coulombIntegral"] and x.get("vertexModel", "pw") != "dw"
                    and (dist == "coulomb"
                         or (dist == "optical" and x.get("opticalAA", "coulomb") != "plane"))):
                raise ValueError(where + "coulombIntegral=1 adds the x-A Coulomb "
                                 "interaction outside the channel radius (C_l) to the "
                                 "vertex, and the a + A Coulomb wave of the distortion "
                                 "factor R(E) already contains it (Z_a = Z_x + Z_s; in the "
                                 "zero range of R the prior operator V_xA + V_sA - U_aA "
                                 "outside the radius vanishes), so it would be counted "
                                 f"twice. Use coulombIntegral=0 with distortion={dist}, or "
                                 "opticalAA=plane.")
            if x.get("vertexModel", "pw") != "dw":
                continue
            if s["coulombIntegral"]:
                raise ValueError(where + "vertexModel=dw is the surface term of the DWBA "
                                 "vertex; its external part (the three-body remnant outside "
                                 "the channel radius, whose plane-wave limit is the Coulomb "
                                 "term C_l) is not computed, so it cannot be combined with "
                                 "coulombIntegral=1.")
            if s["entranceL"] == "coherent":
                raise ValueError(where + "vertexModel=dw sums the entrance partial waves "
                                 "(and their projections) incoherently, as the "
                                 "angle-integrated observable requires; entranceL=coherent "
                                 "cannot be combined with it.")
            for k in x.get("segments", []):
                if k > len(seg_lines):
                    continue
                tok = seg_lines[k - 1].split()
                key = int(float(tok[1])) if len(tok) > 1 and _isnum(tok[1]) else None
                if s["spectatorEnergy"] != 0.0 or s["spectatorByPair"].get(key, 0.0) != 0.0:
                    raise ValueError(where + "with vertexModel=dw the spectator kinematics "
                                     "come from Ebeam and the spectator direction; "
                                     f"spectatorEnergy for entrance pair {key} must be 0.")

    def set_thm_weight(self, segment, path, test=False):
        """``weight[<segment>]=<path>`` (``weightTest`` with ``test=True``): the
        energy-dependent weight w(E) of the THM model of one segment.

        ``segment`` counts every line of ``<segmentsData>`` (``<segmentsTest>``),
        inactive ones included -- the engine's segment key.  The line must
        exist and be a THM segment (isDiff >= 10).  ``path`` is a two-column
        table (E_cm of the THM entrance pair in MeV, w > 0; ``#`` comments;
        at least two rows, E strictly increasing), read here with the engine's
        rules.  A relative path is taken from the directory of the .azr --
        write the edited model next to its source, or give an absolute path.
        """
        segment = int(segment)
        block = "segmentsTest" if test else "segmentsData"
        seg_lines = self._block_lines(block)
        if seg_lines is None:
            raise ValueError(f"no <{block}> block.")
        if not 1 <= segment <= len(seg_lines):
            raise ValueError(f"<{block}> has {len(seg_lines)} lines; no "
                             f"segment {segment}.")
        tok = seg_lines[segment - 1].split()
        is_diff = int(float(tok[7])) if len(tok) > 7 and _isnum(tok[7]) else -1
        if is_diff < 10:
            raise ValueError(f"<{block}> line {segment} is not a THM segment "
                             f"(isDiff {is_diff} < 10); AZURE2 refuses a "
                             "weight on it.")
        key = f"{'weightTest' if test else 'weight'}[{segment}]"
        return self.set_thm_option(key, str(path))

    def clear_thm_option(self, key):
        """Remove an option from the ``<thm>`` block (back to its default).
        Removing the last one removes the block.  Unknown keys raise
        ValueError; a key that is not set is a no-op."""
        s = self._thm_settings()
        key = str(key).strip()
        if key in _THM_GLOBAL_DEFAULTS:
            s[key] = _thm_default_settings()[key]
        else:
            m = re.fullmatch(r"(spectatorEnergy|weightTest|weight)\[(.+)\]", key)
            if not m:
                raise ValueError(f"<thm> {key!r}: not an option.")
            try:
                k = int(m.group(2))
            except ValueError:
                raise ValueError(f"<thm> {key!r}: not an option.") from None
            table = {"spectatorEnergy": "spectatorByPair", "weight": "weight",
                     "weightTest": "weightTest"}[m.group(1)]
            s[table].pop(k, None)
        self._thm_write(s)
        return self

    def clear_thm_weight(self, segment, test=False):
        """Remove ``weight[<segment>]`` (``weightTest`` with ``test=True``)."""
        return self.clear_thm_option(
            f"{'weightTest' if test else 'weight'}[{int(segment)}]")

    # -- THM experiments: experiment[<name>] lines of the <thm> block ---------

    def thm_experiments(self):
        """The THM experiments of the ``<thm>`` block, ``{name: record}``.

        A record has ``segments`` (list of ``<segmentsData>`` line numbers,
        inactive lines counted, like ``weight[k]``), ``background`` (``none``,
        ``const``, ``linear`` or ``quadratic``) and, if given, ``beam``,
        ``target``, ``spectator`` (nuclide names or ``Z,A,mass``) and
        ``Ebeam`` (lab MeV), ``lineshape: True`` when the Coulomb line
        shape of the spectator is on, ``ps`` / ``psNodes`` when given
        (spectator-momentum window), the distortion keys, and ``theta`` when
        given (the angular window of the fixed-angle observable, "all" or
        "thmin-thmax" as written).  All segments of an experiment share one
        profiled norm and the background; see
        docs/source/theory/thm_implementation.rst, "THM experiments" and
        "Coulomb line shape".  Raises ValueError if the block has a line AZURE2
        would refuse.
        """
        s = self._thm_settings()
        return {name: _thm_experiment_record(x) for name, x in s["experiments"].items()}

    def set_thm_experiment(self, name, segments, background="none", beam=None,
                           target=None, spectator=None, Ebeam=None, lineshape=False,
                           ps=None, psNodes=None, distortion=None, opticalAA=None,
                           opticalSF=None, spectatorAngle=None, distortionRef=None,
                           distortionRatio=None, boundState=None, theta=None,
                           vertexModel=None, spectatorAngles=None,
                           spectatorAngleNodes=None, cbackground=None):
        """Define (or replace) ``experiment[<name>]`` in the ``<thm>`` block.

        ``segments`` is a list of ``<segmentsData>`` line numbers (or the
        engine's text form, ``"1,2,4-6"``); each must exist, be a THM segment
        (isDiff >= 10) with a free norm, and belong to no other experiment.
        ``background`` is ``none``, ``const``, ``linear`` or ``quadratic``.
        ``beam``, ``target``, ``spectator`` (a name of the built-in table --
        n p d t 3He 4He 6Li 7Li 9Be 10B 11B 12C 13C 14N 15N 16O 17O 18O 19F
        20Ne 23Na 24Mg -- or ``"Z,A,mass"``, nuclear mass in u) and ``Ebeam``
        (lab MeV) go together: all four or none.  ``lineshape=True`` turns on
        the Coulomb line-shape factor N_C of the spectator (needs the four
        kinematics keys; see :meth:`pyazr.azure2.azure2.thm_lineshape`).
        ``ps`` sets the spectator-momentum window over which the HOES cross
        section is averaged (weight |phi(p_s)|^2 p_s^2, momenta in MeV/c):
        ``"delta"`` (the quasi-free point, the default), ``"hulthen:0-40"``,
        ``"hulthen:a,b:0-40"`` (a, b in fm^-1; default the deuteron's 0.2317,
        1.202), ``"gauss:FWHM:0-40"`` or ``"table:<file>"`` (p_s and the event
        weight per unit p_s; relative to the .azr, read by the engine); a
        window needs the kinematics keys and excludes ``spectatorEnergy`` for
        the experiment's entrance pair.  ``psNodes`` (1-64, default 16) is the
        number of Gauss-Legendre nodes in p_s (see
        :meth:`pyazr.azure2.azure2.thm_vertex`).
        ``distortion`` multiplies the model of every segment by the
        distortion factor R(E) (see :meth:`pyazr.azure2.azure2.thm_distortion`):
        ``"coulomb"`` (point-Coulomb waves in a + A and s + F), ``"optical"``
        (per channel ``opticalAA`` / ``opticalSF``: ``"plane"``,
        ``"coulomb"``, a global optical potential -- ``"ancai06"``,
        ``"daehnick80"`` (d), ``"kd03"`` (n, p), ``"bg71"`` (t, 3He),
        ``"liang09"`` (3He), ``"mcfadden66"``, ``"avrigeanu94"`` (4He), with
        ``":extrapolate"`` to allow it outside its validity range -- or ten
        numbers ``"V,R,a,W,RW,aW,WD,RD,aD,RC"`` in MeV and
        fm, or a sequence of them) -- both need the kinematics keys -- or
        ``"table:<file>"`` (E, w columns as ``weight[k]``).  With coulomb or
        optical: ``spectatorAngle`` (``"qf"``, a lab angle in degrees or
        ``"cm:<deg>"``), ``distortionRef`` (E_ref, MeV), ``distortionRatio``
        (``"dwpw"`` or ``"dw"``), ``boundState`` (``"whittaker"`` or
        ``"yukawa"``, optionally ``":rmin"`` in fm).
        ``theta`` turns the model of every segment into the fixed-angle HOES
        observable: dsigma/dOmega averaged over the c.m. angle of the exit
        pair (particle 1 relative to 2, from p_xA = entrance particle 1
        relative to 2) in a window, ``"50-70"`` or ``(50, 70)`` in degrees
        (0 <= min <= max <= 180; min == max is one angle); ``"all"`` (or None)
        keeps the angle-integrated cross section.  0-180 gives it / 4 pi.  Not
        with ``entranceL=coherent``.
        ``vertexModel`` is ``"pw"`` (the plane-wave vertex M_l, the default) or
        ``"dw"``: the entrance vertex of every segment becomes the surface term
        of the prior-form DWBA built from the experiment's distorted waves
        (needs ``distortion="coulomb"`` or ``"optical"``; R(E) is then not
        applied, so ``distortionRef``/``distortionRatio`` are refused, as are
        ``theta``, ``spectatorAngle`` together with a ``ps`` window,
        ``coulombIntegral=1``, ``entranceL=coherent`` and a spectator energy
        for the entrance pair; see :meth:`pyazr.azure2.azure2.thm_vertex`).
        ``spectatorAngles`` averages R(E) or the DW vertex over the accepted
        spectator directions instead of taking one (``spectatorAngle``, which
        it excludes): ``"10-30"`` or ``(10, 30)`` (lab polar angles of the
        spectator to the beam, degrees), ``"cm:150-180"`` (c.m.),
        ``"table:<file>"`` / ``"cm:table:<file>"`` (angle and acceptance
        columns, relative to the .azr); weight d cos(theta_cm) x acceptance; a
        ``ps`` window then only cuts |p_s| for the DW vertex and for R.
        ``spectatorAngleNodes`` (1-64, default 8) Gauss-Legendre nodes per
        c.m. interval (a lab window can map to two).  Needs ``distortion=
        "coulomb"`` or ``"optical"``; with ``vertexModel="dw"`` not together
        with ``psNodes``.
        ``cbackground`` adds a coherent (interfering) THM-only background
        amplitude c(E) M_l to the resonant HOES amplitude before squaring
        (docs/source/theory/thm_implementation.rst, "Coherent background"):
        a term ``"<J><+|->:<exit pair key>[:<s>,<l>,<s'>,<l'>][:const|:linear]
        [=<Re c0>,<Im c0>[,<Re c1>,<Im c1>]]"`` or a list of them (joined
        with ``;``); without the channels every (entrance (s,l), exit
        (s',l')) combination of the J^pi group gets its own amplitude; start
        values default to 0, a value followed by ``f`` is fixed.  The real
        and imaginary parts are ordinary fit parameters (``cbkg_*`` in the
        parameter vector, kind ``"cbkg"``).  Not with ``entranceL=coherent``;
        the J^pi must be a group of the model and the exit pair that of a
        segment of the experiment (the engine checks the channels when the
        session opens).
        The record replaces every
        earlier line of that name with one line; other lines stay as they
        are.  Raises ValueError (model unchanged) for anything AZURE2 would
        refuse.
        """
        if not isinstance(segments, str):
            try:
                segments = ",".join(str(int(k)) for k in segments)
            except TypeError:
                segments = str(int(segments))
        text = f"experiment[{name}] segments={segments}"
        if background != "none":
            text += f" background={background}"
        for key, value in (("beam", beam), ("target", target), ("spectator", spectator)):
            if value is not None:
                text += f" {key}={value}"
        if Ebeam is not None:
            text += f" Ebeam={_thm_number(float(Ebeam))}"
        if lineshape:
            text += " lineshape=on"
        if ps is not None:
            text += f" ps={ps}"
        if psNodes is not None:
            text += f" psNodes={psNodes}"
        for key, value in (("distortion", distortion), ("opticalAA", opticalAA),
                           ("opticalSF", opticalSF), ("spectatorAngle", spectatorAngle),
                           ("distortionRef", distortionRef),
                           ("distortionRatio", distortionRatio), ("boundState", boundState)):
            if value is None:
                continue
            if key in ("opticalAA", "opticalSF") and not isinstance(value, str):
                value = ",".join(_thm_plain_number(v) for v in value)
            elif key in ("distortionRef", "spectatorAngle") and not isinstance(value, str):
                value = _thm_plain_number(value)
            text += f" {key}={value}"
        if theta is not None:
            if not isinstance(theta, str):
                lo, hi = theta
                theta = f"{_thm_plain_number(lo)}-{_thm_plain_number(hi)}"
            text += f" theta={theta}"
        if vertexModel is not None:
            text += f" vertexModel={vertexModel}"
        if spectatorAngles is not None:
            if not isinstance(spectatorAngles, str):
                lo, hi = spectatorAngles
                spectatorAngles = f"{_thm_plain_number(lo)}-{_thm_plain_number(hi)}"
            text += f" spectatorAngles={spectatorAngles}"
        if spectatorAngleNodes is not None:
            text += f" spectatorAngleNodes={spectatorAngleNodes}"
        if cbackground is not None:
            if not isinstance(cbackground, str):
                cbackground = ";".join(str(t) for t in cbackground)
            text += f" cbackground={cbackground}"
        if any(c in text for c in "#\r\n"):
            raise ValueError(f"<thm> experiment[{name}]: a value cannot contain '#' "
                             "or a line break.")
        s = self._thm_settings()
        others = {k: v for k, v in s["experiments"].items() if k != name}
        trial = dict(others)
        _thm_parse_experiment(text, trial)
        _thm_check_experiments(trial)
        rec = _thm_experiment_record(trial[name])
        if rec.get("theta", "all") != "all" and s["entranceL"] == "coherent":
            raise ValueError(f"<thm> experiment[{name}]: theta= computes the interference "
                             "of the entrance partial waves exactly (at fixed angle they "
                             "interfere); entranceL=coherent is an approximation of the "
                             "angle-integrated observable and cannot be combined with it.")
        if "cbackground" in rec and s["entranceL"] == "coherent":
            raise ValueError(f"<thm> experiment[{name}]: cbackground= adds an amplitude "
                             "per entrance bucket (s, l); entranceL=coherent merges the l "
                             "of a channel spin, so it cannot be combined with it.")
        seg_lines = self._block_lines("segmentsData") or []
        for k in rec["segments"]:
            if k > len(seg_lines):
                raise ValueError(f"<thm> experiment[{name}]: segment {k}: "
                                 f"<segmentsData> has only {len(seg_lines)} line(s).")
            tok = seg_lines[k - 1].split()
            is_diff = int(float(tok[7])) if len(tok) > 7 and _isnum(tok[7]) else -1
            if is_diff < 10:
                raise ValueError(f"<thm> experiment[{name}]: segment {k} is not a "
                                 "THM segment (isDiff < 10).")
            vary = int(float(tok[9])) if len(tok) > 9 and _isnum(tok[9]) else 0
            if not vary:
                raise ValueError(f"<thm> experiment[{name}]: segment {k} has a fixed "
                                 "norm; the segments of an experiment share one free "
                                 "(profiled) norm, so free it.")
            if rec.get("ps", "delta") != "delta":
                key = int(float(tok[1])) if len(tok) > 1 and _isnum(tok[1]) else None
                if s["spectatorEnergy"] != 0.0 or s["spectatorByPair"].get(key, 0.0) != 0.0:
                    raise ValueError(f"<thm> experiment[{name}]: a ps window and "
                                     "spectatorEnergy both set the spectator motion of "
                                     f"entrance pair {key}; use one (ps=delta keeps "
                                     "spectatorEnergy).")
        if "cbackground" in rec:
            # What can be checked without the compound nucleus (EData::
            # BuildThmGroups checks the channels): the exit pair is that of a
            # segment of the experiment, the J^pi a group of <levels>.
            exits = set()
            for k in rec["segments"]:
                tok = seg_lines[k - 1].split()
                if len(tok) > 2 and _isnum(tok[2]):
                    exits.add(int(float(tok[2])))
            groups = {(round(2 * float(lv.J)), 1 if lv.parity > 0 else -1)
                      for lv in self.levels}
            for t in _thm_parse_cbackground(rec["cbackground"]):
                jpi = f"{_thm_spin_text(t['J'])}{'+' if t['parity'] > 0 else '-'}"
                if t["exit"] not in exits:
                    raise ValueError(f"<thm> experiment[{name}]: cbackground {jpi}:"
                                     f"{t['exit']}: no segment of the experiment has exit "
                                     f"pair {t['exit']}.")
                if (round(2 * t["J"]), t["parity"]) not in groups:
                    raise ValueError(f"<thm> experiment[{name}]: cbackground {jpi}:"
                                     f"{t['exit']}: the model has no J^pi = {jpi} group.")
        trial_settings = dict(s)
        trial_settings["experiments"] = trial
        self._thm_check_dw_globals(trial_settings, only=name)
        body = self._thm_body_without_experiment(name)
        body.append(_thm_experiment_line(name, rec))
        self._thm_set_body(body)
        return self

    def clear_thm_experiment(self, name):
        """Remove every ``experiment[<name>]`` line.  A block left with only
        default options (and no experiment) is removed, as the GUI does.
        Unknown names are a no-op."""
        s = self._thm_settings()
        if name not in s["experiments"]:
            return self
        body = self._thm_body_without_experiment(name)
        rest = _thm_default_settings()
        for line in body:
            _thm_parse_line(line, rest)
        if _thm_is_default(rest):
            body = []
        self._thm_set_body(body)
        return self

    def set_thm_cbackground(self, name, cbackground):
        """Replace the ``cbackground=`` value of ``experiment[<name>]`` in
        place (the line that carries it; other keys and lines untouched), e.g.
        with fitted start values.  ``cbackground`` as in
        :meth:`set_thm_experiment`.  KeyError if the experiment has no
        cbackground; ValueError (model unchanged) if AZURE2 would refuse it."""
        if not isinstance(cbackground, str):
            cbackground = ";".join(str(t) for t in cbackground)
        text = _thm_format_cbackground(_thm_parse_cbackground(cbackground))
        body = self._thm_body()
        head = f"experiment[{name}]"
        for i, line in enumerate(body):
            code, hash_, comment = line.partition("#")
            stripped = code.strip()
            if not (stripped.startswith(head) and stripped[len(head):][:1] in ("", " ", "\t")):
                continue
            tokens = code.split()
            for k, tok in enumerate(tokens):
                if tok.startswith("cbackground="):
                    new = re.sub(r"(?<!\S)cbackground=\S+", "cbackground=" + text, code, count=1)
                    trial = list(body)
                    trial[i] = new + hash_ + comment
                    s = _thm_default_settings()
                    for ln in trial:
                        _thm_parse_line(ln, s)
                    _thm_check_experiments(s["experiments"])
                    self._thm_set_body(trial)
                    return self
        raise KeyError(f"<thm> experiment[{name}] has no cbackground= to replace")

    def _thm_body_without_experiment(self, name):
        out = []
        for line in self._thm_body():
            code = line.split("#", 1)[0].strip()
            if code.startswith(f"experiment[{name}]") and \
                    code[len(f"experiment[{name}]"):][:1] in ("", " ", "\t"):
                continue
            out.append(line)
        return out

    def _thm_set_body(self, body):
        """Replace the block's lines by ``body``; no lines: no block."""
        loc = self._thm_locate()
        if not body:
            if loc is not None:
                attr, lines, i, j = loc
                setattr(self, attr, "".join(lines[:i] + lines[j + 1:]))
            return
        if loc is not None:
            attr, lines, i, j = loc
            nl = "\r\n" if lines[i].endswith("\r\n") else "\n"
            setattr(self, attr, "".join(lines[:i + 1] + [ln + nl for ln in body] + lines[j:]))
            return
        block = ["<thm>"] + body + ["</thm>"]
        lines = self._suffix.splitlines()
        if "</targetInt>" in lines:
            end = lines.index("</targetInt>") + 1
            self._set_suffix_lines(lines[:end] + block + lines[end:])
        else:
            self._set_suffix_lines(lines + block)

    def _block_lines(self, tag):
        """Non-blank lines of a block after <levels>, or None if absent."""
        lines = self._suffix.splitlines()
        try:
            start = lines.index(f"<{tag}>") + 1
            end = lines.index(f"</{tag}>")
        except ValueError:
            return None
        return [ln for ln in lines[start:end] if ln.strip()]

    def _thm_check_weight_file(self, key, name):
        name = name.strip()
        if "#" in name:
            raise ValueError(f"<thm> {key}: the path contains '#', which "
                             "starts a comment in the .azr.")
        base = os.path.dirname(os.path.abspath(self.source)) if self.source else os.getcwd()
        path = name if os.path.isabs(name) else os.path.join(base, name)
        why = _thm_read_weight_table(path)
        if why:
            raise ValueError(f"<thm> {key}: {why}")

    def _thm_write(self, s):
        """Write settings ``s`` back into the block, the GUI's way."""
        if s == self._thm_settings():
            return                               # untouched: keep it verbatim
        loc = self._thm_locate()
        if _thm_is_default(s):
            if loc is not None:                  # all default: no block
                attr, lines, i, j = loc
                setattr(self, attr, "".join(lines[:i] + lines[j + 1:]))
            return
        body = _thm_compose(self._thm_body(), s)
        if loc is not None:
            attr, lines, i, j = loc
            nl = "\r\n" if lines[i].endswith("\r\n") else "\n"
            new = [ln + nl for ln in body]
            setattr(self, attr, "".join(lines[:i + 1] + new + lines[j:]))
            return
        block = ["<thm>"] + body + ["</thm>"]
        lines = self._suffix.splitlines()
        if "</targetInt>" in lines:
            end = lines.index("</targetInt>") + 1
            self._set_suffix_lines(lines[:end] + block + lines[end:])
        else:
            self._set_suffix_lines(lines + block)

    # -- rendering ------------------------------------------------------------

    def __str__(self):
        head = f"AzrModel"
        if self.source:
            head += f"  [{os.path.basename(self.source)}]"
        head += f"   {len(self.levels)} levels"
        lines = [head, "=" * len(head)]
        seen = []
        for lv in self.levels:
            if lv.jpi not in seen:
                seen.append(lv.jpi)
                lines.append(f"\nJ^pi = {lv.jpi}")
            fix = "fixed" if lv.fixed else "FREE"
            lines.append(f"  E = {lv.energy:>9.4g} MeV  [{fix}]")
            for c in lv.channels:
                kind = f"{'photon' if c.is_photon else 'particle'} pair{c.pair}"
                cfix = "fixed" if c.channel_fixed else "FREE"
                lines.append(f"      {kind:<16} L={c.L} S={c.S:g}  "
                             f"Gamma={c.gamma:>11.5g}  [{cfix}]")
        return "\n".join(lines)

    def __repr__(self):
        return f"AzrModel({len(self.levels)} levels, source={self.source!r})"
