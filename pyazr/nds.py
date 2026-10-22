"""Fetch nuclear data from IAEA EXFOR and NDS (LiveChart/ENSDF) in
AZURE2-friendly form.

The two web services wrapped here:

- **EXFOR** (``https://nds.iaea.org/exfor/``): experimental cross sections,
  differential cross sections, analyzing powers, yields, and more, with the
  bibliographic record for each entry.  ``search_exfor`` / ``fetch_exfor``
  wrap the ``x4list`` / ``x4get`` Web-API.
- **LiveChart/ENSDF** (``https://www-nds.iaea.org/relnsd/v1/data``): evaluated
  level schemes, gamma transitions and ground-state properties.  ``fetch_levels``
  / ``fetch_gammas`` / ``fetch_ground_state`` wrap it.  ``fetch_levels`` falls
  back to the NNDC ENSDF API (``https://www.nndc.bnl.gov/ensdf-api/``, the same
  ENSDF adopted levels as JSON) when LiveChart cannot be reached.

``ExforData.to_azr`` converts an EXFOR dataset into an AZURE2 data file
(``energy angle crossSection error``, lab frame) and returns the keyword
arguments for :meth:`pyazr.AzrModel.add_data_segment`, so fetched data can be
dropped straight into an evaluation::

    data = fetch_exfor("O2599004")                      # 13C(p,g)14N S-factor
    kw = data.to_azr("run_dir/data", entrance=1, exit=2,
                     observable="total-capture")
    AzrModel.from_file("13N.azr").add_data_segment(**kw).write("13N_new.azr")

When the EXFOR Web-API cannot be reached (it is intermittently behind a
Cloudflare challenge), ``fetch_exfor`` and ``reference`` read the entry from a
local copy of the EXFOR master files instead (https://github.com/IAEA-NDS/exfor_master,
``git clone --depth 1``, about 1.5 GB), found at ``$AZURE2_EXFOR_MASTER`` or
``~/.azure2/exfor_master``.  ``source="local"`` reads only the copy,
``source="web"`` only the server.  ``search_exfor`` has no local fallback.

References resolve to a DOI through the CrossRef API, so the paper behind any
dataset can be pulled up (``reference`` + ``resolve_doi``).

This is the one EXFOR client in the package.  ``gui/src/ExforData.cpp`` is its
counterpart inside AZURESetup; the two are independent implementations of the
same Web-API, and a parsing rule learned by either belongs in both.

Command line, for a quick look without writing a script::

    python -m pyazr.nds search --target C-13 --reaction p,g --quantity SIG
    python -m pyazr.nds download O2599004 -o data/roughton.dat
    python -m pyazr.nds download O2599004 --source local
"""

import csv
import dataclasses
import io
import json
import math
import os
import re
import urllib.error
import urllib.parse
import urllib.request
import warnings

import numpy as np

EXFOR_URL = "https://nds.iaea.org/exfor/"
LIVECHART_URL = "https://www-nds.iaea.org/relnsd/v1/data"
ENSDF_API_URL = "https://www.nndc.bnl.gov/ensdf-api"
CROSSREF_URL = "https://api.crossref.org/works"
UA = {"User-Agent": "pyazr/nds (research tool)"}

#: Projectile shorthand used by EXFOR reaction strings -> mass number.
#: ``(Z, A)`` of the light particles EXFOR names with a shorthand instead of a
#: ``Z-SYM-A`` code.  The charge matters as much as the mass: an S-factor
#: conversion divides by ``exp(2 pi eta)`` with ``eta ~ Z1 Z2``, so guessing
#: ``Z1`` wrong by a factor of two is wrong by ``exp(2 pi eta)`` -- a factor of
#: 5000 on 3He(a,g)7Be at 93 keV, verified against the Costantini 2008 cross
#: sections for the same LUNA measurement.
_PROJ_ZA = {
    "P": (1, 1), "N": (0, 1), "D": (1, 2), "T": (1, 3), "G": (0, 0),
    "A": (2, 4), "ALPHA": (2, 4), "4-HE": (2, 4), "HE4": (2, 4),
    "HE3": (2, 3), "3-HE": (2, 3), "H1": (1, 1), "H2": (1, 2), "H3": (1, 3),
    "E": (-1, 0), "E-": (-1, 0), "E+": (1, 0),
}

_PROJ_MASS = {k: v[1] for k, v in _PROJ_ZA.items()}

#: Assumed relative uncertainty for a point EXFOR gives no error for.  A zero
#: error would be read by AZURE2 as an infinitely precise measurement.
_DEFAULT_REL_ERROR = 0.05

# Constants matching include/Constants.h, so a converted S-factor or Rutherford
# ratio reproduces what AZURE2 computes internally.
_UCONV = 931.4940880
_FSTRUC = 1.00 / 137.0359996790
_HBARC = 197.32696310


def _http_get(url, params=None, timeout=60):
    """GET ``url`` and return the body as text.

    urllib rather than requests: pyazr declares numpy, mpmath and scipy, and a
    data-fetching convenience is not worth adding a dependency the rest of the
    package does not need.
    """
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


@dataclasses.dataclass
class ExforDataset:
    """One dataset found by :func:`search_exfor`."""

    dataset_id: str
    reaction: str        # EXFOR reaction code, e.g. "6-C-13(P,G)7-N-14,,SIG,,SFC"
    npoints: int
    en_min: float        # eV
    en_max: float        # eV
    an_min: float        # deg (angle grid), None if not angular
    an_max: float        # deg
    author: str
    reference: str       # journal citation string
    year: int

    def __repr__(self):  # pragma: no cover - display helper
        return (f"<ExforDataset {self.dataset_id} {self.reaction} "
                f"n={self.npoints} E={self.en_min:.3g}-{self.en_max:.3g} eV"
                + (f" A={self.an_min:.0f}-{self.an_max:.0f}" if self.an_min else "")
                + f" {self.author} {self.year}>")


@dataclasses.dataclass
class ExforData:
    """A parsed EXFOR dataset.

    ``columns`` is a list of :class:`ExforColumn` (name + unit as reported).
    ``arrays`` holds the numeric columns keyed by their plain names (``DATA``,
    ``ERR-S``, ``ERR-SYS``, ``EN``/``EN-CM``, ``ANG``/``ANG-CM``).
    """

    dataset_id: str
    year: int
    author: str
    reaction: str            # Reacode, e.g. "6-C-13(P,G)7-N-14,,SIG,,SFC"
    columns: list
    arrays: dict
    raw: str                 # the raw CSV text

    # -- reaction parsing ----------------------------------------------------
    @property
    def projectile(self):
        """Projectile of the reaction, as EXFOR names it."""
        m = re.search(r"\(([^)]+)\)", self.reaction)
        if not m:
            return ""
        return m.group(1).split(",")[0]

    @property
    def target(self):
        """Target nucleus."""
        m = re.match(r"([0-9]+-[A-Z]+-[0-9]+)", self.reaction)
        return m.group(1) if m else ""

    @property
    def exit(self):
        """Outgoing particle."""
        m = re.search(r"\(([^)]+)\)", self.reaction)
        if not m:
            return ""
        parts = m.group(1).split(",")
        return parts[1] if len(parts) > 1 else ""

    def za(self, name):
        """``(Z, A)`` of a nucleus named like ``6-C-13`` / ``1-H-1``.

        Also resolves the light-particle shorthands EXFOR uses in the reaction
        parentheses -- ``P``, ``N``, ``D``, ``T``, ``A``, ``HE3``, ``G`` --
        which the ``Z-SYM-A`` pattern does not match.  Returns ``None`` only
        for a name that is neither.
        """
        if not name:
            return None
        m = re.match(r"([0-9]+)-([A-Z]+)-([0-9]+)", name)
        if m:
            return int(m.group(1)), int(m.group(3))
        return _PROJ_ZA.get(name.strip().upper())

    def masses(self, target_mass=None, projectile_mass=None):
        """Mass numbers for the entrance pair, in ``(m_target, m_proj)`` u."""
        if target_mass is None:
            za = self.za(self.target)
            target_mass = za[1] if za else None
        if projectile_mass is None:
            proj = self.projectile
            p = self.za(proj)
            if p:
                projectile_mass = p[1]
            else:
                projectile_mass = _PROJ_MASS.get(proj)
        if target_mass is None or projectile_mass is None:
            raise ValueError(
                "cannot determine entrance masses from reaction "
                f"{self.reaction!r}; pass target_mass= / projectile_mass=")
        return float(target_mass), float(projectile_mass)

    # -- conversion ----------------------------------------------------------
    def to_azr(self, data_dir, entrance, exit, observable=None,
               target_mass=None, projectile_mass=None,
               cross_section_scale=1.0, energy_scale=1.0, angle_scale=1.0,
               file_name=None, rutherford=None):
        """Write the dataset as an AZURE2 data file and return
        :meth:`~pyazr.AzrModel.add_data_segment` keyword arguments.

        Returns a dict with ``data_file``, ``entrance``, ``exit``,
        ``observable``, ``energy_min``/``energy_max`` (lab MeV),
        ``angle_min``/``angle_max`` (deg) and ``norm_error`` (the EXFOR
        ``ERR-SYS`` percent when present), ready to pass to
        ``AzrModel.add_data_segment(**kw)``.

        **Frames.** AZURE2 data files are always **lab** energy.  ``EN-CM``
        columns are converted with the standard ``E_lab = E_cm (m_t + m_p)/m_t``
        factor; ``EN`` columns are taken as lab as-is.  The angle column is
        written through unchanged, so pick the ``observable`` that matches its
        frame: EXFOR ``ANG-CM`` data belongs with ``observable="differential-cm"`
        (or ``"analyzing-power"``), EXFOR ``ANG`` (lab) data with
        ``observable="differential"``.

        **S-factor data.** When ``DATA`` is reported as ``B*KEV`` / ``B*EV``
        (an ``SFC`` reaction), the values are converted to barns using
        :func:`sfactor_to_cross_section`, so AZURE2 receives a cross section.
        Set ``observable`` to ``total-capture`` for a radiative-capture
        reaction.

        **Ratio-to-Rutherford data** (EXFOR ``,,RTH``) is dimensionless and
        indistinguishable from an analyzing power once fetched, because
        ``x4get`` drops the quantity suffix that says which.  ``rutherford=True``
        multiplies by the Coulomb cross section to give barn/sr; ``False``
        passes the values through; the default warns when the choice matters.

        Other unit conversions (``MB``, ``NB/SR``, ...) are applied from the
        column unit; ``cross_section_scale`` / ``energy_scale`` / ``angle_scale``
        multiply on top for any manual fix-up.
        """
        os = _import_os()
        data = _azr_columns(self, target_mass=target_mass,
                            projectile_mass=projectile_mass,
                            cross_section_scale=cross_section_scale,
                            energy_scale=energy_scale, angle_scale=angle_scale,
                            rutherford=rutherford)
        os.makedirs(data_dir, exist_ok=True)
        if file_name is None:
            file_name = _safe_name(self.dataset_id, self.author, self.year)
        path = os.path.join(data_dir, file_name + ".dat")
        with open(path, "w") as f:
            for row in data:
                f.write("  ".join(f"{v:.8g}" for v in row) + "\n")
        obs = observable or self._guess_observable()
        ang = data[:, 1]
        en = data[:, 0]
        syspct = self._sys_error_percent()
        kw = dict(data_file=path, entrance=int(entrance), exit=int(exit),
                  observable=obs,
                  energy_min=float(en.min()), energy_max=float(en.max()),
                  angle_min=float(ang.min()), angle_max=float(ang.max()))
        if syspct is not None:
            kw["norm_error"] = float(syspct)
        return kw

    # -- internal helpers ----------------------------------------------------
    def _sys_error_percent(self):
        arr = self.arrays.get("ERR-SYS")
        if arr is None or len(arr) == 0:
            return None
        unit = None
        for c in self.columns:
            if c.name == "ERR-SYS":
                unit = c.unit
        return float(arr[0]) if unit == "PER-CENT" else None

    def _guess_observable(self):
        if self.exit == "G" or "SFC" in self.reaction:
            return "total-capture"
        if "POL" in self.reaction or "NO-DIM" in [c.unit for c in self.columns]:
            return "analyzing-power"
        if any(c.name.startswith("ANG") for c in self.columns):
            return "differential-cm"
        return "angle-integrated"


@dataclasses.dataclass
class ExforColumn:
    """One CSV column: plain ``name`` plus the unit reported in parentheses."""

    name: str
    unit: str = ""

    def __repr__(self):  # pragma: no cover
        return f"<ExforColumn {self.name} ({self.unit})>" if self.unit \
            else f"<ExforColumn {self.name}>"


def _import_os():
    import os
    return os


def _safe_name(dataset_id, author, year):
    stem = re.sub(r"[^A-Za-z0-9]+", "_", author or "").strip("_")
    return f"{stem}_{year}_{dataset_id}" if stem else f"{dataset_id}"


# -- units ------------------------------------------------------------------
#
# EXFOR spells its units many ways for the same quantity (B, MB, MICRO-B, MU-B,
# NB/SR, B*KEV, ...), so these match on prefix rather than looking a fixed
# string up in a table: an unlisted spelling used to fall through to "assume
# it is already in the right unit", which is silently wrong by orders of
# magnitude.

def _energy_factor(unit):
    """EXFOR energy unit -> factor to MeV.  eV is the computational default."""
    unit = (unit or "").upper()
    if "GEV" in unit:
        return 1e3
    if "MEV" in unit:
        return 1.0
    if "KEV" in unit:
        return 1e-3
    if "MILLI-EV" in unit or "MILLIEV" in unit:
        return 1e-9
    return 1e-6


def _barn_factor(unit):
    """EXFOR cross-section unit -> factor to barn (or barn/sr)."""
    unit = (unit or "").upper()
    # Longest-first, so MICRO-B is not read as M-then-B.  A microbarn is 1e-6
    # barn however it is spelled -- the old table had UB/SR at 1e-12, which is
    # a factor of a million.
    for prefix, factor in (("MICRO-B", 1e-6), ("MU-B", 1e-6), ("MUB", 1e-6),
                           ("UB", 1e-6), ("MB", 1e-3), ("KB", 1e3),
                           ("NB", 1e-9), ("PB", 1e-12), ("FB", 1e-15),
                           ("B", 1.0)):
        if unit.startswith(prefix):
            return factor
    return 1.0


def _sfactor_factors(unit):
    """``(factor to barn*MeV, energy factor to MeV)`` for an S-factor unit.

    ``None`` when ``unit`` is not an S-factor (barn times energy) unit.
    """
    unit = (unit or "").upper()
    if "*" not in unit or "EV" not in unit:
        return None
    return _barn_factor(unit) * _energy_factor(unit), _energy_factor(unit)


def _rutherford_barn_per_sr(z1, z2, e_cm_mev, angle_deg):
    """Rutherford differential cross section, c.m. frame, in barn/sr.

    ``dsigma/dOmega = (Z1 Z2 e^2 / 4E)^2 / sin^4(theta/2)`` with
    ``e^2 = alpha hbar c`` in MeV fm, giving fm^2/sr; 1 fm^2 = 1e-2 barn.
    """
    e_cm_mev = np.asarray(e_cm_mev, float)
    sin_half = np.sin(np.radians(np.asarray(angle_deg, float)) / 2.0)
    amplitude = (z1 * z2 * _FSTRUC * _HBARC) / (4.0 * e_cm_mev)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = (amplitude ** 2) / (sin_half ** 4) * 1e-2
    return np.where((e_cm_mev > 0) & (sin_half != 0), out, np.nan)


def _parse_columns(header):
    columns = []
    for raw in header:
        name = raw
        unit = ""
        m = re.search(r"\(([^)]+)\)", raw)
        if m:
            unit = m.group(1)
            name = raw[:m.start()].strip().rstrip()
        columns.append(ExforColumn(name, unit))
    return columns


def search_exfor(target=None, reaction=None, quantity=None, author=None,
                 accnum=None, limit=50):
    """Search EXFOR for datasets.

    ``target`` is the EXFOR nucleus name (``"C-13"``), ``reaction`` the
    projectile/exit pair (``"p,g"``, ``"p,el"``), ``quantity`` one of the
    EXFOR quantity codes (``"SIG"`` total cross section, ``"DA"`` differential
    cross section, ``"pol"`` analyzing power, ``"FY"`` yields, ...).  All
    filters are optional; give at least one.  Returns a list of
    :class:`ExforDataset`.
    """
    params = {"json": ""}
    if target:
        params["Target"] = target
    if reaction:
        params["Reaction"] = reaction
    if quantity:
        params["Quantity"] = quantity
    if author:
        params["Author"] = author
    if accnum:
        params["AccessNumber"] = accnum
    body = _http_get(EXFOR_URL + "x4list", params)
    # x4list emits a missing value for zero-point datasets -- `"enMin":,` --
    # which a strict parser rejects, failing the whole search because of one
    # bad entry.  Substituting null salvages the rest.
    repaired = re.sub(r'":\s*(?=[,}])', '":null', body)
    try:
        payload = json.loads(repaired)
    except ValueError:
        raise ValueError(f"EXFOR x4list returned non-JSON:\n{body[:300]}")
    out = []
    for d in payload.get("x4Datasets", []):
        out.append(ExforDataset(
            dataset_id=str(d.get("id", "")),
            reaction=str(d.get("RC", "")),
            npoints=int(d.get("npts", 0)),
            en_min=float(d.get("enMin", 0) or 0),
            en_max=float(d.get("enMax", 0) or 0),
            an_min=float(d.get("anMin", 0) or 0),
            an_max=float(d.get("anMax", 0) or 0),
            author=str(d.get("A1", "")),
            reference=str(d.get("ref", "")),
            year=int(d.get("year1", 0) or 0),
        ))
        if len(out) >= limit:
            break
    return out


def fetch_exfor(dataset_id, plus=0, source="auto"):
    """Fetch one EXFOR dataset by its ID (``"O2599004"``; a pointer dataset as
    ``"O2599004.1"``).

    ``plus=0`` returns the original units (``B*KEV`` S-factors, ``NB/SR``
    differentials, ...); ``plus=1`` the computational units (eV, ``B*EV``,
    ``B/SR``); ``plus=2`` the universal grid (residual/incident/sec energies,
    angle).  Returns an :class:`ExforData`.

    ``source="auto"`` asks the Web-API and, if it cannot be reached or answers
    with something other than data (a Cloudflare challenge page), reads the
    local EXFOR master files (:func:`exfor_master_path`) with a warning.
    ``"web"`` and ``"local"`` use only one of them.  The local copy gives the
    original units only (``plus=0``); its columns and units are those of the
    EXFOR entry, as ``x4get`` reports them.
    """
    if source not in ("auto", "web", "local"):
        raise ValueError(f"source must be 'auto', 'web' or 'local', not {source!r}")
    if source == "local":
        return _fetch_exfor_local(dataset_id, plus)
    try:
        return _fetch_exfor_web(dataset_id, plus)
    except (OSError, _WebUnavailable) as err:   # URLError and timeouts are OSErrors
        if source == "web" or exfor_master_path() is None:
            raise
        warnings.warn(f"EXFOR Web-API unavailable for {dataset_id} ({err}); "
                      f"reading the local EXFOR master files at "
                      f"{exfor_master_path()}.")
        return _fetch_exfor_local(dataset_id, plus)


class _WebUnavailable(Exception):
    """The server answered, but with a page instead of data (a Cloudflare
    challenge or a maintenance page)."""


def _web_text(url, params):
    text = _http_get(url, params)
    head = text.lstrip()[:200].lower()
    if head.startswith(("<!doctype html", "<html")):
        title = re.search(r"<title>(.*?)</title>", text, re.I | re.S)
        raise _WebUnavailable(f"HTML page instead of data"
                              + (f": {title.group(1).strip()!r}" if title else ""))
    return text


def _fetch_exfor_web(dataset_id, plus=0):
    params = {"DatasetID": dataset_id, "op": "csv"}
    if plus:
        params["plus"] = int(plus)
    text = _web_text(EXFOR_URL + "x4get", params)
    if text.strip().startswith("Error"):
        raise ValueError(f"EXFOR x4get {dataset_id}: {text.strip()[:200]}")
    reader = csv.reader(io.StringIO(text))
    try:
        header = next(reader)
    except StopIteration:
        raise ValueError(f"EXFOR x4get {dataset_id} returned no rows.")
    columns = _parse_columns(header)
    rows = [row for row in reader if row and row[0].strip()]
    arrays, meta = {}, {}
    for i, c in enumerate(columns):
        vals = [row[i] for row in rows if i < len(row)]
        if vals and _is_numeric(vals):
            arrays[c.name] = np.array([float(v) for v in vals])
            meta[c.name] = vals[0]
        else:
            meta[c.name] = vals[0] if vals else ""
    reaction = str(meta.get("Reacode", ""))
    year = int(float(meta.get("year1", 0) or 0))
    author = str(meta.get("author1", ""))
    if "DATA" not in arrays and "DATA-CM" not in arrays:
        server_msg = meta.get("___01d___wrong request___") or ""
        if "wrong request" in str(meta):
            server_msg = str(meta["___01d___wrong request___"])
        raise ValueError(f"EXFOR x4get {dataset_id}: no DATA column; server "
                         f"said: {server_msg or 'unknown'}")
    return ExforData(dataset_id=dataset_id, year=year, author=author,
                     reaction=reaction, columns=columns, arrays=arrays,
                     raw=text)


# -- local EXFOR master files ------------------------------------------------
#
# The EXFOR exchange format (LEXFOR / the EXFOR Formats Manual): 80-column
# records; a BIB keyword in columns 1-10, its pointer in column 11 and the text
# in 12-66; COMMON and DATA sections give a heading line and a unit line of
# 11-column fields, six per record, then the values, a row continuing on the
# next record after six fields.  Columns 67-80 (record identification) are not
# used.  Numbers may be written without the E (``1.23-4``), and a blank field is
# a missing value.

def exfor_master_path():
    """The local EXFOR master files (the ``exfor_master`` repository, whose
    entries are ``exforall/<first 3 characters>/<entry>.x4``):
    ``$AZURE2_EXFOR_MASTER``, else ``~/.azure2/exfor_master``; ``None`` if
    neither holds an ``exforall`` directory."""
    for candidate in (os.environ.get("AZURE2_EXFOR_MASTER"),
                      os.path.join(os.path.expanduser("~"), ".azure2", "exfor_master")):
        if candidate and os.path.isdir(os.path.join(candidate, "exforall")):
            return candidate
    return None


def _local_entry_text(entry):
    """The raw text of EXFOR entry ``entry`` (5 characters) from the local copy."""
    root = exfor_master_path()
    if root is None:
        raise FileNotFoundError(
            "no local EXFOR master files: clone https://github.com/IAEA-NDS/exfor_master "
            "and point $AZURE2_EXFOR_MASTER (or ~/.azure2/exfor_master) at it.")
    entry = entry.upper()
    path = os.path.join(root, "exforall", entry[:3], entry + ".x4")
    if not os.path.isfile(path):
        raise ValueError(f"EXFOR entry {entry} is not in the local copy ({path}).")
    with open(path, encoding="utf-8", errors="replace") as f:
        return f.read()


def _x4_number(field):
    """One EXFOR numeric field; NaN when blank.  Accepts the Fortran-style
    exponent without E (``1.23-4``, ``-5.6+2``) and blanks inside the field,
    which the Fortran reading of the format ignores (``' 1.44   -06'``)."""
    t = "".join(field.split())
    if not t:
        return math.nan
    try:
        return float(t)
    except ValueError:
        pass
    m = re.fullmatch(r"([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)", t)
    if m:
        return float(f"{m.group(1)}e{m.group(2)}")
    raise ValueError(f"not an EXFOR number: {field!r}")


def _x4_subentries(text):
    """{subentry number (3 characters): its records} for one entry."""
    subs, current, key = {}, None, None
    for line in text.splitlines():
        tag = line[:10].strip()
        if tag == "SUBENT":
            m = re.search(r"[0-9A-Z]{8}", line[10:33])
            key = m.group(0)[-3:] if m else None
            current = []
        elif tag in ("ENDSUBENT", "NOSUBENT"):
            if current is not None:
                subs[key] = current
            current = None
        elif current is not None:
            current.append(line.rstrip("\n"))
    return subs


def _x4_bib(records):
    """{keyword: [(pointer, text line), ...]} from a subentry's BIB section."""
    bib, keyword, inside = {}, None, False
    for line in records:
        tag = line[:10].strip()
        if tag == "BIB":
            inside = True
            continue
        if tag == "ENDBIB":
            break
        if not inside:
            continue
        if tag:
            keyword = tag
        if keyword:
            bib.setdefault(keyword, []).append((line[10:11].strip(), line[11:66].rstrip()))
    return bib


def _x4_table(records, section):
    """The COMMON or DATA section of a subentry: ``(headings, units, rows)``,
    headings as (name, pointer) and rows as lists of floats; None if absent."""
    for i, line in enumerate(records):
        if line[:10].strip() != section:
            continue
        nfields = int(line[11:22])
        nrows = int(line[22:33]) if section == "DATA" else 1
        per = (nfields + 5) // 6                    # records per heading, unit or data row
        def fields(start):
            text = "".join(records[start + k][:66].ljust(66) for k in range(per))
            return [text[11 * j:11 * (j + 1)] for j in range(nfields)]
        heads = [(f[:10].strip(), f[10:11].strip()) for f in fields(i + 1)]
        units = [f.strip() for f in fields(i + 1 + per)]
        rows = [[_x4_number(v) for v in fields(i + 1 + per * (2 + r))]
                for r in range(nrows)]
        return heads, units, rows
    return None


def _x4_reaction(bib, pointer):
    """The REACTION code for ``pointer`` ('' = the first), without the outer
    parentheses -- as ``x4get`` reports it in ``Reacode``.  A code continues on
    the next lines while its parentheses are open (a combination of reactions
    can start a continuation line with ``(``); the free text after it is
    dropped."""
    codes, depth = [], 0                            # [pointer, code text]
    for ptr, text in bib.get("REACTION", []):
        t = text.strip()
        if depth > 0 and codes:
            part = t
        elif ptr or t.startswith("("):
            codes.append([ptr, ""])
            part = t
        else:
            continue                                # free text
        for k, ch in enumerate(part):
            depth += (ch == "(") - (ch == ")")
            codes[-1][1] += ch
            if depth == 0:                          # the code closed: the rest is free text
                break
    for ptr, code in codes:
        if ptr == pointer or not pointer:
            return code[1:-1] if code.startswith("(") and code.endswith(")") else code
    return ""


def _fetch_exfor_local(dataset_id, plus=0):
    if plus:
        raise ValueError(f"{dataset_id}: the local EXFOR master files give the "
                         f"original units only (plus=0), not plus={plus}.")
    m = re.fullmatch(r"([0-9A-Z]{5})([0-9]{3})(?:\.?([0-9A-Z]))?", dataset_id.strip().upper())
    if not m:
        raise ValueError(f"not an EXFOR dataset ID: {dataset_id!r} "
                         f"(entry + subentry, e.g. 'O2599004', or 'O2599004.1' for a pointer)")
    entry, number, pointer = m.group(1), m.group(2), m.group(3) or ""
    text = _local_entry_text(entry)
    subs = _x4_subentries(text)
    if number not in subs:
        raise ValueError(f"{dataset_id}: entry {entry} has no subentry {number} "
                         f"(it has {', '.join(sorted(subs))}).")
    records = subs[number]
    data = _x4_table(records, "DATA")
    if data is None:
        raise ValueError(f"{dataset_id}: subentry {number} has no DATA section.")
    heads, units, rows = data
    pointers = sorted({p for name, p in heads if p})
    if pointers and pointer not in pointers:
        raise ValueError(f"{dataset_id}: subentry {number} holds several datasets "
                         f"(pointers {', '.join(pointers)}); ask for one, e.g. "
                         f"'{entry}{number}.{pointers[0]}'.")
    columns, arrays = [], {}
    def add(name, unit, values):
        if name in arrays or not np.any(np.isfinite(values)):
            return
        columns.append(ExforColumn(name, unit))
        arrays[name] = values
    nrow = len(rows)
    for j, ((name, ptr), unit) in enumerate(zip(heads, units)):
        if ptr and ptr != pointer:
            continue
        add(name, unit, np.array([r[j] for r in rows], dtype=float))
    # COMMON values of this subentry, then of subentry 001, hold for every row
    # (x4get repeats them in each line of its CSV).
    for records_common in (records, subs.get("001", [])):
        common = _x4_table(records_common, "COMMON")
        if common is None:
            continue
        for ((name, ptr), unit, value) in zip(common[0], common[1], common[2][0]):
            if ptr and ptr != pointer:
                continue
            add(name, unit, np.full(nrow, value))
    bib = _x4_bib(records)
    bib1 = _x4_bib(subs.get("001", []))
    reaction = _x4_reaction(bib, pointer) or _x4_reaction(bib1, pointer)
    authors = " ".join(t for _, t in bib1.get("AUTHOR", [])).strip()
    author = authors.strip("()").split(",")[0].strip()
    year = 0
    for _, t in bib1.get("REFERENCE", []) + bib.get("REFERENCE", []):
        for part in re.split(r"[,()]", t):
            part = part.strip()
            if re.fullmatch(r"(19|20)\d\d", part):
                year = int(part)
                break
            if re.fullmatch(r"(19|20)\d{4,6}", part):   # a date, 19870515
                year = int(part[:4])
                break
        if year:
            break
    if "DATA" not in arrays and "DATA-CM" not in arrays:
        raise ValueError(f"{dataset_id}: no DATA column in the local entry "
                         f"(headings {[h for h, _ in heads]}).")
    return ExforData(dataset_id=dataset_id, year=year, author=author,
                     reaction=reaction, columns=columns, arrays=arrays,
                     raw="\n".join(records))


def _is_numeric(vals):
    try:
        for v in vals:
            float(v)
    except (TypeError, ValueError):
        return False
    return True


# -- unit/frame conversion helpers ------------------------------------------

def cm_to_lab(E_cm, m_target, m_projectile):
    """Center-of-mass energy (MeV) -> lab energy (MeV), non-relativistic."""
    return E_cm * (m_target + m_projectile) / m_target


def sfactor_to_cross_section(S, E_cm_kev, z1, z2, mu):
    """Convert an S-factor to a cross section: ``sigma = S exp(-2 pi eta) / E``.

    ``S`` in b keV, ``E_cm_kev`` in keV, ``mu`` the reduced mass in u,
    ``z1``/``z2`` the projectile/target charges.  The Sommerfeld parameter uses
    the same constant as AZURE2's ``EPoint::CalcEDependentValues``
    (``pi * sqrt(uconv/2) * fstruct = 0.157488``), so feeding AZURE2 the
    resulting barns reproduces the reported S-factor.
    """
    E_cm_MeV = E_cm_kev / 1e3
    eta = 0.157488 * z1 * z2 * np.sqrt(mu / E_cm_MeV)
    return S * np.exp(-2 * np.pi * eta) / E_cm_kev


def _azr_columns(data, target_mass=None, projectile_mass=None,
                 cross_section_scale=1.0, energy_scale=1.0, angle_scale=1.0,
                 rutherford=None):
    """Build the ``(lab E_MeV, angle_deg, cs, err)`` array for a dataset."""
    m_t, m_p = data.masses(target_mass, projectile_mass)
    cols = {c.name: c.unit for c in data.columns}
    energy = _energy_column(data, cols, m_t, m_p) * energy_scale
    angle = _angle_column(data, cols) * angle_scale
    cs, err = _data_column(data, cols, m_t, m_p,
                           cross_section_scale=cross_section_scale,
                           rutherford=rutherford)
    out = np.column_stack([energy, angle, cs, err])
    # A row with a blank energy, angle or value (EXFOR leaves fields blank, and
    # the local master files read them as NaN) cannot go into an AZURE2 file.
    keep = np.all(np.isfinite(out), axis=1)
    if not np.all(keep):
        warnings.warn(f"{data.dataset_id}: {int(np.sum(~keep))} of {len(out)} rows "
                      f"have a blank energy, angle or value and were left out.")
        out = out[keep]
    return out


def _energy_column(data, cols, m_t, m_p):
    name = None
    for c in data.columns:
        if c.name in ("EN-CM", "EN"):
            name = c.name
    if name is None:
        raise ValueError(f"{data.dataset_id}: no EN / EN-CM energy column in "
                         f"{list(cols)}.")
    factor = _energy_factor(cols[name])
    # An S-factor unit such as B*KEV is a unit of the *product* barn x energy,
    # and its energy part is chosen so the tabulated S values read well.  It
    # says nothing about the unit of the abscissa, so it must not override it.
    # Checked against the independent lab ranges x4list reports: C1610002
    # (Brown 2007, EN-CM/MEV with DATA/B*KEV), O2521002 (Piatti 2020,
    # EN-CM/KEV with DATA/B*EV) and O1590003 (Cruz 2008, EN/KEV) all have the
    # declared column unit right and the S-factor unit different.  Disagreement
    # is still worth a warning, because a genuine EXFOR mislabel is a factor of
    # 1000 either way; pass energy_scale= to correct one.
    sf = _sfactor_factors(cols.get("DATA-CM") or cols.get("DATA"))
    if sf is not None and sf[1] != factor:
        warnings.warn(
            f"{data.dataset_id}: the {name} column says {cols[name]!r} while "
            f"the S-factor unit's energy part implies {sf[1]:g} MeV per unit. "
            f"Using the declared {name} unit, which is the authority for this "
            f"column; check the energies against the reference if the dataset "
            f"looks displaced.")
    E = data.arrays[name].astype(float) * factor
    if name == "EN-CM":
        E = cm_to_lab(E, m_t, m_p)
    return E


def _angle_column(data, cols):
    """Angle of each point in degrees; zeros for angle-integrated data."""
    for c in data.columns:
        if c.name in ("ANG-CM", "ANG"):
            return data.arrays[c.name].astype(float)
        if c.name.startswith("COS"):
            # A cosine grid is the same axis in another variable.
            cosine = np.clip(data.arrays[c.name].astype(float), -1.0, 1.0)
            return np.degrees(np.arccos(cosine))
    return np.zeros(len(_values(data)))


def _values(data):
    return data.arrays["DATA-CM" if "DATA-CM" in data.arrays else "DATA"]


def _error_column(data, cols, raw, data_unit):
    """Uncertainty per point, expressed in the *raw* units of ``raw``.

    Raw, not barns: the caller converts the data column afterwards -- and for
    an S-factor or a Rutherford ratio that conversion is per-point and
    non-linear -- so the error has to go through exactly the same steps.

    EXFOR gives the uncertainty absolutely (``ERR-S``, ``DATA-ERR``) or as a
    percentage (often ``ERR-T``), and sometimes not at all.  A point with no
    error must not reach AZURE2 as zero, which reads as an infinitely precise
    measurement and would dominate the chi-squared.
    """
    for name in ("ERR-S", "DATA-ERR", "ERR-T"):
        arr = data.arrays.get(name)
        if arr is None:
            continue
        unit = (cols.get(name) or "").upper()
        err = arr.astype(float)
        if "PER-CENT" in unit or "PERCENT" in unit:
            err = np.abs(raw) * err / 100.0
        elif unit and unit != data_unit and "NO-DIM" not in unit:
            # Quoted in a different unit from the data; restate it in the
            # data's, so the one conversion below serves both.
            err = err * _barn_factor(unit) / (_barn_factor(data_unit) or 1.0)
        # A blank field (NaN from the local master files) counts as missing too.
        given = np.isfinite(err) & (err != 0.0)
        if np.any(given):
            return np.where(given, err, _DEFAULT_REL_ERROR * np.abs(raw))
    return _DEFAULT_REL_ERROR * np.abs(raw)


def _is_rutherford_ratio(data, cols, rutherford=None):
    """Is this dimensionless differential data a ratio to Rutherford?

    ``rutherford=True``/``False`` answers it outright.  Left at ``None`` the
    reaction code decides -- but only ``x4list`` reports the ``,,RTH`` suffix;
    ``x4get`` drops it, so a dataset fetched by ID alone can be a Rutherford
    ratio with nothing in it to say so.  A dimensionless *differential* is then
    genuinely ambiguous -- an analyzing power looks exactly the same -- and
    guessing either way is silently wrong for the other, so this warns and
    leaves the values alone.
    """
    if rutherford is not None:
        return bool(rutherford)
    unit = (cols.get("DATA-CM") or cols.get("DATA") or "").upper()
    if "NO-DIM" not in unit:
        return False
    if "RTH" in data.reaction.upper():
        return True
    is_polarization = "POL" in data.reaction.upper()
    has_angles = any(c.name.startswith(("ANG", "COS")) for c in data.columns)
    if has_angles and not is_polarization:
        warnings.warn(
            f"{data.dataset_id}: dimensionless differential data, which is "
            f"either a ratio to Rutherford or an analyzing power -- x4get does "
            f"not say which (x4list does: check the reaction code for ',,RTH'). "
            f"Values passed through unconverted; pass rutherford=True to "
            f"multiply by the Coulomb cross section, or rutherford=False to "
            f"silence this.")
    return False


def _data_column(data, cols, m_t, m_p, cross_section_scale=1.0,
                 rutherford=None):
    dname = "DATA-CM" if "DATA-CM" in cols else "DATA"
    if dname not in data.arrays:
        raise ValueError(f"{data.dataset_id}: no DATA / DATA-CM column in "
                         f"{list(cols)}.")
    unit = (cols[dname] or "").upper()
    cs = data.arrays[dname].astype(float)
    stat = _error_column(data, cols, cs, unit)

    sf = _sfactor_factors(unit)
    if sf is not None:
        # S-factor (an SFC reaction) -> cross section.
        E_cm_kev = _cm_kev(data, cols, m_t, m_p)
        to_barn_kev = sf[0] * 1e3           # b*MeV -> b*keV, to match E_cm_kev
        za_t, za_p = data.za(data.target), data.za(data.projectile)
        # No silent default here.  The penetrability factor exp(-2 pi eta) is
        # exponential in Z1 Z2, so an assumed charge is not a small error: the
        # old p + 12C fallback turned the LUNA 3He(a,g)7Be S factor into a
        # cross section 5000 times too large at 93 keV.
        if not (za_t and za_p):
            raise ValueError(
                f"{data.dataset_id}: cannot read the entrance-channel charges "
                f"from {data.reaction!r}, and an S-factor cannot be converted "
                f"without them (sigma = S exp(-2 pi eta)/E, eta ~ Z1 Z2).")
        z2, z1 = za_t[0], za_p[0]
        mu = m_t * m_p / (m_t + m_p)
        cs = sfactor_to_cross_section(cs * to_barn_kev, E_cm_kev, z1, z2, mu)
        stat = sfactor_to_cross_section(stat * to_barn_kev, E_cm_kev, z1, z2, mu)
    elif _is_rutherford_ratio(data, cols, rutherford):
        # Ratio to Rutherford (EXFOR quantity code ...,,RTH): dimensionless,
        # but not a cross section -- passing it through as barns would be
        # wrong by the Coulomb cross section itself.
        za_t, za_p = data.za(data.target), data.za(data.projectile)
        z2 = za_t[0] if za_t else 0
        z1 = za_p[0] if za_p else 0
        if not (z1 and z2):
            raise ValueError(
                f"{data.dataset_id}: ratio-to-Rutherford data needs the "
                f"entrance-channel charges, which could not be read from "
                f"{data.reaction!r}.")
        coulomb = _rutherford_barn_per_sr(
            z1, z2, _cm_kev(data, cols, m_t, m_p) / 1e3,
            _angle_column(data, cols))
        cs = cs * coulomb
        stat = stat * coulomb
    elif "NO-DIM" in unit or "ARB" in unit:
        pass                                 # analyzing power, yield ratio, ...
    else:
        factor = _barn_factor(unit)
        cs = cs * factor
        stat = stat * factor
    return cs * cross_section_scale, stat * cross_section_scale


def _cm_kev(data, cols, m_t, m_p):
    """Center-of-mass energy of each point in keV."""
    for c in data.columns:
        if c.name == "EN-CM":
            E = data.arrays["EN-CM"].astype(float) * _energy_factor(c.unit)
            return E * 1e3
    E_lab = _energy_column(data, cols, m_t, m_p)
    return E_lab * m_t / (m_t + m_p) * 1e3


# -- ENSDF / LiveChart -------------------------------------------------------

@dataclasses.dataclass
class Level:
    """One excited state from :func:`fetch_levels`."""

    z: int
    n: int
    symbol: str
    index: int
    energy: float       # keV
    energy_unc: float   # keV
    jp: str             # "1+", "3/2-", ...
    half_life: str      # "STABLE", "1.23E-3 S", ...
    half_life_sec: float
    width_kev: float = None    # ENSDF total width, keV (NNDC source only; None when not given)
    source: str = "livechart"  # "livechart" or "nndc": which service the record came from

    def __repr__(self):  # pragma: no cover
        return (f"<Level {self.symbol}{self.z + self.n} idx={self.index} "
                f"{self.energy:.3f} keV {self.jp} {self.half_life}>")

    @property
    def energy_mev(self):
        """Level energy in MeV."""
        return self.energy / 1e3

    @property
    def spin(self):
        """Half-integer spin ``J`` as a float (``3/2-`` -> 1.5)."""
        if not self.jp:
            return None
        j = re.match(r"([0-9]+(?:/[0-9]+)?)", self.jp)
        if not j:
            return None
        v = j.group(1)
        if "/" in v:
            num, den = v.split("/")
            return int(num) / int(den)
        return int(v)

    @property
    def parity(self):
        """Parity as +1 or -1, or None when EXFOR gives none."""
        if not self.jp:
            return None
        return +1 if self.jp.endswith("+") else -1


@dataclasses.dataclass
class Gamma:
    """One gamma transition from :func:`fetch_gammas`."""

    start_energy: float    # keV
    start_jp: str
    end_energy: float      # keV
    end_jp: str
    energy: float          # keV
    intensity: float       # relative %
    multipolarity: str


@dataclasses.dataclass
class GroundState:
    """Ground-state properties from :func:`fetch_ground_state`."""

    z: int
    n: int
    symbol: str
    jp: str
    half_life: str
    half_life_sec: float
    spin: float
    parity: int
    mass_excess: float     # keV
    neutron_sep: float     # keV (Sn)
    proton_sep: float      # keV (Sp)
    radius: float          # fm


def _livechart(fields, nuclides):
    if not re.match(r"^[0-9]+[a-z]+$", nuclides):
        raise ValueError(f"nuclides must be like '14n', got {nuclides!r}.")
    text = _http_get(LIVECHART_URL, {"fields": fields, "nuclides": nuclides})
    reader = csv.DictReader(io.StringIO(text))
    # A challenge or error page (HTML, sometimes with HTTP 200) has no CSV header: without this check it
    # reads as a header line and zero rows, i.e. a silent empty level scheme.
    if text.lstrip().startswith("<") or not reader.fieldnames or "z" not in reader.fieldnames:
        raise ValueError(f"LiveChart did not answer with its CSV (starts {text.strip()[:60]!r}).")
    return list(reader)


def fetch_levels(nuclide, source="auto"):
    """Level scheme for ``nuclide`` (e.g. ``"14n"``) from the ENSDF adopted levels.

    Returns a list of :class:`Level` ordered by excitation energy.  Energies
    are in keV (``.energy_mev`` gives MeV); ``.spin`` / ``.parity`` parse the
    ``J^pi`` column.

    ``source`` picks the service: ``"livechart"`` (IAEA LiveChart, CSV),
    ``"nndc"`` (the NNDC ENSDF API at ``ENSDF_API_URL``, JSON) or ``"auto"``
    (the default): LiveChart first, and the NNDC API with a warning when
    LiveChart cannot be reached or does not answer with its CSV -- from some
    networks every IAEA request gets a Cloudflare challenge page (HTTP 403;
    the Notre Dame CRC cluster, 2026-10-15).  Both serve the ENSDF adopted
    levels, and the NNDC records are written in LiveChart's formats, checked
    against it on 14N and 16O (2026-10-15: same levels, energies and
    uncertainties): a multi-valued ``jp`` as ``"2+ 1+"`` / ``"(0+ 2+)"`` (the API
    says ``"2+,1+"`` / ``"(0+),(2+)"``), ``half_life`` as the bare number
    LiveChart gives (``"68"`` for 68 fs, ``"12"`` for a 12 keV width -- the unit
    is not in it; use ``half_life_sec``) or ``"STABLE"``.  A NNDC record also
    carries ``width_kev``, the evaluated total width in keV when ENSDF gives one.
    """
    if source not in ("auto", "livechart", "nndc"):
        raise ValueError(f"source must be 'auto', 'livechart' or 'nndc', got {source!r}.")
    if not re.match(r"^[0-9]+[a-zA-Z]+$", nuclide):
        raise ValueError(f"nuclides must be like '14n', got {nuclide!r}.")
    if source != "nndc":
        try:
            return _levels_livechart(nuclide.lower())
        except (OSError, KeyError, ValueError) as err:   # URLError/HTTPError/timeouts are OSError;
            if source == "livechart":                      # KeyError/ValueError = not LiveChart's CSV
                raise
            warnings.warn(f"LiveChart did not answer for {nuclide!r} ({err.__class__.__name__}: "
                          f"{str(err)[:120]}); using the NNDC ENSDF API ({ENSDF_API_URL}).")
    return _levels_nndc(nuclide)


def _levels_livechart(nuclide):
    out = []
    for row in _livechart("levels", nuclide):
        out.append(Level(
            z=int(row["z"]), n=int(row["n"]), symbol=row["symbol"],
            index=int(row["idx"]),
            energy=float(row["energy"] or 0),
            energy_unc=float(row["unc_e"] or 0),
            jp=row["jp"], half_life=row["half_life"],
            half_life_sec=float(row["half_life_sec"] or 0),
        ))
    return out


_WIDTH_TO_KEV = {"EV": 1e-3, "KEV": 1.0, "MEV": 1e3}


def _levels_nndc(nuclide, page=500):
    """ENSDF adopted levels from the NNDC ENSDF API (``/levels``, paged)."""
    rows, offset = [], 0
    while True:
        d = json.loads(_http_get(f"{ENSDF_API_URL}/levels",
                                 {"nuclide": nuclide, "document-type": "adopted",
                                  "limit": page, "offset": offset}))
        rows += d["results"]
        offset += page
        if not d["meta"].get("next") or offset >= d["meta"].get("total", 0):
            break
    out = []
    for r in rows:
        ds = r.get("datasetSummary") or {}
        e = r.get("energy") or {}
        unc = (e.get("uncertainty") or {}).get("value")
        hl = r.get("halflife") or {}
        w = (r.get("decayWidths") or {}).get("totalWidth") or {}
        if r.get("isStable"):
            half_life = "STABLE"
        elif hl.get("value") is not None:
            half_life = f"{hl['value']:g}"            # LiveChart: the number alone, unit elsewhere
        elif w.get("value") is not None:
            half_life = f"{w['value']:g}"             # a width (keV for broad levels), as LiveChart
        else:
            half_life = ""
        width = None
        if w.get("value") is not None and str(w.get("unit", "")).upper() in _WIDTH_TO_KEV:
            width = float(w["value"]) * _WIDTH_TO_KEV[str(w["unit"]).upper()]
        out.append(Level(
            z=int(ds.get("z", 0)), n=int(ds.get("n", 0)),
            symbol=str(ds.get("elementSymbol", "")),
            index=int(r.get("levelIndex", len(out))),
            energy=float(e.get("value") or 0),
            energy_unc=float(unc or 0),
            jp=((r.get("spinParity") or {}).get("stringRepresentation") or "")
                .replace("),(", " ").replace(",", " "),     # LiveChart's spelling of a J^pi list
            half_life=half_life,
            half_life_sec=float(r.get("halflifeSeconds") or 0),
            width_kev=width, source="nndc",
        ))
    out.sort(key=lambda lv: (lv.energy, lv.index))
    return out


def fetch_gammas(nuclide):
    """Gamma transitions for ``nuclide`` from ENSDF via LiveChart."""
    out = []
    for row in _livechart("gammas", nuclide):
        out.append(Gamma(
            start_energy=float(row["start_level_energy"] or 0),
            start_jp=row["start_level_jp"],
            end_energy=float(row["end_level_energy"] or 0),
            end_jp=row["end_level_jp"],
            energy=float(row["energy"] or 0),
            intensity=float(row["relative_intensity"] or 0),
            multipolarity=row["multipolarity"] or "",
        ))
    return out


def fetch_ground_state(nuclide):
    """Ground-state properties for ``nuclide`` (masses, separations, J^pi)."""
    rows = _livechart("ground_states", nuclide)
    if not rows:
        raise ValueError(f"no ground-state data for {nuclide!r}.")
    r = rows[0]
    jp = r["jp"]
    spin, parity = None, None
    if jp:
        m = re.match(r"([0-9]+(?:/[0-9]+)?)([+-])", jp)
        if m:
            v = m.group(1)
            spin = int(v.split("/")[0]) / int(v.split("/")[1]) if "/" in v \
                else int(v)
            parity = +1 if m.group(2) == "+" else -1
    return GroundState(
        z=int(r["z"]), n=int(r["n"]), symbol=r["symbol"],
        jp=jp, half_life=r["half_life"],
        half_life_sec=float(r["half_life_sec"] or 0),
        spin=spin, parity=parity,
        mass_excess=float(r["massexcess"] or 0),
        neutron_sep=float(r["sn"] or 0),
        proton_sep=float(r["sp"] or 0),
        radius=float(r["radius"] or 0),
    )


# -- references --------------------------------------------------------------

@dataclasses.dataclass
class Reference:
    """The bibliographic record of an EXFOR entry."""

    entry: str            # EXFOR entry ID, e.g. "O2599"
    title: str
    authors: str
    reference: str        # "Jour: Physical Review Letters, Vol.131, p.162701 (2023)"
    x4ref: str            # EXFOR reference code, e.g. "J,PRL,131,162701,2023"
    year: int
    doi: str = ""

    def __repr__(self):  # pragma: no cover
        return (f"<Reference {self.entry}: {self.reference}"
                + (f"  doi:{self.doi}" if self.doi else "") + ">")


def reference(entry_or_dataset_id, source="auto"):
    """Bibliography for an EXFOR **entry** (``"O2599"``) or dataset
    (``"O2599004"`` -> entry ``O2599``).

    Returns a :class:`Reference`.  The ``REFERENCE`` line is parsed from the
    entry's EXFOR BIB section; ``doi`` is left blank (use :func:`resolve_doi`).
    ``source`` as in :func:`fetch_exfor`: the local EXFOR master files hold the
    same BIB text.
    """
    entry = entry_or_dataset_id[:5]
    if source not in ("auto", "web", "local"):
        raise ValueError(f"source must be 'auto', 'web' or 'local', not {source!r}")
    if source == "local":
        text = _local_entry_text(entry)
    else:
        try:
            text = _web_text(EXFOR_URL + "x4get", {"sub": entry})
        except (OSError, _WebUnavailable) as err:
            if source == "web" or exfor_master_path() is None:
                raise
            warnings.warn(f"EXFOR Web-API unavailable for {entry} ({err}); "
                          f"reading the local EXFOR master files.")
            text = _local_entry_text(entry)
    title = _bib_field(text, "TITLE")
    authors = _bib_field(text, "AUTHOR")
    ref_line = _bib_field(text, "REFERENCE")
    x4ref = _bib_field(text, "X4REF") or _parse_x4ref(ref_line)
    year = 0
    for part in (x4ref or "").split(","):
        if part.isdigit() and len(part) == 4:
            year = int(part)
            break
    if year == 0 and ref_line:
        for part in re.split(r"[,() ]+", ref_line):
            if part.isdigit() and len(part) == 4:
                year = int(part)
                break
    return Reference(entry=entry, title=title, authors=authors,
                     reference=ref_line, x4ref=x4ref, year=year)


def _parse_x4ref(ref_line):
    """Recover the ``J,PRL,131,162701,2023`` code from a ``REFERENCE`` line
    like ``(J,PRL,131,162701,2023)`` or a free-text citation."""
    # journal codes may hold "/" and "-" (PR/C, ZP/A, J,NP/A): without them the
    # code of every Phys. Rev. C entry came back empty
    m = re.search(r"\(([A-Z]+,[A-Z0-9+/\-]+,[0-9]+,[0-9A-Z]+(?:,[0-9]{4,6})?)\)",
                  ref_line or "")
    return m.group(1) if m else ""


def _bib_field(text, key):
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if not line.startswith(key):
            continue
        vals = [line[len(key):].strip()]
        j = i + 1
        while j < len(lines) and lines[j][:1] in (" ", "\t", "+") and \
                lines[j].strip() and not lines[j].strip().startswith("DATA"):
            vals.append(lines[j].lstrip("+ \t"))
            j += 1
        return " ".join(v for v in vals if v)
    return ""


def resolve_doi(reference_or_title, year=None, preferred_doi=None):
    """Best-effort DOI lookup for an EXFOR citation.

    Accepts a :class:`Reference` or a plain title string.  When given a
    :class:`Reference`, its structured ``x4ref`` code (``J,PRL,131,162701,2023``)
    is used to construct the DOI for journals with predictable patterns (APS),
    verified against CrossRef; otherwise the title is searched on CrossRef.
    Returns a DOI string or None.
    """
    if preferred_doi:
        return preferred_doi
    if isinstance(reference_or_title, Reference):
        ref = reference_or_title
        x4 = (ref.x4ref or "").split(",")
        if len(x4) >= 4 and x4[0] == "J":
            journal, volume, page = x4[1], x4[2], x4[3]
            doi = aps_doi(journal, volume, page)
            if doi:
                try:                        # confirm the guess resolves
                    _http_get(f"{CROSSREF_URL}/{doi}", timeout=30)
                    return doi
                except OSError:
                    pass
        title = ref.title or ref.reference
        y = ref.year or year
    else:
        title = reference_or_title
        y = year
    if not title:
        return None
    params = {"query.title": title, "rows": 5, "select": "DOI,title,issued"}
    if y:
        params["filter"] = (f"from-pub-date:{y}-01-01,"
                            f"until-pub-date:{y}-12-31")
    body = _http_get(CROSSREF_URL, params)
    items = json.loads(body).get("message", {}).get("items", [])
    if not items:
        return None
    return items[0].get("DOI")


def aps_doi(journal, volume, page):
    """Construct the DOI for an APS article from ``(PRL|PRC|.../vol/page)``."""
    names = {"PRL": "physrevlett", "PRC": "physrevc", "PRA": "physreva",
             "PRD": "physrevd", "PR": "physrev", "P": "physrev"}
    prefix = names.get(journal)
    return f"10.1103/{prefix}.{volume}.{page}" if prefix else None


# -- command line -------------------------------------------------------------

def _main(argv=None):
    """``python -m pyazr.nds`` -- search EXFOR, or fetch one dataset.

    The scripted path is the module API; this is for looking something up
    without opening an editor.
    """
    import argparse
    import sys

    ap = argparse.ArgumentParser(
        prog="python -m pyazr.nds",
        description="Search the IAEA EXFOR database, or convert one dataset "
                    "into an AZURE2 data file.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("search", help="search by target / reaction / quantity")
    s.add_argument("--target", help='EXFOR nucleus, e.g. "C-13"')
    s.add_argument("--reaction", help='projectile,exit -- "p,g", "p,el"')
    s.add_argument("--quantity", help='"SIG", "DA", "POL", ...')
    s.add_argument("--author")
    s.add_argument("--limit", type=int, default=50)

    d = sub.add_parser("download", help="fetch a dataset by its DatasetID")
    d.add_argument("dataset_id")
    d.add_argument("-o", "--output", help="write the AZURE2 data file here "
                                          "(default: stdout)")
    d.add_argument("--source", choices=("auto", "web", "local"), default="auto",
                   help="the Web-API, the local EXFOR master files, or the "
                        "Web-API falling back to them (default)")

    r = sub.add_parser("reference", help="the paper behind an entry, with DOI")
    r.add_argument("entry_or_dataset_id")

    args = ap.parse_args(argv)

    if args.cmd == "search":
        hits = search_exfor(target=args.target, reaction=args.reaction,
                            quantity=args.quantity, author=args.author,
                            limit=args.limit)
        for h in hits:
            print(f"{h.dataset_id}\t{h.npoints:>5} pts\t"
                  f"{h.en_min / 1e6:.4g}-{h.en_max / 1e6:.4g} MeV\t"
                  f"{h.author} {h.year}\t{h.reaction}")
        if not hits:
            print("no datasets matched.", file=sys.stderr)
            return 1
        return 0

    if args.cmd == "download":
        data = fetch_exfor(args.dataset_id, source=args.source)
        rows = _azr_columns(data)
        text = "".join("  ".join(f"{v:.8g}" for v in row) + "\n" for row in rows)
        print(f"# {args.dataset_id}: {len(rows)} points, "
              f"observable {data._guess_observable()}", file=sys.stderr)
        if args.output:
            with open(args.output, "w") as f:
                f.write(text)
        else:
            sys.stdout.write(text)
        return 0

    ref = reference(args.entry_or_dataset_id)
    print(ref)
    doi = resolve_doi(ref)
    print(f"DOI: {doi}" if doi else "DOI: not found")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
