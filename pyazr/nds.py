"""Fetch nuclear data from IAEA EXFOR and NDS (LiveChart/ENSDF) in
AZURE2-friendly form.

The two web services wrapped here:

- **EXFOR** (``https://nds.iaea.org/exfor/``): experimental cross sections,
  differential cross sections, analyzing powers, yields, and more, with the
  bibliographic record for each entry.  ``search_exfor`` / ``fetch_exfor``
  wrap the ``x4list`` / ``x4get`` Web-API.
- **LiveChart/ENSDF** (``https://www-nds.iaea.org/relnsd/v1/data``): evaluated
  level schemes, gamma transitions and ground-state properties.  ``fetch_levels``
  / ``fetch_gammas`` / ``fetch_ground_state`` wrap it.

``ExforData.to_azr`` converts an EXFOR dataset into an AZURE2 data file
(``energy angle crossSection error``, lab frame) and returns the keyword
arguments for :meth:`pyazr.AzrModel.add_data_segment`, so fetched data can be
dropped straight into an evaluation::

    data = fetch_exfor("O2599004")                      # 13C(p,g)14N S-factor
    kw = data.to_azr("run_dir/data", entrance=1, exit=2,
                     observable="total-capture")
    AzrModel.from_file("13N.azr").add_data_segment(**kw).write("13N_new.azr")

References resolve to a DOI through the CrossRef API, so the paper behind any
dataset can be pulled up (``reference`` + ``resolve_doi``).
"""

import csv
import dataclasses
import io
import re
import warnings

import numpy as np
import requests

EXFOR_URL = "https://nds.iaea.org/exfor/"
LIVECHART_URL = "https://www-nds.iaea.org/relnsd/v1/data"
CROSSREF_URL = "https://api.crossref.org/works"
UA = {"User-Agent": "pyazr/nds (research tool)"}

#: EXFOR energy-unit string -> factor to lab MeV (after the frame conversion).
_ENERGY_FACTOR = {"EV": 1e-6, "KEV": 1e-3, "MEV": 1.0, "GEV": 1e3}
#: EXFOR data-unit string -> factor to barns (or barns/sr for differential).
_DATA_FACTOR = {"B": 1.0, "MB": 1e-3, "B/SR": 1.0, "MB/SR": 1e-3,
                "NB/SR": 1e-9, "UB/SR": 1e-12, "NO-DIM": 1.0, "ARB": 1.0}
#: Projectile shorthand used by EXFOR reaction strings -> mass number.
_PROJ_MASS = {"P": 1, "N": 1, "D": 2, "T": 3, "3-HE": 3, "4-HE": 4,
              "ALPHA": 4, "G": 0}


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
        m = re.search(r"\(([^)]+)\)", self.reaction)
        if not m:
            return ""
        return m.group(1).split(",")[0]

    @property
    def target(self):
        m = re.match(r"([0-9]+-[A-Z]+-[0-9]+)", self.reaction)
        return m.group(1) if m else ""

    @property
    def exit(self):
        m = re.search(r"\(([^)]+)\)", self.reaction)
        if not m:
            return ""
        parts = m.group(1).split(",")
        return parts[1] if len(parts) > 1 else ""

    def za(self, name):
        """``(Z, A)`` of a nucleus named like ``6-C-13`` / ``1-H-1``, or None."""
        if not name:
            return None
        m = re.match(r"([0-9]+)-([A-Z]+)-([0-9]+)", name)
        if m:
            return int(m.group(1)), int(m.group(3))
        return None

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
               file_name=None):
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

        Other unit conversions (``MB``, ``NB/SR``, ...) are applied from the
        column unit; ``cross_section_scale`` / ``energy_scale`` / ``angle_scale``
        multiply on top for any manual fix-up.
        """
        os = _import_os()
        data = _azr_columns(self, target_mass=target_mass,
                            projectile_mass=projectile_mass,
                            cross_section_scale=cross_section_scale,
                            energy_scale=energy_scale, angle_scale=angle_scale)
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
    r = requests.get(EXFOR_URL + "x4list", params=params, headers=UA, timeout=60)
    r.raise_for_status()
    try:
        payload = r.json()
    except ValueError:
        raise ValueError(f"EXFOR x4list returned non-JSON:\n{r.text[:300]}")
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


def fetch_exfor(dataset_id, plus=0):
    """Fetch one EXFOR dataset by its ID (``"O2599004"``).

    ``plus=0`` returns the original units (``B*KEV`` S-factors, ``NB/SR``
    differentials, ...); ``plus=1`` the computational units (eV, ``B*EV``,
    ``B/SR``); ``plus=2`` the universal grid (residual/incident/sec energies,
    angle).  Returns an :class:`ExforData`.
    """
    params = {"DatasetID": dataset_id, "op": "csv"}
    if plus:
        params["plus"] = int(plus)
    r = requests.get(EXFOR_URL + "x4get", params=params, headers=UA, timeout=60)
    r.raise_for_status()
    text = r.text
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
                 cross_section_scale=1.0, energy_scale=1.0, angle_scale=1.0):
    """Build the ``(lab E_MeV, angle_deg, cs, err)`` array for a dataset."""
    m_t, m_p = data.masses(target_mass, projectile_mass)
    cols = {c.name: c.unit for c in data.columns}
    energy = _energy_column(data, cols, m_t, m_p) * energy_scale
    angle = _angle_column(data, cols) * angle_scale
    cs, err = _data_column(data, cols, m_t, m_p,
                           cross_section_scale=cross_section_scale)
    return np.column_stack([energy, angle, cs, err])


def _energy_column(data, cols, m_t, m_p):
    name = None
    for c in data.columns:
        if c.name in ("EN-CM", "EN"):
            name = c.name
    if name is None:
        raise ValueError(f"{data.dataset_id}: no EN / EN-CM energy column in "
                         f"{list(cols)}.")
    unit = cols[name]
    E = data.arrays[name].astype(float)
    factor = _ENERGY_FACTOR.get(unit)
    if factor is None:
        warnings.warn(f"{data.dataset_id}: unknown energy unit {unit!r}; "
                      f"assuming MeV.")
        factor = 1.0
    E = E * factor
    if name == "EN-CM":
        E = cm_to_lab(E, m_t, m_p)
    return E


def _angle_column(data, cols):
    for c in data.columns:
        if c.name in ("ANG-CM", "ANG"):
            return data.arrays[c.name].astype(float)
    if "EN" not in cols and "EN-CM" not in cols:
        return np.zeros(len(data.arrays["DATA"]))
    return np.zeros(len(data.arrays["DATA"]))


def _data_column(data, cols, m_t, m_p, cross_section_scale=1.0):
    dname = "DATA-CM" if "DATA-CM" in cols else "DATA"
    if dname not in data.arrays:
        raise ValueError(f"{data.dataset_id}: no DATA / DATA-CM column in "
                         f"{list(cols)}.")
    unit = cols[dname]
    cs = data.arrays[dname].astype(float)
    ename = "ERR-S" if "ERR-S" in data.arrays else \
            ("DATA-ERR" if "DATA-ERR" in data.arrays else None)
    stat = data.arrays[ename].astype(float) if ename else \
        np.zeros_like(cs)
    if unit in ("B*KEV", "B*EV"):
        E_cm_kev = _cm_kev(data, cols, m_t, m_p)
        if unit == "B*EV":              # S in b·eV -> b·keV, E stays in keV
            cs = cs * 1e-3
            stat = stat * 1e-3
        za_t = data.za(data.target)
        za_p = data.za(data.projectile)
        z2 = za_t[0] if za_t else 6
        z1 = za_p[0] if za_p else 1
        mu = m_t * m_p / (m_t + m_p)
        cs = sfactor_to_cross_section(cs, E_cm_kev, z1, z2, mu)
        stat = sfactor_to_cross_section(stat, E_cm_kev, z1, z2, mu)
    elif unit in _DATA_FACTOR:
        cs = cs * _DATA_FACTOR[unit]
        stat = stat * _DATA_FACTOR[unit]
    else:
        warnings.warn(f"{data.dataset_id}: unknown data unit {unit!r}; "
                      f"leaving values unchanged.")
    return cs * cross_section_scale, stat * cross_section_scale


def _cm_kev(data, cols, m_t, m_p):
    """Center-of-mass energy of each point in keV."""
    for c in data.columns:
        if c.name == "EN-CM":
            E = data.arrays["EN-CM"].astype(float)
            unit = c.unit
            if unit in _ENERGY_FACTOR:
                E = E * _ENERGY_FACTOR[unit]
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

    def __repr__(self):  # pragma: no cover
        return (f"<Level {self.symbol}{self.z + self.n} idx={self.index} "
                f"{self.energy:.3f} keV {self.jp} {self.half_life}>")

    @property
    def energy_mev(self):
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
    r = requests.get(LIVECHART_URL,
                     params={"fields": fields, "nuclides": nuclides},
                     headers=UA, timeout=60)
    r.raise_for_status()
    return list(csv.DictReader(io.StringIO(r.text)))


def fetch_levels(nuclide):
    """Level scheme for ``nuclide`` (e.g. ``"14n"``) from ENSDF via LiveChart.

    Returns a list of :class:`Level` ordered by excitation energy.  Energies
    are in keV (``.energy_mev`` gives MeV); ``.spin`` / ``.parity`` parse the
    ``J^pi`` column.
    """
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


def reference(entry_or_dataset_id):
    """Bibliography for an EXFOR **entry** (``"O2599"``) or dataset
    (``"O2599004"`` -> entry ``O2599``).

    Returns a :class:`Reference`.  The ``REFERENCE`` line is parsed from the
    entry's EXFOR BIB section; ``doi`` is left blank (use :func:`resolve_doi`).
    """
    entry = entry_or_dataset_id[:5]
    r = requests.get(EXFOR_URL + "x4get", params={"sub": entry}, headers=UA,
                     timeout=60)
    r.raise_for_status()
    text = r.text
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
    m = re.search(r"\(([A-Z]+,[A-Z0-9+]+,[0-9]+,[0-9]+(?:,[0-9]{4})?)\)",
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
                try:
                    r = requests.get(f"{CROSSREF_URL}/{doi}",
                                     headers=UA, timeout=30)
                    if r.ok:
                        return doi
                except requests.RequestException:
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
    r = requests.get(CROSSREF_URL, params=params, headers=UA, timeout=60)
    r.raise_for_status()
    items = r.json().get("message", {}).get("items", [])
    if not items:
        return None
    return items[0].get("DOI")


def aps_doi(journal, volume, page):
    """Construct the DOI for an APS article from ``(PRL|PRC|.../vol/page)``."""
    names = {"PRL": "physrevlett", "PRC": "physrevc", "PRA": "physreva",
             "PRD": "physrevd", "PR": "physrev", "P": "physrev"}
    prefix = names.get(journal)
    return f"10.1103/{prefix}.{volume}.{page}" if prefix else None
