"""Lightweight client for the IAEA EXFOR Web API.

Python port of gui/src/ExforData.cpp (the "Retrieve Experimental Data from
EXFOR" dialog in AZURESetup), for headless/scripted use -- e.g. from a
notebook or a batch data-gathering script, without launching the Qt GUI.

Mirrors the two calls the GUI makes:

  * ``search_datasets`` -> ``https://nds.iaea.org/exfor/x4list`` (dataset
    search by target/reaction/quantity EXFOR codes, JSON response).
  * ``download_dataset`` -> ``https://nds.iaea.org/exfor/x4get`` (single
    dataset, native-unit computational CSV), parsed into AZURE2's
    (energy[MeV,lab], angle[deg], cross[barn or barn/sr], error) convention
    exactly as ``ExforData::parseCsv`` does, including the astrophysical
    S-factor -> cross-section back-conversion.

Kept dependency-free (urllib + re only) so it works in any environment that
can already run AZURE2/pyazr.
"""

from __future__ import annotations

import csv
import io
import json
import math
import re
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import List, Optional

EXFOR_BASE = "https://nds.iaea.org/exfor"

# Physical constants kept identical to AZURE2's include/Constants.h so the
# S-factor <-> cross-section conversion matches AZURE2's internal one
# (EPoint::CalcEDependentValues / AZUREPlot), and ExforData.cpp's port of it.
# _HBARC matches include/Constants.h, needed for the Rutherford formula below
# (AZURE2 itself has no native "ratio to Rutherford" data mode/formula --
# this conversion is new functionality, not a port of existing C++ logic).
_PI = 3.141592650
_UCONV = 931.4940880
_FSTRUC = 1.00 / 137.0359996790
_HBARC = 197.32696310

_DEFAULT_REL_ERROR = 0.05  # assumed relative uncertainty when none is given


def _rutherford_barn_per_sr(z1: float, z2: float, e_cm_mev: float, angle_deg: float) -> float:
    """Coulomb (Rutherford) differential cross section in the c.m. frame,
    for point charges z1/z2 at c.m. energy e_cm_mev (MeV) and c.m. angle
    angle_deg (degrees). Returns barn/sr.

    dsigma/dOmega = (Z1*Z2*e^2 / (4*Ecm))^2 / sin^4(theta/2), with
    e^2 = fine-structure-constant * hbar*c (MeV*fm), giving fm^2/sr;
    1 fm^2 = 1e-2 barn.
    """
    theta = math.radians(angle_deg)
    sin_half = math.sin(theta / 2.0)
    if e_cm_mev <= 0.0 or sin_half == 0.0:
        return float("nan")
    e2 = _FSTRUC * _HBARC  # MeV*fm
    amplitude_fm = (z1 * z2 * e2) / (4.0 * e_cm_mev)
    dsdo_fm2 = (amplitude_fm ** 2) / (sin_half ** 4)
    return dsdo_fm2 * 1.0e-2


@dataclass
class ExforDataset:
    id: str
    reaction_code: str
    author: str
    reference: str
    npts: int
    en_min: float  # eV
    en_max: float  # eV


@dataclass
class ExforPoint:
    energy: float  # MeV, lab frame
    angle: float = 0.0  # degrees (0 for angle-integrated data)
    cross: float = 0.0  # barn (SIG) or barn/sr (DA)
    error: float = 0.0  # same units as cross


def _fetch(url: str, timeout: float = 30.0) -> bytes:
    request = urllib.request.Request(
        url, headers={"User-Agent": "AZURE2-EXFOR-client/1.0"}
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def search_datasets(
    target: str = "", reaction: str = "", quantity: str = ""
) -> List[ExforDataset]:
    """Search EXFOR for datasets matching EXFOR codes.

    ``target``: e.g. "C-12". ``reaction``: "projectile,emission" EXFOR
    codes, e.g. "P,EL" (elastic scattering) or "P,G" (radiative capture).
    ``quantity``: "SIG" (cross section) or "DA" (angular distribution).
    """
    params = {}
    if target:
        params["Target"] = target
    if reaction:
        params["Reaction"] = reaction
    if quantity:
        params["Quantity"] = quantity
    params["json"] = ""
    url = f"{EXFOR_BASE}/x4list?{urllib.parse.urlencode(params)}"
    body = _fetch(url).decode("utf-8", errors="replace")

    # The live x4list service occasionally emits a syntactically invalid
    # value for zero-point datasets, e.g. `"enMin":,"A1":...` (no value
    # before the comma). This breaks a strict JSON parser -- including,
    # almost certainly, Qt's QJsonDocument::fromJson used by
    # ExforData::parseSearch in the C++ client, meaning any search whose
    # result set happens to include such a dataset would silently fail
    # there too. Repaired here by substituting `null` for the missing value.
    repaired = re.sub(r'":\s*(?=[,}])', '":null', body)

    try:
        doc = json.loads(repaired)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            f"Unexpected response from EXFOR server: {body[:200]!r}"
        ) from exc

    out = []
    for entry in doc.get("x4Datasets", []):
        ds_id = entry.get("id")
        if not ds_id:
            continue
        en_min = entry.get("enMin") or 0.0
        en_max = entry.get("enMax")
        if en_max is None:
            en_max = en_min
        out.append(
            ExforDataset(
                id=ds_id,
                reaction_code=entry.get("RC", ""),
                author=entry.get("A1", ""),
                reference=entry.get("ref", ""),
                npts=entry.get("npts", 0) or 0,
                en_min=en_min,
                en_max=en_max,
            )
        )
    return out


def _extract_unit(header: str) -> str:
    start = header.find("(")
    end = header.find(")", start)
    if start >= 0 and end > start:
        return header[start + 1 : end].strip().upper()
    return ""


def _energy_to_mev(header: str) -> float:
    unit = _extract_unit(header) or header.upper()
    if "GEV" in unit:
        return 1e3
    if "MEV" in unit:
        return 1.0
    if "KEV" in unit:
        return 1.0e-3
    if "MILLI-EV" in unit or "MILLIEV" in unit:
        return 1.0e-9
    if "EV" in unit:
        return 1.0e-6
    return 1.0e-6  # EXFOR computational default is eV


def _cross_to_barn(header: str) -> float:
    unit = _extract_unit(header) or header.upper()
    if unit.startswith("MICRO-B") or unit.startswith("MU-B") or unit.startswith("MUB"):
        return 1.0e-6
    if unit.startswith("MB"):
        return 1.0e-3
    if unit.startswith("KB"):
        return 1.0e3
    if unit.startswith("NB"):
        return 1.0e-9
    if unit.startswith("PB"):
        return 1.0e-12
    if unit.startswith("FB"):
        return 1.0e-15
    if unit.startswith("B"):
        return 1.0
    return 1.0


def _sfactor_barn_mev(header: str):
    """Returns (factor_to_barnMeV, energy_factor_to_MeV) or None if `header`
    is not an astrophysical S-factor column (barn*energy units)."""
    unit = _extract_unit(header) or header.upper()
    if "*" not in unit or "EV" not in unit:
        return None

    barn = 1.0
    if unit.startswith("MICRO-B") or unit.startswith("MU-B") or unit.startswith("MUB"):
        barn = 1.0e-6
    elif unit.startswith("MB"):
        barn = 1.0e-3
    elif unit.startswith("NB"):
        barn = 1.0e-9
    elif unit.startswith("B"):
        barn = 1.0

    energy = 1.0e-6
    if "GEV" in unit:
        energy = 1e3
    elif "MEV" in unit:
        energy = 1.0
    elif "KEV" in unit:
        energy = 1.0e-3
    elif "MILLI-EV" in unit or "MILLIEV" in unit:
        energy = 1.0e-9
    elif "EV" in unit:
        energy = 1.0e-6

    return barn * energy, energy


def _split_csv_line(line: str) -> List[str]:
    return next(csv.reader(io.StringIO(line)))


def _parse_csv(csv_text: str, z1: float, z2: float, m1: float, m2: float):
    lines = [ln for ln in csv_text.splitlines() if ln.strip()]
    if not lines:
        raise RuntimeError("Empty dataset returned by EXFOR server.")
    if "wrong request" in lines[0] or "," not in lines[0]:
        raise RuntimeError(
            f"EXFOR server could not return this dataset: {lines[0].strip()!r}"
        )

    headers = _split_csv_line(lines[0])
    e_col = e_cm_col = data_col = ang_col = -1
    err_abs_col = err_pct_col = -1
    ang_is_cosine = False
    abs_is_total = pct_is_total = False
    is_sfactor = False
    is_rth = False
    err_abs_is_ratio = False
    energy_is_cm = False
    e_factor = e_cm_factor = 1.0e-6
    cross_factor = err_abs_factor = 1.0
    sfactor_factor = 1.0
    sfactor_energy_factor = 0.0

    for i, h in enumerate(headers):
        h = h.strip()
        up = h.upper()
        unit = _extract_unit(up)

        if e_col < 0 and (up.startswith("EN ") or up.startswith("EN(")) and "EV" in unit:
            e_col = i
            e_factor = _energy_to_mev(h)
        elif e_cm_col < 0 and up.startswith("EN-CM") and "EV" in unit:
            e_cm_col = i
            e_cm_factor = _energy_to_mev(h)

        # Data column: "DATA (...)"/"DATA(...)" for angle-integrated data, or
        # "DATA-CM (...)" for a differential cross section already in the
        # centre-of-mass frame (the x4get CSV uses this name whenever a
        # dataset has no companion plain "DATA" column, e.g. most DA
        # datasets) -- but never a "DATA-ERR..." column, which is an error,
        # not a value.
        is_data_header = (
            up.startswith("DATA ") or up.startswith("DATA(") or up.startswith("DATA-CM")
        ) and "ERR" not in up
        if data_col < 0 and is_data_header:
            data_col = i
            sf = _sfactor_barn_mev(h)
            if sf is not None:
                is_sfactor = True
                sfactor_factor, sfactor_energy_factor = sf
            elif unit == "NO-DIM":
                # A dimensionless differential-data column is a ratio to the
                # Rutherford cross section (EXFOR quantity code ...,,RTH),
                # not an absolute cross section already in some unit AZURE2
                # doesn't recognise -- must not be treated as a plain-barn
                # value (cross_factor=1.0 would silently be wrong).
                is_rth = True
            else:
                cross_factor = _cross_to_barn(h)

        if ang_col < 0 and up.startswith("ANG") and "DEG" in unit:
            ang_col = i
            ang_is_cosine = False
        elif ang_col < 0 and up.startswith("COS"):
            ang_col = i
            ang_is_cosine = True

        # Error column: normally "ERR-..." (angle-integrated data), but the
        # x4get CSV names the error of a "DATA-CM" value "DATA-ERR" instead.
        if up.startswith("ERR") or up.startswith("DATA-ERR"):
            is_total = up.startswith("ERR-T")
            is_pct = "PER-CENT" in unit or "PERCENT" in unit or "PER" in unit
            is_abs = bool(unit) and not is_pct and ("B" in unit or unit == "NO-DIM")
            if is_abs and (err_abs_col < 0 or (is_total and not abs_is_total)):
                err_abs_col = i
                err_abs_is_ratio = unit == "NO-DIM"
                err_abs_factor = 1.0 if err_abs_is_ratio else _cross_to_barn(h)
                abs_is_total = is_total
            elif is_pct and (err_pct_col < 0 or (is_total and not pct_is_total)):
                err_pct_col = i
                pct_is_total = is_total

    if e_col < 0 and e_cm_col >= 0:
        e_col = e_cm_col
        e_factor = e_cm_factor
        energy_is_cm = True

    # EXFOR occasionally mislabels the incident-energy unit of S-factor
    # datasets; the S-factor's own energy unit is the self-consistent
    # choice, so trust it over a possibly-wrong energy-column unit.
    if is_sfactor and sfactor_energy_factor > 0.0:
        e_factor = sfactor_energy_factor

    if e_col < 0 or data_col < 0:
        raise RuntimeError(
            f"Could not find energy/cross-section columns. Headers: {headers}"
        )
    if is_rth and (z1 == 0.0 or z2 == 0.0):
        raise RuntimeError(
            "This dataset is a ratio to the Rutherford cross section "
            "(quantity code ...,,RTH) -- pass the entrance-channel charges "
            "z1 (projectile) and z2 (target) to convert it to an absolute "
            "cross section."
        )

    differential = ang_col >= 0
    points: List[ExforPoint] = []
    needed = max(e_col, data_col, err_abs_col, err_pct_col, ang_col)

    for line in lines[1:]:
        fields = _split_csv_line(line)
        if len(fields) <= needed:
            continue
        try:
            e = float(fields[e_col].strip())
            d = float(fields[data_col].strip())
        except ValueError:
            continue  # header echoes / blank rows

        e_raw = e * e_factor
        tot_m = m1 + m2
        if energy_is_cm:
            e_cm = e_raw
            e_lab = e_raw * tot_m / m2 if m2 > 0.0 else e_raw
        else:
            e_lab = e_raw
            e_cm = e_raw * m2 / tot_m if tot_m > 0.0 else e_raw

        # Angle is needed before the cross section itself for RTH data (the
        # Rutherford formula depends on the c.m. angle), so extract it first.
        angle = 0.0
        if ang_col >= 0:
            try:
                a = float(fields[ang_col].strip())
                angle = math.degrees(math.acos(max(-1.0, min(1.0, a)))) if ang_is_cosine else a
            except (ValueError, IndexError):
                pass

        cross_conv = cross_factor
        err_conv = err_abs_factor
        rutherford = None
        if is_sfactor:
            if e_cm <= 0.0:
                continue
            mu = m1 * m2 / tot_m if tot_m > 0.0 else 0.0
            two_pi_eta = (
                2.0 * _PI * math.sqrt(_UCONV / 2.0) * _FSTRUC * z1 * z2 * math.sqrt(mu / e_cm)
            )
            sigma_per_s = math.exp(-two_pi_eta) / e_cm
            cross_conv = sfactor_factor * sigma_per_s
            err_conv = sfactor_factor * sigma_per_s
        elif is_rth:
            rutherford = _rutherford_barn_per_sr(z1, z2, e_cm, angle)
            if math.isnan(rutherford):
                continue
            cross_conv = rutherford
            if err_abs_is_ratio:
                err_conv = rutherford  # error column is itself a ratio to Rutherford

        cross = d * cross_conv

        error = None
        if err_abs_col >= 0:
            try:
                er = float(fields[err_abs_col].strip())
                if er != 0.0:
                    error = er * err_conv
            except (ValueError, IndexError):
                pass
        if error is None and err_pct_col >= 0:
            try:
                pct = float(fields[err_pct_col].strip())
                if pct != 0.0:
                    error = abs(cross) * pct / 100.0
            except (ValueError, IndexError):
                pass
        if error is None:
            error = _DEFAULT_REL_ERROR * abs(cross)

        points.append(ExforPoint(energy=e_lab, angle=angle, cross=cross, error=error))

    if not points:
        raise RuntimeError("No numeric data points could be parsed from this dataset.")
    return points, differential


def download_dataset(
    dataset_id: str, z1: float = 0.0, z2: float = 0.0, m1: float = 0.0, m2: float = 0.0
):
    """Download and parse a single EXFOR dataset by its DatasetID.

    ``z1``/``m1`` are the projectile's charge/mass (amu), ``z2``/``m2`` the
    target's -- needed only to convert astrophysical S-factor data back to
    a cross section. Returns ``(raw_csv, points, differential)``.
    """
    params = {"DatasetID": dataset_id, "op": "csv"}
    url = f"{EXFOR_BASE}/x4get?{urllib.parse.urlencode(params)}"
    raw_csv = _fetch(url).decode("utf-8", errors="replace")
    points, differential = _parse_csv(raw_csv, z1, z2, m1, m2)
    return raw_csv, points, differential


def to_azure_text(points: List[ExforPoint]) -> str:
    """Format parsed points as AZURE2 data-file text (energy angle cross error)."""
    lines = [f"{p.energy:.8g}\t{p.angle:.8g}\t{p.cross:.8g}\t{p.error:.8g}" for p in points]
    return "\n".join(lines) + ("\n" if lines else "")


def _main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_search = sub.add_parser("search", help="search EXFOR by target/reaction/quantity")
    p_search.add_argument("--target", default="", help='e.g. "C-12"')
    p_search.add_argument("--reaction", default="", help='e.g. "P,EL" or "P,G"')
    p_search.add_argument("--quantity", default="", help='"SIG" or "DA"')

    p_dl = sub.add_parser("download", help="download+convert a dataset by DatasetID")
    p_dl.add_argument("dataset_id")
    p_dl.add_argument("--z1", type=float, default=0.0)
    p_dl.add_argument("--z2", type=float, default=0.0)
    p_dl.add_argument("--m1", type=float, default=0.0)
    p_dl.add_argument("--m2", type=float, default=0.0)
    p_dl.add_argument("-o", "--output", default=None, help="write AZURE2-format text here")

    args = parser.parse_args()
    if args.cmd == "search":
        for ds in search_datasets(args.target, args.reaction, args.quantity):
            print(
                f"{ds.id}\t{ds.npts} pts\t{ds.en_min/1e6:.4g}-{ds.en_max/1e6:.4g} MeV"
                f"\t{ds.author}\t{ds.reaction_code}\t{ds.reference}"
            )
    elif args.cmd == "download":
        _, points, differential = download_dataset(
            args.dataset_id, args.z1, args.z2, args.m1, args.m2
        )
        text = to_azure_text(points)
        kind = "differential (b/sr)" if differential else "angle-integrated (barn)"
        print(f"# {args.dataset_id}: {len(points)} point(s), {kind}", file=__import__("sys").stderr)
        if args.output:
            with open(args.output, "w") as f:
                f.write(text)
        else:
            print(text, end="")


if __name__ == "__main__":
    _main()
