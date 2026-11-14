#!/usr/bin/env python3
"""First checks of an AZURE2 project before any fit: format, masses, thresholds, radii, data sets, options.

Reads a format-3 project (``.azr``) and its data files (relative to the project's directory, as a run reads them)
and prints:

1. the engine's validation (``pyazr.validate``): every error and warning with its code;
2. every particle pair against AME2020 (``pyazr.sources.masses.mass_check``): masses in u, S in MeV, and the difference
   of S from the table in keV; a mass that is a mass number (an integer, other than 12C) is flagged;
3. the channel radius of each pair beside 1.4 (A1^1/3 + A2^1/3) fm, a usual starting value (information only);
4. every data set as the engine reads it (``pyazr.describe(path, as_json=True)``): quantity, reaction, rows of
   the file and rows kept inside the windows; a data set that keeps no row is flagged;
5. the options of the calculation with their source (default, project, keyword), and ``includeAllJPi`` when the
   project has elastic-scattering data without it.

Writes nothing. Exit status: 0 when nothing is flagged, 1 when something is (each flag is a line ``FLAG ...``),
2 when the file cannot be read or is not format 3 (the message names the conversion: ``AZURE2 --convert OLD NEW``).
"""
import argparse
import sys
from pathlib import Path


def audit(path, tol_kev=1.0):
    """(lines, flags) of the audit of the project ``path``."""
    import pyazr
    from pyazr.sources import masses
    lines, flags = [], []
    diagnostics = pyazr.validate(path)
    lines.append(f"validation: {sum(d['severity'] == 'error' for d in diagnostics)} error(s), "
                 f"{sum(d['severity'] == 'warning' for d in diagnostics)} warning(s)")
    for d in diagnostics:
        text = f"  {d['severity']} {d['code']} line {d.get('line')}: {d['message']}"
        lines.append(text)
        if d["severity"] == "error":
            flags.append(f"validation {d['code']}: {d['message']}")
    project = pyazr.Project.open(path)
    lines.append("pairs against AME2020:")
    for row in masses.mass_check(project):
        lines.append("  " + row.describe(tol_kev=tol_kev))
        if row.dS_kev is not None and abs(row.dS_kev) > tol_kev:
            flags.append(f"pair {row.pair}: S differs from the masses' table by {row.dS_kev:+.2f} keV")
        for name, m, ident in (("M1", row.M1, row.light), ("M2", row.M2, row.heavy)):
            if ident is not None and ident.A >= 1 and m == round(m) and (ident.Z, ident.A) != (6, 12):
                flags.append(f"pair {row.pair}: {name} = {m} u is a mass number, not an atomic mass")
    lines.append("channel radii (fm), against 1.4 (A1^1/3 + A2^1/3):")
    for q in project.pairs():
        if q.kind == 0:
            usual = 1.4 * (q.light[1] ** (1 / 3) + q.heavy[1] ** (1 / 3))
            lines.append(f"  pair {q.key} {q.label}: a_c = {q.radius:g}, usual start {usual:.2f}")
    described = pyazr.describe(path, as_json=True)
    lines.append("data sets (rows kept / rows in the file):")
    elastic = False
    for d in described.get("datasets", []):
        kept, read = d.get("rows_kept", 0), d.get("rows_read", 0)
        lines.append(f"  {d['name']}: {d['quantity']} {d['tree']} {kept}/{read}" + ("" if d["active"] else " inactive"))
        if d["active"] and kept == 0:
            flags.append(f"data set {d['name']}: no row of {d['data_file']} inside its windows")
        if d["active"] and any(t["entrance_pair"] == t["exit_pair"] for t in d.get("terms", [])):
            elastic = True
    lines.append("options of the calculation:")
    for name, o in described.get("options", {}).items():
        lines.append(f"  {name} = {o['value']} ({o['source']})")
    config = project.config()
    if elastic and str(config.get("includeAllJPi", "0")) not in ("1", "true", "True"):
        lines.append("note: elastic data without includeAllJPi=1: a J^pi without an active level is left out of the "
                     "calculation (its hard-sphere scattering too); see references/azr-levels.md")
    return lines, flags


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--tol-kev", type=float, default=1.0, help="tolerance on S against the mass table (keV)")
    a = ap.parse_args(argv)
    try:
        lines, flags = audit(a.azr, a.tol_kev)
    except (OSError, ValueError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    print("\n".join(lines))
    for f in flags:
        print("FLAG", f)
    print(f"{len(flags)} flag(s)")
    return 1 if flags else 0


if __name__ == "__main__":
    sys.exit(main())
