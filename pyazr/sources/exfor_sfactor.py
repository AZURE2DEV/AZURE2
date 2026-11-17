"""S-factor data of EXFOR as a format-3 data set of quantity ``sfactor`` (``ExforData.to_sfactor_dataset``).

The data are written as S, not converted to a cross section: the engine compares them with E sigma exp(2 pi eta) at
each point's c.m. energy, with the masses of the project's ``<levels>`` (physics/AstrophysicalFactor.h), so that a
fitted energy shift moves the factor and no mass number enters.
"""
import os
import re
import warnings

import numpy as np


def _nds():
    import importlib.util
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nds.py")
    spec = importlib.util.spec_from_file_location("pyazr_nds_for_sfactor", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_LIGHT = {"G": "g", "P": "p", "N": "n", "A": "a", "D": "d", "T": "t"}


def reaction_label(data):
    """The reaction as a format-3 label: ``6-C-12(A,G)8-O-16,,SIG,,SFC`` -> ``12C(a,g)16O``; a ValueError where
    the particles are not single EXFOR letters (give the label by hand)."""
    m = re.match(r"[0-9]+-([A-Z]+)-([0-9]+)\(([A-Z]+),([A-Z]+)\)[0-9]+-([A-Z]+)-([0-9]+)", data.reaction)
    if not m or m.group(3) not in _LIGHT or m.group(4) not in _LIGHT:
        raise ValueError(f"{data.dataset_id}: cannot name the reaction {data.reaction!r} as a format-3 label.")
    return (f"{m.group(2)}{m.group(1).capitalize()}({_LIGHT[m.group(3)]},{_LIGHT[m.group(4)]})"
            f"{m.group(6)}{m.group(5).capitalize()}")


def dataset(data, data_dir, name, norm=None, file_name=None):
    """The file ``E unused S dS`` (E in MeV as the EXFOR energy column gives it: ``EN-CM`` as ``E_cm``, ``EN`` as
    ``E_lab``; S in MeV b) in ``data_dir`` and the text of its ``dataset`` record.  A ValueError for a unit that is not
    barn times energy, a missing energy column, or no ``norm`` and no ERR-SYS in per cent."""
    nds = _nds()
    cols = {c.name: c.unit for c in data.columns}
    unit = (cols.get("DATA-CM") or cols.get("DATA") or "").upper()
    sf = nds._sfactor_factors(unit)
    if sf is None:
        raise ValueError(f"{data.dataset_id}: {unit!r} is not an S-factor unit (barn times energy).")
    energy = "EN-CM" if "EN-CM" in data.arrays else "EN" if "EN" in data.arrays else None
    if energy is None:
        raise ValueError(f"{data.dataset_id}: no EN / EN-CM energy column in {list(cols)}.")
    if norm is None:
        sys_pct = data._sys_error_percent()
        if sys_pct is None:
            raise ValueError(f"{data.dataset_id}: the entry gives no ERR-SYS in per cent: pass norm= (for "
                             f"example '1 fit prior 1 5%', or '1 fixed').")
        norm = f"1 fit prior 1 {sys_pct:g}%"
    S = nds._values(data).astype(float)
    rows = np.column_stack([data.arrays[energy].astype(float) * nds._energy_factor(cols[energy]),
                            np.zeros(len(S)), S * sf[0], nds._error_column(data, cols, S, unit) * sf[0]])
    keep = np.all(np.isfinite(rows), axis=1)
    if not np.all(keep):
        warnings.warn(f"{data.dataset_id}: {int(np.sum(~keep))} of {len(rows)} rows have a blank energy or value "
                      f"and were left out.")
        rows = rows[keep]
    role = "E_cm" if energy == "EN-CM" else "E_lab"
    os.makedirs(data_dir, exist_ok=True)
    path = os.path.join(data_dir, (file_name or nds._safe_name(data.dataset_id, data.author, data.year)) + ".dat")
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# EXFOR {data.dataset_id} {data.reaction}: {role}(MeV) unused S(MeV b) dS(MeV b)\n")
        for row in rows:
            f.write("  ".join(f"{v:.10g}" for v in row) + "\n")
    return (f"dataset {name}\n"
            f"    reaction       {reaction_label(data)}\n"
            f"    quantity       sfactor\n"
            f"    unit           MeV*b\n"
            f"    file           {path}\n"
            f"    columns        {role}[MeV] _ value error\n"
            f"    norm           {norm}\n"
            f"    exfor          {data.dataset_id}\n")
