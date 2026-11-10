"""EXFOR ratios to a Coulomb cross section as format-3 quantities (docs/upgrade/data_kinds.md sections 3.2, 4.3).

The rules ``pyazr.nds.ExforData.format3_quantity`` and ``to_dataset`` follow: which EXFOR quantity codes are a ratio to
Rutherford or Mott, when such an entry is refused, the format-3 label of its reaction and the text of its ``dataset``
record.  Plain functions of an ``ExforData`` (its ``reaction``, ``dataset_id``, ``target``, ``projectile`` and
``za``); no numbers are converted here.  Standard library only (``nds.py`` loads this file by its path).
"""
import re

#: The EXFOR quantity codes (SF5,SF6,SF7,SF8) of a ratio to a Coulomb cross section and the format-3 quantity each
#: maps to: the ``exfor`` column of the quantity table (experiment/src/Quantities.cpp; tests/pyazr/
#: nds_rutherford_test.py checks the two agree).  ``MOT`` is the identical-particle reading of LEXFOR p. S.4;
#: Dict. 236's "relativistic scattering by point nucleus" is refused (data_kinds.md section 3.2, finding 2).
COULOMB_RATIO_CODES = {
    ",DA,,RTH": "ratio_to_rutherford",
    "PAR,DA,,RTH": "ratio_to_rutherford",
    ",DA,,MOT": "ratio_to_mott",
}

#: Format-3 names of the light particles (experiment/src/ReactionParser.cpp), by (Z, A).
LIGHT_LABEL = {(0, 1): "n", (1, 1): "p", (1, 2): "d", (1, 3): "t", (2, 3): "h", (2, 4): "a"}


def quantity_code(reaction):
    """The quantity code (SF5,SF6,SF7,SF8) of an EXFOR reaction string when it is a ratio to a Coulomb cross section,
    else ``None``: ``6-C-12(P,EL)6-C-12,,DA,,RTH`` -> ``,DA,,RTH``."""
    m = re.match(r"^[^(]*\([^)]*\)[^,]*((?:,[^,]*)*)$", reaction.strip())
    if not m:
        return None
    fields = m.group(1).split(",")[1:]  # SF5, SF6, SF7, SF8, SF9
    code = ",".join(fields[:4]).rstrip(",")
    return code if code in COULOMB_RATIO_CODES else None


def is_elastic(data):
    """SF3 ``EL``, or the same ejectile and residual as projectile and target."""
    m = re.search(r"\(([^)]+)\)([^,]*)", data.reaction)
    if not m:
        return False
    parts = m.group(1).split(",")
    if len(parts) < 2:
        return False
    return parts[1].strip().upper() == "EL" or (data.za(parts[1]) == data.za(parts[0]) and
                                                data.za(m.group(2)) == data.za(data.target))


def quantity(data):
    """``ratio_to_rutherford`` for ``,DA,,RTH`` / ``PAR,DA,,RTH`` of an elastic reaction of two different nuclei,
    ``ratio_to_mott`` for ``,DA,,MOT`` of identical particles, ``None`` for any other code.  A ValueError that says
    why for a ratio of an inelastic group (not defined: give dsigma_domega), RTH of identical particles (Rutherford or
    Mott, which the code does not say: decision K2) and MOT of different particles (the point-nucleus cross section)."""
    code = quantity_code(data.reaction)
    if code is None:
        return None
    name = COULOMB_RATIO_CODES[code]
    if not is_elastic(data):
        raise ValueError(f"{data.dataset_id}: {data.reaction!r} is not elastic: a ratio to the Coulomb cross section "
                         f"is not defined for an inelastic group; give dsigma_domega.")
    identical = data.za(data.projectile) == data.za(data.target)
    if name == "ratio_to_rutherford" and identical:
        raise ValueError(f"{data.dataset_id}: {data.reaction!r} is a ratio to Rutherford of identical particles: the "
                         f"reference may be Rutherford or the Mott-symmetrised cross section, which the code does not "
                         f"say (decision K2).  Read the entry and pass quantity='ratio_to_mott' to to_dataset if it is "
                         f"the Mott cross section.")
    if name == "ratio_to_mott" and not identical:
        raise ValueError(f"{data.dataset_id}: {data.reaction!r}: MOT of different particles is the relativistic "
                         f"point-nucleus cross section of Dict. 236, which AZURE2 does not have.")
    return name


def chosen_quantity(data, given=None):
    """The quantity ``to_dataset`` writes: ``given`` (the user's answer to a refusal of :func:`quantity`), else
    :func:`quantity`'s; a ValueError for anything but the two ratio quantities."""
    try:
        found = quantity(data)
    except ValueError:
        if given is None:
            raise
        found = None
    name = given or found
    if name not in ("ratio_to_rutherford", "ratio_to_mott"):
        raise ValueError(f"{data.dataset_id}: {data.reaction!r} is not a ratio to a Coulomb cross section; to_dataset "
                         f"writes ratio_to_rutherford and ratio_to_mott (other data: to_azr).")
    return name


def reaction_label(data):
    """The elastic reaction as a format-3 label: ``6-C-12(P,EL)6-C-12`` -> ``12C(p,p)12C``."""
    za_t, za_p = data.za(data.target), data.za(data.projectile)
    m = re.match(r"[0-9]+-([A-Z]+)-([0-9]+)", data.target)
    if not (za_t and za_p and m):
        raise ValueError(f"{data.dataset_id}: cannot name the nuclei of {data.reaction!r}.")
    target = m.group(2) + m.group(1).capitalize()
    light = LIGHT_LABEL.get(tuple(za_p))
    if light is None:
        pm = re.match(r"[0-9]+-([A-Z]+)-([0-9]+)", data.projectile)
        if not pm:
            raise ValueError(f"{data.dataset_id}: cannot name the projectile of {data.reaction!r}.")
        light = pm.group(2) + pm.group(1).capitalize()
    return f"{target}({light},{light}){target}"


def record(data, name, quantity_name, path, frame, norm):
    """The text of the ``dataset`` record of a ratio data file (``frame``: ``lab`` or ``cm``, of its angle column)."""
    return (f"dataset {name}\n"
            f"    reaction       {reaction_label(data)}\n"
            f"    quantity       {quantity_name}\n"
            f"    unit           1\n"
            f"    file           {path}\n"
            f"    columns        E_lab[MeV] theta_{frame}[deg] value error\n"
            f"    norm           {norm}\n"
            f"    exfor          {data.dataset_id}\n")
