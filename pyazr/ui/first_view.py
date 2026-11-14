"""What the page's first view shows of the lists of pyazr's tables (docs/upgrade/web_ui.md section 11.11): the
columns of the list of data sets (a kind the engine describes: ``pyazr.kinds()["description"]``) and of the list of
particle pairs (``pyazr.kinds()["pairs"]``).

The first view holds the old setup program's columns (its segments table: active, reaction, energy range, angle
range, data type, data file, normalization; its pairs table: the partners, their spins, parities and masses, the
excitation and separation energies, the pair type and the channel radius), so that a list fits a window of 1400
pixels; every other column of the tables is in the expert's depth and in the record's form (the PI's decision of
2026-11-14).  This is the one table
that decides it: a column key here is a key of pyazr's table, or one of the list's own (``name``, ``norm`` the
normalization in use, ``chi2``, ``users``, ``problems`` the validator's messages, ``line``).  Served in
``GET /api/describe`` as ``first_view``.
"""

FIRST_VIEW_COLUMNS = {
    "description": ("name", "tree", "observable", "data_file", "kept_E_lab_MeV", "angle_window_deg", "active", "norm",
                    "problems"),
    "pairs": ("key", "light_name", "light_jpi", "light_mass", "heavy_name", "heavy_jpi", "heavy_mass",
              "heavy_excitation", "separation_energy", "kind_name", "radius"),
}


def first_view_columns():
    """The table as data: ``{list: [column keys]}``."""
    return {name: list(keys) for name, keys in FIRST_VIEW_COLUMNS.items()}
