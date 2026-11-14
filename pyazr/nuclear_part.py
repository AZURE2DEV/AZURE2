"""The nuclear part of a format-3 project in :class:`pyazr.Project`: levels, pairs, ``<config>`` and ``<thm>``.

docs/upgrade/pyazr_on_records.md, "Edits of the nuclear part".  The engine owns each grammar: it reads the
``<levels>`` rows as levels and pairs (``_azure2.level_table``), the ``<config>`` entries (``config_entries``) and the
``<thm>`` lines (``thm_entries``), and makes every edited or new line (``set_level_edits``, ``config_edits``,
``thm_option_edits``, ...).  This module holds no column, position or grammar: it names things (a level by its label
or by J^pi and E_x, a channel by pair key and (l, s), an entry by its name), applies the replacements the engine
returns to the project's document and lets :meth:`Project._apply` validate the project as it would stand.
"""
from dataclasses import dataclass
from typing import List, Optional, Tuple

from .described import _engine as _engine_module

__all__ = ["ProjectLevel", "ProjectChannel", "ProjectPair", "NuclearPart"]


def _text(b):
    return b.decode("utf-8", "surrogateescape") if isinstance(b, bytes) else b


@dataclass(frozen=True)
class ProjectChannel:
    """One row of a level: the channel (pair key, l, s) and its width column.

    ``width`` is the row's value in the row's convention -- a partial width in eV for an open channel or a photon,
    an ANC in fm^-1/2 for a closed channel, a reduced-width amplitude in MeV^1/2 when ``amplitude`` (gammaIsRWA) --
    and ``width_text`` the token as written.  ``channel`` is the engine's index of the channel in its J^pi group (0
    for an inactive row), ``line`` the row's line in the file (1-based)."""
    pair: int
    l: int
    s: float
    width: float
    width_text: str
    fixed: bool
    amplitude: bool
    channel: int
    line: int


@dataclass(frozen=True)
class ProjectLevel:
    """One level of ``<levels>`` (its rows), as :meth:`Project.levels` reads it.

    ``label`` names it in every edit: ``"<J^pi> <E_x as written>"`` (``"3/2- 3.5032"``; ``" #2"`` appended to a later
    level with the same label).  ``energy`` is E_x of the compound nucleus in MeV; ``group`` and ``level`` the
    engine's indices of the fit parameters (:attr:`pyazr.parameters.Parameter.jgroup`, ``.level``), 0 for an inactive
    level, which the engine does not read; ``nucleus`` the compound nucleus (Z, A)."""
    label: str
    jpi: str
    J: float
    parity: int
    energy: float
    energy_text: str
    fixed: bool
    active: bool
    ids: Tuple[int, ...]
    group: int
    level: int
    nucleus: Tuple[int, int]
    channels: Tuple[ProjectChannel, ...]

    def channel(self, pair, l, s):
        """The row of channel (pair key, l, s); KeyError when the level has none."""
        for c in self.channels:
            if c.pair == pair and c.l == l and round(2 * c.s) == round(2 * s):
                return c
        raise KeyError(f"level {self.label!r} has no row for pair {pair}, l = {l}, s = {s}")


@dataclass(frozen=True)
class ProjectPair:
    """A particle pair of ``<levels>`` (the columns of its first row): ``kind`` 0 particles, 10 photon, 20 beta;
    ``light`` and ``heavy`` (Z, A); masses in u, energies in MeV, ``radius`` the channel radius in fm;
    ``binding_energy`` the THM binding energy B (MeV) of the composite projectile whose transferred particle is the
    light partner (0: none)."""
    key: int
    label: str
    kind: int
    light: Tuple[int, int]
    heavy: Tuple[int, int]
    light_mass: float
    heavy_mass: float
    light_excitation: float
    heavy_excitation: float
    separation_energy: float
    threshold: float
    radius: float
    two_j_light: int
    two_j_heavy: int
    heavy_parity: int
    light_parity: int = 1
    binding_energy: float = 0.0
    light_name: str = ""
    heavy_name: str = ""
    light_jpi: str = ""
    heavy_jpi: str = ""
    kind_name: str = ""

    BUILT_BY = "Project.pairs()"
    # The words of each field (pyazr.kinds()["pairs"]): ``title`` (the old setup program's column head where it had
    # the column), ``unit``, ``doc``; ``list``: a column of the list of pairs; ``rare``: a column shown only where a
    # pair sets it (a value other than 0).
    FIELDS = {
        "key": {"title": "Key", "unit": "", "doc": "the pair's key in <levels> (data sets and levels name it)"},
        "label": {"title": "Pair", "unit": "", "doc": "light + heavy partner, as the engine names the pair"},
        "light": {"title": "Light (Z, A)", "unit": "", "doc": "the light partner's charge and mass numbers"},
        "light_name": {"title": "Light particle", "unit": "", "doc": "the light partner", "list": True},
        "light_jpi": {"title": "Light Jπ", "unit": "", "doc": "the light partner's spin and parity",
                      "list": True},
        "light_mass": {"title": "Light M", "unit": "u", "doc": "the light partner's mass", "list": True},
        "heavy": {"title": "Heavy (Z, A)", "unit": "", "doc": "the heavy partner's charge and mass numbers"},
        "heavy_name": {"title": "Heavy particle", "unit": "", "doc": "the heavy partner", "list": True},
        "heavy_jpi": {"title": "Heavy Jπ", "unit": "", "doc": "the heavy partner's spin and parity",
                      "list": True},
        "heavy_mass": {"title": "Heavy M", "unit": "u", "doc": "the heavy partner's mass", "list": True},
        "light_excitation": {"title": "Light excitation energy", "unit": "MeV",
                             "doc": "the light partner's excitation energy", "list": True, "rare": True},
        "heavy_excitation": {"title": "Excitation energy", "unit": "MeV",
                             "doc": "the heavy partner's excitation energy (the final level of a photon pair)",
                             "list": True},
        "separation_energy": {"title": "Separation energy", "unit": "MeV",
                              "doc": "S: the pair opens at E_x = S + the partners' excitation energies",
                              "list": True},
        "kind": {"title": "Pair type number", "unit": "",
                 "doc": "0 particles, 10 photon, 20 beta (the <levels> column)"},
        "kind_name": {"title": "Pair type", "unit": "", "doc": "particle, photon or beta", "list": True},
        "threshold": {"title": "Threshold", "unit": "MeV", "doc": "the compound nucleus's E_x where the pair opens"},
        "radius": {"title": "Channel radius", "unit": "fm", "doc": "the matching radius of the R-matrix",
                   "list": True},
        "two_j_light": {"title": "Light 2J", "unit": "", "doc": "twice the light partner's spin"},
        "two_j_heavy": {"title": "Heavy 2J", "unit": "", "doc": "twice the heavy partner's spin"},
        "heavy_parity": {"title": "Heavy parity", "unit": "", "doc": "+1 or -1"},
        "light_parity": {"title": "Light parity", "unit": "", "doc": "+1 or -1"},
        "binding_energy": {"title": "THM binding energy", "unit": "MeV",
                           "doc": "B of the Trojan horse whose transferred particle is the light partner (0: none)",
                           "list": True, "rare": True},
    }
    # The fields of a new pair's partner (add_pair's dicts ``light`` and ``heavy``), with the same words; ``values``
    # where the field is a choice.
    PARTNER_FIELDS = {
        "J": {"title": "Spin J", "unit": "", "doc": "the partner's spin, a multiple of 1/2 (1/2, 0.5, 1, ...)"},
        "parity": {"title": "Parity", "unit": "", "doc": "the partner's parity", "values": [1, -1],
                   "labels": ["+", "-"]},
        "Z": {"title": "Z", "unit": "", "doc": "the partner's charge number"},
        "mass": {"title": "M", "unit": "u", "doc": "the partner's mass"},
        "g": {"title": "g-factor", "unit": "", "doc": "the partner's g-factor (magnetic moment / (J nuclear "
              "magnetons)); 0 when not needed", "default": 0},
    }
    # add_pair's other arguments that are not a field of a pair, with the same words.
    ARGUMENTS = {
        "external_capture": {"title": "EC multipolarities", "unit": "",
                             "doc": "the external-capture multipole mask of a photon pair (0: none)"},
    }
    # add_pair's kinds, in the order the <levels> numbers them.
    KINDS = {"particle": 0, "photon": 10, "beta": 20}


_KIND_NAMES = {number: name for name, number in ProjectPair.KINDS.items()}


def _partner_names(label):
    """The light and the heavy partner as the engine's pair label names them (``"a+12C[4.43891]"``: the light
    partner, then "+", then the heavy one; no partner's name holds a "+")."""
    light, _, heavy = label.partition("+")
    return light, heavy


def _jpi(two_j, parity):
    """J^pi as a level writes it: ``"1/2+"``, ``"0+"``."""
    return (f"{two_j // 2}" if two_j % 2 == 0 else f"{two_j}/2") + ("+" if parity > 0 else "-")


def _level_name(level):
    """A level as a message names it: its label for a :class:`ProjectLevel` (not its repr), else as given."""
    return level.label if isinstance(level, ProjectLevel) else level


def _level_address(level):
    """What the engine takes as a level's name: a label (str), a :class:`ProjectLevel` (its label), (J^pi, E_x) as
    ``("3/2-", 3.5)`` (the label form "3/2- 3.5": the level of that J^pi within 1 keV), or the engine's indices
    ``(group, level)`` of two integers."""
    if isinstance(level, ProjectLevel):
        return level.label
    if isinstance(level, str):
        return level
    if isinstance(level, tuple) and len(level) == 2:
        if isinstance(level[0], str):
            return f"{level[0]} {level[1]!r}"
        return (int(level[0]), int(level[1]))
    raise TypeError(f"a level is named by its label, a ProjectLevel, (J^pi, E_x) or the engine's (group, level); "
                    f"not {level!r}")


def _two_s(s):
    two = round(2 * float(s))
    if abs(2 * float(s) - two) > 1e-9:
        raise ValueError(f"a channel spin is a multiple of 1/2, not {s}")
    return int(two)


def _flag(name, value):
    """A flag argument as a bool: True, False (or 1, 0); TypeError otherwise."""
    if isinstance(value, bool) or value in (0, 1):
        return bool(value)
    raise TypeError(f"{name} is True or False, not {value!r}")


def _partner(which, d):
    """(2j, parity, Z, mass, g) of a pair partner given as ``{J, parity, Z, mass, g}``."""
    if not isinstance(d, dict) or not {"J", "parity", "Z", "mass"} <= set(d):
        raise ValueError(f"{which} is a dict {{J, parity, Z, mass (u), g (default 0)}}, not {d!r}")
    unknown = set(d) - set(ProjectPair.PARTNER_FIELDS)
    if unknown:
        raise ValueError(f"{which}: unknown keys {sorted(unknown)} (the keys: J, parity, Z, mass, g)")
    return (_two_s(d["J"]), int(d["parity"]), int(d["Z"]), float(d["mass"]), float(d.get("g", 0.0)))


def _replace(doc, replacements):
    for first, last, lines in replacements:
        doc.replace_lines(first, last, list(lines))


def _engine_call(what, call, *args, **kwargs):
    """The engine's replacements, its refusal (ValueError) as FormatError headed by ``what`` (the object edited).
    The engine's refusal of an edit carries its three parts as attributes (``what``, ``allowed``, ``fix``); they
    become the FormatError's, ``what`` headed by the object as the message is."""
    from .project import FormatError
    try:
        return call(*args, **kwargs)
    except ValueError as e:
        part = getattr(e, "what", None)
        raise FormatError(f"{what}: {e}", what=f"{what}: {part}" if part else None,
                          allowed=getattr(e, "allowed", None), fix=getattr(e, "fix", None)) from None


class NuclearPart:
    """The calls of :class:`pyazr.Project` on the nuclear part (mixed into Project).  Every edit is validated by
    the engine on the project as it would stand and returns the new warnings of the validator (a list of the
    engine's diagnostics, empty when there are none); an error raises :class:`pyazr.FormatError` and leaves the
    project unchanged."""

    # -- levels and pairs ---------------------------------------------------------------------------------------

    def _level_table(self, doc=None):
        return _engine_call("<levels>", _engine_module().level_table, doc if doc is not None else self._doc)

    def levels(self) -> List[ProjectLevel]:
        """Every level of ``<levels>`` in the order of its first row (:class:`ProjectLevel`, with its channels)."""
        out = []
        for d in self._level_table()["levels"]:
            channels = tuple(ProjectChannel(c["pair"], c["l"], c["two_s"] / 2, c["width"], c["width_text"], c["fixed"],
                                     c["amplitude"], c["channel"], c["line"]) for c in d["channels"])
            out.append(ProjectLevel(d["label"], d["jpi"], d["two_j"] / 2, d["parity"], d["energy"], d["energy_text"],
                             d["fixed"], d["active"], tuple(d["ids"]), d["group"], d["level"], tuple(d["nucleus"]),
                             channels))
        return out

    def level(self, level) -> ProjectLevel:
        """The :class:`ProjectLevel` named by ``level`` (a label, (J^pi, E_x), a ProjectLevel, or the engine's (group, level));
        FormatError when it names none or several."""
        k = _engine_call(f"level {_level_name(level)!r}", _engine_module().find_level, self._doc, _level_address(level))
        return self.levels()[k]

    def pairs(self) -> List[ProjectPair]:
        """The particle pairs of ``<levels>`` (:class:`ProjectPair`), in the order of their first rows."""
        return [ProjectPair(d["key"], d["label"], d["kind"], tuple(d["light"]), tuple(d["heavy"]), d["light_mass"],
                     d["heavy_mass"], d["light_excitation"], d["heavy_excitation"], d["separation_energy"],
                     d["threshold"], d["radius"], d["two_j_light"], d["two_j_heavy"], d["heavy_parity"],
                     d["light_parity"], d["binding_energy"], *_partner_names(d["label"]),
                     _jpi(d["two_j_light"], d["light_parity"]), _jpi(d["two_j_heavy"], d["heavy_parity"]),
                     _KIND_NAMES.get(d["kind"], str(d["kind"])))
                for d in self._level_table()["pairs"]]

    def reaction(self, entrance, exit):
        """The reaction label from pair ``entrance`` to pair ``exit`` (keys of ``<levels>``; ``exit=-1``: the
        total capture of the entrance pair), as a data set or extrapolation names it: ``"12C(p,g)13N"``."""
        return _engine_call("reaction", _engine_module().reaction_label, self._doc, int(entrance), int(exit))

    def _edit_levels(self, what, call, *args, **kwargs):
        def edit(doc):
            _replace(doc, _engine_call(what, call, doc, *args, **kwargs))
        return self._apply(what, edit)

    def set_level(self, level, *, energy=None, fixed=None, active=None):
        """Set a level's E_x (MeV), whether it is fixed in a fit, whether it is active (an inactive level is not
        read by the engine), on every row of the level; ``None`` leaves a property as it is."""
        return self._edit_levels(f"level {_level_name(level)!r}", _engine_module().set_level_edits, _level_address(level),
                                 energy=None if energy is None else float(energy), fixed=fixed, active=active)

    def set_channel(self, level, pair, l, s, *, width=None, fixed=None, width_is_amplitude=None):
        """Set the width column of channel (pair key, l, s) of a level -- in the row's convention: eV for an open
        channel or a photon, the ANC in fm^-1/2 for a closed channel, MeV^1/2 for a reduced-width amplitude row
        (:attr:`ProjectChannel.amplitude`) -- whether it is fixed in a fit, or ``width_is_amplitude``: whether the
        width column is a reduced-width amplitude (MeV^1/2, column 33 ``gammaIsRWA``) instead of a partial width or
        ANC, a particle channel only (a level known from THM data alone, say).  The value is not converted when the
        flag changes: give ``width`` in the new convention in the same call."""
        return self._edit_levels(f"level {_level_name(level)!r}, channel ({pair}, {l}, {s})", _engine_module().set_channel_edits,
                                 _level_address(level), int(pair), int(l), _two_s(s),
                                 width=None if width is None else float(width), fixed=fixed,
                                 width_is_amplitude=None if width_is_amplitude is None else _flag(
                                     "width_is_amplitude", width_is_amplitude))

    def set_pair(self, pair, *, radius=None, light_mass=None, heavy_mass=None, separation_energy=None,
                 light_excitation=None, heavy_excitation=None, binding_energy=None):
        """Set a pair's channel radius (fm), masses (u), separation energy (MeV), the excitation energies of its
        light and heavy partner (MeV) or its THM binding energy (MeV) on every row of the pair (the separation energy
        also in the entrance-pair copy of the rows whose entrance it is).  The light partner's excitation is the
        optional last column of a row: a row that does not write it gets it only for a value other than 0.

        ``binding_energy`` is B of the Trojan horse a = x + s whose transferred particle x is the pair's light
        partner (column 32): the entrance vertex of a THM experiment whose entrance this pair is reads it.  A particle
        pair only; a number >= 0 (0: none); a row without the THM columns gets them only for a value other than 0.

        The channel radius is the matching surface of the R-matrix: a new radius changes every penetrability, shift
        function and external-capture integral, so the level parameters no longer describe a fit -- refit before
        reading anything off the result."""
        floats = {k: float(v) for k, v in
                  dict(radius=radius, light_mass=light_mass, heavy_mass=heavy_mass, separation_energy=separation_energy,
                       light_excitation=light_excitation, heavy_excitation=heavy_excitation,
                       binding_energy=binding_energy).items() if v is not None}
        return self._edit_levels(f"pair {pair}", _engine_module().set_pair_edits, int(pair), **floats)

    def add_level(self, J, parity, energy, channels, *, fixed=True):
        """Add a level of spin ``J``, parity ``parity`` (+1 or -1) at E_x ``energy`` (MeV) with ``channels``, a list
        of dicts ``{pair, l, s, width (row convention, default 0), fixed (default False)}`` naming existing pairs,
        or ``"allowed"``: the channels of the J^pi group when the project has one, else every channel
        :meth:`allowed_channels` gives with its default limits, each with width 0, fixed.
        The levels of a J^pi share one channel set: a J^pi the project has gets the rows of its group, the channels
        not given with width 0, fixed (a channel outside the set is refused); a new J^pi gets one row per channel,
        the pair's columns copied from a row of the pair.  The level gets the next free level ID and goes at the end
        of ``<levels>``; ``fixed`` fixes its energy.  A pair row (:meth:`add_pair`) of a pair the level names is
        removed: the level now holds the pair."""
        if isinstance(channels, str):
            if channels != "allowed":
                raise ValueError(f"channels is a list of dicts or 'allowed', not {channels!r}")
            two_j, par = _two_s(J), int(parity)
            if any(lv.active and round(2 * lv.J) == two_j and lv.parity == par for lv in self.levels()):
                channels = []  # the group's channels, width 0, fixed
            else:
                channels = [dict(c, width=0.0, fixed=True) for c in self.allowed_channels(J, parity)]
        rows = [(int(c["pair"]), int(c["l"]), _two_s(c["s"]), float(c.get("width", 0.0)), bool(c.get("fixed", False)))
                for c in channels]
        return self._edit_levels(f"new level {J}{'+' if parity > 0 else '-'} at {energy} MeV",
                                 _engine_module().add_level_edits, _two_s(J), int(parity), float(energy),
                                 fixed=bool(fixed), channels=rows)

    def allowed_channels(self, J, parity, *, l_max=2, multipolarity_max=2, multipolarities_max=2, pairs=None):
        """The channels a level of spin ``J`` and parity ``parity`` (+1 or -1) has with the project's pairs (or the
        pair keys ``pairs``), as dicts ``{pair, l, s}`` in pair order, then s and l ascending.  The rule is the Qt
        GUI's (``LevelsTab::calculateChannels``), the engine's ``experiment::AllowedChannels``:

        - a particle pair (spins j1, j2, parities pi1, pi2): every channel spin s = |j1 - j2| ... j1 + j2 and every
          whole l = |s - J| ... s + J with pi1 pi2 (-1)^l = parity and l <= ``l_max``;
        - a photon pair (the final level J_f): s = J_f and the multipolarities L = 1 ... ``multipolarity_max`` with
          |L - J_f| <= J <= L + J_f, the lowest ``multipolarities_max`` of them (E or M follows from the parities);
        - a beta pair (parent j2, pi2), parity = pi2: J = j2 gives s = 0 with l = 1 (Gamow-Teller) and l = 0
          (Fermi), |J - j2| = 1 gives s = 1, l = 1 (Gamow-Teller).

        The pairs of every compound nucleus are taken; a level has the channels of one nucleus (``pairs`` chooses
        them).  :class:`pyazr.FormatError` for a J or parity out of range, a negative limit or an unknown key."""
        found = _engine_call("allowed channels", _engine_module().allowed_channels, self._doc, _two_s(J), int(parity),
                             l_max=int(l_max), multipolarity_max=int(multipolarity_max),
                             multipolarities_max=int(multipolarities_max),
                             pairs=[int(k) for k in (pairs or ())])
        return [{"pair": k, "l": l, "s": two_s / 2} for k, l, two_s in found]

    def add_pair(self, *, kind="particle", light=None, heavy, separation_energy, radius=None, light_excitation=0.0,
                 heavy_excitation=0.0, external_capture=0, binding_energy=0.0, key=None):
        """Add a pair; returns its key (``key``, or the next free key when None).

        ``kind``: ``"particle"``, ``"photon"`` (capture to a level of the compound nucleus: ``heavy`` is that
        nucleus in the final level, ``heavy_excitation`` its E_x) or ``"beta"``.  ``light`` and ``heavy``: dicts
        ``{J, parity, Z, mass (u), g}`` (``g`` the g-factor, default 0; a photon pair takes no ``light``).
        ``separation_energy`` S (MeV), ``radius`` the channel radius (fm, a particle pair), ``light_excitation`` and
        ``heavy_excitation`` (MeV), ``external_capture`` the external-capture multipole mask, ``binding_energy`` the
        THM binding energy (MeV, see :meth:`set_pair`).  The pair opens at E_x = S + e1 + e2 of the compound
        nucleus.

        ``<levels>`` holds a pair only in its rows, so a pair no level names yet is kept by a row of its own (a
        "pair row": inactive, level ID 0), which :meth:`levels` does not list and :meth:`add_level` removes once a
        level names the pair.  Every number is checked by the engine (:class:`pyazr.FormatError` with the reason);
        a key in use is refused."""
        kinds = ProjectPair.KINDS
        if kind not in kinds:
            raise ValueError(f"kind is one of {', '.join(kinds)}, not {kind!r}")
        if kind == "photon" and light is not None:
            raise ValueError("a photon pair has no light partner to give: the photon is implied (light=None)")
        if kind != "photon" and light is None:
            raise ValueError(f"a {kind} pair needs light={{J, parity, Z, mass}}")
        lt = _partner("light", light or {"J": 1, "parity": 1, "Z": 0, "mass": 0.0})
        hv = _partner("heavy", heavy)
        replacements_and_key = []

        def edit(doc):
            replacements, used = _engine_call(
                "new pair", _engine_module().add_pair_edits, doc, key=int(key or 0), kind=kinds[kind],
                two_j_light=lt[0], parity_light=lt[1], z_light=lt[2], mass_light=lt[3], g_light=lt[4],
                two_j_heavy=hv[0], parity_heavy=hv[1], z_heavy=hv[2], mass_heavy=hv[3], g_heavy=hv[4],
                light_excitation=float(light_excitation), heavy_excitation=float(heavy_excitation),
                separation_energy=float(separation_energy), radius=float(radius or 0.0),
                external_capture=int(external_capture), binding_energy=float(binding_energy))
            _replace(doc, replacements)
            replacements_and_key.append(used)
        self._apply("new pair", edit)
        return replacements_and_key[-1]

    def remove_pair(self, pair, *, remove_channels=False):
        """Remove pair ``pair`` (its key).  Refused (:class:`pyazr.FormatError`) while a level has a channel of it,
        naming the levels and channels -- or, with ``remove_channels=True``, those channel rows go with the pair
        (refused when a level would be left without a channel: remove that level first).  The default refuses
        because removing channels removes fit parameters and changes every level that has them; asking for it is
        explicit.  A data set or record that names the pair (its reaction, a potential) is refused by the validator,
        and the project is unchanged."""
        return self._edit_levels(f"remove pair {pair}", _engine_module().remove_pair_edits, int(pair),
                                 remove_channels=_flag("remove_channels", remove_channels))

    def remove_level(self, level):
        """Remove every row of a level."""
        return self._edit_levels(f"remove level {_level_name(level)!r}", _engine_module().remove_level_edits, _level_address(level))

    def _write_level_values(self, values):
        """Write fitted level values ``[(kind "energy" | "width", group, level, pair key, l, 2s, value)]`` into the
        rows they come from (``_azure2.level_value_edits``; pyazr.Session.save_fit).  Returns the values written
        nowhere: ``[(index into values, why)]``."""
        replacements, unplaced = _engine_call("level values", _engine_module().level_value_edits, self._doc, values)
        _replace(self._doc, replacements)
        return list(unplaced)

    # -- <config> -----------------------------------------------------------------------------------------------

    def config(self):
        """The ``<config>`` entries as the engine reads them: ``{name: value}`` (texts), the eleven positional ones
        (``amatrix``, ``output_dir``, ``checks_dir``, ``check_compound_nucleus``, ... ``check_angular_distributions``)
        then the keys given (``relativistic``, ``priorCovariance``, ``format``, ...); the last line of a key wins."""
        entries, _ = _engine_call("<config>", _engine_module().config_entries, self._doc)
        return {e["name"]: _text(e["value"]) for e in entries}

    def set_config(self, **entries):
        """Set ``<config>`` entries by name (see :meth:`config`; the keys and their kinds:
        ``pyazr._azure2.config_keys()``): a value is the word or path the engine reads (``output_dir="out/"``,
        ``check_data="screen"``, ``amatrix="true"``), a number, or a bool (``amatrix`` true / false, a key of kind
        bool 1 / 0); ``None`` removes a key.  A directory is read with a trailing "/" whether or not it is written
        with one; :meth:`pyazr.Session.write_outputs` creates the output directory when it does not exist."""
        engine = _engine_module()

        def edit(doc):
            for name, value in entries.items():
                if isinstance(value, bool):  # the engine's words: amatrix true|false, a key of kind bool 1|0
                    value = ("true" if value else "false") if name == "amatrix" else ("1" if value else "0")
                text = "" if value is None else (value if isinstance(value, str) else repr(value))
                _replace(doc, _engine_call(f"<config> {name}", engine.config_edits, doc, name, text,
                                           remove=value is None))
        return self._apply("<config>", edit)

    # -- <thm> --------------------------------------------------------------------------------------------------

    def thm(self):
        """The ``<thm>`` section as written: ``{"options": {key: value}, "experiments": {name: {key: value}}}``
        (texts; the last line of an option wins); ``None`` without the section.  In format 3 an experiment's data
        sets name it (their ``thm`` field) and a data set's weight table is its ``thm_weight`` field."""
        exists, entries = _engine_call("<thm>", _engine_module().thm_entries, self._doc)
        if not exists:
            return None
        out = {"options": {}, "experiments": {}}
        for experiment, key, value, _ in entries:
            if experiment:
                out["experiments"].setdefault(experiment, {})[key] = value
            else:
                out["options"][key] = value
        return out

    def _edit_thm(self, what, changes, members=()):
        engine = _engine_module()

        def edit(doc):
            for call, args, remove in changes:
                _replace(doc, _engine_call(what, call, doc, *args, remove=remove))
            for dataset, experiment in members:
                self._set_fields(doc, dataset, {"thm": experiment})
            from .project import FormatError
            why = engine.check_format3_thm(doc.write())
            if why:
                raise FormatError(f"{what}: {why}")
        return self._apply(what, edit)

    def set_thm(self, options=None, /, **more):
        """Set global ``<thm>`` options (``vertex``, ``kinematics``, ``entranceL``, ``coulombIntegral``,
        ``spectatorEnergy``; a per-pair spectator energy ``spectatorEnergy[<pair key>]`` in the dict ``options``);
        ``None`` removes one.  The block is then checked as the engine reads it."""
        items = dict(options or {}, **more)
        engine = _engine_module()
        return self._edit_thm("<thm>", [(engine.thm_option_edits, (k, _thm_value(v)), v is None) for k, v in items.items()])

    def set_thm_experiment(self, name, /, *, datasets=None, **keys):
        """Set keys of ``experiment[name]`` (``background``, ``beam``, ``target``, ``spectator``, ``Ebeam``,
        ``cbackground``, ``distortion``, ...; ``None`` removes one); a new name adds the experiment.  Its data sets
        name it in their ``thm`` field; ``datasets`` (names) sets that field of each to ``name`` in the same edit, as
        a new experiment needs (it has no member until a data set names it, and a data set cannot name it before it
        exists)."""
        engine = _engine_module()
        return self._edit_thm(f"<thm> experiment[{name}]",
                              [(engine.thm_experiment_edits, (name, k, _thm_value(v)), v is None) for k, v in keys.items()],
                              [(d, name) for d in (datasets or ())])

    def remove_thm_experiment(self, name):
        """Remove ``experiment[name]`` (every line of it); refused while a data set names it."""
        engine = _engine_module()
        return self._edit_thm(f"remove experiment[{name}]", [(engine.thm_experiment_edits, (name, "", ""), True)])


def _thm_value(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    return v if isinstance(v, str) else repr(v)
