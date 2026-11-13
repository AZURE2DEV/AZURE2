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
    ``light`` and ``heavy`` (Z, A); masses in u, energies in MeV, ``radius`` the channel radius in fm."""
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


def _replace(doc, replacements):
    for first, last, lines in replacements:
        doc.replace_lines(first, last, list(lines))


def _engine_call(what, call, *args, **kwargs):
    """The engine's replacements, its refusal (ValueError) as FormatError."""
    from .project import FormatError
    try:
        return call(*args, **kwargs)
    except ValueError as e:
        raise FormatError(f"{what}: {e}") from None


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
                     d["threshold"], d["radius"], d["two_j_light"], d["two_j_heavy"], d["heavy_parity"])
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

    def set_channel(self, level, pair, l, s, *, width=None, fixed=None):
        """Set the width column of channel (pair key, l, s) of a level -- in the row's convention: eV for an open
        channel or a photon, the ANC in fm^-1/2 for a closed channel, MeV^1/2 for a reduced-width amplitude row
        (:attr:`ProjectChannel.amplitude`) -- or whether it is fixed in a fit."""
        return self._edit_levels(f"level {_level_name(level)!r}, channel ({pair}, {l}, {s})", _engine_module().set_channel_edits,
                                 _level_address(level), int(pair), int(l), _two_s(s),
                                 width=None if width is None else float(width), fixed=fixed)

    def set_pair(self, pair, *, radius=None, light_mass=None, heavy_mass=None, separation_energy=None,
                 light_excitation=None, heavy_excitation=None):
        """Set a pair's channel radius (fm), masses (u), separation energy (MeV) or the excitation energies of its
        light and heavy partner (MeV) on every row of the pair (the separation energy also in the entrance-pair copy
        of the rows whose entrance it is).  The light partner's excitation is the optional last column of a row: a
        row that does not write it gets it only for a value other than 0.

        The channel radius is the matching surface of the R-matrix: a new radius changes every penetrability, shift
        function and external-capture integral, so the level parameters no longer describe a fit -- refit before
        reading anything off the result."""
        floats = {k: float(v) for k, v in
                  dict(radius=radius, light_mass=light_mass, heavy_mass=heavy_mass, separation_energy=separation_energy,
                       light_excitation=light_excitation, heavy_excitation=heavy_excitation).items() if v is not None}
        return self._edit_levels(f"pair {pair}", _engine_module().set_pair_edits, int(pair), **floats)

    def add_level(self, J, parity, energy, channels, *, fixed=True):
        """Add a level of spin ``J``, parity ``parity`` (+1 or -1) at E_x ``energy`` (MeV) with ``channels``, a list
        of dicts ``{pair, l, s, width (row convention, default 0), fixed (default False)}`` naming existing pairs.
        The levels of a J^pi share one channel set: a J^pi the project has gets the rows of its group, the channels
        not given with width 0, fixed (a channel outside the set is refused); a new J^pi gets one row per channel,
        the pair's columns copied from a row of the pair.  The level gets the next free level ID and goes at the end
        of ``<levels>``; ``fixed`` fixes its energy."""
        rows = [(int(c["pair"]), int(c["l"]), _two_s(c["s"]), float(c.get("width", 0.0)), bool(c.get("fixed", False)))
                for c in channels]
        return self._edit_levels(f"new level {J}{'+' if parity > 0 else '-'} at {energy} MeV",
                                 _engine_module().add_level_edits, _two_s(J), int(parity), float(energy),
                                 fixed=bool(fixed), channels=rows)

    def remove_level(self, level):
        """Remove every row of a level."""
        return self._edit_levels(f"remove level {_level_name(level)!r}", _engine_module().remove_level_edits, _level_address(level))

    def _write_level_values(self, values):
        """Write fitted level values ``[(kind "energy" | "width", group, level, pair key, l, 2s, value)]`` into the
        rows they come from (``_azure2.level_value_edits``; pyazr.azure2.save_fit).  Returns the values written
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
        with one; :meth:`pyazr.azure2.write_output_files` creates the output directory when it does not exist."""
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
