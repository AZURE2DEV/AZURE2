"""A prior or a limit on a level parameter, as a :class:`pyazr.Project` edit of ``<parameterSettings>``.

``<parameterSettings>`` (``docs/source/reference/data_formats.rst``) holds one row per parameter: a label, then
``value lower upper error fit_error nuisance category minuit_index``.  The engine reads it at the start of a fit, of a
sampler and of a session's objective (``ParameterLimitsManager::ReadParameterSettings``): a row whose ``nuisance`` is
1 adds the prior ((p - value)/error)^2 to the objective (the term "nuisance priors"), and a row whose limits are not
both 0 limits the parameter.  A level's energy is named ``Level N Energy (MeV)``, a level's channel ``Level N Channel C
Width (eV)``, with N and C the engine's numbers of the level and of its J^pi group's channel
(``ParameterLimitsManager::SettingNameFor``), the numbers of its parameter names ``energy_N`` and ``width_N_C``.  The
value, error and limits of a channel are in the convention its ``<levels>`` row states (a partial width in eV, an ANC
in fm^-1/2 for a level below the pair's threshold, a reduced-width amplitude in MeV^1/2 for a row flagged as one); the
engine converts them to the fit's amplitude.

The engine applies a row to a free parameter only: a prior or a limit on a fixed one is refused here.  N and C follow
the active levels: a level made inactive, added or removed renumbers the levels after it, and a row then names
another parameter (rewrite the rows after such an edit; :meth:`ParameterSettings.parameter_settings` lists them with
the level and channel they name now).

:class:`ParameterSetting` is one row; :data:`FIELDS` its fields with their units (``pyazr.kinds()
["parameter_settings"]``).
"""
import math
from dataclasses import dataclass
from typing import Optional, Tuple

from .errors import FormatError

#: The columns after a row's label, in the order the engine reads them (ParameterLimitsManager::ReadParameterSettings).
COLUMNS = ("value", "lower", "upper", "error", "fit_error", "nuisance", "category", "minuit_index")

#: The fields of a :class:`ParameterSetting`: title, unit, doc (``pyazr.kinds()["parameter_settings"]``).
FIELDS = {
    "label": {"title": "Row", "unit": "", "doc": "the row's label: 'Level N Energy (MeV)' or 'Level N Channel C "
                                                 "Width (eV)', the engine's numbers of the level and the channel"},
    "level": {"title": "Level", "unit": "", "doc": "the label of the level the row names now (None: no active level)"},
    "channel": {"title": "Channel", "unit": "", "doc": "(pair key, l, s) of the channel the row names; None for the "
                                                       "level's energy"},
    "value": {"title": "Prior centre", "unit": "MeV (energy), the channel's convention (width)",
              "doc": "the prior's centre (a row without prior: the value it was written with)"},
    "error": {"title": "Prior sigma", "unit": "as value", "doc": "the prior's standard deviation (> 0)"},
    "lower": {"title": "Lower limit", "unit": "as value", "doc": "the lower limit (0 and 0: no limit)"},
    "upper": {"title": "Upper limit", "unit": "as value", "doc": "the upper limit (0 and 0: no limit)"},
    "nuisance": {"title": "Prior", "unit": "", "doc": "the prior is applied (the column 'nuisance' = 1)"},
    "fit_error": {"title": "Fit error", "unit": "as value", "doc": "a fit's error, a report (not read as input)"},
    "category": {"title": "Category", "unit": "", "doc": "'level' for a level parameter"},
    "minuit_index": {"title": "Index", "unit": "", "doc": "the free index of an older file (the name wins; -1: none)"},
    "line": {"title": "Line", "unit": "", "doc": "the row's line in the file (1-based)"},
}


@dataclass(frozen=True)
class ParameterSetting:
    """One row of ``<parameterSettings>`` that names a level parameter (fields: :data:`FIELDS`).  ``prior`` is
    ``(centre, sigma)`` or None, ``limits`` ``(lower, upper)`` or None."""
    label: str
    level: Optional[str]
    channel: Optional[Tuple[int, int, float]]
    value: float
    lower: float
    upper: float
    error: float
    fit_error: float
    nuisance: bool
    category: str
    minuit_index: int
    line: int

    BUILT_BY = "Project.parameter_settings()"

    @property
    def prior(self):
        return (self.value, self.error) if self.nuisance else None

    @property
    def limits(self):
        return None if self.lower == 0.0 and self.upper == 0.0 else (self.lower, self.upper)


def _label(number, channel=None):
    """The row's label of level ``number``'s energy, or of its group's channel ``channel``."""
    if channel is None:
        return f"Level {number} Energy (MeV)"
    return f"Level {number} Channel {channel} Width (eV)"


def _number(x):
    return repr(float(x))


def _refuse(what, allowed, fix):
    return FormatError(f"{what}; allowed: {allowed}; {fix}", what=what, allowed=allowed, fix=fix)


def _numbered_levels(table):
    """``{(group, level): N}`` of the active levels of a ``level_table``: the engine numbers them 1, 2, ... in the
    order of their J^pi groups and, within a group, of their levels (the numbers of ``energy_N``)."""
    keys = sorted((d["group"], d["level"]) for d in table["levels"] if d["active"])
    return {k: n for n, k in enumerate(keys, start=1)}


class ParameterSettings:
    """The calls of :class:`pyazr.Project` on ``<parameterSettings>`` (mixed into Project): :meth:`parameter_settings`
    and :meth:`set_parameter_setting`."""

    def parameter_settings(self):
        """The rows of ``<parameterSettings>`` that name a level parameter (:class:`ParameterSetting`), in file
        order, each with the level and channel it names in the project as it stands."""
        from .nuclear_part import _text
        table = self._level_table()
        numbers = _numbered_levels(table)
        by_number = {numbers[(d["group"], d["level"])]: d for d in table["levels"] if d["active"]}
        out = []
        for _, first, last in _settings_sections(self._doc):
            for i in range(first + 1, last):
                parts = _text(self._doc.line(i)).split()
                if not parts or parts[0].startswith("#") or parts[0] != "Level":
                    continue
                n_label = 6 if len(parts) == 14 and parts[2] == "Channel" else 4
                if len(parts) != n_label + len(COLUMNS):
                    continue
                d = by_number.get(int(parts[1])) if parts[1].isdigit() else None
                channel = None
                if n_label == 6 and d is not None:
                    c = next((c for c in d["channels"] if c["channel"] == int(parts[3])), None)
                    channel = None if c is None else (c["pair"], c["l"], c["two_s"] / 2)
                v = parts[n_label:]
                out.append(ParameterSetting(" ".join(parts[:n_label]), d["label"] if d else None, channel,
                                            float(v[0]), float(v[1]), float(v[2]), float(v[3]), float(v[4]),
                                            int(v[5]) == 1, v[6], int(v[7]), i + 1))
        return out

    def set_parameter_setting(self, level, channel=None, *, prior=None, limits=None):
        """Put a prior and/or limits on a level parameter: the energy of ``level`` (a label, a ProjectLevel, (J^pi,
        E_x)), or with ``channel=(pair key, l, s)`` the width of that channel.  ``prior=(centre, sigma)`` adds
        ((p - centre)/sigma)^2 to the objective; ``limits=(lower, upper)`` limits the parameter in a fit and the
        sampler's support; both None remove the parameter's row.  Units: MeV for an energy; for a channel the
        convention its ``<levels>`` row states (eV, fm^-1/2 for an ANC, MeV^1/2 for an amplitude).

        The row replaces the parameter's row when it has one, else it is added (the section too).  Refused
        (FormatError with ``what``, ``allowed``, ``fix``; the project unchanged): a level or channel the project has
        not, an inactive level, a fixed parameter (the engine applies a row to a free one only), a sigma that is not
        > 0, limits that are not lower < upper, a number that is not finite.  Returns the row
        (:class:`ParameterSetting`), None when it was removed."""
        from .nuclear_part import _engine_call, _engine_module, _level_address, _level_name, _two_s
        name = _level_name(level)
        table = self._level_table()
        k = _engine_call(f"level {name!r}", _engine_module().find_level, self._doc, _level_address(level))
        d = table["levels"][k]
        what = f"level {d['label']!r}" + ("" if channel is None else f", channel {tuple(channel)!r}")
        if not d["active"]:
            raise _refuse(f"{what}: the level is inactive, so it has no parameters", "an active level",
                          f"Project.set_level({d['label']!r}, active=True) first")
        if channel is None:
            column, fixed, current = None, d["fixed"], d["energy"]
        else:
            try:
                pair, l, s = channel
                two_s = _two_s(s)
            except (TypeError, ValueError):
                raise _refuse(f"{what}: not a channel", "channel=(pair key, l, s), or None for the energy",
                              "Project.levels() lists each level's channels") from None
            c = next((c for c in d["channels"] if c["pair"] == int(pair) and c["l"] == int(l) and c["two_s"] == two_s),
                     None)
            if c is None:
                have = ", ".join(f"({c['pair']}, {c['l']}, {c['two_s'] / 2:g})" for c in d["channels"])
                raise _refuse(f"{what}: the level has no such channel", f"its channels (pair, l, s): {have}",
                              "Project.levels() lists each level's channels")
            column, fixed, current = c["channel"], c["fixed"], c["width"]
        if fixed and (prior is not None or limits is not None):
            fix = (f"Project.set_level({d['label']!r}, fixed=False)" if channel is None else
                   f"Project.set_channel({d['label']!r}, {channel[0]}, {channel[1]}, {channel[2]}, fixed=False)")
            raise _refuse(f"{what}: the parameter is fixed, and the engine applies a prior or a limit to a free "
                          f"parameter only", "a free level parameter", f"free it first: {fix}")
        if prior is not None:
            centre, sigma = (float(v) for v in prior)
            if not (math.isfinite(centre) and math.isfinite(sigma) and sigma > 0):
                raise _refuse(f"{what}: prior {tuple(prior)!r}", "prior=(centre, sigma), finite, sigma > 0",
                              "give the centre and the standard deviation")
        if limits is not None:
            lower, upper = (float(v) for v in limits)
            if not (math.isfinite(lower) and math.isfinite(upper) and lower < upper):
                raise _refuse(f"{what}: limits {tuple(limits)!r}", "limits=(lower, upper), finite, lower < upper",
                              "give the two limits in order")
        numbers = _numbered_levels(table)
        label = _label(numbers[(d["group"], d["level"])], column)
        row = None
        if prior is not None or limits is not None:
            centre, sigma = prior if prior is not None else (current, 0.0)
            lower, upper = limits if limits is not None else (0.0, 0.0)
            row = (f"{label} {_number(centre)} {_number(lower)} {_number(upper)} {_number(sigma)} 0 "
                   f"{1 if prior is not None else 0} level -1")

        def edit(doc):
            _replace_row(doc, label, row)

        self._apply(f"<parameterSettings> {label}", edit)
        if row is None:
            return None
        return next(r for r in self.parameter_settings() if r.label == label)


def _settings_sections(doc):
    return [(s[0], s[1], s[2]) for s in doc.sections() if s[0] == "parameterSettings" and s[2] >= 0]


def _replace_row(doc, label, row):
    """``row`` (text, or None: removed) in place of ``label``'s row of ``<parameterSettings>``; a new row goes before
    the section's closing line, a new section after the last section."""
    from .nuclear_part import _text
    from .project import _bytes, _insert_at
    sections = _settings_sections(doc)
    for _, first, last in sections:
        for i in range(first + 1, last):
            if " ".join(_text(doc.line(i)).split()[:len(label.split())]) == label:
                doc.replace_lines(i, i + 1, [] if row is None else [_bytes(row)])
                return
    if row is None:
        return
    if sections:
        last = sections[0][2]
        doc.replace_lines(last, last, [_bytes(row)])
        return
    at = _insert_at(doc)
    doc.replace_lines(at, at, [_bytes(x) for x in (
        "<parameterSettings>", "# Fitting parameter settings (only non-fixed parameters shown)",
        "# Format: name value lower_limit upper_limit error fit_error nuisance category minuit_index", row,
        "</parameterSettings>")])
