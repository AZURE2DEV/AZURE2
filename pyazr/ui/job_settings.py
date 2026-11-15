"""The settings of a job kind (:class:`Setting`, a row of :class:`pyazr.ui.jobs.JobRow`): their value types
(:data:`SETTING_TYPES`, one widget per type on the page: ``tabs/fitjob.js SETTING_WIDGETS``) and the conversion of a
request's values (:func:`merged_settings`).  A new type is a row here and a widget there."""
from dataclasses import dataclass

import numpy as np

from .errors import Refused


def _as_integer(value):
    if value is None or value == "":
        return None
    number = float(value)
    if not number.is_integer():
        raise ValueError(f"{value!r} is not an integer")
    return int(number)


def _as_boolean(value):
    if isinstance(value, bool):
        return value
    if str(value).lower() in ("true", "1", "yes"):
        return True
    if str(value).lower() in ("false", "0", "no"):
        return False
    raise ValueError(f"{value!r} is not true or false")


def _as_number(value):
    if value is None or value == "":
        return None
    number = float(value)
    if not np.isfinite(number):
        raise ValueError(f"{value!r} is not a finite number")
    return number


def _as_text(value):
    return "" if value is None else str(value).strip()


def _as_names(value):
    if value is None or value == "":
        return []
    words = value.replace(",", " ").split() if isinstance(value, str) else list(value)
    return [str(w) for w in words]


# The value types of a job's settings: type -> (conversion of a request's value, what the type allows).  The page has
# one widget per type (tabs/fit.js, SETTING_WIDGETS); a new type is a row here and a widget there.
SETTING_TYPES = {
    "integer": (_as_integer, "an integer, or empty for the default"),
    "boolean": (_as_boolean, "true or false"),
    "names": (_as_names, "names of free parameters, separated by blanks or commas; empty: every free parameter"),
    "number": (_as_number, "a finite number, or empty for the default"),
    "text": (_as_text, "a word or a path, as the setting's doc says"),
    "choice": (_as_text, "one of the setting's choices"),
}


@dataclass(frozen=True)
class Setting:
    """One setting of a job kind: its name, type (a key of :data:`SETTING_TYPES`), default and doc; a ``choice``
    has ``choices``, each (value, its words)."""
    name: str
    type: str
    default: object
    doc: str
    choices: tuple = ()



def merged_settings(row, given):
    """The row's defaults with the request's values over them, each converted by its type; Refused on an unknown
    setting or a value its type does not allow."""
    by_name = {x.name: x for x in row.settings}
    unknown = sorted(set(given or {}) - set(by_name))
    if unknown:
        raise Refused(f"job {row.kind}/{row.method} has no setting {', '.join(unknown)}", allowed=sorted(by_name))
    out = row.defaults()
    for name, value in (given or {}).items():
        setting = by_name[name]
        convert, allowed = SETTING_TYPES[setting.type]
        try:
            out[name] = convert(value)
        except (TypeError, ValueError) as e:
            raise Refused(f"job {row.kind}/{row.method}, setting {name}: {e}", allowed=allowed) from None
        if setting.choices and out[name] not in [c[0] for c in setting.choices]:
            raise Refused(f"job {row.kind}/{row.method}, setting {name}: {value!r} is not one of its choices",
                          allowed=", ".join(str(c[0]) for c in setting.choices),
                          fix=f"GET /api/jobs/kinds lists the choices of {name}")
    return out
