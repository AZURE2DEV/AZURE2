"""The physics functions of :mod:`pyazr.physics` as requests of the web UI (``POST /api/physics/{fn}``,
docs/upgrade/web_ui.md sections 4.2 and 6.6), and the session's own functions (``POST /api/session/functions/{fn}``).

The arguments of a request are checked against the function's row (``pyazr.kinds()["physics"]``: names, kinds,
domains), so that a refusal says what was wrong, what is allowed and which call shows it; the function is then called
with ``errors="mask"``: an element outside its domain is NaN in the answer with its status, not a refusal.

An argument is a number, a (nested) list of numbers (numpy's broadcasting rules), a choice by its name or a list of
names, or an axis ``{"from", "to", "points", "spacing"}`` that the server turns into an array with numpy
(``linspace``, ``geomspace``, ``arange`` for ``"integer"``), so that the page makes no number.  The answer carries the
arguments as used (``arguments``) and ``python``: the text of the pyazr call that makes the same arrays.

Where it runs: in the server's process, on the request's thread.  These functions need no session, so they must not
wait for the session worker (a build or an evaluation can hold it for seconds); the element-wise C++ loop releases
the GIL (``_azure2.physics_evaluate``), so other requests are served meanwhile.  The work of one request is bounded:
at most :data:`MAX_ELEMENTS` elements of the broadcast arguments, about 1 s for the Coulomb functions (coulsm, 45 us
per element) and 5 s for ``shift_derivative`` (five evaluations per element).

The session's functions (``channel_functions`` and the THM diagnostics) are pyazr's rows
(``Session.describe()["functions"]``): their arguments, their returned fields and where each sits in the answer.  The
page keeps only its own part of each (:data:`PAGE`: the tab, the group, the bound of one request, the status words,
the fields left off the first plot); :func:`answer_parts` walks the session's answer by the row, so that a new
diagnostic of the session needs no line here.
"""
import copy

import numpy as np

from .errors import NotFound, Refused

MAX_ELEMENTS = 20000
# Kinds of argument that take one text per call (not broadcast): a material's formula, a THM experiment's name.
TEXT_KINDS = ("text", "experiment")
SPACINGS = {"linear": "np.linspace", "log": "np.geomspace", "integer": "np.arange"}

# An argument's fields the page reads, where a row does not give them.
ARG = {"unit": "", "lower": None, "upper": None, "lower_open": False, "upper_open": False, "choices": [],
       "default": None, "example": None, "sweep": None}
# The kind of the argument the server fills itself: the fit vector of a row that depends on it (``at_x``), the
# working x of the page.
FIT_VECTOR_KIND = "vector"
NOT_FINITE = {"-1": "not a finite number at this energy"}

# The page's own part of each of the session's functions; the rows themselves are pyazr's
# (``Session.describe()["functions"]``, :mod:`pyazr.session_functions`): arguments, returned fields with their units,
# where each sits in the answer and which fields label a line.  Per row: ``tab`` (the tab that draws it; none: the
# Physics tab), ``group`` (its group in the Physics tab's list), ``max_elements`` (the energies one request may ask
# for, where one costs more than a physics function's), ``axis_field`` (the returned field that repeats the axis;
# default: the field named as the axis argument), ``status`` (the returned field whose False marks an element with
# the code, and the words of the code; without one an element not finite in every field is -1), ``hidden`` (fields
# left off the first plot), ``tables`` (a function that makes the tables of an answer without an axis).  A row of
# pyazr's table not named here is drawn by the Physics tab with the defaults (:data:`PAGE_DEFAULTS`).
PAGE_DEFAULTS = {"tab": None, "group": "this project", "max_elements": MAX_ELEMENTS, "axis_field": None,
                 "status": (None, -1, NOT_FINITE["-1"]), "hidden": (), "tables": None}
PAGE = {
    "channel_functions": {"status": ("open", -2, "closed channel (E <= 0 or the functions undefined there)")},
}


def check_size(shape, what, fix, bound=MAX_ELEMENTS):
    """Refuse a request whose broadcast ``shape`` has more than ``bound`` elements (:data:`MAX_ELEMENTS`, or a row's
    own ``max_elements`` where one element costs more)."""
    n = int(np.prod(shape, dtype=np.int64))
    if n > bound:
        raise Refused(f"{what}: {n} elements (shape {list(shape)}) is above the bound of {bound} per request; {fix}",
                      what=f"{n} elements in one request", allowed=f"at most {bound} elements per request", fix=fix)
    return n


def _number(x):
    """A number as Python text that reads back as the same double (an integral value as an int)."""
    if isinstance(x, bool):
        return repr(x)
    x = float(x)
    if x == int(x) and abs(x) < 2.0 ** 53 and not (x == 0 and np.signbit(x)):
        return repr(int(x))
    return repr(x)


def _axis(fn, name, spec, fix):
    """An axis {"from", "to", "points", "spacing"} as (array, Python expression)."""
    allowed = (f"{name}: {{\"from\": number, \"to\": number, \"points\": 2 to {MAX_ELEMENTS}, \"spacing\": "
               f"{' / '.join(SPACINGS)}}} (log: from and to > 0; integer: whole numbers, every one from to to)")
    unknown = sorted(set(spec) - {"from", "to", "points", "spacing"})
    spacing = spec.get("spacing", "linear")
    try:
        lo, hi = float(spec["from"]), float(spec["to"])
        points = int(spec.get("points", 0)) if spacing != "integer" else 0
    except (KeyError, TypeError, ValueError):
        lo = hi = points = None
    ok = (not unknown and spacing in SPACINGS and lo is not None and np.isfinite(lo) and np.isfinite(hi)
          and (spacing == "integer" or 2 <= points <= MAX_ELEMENTS)
          and (spacing != "log" or (lo > 0 and hi > 0))
          and (spacing != "integer" or (lo == int(lo) and hi == int(hi) and lo <= hi)))
    if not ok:
        raise Refused(f"{fn}: the axis {name} = {spec!r} is not one the server makes", what=f"{name} = {spec!r}",
                      allowed=allowed, fix=fix)
    if spacing == "integer":
        lo, hi = int(lo), int(hi)
        check_size((hi - lo + 1,), fn, "a shorter range")
        return np.arange(lo, hi + 1), f"np.arange({lo}, {hi + 1})"
    make = np.linspace if spacing == "linear" else np.geomspace
    return make(lo, hi, points), f"{SPACINGS[spacing]}({_number(lo)}, {_number(hi)}, {points})"


def _array_text(array):
    flat = ", ".join(_number(v) for v in np.ravel(array))
    return f"np.array([{flat}])" + (f".reshape{tuple(array.shape)}" if array.ndim != 1 else "")


def _all_strings(value):
    if isinstance(value, str):
        return True
    return isinstance(value, list) and len(value) > 0 and all(_all_strings(v) for v in value)


def prepare(fn, row, arguments):
    """The arguments of a request checked against ``row``: (keyword arguments for the call, the arguments as used
    for the answer, the Python lines that make the non-scalar ones, the call's keyword texts, broadcast shape)."""
    names = [a["name"] for a in row["arguments"]]
    signature = ", ".join(f"{a['name']} ({a['unit'] or 'no unit'}, {a['domain']})" for a in row["arguments"])
    fix = (f"GET /api/physics (or /api/stopping) gives the row of {fn}; send {{\"arguments\": {{name: value, list or axis "
           f"{{from, to, points, spacing}}}}}}")
    unknown = [k for k in arguments if k not in names]
    missing = [a["name"] for a in row["arguments"]
               if a.get("required", a["default"] is None) and a["name"] not in arguments]
    if unknown or missing:
        what = "; ".join(x for x in (f"unknown argument(s) {unknown}" if unknown else "",
                                     f"missing argument(s) {missing}" if missing else "") if x)
        raise Refused(f"{fn}: {what}", what=what, allowed=f"{fn}({signature})", fix=fix)
    kwargs, used, lines, texts, shapes = {}, {}, [], [], []
    for a in row["arguments"]:
        name = a["name"]
        if name not in arguments:
            used[name] = a["default"]
            continue
        value = arguments[name]
        if a["kind"] in TEXT_KINDS:
            if not isinstance(value, str):
                raise Refused(f"{fn}: {name} must be a text ({type(value).__name__} given)", what=f"{name} = {value!r}",
                              allowed=f"{name}: one text ({a['doc']})", fix=fix)
            kwargs[name] = used[name] = value
            texts.append(f"{name}={value!r}")
            continue
        if isinstance(value, dict):
            array, expr = _axis(fn, name, value, fix)
            kwargs[name], used[name] = array, array.tolist()
            lines.append(f"{name} = {expr}")
            texts.append(f"{name}={name}")
            shapes.append(array.shape)
            continue
        if a["kind"] == "choice" and _all_strings(value):
            given = np.asarray(value, dtype=object)
            bad = [v for v in np.ravel(given) if v not in a["choices"]]
            if bad:
                raise Refused(f"{fn}: {name} = {bad[0]!r} is not one of {a['choices']}", what=f"{name} = {bad[0]!r}",
                              allowed=f"{name} {a['domain']}", fix=fix)
            used[name] = value
            if given.ndim == 0:
                kwargs[name] = value
                texts.append(f"{name}={value!r}")
                shapes.append(())
                continue
            index = np.vectorize(a["choices"].index, otypes=[float])(given)
            kwargs[name] = index
            lines.append(f"{name} = {_array_text(index.astype(int))}  # {', '.join(map(str, np.ravel(given)))} "
                         f"(indices of the choices of {name})")
            texts.append(f"{name}={name}")
            shapes.append(index.shape)
            continue
        try:
            array = np.asarray(value, dtype=float)
        except (TypeError, ValueError):
            raise Refused(f"{fn}: {name} is not a number or a list of numbers ({type(value).__name__})",
                          what=f"{name} of type {type(value).__name__}",
                          allowed=f"{name}: a number or a (nested) list of numbers, or an axis, {a['domain']}",
                          fix=fix) from None
        kwargs[name], used[name] = value, value
        if array.ndim == 0:
            texts.append(f"{name}={_number(value)}")
        else:
            lines.append(f"{name} = {_array_text(array)}")
            texts.append(f"{name}={name}")
        shapes.append(array.shape)
    try:
        shape = np.broadcast_shapes(*shapes)
    except ValueError:
        raise Refused(f"{fn}: the argument shapes {shapes} do not broadcast against each other",
                      what=f"argument shapes {shapes}", allowed="shapes that broadcast (numpy rules)",
                      fix="give lists of one length, or scalars") from None
    check_size(shape, fn, "split the request", row.get("max_elements", MAX_ELEMENTS))
    return kwargs, used, lines, texts, shape


def _plain(v):
    return np.asarray(v).tolist() if np.ndim(v) else float(v)


def _status_words(row, code):
    """Words of a status code of a row (pyazr.physics.status_text's rule, for any table's rows)."""
    if code > 0:
        numeric = [a for a in row["arguments"] if a["kind"] not in TEXT_KINDS]
        arg = numeric[code - 1]
        return f"{arg['name']} outside its domain ({arg['domain']})"
    import pyazr
    return pyazr.physics.STATUS.get(code, f"unknown status {code}")


def _physics_call(fn, kwargs, texts):
    import pyazr
    result, status = pyazr.physics.functions[fn](**kwargs, errors="mask")
    return result, status, f"result, status = pyazr.physics.{fn}({', '.join(texts)}, errors=\"mask\")"


def _stopping_call(fn, kwargs, texts):
    import pyazr
    formula = kwargs.pop("formula")
    try:
        result, status = pyazr.stopping.evaluate(fn, formula, **kwargs)
    except ValueError as err:  # a formula the library refuses
        raise Refused(f"{fn}: {err}", what=f"formula = {formula!r}", allowed=str(err),
                      fix="a chemical formula: element symbols with counts (Ta2O5, TiN, CH4)") from None
    words = [f"{fn!r}", *texts]
    return result, status, f"result, status = pyazr.stopping.evaluate({', '.join(words)})"


# The function tables served as /api/<table>/{fn}: the rows (pyazr.kinds()[table]) and the call of a row.
TABLES = {"physics": _physics_call, "stopping": _stopping_call}


def call(fn, arguments, table="physics"):
    """A row of ``pyazr.kinds()[table]`` (``physics``: ``pyazr.physics.<fn>(..., errors="mask")``; ``stopping``:
    ``pyazr.stopping.evaluate``) on ``arguments``, as plain data: the returned fields as lists (floats for scalar
    arguments), their units, the status per element and the words of each status that occurs, the arguments as used,
    the Python text of the call, and the row's citation where it has one."""
    import pyazr
    rows = pyazr.kinds()[table]
    if fn not in rows:
        raise NotFound(f"no {table} function {fn!r}", what=f"the {table} function {fn!r}",
                       allowed=", ".join(rows), fix=f"GET /api/{table} lists the functions", known=list(rows))
    row = rows[fn]
    kwargs, used, lines, texts, shape = prepare(fn, row, arguments)
    result, status, line = TABLES[table](fn, kwargs, texts)
    fields = result if isinstance(result, dict) else {row["returns"][0]["name"]: result}
    codes = sorted({int(c) for c in np.ravel(status)})
    cite = [f"# cite: {row['citation']}"] if row.get("citation") else []
    python = "\n".join(["import numpy as np", "import pyazr", "", *lines, line, *cite, ""])
    answer = {"function": fn, "shape": list(shape), "axis": row["axis"],
              "returns": {k: _plain(v) for k, v in fields.items()},
              "units": {r["name"]: r["unit"] for r in row["returns"]},
              "status": np.asarray(status).tolist() if np.ndim(status) else int(status),
              "status_words": {str(c): _status_words(row, c) for c in codes},
              "arguments": used, "python": python}
    if row.get("citation"):
        answer["citation"] = row["citation"]
    return answer


def _pages():
    """The page's part of every session function (:data:`PAGE` and the THM tab's rows, pyazr/ui/thm.py)."""
    from . import thm
    return {**PAGE, **thm.PAGE}


def _pyazr_rows():
    from pyazr.session_functions import session_functions
    return session_functions()


def _top_level(r):
    """A returned field that is a key of the answer itself (not under a list or an experiment's name)."""
    return r["path"] == r["name"]


def _is_fact(r, arguments):
    return _top_level(r) and r["shape"] in ("scalar", "text") and r["name"] not in arguments


def page_row(name, row, page=None):
    """The row the page draws (the shape of a physics row) from pyazr's row of a session function: its arguments but
    the fit vector (the server gives the working x), each marked required where the session's method has no default;
    its returned fields on the axis (``returns``: name, unit, doc, ``lines`` where one field runs over channels or
    levels, ``shown`` False for a field left off the first plot), its scalars (``facts``), and the page's fields."""
    from .engine import required_arguments
    page = dict(PAGE_DEFAULTS, **(page if page is not None else _pages().get(name, {})))
    required = required_arguments(name)
    arguments = [dict(ARG, **a, required=a["name"] in required) for a in row["arguments"]
                 if a["kind"] != FIT_VECTOR_KIND]
    names = {a["name"] for a in row["arguments"]}
    skipped = {page["axis_field"] or row["axis"], page["status"][0]}
    returns = []
    for r in row["returns"]:
        if r["shape"] != "axis" or r["name"] in skipped:
            continue
        field = {"name": r["name"], "unit": r["unit"], "doc": r["doc"]}
        if r.get("lines"):
            field["lines"] = dict(r["lines"])
        if r["name"] in page["hidden"]:
            field["shown"] = False
        returns.append(field)
    facts = [{"name": r["name"], "unit": r["unit"], "doc": r["doc"]} for r in row["returns"] if _is_fact(r, names)]
    code, words = page["status"][1], page["status"][2]
    return dict(row, arguments=arguments, returns=returns, facts=facts, group=page["group"], tab=page["tab"],
                max_elements=page["max_elements"], status_words={"0": "computed", str(code): words})


def session_rows():
    """The rows of the session's functions as the page draws them (:func:`page_row` of each of pyazr's rows)."""
    return {name: page_row(name, row) for name, row in _pyazr_rows().items()}


def session_request(fn, arguments):
    """Check a request of a session function: (row, keyword arguments for the session's method, the arguments as
    used, the Python lines, the call's keyword texts)."""
    rows = session_rows()
    if fn not in rows:
        raise NotFound(f"no session function {fn!r}", what=f"the session function {fn!r}",
                       allowed=", ".join(rows), fix="GET /api/session/functions lists them", known=list(rows))
    row = rows[fn]
    kwargs, used, lines, texts, _ = prepare(fn, row, arguments)
    several = [a["name"] for a in row["arguments"] if a["name"] != row["axis"] and np.ndim(kwargs.get(a["name"], 0))]
    if several:
        raise Refused(f"{fn}: one value each for {several}", what=f"several values of {several}",
                      allowed=f"one value of each argument but the axis {row['axis']}", fix="one request per value")
    for a in row["arguments"]:
        if a["kind"] in ("integer", "pair") and a["name"] in kwargs and a["name"] != row["axis"]:
            kwargs[a["name"]] = int(kwargs[a["name"]])
    if row["axis"] in kwargs:
        kwargs[row["axis"]] = np.asarray(kwargs[row["axis"]], float).tolist()
    return row, kwargs, used, lines, texts


def _chains(out, over):
    """The elements a field's lines run over (``lines.over`` of its row, as ``exits[].levels[]`` or
    ``{experiment}``), each as the chain of objects from the answer down to the element (the element last)."""
    chains = [[out]]
    for part in over.split("."):
        step = []
        for chain in chains:
            node = chain[-1]
            if part.startswith("{"):
                step += [chain + [dict(v, **{part[1:-1]: k})] for k, v in node.items()]
            else:
                step += [chain + [v] for v in node.get(part[:-2], [])]
        chains = step
    return chains


def _label(chain, names, units):
    """The label of one line: each labelling field (``lines.label``) with its value and unit, from the element up."""
    parts = []
    for n in names:
        v = next((c[n] for c in reversed(chain) if isinstance(c, dict) and n in c), None)
        if v is None:
            continue
        text = f"{v:.6g}" if isinstance(v, float) else str(v)
        parts.append(f"{n} {text}" + (f" {units[n]}" if units.get(n) else ""))
    return ", ".join(parts)


def _finite_status(fields, n, code):
    """0 where every field (and every line of it) is finite at an element, ``code`` elsewhere."""
    ok = np.ones(n, bool)
    for v in fields.values():
        a = np.asarray(v, float).reshape(-1, n)
        ok &= np.isfinite(a).all(axis=0)
    return [0 if o else code for o in ok]


def answer_parts(fn, row, out, n, page=None):
    """The parts of the answer of a session function from the session's answer ``out`` (plain data), walked by
    pyazr's row ``row``: the fields on the axis (``returns``; a field under a list one line per element, labelled by
    its row's ``lines.label``), ``lines``, the status of each of the ``n`` elements of the axis, the call's scalars
    (``facts``) and, for a row without an axis, the page's ``tables``."""
    page = dict(PAGE_DEFAULTS, **(page if page is not None else _pages().get(fn, {})))
    names = {a["name"] for a in row["arguments"]}
    units = {r["name"]: r["unit"] for r in row["returns"]}
    skipped = {page["axis_field"] or row["axis"], page["status"][0]}
    returns, lines = {}, {}
    for r in row["returns"]:
        if r["shape"] != "axis" or r["name"] in skipped:
            continue
        if r.get("lines"):
            chains = [c for c in _chains(out, r["lines"]["over"]) if r["name"] in c[-1]]
            if chains:
                returns[r["name"]] = [c[-1][r["name"]] for c in chains]
                lines[r["name"]] = [_label(c, r["lines"]["label"], units) for c in chains]
        elif r["path"] in out:
            returns[r["name"]] = out[r["path"]]
    parts = {"returns": returns,
             "facts": [{"name": r["name"], "unit": r["unit"], "doc": r["doc"], "value": out[r["path"]]}
                       for r in row["returns"] if _is_fact(r, names) and r["path"] in out]}
    if lines:
        parts["lines"] = lines
    field, code, _ = page["status"]
    if not row["axis"]:
        parts["status"] = []
    elif field:
        parts["status"] = [0 if v else code for v in np.ravel(out[field])]
    else:
        parts["status"] = _finite_status(returns, n, code)
    if page["tables"]:
        made = page["tables"](row, out)
        parts["tables"] = made["tables"]
        parts["facts"] += made.get("facts", [])
    if not parts["facts"]:
        del parts["facts"]
    return parts


def session_answer(fn, row, out, used, lines, texts, opening):
    """The answer of a session function from the session's answer ``out``, in the shape of :func:`call`'s: the
    returned fields (a field over several channels or levels with its ``lines``), units, the status per element and
    its words, the arguments as used, the call's facts and tables, and the Python text (``opening``: the line that
    opens a session, ``engine.session_text``)."""
    n = len(used[row["axis"]]) if row["axis"] else 0
    parts = answer_parts(fn, _pyazr_rows()[fn], out, n)
    status = parts["status"]
    words = {int(k): v for k, v in row["status_words"].items()}
    note = ["# the page evaluates at its working x: params=x (the free vector) gives another point"] \
        if row.get("at_x") else []
    python = "\n".join(["import numpy as np", "import pyazr", "", opening, *lines, *note,
                        f"result = session.{fn}({', '.join(texts)})", ""])
    answer = {"function": fn, "shape": [len(status)], "axis": row["axis"], "returns": parts["returns"],
              "units": {r["name"]: r["unit"] for r in row["returns"]}, "status": status,
              "status_words": {str(c): words.get(c, f"status {c}") for c in sorted(set(status))}, "arguments": used,
              "python": python}
    for key in ("lines", "facts", "tables"):
        if key in parts:
            answer[key] = parts[key]
    return answer
