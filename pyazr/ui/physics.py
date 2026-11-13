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
"""
import copy

import numpy as np

from .errors import NotFound, Refused

MAX_ELEMENTS = 20000
SPACINGS = {"linear": "np.linspace", "log": "np.geomspace", "integer": "np.arange"}

# The session's functions as rows of the shape of the physics table (pyazr.kinds()["physics"]): the session's methods
# have no table in the binding (web_ui.md section 12), so the row of each is written here from its docstring.  `op` is
# the adapter's operation (pyazr.ui.engine.Engine) that answers it; `broadcast` False: one value per argument, the
# axis an array.  A pair argument takes a key of the project's pairs (GET /api/project/pairs).
_ARG = {"unit": "", "lower": None, "upper": None, "lower_open": False, "upper_open": False, "choices": [],
        "default": None, "example": None, "sweep": None}
SESSION_ROWS = {
    "channel_functions": {
        "name": "channel_functions", "group": "this project", "op": "channel_functions", "broadcast": False,
        "summary": "The session's channel functions of a pair: Coulomb functions, P, S and the hard-sphere phase.",
        "doc": "As the session computes them (pyazr azure2.channel_functions): the Coulomb table when one is in use "
               "(about 1e-7 from the exact routine of the physics functions), GSL with use_gsl_coul, the hybrid "
               "model's Numerov solution for a pair with a nuclear potential.  Where the channel is closed (E <= 0, "
               "or the functions undefined there) the status is -2: F, G, P and delta_hs are 0 and S is the "
               "closed-channel shift function at the channel radius.",
        "library": "pyazr: azure2.channel_functions (the session's)", "axis": "e_cm",
        "arguments": [
            dict(_ARG, name="pair", kind="pair", doc="particle pair key of <levels>", domain="a pair of the project"),
            dict(_ARG, name="l", kind="integer", doc="orbital angular momentum", domain="integer >= 0", lower=0.0,
                 example=0.0, sweep={"lower": 0.0, "upper": 10.0, "spacing": "linear"}),
            dict(_ARG, name="radius", kind="real", unit="fm", doc="evaluation radius; 0: the pair's channel radius",
                 domain="finite >= 0", lower=0.0, default=0.0, sweep={"lower": 1.0, "upper": 10.0,
                                                                      "spacing": "linear"}),
            dict(_ARG, name="e_cm", kind="real", unit="MeV", doc="c.m. energy of the pair", domain="finite",
                 example=1.0, sweep={"lower": 0.01, "upper": 10.0, "spacing": "log"}),
        ],
        "returns": [{"name": n, "unit": u, "doc": d} for n, u, d in (
            ("F", "", "regular Coulomb function"), ("Fp", "", "dF/drho"), ("G", "", "irregular Coulomb function"),
            ("Gp", "", "dG/drho"), ("P", "", "penetrability"), ("S", "", "shift function"),
            ("delta_hs", "rad", "hard-sphere phase shift"))],
    },
}
SESSION_STATUS = {0: "computed", -2: "closed channel (E <= 0 or the functions undefined there)"}


def check_size(shape, what, fix):
    """Refuse a request whose broadcast ``shape`` has more than :data:`MAX_ELEMENTS` elements."""
    n = int(np.prod(shape, dtype=np.int64))
    if n > MAX_ELEMENTS:
        raise Refused(f"{what}: {n} elements (shape {list(shape)}) is above the bound of {MAX_ELEMENTS} per request; "
                      f"{fix}",
                      what=f"{n} elements in one request", allowed=f"at most {MAX_ELEMENTS} elements per request",
                      fix=fix)
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
    missing = [a["name"] for a in row["arguments"] if a["default"] is None and a["name"] not in arguments]
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
        if a["kind"] == "text":
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
    check_size(shape, fn, "split the request")
    return kwargs, used, lines, texts, shape


def _plain(v):
    return np.asarray(v).tolist() if np.ndim(v) else float(v)


def _status_words(row, code):
    """Words of a status code of a row (pyazr.physics.status_text's rule, for any table's rows)."""
    if code > 0:
        numeric = [a for a in row["arguments"] if a["kind"] != "text"]
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


def session_rows():
    """The rows of the session's functions (:data:`SESSION_ROWS`), a fresh copy."""
    return copy.deepcopy(SESSION_ROWS)


def session_request(fn, arguments):
    """Check a request of a session function: (row, keyword arguments for the adapter's operation, the arguments as
    used, the Python lines, the call's keyword texts)."""
    if fn not in SESSION_ROWS:
        raise NotFound(f"no session function {fn!r}", what=f"the session function {fn!r}",
                       allowed=", ".join(SESSION_ROWS), fix="GET /api/session/functions lists them",
                       known=list(SESSION_ROWS))
    row = SESSION_ROWS[fn]
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


def session_answer(fn, row, out, used, lines, texts, opening):
    """The answer of a session function from the adapter's ``out`` (dict of lists with ``open``), in the shape of
    :func:`call`'s, with the Python text: ``opening`` (the line that opens a session, ``engine.session_text``)."""
    status = [0 if o else -2 for o in np.ravel(out["open"])]
    python = "\n".join(["import numpy as np", "import pyazr", "", opening, *lines, f"result = session.{fn}({', '.join(texts)})", ""])
    return {"function": fn, "shape": [len(status)], "axis": row["axis"],
            "returns": {r["name"]: out[r["name"]] for r in row["returns"]},
            "units": {r["name"]: r["unit"] for r in row["returns"]}, "status": status,
            "status_words": {str(c): SESSION_STATUS[c] for c in sorted(set(status))}, "arguments": used,
            "python": python}
