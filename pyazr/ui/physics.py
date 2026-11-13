"""The physics functions of :mod:`pyazr.physics` as requests of the web UI (``POST /api/physics/{fn}``,
docs/upgrade/web_ui.md section 4.2).

The arguments of a request are checked against the function's row (``pyazr.kinds()["physics"]``: names, kinds,
domains), so that a refusal says what was wrong, what is allowed and which call shows it; the function is then called
with ``errors="mask"``: an element outside its domain is NaN in the answer with its status, not a refusal.

Where it runs: in the server's process, on the request's thread.  These functions need no session, so they must not
wait for the session worker (a build or an evaluation can hold it for seconds); the element-wise C++ loop releases
the GIL (``_azure2.physics_evaluate``), so other requests are served meanwhile.  The work of one request is bounded:
at most :data:`MAX_ELEMENTS` elements of the broadcast arguments, about 1 s for the Coulomb functions (coulsm, 45 us
per element) and 5 s for ``shift_derivative`` (five evaluations per element).
"""
import numpy as np

from .errors import NotFound, Refused

MAX_ELEMENTS = 20000


def check_size(shape, what, fix):
    """Refuse a request whose broadcast ``shape`` has more than :data:`MAX_ELEMENTS` elements."""
    n = int(np.prod(shape, dtype=np.int64))
    if n > MAX_ELEMENTS:
        raise Refused(f"{what}: {n} elements (shape {list(shape)}) is above the bound of {MAX_ELEMENTS} per request; "
                      f"{fix}",
                      what=f"{n} elements in one request", allowed=f"at most {MAX_ELEMENTS} elements per request",
                      fix=fix)
    return n


def call(fn, arguments):
    """``pyazr.physics.<fn>(**arguments, errors="mask")`` as plain data: the returned fields as lists (floats for
    scalar arguments), their units, the status per element and the words of each status that occurs."""
    import pyazr
    rows = pyazr.kinds()["physics"]
    if fn not in rows:
        raise NotFound(f"no physics function {fn!r}", what=f"the physics function {fn!r}",
                       allowed=", ".join(rows), fix="GET /api/physics lists the functions", known=list(rows))
    row = rows[fn]
    names = [a["name"] for a in row["arguments"]]
    signature = ", ".join(f"{a['name']} ({a['unit'] or 'no unit'}, {a['domain']})" for a in row["arguments"])
    fix = f"GET /api/physics gives the row of {fn}; send {{\"arguments\": {{name: value or list}}}}"
    unknown = [k for k in arguments if k not in names]
    missing = [a["name"] for a in row["arguments"] if a["default"] is None and a["name"] not in arguments]
    if unknown or missing:
        what = "; ".join(x for x in (f"unknown argument(s) {unknown}" if unknown else "",
                                     f"missing argument(s) {missing}" if missing else "") if x)
        raise Refused(f"{fn}: {what}", what=what, allowed=f"{fn}({signature})", fix=fix)
    shapes = []
    for a in row["arguments"]:
        if a["name"] not in arguments:
            continue
        value = arguments[a["name"]]
        if a["kind"] == "choice" and isinstance(value, str):
            if value not in a["choices"]:
                raise Refused(f"{fn}: {a['name']} = {value!r} is not one of {a['choices']}",
                              what=f"{a['name']} = {value!r}", allowed=f"{a['name']} {a['domain']}", fix=fix)
            shapes.append(())
            continue
        try:
            array = np.asarray(value, dtype=float)
        except (TypeError, ValueError):
            raise Refused(f"{fn}: {a['name']} is not a number or a list of numbers ({type(value).__name__})",
                          what=f"{a['name']} of type {type(value).__name__}",
                          allowed=f"{a['name']}: a number or a (nested) list of numbers, {a['domain']}",
                          fix=fix) from None
        shapes.append(array.shape)
    try:
        shape = np.broadcast_shapes(*shapes)
    except ValueError:
        raise Refused(f"{fn}: the argument shapes {shapes} do not broadcast against each other",
                      what=f"argument shapes {shapes}", allowed="shapes that broadcast (numpy rules)",
                      fix="give lists of one length, or scalars") from None
    check_size(shape, fn, "split the request")
    result, status = pyazr.physics.functions[fn](**arguments, errors="mask")
    fields = result if isinstance(result, dict) else {row["returns"][0]["name"]: result}
    codes = sorted({int(c) for c in np.ravel(status)})
    return {"function": fn, "shape": list(shape), "axis": row["axis"],
            "returns": {k: np.asarray(v).tolist() if np.ndim(v) else float(v) for k, v in fields.items()},
            "units": {r["name"]: r["unit"] for r in row["returns"]},
            "status": np.asarray(status).tolist() if np.ndim(status) else int(status),
            "status_words": {str(c): pyazr.physics.status_text(fn, c) for c in codes}}
