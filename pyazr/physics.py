"""Physics functions of the engine's libraries as free functions on arrays (docs/upgrade/pyazr3.md section 2.7).

One Python function per row of the binding's table (``api/src/pybind_physics.cpp``), built here at import from
``_azure2.physics_functions()``: its signature, docstring, units and domains are the row's, and each call evaluates
the row's library function element by element in C++ (no session, no project, no physics in Python).

Calling convention, the same for every function:

- arguments broadcast against each other (numpy rules); scalars in give scalars out, arrays give arrays of the
  broadcast shape;
- a function with one returned field returns it (a float or an array); one with several returns a dict of them;
- ``errors="raise"`` (the default): an element outside the domain raises :class:`pyazr.DomainError` naming the
  argument, its value, the element and the allowed range; a computation that does not converge raises
  :class:`pyazr.EngineError`;
- ``errors="mask"``: nothing is raised; such elements are NaN and the call returns ``(result, status)``, ``status`` an
  integer array of the broadcast shape: 0 computed, k > 0 argument k (1-based, in the order of the signature) outside
  its domain, -2 the arguments together outside the domain (each in its own: a lab angle beyond the kinematic
  maximum, an uncharged pair), -1 not converged.  :func:`status_text` turns codes into words.

``functions`` maps each name to its function; ``table()`` is the rows as plain data (``pyazr.kinds()["physics"]``
is the same), with the argument a plot puts on its horizontal axis (``axis``) and every returned field's unit.
A new function is a row of the binding's table and a test naming it (``pyazr/CLAUDE.md``, "Add a physics function").
"""
import copy
import inspect

import numpy as np

from .azure2 import _azure2
from .errors import DomainError, EngineError

STATUS = {0: "computed", -1: "not converged", -2: "the arguments together are outside the domain"}
_ROWS = {row["name"]: row for row in _azure2.physics_functions()}


def table():
    """The rows of the physics table as plain data: a dict by name of ``summary``, ``doc``, ``library`` (the C++
    function called), ``axis``, ``arguments`` (name, unit, doc, kind, domain in words, lower, upper, lower_open,
    upper_open, choices, default) and ``returns`` (name, unit, doc).  A fresh copy at each call."""
    return copy.deepcopy(_ROWS)


def status_text(name, status):
    """Words for status codes of function ``name`` (an int or an array of them): ``"computed"``, ``"e_cm outside
    its domain (finite > 0)"``, ``"not converged"``, ...  Returns a str or a list of str."""
    row = _ROWS[name]

    def one(code):
        code = int(code)
        if code > 0:
            arg = row["arguments"][code - 1]
            return f"{arg['name']} outside its domain ({arg['domain']})"
        return STATUS.get(code, f"unknown status {code}")
    return one(status) if np.ndim(status) == 0 else [one(c) for c in np.ravel(status)]


def _value_text(value, unit):
    return f"{value!r}" + (f" {unit}" if unit else "")


def _raise_for(row, status, columns, shape):
    """The error of the first element whose status is not 0 (columns: the flat arguments)."""
    name = row["name"]
    bad = np.flatnonzero(status != 0)
    i = int(bad[0])
    code = int(status[i])
    where = f"element {i}" + (f" (index {np.unravel_index(i, shape)})" if len(shape) > 1 else "")
    count = f"{bad.size} of {status.size} element(s) are not computed"
    fix = f"help(pyazr.physics.{name}) gives the domain; errors='mask' returns NaN with a status per element instead"
    if code > 0:
        arg = row["arguments"][code - 1]
        value = float(columns[code - 1][i])
        if arg["kind"] == "choice":
            value = arg["choices"][int(value)] if value == int(value) and 0 <= value < len(arg["choices"]) else value
        raise DomainError(f"{name}: {arg['name']} = {_value_text(value, arg['unit'])} at {where} is outside the "
                          f"domain; allowed: {arg['name']} {arg['domain']} ({arg['doc']}); {count}.  {fix}.")
    given = ", ".join(f"{a['name']} = {_value_text(float(columns[k][i]), a['unit'])}"
                      for k, a in enumerate(row["arguments"]))
    if code == -2:
        raise DomainError(f"{name}: the arguments at {where} ({given}) are outside the domain together, though each "
                          f"is in its own: {row['doc']}  {count}.  {fix}.")
    raise EngineError(f"{name}: the computation did not converge at {where} ({given}); {count}.  Please report it "
                      f"with these arguments; errors='mask' gives the other elements.")


def _argument_column(row, k, value):
    """One argument as a float array: a choice by its name, a flag from a bool."""
    arg = row["arguments"][k]
    if arg["kind"] == "choice" and isinstance(value, str):
        if value not in arg["choices"]:
            raise DomainError(f"{row['name']}: {arg['name']} = {value!r} is not known; allowed: {arg['name']} "
                              f"{arg['domain']}.  help(pyazr.physics.{row['name']}).")
        return np.asarray(float(arg["choices"].index(value)))
    try:
        return np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{row['name']}: {arg['name']} must be a number or an array of numbers "
                        f"({arg['unit'] or 'no unit'}, {arg['domain']}), got {type(value).__name__}") from None


def _call(row, values, errors):
    if errors not in ("raise", "mask"):
        raise ValueError(f"{row['name']}: errors = {errors!r}; allowed: 'raise' (DomainError, EngineError) or 'mask' "
                         f"(NaN and a status per element)")
    arrays = [_argument_column(row, k, v) for k, v in enumerate(values)]
    scalar = all(a.ndim == 0 for a in arrays)
    broadcast = np.broadcast_arrays(*arrays)
    shape = broadcast[0].shape
    columns = [np.ascontiguousarray(a, dtype=float).ravel() for a in broadcast]
    out, status = _azure2.physics_evaluate(row["name"], columns)
    if errors == "raise" and np.any(status != 0):
        _raise_for(row, status, columns, shape)
    fields = {}
    for j, ret in enumerate(row["returns"]):
        fields[ret["name"]] = float(out[j][0]) if scalar else out[j].reshape(shape)
    result = fields[row["returns"][0]["name"]] if len(fields) == 1 else fields
    if errors == "mask":
        return result, (int(status[0]) if scalar else status.reshape(shape))
    return result


def _docstring(row):
    lines = [row["summary"], "", row["doc"], "", "Parameters", "----------"]
    for a in row["arguments"]:
        unit = f", {a['unit']}" if a["unit"] else ""
        default = f", default {a['default']!r}" if a["default"] is not None else ""
        lines.append(f"{a['name']} : {a['kind']}{unit}{default}")
        lines.append(f"    {a['doc']}; domain: {a['domain']}.")
    lines.append("errors : 'raise' (default) or 'mask'")
    lines.append("    'mask': NaN outside the domain, and the call returns (result, status) (help(pyazr.physics)).")
    lines += ["", "Returns", "-------"]
    if len(row["returns"]) > 1:
        lines.append("dict of arrays of the broadcast shape (floats for scalar arguments):")
    for r in row["returns"]:
        lines.append(f"{r['name']} : {r['unit'] or 'dimensionless'}")
        lines.append(f"    {r['doc']}")
    lines += ["", "Raises", "------",
              "pyazr.DomainError", "    (a ValueError) an element outside the domain, named with its value.",
              "pyazr.EngineError", "    (a RuntimeError) a computation that does not converge.",
              "", f"Library: {row['library']} (element-wise in C++ through _azure2.physics_evaluate)."]
    return "\n".join(lines)


def _make(row):
    params = [inspect.Parameter(a["name"], inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                default=inspect.Parameter.empty if a["default"] is None else a["default"])
              for a in row["arguments"]]
    params.append(inspect.Parameter("errors", inspect.Parameter.KEYWORD_ONLY, default="raise"))
    signature = inspect.Signature(params)

    def function(*args, **kwargs):
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        return _call(row, [bound.arguments[a["name"]] for a in row["arguments"]], bound.arguments["errors"])

    function.__name__ = function.__qualname__ = row["name"]
    function.__module__ = __name__
    function.__doc__ = _docstring(row)
    function.__signature__ = signature
    return function


functions = {name: _make(row) for name, row in _ROWS.items()}
globals().update(functions)
__all__ = sorted(functions) + ["functions", "table", "status_text", "STATUS"]
