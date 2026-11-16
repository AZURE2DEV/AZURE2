"""The page's prior builder on :mod:`pyazr.priors`: what the screen opens with, the build, and the write that names
the file in ``<config>`` as one history entry.  The arithmetic and the checks are pyazr's (:func:`pyazr.priors.build`);
these functions hand it the session's parameter table and the project's directory.

Where the builder's parts are kept: in the prior file, on its ``#pyazr.prior`` comment line
(:func:`pyazr.priors.read_parts`), so that the builder reopens with them; nothing is kept in the browser or in the
server beyond the project's history.  Removing the prior is the ``<config>`` route's (``priorCovariance`` set to
null), with no call of its own.
"""
import os

from .errors import Refused
from ..paths import is_absolute

#: The ``<config>`` key that names the prior file.
PRIOR_KEY = "priorCovariance"

#: The file name the builder proposes when the project names none.
PROPOSED = "prior_parameters.txt"

def parameter_table(described):
    """The session's parameter table (:class:`pyazr.Parameters`, with its basis) from the server's description of
    the session (``Server.session_describe()``: the engine's records as :meth:`pyazr.ui.engine.Engine._parameters`
    hands them, the session's options)."""
    from pyazr.parameters import Parameter, Parameters
    from pyazr.session import basis_of
    renamed = {"x": "start", "physical_value": "value"}
    rows = []
    for r in described["parameters"]:
        record = {renamed.get(k, k): v for k, v in r.items() if k not in ("free", "x0")}
        record["fixed"] = not r["free"]
        lo, hi = r["bounds"]  # JSON's null where the engine has no limit
        record["bounds"] = (float("-inf") if lo is None else lo, float("inf") if hi is None else hi)
        rows.append(Parameter(record))
    return Parameters(rows, basis=basis_of(described["options"]))


def prior_open(state, files):
    """What the builder opens with: the key's file, its parts (from the file's parts line; ``hand_written`` when it
    has none), the file name to propose, the kinds of part and the project's files (``files``: the relative paths
    of the page's one file listing, ``Server.files``, with the program's ``output``)."""
    from pyazr import priors
    with state.lock:
        named = state.project.config().get(PRIOR_KEY)
        directory = state.directory
    saved, hand_written, missing = None, False, False
    if named:
        path = os.path.join(directory, named)
        if os.path.isfile(path):
            saved = priors.read_parts(path)
            hand_written = saved is None
        else:
            missing = True
    return {"key": PRIOR_KEY, "named": named, "hand_written": hand_written, "missing": missing,
            "parts": saved["parts"] if saved else [], "proposed": named or PROPOSED,
            "kinds": priors.kinds_table(), "files": files, "bases": priors.BASES}


def prior_build(state, described, parts):
    """The prior from ``parts`` on the session of the current text (:func:`pyazr.priors.build`): its summary."""
    from pyazr import priors
    return priors.build(parameter_table(described), parts, directory=state.directory).summary()


def prior_write(state, described, parts, path, rev=None):
    """Build the prior, write it at ``path`` (relative to the project's directory, inside the root) and set the
    ``priorCovariance`` key to it: one history entry, which undo takes back with the file (removed, or restored as it
    was).  Refused when the engine would refuse the matrix."""
    from pyazr import priors
    if not path or is_absolute(str(path)):
        raise Refused(f"{path!r}: a file name relative to the project's directory", allowed="a relative path",
                      fix="give the name the builder proposes, or another inside the project's directory")
    with state.lock:
        state.inside_root(path)
        prior = priors.build(parameter_table(described), parts, directory=state.directory)
        if prior.refusal:
            r = prior.refusal
            raise Refused(f"the prior is not written: {r['what']}", allowed=r["allowed"], fix=r["fix"])
        state.edit_with_files(rev, f"prior on {len(prior.names)} parameter(s) in {path}",
                              lambda p: p.set_config(**{PRIOR_KEY: str(path)}), {path: prior.text()})
        return {"rev": state.rev, "path": path, "names": prior.names}
