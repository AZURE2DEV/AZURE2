"""The options of a session (``docs/upgrade/project_format3.md`` section 3.12).

The options of the calculation (``use_brune``, ``transform``, ``use_coulomb_table``, ...) are keys of the project's
``<config>`` under the same names.  A keyword given to :class:`pyazr.Session` wins over the key, the key over the
default; a keyword left out is not given, even though the signature shows its default.  The engine applies that rule
when the session reads the project (``src/ProjectOptions.cpp``); this module tells it which keywords a call gave and
reads back the values in effect with their source (``"default"``, ``"project"``, ``"keyword"``).
"""
import functools
import inspect
import numbers

#: The keywords of :class:`pyazr.Session` that are not options of the calculation (not project state).
NOT_OPTIONS = ("self", "file", "cwd", "extrapolations", "coulomb_table", "parameter_file", "fold_depth_steps")


def records_given(init):
    """Wrap ``Session.__init__`` so that it knows which keywords the call named: ``self._given_options`` (a set).  The
    signature is unchanged (``inspect.signature`` follows the wrapper to it)."""
    signature = inspect.signature(init)

    @functools.wraps(init)
    def wrapper(self, *args, **kwargs):
        for name in kwargs:
            if name not in signature.parameters:
                removed = _removed_keyword(name)  # the engine's one text for an option that is gone
                if removed:
                    raise TypeError(removed)
        self._given_options = set(signature.bind(self, *args, **kwargs).arguments) - set(NOT_OPTIONS)
        return init(self, *args, **kwargs)
    return wrapper


def _removed_keyword(name):
    """The refusal of a keyword that was an option of the calculation and is one no more (the engine's table
    ``io::RemovedConfigKeys()``), or "" for any other name."""
    from . import _azure2
    return _azure2.removed_option_text(name)


def fold_depth_steps_value(value):
    """The engine's ``RuntimeOptions.fold_depth_steps`` of the keyword: -1 for ``None`` (not given), else the whole
    number from 0 to 1000 it is; anything else (a negative number, a fraction, ``nan``, ``inf``, a bool, text) is
    refused with a ``ValueError`` that says what is taken."""
    if value is None:
        return -1
    if isinstance(value, bool) or not isinstance(value, numbers.Integral) or not 0 <= value <= 1000:
        raise ValueError(f"fold_depth_steps needs a whole number from 0 to 1000 (0: the node intervals), not {value!r}")
    return int(value)


def runtime_options(engine, given, data_mode, coulomb_table, cwd, fold_depth_steps=None):
    """The engine's ``RuntimeOptions``: ``given`` ({keyword: value}, the options of the calculation the call named),
    the run's mode, the Coulomb table file of this machine, the directory of the project's relative paths and the
    panels of the target's depth integral (``None``: not given)."""
    opts = engine.RuntimeOptions()
    opts.data_mode = bool(data_mode)
    opts.coulomb_table = str(coulomb_table or "")
    opts.base_directory = cwd  # the model's relative paths: this directory, not the process's
    opts.fold_depth_steps = fold_depth_steps_value(fold_depth_steps)
    opts.options = {name: bool(value) for name, value in given.items()}
    return opts


def in_use(sess, extrapolations, coulomb_table, parameter_file=None, fold_depth_steps=None):
    """The options in effect, as ``(values, sources)``: ``values`` {keyword: value} with ``extrapolations``,
    ``coulomb_table`` and ``parameter_file`` (what :attr:`pyazr.Session.options` holds: passed back as keywords they rebuild the same
    session), ``sources`` {keyword: "default" | "project" | "keyword"} of the options of the calculation."""
    values = {"extrapolations": bool(extrapolations)}
    sources = {}
    for name, (value, source) in sess.option_sources().items():
        values[name], sources[name] = bool(value), source
    values["coulomb_table"] = str(coulomb_table or "")
    values["parameter_file"] = parameter_file
    values["fold_depth_steps"] = fold_depth_steps  # as given: None leaves it to $AZURE_FOLD_DEPTH_STEPS, else 8
    return values, sources
