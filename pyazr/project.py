"""A format-3 project as records (docs/upgrade/project_format3.md sections 4 and 6): :class:`Project`.

The project is the engine's lossless document (``_azure2.Document``): an edit replaces only the lines of the fields it
changes, and a project opened and saved without an edit is written back byte for byte.  Records are read, written and
resolved by the engine (``_azure2.record_sections``, ``record_lines``, ``record_references``); every edit is checked by
the engine's validator, the one ``AZURE2 --validate`` uses, on the whole project as it would stand after the edit, and
refused with :class:`FormatError` when it adds an error.  The ``add_<kind>`` methods are generated at import from the
engine's table (``_azure2.schema()``): pyazr spells no record kind, field, quantity, unit or grammar of its own.
"""
import inspect
import json
import numbers
import os
import re
import warnings
from collections import Counter

from .described import _engine as _engine_module
from .errors import FormatError, MissingFile, ParameterKeyError, _diagnostic_text, close_names, listed
from .nuclear_part import NuclearPart

__all__ = ["Project", "Record", "Nuisance", "FormatError", "records_from_csv"]

_ENCODING = ("utf-8", "surrogateescape")  # bytes of the file <-> str, every byte kept


def _text(b):
    return b.decode(*_ENCODING) if isinstance(b, bytes) else b


def _bytes(s):
    return s.encode(*_ENCODING) if isinstance(s, str) else s


# FormatError (the engine's refusal of a project or an edit) is pyazr.errors.FormatError, imported above.


class Nuisance:
    """A nuisance value of a format-3 field (docs/upgrade/project_format3.md section 3.3), written by the engine's one
    writer of a nuisance (``experiment::NuisanceText``, the converter's too).

    ``value``: the number (None for a mode without a value, as THM's profiled norm); ``mode``: the mode's word of the
    schema's list (``pyazr.schema()["lists"]["nuisance_modes"]``); ``unit``: a unit token of the field's dimension or
    None (a dimensionless nuisance); ``prior``: ``(centre, width)`` with the width a number in the nuisance's unit or a
    text ending in ``%`` (per cent of the centre); ``grid_max``: the largest value the sub-point grids must span.
    ``str(n)`` is the file's text; the text itself is always accepted wherever a Nuisance is.  Nothing is checked
    here: the record that holds it is validated by the engine.

    >>> str(Nuisance(1, "fit", prior=(1, "5%")))
    '1 fit prior 1 5%'
    """

    def __init__(self, value, mode, *, unit=None, prior=None, grid_max=None):
        self.value, self.mode, self.unit, self.prior, self.grid_max = value, mode, unit, prior, grid_max

    def __str__(self):
        centre, width = ("", "") if self.prior is None else (_number_text(self.prior[0]), _number_text(self.prior[1]))
        return _engine_module().nuisance_text(
            "" if self.value is None else _number_text(self.value), self.unit or "", str(self.mode), centre, width,
            "" if self.grid_max is None else _number_text(self.grid_max))

    def __repr__(self):
        return f"Nuisance({str(self)!r})"


def _number_text(x):
    """A Python number as the file's text (repr: the shortest decimal that reads back as the same double); a text is
    taken as written."""
    if isinstance(x, bool):
        raise TypeError("a boolean is not a value of a format-3 field: give the word of the field's list")
    if isinstance(x, numbers.Integral):
        return str(int(x))
    if isinstance(x, numbers.Real):  # numpy scalars too
        return repr(float(x))
    return str(x)


def _value_text(value):
    """The file's text of one field value (section 4.1): a text as written; a number; a :class:`Nuisance`; a tuple
    of parts joined by blanks, as ``(5e17, "atoms/cm2")``.  One function for every structured value."""
    if isinstance(value, (tuple,)):
        return " ".join(_value_text(v) for v in value)
    if isinstance(value, Nuisance):
        return str(value)
    return _number_text(value)


def _field_texts(fields):
    """(key, text) pairs from keyword arguments: None leaves a field out, a list gives a repeatable field once per
    element, in order."""
    out = []
    for key, value in fields.items():
        if value is None:
            continue
        for v in value if isinstance(value, list) else [value]:
            out.append((key, _value_text(v)))
    return out


class Record:
    """One record of a :class:`Project`, read live from the project's document.

    Fields read as attributes and as a dict: ``r.norm``, ``r["norm"]`` give the field's text as written (a list of
    texts for a field given more than once), None for a field of the record's kind that is not given.  The record's
    own kind and name are ``record_kind`` and ``record_name`` (``kind`` is a field of a target).  ``set(**fields)``
    edits it (:meth:`Project.set`), ``describe()`` shows what the engine makes of it.
    """

    def __init__(self, project, name):
        self._project, self._name = project, name

    @property
    def record_name(self):
        return self._name

    @property
    def record_kind(self):
        return self._raw()["kind"]

    @property
    def line(self):
        """The line of the record's header in the file (1-based)."""
        return self._raw()["line"]

    @property
    def fields(self):
        """``{key: text}`` in file order; a key given more than once maps to the list of its texts."""
        out = {}
        for f in self._raw()["fields"]:
            key, value = _text(f["key"]), _text(f["value"])
            if key in out:
                out[key] = (out[key] if isinstance(out[key], list) else [out[key]]) + [value]
            else:
                out[key] = value
        return out

    def _raw(self):
        return self._project._find(self._name)[1]

    def __getitem__(self, key):
        fields = self.fields
        if key in fields:
            return fields[key]
        if key in _SCHEMA_FIELDS.get(self.record_kind, ()):
            return None
        raise KeyError(f"{self.record_kind} {self._name} has no field {key!r} (the table's fields of a "
                       f"{self.record_kind}: {', '.join(_SCHEMA_FIELDS.get(self.record_kind, ()))})")

    def __getattr__(self, key):
        if key.startswith("_"):
            raise AttributeError(key)
        try:
            return self[key]
        except KeyError as e:
            raise AttributeError(str(e)) from None

    def __contains__(self, key):
        return key in self.fields

    def keys(self):
        return self.fields.keys()

    def items(self):
        return self.fields.items()

    def get(self, key, default=None):
        return self.fields.get(key, default)

    def set(self, /, **fields):
        """Edit fields of this record (:meth:`Project.set`); returns the record."""
        self._project.set(self._name, **fields)
        return self

    def describe(self, as_json=False):
        """What the engine makes of this record (:meth:`Project.describe`)."""
        return self._project.describe(self._name, as_json=as_json)

    @property
    def text(self):
        """The record's lines as they stand in the file."""
        raw = self._raw()
        return "\n".join(self._project._line(i) for i in range(raw["line"] - 1, raw["last_line"]))

    def __repr__(self):
        return f"<Record {self.record_kind} {self._name}>"


class Project(NuclearPart):
    """A format-3 project as records (docs/upgrade/project_format3.md section 4.1).

    Open a file with :meth:`open`, or start from the nuclear part of one with :meth:`new`; read records with
    :meth:`records` or ``p[name]``; add them with ``add_<kind>`` (one method per record kind of the engine's table:
    ``help(p.add_dataset)`` lists the fields with their help and units) or :meth:`add`; edit with :meth:`set`,
    :meth:`rename`, :meth:`remove`; check with :meth:`validate` and :meth:`describe`; write with :meth:`save`; run with
    :meth:`session`.  The nuclear part (:class:`pyazr.nuclear_part.NuclearPart`): :meth:`levels`, :meth:`pairs`,
    :meth:`set_level`, :meth:`set_channel`, :meth:`set_pair`, :meth:`add_level`, :meth:`remove_level`; :meth:`config`,
    :meth:`set_config`; :meth:`thm`, :meth:`set_thm`, :meth:`set_thm_experiment`, :meth:`remove_thm_experiment`.  Data files are named relative to the project's directory (the directory of :attr:`path`), as a
    session reads them.

    A field's value is the text the file holds (``norm="1 fit prior 1 8%"``), a number, a ``(value, "unit")`` tuple, a
    :class:`Nuisance`, or a list of these for a field given more than once; a ``#`` in a text starts the field's
    comment, as in the file.  Every edit is validated by the engine on the project as it would stand and refused with
    :class:`FormatError` (the engine's diagnostics and codes) when it adds an error; the project is then unchanged.
    References are checked as the project stands: add a target before the data set that names it.
    """

    def __init__(self, data, path):
        self._doc = _engine_module().Document.parse(_bytes(data))
        self.path = os.path.abspath(path)
        self._saved = None  # the bytes last read from or written to `path`
        self._side_files = {}  # path relative to the project directory -> bytes, written by save()

    # -- opening and saving ---------------------------------------------------------------------------------------

    @classmethod
    def open(cls, path):
        """Open the format-3 project at ``path``.  A format-1 or format-2 file is refused (FormatError, a ValueError)
        with the engine's conversion command (:func:`require_format3`)."""
        require_format3(path)
        with open(path, "rb") as f:
            data = f.read()
        p = cls(data, path)
        p._saved = data
        return p

    @classmethod
    def new(cls, path, *, nuclear_part):
        """A project to be saved at ``path``, started from the nuclear part of the format-3 project ``nuclear_part``:
        that file's text with every record removed (the bodies of its record sections emptied: data sets,
        normalizations, energy scales, extrapolations, targets, beams, detectors, model data, potentials), every other
        section kept as it is (``<config>``, ``<levels>``, ``<thm>``, ``<parameterSettings>``, ...).  Nothing is
        written until :meth:`save`.  A format-1/2 file is refused with the conversion command."""
        src = cls.open(nuclear_part)
        doc = src._doc
        for s in sorted(src._sections(), key=lambda s: -s["open"]):
            if s["close"] > 0:
                doc.replace_lines(s["open"], s["close"] - 1, [])
        return cls(doc.write(), path)

    @property
    def directory(self):
        """The directory data files are read from: that of :attr:`path`."""
        return os.path.dirname(self.path)

    def to_bytes(self):
        """The project's text as it stands, bytes."""
        return self._doc.write()

    @property
    def text(self):
        return _text(self.to_bytes())

    @property
    def unsaved(self):
        """True when the document differs from what was last read from or written to :attr:`path`."""
        return self._saved != self.to_bytes()

    def save(self, path=None):
        """Write the project to ``path`` (default :attr:`path`; a new path becomes the project's) and the side files
        an edit changed (a prior-covariance file whose parameter names :meth:`rename` changed, at its path relative to
        the project).  Returns the path."""
        if path is not None:
            self.path = os.path.abspath(path)
        data = self.to_bytes()
        with open(self.path, "wb") as f:
            f.write(data)
        for rel, content in self._side_files.items():
            with open(os.path.join(self.directory, rel), "wb") as f:
                f.write(content)
        self._side_files = {}
        self._saved = data
        return self.path

    def session(self, **options):
        """A calculation of the project as it stands: :class:`pyazr.azure2` on :attr:`path` with ``options`` (its
        keyword arguments: ``use_brune``, ``data_mode``, ...).  The engine reads the file, so the project is saved
        first when it has unsaved edits (:attr:`unsaved`); the session's directory is the project's."""
        from .azure2 import azure2
        if self.unsaved or self._side_files:
            self.save()
        return azure2(self.path, **options)

    # -- reading --------------------------------------------------------------------------------------------------

    def _line(self, i):
        return _text(self._doc.line(i))

    def _sections(self, doc=None):
        return _engine_module().record_sections(doc if doc is not None else self._doc)

    def _find(self, name, doc=None):
        """(section dict, record dict) of the record ``name``; KeyError when there is none."""
        for s in self._sections(doc):
            for r in s["records"]:
                if _text(r["name"]) == name:
                    return s, r
        names = [_text(r["name"]) for s in self._sections(doc) for r in s["records"]]
        raise ParameterKeyError(f"no record {name!r} in {self.path}{close_names(name, names)}; its records: "
                                f"{listed(names)}; Project.records() lists them by name.")

    def records(self):
        """``{name: Record}`` of every record, in file order."""
        return {_text(r["name"]): Record(self, _text(r["name"])) for s in self._sections() for r in s["records"]}

    def __getitem__(self, name):
        self._find(name)
        return Record(self, name)

    def __contains__(self, name):
        try:
            self._find(name)
            return True
        except KeyError:
            return False

    def __iter__(self):
        return iter(self.records())

    def __len__(self):
        return len(self.records())

    # -- what the engine says -------------------------------------------------------------------------------------

    def _diagnostics(self, doc=None):
        data = (doc if doc is not None else self._doc).write()
        refusal, diagnostics, skipped = _engine_module().validate_format3_text(data, self.path, self.directory)
        if refusal:
            raise FormatError(refusal)
        return json.loads(diagnostics), skipped

    def validate(self):
        """Every problem of the project as it stands, as ``AZURE2 --validate --json`` reports it (a list of dicts,
        docs/upgrade/project_format3.md section 6.2; empty: valid), data files read relative to :attr:`directory`.
        The checks that need a valid project and its ``<levels>`` are made too; when they cannot be, a warning says
        how many were not."""
        diagnostics, skipped = self._diagnostics()
        if skipped:
            warnings.warn(f"{self.path}: {skipped} build check(s) were not run (they need a valid project and its "
                          f"<levels>)", stacklevel=2)
        return diagnostics

    def describe(self, name=None, as_json=False):
        """What the engine makes of the record ``name`` (a data set, extrapolation, datum or potential), or of every
        record (``name=None``), as ``AZURE2 --describe [--json]`` prints it: the quantity, the observable, the
        reaction with its pair keys, the data file and its rows read and kept, the effects and the fit parameters.
        ``as_json``: a dict (for one record, its entry).  A project with an error raises FormatError; a record the
        engine does not describe (a target, beam, detector, normalization, energy scale) raises ValueError naming the
        records that refer to it."""
        if name is not None:
            self._find(name)
        refusal, diagnostics, description = _engine_module().describe_format3_text(
            self.to_bytes(), self.path, bool(as_json), self.directory, name or "")
        if refusal:
            raise FormatError(refusal)
        errors = [d for d in json.loads(diagnostics) if d["severity"] == "error"]
        if errors:
            raise FormatError(f"{self.path}: not described, the project has errors:", errors)
        if as_json:
            description = json.loads(description)
            if name is not None:
                entries = [e for v in description.values() for e in v]
                if not entries:
                    raise ValueError(self._not_described(name))
                return entries[0]
            return description
        if name is not None and not description.strip():
            raise ValueError(self._not_described(name))
        return description

    def _not_described(self, name):
        users = sorted({f"{_text(r['record'])}.{r['field']}" for r in self._references() if _text(r["name"]) == name})
        kind = self[name].record_kind
        return (f"the engine describes the records that compute something; {name} is a {kind}"
                + (f", used by {', '.join(users)}: describe those" if users else ""))

    def _references(self, doc=None):
        return _engine_module().record_references((doc if doc is not None else self._doc).write())

    # -- edits ----------------------------------------------------------------------------------------------------

    def _apply(self, what, edit):
        """Run ``edit(doc)`` on a copy of the document and keep the copy when the engine's validator finds no error
        that the project did not have before; else raise FormatError with the new errors (``what`` heads the
        message).  New warnings are issued as warnings and returned (a list of the engine's diagnostics)."""
        before, _ = self._diagnostics()
        doc = _engine_module().Document.parse(self._doc.write())
        edit(doc)
        problems = doc.problems()
        if problems:
            raise FormatError(f"{what}: the document's structure: " + "; ".join(problems))
        after, _ = self._diagnostics(doc)
        new = self._new_diagnostics(before, after)
        errors = [d for d in new if d["severity"] == "error"]
        if errors:
            raise FormatError(f"{what}:", errors)
        for d in new:
            warnings.warn(f"{what}: {_diagnostic_text(d)}", stacklevel=3)
        self._doc = doc
        return new

    @staticmethod
    def _new_diagnostics(before, after):
        """The diagnostics of ``after`` that ``before`` does not have, by (severity, code, record, field) and count:
        line numbers and messages that quote them move with any edit."""
        key = lambda d: (d.get("severity"), d.get("code"), d.get("record"), d.get("field"))
        extra = Counter(key(d) for d in after) - Counter(key(d) for d in before)
        out = []
        for d in after:
            if extra[key(d)] > 0:
                extra[key(d)] -= 1
                out.append(d)
        return out

    def add(self, kind, name, /, **fields):
        """Add a record ``kind name`` with ``fields`` (keyword arguments named by the field keys; see ``add_<kind>``)
        at the end of the first section of the table that holds the kind, creating the section when the project has
        none.  Refused (FormatError) when the engine finds an error the project did not have: a missing required
        field, a value of the wrong grammar or unit, a reference to a record that does not exist, a name in use.
        Returns the :class:`Record`."""
        lines, problems = _engine_module().record_lines(kind, name, _field_texts(fields))
        if problems:
            raise FormatError(f"{kind} {name!r}:", [{"severity": "error", "code": c, "record": name, "message": m}
                                                    for c, m in problems])
        section = next((s for s, kinds in _SCHEMA["sections"].items() if kind in kinds), None)

        def edit(doc):
            existing = [s for s in self._sections(doc) if s["section"] == section and s["close"] > 0]
            if existing:  # before the closing tag
                s = existing[0]
                body = [b""] if s["records"] else []
                doc.replace_lines(s["close"] - 1, s["close"] - 1, body + list(lines))
                return
            block = [_bytes(f"<{section}>")] + list(lines) + [_bytes(f"</{section}>")] if section else list(lines)
            doc.replace_lines(*([_insert_at(doc)] * 2), block)

        self._apply(f"{kind} {name!r}", edit)
        return Record(self, name)

    def set(self, name, /, **fields):
        """Edit fields of the record ``name``: each keyword replaces the lines of that field (its first line keeps its
        place, indentation and value column; a comment after the value stays unless the new text brings its own);
        a field the record lacks is added after its last field; None removes the field; a list gives a repeatable
        field once per element.  Only the lines of the fields named change.  Refused (FormatError) when the engine
        finds an error the project did not have, or when there is no record ``name``."""
        self._require(name)
        self._apply(f"{self[name].record_kind} {name!r}", lambda doc: self._set_fields(doc, name, fields))
        return Record(self, name)

    def _require(self, name):
        """FormatError (the refusal of every edit) when the project has no record ``name``."""
        if name not in self:
            raise FormatError(f"{self.path}: no record {name!r} (the records: {', '.join(self.records()) or 'none'})")

    def _set_fields(self, doc, name, fields):
        """The edit of :meth:`set` on ``doc`` (not validated; :meth:`set` and the edits that change a record's fields
        together with another section, as a new THM experiment with its data sets, run it inside ``_apply``)."""
        engine = _engine_module()
        s, r = self._find(name, doc)
        for key, value in fields.items():
            texts = [] if value is None else [_value_text(v) for v in (value if isinstance(value, list) else [value])]
            s, r = self._find(name, doc)
            own = [f for f in r["fields"] if _text(f["key"]) == key]
            new_lines = []
            for i, t in enumerate(texts):
                canonical, comment, problems = engine.value_text(key, t)
                if problems:
                    raise FormatError(f"{r['kind']} {name!r}:", [
                        {"severity": "error", "code": c, "record": name, "field": key, "message": m}
                        for c, m in problems])
                model = own[min(i, len(own) - 1)] if own else (r["fields"][-1] if r["fields"] else None)
                new_lines.append(_field_line(doc, model, key, canonical, comment if comment else
                                             (model["comment"] if own and i < len(own) else b"")))
            if own:
                # the repeats' lines take the new values in order; every other line between them (another field,
                # a whole-line comment) stays where it is; more values follow the last repeat, fewer drop the rest
                first, last = own[0]["line"], own[-1]["line"]
                lines_of_own = {f["line"]: k for k, f in enumerate(own)}
                out = []
                for ln in range(first, last + 1):
                    k = lines_of_own.get(ln)
                    if k is None:
                        out.append(doc.line(ln - 1))
                    elif k < len(new_lines):
                        out.append(new_lines[k])
                doc.replace_lines(first - 1, last, out + new_lines[len(own):])
            elif new_lines:
                doc.replace_lines(r["last_line"], r["last_line"], new_lines)

    def rename(self, old, new):
        """Rename the record ``old`` to ``new``, with every reference to it (the validator's resolution of names) and
        every parameter name made from it (``old.norm``, ``D.value[old]``, ``old@epoch...``) in ``<parameterSettings>``,
        ``<mcmc>`` and the prior-covariance file the project names (written by :meth:`save`).  Refused (FormatError)
        when the engine finds an error the project did not have (a name of the wrong syntax or already in use)."""
        self._require(old)
        refs = [r for r in self._references() if _text(r["name"]) == old]
        engine = _engine_module()
        names = {}  # parameter name before -> after, from the engine (experiment::ParameterNamesOf)

        def edit(doc):
            _, r = self._find(old, doc)
            spots = [(r["line"], r["name_col"])] + [(x["line"], x["col"]) for x in refs]
            for line, col in sorted(spots, reverse=True):  # right to left: earlier columns stay valid
                b = doc.line(line - 1)
                start = col - 1
                if b[start:start + len(_bytes(old))] != _bytes(old):
                    raise RuntimeError(f"line {line}: '{old}' not at column {col}")
                doc.replace_lines(line - 1, line, [b[:start] + _bytes(new) + b[start + len(_bytes(old)):]])
            try:  # the names the engine gives, before and after, inactive records taken as active (their rows apply when
                # they are switched on); the same parameters in the same order
                before = engine.parameter_names(self._doc.write(), True)
                after = engine.parameter_names(doc.write(), True)
            except ValueError:
                warnings.warn(f"rename {old!r}: the project has errors, so the engine names no parameters: "
                              f"<parameterSettings>, <mcmc> and the prior covariance are left as they are", stacklevel=3)
                return
            if len(before) != len(after):
                raise RuntimeError(f"rename {old!r}: {len(before)} parameter names before, {len(after)} after")
            names.update((_bytes(b), _bytes(a)) for b, a in zip(before, after) if b != a)
            unnamed = []  # rows that still hold the old name after the engine's names are renamed
            for s in _sections_named(doc, _PARAMETER_SECTIONS):
                for i in range(s[1] + 1, s[2]):
                    b = doc.line(i)
                    renamed = _rename_parameters(b, names)
                    if renamed != b:
                        doc.replace_lines(i, i + 1, [renamed])
                    if _bytes(old) in renamed:
                        unnamed.append(f"{i + 1}: {_text(renamed).strip()}")
            if unnamed:
                warnings.warn(f"rename {old!r}: these rows contain {old!r} but name no parameter the engine gives the "
                              f"project; left as they are: " + "; ".join(unnamed), stacklevel=3)

        self._apply(f"rename {old!r} to {new!r}", edit)
        cov = self._prior_covariance_renamed(names)
        if cov:
            self._side_files[cov[0]] = cov[1]
        return Record(self, new)

    def _prior_covariance_renamed(self, names):
        """(path relative to the project, renamed bytes) of the prior-covariance file whose names line holds a
        parameter renamed by ``names`` (bytes -> bytes); None when there is none to change."""
        rel = _config_value(self._doc, "priorCovariance")
        if not rel:
            return None
        rel = rel.strip('"')
        path = os.path.join(self.directory, rel)
        content = self._side_files.get(rel)
        if content is None:
            if not os.path.isfile(path):
                return None
            with open(path, "rb") as f:
                content = f.read()
        lines = content.split(b"\n")
        for i, line in enumerate(lines):
            if line.strip() and not line.lstrip().startswith(b"#"):  # the line of names (data_formats.rst)
                renamed = _rename_parameters(line, names)
                if renamed == line:
                    return None
                lines[i] = renamed
                return rel, b"\n".join(lines)
        return None

    def remove(self, name):
        """Remove the record ``name`` (its header, fields and the comments above it).  Refused (FormatError) while
        another record refers to it: the message names the referrers, the diagnostics are the engine's
        (E_UNKNOWN_REFERENCE of each)."""
        self._require(name)

        def edit(doc):
            _, r = self._find(name, doc)
            first, last = r["line"] - 1, r["last_line"]
            while first > 0 and _text(doc.line(first - 1)).lstrip().startswith("#"):
                first -= 1
            if first > 0 and not doc.line(first - 1).strip():
                first -= 1
            elif last < doc.num_lines() and not doc.line(last).strip():
                last += 1
            doc.replace_lines(first, last, [])

        try:
            self._apply(f"remove {name!r}", edit)
        except FormatError as e:
            users = sorted({f"{d.get('record')}.{d.get('field')}" for d in e.diagnostics
                            if d.get("found") == name})
            if not users:
                raise
            raise FormatError(f"remove {name!r}: it is referenced by {', '.join(users)}", e.diagnostics) from None

    def _replace_record_section(self, section, lines):
        """The body of the record section ``section`` replaced by ``lines`` (records as the engine wrote them, e.g.
        the session's potential records): the section removed when ``lines`` is empty, added where :meth:`add` adds
        one when the project has none; validated as every edit."""
        lines = [_bytes(x) for x in lines]

        def edit(doc):
            found = [s for s in self._sections(doc) if s["section"] == section and s["close"] > 0]
            if found:
                s = found[0]
                if lines:
                    doc.replace_lines(s["open"], s["close"] - 1, lines)
                else:
                    doc.replace_lines(s["open"] - 1, s["close"], [])
            elif lines:
                doc.replace_lines(*([_insert_at(doc)] * 2), [_bytes(f"<{section}>")] + lines + [_bytes(f"</{section}>")])

        return self._apply(f"<{section}>", edit)

    # -- parameter values (pyazr.azure2.save_fit) -----------------------------------------------------------------

    def _write_parameter_values(self, values):
        """Write ``{parameter name: value in the engine's unit}`` into the nuisances the engine says they come from
        (``_azure2.parameter_value_edits``: a data set's norm and shift terms, a shared normalization's value, an
        energy scale's offset and function parameters), each in the unit it is written in.  Returns the names written
        nowhere, with the reason."""
        edits, missing, _ = _engine_module().parameter_value_edits(self.to_bytes(), self.directory,
                                                                    {k: float(v) for k, v in values.items()})
        for line, col, old, new in sorted(edits, key=lambda e: (e[0], -e[1])):
            b = self._doc.line(line - 1)
            if b[col - 1:col - 1 + len(old)] != old:
                raise RuntimeError(f"{self.path}:{line}: the value {old!r} is not at column {col}")
            self._doc.replace_lines(line - 1, line, [b[:col - 1] + new + b[col - 1 + len(old):]])
        return list(missing)


def records_from_csv(path, kind="dataset", project=None):
    """Records of ``kind`` declared in a sheet (docs/upgrade/project_format3.md section 6.5): a CSV file whose header
    row holds field keys, a ``name`` column among them, and one row per record (``name, file, reaction, quantity,
    unit, columns, norm, ...``; an empty cell leaves the field out; blank rows and rows starting with ``#`` are
    skipped).  Every row is checked by the engine: the grammar of its values, and with ``project`` (a
    :class:`Project`) the whole project as it would stand with the rows added before it (references, units, required
    fields).  All refusals come at once: FormatError whose diagnostics are the engine's, each with ``row`` (the line of
    the file) and its message starting "row N:".  Returns ``[(kind, name, fields)]`` for :meth:`Project.add`::

        for kind, name, fields in pyazr.records_from_csv("datasets.csv", project=p):
            p.add(kind, name, **fields)
    """
    import csv
    engine = _engine_module()
    with open(path, newline="") as f:
        lines = [(n, row) for n, row in enumerate(csv.reader(f), start=1)
                 if any(c.strip() for c in row) and not row[0].lstrip().startswith("#")]
    if not lines:
        raise FormatError(f"{path}: no header row")
    header = [h.strip() for h in lines[0][1]]
    if "name" not in header:
        raise FormatError(f"{path}: the header row has no 'name' column (the keys: {', '.join(header)})")
    trial = Project(project.to_bytes(), project.path) if project is not None else None
    records, refusals = [], []

    def refuse(n, diagnostics):
        refusals.extend(dict(d, row=n, message=f"row {n}: {d.get('message', '')}") for d in diagnostics)

    for n, row in lines[1:]:
        cells = dict(zip(header, (c.strip() for c in row)))
        name = cells.pop("name", "")
        fields = {k: v for k, v in cells.items() if v != ""}
        _, problems = engine.record_lines(kind, name, list(fields.items()))
        if problems:
            refuse(n, [{"severity": "error", "code": c, "record": name, "message": m} for c, m in problems])
            continue
        if trial is not None:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    trial.add(kind, name, **fields)
            except FormatError as e:
                refuse(n, e.diagnostics)
                continue
        records.append((kind, name, fields))
    if refusals:
        raise FormatError(f"{path}: {len({d['row'] for d in refusals})} row(s) refused", refusals)
    return records


def require_format3(path):
    """Raise :class:`FormatError` unless the project file at ``path`` is format 3: the engine's refusal
    (io::NotFormat3Message; a format-1 or format-2 file, a format-3 header over an old section, a file that is not a
    project), which names the conversion ``AZURE2 --convert OLD NEW`` (:func:`pyazr.convert_azr`).  The one refusal of
    an old file in pyazr: every entry point that takes a path (a session, :meth:`Project.open`, the functions that
    open a project) calls it."""
    if not os.path.isfile(os.fspath(path)):
        raise MissingFile(path, _missing_project_message(path))
    with open(os.fspath(path), "rb") as f:
        data = f.read()
    message = _engine_module().not_format3_message(data, os.fspath(path))
    if message:
        raise FormatError(message)


def _missing_project_message(path):
    """What a missing project file needs: the path as given and as resolved, the .azr files beside it, the calls."""
    given = os.fspath(path)
    full = os.path.abspath(given)
    folder = os.path.dirname(full)
    what = "is a directory, not a project file" if os.path.isdir(full) else "does not exist"
    if os.path.isdir(folder):
        beside = sorted(f for f in os.listdir(folder) if f.endswith(".azr"))
        there = (f"the .azr files in {folder}: {', '.join(beside[:8])}{', ...' if len(beside) > 8 else ''}"
                 if beside else f"{folder} holds no .azr file")
    else:
        there = f"its directory {folder} does not exist either"
    return (f"project file {given!r} {what} (resolved: {full}; the working directory is {os.getcwd()}); a session "
            f"and Project.open need an existing format-3 project file ({there}).  Give its path, or convert an old "
            f"project with pyazr.convert_azr(OLD, NEW).")


def write_parameter_values(path, values):
    """Write ``{parameter name: value in the engine's unit}`` into the format-3 project at ``path``
    (:meth:`Project._write_parameter_values`; pyazr.azure2.save_fit).  Returns the names written nowhere."""
    p = Project.open(path)
    missing = p._write_parameter_values(values)
    p.save()
    return missing


# -- helpers of the edits -----------------------------------------------------------------------------------------

#: Sections whose lines name fit parameters (data_formats.rst: <parameterSettings> rows, <mcmc> priors).
_PARAMETER_SECTIONS = ("parameterSettings", "mcmc")


def _sections_named(doc, names):
    """(name, open, close) of every closed section of ``doc`` named in ``names`` (0-based line indices)."""
    return [(s[0], s[1], s[2]) for s in doc.sections() if s[0] in names and s[2] >= 0]


def _rename_parameters(line, names):
    """``line`` (bytes) with every whitespace-separated token that is a key of ``names`` replaced by its value: the
    parameter names are the engine's, compared whole."""
    return b"".join(names.get(t, t) for t in re.split(rb"(\s+)", line))


def _config_value(doc, key):
    """The value of a ``<config>`` entry as the engine reads it (``_azure2.config_entries``; the last line of a key
    wins), None when absent."""
    value = None
    try:
        entries = _engine_module().config_entries(doc)[0]
    except ValueError:  # no <config> section
        return None
    for e in entries:
        if e["name"] == key:
            value = _text(e["value"])
    return value


def _insert_at(doc):
    """Where a new section goes: after the last top-level section whose closing line has a line ending, else at the
    start of ``<levels>``, else at the end."""
    sections = [s for s in doc.sections() if s[3] == 0 and s[2] >= 0]
    for s in reversed(sections):
        if doc.ending(s[2]):
            return s[2] + 1
    levels = [s for s in sections if s[0] == "levels"]
    return levels[0][1] if levels else doc.num_lines()


def _field_line(doc, model, key, value, comment):
    """A field line ``key value [comment]`` laid out as ``model`` (a field dict of the same record: its indentation
    and value column), or as the record writer writes it when the record has no field to follow."""
    if model is None or not model["value_col"]:
        lines, _ = _engine_module().record_lines("x", "x", [(key, _text(value))])
        line = lines[1]
    else:
        src = doc.line(model["line"] - 1)
        indent = src[:model["col"] - 1]
        pad = max(1, model["value_col"] - model["col"] - len(_bytes(key)))
        line = indent + _bytes(key) + b" " * pad + value
    return line + (b"  " + comment if comment else b"")


# -- the add_<kind> methods, generated from the engine's table ----------------------------------------------------

_SCHEMA = _engine_module().schema()
_SCHEMA_FIELDS = {kind: list(k["fields"]) for kind, k in _SCHEMA["records"].items()}


def _required_always(field):
    return field.get("required") is True and "applies_when" not in field and not field.get("by_quantity")


def _field_doc(key, f):
    parts = [f.get("title", key)]
    if f.get("grammar"):
        parts.append(f["grammar"])
    if _required_always(f):
        parts.append("required")
    elif f.get("required") and f.get("applies_when"):
        parts.append("required when " + f["applies_when"].get("text", ""))
    elif f.get("applies_when"):
        parts.append("only when " + f["applies_when"].get("text", ""))
    if f.get("by_quantity"):
        parts.append("where the quantity takes it")
    if f.get("repeatable"):
        parts.append("repeatable (a list)")
    if f.get("default"):
        parts.append("default " + str(f["default"]))
    if f.get("units"):
        parts.append("units " + ", ".join(f["units"]))
    if f.get("refers_to") or f.get("or_refers_to"):
        parts.append("names a " + " or ".join(f.get("refers_to") or f.get("or_refers_to")))
    text = f"    {key}: " + "; ".join(parts) + "."
    if f.get("doc"):
        text += "\n        " + f["doc"]
    for a in f.get("alternatives", []):
        if a.get("keyword") and a.get("doc"):
            text += f"\n        {a['keyword']}: {a['doc']}"
    return text


def _make_add(kind, spec):
    def add(self, name, /, **fields):
        return self.add(kind, name, **fields)

    params = [inspect.Parameter("self", inspect.Parameter.POSITIONAL_ONLY),
              inspect.Parameter("name", inspect.Parameter.POSITIONAL_ONLY)]
    fields = spec["fields"]
    params += [inspect.Parameter(k, inspect.Parameter.KEYWORD_ONLY) for k, f in fields.items() if _required_always(f)]
    params += [inspect.Parameter(k, inspect.Parameter.KEYWORD_ONLY, default=None)
               for k, f in fields.items() if not _required_always(f)]
    add.__signature__ = inspect.Signature(params)
    add.__name__ = add.__qualname__ = f"add_{kind}"
    sections = [s for s, kinds in _SCHEMA["sections"].items() if kind in kinds]
    add.__doc__ = (f"Add a {kind} record ({spec.get('title', kind)}) to <{sections[0] if sections else '?'}>: "
                   f"Project.add({kind!r}, name, **fields).\n\n{spec.get('doc', '')}\n\n"
                   f"Fields (the engine's table, pyazr.schema()['records'][{kind!r}]; a value is the file's text, a "
                   f"number, a (value, unit) tuple, a pyazr.Nuisance, or a list for a repeatable field):\n"
                   + "\n".join(_field_doc(k, f) for k, f in fields.items())
                   + "\n\nRaises FormatError (the engine's diagnostics) when the record is refused.  Returns the Record.")
    return add


for _kind, _spec in _SCHEMA["records"].items():
    setattr(Project, f"add_{_kind}", _make_add(_kind, _spec))
