"""The state of the web UI's server (docs/upgrade/web_ui.md section 5): the project text with its history.

**The project text is the state.**  Every change of the model is an edit through :class:`pyazr.Project`, validated
by the engine on the project as it would stand; a refused edit leaves the text and the revision as they were.  Each
accepted edit, undo or redo gives a new revision number; an edit carries the revision it was made on and is refused
when that is no longer the current one (two clients).  The history keeps the texts, so undo and redo restore a text
byte for byte.

Files, all inside the project's directory (``<name>`` the project file's name without ``.azr``):

- the user's file: written only by :meth:`State.save`;
- the working file ``.<name>.ui.azr`` beside it: the current text, rewritten on every change while it differs from
  the saved text; it is the crash copy a restarted server offers back (:meth:`State.restore`);
- the session files ``.<name>.ui-session.azr`` (the text the session worker reads) and ``.<name>.ui-job-<id>.azr``
  (a job's) beside it: the engine resolves some relative paths (the prior covariance) against the directory of the
  file it reads, so they cannot live elsewhere;
- the UI's directory ``.<name>.ui/``: ``engine/`` (the output directory of the UI's sessions), ``jobs/`` (job
  checkpoints and logs), ``snapshots/``, ``logs/``, ``url``.

**Session options** (Brune, Park, transform, the Coulomb table, ...): the options of the calculation are keys of the
project's ``<config>`` (``docs/upgrade/project_format3.md`` section 3.12), so the text is the whole state: the
server opens its sessions without keywords, and :func:`session_options` is the engine's table of options
(``pyazr.kinds()["options"]``) with the value in use and its ``source``.  Setting one is the config route's edit
(``Project.set_config``: one history entry, a new session).
"""
import difflib
import os
import re
import threading
import time
from collections import Counter

from .errors import Conflict, NotFound, OutsideRoot, RevisionRequired, StaleRevision

def session_options(project=None, session=None):
    """The engine's table of session options (``pyazr.kinds()["options"]``), each row with ``value`` (the one in use)
    and ``source``: from ``session`` (``Session.describe()``: its ``options`` and ``option_sources``) when given, else
    from ``project``'s ``<config>`` (``"project"`` for a key it states, ``"default"`` otherwise).  An option that is
    not project state (``extrapolations``, ``coulomb_table``, ``cwd``) has its default and source ``"default"``."""
    import pyazr
    config = project.config() if project is not None and session is None else {}
    out = {}
    for keyword, row in pyazr.kinds()["options"].items():
        value, source = row["default"], "default"
        if session is not None and keyword in session.get("option_sources", {}):
            value, source = session["options"][keyword], session["option_sources"][keyword]
        elif row["key"] and row["key"] in config:
            value, source = config[row["key"]] in ("1", "true"), "project"
        out[keyword] = dict(row, value=value, source=source)
    return out


_ENCODING = ("utf-8", "surrogateescape")


class Entry:
    """One history entry: the revision, the time, what was done, the text after it, and the files of the project's
    directory the entry wrote beside the text (``files``: path -> (bytes before or None when absent, bytes after)),
    which undo puts back as they were and redo writes again."""

    def __init__(self, rev, what, text, files=None):
        self.rev, self.what, self.text, self.time = rev, what, text, time.time()
        self.files = files or {}


class State:
    """The open project: its :class:`pyazr.Project`, history, revision, files and root.

    ``root``: the directory every path read or written must resolve into (default the project's directory)."""

    def __init__(self, path, root=None, *, project=None):
        from pyazr import Project
        path = os.path.abspath(path)
        self.root = os.path.realpath(root or os.path.dirname(path))
        path = self.inside_root(path)
        self.project = project if project is not None else Project.open(path)
        self.saved_text = None if project is not None else self.project.text
        self.lock = threading.RLock()
        self.rev = 1
        self.history = [Entry(1, "open" if project is None else "new", self.project.text)]
        self.pos = 0
        self.crash_text = None
        if os.path.isfile(self.working_file):
            with open(self.working_file, "rb") as f:
                text = f.read().decode(*_ENCODING)
            if text != self.project.text:
                self.crash_text = text
        os.makedirs(self.ui_dir, exist_ok=True)

    # -- files ----------------------------------------------------------------------------------------------------

    @property
    def path(self):
        return self.project.path

    @property
    def directory(self):
        return os.path.dirname(self.path)

    @property
    def stem(self):
        return os.path.splitext(os.path.basename(self.path))[0]

    @property
    def working_file(self):
        return os.path.join(self.directory, f".{self.stem}.ui.azr")

    @property
    def ui_dir(self):
        return os.path.join(self.directory, f".{self.stem}.ui")

    def session_file(self, job=None):
        """The file a session of the UI reads: the session worker's, or job ``job``'s (:func:`project_file_of` gives
        the project's file name back)."""
        return os.path.join(self.directory, f".{self.stem}.ui-{'session' if job is None else 'job-' + job}.azr")

    def ui_path(self, *parts):
        p = os.path.join(self.ui_dir, *parts)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        return p

    def inside_root(self, path):
        """``path`` (relative to the project's directory when relative) resolved, symbolic links followed; refused
        (403) when it is outside the root."""
        base = self.directory if hasattr(self, "project") else self.root
        full = os.path.realpath(os.path.join(base, str(path)))
        if os.path.commonpath([full, self.root]) != self.root:
            raise OutsideRoot(f"{path}: outside the directory the server may read and write ({self.root})",
                              allowed=f"paths inside {self.root}", root=self.root,
                              fix=f"name a file inside {self.root}, relative to the project's folder {base}")
        return full

    # -- reading ----------------------------------------------------------------------------------------------------

    @property
    def unsaved(self):
        return self.project.text != self.saved_text

    def describe(self):
        with self.lock:
            crash = None if self.crash_text is None else {"path": self.working_file,
                                                          "modified": os.path.getmtime(self.working_file)}
            return {"path": self.path, "root": self.root, "rev": self.rev, "unsaved": self.unsaved,
                    "can_undo": self.pos > 0, "can_redo": self.pos < len(self.history) - 1,
                    "working_file": self.working_file, "ui_dir": self.ui_dir, "crash_copy": crash}

    def history_view(self):
        with self.lock:
            return [{"rev": e.rev, "time": e.time, "what": e.what, "current": i == self.pos}
                    for i, e in enumerate(self.history)]

    # -- edits ----------------------------------------------------------------------------------------------------

    def check_rev(self, rev):
        if rev is None:
            raise RevisionRequired(f"an edit names the revision it was made on (the current one is {self.rev})",
                                   fix="send rev (in the body, the query or an If-Match header)", rev=self.rev)
        if int(rev) != self.rev:
            raise StaleRevision(f"the edit was made on revision {rev}, the project is at revision {self.rev}",
                                fix="GET /api/project, then redo the edit on the current text", rev=self.rev)

    def edit(self, rev, what, change):
        """Apply ``change(project)`` (an edit of :class:`pyazr.Project`, which validates it and leaves the project
        unchanged when it refuses) made on revision ``rev``; returns the new revision and the change's result."""
        with self.lock:
            self.check_rev(rev)
            before = self.project.text
            result = change(self.project)
            if self.project.text != before:
                self._push(what)
            return result

    def edit_with_files(self, rev, what, change, files):
        """An edit that writes files beside the text (a covariance file the edit names): ``files`` {path inside the
        root: text} are written first (the engine validates the edit with them), then ``change(project)`` is made;
        a refusal puts the files back as they were and leaves the text.  One history entry when the text or a file
        changes (a file rewritten under the name the text already holds is an entry too): undo restores the files
        as they were before (removing those that did not exist) with the text, redo writes them again."""
        with self.lock:
            self.check_rev(rev)
            written = {}
            for path, text in files.items():
                full = self.inside_root(path)
                before = None
                if os.path.isfile(full):
                    with open(full, "rb") as f:
                        before = f.read()
                written[full] = (before, text.encode(*_ENCODING))
            before_text = self.project.text
            _put_files({p: after for p, (_, after) in written.items()})
            try:
                result = change(self.project)
            except BaseException:
                _put_files({p: b for p, (b, _) in written.items()})
                raise
            if self.project.text == before_text and all(b == a for b, a in written.values()):
                return result
            self._push(what, written)
            return result

    def replace_text(self, rev, text, what="replace text"):
        """A whole-text edit: refused when the engine finds an error the current text does not have, or when a path
        named in the text resolves outside the root (checked before the engine reads any file)."""
        from pyazr import Project, FormatError
        with self.lock:
            self.check_rev(rev)
            candidate = Project(text.encode(*_ENCODING), self.path)
            self.check_paths_of(candidate)
            (before, skipped_before), (after, skipped_after) = _diagnostics(self.project), _diagnostics(candidate)
            new = _new_errors(before, after)
            if new:
                raise FormatError(f"{what}: the engine refuses the text", new)
            if skipped_after and not skipped_before:
                raise FormatError(f"{what}: the engine cannot check the text: {skipped_after}")
            candidate.levels(), candidate.config(), candidate.thm()  # the engine's reading of the nuclear part
            self.project = candidate
            self._push(what)

    def _push(self, what, files=None):
        del self.history[self.pos + 1:]
        self.rev += 1
        self.history.append(Entry(self.rev, what, self.project.text, files))
        self.pos = len(self.history) - 1
        self._write_working_file()

    def _write_working_file(self):
        if self.unsaved:
            tmp = self.working_file + ".tmp"
            with open(tmp, "wb") as f:
                f.write(self.project.text.encode(*_ENCODING))
            os.replace(tmp, self.working_file)
            self.crash_text = None
        elif os.path.exists(self.working_file) and self.crash_text is None:
            os.remove(self.working_file)

    def _move(self, step, rev):
        from pyazr import Project
        with self.lock:
            self.check_rev(rev)
            pos = self.pos + step
            if not 0 <= pos < len(self.history):
                raise Conflict(f"nothing to {'undo' if step < 0 else 'redo'}", rev=self.rev)
            # the files of the entry undone go back to how they were before it; those of the entry redone are written
            if step < 0:
                _put_files({p: before for p, (before, _) in self.history[self.pos].files.items()})
            else:
                _put_files({p: after for p, (_, after) in self.history[pos].files.items()})
            self.pos = pos
            self.project = Project(self.history[pos].text.encode(*_ENCODING), self.path)
            self.rev += 1
            self._write_working_file()
            return self.rev

    def undo(self, rev):
        return self._move(-1, rev)

    def redo(self, rev):
        return self._move(1, rev)

    def save(self, rev, path=None):
        """Write the project to its file, or to ``path`` inside the root (which becomes the project's file).  The
        only call that writes the user's file."""
        with self.lock:
            self.check_rev(rev)
            target = None if path is None else self.inside_root(path)
            written = self.project.save(target)
            self.saved_text = self.project.text
            self._write_working_file()
            return written

    def restore(self, rev):
        """Take the crash copy found at start-up as an edit (validated like any text)."""
        with self.lock:
            if self.crash_text is None:
                raise NotFound("no crash copy to restore", rev=self.rev)
            text, self.crash_text = self.crash_text, None
            self.replace_text(rev, text, "restore the working copy")
            return self.rev

    # -- paths named by an edit -------------------------------------------------------------------------------------

    def check_paths(self, values):
        """Refuse (403) a text value of an edit any of whose tokens resolves outside the root: a data file, a
        covariance, a stopping table, a prior covariance, the output directory.  Every token of every value is taken
        as a path relative to the project's directory: a token that is not a path (a number, a unit, a keyword) stays
        inside, only an absolute path, a ``..`` or a symbolic link leading out leaves it."""
        for value in _flatten(values):
            for token in re.split(r"[\s,;=:]+", str(value)):
                if re.search(r"\w", token):  # an operator alone ("/" of a formula) names no file
                    self.inside_root(token)

    def check_paths_of(self, project):
        """:meth:`check_paths` on every record field, ``<config>`` entry and ``<thm>`` value of ``project``, parsed
        without reading any file."""
        values = [r.fields for r in project.records().values()]
        values.append(project.config())
        values.append(project.thm() or {})
        self.check_paths(values)

    # -- the session's text -----------------------------------------------------------------------------------------

    def session_text(self, text=None):
        """The text a session of the UI reads: the project's text with ``<config>``'s output directory pointed at the
        UI's ``engine/`` directory, so that sessions never write into the user's output directory."""
        from pyazr import Project
        p = Project((text or self.project.text).encode(*_ENCODING), self.path)
        engine_dir = os.path.relpath(self.ui_path("engine", "x"), self.directory)
        p.set_config(output_dir=os.path.dirname(engine_dir) + "/")
        return p.text

    def merge_snapshot(self, text, session_text, snapshot_text):
        """The project text ``text`` with the changes a snapshot made to the session's text ``session_text`` (a
        three-way merge by lines): the lines where the session's text differs from the project's (the output
        directory) are the project's; every line the snapshot changed is the snapshot's."""
        base = text.splitlines(keepends=True)
        ours = session_text.splitlines(keepends=True)
        theirs = snapshot_text.splitlines(keepends=True)
        if len(base) != len(ours):
            raise RuntimeError("the session's text does not have the project's lines")
        redirected = {i for i in range(len(base)) if base[i] != ours[i]}
        out = []
        for op, i1, i2, j1, j2 in difflib.SequenceMatcher(None, ours, theirs, autojunk=False).get_opcodes():
            if op == "equal":
                out.extend(base[i1:i2])
                continue
            if redirected & set(range(i1, i2)):
                raise Conflict("the snapshot changed a line of <config> the UI redirects; not merged")
            out.extend(theirs[j1:j2])
        return "".join(out)


def _put_files(contents):
    """Write each file of ``contents`` {path: bytes}, or remove it where the bytes are None."""
    for path, data in contents.items():
        if data is None:
            if os.path.exists(path):
                os.remove(path)
            continue
        tmp = path + ".tmp"
        with open(tmp, "wb") as f:
            f.write(data)
        os.replace(tmp, path)


def _flatten(values):
    if isinstance(values, dict):
        for v in values.values():
            yield from _flatten(v)
    elif isinstance(values, (list, tuple)):
        for v in values:
            yield from _flatten(v)
    elif values is not None:
        yield values


def _diagnostics(project):
    """(the validator's diagnostics, the warning that build checks were skipped or None) of ``project``."""
    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        diagnostics = project.validate()
    skipped = [str(w.message) for w in caught if "not run" in str(w.message)]
    return diagnostics, (skipped[0] if skipped else None)


def _new_errors(before, after):
    """The errors of ``after`` that ``before`` does not have, by (code, record, field) and count."""
    key = lambda d: (d.get("code"), d.get("record"), d.get("field"))
    errors = lambda ds: [d for d in ds if d.get("severity") == "error"]
    extra = Counter(key(d) for d in errors(after)) - Counter(key(d) for d in errors(before))
    out = []
    for d in errors(after):
        if extra[key(d)] > 0:
            extra[key(d)] -= 1
            out.append(d)
    return out


def project_file_of(session_file):
    """The project's file name (``T1.azr``) of a file a session of the UI reads (``.T1.ui-session.azr``,
    ``.T1.ui-job-<id>.azr``: :meth:`State.session_file`); ``None`` for any other name."""
    m = re.match(r"^\.(.+)\.ui-(session|job-[^/]+)\.azr$", os.path.basename(session_file))
    return f"{m.group(1)}.azr" if m else None
