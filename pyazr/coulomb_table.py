"""The Coulomb table: where the engine looks for it, how to build one here, how to copy one you name.

The engine interpolates the Coulomb functions in a table (``CoulombTable``, 85 MB) when it has one and uses its exact
routines otherwise; the numbers agree within the table's accuracy (1e-7), the table is about four times faster.  A
session and the program read the file named by the option ``coulomb_table``, else ``$AZURE2_COULOMB_TABLE``, else
``~/.azure2/coulomb_table.bin`` (:func:`where`).

The table is not published (the owner's decision of 2026-11-18): every user builds it once, locally, with
``python -m pyazr coulomb-table build`` (or ``AZURE2 --build-coulomb-table``): about 5 minutes on 4 threads, 88 MB.
:func:`build` writes the default file, so every later session and run uses the table with nothing set.
:func:`fetch` copies a table from a URL you name, accepted only with the SHA-256 you give (a table built once on
a shared disk or a web server of your own); nothing is downloaded from anywhere by default.
Command line: ``pyazr coulomb-table where|build|fetch|check``.
"""

import hashlib
import os
import urllib.request

_BLOCK = 1 << 20


class ChecksumMismatch(OSError):
    """A copied table whose SHA-256 is not the one given (a damaged or different file); it was removed."""


def _engine():
    from . import _azure2
    return _azure2


def default_path():
    """The file a session and the program read when no table is named: ``$AZURE2_COULOMB_TABLE``, else
    ``~/.azure2/coulomb_table.bin`` (the engine's ``CoulombTable::DefaultPath``)."""
    return _engine().coulomb_table_default_path()


def file_sha256(path):
    """The SHA-256 of a file (hex)."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(_BLOCK), b""):
            h.update(block)
    return h.hexdigest()


def where(path=""):
    """Where the table is (``path``, else :func:`default_path`) and whether a file is there, without loading it:
    ``{"path", "exists", "bytes", "source", "url", "sha256"}``; ``source`` says where the path comes from (the
    argument, ``AZURE2_COULOMB_TABLE`` or the default).  ``url`` and ``sha256`` are always empty (no table is
    published, so there is no address to give); the keys stay so that a script reading them keeps working."""
    if path:
        source = "argument"
    elif os.environ.get("AZURE2_COULOMB_TABLE"):
        source = "AZURE2_COULOMB_TABLE"
    else:
        source = "default"
    file = path or default_path()
    exists = os.path.isfile(file)
    return {"path": file, "exists": exists, "bytes": os.path.getsize(file) if exists else 0, "source": source,
            "url": "", "sha256": ""}


def fetch(path="", *, url=None, sha256=None, overwrite=False):
    """Copy a table you name to ``path`` (else :func:`default_path`) and check it.  No table is published: to get one
    build it here (:func:`build`, ``python -m pyazr coulomb-table build``, about 5 minutes).

    ``url`` (required; ``file://`` and ``http(s)://``) is read into ``path + ".part"``; its SHA-256 must equal
    ``sha256`` (required), else the partial file is removed and :class:`ChecksumMismatch` (an ``OSError``) is raised;
    then it is renamed over ``path`` (a process that has the old file open keeps it).  An existing ``path`` is kept
    unless ``overwrite``.  Returns :func:`where` of the file, with ``fetched`` (True, or False if an existing file was
    kept).  Raises ``ValueError`` when no URL or no hash is given (build the table instead: :func:`build`),
    ``OSError`` when the download fails."""
    if not url or not sha256:
        raise ValueError("the Coulomb table is not published: build it here with pyazr.coulomb_table.build() "
                         "(python -m pyazr coulomb-table build, about 5 minutes), or give the URL and the SHA-256 "
                         "of a copy you have")
    file = path or default_path()
    if os.path.exists(file) and not overwrite:
        return dict(where(path), fetched=False)
    folder = os.path.dirname(os.path.abspath(file))
    os.makedirs(folder, exist_ok=True)
    part = file + ".part"
    h = hashlib.sha256()
    try:
        with urllib.request.urlopen(url) as response, open(part, "wb") as out:
            for block in iter(lambda: response.read(_BLOCK), b""):
                h.update(block)
                out.write(block)
        if h.hexdigest() != sha256.lower():
            raise ChecksumMismatch(f"the file at {url} has SHA-256 {h.hexdigest()}, not the given "
                                   f"{sha256.lower()}: refused and removed (a damaged or different file); try "
                                   "again, or build the table here (pyazr coulomb-table build)")
        os.replace(part, file)
    finally:
        if os.path.exists(part):
            os.remove(part)
    return dict(where(path), fetched=True)


def build(path="", *, threads=None, l_max=None, per_decade=None):
    """Build the table from the engine's exact routines and write it to ``path`` (else :func:`default_path`), as
    ``AZURE2 --build-coulomb-table`` does; ``threads`` (default: the engine's OpenMP threads), ``l_max`` and
    ``per_decade`` (default: the program's table; smaller values give a smaller, coarser table, for tests).  Returns
    :func:`where` of the file."""
    engine = _engine()
    options = {"threads": int(threads if threads is not None else engine.build_info()["max_threads"])}
    if l_max is not None:
        options["l_max"] = int(l_max)
    if per_decade is not None:
        options["per_decade"] = float(per_decade)
    return where(engine.build_coulomb_table(path or "", **options))
