"""The Coulomb table: where the engine looks for it, how to fetch a published copy, how to build one here.

The engine interpolates the Coulomb functions in a table (``CoulombTable``, 85 MB) when it has one and uses its exact
routines otherwise; the numbers agree within the table's accuracy (1e-7), the table is about four times faster.  A
session and the program read the file named by the option ``coulomb_table``, else ``$AZURE2_COULOMB_TABLE``, else
``~/.azure2/coulomb_table.bin`` (:func:`where`).  :func:`fetch` and :func:`build` write that file by default, so every
later session and run uses the table with nothing set.  Command line: ``pyazr coulomb-table where|fetch|build|check``.

:func:`fetch` downloads :data:`TABLE_URL` and accepts the file only if its SHA-256 is :data:`TABLE_SHA256`; both
are empty until the owner publishes the file (``docs/upgrade/FOR_THE_OWNER.md``, section 8), and :func:`build` is
the way until then (and offline): about 5 minutes on 4 threads.
"""

import hashlib
import os
import urllib.request

#: Where the published table is downloaded from: set by the owner when he publishes the file (a release asset).
TABLE_URL = ""
#: The SHA-256 of the published file (hex): set together with :data:`TABLE_URL`.
TABLE_SHA256 = ""

_BLOCK = 1 << 20


class ChecksumMismatch(OSError):
    """A downloaded table whose SHA-256 is not the published one (a damaged or different file); it was removed."""


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
    argument, ``AZURE2_COULOMB_TABLE`` or the default), ``url`` and ``sha256`` are those :func:`fetch` uses."""
    if path:
        source = "argument"
    elif os.environ.get("AZURE2_COULOMB_TABLE"):
        source = "AZURE2_COULOMB_TABLE"
    else:
        source = "default"
    file = path or default_path()
    exists = os.path.isfile(file)
    return {"path": file, "exists": exists, "bytes": os.path.getsize(file) if exists else 0, "source": source,
            "url": TABLE_URL, "sha256": TABLE_SHA256}


def fetch(path="", *, url=None, sha256=None, overwrite=False):
    """Download the table to ``path`` (else :func:`default_path`) and check it.

    ``url`` (default :data:`TABLE_URL`; ``file://`` and ``http(s)://``) is read into ``path + ".part"``; its SHA-256
    must equal ``sha256`` (default :data:`TABLE_SHA256`), else the partial file is removed and
    :class:`ChecksumMismatch` (an ``OSError``) is raised; then it is renamed over ``path`` (a process that has the old
    file open keeps it).  An existing ``path`` is kept unless ``overwrite``.  Returns :func:`where` of the file, with
    ``fetched`` (True, or False if an existing file was kept).  Raises ``ValueError`` when no URL or no hash is known
    (build the table instead: :func:`build`), ``OSError`` when the download fails."""
    url = TABLE_URL if url is None else url
    sha256 = TABLE_SHA256 if sha256 is None else sha256
    if not url or not sha256:
        raise ValueError("no published Coulomb table is known to this pyazr (TABLE_URL and TABLE_SHA256 of "
                         "pyazr.coulomb_table are empty): build it here with pyazr.coulomb_table.build() "
                         "(pyazr coulomb-table build), or give the URL and SHA-256 of a copy")
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
            raise ChecksumMismatch(f"the file at {url} has SHA-256 {h.hexdigest()}, not the published "
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
