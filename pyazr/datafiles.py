"""The files under ``pyazr/`` that are not Python: which the package reads at run time (and so the wheel carries),
which it does not (and why).

The wheel carries the package directory as it is (``pyproject.toml``: ``wheel.packages``) less the patterns of
:data:`NOT_SHIPPED` (``wheel.exclude``, kept equal by consistency check 48).  The readers find their files next to
their own module (``__file__``), so an installed or unpacked package needs no environment variable.  A new run-time
data file is one file under ``pyazr/`` and one row of :data:`DATA_FILES`; the wheel's test
(``tests/pyazr/wheel_test.py``) and check 48 fail until the row is there.

Paths are relative to the package directory, with ``/``; a path ending in ``/`` stands for every file below it.
"""

import fnmatch
import os

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

#: (path, what it is and who reads it) of every non-Python file the wheel carries.
DATA_FILES = (
    ("data/mass_1.mas20", "the AME2020 mass table; pyazr.masses"),
    ("data/README.md", "where the mass table comes from and its terms (provenance, read by people)"),
    ("ui/static/", "the web page: HTML, scripts, style and the vendored plotly with its licence; pyazr.ui.routes"),
)

#: (pattern of pyproject.toml's wheel.exclude, relative to the repository, why it is not shipped).  The engine module
#: is not a file of the source: the wheel's comes from its CMake install rule, and a wheel build refuses to start
#: while an in-tree build's module lies in pyazr/ (api/src/CMakeLists.txt).
NOT_SHIPPED = (
    ("pyazr/CLAUDE.md", "guidance for agents that edit the package"),
    ("pyazr/examples/**", "worked examples that read the repository's test projects (tests/)"),
    ("**/__pycache__/**", "bytecode of this machine"),
)


def path(name):
    """The absolute path of the data file ``name`` (a path of :data:`DATA_FILES`)."""
    return os.path.join(PACKAGE_DIR, *name.rstrip("/").split("/"))


def shipped(rel):
    """True if ``rel`` (a path relative to the package directory, with ``/``) is a file :data:`DATA_FILES` names."""
    return any(rel == p or (p.endswith("/") and rel.startswith(p)) for p, _ in DATA_FILES)


def excluded(rel):
    """True if ``rel`` (relative to the package directory) matches a pattern of :data:`NOT_SHIPPED` (``*`` matches
    across ``/`` here, which is what the patterns mean in ``wheel.exclude``)."""
    return any(fnmatch.fnmatchcase("pyazr/" + rel, pattern) for pattern, _ in NOT_SHIPPED)


def package_files(root=PACKAGE_DIR):
    """The files of the package directory ``root`` that the wheel carries, relative to it, sorted: every ``.py`` file
    and every file :data:`DATA_FILES` names, less :data:`NOT_SHIPPED`.  The wheel's expected list is these, the
    engine module, the program and the wheel's metadata."""
    out = []
    for folder, dirs, files in os.walk(root):
        dirs.sort()
        for name in files:
            rel = os.path.relpath(os.path.join(folder, name), root).replace(os.sep, "/")
            if excluded(rel) or rel.startswith("_azure2"):
                continue
            if rel.endswith(".py") or shipped(rel):
                out.append(rel)
    return sorted(out)
