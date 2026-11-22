"""rmfit -- agent-operable global R-matrix fitting on top of AZURE2 / pyazr.

See README.md for the layout.  Every module imports pyazr from ``AZURE2_ROOT``.
"""
import os
import sys


def _containing_azure2():
    """The AZURE2 checkout this package sits in (rmfit/ at the top of an AZURE2 tree), else None."""
    parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if os.path.isfile(os.path.join(parent, "pyazr", "__init__.py")) and os.path.isfile(os.path.join(parent, "src", "AZURE2.cpp")):
        return parent
    return None


#: The AZURE2 tree pyazr is imported from: $AZURE2_ROOT, else the AZURE2 checkout rmfit sits in, else ~/AZURE2.
AZURE2_ROOT = os.environ.get("AZURE2_ROOT") or _containing_azure2() or os.path.expanduser("~/AZURE2")
# pyazr must come from AZURE2_ROOT even when another directory on the path also has a pyazr package -- the directory
# rmfit was imported from, when that is a different AZURE2 checkout than the campaign's pinned engine.  rmfit itself is
# already imported, so moving AZURE2_ROOT to the front changes only where pyazr (and nothing of rmfit's) is found.
_root = os.path.abspath(AZURE2_ROOT)
sys.path[:] = [p for p in sys.path if not p or os.path.abspath(p) != _root]
sys.path.insert(0, AZURE2_ROOT)
# ...but a spawned worker (multiprocessing "spawn": Sharded shards, SingleSession) starts from THIS process's sys.path and
# imports rmfit afresh, so with AZURE2_ROOT first it finds AZURE2_ROOT/rmfit when that tree has one -- every clone of the
# upgrade branch tracks an old rmfit/ (workers on such a frozen engine copy ran it).
# Then import pyazr from AZURE2_ROOT now and put rmfit's own parent back in front: pyazr is bound to AZURE2_ROOT's copy
# from here on, and the worker finds this rmfit (whose import repeats the same two steps).
_parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.isdir(os.path.join(_root, "rmfit")) and _root != _parent:
    try:
        import pyazr  # noqa: F401
    except Exception:  # noqa: BLE001 -- not importable here: keep the order; rmfit.backend reports the error
        pass
    else:
        sys.path[:] = [p for p in sys.path if not p or os.path.abspath(p) != _parent]
        sys.path.insert(0, _parent)

#: Runtime options the CLI fits of the development archive used
#: (``--no-long-wavelength --ignore-externals --use-brune``).  pyazr's own
#: defaults differ (``use_long_wavelength=True``), so pass these explicitly.
ARCHIVE_FLAGS = dict(use_long_wavelength=False, ignore_externals=True, use_brune=True)

__version__ = "0.1.0"
