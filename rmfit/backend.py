"""Which pyazr rmfit runs on: the one backend switch.

rmfit is one code base on two engine interfaces (the evaluator 2026-11-16: "One rmfit, two backends"):

* **pyazr 3** (AZURE2 ``dev`` since 2026-11-25, when feature/upgrade was merged; the engine of new campaigns): project
  format 3 (records), ``pyazr.Session`` / ``Objective`` / ``Evaluation`` / ``Parameters`` / ``Project``; rmfit's only
  contact with a project is :mod:`rmfit.azr3`.
* **pyazr 2** (the archive's pinned 2.x builds, which older campaigns stay on): the project is the ``.azr`` text format,
  the engine session is ``pyazr.azure2``, the model editor ``pyazr.AzrModel``, the data-set table
  ``pyazr.datasets.SegmentSet``.

Everything that does not touch the engine or the project file (search, ledger, driver logic, review criteria, keys)
is shared.  The contact points branch on :data:`PYAZR3`; the large contact layers have one module per backend
(``_worker_p2`` / ``_worker_p3``, ``_model_p2`` / ``_model_p3``) behind the module the rest of rmfit imports.

The version is read once, from the pyazr that ``rmfit.AZURE2_ROOT`` puts first on the path.
"""
from . import AZURE2_ROOT  # noqa: F401  (puts AZURE2_ROOT's pyazr first on sys.path)


class BackendError(ImportError):
    """pyazr could not be imported, or its version is neither 2 nor 3."""


def _detect():
    try:
        import pyazr
    except Exception as err:  # pragma: no cover - environment problem
        raise BackendError(f"rmfit: cannot import pyazr from AZURE2_ROOT={AZURE2_ROOT!r}: {err}") from err
    ver = str(getattr(pyazr, "__version__", "") or "")
    try:
        major = int(ver.split(".")[0])
    except ValueError:
        major = 3 if hasattr(pyazr, "Session") else (2 if hasattr(pyazr, "AzrModel") else 0)
    if major not in (2, 3):
        raise BackendError(f"rmfit: pyazr {ver or '(no version)'} at {getattr(pyazr, '__file__', '?')} is neither "
                           f"pyazr 2 nor pyazr 3")
    return major, ver


PYAZR_MAJOR, PYAZR_VERSION = _detect()
#: True on pyazr 3 (project format 3), False on pyazr 2 (.azr text format)
PYAZR3 = PYAZR_MAJOR == 3

# campaign.json "flags" keys that the side sessions (init's .sav check, uq) pass on to the engine.  The coulomb_table pair
# was missing until 2026-11-26: a large campaign's init opened its two .sav-check sessions on the exact Coulomb functions
# (the default table path held a format-1 table) and invalidated output/intEC.dat, about a fifth of its init time.
ENGINE_FLAG_KEYS = ("use_brune", "use_park", "ignore_externals", "use_long_wavelength", "use_coulomb_table", "coulomb_table")


def name() -> str:
    """``"pyazr 3.0.0.dev3"`` style: the backend in use, for logs and reports."""
    return f"pyazr {PYAZR_VERSION or PYAZR_MAJOR}"


def require(major: int, what: str):
    """Raise a clear error when ``what`` is not available on the pyazr in use."""
    if PYAZR_MAJOR != major:
        raise NotImplementedError(f"{what} needs pyazr {major}; this process runs {name()}")


def __getattr__(attr):
    """The backend's classes by the names rmfit's shared code uses (imported on first use):

    * ``AzrModel``: the model editor -- ``pyazr.AzrModel`` (pyazr 2) or :class:`rmfit.azr3.Model` (pyazr 3);
    * ``SegmentSet``: the data-set table -- ``pyazr.datasets.SegmentSet`` or :class:`rmfit.azr3.SegmentTable`;
    * ``Nucleus``: ``pyazr.nuclei.Nucleus`` or ``pyazr.system.Nucleus``.
    """
    if attr == "AzrModel":
        if PYAZR3:
            from .azr3 import Model
            return Model
        from pyazr import AzrModel
        return AzrModel
    if attr == "SegmentSet":
        if PYAZR3:
            from .azr3 import SegmentTable
            return SegmentTable
        from pyazr.datasets import SegmentSet
        return SegmentSet
    if attr == "Nucleus":
        if PYAZR3:
            from pyazr.system import Nucleus
            return Nucleus
        from pyazr.nuclei import Nucleus
        return Nucleus
    raise AttributeError(f"module {__name__!r} has no attribute {attr!r}")
