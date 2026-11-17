"""One engine session in one process, driven over a pipe: the worker of the backend in use.

``_worker_p2`` holds a ``pyazr.azure2`` session (pyazr 2), ``_worker_p3`` a ``pyazr.Session`` (pyazr 3); both answer
the commands of :class:`rmfit.engine._Worker`.  The pyazr 3 worker also answers ``model_rows`` and ``save_fit``; the
pyazr 2 worker ``scatlen`` / ``scatlen_jac`` (see :mod:`rmfit.backend`).
"""
from .backend import PYAZR3

if PYAZR3:
    from ._worker_p3 import worker_main  # noqa: F401
else:
    from ._worker_p2 import worker_main  # noqa: F401
