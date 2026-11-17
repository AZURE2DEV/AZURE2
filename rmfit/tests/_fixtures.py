"""The fixture directory of the backend in use: ``fixtures/`` (pyazr 2, .azr text format) or ``fixtures_p3/`` (pyazr 3,
the same models converted to project format 3 with ``pyazr.convert_azr``; data files as in ``fixtures/``)."""
import os

from rmfit.backend import PYAZR3

DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures_p3" if PYAZR3 else "fixtures")


def fixture(*parts):
    return os.path.join(DIR, *parts)
