"""``python -m pyazr COMMAND ...``: the pyazr command line (:mod:`pyazr.cli`; ``python -m pyazr commands``)."""
import sys

from .cli import main

sys.exit(main())
