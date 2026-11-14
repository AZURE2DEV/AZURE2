#!/usr/bin/env python3
"""The published pages build as the workflow builds them (.github/workflows/docs.yml): ``sphinx-build -W
--keep-going -b html docs/source``, where a warning is an error (a broken reference, a page outside every toctree,
text that reStructuredText reads as markup); and nitpicky (``-n``): a ``:py:`` reference (``:class:``, ``:meth:``,
...) needs a target that exists.  The pages have no API reference of pyazr, so a name of pyazr is written as a
literal (``pyazr.Project``).  Without Sphinx and the two packages of docs/requirements.txt: prints
``skip: ...`` and exits 77.
"""
import importlib.util
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NEEDED = ("sphinx", "sphinx_rtd_theme", "sphinx_copybutton")


def main():
    missing = [m for m in NEEDED if importlib.util.find_spec(m) is None]
    if missing:
        print(f"skip: {', '.join(missing)} not installed (pip install -r docs/requirements.txt)")
        return 77
    with tempfile.TemporaryDirectory() as out:
        run = subprocess.run([sys.executable, "-m", "sphinx", "-W", "-n", "--keep-going", "-q", "-b", "html",
                              "-d", os.path.join(out, "doctrees"), os.path.join(ROOT, "docs", "source"),
                              os.path.join(out, "html")], capture_output=True, text=True)
        pages = sum(name.endswith(".html") for _, _, names in os.walk(os.path.join(out, "html")) for name in names)
    if run.returncode != 0:
        print(run.stdout + run.stderr)
        print("FAIL: the published pages do not build without warnings")
        return 1
    print(f"the published pages build without warnings: {pages} pages")
    return 0


if __name__ == "__main__":
    sys.exit(main())
