"""The old projects of the test tree for the tests of pyazr's old file class (pyazr/azrfile.py, pyazr/datasets.py).

The tree's projects are format 3 (docs/upgrade/format3_tree.md); each original old file is kept in
tests/convert/fixtures/<path from the repository root>, without data files.  A test of the old reader reads the
original instead of the tree's file:

    original(path)        the original of a tree file (absolute path), or the file itself when it was not converted
    copytree(src, dst)    shutil.copytree of a tree folder, with the originals of its converted projects laid over

They go with the old reader.
"""
import os
import shutil

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
FIXTURES = os.path.join(ROOT, "tests", "convert", "fixtures")


def _fixture(path):
    rel = os.path.relpath(os.path.abspath(path), ROOT)
    return os.path.join(FIXTURES, rel)


def original(path):
    """The old original of PATH (a project of tests/ or examples/), else PATH."""
    f = _fixture(path)
    return f if os.path.isfile(f) else path


def copytree(src, dst, **kwargs):
    """shutil.copytree(src, dst), then every original kept for a project of src copied over its conversion."""
    out = shutil.copytree(src, dst, **kwargs)
    top = _fixture(src)
    if os.path.isdir(top):
        for dirpath, _, files in os.walk(top):
            for name in files:
                if name.endswith(".azr"):
                    rel = os.path.relpath(os.path.join(dirpath, name), top)
                    target = os.path.join(dst, rel)
                    os.makedirs(os.path.dirname(target), exist_ok=True)
                    shutil.copyfile(os.path.join(dirpath, name), target)
    return out
