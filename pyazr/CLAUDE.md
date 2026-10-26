# pyazr/ — the Python API

Read the root `CLAUDE.md` first. pyazr is how most users script AZURE2, so its
public names are a contract.

## How it is built

- `azure2.py` drives the engine **in-process** through the pybind11 module
  `_azure2` (`api/src/pybind_api.cpp` → `api/src/AZUREAPI.cpp`). It runs no
  subprocess and opens no socket. Rebuilding with CMake (`-DUSE_API=ON` and
  pybind11 found) writes `pyazr/_azure2*.so` here. That file is git-ignored;
  never commit it.
- `azrfile.py` (`AzrModel`) and `datasets.py` (`SegmentSet`) read and write the
  `.azr` file in pure Python, independently of the engine and the GUI. They
  must round-trip every file the engine accepts **unchanged**, including
  inactive lines, unknown keyword tails and comments.
- `transform.py` is a Python port of the engine's `TransformOut`, and
  `widths.py` and `bands.py` re-derive engine quantities. Any such duplicate
  needs a test against the engine (see `tests/pyazr/`). Where possible, call
  through `_azure2` instead of adding another copy.
- `nds.py` (EXFOR/ENSDF client) shares its parsing rules with
  `gui/src/ExforData.cpp`. A rule learned by one belongs in both.

## Rules

- **Public API = `__all__` in `__init__.py`.** Adding a name is a minor version
  bump of `__version__` in the same file. Renaming or removing one needs a
  deprecation period first: keep the old name, have it emit
  `DeprecationWarning`, and note it in the docstring for at least one minor
  release. `tests/consistency` checks that every exported name exists.
- **Units and frames are in the name or the docstring** of every argument and
  return value, for example "lab energy, MeV" or "centre-of-mass, barns/sr".
  Return NumPy arrays, not lists, for numeric columns.
- Every public function has a docstring that says what it computes, its
  units, and what it raises. Match the style of the module you are editing;
  numpydoc sections (Parameters / Returns) are used in `azure2.py` and
  `bands.py`. A new public name also gets an entry in
  `docs/source/user_guide/pyazr.rst`. More than half of `__all__` is not
  documented there yet, so document what you touch.
- Errors: raise a specific exception with a message that names the file, line
  or parameter. Never return `None` or `nan` to signal failure, and never
  `print` a warning where `warnings.warn` belongs.
- Python ≥ 3.8 (`pyproject.toml`). Core dependencies are numpy, scipy and
  mpmath only. Anything else is an optional extra and is imported inside the
  function that needs it.
- Engine semantics (χ², penalties, normalizations, shifts, parameter order) are
  defined by the C++ code. pyazr must reproduce them exactly. Tests compare
  against the engine's `chiSquared.out` and parameter files, not against
  hand-computed expectations.
- `pyazr/examples/*.py` are user documentation. When an API changes, run the
  examples it affects.

## Tests

`tests/pyazr/<name>_test.py`, a plain script registered in
`tests/pyazr/CMakeLists.txt`, with `OMP_NUM_THREADS=2` and `LABELS slow` if it
takes more than a minute. Without the engine module, a test prints `skip:` and
exits 0. Read the output and confirm it says `ok`, not `skip`.
