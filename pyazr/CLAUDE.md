# pyazr/ — the Python API

Read the root `CLAUDE.md` first. pyazr is how most users script AZURE2, so its
public names are a contract.

## Using pyazr as an agent: ask the package first

- `help(pyazr)`: the objects of today's API and a ten-line session (tested by `tests/pyazr/describe_test.py`).
- `pyazr.describe()`: version, engine build, Coulomb table in use (path, content identity or "exact"), stopping and
  mass data, threads.  `pyazr.describe(path)` is still the engine's description of a project.
- `pyazr.kinds()`: record kinds and fields, data kinds (units, frames), parameter kinds (unit, owner fields), session
  options (type, default, CLI flag, doc), bound physics functions, error classes.  Read from the engine's tables
  (`described.py`); never add a hand-typed list of kinds in Python.
- `azure2(...).describe()`: pairs, levels, data sets, parameters (name, kind, unit, value, free), options in effect.
- Errors (`errors.py`): every deliberate error is a `pyazr.Error`; `FormatError`, `ParameterError` (with
  `ParameterKeyError`, `ParameterIndexError`, `ParameterVectorError`), `EngineError` (= `AZURE2Error`, one class,
  which keeps that name), `DomainError` (`DomainZeroDivision`), `SessionClosed`, `MissingFile`; engine notices are
  `EngineWarning`.  A new error derives from the class above and from the builtin the call raised before; its
  message says what was wrong, what is allowed and which call shows or fixes it (`tests/pyazr/errors_test.py`
  holds the table of wrong calls).
  The argument checks of the session's calls live in `session_checks.py`, not in `azure2.py`.
- An engine notice becomes a warning when its line is a row of `api/src/pybind_notices.cpp` (the session's
  `Config::outStream` passes through that buffer); every other line still goes to standard output.
- `pyazr.physics.<name>`: the engine's physics functions on arrays (`help()` gives arguments, units, domains);
  `pyazr.kinds()["physics"]` the same rows as data; `errors="mask"` gives NaN and a status per element instead of
  `DomainError`.  The per-session ones (the session's Coulomb table or potential) are `azure2(...).channel_functions`.

## Recipe: add a physics function

1. The library computes it (`physics/`, `experiment/`, `stopping/`); if it does not, a thin free function there with
   its unit test (as `physics/include/azure2/physics/Widths.h`).  New physics is a design note first, not a row.
2. One row in `BuildTable()` of `api/src/pybind_physics.cpp`: name, summary, doc, library function, axis, arguments
   (`Real`/`Integer`/`Flag`/`Choice` with unit and domain: `.gt`, `.ge`, `.lt`, `.le`, `.otherwise` for a default),
   returns (name, unit, doc), and an adapter that converts the arguments and calls the library (no formula there);
   the adapter returns 0, `kJointly` (-2) for arguments outside the domain together, `kFailed` (-1), or the 1-based
   index of an argument.
3. A `row_test("<name>", ...)` in `tests/pyazr/physics_table_test.py` against mpmath, a published table or the
   engine's own use on a project (consistency check 39 fails without one).
4. `python3 tools/gen_physics_docs.py` (the user-guide page; the test fails until it is run).
Nothing in Python: `pyazr/physics.py` builds the function, its signature and docstring from the row.

## How it is built

- `azure2.py` drives the engine **in-process** through the pybind11 module
  `_azure2` (`api/src/pybind_api.cpp` → `api/src/AZUREAPI.cpp`). It runs no
  subprocess and opens no socket. Rebuilding with CMake (`-DUSE_API=ON` and
  pybind11 found) writes `pyazr/_azure2*.so` here. That file is git-ignored;
  never commit it.
- pyazr reads project format 3 only.  `project.py` (`Project`, records) and
  `nuclear_part.py` (`<levels>`, `<config>`, `<thm>`) edit the file on the
  engine's reading of it: the lossless document (`_azure2.Document`), the
  record grammar and validator, the level edits (`experiment/LevelEdits.h`),
  the `<config>` and `<thm>` edits. They need the compiled module. Never add a
  Python parser or a column position for a section: bind the engine's. A
  project must round-trip **byte for byte** (any layout, any line ending), and
  an edit changes only the lines it edits (`tests/pyazr/project_test.py`,
  `azr_layout_test.py`). A format-1/2 file is refused in one place,
  `project.require_format3` (the engine's message, io::NotFormat3Message, names
  `AZURE2 --convert OLD NEW`); `convert.py` is the converter's binding
  (`convert_azr`, `convert_params`).
- `transform.py` is a Python port of the engine's `TransformOut`, and
  `widths.py` and `bands.py` re-derive engine quantities. Their formulas are
  compared with the rows of `pyazr.physics` that replace them
  (`tests/pyazr/physics_table_test.py`) and go in the rename step. Never add
  another copy: a missing function is a row of the physics table.
- `nds.py` (EXFOR/ENSDF client) is the only EXFOR reader on this branch (the
  Qt GUI's `gui/src/ExforData.cpp` exists on `dev` only).

## Rules

- **Public API = `__all__` in `__init__.py`.** Adding a name is a minor version
  bump of `__version__` in the same file. Renaming or removing one within a
  major version needs a deprecation period first (3.0 removed the old-format
  API at once: the converter handles old files): keep the old name, have it emit
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
