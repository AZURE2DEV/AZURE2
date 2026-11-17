# rmfit on two backends: pyazr 2 and pyazr 3 in one package

Design note, 2026-11-17, for the owner's review before `feature/upgrade-5-dev-rmfit-fold` is merged into
`feature/upgrade-5-dev-rmfit`. The PI's decisions of 2026-11-16: "one rmfit, two backends"; uncertainty
quantification stays on pyazr 3 (`rmfit/sampling.py`); `mixture.py` is pyazr 2 only.

## What it is

`rmfit/` on this branch is the port of `feature/upgrade-5-dev-rmfit` (`543c7754`, HANDOFF entry "rmfit beside
pyazr") folded back into the one rmfit package that `dev` and the R-matrix archive also run. The same tree runs on
pyazr 2 (`dev`, the `.azr` text format, the archive's pinned builds) and on pyazr 3 (this branch, project format 3).
Its public copy is `github.com/rdeboer1/rmfit`; this branch carries that tree unchanged (`cc41f0e`, then
`1997494`: rmfit imports on Python < 3.12 again, and the dataset-off normgroup fix).

Why one package: the search, the ledger, the driver's state machine, the review criteria and the keys are the same
on both engines (the port changed them nowhere), and the archive's running campaigns stay on pyazr 2 until they end.
Two copies would diverge at the first fix; a fix to the shared code now reaches both.

## Where it plugs in

- `rmfit/backend.py` reads `pyazr.__version__` once, from the pyazr that `AZURE2_ROOT` puts first on the path, and
  sets `PYAZR3`. It hands out the backend's model editor, data-set table and nucleus class under one name each
  (`AzrModel`, `SegmentSet`, `Nucleus`: `pyazr.AzrModel` or `rmfit.azr3.Model`, and so on). A process runs one
  backend; nothing switches at run time.
- The port's modules are here unchanged: `azr3.py` (still rmfit's only contact with a format-3 project) and
  `sampling.py` (one comment anonymised by the public export).
- The engine and project layers are one module per backend behind the module the rest imports: `_worker.py`
  dispatches to `_worker_p2.py` (dev's worker, byte for byte) or `_worker_p3.py` (the port's); `model.py` keeps the
  shared shard planning and takes the file work from `_model_p2.py` or `_model_p3.py`.
- Elsewhere a contact point branches on `PYAZR3`; where a whole function is file work it has `_p2` / `_p3` bodies
  behind one dispatcher. Branch points per module: engine 12, audit 13, review 10, uq 8, windows 6, campaign 5,
  covariance 5, drive, local and structure 4 each.
- `mixture.py` refuses to import on pyazr 3, with the reason.

## What it replaces

The port's single-backend `rmfit/` (77 files modified, 57 added, none removed). On pyazr 3 the code paths are the
port's; the difference is the pyazr 2 branches beside them.

## How it is tested

- Tests: where the port left a dev script unchanged it is shared. Where it re-wrote one, `test_<name>.py` is the
  pyazr 2 version and `test_p3_<name>.py` the pyazr 3 one, and each exits 77 on the other backend.
  `rmfit/tests/CMakeLists.txt` registers both, and the fast/slow split matches `p3_` names to their base names.
  On this branch's build the `test_p3_*` scripts run and the pyazr 2 ones skip.
- pyazr 2 (shared dev build `b60aed09`): the 13N driver end to end gives the same verified candidates and check
  lines as the single-backend rmfit (`dev` `527f63de`). The fast set passes except `sanity_same_jpi` and
  `separate_close_levels`, which fail the same way before the fold (the port's HANDOFF entry gives the reason).
- pyazr 3 (build of `543c7754`): fast set 40 of 40; engine scripts 20 of 20; the driver 22 of 22 checks.
- On this branch head (`49adbcf9` + the fold, SGE 1766005): build, then every `rmfit_*` ctest: 80 passed
  (driver end to end 1936 s), 49 pyazr 2 versions skipped, none failed. Not run: the full gate, CI.

## For the owner

- A pyazr 2 branch inside `feature/upgrade` is code this branch's engine never runs. It is there because the same
  tree serves `dev` and the archive. If `rmfit/` is to follow this branch's rule that a phase deletes what it
  replaces, the pyazr 2 halves can be dropped when `dev` stops being used; `backend.py` makes them one `if` each.
- The port's open questions stand: which `rmfit-campaign` skill stays (`FOR_THE_OWNER.md` section 8), the size
  ratchet, no CI job for the rmfit tests, `W_PARAMETER_SETTING_UNKNOWN` at each shard start.
- A worker started with `AZURE2_ROOT` at an upgrade checkout imports that checkout's `rmfit/`, ahead of an rmfit
  elsewhere on the path. In-tree that is the rmfit under test; testing another rmfit against an upgrade build needs
  a view of the build without its `rmfit/`.
