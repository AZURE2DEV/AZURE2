# shared_tables_race

Energy shifts on segments whose points are mapped onto points of another segment
(`EData::MapData`), evaluated with several threads.

Six segments at the same energies (the `shift_fn` test's linear-calibration data at 40 and
90 degrees, three copies each) form one shift group with a fitted shift function, so every
point of segments 2-6 is mapped onto a point of segment 1. On every fit step
`EData::FillEnergyShiftsFromParams` re-energises the points in parallel, and each mapped
point is re-energised by its own segment.

Before the fix, `EPoint::RecalcEDependentValues` on the host point also gave its new tables
to the points mapped onto it (`EPoint::CalcEDependentValues`, `mappedPoint->edep_ = edep_`).
That write came from the host's thread while the mapped point's own segment was re-energising
the point on another thread (`EPoint::MutableEDep`). Two threads assigned the same
`shared_ptr`, and the tables were freed under the other thread's copy-on-write: a
heap-use-after-free. Shared point tables (dev 3264896) turned an old overlapping write into a
crash. AddressSanitizer found it on 2026-11-06. It caused the macOS CI aborts in `shift_fn`
(the fits) and `pyazr_data_covariance`.

The test runs one fit with one thread, then `REPEATS` (default 6) `--use-gradient` Jacobian checks
(`AZURE_GRAD_CHECK`) and `REPEATS` fits with `OMP_NUM_THREADS` threads (default 8). It checks that:

1. every threaded gradient check finishes (exit 0, no sanitizer report or abort) and every shift
   column agrees with finite differences;
2. every threaded fit reaches the one-thread chi2 and shift parameters to 1e-6.

The gradient check is the path that re-energises the whole data set in parallel. A MIGRAD fit
evaluates its copies inside a parallel region, where the shift step runs serially.

**This test does not reliably fail on the old code.** The race needs two threads to collide in a
window of a few instructions. Under AddressSanitizer, `tests/shift_fn` at 8 threads hit it in 2 of
26 runs on 3264896 and 1 of 49 on 02b0001, and this test hit it in 0 of 40 runs on 6982601.
ThreadSanitizer (gcc 15, uninstrumented libgomp) catches the race without needing a collision:

- **6982601:** 888 reports, in a single run of this test, in which both conflicting accesses are in
  the energy-shift step and one of them is the host's `mappedPoint->edep_ = edep_`.
- **Fix:** 0 such reports.

The remaining reports are libgomp false positives. The record is
`R-matrix/AZURE2_code_audit/10-25-26/shared_tables_race/`. The test guards the threaded shift path
against crashes and wrong results.

    ./tests/shared_tables_race/check.sh path/to/AZURE2 [repeats]
