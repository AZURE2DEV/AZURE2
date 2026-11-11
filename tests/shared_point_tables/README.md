# shared_point_tables

Shared per-point tables, pooled copies and sub-point grid re-placement (2026-11-05).

A fit evaluates on `max(4, threads)` copies of the data (`AZURECalc::InitializePools`).  The copies share the
energy-dependent tables of every point and sub-point (copy on write), and a re-placed sub-point grid keeps the
tables of the nodes that did not move, so a copy costs a fraction of the data instead of a full copy
(`tests/natural_neon`: about 0.8 GB per copy instead of 7.7 GB; a fit at four threads 7.6 GB at the start and
9.2 GB after 50 minutes, instead of about 41 GB).  Re-placed nodes are shared between the copies
(`NodeTableCache`), and the tables are packed into flat arrays.  See the developer guide, *EData, Data Points, and Target Effects*, "Memory: shared per-point
tables".

Two races made fits with grid re-placement crash at four threads (`malloc(): unaligned tcache chunk detected`,
segmentation faults), in dev 2d97b1f0 as well:

- `EData::Clone` reset the mapped points only of the points `EDataIterator` visits, which skips segments with
  components: the base points of natural-target and SUM segments kept mapping the ORIGINAL's component points,
  and every pooled copy filled them, from every thread at once.  `Clone` now checks that every mapped point
  belongs to the copy.
- `CoulFuncCache::TransientScope` (used by `CNuc::TransformOut`, which grid re-placement calls on every
  evaluation) counted its depth for all threads together, and the end of a scope on one thread cleared every
  thread's scratch memo while the others were using theirs.  The depth is now per thread;
  `tests/reference/coulfunc_transient_scope_test.cpp` checks it deterministically.

Model: `tests/natural_neon` with every 16th data point, built in a temporary directory by `check.sh`.

Checks (`check.sh path/to/AZURE2`):

- (a) a calculation gives byte-identical output files at 1 and 4 threads;
- (b) MIGRAD with the analytic gradient, capped at 5 calls, completes at 1 and 4 threads, re-placing grids, with
  the recorded chi-squared (`expected/chi2_fit`, the same as dev);
- (c) its output files are byte-identical at 1 and 4 threads;
- (d) the peak RSS of the 4-thread fit is at most 3.0 times that of the calculation (Linux): 2.5 measured,
  4.7 with a full copy of the data per pooled copy.

About a minute at most.
