# shared_point_tables

Shared per-point tables, pooled copies and sub-point grid re-placement (2026-11-05).

A fit evaluates on `max(4, threads)` copies of the data (`AZURECalc::InitializePools`).  The copies share the
energy-dependent tables of every point and sub-point (copy on write), and a re-placed sub-point grid keeps the
tables of the nodes that did not move, so a copy costs a fraction of the data instead of a full copy
(`tests/natural_neon`: about 0.8 GB per copy instead of 7.7 GB; a fit at four threads 7.6 GB at the start and
9.2 GB after 50 minutes, instead of about 41 GB).  The tables are packed into flat arrays.  (Re-placed nodes were
shared between the copies through a cache with locks until 2026-11-12; on the shared mesh only the nodes near a moved
narrow level are re-placed, each copy computes its own.)  See the developer guide, *WorkingCopy, Data Points, and Target Effects*, "Memory: shared per-point
tables".

Two races made fits with grid re-placement crash at four threads (`malloc(): unaligned tcache chunk detected`,
segmentation faults), in dev 2d97b1f0 as well:

- `WorkingCopy::Clone` reset the mapped points only of the points the data iterator of that time visited, which skipped segments with
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
  the recorded chi-squared (`expected/chi2_fit`: 109.263831704 since the shared mesh, 2026-11-12; before it
  thm's 109.250542525, dev: 109.134 -- thm's adaptive sub-point grid, MERGE_NOTES.md section A, row 02ae97b);
- (c) its output files are byte-identical at 1 and 4 threads;
- (d) the peak RSS of the 4-thread fit is at most 3.0 times that of the calculation (Linux): 2.5 measured,
  4.7 with a full copy of the data per pooled copy.

About a minute at most.

2026-11-12, the shared mesh (`docs/upgrade/energy_grid.md` §2.2, §7.1): the folded segments of the cut-down project are
integrated on one mesh per compound nucleus.  At the parameters where the capped fit of the engine before it stopped
(its `param.sav`), the calculation gives 109.250994445 (engine before), 109.252912382 (shared mesh) and 109.25962319
(engine before, every grid refined eightfold): the mesh is the closer one, segment by segment as well.  The capped fit
then ends elsewhere after its five calls: 109.263831704 (shared mesh, 1 and 4 threads), 109.250542525 (before),
109.219622427 (refined eightfold); five calls are not a minimum, so the end point follows the path, not the rule.
