"""rmfit bench: the benchmark suite of DESIGN.md Phase 3.  ``python3 -m rmfit.bench score <case.json> [run ...]`` scores
finished runs against the case's reference with the recovery yardstick (the evaluator 2026-11-09) and appends rows to
``rmfit/bench/results.tsv``; ``python3 -m rmfit.bench run <case.json> <run>`` writes (and with ``--submit`` queues) the SGE
job that builds the reference campaign, the perturbed start and runs ``rmfit drive`` on it.  Cases live in ``cases/``."""
