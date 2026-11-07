# Golden outputs of T2

The outputs of today's engine on `tests/upgrade/T2`, kept so that each phase of the upgrade can show that it
reproduces them or explain where it does not (`docs/upgrade/ROADMAP.md`, the gate of every phase). `ENGINE` records
the engine commit, the build and the run settings.

| directory | run | compared at |
|---|---|---|
| `calc_start/` | menu 1 on `T2.azr` (the fit project at its start values) | 1 × 10⁻⁵ relative, every number |
| `calc_truth/` | menu 1 on `T2_true.azr` (the truth) | 1 × 10⁻⁵ relative |
| `fit/` | menu 2 with `--use-lm` on `T2.azr` (numerical MIGRAD to convergence, 1082 calls) | 1 × 10⁻³ relative |

The 10⁻⁵ is the spread between compilers stated in the repository's `CLAUDE.md`; the 10⁻³ of the fit allows the end
point of MIGRAD, set by its stopping rule (EDM < 2 × 10⁻⁴), to move by a small fraction of a parameter error when the
rounding changes. All runs use exact Coulomb functions (`--no-coulomb-table`). The calculations regenerated at 1 and
at 8 threads agree within 10⁻⁵; a second fit of the same project at 8 threads (the recovery check) agreed with
`fit/param.sav` to 1 × 10⁻⁸ relative.

```bash
AZURE2_BINARY=/path/to/AZURE2 python3 tests/upgrade/golden/T2/regen_compare.py            # compare (10 min at 8 threads)
AZURE2_BINARY=/path/to/AZURE2 python3 tests/upgrade/golden/T2/regen_compare.py --only calc_start --only calc_truth
AZURE2_BINARY=/path/to/AZURE2 python3 tests/upgrade/golden/T2/regen_compare.py --regenerate --commit <hash>
```

A difference outside the tolerance is a change of the engine's numbers; re-baseline only as `tests/CLAUDE.md`
("Re-baselining") says, in a commit of its own that states the old and new values and why they moved.
