# tests/ — how AZURE2 is verified

Read the root `CLAUDE.md` first. This file covers which kind of test to write
and how to keep the references honest.

## Which kind of test

| Kind | Where | Use it for | Registered by |
|---|---|---|---|
| Physics regression | `tests/<name>/<name>.azr` + `data/` + `expected/chiSquared.out` | Pinning a whole evaluation's total χ², per-segment χ² and point counts | Automatic (`tests/run_tests.sh` finds every `.azr`) |
| Behavioural check | `tests/<name>/check.sh` (must be executable) | Sequences of runs, caches, CLI modes, error messages, fits that must recover known parameters, gradients vs finite differences | Automatic, *only if* `chmod +x` (`git update-index --chmod=+x`) |
| Reference | `tests/reference/<x>_test.cpp` | A quantity computed a second way: a closed form or an independent implementation. **The only kind that shows a number is right, not just unchanged.** | `tests/reference/CMakeLists.txt` |
| GUI | `tests/gui/<x>_test.cpp` | `.azr` round trips through the Qt models, dialogs, keyword whitelist | `tests/gui/CMakeLists.txt` (offscreen) |
| pyazr | `tests/pyazr/<x>_test.py` | Python reader/writer round trips; pyazr against the engine | `tests/pyazr/CMakeLists.txt` (`LABELS slow` if > 1 min) |
| Consistency | `tests/consistency/` | Repository drift: CLI flags vs help/docs, test registration, modes | `tests/consistency/CMakeLists.txt` |

A new feature usually needs a regression project or `check.sh`, a pyazr test
and a GUI round trip. A new formula needs a reference test.

## Writing one

- Follow the house pattern. A docstring or header comment says what is
  checked, against what, and which build first passed. A `check NAME OK
  DETAIL` helper prints `ok` or `FAIL` per item, and the test exits non-zero
  if anything failed. Copy from `tests/pyazr/norm_group_test.py`,
  `tests/shift_fn_groups/check.sh` or `tests/reference/*_test.cpp`.
- **Prove the test can fail.** Run it against the build before your fix, or
  break the code on purpose, and see it fail. Say in the commit message which
  you did.
- Use independent truth where possible: `truth.py` with `mpmath`, m-scheme
  sums, brute-force quadrature, published numbers with a citation. A value
  printed by the code under test is a regression reference, not a truth.
- `check.sh` works in `mktemp -d` with `trap 'rm -rf "$WORK"' EXIT`, sets
  `LC_ALL=C` (awk number parsing) and `OMP_NUM_THREADS=2`, and never writes
  into the source tree.
- Keep fixtures small: trim data files to the points the test needs. Keep
  every file under 1 MB.
- pyazr tests print `skip:` and exit 0 when `pyazr`/`_azure2` cannot be
  imported. A green ctest line can therefore be a skip. After touching pyazr,
  read the test output and confirm it says `ok`.
- A test directory gets a `README.md` when the model or data need explaining:
  where the data came from, and why this model.

## Re-baselining (`expected/`, tolerances, reference numbers)

Only for an intended, understood physics change, and never in the same commit
as the code change. The commit states, per project, the old and new totals,
the reason, and the independent evidence that the new number is right. A
`tests/<name>/tolerance` file overrides `TOL=1e-3` for one project and must
contain the reason as a `#` comment. If you are not sure the new number is
right, stop and ask.

## Running

```bash
./tests/run_tests.sh build/src/AZURE2              # all regression projects + check.sh
ctest --test-dir build -LE slow --output-on-failure
ctest --test-dir build -R pyazr_norm_group -V      # one test, full output
python3 tests/pyazr/norm_group_test.py             # pyazr tests also run directly
tests/norm_groups/check.sh build/src/AZURE2        # one check.sh directly
```

`run_tests.sh` leaves `output/`, `checks/` and `run.log` inside each project.
They are git-ignored build products: never commit them, and never copy them
into `expected/` without the re-baselining procedure above.
