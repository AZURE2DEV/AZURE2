# tests/ — how AZURE2 is verified

Read the root `CLAUDE.md` first. This file covers which kind of test to write
and how to keep the references honest.

## Which kind of test

| Kind | Where | Use it for | Registered by |
|---|---|---|---|
| Physics regression | `tests/<name>/<name>.azr` + `data/` + `expected/chiSquared.out` | Pinning a whole evaluation's total χ², per-segment χ² and point counts | Automatic (`tests/run_tests.sh` finds every `.azr`) |
| Behavioural check | `tests/<name>/check.sh` (must be executable) | Sequences of runs, caches, CLI modes, error messages, fits that must recover known parameters, gradients vs finite differences | Automatic, *only if* `chmod +x` (`git update-index --chmod=+x`) |
| Reference | `tests/reference/<x>_test.cpp` | A quantity computed a second way: a closed form or an independent implementation. **The only kind that shows a number is right, not just unchanged.** | `tests/reference/CMakeLists.txt` |
| pyazr | `tests/pyazr/<x>_test.py` | Python reader/writer round trips; pyazr against the engine | `tests/pyazr/CMakeLists.txt` (`LABELS slow` if > 1 min) |
| Consistency | `tests/consistency/` | Repository drift: CLI flags vs help/docs, test registration, modes | `tests/consistency/CMakeLists.txt` |

A new feature usually needs a regression project or `check.sh` and a pyazr
test (with a round trip). A new formula needs a reference test.

A regression reference is not a correctness check: `tests/run_tests.sh` pins
each project's χ² against a number this code produced, so it catches a change
and nothing else. `tests/reference/beam_profile_reference_test.cpp` is the
pattern for an integrator: it drives `FoldPoint` (`src/PointFold.cpp`) with
hand-chosen sub-point cross sections, needs no model, runs in 0.1 s, and rests
on an identity (a constant cross section folds back to itself exactly, whatever
the profile, window and detailed-balance weight, because numerator and
denominator share the kernel). Mutation-test such a check: forcing the window
shift to zero, or dropping the detailed-balance factor, each fails exactly one
case.

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
- An m-scheme truth for capture interference (as `tests/ec_e1m1_phase`) has
  three traps: with a spin-0 target and unpolarized beam, E1 and a pure
  spin-flip M1 (s-wave final state) do not interfere, so such a case cannot
  test the phase; the final state is the |(l_f s)J_f⟩ level, not every product
  state of the orbital; and the ANC must be physical for the channel radius,
  or the cross section stops scaling as ANC².
- `check.sh` works in `mktemp -d` with `trap 'rm -rf "$WORK"' EXIT`, sets
  `LC_ALL=C` (awk number parsing) and `OMP_NUM_THREADS=2`, and never writes
  into the source tree.
- Under `set -o pipefail` never pipe into a reader that leaves early
  (`grep -q`/`-m`, `head -n`, `sed …q`, `awk` with `exit` outside `END`):
  the writer's SIGPIPE makes the pipeline fail at random on a loaded machine
  (a match read as "not found"); use `grep -q PATTERN <<< "$text"` or capture
  the command's output first [`repo_consistency_test.py` check 50].
- Keep fixtures small: trim data files to the points the test needs. Keep
  every file under 1 MB.
- pyazr tests print `skip:` and exit 77 when `pyazr`/`_azure2` cannot be
  imported; ctest shows them as skipped (`SKIP_RETURN_CODE 77`), never as
  passed. A run with skipped tests has not tested pyazr: give the build its
  engine module (`PYAZR_ENGINE_DIR`).
- A test directory gets a `README.md` when the model or data need explaining:
  where the data came from, and why this model.

## Portability (CI runs Linux, macOS and Windows/MinGW)

CI builds and tests with numpy and pybind11 on all three, so the pyazr tests
run there too. Before pushing:

1. No `setenv`/`unsetenv` in C++ (MinGW lacks them); use `setEnvVar` /
   `unsetEnvVar` at the top of `src/AZURECalc.cpp` (`_putenv` on Windows).
2. No GNU-only shell in tests (`sed -i` without a suffix, `sed '0,/re/'`,
   `grep -P`); macOS runs BSD tools. awk and perl are portable.
3. A pyazr test imports numpy and pyazr inside its skip guard
   (`try: ... except: print("skip: ..."); sys.exit(77)`).
4. CI runs `ctest -LE slow`. A test that takes minutes gets `LABELS slow` in
   `tests/pyazr/CMakeLists.txt`; run `ctest -L slow` locally (the longest,
   `park_gradient`, takes about 40 min) when the pyazr or engine paths those
   tests cover change.

Job logs need a GitHub token; annotations are public but capped at 10 per step.
Pushing a change to `.github/workflows/` needs the token's `workflow` scope.

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
