#!/usr/bin/env bash
#
# Physics regression tests.
#
# Runs every evaluation in this directory through AZURE2 and compares the
# resulting chi-squared against a recorded reference. This is what catches a
# change that still compiles and still runs but quietly moves the numbers --
# the smoke test only proves the binary starts.
#
#   ./tests/run_tests.sh [path/to/AZURE2]
#
# Adding a case: create tests/<name>/ containing <name>.azr, its data/, and
# expected/chiSquared.out from a run you trust. Nothing here needs editing --
# every directory holding an .azr is picked up automatically.
#
# TOL is the relative tolerance on chi-squared (default 1e-3). It has to absorb
# genuine cross-platform floating-point differences -- GCC, Clang and MinGW
# vectorise and contract differently, and libm results differ in the last bits --
# while still catching a real physics regression, which moves chi-squared by
# far more than a tenth of a percent.

set -uo pipefail

# Fixed locale: awk parses and prints numbers according to it, and a
# comma-decimal locale turns "24.335" into 24 on some awk implementations. The
# comparison below would then silently succeed against the wrong value.
export LC_ALL=C

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
AZURE2_BIN="${1:-}"
TOL="${TOL:-1e-3}"

# Several build directories may coexist (build/, build-py/, ...) and an old one
# testing green while the source has moved on is worse than no test at all, so
# take the newest binary rather than the first one found.
if [ -z "$AZURE2_BIN" ]; then
  newest=""
  for candidate in "$REPO_ROOT"/build*/src/AZURE2 "$REPO_ROOT"/build*/src/AZURE2.exe; do
    [ -x "$candidate" ] || continue
    if [ -z "$newest" ] || [ "$candidate" -nt "$newest" ]; then
      newest="$candidate"
    fi
  done
  AZURE2_BIN="$newest"
fi

if [ -z "$AZURE2_BIN" ] || [ ! -x "$AZURE2_BIN" ]; then
  echo "ERROR: AZURE2 binary not found. Pass it as the first argument." >&2
  exit 1
fi
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"

echo "AZURE2: $AZURE2_BIN"
echo "relative tolerance: $TOL"
echo

pass=0
fail=0

for project_dir in "$REPO_ROOT"/tests/*/; do
  azr="$(find "$project_dir" -maxdepth 1 -name '*.azr' | head -1)"
  [ -n "$azr" ] || continue

  name="$(basename "$project_dir")"
  expected="$project_dir/expected/chiSquared.out"

  if [ ! -f "$expected" ]; then
    echo "SKIP $name: no expected/chiSquared.out"
    continue
  fi

  echo "=== $name ==="

  # The .azr stores output/, checks/ and data paths relative to itself, so the
  # run has to happen from the project directory.
  (
    cd "$project_dir" || exit 1

    # Start clean: stale external-capture integral caches in output/ are reused
    # silently and would mask a real change.
    rm -rf output checks
    mkdir -p output checks

    # Menu: 1 = Calculate Segments From Data, then the external parameter file
    # and external capture amplitude file prompts (blank = build from the .azr),
    # then 7 = Exit. The prompt loops spin on EOF, so the output is capped as a
    # backstop -- SIGPIPE from head then stops the process.
    printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline "$(basename "$azr")" 2>&1 \
      | head -c 2000000 > run.log
  )

  actual="$project_dir/output/chiSquared.out"
  if [ ! -f "$actual" ]; then
    echo "  FAIL: no output/chiSquared.out produced"
    echo "  --- last 20 lines of run.log ---"
    tail -20 "$project_dir/run.log" 2>/dev/null | sed 's/^/  /'
    fail=$((fail + 1))
    continue
  fi

  # Compare the total, the per-segment chi-squared, and the per-segment point
  # counts. The point counts are integers and must match exactly: a change there
  # means data was dropped or misread, not a rounding difference.
  if awk -v tol="$TOL" '
    function relerr(a, b) { return (b == 0) ? ((a == 0) ? 0 : 1) : ((a - b) / b < 0 ? -(a - b) / b : (a - b) / b) }
    FNR == NR {
      if ($0 ~ /^Total-Chi-Squared:/) { split($0, t, " "); exp_total = t[2]; next }
      if ($0 ~ /^Segment#/) next
      split($0, f, ",")
      if (f[1] != "") { exp_chi[f[1]] = f[2]; exp_n[f[1]] = f[3] }
      next
    }
    {
      if ($0 ~ /^Total-Chi-Squared:/) { split($0, t, " "); act_total = t[2]; next }
      if ($0 ~ /^Segment#/) next
      split($0, f, ",")
      if (f[1] != "") { act_chi[f[1]] = f[2]; act_n[f[1]] = f[3] }
    }
    END {
      bad = 0
      for (s in exp_chi) {
        if (!(s in act_chi)) { printf "  segment %s missing from output\n", s; bad = 1; continue }
        if (act_n[s] != exp_n[s]) {
          printf "  segment %s: N %s, expected %s\n", s, act_n[s], exp_n[s]; bad = 1
        }
        e = relerr(act_chi[s], exp_chi[s])
        if (e > tol) {
          printf "  segment %s: chi2 %s, expected %s (rel err %.3g)\n", s, act_chi[s], exp_chi[s], e
          bad = 1
        }
      }
      for (s in act_chi) if (!(s in exp_chi)) { printf "  unexpected segment %s in output\n", s; bad = 1 }
      e = relerr(act_total, exp_total)
      printf "  total chi2 %s, expected %s (rel err %.3g)\n", act_total, exp_total, e
      if (e > tol) { print "  total chi-squared outside tolerance"; bad = 1 }
      exit bad
    }
  ' "$expected" "$actual"; then
    echo "  PASS"
    pass=$((pass + 1))
  else
    echo "  FAIL"
    fail=$((fail + 1))
  fi
  echo
done

echo "======================================"
echo "passed: $pass   failed: $fail"
[ "$fail" -eq 0 ] || exit 1
