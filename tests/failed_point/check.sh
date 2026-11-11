#!/usr/bin/env bash
#
# A point whose calculation fails is NaN in the CLI too, and the run says so.
#
# The CLI's objective takes the values of the parallel pre-pass (EData::CalculatePointsInParallel), which before
# 2026-11-01 set a failed point to 0 without a word (and ESegment::CalculateTheoreticalCrossSection gave 0 or left
# a SUM term out).  The failure is injected with the test hook AZURE_TEST_FAIL_E=<lab energy> (here the base term
# of tests/sum_lab_frame's SUM segment at one energy; the segment carries a 0.05 MeV energy shift, so 6.15 -> 6.2).
# Checked: mode 1 prints the warning and chiSquared.out's total is nan, the output file has nan at that point only;
# without the hook the same run is finite.  pyazr: tests/pyazr/failed_point_test.py.
#
#   ./tests/failed_point/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/failed_point.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME [VAR=VALUE...]: mode 1 of the SUM model under $WORK/NAME with the given environment.
run() {
  local dir="$WORK/$1"; shift
  mkdir -p "$dir/output" "$dir/checks" "$dir/data"
  cp "$HERE/../sum_lab_frame/models/sum_data.azr" "$dir/"
  cp "$HERE/../sum_lab_frame/models/data/"* "$dir/data/"
  for e in 5.95 6.05 6.15 6.25 6.35; do echo "$e 90 1.0 0.01"; done > "$dir/data/sum.dat"
  (cd "$dir" && printf '1\n\n\n7\n' | env "$@" "$AZURE2_BIN" --no-gui --no-readline sum_data.azr > run.log 2>&1)
}
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2>/dev/null; }

run ref AZURE_UNUSED=1
t="$(total ref)"
case "$t" in ''|nan|-nan|inf) bad "reference run: total chi2 '$t'";; *) ok "reference run: total chi2 $t";; esac

run fail AZURE_TEST_FAIL_E=6.2 AZURE_TEST_FAIL_EXIT=3
grep -q "WARNING: the calculation failed at segment 1 .*injected failure" "$WORK/fail/run.log" \
  && ok "failed point: the run prints the warning" || bad "failed point: no warning in the log"
t="$(total fail)"
case "$t" in nan|-nan) ok "failed point: total chi2 is nan";; *) bad "failed point: total chi2 '$t' (expected nan)";; esac
out="$(ls "$WORK/fail/output/"AZUREOut_*.out 2>/dev/null | head -1)"
n="$(awk '$4 ~ /nan/ { n++ } END { print n + 0 }' "$out" 2>/dev/null)"
[ "$n" = 1 ] && ok "failed point: one nan point in the output file" || bad "failed point: $n nan points in '$out'"

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
