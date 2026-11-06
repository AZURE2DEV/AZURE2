#!/usr/bin/env bash
#
# Shared per-point tables, pooled copies and sub-point grid re-placement (2026-11-05).
#
# A fit evaluates on max(4, threads) copies of the data (AZURECalc::InitializePools).  The copies share the
# energy-dependent tables of every point and sub-point (copy on write, EPoint::MutableEDep), and a re-placed
# sub-point grid keeps the tables of the nodes that did not move (EData::ReplaceSubPointGrids), so a copy costs
# a fraction of the data instead of a full copy (natural neon: under 1 GB when made, instead of ~7.7 GB).  Two races made such
# fits crash at four threads, in dev too: EData::Clone left the base points of segments with components mapping
# the ORIGINAL's component points, and the end of a CoulFuncCache::TransientScope on one thread cleared every
# thread's scratch memo (tests/reference/coulfunc_transient_scope_test.cpp checks that one on its own).
#
# Model: tests/natural_neon with every 16th data point (same segments, components and target effects); the fit
# is MIGRAD with the analytic gradient capped at 5 calls, which re-places sub-point grids at every evaluation.
#   (a) a calculation gives byte-identical output files at 1 and 4 threads;
#   (b) the capped fit completes at 1 and at 4 threads, re-placing grids, with the recorded chi-squared;
#   (c) its output files are byte-identical at 1 and 4 threads;
#   (d) the peak RSS of the 4-thread fit stays below FACTOR x that of the calculation (Linux only): measured
#       1795 MB against 755 MB (2.4x); with a full copy of the data per pooled copy (dev 2d97b1f0) 3513 MB (4.7x).
#
#   ./tests/shared_point_tables/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
SRC="$HERE/../natural_neon"
FACTOR="${FACTOR:-3.0}"  # measured 2.4 (dev 2d97b1f0, full copies: 4.7)
CHI2_FIT="$(cat "$HERE/expected/chi2_fit")"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/shared_point_tables.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# A copy of natural_neon with every 16th data point.
make_project() {  # DIR
  mkdir -p "$1/data" "$1/output" "$1/checks"
  cp "$SRC/natural_neon.azr" "$1/t.azr"
  for f in "$SRC"/data/*.dat; do
    awk '/^[[:space:]]*#/ || NF == 0 { print; next } { n++ } n % 16 == 1' "$f" > "$1/data/$(basename "$f")"
  done
}
# run DIR THREADS MENU [args]: runs AZURE2, prints its peak RSS in kB (0 where /proc is missing)
run() {
  local dir="$1" th="$2" menu="$3"; shift 3
  printf "$menu" | (cd "$dir" && OMP_NUM_THREADS="$th" AZURE_GRID_TIMING=1 exec ${TIMEOUT_CMD:+"$TIMEOUT_CMD" 1800} \
                      "$AZURE2_BIN" --no-gui --no-readline "$@" t.azr > run.log 2>&1) &
  local sub=$! pid="" hwm=0 h
  while kill -0 "$sub" 2> /dev/null; do
    # the AZURE2 process: the subshell itself, or the child of timeout
    [ -z "$pid" ] && { if [ -n "$TIMEOUT_CMD" ]; then pid="$(pgrep -P "$sub" 2> /dev/null | head -1)"; else pid="$sub"; fi; }
    if [ -n "$pid" ] && [ -r "/proc/$pid/status" ]; then
      h="$(awk '/VmHWM/ { print $2 }' "/proc/$pid/status" 2> /dev/null)"
      [ -n "$h" ] && [ "$h" -gt "$hwm" ] && hwm="$h"
    fi
    sleep 0.5
  done
  wait "$sub"; echo "$? $hwm"
}
chi2() { grep -h 'Total Chi-Squared' "$1/run.log" | tail -1 | awk '{ print $NF }'; }

for d in calc1 calc4 fit1 fit4; do make_project "$WORK/$d"; done
read -r rc1 hwm_calc <<< "$(run "$WORK/calc1" 1 '1\n\n\n7\n')"
read -r rc4 _ <<< "$(run "$WORK/calc4" 4 '1\n\n\n7\n')"
if [ "$rc1" -ne 0 ] || [ "$rc4" -ne 0 ]; then
  bad "(a) calculation failed (exit $rc1 / $rc4)"; tail -5 "$WORK/calc4/run.log"
else
  n=0; diff=0
  for f in "$WORK"/calc1/output/*; do
    n=$((n + 1)); cmp -s "$f" "$WORK/calc4/output/$(basename "$f")" || { diff=$((diff + 1)); echo "        differs: $(basename "$f")"; }
  done
  if [ "$n" -gt 0 ] && [ "$diff" -eq 0 ]; then ok "(a) calculation: $n output files identical at 1 and 4 threads"; else bad "(a) $diff of $n output files differ between 1 and 4 threads"; fi
fi

same_outputs() {  # DIR1 DIR2 -> "n d"
  local n=0 d=0 f
  for f in "$1"/output/*; do n=$((n + 1)); cmp -s "$f" "$2/output/$(basename "$f")" || { d=$((d + 1)); echo "        differs: $(basename "$f")" >&2; }; done
  echo "$n $d"
}
replaced() { grep -o '[0-9]* re-placed' "$1/run.log" | awk '{ s += $1 } END { print s + 0 }'; }
FIT='2\nn\n\n\n7\n'
read -r rc1 _ <<< "$(run "$WORK/fit1" 1 "$FIT" --use-gradient --max-calls=5)"
read -r rc4 hwm_fit <<< "$(run "$WORK/fit4" 4 "$FIT" --use-gradient --max-calls=5)"
for t in 1 4; do
  d="$WORK/fit$t"; rc="$( [ "$t" = 1 ] && echo "$rc1" || echo "$rc4")"; c="$(chi2 "$d")"; r="$(replaced "$d")"
  if [ "$rc" -eq 0 ] && [ "$r" -gt 0 ] &&
     awk -v a="$c" -v b="$CHI2_FIT" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(a != "" && d <= 1e-4 * b) }'; then
    ok "(b) capped fit, $t thread(s): chi2 $c (recorded $CHI2_FIT), $r grids re-placed"
  else
    bad "(b) capped fit, $t thread(s): exit $rc, chi2 '$c' (recorded $CHI2_FIT), $r grids re-placed"; tail -3 "$d/run.log"
  fi
done
read -r n d <<< "$(same_outputs "$WORK/fit1" "$WORK/fit4")"
if [ "$rc1" -eq 0 ] && [ "$rc4" -eq 0 ] && [ "$n" -gt 0 ] && [ "$d" -eq 0 ]; then ok "(c) capped fit: $n output files identical at 1 and 4 threads"; else bad "(c) capped fit: $d of $n output files differ between 1 and 4 threads"; fi

if [ "$hwm_calc" -gt 0 ] && [ "$hwm_fit" -gt 0 ]; then
  if awk -v f="$hwm_fit" -v c="$hwm_calc" -v k="$FACTOR" 'BEGIN { exit !(f <= k * c) }'; then
    ok "(d) peak RSS: fit at 4 threads $((hwm_fit / 1024)) MB <= $FACTOR x calculation $((hwm_calc / 1024)) MB"
  else
    bad "(d) peak RSS: fit at 4 threads $((hwm_fit / 1024)) MB > $FACTOR x calculation $((hwm_calc / 1024)) MB"
  fi
else
  echo "  skip  (d) peak RSS: no /proc/<pid>/status here"
fi

[ "$failures" -eq 0 ] && echo "shared_point_tables: all checks passed" || echo "shared_point_tables: $failures check(s) failed"
exit $((failures > 0))
