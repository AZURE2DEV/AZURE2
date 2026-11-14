#!/usr/bin/env bash
#
# A fit does not depend on the number of threads (docs/upgrade/fit_threads.md).  Model: tests/upgrade/T2, MIGRAD with
# Minuit's numerical (OpenMP) gradient, as scripts/bench_upgrade.sh runs it, bounded to 30 calls.  AZURE_TRACE_FCN
# prints every objective value with its parameter vector in hexadecimal floating point.
#
#   1. a pooled working copy gives, at every parameter vector of the fit, the bits a fresh copy gives
#      (AZURE_TRACE_FCN=fresh, 2 threads).  Before the fix the re-placed sub-point grids kept nodes within 1e-12 of
#      an earlier grid, so a copy's value depended on the vectors it had evaluated before (26 of 79 fit evaluations
#      differed in their last bits, the first at evaluation 15);
#   2. fits at 1 and 4 threads evaluate the same parameter vectors and give, at each, the same objective bit for bit.
#      The number of evaluations is not compared: Minuit2's call counter is not atomic (MnUserFcn.cxx), so at 4
#      threads it loses counts and a --max-calls fit may run one more iteration (fit_threads.md section 3).
#
# Shown to fail with the tolerance restored in EData::ReplaceSubPointGrids (2026-11-09): item 1 as above, item 2 with
# 34 of 59 common vectors at different values.  Takes about a minute (three 30-call fits of T2).
#
#   ./tests/fit_threads/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL="$HERE/../upgrade/T2"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/fit_threads.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

fit() {  # NAME THREADS TRACE
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$MODEL/T2.azr" "$MODEL/cov_diff.txt" "$MODEL/prior_cov.f3.txt" "$MODEL/data" "$dir/"
  (cd "$dir" && AZURE_TRACE_FCN="$3" OMP_NUM_THREADS="$2" "$AZURE2_BIN" fit --no-gui \
     --no-coulomb-table --use-lm --max-calls=30 T2.azr > run.log 2> trace.log)
  if [ ! -s "$dir/output/chiSquared.out" ] || ! grep -q '^\[fcn\]' "$dir/trace.log"; then
    bad "$1: no fit or no trace"; tail -5 "$dir/run.log"; return 1
  fi
}

# ---- 1. pooled copy == fresh copy ------------------------------------------------------------------------------
if fit fresh 2 fresh; then
  read -r n nd first <<< "$(awk '/^\[fcn\]/ && /\|\| fresh/ {
      split($0, a, "|"); split(a[1], h, " "); split(a[4], f, " "); n++
      if (h[3] != f[2]) { nd++; if (!first) first = h[2] } }
      END { print n + 0, nd + 0, first + 0 }' "$WORK/fresh/trace.log")"
  if [ "$n" -lt 20 ]; then bad "pooled vs fresh: only $n fit evaluations traced"
  elif [ "$nd" -eq 0 ]; then ok "pooled copy = fresh copy, bit for bit, at all $n fit evaluations"
  else bad "pooled copy differs from a fresh copy at $nd of $n fit evaluations (first: evaluation $first)"; fi
fi

# ---- 2. 1 and 4 threads: the same vectors, the same values ------------------------------------------------------
if fit t1 1 1 && fit t4 4 1; then
  for t in t1 t4; do awk '/^\[fcn\]/ { split($0, a, "|"); split(a[1], h, " "); print a[2] "\t" h[3] }' \
    "$WORK/$t/trace.log" | sort -u > "$WORK/$t.values"; done
  read -r n1 common nd <<< "$(awk -F'\t' 'NR == FNR { v[$1] = $2; n1++; next }
      ($1 in v) { c++; if (v[$1] != $2) nd++ } END { print n1 + 0, c + 0, nd + 0 }' "$WORK/t1.values" "$WORK/t4.values")"
  if [ "$nd" -eq 0 ] && [ "$common" -ge $((n1 - 1)) ] && [ "$n1" -ge 20 ]; then
    ok "1 and 4 threads: $common of the $n1 vectors of the 1-thread fit evaluated at 4 threads too, same bits at each"
  else bad "1 and 4 threads: $common of $n1 vectors in common, $nd with different values"; fi
fi

[ "$failures" -eq 0 ] && echo "fit_threads: all checks passed" || { echo "fit_threads: $failures check(s) FAILED"; exit 1; }
