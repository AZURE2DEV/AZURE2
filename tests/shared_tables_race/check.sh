#!/usr/bin/env bash
#
# Energy shifts on segments whose points are mapped across segments (EData::MapData), threaded.
#
# Six segments at the same energies (the shift_fn test's linear-calibration data at 40 and 90 deg, three
# copies each) share one shift group with a shift function, so every point of segments 2-6 is mapped onto a
# point of segment 1.  On every fit step EData::FillEnergyShiftsFromParams re-energises the points in
# parallel, each mapped point by its own segment.  Before the fix, the host's RecalcEDependentValues also
# assigned its tables to those mapped points (EPoint::CalcEDependentValues), from another thread, while their
# own segment was re-energising them: two threads assigning one shared_ptr, a heap-use-after-free
# (AddressSanitizer; the macOS CI aborts of tests/shift_fn and pyazr_data_covariance, 2026-11-06).
#
# The parallel re-energising runs where shifts are applied to the whole data set outside a parallel region:
# the --use-gradient Jacobian check (AZURE_GRAD_CHECK) does that for every finite-difference step (a MIGRAD
# fit evaluates its copies inside a parallel region, where the shift step runs serially).  The race is
# intermittent: a plain build mostly survives it on Linux.  The checks are
#   1. R gradient checks at OMP_NUM_THREADS (default 8) all finish (no crash, no abort, no sanitizer report)
#      and every shift column agrees with finite differences;
#   2. R threaded fits reach the same chi2 and parameters as a one-thread fit (1e-6 relative).
# Under an AddressSanitizer build (as the 2026-11-06 audit did) the old code fails check 1.
#
#   ./tests/shared_tables_race/check.sh path/to/AZURE2 [repeats]
set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../shift_fn"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2 [repeats]}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
REPEATS="${2:-6}"
THREADS="${OMP_NUM_THREADS:-8}"
[ "$THREADS" -lt 2 ] && THREADS=8
WORK="$(mktemp -d "${TMPDIR:-/tmp}/shared_tables_race.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] && echo "WORK=$WORK" || trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

mkdir -p "$WORK/model"
cp -R "$SRC/data" "$WORK/model/data"
{
  echo '1 1 1 0 100 40 40 4 1 0 0 0 1 1 data/xs_40_lin_mislabelled.dat 0 0 shiftgroup G shiftfn "a0*(x-1)" 1 0 1 1'
  echo '1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90_lin_mislabelled.dat 0 0 shiftgroup G'
  echo '1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40_lin_mislabelled.dat 0 0 shiftgroup G'
  echo '1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90_lin_mislabelled.dat 0 0 shiftgroup G'
  echo '1 1 1 0 100 40 40 4 1 0 0 0 0 0 data/xs_40_lin_mislabelled.dat 0 0 shiftgroup G'
  echo '1 1 1 0 100 90 90 4 1 0 0 0 0 0 data/xs_90_lin_mislabelled.dat 0 0 shiftgroup G'
} > "$WORK/segs.txt"
awk -v segs="$WORK/segs.txt" '
  /<levels>/ { inlev = 1; print; next }
  /<\/levels>/ { inlev = 0 }
  inlev && NF > 12 { $4 = 1; $11 = 1; print; next }   # levels fixed: only the shifts move
  /<segmentsData>/ { print; while ((getline l < segs) > 0) print l; skip = 1; next }
  /<\/segmentsData>/ { skip = 0 }
  skip { next }
  { print }' "$SRC/shift_fn.azr" > "$WORK/model/t.azr"

fit() {  # NAME THREADS
  local dir="$WORK/$1"
  rm -rf "$dir"; cp -R "$WORK/model" "$dir"; mkdir -p "$dir/output" "$dir/checks"
  (cd "$dir" && printf '2\n\n\n\n7\n' | OMP_NUM_THREADS="$2" "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  echo $?
}
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
shifts() { awk '$1 ~ /energy_shift/ { printf "%s ", $2 }' "$1/output/param.sav" 2>/dev/null; }
same() {  # two blank-separated number lists agree to 1e-6 relative
  awk -v a="$1" -v b="$2" 'BEGIN { n = split(a, x, " "); m = split(b, y, " "); if (n != m || n == 0) exit 1
    for (i = 1; i <= n; i++) { d = x[i] - y[i]; if (d < 0) d = -d; s = y[i] < 0 ? -y[i] : y[i]; if (s < 1e-12) s = 1e-12; if (d / s > 1e-6) exit 1 } }'
}

st="$(fit ref 1)"
ref_chi="$(total "$WORK/ref")"; ref_sh="$(shifts "$WORK/ref")"
if [ "$st" = 0 ] && [ -n "$ref_chi" ] && [ -n "$ref_sh" ]; then
  ok "one-thread fit: chi2 $ref_chi, shift parameters $ref_sh"
else
  bad "one-thread fit failed (exit $st)"; tail -15 "$WORK/ref/run.log"
fi
mapped="$(grep -c . "$WORK/model/data/xs_40_lin_mislabelled.dat")"
for r in $(seq 1 "$REPEATS"); do
  dir="$WORK/g$r"; rm -rf "$dir"; cp -R "$WORK/model" "$dir"; mkdir -p "$dir/output" "$dir/checks"
  (cd "$dir" && printf '2\n\n\n\n7\n' | OMP_NUM_THREADS="$THREADS" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=50 \
     "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
  st=$?
  if [ "$st" != 0 ] || grep -qE "AddressSanitizer|ThreadSanitizer|terminate called|Abort" "$dir/run.log"; then
    bad "gradient check $r ($THREADS threads): exit $st"; grep -m1 -A12 -E "ERROR: AddressSanitizer|WARNING: ThreadSanitizer|terminate called" "$dir/run.log"
  elif grep -q "Jacobian path: shift 0/1 shift_fn 0/1" "$dir/run.log"; then
    ok "gradient check $r ($THREADS threads): shift columns agree with finite differences"
  else
    bad "gradient check $r ($THREADS threads): $(grep -m1 'components with relative' "$dir/run.log")"
  fi
done
for r in $(seq 1 "$REPEATS"); do
  st="$(fit "r$r" "$THREADS")"
  log="$WORK/r$r/run.log"
  if [ "$st" != 0 ] || grep -qE "AddressSanitizer|ThreadSanitizer|terminate called|Abort" "$log"; then
    bad "run $r ($THREADS threads): exit $st"; grep -m1 -A12 -E "ERROR: AddressSanitizer|WARNING: ThreadSanitizer|terminate called" "$log"
    continue
  fi
  c="$(total "$WORK/r$r")"; s="$(shifts "$WORK/r$r")"
  if same "$c $s" "$ref_chi $ref_sh"; then ok "run $r ($THREADS threads): same fit (chi2 $c)"
  else bad "run $r ($THREADS threads): chi2 $c, shifts $s; one thread: $ref_chi, $ref_sh"; fi
done
if [ "$failures" -eq 0 ]; then echo "shared_tables_race: all checks passed"; else echo "shared_tables_race: $failures check(s) FAILED"; exit 1; fi
