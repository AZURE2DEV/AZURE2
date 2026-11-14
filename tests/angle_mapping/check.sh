#!/usr/bin/env bash
#
# Points at one energy across segments, angles and observables, and the energy-grouped analytic gradient, on a
# synthetic multi-angle project (angle_mapping.azr): two compound nuclei (p+12C,
# p+13C, the model of tests/natural_target), thirteen data segments on one energy
# grid, most of them one angle each -- cross sections at nine angles, one
# segment with two angles, two A_y segments, two natural-target ("sample 1")
# segments with free abundances under priors, free normalizations, a shift
# group with a free shift, a free sqrt(E) shift, a free shift function and a
# fixed shift of its own.
#
#   (a) the points of all segments at one compound energy share one level-matrix solve (the energy tasks,
#       docs/upgrade/energy_grid.md section 8), and that gives every unit, bit for bit, what its own solve gives
#       (AZURE_PLAN_CHECK=1: 0 units differ), with every output file of "calculate with data" byte-identical to the
#       run without the check;
#   (b) the sharing is real: fewer energy tasks than units (the one-angle segments, the A_y segments and the
#       natural-target segments sit at the energies of others).  (Until 2026-11-11 a point mapping made one point
#       calculate the others; it is deleted.)
#   (c) the analytic gradient (adjoint and residual-Jacobian paths) agrees with
#       central finite differences of chi^2 for every free parameter (levels,
#       norms, shifts, abundances) -- Brune (default), Lane-Thomas
#       (--no-brune) and Park (--use-park);
#   (d) the analytic gradient (adjoint and residual-Jacobian paths) is the same, bit for bit, at one thread and at
#       the thread count of (c): its sums are formed in an order fixed by the energy tasks, not by the threads
#       (src/GradEnergyTasks.cpp).
#
#   ./tests/angle_mapping/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# The time limit of tests/lib/guard.sh: timeout, gtimeout or a shell watchdog.
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/angle_mapping.XXXXXX")"
# Byte-identical files without cmp (not in every MSYS2 image); \r stripped so a
# CRLF-writing Windows build compares like the others.
same_file() { [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MODE "FLAGS" [VAR=VALUE ...] -> a fresh copy of the project in $WORK/NAME, run there
run() {
  local name="$1" mode="$2" flags="$3"
  shift 3
  local dir="$WORK/$name"
  rm -rf "$dir"
  mkdir -p "$dir/output" "$dir/checks"
  cp -r "$HERE/data" "$dir/"
  cp "$HERE/angle_mapping.azr" "$dir/t.azr"
  (cd "$dir" && for v in "$@"; do export "$v"; done &&
   printf "$mode\n\n\n7\n" | run_limited 900 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}

# ---- (a), (b) the shared solves against each unit's own -------------------------------------------
run plain 1 ""
run planchk 1 "" AZURE_PLAN_CHECK=1
n="$(ls "$WORK/plain/output" | wc -l | tr -d ' ')"
if [ "$n" -lt 3 ] || [ ! -s "$WORK/plain/output/chiSquared.out" ]; then
  bad "calculation wrote no output ($n files)"; tail -5 "$WORK/plain/run.log"
else
  diffs=""
  for f in "$WORK/plain/output"/*; do
    same_file "$f" "$WORK/planchk/output/$(basename "$f")" || diffs="$diffs $(basename "$f")"
  done
  if [ -z "$diffs" ]; then ok "(a) outputs ($n files) byte-identical with and without AZURE_PLAN_CHECK"
  else bad "(a) outputs differ under AZURE_PLAN_CHECK:$diffs"; fi
fi
# "AZURE_PLAN_CHECK: <units> units in <tasks> energy tasks (<shared> units at a shared energy) against their own solves:
# <differ> differ", once per pass
line="$(grep -m1 '^AZURE_PLAN_CHECK: .* units in .* energy tasks' "$WORK/planchk/run.log")"
read -r units tasks shared differ <<< "$(echo "$line" | awk '{ print $2, $5, substr($8, 2), $(NF - 1) }')"
if [ -n "$line" ] && [ "$differ" = 0 ] && ! grep -q '^AZURE_PLAN_CHECK: .* [1-9][0-9]* differ' "$WORK/planchk/run.log"; then
  ok "(a) every unit's shared solve = its own solve ($units units, 0 differ)"
else bad "(a) AZURE_PLAN_CHECK: '$line'"; fi
if [ -n "$line" ] && [ "$tasks" -lt "$units" ] && [ "$shared" -gt 0 ]; then
  ok "(b) $units units on $tasks energy tasks ($shared units at a shared energy)"
else bad "(b) no sharing: '$line'"; fi

# ---- (c), (d) gradients ------------------------------------------------------------------------
# rows NAME -> "idx analytic fd jacobian" of every checked parameter
rows() {
  awk '/^\[grad-check\] idx / { idx = $3; a = f = ""; j = "-"
         for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
         print idx, a, f, j }' "$WORK/$1/run.log" | sort -n
}
gradcheck() {  # NAME "FLAGS" [VAR=VALUE ...]
  local name="$1" flags="$2"
  shift 2
  run "$name" 2 "--use-gradient $flags" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
    AZURE_GRAD_CHECK_H=1e-8 "$@"
  rows "$name" > "$WORK/$name.rows"
}
for mode in brune:"" lanethomas:"--no-brune" park:"--use-park"; do
  tag="${mode%%:*}"; flags="${mode#*:}"
  gradcheck "g_$tag" "$flags"
  nrows="$(wc -l < "$WORK/g_$tag.rows" | tr -d ' ')"
  summary="$(grep "components with relative" "$WORK/g_$tag/run.log" | sed -n 1p | sed 's/.*analytic block //')"
  if [ "$nrows" -lt 20 ]; then
    bad "$tag: gradient check listed $nrows parameters"; grep -i "grad" "$WORK/g_$tag/run.log" | sed -n '1,5p'
    continue
  fi
  # finite differences with a 1e-8 relative step (the J = 3/2- level is narrow, 1e-6 is too
  # coarse for its energy); relative to the largest component, 1e-6 (Park: 1e-4, its level
  # matrix carries dS/dE differentiated numerically, see tests/pyazr/park_gradient_test.py)
  tol=1e-6; [ "$tag" = park ] && tol=1e-4
  worst="$(awk -v t="$tol" '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a; g[NR] = $2; f[NR] = $3; j[NR] = $4 }
             END { for (i = 1; i <= NR; i++) { d = g[i] - f[i]; if (d < 0) d = -d; e = j[i] - f[i]; if (e < 0) e = -e
                                               if (d > w) w = d; if (j[i] != "-" && e > w) w = e }
                   printf "%.2e", w / m }' "$WORK/g_$tag.rows")"
  if awk -v w="$worst" -v t="$tol" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 <= t + 0) }'; then
    ok "$tag: analytic and Jacobian-path gradients agree with finite differences for $nrows parameters (worst $worst of the largest)"
  else
    bad "$tag: gradient vs finite differences, worst $worst of the largest (> $tol); $summary"
  fi
  # the same gradient at one thread
  gradcheck "p_$tag" "$flags" OMP_NUM_THREADS=1
  same="$(join "$WORK/g_$tag.rows" "$WORK/p_$tag.rows" | awk '{ n++; if ($2 != $5 || $4 != $7) d++ }
              END { printf "%d %d", d + 0, n + 0 }')"
  if awk -v s="$same" 'BEGIN { split(s, x, " "); exit !(x[2] > 0 && x[1] == 0) }'; then
    ok "$tag: gradient at 1 and $OMP_NUM_THREADS threads identical (${same#* } parameters)"
  else
    bad "$tag: gradient at 1 and $OMP_NUM_THREADS threads: ${same% *} of ${same#* } parameters differ"
  fi
done

echo
if [ "$failures" -eq 0 ]; then echo "angle_mapping: all checks passed"; else echo "angle_mapping: $failures check(s) FAILED"; fi
[ "$failures" -eq 0 ]
