#!/usr/bin/env bash
#
# Point mapping across segments and the energy-grouped analytic gradient, on a
# synthetic multi-angle project (angle_mapping.azr): two compound nuclei (p+12C,
# p+13C, the model of tests/natural_target), thirteen data segments on one energy
# grid, most of them one angle each -- cross sections at nine angles, one
# segment with two angles, two A_y segments, two natural-target ("sample 1")
# segments with free abundances under priors, free normalizations, a shift
# group with a free shift, a free sqrt(E) shift, a free shift function and a
# fixed shift of its own.
#
#   (a) every output file of "calculate with data" is byte-identical with
#       mapping across segments (the default), with AZURE_NO_CROSS_SEGMENT_MAPPING=1
#       (within each segment only) and with AZURE_NO_POINT_MAPPING=1;
#   (b) the data check shows the mapping: points of one-angle segments are
#       mapped onto other segments; the two shift-group members only onto each
#       other; the segments with their own free sqrt shift, free shift function
#       or different fixed shift onto none;
#   (c) the analytic gradient (adjoint and residual-Jacobian paths) agrees with
#       central finite differences of chi^2 for every free parameter (levels,
#       norms, shifts, abundances) -- Brune (default), Lane-Thomas
#       (--no-brune) and Park (--use-park);
#   (d) the energy-grouped gradient equals the point-by-point one
#       (AZURE_NO_GROUPED_GRADIENT=1) to 1e-10 of its largest component.
#
#   ./tests/angle_mapping/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/angle_mapping.XXXXXX")"
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

# ---- (a) model outputs with and without mapping --------------------------------------------------
run map 1 ""
run within 1 "" AZURE_NO_CROSS_SEGMENT_MAPPING=1
run nomap 1 "" AZURE_NO_POINT_MAPPING=1
n="$(ls "$WORK/map/output" | wc -l | tr -d ' ')"
if [ "$n" -lt 3 ] || [ ! -s "$WORK/map/output/chiSquared.out" ]; then
  bad "calculation wrote no output ($n files)"; tail -5 "$WORK/map/run.log"
else
  for other in within nomap; do
    diffs=""
    for f in "$WORK/map/output"/*; do
      cmp -s "$f" "$WORK/$other/output/$(basename "$f")" || diffs="$diffs $(basename "$f")"
    done
    if [ -z "$diffs" ]; then ok "outputs ($n files) byte-identical with mapping across segments and '$other'"
    else bad "outputs differ between mapping across segments and '$other':$diffs"; fi
  done
fi

# ---- (b) what is mapped onto what ---------------------------------------------------------------
dir="$WORK/datachk"
mkdir -p "$dir/output" "$dir/checks"
cp -r "$HERE/data" "$dir/"
sed 's/^none #Data Check/file #Data Check/' "$HERE/angle_mapping.azr" > "$dir/t.azr"
(cd "$dir" && printf '1\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
chk="$dir/checks/data.chk"
if [ ! -s "$chk" ]; then
  bad "no data check written"
else
  # per segment: "host segments" of its mapped points, as a sorted list
  hosts() { awk -v s="$1" '$1 == s && $8 ~ /^\(/ { split(substr($8, 2), h, ","); print h[1] }' "$chk" | sort -un | tr '\n' ' ' | sed 's/ $//'; }
  [ "$(hosts 2)" = "1" ] && ok "segment 2 (90 deg) mapped onto segment 1 (60 deg)" || bad "segment 2 mapped onto segments '$(hosts 2)', want '1'"
  [ "$(hosts 3)" = "1" ] && ok "segment 3 (120 and 150 deg) mapped onto segment 1" || bad "segment 3 mapped onto '$(hosts 3)', want '1'"
  [ "$(hosts 5)" = "4" ] && ok "A_y segment 5 mapped onto A_y segment 4, not onto cross sections" || bad "segment 5 mapped onto '$(hosts 5)', want '4'"
  [ "$(hosts 7)" = "6" ] && [ -z "$(hosts 6)" ] && ok "shift-group members (free shared shift) mapped onto each other only" \
    || bad "shift group: segment 7 onto '$(hosts 7)', segment 6 onto '$(hosts 6)'"
  [ -z "$(hosts 8)$(hosts 9)$(hosts 13)" ] && ok "free sqrt shift, free shift function, other fixed shift: not mapped" \
    || bad "segments 8, 9, 13 mapped onto '$(hosts 8)' '$(hosts 9)' '$(hosts 13)'"
  [ "$(hosts 11)" = "10" ] && [ "$(hosts 12)" = "10" ] && ok "p+13C segments (natural and not) mapped onto segment 10" \
    || bad "p+13C: segment 11 onto '$(hosts 11)', 12 onto '$(hosts 12)'"
fi

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
  summary="$(grep "components with relative" "$WORK/g_$tag/run.log" | head -1 | sed 's/.*analytic block //')"
  if [ "$nrows" -lt 20 ]; then
    bad "$tag: gradient check listed $nrows parameters"; grep -i "grad" "$WORK/g_$tag/run.log" | head -5
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
  if awk -v w="$worst" -v t="$tol" 'BEGIN { exit !(w + 0 <= t + 0) }'; then
    ok "$tag: analytic and Jacobian-path gradients agree with finite differences for $nrows parameters (worst $worst of the largest)"
  else
    bad "$tag: gradient vs finite differences, worst $worst of the largest (> $tol); $summary"
  fi
  # the point-by-point gradient
  gradcheck "p_$tag" "$flags" AZURE_NO_GROUPED_GRADIENT=1
  same="$(join "$WORK/g_$tag.rows" "$WORK/p_$tag.rows" | awk '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a
              d = $2 - $5; if (d < 0) d = -d; if (d > w) w = d
              if ($4 != "-" && $7 != "-") { e = $4 - $7; if (e < 0) e = -e; if (e > w) w = e }; n++ }
              END { if (n == 0 || m == 0) print "nan"; else printf "%.2e %d", w / m, n }')"
  if awk -v s="$same" 'BEGIN { split(s, x, " "); exit !(x[2] > 0 && x[1] + 0 <= 1e-10) }'; then
    ok "$tag: energy-grouped gradient = point-by-point gradient (${same% *} of the largest, ${same#* } parameters)"
  else
    bad "$tag: energy-grouped vs point-by-point gradient: $same"
  fi
done

echo
if [ "$failures" -eq 0 ]; then echo "angle_mapping: all checks passed"; else echo "angle_mapping: $failures check(s) FAILED"; fi
[ "$failures" -eq 0 ]
