#!/usr/bin/env bash
#
# Point mapping inside segments with components (EData::MapComponentSegments): each component's point is mapped onto
# the base segment's point at the same energy, and the base's calculation (one level-matrix solution) builds the
# component's T-matrix and observable (EPoint::Calculate).  Before 2026-10-24 every term was its own calculation.
# The 13N reference model is edited at run time into composites of both exit pairs (p+12C elastic, 1; capture to
# 13N g.s. with external capture, 2):
#
#   seg 2  Seagrave capture, gas-target integration: SUM + 1e-7 x elastic (angle-integrated) + 0.5 x capture
#   seg 3  Meyer 84.3 deg elastic, fixed energy shift 2 keV: SUM + 1000 x capture at the same angle
#   seg 4  Meyer 114.5 deg elastic: RATIO elastic / capture
#   seg 9  Kettner 0 deg capture, sub-point target effect: SUM + capture at 90 deg (same exit, another angle)
#   seg 10 Kettner 55 deg capture: SUM + 1e-7 x elastic at 55 deg
#   test   0.3-3 MeV, 90 deg: SUM elastic + 1000 x capture (Calculate Segments Without Data)
#
#   (a) every output file identical with and without the component mapping (AZURE_NO_COMPONENT_MAPPING=1), for the
#       data (mode 1) and the extrapolation (mode 3);
#   (b) the mapping is used: fewer level-matrix solutions (AZURE_PROFILE);
#   (c) without the ratio: a fit with the analytic gradient starts from the chi2 of the calculation, the gradient
#       (one level-matrix solution for all terms of a composite point, AZUREGrad GradOnePoint) agrees with finite
#       differences of the mapped evaluation, and the Jacobian path with the adjoint (AZURE_GRAD_CHECK).
#
#   ./tests/component_mapping/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# macOS has no timeout(1); run without a limit there
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/component_mapping.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

awk '/<segmentsData>/ { f = 1; n = 0; print; next }
     /<\/segmentsData>/ { f = 0 }
     /<segmentsTest>/ { print; print "1 1 1 0.3 3.0 0.01 90 90 0 1 1 0 -1 2 1 1 -999 1 1 2 -999 1000"; next }
     f { n++
         tail = ""
         if (n == 2)  tail = "1 0 -1 2 1 1 -999 1e-7 1 2 -999 0.5 0"
         if (n == 3)  { $12 = 0.002; tail = "1 0 -1 1 1 2 -999 1000 0" }
         if (n == 4)  tail = "1 1 1 1 2 -999 0"
         if (n == 9)  tail = "1 0 -1 1 1 2 90 1 0"
         if (n == 10) tail = "1 0 -1 1 1 1 -999 1e-7 0"
         if (tail != "") { NF = 15; $0 = $0 " " tail } }
     { print }' "$REF/13N.azr" > "$WORK/t.azr"
[ "$(awk '/<segmentsData>/ { f = 1; next } /<\/segmentsData>/ { f = 0 } f && NF > 17' "$WORK/t.azr" | awk 'END { print NR }')" = 5 ] \
  || { bad "model not built"; exit 1; }
# The fit copy (c): without the ratio, whose elastic/capture values (~1e8 against data of ~1) would dominate chi2 and
# hide the other composites in the gradient check.
awk '/<segmentsData>/ { f = 1; n = 0; print; next } /<\/segmentsData>/ { f = 0 }
     f { n++; if (n == 4) { NF = 15; $0 = $0 " 0 0" } } { print }' "$WORK/t.azr" > "$WORK/g.azr"
# The extrapolation copy: target effects off (they would fold the test segment, see the azure2-eval skill).
awk '/<targetInt>/ { f = 1; print; next } /<\/targetInt>/ { f = 0 } f { $1 = 0 } { print }' "$WORK/t.azr" > "$WORK/x.azr"

# run NAME MODEL MODE [ENV...] -> $WORK/NAME
run() {
  local name="$1" model="$2" mode="$3"; shift 3
  local dir="$WORK/$name" input
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$REF/data" "$dir/data"
  cp "$WORK/$model" "$dir/t.azr"
  if [ "$mode" = 3 ]; then input='3\n\n\n7\n'; else input='1\n\n\n7\n'; fi
  (cd "$dir" && for kv in "$@" AZURE_PROFILE=1; do export "$kv"; done &&
     printf "$input" | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  if [ -z "$(ls "$dir/output" | grep -i 'AZUREOut')" ]; then bad "$name: no output written"; tail -15 "$dir/run.log"; return 1; fi
}
# awk only: the Windows runner has no cmp/diff
same_outputs() {
  local a="$1" b="$2" n=0 d=0 f
  for f in "$a"/output/*.out "$a"/output/*.extrap; do
    [ -e "$f" ] || continue
    n=$((n + 1))
    if [ ! -e "$b/output/$(basename "$f")" ]; then d=$((d + 1)); continue; fi
    awk 'NR == FNR { x[FNR] = $0; m = FNR; next } { k = FNR; if (x[FNR] != $0) e++ } END { exit !(e == 0 && k == m) }' \
      "$f" "$b/output/$(basename "$f")" || d=$((d + 1))
  done
  echo "$n $d"
}
fills() { awk '/^\[profile\]   matrix fill/ { print $5 }' "$1/run.log"; }

run on  t.azr 1 X=0 || exit 1
run off t.azr 1 AZURE_NO_COMPONENT_MAPPING=1 || exit 1
read -r n d <<< "$(same_outputs "$WORK/on" "$WORK/off")"
if [ "$n" -ge 3 ] && [ "$d" = 0 ]; then ok "(a) data: all $n output files identical with and without the mapping"
else bad "(a) data: $d of $n output files differ"; fi
fon="$(fills "$WORK/on")"; foff="$(fills "$WORK/off")"
if [ -n "$fon" ] && [ -n "$foff" ] && [ "$fon" -lt "$foff" ]; then ok "(b) level-matrix solutions $foff -> $fon"
else bad "(b) level-matrix solutions '$foff' -> '$fon'"; fi

run xon  x.azr 3 X=0 || exit 1
run xoff x.azr 3 AZURE_NO_COMPONENT_MAPPING=1 || exit 1
read -r n d <<< "$(same_outputs "$WORK/xon" "$WORK/xoff")"
if [ "$n" -ge 1 ] && [ "$d" = 0 ]; then ok "(a) extrapolation: all $n output files identical"
else bad "(a) extrapolation: $d of $n output files differ"; fi
fon="$(fills "$WORK/xon")"; foff="$(fills "$WORK/xoff")"
if [ -n "$fon" ] && [ -n "$foff" ] && [ "$fon" -lt "$foff" ]; then ok "(b) extrapolation: level-matrix solutions $foff -> $fon"
else bad "(b) extrapolation: level-matrix solutions '$foff' -> '$fon'"; fi

# (c) the gradient check of tests/13N on the composite model
mkdir -p "$WORK/fit/output" "$WORK/fit/checks"; ln -s "$REF/data" "$WORK/fit/data"; cp "$WORK/g.azr" "$WORK/fit/t.azr"
run gcalc g.azr 1 X=0 || exit 1
(cd "$WORK/fit" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
   run_limited 1200 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
chi="$(awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { sub(",", "", $3); print $3 }' "$WORK/fit/run.log")"
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/gcalc/output/chiSquared.out")"
if [ -n "$chi" ] && awk -v a="$chi" -v b="$want" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(d <= 1e-5) }'; then
  ok "(c) chi2 at the start of the fit $chi = calculation $want"
else bad "(c) chi2 at the start of the fit '$chi', calculation $want"; fi
awk '{ sub(/\r$/, "") }
     /^\[grad-check\] idx / { a = f = j = ""
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print a, f, j }' "$WORK/fit/run.log" > "$WORK/fit/rows"
nrow="$(awk 'END { print NR }' "$WORK/fit/rows")"
if [ "$nrow" -lt 5 ]; then bad "(c) gradient check listed $nrow parameters"
else
  fd="$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a; d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d }
             END { printf "%.2e", w / m }' "$WORK/fit/rows")"
  awk -v w="$fd" 'BEGIN { exit !(w + 0 <= 1e-5) }' && ok "(c) gradient = finite differences for $nrow parameters ($fd of the largest)" \
    || bad "(c) gradient vs finite differences: $fd of the largest"
  jac="$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a; e = $1 - $3; if (e < 0) e = -e; if (e > v) v = e }
              END { printf "%.2e", v / m }' "$WORK/fit/rows")"
  awk -v w="$jac" 'BEGIN { exit !(w + 0 <= 1e-10) }' && ok "(c) Jacobian-path gradient = adjoint gradient ($jac of the largest)" \
    || bad "(c) Jacobian-path vs adjoint gradient: $jac of the largest"
fi

if [ "$failures" -eq 0 ]; then echo "PASS: component mapping"; else echo "FAIL: $failures check(s)"; exit 1; fi
