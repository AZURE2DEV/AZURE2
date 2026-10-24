#!/usr/bin/env bash
#
# "total 1" on a segment line: the total cross section of a neutral entrance pair by the optical theorem,
# sigma_tot = (pi/k^2) g sum_J (2J+1) sum_c 2 Re T_cc (GenMatrixFunc::CalculateBareCrossSection), one T-matrix per
# point instead of one per exit pair of a SUM.  Built at run time from tests/natural_neon (n+20Ne, n+22Ne, n+21Ne;
# three compound nuclei, natural-abundance sums, target integration) with a free 1 keV width given to every open
# reaction channel (alpha+17O, n'+20Ne, n'+22Ne, n'+21Ne), which are all zero in the fitted model:
#
#   sum    each total as the SUM of every exit pair of its compound nucleus (the suite's own lines);
#   total  each total as "total 1" (the natural sums over the three nuclei stay SUMs of totals);
#   elast  the same lines without "total 1": elastic only.
#
#   (a) every cross section of `total` equals that of `sum` to 1e-9, data (mode 1) and extrapolation (mode 3), and
#       `elast` differs (the reaction channels matter);
#   (b) on segment 1 (Vaughn, 0-2.33 MeV, inelastic channels open), the reaction widths free: a fit with the
#       analytic gradient starts from the chi2 of the calculation, and the gradient agrees with finite differences
#       and with the Jacobian path (AZURE_GRAD_CHECK);
#   (c) refused, with an ERROR and no output: "total 1" on a charged pair, with an exit pair other than the entrance
#       pair, and on a differential segment.
#
#   ./tests/total_optical/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$HERE/../natural_neon"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
# macOS has no timeout(1); run without a limit there
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/total_optical.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# 1 keV (alternating sign, free) on every reaction-channel line of a level above that channel's threshold; the test
# segments up to 2.5 MeV, where the inelastic channels are open
awk '/<levels>/ { f = 1; print; next } /<\/levels>/ { f = 0 }
     f && NF > 30 && ($6 == 1 || $6 == 3 || $6 == 5 || $6 == 7) && $3 > $23 + $26 + 0.05 { $11 = 0; $12 = (k++ % 2 ? -1000 : 1000) }
     /<segmentsTest>/ { x = 1; print; next } /<\/segmentsTest>/ { x = 0 }
     x && NF > 10 { $5 = 2.5; $6 = 0.01 }
     { print }' "$REF/natural_neon.azr" > "$WORK/sum.azr"
# MODE tot|ela: each data and test line rewritten without the reaction components; "natural" sums keep the three totals
rewrite() {
  awk -v mode="$1" '
    function rw(line, keep,   n, a, i, head, nat, t) {
      n = split(line, a, " "); head = ""
      for (i = 1; i <= keep; i++) head = head (i > 1 ? " " : "") a[i]
      nat = (line ~ /natural 1/) ? " natural 1" : ""
      t = (mode == "tot") ? " total 1" : ""
      if (line ~ / 4 4 -999/) return head " 1 0 2 4 4 -999 6 6 -999" (keep == 15 ? " 0" : "") nat t
      return head (keep == 15 ? " 0 0" : "") nat t
    }
    /<segmentsData>/ { d = 1; print; next } /<\/segmentsData>/ { d = 0 }
    /<segmentsTest>/ { x = 1; print; next } /<\/segmentsTest>/ { x = 0 }
    d && NF > 15 { print rw($0, 15); next }
    x && NF > 10 { print rw($0, 10); next }
    { print }' "$WORK/sum.azr"
}
rewrite tot > "$WORK/total.azr"
rewrite ela > "$WORK/elast.azr"
[ "$(grep -c ' total 1' "$WORK/total.azr")" = 11 ] || { bad "models not built ($(grep -c ' total 1' "$WORK/total.azr") total lines)"; exit 1; }
# The extrapolation copies: target effects off (they would fold the test segments, see the azure2-eval skill).
for m in sum total elast; do
  awk '/<targetInt>/ { f = 1; print; next } /<\/targetInt>/ { f = 0 } f { $1 = 0 } { print }' "$WORK/$m.azr" > "$WORK/x$m.azr"
done

# run NAME MODEL MODE -> $WORK/NAME
run() {
  local dir="$WORK/$1" input
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$REF/data" "$dir/data"
  cp "$WORK/$2" "$dir/t.azr"
  if [ "$3" = 3 ]; then input='3\n\n\n7\n'; else input='1\n\n\n7\n'; fi
  (cd "$dir" && printf "$input" | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}
# max relative difference of the calculated cross section (4th column) over every AZUREOut file of two runs
maxrel() {
  local a="$1" b="$2" f n=0 w=0 r
  for f in "$a"/output/AZUREOut_aa=*."$3"; do
    [ -e "$f" ] || continue
    [ -e "$b/output/$(basename "$f")" ] || { echo "missing"; return; }
    r="$(awk 'NR == FNR { if (NF >= 4) v[FNR] = $4; next }
              NF >= 4 { d = $4 - v[FNR]; if (d < 0) d = -d; s = v[FNR] < 0 ? -v[FNR] : v[FNR]; if (s > 0 && d / s > w) w = d / s; n++ }
              END { printf "%d %.3e", n, w }' "$f" "$b/output/$(basename "$f")")"
    n=$((n + ${r%% *})); w="$(awk -v a="$w" -v b="${r##* }" 'BEGIN { print (b > a ? b : a) }')"
  done
  echo "$n $w"
}

for m in sum total elast; do run "$m" "$m.azr" 1; run "x$m" "x$m.azr" 3; done
for kind in "data out ''" "extrapolation extrap x"; do
  set -- $kind; what="$1"; ext="$2"; p="${3//\'/}"
  read -r n w <<< "$(maxrel "$WORK/${p}sum" "$WORK/${p}total" "$ext")"
  if [ "$n" = missing ] || [ -z "$n" ] || [ "$n" -lt 100 ]; then bad "(a) $what: outputs missing ($n points)"
  elif awk -v w="$w" 'BEGIN { exit !(w <= 1e-9) }'; then ok "(a) $what: total = sum of every exit pair at $n points (max rel $w)"
  else bad "(a) $what: total vs sum of every exit pair: max rel $w at $n points"; fi
  read -r n w <<< "$(maxrel "$WORK/${p}sum" "$WORK/${p}elast" "$ext")"
  if awk -v w="$w" 'BEGIN { exit !(w > 1e-3) }'; then ok "(a) $what: elastic alone differs (max rel $w): the reaction channels matter"
  else bad "(a) $what: elastic alone is the sum (max rel $w): the test does not discriminate"; fi
done

# (b) the gradient check of tests/13N on the total model with only segment 1 active (Vaughn, 0-2.33 MeV: the natural
# sum of the three totals, inelastic channels open, no target integration), the reaction widths free.  (The
# target-integrated Sikkema segments disagree with finite differences by ~1e-2 on a few parameters with or without
# "total 1", from the integration grid.)
awk '/<segmentsData>/ { f = 1; n = 0; print; next } /<\/segmentsData>/ { f = 0 }
     f { n++; if (n != 1) $1 = 0 } { print }' "$WORK/total.azr" > "$WORK/gtotal.azr"
run gcalc gtotal.azr 1
mkdir -p "$WORK/fit/output" "$WORK/fit/checks"; ln -s "$REF/data" "$WORK/fit/data"; cp "$WORK/gtotal.azr" "$WORK/fit/t.azr"
(cd "$WORK/fit" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
   run_limited 1200 "$AZURE2_BIN" --use-gradient --no-gui --no-readline t.azr > run.log 2>&1)
chi="$(awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { sub(",", "", $3); print $3 }' "$WORK/fit/run.log")"
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/gcalc/output/chiSquared.out")"
if [ -n "$chi" ] && awk -v a="$chi" -v b="$want" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(d <= 1e-5) }'; then
  ok "(b) chi2 at the start of the fit $chi = calculation $want"
else bad "(b) chi2 at the start of the fit '$chi', calculation $want"; fi
awk '{ sub(/\r$/, "") }
     /^\[grad-check\] idx / { a = f = j = ""
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print a, f, j }' "$WORK/fit/run.log" > "$WORK/fit/rows"
nrow="$(awk 'END { print NR }' "$WORK/fit/rows")"
if [ "$nrow" -lt 5 ]; then bad "(b) gradient check listed $nrow parameters"
else
  read -r fd jac <<< "$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a
                               d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d
                               e = $1 - $3; if (e < 0) e = -e; if (e > v) v = e }
                             END { printf "%.2e %.2e", w / m, v / m }' "$WORK/fit/rows")"
  awk -v w="$fd" 'BEGIN { exit !(w + 0 <= 1e-5) }' && ok "(b) gradient = finite differences for $nrow parameters ($fd of the largest)" \
    || bad "(b) gradient vs finite differences: $fd of the largest"
  awk -v w="$jac" 'BEGIN { exit !(w + 0 <= 1e-10) }' && ok "(b) Jacobian-path gradient = adjoint gradient ($jac of the largest)" \
    || bad "(b) Jacobian-path vs adjoint gradient: $jac of the largest"
fi

# (c) refusals: segment 3 (n+21Ne, 6 -> 6) rewritten
refuse() {
  local name="$1" edit="$2"
  awk -v edit="$edit" '/<segmentsData>/ { f = 1; n = 0; print; next } /<\/segmentsData>/ { f = 0 }
       f { n++; if (n == 3) { split(edit, e, ":"); $(e[1]) = e[2] } } { print }' "$WORK/total.azr" > "$WORK/$name.azr"
  run "$name" "$name.azr" 1
  if grep -q "ERROR: segment #3 (total 1)" "$WORK/$name/run.log" && [ ! -s "$WORK/$name/output/chiSquared.out" ]; then
    ok "(c) refused: $name"
  else bad "(c) not refused: $name"; fi
}
refuse charged "2:1"          # entrance alpha+17O (and exit 6): charged pair
refuse exit "3:7"             # n+21Ne -> n'+21Ne
refuse differential "8:1"     # data type 1

if [ "$failures" -eq 0 ]; then echo "PASS: total cross section (optical theorem)"; else echo "FAIL: $failures check(s)"; exit 1; fi
