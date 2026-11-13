#!/usr/bin/env bash
#
# A level's radiative decay written as Gamma_gamma, branchings and mixing ratios (<gammaDecays>, gamma_decay;
# docs/upgrade/gamma_parametrisation.md), against the same level written with partial widths:
#   (a) identical calculations: models/partial.azr (the 14N(p,g)15O model, level 5 the 259.7 keV resonance with five
#       radiative widths) and the same project with the record made by hand (models/make_record.py, the note's
#       definitions, not the engine's map), in the Brune and Park bases; the formal basis and external capture on a
#       record channel are refused;
#   (b) the sign of the mixing ratio: a 3/2+ level decaying by M1 and E2 (models/sign_partial.azr, E2 width written
#       negative) and its record (delta = +0.5): the same angular distribution, whose a2/a0 is Krane-Steffen's for
#       delta = +0.5 (the note's section 2.3, case a);
#   (c) fits of the five data sets in both parametrisations, MIGRAD from displaced values and Levenberg-Marquardt:
#       the same chi^2 minimum; a prior on one branching: the pull the prior's formula gives;
#   (d) the gradient in the record's numbers, the level's energy and a proton amplitude against central differences
#       (AZURE_GRAD_CHECK, the gradient and the Jacobian paths), Brune and Park;
#   (e) the boundary: a branching and the remainder driven below 0 end at 0 without NaN, with a note.
#
#   ./tests/gamma_decay/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
M="$HERE/models"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/edit_guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }
PY="${PYTHON:-python3}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/gamma_decay.XXXXXX")"
trap '[ -n "${KEEP_WORK:-}" ] || rm -rf "$WORK"' EXIT  # KEEP_WORK=1 keeps the runs

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
close() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; exit !(a != "" && d <= t * m) }'; }
check() {  # check LABEL VALUE REFERENCE TOL
  if close "$2" "$3" "$4"; then ok "$1: $2 (reference $3)"; else bad "$1: ${2:-missing}, reference $3 (rel. tol. $4)"; fi
}

# project NAME FILE: $WORK/NAME/p.azr = FILE, with the data
project() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  ln -s "$M/data" "$WORK/$1/data"
  cp "$2" "$WORK/$1/p.azr"
}
# run NAME MODE [FLAGS...]: menu MODE (1 calculate, 2 fit; 3 extrapolate) on $WORK/NAME/p.azr; exit status in NAME/exit
run() {
  local name="$1" mode="$2"
  shift 2
  (cd "$WORK/$name" && printf '%s\n\n\n' "$mode" | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline "$@" p.azr > run.log 2>&1
   echo $? > exit)
}
chi2() { sed -n 's/^Total-Chi-Squared: \([^ ]*\).*/\1/p' "$WORK/$1/output/chiSquared.out"; }
sav() { awk -v n="$2" -v c="${3:-2}" '$1 == n { print $c; exit }' "$WORK/$1/output/param.sav"; }
# the parameters.out lines before the record's block (the level parameters of both parametrisations)
levelpart() { awk '/^RADIATIVE DECAYS/ { exit } { print }' "$WORK/$1/output/parameters.out"; }

"$PY" "$M/make_record.py" "$M/partial.azr" "$WORK/record.azr" 5 5.183 fit
if cmp -s "$WORK/record.azr" "$M/record.azr"; then ok "models/record.azr is make_record.py's record of partial.azr"
else bad "models/record.azr differs from make_record.py's record of partial.azr (regenerate it)"; fi

echo "-- (a) the same project written with partial widths and with the record: identical calculations"
for b in "" --use-park; do
  project "partial$b" "$M/partial.azr"
  project "record$b" "$M/record.azr"
  must_differ "$WORK/partial$b/p.azr" "$WORK/record$b/p.azr" || bad "the record project is the partial-width one"
  run "partial$b" 1 $b
  run "record$b" 1 $b
  n=0; d=0
  for f in "$WORK/partial$b"/output/*; do
    base="$(basename "$f")"
    case "$base" in param.par|parameters.out) continue ;; esac
    n=$((n + 1))
    cmp -s "$f" "$WORK/record$b/output/$base" || { d=$((d + 1)); echo "        differs: $base"; }
  done
  if [ "$n" -ge 6 ] && [ "$d" -eq 0 ]; then ok "${b:---use-brune}: $n output files identical (the five data sets, chi^2)"
  else bad "${b:---use-brune}: $d of $n output files differ"; fi
  if [ "$(levelpart "partial$b" | md5sum)" = "$(levelpart "record$b" | md5sum)" ]; then
    ok "${b:---use-brune}: parameters.out's level parameters identical"
  else bad "${b:---use-brune}: parameters.out's level parameters differ"; fi
done
grep -q '^RADIATIVE DECAYS' "$WORK/record/output/parameters.out" && grep -q 'g+15O\[5.183\]: b = 0.168819 (remainder)' "$WORK/record/output/parameters.out" \
  && ok "parameters.out: the record's block with the remainder b = 0.168819 (0.0065 / 0.0385027 eV)" \
  || bad "parameters.out: no record block, or the remainder is not 0.168819"

project formal "$M/record.azr"
run formal 1 --no-brune
if grep -q 'ERROR: <gammaDecays> record .level5. (level 5): the formal basis' "$WORK/formal/run.log" && [ ! -s "$WORK/formal/output/chiSquared.out" ]; then
  ok "--no-brune: refused with the message, nothing calculated"
else bad "--no-brune: not refused"; fi
project external "$M/record.azr"
# external capture into the ground state of 15O (pair 2) by E1: column 31 of its rows
edited "$M/record.azr" "$WORK/external/p.azr" awk '{ n = split($0, a, " "); if (n > 30 && a[6] == "2" && a[27] == "10") { $31 = 1; print " " $0 } else print }' \
  || bad "external-capture variant"
run external 1
if grep -q 'has external capture by this multipolarity' "$WORK/external/run.log"; then ok "external capture on a record channel: refused"
else bad "external capture on a record channel: not refused"; fi

echo "-- (b) the sign of the mixing ratio (Krane-Steffen): M1/E2 of a 3/2+ level, the E2 width written negative"
"$PY" "$M/make_record.py" "$M/sign_partial.azr" "$WORK/sign_record.azr" 1 5.183 fixed
grep -q 'mixing       g+15O\[5.183\] 0.5 fixed' "$WORK/sign_record.azr" && ok "make_record.py: delta = +0.5 (widths of opposite signs)" \
  || bad "make_record.py: delta is not +0.5"
for v in partial record; do
  mkdir -p "$WORK/sign_$v/output" "$WORK/sign_$v/checks"
  [ "$v" = partial ] && cp "$M/sign_partial.azr" "$WORK/sign_$v/p.azr" || cp "$WORK/sign_record.azr" "$WORK/sign_$v/p.azr"
  run "sign_$v" 3
done
ext="AZUREOut_aa=1_R=2.extrap"
if [ -s "$WORK/sign_partial/output/$ext" ] && cmp -s "$WORK/sign_partial/output/$ext" "$WORK/sign_record/output/$ext"; then
  ok "dsigma/dOmega at 19 angles identical"
else bad "dsigma/dOmega differs (or was not written)"; fi
a2="$("$PY" -c "
import numpy as np, sys
t = np.loadtxt(sys.argv[1]); c = np.polynomial.legendre.legfit(np.cos(np.radians(t[:, 2])), t[:, 3], 4)
print('%.6f' % (c[2] / c[0]))" "$WORK/sign_record/output/$ext" 2>/dev/null)"
# Krane, Steffen and Wheeler (the note's table, case a): a2/a0 = -0.99282 for delta = +0.5, +0.39282 for -0.5
check "a2/a0 of the record with delta = +0.5 is Krane-Steffen's" "$a2" -0.99282 1e-4

echo "-- (c) fits in both parametrisations"
# displaced starts: level 5's radiative widths x 1.3 / Gamma_gamma x 1.3 and two branchings moved
edited "$M/partial.azr" "$WORK/partial_d.azr" awk '{ n = split($0, a, " "); if (n > 30 && a[9] == "5" && a[27] == "10" && a[12] != "0") { $12 = a[12] * 1.3; print " " $0 } else print }' \
  || bad "displaced partial-width start"
edited "$M/record.azr" "$WORK/record_d.azr" sed -e 's/^\(    total_width  \)[^ ]*/\10.05/' -e 's/^\(    branching    g+15O\[6.7931\] \)[^ ]*/\10.3/' \
  -e 's/^\(    branching    g+15O\[6.1763\] \)[^ ]*/\10.5/' || bad "displaced record start"
project fit_partial "$WORK/partial_d.azr"
project fit_record "$WORK/record_d.azr"
project fit_record_lm "$WORK/record_d.azr"
run fit_partial 2
run fit_record 2
run fit_record_lm 2 --use-lm
c0="$(chi2 fit_partial)"
check "MIGRAD: the record's chi^2 minimum equals the partial widths'" "$(chi2 fit_record)" "$c0" 1e-5
check "Levenberg-Marquardt: the same minimum" "$(chi2 fit_record_lm)" "$c0" 1e-5
b6="$(sav fit_record 'level5.branching.g+15O[6.7931]')"; s6="$(sav fit_record 'level5.branching.g+15O[6.7931]' 3)"
check "MIGRAD and LM: the same b(6.7931)" "$(sav fit_record_lm 'level5.branching.g+15O[6.7931]')" "$b6" 1e-4
# A prior two errors above the fitted b(6.7931), as wide as its error: to the linearisation the fit moves halfway,
# b' = b + s, and the pull (b' - c) / w is -1.
c="$(awk -v b="$b6" -v s="$s6" 'BEGIN { printf "%.10g", b + 2 * s }')"
edited "$M/record.azr" "$WORK/record_prior.azr" sed "s/^\(    branching    g+15O\[6.7931\] [^ ]* fit\)$/\1 prior $c $s6/" || bad "prior variant"
project fit_prior "$WORK/record_prior.azr"
run fit_prior 2
bp="$(sav fit_prior 'level5.branching.g+15O[6.7931]')"
pull="$(awk -v b="$bp" -v c="$c" -v w="$s6" 'BEGIN { printf "%.6f", (b - c) / w }')"
check "a prior (b + 2 sigma, sigma) on b(6.7931): the pull" "$pull" -1 0.02

echo "-- (d) the gradient against central differences (AZURE_GRAD_CHECK)"
# level 5's energy and its first proton width free as well: their gradients take the photon amplitudes' share
# and priors on a branching and on the remainder, whose gradients go straight to the branchings
edited "$M/record.azr" "$WORK/record_free0.azr" awk '{ n = split($0, a, " "); if (n > 30 && a[9] == "5") { $4 = 0; if (a[27] == "0" && a[7] == "1") $11 = 0; print " " $0 } else print }' \
  || bad "free-energy variant"
edited "$WORK/record_free0.azr" "$WORK/record_free.azr" sed -e 's/^\(    branching    g+15O\[6.1763\] [^ ]* fit\)$/\1 prior 0.6 0.01/' \
  -e 's/^\(    remainder    g+15O\[5.183\]\)$/\1 prior 0.15 0.01/' || bad "prior variant of the gradient check"
for b in "" --use-park; do
  project "grad$b" "$WORK/record_free.azr"
  run "grad$b" 1 $b  # param.par: the index of the record's first number
  i="$(awk '$1 == "level5.total_width" { print NR - 2; exit }' "$WORK/grad$b/output/param.par")"
  (cd "$WORK/grad$b" && printf '2\n\n\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_IDX="$((i - 2))-$((i + 7))" \
     run_limited 900 "$AZURE2_BIN" --no-gui --no-readline --use-gradient $b p.azr > grad.log 2>&1)
  # the energy with the step 1e-8 (|x| + 1): the default 1e-6 is 8.5 eV in E_x, a truncation error of 1.5e-3 across
  # the 1 keV resonance (the same with partial widths)
  (cd "$WORK/grad$b" && printf '2\n\n\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_H=1e-8 AZURE_GRAD_CHECK_IDX="$((i - 3))" \
     run_limited 900 "$AZURE2_BIN" --no-gui --no-readline --use-gradient $b p.azr > gradE.log 2>&1)
  rows="$(grep -c '\[grad-check\] idx' "$WORK/grad$b/grad.log")"
  worst="$(awk '/\[grad-check\] idx/ && $4 != "energy" { for (k = 1; k <= NF; k++) if ($k == "rel" || $k == "relJ") if ($(k + 1) > m) m = $(k + 1) } END { print m + 0 }' "$WORK/grad$b/grad.log")"
  wE="$(awk '/\[grad-check\] idx/ && $4 == "energy" { m = 0; for (k = 1; k <= NF; k++) if ($k == "rel" || $k == "relJ") if ($(k + 1) > m) m = $(k + 1) } END { print m + 0 }' "$WORK/grad$b/gradE.log")"
  kinds="$(awk '/\[grad-check\] idx/ { print $4 }' "$WORK/grad$b/grad.log" | sort -u | tr '\n' ' ')"
  if [ "$rows" -ge 7 ] && awk -v w="$worst" 'BEGIN { exit !(w < 1e-6) }' && [ "$kinds" = "energy gamma radiative " ]; then
    ok "${b:---use-brune}: $rows parameters ($kinds), gradient and Jacobian paths, largest rel. difference $worst"
  else bad "${b:---use-brune}: $rows parameters ($kinds), largest rel. difference $worst (limit 1e-6)"; fi
  # the energy: the photon amplitudes' share is 5e-4 of its gradient here (Gamma_gamma held while E moves)
  if awk -v w="$wE" 'BEGIN { exit !(w > 0 && w < 1e-5) }'; then ok "${b:---use-brune}: level energy, rel. difference $wE"
  else bad "${b:---use-brune}: level energy, rel. difference ${wE:-missing} (limit 1e-5)"; fi
done

echo "-- (f) the record's fit at 1 and at 3 threads: every output file identical, fitStatus.out's call count too"
# The count is the objective's own (AZURECalc::CountCallsInto): Minuit2's NFcn() loses increments under its OpenMP
# numerical gradient, which this fit takes, and came out 65 or 66 at 3 threads from run to run.
for t in 1 3; do
  project "thr$t" "$WORK/record_d.azr"
  OMP_NUM_THREADS=$t run "thr$t" 2 --use-gradient
done
n=0; d=0
for f in "$WORK/thr1"/output/*; do
  n=$((n + 1)); cmp -s "$f" "$WORK/thr3/output/$(basename "$f")" || { d=$((d + 1)); echo "        differs: $(basename "$f")"; }
done
if [ "$n" -ge 9 ] && [ "$d" -eq 0 ] && [ -s "$WORK/thr1/output/fitStatus.out" ]; then
  ok "$n output files identical at 1 and 3 threads ($(awk '/^calls/ { print $2 " calls" }' "$WORK/thr1/output/fitStatus.out"))"
else bad "$d of $n output files differ between 1 and 3 threads"; fi

echo "-- (e) the boundary"
# Without the data sets of 6.8594 and 5.183 MeV, priors below 0 on b(6.8594) and on the remainder (5.183 MeV).
edited "$M/record.azr" "$WORK/record_zero.azr" sed -e 's/^\(    branching    g+15O\[6.8594\] [^ ]* fit\)$/\1 prior -0.01 0.01/' \
  -e 's/^\(    remainder    g+15O\[5.183\]\)$/\1 prior -0.2 0.05/' -e 's/^dataset g6859$/&\n    active         no/' \
  -e 's/^dataset g5183$/&\n    active         no/' || bad "boundary variant"
project zero "$WORK/record_zero.azr"
run zero 2
if [ "$(cat "$WORK/zero/exit")" = 0 ] && ! grep -qiw 'nan' "$WORK/zero/output/parameters.out" "$WORK/zero/output/chiSquared.out" "$WORK/zero/output/param.sav"; then
  ok "a fit pulled below the bounds ends without NaN"
else bad "the fit pulled below the bounds failed or wrote NaN"; fi
grep -q 'note: the branching to g+15O\[6.8594\] is at its bound 0' "$WORK/zero/output/parameters.out" \
  && ok "b(6.8594) ends at its bound 0, with the note" || bad "b(6.8594) not at its bound, or no note"
r="$(awk '/g\+15O\[5.183\]: b = / { sub(/;.*/, ""); print $4; exit }' "$WORK/zero/output/parameters.out")"
if awk -v r="$r" 'BEGIN { exit !(r != "" && r < 2e-3 && r > -2e-3) }'; then ok "the remainder ends at 0 (within the wall's 1e-3): $r"
else bad "the remainder is ${r:-missing}, not at 0"; fi

sum="$(awk '$1 ~ /^level5\.branching\./ { s += $2 } END { printf "%.6f", s }' "$WORK/zero/output/param.sav")"
if awk -v s="$sum" 'BEGIN { exit !(s > 0.99 && s < 1.002) }'; then ok "the wall holds the explicit branchings' sum at 1: $sum"
else bad "the explicit branchings sum to ${sum:-missing}, not 1 (no wall?)"; fi

echo
if [ "$failures" -eq 0 ]; then echo "gamma_decay: all checks passed"; else echo "gamma_decay: $failures check(s) FAILED"; fi
exit $((failures > 0))
