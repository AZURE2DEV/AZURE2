#!/usr/bin/env bash
#
# A level's radiative decay written as Gamma_gamma, branchings and mixing ratios (<gammaDecays>, gamma_decay;
# docs/upgrade/gamma_parametrisation.md), against the same level written with partial widths:
#   (a) identical calculations: models/partial.azr (the 14N(p,g)15O model, level 5 the 259.7 keV resonance with five
#       radiative widths) and the same project with the record made by hand (models/make_record.py, the note's
#       definitions, not the engine's map), in the Brune and Park bases; the formal basis is refused;
#   (b) the sign of the mixing ratio: a 3/2+ level decaying by M1 and E2 (models/sign_partial.azr, E2 width written
#       negative) and its record (delta = +0.5): the same angular distribution, whose a2/a0 is Krane-Steffen's for
#       delta = +0.5 (the note's section 2.3, case a);
#   (c) fits of the five data sets in both parametrisations, MIGRAD from displaced values and Levenberg-Marquardt:
#       the same chi^2 minimum; a prior on one branching: the pull the prior's formula gives;
#   (d) the gradient in the record's numbers, the level's energy and a proton amplitude against central differences
#       (AZURE_GRAD_CHECK, the gradient and the Jacobian paths), Brune and Park;
#   (e) the boundary: a branching and the remainder driven below 0 end at 0 without NaN, with a note;
#   (h) external capture on the record's channels (the twin: models/partial.azr with capture into the ground state
#       and 6.1763 MeV by E1 and 6.7931 MeV by M1): identical calculations in both bases, also after moving the
#       level's proton width, its energy and the ground state's ANC; the record's widths are TransformOut's; the
#       gradient in every free parameter (the ANCs of the final states included) at two points; fits from displaced
#       starts (MIGRAD with the gradient, Levenberg-Marquardt) with the twin's minimum and branchings; 1 against 3
#       threads; a width below the floor 2P (Im gamma^ext)^2 / n as TransformIn handles it;
#   (i) the vector of observed parameters (note section 2.8): the grid copies of a narrow resonance (rate and fit) see
#       the twin's widths; limits on a particle width and on Gamma_gamma; the observed-parameter MCMC moves with the
#       record's numbers (pyazr's calls: tests/pyazr/gamma_decay_test.py part 5);
#   (j) a final state of an external part above a second pair's threshold, its energy free: the gradient against
#       central differences.
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

echo "-- (g) a mixing ratio fitted from 0 to an angular distribution (models/sign_ad.azr: the sign case's dsigma/dOmega)"
# delta = 0 is regular: the L+1 amplitude is linear in delta there, so its derivative is not 0 (the width form
# p / (2 Gamma) dGamma/ddelta is 0/0 at delta = 0 and gave a gradient of exactly 0)
"$PY" "$M/make_record.py" "$M/sign_ad.azr" "$WORK/sad.azr" 1 5.183 fixed
for start in 0 0.3; do
  edited "$WORK/sad.azr" "$WORK/sad$start.azr" sed "s/^\(    mixing       g+15O\[5.183\]\) 0.5 fixed$/\1 $start fit/" || bad "mixing start $start"
  project "sad$start" "$WORK/sad$start.azr"
  (cd "$WORK/sad$start" && printf '2\n\n\n' | AZURE_GRAD_CHECK=exit run_limited 900 "$AZURE2_BIN" --no-gui --no-readline --use-gradient p.azr > grad.log 2>&1)
  line="$(grep '\[grad-check\] idx' "$WORK/sad$start/grad.log" | sed -n 1p)"
  rel="$(echo "$line" | awk '{ for (k = 1; k <= NF; k++) if ($k == "rel" || $k == "relJ") if ($(k + 1) > m) m = $(k + 1) } END { print m + 0 }')"
  an="$(echo "$line" | awk '{ for (k = 1; k <= NF; k++) if ($k == "analytic") print $(k + 1) }')"
  if [ -n "$line" ] && awk -v r="$rel" -v a="$an" 'BEGIN { exit !(r < 1e-6 && a != 0) }'; then
    ok "delta = $start: d chi^2 / d delta $an, gradient and Jacobian paths against central differences, rel. $rel"
  else bad "delta = $start: d chi^2 / d delta ${an:-missing}, rel. ${rel:-missing} (limit 1e-6)"; fi
done
for m in --use-gradient --use-lm; do
  project "sadfit$m" "$WORK/sad0.azr"
  run "sadfit$m" 2 $m
  check "$m from delta = 0: the fitted delta is the twin's +0.5" "$(sav "sadfit$m" 'level1.mixing.g+15O[5.183]')" 0.5 1e-4
done

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

echo "-- (h) external capture on the record's channels: ground state and 6.1763 MeV by E1, 6.7931 MeV by M1"
# The twin: models/partial.azr with capture into the ground state (pair 2) and 6.1763 MeV (pair 6) by E1 and into
# 6.7931 MeV (pair 3) by M1, column 31; the record made from it as models/record.azr is made from partial.azr.
ec() { awk '{ n = split($0, a, " "); if (n > 30 && a[27] == "10" && (a[6] == "2" || a[6] == "6")) { $31 = 1; print " " $0 }
  else if (n > 30 && a[27] == "10" && a[6] == "3") { $31 = 2; print " " $0 } else print }' "$1"; }
# moved: level 5's proton width 974.937 -> 1100 eV, its energy 7.5567 -> 7.5569 MeV, the ground state's ANC 0.23 -> 0.30
moved() { awk '{ n = split($0, a, " "); if (n > 30 && a[9] == "5") { $3 = 7.5569; if (a[27] == "0" && a[7] == "1") $12 = 1100; print " " $0 }
  else if (n > 30 && a[9] == "1" && a[27] == "0" && a[7] == "1") { $12 = 0.30; print " " $0 } else print }' "$1"; }
# free: level 5's energy and proton width, the ground state's first ANC and the two ANCs of 6.1763 MeV
free() { awk '{ n = split($0, a, " "); if (n > 30 && a[9] == "5" ) { $4 = 0; if (a[27] == "0" && a[7] == "1") $11 = 0; print " " $0 }
  else if (n > 30 && a[27] == "0" && ((a[9] == "1" && a[7] == "1") || (a[9] == "9" && a[12] != "0"))) { $11 = 0; print " " $0 }
  else print }' "$1"; }
# same_outputs A B LABEL: every output file but param.par and parameters.out identical, parameters.out's level part too
same_outputs() {
  local n=0 d=0 f base
  for f in "$WORK/$1"/output/*; do
    base="$(basename "$f")"
    case "$base" in param.par|parameters.out) continue ;; esac
    n=$((n + 1))
    cmp -s "$f" "$WORK/$2/output/$base" || { d=$((d + 1)); echo "        differs: $base"; }
  done
  if [ "$n" -ge 6 ] && [ "$d" -eq 0 ] && [ "$(levelpart "$1" | md5sum)" = "$(levelpart "$2" | md5sum)" ]; then
    ok "$3: $n output files and parameters.out's level parameters identical"
  else bad "$3: $d of $n output files differ, or parameters.out's level parameters"; fi
}
# twin_b NAME E_x: Gamma_gamma (eV) and the branchings per photon pair of the level at E_x, from parameters.out's
# observed widths (what TransformOut prints)
twin_b() { awk -v E="$2" '$1 == "J" && $6 == E { on = 1; next } on && NF == 0 { on = 0 }
  on && $4 == "L" { w = $12 * ($13 == "keV" ? 1e3 : $13 == "eV" ? 1 : $13 == "meV" ? 1e-3 : 0); g[$3] += w; t += w }
  END { printf "total %.10g\n", t; for (r in g) printf "R%s %.10g\n", r, g[r] / t }' "$WORK/$1/output/parameters.out"; }
# block_w NAME LABEL: the record block's first width (eV) of a final state; level_w NAME R [E_x]: TransformOut's (eV)
block_w() { awk -v l="$2:" '$1 == l { for (k = 1; k <= NF; k++) if ($k == "Gamma") { print $(k + 2); exit } }' "$WORK/$1/output/parameters.out"; }
level_w() { awk -v r="$2" -v E="${3:-7.5567}" '$1 == "J" && $6 == E { on = 1; next } on && NF == 0 { on = 0 }
  on && $3 == r && $4 == "L" { print $12 * ($13 == "keV" ? 1e3 : $13 == "eV" ? 1 : 1e-3); exit }' "$WORK/$1/output/parameters.out"; }

edited "$M/partial.azr" "$WORK/partial_ec.azr" ec || bad "external-capture twin"
"$PY" "$M/make_record.py" "$WORK/partial_ec.azr" "$WORK/record_ec.azr" 5 5.183 fit
for b in "" --use-park; do
  project "pec$b" "$WORK/partial_ec.azr"
  project "rec$b" "$WORK/record_ec.azr"
  run "pec$b" 1 $b
  run "rec$b" 1 $b
  same_outputs "pec$b" "rec$b" "${b:---use-brune}, external capture"
done
# capture is on and large: g_ext of level 5 into the ground state and 6.1763 MeV
gext="$(awk '$1 == "J" && $6 == "7.5567" { on = 1; next } on && NF == 0 { on = 0 } on && $4 == "L" && ($3 == 2 || $3 == 6) && $6 == "E1" { print $NF == "MeV^(1/2)" ? $(NF - 1) : $NF }' "$WORK/rec/output/parameters.out" | tr '\n' ' ')"
if awk '{ exit !(NF == 2 && $1 !~ /^\(0.000000,0.000000\)$/ && $2 !~ /^\(0.000000,0.000000\)$/) }' <<< "$gext"; then
  ok "level 5's external amplitudes into the ground state and 6.1763 MeV: $gext"
else bad "level 5 has no external amplitude into the ground state or 6.1763 MeV: '$gext'"; fi
for fr in "g+15O 2" "g+15O[6.1763] 6"; do
  set -- $fr
  check "the record's width to $1 is TransformOut's observed width" "$(block_w rec "$1")" "$(level_w rec "$2")" 2e-6
done

edited "$WORK/partial_ec.azr" "$WORK/partial_ec_m.azr" moved || bad "moved twin"
edited "$WORK/record_ec.azr" "$WORK/record_ec_m.azr" moved || bad "moved record"
project pec_m "$WORK/partial_ec_m.azr"
project rec_m "$WORK/record_ec_m.azr"
run pec_m 1
run rec_m 1
same_outputs pec_m rec_m "moved proton width, level energy and ground-state ANC"
check "after the move the record's width to g+15O is still TransformOut's" "$(block_w rec_m g+15O)" "$(level_w rec_m 2 7.5569)" 2e-6

echo "   the gradient against central differences, every free parameter, at three points"
edited "$WORK/record_ec.azr" "$WORK/record_ec_free.azr" free || bad "free variant"
edited "$WORK/record_ec_m.azr" "$WORK/record_ec_free2.azr" free || bad "free variant, second point"
edited "$WORK/record_ec_free2.azr" "$WORK/record_ec_free2.azr" sed -e 's/^\(    total_width  \)[^ ]*/\10.05/' \
  -e 's/^\(    branching    g+15O\[6.1763\] \)[^ ]*/\10.5/' -e 's/^\(    branching    g+15O \)[^ ]*/\10.03/' || bad "second point"
# point 3: also level 6 (8.7502 MeV, 1/2+, Gamma_p = 35.7 keV) as a record; there Im gamma^ext into the ground state
# is 1.3 % of the amplitude, so d gamma_int / d g = |g| / (g^2 - Im^2)^1/2 differs from 1 by 8e-5 (not at level 5)
"$PY" "$M/make_record.py" "$WORK/record_ec_free.azr" "$WORK/record_ec_free3.azr" 6 5.183 fit
edited "$WORK/record_ec_free3.azr" "$WORK/record_ec_free3.azr" sed -e '/^<\/gammaDecays>$/{N;s/^<\/gammaDecays>\n<gammaDecays>$//;}' \
  -e 's/^\(    branching    g+15O\[6.7931\] 0.0\) fit$/\1 fixed/' || bad "the second record"
grep -qx 2 <<< "$(grep -c '^gamma_decay level[56]$' "$WORK/record_ec_free3.azr")" || bad "the second record: not two records in one section"
for pt in "" 2 3; do
  project "gec$pt" "$WORK/record_ec_free$pt.azr"
  (cd "$WORK/gec$pt" && printf '2\n\n\n' | AZURE_GRAD_CHECK=exit run_limited 900 "$AZURE2_BIN" --no-gui --no-readline --use-gradient p.azr > grad.log 2>&1)
  iE="$(awk '/\[grad-check\] idx/ && $4 == "energy" { print $3; exit }' "$WORK/gec$pt/grad.log")"
  (cd "$WORK/gec$pt" && printf '2\n\n\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_H=1e-8 AZURE_GRAD_CHECK_IDX="${iE:-0}" \
     run_limited 900 "$AZURE2_BIN" --no-gui --no-readline --use-gradient p.azr > gradE.log 2>&1)
  rows="$(grep -c '\[grad-check\] idx' "$WORK/gec$pt/grad.log")"
  worst="$(awk '/\[grad-check\] idx/ && $4 != "energy" { for (k = 1; k <= NF; k++) if ($k == "rel" || $k == "relJ") if ($(k + 1) > m) m = $(k + 1) } END { print m + 0 }' "$WORK/gec$pt/grad.log")"
  wE="$(awk '/\[grad-check\] idx/ && $4 == "energy" { for (k = 1; k <= NF; k++) if ($k == "rel" || $k == "relJ") if ($(k + 1) > m) m = $(k + 1) } END { print m + 0 }' "$WORK/gec$pt/gradE.log")"
  # the ANCs: idx 1 (ground state, J1 L1 ch1) and the 6.1763 MeV level (J3 L1 ch1, ch2)
  ancs="$(grep -c '\[grad-check\] idx [0-9]* gamma J[13] L1 ' "$WORK/gec$pt/grad.log")"
  if [ "$rows" -ge 10 ] && [ "$ancs" -eq 3 ] && awk -v w="$worst" 'BEGIN { exit !(w < 1e-6) }'; then
    ok "point ${pt:-1}: $rows parameters (3 ANCs of the final states), gradient and Jacobian paths, largest rel. difference $worst"
  else bad "point ${pt:-1}: $rows parameters, $ancs ANCs, largest rel. difference $worst (limit 1e-6)"; fi
  if awk -v w="$wE" 'BEGIN { exit !(w > 0 && w < 1e-5) }'; then ok "point ${pt:-1}: level energy, rel. difference $wE"
  else bad "point ${pt:-1}: level energy, rel. difference ${wE:-missing} (limit 1e-5)"; fi
done

echo "   fits from displaced starts: the twin's minimum and branchings"
# the twin: level 5's radiative widths x 1.3; the record: Gamma_gamma 0.05, two branchings moved; both: the ground
# state's ANC free from 0.26
anc() { awk '{ n = split($0, a, " "); if (n > 30 && a[9] == "1" && a[27] == "0" && a[7] == "1") { $11 = 0; $12 = 0.26; print " " $0 } else print }' "$1"; }
edited "$WORK/partial_ec.azr" "$WORK/pec_d.azr" awk '{ n = split($0, a, " "); if (n > 30 && a[9] == "5" && a[27] == "10" && a[12] != "0") { $12 = a[12] * 1.3; print " " $0 } else print }' \
  && edited "$WORK/pec_d.azr" "$WORK/pec_d.azr" anc || bad "displaced twin"
edited "$WORK/record_ec.azr" "$WORK/rec_d.azr" sed -e 's/^\(    total_width  \)[^ ]*/\10.05/' -e 's/^\(    branching    g+15O\[6.7931\] \)[^ ]*/\10.3/' \
  -e 's/^\(    branching    g+15O\[6.1763\] \)[^ ]*/\10.5/' && edited "$WORK/rec_d.azr" "$WORK/rec_d.azr" anc || bad "displaced record"
project fit_pec "$WORK/pec_d.azr"
project fit_rec "$WORK/rec_d.azr"
project fit_rec_lm "$WORK/rec_d.azr"
run fit_pec 2 --use-gradient
run fit_rec 2 --use-gradient
run fit_rec_lm 2 --use-lm
c0="$(chi2 fit_pec)"
check "MIGRAD with the gradient: the record's chi^2 minimum equals the twin's" "$(chi2 fit_rec)" "$c0" 1e-5
check "Levenberg-Marquardt: the same minimum" "$(chi2 fit_rec_lm)" "$c0" 1e-5
twin_b fit_pec 7.5567 > "$WORK/twin_b.txt"
check "Gamma_gamma (LM) equals the twin's fitted sum of widths" "$(sav fit_rec_lm level5.total_width)" "$(awk '$1 == "total" { print $2 }' "$WORK/twin_b.txt")" 1e-4
for fr in "g+15O R2" "g+15O[6.7931] R3" "g+15O[6.1763] R6" "g+15O[6.8594] R7"; do
  set -- $fr
  ref="$(awk -v r="$2" '$1 == r { print $2 }' "$WORK/twin_b.txt")"
  check "b($1) (LM) equals the twin's Gamma_f / Gamma_gamma" "$(sav fit_rec_lm "level5.branching.$1")" "$ref" 1e-4
  check "b($1): MIGRAD and LM" "$(sav fit_rec "level5.branching.$1")" "$(sav fit_rec_lm "level5.branching.$1")" 1e-4
done

echo "   the record's fit at 1 and 3 threads"
for t in 1 3; do
  project "ecthr$t" "$WORK/rec_d.azr"
  OMP_NUM_THREADS=$t run "ecthr$t" 2 --use-gradient
done
n=0; d=0
for f in "$WORK/ecthr1"/output/*; do
  n=$((n + 1)); cmp -s "$f" "$WORK/ecthr3/output/$(basename "$f")" || { d=$((d + 1)); echo "        differs: $(basename "$f")"; }
done
if [ "$n" -ge 9 ] && [ "$d" -eq 0 ]; then ok "$n output files identical at 1 and 3 threads"
else bad "$d of $n output files differ between 1 and 3 threads"; fi

echo "   a width below the floor 2P (Im gamma^ext)^2 / n: as TransformIn"
# the ground-state width 3.85e-12 eV, below the floor of about 7e-10 eV
edited "$WORK/partial_ec.azr" "$WORK/pec_floor.azr" awk '{ n = split($0, a, " "); if (n > 30 && a[9] == "5" && a[27] == "10" && a[6] == "2") { $12 = 3.85e-12; print " " $0 } else print }' \
  || bad "floor twin"
"$PY" "$M/make_record.py" "$WORK/pec_floor.azr" "$WORK/rec_floor.azr" 5 5.183 fit
project pec_floor "$WORK/pec_floor.azr"
project rec_floor "$WORK/rec_floor.azr"
run pec_floor 1
run rec_floor 1
same_outputs pec_floor rec_floor "below the floor"
grep -q 'Imaginary portion of the external width is greater' "$WORK/pec_floor/run.log" \
  && grep -q 'below 2P (Im gamma^ext)^2 / n' "$WORK/rec_floor/run.log" \
  && grep -q 'note: the width 3.85e-12 eV to g+15O is below the floor' "$WORK/rec_floor/output/parameters.out" \
  && ok "the twin's TransformIn warning, the record's warning and its note in parameters.out" \
  || bad "the floor: a warning or the note of parameters.out is missing"

echo "-- (i) the vector of observed parameters (note section 2.8): the record's numbers in, the levels' widths out"
# The copies of the compound that place the rate's breakpoints and the meshes' anchors read the levels' observed
# widths: a narrow 14N(p,g) resonance (tests/thick_target_resonance with Gamma_p = 0.01 eV, extrapolations left out,
# a layer yield added) as twin and record; the grid's resonance line (AZR_DEBUG_GRID) gives the photon widths it sums.
TTR="$HERE/../thick_target_resonance"
narrow() { awk '/^<extrapolations>/ { skip = 1; print; next } /^<\/extrapolations>/ { skip = 0 } !skip' "$1" |
  awk '{ n = split($0, a, " "); if (n > 30 && a[9] == "5" && a[27] == "0" && a[12] == "974.937") { $12 = 0.01; print " " $0 } else print }' |
  awk '/^<\/data>/ { print "dataset layer_yield\n    reaction       14N(p,g)15O[total]\n    quantity       yield\n    unit           1/particle\n    file           data/yield_scan.dat\n    columns        E_lab[MeV] _ value error\n    energy_window  0 6 MeV\n    norm           1 fixed\n    target         tgt_targetInt1\n    quadrature     points 20 width_multiplier 5 per_width 50" } { print }'; }
edited "$TTR/thick_target_resonance.azr" "$WORK/nr_partial.azr" narrow || bad "narrow-resonance twin"
"$PY" "$M/make_record.py" "$WORK/nr_partial.azr" "$WORK/nr_record.azr" 5 5.183 fixed
edited "$WORK/nr_partial.azr" "$WORK/nr_partial_fit.azr" sed 's/^\(    norm           360000\) fixed$/\1 fit/' || bad "narrow twin, free norms"
edited "$WORK/nr_record.azr" "$WORK/nr_record_fit.azr" sed 's/^\(    norm           360000\) fixed$/\1 fit/' || bad "narrow record, free norms"
for v in partial record; do
  for k in rate fit; do
    mkdir -p "$WORK/nr_${k}_$v/output" "$WORK/nr_${k}_$v/checks"
    ln -s "$TTR/data" "$WORK/nr_${k}_$v/data"
    cp "$HERE/../reaction_rate/temperatures.dat" "$WORK/nr_${k}_$v/"
  done
  cp "$WORK/nr_$v.azr" "$WORK/nr_rate_$v/p.azr"
  cp "$WORK/nr_${v}_fit.azr" "$WORK/nr_fit_$v/p.azr"
  (cd "$WORK/nr_rate_$v" && printf '5\n\n1\n-1\nyes\ntemperatures.dat\n' | AZR_DEBUG_GRID=1 run_limited 900 "$AZURE2_BIN" --no-gui --no-readline p.azr > run.log 2>&1)
  (cd "$WORK/nr_fit_$v" && printf '2\n\n\n' | AZR_DEBUG_GRID=1 run_limited 900 "$AZURE2_BIN" --no-gui --no-readline p.azr > run.log 2>&1)
done
for k in rate fit; do
  gp="$(grep 'level E=7.556700: input widths' "$WORK/nr_${k}_partial/run.log" | sort -u)"
  gr="$(grep 'level E=7.556700: input widths' "$WORK/nr_${k}_record/run.log" | sort -u)"
  if [ -n "$gp" ] && [ "$gp" = "$gr" ] && [ -s "$WORK/nr_${k}_record/output/chiSquared.out" -o "$k" = rate ]; then
    ok "$k: the grid's copy of the compound has the twin's widths ($(echo "$gp" | sed 's/.*\(gamma=[^ ]*\).*/\1/'))"
  else bad "$k: the grid's widths differ: twin '$(echo "$gp" | sed 's/.*\(gamma=[^ ]*\).*/\1/')', record '$(echo "$gr" | sed 's/.*\(gamma=[^ ]*\).*/\1/')'"; fi
done
cmp -s "$WORK/nr_rate_partial/output/reactionrates.out" "$WORK/nr_rate_record/output/reactionrates.out" \
  && ok "the narrow resonance's rate (total capture) is the twin's" || bad "the narrow resonance's rate differs from the twin's"

# A limit on a particle width of the record's level is converted as the twin's (both end at the 950 eV limit); a limit
# on a record's number is applied as written (Gamma_gamma ends at 0.035 eV).
gp_of() { awk '$1 == "J" && $6 == "7.5567" { on = 1; next } on && NF == 0 { on = 0 } on && $3 == 1 && $6 == "0" { print $12 $13; exit }' "$WORK/$1/output/parameters.out"; }
free_gp() { awk '{ n = split($0, a, " "); if (n > 30 && a[9] == "5" && a[27] == "0" && a[7] == "1") { $11 = 0; print " " $0 } else print }' "$1"; }
for v in partial record; do
  edited "$M/$v.azr" "$WORK/lim_$v.azr" free_gp || bad "free proton width ($v)"
  printf '<parameterSettings>\nLevel 5 Channel 1 Width (eV) 974.937 900 950 0 0 0 level 0\n</parameterSettings>\n' >> "$WORK/lim_$v.azr"
  project "lim_$v" "$WORK/lim_$v.azr"
  run "lim_$v" 2
done
if [ "$(gp_of lim_partial)" = "950.000000eV" ] && [ "$(gp_of lim_record)" = "950.000000eV" ] && [ "$(chi2 lim_partial)" = "$(chi2 lim_record)" ]; then
  ok "a limit on level 5's proton width: both end at $(gp_of lim_record), chi^2 $(chi2 lim_record)"
else bad "a limit on level 5's proton width: twin $(gp_of lim_partial), record $(gp_of lim_record)"; fi
edited "$WORK/lim_record.azr" "$WORK/limq.azr" sed 's/^Level 5 Channel 1 Width (eV) 974.937 900 950 0 0 0 level 0$/level5.total_width 0.0385027 0.03 0.035 0 0 0 radiative 0/' \
  || bad "limit on Gamma_gamma"
project limq "$WORK/limq.azr"
run limq 2
check "a limit (0.03, 0.035) eV on Gamma_gamma: the fit ends at it" "$(sav limq level5.total_width)" 0.035 1e-6

# pyazr (tests/pyazr/gamma_decay_test.py part 5) and the MCMC set the vector in observed form: the record's numbers of
# the vector must reach the engine.  Physical-parameter MCMC with the record's numbers free: the likelihood moves with
# them, param.mcmc carries the best walker's numbers, and its -2 ln L is the chi^2 of a calculation at them.
project mcmc "$M/record.azr"
(cd "$WORK/mcmc" && printf '6\n\n\n4\n100\n0.5\n\n1\nno\nyes\n' | OMP_NUM_THREADS=1 run_limited 900 "$AZURE2_BIN" --no-gui --no-readline p.azr > run.log 2>&1)
if grep -q "Parameter mode: Physical" "$WORK/mcmc/run.log" && [ -s "$WORK/mcmc/output/samples.mcmc" ]; then
  spread="$(awk -F, 'NR > 1 { if (a == "" || $4 < a) a = $4; if (b == "" || $4 > b) b = $4 } END { printf "%.6g", b - a }' "$WORK/mcmc/output/samples.mcmc")"
  awk -v s="$spread" 'BEGIN { exit !(s > 0.1) }' && ok "MCMC (observed parameters): ln L moves with the record's numbers (spread $spread)" \
    || bad "MCMC (observed parameters): ln L does not move with the record's numbers (spread $spread)"
  read -r logL g0 < <(awk -F, 'NR > 1 { s = $1 } END { print s }' "$WORK/mcmc/output/samples.mcmc" | {
    read -r last; awk -F, -v last="$last" 'NR > 1 && $1 == last && ($4 > best || best == "") { best = $4; g = $6 } END { print best, g }' "$WORK/mcmc/output/samples.mcmc"; })
  check "param.mcmc carries the best walker's Gamma_gamma" "$(awk '$1 == "level5.total_width" { print $2 }' "$WORK/mcmc/output/param.mcmc")" "$g0" 1e-6
  mkdir -p "$WORK/mcmc_calc/output" "$WORK/mcmc_calc/checks"; ln -s "$M/data" "$WORK/mcmc_calc/data"
  cp "$M/record.azr" "$WORK/mcmc_calc/p.azr"; cp "$WORK/mcmc/output/param.mcmc" "$WORK/mcmc_calc/best.par"
  (cd "$WORK/mcmc_calc" && printf '1\nbest.par\n\n' | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline p.azr > run.log 2>&1)
  check "-2 ln L of that walker is the chi^2 of a calculation at its parameters" "$(chi2 mcmc_calc)" \
    "$(awk -v l="$logL" 'BEGIN { printf "%.10g", -2 * l }')" 1e-4
else bad "the observed-parameter MCMC did not run"; fi

echo "-- (j) a final state's energy free: the pull-back's term in E_f"
# The external-capture twin with a second particle pair (alpha + 11C, its threshold put at 5.0 MeV, in a 9/2+ group of
# its own, width 0): the 6.1763 MeV final state lies above it, so its energy is a free parameter (CNuc::FillMnParams).
alpha() { awk '/^<\/levels>/ { print " 4.5 1 20 1 1 9 3 6 34 1 1 0 0 1 1.5 -1 0 4.0026 11.0114 2 6 7.297 5.0 0 0 0.0 0 5.5 0 0 0" } { print }' "$1"; }
efree() { awk '{ n = split($0, a, " "); if (n > 30 && a[9] == "9") { $4 = 0; if (a[27] == "0" && a[12] != "0") $11 = 0; print " " $0 }
  else if (n > 30 && a[9] == "5" && a[27] == "0" && a[7] == "1") { $11 = 0; print " " $0 } else print }' "$1"; }
edited "$WORK/partial_ec.azr" "$WORK/pec_alpha.azr" alpha || bad "the second pair"
"$PY" "$M/make_record.py" "$WORK/pec_alpha.azr" "$WORK/rec_alpha.azr" 5 5.183 fit
edited "$WORK/rec_alpha.azr" "$WORK/rec_ef.azr" efree || bad "the free final-state energy"
project ef "$WORK/rec_ef.azr"
(cd "$WORK/ef" && printf '2\n\n\n' | AZURE_GRAD_CHECK=exit run_limited 900 "$AZURE2_BIN" --no-gui --no-readline --use-gradient p.azr > grad.log 2>&1)
rows="$(grep -c '\[grad-check\] idx' "$WORK/ef/grad.log")"
ef_line="$(grep '\[grad-check\] idx [0-9]* energy J3 L1 ' "$WORK/ef/grad.log")"
wEf="$(echo "$ef_line" | awk '{ for (k = 1; k <= NF; k++) if ($k == "rel" || $k == "relJ") if ($(k + 1) > m) m = $(k + 1) } END { print m + 0 }')"
worst="$(awk '/\[grad-check\] idx/ { for (k = 1; k <= NF; k++) if ($k == "rel" || $k == "relJ") if ($(k + 1) > m) m = $(k + 1) } END { print m + 0 }' "$WORK/ef/grad.log")"
if [ -n "$ef_line" ] && [ "$rows" -ge 8 ] && awk -v w="$worst" 'BEGIN { exit !(w < 1e-6) }'; then
  ok "$rows parameters, the 6.1763 MeV state's energy among them (rel. $wEf), largest rel. difference $worst"
else bad "$rows parameters, the final state's energy row '${ef_line:-missing}', largest rel. difference $worst (limit 1e-6)"; fi

echo
if [ "$failures" -eq 0 ]; then echo "gamma_decay: all checks passed"; else echo "gamma_decay: $failures check(s) FAILED"; fi
exit $((failures > 0))
