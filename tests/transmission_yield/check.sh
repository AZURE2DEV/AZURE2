#!/usr/bin/env bash
#
# "transmission <n>" and "yield <n>" on a segment line: the transmission T = exp(-n sigma_t) of a sample n atoms/barn
# thick, and SAMMY's self-shielded capture or reaction yield Y0 = (sigma / sigma_t)(1 - exp(-n sigma_t)), sigma_t the
# total cross section of the (neutral) entrance pair by the optical theorem.  Under a sub-point target effect the
# transform is applied at every sub-point BEFORE the kernel (SAMMY manual III.E.1), never to the broadened cross
# section.  Models (make_models.py) from tests/reaction_rate/models/n12C.azr (n + 12C, s-wave 1/2+ resonance at
# E_cm = 0.425 MeV, Gamma_n = 34 keV; n = 0.15 atoms/b gives n sigma_t ~ 1 at the peak) and tests/natural_neon:
#
#   (a) no spread: T = exp(-n sigma_t) and Y0 = (sg/st)(1 - exp(-n st)) from the plain segments, to 1e-9;
#   (b) Gaussian spread (15 keV lab): T and Y0 equal the kernel average of the transformed dense, unbroadened
#       cross sections (an independent integration of a 0.1 keV grid), to 2e-4; transforming the BROADENED cross
#       sections instead is off by more than 1e-2 at the dip (the test discriminates);
#   (c) a natural (three-nucleus) sum: T = exp(-n sigma_nat), to 1e-9;
#   (d) the analytic gradient equals finite differences and the Jacobian path (AZURE_GRAD_CHECK) on a project with
#       transmission and yield segments with and without the spread, every width free.  --no-grid-replacement: a
#       finite-difference step of the 1/2+ energy or width otherwise re-places the sub-point grids around it, and the
#       difference quotient then differs from the derivative at fixed grids by 1e-2 on those two parameters (as for
#       any sub-point segment, transformed or not);
#   (e) refused, with an ERROR and no output: transmission and yield on one line, yield on the elastic total,
#       transmission on a differential segment, a charged entrance pair, a natural sum under a sub-point spread
#       (its nuclei's terms are on different grids), and yield on a sample segment;
#   (f) test segments (mode 3): the transmission and yield curves are the formula of the plain curves, to 1e-9;
#   (g) a yield's sigma_t is the total, not the elastic: n + 20Ne with 1 keV widths in the reaction channels, a yield
#       to alpha + 17O equals the formula with "total 1" to 1e-9 where the elastic is > 1e-3 off the total.
#
#   ./tests/transmission_yield/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PY="$(command -v python3 || command -v python)"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/transmission_yield.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
# a must be a number: an empty string (analyze.py crashed) or nan/inf is a failure, never 0
num='^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$'
le() { awk -v a="$1" -v b="$2" -v re="$num" 'BEGIN { exit !(a ~ re && a + 0 <= b + 0) }'; }
gt() { awk -v a="$1" -v b="$2" -v re="$num" 'BEGIN { exit !(a ~ re && a + 0 > b + 0) }'; }

N=0.15; NNAT=0.05; SIG=0.015; M1=1.008665; M2=12
"$PY" -I "$HERE/make_models.py" "$HERE/../reaction_rate/models/n12C.azr" "$HERE/../natural_neon/natural_neon.azr" "$WORK/models" ||
  { bad "models not built"; exit 1; }

# run NAME MODEL MODE DATA -> $WORK/NAME (DATA: n12C or neon)
run() {
  local dir="$WORK/$1" input
  mkdir -p "$dir/output" "$dir/checks"
  if [ "$4" = neon ]; then ln -s "$HERE/../natural_neon/data" "$dir/data"; else ln -s "$WORK/models/data" "$dir/data"; fi
  cp "$WORK/models/$2.azr" "$dir/t.azr"
  if [ "$3" = 3 ]; then input='3\n\n\n7\n'; else input='1\n\n\n7\n'; fi
  (cd "$dir" && printf "$input" | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}

for m in bare thin conv bconv; do run $m $m 1 n12C; done
run truth truth 3 n12C

# (a) no spread
for k in T Y; do
  read -r w rest <<< "$("$PY" -I "$HERE/analyze.py" formula $k "$WORK/bare" "$WORK/thin" $N)"
  if le "$w" 1e-9; then ok "(a) $k without a spread = formula of the plain cross sections (max rel $w, $rest)"
  else bad "(a) $k without a spread vs formula: max rel $w ($rest)"; fi
done

# (b) Gaussian spread: transform first, then the kernel
for k in T Y; do
  read -r w rest <<< "$("$PY" -I "$HERE/analyze.py" conv $k "$WORK/truth" "$WORK/conv" $N $SIG $M1 $M2)"
  if le "$w" 2e-4; then ok "(b) $k with a ${SIG} MeV spread = kernel average of the transformed cross section (max rel $w, $rest)"
  else bad "(b) $k with a spread vs independent integration: max rel $w ($rest)"; fi
  read -r w rest <<< "$("$PY" -I "$HERE/analyze.py" wrong $k "$WORK/truth" "$WORK/bconv" $N $SIG $M1 $M2)"
  if gt "$w" 1e-2; then ok "(b) $k of the broadened cross section is off by $w: the order matters here ($rest)"
  else bad "(b) $k of the broadened cross section is only $w off: the test does not discriminate"; fi
done

# (c) natural sum of three nuclei, no target effect
run nat nat 1 neon
run natT natT 1 neon
read -r w rest <<< "$("$PY" -I "$HERE/analyze.py" natural "$WORK/nat" "$WORK/natT" $NNAT)"
if le "$w" 1e-9; then ok "(c) natural sum: T = exp(-n sigma_nat) (max rel $w, $rest)"
else bad "(c) natural sum: T vs exp(-n sigma_nat): max rel $w ($rest)"; fi
grep -q "transmission of n = $NNAT atoms/b" "$WORK/natT/output/AZUREOut_sample_seg=1.out" 2>/dev/null &&
  ok "(c) the sample decomposition says the total column is the transmission" ||
  bad "(c) no transmission note in AZUREOut_sample_seg=1.out"

# (d) gradient check
run gcalc grad 1 n12C
mkdir -p "$WORK/fit/output" "$WORK/fit/checks"; ln -s "$WORK/models/data" "$WORK/fit/data"; cp "$WORK/models/grad.azr" "$WORK/fit/t.azr"
(cd "$WORK/fit" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
   run_limited 1200 "$AZURE2_BIN" --use-gradient --no-grid-replacement --no-gui --no-readline t.azr > run.log 2>&1)
chi="$(awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { sub(",", "", $3); print $3 }' "$WORK/fit/run.log")"
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/gcalc/output/chiSquared.out")"
if [ -n "$chi" ] && awk -v a="$chi" -v b="$want" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(d <= 1e-5) }'; then
  ok "(d) chi2 at the start of the fit $chi = calculation $want"
else bad "(d) chi2 at the start of the fit '$chi', calculation $want"; fi
awk '{ sub(/\r$/, "") }
     /^\[grad-check\] idx / { a = f = j = ""
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print a, f, j }' "$WORK/fit/run.log" > "$WORK/fit/rows"
nrow="$(awk 'END { print NR }' "$WORK/fit/rows")"
if [ "$nrow" -lt 3 ]; then bad "(d) gradient check listed $nrow parameters"
else
  read -r fd jac <<< "$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a
                               d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d
                               e = $1 - $3; if (e < 0) e = -e; if (e > v) v = e }
                             END { printf "%.2e %.2e", w / m, v / m }' "$WORK/fit/rows")"
  le "$fd" 1e-7 && ok "(d) gradient = finite differences for $nrow parameters ($fd of the largest)" \
    || bad "(d) gradient vs finite differences: $fd of the largest"
  le "$jac" 1e-10 && ok "(d) Jacobian-path gradient = adjoint gradient ($jac of the largest)" \
    || bad "(d) Jacobian-path vs adjoint gradient: $jac of the largest"
fi

# (f) test segments (mode 3)
run xbare xbare 3 n12C
run xthin xthin 3 n12C
for k in T Y; do
  read -r w rest <<< "$("$PY" -I "$HERE/analyze.py" xformula $k "$WORK/xbare" "$WORK/xthin" $N)"
  if le "$w" 1e-9; then ok "(f) test segment: $k = formula of the plain test curves (max rel $w, $rest)"
  else bad "(f) test segment: $k vs formula: max rel $w ($rest)"; fi
done

# (g) the yield's sigma_t is the total cross section (elastic + reactions), not the elastic
run ybare ybare 1 neon
run yelast yelast 1 neon
run yreac yreac 1 neon
read -r w npts _ sep <<< "$("$PY" -I "$HERE/analyze.py" reaction "$WORK/ybare" "$WORK/yelast" "$WORK/yreac" $NNAT)"
if le "$w" 1e-9 && gt "$sep" 1e-3; then
  ok "(g) Y0 of n+20Ne -> alpha+17O uses the total (max rel $w, $npts open points; elastic is $sep off the total)"
else bad "(g) Y0 to alpha+17O: max rel $w vs the formula with the total ($npts points; elastic vs total $sep)"; fi

# (e) refusals
refuse() {
  local name="$1" data="$2" want="$3"
  run "$name" "$name" 1 "$data"
  if grep -q "ERROR: segment #[0-9]* ($want" "$WORK/$name/run.log" && [ ! -s "$WORK/$name/output/chiSquared.out" ]; then
    ok "(e) refused: $name"
  else bad "(e) not refused: $name"; grep -m2 -E 'ERROR|WARN' "$WORK/$name/run.log"; fi
}
refuse refuse_both n12C transmission
refuse refuse_yield_total n12C yield
refuse refuse_diff n12C transmission
refuse refuse_charged n12C transmission
refuse refuse_natsub neon transmission
refuse refuse_natyield neon yield

if [ "$failures" -eq 0 ]; then echo "PASS: transmission and self-shielded yield"; else echo "FAIL: $failures check(s)"; exit 1; fi
