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
#
# The fitted sample thickness ("transmission n vary sigma", "yield n vary sigma", "thicknessgroup <label>"):
#
#   (h) data = the model's T (and Y0) at n = 0.18 atoms/b, 1 % errors, every level fixed: from "transmission 0.15 1 0"
#       (and "yield 0.15 1 0") MIGRAD (--use-gradient) and --use-lm recover segment_1_thickness = 0.18 to 1e-4;
#   (i) a thicknessgroup of the transmission and the yield: one parameter from both, the follower's fixed
#       segment_2_thickness written as the leader's value;
#   (j) AZURE_GRAD_CHECK with seven thicknesses free with every width: lines without a target effect, under a
#       Gaussian spread (with a prior), under doppler alone (a thickness group) and under doppler + the Gaussian:
#       analytic = finite differences to 1e-6 of the largest, the Jacobian path to 1e-10;
#   (k) a prior: "transmission 0.15 1 0.01" on 20 % data holds n between 0.15 and 0.18, and chiSquared.out's
#       Total-Thickness-Chi-Squared = ((n - 0.15)/0.01)^2;
#   (l) refused with an ERROR naming the segment: vary 2, sigma < 0, sigma missing, n = 0, thicknessgroup without
#       a label, a thicknessgroup member without its own transmission / yield block;
#   (m) with REF_AZURE2 set (a build from before the fit form, e.g. 7377779): the plain "transmission n" / "yield n"
#       outputs are byte-identical to it.
#
#   ./tests/transmission_yield/check.sh path/to/AZURE2
#   REF_AZURE2=path/to/reference/AZURE2 ./tests/transmission_yield/check.sh path/to/AZURE2

set -uo pipefail
# Byte-identical files without cmp (not in every MSYS2 image: the Windows runner has none); \r stripped so a
# CRLF-writing Windows build compares like the others.
same_file() { [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
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

# ---------------------------------------------------------------------------------------------------------------
# The fitted sample thickness
fit() {  # fit <name> <model> <flags>: mode 2, no parameter file
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$WORK/models/data" "$dir/data"
  cp "$WORK/models/$2.azr" "$dir/t.azr"
  (cd "$dir" && printf '2\nn\n\n\n' | run_limited 900 "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
}
par() { awk -v n="$2" '$1 == n { print $2; exit }' "$WORK/$1/output/param.sav" 2>/dev/null; }
chi() { sed -n "s/.*$2: \([^ ]*\).*/\1/p" "$WORK/$1/output/chiSquared.out" 2>/dev/null | head -1; }
near() { awk -v a="$1" -v b="$2" -v t="$3" -v re="$num" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(a ~ re && d <= t) }'; }

# (h), (i) fits on data made at n = 0.18
run thin18 thin18 1 n12C
if "$PY" -I "$HERE/make_models.py" fitdata "$HERE/../reaction_rate/models/n12C.azr" "$WORK/thin18" "$WORK/models"; then
  for v in "fitT:transmission" "fitY:yield"; do
    m="${v%%:*}"; what="${v#*:}"
    for flags in --use-gradient --use-lm; do
      fit "${m}$flags" "$m" "$flags"
      n="$(par "${m}$flags" segment_1_thickness)"
      if near "$n" 0.18 1e-4; then ok "(h) $what, $flags: segment_1_thickness = $n from 0.15 (data at 0.18)"
      else bad "(h) $what, $flags: segment_1_thickness '${n}' (data at 0.18)"; tail -3 "$WORK/${m}$flags/run.log"; fi
    done
  done
  for flags in --use-gradient --use-lm; do
    fit "fitG$flags" fitG "$flags"
    n1="$(par "fitG$flags" segment_1_thickness)"; n2="$(par "fitG$flags" segment_2_thickness)"
    if near "$n1" 0.18 1e-4 && [ "$n1" = "$n2" ] && grep -q "thickness group 's': segments #1 + 2 share" "$WORK/fitG$flags/run.log"; then
      ok "(i) thicknessgroup, $flags: one thickness $n1 from the transmission and the yield (follower's parameter $n2)"
    else bad "(i) thicknessgroup, $flags: segment_1_thickness '${n1}', segment_2_thickness '${n2}'"; tail -3 "$WORK/fitG$flags/run.log"; fi
  done
  # (k) a prior
  fit prior fitP --use-gradient
  n="$(par prior segment_1_thickness)"; tc="$(chi prior Total-Thickness-Chi-Squared)"
  if [ -n "$n" ] && [ -n "$tc" ] && awk -v n="$n" -v s="$tc" -v re="$num" 'BEGIN { p = ((n - 0.15) / 0.01)^2; d = s - p
         if (d < 0) d = -d; exit !(n ~ re && s ~ re && n > 0.151 && n < 0.179 && d <= 1e-5 * (p + 1)) }'; then
    ok "(k) a 0.01 atoms/b prior holds n at $n (data alone: 0.18); Total-Thickness-Chi-Squared $tc = ((n - 0.15)/0.01)^2"
  else bad "(k) prior fit: segment_1_thickness '${n}', Total-Thickness-Chi-Squared '${tc}'"; tail -3 "$WORK/prior/run.log"; fi
else bad "(h) fit data not built"; fi

# (j) gradient check with the thicknesses free
mkdir -p "$WORK/tfit/output" "$WORK/tfit/checks"; ln -s "$WORK/models/data" "$WORK/tfit/data"
cp "$WORK/models/grad_fit.azr" "$WORK/tfit/t.azr"
(cd "$WORK/tfit" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
   run_limited 1200 "$AZURE2_BIN" --use-gradient --no-grid-replacement --no-gui --no-readline t.azr > run.log 2>&1)
awk '{ sub(/\r$/, "") }
     /^\[grad-check\] idx / { a = f = j = ""
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print a, f, j, $4 }' "$WORK/tfit/run.log" > "$WORK/tfit/rows"
nrow="$(awk 'END { print NR }' "$WORK/tfit/rows")"
nthk="$(awk '$4 == "thickness"' "$WORK/tfit/rows" | grep -c .)"
if [ "$nthk" -ne 7 ] || [ "$nrow" -lt 10 ]; then bad "(j) gradient check listed $nrow parameters, $nthk thicknesses (want 7)"
else
  read -r fd jac <<< "$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a
                               d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d
                               e = $1 - $3; if (e < 0) e = -e; if (e > v) v = e }
                             END { printf "%.2e %.2e", w / m, v / m }' "$WORK/tfit/rows")"
  read -r tfd <<< "$(awk '$4 == "thickness" { d = ($1 - $2) / $2; if (d < 0) d = -d; if (d > w) w = d } END { printf "%.2e", w }' \
                     "$WORK/tfit/rows")"
  le "$fd" 1e-6 && ok "(j) gradient = finite differences for $nrow parameters, $nthk thicknesses among them ($fd of the largest; the thicknesses to $tfd of their own)" \
    || bad "(j) gradient vs finite differences: $fd of the largest (thicknesses $tfd of their own)"
  le "$jac" 1e-10 && ok "(j) Jacobian-path gradient = adjoint gradient ($jac of the largest)" \
    || bad "(j) Jacobian-path vs adjoint gradient: $jac of the largest"
fi

# (l) malformed fit forms and group members
refuse_t() {
  local name="$1" want="$2"
  run "$name" "$name" 1 n12C
  if grep -q "ERROR: [A-Za-z ]*segment #[0-9]*: .*$want" "$WORK/$name/run.log" && [ ! -s "$WORK/$name/output/chiSquared.out" ]; then
    ok "(l) refused: $name ($(grep -m1 -o 'ERROR: .*' "$WORK/$name/run.log" | cut -c1-110))"
  else bad "(l) not refused: $name"; grep -m2 -E 'ERROR|WARN' "$WORK/$name/run.log"; fi
}
refuse_t refuse_t_vary "vary 0 or 1"
refuse_t refuse_t_sigma "sigma must be >= 0"
refuse_t refuse_t_nosigma "needs the prior width sigma"
refuse_t refuse_t_zero "needs the sample thickness n > 0"
refuse_t refuse_t_nolabel "thicknessgroup needs a label"
refuse_t refuse_t_member "thicknessgroup 's' needs the segment's own transmission or yield block"

# (m) plain "transmission n" / "yield n" byte-identical to a reference build (optional)
if [ -n "${REF_AZURE2:-}" ]; then
  ref="$AZURE2_BIN"; AZURE2_BIN="$REF_AZURE2"
  for m in thin conv; do run "ref_$m" "$m" 1 n12C; done
  run ref_xthin xthin 3 n12C
  run ref_natT natT 1 neon
  AZURE2_BIN="$ref"
  same=0; diff=0
  for m in thin conv xthin natT; do
    for f in "$WORK/ref_$m"/output/AZUREOut_*.out "$WORK/ref_$m"/output/AZUREOut_*.extrap "$WORK/ref_$m/output/chiSquared.out"; do
      [ -e "$f" ] || continue
      if same_file "$f" "$WORK/$m/output/$(basename "$f")"; then same=$((same + 1)); else diff=$((diff + 1)); fi
    done
  done
  if [ "$diff" -eq 0 ] && [ "$same" -ge 8 ]; then ok "(m) plain transmission / yield: $same output files byte-identical to $REF_AZURE2"
  else bad "(m) plain transmission / yield: $diff of $((same + diff)) output files differ from $REF_AZURE2"; fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: transmission and self-shielded yield"; else echo "FAIL: $failures check(s)"; exit 1; fi
