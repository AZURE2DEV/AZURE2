#!/usr/bin/env bash
#
# Free-gas Doppler broadening ("doppler T" on a <targetInt> line, T the effective temperature in kelvin): SAMMY's
# free-gas model (manual, "Free-Gas Model of Doppler Broadening"), applied to the cross section of every sub-point
# BEFORE the transmission/yield transform and the line's resolution kernel.  Model (make_models.py): n + 12C from
# tests/reaction_rate/models/n12C.azr with the 1/2+ resonance (E_cm = 0.4247 MeV) narrowed to Gamma_n = 200 eV, so
# that 300 K broadening (Doppler width ~60 eV) is a large effect.  Independent integrals in analyze.py, on a 1 eV
# grid of unbroadened cross sections (mode 3):
#
#   (a) doppler alone: sigma_t and sigma_gamma at the data points = the free-gas integral, to 2e-4; the step
#       changes sigma_t by more than 5 % somewhere (the check discriminates);
#   (b) doppler + Gaussian resolution: the Gaussian average of sigma_D, to 2e-4; and the Gaussian alone differs
#       from it by more than 1 %;
#   (c) transmission and yield: with doppler alone, T = exp(-N sigma_D) and Y0 of sigma_D at the point; with the
#       Gaussian too, the Gaussian average of the transform of sigma_D (Doppler, transform, kernel), to 2e-4; the
#       other order (the transform of the broadened sigma) is more than 1 % off;
#   (d) the analytic gradient equals finite differences and the Jacobian path (AZURE_GRAD_CHECK) on a project with
#       all three kinds, every free parameter varied, --no-grid-replacement (as tests/transmission_yield);
#   (e) refused with an ERROR and no output: "doppler 0" and "doppler" without a temperature.
#
# The fitted temperature, "doppler T vary sigma [Tmax]" (parameter doppler_<n>):
#
#   (f) data = the dop model's sigma_t and sigma_gamma at T = 450 K (1 % errors), every level fixed, the line
#       "doppler 300 1 0 900": MIGRAD (--use-gradient) and --use-lm recover T within 1 K;
#   (g) AZURE_GRAD_CHECK on grad_fit (doppler_1 "300 1 0 900" on the Doppler-only and transmission/yield
#       segments, doppler_2 "300 1 20" on the Doppler + Gaussian transmission/yield ones), free with the widths:
#       the analytic gradient = finite differences (1e-6 of the largest) and the Jacobian path (1e-10);
#   (h) "doppler 300 1 20" against the same data with 20 % errors: T lands between 300 K and the default layout
#       temperature 363 K, and Total-Doppler-Chi-Squared in chiSquared.out = ((T - 300)/20)^2;
#   (i) malformed fit forms refused with an ERROR and no output: vary 2, sigma missing, sigma < 0, Tmax <= T,
#       T < 0;
#   (j) only with REF_AZURE2=<an f373664 build>: the plain "doppler T" models' AZUREOut files and chiSquared.out
#       are byte-identical to that build's (the fit form changed nothing for a fixed line).
#
#   [REF_AZURE2=path/to/reference/AZURE2] ./tests/doppler_fgm/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PY="$(command -v python3 || command -v python)"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/doppler_fgm.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
# a must be a number: an empty string (analyze.py crashed) or nan/inf is a failure, never 0
num='^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$'
le() { awk -v a="$1" -v b="$2" -v re="$num" 'BEGIN { exit !(a ~ re && a + 0 <= b + 0) }'; }
gt() { awk -v a="$1" -v b="$2" -v re="$num" 'BEGIN { exit !(a ~ re && a + 0 > b + 0) }'; }

N=0.15; T_K=300; SIG=0.0001; M1=1.008665; M2=12; TOL=2e-4
"$PY" -I "$HERE/make_models.py" "$HERE/../reaction_rate/models/n12C.azr" "$WORK/models" || { bad "models not built"; exit 1; }

run() {
  local dir="$WORK/$1" input
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$WORK/models/data" "$dir/data"
  cp "$WORK/models/$2.azr" "$dir/t.azr"
  if [ "$3" = 3 ]; then input='3\n\n\n7\n'; else input='1\n\n\n7\n'; fi
  (cd "$dir" && printf "$input" | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}
an() { "$PY" -I "$HERE/analyze.py" "$@"; }

run truth truth 3
for m in bare dop dopg gauss dopT dopgT; do run $m $m 1; done

# (a)
for k in T G; do
  read -r w rest <<< "$(an doppler $k "$WORK/truth" "$WORK/dop" $T_K $M1 $M2)"
  if le "$w" $TOL; then ok "(a) doppler alone, sigma_$k = free-gas integral (max rel $w, $rest)"
  else bad "(a) doppler alone, sigma_$k vs free-gas integral: max rel $w ($rest)"; fi
done
read -r w rest <<< "$(an effect "$WORK/bare" "$WORK/dop")"
if gt "$w" 0.05; then ok "(a) the Doppler step changes sigma_t by up to $w ($rest)"
else bad "(a) the Doppler step changes sigma_t by only $w: the check does not discriminate"; fi

# (b)
for k in T G; do
  read -r w rest <<< "$(an dopgauss $k "$WORK/truth" "$WORK/dopg" $T_K $M1 $M2 $SIG)"
  if le "$w" $TOL; then ok "(b) doppler + Gaussian, sigma_$k = Gaussian average of sigma_D (max rel $w, $rest)"
  else bad "(b) doppler + Gaussian, sigma_$k: max rel $w ($rest)"; fi
done
read -r w rest <<< "$(an effect "$WORK/gauss" "$WORK/dopg")"
if gt "$w" 0.01; then ok "(b) under the Gaussian the Doppler step still changes sigma_t by up to $w"
else bad "(b) under the Gaussian the Doppler step changes sigma_t by only $w"; fi

# (c)
for k in T Y; do
  read -r w rest <<< "$(an dopT $k "$WORK/truth" "$WORK/dopT" $T_K $M1 $M2 $N)"
  if le "$w" $TOL; then ok "(c) doppler alone, $k = transform of sigma_D (max rel $w, $rest)"
  else bad "(c) doppler alone, $k: max rel $w ($rest)"; fi
  read -r w rest <<< "$(an dopT $k "$WORK/truth" "$WORK/dopgT" $T_K $M1 $M2 $N $SIG)"
  if le "$w" $TOL; then ok "(c) doppler + Gaussian, $k = Gaussian average of the transform of sigma_D (max rel $w, $rest)"
  else bad "(c) doppler + Gaussian, $k: max rel $w ($rest)"; fi
done
read -r w rest <<< "$(an order "$WORK/truth" "$WORK/dopgT" $T_K $M1 $M2 $N $SIG)"
if gt "$w" 0.01; then ok "(c) transforming the broadened sigma instead is off by up to $w: the order matters here"
else bad "(c) the other order is only $w off: the check does not discriminate"; fi

# (d) gradient check
run gcalc grad 1
mkdir -p "$WORK/fit/output" "$WORK/fit/checks"; ln -s "$WORK/models/data" "$WORK/fit/data"; cp "$WORK/models/grad.azr" "$WORK/fit/t.azr"
(cd "$WORK/fit" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
   run_limited 1200 "$AZURE2_BIN" --use-gradient --no-grid-replacement --no-gui --no-readline t.azr > run.log 2>&1)
chi="$(awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { sub(",", "", $3); print $3 }' "$WORK/fit/run.log")"
want="$(awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/gcalc/output/chiSquared.out")"
if [ -n "$chi" ] && [ -n "$want" ] && awk -v a="$chi" -v b="$want" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(d <= 1e-5) }'; then
  ok "(d) chi2 at the start of the fit $chi = calculation $want"
else bad "(d) chi2 at the start of the fit '$chi', calculation '$want'"; fi
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
  le "$fd" 1e-6 && ok "(d) gradient = finite differences for $nrow parameters ($fd of the largest)" \
    || bad "(d) gradient vs finite differences: $fd of the largest"
  le "$jac" 1e-10 && ok "(d) Jacobian-path gradient = adjoint gradient ($jac of the largest)" \
    || bad "(d) Jacobian-path vs adjoint gradient: $jac of the largest"
fi

# (e) refusals
for name in refuse_zero refuse_missing; do
  run "$name" "$name" 1
  if grep -q "a doppler block is 'doppler T'" "$WORK/$name/run.log" && [ ! -s "$WORK/$name/output/chiSquared.out" ]; then
    ok "(e) refused: $name"
  else bad "(e) not refused: $name"; grep -m2 -E 'ERROR|WARN' "$WORK/$name/run.log"; fi
done

# (f) a fitted temperature: data at T = 450 K, the fit starts at 300 K
fit() {  # fit <name> <model> <flags>: mode 2, no parameter file
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$WORK/models/data" "$dir/data"
  cp "$WORK/models/$2.azr" "$dir/t.azr"
  (cd "$dir" && printf '2\nn\n\n\n' | run_limited 900 "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
}
par() { awk -v n="$2" '$1 == n { print $2; exit }' "$WORK/$1/output/param.sav"; }
chi() { sed -n "s/.*$2: \([^ ]*\).*/\1/p" "$WORK/$1/output/chiSquared.out" | head -1; }
run dop450 dop450 1
if "$PY" -I "$HERE/make_models.py" fitdata "$HERE/../reaction_rate/models/n12C.azr" "$WORK/dop450" "$WORK/models"; then
  for v in "fit_grad:--use-gradient" "fit_lm:--use-lm"; do
    name="${v%%:*}"; flags="${v#*:}"
    fit "$name" fit "$flags"
    t="$(par "$name" doppler_1)"; c="$(chi "$name" Total-Chi-Squared)"
    if [ -n "$t" ] && awk -v t="$t" -v re="$num" 'BEGIN { d = t - 450; if (d < 0) d = -d; exit !(t ~ re && d < 1) }'; then
      ok "(f) $flags: doppler_1 = $t K from 300 K (truth 450 K), chi2 $c"
    else bad "(f) $flags: doppler_1 '${t}', chi2 '${c}'"; tail -5 "$WORK/$name/run.log"; fi
  done
else bad "(f) fit data not built"; fi

# (g) gradient check with the temperatures free
mkdir -p "$WORK/gfit/output" "$WORK/gfit/checks"; ln -s "$WORK/models/data" "$WORK/gfit/data"
cp "$WORK/models/grad_fit.azr" "$WORK/gfit/t.azr"
(cd "$WORK/gfit" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
   run_limited 1200 "$AZURE2_BIN" --use-gradient --no-grid-replacement --no-gui --no-readline t.azr > run.log 2>&1)
awk '{ sub(/\r$/, "") }
     /^\[grad-check\] idx / { a = f = j = ""
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print a, f, j, $4 }' "$WORK/gfit/run.log" > "$WORK/gfit/rows"
nrow="$(awk 'END { print NR }' "$WORK/gfit/rows")"
ndop="$(awk '$4 == "doppler"' "$WORK/gfit/rows" | grep -c .)"
if [ "$ndop" -ne 2 ] || [ "$nrow" -lt 4 ]; then bad "(g) gradient check listed $nrow parameters, $ndop temperatures"
else
  read -r fd jac <<< "$(awk '{ a = $1 < 0 ? -$1 : $1; if (a > m) m = a
                               d = $1 - $2; if (d < 0) d = -d; if (d > w) w = d
                               e = $1 - $3; if (e < 0) e = -e; if (e > v) v = e }
                             END { printf "%.2e %.2e", w / m, v / m }' "$WORK/gfit/rows")"
  read -r dfd <<< "$(awk '$4 == "doppler" { d = ($1 - $2) / $2; if (d < 0) d = -d; if (d > w) w = d } END { printf "%.2e", w }' \
                     "$WORK/gfit/rows")"
  le "$fd" 1e-6 && ok "(g) gradient = finite differences for $nrow parameters, doppler_1 and doppler_2 among them ($fd of the largest; the temperatures to $dfd of their own)" \
    || bad "(g) gradient vs finite differences: $fd of the largest (temperatures $dfd of their own)"
  le "$jac" 1e-10 && ok "(g) Jacobian-path gradient = adjoint gradient ($jac of the largest)" \
    || bad "(g) Jacobian-path vs adjoint gradient: $jac of the largest"
fi

# (h) a prior
fit prior fit_prior --use-gradient
t="$(par prior doppler_1)"; dc="$(chi prior Total-Doppler-Chi-Squared)"
if [ -n "$t" ] && [ -n "$dc" ] && awk -v t="$t" -v s="$dc" -v re="$num" 'BEGIN { p = ((t - 300) / 20)^2; d = s - p; if (d < 0) d = -d
       exit !(t ~ re && s ~ re && t > 300 && t < 363 && d <= 1e-5 * (p + 1)) }'; then
  ok "(h) a 20 K prior holds T at $t K (data alone: 450 K); Total-Doppler-Chi-Squared $dc = ((T - 300)/20)^2"
else bad "(h) prior fit: doppler_1 '${t}', Total-Doppler-Chi-Squared '${dc}'"; tail -5 "$WORK/prior/run.log"; fi

# (i) malformed fit forms
for name in refuse_vary refuse_nosigma refuse_sigma refuse_tmax refuse_negT; do
  run "$name" "$name" 1
  if grep -q "a doppler block is 'doppler T'" "$WORK/$name/run.log" && [ ! -s "$WORK/$name/output/chiSquared.out" ]; then
    ok "(i) refused: $name ($(grep -m1 -o "a doppler block is.*" "$WORK/$name/run.log" | cut -c1-120))"
  else bad "(i) not refused: $name"; grep -m2 -E 'ERROR|WARN' "$WORK/$name/run.log"; fi
done

# (j) plain "doppler T" byte-identical to a reference build (optional)
if [ -n "${REF_AZURE2:-}" ]; then
  ref="$AZURE2_BIN"; AZURE2_BIN="$REF_AZURE2"
  for m in dop dopg dopT dopgT; do run "ref_$m" "$m" 1; done
  AZURE2_BIN="$ref"
  same=0; diff=0
  for m in dop dopg dopT dopgT; do
    for f in "$WORK/ref_$m"/output/AZUREOut_*.out "$WORK/ref_$m/output/chiSquared.out"; do
      if cmp -s "$f" "$WORK/$m/output/$(basename "$f")"; then same=$((same + 1)); else diff=$((diff + 1)); fi
    done
  done
  if [ "$diff" -eq 0 ] && [ "$same" -ge 12 ]; then ok "(j) plain doppler T: $same output files byte-identical to $REF_AZURE2"
  else bad "(j) plain doppler T: $diff of $((same + diff)) output files differ from $REF_AZURE2"; fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: free-gas Doppler broadening"; else echo "FAIL: $failures check(s)"; exit 1; fi
