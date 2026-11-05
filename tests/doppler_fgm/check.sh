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
#   ./tests/doppler_fgm/check.sh path/to/AZURE2

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

if [ "$failures" -eq 0 ]; then echo "PASS: free-gas Doppler broadening"; else echo "FAIL: $failures check(s)"; exit 1; fi
