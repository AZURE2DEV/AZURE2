#!/usr/bin/env bash
#
# s-wave scattering lengths (ScatteringLength.cpp): the threshold evaluation against independent values,
# against AZURE2's own phase shifts through the effective-range expansion, across parametrizations, the
# <scatteringLengths> chi2 term, its analytic gradient, a fit with it, and the block's error messages.
# See README.md for the models and the reference values.
#
#   ./tests/scattering_length/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PYTHON="$(command -v python3 || command -v python || true)"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/scattering_length.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME MODEL MODE [FLAGS...] -> $WORK/NAME (output/, run.log); env vars pass through
run() {
  local name="$1" model="$2" mode="$3"
  shift 3
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$HERE/data" "$dir/data"
  cp "$HERE/$model" "$dir/m.azr"
  (cd "$dir" && printf '%s\n\n\n' "$mode" | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline "$@" m.azr > run.log 2>&1)
}
# a0 NAME PAIR JPI [im]: Re (or Im) a0 of one row of scatteringLengths.out
a0() {
  awk -v p="$2" -v j="$3" -v col="${4:-re}" '$1 == p && $3 == j { print (col == "im" ? $6 : $5); exit }' \
    "$WORK/$1/output/scatteringLengths.out"
}
# close A B TOL: |A - B| <= TOL |B|
close() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; exit !(a != "" && d <= t * m) }'; }
check() {  # check LABEL VALUE REFERENCE TOL
  if close "$2" "$3" "$4"; then ok "$1: $2 (reference $3)"; else bad "$1: ${2:-missing}, reference $3 (rel. tol. $4)"; fi
}

# 1. Threshold evaluation against independent values (README.md): mpmath Coulomb functions with the
#    Brune level matrix at E = 0 (proto_*.py of the archive), and the neutral closed form a (1 - g^2/E).
run he3he4 he3he4.azr 3
run be7p be7p.azr 3
run n16o n16o.azr 3
check "3He+4He 1/2+ (Paneru 2024: 33.10 fm)" "$(a0 he3he4 1 1/2+)" 33.099574 1e-5
check "7Be+p 1-, a01" "$(a0 be7p 1 1-)" 17.06776438 1e-5
check "7Be+p 2-, a02" "$(a0 be7p 1 2-)" -5.018395912 1e-5
check "n+16O 1/2+ (neutral, a (1 - g^2/E))" "$(a0 n16o 1 1/2+)" "$(awk 'BEGIN { printf "%.10f", 4.0 * (1 - 0.64 / (6.0 - 4.1431)) }')" 1e-7
im="$(a0 be7p 2 1- im)"
if awk -v x="$im" 'BEGIN { exit !(x + 0 < -1e-3) }'; then ok "7Be*+p 1-: complex (Im a0 = $im fm, the 7Be+p channel is open)"
else bad "7Be*+p 1-: Im a0 = ${im:-missing}, expected negative (open 7Be+p channel)"; fi

# 2. The same lengths from AZURE2's own phase shifts: effective-range expansion fitted in k^2 (ere.py).
if [ -n "$PYTHON" ]; then
  ere() { "$PYTHON" "$HERE/ere.py" "$WORK/$1/output/AZUREOut_aa=1_R=1.extrap" "${@:2}"; }
  check "ERE from phase shifts, 3He+4He" "$(ere he3he4 1 4 3.0160293 4.0026032)" "$(a0 he3he4 1 1/2+)" 2e-4
  check "ERE from phase shifts, 7Be+p 2-" "$(ere be7p 1 4 1.00783 7.01693)" "$(a0 be7p 1 2-)" 2e-4
  check "ERE from phase shifts, 7Be+p 1-" "$(ere be7p 2 4 1.00783 7.01693)" "$(a0 be7p 1 1-)" 2e-4
  check "ERE from phase shifts, n+16O" "$(ere n16o 1 0 1.00866492 15.99491462)" "$(a0 n16o 1 1/2+)" 2e-4
else
  echo "  skip  effective-range checks (no python3)"
fi

# 3. Independent of the parametrization: the same observed levels in the standard formalism.
run be7p_std be7p.azr 3 --no-brune
for j in 1- 2-; do check "7Be+p $j, standard parametrization = Brune" "$(a0 be7p_std 1 $j)" "$(a0 be7p 1 $j)" 1e-8; done
check "7Be*+p 1- (Im), standard = Brune" "$(a0 be7p_std 2 1- im)" "$(a0 be7p 2 1- im)" 1e-8

# 4. The <scatteringLengths> chi2: every kind of line (a01, a02, the coherent combination, a complex
#    length with both parts, an inactive line).  Coherent = 3/8 a01 + 5/8 a02 for p (1/2) + 7Be (3/2).
run fit7 be7p_fit.azr 1
SL="$WORK/fit7/output/scatteringLengths.out"
coh="$(awk '/coherent/ { for (i = 1; i <= NF; i++) if ($i == "model") { print $(i + 1); exit } }' "$SL")"
check "coherent combination = 3/8 a01 + 5/8 a02" "$coh" "$(awk -v a="$(a0 fit7 1 1-)" -v b="$(a0 fit7 1 2-)" 'BEGIN { printf "%.8f", 3 / 8 * a + 5 / 8 * b }')" 1e-5
total="$(awk '/total scattering-length chi2/ { print $NF }' "$SL")"
sum="$(awk '/^  line/ { s += $NF } END { printf "%.8g", s }' "$SL")"
check "chi2 total = sum of the lines (inactive line 0)" "$total" "$sum" 1e-6
check "chi2 total in chiSquared.out" "$(sed -n 's/.*Total-Scattering-Length-Chi-Squared: \([^ ]*\).*/\1/p' "$WORK/fit7/output/chiSquared.out")" "$total" 1e-5
check "chi2 total, reference" "$total" 48.3471 1e-4

# 5. Analytic gradient (level-matrix adjoint at the thresholds) against finite differences, in the
#    three parametrizations; AZURE_GRAD_CHECK compares every free energy and width.
for form in "" --no-brune --use-park; do
  AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=40 run "grad${form}" be7p_fit.azr 2 --use-gradient $form
  log="$WORK/grad${form}/run.log"
  worst="$(awk '/\[grad-check\] idx/ { for (i = 1; i <= NF; i++) if ($i == "rel" || $i == "relJ") if ($(i + 1) + 0 > m) m = $(i + 1) + 0 } END { print m + 0 }' "$log")"
  n="$(grep -c '\[grad-check\] idx' "$log")"
  if [ "$n" -ge 8 ] && awk -v w="$worst" 'BEGIN { exit !(w < 1e-6) }'; then ok "gradient ${form:---use-brune}: $n parameters, worst relative difference $worst"
  else bad "gradient ${form:---use-brune}: $n parameters checked, worst relative difference ${worst} (limit 1e-6)"; tail -5 "$log"; fi
done

# 6. A fit with the constraints (least squares: the constraints are residual rows).
run lm be7p_fit.azr 2 --use-lm
fitted="$(awk '/total scattering-length chi2/ { print $NF }' "$WORK/lm/output/scatteringLengths.out")"
data="$(sed -n 's/^Total-Chi-Squared: \([^ ]*\).*/\1/p' "$WORK/lm/output/chiSquared.out")"
if awk -v a="$fitted" -v b="$data" 'BEGIN { exit !(a != "" && b != "" && a + b < 6) }'; then ok "--use-lm fit: data chi2 $data + scattering-length chi2 $fitted (start 227.0 + 48.3)"
else bad "--use-lm fit: data chi2 ${data:-missing} + scattering-length chi2 ${fitted:-missing}, expected below 6"; fi

# 7. Errors in the block stop the run with a message.
refuse() {  # refuse LABEL LINE EXPECTED-TEXT
  local dir="$WORK/refuse_$4"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$HERE/data" "$dir/data"
  sed "s|^1 1 1 17.0 0.5\$|$2|" "$HERE/be7p_fit.azr" > "$dir/m.azr"
  (cd "$dir" && printf '1\n\n\n' | run_limited 120 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
  if grep -q "ERROR: <scatteringLengths> line 1: $3" "$dir/run.log" && [ ! -s "$dir/output/chiSquared.out" ]; then ok "refused: $1"
  else bad "not refused: $1"; grep -i error "$dir/run.log" | head -3; fi
}
refuse "channel spin 3 of p+7Be" "1 1 3 17.0 0.5" "pair 1 has no s-wave channel with channel spin 3" a
refuse "unknown pair key" "1 9 1 17.0 0.5" "no particle pair with key 9" b
refuse "missing sigma" "1 1 1 17.0" "expected" c
refuse "imaginary part without its sigma" "1 1 1 17.0 0.5 -0.3" "an imaginary part needs its uncertainty too" d

if [ "$failures" -eq 0 ]; then echo "PASS scattering_length"; exit 0; fi
echo "FAIL scattering_length ($failures)"; exit 1
