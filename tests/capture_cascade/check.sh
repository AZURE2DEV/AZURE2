#!/usr/bin/env bash
#
# Capture cascades: the secondary gamma ray of a capture final state after an UNOBSERVED primary capture photon
# (the field secondary_gamma on a capture reaction).  The extension of Brune and deBoer, PRC 102, 024628 (2020) named
# in that paper's conclusions: an unobserved photon summed over direction and helicity leaves no interference between
# multipoles and orients the final state like a spin-0 particle with l' = L.  Models: models/ (dev's, written by
# AZURE2 --convert); expectations: expected/.
#
#   (a) resonant capture (13C + alpha): one level and one multipole; two primary multipoles of one level (incoherent)
#       with a mixed secondary transition; two levels of different J through one multipole, both signs;
#   (b) external capture alone (p + a 1/2- target): into an s state (isotropic only as a coherent sum over J) and a
#       d -> p single-particle transition;
#   (c) a resonance and external capture together, at three energies across the resonance;
#   (d) an angle-integrated data set with secondary_gamma on a capture pair = the capture cross section;
#   (e) analytic and Jacobian-path gradients = finite differences (AZURE_GRAD_CHECK), free level and ANC.
#
# Expectations are magnetic-substate sums that use neither the UPOS coefficient nor R_k; for (a) and (b) the PRIMARY
# distribution of the same amplitudes is checked too, which ties the conventions to code that existed before.
#
#   ./tests/capture_cascade/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
PYTHON="$(command -v python3 || command -v python || true)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/capture_cascade.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run DIR MODEL COMMAND [VAR=VALUE ...] -> $WORK/DIR
run() {
  local dir="$WORK/$1" model="$2" command="$3"
  shift 3
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$model.azr" "$dir/t.azr"
  cp -R "$HERE/models/data" "$dir/"
  (cd "$dir" && env "$@" $(guard_command 300) "$AZURE2_BIN" $command t.azr < /dev/null > run.log 2>&1)
}

MODELS="single_e1 multipoles twoj_plus twoj_minus ec_s ec_p ecres_both"
for m in $MODELS; do
  run "$m" "$m" calculate
  ls "$WORK/$m/output/"AZUREOut_* > /dev/null 2>&1 || { bad "$m: no output"; tail -15 "$WORK/$m/run.log"; }
done

# ---- (a)-(d) ---------------------------------------------------------------------------------------------------
if [ -n "$PYTHON" ] && "$PYTHON" -c "import numpy" > /dev/null 2>&1; then
  "$PYTHON" "$HERE/verify.py" "$WORK" || failures=$((failures + $?))
else
  echo "  skip  (a)-(d): no python3 with numpy"
fi

# ---- (e) gradients ---------------------------------------------------------------------------------------------
run grad grad "fit --use-gradient" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_H=1e-7
awk '/^\[grad-check\] idx / { idx = $3; a = f = ""; j = "-"
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print idx, a, f, j }' "$WORK/grad/run.log" | sort -n > "$WORK/grad.rows"
nrows="$(wc -l < "$WORK/grad.rows" | tr -d ' ')"
if grep -q "adjoint bailed" "$WORK/grad/run.log"; then
  bad "(e) the adjoint bailed on a cascade data set"; grep "adjoint bailed" "$WORK/grad/run.log" | sed -n '1,2p'
elif [ "$nrows" -lt 4 ]; then
  bad "(e) gradient check listed $nrows parameters"; grep -i "grad\|error" "$WORK/grad/run.log" | sed -n '1,5p'
else
  worst="$(awk '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a; g[NR] = $2; f[NR] = $3; j[NR] = $4 }
             END { for (i = 1; i <= NR; i++) { d = g[i] - f[i]; if (d < 0) d = -d; e = j[i] - f[i]; if (e < 0) e = -e
                                               if (d > w) w = d; if (j[i] != "-" && e > w) w = e }
                   printf "%.2e", w / m }' "$WORK/grad.rows")"
  if awk -v w="$worst" 'BEGIN { exit !(w ~ /^ *[-+]?\.?[0-9]/ && w + 0 <= 1e-5) }'; then
    ok "(e) adjoint and Jacobian-path gradients = finite differences, $nrows parameters, worst $worst of the largest"
  else
    bad "(e) gradient vs finite differences, worst $worst of the largest (> 1e-5)"
  fi
fi

echo
[ "$failures" -eq 0 ] && echo "capture_cascade: PASS" || echo "capture_cascade: $failures FAILURE(S)"
[ "$failures" -eq 0 ]
