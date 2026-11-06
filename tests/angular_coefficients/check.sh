#!/usr/bin/env bash
#
# Angular-coefficient data segments (isDiff 9, 2026-10-26): the Legendre coefficients a_k = A_k / A_0 of a reaction's
# c.m. angular distribution (or the absolute A_k with "absolute 1") fitted as data.  Models from make_models.py
# (13C + alpha -> n + 16O(3-), two levels, no external capture); data written by verify.py from the model's own
# extrapolation, so the engine's convention is pinned by (a) and the plumbing by the rest:
#
#   (a) the isDiff 3 coefficients equal the Legendre projection of the c.m. differential cross section, and
#       (b) its c_0 equals sigma_int / 4 pi;
#   (c) a normalised and an absolute coefficient data segment reproduce their data (mode 1 chi2 ~ 0);
#   (d) a fixed 0.05 MeV energy shift moves the coefficients to E + 0.05;
#   (e) under a 30 keV beam spread the folded coefficients equal the projection of the folded distribution
#       (the absolute A_k are folded like cross sections, then the ratio is formed);
#   (f) the analytic gradient (AZURE_GRAD_CHECK) equals finite differences with normalised, absolute and
#       folded coefficient segments and free widths;
#   (g) refused: a free normalization on normalised coefficients, an elastic pair, isDiff 9 in <segmentsTest>.
#
#   ./tests/angular_coefficients/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PYTHON="$(command -v python3 || command -v python || true)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/angular_coefficients.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
if [ -z "$PYTHON" ] || ! "$PYTHON" -c "import numpy" > /dev/null 2>&1; then echo "  skip  no python3 with numpy"; echo "angular_coefficients: SKIP"; exit 0; fi

# run MODEL MODE [FLAGS] [ENV...] -> $WORK/MODEL (data/ = $WORK/data)
run() {
  local model="$1" mode="$2" flags="${3:-}"; shift 3 2> /dev/null || shift $#
  local dir="$WORK/$model"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$model.azr" "$dir/t.azr"
  [ -d "$WORK/datafiles" ] && cp -r "$WORK/datafiles" "$dir/data"
  (cd "$dir" && for kv in "$@"; do export "$kv"; done
   printf "$mode\n\n\n7\n" | run_limited 600 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}

# ---- (a), (b): the convention, and the data files -------------------------------------------------------------------
run extrap 3 "--use-brune"
[ -s "$WORK/extrap/output/AZUREOut_aa=2_R=5.extrap" ] || { bad "extrap: no output"; tail -5 "$WORK/extrap/run.log"; exit 1; }
"$PYTHON" "$HERE/verify.py" prepare "$WORK/extrap" "$WORK/datafiles" || failures=$((failures + 1))

# ---- (c) data segments ----------------------------------------------------------------------------------------------
run data 1 "--use-brune"
if [ -s "$WORK/data/output/chiSquared.out" ]; then
  chi="$(grep "Total-Chi" "$WORK/data/output/chiSquared.out" | awk '{print $2}')"
  if awk -v c="$chi" 'BEGIN { exit !(c + 0 <= 1e-8) }'; then ok "(c) both coefficient data segments: total chi2 $chi"; else bad "(c) total chi2 $chi (expected ~0)"; grep -i "error" "$WORK/data/run.log" | head -3; fi
  "$PYTHON" "$HERE/verify.py" data "$WORK/data" "$WORK/datafiles" || failures=$((failures + 1))
else
  bad "(c) the data run wrote no chiSquared.out"; grep -i "error" "$WORK/data/run.log" | head -5
fi

# ---- (d) energy shift -----------------------------------------------------------------------------------------------
run shift 1 "--use-brune"
"$PYTHON" "$HERE/verify.py" shift "$WORK/shift" "$WORK/extrap" || failures=$((failures + 1))

# ---- (e) beam-energy spread -----------------------------------------------------------------------------------------
run conv 1 "--use-brune"
run conv_diff 3 "--use-brune"
"$PYTHON" "$HERE/verify.py" conv "$WORK/conv" "$WORK/conv_diff" || failures=$((failures + 1))

# ---- (f) gradient ---------------------------------------------------------------------------------------------------
# --no-grid-replacement: the finite differences otherwise include the motion of the convolved segment's adaptive
# sub-point grid with the widths (a ~1e-5 effect the analytic gradient leaves out by design, as everywhere else)
run fit 2 "--use-brune --use-gradient --no-grid-replacement" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_H=1e-7
awk '/^\[grad-check\] idx / { idx = $3; a = f = ""; j = "-"
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print idx, a, f, j }' "$WORK/fit/run.log" | sort -n > "$WORK/grad.rows"
nrows="$(wc -l < "$WORK/grad.rows" | tr -d ' ')"
if grep -q "adjoint bailed" "$WORK/fit/run.log"; then
  bad "(f) the adjoint bailed on a coefficient segment"; grep "adjoint bailed" "$WORK/fit/run.log" | head -2
elif [ "$nrows" -lt 5 ]; then
  bad "(f) gradient check listed $nrows parameters"; grep -i "grad\|error" "$WORK/fit/run.log" | head -5
else
  worst="$(awk '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a; g[NR] = $2; f[NR] = $3; j[NR] = $4 }
             END { for (i = 1; i <= NR; i++) { d = g[i] - f[i]; if (d < 0) d = -d; e = j[i] - f[i]; if (e < 0) e = -e
                                               if (d > w) w = d; if (j[i] != "-" && e > w) w = e }
                   printf "%.2e", w / m }' "$WORK/grad.rows")"
  if awk -v w="$worst" 'BEGIN { exit !(w + 0 <= 1e-5) }'; then
    ok "(f) adjoint and Jacobian-path gradients = finite differences, $nrows parameters, worst $worst of the largest"
  else
    bad "(f) gradient vs finite differences, worst $worst of the largest (> 1e-5)"
  fi
fi

# ---- (g) refusals ---------------------------------------------------------------------------------------------------
run bad_norm 1 "--use-brune"
if grep -q "ERROR: segment #1 (angular coefficients): normalised coefficients" "$WORK/bad_norm/run.log" && [ ! -s "$WORK/bad_norm/output/chiSquared.out" ]; then ok "(g) a free normalization on normalised coefficients is refused"; else bad "(g) free norm not refused"; fi
run bad_elastic 1 "--use-brune"
if grep -q "ERROR: segment #1 (angular coefficients): an elastic distribution" "$WORK/bad_elastic/run.log"; then ok "(g) an elastic pair is refused"; else bad "(g) elastic pair not refused"; fi
run bad_test 3 "--use-brune"
if grep -q "ERROR: Test segment #1 uses isDiff 9" "$WORK/bad_test/run.log"; then ok "(g) isDiff 9 in <segmentsTest> is refused"; else bad "(g) isDiff 9 test segment not refused"; fi

echo
if [ "$failures" -eq 0 ]; then echo "angular_coefficients: PASS"; else echo "angular_coefficients: FAIL ($failures)"; exit 1; fi
