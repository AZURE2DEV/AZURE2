#!/usr/bin/env bash
#
# Secondary gamma rays following particle emission (UPOS segments; Brune and deBoer, PRC 102, 024628 (2020)),
# 13C + alpha in 17O (models from make_models.py):
#
#   (a) n + 16O(6.13, 3-): two exit channel spins (5/2, 7/2).  The E3 gamma angular distribution of a single level
#       against a substate calculation that uses neither Eq. (22) nor R_k (truth.py): the 5/2-7/2 interference with
#       both signs of the 7/2 amplitude, a zero 7/2 width, and l' = 1 and 3 together;
#   (b) alpha + 13C(3.68, 3/2-), M1/E2 mixing (delta = -0.094): Eqs. (4) and (7);
#   (c) two levels into one exit channel spin (J1 != J2 interference): unchanged from the engine before the fix
#       of the channel-spin pairing (expected/twolevel_pre_fix.out, dev 6c6932d);
#   (d) analytic gradient and Jacobian-path gradient = finite differences on a two-level model with UPOS data in
#       both channel spins (AZURE_GRAD_CHECK);
#   (e) with a second binary: the non-UPOS segments of the same models (angle-integrated and particle angular
#       distributions) are byte-identical between the two builds.
#
#   ./tests/secondary_gamma/check.sh path/to/AZURE2 [path/to/reference/AZURE2]

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2 [reference AZURE2]}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
REF_BIN="${2:-}"
[ -n "$REF_BIN" ] && REF_BIN="$(cd "$(dirname "$REF_BIN")" && pwd)/$(basename "$REF_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PYTHON="$(command -v python3 || command -v python || true)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/secondary_gamma.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run BIN DIR MODEL MODE [FLAGS] [ENV...] -> $WORK/DIR
run() {
  local bin="$1" dir="$WORK/$2" model="$3" mode="$4" flags="${5:-}"
  shift 5 2> /dev/null || shift $#
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$model.azr" "$dir/t.azr"
  cp -r "$HERE/models/data" "$dir/"
  # shellcheck disable=SC2086
  (cd "$dir" && for kv in "$@"; do export "$kv"; done
   printf "$mode\n\n\n7\n" | run_limited 300 "$bin" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}

MODELS="spin2_plus spin2_minus spin2_zero spin2_multil mixing twolevel"
for m in $MODELS; do
  run "$AZURE2_BIN" "$m" "$m" 1 ""
  ls "$WORK/$m/output/"AZUREOut_* > /dev/null 2>&1 || { bad "$m: no output"; tail -15 "$WORK/$m/run.log"; }
done

# ---- (a)-(c) ---------------------------------------------------------------------------------------------------
if [ -n "$PYTHON" ] && "$PYTHON" -c "import numpy" > /dev/null 2>&1; then
  "$PYTHON" "$HERE/verify.py" "$WORK" || failures=$((failures + $?))
else
  echo "  skip  (a)-(c): no python3 with numpy"
fi

# ---- (d) gradients ---------------------------------------------------------------------------------------------
run "$AZURE2_BIN" grad grad 2 "--use-gradient" AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_H=1e-7
awk '/^\[grad-check\] idx / { idx = $3; a = f = ""; j = "-"
       for (i = 4; i <= NF; i++) { if ($i == "analytic") a = $(i + 1); if ($i == "fd") f = $(i + 1); if ($i == "jacobian") j = $(i + 1) }
       print idx, a, f, j }' "$WORK/grad/run.log" | sort -n > "$WORK/grad.rows"
nrows="$(wc -l < "$WORK/grad.rows" | tr -d ' ')"
if [ "$nrows" -lt 10 ]; then
  bad "(d) gradient check listed $nrows parameters"; grep -i "grad\|error" "$WORK/grad/run.log" | head -5
else
  worst="$(awk '{ a = $2 < 0 ? -$2 : $2; if (a > m) m = a; g[NR] = $2; f[NR] = $3; j[NR] = $4 }
             END { for (i = 1; i <= NR; i++) { d = g[i] - f[i]; if (d < 0) d = -d; e = j[i] - f[i]; if (e < 0) e = -e
                                               if (d > w) w = d; if (j[i] != "-" && e > w) w = e }
                   printf "%.2e", w / m }' "$WORK/grad.rows")"
  if awk -v w="$worst" 'BEGIN { exit !(w + 0 <= 1e-5) }'; then
    ok "(d) adjoint and Jacobian-path gradients = finite differences, $nrows parameters, worst $worst of the largest"
  else
    bad "(d) gradient vs finite differences, worst $worst of the largest (> 1e-5)"
  fi
fi

# ---- (e) non-UPOS segments byte-identical to a reference build ----------------------------------------------
if [ -n "$REF_BIN" ]; then
  for m in spin2_plus spin2_multil mixing; do
    run "$REF_BIN" "ref_$m" "$m" 1 ""
    f="$(basename "$(ls "$WORK/$m/output/"AZUREOut_* | head -1)")"
    if cmp -s <(tail -n +14 "$WORK/$m/output/$f") <(tail -n +14 "$WORK/ref_$m/output/$f"); then
      ok "(e) $m: angle-integrated and particle angular distribution byte-identical to the reference build"
    else
      bad "(e) $m: non-UPOS segments differ from the reference build"
    fi
  done
fi

echo
[ "$failures" -eq 0 ] && echo "secondary_gamma: PASS" || echo "secondary_gamma: $failures FAILURE(S)"
exit "$failures"
