#!/usr/bin/env bash
#
# The capture analyzing power couples the channel spin light particle first, like the rest of AZURE2
# (audit item 9, 2026-10-26).  Seyler and Weller's b_k (CNuc::CalcCaptureAnalyzingPower) couple it target
# first; the (-1)^(s - s') that takes their Eq. (21) to AZURE2's convention was missing, so for a target
# with spin the capture A_y had every channel-spin off-diagonal term sign-flipped.  Models from make_models.py.
#
#   (a) twolevel (n + 1/2- target, channel spins 0 and 1, two 1+ levels decaying by M1): the engine's A_y
#       against a magnetic-substate calculation from the engine's own reduced-width amplitudes (verify.py);
#       the old engine gives exactly -A_y;
#   (b) deuteron: an analyzing-power segment with a spin-1 projectile is refused with an error;
#   (c) tests/13N_capture_ay (spin-0 target, one channel spin) is unchanged: its expected chi2 still holds.
#
#   ./tests/capture_ay_coupling/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PYTHON="$(command -v python3 || command -v python || true)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/capture_ay_coupling.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run DIR MODEL MODE [FLAGS] -> $WORK/DIR
run() {
  local dir="$WORK/$1" model="$2" mode="$3" flags="${4:-}"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$model.azr" "$dir/t.azr"
  [ -d "$HERE/models/data" ] && cp -r "$HERE/models/data" "$dir/"
  # shellcheck disable=SC2086
  (cd "$dir" && printf "$mode\n\n\n7\n" | run_limited 300 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}

# ---- (a) the two-level model --------------------------------------------------------------------------------------
run twolevel twolevel 3 "--use-brune"
if [ -s "$WORK/twolevel/output/AZUREOut_aa=1_R=2.extrap" ] && [ -s "$WORK/twolevel/output/parameters.out" ]; then
  if [ -n "$PYTHON" ] && "$PYTHON" -c "import numpy, sympy" > /dev/null 2>&1; then
    "$PYTHON" "$HERE/verify.py" "$WORK/twolevel" || failures=$((failures + $?))
  else
    echo "  skip  (a): no python3 with numpy and sympy"
  fi
else
  bad "(a) twolevel wrote no output"; tail -8 "$WORK/twolevel/run.log"
fi

# ---- (b) a spin-1 projectile is refused -------------------------------------------------------------------------
run deuteron deuteron 3 "--use-brune"
if grep -q "ERROR: test segment #1: an analyzing power needs a spin-1/2 projectile" "$WORK/deuteron/run.log" &&
   [ ! -s "$WORK/deuteron/output/AZUREOut_aa=1_R=2.extrap" ]; then
  ok "(b) an analyzing-power segment with a spin-1 projectile is refused"
else
  bad "(b) the deuteron analyzing-power segment was not refused"; grep -i "error\|analyzing" "$WORK/deuteron/run.log" | head -3
fi

# ---- (c) the spin-0 case is unchanged ---------------------------------------------------------------------------
if [ -f "$HERE/../13N_capture_ay/expected/chiSquared.out" ]; then
  mkdir -p "$WORK/13N/output" "$WORK/13N/checks"
  cp "$HERE/../13N_capture_ay/13N_capture_ay.azr" "$WORK/13N/t.azr" && cp -r "$HERE/../13N_capture_ay/data" "$WORK/13N/"
  (cd "$WORK/13N" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  exp="$(awk -F, 'NR > 1 && $1 + 0 == 2 { gsub(/ /, "", $2); print $2 }' "$HERE/../13N_capture_ay/expected/chiSquared.out")"
  got="$(awk -F, 'NR > 1 && $1 + 0 == 2 { gsub(/ /, "", $2); print $2 }' "$WORK/13N/output/chiSquared.out" 2>/dev/null)"
  if [ -n "$got" ] && awk -v a="$got" -v b="$exp" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d <= 1e-3 * (b < 0 ? -b : b) + 1e-9) }'; then
    ok "(c) 13N_capture_ay (spin-0 target) unchanged: A_y segment chi2 $got (expected $exp)"
  else
    bad "(c) 13N_capture_ay A_y segment chi2 $got vs expected $exp"
  fi
fi

echo
if [ "$failures" -eq 0 ]; then echo "capture_ay_coupling: PASS"; else echo "capture_ay_coupling: FAIL ($failures)"; exit 1; fi
