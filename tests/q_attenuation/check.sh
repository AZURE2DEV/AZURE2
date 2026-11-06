#!/usr/bin/env bash
#
# Finite-geometry attenuation coefficients (Q_K, <targetInt> line) on a PARTICLE exit (EPoint::CalcLegendreP):
# the detector averages the LAB angular distribution, so Q_K multiplies its lab Legendre coefficients.  Models from
# make_models.py: 13C(a,n)16O*(6.13) near threshold (g = v_cm/v_n* = 0.33-0.37), c.m. differential test segment.
#
#   (a) 40 coefficients all 1 (q1.azr) give output byte-identical to no coefficients (plain.azr);
#   (b) 40 coefficients of a uniform 15-degree cone (q.azr) agree with a brute-force average of the lab distribution
#       over that cone (verify.py, numpy) -- no Legendre algebra of the attenuation in the reference;
#   (c) only Q_0..Q_4 of the cone (q5.azr, as archive fits give them) agree with W* + J sum_K (Q_K - 1) b_K P_K(x_lab),
#       b_K projected in the lab variable.
# Before 2026-10-25 (a) failed (the B = U Q vector changed the result with every Q_K = 1), (b) was off by up to a few
# per cent and more than 20 coefficients read past the end of an array.
#
#   ./tests/q_attenuation/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PYTHON="$(command -v python3 || command -v python || true)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/q_attenuation.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
OUT=output/AZUREOut_aa=2_R=5.extrap

# Byte-identical files without cmp (not in every MSYS2 image);  stripped so a
# CRLF-writing Windows build compares like the others.
same_file() { [ -f "$2" ] && [ "$(tr -d '' < "$1")" = "$(tr -d '' < "$2")" ]; }

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {   # run MODEL: extrapolate (mode 3) in $WORK/MODEL
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$1.azr" "$dir/t.azr"
  cp -r "$HERE/models/data" "$dir/"
  (cd "$dir" && printf '3\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}

echo "q_attenuation: $AZURE2_BIN"
for m in plain q q1 q5; do
  run "$m"
  [ -s "$WORK/$m/$OUT" ] || { bad "$m: no output"; tail -15 "$WORK/$m/run.log"; }
done
[ "$failures" -eq 0 ] || { echo "q_attenuation: FAIL ($failures)"; exit 1; }

if same_file "$WORK/plain/$OUT" "$WORK/q1/$OUT"; then
  ok "(a) every Q_K = 1 leaves the output byte-identical"
else
  bad "(a) every Q_K = 1 changes the output"
fi

if [ -n "$PYTHON" ] && "$PYTHON" -c "import numpy" > /dev/null 2>&1; then
  "$PYTHON" "$HERE/verify.py" "$WORK" || failures=$((failures + $?))
else
  echo "  skip  (b), (c): no python3 with numpy"
fi

if [ "$failures" -eq 0 ]; then echo "q_attenuation: PASS"; else echo "q_attenuation: FAIL ($failures)"; exit 1; fi
