#!/usr/bin/env bash
#
# The upgrade test projects T1 and T3 (tests/upgrade/T1, T3) against today's engine:
#   (a) golden/golden.py: the calculate, extrapolate and rate outputs equal the golden files to 1e-6 relative;
#   (b) recover.py: a MIGRAD fit from the moved start model returns the true parameters (|pull| <= 4 on the
#       committed noisy data, <= 0.05 on noise-free data) with chi2(fit) <= chi2(truth).
# Needs python3 with numpy; without them it says so and passes (as the pyazr tests skip).
#
#   ./tests/upgrade/check.sh path/to/AZURE2
set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
export AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
PY="$(command -v python3 || true)"
if [ -z "$PY" ] || ! "$PY" -c 'import numpy' 2>/dev/null; then
  echo "skip: python3 with numpy not found"; exit 0
fi
failures=0
"$PY" "$HERE/golden/golden.py" || failures=$((failures + 1))
"$PY" "$HERE/recover.py" || failures=$((failures + 1))
if [ "$failures" -eq 0 ]; then echo "upgrade T1/T3: ok"; else echo "upgrade T1/T3: FAIL"; exit 1; fi
