#!/usr/bin/env bash
#
# The final state of an external-capture pair must be a level of the project (user guide, "Gamma-Ray
# Calculations": one level per gamma-ray pair, carrying the ANC and gamma widths).  Before 2026-10-26 a missing
# level was created silently with a fixed energy and FREE particle widths of 0.1 MeV^1/2 that appeared in no
# .azr line, and -- since it was only added when external capture was on (CLI: modes with a capture segment;
# API: always) -- every later level's energy_N / width_N_c name moved by one between runs of the same project,
# so a param.sav landed on the wrong levels without a message (audit 2026-10-25, item 4).
#
# Models (p + 16O -> 17F(0.495 MeV, 1/2+) + gamma, E1 external capture, from tests/ec_m1_exact/models/e1.azr):
#   listed        : the final state is a level -> runs, chi2 as expected/listed_chi2
#   missing_level : the 1/2+ group exists (a 30 MeV level) but the 0.495 MeV level is not listed -> refused
#   missing_group : no 1/2+ level at all -> refused
# "Refused" = the ERROR line naming the pair, J^pi and E_x, and no chiSquared.out.
#
#   ./tests/ec_final_state_listed/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_final_state_listed.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
ERR="ERROR: external capture into pair 3 needs its final state"

run() {  # MODEL -> $WORK/MODEL, mode 1
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$1.azr" "$dir/t.azr"; cp -r "$HERE/data" "$dir/"
  (cd "$dir" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}
refused() {  # MODEL LABEL
  if grep -q "$ERR" "$WORK/$1/run.log" && [ ! -s "$WORK/$1/output/chiSquared.out" ]; then ok "$2"
  else bad "$2: $(grep -c "$ERR" "$WORK/$1/run.log") ERROR lines, chiSquared.out $([ -s "$WORK/$1/output/chiSquared.out" ] && echo written || echo absent)"; fi
}

echo "ec_final_state_listed: $AZURE2_BIN"
run listed; run missing_level; run missing_group
chi2=$(awk '/^Total-Chi-Squared:/ { print $2 }' "$WORK/listed/output/chiSquared.out" 2>/dev/null)
exp=$(cat "$HERE/expected/listed_chi2")
if [ -n "$chi2" ] && awk -v a="$chi2" -v b="$exp" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(d <= 1e-3 * b) }'; then
  ok "(a) final state listed: runs, chi2 $chi2 (expected $exp)"
else bad "(a) final state listed: chi2 '$chi2' vs expected $exp"; fi
refused missing_level "(b) final state not listed although its J^pi group exists: refused with the ERROR line, no output"
refused missing_group "(c) final state's J^pi group absent: refused with the ERROR line, no output"
if [ "$failures" -eq 0 ]; then echo "ec_final_state_listed: PASS"; else echo "ec_final_state_listed: FAIL ($failures)"; exit 1; fi
