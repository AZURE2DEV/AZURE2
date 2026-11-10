#!/usr/bin/env bash
#
# What a beta-delayed spectrum (entrance pair of type 20) may carry (docs/upgrade/data_model.md S2).  Its energies are
# the light decay particle's (lab), converted with the decay pair's masses, so
#   (a) a segment without its decay pair as the exit (exit -1, total capture, data type 3) is refused at load: its energies would
#       be converted with the beta pair's masses, as a beam;
#   (b) a target integration, (c) a beam profile and (d) electron screening are refused: a decay has no beam and the
#       beta pair's charges are not a reaction's;
#   (e) a Gaussian resolution and (f) an energy shift of the measured particle energy run (a detector resolution and
#       a calibration; the sub-points move by the shift times the decay pair's c.m. factor).
#
#   ./tests/beta_delayed_alpha/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/beta_delayed_alpha.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
# variant NAME SEGMENT-LINE [TARGETINT-LINE]: the model with this data segment (and target line), menu 1.
variant() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$HERE/data" "$WORK/$1/"
  awk -v seg="$2" -v tgt="${3:-}" '
    /<segmentsData>/ { print; print seg; skip = 1; next }
    /<\/segmentsData>/ { skip = 0 }
    /<targetInt>/ { print; if (tgt != "") print tgt; next }
    !skip { print }' "$HERE/16N.azr" > "$WORK/$1/m.azr"
  (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
}
refused() { grep -q "$2" "$WORK/$1/run.log" && [ ! -s "$WORK/$1/output/chiSquared.out" ]; }
ran() { grep -q '^Total' "$WORK/$1/output/chiSquared.out" 2> /dev/null && ! grep -qi nan "$WORK/$1/output/chiSquared.out"; }
SEG="1 2 1 0 50 0 180 0 1 0 0 0 0 0 data/bda.dat 0 0"
TAIL='0 0 0 0 "" 0 0 0.04 20 50'

variant a "1 2 -1 0 50 0 180 3 1 0 0 0 0 0 data/bda.dat 0 0"
if refused a "needs its decay pair as the exit pair"; then ok "(a) exit -1: refused"
else bad "(a) exit -1: $(grep -m1 ERROR "$WORK/a/run.log")"; fi

variant b "$SEG" "1 \"1\" 10 0 0 1 1e18 \"10\" $TAIL"
if refused b "takes a detector resolution only"; then ok "(b) target integration: refused"
else bad "(b) target integration: $(grep -m1 ERROR "$WORK/b/run.log")"; fi

variant c "$SEG" "1 \"1\" 10 0 0 0 0 \"\" $TAIL beamprofile 1 1.0 0.1 0 1 0.05 3 0"
if refused c "takes a detector resolution only"; then ok "(c) beam profile: refused"
else bad "(c) beam profile: $(grep -m1 ERROR "$WORK/c/run.log")"; fi

variant d "$SEG" "1 \"1\" 200 0 0 0 0 \"\" $TAIL screening 1000 0 0"
if refused d "takes a detector resolution only"; then ok "(d) electron screening: refused"
else bad "(d) electron screening: $(grep -m1 ERROR "$WORK/d/run.log")"; fi

variant e "$SEG" "1 \"1\" 20 1 0.02 0 0 \"\" $TAIL"
if ran e; then ok "(e) Gaussian resolution: runs ($(grep '^Total' "$WORK/e/output/chiSquared.out"))"
else bad "(e) Gaussian resolution: $(grep -m1 ERROR "$WORK/e/run.log")"; fi

variant f "1 2 1 0 50 0 180 0 1 0 0 0.01 0 0 data/bda.dat 0 0"
if ran f; then ok "(f) energy shift: runs ($(grep '^Total' "$WORK/f/output/chiSquared.out"))"
else bad "(f) energy shift: $(grep -m1 ERROR "$WORK/f/run.log")"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
