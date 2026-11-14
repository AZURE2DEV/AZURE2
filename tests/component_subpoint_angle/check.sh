#!/usr/bin/env bash
#
# The sub-points of a SUM/RATIO component sit at the component's own c.m. angle (audit 2026-10-25, item 32; 2026-10-27).
#
# Component points are copies of the base segment's points.  SetUpTargetEffects laid out a target effect's
# sub-points from the copied point's c.m. angle, before InitializeComponentSegments converted the component point's
# angle with its own pair, and never revisited them: a component with a different exit pair (or a fixed-angle RATIO
# denominator, whose raw lab angle was used) was integrated at the base pair's angle.  Now the sub-points are put at
# the converted angle before their Legendre polynomials are computed.
#
# Model: tests/natural_target/models/M.azr (p + 12C and p + 13C in one project), test segments at 90 deg lab:
#   P:   1 -> 1 and 2 -> 2 alone;
#   Q:   the SUM 12C(p,p)12C + 13C(p,p)13C (1 -> 1 and 2 -> 2), no target effect;
#   Qc:  the same SUM with a beam spread of sigma = 0.1 keV lab (a sub-point effect that changes nothing but the
#        sub-point bookkeeping).
#   (a) Qc equals Q to 1e-6 (the second term's sub-points at the p + 13C angle, 94.4 deg, not the p + 12C 94.8 deg);
#   (b) in the lab output (kinflag.dat = 1) Q equals P1 + P2 to 1e-5 (2e-6 observed): the SUM of lab-differential terms is their lab
#       sum (the default c.m. output of a SUM is the lab sum taken to the BASE pair's c.m. frame, so there the
#       second term carries the Jacobian ratio J_base / J_term, 1.0005 here -- by design, not a defect).
# dev 510e6d1: (a) differs by ~1 %.
#
#   ./tests/component_subpoint_angle/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../natural_target/models/M.azr"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/component_subpoint_angle.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
# model NAME TARGETS EXTRAPOLATION...: M.azr with the given extrapolation records and <targets> records, run in mode 3
model() {
  local name="$1" targets="$2"; shift 2
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"
  [ -d "$(dirname "$SRC")/data" ] && cp -r "$(dirname "$SRC")/data" "$WORK/$name/"
  echo 1 > "$WORK/$name/kinflag.dat"   # lab-frame output files
  # the records go through the environment: BSD awk (macOS) refuses a newline in -v ("newline in string")
  RECORDS="$(printf '%s\n\n' "$@")" TARGETS="$targets" awk '
    /^<extrapolations>$/ { print; print ENVIRON["RECORDS"]; skip = 1; next }
    /^<\/extrapolations>$/ { skip = 0 }
    /^<targets>$/ { print; if (ENVIRON["TARGETS"] != "") print ENVIRON["TARGETS"]; skip = 1; next }
    /^<\/targets>$/ { skip = 0 }
    !skip { print }' "$SRC" > "$WORK/$name/m.azr"
  (cd "$WORK/$name" && run_limited 600 "$AZURE2_BIN" extrapolate --no-gui m.azr > run.log 2>&1; echo "exit $?" >> run.log)
}
# column 4 (the cross section, lab frame here) of a test-segment output file, blank lines dropped
col() { awk 'NF { print $4 }' "$1"; }
maxrel() { paste "$1" "$2" | awk '{ d = $1 - $2; if (d < 0) d = -d; s = ($2 < 0 ? -$2 : $2); r = (s > 0 ? d / s : d); if (r > m) m = r } END { printf "%.3e", m + 0 }'; }

# 1.0 to 3.0 MeV lab in 0.1 MeV steps (accumulated, as the old test lines), dsigma/dOmega at 90 deg lab
E="list $(awk 'BEGIN { for (e = 1.0; e <= 3.0; e += 0.1) printf "%.17g ", e }')MeV lab"
x() {  # x NAME REACTION [MORE FIELDS]: an extrapolation record
  local name="$1" reaction="$2"; shift 2
  printf 'extrapolation %s\n    reaction       %s\n    quantity       dsigma_domega\n    unit           b/sr\n' "$name" "$reaction"
  printf '    energies       %s\n    angles         list 90 deg lab' "$E"
  for f in "$@"; do printf '\n    %s' "$f"; done
}
L1="$(x x1 '12C(p,p)12C')"
L2="$(x x2 '13C(p,p)13C')"
LS="$(x x1 '12C(p,p)12C + 13C(p,p)13C')"
LSc="$(x x1 '12C(p,p)12C + 13C(p,p)13C' 'beam           B' 'quadrature     points 50 width_multiplier 20 per_width 50')"
model P "" "$L1" "$L2"
model Q "" "$LS"
model Qc "$(printf 'beam B\n    spread         0.0001 MeV lab')" "$LSc"
for n in P Q Qc; do grep -q "^exit 0" "$WORK/$n/run.log" || bad "$n did not run: $(tail -1 "$WORK/$n/run.log"); $(grep -m1 ERROR "$WORK/$n/run.log")"; done

col "$WORK/Q/output/AZUREOut_aa=1_R=1.extrap" > "$WORK/q.txt"
col "$WORK/Qc/output/AZUREOut_aa=1_R=1.extrap" > "$WORK/qc.txt"
r="$(maxrel "$WORK/q.txt" "$WORK/qc.txt")"
if [ -s "$WORK/q.txt" ] && [ "$(wc -l < "$WORK/q.txt")" = "$(wc -l < "$WORK/qc.txt")" ] && awk -v r="$r" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && r < 1e-6) }'; then
  ok "(a) the SUM with a 0.1 keV convolution equals the SUM without it (max rel. $r, $(wc -l < "$WORK/q.txt" | tr -d ' ') energies)"
else bad "(a) convolved SUM vs plain SUM: max rel. $r ($(wc -l < "$WORK/q.txt" | tr -d ' ') vs $(wc -l < "$WORK/qc.txt" | tr -d ' ') energies)"; fi

paste <(col "$WORK/P/output/AZUREOut_aa=1_R=1.extrap") <(col "$WORK/P/output/AZUREOut_aa=2_R=2.extrap") | awk '{ print $1 + $2 }' > "$WORK/p.txt"
r="$(maxrel "$WORK/q.txt" "$WORK/p.txt")"
if [ -s "$WORK/p.txt" ] && awk -v r="$r" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && r < 1e-5) }'; then ok "(b) the lab SUM equals the sum of its terms' lab cross sections (max rel. $r)"
else bad "(b) SUM vs P1 + P2: max rel. $r"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
