#!/usr/bin/env bash
# Capture into a BROAD unbound final state: the line-shape integration ("linewindow W" on a capture segment,
# 2026-11-15).  The cross section is integrated over the final-state energy E_f within E_R +- W (keV), with the
# photon energy E - E_x(E_f), the real scattering wave of the final channel at E_f as the exterior of external
# capture, and the one-level line shape rho(E_f)/N_f as the weight (EPoint::CalculateLineWindow).  Without the
# keyword a segment keeps the narrow (stage-1) treatment unchanged.
#
# Models: models/dc39.azr, 12C + p with capture pair 3 into the 13N 1/2+ state at E_x = 2.3655 MeV (E_R = 422 keV,
# Gamma_p = 39 keV as in Rolfs & Azuma, NPA 227, 291 (1974)), a = 4.48 fm, every other width 0 (pure E1 direct
# capture from hard-sphere p waves); models/d52.azr, the same with a 5/2+ level at E_R = 1.60 MeV (Gamma_p 60 keV)
# that decays by E2 (0.05 eV) into the 2.3655 state.
#   (a) direct capture, windowed over unwindowed, at four energies for W = 39 and 117 keV: the independent
#       reference (R-matrix/AZURE2_code_audit/11-14-26_ec_unbound_broad/reference_broad.py: damped radial
#       integrals + Richardson, mpmath/ODE Coulomb functions, quadrature over the window): within 5e-5 / 3e-4;
#   (b) narrow limit: final width 39 eV, W = 20 Gamma: the reference (the window fraction of the line plus the
#       non-resonant continuum term) within 5e-5;
#   (c) internal capture alone (the E2 part of d52, angle-integrated, so E1 and E2 add incoherently): windowed over
#       unwindowed = int rho/N_f ((E - E_f)/(E - E_R))^5 dE_f (scipy quad), within 2e-5;
#   (d) --no-brune gives the Brune result;
#   (e) start-up: a NOTE per windowed segment; "linewindow" on an elastic segment stops with an ERROR.
#
#   ./tests/ec_unbound_broad/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_unbound_broad.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# model NAME BASE SED-EXPR KEYWORD ISDIFF E_CM...: BASE with SED-EXPR applied (portable sed, one substitution),
# one <segmentsTest> line per c.m. energy (entrance 1, exit 3, 90 deg if ISDIFF = 1), each ending in KEYWORD.
model() {
  local name="$1" base="$2" expr="$3" kw="$4" isdiff="$5"; shift 5
  mkdir -p "$WORK/$name/output" "$WORK/$name/checks"
  local segs="" e lab ang=0
  [ "$isdiff" = 1 ] && ang=90
  for e in "$@"; do
    lab=$(awk -v e="$e" 'BEGIN { printf "%.6f", e * 13.00728 / 12.0 }')
    segs="${segs}1 1 3 $lab $lab 0.01 $ang $ang 0 $isdiff$kw
"
  done
  # the segment lines go through the environment: BSD awk (macOS) refuses a newline in a -v value
  sed -e "$expr" "$HERE/models/$base" | EC_BROAD_SEGS="$segs" awk '
    /^<segmentsTest>/ { print; printf "%s", ENVIRON["EC_BROAD_SEGS"]; next }
    { print }' > "$WORK/$name/t.azr"
}
run() {  # NAME FLAGS: mode 3 of $WORK/NAME/t.azr; returns the engine's exit status
  (cd "$WORK/$1" && printf '3\nn\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" $2 --no-gui --no-readline t.azr > run.log 2>&1)
}
CAP3="output/AZUREOut_aa=1_R=3.extrap"
# column 4 of NAME's capture output, one value per line (blank lines between segments dropped)
values() { awk 'NF { print $4 }' "$WORK/$1/$CAP3" 2>/dev/null; }
# compare LABEL NUM DEN TOL REF...: (NUM/DEN)_i against REF_i, relative TOL
compare() {
  local label="$1" num="$2" den="$3" tol="$4"; shift 4
  local refs="$*"
  local out
  out=$(paste -d ' ' <(values "$num") <(values "$den") | awk -v refs="$refs" -v tol="$tol" '
    BEGIN { n = split(refs, r, " ") }
    { i++; v = $1 / $2; d = v / r[i] - 1; if (d < 0) d = -d; if (d > worst) worst = d
      printf "      %d: %.6f (ref %.6f)\n", i, v, r[i] }
    END { if (i != n) { print "      MISSING: " i " of " n " values"; exit 1 }
          printf "      worst relative deviation %.1e (tolerance %s)\n", worst, tol; exit (worst > tol) }')
  local status=$?
  echo "$out"
  if [ $status -eq 0 ]; then ok "$label"; else bad "$label"; fi
}

echo "(a) direct capture through +-1 and +-3 Gamma windows vs reference_broad.py"
EA="0.65 0.90 1.40 2.30"
model a_none dc39.azr 's/^$/&/' "" 1 $EA
model a_w1 dc39.azr 's/^$/&/' " linewindow 39" 1 $EA
model a_w3 dc39.azr 's/^$/&/' " linewindow 117" 1 $EA
for m in a_none a_w1 a_w3; do run $m "" || bad "(a) $m: engine exit status $?"; done
compare "(a) W = 39 keV" a_w1 a_none 5e-5 0.70442 0.70120 0.69904 0.69831
compare "(a) W = 117 keV" a_w3 a_none 3e-4 1.84144 1.03167 0.91742 0.89875

echo "(b) narrow limit: final width 39 eV, W = 20 Gamma"
model b_none dc39.azr 's/ 39000 / 39 /' "" 1 1.20 2.00
model b_w20 dc39.azr 's/ 39000 / 39 /' " linewindow 0.78" 1 1.20 2.00
for m in b_none b_w20; do run $m "" || bad "(b) $m: engine exit status $?"; done
compare "(b) narrow limit" b_w20 b_none 5e-5 1.06495 0.99331

echo "(c) internal E2 capture alone, angle-integrated"
ED="1.20 1.62 2.00"
model d_none d52.azr 's/^$/&/' "" 0 $ED
model d_w1 d52.azr 's/^$/&/' " linewindow 39" 0 $ED
model d_ecnone d52.azr 's/ 0 1 3 1 4 9 1 0 0.05 / 0 1 3 1 4 9 1 0 0 /' "" 0 $ED
model d_ecw1 d52.azr 's/ 0 1 3 1 4 9 1 0 0.05 / 0 1 3 1 4 9 1 0 0 /' " linewindow 39" 0 $ED
for m in d_none d_w1 d_ecnone d_ecw1; do run $m "" || bad "(c) $m: engine exit status $?"; done
# internal part = total - (external alone), windowed over unwindowed
mkdir -p "$WORK/d_int_w/output" "$WORK/d_int_n/output"
paste -d ' ' <(values d_w1) <(values d_ecw1) | awk '{ print 0, 0, 0, $1 - $2 }' > "$WORK/d_int_w/$CAP3"
paste -d ' ' <(values d_none) <(values d_ecnone) | awk '{ print 0, 0, 0, $1 - $2 }' > "$WORK/d_int_n/$CAP3"
compare "(c) internal capture over the line shape" d_int_w d_int_n 2e-5 0.698914 0.697808 0.697589

echo "(d) --no-brune"
model e_formal dc39.azr 's/^$/&/' " linewindow 39" 1 $EA
run e_formal "--no-brune" || bad "(d) engine exit status $?"
compare "(d) --no-brune = Brune" e_formal a_w1 1e-9 1 1 1 1

echo "(e) start-up"
if grep -q "NOTE: segment 1: capture into a broad unbound final state" "$WORK/a_w1/run.log"; then
  ok "(e) NOTE for a windowed segment"
else
  bad "(e) no NOTE for a windowed segment"
fi
mkdir -p "$WORK/f_elastic/output" "$WORK/f_elastic/checks"
awk '/^<segmentsTest>/ { print; print "1 1 1 1.2 1.2 0.01 90 90 0 1 linewindow 39"; next } { print }' \
  "$HERE/models/dc39.azr" > "$WORK/f_elastic/t.azr"
run f_elastic ""
if grep -q "ERROR: segment 1 (entrance 1, exit 1) has \"linewindow 39\", but its exit pair (key 1) is not a capture pair" \
  "$WORK/f_elastic/run.log" && [ ! -s "$WORK/f_elastic/output/AZUREOut_aa=1_R=1.extrap" ]; then
  ok "(e) linewindow on an elastic segment stops with an ERROR"
else
  bad "(e) linewindow on an elastic segment was not refused"
  tail -5 "$WORK/f_elastic/run.log"
fi

if [ $failures -eq 0 ]; then
  echo "PASS"
  exit 0
fi
echo "FAIL ($failures)"
exit 1
