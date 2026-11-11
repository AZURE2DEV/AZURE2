#!/usr/bin/env bash
#
# A lab angle beyond the kinematic maximum is an error, not a NaN (audit 2026-10-25, item 31; 2026-10-27).
#
# EPoint::ConvertLabAngle models only the forward branch of the lab -> c.m. map.  Where two c.m. angles reach one lab
# angle (a reaction near the threshold of an endothermic exit, 7Li(p,n)7Be at E_p = 1.88-1.90 MeV: theta_max 21-45
# deg) or where the projectile is heavier than the target (elastic, theta_max = asin(m2/m1)), a lab angle past
# theta_max gave NaN: the reaction branch printed one line in the log and asserted (compiled out in release builds),
# the elastic branch said nothing; chiSquared.out and the extrapolation files then carried nan and the run exited 0.
# Now the point is reported ("ERROR: ... beyond the kinematic maximum") and the run stops before any output.
#
#   (a) tests/7Li_p_ay + a lab-differential TEST segment 7Li(p,n) at 60 deg for E_p 1.89-1.90 MeV, mode 3: ERROR, no
#       extrap file;
#   (b) the same as a DATA segment (two points), mode 1: ERROR, no chiSquared.out;
#   (c) tests/13N with pair 1 written as 12C + p (m1 > m2, theta_max 4.8 deg), segment 3 at 84.3 deg, mode 1:
#       ERROR, no chiSquared.out;
#   (d) control: the 7Li(p,n) test segment at 10 deg (inside theta_max) runs and gives finite numbers;
#   (e) a lab-differential point BELOW the (p,n) threshold (E_p = 1.88 MeV, threshold 1.885) at 10 deg: refused ("below
#       the threshold"); before, the Marion-Young map went NaN there with no message.
# dev 265c259: (a)-(c) and (e) exit 0 with nan in the output files ((c) and (e) with no message at all).
#
#   ./tests/kinematic_limits/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LI="$HERE/../7Li_p_ay"
N13="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/kinematic_limits.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
prep() { mkdir -p "$WORK/$1/output" "$WORK/$1/checks"; cp -r "$2/data" "$WORK/$1/"; }
run() { (cd "$WORK/$1" && printf '%s\n\n\n7\n' "$2" | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1); }
errors() { grep -c "beyond the kinematic maximum" "$WORK/$1/run.log"; }
below() { grep -c "below the threshold of the exit pair" "$WORK/$1/run.log"; }
has_nan() { grep -il "nan" "$WORK/$1/output"/AZUREOut_* "$WORK/$1/output/chiSquared.out" 2> /dev/null | head -1; }

# (a) test segment past theta_max
prep a "$LI"
awk '{print} /<segmentsTest>/{print "1 5 2 1.89 1.90 0.01 60 60 0 1"}' "$LI/7Li_p_ay.azr" > "$WORK/a/m.azr"
run a 3
if [ "$(errors a)" -ge 1 ] && [ ! -e "$WORK/a/output/AZUREOut_aa=5_R=2.extrap" ]; then
  ok "(a) 7Li(p,n) test segment at 60 deg lab (theta_max 21-45 deg): refused, no extrapolation file"
else bad "(a) test segment past theta_max: errors $(errors a), extrap file $(ls "$WORK/a/output" 2> /dev/null | tr '\n' ' '), nan in $(has_nan a)"; fi

# (b) data segment past theta_max
prep b "$LI"
printf '1.89 60 1.0 0.1\n1.90 60 1.0 0.1\n' > "$WORK/b/data/pn60.dat"
awk '{print} /<segmentsData>/{print "1 5 2 1.8 2.0 60 60 1 1 0 0 0 0 0 data/pn60.dat 0 0"}' "$LI/7Li_p_ay.azr" > "$WORK/b/m.azr"
run b 1
if [ "$(errors b)" -ge 1 ] && [ ! -s "$WORK/b/output/chiSquared.out" ]; then
  ok "(b) the same as a data segment: refused, no chiSquared.out"
else bad "(b) data segment past theta_max: errors $(errors b); chiSquared.out: $(grep '^Total' "$WORK/b/output/chiSquared.out" 2> /dev/null)"; fi

# (c) elastic with m1 > m2: 13N's pair 1 as 12C + p (swap j, pi, m and z of the two particles; the A_y segments off)
prep c "$N13"
awk 'BEGIN{inl=0;ins=0} /<levels>/{inl=1} /<\/levels>/{inl=0} /<segmentsData>/{ins=1} /<\/segmentsData>/{ins=0}
     { if (ins && NF > 14 && $8 == 7) $1 = 0;
       if (inl && NF > 25 && $6 == 1) { t=$13; $13=$15; $15=t; t=$14; $14=$16; $16=t; t=$18; $18=$19; $19=t; t=$20; $20=$21; $21=t }
       print }' "$N13/13N.azr" > "$WORK/c/m.azr"
run c 1
if [ "$(errors c)" -ge 1 ] && [ ! -s "$WORK/c/output/chiSquared.out" ]; then
  ok "(c) elastic 12C + p (m1 > m2, theta_max 4.8 deg) at 84.3 deg lab: refused, no chiSquared.out"
else bad "(c) elastic m1 > m2 past theta_max: errors $(errors c); chiSquared.out: $(grep '^Total' "$WORK/c/output/chiSquared.out" 2> /dev/null)"; fi

# (d) control inside theta_max
prep d "$LI"
awk '{print} /<segmentsTest>/{print "1 5 2 1.89 1.91 0.01 10 10 0 1"}' "$LI/7Li_p_ay.azr" > "$WORK/d/m.azr"
run d 3
f="$WORK/d/output/AZUREOut_aa=5_R=2.extrap"
if [ -s "$f" ] && [ "$(errors d)" -eq 0 ] && ! grep -qi nan "$f"; then
  ok "(d) control: the same segment at 10 deg runs ($(grep -c . "$f") energies, finite)"
else bad "(d) control at 10 deg: errors $(errors d), file $(ls -la "$f" 2> /dev/null | awk '{print $5}') bytes, nan: $(grep -ci nan "$f" 2> /dev/null)"; fi

# (e) below the exit threshold
prep e "$LI"
awk '{print} /<segmentsTest>/{print "1 5 2 1.88 1.88 0 10 10 0 1"}' "$LI/7Li_p_ay.azr" > "$WORK/e/m.azr"
run e 3
if [ "$(below e)" -ge 1 ] && [ ! -e "$WORK/e/output/AZUREOut_aa=5_R=2.extrap" ]; then
  ok "(e) a lab-differential point below the (p,n) threshold (E_p = 1.88 MeV): refused, no extrapolation file"
else bad "(e) below threshold: messages $(below e), nan in $(has_nan e)"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
