#!/usr/bin/env bash
#
# A thin-target integration restricted to energy windows evaluates the points OUTSIDE the windows at their
# depth-averaged energy E_bar, not at the data (beam) energy (2026-11-05).
#
# The data of a target integration carry the beam energy E0.  Since aeb54b7 (2026-08-31) a <targetInt> line can be
# restricted to lab-energy windows ("lo-hi,..." blend [tol]); a point outside every window was left a plain point at
# E0, so its model was sigma(E0) instead of the target average -- off by the centroid shift ~Delta/2.  On real data
# (13C+a ND2021r, 16 lines, windows around the narrow levels) that made the windowed model WORSE than no target
# integration at all: chi2 421,432 against 400,616 integrated everywhere.  Now such a point is sigma(E_bar),
#   E_bar = int E dE/eps / int dE/eps over [E0 - Delta, E0]   (TargetEffect::MeanEnergyLoss),
# exact for a cross section linear across the target; window edges blend with sigma(E_bar) and the automatic
# tolerance falls back to it too.
#
# Model: tests/thick_target_resonance (p + 14N, external capture off), total capture.  Segment 1: 57 lab energies
# 240-520 keV under one thin-target line, epsilon = a0 / x^0.8 (a0 = 7.95e-21 MeV cm^2), N = 1e18 atoms/cm^2
# (Delta = 18 keV at 360 keV), 200 sub-points; the 259.7 keV c.m. (278 keV lab) resonance is inside the window
# 0.27-0.32 MeV.  Segment 2: the bare cross section on a 0.1 keV lab grid, 200-530 keV (the truth table).  For the power law,
# E_exit = (E0^1.8 - 1.8 a0 N)^(1/1.8) and E_bar = (1.8/2.8) (E0^2.8 - E_exit^2.8) / (E0^1.8 - E_exit^1.8).
#   (a) windows, hard edges: every point outside the window equals sigma(E_bar) of the truth table to 2e-5;
#   (b) at the points whose target stays 30 keV clear of the resonance, the full integration and sigma(E_bar) agree
#       within 1.3 x the curvature estimate |sigma''(E_bar)| Delta^2 / 24 (sigma'' from the truth table; up to 11 %
#       at 240 keV, where the Coulomb barrier bends sigma by e every 16 keV across a 25 keV target), and sigma(E0), the
#       old value, is at least 3 x further off at all but two of them (the flat minimum of sigma near 420 keV);
#   (a') the same as (a) on a plain segment (one capture pair) instead of the total-capture composite;
#   (c) the points inside the window are identical to the full integration (all digits);
#   (d) blend width 10 keV (a ramp inside each window edge): every point lies between sigma(E_bar) and the full
#       integration (2e-5 slack, the table's interpolation), the points outside the window are sigma(E_bar) as in (a);
#   (e) automatic tolerance 1e-3 without windows, on the plain segment (the total-capture composite never applies
#       the automatic decision, before or after): every point within 1e-3 of the full integration, some points let
#       off (not identical to it), the resonance points integrated (identical);
#   (f) windows with a fitted density (three fallback nodes, at E_bar of 0.01 N0, N0 and the layout density): at
#       the nominal density within 1e-4 of (a); with a 10 keV blend the gradient check (AZURE_GRAD_CHECK) stays on the
#       analytic block (the blend is now linear inside the kernel) and finds no disagreement, density included;
#   (g) refused with an ERROR and no output: windows on a THICK target (density 0); windows on a line with doppler;
#   (h) the line without windows is unchanged: with REF_AZURE2=<older build> its output is byte-identical.
# Builds before 2026-11-05 fail 7 checks: (a), (a') and (d) (the outside points are sigma(E0), up to 3.7x off),
# (e) (1 point let off: the old estimate charges the slope as well), (f) (blended points bail out to finite
# differences) and both refusals of (g).  (b) and (c) pass on every build: they test the model and the inside points.
#
#   ./tests/target_window_fallback/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../thick_target_resonance"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/target_window_fallback.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
num='^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$'
le() { awk -v a="$1" -v b="$2" -v re="$num" 'BEGIN { exit !(a ~ re && a + 0 <= b + 0) }'; }

A0=7.95e-21
A1=0.8
DENS=1.0e18
WIN="0.27-0.32"
F="$(awk 'BEGIN { printf "%.10f", 14.0031 / (1.00783 + 14.0031) }')"   # c.m. / lab
LINE="1  \"1\"  200  0  0  1  $DENS  \"a0/x^a1\"  2  $A0  $A1  0  0  0  \"\"  0  0  0  20  50"

# model NAME TARGETINT-LINE [EXIT]: the resonance project with the thin-target segment and the truth table, both
# for exit key EXIT (default -1, the total capture of all 12 capture pairs, a composite)
model() {
  local d="$WORK/$1" exit="${3:--1}" typ=3
  [ "$exit" = "-1" ] || typ=0   # data type 3 is the total capture whatever the exit key; 0 = one pair, angle-integrated
  mkdir -p "$d/output" "$d/checks" "$d/data"
  awk 'BEGIN { for (i = 0; i <= 56; i++) printf "%.3f 0 1.0e-06 1.0e-07\n", 0.240 + i * 0.005 }' > "$d/data/thin.dat"
  awk 'BEGIN { for (i = 0; i <= 3300; i++) printf "%.4f 0 1.0e-06 1.0e-07\n", 0.20 + i * 0.0001 }' > "$d/data/fine.dat"
  awk -v tline="$2" \
      -v s1="1 1 $exit 0 6 0 180 $typ 1 0 0 0 0 0 data/thin.dat 0 0" \
      -v s2="1 1 $exit 0 6 0 180 $typ 1 0 0 0 0 0 data/fine.dat 0 0" '
    { sub(/\r$/, "") }
    /^<targetInt>$/ { print; print tline; skip = 1; next }
    /^<segmentsData>$/ { print; print s1; print s2; skip = 1; next }
    /^<segmentsTest>$/ { print; skip = 1; next }
    /^<\/targetInt>$|^<\/segmentsData>$|^<\/segmentsTest>$/ { skip = 0 }
    !skip { print }' "$BASE/thick_target_resonance.azr" > "$d/t.azr"
}
# run NAME MODE FLAGS [ENV...]
run() {
  local d="$WORK/$1" mode="$2" flags="$3"; shift 3
  (cd "$d" && printf "$mode\n\n\n7\n" | env "$@" ${TIMEOUT_CMD:+"$TIMEOUT_CMD" 900} "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}
OUTNAME="AZUREOut_aa=1_TOTAL_CAPTURE.out"
# block K of NAME's output -> "E_cm fit" rows (full precision as written); OUT=<file name> for another exit
block() { awk -v k="$2" 'BEGIN { RS = ""; FS = "\n" } NR == k { print }' "$WORK/$1/output/${OUT:-$OUTNAME}" 2> /dev/null |
          awk '{ sub(/\r$/, "") } NF >= 5 { print $1, $4 }'; }

model full "$LINE";                      run full 1 ""
model win "$LINE  \"$WIN\"  0";          run win 1 ""
model blend "$LINE  \"$WIN\"  0.01";     run blend 1 ""
for m in full win blend; do block $m 1 > "$WORK/$m.y"; done
block full 2 > "$WORK/sigma"
# the plain segment (exit 3, one capture pair): full, windows, automatic tolerance
OUT3="AZUREOut_aa=1_R=3.out"
model full3 "$LINE" 3;                   run full3 1 ""
model win3 "$LINE  \"$WIN\"  0" 3;       run win3 1 ""
model auto "$LINE  \"\"  0  0.001" 3;    run auto 1 ""
for m in full3 win3 auto; do OUT="$OUT3" block $m 1 > "$WORK/$m.y"; done
OUT="$OUT3" block full3 2 > "$WORK/sigma3"
nf="$(wc -l < "$WORK/full.y" | tr -d ' ')"; ns="$(wc -l < "$WORK/sigma" | tr -d ' ')"
if [ "$nf" -eq 57 ] && [ "$ns" -eq 3301 ]; then ok "the full run wrote the 57 yields and the bare cross section on 3301 energies"
else bad "the full run wrote $nf yields and $ns cross sections (expected 57 and 3301)"; tail -3 "$WORK/full/run.log"; fi

# table: E0_lab  E_bar_lab  sigma(E_bar)  sigma(E0)  Y_full  Y_win  Y_blend  Y_auto  inWindow  far  curvature
#   far: the target [E_exit, E0] stays 30 keV clear of the 278 keV (lab) resonance;
#   curvature: |sigma''(E_bar)| Delta^2 / 24, sigma'' a second difference of the truth table with step Delta/2
# table SIGMA-TABLE FULL WIN BLEND AUTO (a missing column: give the same file again)
table() {
paste "$2" "$3" "$4" "$5" | awk -v f="$F" -v a0="$A0" -v a1="$A1" -v n="$DENS" '
  FNR == NR { E[NR] = $1; S[NR] = $2; ns = NR; next }
  function sig(ecm,   i, t) {
    for (i = 1; i < ns; i++) if (ecm >= E[i] && ecm <= E[i + 1]) { t = (ecm - E[i]) / (E[i + 1] - E[i]); return S[i] + t * (S[i + 1] - S[i]) }
    return -1 }
  { e0 = $1 / f; p = 1 + a1
    ex = (e0 ^ p - p * a0 * n) ^ (1 / p)
    eb = (p / (p + 1)) * (e0 ^ (p + 1) - ex ^ (p + 1)) / (e0 ^ p - ex ^ p)
    inw = (e0 >= 0.27 - 1e-9 && e0 <= 0.32 + 1e-9) ? 1 : 0
    far = (e0 - 0.2779 > 0.030 + (e0 - ex) || 0.2779 - e0 > 0.030) ? 1 : 0
    h = 0.5 * (e0 - ex); s0 = sig(f * eb); sp = sig(f * (eb + h)); sm = sig(f * (eb - h))
    curv = (sp > 0 && sm > 0) ? (sp - 2 * s0 + sm) / (h * h) : 0; if (curv < 0) curv = -curv
    printf "%.4f %.8f %.12e %.12e %s %s %s %s %d %d %.6e\n", e0, eb, s0, sig(f * e0), $2, $4, $6, $8, inw, far,
           curv * (e0 - ex) ^ 2 / 24 }' "$1" -
}
table "$WORK/sigma" "$WORK/full.y" "$WORK/win.y" "$WORK/blend.y" "$WORK/full.y" > "$WORK/table"
table "$WORK/sigma3" "$WORK/full3.y" "$WORK/win3.y" "$WORK/full3.y" "$WORK/auto.y" > "$WORK/table3"
[ -n "${KEEP_WORK:-}" ] && echo "  (table: $WORK/table)"

# (a) outside the window = sigma(E_bar)
read na worst <<< "$(awk '$9 == 0 { r = ($6 / $3 - 1); if (r < 0) r = -r; if (r > w) w = r; c++ } END { printf "%d %.3e", c, w }' "$WORK/table")"
if [ "$na" -ge 40 ] && le "$worst" 2e-5; then ok "(a) windows: the $na points outside are sigma(E_bar) of the truth table (max rel difference $worst)"
else bad "(a) windows: the $na points outside differ from sigma(E_bar) by up to $worst (limit 2e-5)"; fi

read na3 worst3 <<< "$(awk '$9 == 0 { r = ($6 / $3 - 1); if (r < 0) r = -r; if (r > w) w = r; c++ } END { printf "%d %.3e", c, w }' "$WORK/table3")"
if [ "$na3" -ge 40 ] && le "$worst3" 2e-5; then ok "(a') the same on a plain segment (one capture pair): $na3 points, max rel difference $worst3"
else bad "(a') plain segment: the $na3 points outside differ from sigma(E_bar) by up to $worst3 (limit 2e-5)"; fi

# (b) far from the resonance: full integration = sigma(E_bar) to 2e-3; sigma(E0) more than 2 % off
read nb nover wfb weak <<< "$(awk '$10 == 1 { d = $5 - $3; if (d < 0) d = -d; q = $5 - $4; if (q < 0) q = -q
                                        if (d > 1.3 * $11 + 1e-5 * $5) o++
                                        r = d / $5; if (r > w) w = r
                                        if (q < 3 * d) k++; c++ }
                               END { printf "%d %d %.3e %d", c, o + 0, w, k + 0 }' "$WORK/table")"
if [ "$nb" -ge 30 ] && [ "$nover" -eq 0 ]; then ok "(b) $nb points clear of the resonance: |full - sigma(E_bar)| within 1.3 x |sigma''| Delta^2/24 everywhere (largest $wfb relative)"
else bad "(b) $nover of $nb points clear of the resonance exceed 1.3 x the curvature estimate (largest $wfb relative)"; fi
if [ "$weak" -le 2 ]; then ok "(b) ... and sigma(E0), the old outside value, is at least 3 x further off at $((nb - weak)) of them"
else bad "(b) sigma(E0) is not 3 x further off at $weak points (allowed 2): the check has no teeth"; fi

# (c) inside the window: identical to the full integration
read nc ndiff <<< "$(awk '$9 == 1 { c++; if ($6 != $5) d++ } END { printf "%d %d", c, d + 0 }' "$WORK/table")"
if [ "$nc" -ge 9 ] && [ "$ndiff" -eq 0 ]; then ok "(c) the $nc points inside the window are the full integration (all digits)"
else bad "(c) $ndiff of the $nc points inside the window differ from the full integration"; fi

# (d) blend: between sigma(E_bar) and the full integration; sigma(E_bar) beyond the blend
read nd nout wfar <<< "$(awk '{ lo = ($3 < $5) ? $3 : $5; hi = ($3 < $5) ? $5 : $3
                                if ($7 < lo * (1 - 2e-5) || $7 > hi * (1 + 2e-5)) o++
                                if ($9 == 0) { r = ($7 / $3 - 1); if (r < 0) r = -r; if (r > w) w = r }
                                c++ } END { printf "%d %d %.3e", c, o + 0, w }' "$WORK/table")"
if [ "$nd" -eq 57 ] && [ "$nout" -eq 0 ] && le "$wfar" 2e-5; then ok "(d) blend 10 keV: all $nd points between sigma(E_bar) and the full integration; outside the window sigma(E_bar) to $wfar"
else bad "(d) blend: $nout points outside [sigma(E_bar), full], beyond the blend up to $wfar from sigma(E_bar)"; fi

# (e) automatic tolerance
read ne wauto nlet nres <<< "$(awk '{ r = ($8 / $5 - 1); if (r < 0) r = -r; if (r > w) w = r; if ($8 != $5) l++
                                      if ($1 >= 0.278 && $1 <= 0.300 && $8 == $5) k++; c++ }
                                 END { printf "%d %.3e %d %d", c, w, l + 0, k + 0 }' "$WORK/table3")"
if le "$wauto" 1e-3 && [ "$nlet" -ge 10 ] && [ "$nres" -ge 4 ]; then ok "(e) automatic tolerance 1e-3: within $wauto of the full integration, $nlet of $ne points let off, the $nres resonance points integrated"
else bad "(e) automatic tolerance: up to $wauto from the full integration (limit 1e-3), $nlet let off (want >= 10), $nres resonance points integrated (want >= 4)"; fi

# (f) fitted density with windows: three fallback nodes, and the analytic gradient
model fitwin "$LINE  \"$WIN\"  0  densityfit 1 0"
run fitwin 1 ""
block fitwin 1 > "$WORK/fitwin.y"
read nfit wfit <<< "$(paste "$WORK/win.y" "$WORK/fitwin.y" | awk '{ r = ($4 / $2 - 1); if (r < 0) r = -r; if (r > w) w = r; c++ } END { printf "%d %.3e", c, w }')"
if [ "$nfit" -eq 57 ] && le "$wfit" 1e-4; then ok "(f) fitted density, windows: at the nominal density within $wfit of the fixed-density run"
else bad "(f) fitted density, windows: $nfit points, up to $wfit from the fixed-density run (limit 1e-4)"; tail -3 "$WORK/fitwin/run.log"; fi
model grad "$LINE  \"$WIN\"  0.01  densityfit 1 0"   # with the blend: linear in the kernel now
run grad 2 "--use-gradient" AZURE_GRAD_CHECK=exit
line="$(grep "components with relative disagreement" "$WORK/grad/run.log" | head -1)"
if echo "$line" | grep -q "analytic block used;.*density 0/1;.*Jacobian path: density 0/1"; then
  ok "(f) gradient check on the analytic path: $line"
else bad "(f) gradient check: ${line:-no [grad-check] line}"; grep -m2 "unsupported\|grad-check" "$WORK/grad/run.log"; fi

# (g) refusals
model thick "1  \"1\"  200  0  0  1  0  \"a0/x^a1\"  2  $A0  $A1  0  0  0  \"\"  0  0  0  20  50  \"$WIN\"  0"
run thick 1 ""
if grep -q "ERROR.*THICK-target" "$WORK/thick/run.log" && [ ! -s "$WORK/thick/output/$OUTNAME" ]; then ok "(g) windows on a thick target: refused, no output"
else bad "(g) windows on a thick target were not refused"; tail -3 "$WORK/thick/run.log"; fi
model dop "$LINE  \"$WIN\"  0  doppler 300"
run dop 1 ""
if grep -q "ERROR.*doppler" "$WORK/dop/run.log" && [ ! -s "$WORK/dop/output/$OUTNAME" ]; then ok "(g) windows on a line with doppler: refused, no output"
else bad "(g) windows on a line with doppler were not refused"; tail -3 "$WORK/dop/run.log"; fi

# (h) optional: the line without windows is unchanged
if [ -n "${REF_AZURE2:-}" ]; then
  model ref "$LINE"
  (cd "$WORK/ref" && printf '1\n\n\n7\n' | run_limited 900 "$REF_AZURE2" --no-gui --no-readline t.azr > run.log 2>&1)
  if "${PY:-python3}" -I -c 'import sys; sys.exit(open(sys.argv[1], "rb").read() != open(sys.argv[2], "rb").read())' \
       "$WORK/ref/output/$OUTNAME" "$WORK/full/output/$OUTNAME"; then ok "(h) without windows: byte-identical to $REF_AZURE2"
  else bad "(h) without windows: output differs from $REF_AZURE2"; fi
fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
