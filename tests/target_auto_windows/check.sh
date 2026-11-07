#!/usr/bin/env bash
#
# Automatic target-integration windows from the level scheme ("auto" in the windows token, 2026-11-06).
#
# A thin-target integration restricted to energy windows evaluates a point outside them once, at the depth-averaged
# energy E_bar of its target (tests/target_window_fallback).  With "auto" the engine places the windows itself, per
# segment and from that segment's entrance pair: around every resonance whose total width Gamma is small enough for
# the fallback to miss its peak by more than the tolerance (Delta^2/(3 Gamma^2) > tol, Delta the target thickness
# plus twice the beam spread), and around every particle threshold of the same compound nucleus (a target that
# straddles a threshold sees the cusp), reaching margin x max(Gamma, Delta) beyond the points whose target holds the
# resonance or threshold.  "auto:<tol>:<margin>" sets them (defaults 0.01 and 3); explicit ranges in the same token
# are added; a Note per segment lists the windows.  Older engines read no range from "auto" and integrate everywhere.
#
# Model: tests/thick_target_resonance (p+14N, the 278 keV lab 1/2+ resonance), 57 points 0.240-0.520 MeV lab on one
# capture pair, density 1e18, stopping a0/x^0.8 (the model of tests/target_window_fallback).
#   (a) "auto": the Note reports tolerance 0.01, margin 5 (max(3, 1/(2 sqrt(tol))), 2026-11-07), one narrow level, no
#       threshold, and some points left out;
#   (b) the output is byte-identical to the same line with the windows the Note lists written out (the windows are
#       applied as reported); the points inside are the full integration to every digit;
#   (c) every point within the tolerance (1e-2) of the full integration;
#   (d) "auto:0.01:1": the Note reports them, and the narrower margin integrates fewer points (the tolerance alone
#       is tested in (j): a tighter one also integrates the sub-barrier points);
#   (e) a particle threshold (a pair whose heavy partner is excited by 0.42 MeV, zero ANC, so the cross section does
#       not change): one threshold and a window over it (margin 3, so that some points stay outside);
#   (f) "auto,0.48-0.50": the explicit range is added to the automatic ones;
#   (g) data far from every level (0.40-0.52 MeV only): no window, every point at E_bar (equal to a line whose only
#       window lies outside the data);
#   (h) with a fitted density and a 10 keV blend the gradient check stays on the analytic block, no disagreement;
#   (i) refused with an ERROR and no output: "auto" on a thick target and on a line without a target integration;
#   (j) the sub-barrier curvature criterion (2026-11-06): data 0.40-0.52 MeV far from the level, a 41 keV target,
#       "auto:0.001:1": the Note reports sub-barrier points, and every point is within 1e-3 of the full integration
#       (without the criterion no window is placed and the E_bar fallback is ~0.4 % low).
#   (k) the default margin follows the tolerance: "auto" (margin 5) integrates more points than "auto:0.01:3" and is at
#       least as close to the full integration.
# Builds before the "auto" token (a20918d) ignore it and integrate every point: they fail 10 checks of (a)-(i), all
# but (c); builds with "auto" but without the curvature criterion and the tolerance-tied margin fail (a) and (j).
#
#   ./tests/target_auto_windows/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../thick_target_resonance"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
PY="$(command -v python3 || command -v python)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/target_auto_windows.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

A0=7.95e-21
A1=0.8
LINE="1  \"1\"  200  0  0  1  1.0e18  \"a0/x^a1\"  2  $A0  $A1  0  0  0  \"\"  0  0  0  20  50"
OUTNAME="AZUREOut_aa=1_R=3.out"

# model NAME TARGETINT-LINE [E_LO E_HI [THRESHOLD]]: the resonance project, one capture pair (exit 3, type 0), points
# every 5 keV from E_LO to E_HI (lab MeV, default 0.240-0.520); THRESHOLD=1 adds a closed channel of a pair whose
# heavy partner is excited by 0.42 MeV (key 9, zero ANC on the ground state: a threshold, no change in sigma)
model() {
  local d="$WORK/$1" lo="${3:-0.240}" hi="${4:-0.520}" thr="${5:-0}"
  mkdir -p "$d/output" "$d/checks" "$d/data"
  awk -v lo="$lo" -v hi="$hi" 'BEGIN { for (e = lo; e <= hi + 1e-9; e += 0.005) printf "%.3f 0 1.0e-06 1.0e-07\n", e }' > "$d/data/thin.dat"
  awk -v tline="$2" -v thr="$thr" \
      -v s1="1 1 3 0 6 0 180 0 1 0 0 0 0 0 data/thin.dat 0 0" '
    { sub(/\r$/, "") }
    /^<targetInt>$/ { print; print tline; skip = 1; next }
    /^<segmentsData>$/ { print; print s1; skip = 1; next }
    /^<segmentsTest>$/ { print; skip = 1; next }
    /^<\/targetInt>$|^<\/segmentsData>$|^<\/segmentsTest>$/ { skip = 0 }
    /^<levels>$/ { inl = 1; print; next }
    /^<\/levels>$/ { inl = 0 }
    inl && thr == 1 && !done && NF >= 31 && $6 == 1 {
      print; $6 = 9; $12 = 0; $17 = 0.42; print; done = 1; next }
    !skip { print }' "$BASE/thick_target_resonance.azr" > "$d/t.azr"
}
run() {  # NAME MODE FLAGS [ENV...]
  local d="$WORK/$1" mode="$2" flags="$3"; shift 3
  (cd "$d" && printf "$mode\n\n\n7\n" | env "$@" ${TIMEOUT_CMD:+"$TIMEOUT_CMD" 900} "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}
note() { grep -m1 "automatic target-integration windows" "$WORK/$1/run.log"; }
# windows NAME: the window list of the Note as a ranges token body ("lo-hi,lo-hi")
windows() { note "$1" | sed -n 's/.*windows (lab MeV) //p' | tr -d ' '; }
inside() { note "$1" | sed -n 's/.* \([0-9]*\) of [0-9]* points integrated.*/\1/p'; }
y() { awk 'NF >= 5 { sub(/\r$/, ""); print $4 }' "$WORK/$1/output/$OUTNAME" 2> /dev/null; }
same() { [ -s "$WORK/$1/output/$OUTNAME" ] && "$PY" -I -c 'import sys; sys.exit(open(sys.argv[1], "rb").read() != open(sys.argv[2], "rb").read())' \
           "$WORK/$1/output/$OUTNAME" "$WORK/$2/output/$OUTNAME"; }

model full "$LINE";                    run full 1 ""
model auto "$LINE  \"auto\"  0";       run auto 1 ""

# (a) the Note
n="$(note auto)"; nin="$(inside auto)"
if echo "$n" | grep -q "tolerance 0.01, margin 5): 1 narrow level(s), 0 threshold(s), [0-9]* of 57 points" && [ "${nin:-57}" -lt 57 ] && [ "${nin:-0}" -gt 5 ]; then
  ok "(a) \"auto\": one narrow level, no threshold, $nin of 57 points integrated, windows $(windows auto)"
else bad "(a) \"auto\": Note '${n:-none}'"; tail -3 "$WORK/auto/run.log"; fi

# (b) the same as the listed windows written out; inside = full to every digit
W="$(windows auto)"
model listed "$LINE  \"${W:-none}\"  0"; run listed 1 ""
if [ -n "$W" ] && same auto listed; then ok "(b) identical to the line with the listed windows \"$W\" written out"
else bad "(b) differs from the line with the listed windows \"${W:-none}\""; fi
read -r nsame nall <<< "$(paste <(y full) <(y auto) | awk '$1 == $2 { s++ } END { print s + 0, NR }')"
if [ "$nall" -eq 57 ] && [ "$nsame" -ge "${nin:-99}" ] && [ "$nsame" -lt 57 ]; then ok "(b) $nsame of 57 points equal the full integration to every digit (the $nin inside the windows)"
else bad "(b) $nsame of $nall points equal the full integration (inside: ${nin:-?})"; fi

# (c) within the tolerance everywhere
read -r worst at <<< "$(paste <(y full) <(y auto) <(awk '{print $1}' "$WORK/auto/data/thin.dat") | awk '{ r = $2 / $1 - 1; if (r < 0) r = -r; if (r > w) { w = r; e = $3 } } END { printf "%.3e %s", w, e }')"
if awk -v w="$worst" 'BEGIN { exit !(w <= 1e-2) }'; then ok "(c) every point within $worst of the full integration (tolerance 1e-2; largest at $at MeV)"
else bad "(c) a point is $worst off the full integration at $at MeV (tolerance 1e-2)"; fi

# (d) tolerance and margin
model tight "$LINE  \"auto:0.01:1\"  0"; run tight 1 ""
nt="$(inside tight)"
if note tight | grep -q "tolerance 0.01, margin 1)" && [ -n "$nt" ] && [ -n "$nin" ] && [ "$nt" -lt "$nin" ]; then
  ok "(d) \"auto:0.01:1\": reported, $nt points integrated (default: $nin)"
else bad "(d) \"auto:0.01:1\": Note '$(note tight)'"; fi

# (e) a threshold of the same compound nucleus
model thr "$LINE  \"auto:0.01:3\"  0" 0.240 0.520 1; run thr 1 ""   # margin 3: with 5 the windows cover every point
F="$(awk 'BEGIN { printf "%.10f", 14.0031 / (1.00783 + 14.0031) }')"
eth="$(awk -v f="$F" 'BEGIN { printf "%.4f", 0.42 / f }')"
covered="$(windows thr | tr ',' '\n' | awk -F- -v e="$eth" '$1 < e && $2 > e { c = 1 } END { print c + 0 }')"
if note thr | grep -q "1 narrow level(s), 1 threshold(s)" && [ "$covered" -eq 1 ] && same thr full 2> /dev/null; then
  bad "(e) threshold: windows found but the output equals the full integration (nothing left out)"
elif note thr | grep -q "1 narrow level(s), 1 threshold(s)" && [ "$covered" -eq 1 ]; then
  ok "(e) the threshold at $eth MeV lab gets a window ($(windows thr))"
else bad "(e) threshold at $eth MeV: Note '$(note thr)'"; fi

# (f) explicit ranges are added
model union "$LINE  \"auto,0.48-0.50\"  0"; run union 1 ""
if windows union | grep -q "0.48-0.5" && [ "$(inside union)" -gt "${nin:-999}" ]; then ok "(f) \"auto,0.48-0.50\": $(windows union)"
else bad "(f) \"auto,0.48-0.50\": Note '$(note union)'"; fi

# (g) no level near the data: every point at E_bar
model far "$LINE  \"auto\"  0" 0.400 0.520; run far 1 ""
model faroff "$LINE  \"0.10-0.11\"  0" 0.400 0.520; run faroff 1 ""
if note far | grep -q "0 narrow level(s), 0 threshold(s), 0 of 25 points" && same far faroff; then ok "(g) data far from every level: no window, every point at E_bar"
else bad "(g) data far from every level: Note '$(note far)'"; fi

# (h) gradient: fitted density and a blend, analytic block
model grad "$LINE  \"auto\"  0.01  densityfit 1 0"; run grad 2 "--use-gradient" AZURE_GRAD_CHECK=exit
line="$(grep "components with relative disagreement" "$WORK/grad/run.log" | head -1)"
if echo "$line" | grep -q "analytic block used;.*density 0/1;.*Jacobian path: density 0/1" && note grad > /dev/null; then
  ok "(h) gradient check on the analytic path: $line"
else bad "(h) gradient check: ${line:-no [grad-check] line}"; fi

# (i) refusals
model thick "1  \"1\"  200  0  0  1  0  \"a0/x^a1\"  2  $A0  $A1  0  0  0  \"\"  0  0  0  20  50  \"auto\"  0"; run thick 1 ""
if grep -q "ERROR.*auto.*THIN-target" "$WORK/thick/run.log" && [ ! -s "$WORK/thick/output/$OUTNAME" ]; then ok "(i) \"auto\" on a thick target: refused, no output"
else bad "(i) \"auto\" on a thick target was not refused"; tail -2 "$WORK/thick/run.log"; fi
model conv "1  \"1\"  50  1  0.002  0  0  \"\"  0  0  0  0  \"\"  0  0  0  20  50  \"auto\"  0"; run conv 1 ""
if grep -q "ERROR.*auto.*THIN-target" "$WORK/conv/run.log" && [ ! -s "$WORK/conv/output/$OUTNAME" ]; then ok "(i) \"auto\" without a target integration: refused, no output"
else bad "(i) \"auto\" on a convolution-only line was not refused"; tail -2 "$WORK/conv/run.log"; fi

# (j) sub-barrier curvature: thicker target, data away from the level, tolerance 1e-3, margin 1
LINED="1  \"1\"  200  0  0  1  2.5e18  \"a0/x^a1\"  2  $A0  $A1  0  0  0  \"\"  0  0  0  20  50"
model steep "$LINED  \"auto:0.001:1\"  0" 0.400 0.520; run steep 1 ""
model steepfull "$LINED" 0.400 0.520; run steepfull 1 ""
read -r worst at <<< "$(paste <(y steepfull) <(y steep) <(awk '{print $1}' "$WORK/steep/data/thin.dat") | awk '{ r = $2 / $1 - 1; if (r < 0) r = -r; if (r > w) { w = r; e = $3 } } END { printf "%.3e %s", w, e }')"
if note steep | grep -q "sub-barrier curvature: [0-9]* point(s)" && awk -v w="$worst" 'BEGIN { exit !(w <= 1e-3) }'; then
  ok "(j) sub-barrier curvature: $(note steep | sed -n 's/.*\(sub-barrier curvature: [^;]*\).*/\1/p'); every point within $worst of the full integration"
else bad "(j) sub-barrier curvature: Note '$(note steep)'; worst $worst at $at MeV (tolerance 1e-3)"; fi

# (k) the default margin follows the tolerance
model m3 "$LINE  \"auto:0.01:3\"  0"; run m3 1 ""
n3="$(inside m3)"
w3="$(paste <(y full) <(y m3) | awk '{ r = $2 / $1 - 1; if (r < 0) r = -r; if (r > w) w = r } END { printf "%.3e", w }')"
wd="$(paste <(y full) <(y auto) | awk '{ r = $2 / $1 - 1; if (r < 0) r = -r; if (r > w) w = r } END { printf "%.3e", w }')"
if [ -n "$n3" ] && [ -n "$nin" ] && [ "$nin" -gt "$n3" ] && awk -v a="$wd" -v b="$w3" 'BEGIN { exit !(a <= b) }'; then
  ok "(k) default margin 5: $nin points integrated, worst $wd; margin 3: $n3 points, worst $w3"
else bad "(k) default margin: $nin points (worst $wd) vs margin 3: ${n3:-?} points (worst $w3)"; fi

if [ "$failures" -ne 0 ]; then echo "target_auto_windows: $failures check(s) FAILED"; exit 1; fi
echo "target_auto_windows: all checks passed"
