#!/usr/bin/env bash
#
# Secondary-gamma multipole mixing at every Legendre order, and a warning for channels beyond the Legendre cap
# (audit 2026-10-25, checklist t3-24; 2026-10-27).
#
# (1) A UPOS segment with a mixed L/L+1 gamma decay uses the radiation parameter
#       R_k = (R_k(L,L) + 2 delta R_k(L,L') + delta^2 R_k(L',L')) / (1 + delta^2)        (Brune & deBoer 2020, Eq. 7)
#     at each even order k.  Up to dev 7d4558b the mixed form was applied only "if (R_LpLp != 0)", so at an order where
#     R_k(L',L') happens to vanish while R_k(L,L) does not -- a 3/2 -> 3/2 M1/E2 decay at k = 2: R_2(2,2) = 0,
#     R_2(1,1) = -0.400, R_2(1,2) = +0.775 -- the 2 delta R_k(L,L') interference and the 1/(1 + delta^2) were dropped
#     for that k, leaving the a_2 coefficient at its pure-M1 value whatever delta.
#     Model: tests/secondary_gamma/models/mixing.azr (5/2+ -> alpha + 13C(3.68, 3/2-), 13 angles at one energy) with the
#     UPOS line's final state set to J = 3/2 and delta = 0.3 (and 0): the distribution is a0 + a2 P2(cos theta) only
#     (k = 4 needs (3/2 3/2 4), impossible), and a2/a0 at delta = 0.3 over a2/a0 at delta = 0 must be
#       F = (R_2(1,1) + 0.6 R_2(1,2)) / (1.09 R_2(1,1)) = -0.148528.
#       (a) both distributions are a0 + a2 P2 at every angle (1e-8 of a0);  (b) the ratio is F to 1e-5.
#     dev 7d4558b: the ratio is exactly 1 (a2 does not change with delta).
# (2) Legendre orders above Config::maxLOrder = 20 are never computed, so a particle channel with l > 10 silently
#     truncated every differential cross section it entered.  Now the run log warns per such channel:
#       (c) tests/13N with an added J = 21/2- level whose p + 12C channel has l = 11: the log carries the WARNING naming
#           l = 11, and the run still completes;  (d) the unmodified 13N prints no such warning.
#     dev 7d4558b: no warning.
#
#   ./tests/mixing_guard_lorder/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SG="$HERE/../secondary_gamma/models"
BASE13="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/mixing_guard_lorder.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }

UPOS='1  2  6  5.569971  6.569971  0  180  4  1  0  0  0  0.02  1  data/mixing_upos.dat  0  1  1  0.5  -0.094'
grep -qF "$UPOS" "$SG/mixing.azr" || { echo "  FAIL  tests/secondary_gamma/models/mixing.azr no longer holds the UPOS line '$UPOS'"; exit 1; }
# mixing NAME DELTA: the mixing model with the UPOS final state 3/2 and the given delta, the other two segments dropped
mixing() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$SG/data" "$WORK/$1/"
  awk -v upos="$UPOS" -v delta="$2" '
    { sub(/\r$/, "") }
    /^<segmentsData>$/ { print; print "1  2  6  5.569971  6.569971  0  180  4  1  0  0  0  0.02  1  data/mixing_upos.dat  0  1  1  1.5  " delta; skip = 1; next }
    /^<\/segmentsData>$/ { skip = 0 }
    !skip { print }' "$SG/mixing.azr" > "$WORK/$1/t.azr"
}
mixing d03 0.3
mixing d00 0
run d03; run d00
# W(theta) rows: angle, fit
dist() { awk '{ sub(/\r$/, "") } NF >= 5 { printf "%.4f %.12e\n", $3, $4 }' "$WORK/$1/output/AZUREOut_aa=2_R=6.out" 2> /dev/null; }
dist d03 > "$WORK/w03"; dist d00 > "$WORK/w00"
n="$(wc -l < "$WORK/w03" | tr -d ' ')"
if [ "$n" -eq 13 ]; then ok "the mixed-decay runs wrote the 13 angles"; else bad "d03 wrote $n angles, expected 13"; tail -3 "$WORK/d03/run.log"; fi
# a2/a0 from W(0) and W(90): a0 = (W0 + 2 W90)/3, a2 = 2 (W0 - W90)/3; then the largest residual of a0 + a2 P2 over all angles
coef() { awk 'BEGIN { pi = atan2(0, -1) } { W[$1 + 0] = $2; A[NR] = $1 + 0 } END { w0 = W[0]; w90 = W[90]; a0 = (w0 + 2 * w90) / 3; a2 = 2 * (w0 - w90) / 3; m = 0;
  for (i = 1; i <= NR; i++) { c = cos(A[i] * pi / 180); p2 = 0.5 * (3 * c * c - 1); r = W[A[i]] - a0 - a2 * p2; if (r < 0) r = -r; if (r > m) m = r }
  printf "%.12e %.12e %.3e\n", a0, a2, m / a0 }' "$1"; }
read a0_03 a2_03 res_03 < <(coef "$WORK/w03")
read a0_00 a2_00 res_00 < <(coef "$WORK/w00")
if awk -v a="$res_03" -v b="$res_00" 'BEGIN { exit !(a <= 1e-8 && b <= 1e-8) }'; then ok "(a) both distributions are a0 + a2 P2 at every angle (residuals $res_03, $res_00 of a0)"; else bad "(a) the distributions are not a0 + a2 P2 (residuals $res_03, $res_00 of a0)"; fi
F="$(awk 'BEGIN { r11 = -0.4; r12 = 0.774596669; d = 0.3; printf "%.6f", (r11 + 2 * d * r12) / ((1 + d * d) * r11) }')"
ratio="$(awk -v a="$a2_03" -v b="$a0_03" -v c="$a2_00" -v d="$a0_00" 'BEGIN { printf "%.6f", (a / b) / (c / d) }')"
if awk -v r="$ratio" -v f="$F" 'BEGIN { d = r - f; if (d < 0) d = -d; exit !(d <= 1e-5) }'; then ok "(b) a2/a0 at delta = 0.3 over delta = 0 is $ratio = F = $F (3/2 -> 3/2 M1/E2, k = 2)"; else bad "(b) a2/a0 at delta = 0.3 over delta = 0 is $ratio, expected F = $F (the guard dropped the mixing at k = 2)"; fi

# ---- (c), (d) the Legendre cap warning -----------------------------------------------------------------------------
LINE1="$(awk '/<levels>/ { f = 1; next } /<\/levels>/ { f = 0 } f && NF > 26 { print; exit }' "$BASE13/13N.azr")"
[ -n "$LINE1" ] || { bad "tests/13N/13N.azr has no level line"; echo "$failures check(s) FAILED"; exit 1; }
# a J = 21/2- level at Ex 8 MeV with the p + 12C channel at l = 11 (the file stores 2s and 2l): J, pi, E, fix, aa, ir, 2s, 2l,
# levelID, active, chFix, width (eV), then the pair columns copied from the first line
HIGHL="$(echo "$LINE1" | awk '{ $1 = "10.5"; $2 = "-1"; $3 = "8.000000"; $4 = "0"; $7 = "1"; $8 = "22"; $9 = "13"; $10 = "1"; $11 = "0"; $12 = "50000"; print }')"
for m in highl plain; do mkdir -p "$WORK/$m/output" "$WORK/$m/checks"; cp -r "$BASE13/data" "$WORK/$m/"; done
awk -v extra="$HIGHL" '{ sub(/\r$/, "") } /^<\/levels>$/ { print extra } { print }' "$BASE13/13N.azr" > "$WORK/highl/t.azr"
sed 's/\r$//' "$BASE13/13N.azr" > "$WORK/plain/t.azr"
run highl; run plain
if grep -q "WARNING: channel l = 11 (J = 10.5, pair 1) exceeds 10" "$WORK/highl/run.log" && [ -f "$WORK/highl/output/chiSquared.out" ]; then ok "(c) a p + 12C channel with l = 11 is flagged: $(grep -o 'WARNING: channel l = 11[^:]*' "$WORK/highl/run.log" | head -1), and the run completes"; else bad "(c) no warning for the l = 11 channel (or the run failed)"; grep -i "warning\|error" "$WORK/highl/run.log" | head -3; fi
if ! grep -q "WARNING: channel l =" "$WORK/plain/run.log" && [ -f "$WORK/plain/output/chiSquared.out" ]; then ok "(d) the unmodified 13N prints no such warning"; else bad "(d) the unmodified 13N printed a channel-l warning or failed"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
