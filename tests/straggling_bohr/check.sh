#!/usr/bin/env bash
#
# Energy straggling: the coefficient is a c.m. quantity, and a coefficient of 0 means Bohr straggling
# (audit 2026-10-25, tier 3 entry 1; 2026-10-27).
#
# The straggling width is sigma = c_s sqrt(dE_cm [keV]) keV with the energy loss in the c.m. and no frame factor: c_s is
# a CENTRE-OF-MASS coefficient (it was never said which frame it belonged to; every other number on the line is lab).
# Straggling on with c_s = 0 now means Bohr straggling, Omega_lab^2 = 4 pi Z1^2 e^4 Z2 N x with N x the areal density
# of target atoms traversed (int dE / epsilon over the sub-points), converted to the c.m. with m2/(m1+m2); the default
# coefficient was 0.04 keV^1/2 (several times below Bohr, essentially no straggling) and is 0 now.
#
# Models are made at run time from tests/thick_target_resonance (p + 14N, the 259.7 keV c.m. resonance, Gamma = 0.97
# keV) with one <targetInt> line for both segments: a 0.1 keV beam sigma, a CONSTANT stopping cross section
# epsilon = 1.8e-20 MeV cm^2 (18 eV / 1e15 atoms cm^-2) and 1.1e18 atoms/cm^2 (a 20 keV lab target), and eight yield
# points on the rising edge of the thick-target yield (lab 276 - 300 keV).  With a constant epsilon Bohr's width is
# exactly c_B sqrt(dE_cm) with c_B = sqrt(1000 f 4 pi Z1^2 Z2 e^4 / epsilon) keV^1/2, f = m2/(m1+m2): 0.3075 here.
#   (a) straggling off, (b) straggling on with coefficient 0 (Bohr), (c) coefficient c_B explicit:
#       (b) == (c) at every point to 1e-6 relative, and (b) differs from (a) by more than 1e-3 at some point;
#   (d) the 13N thick target (density 0, p + 12C down to 1 keV) with Bohr straggling: every yield finite and positive
#       (the areal density accumulates 1/epsilon over the sub-points, ignoring any node with epsilon <= 0).
# dev c6ed589: (b) has no straggling at all, equal to (a), and up to 90 % off (c) (lab 300 keV, the back of the target on the resonance).
#
#   ./tests/straggling_bohr/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../thick_target_resonance"
BASE13="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/straggling_bohr.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# Bohr's coefficient for the constant stopping cross section, from the same constants as Constants.h
EPS=1.8e-20
DENS=1.1e18
CB="$(awk -v eps="$EPS" 'BEGIN { pi = 3.141592650; hbarc = 197.32696310; alpha = 1.0 / 137.0359996790;
  e2 = alpha * hbarc * 1.0e-13; f = 14.0031 / (1.00783 + 14.0031); z1 = 1; z2 = 7;
  printf "%.6f", sqrt(1000.0 * f * 4.0 * pi * z1 * z1 * z2 * e2 * e2 / eps) }')"
echo "  Bohr coefficient for p + 14N, epsilon = $EPS MeV cm^2: c_B = $CB keV^1/2 (c.m.)"
if awk -v c="$CB" 'BEGIN { exit !(c > 0.2 && c < 0.6) }'; then ok "c_B is in the 0.2 - 0.6 keV^1/2 range expected for light ions"; else bad "c_B = $CB is outside 0.2 - 0.6"; fi

# make_model NAME STRAGGLING-FLAG COEFFICIENT: the resonance project with its <targetInt> block replaced by one line
# for segments 1 and 2 (fields: active, segments, sub-points, convolution, sigma, target integration, density,
# equation, n, parameters, Q flag, nQ, conv flag, equation, nconv, straggling, coefficient, widths, points per width)
make_model() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks" "$WORK/$1/data"
  printf '%s 0 5.0e-12 5.0e-13\n' 0.276 0.277 0.278 0.279 0.280 0.282 0.290 0.300 > "$WORK/$1/data/yield_scan.dat"
  awk -v line="1  \"1,2\"  20  1  0.0001  1  $DENS  \"a0\"  1  $EPS  0  0  0  \"\"  0  $2  $3  20  50" '
    { sub(/\r$/, "") }
    /^<targetInt>$/ { print; print line; skip = 1; next }
    /^<\/targetInt>$/ { skip = 0 }
    !skip { print }' "$BASE/thick_target_resonance.azr" > "$WORK/$1/t.azr"
}
run() { (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1); }
# yields NAME -> "E fit" rows of both segments (total-capture output file)
yields() { awk '{ sub(/\r$/, "") } NF >= 5 { print $1, $4 }' "$WORK/$1/output/AZUREOut_aa=1_TOTAL_CAPTURE.out" 2> /dev/null; }
# maxrel A-FILE B-FILE -> the largest |a/b - 1| over paired rows
maxrel() { paste "$1" "$2" | awk 'BEGIN { m = 0 } { if ($4 == 0) { if ($2 != 0) m = 1 } else { r = $2 / $4 - 1; if (r < 0) r = -r; if (r > m) m = r } } END { print m }'; }

make_model off 0 0
make_model bohr 1 0
make_model coef 1 "$CB"
for m in off bohr coef; do run "$m"; yields "$m" > "$WORK/$m.yields"; done
n="$(wc -l < "$WORK/bohr.yields" | tr -d ' ')"
if [ "$n" -eq 16 ]; then ok "the three runs wrote 16 yields each (two segments of eight points)"; else bad "bohr run wrote $n yields, expected 16"; tail -3 "$WORK/bohr/run.log"; fi

r="$(maxrel "$WORK/bohr.yields" "$WORK/coef.yields")"
if awk -v r="$r" 'BEGIN { exit !(r <= 1e-6) }'; then ok "(b) == (c): coefficient 0 gives Bohr's width, equal to the explicit c_B = $CB to $r"; else bad "(b) coefficient 0 differs from the explicit Bohr coefficient c_B = $CB by $r (expected <= 1e-6)"; paste "$WORK/bohr.yields" "$WORK/coef.yields" | head -8; fi
r="$(maxrel "$WORK/bohr.yields" "$WORK/off.yields")"
if awk -v r="$r" 'BEGIN { exit !(r > 1e-3) }'; then ok "(b) != (a): Bohr straggling moves the yield on the edge by up to $r"; else bad "(b) coefficient 0 is within $r of straggling off: no straggling applied"; fi
echo "  E_cm [MeV], yield off / Bohr / c_B on the rising edge (segment 1):"
paste "$WORK/off.yields" "$WORK/bohr.yields" "$WORK/coef.yields" | head -8 | awk '{ printf "    %s  %s  %s  %s\n", $1, $2, $4, $6 }'

# ---- (d) the 13N thick target with Bohr straggling --------------------------------------------------------------------
LINE_PREFIX='1  "1,2"  10  0  0  1  0  '
if grep -qF "$LINE_PREFIX" "$BASE13/13N.azr"; then
  mkdir -p "$WORK/thick/output" "$WORK/thick/checks"
  cp -r "$BASE13/data" "$WORK/thick/"
  awk -v pre="$LINE_PREFIX" '{ sub(/\r$/, ""); if (index($0, pre) == 1) { $23 = 1; $24 = 0 } print }' "$BASE13/13N.azr" > "$WORK/thick/t.azr"
  run thick
  rows="$(awk 'BEGIN { RS = ""; FS = "\n" } NR <= 2 { print }' "$WORK/thick/output/AZUREOut_aa=1_R=2.out" 2> /dev/null | awk 'NF >= 5 { print $4 }')"
  nt="$(printf '%s\n' "$rows" | grep -c .)"
  npos="$(printf '%s\n' "$rows" | awk '$1 > 0 && $1 == $1 + 0 && $1 < 1e300' | wc -l | tr -d ' ')"
  if [ "$nt" -ge 30 ] && [ "$npos" -eq "$nt" ]; then ok "(d) 13N thick target with Bohr straggling: every yield of segments 1 and 2 finite and positive ($nt points)"; else bad "(d) 13N thick target with Bohr straggling: $((nt - npos)) of $nt yields not positive/finite"; tail -3 "$WORK/thick/run.log"; fi
else
  bad "tests/13N/13N.azr no longer holds the thick-target line '$LINE_PREFIX...'"
fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
