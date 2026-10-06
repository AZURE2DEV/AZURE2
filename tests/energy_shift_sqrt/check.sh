#!/usr/bin/env bash
#
# sqrt(E) energy-shift term: (1) shifting a segment's energies with the
# "sqrtshift" block must equal reading the same yields at the shifted energies;
# (2) a MIGRAD fit must recover a known (a, b) from mislabelled energies, both
# with Minuit's numerical gradient and with --use-gradient.
#
#   ./tests/energy_shift_sqrt/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/eshift_sqrt.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
seg() { awk -F, -v k="$2" '$1 == k { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$1/output/chiSquared.out" 2>/dev/null; }
# relclose A B TOL -> 0 when |A-B| <= TOL*max(|B|,1e-300)
relclose() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m == 0) m = 1e-300; exit !(d <= t * m) }'; }

# ---- 1. the shift equals re-labelled data ---------------------------------
mkdir -p "$WORK/calc/output" "$WORK/calc/checks"
cp -R "$HERE/data" "$WORK/calc/data"
cp "$HERE/energy_shift_sqrt.azr" "$WORK/calc/t.azr"
(cd "$WORK/calc" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
if [ ! -f "$WORK/calc/output/chiSquared.out" ]; then
  bad "calculation produced no chiSquared.out"; tail -20 "$WORK/calc/run.log"; exit 1
fi
s1="$(seg "$WORK/calc" 1)"; s2="$(seg "$WORK/calc" 2)"; s3="$(seg "$WORK/calc" 3)"; s4="$(seg "$WORK/calc" 4)"; s5="$(seg "$WORK/calc" 5)"
if relclose "$s1" "$s2" 1e-8; then ok "a + b*sqrt(E) shift == pre-shifted data, cross section (chi2 $s1 vs $s2)"; else bad "cross section: shifted $s1 != pre-shifted $s2"; fi
if relclose "$s3" "$s4" 1e-8; then ok "b*sqrt(E) shift == pre-shifted data, analyzing power (chi2 $s3 vs $s4)"; else bad "analyzing power: shifted $s3 != pre-shifted $s4"; fi
if relclose "$s5" 24 1e-6; then ok "untouched segment unchanged (chi2 $s5)"; else bad "untouched segment moved: $s5 (expected 24)"; fi
if relclose "$s1" 24 1e-2; then bad "shift had no effect on segment 1 (chi2 $s1 = pristine 24)"; else ok "shift changes segment 1 (chi2 $s1 vs pristine 24)"; fi

# ---- 2. recover (a, b) by fitting -----------------------------------------
# Truth: labels L obey L + A0 + B0*sqrt(L) = E_true with A0 = -0.006, B0 = 0.005,
# and the yields are the model's exact values, so the truth is the chi2 = 0
# minimum.  Fitted twice: with Minuit's numerical gradient, and with
# --use-gradient (AZURECalc::Gradient, where both shift terms are finite
# differenced per segment); both must land on the truth.
fit_variant() {  # NAME EXTRA_FLAGS
  local dir="$WORK/fit_$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/data" "$dir/data"
  awk '
    /<levels>/ { inlev = 1; print; next }
    /<\/levels>/ { inlev = 0 }
    inlev && NF > 12 { $4 = 1; $11 = 1; print; next }          # fix every energy and width
    /<segmentsData>/ { print; print "1 1 1 0 100 40 40 4 1 0 0 0 1 1 data/xs_40_mislabelled.dat 0 0 sqrtshift 0 1 1"; skip = 1; next }
    /<\/segmentsData>/ { skip = 0 }
    skip { next }
    { print }' "$HERE/energy_shift_sqrt.azr" > "$dir/t.azr"
  (cd "$dir" && printf '2\n\n\n\n7\n' | "$AZURE2_BIN" $2 --no-gui --no-readline t.azr > run.log 2>&1)
  if [ ! -f "$dir/output/param.sav" ]; then
    bad "$1: fit produced no param.sav"; tail -20 "$dir/run.log"; return
  fi
  local a b t
  a="$(awk '$1 == "segment_1_energy_shift" { print $2 }' "$dir/output/param.sav")"
  b="$(awk '$1 == "segment_1_energy_shift_sqrt" { print $2 }' "$dir/output/param.sav")"
  t="$(total "$dir")"
  if awk -v a="$a" 'BEGIN { d = a + 0.006; if (d < 0) d = -d; exit !(d < 1e-4) }'; then ok "$1: fitted constant shift a = $a (truth -0.006)"; else bad "$1: fitted constant shift a = $a, expected -0.006"; fi
  if awk -v b="$b" 'BEGIN { d = b - 0.005; if (d < 0) d = -d; exit !(d < 1e-4) }'; then ok "$1: fitted sqrt(E) coefficient b = $b (truth 0.005)"; else bad "$1: fitted sqrt(E) coefficient b = $b, expected 0.005"; fi
  if awk -v t="$t" 'BEGIN { exit !(t + 0 < 1e-3) }'; then ok "$1: fit reaches chi2 = 0 at the truth (chi2 $t)"; else bad "$1: fit chi2 $t, expected ~0"; fi
}
fit_variant numerical ""
fit_variant analytic "--use-gradient"

if [ "$failures" -eq 0 ]; then echo "PASS: sqrt(E) energy shift"; else echo "FAIL: $failures check(s)"; exit 1; fi
