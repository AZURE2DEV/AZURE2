#!/usr/bin/env bash
#
# The thin-target energy loss is the exact range integral, not the first-order epsilon(E0) N
# (audit 2026-10-25, tier 3 entry 2; 2026-10-27).
#
# The target-integration yield is Y(E0) = (1/N) int_{E0-Delta}^{E0} sigma(E)/epsilon(E) dE.  Up to dev 5ef20e5 the energy
# loss Delta was epsilon(E0) N; with a stopping cross section that rises through the target (every ion above the Bragg
# peak) that puts the back edge too shallow, and since 1/epsilon sits in the integrand a FLAT cross section came out low
# by (1/N) int dE/epsilon < 1: 2 % for a 5 % target, 8 % for a 20 % one, more at lower energy (a shape distortion of an
# excitation function).  Now TargetEffect::TargetThickness solves int_{E0-Delta}^{E0} dE/epsilon = N.
#
# Model: tests/thick_target_resonance (p + 14N, external capture off) with two data segments of the same observable
# (total capture): segment 1 at lab 320, 360, 400, 450 and 500 keV under ONE <targetInt> line -- no beam spread,
# epsilon = a0 / x^0.8 MeV cm^2 (a0 = 7.95e-21: 18 eV/1e15 at 360 keV), N = 4e18 atoms/cm^2 (a 20 % target at
# 360 keV), 20 sub-points -- and segment 2 the bare cross section on a 0.1 keV lab grid from 200 to 520 keV, no target
# effect.  Truth: segment 2 integrated by hand (awk) over the EXACT window [E_exit, E0] with
# E_exit = (E0^1.8 - 1.8 a0 N)^(1/1.8) (closed form for the power law), trapezoid with the window edges interpolated.
# The 320 keV point has the 259.7 keV c.m. resonance inside its window, the others are on the smooth background above it.
#   (a) every yield within 0.5 % of the truth (1 % for the resonance point);
#   (b) the first-order yields (the truth integrated over [E0 - eps(E0) N, E0]) are more than 3 % below the exact ones
#       at every background point (4 % at 500 keV, 9 % at 400 keV, 65 % at 360 keV where the back edge of the exact
#       window reaches the resonance tail; a flat cross section would give 8 %), so the check has teeth.
# dev 5ef20e5: the background points come out 4-65 % low (the first-order window), (a) fails there.
#
#   ./tests/target_thickness_exact/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../thick_target_resonance"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/target_thickness_exact.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

A0=7.95e-21
A1=0.8
DENS=4.0e18
F="$(awk 'BEGIN { printf "%.10f", 14.0031 / (1.00783 + 14.0031) }')"   # c.m. / lab
ENERGIES="0.320 0.360 0.400 0.450 0.500"

# the resonance project with two data segments (the thin target; the bare cross section on a fine grid) and no test
# segment; the <targetInt> block holds one line for segment 1 only
mkdir -p "$WORK/m/output" "$WORK/m/checks" "$WORK/m/data"
for e in $ENERGIES; do printf '%s 0 1.0e-06 1.0e-07\n' "$e"; done > "$WORK/m/data/thin.dat"
awk 'BEGIN { for (i = 0; i <= 3200; i++) printf "%.4f 0 1.0e-06 1.0e-07\n", 0.20 + i * 0.0001 }' > "$WORK/m/data/fine.dat"
awk -v tline="1  \"1\"  20  0  0  1  $DENS  \"a0/x^a1\"  2  $A0  $A1  0  0  0  \"\"  0  0  0  20  50" \
    -v s1="1 1 -1 0 6 0 180 3 1 0 0 0 0 0 data/thin.dat 0 0" \
    -v s2="1 1 -1 0 6 0 180 3 1 0 0 0 0 0 data/fine.dat 0 0" '
  { sub(/\r$/, "") }
  /^<targetInt>$/ { print; print tline; skip = 1; next }
  /^<segmentsData>$/ { print; print s1; print s2; skip = 1; next }
  /^<segmentsTest>$/ { print; skip = 1; next }
  /^<\/targetInt>$|^<\/segmentsData>$|^<\/segmentsTest>$/ { skip = 0 }
  !skip { print }' "$BASE/thick_target_resonance.azr" > "$WORK/m/t.azr"
(cd "$WORK/m" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
OUT="$WORK/m/output/AZUREOut_aa=1_TOTAL_CAPTURE.out"
# block K of the output file -> "E_cm fit" rows
block() { awk -v k="$1" 'BEGIN { RS = ""; FS = "\n" } NR == k { print }' "$OUT" 2> /dev/null | awk '{ sub(/\r$/, "") } NF >= 5 { printf "%.10f %.12e\n", $1, $4 }'; }
block 1 > "$WORK/yields"
block 2 > "$WORK/sigma"
n="$(wc -l < "$WORK/yields" | tr -d ' ')"; ns="$(wc -l < "$WORK/sigma" | tr -d ' ')"
if [ "$n" -eq 5 ] && [ "$ns" -eq 3201 ]; then ok "the run wrote the 5 yields and the bare cross section on 3201 energies"; else bad "the run wrote $n yields and $ns cross sections (expected 5 and 3201)"; tail -3 "$WORK/m/run.log"; fi

# truth E0 EXIT-MODE -> (1/N) int_{E_exit}^{E0} sigma / epsilon dE_lab over the c.m. table; EXIT-MODE exact | first
truth() {
  awk -v e0="$1" -v mode="$2" -v a0="$A0" -v a1="$A1" -v dens="$DENS" -v f="$F" '
    { E[NR] = $1; S[NR] = $2 }
    END {
      if (mode == "exact") ex = (e0 ^ (1 + a1) - (1 + a1) * a0 * dens) ^ (1 / (1 + a1)); else ex = e0 - a0 / e0 ^ a1 * dens;
      lo = f * ex; hi = f * e0;              # c.m. window
      # integrand per c.m. energy: sigma(Ecm) / eps_lab(Ecm/f) / f  (dE_lab = dE_cm / f)
      sum = 0;
      for (i = 1; i < NR; i++) {
        a = E[i]; b = E[i + 1]; if (b <= lo || a >= hi) continue;
        ya = S[i]; yb = S[i + 1];
        if (a < lo) { ya = ya + (yb - ya) * (lo - a) / (b - a); a = lo }
        if (b > hi) { yb = S[i] + (S[i + 1] - S[i]) * (hi - E[i]) / (E[i + 1] - E[i]); b = hi }
        ga = ya * (a / f) ^ a1 / a0 / f; gb = yb * (b / f) ^ a1 / a0 / f;
        sum += 0.5 * (ga + gb) * (b - a);
      }
      printf "%.12e", sum / dens;
    }' "$WORK/sigma"
}

i=0
teeth=0
for e in $ENERGIES; do
  i=$((i + 1))
  y="$(sed -n "${i}p" "$WORK/yields" | awk '{ print $2 }')"
  [ -n "$y" ] || { bad "(a) E_lab = $e MeV: no yield"; continue; }
  tx="$(truth "$e" exact)"; tf="$(truth "$e" first)"
  tol=0.005; [ "$e" = "0.320" ] && tol=0.01
  rel="$(awk -v a="$y" -v b="$tx" 'BEGIN { printf "%.4f", (b > 0) ? a / b - 1 : 9 }')"
  rf="$(awk -v a="$tf" -v b="$tx" 'BEGIN { printf "%.4f", (b > 0) ? a / b - 1 : 9 }')"
  if awk -v r="$rel" -v t="$tol" 'BEGIN { exit !(r < t && r > -t) }'; then ok "(a) E_lab = $e MeV: yield $y vs exact truth $tx ($rel; the first-order window would give $rf)"; else bad "(a) E_lab = $e MeV: yield $y vs exact truth $tx ($rel, tolerance $tol; the first-order window gives $rf)"; fi
  if [ "$e" != "0.320" ]; then
    if awk -v r="$rf" 'BEGIN { exit !(r < -0.03) }'; then teeth=$((teeth + 1)); else bad "(b) E_lab = $e MeV: the first-order window is only $rf off the exact one, expected below -3 %"; fi
  fi
done
[ "$teeth" -eq 4 ] && ok "(b) the first-order window is more than 3 % below the exact one at the four background points"

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
