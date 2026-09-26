#!/usr/bin/env bash
#
# A Gaussian fold across a resonance far narrower than the Gaussian must give
# the fold of the unfolded curve -- pointwise, and with the resonance's area
# conserved -- and must keep doing so when a fit moves the resonance.
#
# Model: tests/17O (a THM segment, 21 keV lab Gaussian; its 5- level at
# E_cm = 74 keV is 46 eV wide).  Three checks:
#
#  A. The <targetInt> fold at 38 energies across the resonance against a
#     trapezoid fold, done here in awk, of AZURE2's own UNFOLDED curve on a
#     0.5 eV / 20 eV grid, over the same +-5 sigma window: once with the
#     default grid settings and once with the "5 50" (resonanceWidthMultiplier,
#     pointsPerWidth) of tests/17O.  Before the lattice had geometric tails the
#     second was 16 % off at 150 sub-points -- the chord over the resonance's
#     1/(E-E_R)^2 wings, bridged with the smooth step, lies above them.
#  B. The same energies with the level energy supplied through an external
#     parameter file 1.5 keV (33 widths) above the .azr value -- the grid is
#     built at Fill time from the .azr, as it is when a fit moves a level --
#     against a run whose .azr carries the moved energy.  The stale grid put
#     the peak between coarse sub-points (68 % off); the grid is now rebuilt
#     from the current parameters (EPoint::RefreshSubPointGrid) and the two
#     agree to rounding.
#
#   ./tests/thm_narrow_fold/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../17O"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-900}"; else RUN=""; fi
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_narrow_fold.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
fail=0

SIGMA_LAB=0.021            # the targetInt sigma of tests/17O (lab)
SIGMA_CM=$(awk -v s=$SIGMA_LAB 'BEGIN { printf "%.12g", s * 17 / 18 }')  # n + 17O, masses 1 and 17
LEVEL_EX=1.891             # the narrow 5- level (E_cm = 74 keV)
TOL_POINT=1e-3             # pointwise, relative
TOL_AREA=5e-4              # area of the folded curve across the resonance

# project DIR DATAFILE TARGETINT-LINE LEVEL-EX: a copy of tests/17O reading
# DATAFILE (lab energies) with the given targetInt line (empty = unfolded),
# the narrow level at LEVEL-EX with a free energy, and no <parameterSettings>
# overrides (so an external parameter file is not overridden in turn).
project() {
  local d="$1" data="$2" tline="$3" ex="$4"
  mkdir -p "$d/data" "$d/output" "$d/checks"
  cp "$data" "$d/data/pts.dat"
  awk -v tline="$tline" -v ex="$ex" -v lev="$LEVEL_EX" '
    /^<segmentsData>/ { print; print "1 1 2 -1 1 0 180 10 1.0 0 0 0 0 0 data/pts.dat 0 0"; skip = 1; next }
    /^<\/segmentsData>/ { skip = 0 }
    /^<targetInt>/ { print; if (tline != "") print tline; skip = 1; next }
    /^<\/targetInt>/ { skip = 0 }
    /^<parameterSettings>/ { print; skip = 1; next }
    /^<\/parameterSettings>/ { skip = 0 }
    skip { next }
    NF > 30 && $3 == lev { $3 = ex; $4 = 0 }
    { print }' "$SRC/17O.azr" > "$d/run.azr"
}

# run DIR [PARFILE]: calculate; prints nothing, leaves output/
run() {
  (cd "$1" && printf '1\n%s\n\n7\n' "${2:-}" | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 \
     | head -c 1000000 > log)
  [ -s "$1/output/AZUREOut_aa=1_R=2.out" ] || { echo "  FAIL  no output in $(basename "$1")"; tail -3 "$1/log" | sed 's/^/        /'; exit 1; }
}
curve() { awk 'NF { print $1, $4 }' "$1/output/AZUREOut_aa=1_R=2.out"; }   # E_cm, fit

# Evaluation energies, E_cm 20-131 keV every 3 keV, and the dense grid (0.5 eV
# within 2 keV of the level, 20 eV elsewhere, covering +-5 sigma beyond them);
# data files hold lab energies, E_lab = 18/17 E_cm.
awk 'BEGIN { for (e = 0.020; e < 0.1315; e += 0.003) printf "%.9f 0 1 0.1\n", e * 18 / 17 }' > "$WORK/eval.dat"
awk -v s="$SIGMA_CM" 'BEGIN {
  lo = 0.020 - 5.2 * s; hi = 0.131 + 5.2 * s; er = 0.074
  for (e = lo; e <= hi; e += 2.0e-5) if (e < er - 0.002 || e > er + 0.002) printf "%.10f 0 1 0.1\n", e * 18 / 17
  for (e = er - 0.002; e <= er + 0.002; e += 5.0e-7) printf "%.10f 0 1 0.1\n", e * 18 / 17
}' | sort -g > "$WORK/dense.dat"

TL_DEFAULT="1 \"1\" 150 1 $SIGMA_LAB 0 0 \"\" 0 0 0 0 \"\" 0"
TL_17O="$TL_DEFAULT 0 0.04 5 50"

project "$WORK/ref" "$WORK/dense.dat" "" "$LEVEL_EX";         run "$WORK/ref"
project "$WORK/def" "$WORK/eval.dat" "$TL_DEFAULT" "$LEVEL_EX"; run "$WORK/def"
project "$WORK/x5"  "$WORK/eval.dat" "$TL_17O" "$LEVEL_EX";     run "$WORK/x5"

# compare NAME FOLDED: pointwise and area comparison with the reference fold
compare() {
  curve "$WORK/ref" | sort -g > "$WORK/ref.txt"
  curve "$2" | awk -v s="$SIGMA_CM" -v tp="$TOL_POINT" -v ta="$TOL_AREA" -v name="$1" '
    FNR == NR { re[++n] = $1; rs[n] = $2; next }
    { x[++m] = $1; y[m] = $2 }
    END {
      c = 1 / (sqrt(2 * 3.14159265358979) * s)
      for (i = 1; i <= m; i++) {
        f = 0; prev = 0
        for (j = 1; j <= n; j++) {
          if (re[j] < x[i] - 5 * s || re[j] > x[i] + 5 * s) { prev = 0; continue }
          g = rs[j] * c * exp(-0.5 * ((re[j] - x[i]) / s) ^ 2)
          if (prev) f += 0.5 * (g + gp) * (re[j] - ep)
          prev = 1; gp = g; ep = re[j]
        }
        r = y[i] / f - 1; if (r < 0) r = -r
        if (r > worst) { worst = r; at = x[i] }
        ref[i] = f
      }
      for (i = 2; i <= m; i++) { a += 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1]); b += 0.5 * (ref[i] + ref[i-1]) * (x[i] - x[i-1]) }
      ar = a / b - 1; aa = ar < 0 ? -ar : ar
      printf "  %-34s max |rel err| %.2e (at E_cm %.1f keV), area ratio - 1 = %+.2e\n", name, worst, at * 1000, ar
      exit !(worst < tp && aa < ta)
    }' "$WORK/ref.txt" - || { echo "  FAIL  $1: the fold does not match the reference"; fail=1; }
}
echo "  reference: unfolded curve on $(wc -l < "$WORK/dense.dat") energies, sigma_cm = $SIGMA_CM MeV"
compare "A. default grid" "$WORK/def"
compare "A. tests/17O grid (5 widths, 50/width)" "$WORK/x5"

# B. level moved through an external parameter file
MOVED=$(awk -v e="$LEVEL_EX" 'BEGIN { printf "%.6f", e + 0.0015 }')
project "$WORK/moved" "$WORK/eval.dat" "$TL_DEFAULT" "$MOVED"; run "$WORK/moved"
project "$WORK/par" "$WORK/eval.dat" "$TL_DEFAULT" "$LEVEL_EX"; run "$WORK/par" "$WORK/moved/output/param.par"
paste <(curve "$WORK/par") <(curve "$WORK/moved") | awk '
  { r = $2 / $4 - 1; if (r < 0) r = -r; if (r > w) w = r }
  END { printf "  %-34s max |rel diff| %.2e against a grid built at the moved energy\n", "B. level moved 1.5 keV by param file", w
        exit !(w < 1e-7) }' || { echo "  FAIL  B: the sub-point grid did not follow the moved level"; fail=1; }

if [ "$fail" = 0 ]; then echo "  ok    narrow-resonance fold matches the dense reference and follows a moved level"; fi
exit "$fail"
