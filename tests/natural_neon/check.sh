#!/usr/bin/env bash
#
# Natural neon: a real three-nucleus project (n+20Ne -> 21Ne, n+22Ne -> 23Ne,
# n+21Ne -> 22Ne in one file, natural abundances) taken from refit1, the model
# of record of the neon total cross sections (2026-10-18).  natural_neon.azr itself is a
# chi-squared suite (expected/chiSquared.out, run by run_tests.sh): natural-neon,
# enriched and folded data segments at the fitted parameters.  This script adds
# what a chi-squared cannot show, on the test segments (one lab-energy grid):
#   - the natural SUM over the three nuclei equals
#     0.9048 sigma(21Ne) + 0.0925 sigma(23Ne) + 0.0027 sigma(22Ne), each
#     nucleus's total taken from its own segment, to 1e-9;
#   - the decomposition file carries those three weighted terms and they add
#     up to the total.
# The <targetInt> lines are switched off for this run: they fold the test
# segments with the same numbers as the data segments they name.
#
#   ./tests/natural_neon/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/natural_neon.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# block FILE N -> the N-th blank-line-separated block of an output file (segments with the same pair keys share a file)
block() { awk -v blk="$2" 'BEGIN { RS = ""; FS = "\n" } NR == blk { print }' "$1"; }
# xs FILE N -> energy and cross section (the column before the S-factor) of block N
xs() { block "$1" "$2" | awk 'NF >= 4 { printf "%.12e %.12e\n", $1, $(NF - 1) }'; }
# maxrel A B -> largest relative difference between two two-column lists, row by row (the energies are c.m.
# energies of each entrance pair, so they differ between nuclei; the rows share one lab-energy grid)
maxrel() {
  paste "$1" "$2" | awk '
    NF >= 4 { d = $2 - $4; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300;
              r = d / m; if (r > worst) worst = r; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t + 0) }'; }
nlines() { awk 'END { print NR }' "$1"; }

# Test segments only, target integration off (first field of each <targetInt> line set to 0)
dir="$WORK/run"
mkdir -p "$dir/output" "$dir/checks"
ln -s "$HERE/data" "$dir/data"
awk '/<targetInt>/ { t = 1; print; next } /<\/targetInt>/ { t = 0 } t { sub(/^1 /, "0 ") } { print }' \
  "$HERE/natural_neon.azr" > "$dir/t.azr"
(cd "$dir" && printf '3\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)

f20="$dir/output/AZUREOut_aa=2_R=2.extrap"; f22="$dir/output/AZUREOut_aa=4_R=4.extrap"; f21="$dir/output/AZUREOut_aa=6_R=6.extrap"
if [ ! -s "$f20" ] || [ ! -s "$f22" ] || [ ! -s "$f21" ]; then
  bad "an extrapolation is missing"; tail -8 "$dir/run.log"
  echo "FAIL: $failures check(s)"; exit 1
fi
xs "$f20" 1 > "$WORK/nat"; xs "$f20" 2 > "$WORK/s21"; xs "$f22" 1 > "$WORK/s23"; xs "$f21" 1 > "$WORK/s22"
n="$(nlines "$WORK/nat")"
if [ "$n" -gt 300 ] && [ "$(nlines "$WORK/s21")" = "$n" ] && [ "$(nlines "$WORK/s23")" = "$n" ] && [ "$(nlines "$WORK/s22")" = "$n" ]; then
  ok "four test segments with $n points each"
else
  bad "test segments have $n / $(nlines "$WORK/s21") / $(nlines "$WORK/s23") / $(nlines "$WORK/s22") points"
fi

paste "$WORK/s21" "$WORK/s23" "$WORK/s22" | awk '{ printf "%.12e %.12e\n", $1, 0.9048 * $2 + 0.0925 * $4 + 0.0027 * $6 }' > "$WORK/sum"
r="$(maxrel "$WORK/nat" "$WORK/sum")"
if le "$r" 1e-9; then ok "natural SUM == 0.9048 sigma(21Ne) + 0.0925 sigma(23Ne) + 0.0027 sigma(22Ne) (max rel $r)"
else bad "natural SUM vs the weighted sum of the three nuclei (max rel $r)"; fi

dec="$dir/output/AZUREOut_sample_seg=1.extrap"
if [ -s "$dec" ]; then
  awk '!/^#/ && NF >= 5 { printf "%.12e %.12e\n", $1, $3 }' "$dec" > "$WORK/d21"
  awk '!/^#/ && NF >= 5 { printf "%.12e %.12e\n", $1, $4 }' "$dec" > "$WORK/d23"
  awk '!/^#/ && NF >= 5 { printf "%.12e %.12e\n", $1, $5 }' "$dec" > "$WORK/d22"
  awk '!/^#/ && NF >= 5 { t = $2; r = $2 - $3 - $4 - $5; if (r < 0) r = -r; if (t < 0) t = -t; if (r / t > w) w = r / t } END { printf "%.3e", w }' "$dec" > "$WORK/res"
  awk '{ printf "%.12e %.12e\n", $1, 0.9048 * $2 }' "$WORK/s21" > "$WORK/w21"
  awk '{ printf "%.12e %.12e\n", $1, 0.0925 * $2 }' "$WORK/s23" > "$WORK/w23"
  awk '{ printf "%.12e %.12e\n", $1, 0.0027 * $2 }' "$WORK/s22" > "$WORK/w22"
  r1="$(maxrel "$WORK/d21" "$WORK/w21")"; r2="$(maxrel "$WORK/d23" "$WORK/w23")"; r3="$(maxrel "$WORK/d22" "$WORK/w22")"
  res="$(cat "$WORK/res")"
  # the decomposition file is written with 10 significant digits
  if le "$r1" 1e-9 && le "$r2" 1e-9 && le "$r3" 1e-9 && le "$res" 1e-9; then
    ok "decomposition file carries the three weighted terms and they add up to the total"
  else
    bad "decomposition file: terms vs weighted totals $r1, $r2, $r3; total minus terms $res"
  fi
else
  bad "no decomposition file AZUREOut_sample_seg=1.extrap"; ls "$dir/output"
fi

if [ "$failures" -eq 0 ]; then echo "PASS: natural neon"; else echo "FAIL: $failures check(s)"; exit 1; fi
