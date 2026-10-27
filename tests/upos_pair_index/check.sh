#!/usr/bin/env bash
#
# The UPOS secondary-gamma distribution takes the residual nucleus's spin from the EXIT PAIR, not from the pair whose
# number happens to equal the decay index (audit item 10, 2026-10-26).  GenMatrixFunc::CalculateBareCrossSection and
# the AMatrixFunc gradient mirror looked up GetPair(ir) with ir the decay index within the entrance pair's list; the two
# agree only while every earlier pair has a decay from the entrance pair.  Prepending a pair of another compound nucleus
# (p + 12C -> 13N, key 1, one 1/2- level) to tests/secondary_gamma/models/spin2_plus.azr (alpha + 13C -> n + 16O(3-),
# keys 2 and 5) shifts the pair numbers by one while the decay indices stay: the old engine then took j2f = 1/2 (13C, the
# alpha pair's own target) instead of 3 (16O*), and the 13C(a,n'g) distribution changed.
#
#   (a) the UPOS angular distribution (segment 1, AZUREOut_aa=2_R=5.out block 1) of the two models is identical;
#   (b) so are the angle-integrated and the particle angular distributions (blocks 2 and 3), which never saw j2f;
#   (c) the added nucleus's own elastic output exists (the model really has two compound nuclei).
#
#   ./tests/upos_pair_index/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SG="$HERE/../secondary_gamma"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/upos_pair_index.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# the extra compound nucleus: p + 12C (key 1), one 1/2- level at E_x 3.0 MeV with a 1 keV proton width, no capture
EXTRA='0.5  -1  3.0  1  1  1  1  2  9  1  1  1000.0  0.5  1  0  1  0  1.00728  12  1  6  1.94351  1.94351  0  0  0.0  0  3.6  0  0  0'
mkdir -p "$WORK/base/output" "$WORK/base/checks" "$WORK/extra/output" "$WORK/extra/checks"
cp "$SG/models/spin2_plus.azr" "$WORK/base/t.azr" && cp -r "$SG/models/data" "$WORK/base/"
awk -v extra="$EXTRA" '{ print } /^<levels>$/ { print extra; print "" }' "$SG/models/spin2_plus.azr" > "$WORK/extra/t.azr"
cp -r "$SG/models/data" "$WORK/extra/"
grep -q "^0.5  -1  3.0  1  1  1" "$WORK/extra/t.azr" || { bad "could not build the two-nucleus model"; exit 1; }
for m in base extra; do
  (cd "$WORK/$m" && printf '1\n\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  [ -s "$WORK/$m/output/AZUREOut_aa=2_R=5.out" ] || { bad "$m: no output"; tail -8 "$WORK/$m/run.log"; }
done

# block FILE N -> the N-th blank-line-separated block (the three segments of the file, in order)
block() { awk -v blk="$2" 'BEGIN { RS = ""; FS = "\n" } NR == blk { print }' "$1"; }
# maxrel of the fit column (4) between two blocks
maxrel() {
  paste <(block "$1" "$3" | awk '{ print $1, $3, $4 }') <(block "$2" "$3" | awk '{ print $1, $3, $4 }') | awk '
    { if ($1 != $4 || $2 != $5) { print "grid"; exit }
      d = $3 - $6; if (d < 0) d = -d; m = $6 < 0 ? -$6 : $6; if (m < 1e-300) m = 1e-300; r = d / m; if (r > w) w = r; n++ }
    END { if (n == 0) print "nan"; else printf "%.3e\n", w }'
}
le() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r + 0 <= t + 0) }'; }
f1="$WORK/base/output/AZUREOut_aa=2_R=5.out"; f2="$WORK/extra/output/AZUREOut_aa=2_R=5.out"
if [ -s "$f1" ] && [ -s "$f2" ]; then
  n="$(block "$f1" 1 | grep -c .)"
  r="$(maxrel "$f1" "$f2" 1)"
  if le "$r" 1e-10; then ok "(a) UPOS angular distribution unchanged by a preceding pair of another nucleus (max rel $r, $n points)"
  else bad "(a) UPOS angular distribution changed by a preceding pair of another nucleus (max rel $r)"; fi
  r2="$(maxrel "$f1" "$f2" 2)"; r3="$(maxrel "$f1" "$f2" 3)"
  if le "$r2" 1e-10 && le "$r3" 1e-10; then ok "(b) angle-integrated and particle distributions unchanged ($r2, $r3)"
  else bad "(b) angle-integrated / particle distributions changed ($r2, $r3)"; fi
  if grep -q "ERROR" "$WORK/extra/run.log"; then bad "(c) the two-nucleus model reported an error"; grep ERROR "$WORK/extra/run.log" | head -2
  elif grep -q "Z=7 A=13\|13N\|two compound\|compound nuclei" "$WORK/extra/run.log" || [ "$(grep -c "" "$WORK/extra/output/parameters.out")" -gt "$(grep -c "" "$WORK/base/output/parameters.out")" ]; then
    ok "(c) the model has the extra nucleus (parameters.out longer: $(grep -c "" "$WORK/extra/output/parameters.out") vs $(grep -c "" "$WORK/base/output/parameters.out") lines)"
  else bad "(c) the extra nucleus is not in the model"; fi
fi

echo
if [ "$failures" -eq 0 ]; then echo "upos_pair_index: PASS"; else echo "upos_pair_index: FAIL ($failures)"; exit 1; fi
