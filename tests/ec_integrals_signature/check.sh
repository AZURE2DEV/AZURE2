#!/usr/bin/env bash
#
# Saved external-capture integrals must not be reused on a different grid.
#
# intEC.dat holds the EC amplitudes at the sub-point energies of one
# integration grid, for one set of inputs.  A channel radius, the grid
# (straggling, thickness, adaptive settings) or --gsl-coul can change the
# amplitudes without changing how many there are, so the amplitude count alone
# cannot tell the files apart, and a stale file used to be read back silently.
# The signature written beside it (intEC.dat.sig) is what catches that.
#
# Model: tests/13N cut down to one gas-target capture segment (Artemov), which
# carries a target integration with external capture.
#
#   ./tests/ec_integrals_signature/check.sh path/to/AZURE2
#
# Exit status 0 when every check passes.

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
# Absolute, since every run happens from inside its own project directory.
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_signature.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# make_variant NAME STRAGGLING RADIUS: a copy of the model with only segment 1
# (Artemov) active, the straggling flag of the "1,2" target effect set, and the
# 12C+p channel radius set.  The trailing fields of a <targetInt> line are
# isStraggling stragglingCoefficient widthMultiplier pointsPerWidth; the last
# is lowered from 50 to 10 to keep the runs short -- the test compares the
# model with itself, so the grid need not be converged.  In a <levels> line,
# field 6 is the pair key and field 28 the channel radius.
make_variant() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$SRC/data" "$dir/data"
  awk -v strag="$2" -v radius="$3" '
    /<levels>/ { inlev = 1; print; next }
    /<\/levels>/ { inlev = 0 }
    inlev && NF > 27 && $6 == 1 { $28 = radius }
    /<segmentsData>/ { inseg = 1; n = 0; print; next }
    /<\/segmentsData>/ { inseg = 0 }
    inseg && NF { n++; if (n > 1) $1 = 0 }
    /<targetInt>/ { inte = 1; print; next }
    /<\/targetInt>/ { inte = 0 }
    inte && $2 == "\"1,2\"" { $(NF - 3) = strag; $NF = 10 }
    { print }' "$SRC/13N.azr" > "$dir/13N.azr"
}

# run NAME [EC_FILE] [FLAGS...]: mode 1, with or without a saved EC file.
run() {
  local dir="$WORK/$1" ec="${2:-}"
  shift 2 || shift $#
  (cd "$dir" && printf '1\n\n%s\n7\n' "$ec" |
    "$AZURE2_BIN" --no-gui --no-readline "$@" 13N.azr > run.log 2>&1)
}

chi2() { awk '/^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2>/dev/null; }

# Reading a file back is not bit-exact: intEC.dat stores each amplitude to six
# significant digits, which moves chi2 in the fifth or sixth.  A recomputed
# file is compared exactly; a reused one to this relative tolerance.
close() { awk -v a="$1" -v b="$2" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; exit !(a != "" && d < 1e-4) }'; }
warned() { grep -q "$2" "$WORK/$1/run.log"; }
SIGWARN="different integration grid or different inputs"

make_variant plain 0 3.6
make_variant straggling 1 3.6
make_variant radius 0 4.0

echo "fresh calculations"
run plain
run straggling
run radius
A="$(chi2 plain)"
B="$(chi2 straggling)"
C="$(chi2 radius)"
echo "  chi2: as is $A, with straggling $B, radius 4.0 fm $C"
[ -n "$A" ] && [ -n "$B" ] && [ -n "$C" ] || { echo "  FAIL  a fresh run produced no chi-squared"; exit 1; }
[ -f "$WORK/plain/output/intEC.dat.sig" ] && ok "the signature file is written" || bad "no output/intEC.dat.sig after a fresh run"
[ "$A" != "$C" ] && ok "the radius changes the result" || bad "the radius does not change chi2, so this test cannot tell a stale file from a good one"
if [ "$(wc -l < "$WORK/plain/output/intEC.dat")" -eq "$(wc -l < "$WORK/radius/output/intEC.dat")" ]; then
  ok "both radii give the same number of amplitudes (the count check cannot see the change)"
else
  bad "the two radii give files of different size; the signature is not what is being tested"
fi

cp "$WORK/plain/output/intEC.dat" "$WORK/intEC.dat"
cp "$WORK/plain/output/intEC.dat.sig" "$WORK/intEC.dat.sig"

echo "reusing the file on the grid it was built for"
run plain "$WORK/intEC.dat"
warned plain "WARNING: .*intEC" && bad "a matching file was rejected" || ok "a matching file is reused without a warning"
close "$(chi2 plain)" "$A" && ok "and gives the same chi2 ($(chi2 plain), to the file's precision)" || bad "reused file gives chi2 $(chi2 plain), fresh $A"

echo "reusing the file with a different channel radius"
run radius "$WORK/intEC.dat"
warned radius "$SIGWARN" && ok "the signature mismatch is reported" || bad "no signature warning; the stale file was used"
[ "$(chi2 radius)" = "$C" ] && ok "and the recomputed chi2 equals the fresh one ($C)" || bad "chi2 $(chi2 radius), fresh calculation $C"

# Straggling widens the kernel, and here also the depth the grid spans only
# where the target is thick; whether it moves the sub-point energies depends
# on the model.  Either way the answer must be the fresh one: reused when the
# grid is unchanged, recomputed when it is not.
echo "reusing the file with straggling switched on"
run straggling "$WORK/intEC.dat"
if warned straggling "$SIGWARN"; then
  echo "  note  the grid moved, the file was recomputed"
  [ "$(chi2 straggling)" = "$B" ] && ok "chi2 equals the fresh one" || bad "chi2 $(chi2 straggling), fresh calculation $B"
else
  echo "  note  the grid did not move, the file was reused"
  close "$(chi2 straggling)" "$B" && ok "chi2 equals the fresh one to the file's precision" || bad "chi2 $(chi2 straggling), fresh calculation $B"
fi

echo "a file with no signature (written before signatures existed)"
rm "$WORK/intEC.dat.sig"
run plain "$WORK/intEC.dat"
warned plain "has no signature file" && ok "the missing signature is reported" || bad "no warning for a file without a signature"
[ "$(chi2 plain)" = "$A" ] && ok "and the recomputed chi2 equals the fresh one" || bad "chi2 $(chi2 plain), fresh $A"

echo "the same grid with the other Coulomb-function routine"
cp "$WORK/plain/output/intEC.dat.sig" "$WORK/intEC.dat.sig"
run plain "$WORK/intEC.dat" --gsl-coul
warned plain "$SIGWARN" && ok "switching to --gsl-coul invalidates the file" || bad "--gsl-coul reused integrals computed with the other routine"

[ "$failures" -eq 0 ] && echo "  PASS" || echo "  FAIL ($failures)"
[ "$failures" -eq 0 ]
