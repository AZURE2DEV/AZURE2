#!/usr/bin/env bash
#
# Saved external-capture integrals must not be reused on a different grid.
#
# intEC.dat holds the EC amplitudes at the sub-point energies of one
# integration grid, for one set of inputs.  A channel radius, the grid
# (straggling, thickness, adaptive settings) can change the
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
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-900}")"

failures=0
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# make_variant NAME STRAGGLING RADIUS: a copy of the model with only data set 1
# (Artemov) active, the straggling of its thick target tgt_targetInt2 set (0:
# none, 1: the coefficient 0.04 keV^1/2 cm the old line carried), and the 12C+p channel radius set.  The points per width of
# the data sets on that target are lowered from 50 to 10 to keep the runs
# short -- the test compares the model with itself, so the grid need not be
# converged.  In a <levels> line, field 6 is the pair key and field 28 the
# channel radius.
make_variant() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$SRC/data" "$dir/data"
  awk -v strag="$2" -v radius="$3" '
    /<levels>/ { inlev = 1; print; next }
    /<\/levels>/ { inlev = 0 }
    inlev && NF > 27 && $6 == 1 { $28 = radius }
    /^[a-z_]+ [A-Za-z0-9_-]+$/ { rec = $0 }
    /^dataset / { n++; if (n > 1) { print; print "    active         no"; next } }
    $1 == "quadrature" && $NF == 50 && $(NF - 1) == "per_width" { $NF = 10; $0 = "    " $0 }
    rec == "target tgt_targetInt2" && $1 == "straggling" { $0 = "    straggling     " (strag == 1 ? "coefficient 0.04 keV^1/2 cm" : "none") }
    { print }' "$SRC/13N.azr" > "$dir/13N.azr"
}

# run NAME [EC_FILE] [FLAGS...]: mode 1, with or without a saved EC file.
# The saved file is passed relative to the run directory: the path is read
# from stdin, where MSYS does not translate POSIX paths, so an absolute
# /tmp/... path cannot be opened by the Windows binary.
run() {
  local dir="$WORK/$1" ec="${2:-}"
  [ -n "$ec" ] && ec="../$(basename "$ec")"
  shift 2 || shift $#
  (cd "$dir" && $RUN "$AZURE2_BIN" calculate ${ec:+--capture-amplitudes="$ec"} --no-gui "$@" 13N.azr > run.log 2>&1)
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
warned plain "WARNING: .*intEC\\|Could not open external capture file" && bad "a matching file was rejected or not found" || ok "a matching file is reused without a warning"
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

# Until 2026-11-16 --gsl-coul chose GSL's Coulomb functions, and this case showed that a file of the default routine
# was not reused with it.  The option is retired (docs/upgrade/coulomb.md section 9.1): the run is refused before
# anything is read, the saved file stays as it is, and the next run takes it as before (the signature of a run
# without the option did not change with the retirement).
echo "the retired option --gsl-coul"
cp "$WORK/plain/output/intEC.dat.sig" "$WORK/intEC.dat.sig"
cp "$WORK/intEC.dat" "$WORK/intEC.before"
run plain "$WORK/intEC.dat" --gsl-coul
warned plain "is no longer an option of AZURE2: leave it out" && warned plain "wrong by up to 5 % in the penetrability" \
  && ok "--gsl-coul is refused with the reason" || bad "--gsl-coul was not refused with the reason: $(sed -n 1,2p "$WORK/plain/run.log")"
cmp -s "$WORK/intEC.dat" "$WORK/intEC.before" && cmp -s "$WORK/plain/output/intEC.dat.sig" "$WORK/intEC.dat.sig" \
  && ok "and the saved file and its signature are untouched" || bad "the refused run changed the saved file"
run plain "$WORK/intEC.dat"
warned plain "WARNING: .*intEC\\|Could not open external capture file" && bad "after the refusal the matching file was rejected" || ok "and the next run reuses the file without a warning"

[ "$failures" -eq 0 ] && echo "  PASS" || echo "  FAIL ($failures)"
[ "$failures" -eq 0 ]
