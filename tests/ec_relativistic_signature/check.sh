#!/usr/bin/env bash
#
# Relativistic kinematics are part of the external-capture integral signature (audit 2026-10-25, checklist t3-25;
# 2026-10-27).
#
# With "relativistic=1" in <config> every channel's eta, rho and flux factor change, so the external-capture amplitudes
# change too (0.1-1 % through the penetrabilities).  Up to dev b80b857 the switch was not part of the signature written
# beside intEC.dat, so a relativistic run handed a file written without the switch read it back silently (the
# long-wavelength pattern of 41a468b), and the process-wide ECTable key lacked it too (tests/pyazr/ec_relativistic_test.py).
#
# Model: tests/13N cut down to its gas-target capture segment, as in tests/ec_integrals_signature.
#   (a) the relativistic run differs from the non-relativistic one (teeth, > 1e-4 relative);
#   (b) the relativistic run given the NON-relativistic intEC.dat warns that the file was computed for different inputs,
#       recomputes, and equals the fresh relativistic chi2 (1e-6);
#   (c) the non-relativistic run given its own file reuses it without a warning and equals its fresh chi2 (1e-4: the
#       amplitudes went through the text file and back, 1e-5 here).
# dev b80b857: (b) already held, because the signature hashes the c.m. energy of every sub-point and the switch moves
# those energies; the field added to the signature is belt and braces for a grid whose energies coincide.  The
# regression the old build DOES fail is the in-process node table, tests/pyazr/ec_relativistic_test.py.
#
#   ./tests/ec_relativistic_signature/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_relativistic.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# make_variant NAME REL: the 13N project with only its first data segment active, and relativistic=1 added when REL = 1
make_variant() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$SRC/data" "$dir/data"
  awk -v rel="$2" '
    { sub(/\r$/, "") }
    /<segmentsData>/ { inseg = 1; n = 0; print; next }
    /<\/segmentsData>/ { inseg = 0 }
    inseg && NF { n++; if (n > 1) $1 = 0 }
    /^<\/config>$/ && rel == 1 { print "relativistic=1" }
    { print }' "$SRC/13N.azr" > "$dir/13N.azr"
}
# run NAME [INTEGRALS-FILE]: mode 1; the second prompt is the parameter file (blank), the third the integrals file
run() {
  local dir="$WORK/$1" ec="${2:-}"
  (cd "$dir" && printf '1\n\n%s\n7\n' "$ec" | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline 13N.azr > run.log 2>&1)
}
chi2() { awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2> /dev/null; }
rel() { awk -v a="$1" -v b="$2" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; printf "%.3e", d }'; }
SIGWARN="was computed for a different integration grid or different inputs"

make_variant plain 0
make_variant relfresh 1
make_variant relstale 1
make_variant plainreuse 0
run plain; run relfresh
A="$(chi2 plain)"; B="$(chi2 relfresh)"
echo "  chi2: non-relativistic $A, relativistic $B"
[ -n "$A" ] && [ -n "$B" ] || { bad "a fresh run produced no chi-squared"; tail -3 "$WORK/relfresh/run.log"; echo "$failures check(s) FAILED"; exit 1; }
[ -f "$WORK/plain/output/intEC.dat.sig" ] || { bad "no output/intEC.dat.sig after the non-relativistic run"; echo "$failures check(s) FAILED"; exit 1; }
d="$(rel "$B" "$A")"
if awk -v d="$d" 'BEGIN { exit !(d > 1e-4) }'; then ok "(a) relativistic kinematics move the capture chi2 by $d (so a stale file would be visible)"; else bad "(a) the relativistic chi2 differs from the non-relativistic one by only $d"; fi

cp "$WORK/plain/output/intEC.dat" "$WORK/plain/output/intEC.dat.sig" "$WORK/relstale/"
run relstale "intEC.dat"
C="$(chi2 relstale)"
d="$(rel "$C" "$B")"
if grep -q "$SIGWARN" "$WORK/relstale/run.log" && awk -v d="$d" 'BEGIN { exit !(d <= 1e-6) }'; then ok "(b) the relativistic run refuses the non-relativistic intEC.dat (warned, recomputed) and gives the fresh relativistic chi2 ($C, $d)"; else bad "(b) the relativistic run given the non-relativistic intEC.dat gives $C vs fresh $B ($d); warning present: $(grep -c "$SIGWARN" "$WORK/relstale/run.log")"; fi

cp "$WORK/plain/output/intEC.dat" "$WORK/plain/output/intEC.dat.sig" "$WORK/plainreuse/"
run plainreuse "intEC.dat"
D="$(chi2 plainreuse)"
d="$(rel "$D" "$A")"
if ! grep -q "$SIGWARN" "$WORK/plainreuse/run.log" && awk -v d="$d" 'BEGIN { exit !(d <= 1e-4) }'; then ok "(c) the non-relativistic run reuses its own file without a warning ($D)"; else bad "(c) the non-relativistic run with its own file: $D vs $A ($d); warning present: $(grep -c "$SIGWARN" "$WORK/plainreuse/run.log")"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
