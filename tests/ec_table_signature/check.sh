#!/usr/bin/env bash
#
# The Coulomb table a run reads is part of the external-capture integral signature written beside intEC.dat.
#
# The integrals take the Coulomb waves at their radial nodes from the table when one is used (CoulFunc), so a file
# computed with one table (or with the exact functions) holds other digits than one computed with another.  The
# signature names the table by what it holds (CoulombTable::ContentIdentity: its header and a hash of its values),
# not by its file: the same table at another path is the same table.
#
# Model: tests/13N cut down to its gas-target capture segment, as in tests/ec_relativistic_signature.  The table is
# $AZURE2_COULOMB_TABLE (required); the second table, with other values, is built here with the table tool
# (azure2-coulomb-table beside the engine, or $AZURE2_COULOMB_TABLE_TOOL; without it case (c) is skipped).
#   (a) a run with the table given the intEC.dat of a run with the same table reuses it without a warning;
#   (b) a run with a byte copy of the table at another path (another modification time) reuses it too;
#   (c) a run with a coarser table (other values) given that file refuses it (signature warning) and recomputes;
#   (d) a run with --no-coulomb-table given that file refuses it and recomputes;
#   (e) a run from scratch with the copy writes the same output/intEC.dat.sig as the first run;
#   (f) a copy with one tabulated value changed in its last bit (the header the same) refuses it.
#
#   ./tests/ec_table_signature/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
if [ -z "${AZURE2_COULOMB_TABLE:-}" ] || [ ! -f "$AZURE2_COULOMB_TABLE" ]; then
  echo "SKIP: AZURE2_COULOMB_TABLE is not set to a table file"
  exit 0
fi
TABLE="$AZURE2_COULOMB_TABLE"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_table_sig.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

make_variant() {
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$SRC/data" "$dir/data"
  awk '{ sub(/\r$/, "") } /^dataset / { n++; if (n > 1) { print; print "    active         no"; next } } { print }' \
    "$SRC/13N.azr" > "$dir/13N.azr"
}
# run NAME TABLE [INTEGRALS-FILE] [FLAGS...]: mode 1, the integrals file at the third prompt
run() {
  local dir="$WORK/$1" table="$2" ec="${3:-}"
  shift 3
  (cd "$dir" && AZURE2_COULOMB_TABLE="$table" run_limited 900 "$AZURE2_BIN" calculate ${ec:+--capture-amplitudes="$ec"} --no-gui "$@" 13N.azr > run.log 2>&1)
}
SIGWARN="was computed for a different integration grid or different inputs"

TOOL="${AZURE2_COULOMB_TABLE_TOOL:-$(dirname "$AZURE2_BIN")/azure2-coulomb-table}"
make_variant first; run first "$TABLE" ""
[ -f "$WORK/first/output/intEC.dat.sig" ] || { bad "no output/intEC.dat.sig after the first run"; tail -3 "$WORK/first/run.log"; echo "$failures check(s) FAILED"; exit 1; }
mkdir -p "$WORK/elsewhere"
cp "$TABLE" "$WORK/elsewhere/copy_of_table.bin"   # no -p: another modification time as well as another path
COARSE=""
if [ -x "$TOOL" ]; then
  # A coarse table (l <= 4, 24 nodes per decade), as tests/pyazr/concurrent_sessions_test.py builds one.
  if "$TOOL" build "$WORK/coarse_table.bin" --lmax 4 --per-decade 24 --threads 1 > "$WORK/coarse_build.log" 2>&1; then
    COARSE="$WORK/coarse_table.bin"
  else
    bad "the table tool could not build the coarse table"; tail -3 "$WORK/coarse_build.log"
  fi
else
  echo "  note  no Coulomb table tool at $TOOL: case (c) skipped"
fi

# The same header with one tabulated value changed (the last bit of the file's last double: a value of the last sheet).
cp "$TABLE" "$WORK/elsewhere/one_value_changed.bin"
python3 - "$WORK/elsewhere/one_value_changed.bin" <<'PY'
import sys
with open(sys.argv[1], "r+b") as f:
    f.seek(-8, 2)
    b = bytearray(f.read(8))
    b[0] ^= 1
    f.seek(-8, 2)
    f.write(b)
PY

variants="same copy none changed"; [ -n "$COARSE" ] && variants="$variants coarse"
for v in $variants; do
  make_variant "$v"
  cp "$WORK/first/output/intEC.dat" "$WORK/first/output/intEC.dat.sig" "$WORK/$v/"
done
run same "$TABLE" intEC.dat
run copy "$WORK/elsewhere/copy_of_table.bin" intEC.dat
[ -n "$COARSE" ] && run coarse "$COARSE" intEC.dat
run changed "$WORK/elsewhere/one_value_changed.bin" intEC.dat
run none "$TABLE" intEC.dat --no-coulomb-table
warned() { grep -c "$SIGWARN" "$WORK/$1/run.log"; }
if [ "$(warned same)" -eq 0 ]; then ok "(a) the same table reuses its file"; else bad "(a) the same table refused its own file"; fi
if [ "$(warned copy)" -eq 0 ]; then ok "(b) a copy of the table at another path reuses the file"; else bad "(b) a copy of the table at another path refused the file"; fi
if [ -n "$COARSE" ]; then
  if [ "$(warned coarse)" -gt 0 ]; then ok "(c) a table with other values refuses the file and recomputes"; else bad "(c) a run with a table of other values read the file back"; fi
fi
if [ "$(warned changed)" -gt 0 ]; then ok "(f) a table with one value changed refuses the file"; else bad "(f) a table with one value changed read the file back"; fi
if [ "$(warned none)" -gt 0 ]; then ok "(d) --no-coulomb-table refuses the file and recomputes"; else bad "(d) a run without the table read back a file computed with it"; fi
# (e) a run from scratch with the copy writes the first run's signature byte for byte.
make_variant fresh_copy; run fresh_copy "$WORK/elsewhere/copy_of_table.bin" ""
if cmp -s "$WORK/first/output/intEC.dat.sig" "$WORK/fresh_copy/output/intEC.dat.sig"; then
  ok "(e) a run with the copy writes the same signature"
else
  bad "(e) a run with the copy wrote another signature"
fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
