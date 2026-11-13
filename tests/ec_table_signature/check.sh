#!/usr/bin/env bash
#
# The Coulomb table a run reads is part of the external-capture integral signature written beside intEC.dat.
#
# The integrals take the Coulomb waves at their radial nodes from the table when one is used (CoulFunc), so a file
# computed with one table (or with the exact functions) holds other digits than one computed with another.  Until
# branch feature/upgrade-6-evaluation-3 the signature held neither whether a table was used nor which file, so a run
# handed the file of a run with another table, or without one, read it back silently.
#
# Model: tests/13N cut down to its gas-target capture segment, as in tests/ec_relativistic_signature.  The table is
# $AZURE2_COULOMB_TABLE (required); a copy of it at another path is "another table" (the signature names the file by
# its path, size and modification time, so a rebuilt table at the same path is another table too).
#   (a) a run with the table given the intEC.dat of a run with the same table reuses it without a warning;
#   (b) a run with the copy at another path given that file refuses it (signature warning) and recomputes;
#   (c) a run with --no-coulomb-table given that file refuses it and recomputes.
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
  (cd "$dir" && printf '1\n\n%s\n7\n' "$ec" |
    AZURE2_COULOMB_TABLE="$table" run_limited 900 "$AZURE2_BIN" --no-gui --no-readline "$@" 13N.azr > run.log 2>&1)
}
SIGWARN="was computed for a different integration grid or different inputs"

make_variant first; run first "$TABLE" ""
[ -f "$WORK/first/output/intEC.dat.sig" ] || { bad "no output/intEC.dat.sig after the first run"; tail -3 "$WORK/first/run.log"; echo "$failures check(s) FAILED"; exit 1; }
cp -p "$TABLE" "$WORK/other_table.bin"

for v in same other none; do
  make_variant "$v"
  cp "$WORK/first/output/intEC.dat" "$WORK/first/output/intEC.dat.sig" "$WORK/$v/"
done
run same "$TABLE" intEC.dat
run other "$WORK/other_table.bin" intEC.dat
run none "$TABLE" intEC.dat --no-coulomb-table
warned() { grep -c "$SIGWARN" "$WORK/$1/run.log"; }
if [ "$(warned same)" -eq 0 ]; then ok "(a) the same table reuses its file"; else bad "(a) the same table refused its own file"; fi
if [ "$(warned other)" -gt 0 ]; then ok "(b) another table file refuses the file and recomputes"; else bad "(b) a run with another table file read the file back"; fi
if [ "$(warned none)" -gt 0 ]; then ok "(c) --no-coulomb-table refuses the file and recomputes"; else bad "(c) a run without the table read back a file computed with it"; fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
