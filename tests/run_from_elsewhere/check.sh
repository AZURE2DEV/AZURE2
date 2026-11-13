#!/usr/bin/env bash
#
# A project's relative paths are the project file's directory's, not the working directory's of the process.
#
# The data files, the output and checks directories, intEC.dat and a file named at a prompt are opened relative to the
# directory of the project file (Config::baseDirectory, Config::OpenPath); the printed names stay as the project writes
# them.  So the command-line program started elsewhere reads and writes the same files as started in the project's
# directory.  Model: tests/13N cut down to its gas-target capture segment (external capture: intEC.dat and its
# signature), as in tests/ec_table_signature.
#   (a) a calculation started in the project's directory (the reference);
#   (b) the same calculation started from the parent directory (AZURE2 proj/13N.azr) writes the same output and checks
#       files byte for byte, and nothing in the directory it was started from;
#   (c) the same from an unrelated directory with the project's absolute path;
#   (d) from the parent directory, the integrals file typed at the prompt as the project names it (output/intEC.dat)
#       is read back (no "Could not open", no signature warning).
# --validate and --describe are not tested here: they read data files relative to the working directory, as
# pyazr.validate and pyazr.describe do (tests/pyazr/validate_test.py holds the two to the same text).
#
#   ./tests/run_from_elsewhere/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/run_from_elsewhere.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# make_project DIR: the cut-down project in DIR/proj
make_project() {
  local dir="$WORK/$1/proj"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$SRC/data" "$dir/data"
  awk '{ sub(/\r$/, "") } /^dataset / { n++; if (n > 1) { print; print "    active         no"; next } } { print }' \
    "$SRC/13N.azr" > "$dir/13N.azr"
}
# calc START-DIR PROJECT-ARGUMENT [INTEGRALS-FILE]: mode 1, the integrals file at the third prompt
calc() {
  (cd "$1" && printf '1\n\n%s\n7\n' "${3:-}" | run_limited 900 "$AZURE2_BIN" --no-gui --no-readline "$2" > run.log 2>&1)
}
same_files() {  # same_files DIR-A DIR-B: the output and checks files of two projects equal byte for byte
  diff -r "$1/output" "$2/output" > /dev/null && diff -r "$1/checks" "$2/checks" > /dev/null
}

make_project a; calc "$WORK/a/proj" 13N.azr
[ -s "$WORK/a/proj/output/intEC.dat" ] || { bad "(a) the reference run wrote no output/intEC.dat"; tail -5 "$WORK/a/proj/run.log"; echo "$failures check(s) FAILED"; exit 1; }
ok "(a) the reference run in the project's directory ($(ls "$WORK/a/proj/output" | wc -l) output files)"

make_project b; calc "$WORK/b" proj/13N.azr
if same_files "$WORK/a/proj" "$WORK/b/proj"; then ok "(b) started from the parent directory: the same files"; else bad "(b) started from the parent directory: the files differ ($(diff -rq "$WORK/a/proj/output" "$WORK/b/proj/output" 2>&1 | head -3 | tr '\n' ' '))"; fi
stray="$(ls "$WORK/b" | grep -v -x -e proj -e run.log)"
if [ -z "$stray" ]; then ok "(b) nothing written where it was started"; else bad "(b) written where it was started: $stray"; fi

make_project c; mkdir -p "$WORK/elsewhere"; calc "$WORK/elsewhere" "$WORK/c/proj/13N.azr"
if same_files "$WORK/a/proj" "$WORK/c/proj"; then ok "(c) started from an unrelated directory: the same files"; else bad "(c) started from an unrelated directory: the files differ"; fi

make_project d; cp "$WORK/a/proj/output/intEC.dat" "$WORK/a/proj/output/intEC.dat.sig" "$WORK/d/proj/output/"
calc "$WORK/d" proj/13N.azr output/intEC.dat
if grep -q "Could not open external capture file\|Cannot Read From\|was computed for a different" "$WORK/d/run.log"; then
  bad "(d) the integrals file named at the prompt was not read back"
elif same_files "$WORK/a/proj" "$WORK/d/proj"; then
  ok "(d) the integrals file named at the prompt is the project's (read back, the same files)"
else
  bad "(d) the run that read the integrals back wrote other files"
fi

if [ "$failures" -ne 0 ]; then echo "$failures check(s) FAILED"; exit 1; fi
echo "all checks passed"
