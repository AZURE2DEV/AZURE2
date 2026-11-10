#!/usr/bin/env bash
#
# The four data check files (lines 6-9 of <config>: data, energy-dependent values, Legendre polynomials, Coulomb
# amplitudes; EData::PrintData, PrintEDependentValues, PrintLegendreP, PrintCoulombAmplitude in
# src/EDataCheckFiles.cpp) are written.  No other test project sets these checks, so without this test the suite
# does not reach the four writers.  Model: tests/upgrade/T1 calculated with data (mode 1) with the four checks set
# to "file"; each of checks/data.chk, lomatrixandpene.chk, legendre.chk, coulombamplitudes.chk must exist and have
# the line count recorded below (counted with the build of feature/upgrade-4-data at dcc1fb0a).  A changed count
# means the writer or the points it walks changed: find out which before updating the number.
#
#   ./tests/check_files_data/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/check_files_data.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

cp -R "$HERE/../upgrade/T1/." "$WORK/"
rm -rf "$WORK/output" "$WORK/checks"
mkdir -p "$WORK/output" "$WORK/checks"
# lines 6-9 of <config> (the data checks): none -> file
awk '/^<config>/ { inside = 1; n = 0; print; next }
     /^<\/config>/ { inside = 0 }
     inside { n++; if (n >= 6 && n <= 9) sub(/^none/, "file") }
     { print }' "$HERE/../upgrade/T1/T1.azr" > "$WORK/T1.azr"
(cd "$WORK" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline T1.azr > run.log 2>&1)
[ -f "$WORK/output/chiSquared.out" ] || { bad "the calculation wrote no chiSquared.out"; tail -20 "$WORK/run.log"; }

# file and its line count
for entry in data.chk:472 lomatrixandpene.chk:4425 legendre.chk:9287 coulombamplitudes.chk:319; do
  f=${entry%%:*}; want=${entry##*:}
  if [ ! -f "$WORK/checks/$f" ]; then bad "checks/$f not written"; continue; fi
  got=$(wc -l < "$WORK/checks/$f" | tr -d ' ')
  if [ "$got" = "$want" ]; then ok "checks/$f: $got lines"; else bad "checks/$f: $got lines, expected $want"; fi
done

if [ "$failures" -eq 0 ]; then echo "PASS: the four data check files of T1"; exit 0; fi
echo "FAIL: $failures check(s)"; exit 1
