#!/usr/bin/env bash
#
# The --ignore-externals mask and a parameter file (docs/upgrade/external_capture.md section 3).
#
# The mask is a property of the input: a photon channel whose width is given as exactly 0 has no gamma decay
# and no external capture.  A parameter file (--parameter-file) is input too: a photon amplitude it gives as
# exactly 0 is masked and fixed whatever <levels> says, and one it gives as non-zero is not masked.  Before
# 2026-11-25 only <levels> decided, so a project whose <levels> holds a photon width that the fit's param.sav
# holds as 0 kept that channel's external capture (12C(p,g) angular distributions moved by 1.25e-3 against
# the engine before the restructuring, which tested the value in use).
#
# Model: tests/13N (12C(p,g)), the E1 width of the 3/2- level at 3.5032 MeV (width_5_2).
#
#   ./tests/ec_mask_parameter_file/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../13N"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ec_mask_pf.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME PROJECT [PARAMETER_FILE]
run() {
  local dir="$WORK/$1" option=""
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$SRC/data" "$dir/data"
  cp "$2" "$dir/t.azr"
  if [ -n "${3:-}" ]; then cp "$3" "$dir/start.sav"; option="--parameter-file=start.sav"; fi
  (cd "$dir" && run_limited 600 "$AZURE2_BIN" calculate $option t.azr > run.log 2>&1)
  if [ ! -s "$dir/output/AZUREOut_aa=1_R=2.out" ]; then
    bad "$1: no output written"; tail -15 "$dir/run.log"; return 1
  fi
}
# maxrel A B: largest relative difference of the calculated capture cross sections (column 4), rows in order
maxrel() {
  paste "$WORK/$1/output/AZUREOut_aa=1_R=2.out" "$WORK/$2/output/AZUREOut_aa=1_R=2.out" | awk '
    NF >= 8 { n = NF / 2; a = $4; b = $(n + 4); d = a - b; if (d < 0) d = -d; m = b < 0 ? -b : b; if (m < 1e-300) m = 1e-300
              r = d / m; if (r > worst) worst = r; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
within() { awk -v r="$1" -v t="$2" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && r != "nan" && r + 0 <= t) }'; }
noted() { awk -v want="$2" 'index($0, "Note:") == 1 && index($0, want) { found = 1 } END { exit !found }' "$WORK/$1/run.log"; }

# The project, and the project with the E1 width of the 3.5032 MeV level given as 0.
cp "$SRC/13N.azr" "$WORK/with.azr"
awk '{ sub(/\r$/, "") }
     /^<levels>$/ { inlevels = 1 } /^<\/levels>$/ { inlevels = 0 }
     inlevels && $3 == "3.5032" && $12 == "-0.49" { sub(/-0\.49/, "    0"); changed++ }
     { print } END { exit changed != 1 }' "$SRC/13N.azr" > "$WORK/without.azr" || bad "the E1 width line of the 3.5032 MeV level was not found once"

# The project's own parameters (param.par of a run without a parameter file), and the same with width_5_2 = 0.
run plain "$WORK/with.azr" || { echo "FAILED"; exit 1; }
cp "$WORK/plain/output/param.par" "$WORK/full.sav"
awk '$1 == "width_5_2" { printf "%20s %19s %19s\n", $1, "0.0000000e+00", "0.0000000e+00"; changed++; next }
     { print } END { exit changed != 1 }' "$WORK/full.sav" > "$WORK/zero.sav" || bad "width_5_2 was not found once in param.par"
if awk '$1 == "width_5_2" && $2 + 0 != 0 { found = 1 } END { exit !found }' "$WORK/full.sav"; then ok "param.par holds a non-zero width_5_2"
else bad "param.par holds no non-zero width_5_2"; fi

run with_zero "$WORK/with.azr" "$WORK/zero.sav"
run without_zero "$WORK/without.azr" "$WORK/zero.sav"
run with_full "$WORK/with.azr" "$WORK/full.sav"
run without_full "$WORK/without.azr" "$WORK/full.sav"
run without_plain "$WORK/without.azr"

# 1. A 0 in the parameter file is a 0 in <levels>.
r="$(maxrel with_zero without_zero)"
if within "$r" 1e-12; then ok "0 in the parameter file = 0 in <levels> (max rel diff $r)"
else bad "0 in the parameter file differs from 0 in <levels> by $r (limit 1e-12)"; fi
r="$(maxrel with_zero without_plain)"
if within "$r" 1e-6; then ok "and = the project with the width given as 0, without a file (max rel diff $r)"
else bad "0 in the parameter file differs from the project with the width 0 by $r (limit 1e-6)"; fi
if noted with_zero "exactly 0 in the parameter"; then ok "the run says that the file's 0 switched the external capture off"
else bad "no note about the file's 0"; fi

# 2. A non-zero value in the parameter file is not masked, whatever <levels> says.
r="$(maxrel without_full with_full)"
if within "$r" 1e-12; then ok "non-zero in the parameter file, 0 in <levels> = non-zero in both (max rel diff $r)"
else bad "non-zero in the parameter file with 0 in <levels> differs from non-zero in both by $r (limit 1e-12)"; fi
if noted without_full "non-zero in the parameter file"; then ok "the run says that the file's value switched the external capture on"
else bad "no note about the file's non-zero value"; fi
r="$(maxrel with_full plain)"
if within "$r" 1e-6; then ok "the project's own param.par reproduces the project, to its 8 digits (max rel diff $r)"
else bad "the project's own param.par differs from the project by $r (limit 1e-6)"; fi
if noted with_full "photon width"; then bad "a note although file and <levels> agree"; else ok "no note when file and <levels> agree"; fi

# 3. The channel matters (else the comparisons above prove nothing).
r="$(maxrel with_zero with_full)"
if [ "$r" = nan ] || within "$r" 1e-3; then bad "the width's 0 changed the capture cross sections by only $r"
else ok "the width's 0 moves the capture cross sections (max rel change $r)"; fi

if [ "$failures" -eq 0 ]; then echo "PASSED"; else echo "FAILED: $failures"; fi
exit "$failures"
