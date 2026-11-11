#!/usr/bin/env bash
# The CLI's route selection and exit status (2026-11-01).
#
# Route.  The R-matrix route (<config> first line false) does not implement ignore-externals, which is on by default.
# Before 2026-11-01 the CLI forced the A-matrix whenever Brune OR ignore-externals was on, with "R-Matrix specified but
# --ignore-externals and --use-brune options require A-Matrix" (out of date: RMatrixFunc implements Brune), and had no
# flag to turn ignore-externals off, so it could never run the R-matrix route.  Now --no-ignore-externals turns it off;
# the R-matrix route then runs (with Brune, the default) and gives the A-matrix route's numbers.  With ignore-externals
# on, an R-matrix request still runs the A-matrix, with a warning that says why; so does one with --use-park
# (RMatrixFunc has no Park matrices; pyazr and the GUI ran it, reading Park's amplitudes as Brune's).
#
# Exit status.  Before 2026-11-01 AZURE2 exited 0 whatever happened, also after "Calculation was aborted." or a
# refused project, so a script or job chain could not tell.  Now 1, also for a format-3 project refused when it is read
# or by its set-up (a covariance file: the refusal is the last line),
# or for its data files; a project file that cannot be opened exits 255 with its own message.
#
# Checked on tests/13N (with its capture segment and target integration) and tests/7Li_p_ay --no-brune (refused: no
# standard R-matrix corresponds to its 0+ levels, tests/no_brune_transform), and a SUM data set with a term scaled by
# something that is not a number (tests/sum_lab_frame/models/sum_data.azr; also tests/pyazr/composite_terms_test.py).
#
# Usage:
#   ./tests/cli_route_exit/check.sh path/to/AZURE2
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/cli_route_exit.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

. "$HERE/../lib/edit_guard.sh"
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run NAME PROJECT ROUTE(true|false|keep) [FLAGS...]: mode 1 of tests/PROJECT under $WORK/NAME; sets $rc.
run() {
  local dir="$WORK/$1" proj="$2" route="$3"; shift 3
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/../$proj/data" "$dir/data"
  if [ "$route" = keep ]; then
    cp "$HERE/../$proj/$proj.azr" "$dir/"
  else
    awk -v r="$route" 'p { sub(/^[a-z]+/, r); p = 0 } /^<config>$/ { p = 1 } { print }' "$HERE/../$proj/$proj.azr" > "$dir/$proj.azr"
  fi
  (cd "$dir" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline "$@" "$proj.azr" > run.log 2>&1)
  rc=$?
}
total() { awk '/^Total-Chi-Squared:/ { print $2 }' "$WORK/$1/output/chiSquared.out" 2>/dev/null; }
said() { grep -q "$2" "$WORK/$1/run.log"; }
ROUTE_WARNING="WARNING: the R-matrix route (<config> first line false) does not implement ignore-externals"

run a 13N true
run r 13N false
run a_noi 13N true --no-ignore-externals
run r_noi 13N false --no-ignore-externals
run a_park 13N true --no-ignore-externals --use-park
run r_park 13N false --no-ignore-externals --use-park

# Route selection
said r "$ROUTE_WARNING" && said r "Performing A-Matrix Calculation" \
  && ok "R-matrix asked, ignore-externals on (default): A-matrix with the warning" \
  || bad "R-matrix asked, ignore-externals on: no warning or not the A-matrix"
[ "$(total r)" = "$(total a)" ] && ok "  ... chi2 $(total r) = A-matrix's" || bad "  ... chi2 '$(total r)' vs '$(total a)'"
if grep -q "require A-Matrix" "$WORK"/*/run.log; then bad "the out-of-date 'require A-Matrix' warning is still printed"
else ok "no 'options require A-Matrix' warning"; fi
said r_noi "Unknown option" && bad "--no-ignore-externals is unknown" || ok "--no-ignore-externals is an option"
! said r_noi "$ROUTE_WARNING" && said r_noi "Performing R-Matrix Calculation" \
  && ok "R-matrix asked, --no-ignore-externals: the R-matrix route runs (Brune), no warning" \
  || bad "R-matrix asked, --no-ignore-externals: not the R-matrix route"
said a_noi "Performing A-Matrix Calculation" && ok "A-matrix asked, --no-ignore-externals: A-matrix" \
  || bad "A-matrix asked, --no-ignore-externals: not the A-matrix"
t1="$(total a_noi)"; t2="$(total r_noi)"
[ -n "$t1" ] && [ "$t1" = "$t2" ] && [ "$t1" != "$(total a)" ] \
  && ok "--no-ignore-externals: R-matrix chi2 $t2 = A-matrix's (and differs from ignore-externals on, $(total a))" \
  || bad "--no-ignore-externals: chi2 R '$t2', A '$t1', A with ignore-externals '$(total a)'"
worst="$(python3 - "$WORK/a_noi/output" "$WORK/r_noi/output" <<'EOF'
import glob, os, sys
w, n = 0.0, 0
for f in sorted(glob.glob(os.path.join(sys.argv[1], "AZUREOut_*.out"))):
    g = os.path.join(sys.argv[2], os.path.basename(f))
    for la, lb in zip(open(f), open(g)):
        a, b = la.split(), lb.split()
        if len(a) < 4: continue
        x, y = float(a[3]), float(b[3]); n += 1
        if x != y: w = max(w, abs(x - y) / max(abs(x), 1e-300))
print(f"{w:.1e} {n}")
EOF
)"
set -- $worst
[ "${2:-0}" -gt 0 ] && python3 -c "import sys; sys.exit(0 if float('$1') <= 1e-9 else 1)" \
  && ok "  ... every calculated point equal (worst $1 over $2 points, 11-digit output)" \
  || bad "  ... calculated points differ: worst ${1:-?} over ${2:-0} points"

said r_park "Park's parametrization.  The A-matrix route will be used" && said r_park "Performing A-Matrix Calculation" \
  && [ -n "$(total r_park)" ] && [ "$(total r_park)" = "$(total a_park)" ] \
  && ok "R-matrix asked, --use-park: A-matrix with the warning, chi2 $(total r_park) = A-matrix's" \
  || bad "R-matrix asked, --use-park: chi2 '$(total r_park)' vs '$(total a_park)', warning/route wrong"

# Exit status
[ "$rc" = 0 ] && ok "a calculation that runs: exit 0" || bad "a calculation that runs: exit $rc"
run li7 7Li_p_ay keep --no-brune
said li7 "Calculation was aborted" && [ "$rc" = 1 ] && ok "7Li_p_ay --no-brune (aborted): exit 1" \
  || bad "7Li_p_ay --no-brune: exit $rc (aborted: $(said li7 'Calculation was aborted' && echo yes || echo no))"
dir="$WORK/badscale"; mkdir -p "$dir/output" "$dir/checks" "$dir/data"
for e in 5.95 6.05 6.15; do echo "$e  90  1.0  0.01"; done > "$dir/data/sum.dat"
sed 's/^\(    reaction       13C(a,n)16O\[6.0494\] + \)\(13C(a,n)16O\[6.12989\]\)$/\1abc*\2/' \
  "$HERE/../sum_lab_frame/models/sum_data.azr" > "$dir/sum_data.azr"
grep -q -- "+ abc\*13C" "$dir/sum_data.azr" || bad "badscale: the model edit did not apply"
(cd "$dir" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline sum_data.azr > run.log 2>&1); rc=$?
# A format-3 project is refused when it is read (Config::ReadFormat3Project): its diagnostics, then the line that says
# it is refused as its last line, exit 1 (as a refused set-up).
REFUSED="ERROR: sum_data.azr (format 3) is refused; AZURE2 --validate sum_data.azr lists the same diagnostics."
grep -q "E_REACTION_SYNTAX: dataset 'sum', field 'reaction': the scale of a term is a number" "$dir/run.log" && [ "$rc" = 1 ] \
  && [ "$(tail -n 1 "$dir/run.log")" = "$REFUSED" ] && ! grep -q "Could not open" "$dir/run.log" \
  && [ ! -s "$dir/output/chiSquared.out" ] \
  && ok "SUM term scaled by 'abc': refused with the data set named, the refusal line last, exit 1" \
  || bad "SUM term scaled by 'abc': exit $rc (want 1), last line '$(tail -n 1 "$dir/run.log")'"

# A data file of a format-3 data run that cannot be read refuses the project the same way: E_DATA_FILE, the refusal
# line (before the closing "Thanks for using AZURE2."), exit 1.
dir="$WORK/nodata"; mkdir -p "$dir/output" "$dir/checks"
cp "$HERE/../13N/13N.azr" "$dir/"
(cd "$dir" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline 13N.azr > run.log 2>&1); rc=$?
REFUSED="ERROR: 13N.azr (format 3) is refused; AZURE2 --validate 13N.azr lists the same diagnostics."
last="$(grep -v '^ *$' "$dir/run.log" | grep -v '^Thanks for using AZURE2' | tail -n 1)"
grep -q "E_DATA_FILE: dataset baumann_1658, field file: cannot read data/baumann_1658.dat" "$dir/run.log" \
  && [ "$rc" = 1 ] && [ "$last" = "$REFUSED" ] && ! grep -q "Could not fill data object" "$dir/run.log" \
  && ok "a data file that cannot be read: E_DATA_FILE, the refusal line, exit 1" \
  || bad "a data file that cannot be read: exit $rc (want 1), last line '$last'"

# A refusal of the set-up other than a data file (here a covariance file that cannot be read) ends the run the same
# way: its ERROR line is the last before "Thanks for using AZURE2.", without "Could not fill data object from file."
# after it, exit 1.
dir="$WORK/nocov"; mkdir -p "$dir/output" "$dir/checks"
cp -R "$HERE/../13N/data" "$dir/data"
edited "$HERE/../13N/13N.azr" "$dir/13N.azr" \
  awk '$1 == "dataset" { ds = $2 } { print } ds == "artemov" && $1 == "file" { print "    covariance     no_such_cov.txt" }' \
  || bad "nocov: the covariance field was not added"
(cd "$dir" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline 13N.azr > run.log 2>&1); rc=$?
last="$(grep -v '^ *$' "$dir/run.log" | grep -v '^Thanks for using AZURE2' | tail -n 1)"
case "$last" in "ERROR: Covariance file: "*no_such_cov.txt*) last_ok=1 ;; *) last_ok=0 ;; esac
[ "$rc" = 1 ] && [ "$last_ok" = 1 ] && ! grep -q "Could not fill data object" "$dir/run.log" \
  && ok "a covariance file that cannot be read: its ERROR line last, exit 1" \
  || bad "a covariance file that cannot be read: exit $rc (want 1), last line '$last'"

# A project file that does not exist keeps its own message (and exit status 255).
dir="$WORK/nofile"; mkdir -p "$dir"
(cd "$dir" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline absent.azr > run.log 2>&1); rc=$?
grep -q "Could not open absent.azr.  Check that file exists." "$dir/run.log" && [ "$rc" = 255 ] \
  && ! grep -q "is refused" "$dir/run.log" \
  && ok "a missing project file: 'Could not open', exit 255" \
  || bad "a missing project file: exit $rc, $(tail -n 1 "$dir/run.log")"

echo
if [ "$failures" -gt 0 ]; then echo "FAIL: $failures check(s)"; exit 1; fi
echo "all checks passed"
