#!/usr/bin/env bash
#
# The modified S factor of heavy-ion fusion, S*(E) = E sigma(E) exp(2 pi eta + g E) (quantity modified_sfactor, the
# record field 'shape' = g; docs/upgrade/data_kinds.md section 4.1), on 15N(p,a)12C (the levels of tests/15N_p_a):
#   1. models/pa_extrap.azr: S* with g = 0.46 /MeV equals the same run's sigma times E exp(2 pi eta + g E) of truth.py
#      (mpmath, eta from the masses <levels> prints) to 1e-10, its sigma column the sigma block's, and S* with g = 0
#      the sfactor block character for character;
#   2. a fit of S* data (models/fit_sstar.azr, data/sstar.dat) and of the same rows as sigma (fit_sigma.azr,
#      data/sigma.dat): the same chi^2 as printed and every number of parameters.out within 1e-6 relative;
#   3. --describe names the modified S factor and g;
#   4. the S* fit at 1 and 2 threads: every output identical;
#   5. refusals: no 'shape' (E_FIELD_REQUIRED), a 'shape' in MeV (E_UNIT_NOT_ALLOWED), a neutral entrance
#      (E_REACTION_NOT_ALLOWED, 16O(n,g)17O of tests/ec_threshold).
#
#   ./tests/modified_sfactor/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/edit_guard.sh"
PY="${PYTHON:-python3}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/modified_sfactor.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# setup NAME MODEL: $WORK/NAME with the data and MODEL as m.azr
setup() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  ln -s "$HERE/models/data" "$WORK/$1/data"
  cp "$2" "$WORK/$1/m.azr"
}
# run NAME MENU THREADS: menu 2 (fit) or 3 (extrapolation) of $WORK/NAME/m.azr
run() {
  local name="$1" menu="$2" threads="$3" answers
  answers="$([ "$menu" = 2 ] && printf '2\nn\n\n\n' || printf '3\n\n\n7\n')"
  (cd "$WORK/$name" && rm -rf output/* checks/* && printf '%s\n' "$answers" |
     env OMP_NUM_THREADS="$threads" $(guard_command 300) "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
}

# ---- 1. the extrapolation against truth.py
setup ex "$HERE/models/pa_extrap.azr"
run ex 3 "$OMP_NUM_THREADS"
out="$WORK/ex/output/AZUREOut_aa=1_R=2.extrap"
if [ ! -s "$out" ]; then
  bad "no extrapolation written ($(grep -m1 ERROR "$WORK/ex/run.log"))"
elif res="$("$PY" "$HERE/truth.py" check "$out")"; then
  ok "extrapolation: $res"
else
  bad "extrapolation against truth.py: $res"
fi

# ---- 2. the fit of S* data against the fit of the same rows as sigma
for q in sigma sstar; do setup "fit_$q" "$HERE/models/fit_$q.azr"; run "fit_$q" 2 "$OMP_NUM_THREADS"; done
cs="$WORK/fit_sigma/output/chiSquared.out" cx="$WORK/fit_sstar/output/chiSquared.out"
ps="$WORK/fit_sigma/output/parameters.out" px="$WORK/fit_sstar/output/parameters.out"
if [ ! -s "$cs" ] || [ ! -s "$cx" ] || [ ! -s "$ps" ] || [ ! -s "$px" ]; then
  bad "fits: no output ($(grep -h -m1 ERROR "$WORK"/fit_*/run.log))"
elif ! grep -q "MIGRAD converged" "$WORK/fit_sstar/run.log"; then
  bad "the S* fit did not converge: $(grep -m1 MIGRAD "$WORK/fit_sstar/run.log")"
elif ! diff -q "$cs" "$cx" > /dev/null; then
  bad "fits: chi^2 differ: $(grep Total "$cs") / $(grep Total "$cx")"
elif worst="$(awk 'NR == FNR { n1++; for (i = 1; i <= NF; i++) a[n1, i] = $i; nf[n1] = NF; next }
                   { n2++; if (NF != nf[n2]) { bad = 1; next }
                     for (i = 1; i <= NF; i++) {
                       x = a[n2, i]; y = $i
                       if (x ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ && y ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/) {
                         d = x - y; if (d < 0) d = -d; s = x < 0 ? -x : x; if (s > 0 && d / s > w) w = d / s
                       } else if (x != y) bad = 1 } }
                   END { if (bad || n1 != n2) { print "layout differs"; exit 1 } printf "%.1e", w; exit !(w <= 1e-6) }' \
                  "$ps" "$px")"; then
  ok "fit of S* data = fit of the same rows as sigma: $(grep Total "$cx" | cut -d' ' -f1-2), parameters within $worst"
else
  bad "fits: parameters.out differ ($worst)"
fi

# ---- 3. --describe
desc="$(cd "$WORK/fit_sstar" && "$AZURE2_BIN" --describe m.azr 2>&1)"
if echo "$desc" | grep -q "modified_sfactor: angle-integrated" &&
   echo "$desc" | grep -q "the model is the modified S factor E sigma exp(2 pi eta + g E), g = 0.46 1/MeV"; then
  ok "--describe names the modified S factor and g"
else
  bad "--describe: $(echo "$desc" | head -5)"
fi

# ---- 4. 1 and 2 threads
run fit_sstar 2 1 && cp -r "$WORK/fit_sstar/output" "$WORK/sstar_1"
run fit_sstar 2 2
# fitStatus.out's "calls" line is left out: MIGRAD's call count was seen to differ between two 2-thread runs of
# the same fit (1142 against 1146) with every other byte equal, also for fit_sigma.azr and with the engine before
# modified_sfactor existed (a defect outside S*).
if [ -s "$WORK/sstar_1/chiSquared.out" ] && diff -r -x fitStatus.out "$WORK/sstar_1" "$WORK/fit_sstar/output" > /dev/null &&
   diff <(grep -v '^calls ' "$WORK/sstar_1/fitStatus.out") <(grep -v '^calls ' "$WORK/fit_sstar/output/fitStatus.out") > /dev/null; then
  ok "S* fit at 1 and 2 threads: every output identical (but MIGRAD's call count)"
else
  bad "S* fit: outputs differ between 1 and 2 threads"
fi

# ---- 5. refusals
cd "$WORK/fit_sstar" || exit 1
refused() {  # refused WHAT FILE CODE TEXT
  local msg
  msg="$("$AZURE2_BIN" --validate "$2" 2>&1)"
  if echo "$msg" | grep -q "$3" && echo "$msg" | grep -q "$4"; then ok "refused: $1 ($3)"; else bad "$1 not refused: $msg"; fi
}
edited m.azr noshape.azr sed '/^    shape          0.46 1\/MeV$/d' || bad "edit: no shape"
refused "no shape" noshape.azr E_FIELD_REQUIRED "field 'shape' is missing"
edited m.azr unit.azr sed 's|^    shape          0.46 1/MeV$|    shape          0.46 MeV|' || bad "edit: unit"
refused "shape in MeV" unit.azr E_UNIT_NOT_ALLOWED "inverse_energy"
edited "$HERE/../ec_threshold/models/ng.azr" ng.azr awk '/^<extrapolations>$/ { print; print "extrapolation ng_sstar"
  print "    reaction       16O(n,g)17O"; print "    quantity       modified_sfactor"; print "    shape          0.46 1/MeV"
  print "    unit           MeV*b"; print "    energies       range 0.1 1.0 MeV cm step 0.1"; next } { print }' ||
  bad "edit: neutral"
refused "neutral entrance 16O(n,g)17O" ng.azr E_REACTION_NOT_ALLOWED "modified_sfactor"

echo
if [ "$failures" -ne 0 ]; then
  echo "modified_sfactor: $failures failure(s)"
  exit 1
fi
echo "modified_sfactor: all checks passed"
