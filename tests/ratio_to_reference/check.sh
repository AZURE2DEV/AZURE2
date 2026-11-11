#!/usr/bin/env bash
#
# A cross section measured relative to a tabulated standard (quantity ratio_to_reference, field 'reference <file>
# [unit <u>]'; docs/upgrade/data_kinds.md section 4.9), on 15N(p,a)12C (the levels of tests/15N_p_a):
#   1. models/extrap.azr: the ratio to data/ref_cm.dat (E_cm, mb) equals the same run's sigma over the table as
#      truth.py interpolates it (linear in ln sigma against E_cm, mpmath) to 1e-10, and the ratio to a constant 1 b
#      equals sigma character for character;
#   2. a fit of the ratio data (fit_ratio.azr, data/ratio.dat = sigma.dat / table) and of sigma.dat (fit_sigma.azr):
#      the same chi^2 as printed and every number of parameters.out within 1e-6 relative;
#   3. folded.azr: sigma.dat as sigma and as the ratio to 1 b through a 10 keV beam spread (the reference folded on
#      the point's own sub-points): the model columns equal character for character;
#   4. fit_ratio.azr with that beam spread: the analytic gradient (AZURE_GRAD_CHECK) equals finite differences;
#   5. --describe names the table; 1 against 2 threads identical (calculation of folded.azr);
#   6. refusals: no 'reference' (E_FIELD_REQUIRED), a grid point outside the table, a table value <= 0 and a data
#      row outside the table (E_REFERENCE_TABLE).
#
#   ./tests/ratio_to_reference/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/edit_guard.sh"
PY="${PYTHON:-python3}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ratio_to_reference.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

setup() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  ln -s "$HERE/models/data" "$WORK/$1/data"
  cp "$2" "$WORK/$1/m.azr"
}
# run NAME MENU THREADS [env ...] [-- flags]: menu 1, 2 (fit, no MINOS) or 3
run() {
  local name="$1" menu="$2" threads="$3" answers; shift 3
  answers="$([ "$menu" = 2 ] && printf '2\nn\n\n\n' || printf '%s\n\n\n7\n' "$menu")"
  (cd "$WORK/$name" && rm -rf output/* checks/* && printf '%s\n' "$answers" |
     env OMP_NUM_THREADS="$threads" "$@" $(guard_command 300) "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
}
# block N FILE: columns 1-4 of the N-th block (0-based) of an AZUREOut file
block() { awk -v b="$1" 'NF == 0 { if (seen) n++; seen = 0; next } { seen = 1; if (n == b) print $1, $2, $3, $4 }' "$2"; }

# ---- 1. the extrapolation against truth.py
setup ex "$HERE/models/extrap.azr"
run ex 3 "$OMP_NUM_THREADS"
out="$WORK/ex/output/AZUREOut_aa=1_R=2.extrap"
if [ ! -s "$out" ]; then bad "no extrapolation written ($(grep -m1 ERROR "$WORK/ex/run.log"))"
elif res="$("$PY" "$HERE/truth.py" check "$out")"; then ok "extrapolation: $res"
else bad "extrapolation against truth.py: $res"; fi

# ---- 2. the fit of ratio data against the fit of sigma data
for q in sigma ratio; do setup "fit_$q" "$HERE/models/fit_$q.azr"; run "fit_$q" 2 "$OMP_NUM_THREADS"; done
cs="$WORK/fit_sigma/output/chiSquared.out" cx="$WORK/fit_ratio/output/chiSquared.out"
ps="$WORK/fit_sigma/output/parameters.out" px="$WORK/fit_ratio/output/parameters.out"
if [ ! -s "$cs" ] || [ ! -s "$cx" ] || [ ! -s "$ps" ] || [ ! -s "$px" ]; then
  bad "fits: no output ($(grep -h -m1 ERROR "$WORK"/fit_*/run.log))"
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
  ok "fit of ratio data = fit of the same rows as sigma: $(grep Total "$cx" | cut -d' ' -f1-2), parameters within $worst"
else
  bad "fits: parameters.out differ ($worst)"
fi

# ---- 3. folded: the ratio to 1 b equals sigma through the same beam spread
setup fold "$HERE/models/folded.azr"
run fold 1 1 && cp -r "$WORK/fold/output" "$WORK/fold_1"
f="$WORK/fold_1/AZUREOut_aa=1_R=2.out"
if [ -s "$f" ] && [ "$(block 0 "$f" | wc -l)" = 16 ] && diff <(block 0 "$f") <(block 1 "$f") > /dev/null; then
  ok "folded through a 10 keV beam spread: the ratio to 1 b equals sigma character for character (16 points)"
else
  bad "folded: the ratio to 1 b differs from sigma ($(grep -m1 ERROR "$WORK/fold/run.log"))"
fi

# ---- 4. the gradient with the beam spread
edited "$HERE/models/fit_ratio.azr" "$WORK/grad.azr" sed -e 's|^    norm           1 fixed$|&\n    energy_window  0.3 1.6 MeV\n    beam           B\n    quadrature     points 30|' \
  -e 's|^<targets>$|&\nbeam B\n    spread         0.01 MeV lab|' || bad "edit: beam"
setup grad "$WORK/grad.azr"
(cd "$WORK/grad" && printf '2\nn\n\n\n' | AZURE_GRAD_CHECK=exit $(guard_command 300) "$AZURE2_BIN" --use-gradient --no-gui \
   --no-readline m.azr > run.log 2>&1)
line="$(grep "components with relative disagreement" "$WORK/grad/run.log" | head -1)"
worst="$(awk '/\[grad-check\] idx/ { for (i = 1; i < NF; i++) if ($i == "rel") { r = $(i + 1); if (r > w) w = r } }
              END { printf "%.2e", w + 0 }' "$WORK/grad/run.log")"
if echo "$line" | grep -q "energy 0/4 gamma 0/10; Jacobian path: energy 0/4 gamma 0/10" && awk -v w="$worst" 'BEGIN { exit !(w < 1e-5) }'; then
  ok "folded ratio: analytic gradient = finite differences (14 parameters, largest relative $worst)"
else
  bad "gradient check: ${line:-no [grad-check] line} (largest relative $worst)"
fi

# ---- 5. --describe, threads
desc="$(cd "$WORK/fit_ratio" && "$AZURE2_BIN" --describe m.azr 2>&1)"
if echo "$desc" | grep -q "divided by the reference cross section of data/ref_cm.dat (37 rows"; then
  ok "--describe names the table"
else
  bad "--describe: $(echo "$desc" | head -4)"
fi
run fold 1 2
if diff -r "$WORK/fold_1" "$WORK/fold/output" > /dev/null; then ok "folded.azr at 1 and 2 threads: every output identical"
else bad "folded.azr: outputs differ between 1 and 2 threads"; fi

# ---- 6. refusals
cd "$WORK/fit_ratio" || exit 1
refused() {  # refused WHAT FILE CODE TEXT
  local msg
  msg="$("$AZURE2_BIN" --validate "$2" 2>&1)"
  if echo "$msg" | grep -q "$3" && echo "$msg" | grep -q "$4"; then ok "refused: $1 ($3)"; else bad "$1 not refused: $msg"; fi
}
edited m.azr noref.azr sed '/^    reference      data\/ref_cm.dat unit mb$/d' || bad "edit: no reference"
refused "no reference" noref.azr E_FIELD_REQUIRED "field 'reference' is missing"
edited "$HERE/models/extrap.azr" far.azr sed 's|range 0.1 1.6 MeV cm step 0.05|range 0.1 2.5 MeV cm step 0.05|' || bad "edit: far grid"
ln -s "$HERE/models/data" "$WORK/fit_ratio/data" 2> /dev/null
refused "a grid point outside the table" far.azr E_REFERENCE_TABLE "is outside the reference table"
mkdir -p bad
edited "$HERE/models/data/ref_cm.dat" bad/ref0.dat awk 'NR == 5 { $2 = "0" } { print }' || bad "edit: table"
edited m.azr zero.azr sed 's|data/ref_cm.dat unit mb|bad/ref0.dat unit mb|' || bad "edit: zero"
refused "a table value 0" zero.azr E_REFERENCE_TABLE "is not positive"
edited "$HERE/models/data/ratio.dat" bad/far.dat awk '{ print } END { print "1.95  0  1.0  0.1" }' || bad "edit: row"
edited m.azr rowfar.azr sed 's|^    file           data/ratio.dat$|    file           bad/far.dat|' || bad "edit: rowfar"
refused "a data row outside the table" rowfar.azr E_REFERENCE_TABLE "the row at E_lab = 2.08098 MeV is outside"

echo
if [ "$failures" -ne 0 ]; then
  echo "ratio_to_reference: $failures failure(s)"
  exit 1
fi
echo "ratio_to_reference: all checks passed"
