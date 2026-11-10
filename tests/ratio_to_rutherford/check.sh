#!/usr/bin/env bash
#
# Ratio to Rutherford and to Mott (docs/upgrade/data_kinds.md section 4.3; the leaf `rutherford` of step S14):
#   1. alpha + 12C, every partial wave l <= 4 pure hard-sphere scattering (models/aC_hs.azr): dsigma/dOmega and the
#      ratio to Rutherford at c.m. and at lab angles equal |f_C + f_HS|^2 and |f_C + f_HS|^2 / |f_C|^2 of truth.py
#      (mpmath, expected/aC_hs.dat) to 1e-9;
#   2. alpha + alpha, every allowed partial wave l <= 6 hard-sphere (models/aa_hs.azr): the engine's elastic
#      dsigma/dOmega is the symmetrised |f_C(theta) + f_C(pi - theta) + 2 f_HS|^2 and the ratio to Mott that over
#      |f_C(theta) + f_C(pi - theta)|^2 (expected/aa_hs.dat) to 1e-9;
#   3. the worked example of project_format3.md section 9.2 (models/aC_story.azr: lab angles, five angles in one file,
#      3 % common normalization): --describe and --print-predictions as committed (expected/story_*.txt);
#   4. the same project at 1 and 2 threads: every output file identical;
#   5. its analytic gradient (AZURE_GRAD_CHECK, with a beam spread: the Coulomb leaf folded on the sub-points) equals
#      finite differences;
#   6. ratio_to_rutherford with identical particles is refused, naming ratio_to_mott.
# The chi^2 against the same data multiplied by sigma_R, the 5 keV shift and the ratio of folds: tests/pyazr/
# ratio_to_rutherford_test.py (they need the chi^2 to more digits than chiSquared.out prints).
#
#   ./tests/ratio_to_rutherford/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ratio_to_rutherford.XXXXXX")"
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
# run NAME MENU THREADS [env ...] [-- flags]: menu 1 (calculation) or 3 (extrapolation) of $WORK/NAME/m.azr
run() {
  local name="$1" menu="$2" threads="$3"; shift 3
  (cd "$WORK/$name" && rm -rf output/* checks/* && printf '%s\n\n\n7\n' "$menu" |
     env OMP_NUM_THREADS="$threads" "$@" $(guard_command 300) "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
}

# compare EXPECTED EXTRAP: column 4 of the .extrap rows against the truth, row by row (and the grid itself)
compare() {
  awk 'NR == FNR { if ($1 !~ /^#/ && NF) { n++; e[n] = $3; t[n] = $4; v[n] = $5; k[n] = $2 } next }
       NF >= 4 { m++; if (m > n) { extra = 1; next }
                 d = $1 - e[m]; if (d < 0) d = -d; if (d > 1e-9 * e[m]) grid = 1
                 d = $3 - t[m]; if (d < 0) d = -d; if (d > 1e-9 * t[m]) grid = 1
                 r = ($4 - v[m]) / v[m]; if (r < 0) r = -r; if (r > w[k[m]]) w[k[m]] = r }
       END { if (m != n || extra || grid) { printf "rows %d of %d, grid differs %d\n", m, n, grid; exit 1 }
             for (x in w) printf "%s %.2e ", x, w[x]; printf "\n"
             for (x in w) if (w[x] > 1e-9) exit 1 }' "$1" "$2"
}

# ---- 1, 2. hard-sphere scattering against truth.py
for m in aC aa; do
  setup "hs_$m" "$HERE/models/${m}_hs.azr"
  run "hs_$m" 3 "$OMP_NUM_THREADS"
  out="$WORK/hs_$m/output/AZUREOut_aa=1_R=1.extrap"
  if [ ! -s "$out" ]; then bad "$m: no extrapolation written ($(grep -m1 ERROR "$WORK/hs_$m/run.log"))"; continue; fi
  if res="$(compare "$HERE/expected/${m}_hs.dat" "$out")"; then
    ok "$m hard sphere: dsigma/dOmega and the ratio equal truth.py to 1e-9 (largest: $res)"
  else
    bad "$m hard sphere against truth.py: $res"
  fi
done

# ---- 3. the worked example: --describe and --print-predictions
setup story "$HERE/models/aC_story.azr"
(cd "$WORK/story" && "$AZURE2_BIN" --describe m.azr 2> describe.err | sed "s|$WORK/story|<project>|g" > describe.txt)
if diff -u "$HERE/expected/story_describe.txt" "$WORK/story/describe.txt" > "$WORK/describe.diff"; then
  ok "story: --describe as committed"
else
  bad "story: --describe differs"; head -20 "$WORK/describe.diff"
fi
(cd "$WORK/story" && "$AZURE2_BIN" --print-predictions m.azr > predictions.txt 2> predictions.err)
if diff -u "$HERE/expected/story_predictions.txt" "$WORK/story/predictions.txt" > "$WORK/predictions.diff"; then
  ok "story: --print-predictions as committed ($(head -1 "$WORK/story/predictions.txt"))"
else
  bad "story: --print-predictions differs"; head -20 "$WORK/predictions.diff"
fi

# ---- 4. 1 and 2 threads
run story 1 1 && cp -r "$WORK/story/output" "$WORK/story_1"
run story 1 2
if [ -s "$WORK/story_1/chiSquared.out" ] && diff -r "$WORK/story_1" "$WORK/story/output" > /dev/null; then
  ok "story: calculation at 1 and 2 threads, every output identical ($(grep Total "$WORK/story_1/chiSquared.out"))"
else
  bad "story: outputs differ between 1 and 2 threads (or none: $(grep -m1 ERROR "$WORK/story/run.log"))"
fi

# ---- 5. the gradient, with a beam spread
sed -e 's|^    norm           1 fit prior 1 3%$|&\n    beam           B\n    quadrature     points 60|' \
    -e 's|^<targets>$|&\nbeam B\n    spread         0.05 MeV lab|' "$HERE/models/aC_story.azr" > "$WORK/grad.azr"
setup grad "$WORK/grad.azr"
(cd "$WORK/grad" && printf '2\nn\n\n\n' | AZURE_GRAD_CHECK=exit run_limited 300 "$AZURE2_BIN" --use-gradient --no-gui \
   --no-readline m.azr > run.log 2>&1)
line="$(grep "components with relative disagreement" "$WORK/grad/run.log" | head -1)"
worst="$(awk '/\[grad-check\] idx/ { for (i = 1; i < NF; i++) if ($i == "rel") { r = $(i + 1); if (r > w) w = r } }
              END { printf "%.2e", w + 0 }' "$WORK/grad/run.log")"
if echo "$line" | grep -q "energy 0/3 gamma 0/3 norm 0/1; Jacobian path: energy 0/3 gamma 0/3 norm 0/1" &&
   awk -v w="$worst" 'BEGIN { exit !(w < 1e-5) }'; then
  ok "story with a beam spread: analytic gradient = finite differences (7 parameters, largest relative $worst)"
else
  bad "gradient check: ${line:-no [grad-check] line} (largest relative $worst)"; tail -5 "$WORK/grad/run.log"
fi

# ---- 6. a refusal
sed -e 's|quantity       ratio_to_mott|quantity       ratio_to_rutherford|' "$HERE/models/aa_hs.azr" > "$WORK/aa_rth.azr"
msg="$("$AZURE2_BIN" --validate "$WORK/aa_rth.azr" 2>&1)"
if echo "$msg" | grep -q "E_REACTION_NOT_ALLOWED" && echo "$msg" | grep -q "ratio_to_mott"; then
  ok "ratio_to_rutherford of alpha + alpha refused, naming ratio_to_mott"
else
  bad "the identical-particle refusal: $msg"
fi

echo
if [ "$failures" -ne 0 ]; then
  echo "ratio_to_rutherford: $failures failure(s)"
  exit 1
fi
echo "ratio_to_rutherford: all checks passed"
