#!/usr/bin/env bash
#
# The five Legendre expansions of quantity legendre (expansion a_over_a0 | absolute | rsl | 2l2 | 4pi; docs/upgrade/
# data_kinds.md section 4.4), on 13C(a,n)16O[6.13] (the levels of tests/angular_coefficients and a 5/2- level that
# makes the odd orders):
#   1. models/extrap.azr, the CLI: every order of every expansion equals the Gauss-Legendre projection of the same
#      run's dsigma/dOmega at the 24 nodes (truth.py cli) to 1e-8;
#   2. the same through pyazr at full precision (truth.py check): to 1e-8, and the five expansions against each other
#      by their defining relations to 1e-12 (skipped when pyazr._azure2 cannot be imported);
#   3. models/data.azr, the same rows in the five expansions (truth.py data): five equal chi^2;
#   4. that run at 1 and 2 threads: every output identical;
#   5. its analytic gradient with the 5/2- widths free (AZURE_GRAD_CHECK) equals finite differences;
#   6. --describe names each expansion's formula; refusals: rsl in b/sr (E_UNIT_NOT_ALLOWED), a fitted norm on rsl
#      (E_MODE_NOT_ALLOWED).
#
#   ./tests/legendre_expansions/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/edit_guard.sh"
PY="${PYTHON:-python3}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/legendre_expansions.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

setup() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  ln -s "$HERE/models/data" "$WORK/$1/data"
  cp "$2" "$WORK/$1/m.azr"
}
# run NAME MENU THREADS [flags]: menu 1 (calculation) or 3 (extrapolation)
run() {
  local name="$1" menu="$2" threads="$3"; shift 3
  (cd "$WORK/$name" && rm -rf output/* checks/* && printf '%s\n\n\n7\n' "$menu" |
     env OMP_NUM_THREADS="$threads" $(guard_command 300) "$AZURE2_BIN" --no-gui --no-readline "$@" m.azr > run.log 2>&1)
}

# ---- 1. the CLI's coefficients against the projection
setup ex "$HERE/models/extrap.azr"
run ex 3 "$OMP_NUM_THREADS"
if res="$("$PY" "$HERE/truth.py" cli "$WORK/ex/output" 2>&1)"; then ok "$res"; else bad "CLI coefficients: $res"; fi

# ---- 2. pyazr, full precision
if (cd /tmp && PYTHONPATH="$ROOT" "$PY" -c "import pyazr._azure2" > /dev/null 2>&1); then
  if res="$(PYTHONPATH="$ROOT" "$PY" "$HERE/truth.py" check "$WORK/ex" 2>&1 | grep -v 'WARNING\|^Coulomb')"; then
    ok "pyazr: $res"
  else
    bad "pyazr coefficients: $res"
  fi
else
  echo "  skip  pyazr (pyazr._azure2 not importable)"
fi

# ---- 3. the same rows in five expansions: five equal chi^2
setup data "$HERE/models/data.azr"
run data 1 1 && cp -r "$WORK/data/output" "$WORK/data_1"
chi="$WORK/data_1/chiSquared.out"
n="$(awk -F, 'NR > 1 && $1 ~ /^[0-9]+$/ { print $2 }' "$chi" 2>/dev/null | sort -u | wc -l)"
segs="$(awk -F, 'NR > 1 && $1 ~ /^[0-9]+$/' "$chi" 2>/dev/null | wc -l)"
if [ "$segs" = 5 ] && [ "$n" = 1 ]; then
  ok "five expansions of the same rows: five equal chi^2 ($(awk -F, 'NR == 2 { print $2 }' "$chi"))"
else
  bad "chi^2 of the five expansions: $(cat "$chi" 2>/dev/null || grep -m1 ERROR "$WORK/data/run.log")"
fi

# ---- 4. 1 and 2 threads
run data 1 2
if [ -s "$chi" ] && diff -r "$WORK/data_1" "$WORK/data/output" > /dev/null; then
  ok "data at 1 and 2 threads: every output identical"
else
  bad "data: outputs differ between 1 and 2 threads"
fi

# ---- 5. the gradient with the 5/2- widths free
edited "$HERE/models/data.azr" "$WORK/grad.azr" sed -E 's/^(2\.5  -1  11\.05  1  1  [0-9]  [0-9]  [0-9]  3  1  )1 /\10 /' ||
  bad "edit: free widths"
setup grad "$WORK/grad.azr"
(cd "$WORK/grad" && printf '2\nn\n\n\n' | AZURE_GRAD_CHECK=exit $(guard_command 300) "$AZURE2_BIN" --use-gradient --no-gui \
   --no-readline m.azr > run.log 2>&1)
line="$(grep "components with relative disagreement" "$WORK/grad/run.log" | sed -n 1p)"
worst="$(awk '/\[grad-check\] idx/ { for (i = 1; i < NF; i++) if ($i == "rel") { r = $(i + 1); if (r > w) w = r } }
              END { printf "%.2e", w + 0 }' "$WORK/grad/run.log")"
if grep -q "gamma 0/4; Jacobian path: gamma 0/4" <<< "$line" && awk -v w="$worst" 'BEGIN { exit !(w < 1e-6) }'; then
  ok "analytic gradient = finite differences over the five expansions (4 widths, largest relative $worst)"
else
  bad "gradient check: ${line:-no [grad-check] line} (largest relative $worst)"
fi

# ---- 6. --describe and refusals
desc="$(cd "$WORK/data" && "$AZURE2_BIN" --describe m.azr 2>&1)"
missing=""
for f in "a_over_a0: the coefficient of order L is A_L / A_0" "absolute: the coefficient of order L is A_L (b/sr)" \
         "rsl: the coefficient of order L is A_L / ((2L+1) A_0)" "2l2: the coefficient of order L is 2 A_L / (2L+1) (b/sr)" \
         "4pi: the coefficient of order L is 4 pi A_L (b)"; do
  grep -qF "Legendre expansion $f" <<< "$desc" || missing="$missing [$f]"
done
if [ -z "$missing" ]; then ok "--describe names the five expansions"; else bad "--describe lacks$missing"; fi
cd "$WORK/data" || exit 1
edited m.azr unit.azr awk '/^dataset leg_rsl$/ { r = 1 } r && /^    unit           1$/ { print "    unit           b/sr"; r = 0; next } { print }' ||
  bad "edit: unit"
msg="$("$AZURE2_BIN" --validate unit.azr 2>&1)"
if grep -q "E_UNIT_NOT_ALLOWED" <<< "$msg"; then ok "refused: rsl in b/sr (E_UNIT_NOT_ALLOWED)"; else bad "rsl in b/sr: $msg"; fi
edited m.azr norm.azr awk '/^dataset leg_rsl$/ { r = 1 } r && /^    norm           1 fixed$/ { print "    norm           1 fit"; r = 0; next } { print }' ||
  bad "edit: norm"
msg="$("$AZURE2_BIN" --validate norm.azr 2>&1)"
if grep -q "E_MODE_NOT_ALLOWED" <<< "$msg"; then ok "refused: a fitted norm on rsl (E_MODE_NOT_ALLOWED)"; else bad "rsl norm fit: $msg"; fi

echo
if [ "$failures" -ne 0 ]; then
  echo "legendre_expansions: $failures failure(s)"
  exit 1
fi
echo "legendre_expansions: all checks passed"
