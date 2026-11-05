#!/usr/bin/env bash
#
# The working copies a CLI fit (and MCMC) evaluates chi2 on are made by EData::Clone.  Before 2026-11-05 Clone
# skipped the segments that hold components (SUM, RATIO, total capture, natural sums with components, a yield's
# hidden total term):
#   - their points kept the ORIGINAL data as their parent, so they read the original's target effects: a fitted
#     areal density (densityfit) -- or a fitted Doppler temperature -- never moved for them inside a fit;
#   - their hosts kept the original's mapped-point list next to the copy's, so every copy also computed and wrote
#     the original's component points (twice the component work, unsynchronised writes between threads).
#
# Model (make_models.py): tests/target_density with segment 3 (12C(p,p), 84.3 deg) made a SUM of its own term and an
# identical elastic component, the only active segment, under a thin target with a fitted density:
#
#   (a) MIGRAD (--use-gradient) recovers density_1 = 1.1e18 from 1e18 to 1e-4 with chi2 < 1e-2 (data made at 1.1e18);
#       before the fix the copies saw the SUM segment at the starting density, so its chi2 did not respond;
#   (b) the gradient check (AZURE_GRAD_CHECK): d chi2/d N from finite differences on the copies agrees with the
#       analytic and Jacobian paths.
#
#   ./tests/clone_composite/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PY="$(command -v python3 || command -v python)"
SRC="$HERE/../13N"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/clone_composite.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

"$PY" -I "$HERE/make_models.py" "$HERE/../target_density" "$WORK" || { bad "models not built"; exit 1; }

# run NAME AZR MODE FLAGS [env...] -> $WORK/NAME (mode 1 calculate, 2 fit)
run() {
  local name="$1" azr="$2" mode="$3" flags="$4"; shift 4
  local dir="$WORK/$name"
  mkdir -p "$dir/output" "$dir/checks"
  ln -s "$SRC/data" "$dir/data"
  [ -d "$WORK/fitdata" ] && ln -s "$WORK/fitdata" "$dir/fitdata"
  cp "$azr" "$dir/t.azr"
  local input='1\n\n\n7\n'
  [ "$mode" = 2 ] && input='2\nn\n\n\n'
  (cd "$dir" && for kv in "$@"; do export "$kv"; done
   printf "$input" | run_limited 900 "$AZURE2_BIN" $flags --no-gui --no-readline t.azr > run.log 2>&1)
}
out() { echo "$WORK/$1/output/AZUREOut_aa=1_R=1.out"; }
par() { awk -v n="$2" '$1 == n { print $2; exit }' "$WORK/$1/output/param.sav"; }
chi() { sed -n "s/.*$2: \([^ ]*\).*/\1/p" "$WORK/$1/output/chiSquared.out" | head -1; }
near() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a / b - 1; if (d < 0) d = -d; exit !(d < t) }'; }

# data from the truth, N = 1.1e18 (lab value = model_cm x data_lab / data_cm, as tests/target_density)
run true "$WORK/true.azr" 1 ""
if [ -s "$(out true)" ]; then
  mkdir -p "$WORK/fitdata"
  awk 'NR == FNR { if (NF >= 7) { k++; m[k] = $4; d[k] = $6 } next }
       NF >= 4 && $1 >= 0.3 && $1 <= 0.9 { j++; v = m[j] * $3 / d[j]; printf "%s %s %.12e %.12e\n", $1, $2, v, 0.01 * v }' \
      "$(out true)" "$SRC/data/meyer_84.dat" > "$WORK/fitdata/meyer_84.dat"
else
  bad "true.azr: no output"; tail -5 "$WORK/true/run.log"
fi

# (a)
run fit "$WORK/fit.azr" 2 "--use-gradient"
n="$(par fit density_1)"; c="$(chi fit Total-Chi-Squared)"
if [ -n "$n" ] && near "$n" 1.1e18 1e-4 && awk -v c="${c:-1e9}" 'BEGIN { exit !(c < 1e-2) }'; then
  ok "(a) a SUM segment under a fitted density: density_1 $n from 1e18 (truth 1.1e18), chi2 $c"
else bad "(a) a SUM segment under a fitted density: density_1 ${n:-missing}, chi2 ${c:-missing}"; tail -3 "$WORK/fit/run.log"; fi

# (b)
run grad "$WORK/fit.azr" 2 "--use-gradient" AZURE_GRAD_CHECK=exit
line="$(grep "components with relative disagreement" "$WORK/grad/run.log" | head -1)"
if echo "$line" | grep -q "density 0/1;.*density 0/1"; then ok "(b) d chi2/d N on the working copies = the analytic and Jacobian paths"
else bad "(b) gradient check: ${line:-no [grad-check] line}"; grep -m3 "density" "$WORK/grad/run.log"; fi

if [ "$failures" -eq 0 ]; then echo "PASS: clone_composite"; else echo "FAIL: $failures check(s)"; exit 1; fi
