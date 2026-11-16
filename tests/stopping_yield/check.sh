#!/usr/bin/env bash
#
# The yield of a layer whose stopping is SRIM-2013's (stopping evaluation, the engine's stopping/ library) against
# LunaSRIM's numbers (README.md of this directory; docs/upgrade/energy_grid.md section 6):
#   1. --validate prints the source of each target's stopping, with SRIM's citation;
#   2. a calculation (menu 1) of models/stopping_yield.azr at 1 and 2 threads: every output identical;
#   3. check.py: the plateau of a thick TiN layer against (lambda_r^2/2) omega gamma (m_p + m_N)/(m_N eps(E_r)) with
#      LunaSRIM's eps; the width of the step of a thinner layer against the scan with LunaSRIM's E_out; the stopping
#      as a table of the library's own values (azure2-stopping, written here) against the evaluation;
#   4. the fit (models/fit.azr: two TiN layers, N1, N2 and the one stopping scale of (TiN, p) free, its prior 1 4 %)
#      recovers N1, N2 and s within two standard errors (one thread), with one scale parameter for both layers;
#   5. the refusals: an evaluation without composition, an active element not in it, a neutral projectile, two
#      targets that share a scale and write it differently, no SRIM data, a table that does not reach the energies.
#
#   ./tests/stopping_yield/check.sh path/to/AZURE2      (azure2-stopping is taken from <build>/stopping/)

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
STOPPING_BIN="$(dirname "$AZURE2_BIN")/../stopping/azure2-stopping"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/stopping_yield.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

if [ ! -x "$STOPPING_BIN" ]; then
  echo "stopping_yield: no azure2-stopping beside the engine ($STOPPING_BIN)"
  exit 1
fi
# table FROM: the stopping table of the library's values (p in TiN, from FROM to 1.5 MeV), as models/data/tin_p.tab
table() {
  "$STOPPING_BIN" stopping --ion-z 1 --ion-mass 1.007825 --material TiN --from "$1" --to 1.5 --points 5000 |
    awk '!/^#/ && $1 + 0 == $1 { print $1, $4 }'
}
# setup NAME FROM MODEL: $WORK/NAME with a copy of the data (the table from FROM MeV) and MODEL as m.azr
setup() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -r "$HERE/models/data" "$WORK/$1/data"
  [ "$2" = 0.5 ] || table "$2" > "$WORK/$1/data/tin_p.tab"
  cp "$3" "$WORK/$1/m.azr"
}
# the committed table is the library's present values
if diff <(grep -v '^#' "$HERE/models/data/tin_p.tab") <(table 0.5) > /dev/null; then
  ok "models/data/tin_p.tab is azure2-stopping's table of the library's values"
else bad "models/data/tin_p.tab differs from azure2-stopping's output (the library changed: rewrite it, README.md)"; fi
run() {
  local name="$1" command="$2" threads="$3"
  (cd "$WORK/$name" && rm -rf output/* checks/* &&
     OMP_NUM_THREADS="$threads" $(guard_command 600) "$AZURE2_BIN" "$command" --no-gui m.azr > run.log 2>&1)
}
# refused WHAT CODE TEXT NAME: --validate, or a calculation (the set-up's refusals), of $WORK/NAME/m.azr prints CODE
refused() {
  local msg
  msg="$(cd "$WORK/$4" && "$AZURE2_BIN" --validate m.azr 2>&1; $(guard_command 120) "$AZURE2_BIN" calculate --no-gui m.azr 2>&1)"
  if grep -q "$2" <<< "$msg" && grep -q "$3" <<< "$msg"; then ok "$1 is refused ($2)"
  else bad "$1: no $2 ($(grep -m2 -E 'ERROR|error' <<< "$msg"))"; fi
}

# ---- 1. --validate names the source of the stopping
setup calc 0.5 "$HERE/models/stopping_yield.azr"
(cd "$WORK/calc" && "$AZURE2_BIN" --validate m.azr > validate.txt 2>&1)
n="$(grep -c '^Stopping of target' "$WORK/calc/validate.txt")"
c="$(grep -c 'cite SRIM: J. F. Ziegler' "$WORK/calc/validate.txt")"
if [ "$n" = 3 ] && [ "$c" = 2 ] && grep -q 'valid\.' "$WORK/calc/validate.txt"; then
  ok "--validate: the stopping of the 3 targets, SRIM cited for the 2 evaluations"
else bad "--validate: $n stopping sources (3 expected), $c SRIM citations (2): $(tail -1 "$WORK/calc/validate.txt")"; fi

# ---- 2. 1 and 2 threads
run calc calculate 1 && cp -r "$WORK/calc/output" "$WORK/calc_1"
run calc calculate 2
if [ -s "$WORK/calc_1/chiSquared.out" ] && diff -r "$WORK/calc_1" "$WORK/calc/output" > /dev/null; then
  ok "calculation at 1 and 2 threads: every output identical"
else
  bad "calculation: outputs differ between 1 and 2 threads (or none: $(grep -m1 ERROR "$WORK/calc/run.log"))"
fi
if grep -q 'cite SRIM' "$WORK/calc/run.log"; then ok "the run's log cites SRIM"; else bad "the run's log does not cite SRIM"; fi

# ---- 3. against LunaSRIM
if python3 -I "$HERE/check.py" "$HERE/models/truth.json" "$WORK/calc/output" > "$WORK/limits.txt" 2>&1; then
  ok "every check within its tolerance:"
else
  bad "a check outside its tolerance:"
fi
sed 's/^/        /' "$WORK/limits.txt"

# ---- 4. the fit
setup fit 0.5 "$HERE/models/fit.azr"
# at 1 thread: a fit at several threads is not reproducible to the last digit (Minuit's parallel gradient sums)
run fit fit 1
# (the text goes to a file: bash 3.2 of macOS reads an apostrophe of a here-document inside $( ) as a quote)
python3 -I - "$HERE/models/truth.json" "$WORK/fit/output" "$WORK/fit/run.log" > "$WORK/fit_result.txt" <<'PY'
import json, sys
truth = json.load(open(sys.argv[1]))["fit"]
out = sys.argv[2]
sav = {l.split()[0]: (float(l.split()[1]), float(l.split()[2])) for l in open(out + "/param.sav") if not l.startswith("#")}
scales = [k for k in sav if k.startswith("stopping_scale.")]
s, ds = sav["stopping_scale.TiN.1H"]
pulls = {"s": (s - truth["s_true"]) / ds}
text = [f"s = {s:.5f} +- {ds:.5f} (truth {truth['s_true']}, pull {pulls['s']:+.2f})"]
for name, d in truth["sets"].items():
    n, dn = sav[d["target"] + ".areal_density"]
    pulls[name] = (n - d["N_true"]) / dn
    text.append(f"{d['target']} N = {n:.6e} +- {dn:.1e} (truth {d['N_true']:.3e}, pull {pulls[name]:+.2f})")
chi2 = float([l.split()[1] for l in open(out + "/chiSquared.out") if l.startswith("Total-Chi-Squared:")][0])
nu = sum(len(d["E0"]) for d in truth["sets"].values()) - 3
# the objective MIGRAD minimised (the run's log) is the data's chi^2 plus the scale's prior ((s - 1)/0.04)^2
total = float([l.split(":")[1] for l in open(sys.argv[3]) if l.startswith("Total Chi-Squared:")][-1])
centre, width = truth["s_prior"]
prior = ((s - centre) / width) ** 2
prior_ok = abs(total - chi2 - prior) <= 1e-3
good = len(scales) == 1 and all(abs(p) <= 2 for p in pulls.values()) and abs(chi2 / nu - 1) <= 3 * (2.0 / nu) ** 0.5 and prior_ok
print(int(good), f"one scale parameter of {len(scales)}; " + "; ".join(text) + f"; chi2/nu = {chi2:.2f}/{nu}; "
      f"objective - chi2 = {total - chi2:.4f}, the prior {prior:.4f}")
PY
res="$(cat "$WORK/fit_result.txt")"
if [ "${res%% *}" = 1 ]; then ok "fit: ${res#* }"; else bad "fit: ${res:-no result ($(grep -m1 ERROR "$WORK/fit/run.log"))}"; fi

# ---- 5. refusals
refusal() {  # refusal NAME PYTHON-EDIT: a copy of the calculation's project with one edit
  setup "$1" "${3:-0.5}" "$HERE/models/stopping_yield.azr"
  python3 -I -c "import sys; p = sys.argv[1]; t = open(p).read(); $2; open(p, 'w').write(t)" "$WORK/$1/m.azr"
}
refusal no_composition "t = t.replace('    composition    Ti 1 N 1\n', '', 1)"
refused "an evaluation without composition" E_FIELD_REQUIRED "composition" no_composition
refusal no_active "t = t.replace('Ti 1 N 1', 'Ti 1 O 1', 1)"
refused "an active element not in the composition" E_STOPPING_ACTIVE "not in the composition" no_active
refusal neutral "t = t.replace('</data>', 'dataset nn\n    reaction       12C(n,n)12C\n    quantity       yield\n    unit           1/particle\n    file           data/plateau.dat\n    columns        E_lab[MeV] _ value error\n    norm           1 fixed\n    target         TC\n    quadrature     points 10\n</data>').replace('</targets>', 'target TC\n    kind           layer\n    active         12C\n    composition    C 1\n    areal_density  1e18 atoms/cm2 fixed\n    stopping       evaluation\n    straggling     none\n</targets>')"
refused "a neutral projectile" E_STOPPING_PROJECTILE "charged projectile" neutral
refusal conflict "t = t.replace('target TS\n    kind           layer\n    active         15N\n    composition    Ti 1 N 1\n    areal_density  5e+17 atoms/cm2 fixed\n    stopping       evaluation', 'target TS\n    kind           layer\n    active         15N\n    composition    Ti 1 N 1\n    areal_density  5e+17 atoms/cm2 fixed\n    stopping       evaluation scale 1.02 fixed')"
refused "two targets of one scale written differently" E_STOPPING_SCALE_CONFLICT "stopping_scale.TiN.1H" conflict
refusal range "pass" 0.93
refused "a table that does not reach the layer's back" E_STOPPING_RANGE "outside its range" range
mkdir -p "$WORK/empty"
msg="$(cd "$WORK/calc" && AZURE2_SRIM_DATA_DIR="$WORK/empty" "$AZURE2_BIN" --validate m.azr 2>&1; AZURE2_SRIM_DATA_DIR="$WORK/empty" $(guard_command 120) "$AZURE2_BIN" calculate --no-gui m.azr 2>&1)"
if [ "$(echo "$msg" | grep -c E_STOPPING_DATA)" -ge 2 ]; then ok "no SRIM data (AZURE2_SRIM_DATA_DIR an empty directory): refused by --validate and the run (E_STOPPING_DATA)"
else bad "no SRIM data: no E_STOPPING_DATA from both ($(grep -m2 -E 'ERROR|error' <<< "$msg"))"; fi

echo
if [ "$failures" -ne 0 ]; then
  echo "stopping_yield: $failures failure(s)"
  exit 1
fi
echo "stopping_yield: all checks passed"
