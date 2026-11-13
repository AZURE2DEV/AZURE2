#!/usr/bin/env bash
#
# The yield of a layer of finite thickness (quantity yield on a target of kind layer; docs/upgrade/energy_grid.md
# section 5, README.md of this directory):
#   1. --describe names the layer kernel on every yield on a layer;
#   2. a calculation (menu 1) of models/layer_yield.azr at 1 and 2 threads: every output identical;
#   3. its model values against the closed forms check.py computes from models/truth.json (constant sigma, an isotope
#      abundance and a beam spread; thin layer against N sigma(E0 - Delta/2) and against the thin-target sigma; thick
#      layer against the infinite target; the narrow resonance against the arctan form; per uC for charge states 1 and
#      2), each to its stated tolerance;
#   3b. the same with every sub-point grid refined twofold (AZURE_REFINE_QUADRATURE=2): the unfolded data sets do not
#      move, the folded ones move by less than 1/p^2, the thick layer comes closer to its independent integral;
#   4. the fit (models/fit.azr: a scan of the resonance with 2 % noise, N and Gamma_p free) recovers N and Gamma_p
#      within two standard errors, with chi^2/nu near 1, and gives the same output at 1 and 2 threads;
#   5. the refusals that remain: a layer without an areal density, a yield per uC without charge_state, the deleted
#      fold_windows field, a stopping not positive in the layer, a layer the beam does not leave;
#   6. pyazr (pyazr_layer_test.py: the records built in Python, the arctan form, dY/dN against 1e-24 sigma(E_out))
#      when the compiled module can be imported (PYAZR_ENGINE_DIR), else a skip.
#
#   ./tests/layer_yield/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/layer_yield.XXXXXX")"
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
# run NAME MENU THREADS: menu 1 (calculation) or 2 (fit) of $WORK/NAME/m.azr
run() {
  local name="$1" menu="$2" threads="$3" input='1\n\n\n7\n'
  [ "$menu" = 2 ] && input='2\nn\n\n\n'
  (cd "$WORK/$name" && rm -rf output/* checks/* && printf "$input" |
     OMP_NUM_THREADS="$threads" $(guard_command 600) "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
}
# refused NAME CODE TEXT FILE: --describe of FILE prints CODE and TEXT
refused() {
  local msg
  # the validator's and Build's refusals by --describe (data files relative to the working directory); the set-up's by
  # a calculation (menu 1) of a copy of FILE beside the data (a run reads them relative to the project file)
  cp "$WORK/$4" "$WORK/calc/$4"
  msg="$(cd "$WORK/calc" && "$AZURE2_BIN" --describe "$WORK/$4" 2>&1; printf '1\n\n\n7\n' |
         $(guard_command 120) "$AZURE2_BIN" --no-gui --no-readline "$4" 2>&1)"
  if echo "$msg" | grep -q "$2" && echo "$msg" | grep -q "$3"; then ok "$1 is refused ($2)"
  else bad "$1: no $2 ($(echo "$msg" | grep -m2 -E 'ERROR|error'))"; fi
}

# ---- 1. --describe
setup calc "$HERE/models/layer_yield.azr"
(cd "$WORK/calc" && "$AZURE2_BIN" --describe m.azr > describe.txt 2> describe.err)
n="$(grep -c 'yield per incident particle on the layer' "$WORK/calc/describe.txt")"
if [ "$n" = 9 ]; then ok "--describe: the layer kernel on each of the 9 yields on a layer"
else bad "--describe names the layer kernel on $n data sets (9 expected): $(head -3 "$WORK/calc/describe.err")"; fi

# ---- 2. 1 and 2 threads
run calc 1 1 && cp -r "$WORK/calc/output" "$WORK/calc_1"
run calc 1 2
if [ -s "$WORK/calc_1/chiSquared.out" ] && diff -r "$WORK/calc_1" "$WORK/calc/output" > /dev/null; then
  ok "calculation at 1 and 2 threads: every output identical"
else
  bad "calculation: outputs differ between 1 and 2 threads (or none: $(grep -m1 ERROR "$WORK/calc/run.log"))"
fi

# ---- 3. the limits
if python3 -I "$HERE/check.py" "$HERE/models/truth.json" "$WORK/calc/output" > "$WORK/limits.txt" 2>&1; then
  ok "every limit within its tolerance:"
else
  bad "a limit outside its tolerance:"
fi
sed 's/^/        /' "$WORK/limits.txt"

# ---- 3b. the refined quadrature (AZURE_REFINE_QUADRATURE=2: every sub-point grid twice as fine, the diagnostic of the
#      quadrature gate, docs/upgrade/energy_grid.md section 7): the data sets without a target (thin_bare, deriv_bare,
#      res_bare) do not move; every folded yield moves by less than the grid's own error estimate, 1/p^2 of the largest
#      value of its data set (p = 50 points per resonance width: the linear interpolation over a width, section 2.4;
#      measured 1.2e-4, thick_inf), and at least one moves;
#      the thick layer comes closer to the independent integral (models/thick_reference.json).
setup fine "$HERE/models/layer_yield.azr"
(export AZURE_REFINE_QUADRATURE=2; run fine 1 2)
res="$(python3 -I - "$WORK/calc/output" "$WORK/fine/output" "$HERE/models/thick_reference.json" <<'PY'
import json, sys
def blocks(path):
    out, cur = [], []
    for line in open(path):
        if line.strip():
            cur.append([float(x) for x in line.split()])
        elif cur:
            out.append(cur); cur = []
    if cur:
        out.append(cur)
    return out
names = ["thin", "thin_sigma", "thin_bare", "thick", "thick_inf", "res", "res_q1", "res_q2", "deriv", "deriv_bare", "res_bare"]
unfolded = {"thin_bare", "deriv_bare", "res_bare"}
a = dict(zip(names, blocks(sys.argv[1] + "/AZUREOut_aa=1_R=2.out")))
b = dict(zip(names, blocks(sys.argv[2] + "/AZUREOut_aa=1_R=2.out")))
na, nb = blocks(sys.argv[1] + "/AZUREOut_aa=3_R=3.out"), blocks(sys.argv[2] + "/AZUREOut_aa=3_R=3.out")
for k, (x, y) in enumerate(zip(na, nb)):
    a["n%d" % k], b["n%d" % k] = x, y
bound = 1.0 / 50 ** 2
moved_unfolded, worst, moved = [], 0.0, 0
for name in a:
    scale = max(abs(r[3]) for r in a[name])  # a yield relative to the largest of its data set (a plateau, a peak)
    for ra, rb in zip(a[name], b[name]):
        d = abs(rb[3] - ra[3]) / scale
        if name in unfolded:
            if ra[3] != rb[3]:
                moved_unfolded.append(name)
        else:
            worst = max(worst, d)
            moved += d > 0
ref = json.load(open(sys.argv[3]))
closer = all(abs(yb / r - 1) < abs(ya / r - 1)
             for ya, yb, r in zip([x[3] for x in a["thick"]], [x[3] for x in b["thick"]], ref["Y_reference"]))
dev = max(abs(x[3] / r - 1) for x, r in zip(b["thick"], ref["Y_reference"]))
good = not moved_unfolded and moved > 0 and worst <= bound and closer
print(int(good), f"unfolded moved: {sorted(set(moved_unfolded)) or 'none'}; {moved} folded values moved, the largest by "
      f"{worst:.2e} of its data set's largest value (bound {bound:.0e}); thick layer {dev:.2e} from the independent integral ({'closer' if closer else 'NOT closer'})")
PY
)"
if [ "${res%% *}" = 1 ]; then ok "refined quadrature: ${res#* }"; else bad "refined quadrature: ${res#* }"; fi

# ---- 4. the fit
setup fit "$HERE/models/fit.azr"
run fit 2 1 && cp -r "$WORK/fit/output" "$WORK/fit_1"
run fit 2 2
# fitStatus.out's call count is Minuit's unsynchronised counter (docs/upgrade/fit_threads.md section 3): left out
if [ -s "$WORK/fit_1/param.sav" ] && diff -r -x fitStatus.out "$WORK/fit_1" "$WORK/fit/output" > /dev/null &&
   diff <(grep -v '^calls' "$WORK/fit_1/fitStatus.out") <(grep -v '^calls' "$WORK/fit/output/fitStatus.out") > /dev/null; then
  ok "fit at 1 and 2 threads: every output identical (but Minuit's call count)"
else
  bad "fit: outputs differ between 1 and 2 threads (or none: $(grep -m1 ERROR "$WORK/fit/run.log"))"
fi
res="$(python3 -I - "$HERE/models/truth.json" "$WORK/fit/output" <<'PY'
import json, sys
truth = json.load(open(sys.argv[1]))["fit"]
out = sys.argv[2]
sav = {l.split()[0]: (float(l.split()[1]), float(l.split()[2])) for l in open(out + "/param.sav") if not l.startswith("#")}
N, dN = sav["TF.areal_density"]
g, dg = sav["width_1_1"]
gp = None
for line in open(out + "/parameters.out"):
    if "R =  1" in line and "G  =" in line and gp is None:
        gp = float(line.split("G  =")[1].split()[0])
dgp = 2 * gp * dg / g          # Gamma_p = 2 P g^2 at the resonance energy, which does not move
chi2 = float([l.split()[1] for l in open(out + "/chiSquared.out") if l.startswith("Total-Chi-Squared:")][0])
nu = len(truth["E0"]) - 2
pN, pG = (N - truth["N_true"]) / dN, (gp - truth["Gp_true"]) / dgp
good = abs(pN) <= 2 and abs(pG) <= 2 and abs(chi2 / nu - 1) <= 3 * (2.0 / nu) ** 0.5
print(int(good), f"N = {N:.7e} +- {dN:.2e} (truth {truth['N_true']:.3e}, pull {pN:+.2f}); Gamma_p = {gp:.4f} +- "
      f"{dgp:.4f} eV (truth {truth['Gp_true']}, pull {pG:+.2f}); chi2/nu = {chi2:.2f}/{nu}")
PY
)"
if [ "${res%% *}" = 1 ]; then ok "fit: ${res#* }"; else bad "fit: ${res:-no result ($(grep -m1 ERROR "$WORK/fit/run.log"))}"; fi

# ---- 5. refusals
python3 -I - "$HERE/models/layer_yield.azr" "$WORK" <<'PY'
import re, sys
src, work = sys.argv[1:3]
text = open(src).read()
def target(name):
    return re.search(r"target %s\n(    .*\n)+" % name, text).group(0)
def write(name, new):
    open(f"{work}/{name}.azr", "w").write(new)
t = target("TR")
write("no_density", text.replace(t, t.replace("    areal_density  1e+18 atoms/cm2 fixed\n", "")))
write("no_charge", text.replace("    charge_state   1\n", "", 1))
write("windows", text.replace("dataset res\n", "dataset res\n    fold_windows   0.95-0.97 MeV\n", 1))
t = target("LN")
write("stopping", text.replace(t, t.replace('"a0*x^a1" a0=1.2e-24 a1=-1', '"a0*(x-0.019)" a0=1e-20')))
write("stopping_deep", text.replace(t, t.replace('"a0*x^a1" a0=1.2e-24 a1=-1', '"a0*(x-0.0125)" a0=1e-19')))
write("stops", text.replace(t, t.replace("1e+20 atoms/cm2", "1.7e+20 atoms/cm2")))
PY
refused "a layer without areal_density" E_FIELD_REQUIRED "areal_density" no_density.azr
refused "a yield per uC without charge_state" E_FIELD_REQUIRED "charge_state" no_charge.azr
refused "fold_windows (deleted from format 3)" E_FIELD_DELETED "fold_windows" windows.azr
refused "a stopping not positive at the beam energy" E_LAYER_STOPPING "not positive and finite" stopping.azr
refused "a stopping not positive below the beam energy" E_LAYER_STOPPING "not positive and finite" stopping_deep.azr
refused "a layer the beam does not leave" E_LAYER_EXIT_BELOW_MODEL "use a target of kind infinite" stops.azr

# ---- 6. pyazr
if [ -n "${PYAZR_ENGINE_DIR:-}" ]; then
  AZURE2_BIN="$AZURE2_BIN" python3 "$HERE/pyazr_layer_test.py" > "$WORK/pyazr.txt" 2>&1
  status=$?
  if [ "$status" = 0 ]; then ok "pyazr: $(grep -c '^  ok' "$WORK/pyazr.txt") checks of pyazr_layer_test.py"
  elif [ "$status" = 77 ]; then echo "  skip  pyazr: $(head -1 "$WORK/pyazr.txt")"
  else bad "pyazr_layer_test.py:"; grep -v '^Coulomb\|EC Amplitude' "$WORK/pyazr.txt" | sed 's/^/        /'; fi
else
  echo "  skip  pyazr (PYAZR_ENGINE_DIR not set)"
fi

echo
if [ "$failures" -ne 0 ]; then
  echo "layer_yield: $failures failure(s)"
  exit 1
fi
echo "layer_yield: all checks passed"
