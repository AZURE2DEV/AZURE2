#!/usr/bin/env bash
#
# Extrapolation grids of format 3 (docs/upgrade/project_format3.md section 3.8, brief K), in the engine's
# extrapolation run (menu 3) of narrow.azr.in (a template, @GRID@ the grid; not .azr, so that the tree's runs skip
# it): 12C(p,g)13N through one 1/2+ level at E_cm = E_r = 0.3001234 MeV with Gamma_p = 1 eV and Gamma_gamma = 0.1 eV
# (E1 to the ground state).  The level width in the R-matrix is Gamma_p (a photon channel does not enter the level
# matrix), so the closed form is the single-level Breit-Wigner
# sigma(E) = (pi/k^2) omega Gamma_p Gamma_g / ((E - E_r)^2 + Gamma_p^2/4), omega = 2/(2 x 1) = 1.
#   1. `adaptive 0.2 0.4 MeV cm tolerance 1e-3 max_points 2000 start_step 0.01`: the peak value to 1e-3 of
#      (pi/k^2) 4 Gamma_g/Gamma_p; the trapezoid integral of the written cross section over the grid to 1e-3 of the
#      same model's integral on a dense grid (a `file` grid of 41800 c.m. energies: step 1e-8 MeV within 1e-4 MeV of
#      E_r, 1e-6 within 1e-2, 1e-4 elsewhere), with fewer than 400 points; `range 0.2 0.4 MeV cm points <the same
#      count>` misses the peak by more than 50 %.  The integral of the closed form over [0.2, 0.4],
#      (pi/k^2) Gamma_g 2 [atan((0.4 - E_r)/h) + atan((E_r - 0.2)/h)], h = Gamma_p/2, is printed beside it: the
#      model differs from it by 6e-3, since Gamma_p and k vary across the line and the bound 1/2- level adds a
#      non-resonant term (the written sigma over the fixed-width line is 0.96 at 50 eV below E_r, 134 at 0.2 MeV);
#   2. the adaptive grid at 1 and 3 threads: the output files identical; with AZURE_PLAN_CHECK=1 (the calculation plan
#      checked for every sweep's data set and for the refined one);
#   3. max_points 60: the warning with the largest remaining deviation, 60 points written;
#   4. `list` in c.m. and `file` in keV c.m. (relative to the project) give the energies asked for (E_cm column), and
#      --describe prints the grid file's resolved path and its number of energies;
#   5. an angle grid (dsigma/dOmega at 45 and 135 deg c.m.): one adaptive energy grid for both angles (the union),
#      the same energies at each; the 1/2+ level decays by E1 to a 1/2- state, so the distribution is isotropic, every
#      angle's curve is sigma/(4 pi) and the grid is the angle-integrated one, energy for energy;
#   6. the data check file (`file` on the Data Check line of <config>) of an adaptive grid lists the refined grid the
#      output has, point for point, not the start grid.
#
#   ./tests/extrapolation_grids/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"
run_limited() { local limit="$1"; shift; $(guard_command "$limit") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/extrapolation_grids.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# setup NAME GRID [ANGLES]: $WORK/NAME/m.azr, narrow.azr.in with GRID as the extrapolation's energies; with ANGLES, a
# differential cross section (b/sr) on that angle grid
setup() {
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  python3 - "$HERE/narrow.azr.in" "$WORK/$1/m.azr" "$2" "${3:-}" <<'PY'
import sys
src, dst, grid, angles = sys.argv[1:5]
text = open(src).read().replace("@GRID@", grid)
if angles:
    text = text.replace("    energies       " + grid, "    energies       " + grid + "\n    angles         " + angles)
    text = text.replace("quantity       sigma\n", "quantity       dsigma_domega\n").replace("unit           b\n", "unit           b/sr\n")
open(dst, "w").write(text)
PY
}
# run NAME THREADS: the extrapolation run (menu 3) of $WORK/NAME/m.azr
run() {
  (cd "$WORK/$1" && rm -rf output/* checks/* && printf '3\nn\n\n\n7\n' |
     OMP_NUM_THREADS="$2" run_limited 600 "$AZURE2_BIN" --no-gui --no-readline m.azr > run.log 2>&1)
}
OUT=output/AZUREOut_aa=1_R=2.extrap

# ---- 1. the narrow line on an adaptive grid
nint=0
setup adaptive "adaptive 0.2 0.4 MeV cm tolerance 1e-3 max_points 2000 start_step 0.01"
AZURE_PLAN_CHECK=1 run adaptive 1
if [ ! -s "$WORK/adaptive/$OUT" ]; then
  bad "adaptive: no output ($(grep -m1 -E '^ERROR|error' "$WORK/adaptive/run.log"))"
else
  n=$(awk 'NF && $1 !~ /^#/' "$WORK/adaptive/$OUT" | wc -l)
  nint=$n
  setup uniform "range 0.2 0.4 MeV cm points $n"
  run uniform 1
  setup dense 'file "dense.txt" MeV cm'
  python3 -c "
er = 0.3001234
g = {round(0.2 + k * 1e-4, 12) for k in range(2001)}
g |= {round(er + k * 1e-6, 12) for k in range(-10000, 10001)}
g |= {round(er + k * 1e-8, 14) for k in range(-10000, 10001)}
open('$WORK/dense/dense.txt', 'w').write(''.join('%.14g\\n' % x for x in sorted(x for x in g if 0.2 <= x <= 0.4)))"
  run dense 3
  if python3 - "$WORK/adaptive/$OUT" "$WORK/uniform/$OUT" "$WORK/dense/$OUT" <<'PY'
import math, sys
def read(path):
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    return [float(r[0]) for r in rows], [float(r[3]) for r in rows]
e, s = read(sys.argv[1])
eu, su = read(sys.argv[2])
ed, sd = read(sys.argv[3])
# closed form: m(p) = 1.00728 u, m(12C) = 12 u (the <levels> masses), u = 931.49410242 MeV, hbar c = 197.3269804 MeV fm
mu = 1.00728 * 12.0 / 13.00728 * 931.49410242
er = 0.3001234
k2 = 2 * mu * er / 197.3269804**2           # fm^-2 at E_r
gp, gg = 1.0e-6, 0.1e-6                     # MeV
h = gp / 2
peak = math.pi / k2 * 4 * gg / gp / 100.0   # b (1 b = 100 fm^2)
integral = math.pi / k2 * gg * 2 * (math.atan((0.4 - er) / h) + math.atan((er - 0.2) / h)) / 100.0
trapezoid = lambda x, y: sum(0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1]) for i in range(1, len(x)))
trap, dense = trapezoid(e, s), trapezoid(ed, sd)
fails = 0
def check(name, good, detail):
    global fails
    print(("  ok    " if good else "  FAIL  ") + name + ": " + detail)
    fails += not good
check("adaptive: the peak to 1e-3 of the closed form", abs(max(s) / peak - 1) < 1e-3,
      f"{max(s):.6e} b against {peak:.6e} b ({max(s) / peak - 1:+.1e})")
check("adaptive: the trapezoid integral to 1e-3 of the dense grid's", len(ed) == 41800 and abs(trap / dense - 1) < 1e-3,
      f"{trap:.6e} against {dense:.6e} MeV b ({trap / dense - 1:+.1e}; the closed form's {integral:.6e})")
check("adaptive: fewer than 400 points, ascending, from 0.2 to 0.4 MeV c.m.",
      len(e) < 400 and e == sorted(e) and abs(e[0] - 0.2) < 1e-9 and abs(e[-1] - 0.4) < 1e-9, f"{len(e)} points")
check("the uniform grid of the same count misses the peak by more than 50 %", len(eu) == len(e) and max(su) < 0.5 * peak,
      f"{len(eu)} points, largest {max(su):.3e} b")
sys.exit(1 if fails else 0)
PY
  then :; else failures=$((failures + 1)); fi
  grep -q 'Adaptive grid of extrapolation #1' "$WORK/adaptive/run.log" &&
    ok "adaptive: the run log reports the grid" || bad "adaptive: no report of the grid in the run log"
fi

# ---- 2. 1 and 3 threads
setup threads "adaptive 0.2 0.4 MeV cm tolerance 1e-3 max_points 2000 start_step 0.01"
AZURE_PLAN_CHECK=1 run threads 3
same=1
for f in $(cd "$WORK/adaptive/output" && ls); do
  case "$f" in intEC*) continue ;; esac
  cmp -s "$WORK/adaptive/output/$f" "$WORK/threads/output/$f" || { same=0; bad "1 and 3 threads: output/$f differs"; }
done
[ "$same" = 1 ] && [ -s "$WORK/threads/$OUT" ] && ok "1 and 3 threads (plan checked): every output file identical"

# ---- 3. max_points
setup budget "adaptive 0.2 0.4 MeV cm tolerance 1e-3 max_points 60 start_step 0.01"
run budget 1
n=$(awk 'NF && $1 !~ /^#/' "$WORK/budget/$OUT" 2>/dev/null | wc -l)
if grep -q 'WARNING: extrapolation #1: the adaptive grid reached max_points = 60' "$WORK/budget/run.log" && [ "$n" = 60 ]; then
  ok "max_points 60: the warning ($(grep -o 'the largest deviation.*' "$WORK/budget/run.log")), 60 points"
else
  bad "max_points 60: $n points, warning: $(grep -m1 'WARNING: extrapolation' "$WORK/budget/run.log")"
fi

# ---- 4. list and file grids in c.m.
setup list "list 0.25 0.3 0.35 MeV cm"
run list 1
setup file "file \"grid.txt\" keV cm"
printf '# E_cm (keV)\n250\n300  # the line\n350\n' > "$WORK/file/grid.txt"
run file 1
for g in list file; do
  got=$(awk 'NF && $1 !~ /^#/ {printf "%.9f ", $1}' "$WORK/$g/$OUT" 2>/dev/null)
  [ "$got" = "0.250000000 0.300000000 0.350000000 " ] && ok "$g grid in c.m.: E_cm = 0.25, 0.3, 0.35 MeV" ||
    bad "$g grid in c.m.: E_cm = $got"
done
desc=$(cd "$WORK/file" && "$AZURE2_BIN" --describe m.azr 2>&1)
echo "$desc" | grep -q 'grid file: grid.txt (3 energies)' && ok "--describe: the grid file and its 3 energies" ||
  bad "--describe: $(echo "$desc" | grep -i grid)"

# ---- 5. an angle grid: the union of the angles' needs
setup angles "adaptive 0.2 0.4 MeV cm tolerance 1e-3 max_points 2000 start_step 0.01" "list 45 135 deg cm"
run angles 1
a45=$(awk 'NF && $3 + 0 == 45 {print $1}' "$WORK/angles/$OUT" 2>/dev/null)
a135=$(awk 'NF && $3 + 0 == 135 {print $1}' "$WORK/angles/$OUT" 2>/dev/null)
na=$(echo "$a45" | grep -c .)
aint=$(awk 'NF && $1 !~ /^#/ {print $1}' "$WORK/adaptive/$OUT" 2>/dev/null)
if [ "$nint" -gt 0 ] && [ "$a45" = "$aint" ] && [ "$a45" = "$a135" ]; then
  ok "an angle grid: the $na energies of the angle-integrated grid at 45 and at 135 deg (isotropic)"
else
  bad "an angle grid: $na energies at 45 deg, $(echo "$a135" | grep -c .) at 135 deg ($(grep -m1 ERROR "$WORK/angles/run.log"))"
fi

# ---- 6. the data check file lists the refined grid
setup datachk "adaptive 0.2 0.4 MeV cm tolerance 1e-3 max_points 2000 start_step 0.01"
sed -i 's/^none\( *#Data Check\)/file\1/' "$WORK/datachk/m.azr"
run datachk 1
# data.chk: the rows of the "Data" table (segment 1), their c.m. energies (column 4)
chk=$(awk '/^\*  *Data  *\*/ {d = 1; next} d && $1 == 1 && NF >= 7 {print $4}' "$WORK/datachk/checks/data.chk" 2>/dev/null)
nchk=$(echo "$chk" | grep -c .)
nout=$(awk 'NF && $1 !~ /^#/' "$WORK/datachk/$OUT" 2>/dev/null | wc -l)
if grep -q '^file *#Data Check' "$WORK/datachk/m.azr" && [ "$nout" -gt 21 ] && [ "$nchk" = "$nout" ] &&
   [ "$(echo "$chk" | head -1)" = 0.2 ] && [ "$(echo "$chk" | tail -1)" = 0.4 ]; then
  ok "data check file: the $nchk points of the refined grid (the start grid has 21)"
else
  bad "data check file: $nchk points listed, $nout calculated"
fi

echo
[ "$failures" -eq 0 ] && echo "extrapolation_grids: all checks passed" || echo "extrapolation_grids: $failures failure(s)"
exit $((failures > 0))
