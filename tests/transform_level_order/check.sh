#!/usr/bin/env bash
#
# The standard (non-Brune) parameter transformations hand the eigenvalues to the levels by energy, not by
# their position in the file (audit item 13, 2026-10-27).
#
# CNuc::TransformIn solves a generalized eigenproblem whose eigenvalues come back in ascending order and
# assigned them to the levels of a J^pi group in FILE order; TransformOut followed the same positional
# index.  A group listed out of energy order (a background pole first, say) therefore had its energies,
# widths and fix flags permuted: the first-listed level carried the lowest formal energy, and a fit that
# fixed the pole varied the low level instead.  The cross sections were never affected (the set of levels
# is the same whatever the labels), which is why nothing warned.  Only runs with --no-brune (or the Brune
# box unchecked in the GUI) go through the transformation.
#
# Models (make_models.py): a 1+ group of n + 15N with levels 3.491 (energy free) and 3.791 (energy fixed),
# listed ascending (asc) and with the 3.791 level first (desc); truth = asc with the free level at 3.471.
#
#   (a) mode 3 with --no-brune: asc and desc give the same cross sections (labels do not matter);
#   (b) desc: param.par's energy_1 belongs to the first-listed level (3.791, the boundary-condition
#       reference, so its formal energy is 3.791 exactly) and energy_2 is the formal energy of 3.491;
#   (c) desc: parameters.out lists the levels in file order (3.791 then 3.491);
#   (d) a --no-brune fit of desc_fit (data from truth, one free parameter: the 3.491 level's energy)
#       moves the free level to 3.471 and leaves the fixed one at 3.791, chi2 ~ 0.
#   Old build: (b) energy_1 = 3.4888, (c) 3.491 first, (d) the fixed flag sits on the 3.491 level and
#   the fit moves the 3.791 level instead.
#
#   ./tests/transform_level_order/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/transform_level_order.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }
near() { awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(a != "" && d <= t) }'; }

prep() { mkdir -p "$WORK/$1/output" "$WORK/$1/checks"; cp "$HERE/models/$2" "$WORK/$1/t.azr"; }
extrap() { (cd "$WORK/$1" && printf '3\nn\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-brune --no-gui --no-readline t.azr > run.log 2>&1); }
energies() { awk '/E_level/ { for (i = 1; i <= NF; i++) if ($i == "E_level") print $(i + 2) }' "$WORK/$1/output/parameters.out" | tr '\n' ' '; }
param() { awk -v n="$2" '$1 == n { print $2 }' "$WORK/$1/output/param.par"; }

# ---- (a) the cross sections do not depend on the listing order
prep asc asc.azr; extrap asc
prep desc desc.azr; extrap desc
f="$WORK/asc/output/AZUREOut_aa=1_TOTAL_CAPTURE.extrap"; g="$WORK/desc/output/AZUREOut_aa=1_TOTAL_CAPTURE.extrap"
[ -s "$f" ] && [ -s "$g" ] || { bad "no extrapolation output"; tail -5 "$WORK/desc/run.log"; exit 1; }
w="$(paste "$f" "$g" | awk 'NF >= 10 { d = $4 - $9; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300; r = d / m; if (r > w) w = r; n++ } END { printf "%.2e %d", w + 0, n + 0 }')"
if [ "${w##* }" -gt 5 ] && awk -v r="${w%% *}" 'BEGIN { exit !(r + 0 <= 1e-10) }'; then ok "(a) asc and desc give the same total capture cross section (${w##* } energies, worst ${w%% *})"
else bad "(a) asc and desc differ: $w"; fi

# ---- (b) the formal energies belong to the listed levels
e1="$(param desc energy_1)"; e2="$(param desc energy_2)"
if near "$e1" 3.791 1e-6 && awk -v e="$e2" 'BEGIN { exit !(e + 0 > 3.4 && e + 0 < 3.6) }'; then ok "(b) desc: energy_1 = $e1 is the first-listed (fixed) 3.791 level, energy_2 = $e2 the 3.491 level"
else bad "(b) desc: energy_1 = $e1, energy_2 = $e2 (expected 3.791 and ~3.489: the first-listed level got the lowest eigenvalue)"; fi
a1="$(param asc energy_1)"; a2="$(param asc energy_2)"
if near "$a1" 3.491 1e-6 && near "$a2" 3.7921077 1e-5; then ok "(b) asc unchanged: energy_1 = $a1, energy_2 = $a2"; else bad "(b) asc: energy_1 = $a1, energy_2 = $a2"; fi

# ---- (c) parameters.out in file order
es="$(energies desc)"
if [ "$es" = "3.7910 3.4910 " ]; then ok "(c) desc: parameters.out lists the levels in file order ($es)"; else bad "(c) desc: parameters.out order '$es', expected '3.7910 3.4910 '"; fi

# ---- (d) a fit varies the level whose energy is free
# Truth: desc itself with the 3.491 level's formal energy (3.4888353 under the boundary condition of the
# first-listed level, whichever name a build gives it) lowered by 0.02 MeV through a parameter file, so the
# truth lies in the fit's own one-parameter space on every build.  The fit has one free parameter, the
# energy of the 3.491 level, and must reach it.
prep truth desc.azr
awk '$1 ~ /^energy_/ && $2 + 0 > 3.48 && $2 + 0 < 3.50 { $2 = "3.4688353e+00"; print; next } { print }' "$WORK/desc/output/param.par" > "$WORK/truth/truth.sav"
grep -q "3.4688353e+00" "$WORK/truth/truth.sav" || { bad "(d) no formal energy near 3.4888 in desc's param.par"; cat "$WORK/desc/output/param.par"; exit 1; }
(cd "$WORK/truth" && printf '3\nn\ntruth.sav\n\n7\n' | run_limited 300 "$AZURE2_BIN" --no-brune --no-gui --no-readline t.azr > run.log 2>&1)
te="$(energies truth | tr ' ' '\n' | sort -g | head -1)"       # the moved level's physical energy
prep fit desc_fit.azr; mkdir -p "$WORK/fit/data"
awk 'NF >= 4 && $4 > 0 { printf "%.6f  0  %.8e  %.8e\n", $1 * 16.00876 / 15.0001, $4, 0.01 * $4 }' "$WORK/truth/output/AZUREOut_aa=1_TOTAL_CAPTURE.extrap" > "$WORK/fit/data/cap.dat"
n="$(grep -c . "$WORK/fit/data/cap.dat")"
[ "$n" -ge 10 ] && near "$te" 3.471 5e-3 || { bad "(d) truth: $n data points, moved level at '$te'"; tail -3 "$WORK/truth/run.log"; exit 1; }
(cd "$WORK/fit" && printf '2\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-brune --no-gui --no-readline t.azr > run.log 2>&1)
chi2="$(sed -n 's/^Total-Chi-Squared: \([^ ]*\).*/\1/p' "$WORK/fit/output/chiSquared.out" | tail -1)"
es="$(energies fit)"; fe="${es%% *}"; me="$(echo "$es" | awk '{ print $2 }')"
if near "$fe" 3.791 1e-4 && near "$me" "$te" 1e-3 && awk -v c="$chi2" 'BEGIN { exit !(c != "" && c + 0 < 1) }'; then ok "(d) the fit moved the free level to $me (truth $te) and kept the fixed one at $fe (chi2 $chi2, $n points)"
else bad "(d) fit: levels '$es', chi2 '$chi2' (expected 3.791 fixed, $te fitted, chi2 ~ 0)"; fi

echo
if [ "$failures" -eq 0 ]; then echo "transform_level_order: PASS"; else echo "transform_level_order: FAIL ($failures)"; exit 1; fi
