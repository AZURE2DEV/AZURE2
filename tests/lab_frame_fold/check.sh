#!/usr/bin/env bash
#
# A point given at a lab angle is folded in the lab (docs/upgrade/energy_grid.md section 4.3, decision D5): at every
# node E_n of its fold the c.m. angle is theta_cm(E_n; theta_lab) and the value is divided by J(E_n; theta_lab) =
# dOmega_lab/dOmega_cm (include/RowSet.h).  models/lab_fold.azr (README.md): 7Li(p,n)7Be at 30 deg lab, E_lab = 2.3 MeV,
# 0.4 MeV above threshold, where theta_cm moves by 3 deg and J by 14 % across a 300 keV rectangular resolution.
#   folded : the point, `detector det_wide` (rectangular, width 300 keV lab);
#   fine   : the same lab angle unfolded at 301 lab energies, 1 keV apart, across the window (value 1: the c.m. data
#            column is J);
#   frozen : the same energies at the folded point's own c.m. angle, as c.m. data (the fold before D5).
# Checked, with the unfolded engine as the reference (no second formula of the cross section):
#   (a) the folded point's lab value (its c.m. model over its J) equals the brute-force average (1/w) int sigma_cm(E,
#       theta_cm(E)) / J(E) dE_lab (Simpson on `fine`) within 1e-4 (the fold's quadrature: 2.2e-5 at 200 points, 3.6e-4
#       at 50, 1.4e-6 at 800 -- second order; Simpson's error on the 1 keV grid is below 1e-9);
#   (b) the fold with the angle frozen at the nominal energy, (1/w) int sigma_cm(E, theta_i) dE_lab / J_i, differs from
#       the brute force by more than 1e-2 (measured 2.2e-2), so (a) separates the two; the build before D5 folded at the
#       frozen angle and fails (a) by 2.2e-2;
#   (c) frozen.dat's c.m. angle is the folded point's own (1e-8 deg), and the window does change the angle and J (the
#       case bites);
#   (d) a fitted energy offset with a wide prior (0.3 MeV) gives the fold a reach margin that holds nodes where 30 deg
#       lab is beyond the kinematic maximum (E_c.m. below 1.678 MeV): at offset 0 the run completes and the point is
#       finite (the margin's nodes are not read; the milestone-4 build aborted at set-up), at -0.32 MeV the window
#       reads them and the point is NaN with the failed-point warning that names the cause.
#
#   ./tests/lab_frame_fold/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/lab_frame_fold.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

echo "lab_frame_fold: $AZURE2_BIN"
cp -r "$HERE/models/." "$WORK/"
mkdir -p "$WORK/output" "$WORK/checks"
(cd "$WORK" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline lab_fold.azr > run.log 2>&1)
OUT="$WORK/output/AZUREOut_aa=1_R=2.out"
[ -s "$OUT" ] || { echo "  FAIL  no output"; tail -5 "$WORK/run.log"; echo "lab_frame_fold: FAIL"; exit 1; }

. "$HERE/../lib/edit_guard.sh"
REACH_FAILS=0
reach() {  # NAME OFFSET: the folded data set with a fitted offset (prior 0.3 MeV) at OFFSET MeV; prints its row 1
  edited "$WORK/lab_fold.azr" "$WORK/$1.azr" sed \
    "s/^    detector       det_wide\$/    detector       det_wide\n    shift          offset $2 MeV fit prior 0 0.3/" ||
    { echo "  FAIL  (d) the variant $1 was not made"; REACH_FAILS=$((REACH_FAILS + 1)); return; }
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"; cp -r "$WORK/data" "$WORK/$1/"; cp "$WORK/$1.azr" "$WORK/$1/t.azr"
  (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
}
reach reach0 0
reach reach32 -0.32
row0="$(head -1 "$WORK/reach0/output/AZUREOut_aa=1_R=2.out" 2>/dev/null | awk '{print $4}')"
row32="$(head -1 "$WORK/reach32/output/AZUREOut_aa=1_R=2.out" 2>/dev/null | awk '{print $4}')"
if [ -n "$row0" ] && [ "$row0" != nan ] && ! grep -q 'calculation failed' "$WORK/reach0/run.log"; then
  echo "  ok    (d) a margin beyond the kinematic maximum, not read: the point is $row0"
else echo "  FAIL  (d) offset 0: '$row0' (run.log: $(grep -m1 -i 'terminate\|what()\|failed' "$WORK/reach0/run.log"))"
  REACH_FAILS=$((REACH_FAILS + 1)); fi
if [ "$row32" = nan ] && grep -q 'beyond the kinematic maximum at an energy its fold reads' "$WORK/reach32/run.log"; then
  echo "  ok    (d) the window reads it at -0.32 MeV: NaN with the warning"
else echo "  FAIL  (d) offset -0.32 MeV: '$row32', no failed-point warning naming the kinematic maximum"
  REACH_FAILS=$((REACH_FAILS + 1)); fi

python3 - "$OUT" <<'PY'
import sys

rows = [list(map(float, l.split())) for l in open(sys.argv[1]) if len(l.split()) >= 7]
fails = 0


def check(name, ok, detail):
    global fails
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}: {detail}")
    fails += 0 if ok else 1


# columns: E_cm, E_x, theta_cm, model (c.m.), ..., data (c.m.: J times the lab value 1), ...
if len(rows) != 603:
    print(f"  FAIL  {len(rows)} rows, expected 1 + 301 + 301")
    sys.exit(1)
fold, fine, frozen = rows[0], rows[1:302], rows[302:603]
w, h = 0.300, 0.001


def simpson(y):
    return h / 3 * (y[0] + y[-1] + 4 * sum(y[1:-1:2]) + 2 * sum(y[2:-1:2]))


brute = simpson([r[3] / r[5] for r in fine]) / w
j_i = fold[5]
model = fold[3] / j_i
frozen_avg = simpson([r[3] for r in frozen]) / w / j_i
rel = model / brute - 1
rel_frozen = frozen_avg / brute - 1
check("(a) folded = brute force", abs(rel) < 1e-4, f"{model:.10e} against {brute:.10e}, {rel:.2e} (limit 1e-4)")
check("(b) frozen angle is worse", abs(rel_frozen) > 1e-2, f"{frozen_avg:.10e}, {rel_frozen:.2e} (must exceed 1e-2)")
check("(c) frozen at the point's angle", all(abs(r[2] - fold[2]) < 1e-8 for r in frozen),
      f"theta_i {fold[2]:.9f} deg, frozen rows {frozen[0][2]:.9f} deg")
dj = (max(r[5] for r in fine) - min(r[5] for r in fine)) / j_i
dth = max(r[2] for r in fine) - min(r[2] for r in fine)
check("(c) the window moves the angle and J", dj > 0.1 and dth > 2.0, f"J by {dj:.3f} of J_i, theta_cm by {dth:.2f} deg")
sys.exit(1 if fails else 0)
PY
status=$?
[ $REACH_FAILS -eq 0 ] || status=1
[ $status -eq 0 ] && echo "lab_frame_fold: PASS" || echo "lab_frame_fold: FAIL"
exit $status
