#!/usr/bin/env bash
#
# Relativistic two-body kinematics ("relativistic=1" in <config>, include/Kinematics.h):
#
#   (a) n + 4He s-wave hard sphere: sigma = 4 pi sin^2(k a)/k^2 with k from the invariant mass, and the
#       c.m. energies of the lab grid;
#   (b) p + 4He s-wave hard sphere + Coulomb at 90/150 deg against mpmath Coulomb functions at the
#       relativistic eta and rho;
#   (c) 5He: t(d,n)a and 4He(n,d)t at equal invariant mass obey detailed balance with relativistic k
#       (and, in the non-relativistic models, with the non-relativistic mapping and k);
#   (d) relativistic.out: the levels in the E(s) convention of the LANL code EDA;
#   (e) lab-frame differential data: the relativistic c.m. angle and dOmega_lab/dOmega_cm against an
#       explicit four-vector boost (t(d,n)a, p + 4He elastic, and 7Li(p,n)7Be just above threshold, forward branch);
#   (h) a lab angle between the relativistic and the non-relativistic maximum is refused with relativistic=1;
#   (g) relativistic=1 with the hybrid nuclear-potential Coulomb functions (<potential> block) stops with an ERROR.
#   (f) with a second binary (a build without relativistic kinematics), every non-relativistic model's
#       output is byte-identical between the two.
#
#   ./tests/relativistic/check.sh path/to/AZURE2 [path/to/reference/AZURE2]

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2 [reference AZURE2]}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
REF_BIN="${2:-}"
[ -n "$REF_BIN" ] && REF_BIN="$(cd "$(dirname "$REF_BIN")" && pwd)/$(basename "$REF_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }
PYTHON="$(command -v python3 || command -v python || true)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/relativistic.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run BIN NAME MODEL MODE -> $WORK/NAME
run() {
  local bin="$1" dir="$WORK/$2"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$HERE/models/$3.azr" "$dir/t.azr"
  cp -r "$HERE/models/data" "$dir/"
  (cd "$dir" && printf "$4\n\n\n7\n" | run_limited 300 "$bin" --no-coulomb-table --no-gui --no-readline t.azr > run.log 2>&1)
}
for m in na_hs pa_hs; do run "$AZURE2_BIN" "$m" "$m" 3; run "$AZURE2_BIN" "${m}_rel" "${m}_rel" 3; done
for m in he5 he5_ang pa_ang li7pn_ang; do run "$AZURE2_BIN" "$m" "$m" 1; run "$AZURE2_BIN" "${m}_rel" "${m}_rel" 1; done
for d in na_hs na_hs_rel pa_hs pa_hs_rel he5 he5_rel he5_ang he5_ang_rel pa_ang pa_ang_rel li7pn_ang li7pn_ang_rel; do
  if ! ls "$WORK/$d/output/"AZUREOut_* > /dev/null 2>&1; then bad "$d: no output"; tail -15 "$WORK/$d/run.log"; fi
done

if [ -n "$PYTHON" ] && "$PYTHON" -c "import numpy" > /dev/null 2>&1; then
  "$PYTHON" "$HERE/verify.py" "$WORK" || failures=$((failures + $?))
else
  echo "  skip  (a)-(d): no python3 with numpy"
fi

# (h) a 7Li(p,n) lab angle between the relativistic (44.60 deg) and non-relativistic maxima at 1.90 MeV (these masses):
#     refused with relativistic=1 (until 2026-10-29 converted with the non-relativistic angle and Jacobian, 9.5 % off,
#     silently), converted without it
run "$AZURE2_BIN" li7pn_edge li7pn_edge 1; run "$AZURE2_BIN" li7pn_edge_rel li7pn_edge_rel 1
if grep -q "beyond the kinematic maximum 44\.[0-9]*.*relativistic kinematics" "$WORK/li7pn_edge_rel/run.log" &&
   ! ls "$WORK/li7pn_edge_rel/output/"AZUREOut_* > /dev/null 2>&1 && ls "$WORK/li7pn_edge/output/"AZUREOut_* > /dev/null 2>&1; then
  ok "(h) 7Li(p,n) at 1.90 MeV, 44.9 deg: refused at the relativistic maximum (44.60 deg with these masses) with relativistic=1, converted without (Marion-Young 44.93 deg)"
else
  bad "(h) 7Li(p,n) 44.9 deg edge: $(grep -m1 ERROR "$WORK/li7pn_edge_rel/run.log"); rel outputs: $(ls "$WORK/li7pn_edge_rel/output" | grep -c AZUREOut)"
fi

# (e) the hybrid model is refused
mkdir -p "$WORK/hyb/output" "$WORK/hyb/checks"
# the hybrid model is switched on by a <potential> block after </config>
awk '{ print } /^<\/config>$/ { print "<potential>\nuseHybridPotential=1\npotentialType=0\nV0=40\nR=2.5\na=0.6\n</potential>" }' \
  "$HERE/models/na_hs_rel.azr" > "$WORK/hyb/t.azr"
(cd "$WORK/hyb" && printf "3\n\n\n7\n" | run_limited 120 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
if grep -q "ERROR: relativistic=1 cannot be combined" "$WORK/hyb/run.log" && ! ls "$WORK/hyb/output/"AZUREOut_* > /dev/null 2>&1; then
  ok "(g) relativistic=1 with the hybrid Coulomb functions stops with an ERROR"
else
  bad "(g) relativistic=1 with a <potential> block did not stop"; tail -5 "$WORK/hyb/run.log"
fi

# (f) non-relativistic models byte-identical to a reference build
if [ -n "$REF_BIN" ]; then
  for m in na_hs pa_hs he5 he5_ang pa_ang; do
    mode=1; case "$m" in na_hs|pa_hs) mode=3;; esac
    run "$REF_BIN" "ref_$m" "$m" "$mode"
    same=1
    for f in "$WORK/$m/output/"AZUREOut_*; do  # no cmp on the Windows runners
      g="$WORK/ref_$m/output/$(basename "$f")"; [ -f "$g" ] && [ "$(tr -d '\r' < "$f")" = "$(tr -d '\r' < "$g")" ] || same=0; done
    [ "$same" = 1 ] && ok "(f) $m: non-relativistic output byte-identical to the reference build" \
                    || bad "(f) $m: non-relativistic output differs from the reference build"
  done
fi

echo
[ "$failures" -eq 0 ] && echo "relativistic: PASS" || echo "relativistic: $failures FAILURE(S)"
exit "$failures"
