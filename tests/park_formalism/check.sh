#!/usr/bin/env bash
#
# --use-park must reproduce the default (Brune) calculation.
#
# Park's level matrix (Phys. Rev. C 104, 064612, Eqs. 9, 21, 22, 26, 28) is
# Brune's alternative level matrix (Phys. Rev. C 66, 044611, Eq. 33) with each
# basis state rescaled: gamma_Park = gamma_Brune * sqrt(J), J = 1 - sum_c
# gamma_Park^2 dS_c/dE.  The collision matrix is invariant under that
# rescaling, so a project read from the same .azr (observed energies and
# widths) has to give the same cross sections in both modes -- while the two
# modes build their level matrices from different equations and different
# reduced width amplitudes.
#
# Models: tests/13N (12C+p elastic scattering and capture, external capture
# included) and tests/identical_pp_res (a two-channel 3P2-3F2 group).
#
#   ./tests/park_formalism/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TOL="${TOL:-1e-5}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/park.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# run PROJECT MODE [FLAG]: mode-1 calculate in a private copy, print total chi2.
run() {
  local dir="$WORK/$1.$2"
  mkdir -p "$dir/output" "$dir/checks"
  cp -R "$HERE/../$1/data" "$dir/data"
  cp "$HERE/../$1/$1.azr" "$dir/"
  (cd "$dir" && printf '1\n\n\n7\n' | "$AZURE2_BIN" --no-gui --no-readline ${3:-} "$1.azr" > run.log 2>&1)
  awk '/^Total-Chi-Squared:/ { print $2 }' "$dir/output/chiSquared.out" 2>/dev/null
}

for project in 13N identical_pp_res; do
  b="$(run "$project" brune)"
  p="$(run "$project" park --use-park)"
  if [ -z "$b" ] || [ -z "$p" ]; then
    bad "$project: no chi-squared (Brune '$b', Park '$p')"
  elif awk -v a="$b" -v c="$p" -v t="$TOL" 'BEGIN { d = a - c; if (d < 0) d = -d; exit !(d <= t * a) }'; then
    ok "$project: Brune $b = Park $p"
  else
    bad "$project: Brune $b, Park $p"
  fi
  # The two modes must not be the same calculation in disguise: the fitted
  # reduced width amplitudes differ by the factor sqrt(J).
  if cmp -s "$WORK/$project.brune/output/param.par" "$WORK/$project.park/output/param.par"; then
    bad "$project: param.par identical in both modes (--use-park had no effect)"
  else
    ok "$project: reduced width amplitudes differ between the modes"
  fi
done

[ "$failures" -eq 0 ] && echo "  PASS" || echo "  FAIL ($failures)"
[ "$failures" -eq 0 ]
