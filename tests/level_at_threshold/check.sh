#!/usr/bin/env bash
#
# A level energy exactly at, or within round-off of, a particle threshold.
#
# The shift function, penetrability and dS/dE at a level energy (CNuc's
# boundary conditions, shift functions and parameter transformations, the THM
# vertex=constant shift, the external-capture integral) all go through
# ChannelFunc now, which continues S through threshold and treats a channel
# with Sommerfeld parameter above 100 as closed.  Before, they called the
# Coulomb functions at E_res >= 0 (NaN at 0, NaN/garbage or seconds per call
# just above for large eta) and ShftFunc's dS/dE, which mirrored its probe above
# threshold back below it (dS/dE = 0 at a level exactly at threshold).
#
# Three one-level projects, each with the level at E_res = -1e-9, 0 and +1e-9
# MeV in the entrance channel, whose reduced width is given as an amplitude
# (gammaIsRWA) so that it means the same thing on both sides:
#
#   c12c12  12C+12C -> a+20Ne and 12C+12C elastic, 2+ (eta ~ 1e4 at 1e-9 MeV)
#   be8     p+7Li -> a+a, 2+
#   thm     tests/18O_p_a_thm (THM, vertex=constant takes B_c from this level)
#
# Every run must finish within 60 s with a finite chi2 and finite output, and
# the model cross sections of the three runs of a project (every AZUREOut
# file, 10 digits) must agree to 1e-6 relative: the model is continuous
# through threshold.
#
#   ./tests/level_at_threshold/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/level_at_threshold.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-60}"; else RUN=""; fi
IDENT="$HERE/../identical_entrance_reaction/projects"
fail=0

# make NAME SRC_AZR DATA_DIR SEP EXPR -- copy SRC_AZR, then on every level line
# of the first level (field 9 == 1): field 3 = SEP + EXPR (printed exactly),
# and on the line whose field 6 is the entrance pair key 1: gamma (field 12) =
# the amplitude RWA, fields 32 and 33 = 0 1 (binding untouched on a THM line).
make() {
  local d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$3" "$d/data"
  awk -v sep="$4" -v eres="$5" -v rwa="$6" '
    /<levels>/ { L = 1; print; next } /<\/levels>/ { L = 0 }
    L && NF >= 31 && $9 == 1 {
      $3 = sprintf("%.12f", sep + eres)
      if ($6 == 1) { $12 = rwa; if (NF < 32) $32 = 0; $33 = 1 }
    }
    { print }' "$2" > "$d/run.azr"
}

# run NAME -> total chi2, or empty on failure (reason printed)
run() {
  local d="$WORK/$1" t0 t1 status
  t0=$(date +%s)
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 \
     | head -c 1000000 > log)
  status=${PIPESTATUS[0]}
  t1=$(date +%s)
  if [ ! -f "$d/output/chiSquared.out" ]; then
    echo "  FAIL  $1: no chiSquared.out ($((t1 - t0)) s)" >&2; tail -3 "$d/log" | sed 's/^/        /' >&2; return
  fi
  if [ $((t1 - t0)) -gt 60 ]; then echo "  FAIL  $1: took $((t1 - t0)) s" >&2; return; fi
  if grep -qiE 'nan|inf' "$d/output/chiSquared.out" "$d"/output/AZUREOut_*.out; then
    echo "  FAIL  $1: non-finite output" >&2; grep -iE 'nan|inf' "$d"/output/*.out | head -3 >&2; return
  fi
  grep -oE 'Total-Chi-Squared: [0-9.eE+-]+' "$d/output/chiSquared.out" | awk '{ print $2 }'
}

case_set() {  # case_set LABEL SRC_AZR DATA_DIR SEP RWA
  local label="$1" c
  local -a chi=()
  for eres in -1e-9 0 1e-9; do
    make "$label$eres" "$2" "$3" "$4" "$eres" "$5"
    c="$(run "$label$eres")"
    [ -n "$c" ] || { fail=1; chi+=(""); continue; }
    chi+=("$c")
    echo "  $label E_res = $eres MeV: chi2 $c"
  done
  [ -n "${chi[0]}" ] && [ -n "${chi[1]}" ] && [ -n "${chi[2]}" ] || { fail=1; return; }
  # The model (column 4 of every AZUREOut file, 10 digits) at E_res = +-1e-9
  # against E_res = 0.
  local worst
  worst="$(for f in "$WORK/${label}0"/output/AZUREOut_*.out; do
      b="$(basename "$f")"
      paste "$WORK/${label}-1e-9/output/$b" "$f" "$WORK/${label}1e-9/output/$b"
    done | awk '{ n = NF / 3; if (n < 4) next
      m = $(n + 4); if (m == 0) next
      for (k = 0; k <= 2; k += 2) { d = ($(k * n + 4) - m) / m; if (d < 0) d = -d; if (d > w) w = d }
    } END { printf "%.3g\n", w + 0 }')"
  if awk -v w="$worst" 'BEGIN { exit !(w < 1e-6) }'; then
    echo "  ok    $label: model continuous through threshold (max rel diff $worst)"
  else
    echo "  FAIL  $label: model jumps across threshold (max rel diff $worst)"; fail=1
  fi
}

# Data: a few points away from threshold (lab energies), 10 % errors.
mkdir -p "$WORK/d_c12" "$WORK/d_be8"
printf '%s 0 1e-10 1e-11\n' 3.0 4.0 5.0 6.0 > "$WORK/d_c12/c12_a_int.dat"
printf '%s 0 1e-12 1e-13\n' 3.0 4.0 5.0 6.0 > "$WORK/d_c12/c12_el_int.dat"
printf '%s 0 1e-3 1e-4\n' 0.3 0.6 1.0 1.5 > "$WORK/d_be8/p7li_aa_int.dat"
printf '%s 0 1e-3 1e-4\n' 17.0 18.0 19.0 20.0 > "$WORK/d_be8/aa_p7li_int.dat"
grep -v 'c12_a_diff' "$IDENT/c12c12.azr" > "$WORK/c12c12.azr"

case_set c12c12 "$WORK/c12c12.azr" "$WORK/d_c12" 13.9336 0.3
case_set be8 "$IDENT/be8.azr" "$WORK/d_be8" 17.2551 0.5
case_set thm "$HERE/../18O_p_a_thm/18O_p_a_thm.azr" "$HERE/../18O_p_a_thm/data" 7.993600 -0.407

if [ "$fail" -ne 0 ]; then echo "FAIL"; exit 1; fi
echo "PASS: levels at E_res = -1e-9, 0, +1e-9 MeV: finite, fast, continuous (12C+12C, p+7Li, THM p+18O)"
