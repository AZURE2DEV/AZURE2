#!/usr/bin/env bash
#
# The reaction rate (mode 5) equals N_A <sigma v> integrated from a dense extrapolation of
# the same model (expected/rates_*.dat, from make_truth.sh; mode 3 does not go through
# ReactionRate.cpp).  Five temperatures, T9 = 0.01 to 1, relative tolerance 1e-3:
#
#   (a) p12C       charged entrance (the Gamow tail): the spline of sigma against E lay
#                  0.9 % low at T9 = 0.01 (now ln sigma against ln E; audit F1);
#   (b) n12C       neutron entrance, sigma ~ 1/v: sigma was zero below 1 keV, 49 % low at
#                  T9 = 0.01 (now from 0.1 eV; audit F2);
#   (c) n12C_pert  the same through the per-temperature quadrature (useAdaptiveGrid=0);
#   (d) moved      a 22 eV wide 1/2+ level at E_x = 2.0935 MeV in the .azr, evaluated with the
#                  parameter file of the same model at 2.0985 MeV (expected/moved.par): the grid
#                  clusters followed the .azr value, the resonance fell between backbone nodes
#                  and the rate lost it (99.5 % low at T9 = 0.1 to 1; audit F3);
#   (e) moved_azr  the model with 2.0985 MeV in the .azr and no parameter file gives the same
#                  rates as (d) to 1e-9: the grid depends on the parameters, not on the file.
#
# The Boltzmann constant moved from 8.6171e-2 to the CODATA 8.617333262e-2 MeV/GK in the
# same change (audit F6); the truth uses CODATA throughout.
#
#   ./tests/reaction_rate/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/reaction_rate.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

# rate <model> <parameter file or ""> -> $WORK/<model>/output/reactionrates.out
rate() {
  local model="$1" par="$2" d="$WORK/$1"
  mkdir -p "$d/output" "$d/checks"
  cp "$HERE/models/$model.azr" "$d/t.azr"
  cp "$HERE/temperatures.dat" "$d/"
  [ -n "$par" ] && cp "$HERE/expected/$par" "$d/"
  (cd "$d" && printf '5\n%s\n1\n2\nyes\ntemperatures.dat\n' "$par" |
     run_limited 600 "$AZURE2_BIN" --no-gui --no-readline t.azr > run.log 2>&1)
  local status=$?
  if [ "$status" -ne 0 ] || [ ! -s "$d/output/reactionrates.out" ]; then
    bad "$model: the rate run failed (exit status $status)"; tail -5 "$d/run.log"; return 1
  fi
  return 0
}

# compare <label> <rates.out> <truth> <tol>: every temperature within tol
compare() {
  local label="$1" got="$2" want="$3" tol="$4"
  local worst
  worst="$(paste <(awk 'NR > 1 && NF >= 2 { print $1, $2 }' "$got") <(awk 'NR > 1 && NF >= 2 { print $2 }' "$want") |
    awk -v tol="$tol" 'BEGIN { w = 0; n = 0 }
         NF == 3 { d = ($2 - $3) / $3; if (d < 0) d = -d; if (d > w) { w = d; wt = $1; wg = $2; ww = $3 }; n++ }
         END { if (n != 5) { print "nrows", n; exit } printf "%.2e T9=%s got %s want %s\n", w, wt, wg, ww }')"
  case "$worst" in
    nrows*) bad "$label: $worst temperature rows (want 5)"; return ;;
  esac
  if awk -v w="${worst%% *}" -v tol="$tol" 'BEGIN { exit !(w <= tol) }'; then
    ok "$label: largest deviation $worst (tol $tol)"
  else
    bad "$label: largest deviation $worst (tol $tol)"
  fi
}

TOL=1e-3
rate p12C ""      && compare "(a) p12C, charged entrance"              "$WORK/p12C/output/reactionrates.out"      "$HERE/expected/rates_p12C.dat" $TOL
rate n12C ""      && compare "(b) n12C, 1/v neutron entrance"          "$WORK/n12C/output/reactionrates.out"      "$HERE/expected/rates_n12C.dat" $TOL
rate n12C_pert "" && compare "(c) n12C, per-temperature quadrature"    "$WORK/n12C_pert/output/reactionrates.out" "$HERE/expected/rates_n12C.dat" $TOL
rate moved moved.par && compare "(d) moved, narrow level from the parameter file" "$WORK/moved/output/reactionrates.out" "$HERE/expected/rates_moved_azr.dat" $TOL
rate moved_azr ""    && compare "(e) moved_azr, the same level in the .azr"      "$WORK/moved_azr/output/reactionrates.out" "$WORK/moved/output/reactionrates.out" 1e-9

if [ "$failures" -eq 0 ]; then echo "reaction_rate: PASS"; else echo "reaction_rate: FAIL ($failures)"; exit 1; fi
