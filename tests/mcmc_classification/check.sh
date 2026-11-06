#!/usr/bin/env bash
#
# MCMC parameter classification: every sampled parameter gets its kind, so the
# level-energy spread and the automatic norm/shift priors are applied.
#
# AZURECalcMCMC::BuildAutoPriors gives each parameter a kind (level energy,
# width, norm, energy shift, sqrt(E) coefficient, THM coherent background).  It
# used to re-walk the levels and segments, which gave an entry to a free THM
# norm (profiled, not a parameter) and none to a coherent background
# (cbkg_*): the classification and the parameter vector then had different
# lengths ("parameter classification produced 6 entries but there are 7
# parameters"), and the run dropped the automatic priors AND the level-energy
# spread, scattering every energy by the percentage spread instead (5 % of a
# 4.6 MeV level is 230 keV).  Kinds now come from the parameter names.
#
# Runs (24 walkers, 100 steps, level-energy spread 0.5 keV, physical
# parameters, a few seconds each):
#   pp_sqrt     tests/energy_shift_sqrt with segment 1's norm (5 %), shift
#               (2 keV) and sqrt(E) coefficient (0.001) free: 1 automatic norm
#               prior, 2 shift priors;
#   pp_plain    the same with the sqrt(E) coefficient fixed: 1 and 1;
#   thm         tests/thm_coherent_background/toy with its level energy free:
#               the free THM norm is profiled (no parameter);
#   thm_cbkg    the same with a free coherent background (Re, Im) and a prior
#               on Re defined in <mcmc>: no automatic prior for cbkg, the
#               user's kept.
# Each must classify without the warning, report the level energies as
# limited to 0.5 keV, and start its walkers within a few keV of each level
# (range of each energy over the walkers of the first step < 10 keV, > 0).
#
#   ./tests/mcmc_classification/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-600}")"
export OMP_NUM_THREADS=1

WORK="$(mktemp -d "${TMPDIR:-/tmp}/mcmc_class.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
fail=0

# The CLI's MCMC questions: parameter file (blank), walkers, steps, spread %
# (default), level-energy spread keV, threads, RWA, overwrite.
INPUT='6\n\n24\n100\n\n0.5\n1\nno\nyes\n'

# mcmc NAME: run the sampler on $WORK/NAME/run.azr; status in NAME/status.
mcmc() {
  local d="$WORK/$1"
  (cd "$d" && printf '%b' "$INPUT" | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 \
     | head -c 5000000 | tr '\r' '\n' > log
   echo "${PIPESTATUS[1]}" > status)
}
stage() {  # NAME SOURCE_DIR AZR
  mkdir -p "$WORK/$1/output" "$WORK/$1/checks"
  cp -R "$2/data" "$WORK/$1/data"
  tr -d '\r' < "$2/$3" > "$WORK/$1/source.azr"
}
# energy_ranges NAME: the range (max - min) over the walkers of the first
# recorded step of every level-energy column, one per line.  The columns are
# found from the "Free Parameters for MCMC" table of the log, in order.
energy_ranges() {
  local d="$WORK/$1"
  local cols
  # \r stripped: a Windows build writes the log with CRLF, which hides "^-+$".
  cols="$(tr -d '\r' < "$d/log" | awk '/^Free Parameters for MCMC:/ { t = 1; next }
               t == 1 && /^-+$/ { t = 2; next }
               t == 2 && NF == 0 { next }   # the Windows build puts a blank line after every row
               t == 2 && !/(Yes|No) \(/ { exit }
               t == 2 { k++; if (index($0, "Energy (MeV)")) printf "%d ", 5 + k }')"
  [ -n "$cols" ] || return
  tr -d '\r' < "$d/output/samples.mcmc" | awk -F, -v cols="$cols" '
    NR == 1 { n = split(cols, c, " "); next }
    first == "" { first = $1 }
    $1 != first { exit }
    { for (i = 1; i <= n; i++) { v = $(c[i]); if (!(i in lo) || v < lo[i]) lo[i] = v; if (!(i in hi) || v > hi[i]) hi[i] = v } }
    END { for (i = 1; i <= n; i++) printf "%.6g\n", hi[i] - lo[i] }'
}
# check NAME N_ENERGIES AUTOPRIOR_LINE [EXTRA_LINE]
check() {
  local d="$WORK/$1"
  if [ "$(cat "$d/status")" != 0 ] || ! grep -q "MCMC sampling completed successfully" "$d/log"; then
    bad "$1: the run failed (status $(cat "$d/status"))"; tail -5 "$d/log" | sed 's/^/        /'; return
  fi
  if grep -q "parameter classification produced" "$d/log"; then
    bad "$1: $(grep "parameter classification produced" "$d/log" | head -1)"
  else ok "$1: classification matches the parameter vector"; fi
  local want="Initial spread: $2 level-energy parameter"
  if grep "$want" "$d/log" | grep -qE "limited to (0\.5|5\.000000e-01) keV"; then
    ok "$1: $2 level energies limited to 0.5 keV"
  else bad "$1: no '$want... limited to 0.5 keV' line"; fi
  if grep -qF "$3" "$d/log"; then ok "$1: $3"; else bad "$1: no '$3'"; grep "Automatic priors" "$d/log" | sed 's/^/        /'; fi
  if [ -n "${4:-}" ]; then
    if grep -qF "$4" "$d/log"; then ok "$1: $4"; else bad "$1: no '$4'"; fi
  fi
  local ranges nr
  ranges="$(energy_ranges "$1")"
  nr="$(printf '%s\n' "$ranges" | awk 'NF { n++ } END { print n + 0 }')"
  if [ "$nr" != "$2" ]; then
    bad "$1: $nr level-energy columns in samples.mcmc, expected $2"
    # What the parser saw (sed l shows \r and other invisible characters).
    awk '/Free Parameters for MCMC/ { t = 14 } t > 0 { print; t-- }' "$d/log" | sed -n 'l' | sed 's/^/        log: /'
    head -2 "$d/output/samples.mcmc" 2>&1 | cut -c1-200 | sed 's/^/        samples: /'
    return
  fi
  if printf '%s\n' "$ranges" | awk '{ if (!($1 > 0 && $1 < 0.01)) bad = 1 } END { exit bad }'; then
    ok "$1: walkers start within $(printf '%s\n' "$ranges" | sort -g | tail -1 | awk '{ printf "%.2g", $1 * 1000 }') keV of each level"
  else bad "$1: level-energy ranges over the first step (MeV): $(echo $ranges)"; fi
}

ESQ="$HERE/../energy_shift_sqrt"
TOY="$HERE/../thm_coherent_background/toy"

# Segment 1 of energy_shift_sqrt: norm free (5 %), shift free (0.004 +- 0.002
# MeV), sqrt(E) coefficient free (0.003 +- 0.001) or fixed.
stage pp_sqrt "$ESQ" energy_shift_sqrt.azr
awk '!done && index($0, "data/xs_40.dat") { $10 = 1; $11 = 5; $13 = 0.002; $14 = 1; sub(/sqrtshift 0.003 0 0/, "sqrtshift 0.003 0.001 1"); done = 1 } { print }' \
  "$WORK/pp_sqrt/source.azr" > "$WORK/pp_sqrt/run.azr"
stage pp_plain "$ESQ" energy_shift_sqrt.azr
awk '!done && index($0, "data/xs_40.dat") { $10 = 1; $11 = 5; $13 = 0.002; $14 = 1; done = 1 } { print }' \
  "$WORK/pp_plain/source.azr" > "$WORK/pp_plain/run.azr"
# The toy's level energy free (field 4 of both channel lines).
free_toy() { awk '/<levels>/ { l = 1; print; next } /<\/levels>/ { l = 0 } l && NF > 12 { $4 = 0 } { print }' "$1"; }
stage thm "$TOY" toy.azr
free_toy "$WORK/thm/source.azr" > "$WORK/thm/run.azr"
stage thm_cbkg "$TOY" toy.azr
{ free_toy "$WORK/thm_cbkg/source.azr"
  printf '<thm>\nexperiment[E] segments=1 cbackground=1/2+:2=0.05,-0.08\n</thm>\n'
  printf '<mcmc>\n<parameters>\ncbkg_E_1/2+_2_1/2,0,1/2,0_re0 0.05 0.05 0.01 1\n</parameters>\n</mcmc>\n'
} > "$WORK/thm_cbkg/run.azr"

for v in pp_sqrt pp_plain thm thm_cbkg; do mcmc "$v"; done

echo "== sqrt(E) coefficient free"
check pp_sqrt 5 "Automatic priors: 1 normalization, 2 energy shifts"
echo "== sqrt(E) coefficient fixed"
check pp_plain 5 "Automatic priors: 1 normalization, 1 energy shift "
echo "== THM, profiled norm"
check thm 1 "Automatic priors: 0 normalizations, 0 energy shifts"
echo "== THM, coherent background"
check thm_cbkg 1 "2 THM coherent-background parameters: 1 with a user-defined prior, the others uniform (no automatic prior)" \
  "Found 1 Gaussian prior(s) in AZURE2 file."

[ "$fail" -eq 0 ] && echo "  PASS" || echo "  FAIL"
[ "$fail" -eq 0 ]
