#!/usr/bin/env bash
#
# A nuisance prior (a <parameterSettings> row with the nuisance flag) must be
# evaluated on the parameter it names, also when fixed parameters precede it.
#
# The chi-squared function receives the FULL parameter vector (fixed parameters
# included, as CNuc::FillCompoundFromParams reads it), while the limits manager
# numbers the parameters it knows by their position among the non-fixed ones.
# AZURECalc used to index the full vector with that non-fixed index, so the
# prior on energy_2 of tests/7Li_p_a (preceded by the fixed zero widths
# width_1_6 and width_1_7) was evaluated on width_1_6 = 0: a penalty of
# ((0 - 27.494)/0.01)^2 = 7.6e6 for a prior centred on the parameter's own
# value.
#
# Model: tests/7Li_p_a, "Level 2 Energy (MeV)" = energy_2 = 27.494 MeV (free).
#   none  no <parameterSettings>                        -> chi2_0
#   (a)   prior 27.494 +- 0.01                          -> chi2_0 exactly
#   (b)   prior 27.504 +- 0.01 (value is 1 sigma off)   -> chi2_0 + 1 (rel 1e-6)
#   (c)   (b) with an external param.par that fixes energy_1, so the fixed
#         flags the limits manager was applied with differ from a freshly
#         filled parameter list                         -> chi2 of the same
#         run without the prior + 1 (param.par rounds the values, so its own
#         baseline)
# The "Total Chi-Squared:" line on stdout includes the priors.  Since 2026-11-09 (docs/upgrade/objective.md D6,
# step 5) the nuisance priors are in every path:
#   (d)   chiSquared.out of (b) reports Total-Nuisance-Chi-Squared = 1 (rel 1e-6);
#   (e)   an MCMC on the model of (b) (24 walkers, 100 steps, physical parameters; its set-up only) applies the prior
#         ("Nuisance priors from <parameterSettings>: 1").
#
#   ./tests/nuisance_prior/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"
SRC="$HERE/../7Li_p_a"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/nuisance_prior.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

# run NAME "SETTINGS ROW or empty" [param file] -> prints the total chi2
run() {
  local d="$WORK/$1" row="$2" par="${3:-}"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$SRC/data" "$d/"
  awk -v row="$row" '
    {print}
    /<\/targets>/ && row != "" {print "<parameterSettings>"; print row; print "</parameterSettings>"}
  ' "$SRC/7Li_p_a.azr" > "$d/run.azr"
  # Relative to the run directory: read from stdin, a POSIX /tmp/... path is
  # not translated by MSYS and the native Windows binary cannot open it.
  [ -z "$par" ] || { cp "$par" "$d/start.par"; par="start.par"; }
  (cd "$d" && printf '1\n%s\n\n7\n' "$par" | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 \
     | head -c 2000000 > run.log)
  grep -oE 'Total Chi-Squared: [0-9.eE+-]+' "$d/run.log" | awk '{print $3}'
}

row_a="Level 2 Energy (MeV) 27.494 0 0 0.01 0 1 level 6"
row_b="Level 2 Energy (MeV) 27.504 0 0 0.01 0 1 level 6"

c0="$(run none "")"
ca="$(run a "$row_a")"
cb="$(run b "$row_b")"
# param.par of the plain run with energy_1 marked fixed.
awk '{ sub(/\r$/, "") } $1 == "energy_1" {print $0 "  fixed"; next} {print}' "$WORK/none/output/param.par" > "$WORK/fixed.par"
c0p="$(run none_par "" "$WORK/fixed.par")"
cc="$(run c "$row_b" "$WORK/fixed.par")"

echo "  chi2: no prior $c0   (a) $ca   (b) $cb"
echo "  chi2 with param.par (energy_1 fixed): no prior $c0p   (c) $cc"
for v in "$c0" "$ca" "$cb" "$c0p" "$cc"; do
  if [ -z "$v" ]; then echo "  FAIL  a run produced no chi-squared"; exit 1; fi
done

bad=0
if awk -v a="$ca" -v z="$c0" 'BEGIN{d=(a-z)/z; if(d<0)d=-d; exit !(d<1e-12)}'; then
  echo "  ok    (a) prior at the parameter's value adds nothing"
else
  echo "  FAIL  (a) prior at the parameter's value changed chi2 by $(awk -v a="$ca" -v z="$c0" 'BEGIN{print a-z}')"; bad=1
fi
for case in b c; do
  v="$cb"; z="$c0"
  [ "$case" = c ] && { v="$cc"; z="$c0p"; }
  if awk -v a="$v" -v z="$z" 'BEGIN{d=(a-z)-1; if(d<0)d=-d; exit !(d<1e-6)}'; then
    echo "  ok    ($case) a 1-sigma offset adds 1"
  else
    echo "  FAIL  ($case) a 1-sigma offset added $(awk -v a="$v" -v z="$z" 'BEGIN{print a-z}'), expected 1"; bad=1
  fi
done
[ "$bad" -eq 0 ] && echo "  ok    nuisance prior evaluated on the parameter it names"

# (d) the writer reports the nuisance term (D6)
nu="$(grep -oE 'Total-Nuisance-Chi-Squared: [0-9.eE+-]+' "$WORK/b/output/chiSquared.out" | awk '{print $2}')"
if [ -n "$nu" ] && awk -v a="$nu" 'BEGIN{d=a-1; if(d<0)d=-d; exit !(d<1e-6)}'; then
  echo "  ok    (d) chiSquared.out: Total-Nuisance-Chi-Squared $nu (the 1-sigma offset)"
else
  echo "  FAIL  (d) chiSquared.out: Total-Nuisance-Chi-Squared '${nu}', expected 1"; bad=1
fi

# (e) the MCMC applies the same prior (the shared configuration, step 5)
d="$WORK/b"
# Only the set-up is needed (the priors are reported, and flushed, before sampling starts): stop after 30 s.
(cd "$d" && printf '6\n\n24\n100\n\n0.5\n1\nno\nyes\n' | timeout 30 "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 \
   | head -c 5000000 | tr '\r' '\n' > mcmc.log)
if grep -q "Nuisance priors from <parameterSettings>: 1 level/background parameter" "$d/mcmc.log"; then
  echo "  ok    (e) MCMC: the nuisance prior of <parameterSettings> is applied"
else
  echo "  FAIL  (e) MCMC: no nuisance prior applied"; grep -i "prior" "$d/mcmc.log" | sed -n '1,5p'; bad=1
fi
exit "$bad"
