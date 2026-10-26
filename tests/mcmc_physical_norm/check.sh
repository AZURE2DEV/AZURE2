#!/bin/bash
# Physical-parameter MCMC (mode 6, "Use Reduced Width Amplitudes: no") must apply the sampled
# normalizations.  Before 2026-10-26 the likelihood never copied the norms out of the physical
# vector: the sampler moved them under their prior while every evaluation used the starting norm
# (norm posterior = prior, R-matrix posterior at frozen norms), and param.mcmc recorded the start.
# The chain itself is random (no seed control); the checks are relations that hold for any chain:
#   (a) param.mcmc holds the norm of the best walker of the last step;
#   (b) the log-likelihood the sampler recorded for that walker is -chi2/2 of its parameters
#       with the SAMPLED norm (a mode-1 calculation from param.mcmc with that norm), not with the
#       start norm.
# The 2026-10-26 build passes both; dev dfc64e1 fails both (param.mcmc norm = 1, chi2 off by the
# norm's effect).  Model: 13N_fast.azr (make_model.py), 4 walkers x 100 steps, ~15 s.
set -u
[ $# -ge 1 ] || { echo "usage: check.sh path/to/AZURE2"; exit 2; }
AZURE2_BIN="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail=0
ok()  { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; fail=1; }

mkdir -p "$WORK/output" "$WORK/checks"
cp "$HERE/13N_fast.azr" "$WORK/t.azr"
cp -r "$HERE/../13N/data" "$WORK/data"
cd "$WORK"
# mode 6: param file (blank), EC file (blank), walkers, steps, spread %, energy spread (default), threads, RWA?, overwrite?
printf '6\n\n\n4\n100\n0.5\n\n1\nno\nyes\n' | "$AZURE2_BIN" --no-gui --no-readline t.azr > mcmc.log 2>&1
grep -q "Parameter mode: Physical" mcmc.log || { bad "the run was not in physical-parameter mode"; tail -5 mcmc.log; exit 1; }
[ -s output/samples.mcmc ] && [ -s output/param.mcmc ] || { bad "no samples.mcmc / param.mcmc written"; tail -5 mcmc.log; exit 1; }
nfree="$(awk '/Number of free parameters:/ { print $NF }' mcmc.log)"
# the free vector is level parameters then the one free norm: the norm is the last column
read -r step walker logL norm < <(awk -F, -v n="$nfree" 'NR > 1 { s = $1 } END { print s }' output/samples.mcmc | {
  read -r last
  awk -F, -v last="$last" 'NR > 1 && $1 == last && ($3 > best || best == "") { best = $3; line = $0 }
       END { n = split(line, f, ","); print f[1], f[2], f[4], f[n] }' output/samples.mcmc; })
pnorm="$(awk '$1 ~ /_norm$/ { print $2 }' output/param.mcmc)"
echo "  last step $step, best walker $walker: sampled norm $norm, param.mcmc norm $pnorm, logL $logL"
rel() { awk -v a="$1" -v b="$2" 'BEGIN { d = a - b; if (d < 0) d = -d; m = (a < 0 ? -a : a); if (m == 0) m = 1; printf "%.2e", d / m }'; }
d="$(rel "$pnorm" "$norm")"
if awk -v d="$d" 'BEGIN { exit !(d + 0 <= 1e-6) }'; then ok "(a) param.mcmc carries the sampled norm (rel $d)"
else bad "(a) param.mcmc norm $pnorm is not the sampled norm $norm (rel $d)"; fi
# (b) re-evaluate the sampled point two ways: param.mcmc (the full RWA set, 7 digits) with the norm set
# to the sampled value, and with the norm at its start value 1.  The likelihood the sampler recorded must
# be the first (to the 7-digit rounding of the file, ~1e-5 of chi2) and not the second.
norm0="$(awk '/<segmentsData>/ { on = 1; next } /<\/segmentsData>/ { on = 0 } on && $15 ~ /meyer_84/ { print $8 }' t.azr)"
awk -v n="$norm" '$1 ~ /_norm$/ { $2 = n } { printf "%20s %20s %20s\n", $1, $2, $3 }' output/param.mcmc > sampled.par
awk -v n="$norm0" '$1 ~ /_norm$/ { $2 = n } { printf "%20s %20s %20s\n", $1, $2, $3 }' output/param.mcmc > start.par
chi2_at() {  # chi2_at <parfile> -> data chi2 of a mode-1 calculation from that parameter file
  local d="calc_$(basename "$1" .par)"
  mkdir -p "$d/output" "$d/checks" && cp t.azr "$d/" && cp -r data "$d/"
  ( cd "$d" && printf '1\n../%s\n\n' "$1" | "$AZURE2_BIN" --no-gui --no-readline t.azr > calc.log 2>&1 )
  grep -q "Reading User Parameter File" "$d/calc.log" || { echo "nan"; return; }
  awk '/Total-Chi-Squared:/ { for (i = 1; i <= NF; i++) if ($i == "Total-Chi-Squared:") print $(i + 1) }' "$d/output/chiSquared.out" | tail -1
}
cs="$(chi2_at sampled.par)"; c0="$(chi2_at start.par)"
m2="$(awk -v l="$logL" 'BEGIN { printf "%.10g", -2 * l }')"
echo "  -2 logL = $m2; mode-1 chi2 with the sampled norm $norm: $cs, with the start norm $norm0: $c0"
ds="$(rel "$cs" "$m2")"; d0="$(rel "$c0" "$m2")"
if awk -v a="$ds" -v b="$d0" 'BEGIN { exit !(a + 0 <= 1e-4 && a + 0 < b + 0) }'; then
  ok "(b) the recorded log-likelihood is the chi2 of the sampled parameters (rel $ds; at the start norm $d0)"
else
  bad "(b) -2 logL $m2 vs chi2 at the sampled norm $cs (rel $ds) and at the start norm $c0 (rel $d0): the likelihood did not use the sampled norm"
fi
[ $fail -eq 0 ] && echo "mcmc_physical_norm: PASS" || echo "mcmc_physical_norm: FAIL"
exit $fail
