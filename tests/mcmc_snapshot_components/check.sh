#!/bin/bash
# The MCMC snapshot (written every 10 steps: AZUREOut_*.out, param.mcmc) carries each segment's model value -- the
# value the likelihood compares with the data -- also for an E1-only or E2-only capture segment (isDiff 5, 6), where
# it is the E1 or E2 part of the angle-integrated cross section, not the total (docs/upgrade/data_model.md S4).
# Model: tests/mcmc_physical_norm/13N_fast.azr with three copies of a 12C(p,g) capture segment (total, E1 only,
# E2 only; tests/13N/data/skowronski_luna_bgo.dat); 4 walkers x 100 steps, physical parameters.  Checks:
#   (a) each of the three capture blocks of the snapshot equals a mode-1 calculation from param.mcmc (the sampled
#       point, 7 digits) to 1e-5 relative;
#   (b) the E2-only block is not the total (they differ by orders of magnitude here: E2 is ~1e-4 of the total).
# Fails if the snapshot writes the base point's leaf cross section for an E1/E2-only segment (the total): shown by
# injecting that into the snapshot once.  About 15 s.
set -u
export LC_ALL=C OMP_NUM_THREADS=2
[ $# -ge 1 ] || { echo "usage: check.sh path/to/AZURE2"; exit 2; }
AZURE2_BIN="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fail=0
ok()  { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; fail=1; }

mkdir -p "$WORK/output" "$WORK/checks"
cp -r "$HERE/../13N/data" "$WORK/data"
# three capture data sets on one file: the total, E1 only, E2 only
export RECS="$(for m in total E1 E2; do
  printf 'dataset bgo_%s\n    reaction       12C(p,g)13N\n    quantity       sigma\n    unit           b\n' "$m"
  [ "$m" = total ] || printf '    multipole      %s\n' "$m"
  printf '    file           data/skowronski_luna_bgo.dat\n    columns        E_lab[MeV] _ value error\n'
  printf '    energy_window  0.05 0.5 MeV\n    norm           1 fixed\n'
done)"
awk '/^<\/data>$/ { print ENVIRON["RECS"] } { print }' "$HERE/../mcmc_physical_norm/13N_fast.azr" > "$WORK/t.azr"
cd "$WORK"
# mode 6: param file (blank), EC file (blank), walkers, steps, spread %, energy spread (default), threads, RWA?, overwrite?
"$AZURE2_BIN" mcmc --walkers=4 --steps=100 --parameter-spread=0.5 --sampler-threads=1 --basis=physical --overwrite-samples --no-gui t.azr > mcmc.log 2>&1
[ -s output/param.mcmc ] && [ -s "output/AZUREOut_aa=1_R=2.out" ] || { bad "no snapshot written"; tail -5 mcmc.log; exit 1; }
cp "output/AZUREOut_aa=1_R=2.out" snapshot.out
cp output/param.mcmc sampled.par
mkdir -p calc/output calc/checks && cp t.azr calc/ && cp -r data calc/
( cd calc && "$AZURE2_BIN" calculate --parameter-file=../sampled.par --no-gui t.azr > calc.log 2>&1 )
grep -q "Reading User Parameter File" calc/calc.log || { bad "the calculation did not read param.mcmc"; tail -5 calc/calc.log; exit 1; }
# the blocks of a file (blank-line separated) -> "block point value" with the calculated cross section (column 4)
blocks() { awk 'NF == 0 { if (n) b++; n = 0; next } { n++; print b + 0, n, $4 }' "$1"; }
blocks snapshot.out > snap.txt
blocks "calc/output/AZUREOut_aa=1_R=2.out" > calc.txt
nb="$(awk 'END { print $1 + 1 }' snap.txt)"
[ "$nb" = 3 ] && [ "$(wc -l < snap.txt)" = "$(wc -l < calc.txt)" ] || { bad "expected 3 capture blocks of equal length ($nb)"; exit 1; }
for b in 0 1 2; do
  name=$(echo "total E1-only E2-only" | cut -d' ' -f$((b + 1)))
  rel="$(paste snap.txt calc.txt | awk -v b="$b" '$1 == b { d = $3 - $6; if (d < 0) d = -d; m = ($6 < 0 ? -$6 : $6);
         r = (m > 0 ? d / m : d); if (r > w) w = r } END { printf "%.2e", w }')"
  if awk -v r="$rel" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && r + 0 <= 1e-5) }'; then ok "(a) $name: snapshot = calculation at param.mcmc (max rel $rel)"
  else bad "(a) $name: snapshot differs from the calculation at param.mcmc (max rel $rel)"; fi
done
ratio="$(awk '$1 == 0 { t[$2] = $3 } $1 == 2 { e[$2] = $3 } END { m = 0; for (k in e) if (t[k] != 0) { r = e[k] / t[k];
          if (r < 0) r = -r; if (r > m) m = r } printf "%.2e", m }' snap.txt)"
if awk -v r="$ratio" 'BEGIN { exit !(r ~ /^ *[-+]?\.?[0-9]/ && r + 0 < 0.5) }'; then ok "(b) the E2-only block is not the total (max E2/total $ratio)"
else bad "(b) the E2-only block looks like the total (max E2/total $ratio)"; fi
[ $fail -eq 0 ] && echo "mcmc_snapshot_components: PASS" || echo "mcmc_snapshot_components: FAIL"
exit $fail
