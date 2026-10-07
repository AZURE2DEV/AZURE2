#!/usr/bin/env bash
#
# Identical particles in the ENTRANCE channel of a reaction carry (1+delta_12).
#
# Three one-level projects (projects/*.azr), each evaluated exactly at the level
# energy E_R, where the one-level R-matrix cross section is, in any formalism
# whose boundary condition is S_c(E_R) (the default, and --use-brune),
#
#   sigma_ab(E_R) = (pi/k^2) g_J (1+delta_12) 4 Gamma_a Gamma_b / Gamma^2 ,
#   g_J = (2J+1) / [(2i_1+1)(2i_2+1)] ,
#
# and the observed partial widths of the .azr enter only through their ratio.
# Checked against that expression (awk, AZURE2's own constants):
#
#   c12c12  12C+12C -> a+20Ne, 2+ (i=0):  x2.  dsigma/dOmega = sigma 5 P_2^2 / 4pi,
#           i.e. the differential is x2 too.  Elastic 12C+12C: x2 (7c34992).
#   dd_pt   d+d -> p+t, 2+ (i=1, g = 5/9): x2; dsigma/dOmega isotropic = sigma/4pi.
#   be8     p+7Li -> a+a: NO factor (identical EXIT pair: one count per event);
#           a+a -> p+7Li: x2; and the two obey reciprocity
#           8 k_p^2 sigma(p7Li->aa) = 1 k_a^2 sigma(aa->p7Li) / (1+delta_aa),
#           the reciprocity theorem with the (2i_1+1)(2i_2+1) weights.
#
# Before the fix every "x2" above came out x1 except the elastic one.
#
#   ./tests/identical_entrance_reaction/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/identical_entrance.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if command -v timeout >/dev/null 2>&1; then RUN="timeout ${TEST_TIMEOUT:-600}"; else RUN=""; fi
TOL=1e-6
fail=0

# run NAME PROJECT [options] -> project evaluated in $WORK/NAME
run() {
  local d="$WORK/$1" p="$2"
  shift 2
  mkdir -p "$d/output" "$d/checks"
  cp -r "$HERE/projects/data" "$d/"
  cp "$HERE/projects/$p.azr" "$d/"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline "$@" "$p.azr" > log 2>&1)
}
# row FILE N COL -> column COL of the N-th non-blank line
# (\r stripped first: Windows writes CRLF, and a bare "\r" line has NF = 1)
row() { awk -v n="$2" -v c="$3" '{ sub(/\r$/, "") } NF { if (++i == n) { print $c; exit } }' "$1"; }

# peak M1 M2 ECM G SYM GA GB -> (pi/k^2)/100 g sym 4 GA GB/(GA+GB)^2 in barn
peak() {
  awk -v m1="$1" -v m2="$2" -v e="$3" -v g="$4" -v sym="$5" -v ga="$6" -v gb="$7" 'BEGIN {
    pi = 3.141592650; hbarc = 197.32696310; uconv = 931.4940880
    mu = m1 * m2 / (m1 + m2)
    printf "%.12e\n", pi * hbarc * hbarc / (2 * mu * uconv * e) / 100 * g * sym * 4 * ga * gb / ((ga + gb) * (ga + gb))
  }'
}
expect() {  # expect LABEL ACTUAL EXPECTED
  if [ -z "$2" ]; then echo "  FAIL  $1: no value"; fail=1; return; fi
  if awk -v a="$2" -v b="$3" -v t="$TOL" 'BEGIN{d=(a-b)/b; if(d<0)d=-d; exit !(d<=t)}'; then
    printf '  ok    %-44s %s (analytic %s)\n' "$1:" "$2" "$3"
  else
    printf '  FAIL  %-44s %s, analytic %s (ratio %s)\n' "$1:" "$2" "$3" "$(awk -v a="$2" -v b="$3" 'BEGIN{printf "%.6f", a/b}')"
    fail=1
  fi
}

# ---- 12C+12C (default formalism and Brune) -----------------------------------
for mode in default brune; do
  opt=""; [ "$mode" = brune ] && opt="--use-brune"
  run "c12_$mode" c12c12 $opt
  o="$WORK/c12_$mode/output"
  e=$(row "$o/AZUREOut_aa=1_R=2.out" 1 1)
  s=$(row "$o/AZUREOut_aa=1_R=2.out" 1 4)
  sr=$(peak 12.0 12.0 "$e" 5 2 1.5e-6 80000)
  expect "12C(12C,a0)20Ne sigma ($mode)" "$s" "$sr"
  for n in 2 3 4 5; do
    th=$(row "$o/AZUREOut_aa=1_R=2.out" $n 3)
    ds=$(row "$o/AZUREOut_aa=1_R=2.out" $n 4)
    x=$(awk -v s="$sr" -v t="$th" 'BEGIN{c=cos(t*3.141592650/180); p=0.5*(3*c*c-1); printf "%.12e", s*5*p*p/(4*3.141592650)}')
    expect "  dsigma/dOmega at $(printf %.1f "$th") deg ($mode)" "$ds" "$x"
  done
  # Elastic: |1 - U|^2 with U = exp(2i phi)(1 - 2 Gamma_c/Gamma) at E_R, so the
  # hard-sphere phase enters; phi = -atan(F_2/G_2) at E = 2 MeV, rho = k 6.4 fm
  # (mpmath) -- the one pinned input.  AZURE2's Coulomb functions give phi to
  # ~1e-6 relative here, hence the looser tolerance.
  el=$(peak 12.0 12.0 "$(row "$o/AZUREOut_aa=1_R=1.out" 1 1)" 5 2 1 1 |
       awk -v phi=4.4955151025e-12 '{x = 1.5e-6 / (1.5e-6 + 80000); a = 1 - 2 * x
         re = 1 - a * cos(2 * phi); im = a * sin(2 * phi)
         printf "%.12e", $1 * (re * re + im * im)}')
  TOL=1e-4 expect "12C+12C elastic sigma ($mode)" "$(row "$o/AZUREOut_aa=1_R=1.out" 1 4)" "$el"
done

# ---- d+d -> p+t: spin-1 identical pair --------------------------------------
run dd dd_pt
o="$WORK/dd/output/AZUREOut_aa=1_R=2.out"
s=$(row "$o" 1 4)
sr=$(peak 2.01410178 2.01410178 "$(row "$o" 1 1)" 0.5555555555555556 2 100000 200000)
expect "d(d,p)t sigma, g=5/9" "$s" "$sr"
for n in 2 3 4; do
  expect "  dsigma/dOmega at $(printf %.1f "$(row "$o" $n 3)") deg" "$(row "$o" $n 4)" \
    "$(awk -v s="$sr" 'BEGIN{printf "%.12e", s/(4*3.141592650)}')"
done

# ---- p+7Li <-> a+a: identical exit pair (no factor) and reciprocity -----------
run be8 be8
f="$WORK/be8/output/AZUREOut_aa=1_R=2.out"
r="$WORK/be8/output/AZUREOut_aa=2_R=1.out"
ep=$(row "$f" 1 1); sf=$(row "$f" 1 4)
ea=$(row "$r" 1 1); sb=$(row "$r" 1 4)
expect "7Li(p,a)a sigma, per event (no factor)" "$sf" "$(peak 1.00782503 7.01600344 "$ep" 0.625 1 50000 100000)"
expect "4He(a,p)7Li sigma, x(1+delta)" "$sb" "$(peak 4.00260325 4.00260325 "$ea" 5 2 50000 100000)"
# w k^2 sigma / (1+delta_entrance) is the same both ways, w = (2i_1+1)(2i_2+1)
# (8 for p+7Li, 1 for a+a); k^2 is proportional to mu E.
lhs=$(awk -v s="$sf" -v e="$ep" 'BEGIN{mu=1.00782503*7.01600344/(1.00782503+7.01600344); printf "%.12e", 8*mu*e*s}')
rhs=$(awk -v s="$sb" -v e="$ea" 'BEGIN{mu=4.00260325/2; printf "%.12e", 1*mu*e*s/2}')
TOL=1e-5 expect "reciprocity w k^2 sigma/(1+delta)" "$lhs" "$rhs"

if [ "$fail" -eq 0 ]; then
  echo "  ok    identical entrance pairs carry (1+delta_12); identical exit pairs count events"
  exit 0
fi
exit 1
