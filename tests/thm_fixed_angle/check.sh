#!/usr/bin/env bash
#
# Fixed-angle (differential) HOES observable of a THM experiment: `theta=` on
# an `experiment[<name>]` line of the <thm> block (ThmAngular.h;
# docs/source/theory/thm_implementation.rst, "Fixed-angle observable").  With
# theta=thmin-thmax the model of every segment of the experiment is
#   <dsigma/dOmega> = (1/pi) sum_L b_L <P_L(cos theta)>_window,
# the Blatt-Biedenharn sum over the HOES partial amplitudes of all J groups
# (m_l = 0 along p_xA), averaged over the window in cos theta; its 4 pi
# integral is the angle-integrated HOES cross section.
#
# Models: tests/7Li_p_a (7Li+p -> a+a: two entrance l per channel spin,
# identical-boson exit, 30 keV folding), tests/17O (17O(n,a)14C: J^pi groups
# of both parities -> odd L), tests/6Li_d (6Li(d,a)a), tests/18O_p_a_thm
# (one 1/2+ group: isotropic; two segments, the 2H(18O,a15N)n kinematics).
#
#   (a) theta=all changes nothing: output files byte-identical to the same
#       line without the key, and the model to no <thm> block at all;
#   (b) a 0-180 window is the angle-integrated model / 4 pi (1.2e-10: the
#       11 significant digits of the output file) and gives the same chi2;
#       for the one-group 18O model any window does, also combined with a
#       spectator-momentum window and a linear background;
#   (c) theta -> 0: windows 0-t approach the single angle 0-0 as t^2; the
#       0-deg observable is not the entranceL=coherent one;
#   (d) (tests/reference thm_fixed_angle: the formula against an independent
#       M-sum -- not repeated here);
#   (e) identical exit particles: 30-60 == 120-150 (7Li+p, 6Li+d); 17O's
#       distribution is not symmetric;
#   (f) refusals;
#   and the effect size of the Tumino 2006 window theta_cm = 50-70 deg on
#   7Li(p,a) (parameters not refitted).
#
#   ./tests/thm_fixed_angle/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_fixed_angle.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
. "$HERE/../lib/check_common.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-600}")"
KIN18="beam=18O target=d spectator=n Ebeam=54"
FOURPI="12.566370614359172"
fail=0

# run PROJECT NAME BLOCK -- a "calculate with data" run of a copy of
# tests/PROJECT with BLOCK as its <thm> block ("" = none).  Status in NAME/status.
run() {
  local d="$WORK/$2"
  rm -rf "$d"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$HERE/../$1/data" "$d/"
  cp "$HERE/../$1/$1.azr" "$d/run.azr"
  [ -z "$3" ] || printf '\n<thm>\n%s\n</thm>\n' "$3" >> "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 | head -c 1000000 > log;
   echo "${PIPESTATUS[1]}" > status)
}
out() { ls "$WORK/$1/output/" | grep '^AZUREOut_aa=' | head -1; }
ran() { [ "$(cat "$WORK/$1/status")" = 0 ] && [ -n "$(out "$1")" ] || { bad "run $1 failed"; tail -5 "$WORK/$1/log" | sed 's/^/        /'; return 1; }; }
model() { tr -d '\r' < "$WORK/$1/output/$(out "$1")" | awk 'NF > 4 { print $1, $4 }'; }
chi2() { tr -d '\r' < "$WORK/$1/log" | awk '/Total Chi-Squared:/ { v = $NF } END { print v }'; }
# worst |model(A) * F / model(B) - 1| (energies must agree): "n worst"
worst() {
  paste -d' ' <(model "$1") <(model "$2") |
    awk -v f="${3:-1}" '{ if ($1 != $3) { e = 1; exit } d = $2 * f / $4 - 1; if (d < 0) d = -d; if (d > w) w = d; n++ }
         END { if (e) print "0 ENERGY-MISMATCH"; else printf "%d %.3e\n", n, w }'
}
# check_worst A B FACTOR TOL LABEL
check_worst() {
  local n w
  read -r n w <<< "$(worst "$1" "$2" "$3")"
  awk -v w="$w" -v n="$n" -v t="$4" 'BEGIN { exit !(n > 10 && w <= t) }' && ok "$5: $n points, worst rel $w" \
    || bad "$5: $n points, worst rel $w (tol $4)"
}
same_chi2() {  # same_chi2 A B LABEL
  local a b
  a="$(chi2 "$1")"; b="$(chi2 "$2")"
  awk -v a="$a" -v b="$b" 'BEGIN { d = a / b - 1; if (d < 0) d = -d; exit !(a != "" && d < 1e-9) }' \
    && ok "$3: chi2 $a == $b" || bad "$3: chi2 $a vs $b"
}

# (a) ----------------------------------------------------------------------
echo "(a) theta=all == no key"
run 7Li_p_a none ""
run 7Li_p_a nokey "experiment[A] segments=1"
run 7Li_p_a all "experiment[A] segments=1 theta=all"
if ran none && ran nokey && ran all; then
  for f in chiSquared.out normalizations.out parameters.out thm_experiments.out "$(out all)"; do
    same "$WORK/nokey/output/$f" "$WORK/all/output/$f" && ok "$f byte-identical" || bad "$f differs"
  done
  for f in chiSquared.out "$(out all)"; do
    same "$WORK/none/output/$f" "$WORK/all/output/$f" && ok "no <thm> block: $f byte-identical" || bad "no <thm> block: $f differs"
  done
fi

# (b) ----------------------------------------------------------------------
echo "(b) a 0-180 window is the angle-integrated model / 4 pi"
run 7Li_p_a full "experiment[A] segments=1 theta=0-180"
if ran full; then
  check_worst full nokey "$FOURPI" 1.2e-10 "7Li(p,a): 4 pi <dsigma/dOmega>_{0-180} vs angle-integrated"
  same_chi2 full nokey "7Li(p,a) 0-180"
fi
run 17O o17_all ""
run 17O o17_full "experiment[A] segments=1 theta=0-180"
if ran o17_all && ran o17_full; then
  check_worst o17_full o17_all "$FOURPI" 1.2e-10 "17O(n,a): 4 pi <dsigma/dOmega>_{0-180} vs angle-integrated"
  same_chi2 o17_full o17_all "17O(n,a) 0-180"
fi
echo "    18O(p,a): one 1/2+ group, one entrance and one exit wave -- isotropic"
run 18O_p_a_thm o18_ps "experiment[A] segments=1,2 $KIN18 ps=hulthen:0-40 background=linear"
run 18O_p_a_thm o18_ps_theta "experiment[A] segments=1,2 $KIN18 ps=hulthen:0-40 background=linear theta=50-70"
if ran o18_ps && ran o18_ps_theta; then
  check_worst o18_ps_theta o18_ps "$FOURPI" 1.2e-10 "18O(p,a) theta=50-70 with ps window and background: 4 pi <.> vs angle-integrated"
  same_chi2 o18_ps_theta o18_ps "18O(p,a) theta=50-70, ps window, linear background"
fi

# (c) ----------------------------------------------------------------------
echo "(c) theta -> 0"
run 7Li_p_a p0 "experiment[A] segments=1 theta=0-0"
run 7Li_p_a p05 "experiment[A] segments=1 theta=0-0.5"
run 7Li_p_a p1 "experiment[A] segments=1 theta=0-1"
if ran p0 && ran p05 && ran p1; then
  # (m(0-1) - m(0)) / (m(0-0.5) - m(0)) -> 4 (the window mean is m(0) + c t^2)
  read -r n q lo hi <<< "$(paste -d' ' <(model p0) <(model p05) <(model p1) |
    awk '{ d1 = $4 - $2; d2 = $6 - $2; if (d1 * d1 < 1e-14 * $2 * $2) next; q = d2 / d1; s += q; n++;
           if (n == 1 || q < lo) lo = q; if (n == 1 || q > hi) hi = q } END { printf "%d %.6f %.6f %.6f\n", n, s / n, lo, hi }')"
  awk -v lo="$lo" -v hi="$hi" -v n="$n" 'BEGIN { exit !(n > 10 && lo > 3.99 && hi < 4.01) }' \
    && ok "windows 0-1 and 0-0.5 deg approach 0-0 as t^2: ratio of differences $lo..$hi over $n points" \
    || bad "not t^2: ratio $lo..$hi over $n points"
  run 7Li_p_a coherent "entranceL=coherent"
  if ran coherent; then
    read -r n w <<< "$(worst p0 coherent)"
    ok "(size) theta = 0 vs entranceL=coherent (angle-integrated): chi2 $(chi2 p0) vs $(chi2 coherent) -- different observables"
  fi
fi

# (e) ----------------------------------------------------------------------
echo "(e) identical exit particles: symmetric about 90 deg"
run 7Li_p_a f3060 "experiment[A] segments=1 theta=30-60"
run 7Li_p_a b120150 "experiment[A] segments=1 theta=120-150"
if ran f3060 && ran b120150; then
  check_worst f3060 b120150 1 1e-12 "7Li(p,a)a: 30-60 == 120-150"
  read -r n w <<< "$(worst f3060 nokey "$FOURPI")"
  ok "(size) 7Li(p,a) 30-60 vs angle-integrated: worst rel $w"
fi
run 6Li_d li6_f "experiment[A] segments=1 theta=10-40"
run 6Li_d li6_b "experiment[A] segments=1 theta=140-170"
if ran li6_f && ran li6_b; then
  check_worst li6_f li6_b 1 1e-12 "6Li(d,a)a: 10-40 == 140-170"
fi
run 17O o17_f "experiment[A] segments=1 theta=20-60"
run 17O o17_b "experiment[A] segments=1 theta=120-160"
if ran o17_f && ran o17_b; then
  read -r n w <<< "$(worst o17_f o17_b)"
  counted "$n" && awk -v w="$w" 'BEGIN { exit !(w > 1e-3) }' && ok "17O(n,a)14C (1-, 2+, 3-, 5- interfere): 20-60 vs 120-160 differ by up to $w" \
    || bad "17O(n,a): forward and backward windows agree ($w): no odd L"
fi

# (f) ----------------------------------------------------------------------
echo "(f) refused"
refuse() {  # refuse NAME MESSAGE-FRAGMENT BLOCK
  run 7Li_p_a "$1" "$3"
  if [ "$(cat "$WORK/$1/status")" != 0 ] && grep -q "ERROR: <thm> experiment\[" "$WORK/$1/log" &&
     grep -q -- "$2" "$WORK/$1/log"; then
    ok "$1: $(grep -m1 'ERROR' "$WORK/$1/log" | tr -d '\r')"
  else
    bad "$1: not refused with '$2' (status $(cat "$WORK/$1/status"))"; tail -3 "$WORK/$1/log" | sed 's/^/        /'
  fi
}
USAGE="expected all or thmin-thmax"
refuse word "$USAGE" "experiment[A] segments=1 theta=forward"
refuse single "$USAGE" "experiment[A] segments=1 theta=60"
refuse reversed "$USAGE" "experiment[A] segments=1 theta=70-50"
refuse negative "$USAGE" "experiment[A] segments=1 theta=-10-20"
refuse beyond "$USAGE" "experiment[A] segments=1 theta=170-190"
refuse empty_end "$USAGE" "experiment[A] segments=1 theta=50-"
refuse twice "given twice" "experiment[A] segments=1 theta=50-70
experiment[A] theta=all"
refuse coherent "entranceL=coherent" "entranceL=coherent
experiment[A] segments=1 theta=50-70"

# effect size ----------------------------------------------------------------
echo "(size) Tumino et al. EPJA 27 (2006) 243: theta_cm = 50-70 deg, same parameters"
run 7Li_p_a tumino "experiment[A] segments=1 theta=50-70"
if ran tumino; then
  # shape: the two models scaled to the same mean over the points
  paste -d' ' <(model tumino) <(model nokey) |
    awk '{ a[NR] = $2; b[NR] = $4; sa += $2; sb += $4; n++ }
         END { for (i = 1; i <= n; i++) { r = (a[i] / sa) / (b[i] / sb) - 1; s2 += r * r; if (r < 0) r = -r; if (r > m) m = r }
               printf "  (size) shape 50-70 vs angle-integrated (equal means): rms %.3f, max %.3f over %d points\n", sqrt(s2 / n), m, n }'
  echo "  (size) chi2 with the profiled norm: $(chi2 tumino) (50-70) vs $(chi2 nokey) (angle-integrated)"
  grep -q "Angular window: theta_cm = 50-70 deg" "$WORK/tumino/log" && ok "the log states the window" || bad "no window line in the log"
fi

echo
if [ "$fail" -eq 0 ]; then echo "PASS: THM fixed-angle observable"; else echo "FAIL"; exit 1; fi
