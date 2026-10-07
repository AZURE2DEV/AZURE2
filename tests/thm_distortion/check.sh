#!/usr/bin/env bash
#
# Distortion factor R(E) of a THM experiment: `distortion=` on an
# `experiment[<name>]` line of the <thm> block (ThmDistortion.h;
# docs/source/theory/thm_implementation.rst, "Distortion factor R(E)").  The
# model of every segment is multiplied, before the folding, by
#   R(E) = rho(E)/rho(E_ref),  rho = |M|^2/|M_PW|^2 (dwpw) or |M|^2 (dw),
# M the zero-range prior-form DWBA transfer amplitude
# <chi(-)_sF phi_sx chi(+)_aA(beta r)> and M_PW its plane-wave limit.
#
# Models: tests/18O_p_a_thm (19F, two THM segments with free norms; its real
# reaction 2H(18O,a15N)n has a neutron spectator, a charged one is made up as
# 18O(3He,a15N)d at 115 MeV) and a one-segment copy of
# examples/c12c12_tumino2018 (12C(14N,a/p)d at 30 MeV) on points at chosen
# energies, without folding.
#
#   (a) distortion=none is byte-identical to the same line without the key;
#   (b) no distortion left: a neutral spectator with the a + A wave a plane
#       wave (distortion=optical opticalAA=plane) gives R = 1 at every point
#       (1e-6); with point Coulomb in a + A it stays within 12 % of 1;
#   (c) the model ratio on/off at the lowest and highest point equals the R
#       that thm_experiments.out evaluates there directly (grid interpolation,
#       1e-5); R(E_ref) = 1 at a point chosen as E_ref.  The amplitude itself
#       is checked against mpmath/scipy in tests/reference (ctest
#       thm_distortion);
#   (d) 12C+12C: the published factors that multiply the PWA S*, normalized at
#       2.664 MeV -- Mukhamedzhanov & Pang PRC 99 (2019) Fig. 10 (30 MeV,
#       Coulomb, FRESCO) with the default conventions (qf, whittaker, dwpw),
#       and Mukhamedzhanov arXiv:2609.04498 Fig. 9 (forward angles) with
#       distortionRatio=dw and E_sF 50 keV higher (Ebeam=30.11) -- within
#       the digitisation (rms <= 0.035 dex, max <= 0.06 dex);
#   (e) an optical potential with the nuclear part off (ten zeros) is
#       byte-identical to distortion=coulomb in the model; a real one is not;
#   (f) refusals, and the warning when the masses and field 32 disagree on
#       B(x+s);
#   (g) global optical potentials (ThmOptical.h): opticalAA=ancai06 equals
#       its ten numbers written out at E_aA (1e-6), kd03:extrapolate for
#       n + 19F is warned and reported at both data ends, refusals (outside
#       the range, a projectile the model does not describe, a bad name).
#
#   ./tests/thm_distortion/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
O18="$HERE/../18O_p_a_thm"
C12="$HERE/../../examples/c12c12_tumino2018"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/thm_distortion.XXXXXX")"
# The engine is a native binary: on MSYS2 give it a Windows path for files it reads.
WP="$WORK"; if command -v cygpath > /dev/null 2>&1; then WP="$(cygpath -m "$WORK")"; fi
trap 'rm -rf "$WORK"' EXIT
. "$HERE/../lib/guard.sh"
RUN="$(guard_command "${TEST_TIMEOUT:-300}")"
OUT="AZUREOut_aa=1_R=2.out"
KIN="beam=18O target=3He spectator=d Ebeam=115"
KINN="beam=18O target=d spectator=n Ebeam=54"
fail=0

# run NAME BLOCK [EDITS [PROJECT]] -- a "calculate with data" run of a copy of
# PROJECT (o18, the default, or c12) with BLOCK as its <thm> block ("" = none),
# the .azr first passed through the awk programs named in EDITS.  Status in
# NAME/status.
run() {
  local d="$WORK/$1" src azr
  if [ "${4:-o18}" = c12 ]; then src="$C12"; azr=c12c12_tumino2018.azr; else src="$O18"; azr=18O_p_a_thm.azr; fi
  rm -rf "$d"
  mkdir -p "$d/output" "$d/checks"
  cp -r "$src/data" "$d/"
  [ "${4:-o18}" = c12 ] && cp "$WORK/pts.dat" "$d/data/pts.dat"
  cp "$src/$azr" "$d/run.azr"
  local edit
  for edit in ${3:-}; do awk "${!edit}" "$d/run.azr" > "$d/tmp.azr" && mv "$d/tmp.azr" "$d/run.azr"; done
  [ -z "$2" ] || printf '\n<thm>\n%s\n</thm>\n' "$2" >> "$d/run.azr"
  (cd "$d" && printf '1\n\n\n7\n' | $RUN "$AZURE2_BIN" --no-gui --no-readline run.azr 2>&1 | head -c 1000000 > log;
   echo "${PIPESTATUS[1]}" > status)
}
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; fail=1; }
same() { [ -f "$1" ] && [ -f "$2" ] && [ "$(tr -d '\r' < "$1")" = "$(tr -d '\r' < "$2")" ]; }
ran() { [ "$(cat "$WORK/$1/status")" = 0 ] && [ -f "$WORK/$1/output/$OUT" ] || { bad "run $1 failed"; tail -5 "$WORK/$1/log" | sed 's/^/        /'; return 1; }; }
# "E ratio" per point: model(ON)/model(OFF), energies must agree.
ratios() {
  paste -d' ' <(tr -d '\r' < "$WORK/$1/output/$OUT" | awk 'NF > 4 { print $1, $4 }') \
              <(tr -d '\r' < "$WORK/$2/output/$OUT" | awk 'NF > 4 { print $1, $4 }') |
    awk '{ if ($1 != $3) { print "ENERGY", $1, $3; exit 1 } printf "%.12e %.12e\n", $1, $2 / $4 }'
}

# awk programs editing the .azr
NOFOLD='/<targetInt>/ { print; T = 1; next } /<\/targetInt>/ { T = 0 } T { next } { print }'
# examples/c12c12_tumino2018 reduced to one THM segment on data/pts.dat, no
# folding, its <thm> block dropped (the runs append their own).
C12MIN='/<segmentsData>/ { print; print "1  1  2  0  7  0  180  10  1  1  0  0  0  0  data/pts.dat  0  0"; S = 1; next }
/<\/segmentsData>/ { S = 0 } /<targetInt>/ { print; S = 1; next } /<\/targetInt>/ { S = 0 }
/<thm>/ { T = 1; next } /<\/thm>/ { T = 0; next } S || T { next } { print }'

# (a) ----------------------------------------------------------------------
echo "(a) distortion=none == no key"
run nokey "experiment[A] segments=1,2 $KIN"
run none "experiment[A] segments=1,2 $KIN distortion=none"
ran nokey && ran none
for f in chiSquared.out normalizations.out parameters.out thm_experiments.out "$OUT"; do
  same "$WORK/nokey/output/$f" "$WORK/none/output/$f" && ok "$f byte-identical" || bad "$f differs"
done

# (b) ----------------------------------------------------------------------
echo "(b) neutral spectator"
run n_off "experiment[A] segments=1,2 $KINN"
run n_plane "experiment[A] segments=1,2 $KINN distortion=optical opticalAA=plane"
run n_coul "experiment[A] segments=1,2 $KINN distortion=coulomb"
if ran n_off && ran n_plane && ran n_coul; then
  w="$(ratios n_plane n_off | awk '{ d = $2 - 1; if (d < 0) d = -d; if (d > w) w = d; n++ } END { printf "%d %.3e", n, w }')"
  awk -v w="${w#* }" -v n="${w% *}" 'BEGIN { exit !(n > 100 && w < 1e-6) }' \
    && ok "no distortion left (plane a + A, eta_sF = 0): R = 1 at ${w% *} points, worst |R - 1| = ${w#* }" \
    || bad "plane waves: R != 1 ($w)"
  r="$(ratios n_coul n_off | awk 'NR == 1 { lo = $2; hi = $2 } { if ($2 < lo) lo = $2; if ($2 > hi) hi = $2 } END { printf "%.4f %.4f", lo, hi }')"
  awk -v lo="${r% *}" -v hi="${r#* }" 'BEGIN { exit !(lo > 0.88 && hi < 1.12 && hi - lo > 0.01) }' \
    && ok "point Coulomb in d + 18O only (eta_aA = 0.73): R = ${r% *} to ${r#* } over the data" \
    || bad "neutral spectator, Coulomb a + A: R = $r"
fi

# (c) ----------------------------------------------------------------------
echo "(c) the model carries R(E)"
run c_off "experiment[A] segments=1,2 $KIN" "NOFOLD"
run c_on "experiment[A] segments=1,2 $KIN distortion=coulomb spectatorAngle=30" "NOFOLD"
if ran c_off && ran c_on; then
  ratios c_on c_off | sort -g > "$WORK/c.ratio"
  # distortion_point rows: E E_sF eta_sF theta |M|^2 |M_PW|^2 R lmax (lowest point, E_ref, highest point)
  tr -d '\r' < "$WORK/c_on/output/thm_experiments.out" | awk '$1 == "distortion_point" { print $2, $8 }' > "$WORK/c.rows"
  lo="$(head -1 "$WORK/c.ratio")"; hi="$(tail -1 "$WORK/c.ratio")"
  rlo="$(sed -n 1p "$WORK/c.rows")"; rhi="$(sed -n 3p "$WORK/c.rows")"
  awk -v a="$lo" -v b="$rlo" -v c="$hi" -v d="$rhi" 'BEGIN { split(a, x, " "); split(b, y, " "); split(c, u, " "); split(d, v, " ")
      e1 = x[2] / y[2] - 1; e2 = u[2] / v[2] - 1; if (e1 < 0) e1 = -e1; if (e2 < 0) e2 = -e2
      exit !(x[1] - y[1] < 1e-9 && y[1] - x[1] < 1e-9 && e1 < 1e-5 && e2 < 1e-5 && y[2] != 1 && v[2] != 1) }' \
    && ok "model ratio == R of thm_experiments.out at E = ${lo% *} (${lo#* }) and ${hi% *} (${hi#* }), lab 30 deg" \
    || bad "model ratio vs distortion_point: '$lo' '$rlo' / '$hi' '$rhi'"
  eref="$(awk 'NR == 10 { print $1 }' "$WORK/c.ratio")"
  run c_ref "experiment[A] segments=1,2 $KIN distortion=coulomb spectatorAngle=30 distortionRef=$eref" "NOFOLD"
  if ran c_ref; then
    r="$(ratios c_ref c_off | awk -v e="$eref" '$1 == e { print $2 }')"
    awk -v r="$r" 'BEGIN { d = r - 1; if (d < 0) d = -d; exit !(r != "" && d < 1e-6) }' \
      && ok "distortionRef=$eref: R = $r there" || bad "distortionRef=$eref: R = '$r'"
  fi
fi

# (d) ----------------------------------------------------------------------
# Digitized factors that multiply S* (1/R here, normalized at 2.664 MeV): E,
# 2019 Fig. 10 (E_14N = 30 MeV, pure Coulomb), 2026 Fig. 9 (forward angles).
DIG="0.80 2.2944e-03 3.7018e-03
0.90 2.6269e-03 4.2111e-03
1.00 3.0237e-03 4.8286e-03
1.20 4.1329e-03 6.5185e-03
1.40 5.9919e-03 9.1858e-03
1.60 9.1856e-03 1.3680e-02
1.80 1.5383e-02 2.1925e-02
2.00 2.8743e-02 3.8695e-02
2.20 6.2133e-02 7.8146e-02
2.40 1.6601e-01 1.9134e-01
2.55 4.3032e-01 4.5005e-01"
echo "$DIG" | awk '{ printf "%.4f 0 1.0 0.1\n", 2 * $1 }' > "$WORK/pts.dat"   # E_lab = 2 E_cm
echo "(d) 12C(14N,d): Mukhamedzhanov's published factors"
K12="beam=14N target=12C spectator=d"
run d_off "experiment[T] segments=1 $K12 Ebeam=30" "C12MIN" c12
run d_19 "experiment[T] segments=1 $K12 Ebeam=30 distortion=coulomb distortionRef=2.664" "C12MIN" c12
run d_26 "experiment[T] segments=1 $K12 Ebeam=30.11 distortion=coulomb distortionRatio=dw distortionRef=2.664" "C12MIN" c12
compare() {  # compare RUN COLUMN LABEL
  paste -d' ' <(ratios "$1" d_off) <(echo "$DIG") |
    awk -v c="$2" '{ if (($1 - $3) > 1e-6 || ($3 - $1) > 1e-6) { print "ENERGY", $1, $3; exit 1 }
      d = log(1 / $2) / log(10) - log($c) / log(10); s += d * d; n++; if (d < 0) d = -d; if (d > m) m = d }
      END { printf "%d %.4f %.4f\n", n, sqrt(s / n), m }' > "$WORK/$1.cmp"
  read -r n rms mx < "$WORK/$1.cmp"
  awk -v n="$n" -v r="$rms" -v m="$mx" 'BEGIN { exit !(n == 11 && r <= 0.035 && m <= 0.06) }' \
    && ok "$3: $n energies 0.8-2.55 MeV, rms $rms dex, max $mx dex" || bad "$3: $(cat "$WORK/$1.cmp")"
}
if ran d_off && ran d_19 && ran d_26; then
  compare d_19 4 "2019 Fig. 10 (FRESCO, Coulomb, 30 MeV) with qf, whittaker, dwpw"
  compare d_26 5 "2026 Fig. 9 (forward angles) with qf, whittaker, dw, E_sF + 50 keV"
fi

# (e) ----------------------------------------------------------------------
echo "(e) optical potential, nuclear part off == coulomb"
run e_coul "experiment[A] segments=1,2 $KIN distortion=coulomb"
run e_zero "experiment[A] segments=1,2 $KIN distortion=optical opticalAA=0,0,0,0,0,0,0,0,0,0 opticalSF=0,0,0,0,0,0,0,0,0,0"
run e_ws "experiment[A] segments=1,2 $KIN distortion=optical opticalAA=100,4.2,0.7,20,4.2,0.7,0,0,0,4.2 opticalSF=90,3.5,0.75,0,0,0,10,3.8,0.65,3.5"
if ran e_coul && ran e_zero && ran e_ws; then
  for f in chiSquared.out normalizations.out "$OUT"; do
    same "$WORK/e_coul/output/$f" "$WORK/e_zero/output/$f" && ok "$f byte-identical" || bad "$f differs"
  done
  same "$WORK/e_coul/output/$OUT" "$WORK/e_ws/output/$OUT" && bad "a Woods-Saxon potential changes nothing" \
    || ok "a Woods-Saxon potential changes the model"
fi

# (f) ----------------------------------------------------------------------
echo "(f) refused, and warned"
refuse() {  # refuse NAME MESSAGE-FRAGMENT BLOCK [EDITS PROJECT]
  run "$1" "$3" "${4:-}" "${5:-o18}"
  if [ "$(cat "$WORK/$1/status")" != 0 ] && grep -q "ERROR: <thm> experiment\[" "$WORK/$1/log" &&
     grep -q -- "$2" "$WORK/$1/log"; then
    ok "$1: $(grep -m1 'ERROR' "$WORK/$1/log" | tr -d '\r')"
  else
    bad "$1: not refused with '$2' (status $(cat "$WORK/$1/status"))"; tail -3 "$WORK/$1/log" | sed 's/^/        /'
  fi
}
printf '0.6 1.0\n0.7 2.0\n' > "$WORK/short.dat"
refuse no_kinematics "distortion=coulomb needs the kinematics" "experiment[A] segments=1,2 distortion=coulomb"
refuse bad_kind "expected none, coulomb, optical or table" "experiment[A] segments=1,2 $KIN distortion=dwba"
refuse optical_key "need distortion=optical" "experiment[A] segments=1,2 $KIN distortion=coulomb opticalAA=plane"
refuse optical_count "ten numbers" "experiment[A] segments=1,2 $KIN distortion=optical opticalSF=50,1.2,0.6"
refuse optical_radius "ten numbers" "experiment[A] segments=1,2 $KIN distortion=optical opticalAA=50,0,0.6,0,0,0,0,0,0,0"
refuse angle_value "expected qf, a lab angle" "experiment[A] segments=1,2 $KIN distortion=coulomb spectatorAngle=200"
refuse angle_alone "needs distortion=coulomb or distortion=optical" "experiment[A] segments=1,2 $KIN spectatorAngle=10"
refuse ratio_value "expected dwpw or dw" "experiment[A] segments=1,2 $KIN distortion=coulomb distortionRatio=pw"
refuse bound_value "expected whittaker or yukawa" "experiment[A] segments=1,2 $KIN distortion=coulomb boundState=hulthen"
refuse theta "expected all or thmin-thmax" "experiment[A] segments=1,2 $KIN theta=5"
refuse twice "is given twice" "experiment[A] segments=1,2 $KIN distortion=coulomb distortion=none"
refuse table_missing "cannot read" "experiment[A] segments=1,2 distortion=table:$WP/none_such.dat"
refuse table_short "beyond the table" "experiment[A] segments=1,2 distortion=table:$WP/short.dat"
refuse no_energy "the spectator has no energy left" "experiment[A] segments=1,2 beam=18O target=3He spectator=d Ebeam=40 distortion=coulomb"
# 12C(14N,d): above E = 2.3 MeV the deuteron cannot reach 75 deg in the lab.
refuse lab_reach "beyond the reach of the spectator" "experiment[T] segments=1 $K12 Ebeam=30 distortion=coulomb spectatorAngle=75" "C12MIN" c12
refuse ref_energy "distortionRef: at E = 20" "experiment[A] segments=1,2 $KIN distortion=coulomb distortionRef=20"
# B(x+s): masses 5.493 MeV (3He = d + p) against field 32 = 2.2246 MeV; the real d kinematics agree.
grep -q "WARNING: <thm> experiment\[A\]: B(x+s) from the masses of 3He = x + d is 5.49343 MeV, but the entrance pair 1 carries B = 2.22457 MeV (field 32" "$WORK/c_off/log" \
  && ok "B(x+s) warning: masses 5.49343 MeV vs field 32 2.22457 MeV" || bad "no B(x+s) warning for 3He"
grep -q "B(x+s) from the masses" "$WORK/n_off/log" && bad "B(x+s) warning for the d kinematics (they agree)" \
  || ok "no B(x+s) warning when masses and field 32 agree (d)"

# (g) ----------------------------------------------------------------------
echo "(g) global optical potentials (2H(18O,a15N)n: d + 18O, n + 19F)"
# ancai06 for d + 18O at E_d = 6.0424 MeV, the ten numbers written out (An & Cai 2006).
ANCAI="92.3094248915,3.01133408613,0.752021341567,1.47983816455,3.51976751944,0.592925860579,10.6451021248,3.64806901987,0.693485966441,3.41482603665"
run g_name "experiment[A] segments=1,2 $KINN distortion=optical opticalAA=ancai06 opticalSF=plane"
run g_ten "experiment[A] segments=1,2 $KINN distortion=optical opticalAA=$ANCAI opticalSF=plane"
if ran g_name && ran g_ten; then
  ratios g_name g_ten | awk '{ d = $2 - 1; if (d < 0) d = -d; if (d > m) m = d; n++ }
    END { printf "%d %.2e\n", n, m; exit !(n > 0 && m <= 1e-6) }' > "$WORK/g.cmp" \
    && ok "opticalAA=ancai06 == its ten numbers at E_aA: $(cat "$WORK/g.cmp") (points, max |ratio - 1|)" \
    || bad "ancai06 against its ten numbers: $(cat "$WORK/g.cmp")"
  grep -q "WARNING" "$WORK/g_name/log" && bad "a warning for ancai06 inside its range" || ok "no warning inside the range"
fi
run g_extra "experiment[A] segments=1,2 $KINN distortion=optical opticalAA=ancai06 opticalSF=kd03:extrapolate"
if ran g_extra; then
  grep -q "WARNING: <thm> experiment\[A\]: opticalSF=kd03 (Koning & Delaroche, NPA 713 (2003) 231 (global)) is outside its validity range for n + 19F: target A = 19 (valid 24-209); extrapolated" "$WORK/g_extra/log" \
    && ok "kd03:extrapolate for n + 19F: warned" || bad "no extrapolation warning for kd03"
  grep -q "s + F: kd03:extrapolate (n on 19F) Woods-Saxon V=.* at E_lab = .* MeV (E = .* MeV) to V=" "$WORK/g_extra/log" \
    && ok "the s + F potential is reported at both ends of the data" || bad "no s + F potential at both ends"
fi
refuse g_range "Write kd03:extrapolate" "experiment[A] segments=1,2 $KINN distortion=optical opticalSF=kd03"
refuse g_species "which it does not describe" "experiment[A] segments=1,2 $KINN distortion=optical opticalAA=kd03"
refuse g_mass "outside its validity range for d + 18O: target A = 18 (valid 27-238)" "experiment[A] segments=1,2 $KINN distortion=optical opticalAA=daehnick80"
refuse g_name_bad "a global optical potential (ancai06, daehnick80" "experiment[A] segments=1,2 $KINN distortion=optical opticalAA=ancai"

echo
if [ "$fail" -eq 0 ]; then echo "PASS: THM distortion factor"; else echo "FAIL"; exit 1; fi
