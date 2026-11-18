#!/usr/bin/env bash
#
# External capture to a final state that is not the first level of its J^pi group, in the formal
# parametrization (--no-brune).
#
# The external-capture terms need the final state's OBSERVED energy and its on-resonance (Brune) reduced
# width amplitude, which is what the ANC is made of.  In the formal parametrization the working energy
# and amplitude of a level are those only for the level the boundary conditions were taken from, the first
# of its J^pi group.  Before the fix the formal values were used for every final state, so the capture
# cross section depended on the ORDER of the <levels> lines: with the 1/2- ground state of 13N listed
# after the 1/2- background pole (models/A_swap.azr) the --no-brune 12C(p,g0) cross section was up to
# 9.7 times the one of the same levels in the other order (models/A_orig.azr), while --use-brune and the
# elastic scattering of both runs agreed exactly.
#
# Checked here, for a two-level (A) and a three-level (B) 1/2- group:
#   - capture and elastic scattering do not depend on the order of the levels, in either parametrization;
#   - the formal and the Brune parametrization give the same cross sections for both orders;
#   - the sign of the entered ANC: a negative ANC changes the capture cross section, and gives the same one
#     in both orders and both parametrizations;
#   - a static moment on the ground state (the M1 line to its own pair) is reported as a moment in both
#     orders (the test compared the level's FORMAL energy with the pair's before);
#   - the gamma widths reported in parameters.out by a formal run are those of the Brune run (with the final
#     state's row returned in the wrong sign their external part had the wrong relative sign);
#   - --no-transform (formal input parameters): the formal parameters of a --no-brune run, written into the
#     level lines, reproduce that run in both orders; a Note says that the boundary condition is at the lowest
#     level when that is not the first listed one; a final state moved off its pair's energy is refused;
#   - with capture pseudo-data: chi^2 is the same for both orders and both parametrizations, the analytic
#     gradient is used and agrees with finite differences for both orders (with every parameter of the
#     final state's group free: by finite differences, announced), and a formal fit comes back to the data.
#
# From dev's tests/formal_ec_final_state (081dddf9, 630dff75), models converted with --convert; README.md says
# what differs from dev's script and why.
#
#   ./tests/formal_ec_final_state/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
. "$HERE/../lib/guard.sh"
run_limited() { local t="$1"; shift; $(guard_command "$t") "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/formal_ec_final_state.XXXXXX")"
[ -n "${KEEP_WORK:-}" ] || trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {  # NAME MODEL FLAGS
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$2" "$dir/t.azr"
  (cd "$dir" && run_limited 600 "$AZURE2_BIN" extrapolate $3 t.azr < /dev/null > run.log 2>&1)
}
xs() { awk 'NF >= 4 { printf "%.12e %.12e\n", $1, $(NF - 1) }' "$1"; }
maxrel() {
  paste "$1" "$2" | awk '
    NF >= 4 { d = $2 - $4; if (d < 0) d = -d; m = $4 < 0 ? -$4 : $4; if (m < 1e-300) m = 1e-300;
              r = d / m; if (r != r || r > worst) worst = r; rows++ }
    END { if (rows == 0) print "nan"; else printf "%.3e\n", worst }'
}
# both must be plain numbers: "nan" (no rows, a failed run) compared as a string and passed before 2026-11-12
le() { awk -v r="$1" -v t="$2" 'BEGIN { n = "^[0-9]+\\.?[0-9]*(e[-+]?[0-9]+)?$"; exit !(r ~ n && t ~ n && r + 0 <= t + 0) }'; }
same() {  # LABEL FILE1 FILE2 TOLERANCE
  local r; r="$(maxrel "$2" "$3")"
  if le "$r" "$4"; then ok "$1 (max rel $r)"; else bad "$1: max rel $r > $4"; fi
}

CAP="output/AZUREOut_aa=1_R=2.extrap"; ELA="output/AZUREOut_aa=1_R=1.extrap"
M="$HERE/models"
# The formal run reaches the Brune one through an eigenvalue problem (TransformIn) and its inverse
# (the final state's on-resonance amplitude): agreement to rounding of those, not bit for bit.
TOL=1e-7
for pair in "A_orig A_swap" "B_orig B_mid"; do
  set -- $pair; first="$1"; second="$2"
  missing=0
  for m in "$first" "$second"; do
    for f in use-brune no-brune; do
      run "${m}_$f" "$M/$m.azr" "--$f" || true
      for o in "$CAP" "$ELA"; do
        if [ ! -s "$WORK/${m}_$f/$o" ]; then
          bad "$m --$f: no $(basename "$o") ($(tail -1 "$WORK/${m}_$f/run.log" | cut -c1-80))"; missing=1
        else
          xs "$WORK/${m}_$f/$o" > "$WORK/${m}_$f.$( [ "$o" = "$CAP" ] && echo cap || echo ela )"
        fi
      done
    done
  done
  [ "$missing" -eq 0 ] || continue
  n="$(wc -l < "$WORK/${first}_use-brune.cap")"
  [ "$n" -gt 100 ] || bad "$first: only $n capture points"
  if grep -q "nan" "$WORK/${second}_no-brune.cap"; then bad "$second --no-brune: NaN in the capture cross section"; fi
  same "$second vs $first, --use-brune, capture" "$WORK/${second}_use-brune.cap" "$WORK/${first}_use-brune.cap" 1e-10
  same "$second vs $first, --use-brune, elastic" "$WORK/${second}_use-brune.ela" "$WORK/${first}_use-brune.ela" 1e-10
  same "$second vs $first, --no-brune, capture" "$WORK/${second}_no-brune.cap" "$WORK/${first}_no-brune.cap" "$TOL"
  same "$second vs $first, --no-brune, elastic" "$WORK/${second}_no-brune.ela" "$WORK/${first}_no-brune.ela" "$TOL"
  same "$first, formal vs Brune, capture" "$WORK/${first}_no-brune.cap" "$WORK/${first}_use-brune.cap" "$TOL"
  same "$second, formal vs Brune, capture" "$WORK/${second}_no-brune.cap" "$WORK/${second}_use-brune.cap" "$TOL"
  same "$second, formal vs Brune, elastic" "$WORK/${second}_no-brune.ela" "$WORK/${second}_use-brune.ela" "$TOL"
  # the corrected final state must not be announced as displaced: nothing was fitted
  if grep -q "external-capture final state in the formal parametrization" "$WORK/${second}_no-brune/run.log"; then
    bad "$second --no-brune: a pole-offset warning in an unfitted run"
  else ok "$second --no-brune: no pole-offset warning"; fi
done

# The sign of the entered ANC.  The external-capture amplitude is linear in the final state's reduced width
# amplitude, so its sign matters (interference with the internal capture); it must follow the ANC of the
# input in every order and parametrization.
mkdir -p "$WORK/models"
for m in A_orig A_swap B_mid; do
  awk '$1 == "0.5" && $2 == "-1" && $12 == "1.62" { sub(/ 1\.62 /, " -1.62 ") } { print }' "$M/$m.azr" > "$WORK/models/${m}_neg.azr"
  for f in use-brune no-brune; do
    run "${m}_neg_$f" "$WORK/models/${m}_neg.azr" "--$f" || true
    xs "$WORK/${m}_neg_$f/$CAP" > "$WORK/${m}_neg_$f.cap" 2>/dev/null
  done
done
r="$(maxrel "$WORK/A_orig_neg_use-brune.cap" "$WORK/A_orig_use-brune.cap")"
if le 1e-3 "$r"; then ok "negative ANC changes the capture cross section (max rel $r)"; else bad "negative ANC: capture unchanged (max rel $r)"; fi
same "negative ANC, A_orig, formal vs Brune, capture" "$WORK/A_orig_neg_no-brune.cap" "$WORK/A_orig_neg_use-brune.cap" "$TOL"
same "negative ANC, A_swap --no-brune vs A_orig Brune, capture" "$WORK/A_swap_neg_no-brune.cap" "$WORK/A_orig_neg_use-brune.cap" "$TOL"
same "negative ANC, B_mid --no-brune vs A_orig Brune, capture" "$WORK/B_mid_neg_no-brune.cap" "$WORK/A_orig_neg_use-brune.cap" "$TOL"

# A static moment on the ground state: its M1 line to its own pair with a non-zero value.  parameters.out
# must print it as a moment ("mu =") with the same numbers in both level orders.
# A printed zero carries no sign here: with the ground state second its energy comes back from the formal
# transform as round-off (7.4e-18 MeV on Linux gcc, where the first order gives exactly 0), and a build that rounds
# the other way prints -0.0000 (the Windows runner failed this comparison, 2026-11-17).  The fields are rejoined by
# single blanks in both files.
gsblock() { awk '{ sub(/\r$/, "") } /^J = / { on = ($1 == "J" && $3 == "0.5-" && $6 + 0 < 1 && $6 + 0 > -1) }
  on && NF { gsub(/\(-0\.0+,/, "(0.000000,"); gsub(/,-0\.0+\)/, ",0.000000)")
             for (i = 1; i <= NF; i++) if ($i ~ /^-0\.0+$/) $i = substr($i, 2)
             if ($1 == "J") print $1, $2, $3, $4, $5, $6; else print }' "$1"; }
for m in A_orig A_swap; do
  awk '$1 == "0.5" && $2 == "-1" && $3 == "0" && $6 == "2" { $12 = "-0.3222"; print " " $0; next } { print }' "$M/$m.azr" > "$WORK/models/${m}_mom.azr"
  run "${m}_mom" "$WORK/models/${m}_mom.azr" "--no-brune" || true
  gsblock "$WORK/${m}_mom/output/parameters.out" > "$WORK/${m}_mom.gs" 2>/dev/null
done
if grep -q "mu = *-0.3222" "$WORK/A_orig_mom.gs"; then ok "static moment, ground state first: reported as mu"; else bad "static moment, ground state first: no 'mu = -0.3222' in parameters.out"; fi
# Byte-identical files without cmp (the Windows MSYS2 runner has none); gsblock has already dropped any \r.
same_file() { [ -f "$2" ] && [ "$(cat "$1")" = "$(cat "$2")" ]; }
if [ -s "$WORK/A_swap_mom.gs" ] && same_file "$WORK/A_orig_mom.gs" "$WORK/A_swap_mom.gs"; then ok "static moment, ground state second: same ground-state lines in parameters.out"
else
  bad "static moment, ground state second: parameters.out differs ($(grep "L = " "$WORK/A_swap_mom.gs" 2>/dev/null | tail -1 | cut -c1-60))"
  # both blocks, so that the log of a runner nobody can log into shows the difference
  for m in A_orig A_swap; do sed "s/^/    $m: /" "$WORK/${m}_mom.gs" 2>/dev/null; done
fi

# The gamma widths REPORTED in parameters.out.  A reported width holds the internal and the external part of
# the capture amplitude; with the final state's transformed row returned in the wrong sign the external part
# entered with the wrong relative sign, and B_mid --no-brune reported 1.228 eV for the entered 480 meV of the
# 1/2+ level (the cross sections were right).  The formal run must report what the Brune run does.  Zero
# widths are left out, so that the ground state's M1 line (a moment or a width, the check above) does not
# count here.
gwidths() { awk '{ sub(/\r$/, "") } $1 == "J" { j = $3 " " $6 } $1 == "R" && $4 == "L" && $10 == "G" && $12 + 0 != 0 { u = $13 == "meV" ? 1e-3 : $13 == "keV" ? 1e3 : $13 == "MeV" ? 1e6 : 1; printf "%s %s %.9e eV\n", j, $6, $12 * u }' "$1"; }
for m in A_swap B_mid; do
  gwidths "$WORK/${m}_no-brune/output/parameters.out" > "$WORK/${m}_no-brune.gw" 2>/dev/null
  gwidths "$WORK/${m}_use-brune/output/parameters.out" > "$WORK/${m}_use-brune.gw" 2>/dev/null
  paste "$WORK/${m}_no-brune.gw" "$WORK/${m}_use-brune.gw" > "$WORK/${m}.gw2"
  w="$(awk '
    bad { next }
    NF != 10 || $1 != $6 || $3 != $8 || $5 != $10 { print "layout"; bad = 1; next }
    { d = $4 - $9; if (d < 0) d = -d; s = $9 < 0 ? -$9 : $9; if (d > 1e-6 * s + 1e-9) { print $1, $2, $3, $4, $5, "for", $9, $10; bad = 1; next } rows++ }
    END { if (!bad && rows < 4) print "only " rows + 0 " gamma widths" }' "$WORK/${m}.gw2")"
  if [ -z "$w" ]; then ok "$m --no-brune: reported gamma widths are those of the Brune run"; else bad "$m --no-brune: reported gamma width differs from the Brune run ($w)"; fi
done

# Formal input parameters (--no-transform).  The level lines then hold formal energies (MeV) and formal
# reduced width amplitudes (MeV^1/2).  Round trip: the formal parameters of the --no-brune runs above
# (output/param.par, 8 digits, hence the 1e-5) written into the level lines must give the same cross
# sections.  The final state is found by its pair's energy; its formal energy is 3 MeV from it when the
# ground state is listed second, and a NOTE says which level was taken.  A WARNING is due only when the
# pole of the group is not at the pair's energy: here with the ground state's formal energy moved by 50 keV.
formal_model() {  # PARAM.PAR MODEL [LEVEL SHIFT_MEV]: the model with the formal parameters in its level lines
  awk -v ps="$1" -v sl="${3:-0}" -v sh="${4:-0}" '
    BEGIN { while ((getline line < ps) > 0) { sub(/\r$/, "", line); split(line, f, " "); p[f[1]] = f[2] } }
    /^<levels>/ { lev = 1; print; next } /^<\/levels>/ { lev = 0; print; next }
    !lev { print; next }
    NF == 0 { open = 0; print; next }
    { if (!open) { n++; c = 0; open = 1 } c++
      e = p["energy_" n] + 0; if (n == sl) e += sh
      $3 = sprintf("%.10g", e); $12 = sprintf("%.10g", p["width_" n "_" c] + 0); print " " $0 }' "$2"
}
# This branch sets the boundary condition of a group at its LOWEST level, whatever the order of the lines, and says
# so in a Note when that is not the first listed one; a final state that is not within 1 keV of its pair's energy
# is refused (dev: the nearest level is taken, with a NOTE, and a WARNING when the pole is off the pair's energy).
W="needs its final state"
NOTE="the boundary condition is set at the lowest level"
for m in A_orig A_swap A_swap_neg; do
  src="$M/$m.azr"; [ -f "$src" ] || src="$WORK/models/$m.azr"
  par="$WORK/${m}_no-brune/output/param.par"
  if [ ! -s "$par" ]; then bad "--no-transform $m: no param.par of the --no-brune run"; continue; fi
  formal_model "$par" "$src" > "$WORK/models/${m}_formal.azr"
  run "${m}_nt" "$WORK/models/${m}_formal.azr" "--no-brune --no-transform" || true
  xs "$WORK/${m}_nt/$CAP" > "$WORK/${m}_nt.cap" 2>/dev/null; xs "$WORK/${m}_nt/$ELA" > "$WORK/${m}_nt.ela" 2>/dev/null
  [ -s "$WORK/${m}_no-brune.cap" ] || xs "$WORK/${m}_no-brune/$CAP" > "$WORK/${m}_no-brune.cap"
  [ -s "$WORK/${m}_no-brune.ela" ] || xs "$WORK/${m}_no-brune/$ELA" > "$WORK/${m}_no-brune.ela"
  same "--no-transform $m: capture of the formal parameters" "$WORK/${m}_nt.cap" "$WORK/${m}_no-brune.cap" 1e-5
  same "--no-transform $m: elastic of the formal parameters" "$WORK/${m}_nt.ela" "$WORK/${m}_no-brune.ela" 1e-5
  if grep -q "$W" "$WORK/${m}_nt/run.log"; then bad "--no-transform $m: warned although the pole is at the pair's energy"; else ok "--no-transform $m: no warning"; fi
done
if grep -q "$NOTE" "$WORK/A_swap_nt/run.log"; then ok "--no-transform, ground state second: Note says the boundary condition is at the lowest level"; else bad "--no-transform, ground state second: no Note on the boundary condition's level"; fi
if grep -q "$NOTE" "$WORK/A_orig_nt/run.log"; then bad "--no-transform, ground state first: Note without cause"; else ok "--no-transform, ground state first: no Note"; fi
if grep -q "$W" "$WORK/A_swap_no-brune/run.log" 2>/dev/null; then bad "transformed run carries the --no-transform warning"; fi
# displaced: ground-state formal energy 50 keV lower.  The lowest level is the pole in both orders, and it is no
# longer at the pair's energy: the run is refused with a message that names the pair, and writes no cross section.
for c in "A_orig 1" "A_swap 2"; do
  set -- $c
  formal_model "$WORK/$1_no-brune/output/param.par" "$M/$1.azr" "$2" -0.05 > "$WORK/models/$1_shift.azr"
  if run "$1_shift" "$WORK/models/$1_shift.azr" "--no-brune --no-transform"; then
    bad "--no-transform $1, formal energy moved 50 keV: the run was not refused"
  elif grep -q "$W" "$WORK/$1_shift/run.log" && [ ! -s "$WORK/$1_shift/$CAP" ]; then
    ok "--no-transform $1, formal energy moved 50 keV: refused, no final state within 1 keV of the pair's energy"
  else bad "--no-transform $1, formal energy moved 50 keV: failed without the message ($(tail -1 "$WORK/$1_shift/run.log" | cut -c1-80))"; fi
done

# An ANC of exactly zero gives the sign rule above no direction.  The final state's formal row must then be
# its input row (DeBoer 2026-11-13), which is also the row the rule picks for an ANC of +1e-6 or of -1e-6
# (the amplitude used is a_f times the entered one), so zero is the continuous point.  The ground state
# carries a static moment (its M1 line), so its row is not zero and its sign can be seen: width_<gs>_2 in
# param.par, and the moment printed in parameters.out.  Before, the solver's sign decided; with the ground
# state second of three levels (B_mid) its row came back reversed and parameters.out printed mu = +0.3222
# for the entered -0.3222 (the cross sections were right).
zero_model() {  # MODEL ANC
  awk -v a="$2" '$1 == "0.5" && $2 == "-1" && $3 == "0" && $6 == "1" && $12 == "1.62" { $12 = a; print " " $0; next }
                 $1 == "0.5" && $2 == "-1" && $3 == "0" && $6 == "2" { $12 = "-0.3222"; print " " $0; next } { print }' "$1"
}
parval() { awk -v k="$2" '{ sub(/\r$/, "") } $1 == k { print $2; exit }' "$1"; }
for c in "A_swap 2" "B_mid 2"; do
  set -- $c
  for a in 0 1e-6 -1e-6; do
    zero_model "$M/$1.azr" "$a" > "$WORK/models/$1_anc$a.azr"
    run "$1_anc$a" "$WORK/models/$1_anc$a.azr" "--no-brune" || true
  done
  z="$(parval "$WORK/$1_anc0/output/param.par" "width_$2_2")"; p="$(parval "$WORK/$1_anc1e-6/output/param.par" "width_$2_2")"
  n="$(parval "$WORK/$1_anc-1e-6/output/param.par" "width_$2_2")"; z1="$(parval "$WORK/$1_anc0/output/param.par" "width_$2_1")"
  if [ -n "$z" ] && [ -n "$p" ] && [ -n "$n" ] && awk -v z="$z" -v p="$p" -v n="$n" -v z1="$z1" '
       function abs(x) { return x < 0 ? -x : x }
       BEGIN { exit !(z + 0 != 0 && abs(z - p) <= 1e-4 * abs(z) && abs(z - n) <= 1e-4 * abs(z) && z1 !~ /^-/ && abs(z1) < 1e-12) }'; then
    ok "zero ANC, $1: ground-state formal row = input row = that of ANC +-1e-6 (moment amplitude $z)"
  else bad "zero ANC, $1: ground-state moment amplitude '$z', with ANC +1e-6 '$p', -1e-6 '$n' (particle '$z1')"; fi
  if grep -q "mu = *-0.3222" "$WORK/$1_anc0/output/parameters.out" 2>/dev/null; then ok "zero ANC, $1: static moment reported with the entered sign"
  else bad "zero ANC, $1: no 'mu = -0.3222' in parameters.out ($(awk 'match($0, /mu = *[-0-9.]*/) && !seen { seen = 1; found = substr($0, RSTART, RLENGTH) } END { print found }' "$WORK/$1_anc0/output/parameters.out" 2>/dev/null))"; fi
done

# Data mode.  Pseudo-data: the Brune capture cross section of A_orig (every fourth point, 5 %), entered in
# the laboratory frame; then the 1/2+ width is moved from 34 to 30 keV so that chi^2 is far from zero.
# chiSquared.out carries six digits, hence the 2e-5.
SRC="$WORK/A_orig_use-brune/$CAP"
if [ ! -s "$SRC" ]; then
  bad "no Brune capture cross section to make the pseudo-data from"
else
  mkdir -p "$WORK/fitdata"
  awk 'NF >= 5 && (++n % 4) == 1 { printf "%.10e 0 %.10e %.10e\n", $1 * (1.00728 + 12) / 12, $(NF - 1), 0.05 * $(NF - 1) }' \
    "$SRC" > "$WORK/fitdata/pseudo.dat"
  fitmodel() {  # NAME MODEL
    local dir="$WORK/$1"
    mkdir -p "$dir/output" "$dir/checks" "$dir/data"
    cp "$WORK/fitdata/pseudo.dat" "$dir/data/"
    awk '
      /^<data>/ { print; print "dataset pseudo"; print "    reaction       12C(p,g)13N"; print "    quantity       sigma"
                  print "    unit           b"; print "    file           data/pseudo.dat"
                  print "    columns        E_lab[MeV] _ value error"; print "    energy_window  0 5 MeV"
                  print "    norm           1 fixed"; next }
      !moved && $3 == "2.3682" && sub(/ 34000 /, " 30000 ") { moved = 1 }
      { print }' "$2" > "$dir/t.azr"
  }
  total() { awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$1" 2>/dev/null; }
  calc() {  # NAME FLAGS -> total chi^2 of a calculation
    rm -f "$WORK/$1/output/chiSquared.out"
    (cd "$WORK/$1" && run_limited 600 "$AZURE2_BIN" calculate $2 t.azr < /dev/null > calc.log 2>&1)
    total "$WORK/$1/output/chiSquared.out"
  }
  close() { awk -v a="$1" -v b="$2" 'BEGIN { if (a == "" || b == "" || b + 0 == 0) exit 1; d = (a - b) / b; if (d < 0) d = -d; exit !(d <= 2e-5) }'; }
  fitmodel fit_orig "$M/A_orig.azr"; fitmodel fit_swap "$M/A_swap.azr"
  ref="$(calc fit_orig --use-brune)"
  if [ -n "$ref" ] && awk -v c="$ref" 'BEGIN { exit !(c > 100) }'; then ok "pseudo-data: chi^2 $ref at the moved start (Brune, ground state first)"
  else bad "pseudo-data: chi^2 '$ref' at the moved start, expected > 100"; fi
  for c in "fit_orig --no-brune" "fit_swap --use-brune" "fit_swap --no-brune"; do
    set -- $c; got="$(calc "$1" "$2")"
    if close "$got" "$ref"; then ok "chi^2 $1 $2 = $got"; else bad "chi^2 $1 $2 = '$got', expected $ref"; fi
  done

  gradcheck() {  # NAME -> the [grad-check] summary line of a formal fit start
    (cd "$WORK/$1" && AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
       run_limited 600 "$AZURE2_BIN" fit --no-brune --use-gradient t.azr < /dev/null > grad.log 2>&1)
    awk '{ sub(/\r$/, "") } /^\[grad-check\] chi2 / { print; exit }' "$WORK/$1/grad.log"
  }
  # fit_free: every energy and proton width of the final state's 1/2- group free as well, so that the
  # derivative through the on-resonance amplitude (CNuc::FormalECFinalPullback) carries weight: the
  # background pole's parameters move the ground state's amplitude only through the eigenvector.
  mkdir -p "$WORK/fit_free/output" "$WORK/fit_free/checks" "$WORK/fit_free/data"
  cp "$WORK/fitdata/pseudo.dat" "$WORK/fit_free/data/"
  awk '/^<levels>/ { lev = 1 } /^<\/levels>/ { lev = 0 }
       lev && $1 == "0.5" && $2 == "-1" { $4 = 0; if ($6 == 1) $11 = 0; print " " $0; next } { print }' \
    "$WORK/fit_swap/t.azr" > "$WORK/fit_free/t.azr"
  for c in "fit_swap:ground state second" "fit_orig:ground state first"; do
    line="$(gradcheck "${c%%:*}")"
    case "$line" in
      *"analytic block used"*"energy 0/"*"gamma 0/"*) ok "${c#*:}: analytic gradient used, no disagreement with finite differences" ;;
      *) bad "${c#*:}: analytic gradient ('$(echo "$line" | cut -c1-140)')" ;;
    esac
  done
  # With the final state's whole group free this branch does not differentiate the final state's amplitude
  # analytically (AZUREGrad.cpp, HasFormalECFinalBlocks): the gradient is by finite differences and the run says so
  # (dev: analytic, CNuc::FormalECFinalPullback; open in README.md).  Either is accepted; a silent third is not.
  line="$(gradcheck fit_free)"
  case "$line" in
    *"analytic block used"*"energy 0/"*"gamma 0/"*) ok "its whole group free: analytic gradient used, no disagreement with finite differences" ;;
    *"analytic block BAILED (all finite differences)"*) ok "its whole group free: gradient by finite differences, announced" ;;
    *) bad "its whole group free: gradient ('$(echo "$line" | cut -c1-140)')" ;;
  esac
  # ... and a formal fit of that model comes back to the data as well
  rm -f "$WORK/fit_free/output/chiSquared.out"
  (cd "$WORK/fit_free" && run_limited 900 "$AZURE2_BIN" fit --no-brune --use-gradient t.azr < /dev/null > fit.log 2>&1)
  got="$(total "$WORK/fit_free/output/chiSquared.out")"
  if [ -n "$got" ] && awk -v c="$got" 'BEGIN { exit !(c ~ /^ *[-+]?\.?[0-9]/ && c + 0 < 1) }'; then ok "formal fit, its whole group free: chi^2 $ref -> $got"
  else bad "formal fit, its whole group free: chi^2 '$got' ($(tail -1 "$WORK/fit_free/fit.log" | cut -c1-80))"; fi

  # a formal fit from the corrected order: back to the data (the true parameters give chi^2 = 0)
  rm -f "$WORK/fit_swap/output/chiSquared.out"
  (cd "$WORK/fit_swap" && run_limited 900 "$AZURE2_BIN" fit --no-brune --use-gradient t.azr < /dev/null > fit.log 2>&1)
  got="$(total "$WORK/fit_swap/output/chiSquared.out")"
  if [ -n "$got" ] && awk -v c="$got" 'BEGIN { exit !(c ~ /^ *[-+]?\.?[0-9]/ && c + 0 < 1) }'; then ok "formal fit, ground state second: chi^2 $ref -> $got"
  else bad "formal fit, ground state second: chi^2 '$got' ($(tail -1 "$WORK/fit_swap/fit.log" | cut -c1-80))"; fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: formal-parametrization external-capture final states"; else echo "FAIL: $failures check(s)"; exit 1; fi
