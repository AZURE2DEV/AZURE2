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
#     level lines, reproduce that run in both orders; the final state is the level nearest its pair's energy
#     (a NOTE says which when none is within 1 keV), and a WARNING appears only when the pole of the group
#     is not at the pair's energy;
#   - with capture pseudo-data: chi^2 is the same for both orders and both parametrizations, the analytic
#     gradient is used and agrees with finite differences for both orders (also with every parameter of the
#     final state's group free), and a formal fit from the corrected order comes back to the data.
#
#   ./tests/formal_ec_final_state/check.sh path/to/AZURE2

set -uo pipefail
export LC_ALL=C

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AZURE2_BIN="${1:?usage: check.sh path/to/AZURE2}"
AZURE2_BIN="$(cd "$(dirname "$AZURE2_BIN")" && pwd)/$(basename "$AZURE2_BIN")"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
run_limited() { if [ -n "$TIMEOUT_CMD" ]; then "$TIMEOUT_CMD" "$@"; else shift; "$@"; fi; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/formal_ec_final_state.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

failures=0
ok() { echo "  ok    $1"; }
bad() { echo "  FAIL  $1"; failures=$((failures + 1)); }

run() {  # NAME MODEL FLAGS
  local dir="$WORK/$1"
  mkdir -p "$dir/output" "$dir/checks"
  cp "$2" "$dir/t.azr"
  (cd "$dir" && printf '3\nn\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" $3 --no-gui --no-readline t.azr > run.log 2>&1)
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
gsblock() { awk '{ sub(/\r$/, "") } /^J = / { on = ($1 == "J" && $3 == "0.5-" && $6 + 0 < 1 && $6 + 0 > -1) } on && NF { if ($1 == "J") print $1, $2, $3, $4, $5, $6; else print }' "$1"; }
for m in A_orig A_swap; do
  awk '$1 == "0.5" && $2 == "-1" && $3 == "0" && $6 == "2" { $12 = "-0.3222"; print " " $0; next } { print }' "$M/$m.azr" > "$WORK/models/${m}_mom.azr"
  run "${m}_mom" "$WORK/models/${m}_mom.azr" "--no-brune" || true
  gsblock "$WORK/${m}_mom/output/parameters.out" > "$WORK/${m}_mom.gs" 2>/dev/null
done
if grep -q "mu = *-0.3222" "$WORK/A_orig_mom.gs"; then ok "static moment, ground state first: reported as mu"; else bad "static moment, ground state first: no 'mu = -0.3222' in parameters.out"; fi
if [ -s "$WORK/A_swap_mom.gs" ] && cmp -s "$WORK/A_orig_mom.gs" "$WORK/A_swap_mom.gs"; then ok "static moment, ground state second: same ground-state lines in parameters.out"
else bad "static moment, ground state second: parameters.out differs ($(grep "L = " "$WORK/A_swap_mom.gs" 2>/dev/null | tail -1 | cut -c1-60))"; fi

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
  w="$(paste "$WORK/${m}_no-brune.gw" "$WORK/${m}_use-brune.gw" | awk '
    NF != 10 || $1 != $6 || $3 != $8 || $5 != $10 { print "layout"; bad = 1; exit }
    { d = $4 - $9; if (d < 0) d = -d; s = $9 < 0 ? -$9 : $9; if (d > 1e-6 * s + 1e-9) { print $1, $2, $3, $4, $5, "for", $9, $10; bad = 1; exit } rows++ }
    END { if (!bad && rows < 4) print "only " rows + 0 " gamma widths" }')"
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
W="external capture in the formal parametrization without the input parameter transformation"
NOTE="is taken as the final state"
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
if grep -q "$NOTE" "$WORK/A_swap_nt/run.log"; then ok "--no-transform, ground state second: NOTE names the level taken as the final state"; else bad "--no-transform, ground state second: no NOTE on the final-state level"; fi
if grep -q "$NOTE" "$WORK/A_orig_nt/run.log"; then bad "--no-transform, ground state first: NOTE without cause"; else ok "--no-transform, ground state first: no NOTE"; fi
if grep -q "$W" "$WORK/A_swap_no-brune/run.log" 2>/dev/null; then bad "transformed run carries the --no-transform warning"; fi
# displaced: ground-state formal energy 50 keV lower.  With the ground state first that IS the pole (the
# boundary conditions are taken there), so the distance reported must be -50 keV.
for c in "A_orig 1" "A_swap 2"; do
  set -- $c
  formal_model "$WORK/$1_no-brune/output/param.par" "$M/$1.azr" "$2" -0.05 > "$WORK/models/$1_shift.azr"
  run "$1_shift" "$WORK/models/$1_shift.azr" "--no-brune --no-transform" || true
  if grep -q "$W" "$WORK/$1_shift/run.log"; then ok "--no-transform $1, formal energy moved 50 keV: warned"; else bad "--no-transform $1, formal energy moved 50 keV: no warning"; fi
done
d="$(awk '{ sub(/\r$/, "") } /energy is .* keV away/ { for (i = 1; i < NF; i++) if ($i == "is" && $(i + 2) == "keV") { print $(i + 1); exit } }' "$WORK/A_orig_shift/run.log")"
if [ -n "$d" ] && awk -v d="$d" 'BEGIN { exit !(d < -49 && d > -51) }'; then ok "--no-transform A_orig: pole reported $d keV from the pair's energy (true -50)"
else bad "--no-transform A_orig: pole distance '$d' keV, expected -50"; fi

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
    awk -v seg="1 1 2 0 5 0 180 0 1 0 0 0 0 0 data/pseudo.dat 0 0" '
      /^<segmentsData>/ { print; print seg; next }
      !moved && $3 == "2.3682" && sub(/ 34000 /, " 30000 ") { moved = 1 }
      { print }' "$2" > "$dir/t.azr"
  }
  total() { awk '{ sub(/\r$/, "") } /^Total-Chi-Squared:/ { print $2 }' "$1" 2>/dev/null; }
  calc() {  # NAME FLAGS -> total chi^2 of a calculation
    rm -f "$WORK/$1/output/chiSquared.out"
    (cd "$WORK/$1" && printf '1\n\n\n7\n' | run_limited 600 "$AZURE2_BIN" $2 --no-gui --no-readline t.azr > calc.log 2>&1)
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
    (cd "$WORK/$1" && printf '2\n\n\n7\n' | AZURE_GRAD_CHECK=exit AZURE_GRAD_CHECK_N=100000 AZURE_GRAD_CHECK_ROWS=100000 \
       run_limited 600 "$AZURE2_BIN" --no-brune --use-gradient --no-gui --no-readline t.azr > grad.log 2>&1)
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
  for c in "fit_swap:ground state second" "fit_orig:ground state first" "fit_free:ground state second, its whole group free"; do
    line="$(gradcheck "${c%%:*}")"
    case "$line" in
      *"analytic block used"*"energy 0/"*"gamma 0/"*) ok "${c#*:}: analytic gradient used, no disagreement with finite differences" ;;
      *) bad "${c#*:}: analytic gradient ('$(echo "$line" | cut -c1-140)')" ;;
    esac
  done
  # ... and the background pole's width does have a gradient there (it is zero with the ground state first)
  g="$(awk '{ sub(/\r$/, "") } /^\[grad-check\] idx .* gamma J1 L1 ch1 / { for (i = 1; i < NF; i++) if ($i == "analytic") { print $(i + 1); exit } }' "$WORK/fit_free/grad.log")"
  if [ -n "$g" ] && awk -v g="$g" 'BEGIN { exit !(g > 1 || g < -1) }'; then ok "its whole group free: d chi^2 / d gamma(background pole) = $g"
  else bad "its whole group free: no gradient for the background pole's width ('$g')"; fi

  # a formal fit from the corrected order: back to the data (the true parameters give chi^2 = 0)
  rm -f "$WORK/fit_swap/output/chiSquared.out"
  (cd "$WORK/fit_swap" && printf '2\n\n\n7\n' | run_limited 900 "$AZURE2_BIN" --no-brune --use-gradient --no-gui --no-readline t.azr > fit.log 2>&1)
  got="$(total "$WORK/fit_swap/output/chiSquared.out")"
  if [ -n "$got" ] && awk -v c="$got" 'BEGIN { exit !(c == c && c < 1) }'; then ok "formal fit, ground state second: chi^2 $ref -> $got"
  else bad "formal fit, ground state second: chi^2 '$got' ($(tail -1 "$WORK/fit_swap/fit.log" | cut -c1-80))"; fi
fi

if [ "$failures" -eq 0 ]; then echo "PASS: formal-parametrization external-capture final states"; else echo "FAIL: $failures check(s)"; exit 1; fi
