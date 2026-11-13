# Runtime options that change the physics

Each option is a CLI flag, a keyword of `Session()` and a `<config>` key of the same name as the keyword (`use_brune=0`;
a flag or keyword wins over the key, the key over the default; `AZURE2 --describe` shows the source). Full list:
`docs/source/reference/command_line.rst`. Take the flags of an established project from its job scripts; a calculation
under different flags is a different model.

## Route: A-matrix or R-matrix (`<config>` first line)

- The R-matrix route implements neither ignore-externals (on by default) nor Park. With either on, every front end
  (CLI, pyazr) runs the A-matrix and warns "the R-matrix route ... does not implement ...".
- To run the R-matrix route: `--no-ignore-externals` (pyazr `ignore_externals=False`) and no `--use-park`. It equals the
  A-matrix route to ~10⁻¹³ in the Brune or standard basis.
- Builds before 2026-11-01 (expires: 2027-02-01): the CLI never ran the R-matrix route; pyazr and the GUI ran it with
  Park amplitudes read as Brune's (up to 9× off on 13N) and with ignore-externals on (up to 4.5 % off). Exit status was
  always 0.

## External-capture operator

- `--no-long-wavelength` changes external capture (EC) only; resonant photon channels are untouched.
  `--ignore-externals` does not switch EC off.
- The exact operator is the extended-Siegert charge term (Friar and Fallieros, Phys. Rev. C 29 (1984) 1645) plus the
  convection current of the relative motion; `--no-ec-convection-current` (pyazr `ec_convection_current=False`) keeps
  the charge term alone. Derivation: `docs/source/theory/external_capture_operator.rst`; independent check:
  `tests/ec_m1_exact`.
- Size: about 0.1 % on ordinary E1 and E2 EC. Where the long-wavelength E1 effective charge Z₁m₂/M − Z₂m₁/M nearly
  vanishes (¹²C+α, ¹⁰B+α), the exact E1 EC amplitude is 3.5–4.7 × the long-wavelength one, plus 10 % from the
  convection current (¹²C(α,γ₀): 0.1–0.4 % in S at fixed parameters).
- Keep the photon pair's final-state energy (`e2`) equal to the EC level's energy: the ANC conversion uses the former
  and the radial integral the latter; a 0.03 keV mismatch moves EC by 0.16 %.
- Delete `output/intEC.*` when toggling the operator flags.
- Older-build errors of the M1 EC term (6.25 % low; wrong spin-term sign for half-integer partners; spurious orbital
  term on channel-spin-flip pathways): `older-builds.md`.

## Hybrid nuclear potential (`potential` record; older projects: `<potential>` block)

- Format 3: one `potential` record per particle pair in `<potentials>` (`pair`, required `coulomb_radius`, `central
  woods_saxon|gaussian ...`, optional `matching_radius` 15 fm and `step` 0.0275 fm); no record, no potential.
  `AZURE2 --convert OLD NEW` writes the records from an old block (`coulomb_radius 3.5 fm`). Fields and
  refusals: `docs/source/reference/format3.rst`.
- The potential changes the radial solution between a_c and the matching radius, hence P, S, the phases and the EC
  entrance wave, and below threshold the decaying solution (S, ANC conversion, bound-state tail of EC). Behaviour:
  `docs/source/user_guide/pyazr.rst`.
- Optional `spin_orbit woods_saxon depth .. radius .. diffuseness .. on light|heavy operator 2l.s` for a pair with
  one spinning particle (FRESCO's convention: a strength quoted for l.s without the 2 is twice V_so); the channel
  functions then depend on J (`Session.channel_functions(..., J=)` needs it).  Pairs with two spins are refused
  (`E_SO_TWO_SPINS`).  Theory: `docs/source/theory/hybrid_potential.rst`.
- A pyazr session owns its potential: several sessions with different potentials in one process stay independent
  (`tests/pyazr/potential_sessions_test.py`). Builds before 2026-10-27 (expires: 2027-02-01) shared one potential per
  process (the last session opened set it for all) and had an O(dr) derivative of the Numerov solution at a_c (shift
  function off by 4 × 10⁻³ to 4 × 10⁻² for l = 0–3).

## Coulomb functions

- The universal Coulomb table is on by default (`--no-coulomb-table`, pyazr `use_coulomb_table=False` turns it off):
  P, S and the hard-sphere phase of charged channels, open and closed. Build it once with `AZURE2 --build-coulomb-table`
  into `~/.azure2/coulomb_table.bin` or `$AZURE2_COULOMB_TABLE`; without the file
  the log says "Coulomb table: none at ..." and the exact routines run. Not tabulated: neutral pairs, pairs with a
  nuclear potential, the EC integrals. Design and validation: `docs/source/developer/coulomb_table.rst`.
- Closed neutral channels use the closed form S₀ = −ρ, S_l = −ρ²/(l − S_{l−1}) − l, ρ = κa (agrees with the numerical
  route to 3 × 10⁻⁹; `AZURE_NO_NEUTRAL_CLOSED_FORM=1` restores the old route for comparison).
- Charged-particle fits with bound or subthreshold levels differ from builds before 2026-10-19 (expires: 2027-02-01)
  by up to ~10⁻³ in χ² (exact-route dS/dE corrected; `AZURE_BOUND_SHIFT_NUMERICAL=1` and `AZURE_SHIFT_DE_1EV=1` restore
  the old routes).
- **`--gsl-coul` trades accuracy for speed.** The default routine (Michel 2007) agrees with mpmath to 10⁻⁹ and is
  slower; GSL's is faster but wrong by 3–6 % in P and S for 0.1 ≤ η < 10, 1 ≤ ρ < 10. Use it only for timing or
  exploratory runs where that error is acceptable, never for results; the Coulomb table gives the speed without the
  error.

## External-capture energy table

On by default (`--no-ec-table`, pyazr `use_ec_table=False`): each capture term's radial integral is tabulated in
ln E_out, 50 nodes per decade, nodes converged to 10⁻¹⁰. `AZURE2_EC_TABLE_CHECK=1` compares every tabulated integral
with a direct one; `AZURE2_EC_TABLE_NODES=N` changes the density; the log prints "External capture table: ...".
Design: `docs/source/developer/ec_table.rst`. Shifted capture segments are evaluated at the shifted energies.
Builds before 2026-10-21 (expires: 2027-02-01) interpolated EC amplitudes linearly between a segment's own points and
held them flat beyond its ends (0.4–15 % errors, a factor 57 at the low end of a shifted data set): fitted energy
shifts of capture segments from those builds are biased and need a refit.

## Relativistic kinematics

`relativistic=1` in `<config>` changes η, ρ and the flux factor of every EC amplitude (0.1–1 %). It is part of the
`intEC` signature and of the EC table key. A relativistic run from builds before 2026-10-27 (expires: 2027-02-01)
could reuse a non-relativistic `intEC.dat`: rerun once with fresh caches.

## Cache signatures

`intEC.dat` written by an older build is recomputed once, with a WARNING, when its signature differs. Files from builds
before 2026-10-30 (expires: 2027-02-01) held 6 significant digits, so a run that reused one differed from a run that
computed it by ~2 × 10⁻⁶ in χ² (up to ~10⁻⁵ on single capture points). When comparing older runs to better than
10⁻⁵, check whether one computed and the other reused `intEC.dat`.

A pyazr session that reloads an existing `intEC.dat` has returned a wrong χ² on at least one large model, while the CLI
reloading the same file was exact. Until that is understood, pyazr scripts delete `output/intEC.*` first or use a fresh
output directory.

## Parametrization and formalism

Brune (default), Park (`--use-park`), standard Lane–Thomas (`--no-brune`), Reich–Moore (`--use-rmc`, (n,γ) only),
`--no-transform`: what each basis means and how to convert between them is in azure2-exchange.
