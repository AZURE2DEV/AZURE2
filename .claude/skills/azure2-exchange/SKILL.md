---
name: azure2-exchange
description: Move AZURE2 R-matrix parameters between conventions and codes - the level-parameter basis (standard/formal, Brune, Park; --no-brune, --use-park, the #parametrization tag), export to and import from GNDS/ENDF with LLNL Ferdinand and FUDGE (also FRESCO, EDA, AMUR, RAC, LaTeX), comparison with ENDF/B evaluations, parameters for rmfit and SAMMY. Use when the task names Ferdinand, FUDGE, GNDS, ENDF, a basis conversion or another code's parameters. Fitting in a basis: azure2-fit; setting up a project: azure2-setup.
---

# Conventions and exchange of R-matrix parameters

AZURE2 describes a level by its energy E_λ (MeV, excitation energy E_x of the compound nucleus in `<levels>`) and its
partial widths Γ_λc (eV) or, for closed channels, its ANC (fm^-½); internally the minimizer works on reduced-width
amplitudes γ_λc (MeV^½) in one of three bases. Exchanging parameters with another code, or with an evaluation in
GNDS or ENDF, means stating the basis, the boundary condition, the masses and the energy frame, and checking the
round trip by numbers. The behaviour of the bases is in `docs/source/reference/command_line.rst` ("Brune and Park
parametrizations"); the Brune transformation and its failure modes in `docs/source/reference/troubleshooting.rst`
and `docs/source/user_guide/levels_channels.rst`.

## When to use, and when a sibling skill is the owner

- Use this skill to choose or change the basis, to export to or import from GNDS/ENDF and other codes, and to
  compare an AZURE2 fit with an evaluation.
- Fitting strategy within a basis (fixing, bounds, priors, staged fits): `azure2-fit`.
- Dimensionless widths and Wigner limits as a trust check: `azure2-check`.
- Building a new `.azr` (masses, separation energies, channels): `azure2-setup`.
- rmfit campaigns (`"use_park": true`): `rmfit-campaign`.

## 1. The three bases

| Basis | CLI / pyazr | Parameters | Boundary condition |
|---|---|---|---|
| standard (formal, Lane and Thomas) | `--no-brune`, `use_brune=False` | formal E_λ, γ_λc; transformed from the observed inputs | constant B_c chosen at the level energies |
| Brune (default) | `--use-brune` | observed E_λ, γ_λc^Brune | none needed (Brune, Phys. Rev. C 66 (2002) 044611) |
| Park | `--use-park`, `use_park=True` | observed E_λ, γ_λc^Park = γ_λc^Brune √J_λ | none needed (Park, Phys. Rev. C 104 (2021) 064612) |

with J_λ = 1 − Σ_c (γ_λc^Park)² dS_c/dE |_{E_λ}. Brune and Park give the same collision matrix and read the same
`.azr`; what changes is the space the minimizer works in.

**Choose Park** when the purpose of the fit is to constrain or report widths, for MCMC, and when exchanging with
codes that quote observed amplitudes: a Park amplitude is the observed width alone, Γ_c = 2P_cγ_c², so a fix, bound,
prior or Wigner-limit bound acts on one parameter (on the observed θ²), and the parameters of a level are less
correlated (one p + p test: parabolic errors 0.4 % against 1.3 % in Brune).

**Rules for Park fits.**
- **J_λ must stay positive.** A Park set with J_λ ≤ 0 has no R-matrix meaning (one channel: Γ < 2P/(dS/dE)). The
  objective adds Σ(J_λ/10⁻³)² over J_λ < 0 (`Total-Park-Chi-Squared` in `chiSquared.out`), MCMC rejects such points,
  and a run ending there prints `**WARNING: Park norm J = ... is not positive`. pyazr: `m.park_norms(x)` (one J
  per level, `physical_levels()` order), `m.penalties(x)["park"]`. Levels that saturate in Brune (amplitudes of
  10²–10³ MeV^½) appear in Park as J_λ → 0.
- **Do not expect a different minimum.** The p + p test reaches χ² 99.6121 (Park) and 99.6123 (Brune), widths equal
  to 2 × 10⁻⁵. Near J_λ → 0 the Park form loses precision by cancellation (8 × 10⁻⁷ at the worst point of one
  ¹⁷O model).
- **Gradient checks resolve about 10⁻⁵ in Park, 10⁻⁷ in Brune**: dS/dE in the level matrix is differentiated
  numerically. This is expected, not a defect.
- Park (`--use-park`, pyazr `use_park=True`) implies Brune and excludes RMC. Physical widths of a Park
  vector: `pyazr.transform.transform_out(..., park=True)`. Park-mode MCMC and covariance bands share Brune's code
  paths and are not separately tested.

**Rules for the standard basis.**
- **A Brune fit with very broad background poles may have no standard equivalent**: the observed → formal
  transformation N a = E M a needs M = 1 − Σ γγ′ dS/dE positive definite. The run stops with "No standard R-matrix
  ... corresponds to the levels of J^π = ..." (pyazr `AZURE2Error`); the preceding "Denominator less than zero while
  transforming" warnings name the levels (`tests/7Li_p_ay` is an example). Stay in Brune for such files.
  builds before 2026-11-01 (expires: 2027-02-01) crashed instead ("z is not finite in log_Gamma").
- builds before 2026-10-27 (expires: 2027-02-01): with `--no-brune` the formal eigenvalues went to a J^π group's
  levels in file order, so a group listed out of energy order had energies, widths and fix flags permuted (cross
  sections unaffected). Suspect any old non-Brune fit whose `param.par` `energy_N` values are out of step with the
  file's levels.

## 2. Parameter files carry their basis

`param.par`, `param.sav`, `param.fit` and `param.mcmc` start with `#parametrization <0|1|2>` (standard, Brune,
Park), a comment line so positional readers are unaffected. A file in the other alternative basis is converted on
read (log: `Converted the parameter file from Brune to Park (observed) amplitudes.`). Standard ↔ Brune/Park cannot
be converted level by level: the run stops. An untagged file is read as the running basis's own (Brune for a Park
run).

- **Read `.sav` files by parameter name** (first column), never by row position.
- builds before 2026-11-04 (expires: 2027-02-01) wrote `param.fit` and `param.mcmc` without the tag. Before resuming
  a Park fit from such a file, prepend
  `    #parametrization       2.0000000e+00       0.0000000e+00`; otherwise every amplitude, fixed ones included, is
  multiplied by √J_λ again. Brune and standard runs are unaffected.
- builds before 2026-10-08 (expires: 2027-02-01): a few builds wrote the tag as a plain `parametrization` row, which
  shifts any positional reader by one line; strip it or prefix it with `#`.

## 3. GNDS, ENDF and other codes with Ferdinand and FUDGE

LLNL's Ferdinand (github.com/llnl/ferdinand) translates between AZURE2 `.azr`, GNDS (XML), ENDF, FRESCO/SFRESCO,
EDA, AMUR, RAC and LaTeX; FUDGE (github.com/LLNL/fudge) reads GNDS. Set `$FERDINAND_HOME` to the Ferdinand checkout
(holds `ferdinand.py`) and put the FUDGE checkout on `PYTHONPATH`. Only the AZURE2 ↔ GNDS path and the LaTeX output
have been exercised; options: `python3 $FERDINAND_HOME/ferdinand.py -h`.

```bash
F=$FERDINAND_HOME/ferdinand.py
python3 $F fit.azr gnds -o fit.gnds          # AZURE2 -> GNDS (add -g to drop photon channels, see below)
python3 $F fit.gnds azure -o back.azr        # GNDS -> AZURE2
python3 $F fit.azr tex -o fit.tex            # RevTeX document with particle and parameter tables
```

- The input format comes from the suffix (`.azr` → azure) or `-i <format>`. Without `-o` the output is named
  `<input>+.<format>`, next to the input. Work on a copy of the `.azr`, never in a live fit directory.
- Ferdinand reads the `<levels>` block only. To export fitted values, write the fit back into `<levels>` first
  (`m.save_fit`, `azure2-fit`).
- When writing AZURE2, the log prints "Setting CN zero energy from eSepE = ...": **check that eSepE equals the
  separation energy in the original `.azr`.**

**Known problems.**
1. **Mixed nuclear and atomic masses shift every level energy (GNDS → AZURE2).** Ferdinand recomputes the
   separation energy from the pair masses. A file with a nuclear proton mass (1.00728 u) and an atomic ¹²C mass got
   S_p = 1.43580 MeV instead of 1.94351 MeV, and every level came back 0.5077 MeV low; with both masses atomic
   (p 1.007825 u) the round trip agreed to about 50 eV. Check the projectile mass (level-line token 17) before
   converting; AZURE2 expects atomic masses (`azure2-setup`).
2. **AZURE2-only settings are lost.** After `.azr → GNDS → .azr`: every `levelFix` (token 3) and `channelFix`
   (token 10) is reset to 0, so every parameter becomes free; on capture-pair lines `pType` (token 26) goes 10 → 0,
   `ecMultMask` (token 30) 5 → 0, `chRad` (token 27) 0 → 1.0 and `e2`/`sepE` (tokens 16/22) 0 → ±Q. The γ channel is
   no longer a capture channel and external capture is off. Use GNDS for exchange and archiving; when importing,
   restore fix flags and capture columns from a known-good `.azr`, recalculate, and compare χ² before trusting it.
   Token order: `docs/source/reference/data_formats.rst`; check a round trip token by token with
   `references/ferdinand-checks.md`.
3. **Radiative capture is not handled.** Ferdinand treats the photon pair as a particle channel; with γ channels a
   ¹²C + p Brune fit had its elastic and (p,γ) peaks at E_c.m. ≈ 0.54 MeV instead of 0.424 MeV. Convert particle
   channels with `-g` (log: "Channel photon excluded"); γ widths are eV-scale, so elastic and particle-reaction
   results do not change. Do not use Ferdinand for capture cross sections.
4. **Angular distributions do not reconstruct**: `reconstructLegendre.py` failed with a TensorFlow `MatrixInverse ...
   Input is not invertible`, cause unknown. Compare resonance positions and widths only.

**Reading GNDS and reconstructing cross sections**: `references/ferdinand-checks.md` (FUDGE snippet, pole energies
in the lab frame relative to the entrance threshold, the placeholder compound mass, TensorFlow reconstruction).

## 4. Comparing with an ENDF/B evaluation

Fetching, parsing with the `endf` package, the four LAW=5 representations of charged-particle elastic scattering
and the Rutherford term: `references/endf-comparison.md`. Two rules first:

- **Check the tabulated energy grid before building the comparison.** Light-target proton sublibraries often start
  at E_lab = 1 MeV; the p + ¹²C ENDF/B-VIII.1 file has nothing below it, so a resolved-resonance fit at 0.3–0.6 MeV has
  no ENDF curve to compare with.
- **σ_NI in MF=3/MT=2 of a charged-particle file can be negative**: it is the nuclear-plus-interference integral,
  a difference of cross sections, not corrupted data.

## 5. rmfit and SAMMY

- rmfit uses AZURE2's bases (`"use_park": true` in the campaign settings) and writes the `#parametrization` tag
  when it bakes a fit; procedure in `rmfit-campaign`.
- SAMMY conventions that AZURE2 follows: the free-gas Doppler model, transmission and self-shielded yield, the UDR
  resolution format (`azure2-targets`), and the prior-covariance (PUP) sequential update (`azure2-uncertainty`).
  Parameter files are not exchanged directly; go through GNDS/ENDF or rewrite the levels by hand, then compare χ².

## Pointers

- `docs/source/reference/command_line.rst`: `--no-brune`, `--use-brune`, `--use-park`; "Brune and Park
  parametrizations".
- `docs/source/reference/output_files.rst`: the `#parametrization` line of parameter files.
- `docs/source/reference/troubleshooting.rst`: "No Standard R-Matrix Corresponds to the Levels".
- `docs/source/user_guide/particle_pairs.rst`: EDA and ENDF `KRL=1` conventions.
- `references/ferdinand-checks.md`, `references/endf-comparison.md`.
- Tests: `tests/park_formalism`, `tests/transform_level_order`, `tests/7Li_p_ay`.
