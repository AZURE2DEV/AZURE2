# Special observables

How to compute and fit the observables that are not a plain σ or dσ/dΩ. Segment syntax:
`docs/source/user_guide/segments.rst` and `docs/source/reference/data_formats.rst`.

## Analyzing power A_y (observable 7)

The vector analyzing power of a spin-½ projectile on a target of any spin, in particle and capture exits. Declare
the segment with `observable="analyzing-power"` (code 7 in the raw `.azr`); A_y is reported in place of σ, so χ²,
fitting, plotting and output files need nothing special.

```python
ay = m.calculate_analyzing_power_rwa(m.params_rwa)
```

- **Two routes.** Particle exits use the channel-spin amplitude matrix (Seyler), with the entrance index decomposed
  into projectile and target projections before σ_y acts on the projectile. Capture exits use the Legendre
  coefficients of Seyler and Weller, Phys. Rev. C 20 (1979) 453, Eqs. (20)–(21):
  A_y = Σ_k b_k P_k¹(cos θ) / Σ_k a_k P_k(cos θ), external capture included. `AZURE2_CAPPOL_DEBUG=1` prints, per
  point, the check that Σ a_k P_k reproduces the differential capture cross section. Formalism:
  `docs/source/theory/polarization_theory.rst`, `docs/source/theory/polarization_implementation.rst`.
- **Angles are θ_c.m.** in the data file, unlike an ordinary differential segment: columns
  `E_lab θ_c.m. A_y ΔA_y`.
- **A_y lies in [−1, 1] and is often negative.** Take the experimenters' absolute point-to-point uncertainties
  (statistics, background subtraction). Never convert them to relative errors: a point near a zero crossing does not
  have a small uncertainty, and treating it so gives it enormous weight.
- **Free a normalization.** The beam-polarization calibration P scales A_y multiplicatively
  (A_y ∝ (N⁺ − N⁻)/(P (N⁺ + N⁻))), which is exactly a normalization:
  `Project.set(name, norm="1 fit prior 1 <polarization uncertainty>%")`.
- **No target integration** on A_y segments compared with thin-target data. The target average
  ⟨A_y⟩ = ∫A_y σ dE / ∫σ dE is dominated by the low-energy Rutherford σ where A_y ≈ 0 (≈ 10⁻⁶ on a thick target);
  that is physics. A target-integrated A_y point also switches the whole fit to numerical derivatives, and for
  capture the attenuation coefficients have no counterpart for the P_k¹ numerator.
- **Projectile spin ½ only.** Other projectiles are refused.
- **Outgoing polarization P of A(a,b)B** is the analyzing power of the inverse reaction B(b,a)A by time reversal.
  Enter P (divide published P·dσ/dΩ by dσ/dΩ at the same energies and angles) as an A_y segment on the inverse
  channel: entrance = the measured exit pair, exit = the measured entrance pair. The former observable 8
  (P·dσ/dΩ) is refused at read time; old files with such lines must be converted this way.
- Worked cases: `tests/13N` segments 11–16 (particle A_y), `tests/13N_capture_ay` (capture).

## Legendre coefficients as data (observable 9)

For papers that publish only a_k = A_k/A_0 of dσ/dΩ_c.m. = A_0 [1 + Σ_k a_k P_k(cos θ_c.m.)]:

- Data segment type 9, `observable="angular-coefficients"`; rows `E_lab k a_k Δa_k`, the order k in the angle column
  (several orders at one energy are several rows). `absolute 1` (`absolute=True`) makes the rows absolute A_k in
  b/sr (A_0 = σ/4π).
- The coefficients are those an `angular-distribution` test segment writes to `AZUREOut_*.acoeff` (a_0 = 1).
- Target integration and resolution fold the absolute A_k, then divide; energy shifts and the analytic gradient work.
- Refused: a free norm on normalised coefficients, an elastic pair (no Coulomb term in the expansion), total capture,
  SUM/RATIO, samples, `total`, UPOS, and type 9 in `<segmentsTest>` (use 3 there).
- Not implemented: emulating a finite-angle Legendre fit (only P_0…P_Lmax at the measured angles). Compare the
  model's A_5+ content first; if it is not negligible, the published a_k are not comparable.
- A gradient check on a folded coefficient segment needs `--no-grid-replacement`.

## Phase shifts

A phase-shift segment is one (J, l, s) wave. For a target with spin, give the channel spin with `spin <s>` after the
shift fields (pyazr `Segment.phase_S`). A run stops with "phase-shift segment #k (J, l) has N channel
spins" when the scheme has several spins for that (J, l) and none is given. Spin-0 targets need nothing. Older
engines ignore the keyword and report half the argument of the summed diagonal elements, which is the phase of no
single wave.

## Composite segments: SUM and RATIO

- A lab-frame differential SUM is compared in the c.m. frame of its base pair; a term to another exit pair (another
  Q value, another compound nucleus of a natural target) carries J_base/J_term. The default output of such a test
  segment is in the base pair's c.m. frame; `kinflag.dat` = 1 gives the lab sum.
- A lab-frame differential RATIO is base/denominator in the lab, each term with its own Jacobian.
- A RATIO whose denominator is exactly 0 is NaN, with one warning per run.
- A SUM component's `Scaling` keeps all digits; a non-numeric scaling refuses the project.
- Component lists from old Qt GUI saves, with entrance/exit keys only (no angle), are read in that layout; too few numbers for
  the count gives a WARNING and the line is read as a plain segment.
- Two segment lines that share a data file are not summed (use SUM; azure2-setup).

## Total capture

Total capture is exit key −1 together with the total-capture data type (3 on a data line, 4 on a test line). Key −1
with type 0 is refused ("the exit key −1 (total capture) needs the total-capture data type"); type 3/4 with a real exit
key means total capture with a WARNING. A type-4 TOTAL_CAPTURE extrapolation equals the sum over the capture pairs
(`tests/total_capture_code`).

## Excitation energy column

The E_x column of `AZUREOut_*` and pyazr's `calculate_excitation_energy` include the entrance pair's excitation
energy E* (a segment on an excited-state pair, e.g. n + ⁷Be*). Builds before 2026-10-28 (expires: 2027-02-01)
omitted E* in the column only; the physics used the full value.

## UPOS secondary γ rays and high partial waves

- UPOS multipole mixing (R_k^LL + 2δ R_k^LL' + δ² R_k^L'L')/(1 + δ²) applies at every Legendre order.
- Legendre orders stop at 20 (`Config::maxLOrder`): a particle channel with l > 10 truncates every differential
  cross section it enters (angle-integrated values are complete); the log prints "WARNING: channel l = ... exceeds
  10". Keep l ≤ 10 for differential data.
- Extrapolation cannot compute UPOS / secondary-γ curves; compute them in data mode on a pseudo-data grid.

## S-wave scattering lengths

Every calculation and fit writes `output/scatteringLengths.out`: a_0 (fm) of every channel spin of every particle
pair at threshold, in the standard, Brune and Park bases; complex when a lower channel is open; "hard sphere" when no
level has a width in that channel spin. Formalism: `docs/source/theory/scattering_length.rst`.

- Measured values enter as an optional block
  `<scatteringLengths>` / `active pairKey s value sigma [imValue imSigma]` / `</scatteringLengths>`
  (s = −1: the coherent sum Σ_s (2s+1) a_s / ((2j₁+1)(2j₂+1)); neutron tables give the bound b, enter b·A/(A+1)).
  Each line is a χ² term in every fitter, the MCMC likelihood and `calculate_chi2_rwa`; `chiSquared.out` gains
  `Total-Scattering-Length-Chi-Squared`.
- pyazr: `m.scattering_lengths(x)`, `m.penalties()["scattering_length"]`, `m.scattering_length_chi2(x)`; the
  measured values are `datum` records of quantity `scattering_length` (`Project.add_datum`, `set`, `remove`).
- Trap when entering published Brune widths γ̃: the `.azr` width is the observed width
  Γ = 2P γ̃² / (1 + Σ γ̃² dS/dE), not 2P γ̃²; dropping the denominator moved a ³He+⁴He a_0 by 0.6 %
  (check `g_int` in `parameters.out`).
- A particle channel exactly at threshold (|E| < 10⁻¹⁰ MeV) takes the exact E → 0 limit (P = 0).
