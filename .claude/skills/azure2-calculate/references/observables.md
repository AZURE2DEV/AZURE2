# Special observables

How to compute and fit the observables that are not a plain σ or dσ/dΩ. The fields of each quantity:
`docs/source/reference/format3.rst`; what each computes: `docs/source/user_guide/segments.rst`; the chooser:
azure2-setup (its data-kinds reference).

## Analyzing power A_y (`quantity analyzing_power`)

The vector analyzing power of a spin-½ projectile on a target of any spin, in particle and capture exits. Declare
a data set of `quantity analyzing_power`, unit `1`; A_y is reported in place of σ, so χ², fitting, plotting and
output files need nothing special (`s.evaluate(x).data[i].model` is A_y).

- **Two routes.** Particle exits use the channel-spin amplitude matrix (Seyler), with the entrance index decomposed
  into projectile and target projections before σ_y acts on the projectile. Capture exits use the Legendre
  coefficients of Seyler and Weller, Phys. Rev. C 20 (1979) 453, Eqs. (20)–(21):
  A_y = Σ_k b_k P_k¹(cos θ) / Σ_k a_k P_k(cos θ), external capture included. `AZURE2_CAPPOL_DEBUG=1` prints, per
  point, the check that Σ a_k P_k reproduces the differential capture cross section. Formalism:
  `docs/source/theory/polarization_theory.rst`, `docs/source/theory/polarization_implementation.rst`.
- **Angles are θ_c.m.**: columns `E_lab[MeV] theta_cm[deg] value error`.
- **A_y lies in [−1, 1] and is often negative.** Take the experimenters' absolute point-to-point uncertainties
  (statistics, background subtraction). Never convert them to relative errors: a point near a zero crossing does not
  have a small uncertainty, and treating it so gives it enormous weight.
- **Free a normalization.** The beam-polarization calibration P scales A_y multiplicatively
  (A_y ∝ (N⁺ − N⁻)/(P (N⁺ + N⁻))), which is exactly a normalization:
  `Project.set(name, norm="1 fit prior 1 <polarization uncertainty>%")`.
- **No target integration** on A_y data sets compared with thin-target data. The target average
  ⟨A_y⟩ = ∫A_y σ dE / ∫σ dE is dominated by the low-energy Rutherford σ where A_y ≈ 0 (≈ 10⁻⁶ on a thick target);
  that is physics. A target-integrated A_y point also switches the whole fit to numerical derivatives, and for
  capture the attenuation coefficients have no counterpart for the P_k¹ numerator.
- **Projectile spin ½ only.** Other projectiles are refused.
- **Outgoing polarization P of A(a,b)B** is the analyzing power of the inverse reaction B(b,a)A by time reversal.
  Enter P (divide published P·dσ/dΩ by dσ/dΩ at the same energies and angles) as an A_y data set of the inverse
  reaction (B(b,a)A).
- Worked cases: `tests/upgrade/T1` (`ay_cm090`, `ay_cm150`), `tests/13N`, `tests/13N_capture_ay` (capture).

## Legendre coefficients as data (`quantity legendre`)

For papers that publish only coefficients of dσ/dΩ_c.m. = A_0 [1 + Σ_k a_k P_k(cos θ_c.m.)]:

- Columns `E_lab[MeV] order value error` (several orders at one energy are several rows); `expansion` says which
  coefficients: `a_over_a0` (a_k), `absolute` (A_k in b/sr, A_0 = σ/4π), `rsl`, `2l2` (EXFOR's LEG/2L2), `4pi`.
- Target integration and resolution fold the absolute A_k, then divide; energy shifts and the analytic gradient work.
- Refused (the engine says why): a free norm on normalised coefficients, an elastic pair (no Coulomb term in the
  expansion), sums and ratios, a secondary γ ray.
- Not implemented: emulating a finite-angle Legendre fit (only P_0…P_Lmax at the measured angles). Compare the
  model's high-order content first; if it is not negligible, the published a_k are not comparable.

## Phase shifts (`quantity phase_shift`)

A phase-shift data set is one (J, l, s) wave: the fields `J`, `l` and `s`. For a target with spin, `s` chooses the
channel spin; several channel spins for one (J, l) without `s` are refused.

## Sums and ratios (a reaction `A + B` or `A / B`)

- A lab-frame differential sum is compared in the c.m. frame of its first term; a term to another exit pair (another
  Q value, another compound nucleus of a natural target) carries J_first/J_term. A lab-frame differential ratio is
  first/denominator in the lab, each term with its own Jacobian. `AZURE2 --print-predictions` shows the tree.
- A ratio whose denominator is exactly 0 is NaN, with one warning per run.
- Two data sets on one data file are not summed (azure2-setup, its data-kinds reference, §3).

## Total capture

A reaction summed over every capture pair of the entrance (`p.reaction(entrance, -1)` writes it) is the total
capture cross section; the rate's exit `-1` is the same sum. One level-matrix solution per energy serves all terms.

## Excitation energy column

The E_x column of `AZUREOut_*` and pyazr's `evaluate(x).data[i].excitation_energies` include the entrance pair's
excitation energy E* (a data set on an excited-state pair, e.g. n + ⁷Be*).

## Secondary γ rays and high partial waves

- A secondary γ ray (`secondary_gamma L <l> Ic <J> delta <δ>`): multipole mixing (R_k^LL + 2δ R_k^LL' + δ² R_k^L'L')/(1 + δ²) applies at every Legendre order.
- Legendre orders stop at 20 (`Config::maxLOrder`): a particle channel with l > 10 truncates every differential
  cross section it enters (angle-integrated values are complete); the log prints "WARNING: channel l = ... exceeds
  10". Keep l ≤ 10 for differential data.
- An extrapolation of the primary reaction is not the secondary-γ curve; evaluate a copy of the data set on a
  pseudo-data file of the grid energies.

## S-wave scattering lengths

Every calculation and fit writes `output/scatteringLengths.out`: a_0 (fm) of every channel spin of every particle
pair at threshold, in the standard, Brune and Park bases; complex when a lower channel is open; "hard sphere" when no
level has a width in that channel spin. Formalism: `docs/source/theory/scattering_length.rst`.

- Measured values are `datum` records of `quantity scattering_length` in `<modelData>` (`Project.add_datum`, `set`,
  `remove`; the coherent sum Σ_s (2s+1) a_s / ((2j₁+1)(2j₂+1)) is a choice of the record; neutron tables give the
  bound b, enter b·A/(A+1)). Each is a term of the objective in every fitter and sampler:
  `s.evaluate(x).terms["model data"]`, per datum `s.evaluate(x).model_data`; `s.scattering_lengths(x)` gives every
  a_0 at x. Comparing a_0 and the effective range with other models: the user-level effective-range skill.
- Trap when entering published Brune widths γ̃: the `.azr` width is the observed width
  Γ = 2P γ̃² / (1 + Σ γ̃² dS/dE), not 2P γ̃²; dropping the denominator moved a ³He+⁴He a_0 by 0.6 %
  (check `g_int` in `parameters.out`).
- A particle channel exactly at threshold (|E| < 10⁻¹⁰ MeV) takes the exact E → 0 limit (P = 0).
